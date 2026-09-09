#!/usr/bin/env nextflow
// Recoding-alignment pipeline.
// Dual alignment (WT + recoded) -> competitive assignment -> calmd + merge -> per-codon recoding.
// WT and recoded references must share identical coordinates and contig name.
//
// Usage:
//   nextflow run main.nf -params-file params.yaml -resume

nextflow.enable.dsl = 2

// A reference as BWA_ALIGN_NSORT wants it: tuple(fasta name, [fasta + indexes]).
// Everything matching '<fasta>*' is picked up, i.e. the bwa/samtools indexes the
// README asks you to build once per reference.
def ref_bundle(fasta, tag) {
    def f = file(fasta).toAbsolutePath()
    if (!f.exists()) { error "${tag} reference not found: ${f}" }
    def hits = file("${f}*")
    hits = (hits instanceof List) ? hits : [hits]
    return tuple(f.name, hits)
}

include { TRIM }                       from './modules/fastp/main.nf'
include { BWA_ALIGN_NSORT; INDEX_REF } from './modules/bwa_samtools/main.nf'
include { BIGWIG }                     from './modules/bedtools/main.nf'
include { COMPETITIVE_ASSIGN; VARIANTS;
          RECODING_LANDSCAPE;
          AGGREGATE_ANNDATA; MAKE_WT_REF } from './modules/recoding-alignment-utils/main.nf'


workflow {
    def outdir = new File(params.alignment.outdir)
    if (!outdir.exists()) { outdir.mkdirs() }

    rec_ref     = file(params.alignment.recoded_ref).toAbsolutePath().toString()
    chrom_sizes = file("${params.alignment.recoded_ref}.chrom.sizes").toAbsolutePath().toString()

    // WT reference: take the one named in the params file, or build it from the
    // recoded ref + GenBank when `wt_ref` is unset (flipping every 'XXX to YYY'
    // misc_feature back to its WT codon, so coordinates stay identical). The
    // built ref lands in ${outdir}/reference/ — pin it as `wt_ref` in the params
    // file afterwards to skip the rebuild.
    if (params.alignment.wt_ref) {
        wt_ref = Channel.value(ref_bundle(params.alignment.wt_ref, 'WT'))
    } else {
        wt_name  = file(params.alignment.recoded_ref).baseName + '_wt'
        wt_fasta = MAKE_WT_REF(file(params.alignment.recoded_ref),
                               file(params.alignment.genbank),
                               wt_name).fasta
        wt_ref = INDEX_REF(wt_fasta.map { f -> tuple(f.name, f) }).ref
    }

    def in_path = file(params.alignment.fastq_dir)
    samples = Channel.fromFilePairs("${in_path}/*_R{1,2}*.fastq*", flat: true)

    trimmed = TRIM(samples).map { sample, r1, r2, _html, _json -> tuple(sample, r1, r2) }

    // Dual alignment: same trimmed reads against both refs, name-sorted in one step.
    def (rec_name, rec_files) = ref_bundle(params.alignment.recoded_ref, 'Recoded')
    wt_in  = trimmed.combine(wt_ref)
                    .map { sample, r1, r2, name, files -> tuple(sample, r1, r2, 'wt', name, files) }
    rec_in = trimmed.map { sample, r1, r2 -> tuple(sample, r1, r2, 'rec', rec_name, rec_files) }
    aligned = BWA_ALIGN_NSORT(wt_in.mix(rec_in))

    wt_ns  = aligned.filter { _s, rg_id, _b -> rg_id == 'wt'  }.map { s, _rg, b -> tuple(s, b) }
    rec_ns = aligned.filter { _s, rg_id, _b -> rg_id == 'rec' }.map { s, _rg, b -> tuple(s, b) }

    // Competitive assignment + calmd + merge in one step.
    paired_ns = wt_ns.join(rec_ns)
    final_bam = COMPETITIVE_ASSIGN(paired_ns, rec_ref).final_bam

    BIGWIG(final_bam, chrom_sizes)
    VARIANTS(final_bam, rec_ref)

    recoding = RECODING_LANDSCAPE(
        final_bam.map { sample, _label, bam, bai -> tuple(sample, bam, bai) },
        file(params.alignment.genbank),
        file(params.alignment.recoded_ref))

    AGGREGATE_ANNDATA(recoding.csv.map { _sample, csv -> csv }.collect())
}
