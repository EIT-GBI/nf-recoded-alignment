#!/usr/bin/env nextflow
// Recoding-alignment pipeline.
// Dual alignment (WT + recoded) -> competitive assignment -> calmd + merge -> per-codon recoding.
// WT and recoded references must share identical coordinates and contig name.
//
// Usage:
//   nextflow run main.nf -params-file params.yaml -resume

nextflow.enable.dsl = 2

include { TRIM }                       from './modules/fastp/main.nf'
include { BWA_ALIGN_NSORT }            from './modules/bwa_samtools/main.nf'
include { BIGWIG }                     from './modules/bedtools/main.nf'
include { COMPETITIVE_ASSIGN; VARIANTS;
          RECODING_LANDSCAPE;
          AGGREGATE_ANNDATA }          from './modules/recoding-alignment-utils/main.nf'


workflow {
    def outdir = new File(params.alignment.outdir)
    if (!outdir.exists()) { outdir.mkdirs() }

    wt_ref      = file(params.alignment.wt_ref).toAbsolutePath().toString()
    rec_ref     = file(params.alignment.recoded_ref).toAbsolutePath().toString()
    chrom_sizes = file("${params.alignment.recoded_ref}.chrom.sizes").toAbsolutePath().toString()

    def in_path = file(params.alignment.fastq_dir)
    samples = Channel.fromFilePairs("${in_path}/*_R{1,2}*.fastq*", flat: true)

    trimmed = TRIM(samples).map { sample, r1, r2, _html, _json -> tuple(sample, r1, r2) }

    // Dual alignment: same trimmed reads against both refs, name-sorted in one step.
    wt_in  = trimmed.map { sample, r1, r2 -> tuple(sample, r1, r2, 'wt',  wt_ref)  }
    rec_in = trimmed.map { sample, r1, r2 -> tuple(sample, r1, r2, 'rec', rec_ref) }
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
