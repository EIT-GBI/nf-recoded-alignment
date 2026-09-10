#!/usr/bin/env nextflow
// Recoding-alignment pipeline.
// Dual alignment (WT + recoded) -> competitive assignment -> calmd + merge -> per-codon recoding.
// WT and recoded references must share identical coordinates and contig name.
//
// Usage:
//   nextflow run main.nf -params-file params.yaml -resume

nextflow.enable.dsl = 2

// Resolve exactly one FASTQ for a glob; fail loudly on 0 or >1 matches.
def one_fastq(pattern, id, tag) {
    def hits = file(pattern)
    hits = (hits instanceof List) ? hits : [hits]
    if (hits.size() != 1) {
        error "Expected exactly 1 ${tag} FASTQ for ${id} matching ${pattern}, " +
              "found ${hits.size()}: ${hits}"
    }
    return hits[0]
}

// Index files a reference needs on disk. bwa mem needs the first five; samtools
// calmd and both mpileups need .fai; bedGraphToBigWig needs .chrom.sizes, so
// only the recoded reference is required to have that one. A function, not a
// top-level constant: Nextflow 26's parser rejects statements at script level.
def bwa_sidecars() {
    return ['.amb', '.ann', '.bwt', '.pac', '.sa', '.fai']
}

// Sidecars that aren't on disk yet. Empty list => the reference is ready to use.
def missing_sidecars(fasta, suffixes, tag) {
    def f = file(fasta).toAbsolutePath()
    if (!f.exists()) { error "${tag} reference not found: ${f}" }
    return suffixes.findAll { !file("${f}${it}").exists() }
}

// A reference as the processes want it: tuple(fasta name, [fasta + sidecars]),
// staged side by side. Everything matching '<fasta>*' is picked up, so a
// hand-indexed reference and one INDEX_REF built here look the same downstream.
def ref_bundle(fasta) {
    def f = file(fasta).toAbsolutePath()
    def hits = file("${f}*")
    return tuple(f.name, (hits instanceof List) ? hits : [hits])
}

// The chrom.sizes out of a reference bundle, as an absolute path string: BIGWIG
// lives in a shared submodule and takes it as a `val`, so it is not staged. That
// path is either next to the reference or in INDEX_REF's work dir, and
// nextflow.config binds both into the containers.
def chrom_sizes_of(name, files) {
    def hits = (files instanceof List) ? files : [files]
    def hit = hits.find { it.name == "${name}.chrom.sizes" }
    if (!hit) { error "No ${name}.chrom.sizes in the recoded reference bundle" }
    return hit.toAbsolutePath().toString()
}

// INDEX_REF is invoked for both references, so it needs one alias per call site.
include { TRIM }                       from './modules/fastp/main.nf'
include { BWA_ALIGN_NSORT;
          INDEX_REF as INDEX_REC_REF;
          INDEX_REF as INDEX_WT_REF }  from './modules/bwa_samtools/main.nf'
include { BIGWIG }                     from './modules/bedtools/main.nf'
include { COMPETITIVE_ASSIGN; VARIANTS;
          RECODING_LANDSCAPE;
          AGGREGATE_ANNDATA; MAKE_WT_REF } from './modules/recoding-alignment-utils/main.nf'


workflow {
    def outdir = new File(params.alignment.outdir)
    if (!outdir.exists()) { outdir.mkdirs() }

    // Recoded reference. Index it here when any sidecar is missing, rather than
    // failing halfway through the run; the indexes are published to
    // ${outdir}/reference/ so the next run can pick them up from there.
    rec_fasta = file(params.alignment.recoded_ref).toAbsolutePath()
    rec_todo  = missing_sidecars(rec_fasta, bwa_sidecars() + ['.chrom.sizes'], 'Recoded')
    if (rec_todo) {
        log.info "Indexing recoded reference ${rec_fasta.name} (missing: ${rec_todo.join(' ')})"
        rec_ref = INDEX_REC_REF(Channel.value(tuple(rec_fasta.name, rec_fasta))).ref
    } else {
        rec_ref = Channel.value(ref_bundle(rec_fasta))
    }
    chrom_sizes = rec_ref.map { name, files -> chrom_sizes_of(name, files) }

    // WT reference: take the one named in the params file, or build it from the
    // recoded ref + GenBank when `wt_ref` is unset (flipping every 'XXX to YYY'
    // misc_feature back to its WT codon, so coordinates stay identical). The
    // built ref lands in ${outdir}/reference/ — pin it as `wt_ref` in the params
    // file afterwards to skip the rebuild. A WT ref that is named but not
    // indexed gets indexed here too. No .chrom.sizes needed: BIGWIG only ever
    // runs off the recoded reference.
    if (params.alignment.wt_ref) {
        wt_fasta_path = file(params.alignment.wt_ref).toAbsolutePath()
        wt_todo = missing_sidecars(wt_fasta_path, bwa_sidecars(), 'WT')
        if (wt_todo) {
            log.info "Indexing WT reference ${wt_fasta_path.name} (missing: ${wt_todo.join(' ')})"
            wt_ref = INDEX_WT_REF(Channel.value(tuple(wt_fasta_path.name, wt_fasta_path))).ref
        } else {
            wt_ref = Channel.value(ref_bundle(wt_fasta_path))
        }
    } else {
        wt_name  = rec_fasta.baseName + '_wt'
        wt_fasta = MAKE_WT_REF(rec_fasta, file(params.alignment.genbank), wt_name).fasta
        wt_ref = INDEX_WT_REF(wt_fasta.map { f -> tuple(f.name, f) }).ref
    }

    // Samples come either from a samplesheet (many runs in one go) or from a
    // single flat FASTQ directory. Exactly one of the two must be in the params
    // file; the samplesheet wins if somebody sets both.
    if (params.alignment.samplesheet) {
        // One row per sample across all runs. `unique_id` is the sample key end to
        // end (fastq_prefix repeats across runs, so it alone would collide). R1/R2
        // are resolved separately — one R{1,2} glob does NOT reliably order R1 first.
        samples = Channel.fromPath(params.alignment.samplesheet)
            .splitCsv(header: true)
            .map { row ->
                def base = "${row.run_path}/${row.fastq_prefix}"
                def r1 = one_fastq("${base}_R1*.fastq*", row.unique_id, 'R1')
                def r2 = one_fastq("${base}_R2*.fastq*", row.unique_id, 'R2')
                tuple(row.unique_id, r1, r2)
            }
    } else if (params.alignment.fastq_dir) {
        def in_path = file(params.alignment.fastq_dir)
        samples = Channel.fromFilePairs("${in_path}/*_R{1,2}*.fastq*", flat: true)
    } else {
        error "Set either alignment.samplesheet or alignment.fastq_dir in the params file"
    }

    // Undetermined_* is the demultiplexer's leftover bucket, not a sample. A flat
    // fastq_dir glob picks it up like anything else, and because it holds every read
    // that failed indexing it is large and aligns all over the genome — which is what
    // OOM-kills RECODING_LANDSCAPE. Drop it however the sample list was built.
    samples = samples.filter { row ->
        def keep = !(row[0] ==~ /(?i)undetermined.*/)
        if (!keep) { log.info "Skipping non-sample FASTQ pair: ${row[0]}" }
        return keep
    }

    trimmed = TRIM(samples).map { sample, r1, r2, _html, _json -> tuple(sample, r1, r2) }

    // Dual alignment: same trimmed reads against both refs, name-sorted in one step.
    wt_in  = trimmed.combine(wt_ref)
                    .map { sample, r1, r2, name, files -> tuple(sample, r1, r2, 'wt', name, files) }
    rec_in = trimmed.combine(rec_ref)
                    .map { sample, r1, r2, name, files -> tuple(sample, r1, r2, 'rec', name, files) }
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
        rec_ref)

    // The samplesheet doubles as obs metadata; without one, AGGREGATE_ANNDATA
    // falls back to the well info it derives from the sample names.
    AGGREGATE_ANNDATA(
        recoding.csv.map { _sample, csv -> csv }.collect(),
        params.alignment.samplesheet ? file(params.alignment.samplesheet) : [])
}
