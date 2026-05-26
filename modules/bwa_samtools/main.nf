// Align paired reads to a single FASTA and emit a name-sorted BAM.
// Uses the nf-mod-bwa image, which bundles both bwa and samtools.

process BWA_ALIGN_NSORT {
    tag "${sample}/${rg_id}"
    cpus params.alignment.threads

    input:
    tuple val(sample), path(r1), path(r2), val(rg_id), val(ref)

    output:
    tuple val(sample), val(rg_id), path("${sample}.${rg_id}.ns.bam")

    script:
    """
    bwa mem -t ${task.cpus} \\
      -R "@RG\\tID:${rg_id}\\tSM:${sample}\\tLB:${rg_id}\\tPL:ILLUMINA" \\
      ${ref} ${r1} ${r2} \\
      | samtools sort -n -@${task.cpus} -o ${sample}.${rg_id}.ns.bam -
    """
}
