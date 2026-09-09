// Align paired reads to a single FASTA and emit a name-sorted BAM.
// Uses the nf-mod-bwa image, which bundles both bwa and samtools.

process BWA_ALIGN_NSORT {
    tag "${sample}/${rg_id}"
    cpus params.alignment.threads

    // `ref_files` is the fasta plus its bwa/samtools indexes, all staged side by
    // side; `ref` names the fasta among them. Staging (rather than a bare host
    // path) is what lets a reference built by MAKE_WT_REF/INDEX_REF be used here.
    input:
    tuple val(sample), path(r1), path(r2), val(rg_id), val(ref), path(ref_files)

    output:
    tuple val(sample), val(rg_id), path("${sample}.${rg_id}.ns.bam")

    script:
    """
    bwa mem -t ${task.cpus} \\
      -R "@RG\\tID:${rg_id}\\tSM:${sample}\\tLB:${rg_id}\\tPL:ILLUMINA" \\
      ${ref} ${r1} ${r2} \\
      | samtools sort -n -@${task.cpus} -o ${sample}.${rg_id}.ns.bam -
    """

    stub:
    """
    which bwa && bwa 2>&1 | head -3
    which samtools && samtools --version | head -1
    echo "ref staged as: ${ref}" && ls -l ${ref}*
    touch ${sample}.${rg_id}.ns.bam
    """
}


// bwa-index + faidx a FASTA that the pipeline built itself (a user-supplied
// reference is expected to come pre-indexed, see the README).
process INDEX_REF {
    tag "$name"
    publishDir "${params.alignment.outdir}/reference", mode: 'copy'

    // The fasta is staged in a subdirectory and copied back out under `name`, so
    // the indexes land next to a real file that this task also owns as an output.
    // (`name` is passed in because a staged path's .name keeps that subdirectory.)
    input:
    tuple val(name), path(fasta, stageAs: 'staged/*')

    output:
    tuple val(name), path("${name}*"), emit: ref

    script:
    """
    cp -L ${fasta} ${name}
    bwa index ${name}
    samtools faidx ${name}
    cut -f1,2 ${name}.fai > ${name}.chrom.sizes
    """

    stub:
    """
    which bwa && which samtools
    cp -L ${fasta} ${name}
    touch ${name}.amb ${name}.ann ${name}.bwt \\
          ${name}.pac ${name}.sa ${name}.fai ${name}.chrom.sizes
    """
}
