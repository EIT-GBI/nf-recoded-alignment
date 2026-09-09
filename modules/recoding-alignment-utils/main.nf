// Pipeline-specific processes that all share the recoding-alignment-utils image
// (python + pysam + biopython + pandas + anndata + samtools + bcftools).

process COMPETITIVE_ASSIGN {
    tag "$sample"
    cpus params.alignment.threads
    publishDir "${params.alignment.outdir}/variants", mode: 'link', pattern: "*.assignment.tsv"
    publishDir "${params.alignment.outdir}/bam",      mode: 'link', pattern: "*.final.sorted.bam*"

    input:
    tuple val(sample), path(wt_ns_bam), path(rec_ns_bam)
    val recoded_ref

    output:
    tuple val(sample), val('final'),
          path("${sample}.final.sorted.bam"),
          path("${sample}.final.sorted.bam.bai"), emit: final_bam
    path "${sample}.assignment.tsv",              emit: assignment

    script:
    """
    set -euo pipefail

    # Competitive assignment via pysam: name-sorted BAMs in, name-sorted BAMs out (no SAM intermediates).
    python ${workflow.projectDir}/scripts/competitive_assign.py \\
      --wt-body ${wt_ns_bam} \\
      --rec-body ${rec_ns_bam} \\
      --wt-out wt.bam \\
      --rec-out rec.bam \\
      --summary ${sample}.assignment.tsv \\
      --sample ${sample} \\
      --threads ${task.cpus}

    # Coord-sort, calmd wt, merge with rec, sort, index — no merged.bam intermediate.
    samtools sort -@${task.cpus} -o wt.sorted.bam wt.bam
    samtools sort -@${task.cpus} -o rec.sorted.bam rec.bam
    samtools calmd -b wt.sorted.bam ${recoded_ref} 2>/dev/null > wt.calmd.bam
    samtools merge -@${task.cpus} -f - rec.sorted.bam wt.calmd.bam \\
      | samtools sort -@${task.cpus} -o ${sample}.final.sorted.bam -
    samtools index ${sample}.final.sorted.bam

    rm wt.bam rec.bam wt.sorted.bam rec.sorted.bam wt.calmd.bam
    """

    stub:
    """
    which samtools && samtools --version | head -1
    which python && python -c 'import pysam; print("pysam", pysam.__version__)'
    echo "ref path visible: ${recoded_ref}" && ls "${recoded_ref}"
    touch ${sample}.final.sorted.bam ${sample}.final.sorted.bam.bai ${sample}.assignment.tsv
    """
}


process VARIANTS {
    tag "$sample"
    cpus params.alignment.threads
    publishDir "${params.alignment.outdir}/variants", mode: 'link'

    input:
    tuple val(sample), val(label), path(final_bam), path(final_bai)
    val recoded_ref

    output:
    tuple path("${sample}.recoding_state.vcf.gz"),
          path("${sample}.recoding_state.vcf.gz.csi")

    script:
    """
    bcftools mpileup -f ${recoded_ref} -a AD,DP ${final_bam} \\
      | bcftools call --ploidy 1 -mv -Oz -o ${sample}.recoding_state.vcf.gz
    bcftools index ${sample}.recoding_state.vcf.gz
    """

    stub:
    """
    which bcftools && bcftools --version | head -1
    touch ${sample}.recoding_state.vcf.gz ${sample}.recoding_state.vcf.gz.csi
    """
}


process RECODING_LANDSCAPE {
    tag "$sample"
    cache 'lenient'
    publishDir "${params.alignment.outdir}/recoding", mode: 'link'

    input:
    tuple val(sample), path(bam), path(bai)
    path genbank
    path ref_fasta

    output:
    tuple val(sample), path("csv/${sample}_recoding_analysis.csv"), emit: csv
    path "plots/${sample}_recoding_landscape.png",                  emit: plot, optional: true

    script:
    """
    python ${workflow.projectDir}/scripts/extract_recoded_codons.py \\
        --bam ${bam} \\
        --genbank ${genbank} \\
        --ref-fasta ${ref_fasta} \\
        --sample-name ${sample} \\
        --output-dir .
    """

    stub:
    """
    which python && python -c 'import pysam, Bio, pandas; print("imports ok")'
    ls "${genbank}" "${ref_fasta}"
    mkdir -p csv plots
    touch csv/${sample}_recoding_analysis.csv
    """
}


process AGGREGATE_ANNDATA {
    tag 'aggregate'
    cache 'lenient'
    publishDir "${params.alignment.outdir}/recoding/anndata", mode: 'link'

    input:
    path csvs
    path metadata

    output:
    path 'recoding_landscape.h5ad'

    script:
    """
    python ${workflow.projectDir}/scripts/aggregate_anndata.py \\
        --csvs ${csvs} \\
        --metadata ${metadata} \\
        --output recoding_landscape.h5ad
    """

    stub:
    """
    which python && python -c 'import anndata; print("anndata", anndata.__version__)'
    echo "received \$(ls *.csv 2>/dev/null | wc -l) csv(s); metadata=${metadata}"
    touch recoding_landscape.h5ad
    """
}


process MAKE_WT_REF {
    tag "$wt_name"
    cache 'lenient'
    // The fasta itself is published by INDEX_REF (with its indexes alongside).
    publishDir "${params.alignment.outdir}/reference", mode: 'copy', pattern: "*.gbk"

    input:
    path recoded_fasta
    path genbank
    val wt_name

    output:
    path "${wt_name}.fasta", emit: fasta
    path "${wt_name}.gbk",   emit: gbk

    script:
    """
    python ${workflow.projectDir}/scripts/make_wt_ref.py \\
        --recoded-fasta ${recoded_fasta} \\
        --genbank ${genbank} \\
        --out-fasta ${wt_name}.fasta \\
        --out-gbk ${wt_name}.gbk
    """

    stub:
    """
    which python && python -c 'import Bio; print("biopython", Bio.__version__)'
    ls "${recoded_fasta}" "${genbank}"
    touch ${wt_name}.fasta ${wt_name}.gbk
    """
}
