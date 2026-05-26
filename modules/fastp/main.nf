process TRIM {
    tag "$sample"
    cpus params.alignment.threads
    publishDir "${params.alignment.outdir}/trimmed", mode: 'link', pattern: '*.{html,json}'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample),
          path("${sample}_R1.trimmed.fastq.gz"),
          path("${sample}_R2.trimmed.fastq.gz"),
          path("${sample}.fastp.html"),
          path("${sample}.fastp.json")

    script:
    """
    fastp -w ${task.cpus} -i ${r1} -I ${r2} \\
      -o ${sample}_R1.trimmed.fastq.gz -O ${sample}_R2.trimmed.fastq.gz \\
      -h ${sample}.fastp.html -j ${sample}.fastp.json
    """
}
