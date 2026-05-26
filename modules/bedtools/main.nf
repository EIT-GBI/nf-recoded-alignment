process BIGWIG {
    tag "${sample}/${label}"
    publishDir "${params.alignment.outdir}/bigwig", mode: 'link'

    input:
    tuple val(sample), val(label), path(bam), path(bai)
    val chrom_sizes

    output:
    path "${sample}.${label}.bw", optional: true

    script:
    """
    set -euo pipefail
    bedtools genomecov -ibam ${bam} -bg | sort -k1,1 -k2,2n > tmp.bg
    if [ -s tmp.bg ]; then
        bedGraphToBigWig tmp.bg ${chrom_sizes} ${sample}.${label}.bw
    else
        echo "No coverage records for ${sample} — skipping bigwig"
    fi
    """
}
