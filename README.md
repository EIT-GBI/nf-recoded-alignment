# nf-recoded-alignment

Nextflow pipeline: dual alignment (WT + recoded) → competitive assignment → per-codon recoding landscape.

## Requirements

- Nextflow (25.10+)
- Docker (or Singularity via `-profile cluster`)

```bash
git clone --recurse-submodules <repo>
```

## Inputs

- Paired FASTQs named `*_R{1,2}*.fastq*`
- Recoded FASTA (bwa-indexed, `samtools faidx`'d, with `<ref>.chrom.sizes`)
- WT FASTA at identical coordinates
- Recoded GenBank with `misc_feature` codon annotations

Per-reference setup, run once:

```bash
bwa index <ref> && samtools faidx <ref>
cut -f1,2 <ref>.fai > <ref>.chrom.sizes
```

## Configure

Edit [`params.yaml`](params.yaml) — paths to FASTQs, refs, outdir, and `bind_mounts` (host dirs to expose to containers).

## Run

```bash
nextflow run main.nf -params-file params.yaml -resume
# cluster:
nextflow run main.nf -params-file params.yaml -profile cluster -resume
```

Use `-params-file` (single dash).

## Outputs

Published under `${alignment.outdir}/`:

- `bam/*.final.sorted.bam` — merged competitively-assigned BAMs
- `variants/*.assignment.tsv`, `*.recoding_state.vcf.gz`
- `bigwig/*.final.bw`
- `recoding/csv/*_recoding_analysis.csv`, `recoding/anndata/recoding_landscape.h5ad`

## Graph

```
TRIM → BWA_ALIGN_NSORT (×2: wt, rec) → COMPETITIVE_ASSIGN ┬→ BIGWIG
                                                          ├→ VARIANTS
                                                          └→ RECODING_LANDSCAPE → AGGREGATE_ANNDATA
```
