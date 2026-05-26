# nf-mod-bedtools

Nextflow module for bedtools + ucsc-bedGraphToBigWig. Used as a git submodule by pipelines. 

Image: `ghcr.io/eit-gbi/nf-mod-bedtools:latest`

## Processes

- `BIGWIG` — inputs: `tuple(sample, label, bam, bai)` + `val chrom_sizes` → output: `path ${sample}.${label}.bw` (optional)

Generates a coverage bigWig from an indexed BAM via `bedtools genomecov -bg` piped into `bedGraphToBigWig`. Empty-coverage samples are skipped.

## Use as submodule
```bash
git submodule add https://github.com/eit-gbi/nf-mod-bedtools.git modules/bedtools
```

Then in your pipeline:
```
include { BIGWIG } from './modules/bedtools/main.nf'
```
