# nf-mod-fastp

Nextflow module for fastp (paired-end read trimming/QC). Used as a git submodule by pipelines.

Image: `ghcr.io/eit-gbi/nf-mod-fastp:latest`

## Processes

- `TRIM` — input: `tuple(sample, r1, r2)` → output: `tuple(sample, r1_trimmed, r2_trimmed, html_report, json_report)`

## Use as submodule
```bash
git submodule add https://github.com/eit-gbi/nf-mod-fastp.git modules/fastp
```

Then in your pipeline:
```
include { TRIM } from './modules/fastp/main.nf'
```
