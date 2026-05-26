# recoding_landscape

Nextflow pipeline for processing dual-aligned sequencing data and producing per-sample recoding landscapes.

The pipeline trims reads, aligns each pair against a WT and a recoded reference, performs per-read competitive assignment between the two alignments, merges into one BAM per sample, and extracts per-codon recoding frequencies. All per-sample CSVs are aggregated into a single AnnData (`.h5ad`).

## Requirements

- [Nextflow](https://www.nextflow.io/) (tested with 25.10.4)
- Docker (or Singularity, via the `cluster` profile). Every tool runs inside a container — nothing else needs to be installed on the host.

Clone with submodules so the bwa and samtools container modules come along:

```bash
git clone --recurse-submodules <this repo>
# or, if you've already cloned:
git submodule update --init --recursive
```

## Inputs

- A directory of paired FASTQs named `*_R{1,2}*.fastq*`
- A recoded reference FASTA (used for alignment) — bwa-indexed and `samtools faidx`ed, with a `<ref>.chrom.sizes` file alongside
- A WT reference FASTA at the same coordinates as the recoded reference
- A recoded GenBank file with `misc_feature` annotations marking recoded codons

First-time per-reference setup (run once, on the host or inside the bwa/samtools container):

```bash
bwa index <ref>
samtools faidx <ref>
cut -f1,2 <ref>.fai > <ref>.chrom.sizes
```

## Configuration

All paths and thread counts live in a YAML params file. See [`params.yaml`](params.yaml) and [`params_test.yaml`](params_test.yaml). Keys the pipeline reads:

```yaml
# Directories bind-mounted into every Docker container (refs, fastqs).
bind_mounts:
  - '/path/that/contains/your/refs/and/fastqs'

alignment:
  fastq_dir:    '/path/to/fastq/'
  recoded_ref:  '/path/to/recoded.fasta'
  wt_ref:       '/path/to/wt.fasta'
  genbank:      '/path/to/recoded.gbk'
  outdir:       '/path/to/results/'
  threads:      8
```

`bind_mounts` is the list of host directories made visible inside the containers. Refs and fastqs are passed as paths (not staged), so their host directory must be in this list.

## Container images

The pipeline expects these images on `ghcr.io/eit-gbi/` (or wherever you push them). The `bwa` and `samtools` images come from the submodules; the other four are built from Dockerfiles in [`modules/`](modules/):

| Process | Image |
|---|---|
| TRIM | `ghcr.io/eit-gbi/nf-mod-fastp:latest` |
| BWA_MEM | `ghcr.io/eit-gbi/nf-mod-bwa-samtools:latest` |
| COMPETITIVE_ASSIGN, RECODING_LANDSCAPE, AGGREGATE_ANNDATA | `ghcr.io/eit-gbi/nf-mod-pybio:latest` |
| BIGWIG | `ghcr.io/eit-gbi/nf-mod-bedtools:latest` |
| VARIANTS | `ghcr.io/eit-gbi/nf-mod-bcftools:latest` |

Build them locally for testing:

```bash
docker build -t ghcr.io/eit-gbi/nf-mod-fastp:latest          modules/fastp/
docker build -t ghcr.io/eit-gbi/nf-mod-bcftools:latest       modules/bcftools/
docker build -t ghcr.io/eit-gbi/nf-mod-bedtools:latest      modules/bedtools/
docker build -t ghcr.io/eit-gbi/nf-mod-bwa-samtools:latest   modules/bwa_samtools/
docker build -t ghcr.io/eit-gbi/nf-mod-pybio:latest          modules/pybio/
```

Then `docker push` each one to GHCR (or change the container URIs in [`nextflow.config`](nextflow.config)).

## Run

```bash
nextflow run main.nf -params-file params.yaml -resume
```

Always use `-params-file` (single dash) — `--params-file` would be parsed as a pipeline parameter and silently ignored.

On a SLURM/Singularity cluster:

```bash
nextflow run main.nf -params-file params.yaml -profile cluster -resume
```

## Outputs

Published under `${alignment.outdir}/` (hardlinked from `work/` — zero extra disk):

- `bam/*.final.sorted.bam(.bai)`         — merged competitively-assigned BAMs
- `variants/*.assignment.tsv`             — per-sample wt/rec/tie/unmapped counts
- `variants/*.recoding_state.vcf.gz`      — per-sample variant calls
- `bigwig/*.final.bw`                     — coverage tracks
- `recoding/csv/*_recoding_analysis.csv`  — per-sample per-codon depth + frequency
- `recoding/plots/*_recoding_landscape.png`
- `recoding/anndata/recoding_landscape.h5ad` — aggregated AnnData (X = frequency; layers = depth, recoded_count)

## Cache and disk

- `cache 'lenient'` is set on `RECODING_LANDSCAPE` and `AGGREGATE_ANNDATA` to ignore reference-file mtime jitter.
- Published outputs are hardlinked from `work/`, so they survive when you reclaim the work directory.
- Clean older work dirs (keeping the latest run for `-resume`): `nextflow clean -f -keep-last 1`.

## Pipeline graph

```
TRIM → BWA_MEM (×2: wt, rec) → COMPETITIVE_ASSIGN ┬→ BIGWIG
                                                  ├→ VARIANTS
                                                  └→ RECODING_LANDSCAPE → AGGREGATE_ANNDATA
```

`COMPETITIVE_ASSIGN` runs the assignment in pysam (see [`scripts/competitive_assign.py`](scripts/competitive_assign.py)) and folds the calmd+merge+sort steps into the same task.
