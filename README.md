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
- Recoded GenBank with `misc_feature` codon annotations
- WT FASTA at identical coordinates — **optional**, see below

Per-reference setup, run once:

```bash
bwa index <ref> && samtools faidx <ref>
cut -f1,2 <ref>.fai > <ref>.chrom.sizes
```

### WT reference

Leave `alignment.wt_ref` out of the params file and the pipeline builds it
(`MAKE_WT_REF` → `INDEX_REF`) from the recoded FASTA + GenBank, flipping every
`XXX to YYY` `misc_feature` back to its WT codon in place — so contig name,
length and coordinates stay identical to the recoded ref. Features whose
sequence matches neither the labelled recoded codon nor its reverse complement
are skipped with a warning (see the task's `.command.out` for the counts).

The result is published, bwa-indexed, to `${alignment.outdir}/reference/`
alongside a WT GenBank. Point `wt_ref` at it on later runs to skip the rebuild:

```yaml
alignment:
  wt_ref: '<outdir>/reference/<recoded basename>_wt.fasta'
```

The same thing standalone:

```bash
scripts/make_wt_ref.py --recoded-fasta <recoded.fasta> --genbank <recoded.gbk> \
    --out-fasta <wt.fasta> --out-gbk <wt.gbk>
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
- `reference/*_wt.fasta*`, `*_wt.gbk` — only when the WT reference was built here

## Tests

```bash
uv run --group dev pytest tests/ -v
```

`tests/test_make_wt_ref.py` covers the WT-reference build, including a check that
regenerating the WT FASTA from the reference pair in `params.yaml` reproduces the
existing `wt_ref` byte-for-byte (skipped when those paths aren't reachable).

## Graph

```
[MAKE_WT_REF → INDEX_REF]  (only when `wt_ref` is unset)

TRIM → BWA_ALIGN_NSORT (×2: wt, rec) → COMPETITIVE_ASSIGN ┬→ BIGWIG
                                                          ├→ VARIANTS
                                                          └→ RECODING_LANDSCAPE → AGGREGATE_ANNDATA
```
