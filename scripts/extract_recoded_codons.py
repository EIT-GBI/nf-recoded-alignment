#!/usr/bin/env python3
import argparse
import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq

# Syn57 non-recoded codons (codons to exclude)
NON_RECODED = {"TCT", "TCC", "TCA", "TCG", "GCG", "GCA", "TAG"}


def get_recoded_codons_and_bed(genbank_file, bed_file=None, chrom=None):
    '''Extract all the recoded codon positions from genbank file and create a bed file.
    This assumes that the genbank file has misc_feature annotations for recoded codons.
    chrom: contig name to write in the BED. If None, falls back to the genbank record id.
    '''
    if bed_file is None:
        bed_file = genbank_file.replace('.gb', '.bed')

    record = SeqIO.read(genbank_file, format='genbank')
    if chrom is None:
        chrom = record.id
    recoded_codons = []

    cds_features = [feature for feature in record.features if feature.type == 'CDS']
    recoding_features = [feature for feature in record.features
                         if feature.type == 'misc_feature'
                         and 'to' in str(feature.qualifiers)]

    for cds in cds_features:
        for recoding_feature in recoding_features:
            if cds.location.start <= recoding_feature.location.start and recoding_feature.location.end <= cds.location.end:
                strand = 'f' if cds.location.strand == 1 else 'r'
                codon_info = {
                    'start': recoding_feature.location.start,
                    'end': recoding_feature.location.end,
                    'position': sorted(list(range(recoding_feature.location.start, recoding_feature.location.end))),
                    'strand': strand
                }
                recoded_codons.append(codon_info)

    recoded_codons = sorted(recoded_codons, key=lambda x: x['start'])
    print(f"Total recoded codons found: {len(recoded_codons)}")

    print(f"Writing recoded codons to bed file: {bed_file}")
    with open(bed_file, 'w') as bed:
        for codon in recoded_codons:
            bed.write(f"{chrom}\t{codon['start']}\t{codon['end']}\t{codon['strand']}\n")

    return recoded_codons


def parse_mpileup_bases(bases, ref_base):
    '''Parse the mpileup base string into a list of actual bases'''
    parsed = []
    i = 0
    while i < len(bases):
        base = bases[i]
        if base in '.,':
            parsed.append(ref_base)
            i += 1
        elif base == '^':
            i += 2
        elif base == '$':
            i += 1
        elif base == '+':
            i += 1
            num_str = ''
            while i < len(bases) and bases[i].isdigit():
                num_str += bases[i]
                i += 1
            i += int(num_str)
        elif base == '-':
            i += 1
            num_str = ''
            while i < len(bases) and bases[i].isdigit():
                num_str += bases[i]
                i += 1
            i += int(num_str)
        else:
            parsed.append(base.upper())
            i += 1
    return parsed


def analyze_all_codons_mpileup(bam_file, ref_fasta, bed_file, recoded_codons):
    '''Analyze all codons efficiently using samtools mpileup'''
    cmd = [
        'samtools', 'mpileup',
        '-f', ref_fasta,
        '-l', bed_file,
        '-Q', '0',
        '--output-QNAME',
        bam_file
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("Error running samtools mpileup:", result.stderr)
        return []

    pos_to_codon_start = {}
    for codon_info in recoded_codons:
        for pos in codon_info['position']:
            pos_to_codon_start[pos] = codon_info['start']

    codon_reads = {}
    for codon_info in recoded_codons:
        codon_reads[codon_info['start']] = {'info': codon_info, 'reads': {}}

    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        pos = int(parts[1]) - 1
        ref_base = parts[2]
        bases = parts[4]
        qnames = parts[6].split(',')

        if pos in pos_to_codon_start:
            codon_start = pos_to_codon_start[pos]
            parsed_bases = parse_mpileup_bases(bases, ref_base)
            for base, qname in zip(parsed_bases, qnames):
                if qname not in codon_reads[codon_start]['reads']:
                    codon_reads[codon_start]['reads'][qname] = {}
                codon_reads[codon_start]['reads'][qname][pos] = base

    results = []
    for codon_key in sorted(codon_reads.keys()):
        codon_data = codon_reads[codon_key]
        positions = codon_data['info']['position']
        strand = codon_data['info']['strand']

        total_codons = 0
        recoded_count = 0
        for read_bases in codon_data['reads'].values():
            if all(pos in read_bases for pos in positions):
                total_codons += 1
                codon_seq = ''.join(read_bases[pos] for pos in positions)
                if strand == 'r':
                    codon_seq = str(Seq(codon_seq).reverse_complement())
                if codon_seq not in NON_RECODED:
                    recoded_count += 1

        frequency = round(recoded_count / total_codons, 2) if total_codons > 0 else 0

        results.append({
            'position': codon_key + 1,
            'depth': total_codons,
            'recoded_codons': recoded_count,
            'frequency': frequency
        })

    return results


def process_sample(bam_file, ref_fasta, bed_file, recoded_codons, output_dir, sample_name):
    '''Process a single BAM. Returns a DataFrame; rows are NaN if mpileup yielded nothing.'''
    results = analyze_all_codons_mpileup(bam_file, ref_fasta, bed_file, recoded_codons)

    if not results:
        print(f"No results for sample {sample_name}")
        positions = [c['start'] + 1 for c in recoded_codons]
        df = pd.DataFrame({
            'position': positions,
            'frequency': [np.nan] * len(positions),
            'depth': [0] * len(positions),
            'recoded_codons': [0] * len(positions),
        })
    else:
        df = pd.DataFrame(results)

    csv_dir = os.path.join(output_dir, 'csv')
    os.makedirs(csv_dir, exist_ok=True)
    csv_path = os.path.join(csv_dir, f"{sample_name}_recoding_analysis.csv")
    df.to_csv(csv_path, index=False)
    print(f"Wrote {csv_path}")

    return df


def plot_sample(df, output_dir, sample_name, genome_length=None):
    """Plot the recoding landscape for a sample."""
    df = df.sort_values("position").copy()

    plt.figure(figsize=(15, 5))
    sns.scatterplot(data=df, x='position', y='frequency', s=10, alpha=0.5)
    plt.title(f"Recoding Landscape for {sample_name}")
    plt.xlabel("Genome Position")
    plt.ylabel("Frequency of Non-Recoded Codons")
    plt.ylim(0, 1)
    if genome_length:
        plt.xlim(0, genome_length)

    plots_dir = os.path.join(output_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    output_file = os.path.join(plots_dir, f"{sample_name}_recoding_landscape.png")
    plt.savefig(output_file, bbox_inches="tight")
    plt.close()
    print(f"Wrote {output_file}")


def main():
    print(f'\n[INFO] Running {os.path.basename(__file__)}\n')

    parser = argparse.ArgumentParser(description='Per-sample recoding landscape from one BAM file')
    parser.add_argument('--bam', required=True, help='Path to a single BAM file')
    parser.add_argument('--genbank', required=True, help='Path to GenBank file with recoded codon annotations')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--sample-name', help='Sample name (default: derived from BAM filename)')
    parser.add_argument('--ref-fasta', help='Pre-built reference FASTA (default: derived from --genbank)')
    parser.add_argument('--bed', help='Pre-built BED of recoded codons (default: derived from --genbank)')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    sample_name = args.sample_name or os.path.basename(args.bam).replace('.sorted.bam', '').replace('.bam', '')

    bed_file = args.bed or args.genbank.replace('.gb', '.bed')
    ref_fasta = args.ref_fasta or args.genbank.replace('.gb', '.fasta')

    if not args.ref_fasta:
        SeqIO.write(SeqIO.read(args.genbank, 'genbank'), ref_fasta, 'fasta')

    ref_record = SeqIO.read(ref_fasta, 'fasta')
    chrom = ref_record.id
    genome_length = len(ref_record.seq)

    print(f"Extracting recoded codons from {args.genbank} (chrom={chrom})")
    recoded_codons = get_recoded_codons_and_bed(args.genbank, bed_file, chrom=chrom)

    df = process_sample(args.bam, ref_fasta, bed_file, recoded_codons, args.output_dir, sample_name)

    if not df['frequency'].isna().all():
        plot_sample(df, args.output_dir, sample_name, genome_length)

    print(f'\n[INFO] End of {os.path.basename(__file__)}.\n')


if __name__ == '__main__':
    main()
