#!/usr/bin/env python3
"""
Build a 'WT' reference from the recoded reference by flipping every
'XXX to YYY' misc_feature in the recoded GenBank back to its WT codon at the
annotated position.

The sequence is taken from the recoded FASTA (not from the GenBank), and the
contig name/description are carried over unchanged, so the two references share
identical contig name, length and coordinates and downstream alignments are
interchangeable position-wise. The GenBank is used only for the annotations,
and its sequence must match the FASTA's length.

Sanity checks per feature:
  - codon lengths match (len(WT) == len(recoded) == feature span)
  - the current sequence at the feature position equals the labelled
    recoded codon, either directly (CDS on + strand) or as its reverse
    complement (CDS on - strand; the label is in CDS reading direction)
Length mismatches abort the run. Sequence mismatches print a warning and
the feature is skipped (the gbk can carry stale or overlapping recoding
annotations); the script reports the total skip count at the end.

Usage:
    scripts/make_wt_ref.py --recoded-fasta <recoded.fasta> \\
        --genbank <recoded.gbk> --out-fasta <wt.fasta> [--out-gbk <wt.gbk>]
"""
import argparse
import re
import sys
from collections import Counter

from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq

LABEL_RE = re.compile(r'^\s*([ACGT]+)\s+to\s+([ACGT]+)\s*$', re.IGNORECASE)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--recoded-fasta', required=True,
                   help='recoded reference FASTA (single contig)')
    p.add_argument('--genbank', required=True,
                   help='recoded GenBank carrying the "XXX to YYY" misc_features')
    p.add_argument('--out-fasta', required=True, help='WT FASTA to write')
    p.add_argument('--out-gbk', help='optional WT GenBank to write')
    return p.parse_args(argv)


def build_wt(recoded_fasta, genbank):
    """Return (wt_record, stats). wt_record carries the recoded contig's
    name/description with every annotated codon flipped back to WT."""
    record = SeqIO.read(recoded_fasta, 'fasta')
    gbk = SeqIO.read(genbank, 'genbank')

    if len(gbk.seq) != len(record.seq):
        sys.exit(
            f"[ERROR] length mismatch: {recoded_fasta} is {len(record.seq)} bp "
            f"but {genbank} is {len(gbk.seq)} bp — coordinates would not line up"
        )

    seq = MutableSeq(str(record.seq))
    stats = {
        'flipped': 0,
        'rc': 0,            # features flipped on the reverse strand
        'skipped': 0,
        'bp_changed': 0,
        'diff_dist': Counter(),   # 1/2/3 bases changed within the codon
        'swap_dist': Counter(),
    }

    for feat in gbk.features:
        if feat.type != 'misc_feature':
            continue
        wt_codon = rec_codon = None
        for label in feat.qualifiers.get('label', []):
            m = LABEL_RE.match(label)
            if m:
                wt_codon = m.group(1).upper()
                rec_codon = m.group(2).upper()
                break
        if wt_codon is None:
            continue

        start = int(feat.location.start)   # BioPython is 0-based
        end = int(feat.location.end)
        span = end - start

        if len(wt_codon) != span or len(rec_codon) != span:
            sys.exit(
                f"[ERROR] length mismatch at {start + 1}..{end}: "
                f"span={span}  WT='{wt_codon}' ({len(wt_codon)})  "
                f"rec='{rec_codon}' ({len(rec_codon)})"
            )

        actual = str(seq[start:end]).upper()
        rec_rc = str(Seq(rec_codon).reverse_complement())
        if actual == rec_codon:
            new_codon = wt_codon
        elif actual == rec_rc:
            new_codon = str(Seq(wt_codon).reverse_complement())
            stats['rc'] += 1
        else:
            print(
                f"[WARN] skipped {start + 1}..{end}: "
                f"actual='{actual}', label rec='{rec_codon}' or rc='{rec_rc}'",
                file=sys.stderr,
            )
            stats['skipped'] += 1
            continue

        for i, base in enumerate(new_codon):
            seq[start + i] = base

        diffs = sum(1 for a, b in zip(wt_codon, rec_codon) if a != b)
        stats['flipped'] += 1
        stats['bp_changed'] += diffs
        stats['diff_dist'][diffs] += 1
        stats['swap_dist'][f"{rec_codon}->{wt_codon}"] += 1

    record.seq = Seq(str(seq))
    return record, stats


def main(argv=None):
    args = parse_args(argv)

    wt_record, stats = build_wt(args.recoded_fasta, args.genbank)
    SeqIO.write(wt_record, args.out_fasta, 'fasta')

    if args.out_gbk:
        # Re-read the gbk and swap in the WT sequence so the WT GenBank keeps
        # every annotation (the codon labels included) at the same coordinates.
        gbk = SeqIO.read(args.genbank, 'genbank')
        gbk.seq = wt_record.seq
        SeqIO.write(gbk, args.out_gbk, 'genbank')

    print(f"flipped {stats['flipped']} codon features "
          f"({stats['bp_changed']} bp changed total)")
    print(f"  on + strand: {stats['flipped'] - stats['rc']}, on - strand: {stats['rc']}")
    print(f"skipped {stats['skipped']} features (sequence didn't match label or its RC)")
    print(f"changes-per-codon distribution: {dict(sorted(stats['diff_dist'].items()))}")
    print("top swaps (recoded -> WT):")
    for swap, n in sorted(stats['swap_dist'].items(), key=lambda x: -x[1])[:10]:
        print(f"  {swap}: {n}")
    print(f"wrote: {args.out_fasta}" + (f", {args.out_gbk}" if args.out_gbk else ""))


if __name__ == '__main__':
    main()
