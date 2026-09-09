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

Nothing is rewritten unless the sequence agrees with the label: the bases at
the feature position must equal the labelled recoded codon, either directly
(CDS on + strand) or as its reverse complement (CDS on - strand; the label is
in CDS reading direction). Features that match neither are skipped with a
warning -- a gbk can carry stale or overlapping recoding annotations -- and
the totals are reported at the end.

Two quirks of real gbks are handled:
  - a feature may carry several 'XXX to YYY' labels (alternative recodings);
    the first label whose recoded codon is actually present at the position is
    the one that gets flipped
  - a feature's span may disagree with its codon length by a base at one end.
    The codon-sized windows anchored at the feature's start and at its end are
    then both tried, so the sequence decides which end is off; each such repair
    is warned about individually.

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
        'repaired': 0,      # span disagreed with the codon length, resolved by sequence
        'bp_changed': 0,
        'diff_dist': Counter(),   # 1/2/3 bases changed within the codon
        'swap_dist': Counter(),
    }

    for feat in gbk.features:
        if feat.type != 'misc_feature':
            continue
        labels = [(m.group(1).upper(), m.group(2).upper())
                  for label in feat.qualifiers.get('label', [])
                  if (m := LABEL_RE.match(label))]
        if not labels:
            continue

        start = int(feat.location.start)   # BioPython is 0-based
        end = int(feat.location.end)
        span = end - start

        # A feature can carry several 'XXX to YYY' labels (alternative recodings);
        # the first one whose recoded codon is actually present at the position
        # wins. Where the annotated span disagrees with the codon length -- a
        # handful of gbks have a feature truncated or extended by a base -- the
        # codon-sized windows anchored at the feature's start and at its end are
        # tried too, so the sequence itself decides which end is off.
        hit = None
        for wt_codon, rec_codon in labels:
            if len(wt_codon) != len(rec_codon):
                continue
            size = len(rec_codon)
            windows = [(start, start + size)] if size == span else [
                (start, start + size), (end - size, end)]
            rec_rc = str(Seq(rec_codon).reverse_complement())
            for lo, hi in windows:
                if lo < 0 or hi > len(seq):
                    continue
                actual = str(seq[lo:hi]).upper()
                if actual == rec_codon:
                    hit = (lo, hi, wt_codon, rec_codon, wt_codon, False)
                elif actual == rec_rc:
                    hit = (lo, hi, wt_codon, rec_codon,
                           str(Seq(wt_codon).reverse_complement()), True)
                if hit:
                    break
            if hit:
                break

        if hit is None:
            wt_codon, rec_codon = labels[0]
            print(
                f"[WARN] skipped {start + 1}..{end}: actual="
                f"'{str(seq[start:end]).upper()}', label rec='{rec_codon}' or "
                f"rc='{str(Seq(rec_codon).reverse_complement())}'",
                file=sys.stderr,
            )
            stats['skipped'] += 1
            continue

        lo, hi, wt_codon, rec_codon, new_codon, is_rc = hit
        if (lo, hi) != (start, end):
            print(
                f"[WARN] span {start + 1}..{end} ({span} bp) disagrees with the "
                f"{len(rec_codon)} bp label '{wt_codon} to {rec_codon}'; used "
                f"{lo + 1}..{hi}, where the sequence matches",
                file=sys.stderr,
            )
            stats['repaired'] += 1
        if is_rc:
            stats['rc'] += 1

        for i, base in enumerate(new_codon):
            seq[lo + i] = base

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
    print(f"  of which {stats['repaired']} had a span/codon-length mismatch "
          f"resolved from the sequence")
    print(f"skipped {stats['skipped']} features (sequence didn't match label or its RC)")
    print(f"changes-per-codon distribution: {dict(sorted(stats['diff_dist'].items()))}")
    print("top swaps (recoded -> WT):")
    for swap, n in sorted(stats['swap_dist'].items(), key=lambda x: -x[1])[:10]:
        print(f"  {swap}: {n}")
    print(f"wrote: {args.out_fasta}" + (f", {args.out_gbk}" if args.out_gbk else ""))


if __name__ == '__main__':
    main()
