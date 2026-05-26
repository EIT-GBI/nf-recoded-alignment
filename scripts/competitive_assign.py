#!/usr/bin/env python3
"""Per-pair competitive assignment between two name-sorted BAMs.

Faithful port of tools/competitive_assign.awk. Reads two name-sorted BAMs in
lockstep, sums the AS tag of primary alignments per side per qname, and routes
each qname's records to wt-out, rec-out, or rec-out on ties.
"""
import argparse
import sys
import pysam


def is_primary(flag: int) -> bool:
    return (flag & 0x100) == 0 and (flag & 0x800) == 0


def get_AS(read) -> int:
    try:
        return int(read.get_tag('AS'))
    except KeyError:
        return -1


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--wt-body',  required=True)
    p.add_argument('--rec-body', required=True)
    p.add_argument('--wt-out',   required=True)
    p.add_argument('--rec-out',  required=True)
    p.add_argument('--summary',  required=True)
    p.add_argument('--sample',   required=True)
    p.add_argument('--threads', type=int, default=4)
    args = p.parse_args()

    wt_in  = pysam.AlignmentFile(args.wt_body,  'rb', threads=args.threads)
    rec_in = pysam.AlignmentFile(args.rec_body, 'rb', threads=args.threads)
    wt_out  = pysam.AlignmentFile(args.wt_out,  'wb', template=wt_in,  threads=args.threads)
    rec_out = pysam.AlignmentFile(args.rec_out, 'wb', template=rec_in, threads=args.threads)

    wt_iter  = iter(wt_in)
    rec_iter = iter(rec_in)
    wt_pending  = next(wt_iter,  None)
    rec_pending = next(rec_iter, None)

    n_wt = n_rec = n_tie = n_unmapped = 0

    while wt_pending is not None or rec_pending is not None:
        if wt_pending is not None and rec_pending is not None:
            if wt_pending.query_name != rec_pending.query_name:
                sys.stderr.write(
                    f"ERROR: qname order divergence: "
                    f"{wt_pending.query_name} vs {rec_pending.query_name}\n"
                )
                sys.exit(1)
            cur = wt_pending.query_name
        elif wt_pending is not None:
            cur = wt_pending.query_name
        else:
            cur = rec_pending.query_name

        wt_buf, rec_buf = [], []
        sum_wt = sum_rec = 0
        any_wt = any_rec = False

        while wt_pending is not None and wt_pending.query_name == cur:
            wt_buf.append(wt_pending)
            if is_primary(wt_pending.flag):
                as_score = get_AS(wt_pending)
                if as_score >= 0:
                    sum_wt += as_score
                    any_wt = True
            wt_pending = next(wt_iter, None)

        while rec_pending is not None and rec_pending.query_name == cur:
            rec_buf.append(rec_pending)
            if is_primary(rec_pending.flag):
                as_score = get_AS(rec_pending)
                if as_score >= 0:
                    sum_rec += as_score
                    any_rec = True
            rec_pending = next(rec_iter, None)

        if not any_wt and not any_rec:
            n_unmapped += 1
        elif any_wt and (not any_rec or sum_wt > sum_rec):
            for r in wt_buf:
                wt_out.write(r)
            n_wt += 1
        elif any_rec and (not any_wt or sum_rec > sum_wt):
            for r in rec_buf:
                rec_out.write(r)
            n_rec += 1
        else:
            for r in rec_buf:
                rec_out.write(r)
            n_tie += 1

    wt_out.close()
    rec_out.close()
    wt_in.close()
    rec_in.close()

    with open(args.summary, 'w') as f:
        f.write("sample\twt_best\trec_best\ttie\tunmapped\n")
        f.write(f"{args.sample}\t{n_wt}\t{n_rec}\t{n_tie}\t{n_unmapped}\n")


if __name__ == '__main__':
    main()
