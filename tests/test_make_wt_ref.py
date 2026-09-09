"""Tests for scripts/make_wt_ref.py.

The synthetic tests always run. `test_reproduces_reference_wt_fasta` regenerates
the WT FASTA from the recoded FASTA + GenBank named in params.yaml and asserts it
is byte-identical to the `wt_ref` in the same file; it skips when those paths
aren't present on this machine (they're host paths, not repo data).

Run with:  uv run --group dev pytest tests/ -v
"""
import re
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / 'scripts' / 'make_wt_ref.py'
PARAMS = REPO / 'params.yaml'

sys.path.insert(0, str(REPO / 'scripts'))
from make_wt_ref import build_wt  # noqa: E402


def alignment_params():
    """The `alignment:` block of params.yaml, without needing PyYAML."""
    text = PARAMS.read_text()
    block = text.split('alignment:', 1)[1]
    out = {}
    for line in block.splitlines():
        if line and not line[0].isspace() and not line.startswith('#'):
            break                       # next top-level key
        m = re.match(r"\s+(\w+):\s*'([^']*)'", line)
        if m:
            out[m.group(1)] = m.group(2)
    return out


def run_script(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *map(str, args)],
        capture_output=True, text=True, check=True,
    )


# --- synthetic fixtures -------------------------------------------------------

def write_pair(tmp_path, recoded_seq, features):
    """Write a 1-contig recoded FASTA + GenBank. `features` is a list of
    (1-based start, end, label) for misc_features; `label` may be a list to put
    several /label qualifiers on one feature, as real gbks do."""
    fasta = tmp_path / 'rec.fasta'
    gbk = tmp_path / 'rec.gbk'
    fasta.write_text('>chr test contig\n' + recoded_seq + '\n')

    parts = []
    for s, e, label in features:
        parts.append(f"     misc_feature    {s}..{e}\n")
        labels = [label] if isinstance(label, str) else label
        parts.extend(f'                     /label="{one}"\n' for one in labels)
    feature_lines = ''.join(parts)
    seq_lines = ''.join(
        f"{i + 1:>9} {recoded_seq[i:i + 60].lower()}\n"
        for i in range(0, len(recoded_seq), 60)
    )
    gbk.write_text(
        f"LOCUS       chr                    {len(recoded_seq)} bp    DNA     linear   01-JAN-2026\n"
        "FEATURES             Location/Qualifiers\n"
        f"{feature_lines}"
        "ORIGIN\n"
        f"{seq_lines}"
        "//\n"
    )
    return fasta, gbk


def test_flips_plus_strand_codon(tmp_path):
    #                 0-based 3..5 -> AGC (recoded), WT is TCG
    recoded = 'ATG' + 'AGC' + 'GGGTTT'
    fasta, gbk = write_pair(tmp_path, recoded, [(4, 6, 'TCG to AGC')])

    rec, stats = build_wt(fasta, gbk)

    assert str(rec.seq) == 'ATG' + 'TCG' + 'GGGTTT'
    assert (stats['flipped'], stats['rc'], stats['skipped']) == (1, 0, 0)
    assert stats['bp_changed'] == 3          # TCG vs AGC differs at all 3 positions


def test_flips_minus_strand_codon_as_reverse_complement(tmp_path):
    # The label is in CDS reading direction, so a - strand CDS carries the
    # reverse complement of 'AGC' (= GCT) on the plus strand of the reference.
    recoded = 'ATG' + 'GCT' + 'GGGTTT'
    fasta, gbk = write_pair(tmp_path, recoded, [(4, 6, 'TCG to AGC')])

    rec, stats = build_wt(fasta, gbk)

    assert str(rec.seq) == 'ATG' + 'CGA' + 'GGGTTT'   # revcomp of TCG
    assert (stats['flipped'], stats['rc'], stats['skipped']) == (1, 1, 0)


def test_skips_feature_whose_sequence_does_not_match_the_label(tmp_path):
    recoded = 'ATG' + 'TTT' + 'GGGTTT'                # neither AGC nor its RC
    fasta, gbk = write_pair(tmp_path, recoded, [(4, 6, 'TCG to AGC')])

    rec, stats = build_wt(fasta, gbk)

    assert str(rec.seq) == recoded                    # untouched
    assert (stats['flipped'], stats['skipped']) == (0, 1)


def test_keeps_contig_name_length_and_coordinates(tmp_path):
    recoded = 'ATG' + 'AGC' + 'GGGTTT'
    fasta, gbk = write_pair(tmp_path, recoded, [(4, 6, 'TCG to AGC')])

    rec, _ = build_wt(fasta, gbk)

    assert rec.id == 'chr'
    assert len(rec.seq) == len(recoded)


def test_repairs_span_truncated_at_the_end(tmp_path):
    # As in Syn57_evo2 at 187785..187786: a 3-base codon annotated on 2 bases,
    # where the sequence at start..start+3 is the labelled recoded codon.
    recoded = 'ATG' + 'AGC' + 'GGGTTT'
    fasta, gbk = write_pair(tmp_path, recoded, [(4, 5, 'TCG to AGC')])

    rec, stats = build_wt(fasta, gbk)

    assert str(rec.seq) == 'ATG' + 'TCG' + 'GGGTTT'
    assert (stats['flipped'], stats['repaired'], stats['skipped']) == (1, 1, 0)


def test_repairs_span_extended_at_the_start(tmp_path):
    # As in Syn57_evo2 at 675061..675064: the start is a base too early, so the
    # codon sits at the feature's END. Only the end-anchored window matches.
    recoded = 'ATG' + 'AGC' + 'GGGTTT'
    fasta, gbk = write_pair(tmp_path, recoded, [(3, 6, 'TCG to AGC')])

    rec, stats = build_wt(fasta, gbk)

    assert str(rec.seq) == 'ATG' + 'TCG' + 'GGGTTT'
    assert (stats['flipped'], stats['repaired'], stats['skipped']) == (1, 1, 0)


def test_skips_rather_than_aborts_when_no_window_matches(tmp_path):
    # Span disagrees with the label AND neither anchoring matches the sequence:
    # warn and skip, so one bad annotation can't abort a whole run.
    recoded = 'ATG' + 'TTT' + 'GGGTTT'
    fasta, gbk = write_pair(tmp_path, recoded, [(4, 5, 'TCG to AGC')])

    rec, stats = build_wt(fasta, gbk)

    assert str(rec.seq) == recoded                    # untouched
    assert (stats['flipped'], stats['repaired'], stats['skipped']) == (0, 0, 1)


def test_uses_the_first_label_whose_codon_is_actually_present(tmp_path):
    # Syn57_evo2 has 108 features with two alternative recoding labels. The
    # sequence carries GCT, which is the second label's recoded codon.
    recoded = 'ATG' + 'GCT' + 'GGGTTT'
    fasta, gbk = write_pair(tmp_path, recoded,
                            [(4, 6, ['GCG to GCC', 'GCG to GCT'])])

    rec, stats = build_wt(fasta, gbk)

    assert str(rec.seq) == 'ATG' + 'GCG' + 'GGGTTT'
    assert (stats['flipped'], stats['skipped']) == (1, 0)


def test_earlier_matching_label_takes_precedence(tmp_path):
    # When more than one label matches, the first one wins -- this is what keeps
    # the output stable for gbks whose labels overlap in opposite orientations.
    recoded = 'ATG' + 'GCC' + 'GGGTTT'
    fasta, gbk = write_pair(tmp_path, recoded,
                            [(4, 6, ['GCG to GCC', 'AAA to GGC'])])

    rec, _ = build_wt(fasta, gbk)

    assert str(rec.seq) == 'ATG' + 'GCG' + 'GGGTTT'   # not the AAA/GGC label


def test_aborts_when_fasta_and_genbank_lengths_differ(tmp_path):
    fasta, gbk = write_pair(tmp_path, 'ATGAGCGGGTTT', [])
    (tmp_path / 'short.fasta').write_text('>chr\nATGAGC\n')

    with pytest.raises(SystemExit, match='coordinates would not line up'):
        build_wt(tmp_path / 'short.fasta', gbk)


def test_cli_writes_fasta_and_genbank(tmp_path):
    recoded = 'ATG' + 'AGC' + 'GGGTTT'
    fasta, gbk = write_pair(tmp_path, recoded, [(4, 6, 'TCG to AGC')])
    out_fa, out_gbk = tmp_path / 'wt.fasta', tmp_path / 'wt.gbk'

    proc = run_script('--recoded-fasta', fasta, '--genbank', gbk,
                      '--out-fasta', out_fa, '--out-gbk', out_gbk)

    assert 'flipped 1 codon features' in proc.stdout
    assert str(SeqIO.read(out_fa, 'fasta').seq) == 'ATGTCGGGGTTT'
    # the WT gbk keeps the annotations at the same coordinates, on the WT sequence
    wt_gbk = SeqIO.read(out_gbk, 'genbank')
    assert str(wt_gbk.seq) == 'ATGTCGGGGTTT'
    assert [f.type for f in wt_gbk.features] == ['misc_feature']


# --- the real reference pair from params.yaml ---------------------------------

def test_reproduces_reference_wt_fasta(tmp_path):
    p = alignment_params()
    paths = {k: Path(p[k]) for k in ('recoded_ref', 'wt_ref', 'genbank') if k in p}
    missing = [k for k in ('recoded_ref', 'wt_ref', 'genbank')
               if k not in paths or not paths[k].exists()]
    if missing:
        pytest.skip(f"params.yaml reference(s) not available here: {missing}")

    out_fa = tmp_path / 'wt.fasta'
    run_script('--recoded-fasta', paths['recoded_ref'],
               '--genbank', paths['genbank'], '--out-fasta', out_fa)

    assert out_fa.read_bytes() == paths['wt_ref'].read_bytes()
