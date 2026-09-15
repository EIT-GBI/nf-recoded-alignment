"""Tests for scripts/extract_recoded_codons.py.

The codon caller reads the reference base straight out of mpileup, which echoes
the FASTA verbatim. A soft-masked (lowercase) reference — Syn57_evo2.fa is one —
used to make every match-to-reference codon lowercase, so it never equalled a
codon in NON_RECODED and every read was scored as recoded. These tests pin the
caller to be case-insensitive about the reference.

Run with:  uv run --group dev pytest tests/ -v
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'scripts'))
from extract_recoded_codons import NON_RECODED, parse_mpileup_bases  # noqa: E402


@pytest.mark.parametrize('ref_base', ['g', 'G'])
def test_matches_are_uppercase_whatever_the_reference_case(ref_base):
    """'.' and ',' mean 'same as reference'; both must come back uppercase."""
    assert parse_mpileup_bases('.,.', ref_base) == ['G', 'G', 'G']


def test_mismatches_are_uppercase_too():
    assert parse_mpileup_bases('aA,', 'g') == ['A', 'A', 'G']


def test_indel_and_read_edge_markers_are_skipped():
    # ^K opens a read (with its mapping quality), $ closes one, +2ag is an
    # insertion of 2 bases, -1a a deletion: none of them contribute a base.
    assert parse_mpileup_bases('^K.$.+2ag.-1a', 'c') == ['C', 'C', 'C']


def test_wt_codon_is_recognised_against_a_lowercase_reference():
    """A read carrying the WT codon at a GCG -> GCT site, on a lowercase ref.

    The read matches the recoded reference at the first two bases and differs at
    the third, so the reconstructed codon is GCG — which is in NON_RECODED, i.e.
    not recoded. Before the fix this came out 'gcG' and was scored as recoded.
    """
    codon = ''.join(parse_mpileup_bases(b, r)[0] for b, r in [('.', 'g'), ('.', 'c'), ('G', 't')])
    assert codon == 'GCG'
    assert codon in NON_RECODED


def test_recoded_codon_is_still_recoded():
    """The same site in a fully recoded clone: every base matches the ref."""
    codon = ''.join(parse_mpileup_bases(b, r)[0] for b, r in [('.', 'g'), ('.', 'c'), ('.', 't')])
    assert codon == 'GCT'
    assert codon not in NON_RECODED
