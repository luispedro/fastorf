from fastorf import genetic_code
from fastorf import find_orfs
from fastorf import fasta

from os import path
from hypothesis import given, strategies as st

DATA_DIR = path.join(
        path.dirname(__file__),
        'data')

def _get_file(fname):
    return path.join(DATA_DIR, fname)

def test_find_orfs():
    orfs = find_orfs.find_orfs('ATG' + 'AAG' * 100 + 'TAA')
    assert orfs == [(0, 306, False)]

def test_ecoli_10k():
    '''Test that we can find the ORFs in the E. coli genome (first 10kbp).'''
    [(_, seq)] = list(fasta.fasta_iter(_get_file('NC_000913_10kbp.fna.xz')))
    orfs = find_orfs.find_orfs(seq)
    with open(_get_file('NC_000913_10kbp.orfs'), 'rt') as f:
        expected = [tuple(map(int, line.split())) for line in f]
    orfs.sort()
    expected.sort()
    assert orfs == expected

@given(seq=st.text(min_size=203, max_size=903, alphabet='ATGC'))
def test_reverse_comp(seq):
    rc_seq = genetic_code.reverse_complement(seq)

    orfs = find_orfs.find_orfs(seq)
    rc_orfs = find_orfs.find_orfs(rc_seq)
    assert len(orfs) == len(rc_orfs)
    assert sum(orf.rc for orf in orfs) == sum((not orf.rc) for orf in rc_orfs)


@given(seq=st.text(min_size=203, max_size=903, alphabet='ATGC'),
        pats=st.lists(min_size=1, max_size=7,
            elements=st.text(min_size=1, max_size=5, alphabet='ATGC')))
def test_findall(seq, pats):
    matches = find_orfs.findall(seq, pats)
    matches = set(matches)
    for i in range(len(seq)):
        is_match = any(seq[i:].startswith(p) for p in pats)
        assert (i in matches) == is_match
