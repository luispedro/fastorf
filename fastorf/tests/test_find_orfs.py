import genetic_code
import find_orfs
import fasta
from hypothesis import given, strategies as st

def test_find_orfs():
    orfs = find_orfs.find_orfs('ATG' + 'AAG' * 100 + 'TAA')
    assert orfs == [(0, 306, False)]

def test_ecoli_10k():
    '''Test that we can find the ORFs in the E. coli genome (first 10kbp).'''
    [seq] = list(fasta.fasta_iter('./tests/data/NC_000913_10kbp.fna.xz'))
    orfs = find_orfs.find_orfs(seq)
    with open('./tests/data/NC_000913_10kbp.orfs', 'rt') as f:
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

