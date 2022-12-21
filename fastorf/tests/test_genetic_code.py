import genetic_code
from hypothesis import given, strategies as st

def test_genetic_code():
    assert genetic_code.is_start('ATG')
    assert not genetic_code.is_start('AAG')
    assert not genetic_code.is_stop('ATG')
    assert  genetic_code.is_stop('TAA')

    assert genetic_code.is_start_reverse_complement(reverse_complement('ATG'))


@given(codon=st.lists(min_size=3, max_size=3, elements=st.sampled_from('ATGC')))
def test_reverse_comp(codon):
    codon = ''.join(codon)
    assert genetic_code.is_start(codon) == genetic_code.is_start_reverse_complement(reverse_complement(codon))
    assert genetic_code.is_stop(codon) == genetic_code.is_stop_reverse_complement(reverse_complement(codon))


