from epimap import map
import pandas as pd
import pytest

def test_float_peptides():
    sequence = "abcdefghi"
    epitopes = pd.DataFrame(
        {
            'START': [1.0,2,3,4,6],
            'END': [4,5,6,7,9],
            'SEQ': ['abc','bcd','cde','def','fgh']
        }
    )
    floating_peptides = map.float_peptides(epitopes, index=1, start_col="START", seq_col="SEQ")
    assert floating_peptides == ['abc', '-bcd', '--cde', '---def', '-----fgh']


sequence = "ABCDEFGHI"
epitopes = pd.DataFrame(
    {
        'START': [1.0,2,3,4,6],
        'END': [4,5,6,7,9],
        'SEQ': ['abc','xcd','cxe','dex','xxh']
    }
)

def test_score_1_toupper():
    floating_peptides = map.float_peptides(
        epitopes,
        index=1,
        start_col="START",
        seq_col="SEQ"
    )
    scores = [map.score_epitope_alignment(epi, sequence, toupper=True)[0] for epi in floating_peptides]
    true_scores = [1,2/3,2/3,2/3,1/3]
    assert scores == pytest.approx(true_scores)

def test_score_length_assertion():
    seq = "ABC"
    epi = "-BCD"
    with pytest.raises(AssertionError):
        map.score_epitope_alignment(epi, seq)

# To write...
# Test score with/out upper
# 0/1 index
# and combinations

# sequence = "ABCDEFGHI"
# epitopes = pd.DataFrame(
#     {
#         'START': [1.0,2,3,4,6],
#         'END': [4,5,6,7,9],
#         'SEQ': ['abc','bcd','cde','def','fgh']
#     }
# )
# floating_peptides = map.float_peptides(epitopes, index=1, start_col="START", seq_col="SEQ")
# [map.score_epitope_alignment(epi, sequence)[0] for epi in floating_peptides]
