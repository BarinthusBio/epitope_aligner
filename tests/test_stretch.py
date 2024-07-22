
import pytest
import numpy as np
import pandas as pd
from epimap import map, stretch, utils

# missing first residue?
# how to get length
# check residues with indexes

sequence = "abcdefghi"
epitopes = pd.DataFrame(
    {
        'start': [1.0,2,2,6],
        'end': [4,5,5,9],
        'seq': ['abc','bcd','bcd','fgh'],
        'mhc_allele': ["x","x","y","z"]
    }
)
# Each epitope needs a length column to indicate how
# much it needs to be stretched
epitopes['length'] = epitopes.seq.apply(len)

@pytest.fixture(params=[0,1], ids=["0index","1index"])
def index(request):
    return request.param

@pytest.fixture(params=[True, False], ids=["endin","endout"])
def includeend(request):
    return request.param

@pytest.fixture
def sequence():
    return utils.random_seq(100)

@pytest.fixture
def epitopes(sequence, index, includeend):
    return utils.random_epitopes(
        sequence,
        n=100,
        epitope_lengths=(5,15),
        index=index,
        includeend=includeend
    )

def test_stretch_size(epitopes):
    stretched_epitopes = stretch.stretch(epitopes)
    total_stretch_size = sum(stretched_epitopes.groupby(["seq", "start"]).size())
    assert sum(epitopes.length) == total_stretch_size

