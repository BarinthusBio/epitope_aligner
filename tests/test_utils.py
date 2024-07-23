import pytest

import numpy as np
import pandas as pd
from epimap import map, stretch, utils

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

def test_seq_len(epitopes):
    len_matches = epitopes.apply(lambda row: len(row.seq)==row.length, axis=1)
    assert all(len_matches)

def test_random_epitopes(sequence, epitopes, index):
    floating_epitopes = map.float_peptides(epitopes, index)
    scores = ([map.score_epitope_alignment(epi, sequence)[0] for epi in floating_epitopes])
    assert sum(scores) == len(scores)


