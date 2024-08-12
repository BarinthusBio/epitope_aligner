from epimap import map, utils
import pandas as pd
import pytest

@pytest.fixture(params=[0, 1], ids=['0index','1index'])
def index(request):
    return request.param


@pytest.fixture(params=[True, False], ids=["endin", "endout"])
def includeend(request):
    return request.param


@pytest.fixture
def sequence():
    return utils.random_seq(100)


@pytest.fixture
def epitopes(sequence, index, includeend):
    return utils.random_epitopes(
        sequence=sequence,
        n=100,
        epitope_lengths=(5, 15),
        index=index,
        includeend=includeend
    )

def test_locate_peptide(epitopes, index, includeend):
    floating_epitopes = map.float_peptides(epitopes, index=index)
    located_pos = [map.locate_peptide(floater, index, includeend) for floater in floating_epitopes]
    epitopes['located_start'] = [pos[0] for pos in located_pos]
    epitopes['located_end'] = [pos[1] for pos in located_pos]
    assert all((epitopes['start'] == epitopes['located_start']) & (epitopes['end'] == epitopes['located_end']))
