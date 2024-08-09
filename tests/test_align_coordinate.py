from epimap import map, utils
import pytest

@pytest.fixture(params=[0,1], ids=["0index","1index"])
def index(request):
    return request.param

@pytest.fixture(params=[True, False], ids=["endin","endout"])
def includeend(request):
    return request.param

@pytest.fixture
def sequence():
    return utils.random_seq(50)

@pytest.fixture
def epitopes(sequence, index, includeend):
    return utils.random_epitopes(
        sequence,
        n=100,
        epitope_lengths=(5,15),
        index=index,
        includeend=includeend
    )

@pytest.fixture
def aligned_seq(sequence):
    return utils.random_gaps(
        seq=sequence,
        gap_prob=0.2,
        gap_size_interval=(1,5)
    )

def test_align_coordinate(epitopes, aligned_seq, index, includeend):
    epitopes['newstart'] = epitopes.start.apply(map.align_coords, aligned_seq=aligned_seq, index=index)
    epitopes['newend'] = epitopes.end.apply(map.align_coords, aligned_seq=aligned_seq, index=index)
    epitopes['aligned_seq'] = epitopes.apply(lambda x: aligned_seq[x.newstart-index:x.newend-index+includeend], axis=1)
    assert all(epitopes.seq == epitopes.aligned_seq.str.replace("-",""))

def test_unalign_coordinate(epitopes, aligned_seq, index, includeend, sequence):
    epitopes['newstart'] = epitopes.start.apply(map.align_coords, aligned_seq=aligned_seq, index=index)
    epitopes['newend'] = epitopes.end.apply(map.align_coords, aligned_seq=aligned_seq, index=index)
    epitopes['oldstart'] = epitopes.newstart.apply(map.unalign_coordinate, aligned_seq=aligned_seq, index=index)
    epitopes['oldend'] = epitopes.newend.apply(map.unalign_coordinate, aligned_seq=aligned_seq, index=index)
    epitopes['unaligned_seq'] = epitopes.apply(lambda x: sequence[x.oldstart-index:x.oldend-index+includeend], axis=1)
    assert all(epitopes.seq == epitopes.unaligned_seq)
