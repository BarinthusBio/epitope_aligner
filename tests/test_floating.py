from epitope_aligner import map, utils
import pandas as pd
import pytest

def test_align_float_peptide():
    assert map._align_float(
        start=2,
        seq="CDEF",
        parent_seq="ABC--DE-FG",
        index=0
    ) == "C--DE-F"

def test_float_peptides():
    sequence = "abcdefghi"
    epitopes = pd.DataFrame(
        {
            'START': [1.0,2,3,4,6],
            'END': [4,5,6,7,9],
            'SEQ': ['abc','bcd','cde','def','fgh']
        }
    )
    floating_peptides = map.float_peptides(epitopes, parent_seq=sequence, index=1, start_col="START", seq_col="SEQ")
    assert floating_peptides == ['abc', '-bcd', '--cde', '---def', '-----fgh']


def test_float_peptides_none_parent():
    sequence = "abcdefghi"
    epitopes = pd.DataFrame(
        {
            'START': [1.0,2,3,4,6],
            'END': [4,5,6,7,9],
            'SEQ': ['abc','bcd','cde','def','fgh']
        }
    )
    floating_peptides = map.float_peptides(epitopes, parent_seq=None, index=1, start_col="START", seq_col="SEQ")
    assert floating_peptides == ['abc', '-bcd', '--cde', '---def', '-----fgh']


def test_score_1_toupper():
    sequence = "ABCDEFGHI"
    epitopes = pd.DataFrame(
        {
            'START': [1.0,2,3,4,6],
            'END': [4,5,6,7,9],
            'SEQ': ['abc','xcd','cxe','dex','xxh']
        }
    )
    floating_peptides = map.float_peptides(
        epitopes,
        parent_seq=sequence,
        index=1,
        start_col="START",
        seq_col="SEQ"
    )
    scores = [map.score_epitope_alignment(epi, sequence, toupper=True)[0] for epi in floating_peptides]
    true_scores = [1, 2/3, 2/3, 2/3, 1/3]
    assert scores == pytest.approx(true_scores)

def test_score_1_noupper():
    sequence = "ABCDEFGHI"
    epitopes = pd.DataFrame(
        {
            'START': [1.0,2,3,4,6],
            'END': [4,5,6,7,9],
            'SEQ': ['abc','xcd','cxe','dex','xxh']
        }
    )
    floating_peptides = map.float_peptides(
        epitopes,
        parent_seq=sequence,
        index=1,
        start_col="START",
        seq_col="SEQ"
    )
    scores = [map.score_epitope_alignment(epi, sequence, toupper=False)[0] for epi in floating_peptides]
    true_scores = [0, 0, 0, 0, 0]
    assert scores == pytest.approx(true_scores)

def test_score_length_assertion():
    seq = "ABC"
    epi = "-BCD"
    with pytest.raises(AssertionError):
        map.score_epitope_alignment(epi, seq)


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

def test_float_random_peptides(epitopes, sequence, index):
    floating_peptides = map.float_peptides(table=epitopes, parent_seq=sequence, index=index)
    scores = [map.score_epitope_alignment(epi, sequence) for epi in floating_peptides]
    perfect_score = [1==x[0] for x in scores]
    all_match = [sum(x[1])==len(x[1]) for x in scores]
    assert all(perfect_score) and all(all_match)


def test_score_matches_length(epitopes, sequence, index):
    floating_peptides = map.float_peptides(table=epitopes, parent_seq=sequence, index=index)
    scores = [map.score_epitope_alignment(epi, sequence) for epi in floating_peptides]
    match_lengths = pd.Series([len(m) for s,m in scores])
    assert all(epitopes.length == match_lengths)


def test_floating_seqrecord(epitopes, index, sequence):
    floating_seqs = map.float_peptides(
        table=epitopes,
        parent_seq=sequence,
        index=index,
        id_col="seq"
    )
    seqrecord_ids = [r.id for r in floating_seqs]
    assert all(seqrecord_ids == epitopes.seq)

@pytest.fixture
def aligned_seq(sequence):
    return utils.random_gaps(
        seq=sequence,
        gap_prob=0.2,
        gap_size_interval=(1,5)
    )

@pytest.fixture
def aligned_epitopes(epitopes, aligned_seq, index):
    epitopes['newstart'] = epitopes.start.apply(
        map.align_coords,
        aligned_seq=aligned_seq,
        index=index
    )
    epitopes['newend'] = epitopes.end.apply(
        map.align_coords,
        aligned_seq=aligned_seq,
        index=index
    )
    return epitopes

def test_align_float(aligned_epitopes, aligned_seq, index, includeend):
    aligned_epitopes['float'] = map.float_peptides(
        table=aligned_epitopes,
        parent_seq=aligned_seq,
        start_col="newstart",
        index=index
    )
    trimmed_epitope_floats = aligned_epitopes.float.str.replace("^-*", "", regex=True)
    aligned_seq_epitope_slices = aligned_epitopes.apply(
        lambda x: aligned_seq[x.newstart-index:x.newend-index+includeend], axis=1
    )
    aligned_seq_epitope_slices = aligned_seq_epitope_slices.str.replace(
        "-*$", "", regex=True
    )
    assert all(trimmed_epitope_floats == aligned_seq_epitope_slices)
