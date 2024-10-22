from typing import Literal
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class ParentSeqSerialiser(object):
    """A class to get the correct parent sequence for each epitope

    The method `get_parent_seq()` takes a single argument `parent` and returns
    the parent sequence.

    How to get the parent sequence is based on the value of `parent_seq_object`.

    Attributes:
        parent_seq_object (dict | str): If `parent_seq_object` is a dictionary,
        parent sequence names should be keys and the values should be sequences
        as strings.

        If `parent_seq_object` is a string it is usually a single sequence used
        for all epitopes.
    
        However, if the string is "parent_seq_column" the argument `parent` passed
        to `get_parent_seq() is assumed to be the parent sequence and is returned
        as is.
    
    Methods:
        get_parent_seq(parent:str)->str:
            Returns the parent sequence.

            If `paret_seq_object` is a dictionary, `parent` should be a key in that
                dictionary.
            If `parent_seq_object` is the string "parent_seq_column", `parent` should
                be the parent sequence as is returned as is.
            If `parent` is anything else, the value of `parent` doesn't matter,
                `parent_seq_object` is returned regardless.
    """
    def __init__(self, parent_seq_object:dict|str):
        self.parent_seq_object = parent_seq_object
        self.get_parent_seq = self.get_serialiser()

    def serialise_single_seq(self, parent):
        return self.parent_seq_object
    
    def serialise_parent_seq_dict(self, parent):
        try:
            return self.parent_seq_object[parent]
        except KeyError as e:
            raise KeyError(f"{parent} not a key in parent_seq_object. Did you set `parent_col`?")

    def serialise_parent_seq_column(self, parent):
        return parent

    def get_serialiser(self):
        """Determines which function to use for the get_parent_seq() method"""
        if isinstance(self.parent_seq_object, dict):
            serialiser = self.serialise_parent_seq_dict
        elif self.parent_seq_object == "parent_seq_column":
            serialiser = self.serialise_parent_seq_column
        else:
            serialiser = self.serialise_single_seq
        return serialiser


def _align_float(start:int, seq:str, parent_seq:str, index:int, gap:str="-")->str:
    """Match spacing between peptide and parent seq

    If the peptide spans an insertion in the parent sequence,
    indicated by gaps, these gaps will be added to the floating
    peptide. E.g. 'LOPS' in 'ABCD--FGHIJKL--OP--STUVWXYZ'
    becomes 'L--OP--S'.

    Note, this does not add gaps the start of the peptide.
    Instead this function is called by _float_peptide() before
    adding gaps to the start of the peptide.

    Args:
        start (int): The start position of the peptide
        seq (str): The sequence to float
        parent_seq (str): The parent seq the peptide is derived from
        index (int): Counting index, i.e. do the position counts start at
            0 or 1?
        gap (str, optional): The gap character. Defaults to "-".

    Returns:
        str: The "floating" sequence
    """
    if isinstance(start, float):
        start = int(start)
    if index == 1:
        start = start - 1
    i = 0
    floating_peptide = []
    for aa in parent_seq[start:]:
        if aa == gap:
            floating_peptide.append(gap)
        else:
            floating_peptide.append(seq[i])
            i += 1
        if i >= len(seq):
            break
    floating_peptide = "".join(floating_peptide)
    return floating_peptide


def _float_epitope(start, seq, parent_seq, index):
    """Add gaps to sequence to it aligns with its parent.


    Args:
        start (int): The start position of the peptide
        seq (str): The sequence to float
        parent_seq (str): The parent seq the peptide is derived from
        index (int): Counting index, i.e. do the position counts start at
            0 or 1?

    Returns:
        str: The "floating" sequence
    """
    if not parent_seq is None:
        seq = _align_float(
            start=start,
            seq=seq,
            parent_seq=parent_seq,
            index=index
        )
    if index == 1:
        start = start - 1
    start = int(start)
    floating_peptide = "-" * start + seq
    return floating_peptide


def float_epitopes(table, parent_seq:str|dict, index, start_col="start", seq_col="seq", parent_col=None, id_col=None):
    """Add gaps to sequences so they align to their parent

    Args:
        table (pd.DataFrame): Dataframe with sequences and their start
            position as columns
        parent_seq (str|dict): The parent seq the peptide is derived from.
        index (int): Counting index, i.e. do the positions start at 0 or 1?
        start_col (str, optional): Name of the column with start positions.
            Defaults to "start".
        seq_col (str, optional): Name of column with sequences. Defaults to
            "seq".
        id_col (str, optional): If provided, this column is used as the id
            for sequence records. Defaults to None.

    Returns:
        list: List of floating sequences or (if `id_col` provided) a list
            of SeqRecords
    """
    pss = ParentSeqSerialiser(parent_seq_object=parent_seq)
    if parent_col is None:
        parent_col = start_col
    floating_peptides = table.apply(
        lambda row: _float_epitope(
            start = row[start_col],
            seq = row[seq_col],
            parent_seq = pss.get_parent_seq(row[parent_col]),
            index = index
        ), axis=1
    )
    # floating_peptides = floating_peptides.apply(Seq.Seq)
    floating_peptides = floating_peptides.tolist()
    if id_col is not None:
        ids = table[id_col]
        floating_peptides = [
            SeqRecord(Seq(seq), id=id) for id, seq in zip(ids, floating_peptides)
        ]
    return floating_peptides


def _score_epitope_alignment(seq, parent_seq, gap="-", toupper=True):
    """Proportion of aligned peptide positions that match the sequence

    Args:
        seq (str): Aligned peptide sequence
        parent_seq (str): Sequence peptides are aligned to
        gap (str, optional): Gap characters to ignore.
            Use a list of strings to ignore multiple gap types.
            Defaults to "-".
        toupper (bool, optional): Convert peptide and sequence to upper
            case before
        comparison. Defaults to True.

    Returns:
        tuple: (score, matches)

        - score (float): The proportion of non-gap peptide positions that
            match
        the sequence.
        - matches (list): List of booleans for matches of each non-gap position.
    """
    if not len(parent_seq) >= len(seq):
        raise AssertionError(f"The peptide ({seq}) is longer than the parent sequence.")
    if toupper:
        parent_seq = parent_seq.upper()
        seq = seq.upper()
    matches = []
    for parent_seqaa, epiaa in zip(parent_seq, seq):
        if epiaa not in gap:
            matches.append(parent_seqaa == epiaa)
    score = sum(matches) / len(matches)
    return score, matches


def score_epitope_alignments(table, parent_seq, seq_col, parent_col=None, gap="-", toupper=True):
    pss = ParentSeqSerialiser(parent_seq_object=parent_seq)
    if parent_col is None:
        parent_col=seq_col
    scores_matches = table.apply(
        lambda row: _score_epitope_alignment(
            seq=row[seq_col],
            parent_seq = pss.get_parent_seq(row[parent_col]),
            gap=gap,
            toupper=toupper
        ),
        axis=1,
        result_type="expand"
    )
    scores_matches.columns = ["score", "matches"]
    return scores_matches


def locate_epitope(aligned_seq, index, includeend):
    """Get start and end position of peptide in aligned sequence

    Args:
        aligned_seq (str): The aligned epitope sequence
        index (int): Counting index, i.e. do positions start at 0 or 1?
        includeend (bool): Should the end position be included in the peptide?

    Returns:
        tuple: Start and end positions of the peptide in the provided sequence
    """
    # This regex pattern matches a string with non hypher at start and end and
    # anything inbetween (including hyphens)
    pattern = "[^-].*[^-]"
    aligned_seq = str(aligned_seq)
    start, end = re.search(pattern, aligned_seq).span()
    if index == 1:
        start += 1
        end += 1
    if includeend:
        end -= 1
    return start, end


def _align_coord(coordinate: int, aligned_parent_seq: str, index: Literal[0, 1], gap="-") -> int:
    """Convert coordinate from unaligned to aligned position

    The position in an unaligned antigen sequence is converted to the
    equivalent position in an aligned version of the antigen sequence.

    Note though that if you use it for the end of a slice (:x) it may
    run to the end of the next gaps.

    Args:
        coordinate (int): Position in an unaligned sequence.
        aligned_parent_seq (str): The aligned version of the sequence.
        index (Literal[0,1]): Counting index, i.e. do the position counts start at
            0 or 1?
        gap (str, optional): Character used for alignment gaps. Defaults to "-".

    Raises:
        Exception: If the new coordinate raises an IndexError and is
            not the length of the sequence (accounting for index).

    Returns:
        int: The new coordinate
    """
    try:
        new_coord = [i for i, aa in enumerate(aligned_parent_seq, index) if aa != gap][coordinate-index]
    except IndexError as e:
        if coordinate - index == len(aligned_parent_seq.replace(gap, "")):
            new_coord = len(aligned_parent_seq) + index
        else:
            raise Exception(f"{coordinate} not a valid position in ungapped aligned_seq") from e
    return new_coord


def align_coords(table, aligned_parent_seq, coord_col, index, parent_col=None, gap="-"):
    pss = ParentSeqSerialiser(parent_seq_object=aligned_parent_seq)
    if parent_col is None:
        parent_col = coord_col
    new_coords = table.apply(
        lambda row: _align_coord(
            coordinate=row[coord_col],
            aligned_parent_seq=pss.get_parent_seq(row[parent_col]),
            index=index,
            gap=gap
        ),
        axis=1
    )
    return new_coords
    

def _unalign_coord(coordinate: int, aligned_parent_seq: str, index: Literal[0, 1], gap="-") -> int:
    """Convert aligned coordinate to unaligned

    Convert a position in an anligned sequence to the equivalent in
    the unaligned sequence.

    Args:
        coordinate (int): the zero-indexed python coordinate
        aligned_parent_seq (str): Aligned sequence the coordinate refers to.
        index (Literal[0, 1]): Counting index, i.e. do the position
            counts start at 0 or 1?
        gap (str, optional):  Character used for alignment gaps. Defaults to "-".

    Returns:
        int: The equivalent coordinate in an unaligned sequence.
    """
    gaps = aligned_parent_seq[: (coordinate - index + 1)].count(gap)
    new_coord = coordinate - gaps
    return new_coord

def unalign_coords(table, aligned_parent_seq, coord_col, index, parent_col=None, gap="-"):
    pss = ParentSeqSerialiser(parent_seq_object=aligned_parent_seq)
    if parent_col is None:
        parent_col = coord_col
    new_coords = table.apply(
        lambda row: _unalign_coord(
            row[coord_col],
            aligned_parent_seq=pss.get_parent_seq(row[parent_col]),
            index=index,
            gap=gap
        ),
        axis=1
    )
    return new_coords
