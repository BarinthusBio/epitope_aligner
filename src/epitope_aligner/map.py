from typing import Literal
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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


def _float_peptide(start, seq, parent_seq, index):
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


def float_peptides(table, parent_seq, index, start_col="start", seq_col="seq", id_col=None):
    """Add gaps to sequences so they align to their parent

    Args:
        table (pd.DataFrame): Dataframe with sequences and their start
            position as columns
        parent_seq (str): The parent seq the peptide is derived from.
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
    floating_peptides = table.apply(
        lambda row: _float_peptide(
            start = row[start_col],
            seq = row[seq_col],
            parent_seq=parent_seq,
            index=index
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


def score_epitope_alignment(epitope, sequence, gap="-", toupper=True):
    """Proportion of aligned epitope positions that match the sequence

    Args:
        epitope (str): Aligned epitope sequence
        sequence (str): Sequence epitopes are aligned to
        gap (str, optional): Gap characters to ignore.
            Use a list of strings to ignore multiple gap types.
            Defaults to "-".
        toupper (bool, optional): Convert epitope and sequence to upper
            case before
        comparison. Defaults to True.

    Returns:
        tuple: (score, matches)

        - score (float): The proportion of non-gap epitope positions that
            match
        the sequence.
        - matches (list): List of booleans for matches of each non-gap position.
    """
    if not len(sequence) >= len(epitope):
        raise AssertionError(f"The epitope ({epitope}) is longer than the sequence.")
    if toupper:
        sequence = sequence.upper()
        epitope = epitope.upper()
    matches = []
    for seqaa, epiaa in zip(sequence, epitope):
        if epiaa not in gap:
            matches.append(seqaa == epiaa)
    score = sum(matches) / len(matches)
    return score, matches


def locate_peptide(seq, index, includeend):
    """Get start and end position of peptide in aligned sequence

    Args:
        seq (str): The aligned sequence containing the peptide
        index (int): Counting index, i.e. do positions start at 0 or 1?
        includeend (bool): Should the end position be included in the peptide?

    Returns:
        tuple: Start and end positions of the peptide in the provided sequence
    """
    # This regex pattern matches a string with non hypher at start and end and
    # anything inbetween (including hyphens)
    pattern = "[^-].*[^-]"
    seq = str(seq)
    start, end = re.search(pattern, seq).span()
    if index == 1:
        start += 1
        end += 1
    if includeend:
        end -= 1
    return start, end


def align_coords(coordinate: int, aligned_seq: str, index: Literal[0, 1], gap="-") -> int:
    """Convert coordinate from unaligned to aligned position

    The position in an unaligned sequence is converted to the
    equivalent position in an aligned version of the sequence.

    Note though that if you use it for the end of a slice (:x) it may
    run to the end of the next gaps.

    Args:
        coordinate (int): Position in an unaligned sequence.
        aligned_seq (str): The aligned version of the sequence.
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
        new_coord = [i for i, aa in enumerate(aligned_seq, index) if aa != gap][coordinate-index]
    except IndexError as e:
        if coordinate - index == len(aligned_seq.replace(gap, "")):
            new_coord = len(aligned_seq) + index
        else:
            raise Exception(f"{coordinate} not a valid position in ungapped aligned_seq") from e
    return new_coord


def unalign_coordinate(coordinate: int, aligned_seq: str, index: Literal[0, 1], gap="-") -> int:
    """Convert aligned coordinate to unaligned

    Convert a position in an anligned sequence to the equivalent in
    the unaligned sequence.

    Args:
        coordinate (int): the zero-indexed python coordinate
        aligned_seq (str): Aligned sequence the coordinate refers to.
        index (Literal[0, 1]): Counting index, i.e. do the position
            counts start at 0 or 1?
        gap (str, optional):  Character used for alignment gaps. Defaults to "-".

    Returns:
        int: The equivalent coordinate in an unaligned sequence.
    """
    gaps = aligned_seq[: (coordinate - index + 1)].count(gap)
    new_coord = coordinate - gaps
    return new_coord
