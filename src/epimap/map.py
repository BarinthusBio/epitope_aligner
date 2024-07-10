
from Bio import Seq
import pandas as pd

def _float_peptide(start, seq, index):
    if index == 1:
        start = start -1
    start = int(start)
    floating_peptide = "-" * start + seq
    return floating_peptide

def float_peptides(table, index, start_col="start", seq_col="seq"):
    floating_peptides = table.apply(
        lambda row: _float_peptide(row[start_col], row[seq_col], index=index),
        axis=1
    )
    # floating_peptides = floating_peptides.apply(Seq.Seq)
    floating_peptides = floating_peptides.tolist()
    return floating_peptides

def score_epitope_alignment(epitope, sequence, gap="-", toupper=True):
    """Proportion of aligned epitope positions that match the sequence

    Args:
        epitope (str): Aligned epitope sequence
        sequence (str): Sequence epitopes are aligned to
        gap (str, optional): Gap characters to ignore.
            Use a list of strings to ignore multiple gap types. Defaults to "-".
        toupper (bool, optional): Convert epitope and sequence to upper case before
        comparison. Defaults to True.

    Returns:
        tuple: (score, matches)

        - score (float): The proportion of non-gap epitope positions that match
        the sequence.
        - matches (list): List of booleans for matches of each non-gap position.
    """
    assert len(sequence) > len(epitope), f"The epitope ({epitope}) is longer than the sequence."
    if toupper:
        sequence = sequence.upper()
        epitope = epitope.upper()
    matches = []
    for seqaa, epiaa in zip(sequence, epitope):
        if epiaa not in gap:
            matches.append(seqaa == epiaa)
    score = sum(matches)/len(matches)
    return score, matches
