
from Bio import Seq
import pandas as pd

def _float_peptide(start, seq):
    floating_peptide = "-" * start + seq
    return floating_peptide

def float_peptides(table, start_col="start", seq_col="seq"):
    floating_peptides = table.apply(
        lambda row: _float_peptide(row[start_col], row[seq_col]),
        axis=1
    )
    floating_peptides = floating_peptides.apply(Seq.Seq)
    floating_peptides = floating_peptides.tolist()
    return floating_peptides

