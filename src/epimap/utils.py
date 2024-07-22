"""Utility functions such as random epitopes
"""
import numpy as np
import pandas as pd
from epimap import map

def random_seq(seq_length, AAs="ARNDCQEGHILKMFPSTWYV"):
    AAs = [a for a in AAs]
    seq = np.random.choice(AAs, seq_length)
    seq = "".join(seq)
    return seq

def random_epitopes(sequence, n, epitope_lengths, index, includeend):
    lengths = np.random.randint(
        epitope_lengths[0],
        epitope_lengths[1],
        n
    )
    starts = np.random.randint(0, len(sequence)-epitope_lengths[0], n)
    ends = starts + lengths
    epitopes = pd.DataFrame({
        "start":starts,
        "end":ends,
        "length":lengths
    })
    epitopes.end = epitopes.apply(lambda row: row.end if row.end<len(sequence) else len(sequence), axis=1)
    epitopes['seq'] = epitopes.apply(lambda row:sequence[row.start:row.end+includeend], axis=1)
    epitopes.start = epitopes.start + index
    epitopes.end = epitopes.end + index
    return epitopes

