"""Utility functions such as random epitopes
"""
import numpy as np
import pandas as pd
from epimap import map
import matplotlib.pyplot as plt

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
    })
    epitopes.end = epitopes.apply(lambda row: row.end if row.end<len(sequence) else len(sequence), axis=1)
    epitopes['seq'] = epitopes.apply(lambda row:sequence[row.start:row.end+includeend], axis=1)
    epitopes['length'] = epitopes.seq.apply(len)
    epitopes.start = epitopes.start + index
    epitopes.end = epitopes.end + index
    return epitopes

def plot_line(dataframe, y, start_col="start", end_col="end", jitter=0, ax=None, **kwargs):
    if ax == None:
        ax = plt.gca()
    dataframe = dataframe.copy()
    jitter = np.random.uniform(0, jitter, dataframe.shape[0])
    dataframe['jitter'] = jitter
    for i,row in dataframe.iterrows():
        ax.plot(
            (row.start, row.end),
            (row[y]+row['jitter'], row[y]+row['jitter']),
            **kwargs
        )
    return ax

