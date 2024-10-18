"""Utility functions such as random epitopes
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def random_seq(seq_length, AAs="ARNDCQEGHILKMFPSTWYV"):
    """Generate a random sequence, useful for tests.

    Args:
        seq_length (int): Length of sequence to generate.
        AAs (str, optional): Possible characters for seq.
            Defaults to "ARNDCQEGHILKMFPSTWYV".

    Returns:
        str: Random sequence
    """
    AAs = [a for a in AAs]
    seq = np.random.choice(AAs, seq_length)
    seq = "".join(seq)
    return seq


def random_epitopes(sequence, n, epitope_lengths, index, includeend):
    """Generate random epitopes for a sequence, useful for tests.

    Args:
        sequence (str): The parent sequence
        n (int): Number of epitopes to generate.
        epitope_lengths (tuple): Tuple of min and max epitope length.
        index (int): Counting index, Do positions start at 0 or 1?
        includeend (bool): Is the end position included in the epitope?

    Returns:
        pd.DataFrame: Table of randomly generated epitopes.
    """
    lengths = np.random.randint(epitope_lengths[0], epitope_lengths[1], n)
    starts = np.random.randint(0, len(sequence) - epitope_lengths[0], n)
    ends = starts + lengths
    epitopes = pd.DataFrame(
        {
            "start": starts,
            "end": ends,
        }
    )
    epitopes.end = epitopes.apply(
        lambda row: row.end if row.end < len(sequence) else len(sequence)-1, axis=1
    )
    epitopes["seq"] = epitopes.apply(
        lambda row: sequence[row.start:row.end + includeend], axis=1
    )
    epitopes["length"] = epitopes.seq.apply(len)
    epitopes.start = epitopes.start + index
    epitopes.end = epitopes.end + index
    return epitopes


def random_gaps(seq: str, gap_prob: float, gap_size_interval: tuple[int, int]) -> str:
    """Introduce gaps to a sequence to simulate alignment

    Args:
        seq (str): Sequence to introduce gaps to
        gap_prob (float): Probability of introducing a gap at each
            position
        gap_size_interval (tuple[int, int]): Minimum (inclusive) and maximum (exclusive)
            size of gaps introduced.

    Returns:
        str: The sequence with gaps.
    """
    aligned_seq = []
    gap_prob = 0.1
    for a in seq:
        if gap_prob > np.random.random():
            gap_size = np.random.randint(gap_size_interval[0], gap_size_interval[1])
            aligned_seq.append("-" * gap_size)
        aligned_seq.append(a)
    aligned_seq = "".join(aligned_seq)
    return aligned_seq


def plot_line(
    dataframe, y, start_col="start", end_col="end", jitter=0, ax=None, **kwargs
):
    """Plot epitope values as a line showing the epitope's location.

    Args:
        dataframe (pd.DataFrame): Table of epitopes
        y (str): Column of values to plot
        start_col (str, optional): Column of epitope start positions. Defaults to "start".
        end_col (str, optional): Column of epitope end positions. Defaults to "end".
        jitter (int, optional): How much to jitter y values to avoid overlapping.
            Defaults to 0.
        ax (matplotlib.Axes, optional): Axes to plot to,
            if None uses `matplotlib.pyplot.gca()`. Defaults to None.

    Returns:
        matplotlib.Axes: Axes with the generated plot
    """
    if ax is None:
        ax = plt.gca()
    dataframe = dataframe.copy()
    jitter = np.random.uniform(0, jitter, dataframe.shape[0])
    dataframe["jitter"] = jitter
    for i, row in dataframe.iterrows():
        ax.plot(
            (row[start_col], row[end_col]),
            (row[y] + row["jitter"], row[y] + row["jitter"]),
            **kwargs,
        )
    return ax