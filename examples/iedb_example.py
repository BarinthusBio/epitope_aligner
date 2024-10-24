"""
Investigate eptiopes from related sequences in a shared alignment

Map epitopes from related sequences onto a shared alignment,
quantify similarity to a specific sequence,
plot epitopes and overall epitope density in the alignment.
"""

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
from matplotlib import cm
import re

import epitope_aligner.map
import epitope_aligner.stretch
import epitope_aligner.utils

# def acc2seq(acc, seqs):
#     seq = [r for r in seqs if acc in r.id][0]
#     return seq

seqs = list(SeqIO.parse("examples/antigens_al.fa", "fasta"))

def seqid2acc(id):
    """strip seq.id to get accession number"""
    acc = id
    pipes = acc.count("|")
    if pipes == 2:
        parts = acc.split("|")
        acc = parts[1]
        if acc == "":
            acc = parts[2]
    acc = re.sub("\\.\\d+", "", acc)
    return acc

# Create dictionary of sequences with accession numbers as keys.
aligned_parent_seqs = {}
for r in seqs:
    acc = seqid2acc(r.id)
    aligned_parent_seqs[acc] = r

epitopes = pd.read_csv("examples/epitopes.csv")

epitopes['al_start'] = epitope_aligner.map.align_coords(
    table = epitopes,
    aligned_parent_seq = aligned_parent_seqs,
    index=1,
    coord_col="start",
    parent_col="antigen_acc"
)
epitopes['al_end'] = epitope_aligner.map.align_coords(
    table = epitopes,
    aligned_parent_seq = aligned_parent_seqs,
    index=1,
    coord_col="end",
    parent_col="antigen_acc"
)

epitopes = epitopes.sort_values("al_start")

epitopes['float'] = epitope_aligner.map.float_epitopes(
    table = epitopes,
    parent_seq = aligned_parent_seqs,
    start_col = "al_start",
    seq_col = "linear_sequence",
    parent_col = "antigen_acc",
    index = 1
)

reference_sequence = "P03452"
reference_sequence = aligned_parent_seqs[reference_sequence]
epitopes[['reference_antigen_score','reference_antigen_matches']] = epitope_aligner.map.score_epitope_alignments(
    table = epitopes,
    parent_seq = reference_sequence.seq,
    seq_col = "float",
)

epitopes[['source_antigen_score','source_antigen_matches']] = epitope_aligner.map.score_epitope_alignments(
    table = epitopes,
    parent_seq = aligned_parent_seqs,
    seq_col = "float",
    parent_col = "antigen_acc"
)

# Save out epitopes and their source sequences as alignment
epitope_seqs = epitopes.apply(lambda x: f">{x.antigen_acc}\n{x.float}\n", axis=1)
with open("examples/antigen_epitope_alignment.fa", "w") as f:
    for r in seqs:
        f.write(f">{r.id}\n{r.seq}\n")
        antigen_mask = epitopes.antigen_acc.apply(lambda x: x in r.id)
        f.writelines(epitope_seqs[antigen_mask].tolist())

stretched_epitopes = epitope_aligner.stretch.stretch(
    epitopes,
    length_col="length",
    start_col="al_start",
    seq_col="linear_sequence"
)

positional_count = stretched_epitopes.groupby("position").size()
positional_count = epitope_aligner.stretch.add_empty_positions(
    series=positional_count,
    parent_seq_length=len(reference_sequence),
    index=1,
    empty_value=0
)

antigen_position_count = stretched_epitopes.groupby(["antigen_acc", "position"]).size()
grid = epitope_aligner.stretch.make_grid(
    grid_values = antigen_position_count,
    index=1,
    parent_seq_length = len(reference_sequence),
    empty_value=0
)

fig, ax = plt.subplots(3, sharex=True)
for i,acc in enumerate(epitopes.antigen_acc.unique()):
    mask = epitopes.antigen_acc == acc
    epitope_aligner.utils.plot_line(
        epitopes[mask],
        start_col="al_start",
        end_col="al_end",
        y="reference_antigen_score",
        jitter=0.1,
        c=cm.tab10(i),
        ax=ax[0]
    )
ax[0].set_ylabel("Reference alignment score")
ax[1].plot(positional_count)
ax[1].set_ylabel("Epitope count")
ax[2].matshow(grid, aspect="auto")
ax[2].set_title("Epitope count by antigen sequence")
ax[2].set_yticks(np.arange(len(grid.index)), labels=grid.index)
ax[2].set_xlabel("Position in alignment")
plt.show()
