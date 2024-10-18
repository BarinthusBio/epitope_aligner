"""
Investigate eptiopes from related sequences in a shared alignment

Map epitopes from related sequences onto a shared alignment,
quantify similarity to a specific sequence,
plot epitopes and overall epitope density in the alignment.
"""

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
from matplotlib import cm

import epitope_aligner.map
import epitope_aligner.stretch
import epitope_aligner.utils

def acc2seq(acc, seqs):
    seq = [r for r in seqs if acc in r.id][0]
    return seq

seqs = list(SeqIO.parse("examples/antigens_al.fa", "fasta"))

epitopes = pd.read_csv("examples/epitopes.csv")

epitopes['al_start'] = epitopes.apply(lambda x: epitope_aligner.map.align_coords(
    x.start,
    acc2seq(x.antigen_acc, seqs).seq.__str__(),
    index=1), axis=1)
epitopes['al_end'] = epitopes.apply(lambda x: epitope_aligner.map.align_coords(
    x.end,
    acc2seq(x.antigen_acc, seqs).seq.__str__(),
    index=1), axis=1)

epitopes = epitopes.sort_values("al_start")

epitopes['float'] = epitope_aligner.map.float_peptides(
    epitopes,
    start_col="al_start",
    seq_col="linear_sequence",
    index=1
)

reference_sequence = "P03452"
reference_sequence = acc2seq(reference_sequence, seqs).seq.__str__()
ref_scores = epitopes.float.apply(
    epitope_aligner.map.score_epitope_alignment,
    sequence = reference_sequence
)
ref_scores = ref_scores.apply(pd.Series, index=['score', 'bool'])
epitopes['reference_antigen_score'] = ref_scores.score

scores = epitopes.apply(
    lambda x: epitope_aligner.map.score_epitope_alignment(
        x.float,
        acc2seq(
            x.antigen_acc,
            seqs
        ).seq.__str__()
    ),
    axis=1
).apply(pd.Series, index=["score", "bool"])
epitopes['source_antigen_score'] = scores.score

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

fig, ax = plt.subplots(2)
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

ax[1].plot(positional_count)
plt.show()

