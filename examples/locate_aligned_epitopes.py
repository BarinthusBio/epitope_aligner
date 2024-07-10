"""Locate aligned epitopes

Given epitopes aligned to sequences, what are the positions of
the epitopes in that alignment"""

import re
import pandas as pd
from epimap import map
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord

seq1 = SeqRecord(Seq("abcdefghi"), description="seq1")
seq2 = SeqRecord(Seq("ABCDEFGHI"), description="seq2")
epi1 = SeqRecord(Seq("abc"), description="epi1")
epi2 = SeqRecord(Seq("--cdef"), description="epi2")
epi3 = SeqRecord(Seq("-----fghi"), description="epi3")
epi4 = SeqRecord(Seq("--cdef---"), description="epi4")
epi5 = SeqRecord(Seq("--cd-f---"), description="epi5")

al = [
    seq1, seq2,
    epi1, epi2, epi3, epi4, epi5
]
# Given aligned sequences you can get their non-gap start and end
al_df = pd.DataFrame({
    'description': [r.description for r in al],
    'peptide': [str(r.seq.replace("-", "")) for r in al],
    'seq': [str(r.seq) for r in al]
})

pos = al_df.seq.apply(map.locate_peptide, index=0, includeend=False)
al_df[['start','end']] = pd.DataFrame(pos.tolist())

# Python uses zero indexes and slicing that doesn't include
# the end position. However, you can use the IEDB
# format if you'd like
for r in al:
    map.locate_peptide(r.seq, index=1, includeend=True)

