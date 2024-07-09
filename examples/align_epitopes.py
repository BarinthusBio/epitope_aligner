import pandas as pd
from epimap import map

#Define an example sequence and set of epitopes
sequence = "abcdefghi"
epitopes = pd.DataFrame(
    {
        'start': [1.0,2,3,4,6],
        'end': [4,5,6,7,9],
        'seq': ['abc','bcd','cde','def','fgh']
    }
)

# Float the epitope sequences so they would be in the correct location in
# an alignment
map.float_peptides(epitopes)

# If your epitope table doesn't use "start" and "seq" columns
# you can set custom column names
epitopes2 = epitopes.copy()
epitopes2.columns = [c.upper() for c in epitopes2.columns]
map.float_peptides(epitopes2, start_col="START", seq_col="SEQ")

# Assess epitope alignment accuracy
sequence = "abcdefghi"
epitopes = pd.DataFrame(
    {
        'start': [1.0,2,3,4,6],
        'end': [4,5,6,7,9],
        'seq': ['abc','bcd','cde','def','fgh']
    }
)
floating_epitopes = map.float_peptides(epitopes)
epitope = floating_epitopes[-1]

map.score_epitope_alignment(epitope, sequence)
map.score_epitope_alignment(epitope, sequence, [" ", "-"])

[map.score_epitope_alignment(epi, sequence)[0] for epi in floating_epitopes]