# How to align epitope sequences to parent sequences
# and assess how well they match

# Different people use different numbering systems for their epitope
# positions, and we can accomidate that. IEDB uses a 1-index and
# the end position is part of the epitope, compared to python
# which uses 0-indexing and when slicing iterables the end position
# is not included.

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
map.float_peptides(epitopes, index=1)

# If your epitope positions use a zero index eg
epitopes0 = epitopes.copy()
epitopes0.start = epitopes0.start-1
# set the index as 0 when floating
map.float_peptides(epitopes0, index=0)

# If your epitope table doesn't use "start" and "seq" columns
# you can set custom column names
epitopes2 = epitopes.copy()
epitopes2.columns = [c.upper() for c in epitopes2.columns]
map.float_peptides(epitopes2, start_col="START", seq_col="SEQ", index=1)

# Assess epitope alignment accuracy
sequence = "abcdefghi"
epitopes = pd.DataFrame(
    {
        'start': [1.0,2,3,4,6],
        'end': [4,5,6,7,9],
        'seq': ['abc','bcd','cde','def','fgh']
    }
)
floating_epitopes = map.float_peptides(epitopes, index=1)
for epitope in floating_epitopes:
    print(sequence)
    print(epitope)

# Calculate the score of a single epitope
epitope = floating_epitopes[-1]
map.score_epitope_alignment(epitope, sequence)
# Returns a tuple, first element is the proportion of matches.
# Second element is a list of booleans for matches

# Speicify multiple gap characters to skip
map.score_epitope_alignment(epitope, sequence, gap=[" ", "-"])

# Use list comprehension to calculate the score for all epitopes
[map.score_epitope_alignment(epi, sequence) for epi in floating_epitopes]

# Low scores can indicate using the wrong index
floating_epitopes = map.float_peptides(epitopes0, index=1)
# selecting only the scores and not boolean lists
[map.score_epitope_alignment(epi, sequence)[0] for epi in floating_epitopes]
for epitope in floating_epitopes:
    print(sequence)
    print(epitope)

