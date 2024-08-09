
import numpy as np
import pandas as pd
from epimap import map, stretch
import itertools

#Define an example sequence and set of epitopes
sequence = "abcdefghi"
epitopes = pd.DataFrame(
    {
        'start': [1.0,2,2,6],
        'end': [3,4,4,8],
        'seq': ['abc','bcd','bcd','fgh'],
        'mhc_allele': ["x","x","y","z"]
    }
)
# Each epitope needs a length column to indicate how
# much it needs to be stretched
epitopes['length'] = epitopes.seq.apply(len)


# Stretch epitopes so each position they occupy is recorded.
# This means that instead of each row being an epitope,
# each row is a position occupied by a given epitope.
# Effectively duplicate each row for each residue in its epitope
# and adding a column indicating the position of that residue.
stretched_epitopes = stretch.stretch(epitopes)
# This allows us to assess all the epitopes at a given position
# Use groupby on the new position column and then apply whatever function you like

# For example size counts the number of epitopes overlapping a given position
positional_count = stretched_epitopes.groupby("position").size()
positional_count

# Or you can apply any function you like to different columns using agg
stretched_epitopes.groupby("position").agg(
    # Average start position of epitopes overlapping this position
    mean_start=('start', np.mean),
    # number of unique mhc alleles with epitopes at this position
    n_alleles = ('mhc_allele', lambda x: len(set(x)))
)

# Note that there are not records for positions without any epitopes
# But we can add them with add_empty_positions()
stretch.add_empty_positions(
    positional_count,
    seq_length=len(sequence),
    index=1,
    empty_value=0
)


# You can also group by position and some other value.
# This is useful for populating a grid with values.

# Number of epitopes for a given mhc allele at a given position
allele_position_count = stretched_epitopes.groupby(["mhc_allele", "position"]).size()
stretch.make_grid(
    allele_position_count,
    index=1,
    seq_length=len(sequence),
    empty_value=0
)

