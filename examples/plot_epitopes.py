"""
Plot epitopes by some value
"""

import numpy as np
import pandas as pd
from epimap import stretch, utils
import itertools
import matplotlib.pyplot as plt

#Define an example sequence and set of epitopes
sequence = "abcdefghi"
epitopes = pd.DataFrame(
    {
        'start': [1.0,2,2,6],
        'end': [3,5,3,8],
        'seq': ['abc','bcde','bc','fgh'],
        'mhc_allele': ["x","x","y","z"]
    }
)
# Each epitope needs a length column to indicate how
# much it needs to be stretched
epitopes['length'] = epitopes.seq.apply(len)


# Stretch epitopes so each position they occupy is recorded.
stretched_epitopes = stretch.stretch(epitopes)
# Count how many epitopes overlap each position
positional_count = stretched_epitopes.groupby("position").size()
# add in missing positions
positional_count = stretch.add_empty_positions(
    positional_count,
    seq_length=len(sequence),
    index=1,
    empty_value=0
)
# Plot the number of epitopes at each position
plt.plot(positional_count)
plt.show()

utils.plot_line(epitopes, y="length", c="black")
plt.show()

# If you define some matplotlib axes you can pass those to make a multi panel
# figure. Set `sharex=True` to ensure positions line up between figures.
fig,ax = plt.subplots(3, sharex=True, layout="constrained")
ax[0].plot(positional_count)
utils.plot_line(epitopes, y="length", c="black", ax=ax[1])
plt.show()

# Plot an epitope grid
# Note that plt.matshow() discards the names of columns and rows
# But we can set these as shown in the next example
allele_position_count = stretched_epitopes.groupby(["mhc_allele", "position"]).size()
grid = stretch.make_grid(
    allele_position_count,
    index=1,
    seq_length=len(sequence),
    empty_value=0
)
plt.matshow(grid)
plt.show()


# It is also important to set the aspect ration to "auto" for matshow()
# We use .set_xticks to set the tick labels, it takes the tick number and
# label.

# Be careful when plotting grids with other figures. Even though we have
# relabelled the ticks they are actually still zero indexed and so will
# not align correctly with 1 indexed plots.
fig,ax = plt.subplots()
ax.matshow(grid)
ax.set_aspect("auto")
ax.set_xticks(np.arange(len(grid.columns)), grid.columns.get_level_values('position'))
ax.set_yticks(np.arange(len(grid.index)), labels=grid.index)
plt.show()

