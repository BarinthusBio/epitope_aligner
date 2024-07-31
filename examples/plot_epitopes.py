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
        'seq': ['abc','bcde','bc','fgh']
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
positional_count = stretch.add_empty_positions(positional_count, len(sequence), 1, 0)
# Plot the number of epitopes at each position
plt.plot(positional_count)
plt.show()

utils.plot_line(epitopes, y="length", c="black")
plt.show()

# If you define some matplotlib axes you can pass those
fig,ax = plt.subplots(2)
ax[0].plot(positional_count)
utils.plot_line(epitopes, y="length", c="black", ax=ax[1])
plt.show()
