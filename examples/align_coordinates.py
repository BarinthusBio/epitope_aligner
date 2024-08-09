"""Align coordinates to a sequence with gaps in it

Epitope positions are usually given for the parent sequence
without any gaps in it. However, if you are looking at epitopes
from different sequences in an alignment it can be useful to convert
positions to that alignment.

`align_coordinate()` can convert epitope positions to positions in
an alignment if we have the aligned version of the parent sequence.
"""

import pandas as pd
from epimap import map

# The aligned version of our sequence
aligned_seq = "ABCD--FGHIJKL--OP--STUVWXYZ"
# The unaligned version
seq = aligned_seq.replace("-","")

# A table of epitopes in sequence, positions are relative to the
# unaligned version
epitopes = pd.DataFrame({
        'start':[1,     4,      7,      10,     9,      11,     13],
        'end':  [4,     8,      11,     14,     13,     15,     17],
        'seq':  ["BCD", "FGHI", "IJKL", "LOPS", "KLOP", "OPST", "STUV"]
})

# epitope coords are correct
assert all(epitopes.apply(lambda x: seq[x.start:x.end]==x.seq, axis=1))

# `align_coords()` takes an ungapped position, the sequence you would
# like it aligned in, and a counting index. It returns the equivalent position in that aligned sequence.
# Here position 4 in the unaligned sequence is equivalent to position 6
# in the aligned sequence
i = 4
j = map.align_coords(i, aligned_seq, index=0)
print(f"Position {i} in seq is {seq[i]}")
print(f"Position {j} in aligned_seq is {aligned_seq[j]}")

# A much faster way to use `align_coords()` is to `apply` it to several epitopes at once.
epitopes['newstart'] = epitopes.start.apply(map.align_coords, aligned_seq=aligned_seq, index=0)
epitopes['newend'] = epitopes.end.apply(map.align_coords, aligned_seq=aligned_seq, index=0)

# We can verify that the new positions are correct by slicing these
# new coordinates from the `aligned_seq`. Note that python slicing
# uses a zero index, and does not include the end index. Our example
# positions follow this form, but if you are using a different index
# and/or including the end index, you will have to account for this
# in your slicing.
epitopes['alseq'] = epitopes.apply(lambda x: aligned_seq[x.newstart:x.newend], axis=1)
epitopes[['seq','alseq']]

