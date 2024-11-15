# Epitope aligner
Easily map epitope coordinates between sequences in an alignment,
regardless of which coordinate system you are using.
This lets you combine epitopes from different sources and calculate
things like epitope density in a set of proteins.

# Contents
- [Install](#install)
- [Quick start](#quick-start)
- [Examples](#examples)
- [Docs](#docs)
- [Dev](#dev)

# Install
Install directly from github using one of:
- `pip install git+https://github.com/Vaccitech/epitope_aligner.git`
- `pip install git+git@github.com:Vaccitech/epitope_aligner.git`

# Quickstart
The full quickstart example is [here](examples/quickstart.html) which analyses and plots the epitopes from different strains of the influenza virus.

In the current minimal example we'll:
- convert epitope coordinates to an aligned antigen
- float the epitope sequences to match it
- calculate the number of epitopes at each position in the antigen

For the inverse of these aligning and floating operations
see the [cookbook](examples/cookbook.html).

Import functions from `epitope_aligner` modules and pandas to create
an example dataset.
```
from epitope_aligner import map, stretch, utils
import pandas as pd
```

We'll define a short example antigen sequence, with an aligned
and unaligned version.
```
aligned_seq = "ABC---DEFGH-IJK--LM"
seq = aligned_seq.replace("-","")
```

We'll define some exmple epitopes with positions in the unaligned antigen sequence.
```
epitopes = pd.DataFrame({
        'start':  [2,      6,      9],
        'end':    [4,      9,      12],
        'seq':    ["BCD",  "FGHI", "IJKL"],
        "length": [3,       4,     4]
})
epitopes
```
```
#    start  end   seq  length
# 0      2    4   BCD       3
# 1      6    9  FGHI       4
# 2      9   12  IJKL       4
```

Let's calculate the start positions of these epitopes in the aligned
antigen sequence.
```
epitopes['newstart'] = map.align_coords(
    table = epitopes,
    aligned_parent_seq = aligned_seq,
    coord_col = "start",
    index = 1
)
epitopes
```
```
#    start  end   seq  length  newstart
# 0      2    4   BCD       3         2
# 1      6    9  FGHI       4         9
# 2      9   12  IJKL       4        13
```

Now we can "float" an epitope to line up with its antigen based on a start position and antigen sequence.
```
epitopes['float'] = map.float_epitopes(
    table=epitopes,
    parent_seq=aligned_seq,
    start_col="newstart",
    index=1,
)
epitopes
```
```
# Aligned antigen
ABC---DEFGH-IJK--LM

# Aligned epitopes
-BC---D
--------FGH-I
------------IJK--L
```

We can easily count the number of epitopes overlapping each position
by "stretching" them. For plotting, it is often helpful to add zeros
for positions with no epitopes.
```
stretched_epitopes = stretch.stretch(epitopes)
positional_count = stretched_epitopes.groupby("position").size()
positional_count = stretch.add_empty_positions(
    positional_count,
    parent_seq_length=len(seq),
    index=1,
    empty_value=0
)
positional_count
```
```
# position
# 1     0.0
# 2     1.0
# 3     1.0
# 4     1.0
# 5     0.0
# 6     1.0
# 7     1.0
# 8     1.0
# 9     2.0
# 10    1.0
# 11    1.0
# 12    1.0
# 13    0.0
# dtype: float64
```
Read the [cookbook](examples/cookbook.html) for tips on calculating more interesting measures than counts.

# Examples
Detailed examples are in `examples/` as jupyter notebooks,
and can be converted to html with `jupyter nbconvert --to html --output-dir docs/epitope_aligner/examples examples/*.ipynb --TagRemovePreprocessor.remove_cell_tags=hidden`.

# Docs
## Make docs
`pdoc -d google -o docs/ epitope_aligner`.

# Dev

## Set up
Create a virtual environment with `python3 -m venv .venv`.
Activate that environment with `. .venv/bin/activate`.
Install in editable mode with `pip install -e .`.
Deactivate it with `deactivate`.

## Nox
Linting, documentation, examples, and testing can all be run with
`nox`.
