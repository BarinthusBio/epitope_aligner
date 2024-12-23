# Epitope aligner
Easily map epitope coordinates between sequences in an alignment,
regardless of which coordinate system you are using.
This lets you combine epitopes from different sources and calculate
things like epitope density in a set of proteins.

`epitope_aligner` is a python package hosted on Github at
[BarinthusBio/epitope_aligner](https://github.com/BarinthusBio/epitope_aligner).

Full documentation at [barinthusbio.github.io/epitope_aligner](https://barinthusbio.github.io/epitope_aligner).

If you have any suggestions or problems, please open an [issue](https://github.com/BarinthusBio/epitope_aligner/issues).

# Contents
- [Install](#install)
- [Examples](#examples)
    - [Minimal](#minimal-examples)
    - [Quickstart](#quickstart)
    - [Cookbook](#cookbook)
- [Docs](#docs)
- [Dev](#dev)

# Install
Install from `pip` with:
```
pip install epitope-aligner
```
<!-- Install directly from github using one of:
- `pip install git+https://github.com/BarinthusBio/epitope_aligner.git`
- `pip install git+git@github.com:BarinthusBio/epitope_aligner.git` -->

# Examples
## Minimal examples
In the current minimal example we'll:
- convert epitope coordinates to an aligned antigen
- float the epitope sequences to match it
- calculate the number of epitopes at each position in the antigen

For the inverse of these aligning and floating operations
see the [cookbook](https://barinthusbio.github.io/epitope_aligner/epitope_aligner/examples/cookbook.html).

Import functions from `epitope_aligner` modules and pandas to create
an example dataset.
```python
from epitope_aligner import map, stretch, utils
import pandas as pd
```

We'll define a short example antigen sequence, with an aligned
and unaligned version.
```python
aligned_seq = "ABC---DEFGH-IJK--LM"
seq = aligned_seq.replace("-","")
```

We'll define some exmple epitopes with positions in the unaligned antigen sequence.
```python
epitopes = pd.DataFrame({
        'start':  [2,      6,      9],
        'end':    [4,      9,      12],
        'seq':    ["BCD",  "FGHI", "IJKL"],
        "length": [3,       4,     4]
})
epitopes
```
```python
#    start  end   seq  length
# 0      2    4   BCD       3
# 1      6    9  FGHI       4
# 2      9   12  IJKL       4
```

Let's calculate the start positions of these epitopes in the aligned
antigen sequence.
```python
epitopes['newstart'] = map.align_coords(
    table = epitopes,
    aligned_parent_seq = aligned_seq,
    coord_col = "start",
    index = 1
)
epitopes
```
```python
#    start  end   seq  length  newstart
# 0      2    4   BCD       3         2
# 1      6    9  FGHI       4         9
# 2      9   12  IJKL       4        13
```

Now we can "float" an epitope to line up with its antigen based on a start position and antigen sequence.
```python
epitopes['float'] = map.float_epitopes(
    table=epitopes,
    parent_seq=aligned_seq,
    start_col="newstart",
    index=1,
)
epitopes
```
```python
# Aligned antigen
# ABC---DEFGH-IJK--LM

# Aligned epitopes
# -BC---D
# --------FGH-I
# ------------IJK--L
```

We can easily count the number of epitopes overlapping each position
by "stretching" them. For plotting, it is often helpful to add zeros
for positions with no epitopes.
```python
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
```python
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
Read the [cookbook](https://barinthusbio.github.io/epitope_aligner/epitope_aligner/examples/cookbook.html#Stretch-epitopes) for tips on calculating more interesting measures than counts.

## Quickstart
A real world example is demonstrated in the [quickstart](https://barinthusbio.github.io/epitope_aligner/epitope_aligner/examples/quickstart.html) which analyses and plots the epitopes from different strains of the influenza virus.

## Cookbook
The [cookbook](https://barinthusbio.github.io/epitope_aligner/epitope_aligner/examples/cookbook.html) provides a detailed description and example of all functions.

# Docs
The [full documentation](https://barinthusbio.github.io/epitope_aligner/epitope_aligner.html) includes function APIs under the submodules:
- [map](https://barinthusbio.github.io/epitope_aligner/epitope_aligner/map.html)
- [stretch](https://barinthusbio.github.io/epitope_aligner/epitope_aligner/stretch.html)
- [utils](https://barinthusbio.github.io/epitope_aligner/epitope_aligner/utils.html)

# Dev
Details on testing, creating docs, and virtual envinments.

### Dev: Set up
Create a virtual environment with `python3 -m venv .venv`.
Activate that environment with `. .venv/bin/activate`.
Install in editable mode with `pip install -e .`.
Deactivate it with `deactivate`.

### Dev: Nox
Linting, bandit, documentation, examples, and testing can all be run with
`nox` based on `noxfile.py`. This is also run by github actions.

### Dev: Make docs
The full guide is `docs/README.md` but in short pdoc generates the
api documentation and renders the read me, jupyter notebook examples
are converted to html, and the complete docs are hosted at [barinthusbio.github.io/epitope_aligner/index.html](https://barinthusbio.github.io/epitope_aligner).

Generating the docs and hosting them is handled by the github actions, but
if you want to produce them locally just run `nox`.

### Dev: Publish to PyPI
Uploading requires the `build` and `twine` packages, 
`pip install --upgrade twine build`.

`python -m build` will create both the `--sdist` and `--wheel`.
`twine check dist/*` will check the package is ready for uploading.
`twine upload dist/*` will actually upload to pypi.
