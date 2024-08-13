# epimap
Easily map epitope coordinates onto sequences and back,
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
- `pip install git+https://github.com/Vaccitech/epimap.git`
- `pip install git+git@github.com:Vaccitech/epimap.git`

# Quick start
Import functions `epimap` modules.
```
from epimap import map, stretch, utils
```
Let use a random sequence and random epitopes for examples.
```
sequence = utils.random_seq(10)
sequence
# KHPRILPWVF
```
```
epitopes = utils.random_epitopes(sequence, n=5, epitope_lengths=(3,6), index=0, includeend=False)
epitopes
#    start  end   seq  length
# 0      4    8  ILPW       4
# 1      1    5  HPRI       4
# 2      6    9   PWV       3
# 3      1    5  HPRI       4
# 4      3    7  RILP       4
```

You can turn an epitope start and end position into a sequence.
```
map.float_peptides(epitopes, index=0)
# ['----ILPW', '-HPRI', '------PWV', '-HPRI', '---RILP']
```
Given a floating epitope you can get the start and end positions.
```
map.locate_peptide("---ABC--", index=0, includeend=False)
# (3,6)
```
You can find the equivalent position in aligned an unaligned sequences
```
sequence = "ABCDE"
aligned_sequence = "A---BC-DE-"
map.align_coords(2, aligned_sequence, index=0)
# 5

map.unalign_coordinate(5, aligned_sequence, index=0)
# 2
```

To analyse epitopes overlapping a given position you can convert epitope
coordinates into a table with records for each position the epitope overlaps.
```
stretch.stretch(epitopes)
```
Read the examples to see how to calculate epitope density with this.

# Examples
Detailed examples are in `examples/` as jupyter notebooks,
and can be converted to html with `jupyter nbconvert --to html --output-dir docs/epimap/examples examples/*.ipynb --TagRemovePreprocessor.remove_cell_tags=hidden`.

# Docs
## Make docs
`pdoc -d google -o docs/ epimap`.

# Dev

## Set up
Create a virtual environment with `python3 -m venv .venv`.
Install in editable mode with `pip install -e .`.
Activat that environment with `. .venv/bin/activate`.
Deactivate it with `deactivate`.

## Nox
Linting, documentation, examples, and testing can all be run with
`nox`.
