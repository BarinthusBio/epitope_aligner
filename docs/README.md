# Docs

The documentation is autogenerated using pdoc.

`docs/` also includes detailed examples generated
with:

`jupyter nbconvert --to html --output-dir docs/epitope_aligner/examples examples/*.ipynb --TagRemovePreprocessor.remove_cell_tags=hidden`

`pdoc -d google -o docs/ epitope_aligner`

The readme is included in the pdoc doc for `epitope_aligner/__init__.py`.
In fact it is included twice so that docs can be generated by calling
the pdoc command directly or by nox.

The `noxfile.py` generates the documentation, use this for local docs
generation. Although it does generate the docs when run as a github,
action it doesn't host the docs page and they are not kept. The docs
github action does both docs generation and hosting.
