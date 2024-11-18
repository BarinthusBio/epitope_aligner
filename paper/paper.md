---
title: 'epitope_aligner: placing epitopes in the context of sequence alignments'
tags:
    - Python
    - epitopes
    - vaccine design
    - sequence alignment
authors:
    - name: David A. Wells
      orcid: 0000-0002-4531-5968
      corresponding: true
affiliations:
    - name: Barinthus Biotherapeutics, UK
      index: 1
date: 13 November 2024
bibliography: paper.bib
---

# Summary
The immune system recognises specific regions of protein sequences, these regions are called epitopes. It is important to know these epitope locations when designing vaccines or immunotherapies. However, the location of equivalent epitopes can be different in related proteins. These differences are common because protein antigens are often highly variable; in some cases this is because the pathogen is evolving to evade the immune system. `epitope_aligner` is a python package designed to bring together epitope information from multiple related proteins to help design vaccines which provide broader protection.

# Statement of need
The location of epitopes in different proteins is not directly comparable unless the proteins have been aligned, but epitope locations are usually reported in unaligned sequences. This impedes our ability to generalise epitope information to multiple sequences. For example, viral proteins frequently contain insertions and deletions so position $i$ in one strain of influenza virus is not necessarily equivalent to position $i$ in a different strain. To identify common epitope positions across multiple viral strains, we must first account for these differences. One way to do this is to align the epitopes to a common reference sequence(s) (for example with MAFFT [@katoh2012adding] or QuickAlign [@quickalign]); however, epitopes sequences are short and often align poorly which requires unnecessary manual curation. `epitope_aligner` uses the aligned parent antigens to correctly convert epitope locations to a common reference frame, enabling analyse of epitopes from different but related sequences. Identifying epitope hotspots and conservation is key for designing broadly protective vaccines or immunotherapies targeting diverse proteins. This tool has been used internally by Barinthus Biotherapeutics to design several antigens to treat viral and autoimmune diseases.

# Availability
`epitope_aligner` is available at [github.com/BarinthusBio/epitope_aligner](https://github.com/BarinthusBio/epitope_aligner) with installation instructions and detailed examples. The [quickstart](https://barinthusbio.github.io/epitope_aligner/epitope_aligner/examples/quickstart.html) demonstrates combining epitopes from different influenza virus strains. All of the available functions are described in detail in the [cookbook](https://barinthusbio.github.io/epitope_aligner/epitope_aligner/examples/cookbook.html) and [submodule API docs](https://barinthusbio.github.io/epitope_aligner/epitope_aligner.html).

# References
