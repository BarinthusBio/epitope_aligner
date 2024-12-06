---
title: '`epitope_aligner`: placing epitopes in the context of sequence alignments'
tags:
    - Python
    - epitopes
    - vaccine design
    - sequence alignment
authors:
    - name: David A. Wells
      orcid: 0000-0002-4531-5968
affiliations:
    - name: Barinthus Biotherapeutics, UK
      index: 1
date: 13 November 2024
bibliography: paper.bib
---

# Summary
The adaptive immune system recognises specific regions of protein sequences, called epitopes. These epitopes are between 8-15 amino acids in length. When designing vaccines or immunotherapies, it is important to know the location of epitopes, both within the protein as a whole and relative to other epitopes and protein subdomains [@martinelli_silico_2022;@huang_antigenic_2012]. However, the precise location of an epitope can vary between related proteins because of insertions, deletions, and alternative splicing [REF]. These differences are common because protein antigens are often highly variable; in some cases this is because the pathogen is evolving to evade the immune system[@huang_antigenic_2012;@guo_evolutionary_2020]. `epitope_aligner` is a python package which harmonises epitope locations across related proteins, overcoming inconsistent epitope locations to identify epitope hotspots and help design vaccines and immunotherapies. 

# Statement of need
The location of epitopes in different proteins is not directly comparable unless the proteins have been aligned, but epitope locations are usually reported in unaligned sequences. This impedes our ability to generalise epitope information to multiple sequences. For example, viral proteins frequently contain insertions and deletions so position $i$ in one pathogen strain is not necessarily equivalent to position $i$ in a different strain. To define a common epitope coordinate system across multiple viral strains of a protein, we must first account for these differences. One way to do this is to align the epitopes to a common reference sequence(s) (for example with MAFFT [@katoh2012adding] or QuickAlign [@quickalign]); however, epitopes sequences are short and often align poorly which requires manual curation. `epitope_aligner` uses the aligned parent antigens to correctly convert epitope locations to a common reference frame, enabling analyse of epitopes from different but related sequences. This allows users to identify conserved epitopes and epitope hotspots, and to design effective vaccines and immunotherapies targeting diverse proteins. This tool has been used internally by Barinthus Biotherapeutics to design several antigens for candidate immunotherapies targeting infectious diseases, cancer and autoimmune diseases.

# Availability
`epitope_aligner` is available at [github.com/BarinthusBio/epitope_aligner](https://github.com/BarinthusBio/epitope_aligner) with installation instructions and detailed examples. The [quickstart](https://barinthusbio.github.io/epitope_aligner/epitope_aligner/examples/quickstart.html) demonstrates combining epitopes from different influenza virus strains. All of the available functions are described in detail in the [cookbook](https://barinthusbio.github.io/epitope_aligner/epitope_aligner/examples/cookbook.html) and [submodule API docs](https://barinthusbio.github.io/epitope_aligner/epitope_aligner.html).

# Acknowledgements
DW is employed by Barinthus Biotherapeutics PLC. I am grateful to Hugh Welles for comments on the manuscript.

# References
