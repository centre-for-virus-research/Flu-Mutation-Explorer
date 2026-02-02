# Flu-GDB Mutation Explorer

## Overview

The in**FLU**enza Virus **G**enome sequence **D**ata**B**ase (*Flu-GDB*) Mutation Explorer is a web resource which allows users to search for and visualise sequence variations in the genomes of influenza A viruses (IAVs).

Variation can be assessed using two tools:

### Tree

The Tree tool displays an interactive phylogenetic tree, based on nucleotide sequences, for each of the eight segments of the IAV genome. The tree can be coloured to illustrate features including subtype, host taxonomy and amino acid usage at specified positions, and publication-quality images can be downloaded.

### Adaptation Mutations

The Adaptation Mutations tool takes a user-supplied IAV protein sequence and identifies positions previously associated with changes in host species. For these positions, the amino acid usage in different host taxa is illustrated, in a graph that can be downloaded as a publication-quality image, and a table provides more detail on sites of adaptation mutations in the sequence, including links to supporting literature.

## Methods Summary

### Tree visualisation

The IAV phylogenetic trees are visualized using Taxonium, a tool for exploring large phylogenetic trees.

### Influenza virus Genome sequence Database

To assemble the sequence database, influenza A virus nucleotide sequences were retrieved from the NCBI Entrez databases using the influenza A virus taxonomic ID, using an in-house Python tool to access the E-utilities API. A BLAST search against a curated IAV reference set was performed to identify the closest reference for each sequence and validate the genome segment number.

Alignments for each genome segment were then obtained as follows:

1\. Nextalign was used to obtain sub-alignments of the sequences which shared the same closest reference in the IAV reference set

2\. the reference alignment of the IAV reference set was used to guide the insertion of gaps in individual sub-alignments, ensuring consistency with the corresponding reference sequences;

3\. sub-alignments were concatenated to obtain the full alignment per segment.

Several curation steps were applied, including removing non-IAV genomic sequences, eliminating sequences with non-significant BLAST hits, removing sequences shorter than a pre-established length threshold for each segment, and removing sequences that could not be aligned to a reference sequence with Nextalign. Nucleotide alignments were then trimmed to the coding sequences (CDS) for the IAV proteins, which were translated *in silico* to give protein sequence alignments.

For visualization, the IAV nucleotide sequence database was clustered, and phylogenetic trees were constructed. Clustering was performed on the nucleotide sequences using MMseqs2 with a 95% sequence identity threshold ). Representative sequences for each cluster and segment were aligned using MAFFT with default parameters. Maximum likelihood trees were inferred with IQ-TREE using the best-fit model determined by the Bayesian Information Criterion (BIC), and were midpoint-rooted for visualization purposes.

### Adaptation mutations

The database of mammalian adaptations is curated and maintained by Daniel Goldhill (Royal Veterinary College).

## Media

The virus image on the Home tab is attributed to Naina Nair and Ed Hutchinson, MRC-University of Glasgow Centre for Virus Research (CC-BY 2022).

## GitHub

Source code is available on <a href="https://github.com/wrightdw/flu-gdb" target="_blank" title="Flu-GDB GitHub repository">GitHub</a>.
Bug reports and feature requests may be submitted via the <a href="https://github.com/wrightdw/flu-gdb/issues" target="_blank" title="Flu-GDB issue tracker">issue tracker</a>.

## Software

Taxonium Sanderson, T. (2022). 
Taxonium, a web-based tool for exploring large phylogenetic trees. 
eLife, 11:e82392. 
DOI: <a href="https://doi.org/10.7554/eLife.82392" target="_blank">10.7554/eLife.82392</a>.

E-utilities Kans Jonathan. 
Entrez Direct: E-utilities on the Unix Command Line. 2013 Apr 23. 
In: Entrez Programming Utilities Help [Internet]. Bethesda (MD): National Center for Biotechnology Information (US); 2010-. 
<a href="https://www.ncbi.nlm.nih.gov/books/NBK179288/" target="_blank">https://www.ncbi.nlm.nih.gov/books/NBK179288/</a>.

MAFFT Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. 
Mol Biol Evol. 2013 Apr;30(4):772-80. 
DOI: <a href="https://doi.org/10.1093/molbev/mst010" target="_blank">10.1093/molbev/mst010</a>. 
Epub 2013 Jan 16. PMID: 23329690; PMCID: PMC3603318.

MMSEQ2 Steinegger M and Söding J. 
MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. 
Nature Biotechnology, 35, 1026–1028 (2017). 
DOI: <a href="https://doi.org/10.1038/nbt.3988" target="_blank">10.1038/nbt.3988</a>.

IQ-TREE Nguyen, L.-T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2015). 
IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies. 
Molecular Biology and Evolution, 32(1), 268–274. 
DOI: <a href="https://doi.org/10.1093/molbev/msu300" target="_blank">10.1093/molbev/msu300</a>. 

Nextalign Aksamentov, I., Roemer, C., Hodcroft, E. B., & Neher, R. A. (2021). 
Nextclade: clade assignment, mutation calling and quality control for viral genomes. 
Journal of Open Source Software, 6(67), 3773. 
DOI: <a href="https://doi.org/10.21105/joss.03773" target="_blank">10.21105/joss.03773</a>.

Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Xie Y, Allen J, McPherson J, Dipert A, Borges B (2025). 
*shiny: Web Application Framework for R*. R package version 1.11.0, 
<a href="https://shiny.posit.co/" target="_blank">Posit</a>.