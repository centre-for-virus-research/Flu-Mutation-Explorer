## Overview

The Flu-GDB project provides influenza phylogenies and associated metadata, plus a database of mammalian adaptation amino acid replacements. Phylogenies are created by clustering Genbank sequences for the 8 major genomic segments into representatives, which are then used to construct phylogenetic trees. The *Home* tab provides a summary table of the numbers of original sequences and clusters for each segment.

## Tree

The *Tree* tab provides an interactive phylogenic tree for each segment. Switch between segments by selecting *Segment*. Within a tree, the tips are named by strain and coloured by various metadata, which may be chosen by selecting *Colour by*. The tree may be searched for tip names or metadata by entering text into the *Search* box.

A query amino acid sequence may be placed into a tree. Paste into the *Query sequence* box and click *SUBMIT*. The closest segment will be determined and then the query sequence will be aligned with a database of represetative sequences, corresponding to the sequences in the tree. A table of hits is then displayed, sorted by closest match. Table rows may be clicked (ctrl/command-click for multiple) and corresponding sequences will be coloured in the tree.

Amino acids aligned to a position in a reference amino acid sequence may be visualised. Select a reference strain in *Subtype* and enter a position number in *Position* and click *SEARCH*. Amino acids aligned to this reference position will be coloured in the tree, using the ClustalX colouring scheme for amino acid chemistry.

## Adaptation Mutations

The *Adaptation Mutations* tab provides detection and analysis of mammalian adaptions in an amino acid sequence. Paste an amino acid sequence into the *Query sequence* box, and click *SUBMIT*. The closest segment will be determined and a table of any adaption mutations detected displayed. Click on a table row to visualise an interactive chart of the frequencies of amino acids aligned at that position across our dataset of cluster representative sequences.
