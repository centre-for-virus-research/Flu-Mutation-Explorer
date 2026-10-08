# **Flu Mutation Explorer**

## **Overview**

The in**FLU**enza Virus Mutation Explorer is a web resource which allows users to search for and visualise sequence variations in the genomes of influenza A viruses (IAVs).

Variation can be assessed using two tools:

### **Tree**

The Tree tool displays an interactive phylogenetic tree, based on nucleotide or protein sequences, for each of the eight segments of the IAV genome. The tree can be coloured to illustrate features including subtype, host taxonomy and amino acid usage at specified positions, and publication-quality images can be downloaded.

### **Adaptation Mutations**

The Adaptation Mutations tool takes a nucleotide or protein protein sequence, and identifies positions previously associated with changes in host species. Segments 7 and 8 each encode two products, and M2 and NEP are listed separately from M1 and NS1 because their catalogued positions are numbered against those proteins. For these positions, the amino acid usage in different host taxa is illustrated, in a graph that can be downloaded as a publication-quality image, and a table provides more detail on sites of adaptation mutations in the sequence, including links to supporting literature.

### **Batch screening**

The Batch screening tab allows for the submission of up to 2000 sequences. The sequences have to be either all nucletoide or all protein sequences. The sequences are scanned against the database of adaptation mutations and for each site with a match to an adaptation mutation, the user can view the barplot frequencies of residues in different hosts or the distribution of residues at that site on the tree.

## **Methods Summary**

### **Tree visualisation**

The IAV phylogenetic trees are visualized using <a href="https://taxonium.org/">Taxonium</a>, a tool for exploring large phylogenetic trees.

### **Influenza virus Genome sequence Database**

To assemble the sequence database, influenza A virus nucleotide sequences were retrieved from the NCBI Entrez databases using the influenza A virus taxonomic ID, using an in-house Python tool to access the E-utilities API. A BLAST search against a curated IAV reference set was performed to identify the closest reference for each sequence and validate the genome segment number.

Alignments for each genome segment were then obtained as follows:

1\. Nextalign was used to obtain sub-alignments of the sequences which shared the same closest reference in the IAV reference set

2\. the reference alignment of the IAV reference set was used to guide the insertion of gaps in individual sub-alignments, ensuring consistency with the corresponding reference sequences;

3\. sub-alignments were concatenated to obtain the full alignment per segment.

#### **Subtype reference strains**

The IAV reference set is made up of the strains below, chosen to give at least one representative of every HA and every NA subtype, with additional strains kept where a subtype is more prevalent or more relevant for research (e.g. H1N1). It is not the full set of subtype combinations reported in GenBank. Strains were taken, in order of priority, from those proposed by <a href="https://doi.org/10.1371/journal.pone.0112302" target="_blank">Burke & Smith (2014)</a>, from the NCBI RefSeq influenza A reference genomes, and from a manual GenBank search where a subtype or a relevant subtype combination was still unrepresented. A dash indicates that no sequence for that segment is included in the reference set.

<div class="table-scroll">

| Strain | Subtype | PB2 (seg 1) | PB1 (seg 2) | PA (seg 3) | HA (seg 4) | NP (seg 5) | NA (seg 6) | M (seg 7) | NS (seg 8) | Source |
|---|---|---|---|---|---|---|---|---|---|---|
| A/California/07/2009 | H1N1 | NC_026438 | NC_026435 | NC_026437 | NC_026433 | NC_026436 | NC_026434 | NC_026431 | NC_026432 | RefSeq |
| A/New Caledonia/20/1999 | H1N1 | CY033629 | CY033628 | CY033627 | CY033622 | CY033625 | CY033624 | CY033623 | CY033626 | Burke & Smith |
| A/Puerto Rico/8/1934 | H1N1 | NC_002023 | NC_002021 | NC_002022 | NC_002017 | NC_002019 | NC_002018 | NC_002016 | NC_002020 | Burke & Smith |
| A/Korea/426/1968 | H2N2 | NC_007378 | NC_007375 | NC_007376 | NC_007374 | NC_007381 | NC_007382 | NC_007377 | NC_007380 | RefSeq |
| A/New York/392/2004 | H3N2 | NC_007373 | NC_007372 | NC_007371 | NC_007366 | NC_007369 | NC_007368 | NC_007367 | NC_007370 | RefSeq |
| A/Swine/Ontario/01911-1/99 | H4N6 | AF285892 | AF285891 | AF285890 | AF285885 | AF285888 | AF285887 | AF285886 | AF285889 | Burke & Smith |
| A/goose/Guangdong/1/1996 | H5N1 | NC_007357 | NC_007358 | NC_007359 | NC_007362 | NC_007360 | NC_007361 | NC_007363 | NC_007364 | RefSeq |
| A/chicken/Taiwan/0705/99 | H6N1 | DQ376876 | DQ376840 | DQ376803 | DQ376624 | DQ376732 | DQ376696 | DQ376659 | DQ376768 | Burke & Smith |
| A/turkey/Italy/8912/2002 | H7N3 | CY020612 | CY020611 | CY020610 | CY020605 | CY020608 | CY020607 | CY020606 | CY020609 | Manual&nbsp;† |
| A/Netherlands/219/03 | H7N7 | AY342413 | AY340083 | AY342418 | AY338459 | AY342425 | AY340079 | AY340089 | AY342422 | Burke & Smith |
| A/Shanghai/02/2013 | H7N9 | NC_026422 | NC_026423 | NC_026424 | NC_026425 | NC_026426 | NC_026429 | NC_026427 | NC_026428 | RefSeq |
| A/turkey/Ontario/6118/1968 | H8N4 | CY005831 | CY014662 | CY005830 | CY014659 | CY005829 | CY014660 | CY005828 | CY014661 | Burke & Smith |
| A/Hong Kong/1073/99 | H9N2 | NC_004910 | NC_004911 | NC_004912 | NC_004908 | NC_004905 | NC_004909 | NC_004907 | NC_004906 | RefSeq |
| A/Rousettus aegyptiacus/Egypt/381OP/2017 | H9N2 | PP273148 | PP273149 | PP273150 | PP273151 | PP273152 | PP273153 | PP273154 | PP273155 | Manual |
| A/mallard/Bavaria/3/2006 | H10N7 | – | – | – | FJ183474 | DQ792927 | FJ183475 | FJ743478 | – | Burke & Smith |
| A/duck/England/1/1956 | H11N6 | GU052209 | GU052208 | GU052207 | GU052203 | GU052205 | EU429795 | GU052204 | GU052206 | Burke & Smith |
| A/duck/Alberta/60/1976 | H12N5 | CY130085 | CY130084 | CY130083 | CY130078 | CY130081 | CY130080 | CY130079 | CY130082 | Burke & Smith |
| A/gull/Maryland/704/1977 | H13N6 | CY014701 | CY014700 | CY014699 | CY014694 | CY014697 | CY014696 | CY014695 | CY014698 | Burke & Smith |
| A/mallard/Astrakhan/263/1982 | H14N5 | CY130101 | CY130100 | CY130099 | CY130094 | CY130097 | CY130096 | CY130095 | CY130098 | Burke & Smith |
| A/duck/Australia/341/1983 | H15N8 | – | – | – | AB295613 | – | AB295614 | – | – | Burke & Smith |
| A/wedge-tailed shearwater/Western Australia/2576/1979 | H15N9 | CY005412 | CY005411 | CY005410 | CY006010 | CY005408 | CY005407 | CY005406 | CY005409 | Manual |
| A/black-headed gull/Turkmenistan/13/76 | H16N3 | – | – | – | EU293864 | – | – | – | – | Burke & Smith |
| A/little yellow-shouldered bat/Guatemala/060/2010 | H17N10 | CY103889 | CY103890 | CY103891 | CY103892 | CY103893 | CY103894 | CY103895 | CY103896 | Burke & Smith |
| A/flat-faced bat/Peru/033/2010 | H18N11 | CY125942 | CY125943 | CY125944 | CY125945 | CY125946 | CY125947 | CY125948 | CY125949 | Burke & Smith |
| A/lesser scaup/CA/1742/2013 | H19 | – | – | – | OR611723 | OR479766 | – | OR479890 | OR602902 | Manual |

</div>

† Found by manual search and used in place of A/turkey/Italy/220158/02 (Burke & Smith 2014), whose sequences are partial.

#### **Data curation in the pipeline**

Several curation steps were applied, including removing non-IAV genomic sequences, eliminating sequences with non-significant BLAST hits, removing sequences shorter than a pre-established length threshold for each segment, and removing sequences that could not be aligned to a reference sequence with Nextalign. Nucleotide alignments were then trimmed to the coding sequences (CDS) for the IAV proteins, which were translated *in silico* to give protein sequence alignments.

For visualization, the IAV nucleotide sequence database was clustered, and phylogenetic trees were constructed. Clustering was performed on the nucleotide sequences using MMseqs2 with a 95% sequence identity threshold. Representative sequences for each cluster and segment were aligned using MAFFT with default parameters. Maximum likelihood trees were inferred with IQ-TREE using the best-fit model determined by the Bayesian Information Criterion (BIC), and were midpoint-rooted for visualization purposes.

### **Adaptation mutations**

The database of mammalian adaptations is curated and maintained by Daniel Goldhill (Royal Veterinary College).

## **Media**

The virus image on the Home tab is from Naina Nair and Ed Hutchinson, MRC-University of Glasgow Centre for Virus Research (CC-BY 2022).

## **GitHub**

Source code is available on <a href="https://github.com/centre-for-virus-research/Flu-Mutation-Explorer" target="_blank" title="Flu Mutation Explorer GitHub repository">GitHub</a>.
Bug reports and feature requests may be submitted via the <a href="https://github.com/centre-for-virus-research/Flu-Mutation-Explorer/issues" target="_blank" title="Flu Mutation Explorer issue tracker">issue tracker</a>.

## **Software**

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

## **Funding**

We acknowledge funding from the Medical Research Council (MRC) and Department for Environment, Food and Rural Affairs (Defra, UK) as FluTrailMap-One Health (MR/Y03368X/1); and MRC funding to the MRC-University of Glasgow Centre for Virus Research: MC_UU_00034/5, MC_UU_00034/6, MC_UU_00034/1.

