Using a metabarcoding analysis of 16SrRNA and ITS sequence, it was calculated the diversity of the endophytic microbiome (fungi and bacteria) in 3 groups if plants: a) symptomatic, b) asymptomatic and c) wild plants of Vanilla planifolia, and make their characterization taxonomic. We sequenced the bacterial 16S rRNA region V3-V4 with primers Bakt_341F (CCTACGGGNGGCWGCAG) and Bakt_805R (GACTACHVGGGTATCTAATCC), and the fungal region ITS2 with primers 3F (GCATCGATGAAGAACGCAGC) and 4R (TCCTCCGCTTATTGATATGC). Sequencing was performed using an Illumina MiSeq 2 × 300 bp platform. 
# Organization
This repository is divided in 4 directories: 
1.	bin. This directory has scripts necessary to bioinformatics processing and diversity analysis. It is divided in two subdirectories: a) Sequence processing, with necessary scripts in Shell for merge pairs sequences and their quality filter, identify chimeras sequence, OTU assignment and their taxonomy assignment and make OTU table; b) Diversity analysis, with scripts in R for normalization of data set, calculate alpha diversity and beta diversity and their visualizations and statistical analysis, and graphics of binning taxonomic. 
2.	Objets. This directory has reference data base and parameters files necessary for sequence processing scripts. 
3.	Meta. It contains the OTU table of 16SrRNA data and other of ITS data.  It also contains the meta file with information of each sample and groups. 
4.	Figures. Graphics and visualizations that are files out of scripts in Diversity analysis directory.
