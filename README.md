# Introduction
The susceptibility of _Vanilla planifolia_ to certain fungal agents may be mediated by the interaction of the plant with its microbiomes. One of these fungi, Fusarium oxysporum, causes stem and root rot disease, which is particularly devastating in vanilla plantations. It has been observed that the diversity of microbiomes in the cultivation soil and the presence of certain microorganisms in the microbiomes influence the development of this disease. However, the role that the endophyte microbiome may be playing in the development of the disease has not been fully explored. Using a metabarcoding analysis of 16SrRNA and ITS sequence, it was calculated the diversity of the endophytic microbiome (fungi and bacteria) in 3 groups if plants: a) symptomatic, b) asymptomatic and c) wild plants of _Vanilla planifolia._ Sequencing was performed using an Illumina MiSeq 2 × 300 bp platform. Analysing of the microbiomes of wild plants is a strategy to understand the role of the microbiome of wild plants in disease tolerance. In this repository, the procedure by which the bioinformatic analysis of the data was carried out is detailed. First step consist in quality control of the reads, OTU Clustering and Taxonomic Annotation. Second step consist in analysis of alpha diversity and beta diversity, with their respective visualizations and hypothesis tests. Third step consisted of visualizing taxonomic diversity and comparing the abundance of taxa between groups. 
## Summary
Our results show a high species diversity, with 889 OTUs for ITS data, and for 16SrRNA data we have obtained 1670 OTUs in dataset. In ITS data, we obtained 8 phylum inside Fungal kingdom. In all sample groups, the principal phylum present was Ascomycota, with 607 OTU representing 68.27% of the total in dataset. In 16SrRNA data, we obtained 20 bacterial phylums, where 5 phylums represented >80% of total OTU in dataset. This phylums were Proteobacteria, with 573 OTU (34.3%), Actinobacteria with 277 OTU (16.6%), Bacteroidetes with 245 OTU (14.7%), Acidobacteria with 125 OTU (7.5%) and Planctomycetes with 124 OTU (7.4%). We find a higher diversity in root samples than stems in both of them amplicons. In Bacterial data, for root samples, symptomatic values are the highest values(R′=1009.00, H′=6.91), followed by asymptomatic samples (R′=915.00, H′=6.82) and wild samples (R′=494.50, H′=6.14). In stem samples, wild samples are the highest values (R′=442.00, H′=6.08). In Fungi data, symptomatic samples have higher values (R′=292.33, H′=5.64), followed by wild samples (R′=187.67, H′=4.98) and asymptomatic samples (R′=178.00, H′=5.13) in root samples. We don't observe the same effect in stem samples. Beta diversity values shows a important difference between sample groups (roots and stem; symptomatic, asymptomatic ans wild plants), accourding PERMANOVA test (p = 0.002 for fungi) and (p = 0.001 for bacteria). These results suggest an important difference in community structure relationed with groups, being the most important between stems and roots. In stems group, we observed a lower clustering into samples when compared with roots. 
# Description
## Organization
This repository is divided in 3 directories: 
1.	**bin.** This directory is subdivided in three directories: "1_Sequence_processing", "2_Normalizing_and_diversity_analysis", "3_Taxonomic_abundance_analysis". Inside each repositories, there are two scripts. one to 16S data and another to ITS data. The number of the directories correspond to order of the scripts in the process. 
2.	**Data.** It contains the OTU table of 16SrRNA data and other of ITS data. They are output file in the script "1_Sequence_processing", and they are input files in the another scripts. Also, here are contined a "mapping_file.txt", that describe the sample groups.  
3.	**Figures.** Graphics and visualizations that are output files of scripts of Diversity analysis.
## Process structure
### Quality filter, OTU clustering and assign taxonomy
The initial sequencing data processing step is filtering based on read quality scores and the presence of the expected primer and adapter sequences, and the removal of these non-biological sequences. The next step is OTU clustering. After, its necessary assign a taxonomic classification to OTUs. For all these steps, we use the AMPtk pipeline(https://amptk.readthedocs.io/en/latest/index.html). Its scripts are conteined in the directory "1_Sequence_processing", and content four comands:  
1. `amptk illumina` Demultiplexing Illumina data that has been delivered from sequencing center in a folder of FASTQ files, one set for each sample.
2. `amptk cluster` OTU picking using UPARSE algorithm and chimera filtering with reference based.  
3. `amptk filter` Removing index-bleed or sample cross-over from datasets. 
4. `amptk taxonomy` Assign taxonomy information to OTUs. 
### Estimating Alpha-diversity and Beta-diversity
In the scripts content in the directory "2_Normalizing_and_diversity_analysis", we calculate alpha-diversity and beta-diversity. 
Unequal sampling depth between samples makes necessary normalization. These normalization was made transform read counts into  % relative abundances. Next, was calculated to presence-absence of OTU in samples(binary tables) and we use this data for diversity analysis. 
The species richness was calculated in based number of OTU(observed), and we calculated Shannon’s Diversity Index. We aplicated Kruskall-Wallis test for evaluated distribution of these index. 
For Beta-diversity we use Bray-Curtis distance and Raup-Crick similarity index. For Visualization of this metrics, we use a NMDS plot. By hypothesis testing for this data, we use a PERMANOVA test.  
### Plotting taxonomic data and comparissons tax abund
For this visualizations, we use `metacoder` library. It's necessary to convert our data to dataset that they can use for this library. This dataset consist in an abundance matrix called `hmp_otus`, and a sample data table called `hmp_samples`. This library offers its own data transformation (for more information visit [here](https://grunwaldlab.github.io/metacoder_documentation/index.html)), which consists in filter OTUs with low count reads and divide each sample’s counts by the total number of counts observed for each sample, resulting in a proportion.
Taxonomic plot of the sample set its made with `heat_tree` command.
For comparisson of taxonomic abund between groups we use `compare_grops` funtion and this data its visualated with plot using `heat_tree_matrix` funtion.
## Installation 
For information on the requirements for installing AMPtk dependencies, please visit their webside [here](https://amptk.readthedocs.io/en/latest/index.html).

The libraries use in R you can install with this command: 
```
packages <-c("ape","dplyr","ggplot2","gplots","lme4","miLineage","phangorn","plotly","tidyr","vegan","VennDiagram","metacoder")
InsPack <- function(pack)
{
  if ((pack %in% installed.packages()) == FALSE) {
    install.packages(pack,repos ="http://cloud.r-project.org/")
  } 
}
lapply(packages, InsPack)
```
For installation of `phyloseq` library you can use the next form: 
```
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
```
