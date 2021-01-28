# Introduction
The susceptibility of _Vanilla planifolia_ to certain fungal agents may be mediated by the interaction of the plant with its microbiomes (Ramirez-Mosqueda J. *et al.* 2018). One of these fungi, *Fusarium oxysporum*, causes stem and root rot disease, which is particularly devastating in vanilla plantations (Hernández-Hernández, J. *et al.* 2011; Loredo, 1990). It has been observed that the diversity of microbiomes in the cultivation soil and the presence of certain microorganisms in the microbiomes influence the development of this disease (Xiong *et al.* 2017). However, the role that the endophyte microbiome may be playing in the development of the disease has not been fully explored. Using a metabarcoding analysis of 16SrRNA (V3-V4 region, primers Bakt_341F (CCTACGGGNGGCWGCAG) and Bakt_805R (GACTACHVGGGTATCTAATCC)), and ITS sequence (ITS2 region, primers ITS3 (GCATCGATGAAGAACGCAGC) and ITS4 (TCCTCCGCTTATTGATATGC)). It was calculated the diversity of the endophytic microbiome (fungi and bacteria) in 3 groups if plants: a) symptomatic, b) asymptomatic and c) wild plants of _Vanilla planifolia._ Sequencing was performed using an Illumina MiSeq 2 × 300 bp platform. Analysing of the microbiomes of wild plants is a strategy to understand the role of the microbiome of wild plants in disease tolerance. In this repository, the procedure by which the bioinformatic analysis of the data was carried out is detailed. First step consist in quality control of the reads, OTU Clustering and Taxonomic Annotation. Second step consist in analysis of alpha diversity and beta diversity, with their respective visualizations and hypothesis tests. Third step consisted of visualizing taxonomic diversity and comparing the abundance of taxa between groups. 
## Summary report
Our results show a high species diversity, with 889 OTUs for ITS data, and for 16SrRNA data we have obtained 1670 OTUs in dataset. In ITS data, we obtained 8 phylum inside Fungal kingdom. In all sample groups, the principal phylum present was Ascomycota, with 607 OTU representing 68.27% of the total in dataset. In 16SrRNA data, we obtained 20 bacterial phylums, where 5 phylums represented >80% of total OTU in dataset. This phylums were Proteobacteria, with 573 OTU (34.3%), Actinobacteria with 277 OTU (16.6%), Bacteroidetes with 245 OTU (14.7%), Acidobacteria with 125 OTU (7.5%) and Planctomycetes with 124 OTU (7.4%). We find a higher diversity in root samples than stems in both of them amplicons. In Bacterial data, for root samples, symptomatic values are the highest values(R′=1009.00, H′=6.91), followed by asymptomatic samples (R′=915.00, H′=6.82) and wild samples (R′=494.50, H′=6.14). In stem samples, wild samples are the highest values (R′=442.00, H′=6.08). In Fungi data, symptomatic samples have higher values (R′=292.33, H′=5.64), followed by wild samples (R′=187.67, H′=4.98) and asymptomatic samples (R′=178.00, H′=5.13) in root samples. We don't observe the same effect in stem samples. Beta diversity values shows a important difference between sample groups (roots and stem; symptomatic, asymptomatic ans wild plants), accourding PERMANOVA test (p = 0.002 for fungi) and (p = 0.001 for bacteria). These results suggest an important difference in community structure relationed with groups, being the most important between stems and roots. In stems group, we observed a lower clustering into samples when compared with roots. 
# Description
## Organization
This repository is divided in 3 directories, and these are organized in subdirectories: 

```
###########################################*DIRECTORY STRUCTURE*######################################################
|
├── Data <- It contains the OTU table of 16SrRNA data and other of ITS data, and with results of diversity analysis.
|    ├── Alpha_diversity <- Output tables with values of Alpha divsersity.
|    ├── Beta_diversity <- Output tables with values of PERMANOVA test applied to dissimilarity matrices.
|    ├── Data_transformation <- Output OTU Tables of the scripts Data_transformation with relative abundance and presence-absence values. 
|
├── Figures <- Graphics and visualizations that are output files of scripts of Diversity analysis.
|      ├── Abundance_comparation <- Output plots of scripts of Abundance_comparison.
|      ├── Alpha_diversity <- Output plots with Alpha diversity values.
|      ├── Beta_diversity <- Output plots with Beta diversity values.
|      ├── Taxonomic_plots <- Output plots of Taxonomic anotation between sample groups. 
|
├── bin <- It contains scripts to run analysis 
     ├── 16S_analysis <- Pipeline of flow process using 16SrRNA amplicons.
     ├── ITS_analysis <- Pipeline of flow process using ITS amplicons.
```
## Process structure
### Quality filter, OTU clustering and assign taxonomy
The initial sequencing data processing step is filtering based on read quality scores and the presence of the expected primer and adapter sequences, and the removal of these non-biological sequences. The next step is OTU clustering. After, its necessary assign a taxonomic classification to OTUs. For all these steps, we use the [AMPtk pipeline](https://amptk.readthedocs.io/en/latest/index.html). Its scripts are conteined in the directory "1_Sequence_processing", and content four comands:  
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

We use the `R version 3.6.1 (2019-07-05)`

The libraries use in R you can install with this command: 
```
packages <-c("dplyr","ggplot2","taxa","miLineage","tidyr","vegan","metacoder")
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
Package versions in R that we use are next: 

- `dplyr` 1.02

- `ggplot2` 3.3.2

- `taxa` 0.3.4

- `miLineage` 2.1

- `tidyr` 1.1.2

- `vegan` 2.5.6

- `phyloseq` 1.30.0

- `metacoder` 0.3.4

# References 
- Hernández H., J . (2011). Paquete tecnológico vainilla (Vainilla planifolia Jackson) establecimiento y mantenimiento. Centro de Investigación Regional Golfo Centro. Campo Experimental Ixtacuaco. Instituto Nacional de Investigaciones Forestales, Agrícolas y Pecuarias.  
- Ramírez-Mosqueda, Marco Antonio; Iglesias Andreu, Lourdes Georgina; Noa Carrazana, Juan Carlos; Armas-Silva, Arturo Alonso (2018). Técnicas biotecnológicas para la obtención de genotipos de Vanilla planifolia Jacks. resistentes a Fusarium oxysporum f. sp. vanillae. Cuadernos de biodiversidad. 54. 9-14. 
- Xiong, Wu; Li, Rong; Ren, Yi; Liu, Chen; Zhao, Qingyun; Wu, Huasong; Jousset, Alexandre; Shen, Qirong. (2017) Distinct roles for soil fungal and bacterial communities associated with the suppression of vanilla Fusarium wilt disease. Soil Biology & Biochemistry 107. 198-207. 
- Loredo, S.X. (1990) Etiología de la necrosis del tallo de vainilla (Vanilla planifolia Andrews) en Papantla, Veracruz. Tesis de Maestría en Ciencias. Colegio de Postgraduados. Montecillos, México.
