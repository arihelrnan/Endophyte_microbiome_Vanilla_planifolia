# ITS flow procces

## Purpose

Here are contained scripts that are part of flow process to ITS sequences. They consist on the following: a) Processing of ITS amplicons since raw data: merged paired-end sequences, quality filter, assign OTU, remove chimeras, clustering sequence in OTU and its taxonomy assign b) Data transformation in relative abundance and transform OTU table in presence/absence of OTU through binary table, and use binary data analysis c) Calculate alpha diversity using observed OTU values and Shannon´s Index and explain results through a table with this values and boxplots plots and statical analysis d) Calculate beta diversity using Raup-Crick Dissimilarity and explain results in an NMDS plot and use PERMANOVA test to evaluated distance f) Obtain a heat tree plot with taxonomic description across organ and plant state g) Using data in proportions, make a heat tree plot with significant differences using Wilcoxon test for multiple comparisons.       

## Backgraund 

The number in the file names correspond to order of the scripts in the process. 

### `1_ITS_paired.sh`

The firts file is the `1_ITS_paired.sh`. This file consist in the four scripts: 

1. `amptk illumina` Take a folder with Forward and reverse sequence hat is already de-multiplexed and processes it for clustering using AMPtk. Raw data sequences are not contained in this repository. In this step, reads Pared-end are merged using USEARCH algorithm. We use next arguments: 

```
-i Folder name that contained raw data sequences in Paired-end. This data is not contained in this repository, but I named the forder as `Illumina_ITS_data/` 
-o Output folder name with merged reads. Was name as `ITS_merged/`, but is not contained in this repository. 
-f Name of forwrad primer. We use primer ITS3 
-r Name of reverse primer. We use primer ITS4
-l Length to trim/pad reads. We trim our reads in 300 pb
-min_len Minimum length read to keep. The minimum length that we keep is 150 pb
```
2. `amptk cluster` Cluster reads in OTU using UPARSE algorithm. Also, make a chimera filtering. We use next arguments: 

```
-i
-o
--uchime_ref
-m
```

3. `amptk filter`

```

```

4. `amptk taxonomy`

```

```
