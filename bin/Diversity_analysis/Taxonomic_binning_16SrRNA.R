packages <-c("ape","dplyr","ggplot2","gplots","lme4","miLineage","phangorn","plotly","tidyr","vegan","VennDiagram","DESeq2","phyloseq")
lib <- lapply(packages, require, character.only = TRUE)
MetaFile <- "mapping_file.txt"
MetaFile <- read.table(file=MetaFile,header=TRUE,sep="\t",comment.char = "",row.names = 1,check.names = F)
MetaFile <- MetaFile[!apply(is.na(MetaFile) | MetaFile=="",1,all),]
otu_file<-"OTUs_Table-norm-rel-tax.tab"
otu_table <-  read.table (otu_file,
                          check.names = FALSE,
                          header = TRUE,
                          dec = ".",
                          sep = "\t",
                          row.names = 1,
                          comment.char = "")
otu_table <- otu_table[!apply(is.na(otu_table) | otu_table=="",1,all),]
taxonomy=otu_table[,c(1,14)]
taxonomy$V3.16s=NULL
taxonomy = separate(taxonomy, taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";")
otu_mat=as.matrix(otu_file)
otu_mat=as.matrix(otu_table)
tax_mat=as.matrix(taxonomy)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
OTU = otu_table(otu_mat, taxa_are_rows = FALSE)
otu_table$taxonomy=NULL
otu_mat=as.matrix(otu_table)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
samples = sample_data(MetaFile)
TAX = tax_table(tax_mat)
physeq <- phyloseq(OTU, TAX, samples)
library(microbiome)
library(knitr)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(microbiomeutilities)
library(viridis)
library(tibble)
plot_composition(physeq, "Phylum", top = 11)
View(taxonomy)
lib <- lapply(packages, require, character.only = TRUE)
library(microbiome)
library(knitr)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(microbiomeutilities)
library(viridis)
library(tibble)
phylum_16S=physeq %>%
  tax_glom(taxrank = "Phylum") %>% psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)
plot_bar(phylum_16S, fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
View(taxonomy)
View(taxonomy)
View(phylum_16S)
View(phylum_16S)
plot_bar(phylum_16S, fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
top10P.names = sort(tapply(taxa_sums(physeq), tax_table(physeq)[, "Phylum"], sum), TRUE)[1:10]
top10P = subset_taxa(physeq, Phylum %in% names(top10P.names))
plot_bar(top10P, fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
library(taxa)
taxonomy <- filter_taxa(taxonomy, taxon_names != "")
taxonomy=as.data.frame(taxonomy)
library(taxa)
taxonomy <- filter_taxa(taxonomy, taxon_names != "")
physeq <- filter_taxa(physeq, taxon_names != "")
lib <- lapply(packages, require, character.only = TRUE)
sample.names(top10P)
library(forcats)
mutate(physeq = fct_relevel(physeq,
                            "V1_16S","V2_16S","V3_16S",
                            "V4_16S","V5_16S","V6_16S",
                            "V7_16S","V8_166S","V9_16S","V10_16S","V11_16S","V12_16S","V13_16S","V14_16S","V15_16S","V16_16S"))
View(samples)
View(samples)
plot_bar(top10P, x= fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
Muestras_16S <- c(rep("V1_16S" , 6) , rep("V2_16S" , 6) , rep("V3_16S" , 6) , rep("V4_16S" , 6), rep("V5_16S" , 6) , rep("V6_16S" , 6) , rep("V7_16S" , 6) , rep("V8_16S" , 6) , rep("V9_16S" , 6) , rep("V10_16S" , 6) , rep("V11_16S" , 6) , rep("V12_16S" , 6) , rep("V13_16S" , 6) , rep("V14_16S" , 6) , rep("V15_16S" , 6) , rep("V16_16S" , 6))
mutate(Muestras_16S = fct_relevel(Muestras_16S,
                                  "V1_16S","V2_16S","V3_16S",
                                  "V4_16S","V5_16S","V6_16S",
                                  "V7_16S","V8_16S","V9_16S","V10_16S","V11_16S","V12_16S","V13_16S","V14_16S","V15_16S","V16_16S"))
mutate(Muestras_16S = fct_relevel(Muestras_16S,
                                  "V1_16S" "V2_16S" "V3_16S"
                                  "V4_16S" "V5_16S" "V6_16S"
                                  "V7_16S" "V8_16S" "V9_16S" "V10_16S" "V11_16S" "V12_16S" "V13_16S" "V14_16S" "V15_16S" "V16_16S"))
mutate(Muestras_16S = fct_relevel(Muestras_16S,"V1_16S","V2_16S","V3_16S","V4_16S","V5_16S","V6_16S","V7_16S","V8_16S","V9_16S","V10_16S","V11_16S","V12_16S","V13_16S","V14_16S","V15_16S","V16_16S"))
View(samples)
View(MetaFile)
Muestras_16S <- c(rep("V3_16S" , 1) , rep("V5_16S" , 1) , rep("V6_16S" , 1) , rep("V7_16S" , 1) , rep("V8_16S" , 1) , rep("V9_16S" , 1) , rep("V10_16S" , 1) , rep("V11_16S" , 1) , rep("V12_16S" , 1) , rep("V13_16S" , 1) , rep("V14_16S" , 1) , rep("V15_16S" , 1) , rep("V16_16S" , 1))
mutate(Muestras_16S = fct_relevel(Muestras_16S,"V3_16S","V5_16S","V6_16S","V7_16S","V8_16S","V9_16S","V10_16S","V11_16S","V12_16S","V13_16S","V14_16S","V15_16S","V16_16S"))
library(dplyr)
mutate(Muestras_16S = fct_relevel(Muestras_16S,"V3_16S","V5_16S","V6_16S","V7_16S","V8_16S","V9_16S","V10_16S","V11_16S","V12_16S","V13_16S","V14_16S","V15_16S","V16_16S"))
mutate(Muestras_16S = fct_relevel(Muestras_16S,"V3_16S","V5_16S","V6_16S","V7_16S","V8_16S","V9_16S","V10_16S","V11_16S","V12_16S","V13_16S","V14_16S","V15_16S","V16_16S"))
packages <-c("ape","dplyr","ggplot2","gplots","lme4","miLineage","phangorn","plotly","tidyr","vegan","VennDiagram","DESeq2","phyloseq")
lib <- lapply(packages, require, character.only = TRUE)
MetaFile <- "mapping_file.txt"
MetaFile <- read.table(file=MetaFile,header=TRUE,sep="\t",comment.char = "",row.names = 1,check.names = F)
MetaFile <- MetaFile[!apply(is.na(MetaFile) | MetaFile=="",1,all),]
packages <-c("ape","dplyr","ggplot2","gplots","lme4","miLineage","phangorn","plotly","tidyr","vegan","VennDiagram","DESeq2","phyloseq")
lib <- lapply(packages, require, character.only = TRUE)
MetaFile <- "mapping_file.txt"
MetaFile <- read.table(file=MetaFile,header=TRUE,sep="\t",comment.char = "",check.names = F)
MetaFile <- MetaFile[!apply(is.na(MetaFile) | MetaFile=="",1,all),]
View(MetaFile)
otu_file<-"OTUs_Table-norm-rel-tax.tab"
otu_table <-  read.table (otu_file,
                          check.names = FALSE,
                          header = TRUE,
                          dec = ".",
                          sep = "\t",
                          row.names = 1,
                          comment.char = "")
otu_table <- otu_table[!apply(is.na(otu_table) | otu_table=="",1,all),]
taxonomy=otu_table[,c(1,14)]
taxonomy$V3.16s=NULL
taxonomy = separate(taxonomy, taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";")
otu_mat=as.matrix(otu_file)
otu_mat=as.matrix(otu_table)
tax_mat=as.matrix(taxonomy)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
OTU = otu_table(otu_mat, taxa_are_rows = FALSE)
otu_table$taxonomy=NULL
otu_mat=as.matrix(otu_table)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
samples = sample_data(MetaFile)
TAX = tax_table(tax_mat)
physeq <- phyloseq(OTU, TAX, samples)
packages <-c("ape","dplyr","ggplot2","gplots","lme4","miLineage","phangorn","plotly","tidyr","vegan","VennDiagram","DESeq2","phyloseq")
lib <- lapply(packages, require, character.only = TRUE)
MetaFile <- "mapping_file.txt"
MetaFile <- read.table(file=MetaFile,header=TRUE,sep="\t",comment.char = "",row.names = 1,check.names = F)
MetaFile <- MetaFile[!apply(is.na(MetaFile) | MetaFile=="",1,all),]
otu_file<-"OTUs_Table-norm-rel-tax.tab"
otu_table <-  read.table (otu_file,
                          check.names = FALSE,
                          header = TRUE,
                          dec = ".",
                          sep = "\t",
                          row.names = 1,
                          comment.char = "")
otu_table <- otu_table[!apply(is.na(otu_table) | otu_table=="",1,all),]
taxonomy=otu_table[,c(1,14)]
taxonomy$V3.16s=NULL
taxonomy = separate(taxonomy, taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";")
otu_mat=as.matrix(otu_file)
otu_mat=as.matrix(otu_table)
tax_mat=as.matrix(taxonomy)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
OTU = otu_table(otu_mat, taxa_are_rows = FALSE)
otu_table$taxonomy=NULL
otu_mat=as.matrix(otu_table)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
samples = sample_data(MetaFile)
TAX = tax_table(tax_mat)
physeq <- phyloseq(OTU, TAX, samples)
top10P.names = sort(tapply(taxa_sums(physeq), tax_table(physeq)[, "Phylum"], sum), TRUE)[1:10]
top10P = subset_taxa(physeq, Phylum %in% names(top10P.names))
plot_bar(top10P, x="Variables" fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
plot_bar(top10P, x="Variables", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
MetaFile$Muestra=("V1_16S","V2_16S","V3_16S",
                  "V4_16S","V5_16S","V6_16S",
                  "V7_16S","V8_166S","V9_16S","V10_16S","V11_16S","V12_16S","V13_16S","V14_16S","V15_16S","V16_16S"))
MetaFile$Muestra <- c(rep("V3_16S" , 1) , rep("V5_16S" , 1) , rep("V6_16S" , 1) , rep("V7_16S" , 1) , rep("V8_16S" , 1) , rep("V9_16S" , 1) , rep("V10_16S" , 1) , rep("V11_16S" , 1) , rep("V12_16S" , 1) , rep("V13_16S" , 1) , rep("V14_16S" , 1) , rep("V15_16S" , 1) , rep("V16_16S" , 1))
View(MetaFile)
View(MetaFile)
samples = sample_data(MetaFile)
TAX = tax_table(tax_mat)
physeq <- phyloseq(OTU, TAX, samples)
View(physeq)
View(physeq)
top10P.names = sort(tapply(taxa_sums(physeq), tax_table(physeq)[, "Phylum"], sum), TRUE)[1:10]
top10P = subset_taxa(physeq, Phylum %in% names(top10P.names))
plot_bar(top10P, x="Muestra", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
phyloseq%>%mutate(sam_data$Muestra = fct_relevel(sam_data$Muestra,
                                                 "V3_16S","V5_16S","V6_16S",
                                                 "V7_16S","V8_16S","V9_16S","V10_16S","V11_16S","V12_16S","V13_16S","V14_16S","V15_16S","V16_16S"))
)
)
)
phyloseq%>%mutate(sam_data$Muestra = fct_relevel(sam_data$Muestra,"V3_16S","V5_16S","V6_16S","V7_16S","V8_16S","V9_16S","V10_16S","V11_16S","V12_16S","V13_16S","V14_16S","V15_16S","V16_16S"))
MetaFile$Muestra <- c(rep("V03_16S" , 1) , rep("V05_16S" , 1) , rep("V06_16S" , 1) , rep("V07_16S" , 1) , rep("V08_16S" , 1) , rep("V09_16S" , 1) , rep("V10_16S" , 1) , rep("V11_16S" , 1) , rep("V12_16S" , 1) , rep("V13_16S" , 1) , rep("V14_16S" , 1) , rep("V15_16S" , 1) , rep("V16_16S" , 1))
samples = sample_data(MetaFile)
TAX = tax_table(tax_mat)
physeq <- phyloseq(OTU, TAX, samples)
top10P.names = sort(tapply(taxa_sums(physeq), tax_table(physeq)[, "Phylum"], sum), TRUE)[1:10]
top10P = subset_taxa(physeq, Phylum %in% names(top10P.names))
plot_bar(top10P, x="Muestra", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)
plot_bar(top10P, x="Muestra", fill="Phylum") +
  facet_grid(Station~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) + position="stack")
plot_bar(top10P, x="Muestra", fill="Phylum") +
  facet_grid(Station~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors)
plot_bar(top10P, x="Muestra", fill="Phylum") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors)
plot_bar(top10P, x="Muestra", fill="Phylum") +
  theme_ipsum + geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors)
plot_bar(top10P, x="Muestra", fill="Phylum") +
  theme_ipsum() + geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors)
plot_bar(top10P, x="Muestra", fill="Phylum") + geom_bar(stat = "identity") + theme_ipsum() +
  scale_fill_manual(values = phylum_colors)
plot_bar(top10P, x="Muestra", fill="Phylum") + geom_bar(stat = "identity") +   theme_minimal() +
  scale_fill_manual(values = phylum_colors)
plot_bar(top10P, x="Muestra", fill="Phylum") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +   ylab("Abundancia relativa")
plot_bar(top10P, x="Muestra", fill="Phylum") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +   ylab("Abundancia relativa") +
  ggtitle("Composición a nivel Phylum de las comunidades bacterianas")
View(taxonomy)
top10P.names = sort(tapply(taxa_sums(physeq), tax_table(physeq)[, "Class"], sum), TRUE)[1:10]
top10P.names = sort(tapply(taxa_sums(physeq), tax_table(physeq)[, "Class"], sum), TRUE)[1:10]
top10P = subset_taxa(physeq, Phylum %in% names(top10P.names))
top10P.names = sort(tapply(taxa_sums(physeq), tax_table(physeq)[, "Class"], sum), TRUE)[1:10]
top10P = subset_taxa(physeq, Class %in% names(top10P.names))
top10P.names = sort(tapply(taxa_sums(physeq), tax_table(physeq)[, "Class"], sum), TRUE)[1:15]
top10P = subset_taxa(physeq, Class %in% names(top10P.names))
plot_bar(top10P, x="Muestra", fill="Phylum") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +   ylab("Abundancia relativa") +
  ggtitle("Composición a nivel Clase de las comunidades bacterianas")
plot_bar(top10P, x="Muestra", fill="Class") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +   ylab("Abundancia relativa") +
  ggtitle("Composición a nivel Clase de las comunidades bacterianas")
top10P.names = sort(tapply(taxa_sums(physeq), tax_table(physeq)[, "Order"], sum), TRUE)[1:15]
top10P = subset_taxa(physeq, Order %in% names(top10P.names))
plot_bar(top10P, x="Muestra", fill="Order") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +   ylab("Abundancia relativa") +
  ggtitle("Composición a nivel Orden de las comunidades bacterianas")
