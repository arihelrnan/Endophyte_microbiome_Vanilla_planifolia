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
View(otu_table)
taxonomy=otu_table[,c(1,15)]
View(taxonomy)
taxonomy$V2.ITS=NULL
View(taxonomy)
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
plot_bar(top10P, x="Muestra", fill="Order") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +   ylab("Abundancia relativa") +
  ggtitle("ComposiciÃƒÂ³n a nivel Orden de las comunidades bacterianas")
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)
plot_bar(top10P, x="Muestra", fill="Order") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +   ylab("Abundancia relativa") +
  ggtitle("ComposiciÃƒÂ³n a nivel Orden de las comunidades bacterianas")
MetaFile$Muestra <- c(rep("V03_16S" , 1) , rep("V05_16S" , 1) , rep("V06_16S" , 1) , rep("V07_16S" , 1) , rep("V08_16S" , 1) , rep("V09_16S" , 1) , rep("V10_16S" , 1) , rep("V11_16S" , 1) , rep("V12_16S" , 1) , rep("V13_16S" , 1) , rep("V14_16S" , 1) , rep("V15_16S" , 1) , rep("V16_16S" , 1))
MetaFile$Muestra <- c(rep("V02_ITS" , 1) , rep("V03_ITS" , 1) , rep("V04_ITS" , 1) , rep("V06_ITS" , 1) , rep("V07_ITS" , 1) , rep("V08_ITS" , 1) , ("V09_ITS" , 1) , rep("V10_ITS" , 1) , rep("V11_ITS" , 1) , rep("V12_ITS" , 1) , rep("V13_ITS" , 1) , rep("V14_ITS" , 1) , rep("V15_ITS" , 1) , rep("V16_ITS" , 1))
MetaFile$Muestra <- c(rep("V02_ITS" , 1) , rep("V03_ITS" , 1) , rep("V04_ITS" , 1) , rep("V06_ITS" , 1) , rep("V07_ITS" , 1) , rep("V08_ITS" , 1) , rep("V09_ITS" , 1) , rep("V10_ITS" , 1) , rep("V11_ITS" , 1) , rep("V12_ITS" , 1) , rep("V13_ITS" , 1) , rep("V14_ITS" , 1) , rep("V15_ITS" , 1) , rep("V16_ITS" , 1))
samples = sample_data(MetaFile)
TAX = tax_table(tax_mat)
physeq <- phyloseq(OTU, TAX, samples)
plot_bar(top10P, x="Muestra", fill="Order") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +   ylab("Abundancia relativa") +
  ggtitle("ComposiciÃƒÂ³n a nivel Orden de las comunidades bacterianas")
plot_bar(top10P, x="Muestra", fill="phylum") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +   ylab("Abundancia relativa") +
  ggtitle("ComposiciÃ³n a nivel Phylum de las comunidades fÃºngicas")
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
View(otu_table)
taxonomy=otu_table[,c(1,15)]
View(taxonomy)
taxonomy$V2.ITS=NULL
View(taxonomy)
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
plot_bar(top10P, x="Muestra", fill="phylum") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +   ylab("Abundancia relativa") +
  ggtitle("ComposiciÃ³n a nivel Phylum de las comunidades fÃºngicas")
View(MetaFile)
MetaFile$Muestra <- c(rep("V02_ITS" , 1) , rep("V03_ITS" , 1) , rep("V04_ITS" , 1) , rep("V06_ITS" , 1) , rep("V07_ITS" , 1) , rep("V08_ITS" , 1) , rep("V09_ITS" , 1) , rep("V10_ITS" , 1) , rep("V11_ITS" , 1) , rep("V12_ITS" , 1) , rep("V13_ITS" , 1) , rep("V14_ITS" , 1) , rep("V15_ITS" , 1) , rep("V16_ITS" , 1))
samples = sample_data(MetaFile)
TAX = tax_table(tax_mat)
physeq <- phyloseq(OTU, TAX, samples)
View(physeq)
View(physeq)
top10P.names = sort(tapply(taxa_sums(physeq), tax_table(physeq)[, "Phylum"], sum), TRUE)[1:10]
top10P = subset_taxa(physeq, Phylum %in% names(top10P.names))
plot_bar(top10P, x="Muestra", fill="phylum") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +   ylab("Abundancia relativa") +
  ggtitle("ComposiciÃ³n a nivel Phylum de las comunidades fÃºngicas")
View(taxonomy)
View(taxonomy)
plot_bar(top10P, x="Muestra", fill="Phylum") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +   ylab("Abundancia relativa") +
  ggtitle("ComposiciÃ³n a nivel Phylum de las comunidades fÃºngicas")
top10P.names = sort(tapply(taxa_sums(physeq), tax_table(physeq)[, "Class"], sum), TRUE)[1:15]
top10P = subset_taxa(physeq, Class %in% names(top10P.names))
plot_bar(top10P, x="Muestra", fill="Class") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +   ylab("Abundancia relativa") +
  ggtitle("ComposiciÃ³n a nivel Clase de las comunidades fÃºngicas")
top10P.names = sort(tapply(taxa_sums(physeq), tax_table(physeq)[, "Order"], sum), TRUE)[1:15]
top10P = subset_taxa(physeq, Order %in% names(top10P.names))
plot_bar(top10P, x="Muestra", fill="Order") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +   ylab("Abundancia relativa") +
  ggtitle("ComposiciÃ³n a nivel Orden de las comunidades fÃºngicas")
Fusarium = subset_taxa(physeq, Genus == "g__Fusarium")
Fusarium = subset_taxa(physeq, Genus == "g__Fusarium")
Fusarium = subset_taxa(physeq, Genus == "g__Fusarium")
subset_taxa(physeq, Genus == "g__Fusarium")
subset_taxa(physeq, Genus = "g__Fusarium")
subset_taxa(physeq, Phylum == "p__Ascomycota")
subset_taxa(physeq, Phylum == "p__Ascomycota")
subset_taxa(physeq, Phylum == p__Ascomycota)
Fusarium = subset_taxa(physeq, Genus == " g__Fusarium")
plot_bar(Fusarium, x="Muestra", fill="Genus") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +   ylab("Abundancia relativa") +
  ggtitle("Abundancia del gÃ©nero Fusarium en las muestras")
