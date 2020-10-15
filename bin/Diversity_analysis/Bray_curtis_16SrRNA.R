packages <-c("ape","dplyr","ggplot2","gplots","lme4","miLineage","phangorn","plotly","tidyr","vegan","VennDiagram","DESeq2","phyloseq")
lib <- lapply(packages, require, character.only = TRUE)
MetaFile <- "mapping_file.txt"
MetaFile <- read.table(file=MetaFile,header=TRUE,sep="\t",comment.char = "",row.names = 1,check.names = F)
MetaFile <- MetaFile[!apply(is.na(MetaFile) | MetaFile=="",1,all),]
otu_file<-"OTUs_Table-norm-tax.tab"
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
View(physeq)
physeq
physeq_rel <- microbiome::transform(physeq, "compositional")
bx.ord_pcoa_bray <- ordinate(physeq_rel, "NMDS", "bray")
plot_scree(bx.ord_pcoa_bray) + theme_bw()
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable",
                            label = "Sites") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 12))
beta.ps3 <- beta.ps3 + theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 + scale_color_brewer(palette = "Dark2") + stat_ellipse()
beta.ps3
beta.ps3 + stat_ellipse()
beta.ps3 +  theme_bw() +
  stat_ellipse()
metadf.bx <- data.frame(sample_data(physeq_rel))
bray_ps.bxn <- phyloseq::distance(physeq = physeq_rel, method = "bray")
library(vegan)
adonis.test <- adonis(bray_ps.bxn ~ Variables, data = metadf.bx)
adonis.test
dist <- vegdist(t(abundances(ps4.rel)))
anova(betadisper(dist, metadf.bx$Variables))
library(vegan)
dist <- vegdist(t(abundances(physeq_rel)))
anova(betadisper(dist, metadf.bx$Variables))
View(metadf.bx)
distancia_organo=as.matrix(bray_ps.bxn)
heatmap.2(distancia_organo, key=T, trace="none", ColSideColors = c(rep("aquamarine",2),rep("maroon1",2),rep("blue",2),rep("red",3),rep("green",3),rep("darkgoldenrod3",2)), RowSideColors = c(rep("aquamarine",2),rep("maroon1",2),rep("blue",2),rep("red",3),rep("green",3),rep("darkgoldenrod3",2)))
