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
physeq
View(samples)
physeq_rel <- microbiome::transform(physeq, "compositional")
bx.ord_pcoa_bray <- ordinate(physeq_rel, "NMDS", "bray")
plot_scree(bx.ord_pcoa_bray) + theme_bw()
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 12))
View(samples)
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 12))
beta.ps3 <- beta.ps3 + theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 + scale_color_brewer(palette = "Dark2") + stat_ellipse()
metadf.bx <- data.frame(sample_data(physeq_rel))
bray_ps.bxn <- phyloseq::distance(physeq = physeq_rel, method = "bray")
library(vegan)
adonis.test <- adonis(bray_ps.bxn ~ Variable, data = metadf.bx)
adonis.test
physeq_rel <- microbiome::transform(physeq, "compositional")
bx.ord_pcoa_bray <- ordinate(physeq_rel, "NMDS", "bray")
plot_scree(bx.ord_pcoa_bray) + theme_bw()
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 12))
View(samples)
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 20))
beta.ps3 <- beta.ps3 + theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 + scale_color_brewer(palette = "Dark2") + stat_ellipse()
physeq_rel <- microbiome::transform(physeq, "compositional")
bx.ord_pcoa_bray <- ordinate(physeq_rel, "NMDS", "bray")
plot_scree(bx.ord_pcoa_bray) + theme_bw()
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 12))
View(samples)
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 12))
beta.ps3 <- beta.ps3 + theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 + scale_color_brewer(palette = "Dark1")
physeq_rel <- microbiome::transform(physeq, "compositional")
bx.ord_pcoa_bray <- ordinate(physeq_rel, "NMDS", "bray")
plot_scree(bx.ord_pcoa_bray) + theme_bw()
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 12))
View(samples)
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 12))
beta.ps3 <- beta.ps3 + theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 + scale_color_brewer(palette = "Dark2")
View(bx.ord_pcoa_bray)
View(bx.ord_pcoa_bray)
physeq_rel <- microbiome::transform(physeq, "compositional")
bx.ord_pcoa_bray <- ordinate(physeq_rel, "NMDS", "bray")
plot_scree(bx.ord_pcoa_bray) + theme_bw()
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 12))
View(samples)
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 12))
beta.ps3 <- beta.ps3 + theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 + scale_color_brewer(palette = "Dark2")
beta.ps3 + StatEllipse
beta.ps3 + stat_ellipse()
plot_ordination(physeq_rel, bx.ord_pcoa_bray, type="sites", color="Variable") +
  theme_bw() +
  stat_ellipse()
plot_ordination(physeq_rel, bx.ord_pcoa_bray, type="sites", color="Variable") +
  theme_bw() +
  stat_ellipse()
plot_ordination(physeq_rel, bx.ord_pcoa_bray, type="sites", color="Estado") +
  theme_bw() +
  stat_ellipse()
plot_ordination(physeq_rel, bx.ord_pcoa_bray, type="sites", color="Variable") +
  theme_bw() +
  stat_ellipse()
lib <- lapply(packages, require, character.only = TRUE)
plot_ordination(physeq_rel, bx.ord_pcoa_bray, type="sites", color="Variable") +
  theme_bw() +
  stat_ellipse()
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="SampleType",
                            label = "sites") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 10))
View(samples)
beta.ps3 <- beta.ps3 + theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 + scale_color_brewer(palette = "Dark2") + stat_ellipse()
beta.ps3
beta.ps3 + stat_ellipse()
beta.ps3 +  theme_bw() +
  stat_ellipse()
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable",
                            label = "sites") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 10))
View(samples)
beta.ps3 <- beta.ps3 + theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 + scale_color_brewer(palette = "Dark2") + stat_ellipse()
beta.ps3
beta.ps3 + stat_ellipse()
beta.ps3 +  theme_bw() +
  stat_ellipse()
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable",
                            label = "sites") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 10))
View(samples)
beta.ps3 <- beta.ps3 + theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 + scale_color_brewer(palette = "Dark2") + stat_ellipse()
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable",
                            label = "sites") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 14))
View(samples)
beta.ps3 <- beta.ps3 + theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 + scale_color_brewer(palette = "Dark2") + stat_ellipse()
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable",
                            label = "sites") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 20))
View(samples)
beta.ps3 <- beta.ps3 + theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 + scale_color_brewer(palette = "Dark2") + stat_ellipse()
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable",
                            label = "sites") +
  geom_point(size= 2) +
  theme(plot.title = element_text(hjust = 0, size = 20))
View(samples)
beta.ps3 <- beta.ps3 + theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 + scale_color_brewer(palette = "Dark2") + stat_ellipse()
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable",
                            label = "sites") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 20))
View(samples)
beta.ps3 <- beta.ps3 + theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 + scale_color_brewer(palette = "Dark2") + stat_ellipse()
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable",
                            label = "Samples") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 20))
View(samples)
beta.ps3 <- beta.ps3 + theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 + scale_color_brewer(palette = "Dark2") + stat_ellipse()
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable",
                            label = "row_names") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 20))
View(samples)
beta.ps3 <- beta.ps3 + theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 + scale_color_brewer(palette = "Dark2") + stat_ellipse()
beta.ps3 <- plot_ordination(physeq_rel,
                            bx.ord_pcoa_bray,
                            color="Variable",
                            label = "row_names") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 20))
View(samples)
beta.ps3 <- beta.ps3 + theme_bw(base_size = 9) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 + scale_color_brewer(palette = "Dark2") + stat_ellipse()
bx.ord_pcoa_bray
View(otu_mat)
View(otu_mat)
lib <- lapply(packages, require, character.only = TRUE)
beta.ps3
