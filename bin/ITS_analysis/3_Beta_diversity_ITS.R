#Load packpage 
packages <- c("dplyr","ggplot2","plotly","tidyr","vegan","phyloseq")
lib <- lapply(packages, require, character.only = TRUE)

#Load input files
MetaFile <- "../../Data/mapping_file.txt"
MetaFile <- read.table(file=MetaFile,header=TRUE,sep="\t",comment.char = "",row.names = 1,check.names = F)
MetaFile <- MetaFile[!apply(is.na(MetaFile) | MetaFile=="",1,all),]
otu_file_binary <- "../../Data/Data_transformation/binary_table_ITS.txt"
otu_table_binary <-  read.table (otu_file_binary,
                                 check.names = FALSE,
                                 header = TRUE,
                                 dec = ".",
                                 sep = "\t",
                                 row.names = 1,
                                 comment.char = "")
otu_file_relative <- "../../Data/Data_transformation/relative_table_ITS.txt"
otu_table_relative <-  read.table (otu_file_relative,
                                   check.names = FALSE,
                                   header = TRUE,
                                   dec = ".",
                                   sep = "\t",
                                   row.names = 1,
                                   comment.char = "")

####################Make phyloseq objects##################################
#Make taxonomy object
taxonomy=otu_table_binary[,c(1,17)]
taxonomy$V1=NULL
taxonomy = separate(taxonomy, taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Specie"), sep=",")
otu_table_binary$taxonomy=NULL
otu_table_relative$taxonomy=NULL
#Make phyloseq object with relative data
otu_mat_relative=as.matrix(otu_table_relative)
tax_mat=as.matrix(taxonomy)
OTU = otu_table(otu_mat_relative, taxa_are_rows = TRUE)
samples = sample_data(MetaFile)
TAX = tax_table(tax_mat)
relative_data <- phyloseq(OTU, TAX, samples)
#Make phyloseq object with binary data
otu_mat_binary=as.matrix(otu_table_binary)
tax_mat=as.matrix(taxonomy)
OTU = otu_table(otu_mat_binary, taxa_are_rows = TRUE)
samples = sample_data(MetaFile)
TAX = tax_table(tax_mat)
binary_data <- phyloseq(OTU, TAX, samples)

#Calculate Bray Curtis Dissimilarity with abundance relative data and make NMDS with matrix distance
bx.ord_pcoa_bray <- ordinate(relative_data, "NMDS", "bray")
plot_scree(bx.ord_pcoa_bray) + theme_bw()
beta.ps3 <- plot_ordination(relative_data,
                            bx.ord_pcoa_bray,
                            color="Organ",
                            shape ="State",
                            title="NDMS using Bray-Curtis Dissimilarity with ITS data",
                            label = "Sites") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 12))
beta.ps3 <- beta.ps3 + theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 <- beta.ps3 +  theme_bw() +
  stat_ellipse() + scale_color_manual(values = c("orange2", "green4"))
pdf("../../Figures/Beta_diversity/beta_bc_ITS.pdf") 
print(beta.ps3)
dev.off()
png("../../Figures/Beta_diversity/beta_bc_ITS.png") 
print(beta.ps3)
dev.off()

##Calculate Raup-Crick Dissimilarity with binary data (presence-absence) and make NMDS with matrix distance  
bx.ord_pcoa_raup <- ordinate(binary_data, "NMDS", "raup")
plot_scree(bx.ord_pcoa_raup) + theme_bw()
beta.raup <- plot_ordination(binary_data,
                             bx.ord_pcoa_raup,
                             color="Organ",
                             shape ="State", 
                             title="NDMS using Raup-Crick Dissimilarity with ITS data",
                             label = "Sites") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 12))
beta.raup <- beta.raup + theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.raup <- beta.raup +  theme_bw() +
  stat_ellipse() + scale_color_manual(values = c("orange2", "green4"))
pdf("../../Figures/Beta_diversity/beta_raup_ITS.pdf") 
print(beta.raup)
dev.off()
png("../../Figures/Beta_diversity/beta_raup_ITS.png") 
print(beta.raup)
dev.off()

#Make PERMANOVA test to Bray Curtis Dissimilarity
metadf.bx <- data.frame(sample_data(relative_data))
bray_ps.bxn <- phyloseq::distance(physeq = relative_data, method = "bray")
adonis.test <- adonis(bray_ps.bxn ~ Organ*State, data = metadf.bx)
adonis.test
write.table(as.data.frame(adonis.test[[1]]), "../../Data/Beta_diversity/PERMANOVA_Relative_data_ITS.tab", sep = "\t")

#Make PERMANOVA test to Raup-Crick Dissimilarity
metadf.rp <- data.frame(sample_data(binary_data))
raup_ps.bxn <- phyloseq::distance(physeq = binary_data, method = "raup")
adonis.test_2 <- adonis(raup_ps.bxn ~ Organ*State, data = metadf.rp)
adonis.test_2
write.table(as.data.frame(adonis.test_2[[1]]), "../../Data/Beta_diversity/PERMANOVA_Binary_data_ITS.tab", sep = "\t")

