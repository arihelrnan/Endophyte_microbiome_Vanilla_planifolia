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

##Calculate alpha diversity for relative abundance data
MetaFile$Richness <- estimate_richness(relative_data, split = TRUE, measures = "Observed")
MetaFile$Alpha <- estimate_richness(relative_data, split = TRUE, measures = "Shannon")
Alpha_diversity <- estimate_richness(relative_data, split = TRUE, measures = c("Observed","Chao1","Shannon","InvSimpson"))
write.table(Alpha_diversity, "../../Data/Alpha_diversity/Alpha_ITS_Relative_data.tab", sep = "\t")

#Use Kruskal-Wallis test for evaluate diference diversity between symptomatic, asymptomatic and wild samples groups
kruskal_rel_variables_obs <- kruskal.test(Alpha_diversity$Observed~State, MetaFile)
kruskal_rel_variables_shan <- kruskal.test(Alpha_diversity$Shannon~State, MetaFile)
#Use Kruskal-Wallis test for evaluate diference diversity between organ groups (stem and root)
kruskal_rel_organ_obs <- kruskal.test(Alpha_diversity$Observed~Organ, MetaFile)
kruskal_rel_organ_shan <- kruskal.test(Alpha_diversity$Shannon~Organ, MetaFile)

#Calculate alpha diversity for binary data
MetaFile$Richness.bin <- estimate_richness(binary_data, split = TRUE, measures = "Observed")
MetaFile$Alpha.bin <- estimate_richness(binary_data, split = TRUE, measures = "Shannon")
Alpha_diversity.bin <- estimate_richness(binary_data, split = TRUE, measures = c("Observed","Chao1","Shannon","InvSimpson"))
write.table(Alpha_diversity.bin, "../../Data/Alpha_diversity/Alpha_ITS_Binary_data.tab", sep = "\t")

#Use Kruskal-Wallis test for evaluate diference diversity between symptomatic, asymptomatic and wild samples groups
kruskal_bin_variables_obs <- kruskal.test(Alpha_diversity.bin$Observed~State, MetaFile)
kruskal_bin_variables_shan <- kruskal.test(Alpha_diversity.bin$Shannon~State, MetaFile)
#Use Kruskal-Wallis test for evaluate diference diversity between organ groups (stem and root)
kruskal_bin_organ_obs <- kruskal.test(Alpha_diversity.bin$Observed~Organ, MetaFile)
kruskal_bin_organ_shan <- kruskal.test(Alpha_diversity.bin$Shannon~Organ, MetaFile)

#Male plot with  alpha diversity (Shannon?s diversity and Observed Richness)
alpha_meas = c("Observed", "Shannon")
Alpha_diversity_plot <- plot_richness(binary_data, color = "State", measures=alpha_meas, x = "Organ", title = "Alpha diversity with ITS data")
Alpha_diversity_plot <- Alpha_diversity_plot + geom_point(size=3, alpha=0.7)
Alpha_diversity_plot <- Alpha_diversity_plot + theme_light()

Alpha_diversity_plot

png("../../Figures/Alpha_diversity/Alpha_diversity_with_ITS.png") 
print(Alpha_diversity_plot)
dev.off()


kruskal_bin_variables_obs
kruskal_bin_variables_shan
kruskal_bin_organ_obs
kruskal_bin_organ_shan