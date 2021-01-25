#Make phyloseq object
packages <-c("ape","dplyr","ggplot2","gplots","lme4","miLineage","phangorn","plotly","tidyr","vegan","VennDiagram","metacoder","phyloseq","dunn.test")
lib <- lapply(packages, require, character.only = TRUE)
MetaFile <- "../../Data/mapping_file.txt"
MetaFile <- read.table(file=MetaFile,header=TRUE,sep="\t",comment.char = "",row.names = 1,check.names = F)
MetaFile <- MetaFile[!apply(is.na(MetaFile) | MetaFile=="",1,all),]
otu_file<-"../../Data/16S.otu_table.taxonomy.txt"
otu_table_2 <-  read.table (otu_file,
                            check.names = FALSE,
                            header = TRUE,
                            dec = ".",
                            sep = "\t",
                            row.names = 1,
                            comment.char = "")
taxonomy=otu_table_2[,c(1,17)]
taxonomy$V1=NULL
taxonomy_2 <- taxonomy
taxonomy = separate(taxonomy, taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Specie"), sep=";")
otu_table_2$taxonomy=NULL
otu_mat=as.matrix(otu_table_2)
tax_mat=as.matrix(taxonomy)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
samples = sample_data(MetaFile)
TAX = tax_table(tax_mat)
physeq <- phyloseq(OTU, TAX, samples)

# Step 2. DATA TRANSFORMATION

# Transform read counts into  % relative abundances: multiply by 1000 and transform to next integer so it looks like read count
phyloseq.rel = transform_sample_counts(physeq, function(x) 100000 * x/sum(x))
otu_table(phyloseq.rel) = ceiling(otu_table(phyloseq.rel, "matrix")) # transform to next integer so it looks like read count
otu_table(phyloseq.rel) # check otu table
phyloseq.rel # check project

#check if any OTUs are still counted as relative abundance and not integer
any(taxa_sums(phyloseq.rel) < 1)
ntaxa(phyloseq.rel)

#check distribution of how many reads/OTU, reads/sample 

sum(taxa_sums(phyloseq.rel))

# Step 3. check distribution of how many reads/OTU, reads/sample: Plot number of reads per OTU / samples 


readsumsdf = data.frame(no.reads = sort(taxa_sums(physeq), TRUE), sorted = 1:ntaxa(physeq), type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(no.reads = sort(sample_sums(physeq), TRUE), sorted = 1:nsamples(physeq), type = "Samples"))
title = ""
p = ggplot(readsumsdf, aes(x = sorted, y = no.reads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

readsumsdf = data.frame(no.reads = sort(taxa_sums(phyloseq.rel), TRUE), sorted = 1:ntaxa(phyloseq.rel), type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(no.reads = sort(sample_sums(phyloseq.rel), TRUE), sorted = 1:nsamples(phyloseq.rel), type = "Samples"))
title = ""
p = ggplot(readsumsdf, aes(x = sorted, y = no.reads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")


#Make Binomial table: how many OTUS per samples. Filter OTUs from only one sample 

binary_table = transform_sample_counts(phyloseq.rel, function(x, minthreshold=0){
  x[x > minthreshold] <- 1
  return(x)})

head(otu_table(binary_table))


# Remove OTUs that appear only in 1 sample (using presence/absence)

any(taxa_sums(binary_table) == 1)
otu_table(prune_taxa(taxa_sums(binary_table) <= 1, binary_table))
binary_table_OTU2 <- prune_taxa(taxa_sums(binary_table) > 1, binary_table)

binary_table
binary_table_OTU2
phyloseq.rel

#Save binary table as data frame
binary_table_16S <- as.data.frame(otu_table(binary_table))
binary_table_16S <- cbind(binary_table_16S, taxonomy_2)
write.table(binary_table_16S, "../../Data/binary_table_16S.txt", sep = "\t")

#Make bray curtis plot

bx.ord_pcoa_bray <- ordinate(phyloseq.rel, "NMDS", "bray")
plot_scree(bx.ord_pcoa_bray) + theme_bw()
beta.ps3 <- plot_ordination(phyloseq.rel,
                            bx.ord_pcoa_bray,
                            color="Organ",
                            shape ="State",
                            title="NDMS using Bray-Curtis Dissimilarity with 16SrRNA data",
                            label = "Sites") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 12))
beta.ps3 <- beta.ps3 + theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.ps3 <- beta.ps3 +  theme_bw() +
  stat_ellipse() + scale_color_manual(values = c("orange2", "green4"))
pdf("../../Figures/beta_bc_16S.pdf") 
print(beta.ps3)
dev.off()
png("../../Figures/beta_bc_16S.png") 
print(beta.ps3)
dev.off()

#Make Raup-Crick plot
bx.ord_pcoa_raup <- ordinate(binary_table, "NMDS", "raup")
plot_scree(bx.ord_pcoa_raup) + theme_bw()
beta.raup <- plot_ordination(binary_table,
                             bx.ord_pcoa_raup,
                             color="Organ",
                             shape ="State", 
                             title="NDMS using Raup-Crick Dissimilarity with 16SrRNA data",
                             label = "Sites") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 12))
beta.raup <- beta.raup + theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.raup <- beta.raup +  theme_bw() +
  stat_ellipse() + scale_color_manual(values = c("orange2", "green4"))
pdf("../../Figures/beta_raup_16S.pdf") 
print(beta.raup)
dev.off()
png("../../Figures/beta_raup_16S.png") 
print(beta.raup)
dev.off()

#Make PERMANOVA test
library(vegan)
metadf.bx <- data.frame(sample_data(phyloseq.rel))
bray_ps.bxn <- phyloseq::distance(physeq = phyloseq.rel, method = "bray")
adonis.test <- adonis(bray_ps.bxn ~ Organ*State, data = metadf.bx)
adonis.test
write.table(as.data.frame(adonis.test[[1]]), "../../Data/PERMANOVA_Relative_data_16S.tab", sep = "\t")


metadf.rp <- data.frame(sample_data(binary_table))
raup_ps.bxn <- phyloseq::distance(physeq = binary_table, method = "raup")
adonis.test_2 <- adonis(raup_ps.bxn ~ Organ*State, data = metadf.rp)
adonis.test_2
write.table(as.data.frame(adonis.test_2[[1]]), "../../Data/PERMANOVA_Binary_data_16S.tab", sep = "\t")


#Calculate alpha diversity for relative abundance data
MetaFile$Richness <- estimate_richness(phyloseq.rel, split = TRUE, measures = "Observed")
MetaFile$Alpha <- estimate_richness(phyloseq.rel, split = TRUE, measures = "Shannon")
Alpha_diversity <- estimate_richness(phyloseq.rel, split = TRUE, measures = c("Observed","Chao1","Shannon","InvSimpson"))
write.table(Alpha_diversity, "../../Data/Alpha_16S_Relative_data.tab", sep = "\t")

#Use Kruskal-Wallis test for evaluate diference diversity between symptomatic, asymptomatic and wild samples groups
kruskal_rel_variables_obs <- kruskal.test(Alpha_diversity$Observed~Variables, MetaFile)
kruskal_rel_variables_shan <- kruskal.test(Alpha_diversity$Shannon~Variables, MetaFile)
#Use Kruskal-Wallis test for evaluate diference diversity between organ groups (stem and root)
kruskal_rel_organ_obs <- kruskal.test(Alpha_diversity$Observed~Organ, MetaFile)
kruskal_rel_organ_shan <- kruskal.test(Alpha_diversity$Shannon~Organ, MetaFile)

#Calculate alpha diversity for binary data
MetaFile$Richness.bin <- estimate_richness(binary_table, split = TRUE, measures = "Observed")
MetaFile$Alpha.bin <- estimate_richness(binary_table, split = TRUE, measures = "Shannon")
Alpha_diversity.bin <- estimate_richness(binary_table, split = TRUE, measures = c("Observed","Chao1","Shannon","InvSimpson"))
write.table(Alpha_diversity.bin, "../../Data/Alpha_16S_Binary_data.tab", sep = "\t")

#Use Kruskal-Wallis test for evaluate diference diversity between symptomatic, asymptomatic and wild samples groups
kruskal_bin_variables_obs <- kruskal.test(Alpha_diversity.bin$Observed~, MetaFile)
kruskal_bin_variables_shan <- kruskal.test(Alpha_diversity.bin$Shannon~Variables, MetaFile)
#Use Kruskal-Wallis test for evaluate diference diversity between organ groups (stem and root)
kruskal_bin_organ_obs <- kruskal.test(Alpha_diversity.bin$Observed~Organ, MetaFile)
kruskal_bin_organ_shan <- kruskal.test(Alpha_diversity.bin$Shannon~Organ, MetaFile)

#Male plot with  alpha diversity (Shannon?s diversity and Observed Richness)
alpha_meas = c("Observed", "Shannon")
Alpha_diversity_plot <- plot_richness(binary_table, color = "State", measures=alpha_meas, x = "Organ", title = "Alpha diversity with 16SrRNA data")
Alpha_diversity_plot <- Alpha_diversity_plot + geom_point(size=3, alpha=0.7)
Alpha_diversity_plot <- Alpha_diversity_plot + theme_light()

Alpha_diversity_plot

png("../../Figures/Alpha_diversity_with_16SrRNA.png") 
print(Alpha_diversity_plot)
dev.off()


kruskal_bin_variables_obs
kruskal_bin_variables_shan
kruskal_bin_organ_obs
kruskal_bin_organ_shan

#Analize the number of OTU across samples groups
ntaxa(sample_variables(binary_table))
sample_variables(binary_table)
rowSums(binary_table)
Stem_OTUs <- subset_samples(binary_table, Organ == "Stem")
Root_OTUs <- subset_samples(binary_table, Organ == "Root")
Stem_OTUs
Root_OTUs
head(otu_table(Stem_OTUs))
head(otu_table(Root_OTUs))
Asymptomatic_OTUs <- subset_samples(binary_table, Variables == "Asymptomatic")
Symptomatic_OTUs <- subset_samples(binary_table, Variables == "Symptomatic")
Wild_OTUs <- subset_samples(binary_table, Variables == "Wild")
head(otu_table(Asymptomatic_OTUs))
head(otu_table(Symptomatic_OTUs))
head(otu_table(Wild_OTUs))

#Taxonomic analysis in phylum and family class

plot_bar(Asymptomatic_OTUs, x="Variables", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
plot_bar(Symptomatic_OTUs, x="Variables", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
plot_bar(Wild_OTUs, x="Variables", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

plot_bar(Stem_OTUs, x="Variables", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
plot_bar(Root_OTUs, x="Variables", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

plot_bar(binary_table, x="Variables", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

subset_taxa(binary_table, Order == " o__Actinomycetales")
subset_taxa(binary_table, Order == " o__Rhizobiales")
subset_taxa(binary_table, Order == " o__Planctomycetales")
subset_taxa(binary_table, Order == " o__\"Sphingobacteriales\"")

top5P.phylum = sort(tapply(taxa_sums(binary_table), tax_table(binary_table)[, "Phylum"], sum), TRUE)[1:5]
top5P.family = sort(tapply(taxa_sums(binary_table), tax_table(binary_table)[, "Order"], sum), TRUE)[1:5]
top5P.family
