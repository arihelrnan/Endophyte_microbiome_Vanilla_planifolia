#Make phyloseq object
packages <-c("ape","dplyr","ggplot2","gplots","lme4","miLineage","phangorn","plotly","tidyr","vegan","VennDiagram","metacoder","phyloseq")
lib <- lapply(packages, require, character.only = TRUE)
MetaFile <- "mapping_file.txt"
MetaFile <- read.table(file=MetaFile,header=TRUE,sep="\t",comment.char = "",row.names = 1,check.names = F)
MetaFile <- MetaFile[!apply(is.na(MetaFile) | MetaFile=="",1,all),]
otu_file<-"ITS_taxonomy.otu_table.taxonomy.txt"
otu_table_2 <-  read.table (otu_file,
                          check.names = FALSE,
                          header = TRUE,
                          dec = ".",
                          sep = "\t",
                          row.names = 1,
                          comment.char = "")
taxonomy=otu_table_2[,c(1,17)]
taxonomy$V1=NULL
taxonomy = separate(taxonomy, Taxonomy, into = c("Details", "Clasification"), sep=";")
taxonomy$Details=NULL
taxonomy = separate(taxonomy, Clasification, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Specie"), sep=",")
otu_table_2$Taxonomy=NULL
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

#Check for distribution of OTUs
hist(log10(taxa_sums(phyloseq.rel)))


#################FOR BETA DIVERSITY:  Binomial table: how many OTUS per samples. Filter OTUs from only one sample 

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


#Make bray curtis plot

#Run : 
#1) binary_OTU1 con Raup
#2)rel con Bray
#3)permanova para ambos y para los dos genes

## Not run= physeq_rel <- microbiome::transform(physeq, "compositional")
#Change the matrix ordenation with Raup-Crick
bx.ord_pcoa_bray <- ordinate(phyloseq.rel, "NMDS", "bray")
plot_scree(bx.ord_pcoa_bray) + theme_bw()
beta.ps3 <- plot_ordination(phyloseq.rel,
                            bx.ord_pcoa_bray,
                            color="Organ",
                            shape ="Variables", 
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

#Make Raup-Crick plot
bx.ord_pcoa_raup <- ordinate(binary_table, "NMDS", "raup")
plot_scree(bx.ord_pcoa_raup) + theme_bw()
beta.raup <- plot_ordination(binary_table,
                            bx.ord_pcoa_raup,
                            color="Organ",
                            shape ="Variables", 
                            label = "Sites") +
  geom_point(size= 4) +
  theme(plot.title = element_text(hjust = 0, size = 12))
beta.raup <- beta.raup + theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beta.raup + scale_color_brewer(palette = "Dark2") + stat_ellipse()
beta.raup
beta.raup + stat_ellipse()
beta.raup +  theme_bw() +
  stat_ellipse()

#Make PERMANOVA test
library(vegan)
metadf.bx <- data.frame(sample_data(phyloseq.rel))
bray_ps.bxn <- phyloseq::distance(physeq = phyloseq.rel, method = "bray")
adonis.test <- adonis(bray_ps.bxn ~ Organ*Variables, data = metadf.bx)
adonis.test

metadf.rp <- data.frame(sample_data(binary_table))
raup_ps.bxn <- phyloseq::distance(physeq = binary_table, method = "raup")
adonis.test_2 <- adonis(raup_ps.bxn ~ Organ*Variables, data = metadf.rp)
adonis.test_2

#Make dendogram and heat plot
View(metadf.bx)
distancia_organo=as.matrix(bray_ps.bxn)
heatmap.2(distancia_organo, key=T, trace="none", ColSideColors = c(rep("aquamarine",2),rep("maroon1",2),rep("blue",2),rep("red",3),rep("green",3),rep("darkgoldenrod3",2)), RowSideColors = c(rep("aquamarine",2),rep("maroon1",2),rep("blue",2),rep("red",3),rep("green",3),rep("darkgoldenrod3",2)))

#Calculate alpha diversity for relative abundance data
MetaFile$Richness <- estimate_richness(phyloseq.rel, split = TRUE, measures = "Observed")
MetaFile$Alpha <- estimate_richness(phyloseq.rel, split = TRUE, measures = "Shannon")
Alpha_diversity <- estimate_richness(phyloseq.rel, split = TRUE, measures = c("Observed","Chao1","Shannon","InvSimpson"))

#Use Kruskal-Wallis test for evaluate diference diversity between symptomatic, asymptomatic and wild samples groups
kruskal_rel_variables_obs <- kruskal.test(Alpha_diversity$Observed~Variables, MetaFile)
kruskal_rel_variables.shan <- kruskal.test(Alpha_diversity$Shannon~Variables, MetaFile)
#Use Kruskal-Wallis test for evaluate diference diversity between organ groups (stem and root)
kruskal_rel_organ_obs <- kruskal.test(Alpha_diversity$Observed~Organ, MetaFile)
kruskal_rel_organ_shan <- kruskal.test(Alpha_diversity$Shannon~Organ, MetaFile)

#Calculate alpha diversity for binary data
MetaFile$Richness.bin <- estimate_richness(binary_table, split = TRUE, measures = "Observed")
MetaFile$Alpha.bin <- estimate_richness(binary_table, split = TRUE, measures = "Shannon")
Alpha_diversity.bin <- estimate_richness(binary_table, split = TRUE, measures = c("Observed","Chao1","Shannon","InvSimpson"))

#Use Kruskal-Wallis test for evaluate diference diversity between symptomatic, asymptomatic and wild samples groups
kruskal_bin_variables_obs <- kruskal.test(Alpha_diversity.bin$Observed~Variables, MetaFile)
kruskal_bin_variables_shan <- kruskal.test(Alpha_diversity.bin$Shannon~Variables, MetaFile)
#Use Kruskal-Wallis test for evaluate diference diversity between organ groups (stem and root)
kruskal_bin_organ_obs <- kruskal.test(Alpha_diversity.bin$Observed~Organ, MetaFile)
kruskal_bin_organ_shan <- kruskal.test(Alpha_diversity.bin$Shannon~Organ, MetaFile)

shapiro.test(Alpha_diversity.bin)
