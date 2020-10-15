file_name <- "OTUs_Table-norm.tab"                #<--- CHANGE ACCORDINGLY
######                  NO CHANGES ARE NEEDED BELOW THIS LINE               ######
##################################################################################
######                        Diversity Functions                           ######
##################################################################################
# Calculate the species richness in a sample
Species.richness <- function(x)
{
  # Count only the OTUs that are present >0.5 normalized counts (normalization produces real values for counts)
  count=sum(x[x>0.5]^0)
  return(count)
}
# Calculate the Shannon diversity index
Shannon.entropy <- function(x)
{
  total=sum(x)
  se=-sum(x[x>0]/total*log(x[x>0]/total))
  return(se)
}
# Calculate the effective number of species for Shannon
Shannon.effective <- function(x)
{
  total=sum(x)
  se=round(exp(-sum(x[x>0]/total*log(x[x>0]/total))),digits =2)
  return(se)
}
# Calculate the Simpson diversity index
Simpson.concentration <- function(x)
{
  total=sum(x)
  si=sum((x[x>0]/total)^2)
  return(si)
}
# Calculate the effective number of species for Simpson
Simpson.effective <- function(x)
{
  total=sum(x)
  si=round(1/sum((x[x>0]/total)^2),digits =2)
  return(si)
}
##################################################################################
######                             Main Script                              ######
##################################################################################
# Read a normalized OTU-table without taxonomy
otu_table <- read.table (file_name,
                         check.names = FALSE,
                         header=TRUE,
                         dec=".",
                         sep = "\t",
                         row.names = 1)
# Clean table from empty lines
otu_table <- otu_table[!apply(is.na(otu_table) | otu_table=="",1,all),]
# Order and transpose OTU-table
my_otu_table <- otu_table[,order(names(otu_table))]
my_otu_table <-data.frame(t(my_otu_table))
# Apply diversity functions to table
otus_div_stats<-data.frame(my_otu_table[,0])
otus_div_stats$Richness<-apply(my_otu_table,1,Species.richness)
otus_div_stats$Shannon<-apply(my_otu_table,1,Shannon.entropy)
otus_div_stats$Shannon.effective<-apply(my_otu_table,1,Shannon.effective)
otus_div_stats$Simpson<-apply(my_otu_table,1,Simpson.concentration)
otus_div_stats$Simpson.effective<-apply(my_otu_table,1,Simpson.effective)
otus_div_stats$Evenness <- otus_div_stats$Shannon/log(otus_div_stats$Richness,2)
# Write the results in a file and copy in directory "Serial-Group-Comparisons" if existing
write.table(otus_div_stats, "alpha-diversity.tab", sep="\t", col.names=NA, quote=FALSE)
suppressWarnings (try(write.table(otus_div_stats[c(1,3,5)], "../5.Serial-Group-Comparisons/alpha-diversity.tab", sep = "\t",col.names = NA, quote = FALSE), silent =TRUE))
View(my_otu_table)
View(my_otu_table)
View(otus_div_stats)
library(vegan)
pairwise.wilcox.test(my_otu_table, meta$AgeGroup, p.adjust.method="fdr")
MetaFile <- "mapping_file.txt"
MetaFile <- read.table(file=MetaFile,header=TRUE,sep="\t",comment.char = "",row.names = 1,check.names = F)
MetaFile <- MetaFile[!apply(is.na(MetaFile) | MetaFile=="",1,all),]
pairwise.wilcox.test(my_otu_table, MetaFile$Variables, p.adjust.method="fdr")
pairwise.wilcox.test(my_otu_table$Shannon.effective, MetaFile$Variables, p.adjust.method="fdr")
my_otu_table
otus_div_stats
pairwise.wilcox.test(otus_div_stats$Shannon.effective, MetaFile$Variables, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Simpson.effective, MetaFile$Variables, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Richness, MetaFile$Variables, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Shannon, MetaFile$Variables, p.adjust.method="fdr")
View(MetaFile)
pairwise.wilcox.test(otus_div_stats$Shannon, MetaFile$Tejido, p.adjust.method="fdr")
(otus_div_stats$Shannon, MetaFile$Variables)
kruskal.test(otus_div_stats$Shannon, MetaFile$Variables)
kruskal.test(otus_div_stats$Shannon, MetaFile$Tejido)
kruskal.test(otus_div_stats$Shannon, MetaFile$Estado)
kruskal.test(otus_div_stats$Shannon, MetaFile$Estado)
kruskal.test(otus_div_stats$Shannon, MetaFile$Estado)
?kruskal.test
pairwise.wilcox.test(otus_div_stats$Shannon, MetaFile$Tejido, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Simpson.effective, MetaFile$Tejido, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Shannon.effective, MetaFile$Tejido, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Richness, MetaFile$Tejido, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Richness, MetaFile$Variables, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Shannon.effective, MetaFile$Variables, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Simpson.effective, MetaFile$Variables, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Simpson.effective, MetaFile$Estadi, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Simpson.effective, MetaFile$Estado, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Simpson.effective, MetaFile$Tejido, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Shannon.effective, MetaFile$Tejido, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Richness, MetaFile$Variable, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Richness, MetaFile$Variables, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Simpson.effective, MetaFile$Variables, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Shannon.effective, MetaFile$Variables, p.adjust.method="fdr")
> kruskal.test(otus_div_stats$Simpson.Effective, MetaFile$Tejido)
kruskal.test(otus_div_stats$Simpson.Effective, MetaFile$Tejido)
kruskal.test(otus_div_stats$Simpson.Effective, MetaFile$Tejido)
kruskal.test(otus_div_stats$Simpson.Effective, MetaFile$Tejido)
kruskal.test(otus_div_stats$Simpson.Effective, MetaFile$Tejido)
View(otus_div_stats)
kruskal.test(otus_div_stats$Simpson.effective, MetaFile$Tejido)
kruskal.test(otus_div_stats$Simpson, MetaFile$Tejido)
kruskal.test(otus_div_stats$Shannon.effective, MetaFile$Tejido)
kruskal.test(otus_div_stats$Shannon, MetaFile$Tejido)
kruskal.test(otus_div_stats$Richness, MetaFile$Tejido)
kruskal.test(otus_div_stats$Evenness, MetaFile$Tejido)
kruskal.test(otus_div_stats$Richness, MetaFile$Tejido)
kruskal.test(otus_div_stats$Richness, MetaFile$Variables)
kruskal.test(otus_div_stats$Shannon, MetaFile$Variables)
kruskal.test(otus_div_stats$Simpson, MetaFile$Variables)
pairwise.wilcox.test(otus_div_stats$Simpson.effective, MetaFile$Variables, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Simpson, MetaFile$Variables, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Simpson, MetaFile$Variables, p.adjust.method="fdr")
kruskal.test(otus_div_stats$Simpson, MetaFile$Variables)
kruskal.test(otus_div_stats$Simpson, MetaFile$Variables)
kruskal.test(otus_div_stats$Simpson, MetaFile$TeTejido)
kruskal.test(otus_div_stats$Simpson, MetaFile$Tejido)
kruskal.test(otus_div_stats$Shannon.effective, MetaFile$Tejido)
kruskal.test(otus_div_stats$Richness, MetaFile$Tejido)
kruskal.test(otus_div_stats$Shannon.effective, MetaFile$Tejido)
pairwise.wilcox.test(otus_div_stats$Shannon.effective, MetaFile$Variables, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Shannon.effective, MetaFile$Tejido, p.adjust.method="fdr")
kruskal.test(otus_div_stats$Richness, MetaFile$Tejido)
kruskal.test(otus_div_stats$Shannon.effective, MetaFile$Tejido)
kruskal.test(otus_div_stats$Simpson.effective, MetaFile$Tejido)
kruskal.test(otus_div_stats$Richness, MetaFile$Variables)
kruskal.test(otus_div_stats$Shannon.effective, MetaFile$Variables)
kruskal.test(otus_div_stats$Simpson.effective, MetaFile$Variables)
