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
library(vegan)
#' Please give the file name of the normalized OTU-table (without taxonomic classification)
file_name <- "OTUs_Table-norm.tab"                #<--- CHANGE ACCORDINGLY
#' The abundance filtering cutoff
eff.cutoff <- 0.0025 # this is the default value for Effective Richness (0.25%)
#' The normalized depth cutoff
norm.cutoff <- 1000 # this is the default value for Standard Richness (1000)
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
# Calculate the Effective species richness in each individual sample
Eff.Species.richness <- function(x)
{
  # Count only the OTUs that are present more than the set proportion
  total=sum(x)
  count=sum(x[x/total>eff.cutoff]^0)
  return(count)
}
# Calculate the Normalized species richness in each individual sample
Norm.Species.richness <- function(x)
{
  # Count only the OTUs that are present >0.5 normalized counts (normalization produces real values for counts)
  # Given a fixed Normalization reads depth
  total=sum(x)
  count=sum(x[norm.cutoff*x/total>0.5]^0)
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
otus_div_stats$Normalized.Richness<-apply(my_otu_table,1,Norm.Species.richness)
otus_div_stats$Effective.Richness<-apply(my_otu_table,1,Eff.Species.richness)
otus_div_stats$Shannon.Index<-apply(my_otu_table,1,Shannon.entropy)
otus_div_stats$Shannon.Effective<-apply(my_otu_table,1,Shannon.effective)
otus_div_stats$Simpson.Index<-apply(my_otu_table,1,Simpson.concentration)
otus_div_stats$Simpson.Effective<-apply(my_otu_table,1,Simpson.effective)
otus_div_stats$Evenness <- otus_div_stats$Shannon.Index/log(otus_div_stats$Richness,2)
# Write the results in a file and copy in directory "Serial-Group-Comparisons" if existing
write.table(otus_div_stats, "alpha-diversity.tab", sep="\t", col.names=NA, quote=FALSE)
suppressWarnings (try(write.table(otus_div_stats[c(1,3,5)], "../5.Serial-Group-Comparisons/alpha-diversity.tab", sep = "\t",col.names = NA, quote = FALSE), silent =TRUE))
##################################################################################
######                          End of Script                               ######
##################################################################################
MetaFile <- "mapping_file.txt"
MetaFile <- read.table(file=MetaFile,header=TRUE,sep="\t",comment.char = "",row.names = 1,check.names = F)
MetaFile <- MetaFile[!apply(is.na(MetaFile) | MetaFile=="",1,all),]
View(MetaFile)
View(my_otu_table)
View(my_otu_table)
View(otus_div_stats)
pairwise.wilcox.test(otus_div_stats$Shannon.Effective, MetaFile$Variable, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Shannon.Effective, MetaFile$Variable, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Simpson.Effective, MetaFile$Variable, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Richness, MetaFile$Variable, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Richness, MetaFile$Tejido, p.adjust.method="fdr")
pairwise.wilcox.test(otus_div_stats$Richness, MetaFile$Estado, p.adjust.method="fdr")
t.test(otus_div_stats$Richness, MetaFile$Estado, p.adjust.method="fdr")
t.test(otus_div_stats$Richness, MetaFile$Estado)
t.test(otus_div_stats$Richness, MetaFile$Variable)
Wilconxon_tejido_richness=pairwise.wilcox.test(otus_div_stats$Richness, MetaFile$Tejido, p.adjust.method="fdr")
Wilconxon_tejido_shannon=pairwise.wilcox.test(otus_div_stats$Shannon.Effective, MetaFile$Tejido, p.adjust.method="fdr")
Wilconxon_tejido_simpson=pairwise.wilcox.test(otus_div_stats$Simpson.Effective, MetaFile$Tejido, p.adjust.method="fdr")
Wilconxon_estado_richness=pairwise.wilcox.test(otus_div_stats$Richness, MetaFile$Estado, p.adjust.method="fdr")
Wilconxon_Estado_shannon=pairwise.wilcox.test(otus_div_stats$Shannon.Effective, MetaFile$Estado, p.adjust.method="fdr")
Wilconxon_estado_simpson=pairwise.wilcox.test(otus_div_stats$Simpson.Effective, MetaFile$Estado, p.adjust.method="fdr")
library(vegan)
Wilconxon_variable_simpson=pairwise.wilcox.test(otus_div_stats$Simpson.Effective, MetaFile$Variable, p.adjust.method="fdr")
Wilconxon_variable_Shannon=pairwise.wilcox.test(otus_div_stats$Shannon.Effective, MetaFile$Variable, p.adjust.method="fdr")
Wilconxon_variable_Richness=pairwise.wilcox.test(otus_div_stats$Richness, MetaFile$Variable, p.adjust.method="fdr")
Wilconxon_variable_Shannon
Wilconxon_Estado_Shannon
Wilconxon_estado_Shannon
Wilconxon_Estado_shannon
Wilconxon_tejido_shannon
Wilconxon_tejido_simpson
Wilconxon_tejido_shannon
Wilconxon_tejido_richness
kruskal.test(otus_div_stats$Shannon, MetaFile$Estado)
kruskal.test(otus_div_stats$Shannon.Effective, MetaFile$Estado)
kruskal.test(otus_div_stats$Shannon.Effective, MetaFile$Tejido)
kruskal.test(otus_div_stats$Simpson.Effective, MetaFile$Tejido)
kruskal.test(otus_div_stats$Richness, MetaFile$Tejido)
kruskal.test(otus_div_stats$Richness, MetaFile$Variable)
kruskal.test(otus_div_stats$Shannon.Effective, MetaFile$Variable)
kruskal.test(otus_div_stats$Simpson.Effective, MetaFile$Variable)
Wilconxon_variable_Richness
Wilconxon_variable_simpson
Wilconxon_variable_Shannon
kruskal.test(otus_div_stats$Richness, MetaFile$Variable)
kruskal.test(otus_div_stats$Richness, MetaFile$Tejido)
kruskal.test(otus_div_stats$Shannon.Effective, MetaFile$Tejido)
kruskal.test(otus_div_stats$Simpson.Effective, MetaFile$Tejido)
kruskal.test(otus_div_stats$Simpson.Effective, MetaFile$Tejido)
View(Wilconxon_Estado_shannon)
View(Wilconxon_Estado_shannon)
kruskal.test(otus_div_stats$Richness, MetaFile$Tejido)
kruskal.test(otus_div_stats$Shannon.Effective, MetaFile$Tejido)
kruskal.test(otus_div_stats$Simpson.Effective, MetaFile$Tejido)
kruskal.test(otus_div_stats$Richness, MetaFile$Variable)
kruskal.test(otus_div_stats$Shannon.Effective, MetaFile$Variable)
kruskal.test(otus_div_stats$Simpson.Effective, MetaFile$Variable)
