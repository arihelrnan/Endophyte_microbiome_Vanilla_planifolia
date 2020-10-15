file_name<-"Otu_table_ITS-tab.txt"                   #<--- CHANGE ACCORDINGLY
method <- 0                                   #<--- CHANGE ACCORDINGLY
#' Please choose the number of samples with the stipest rarefaction curves to be selectively plotted
#' The default number of samples presented seperately is 5
labelCutoff <- 5                              #<--- CHANGE ACCORDINGLY
method <- 0                                   #<--- CHANGE ACCORDINGLY
#' Please choose the number of samples with the stipest rarefaction curves to be selectively plotted
#' The default number of samples presented seperately is 5
labelCutoff <- 13                              #<--- CHANGE ACCORDINGLY
packages <-c("GUniFrac","vegan")
# Function to check whether the package is installed
InsPack <- function(pack)
{
  if ((pack %in% installed.packages()) == FALSE) {
    install.packages(pack,repos = "http://cloud.r-project.org/")
  }
}
# Applying the installation on the list of packages
lapply(packages, InsPack)
# Make the libraries
lib <- lapply(packages, require, character.only = TRUE)
# Check if it was possible to install all required libraries
flag <- all(as.logical(lib))
otu_table <-  read.table (file_name,
                          check.names = FALSE,
                          header = TRUE,
                          dec = ".",
                          sep = "\t",
                          row.names = 1,
                          comment.char = "")
# Clean table from empty lines
otu_table <- otu_table[!apply(is.na(otu_table) | otu_table=="",1,all),]
View(otu_table)
View(otu_table)
otu_table$V1.ITS=NULL
otu_table$V5.ITS=NULL
taxonomy <- as.vector(otu_table$taxonomy)
# Delete column with taxonomy information in dataframe
otu_table$taxonomy <- NULL
# Calculate the minimum sum of all columns/samples
min_sum <- min(colSums(otu_table))
if (method == 0) {
  # Divide each value by the sum of the sample and multiply by the minimal sample sum
  norm_otu_table <- t(min_sum * t(otu_table) / colSums(otu_table))
} else {
  # Rarefy the OTU table to an equal sequencing depth
  norm_otu_table <- Rarefy(t(otu_table),depth = min_sum)
  norm_otu_table <- t(as.data.frame(norm_otu_table$otu.tab.rff))
}
# Calculate relative abundances for all OTUs over all samples
# Divide each value by the sum of the sample and multiply by 100
rel_otu_table <- t(100 * t(otu_table) / colSums(otu_table))
# Re-insert the taxonomy information in normalized counts table
norm_otu_table_tax <- cbind(norm_otu_table,taxonomy)
# Reinsert the taxonomy information in relative abundance table
rel_otu_table_tax <- cbind(rel_otu_table,taxonomy)
################################################################################
# Generate a twosided pdf with a rarefaction curve for all samples and a curve
pdf(file = "RarefactionCurve.pdf")
# Plot the rarefaction curve for all samples
rarefactionCurve <- rarecurve(data.frame(t(otu_table)),
                              step = 20,
                              col = "black",
                              lty = "solid",
                              label=F,
                              xlab = "Number of Reads",
                              ylab = "Number of Species",
                              main = "Rarefaction Curves of All Samples")
# Generate empy vectors for the analysis of the rarefaction curve
slope=vector()
SampleID=vector()
# Iterate through all samples
for(i in seq_along(rarefactionCurve)) {
  # If the sequencing depth is greater 100 the difference between the last and last-100 richness is calcualted
  richness <- ifelse(length(rarefactionCurve[[i]])>=100,rarefactionCurve[[i]][length(rarefactionCurve[[i]])] - rarefactionCurve[[i]][length(rarefactionCurve[[i]])-100],1000)
  slope<- c(slope,richness)
  SampleID <- c(SampleID,as.character(names(otu_table)[i]))
}
# Generate the output table for rarefaction curve
curvedf <- cbind(SampleID,slope)
order <- order(curvedf[,2],decreasing = TRUE)
# Order the table
curvedf <- curvedf[order(curvedf[,2],decreasing = TRUE),]
# Generates a graph with all samples
# Underestimated cases are shown in red
for ( i in 1:labelCutoff) {
  N <- attr(rarefactionCurve[[order[i]]], "Subsample")
  lines(N, rarefactionCurve[[order[i]]],col="red")
}
# Determine the plotting width and height
Nmax <- sapply(rarefactionCurve, function(x) max(attr(x, "Subsample")))
Smax <- sapply(rarefactionCurve, max)
# Creates an empty plot for rarefaction curves of underestimated cases
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Number of Reads",
     ylab = "Number of Species", type = "n", main=paste(labelCutoff,"- most undersampled cases"))
for (i in 1:labelCutoff) {
  N <- attr(rarefactionCurve[[order[i]]], "Subsample")
  lines(N, rarefactionCurve[[order[i]]],col="red")
  text(max(attr(rarefactionCurve[[order[i]]],"Subsample")),max(rarefactionCurve[[order[i]]]), curvedf[i,1],cex=0.6)
}
dev.off()
#################################################################################
######                        Write Output Files                           ######
#################################################################################
# Write the normalized table in a file and copy in directories alpha-diversity and beta-diversity if existing
write.table(norm_otu_table, "OTUs_Table-norm.tab", sep = "\t",col.names = NA, quote = FALSE)
suppressWarnings (try(write.table(norm_otu_table, "../2.Alpha-Diversity/OTUs_Table-norm.tab", sep = "\t",col.names = NA, quote = FALSE), silent =TRUE))
suppressWarnings (try(write.table(norm_otu_table, "../3.Beta-Diversity/OTUs_Table-norm.tab", sep = "\t",col.names = NA, quote = FALSE), silent =TRUE))
# Write the normalized table with taxonomy in a file
write.table(norm_otu_table_tax, "OTUs_Table-norm-tax.tab", sep = "\t",col.names = NA, quote = FALSE)
# Write the normalized relative abundance table in a file and copy in directory Serial-Group-Comparisons if existing
write.table(rel_otu_table, "OTUs_Table-norm-rel.tab", sep = "\t",col.names = NA, quote = FALSE)
suppressWarnings (try(write.table(rel_otu_table, "../5.Serial-Group-Comparisons/OTUs_Table-norm-rel.tab", sep = "\t",col.names = NA, quote = FALSE), silent =TRUE))
# Write the normalized relative abundance with taxonomy table in a file and copy in directory Taxonomic-Binning if existing
write.table(rel_otu_table_tax, "OTUs_Table-norm-rel-tax.tab", sep ="\t",col.names = NA, quote = FALSE)
suppressWarnings (try(write.table(rel_otu_table_tax, "../4.Taxonomic-Binning/OTUs_Table-norm-rel-tax.tab", sep ="\t",col.names = NA, quote = FALSE), silent =TRUE))
# Write the rarefaction table
write.table(curvedf, "RarefactionCurve.tab", sep ="\t", quote = FALSE, row.names = FALSE)
# Error message
if(!flag) { stop("
It was not possible to install all required R libraries properly.
Please check the installation of all required libraries manually.\n
Required libaries: GUniFrac")
}