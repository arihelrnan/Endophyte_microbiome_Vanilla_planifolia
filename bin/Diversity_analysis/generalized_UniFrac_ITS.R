input_otu = "OTUs_Table-norm.tab"              #<--- CHANGE ACCORDINGLY !!!
#' Please give the name of the meta-file that contains individual sample information
input_meta = "mapping_file.txt"                #<--- CHANGE ACCORDINGLY !!!
#' Please give the name of the phylogenetic tree constructed from the OTU sequences
input_tree = "OTUs-Tree.tre"                   #<--- CHANGE ACCORDINGLY !!!
#' Please give the column name (in the mapping file) of the categorical variable to be used for comparison (e.g. Genotype)
group_name = "Variable"
label_samples = 1
label_id =c("")
kmers_limit=100
packages <-c("ade4","GUniFrac","phangorn","cluster","fpc")
# Function to check whether the package is installed
InsPack <- function(pack)
{
  if ((pack %in% installed.packages()) == FALSE) {
    install.packages(pack,repos ="http://cloud.r-project.org/")
  }
}
# Applying the installation on the list of packages
lapply(packages, InsPack)
# Make the libraries
lib <- lapply(packages, require, character.only = TRUE)
# Check if it was possible to install all required libraries
flag <- all(as.logical(lib))
otu_file <- read.table (file = input_otu, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "")
# Clean table from empty lines
otu_file <- otu_file[!apply(is.na(otu_file) | otu_file =="",1,all),]
# Load the mapping file containing individual sample information (sample names in the first column)
meta_file <- read.table (file = input_meta, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "")
# Save the column names of the mapping file
mappingVar <- names(meta_file)
# Clean table from empty lines
meta_file <- data.frame(meta_file[!apply(is.na(meta_file) | meta_file=="",1,all),],row.names=row.names(meta_file))
tree_file <- read.tree(input_tree)
# Create the directory where all output files are saved (is named after the target group name set above for comparisons)
dir.create(group_name)
otu_file <- otu_file[,order(names(otu_file))]
# Transpose OTU-table and convert format to a data frame
otu_file <- data.frame(t(otu_file))
# Order the mapping file by sample names (ascending)
meta_file <- data.frame(meta_file[order(row.names(meta_file)),],row.names=row.names(meta_file))
# Assign the column names to the modified mapping file
names(meta_file) <- mappingVar
# Save the position of the target group name in the mapping file
meta_file_pos <- which(colnames(meta_file) == group_name)
# Select metadata group based on the pre-set group name
all_groups <- as.factor(meta_file[,meta_file_pos])
# Root the OTU tree at midpoint
rooted_tree <- midpoint(tree_file)
absent <- rooted_tree$tip.label[!(rooted_tree$tip.label %in% colnames(otu_file))]
if (length(absent) != 0) {
  rooted_tree <- drop.tip(rooted_tree, absent)}
tip.label <- rooted_tree$tip.label
otu_file <- otu_file[, tip.label]
ntip <- length(tip.label)
nbr <- nrow(rooted_tree$edge)
edge <- rooted_tree$edge
edge2 <- edge[, 2]
br.len <- rooted_tree$edge.length
unifracs <- GUniFrac(otu_file, rooted_tree, alpha = c(0.0,0.5,1.0))$unifracs
# Weight on abundant lineages so the distance is not dominated by highly abundant lineages with 0.5 having the best power
unifract_dist <- unifracs[, , "d_0.5"]
all_dist_matrix <- as.dist(unifract_dist)
# Apply a hierarchical cluster analysis on the distance matrix based on the Ward's method
all_fit <- hclust(all_dist_matrix, method = "ward.D2")
# Generates a tree from the hierarchically generated object
tree <- as.phylo(all_fit)
my_tree_file_name <- paste(group_name,"/phylogram.pdf",sep="")
plot_color<-rainbow(length(levels(all_groups)))[all_groups]
# Save the generated phylogram in a pdf file
pdf(my_tree_file_name)
# The tree is visualized as a Phylogram color-coded by the selected group name
plot(tree, type = "phylogram",use.edge.length = TRUE, tip.color = (plot_color), label.offset = 0.01)
print.phylo(tree)
axisPhylo()
tiplabels(pch = 16, col = plot_color)
dev.off()
#################            Build NMDS plot           ########################
# Generated figures are saved in a pdf file
file_name <- paste(group_name,"beta-diversity.pdf",sep="_")
pdf(paste(group_name,"/",file_name,sep=""))
# Calculate the significance of variance to compare multivariate sample means (including two or more dependent variables)
# Omit cases where there isn't data for the sample (NA)
all_groups_comp <- all_groups[!is.na(all_groups)]
unifract_dist_comp <- unifract_dist[, !is.na(all_groups)]
adonis<-adonis(as.dist(unifract_dist_comp) ~ all_groups_comp)
all_groups_comp<-factor(all_groups_comp,levels(all_groups_comp)[unique(all_groups_comp)])
# Calculate and display the MDS plot (Multidimensional Scaling plot)
s.class(
  cmdscale(unifract_dist_comp, k = 2), col = unique(plot_color), cpoint =
    2, fac = all_groups_comp, sub = paste("MDS plot of Microbial Profiles\n(p-value ",adonis[[1]][6][[1]][1],")",sep="")
)
if (label_samples==1) {
  lab_samples <- row.names(cmdscale(unifract_dist_comp, k = 2))
  ifelse (label_id != "",lab_samples <- replace(lab_samples, !(lab_samples %in% label_id), ""), lab_samples)
  text(cmdscale(unifract_dist_comp, k = 2),labels=lab_samples,cex=0.7,adj=c(-.1,-.8))
}
# Calculate and display the NMDS plot (Non-metric Multidimensional Scaling plot)
meta <- metaMDS(unifract_dist_comp,k = 2)
s.class(
  meta$points, col = unique(plot_color), cpoint = 2, fac = all_groups_comp,
  sub = paste("metaNMDS plot of Microbial Profiles\n(p-value ",adonis[[1]][6][[1]][1],")",sep="")
)
if (label_samples==1){
  lab_samples <- row.names(meta$points)
  ifelse (label_id != "",lab_samples <- replace(lab_samples, !(lab_samples %in% label_id), ""), lab_samples)
  text(meta$points,labels=lab_samples,cex=0.7,adj=c(-.1,-.8))
}
#close the pdf file
dev.off()
###############          NMDS for pairwise analysis        ###################
# This plot is only generated if there are more than two groups included in the comparison
# Calculate the pairwise significance of variance for group pairs
# Get all groups contained in the mapping file
unique_groups <- levels(all_groups_comp)
if (dim(table(unique_groups)) > 2) {
  # Initialise vector and lists
  pVal = NULL
  pairedMatrixList <- list(NULL)
  pair_1_list <- NULL
  pair_2_list <- NULL
  for (i in 1:length(combn(unique_groups,2)[1,])) {
    # Combine all possible pairs of groups
    pair_1 <- combn(unique_groups,2)[1,i]
    pair_2 <- combn(unique_groups,2)[2,i]
    # Save pairs information in a vector
    pair_1_list[i] <- pair_1
    pair_2_list[i] <- pair_2
    # Generate a subset of all samples within the mapping file related to one of the two groups
    inc_groups <-
      rownames(subset(meta_file, meta_file[,meta_file_pos] == pair_1
                      |
                        meta_file[,meta_file_pos] == pair_2))
    # Convert UniFrac distance matrix to data frame
    paired_dist <- as.data.frame(unifract_dist_comp)
    # Save all row names of the mapping file
    row_names <- rownames(paired_dist)
    # Add row names to the distance matrix
    paired_dist <- cbind(row_names,paired_dist)
    # Generate distance matrix with samples of the compared groups (column-wise)
    paired_dist <- paired_dist[sapply(paired_dist[,1], function(x) all(x %in% inc_groups)),]
    # Remove first column with unnecessary group information
    paired_dist[,1] <- NULL
    paired_dist <- rbind(row_names,paired_dist)
    # Generate distance matrix with samples of the compared group (row-wise)
    paired_dist <- paired_dist[,sapply(paired_dist[1,], function(x) all(x %in% inc_groups))]
    # Remove first row with unnecessary group information
    paired_dist <- paired_dist[-1,]
    # Convert generated distance matrix to data type matrix (needed by multivariate analysis)
    paired_matrix <- as.matrix(paired_dist)
    class(paired_matrix) <- "numeric"
    # Save paired matrix in list
    pairedMatrixList[[i]] <- paired_matrix
    # Applies multivariate analysis to a pair out of the selected groups
    adonis <- adonis(paired_matrix ~ all_groups_comp[all_groups_comp == pair_1 |
                                                       all_groups_comp == pair_2])
    # List p-values
    pVal[i] <- adonis[[1]][6][[1]][1]
  }
  # Adjust p-values for multiple testing according to Benjamini-Hochberg method
  pVal_BH <- p.adjust(pVal,method="BH", n=length(pVal))
  # Generated NMDS plots are stored in one pdf file called "pairwise-beta-diversity-nMDS.pdf"
  file_name <- paste(group_name,"pairwise-beta-diversity.pdf",sep="_")
  pdf(paste(group_name,"/",file_name,sep=""))
  for(i in 1:length(combn(unique_groups,2)[1,])){
    meta <- metaMDS(pairedMatrixList[[i]], k = 2)
    s.class(
      meta$points,
      col = rainbow(length(levels(all_groups_comp))), cpoint = 2,
      fac = as.factor(all_groups_comp[all_groups_comp == pair_1_list[i] |
                                        all_groups_comp == pair_2_list[i]]),
      sub = paste("NMDS plot of Microbial Profiles\n ",pair_1_list[i]," - ",pair_2_list[i], "\n(p-value ",pVal[i],","," corr. p-value ", pVal_BH[i],")",sep="")
    )
  }
  dev.off()
  # Generated MDS plots are stored in one pdf file called "pairwise-beta-diversity-MDS.pdf"
  file_name <- paste(group_name,"pairwise-beta-diversity-MDS.pdf",sep="_")
  pdf(paste(group_name,"/",file_name,sep=""))
  for(i in 1:length(combn(unique_groups,2)[1,])){
    # Calculate and display the MDS plot (Multidimensional Scaling plot)
    s.class(
      cmdscale(pairedMatrixList[[i]], k = 2), col = rainbow(length(levels(all_groups_comp))), cpoint =
        2, fac = as.factor(all_groups_comp[all_groups_comp == pair_1_list[i] |
                                             all_groups_comp == pair_2_list[i]]), sub = paste("MDS plot of Microbial Profiles\n ",pair_1_list[i]," - ",pair_2_list[i], "\n(p-value ",pVal[i],","," corr. p-value ", pVal_BH[i],")",sep="")
    )
  }
  dev.off()
}
######                        Determine number of clusters                           ######
nclusters=NULL
if (dim(otu_file)[1]-1 <= kmers_limit) {
  kmers_limit=dim(otu_file)[1]-1
}
for (k in 1:kmers_limit) {
  if (k==1) {
    nclusters[k]=NA
  } else {
    # Partitioning the data into k clusters (max k is number of samples within the dataset)
    data_cluster=as.vector(pam(as.dist(unifract_dist_comp), k, diss=TRUE)$clustering)
    # Calculate Calinski-Harabasz Index
    nclusters[k]=calinhara(otu_file,data_cluster,k)
    print(k)
  }
}
# Generated plot showing the optimal number of clusters
pdf("de-novo-clustering.pdf")
plot(nclusters, type="h", xlab="k clusters", ylab="CH index",main="Optimal number of clusters")
dev.off()
#################################################################################
######                        Write Output Files                           ######
#################################################################################
# Write the distance matrix table in a file
file_name <- paste(group_name,"distance-matrix-gunif.tab",sep="_")
write.table( unifract_dist_comp, paste(group_name,"/",file_name,sep=""), sep = "\t", col.names = NA, quote = FALSE)
write.table( unifract_dist_comp, "distance-matrix-gunif.tab", sep = "\t", col.names = NA, quote = FALSE)
write.tree(tree,"samples-Tree.nwk",tree.names = FALSE)
# Graphical output files are generated in the main part of the script
if(!flag) { stop("
It was not possible to install all required R libraries properly.
Please check the installation of all required libraries manually.\n
Required libaries:ade4, GUniFrac, phangorn")
}
input_otu = "OTUs_Table-norm.tab"              #<--- CHANGE ACCORDINGLY !!!
#' Please give the name of the meta-file that contains individual sample information
input_meta = "mapping_file.txt"                #<--- CHANGE ACCORDINGLY !!!
#' Please give the name of the phylogenetic tree constructed from the OTU sequences
input_tree = "OTUs-Tree.tre"                   #<--- CHANGE ACCORDINGLY !!!
#' Please give the column name (in the mapping file) of the categorical variable to be used for comparison (e.g. Genotype)
group_name = "Variable"
label_samples = 1
label_id =c("")
kmers_limit=100
packages <-c("ade4","GUniFrac","phangorn","cluster","fpc")
# Function to check whether the package is installed
InsPack <- function(pack)
{
  if ((pack %in% installed.packages()) == FALSE) {
    install.packages(pack,repos ="http://cloud.r-project.org/")
  }
}
# Applying the installation on the list of packages
lapply(packages, InsPack)
# Make the libraries
lib <- lapply(packages, require, character.only = TRUE)
# Check if it was possible to install all required libraries
flag <- all(as.logical(lib))
otu_file <- read.table (file = input_otu, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "")
# Clean table from empty lines
otu_file <- otu_file[!apply(is.na(otu_file) | otu_file =="",1,all),]
# Load the mapping file containing individual sample information (sample names in the first column)
meta_file <- read.table (file = input_meta, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "")
# Save the column names of the mapping file
mappingVar <- names(meta_file)
# Clean table from empty lines
meta_file <- data.frame(meta_file[!apply(is.na(meta_file) | meta_file=="",1,all),],row.names=row.names(meta_file))
tree_file <- read.tree(input_tree)
# Create the directory where all output files are saved (is named after the target group name set above for comparisons)
dir.create(group_name)
otu_file <- otu_file[,order(names(otu_file))]
# Transpose OTU-table and convert format to a data frame
otu_file <- data.frame(t(otu_file))
# Order the mapping file by sample names (ascending)
meta_file <- data.frame(meta_file[order(row.names(meta_file)),],row.names=row.names(meta_file))
# Assign the column names to the modified mapping file
names(meta_file) <- mappingVar
# Save the position of the target group name in the mapping file
meta_file_pos <- which(colnames(meta_file) == group_name)
# Select metadata group based on the pre-set group name
all_groups <- as.factor(meta_file[,meta_file_pos])
# Root the OTU tree at midpoint
rooted_tree <- midpoint(tree_file)
absent <- rooted_tree$tip.label[!(rooted_tree$tip.label %in% colnames(otu_file))]
if (length(absent) != 0) {
  rooted_tree <- drop.tip(rooted_tree, absent)}
tip.label <- rooted_tree$tip.label
otu_file <- otu_file[, tip.label]
ntip <- length(tip.label)
nbr <- nrow(rooted_tree$edge)
edge <- rooted_tree$edge
edge2 <- edge[, 2]
br.len <- rooted_tree$edge.length
unifracs <- GUniFrac(otu_file, rooted_tree, alpha = c(0.0,0.5,1.0))$unifracs
# Weight on abundant lineages so the distance is not dominated by highly abundant lineages with 0.5 having the best power
unifract_dist <- unifracs[, , "d_0.5"]
all_dist_matrix <- as.dist(unifract_dist)
# Apply a hierarchical cluster analysis on the distance matrix based on the Ward's method
all_fit <- hclust(all_dist_matrix, method = "ward.D2")
# Generates a tree from the hierarchically generated object
tree <- as.phylo(all_fit)
my_tree_file_name <- paste(group_name,"/phylogram.pdf",sep="")
plot_color<-rainbow(length(levels(all_groups)))[all_groups]
# Save the generated phylogram in a pdf file
pdf(my_tree_file_name)
# The tree is visualized as a Phylogram color-coded by the selected group name
plot(tree, type = "phylogram",use.edge.length = TRUE, tip.color = (plot_color), label.offset = 0.01)
print.phylo(tree)
axisPhylo()
tiplabels(pch = 16, col = plot_color)
dev.off()
#################            Build NMDS plot           ########################
# Generated figures are saved in a pdf file
file_name <- paste(group_name,"beta-diversity.pdf",sep="_")
pdf(paste(group_name,"/",file_name,sep=""))
# Calculate the significance of variance to compare multivariate sample means (including two or more dependent variables)
# Omit cases where there isn't data for the sample (NA)
all_groups_comp <- all_groups[!is.na(all_groups)]
unifract_dist_comp <- unifract_dist[, !is.na(all_groups)]
adonis<-adonis(as.dist(unifract_dist_comp) ~ all_groups_comp)
all_groups_comp<-factor(all_groups_comp,levels(all_groups_comp)[unique(all_groups_comp)])
# Calculate and display the MDS plot (Multidimensional Scaling plot)
s.class(
  cmdscale(unifract_dist_comp, k = 2), col = unique(plot_color), cpoint =
    2, fac = all_groups_comp, sub = paste("MDS plot of Microbial Profiles\n(p-value ",adonis[[1]][6][[1]][1],")",sep="")
)
if (label_samples==1) {
  lab_samples <- row.names(cmdscale(unifract_dist_comp, k = 2))
  ifelse (label_id != "",lab_samples <- replace(lab_samples, !(lab_samples %in% label_id), ""), lab_samples)
  text(cmdscale(unifract_dist_comp, k = 2),labels=lab_samples,cex=0.7,adj=c(-.1,-.8))
}
# Calculate and display the NMDS plot (Non-metric Multidimensional Scaling plot)
meta <- metaMDS(unifract_dist_comp,k = 2)
s.class(
  meta$points, col = unique(plot_color), cpoint = 2, fac = all_groups_comp,
  sub = paste("metaNMDS plot of Microbial Profiles\n(p-value ",adonis[[1]][6][[1]][1],")",sep="")
)
if (label_samples==1){
  lab_samples <- row.names(meta$points)
  ifelse (label_id != "",lab_samples <- replace(lab_samples, !(lab_samples %in% label_id), ""), lab_samples)
  text(meta$points,labels=lab_samples,cex=0.7,adj=c(-.1,-.8))
}
#close the pdf file
dev.off()
###############          NMDS for pairwise analysis        ###################
# This plot is only generated if there are more than two groups included in the comparison
# Calculate the pairwise significance of variance for group pairs
# Get all groups contained in the mapping file
unique_groups <- levels(all_groups_comp)
if (dim(table(unique_groups)) > 2) {
  # Initialise vector and lists
  pVal = NULL
  pairedMatrixList <- list(NULL)
  pair_1_list <- NULL
  pair_2_list <- NULL
  for (i in 1:length(combn(unique_groups,2)[1,])) {
    # Combine all possible pairs of groups
    pair_1 <- combn(unique_groups,2)[1,i]
    pair_2 <- combn(unique_groups,2)[2,i]
    # Save pairs information in a vector
    pair_1_list[i] <- pair_1
    pair_2_list[i] <- pair_2
    # Generate a subset of all samples within the mapping file related to one of the two groups
    inc_groups <-
      rownames(subset(meta_file, meta_file[,meta_file_pos] == pair_1
                      |
                        meta_file[,meta_file_pos] == pair_2))
    # Convert UniFrac distance matrix to data frame
    paired_dist <- as.data.frame(unifract_dist_comp)
    # Save all row names of the mapping file
    row_names <- rownames(paired_dist)
    # Add row names to the distance matrix
    paired_dist <- cbind(row_names,paired_dist)
    # Generate distance matrix with samples of the compared groups (column-wise)
    paired_dist <- paired_dist[sapply(paired_dist[,1], function(x) all(x %in% inc_groups)),]
    # Remove first column with unnecessary group information
    paired_dist[,1] <- NULL
    paired_dist <- rbind(row_names,paired_dist)
    # Generate distance matrix with samples of the compared group (row-wise)
    paired_dist <- paired_dist[,sapply(paired_dist[1,], function(x) all(x %in% inc_groups))]
    # Remove first row with unnecessary group information
    paired_dist <- paired_dist[-1,]
    # Convert generated distance matrix to data type matrix (needed by multivariate analysis)
    paired_matrix <- as.matrix(paired_dist)
    class(paired_matrix) <- "numeric"
    # Save paired matrix in list
    pairedMatrixList[[i]] <- paired_matrix
    # Applies multivariate analysis to a pair out of the selected groups
    adonis <- adonis(paired_matrix ~ all_groups_comp[all_groups_comp == pair_1 |
                                                       all_groups_comp == pair_2])
    # List p-values
    pVal[i] <- adonis[[1]][6][[1]][1]
  }
  # Adjust p-values for multiple testing according to Benjamini-Hochberg method
  pVal_BH <- p.adjust(pVal,method="BH", n=length(pVal))
  # Generated NMDS plots are stored in one pdf file called "pairwise-beta-diversity-nMDS.pdf"
  file_name <- paste(group_name,"pairwise-beta-diversity.pdf",sep="_")
  pdf(paste(group_name,"/",file_name,sep=""))
  for(i in 1:length(combn(unique_groups,2)[1,])){
    meta <- metaMDS(pairedMatrixList[[i]], k = 2)
    s.class(
      meta$points,
      col = rainbow(length(levels(all_groups_comp))), cpoint = 2,
      fac = as.factor(all_groups_comp[all_groups_comp == pair_1_list[i] |
                                        all_groups_comp == pair_2_list[i]]),
      sub = paste("NMDS plot of Microbial Profiles\n ",pair_1_list[i]," - ",pair_2_list[i], "\n(p-value ",pVal[i],","," corr. p-value ", pVal_BH[i],")",sep="")
    )
  }
  dev.off()
  # Generated MDS plots are stored in one pdf file called "pairwise-beta-diversity-MDS.pdf"
  file_name <- paste(group_name,"pairwise-beta-diversity-MDS.pdf",sep="_")
  pdf(paste(group_name,"/",file_name,sep=""))
  for(i in 1:length(combn(unique_groups,2)[1,])){
    # Calculate and display the MDS plot (Multidimensional Scaling plot)
    s.class(
      cmdscale(pairedMatrixList[[i]], k = 2), col = rainbow(length(levels(all_groups_comp))), cpoint =
        2, fac = as.factor(all_groups_comp[all_groups_comp == pair_1_list[i] |
                                             all_groups_comp == pair_2_list[i]]), sub = paste("MDS plot of Microbial Profiles\n ",pair_1_list[i]," - ",pair_2_list[i], "\n(p-value ",pVal[i],","," corr. p-value ", pVal_BH[i],")",sep="")
    )
  }
  dev.off()
}
######                        Determine number of clusters                           ######
nclusters=NULL
if (dim(otu_file)[1]-1 <= kmers_limit) {
  kmers_limit=dim(otu_file)[1]-1
}
for (k in 1:kmers_limit) {
  if (k==1) {
    nclusters[k]=NA
  } else {
    # Partitioning the data into k clusters (max k is number of samples within the dataset)
    data_cluster=as.vector(pam(as.dist(unifract_dist_comp), k, diss=TRUE)$clustering)
    # Calculate Calinski-Harabasz Index
    nclusters[k]=calinhara(otu_file,data_cluster,k)
    print(k)
  }
}
# Generated plot showing the optimal number of clusters
pdf("de-novo-clustering.pdf")
plot(nclusters, type="h", xlab="k clusters", ylab="CH index",main="Optimal number of clusters")
dev.off()
#################################################################################
######                        Write Output Files                           ######
#################################################################################
# Write the distance matrix table in a file
file_name <- paste(group_name,"distance-matrix-gunif.tab",sep="_")
write.table( unifract_dist_comp, paste(group_name,"/",file_name,sep=""), sep = "\t", col.names = NA, quote = FALSE)
write.table( unifract_dist_comp, "distance-matrix-gunif.tab", sep = "\t", col.names = NA, quote = FALSE)
write.tree(tree,"samples-Tree.nwk",tree.names = FALSE)
# Graphical output files are generated in the main part of the script
if(!flag) { stop("
It was not possible to install all required R libraries properly.
Please check the installation of all required libraries manually.\n
Required libaries:ade4, GUniFrac, phangorn")
}
library(gplots)
distancia=as.matrix(all_dist_matrix)
heatmap.2(distancia, key=T, trace="none", ColSideColors = c(rep("aquamarine",2),rep("yellow",2),rep("blue",2),rep("red",3),rep("green",3),rep("darkgoldenrod3",2)), RowSideColors = c(rep("aquamarine",2),rep("yellow",2),rep("blue",2),rep("red",3),rep("green",3),rep("darkgoldenrod3",2)))
View(distancia)
View(distancia)
heatmap.2(distancia, key=T, trace="none", ColSideColors = c(rep("red",2),rep("green",3),rep("yellow",2),rep("aquamarine",2),rep("maroon1",2),rep("blue",2),rep("red",1)), RowSideColors = c(rep("red",2),rep("green",3),rep("yellow",2),rep("aquamarine",2),rep("maroon1",2),rep("blue",2),rep("red",1)))
