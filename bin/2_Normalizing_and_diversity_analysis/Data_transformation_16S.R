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

#########Save tables as data frame and convert in a .txt file#########
#Save binary table as .txt file
binary_table_16S <- as.data.frame(otu_table(binary_table))
binary_table_16S <- cbind(binary_table_16S, taxonomy_2)
write.table(binary_table_16S, "../../Data/binary_table_16S.txt", sep = "\t")
#Save relative table as .txt file
relative_table_16S <- as.data.frame(otu_table(phyloseq.rel))
relative_table_16S <- cbind(relative_table_16S, taxonomy_2)
write.table(relative_table_16S, "../../Data/relative_table_16S.txt", sep = "\t")
