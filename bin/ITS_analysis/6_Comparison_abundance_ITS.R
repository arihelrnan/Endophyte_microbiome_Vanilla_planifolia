#Load programs
packages <-c("metacoder","tidyr")
lib <- lapply(packages, require, character.only = TRUE)

#Load tab delimited file with OTU table information
otu_table <-  read.table ("../../Data/ITS_taxonomy.otu_table.taxonomy.txt",
                          check.names = FALSE,
                          header = TRUE,
                          dec = ".",
                          sep = "\t",
                          row.names = 1,
                          comment.char = "")
#Load tab delimited file with information of each sample
sample_table <-  read.table ("../../Data/mapping_file.txt",
                             check.names = FALSE,
                             header = TRUE,
                             dec = ".",
                             sep = "\t",
                             comment.char = "")

#Change the taxonomy colum with taxonomic levels
otu_table = separate(otu_table, Taxonomy, into = c("Details", "taxonomy"), sep=";")
otu_table$Details=NULL

#Make taxmap object
metamapa <- parse_tax_data(otu_table, #Make taxmap object.
                           class_cols = "taxonomy", # the column that contains taxonomic information
                           class_sep = ",", # The character used to separate taxa in the classification
                           class_regex = "^([a-z]{0,1}):{0,2}(.*)$", # Regex identifying where the data for each taxon is
                           class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))  
print(metamapa)

#Filter taxa with less 2 OTU
metamapa=filter_taxa(metamapa, n_obs> 2)

#
metamapa$data$tax_abund <- calc_taxon_abund(metamapa, "tax_data")

metamapa$data$type_abund <- calc_group_mean(metamapa, "tax_abund",
                                            cols = sample_table$`#SampleID`,
                                            groups = sample_table$State)

#Removing low-abundance counts
metamapa$data$tax_data <- zero_low_counts(metamapa, dataset = "tax_data", min_count = 10)

#Comprobate OTU remove 
no_reads <- rowSums(metamapa$data$tax_data[, sample_table$`#SampleID`]) == 0
sum(no_reads)

#Remove OTU without reads
metamapa <- filter_obs(metamapa, target = "tax_data", ! no_reads, drop_taxa = TRUE)
print(metamapa)

#Normalize dividing each samples counts by the total number of counts observed for each sample
metamapa$data$tax_data <- calc_obs_props(metamapa, "tax_data")
print(metamapa)

#Calculate value of abundance per taxon
metamapa$data$tax_abund <- calc_taxon_abund(metamapa, "tax_data",
                                            cols = sample_table$`#SampleID`)

#Calculate number of samples for each taxon
metamapa$data$tax_occ <- calc_n_samples(metamapa, "tax_abund", groups = sample_table$State, cols = sample_table$`#SampleID`)

#Comparing taxon abundance in variables groups
metamapa$data$diff_table <- compare_groups(metamapa, data = "tax_abund",
                                           cols = sample_table$`#SampleID`,
                                           groups = paste(sample_table$State))


print(metamapa$data$diff_table)

heat_tree_matrix(metamapa,
                 data = "diff_table",
                 node_label = taxon_names,
                 node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                 node_color_range = diverging_palette(), # The built-in palette for diverging data
                 node_color_trans = "linear", # The default is scaled by circle area
                 node_color_interval = c(-10, 10), # The range of `log2_median_ratio` to display
                 edge_color_interval = c(-10, 10), # The range of `log2_median_ratio` to display
                 node_label_size_range = c(0.015, 0.03),
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "../../Figures/differential_heat_tree_ITS.png") # Saves the plot as a pdf file
