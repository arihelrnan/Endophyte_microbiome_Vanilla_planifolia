#Load programs
packages <-c("metacoder")
lib <- lapply(packages, require, character.only = TRUE)

#Load tab delimited file with OTU table information
otu_table <-  read.table ("../../Data/binary_table_16S.txt",
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
#
sample_table <- sample_table[-c(1,2,3,4,5,6,7,8),]
otu_table <- otu_table[,-c(1,2,3,4,5,6,7,8)]

#Make taxmap object
metamapa <- parse_tax_data(otu_table, #Make taxmap object.
                           class_cols = "taxonomy", # the column that contains taxonomic information
                           class_sep = ";", # The character used to separate taxa in the classification
                           class_regex = "^([a-z]{0,1})__{0,2}(.*)$", # Regex identifying where the data for each taxon is
                           class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))  
print(metamapa)

#Filter taxa with less 10 OTU
metamapa=filter_taxa(metamapa, n_obs> 10)

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

#Plot taxonomic data between organs
heat_tree(metamapa, 
          node_label = taxon_names,
          node_size = metamapa$data$tax_occ[['Root']],
          node_color = Root, 
          initial_layout = "re", layout = "da",
          title = "Root samples",
          node_label_size_range = c(0.01, 0.025),
          node_color_axis_label = "Mean of OTUs",
          node_size_axis_label = "Number of samples",
          output_file = "../../Figures/Root_samples_16S.png")

heat_tree(metamapa, 
          node_label = taxon_names,
          node_size = metamapa$data$tax_occ[['Stem']],
          node_color = Stem, 
          initial_layout = "re", layout = "da",
          title = "Stem samples",
          node_label_size_range = c(0.01, 0.025),
          node_color_axis_label = "Mean of OTUs",
          node_size_axis_label = "Number of samples",
          output_file = "../../Figures/Stem_samples_16S.png")

#Plot taxonomic data in Asymptomatic, Symptomatic and Wild
heat_tree(metamapa, 
          node_label = taxon_names,
          node_size = metamapa$data$tax_occ[['Symptomatic']],
          node_color =  Symptomatic, 
          initial_layout = "re", layout = "da",
          title = "Symptomatic",
          node_label_size_range = c(0.01, 0.025),
          node_color_axis_label = "Mean of OTUs",
          node_size_axis_label = "Number of samples",
          output_file = "../../Figures/Symptomatic_samples_16S.png")

heat_tree(metamapa, 
          node_label = taxon_names,
          node_size = metamapa$data$tax_occ[['Symptomatic']],
          node_color =  Symptomatic, 
          initial_layout = "re", layout = "da",
          title = "Symptomatic",
          node_label_size_range = c(0.01, 0.025),
          node_color_axis_label = "Mean of OTUs",
          node_size_axis_label = "Number of samples",
          output_file = "../../Figures/Symptomatic_samples_16S.png")

heat_tree(metamapa, 
          node_label = taxon_names,
          node_size = metamapa$data$tax_occ[['Wild']],
          node_color =  Wild, 
          initial_layout = "re", layout = "da",
          title = "Wild",
          node_label_size_range = c(0.01, 0.025),
          node_color_axis_label = "Mean of OTUs",
          node_size_axis_label = "Number of samples",
          output_file = "../../Figures/Wild_samples_16S.png")


#Comparing taxon abundance in variables groups
metamapa$data$diff_table <- compare_groups(metamapa, data = "tax_abund",
                                      cols = sample_table$`#SampleID`,
                                      groups = sample_table$State)
print(metamapa$data$diff_table)

metamapa <- mutate_obs(metamapa, "diff_table",
                  wilcox_p_value = p.adjust(wilcox_p_value, method = "BH"))
range(metamapa$data$diff_table$wilcox_p_value, finite = TRUE) 

metamapa$data$diff_table$log2_median_ratio[metamapa$data$diff_table$wilcox_p_value > 0.05] <- 0

heat_tree_matrix(metamapa,
                 data = "diff_table",
                 node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                 node_label = taxon_names,
                 node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                 node_color_range = diverging_palette(), # The built-in palette for diverging data
                 node_color_trans = "linear", # The default is scaled by circle area
                 node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                 edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "../../Figures/differential_heat_tree_ITS.png") # Saves the plot as a pdf file

heat_tree(metamapa, 
          node_label = taxon_names,
          node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
          node_color = log2_median_ratio, # A column from `obj$data$diff_table`
          node_color_interval = c(-2, 2), # The range of `log2_median_ratio` to display
          node_color_range = diverging_palette(), # The color palette used
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Log 2 ratio of median proportions",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations
