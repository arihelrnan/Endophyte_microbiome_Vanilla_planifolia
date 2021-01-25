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
#Make taxmap object
metamapa <- parse_tax_data(otu_table, #Make taxmap object.
                           class_cols = "taxonomy", # the column that contains taxonomic information
                           class_sep = ";", # The character used to separate taxa in the classification
                           class_regex = "^([a-z]{0,1})__{0,2}(.*)$", # Regex identifying where the data for each taxon is
                           class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))  
print(metamapa)

#Filter taxa with less 10 OTU
metamapa=filter_taxa(metamapa, n_obs> 10)

#Calculate OTUs in organ groups
metamapa$data$tax_abund <- calc_taxon_abund(metamapa, "tax_data")

metamapa$data$type_abund <- calc_group_mean(metamapa, "tax_abund",
                                            cols = sample_table$`#SampleID`,
                                            groups = sample_table$Organ)

metamapa$data$tax_occ <- calc_n_samples(metamapa, "tax_abund", groups = sample_table$Organ, cols = sample_table$`#SampleID`)

#Plot taxonomic data between organs
heat_tree(metamapa, 
          node_label = taxon_names,
          node_size = metamapa$data$tax_occ[['Root']],
          node_color = Root, 
          initial_layout = "re", layout = "da",
          title = "Root samples",
          node_label_size_range = c(0.01, 0.027),
          node_color_axis_label = "Mean of OTUs",
          node_size_axis_label = "Number of samples",
          output_file = "../../Figures/Taxonomic_plots/Root_samples_16S.png")

heat_tree(metamapa, 
          node_label = taxon_names,
          node_size = metamapa$data$tax_occ[['Stem']],
          node_color = Stem, 
          initial_layout = "re", layout = "da",
          title = "Stem samples",
          node_label_size_range = c(0.01, 0.027),
          node_color_axis_label = "Mean of OTUs",
          node_size_axis_label = "Number of samples",
          output_file = "../../Figures/Taxonomic_plots/Stem_samples_16S.png")

#Calculate OTUs in state groups
metamapa$data$tax_abund <- calc_taxon_abund(metamapa, "tax_data")

metamapa$data$type_abund <- calc_group_mean(metamapa, "tax_abund",
                                            cols = sample_table$`#SampleID`,
                                            groups = sample_table$State)

metamapa$data$tax_occ <- calc_n_samples(metamapa, "tax_abund", groups = sample_table$State, cols = sample_table$`#SampleID`)


#Plot taxonomic data in Asymptomatic, Symptomatic and Wild
heat_tree(metamapa, 
          node_label = taxon_names,
          node_size = metamapa$data$tax_occ[['Asymptomatic']],
          node_color =  Asymptomatic, 
          initial_layout = "re", layout = "da",
          title = "Asymptomatic samples",
          node_label_size_range = c(0.01, 0.027),
          node_color_axis_label = "Mean of OTUs",
          node_size_axis_label = "Number of samples",
          output_file = "../../Figures/Taxonomic_plots/Asymptomatic_samples_16S.png")

heat_tree(metamapa, 
          node_label = taxon_names,
          node_size = metamapa$data$tax_occ[['Symptomatic']],
          node_color =  Symptomatic, 
          initial_layout = "re", layout = "da",
          title = "Symptomatic samples",
          node_label_size_range = c(0.01, 0.027),
          node_color_axis_label = "Mean of OTUs",
          node_size_axis_label = "Number of samples",
          output_file = "../../Figures/Taxonomic_plots/Symptomatic_samples_16S.png")

heat_tree(metamapa, 
          node_label = taxon_names,
          node_size = metamapa$data$tax_occ[['Wild']],
          node_color =  Wild, 
          initial_layout = "re", layout = "da",
          title = "Wild samples",
          node_label_size_range = c(0.01, 0.027),
          node_color_axis_label = "Mean of OTUs",
          node_size_axis_label = "Number of samples",
          output_file = "../../Figures/Taxonomic_plots/Wild_samples_16S.png")

