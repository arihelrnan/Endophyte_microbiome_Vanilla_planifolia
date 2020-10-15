library(metacoder)
library(dplyr)
library(ggplot2)
input_otu = "OTUs_Table-norm-reÃ±_tax.tab"
input_meta = "mapping_file.tab"
otu_file <- read.table (file = input_otu, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "")
input_otu = "OTUs_Table-norm-rel_tax.tab"
input_meta = "mapping_file.tab"
otu_file <- read.table (file = input_otu, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "")
input_otu = "OTUs_Table-norm-rel-tax.tab"
input_meta = "mapping_file.tab"
otu_file <- read.table (file = input_otu, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "")
otu_file <- otu_file[!apply(is.na(otu_file) | otu_file =="",1,all),]
View(otu_file)
meta_file <- read.table (file = input_meta, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "")
meta_file <- data.frame(meta_file[!apply(is.na(meta_file) | meta_file=="",1,all),])
View(otu_file)
View(meta_file)
sample_table <- read.table (file = input_meta, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "")
sample_table <- data.frame(meta_file[!apply(is.na(meta_file) | meta_file=="",1,all),])
View(sample_table)
metamapa$data$diff_table <- compare_groups(metamapa, dataset = "tax_abund",
                                           cols = sample_table$SampleID,
                                           groups = sample_table$Estado)
metamapa <- parse_tax_data(otu_table,
                           class_cols = "taxonomy",
                           class_sep = ";",
                           class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                           class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))
metamapa=filter_taxa(metamapa, n_obs> 50)
metamapa$data$tax_abund <- calc_taxon_abund(metamapa, "tax_data", cols=sample_table$SampleID)
metamapa$data$tax_occ <- calc_n_samples(metamapa, "tax_abund", groups = sample_table$Estado)
metamapa$data$diff_table <- compare_groups(metamapa, dataset = "tax_abund",
                                           cols = sample_table$SampleID,
                                           groups = sample_table$Estado)
metamapa <- mutate_obs(metamapa, "diff_table",
                       wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))
metamapa$data$diff_table$log2_median_ratio[metamapa$data$diff_table$wilcox_p_value > 0.05] <- 0
heat_tree_matrix(data = "diff_table",
                 node_label = cleaned_names,
                 node_size = n_obs, # number of OTUs
                 node_color = log2_median_ratio, # difference between groups
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3), # symmetric interval
                 edge_color_interval = c(-3, 3), # symmetric interval
                 node_color_range = diverging_palette(), # diverging colors
                 node_size_axis_label = "OTU count",
                 node_color_axis_label = "Log 2 ratio of median counts",
                 layout = "da", initial_layout = "re",
                 key_size = 0.67,
                 seed = 2)
heat_tree_matrix(metamapa, data = "diff_table",
                 node_label = cleaned_names,
                 node_size = n_obs, # number of OTUs
                 node_color = log2_median_ratio, # difference between groups
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3), # symmetric interval
                 edge_color_interval = c(-3, 3), # symmetric interval
                 node_color_range = diverging_palette(), # diverging colors
                 node_size_axis_label = "OTU count",
                 node_color_axis_label = "Log 2 ratio of median counts",
                 layout = "da", initial_layout = "re",
                 key_size = 0.67,
                 seed = 2)
heat_tree_matrix(metamapa, data = "diff_table",
                 node_size = n_obs, # number of OTUs
                 node_color = log2_median_ratio, # difference between groups
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3), # symmetric interval
                 edge_color_interval = c(-3, 3), # symmetric interval
                 node_color_range = diverging_palette(), # diverging colors
                 node_size_axis_label = "OTU count",
                 node_color_axis_label = "Log 2 ratio of median counts",
                 layout = "da", initial_layout = "re",
                 key_size = 0.67,
                 seed = 2)
heat_tree_matrix(metamapa,
                 dataset = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions")
metamapa %>%
  taxa::filter_taxa(taxon_ranks == "o", supertaxa = TRUE, reassign_obs = FALSE) %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  taxa::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$"), reassign_obs = FALSE) %>%
  heat_tree_matrix(data = "diff_table",
                   node_label = cleaned_names,
                   node_size = n_obs, # number of OTUs
                   node_color = log2_median_ratio, # difference between groups
                   node_color_trans = "linear",
                   node_color_interval = c(-3, 3), # symmetric interval
                   edge_color_interval = c(-3, 3), # symmetric interval
                   node_color_range = diverging_palette(), # diverging colors
                   node_size_axis_label = "OTU count",
                   node_color_axis_label = "Log 2 ratio of median counts",
                   layout = "da", initial_layout = "re",
                   key_size = 0.67,
                   seed = 2)
metamapa <- parse_tax_data(otu_table,
                           class_cols = "taxonomy",
                           class_sep = ";",
                           class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                           class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))
metamapa$data$tax_abund <- calc_taxon_abund(metamapa, "tax_data", cols=sample_table$SampleID)
metamapa$data$tax_occ <- calc_n_samples(metamapa, "tax_abund", groups = sample_table$Estado)
metamapa <- parse_tax_data(otu_table,
                           class_cols = "taxonomy",
                           class_sep = ";",
                           class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                           class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))
metamapa$data$tax_abund <- calc_taxon_abund(metamapa, "tax_data",
                                            cols = sample_table$SampleID,
                                            groups = sample_table$Estado)
metamapa$data$tax_occ <- calc_n_samples(metamapa, "tax_abund", groups = sample_table$Estado)
metamapa <- parse_tax_data(otu_table,
                           class_cols = "taxonomy",
                           class_sep = ";",
                           class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                           class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))
metamapa$data$tax_abund <- calc_taxon_abund(metamapa, "tax_data", cols=sample_table$SampleID)
metamapa$data$tax_occ <- calc_n_samples(metamapa, "tax_abund", groups = sample_table$Estado)
metamapa$data$diff_table <- compare_groups(metamapa, dataset = "tax_abund",
                                           cols = sample_table$SampleID,
                                           groups = sample_table$Estado)
metamapa <- mutate_obs(metamapa, "diff_table",
                       wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))
metamapa$data$diff_table$log2_median_ratio[metamapa$data$diff_table$wilcox_p_value > 0.05] <- 0
metamapa %>%
  taxa::filter_taxa(taxon_ranks == "o", supertaxa = TRUE, reassign_obs = FALSE) %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  taxa::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$"), reassign_obs = FALSE) %>%
  heat_tree_matrix(data = "diff_table",
                   node_label = cleaned_names,
                   node_size = n_obs, # number of OTUs
                   node_color = log2_median_ratio, # difference between groups
                   node_color_trans = "linear",
                   node_color_interval = c(-3, 3), # symmetric interval
                   edge_color_interval = c(-3, 3), # symmetric interval
                   node_color_range = diverging_palette(), # diverging colors
                   node_size_axis_label = "OTU count",
                   node_color_axis_label = "Log 2 ratio of median counts",
                   layout = "da", initial_layout = "re",
                   key_size = 0.67,
                   seed = 2)
metamapa <- parse_tax_data(otu_table,
                           class_cols = "taxonomy",
                           class_sep = ";",
                           class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                           class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))
metamapa$data$tax_abund <- calc_taxon_abund(metamapa, "tax_data", cols=sample_table$SampleID)
metamapa$data$tax_occ <- calc_n_samples(metamapa, "tax_abund", groups = sample_table$Estado)
metamapa$data$diff_table <- compare_groups(metamapa, dataset = "tax_abund",
                                           cols = sample_table$SampleID,
                                           groups = sample_table$Estado)
metamapa %>%
  taxa::filter_taxa(taxon_ranks == "o", supertaxa = TRUE, reassign_obs = FALSE) %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  taxa::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$"), reassign_obs = FALSE) %>%
  heat_tree_matrix(data = "diff_table",
                   node_label = cleaned_names,
                   node_size = n_obs, # number of OTUs
                   node_color = log2_median_ratio, # difference between groups
                   node_color_trans = "linear",
                   node_color_interval = c(-3, 3), # symmetric interval
                   edge_color_interval = c(-3, 3), # symmetric interval
                   node_color_range = diverging_palette(), # diverging colors
                   node_size_axis_label = "OTU count",
                   node_color_axis_label = "Log 2 ratio of median counts",
                   layout = "da", initial_layout = "re",
                   key_size = 0.67,
                   seed = 2)
metamapa %>%
  taxa::filter_taxa(taxon_ranks == "o", supertaxa = TRUE, reassign_obs = FALSE) %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  taxa::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$"), reassign_obs = FALSE) %>%
  heat_tree_matrix(data = "diff_table",
                   node_label = cleaned_names,
                   node_size = n_obs, # number of OTUs
                   node_color = log2_median_ratio, # difference between groups
                   node_color_trans = "linear",
                   node_color_interval = c(-3, 3), # symmetric interval
                   edge_color_interval = c(-3, 3), # symmetric interval
                   node_color_range = diverging_palette(), # diverging colors
                   node_size_axis_label = "OTU count",
                   node_color_axis_label = "Log 2 ratio of median counts",
                   layout = "da", initial_layout = "re",
                   key_size = 0.67,
                   seed = 2)
heat_tree_matrix(help)
heat_tree_matrix(?)
metamapa %>%
  taxa::filter_taxa(taxon_ranks == "o", supertaxa = TRUE, reassign_obs = FALSE) %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  taxa::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$"), reassign_obs = FALSE) %>%
  heat_tree_matrix(data = "diff_table",
                   node_label = cleaned_names,node_label_size_range = c(0.02, 0.04),
                   node_size = n_obs, # number of OTUs
                   node_color = log2_median_ratio, # difference between groups
                   node_color_trans = "linear",
                   node_color_interval = c(-3, 3), # symmetric interval
                   edge_color_interval = c(-3, 3), # symmetric interval
                   node_color_range = diverging_palette(), # diverging colors
                   node_size_axis_label = "OTU count",
                   node_color_axis_label = "Log 2 ratio of median counts",
                   layout = "da", initial_layout = "re",
                   key_size = 0.67,
                   seed = 2)
metamapa %>%
  taxa::filter_taxa(taxon_ranks == "o", supertaxa = TRUE, reassign_obs = FALSE) %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  taxa::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$"), reassign_obs = FALSE) %>%
  heat_tree_matrix(data = "diff_table",
                   node_label = cleaned_names,node_label_size_range = c(0.03, 0.04),
                   node_size = n_obs, # number of OTUs
                   node_color = log2_median_ratio, # difference between groups
                   node_color_trans = "linear",
                   node_color_interval = c(-3, 3), # symmetric interval
                   edge_color_interval = c(-3, 3), # symmetric interval
                   node_color_range = diverging_palette(), # diverging colors
                   node_size_axis_label = "OTU count",
                   node_color_axis_label = "Log 2 ratio of median counts",
                   layout = "da", initial_layout = "re",
                   key_size = 0.67,
                   seed = 2)
metamapa %>%
  taxa::filter_taxa(taxon_ranks == "o", supertaxa = TRUE, reassign_obs = FALSE) %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  taxa::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$"), reassign_obs = FALSE) %>%
  heat_tree_matrix(data = "diff_table",
                   node_label = cleaned_names,node_label_size_range = c(0.025, 0.04),
                   node_size = n_obs, # number of OTUs
                   node_color = log2_median_ratio, # difference between groups
                   node_color_trans = "linear",
                   node_color_interval = c(-3, 3), # symmetric interval
                   edge_color_interval = c(-3, 3), # symmetric interval
                   node_color_range = diverging_palette(), # diverging colors
                   node_size_axis_label = "OTU count",
                   node_color_axis_label = "Log 2 ratio of median counts",
                   layout = "da", initial_layout = "re",
                   key_size = 0.67,
                   seed = 2)
