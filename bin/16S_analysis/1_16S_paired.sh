#!/bin/bash
# -*- ENCODING: UTF-8 -*-

#R1 and R2 sequence are contained in the illumina_16S_data/ directory. The output files are content in the 16S_merged/ directory. 
amptk illumina -i Illumina_16S_data/ -o 16S_merged/ -f CCTACGGGNGGCWGCAG -r GACTACHVGGGTATCTAATCC -l 300 --min_len 150 --cleanup 

#OTU clustering using UPARSE with 97% of identy. I run this command of each sample, and I specify the name of output file. 
amptk cluster -i 16S_merged.demux.fq -o 16S_cluster --uchime_ref 16S -m 10

#Filter OTU Table acourding 
amptk filter -i 16S_cluster.otu_table.txt -f 16S_cluster.cluster.otus.fa -o 16S_filter -p 0.005 --min_reads_otu 10

#Assign taxonomy to OTUs
amptk taxonomy -f 16S_filter.filtered.otus.fa -i 16S_filter.final.txt -o 16S_taxonomy -m 16S_merged.mapping_file.txt -d 16S --tax_filter Bacteria