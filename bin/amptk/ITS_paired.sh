#!/bin/bash
# -*- ENCODING: UTF-8 -*-

#Generate merge sequence. R1 and R2 sequence are contained in Illumina_ITS_data/ directory.  
amptk illumina -i Illumina_ITS_data/ -o ITS_merged/ -f ITS3 -r ITS4 -l 300 --min_len 150 --cleanup 

#OTU clustering using UPARSE with 97% of identy. I run this command of each sample, and I specify the name of output file. 
amptk cluster -i ITS_merged.demux.fq -o ITS_cluster --uchime_ref ITS -m 10

#OTU clustering using UPARSE with 97% of identy. I run this command of each sample, and I specify the name of output file. 
amptk filter -i ITS_cluster.otu_table.txt -f ITS_cluster.cluster.otus.fa -o ITS_filter -p 0.005 --min_reads_otu 10

#Assign taxonomy to OTUs
amptk taxonomy -f ITS_filter.filtered.otus.fa -i ITS_filter.final.txt -o ITS_taxonomy -m mapping_file.txt -d ITS2 --tax_filter Fungi(base) arihel@pavilion:/hdd/arihel/se
