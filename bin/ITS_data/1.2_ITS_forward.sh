#!/bin/bash
# -*- ENCODING: UTF-8 -*-

#R1 and R2 sequence are contained in the illumina_ITS_data/ directory. The output files with only forwrad sequence are content in the ITS_forward/ directory. 
amptk illumina -i Illumina_ITS_data/ -o Fwd_ITS/ -f ITS3 -r ITS4 -l 300 --reads forward --min_len 150 --cleanup 

#OTU clustering using UPARSE with 97% of identy. I run this command of each sample, and I specify the name of output file.
amptk cluster -i Fwd_ITS.demux.fq -o Fwd_ITS_cluster --uchime_ref ITS -m 10

#Filter OTU Table acourding 
amptk filter -i Fwd_ITS_cluster.otu_table.txt -f Fwd_ITS_cluster.cluster.otus.fa -o Fwd_filter -p 0.005 --min_reads_otu 10

#Assign taxonomy to OTUs
amptk taxonomy -f Fwd_filter.filtered.otus.fa -i Fwd_filter.final.txt -o Fwd_taxonomy -m mapping_file.txt -d ITS2 --tax_filter Fungi
