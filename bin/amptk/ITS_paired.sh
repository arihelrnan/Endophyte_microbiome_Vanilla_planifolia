#!/bin/bash
# -*- ENCODING: UTF-8 -*-

#R1 and R2 sequence are contained in the illumina_ITS_data/ directory. The output files are content in the ITS_merged/ directory. 
amptk illumina -i illumina_ITS_data/ -o ITS_merged/ -f ITS3 -r ITS4 -l 300 --min_len 150 --cleanup 

#Make new folder in to ITS_merged/ directory and copy .demux files to the new folder
cd ITS_merged
mkdir cluster_data
cp *.demux.fq cluster_data/ 
cd cluster_data/

#OTU clustering using UPARSE with 97% of identy. I run this command of each sample, and I specify the name of output file. 
amptk cluster -i V1.demux.fq -o V1_cluster --uchime_ref ITS -m 10
amptk cluster -i V2.demux.fq -o V2_cluster --uchime_ref ITS -m 10
amptk cluster -i V3.demux.fq -o V3_cluster --uchime_ref ITS -m 10
amptk cluster -i V4.demux.fq -o V4_cluster --uchime_ref ITS -m 10
amptk cluster -i V5.demux.fq -o V5_cluster --uchime_ref ITS -m 10
amptk cluster -i V6.demux.fq -o V6_cluster --uchime_ref ITS -m 10
amptk cluster -i V7.demux.fq -o V7_cluster --uchime_ref ITS -m 10
amptk cluster -i V8.demux.fq -o V8_cluster --uchime_ref ITS -m 10
amptk cluster -i V9.demux.fq -o V9_cluster --uchime_ref ITS -m 10
amptk cluster -i V10.demux.fq -o V10_cluster --uchime_ref ITS -m 10
amptk cluster -i V11.demux.fq -o V11_cluster --uchime_ref ITS -m 10
amptk cluster -i V12.demux.fq -o V12_cluster --uchime_ref ITS -m 10
amptk cluster -i V13.demux.fq -o V13_cluster --uchime_ref ITS -m 10
amptk cluster -i V14.demux.fq -o V14_cluster --uchime_ref ITS -m 10
amptk cluster -i V15.demux.fq -o V15_cluster --uchime_ref ITS -m 10
amptk cluster -i V16.demux.fq -o V16_cluster --uchime_ref ITS -m 10

#Move directory cluster_data/ with OTU Tables and fasta files
cd ..
mv ITS_merged/cluster_data cluster_data
mkdir filter
cp cluster_data/*.txt filter/
cp cluster_data/*.fa filter/
cp cluster_data/*.log filter/
cd filter

#Filter OTU Table acourding 
amptk filter -i V1_cluster.otu_table.txt -f V1_cluster.cluster.otus.fa -o V1_filter -p 0.005 --min_reads_otu 10
amptk filter -i V2_cluster.otu_table.txt -f V2_cluster.cluster.otus.fa -o V2_filter -p 0.005 --min_reads_otu 10
amptk filter -i V3_cluster.otu_table.txt -f V3_cluster.cluster.otus.fa -o V3_filter -p 0.005 --min_reads_otu 10
amptk filter -i V4_cluster.otu_table.txt -f V4_cluster.cluster.otus.fa -o V4_filter -p 0.005 --min_reads_otu 10
amptk filter -i V5_cluster.otu_table.txt -f V5_cluster.cluster.otus.fa -o V5_filter -p 0.005 --min_reads_otu 10
amptk filter -i V6_cluster.otu_table.txt -f V6_cluster.cluster.otus.fa -o V6_filter -p 0.005 --min_reads_otu 10
amptk filter -i V7_cluster.otu_table.txt -f V7_cluster.cluster.otus.fa -o V7_filter -p 0.005 --min_reads_otu 10
amptk filter -i V8_cluster.otu_table.txt -f V8_cluster.cluster.otus.fa -o V8_filter -p 0.005 --min_reads_otu 10
amptk filter -i V9_cluster.otu_table.txt -f V9_cluster.cluster.otus.fa -o V9_filter -p 0.005 --min_reads_otu 10
amptk filter -i V10_cluster.otu_table.txt -f V10_cluster.cluster.otus.fa -o V10_filter -p 0.005 --min_reads_otu 10
amptk filter -i V11_cluster.otu_table.txt -f V11_cluster.cluster.otus.fa -o V11_filter -p 0.005 --min_reads_otu 10
amptk filter -i V12_cluster.otu_table.txt -f V12_cluster.cluster.otus.fa -o V12_filter -p 0.005 --min_reads_otu 10
amptk filter -i V13_cluster.otu_table.txt -f V13_cluster.cluster.otus.fa -o V13_filter -p 0.005 --min_reads_otu 10
amptk filter -i V14_cluster.otu_table.txt -f V14_cluster.cluster.otus.fa -o V14_filter -p 0.005 --min_reads_otu 10
amptk filter -i V15_cluster.otu_table.txt -f V15_cluster.cluster.otus.fa -o V15_filter -p 0.005 --min_reads_otu 10
amptk filter -i V16_cluster.otu_table.txt -f V16_cluster.cluster.otus.fa -o V16_filter -p 0.005 --min_reads_otu 10

#Make taxonomy/ directory and copy OTU tables and fasta files filtered
cd .. 
mkdir taxonomy
cp filter/*.filtered.otus.fa taxonomy/
cp filter/*.final.txt taxonomy/
cp ITS_merged.mapping_file.txt taxonomy/
cd taxonomy

#Assign taxonomy to OTUs
amptk taxonomy -f V1_filter.filtered.otus.fa -i V1_filter.final.txt -o V1_taxonomy -m ITS_merged.mapping_file.txt -d ITS2 --tax_filter Fungi
amptk taxonomy -f V2_filter.filtered.otus.fa -i V2_filter.final.txt -o V2_taxonomy -m ITS_merged.mapping_file.txt -d ITS2 --tax_filter Fungi
amptk taxonomy -f V3_filter.filtered.otus.fa -i V3_filter.final.txt -o V3_taxonomy -m ITS_merged.mapping_file.txt -d ITS2 --tax_filter Fungi
amptk taxonomy -f V4_filter.filtered.otus.fa -i V4_filter.final.txt -o V4_taxonomy -m ITS_merged.mapping_file.txt -d ITS2 --tax_filter Fungi
amptk taxonomy -f V5_filter.filtered.otus.fa -i V5_filter.final.txt -o V5_taxonomy -m ITS_merged.mapping_file.txt -d ITS2 --tax_filter Fungi
amptk taxonomy -f V6_filter.filtered.otus.fa -i V6_filter.final.txt -o V6_taxonomy -m ITS_merged.mapping_file.txt -d ITS2 --tax_filter Fungi
amptk taxonomy -f V7_filter.filtered.otus.fa -i V7_filter.final.txt -o V7_taxonomy -m ITS_merged.mapping_file.txt -d ITS2 --tax_filter Fungi
amptk taxonomy -f V8_filter.filtered.otus.fa -i V8_filter.final.txt -o V8_taxonomy -m ITS_merged.mapping_file.txt -d ITS2 --tax_filter Fungi
amptk taxonomy -f V9_filter.filtered.otus.fa -i V9_filter.final.txt -o V9_taxonomy -m ITS_merged.mapping_file.txt -d ITS2 --tax_filter Fungi
amptk taxonomy -f V10_filter.filtered.otus.fa -i V10_filter.final.txt -o V10_taxonomy -m ITS_merged.mapping_file.txt -d ITS2 --tax_filter Fungi
amptk taxonomy -f V11_filter.filtered.otus.fa -i V11_filter.final.txt -o V11_taxonomy -m ITS_merged.mapping_file.txt -d ITS2 --tax_filter Fungi
amptk taxonomy -f V12_filter.filtered.otus.fa -i V12_filter.final.txt -o V12_taxonomy -m ITS_merged.mapping_file.txt -d ITS2 --tax_filter Fungi
amptk taxonomy -f V13_filter.filtered.otus.fa -i V13_filter.final.txt -o V13_taxonomy -m ITS_merged.mapping_file.txt -d ITS2 --tax_filter Fungi
amptk taxonomy -f V14_filter.filtered.otus.fa -i V14_filter.final.txt -o V14_taxonomy -m ITS_merged.mapping_file.txt -d ITS2 --tax_filter Fungi
amptk taxonomy -f V15_filter.filtered.otus.fa -i V15_filter.final.txt -o V15_taxonomy -m ITS_merged.mapping_file.txt -d ITS2 --tax_filter Fungi
amptk taxonomy -f V16_filter.filtered.otus.fa -i V16_filter.final.txt -o V16_taxonomy -m ITS_merged.mapping_file.txt -d ITS2 --tax_filter Fungi
