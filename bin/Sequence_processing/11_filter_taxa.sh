#Remove OTU classify as Chloroplast from OTU table

filter_taxa_from_otu_table.py -i master_table_16S.biom -o master_table_16S_2 -n c__Chloroplast

#Remove OTU classify as Mitochondria from OTU table

filter_taxa_from_otu_table.py -i master_table_16S_2 -o Otu_table_16S.biom -n k__Mitochondria

#Remove OTU classify in the Kingdom Viridiplantae from OTU table

filter_taxa_from_otu_table.py -i master_table_ITS.biom -o master_table_ITS.biom_2 -n k__Viridiplantae

#Remove OTU classify as Chloroplast and Mitochondria from OTU table

filter_taxa_from_otu_table.py -i master_table_ITS.biom_2 -o Otu_table_ITS.biom -n k__Metazoa
 
#Convert OTU tables in format Biom to .txt

biom convert -i Otu_table_16S.biom -o Otu_table_16S.txt --to-tsv --header-key taxonomy --output-metadata-id "ConsensusLineage"
biom convert -i Otu_table_ITS.biom -o Otu_table_ITS.txt --to-tsv --header-key taxonomy --output-metadata-id "ConsensusLineage"
