#Copy representative sequence of OTU of 16SrRNA data to this directory

cp V1_16S/rep_set.fna repset_V1_16S.fna
cp V2_16S/rep_set.fna repset_V2_16S.fna
cp V3_16S/rep_set.fna repset_V3_16S.fna
cp V4_16S/rep_set.fna repset_V4_16S.fna
cp V5_16S/rep_set.fna repset_V5_16S.fna
cp V6_16S/rep_set.fna repset_V6_16S.fna
cp V7_16S/rep_set.fna repset_V7_16S.fna
cp V8_16S/rep_set.fna repset_V8_16S.fna
cp V9_16S/rep_set.fna repset_V9_16S.fna
cp V10_16S/rep_set.fna repset_V10_16S.fna
cp V11_16S/rep_set.fna repset_V11_16S.fna
cp V12_16S/rep_set.fna repset_V12_16S.fna
cp V13_16S/rep_set.fna repset_V13_16S.fna
cp V14_16S/rep_set.fna repset_V14_16S.fna
cp V15_16S/rep_set.fna repset_V15_16S.fna
cp V16_16S/rep_set.fna repset_V16_16S.fna

#Copy representative sequence of OTU of ITS data to this directory

cp V1_ITS/rep_set.fna repset_V1_ITS.fna
cp V2_ITS/rep_set.fna repset_V2_ITS.fna
cp V3_ITS/rep_set.fna repset_V3_ITS.fna
cp V4_ITS/rep_set.fna repset_V4_ITS.fna
cp V5_ITS/rep_set.fna repset_V5_ITS.fna
cp V6_ITS/rep_set.fna repset_V6_ITS.fna
cp V7_ITS/rep_set.fna repset_V7_ITS.fna
cp V8_ITS/rep_set.fna repset_V8_ITS.fna
cp V9_ITS/rep_set.fna repset_V9_ITS.fna
cp V10_ITS/rep_set.fna repset_V10_ITS.fna
cp V11_ITS/rep_set.fna repset_V11_ITS.fna
cp V12_ITS/rep_set.fna repset_V12_ITS.fna
cp V13_ITS/rep_set.fna repset_V13_ITS.fna
cp V14_ITS/rep_set.fna repset_V14_ITS.fna
cp V15_ITS/rep_set.fna repset_V15_ITS.fna
cp V16_ITS/rep_set.fna repset_V16_ITS.fna

#Assign taxonomy to sequence OTU using default data base with the RDP classifier

assign_taxonomy.py -i repset_V1_16S.fna -m rdp -c 0.8 -o Taxonomy_V1_16S
assign_taxonomy.py -i repset_V2_16S.fna -m rdp -c 0.8 -o Taxonomy_V2_16S
assign_taxonomy.py -i repset_V3_16S.fna -m rdp -c 0.8 -o Taxonomy_V3_16S
assign_taxonomy.py -i repset_V4_16S.fna -m rdp -c 0.8 -o Taxonomy_V4_16S
assign_taxonomy.py -i repset_V5_16S.fna -m rdp -c 0.8 -o Taxonomy_V5_16S
assign_taxonomy.py -i repset_V6_16S.fna -m rdp -c 0.8 -o Taxonomy_V6_16S
assign_taxonomy.py -i repset_V7_16S.fna -m rdp -c 0.8 -o Taxonomy_V7_16S
assign_taxonomy.py -i repset_V8_16S.fna -m rdp -c 0.8 -o Taxonomy_V8_16S
assign_taxonomy.py -i repset_V9_16S.fna -m rdp -c 0.8 -o Taxonomy_V9_16S
assign_taxonomy.py -i repset_V10_16S.fna -m rdp -c 0.8 -o Taxonomy_V10_16S
assign_taxonomy.py -i repset_V11_16S.fna -m rdp -c 0.8 -o Taxonomy_V11_16S
assign_taxonomy.py -i repset_V12_16S.fna -m rdp -c 0.8 -o Taxonomy_V12_16S
assign_taxonomy.py -i repset_V13_16S.fna -m rdp -c 0.8 -o Taxonomy_V13_16S
assign_taxonomy.py -i repset_V14_16S.fna -m rdp -c 0.8 -o Taxonomy_V14_16S
assign_taxonomy.py -i repset_V15_16S.fna -m rdp -c 0.8 -o Taxonomy_V15_16S
assign_taxonomy.py -i repset_V16_16S.fna -m rdp -c 0.8 -o Taxonomy_V16_16S

#Assign taxonomy to sequence OTU of ITS using Unite data base with the RDP classifier

assign_taxonomy.py -i rep_set_V1_ITS.fna -t qiime3/sh_taxonomy_qiime_ver8_97_s_all_02.02.2019.txt -r qiime3/sh_refs_qiime_ver8_97_s_all_02.02.2019.fasta -m rdp -c 0.8 -o Taxonomy_V1_ITS --rdp_max_memory=20000
assign_taxonomy.py -i rep_set_V2_ITS.fna -t qiime3/sh_taxonomy_qiime_ver8_97_s_all_02.02.2019.txt -r qiime3/sh_refs_qiime_ver8_97_s_all_02.02.2019.fasta -m rdp -c 0.8 -o Taxonomy_V2_ITS --rdp_max_memory=20000
assign_taxonomy.py -i rep_set_V3_ITS.fna -t qiime3/sh_taxonomy_qiime_ver8_97_s_all_02.02.2019.txt -r qiime3/sh_refs_qiime_ver8_97_s_all_02.02.2019.fasta -m rdp -c 0.8 -o Taxonomy_V3_ITS --rdp_max_memory=20000
assign_taxonomy.py -i rep_set_V4_ITS.fna -t qiime3/sh_taxonomy_qiime_ver8_97_s_all_02.02.2019.txt -r qiime3/sh_refs_qiime_ver8_97_s_all_02.02.2019.fasta -m rdp -c 0.8 -o Taxonomy_V4_ITS --rdp_max_memory=20000
assign_taxonomy.py -i rep_set_V5_ITS.fna -t qiime3/sh_taxonomy_qiime_ver8_97_s_all_02.02.2019.txt -r qiime3/sh_refs_qiime_ver8_97_s_all_02.02.2019.fasta -m rdp -c 0.8 -o Taxonomy_V5_ITS --rdp_max_memory=20000
assign_taxonomy.py -i rep_set_V6_ITS.fna -t qiime3/sh_taxonomy_qiime_ver8_97_s_all_02.02.2019.txt -r qiime3/sh_refs_qiime_ver8_97_s_all_02.02.2019.fasta -m rdp -c 0.8 -o Taxonomy_V6_ITS --rdp_max_memory=20000
assign_taxonomy.py -i rep_set_V7_ITS.fna -t qiime3/sh_taxonomy_qiime_ver8_97_s_all_02.02.2019.txt -r qiime3/sh_refs_qiime_ver8_97_s_all_02.02.2019.fasta -m rdp -c 0.8 -o Taxonomy_V7_ITS --rdp_max_memory=20000
assign_taxonomy.py -i rep_set_V8_ITS.fna -t qiime3/sh_taxonomy_qiime_ver8_97_s_all_02.02.2019.txt -r qiime3/sh_refs_qiime_ver8_97_s_all_02.02.2019.fasta -m rdp -c 0.8 -o Taxonomy_V8_ITS --rdp_max_memory=20000
assign_taxonomy.py -i rep_set_V9_ITS.fna -t qiime3/sh_taxonomy_qiime_ver8_97_s_all_02.02.2019.txt -r qiime3/sh_refs_qiime_ver8_97_s_all_02.02.2019.fasta -m rdp -c 0.8 -o Taxonomy_V9_ITS --rdp_max_memory=20000
assign_taxonomy.py -i rep_set_V10_ITS.fna -t qiime3/sh_taxonomy_qiime_ver8_97_s_all_02.02.2019.txt -r qiime3/sh_refs_qiime_ver8_97_s_all_02.02.2019.fasta -m rdp -c 0.8 -o Taxonomy_V10_ITS --rdp_max_memory=20000
assign_taxonomy.py -i rep_set_V11_ITS.fna -t qiime3/sh_taxonomy_qiime_ver8_97_s_all_02.02.2019.txt -r qiime3/sh_refs_qiime_ver8_97_s_all_02.02.2019.fasta -m rdp -c 0.8 -o Taxonomy_V11_ITS --rdp_max_memory=20000
assign_taxonomy.py -i rep_set_V12_ITS.fna -t qiime3/sh_taxonomy_qiime_ver8_97_s_all_02.02.2019.txt -r qiime3/sh_refs_qiime_ver8_97_s_all_02.02.2019.fasta -m rdp -c 0.8 -o Taxonomy_V12_ITS --rdp_max_memory=20000
assign_taxonomy.py -i rep_set_V13_ITS.fna -t qiime3/sh_taxonomy_qiime_ver8_97_s_all_02.02.2019.txt -r qiime3/sh_refs_qiime_ver8_97_s_all_02.02.2019.fasta -m rdp -c 0.8 -o Taxonomy_V13_ITS --rdp_max_memory=20000
assign_taxonomy.py -i rep_set_V14_ITS.fna -t qiime3/sh_taxonomy_qiime_ver8_97_s_all_02.02.2019.txt -r qiime3/sh_refs_qiime_ver8_97_s_all_02.02.2019.fasta -m rdp -c 0.8 -o Taxonomy_V14_ITS --rdp_max_memory=20000
assign_taxonomy.py -i rep_set_V15_ITS.fna -t qiime3/sh_taxonomy_qiime_ver8_97_s_all_02.02.2019.txt -r qiime3/sh_refs_qiime_ver8_97_s_all_02.02.2019.fasta -m rdp -c 0.8 -o Taxonomy_V15_ITS --rdp_max_memory=20000
assign_taxonomy.py -i rep_set_V16_ITS.fna -t qiime3/sh_taxonomy_qiime_ver8_97_s_all_02.02.2019.txt -r qiime3/sh_refs_qiime_ver8_97_s_all_02.02.2019.fasta -m rdp -c 0.8 -o Taxonomy_V16_ITS --rdp_max_memory=20000
