#Copy sequence with primers sequence remove to this directory
cp no_sqprimer_16S/V1_16S/reads.fastq V1_16S.fastq
cp no_sqprimer_16S/V2_16S/reads.fastq V2_16S.fastq
cp no_sqprimer_16S/V3_16S/reads.fastq V3_16S.fastq
cp no_sqprimer_16S/V4_16S/reads.fastq V4_16S.fastq
cp no_sqprimer_16S/V5_16S/reads.fastq V5_16S.fastq
cp no_sqprimer_16S/V6_16S/reads.fastq V5_16S.fastq
cp no_sqprimer_16S/V7_16S/reads.fastq V6_16S.fastq
cp no_sqprimer_16S/V8_16S/reads.fastq V7_16S.fastq
cp no_sqprimer_16S/V9_16S/reads.fastq V8_16S.fastq
cp no_sqprimer_16S/V10_16S/reads.fastq V10_16S.fastq
cp no_sqprimer_16S/V11_16S/reads.fastq V11_16S.fastq
cp no_sqprimer_16S/V12_16S/reads.fastq V12_16S.fastq
cp no_sqprimer_16S/V13_16S/reads.fastq V13_16S.fastq
cp no_sqprimer_16S/V14_16S/reads.fastq V14_16S.fastq
cp no_sqprimer_16S/V15_16S/reads.fastq V15_16S.fastq
cp no_sqprimer_16S/V16_16S/reads.fastq V16_16S.fastq

cp no_sqprimer_ITS/V1_ITS/reads.fastq V1_ITS.fastq
cp no_sqprimer_ITS/V2_ITS/reads.fastq V2_ITS.fastq
cp no_sqprimer_ITS/V3_ITS/reads.fastq V3_ITS.fastq
cp no_sqprimer_ITS/V4_ITS/reads.fastq V4_ITS.fastq
cp no_sqprimer_ITS/V5_ITS/reads.fastq V5_ITS.fastq
cp no_sqprimer_ITS/V6_ITS/reads.fastq V5_ITS.fastq
cp no_sqprimer_ITS/V7_ITS/reads.fastq V6_ITS.fastq
cp no_sqprimer_ITS/V8_ITS/reads.fastq V7_ITS.fastq
cp no_sqprimer_ITS/V9_ITS/reads.fastq V8_ITS.fastq
cp no_sqprimer_ITS/V10_ITS/reads.fastq V10_ITS.fastq
cp no_sqprimer_ITS/V11_ITS/reads.fastq V11_ITS.fastq
cp no_sqprimer_ITS/V12_ITS/reads.fastq V12_ITS.fastq
cp no_sqprimer_ITS/V13_ITS/reads.fastq V13_ITS.fastq
cp no_sqprimer_ITS/V14_ITS/reads.fastq V14_ITS.fastq
cp no_sqprimer_ITS/V15_ITS/reads.fastq V15_ITS.fastq
cp no_sqprimer_ITS/V16_ITS/reads.fastq V16_ITS.fastq



#16SrRna_data

usearch -fastq_filter V1_16S.fastq -fastq_maxee_rate 0.001 -fastaout V1_16S.fna 
usearch -fastq_filter V2_16S.fastq -fastq_maxee_rate 0.001 -fastaout V2_16S.fna 
usearch -fastq_filter V3_16S.fastq -fastq_maxee_rate 0.001 -fastaout V3_16S.fna 
usearch -fastq_filter V4_16S.fastq -fastq_maxee_rate 0.001 -fastaout V4_16S.fna 
usearch -fastq_filter V5_16S.fastq -fastq_maxee_rate 0.001 -fastaout V5_16S.fna 
usearch -fastq_filter V6_16S.fastq -fastq_maxee_rate 0.001 -fastaout V6_16S.fna 
usearch -fastq_filter V7_16S.fastq -fastq_maxee_rate 0.001 -fastaout V7_16S.fna 
usearch -fastq_filter V8_16S.fastq -fastq_maxee_rate 0.001 -fastaout V8_16S.fna 
usearch -fastq_filter V9_16S.fastq -fastq_maxee_rate 0.001 -fastaout V9_16S.fna 
usearch -fastq_filter V10_16S.fastq -fastq_maxee_rate 0.001 -fastaout V10_16S.fna 
usearch -fastq_filter V11_16S.fastq -fastq_maxee_rate 0.001 -fastaout V11_16S.fna 
usearch -fastq_filter V12_16S.fastq -fastq_maxee_rate 0.001 -fastaout V12_16S.fna 
usearch -fastq_filter V13_16S.fastq -fastq_maxee_rate 0.001 -fastaout V13_16S.fna 
usearch -fastq_filter V14_16S.fastq -fastq_maxee_rate 0.001 -fastaout V14_16S.fna 
usearch -fastq_filter V15_16S.fastq -fastq_maxee_rate 0.001 -fastaout V15_16S.fna 
usearch -fastq_filter V16_16S.fastq -fastq_maxee_rate 0.001 -fastaout V16_16S.fna 

#ITS_data

usearch -fastq_filter V1_ITS.fastq -fastq_maxee_rate 0.001 -fastaout V1_ITS.fna
usearch -fastq_filter V2_ITS.fastq -fastq_maxee_rate 0.001 -fastaout V2_ITS.fna
usearch -fastq_filter V3_ITS.fastq -fastq_maxee_rate 0.001 -fastaout V3_ITS.fna
usearch -fastq_filter V4_ITS.fastq -fastq_maxee_rate 0.001 -fastaout V4_ITS.fna
usearch -fastq_filter V5_ITS.fastq -fastq_maxee_rate 0.001 -fastaout V5_ITS.fna
usearch -fastq_filter V6_ITS.fastq -fastq_maxee_rate 0.001 -fastaout V6_ITS.fna
usearch -fastq_filter V7_ITS.fastq -fastq_maxee_rate 0.001 -fastaout V7_ITS.fna
usearch -fastq_filter V8_ITS.fastq -fastq_maxee_rate 0.001 -fastaout V8_ITS.fna
usearch -fastq_filter V9_ITS.fastq -fastq_maxee_rate 0.001 -fastaout V9_ITS.fna
usearch -fastq_filter V10_ITS.fastq -fastq_maxee_rate 0.001 -fastaout V10_ITS.fna
usearch -fastq_filter V11_ITS.fastq -fastq_maxee_rate 0.001 -fastaout V11_ITS.fna
usearch -fastq_filter V12_ITS.fastq -fastq_maxee_rate 0.001 -fastaout V12_ITS.fna
usearch -fastq_filter V13_ITS.fastq -fastq_maxee_rate 0.001 -fastaout V13_ITS.fna
usearch -fastq_filter V14_ITS.fastq -fastq_maxee_rate 0.001 -fastaout V14_ITS.fna
usearch -fastq_filter V15_ITS.fastq -fastq_maxee_rate 0.001 -fastaout V15_ITS.fna
usearch -fastq_filter V16_ITS.fastq -fastq_maxee_rate 0.001 -fastaout V16_ITS.fna
