#Copy sequence fasta to Sequence_fasta directory
cp *.fna Sequence_fasta/

#Add Sample ID to fasta files according mapping file
add_qiime_labels.py -i etiquetas -m Mapilla.txt -c InputFileName -o Sequence_labelled/

