# Organization
This directory content 3 scripts: 

**1. ITS_paired.sh.** This script run merged pair-end sequence for each sample using ITS data. 
  
**2. ITS_forward.sh.** This script use only forward sequence to the sequence processing using ITS data.

**3. 16S_paired.sh.**  This script run merged pair-end sequence for each sample using 16SrRNA data.

# Structure
Input sequence forqrd and reverse are in format fastq. AMPtk scripts consist in four steps: Pre-Processing, Clustering, OTU table filtering and Taxonomy. The parameters of each steps its variable acourding scripts: 

**Pre-Processing**

For `ITS_paired.sh` and `16S_paired.sh`, I use merged sequence, using option default of `--reads`. For script `ITS_forward-sh` I use only forward sequence indicate in option `--reads forward`. For ITS merged adn forward sequence, minimum length read to keep its 150pb and trim length in 300pb. For 16SrRNA sequence I trim length in a 450 pb, and min length is a 400 pb. 

**Clustering**

I use the uchime reference chimeras for filter according to sequence: If its 16SrRNA use `--uchime_ref 16S` and if its ITS use `--uchime_ref ITS`. The min size to select a OTU is 10 for two cases. 

**Filtering**

The filter index bleed between samples is `-p 0.005` and OTU with less 10 counts are filtering. 

**Assign taxonomy**

I install database for 16SrRNA and ITS with command `amptk install -i ITS 16S`, preciously to run scipts. I run script with each databases accourding sequence. I select only OTU classified in the kingdom Bacteria and Fungi. For 16SrRNA data, I will use QIIME for delete sequence classified as Mitchondria and Chloropast.  
