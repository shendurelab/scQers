This is an example set of data and scripts to process two barcodes data (for example oBC-mBC paired amplicon, such as [this amplicon](https://github.com/shendurelab/scQers/blob/main/custom_amplicon_structures/oBC_mBC_subassembly.gbk) and obtain list of bona fide pairs of barcodes (for example characterizing a plasmid library for the oBC-mBC composition).

First, you can grab the sequences from the fastq files as follows (you will need to have the requisit dependencies installed, such as Biopython). 
```
python3 parse_combine_two_fastqs.py \
	--out_file p25_paired_oBC_mBC.txt.gz \
	--in_R1 head_p25_S1_R1.fastq.gz \
	--in_R2 head_p25_S1_R2.fastq.gz \
	--R1_name RC_oBC \
	--R2_name mBC	
```

This takes in two fastqs (typically from paired end sequencing), and output the two BC sequences as a compressed text file (p25_paired_oBC_mBC.txt.gz).

These oBC-mBC pairs of sequences can then be counted (i.e., determining the number of reads per pair) as follows. The function call also takes as input a threshold read count (here 300), outputs the pairs of BCs above the count threshold (hi_counts_p25_paired_oBC_mBC.txt.gz), and prints a figure of the count distribution:
```
Rscript count_2BC_file.R \
	p25_paired_oBC_mBC.txt.gz \
	all_counts_p25_paired_oBC_mBC.txt.gz \
	hi_counts_p25_paired_oBC_mBC.txt.gz \
	RC_oBC \
	mBC \
	plot_dist_p25_paired_oBC_mBC.pdf \
	300
```

This repository contains the expected input and outputs of both commands. 
