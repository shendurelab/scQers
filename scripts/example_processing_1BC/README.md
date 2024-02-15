This is an example set of data and scripts to process simple one barcode data and obtain list of bona fide sequenced barcodes (for example characterizing a plasmid library for the mBC composition).

First, you can grab the sequences from the fastq file as follows (you will need to have the requisit dependencies installed, such as Biopython). 

```python3 parse_1read_fastq.py --in_R1 mBC_p29_R1.fastq.gz --R1_name BC --out_file p29_mBC_list.txt.gz```

This takes in a fastq, and output the mBC sequences as a compressed text file (p29_mBC_list.txt.gz).
These mBC sequences can then be counted (i.e., determining the number of reads per mBC) as follows. The function call also takes as input a threshold read count (here 300), outputs the mBCs above the threshold, and prints a figure of the count distribution:

```Rscript count_reads_per_mBC.R p29_mBC_list.txt.gz table_high_count_mBC_p29.txt 300 plot_dist_mBC_p29.pdf```

This repository contains the expected input and outputs of both commands. 
