This is an example set of data and scripts to process oBC-CRE amplicons (structure of the amplicon [here](https://github.com/shendurelab/scQers/blob/main/custom_amplicon_structures/oBC_CRE_subassembly.gbk) and obtain list of bona fide paired oBC-BC sequenced barcodes.

You will need a working installation of bowtie2 (v2.4.4 here), samtools (v1.14 here), and pysam python module. 

Since the structure of the read is the CRE (from tagmented end), the first step is alignment of the reads to a reference set. 

The bowtie index (for the expected sequences in the library) can be generated with:

```bowtie2-build example_CRE.fasta example_CRE_idx```

Then, the sequences from read1 can be aligned to the CREs: 

```bowtie2 --threads 8 -k 2 -x example_CRE_idx -U example_CRE_oBC__R1_001.fastq.gz --un-gz unaligned_example_CRE_oBC__R1_001.fastq.gz -S aln_example_CRE_oBC__R1_001.sam```

Followed compression and index: 
```
samtools view -hb -S aln_example_CRE_oBC__R1_001.sam | samtools sort - -o aln_example_CRE_oBC__R1_001.bam
samtools index aln_example_CRE_oBC__R1_001.bam
```

The reads aligned to the expected CREs can then be merged with the barcodes using the script below, which also includes a sequence check downstream of the barcode (possible if your read is longer then the barcode length and extends in the constant region). 
```
python merge_bam_CRE_w_BC_w_multimapper_20220818.py \
    --input_bam aln_example_CRE_oBC__R1_001.bam \
    --input_BC_fastq example_CRE_oBC__R3_001.fastq.gz \
    --output_file example_CRE_oBC_merged.txt.gz \
    --barcode_start 0 \
    --barcode_end 15 \
    --search_seq GCT
```

The list of counted oBC-CRE pairs can be obtained as follows: 
```
Rscript pileup_BC_CRE_20220404.R \
	example_CRE_oBC_merged.txt.gz \
	example_CRE_oBC_pileup.txt \
	example_CRE_oBC_full_pileup.txt
```

Finally, the final list of bona fide oBC-CRE pairs is obtained by running (note that there are internal parameters in the script that relate to things such as size selection during library preparation that can be modulated according to the specific application, the script also takes in a read count threshold, dependent on the level of sequencing saturation, here threshold of 3 counts). The script also takes a boolean for plotting some QC figures. 
```
Rscript pileup_BC_CRE_20240209.R example_CRE_oBC_full_pileup.txt.gz example_CRE_oBC 3 p025_recloned_complex_mBC_oBC_subassembly.txt.gz TRUE
```

This repository contains the expected input and outputs of all above commands as examples.  
