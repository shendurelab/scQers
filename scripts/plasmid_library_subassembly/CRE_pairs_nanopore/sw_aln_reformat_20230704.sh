#!/usr/bin/env bash

fastq_input=$1
fasta_signpost=$2
out_aln_file=$3

# alignment
echo "performing sw alignment"
./ssw_test -r $fastq_input $fasta_signpost > temp_aln0.txt
echo "done sw alignment"

# reformat
echo "reformatting output file"
sed 'H;1h;$!d;x; s/\n\n/µ/g' temp_aln0.txt > temp_aln1.txt
tr "\n" "\t" < temp_aln1.txt > temp_aln2.txt
tr "µ" "\n" < temp_aln2.txt > temp_aln2_fmt.txt
rm temp_aln0.txt
rm temp_aln1.txt
rm temp_aln2.txt

#clean up in case no good alignment to ensure 7 fields per read
echo "cleaning up reads with no alignments"
awk -F'\t' 'NF==7 {print}' temp_aln2_fmt.txt  > $out_aln_file
rm temp_aln2_fmt.txt
