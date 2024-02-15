

from Bio import SeqIO
import sys, gzip, argparse


# path to pear file
# e.g., p25_full_pear_BC_seqs.txt.gz
# see https://docs.python.org/3/library/argparse.html#module-argparse for argument parsing

parser = argparse.ArgumentParser('Script collating PEAR corrected barcodes with original reads for BC associations.')
parser.add_argument('--in_R1')
parser.add_argument('--R1_name')
parser.add_argument('--out_file')


args = parser.parse_args()
IN_R1_file = args.in_R1
R1_name = args.R1_name
OUT_file = args.out_file


# open all  fastq file and output read content

with gzip.open(IN_R1_file, 'rt') as R1_fq, \
         gzip.open(OUT_file, 'wt') as out_file:

    # write header for output file
    out_file.write(R1_name + '\n')
    
    # handles to read files in
    h_R1_f = SeqIO.parse(R1_fq, 'fastq')

    # loop through reads and write everything to files 
    for r_R1_f in h_R1_f:
        out_file.write(str(r_R1_f.seq) + '\n')