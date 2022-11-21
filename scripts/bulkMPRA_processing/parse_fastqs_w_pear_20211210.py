
# inputs: 
# pear error corrected BC file and fastqs
# assumes R1 and R3 are RC of each other 


from Bio import SeqIO
import sys, gzip, argparse


# path to pear file
# e.g., p25_full_pear_BC_seqs.txt.gz
# see https://docs.python.org/3/library/argparse.html#module-argparse for argument parsing

parser = argparse.ArgumentParser('Script collating PEAR corrected barcodes with original reads for BC associations.')
parser.add_argument('--pear_file')
parser.add_argument('--out_valid')
parser.add_argument('--out_no_pear_BC')
parser.add_argument('--in_R1')
parser.add_argument('--in_R2')
parser.add_argument('--in_R3')
parser.add_argument('--R1_name')
parser.add_argument('--R2_name')
parser.add_argument('--R3_name')
parser.add_argument('--simple_out_bool',type=int)

args = parser.parse_args()
pear_BC_file = args.pear_file
OUT_file_valid = args.out_valid
OUT_no_pear_BC = args.out_no_pear_BC
IN_R1_file = args.in_R1
IN_R2_file = args.in_R2
IN_R3_file = args.in_R3
R1_name = args.R1_name
R2_name = args.R2_name
R3_name = args.R3_name
simple_out_bool = args.simple_out_bool

print(simple_out_bool)


# get error corrected bc from pear output (see e.g., seq_pear_20210504.py)
# pear_barcodes is an error corrected barcode for each read (or NA is no possible to error correct)
print('reading the pear barcodes in:',pear_BC_file)
pear_barcodes = {}
with gzip.open(pear_BC_file, 'rt') as bc: # reads are unique
    next(bc) # skips the header
    for line in bc:
        read, barcode = line.strip().split()
        pear_barcodes[read] = barcode

print('starting the processing.')

# open all the fastq files and output files printed to:
# examples
# R1: oBC
# R3: RC_oBC
# R2: mBC

#R1: mBC
#R2: UMI
#R3: mBC_RC


with gzip.open(IN_R1_file, 'rt') as R1_fq, \
        gzip.open(IN_R3_file, 'rt') as R3_R1_rev_fq, \
        gzip.open(IN_R2_file, 'rt') as R2_fq, \
        gzip.open(OUT_file_valid, 'wt') as out_valid, \
        gzip.open(OUT_no_pear_BC, 'wt') as out_no_pear_BC:

    # write header for output files and parse handles to get iterators
    if simple_out_bool:
        out_valid.write("\t".join([R2_name, R1_name+'_pear']) + '\n')
        out_no_pear_BC.write("\t".join([R2_name, R1_name+'_pear']) + '\n')
    else:
        out_valid.write("\t".join(['read_id', R2_name, R1_name, R3_name, R1_name+'_pear']) + '\n')
        out_no_pear_BC.write("\t".join(['read_id', R2_name, R1_name, R3_name, R1_name+'_pear']) + '\n')
    
    # handles to read files in
    h_R1_f, h_R3_R1_rc, h_R2 = [SeqIO.parse(e, 'fastq') for e in [R1_fq, R3_R1_rev_fq, R2_fq]]

    # loop through reads and write everything to files 
    for r_R1_f in h_R1_f:
    
        # reads for the other files
        r_R3_R1_rc, r_R2 = next(h_R3_R1_rc), next(h_R2)
        
        # looking for the read id in the pear output (because of error correction, cannot look for BC directly!)
        if r_R1_f.id in pear_barcodes:
            pear_barcode = pear_barcodes[r_R1_f.id]
        else:
            pear_barcode = 'none'
        

        if (pear_barcode!='none'):
            # w/ error corrected BC, output to main file
            if (simple_out_bool==1):
                out_valid.write("\t".join([str(r_R2.seq)] + \
                    [str(e) for e in [pear_barcode]]) + '\n')
            else:
                out_valid.write("\t".join([r_R1_f.id, str(r_R2.seq)] + \
                    [str(e) for e in [r_R1_f.seq, r_R3_R1_rc.seq, pear_barcode]]) + '\n')
        else:
            # w/o error corrected BC, output to no pear file
            if (simple_out_bool==1):
                out_no_pear_BC.write("\t".join([str(r_R2.seq)] + \
                    [str(e) for e in [pear_barcode]]) + '\n')
            else:
                out_no_pear_BC.write("\t".join([r_R1_f.id, str(r_R2.seq)] + \
                    [str(e) for e in [r_R1_f.seq, r_R3_R1_rc.seq, pear_barcode]]) + '\n')






