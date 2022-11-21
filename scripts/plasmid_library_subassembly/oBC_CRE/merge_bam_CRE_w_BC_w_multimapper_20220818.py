# # Opens sequencing files then outputs list of barocdes and ORFs
# # @author Diego 07/22/2021
# # modified by JB 04/03/2022: 
# added all in/out as paths as arguments
# pull additional information from bam file (alignment flag, positions-->serves as "UMI")
# changed to fixed sequence position search after barcode

from Bio import SeqIO
import sys, gzip, pysam, re, argparse


parser = argparse.ArgumentParser('.')
parser.add_argument('--input_bam', '-i', help='Position sorted BAM (or list of bams), from bowtie to list of regions.')
parser.add_argument('--input_BC_fastq', '-ibam', help='fastq.gz file corresponding to BC read.')
parser.add_argument('--output_file', '-o', help='Tab delimited file with cell, mutation barcode, read count, umi count. All observed barcodes correctable to a whitelist are reported.')
parser.add_argument('--barcode_start', type=int, help='start position of barcode in read')
parser.add_argument('--barcode_end', type=int, help='length of barcode')
parser.add_argument('--search_seq', help='Sequence to search for immediately after barcode.')

args = parser.parse_args()

input_bam=args.input_bam
input_BC_fastq=args.input_BC_fastq
output_file=args.output_file
barcode_end=args.barcode_end
search_seq=args.search_seq
barcode_start=args.barcode_start

search_seq_start = barcode_start + barcode_end + 1 
search_seq_end = barcode_start + barcode_end + 1 + len(search_seq)


print('reading in ' + input_BC_fastq)

# to deal with multi-mapper, easier to START by reading in the BC fastq and create a READID to BC dictionary. 
read_BC=dict()
with gzip.open(input_BC_fastq, 'rt') as BC_fq:
    bh = SeqIO.parse(BC_fq, 'fastq') 
    for read in bh:
        seq = read.seq.upper()
                
        #print(seq)
        #print(seq[search_seq_start: search_seq_end])
        
        # searching for the constant sequence after the BC for valid read
        is_ok=False
        if seq[search_seq_start: search_seq_end]==search_seq:
            is_ok=True
        
        if is_ok:
            barcode = seq[barcode_start: barcode_start+barcode_end+1]
            read_BC[read.id]=barcode
            #print(barcode)
            #print(read_alignments[read.id])
            #out.write("\t".join([str(read.id), str(barcode), str(read_alignments[read.id])]) + '\n')

print("done reading in BC fastq file")
  

print('reading in ' + input_bam)

# extracting some information from the bam file
with gzip.open(output_file, 'wt') as out:

    out.write("\t".join(['read_id', 'BC', 'CRE_id', 'mapq', 'CRE_read', 'read_forward_in_ref', 'read_start', 'read_end']) + '\n')
    
    #read_alignments=dict()
    with pysam.AlignmentFile(input_bam, "rb") as samfile:
        for read in samfile:
            # if read.query_name not in read_alignments: # REMOVE THIS AS THIS ELIMINATES MULTIMAPPERS
            if read.reference_name is not None and read.mapping_quality is not None and read.query_name in read_BC:                    
                summary="\t".join([read.reference_name, str(read.mapping_quality), str(read.query_sequence), str(read.is_reverse), str(read.reference_start), str(read.reference_end)])
                #read_alignments[read.query_name]=summary
                out.write("\t".join([str(read.query_name), str(read_BC[read.query_name]), str(summary)]) + '\n')
                    

# read through barcode fastq and pair with aligned reads
#with gzip.open(input_BC_fastq, 'rt') as BC_fq, \
#        gzip.open(output_file, 'wt') as out:
        
#    out.write("\t".join(['read_id', 'BC', 'CRE_id', 'mapq', 'CRE_read', 'read_forward_in_ref', 'read_start', 'read_end']) + '\n')
    
    
#    bh = SeqIO.parse(BC_fq, 'fastq') 
#    for read in bh:
#        seq = read.seq.upper()
                
#        print(seq)
#        print(seq[search_seq_start: search_seq_end])
#        # searching for the constant sequence after the BC for valid read
#        is_ok=False
#        if seq[search_seq_start: search_seq_end]==search_seq:
#            is_ok=True
        
#        if read.id in read_alignments and is_ok:
#            barcode = seq[barcode_start: barcode_start+barcode_end+1]
#            print(barcode)
#            print(read_alignments[read.id])
#            out.write("\t".join([str(read.id), str(barcode), str(read_alignments[read.id])]) + '\n')
                

