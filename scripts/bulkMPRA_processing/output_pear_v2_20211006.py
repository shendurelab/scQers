# opens sequencing files and then outputs a list error corrected barcodes

from Bio import SeqIO
import sys, gzip, getopt

# argument parsing based on: https://www.tutorialspoint.com/python/python_command_line_arguments.htm

argv=sys.argv[1:]
opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
for opt, arg in opts:
    if opt in ("-i", "--ifile"):
        pear_input_file = arg
    elif opt in ("-o", "--ofile"):
        output_file = arg

# open read files, loop through reads, and write to file
with gzip.open(pear_input_file, 'rt') as pear_file, \
	gzip.open(output_file, 'wt') as out:

	# write header for output file
	out.write("\t".join(['read', 'barcode']) + '\n')
	for r1 in SeqIO.parse(pear_file, 'fastq'):
		out.write("\t".join([r1.id, str(r1.seq)]) + '\n')
			
   


          