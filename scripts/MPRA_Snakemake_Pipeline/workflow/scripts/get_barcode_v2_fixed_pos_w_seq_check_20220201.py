
import argparse
import pysam
import sys
import collections


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Script to generate a file with counts of all mutation barcodes found within each cell in a CROP-seq or similar experiment.')
    parser.add_argument('--input_bam', '-i', help='Position sorted BAM (or list of bams) from 10X pipestance.')
    parser.add_argument('--output_file', '-o', help='Tab delimited file with cell, mutation barcode, read count, umi count. All observed barcodes correctable to a whitelist are reported.')
    parser.add_argument('--whitelist', required=False, default=None, help='Optional mutation barcode whitelist.')
    parser.add_argument('--seq_start',required=False, type=int, help='position of start of bc if known (python 0 indexing applies!)')
    parser.add_argument('--search_seq', required=False, help='Sequence to search for immediately upstream of the mutation barcode in transcript.')
    parser.add_argument('--barcode_length', required=False, type=int, help='If you are not providing a whitelist, you must provide the length of your barcodes in bp.')
    parser.add_argument('--chimeric_threshold', type=float, default=0.2, help='Threshold for calling a UMI non-chimeric.')
    parser.add_argument('--force_correction', type=int, help='Force correction to a specified edit distance. Mismatches that can map to different sequences will be ignored and left uncorrected.')
    parser.add_argument('--barcode_type', required=False, help='Barcode type for final output file')

    args = parser.parse_args()

    search_seq=args.search_seq
    barcode_start=args.seq_start
    barcode_length=args.barcode_length
    input_bam=args.input_bam
    chimeric_threshold=args.chimeric_threshold
    output_file_name=args.output_file
    barcode_type=args.barcode_type
    
    print("search_seq:", search_seq)
    print("barcode_start:", barcode_start)
    print("barcode_length:", barcode_length)
    print("input_bam:", input_bam)
    print("chimeric_threshold:", chimeric_threshold)
    print("output_file_name:", output_file_name)
    print("barcode_type:", barcode_type)
    

    
    ############################################################
    # tallying up UMIs observed in cells for each mBC
    ############################################################
    
    mutation_barcodes = {}
    
    read_number=0
    read_mapped=0
    read_no_corrected_cBC_umi=0
    read_no_search_seq_match=0
    
    print("read number:",read_number)
    print("mapped reads (not expected):",read_mapped)
    print("read_no_corrected_cBC_umi:", read_no_corrected_cBC_umi)
    print("read_no_search_seq_match:", read_no_search_seq_match)


    for read in pysam.Samfile(input_bam):
        read_number += 1
        #print("read_number:", read_number)

        if not read.is_unmapped:
            read_mapped += 1
            continue

        # get various part of the read in the bam file
        seq = read.seq.upper()
        tags = dict(read.tags)
        cell = tags.get('CB', None)
        umi = tags.get('UB', None)

        # skip read if no corrected UMI or cBC from cell ranger
        if not cell or not umi:
            read_no_corrected_cBC_umi += 1
            continue

        # only retain read if perfect match for sequence downstream of BC from constant region
        search_seq_start = barcode_start + barcode_length
        search_seq_end = barcode_start + barcode_length + len(search_seq)


        if seq[search_seq_start: search_seq_end]==search_seq:
            # based on the alignment to search sequence, extract barcode from the read in bam file print("sequence:", seq)
            barcode = seq[barcode_start: barcode_start+barcode_length]

            # leave barcode error correction out for now, can implement once we have accepted lists.
            corrected_barcode=barcode

            # list of all barcodes already listed w/ cBC=cell
            barcodes_in_cell = mutation_barcodes.get(cell, dict())

            # add the corresponding UMI to cBC
            if corrected_barcode not in barcodes_in_cell:
                barcodes_in_cell[corrected_barcode] = []
            barcodes_in_cell[corrected_barcode].append(umi)

            # update full dictionary with updated list of mBC/umis
            mutation_barcodes[cell] = barcodes_in_cell
            
        else:
            read_no_search_seq_match += 1
      
      
    print("read number:",read_number)
    print("mapped reads (not expected):",read_mapped)
    print("reads without corrected cBC+UMI:",read_no_corrected_cBC_umi)
    print("non mapped reads w/ corrected cBC+UMI:",read_number-read_mapped-read_no_corrected_cBC_umi)
    print("reads without exact match in search seq:", read_no_search_seq_match, read_no_search_seq_match/(read_number-read_mapped-read_no_corrected_cBC_umi))
    
    
    ############################################################
    # generate output file (filter out chimeric UMIs)
    ############################################################
    print("output_file_name:", output_file_name)
    barcode = barcode_type + "BC" 
    output_file=open(output_file_name, 'w')
    original_stdout = sys.stdout
    output_file.write('\t'.join(['cBC', barcode,'n_reads','n_UMI','n_reads_filtered','n_UMI_filtered','list_reads_per_UMIs','list_reads_per_UMIs_filtered'])+'\n')
    
    for cell in mutation_barcodes:

        # get the list of all UMI (each appearance corresponding to a read), irrespective of mBC, in a given cell (for a cBC)
        all_umis = []
        for barcode in mutation_barcodes[cell]:
            all_umis.extend(mutation_barcodes[cell][barcode])

        # tally up read counts for each UMI (again, irrespective of associated mBC)
        umi_counts_all_bc = collections.Counter(all_umis)    

        # loop through barcodes
        for barcode in mutation_barcodes[cell]:

            # counts of reads and UMIs for a given mBC, generate the list of read counts per UMI
            umi_counts_bc=collections.Counter(mutation_barcodes[cell][barcode])
            n_UMIs=len(umi_counts_bc.keys())
            n_reads=sum(umi_counts_bc.values())
            list_umi_counts=list(umi_counts_bc.items())
            str_list_umi_counts=[umi+"_"+str(count) for umi,count in list_umi_counts]

            # Get a list of UMIs that fall below the threshold for being likely non-chimeric.
            # the idea here: for a given cBC, each UMI should in principle be associated with a unique mBC. 
            # the proportion below is the fraction of reads from a UMI associated with the mBC of interest (over all other). If this fraction is too low
            chimeric_umis = set()
            for umi in umi_counts_bc:
                if float(umi_counts_bc[umi]) / umi_counts_all_bc[umi] < chimeric_threshold:
                    chimeric_umis.add(umi)

            # Filter any chimeric UMIs before outputing final stats
            filtered_umis = [umi for umi in mutation_barcodes[cell][barcode] if umi not in chimeric_umis]
            filtered_umi_counts_bc = collections.Counter(filtered_umis)
            n_UMIs_filtered=len(filtered_umi_counts_bc.keys())
            n_reads_filtered=sum(filtered_umi_counts_bc.values())
            list_umi_counts_filtered=list(filtered_umi_counts_bc.items())
            str_list_umi_counts_filtered=[umi+"_"+str(count) for umi,count in list_umi_counts_filtered]

            # output to file: only print if non 0 filtered umis>0. 
            if n_UMIs_filtered>0:
               output_file.write('\t'.join([cell,barcode,str(n_reads),str(n_UMIs),str(n_reads_filtered),str(n_UMIs_filtered)]))
               output_file.write('\t')
               sys.stdout = output_file
               print(*str_list_umi_counts,sep=",",end="\t")
               print(*str_list_umi_counts_filtered,sep=",")
               sys.stdout = original_stdout


        
        
        
        
