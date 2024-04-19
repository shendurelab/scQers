rule get_barcodes:
    input:
        input_file = (config["output_directory"] + "/" + "{ID}" + "/" + "outs/possorted_genome_bam.bam")
    params:
        input_directory = (config["input_directory"] + "{ID}"),
        barcode = lambda wildcards: filepaths_dictionary[wildcards.ID]
    output:
        output_file = (config["output_directory"] + "/" + "{ID}" + "/outs/" + "{ID}" + "_get_bc_v3.txt")
    shell:""" 
 
    if [ {params.barcode} == "m" ]; then
        barcode_length=15
        search_seq="TCGACAA"
        barcode_type="m"
    fi

    if [ {params.barcode} == "o" ]; then
        barcode_length=16
        search_seq="GCTT"
        barcode_type="o"
    fi
    
    #barcode_length and search_seq in dataframe accessing
    
    module load python/3.7.8
    module load pysam

    echo "getting barcodes from cellRanger bam"
        
    python scripts/get_barcode_v2_fixed_pos_w_seq_check_20220201.py \
        --input_bam {input.input_file} \
        --output_file {output.output_file} \
        --barcode_length $barcode_length \
        --seq_start 0 \
        --chimeric_threshold 0.20 \
        --search_seq $search_seq \
        --barcode_type $barcode_type
        
        #add these to params
        

    echo "barcode retrieval complete"
    
    ls {output.output_file}
    """
