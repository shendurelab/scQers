rule get_barcodes:
    input:
        input_file = (config["output_directory"] + "/" + "{ID}" + "/" + "outs/possorted_genome_bam.bam")
    params:
        barcode = lambda wildcards: samples_dictionary[wildcards.ID]
    output:
        output_file = (config["output_directory"] + "/" + "{ID}" + "/outs/" + "{ID}" + "_get_bc_v3.txt")
    conda:
        '/home/maurertm/micromamba/envs/scqers'
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

    echo "getting barcodes from cellRanger bam"
        
    get_barcode_fixed_pos \
        --input_bam {input.input_file} \
        --output_file {output.output_file} \
        --barcode_length $barcode_length \
        --seq_start 0 \
        --chimeric_threshold 0.20 \
        --search_seq $search_seq \
        --barcode_type $barcode_type

    #email about params- chimeric threshold in particular

    echo "barcode retrieval complete"
    
    #Snakemake will error out if you do not do something with the output file in the rule. 
    #To get around this, I ls the output file
    ls {output.output_file}

    """
