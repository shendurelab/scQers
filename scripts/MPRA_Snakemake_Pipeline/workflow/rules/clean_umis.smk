rule clean_umis:
    input:
        input_file = (config["output_directory"] + "/" + "{ID}" + "/outs/" + "{ID}" + "_get_bc_v3.txt")
    params:
        barcode = lambda wildcards: samples_dictionary[wildcards.ID]
    output:
        output_file = (config["output_directory"] + "/" + "{ID}" + "/outs/" + "{ID}" + "_get_bc_v3_no_G_cleaned_UMI.txt"),
        input_file_no_G = (config["output_directory"] + "/" + "{ID}" + "/outs/" + "{ID}" + "_get_bc_v3_no_G.txt")
    conda:
        '/home/maurertm/micromamba/envs/scqers'
    shell:"""
    
    if [ {params.barcode} == "o" ]; then
        awk '$2!="GGGGGGGGGGGGGGGG" {{ print $0 }}' {input.input_file} > {output.input_file_no_G}
        barcode="oBC"
    fi
    if [ {params.barcode} == "m" ]; then
        awk '$2!="GGGGGGGGGGGGGGG" {{ print $0 }}' {input.input_file} > {output.input_file_no_G}
        barcode="mBC"
    fi

    #email about G's
    
    echo "cleaning UMIs"
    clean_up_umi_counts --file_name {output.input_file_no_G} --file_name_corrected {output.output_file} --variable_name ${{barcode}}  

    echo "UMI cleaning complete"

    #Snakemake will error out if you do not do something with the output file in the rule. 
    #To get around this, I ls the output file
    ls {output.output_file}

    """
