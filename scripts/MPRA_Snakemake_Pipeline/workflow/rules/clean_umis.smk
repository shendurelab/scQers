rule clean_umis:
    input:
        input_file = (config["output_directory"] + "/" + "{ID}" + "/outs/" + "{ID}" + "_get_bc_v3.txt")
    params:
        barcode = lambda wildcards: filepaths_dictionary[wildcards.ID]
    output:
        output_file = (config["output_directory"] + "/" + "{ID}" + "/outs/" + "{ID}" + "_get_bc_v3_no_G_cleaned_UMI.txt"),
        input_file_no_G = (config["output_directory"] + "/" + "{ID}" + "/outs/" + "{ID}" + "_get_bc_v3_no_G.txt")
    conda:
        'maggie_python'
    shell:"""
    
    if [ {params.barcode} == "o" ]; then
        awk '$2!="GGGGGGGGGGGGGGGG" {{ print $0 }}' {input.input_file} > {output.input_file_no_G}
        barcode="oBC"
    fi
    if [ {params.barcode} == "m" ]; then
        awk '$2!="GGGGGGGGGGGGGGG" {{ print $0 }}' {input.input_file} > {output.input_file_no_G}
        barcode="mBC"
    fi
    
    echo "cleaning UMIs"
    Rscript --vanilla scripts/clean_up_UMI_counts_v3_20220126.R --file_name {output.input_file_no_G} --file_name_corrected {output.output_file} --variable_name ${{barcode}}  

    echo "UMI cleaning complete"

    #Snakemake will error out if you do not do something with the output file in the rule. 
    #To get around this, I ls the output file
    ls {output.output_file}

    """
