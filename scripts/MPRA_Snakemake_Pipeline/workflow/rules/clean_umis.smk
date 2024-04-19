rule clean_umis:
    input:
        input_file = (config["output_directory"] + "/" + "{ID}" + "/outs/" + "{ID}" + "_get_bc_v3.txt")
    params:
        input_directory = (config["input_directory"] + "{ID}"),
        barcode = lambda wildcards: filepaths_dictionary[wildcards.ID]
    output:
        output_file = (config["output_directory"] + "/" + "{ID}" + "/outs/" + "{ID}" + "_get_bc_v3_no_G_cleaned_UMI.txt"),
        input_file_no_G = (config["output_directory"] + "/" + "{ID}" + "/outs/" + "{ID}" + "_get_bc_v3_no_G.txt")
    shell:"""
    echo "cleaning up UMIs"
    module load R/4.0

    if [ {params.barcode} == "o" ]; then
        awk '$2!="GGGGGGGGGGGGGGGG" {{ print $0 }}' {input.input_file} > {output.input_file_no_G}
        barcode="oBC"
    fi
    #email and find out if extra G is an accident
    if [ {params.barcode} == "m" ]; then
        awk '$2!="GGGGGGGGGGGGGGG" {{ print $0 }}' {input.input_file} > {output.input_file_no_G}
        barcode="mBC"
    fi
    
    # cleaning up the UMIs
    echo "cleaning UMIs"
    #google --vanilla
    Rscript --vanilla scripts/clean_up_UMI_counts_v3_20220126.R \
	    {output.input_file_no_G} \
	    {output.output_file} \
        ${{barcode}}  

    echo "UMI cleaning complete"

    ls {output.output_file}
    """