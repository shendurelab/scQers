rule cellranger:
    input:
        input_file = config["input_directory"] + "/" + "{ID}" + "/" + "{ID}" + "_S1_L001_R1_001.fastq" 
    resources:
        mem = "256G", 
        time = "48:00:00"
    params:
        reference_data = config["reference_data"],
        prefix = "{ID}",
        input_directory = config["input_directory"] + "/" + "{ID}" + "/",
        output_directory = config["output_directory"],
        output_directory_ids = config["output_directory"] + "/" + "{ID}"
    output:
        output_file = (config["output_directory"] + "/" + "{ID}" + "/outs/possorted_genome_bam.bam"),
    shell:"""
        
        #Snakemake automatically makes missing output directories. CellRanger, however, errors out if it does not make the output directory.
        #Thus, the best way to get around this is to remove the output directory Snakemake automatically makes before running CellRanger.
        rm -rf {params.output_directory_ids}
        
        #CellRanger cannot be put into a conda environment. Thus, it must be module loaded before running. 
        module load cellranger

        #CellRanger does not allow you to specify an output directory. It just automatically otuputs in the directory it is located in.
        #Thus, in order to specify an output directory, you must (1) save the code directory, (2) cd to the output directory of your choice
        #(3) run CellRanger and (4) cd back into the code directory
        
        code_directory=$(pwd)

        cd "{params.output_directory}" 

        echo "running CellRanger"
        cellranger count --id={params.prefix} \
             --transcriptome={params.reference_data} \
             --fastqs={params.input_directory} \
             --localmem=64 \
             --localcores=8
        #add 64 and 8 into the params

        cd $code_directory

        echo "CellRanger complete"
        
        #Snakemake will error out if you do not do something with the output file in the rule. 
        #To get around this, I ls the output file
        ls {output.output_file}

        """
