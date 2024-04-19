rule cellranger:
    input:
        input_file = config["input_directory"] + "{ID}" + "/" + "{ID}" + "_S1_L001_R1_001.fastq" 
    resources:
        mem = "64G", 
        time = "2:00:00"
    params:
        reference_data = config["reference_data"],
        prefix = "{ID}",
        input_directory = config["input_directory"] + "{ID}" + "/",
        output_directory = config["output_directory"],
        output_directory_ids = config["output_directory"] + "/" + "{ID}"
    output:
        output_file = (config["output_directory"] + "/" + "{ID}" + "/" + "outs/possorted_genome_bam.bam"),
    shell:"""
        rm -rf {params.output_directory_ids}
        #could do different directory for the outputs so you don't have to remove
        
        module load cellranger
                
        code_directory=$(pwd)

        cd "{params.output_directory}" 

        echo "running cell ranger"
        cellranger count --id={params.prefix} \
             --transcriptome={params.reference_data} \
             --fastqs={params.input_directory} \
             --localmem=64 \
             --localcores=8
        #add 64 and 8 into the params

        cd $code_directory

        echo "cell ranger complete"
        
        ls {output.output_file}
        """