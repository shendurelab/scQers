#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --job-name=single_cell_mpra

module load R/4.0
module load cellranger


input_directory=$1 
output_directory=$2 
barcode_type=$3 
reference_data=$4 

parent_directory=$(dirname $input_directory)
id_file="${parent_directory}/directories"

echo "input_directory: $input_directory"
echo "output_directory: $output_directory"
echo "barcode_type: $barcode_type"
echo "reference_data: $reference_data"
echo "id file: $id_file"

line_number=$SLURM_ARRAY_TASK_ID #get index of which file to process from $SLURM_ARRAY_TASK_ID provided by SLURM
echo "line number: $line_number"
ID="$(basename $(sed "${line_number}q; d" "${id_file}"))"
#ID=`awk NR=${line_number} $id_file`
#ID="$(sed "${line_number}q; d" "${id_file}")" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID
echo "ID: $ID"

input_directory=$input_directory"/"$ID"/"
output_directory=$output_directory"/"$ID"/"

if [[ ${output_directory} != "none" ]] && [[ ! -d ${output_directory} ]]; then
        mkdir -p ${output_directory}
fi

#LOOKUP="SRR22253239down"
#look up to make sure is m barcode
#dir_get_bc_out="get_barcode_mBC/"


# # # # # # # # # # # # # # # # # # # # 
# Calling cellRanger 
# # # # # # # # # # # # # # # # # # # # 

echo "running cellRanger"
#path_ori_fastq="/oak/stanford/groups/smontgom/maurertm/MPRA/"

#path_fastq=${path_ori_fastq}${LOOKUP}
#date_str=$(date "+20%y%m%d")


#path_cellRanger_out="/net/shendure/vol1/home/lalannej/vol10_projects_JB/seq016_scMPRA_w_TAP_miniPilot_crisprQTL/nobackup/CR_out_oBC/"${LOOKUP}"_oBC_CRv6_"${date_str}"/"
#echo $path_cellRanger_out
#output="/oak/stanford/groups/smontgom/maurertm/MPRA/scQers/scripts/scQer_analysis/pre_processing/SRR22253239down_CRv6_20240411/"

code_directory=$(pwd)
echo $code_directory

output_directory_dirname="$(dirname "$output_directory")"
output_directory_basename="$(basename "$output_directory")"
echo "$output_directory_basename"
echo "$output_directory_dirname"

cd $output_directory_dirname

file_path=${output_directory}${output_directory_basename}".mri.tgz"
echo "$file_path"

if [[ ! -f ${file_path} ]]
then
    echo "running cell ranger"
    cellranger count --id=$output_directory_basename \
                 --transcriptome=$reference_data \
                 --fastqs=$input_directory \
                 --localmem=64 \
                 --localcores=8
else
    echo "cell ranger run previously completed"
fi

cd $code_directory

# # # # # # # # # # # # # # # # # # # # 
# get barcode (processing of cellRanger bam)
# # # # # # # # # # # # # # # # # # # # 
echo "get_barcode from cellRanger bam"

#source ~/virt_env/bin/activate
#module load modules modules-init modules-gs
module load python/3.7.8
module load pysam



#output_file1_name=${dir_get_bc_out}${LOOKUP}"_get_bc_v3_"${date_str}".txt"


# parsing the cellRanger bam file

output_outs=${output_directory}"/outs/"
input_bam_file=${output_outs}"possorted_genome_bam.bam"
output_file1_name=${output_outs}${ID}"_get_bc_v3_"${date_str}".txt"

if [[ $barcode_type == "m" ]]
then
    barcode_length=15
    search_seq="TCGACAA"
elif [[ $barcode_type == "o" ]]
then
    barcode_length=16
    search_seq="GCTT"
else
    echo "barcode type not recognized"
    exit 1
fi
    

if [[ ! -f ${output_file1_name} ]]
then
   if [[ -f ${input_bam_file} ]]
   then
       echo "input exists"
    fi 
    python get_barcode_v2_fixed_pos_w_seq_check_20220201.py \
        --input_bam ${input_bam_file} \
        --output_file ${output_file1_name} \
        --barcode_length ${barcode_length} \
        --seq_start 0 \
        --chimeric_threshold 0.20 \
        --search_seq ${search_seq}
else
    echo "parsing the cellRanger bam file previously completed"
fi


# # # # # # # # # # # # # # # # # # # # 
# clean up UMI (Hamming distance correction)
# # # # # # # # # # # # # # # # # # # # 
echo "cleaning up UMIs"

# removing the G only mBCs (spurious and slow down UMI clean up)
#directory="/oak/stanford/groups/smontgom/maurertm/MPRA/scQers/scripts/scQer_analysis/pre_processing/SRR22253239down_CRv6_20240411/outs/"
output_file1_name_no_G=${output_outs}${ID}"_get_bc_v3_no_G_"${date_str}".txt"

if [[ $barcode_type == "m" ]]
then
    awk '$2!="GGGGGGGGGGGGGGG" { print $0 }' ${output_file1_name} > ${output_file1_name_no_G}
else
    awk '$2!="GGGGGGGGGGGGGGGG" { print $0 }' ${output_file1_name} > ${output_file1_name_no_G}
fi

# cleaning up the UMIs
#module load gcc/8.1.0
#module load R/3.5.3

output_file2_name=${output_outs}${ID}"_get_bc_v3_no_G_cleaned_UMI_"${date_str}".txt"

if [[ ! -f $output_file2_name ]]
then
    Rscript --vanilla clean_up_UMI_counts_v3_20220126.R \
	    ${output_file1_name_no_G} \
	    ${output_file2_name} \
	    mBC
else
    echo "UMIs have already been cleaned"
fi


