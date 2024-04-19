#$ -q shendure-short.q
#$ -cwd
#$ -S /bin/bash
#$ -o /net/shendure/vol1/home/lalannej/sge_logs
#$ -e /net/shendure/vol1/home/lalannej/sge_logs
#$ -m ae
#$ -l mfree=75G,h_rt=6:00:00:00
#$ -pe serial 1
#$ -tc 4
#$ -t 1-4:1

# -M lalannej@uw.edu


module load modules modules-init modules-gs
module load pcre2/10.35
module load hdf5/1.10.1
module load R/4.0.0
module load gcc/8.2.0
export PATH=/net/shendure/vol10/projects/Samuel/nobackup/10X/cellranger-6.0.1:$PATH
module load seqtk/1.3

set -e

LOOKUP_FILE="/net/shendure/vol10/projects/JBL/seq030_bulk_MPRA_mEB_v2__oBC_10x_scREA_reseq/nobackup/idx1_dmux/lookup_oBC_seq028_seq030_20220604.txt"


path_ori_fastq="/net/shendure/vol10/projects/JBL/seq030_bulk_MPRA_mEB_v2__oBC_10x_scREA_reseq/nobackup/idx1_dmux/"
dir_get_bc_out="/net/shendure/vol10/projects/JBL/seq030_bulk_MPRA_mEB_v2__oBC_10x_scREA_reseq/nobackup/idx1_dmux/get_barcode_oBC/"



if [[ ! -r "${LOOKUP_FILE}" ]]; then
    echo "Cannot find ${LOOKUP_FILE}" >&2
    exit 1
fi

LOOKUP="$(awk -v SGE_TASK_ID="${SGE_TASK_ID}" '$1 == SGE_TASK_ID {print $2}' < "${LOOKUP_FILE}")"

if [[ "x${LOOKUP}" = "x" ]]; then
    echo "Task ${SGE_TASK_ID} failed to lookup task ID" >&2
    exit 1
fi



# # # # # # # # # # # # # # # # # # # # 
# Calling cellRanger 
# # # # # # # # # # # # # # # # # # # # 

echo "running cellRanger"

path_fastq=${path_ori_fastq}${LOOKUP}
date_str=$(date "+20%y%m%d")


cellranger count --id=${LOOKUP}"_CRv6_"${date_str} \
                 --transcriptome=/net/shendure/vol10/projects/Samuel/nobackup/10X/refdata-cellranger-mm10-3.0.0 \
                 --fastqs=$path_fastq \
                 --localmem=64 \
                 --localcores=8



# # # # # # # # # # # # # # # # # # # # 
# get barcode (processing of cellRanger bam)
# # # # # # # # # # # # # # # # # # # # 

echo "get_barcode from cellRanger bam"

source ~/virt_env/bin/activate
module load modules modules-init modules-gs
module load python/3.7.7
module load pysam


input_bam_file=${LOOKUP}"_CRv6_"${date_str}"/outs/possorted_genome_bam.bam"

output_file1_name=${dir_get_bc_out}${LOOKUP}"_get_bc_v3_"${date_str}".txt"

# parsing the cellRanger bam file
python /net/shendure/vol1/home/lalannej/vol10_projects_JB/seq016_scMPRA_w_TAP_miniPilot_crisprQTL/nobackup/get_barcode_v2_fixed_pos_w_seq_check_20220201.py \
    --input_bam ${input_bam_file}\
    --output_file ${output_file1_name} \
    --barcode_length 16 \
    --seq_start 0 \
    --chimeric_threshold 0.20 \
    --search_seq GCTT



# # # # # # # # # # # # # # # # # # # # 
# clean up UMI (Hamming distance correction)
# # # # # # # # # # # # # # # # # # # # 
echo "cleaning up UMIs"

# removing the G only mBCs (spurious and slow down UMI clean up)
output_file1_name_no_G=${dir_get_bc_out}${LOOKUP}"_get_bc_v3_no_G_"${date_str}".txt"
awk '$2!="GGGGGGGGGGGGGGGG" { print $0 }' ${output_file1_name} > ${output_file1_name_no_G}

# cleaning up the UMIs
module load gcc/8.1.0
module load R/3.5.1
output_file2_name=${dir_get_bc_out}${LOOKUP}"_get_bc_v3_no_G_cleaned_UMI_"${date_str}".txt"
Rscript --vanilla /net/shendure/vol1/home/lalannej/vol10_projects_JB/seq015_mEB_crisprQTL_miniPilot_GFP_amplicons/nobackup/demux_fastq/clean_up_UMI_counts_v3_20220126.R \
	${output_file1_name_no_G} \
	${output_file2_name} \
	BC


