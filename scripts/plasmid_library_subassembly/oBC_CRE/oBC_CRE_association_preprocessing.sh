#$ -q shendure-short.q
#$ -cwd
#$ -S /bin/bash
#$ -o /net/shendure/vol1/home/lalannej/sge_logs
#$ -e /net/shendure/vol1/home/lalannej/sge_logs
#$ -M lalannej@uw.edu
#$ -l mfree=8G,h_rt=6:00:00:00
#$ -pe serial 2
#$ -tc 2
#$ -t 1-2:1


source ~/virt_env/bin/activate

module load modules modules-init modules-gs
module load python/3.7.7
module load pysam

module load tbb/2020_U2 
module load bowtie2/2.4.4
module load samtools/1.14

module load gcc/8.1.0
module load R/3.5.1


# need to get bowtie index file from fasta reference prior to that, here:
# bowtie2-build  reference_files/devCRE_w_plasmid_and_promoters_20220404.fa bwt_idx/devCRE_PCR2_products_v2



# list of relevant in/out directories
fastq_dir="path_fastq/"
bowtie_idx="bwt_idx/devCRE_PCR2_products_v2"
align_out_dir="alignment_outs/"
merged_BC_CRE_out_dir="merge_BC_CRE_outs/"
pileup_outs_dir="pileup_outs/"

LOOKUP_FILE="input_file_names.txt"


# section for arrayed job lookup file parsing

set -e
if [[ ! -r "${LOOKUP_FILE}" ]]; then
    echo "Cannot find ${LOOKUP_FILE}" >&2
    exit 1
fi

LOOKUP="$(awk -v SGE_TASK_ID="${SGE_TASK_ID}" '$1 == SGE_TASK_ID {print $2}' < "${LOOKUP_FILE}")"
if [[ "x${LOOKUP}" = "x" ]]; then
    echo "Task ${SGE_TASK_ID} failed to lookup task ID" >&2
    exit 1
fi

date_str=$(date "+20%y%m%d")



# # # # # # # # # # # # # # # # # # 
# aligning with bowtie
# # # # # # # # # # # # # # # # # # 

out_unaligned=${align_out_dir}${LOOKUP}"_unaligned_k2_"${date_str}".fastq.gz"
out_sam=${align_out_dir}${LOOKUP}"_k2_"${date_str}".sam"

bowtie2 --threads 8 -k 2 \
	-x ${bowtie_idx} \
	-U ${fastq_dir}${LOOKUP}_R1_001.fastq.gz \
 	--un-gz ${out_unaligned} \
 	-S ${out_sam}



# # # # # # # # # # # # # # # # # # 
# convert to bam, order and index
# # # # # # # # # # # # # # # # # # 

out_bam=${align_out_dir}${LOOKUP}"_k2_"${date_str}".bam"
samtools view -hb -S ${out_sam} | samtools sort - -o ${out_bam}
samtools index ${out_bam}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# parsing the bam and fastq files to combine the alignment info with 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

input_bam_file=${out_bam}
input_BC_fastq=${fastq_dir}${LOOKUP}"_R3_001.fastq.gz"
output_merged_file_name=${merged_BC_CRE_out_dir}${LOOKUP}"_merged_k2_"${date_str}".txt.gz"

python /net/shendure/vol10/projects/JBL/seq023_devCRE_suba_sci_mBC/nobackup/suba_devCRE/merge_bam_CRE_w_BC_w_multimapper_20220818.py \
    --input_bam ${input_bam_file} \
    --input_BC_fastq ${input_BC_fastq} \
    --output_file ${output_merged_file_name} \
    --barcode_start 0 \
    --barcode_end 15 \
    --search_seq GCTTTA
  # 0 to 15 is 16 long here.
  
    
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# pile up of counts and positions across BC and CREs 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

input_merged_file=${output_merged_file_name}
out_pileup=$pileup_outs_dir${LOOKUP}"_pileup_k2_"${date_str}".txt"
out_pileup_full=$pileup_outs_dir${LOOKUP}"_full_pileup_k2_"${date_str}".txt"

Rscript --vanilla /net/shendure/vol10/projects/JBL/seq023_devCRE_suba_sci_mBC/nobackup/suba_devCRE/pileup_BC_CRE_20220404.R \
	$input_merged_file \
	$out_pileup \
	$out_pileup_full
	
gzip $out_pileup
gzip $out_pileup_full