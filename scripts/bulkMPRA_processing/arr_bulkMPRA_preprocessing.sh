#$ -q shendure-short.q
#$ -cwd
#$ -S /bin/bash
#$ -o /net/shendure/vol1/home/lalannej/sge_logs
#$ -e /net/shendure/vol1/home/lalannej/sge_logs
#$ -m ae
#$ -l mfree=15G,h_rt=6:00:00:00
#$ -pe serial 1
#$ -tc 16
#$ -t 1-32:1

# -M lalannej@uw.edu

module load modules modules-init modules-gs
module load samtools/1.10
module load bwa/0.7.17
module load bedtools/2.29.2
module load gcc/8.1.0
module load R/3.5.1
module load python/3.6.5

module load pear/0.9.11 
module load seqtk/1.3

set -e

head_dir="data_directory/"

LOOKUP_FILE=${head_dir}"lookup_fastq_file_names.txt"
dir_pear_outs=${head_dir}"pear_outs/"
dir_merged_file=${head_dir}"merged_files/"
dir_final_outs=${head_dir}"BC_quant/"

fastq_dir=${head_dir}"idx2_demux/"


# from seq030: 
#R1: 18 cycles o369: mBC forward (15 bp) 
#R2: 10 cycles o435: UMI (10 bp) 
#R3: 20 cycles o371: RC_mBC (15 bp) 
trim_R1=3
trim_R2=0
trim_R3=5


# # # # # # # # # # # # # # # 
# NOTHING TO MODIFY BELOW THIS POINT
# # # # # # # # # # # # # # # 


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



# # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # demultiplexing with index read as a fastq: 12/10/2021
# # # # # # # # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # 
# #  trimming the end of the BC read. 
# # # # # # # # # # # # # # # # # # # #

echo "trimming"


len_BC=15

new_suffix_R1="_R1_001_e"$trim_R1".fastq"
new_suffix_R2="_R2_001_e"$trim_R2".fastq"
new_suffix_R3="_R3_001_e"$trim_R3".fastq"


seqtk trimfq -e $trim_R1 ${fastq_dir}${LOOKUP}_R1_001.fastq.gz > ${fastq_dir}${LOOKUP}${new_suffix_R1}
seqtk trimfq -e $trim_R2 ${fastq_dir}${LOOKUP}_R2_001.fastq.gz > ${fastq_dir}${LOOKUP}${new_suffix_R2}
seqtk trimfq -e $trim_R3 ${fastq_dir}${LOOKUP}_R3_001.fastq.gz > ${fastq_dir}${LOOKUP}${new_suffix_R3}

# note: command does not preserve compressed format. 


# # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # read merging/error correction of BC reads
# # # # # # # # # # # # # # # # # # # # # # # # # # # 

# "merge" the reads, which effectively corresponds to error correction w/ PEAR
# need to combine the two reads that are reverse complement of each other, here R1 and R3. 
 
echo "using pear to combine BC reads"

out_file=${dir_pear_outs}${LOOKUP}"_pear_BC_"$date_str

parallel_nodes=8
pear -j $parallel_nodes -v $len_BC -m $len_BC -n $len_BC -t $len_BC \
	-f ${fastq_dir}${LOOKUP}${new_suffix_R1} \
	-r ${fastq_dir}${LOOKUP}${new_suffix_R3} \
	-o $out_file
	
# see for pear option description: https://cme.h-its.org/exelixis/web/software/pear/doc.html



# # # # # # # # # # # # # # # # # # # # # #
# # # # processing the pear output fastq
# # # # # # # # # # # # # # # # # # # # # #
echo "processing pear out"

# pear output file must be compressed for SeqIO.
gzip $out_file".assembled.fastq"

# run python script (reformatting the fastq file as: read_id \t pear_corrected_BC)
out_file2=${dir_pear_outs}${LOOKUP}"_pear_BC_seqs_"$date_str".txt.gz"

python /net/shendure/vol1/home/lalannej/vol10_projects_JB/seq007_oBC_mBC_asso_p25_p26/nobackup/demux_fastqs/output_pear_v2_20211006.py \
	-i $out_file".assembled.fastq.gz" \
	-o $out_file2



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # generating final output file w/ gRNA and BC connected in a single text file
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

echo "compressing fastqs"

# compress the fastq's:
gzip ${fastq_dir}${LOOKUP}${new_suffix_R1}
gzip ${fastq_dir}${LOOKUP}${new_suffix_R2} 
gzip ${fastq_dir}${LOOKUP}${new_suffix_R3}

# launch python script to generate the final read/BC association file (need to have 3.6.5 for the biopython module)

echo "generating merged BC-UMI file"

#simple output file
out_BC_valid_simple=${dir_merged_file}${LOOKUP}"_BC_valid_pear_simple_"$date_str".txt.gz"
out_BC_no_pear_simple=${dir_merged_file}${LOOKUP}"_BC_no_pear_simple_"$date_str".txt.gz"

R1_name="BC"
R2_name="UMI"
R3_name="RC_BC"

python /net/shendure/vol10/projects/JBL/seq012_bulk_MPRA_K562/nobackup/parse_fastqs_w_pear_20211210.py \
	--pear_file $out_file2 \
	--out_valid $out_BC_valid_simple \
	--out_no_pear_BC $out_BC_no_pear_simple \
	--in_R1 ${fastq_dir}${LOOKUP}${new_suffix_R1}".gz" \
	--in_R2 ${fastq_dir}${LOOKUP}${new_suffix_R2}".gz" \
	--in_R3 ${fastq_dir}${LOOKUP}${new_suffix_R3}".gz" \
	--R1_name $R1_name \
	--R2_name $R2_name \
	--R3_name $R3_name \
	--simple_out_bool 1



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # pile-up of reads for downstream processing
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

echo "condensing outs"

# condensing the output to a pile-up for downstream processing
out_condensed=${dir_final_outs}${LOOKUP}"_BC_pear_UMI_condensed_"$date_str".txt"

Rscript --vanilla /net/shendure/vol10/projects/JBL/seq012_bulk_MPRA_K562/nobackup/condense_BC_file_v2_20211211.R \
	$out_BC_valid_simple \
	$out_condensed \
	"BC_pear" \
	"UMI"

gzip $out_condensed



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # further pile-up of reads for downstream processing
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# condensing the output to a pile-up for downstream processing
out_condensed2=${dir_final_outs}${LOOKUP}"_BC_pear_UMI_condensed_"$date_str".txt.gz"
out_mBC_only=${dir_final_outs}${LOOKUP}"_BC_quant_"$date_str".txt"

Rscript --vanilla /net/shendure/vol10/projects/JBL/seq013_piggyF_mBC_QTL_EB/nobackup/demux_fastqs/pileup_mBC_20211213.R \
	$out_condensed2 \
	$out_mBC_only

gzip $out_mBC_only


