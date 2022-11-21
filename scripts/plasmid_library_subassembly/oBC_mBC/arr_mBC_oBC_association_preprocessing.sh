#$ -q shendure-short.q
#$ -cwd
#$ -S /bin/bash
#$ -o /net/shendure/vol1/home/lalannej/sge_logs
#$ -e /net/shendure/vol1/home/lalannej/sge_logs
#$ -m ae
#$ -l mfree=8G,h_rt=6:00:00:00
#$ -pe serial 5
#$ -tc 4
#$ -t 1-4:1

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

LOOKUP_FILE="/net/shendure/vol1/home/lalannej/vol10_projects_JB/seq020_p25_mBC_oBC_p22_mBC_oBC_MPRA_scRep_CS2/nobackup/idx1_demux/p25_mBC_oBC_asso/lookup_p25_mBC_oBC_asso_seq020_20220319.txt"
dir_pear_outpear_outs="/net/shendure/vol1/home/lalannej/vol10_projects_JB/seq020_p25_mBC_oBC_p22_mBC_oBC_MPRA_scRep_CS2/nobackup/idx1_demux/p25_mBC_oBC_asso/pear_outs/"




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


# from seq020: 
#R1: 30 cycles o346: oBC forward (16 bp) —> requires 14 trimming from end (this is read1)
#R2: 18 cycles o433: RC_oBC (16 bp) —> requires 2 trimming from end (this is index2). 
#R3: 18 cycles o334: RC_mBC (15 bp) —> requires 3 trimming from end (this is read2)

trim_R1=14
trim_R2=2
trim_R3=3

new_suffix_R1="_R1_001_e"$trim_R1".fastq"
new_suffix_R2="_R2_001_e"$trim_R2".fastq"
new_suffix_R3="_R3_001_e"$trim_R3".fastq"


seqtk trimfq -e $trim_R1 ${LOOKUP}_R1_001.fastq.gz > ${LOOKUP}${new_suffix_R1}
seqtk trimfq -e $trim_R2 ${LOOKUP}_R2_001.fastq.gz > ${LOOKUP}${new_suffix_R2}
seqtk trimfq -e $trim_R3 ${LOOKUP}_R3_001.fastq.gz > ${LOOKUP}${new_suffix_R3}


# note: command does not preserve compressed format. 


# # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # read merging/error correction of BC reads
# # # # # # # # # # # # # # # # # # # # # # # # # # # 

# "merge" the reads, which effectively corresponds to error correction w/ PEAR
# need to combine the two reads that are reverse complement of each other, here R1 and R3. 
 

out_file=$dir_pear_outpear_outs${LOOKUP}"_pear_BC_"$date_str

len_BC=16
parallel_nodes=8
pear -j $parallel_nodes -v $len_BC -m $len_BC -n $len_BC -t $len_BC \
	-f ${LOOKUP}${new_suffix_R1} \
	-r ${LOOKUP}${new_suffix_R2} \
	-o $out_file
	
# see for pear option description: https://cme.h-its.org/exelixis/web/software/pear/doc.html



# # # # # # # # # # # # # # # # # # # # # #
# # # # processing the pear output fastq
# # # # # # # # # # # # # # # # # # # # # #

# pear output file must be compressed for SeqIO.
gzip $out_file".assembled.fastq"

# run python script (reformatting the fastq file as: read_id \t pear_corrected_BC)
out_file2=$dir_pear_outpear_outs${LOOKUP}"_pear_BC_seqs_"$date_str".txt.gz"

python /net/shendure/vol1/home/lalannej/vol10_projects_JB/seq007_oBC_mBC_asso_p25_p26/nobackup/demux_fastqs/output_pear_v2_20211006.py \
	-i $out_file".assembled.fastq.gz" \
	-o $out_file2



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # generating final output file w/ gRNA and BC connected in a single text file
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# compress the fastq's:
gzip ${LOOKUP}${new_suffix_R1}
gzip ${LOOKUP}${new_suffix_R2} 
gzip ${LOOKUP}${new_suffix_R3}

# launch python script to generate the final read/BC association file (need to have 3.6.5 for the biopython module)

#simple output file
out_BC_valid_simple=${LOOKUP}"_BC_valid_pear_simple_"$date_str".txt.gz"
out_BC_no_pear_simple=${LOOKUP}"_BC_no_pear_simple_"$date_str".txt.gz"

# even though they are not the same reads, R1 and R3 need to be the RC of each other here.
#R1_name="mBC"
#R2_name="UMI"
#R3_name="mBC_RC"

#R1: 30 cycles o346: oBC forward (16 bp) —> requires 14 trimming from end (this is read1)
#R2: 18 cycles o433: RC_oBC (16 bp) —> requires 2 trimming from end (this is index2). 
#R3: 18 cycles o334: RC_mBC (15 bp) —> requires 3 trimming from end (this is read2)


element1_fwd="oBC" # called read1 below
element1_rev="RC_oBC" # called read3 below
element2="RC_mBC" # called read2 below

python /net/shendure/vol10/projects/JBL/seq012_bulk_MPRA_K562/nobackup/parse_fastqs_w_pear_20211210.py \
	--pear_file $out_file2 \
	--out_valid $out_BC_valid_simple \
	--out_no_pear_BC $out_BC_no_pear_simple \
	--in_R1 ${LOOKUP}${new_suffix_R1}".gz" \
	--in_R2 ${LOOKUP}${new_suffix_R3}".gz" \
	--in_R3 ${LOOKUP}${new_suffix_R2}".gz" \
	--R1_name $element1_fwd \
	--R2_name $element2 \
	--R3_name $element1_rev \
	--simple_out_bool 1



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # pile-up of reads for downstream processing
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# condensing the output to a pile-up for downstream processing
out_condensed=${LOOKUP}"_pear_condensed_"$date_str".txt"

Rscript --vanilla /net/shendure/vol10/projects/JBL/seq012_bulk_MPRA_K562/nobackup/condense_BC_file_v2_20211211.R \
	$out_BC_valid_simple \
	$out_condensed \
	${element1_fwd}"_pear" \
	$element2

gzip $out_condensed

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # further pile-up of reads for downstream processing
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# condensing the output to a pile-up for downstream processing
# out_condensed2=${LOOKUP}"_pear_condensed_"$date_str".txt.gz"
# out_mBC_only=${LOOKUP}"_quant_"$date_str".txt"

#Rscript --vanilla /net/shendure/vol10/projects/JBL/seq013_piggyF_mBC_QTL_EB/nobackup/demux_fastqs/pileup_mBC_20211213.R \
#	$out_condensed2 \
#	$out_mBC_only

#gzip $out_mBC_only


