#$ -q shendure-short.q
#$ -cwd
#$ -S /bin/bash
#$ -o /net/shendure/vol1/home/lalannej/sge_logs
#$ -e /net/shendure/vol1/home/lalannej/sge_logs
#$ -M lalannej@uw.edu
#$ -m ae
#$ -l mfree=8G,h_rt=6:00:00:00
#$ -pe serial 2
#$ -tc 3
#$ -t 1-3:1

module load modules modules-init modules-gs
module load python/3.6.5
module load numpy/1.18.4
module load pandas/1.0.3
module load scipy/1.5.1
module load statsmodels

set -e

LOOKUP_FILE="$(pwd)/array_job_lookup_table_scATAC_files_20210207.txt"
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



mkdir out_scATAC_${LOOKUP}_${date_str}/

./ATACDoubletDetector.sh \
/net/shendure/vol10/projects/Samuel/nobackup/10X/2020_ATAC/EBD21_crispri/cellranger_atac_count_combined_across_flowcells/EBD21_crispri_${LOOKUP}/outs/possorted_bam.bam \
/net/shendure/vol10/projects/Samuel/nobackup/10X/2020_ATAC/EBD21_crispri/cellranger_atac_count_combined_across_flowcells/EBD21_crispri_${LOOKUP}/outs/singlecell.csv \
mm10_autosomes.txt mm10_blacklist_microsat_segdups_simplereps_rmsk_20210124.bed ./out_scATAC_${LOOKUP}_${date_str}/ .

