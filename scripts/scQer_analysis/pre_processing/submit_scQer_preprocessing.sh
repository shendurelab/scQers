#!/bin/bash

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
#set -o xtrace -o nounset -o errexit

check_for_input_directory() {
    argument_name="${1}"
    directory_path="${2}"
    if [[ ${directory_path} != "none" ]] && [[ ! -d ${directory_path} ]]; then
        echo "Error: directory ${directory_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

check_for_output_directory() {
    argument_name="${1}"
    directory_path="${2}"
    if [[ ${directory_path} != "none" ]] && [[ ! -d ${directory_path} ]]; then
        mkdir -p ${directory_path}
    fi
}

##################################################################################################################################
#############################################---Step 1: set up parameters---######################################################
##################################################################################################################################
options_array=(
    input_directory
    output_directory
    barcode_type
    reference_data
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g'):

arguments=$(getopt --options a --longoptions "${longoptions}" --name 'MPRA pipeline' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --input_directory )
            input_directory="${2}"; check_for_input_directory "${1}" "${2}"; shift 2 ;;
        --output_directory )
            output_directory="${2}"; check_for_output_directory "${1}" "${2}"; shift 2 ;;
        --barcode_type )
            barcode_type="${2}"; shift 2 ;;
        --reference_data )
            reference_data="${2}"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done


if [[ $reference_data == "" ]]
then
    echo "using default mouse reference data"
    reference_data="/oak/stanford/groups/smontgom/maurertm/MPRA/reference_data/refdata-gex-GRCm39-2024-A"
fi

##################################################################################################################################
#############################################--STEP 4: GET ARRAY LENGTHS---#######################################################
##################################################################################################################################
code_directory=$(realpath $(dirname ${BASH_SOURCE[0]}))
parent_directory=$(dirname $input_directory)
input_file="${parent_directory}/directories"
find "${input_directory}/" -mindepth 1 -type d `#list all files in ${fastq_directory}`  > ${input_file}
        #add something to get a string at the end of the thign
       
array_length=$(wc -l < "${input_file}") #get the number of FASTQs
chmod 777 $input_file
echo $array_length


 mkdir -p "${output_directory}/logs"
    sbatch --output "${output_directory}/logs/%A_%a.log" `#put into log` \
        --error "${output_directory}/logs/%A_%a.log" `#put into log` \
        --array "1-${array_length}" `#initiate job array equal to the number of fastq files` \
        --wait `#indicates to the script not to move on until the sbatch operation is complete` \
        "${code_directory}/scQer_preprocessing.sh" \
        ${input_directory} \
        ${output_directory} \
        ${barcode_type} \
        ${reference_data} 
