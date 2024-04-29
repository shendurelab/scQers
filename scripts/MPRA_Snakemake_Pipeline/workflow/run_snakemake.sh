#!/bin/bash

#################################---Step 1: Load Conda Snakemake Environemnt---###################################################
if conda info --envs | grep -q maggie_snakemake; then echo "maggie_snakemake environment already exists"; else conda env create -f environment.yaml; fi

. "$HOME/micromamba/etc/profile.d/mamba.sh"
. "$HOME/micromamba/etc/profile.d/conda.sh"

mamba activate maggie_snakemake

#########################################---Step 2: Define Functions---###########################################################
check_for_config_file() {
    argument_name="${1}"
    file_path="${2}"
    if [[ ! -f ${file_path} ]]; then
        echo "Error: filepath ${file_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

check_for_profile() {
    argument_name="${1}"
    directory_path="${2}"
    if [[ ! -d ${directory_path} ]]; then
        echo "Error: filepath ${directory_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

#############################################---Step 3: Set Up Parameters---######################################################
options_array=(
    config_file
    profile
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g'):

arguments=$(getopt --options a --longoptions "${longoptions}" --name 'MPRA pipeline' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --config_file )
            config_file="${2}"; check_for_config_file "${1}" "${2}"; shift 2 ;;
        --profile )
            profile="${2}"; check_for_profile "${1}" "${2}"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done
    
########################################---Step 4: Get Number of Cores---#########################################################
#This is based on the number of directories in the input_directory (as defined by the config file)
no_cores_filepath=$(awk -F'"' '/input_directory/{print $2}' ${config_file})
no_cores=$(find ${no_cores_filepath}/* -maxdepth 0 -type d | wc -l)

##########################################---Step 4: Clean yaml file---############################################################
#The script will produce an error if there are double slashes in filepaths. Thus, for all directories, they must NOT end in a /
#This code replaces all instances of /" with " in the yaml file.
sed -i 's/\/"/"/g' ${config_file}

###########################################---Step 4: Run Snakemake---############################################################
if [ -z $config_file ]; then
    echo "no config_file selected; please select config_file to run snakemake on"
    exit 1
elif [ -z $profile ]; then
    echo "no profile selected"
    snakemake --configfile ${config_file} -c ${no_cores} --use-conda
else 
    snakemake --configfile ${config_file} -c ${no_cores} --profile ${profile} --use-conda
fi
