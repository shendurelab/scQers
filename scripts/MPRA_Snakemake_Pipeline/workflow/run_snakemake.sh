#!/bin/bash

if conda info --envs | grep -q maggie_snakemake; then echo "maggie_snakemake environment already exists"; else conda env create -f environment.yaml; fi

. "$HOME/micromamba/etc/profile.d/mamba.sh"
. "$HOME/micromamba/etc/profile.d/conda.sh"

mamba activate maggie_snakemake

check_for_config_file() {
    argument_name="${1}"
    file_path="${2}"
    if [[ ! -f ${file_path} ]]; then
        echo "Error: filepath ${file_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

##################################################################################################################################
#############################################---Step 1: set up parameters---######################################################
##################################################################################################################################
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
            profile="${2}"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

echo "input files read"

if [ -z $profile ]; then
    echo "no profile selected"
    echo "job will be run sequentially"
    snakemake --configfile ${config_file} -c 1
else
    if [ ! -d $profile ]; then
        echo "profile directory (${profile}) cannot be found"
    else
        snakemake --configfile ${config_file} -c 1 --profile ${profile}
    fi
fi
