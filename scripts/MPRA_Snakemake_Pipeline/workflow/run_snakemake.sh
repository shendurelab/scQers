#!/bin/bash

if conda info --envs | grep -q snakemake; then echo "snakemake environment already exists"; else conda env create -f environment.yaml; fi
conda init bash 
conda activate snakemake

snakemake --configfile "/home/maurertm/smontgom/maurertm/MPRA/MPRA_snakemake_pipeline/config/config.yaml" -c 2 --profile "/home/maurertm/smontgom/maurertm/MPRA/MPRA_snakemake_pipeline/config/slurm_scg/"
