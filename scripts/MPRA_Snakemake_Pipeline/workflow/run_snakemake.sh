#!/bin/bash
conda init bash 
conda activate snakemake

snakemake --configfile "/home/maurertm/smontgom/maurertm/MPRA/MPRA_snakemake_pipeline/config/config.yaml" -c 2 --profile "/home/maurertm/smontgom/maurertm/MPRA/MPRA_snakemake_pipeline/config/slurm_scg/"
