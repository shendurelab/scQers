#!/usr/bin/env bash

# conda activate snakemake
# snakemake --configfile config/pcqtl.yaml --cores 10 --use-conda --profile scg --resources -j 100 --cluster-cancel scancel

# conda activate snakemake_7.32
snakemake --config config/config.yaml --profile config/slurm_scg -nfp