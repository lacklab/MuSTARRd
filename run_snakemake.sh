#!/bin/bash
#SBATCH -c 64
#SBATCH --mem 720GB
#SBATCH -p long,big-mem,normal,express
source /home/etekoglu/.bashrc
conda activate /home/etekoglu/miniconda3/envs/snakemake

snakemake --profile profile/