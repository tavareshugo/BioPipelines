#!/bin/bash
#SBATCH -A <ACCOUNT>
#SBATCH -D <WORKDIR>
#SBATCH -o logs/00-fastqdump_%a.log
#SBATCH -p <PARTITION>
#SBATCH -c 4
#SBATCH -t 02:00:00
#SBATCH -a 2-<EDIT>

source $(conda info --base)/etc/profile.d/conda.sh
conda activate sra

# fetch current SRA id
sra=$(cat sra_accessions.csv | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

# prefetch
prefetch ${sra}

# validate
vdb-validate ${sra}

# convert
fasterq-dump --outdir  reads/ ${sra}

# gzip
cd reads/
gzip ${sra}.fastq

