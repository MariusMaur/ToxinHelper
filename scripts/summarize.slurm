#!/bin/bash

# Variables for resources; default to previous values if not specified
JOB_NAME=${JOB_NAME:-summarize}
TIME=${TIME:-01:00:00}
NTASKS=${NTASKS:-1}
CPUS=${CPUS:-2}
MEM=${MEM:-8}
PARTITION=${PARTITION:-normal}

#SBATCH --job-name=$JOB_NAME
#SBATCH --time=$TIME
#SBATCH --ntasks=$NTASKS
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem=$MEM
#SBATCH --partition=$PARTITION



# Variables
SIGNALP=$(realpath $1)
TOXPROT=$(realpath $2)
IPS=$(realpath $3)
FASTA=$(realpath $4)
TOXIN_GROUPS=$(realpath $5)
CLUSTER=$(realpath $6)
OUTPUT=$(basename ${FASTA%.fa*})_summary.tsv


python scripts/summarize_stats.py \
	--cluster ${CLUSTER} \
	--signalp ${SIGNALP} --ips ${IPS} \
	--toxprot ${TOXPROT} --fasta ${FASTA} \
	--groups ${TOXIN_GROUPS} --out ${OUTPUT}