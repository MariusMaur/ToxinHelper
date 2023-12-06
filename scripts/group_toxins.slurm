#!/bin/bash

# Variables for resources; default to previous values if not specified
JOB_NAME=${JOB_NAME:-group_toxins}
TIME=${TIME:-01:00:00}
NTASKS=${NTASKS:-1}
CPUS=${CPUS:-4}
MEM=${MEM:-16}
PARTITION=${PARTITION:-normal}

#SBATCH --job-name=$JOB_NAME
#SBATCH --time=$TIME
#SBATCH --ntasks=$NTASKS
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem=$MEM
#SBATCH --partition=$PARTITION

# Load packages 
ml BLAST+/2.14.0-gompi-2022b
ml Biopython/1.81-foss-2022b

# Variables
MATURE_PEPTIDES=$(realpath $1)
OUTDIR="4_toxin_groups"

# Get the groups
python scripts/group_toxins.py --fasta ${MATURE_PEPTIDES} \
	--database ${MATURE_PEPTIDES} \
	--outdir ${OUTDIR}

# get the list for the summary
grep ">" ${OUTDIR}/group_files/* | awk -F: '{print $2 "\t" $1}' | sed 's/>//g' | \
sed 's/${OUTDIR}\/group_files\/group//g' | sed  's/\.fasta//g' > ${OUTDIR}/group_summary.tsv