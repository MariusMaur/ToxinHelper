#!/bin/bash

# Variables for resources; default to previous values if not specified
JOB_NAME=${JOB_NAME:-toxprot}
TIME=${TIME:-04:00:00}
NTASKS=${NTASKS:-1}
CPUS=${CPUS:-4}
MEM=${MEM:-16G}
PARTITION=${PARTITION:-normal}

#SBATCH --job-name=$JOB_NAME
#SBATCH --time=$TIME
#SBATCH --ntasks=$NTASKS
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem=$MEM
#SBATCH --partition=$PARTITION

# Load environment
module purge
module load BLAST+/2.13.0-gompi-2022a
module load SeqKit/2.3.1

# Variables
MATURE_PEPTIDES=$(realpath $1)
OUTPUT_DIR=2_toxprot_hits_and_alignments

# Make output dir
mkdir -p ${OUTPUT_DIR}

# Run blastp on mature peptides
blastp -query ${MATURE_PEPTIDES} -db toxprot/toxprot -out ${OUTPUT_DIR}/toxprot_hits.txt -evalue 1e-3 -outfmt 6

# Get alignments
OUTPUT_ALN=${OUTPUT_DIR}/alignments

# Make output dir for alignments
mkdir -p "${OUTPUT_ALN}"

# Make output dir for alignments
mkdir -p "${OUTPUT_ALN}"

# Declare an associative array to hold hits for each query
declare -A query_hits

while read -r line; do
    query_out=$(echo "$line" | cut -f1 | tr '/' '_')  # Modified line
    query=$(echo "$line" | cut -f1)
    hit=$(echo "$line" | cut -f2)

    echo "Debug: OUTPUT_ALN=$OUTPUT_ALN"  # Debug statement
    echo "Debug: query=$query"  # Debug statement

    # Check if this query has been seen before. If not, fetch its sequence.
    if [[ -z "${query_hits[$query]}" ]]; then
        seqkit grep -p "$query" "$MATURE_PEPTIDES" > "${OUTPUT_ALN}/${query_out}_with_hits.fasta"
        query_hits["$query"]=1  # Mark the query as seen
    fi

    # Add hit to the associative array
    seqkit grep -p "$hit" toxprot/uniprotkb_taxonomy_id_33208_AND_cc_tiss_2023_10_19_NO_SPACES.fasta >> "${OUTPUT_ALN}/${query_out}_with_hits.fasta"

done < "${OUTPUT_DIR}/toxprot_hits.txt"

