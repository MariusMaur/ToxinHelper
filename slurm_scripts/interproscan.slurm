#!/bin/bash

# Variables for resources; default to previous values if not specified
JOB_NAME=${JOB_NAME:-interproscan}
TIME=${TIME:-08:00:00}
NTASKS=${NTASKS:-1}
CPUS=${CPUS:-32}
MEM=${MEM:-64}
PARTITION=${PARTITION:-bigmem}

#SBATCH --job-name=$JOB_NAME
#SBATCH --time=$TIME
#SBATCH --ntasks=$NTASKS
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem=$MEM
#SBATCH --partition=$PARTITION

# Load environment
module purge
module use /cluster/projects/nn9825k/software/modules/all/
module load SeqKit/2.3.1
module load InterProScan/5.62-94.0-foss-2022a
module load parallel/20220722-GCCcore-11.3.0

# Set variables
fasta_file=$(realpath $1)

# Split the file for more speed
cd tmp
seqkit split $fasta_file -s 100 -O interpro_seqkit_split_sequences
SPLIT_FA=$(realpath interpro_seqkit_split_sequences)


# Func call interpro
call_interpro (){
        echo "Running interpro for: " ${1} &> $(basename ${1}).logfile
        interproscan.sh -i ${1} --cpu 16 &>> $(basename ${1}).logfile
        echo "Finished interpro for: " ${1} &>> $(basename ${1}).logfile
}

export -f call_interpro

# Call interpro function
cd ..
mkdir -p 3_interproscan && cd 3_interproscan
parallel -j 2 call_interpro ::: ${SPLIT_FA}/*.fa*
# Combine tsv
cat $(basename ${fasta_file%.fa*})*.tsv > $(basename ${fasta_file%.fa*})_IPR.tsv
rm -r *part_00* temp