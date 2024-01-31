#!/bin/bash

# Variables for resources; default to previous values if not specified
JOB_NAME=${JOB_NAME:-peptides}
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

# Set up environment
module purge
module use /cluster/projects/nn9825k/software/modules/all/
module load SeqKit/2.3.1
module load SignalP/5.0b-Linux
module load CD-HIT/4.8.1-GCC-11.3.0

# Variables
fasta_file=$(realpath ${1})

# Make output dir
mkdir -p 1_SignalP5

# Make work environment
mkdir -p tmp && cd tmp 

# Remove dups by name and stop codon (*), create list of dups, and file of removed seqs
seqkit rmdup -n \
-D $(basename ${fasta_file%.fa*})_list_of_dups.txt \
-d $(basename ${fasta_file%.fa*})_dups.fa \
${fasta_file} | sed 's/*//g' | sed 's/\s//g' > $(basename ${fasta_file%.fa*})_nodups.fasta

## CD-HIT ##
cd-hit -i $(basename ${fasta_file%.fa*})_nodups.fasta -o CD1_$(basename ${fasta_file%.fa*})_nodups.fa -c 1 -n 5 -M 16000 -T 4 -d 0

## SignalP5 ##
# Split the fasta file, signalp only takes 5K seqs at a time
seqkit split CD1_$(basename ${fasta_file%.fa*})_nodups.fa -s 1000 -O seqkit_split_sequences
SPLIT_FA=$(realpath seqkit_split_sequences)

# Run signalP5
for FASTA in $(ls ${SPLIT_FA}/*.fa* | sort -rV);
        do
                echo "Running signalp for: " ${FASTA%.fa*}
                signalp -fasta  ${FASTA} -org euk -format short -mature
                echo "Finished signalp for: " ${FASTA%.fa*};
        done

# Combine signalp out files into one
cat *part*signalp5 | sed '/#/d' >> ../1_SignalP5/CD1_SP_$(basename ${fasta_file%.fa*}).signalp5

# Get sequences with signal peptides
seqkit grep -f <(cat *part*signalp5 | sed '/#/d' | grep "SP" | cut -f 1) CD1_$(basename ${fasta_file%.fa*})_nodups.fa > ../1_SignalP5/CD1_SP_$(basename ${fasta_file%.fa*}).fa

# Get mature sequences
cat *mature* > ../1_SignalP5/Mature_CD1_SP_$(basename ${fasta_file%.fa*}).fa

# Get unique sequences
cp CD1_$(basename ${fasta_file%.fa*})_nodups.fa  ../1_SignalP5/.