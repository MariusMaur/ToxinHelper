# Toxin-Annotation-Helper <!-- omit in toc -->

### Table of contents <!-- omit in toc -->
- [Script overview](#script-overview)
		- [](#)
- [Usage](#usage)
	- [Prepare environment](#prepare-environment)
	- [How to run:](#how-to-run)
	- [Output](#output)
- [Program versions](#program-versions)



## Script overview


**This is written specifically for the venom group at UiO and is integrated with the supercomputer Saga at Sigma2.**

*This Snakemake pipeline takes a file of peptides/proteins and produces a summary table that can be used for toxin annotation.* Specifically, it:

- Runs CD-HIT to cluster sequences with 100% similarity to remove redundancy.
- Annotates signal peptides using SignalP5.
- Annotates sequences that have signal peptides with:
  - InterProScan
  - Tox-Prot (with an e-value of 1e-3), producing a multi-FASTA file for each hit.
- Executes a summary script which creates the output table, with the columns:

###

| Query ID | Signal peptide cleavage sites | Length | C's after Signal Peptide | Panther Annotation | IPR annotations | Best Tox-Prot Hit | 
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|


## Usage

### Prepare environment

No installation required, everything you need is already prepared. Before using the script, make sure the required environment is set up.

First, set up a tmux or screen environment (to run the pipeline while offline), then run:

```bash
module load Mamba/4.14.0-0
source ${EBROOTMAMBA}/bin/activate
conda activate /cluster/projects/nn9825k/admin/mamba/snakemake
```

### How to run:

```{bash}
python master.py --account ACCOUNT --fasta FASTA
```

```{bash}
usage: master.py [-h] --account ACCOUNT --fasta FASTA

Run snakemake pipeline

options:
  -h, --help         show this help message and exit
  --account ACCOUNT  Input slurm account
  --fasta FASTA      Path to the input FASTA file
```

### Output

```{bash}
.
├── 0_slurm_logs
├── 1_SignalP5
├── 2_toxprot_hits_and_alignments <-- TOX-PROT MULTI-FASTA's
├── 3_interproscan
├── CD1_SP_Your_File_summary.tsv <-- OUTPUT TABLE
├── master.py
├── README.md
├── slurm_scripts
├── snakefiles
├── tmp
└── toxprot
```

## Program versions
- SignalP v5.0
- CD-HIT v4.8.1
- SeqKit v2.3.1
- InterProScan v5.62-94.0
- parallel v20220722
- BLAST+ v2.13.0