import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Run snakemake pipeline")
    parser.add_argument('--account', required=True, help='Input slurm account')
    parser.add_argument('--fasta', required=True, help='Path to the input FASTA file')
    args = parser.parse_args()
        
	# Function to create slurm dir outputs
    if not os.path.exists("0_slurm_logs"):
        os.mkdir("0_slurm_logs")
    
	# Snakemake command
    os.system(f"""snakemake --snakefile snakefiles/Snakefile --config fasta_file={args.fasta} \
                --cluster 'sbatch --account={args.account} \
                --output=0_slurm_logs/{{params.rule_name}}_slurm_%j.log \
                --error=0_slurm_logs/{{params.rule_name}}_slurm_%j.log \
                --time={{resources.time}} --mem={{resources.mem}}MB --cpus-per-task={{resources.cpu}}' --jobs 9999""")
if __name__ == '__main__':
    main()

