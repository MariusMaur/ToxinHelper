import argparse
import csv

def process_line(group_id, protein_ids, majority_protein_ids, razor_unique_counts, unique_counts, mol_weight, q_value, score, ibaq):
    protein_ids = protein_ids.split(';')
    razor_unique_counts = razor_unique_counts.split(';')
    unique_counts = unique_counts.split(';')

    if len(protein_ids) != len(razor_unique_counts) or len(protein_ids) != len(unique_counts):
        raise ValueError("Number of protein IDs and counts do not match.")

    for i, protein_id in enumerate(protein_ids):
        yield [protein_id, group_id, majority_protein_ids, razor_unique_counts[i], unique_counts[i], mol_weight, q_value, score, ibaq]

def main(input_file):
    columns = [0, 1, 3, 4, 31, 34, 35, 50]  # Corresponding to 1, 2, 4, 5, 32, 35, 36, 51 in 1-indexed format
    group_id = 1  # Starting Group ID
    custom_headers = ["Protein IDs", "MQProteinGroup", "Majority protein IDs", "Razor + unique peptides", "Unique peptides", "Mol. weight [kDa]", "Q-value", "Score", "iBAQ"]

    with open(input_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        is_header = True
        for line in reader:
            if is_header:
                # Print custom header
                print('\t'.join(custom_headers))
                is_header = False
                continue

            selected_line = [line[i] for i in columns]
            for row in process_line(str(group_id), *selected_line):
                print('\t'.join(row))
            group_id += 1  # Increment Group ID for each new line

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process MQ ProteinGroup data.')
    parser.add_argument('--MQ_proteinGroup', type=str, required=True, help='Tab-separated MQ ProteinGroup file')
    args = parser.parse_args()

    main(args.MQ_proteinGroup)
