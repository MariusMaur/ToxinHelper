import argparse
from collections import defaultdict

def read_cluster_file(cluster_file):
    cluster_counts = defaultdict(int)
    protein_to_cluster_count = {}

    with open(cluster_file, 'r') as f:
        current_cluster = None
        for line in f:
            if line.startswith('>Cluster'):
                if current_cluster is not None:
                    # Update count for proteins in the previous cluster
                    for protein_id in protein_to_cluster_count:
                        if protein_to_cluster_count[protein_id] == current_cluster:
                            protein_to_cluster_count[protein_id] = cluster_counts[current_cluster]
                current_cluster = line.split()[1]
            else:
                # Increment count for current cluster
                cluster_counts[current_cluster] += 1

                # Extract full protein ID including '.p' and number
                protein_id = line.split('>')[1].split('...')[0]

                # Map this protein ID to the current cluster
                protein_to_cluster_count[protein_id] = current_cluster

        # Update for the last cluster
        for protein_id in protein_to_cluster_count:
            if protein_to_cluster_count[protein_id] == current_cluster:
                protein_to_cluster_count[protein_id] = cluster_counts[current_cluster]

    return protein_to_cluster_count




def read_signalp(signalp_file):
    ips_hash = {}
    with open(signalp_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                tokens = line.split()
                if len(tokens) < 7:
                    print(f"WARNING: Skipping incomplete line in signalp file: {line}")
                    continue
                try:
                    # Get both the start and end of the range
                    signalp_range = tokens[6].split('-')
                    signalp_start = int(signalp_range[0])
                    signalp_end = int(signalp_range[-1].rstrip('.'))
                except ValueError:
                    print(f"WARNING: Unable to process signalp position: {tokens[6]}")
                    continue
                ips_hash[tokens[0]] = {'signalp_start': signalp_start, 'signalp_end': signalp_end}
    return ips_hash

def extract_from_fasta(fasta_file, ips_hash):
    fasta_lengths = {}
    fasta_cys_counts = {}
    with open(fasta_file, 'r') as f:
        seq_id = None
        sequence = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    fasta_lengths[seq_id] = len(sequence)
                    if seq_id in ips_hash and 'signalp_end' in ips_hash[seq_id]:
                        # Extract part of the sequence after signalp position
                        post_signalp_seq = sequence[ips_hash[seq_id]['signalp_end']:].upper()
                        fasta_cys_counts[seq_id] = post_signalp_seq.count("C")
                    else:
                        fasta_cys_counts[seq_id] = "-"  
                seq_id = line[1:].split()[0]
                sequence = ""
            else:
                sequence += line
        if seq_id:
            fasta_lengths[seq_id] = len(sequence)
            if seq_id in ips_hash and 'signalp_end' in ips_hash[seq_id]:
                # Extract part of the sequence after signalp position for the last sequence
                post_signalp_seq = sequence[ips_hash[seq_id]['signalp_end']:].upper()
                fasta_cys_counts[seq_id] = post_signalp_seq.count("C")
            else:
                fasta_cys_counts[seq_id] = "-"  
    return fasta_lengths, fasta_cys_counts


def read_groups(groups_file):
    groups_hash = {}
    with open(groups_file, 'r') as f:
        for line in f:
            protein, group = line.strip().split()
            groups_hash[protein] = group
    return groups_hash

def main():
    parser = argparse.ArgumentParser(description='Summarizes a bunch of files for manual toxin annotation.')
    parser.add_argument('--cluster', required=True, help='Path to cluster file')
    parser.add_argument('--ips', required=True, help='Path to interproscan file')
    parser.add_argument('--signalp', required=True, help='Path to signalp file')
    parser.add_argument('--toxprot', required=True, help='Path to toxprot file')
    parser.add_argument('--fasta', required=True, help='Path to fasta file')
    parser.add_argument('--out', required=True, help='Path to output file')
    parser.add_argument('--groups_file', required=True, help='Path to groups file')

    args = parser.parse_args()

    ips_hash = read_signalp(args.signalp)
    fasta_lengths, fasta_cys_counts = extract_from_fasta(args.fasta, ips_hash)
    groups_hash = read_groups(args.groups_file)
    cluster_info = read_cluster_file(args.cluster)

    # Read ipscan file into hash
    with open(args.ips, 'r') as f:
        for line in f:
            line = line.strip()
            tokens = line.split('\t')
        
            protein_id = tokens[0]
            if protein_id not in ips_hash:
                ips_hash[protein_id] = {} 

            ips_hash[protein_id]['length'] = tokens[2]

            if tokens[3] == "PANTHER":
                description = tokens[5]
                if description in ['-', 'UNCHARACTERIZED']:
                    description = f"{tokens[4]} ({description})"
                if ":" in tokens[4]:
                    ips_hash[tokens[0]].setdefault('PANTHER_SUB', {})[tokens[4]] = description
                else:
                    ips_hash[tokens[0]].setdefault('PANTHER', {})[tokens[4]] = description
            elif len(tokens) == 11:
                ips_hash[tokens[0]].setdefault('OTHER', {}).setdefault(tokens[3], {})[tokens[4]] = tokens[5]
            elif len(tokens) == 13:
                if tokens[11] != '-' or tokens[12] != '-':
                    if tokens[11] == '-':
                        ips_hash[tokens[0]].setdefault('INTERPRO', {})[tokens[12]] = tokens[12]
                    elif tokens[12] == '-':
                        ips_hash[tokens[0]].setdefault('INTERPRO', {})[tokens[11]] = tokens[11]
                    else:
                        ips_hash[tokens[0]].setdefault('INTERPRO', {})[tokens[11]] = f"{tokens[11]}:{tokens[12]}"
            else:
                print(f"{args.ips} has the wrong number of columns")
                exit()

    # Parse toxprot to get best hit based on e-values for each ID
    best_hits = {}
    with open(args.toxprot, 'r') as f:
        for line in f:
            data = line.strip().split('\t')
            transcript_id = data[0]
            hit = data[1]
            evalue = float(data[10])  # E-value is the 11th column
            # Update the dictionary if needed
            if transcript_id not in best_hits or evalue < best_hits[transcript_id]['evalue']:
                best_hits[transcript_id] = {'hit': hit, 'evalue': evalue}

    # Produce output
    with open(args.out, 'w') as out_f:
    # Writing headers
        headers = ["Protein ID", "N Prots in CD100 cluster", "Group", "SP position", "Prot Length", "Cysteines after SP", "Panther", "InterPro", "Others", "Toxprot best hit (<1e-3)"]
        out_f.write('\t'.join(headers) + '\n')
        for protein in sorted(ips_hash.keys()):
            signalp_info = ips_hash[protein].get('signalp_start', None)
            if signalp_info is not None:
                signalp_start = signalp_info
                signalp_end = ips_hash[protein].get('signalp_end', 'no end')
                sp_column = f"{signalp_start}-{signalp_end}"
            else:
                sp_column = '-'
            seq_length = fasta_lengths.get(protein, 'unknown length')
            cys_count = fasta_cys_counts.get(protein, 'unknown cysteine count')
            panther = []
            interpro = []
            others = []
            toxprot_hit = best_hits.get(protein, {}).get('hit', 'no-toxprot')
            panther_sub = 0

            if 'PANTHER' in ips_hash[protein] or 'PANTHER_SUB' in ips_hash[protein] or 'INTERPRO' in ips_hash[protein]:
                if 'PANTHER' in ips_hash[protein]:
                    for cur_panther, current in sorted(ips_hash[protein]['PANTHER'].items()):
                        panther.append(current)
                        if current != '-' and current != 'UNCHARACTERIZED':
                            panther_sub = 1
                    if panther_sub == 0 and 'PANTHER_SUB' in ips_hash[protein]:
                        for _, current in sorted(ips_hash[protein]['PANTHER_SUB'].items()):
                            panther.append(current)
                if 'INTERPRO' in ips_hash[protein]:
                    for cur_interpro, current in sorted(ips_hash[protein]['INTERPRO'].items()):
                        interpro.append(current)
            elif 'OTHER' in ips_hash[protein]:
                for other_key, val in sorted(ips_hash[protein]['OTHER'].items()):
                    for sub_key, sub_val in sorted(val.items()):
                        others.append(f"{other_key}:{sub_key}:{sub_val}")

            if not panther:
                panther = ['no-hits']
            if not interpro:
                interpro = ['no-domains']

            group = groups_hash.get(protein, 'no-group')  # Fetch group from groups_hash
            cluster_count = cluster_info.get(protein, 'no-cluster-info')


            out_f.write(f"{protein}\t{cluster_count}\t{group}\t{sp_column}\t{seq_length}\t{cys_count}\t{'|'.join(panther)}\t{'|'.join(interpro)}\t{'|'.join(others)}\t{toxprot_hit}\n")

if __name__ == "__main__":
    main()


