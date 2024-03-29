import argparse
import os
import logging
import uuid
from Bio import SeqIO
from tempfile import TemporaryDirectory
import subprocess
import multiprocessing

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

class UnionFind:
    """
    A class for the Union-Find data structure used for grouping elements
    into disjoint sets and merging these sets efficiently.
    """
    def __init__(self):
        self.parent = {}
        self.size = {}

    def find(self, item):
        if item not in self.parent:
            self.parent[item] = item
            self.size[item] = 1
            return item
        if self.parent[item] == item:
            return item
        self.parent[item] = self.find(self.parent[item])
        return self.parent[item]

    def union(self, group1, group2):
        root1 = self.find(group1)
        root2 = self.find(group2)
        if root1 != root2:
            if self.size[root1] < self.size[root2]:
                self.parent[root1] = root2
                self.size[root2] += self.size[root1]
            else:
                self.parent[root2] = root1
                self.size[root1] += self.size[root2]

def run_psiblast(temp_dir, sequence, db_name, e_value_thresh, output_file, iterations=3, inclusion_ethresh=1e-10):
    """ Run PSI-BLAST on a given sequence against a database. """
    unique_file_name = f"temp_sequence_{uuid.uuid4().hex}.fasta"
    temp_file_path = os.path.join(temp_dir, unique_file_name)

    with open(temp_file_path, 'w') as temp_file:
        SeqIO.write(sequence, temp_file, 'fasta')

    psi_blast_command = [
        "psiblast",
        "-query", temp_file_path,
        "-db", db_name,
        "-evalue", str(e_value_thresh),
        "-num_iterations", str(iterations),
        "-inclusion_ethresh", str(inclusion_ethresh),
        "-outfmt", "6 qseqid sseqid",
        "-out", output_file
    ]

    try:
        subprocess.run(psi_blast_command, check=True, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        logging.error(f"PSI-BLAST command failed: {e}")
        raise
    finally:
        os.unlink(temp_file_path)

    return output_file

def parallel_psiblast(record, temp_db_filename, e_value_thresh, blast_out_dir, temp_dir):
    logging.info(f"Processing {record.id}")
    result_file = os.path.join(blast_out_dir, f"{record.id}_psi_blast.out")
    return run_psiblast(temp_dir, record, temp_db_filename, e_value_thresh, result_file)

def get_full_sequences(hit_ids, db_file, query_records):
    db_sequences = {record.id: record for record in SeqIO.parse(db_file, "fasta")}
    query_sequences = {record.id: record for record in query_records}

    sequences = {}
    for hit_id in hit_ids:
        sequences[hit_id] = db_sequences.get(hit_id, query_sequences.get(hit_id))

    return sequences

def group_queries(all_hits):
    uf = UnionFind()
    for query_id, hits in all_hits.items():
        for other_query_id, other_hits in all_hits.items():
            if query_id != other_query_id and hits.intersection(other_hits):
                uf.union(query_id, other_query_id)

    grouped = {}
    for query_id in all_hits:
        root = uf.find(query_id)
        grouped.setdefault(root, set()).add(query_id)

    return grouped

def main(fasta_file, database, output_dir=None, num_processes=4):
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    tmp_dir_name = f"tmp_{uuid.uuid4().hex}"  # Using uuid for unique directory name
    with TemporaryDirectory(dir=output_dir, prefix=tmp_dir_name) as temp_dir:

        blast_out_dir = os.path.join(output_dir, "blast_out_files") if output_dir else "blast_out_files"
        group_files_dir = os.path.join(output_dir, "group_files") if output_dir else "group_files"

        os.makedirs(blast_out_dir, exist_ok=True)
        os.makedirs(group_files_dir, exist_ok=True)

        temp_db_filename = os.path.join(temp_dir, 'temp_db.fasta')
        with open(temp_db_filename, 'w+') as temp_db_file:
            with open(database, 'r') as db_handle:
                temp_db_file.write(db_handle.read())

        try:
            subprocess.run(["makeblastdb", "-in", temp_db_filename, "-dbtype", "prot"], check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"makeblastdb command failed: {e}")
            raise

        query_records = list(SeqIO.parse(fasta_file, "fasta"))
        all_hits = {}
        unique_groups = set()

        with multiprocessing.Pool(processes=num_processes) as pool:
            results = pool.starmap(parallel_psiblast, [(record, temp_db_filename, "1e-10", blast_out_dir, temp_dir) for record in query_records])

        for record, result_file in zip(query_records, results):
            hit_ids = set()
            with open(result_file) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        hit_ids.add(parts[1])

            full_sequences = get_full_sequences(hit_ids, database, query_records)
            all_hits[record.id] = full_sequences
            unique_groups.add(frozenset(hit_ids))

    groups_before_union_find = len(unique_groups)
    print(f"Number of unique groups based on initial BLAST hits: {groups_before_union_find}")

    query_groups = group_queries({k: set(v.keys()) for k, v in all_hits.items()})
    groups_after_union_find = len(query_groups)
    print(f"Number of groups after applying Union-Find: {groups_after_union_find}")

    group_count = 0
    for root, queries in query_groups.items():
        group_count += 1
        group_name = f"group{group_count}"
        written_hits = set()
        with open(os.path.join(group_files_dir, f"{group_name}.fasta"), "w") as group_file:
            for query_id in queries:
                for seq_id, seq in all_hits[query_id].items():
                    if seq_id not in written_hits:
                        group_file.write(f">{seq_id}\n{str(seq.seq)}\n")
                        written_hits.add(seq_id)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PSI-BLAST automation script")
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--database", required=True, help="Database in FASTA format")
    parser.add_argument("--outdir", help="Output directory for blast and group files", default=".")
    parser.add_argument("--parallel", type=int, default=4, help="Number of parallel processes (default: 4)")
    args = parser.parse_args()

    main(args.fasta, args.database, args.outdir)

