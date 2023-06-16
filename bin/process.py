from process_antigen_fasta import process_antigen_fasta
from process_antigen_pdb import process_antigen_pdb

import csv
from tqdm import tqdm

def all_pdb_codes(tsv_files):
    pdb_codes = []
    for tsv in tsv_files:
        with open(tsv, 'r') as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter='\t')
            for row in reader:
                pdb_codes.append(row['pdb'])

    pdb_codes = list(set(pdb_codes))
    return pdb_codes

def process(pdb_code):
    # download(pdb_code)
    process_antigen_fasta(pdb_code)
    process_antigen_pdb(pdb_code)
    


tsv_files = ["./bio-data/peptide.tsv", "./bio-data/protein.tsv"]
pdb_codes = sorted(all_pdb_codes(tsv_files))
for pdb_code in tqdm(pdb_codes):
    process(pdb_code)