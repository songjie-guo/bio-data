import pandas as pd
import csv
from tqdm import tqdm

from config import file_path, tsv_files
from process import process
from write import write_all_you_want

"""
    The main function
"""
def main():
    pdb_codes = get_pdb_codes(tsv_files)
    
    for pdb_code in tqdm(pdb_codes):
        try: 
            protein = process(pdb_code) # you process the data, get an object(Protein), you may store it for future use
            write_all_you_want(protein_writer,chain_writer,protein) # you write all fasat/pdb files, and summary files
        except Exception as e:
            print(pdb_code, e)
    
"""
    To get all pdb_codes from your search on SAbDab
"""
def get_pdb_codes(tsv_files):
    df = pd.DataFrame()
    for tsv_file in tsv_files:
        tmp = pd.read_csv(tsv_file,sep='\t')
        df = pd.concat([df, tmp], axis=0)
    pdb_codes = sorted(set(df['pdb']))
    return pdb_codes

"""
    The main part, including the two csv writer for the two summary files.
"""
summary_protein = open(f'{file_path}/summary_protein.csv', 'a', newline='')
protein_writer = csv.writer(summary_protein)
protein_writer.writerow(['time','pdb_code','error'])

summary_chain = open(f'{file_path}/summary_chain.csv', 'a', newline='')
chain_writer = csv.writer(summary_chain)
chain_writer.writerow(['pdb_chain','type','sequence','CDR'])

main()
