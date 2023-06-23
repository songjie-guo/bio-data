import pandas as pd
from tqdm import tqdm
from process import process
from write import write
import pdb

# you may change below paths:
file_path = './bio-data-test-ah'
tsv_files = ["./bio-data/peptide.tsv", "./bio-data/protein.tsv"]


def main():
    pdb_codes,paringDF = get_pdb_codes(tsv_files)

    # proteins = []
    # ['7ufn']
    for pdb_code in tqdm(pdb_codes):
        try: 
            protein = process(pdb_code)
            write(protein,paringDF)
            # proteins.append(protein)
        except Exception as e:
            print(pdb_code, e)
            



def get_pdb_codes(tsv_files):
    df = pd.DataFrame()
    for tsv_file in tsv_files:
        tmp = pd.read_csv(tsv_file,sep='\t')
        df = pd.concat([df, tmp], axis=0)
    pdb_codes = sorted(set(df['pdb']))
    paringDF = df[['pdb','Hchain','Lchain','antigen_chain']]
    return pdb_codes,paringDF

##main

main()
