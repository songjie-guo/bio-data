import pandas as pd
from tqdm import tqdm
from process import process
from write import write_all_you_want
import pdb
# you may change below paths:
tsv_files = ["./bio-data/peptide.tsv", "./bio-data/protein.tsv"]


def main():
    pdb_codes = get_pdb_codes(tsv_files)

    # proteins = []
    
    for pdb_code in tqdm(pdb_codes[1000:]):
        try: 
            protein = process(pdb_code)
            write_all_you_want(protein)
            # proteins.append(protein)
        except Exception as e:
            print(pdb_code, e)
    
    

def get_pdb_codes(tsv_files):
    df = pd.DataFrame()
    for tsv_file in tsv_files:
        tmp = pd.read_csv(tsv_file,sep='\t')
        df = pd.concat([df, tmp], axis=0)
    pdb_codes = sorted(set(df['pdb']))
    return pdb_codes

main()