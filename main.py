import pandas as pd
from tqdm import tqdm
import pickle
import datetime

from antigen import Antigen
from process import process_fasta, process_pdb

# you may change below paths:
tsv_files = ["./bio-data/peptide.tsv", "./bio-data/protein.tsv"]
log_file = './bio-data-test/error.log'

"""
    main part of the processing
"""

def main():
    SABdf = get_SAB_df(tsv_files)

    pdb_codes= sorted(set(SABdf['pdb']))
 
    antigens = {}
    for pdb_code in tqdm(pdb_codes):
        # 1. process it
        antigen = process(pdb_code)
        # 2. add_SAB_info into it, for further checking
        antigen = add_SAB_info(antigen,SABdf)
    
        antigens[antigen.code] = antigen
    
    with open('./bio-data-test/all_antigen.pkl', 'wb') as f:
        pickle.dump(antigens, f)

"""
    To process antigen (by its pdb_code, e.g. '1a14')
    And record all errors in error.log
"""

def process(code):
    antigen = Antigen(code)
    if antigen.ori_fasta == None:
        log(f'{code}: No valid ori_fasta.')
    else:
        antigen = process_fasta(antigen)
        
        if antigen.ori_pdb == None:
            log(f'{code}: No valid ori_pdb.')
        elif len(antigen.fasta_list) == 0:
            log(f'{antigen.code}: No fasta generated.')
        else: # process_pdb
            antigen = process_pdb(antigen)
            if len(antigen.pdb_list) == 0:
                log(f'{antigen.code}: No pdb generated.')
    
    return antigen

def get_SAB_df(tsv_files):
    df = pd.DataFrame()
    for tsv_file in tsv_files:
        tmp = pd.read_csv(tsv_file,sep='\t')
        df = pd.concat([df, tmp], axis=0)
    return df

def add_SAB_info(antigen,SABdf):
    for chain_id in SABdf[SABdf['pdb'] == antigen.code]['antigen_chain'].unique():
        antigen.chain_id_SAB.append(chain_id)
    return antigen

def log(message):
    with open(log_file, 'a') as f:
        # current_time = datetime.datetime.now().strftime('%Y/%m/%d %H:%M')
        # f.write(f'{current_time}\t{message}\n')
        f.write(f'{message}\n')

main()