import pandas as pd
from tqdm import tqdm
import csv
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
 
    # Open a new CSV file for error_checking
    csvfile = open('./bio-data-test/master.csv', 'w', newline='')
    writer = csv.writer(csvfile)
    writer.writerow(['code', 'SAb_chain', 'PDB_chain', 'PDB_auth', 'error']) 
    
    for pdb_code in tqdm(pdb_codes):
        line = [None,None,[],[],None] 
        # 1. process it
        antigen = process(pdb_code,line)
        # 2. add_SAB_info into it, for further checking
        antigen = add_SAB_info(antigen, SABdf)

        line[0] = antigen.code
        line[1] = str(antigen.chain_id_SAB)
        line[2] = str(antigen.chain_id_PDB)
        line[3] = str(antigen.auth_id_PDB)

        writer.writerow(line)
    
    csvfile.close()

"""
    To process antigen (by its pdb_code, e.g. '1a14')
    And record all errors in error.log
"""

def process(code,line):
    antigen = Antigen(code)
    if antigen.ori_fasta == None:
        line[4] = 'No valid ori_fasta.'
    else:
        antigen = process_fasta(antigen)
        
        if antigen.ori_pdb == None:
            line[4] = 'No valid ori_pdb.'
        elif len(antigen.fasta_list) == 0:
            line[4] = 'No fasta generated.'
        
        else: # process_pdb
            antigen = process_pdb(antigen)
            if len(antigen.pdb_list) == 0:
                line[4] = 'No pdb generated.'
    
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

main()