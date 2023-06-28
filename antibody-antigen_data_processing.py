import os
import re
import requests
import pandas as pd
from bs4 import BeautifulSoup, Comment
from tqdm import tqdm

"""
NOTE: 
By using this python file to generate fasta/pdb files for antibody-antigen pairs,
you ONLY need to create a folder, with four sub-folders: 
    antibody_fasta, antibody_pdb, 
    antigen_fasta, antigen_pdb,(if antigen exists)
and change the file_path below.
NO original files from SAbDab/PDB website will be downloaded or needed to download.
"""
sab = pd.read_csv("./bound-data/bound.tsv", sep='\t').sort_values(by='pdb') # summary file of SAbDab
pdb_codes = sab.iloc[:,0].unique() # all pdb_codes you process
file_path = './unbound-data'
if_antigen = False # for bound structures: T, undbound: F

"""
    Find the info of fv part we want from the result of Beautful Soup
"""
def extract_comments(comments):
    contents = []
    for comment in comments:
        next_element = comment.next_sibling
        while next_element and not isinstance(next_element, Comment):
            contents.append(next_element)
            next_element = next_element.next_sibling
    return contents

"""
    Extract info from Beautful Soup's <table_result> and <table_alignment>
"""
def process_table_result(table_element):
    table_data = {}
    for row in table_element.find_all('tr'):
        cells = row.find_all('td')
        if len(cells) == 2:
            key = cells[0].text.strip()
            value = cells[1].text.strip()
            table_data[key] = value
    return table_data

def process_table_alignment(table_element):
    for row in table_element.find_all('tr'):
        cells = [cell.text.strip() for cell in row.find_all('td')]
    chain = ''.join(cells)
    return chain

"""
    The main function to process fasta files for a pdb_code
"""
def process_fasta(pdb_code):
    url = f'https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/structureviewer/?pdb={pdb_code}'
    response = requests.get(url, stream=True)
    soup = BeautifulSoup(response.text, features="html.parser")

    comments = soup.find_all(string=lambda text: isinstance(text, Comment))
    fv_contents = extract_comments([c for c in comments if c == ' FV DETAILS '])
    sequence_contents = extract_comments([c for c in comments if c == ' SEQUENCES ' ])
    antigen_contents = extract_comments([c for c in comments if c == ' ANTIGEN DETAILS '])
    cdr_contents = extract_comments([c for c in comments if c == ' CDR SEQUENCES '])

    # get all info you need for processing
    fv_info = [process_table_result(table_element) for table_element in fv_contents[1::5]]
    H_sequence_info = [process_table_alignment(table_element) for table_element in sequence_contents[5::13]]
    L_sequence_info = [process_table_alignment(table_element) for table_element in sequence_contents[9::13]]
    antigen_info = [process_table_result(table_element) for table_element in antigen_contents[1::5]]
    cdr_info = [process_table_result(table_element) for table_element in cdr_contents[1::5]]

    # process fasta one by one from fv regions
    fv_num = len(fv_info)
    for i in range(fv_num):
        H_chain_name = fv_info[i]['Heavy chain']
        L_chain_name = fv_info[i]['Light chain']
        H_sequence = H_sequence_info[i]
        L_sequence = L_sequence_info[i]
        antigen_chain_name = antigen_info[i]['Antigen chains']
        antigen_sequence = antigen_info[i]['Antigen sequence']
        cdrs = cdr_info[i]

        # write antibody's HL fasta, regardless of the antigens
        with open(f'{file_path}/antibody_fasta/{pdb_code}_{H_chain_name}{L_chain_name}.fasta', 'w') as f:
            f.write(f'>:H\n')
            f.write(H_sequence+'\n')
            f.write(f'>:L\n')
            f.write(L_sequence+'\n')
            f.write(f'>:H1\n')
            f.write(cdrs['CDRH1']+'\n')
            f.write(f'>:H2\n')
            f.write(cdrs['CDRH2']+'\n')
            f.write(f'>:H3\n')
            f.write(cdrs['CDRH3']+'\n')
            f.write(f'>:L1\n')
            f.write(cdrs['CDRL1']+'\n')
            f.write(f'>:L2\n')
            f.write(cdrs['CDRL2']+'\n')
            f.write(f'>:L3\n')
            f.write(cdrs['CDRL3']+'\n')

        # if_antigen is False, stop writing antigen fasta
        if not if_antigen:
            return
        
        # if only one antigen pairs to that fv region
        if len(list(antigen_chain_name))==1:
            # write the antigen_fasta directly
            with open(f'{file_path}/antigen_fasta/{pdb_code}_{H_chain_name}{L_chain_name}_{antigen_chain_name}.fasta', 'w') as f:
                f.write(f'>:antigen\n')
                f.write(antigen_sequence)
        # if more than two antigens pair to that fv region
        # use another function to get fasta sequences from PDB website
        else:
            if antigen_chain_name != '': # e.g. idee's antigen_chain_name is ''.
                process_multiple_fasta(pdb_code,H_chain_name,L_chain_name,antigen_chain_name)

"""
    The function to get a list of chain names from PDB website's fasta file
"""
def get_chain_names(fastaF):
    chain_names = []
    for i in range(int(len(fastaF)/2)):
        chain_info = fastaF[i*2].split("|")[1]
        chain_info = chain_info.replace(" ", "").replace('Chains','').replace('Chain','')
        chain_list = chain_info.split(',')
        chain_list_new = []
        for chain in chain_list:
            if 'auth' in chain:
                pattern = r"(?<=auth)[A-Za-z]+"
                chain = re.findall(pattern, chain)
                chain = chain[0]
            chain_list_new.append(chain)
        chain_names.append(''.join(chain_list_new).replace(" ",""))
    return chain_names

"""
    The function to get fasta sequences from PDB website
"""
def process_multiple_fasta(pdb_code,H_chain_name,L_chain_name,antigen_chain_name):
    antigen_chain_name_list = antigen_chain_name.split(',')
    PDB_CODE = pdb_code.upper()
    url = f'https://www.rcsb.org/fasta/entry/{PDB_CODE}/display'
    response = requests.get(url, stream=True)
    soup = BeautifulSoup(response.text, features="html.parser")

    fastaF = [line.replace('&gt;', '') for line in str(soup).split('\n')]
    fasta_chain_names = get_chain_names(fastaF)

    for antigen_chain_name in antigen_chain_name_list:
        for i in range(len(fasta_chain_names)):
            if antigen_chain_name in fasta_chain_names[i]:
                antigen_sequence = fastaF[i*2+1]
        with open(f'{file_path}/antigen_fasta/{pdb_code}_{H_chain_name}{L_chain_name}_{antigen_chain_name}.fasta', 'w') as f:
            f.write(f'>:antigen\n')
            f.write(antigen_sequence)

"""
    Get indices of "ATOM" lines with a specific chain_name from pdb file
"""
def get_indices(chain_name,pdbF):
    indices = []
    for i in range(len(pdbF)):
        row = pdbF[i]
        if row.split() == []:
            continue
        if row.split()[0] == "ATOM":
            if row.split()[4] in list(chain_name):
                indices.append(i)
    return indices

"""
    The main function to process pdb files
    NOTE: We only process pdb based on all fasta files we can get
"""
def process_pdb():
    directory_path = "./bio-data/antigen_fasta"
    for file in tqdm(sorted(os.listdir(directory_path))):
        try:
            pdb_code,HL,antigen = file.replace('.fasta','').split('_')
            
            url = f'https://opig.stats.ox.ac.uk//webapps/sabdab-sabpred/sabdab/pdb/{pdb_code}/?scheme=chothia' # NOTE: We use "Chothia" numbering and CDR definition, same as SAbDab's default setting
            response = requests.get(url, stream=True)
            soup = BeautifulSoup(response.text, features="html.parser")
            pdbF = str(soup).split('\n')
    
            HL_indices = get_indices(HL,pdbF)
            antigen_indices = get_indices(antigen,pdbF)
    
    
            with open(f'{file_path}/antibody_pdb/{pdb_code}_{HL}.pdb','w') as f:
                for index in HL_indices:
                    f.write(pdbF[index]+'\n')
    
            if not if_antigen:
                return
            with open(f'{file_path}/antigen_pdb/{pdb_code}_{HL}_{antigen}.pdb','w') as f:
                for index in antigen_indices:
                    f.write(pdbF[index]+'\n')
        except:
            continue # you may change to see errors


# the main function
print("Start to process fasta.")
for pdb_code in tqdm(pdb_codes):
    try:
        process_fasta(pdb_code)
    except:
        continue # you may change to see errors

print("Start to process pdb.")
process_pdb()


