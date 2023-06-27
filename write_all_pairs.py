import pandas as pd
import requests
from bs4 import BeautifulSoup, Comment
import re
from tqdm import tqdm
import os
import pdb

sab = pd.read_csv("./bio-data/all.tsv", sep='\t').sort_values(by='pdb')
# sab = pd.read_csv("./bu.txt", sep=' ')
pdb_codes = sab.iloc[:,0].unique()
file_path = './bio-data'

def extract_comments(comments):
    contents = []
    for comment in comments:
        next_element = comment.next_sibling
        while next_element and not isinstance(next_element, Comment):
            contents.append(next_element)
            next_element = next_element.next_sibling
    return contents

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

def process_fasta(pdb_code):
    url = f'https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/structureviewer/?pdb={pdb_code}'
    response = requests.get(url, stream=True)
    soup = BeautifulSoup(response.text, features="html.parser")

    comments = soup.find_all(string=lambda text: isinstance(text, Comment))
    fv_contents = extract_comments([c for c in comments if c == ' FV DETAILS '])
    sequence_contents = extract_comments([c for c in comments if c == ' SEQUENCES ' ])
    antigen_contents = extract_comments([c for c in comments if c == ' ANTIGEN DETAILS '])
    cdr_contents = extract_comments([c for c in comments if c == ' CDR SEQUENCES '])

    fv_info = [process_table_result(table_element) for table_element in fv_contents[1::5]]
    H_sequence_info = [process_table_alignment(table_element) for table_element in sequence_contents[5::13]]
    L_sequence_info = [process_table_alignment(table_element) for table_element in sequence_contents[9::13]]
    antigen_info = [process_table_result(table_element) for table_element in antigen_contents[1::5]]
    cdr_info = [process_table_result(table_element) for table_element in cdr_contents[1::5]]
    
    fv_num = len(fv_info)
    for i in range(fv_num):
        H_chain_name = fv_info[i]['Heavy chain']
        L_chain_name = fv_info[i]['Light chain']
        H_sequence = H_sequence_info[i]
        L_sequence = L_sequence_info[i]
        antigen_chain_name = antigen_info[i]['Antigen chains']
        antigen_sequence = antigen_info[i]['Antigen sequence']
        cdrs = cdr_info[i]
        
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
        
        if len(list(antigen_chain_name))==1:
            # write the antigen_fasta
            with open(f'{file_path}/antigen_fasta/{pdb_code}_{H_chain_name}{L_chain_name}_{antigen_chain_name}.fasta', 'w') as f:
                f.write(f'>:antigen\n')
                f.write(antigen_sequence)
        else:
            if antigen_chain_name != '':
                process_multiple_fasta(pdb_code,H_chain_name,L_chain_name,antigen_chain_name)
            else:
                print(pdb_code)

def get_chain_names(fastaF):
    chain_names = []
    # pdb.set_trace()
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

def process_pdb():
    directory_path = "./bio-data/antigen_fasta"
    for file in tqdm(sorted(os.listdir(directory_path))[:10]):
        pdb_code,HL,antigen = file.replace('.fasta','').split('_')
        
        url = f'https://opig.stats.ox.ac.uk//webapps/sabdab-sabpred/sabdab/pdb/{pdb_code}/?scheme=chothia'
        response = requests.get(url, stream=True)
        soup = BeautifulSoup(response.text, features="html.parser")
        pdbF = str(soup).split('\n')

        HL_indices = get_indices(HL,pdbF)
        antigen_indices = get_indices(antigen,pdbF)


        with open(f'{file_path}/antibody_pdb/{pdb_code}_{HL}.pdb','w') as f:
            for index in HL_indices:
                f.write(pdbF[index]+'\n')


        with open(f'{file_path}/antigen_pdb/{pdb_code}_{HL}_{antigen}.pdb','w') as f:
            for index in antigen_indices:
                f.write(pdbF[index]+'\n')


# print("start to process fasta. 🤣")
# for pdb_code in tqdm(pdb_codes):
#     try:
#         process_fasta(pdb_code)
#     except:
#         continue
print("start to process pdb. 🤣")

process_pdb()

