from anarci import anarci
from biopandas.pdb import PandasPdb
import re
from io import StringIO


"""
    This part is to process the original files for the "fasta" files of antigens
    First, get the indices of antigen in the fasta file
    Then, generate antigen fasta files according to the indices
"""
def process_fasta(antigen):
    antigen_indices, antigen = get_indices(antigen)
    generate_fasta_file(antigen_indices, antigen)
    return antigen

"""
    get indices of the beginning of antigens in the fasta file
    return:
    antigen_indices: a list of indices, for example, [0,1]
"""
def get_indices(antigen):
    fastaF = antigen.ori_fasta
    sequences = []  
    for i in range(0,len(fastaF),2):
        sequences.append(('result', fastaF[i+1].replace('\n','')))
    results = anarci(sequences, output=False)
    alignment = results[1]
    antigen_indices = [i for i, item in enumerate(alignment) if item is None]

    antigen = get_chain_indices(antigen_indices, antigen)

    return antigen_indices, antigen

"""
    get properties: chain_info_PDB, chain_id_PDB and auth_id_PDB
    chain_info_PDB, chain_id_PDB and auth_id_PDB are a list of items
"""
def get_chain_indices(antigen_indices, antigen):
    for i in (antigen_indices):
        row = antigen.ori_fasta[i*2] 
        sequence_info_list = row.split("|")

        chain_info = sequence_info_list[1]
        antigen.chain_info_PDB.append(chain_info)

        # get chain_list
        chain_list = chain_info.replace(" ", "")
        chain_list = re.sub(r'\[(.*?)\]','', chain_info).split(',')
        chain_list[0] = chain_list[0][-1]
        
        antigen.chain_id_PDB.append(chain_list)

        # get auth_list
        auth_list = chain_info.replace(" ", "").replace("auth","")
        auth_list = re.findall(r'\[(.*?)\]', auth_list)
        auth_list = [auth[0] if len(set(auth)) == 1 else auth for auth in auth_list]

        antigen.auth_id_PDB.append(auth_list)

    return antigen

"""
    generate antigen fasta files according to the indices
"""
def generate_fasta_file(antigen_indices, antigen):
    count = 1
    for i in (antigen_indices):
        sequence_info = antigen.ori_fasta[i*2] 
        sequence = antigen.ori_fasta[i*2+1]
        # for Antigen
        antigen.fasta_list.append([sequence_info,sequence])
        # for file
        with open(f'./bio-data-test/antigen_fasta/{antigen.code}_{count}.fasta', 'w') as f:
            f.write(sequence_info)  
            f.write(sequence)
        count+=1
"""
    This part is to process the original files for the "pdb" files of antigens
    First, check how many fasta files are there in antigen.fasta_list
    Then, for each antigen.fasta_list, generate a pdb file
"""
def process_pdb(antigen):
    atoms = antigen.ori_pdb['ATOM'] # DataFrame
    total_count = len(antigen.fasta_list)
    for i in range(total_count):
        antigen = generate_pdb_file(atoms,antigen,i)
    
    return antigen

def generate_pdb_file(atoms,antigen,i):
    # pdb.set_trace()
    count = i+1
    # get the chain's id we want
    chain_ids = atoms['chain_id'].unique()
    list1 = set(chain_ids).intersection(set(antigen.chain_id_PDB[i]))
    list2 = set(chain_ids).intersection(set(antigen.auth_id_PDB[i]))
    if len(list1) != 0:
        chain_selected = list1
    elif len(list2) != 0 :
        chain_selected = list2
    else:
        return antigen
    
    atoms_selected = atoms[atoms['chain_id'].isin(chain_selected)]

    # for Antigen
    antigen.pdb_list.append(atoms_selected)

    # for file
    antigen_pdb = PandasPdb()
    antigen_pdb._df = {'ATOM':atoms_selected}
    antigen_pdb.to_pdb(path=f'./bio-data-test/antigen_pdb/{antigen.code}_{count}.pdb')
    
    return antigen

    

