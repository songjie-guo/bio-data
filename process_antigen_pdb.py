import glob
import re
from logger import log
from biopandas.pdb import PandasPdb


def process_antigen_pdb(pdb_code):
    if f"./bio-data/antigen_pdb/{pdb_code}_1.fasta" in glob.glob('./bio-data/antigen_fasta/*'):
        return 
    
    pdbF = get_ori_pdb(pdb_code)
    if pdbF is None:
        log(f"{pdb_code}: No valid ori_pdb.")
        return
    
    # if there is at least an antigen-fasta file
    # also pdbF exist
    count = 1
    antigen_fasta_path = f'./bio-data/antigen_fasta/{pdb_code}_{count}.fasta'
    while (antigen_fasta_path) in glob.glob('./bio-data/antigen_fasta/*'):
        chain_list,auth_list = get_chain_indices(antigen_fasta_path)
        if f"./bio-data/antigen_pdb/{pdb_code}_{count}.pdb" not in glob.glob('./bio-data/antigen_pdb/*'):
            generate_antigen_pdb(chain_list,auth_list,pdbF,pdb_code,count)
        count +=1
        antigen_fasta_path = f'./bio-data/antigen_fasta/{pdb_code}_{count}.fasta'
 

        

def get_ori_pdb(pdb_code):
    ori_pdb_path = f"./bio-data/ori_pdb/{pdb_code}.pdb"
    # pdb.set_trace()
    with open(ori_pdb_path, 'r') as f:
        first_line = f.readline()
    if first_line[0:6] != 'HEADER':
        return None
    else:
        pdbF = PandasPdb()
        pdbF = pdbF.read_pdb(ori_pdb_path).df  #a dictionary
        return pdbF

def get_chain_indices(antigen_fasta_path):
    with open(antigen_fasta_path, "r") as f:
        fastaF = f.readlines()
    sequence_info = fastaF[0].split("|")
    row = sequence_info[1]

    # get chain_list
    chain_list = row.replace(" ", "")
    chain_list = re.sub(r'\[(.*?)\]','', row).split(',')
    chain_list[0] = chain_list[0][-1]

    # get auth_list
    auth_list = row.replace(" ", "").replace("auth","")
    auth_list = re.findall(r'\[(.*?)\]', auth_list)
    auth_list = [auth[0] if len(set(auth)) == 1 else auth for auth in auth_list]

    return chain_list,auth_list

def generate_antigen_pdb(chain_list,auth_list,pdbF,pdb_code,count):
    atoms = pdbF['ATOM'] # a pd.DataFrame
    # get the chain's id we want
    chain_ids = atoms['chain_id'].unique()
    list1 = set(chain_ids).intersection(set(chain_list))
    list2 = set(chain_ids).intersection(set(auth_list))
    
    if len(list1) != 0:
        chain_selected = list1
    elif len(list2) != 0 :
        chain_selected = list2
    else:
        log(f'{pdb_code}: Finding no chains.')
        return
    atoms_selected = atoms[atoms['chain_id'].isin(chain_selected)]
    
    antigen_pdb = PandasPdb()
    antigen_pdb._df = {'ATOM':atoms_selected}
    # pdb.set_trace()
    try:
        antigen_pdb.to_pdb(path=f'./bio-data/antigen_pdb/{pdb_code}_{count}.pdb')
    except Exception as e:
        log(f'{pdb_code}: {e}.')