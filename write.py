from biopandas.pdb import PandasPdb

file_path = './bio-data-test'

def write_all_you_want(protein):
    for chain in protein.chains.values():
        generate_fasta_file(chain,protein)
        generate_pdb_file(chain,protein)

def generate_fasta_file(chain, protein):
    with open(f'{file_path}/fasta/{protein.code}_{chain.type}_{chain.name}.fasta', 'w') as f:
        f.write(f'>{chain.type}:{chain.name}\n')
        f.write(chain.sequence)
        if chain.type != 'antigen':
            for i in range(3):
                f.write(f'>{chain.type}{i+1}\n')
                f.write(chain.CDR[i]+'\n')

def generate_pdb_file(chain,protein):
    chain_selected = chain.name.split("|")
    atoms = protein.ori_pdb['ATOM'] 
    atoms_selected = atoms[atoms['chain_id'].isin(chain_selected)]
    antigen_pdb = PandasPdb()
    antigen_pdb._df = {'ATOM':atoms_selected}
    antigen_pdb.to_pdb(path = f'{file_path}/pdb/{protein.code}_{chain.type}_{chain.name}.pdb')