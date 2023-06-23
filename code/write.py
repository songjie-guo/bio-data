from biopandas.pdb import PandasPdb
import datetime

from config import file_path

current_time = datetime.datetime.now().strftime('%Y/%m/%d %H:%M')

"""
   If you want / donnot want to Write fasta files/ Write PDB files/ Write summary file for all chains/ Write summary file for all protein (including all errors.)
   Please change the parameters: WrtFasta, WrtPDB, SumChain, SumProtein
"""
def write_all_you_want(protein_writer, chain_writer, protein, 
            WrtFasta=True, WrtPDB=True, SumChain=True, SumProtein=True):

    for chain in protein.chains.values():
        if protein.log != 'No valid ori_fasta' and WrtFasta:
            generate_fasta_file(chain,protein)
        if protein.log != 'No valid ori_pdb' and WrtPDB:
            generate_pdb_file(chain,protein)
    if SumProtein:
        protein_writer.writerow([current_time, protein.code, protein.log])
    if SumChain:
        if protein.chains != {}:
            for chain in protein.chains.values():
                chain_writer.writerow([f'{protein.code}_{chain.name}', chain.type, chain.sequence, str(chain.CDR)])

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
    
    if len(atoms_selected) == 0:
        chain_selected_first = []
        for chain_id in chain_selected:
            chain_id = chain_id[0]
            chain_selected_first.append(chain_id)
        atoms_selected = atoms[atoms['chain_id'].isin(chain_selected_first)]
    
    #else
    antigen_pdb = PandasPdb()
    antigen_pdb._df = {'ATOM':atoms_selected}
    antigen_pdb.to_pdb(path = f'{file_path}/pdb/{protein.code}_{chain.type}_{chain.name}.pdb')
