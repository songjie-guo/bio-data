from biopandas.pdb import PandasPdb
import pdb
file_path = './bio-data-test-ah'

def write(protein,paringDF):
    all_chains = get_chains(protein)
    write_fasta(protein,all_chains)
    write_pdb(protein,paringDF)

def get_chains(protein):
    H_chains = []
    L_chains = []
    antigen_chains = []
    for chain in protein.chains.values():
        if chain.type == "H":
            H_chains.append(chain)
        elif chain.type == "L":
            L_chains.append(chain)
        elif chain.type == "antigen":
            antigen_chains.append(chain)

    return [H_chains,L_chains,antigen_chains]

def write_fasta(protein,all_chains):
    H_chains,L_chains,antigen_chains = all_chains
    # pdb.set_trace()
    if len(H_chains)==1 and len(L_chains)==1 and len(antigen_chains)==1:
        H_chain = H_chains[0]
        L_chain = L_chains[0]
        antigen_chain = antigen_chains[0]
        # write the HL_fasta
        with open(f'{file_path}/antibody_fasta/{protein.code}_{H_chain.name}{L_chain.name}.fasta', 'w') as f:
            f.write(f'>{H_chain.type}\n')
            f.write(H_chain.sequence)
            f.write(f'>{L_chain.type}\n')
            f.write(L_chain.sequence)
            for i in range(3):
                f.write(f'>{H_chain.type}{i+1}\n')
                f.write(H_chain.CDR[i]+'\n')
            for i in range(3):
                f.write(f'>{L_chain.type}{i+1}\n')
                f.write(L_chain.CDR[i]+'\n')
        HL_name = H_chain.name+'-'+L_chain.name

        # write the antigen_fasta
        with open(f'{file_path}/antigen_fasta/{protein.code}_{HL_name}_{antigen_chain.name}.fasta', 'w') as f:
            f.write(f'>{antigen_chain.type}\n')
            f.write(antigen_chain.sequence)

def write_pdb(protein,paringDF):
    paringDF = paringDF[paringDF["pdb"] == protein.code].iloc[:, 1:]
    paring_lists = paringDF.values.tolist()
    for a_paring in paring_lists:
        H_name,L_name,antigen_name  = a_paring
        antigen_names = antigen_name.replace(" ",'').split('|')
        atoms = protein.ori_pdb['ATOM']

        # write the HL_pdb
        HL_selected = atoms[atoms['chain_id'].isin([H_name,L_name])]
        
        HL_pdb = PandasPdb()
        HL_pdb._df = {'ATOM':HL_selected}
        HL_pdb.to_pdb(path = f'{file_path}/antibody_pdb/{protein.code}_{H_name}-{L_name}.pdb')

        # write the antigen_pdb
        antigen_selected = atoms[atoms['chain_id'].isin(antigen_names)]

        antigen_pdb = PandasPdb()
        antigen_pdb._df = {'ATOM':antigen_selected}
        antigen_name = '_'.join(antigen_names)
        antigen_pdb.to_pdb(path = f'{file_path}/antigen_pdb/{protein.code}_{H_name}-{L_name}_{antigen_name}.pdb')
