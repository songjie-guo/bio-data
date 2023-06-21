from anarci import anarci
import re

from protein import Protein, Chain
from config import num_schm, cdr_schm


def process(pdb_code):
    protein = Protein(pdb_code)
    # pdb.set_trace()
    if protein.ori_fasta is not None:
        fastaF = protein.ori_fasta
        chain_names = get_chain_names(fastaF)
        for i, chain_name in enumerate(chain_names):
            chain = Chain(chain_name)
            chain.entryID = i
            chain = process_chain(chain, fastaF)
            protein.chains[chain.name] = chain
    return protein

def get_chain_names(fastaF):
    chain_names = []
    for i in range(int(len(fastaF)/2)):
        chain_info = fastaF[i*2].split("|")[1]
        if 'auth' in chain_info:
            pattern = r"\[auth\s+(\w+)\]"
            chain_list = re.findall(pattern, chain_info)
        else:
            chain_list = chain_info.replace(" ", "")
            chain_list = re.sub(r'\[(.*?)\]','', chain_info).split(',')
            chain_list[0] = chain_list[0][-1]
        chain_names.append('|'.join(chain_list).replace(" ",""))
    return chain_names


def process_chain(chain, fastaF):
    chain.sequence = fastaF[chain.entryID*2+1]
    sequence = [('result', chain.sequence.replace('\n',''))]
    result = anarci(sequence, scheme = num_schm, output=False)
    alignment = result[1][0]
    if alignment is None:
        chain.type = 'antigen'
        return chain
    
    chain_type = result[1][0][0]['chain_type']
    if chain_type in ["H", "L", "K"]: 
        # numbering = [(str(x[0][0])+str(x[0][1]).replace(' ',''), x[1]) for x in result[0][0][0][0]]
        keys = [str(x[0][0])+str(x[0][1]).replace(' ','') for x in result[0][0][0][0]]
        values = [x[1] for x in result[0][0][0][0]]
        alignment = result[1][0][0]
        chain = process_HL(chain, keys, values, alignment)
    else:
        chain.type = 'antigen' 
    return chain

def get_cutting(cdr_schm):
    if cdr_schm == 'kabat':
        H_cutting = [31, 35, 50, 65, 95, 102]
        L_cutting = [24, 34, 50, 56, 89, 97]

    return H_cutting, L_cutting

def process_HL(chain, keys, values, alignment):
    H_cutting, L_cutting = get_cutting(cdr_schm)
    if alignment['chain_type'] == 'H':
        chain.type = 'H'
        H_indices = []
        for H_cut in H_cutting:
            for i, key in enumerate(keys):
                if key == str(H_cut):
                    H_indices.append(i)
                    break
            else:
                print(f"No H_cut found with {H_cut}")
        H1 = ''.join(values[H_indices[0]:H_indices[1]+1])
        H2 = ''.join(values[H_indices[2]:H_indices[3]+1])
        H3 = ''.join(values[H_indices[4]:H_indices[5]+1])
        chain.CDR = [H1,H2,H3]
    
    elif alignment['chain_type'] in ['L','K']:
        chain.type = 'L'
        L_indices = []
        for L_cut in L_cutting:
            for i, key in enumerate(keys):
                if key == str(L_cut):
                    L_indices.append(i)
                    break
            else:
                print(f"No H_cut found with {L_cut}")
        L1 = ''.join(values[L_indices[0]:L_indices[1]+1])
        L2 = ''.join(values[L_indices[2]:L_indices[3]+1])
        L3 = ''.join(values[L_indices[4]:L_indices[5]+1])
        chain.CDR = [L1,L2,L3]

    return chain
