from anarci import anarci
import glob
from logger import log

def process_antigen_fasta(pdb_code):
    fastaF = get_ori_fasta(pdb_code)
    if f"./bio-data/antigen_fasta/{pdb_code}_1.fasta" in glob.glob('./bio-data/antigen_fasta/*'):
        return 
    if fastaF is None:
        log(f"{pdb_code}: No valid ori_fasta.")
        return

    # if fastaF is valid
    antigen_indices = get_antigen_indices(fastaF)
    generate_antigen_fasta(antigen_indices, fastaF, pdb_code)
    

def get_ori_fasta(pdb_code):
    ori_fasta_path = f"./bio-data/ori_fasta/{pdb_code}.fasta"
    with open(ori_fasta_path, "r") as f:
        fastaF = f.readlines()
    if fastaF[0] == "No fasta files were found.":
        return None
    else:
        return fastaF

def get_antigen_indices(fastaF):
    sequences = []  
    for i in range(0,len(fastaF),2):
        sequences.append(('result', fastaF[i+1].replace('\n','')))
    results = anarci(sequences, output=False)
    alignment = results[1]
    antigen_indices = [i for i, item in enumerate(alignment) if item is None]
    return antigen_indices

def generate_antigen_fasta(antigen_indices, fastaF, pdb_code):
    count = 1
    for i in (antigen_indices):
        sequence_info = fastaF[i*2] 
        sequence = fastaF[i*2+1]
        with open(f'./bio-data/antigen_fasta/{pdb_code}_{count}.fasta', 'w') as f:
            f.write(sequence_info)  
            f.write(sequence)
        count+=1