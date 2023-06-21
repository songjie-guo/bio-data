from biopandas.pdb import PandasPdb
from config import file_path

class Protein:
    def __init__(self,code):
        self.code = code
        self.log = None
        self.ori_fasta = self.get_ori_fasta() # list
        self.ori_pdb = self.get_ori_pdb() # Dataframe
        self.chains = {}
    
    def get_ori_fasta(self):
        ori_fasta_path = f"{file_path}/ori_fasta/{self.code}.fasta"
        with open(ori_fasta_path, "r") as f:
            fastaF = f.readlines()
        if fastaF[0] == "No fasta files were found.":
            self.log = "No valid ori_fasta"
            return None
        else:
            return fastaF

    def get_ori_pdb(self):
        ori_pdb_path = f"{file_path}/ori_pdb/{self.code}.pdb"
        with open(ori_pdb_path, 'r') as f:
            first_line = f.readline()
        if first_line[0:6] != 'HEADER':
            self.log = "No valid ori_pdb"
            return None
        else:
            pdbF = PandasPdb()
            pdbF = pdbF.read_pdb(ori_pdb_path).df  #a dictionary
            return pdbF

class Chain:
    def __init__(self,name):
        self.name = name # chain name
        self.entryID = None
        self.type = '' # H/L/antigen
        self.sequence = '' # long string
        self.CDR = [] # H1,H2,H3/L1,L2,L3/[]
