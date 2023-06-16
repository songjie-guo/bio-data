from biopandas.pdb import PandasPdb

# it is unique antigen thing
class Antigen:
    def __init__(self, pdb_code):
        self.code = pdb_code
        self.ori_file_path = './bio-data/ori_fasta'
        self.ori_fasta = self.get_ori_fasta() # list
        self.ori_pdb = self.get_ori_pdb() # Dataframe
        self.chain_id_SAB = [] # ?
        self.chain_info_PDB = None # string
        self.chain_id_PDB = [] #[[]]
        self.auth_id_PDB = [] # [[]]
        self.fasta_list = [] # [[]]
        self.pdb_list = [] # [Dataframe]
    
    def get_ori_fasta(self):
        ori_fasta_path = f"{self.ori_file_path}/{self.code}.fasta"
        with open(ori_fasta_path, "r") as f:
            fastaF = f.readlines()
        if fastaF[0] == "No fasta files were found.":
            return None
        else:
            return fastaF
    
    def get_ori_pdb(self):
        ori_pdb_path = f"./bio-data/ori_pdb/{self.code}.pdb"
        with open(ori_pdb_path, 'r') as f:
            first_line = f.readline()
        if first_line[0:6] != 'HEADER':
            return None
        else:
            pdbF = PandasPdb()
            pdbF = pdbF.read_pdb(ori_pdb_path).df  #a dictionary
            return pdbF