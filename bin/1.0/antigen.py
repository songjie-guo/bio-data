from biopandas.pdb import PandasPdb

# it is unique antigen thing
"""
    For an object (Class: Antigen), it has following properties:
        code: its pdb code, like '1a14'
        ori_fasta.ori_pdb: the original fasta/pdb file (you downloaded from website)
        chain_id_SAB: The antigen chain id in the tsv file from SAbDAb
        chain_info_PDB: The antigen chain's information in the fasta file from PDB
        chain_id_PDB: The antigen chain id in the fasta file from PDB (e.g. ['A'], if the chain_info_PDB is "Chain A [auth E]")
        auth_id_PDB: The antigen auth id in the fasta file from PDB (e.g. ['E'], if the chain_info_PDB is "Chain A [auth E]")
        fasta_list: a list of antigen fasta files (list)
        pdb_list: a list of antigen pdb files (string)
"""
class Antigen:
    def __init__(self, pdb_code):
        self.code = pdb_code
        self.ori_fasta = self.get_ori_fasta() # list
        self.ori_pdb = self.get_ori_pdb() # Dataframe
        self.chain_id_SAB = []
        self.chain_info_PDB = [] # [string]
        self.chain_id_PDB = [] #[[]]
        self.auth_id_PDB = [] # [[]]
        self.fasta_list = [] # [[]]
        self.pdb_list = [] # [string]
    
    """
        get the original fasta file
        return: "fastaF": a list
    """

    def get_ori_fasta(self):
        ori_fasta_path = f"./bio-data/ori_fasta/{self.code}.fasta"
        with open(ori_fasta_path, "r") as f:
            fastaF = f.readlines()
        if fastaF[0] == "No fasta files were found.":
            return None
        else:
            return fastaF
    """
        get the original pdb file
        return: "pdbF": a dictionary
    """
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