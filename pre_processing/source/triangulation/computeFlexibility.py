from Bio.PDB import *
from Bio import SeqIO
import numpy as np
from sklearn.neighbors import KDTree
#from source.input_output.ESMFold import parse_pdb_b_factors_mean
from Bio import SeqIO
from source.input_output.PDB import change_res_id

path = "/disk1/fingerprint/provaESMFold"
parser = PDBParser(QUIET=True)

def parse_pdb_b_factors_mean(structure, b_factors = None):
    if b_factors is None:
      b_factors = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                residue_id = residue.id[1]
                b_factor_sum = 0.0
                num_atoms = 0
                for atom in residue:
                    b_factor_sum += atom.get_bfactor()
                    num_atoms += 1
                if num_atoms > 0:
                    mean_b_factor = b_factor_sum / num_atoms
                    key = (chain.id, residue_id, residue.resname)
                    b_factors[key] = mean_b_factor
    return b_factors

def dictionary_bfactors(pdb_file, chains):
    #model = ESMModel()
    pdb_name = pdb_file.split('/')[-1].split('_')[0]
    b_factor_dict = {}
    for record in SeqIO.parse(pdb_file + '.pdb', "pdb-atom"):
        chain = record.id[-1]
        sequence = str(record.seq)
        out_file = pdb_name + '_' + chain + '_ESMFold.pdb'
        path = '/'.join(["/disk1/fingerprint/provaESMFold", out_file])
        #output = model.generate_model(chain = chain, data=sequence, pdb_write=True, model_path=path)
        #change_res_id(path)
        structure = parser.get_structure(path.split('.')[0], path)
        b_factor_dict = parse_pdb_b_factors_mean(structure, b_factor_dict)
    
    return b_factor_dict


def computeFlexibility(pdb_filename, vertices, names):
    '''
    The function takes three parameters: pdb_filename (filename of the protonated protein in PDB format), vertices (surface vertices of the protonated protein), and names (names of each vertex).
    '''
    # The PDB file is parsed using PDBParser from Biopython, and the protein structure is obtained
    struct = parser.get_structure(pdb_filename, pdb_filename + ".pdb")
    chains = pdb_filename.split('_')[-1]

    # A dictionary (residues) is created to store residues based on their chain and residue ID.
    residues = {}
    for res in struct.get_residues():
        chain_id = res.get_parent().get_id()
        if chain_id == "":
            chain_id = " "
        residues[(chain_id, res.get_id())] = res

    # All atoms are retrieved from the protein structure.
    atoms = Selection.unfold_entities(struct, "A")

    # An array (charge) is initialized to store charges for each vertex.
    plddt = np.array([0.0] * len(vertices))
    plddt_dict = dictionary_bfactors(pdb_filename, chains)
    
    # Go over every vertex  --> The function iterates over every vertex, extracts relevant information from the vertex name, 
    # and checks whether the atom is a backbone atom and whether it is already satisfied. If satisfied, the atom is ignored. 
    # Otherwise, the charge is computed using a helper function (computeChargeHelper) and stored in the charge array.
    for ix, name in enumerate(names):
        fields = name.split("_")
        chain_id = fields[0]
        if chain_id == "":
            chain_id = " "
        if fields[2] == "x":
            fields[2] = " "
        res_id = (" ", int(fields[1]), fields[2])
        aa = fields[3]
        atom_name = fields[4]
        # Compute the charge of the vertex
        #key = chain_id + '_' + str(res_id[1]) + '_' + aa + '_' + atom_name
        key = (chain_id, res_id[1], aa)
        plddt[ix] = plddt_dict[key]

    return plddt
    