import numpy as np
from pathlib import Path
from tqdm import tqdm
from Bio.PDB import *
from pathlib import Path
from Bio import SeqIO
from ESMFold import ESMModel, dictionary_pdb_bfactors
from PDB import change_res_id
import os

path_esm = "/disk1/fingerprint/provaESMFold"
parser = PDBParser(QUIET=True)

ele2num = {"C": 0, "H": 1, "O": 2, "N": 3, "S": 4, "SE": 5}

def computeFlexibility(fname):
    fname = str(fname)
    model = ESMModel()
    pdb_name = fname.split('/')[-1].split('_')[0]
    b_factor_dict = {}
    for record in SeqIO.parse(fname, "pdb-atom"):
        chain = record.id[-1]
        print(chain)
        sequence = str(record.seq)
        out_file = pdb_name.split('.')[0] + '_' + chain + '_ESMFold.pdb'
        path = '/'.join([path_esm, out_file])
        output = model.generate_model(chain = chain, data=sequence, pdb_write=True, model_path=path)
        change_res_id(path)
        structure = parser.get_structure(path.split('.')[0], path)
        b_factor_dict = dictionary_pdb_bfactors(structure, b_factor_dict)
    
    return b_factor_dict



def load_structure_np(fname, center):
    """Loads a .ply mesh to return a point cloud and connectivity."""
    # Load the data
    change_res_id(fname)
    structure = parser.get_structure("structure", fname)
    atoms = structure.get_atoms()
    b_factor_dict = computeFlexibility(fname)

    coords = []
    types = []
    plddt = []
    for atom in atoms:
        coords.append(atom.get_coord())
        types.append(ele2num[atom.element])
        key = atom.get_parent().get_parent().id + '_' + str(atom.get_parent().get_id()[1]) + '_' + atom.get_parent().get_resname() + '_' + atom.get_name()
        plddt.append(b_factor_dict[key])

    coords = np.stack(coords)
    types_array = np.zeros((len(types), len(ele2num)))
    for i, t in enumerate(types):
        types_array[i, t] = 1.0

    # Normalize the coordinates, as specified by the user:
    if center:
        coords = coords - np.mean(coords, axis=0, keepdims=True)

    return {"xyz": coords, "types": types_array, "flexibility": plddt}


def convert_pdbs(pdb_dir, npy_dir):
    print("Converting PDBs")
    os.makedirs(npy_dir, exist_ok=True)
    for p in tqdm(pdb_dir.glob("*.pdb")):
        protein = load_structure_np(p, center=False)
        np.save(npy_dir / (p.stem + "_atomxyz.npy"), protein["xyz"])
        np.save(npy_dir / (p.stem + "_atomtypes.npy"), protein["types"])
        np.save(npy_dir / (p.stem + "_atomflex.npy"), protein["flexibility"])

if __name__=="__main__":
    convert_pdbs(Path("/disk1/fingerprint/dMaSIF/z_prova"), Path("/disk1/fingerprint/dMaSIF/z_prova_copy"))