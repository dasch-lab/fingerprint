import numpy as np
from pathlib import Path
from tqdm import tqdm
from Bio.PDB import *
from pathlib import Path
from Bio import SeqIO
import os
#from ImmuneBuilder import ABodyBuilder2
#from ImmuneBuilder.util import sequence_dict_from_fasta
#from ImmuneBuilder.constants import restypes

try:
    from data_preprocessing.ESMFold import ESMModel, parse_pdb_b_factors_mean
    from data_preprocessing.PDB import change_res_id
except Exception:
    from ESMFold import ESMModel, parse_pdb_b_factors_mean
    from PDB import change_res_id

path_esm = "/disk1/fingerprint/provaESMFold"
path_abb2 = "/disk1/fingerprint/provaABodyBuilder2"
parser = PDBParser(QUIET=True)

ele2num = {"C": 0, "H": 1, "O": 2, "N": 3, "S": 4, "SE": 5}
count = 5

def computeFlexibility(fname, abb2 = False, list_chian = None):
    fname = str(fname)
    model = ESMModel()
    pdb_name = fname.split('/')[-1].split('_')[0]
    b_factor_dict = {}
    sequence_list = []
    #predictor = ABodyBuilder2()

    for record in SeqIO.parse(fname, "pdb-atom"):
        chain = record.id[-1]
        print(chain)
        sequence_list.append(str(record.seq))
        print(len(str(record.seq)))
    
    
    """if abb2 and list_chian is not None:
        out_file = pdb_name.split('.')[0] + '_' + ''.join(list_chian) + '_ABB2.pdb'
        pdb_path = '/'.join([path_abb2, out_file])
        if len(sequence_list)==2:
            sequences_dict = {'H': sequence_list[0], 'L': sequence_list[1]}
            try:
                antibody = predictor.predict(sequences_dict)
                try:
                    antibody.save(pdb_path)
                except:
                    print("No template found for residue")
            except:
                print('Error in prediction')
            change_res_id(pdb_path)
            structure = parser.get_structure(pdb_path.split('.')[0], pdb_path)
            b_factor_dict = parse_pdb_b_factors_mean(structure, b_factor_dict)
            b_factor_dict_new = {}
            for key, value in b_factor_dict.items():
                if key[0] == 'H':
                    new_key = (list_chian[0], key[1], key[2])
                elif key[0] == 'L':
                    new_key = (list_chian[1], key[1], key[2])
                b_factor_dict_new[new_key] = value

            print(b_factor_dict_new)
            return(b_factor_dict_new)
    else:"""
    for element in sequence_list:
        out_file = pdb_name.split('.')[0] + '_' + chain + '_ESMFold.pdb'
        path = '/'.join([path_esm, out_file])
        output = model.generate_model(chain = chain, data=element, pdb_write=True, model_path=path)
        change_res_id(path)

        structure = parser.get_structure(path.split('.')[0], path)
        b_factor_dict = parse_pdb_b_factors_mean(structure, b_factor_dict)
        
    
    return b_factor_dict



def load_structure_np(fname, center, abb2 = False, chains = None):
    """Loads a .ply mesh to return a point cloud and connectivity."""
    # Load the data
    print(f'Processing {fname}')
    change_res_id(fname)
    structure = parser.get_structure("structure", fname)
    atoms = structure.get_atoms()
    b_factor_dict = computeFlexibility(fname, abb2, chains)

    coords = []
    types = []
    plddt = []
    #plddt = 0
    fine = True
    for atom in atoms:
        coords.append(atom.get_coord())
        try:
            types.append(ele2num[atom.element])
        except KeyError:
            count = count +1
            ele2num[atom.element] = count
            types.append(ele2num[atom.element])
            print(ele2num)
        key = (atom.get_parent().get_parent().id, atom.get_parent().get_id()[1], atom.get_parent().get_resname())
        try:
            plddt.append(b_factor_dict[key])
        except KeyError:
            fine = False
            print(f'Issues in {fname} for {key}')
            continue

    coords = np.stack(coords)
    types_array = np.zeros((len(types), len(ele2num)))
    for i, t in enumerate(types):
        types_array[i, t] = 1.0

    # Normalize the coordinates, as specified by the user:
    if center:
        coords = coords - np.mean(coords, axis=0, keepdims=True)
    
    if fine:
        if abb2:
            return {"xyz": coords, "types": types_array, "flexibility_ab": plddt}
        else:
            return {"xyz": coords, "types": types_array, "flexibility": plddt}
    else:
        return {"xyz": coords, "types": types_array, "flexibility": None}


def convert_pdbs(pdb_dir, npy_dir):
    print("Converting PDBs")
    os.makedirs(npy_dir, exist_ok=True)
    for p in tqdm(pdb_dir.glob("*.pdb")):
        protein = load_structure_np(p, center=False)
        np.save(npy_dir / (p.stem + "_atomxyz.npy"), protein["xyz"])
        np.save(npy_dir / (p.stem + "_atomtypes.npy"), protein["types"])
        if protein["flexibility"] is not None:
            np.save(npy_dir / (p.stem + "_atomflex.npy"), protein["flexibility"])

"""
if __name__=="__main__":
    convert_pdbs(Path("/disk1/fingerprint/dmasif/z_prova"), Path("/disk1/fingerprint/dmasif/z_prova_copy"))
"""
