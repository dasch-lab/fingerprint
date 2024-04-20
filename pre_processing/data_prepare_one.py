import os
import subprocess
import sys

import os
from tqdm import tqdm
# Path to the directory you want to add to PATH


def export_lib(new_path, name):
    current_path = os.environ.get("PATH", "")
    if current_path:
        new_path = current_path + ":" + new_path
    os.environ["PATH"] = new_path
    print('Added path ' + name)

os.environ["APBS_BIN"] = "/disk1/fingerprint/librerie/APBS-3.4.1.Linux/bin/apbs"
os.environ["MULTIVALUE_BIN"] = "/disk1/fingerprint/librerie/APBS-3.4.1.Linux/share/apbs/tools/bin/multivalue"
#os.environ["PDB2PQR_BIN"] = "/disk1/fingerprint/masif/APBS-3.4.1.Linux/share/apbs/tools/conversion/param/pdb2pqr"
os.environ["PDB2PQR_BIN"] = "/home/s.joubbi/miniconda3/envs/surfaceid2/bin/pdb2pqr"
export_lib("/disk1/fingerprint/librerie/reduce/build/reduce/reduce_src/", 'reduce')
os.environ["REDUCE_HET_DICT"] = "/disk1/fingerprint/librerie/reduce/reduce_wwPDB_het_dict.txt"
os.environ["MSMS_BIN"] = "/disk1/fingerprint/librerie/msms/msms.x86_64Linux2.2.6.1"
os.environ["PDB2XYZRN"] = "/disk1/fingerprint/librerie/msms/pdb_to_xyzrn"
pdb_path = "/disk1/fingerprint/SAbDab_preparation/all_structures"
pdb_names = "/disk1/fingerprint/SAbDab_preparation/sabdab_summary_output_filtered.txt"


def main():
    # Get the root directory of the Git repository
    masif_root = subprocess.run(['git', 'rev-parse', '--show-toplevel'], capture_output=True, text=True).stdout.strip()
    masif_source = os.path.join(masif_root, 'pre_processing')

    # Add masif_source to PYTHONPATH
    #os.environ['PYTHONPATH'] = os.pathsep.join([os.environ.get('PYTHONPATH', ''), masif_source])
    with open(pdb_names, 'r') as f:
        pdb_name_list = f.read().splitlines()
    

    for name in tqdm(pdb_name_list):

        # Extract PDB_ID, CHAIN1, and CHAIN2 from the input argument
        #pdb_id, chain1, chain2 = sys.argv[1].split('_')
        pdb_id, chain1, chain2 = name.split('_')

        # Load environment if necessary

        # Execute the Python scripts
        subprocess.run(['python', os.path.join(masif_source, '00-pdb_download.py'), name])
        print('Done 00-pdb_download.py')
        subprocess.run(['python', os.path.join(masif_source, '01-pdb_extract_and_triangulate.py'), f'{pdb_id}_{chain1}'])
        print(f'Done 01-pdb_extract_and_triangulate.py for {pdb_id}_{chain1}')
        subprocess.run(['python', os.path.join(masif_source, '01-pdb_extract_and_triangulate.py'), f'{pdb_id}_{chain2}'])
        print(f'Don 01-pdb_extract_and_triangulate.py for {pdb_id}_{chain1}')
        #subprocess.run(['python', os.path.join(masif_source, '04-masif_precompute.py'), 'masif_site', name])
        #print(f'Done 04-masif_precompute.py for masif_site')
        subprocess.run(['python', os.path.join(masif_source, '04-masif_precompute.py'), 'masif_ppi_search', name])
        print(f'Done 04-masif_precompute.py for masif_ppi_search')

if __name__ == "__main__":
    main()