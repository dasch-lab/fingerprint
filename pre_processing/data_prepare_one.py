import os
import subprocess
import sys

import os
from tqdm import tqdm
import multiprocessing
import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=FutureWarning)

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
pdb_names = "/disk1/fingerprint/SAbDab_preparation/SAbDab_30_resolution.txt"

folder_path = "/disk1/fingerprint/provaESMFold"
precomputation = "/disk1/fingerprint/data_preparation/04b-precomputation_12A/precomputation/"

# Check if the folder already exists
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

def remove_files_raw(folder_path, pdb_id):
    # List all files in the folder
    files = os.listdir(folder_path)
    # Iterate through the files
    for file in files:
        # Check if the file starts with the PDB ID
        if file.startswith(pdb_id):
            # Construct the full file path
            file_path = os.path.join(folder_path, file)
            # Remove the file
            os.remove(file_path)


def remove_files(folder_path, extension=None):
    # Iterate over the files in the folder
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        try:
            if os.path.isfile(file_path):
                # Remove the file
                os.remove(file_path)
                print(f"Removed file: {file_path}")
        except Exception as e:
            print(f"Error deleting file {file_path}: {e}")

def main():
    # Get the root directory of the Git repository
    masif_root = subprocess.run(['git', 'rev-parse', '--show-toplevel'], capture_output=True, text=True).stdout.strip()
    masif_source = os.path.join(masif_root, 'pre_processing')

    # Add masif_source to PYTHONPATH
    #os.environ['PYTHONPATH'] = os.pathsep.join([os.environ.get('PYTHONPATH', ''), masif_source])
    with open(pdb_names, 'r') as f:
        pdb_name_list = f.read().splitlines()
    
    pdb_name_list = pdb_name_list[4209:]

    for name in tqdm(pdb_name_list):

        # Extract PDB_ID, CHAIN1, and CHAIN2 from the input argument
        #pdb_id, chain1, chain2 = sys.argv[1].split('_')
        pdb_folder = precomputation + name
        if not os.path.exists(folder_path) or len(os.listdir(folder_path)) == 0:

            pdb_id, chain1, chain2 = name.split('_')

            # Load environment if necessary

            try:
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
                #remove_files_raw("/disk1/fingerprint/data_preparation/00-raw_pdbs", pdb_id)
                #remove_files("/disk1/fingerprint/data_preparation/01-benchmark_pdbs")
                #remove_files("/disk1/fingerprint/tmp")
                remove_files_raw("/disk1/fingerprint/data_preparation/00-raw_pdbs", pdb_id)
                remove_files_raw("/disk1/fingerprint/tmp", pdb_id)
                remove_files_raw("/disk1/fingerprint/tmp", "msms")
                remove_files_raw("/disk1/fingerprint/data_preparation/01-benchmark_pdbs", pdb_id)
            except Exception:
                print(f'Issue with {name}')
                remove_files_raw("/disk1/fingerprint/data_preparation/00-raw_pdbs", pdb_id)
                remove_files_raw("/disk1/fingerprint/tmp", pdb_id)
                remove_files_raw("/disk1/fingerprint/data_preparation/01-benchmark_pdbs", pdb_id)
        else:
            print(f'Already processed')

        

if __name__ == "__main__":
    main()
    warnings.filterwarnings("default")