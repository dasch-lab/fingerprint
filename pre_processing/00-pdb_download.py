#!/usr/bin/python
from Bio.PDB import *
import sys
import os

# Local includes
from source.default_config.masif_opts import masif_opts
from source.input_output.protonate import protonate
import shutil

if len(sys.argv) <= 1:
    print("Usage: " + sys.argv[0] + " PDBID_A_B")
    print("A or B are the chains to include in this pdb.")
    sys.exit(1)

if not os.path.exists(masif_opts['raw_pdb_dir']):
    os.makedirs(masif_opts['raw_pdb_dir'])

if not os.path.exists(masif_opts['tmp_dir']):
    os.mkdir(masif_opts['tmp_dir'])

in_fields = sys.argv[1].split('_')
pdb_id = in_fields[0]

# Download pdb
# pdbl = PDBList(server='http://files.wwpdb.org')


pdb_filename = os.path.join(masif_opts['pdb_folder'], pdb_id + ".pdb")
if os.path.exists(pdb_filename):
    shutil.copy(pdb_filename, masif_opts['raw_pdb_dir'])

pdb_filename = os.path.join(masif_opts['raw_pdb_dir'], pdb_id + ".pdb")
if not os.path.exists(pdb_filename):
    pdbl = PDBList()
    pdb_filename = pdbl.retrieve_pdb_file(pdb_id, pdir=masif_opts['tmp_dir'], file_format='pdb')

# Protonate with reduce, if hydrogens included.
# - Always protonate as this is useful for charges. If necessary ignore hydrogens later.
protonated_file = os.path.join(masif_opts['raw_pdb_dir'], pdb_id + ".pdb")
protonate(pdb_filename, protonated_file)
pdb_filename = protonated_file


"""
'''
00-pdb_download.py: Download pdb file from Protein DataBank
'''

#!/usr/bin/python
import Bio
from Bio.PDB import * 
import sys
import importlib
import os

#from default_config.masif_opts import masif_opts
from source.default_config.masif_opts import masif_opts
# Local includes
#from input_output.protonate import protonate
from source.input_output.protonate import protonate
from source.input_output.PDB import change_res_id

if len(sys.argv) <= 1: 
    print("Usage: "+sys.argv[0]+" PDBID_A_B")
    print("A or B are the chains to include in this pdb.")
    sys.exit(1)

if not os.path.exists(masif_opts['raw_pdb_dir']):
    os.makedirs(masif_opts['raw_pdb_dir'])

if not os.path.exists(masif_opts['tmp_dir']):
    os.mkdir(masif_opts['tmp_dir'])

#Extraction of PDB ID from the command line argument (e.g. PDBID_A_B)
in_fields = sys.argv[1].split('_')
pdb_id = in_fields[0]

# Download pdb 
pdbl = PDBList()
pdb_filename = pdbl.retrieve_pdb_file(pdb_id, pdir=masif_opts['tmp_dir'],file_format='pdb')


##### Protonate with reduce, if hydrogens included.
# - Always protonate as this is useful for charges. If necessary ignore hydrogens later.
protonated_file = masif_opts['raw_pdb_dir']+pdb_id+".pdb"
protonate(pdb_filename, protonated_file)
pdb_filename = protonated_file

out = change_res_id(pdb_filename)
if out:
    print(f'Changed res_id of {pdb_filename}')
else:
    print(f'Failed to change res_id of {pdb_filename}')


"""