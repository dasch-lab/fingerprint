'''
Extract the PDB chains analyzed and triangulate them
'''


#!/usr/bin/python
import numpy as np
import os
import Bio
import shutil
from Bio.PDB import * 
import sys
import importlib
#from IPython.core.debugger import set_trace

# Local includes
from source.default_config.masif_opts import masif_opts
from source.triangulation.computeMSMS import computeMSMS
from source.triangulation.fixmesh import fix_mesh
import pymesh
from source.input_output.extractPDB import extractPDB
from source.input_output.save_ply import save_ply
from source.input_output.read_ply import read_ply
from source.input_output.protonate import protonate
from source.triangulation.computeHydrophobicity import computeHydrophobicity
from source.triangulation.computeCharges import computeCharges, assignChargesToNewMesh
from source.triangulation.computeAPBS import computeAPBS
from source.triangulation.compute_normal import compute_normal
from source.triangulation.computeFlexibility import computeFlexibility
from sklearn.neighbors import KDTree

esm_folder = "/disk1/fingerprint/provaESMFold"

# Check if the correct number of command-line arguments is provided. If not, it prints usage information and exits
if len(sys.argv) <= 1: 
    print("Usage: {config} "+sys.argv[0]+" PDBID_A")
    print("A or AB are the chains to include in this surface.")
    sys.exit(1)


# The script processes the command-line arguments, extracts the PDB ID and chain IDs, and determines the source of the PDB 
# file based on additional conditions. It then protonates the PDB file using the protonate function. 
in_fields = sys.argv[1].split("_")
pdb_id = in_fields[0]
chain_ids1 = in_fields[1]

if (len(sys.argv)>2) and (sys.argv[2]=='masif_ligand'):
    pdb_filename = os.path.join(masif_opts["ligand"]["assembly_dir"],pdb_id+".pdb")
else:
    pdb_filename = masif_opts['raw_pdb_dir']+pdb_id+".pdb"
tmp_dir= masif_opts['tmp_dir']
protonated_file = tmp_dir+"/"+pdb_id+".pdb"
protonate(pdb_filename, protonated_file)
pdb_filename = protonated_file

# Extract chains of interest: the script extracts the chains of interest from the protonated PDB file and saves them to a new file.
out_filename1 = tmp_dir+"/"+pdb_id+"_"+chain_ids1
extractPDB(pdb_filename, out_filename1+".pdb", chain_ids1)

# Compute MSMS of surface w/hydrogens: triangulate the surface of the extracted chains, including hydrogen 
try:
    vertices1, faces1, normals1, names1, areas1 = computeMSMS(out_filename1+".pdb",\
        protonate=True)
except:
    print("Error")
    #set_trace()

# Compute "charged" vertices: if hydrogen bonding information is used (use_hbond is true), the script computes "charged" vertices using the computeCharges function
if masif_opts['use_hbond']:
    vertex_hbond = computeCharges(out_filename1, vertices1, names1)

# For each surface residue, assign the hydrophobicity of its amino acid. 
if masif_opts['use_hphob']:
    vertex_hphobicity = computeHydrophobicity(names1)

"""
if masif_opts['use_flexibility']:
    vertex_flexibility = computeFlexibility(out_filename1, vertices1, names1)
"""

# If protonate = false, recompute MSMS of surface, but without hydrogens (set radius of hydrogens to 0).
vertices2 = vertices1
faces2 = faces1

# Fix the mesh: The script forms a mesh from the vertices and faces, then applies mesh fixing with a specified resolution.
mesh = pymesh.form_mesh(vertices2, faces2)
regular_mesh = fix_mesh(mesh, masif_opts['mesh_res'])

# Compute the normals: The script computes normals for the vertices of the regularized mesh.
vertex_normal = compute_normal(regular_mesh.vertices, regular_mesh.faces)
# Assign charges on new vertices based on charges of old vertices (nearest
# neighbor)

# If hydrogen bonding, hydrophobicity, or APBS information is used, the script assigns these values to the new vertices of the regularized mesh.
if masif_opts['use_hbond']:
    vertex_hbond = assignChargesToNewMesh(regular_mesh.vertices, vertices1,\
        vertex_hbond, masif_opts)

if masif_opts['use_hphob']:
    vertex_hphobicity = assignChargesToNewMesh(regular_mesh.vertices, vertices1,\
        vertex_hphobicity, masif_opts)

if masif_opts['use_apbs']:
    vertex_charges = computeAPBS(regular_mesh.vertices, out_filename1+".pdb", out_filename1)


# This part then checks if the option to compute the interface (compute_iface) is enabled. If enabled, it computes the surface of the entire complex using computeMSMS, regularizes the mesh, finds vertices on the interface, and assigns values to the iface array. Finally, it saves the regularized mesh along with computed properties to a PLY file.
iface = np.zeros(len(regular_mesh.vertices))
if 'compute_iface' in masif_opts and masif_opts['compute_iface']:
    # Compute the surface of the entire complex and from that compute the interface.
    v3, f3, _, _, _ = computeMSMS(pdb_filename,\
        protonate=True)
    # Regularize the mesh
    mesh = pymesh.form_mesh(v3, f3)
    # I believe It is not necessary to regularize the full mesh. This can speed up things by a lot.
    full_regular_mesh = mesh
    # Find the vertices that are in the iface.
    v3 = full_regular_mesh.vertices
    # Find the distance between every vertex in regular_mesh.vertices and those in the full complex.
    kdt = KDTree(v3)
    d, r = kdt.query(regular_mesh.vertices)
    d = np.square(d) # Square d, because this is how it was in the pyflann version.
    assert(len(d) == len(regular_mesh.vertices))
    iface_v = np.where(d >= 2.0)[0]
    iface[iface_v] = 1.0
    # Convert to ply and save.
    save_ply(out_filename1+".ply", regular_mesh.vertices,\
                        regular_mesh.faces, normals=vertex_normal, charges=vertex_charges,\
                        normalize_charges=True, hbond=vertex_hbond, hphob=vertex_hphobicity,\
                        iface=iface)

else:
    # Convert to ply and save.
    save_ply(out_filename1+".ply", regular_mesh.vertices,\
                        regular_mesh.faces, normals=vertex_normal, charges=vertex_charges,\
                        normalize_charges=True, hbond=vertex_hbond, hphob=vertex_hphobicity)
if not os.path.exists(masif_opts['ply_chain_dir']):
    os.makedirs(masif_opts['ply_chain_dir'])
if not os.path.exists(masif_opts['pdb_chain_dir']):
    os.makedirs(masif_opts['pdb_chain_dir'])
shutil.copy(out_filename1+'.ply', masif_opts['ply_chain_dir']) 
shutil.copy(out_filename1+'.pdb', masif_opts['pdb_chain_dir']) 

# Finally, the script saves the processed results, including the mesh in PLY format and the extracted chain in PDB format, to specified directories.
