from Bio.PDB import *
from Bio import SeqIO
import numpy as np
from sklearn.neighbors import KDTree

"""
computeCharges.py: Wrapper function to compute hydrogen bond potential (free electrons/protons) in the surface
Pablo Gainza - LPDI STI EPFL 2019
This file is part of MaSIF.
Released under an Apache License 2.0
"""

from source.default_config.chemistry import (
    polarHydrogens,
    radii,
    acceptorAngleAtom,
    acceptorPlaneAtom,
    hbond_std_dev,
    donorAtom,
)

def getSequence(pdb_filename):
    sequence = {}
    for record in SeqIO.parse(pdb_filename + ".pdb","pdb-atom"):
        chain = record.id[-1]
        sequence[chain] = str(record.seq)
    return sequence


# Compute vertex charges based on hydrogen bond potential.
# pdb_filename: The filename of the protonated protein.
# vertices: The surface vertices of the protonated protein
# The name of each vertex in the format, example: B_125_x_ASN_ND2_Green
# where B is chain, 125 res id, x the insertion, ASN aatype, ND2 the name of the
# atom, and green is not used anymore.
def computeFlexibility(pdb_filename, vertices, names):
    '''
    The function takes three parameters: pdb_filename (filename of the protonated protein in PDB format), vertices (surface vertices of the protonated protein), and names (names of each vertex).
    '''
    # The PDB file is parsed using PDBParser from Biopython, and the protein structure is obtained
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure(pdb_filename, pdb_filename + ".pdb")
    sequence = getSequence(pdb_filename)


    # A dictionary (residues) is created to store residues based on their chain and residue ID.
    residues = {}
    for res in struct.get_residues():
        chain_id = res.get_parent().get_id()
        if chain_id == "":
            chain_id = " "
        residues[(chain_id, res.get_id())] = res

    # All atoms are retrieved from the protein structure.
    atoms = Selection.unfold_entities(struct, "A")
    # Satisfied backbone C=O:H-N pairs are computed using the helper function computeSatisfied_CO_HN
    satisfied_CO, satisfied_HN = computeSatisfied_CO_HN(atoms)

    # An array (charge) is initialized to store charges for each vertex.
    charge = np.array([0.0] * len(vertices))
    
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
        # Ignore atom if it is BB and it is already satisfied.
        if atom_name == "H" and res_id in satisfied_HN:
            continue
        if atom_name == "O" and res_id in satisfied_CO:
            continue
        # Compute the charge of the vertex
        charge[ix] = computeChargeHelper(
            atom_name, residues[(chain_id, res_id)], vertices[ix]
        )

    return charge


# Compute the charge of a vertex in a residue.
def computeChargeHelper(atom_name, res, v):
    '''
    atom_name: The name of the atom for which the charge is being computed.
    res: The residue containing the atom.
    v: The coordinates of the vertex (presumably in 3D space).
    '''
    res_type = res.get_resname()
    # Check if it is a polar hydrogen: if it is, it calculates the angle deviation between the donor atom (a) and the hydrogen atom (b) with the vertex (v). 
    # It then calculates the penalty for this angle deviation using computeAnglePenalty. The result is multiplied by 1.0, and this value is returned as the charge.
    if isPolarHydrogen(atom_name, res):
        donor_atom_name = donorAtom[atom_name]
        a = res[donor_atom_name].get_coord()  # N/O
        b = res[atom_name].get_coord()  # H
        # Donor-H is always 180.0 degrees, = pi
        angle_deviation = computeAngleDeviation(a, b, v, np.pi)
        angle_penalty = computeAnglePenalty(angle_deviation)
        return 1.0 * angle_penalty
    # Check if it is an acceptor oxygen or nitrogen
    elif isAcceptorAtom(atom_name, res): # --> If the atom is not a polar hydrogen, it checks if it is an acceptor oxygen or nitrogen by calling the isAcceptorAtom function. 
        acceptor_atom = res[atom_name]
        b = acceptor_atom.get_coord()
        try:
            a = res[acceptorAngleAtom[atom_name]].get_coord()
        except:
            return 0.0
        # 120 degress for acceptor
        angle_deviation = computeAngleDeviation(a, b, v, 2 * np.pi / 3)
        # TODO: This should not be 120 for all atoms, i.e. for HIS it should be
        #       ~125.0
        angle_penalty = computeAnglePenalty(angle_deviation)
        plane_penalty = 1.0
        if atom_name in acceptorPlaneAtom:
            try:
                d = res[acceptorPlaneAtom[atom_name]].get_coord()
            except:
                return 0.0
            plane_deviation = computePlaneDeviation(d, a, b, v)
            plane_penalty = computeAnglePenalty(plane_deviation)
        return -1.0 * angle_penalty * plane_penalty # The final charge is a combination of angle and plane penalties, and it is multiplied by -1.0.
        # Compute the
    return 0.0 # If the atom doesn't fall into the categories above, the function returns 0.0 as the charge.


# Compute the absolute value of the deviation from theta
def computeAngleDeviation(a, b, c, theta):
    return abs(calc_angle(Vector(a), Vector(b), Vector(c)) - theta)


# Compute the angle deviation from a plane
def computePlaneDeviation(a, b, c, d):
    dih = calc_dihedral(Vector(a), Vector(b), Vector(c), Vector(d))
    dev1 = abs(dih)
    dev2 = np.pi - abs(dih)
    return min(dev1, dev2)


# angle_deviation from ideal value. TODO: do a more data-based solution
def computeAnglePenalty(angle_deviation):
    # Standard deviation: hbond_std_dev
    return max(0.0, 1.0 - (angle_deviation / (hbond_std_dev)) ** 2)


def isPolarHydrogen(atom_name, res):
    if atom_name in polarHydrogens[res.get_resname()]:
        return True
    else:
        return False


def isAcceptorAtom(atom_name, res):
    if atom_name.startswith("O"):
        return True
    else:
        if res.get_resname() == "HIS":
            if atom_name == "ND1" and "HD1" not in res:
                return True
            if atom_name == "NE2" and "HE2" not in res:
                return True
    return False


# Compute the list of backbone C=O:H-N that are satisfied. These will be ignored.
def computeSatisfied_CO_HN(atoms):

    # atoms --> list of atoms
    ns = NeighborSearch(atoms) # object for efficient spatial queries on the atoms
    # initialization of two sets to store satisfied C=O and H-N pairs
    satisfied_CO = set()
    satisfied_HN = set()
    for atom1 in atoms: # iterates over each atom in the list of atoms
        res1 = atom1.get_parent()
        if atom1.get_id() == "O": # If the atom is an oxygen (O), it uses NeighborSearch to find neighboring atoms within a distance of 2.5 Angstroms
            neigh_atoms = ns.search(atom1.get_coord(), 2.5, level="A")
            for atom2 in neigh_atoms: # It iterates over the neighboring atoms (atom2) and checks if the neighboring atom is a hydrogen (H).
                if atom2.get_id() == "H":
                    res2 = atom2.get_parent()
                    # Ensure they belong to different residues.
                    if res2.get_id() != res1.get_id():
                        # Compute the angle N-H:O, ideal value is 180 (but in
                        # helices it is typically 160) 180 +-30 = pi
                        # It computes the angles angle_N_H_O_dev and angle_H_O_C_dev representing the deviations from the ideal angles.
                        angle_N_H_O_dev = computeAngleDeviation(
                            res2["N"].get_coord(),
                            atom2.get_coord(),
                            atom1.get_coord(),
                            np.pi,
                        )
                        # Compute angle H:O=C, ideal value is ~160 +- 20 = 8*pi/9
                        # It checks whether these deviations are within the allowed limits (30 degrees for angle_N_H_O_dev and 20 degrees for angle_H_O_C_dev).
                        angle_H_O_C_dev = computeAngleDeviation(
                            atom2.get_coord(),
                            atom1.get_coord(),
                            res1["C"].get_coord(),
                            8 * np.pi / 9,
                        )
                        ## Allowed deviations: 30 degrees (pi/6) and 20 degrees
                        #       (pi/9)
                        # If the deviations are within limits, the residue IDs of the oxygen and hydrogen are added to the corresponding sets.
                        if (
                            angle_N_H_O_dev - np.pi / 6 < 0
                            and angle_H_O_C_dev - np.pi / 9 < 0.0
                        ):
                            satisfied_CO.add(res1.get_id())
                            satisfied_HN.add(res2.get_id())
    # The function returns the sets of satisfied C=O and H-N pairs.
    return satisfied_CO, satisfied_HN


# Compute the charge of a new mesh, based on the charge of an old mesh.
# Use the top vertex in distance, for now (later this should be smoothed over 3
# or 4 vertices)
def assignChargesToNewMesh(new_vertices, old_vertices, old_charges, seeder_opts):
    dataset = old_vertices
    testset = new_vertices
    new_charges = np.zeros(len(new_vertices))

    # The function uses a condition based on seeder_opts["feature_interpolation"]. 
    # If this option is True, it employs feature interpolation; otherwise, it directly assigns charges based on the nearest old vertices.
    if seeder_opts["feature_interpolation"]:
        num_inter = 4  # Number of interpolation features
        # The code uses a k-d tree (KDTree) to find the num_inter nearest neighbors for each new vertex in the old vertices dataset.
        # Assign k old vertices to each new vertex.
        kdt = KDTree(dataset)
        dists, result = kdt.query(testset, k=num_inter)
        # Square the distances (as in the original pyflann) --> The distances are squared, and for each new vertex, the charges are interpolated based on the inverse distances of the nearby old vertices.
        dists = np.square(dists)
        # The size of result is the same as new_vertices
        for vi_new in range(len(result)):
            vi_old = result[vi_new]
            dist_old = dists[vi_new]
            # If one vertex is right on top, ignore the rest.
            if dist_old[0] == 0.0:
                new_charges[vi_new] = old_charges[vi_old[0]]
                continue

            total_dist = np.sum(1 / dist_old)
            for i in range(num_inter):
                new_charges[vi_new] += (
                    old_charges[vi_old[i]] * (1 / dist_old[i]) / total_dist
                )
    else:
        # Assign k old vertices to each new vertex.
        kdt = KDTree(dataset)
        dists, result = kdt.query(testset)
        # The charges for new vertices are directly assigned based on the charges of the nearest old vertices.
        new_charges = old_charges[result]
    return new_charges

