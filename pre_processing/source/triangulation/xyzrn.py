from Bio.PDB import *
from source.default_config.chemistry import radii, polarHydrogens

"""
xyzrn.py: Read a pdb file and output it is in xyzrn for use in MSMS
Pablo Gainza - LPDI STI EPFL 2019
This file is part of MaSIF.
Released under an Apache License 2.0
"""

#  The XYZRN format includes information about the atomic coordinates, radii, and additional properties 
# required for input to programs like MSMS (Michel Sanner's Molecular Surface).

def output_pdb_as_xyzrn(pdbfilename, xyzrnfilename):
    """
        pdbfilename: input pdb filename
        xyzrnfilename: output in xyzrn format.
    """
    parser = PDBParser()
    struct = parser.get_structure(pdbfilename, pdbfilename)
    outfile = open(xyzrnfilename, "w")

    # The script iterates over atoms in the structure. It extracts information such as atom name, residue, residue name, residue key, chain, and atom type.
    for atom in struct.get_atoms():
        name = atom.get_name()
        residue = atom.get_parent()
        # Ignore hetatms.
        if residue.get_id()[0] != " ":
            continue
        resname = residue.get_resname()
        reskey = residue.get_id()[1]
        chain = residue.get_parent().get_id()
        atomtype = name[0]

        # It assigns a default color ("Green") and checks for specific conditions to modify the color based on atom type and whether it is a polar hydrogen. 
        # It also extracts and formats the atomic coordinates.
        color = "Green"
        coords = None

        # This conditional statement checks if the atomtype (the first character of the atom name) is present in the radii dictionary and if the resname (residue name) 
        # is present in the polarHydrogens dictionary. This check is important because not all atoms may have associated radius or polarity information.
        if atomtype in radii and resname in polarHydrogens:
            if atomtype == "O":
                color = "Red"
            if atomtype == "N":
                color = "Blue"
            if atomtype == "H":
                if name in polarHydrogens[resname]:
                    color = "Blue"  # Polar hydrogens
            # If the conditions are met, the coords variable is assigned a formatted string containing the atomic coordinates. 
            # The get_coord() method retrieves the X, Y, and Z coordinates of the atom
            coords = "{:.06f} {:.06f} {:.06f}".format(
                atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]
            )

            # The insertion variable is assigned the insertion code of the residue. If the insertion code is not a space, it 
            # is assigned the actual insertion code. The full_id variable is then assigned a formatted string containing various 
            # identifiers, including chain, residue number, insertion code, residue name, atom name, and color.
            insertion = "x"
            if residue.get_id()[2] != " ":
                insertion = residue.get_id()[2]
            full_id = "{}_{:d}_{}_{}_{}_{}".format(
                chain, residue.get_id()[1], insertion, resname, name, color
            )

        # If valid coordinates are obtained, it writes the information to the output file in the XYZRN format, 
        # including the coordinates, atomic radius, a count (1), and a full identifier.
        if coords is not None:
            outfile.write(coords + " " + radii[atomtype] + " 1 " + full_id + "\n")

