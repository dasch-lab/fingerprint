"""
protonate.py: Wrapper method for the reduce program: protonate (i.e., add hydrogens) a pdb using reduce 
                and save to an output file.
Pablo Gainza - LPDI STI EPFL 2019
Released under an Apache License 2.0
"""

# Comment SJ:
# 1) By removing existing protons first, you ensure that you are starting with a clean slate. This is important because the original PDB 
# file might already contain hydrogen atoms, and you want to avoid any potential inconsistencies or conflicts when adding new hydrogens.
# 2) This helps ensure that the resulting protonation is consistent and meets the criteria or assumptions of the downstream analysis.
# 3) If you skip this step and directly add hydrogens using -HIS, you might end up with redundant or conflicting hydrogen atom positions.

from subprocess import Popen, PIPE
from IPython.core.debugger import set_trace
import os


def protonate(in_pdb_file, out_pdb_file):
    # protonate (i.e., add hydrogens) a pdb using reduce and save to an output file.
    # in_pdb_file: file to protonate.
    # out_pdb_file: output file where to save the protonated pdb file. 
    
    # Remove protons first, in case the structure is already protonated
    # Use Propen class to execute the 'reduce' program with the options '-Trim' and the input PDB file ('in_pdb_file')
    # The '-Trim' option removes existing protons from the structure
    # The output from the execution is captured in 'stdout' and 'stderr'
    args = ["reduce", "-Trim", in_pdb_file] 
    p2 = Popen(args, stdout=PIPE, stderr=PIPE) 
    stdout, stderr = p2.communicate()

    # The stdout output from the reduce command is written to the output PDB file (out_pdb_file) after decoding it from 
    # bytes to a UTF-8 string and removing trailing whitespaces.
    outfile = open(out_pdb_file, "w")
    outfile.write(stdout.decode('utf-8').rstrip())
    outfile.close()

    # Now add them again: After removing existing protons, this section uses reduce again, but this time with the -HIS option. 
    # This step adds hydrogens to the structure.
    args = ["reduce", "-HIS", out_pdb_file]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p2.communicate()

    # The stdout output from this second execution of reduce is again written to the output PDB file. The file is opened in 
    # write mode ("w"), and the decoded output is written.
    outfile = open(out_pdb_file, "w")
    outfile.write(stdout.decode('utf-8'))
    outfile.close()

