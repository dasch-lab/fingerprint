import os
import numpy
from subprocess import Popen, PIPE
import pymesh


# from default_config.global_vars import apbs_bin, pdb2pqr_bin, multivalue_bin
from source.default_config.global_vars import apbs_bin, multivalue_bin
# pdb2pqr bin not found --> installed with pip install pdb2pqr
import random
import pdb2pqr as pdb2pqr_bin

"""
computeAPBS.py: Wrapper function to compute the Poisson Boltzmann electrostatics for a surface using APBS.
Pablo Gainza - LPDI STI EPFL 2019
This file is part of MaSIF.
Released under an Apache License 2.0
"""


def computeAPBS(vertices, pdb_file, tmp_file_base):
    """
        Calls APBS, pdb2pqr, and multivalue and returns the charges per vertex
    """
    # ---------- PD2PQR ----------
    filename_base = pdb_file.split("/")[-1].split(".")[0]
    pdb_id = filename_base.split("_")[0]
    directory = '/disk1/fingerprint/data_preparation/00-raw_pdbs/'
    args = [
        "python",
        "-m",
        "pdb2pqr",
        "--ff=PARSE",
        "--whitespace",
        "--noopt",
        "--apbs-input",
        f"{filename_base}.in",
        f"{pdb_id}.pdb",
        filename_base
    ]
    print(f"Running:\n{' '.join(args)}")
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)
    stdout, stderr = p2.communicate()
    print(stderr.decode("utf-8"))
    print(stdout.decode("utf-8"))

    # ---------- APBS ----------
    args = [apbs_bin, f"{filename_base}.in", '--output-file=test.out']
    print(f"Running:\n{' '.join(args)}")
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)
    stdout, stderr = p2.communicate()

    print(stderr.decode("utf-8"))
    print(stdout.decode("utf-8"))

    print(f'Saving vertices to file {os.path.join(directory, filename_base + ".csv")}')
    vertfile = open(os.path.join(directory, filename_base + ".csv"), "w")
    for vert in vertices:
        vertfile.write("{},{},{}\n".format(vert[0], vert[1], vert[2]))
    vertfile.close()

    # ---------- MULTIVALUE ----------
    args = [
        multivalue_bin,
        filename_base + ".csv",
        filename_base + ".dx",
        filename_base + "_out.csv",
    ]
    print(f"Running:\n{' '.join(args)}")
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)
    stdout, stderr = p2.communicate()

    print(stderr.decode("utf-8"))
    print(stdout.decode("utf-8"))

    # Read the charge file
    chargefile = open(os.path.join(directory, filename_base + "_out.csv"))
    charges = numpy.array([0.0] * len(vertices))
    for ix, line in enumerate(chargefile.readlines()):
        charges[ix] = float(line.split(",")[3])

    # ---------- CLEANUP ----------
    # remove_fn = os.path.join(directory, filename_base)
    # print(f"Removing files: {remove_fn}")
    # os.remove(remove_fn)
    # os.remove(remove_fn + '.csv')
    # os.remove(remove_fn + '.dx')
    # os.remove(remove_fn + '.in')
    # os.remove(remove_fn + '-input.p')
    # os.remove(remove_fn + '_out.csv')

    return charges