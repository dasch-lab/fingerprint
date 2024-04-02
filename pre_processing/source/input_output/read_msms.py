import numpy as np
"""
read_msms.py: Read an msms output file that was output by MSMS (MSMS is the program we use to build a surface) 
Pablo Gainza - LPDI STI EPFL 2019
Released under an Apache License 2.0
"""

def read_msms(file_root):
    # read the surface from the msms output. MSMS outputs two files: {file_root}.vert and {file_root}.face
    
    vertfile = open(file_root + ".vert")
    meshdata = (vertfile.read().rstrip()).split("\n")
    vertfile.close()

    # Read number of vertices.
    count = {}
    header = meshdata[2].split()
    count["vertices"] = int(header[0])
    ## Data Structures
    # This section parses the vertex data. It extracts the number of vertices from the file, initializes data structures 
    # (NumPy arrays and lists) to store vertex information, and iterates through the lines of the file to populate the 
    # arrays with vertex coordinates, normals, atom IDs, and residue IDs.
    vertices = np.zeros((count["vertices"], 3))
    normalv = np.zeros((count["vertices"], 3))
    atom_id = [""] * count["vertices"]
    res_id = [""] * count["vertices"]
    for i in range(3, len(meshdata)):
        fields = meshdata[i].split()
        vi = i - 3
        vertices[vi][0] = float(fields[0])
        vertices[vi][1] = float(fields[1])
        vertices[vi][2] = float(fields[2])
        normalv[vi][0] = float(fields[3])
        normalv[vi][1] = float(fields[4])
        normalv[vi][2] = float(fields[5])
        atom_id[vi] = fields[7]
        res_id[vi] = fields[9]
        count["vertices"] -= 1

    # Read faces_ The script opens and reads the {file_root}.face file, storing the content in the meshdata list.
    facefile = open(file_root + ".face")
    meshdata = (facefile.read().rstrip()).split("\n")
    facefile.close()

    # Read number of vertices: This part parses the face data. It extracts the number of faces from the file, 
    # initializes data structures for storing face information, and iterates through the lines of the file to populate the array with face vertex indices.
    header = meshdata[2].split()
    count["faces"] = int(header[0])
    faces = np.zeros((count["faces"], 3), dtype=int)
    normalf = np.zeros((count["faces"], 3))

    for i in range(3, len(meshdata)):
        fi = i - 3
        fields = meshdata[i].split()
        faces[fi][0] = int(fields[0]) - 1
        faces[fi][1] = int(fields[1]) - 1
        faces[fi][2] = int(fields[2]) - 1
        count["faces"] -= 1

    # Finally, the function performs assertions to ensure that all vertices and faces have been processed, 
    # and it returns the vertices, faces, vertex normals, and residue IDs as NumPy arrays.
    assert count["vertices"] == 0
    assert count["faces"] == 0

    return vertices, faces, normalv, res_id

