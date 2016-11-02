import numpy as np
import warnings
from binascii import hexlify


# Reads one word from a binary file and returns it as an integer
# with little-endian convention
#
# Input: f is an open file
# Returns: an integer which is the value of next word
def read_word(f):
    b1 = f.read(1)
    if b1 == "":
        warnings.warn("read_word(): Something might be wrong with reading b1")
    b2 = f.read(1)
    if b2 == "":
        warnings.warn("read_word(): Something might be wrong with reading b2")
    hex = hexlify(b1)
    bb1 = int(hex, 16)
    hex = hexlify(b2)
    bb2 = int(hex, 16)
    return bb2 * (2 ** 8) + bb1


# Reads a complete graph, assuming that it's stored as unlabeled
# Please refer to
# http://mivia.unisa.it/datasets/graph-database/arg-database/documentation/
#
# Input: myfile is a string that contains the address of the file
# Returns: the adjacency matrix.
def read_graph_unlabeled_from_file(myfile):
    with open(myfile, "rb") as f:
        nodes = read_word(f)
        adj_mat = np.zeros([nodes, nodes])
        for i in range(nodes):
            edges = read_word(f)
            for j in range(edges):
                target = read_word(f)
                adj_mat[i][target] = 1
        b1 = f.read(1)
        if b1 != "":
            warnings.warn("Something might be " +
                          "wrong with reading the file as unlabeled")
    return adj_mat


# Reads the adjacency matrix of a graph and verifies that it's a symmetric
# binary matrix.
#
# Input: myfile is a string that contains the address of the file
# Returns: the adjacency matrix (symmetric and verified).
def read_unlabled(myfile):
    A = read_graph_unlabeled_from_file(myfile)
    if np.max(A + A.T) > 1:
        warnings.warn("something unexpected in A = A + A.T. A" +
                      " is probably not antisymmetric.")
    A = A + A.T
    return A


# Reads a graph from the iso library. Assumes that the iso_m2D folder
# is located in the current directory.
#
# Input: size is an integer that denotes the size of the graph
#        from the following numbers
#				(16, 36, 64, 81, 100, 196, 400, 576, 784, 1024)
#		 index is a number between 0 and 99
# Returns: the adjacency matrice of A,B (symmetrix and verified).
def read_iso_m2D(size, index):
    if (size in [16, 36, 64, 81, 100]):
        prefix = 'iso_m2D/iso_m2D_s' + str(size)
    elif (size in [196, 400, 576, 784, 1024]):
        prefix = 'iso_m2D/iso_m2D_m' + str(size)
    else:
        warnings.warn("Graph size is not correct")

    if (index < 10):
        extension = "0" + str(index)
    else:
        extension = str(index)

    filename_A = prefix + ".A" + extension
    filename_B = prefix + ".B" + extension

    return (read_unlabled(filename_A), read_unlabled(filename_B))


def read_graph_hadamard_from_file(myfile):
    with open(myfile, "rb") as f:
        line = f.readline()
        header = line.split()
        if ((len(header) != 4) | (header[0] != "p") | (header[1] != "edge")):
            warnings.warn("first line is not as expected")
        nodes = int(header[2])
        adj_mat = np.zeros([nodes, nodes])
        while True:
            line = f.readline()
            if not line:
                break
            line = line.split()
            if ((line[0] != "e") | (len(line) != 3)):
                warnings.warn("line does not start with e")
            adj_mat[int(line[1]) - 1, int(line[2]) - 1] = 1
            adj_mat[int(line[2]) - 1, int(line[1]) - 1] = 1
    return adj_mat


def read_hadamard(size):
    filename = "had/had-" + str(size)
    return read_graph_hadamard_from_file(filename)


def read_sts(size):
    dict = {57: 19, 100: 25, 155: 31, 222: 37, 301: 43, 392: 49, 495: 55,
            610: 61, 737: 67, 876: 73, 1027: 79}
    if size in dict:
        filename = "sts/sts-" + str(dict[size]) + "_" + str(size)
        print filename
        return read_graph_unlabeled_from_file(filename)
    else:
        warnings.warn("wrong size for sts")


def read_latin(size):
    if size not in np.arange(2, 31):
        warnings.warn("wrong size for latin")
    filename = "latin/latin-" + str(size)
    return read_graph_hadamard_from_file(filename)


def read_lattice(size):
    if size not in [4, 6, 8, 9, 11, 15, 20, 25, 29, 32]:
        warnings.warn("wrong size for LAT")
    filename = "LAT/Lattice(" + str(size) + ")_" + str(size ** 2)
    return read_graph_unlabeled_from_file(filename)


def read_TRI(size):
    dict = {21: 7, 45: 10, 66: 12, 91: 14, 120: 16, 210: 21,
            406: 29, 595: 35, 780: 40, 990: 45}
    if size in dict:
        filename = "data/TRI/Triangular(" + str(dict[size]) + ")_" + str(size)
        return read_graph_unlabeled_from_file(filename)
    else:
        warnings.warn("wrong size for TRI")


def read_R1N(size):
    if size in [20, 40, 60, 80, 100]:
        filename = "data/R1N/iso_r01N_s" + str(size)
    elif size in [200, 400, 600, 800, 1000]:
        filename = "data/R1N/iso_r01N_m" + str(size)
    else:
        warnings.warn("wrong size for R1D")
    return read_graph_unlabeled_from_file(filename).astype(int)


def read_USR(size, one=True):
    dict = {29: 1, 58: 2, 87: 3, 116: 4, 203: 7, 406: 14,
            609: 21, 812: 28, 986: 34}
    if size in dict:
        filename = "data/USR/usr(" + str(dict[size]) + ")_" + str(size)
        if one:
            filename += "-1"
        else:
            filename += "-2"
        return read_graph_unlabeled_from_file(filename)
    else:
        warnings.warn("wrong size for USR")


def read_AG2(size):
    dict = {10: 2, 21: 3, 36: 4, 55: 5, 105: 7, 136: 8, 171: 9, 253: 11,
            351: 13, 528: 16, 595: 17, 741: 19, 1081: 23}
    if size in dict:
        filename = "data/AG2/ag2-" + str(dict[size]) + "_" + str(size)
        return read_graph_unlabeled_from_file(filename)
    else:
        warnings.warn("wrong size for AG2")


def read_G2N(size):
    if size in [16, 36, 64, 81, 100]:
        filename = "data/G2N/iso_m2N_s" + str(size)
    elif size in [196, 400, 576, 784, 1024]:
        filename = "data/G2N/iso_m2N_m" + str(size)
    else:
        warnings.warn("wrong size for G2N")
    return read_graph_unlabeled_from_file(filename)


def read_KEF(size, one=True):
    dict = {32: 4, 50: 5, 72: 6, 98: 7, 128: 8, 242: 11, 392: 14, 578: 17,
            800: 20, 1058: 23}
    if size in dict:
        filename = "data/KEF/kef-" + str(dict[size]) + "_" + str(size)
        if one:
            filename += "-1"
        else:
            filename += "-2"
        return read_graph_unlabeled_from_file(filename)
    else:
        warnings.warn("wrong size for KEF")


def read_CHH(size, one=True):
    dict = {22: "1-1", 44: "2-1", 88: "2-2", 132: "3-2", 198: "3-3",
            264: "4-3", 352: "4-4", 440: "5-4", 550: "5-5",
            660: "6-5", 792: "6-6", 924: "7-6", 1078: "7-7"}
    if size in dict:
        filename = "data/CHH/CHH_cc(" + dict[size] + ")_" + str(size)
        if one:
            filename += "-1"
        else:
            filename += "-2"
        return read_graph_unlabeled_from_file(filename)
    else:
        warnings.warn("wrong size for CHH")


def read_MZA(size):
    filename = "data/MZA/mz-aug-" + str(size / 20) + "_" + str(size)
    return read_graph_unlabeled_from_file(filename)


def read_PG2(size):
    dict = {14: 2, 26: 3, 42: 4, 62: 5, 114: 7, 146: 8, 182: 9, 266: 11,
            366: 13, 546: 16, 614: 17, 762: 19, 1106: 23}
    if size in dict:
        filename = "data/PG2/PG2(" + str(dict[size]) + ")_" + str(size)
        return read_graph_unlabeled_from_file(filename)
    else:
        warnings.warn("wrong size for AG2")
