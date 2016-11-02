import numpy as np
import warnings


def round(P, threshold=1e-2):
    for i in range(P.shape[0]):
        for j in range(P.shape[1]):
            if P[i, j] < threshold:
                P[i, j] = 0
            elif P[i, j] > 1 - threshold:
                P[i, j] = 1
    if (np.sum(P * (1 - P))):
        warnings.warn("The final solution has entries in [eps, 1-eps].")
    return np.round(P, 2)


def print_header():
    header = ["init", "big entries", "maxP(1-P)", "norm(P1-1)", "max(P.T1)",
              "sum|PART-B|"]
    row_format = "{:>15}" * (len(header))
    print row_format.format(*header)


def print_error(P, A, B, init=""):
    n = A.shape[0]
    row_format = "{:>15}" * 6

    big_ent = np.sum(P > .3)
    max_dev = np.max(np.multiply(P, 1 - P))
    P1dev = np.linalg.norm(np.dot(P, np.ones([n, 1])) - np.ones([n, 1]))
    PT1dev = np.linalg.norm(np.dot(P.T, np.ones([n, 1])) - np.ones([n, 1]))
    APdev = np.linalg.norm(np.dot(P, A) - np.dot(B, P))

    print row_format.format(init, big_ent, np.round(max_dev, 4),
                    np.round(P1dev, 4), np.round(PT1dev, 4), APdev)


def matrix_to_lists(matrix):
    """Convert a matrix to a list of lists.
    """
    rows, cols = matrix.shape
    lists = []
    for i in range(rows):
        lists.append(matrix[i, :].tolist()[0])
    print lists
    return lists


# works for .5 now only
def print_cycles(Z):
    row_format = "{:>5}" * 3
    for i in range(6):
        j = -1
        k = -1
        for l in range(Z.shape[0]):
            if Z[i, l] == .5:
                if j == -1:
                    j = l
                else:
                    k = l
                    break
    print row_format.format(i, j, k)
    print "-----------------"
