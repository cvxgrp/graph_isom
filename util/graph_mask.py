import numpy as np


# Creates a sparsity mask based on two input vectors.
# If v_j != w_i, then the (i,j) entery of the out is false (zero)
#
# Input: two vectors of equal length
# Output: the mask made from these two vector
def create_one_mask(v, w):
    n = len(v)
    outer_prod = np.outer(v, np.ones(n)) - np.outer(np.ones(n), w)
    return (np.abs(outer_prod).T < 1e-5)


# Returns the sparsity mask built based on the degrees invariant
#
# Input: adjacency matrices A,B
# Returns: the sparsity mask
def create_degree_invariants(A, B, powers=1):
    n = A.shape[0]
    sparsity_mask = np.ones([n, n])

    dA = np.ones(n)
    dB = np.ones(n)

    for i in range(powers):
        dA = np.dot(A, dA)
        dB = np.dot(B, dB)
        dA /= np.linalg.norm(dA)
        dB /= np.linalg.norm(dB)
        sparsity_mask *= create_one_mask(dA, dB)
    return sparsity_mask


def create_spectral_invariants(eigenvectors_A, eigenvectors_B,
                               block_sizes, diagons=0, row_sums=0):
    n = eigenvectors_A.shape[0]
    sparsity_mask = np.ones([n, n])
    partition_eigvals = np.cumsum(np.append(0, block_sizes))

    if diagons:
        for i in range(len(block_sizes)):
            diag1 = np.sum(eigenvectors_A[:, partition_eigvals[i]:
                        partition_eigvals[i + 1]] ** 2, axis=1)
            diag2 = np.sum(eigenvectors_B[:, partition_eigvals[i]:
                        partition_eigvals[i + 1]] ** 2, axis=1)
            sparsity_mask *= create_one_mask(diag1, diag2)

    if row_sums:
        for i in range(len(block_sizes)):
            vec1 = eigenvectors_A[:,
                    partition_eigvals[i]:partition_eigvals[i + 1]]
            vec2 = eigenvectors_B[:,
                    partition_eigvals[i]:partition_eigvals[i + 1]]
            row_sum1 = np.dot(vec1, np.dot(vec1.T, np.ones(n)))
            row_sum2 = np.dot(vec2, np.dot(vec2.T, np.ones(n)))
            sparsity_mask *= create_one_mask(row_sum1, row_sum2)
    return sparsity_mask


# XXX not final yet
def prune(R, A, B):
    n = A.shape[0]
    for i in range(n):
        for j in range(n):
            if R[j, i] == 1:
                # now i,j are potential
                for k in range(n):
                    if A[i, k] > 0:
                        # now k is a neighbor of i
                        possible = 0
                        for l in range(n):
                            if B[j, l] == A[i, k]:
                                # now l is a neighbor
                                if R[l, k] == 1:
                                    possible = 1
                        if possible == 0:
                            R[i, j] = 0
                            continue


def prune2(R, A, B):
    n = A.shape[0]
    for i in range(n):
        found = 0
        for j in range(n):
            if R[j, i] == 1:
                found += 1
                ind = j
        if found == 1:
            R[ind, :] = 0
            R[ind, i] = 1
    for i in range(n):
        found = 0
        for j in range(n):
            if R[i, j] == 1:
                found += 1
                ind = j
        if found == 1:
            R[:, ind] = 0
            R[i, ind] = 1
