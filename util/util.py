import numpy as np


# Finds the block sizes of eigenvalues. It assumes that
# the input is soted and decreasing.
#
# Input: eigvals is a vector of sorted eigenvalues
#		 accuracy is the threshold that we consider two value the same
# Returns: an integer vector denoting the size of blocks
def find_block_sizes(eigvals, accuracy=1e-4):
    block_sizes = np.array(1)
    l = 0  # length-1
    for i in range(len(eigvals) - 1):
        if (eigvals[i] - eigvals[i + 1] < accuracy):  # correct later
            if l == 0:
                block_sizes += 1
            else:
                block_sizes[l] += 1
        else:
            block_sizes = np.append(block_sizes, 1)
            l += 1
    return block_sizes


def eigenvalue_info(A):
    eigenvalues_A, eigenvectors_A = np.linalg.eigh(A)
    idx = eigenvalues_A.argsort()[::-1]
    eigenvalues_A = eigenvalues_A[idx]
    eigenvectors_A = eigenvectors_A[:, idx]
    block_sizes = find_block_sizes(eigenvalues_A)
    return eigenvalues_A, eigenvectors_A, block_sizes


def check_eigs(eigenvalues_A, eigenvalues_B, accuracy=1e-4):
    if np.linalg.norm(eigenvalues_A - eigenvalues_B) > accuracy:
        return -1
    return 0


def random_permute(A):
    n = A.shape[0]
    pi = np.random.permutation(n)
    B = np.zeros([n, n])
    for i in range(n):
        for j in range(n):
            B[i, j] = A[pi[i], pi[j]]
    return B.astype(int)
