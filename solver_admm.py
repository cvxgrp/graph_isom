import numpy as np
from munkres import Munkres
from util import util
from util import graph_mask

np.random.seed(0)

def find_projection_mask(block_sizes):
	n = sum(block_sizes)
	ind = np.append(0, np.cumsum(block_sizes))
	M = np.zeros([n, n])
	for i in range(len(block_sizes)):
		M[ind[i]:ind[i + 1], ind[i]:ind[i + 1]] = 1
	return M


def spectral_project(Ptilde, eigenvectors_A, eigenvectors_B, projection_mask):
	pre_project = np.dot(np.dot(eigenvectors_B.T, Ptilde), eigenvectors_A)
	post_project = pre_project * projection_mask
	return np.dot(np.dot(eigenvectors_B, post_project), eigenvectors_A.T)


# this is hungarian algorithm
def hungarian_project(matrix):
	n = matrix.shape[0]
	m = Munkres()
	indexes = m.compute(-matrix)
	result = np.zeros([n, n])
	for row, column in indexes:
		result[row, column] = 1
	return result


def solve_basic(W, eigenvectors_A, eigenvectors_B, projection_mask,
				sparsity_mask=None, max_iter=100):
	n = eigenvectors_A.shape[0]
	if sparsity_mask is None:
		sparsity_mask = np.ones([n, n])

	P1 = np.zeros([n, n])
	P2 = np.zeros([n, n])
	P3 = np.zeros([n, n])
	Z = np.ones([n, n]) / n
	Y1 = np.zeros([n, n])
	Y2 = np.zeros([n, n])
	Y3 = np.zeros([n, n])

	W = W / 2 / (sum(sum(np.abs(W))) / n ** 2)

	for iter in range(max_iter):
		P1 = Z - W - Y1 - np.outer(np.dot(Z - W - Y1, np.ones(n)) -
									np.ones(n), np.ones(n).T) / n
		P2 = Z - W - Y2 - np.outer(np.ones(n), np.dot(np.ones(n).T, Z - W - Y2) -
									np.ones(n).T) / n
		P3 = np.maximum(Z - Y3, 0) * sparsity_mask
		Ptilde = (P1 + P2 + P3) / 3.0

		Z = spectral_project(Ptilde, eigenvectors_A, eigenvectors_B, projection_mask)

		Y1 += (P1 - Z)
		Y2 += (P2 - Z)
		Y3 += (P3 - Z)
	return Z


def solve(A, B, max_rep, verbose=True, use_hungarian=True, TOL=1e-4):
	eigenvalues_A, eigenvectors_A, block_sizes = util.eigenvalue_info(A)
	eigenvalues_B, eigenvectors_B, block_sizes = util.eigenvalue_info(B)

	if util.check_eigs(eigenvalues_A, eigenvalues_B):
		print "Two matrices have different spectra, hence not isomorphic."
		return None

	M1 = graph_mask.create_degree_invariants(A, B, A.shape[0])
	M2 = graph_mask.create_spectral_invariants(eigenvectors_A,
					eigenvectors_B, block_sizes, diagons=1, row_sums=1)
	mask = M1 * M2

	n = A.shape[0]
	if verbose:
		"Sparsity mask found. nnz = ", np.round(np.sum(mask) / (n ** 2), 3)
	graph_mask.prune(mask, A, B)
	projection_mask = find_projection_mask(block_sizes)

	if verbose:
		# should fix this
		row_format = "{:>15}" * 5
		print row_format.format("instance", "||ZA-BZ||", "max(|Z1-1|)",
								"max(|Z.T1-1|)", "max(Z(1-Z))")
	for init in range(max_rep):
		W = np.random.randn(n, n)
		Z = solve_basic(W, eigenvectors_A, eigenvectors_B,
						projection_mask, sparsity_mask=mask)
		if use_hungarian:
			Z = hungarian_project(Z)
		r1 = np.linalg.norm(np.dot(Z, A) - np.dot(B, Z))
		r2 = np.max(np.abs(np.dot(Z, np.ones(n)) - np.ones(n)))
		r3 = np.max(np.abs(np.dot(Z.T, np.ones(n)) - np.ones(n)))
		r4 = np.max(Z * (1 - Z))

		if verbose:
			print row_format.format(init + 1, r1, r2, r3, r4)

		if ((r1 < TOL) and (r2 < TOL) and (r3 < TOL) and (r4 < TOL)):
			if verbose:
				print "Isomorphy found in ", init + 1, " instances."
			return Z

	if verbose:
		print "No isomorphy found in ", max_rep, " instances."
	return None
