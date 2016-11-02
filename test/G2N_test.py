import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from util import graph_read
from util import graph_mask
from util import util

import solver_admm

for i in [16, 36, 64, 81, 100, 196, 400, 576, 784, 1024]:
	A = graph_read.read_G2N(i)
	B = util.random_permute(A)

	print "------------------------------------"
	print "G2N graph read from file with n = ", A.shape[0]

	solver_admm.solve(A, B, 10)