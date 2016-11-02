import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from util import graph_read
from util import graph_mask
from util import util

import solver_admm

for i in [10, 21, 36, 55, 105, 136]:#, 171, 253, 351, 528, 595, 741, 1081]:
	A = graph_read.read_AG2(i)
	B = util.random_permute(A)

	print "------------------------------------"
	print "AG2 graph read from file with n = ", A.shape[0]

	solver_admm.solve(A, B, 10)