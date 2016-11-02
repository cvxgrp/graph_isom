import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from util import graph_read
from util import graph_mask
from util import util

import solver_admm

for i in [26, 42]:
	A = graph_read.read_PG2(i)
	B = util.random_permute(A)

	print "------------------------------------"
	print "PG2 graph read from file with n = ", A.shape[0]

	solver_admm.solve(A, B, 10)