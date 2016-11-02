import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from util import graph_read
from util import graph_mask
from util import util

import solver_admm

for i in [32, 50, 72, 98, 128, 242, 392, 578, 800, 1058]:
	A = graph_read.read_KEF(i)
	B = util.random_permute(A)

	print "------------------------------------"
	print "KEF graph read from file with n = ", A.shape[0]

	solver_admm.solve(A, B, 10)