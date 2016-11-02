import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from util import graph_read
from util import graph_mask
from util import util

import solver_admm

for i in [20, 40, 60, 80, 100, 200, 400, 600, 800, 1000]:
	A = graph_read.read_R1N(i)
	B = util.random_permute(A)

	print "------------------------------------"
	print "R1N graph read from file with n = ", A.shape[0]

	solver_admm.solve(A,B,100)
