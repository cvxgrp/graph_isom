import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from util import graph_read
from util import graph_mask
from util import util

import solver_admm

for i in [22, 44, 88, 132, 198, 264, 352, 440, 550, 660, 792, 924, 1078]:
	A = graph_read.read_CHH(i)
	B = util.random_permute(A)

	print "------------------------------------"
	print "CHH graph read from file with n = ", A.shape[0]

	solver_admm.solve(A, B, 10)