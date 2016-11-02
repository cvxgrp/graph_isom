from util import graph_examples
from util import util
import solver_admm

A = graph_examples.read_petersen()
B = util.random_permute(A)

Z = solver_admm.solve(A, B, 10, verbose=True)
if (Z is not None):
    # print solution
    print Z
