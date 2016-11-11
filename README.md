# graph_isom

A heuristic for the graph isomorphism problem which is based on linear programming. For more information, please see [our associated paper](http://stanford.edu/~boyd/papers/graph_isom.html).

For given symmetric matrices A,B the following line runs the algorithm for N random instances.
```python
Z = solver_admm.solve(A, B, N)
```
If no solution is found Z=None will be returned. For more examples, see test subdirectory.
