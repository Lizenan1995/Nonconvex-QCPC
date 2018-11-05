This is an example code for non-convex QCQP problem.

The optimization problem is shown as follows:

min  x^TPx + p^T x + r

s.t.   x^TQx + q^T x + c <= 0

where inequality constraints can be more than one, and equality constraints can be written as inequality constraints.

The `document/` is the optimization algorithm to solve this problem.

The `code/` is matlab code reproduce.