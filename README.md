Introduction
----
This is an implementation of the successive over relaxation (SOR) linear system solver, which is well suited for large sparse systems of linear equations. Taking a large sparse matrix $A$ (in <i>scipy.sparse.csr_matrix</i> format) and a (<i>numpy.array</i>) vector $b$ it iteratively approaches the solution  to the $Ax=b$ system. Convergence is ensured for symmetric positive-definite matrices for any $\omega \in (0, 2)$.

[See SOR on Wikipedia](https://en.wikipedia.org/wiki/Successive_over-relaxation).

Instalation
----
You need python and cython installed on your system. You additionally need the numpy, scipy and optionally timeit modules installed. If these requirements are fulfilled just run

	$ python setup.py build_ext --inplace

to compile the code into a .so/.dll dynamic library.

Testing
----
After the compilation is done, you can run

	$ python test_sor.py

Which should output something like

	SOR: 12.350s
	PCG: 15.085s

This means that for the example system in <code>matrices/A.npz</code> and <code>matrices/b.npy</code> the SOR reaches the specified relative residual $1.22 \times$ faster than the conjugate gradient algorithm of <i>scipy.sparse.linalg</i> with Jacobi preconditioning, on my computer.

Usage
----
To use this in your project just copy the generated .so/.dll library into your project and import it via <code>from sor_cy import sor</code>. Now you can call <code>sor(A, b, omega, tol=1e-9, maxiter=-1, x0=None)</code>. The required arguments are the system matrix $A$, the right-hand-side vector $b$ and the relaxation parameter $\omega\in (0, 2)$ ('omega'). The optional arguments are 'tol': the relative tolerance of the residual error, 'maxiter': the maximum number of iterations and 'x0': an initial approximation to the solution vector. Setting 'maxiter' to $-1$ means that the algorithm stops when $|Ax-b|_2 < tol\cdot|b|_2$. Setting 'tol' to $-1$ means that the algorithm stops after 'maxiter' iterations. If none of these parameters is set to $-1$, the algorithm stops as soon as one of these conditions is fulfilled.

Note that if 'tol' is not $-1$ the algorithm calculates the residual $|Ax-b|_2$ in each iteration, which might be costly. To obtain a good performance the parameter $\omega$ needs to be tuned for the particular linear system at hand.

Any suggestions on how to improve the code are welcome.