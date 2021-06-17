#   Copyright [2021] [Jan Dorazil]
#  #
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#  #
#       http://www.apache.org/licenses/LICENSE-2.0
#  #
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import numpy as np
from scipy.sparse import load_npz, diags
from sor_cy import sor
from scipy.sparse.linalg import cg
from timeit import timeit

A = load_npz('matrices/A.npz')
b = np.load('matrices/b.npy')

# SOR
print('SOR: {:.3f}s'.format(timeit('sor(A, b, omega=1.75, tol=1e-4)', globals=globals(), number=10)))

# CG with Jacobi preconditioning
M = diags([1 / A.diagonal()], [0], format='csr')
print('PCG: {:.3f}s'.format(timeit('cg(A, b, M=M, tol=1e-4)', globals=globals(), number=10)))
