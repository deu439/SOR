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

from scipy.sparse import isspmatrix_csr, SparseEfficiencyWarning, csr_matrix
from warnings import warn
import numpy as np
import cython
from cython cimport floating


def sor(A, floating [:] b, double omega, double tol=1e-9, int maxiter=-1, x0=None):
    assert(A.shape[0] == A.shape[1] == len(b))
    assert(0 < omega < 2)
    cdef int ln = len(b)
    
    if floating is float:
        dtype = np.float32
    else:
        dtype = np.float64

    # Check that A is CSR
    if not isspmatrix_csr(A):
        A = csr_matrix(A)
        warn('SOR requires A be CSR matrix format', SparseEfficiencyWarning)

    # Construct array views
    cdef floating [:] data = A.data
    cdef int [:] indices = A.indices
    cdef int [:] indptr = A.indptr

    # Initialize the solution vector
    cdef floating [:] x = np.zeros(ln, dtype=dtype)
    if x0 is not None:
        x = x0

    # Start iterations
    cdef floating bnorm = np.linalg.norm(b)
    cdef floating res = np.linalg.norm(A @ x - b)
    cdef bint test = True
    cdef int n = 0

    while test:
        n += 1
        x = c_sor(data, indices, indptr, b, x, omega, ln)

        if maxiter > 0 and n >= maxiter:
            test = False
        if tol > 0:
            res = np.linalg.norm(A @ x - b)
            if res < tol*bnorm:
                test = False

    return x, n


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
@cython.profile(True)

cdef floating [:] c_sor(floating [:] data, int [:] indices, int [:] indptr, floating [:] b, floating [:] x, double omega, int ln) nogil:
    cdef int i
    cdef int j
    cdef floating sigma
    cdef floating aii
    cdef int k

    for i in range(ln):  # Row index
        sigma = 0
        aii = 0
        for j in range(indptr[i], indptr[i+1]):  # data array index
            k = indices[j]  # column index
            if k != i:
                sigma += data[j] * x[k]
            else:
                aii = data[j]

        x[i] = x[i] + omega*((b[i] - sigma)/aii - x[i])

    return x
