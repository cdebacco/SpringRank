"""
New version developed by Nicolò Ruggeri, Max Planck Institute for Intelligent Systems, Tuebingen, Germany, March-2020
It forces to use sparse matrices when possible, results in much more efficent implementation, especially for large matrices
"""

import warnings

import scipy.sparse

from ._core import build_from_dense, build_from_sparse, solve_linear_system


def SpringRank(A, alpha=0., l0=1., l1=1., solver='bicgstab', verbose=False, force_dense=False):
    """
        Main routine to calculate SpringRank by a solving linear system.

        Parameters
        ----------
        A : numpy.ndarray or scipy.sparse.spmatrix
            Has tobe  2 dimensional and with same dimensions.
        alpha, l0, l1: float
            Defined as in the SpringRank paper
            https://arxiv.org/abs/1709.09002
        solver: str
            One between 'spsolve' (direct, slower) and 'bicgstab' (iterative, faster).
            The solver to be used for the linear system returning the ranks.
        verbose: bool
        force_dense: bool
            By default A is converted to a sparse matrix scipy.sparse.csr, if it is not already sparse.
            If force_dense is set to True and a dense ndarray A is input, then it is not converted to sparse.

        Returns
        -------
        rank
            numpy.ndarray of ranks. Indices represent the nodes' indices used in the matrix A.

    """

    # check if input is sparse or can be converted to sparse.
    use_sparse = True
    if force_dense and not scipy.sparse.issparse(A):
        try:
            A = scipy.sparse.csr_matrix(A)
        except:
            warnings.warn('The input parameter A could not be converted to scipy.sparse.csr_matrix. '
                          'Using a dense representation.')
            use_sparse = False
    elif force_dense:
        use_sparse = False

    # build array to feed linear system solver
    if use_sparse:
        A, B = build_from_sparse(A, alpha, l0, l1)
    else:
        A, B = build_from_dense(A, alpha, l0, l1)

    rank = solve_linear_system(A, B, solver, verbose)

    return rank


