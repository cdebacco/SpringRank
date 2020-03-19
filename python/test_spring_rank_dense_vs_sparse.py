"""
Test developed by Nicolò Ruggeri, Max Planck Institute for Intelligent Systems, Tuebingen, Germany, March-2020
It shows that working with sparse matrices it's much faster, this can be noticed for Adjacency matrices of size ~10^4
"""
import scipy.sparse
import numpy as np
import networkx as nx
import timeit

from SpringRank import SpringRank


# test correctness (outputs are the same for old and new function)
def test_correctness(alpha, rtol=1.e-10):
	print('\nTest correctness with alpha = ', alpha)
	res_dense = SpringRank(B, alpha=alpha, verbose=verbose, force_dense=True)
	res_sparse = SpringRank(A, alpha=alpha, verbose=verbose)
	assert np.allclose(res_dense, res_sparse, rtol=rtol)
	print('Passed.')


# test speed (having the adjacency matrix already instanciated)
def test_speed(alpha, rep):
	# print(f'\nTest speed with alpha = {alpha} and {rep} repetitions')
	print("\nTest speed with alpha={alpha} and {rep} repetitions".format(alpha=alpha, rep=rep))
	sr_dense = lambda: SpringRank(B, alpha=alpha, verbose=False, force_dense=True)
	sr_sparse = lambda: SpringRank(A, alpha=alpha, verbose=False)

	print('Time for SpringRank with dense arrays: ', timeit.timeit(sr_dense, number=rep))
	print('Time for SpringRank with sparse arrays:', timeit.timeit(sr_sparse, number=rep))


if __name__ == '__main__':
	# define experiment instance
	experiment = 4

	if experiment == 1:
		A = scipy.sparse.csr_matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float)
		B = A.toarray()
	elif experiment == 2:
		A = scipy.sparse.csr_matrix(np.random.rand(100, 100))
		B = A.toarray()
	elif experiment == 3:
		G = nx.scale_free_graph(1000)
		A = nx.to_scipy_sparse_matrix(G, dtype=float)
		B = A.toarray()
	elif experiment == 4:
		print('Generating synsthetic scale free graph...')
		G = nx.scale_free_graph(10000)
		print('Done.')
		A = nx.to_scipy_sparse_matrix(G, dtype=float)
		B = A.toarray()

	verbose = False
	solver = 'spsolve'  # 'bicgstab'
	rep = 10

	# test correctness
	test_correctness(alpha=0.)
	test_correctness(alpha=1.)

	# test speed
	test_speed(alpha=0., rep=rep)
	test_speed(alpha=1., rep=rep)
