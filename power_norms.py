"""
Functions to estimate (rigorously) norms of powers of a matrix, given a nearby matrix.
"""

from numpy import zeros, ones, inf, minimum
from sage.all import RealField, RR, RealNumber, vector, matrix, identity_matrix

from sparse import norm1, gamma, max_nonzeros_per_row, matrix_norm1, matrix_norminf, sage_sparse_to_scipy_sparse, norminf

import unittest

def basis_of_zero_sum_subspace(n):
	"""
	Gives a matrix V whose column space is the zero-sum subspace.

	Used only for testing purposes.

	"""
	V = -identity_matrix(n)
	V = V[:, 1:]
	V[0,:] = 1.
	return V

def norm1_of_powers(P, M, m=20, K=RealField(rnd='RNDU').zero()):
	"""
	Estimates :math:`\| A^k |V \|_1` for k in range(m).

	:math:`A` is an operator such that :math:`AV \subseteq V`
	:math:`|V` means that the norm is restricted to the set of zero-sum vectors.

	Args:
		P (scipy matrix): a matrix such that :math:`\|P-A\| \leq K`
		m (integer): how many norms to compute
		M, K (reals with RNDU): constants such that :math:`\|P-A\| \leq K` and :math:`\|A^i\| \leq M` for each i in `range(m)`.
	
	Returns:
		a lists of norm bounds, starting from `k=0` (which is simply 1)
    """

	n = P.shape[1]
	Ci = zeros(m)

	prop = (gamma(max_nonzeros_per_row(P)) * matrix_norm1(P, rnd='RNDU') + K) * M
	zero = RealField(rnd='RNDU').zero()

	Ci[0] = 1.
	for j in range(1, n):
		# We test the norm on these vectors: v_j = [1, 0, 0, ... , 0, -1 (in position j), 0, ... , 0]
		# The matrix that has them as columns satisfies ||Vx|| \geq ||x||.
		# So for each Vx \in span V, we have ||P^i V x|| / ||Vx|| \leq ||P^i V||*||x|| / ||x|| = ||P^i V||.

		v = zeros(n)
		v[0] = 1.
		v[j] = -1.

		# Error propagation: we use the bound
		#||(F^k-A^k)v|| \leq \sum ||A^(j-i-1) (F-A) F^i v|| \leq M ||F-A|| \sum ||v_i||,
		# where F is the operator that computes the product in floating point
		current_norm = 2.

		error = zero
		for i in range(1,m):
			v = P * v
			error += prop * current_norm
			current_norm = norm1(v, rnd='RNDU')
			Ci[i] = max(Ci[i], current_norm + error)

	for i in range(1, m): #if the bounds are larger than M, replace them with M
		Ci[i] = min(Ci[i], M)

	return Ci

def norm1_of_powers_alternate(A, M, m=20):
	"""
	Estimates :math:`\| A^k |V \|_1` for k in range(m), using an alternate algorithm

	:math:`A` is an operator such that :math:`AV \subseteq V`
	:math:`|V` means that the norm is restricted to the set of zero-sum vectors.

	Args:
		A (Sage interval matrix): a matrix :math:`ones((1,n))*A = ones((1,n))`
		m (integer): how many norms to compute
		M (real with RNDU): constant such that :math:`\|A^i\| \leq M` for each i in `range(m)`.
	
	Returns:
		a lists of norm bounds, starting from `k=0` (which is simply 1)

	Uses the method in our notes `algorithm_for_computing_norms.pdf`.
    """

    # we take the smallest matrix in the interval A, since we want a lower bound to mins
	P = sage_sparse_to_scipy_sparse(A, op=lambda x: RR(x.endpoints(rnd='RNDD')[0]))

	if P.min() < 0:
		raise NotImplementedError, "This method works only if A is nonnegative"

	n = P.shape[1]
	Ci = zeros(m)
	Ci[0] = 1.

	# Same bound as above, but ignoring the error ||A-P|| which is not relevant
	# Can be improved to an elementwise one, probably?

	prop = (gamma(max_nonzeros_per_row(P)) * matrix_norm1(P, rnd='RNDU')) * M
	zero = RealField(rnd='RNDU').zero()


	# We first compute r_l = min_i (M^l e_i), where the min is componentwise.

	mins = ones((n, m)) * inf
	for k in range(n):
		v = zeros(n)
		v[k] = 1.
		current_norm = 1.

		error = zero
		for i in range(1, m):
			v = P * v
			error += prop * current_norm
			current_norm = norm1(v, rnd='RNDU')

			w = v - error

			mask = w < mins[:, i]
			mins[mask, i] = w[mask]

	for i in range(1, m):
		Ci[i] = RealField(rnd='RNDU').one() - norm1(mins[:, i], rnd='RNDD')
	return Ci

def norm1_of_powers_alternate2(P, M, u, m=20, K=RealField(rnd='RNDU').zero()):
	"""
	Estimates :math:`\| A^k |V \|_1` for k in range(m), using an alternate algorithm

	:math:`A` is an operator such that :math:`AV \subseteq V`
	:math:`|V` means that the norm is restricted to the set of zero-sum vectors.

	Args:
		A (Sage interval matrix): a matrix :math:`ones((1,n))*A = ones((1,n))`
		u (numpy vector): an approximation of the invariant measure of A.
		m (integer): how many norms to compute
		M (real with RNDU): constant such that :math:`\|A^i\| \leq M` for each i in `range(m)`.
	
	Returns:
		a lists of norm bounds, starting from `k=0` (which is simply 1)

	Estimates differences wrt an approximation of v
    """

	n = P.shape[1]
	Ci = zeros(m)

	prop = (gamma(max_nonzeros_per_row(P)) * matrix_norm1(P, rnd='RNDU') + K) * M
	zero = RealField(rnd='RNDU').zero()

	Ci[0] = 1.
	for j in range(n):
		# We test the norm on these vectors: v_j = [1, 0, 0, ... , 0, -1 (in position j), 0, ... , 0]
		# The matrix that has them as columns satisfies ||Vx|| \geq ||x||.
		# So for each Vx \in span V, we have ||P^i V x|| / ||Vx|| \leq ||P^i V||*||x|| / ||x|| = ||P^i V||.

		v = zeros(n)
		v[j] = 1.

		# Error propagation: we use the bound
		#||(F^k-A^k)v|| \leq \sum ||A^(j-i-1) (F-A) F^i v|| \leq M ||F-A|| \sum ||v_i||,
		# where F is the operator that computes the product in floating point
		current_norm = 1.

		error = zero
		for i in range(1,m):
			v = P * v
			error += prop * current_norm
			current_norm = norm1(v, rnd='RNDU')
			Ci[i] = max(Ci[i], norm1(v - u, rnd='RNDU') + error)

#	for i in range(1, m): #if the bounds are larger than M, replace them with M
#		Ci[i] = min(Ci[i], M)

	return Ci

def norminf_of_powers(P, M, m=20, K=RealField(rnd='RNDU').zero()):
	"""
	Estimates :math:`\| A^k |V \|_Linf` for k in range(m).

	:math:`A` is an operator such that :math:`AV \subseteq V`
	:math:`|V` means that the norm is restricted to the set of zero-sum vectors.

	Args:
		P (scipy matrix): a matrix such that :math:`\|P-A\| \leq K`
		m (integer): how many norms to compute
		M, K (reals with RNDU): constants such that :math:`\|P-A\| \leq K` and :math:`\|A^i\| \leq M` for each i in `range(m)`.
	
	Returns:
		a lists of norm bounds, starting from `k=0` (which is simply 1)
    """

	n = P.shape[1]
	prop = (gamma(max_nonzeros_per_row(P)) * matrix_norminf(P, rnd='RNDU') + K) * M
	zero = RealField(rnd='RNDU').zero()

	maxes = zeros(m)
	for j in range(1, n):
		v = zeros(n)
		v[0] = 1.
		v[j] = -1.

		current_norm = 1.

		error = zero
		for i in range(1,m):
			v = P * v
			error += prop * current_norm
			current_norm = norminf(v, rnd='RNDU')
			maxes[i] = max(maxes[i], current_norm + error)

	# A generic extremal point in the Linf ball is the sum (with signs) of (n-1) elementary vectors 

	Ci = vector(RealField(rnd='RNDU'), maxes) * (n-1)
	Ci[0] = 1.

	for i in range(1, m): #if the bounds are larger than M, replace them with M
		Ci[i] = min(Ci[i], M)

	return Ci

def norminf_of_powers_alternate(P, M, m=20, K=RealField(rnd='RNDU').zero()):
	"""
	Estimates :math:`\| A^k |V \|_Linf` for k in range(m).

	:math:`A` is an operator such that :math:`AV \subseteq V`
	:math:`|V` means that the norm is restricted to the set of zero-sum vectors.

	Args:
		P (scipy matrix): a matrix such that :math:`\|P-A\| \leq K`
		m (integer): how many norms to compute
		M, K (reals with RNDU): constants such that :math:`\|P-A\| \leq K` and :math:`\|A^i\| \leq M` for each i in `range(m)`.
	
	Returns:
		a lists of norm bounds, starting from `k=0` (which is simply 1)
    """

	n = P.shape[1]
	prop = (gamma(max_nonzeros_per_row(P)) * matrix_norminf(P, rnd='RNDU') + K) * M
	zero = RealField(rnd='RNDU').zero()

	sums = zeros((n, m))
	for j in range(1, n):
		v = zeros(n)
		v[0] = 1.
		v[j] = -1.

		current_norm = 1.

		error = zero
		for i in range(1,m):
			v = P * v
			error += prop * current_norm
			current_norm = norminf(v, rnd='RNDU')
			sums[:,i] = sums[:,i] + v + error

	Ci = zeros(m)
	Ci[0] = 1.
	for i in range(1, m):
		# the term 1+gamma(2*i) keeps track of the errors obtained when summing
		# sums[:,i] + v + error (summation error on each component)
		# TODO: if proper rounding is used, this term can be omitted.
		# Not that it matters much, though.
		Ci[i] = norminf(sums[:,i]) * (1+gamma(2*i))

	for i in range(1, m): #if the bounds are larger than M, replace them with M
		Ci[i] = min(Ci[i], M)

	return Ci

class BasicTest(unittest.TestCase):
	"""
	Some tests.
	"""
	def test_small(self):
		import scipy.sparse as sp
		import numpy as np
		P = sp.csr_matrix([[0,1],[0,0]])
		expected = np.array([  1.00000000e+00,   1.00000000e+00,   6.66133815e-16, 6.66133815e-16,   6.66133815e-16])
		np.testing.assert_almost_equal(norm1_of_powers(P, 1, m=5), expected)

if __name__ == '__main__':
		unittest.main()
