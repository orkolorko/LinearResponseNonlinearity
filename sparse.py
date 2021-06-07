"""
Utility functions on sparse matrices.
Mostly created to work around Sage inefficiencies.
"""

from __future__ import division
from sage.rings.real_mpfr import RealField
from sage.rings.real_mpfi import is_RealIntervalFieldElement
from sage.all import VectorSpace, RR, vector
import numpy as np
from scipy.sparse import csr_matrix
from collections import Counter
from itertools import izip
from sage.all import RealNumber,load,VectorSpace
from time_profile import *
load('cython_sparse_mat_vec.spyx')
from sage.matrix.matrix import Matrix

__all__ = ['sparse_matvec', 'sparse_vecmat', 'interval_norm1_error', 'sage_sparse_to_scipy_sparse', 'max_nonzeros_per_row', 'norm1']

def sparse_matvec(P,v):
	"""
	Compute Sage sparse matrix * vector product.
	
	Apparently, M*v fills the matrix in Sage, so this should be quicker
	"""
	
	w = VectorSpace(P.base_ring(),P.nrows())(0)
	for (ij, val) in P.dict().iteritems():
		w[ij[0]] += val*v[ij[1]]
	#cython_sparse_matvec(P,v,w)
	return w
	
def sparse_vecmat(v,P):
	"""
	Compute Sage sparse vector * matrix product.
	
	Apparently, v*M fills the matrix in Sage, so this should be quicker
	"""
	w = VectorSpace(P.base_ring(),P.ncols())(0)
	w = P.parent().row_space()(0)
	for (ij, val) in P.dict().iteritems():
		w[ij[1]] += v[ij[0]] * val
	#cython_sparse_vecmat(P,v,w)
	return w
	
	
def interval_norm1_error(P):
	"""
	Compute the "norm-1 width of the interval matrix P", i.e., an upper bound for :math:`\|P_1-P_2\|_1`, where the :math:`P_i` are matrices inside the interval :math:`P`.
	
	This is essentially the norm of P.upper()-P.lower().
	"""
	
	w = VectorSpace(RealField(rnd='RNDU'),P.dimensions()[1])(0)
	for (ij, val) in P.dict().iteritems():
		w[ij[1]] += val.absolute_diameter() #operations are rounded up because the LHS is so
	return max(w)

def interval_norminf_error(P):
	"""
	Compute the "norm-infinity width of the interval matrix P", i.e., an upper bound for :math:`\|P_1-P_2\|_inf`, where the :math:`P_i` are matrices inside the interval :math:`P`.
	
	This is essentially the norm of P.upper()-P.lower().
	"""
	
	w = VectorSpace(RealField(rnd='RNDU'),P.dimensions()[0])(0)
	for (ij, val) in P.dict().iteritems():
		w[ij[0]] += val.absolute_diameter() #operations are rounded up because the LHS is so
	return max(w)

def sage_sparse_to_scipy_sparse(P, op = lambda x: RR(x.center())):
	"""
	Convert a Sage sparse matrix to a scipy CSR one without going 
	through the full matrix.
	If the matrix is an interval one, converts the
	intervals to reals using op (default: center).
	Otherwise, op is ignored.
	
	Method taken from http://www.sagenb.org/src/plot/matrix_plot.py.
	"""
	
	if not is_RealIntervalFieldElement(P[0,0]):
		op = lambda x: x
	entries = list(P._dict().items())
	data = np.asarray([op(d) for _,d in entries], dtype=float)
	positions = np.asarray([[row for (row,col),_ in entries],
		[col for (row,col),_ in entries]], dtype=int)
	return csr_matrix((data, positions), shape=(P.nrows(), P.ncols()))


def max_nonzeros_per_row(P):
	"""
	Max number of nonzeros in a row of the matrix P
	"""
	if isinstance(P, Matrix): #Sage matrix
		c = Counter(i for (i,j) in P.dict())
		return max(c.values())
	else: #assuming scipy matrix
		c = Counter(P.nonzero()[0])
		return max(c.values())

def norm1(v, rnd='RNDN'):
	"""
	Compute the 1-norm of a vector, with configurable rounding mode.
	
	Returns a real with the correct rounding mode.
	"""
	
	nm = RealField(rnd=rnd)(0)
	for vi in v:
		nm += abs(vi) #operations are coherced to the ring of the LHS, so they are done with the correct rounding. Abs is exact in IEEE.
	return nm

def matrix_norm1(PP, rnd='RNDN'):
	"""
	Compute the 1-norm of a matrix, with configurable rounding mode
	
	Returns a real with the correct rounding mode.
	
	Args:
		PP (scipy sparse matrix):
	"""
	if not PP.format == 'csr':
		raise NotImplementedError
	column_norms = vector(RealField(rnd=rnd), PP.shape[1])
	for j, Pij in izip(PP.indices, PP.data):
		column_norms[j] += abs(Pij)
	return max(column_norms)

def norminf(v, rnd='RNDN'):
	"""
	Infinity-norm of a vector.

	The `rnd` parameter is there for consistency with norm1, but it is not used in computation.
	It affects the parent() field of the return type, though.
	"""

	return RealNumber(np.linalg.norm(v, np.inf), rnd=rnd)
	
def matrix_norminf(PP, rnd='RNDN'):
	"""
	Compute the infinity norm of a matrix, with configurable rounding mode.

	Args:
		PP (scipy sparse matrix):
	"""

	if not PP.format == 'csr':
		raise NotImplementedError

	return max(norm1(PP.data[PP.indptr[i]:PP.indptr[i+1]], rnd=rnd) for i in xrange(PP.shape[0]))
	

def gamma(n):
	"""
	gamma constants for error propagation; see [Higham, Accuracy and stability, p. 63]
	"""
	Rup = RealField(rnd='RNDU')
	Rdown = RealField(rnd='RNDD')
		
	gamma = Rup(np.finfo(float).eps) * n;
	if gamma >= 1:
		raise ValueError('n too large for error propagation bounds')

	gamma = gamma / (Rdown(1)-gamma)
	return gamma

def count_nnz_per_row(A):
	"""
	Computes the number of nonzeros per row, needs A to be in csr form
	"""
	return np.diff(A.indptr)