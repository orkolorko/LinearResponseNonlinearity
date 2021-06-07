"""
Generic functions to estimate decay times
"""

from sage.all import RealField,VectorSpace
from joblib import Parallel, delayed
from sparse import sage_sparse_to_scipy_sparse

def vector_decay_time(basis, PP, v, my_norm, M=0, K=0, target = 0.5):
	"""
	Generic function to estimate decay time on a single vector
	
	Args:
		PP (scipy matrix): approximation of the discretized matrix to use
		v (numpy vector):
		my_norm (function):
		M, K (reals with RNDU): constants such that :math:`\|P-PP\| \leq K` and :math:`\|P^n\| \leq M`, where :math:`P` is the exact discretized matrix.
	
	Returns:
		n such that :math:`\|PP^n v\| \leq target`.
		
	Raises:
		ValueError if insufficient precision is detected
	"""
	Rup = RealField(rnd='RNDU')
	
	current_norm = my_norm(v)
	error_on_computed_norm = M.parent().zero()
	n = 0
	MK = M * K
	error_propagation_constant = K #this constant is K at the first step and MK afterwards; see notes
	while current_norm + error_on_computed_norm >= target:
		n += 1
		v = basis.mat_vec(PP,v)
		current_norm = my_norm(v)
		error_on_computed_norm += error_propagation_constant * current_norm
		error_propagation_constant = MK
		#constants[n]=max(constants[n],(current_norm + error_on_computed_norm).upper())
		if error_on_computed_norm > target:
			raise ValueError, 'Insufficient precision'
	return n

def decay_time(dynamic,basis, P, alpha = 0.5, n_jobs = 1):
	r"""
	Number of iterations needed to contract all vectors in `basis.contracting_pairs()` to a given target alpha
	
	If (u_i,s_i) are the output of basis.contracting_pairs(),
	then we choose :math:`k` such that :math:`\|P^k u_i\| \leq \alpha s_i`
	so that the inequality :math:`\|P^k u\| \leq \sum |b_i| \|P^k u_i\| \leq \sum alpha |b_i| s_i \leq alpha \|u\|` holds for each vector :math:`u`.
	
	(The last equality is required in basis.contracting_pairs())
	"""
	
	Rup = RealField(rnd='RNDU')
	Rdown = RealField(rnd='RNDU')
	
	PP = sage_sparse_to_scipy_sparse(P)
	
	M = basis.bound_on_norms_of_powers(dynamic, project_left=True, project_right=True)
	K = basis.numerical_error(dynamic, P, PP)
	my_norm = lambda v: basis.norm_estimate(v)
	
	print 'bound on norms of powers, discretized operator', M
	print 'basis numerical error' , K
	
	alpha = Rdown(alpha) # we need a lower bound to alpha*s to make sure we contract below it
	decay_times = Parallel(n_jobs=n_jobs, verbose=1)(delayed(vector_decay_time)(basis,PP, v, my_norm, M, K, alpha*s ) for v, s in basis.contracting_pairs())
	
	return max(decay_times)

def estimate_tail(C, rho, n ,norm, tau):
	"""
	This function estimates the starting index of the tail such that
	the sum of the tail is less than tau.
	From ||L^{nk} g||_{s/w}\leq C\rho^k ||g||_s
		Args:
			C: the constant in front of rho
			rho: the estimated decay ratio
			n: the decay time
			norm: the norm of g
			tau: the desired error
	"""
	
	l = ((((1-rho)*tau)/(norm*n*C)).log()/rho.log()).ceil()
	return int((l*n).center())
	
	
