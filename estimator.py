"""
Error estimates
"""

from __future__ import division
from sage.all import ceil, log, vector, RealField
from partition import partition_diameter
from refiner import sum_on_same_length

def global_error(D, partition, N, alpha, residual, type='L1'):
	r"""
	Estimate the global error using decay time and Perron vector residual.
	
	Args:
		N (int), alpha (real): so that :math:`\|L^N\| \leq \alpha`.
		residual (interval): the norm-1 residual :math:`|Pv-v|` of the computed approximation of the eigenvector
		partition (int or sequence):
		type (string): 'L1' or 'Linf'/'C0'. The norms in the other inputs above have match with this.
	"""
	
	delta = partition_diameter(partition)
	l, B = D.lasota_yorke_constants(type)
	M = D.bound_on_norms_of_powers(type)
	intervalfield = B.parent()
	
	# (1-alpha) generalizes the factor 2 in [Galatolo, Nisoli, Theorem 1].
	# The rest is just their Theorem 1.
	amplification_factor = intervalfield(N) * M / (intervalfield(1)-alpha)
	
	if not l < 1:
		raise ValueError, 'The Lasota-Yorke inequality does not have lambda<1, something unexpected'
	if not alpha < 1:
		raise ValueError, 'alpha has to be smaller than 1 (true contraction)'
	
	# note that all the arithmetic operations done here have to propagate to intervals correctly
	# i.e., we can't write stuff like float * float + interval
	if type == 'L1':
		# 2* delta* varBound bounds ||Lf-L_delta f||_1 [Galatolo, Nisoli, Lemma 10].
		varBound = B / (1-l) 
		approximationError = varBound * 2 * delta
	elif type == 'Linf' or type == 'C0':
		weakNormBound = B + 1 # last remark in sec 6.1 of GalN
		lipBound = B / (1-l) * weakNormBound
		# Let p be the projection; we have ||Lf-pLpf|| \leq ||(id-p)Lf||+||pL(id-p)f|| \leq (1+M)*Lip f*delta, if ||p||
		approximationError = lipBound * (M+1) * delta # the 2M is here because we have two "projections"
	else:
		raise ValueError, 'Unknown norm type'
	
	# this is the difference between the exact invariant measure and the exact eigenvector of the Ulam discretization
	error1 = approximationError * amplification_factor
	
# this is the difference between the exact eigenvector of the Ulam discretization and the vector we computed.
# proof for this formula is on our Dropbox notes.
	error2 = residual * amplification_factor
	print "Detailed errors: discretization %g, eigenvector residual %g" % (error1.upper(), error2.upper())
	return (error1 + error2).magnitude()

def twogrid_decay_time(D, coarsePartition, finePartition, coarseN, coarseAlpha):
	"""
	Estimate the decay time on the finer grid finePartition using decay time on coarsePartition.
	
	coarsePartition, finePartition may be either numbers (equal intervals) or grid arrays on [0,1].
	N, alpha are such that ||LC^(coarseN)g||_1 <= coarseAlpha*||g||_1 for each g.
	
	Returns (fineN, fineAlpha) so that ||LF^(fineN)||_1 <= fineAlpha. fineAlpha will be 1/2 in this implementation, it is returned only for shape compatibility.
	
	Based on the formula ||LF^(coarseN+K)g|| = (coarseAlpha + deltaC*(coarseN+ll)*BB)*||g||_1 + ...
		... + deltaC*lambda^K*ll*var(g), see notes
	
	BB = B/(1-lambda), ll = lambda/(1-lambda)
	
	"""

	l, B = D.lasota_yorke_constants('L1')
	
	deltaC = partition_diameter(coarsePartition)
	deltaF = partition_diameter(finePartition)
	
	ll = 1 / (1-l)
	BB = B / (1-l)
	
	C1 = coarseAlpha + BB * deltaC*(coarseN+1+ll)
	
	if C1.magnitude() >= 0.5:
		raise ValueError, "Cannot go on, the decay time estimate is insufficient. Retry with a finer CoarsePartition or a smaller coarseAlpha."
	
	# let v_i be the vectors in the fine partition with one component 0.5, one -0.5 and the rest zero.
	# ||v_i||_1 =  deltaF and Var v_i = 2, so for those vectors var(g) = 2/deltaF * ||g||_1.
	# Note that the norm-1 maximum is attained on one of those (without the additional 1/2 factor).
	
	# Hence we need C1 + C2*lambda^K <= 1/2, where C2 = deltaC*ll * 2/deltaF
	
	C2 = deltaC*ll*2/deltaF

	K = log((0.5 - C1)/C2) / log(l)
	K = ceil(K.magnitude())
	return (coarseN + K, 1/2)

def minimize_on_integers(f, xmin, xmax):
	"""
	Find argmin(f) in range(xmin,xmax).
	
	Brutal implementation by trial and error.
	"""
	minvalue = float('inf')
	for x in range(xmin,xmax):
		fx = f(x)
		if fx < minvalue:
			argmin = x
			minvalue = fx
	
	return argmin

def optimal_twogrid_decay_time(D, coarsePartition, finePartition, coarseN, coarseAlpha):
	"""
	Estimate the decay time on the finer grid FinePartition using decay time on CoarsePartition.
	
	As twogrid_decay_time, but returns an "optimal" value.
	"""

	l, B = D.lasota_yorke_constants('L1')

	ll = 1 / (1-l)
	BB = B / (1-l)

	deltaC = partition_diameter(coarsePartition)
	deltaF = partition_diameter(finePartition)
	
	C1 = coarseAlpha + BB * deltaC*(coarseN+1+ll)
	C2 = deltaC*ll*2/deltaF
	
	if C1.magnitude() >= 1:
		raise ValueError, "Cannot go on, the decay time estimate is insufficient. Retry with a finer CoarsePartition or a smaller coarseAlpha."
	
	#the minimization needs not be rigorous, so it's not.
	objective_function = lambda(K): (coarseN+K)/(1-C1.center()-C2.center()*(l.center()**K))
	minK = ceil(log((1-C1.center())/C2.center()) / log(l.center()))
	minK = max(minK,1) #ensures K is never smaller than 1
	K = minimize_on_integers(objective_function,minK,max(100,5*minK)) #magic "reasonable" number for the maximum
	# ok, K is the (non-rigorous) optimal value
	
	fineAlpha = C1 + C2*(l**K)
	assert fineAlpha < 1
	return (coarseN + K, fineAlpha.upper())
	
	
def localized_error(D, v, partition, epsilon, newp, N, alpha, residual):
	"""
	Given an approximation ``v`` of the invariant measure of ``D`` on a partition ``partition``,
	accurate within a factor ``epsilon``, construct a localized error estimate on the unequal
	partition ``newp``.
	"""
	
	l, B = D.lasota_yorke_constants('L1')
	BB = B / (1-l)
	intervalfield = B.parent()
	amplification_factor = intervalfield(N) / (intervalfield(1)-alpha)
	
	v = vector(v).change_ring(RealField(rnd='RNDU')) #so that operations with v are rounded up. There seems to be a Sage bug with the syntax vector(RealField(rnd='RNDU'),v), so we use this one.
	assert v.parent().base().rounding_mode() == 'RNDU'
	d = sum_on_same_length(v, partition, newp)
	
	# this bounds the distance ||L_delta*f-f||
	projBound = 2*BB*sum(delta*(mJ+epsilon) for (delta, mJ) in d.iteritems())
	
	error1 = projBound * amplification_factor
	error2 = residual * amplification_factor

	print "Errors: discretization %g, eigenvector residual %g" % (error1.upper(), error2.upper())
	return (error1 + error2).magnitude()

