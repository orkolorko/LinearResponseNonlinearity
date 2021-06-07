"""
Generic functions to estimate decay times
"""

from sage.all import matrix,RealField,VectorSpace,MatrixSpace,RIF
from joblib import Parallel, delayed
from sparse import *
from sparse import matrix_norminf,interval_norminf_error
from time_profile import *

import numpy as np
import ctypes

from sage.all import RIF,load

load('decay_iterate_cython.spyx')

FE_TONEAREST = 0x0000
FE_DOWNWARD = 0x0400
FE_UPWARD = 0x0800
FE_TOWARDZERO = 0x0c00
libc = ctypes.CDLL('libm.so.6')



def decay_for_iterate_better_infinity(dynamic, basis, N , P, n_jobs = 1, output_rate=16, debug=False):
	"""
	This function bounds the norm in the following way
	let *v the componentwise absolute value of a vector
	Q[:,i]+=\sum *(P e_i), where the e_i's are a(ny) basis of the average 0 space
	This version is made using interval arithmetics in the matrix vector product
	"""
	Rup = RealField(rnd='RNDU')
	ring = dynamic.field
	V = VectorSpace(ring,len(basis))
	MS = MatrixSpace(ring,len(basis),1)
	
	count=0
	Norms = VectorSpace(Rup,N)(0)
	
	Q = matrix(Rup, len(basis), N)
	
	for v, s in basis.contracting_pairs(vectorspace=V):
		count=count+1
		for i in range(N):
			v=basis.interval_mat_vec(P,v,vectorspace=V)
			Q[:,i]+=MS(v.apply_map(abs))
		if (count%output_rate)==0: print count
			
	return Q

@timefunc
def decay_for_iterate_better_infinity_double(dynamic, basis, N , P, n_jobs = 1, block_size=1, output_rate=16, debug=False, interval_matrix=True, absolute_diameter=None):
	"""
	This function bounds the norm in the following way
	let *v the componentwise absolute value of a vector
	Q[:,i]+=\sum *(P e_i), where the e_i's are a(ny) basis of the average 0 space
	This version is made using interval arithmetics in the matrix vector product
	inputs:
	dynamic: 
	basis:
	N: the function bound the norms of P^0, ..., P^N
	n_jobs: number of parallel jobs (seems like it does not give a speedup)
	block_size: working on matrix instead than vectors (seems like it does not give a speedup)
	output_rate: 
	interval_matrix: if True, it expects a matrix of intervals
	absolute_diameter: the error on the coefficients of the matrix if the input is a double matrix
	
	Returns two matrices:
	Q: its j-th column contains \sum_i *P^j e_i 
	Norms: the row 0 contains the max(Q[:,i]), the row 1 contains the max(Q[:,i])+the numerical error, row 2 contains max(Q[:,i])+the numerical error+the error that comes from the 
	fact that we are using a double matrix instead than a interval matrix  
	"""
	
	if interval_matrix is True:
		PP=sage_sparse_to_scipy_sparse(P)
	else:
		PP=P
	
	Rup = RealField(rnd='RNDU')
	Rdown = RealField(rnd='RNDN')
	
	# we first bound from above the norm of the matrix P
	if interval_matrix is True:
		diameter = interval_norminf_error(P)
	else:
		diameter=Rup(absolute_diameter)
		if absolute_diameter is None:
			raise ValueError, "No diameter for the matrix provided"
	
	print "diameter ", diameter
	
	norm_inf_precompose=basis.norm_precompose()
	norm_inf_postcompose=basis.norm_postcompose()
	
	
	# find gamma
	if interval_matrix is True:
		gamma = Rup(np.finfo(float).eps) * max_nonzeros_per_row(P);
	else:
        nnz = max(np.diff(P.indptr))
		gamma = Rup(np.finfo(float).eps) * P.nnz
	gamma = gamma / (Rdown(1)-gamma)
	
	norm_of_the_matrix=matrix_norminf(PP, rnd='RNDU')
	print "Norm of the matrix ", norm_of_the_matrix
	norm_of_linear_application=norm_inf_postcompose*norm_of_the_matrix*norm_inf_precompose
	print "Norm of the linear application, with precompose and postcompose ", norm_of_linear_application
	print "Norm of the matrix times gamma", norm_of_linear_application*gamma
	diameter=norm_inf_postcompose*diameter*norm_inf_precompose
	print "diameter with precompose and postcompose ", diameter
	
	
	# initializes the matrix storing the sum of |Pv|, for v in contracting pairs
	Q = np.zeros([len(basis), N])
	
	#estimation vector by vector
	count=0
	if block_size==1:
		for v, s in basis.contracting_pairs():
			count=count+1
			for i in range(N):
				v=basis.mat_vec(PP,v)
				old_rnd=libc.fegetround()
				libc.fesetround(FE_UPWARD)
				Q[:,i]+=np.absolute(v)
				libc.fesetround(old_rnd)
			if (count%output_rate)==0: print count
	
	#blockwise estimation
	else:
		# the matrix V is a block matrix that contains block_size elements of the basis for the 
		# average zero measure space
		
		V=np.zeros(shape=(len(basis),block_size)) #matrix with blocksize columns
		h=(len(basis)-1)//block_size # the number of blocks we are going to contract
		r=(len(basis)-1)%block_size # the remaining vectors
		
		# we contract matrices of blocksize vectors
		for j in range(h):
			print "basis vectors from ", j*block_size, " to ",j*block_size+block_size-1
			for i in range(block_size):
				v,s=basis.contracting_pairs().next()
				V[:,i]=v
				for it in range(N):
					V=basis.mat_mat(PP,V)
					old_rnd=libc.fegetround()
					libc.fesetround(FE_UPWARD)
					Q[:,it]+=(np.absolute(V)).sum(axis=1)
					libc.fesetround(old_rnd)
		V=np.zeros(shape=(len(basis),block_size))
		for j in range(r):
			print "the last ", r, " vectors"
			v,s=basis.contracting_pairs().next()
			V[:,j]=v
			for it in range(N):
				V=basis.mat_mat(PP,V)
				old_rnd=libc.fegetround()
				libc.fesetround(FE_UPWARD) 
				Q[:,it]+=(np.absolute(V)).sum(axis=1)
				libc.fesetround(old_rnd)

	# we estimate now the infinity norm, componentwise
	
	old_rnd=libc.fegetround()
	libc.fesetround(FE_UPWARD)
	Norms=np.zeros([3,N+1])
	# the norm matrix, the first column contains
	# row 0, the infinity norm of the starting vectors
	# row 1, the infinity norm of the numerical matrix
	# row 2, the infinity norm of P
	Norms[0,0]=1
	Norms[1,0]=1
	Norms[2,0]=1	
	
	# the infinity norm of the 1-st iteration
	Norms[0,1]=max(Q[:,0])
	# keeping track of the numerical error
	Norms[1,1]=Norms[0,1]+norm_of_the_matrix*gamma
	#keeping track of the error that comes from the fact that we passed from interval to double
	Norms[2,1]=Norms[1,1]+diameter
	
	for i in range(2,N+1):
		Norms[0,i]=max(Q[:,i-1])
		add_num=Rup(0)
		add_int=Rup(0)
		for j in range(i):
			add_num+=Norms[1,j]*norm_of_the_matrix*gamma*Norms[0,i-1-j]
			add_int+=Norms[1,j]*diameter*Norms[2,i-1-j]
		Norms[1,i]=Norms[0,i]+add_num
		Norms[2,i]=Norms[1,i]+add_int
	libc.fesetround(old_rnd)
	
	return Q,Norms

