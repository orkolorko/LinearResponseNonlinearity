"""
In this file we put some functions specific for the study of the linear response
of the LSV family at the boundary
"""

import numpy as np
from sage.all import ln,sqrt,exp,RR,pi
from ulam import *
from sparse import *

def compute_f_eta_lsv_at_boundary(b_Ulam):
	"""
	This function gives us the approximation of \hat{L}h on the Ulam basis
	Input:
	b_Ulam: the Ulam basis
	"""
	size=len(b_Ulam)
	eta=1.0/size
	w=np.zeros(size) 
	w[0]=-eta*ln(eta)/4.0
	for i in range(1,size):
		w[i]=-((i+1)*eta*ln((i+1)*eta)-i*eta*ln(i*eta))/4.0
	return w*size

def bound_size_tail_lsv_at_boundary(n):
	"""
	This function bounds the size of the tail of \sum_{n}^{+\infty}L^i \hat{L}h
	for the lsv family at the boundary
	Input:
	n: the first index of the tail
	"""
	return 1.0/2**(n+1)*(1+ln(sqrt(2.0*RR(pi))))+(n+1)/2**(n+1)*ln(2.0)/2+1.0/9*1/2**(2*n+2)

def compute_variation_f_eta(b_Ulam):
	"""
	This function bounds the variation of f_eta
	Input:
	w: the vector containing the coefficients of f_eta
	"""
	size=len(b_Ulam)
	eta=1.0/size
	return (-ln(eta)-(size-1)*ln((size-1)*eta))/4.0

def compute_error_ulam_approximation_f_eta(b_Ulam):
	"""
	This function bounds the L1 error due to the discretization of \hat{L}h
	Input:
	b_Ulam:
	"""
	size=len(b_Ulam)
	eta=1.0/size
	return eta/(2*exp(1.0))-eta*ln(eta)/4
	
def compute_error_sum_L_i_f_eta_minus_hat_L_h(b_Ulam,l):
	"""
	This function computes the size of ||\sum_{i=0}^{l-1}L^i(f_eta-\hat{L}h)||_{L^i}
	Input:
	b_Ulam:
	l: the first index of the tail
	"""
	return l*compute_error_ulam_approximation_f_eta(b_Ulam)
	
def compute_sum_L_i_f_eta_minus_L_eta_i_f_eta(b_Ulam,l,w):
	"""
	This function computes the size of ||\sum_{i=0}^{l-1}L^i f_eta-L^i_{\eta}f_{\eta}||_{L^i}
	Input:
	b_Ulam:
	l: the first index of the tail
	"""
	size=len(b_Ulam)
	eta=1.0/size
	
	return 3*eta/2*(2*l-((2**l-1)/2**(l-1)))*compute_variation_f_eta(b_Ulam)

def compute_linear_response(b_Ulam,P,f_eta,l, interval_matrix=True):
	"""
	This function computes the approximation of the linear response
	Inputs:
	b_Ulam:
	PP:
	f_eta:
	l:
	"""
	w=f_eta.copy()
	
	if interval_matrix:
		PP=sage_sparse_to_scipy_sparse(P)
	else:
		PP=P
	
	for i in range(l-1):
		w=f_eta+b_Ulam.mat_vec(PP,w)
	
	return w
	
	
