"""
In this file is contained the function that
permits to estimate the speed of decay of correlations
using the small matrix method
"""

from sage.all import sqrt


def rho(a,b,c,d):	
	return 0.5*(a+d+sqrt((a-d)**2+4*b*c))

def autonorm(a,b,c,d):
	v_1=a-d+sqrt((a-d)**2+4*b*c)
	v_2=2*b
	lam=v_1+v_2
	
	return [v_1/lam,v_2/lam]


def spectral_radius_and_autonorm_uniform(A, lam_1, B, n, lam_2, C, D):
	""" This function estimates the spectral radius $\rho$ and the coefficients
	(a,b) such that, if :math:'||g||_{a,b}=a||g||_s+b||g||_w'
	Galatolo-Nisoli-Saussol, we have that
	:math:'||L^n_{delta} g||_{a,b}\leq \rho ||g||_{a,b}', for any g of zero average
	
	Please note that this function differs from the one above since it works for the 
	discretized operator.
	
	Args:
	A: coeff L-Y
	lam: coeff L-Y
	B: coeff L-Y
	n: the number of iterates
	lam_2: the norm of L^n_{\delta}
	C: coefficient of the strong norm (multiplied by eta)
	D: coefficient of the weak norm (multiplied by eta)
	"""

	# building the matrix
	
	a=A*lam_1**n
	b=B
	c=C
	d=D+lam_2
	
	print "[",a, b,"]\n[", c, d,"]"
	
	return [rho(a,b,c,d),autonorm(a,b,c,d)]
