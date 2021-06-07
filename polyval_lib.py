from sage.all import load
from time_profile import *



def polyval(polynomial,x):
	"""
	Evaluates a polynomial p (given by its coefficients, e.g., [0, 1, 2] represents :math:`x+2x^2`) in an interval point x.
	Args:
		p (collection):
		x (number):
	Returns:
		:math:`p(x)`, evaluated using the operations of x.
	"""
	
	result = x.parent()(0)
	for coefficient in reversed(polynomial):
		result = x*result + coefficient
	
	return result

load('polyval_lib.spyx')
def polyval_int(polynomial,x):
	return polyval_interval(polynomial,x)

