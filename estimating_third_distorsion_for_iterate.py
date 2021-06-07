from dynamic_wael import *
from sage.all import RIF, VectorSpace

def orbit_for_iterate(D,n_iter,x):
	V=VectorSpace(D.field,n_iter+1)(0)
	
	V[0]=x
	
	f=D.f(x)
	f_prime=D.fprime(x)
	
	for i in range(1,n_iter+1):
		V[i]=D.f(V[i-1])
	
	return V

def first_derivative_for_iterate(D,n_iter,x):
	W=VectorSpace(D.field,n_iter+1)(0)
	
	V=orbit_for_iterate(D,n_iter,x) # v contains x_i=T^i(x)
	
	f_prime=D.fprime(x)
	
	W[0]=1
	
	for i in range(1,n_iter+1):
		W[i]=D.fprime(V[i-1])*W[i-1] #T'(T^{i-1}(x))*(T^{i-1})'(x)
	
	return W

def second_derivative_for_iterate(D,n_iter,x):
	Z=VectorSpace(D.field,n_iter+1)(0)
	
	V=orbit_for_iterate(D,n_iter,x) # v contains x_i=T^i(x)
	W=first_derivative_for_iterate(D,n_iter,RIF(0,0.1))
	
	for i in range(1,n_iter+1):
		Z[i]=D.fsecond(V[i-1])*W[i-1]**2+D.fprime(V[i-1])*Z[i-1] #T''(T^{i-1}(x))*[(T^{i-1})'(x)]^2+T'(T^{i-1}(x))*T^{i-1}''(x)
	
	return Z

def third_derivative_for_iterate(D,n_iter,x):
	Y=VectorSpace(D.field,n_iter+1)(0)
	
	V=orbit_for_iterate(D,n_iter,x) # v contains x_i=T^i(x)
	W=first_derivative_for_iterate(D,n_iter,x)
	Z=second_derivative_for_iterate(D,n_iter,x)
	
	for i in range(1,n_iter+1):
		Y[i]=D.fthird(V[i-1])*W[i-1]**3+3*D.fsecond(V[i-1])*W[i-1]*Z[i-1]+D.fprime(V[i-1])*Y[i-1] 
		#T'''(T^{i-1}(x))*[(T^{i-1})'(x)]^3+3T''(T^{i-1}(x))*T^{i-1}'(x)*T^{i-1}''(x)+T'(T^{i-1}(x))*T^{i-1}'''(x)
	return Y

def estimate_distorsion(D,n_iter,range_prec=0.0001):
	f=lambda x : second_derivative_for_iterate(D,n_iter,x)[n_iter]
	g=lambda x : first_derivative_for_iterate(D,n_iter,x)[n_iter]
	h=lambda x : f(x)/(g(x)**2)
	return func_range(h,D.field(0,1),range_prec)

def third_distorsion(D,n_iter,range_prec=0.0001):
	f=lambda x : third_derivative_for_iterate(D,n_iter,x)[n_iter]
	g=lambda x : first_derivative_for_iterate(D,n_iter,x)[n_iter]
	h=lambda x : f(x)/(g(x)**3)
	plot(h,0,1)
	return func_range(h,D.field(0,1),range_prec)


def test_iter(n,range_prec=0.0001):
	D=WaelDynamic_parametric(a=2, k=2,prec=53)
	return estimate_distorsion(D,n,range_prec=0.0001)
