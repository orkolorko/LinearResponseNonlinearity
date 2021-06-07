"""
This program is used to get the elementary convergence estimate from 

"""

from sage.rings.real_mpfr import RealField
from sage.all import RealNumber, NaN
from sage.rings.real_mpfi import RealIntervalField
from sage.matrix.matrix_space import MatrixSpace
from joblib import Parallel, delayed
import sys

def rho(lambda_1,lambda_2,delta,n_1,A,B,C,D):
	RI = RealInterval
	mu = RI(A)*RI(lambda_1)**n_1+delta*n_1*D+RI(lambda_2)
	nu = RI(A)*RI(lambda_1)**n_1-delta*n_1*D-RI(lambda_2)
	
	tau = 4*delta*RI(B)*C
	discriminant = sqrt(nu**2+tau)
	dom_eig = (mu+discriminant)/2
	small_eig = (mu-discriminant)/2
	
	a = 1/(2*RI(B))*(RI(A)*RI(lambda_1)**n_1-delta*n_1*D-RI(lambda_2)+discriminant)
	b = RI(1)

	normalize=a+b
	
	a/=normalize
	b/=normalize

	return [dom_eig,small_eig,a,b]

def compute_C_D(lambda_1,lambda_2,delta,n_1,A,B):
	C=(RI(A)*RI(lambda_1)+1)*RI(A)/(1-RI(lambda_1))
	D=RI(B)*(RI(A)*RI(lambda_1)+2)
	return [C,D]

prec=53
RI=RealInterval
A = input('Input A: ')
lambda_1 = input ('Input lambda_1: ')
lambda_2 = input ('Input lambda_2: ')
k = input('Input k: ')
delta = RI(1)/k
n_1 = input('Input n_1: ')
B = input('Input B: ')

[C,D]=compute_C_D(lambda_1,lambda_2,delta,n_1,A,B)

M = MatrixSpace(RealIntervalField(prec),2,2)
V = VectorSpace(RealIntervalField(prec),2)

P = M(0)

P[0,0]=RI(A)*RI(lambda_1)**n_1
P[0,1]=B
P[1,0]=delta*C
P[1,1]=delta*n_1*D+lambda_2

[dom_eig,small_eig,a,b] = rho(lambda_1,lambda_2,delta,n_1,A,B,C,D)

v = V(0)
v[0] = a
v[1] = b

PP=P.transpose()

print(dom_eig.endpoints())
print(a.endpoints())
print(b.endpoints())
s= a.str(style = 'brackets',error_digits=3) + '||g||_{BV}+' + b.str(style = 'brackets',error_digits=3) + '||g||_1'
print(s)

Trunc=RealField(15,True,'RNDU')

for i in range(1,5):
	j=2*i
	PP=P**j
	s = 'h=' + repr(j*n_1) + ' & ' + (Trunc((PP[0,0]).upper())).str() + '||g||_{BV}+' + repr(Trunc((PP[0,1]).upper())) + '||g||_1' + ' & ' + repr(Trunc((PP[1,0]).upper())) + '||g||_{BV}+' + repr(Trunc((PP[1,1]).upper())) + '||g||_1'  
	print(s)

coeff_1=A/v[0]+B/v[1]
coeff_2=B/v[1]

