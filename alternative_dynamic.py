from dynamic import *
from polyval_lib import polyval
from interval import func_range,func_range_with_der

from vandermonde import *


def C_3_bump_coeffs():
	M=Vandermonde_first_second_third(0,1,1,1)
	
	phi_0=SR(1)
	phi_1=SR(0)
		
	phi_prime_0=SR(0)
	phi_prime_1=SR(0)
	
	phi_second_0=SR(0)
	phi_second_1=SR(0)
	
	phi_third_0=SR(0)
	phi_third_1=SR(0)
		
	v=M.inverse()*vector(SR,[phi_0,phi_1,phi_prime_0,phi_prime_1,phi_second_0,phi_second_1,phi_third_0,phi_third_1])
	
	return v

def C_3_bump(x):
	"""
	Evaluates an element phi
	Args:
	x (real number) : the point on which we evaluate the function 
	"""
	tilde_x=x
	Interval_Field=x.parent()
	try:
		tilde_x=tilde_x.intersection(Interval_Field(-1,1))
	except ValueError:
		return Interval_Field(0)
	try:
		tilde_x_neg=tilde_x.intersection(Interval_Field(-1,0))
		val_neg = polyval([1,0,0,0,-35,-84,-70,-20], tilde_x_neg)
	except ValueError:
		val_neg=None
	try:
		tilde_x_pos=tilde_x.intersection(Interval_Field(0,1))
		val_pos = polyval([1,0,0,0,-35,84,-70,20], tilde_x_pos)
	except ValueError:
		val_pos=None
	if val_pos is None:
		return val_neg
	if val_neg is None:
		return val_pos
	return val_neg.union(val_pos)

def C_3_bump_prime(x):
	"""
	Evaluates an element phi
	Args:
	x (real number) : the point on which we evaluate the function 
	"""
	tilde_x=x
	Interval_Field=x.parent()
	try:
		tilde_x=tilde_x.intersection(Interval_Field(-1,1))
	except ValueError:
		return Interval_Field(0)
	try:
		tilde_x_neg=tilde_x.intersection(Interval_Field(-1,0))
		val_neg = polyval([0,0,0,-35*4,-84*5,-70*6,-20*7], tilde_x_neg)
	except ValueError:
		val_neg=None
	try:
		tilde_x_pos=tilde_x.intersection(Interval_Field(0,1))
		val_pos = polyval([0,0,0,-35*4,84*5,-70*6,20*7], tilde_x_pos)
	except ValueError:
		val_pos=None
	if val_pos is None:
		return val_neg
	if val_neg is None:
		return val_pos
	return val_neg.union(val_pos)

def C_3_bump_second(x):
	"""
	Evaluates an element phi
	Args:
	x (real number) : the point on which we evaluate the function 
	"""
	tilde_x=x
	Interval_Field=x.parent()
	try:
		tilde_x=tilde_x.intersection(Interval_Field(-1,1))
	except ValueError:
		return Interval_Field(0)
	try:
		tilde_x_neg=tilde_x.intersection(Interval_Field(-1,0))
		val_neg = polyval([0,0,-35*4*3,-84*5*4,-70*6*5,-20*7*6], tilde_x_neg)
	except ValueError:
		val_neg=None
	try:
		tilde_x_pos=tilde_x.intersection(Interval_Field(0,1))
		val_pos = polyval([0,0,-35*4*3,84*5*4,-70*6*5,20*7*6], tilde_x_pos)
	except ValueError:
		val_pos=None
	if val_pos is None:
		return val_neg
	if val_neg is None:
		return val_pos
	return val_neg.union(val_pos)

def C_3_bump_third(x):
	"""
	Evaluates an element phi
	Args:
	x (real number) : the point on which we evaluate the function 
	"""
	tilde_x=x
	Interval_Field=x.parent()
	try:
		tilde_x=tilde_x.intersection(Interval_Field(-1,1))
	except ValueError:
		return Interval_Field(0)
	try:
		tilde_x_neg=tilde_x.intersection(Interval_Field(-1,0))
		val_neg = polyval([0,-35*4*3*2,-84*5*4*3,-70*6*5*4,-20*7*6*5], tilde_x_neg)
	except ValueError:
		val_neg=None
	try:
		tilde_x_pos=tilde_x.intersection(Interval_Field(0,1))
		val_pos = polyval([0,-35*4*3*2,84*5*4*3,-70*6*5*4,20*7*6*5], tilde_x_pos)
	except ValueError:
		val_pos=None
	if val_pos is None:
		return val_neg
	if val_neg is None:
		return val_pos
	return val_neg.union(val_pos)


def smooth_bump(x):
	"""
	A smooth bump function phi
	"""
	IntervalField=x.parent()
	try:
		tilde_x=x.intersection(IntervalField(-1,1))
		y=IntervalField(-1)/(IntervalField(1)-tilde_x**2)
	except ValueError:
		return IntervalField(0)
	return y.exp()

#problem, this is not numerically well-behaved
def smooth_bump_prime(x):
	"""
	A smooth bump function phi
	"""
	IntervalField=x.parent()
	try:
		tilde_x=x.intersection(IntervalField(-1,1))
		y=IntervalField(-1)/(IntervalField(1)-tilde_x**2)
		y_prime=IntervalField(-2)*tilde_x/((IntervalField(1)-tilde_x**2)**2)
	except ValueError:
		return IntervalField(0)
	return y.exp()*y_prime

def bump(x):
	return C_3_bump(x)

def bump_prime(x):
	return C_3_bump_prime(x)

def bump_second(x):
	return C_3_bump_second(x)

def bump_third(x):
	return C_3_bump_third(x)

class BumpPerturbedLinearDynamic(Mod1Dynamic):
	"""
	Dynamic :math:`k x+c\phi(m*(x-a)) \mod 1`, for tunable parameters ``c`` and ``j``.
	"""
	def __init__(self, prec=53, k=4 , c='0.5', m=4,a='0.5'):
		"""
		Args:
			
		"""
		Mod1Dynamic.__init__(self, k, prec)
		self._k = self.field(k)
		self._c = self.field(c)
		self._a=self.field(a)
		self._m=m
		
		self._der_range=func_range(self.fprime_aux,x=self.field(0,1),epsilon=0.001)
		infTprime = (self._der_range).lower()
		
		self._distorsion=func_range(lambda x : (self.fsecond_aux(x)/(self.fprime_aux(x)**2)).abs(),x=self.field(0,1),epsilon=0.001)
		
		if not infTprime > 1:
			warn("The dynamic is not expanding with this choice of the parameters")
	
	def f_unquotiented(self, x):
		"""
		"""
		warn_if_not_interval(x)
		warn_if_not_in_01(x)
		y=self._k * x+self._c * bump(self._m*(x-self._a))
		return y 
		
	def fprime_aux(self,x):
		warn_if_not_interval(x)
		warn_if_not_in_01(x)
		return self._k + self._c * self._m * bump_prime(self._m*(x-self._a))
		
	def fprime(self, x):
		warn_if_not_interval(x)
		warn_if_not_in_01(x)
		return self.fprime_aux(x).intersection(self._der_range)
	
	def fsecond_aux(self, x):
		warn_if_not_interval(x)
		warn_if_not_in_01(x)
		return self._k + self._c * (self._m**2) * bump_second(self._m*(x-self._a))
	
	def fsecond(self,x):
		warn_if_not_interval(x)
		warn_if_not_in_01(x)
		return func_range(self.fsecond_aux,x,epsilon=0.001)


	def expansivity(self):
		return (self._der_range).lower()

	def distorsion(self):
		return self._c * (self._jpi**2) / self.expansivity()**2
	
	def third_distorsion(self):
		return self._c * self._jpi * self._jpi* self._jpi / self.expansivity()**3

	def is_continuous_on_torus(self):
		return True
	
	def fthird(self, x):
		warn_if_not_interval(x)
		return -self._c * (self._jpi**3) * (self._jpi * x).cos()

class TwoBumpPerturbedLinearDynamic(Mod1Dynamic):
	"""
	Dynamic :math:`k x+c\phi(m*(x-a)) \mod 1`, for tunable parameters ``c`` and ``j``.
	"""
	def __init__(self, prec=53, k=4 , c='0.5', m=4,a_1='0.25',a_2='0.75'):
		"""
		Args:
			
		"""
		Mod1Dynamic.__init__(self, k, prec)
		self._k = self.field(k)
		self._c = self.field(c)
		self._a_1=self.field(a_1)
		self._a_2=self.field(a_2)
		self._m=m
		
		self._der_range=func_range(self.fprime_aux,x=self.field(0,1),epsilon=0.001)
		infTprime = (self._der_range).lower()
		
		self._distorsion=func_range(lambda x : (self.fsecond_aux(x)/(self.fprime_aux(x)**2)).abs(),x=self.field(0,1),epsilon=0.001)
		self._third_distorsion=func_range(lambda x : (self.fthird_aux(x)/(self.fprime_aux(x)**3)).abs(),x=self.field(0,1),epsilon=0.001)
		
		
		if not infTprime > 1:
			warn("The dynamic is not expanding with this choice of the parameters")
	
	def f_unquotiented(self, x):
		"""
		"""
		warn_if_not_interval(x)
		warn_if_not_in_01(x)
		y=self._k * x+self._c * bump(self._m*(x-self._a_1))+self._c * bump(self._m*(x-self._a_2))
		return y 
		
	def fprime_aux(self,x):
		warn_if_not_interval(x)
		warn_if_not_in_01(x)
		return self._k + self._c * self._m * bump_prime(self._m*(x-self._a_1))+ self._c * self._m * bump_prime(self._m*(x-self._a_2))
		
	def fprime(self, x):
		warn_if_not_interval(x)
		warn_if_not_in_01(x)
		return self.fprime_aux(x).intersection(self._der_range)
	
	def fsecond_aux(self, x):
		warn_if_not_interval(x)
		warn_if_not_in_01(x)
		return self._c * (self._m**2) * bump_second(self._m*(x-self._a_1))+ self._c * (self._m**2) * bump_second(self._m*(x-self._a_2))
	
	def fsecond(self,x):
		warn_if_not_interval(x)
		warn_if_not_in_01(x)
		return func_range(self.fsecond_aux,x,epsilon=0.001)


	def expansivity(self):
		return (self._der_range).lower()

	def distorsion(self):
		return (self._distorsion).upper()
	
	def third_distorsion(self):
		return (self._third_distorsion).upper()

	def is_continuous_on_torus(self):
		return True
	
	def fthird_aux(self, x):
		warn_if_not_interval(x)
		return self._c * (self._m**3) * (bump_third(self._m*(x-self._a_1))+ bump_third(self._m*(x-self._a_2)))

class MorePerturbedLinearDynamic(Mod1Dynamic):
	"""
	Dynamic :math:`k x+c*(cos(2kpi x)+0.25 cos(4kpi x)) \mod 1`, for tunable parameter ``c``.
	"""
	def __init__(self, prec=53, k=4, c='0.001'):
		"""
		Args:
			
		"""
		Mod1Dynamic.__init__(self, k, prec)
		self._k = self.field(k)
		self._c = self.field(c)
		self._two_k_pi = self.field.pi() *self.field(2)*self._k
		
		self._der_range=func_range(self.fprime_aux,x=self.field(0,1),epsilon=0.001)
		infTprime = (self._der_range).lower()
		
		self._distorsion=func_range(lambda x : (self.fsecond_aux(x)/(self.fprime_aux(x)**2)).abs(),x=self.field(0,1),epsilon=0.001)
		self._third_distorsion=func_range(lambda x : (self.fthird_aux(x)/(self.fprime_aux(x)**3)).abs(),x=self.field(0,1),epsilon=0.001)
		
		
		if not infTprime > 1:
			warn("The dynamic is not expanding with this choice of the parameters")
	
	def f_unquotiented(self, x):
		"""
		"""
		warn_if_not_interval(x)
		warn_if_not_in_01(x)
		y=self._k * x+self._c * ((self._two_k_pi*x).sin()+self.field(0.25)*(self.field(2)*self._two_k_pi*x).sin())
		return y 
		
	def fprime_aux(self,x):
		warn_if_not_interval(x)
		warn_if_not_in_01(x)
		return self._k + self._c * self._two_k_pi*((self._two_k_pi*x).cos()+self.field(0.5)*(self.field(2)*self._two_k_pi*x).cos())
		
	def fprime(self, x):
		warn_if_not_interval(x)
		warn_if_not_in_01(x)
		return self.fprime_aux(x).intersection(self._der_range)
	
	def fsecond_aux(self, x):
		warn_if_not_interval(x)
		warn_if_not_in_01(x)
		return - self._c * (self._two_k_pi**2)*((self._two_k_pi*x).sin()+(self.field(2)*self._two_k_pi*x).sin())
	
	def fsecond(self,x):
		warn_if_not_interval(x)
		warn_if_not_in_01(x)
		return func_range(self.fsecond_aux,x,epsilon=0.001)


	def expansivity(self):
		return self.field((self._der_range).lower())

	def distorsion(self):
		return self.field((self._distorsion).upper())
	
	def third_distorsion(self):
		return self.field((self._third_distorsion).upper())

	def is_continuous_on_torus(self):
		return True
	
	def fthird_aux(self, x):
		warn_if_not_interval(x)
		return - self._c * (self._two_k_pi**3)*((self._two_k_pi*x).cos()+self.field(2)*(self.field(2)*self._two_k_pi*x).cos())

