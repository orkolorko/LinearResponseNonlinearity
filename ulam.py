"""
Ulam basis written in the generic assembler framework
"""

from __future__ import division
from basis import Basis
from partition import check_partition, partition_diameter, partition_minimum_diameter, is_refinement, is_inside, overlap, nonzero_on

from sage.all import load, RealNumber, RIF, floor, VectorSpace

import numpy as np
from spectral_radius import *

from sage.all import vector, RealField, RealIntervalField
from sparse import max_nonzeros_per_row, norm1, matrix_norm1, sparse_matvec
from warnings import warn
from interval import interval_newton, NoZeroException

from sage.rings.real_mpfi import RealIntervalField, is_RealIntervalFieldElement
from dynamic import normalized, Mod1Dynamic, QuotientedLinearDynamic, PiecewiseExpandingDynamic
load('binsearch2.spyx')

class Ulam(Basis):
	r"""
	Basis for the Ulam method
	
	Basis elements are the functions :math:`u_j = \frac{1}{|x_{j+1}-x_j|} 1_{[x_j,x_{j+1}]}`
	
	Dual basis elements are the integrals :math:`V_i = \int_{[x_i,x_{i+1}]}`.
	
	We have to normalize in this way because we want stochasticity in the resulting matrix, so we require \sum_i V_i to be the integral over the whole domain. Any other normalization would get into trouble when using unequal intervals.
	
	There is an important assumption that we are using at the moment: that all points have the same number of preimages and all branches of the function are monotonic.
	
	"""
	
	def __init__(self, partition):
		check_partition(partition)
		self.partition = partition
	
	def __len__(self):
		return len(self.partition) - 1
	
	def evaluate(self, i, x):
		# updated to support canonical representatives
		field = x.parent()
		if is_inside(x, self.partition, i):
			return field(1)
		elif overlap(x, self.partition, i):
			return field(0, 1)
		else:
			return field(0)
		
	def nonzero_on(self, I):
		return nonzero_on(I, self.partition)

			
	def dual_composed_with_dynamic(self, dynamic, epsilon):
		r"""
		For each `i`, yields the set of `dynamic.nbranches` pairs :math:`\{(f_k^{-1}(x_i),f_k^{-1}(x_{i+1})) \colon f_k \in \text{branch inverses}\}`.
		Note that this is a *pair* and not an *interval*. This is necessary, because we have
		uncertainty on the values of the lower and upper bound of the interval, and
		interval arithmetic is used to represent that. If we were to "flatten" the interval,
		we'd lose the interpretation of its minimum and maximum length.

		Repeat after me: interval arithmetic is *not* for representing Ulam intervals.
		"""

		if isinstance(dynamic, Mod1Dynamic):

			# We rely on the fact that the preimages are returned in increasing order here,
			# so be careful if you wish to implement this for IterateDynamic() as well.

			# TODO: this would probably need to be done in a completely different way
			# to be more general, for instance using the multi-interval generalization
			# of interval_newton.
			# It is probably going to be simpler to deal with these issues when we
			# generalize this thing to domains different from [0,1].

			n = len(self)
			preimages_0 = list(dynamic.preimages(self.partition[0], epsilon))
			preimages_xi = preimages_0
			for i in range(len(self) - 1):
				preimages_xi_plus_one = list(dynamic.preimages(self.partition[i+1], epsilon))
				for K in zip(preimages_xi, preimages_xi_plus_one):
					yield (i, K)
				preimages_xi = preimages_xi_plus_one

			# The last interval is special because we need the preimages of "1 unquotiented", not 0, 
			# and simply evaluating f would quotient incorrectly and give again preimages_0.
			preimages_1 = preimages_0[1:]
			preimages_1.append(dynamic.field(1))
			for K in zip(preimages_xi, preimages_1):
				yield (n-1, K)

		elif isinstance(dynamic, PiecewiseExpandingDynamic):
			for k in range(dynamic.nbranches):
				f = dynamic.functions[k]
				fprime = dynamic.derivatives[k]
				if not fprime((dynamic.grid[k]+dynamic.grid[k+1])/2) > 0:
					raise NotImplementedError, 'Only increasing functions for now'
				fmin = f(dynamic.field(dynamic.grid[k]))
				fmax = f(dynamic.field(dynamic.grid[k+1]))
				fdomain = dynamic.field(dynamic.grid[k], dynamic.grid[k+1])
				frange = fmin.union(fmax)
				# integer_shift must be an interval because we need to make arithmetic with it
				integer_shift = dynamic.field(fmin.floor().lower())
				lowest_i = binsearch2((fmin - integer_shift).lower(), self.partition)
				assert lowest_i > -1
				a = dynamic.field(dynamic.grid[k])
				i = lowest_i
				while True:
					next_i = i + 1
					if next_i == len(self):
						next_i = 0
						integer_shift = integer_shift + 1
					try:
						b =	interval_newton(f, fprime, fdomain, integer_shift + self.partition[next_i], epsilon)
					except NoZeroException:
						yield (i, (a, dynamic.field(dynamic.grid[k+1])))
						break
					yield (i, (a, b))
					i = next_i
					a = b
		else:
			raise NotImplementedError, 'We need guaranteed monotonicity intervals in the dynamic for this to work'

	def project_dual_element(self, dual_element):
		a, b = dual_element
		for j in self.nonzero_on(a.union(b)):
			# ensure that they are sorted
			A = a.min(b)
			B = a.max(b)
			
			x = A.parent()(self.partition[j])
			y = A.parent()(self.partition[j+1])
			
			# compute endpoints of the intersection
			lower = A.max(x)
			upper = B.min(y)
			yield (j, (upper - lower).max(0) / (y-x))

	def sanitize_perron_vector(self, v):
		"""
		Uses theoretical knowledge on the Perron vector (e.g., realness, positivity...) to correct numerical output.
		
		This gets called only before computing the rigorous residual, so don't worry, we're not cheating.
		
		Args:
			v (numpy vector):
		Returns:
			v (numpy vector): the same vector but "adjusted". In particular, we require that it is scaled so that the invariant measure is positive and has mass 1 (numerically)
		"""
		
		v = np.real(v)
		v = v / sum(v) # makes sure that the components are (almost all) positive
		v[v<0] = 0
		v = v / np.linalg.norm(v, ord=1) #probably faster than our norm1(v), and than self.integral(v)
		return v

	def evaluate_integral(self, i, ring=RIF):
		"""
		Integral of the `i`th basis function
		"""
		return ring(1)
	
	def mat_vec(self,PP,v):
		return PP*v

class UlamL1(Ulam):
	"""
	Combines Ulam with norm-1 bound estimation.
	"""
	
	def contracting_pairs(self):
		"""		
		Return vectors with u_i[0]=u_i[i+1]=1/2 and the rest 0, and s_i=1/2.
		"""
		n = len(self)
		
		for i in range(n-1):
			u = np.zeros(n)
			u[0] = 1/2
			u[i+1] = -1/2
			yield (u, 1/2)
		
	def norm_estimate(self, v):
		return norm1(v, rnd='RNDU')
	
	def matrix_norm_estimate(self, PP):
		return matrix_norm1(PP, rnd='RNDU')
	
	def rigorous_norm(self, v):
		return v.norm(1)
	
	def bound_on_norms_of_powers(self, dynamic, project_left = False, project_right = False):
		return RealField(rnd='RNDU')(1)
			
	def matrix_norm_diameter(self, P):		
		w = vector(RealField(rnd='RNDU'),P.dimensions()[1])
		for (ij, val) in P.dict().iteritems():
			w[ij[1]] += val.absolute_diameter() #operations are rounded up because the LHS is so
		return max(w)
	
	def numerical_error(self, dynamic, P, PP):
		"""
		Return a bound for the error amplification coefficients.
		
		Args:
			dynamic
			P (sage interval matrix): the computed discretized matrix
			PP (scipy sparse matrix): the matrix that we will use for the floating-point arithmetic
		
		Returns:
			K (float): a constant such that :math:`\|Fl(Pv)-Pv\| \leq K\|v\|` for all vectors :math:`v`
		"""
		Rup = RealField(rnd='RNDU')
		Rdown = RealField(rnd='RNDD')
		
		D = self.matrix_norm_diameter(P)
		
		gamma = Rup(np.finfo(float).eps) * max_nonzeros_per_row(P);
		gamma = gamma / (Rdown(1)-gamma)
		
		K = gamma * self.matrix_norm_estimate(PP) + D
		return K

	def lasota_yorke_constants(self, dynamic):
		"""
		Return Lasota-Yorke constants
		
		This is meant to replace `dynamic.lasota_yorke_constants()` with a more norm-agnostic packaging
		
		Returns:
			(lambda, B) (pair of real intervals): such that :math:`\|Lf\|_s \leq \lambda \|f\|_s + B\|f\|`
		"""
		warn('Deprecated: we are trying to move to dfly')

		if not dynamic.expansivity() > 1:
			raise ValueError, 'Your map is not expansive. A Lasota-Yorke inequality probably doesn''t exist.'

		if dynamic.is_continuous_on_torus():
			return (1/dynamic.expansivity(), dynamic.distorsion()+1)
		else:
			if dynamic.expansivity() > 2:
				raise NotImplementedError, 'since this is deprecated, not implementing it'
			else:
				raise ValueError, "We don't know how to make a LY inequality if lambda > 1/2. Try with 'iteratewithlasotayorke'."
			raise NotImplementedError, "L-Y constants are implemented only for functions which are continuous on the torus. Define a method is_continuous_on_torus() if your map is so."

	def semidiscrete_lasota_yorke_constants(self, dynamic):
		"""
		Return Lasota-Yorke constants of :math:`L\Pi`

		Returns:
			(lambda, B) (pair of real intervals): such that :math:`\|L\Pi f\|_s \leq \lambda \|f\|_s + B\|f\|`
		"""
		return self.lasota_yorke_constants(dynamic)

	def mixed_matrix_norm_approximation_error(self, dynamic):
		r"""
		Return a bound on :math:`\|L-L_h\|_{s \to w}`

		That is, :math:`\sup \frac{\|(L-L_h)f\|_w}{\|f\|_s}`, :math:`\|\cdot\|_w` and :math:`\|\cdot\|_s` 
		being the weak and strong norm, respectively.
		"""

		return partition_diameter(self.partition)*2

	def mixed_matrix_norm_distance(self, other):
		"""
		Distance (in strong->weak norm) between the projection associated to this basis and to another one.

		Args:
			other: another basis
		Returns:
			a bound to :math:`\|\Pi_{self} - \Pi_{other}\|_{s\to w}`, where :math:`\Pi` denotes the projection operator
			If one partition is a refinement of the other, returns delta_C, else delta_C + delta_F.
		"""

		if not type(other) == type(self):
			raise NotImplementedError

		fine = self.partition
		coarse = other.partition
		if len(fine) < len(coarse):
			fine, coarse = coarse, fine

		if is_refinement(fine, coarse):
			return partition_diameter(coarse)
		else:
			return partition_diameter(coarse) + partition_diameter(fine)

	def strong_to_weak_norm_equivalence(self):
		"""
		Return the constant required to estimate a strong norm with a weak norm (in the discretized space)

		Since we have a finite-dimensional discretized space, there is a constant
		such that :math:`\|f\|_s \leq C \|f\|_w` for each `f` in the discrete space.

		Returns:
			C (real with RNDU):
		"""
		
		return RealField(rnd='RNDU')(2) / partition_minimum_diameter(self.partition)
		
	def iterate_with_lasota_yorke(self, D):
		if not D.is_continuous_on_torus():
			raise NotImplementedError
		if D.expansivity() > 1:
			return D
		else:
			raise ValueError, "No iterate of this map admits a Lasota-Yorke inequality"

	def dfly(self, dynamic, discrete=False, n=None):
		"""
		Return constants of the dfly inequality :math:`\|L^n f\|_s \leq A \lambda^n \|f\|_s + B\|f\|_w` 

		This function should make `lasota_yorke_constants` and `iterate_with_lasota_yorke` obsolete.
		
		Input:
			Dynamic:
			discrete: if True, returns constants for the projected operator instead
			n: if non-None, may return a (possibly stronger) inequality that holds only for
			   the given value of n

		Returns:
			(A, B, lambda): tuple of intervals
		"""
		if not dynamic.expansivity() > 1:
			raise ValueError, 'Your map is not expansive. A Lasota-Yorke inequality probably doesn''t exist.'

		if dynamic.is_continuous_on_torus():
			# note that this bound is ok for both continuous and discrete, since all projections have norm 1.
			lam = 1/dynamic.expansivity()
			if n==1:
				return (dynamic.field(1), dynamic.distorsion(), lam, )
			else:
				return (dynamic.field(1), dynamic.distorsion() / (1 - lam), lam, )
		elif isinstance(dynamic, QuotientedLinearDynamic):
			# Galatolo-Nisoli, Remark 8 with d_0=0, d_1=1. TODO: check with Isaia.
			lam = 1/dynamic.expansivity()
			if not 2*lam < 1:
				raise ValueError, 'Try with iterate_with_lasota_yorke(); maybe there is a sharper bound though'
			if n==1:
				return (dynamic.field(1), 2*dynamic.distorsion() + 2, 2*lam, )
			else:
				return (dynamic.field(1), (2*dynamic.distorsion() + 2) / (1-2*lam), 2*lam, ) 
		elif isinstance(dynamic, PiecewiseExpandingDynamic):
			# Galatolo-Nisoli, Remark 8 with d_0=0, d_1=1. 
			lam = 1/dynamic.expansivity()
			t = dynamic.field(0)
			for i in range(dynamic.nbranches):
				t = t.max(2/(dynamic.field(dynamic.grid[i+1]) - dynamic.field(dynamic.grid[i])))
			if n==1:
				return (dynamic.field(1), 2*dynamic.distorsion() + t, 2*lam, ) 
			else:
				return (dynamic.field(1), (2*dynamic.distorsion() + t) / (1-2*lam), 2*lam, ) 
		else:
			raise NotImplementedError, "L-Y constants are implemented only for functions which are continuous on the torus. Define a method is_continuous_on_torus() if your map is so."

	def projection_error(self):
		"""
		Returns a bound on the projection error :math:`\|\Pi-I\|_{s\to w}`
		"""
		return partition_diameter(self.partition)

	def coefficients_approximation_inequality(self,dynamic):
		""" This function returns the coefficients C,D such that
		||L^n -L^n _{\delta}f||\leq \delta C ||f||_s + n\delta D ||f||_w
				args:
					dynamic: the dynamic we are studying
		"""
		[lambda_m, C_m]=self.lasota_yorke_constants(dynamic)
		
		Interval=lambda_m.parent()
		lam= lambda_m
		B= C_m
		
		A= 1
		M= self.bound_on_norms_of_powers(dynamic)
		
		P= Interval(1)
		K= Interval(1)
		
		C=K*M*(A*lam+P)*A/(1-lam)
		D=K*M*B*(A*lam+P+M)
		return C, D

	def spectral_radius_and_autonorm(self,dynamic,n,lambda_2):
		"""
		This function returns bounds on the spectral radius of L
		and the coefficients of the (a,b)-autonorm as presented in [GaNiSa]
		used to estimate in an elementary way the speed of convergence to equilibrium
		by knowing that ||L^n_{\delta}||_w \leq \lambda_2
				args:
					dynamic: the dynamic
					lambda_2 :
					n :
		"""
		
		[lambda_m, C_m]=self.lasota_yorke_constants(dynamic)
		
		Interval=lambda_m.parent()
		lam= lambda_m
		B= C_m/(1-lambda_m)
		
		A= 1
		delta = Interval(1)/(len(self))
		
		C,D=self.coefficients_approximation_inequality(dynamic)
		
		# building the matrix
		
		a=A*lam**n
		b=B
		c=delta*C
		d=delta*n*D+lambda_2
		
		print a, b, c, d
		
		return [rho(a,b,c,d),autonorm(a,b,c,d)]

from rigorous_integral import riemann_integral_on_interval

class UlamL1BasisUtilities(UlamL1):
	""" This class contains some utilities needed if we want to do
	some computation on function when projected on the Ulam basis
	it must be considered a superset of UlamL1and it is not needed for 
	the assembly of the discretized operator nor the estimates.
	It is only used for finer study of the system, as computing the
	diffusion coefficients or integrating observables against the
	computed density"""
		
	def project_on_ulam_basis(self,observable,err, report_frequence=65536, numpy_array_output = False):
		""" This function projects an observable on the Ulam basis,
		with a rigorous error bound
			args:
				obs: the observable
				err: the error
		"""
		n=len(self)
		error=err/n
		if numpy_array_output is False:
			v=VectorSpace(RIF,n)(0)
			final_err=RealNumber(0,rnd='RNDU')
			i_min=(((observable.support()).lower())*n).floor()
			i_max=(((observable.support()).upper())*n).ceil()
			for i in range(i_min,i_max):
				if (i%report_frequence==0):
					print "Projected the first ", i, " coordinates"
				x=RIF(self.partition[i],self.partition[i+1])
				coeff=riemann_integral_on_interval(observable,x,error)
				#final_err=max(final_err,coeff.absolute_diameter())
				v[i]=coeff
			return v
		else:
			v=np.zeros(n)
			final_err=RealNumber(0,rnd='RNDU')
			i_min=(((observable.support()).lower())*n).floor()
			i_max=(((observable.support()).upper())*n).ceil()
			for i in range(i_min,i_max):
				if (i%report_frequence==0):
					print "Projected the first ", i, " coordinates"
				x=RIF(self.partition[i],self.partition[i+1])
				coeff=riemann_integral_on_interval(observable,x,error)
				final_err=max(final_err,coeff.absolute_diameter())
				v[i]=coeff.center()
			return v,final_err

	def integral_observable_against_density(self,observable,density,err,report_frequence=65536):
		""" This function computes the integral of an observable againsta a piecewise constant density
			with a rigorous error bound
			args:
				obs: the observable
				debsity: the density
				err: the error
		"""
		n=len(self)
		error=err/n
		result=RIF(0)
		
		i_min=(((observable.support()).lower())*n).floor()
		i_max=(((observable.support()).upper())*n).ceil()
		
		for i in range(i_min,i_max):
			if (i%report_frequence==0):
				print "Integrated the first ", i, " coordinates"
			x=RIF(self.partition[i],self.partition[i+1])
			result+=riemann_integral_on_interval(observable,x,error)*density[i]/x.absolute_diameter()
			
		return result
		
	def pointwise_sum(self,v,w):
		return v+w
		
	def pointwise_product(self,v,w):
		n=len(self)
		z=VectorSpace(RIF,n)(0)
		for i in range(n):
			x=RIF(self.partition[i],self.partition[i+1])
			z[i]=v[i]*w[i]/x.absolute_diameter()
		return z
			
	def pointwise_product_in_place_of_first_argument_mutable(self,v,w):
		"""
		Pointwise product, please note it is made in place v
		therefore saving memory
		args:
			v: the mutable vector, the projected observable
			w: the const density
		"""
		n=len(self)
		for i in range(n):
			x=RIF(self.partition[i],self.partition[i+1])
			v[i]=(v[i]*w[i]/x.absolute_diameter())

	def integral_given_projection_and_density(self,proj,density):
		"""
		This function computes the integral of an observable given 
		its projection on the basis and an approximation of the density
		Tries to use at least memory as possible
		args:
			proj: projection of observable on the Ulam basis
			density: approximation of the density
		"""
		result=RIF(0)
		n=len(self)
		for i in range(n):
			x=RIF(self.partition[i],self.partition[i+1])
			result+=(proj[i]*density[i]/x.absolute_diameter())
		return result
		
	def pointwise_product_observable_density_outputs_vector_of_double_and_error(self,observable,density):
		n=len(self)
		w=np.zeros(n)
		error=0
		for i in range(n):
			x=RIF(self.partition[i],self.partition[i+1])
			val=riemann_integral_on_interval(observable,x,error)*density[i]/x.absolute_diameter()
			v[i]=val.center()
		return w,error
		
	def Lebesgue_integral_pointwise_product(self,projected_1,projected_2, equispaced_partition = False, projection_np_array = False):
		"""
		Pointwise product, please note it is made in place v
		therefore saving memory
		args:
			projected_1: projected vector 1
			projected_2: projected vector 2
		"""
		n=len(self)
		integral=RIF(0)
		if equispaced_partition is True:
			if projection_np_array is True:
				integral=projected_2.dot(projected_1)
				return integral*n
			if projection_np_array is False:
				integral=projected_2.dot_product(projected_1)
				return integral*n
			
		if equispaced_partition is False:
			for i in range(n):
				x=RIF(self.partition[i],self.partition[i+1])
				integral+=projected_1[i]*projected_2[i]/x.absolute_diameter()
		return integral
	
	def test_Lebesgue_integral_pointwise(self):
		proj_1=np.random.uniform(0,1,len(self))
		proj_2=np.random.uniform(0,1,len(self))
		int_1=self.Lebesgue_integral_pointwise_product(proj_1,proj_2)
		int_2=self.Lebesgue_integral_pointwise_product(proj_1,proj_2,equispaced_partition = True , projection_np_array = True)
		return int_1,int_2
