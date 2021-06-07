from sage.all import RIF

__all__=["riemann_integral_on_interval"]

def riemann_integral_on_interval(observable,interval,error):
	""" This function returns the integral of an observable on an interval
		This function masks the riemann_sum_on_interval_C* functions,
		it gets the information about the regularity of the observable 
		from the regularity method of the class observable and choses 
		the fastest integrator.
		Please note that the integrators are not imported by default,
		and should be used only through this function.
		args:
			observable:
			interval:
			error: please note that this error is an euristic, the final error may be different
	"""
	regularity=observable.regularity()
	if (regularity is 1):
		return riemann_sum_on_interval_C1(observable,interval,error)
	if (regularity is 2):
		return riemann_sum_on_interval_C2(observable,interval,error)
	if (regularity is 3):
		return riemann_sum_on_interval_C3(observable,interval,error)
	
	

def riemann_sum_on_interval_C1(observable,interval,error):
	""" This function computes the integral of a C1 observable on an interval
		Please note that it uses the fact that the observable is C1
		args:
			observable : the observable
			interval : the interval
			error : a bound on the size of the integral, please note that it is not the output error
	"""
	err=error/16
	x_mid=interval.center()
	rad=interval.absolute_diameter()/2

	f_mid=observable.value(RIF(x_mid))
	f_prime_mid=observable.derivative(RIF(x_mid))
	epsilon=(observable.derivative(interval)-f_prime_mid).magnitude()
	
	val=2*(f_mid*rad+RIF(-epsilon,epsilon)*((rad**2)/2))
	
	if (val.absolute_diameter()<error):
		return val
	else:
		interval1,interval2=interval.bisection()
		sum_left=riemann_sum_on_interval_C1(observable,interval1,error)
		sum_right=riemann_sum_on_interval_C1(observable,interval2,error)
		return (sum_left+sum_right)
		

def riemann_sum_on_interval_C2(observable,interval,error):
	""" This function computes the integral of a C1 observable on an interval
		Please note that it uses the fact that the observable is C1
		args:
			observable : the observable
			interval : the interval
			error : a bound on the size of the integral, please note that it is not the output error
	"""
	
	
	err=error/16
	x_mid=interval.center()
	rad=interval.absolute_diameter()/2

	f_mid=observable.value(RIF(x_mid))
	f_second_mid=observable.second_derivative(RIF(x_mid))/2
	epsilon=(observable.second_derivative(interval)/2-f_second_mid).magnitude()
	
	val=2*(f_mid*rad+f_second_mid*((rad**3)/3)+RIF(-epsilon,epsilon)*((rad**3)/3))
	
	if (val.absolute_diameter()<error):
		return val
	else:
		interval1,interval2=interval.bisection()
		left=riemann_sum_on_interval_C2(observable,interval1,error)
		right=riemann_sum_on_interval_C2(observable,interval2,error)
		return (left+right)

def riemann_sum_on_interval_C3(observable,interval,error):
	""" This function computes the integral of a C1 observable on an interval
		Please note that it uses the fact that the observable is C1
		args:
			observable : the observable
			interval : the interval
			error : a bound on the size of the integral, please note that it is not the output error
	"""
	
	
	err=error/16
	x_mid=interval.center()
	rad=interval.absolute_diameter()/2

	f_mid=observable.value(RIF(x_mid))
	f_second_mid=observable.second_derivative(RIF(x_mid))/2
	f_third_mid=observable.third_derivative(RIF(x_mid))/6
	epsilon=(observable.third_derivative(interval)/6-f_third_mid).magnitude()
	
	val=2*(f_mid*rad+f_second_mid*((rad**3)/3)+RIF(-epsilon,epsilon)*((rad**4)/4))
	
	if (val.absolute_diameter()<error):
		return val
	else:
		interval1,interval2=interval.bisection()
		left=riemann_sum_on_interval_C2(observable,interval1,error)
		right=riemann_sum_on_interval_C2(observable,interval2,error)
		return (left+right)

