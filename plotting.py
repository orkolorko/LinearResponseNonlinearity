"""
Istructions to plot invariant measures and the like.

These functions return plots. One can add plots with `p1 + p2` and show them with `show(p)`
"""

from scipy.sparse.linalg import eigs
from sparse import *
from sage.all import scatter_plot, real, imag, RIF, bar_chart, plot_step_function, sqrt, plot, line
import numpy as np
from partition import step_function


def plot_spectrum(P, neigs=30, **kwargs):
	"""
	Plots (selected eigenvalues of) the spectrum of a matrix.
	
	Args:
		P (sage interval matrix, sage matrix or scipy matrix):
		neigs (int, default 30): number of eigenvalues to plot. If neigs==0, plots all eigenvalues (very slow for matrices larger than ca. 1024).
		
	Returns:
		plot
	"""
	
	PP = sage_sparse_to_scipy_sparse(P)
	if neigs == 0:
		L, V = np.linalg.eig(PP.todense())
	else:
		L, V = eigs(PP, neigs)

	return scatter_plot([(real(x), imag(x)) for x in L], **kwargs)
	
def plot_dynamic(D, **kwargs):
	"""
	Plot a dynamic as a map [0,1] -> [0,1]
	
	Args:
		D (Dynamic):
		
	Returns:
		plot
	"""
	
	if not 'legend_label' in kwargs:
		kwargs['legend_label'] = 'Map T'
	return plot(lambda x: D.f(RIF(x)).center(), 0, 1, **kwargs)

def plot_measure(v, partition, error=None, basis='Ulam', norm='L1', **kwargs):
	"""
	Plot an invariant measure
	
	Args:
		v (numpy array): basis elements of the invariant measure in the specified basis
		partition: the partition of [0,1] under which v has been produced
		basis (string, default 'Ulam'): see above
		error (real or None (default)): error to show in the plot, in the specified norm
		norm (string, default 'L1'): norm to use for the error
	"""

	if not 'legend_label' in kwargs:
		kwargs['legend_label'] = 'Invariant measure'
	if not 'color' in kwargs:
		kwargs['color'] = 'red'
	
	if basis == 'Ulam' or basis == 'ulam':
		steps = tuple(step_function(v, partition))
		p = plot_step_function(steps, **kwargs)
	elif basis == 'hat':
		points = zip(partition, np.append(v, v[0]))
		p = line(points, **kwargs)
	else:
		raise ValueError, "invalid basis type"
		
	if error is not None:
		if norm == 'L1':
			p += bar_chart([sqrt(error)], width=sqrt(error), color='green', legend_label='Area of the total error')
		elif norm == 'Linf' or norm == 'C0':
			p += line([(x, y-error) for (x, y) in points], color='green', legend_label='error bounds')
			p += line([(x, y+error) for (x, y) in points], color='green')
		else:
			raise ValueError, "invalid norm type"
	
	return p

def plot_basis_element(basis, i, **kwargs):
	if not 'legend_label' in kwargs:
		kwargs['legend_label'] = 'phi'

	return plot(lambda x: basis.evaluate(i, RIF(x)).center(), 0,1, **kwargs)

def plot_basis_element_composed_with_dynamic(basis, i, dynamic, **kwargs):
	if not 'legend_label' in kwargs:
		kwargs['legend_label'] = 'phi \circ T'
	return plot(lambda x: basis.evaluate(i, dynamic.f(RIF(x))).center(), 0,1, **kwargs)

def plot_perron_operator(basis, i, dynamic, **kwargs):
	f = lambda x: basis.evaluate(i, RIF(x))
	if not 'legend_label' in kwargs:
		kwargs['legend_label'] = 'L\phi'
	return plot(lambda x: dynamic.perron_operator(f, x, epsilon=1e-10).center(), 0, 1, **kwargs)
