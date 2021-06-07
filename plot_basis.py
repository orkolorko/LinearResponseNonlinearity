from quintic_basis import phi,phi_prime,phi_second,psi,psi_prime,psi_second
from sage.all import plot, RIF
	
def test_plot_basis_phi(i,n):
	pl=plot(lambda x: phi(RIF(x),i,n).center(), 0, 1, legend_label='phi')+\
	plot(lambda x: phi_prime(RIF(x),i,n).center(), 0, 1, legend_label='phi_prime',color='red')#+\
	#plot(lambda x: phi_second(RIF(x),i,n).center(), 0, 1, legend_label='phi_second',color='green')
	pl.show()

def test_plot_basis_psi(i,n):
	pl=plot(lambda x: psi(RIF(x),i,n).center(), 0, 1, legend_label='psi')+\
	plot(lambda x: psi_prime(RIF(x),i,n).center(), 0, 1, legend_label='psi_prime',color='red')#+\
	#plot(lambda x: psi_second(RIF(x),i,n).center(), 0, 1, legend_label='psi_second',color='green')
	pl.show()

def test_plot_basis_phi_composed_dynamic(i,n,dynamic):
	pl=plot(lambda x: phi(dynamic.f(RIF(x)),i,n).center(), 0, 1, legend_label='phi\circ T')+\
	plot(lambda x: (phi_prime(dynamic.f(RIF(x)),i,n)*dynamic.fprime(RIF(x))).center(), 0, 1, legend_label='phi_prime',color='red')
	pl.show()

def test_plot_basis_phi_composed_unquotiented_dynamic(i,n,dynamic):
	pl=plot(lambda x: phi(dynamic.f_unquotiented(RIF(x)),i,n).center()+phi(dynamic.f_unquotiented(RIF(x)),i+n,n).center(), 0, 1, legend_label='phi\circ T')+\
	plot(lambda x: (phi_prime(dynamic.f_unquotiented(RIF(x)),i,n)*dynamic.fprime(RIF(x))).center(), 0, 1, legend_label='phi_prime',color='red')
	pl.show()

def test_plot_basis_psi_composed_unquotiented_dynamic(i,n,dynamic):
	pl=plot(lambda x: phi(dynamic.f_unquotiented(RIF(x)),i,n).center(), 0, 4, legend_label='phi\circ T')+\
	plot(lambda x: (psi_prime(dynamic.f_unquotiented(RIF(x)),i,n)*dynamic.fprime(RIF(x))).center(), 0, 1, legend_label='phi_prime',color='red')
	pl.show()

def test_plot_basis_psi_composed_dynamic(i,n,dynamic):
	pl=plot(lambda x: psi(dynamic.f(RIF(x)),i,n).center(), 0, 1, legend_label='psi\circ T')+\
	plot(lambda x: (psi_prime(dynamic.f(RIF(x)),i,n)*dynamic.fprime(RIF(x))).center(), 0, 1, legend_label='psi\circ T_prime',color='red')
	pl.show()

def test_plot_basis_phi_circ_T(i,n,dynamic):
	pl=plot(lambda x: phi_circ_T(RIF(x),i,n,dynamic).center(), 0, 1, legend_label='phi\circ T')+\
	plot(lambda x: phi_circ_T_prime((RIF(x)),i,n,dynamic).center(), 0, 1, legend_label='phi_prime',color='red')
	pl.show()

def test_plot_basis_tau(n):
	pl=plot(lambda x: (tau(RIF(x),n)).center(), 0, 1, legend_label='phi')+\
	plot(lambda x: (tau_prime(RIF(x),n)).center(), 0, 1, legend_label='phi_prime',color='red')
	pl.show()
