from partition import *
from basis import Basis

from sage.all import RealIntervalField, RIF, RealNumber, RealField, vector, floor, plot
from sparse import max_nonzeros_per_row, matrix_norminf, sage_sparse_to_scipy_sparse, sparse_matvec
from scipy.sparse import csr_matrix
from spectral_radius import *
import numpy as np
from polyval_lib import polyval
from rump_implementation import rump_sparse_mat_vec, rump_dot, rump_sparse_mat_vec_unif_bound
from generic_estimator import extend_power_norms

from sage.all import load
load('binsearch2.spyx')

RUP=RealField(rnd='RNDU')

class Spline_Basis_C1(Basis):
    """
    Abstract class for the projection from C^1 to C^0

    :math:`\phi_j` is the function so that :math:`\phi_j(x_i)=\delta_{ij}`, so :math:`\phi_0` has support :math:`[x_{n-1},1] \cup [0,x_1]`.

    """
    
    def __init__(self, partition):
        check_partition(partition)
        self.partition = partition
        #self._high_prec_integrals_generated=False
        
        n=len(self)
        self._integrals_double=np.zeros(n)
        self._integrals_double_radius=np.zeros(n)
        for i in range(n):
            self._integrals_double[i]=self.evaluate_integral(i).center()
            self._integrals_double_radius[i]=self.evaluate_integral(i).absolute_diameter()/2
        self._integrals_double[n-1]=0
        
    def __len__(self):
        return 2*len(self.partition)+1

    def evaluate_phi(self,i,x,partsize=None):
        """
        This is the implementation of the elements of the elements of the partition of unity
        i: identifies the element of the partition phi_i(a_i)=1
        """
        if partsize is None:
            partsize = len(self.partition) - 1
        
        Interval_Field=x.parent()
        tilde_x = x*partsize-i
        #if x is not in [-1,1], return 0
        
        try:
            tilde_x=tilde_x.intersection(Interval_Field(-1,1))
        except ValueError:
            return Interval_Field(0)
        try:
            tilde_x_neg=tilde_x.intersection(Interval_Field(-1,0))
            #val_neg = polyval([Interval_Field(1),Interval_Field(0),Interval_Field(-3),Interval_Field(-2)], tilde_x_neg)
            val_neg = polyval([Interval_Field(1),Interval_Field(0),Interval_Field(0),Interval_Field(10),Interval_Field(15),Interval_Field(6)], tilde_x_neg)
        except ValueError:
            val_neg=None
        try:
            tilde_x_pos=tilde_x.intersection(Interval_Field(0,1))
            #val_pos = polyval([Interval_Field(1),Interval_Field(0),Interval_Field(-3),Interval_Field(2)], tilde_x_pos)
            val_pos = polyval([Interval_Field(1),Interval_Field(0),Interval_Field(0),Interval_Field(-10),Interval_Field(15),Interval_Field(-6)], tilde_x_pos)
        except ValueError:
            val_pos=None
        if val_pos is None:
            return val_neg
        if val_neg is None:
            return val_pos
        return val_neg.union(val_pos)
    
    def evaluate_phi_prime(self,i,x,partsize=None):
        """
        This is the implementation of the elements of the elements of the partition of unity
        i: identifies the element of the partition phi_i(a_i)=1
        """
        if partsize is None:
            partsize = len(self.partition) - 1
        
        Interval_Field=x.parent()
        tilde_x = x*partsize-i
        #if x is not in [-1,1], return 0
        
        try:
            tilde_x=tilde_x.intersection(Interval_Field(-1,1))
        except ValueError:
            return Interval_Field(0)
        try:
            tilde_x_neg=tilde_x.intersection(Interval_Field(-1,0))
            #val_neg = polyval([Interval_Field(0),Interval_Field(-6),Interval_Field(-6)], tilde_x_neg)
            val_neg = polyval([Interval_Field(0),Interval_Field(0),Interval_Field(30),Interval_Field(60),Interval_Field(30)], tilde_x_neg)
        except ValueError:
            val_neg=None
        try:
            tilde_x_pos=tilde_x.intersection(Interval_Field(0,1))
            #val_pos = polyval([Interval_Field(0),Interval_Field(-6),Interval_Field(6)], tilde_x_pos)
            val_pos = polyval([Interval_Field(0),Interval_Field(0),Interval_Field(-30),Interval_Field(60),Interval_Field(-30)], tilde_x_pos)
        except ValueError:
            val_pos=None
        if val_pos is None:
            return val_neg*partsize
        if val_neg is None:
            return val_pos*partsize
        return val_neg.union(val_pos)*partsize
    
    def evaluate_nu(self,i,x,partsize=None):
        """
        This is the implementation of the elements of the elements of the partition of unity
        i: identifies the element of the partition phi_i(a_i)=1
        """
        if partsize is None:
            partsize = len(self.partition) - 1
        
        Interval_Field=x.parent()
        tilde_x = x*partsize-i
        #if x is not in [-1,1], return 0
        
        try:
            tilde_x=tilde_x.intersection(Interval_Field(-1,1))
        except ValueError:
            return Interval_Field(0)
        try:
            tilde_x_neg=tilde_x.intersection(Interval_Field(-1,0))
            #val_neg = polyval([Interval_Field(0),Interval_Field(1),Interval_Field(2),Interval_Field(1)], tilde_x_neg)
            val_neg = polyval([Interval_Field(0),Interval_Field(1),Interval_Field(0),Interval_Field(-6),Interval_Field(-8),Interval_Field(-3)], tilde_x_neg)
        except ValueError:
            val_neg=None
        try:
            tilde_x_pos=tilde_x.intersection(Interval_Field(0,1))
            #val_pos = polyval([Interval_Field(0),Interval_Field(1),Interval_Field(-2),Interval_Field(1)], tilde_x_pos)
            val_pos = polyval([Interval_Field(0),Interval_Field(1),Interval_Field(0),Interval_Field(-6),Interval_Field(8),Interval_Field(-3)], tilde_x_pos)
        except ValueError:
            val_pos=None
        if val_pos is None:
            return val_neg/partsize
        if val_neg is None:
            return val_pos/partsize
        return val_neg.union(val_pos)/partsize
    
    def evaluate_nu_prime(self,i,x,partsize=None):
        """
        This is the implementation of the elements of the elements of the partition of unity
        i: identifies the element of the partition phi_i(a_i)=1
        """
        if partsize is None:
            partsize = len(self.partition) - 1
        
        Interval_Field=x.parent()
        tilde_x = x*partsize-i
        #if x is not in [-1,1], return 0
        
        try:
            tilde_x=tilde_x.intersection(Interval_Field(-1,1))
        except ValueError:
            return Interval_Field(0)
        try:
            tilde_x_neg=tilde_x.intersection(Interval_Field(-1,0))
            #val_neg = polyval([Interval_Field(1),Interval_Field(4),Interval_Field(3)], tilde_x_neg)
            val_neg = polyval([Interval_Field(1),Interval_Field(0),Interval_Field(-18),Interval_Field(-32),Interval_Field(-15)], tilde_x_neg)
        except ValueError:
            val_neg=None
        try:
            tilde_x_pos=tilde_x.intersection(Interval_Field(0,1))
            #val_pos = polyval([Interval_Field(1),Interval_Field(-4),Interval_Field(3)], tilde_x_pos)
            val_pos = polyval([Interval_Field(1),Interval_Field(0),Interval_Field(-18),Interval_Field(32),Interval_Field(-15)], tilde_x_pos)
        except ValueError:
            val_pos=None
        if val_pos is None:
            return val_neg
        if val_neg is None:
            return val_pos
        return val_neg.union(val_pos)
    
    def evaluate_tau(self,x):
        """
        Evaluates the function tau on a point x
        Args:
            partsize (integer) : is the size of the partition. This parameter is now computed automatically from the length of the basis, it is not
                          needed, I kept it only to have similar interfaces for tau and psi
            x (real number) :
        """
        RI=x.parent()
        
        if (x-0.5).contains_zero():
            a=RI(x.lower())
            b=RI(x.upper())
            y_1=a*(1-a)
            y_2=b*(1-b)
            y=y_1.union(y_2)
            return RI(6)*y.union(0.25)
        else:
            a=RI(x.lower())
            b=RI(x.upper())
            y_1=a*(1-a)
            y_2=b*(1-b)
            return RI(6)*y_1.union(y_2)
    
    def evaluate_tau_prime(self,x):
        """
        Evaluates the function tau on a point x
        Args:
            partsize (integer) : is the size of the partition. This parameter is now computed automatically from the length of the basis, it is not
                          needed, I kept it only to have similar interfaces for tau and psi
            x (real number) :
        """
        RI=x.parent()
        return RI(6)-RI(12)*x
    
    def evaluate(self,i,x):
        """
        This evaluates the i-th element on the range, i.e., tau_i
        """
        n=len(self.partition)
        if (i<n):
            return self.evaluate_phi(i,x)
        if (n<=i) and (i<2*n):
            return self.evaluate_nu(i%n,x)
        if i==2*n:
            return self.evaluate_tau(x)
    
    def evaluate_prime(self,i,x):
        """
        This evaluates the derivative of the i-th element on the range, i.e., tau_i
        """
        n = len(self.partition)
        if (i<n):
            return self.evaluate_phi_prime(i,x)
        if (n<=i) and (i<2*n):
            return self.evaluate_nu_prime(i%n,x)
        if i==2*n:
            return self.evaluate_tau(x)
        
        
    def nonzero_on(self, I):
        """
        This function gives us the index of the elements of the basis of the domain
        that are nonzero on the interval I
        """
        n = len(self.partition)
        partsize = n - 1
        
        jmin = binsearch2(I.lower(), self.partition)-1
        jmax = binsearch2(I.upper(), self.partition)+1
        
        #jmin = 0
        #jmax = partsize
        jmin = max ( 0, jmin )
        jmax = min (jmax,partsize)
        for i in range(jmin, jmax+1):
            yield i
        
        for i in range(jmin, jmax+1):
            yield i+n
        yield 2*n

    def evaluate_integral(self,i,ring=RIF):
        """
        This evaluates the i-th element on the range, i.e., tau_i
        """
        n=len(self.partition)
        partsize = len(self.partition) - 1
        if (i<n):
            if i==0 or i==n-1:
                return ring(1)/(2*partsize)
            else:
                return ring(1)/partsize
        if (n<=i) and (i<2*n):
            if i==n:
                return ring(1)/(10*partsize**2)
            elif i==2*n-1:
                return ring(-1)/(10*partsize**2)
            else:
                return ring(0)
        if i==2*n:
            return ring(1)
        
    def evaluate_linear_combination(self,w,x):
        """
        Evaluates a function written as a linear combination of elements of the basis on the point x
            Args:
                w : the coordinate vector
                x : the point
        """
        n=len(self)
        y=x.parent()(0)
        
        for i in self.nonzero_on(x):
            y+=w[i]*self.evaluate(i,x)
        y+=w[n-1]*self.evaluate(n-1,x)
        return y
    
    def evaluate_linear_prime(self,w,x):
        """
        Evaluates a function written as a linear combination of elements of the basis on the point x
            Args:
                w : the coordinate vector
                x : the point
        """
        n=len(self)
        y=x.parent()(0)
        
        for i in self.nonzero_on(x):
            y+=w[i]*self.evaluate_prime(i,x)
        return y
    
    def evaluate_function(self,w,x):
        """
        Alias for evaluate_linear_combination
        """
        return self.evaluate_linear_combination(w,x)
    
    def plot_function(self,w):
        """plots the function whose coefficients are given by w"""
        pl=plot(lambda x: self.evaluate_function(w,RIF(x)).center(),0,1)
        return pl
    
    def mat_vec(self,PP,w,debug=False):
        """ 
        Matrix vector product associated to the basis
        """
        
        n = len(self)
        v = w.copy()
        
        v = PP.dot(v)
        
        # we compute the error in the integral of the image
        v[n-1]-=self._integrals_double.dot(v)
        
        return v

    def mat_vec_with_radius(self, PP, PP_abs_plus_radius, PP_radius,w,w_radius, debug=False):
        """ 
        Matrix vector product associated to the basis
        """
        
        n = len(self)
        v = w.copy()
        v_radius = w_radius.copy()
        
        v, v_radius = rump_sparse_mat_vec(PP, PP_abs_plus_radius, PP_radius, v, v_radius)
        
        # we compute the error in the integral of the image
        integral, integral_radius = rump_dot(self._integrals_double,self._integrals_double_radius,v,v_radius)
        
        v[n-1]-=integral
        v_radius[n-1]+=integral_radius
        
        return v, v_radius
    
    def dual_composed_with_dynamic(self, dynamic, epsilon):
        """
        For each `i`, yields (successively) the sequence :math:`\{(y,T'(y),T''(y)) \colon y\in T^{-1}(x_i)\}`.
        """
        n=len(self.partition)
        for i in range(len(self)-1):
            x=self.partition[i%n]
            yield (i,[(i, y,dynamic.fprime(y),dynamic.fsecond(y),dynamic.field) for y in dynamic.preimages(x, epsilon) ])
        #Now, to ensure that we compute the integrals
        yield (len(self)-1,[(len(self)-1,len(self),1,0,dynamic.field)])
        
    def project_dual_element(self, dual_element):
        n=len(self.partition)
        partsize = n - 1
        for i, y, Tprimey, Tsecondy, dynamic_ring in dual_element:
            if (y==len(self)) and (Tprimey==1):
                for j in range(len(self)):
                    yield (j,self.evaluate_integral(j))
            else:
                if (i<n):
                    for j in self.nonzero_on(y):
                        val = self.evaluate(j,y)/(Tprimey.abs())
                        yield (j, val)
                if (n<=i) and (i<2*n):
                    for j in self.nonzero_on(y):
                        val_tau_prime = self.evaluate_prime(j,y)
                        T_prime_square = (Tprimey.abs()**2)
                        val_prime=  val_tau_prime/ T_prime_square  
                        val_const= - self.evaluate(j,y)*Tsecondy/(Tprimey.abs()**3)
                        val= val_prime + val_const
                        yield (j, val)

    def contracting_pairs_with_radius(self,vectorspace=None, start=0, stop=None):
        """
        """
        n=len(self.partition)
        if stop==None:
            stop=len(self)-1
            assert 2*n==stop
        
        for i in range(start, stop):
            if (i<n):
                v=np.zeros(len(self))
                v_radius=np.zeros(len(self))

                #we have to normalize, we bound the C1 norm from below
                norm_lower_bound=abs(RIF(1)+RIF(15)*RIF(n-1)/RIF(8)-self.evaluate_integral(i)*RIF(7.5))
                
                #norm_lower_bound=abs(RIF(1)+RIF(15)*RIF(n-1)/RIF(8)-self.evaluate_integral(i)*(RIF(1)+RIF(15)*RIF(n-1)/RIF(4)))
                
                val_phi=RIF(1)/norm_lower_bound
                v[i]=val_phi.center()
                v_radius[i]=val_phi.absolute_diameter()/2
            
                val_tau=-self.evaluate_integral(i)/norm_lower_bound
            
                v[len(self)-1] = val_tau.center()
                v_radius[len(self)-1] = val_tau.absolute_diameter()/2
            
                yield (v,v_radius)
            
            if (n<=i) and (i<2*n):
                v=np.zeros(len(self))
                v_radius=np.zeros(len(self))
                norm_lower_bound=RIF(1)+RIF(n-1)/RIF(5)
            
                val_nu=RIF(1)/norm_lower_bound
                v[i]=val_nu.center()
                v_radius[i]=val_nu.absolute_diameter()/2
            
                val_tau=-self.evaluate_integral(i)/norm_lower_bound
            
                v[len(self)-1] = val_tau.center()
                v_radius[len(self)-1] = val_tau.absolute_diameter()/2
            
                yield (v,v_radius)
    
    def save_csr_matrix_with_radius(self, P, P_radius, name='default'):
        """
        This function saves the csr matrix containing the values
        and the radius
        """
        size=len(self.partition)-1
        np.save('P_C1_'+name+'_'+str(size)+'_data.npy',P.data)
        np.save('P_C1_'+name+'_'+str(size)+'_indices.npy',P.indices)
        np.save('P_C1_'+name+'_'+str(size)+'_indptr.npy',P.indptr)

        np.save('P_C1_radius_'+name+'_'+str(size)+'_data.npy',P_radius.data)
        np.save('P_C1_radius_'+name+'_'+str(size)+'_indices.npy',P_radius.indices)
        np.save('P_C1_radius_'+name+'_'+str(size)+'_indptr.npy',P_radius.indptr)
        
        P_abs_plus_radius=abs(P)+P_radius
        
        np.save('P_C1_abs_plus_radius_'+name+'_'+str(size)+'_data.npy',P_abs_plus_radius.data)
        np.save('P_C1_abs_plus_radius_'+name+'_'+str(size)+'_indices.npy',P_abs_plus_radius.indices)
        np.save('P_C1_abs_plus_radius_'+name+'_'+str(size)+'_indptr.npy',P_abs_plus_radius.indptr)
    
    def load_csr_matrix_with_radius_memmapped(self, name='default'):
        """
        This function loads the csr matrix containing the values
        and the radius as a memmap
        """
        size=len(self.partition)-1
        memmap_P_data = np.load('P_C1_'+name+'_'+str(size)+'_data.npy', mmap_mode='r')
        memmap_P_indices = np.load('P_C1_'+name+'_'+str(size)+'_indices.npy', mmap_mode='r')
        memmap_P_indptr = np.load('P_C1_'+name+'_'+str(size)+'_indptr.npy', mmap_mode='r')

        P_memmapped = csr_matrix((memmap_P_data,memmap_P_indices,memmap_P_indptr),(len(self),len(self)))

        #loading the matrix of the radiuses in memory map
        memmap_P_radius_data = np.load('P_C1_radius_'+name+'_'+str(size)+'_data.npy', mmap_mode='r')
        memmap_P_radius_indices = np.load('P_C1_radius_'+name+'_'+str(size)+'_indices.npy', mmap_mode='r')
        memmap_P_radius_indptr = np.load('P_C1_radius_'+name+'_'+str(size)+'_indptr.npy', mmap_mode='r')

        P_radius_memmapped = csr_matrix((memmap_P_radius_data,memmap_P_radius_indices,memmap_P_radius_indptr),(len(self),len(self)))
        
        #loading the matrix of the |P|+radiuses in memory map
        memmap_P_abs_plus_radius_data = np.load('P_C1_abs_plus_radius_'+name+'_'+str(size)+'_data.npy', mmap_mode='r')
        memmap_P_abs_plus_radius_indices = np.load('P_C1_abs_plus_radius_'+name+'_'+str(size)+'_indices.npy', mmap_mode='r')
        memmap_P_abs_plus_radius_indptr = np.load('P_C1_abs_plus_radius_'+name+'_'+str(size)+'_indptr.npy', mmap_mode='r')

        P_abs_plus_radius_memmapped = csr_matrix((memmap_P_abs_plus_radius_data,memmap_P_abs_plus_radius_indices,memmap_P_abs_plus_radius_indptr),(len(self),len(self)))
        
        return P_memmapped, P_abs_plus_radius_memmapped, P_radius_memmapped
    
    def save_Q(self, Q, name='default'):
        """
        This function saves the output of decay_rump
        """
        size=len(self.partition)-1
        np.save('Q_C1_'+name+'_'+str(size)+'.npy',Q)
    
    def load_Q(self, name='default'):
        """
        This function loads the output of decay_rump
        """
        size=len(self.partition)-1
        Q = np.load('Q_C1_'+name+'_'+str(size)+'.npy', mmap_mode='r')
        return Q
    
    def save_vector_with_radius(self, v, v_radius, name_of_vector, name='default'):
        """
        This function saves the output of decay_rump
        """
        size=len(self.partition)-1
        np.save(name_of_vector+'_C1_'+name+'_'+str(size)+'.npy',v)
        np.save(name_of_vector+'_radius_C1_'+name+'_'+str(size)+'.npy',v_radius)
    
        
    def load_vector_with_radius(self, name_of_vector, name='default'):
        """
        This function saves the output of decay_rump
        """
        size=len(self.partition)-1
        v = np.load(name_of_vector+'_C1_'+name+'_'+str(size)+'.npy', mmap_mode='r')
        v_radius = np.load(name_of_vector+'_radius_C1_'+name+'_'+str(size)+'.npy', mmap_mode='r')
        
        return v, v_radius
    
    def bound_C0_norm(self,v):
        """
        Estimates the C0 norm of v
        """
        n=len(self.partition)
        size=n-1
        phi_part=v[0:n]
        nu_part=v[n:2*n]
        kappa_part=v[2*n]
        
        max_phi_part=max(abs(phi_part))
        max_nu_part=(max(abs(nu_part))*5)/(16*size)
        max_kappa_part = abs(kappa_part)*1.5
        
        return max_phi_part+max_nu_part+max_kappa_part
        
    def bound_weak_norm(self,v):
        """
        Estimates the weak norm of v
        """
        n=len(self.partition)
        size=n-1
        
        phi_part=v[0:n]
        nu_part=v[n:2*n]
        kappa_part = v[2*n]
        
        w = np.diff(phi_part)*15*size/8
        
        max_der_phi_part = max(abs(w))
        max_der_nu_part = max(abs(nu_part))
        max_der_kappa_part = abs(v[2*n])*6
        
        return max_der_phi_part+max_der_nu_part+max_der_kappa_part+self.bound_C0_norm(v)
        
    def bound_strong_norm(self, v):
        """
        Bounds the strong norm of a vector written in the basis
        """
        raise NotImplementedError
        
    def project_function_on_basis(self, f, f_prime, value_of_integral, output_rate = None):
        """
        This function projects a function f on the basis
        """
        n=len(self.partition)
        v=np.zeros(len(self))
        v_radius=np.zeros(len(self))
        
        for i in range(n):
            val = f(RIF(self.partition[i]))
            val_prime = f_prime(RIF(self.partition[i]))
            
            v[i]= val.center()
            v_radius[i]= val.absolute_diameter()/2
            
            v[n+i]=val_prime.center()
            v_radius[n+i]= val_prime.absolute_diameter()/2
            
            if output_rate:
                if (i%output_rate==0):
                    print i
        
        integral, integral_radius = rump_dot(self._integrals_double[0:2*n],self._integrals_double_radius[0:2*n],v[0:2*n],v_radius[0:2*n])
        
        v[2*n]=value_of_integral-integral
        v_radius[2*n]=integral_radius
        
        return v, v_radius
        
    
    def norms_from_Q(self, Q, N):
        """
        Takes the result of decay_rump and returns bounds for the norms
        """
        Norms=np.zeros(N+1)

        for i in range(N):
            Norms[i+1]=self.bound_weak_norm(Q[:,i])
        
        Norms[0]=1
        
        return Norms
    
    def bound_on_norms_of_powers_C0(self, dynamic, project_left = False, project_right = False):
        """
        Uniform bound on :math:`\|L^i\|`.
        
        Gives a bound to the norms of :math:`\|(\Pi_1 L \Pi_2)^i\|` for each `i`, where
        :math:`\Pi_1` is the identity if `project_left==False`
        and the discretization projection if otherwise, and similarly for :math:`\Pi_2` and `project_right`.

        We do not use the terms "before" and "after" because they are ambiguous: would they refer to the order
        of function evaluations, or to the order in writing?

        Args:
            project_left, project_right (bools): tell whether to discretize.

        Returns:
            M (real with RNDU): such that :math:`\|(\Pi_1 L \Pi_2)^i\| \leq M` for each i
        """
        lam=1/dynamic.expansivity()
        B=dynamic.distorsion()
        
        if (project_left == False) and (project_right==False):
            return (1+B/(1-lam)).upper()
        if (project_left == True) and (project_right==False):
            raise NotImplementedError
        if (project_left == False) and (project_right==True):
            raise NotImplementedError
        if (project_left == True) and (project_right==True):
            raise NotImplementedError
    
    def dfly_weak(self, dynamic, discrete=False):
        """
        Return constants of the dfly inequality :math:`\|L^n f\|_s \leq A \lambda^n \|f\|_s + B\|f\|_w` 

        This function should make `lasota_yorke_constants` and `iterate_with_lasota_yorke` obsolete.
        
        Input:
            Dynamic:
            discrete: if True, returns constants for the projected operator instead (TODO: experimental)

        Returns:
            (A, B, lambda): tuple of intervals
        """
        n=len(self.partition)-1
        
        M=self.bound_on_norms_of_powers_C0(dynamic, project_left = False, project_right = False)
        lam=1/dynamic.expansivity()
        B=dynamic.distorsion()
        
        if discrete:
            P_s= RUP(23)/8+RUP(45)/(4*n)
            A=2*P_s
            theta=P_s*(M+B)*lam
            if theta>=1:
                print "Not enough expansion"
                raise NotImplementedError
            return (A, P_s*B*(B+1)/(1-theta)+1, theta)
        return (M, M**2, lam)
    
    def bound_on_norms_of_powers(self, dynamic, project_left = False, project_right = False):
        """
        Uniform bound on :math:`\|L^i\|`.
        
        Gives a bound to the norms of :math:`\|(\Pi_1 L \Pi_2)^i\|` for each `i`, where
        :math:`\Pi_1` is the identity if `project_left==False`
        and the discretization projection if otherwise, and similarly for :math:`\Pi_2` and `project_right`.

        We do not use the terms "before" and "after" because they are ambiguous: would they refer to the order
        of function evaluations, or to the order in writing?

        Args:
            project_left, project_right (bools): tell whether to discretize.

        Returns:
            M (real with RNDU): such that :math:`\|(\Pi_1 L \Pi_2)^i\| \leq M` for each i
        """
        
        if (project_left == False) and (project_right==False):
            (A, B, lam) = self.dfly_weak(dynamic, discrete=False)
            return A*lam+B
        if (project_left == True) and (project_right==False):
            raise NotImplementedError
        if (project_left == False) and (project_right==True):
            raise NotImplementedError
        if (project_left == True) and (project_right==True):
            (A, B, lam) = self.dfly_weak(dynamic, discrete=True)
            return A*lam+B
    
    def dfly(self, dynamic, discrete=False):
        """
        Return constants of the dfly inequality :math:`\|L^n f\|_s \leq A \lambda^n \|f\|_s + B\|f\|_w` 

        This function should make `lasota_yorke_constants` and `iterate_with_lasota_yorke` obsolete.
        
        Input:
            Dynamic:
            discrete: if True, returns constants for the projected operator instead (TODO: experimental)

        Returns:
            (A, B, lambda): tuple of intervals
        """
        n=len(self.partition)-1
        M=self.bound_on_norms_of_powers_C0(dynamic, project_left = False, project_right = False)
        lam=1/dynamic.expansivity()
        B=dynamic.distorsion()
        Z=1/(1-lam**2)*(dynamic.third_distorsion().upper()+3*lam/(1-lam)*dynamic.distorsion().upper())
        B_C1=max(3*(lam)*B*M/(1-lam),3*M*(B/1-lam)**2+M*Z)+M*(lam**(1))*(1-lam**(1)) + M**2
        
        if discrete:
            A_2=self.projection_continuity_constants_with_respect_to_strong_norm()
            A_1, A_3 = self.projection_continuity_constants_with_respect_to_weak_norm(norm='strong')
            tilde_M=self.bound_on_norms_of_powers(dynamic, project_left = True, project_right = True)
            
            theta=A_2*M*lam**2+A_3/n #please note that in the function that returns A_3 we already divide once by eta
            
            bound = (A_2*B_C1*A_1*tilde_M)/(1-theta)
            
            if theta>=1:
                print "Not enough expansion"
                raise NotImplementedError
            return (A_2, bound, theta)
        return (M, B_C1 , lam**2)
    
    def compute_error(self,dynamic, norms, N_approx, P, v, n_iter=5):
        """
        This function bounds the error on the approximation of the invariant measure using
        a bootstrap argument; the number of the iterations of the bootstrap is 
        controlled by n_iter
        
        input:
        dynamic:
        norms: the bounds for the norms of the discretized operator
        v: a candidate for the approximation of the invariant measure
        """
        if len(norms)<(N_approx+1):
            norms_extended = extend_power_norms(norms, N_approx)
            norms = norms_extended
        
        residual = self.bound_weak_norm(self.mat_vec(P,v)-v)
        
        C_N=norms[N_approx-1]
        C, D = self.approximation_coefficients_from_norms(dynamic, norms, N_approx)
        #start the bootstrap process
        M = self.bound_on_norms_of_powers_C0(dynamic, project_left = False, project_right = False)
        (A_w, B_w, lam_W) =self.dfly_weak(dynamic, discrete=False )
        
        # Bound for the C^1 norm of the fixed point
        f_w=B_w*M
        
        # Bound for the C^2 norm of the fixed point
        (A_s, B_s, lam_s) =self.dfly(dynamic, discrete=False )
        f_s=B_s*f_w
        
        # first bound for the C1 error
        err=1/(1-C_N)*(C*f_s+D*f_w+residual)
        norm_numerical=self.bound_weak_norm(v)
        
        for i in range(n_iter):
            f_w=err+norm_numerical

            err=1/(1-C_N)*(C*f_s+D*f_w+residual)

        return err
    
    def bound_on_norms_of_powers_strong(self, dynamic, project_left=False, project_right = False):
        """
        This function uses the Lasota-Yorke to bound the strong norm of the operator and
        of the discretized operator
        """
        raise NotImplementedError
    
    def projection_continuity_constants_with_respect_to_weak_norm(self, norm ='strong'):
        """
        This function returns constant P_w,Q_w such that ||\Pi f||_w\leq P_w ||f||_w+Q_w\cdot \delta ||f||_s,
        if the variable norm is 'strong', and the constant P such that ||\Pi f||_w\leq P ||f||_w
        if norm is 'weak'. Default is norm='strong'
        """
        n=len(self.partition)-1
        if (norm=='weak'):
            return RUP(23)/8+(RUP(45)/4+RUP(32)/81+RUP(45)/16)/n
        if (norm=='strong'):
            P_w = RUP(23)/8+RUP(32)/(81*n)
            Q_w = (RUP(9)/4+RUP(9))/n
            return (P_w, Q_w)

    def projection_continuity_constants_with_respect_to_strong_norm(self):
        """
        This function returns constant P_w,Q_w such that ||\Pi f||_s\leq P_s ||f||_s
        """
        n=len(self.partition)-1
        return 4+RUP(32)/(81*n)+RUP(18)/(n**2)  
    
    def projection_coefficients(self,norm):
        """
        This function returns the bound coefficients such that
        ||\Pi f||_w \leq P_w ||f||_w+\delta Q_w ||f||_s, ||\Pi f||_s \leq P_s ||f||_s
        
        Args:
            norm: it can have value 'strong' or 'weak', if it is 'strong' it will return P_s
            if it is 'weak' it will return P_w and Q_w
        """
        if (norm=='strong'):
            return self.projection_continuity_constants_with_respect_to_strong_norm(norm=norm)
        if (norm=='weak'):
            return self.projection_continuity_constants_with_respect_to_weak_norm()

    def projection_error(self):
        """
        This function returns the coefficient Kdelta such that
        ||\Pi f-f||_w \leq Kdelta ||f||_s
        
        """
        n=len(self.partition)-1
        return RUP(5)/(2*n)
    
    def candidate(self,P,interval_matrix=True,n_iter=10):
        """ 
        This function computes a candidate for the approximation of the invariant 
        measure using the power method
        P: discretized matrix
        scipy_matrix: if the matrix is already in double 
        """
        if interval_matrix:
            PP=sage_sparse_to_scipy_sparse(P)
        else:
            PP=P

        v=np.zeros(len(self))
        v[len(self)-1]=1
        
        for i in range(n_iter):
            v=self.mat_vec(PP,v)
        
        return v
    
    def approximation_coefficients_from_norms(self, dynamic, norms, N):
        """
        This function computes the approximation coefficients for
        the N-th iteration, given the C_i
        """
        if len(norms)<(N+1):
            norms_extended = extend_power_norms(norms, N)
            norms = norms_extended
        
        (A,B,lam) = self.dfly(dynamic,discrete=False)
        M = self.bound_on_norms_of_powers(dynamic, project_left = False, project_right = False)
        Kdelta=self.projection_error()
        P= self.projection_continuity_constants_with_respect_to_weak_norm(norm='weak')
        
        C = Kdelta*sum(norms[i]*(lam**(N-i-1)) for i in range(N))*(A*lam+P*M)
        D = Kdelta*sum(norms[i] for i in range(N))*(A*lam+P*M+M)*B
        
        return C, D
    
    def parallel_decay_rump_argument(self, P, P_abs_plus_radius, P_radius, N, k):
        """
        This routine generates the input of the parallel decay rump on k processes
        """
        size = len(self.partition)-1
        step = len(self)//k
        rest = len(self)%k

        input_data_parallel=[(self, P, P_abs_plus_radius, P_radius, N, size*8/256, i*step,(i+1)*step) for i in range(0,k)]
        # the last process takes care of the last vectors
        input_data_parallel[k-1]=(self, P, P_abs_plus_radius, P_radius, N, size*8/256, (k-1)*step, len(self)-1)

        return input_data_parallel