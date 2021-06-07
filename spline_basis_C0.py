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

class Spline_Basis_C0(Basis):
    """
    Abstract class for the projection from C^1 to C^0

    :math:`\phi_j` is the function so that :math:`\phi_j(x_i)=\delta_{ij}`, so :math:`\phi_0` has support :math:`[x_{n-1},1] \cup [0,x_1]`.

    """
    
    def __init__(self, partition):
        check_partition(partition)
        self.partition = partition

        n=len(self)

        self._integrals_double=np.zeros(n)
        self._integrals_double_radius=np.zeros(n)
        
        for i in range(n):
            self._integrals_double[i]=self.evaluate_integral(i).center()
            self._integrals_double_radius[i]=self.evaluate_integral(i).absolute_diameter()/2
        

    def __len__(self):
        """
        Returns the number of elements of the basis
        """
        return len(self.partition)+1

    def evaluate_phi(self,i,partsize,x):
        """
        Evaluates an element phi
        Args:
            i (integer) : the index of the element of the basis
            partsize : the size of the partition (please note that since we need to build a finer partition of unity to define tau,
            this argument is needed)
            x (real number) : the point on which we evaluate the function 
        """
    
        Interval_Field=x.parent()
        tilde_x=x*partsize-i
        #if x is not in [-1,1], return 0
        try:
            tilde_x=tilde_x.intersection(Interval_Field(-1,1))
        except ValueError:
            return Interval_Field(0)
        try:
            tilde_x_neg=tilde_x.intersection(Interval_Field(-1,0))
            val_neg = polyval([1,0,-3,-2], tilde_x_neg)
        except ValueError:
            val_neg=None
        try:
            tilde_x_pos=tilde_x.intersection(Interval_Field(0,1))
            val_pos = polyval([1,0,-3,2], tilde_x_pos)
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
    
    def evaluate(self,i,x):
        """
        A function that evaluates the i-th element of the basis on the point x
            Args:
                i : index of the element
                x : point on which we evaluate
        """
        n = len( self )
        partsize = n - 2
        if (i==partsize+1):
            return self.evaluate_tau(x)
        else:
            return self.evaluate_phi(i,partsize,x)
    
    def nonzero_on(self, I):
        """
        Returns which elements of the basis are nonzero on the interval I
            Args:
                I (interval) :
        """
        n = len( self )
        partsize = n - 2
        
        jmin = binsearch2(I.lower(), self.partition)
        jmax = binsearch2(I.upper(), self.partition)+1
        
        jmin = max ( 0, jmin )
        jmax = min (jmax,partsize)
        
        # please note that in this case, the index may be equal to partsize,
        # we are not on the torus
        
        for i in range(jmin, jmax+1):
            yield i

    def evaluate_integral(self, j, ring=RIF):
        """
        Returns the integral of the index j element of the basis
            Args:
                j : index of the element
                ring : the ring in which we return the value of the integral
        """
        n = len( self )
        partsize = n - 2
        if (j==partsize+1):
            return ring(1)
        elif (j==0) or (j==partsize):
            return ring(1)/(2*partsize)
        else:
            return ring(1)/partsize
        
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
    
    def evaluate_function(self,w,x):
        """
        Alias for evaluate_linear_combination
        """
        return self.evaluate_linear_combination(w,x)
    
    def plot_function(self,w):
        """plots the function whose coefficients are given by w"""
        pl=plot(lambda x: self.evaluate_function(w,RIF(x)).center(),0,1)
        return pl
    
    def mat_vec(self,PP,v):
        """
        The matrix-vector product (on the right) specialized for this basis
        Args:
            PP : the matrix
            v : the vector
        """
        
        n = len(self)
        
        # Matrix multiplication
        v=PP*v
        
        v[n-1]-=v[0:n-1].dot(self._integrals_double[0:n-1])

        return v
    
    def mat_vec_with_radius(self, PP, PP_abs_plus_radius, PP_radius,w, w_radius, known_integral = None, debug=False):
        """ 
        Matrix vector product associated to the basis
        """
        
        n = len(self)
        v = w.copy()
        v_radius = w_radius.copy()
        
        v, v_radius = rump_sparse_mat_vec(PP, PP_abs_plus_radius, PP_radius, v, v_radius)
        
        # we compute the error in the integral of the image
        integral, integral_radius = rump_dot(self._integrals_double[0:n-1],self._integrals_double_radius[0:n-1],v[0:n-1],v_radius[0:n-1])
        
        if known_integral is not None:
            v[n-1]=known_integral-integral
            v_radius[n-1]=integral_radius
        else:
            v[n-1]-=integral
            v_radius[n-1]+=integral_radius
        
        return v, v_radius
    
    def dual_composed_with_dynamic(self, dynamic, epsilon):
        """
        For each `i`, yields (successively) the sequence :math:`\{(y,T'(y)) \colon y\in T^{-1}(x_i)\}`.
        Please note that, for the last element (tau) it yelds (len(basis),1)
        """
        for i in range(len(self) - 1):
            x=self.partition[i]
            yield (i,[ (y,dynamic.fprime(y)) for y in dynamic.preimages(x, epsilon) ])
        yield (len(self)-1,[(len(self),1)]) 
        
    def project_dual_element(self, dual_element):
        """
        Yields the coefficients of the matrix.
        Depending on the dual element returns different things.
        If dual element is of the form :math:`\{(y,T'(y)) \colon y\in T^{-1}(x_i)\}` it computes
        the value of :math:'\phi' and :math:'\tau' through :math:'L' 
        If the dual element is of the form  (len(basis),1) it returns the integral of the element j
        """
        for y, Tprimey in dual_element:
            if (y==len(self)) and (Tprimey==1):
                for j in range(len(self)):
                    yield (j,self.evaluate_integral(j))
            else:
                for j in self.nonzero_on(y):
                    yield (j, self.evaluate(j,y) / Tprimey.abs())
                yield(len(self)-1, self.evaluate(len(self)-1,y) / Tprimey.abs())

    def contracting_pairs_with_radius(self,vectorspace=None, start=0, stop=None):
        """
        This function returns the basis of the space V_0.
        It returns a couple of vectors:
        v contains the midpoint of the element of the basis
        v_radius contains the radiuses of the elements of v, i.e.
        the real vector of the basis \tilde{v} is such that |(\tilde{v})_j-(v)_j|\leq (v_radius)_j
        where j is the j-th component
        """
        n=len(self.partition)
        if stop==None:
            stop=len(self)-1
            
        for i in range(start, stop):
            v=np.zeros(len(self))
            v_radius=np.zeros(len(self))

            #we have to normalize, we bound the C0 norm from below
            norm_lower_bound=abs(RIF(1)-self.evaluate_integral(i)*RIF(1.5))
            val_phi=RIF(1)/norm_lower_bound
            v[i]=val_phi.center()
            v_radius[i]=val_phi.absolute_diameter()/2
            
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
        np.save('P_C0_'+name+'_'+str(size)+'_data.npy',P.data)
        np.save('P_C0_'+name+'_'+str(size)+'_indices.npy',P.indices)
        np.save('P_C0_'+name+'_'+str(size)+'_indptr.npy',P.indptr)

        np.save('P_C0_radius_'+name+'_'+str(size)+'_data.npy',P_radius.data)
        np.save('P_C0_radius_'+name+'_'+str(size)+'_indices.npy',P_radius.indices)
        np.save('P_C0_radius_'+name+'_'+str(size)+'_indptr.npy',P_radius.indptr)
        
        P_abs_plus_radius=abs(P)+P_radius
        
        np.save('P_C0_abs_plus_radius_'+name+'_'+str(size)+'_data.npy',P_abs_plus_radius.data)
        np.save('P_C0_abs_plus_radius_'+name+'_'+str(size)+'_indices.npy',P_abs_plus_radius.indices)
        np.save('P_C0_abs_plus_radius_'+name+'_'+str(size)+'_indptr.npy',P_abs_plus_radius.indptr)
    
    def load_csr_matrix_with_radius_memmapped(self, name='default'):
        """
        This function loads the csr matrix containing the values
        and the radius as a memmap
        """
        size=len(self.partition)-1
        memmap_P_data = np.load('P_C0_'+name+'_'+str(size)+'_data.npy', mmap_mode='r')
        memmap_P_indices = np.load('P_C0_'+name+'_'+str(size)+'_indices.npy', mmap_mode='r')
        memmap_P_indptr = np.load('P_C0_'+name+'_'+str(size)+'_indptr.npy', mmap_mode='r')

        P_memmapped = csr_matrix((memmap_P_data,memmap_P_indices,memmap_P_indptr),(len(self),len(self)))

        #loading the matrix of the radiuses in memory map
        memmap_P_radius_data = np.load('P_C0_radius_'+name+'_'+str(size)+'_data.npy', mmap_mode='r')
        memmap_P_radius_indices = np.load('P_C0_radius_'+name+'_'+str(size)+'_indices.npy', mmap_mode='r')
        memmap_P_radius_indptr = np.load('P_C0_radius_'+name+'_'+str(size)+'_indptr.npy', mmap_mode='r')

        P_radius_memmapped = csr_matrix((memmap_P_radius_data,memmap_P_radius_indices,memmap_P_radius_indptr),(len(self),len(self)))
        
        #loading the matrix of the |P|+radiuses in memory map
        memmap_P_abs_plus_radius_data = np.load('P_C0_abs_plus_radius_'+name+'_'+str(size)+'_data.npy', mmap_mode='r')
        memmap_P_abs_plus_radius_indices = np.load('P_C0_abs_plus_radius_'+name+'_'+str(size)+'_indices.npy', mmap_mode='r')
        memmap_P_abs_plus_radius_indptr = np.load('P_C0_abs_plus_radius_'+name+'_'+str(size)+'_indptr.npy', mmap_mode='r')

        P_abs_plus_radius_memmapped = csr_matrix((memmap_P_abs_plus_radius_data,memmap_P_abs_plus_radius_indices,memmap_P_abs_plus_radius_indptr),(len(self),len(self)))
        
        return P_memmapped, P_abs_plus_radius_memmapped, P_radius_memmapped
    
    def save_Q(self, Q, name='default'):
        """
        This function saves the output of decay_rump
        """
        size=len(self.partition)-1
        np.save('Q_C0_'+name+'_'+str(size)+'.npy',Q)
    
    def load_Q(self, name='default'):
        """
        This function loads the output of decay_rump
        """
        size=len(self.partition)-1
        Q = np.load('Q_C0_'+name+'_'+str(size)+'.npy', mmap_mode='r')
        return Q
    
    def save_vector_with_radius(self, v, v_radius, name_of_vector, name='default'):
        """
        This function saves the output of decay_rump
        """
        size=len(self.partition)-1
        np.save(name_of_vector+'_C0_'+name+'_'+str(size)+'.npy',v)
        np.save(name_of_vector+'_radius_C0_'+name+'_'+str(size)+'.npy',v_radius)
    
        
    def load_vector_with_radius(self, name_of_vector, name='default'):
        """
        This function saves the output of decay_rump
        """
        size=len(self.partition)-1
        v = np.load(name_of_vector+'_C0_'+name+'_'+str(size)+'.npy', mmap_mode='r')
        v_radius = np.load(name_of_vector+'_radius_C0_'+name+'_'+str(size)+'.npy', mmap_mode='r')
        
        return v, v_radius
        
    
    def bound_weak_norm(self,v):
        """
        Estimates the weak norm of v
        """
        n=len(self)
        return max(abs(v[0:n-1]))+1.5*abs(v[n-1])
    
    def bound_strong_norm(self, v):
        """
        Bounds the strong norm of a vector written in the basis
        """
        n = len(self)
        size = n-2
        w = np.diff(v[0:n-1])*3*size/2
        
        max_der = max(abs(w))+6*abs(v[n-1])
        
        return max_der+self.bound_weak_norm(v)
    
    def project_function_on_basis(self, f, value_of_integral, output_rate = None):
        """
        This function projects a function f on the basis
        """
        n=len(self)
        v=np.zeros(n)
        v_radius=np.zeros(n)
        
        for i in range(len(self.partition)):
            val = f(RIF(self.partition[i]))
            v[i]= val.center()
            v_radius[i]= val.absolute_diameter()/2
            if output_rate:
                if (i%output_rate==0):
                    print i
        
        integral, integral_radius = rump_dot(self._integrals_double[0:n-1],self._integrals_double_radius[0:n-1],v[0:n-1],v_radius[0:n-1])
        
        v[n-1]=value_of_integral-integral
        v_radius[n-1]=integral_radius
        
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
        M=self.bound_on_norms_of_powers(dynamic, project_left = False, project_right = False)
        lam=1/dynamic.expansivity()
        B=dynamic.distorsion()
        
        if discrete:
            P_s = RUP(3)/2+RUP(6)/n
            A=2*P_s
            theta=P_s*(M+B)*lam
            if theta>=1:
                print "Not enough expansion"
                raise NotImplementedError
            return (A, P_s*B*(B+1)/(1-theta)+1, theta)
        return (M, M**2, lam)
    
    def bound_on_norms_of_powers_strong(self, dynamic, project_left=False, project_right = False):
        """
        This function uses the Lasota-Yorke to bound the strong norm of the operator and
        of the discretized operator
        """
        if (project_left==False) and (project_right == False):
            (A, B, lam) = self.dfly(dynamic , discrete=False)
            return (A*lam+B).upper()
        if (project_left==True) and (project_right == True):
            (A, B, lam) = self.dfly(dynamic , discrete = True)
            return (A*lam+B).upper()
        if (project_left==False) and (project_right == True):
            raise NotImplementedError
        if (project_left==True) and (project_right == False):
            raise NotImplementedError
    
    def projection_continuity_constants_with_respect_to_weak_norm(self, norm ='strong'):
        """
        This function returns constant P_w,Q_w such that ||\Pi f||_w\leq P_w ||f||_w+Q_w\cdot \delta ||f||_s,
        if the variable norm is 'strong', and the constant P such that ||\Pi f||_w\leq P ||f||_w
        if norm is 'weak'. Default is norm='strong'
        """
        n=len(self.partition)-1
        if (norm=='weak'):
            return RUP(4)
        if (norm=='strong'):
            P_w = RUP(1)
            Q_w = RUP(3)/(2)
            return (P_w, Q_w)

    def projection_continuity_constants_with_respect_to_strong_norm(self):
        """
        This function returns constant P_w,Q_w such that ||\Pi f||_s\leq P_s ||f||_s
        """
        n=len(self.partition)-1
        return RUP(3)/2+RUP(6)/n+RUP(3)/(2*n)
    
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
        
    def bound_difference_of_iteration_on_f_eta_radius(self, dynamic, P_C0, P_C0_abs_plus_radius, P_C0_radius,f_eta, f_eta_radius, l,debug=False):
        """
        This function bounds the error made computing \sum_{i=0}^{l-1}L^i f_{\eta}-\sum_{i=0}^{l-1}L_{\eta}^i f_{\eta}
        Inputs:
        dynamic:
        P_C0: the discretized operator
        f_eta: the approximation of \hat{L}h
        l: the first index of the tail
        """
        
        size=len(self)-2
        
        M = self.bound_on_norms_of_powers(dynamic, project_left = False, project_right = False)
        M_discretized_strong = self.bound_on_norms_of_powers_strong(dynamic, project_left=True, project_right = True)
        
        Kdelta=self.projection_error()
        P= self.projection_continuity_constants_with_respect_to_weak_norm(norm='weak')
        
        # Lasota-Yorke coefficients C0 space
        (A_C0, B_C0, lam_C0)=self.dfly(dynamic,discrete = False)
        
        difference_L_L_eta=Kdelta*(A_C0*lam_C0+P*M+B_C0)

        v=f_eta.copy()
        v_radius=f_eta_radius.copy()

        norms_L_i_f_eta = np.zeros(l)
        sums_L_i_f_eta = 0
        
        # the vector norms_L_i_f_eta contains a bound of the strong norm of f_eta, Q[i]=L_{\eta}^i f_{\eta}
        for i in range(l):
            norms_L_i_f_eta[i] = self.bound_strong_norm(v)+self.bound_strong_norm(v_radius)
            if i>0:
                if norms_L_i_f_eta[i]>norms_L_i_f_eta[i-1]:
                    norms_L_i_f_eta[i:l] = np.full(len(norms_L_i_f_eta[i:l]), min(M_discretized_strong * norms_L_i_f_eta[i-1], norms_L_i_f_eta[i]))
                    break
            v, v_radius = self.mat_vec_with_radius(P_C0,P_C0_abs_plus_radius,P_C0_radius,v,v_radius,known_integral=0) 
        
        # the vector sums_L_i_f_eta contains a bound of the sums of the strong norms up to i
        for i in range(l):
            if debug:
                print i
            sums_L_i_f_eta +=sum(norms_L_i_f_eta[0:i+1]) #M[0]=Q[0], M[1]=Q[0]+Q[1], M[2]=Q[0]+Q[1]+Q[2], ..., Q[l-1]=Q[0]+Q[1]+...+Q[l-1]
        
        return sums_L_i_f_eta*difference_L_L_eta*M, norms_L_i_f_eta
    
    def bound_sum_iterations_L_on_the_difference(self,dynamic,l,err_C0):
        """
        This function bounds the error made computing \sum_{i=0}^{l-1} L^i(f_{\eta}-\hat{L}h)
        Inputs:
        dynamic:
        l: the first index of the tail
        err_C1: the C0 error on f_{\eta}-\hat{L}h
        """
        M = self.bound_on_norms_of_powers(dynamic, project_left = False, project_right = False)
        return M*l*err_C0
    
    def bound_convergence_to_equilibrium(self, dynamic, N, C, D, norms):
        """
        This function returns the rho and the C_1, C_2 used in the estimating the size of the tail
        returns
        (rho, autonorm)
        """
        (A, B, lam_1) = self.dfly(dynamic,discrete = False)
        lam_2 = norms[N]
        
        return spectral_radius_and_autonorm_uniform(A, lam_1, B, N, lam_2, C, D)
        
    
    def estimate_error_tail_C0(self, dynamic, rho, autonorm, N, l, bound_norm):
        """
        This function estimates the size of the tail ||\sum_{i=l}^{+\infty}L^i \hat{L}h||_{\infty}
        dynamic:
        rho: output of the function spectral_radius_and_autonorm
        autonorm: output of the function spectral_radius_and_autonorm
        N: the N used in spectral_radius_and_autonorm
        l: the starting index for the tail (please note that it must be a multiple of N)
        bound_norm: a bound for the strong norm of \hat{L}h
        """
        M = self.bound_on_norms_of_powers(dynamic, project_left = False, project_right = False)
        
        assert ((l%N)==0)
        k=l//N
        
        C=1/autonorm[1]
        
        #This M is due to the fact that we can estimate L^k, to estimate L^{k+i} we use ||L^k||<=
        return bound_norm*(C*rho**k/(1-rho)+M*(N-1)*C*rho**k/(1-rho))