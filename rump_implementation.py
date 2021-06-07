import ctypes
import numpy as np

from sage.all import RIF,load

FE_TONEAREST = 0x0000
FE_DOWNWARD = 0x0400
FE_UPWARD = 0x0800
FE_TOWARDZERO = 0x0c00
libc = ctypes.CDLL('libm.so.6')


def rump_dot(A, A_radiuses, B, B_radiuses):
    """
    This is a Rump[99] dot product
    """
    old_rnd=libc.fegetround()
    
    libc.fesetround(FE_DOWNWARD)
    
    C_1 = np.dot(A,B)
    
    libc.fesetround(FE_UPWARD)
    C_2 = np.dot(A,B)
    
    C = C_1+0.5*(C_2-C_1)
    
    gamma = C-C_1+ np.dot(np.absolute(A)+A_radiuses,B_radiuses)+np.dot(A_radiuses,np.absolute(B))
    
    libc.fesetround(old_rnd)
    
    return C, gamma

def rump_dot_rad(A, A_radius, B, B_radius):
    """
    This is a modified version of Rump[99] dot product,
    with, instead of vectors of radiuses, upper bound on the radiuses
    """
    old_rnd=libc.fegetround()
    
    n=A.shape[1]
    
    libc.fesetround(FE_DOWNWARD)
    
    C_1 = np.dot(A,B)
    
    libc.fesetround(FE_UPWARD)
    C_2 = np.dot(A,B)
    
    C = C_1+0.5*(C_2-C_1)
    
    gamma = np.amax(C-C_1)+ n*(np.amax(np.fabs(A))+A_radius)*B_radius+n*A_radius*np.amax(np.fabs(B))
    
    libc.fesetround(old_rnd)
    
    return C, gamma

from scipy.sparse import *

def rump_sparse_mat_vec(A, A_abs_plus_radius, A_radius, v, v_radius):
    """
    This is a sparse matrix vector product, we assume A in CSR form, A_radius is a bound on the
    radius of the entries of A, v_radius is a bound on the radiuses of v
    input:
        A, the matrix
        A_abs_plus_radius, a matrix containing the absolute values of the elements of A plus a matrix containing the radius of the elements of A
        v, a vector
        v_radius, a vector containing the radiuses of the elements of v 
    """
    old_rnd=libc.fegetround()
    
    #nnz = max(np.diff(A.indptr))
    libc.fesetround(FE_DOWNWARD)
    
    C_1 = A.dot(v)
    
    libc.fesetround(FE_UPWARD)
    C_2 = A.dot(v)
    
    C = C_1+0.5*(C_2-C_1)
    
    #rump formula
    # c-c1+(|A|+rA)*rv+rA*v
    
    gamma = C-C_1+A_abs_plus_radius.dot(v_radius)+A_radius.dot(np.fabs(v))
    
    libc.fesetround(old_rnd)
    
    return C, gamma

def rump_sparse_mat_vec_unif_bound(A, norm_A, A_abs_plus_radius_bound, A_radius_bound, nnz, gamma_nnz, v,v_radius):
    """
    This is a modified version of Rump[99] dot product,
    with, instead of vectors of radiuses, upper bound on the radiuses
    """
    
    C = A.dot(v)
    max_abs_v=max(np.fabs(v))
    fl_error = norm_A*max_abs_v*gamma_nnz
    
    #each component of gamma in rump_sparse_mat_vec is bound by
    radius_bound = fl_error+nnz*A_abs_plus_radius_bound*max(v_radius)+nnz*A_radius_bound*max_abs_v
    
    radius=np.full(len(v_radius),radius_bound)
    return C, radius



def kahan_with_radius(v, v_radius, z):
    """ Kahan summation with control on the error, supposes that z is positive"""
    
    
    # This is the cycle for Kahan summation
    N=len(v)
    cumulator = np.zeros(N)
    compensation = 0.0
    cumulator[0] = v[0]*z[0]  
    for i in range(1, N):
        y = z[i]*v[i] - compensation 
        t = cumulator[i - 1] + y 
        compensation = (t - cumulator[i - 1]) - y 
        cumulator[i] = t
    
    # Now we compute the radius of the cumulative sums
    old_rnd=libc.fegetround()
    # We set the rounding upward, to take into account all the possible floating point errors
    libc.fesetround(FE_UPWARD)
    
    radius = np.zeros(N)
    
    abs_v = np.absolute(v)
    
    acc_radius = 0.0
    acc_abs = 0.0
    
    acc_radius = acc_radius+v_radius[0]
    acc_abs = acc_abs + abs_v[0]
    eps = np.finfo(float).eps 
    
    radius[0] = v_radius[0]
    for i in range(1, N):
        # we accumulate in this the sum with absolute value of the first i terms
        acc_abs = acc_abs + abs_v[i] * z[i]
        # This is the accumulated radius, please note the last term which is an upper bound on the numerical error of the midpoint, following Higham, Kahan, Knuth
        acc_radius = acc_radius + v_radius[i]*z[i] + 3 * eps * acc_abs
        radius[i] = acc_radius
    
    # back to the former rounding
    libc.fesetround(old_rnd)
    
    return cumulator, radius


