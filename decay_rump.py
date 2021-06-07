"""
Generic functions to estimate decay times, using Rump99
"""

from sage.all import matrix,RealField,VectorSpace,MatrixSpace,RIF
#from joblib import Parallel, delayed
from sparse import *
from sparse import matrix_norminf,interval_norminf_error
from time_profile import *
from rump_implementation import *

from sage.parallel.decorate import parallel

import numpy as np
import ctypes

from sage.all import RIF,load

#load('kahan_with_radius.spyx')

FE_TONEAREST = 0x0000
FE_DOWNWARD = 0x0400
FE_UPWARD = 0x0800
FE_TOWARDZERO = 0x0c00
libc = ctypes.CDLL('libm.so.6')

def decay_rump_manager(basis, A, A_abs_plus_radius, A_radius, N, num_process=1, output_rate=1024, start=0, stop=None):
    if num_process>1:
        print "Parallel"
        input_data_parallel=basis.parallel_decay_rump_argument(A, A_abs_plus_radius, A_radius, N, num_process)
        Q_tot = np.zeros([len(basis), N])
        for data, output in decay_rump_parallel(input_data_parallel):
            old_rnd=libc.fegetround()
            libc.fesetround(FE_UPWARD)
            Q_tot=Q_tot + output
            libc.fesetround(old_rnd)
    else:
        print "Not parallel"
        Q_tot = decay_rump(basis, A, A_abs_plus_radius, A_radius, N, output_rate=1024, start=0, stop=None)
    return Q_tot

def decay_rump(basis, A, A_abs_plus_radius, A_radius, N, output_rate=1024, start=0, stop=None):
    print "start", start
    print "stop", stop
    
    # initializes the matrix storing the sum of |Pv|, for v in contracting pairs
    Q = np.zeros([len(basis), N])
    n=len(basis)
    
    #estimation vector by vector
    count=0

    for v, v_radius in basis.contracting_pairs_with_radius(start=start,stop=stop):
        assert len(v)==len(v_radius)
        
        
        count=count+1
        for i in range(N):
            
            v,v_radius = basis.mat_vec_with_radius(A, A_abs_plus_radius, A_radius, v, v_radius, debug=False)
            
            old_rnd=libc.fegetround()
            
            libc.fesetround(FE_UPWARD)
            Q[:,i]+=np.absolute(v)+v_radius
            libc.fesetround(old_rnd)
        if (count%output_rate)==0: print count
    
    return Q

@parallel
def decay_rump_parallel(basis, A, A_abs_plus_radius, A_radius, N, output_rate=1024, start=0, stop=None):
    print "start", start
    print "stop", stop
    
    # initializes the matrix storing the sum of |Pv|, for v in contracting pairs
    Q = np.zeros([len(basis), N])
    n=len(basis)
    
    #estimation vector by vector
    count=0

    for v, v_radius in basis.contracting_pairs_with_radius(start=start,stop=stop):
        assert len(v)==len(v_radius)
        
        
        count=count+1
        for i in range(N):
            
            v,v_radius = basis.mat_vec_with_radius(A, A_abs_plus_radius, A_radius, v, v_radius, debug=False)
            
            old_rnd=libc.fegetround()
            
            libc.fesetround(FE_UPWARD)
            Q[:,i]+=np.absolute(v)+v_radius
            libc.fesetround(old_rnd)
        if (count%output_rate)==0: print count
    
    return Q
