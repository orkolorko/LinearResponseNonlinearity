"""
Rewrite of the various assembler functions in a more generic fashion.
"""

from sage.all import matrix,RDF
from scipy.sparse import *
from huge_matrix import *


def assemble(dynamic, basis, epsilon, prec=53, scipy_matrix=False, scipy_radius=False, huge_matrix_flag=False, output_rate=1024):
    """
    Very generic assembler function
    """
    
    print "Epsilon ", epsilon
    
    absolute_diameter=0
    
    n = len(basis)
    
    if scipy_matrix:
        print "sparse matrix,double, csr,scipy"
        P = lil_matrix((n,n))
        if scipy_radius:
            print "sparse matrix, with radius, double, csr, scipy"
            P_radius = lil_matrix((n,n))
    elif huge_matrix_flag:
        print "huge block matrix"
        P = huge_matrix(n,n)
    else:
        print "interval matrix"
        P = matrix(dynamic.field, n, n, sparse=True)

    i = 0
    print "start"
    for i, dual_element in basis.dual_composed_with_dynamic(dynamic, epsilon):
        for j,x in basis.project_dual_element(dual_element):
            #if (x.absolute_diameter()>1024*epsilon):
            #	print i, j, x.absolute_diameter()
            if scipy_matrix:
                P[i,j] += RDF(x.center())
                if scipy_radius:
                    P_radius [i,j] +=x.absolute_diameter()/2
                absolute_diameter=max(absolute_diameter,x.absolute_diameter())
            elif huge_matrix_flag:
                P.sum_coefficients(i,j,RDF(x.center()))
                absolute_diameter=max(absolute_diameter,x.absolute_diameter())
            else:
                P[i,j] += x	

        if (output_rate is not 0) and (i%output_rate==0):
            print i

    if scipy_matrix:
        P=P.tocsr()
        if scipy_radius:
            P_radius=P_radius.tocsr()
            return P, P_radius
        print "absolute diameter", absolute_diameter
        return P,absolute_diameter

    if huge_matrix_flag:
        P.finished_assembling()
        print "absolute diameter", absolute_diameter
        return P,absolute_diameter

    return P
