︠3e62df91-fca2-40a5-a1e0-07cd5177e879s︠
from spline_basis_C0 import *

from alternative_dynamic import *
name='Deterministic'
# this determines whether it should assemble the matrix or load it from a file
assemble = False

# this determines whether we have to estimate the decay time or load it from a file
decay = False

compute_f_eta = False

prec=53
N=20

dynamic=QuotientedLinearDynamic(k=2)

size=4194304

plot(lambda x: dynamic.f(RIF(x)).center(),0,1)

b_C0=Spline_Basis_C0(equispaced(size))
(A,B,lam)=b_C0.dfly(dynamic)
print "||L^n f||_{C^1}\leq ",A,"*",lam,"^n+",B,"||f||_{C^0}"
(A_unif,B_unif,lam_unif)=b_C0.dfly(dynamic, discrete=True)
print "||L_{\eta}^n f||_{C^1}\leq ",A_unif,"*",lam_unif,"^n+",B_unif,"||f||_{C^0}"

︡983e8ff2-23c2-4f44-ae7b-ef4c4d48632e︡{"file":{"filename":"/projects/f414be55-bd9f-499e-8370-d3ca3e1c83cd/.sage/temp/9aca8c334dfe/17011/tmp_1C7tgT.svg","show":true,"text":null,"uuid":"ee95bf70-b8e4-4585-b7f2-c031703b6bed"},"once":false}︡{"stdout":"||L^n f||_{C^1}\\leq  1.00000000000000 * 0.50000000000000000? ^n+ 1.00000000000000 ||f||_{C^0}\n"}︡{"stdout":"||L_{\\eta}^n f||_{C^1}\\leq  1.50000143051148 * 0.75000071525573731? ^n+ 1 ||f||_{C^0}\n"}︡{"done":true}︡
︠1fcb9df9-5b34-4c71-8676-2a7567fbc10f︠
if assemble:
    from generic_assembler import *
    P_C0, P_C0_radius=assemble(dynamic,b_C0,2**(-prec+8),output_rate=8*size/256,scipy_matrix=True, scipy_radius=True)
    b_C0.save_csr_matrix_with_radius(P_C0,P_C0_radius,name)
︡47e26c9d-5ead-44eb-a251-5a8dbf4c0cce︡{"done":true}︡
︠05fa90cb-82f1-4ec2-b0d9-cbdc0d0cf2ees︠
P_C0, P_C0_abs_plus_radius, P_C0_radius = b_C0.load_csr_matrix_with_radius_memmapped(name)
︡d016abbf-7e0b-4381-b941-fed6327c8617︡{"done":true}︡
︠cc10ea5f-7e71-4b24-81f9-44f710986aa7︠
f = lambda x : float(pi)/8*sin(2*float(pi)*x)+float(pi)/16*sin(4*float(pi)*x)
plot(f,0,1)
︡c3c77dfc-c521-41bb-b15b-39c434413e3e︡{"file":{"filename":"/projects/f414be55-bd9f-499e-8370-d3ca3e1c83cd/.sage/temp/9aca8c334dfe/16547/tmp_547dCs.svg","show":true,"text":null,"uuid":"bbe48d31-57a3-43c8-a0fb-0995c96dbedc"},"once":false}︡{"done":true}︡
︠44ad22f5-5d19-46a0-a203-abe388c93e9e︠
if compute_f_eta is True:
    f_eta, f_eta_radius=b_C0.project_function_on_basis(f, value_of_integral = 0)
    b_C0.save_vector_with_radius(f_eta, f_eta_radius,'f_eta_deterministic', name)
︡8329d7fc-9f9c-4334-8f6c-c913fb896be2︡{"done":true}︡
︠4ecc3756-c0d6-421f-bfb9-49bf356cb3a2s︠
f_eta, f_eta_radius = b_C0.load_vector_with_radius('f_eta_deterministic', name)
print b_C0.bound_strong_norm(f_eta)
︡ee1486f3-1b37-4f0b-9b77-79f9ac555cf6︡{"stdout":"7.91233437391"}︡{"stdout":"\n"}︡{"done":true}︡
︠f37a91d2-3a6b-46cd-b8c4-8dd4335ab342s︠
l = 57
err_diff_iter, norms_L_i_f_eta = b_C0.bound_difference_of_iteration_on_f_eta_radius(dynamic, P_C0, P_C0_abs_plus_radius, P_C0_radius,f_eta,f_eta_radius, l,debug=False)
print err_diff_iter
print norms_L_i_f_eta
︡7bd761e0-1819-420e-af05-9f6f33688643︡{"stdout":"0.001854277934518261?\n"}︡{"stdout":"[  7.91233438e+00   2.04690037e+00   6.41958290e-10   3.23169626e-10\n   3.30094555e-10   3.30094555e-10   3.30094555e-10   3.30094555e-10\n   3.30094555e-10   3.30094555e-10   3.30094555e-10   3.30094555e-10\n   3.30094555e-10   3.30094555e-10   3.30094555e-10   3.30094555e-10\n   3.30094555e-10   3.30094555e-10   3.30094555e-10   3.30094555e-10\n   3.30094555e-10   3.30094555e-10   3.30094555e-10   3.30094555e-10\n   3.30094555e-10   3.30094555e-10   3.30094555e-10   3.30094555e-10\n   3.30094555e-10   3.30094555e-10   3.30094555e-10   3.30094555e-10\n   3.30094555e-10   3.30094555e-10   3.30094555e-10   3.30094555e-10\n   3.30094555e-10   3.30094555e-10   3.30094555e-10   3.30094555e-10\n   3.30094555e-10   3.30094555e-10   3.30094555e-10   3.30094555e-10\n   3.30094555e-10   3.30094555e-10   3.30094555e-10   3.30094555e-10\n   3.30094555e-10   3.30094555e-10   3.30094555e-10   3.30094555e-10\n   3.30094555e-10   3.30094555e-10   3.30094555e-10   3.30094555e-10\n   3.30094555e-10]\n"}︡{"done":true}︡
︠ac1d2bc0-e7e8-4fbd-ab83-86a291999f46s︠
#This is the error due to the C0 approximation
err_approx_C1 = b_C0.projection_error()*(3.15**2)/2 #the derivative of \hat{L}h is bounded by pi**2/2<(3.15)**2/2
print err_approx_C1
err_it_on_diff = b_C0.bound_sum_iterations_L_on_the_difference(dynamic,l,err_approx_C1)
print err_it_on_diff
︡320997fc-fb66-4bd9-ac2b-a692bb42a191︡{"stdout":"2.95713543891907e-6\n"}︡{"stdout":"0.000168556720018387\n"}︡{"done":true}︡
︠9ae02fc7-7f3c-4e42-ae2b-40645c204528s︠
#The total error depends on the error on the tail, from the C0 part of the approximation
err_tail = 0.00054
total_error=err_diff_iter+err_it_on_diff+err_tail
print total_error
︡68e9c773-a8e0-4d1a-a1d8-86617ef3133c︡{"stdout":"0.002562834654536647?\n"}︡{"done":true}︡









