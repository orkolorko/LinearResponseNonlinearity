︠1188fb60-538b-4877-8f30-d79c7244dc6ds︠
print "In this file we do the final estimates needed to compute the linear response in the stochastic perturbed case"
︡780d4c85-67a4-4ea0-a60e-40ccf9031a89︡{"stdout":"In this file we do the final estimates needed to compute the linear response in the stochastic perturbed case\n"}︡{"done":true}︡
︠e5e61c44-b26e-4205-8e1e-679ed357869as︠

from spline_basis_C0 import *
from alternative_dynamic import *
name='Stochastic'
# this determines whether it should assemble the matrix or load it from a file
assemble = True

# this determines whether we have to estimate the decay time or load it from a file
decay = False

prec=53
N=15

dynamic=MorePerturbedLinearDynamic(k=8,c='0.0025')

size=524288

from generic_assembler import *

b_C0=Spline_Basis_C0(equispaced(size))
(A,B,lam)=b_C0.dfly(dynamic)
print "||L^n f||_{C^1}\leq ",A,"*",lam,"^n+",B,"||f||_{C^0}"
(A_unif,B_unif,lam_unif)=b_C0.dfly(dynamic, discrete=True)
print "||L_{\eta}^n f||_{C^1}\leq ",A_unif,"*",lam_unif,"^n+",B_unif,"||f||_{C^0}"

plot(lambda x: dynamic.f(RIF(x)).center(),0,1)
︡1c197af7-e0d3-475c-9dc0-193b08a62a48︡{"stdout":"||L^n f||_{C^1}\\leq  1.19979547394928 * 0.1265744527098697? ^n+ 1.43950917930917 ||f||_{C^0}\n"}︡{"stdout":"||L_{\\eta}^n f||_{C^1}\\leq  1.50001144409180 * 0.260929265561440? ^n+ 1.415982937409840? ||f||_{C^0}\n"}︡{"file":{"filename":"/projects/f414be55-bd9f-499e-8370-d3ca3e1c83cd/.sage/temp/9aca8c334dfe/17425/tmp_1yYacq.svg","show":true,"text":null,"uuid":"78207619-b5f8-4a32-bc1a-23da4ecc1d86"},"once":false}︡{"done":true}︡
︠ba0afc93-6eef-4467-810e-c48f6061e0da︠
if assemble is True:
    from generic_assembler import *
    P_C0, P_C0_radius=assemble(dynamic,b_C0,2**(-prec+8),output_rate=8*size/256,scipy_matrix=True, scipy_radius=True)
    b_C0.save_csr_matrix_with_radius(P_C0,P_C0_radius,name)
︡79f8d97b-ff2b-45ff-8aba-ba67c07f68c2︡{"stdout":"Epsilon  1/35184372088832\nsparse matrix,double, csr,scipy\nsparse matrix, with radius, double, csr, scipy"}︡{"stdout":"\nstart"}︡{"stdout":"\n0\n16384"}︡{"stdout":"\n32768"}︡{"stdout":"\n49152"}︡{"stdout":"\n65536"}︡{"stdout":"\n81920"}︡{"stdout":"\n98304"}︡{"stdout":"\n114688"}︡{"stdout":"\n131072"}︡{"stdout":"\n147456"}︡{"stdout":"\n163840"}︡{"stdout":"\n180224"}︡{"stdout":"\n196608"}︡{"stdout":"\n212992"}︡{"stdout":"\n229376"}︡{"stdout":"\n245760"}︡{"stdout":"\n262144"}︡{"stdout":"\n278528"}︡{"stdout":"\n294912"}︡{"stdout":"\n311296"}︡{"stdout":"\n327680"}︡{"stdout":"\n344064"}︡{"stdout":"\n360448"}︡{"stdout":"\n376832"}︡{"stdout":"\n393216"}︡{"stdout":"\n409600"}︡{"stdout":"\n425984"}︡{"stdout":"\n442368"}︡{"stdout":"\n458752"}︡{"stdout":"\n475136"}︡{"stdout":"\n491520"}︡{"stdout":"\n491520"}︡{"stdout":"\n507904"}︡{"stdout":"\n524288"}︡{"stdout":"\n"}︡{"done":true}
︠69d480f4-2cf4-4d4b-8a5e-8f3a94810766s︠
P_C0, P_C0_abs_plus_radius, P_C0_radius = b_C0.load_csr_matrix_with_radius_memmapped(name)
︡1bf05e19-d007-4b7e-9478-c5201dfe47ea︡{"done":true}︡
︠8663e79e-a769-4bb8-a473-43044b76c8c1s︠
z = np.zeros(len(b_C0))
z[0:len(b_C0)-1] = np.load('derivative_'+str(size)+'.npy') 
int_error = -b_C0.integral(z)
print int_error
z[len(b_C0)-1]=int_error
f_eta = z
print b_C0.bound_strong_norm(f_eta)
︡2974093d-c4e7-487c-978f-5b10cc56a6bd︡{"stdout":"-1.71496427348e-14\n"}︡{"stdout":"2.71278129981\n"}︡{"done":true}︡
︠c2b2c0d2-1637-4e2a-a9fd-7dc74dd4c013s︠
pl=b_C0.plot_function(f_eta)
pl.show()
︡31faaf08-b127-4a35-9e3f-4b45283a6777︡{"file":{"filename":"/projects/f414be55-bd9f-499e-8370-d3ca3e1c83cd/.sage/temp/9aca8c334dfe/17425/tmp_KMr0Fw.svg","show":true,"text":null,"uuid":"8c334c19-7d1b-427a-b036-08c6be8aecbc"},"once":false}︡{"done":true}︡
︠9a7a79ef-81db-425c-bb72-a44fcf1c69dfs︠
f_eta_radius=np.zeros(len(b_C0))
︡a24f8385-68af-41db-b336-285afb0e0816︡{"done":true}︡
︠115b8918-50fd-4215-bf7a-de0316586c8es︠
l = 18
err_diff_iter, norms_L_i_f_eta = b_C0.bound_difference_of_iteration_on_f_eta_radius(dynamic, P_C0, P_C0_abs_plus_radius, P_C0_radius,f_eta,f_eta_radius, l,debug=False)
print err_diff_iter
for i in range(l):
    print i, norms_L_i_f_eta[i]
︡915735e4-9882-4524-b1f4-d1377790aa64︡{"stdout":"0.001785652696123013?\n"}︡{"stdout":"0 2.71278129981\n1 0.000606933204692\n2 1.43421475589e-05\n3 3.34073202251e-07\n4 1.06736166367e-08\n5 4.79289879883e-09\n6 8.66258857717e-09\n7 8.66258857717e-09\n8 8.66258857717e-09\n9 8.66258857717e-09\n10 8.66258857717e-09\n11 8.66258857717e-09\n12 8.66258857717e-09\n13 8.66258857717e-09\n14 8.66258857717e-09\n15 8.66258857717e-09\n16 8.66258857717e-09\n17 8.66258857717e-09\n"}︡{"done":true}︡
︠b0933fd1-66ff-447e-a133-25e1d8fe9e26s︠
#This is a bound on the error obtained in the C1 part of the approximation with stochastic perturbation
err_approx_C1 = 0.0002174
err_it_on_diff = b_C0.bound_sum_iterations_L_on_the_difference(dynamic,l,err_approx_C1)
print err_it_on_diff
︡ce90f7dd-4ece-4dcb-b7f5-0db0bb6f5160︡{"stdout":"0.00469503964865831\n"}︡{"done":true}︡
︠df3b5bef-9067-45b6-86b6-efb00471434cs︠
#The total error depends on the error on the tail, from the C0 part of the approximation
err_tail = 0.0003679
total_error=err_diff_iter+err_it_on_diff+err_tail
print total_error
︡1dab4f83-b58b-4c74-8e64-04c16d3d3bc1︡{"stdout":"0.006848592344781317?\n"}︡{"done":true}︡
︠1e2bbae3-34ce-442a-830c-63ada6107c28︠
v=b_C0.linear_response(P_C0,f_eta,l,interval_matrix=False)
pl=b_C0.plot_function(v)
pl.show()
︡19ccc398-4676-4597-b46a-5b8711fced27︡{"file":{"filename":"/projects/dc74f09d-fa5c-46f2-8b80-867694094c12/.sage/temp/compute8-us/24780/tmp_GPriFP.svg","show":true,"text":null,"text":null,"uuid":"96fd1783-c1bb-467d-b1ba-b2fe75c41dd9"},"once":false}︡{"done":true}︡︡{"done":true}︡









