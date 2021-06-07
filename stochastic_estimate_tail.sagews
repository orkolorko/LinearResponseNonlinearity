︠82fbc4c4-8240-4ca7-9fa2-87e80e24a19f︠
from spline_basis_C0 import *
from alternative_dynamic import *
name='Stochastic'

# this determines whether it should assemble the matrix or load it from a file
assemble_matrix = True

# this determines whether we have to estimate the decay time or load it from a file
compute_decay = True
num_process = 20

prec=53
N=15

dynamic=MorePerturbedLinearDynamic(k=8,c='0.0025')

size=65536*4

plot(lambda x: dynamic.f(RIF(x)).center(),0,1)

from generic_assembler import *

b_C0=Spline_Basis_C0(equispaced(size))
(A,B,lam)=b_C0.dfly(dynamic)
print "||L^n f||_{C^1}\leq ",A,"*",lam,"^n+",B,"||f||_{C^0}"
(A_unif,B_unif,lam_unif)=b_C0.dfly(dynamic, discrete=True)
print "||L_{\eta}^n f||_{C^1}\leq ",A_unif,"*",lam_unif,"^n+",B_unif,"||f||_{C^0}"
︡0d111185-11a5-4d47-9455-d3cbf2ecdd41︡{"stderr":"Compiling ./binsearch2.spyx...\nCompiling ./cython_sparse_mat_vec.spyx..."}︡{"stderr":"\nCompiling ./polyval_lib.spyx..."}︡{"stderr":"\n"}︡{"file":{"filename":"/projects/f414be55-bd9f-499e-8370-d3ca3e1c83cd/.sage/temp/9aca8c334dfe/23978/tmp_GREY00.svg","show":true,"text":null,"uuid":"674f1d7f-8ce5-4073-8966-1203be949795"},"once":false}︡{"stdout":"||L^n f||_{C^1}\\leq  1.19979547394928 * 0.1265744527098697? ^n+ 1.43950917930917 ||f||_{C^0}\n"}︡{"stdout":"||L_{\\eta}^n f||_{C^1}\\leq  3.00004577636719 * 0.260931256278564? ^n+ 1.415987231562125? ||f||_{C^0}\n"}︡{"done":true}︡
︠388aa4fb-25af-4ed5-8b8c-4aa174b3633c︠
if assemble_matrix is True:    
    #Assembling the dis cretized operator
    P_C0, P_C0_radius=assemble(dynamic,b_C0,2**(-prec+8),output_rate=8*size/256,scipy_matrix=True, scipy_radius=True)
    b_C0.save_csr_matrix_with_radius(P_C0,P_C0_radius,name)
︡84dd56b9-f362-4867-83df-0fa82e937da5︡{"stdout":"Epsilon  1/35184372088832\nsparse matrix,double, csr,scipy\nsparse matrix, with radius, double, csr, scipy"}︡{"stdout":"\nstart"}︡{"stdout":"\n0\n8192"}︡{"stdout":"\n16384"}︡{"stdout":"\n24576"}︡{"stdout":"\n32768"}︡{"stdout":"\n40960"}︡{"stdout":"\n49152"}︡{"stdout":"\n57344"}︡{"stdout":"\n65536"}︡{"stdout":"\n73728"}︡{"stdout":"\n81920"}︡{"stdout":"\n90112"}︡{"stdout":"\n98304"}︡{"stdout":"\n106496"}︡{"stdout":"\n114688"}︡{"stdout":"\n122880"}︡{"stdout":"\n131072"}︡{"stdout":"\n139264"}︡{"stdout":"\n147456"}︡{"stdout":"\n155648"}︡{"stdout":"\n163840"}︡{"stdout":"\n172032"}︡{"stdout":"\n180224"}︡{"stdout":"\n188416"}︡{"stdout":"\n196608"}︡{"stdout":"\n204800"}︡{"stdout":"\n212992"}︡{"stdout":"\n221184"}︡{"stdout":"\n229376"}︡{"stdout":"\n237568"}︡{"stdout":"\n245760"}︡{"stdout":"\n253952"}︡{"stdout":"\n114688"}︡{"stdout":"\n122880"}︡{"stdout":"\n131072"}︡{"stdout":"\n139264"}︡{"stdout":"\n147456"}︡{"stdout":"\n155648"}︡{"stdout":"\n163840"}︡{"stdout":"\n172032"}︡{"stdout":"\n180224"}︡{"stdout":"\n188416"}︡{"stdout":"\n196608"}︡{"stdout":"\n204800"}︡{"stdout":"\n212992"}︡{"stdout":"\n221184"}︡{"stdout":"\n229376"}︡{"stdout":"\n237568"}︡{"stdout":"\n245760"}︡{"stdout":"\n253952"}︡{"stdout":"\n262144"}︡{"stdout":"\n"}︡{"done":true}
︠262dafba-dcfb-4702-8191-c44c164ba271︠
P_C0, P_C0_abs_plus_radius, P_C0_radius = b_C0.load_csr_matrix_with_radius_memmapped(name)
︡f163d8a3-fdff-42a3-ade9-6b93722e2377︡{"done":true}︡
︠ecfef5a1-f211-4b3d-9c3a-83de1ba7146b︠
if compute_decay is True:
    from decay_rump import *
    Q_tot=decay_rump_manager(b_C0, P_C0, P_C0_abs_plus_radius, P_C0_radius, N, num_process=num_process, output_rate=1024, start=0, stop=None)
    b_C0.save_Q(Q_tot, name)
︡c256adf3-b003-48c6-b4a4-4db63cff9542︡{"stderr":"Error in lines 1-11\nTraceback (most recent call last):\n  File \"/usr/local/sage/local/lib/python2.7/site-packages/smc_sagews/sage_server.py\", line 995, in execute\n    exec compile(block+'\\n', '', 'single') in namespace, locals\n  File \"\", line 1, in <module>\nNameError: name 'decay' is not defined\n"}︡{"done":true}︡
︠85b0631f-ae5f-46eb-8de5-f24aa4e7bcf9s︠
Q_alt = b_C0.load_Q(name)

Norms = b_C0.norms_from_Q(Q_alt, N)
print Norms
︡1a93b255-acc1-40b3-8b8b-2cf6f6d3bdc0︡{"stdout":"[  1.00000000e+00   4.84662819e+00   4.88228359e+00   4.86549054e+00\n   4.73642896e+00   3.83080641e+00   2.86082069e-02   1.51304716e-02\n   1.20797333e-03   9.38640834e-05   3.50135275e-05   9.22571466e-05\n   2.76426887e-04   8.29257878e-04   2.48777259e-03   7.46331770e-03]\n"}︡{"done":true}︡
︠77a850ca-6960-4dc3-a332-5da0984be53ds︠
for N_approx in range(7,10):
    print "Norm of the ", N_approx, " iterate ", Norms[N_approx]
    C, D = b_C0.approximation_coefficients_from_norms(dynamic, Norms, N_approx)
    print "approximation coefficients ", C,"||f||_s", D,"||f||_w"
    rho, autonorm = b_C0.bound_convergence_to_equilibrium(dynamic, N_approx, C, D, Norms)
    print "rho", rho
    print "autonorm", autonorm
    print "C_s", 1/autonorm[0], "C_w", 1/autonorm[1]
    l=2*N_approx
    err_tail = b_C0.estimate_error_tail_C0(dynamic, rho, autonorm, N_approx, l, 5.6)
    print "l", l, "size of tail", err_tail
︡5788532b-8ffc-4fb1-bc25-ea64dc07f6f0︡{"stdout":"Norm of the  7  iterate  0.0151304715982\napproximation coefficients  0.0000283610106762083? ||f||_s 0.002042628073689518? ||f||_w\n[ 6.2449595102974?e-7 1.43950917930917 ]\n[ 0.0000283610106762083? 0.01717309967189476? ]\nrho 0.01928963794441759?\nautonorm [0.00146816073458602?, 0.998531839265414?]\nC_s 681.12433226324? C_w 1.001470319399796?\nl 14 size of tail 0.0174454570563118?\nNorm of the  8  iterate  0.00120797332816\napproximation coefficients  4.30419258087988?e-6 ||f||_s 0.002043905693032173? ||f||_w\n[ 7.9045233221118?e-8 1.43950917930917 ]\n[ 4.30419258087988?e-6 0.00325187902119108? ]\nrho 0.00459910457724841?\nautonorm [0.00093501722413716?, 0.999064982775863?]\nC_s 1069.49901476180? C_w 1.000935892299558?\nl 16 size of tail 0.00111944957350636?\nNorm of the  9  iterate  9.38640833672e-05\napproximation coefficients  6.0183751319626?e-7 ||f||_s 0.002044007694485297? ||f||_w\n[ 1.00051071342870?e-8 1.43950917930917 ]\n[ 6.0183751319626?e-7 0.002137871777852514? ]\nrho 0.00248632011544923?\nautonorm [0.00024200193887656?, 0.999757998061124?]\nC_s 4132.1982982546? C_w 1.000242060517992?\nl 18 size of tail 0.000367897798886193?\n"}︡{"done":true}︡









