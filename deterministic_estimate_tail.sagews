︠c42cd231-117e-4471-b0de-a26c9cf75305s︠
from spline_basis_C0 import *
from alternative_dynamic import *
name='Deterministic'
# this determines whether it should assemble the matrix or load it from a file
assemble_matrix = True

# this determines whether we have to estimate the decay time or load it from a file
compute_decay = True
num_process = 20

prec=53
N=20

dynamic=QuotientedLinearDynamic(k=2)

size=16384

plot(lambda x: dynamic.f(RIF(x)).center(),0,1)

from generic_assembler import *

b_C0=Spline_Basis_C0(equispaced(size))
(A,B,lam)=b_C0.dfly(dynamic)
print "||L^n f||_{C^1}\leq ",A,"*",lam,"^n+",B,"||f||_{C^0}"
(A_unif,B_unif,lam_unif)=b_C0.dfly(dynamic, discrete=True)
print "||L_{\eta}^n f||_{C^1}\leq ",A_unif,"*",lam_unif,"^n+",B_unif,"||f||_{C^0}"
︡6b6f2d6a-3cd5-4d61-aaa2-7ca120f8d8d5︡{"file":{"filename":"/projects/f414be55-bd9f-499e-8370-d3ca3e1c83cd/.sage/temp/9aca8c334dfe/1956/tmp_ctJG_f.svg","show":true,"text":null,"uuid":"a302dc5d-0f51-4b2d-ac43-3cd92fe173b2"},"once":false}︡{"stdout":"||L^n f||_{C^1}\\leq  1.00000000000000 * 0.50000000000000000? ^n+ 1.00000000000000 ||f||_{C^0}\n"}︡{"stdout":"||L_{\\eta}^n f||_{C^1}\\leq  3.00073242187500 * 0.75018310546875000? ^n+ 1 ||f||_{C^0}\n"}︡{"done":true}︡
︠388aa4fb-25af-4ed5-8b8c-4aa174b3633cs︠
if assemble_matrix is True:    
    #Assembling the dis cretized operator
    P_C0, P_C0_radius=assemble(dynamic,b_C0,2**(-prec+8),output_rate=8*size/256,scipy_matrix=True, scipy_radius=True)
    b_C0.save_csr_matrix_with_radius(P_C0,P_C0_radius,name)
︡160140e3-065a-487b-90d5-9616d533492e︡{"stdout":"Epsilon  1/35184372088832\nsparse matrix,double, csr,scipy\nsparse matrix, with radius, double, csr, scipy\nstart"}︡{"stdout":"\n0\n512"}︡{"stdout":"\n1024"}︡{"stdout":"\n1536"}︡{"stdout":"\n2048"}︡{"stdout":"\n2560"}︡{"stdout":"\n3072"}︡{"stdout":"\n3584"}︡{"stdout":"\n4096"}︡{"stdout":"\n4608"}︡{"stdout":"\n5120"}︡{"stdout":"\n5632"}︡{"stdout":"\n6144"}︡{"stdout":"\n6656"}︡{"stdout":"\n7168"}︡{"stdout":"\n7680"}︡{"stdout":"\n8192"}︡{"stdout":"\n8704"}︡{"stdout":"\n9216"}︡{"stdout":"\n9728"}︡{"stdout":"\n10240"}︡{"stdout":"\n10752"}︡{"stdout":"\n11264"}︡{"stdout":"\n11776"}︡{"stdout":"\n12288"}︡{"stdout":"\n12800"}︡{"stdout":"\n13312"}︡{"stdout":"\n13824"}︡{"stdout":"\n14336"}︡{"stdout":"\n14848"}︡{"stdout":"\n15360"}︡{"stdout":"\n15872"}︡{"stdout":"\n16384"}︡{"stdout":"\n"}︡{"done":true}︡
︠262dafba-dcfb-4702-8191-c44c164ba271s︠
P_C0, P_C0_abs_plus_radius, P_C0_radius = b_C0.load_csr_matrix_with_radius_memmapped(name)
︡88eb2307-e2de-42ce-bc18-5f7b356932f8︡{"done":true}︡
︠ecfef5a1-f211-4b3d-9c3a-83de1ba7146bs︠
if compute_decay is True:
    from decay_rump import *
    Q_tot=decay_rump_manager(b_C0, P_C0, P_C0_abs_plus_radius, P_C0_radius, N, num_process=num_process, output_rate=1024, start=0, stop=None)
    b_C0.save_Q(Q_tot, name)
︡75797689-e0d6-431d-a67a-cad5043fff4c︡{"stdout":"Parallel\n"}︡{"stdout":"start 15561\nstop 16385\n512\n"}︡{"stdout":"start 14742\nstop 15561\n512\n"}︡{"stdout":"start 13104\nstop 13923\n512\n"}︡{"stdout":"start 10647\nstop 11466\n512\n"}︡{"stdout":"start 1638\nstop 2457\n512\n"}︡{"stdout":"start 8190\nstop 9009\n512\n"}︡{"stdout":"start 5733\nstop 6552\n512\nstart 9828\nstop 10647\n512\n"}︡{"stdout":"start 12285\nstop 13104\n512\nstart 13923\nstop 14742\n512\nstart 3276\nstop 4095\n512\n"}︡{"stdout":"start 11466\nstop 12285\n512\nstart 819\nstop 1638\n512\nstart 7371\nstop 8190\n512\n"}︡{"stdout":"start 4914\nstop 5733\n512\nstart 9009\nstop 9828\n512\n"}︡{"stdout":"start 6552\nstop 7371\n512\nstart 0\nstop 819\n512\nstart 2457\nstop 3276\n512\n"}︡{"stdout":"start 4095\nstop 4914\n512\n"}︡{"done":true}︡
︠313fd237-ffdb-41b6-8a0b-bcc7a127601cs︠
Q_alt = b_C0.load_Q(name)

Norms = b_C0.norms_from_Q(Q_alt, N)
print Norms
︡51322267-750d-4372-845e-2522ca77fdd6︡{"stdout":"[  1.00000000e+00   2.12496566e+00   2.03093240e+00   2.00701205e+00\n   2.00017931e+00   1.99676293e+00   1.99249103e+00   1.98458673e+00\n   1.96894582e+00   1.93770426e+00   1.87520133e+00   1.75019058e+00\n   1.50016790e+00   1.00012232e+00   3.13059698e-05   2.35260089e-06\n   7.05780234e-06   2.11734070e-05   6.35202209e-05   1.90560662e-04\n   5.71681987e-04]\n"}︡{"done":true}︡
︠77a850ca-6960-4dc3-a332-5da0984be53ds︠
bound_norm= float(3*pi/16+pi**2/2) #bounds the norm of \hat{L}h
for N_approx in range(17,20):
    print "Norm of the ", N_approx, " iterate ", Norms[N_approx]
    C, D = b_C0.approximation_coefficients_from_norms(dynamic, Norms, N_approx)
    print "approximation coefficients ", C,"||f||_s", D,"||f||_w"
    rho, autonorm = b_C0.bound_convergence_to_equilibrium(dynamic, N_approx, C, D, Norms)
    print "rho", rho
    print "autonorm", autonorm
    print "C_s", 1/autonorm[0], "C_w", 1/autonorm[1]
    l=3*N_approx
    err_tail = b_C0.estimate_error_tail_C0(dynamic, rho, autonorm, N_approx, l, 5.6)
    print "l", l, "size of tail", err_tail
︡5c6112ba-6a15-4482-8e6f-72d2e1bb6228︡{"stdout":"Norm of the  17  iterate  2.11734069764e-05\napproximation coefficients  0.0002288926879375552? ||f||_s 0.02112291973300098? ||f||_w\n[ 7.6293945312500000?e-6 1.00000000000000 ]\n[ 0.0002288926879375552? 0.02114409313997740? ]\nrho 0.0290306767053235?\nautonorm [0.0078248720579728?, 0.992175127942028?]\nC_s 127.797616701107? C_w 1.007886583565346?\nl 51 size of tail 0.00241776775667970?\nNorm of the  18  iterate  6.35202208669e-05\napproximation coefficients  0.000114460882593563? ||f||_s 0.02112293750243127? ||f||_w\n[ 3.8146972656250000?e-6 1.00000000000000 ]\n[ 0.000114460882593563? 0.02118645772329814? ]\nrho 0.02564960312165553?\nautonorm [0.0044433142408499?, 0.995556685759151?]\nC_s 225.057231110606? C_w 1.004463145398358?\nl 54 size of tail 0.00175356279008919?\nNorm of the  19  iterate  0.000190560662487\napproximation coefficients  0.0000572740571710950? ||f||_s 0.02112299081072210? ||f||_w\n[ 1.9073486328125000?e-6 1.00000000000000 ]\n[ 0.0000572740571710950? 0.02131355147320866? ]\nrho 0.02372756494906876?\nautonorm [0.00240820004849048?, 0.997591799951510?]\nC_s 415.24789463684? C_w 1.002414013475860?\nl 57 size of tail 0.00145940974764983?\n"}︡{"done":true}︡
︠54517590-487a-4c66-b316-70a3fac68e37︠









