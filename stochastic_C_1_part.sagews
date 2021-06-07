︠0314b26b-6588-42d8-bf80-f79de32b7745︠
print "In this file we compute the derivative of the invariant density, to compute the linear response in the stocastich perturbation case"
︡4c083725-02f4-420c-a308-f4f691ef8df1︡{"stdout":"In this file we compute the derivative of the invariant density, to compute the linear response in the stocastich perturbation case\n"}︡{"done":true}︡
︠2753d58f-e2a3-471b-9093-24bba98cc33es︠
from spline_basis_C1 import *
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
size = 16384

plot(lambda x: dynamic.f(RIF(x)).center(),0,1)

from generic_assembler import *

b_C1=Spline_Basis_C1(equispaced(size))
(A,B,lam)=b_C1.dfly(dynamic)
print "||L^n f||_{C^2}\leq ",A,"*",lam,"^n+",B,"||f||_{C^1}"
(A_unif,B_unif,lam_unif)=b_C1.dfly(dynamic, discrete=True)
print "||L_{\eta}^n f||_{C^2}\leq ",A_unif,"*",lam_unif,"^n+",B_unif,"||f||_{C^1}"
︡3415c5b5-df94-4466-8c21-e4ae634c97e5︡{"stderr":"Compiling ./binsearch2.spyx...\nCompiling ./cython_sparse_mat_vec.spyx..."}︡{"stderr":"\nCompiling ./polyval_lib.spyx..."}︡{"stderr":"\n"}︡{"file":{"filename":"/projects/f414be55-bd9f-499e-8370-d3ca3e1c83cd/.sage/temp/9aca8c334dfe/2320/tmp_TAdNIY.svg","show":true,"text":null,"uuid":"3cd73960-f240-43f5-8e1f-1031342afae9"},"once":false}︡{"stdout":"||L^n f||_{C^2}\\leq  1.19979547394928 * 0.01602109207880304? ^n+ 3.78861754934280? ||f||_{C^1}\n"}︡{"stdout":"||L_{\\eta}^n f||_{C^2}\\leq  4.00002417970955 * 0.0768886417481990? ^n+ 238.653460146613? ||f||_{C^1}\n"}︡{"done":true}︡
︠497e384e-c1b0-4bfd-a1ab-f1c4225f0f66s︠
if assemble_matrix is True:
    from generic_assembler import *
    P_C1,P_C1_radius=assemble(dynamic,b_C1,2**(-prec+8),output_rate=size*8/256,scipy_matrix=True,scipy_radius=True)
    b_C1.save_csr_matrix_with_radius(P_C1,P_C1_radius,name)
︡30c39bf4-bae1-4d51-89ef-cd9f0bd309c3︡{"stdout":"Epsilon  1/35184372088832\nsparse matrix,double, csr,scipy\nsparse matrix, with radius, double, csr, scipy\nstart"}︡{"stdout":"\n0\n512"}︡{"stdout":"\n1024"}︡{"stdout":"\n1536"}︡{"stdout":"\n2048"}︡{"stdout":"\n2560"}︡{"stdout":"\n3072"}︡{"stdout":"\n3584"}︡{"stdout":"\n4096"}︡{"stdout":"\n4608"}︡{"stdout":"\n5120"}︡{"stdout":"\n5632"}︡{"stdout":"\n6144"}︡{"stdout":"\n6656"}︡{"stdout":"\n7168"}︡{"stdout":"\n7680"}︡{"stdout":"\n8192"}︡{"stdout":"\n8704"}︡{"stdout":"\n9216"}︡{"stdout":"\n9728"}︡{"stdout":"\n10240"}︡{"stdout":"\n10752"}︡{"stdout":"\n11264"}︡{"stdout":"\n11776"}︡{"stdout":"\n12288"}︡{"stdout":"\n12800"}︡{"stdout":"\n13312"}︡{"stdout":"\n13824"}︡{"stdout":"\n14336"}︡{"stdout":"\n14848"}︡{"stdout":"\n15360"}︡{"stdout":"\n15872"}︡{"stdout":"\n16384"}︡{"stdout":"\n16896"}︡{"stdout":"\n17408"}︡{"stdout":"\n17920"}︡{"stdout":"\n18432"}︡{"stdout":"\n18944"}︡{"stdout":"\n19456"}︡{"stdout":"\n19968"}︡{"stdout":"\n20480"}︡{"stdout":"\n20992"}︡{"stdout":"\n21504"}︡{"stdout":"\n22016"}︡{"stdout":"\n22528"}︡{"stdout":"\n23040"}︡{"stdout":"\n23552"}︡{"stdout":"\n24064"}︡{"stdout":"\n24576"}︡{"stdout":"\n25088"}︡{"stdout":"\n25600"}︡{"stdout":"\n26112"}︡{"stdout":"\n26624"}︡{"stdout":"\n27136"}︡{"stdout":"\n27648"}︡{"stdout":"\n28160"}︡{"stdout":"\n28672"}︡{"stdout":"\n29184"}︡{"stdout":"\n29696"}︡{"stdout":"\n30208"}︡{"stdout":"\n30720"}︡{"stdout":"\n31232"}︡{"stdout":"\n31744"}︡{"stdout":"\n32256"}︡{"stdout":"\n32768"}︡{"stdout":"\n"}︡{"done":true}︡
︠f92922c8-2a36-4bd6-8dcf-e6a36bff80f1s︠
P_C1, P_C1_abs_plus_radius, P_C1_radius = b_C1.load_csr_matrix_with_radius_memmapped(name)
︡a80d50f3-aee4-4f1e-b27f-5dc970f4f7c2︡{"done":true}︡
︠c1ebda38-d6fe-4812-8047-81c2c9585e2ds︠
if compute_decay is True:
    from decay_rump import *
    Q_tot=decay_rump_manager(b_C1, P_C1, P_C1_abs_plus_radius, P_C1_radius, N, num_process=num_process, output_rate=1024, start=0, stop=None)
    b_C1.save_Q(Q_tot, name)     
︡78643edf-d917-4120-82db-64e466aaaf75︡{"stdout":"Parallel\n"}︡{"stdout":"start 16380\nstop 18018\n512\n1024\n1536\n"}︡{"stdout":"start 14742\nstop 16380\n512\n1024\n1536\n"}︡{"stdout":"start 11466\nstop 13104\n512\n1024\n1536\n"}︡{"stdout":"start 27846\nstop 29484\n512\n1024\n1536\n"}︡{"stdout":"start 31122\nstop 32770\n512\n1024\n1536\n"}︡{"stdout":"start 4914\nstop 6552\n512\n1024\n1536\n"}︡{"stdout":"start 6552\nstop 8190\n512\n1024\n1536\n"}︡{"stdout":"start 24570\nstop 26208\n512\n1024\n1536\n"}︡{"stdout":"start 8190\nstop 9828\n512\n1024\n1536\n"}︡{"stdout":"start 19656\nstop 21294\n512\n1024\n1536\n"}︡{"stdout":"start 26208\nstop 27846\n512\n1024\n1536\n"}︡{"stdout":"start 3276\nstop 4914\n512\n1024\n1536\n"}︡{"stdout":"start 18018\nstop 19656\n512\n1024\n1536\n"}︡{"stdout":"start 13104\nstop 14742\n512\n1024\n1536\n"}︡{"stdout":"start 1638\nstop 3276\n512\n1024\n1536\n"}︡{"stdout":"start 22932\nstop 24570\n512\n1024\n1536\n"}︡{"stdout":"start 9828\nstop 11466\n512\n1024\n1536\n"}︡{"stdout":"start 21294\nstop 22932\n512\n1024\n1536\n"}︡{"stdout":"start 0\nstop 1638\n512\n1024\n1536\n"}︡{"stdout":"start 29484\nstop 31122\n512\n1024\n1536\n"}︡{"done":true}︡
︠86def653-0a4b-440c-a14a-88269351ba71s︠
N = 2
Q_alt = b_C1.load_Q(name)

Norms = b_C1.norms_from_Q(Q_alt, N)
print Norms
︡a77e3dc1-a973-4a8a-934a-90fdb2e84422︡{"stdout":"[ 1.          0.25703015  0.02827165]\n"}︡{"done":true}︡
︠2d3002d6-4ad4-49d9-8f46-d0153f0ce3das︠
N_approx = 60
C, D = b_C1.approximation_coefficients_from_norms(dynamic, Norms, N_approx)
print C, D
︡629d415a-1459-45a3-a5a6-62b68e0367b1︡{"stdout":"2.37287013062653?e-49 0.00462695717145123?\n"}︡{"done":true}︡
︠4c88c465-4a8d-4fa2-a3a9-4203e58ed3eas︠
v = b_C1.candidate(P_C1,interval_matrix=False)

plot(lambda x: b_C1.evaluate_linear_combination(v,RIF(x)).center(),0,1)
plot(lambda x: b_C1.evaluate_linear_prime(v,RIF(x)).center(),0,1)

err = b_C1.compute_error(dynamic, Norms, N_approx, P_C1, v, n_iter=5)
print err
︡f87403b7-473f-4e9a-87a1-245fb05616ca︡{"file":{"filename":"/projects/f414be55-bd9f-499e-8370-d3ca3e1c83cd/.sage/temp/9aca8c334dfe/23335/tmp_uMVz8f.svg","show":true,"text":null,"uuid":"7a38a188-aaab-4971-9aee-2acf91e9b123"},"once":false}︡{"file":{"filename":"/projects/f414be55-bd9f-499e-8370-d3ca3e1c83cd/.sage/temp/9aca8c334dfe/23335/tmp_G67yGu.svg","show":true,"text":null,"uuid":"12183fea-57f8-494e-9c44-a4ac32645334"},"once":false}︡{"stdout":"0.000217388087787333?\n"}︡{"done":true}︡
︠dbfb670e-fc16-483a-b7d6-dc548e1e47b7s︠
# we take now the derivative of the density and put it in a form readable to the next program
z = [b_C1.evaluate_linear_prime(v,RIF(x)).center() for x in b_C1.partition ]
np.save('derivative_'+str(size)+'.npy',z)
︡23eefc72-53ce-4c39-8a00-decefa47c1e3︡{"done":true}︡









