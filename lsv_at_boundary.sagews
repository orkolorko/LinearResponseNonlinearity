︠33110740-d649-4c37-9b4b-c9feb51bcee6s︠
size=32768
eta=1.0/size

from partition import *
from ulam import *

b_Ulam=UlamL1(equispaced(size))
︡2bf100ef-8daa-4c5b-8773-412fc62bb55b︡{"done":true}︡
︠f4bc0215-576b-49d6-8648-df4a3ee4e5ffs︠
from dynamic import *
from generic_assembler import *

prec=53

D=QuotientedLinearDynamic(k=2,prec=prec)
P_Ulam=assemble(D,b_Ulam,2**(-prec+12),output_rate=65536)
︡49ae0ae9-d6f3-4653-9992-01ca6a32bd22︡{"stdout":"Epsilon  1/2199023255552\ninterval matrix\nstart\n0\n0\n"}︡{"done":true}︡
︠435bf52a-a84e-44a7-8c30-6970e5ecd6d8s︠

from linear_lsv_at_the_boundary import *
︡23cecab9-0592-4901-84e1-0d4938901c6a︡{"done":true}︡
︠2592e895-a315-4480-aedb-f9b7fdf1ebe2s︠

f_eta=compute_f_eta_lsv_at_boundary(b_Ulam)
︡8b535b2b-06a0-452a-9791-f8633f460b02︡{"done":true}︡
︠493eb1ac-c33a-49aa-877b-4203f11d4aaas︠
pl=plot(lambda x: b_Ulam.evaluate_linear_combination(f_eta,RIF(x)).center(),0,1)
pl.show()
︡9b54d311-a04f-4297-ad2d-9e8e410c70d8︡{"file":{"filename":"/projects/f414be55-bd9f-499e-8370-d3ca3e1c83cd/.sage/temp/9aca8c334dfe/2600/tmp_4R0ydl.svg","show":true,"text":null,"uuid":"a19d2c77-6d49-42ff-8e2e-bfa7eb9186dd"},"once":false}︡{"done":true}︡
︠5ca5c47e-74c6-4311-8b0c-77b108779d12s︠
#The peak at the end is a problem with the plot function, as it is possible to see with the coefficients of f_eta
print f_eta
︡a71405e8-320b-4401-8333-d26632b9e054︡{"stdout":"[ 2.59930193  2.25272834  2.1219163  ..., -0.24998093 -0.24998856\n -0.24999619]\n"}︡{"done":true}︡
︠c42d661f-23aa-4277-bf06-9665876e3c43s︠
print compute_variation_f_eta(b_Ulam)
︡4f5d66b8-f9b1-471d-8c33-933807d11d5a︡{"stdout":"2.84929811236\n"}︡{"done":true}︡
︠3714c400-b7cb-4f0a-b7b0-cced180bf951s︠
l=20
err_tail = bound_size_tail_lsv_at_boundary(l)
print err_tail
︡7e810047-d2a2-4f1d-89bb-ca7ac0c023f2︡{"stdout":"4.38546370604801e-6\n"}︡{"done":true}︡
︠045a88f7-e734-404a-a1ca-6927ed171e8ds︠
err_sum_L_i = compute_error_sum_L_i_f_eta_minus_hat_L_h(b_Ulam,l)
print err_sum_L_i
︡b8cd0aed-b891-4f92-a9d1-9059fb7f5dee︡{"stdout":"0.00169875588848\n"}︡{"done":true}︡
︠b9b97dc7-f578-4be6-88c9-58909c1961fes︠
sum_iter_error = compute_sum_L_i_f_eta_minus_L_eta_i_f_eta(b_Ulam,l,f_eta)
print sum_iter_error
︡a8ad6da5-011c-429e-8b85-0beae89b18f1︡{"stdout":"0.00495635988027\n"}︡{"done":true}︡
︠3bb738e9-1f6c-4ef8-9961-015e5e83e6e1s︠
tot_error = err_tail+err_sum_L_i+sum_iter_error
print tot_error
︡35d1c81a-00a0-430a-ae96-c8659acf60ca︡{"stdout":"0.00665950123245452\n"}︡{"done":true}︡
︠57c30fa9-e656-4c10-ab72-ff5eafb98489s︠
z = compute_linear_response(b_Ulam,P_Ulam,f_eta,l, interval_matrix=True)
︡caadd8bf-9078-45be-ba59-46c6653a1010︡{"done":true}︡
︠27d43a2e-8b96-42ab-b6a0-3fd081666f34s︠
pl=plot(lambda x: b_Ulam.evaluate_linear_combination(z,RIF(x)).center(),0,1)
pl.show()
︡42fe8177-4c69-4165-bc7c-4a6a0500a97c︡{"file":{"filename":"/projects/f414be55-bd9f-499e-8370-d3ca3e1c83cd/.sage/temp/9aca8c334dfe/2600/tmp_Dm1wKQ.svg","show":true,"text":null,"uuid":"af5dcd8b-9e6e-4129-874d-d1674a028636"},"once":false}︡{"done":true}︡
︠c22544d4-01db-4a93-bfc3-63798c496d90s︠
#As before, the step at the end is an artifact
print z
︡cfde90b6-ebfa-4a1d-acec-743b991803c5︡{"stdout":"[ 5.03830173  4.69172814  4.38762168 ..., -0.65948043 -0.6594995\n -0.65950713]\n"}︡{"done":true}︡









