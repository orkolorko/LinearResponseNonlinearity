from dynamic import *
from assembler import *
from plotting import *

D = MannevillePomeauDynamic()
P = assemble(D, 2**9)
p = plot_spectrum(P, 0)
p.show()

