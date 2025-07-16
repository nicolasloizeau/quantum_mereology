import numpy as np
from numpy.linalg import norm
from quantum_mereology import *


N = 8
D = 2**N
H = GOE(D)
taus = local2(N)
result = localize(H, taus, maxiter=50, verbose=True)
U = result.U
H2 = U@H@U.T.conj() # localized Hamiltonian
Hloc = result.Hloc # local Hamiltonian

print(norm(H2 - Hloc)/norm(H2))
