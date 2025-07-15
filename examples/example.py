import numpy as np
from numpy.linalg import norm
from quantum_mereology import *


D = 2**10
DA = 2 # dimention of subsystem A
H = GOE(D) # intial Hamiltonian
DB = D//DA

# we partition H into diagonal sectors A and B
# U is the unitary matrix that partitions H
res = partition(H, DA, verbose=True)
U = res["U"]
EA = res["EA"]
EB = res["EB"]
Eint = res["Eint"]

# we build the diagonal matrices correspondig to sectors A and B :
# HA=EA⊗I, HB=I⊗EB
HA = np.kron(np.diag(EA), np.identity(DB))
HB = np.kron(np.identity(DA), np.diag(EB))

# interaction Hamiltonian
Hint = np.diag(Eint)

# Check that U partitions H into HA, HB and Hint
Hp = U@H@U.T.conj()
print(norm(Hp - HA-HB-Hint))

# Check how good the partition is
print(norm(Hint)/norm(H))
