import numpy as np
from numpy.linalg import norm
from quantum_mereology import *



N = 10
D = 2**N
H = GOE(D)
taus = local2(N)
result = localize(H, taus, maxiter=50, verbose=True)
