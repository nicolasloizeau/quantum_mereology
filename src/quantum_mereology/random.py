
import numpy as np



def GOE(D):
    H = np.random.normal(size=(D,D), scale=1, loc=0)
    return H+H.T
    return H


def GUE(D):
    H = np.random.normal(size=(D,D), scale=1, loc=0)
    H = H+ np.random.normal(size=(D,D), scale=1, loc=0)*1j
    return H+H.conj().T
