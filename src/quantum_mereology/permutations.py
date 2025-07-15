import numpy as np



def inversep(permutation):
    x = np.empty_like(permutation)
    x[permutation] = np.arange(len(permutation))
    return x

def composep(p1,p2):
    """
    compose by the left side: composep(p1,p2)(x) = p1(p2(x))
    """
    return p2[p1]

def permutation_matrix(p):
    """constructs a permutation matrix from a permutation p"""
    n = len(p)
    mat = np.zeros((n, n))
    for i in range(n):
        mat[p[i], i] = 1
    return mat
