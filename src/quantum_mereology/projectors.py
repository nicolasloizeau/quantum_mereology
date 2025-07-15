import os
from scipy.sparse import coo_matrix, kron
import numpy as np
from itertools import product, combinations
import random
from scipy import sparse

paulis = [np.eye(2), np.array([[0,1],[1,0]]), 1j*np.array([[0,-1],[1,0]]), np.array([[1,0],[0,-1]])]
paulis_sparse = [coo_matrix(p, dtype='complex128') for p in paulis]

def proj_from_indexes(indexes, dtype='float64'):
    proj = paulis_sparse[indexes[0]]
    for i in indexes[1:]:
        proj = sparse.kron(proj, paulis_sparse[i], format='coo')
    if dtype=='float64':
        proj = proj.real
    return coo_matrix(proj, dtype=dtype)


def majorana(N, k):
    pauli_indexes = np.zeros(N, dtype=int)
    pauli_indexes[k//2] = 1 if k%2==0 else 2
    pauli_indexes[k//2+1:]+=3
    return proj_from_indexes(pauli_indexes, dtype='complex128')

def majoranas(N):
    return [majorana(N, k) for k in range(2*N)]


def majoranas4(N):
    psis = majoranas(N)
    for i in range(2*N):
        for j in range(i):
            for k in range(j):
                for l in range(k):
                    yield coo_matrix(psis[l]@psis[k]@psis[j]@psis[i])




def local2(N, dtype='float64', sites=[]):
    sites = range(N) if len(sites) == 0 else sites
    projectors = []
    # 1-local
    for i in sites:
        pauli_indexes = np.zeros(N, dtype=int)
        pauli_indexes[i] = 3
        projectors.append(proj_from_indexes(pauli_indexes, dtype=dtype))
    # 2-local
    for i, j in list(combinations(sites, 2)):
        for k, l in product(range(1,4), repeat=2):
            if dtype=='complex128' or not((k==2 and l!=2) or (l==2 and k!=2)):
                pauli_indexes = np.zeros(N, dtype=int)
                pauli_indexes[i] = k
                pauli_indexes[j] = l
                projectors.append(proj_from_indexes(pauli_indexes, dtype=dtype))
    return projectors



def local2only(N, dtype='float64', sites=[]):
    sites = range(N) if len(sites) == 0 else sites
    projectors = []
    # 2-local
    for i, j in list(combinations(sites, 2)):
        for k, l in product(range(1,4), repeat=2):
            if dtype=='complex128' or not((k==2 and l!=2) or (l==2 and k!=2)):
                pauli_indexes = np.zeros(N, dtype=int)
                pauli_indexes[i] = k
                pauli_indexes[j] = l
                projectors.append(proj_from_indexes(pauli_indexes, dtype=dtype))
    return projectors



def localk(N, k, dtype='float64', progress=list):
    projectors = []
    for sites in progress(list(combinations(range(N), k))):
        for ptypes in product(range(1,4), repeat=k):
            if dtype=='complex128' or ptypes.count(2)%2 == 0:
                pauli_indexes = np.zeros(N, dtype=int)
                pauli_indexes[list(sites)] = list(ptypes)
                projectors.append(proj_from_indexes(pauli_indexes, dtype=dtype))
    return projectors



def chain1d(N, dtype='float64'):
    projectors = []
    # 1-local
    for i in range(N):
        pauli_indexes = np.zeros(N, dtype=int)
        pauli_indexes[i] = 3
        projectors.append(proj_from_indexes(pauli_indexes, dtype=dtype))
    # 2-local
    for i in range(N-1):
        for k, l in product(range(1,4), repeat=2):
            if not('float64') or not((k==2 and l!=2) or (l==2 and k!=2)):
                pauli_indexes = np.zeros(N, dtype=int)
                pauli_indexes[i] = k
                pauli_indexes[i+1] = l
                projectors.append(proj_from_indexes(pauli_indexes, dtype=dtype))
    return projectors



def local1(N, dtype='float64'):
    projectors = []
    # 1-local
    for i in range(N):
        pauli_indexes = np.zeros(N, dtype=int)
        pauli_indexes[i] = 3
        projectors.append(proj_from_indexes(pauli_indexes, dtype=dtype))
    return projectors


def allstrings(N, dtype='float64'):
    strings = product([0,1,2,3], repeat=N)
    strings2 = []
    for string in strings:
        if dtype=='complex128' or string.count(2)%2 == 0:
            strings2.append(string)
    return np.array(strings2)



def partition(N, dtype='float64'):
    strings = allstrings(N//2, dtype=dtype)[1:]
    projectors = []
    # A side
    for string in strings:
        pauli_indexes = np.concatenate((np.zeros(N//2, dtype=int), string))
        projectors.append(proj_from_indexes(pauli_indexes, dtype=dtype))
    # B side
    for string in strings:
        pauli_indexes = np.concatenate((string, np.zeros(N//2, dtype=int)))
        projectors.append(proj_from_indexes(pauli_indexes, dtype=dtype))
    one = np.zeros(N, dtype=int)
    projectors.append(proj_from_indexes(one, dtype=dtype))
    return projectors


def checktrace(ps):
    for i in range(len(ps)):
        for j in range(len(ps)):
            if i != j:
                p1 = ps[i].todense()
                p2 = ps[j].todense()
                if np.trace(p1@p2) > 1e-5:
                    return False
    return True

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    # taus = regular5(10)
    # G = nx.random_regular_graph(2, 10, seed=1)
    # print(G.edges)
    # nx.draw(G)
    # plt.show()
    strings = allstrings(8)
    print(len(strings))
