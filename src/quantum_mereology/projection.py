import os
import numpy as np

def projection0(H, ps):
    """
    basic unoptimized projection
    """
    Hp = 0
    for p in ps:
        p = p.toarray()
        Hp += np.trace(H@p)*p
    return Hp/len(H)

def projection(H, ps, dtype='float64'):
    """
    projection optimized for sparse matrices
    """
    d = len(H)
    Hp = np.zeros((d,d), dtype=dtype)
    for p in ps:
        c = inner_sparse(H, p, dtype=dtype)
        Hp[p.row, p.col] += c*p.data
    return Hp/d


def projection_diag(H, ps):
    d = len(H)
    Hp = np.zeros(d, dtype=float)
    for p in ps:
        c = np.sum(H*p)
        Hp += c*p
    return Hp/d


def projection_coefs(H, ps, dtype='float64'):
    """
    projection optimized for sparse matrices
    """
    d = len(H)
    hs = []
    for p in ps:
        hs.append(inner_sparse(H, p, dtype=dtype).real)
    return np.array(hs)



def inner_sparse(H, P, dtype='float64'):
    p = P.data.conj() if dtype=='complex128' else P.data
    # p = P.data
    h = np.array(H[P.row, P.col])
    return np.sum(h*p)



if __name__ == '__main__':
    from scipy.linalg import norm
    from projectors import*
    from randq import*
    N = 8
    for dtype in ('float64', 'complex128'):
        ps = get_1_2_body(N, dtype=dtype)
        H = randH_real(N)
        print(norm(projection(H, ps, dtype=dtype)-projection0(H, ps)))

    N = 2
    p = get_1_2_body(N, dtype='complex128')[1]
    print(norm(projection0(p.toarray(), [p])-p))
    print(norm(projection(p.toarray(), [p], dtype='complex128')-p))
