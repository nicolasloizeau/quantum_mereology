import numpy as np
from numpy.linalg import norm, eigh, eigvalsh
from scipy.optimize import minimize
from time import time
from .permutations import *

def buildHA(A, DB):
    """eqivalent to np.diag(np.kron(np.diag(A), np.identity(len(B))))"""
    return np.repeat(A, DB)

def buildHB(B, DA):
    return np.tile(B, DA)


class Partitioner:
    def __init__(self, spectrum, DA, disp=False, x0=None, maxiter=1000):
        self.p0 = np.argsort(spectrum)
        self.spectrum = spectrum[self.p0]
        self.D = len(spectrum)
        self.history = []
        self.disp = disp
        self.t0 = time()
        self.times = []
        self.iter = 0
        self.maxiter = maxiter
        self.DA = DA
        self.DB = self.D//DA
        self.x0 = x0
        if self.x0 is None:
            self.x0 = np.random.random(self.DA+self.DB)-0.5


    def callback(self, x):
        c = self.cost(x)
        self.history.append(c)
        self.iter += 1
        self.times.append(time()-self.t0)
        if self.disp:
            print(self.iter, c)

    def get_AB(self,x):
        return x[:self.DA],x[self.DA:]

    def cost(self, x):
        A,B = self.get_AB(x)
        AB = np.sort(buildHA(A, self.DB)+buildHB(B, self.DA))
        c = norm(AB-self.spectrum)**2
        self.history.append(c)
        return c

    def gradient(self, x):
        A,B = self.get_AB(x)
        DB = np.full(self.DA,self.DB)
        DA = np.full(self.DB,self.DA)
        D = np.concatenate((DB,DA))
        G1 = D*2*x
        G2 = 2*np.concatenate((np.full(self.DA, np.sum(B)), np.full(self.DB, np.sum(A))))
        AB = buildHA(A,self.DB)+buildHB(B, self.DA)
        p = np.argsort(AB)
        p = inversep(p)
        spectrum = np.reshape(self.spectrum[p], (self.DA, self.DB))
        G3 = -2*np.concatenate([spectrum.sum(axis=1), spectrum.sum(axis=0)])
        return G1+G2+G3

    def gradient_fd(self, h, d=1e-6):
        grad = np.zeros(len(h))
        for i in range(len(h)):
            dh = np.zeros(len(h))
            dh[i] = d
            grad[i] = (self.cost(h+dh)-self.cost(h-dh))/(2*d)
        return grad

    def partition(self):
        res = minimize(self.cost, self.x0, tol=1e-10, callback=self.callback, jac=self.gradient, method='BFGS',
        options={'gtol':1e-16, 'maxiter':self.maxiter, 'disp':self.disp})
        A,B = self.get_AB(res.x)
        AB = buildHA(A, self.DB)+buildHB(B, self.DA)
        p = inversep(np.argsort(AB))
        return A,B, composep(p, self.p0)


def partition(H, DA, verbose=False, maxiter=1000, x0=None):
    """
    Partitions the Hamiltonian H into two diagonal sectors A and B
    H : a numpy array
    DA : dimension of subsystem A
    """
    assert len(H)%DA == 0, "H must be a square matrix with dimension divisible by DA"
    E, V = eigh(H)
    part = Partitioner(E, DA, disp=verbose, maxiter=maxiter, x0=x0)
    EA,EB, p = part.partition()
    U = permutation_matrix(p).T.conj()@V.T.conj()
    H2 = U@H@U.T.conj()
    HA = buildHA(EA, len(EB))
    HB = buildHB(EB, len(EA))
    Eint = np.diag(H2)-HA-HB
    result = {}
    result["U"] = U
    result["EA"] = EA
    result["EB"] = EB
    result["Eint"] = Eint
    return result
