
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.optimize import minimize
from numpy.linalg import norm
from scipy.linalg import eigh, eigvalsh
from .pauli_strings import *

class Localizer:

    def __init__(self, E0, taus, maxiter=1000, gtol=1e-6, x0=[], verbose=True):
        self.N = int(np.log2(len(E0)))
        self.taus = taus
        self.iter = 0
        self.history = []
        self.maxiter = maxiter
        self.gtol = gtol
        self.E0 = E0
        self.h0 = x0
        if len(self.h0) == 0:
            self.h0 = (np.random.random(len(self.taus))-0.5)/len(self.taus)
        self.disp = verbose
        self.current_cost = 1


    def cost(self, h):
        E = eigvalsh(buildH(self.taus, h))
        c = norm(E-self.E0)**2
        self.current_cost = c
        return c


    def callback(self, h):
        c = self.current_cost
        self.history.append(c)
        self.iter += 1
        if self.disp:
            print(self.iter, c)

    def localize(self):
        res = minimize(self.cost, self.h0, method='BFGS',
        callback=self.callback, jac=self.gradient,
        options={'gtol':1e-16, 'maxiter':self.maxiter, 'disp':True}, tol=1e-16)
        return res.x


    def gradient(self, h):
        e,v = eigh(buildH(self.taus, h))
        H2 = v@np.diag(self.E0)@v.T
        grad = np.zeros(len(h))
        for i in range(len(h)):
            grad[i] = 2*(h[i]*2**self.N-np.sum(H2[self.taus[i].row, self.taus[i].col]*self.taus[i].data))
        return grad




def localize(H, taus, maxiter=200, verbose=False):
    E, V = eigh(H)
    loc = Localizer(E, taus, maxiter=maxiter, verbose=verbose)
    h = loc.localize()
    res = {}
    res['h'] = h
    res['H'] = buildH(taus, h)
    return res
