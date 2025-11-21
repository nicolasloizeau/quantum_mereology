#!/usr/bin/env python
import pytest

"""Tests for `quantum_mereology` package."""

from quantum_mereology.partitioning import *
from quantum_mereology.random import GOE


def test_gradient():
    """test that the gradient aggrees with the finite difference"""
    N = 8
    E = eigvalsh(GOE(N))
    part = Partitioner(E, 2)
    x0 = np.random.random(part.DA + part.DB)
    g1 = part.gradient(x0)
    g2 = part.gradient_fd(x0, d=1e-6)
    print(norm(g1 - g2) / norm(g1))
    assert norm(g1 - g2) / norm(g1) < 1e-8, "gradient and finite difference do not agree"


def test_gradientB():
    """test that the gradient aggrees with the finite difference"""
    N = 8
    E = eigvalsh(GOE(N))
    part = PartitionerB(E, 2)
    x0 = np.random.random(part.DB)
    g1 = part.gradient(x0)
    g2 = part.gradient_fd(x0, d=1e-6)
    print(norm(g1 - g2) / norm(g1))
    assert norm(g1 - g2) / norm(g1) < 1e-8, "gradient and finite difference do not agree"


def run_partitioner(optimizeA=True):
    D = 2**8
    DA = 2  # dimention of subsystem A
    H = GOE(D)  # intial Hamiltonian
    DB = D // DA

    # we partition H into diagonal sectors A and B
    # U is the unitary matrix that partitions H
    res = partition(H, DA, verbose=True, optimizeA=optimizeA)
    U = res.U
    EA = res.EA
    EB = res.EB
    Eint = res.Eint

    # we build the diagonal matrices correspondig to sectors A and B :
    # HA=EA⊗I, HB=I⊗EB
    HA = np.kron(np.diag(EA), np.identity(DB))
    HB = np.kron(np.identity(DA), np.diag(EB))

    # interaction Hamiltonian
    Hint = np.diag(Eint)

    # Check that U partitions H into HA, HB and Hint
    Hp = U @ H @ U.T.conj()
    assert norm(Hp - HA - HB - Hint) < 1e-8

    # Check how good the partition is
    assert (norm(Hint) / norm(H)) < 5e-2


def test_partitioner():
    run_partitioner(optimizeA=True)


def test_partitionerB():
    run_partitioner(optimizeA=False)
