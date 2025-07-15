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
    x0 = np.random.random(part.DA+part.DB)
    g1 = part.gradient(x0)
    g2 = part.gradient_fd(x0, d=1e-6)
    print(norm(g1-g2)/norm(g1))
    assert norm(g1-g2)/norm(g1) < 1e-8, "gradient and finite difference do not agree"
