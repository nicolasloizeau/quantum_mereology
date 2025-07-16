Quantum Mereology
================

.. image:: https://github.com/nicolasloizeau/quantum_mereology/actions/workflows/test.yml/badge.svg
   :target: https://github.com/nicolasloizeau/quantum_mereology/actions/workflows/test.yml
   :alt: Tests Status

.. image:: https://img.shields.io/badge/docs-blue.svg
   :target: https://nicolasloizeau.github.io/quantum_mereology/
   :alt: Documentation

Numerical methods for quantum mereology

This package focuses on the following problem :
Given a many body spin 1/2 Hamiltonian :math:`H`, find a unitary :math:`U` such that
:math:`U H U^\dagger` has a particular tensor product structure.
In particular, minimize the cost function
:math:`C(h) = \sum_{n} |E_n-\varepsilon_n|^2`
where :math:`E_n` are the eigenvalues of :math:`H` and
:math:`\varepsilon_n` are the eigenvalues of
:math:`H'=\sum_i h_i \tau_i` and :math:`\tau_i` are a set of Pauli strings
that represent a particular tensor product structure.


Installation
-----------

From github::

    pip install git+https://github.com/nicolasloizeau/quantum_mereology.git

Example
-------
Given H a GOE matrix, find a U such that UHU^\dagger is as 2-local as possible :

.. code-block:: python

    import numpy as np
    from numpy.linalg import norm
    from quantum_mereology import *
    N = 8 #number of qubits
    H = GOE(2**N)
    taus = local2(N) # a list of 2-local pauli strings
    result = localize(H, taus, maxiter=50, verbose=True) # optimization
    H2 = result.U@H@result.U.T.conj() # localized Hamiltonian
    print(norm(H2 - result.Hloc)/norm(H2)) # chech how good the results is


Citation
--------

.. code-block:: bibtex

    @Article{Loizeau2024,
        author={Loizeau, Nicolas
        and Sels, Dries},
        title={Quantum Mereology and Subsystems from the Spectrum},
        journal={Foundations of Physics},
        year={2024},
        month={Dec},
        day={21},
        volume={55},
        number={1},
        pages={3},
        issn={1572-9516},
        doi={10.1007/s10701-024-00813-2},
        url={https://doi.org/10.1007/s10701-024-00813-2}
    }

    @article{Loizeau2023,
        author = {Nicolas Loizeau  and Flaviano Morone  and Dries Sels },
        title = {Unveiling order from chaos by approximate 2-localization of random matrices},
        journal = {Proceedings of the National Academy of Sciences},
        volume = {120},
        number = {39},
        pages = {e2308006120},
        year = {2023},
        doi = {10.1073/pnas.2308006120},
        URL = {https://www.pnas.org/doi/abs/10.1073/pnas.2308006120},
        eprint = {https://www.pnas.org/doi/pdf/10.1073/pnas.2308006120},
    }
