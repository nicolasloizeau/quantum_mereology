Localization example
====================

Given a random matrix H, we want to find a Unitary U that transform it to a 2-local form

Generate and GOE matrix of 10 spins :

.. code-block:: python

    import numpy as np
    from numpy.linalg import norm
    from quantum_mereology import *
    N = 8
    D = 2**N
    H = GOE(D)

Construct a list of 2-local Pauli strings

.. code-block:: python

    taus = local2(N)

Solve the localization problem

.. code-block:: python

    result = localize(H, taus, maxiter=50, verbose=True)

Check the result

.. code-block:: python

    U = result.U
    H2 = U@H@U.T.conj() # localized Hamiltonian
    Hloc = result.Hloc # local Hamiltonian
    print(norm(H2 - Hloc)/norm(H2))
