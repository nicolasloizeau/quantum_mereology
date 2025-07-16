Partitioning example
====================

Lets partition a random matrix into two sectors

First we construct a random GOE matrix H of dimentions D and specify DA the dimension of subsystem A.

.. code-block:: python

    import numpy as np
    from numpy.linalg import norm
    from quantum_mereology import *
    D = 2**10
    DA = 2 # dimention of subsystem A
    H = GOE(D) # intial Hamiltonian
    DB = D//DA

Partition the Hamiltonian :

.. code-block:: python

    res = partition(H, DA, verbose=True)
    U = res.U
    EA = res.EA
    EB = res.EB
    Eint = res.Eint

Now we can construct HA and HB and check how good the partitioning is

.. code-block:: python

    HA = np.kron(np.diag(EA), np.identity(DB))
    HB = np.kron(np.identity(DA), np.diag(EB))

    # interaction Hamiltonian
    Hint = np.diag(Eint)

    # Check that U partitions H into HA, HB and Hint
    Hp = U@H@U.T.conj()
    print(norm(Hp - HA-HB-Hint))

    # Check how good the partition is
    print(norm(Hint)/norm(H))
