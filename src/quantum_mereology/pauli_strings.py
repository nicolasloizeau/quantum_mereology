import os
from scipy.sparse import coo_matrix, kron
import numpy as np
from itertools import product, combinations



paulis = [np.eye(2), np.array([[0,1],[1,0]]), 1j*np.array([[0,-1],[1,0]]), np.array([[1,0],[0,-1]])]
paulis_sparse = [coo_matrix(p, dtype='complex128') for p in paulis]

def operator_from_indexes(indexes, dtype='float64'):
    """
    indexes : list of pauli string indexes (eg [0,1,2,0,3])
    return : coo_matrix representing a pauli string (eg 1XY1Z)
    """
    op = paulis_sparse[indexes[0]]
    for i in indexes[1:]:
        op = kron(op, paulis_sparse[i], format='coo')
    if dtype=='float64':
        op = op.real
    return coo_matrix(op, dtype=dtype)



def local2(N, dtype='float64'):
    """
    Generate a list of 2-local Pauli strings of N qubits.

    Parameters
    ----------
    N : int
        Number of qubits.
    dtype : str, optional
        Data type for the operators. Use 'complex128' for complex operators,
        'float64' for real operators. Default is 'float64'.

    Returns
    -------
    list
        List of 2-local Pauli string operators as sparse matrices.

    Notes
    -----
    If dtype='float64', complex strings (those with a single Y operator) are excluded.
    """
    taus = []
    # 1-local
    for i in range(N):
        for k in range(1,4):
            pauli_indexes = np.zeros(N, dtype=int)
            pauli_indexes[i] = k
            taus.append(operator_from_indexes(pauli_indexes, dtype=dtype))
    # 2-local
    for i, j in list(combinations(range(N), 2)):
        for k, l in product(range(1,4), repeat=2):
            # exclude complex strings if dtype=float64
            if dtype=='complex128' or not((k==2 and l!=2) or (l==2 and k!=2)):
                pauli_indexes = np.zeros(N, dtype=int)
                pauli_indexes[i] = k
                pauli_indexes[j] = l
                taus.append(operator_from_indexes(pauli_indexes, dtype=dtype))
    return taus



def local1Z(N):
    """
    Generate a list of Z operators N qubits.

    Parameters
    ----------
    N : int
        Number of qubits.

    Returns
    -------
    list
        List of Zi Pauli string operators as sparse matrices.

    """
    taus = []
    for i in range(N):
        pauli_indexes = np.zeros(N, dtype=int)
        pauli_indexes[i] = 3
        taus.append(operator_from_indexes(pauli_indexes, dtype='float64'))
    return taus



def local1(N, dtype='float64'):
    """
    Generate a list of 1-local Pauli strings of N qubits.

    Parameters
    ----------
    N : int
        Number of qubits.
    dtype : str, optional
        Data type for the operators. Use 'complex128' for complex operators,
        'float64' for real operators. Default is 'float64'.

    Returns
    -------
    list
        List of 1-local Pauli string operators as sparse matrices.

    Notes
    -----
    If dtype='float64', complex strings (those with a Y operator) are excluded.
    """
    taus = []
    for i in range(N):
        for k in range(1,4):
            pauli_indexes = np.zeros(N, dtype=int)
            pauli_indexes[i] = k
            taus.append(operator_from_indexes(pauli_indexes, dtype=dtype))
    return taus


def buildH(taus, h):
    """
    Build a Hamiltonian from a list of pauli strings and coefficients.
    Return sum_i h_i * tau_i

    Parameters
    ----------
    taus : list
        List of operators (as scipy.sparse.coo_matrix).
    h : array_like
        Coefficients for each operator.

    Returns
    -------
    numpy.ndarray
        The Hamiltonian matrix constructed as the weighted sum of operators.

    Notes
    -----
    This function efficiently computes the Hamiltonian as:
    H = sum_i h_i * tau_i
    where h_i are the coefficients and tau_i are the operators.
    """
    d = taus[0].shape[0]
    H = np.zeros((d,d), dtype=taus[0].dtype)
    for i in range(len(taus)):
        H[taus[i].row, taus[i].col] += h[i]*taus[i].data
    return H
