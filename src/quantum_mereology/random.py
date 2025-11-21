import numpy as np


def GOE(D):
    """
    Generates a random matrix from the Gaussian Orthogonal Ensemble (GOE).

    Parameters:
    -----------
    D : int
        Dimension of the square matrix.

    Returns:
    --------
    H : ndarray
    """
    H = np.random.normal(size=(D, D), scale=1, loc=0)
    return H + H.T
    return H


def GUE(D):
    """
    Generates a random matrix from the Gaussian Unitary Ensemble (GUE).

    Parameters:
    -----------
    D : int
        Dimension of the square matrix.

    Returns:
    --------
    H : ndarray
    """
    H = np.random.normal(size=(D, D), scale=1, loc=0)
    H = H + np.random.normal(size=(D, D), scale=1, loc=0) * 1j
    return H + H.conj().T
