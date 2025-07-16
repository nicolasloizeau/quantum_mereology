.. Quantum Mereology documentation master file

Quantum Mereology
=================

This package focuses on the following problem :
Given a many body spin 1/2 Hamiltonian :math:`H`, find a unitary :math:`U` such that
:math:`U H U^\dagger` has a particular tensor product structure.
In particular, minimize the cost function
:math:`C(h) = \sum_{n} |E_n-\varepsilon_n|^2`
where :math:`E_n` are the eigenvalues of :math:`H` and
:math:`\varepsilon_n` are the eigenvalues of
:math:`H'=\sum_i h_i \tau_i` and :math:`\tau_i` are a set of Pauli strings
that represent a particular tensor product structure.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   example_partitioning
   example_localization
   reference
