==================
Quantum Mereology
==================



.. image:: https://github.com/nicolasloizeau/quantum_mereology/actions/workflows/test.yml/badge.svg
   :target: https://github.com/nicolasloizeau/quantum_mereology/actions/workflows/test.yml
   :alt: Tests Status
.. image:: https://img.shields.io/badge/docs-blue.svg
   :target: https://nicolasloizeau.github.io/quantum_mereology/
   :alt: Documentation

Numerical methods for quantum mereology

This is a Python package that focuses on the following problem:
Given a many body spin 1/2 Hamiltonian *H*, find a unitary *U* such that *U H U*\ :sup:`†` has a particular tensor product structure.
In particular, minimize the cost function
*C(h) = Σ_n |E_n-ε_n|*\ :sup:`2`
where *E_n* are the eigenvalues of *H* and *ε_n* are the eigenvalues of *H'=Σ_i h_i τ_i* and *τ_i* are a set of Pauli strings that represent a particular tensor product structure.

Installation
------------

From github:

.. code-block:: bash

    pip install git+https://github.com/nicolasloizeau/quantum_mereology.git

Example
-------

.. code-block:: python

    # Example code will go here

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
