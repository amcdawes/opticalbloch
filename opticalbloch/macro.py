# -*- coding: utf-8 -*-

""" Functions for calculating macroscopic properties of atomic ensembles from
solutions of the optical Bloch equations.

Thomas Ogden <t@ogden.eu>
"""

def calc_susceptibility(tdme, E, N, coh):
    """
    Args:
        tdme: Transition dipole matrix element [e a_0]
        E: Electric field amplitude [V/m]
        N: Number density [/m3]
        coh: Coherence between levels []

    Returns:
        Susceptibility []

    """

    a_0 = si.physical_constants['Bohr radius'][0]

    tdme_Cm = tdme*si.e*a_0 # [C m]

    return (2.*N*tdme_Cm/si.epsilon_0/E)*coh # []