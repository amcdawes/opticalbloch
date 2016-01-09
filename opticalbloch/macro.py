# -*- coding: utf-8 -*-

""" Functions for calculating macroscopic properties of atomic ensembles from
solutions of the optical Bloch equations.

Thomas Ogden <t@ogden.eu>
"""

from scipy import exp, log10, sqrt, pi
from scipy import constants as si

### Constants

a_0 = si.physical_constants['Bohr radius'][0] # [m] Bohr radius 

def calc_wavenumber(omega_MHz):
    """ Returns the wavenumber for light at a angular freq. 

    Args:
        omega_MHz: angular freq [2π MHz]

    Returns:
        k: wavenumber [/m]
    """

    return omega_MHz*1.e6/si.c # [/m]

def calc_rabi_freq(electric_field_ampl, tdme):
    """ Returns a Rabi frequnecy for a given E field and tdme.

    Args:
        electric_field_ampl [V /m]
        tdme: transition dipole matrix element [e a_0]

    Returns
        rabi_freq [2π MHz]
    """

    a_0 = si.physical_constants['Bohr radius'][0]

    coeff = si.e*a_0/si.hbar

    # tdme_Cm = tdme*si.e*a_0
    # rabi_freq = electric_field_ampl/si.hbar/(2*pi*1e6)*tdme_Cm # [2π MHz]

    rabi_freq = coeff*electric_field_ampl*tdme/1.e6 # [2π MHz]

    return rabi_freq

def calc_electric_field_ampl(rabi_freq, tdme):
    """ TODO doc 

        rabi_freq [2π MHz]

    Note: Make sure the 2π factor IS on the rabi freq.
    """

    a_0 = si.physical_constants['Bohr radius'][0]

    coeff = si.e*a_0/si.hbar

    electric_field_ampl = 1./coeff*rabi_freq*1.e6/tdme # [V /m]

    return electric_field_ampl

def calc_electric_field_ampl_from_intensity(I, refr_index=1.):
    """
    Args:
        I: intensity [W/m2]
        refr_index: refractive index []

    Returns:
        electric field amplitude [V/m]
    """

    return sqrt(I/(refr_index*si.c*si.epsilon_0)) # [V/m]


def calc_susceptibility(tdme, E, N, coh):
    """ Returns the linear susceptibility (χ) for an ensemble of atoms given
    that a field E generates the given coherence (coh).

    Args:
        tdme: Transition dipole matrix element [e a_0]
        E: Electric field amplitude [V/m]
        N: Number density [/m3]
        coh: Coherence between levels []

    Returns:
        susceptibility []
    """

    tdme_Cm = tdme*si.e*a_0 # [C m]

    return (2.*N*tdme_Cm/si.epsilon_0/E)*coh # []

def calc_absorption_coeff(k, chi):
    """ Returns the absorption coefficient as per the Beer-Lambert law, 
    (alpha where I = I_0*exp(-alpha*z)).

    Args:
        k: Wavenumber [/m]
        chi: Susceptibility []

    Returns:
        alpha: Absorption coefficient [/m]
    """

    n = sqrt(1. + chi)
    alpha = 2.*k*n.imag # [/m]

    return alpha

def calc_transmission(alpha, L):
    """ 
    Args:
        alpha: Absorption coefficient [/m]
        L: Length of medium [m]

    Returns:
        Transmission []

    """

    return exp(-alpha*L)

def maxwell_boltzmann(v, width=1.):
    """ Maxwell Boltzmann probability distribution function. """

    return 1./(width*sqrt(pi))*exp(-(v/width)**2)

def calc_N_Rb(T):
    """ Calculates the Rubidium number density [/m3].
        Taken from ElecSus. 

    Args:
        T: Temperature [K]
    """
    
    if T<312.46:
        p=10.0**(4.857-4215./T)
    else:
        p=10.0**(8.316-4275./T-1.3102*log10(T))
    N=101325.0*p/(si.k*T)
    return N