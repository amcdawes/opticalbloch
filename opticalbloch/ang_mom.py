# -*- coding: utf-8 -*-

import sys

import numpy as np

from wigner import Wigner3j, Wigner6j

def calc_clebsch_hf(J_a, I_a, F_a, mF_a, J_b, I_b, F_b, mF_b, q):
    """ Clebsch-Gordan coefficient for the hyperfine transition dipole matrix
        element """

    coeff_F = ((-1)**(F_b+J_a+1+I_a)*
              np.sqrt((2*F_b+1)*(2*J_a+1))*
              Wigner6j(J_a,J_b,1,F_b,F_a,I_a)) # Steck Rb87, eqn 36

    coeff_hf = ((-1)**(F_b-1+mF_a)*
               np.sqrt(2*F_a+1)*
               Wigner3j(F_b,1,F_a,mF_b,q,-mF_a)) # Steck Rb87, eqn 35

    return coeff_hf*coeff_F

def main():
    pass

if __name__ == '__main__':
    status = main()
    sys.exit(status)