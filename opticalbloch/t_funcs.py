# -*- coding: utf-8 -*-

""" These are time functions provided for using the time-dependent solver.

Q: Why are there multiple versions of each?
A: The solver will want one list of arguments even if there are multiple 
time-dependent parts to the Hamiltonian. (Say one laser is ramped on then CW
    and another is a Gaussian pulse.) To distinguish the arguments we've got 
    multiple versions. Yes this is wasteful, I would like a better way to do
    it but this works.

Thomas Ogden <t@ogden.eu>
"""

import sys
from numpy import exp, log

def square_1(t, args):

    on_1 = args['on_1']
    off_1 = args['off_1']
    ampl_1 = args['ampl_1']

    return ampl_1*(t >= on_1)*(t <= off_1)

def square_2(t, args):

    on_2 = args['on_2']
    off_2 = args['off_2']
    ampl_2 = args['ampl_2']

    return ampl_2*(t >= on_2)*(t <= off_2)

def gaussian_1(t, args):

    ampl_1 = args['ampl_1']
    width_1 = args['width_1']
    centre_1 = args['centre_1']

    return ampl_1*exp(-2*log(2)*((t - centre_1)/width_1)**2)

def gaussian_2(t, args):

    ampl_2 = args['ampl_2']
    width_2 = args['width_2']
    centre_2 = args['centre_2']

    return ampl_2*exp(-2*log(2)*((t - centre_2)/width_2)**2)

def ramp_on_1(t, args):

    ampl_1 = args['ampl_1']
    width_1 = args['width_1']
    centre_1 = args['centre_1']

    return ampl_1*(exp(-2*log(2)*((t - centre_1)/width_1)**2)*(t <= centre_1) +
                   (t > centre_1))

def ramp_on_2(t, args):

    ampl_2 = args['ampl_2']
    width_2 = args['width_2']
    centre_2 = args['centre_2']

    return ampl_2*(exp(-2*log(2)*((t - centre_2)/width_2)**2)*(t <= centre_2) +
                   (t > centre_2))

def ramp_on_3(t, args):

    ampl_3 = args['ampl_3']
    width_3 = args['width_3']
    centre_3 = args['centre_3']

    return ampl_3*(exp(-2*log(2)*((t - centre_3)/width_3)**2)*(t <= centre_3) +
                   (t > centre_3))

def ramp_off_1(t, args):

    ampl_1 = args['ampl_1']
    width_1 = args['width_1']
    centre_1 = args['centre_1']

    return ampl_1*(exp(-2*log(2)*((t - centre_1)/width_1)**2)*(t >= centre_1) +
                   (t < centre_1))

def ramp_off_2(t, args):

    ampl_2 = args['ampl_2']
    width_2 = args['width_2']
    centre_2 = args['centre_2']

    return ampl_2*(exp(-2*log(2)*((t - centre_2)/width_2)**2)*(t >= centre_2) +
                   (t < centre_2))

def ramp_off_3(t, args):

    ampl_3 = args['ampl_3']
    width_3 = args['width_3']
    centre_3 = args['centre_3']

    return ampl_3*(exp(-2*log(2)*((t - centre_3)/width_3)**2)*(t >= centre_3) +
                   (t < centre_3))

def ramp_onoff_1(t, args):

    ampl_1 = args['ampl_1']
    width_1 = args['width_1']
    on_1 = args['on_1']
    off_1 = args['off_1']

    ramp_on = (exp(-2*log(2)*((t - on_1)/width_1)**2)*
                        (t <= on_1) + (t > on_1))

    ramp_off = (exp(-2*log(2)*((t - off_1)/width_1)**2)*
                        (t >= off_1) + (t < off_1))

    return ampl_1*(ramp_on + ramp_off - 1.)

def ramp_onoff_2(t, args):

    ampl_2 = args['ampl_2']
    width_2 = args['width_2']
    on_2 = args['on_2']
    off_2 = args['off_2']

    ramp_on = (exp(-2*log(2)*((t - on_2)/width_2)**2)*
                        (t <= on_2) + (t > on_2))

    ramp_off = (exp(-2*log(2)*((t - off_2)/width_2)**2)*
                        (t >= off_2) + (t < off_2))

    return ampl_2*(ramp_on + ramp_off - 1.)

def ramp_onoff_3(t, args):

    ampl_3 = args['ampl_3']
    width_3 = args['width_3']
    on_3 = args['on_3']
    off_3 = args['off_3']

    ramp_on = (exp(-2*log(2)*((t - on_3)/width_3)**2)*
                        (t <= on_3) + (t > on_3))

    ramp_off = (exp(-2*log(2)*((t - off_3)/width_3)**2)*
                        (t >= off_3) + (t < off_3))

    return ampl_3*(ramp_on + ramp_off - 1.)