# -*- coding: utf-8 -*-

import os
import sys

import numpy as np
import qutip as qu

import matplotlib.pyplot as plt

import ob_two

pi = np.pi

base = os.path.basename(__file__)
savename = os.path.splitext(base)[0]

def solve():
    """ Light interaction with a Two-level Atom. 
        Rabi Flopping, No Decay. 
    """

### TwoOB

    two_obj = ob_two.OBTwo(gammas=[2*pi*1.])

    two_obj.set_H_Delta([2*pi*0.])
    two_obj.set_H_Omega([2*pi*10.])

### Solve with QuTiP

    tlist = np.linspace(0, 1, 201)

    two_obj.mesolve(tlist)

    return two_obj

def plot(two_obj):

    t = two_obj.result.times
    states_t = two_obj.states_t()

    pop_1 = states_t[:,1,1].real
    coh_im = states_t[:,1,0].imag

### Plot

    fig = plt.figure()
    ax = fig.add_subplot(2,1,1)
    ax.plot(t, pop_1)
    ax.set_ylabel(r'Excited Population')
    ax.set_ylim([0.,1.])
    ax = fig.add_subplot(2,1,2)
    ax.plot(t, coh_im, c='r')
    ax.set_xlabel(r'Time ($\tau$)')
    ax.set_ylabel(r'Imag Coherence')
    ax.set_ylim([-.5,.5])
    plt.draw()

def main():
    two_obj = solve()
    plot(two_obj)

### Save and show figures

    for n in plt.get_fignums():
        plt.figure(n).savefig('png/'+savename+'_fig'+str(n)+'.png')
        plt.figure(n).savefig('pdf/'+savename+'_fig'+str(n)+'.pdf')    

    plt.show()

if __name__ == "__main__":
    status = main()
    sys.exit(status)