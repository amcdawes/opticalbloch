# -*- coding: utf-8 -*-

import os
import sys

import numpy as np
import qutip as qu

import matplotlib.pyplot as plt

import ob_two
import opticalbloch.t_funcs as tf

pi = np.pi

base = os.path.basename(__file__)
savename = os.path.splitext(base)[0]

def square_1(t, args):

    on_1 = args['on_1']
    off_1 = args['off_1']
    ampl_1 = args['ampl_1']

    return ampl_1*(t >= on_1)*(t <= off_1)

def solve():

### TwoOB

    two_obj = ob_two.OBTwo(gammas=[2*pi*0.])

    two_obj.set_H_Delta([2*pi*0.])
    two_obj.set_H_Omega([2*pi*10.], [tf.square_1])

    args = {'on_1':0.,
            'off_1':.5, 'ampl_1':1.}

### Solve with QuTiP

    tlist = np.linspace(0., 1., 201)

    two_obj.mesolve(tlist, td=True, args=args)

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

    # for n in plt.get_fignums():
    #     plt.figure(n).savefig('png/'+savename+'_fig'+str(n)+'.png')
    #     plt.figure(n).savefig('pdf/'+savename+'_fig'+str(n)+'.pdf')    

    plt.show()

if __name__ == "__main__":
    status = main()
    sys.exit(status)