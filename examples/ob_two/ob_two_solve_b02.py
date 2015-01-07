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

def solve():

### TwoOB

    two_obj = ob_two.OBTwo(gammas=[2*pi*0.])

    two_obj.set_H_Delta([2*pi*0.])

    squares = lambda t, args: (tf.square_1(t, args) +
                               tf.square_2(t, args))
    two_obj.set_H_Omega([2*pi*10.], [squares])

    args = {'on_1':.3, 'off_1':.325, 'ampl_1':1.,
            'on_2':.6, 'off_2':.625, 'ampl_2':1.}

### Solve with QuTiP

    tlist = np.linspace(0., 1., 201)

    # Need to set a max step or it won't see the switch on
    opts = qu.Options(max_step=.01)
    two_obj.mesolve(tlist, td=True, args=args, opts=opts, show_pbar=False)

    return two_obj, args

def plot(two_obj, args):

    t = two_obj.result.times
    states_t = two_obj.states_t()

    pop_1 = states_t[:,1,1].real
    coh_im = states_t[:,1,0].imag

    squares = lambda t, args: (tf.square_1(t, args) +
                               tf.square_2(t, args))

### Plot

    fig = plt.figure()
    ax = fig.add_subplot(3,1,1)
    ax.plot(t, squares(t, args), c='g', zorder=10)
    ax.set_ylabel(r'Rabi Frequency')
    ax = fig.add_subplot(3,1,2)
    ax.plot(t, pop_1)
    ax.set_ylabel(r'Excited Population')
    ax.set_ylim([0.,1.])
    ax = fig.add_subplot(3,1,3)
    ax.plot(t, coh_im, c='r')
    ax.set_xlabel(r'Time ($\tau$)')
    ax.set_ylabel(r'Imag Coherence')
    ax.set_ylim([-.5,.5])
    plt.draw()

def main():
    two_obj, args = solve()
    plot(two_obj, args)

### Save and show figures

    # for n in plt.get_fignums():
    #     plt.figure(n).savefig('png/'+savename+'_fig'+str(n)+'.png')
    #     plt.figure(n).savefig('pdf/'+savename+'_fig'+str(n)+'.pdf')    

    plt.show()

if __name__ == "__main__":
    status = main()
    sys.exit(status)