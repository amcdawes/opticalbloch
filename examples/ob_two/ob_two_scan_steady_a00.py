# -*- coding: utf-8 -*-

import os
import sys
import time

import numpy as np
import qutip as qu

import matplotlib.pyplot as plt

import ob_two
from opticalbloch import ob_scan

pi = np.pi

base = os.path.basename(__file__)
savename = os.path.splitext(base)[0]

def scan_steady():

    # Detuning range
    Delta_range = np.linspace(-2*pi*10., 2*pi*10., 51)

### TwoOB

    two_obj = ob_two.OBTwo(gammas=[2*pi*1.])

    two_obj.set_H_Delta([2*pi*0.])
    two_obj.set_H_Omega([2*pi*0.001])

    two_scan = ob_scan.OBScan(two_obj, Delta_range)

    two_scan.steadystate(Deltas=[None])

    return two_scan

def plot(scan_obj):

    from scipy.ndimage import zoom
    from scipy.interpolate import UnivariateSpline

    import analytic

### Data

    Delta_range = scan_obj.Delta_range/(2*pi)

    pop_1 = scan_obj.get_rho_Delta()[:,1,1].real

    coh_im = scan_obj.get_rho_Delta()[:,0,1].imag
    coh_re = -scan_obj.get_rho_Delta()[:,0,1].real

    Omega = 2*pi*.001
    Gamma = 2*pi*1.

    I_over_I_sat = analytic.I_over_I_sat(Omega, Gamma)
    fwhm = analytic.fwhm(Omega, Gamma)/(2.*np.pi)

    # Interpolate
    z = 4
    Delta_range = zoom(Delta_range, z)

    pop_1 = zoom(pop_1, z)

    coh_im = zoom(coh_im, z)
    coh_re = zoom(coh_re, z)

### Plot

### Fig 1: Populations and Coherences

    fig = plt.figure()#, figsize=(12.6, 9.4))

    file_text = fig.text(0.01, 0.02, savename,
                           fontsize='smaller', color='darkgrey')

    # Plot populations

    ax = fig.add_subplot(211)

    ax.plot(zoom(Delta_range,z), zoom(pop_1,z),
                         label=r'$\rho_{11}$',
                         c='r',
                         clip_on=False, zorder=10)

    ax.set_ylabel(r'Population')

    I_text = fig.text(0.1, 0.9,
                        r'$I/I_{sat} = %(I)0.3f$' % {'I': I_over_I_sat},
                           fontsize='medium', color='black')
    fwhm_text = fig.text(0.1, 0.45,
                           r'$\mathrm{FWHM} = %(fwhm)0.3f$' % {'fwhm': fwhm},
                           fontsize='medium', color='black')

    # Plot coherences

    ax = fig.add_subplot(212)

    ax.plot(zoom(Delta_range,z), zoom(coh_im,z),
                         label=r'$\mathrm{Im}{\rho_{01}}$',
                         c='g',
                         clip_on=False, zorder=10)
    ax.plot(zoom(Delta_range,z), zoom(coh_re,z),
                         label=r'$\mathrm{Re}{\rho_{01}}$',
                         c='b',
                         clip_on=False, zorder=10)

    ax.set_xlabel(r'Detuning ($\Gamma$)')
    ax.set_ylabel('Coherence')

    for ax in fig.axes:

        ax.set_xlim([Delta_range[0], Delta_range[-1]])
        ax.set_xticks(np.linspace(Delta_range[0],Delta_range[-1],4+1))

        leg = ax.legend()
        leg.get_frame().set_alpha(0.5)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

### FWHM Marker

    half_max = np.max(coh_im)/2

    # Create a spline of Delta and coh_im-half_max
    spline = UnivariateSpline(Delta_range, coh_im-half_max, s=0)
    r1, r2 = spline.roots() # find the roots

    # Draw line at FWHM
    ax.hlines(y=half_max, xmin=r1, xmax=r2,
                color='g', linestyle='dotted', linewidth=1.0)

    # Annotate with arrow and FWHM value
    connectionstyle = "angle3,angleA=0,angleB=110"
    ax.annotate('FWHM: ' + '%0.3f'%(r2 - r1), xy=((r2+r1)/2, half_max),
                  xycoords='data',
                  xytext=(-120, 50), textcoords='offset points',
                  zorder=20,
                  arrowprops=dict(arrowstyle="-|>",
                                  color="0.5",
                                  shrinkA=5, shrinkB=5,
                                  patchA=None,
                                  patchB=None,
                                  relpos=(1., .5),
                                  connectionstyle=connectionstyle))

def main():

    start = time.time()
    two_scan = scan_steady()
    print("Elapsed time: {:.2f}").format(time.time()-start)

    plot(two_scan)

### Save and show figures

    # for fig_num in plt.get_fignums():
    #     plt.figure(fig_num).savefig('png/'+savename+'_fig'+str(fig_num)+'.png')
    #     plt.figure(fig_num).savefig('pdf/'+savename+'_fig'+str(fig_num)+'.pdf')

    plt.show()

if __name__ == "__main__":
    status = main()
    sys.exit(status)