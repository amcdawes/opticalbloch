# -*- coding: utf-8 -*-

import os
import sys

import numpy as np 
import qutip as qu

"""
Contains the class for making optical bloch scan objects.

by Tommy Ogden <t.p.ogden@durham.ac.uk>
"""

class OBScan(object):

    def __init__(self, ob_obj, Delta_range, Delta_idx=0):

        self.ob_obj = ob_obj
        self.Delta_range = Delta_range
        self.Delta_idx = Delta_idx

        self.rho_Delta = [None]*len(self.Delta_range)
        self.result_Delta = [None]*len(self.Delta_range)

    def essolve(self, Deltas, tlist, rho0=None, recalc=True, savefile=None):

        savefile_exists = os.path.isfile(str(savefile) + '.qu')

        # If 1) we ask for it to be recalculated or 2) it *must* be calculated
        # because no savefile exists: do the steady state calc. Else load file.
        if (recalc or not savefile_exists):

            # Reset rho_Delta in case of multiple calcs
            self.rho_Delta = [None]*len(self.Delta_range)
            self.result_Delta = [None]*len(self.Delta_range)

            Deltas_i = list(Deltas) # Make copy, don't modify in place                

            Delta_steps = len(self.Delta_range)-1

            for i, Delta_i in enumerate(self.Delta_range):



                print("\rDelta: " + str(Delta_i) + ", " + str(i) + "/" +
                      str(Delta_steps))

                # Set the omega of the chosen beam to the current delta step.
                Deltas_i[self.Delta_idx] = Delta_i
                self.ob_obj.set_H_Delta(Deltas_i)

                try:
                    result = self.ob_obj.essolve(tlist, rho0=rho0)

                except ValueError:
                    print ("Failed to solve.")

                self.result_Delta[i] = result
                self.rho_Delta[i] = result.states[-1]

            # Only save the file if we have a place to save it
            if (savefile != None):
                qu.qsave((self.rho_Delta, self.result_Delta), savefile)

        # Otherwise load the steady state rho_v_delta from file
        else:
            (self.rho_Delta, self.result_Delta) = qu.qload(savefile)

        return self.rho_Delta, self.result_Delta

    def steadystate(self, Deltas, recalc=True, savefile=None):# , **kwargs):
        
        savefile_exists = os.path.isfile(str(savefile) + '.qu')

        # If 1) we ask for it to be recalculated or 2) it *must* be calculated
        # because no savefile exists: do the steady state calc. Else load file.
        if (recalc or not savefile_exists):

            # Reset rho_Delta in case of multiple calcs
            self.rho_Delta = [None]*len(self.Delta_range)
            Deltas_i = list(Deltas) # Make copy, don't modify in place

            Delta_steps = len(self.Delta_range)-1

            for i, Delta_i in enumerate(self.Delta_range):

                print("\rDelta: " + str(Delta_i) + ", " + str(i) + "/" +
                      str(Delta_steps))

                # Set the omega of the chosen beam to the current delta step.
                Deltas_i[self.Delta_idx] = Delta_i
                self.ob_obj.set_H_Delta(Deltas_i)

                try:
                    self.ob_obj.steadystate()
                except ValueError:
                    print ("Failed to find steadystate.")

                self.rho_Delta[i] = self.ob_obj.rho

            # Only save the file if we have a place to save it
            if (savefile != None):
                qu.qsave(self.rho_Delta, savefile)

        # Otherwise load the steady state rho_v_delta from file
        else:
            self.rho_Delta = qu.qload(savefile)

        return self.rho_Delta

    def steadystate_parfor(self, Deltas, num_cpus=None, recalc=True, 
                           savefile=None):

        savefile_exists = os.path.isfile(str(savefile) + '.qu')

        # If 1) we ask for it to be recalculated or 2) it *must* be calculated
        # because no savefile exists: do the steady state calc. Else load file.
        if (recalc or not savefile_exists):

            # Reset rho_Delta in case of multiple calcs
            self.rho_Delta = [None]*len(self.Delta_range)
            Deltas_i = list(Deltas) # Make copy, don't modify in place

            Delta_steps = len(self.Delta_range)-1

            self.rho_Delta = qu.parfor(steadystate_parfor_func, 
                                       self.Delta_range, 
                                       Delta_idx=self.Delta_idx, 
                                       Deltas=Deltas_i, ob_obj=self.ob_obj,
                                       num_cpus=num_cpus)

            # Only save the file if we have a place to save it
            if (savefile != None):
                qu.qsave(self.rho_Delta, savefile)

        # Otherwise load the steady state rho_v_delta from file
        else:
            self.rho_Delta = qu.qload(savefile)

        return self.rho_Delta

    def get_rho_Delta(self):

        rho_Delta = np.zeros((len(self.Delta_range), self.ob_obj.num_states, 
                              self.ob_obj.num_states), dtype=np.complex)

        for i, rho_i in enumerate(self.rho_Delta):
            rho_Delta[i] = self.rho_Delta[i].full()

        return rho_Delta

    def get_rho_Delta_element_sum(self, rho_list):

        rho_sum = np.zeros(len(self.Delta_range), dtype=np.complex)

        rho_Delta = self.get_rho_Delta()

        for i, Delta in enumerate(self.Delta_range):
            for j in rho_list:         
                rho_sum[i] += rho_Delta[i].item(j)

        return rho_sum

    def get_rho_Delta_element_sum_cutoff(self, rho_list, cutoff):

        rho_sum = np.zeros(len(self.Delta_range), dtype=np.complex)

        rho_Delta = self.get_rho_Delta()

        for j in rho_list:
            if (j[0] == j[1]): # If we're looking at a population
                pass
            else: # we're looking at a coherence
                omega_0 = (self.ob_obj.H_0[j[1],j[1]] - 
                           self.ob_obj.H_0[j[0],j[0]])
                for i, Delta in enumerate(self.Delta_range):
                    # rho_sum[i] += rho_Delta[i].item(j)
                    rho_sum[i] += (rho_Delta[i].item(j)*
                                   (abs(Delta - omega_0) < cutoff))

        return rho_sum

def steadystate_parfor_func(Delta, Delta_idx=0, Deltas=[], ob_obj=None):

    print("Delta: " + str(Delta) + "\r")

    Deltas[Delta_idx] = Delta
    
    ob_obj.set_H_Delta(Deltas)
    H = ob_obj.H_0 + ob_obj.H_Delta + ob_obj.H_I_sum()

    return qu.steadystate(H, ob_obj.c_ops)