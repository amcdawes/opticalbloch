import numpy as np

class Atom1e(object):

    def __init__(self, element, isotope, I, S):

        self.element = element
        self.isotope = isotope
        self.I = I
        self.S = S
        self.nL_levels = []

    def __repr__(self):

        return "<Atom1e :: %s>" % self.__dict__

    def add_nL_level(self, nL_level):

        self.nL_levels.append(nL_level)

    def get_mF_list(self):

        mF_list = []

        for nL_level in self.nL_levels:

            item_n = nL_level.n
            item_L = nL_level.L

            for J_level in nL_level.J_levels:

                item_J = J_level.J

                for F_level in J_level.F_levels:

                    item_F = F_level.F

                    for mF_level in F_level.mF_levels:

                        item_mF = mF_level.mF
                        item_energy = mF_level.energy

                        mF_dict = {'I':self.I, 'S':self.S, 'n':item_n,
                                   'L':item_L, 'J':item_J, 'F':item_F,
                                   'mF':item_mF, 'energy':item_energy}

                        mF_list.append(mF_dict)

        return mF_list

    def get_num_mF_levels(self):

        return len(self.get_mF_list())

    def get_mF_energies(self):

        mF_list = self.get_mF_list()

        mF_energies = []

        for mF_level in mF_list:

            mF_energies.append(mF_level['energy'])

        return mF_energies

class LevelNL(object):

    def __init__(self, n, L, energy):

        self.n = n
        self.L = L
        self.J_levels = []
        self.energy = energy

    def __repr__(self):

        return "<LevelNL :: %s>" % self.__dict__

    def build_J_levels(self, S, I, J_energies, F_energies, mF_energies):

        for i, J in enumerate(self.get_J_range(S)):

            self.add_J_level(LevelJ(J, J_energies[i]))

            # I think this check should be outside loop
            if (mF_energies == None):

                self.J_levels[i].build_F_levels(I, F_energies[i])

            elif (len(mF_energies) == self.get_num_J_levels()):

                 self.J_levels[i].build_F_levels(I, F_energies[i], mF_energies[i])

            else:

                pass # should raise an error here

    def add_J_level(self, J_level):

        self.J_levels.append(J_level)

    def get_J_range(self, S):

        return np.arange(abs(self.L-S),self.L+S+1)

    def get_num_J_levels(self):

        return self.get_J_range().size


class LevelJ(object):

    def __init__(self, J, energy):

        self.J = J
        self.F_levels = []
        self.energy = energy

    def __repr__(self):

        return "<LevelJ :: %s>" % self.__dict__

    def build_F_levels(self, I, F_energies, mF_energies=None):

        for i, F in enumerate(self.get_F_range(I)):

            self.add_F_level(LevelF(F, F_energies[i]))

            # I think this check should be outside loop
            if (mF_energies == None):

                self.F_levels[i].build_mF_levels()

            elif (len(mF_energies) == self.get_num_F_levels(I)):

                self.F_levels[i].build_mF_levels(mF_energies[i])

            else:

                pass # should raise an error here

    def add_F_level(self, F_level):

        self.F_levels.append(F_level)

    def get_F_range(self, I):

        return np.arange(abs(self.J-I),self.J+I+1)

    def get_num_F_levels(self, I):

        return self.get_F_range(I).size


class LevelF(object):

    def __init__(self, F, energy):

        self.F = F
        self.mF_levels = []
        self.energy = energy

    def __repr__(self):

        return "<LevelF :: %s>" % self.__dict__

    def build_mF_levels(self, mF_energies=None):

        # This could be improved. Use the if to make a np array, then do the
        # for loop.

        # if (len(mF_energies) == 0):

        if (mF_energies == None):

            for i, mF in enumerate(self.get_mF_range()):

                self.add_mF_level(LevelMF(mF, self.energy))

        elif (len(mF_energies) == self.get_num_mF_levels()):

            for i, mF in enumerate(self.get_mF_range()):

                self.add_mF_level(LevelMF(mF, mF_energies[i]))

        else:
            pass # Should raise an error here

    def add_mF_level(self, mF_level):

        self.mF_levels.append(mF_level)

    def get_mF_range(self):

        return np.arange(-self.F,self.F+1)

    def get_num_mF_levels(self):

        return self.get_mF_range().size


class LevelMF(object):

    def __init__(self, mF, energy):

        self.mF = mF
        self.energy = energy

    def __repr__(self):

        return "<LevelMF :: %s>" % self.__dict__