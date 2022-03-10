""" Local environment-based atomic features
"""

import os
import sys
import numpy as np
import pandas as pd
from ase.data import chemical_symbols
from pymatgen.analysis.local_env import *
from pymatgen.io.cif import CifParser
from matminer.featurizers.site.fingerprint import OPSiteFingerprint, VoronoiFingerprint

class Leaf:
    """
    Creates a matrix of one-hot encoded atomic representation
    Values: lists of local environment features (lostops, voronoi tes.)
    Keys:   atomic elements
    """

    def __init__(self, featurizer=OPSiteFingerprint()):
        self.featurizer = featurizer
        self.onehot = {atom: {} for atom in chemical_symbols}
        self.features_names = self.get_features_names()

    def get_features_names(self):
        lostops = [i for i in OPSiteFingerprint().feature_labels()]
        voronoi = [i for i in VoronoiFingerprint().feature_labels()]
        voronoi = [v for v in voronoi if 'std' not in v]
        return lostops + voronoi

    @staticmethod
    def get_species(site):
        species = str(site.species).split()
        return [''.join([a for a in s if a.isalpha()]) for s in species \
                if ''.join([a for a in s if a.isalpha()]) in chemical_symbols]

    @staticmethod
    def readfile(list_cifs):
        """ read file with list of cifs e.g. 1.dat """

        cifs = open(list_cifs, 'r').readlines()
        cifs = [i.strip() for i in cifs]
        order = list_cifs.split('.')[0]
        return cifs, order

    @staticmethod
    def select_features(i, s):
        """
        from all voronoi features select:
        'Voro_vol_sum', 'Voro_area_sum', 'Voro_vol_mean', 'Voro_vol_minimum',
        'Voro_vol_maximum', 'Voro_area_mean', 'Voro_area_minimum', 'Voro_area_maximum',
        'Voro_dist_mean', 'Voro_dist_minimum', 'Voro_dist_maximum'
        """
        l = list(VoronoiFingerprint().featurize(s, i))
        b = l[16:19] + l[20:23] + l[24:27] + l[28:31]
        return np.array(b)

    def average_features_cifs(self, list_cifs):
        """ process a list of cifs  (e.g. 1.dat)
        average features over a number of occurences of the elements in icsd
        write results into a dictionary """

        cifs, order = self.readfile(list_cifs)
        # structures dict
        icsd = {'composition':[],'structure':[]}

        # lostops
        # leaf = {atom: np.zeros(37) for atom in chemical_symbols}
        # nleaf = {atom: 0 for atom in chemical_symbols}

        # voronoi
        leaf = {atom: np.zeros(11) for atom in chemical_symbols}
        nleaf = {atom: 0 for atom in chemical_symbols}
        for cif in cifs:
            try:
                structure = CifParser(cif).get_structures()[0]
            except:
                continue
            icsd['composition'].append(str(structure.composition))
            icsd['structure'].append(structure)
            features = []
            for i, s in enumerate(structure):
                species = self.get_species(s)
                for element in species:
                    # lostops
                    # leaf[element] += self.featurizer.featurize(structure, i)
                    # nleaf[element] += 1

                    # voronoi
                    try:
                        features = self.select_features(self.featurizer.featurize(structure, i))
                        leaf[element] += features
                        nlean[element] += 1
                    except:
                        break

        df = pd.DataFrame(leaf)
        df.to_pickle(f'Vleaf_{order}.pickle')
        df = pd.DataFrame(icsd)
        df.to_pickle(f'icsd_{order}.pickle')
        df = pd.DataFrame(nleaf, index=[order])
        df.to_pickle(f'vnleaf_{order}.pickle')

    def get_features(self, list_cifs):
        cifs, order = self.readfile(list_cifs)
        features = []
        elements = []
        for cif in cifs:
            try:
                structure = CifParser(cif).get_structures()[0]
            except:
                continue
            print(structure.composition)
            for i, s in enumerate(structure):
                species = self.get_species(s)
                for element in species:
                    vfeat = []
                    try:
                        vfeat = [round(f,3) for i in select_features(i,structure)]
                    except:
                        pass

                    feats = [round(f,3) for f in self.featurizer.featurize(structure, i)] + \
                            vfeat
                    features.append(feats)
                    elements.append(element)
        df = pd.DataFrame({'elements': elements, 'features': features})
        df.sort_values(by=['elements'])
        df.to_pickle(f'test_features_{order}.pickle')
    def expand_onehot(self, element, features):
        """ For all features 'feature_name' with values 'v'
        create one-hot columns 'feature_name_v'
        fill the number of occurences into self.onehot
        """
        for i,v in enumerate(features):
            feature = self.features_names + '_' + str(v)
            if feature in self.onehot[element]:
                self.onehot[element][feature] += 1
            else:
                for atom in chemical_elements:
                    self.onehot[atom][feature] = 0

     def sort_features(self):
         """ sort by features names
         return dictionary[atoms] = 'features values' """
         for atom in self.onehot:
             self.onehot[atom] = {k: v for k,v in sorted(self.onehot[atom].items(), key=lambda x: x[0])}

         return {a: np.array(self.onehot[a].values()) for a in self.onehot}



if __name__ == "__main__":
    fname = sys.argv[1]
    #leaf = Leaf(VoronoiFingerprint())
    leaf = Leaf()
    #leaf.average_features_cifs(fname)
    leaf.get_features(fname)
