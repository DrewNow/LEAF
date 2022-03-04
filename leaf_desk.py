""" Local environment-based atomic features
"""

import os
import sys
import numpy as np
import pandas as pd
from ase.data import chemical_symbols
from pymatgen.analysis.local_env import *
from pymatgen.io.cif import CifParser
#import pickle

try:
   # Python>=3.7
   from matminer.featurizers.site.fingerprint import OPSiteFingerprint, VoronoiFingerprint
except: pass
#  # Python 3.6
#  from matminer.featurizers.site import OPSiteFingerprint, VoronoiFingerprint

class Leaf:
    """
    Creates a matrix of one-hot encoded atomic representation
    Values: lists of local environment features (lostops, voronoi tes.)
    Keys:   atomic elements
    """

    def __init__(self, featurizer=OPSiteFingerprint()):
        self.featurizer = featurizer

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

    def average_features_cifs(self, list_cifs):
        """ process a list of cifs  (e.g. 1.dat)
        average features over a number of occurences of the elements in icsd
        write results into a dictionary """

        cifs, order = self.readfile(list_cifs)
        leaf = {atom: np.zeros(37) for atom in chemical_symbols}
        nleaf = {atom: 0 for atom in chemical_symbols}

        for cif in cifs:
            # print(f"Processing {cif}")
            try:
                structure = CifParser(cif).get_structures()[0]
            except:
                continue
            features = []
            for i, s in enumerate(structure):
                species = self.get_species(s)
                for element in species:
                    leaf[element] += self.featurizer.featurize(structure, i)
                    nleaf[element] += 1

        df = pd.DataFrame(leaf)
        df.to_pickle(f'OPleaf_{order}.pickle')
        df = pd.DataFrame(nleaf, index=[order])
        df.to_pickle(f'nleaf_{order}.pickle')

if __name__ == "__main__":
    fname = sys.argv[1]
    leaf = Leaf()
    leaf.average_features_cifs(fname)
