""" Local environment-based atomic features
"""

import os
import sys
import numpy as np
import pandas as pd
from ase.data import chemical_symbols
from pymatgen.analysis.local_env import *
from pymatgen.io.cif import CifParser
import pickle

try:
   # Python>=3.7
   from matminer.featurizers.site.fingerprint import OPSiteFingerprint, VoronoiFingerprint
except: 
   # Python 3.6
   from matminer.featurizers.site import OPSiteFingerprint, VoronoiFingerprint


class LEAF:
    """
    Creates a matrix of one-hot encoded atomic representation
    Values: lists of local environment features (lostops, voronoi tes.)
    Keys:   atomic elements
    """

    def __init__(self, leaf=None, featurizer=OPSiteFingerprint, icsd_dic=None):
        if leaf:
            self.leaf = leaf
        else:
            self.leaf = {atom: [] for atom in chemical_symbols}
            self.nleaf = {atom: [] for atom in chemical_symbols}
        self.featurizer = featurizer
        self.icsd = self.parse_icsd(icsd_dic)

    @staticmethod
    def get_species(structure):
        species = [str(s).split()[-1] for s in structure]
        return [''.join([a for a in s if a.isalpha()]) for s in species]

    def get_features(self, structure):
        features = []
        for i in range(len(structure)):
            features.append(self.featurizer.featurize(structure, i))
        return features

    @staticmethod
    def read_file(list_cifs):
        """ read file with list of cifs e.g. 1.dat """

        cifs = open(list_cifs, 'r').readlines()
        cifs = [i.strip() for i in cifs]
        order = list_cifs.split('.')[0]
        return cifs, order

    def average_features_cifs(self, list_cifs):
        """ process a list of cifs  (e.g. 1.dat)
        average features over a number of occurences of the elements in icsd
        write results into a dictionary """

        cifs, order = readfile(list_cifs)
        leaf = {atom: np.zeros(37) for atom in chemical_symbols}
        nleaf = {atom: 0 for atom in chemical_symbols}

        for cif in cifs:
            structure = CifParser(fname).get_structures()[0]
            species = self.get_species(structure)
            features = []
            for i in range(len(structure)): 
                leaf[species[i]] += self.featurizer.featurize(structure, i) 
                nleaf[species[i]] += 1

        pd.DataFrame(leaf).to_pickle(f'leaf_{order}.pickle') 
        pd.DataFrame(nleaf).to_pickle(f'nleaf_{order}.pickle') 
             

    @staticmethod 
    def parse_icsd(icsd_dic=None, f='./ICSD/full_icsd'):
       """ load icsd dictionary 
           or create one (process all cif files) """

        if icsd_dic:
           return icds_dic

        structures = []
        for path, dirs, files in os.walk(f):
            for i in files:
                fname = os.path.join(path, i)
                structure = CifParser(fname).get_structures()[0]
                structures.append([structure])

        icsd_dic = pd.DataFrame({'structure': structures})
        icsd_dic.to_pickle('icsd_structures.pickle')
        return icsd_dic


