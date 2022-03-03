""" Local environment-based atomic features
"""

import os
import sys
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
        self.featurizer = featurizer
        self.icsd = self.parse_icsd(icsd_dic)

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
