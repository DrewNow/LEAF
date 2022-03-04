import sys
import pandas as pd
import pymatgen
#from pymatgen.analysis.local_env import *
from pymatgen.io.cif import CifParser
#from matminer.data_retrieval.retrieve_MP import MPDataRetrieval
#from matminer.featurizers.base import MultipleFeaturizer
from matminer.featurizers.site.fingerprint import OPSiteFingerprint, VoronoiFingerprint

#from ase.data import chemical_symbols

featurizer = OPSiteFingerprint() 
featurizer = VoronoiFingerprint()
#f = featurizer.feature_labels()
#print(f)

#atoms = chemical_symbols[1:3]
#print(atoms)


# get data

# from cifs
parser = CifParser("icsd_083501.cif")
structure = parser.get_structures()[0]


# MP retrieval
#mp = MPDataRetrieval(api_key="4r0hw9Mm0wjGgkSeIE").get_dataframe(criteria='Ta-S', properties=['structure','formula'])
#print(mp.head())
#print(mp.shape)
#for com in mp.formula.values:
#    print([el for el in com.keys()])

# species = [str(s).split()[-1] for s in structure]
# species = [''.join([a for a in s if a.isalpha()]) for s in species]
# print(len(species))

#mp = pd.DataFrame({'structure': [structure]})

features = []
for i in range(len(structure)): features.append(featurizer.featurize(structure, i))
for i in features: print(i)
sys.exit()

features = featurizer.featurize_dataframe(mp, col_id='structure',ignore_errors=True)
features = features.dropna()

print(features.shape)
print(features.head())
#for i in features.iloc[1]: print(i)
