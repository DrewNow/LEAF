import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.metrics.pairwise import cosine_similarity 

#phase  = get_vectors(elements, vectors) # arrays
#attention = normalize(phase @ phase.T)*0.74 # matrix

def similarity(a,b=None):
    a = a.values.reshape(1,-1)
    if any(b):
        b = b.values.reshape(1,-1)
    return cosine_similarity(a,b)[0][0]

def plot_map(matrix, elements):
    fig, ax = plt.subplots()
    axp = ax.pcolormesh(np.arange(len(matrix[0])), np.arange(len(matrix[0])), matrix, shading='auto', cmap='YlOrRd')#'magma')
    cb = plt.colorbar(axp,ax=[ax],location='right')
    #ax.yaxis.tick_right()
    #ax.set_yticks(ticks=[i for i in range(len(elements))])
    #ax.set_yticklabels(labels = elements)
    #ax.set_xticks(ticks=[i for i in range(len(elements))])
    #ax.set_xticklabels(labels = elements)
    #ax.tick_params(axis=u'both', which=u'both',length=0)
    #plt.xlabel('Atomic vector dimension', fontsize=13)
    plt.show()

if __name__ == '__main__':

    m = pd.read_pickle('LEAF_L.pickle')
    m = m[m.columns[1:]] # remove X
    m = m.dropna(axis=1) # remove NaN
    c = np.zeros((m.shape[1], m.shape[1]))

    #print(m.columns)

    #sys.exit()
    for i, a in enumerate(m.columns):
        for j, b in enumerate(m.columns):
            c[i, j] = similarity(m[a], m[b])

    plot_map(c, m.columns.values)
