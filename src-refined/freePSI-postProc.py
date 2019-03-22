#! /bin/env python3

import json
import math
import numpy as np
import scipy.spatial.distance as dis
import scipy.cluster.hierarchy as hac
import sys

inFile = open(sys.argv[1], 'r')
outFile = open(sys.argv[2], 'w')

psi = json.load(inFile)
newPsi = [[0.0 for e in range(len(psi[g]))] for g in range(len(psi))]

#Hierarchical clustering
for g in range(len(psi)):
    b = psi[g]
    newPsi[g] = b
    if len(b) > 1:
        a = np.array([b]).T
        b_dist = dis.pdist(a)
        b_link = hac.linkage(b_dist)
        for k in range(1, len(a) + 1):
            clust = hac.cut_tree(b_link, k)
            avgSD = []
            mu = []
            for i in range(k):
                mu.append(np.mean(a[clust == i, ]))
                avgSD.append(np.std(a[clust == i, ]))
            if np.max(avgSD) < 0.05:
                for e in range(len(a)):
                    newPsi[g][e] = mu[clust[e, 0]]
                maxv = max(newPsi[g])
                if maxv > 0.05:
                    for e in range(len(a)):
                        newPsi[g][e] /= maxv
                else:
                    for e in range(len(a)):
                        newPsi[g][e] = 0
                break

json.dump(newPsi, outFile, indent=4)
