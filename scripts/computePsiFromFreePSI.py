import json
import math
import numpy as np
import scipy.spatial.distance as dis
import scipy.cluster.hierarchy as hac
import scipy.stats as stats
import sys

inFile = open(sys.argv[1], 'r')
outFile = open(sys.argv[2], 'w')

psi = json.load(inFile)
newPsi = [[0.0 for e in range(len(psi[g]))] for g in range(len(psi))]
klist = []
ptplist = []
sdlist = []
skewlist = []
dlist = []

#Hierarchical clustering
for g in range(len(psi)):
    b = psi[g]
    newPsi[g] = b
    if len(b) > 1:
        maxb = max(b)
        a = np.array([b]).T
        b_dist = dis.pdist(a)
        b_link = hac.linkage(b_dist, method='average')
        for k in range(1, len(a)+1):
            clust = hac.cut_tree(b_link, k)
            Mu = []
            SD = []
            Skew = []
            PTP = []
            Dist = []
            for i in range(k):
                data = a[clust == i, ]
                Mu.append(np.mean(data))
                SD.append(np.std(data))
                Skew.append(stats.skew(data))
                PTP.append(np.ptp(data))
                dist = 0
                if len(data) > 1:
                    for p in range(len(data)):
                        for q in range(p+1, len(data)):
                            dist += abs(data[p] - data[q])
                    dist /= (len(data) * (len(data) - 1) / 2)
                Dist.append(dist)

            if np.max(SD) < 0.06 or np.mean(SD) < 0.05:
                klist.append(k)
                sdlist.append(np.max(SD))
                skewlist.append(np.max(Skew))
                ptplist.append(np.max(PTP))
                dlist.append(np.max(Dist))
                for e in range(len(a)):
                    newPsi[g][e] = Mu[clust[e, 0]]
                maxv = max(newPsi[g])
                if maxv > 0.05:
                    for e in range(len(a)):
                        newPsi[g][e] /= maxv
                else:
                    for e in range(len(a)):
                        newPsi[g][e] /= 5.0 
                break

json.dump(newPsi, outFile, indent = 4)
print("AVG cluster number = " + str(np.mean(klist)))
print("max max sd = " + str(np.max(sdlist)))
print("max max skew = " + str(np.max(skewlist)))
print("max max range = " + str(np.max(ptplist)))
print("max max dist = " + str(np.max(dlist)))
