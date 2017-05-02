import random
import json
import pandas 
import scipy.stats as stats
import sys

print('\n=== %s start... ===' % sys.argv[0])

nargv = 1
isoTpmPath = sys.argv[nargv];nargv += 1
profPath = sys.argv[nargv];nargv += 1
namePath = sys.argv[nargv];nargv += 1
psiPath = sys.argv[nargv];nargv += 1
geneTpmPath = sys.argv[nargv];nargv += 1

isoTpmFile = open(isoTpmPath, 'r')
profFile = open(profPath, 'r')
nameFile = open(namePath, 'r')
psiFile = open(psiPath, 'w')
geneTpmFile = open(geneTpmPath, 'w')

prof = json.load(profFile)
NG = prof['NG']
NE = prof['NE']
iso = prof['Iso']

nameMap = json.load(nameFile)

estPsi = [[0.0 for e in range(NE[g])] for g in range(NG)]
estTot = [0.0 for g in range(NG)]

k = 0
for line in isoTpmFile:
    k += 1
    if k == 1:
        continue
    substr = line.strip().split('\t')
    txName = substr[0]
    tpm = float(substr[3])

    if not txName in nameMap:
        continue

    isoId = nameMap[txName]
    geneId = isoId.split('Iso')[0]
    geneNum = int(geneId.replace('Gene', ''))

    estTot[geneNum] += tpm
    for e in iso[geneId][isoId]:
        estPsi[geneNum][e] += tpm

for g in range(NG):
    for e in range(NE[g]):
        if estTot[g] > 0:
            estPsi[g][e] /= estTot[g]
        else:
            estPsi[g][e] = 0.0

json.dump(estPsi, psiFile, indent = 4)
json.dump(estTot, geneTpmFile, indent = 4)

print('=== %s finished! ===\n' % sys.argv[0])
