import random
import json
import pandas 
import scipy.stats as stats
import sys

print('\n=== %s start... ===' % sys.argv[0])

nargv = 1
expPath = sys.argv[1];nargv += 1
profPath = sys.argv[nargv];nargv += 1
namePath = sys.argv[nargv];nargv += 1
psiPath = sys.argv[nargv];nargv += 1
tpmPath = sys.argv[nargv];nargv += 1

expFile = open(expPath, 'r')
profFile = open(profPath, 'r')
nameFile = open(namePath, 'r')
psiFile = open(psiPath, 'w')
tpmFile = open(tpmPath, 'w')

prof = json.load(profFile)
NG = prof['NG']
NE = prof['NE']
iso = prof['Iso']

nameMap = json.load(nameFile)

truePsi = [[0.0 for e in range(NE[g])] for g in range(NG)]
trueTpm = [0.0 for g in range(NG)]

trueIsoCount = [0 for g in range(NG)]

TPMlist = []
NameIndex = {}
for line in expFile:
    substr = line.strip().split('\t')
    txName = substr[1]
    tpm = float(substr[8]) / float(substr[3])

    NameIndex[txName] = len(TPMlist)
    TPMlist.append(tpm)

totTPM = sum(TPMlist)

for txName in NameIndex:
    i = NameIndex[txName]
    tpm = TPMlist[i] / totTPM * 1e6

    if tpm <= 1.0:
        continue

    isoId = nameMap[txName]
    geneId = isoId.split('Iso')[0]
    geneNum = int(geneId.replace('Gene', ''))

    trueTpm[geneNum] += tpm
    trueIsoCount[geneNum] += 1

    for e in iso[geneId][isoId]:
        truePsi[geneNum][e] += tpm

for g in range(NG):
    if trueIsoCount[g] > 0:
        for e in range(NE[g]):
            truePsi[g][e] /= trueTpm[g]

json.dump(truePsi, psiFile, indent = 4)
json.dump(trueTpm, tpmFile, indent = 4)

print('=== %s finished! ===\n' % sys.argv[0])
