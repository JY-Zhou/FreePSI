import os
import math
import sys
import json
import numpy as np
import scipy as scp
import scipy.stats as stats
from outliers import smirnov_grubbs as grubbs
from sklearn import metrics

def procGeneFamily():
    geneList = {}
    for line in geneMapFile:
        substr = line.split('\t')
        isoName = substr[3]
        geneName = substr[4]
        if not geneName in geneList:
            geneList[geneName] = []
        geneList[geneName].append(isoName)

    geneId = {}
    for x in geneList:
        geneId[x] = []
        for y in geneList[x]:
            if y in nameMap:
                z = int(nameMap[y].split('Iso')[0].replace('Gene', ''))
                if not z in geneId[x]:
                    geneId[x].append(z)

    family = {}
    for line in geneFamilyFile:
        substr = line.split('\t')
        if not substr[-1] in family:
            family[substr[-1]] = []
        family[substr[-1]].append((substr[1], substr[8], substr[9], substr[10]))

    familyList = []
    for x in family:
        familyList.append((len(family[x]), family[x]))

    for x in familyList:
        if x[0] > 20 and len(x[1][0][2]) > 0:
            ids = []
            s = 0
            for y in x[1]:
                if y[0] in geneId:
                    ids.extend(geneId[y[0]])
                    s += len(geneList[y[0]])
            #info = '>>> Gene_number=%d\tFamily_name=%s\tDescription=%s\n++> Mapped_gene_number=%d\tTotal_isoforms=%d\t' % (x[0], x[1][0][2], x[1][0][3], len(ids), s)
            info = '%d\t%s\t%s\t%d\t%d\t' % (x[0], x[1][0][2], x[1][0][3], len(ids), s)

            if len(ids) > 20 and s > 2*len(ids):
                evaluateList.append((info, ids))

TPMFILTER = 10
def filter(i):
    if not i in fam[1]:
        return False
    if len(truePsi[i]) < 40:
        if trueTpm[i] >= TPMFILTER:
            return True
    return False

def statFilter():
    print(fam[0], end = '')
    #print("+++ Filter: true TPM >= " + str(TPMFILTER))
    reserveGene = 0
    for i in range(len(estPsi)):
        if filter(i):
            reserveGene += 1
    #print('+++ # genes passed filter = \t' + str(reserveGene))
    #print('+++ # genes failed (low expression level) = \t' + str(len(fam[0]) - reserveGene))
    print(reserveGene, end = '\t')

def globalCorrelation():
    estFlatPsi = []
    trueFlatPsi = []
    totreads = 0
    totmulti = 0
    for i in range(len(estPsi)):
        if filter(i):
            estFlatPsi.extend(estPsi[i])
            trueFlatPsi.extend(truePsi[i])
            totreads += readcov[i]
            totmulti += multicov[i]
    #print('--- Global correlation = ', end = '\t')
    print(stats.pearsonr(trueFlatPsi, estFlatPsi)[0], end = '\t')
    print(totreads, end = '\t')
    print(totmulti)

def genelevel():
    geneList = {}
    for line in geneMapFile:
        substr = line.split('\t')
        isoName = substr[3]
        geneName = substr[4]
        if not geneName in geneList:
            geneList[geneName] = []
        geneList[geneName].append(isoName)

    geneId = {}
    for x in geneList:
        geneId[x] = []
        for y in geneList[x]:
            if y in nameMap:
                z = int(nameMap[y].split('Iso')[0].replace('Gene', ''))
                if not z in geneId[x]:
                    geneId[x].append(z)
    for i in range(len(estPsi)):
        if trueTpm[i] > TPMFILTER:
            print(i)


eps = 1e-6

nargv = 1
geneFamilyPath = sys.argv[nargv];nargv += 1
readsDisPath = sys.argv[nargv];nargv += 1
geneMapPath = sys.argv[nargv];nargv += 1
nameMapPath = sys.argv[nargv];nargv += 1
truePsiPath = sys.argv[nargv];nargv += 1
trueTpmPath = sys.argv[nargv];nargv += 1
estPsiPath = sys.argv[nargv];nargv += 1
estTpmPath = sys.argv[nargv];nargv += 1
outputPath = sys.argv[nargv];nargv += 1

geneFamilyFile = open(geneFamilyPath, 'r')
readsDisFile = open(readsDisPath, 'r')
geneMapFile = open(geneMapPath, 'r')
nameMapFile = open(nameMapPath, 'r')
truePsiFile = open(truePsiPath, 'r')
trueTpmFile = open(trueTpmPath, 'r')
estPsiFile = open(estPsiPath, 'r')
estTpmFile = open(estTpmPath, 'r')
outputFile = open(outputPath, 'w')

nameMap = json.load(nameMapFile)
truePsi = json.load(truePsiFile)
trueTpm = json.load(trueTpmFile)
estPsi = json.load(estPsiFile)
estTpm = json.load(estTpmFile)

readcov = []
multicov = []
for line in readsDisFile:
    substr = line.split('\t')
    readcov.append(float(substr[1]))
    multicov.append(float(substr[2]))

sys.stdout = outputFile

evaluateList = []
procGeneFamily()
evaluateList = sorted(evaluateList, key = lambda x: x[0])

print(outputPath.split('/')[-1])
print('Total reads in experiment = %d, total multireads = %d' % (sum(readcov), sum(multicov)))
print('Gene_Number\tFamily_name\tDescription\tMapped_Gene_Number\tTotal_isoform\tPassed_Gene_Number\tCorrelation\t# reads\t# multi-reads')

for fam in evaluateList:
    statFilter()
    globalCorrelation()
print()

#genelevel()
