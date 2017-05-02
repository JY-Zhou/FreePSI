import os
import math
import sys
import json
import numpy as np
import scipy as scp
import scipy.stats as stats
from sklearn import metrics

TPMFILTER = 1

def filter(i, j):
    if len(truePsi[i]) < 40:
        if trueTpm[i] >= TPMFILTER and mask[i][j] > 0:
            return True
    return False

def statFilter():
    reserveGene = 0
    reserveExon = 0
    for i in range(len(mask)):
        add = False
        for j in range(len(mask[i])):
            if filter(i, j):
                reserveExon += 1
                add = True
        if add:
            reserveGene += 1
    print('+++ # spliced genes = \t' + str(reserveGene))
    print('+++ # spliced exons = \t' + str(reserveExon))

def globalRSME():
    estFlatPsi = []
    trueFlatPsi = []
    for i in range(len(estPsi)):
        for j in range(len(estPsi[i])):
            if filter(i, j):
                estFlatPsi.append(estPsi[i][j])
                trueFlatPsi.append(truePsi[i][j])
    print('--- Global RMSE = ', end = '\t')
    print(np.sqrt(metrics.mean_squared_error(estFlatPsi, trueFlatPsi)))

def globalCorrelation():
    estFlatPsi = []
    trueFlatPsi = []
    for i in range(len(estPsi)):
        for j in range(len(estPsi[i])):
            if filter(i, j):
                estFlatPsi.append(estPsi[i][j])
                trueFlatPsi.append(truePsi[i][j])
    print('--- Global correlation = ', end = '\t')
    print(stats.pearsonr(trueFlatPsi, estFlatPsi)[0])

eps = 1e-6

nargv = 1
truePsiPath = sys.argv[nargv];nargv += 1
trueTpmPath = sys.argv[nargv];nargv += 1
estPsiPath = sys.argv[nargv];nargv += 1
maskPath = sys.argv[nargv];nargv += 1
outputPath = sys.argv[nargv];nargv += 1

truePsiFile = open(truePsiPath, 'r')
trueTpmFile = open(trueTpmPath, 'r')
estPsiFile = open(estPsiPath, 'r')
maskFile = open(maskPath, 'r')
outputFile = open(outputPath, 'w')

truePsi = json.load(truePsiFile)
trueTpm = json.load(trueTpmFile)
estPsi = json.load(estPsiFile)
mask = json.load(maskFile)

sys.stdout = outputFile

print("********** " + outputPath)
statFilter()
print("======")
globalRSME()
print("======")
globalCorrelation()
print("======")
print("\n")
