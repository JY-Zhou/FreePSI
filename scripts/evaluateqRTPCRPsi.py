import os
import math
import sys
import json
import numpy as np
import scipy as scp
import scipy.stats as stats
from outliers import smirnov_grubbs as grubbs
from sklearn import metrics


def filter(i, j):
    if len(truePsi[i]) < 40:
        if mask[i][j] > 0:
            return True
    return False

def statFilter():
    reserveExon = 0
    for i in range(len(mask)):
        for j in range(len(mask[i])):
            if filter(i, j):
                reserveExon += 1
                print(truePsi[i][j], end = '\t')
                print(estPsi[i][j])
    print('+++ # spliced exons = \t' + str(reserveExon))

def globalRSME():
    estFlatPsi = []
    trueFlatPsi = []
    for i in range(len(estPsi)):
        for j in range(len(estPsi[i])):
            if filter(i, j):
                estFlatPsi.append(estPsi[i][j])
                trueFlatPsi.append(truePsi[i][j])
    #print('--- Global RMSE = ', end = '\t')
    print(np.sqrt(metrics.mean_squared_error(estFlatPsi, trueFlatPsi)), end = '\t')

def globalCorrelation():
    estFlatPsi = []
    trueFlatPsi = []
    for i in range(len(estPsi)):
        for j in range(len(estPsi[i])):
            if filter(i, j):
                estFlatPsi.append(estPsi[i][j])
                trueFlatPsi.append(truePsi[i][j])
    #print('--- Global Pearson correlation = ', end = '\t')
    print(stats.pearsonr(trueFlatPsi, estFlatPsi)[0], end = '\t')
    #print('--- Global Spearman correlation = ', end = '\t')
    print(stats.spearmanr(trueFlatPsi, estFlatPsi)[0], end = '\t')

eps = 1e-6

nargv = 1
truePsiPath = sys.argv[nargv];nargv += 1
estPsiPath = sys.argv[nargv];nargv += 1
maskPath = sys.argv[nargv];nargv += 1
outputPath = sys.argv[nargv];nargv += 1

truePsiFile = open(truePsiPath, 'r')
estPsiFile = open(estPsiPath, 'r')
maskFile = open(maskPath, 'r')
outputFile = open(outputPath, 'w')

truePsi = json.load(truePsiFile)
estPsi = json.load(estPsiFile)
mask = json.load(maskFile)
sys.stdout = outputFile

print("********** " + outputPath)
print("# exon\tRMSE\tPearson\tSpearman")
statFilter()
#print("======")
globalRSME()
#print("======")
globalCorrelation()
#print("======")
print("\n")
