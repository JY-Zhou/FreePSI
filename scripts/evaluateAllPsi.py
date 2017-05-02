import os
import math
import sys
import json
import numpy as np
import scipy as scp
import scipy.stats as stats
from sklearn import metrics

TPMFILTER = 10

def filter(i):
    if len(truePsi[i]) < 40:
        if trueTpm[i] >= TPMFILTER:
            return True
    return False

def statFilter():
    print("+++ Filter: true TPM >= " + str(TPMFILTER))
    reserveGene = 0
    for i in range(len(estPsi)):
        if filter(i):
            reserveGene += 1
    print('+++ # genes passed filter = \t' + str(reserveGene))

def globalRSME():
    estFlatPsi = []
    trueFlatPsi = []
    for i in range(len(estPsi)):
        if filter(i):
            estFlatPsi.extend(estPsi[i])
            trueFlatPsi.extend(truePsi[i])
    print('--- Global RMSE = ', end = '\t')
    print(np.sqrt(metrics.mean_squared_error(estFlatPsi, trueFlatPsi)))

def localRSMEMean():
    rmse = []
    for i in range(len(truePsi)):
        if filter(i):
            rmse.append(np.sqrt(metrics.mean_squared_error(truePsi[i], estPsi[i])))
    print('--- Local RMSE = ', end = '\t')
    print(np.mean(rmse))

def globalCorrelation():
    estFlatPsi = []
    trueFlatPsi = []
    for i in range(len(estPsi)):
        if filter(i):
            estFlatPsi.extend(estPsi[i])
            trueFlatPsi.extend(truePsi[i])
    print('--- Global correlation = ', end = '\t')
    print(stats.pearsonr(trueFlatPsi, estFlatPsi)[0])

def localCorrelationMean():
    cor = []
    for i in range(len(truePsi)):
        if filter(i):
            cor.append(stats.pearsonr(truePsi[i], estPsi[i])[0])
    print('--- Local correlation = ', end = '\t')
    print(np.mean(cor))

def splicedCorrelation():
    estSplicedLabel = [0 for i in range(len(estPsi))]
    estSplicedScore = [0 for i in range(len(estPsi))]
    for i in range(len(estPsi)):
        if filter(i):
            #if np.std(estPsi[i]) > 0.01:
            if np.std(estPsi[i]) > 0.05:
                estSplicedLabel[i] = 1
            estSplicedScore[i] = np.std(estPsi[i])
    
    trueSplicedLabel = [0 for i in range(len(truePsi))]
    for i in range(len(truePsi)):
        if filter(i):
            #if abs(min(truePsi[i]) - max(truePsi[i])) > 0.05:
            if np.std(truePsi[i]) > 0.05:
                trueSplicedLabel[i] = 1

    print('--- Precision score of spliced genes = ', end = '\t')
    print(metrics.precision_score(trueSplicedLabel, estSplicedLabel))
    print('--- Recall score of spliced genes = ', end = '\t')
    print(metrics.recall_score(trueSplicedLabel, estSplicedLabel))
    print('--- F1 score of spliced genes = ', end = '\t')
    print(metrics.f1_score(trueSplicedLabel, estSplicedLabel))
    print('!*!*!*!\n--- ROC-AUC score of spliced genes = ', end = '\t')
    print(metrics.roc_auc_score(trueSplicedLabel, estSplicedScore))
    print('!*!*!*!')

    print('------')

    cor = []
    rsme = []
    estSplicedFlat = []
    trueSplicedFlat = []
    err = []
    tp = 0

    for i in range(len(estPsi)):
        if filter(i):
            #if estSplicedLabel[i] == 1 and trueSplicedLabel[i] == 1:
            if trueSplicedLabel[i] == 1:
                tp += 1
                estSplicedFlat.extend(estPsi[i])
                trueSplicedFlat.extend(truePsi[i])
                cor.append(stats.pearsonr(estPsi[i], truePsi[i])[0])
                rsme.append(np.sqrt(metrics.mean_squared_error(estPsi[i], truePsi[i])))
            elif estSplicedLabel[i] != trueSplicedLabel[i]:
                err.append(i)

    print('--- # genes in ground truth splicing set = \t' + str(tp))
    #print('--- # genes in true positive set = \t' + str(tp))
    print('--- Global RMSE = \t\t', end = '')
    print(np.sqrt(metrics.mean_squared_error(estSplicedFlat, trueSplicedFlat)))
    print('--- Local RMSE = \t\t', end = '')
    print(np.mean(rsme))
    print('!*!*!*!\n--- Global correlation = \t', end = '')
    print(stats.pearsonr(estSplicedFlat, trueSplicedFlat)[0])
    print('!*!*!*!')
    print('--- Local correlation = \t', end = '')
    print(np.mean(cor))

eps = 1e-6

nargv = 1
truePsiPath = sys.argv[nargv];nargv += 1
trueTpmPath = sys.argv[nargv];nargv += 1
estPsiPath = sys.argv[nargv];nargv += 1
estTpmPath = sys.argv[nargv];nargv += 1
outputPath = sys.argv[nargv];nargv += 1


truePsiFile = open(truePsiPath, 'r')
trueTpmFile = open(trueTpmPath, 'r')
estPsiFile = open(estPsiPath, 'r')
estTpmFile = open(estTpmPath, 'r')
outputFile = open(outputPath, 'w')

truePsi = json.load(truePsiFile)
trueTpm = json.load(trueTpmFile)
estPsi = json.load(estPsiFile)
estTpm = json.load(estTpmFile)
sys.stdout = outputFile

print("********** " + outputPath)
statFilter()
print("======")
globalRSME()
#localRSMEMean()
print("======")
globalCorrelation()
#localCorrelationMean()
print("======")
#splicedCorrelation()
#print("======")

print("\n")
