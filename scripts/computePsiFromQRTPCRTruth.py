import random
import json
import pandas 
import scipy.stats as stats
import sys

def add_exon(chrName, strand, st, ed, psi):
    if abs(ed - st) < 30:
        return False
    if chrName + strand in refGeneBnd:
        refGeneId = -1
        refGeneInfo = None
        maxRate = 0
        for geneInfo in refGeneBnd[chrName + strand]:
            geneSt = geneInfo[1][0]
            geneEd = geneInfo[1][1]
            overlap = float(min(ed, geneEd) - max(st, geneSt))
            rate = max(overlap / (geneEd - geneSt), overlap / (ed - st))
            if rate > maxRate:
                maxRate = rate
                refGeneId = int(geneInfo[0].replace('Gene', ''))
                refGeneInfo = geneInfo
        if refGeneId == -1:
            return False

        mapped = False
        for e in range(refGeneInfo[2]):
            exonSt = refGeneInfo[3][e][0]
            exonEd = refGeneInfo[3][e][1]
            overlap = float(min(ed, exonEd) - max(st, exonSt))
            rate = overlap / (exonEd - exonSt)
            if rate > 0.8:
                mapped = True
                estPsi[refGeneId][e] += psi
                estMask[refGeneId][e] += 1
        if not mapped:
            return False
    else:
        print('Unknown strand info')
        return False
    return True

print('\n=== %s start... ===' % sys.argv[0])

nargv = 1
qRTPCRPath = sys.argv[nargv];nargv += 1
bndPath = sys.argv[nargv];nargv += 1
profPath = sys.argv[nargv];nargv += 1
namePath = sys.argv[nargv];nargv += 1
psiPath = sys.argv[nargv];nargv += 1
maskPath = sys.argv[nargv];nargv += 1

qRTPCRFile = open(qRTPCRPath, 'r')
bndFile = open(bndPath, 'r')
profFile = open(profPath, 'r')
nameFile = open(namePath, 'r')
psiFile = open(psiPath, 'w')
maskFile = open(maskPath, 'w')

prof = json.load(profFile)
NG = prof['NG']
NE = prof['NE']
iso = prof['Iso']

nameMap = json.load(nameFile)

refGeneBnd = {}
for line in bndFile:
    substr = line.strip().split('\t')
    chrName = substr[0]
    geneBnd = (int(substr[1]), int(substr[2]))
    geneName = substr[3]
    geneStrand = substr[5]
    exonNum = int(substr[9])
    exonLenStr = substr[10].split(',')
    exonStStr = substr[11].split(',')
    exonBnd = []
    for i in range(exonNum):
        l = geneBnd[0] + int(exonStStr[i])
        r = l + int(exonLenStr[i])
        exonBnd.append((l, r))

    geneInfo = (geneName, geneBnd, exonNum, exonBnd)
    if chrName + geneStrand in refGeneBnd:
        refGeneBnd[chrName + geneStrand].append(geneInfo)
    else:
        refGeneBnd[chrName + geneStrand] = [geneInfo]

estPsi = [[0.0 for e in range(NE[g])] for g in range(NG)]
estMask = [[0 for e in range(NE[g])] for g in range(NG)]
refGeneId = -1
refGeneInfo = None
spamGene = {}
lossExon = []
k = 0
failexon = 0
for line in qRTPCRFile:
    substr = line.strip().split('\t')
    k += 1
    chrName = substr[2]
    psi = float(substr[6])
    st = int(substr[4])
    ed = int(substr[5])
    strand = substr[3]
    if not add_exon(chrName, strand, min(st, ed), max(st, ed), psi):
        failexon += 1

print('Total %d, %d failed' % (k, failexon))

totExon = 0
for g in range(NG):
    for e in range(NE[g]):
        if estMask[g][e] > 0:
            totExon += 1
            estPsi[g][e] /= estMask[g][e]
    
print('Total %d spliced exons' % (totExon))
json.dump(estPsi, psiFile, indent = 4)
json.dump(estMask, maskFile, indent = 4)

print('=== %s finished! ===\n' % sys.argv[0])
