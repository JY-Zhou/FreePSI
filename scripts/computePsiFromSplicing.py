import random
import json
import pandas 
import scipy.stats as stats
import sys

def add_exon(chrName, strand, st, ed, psi):
    if abs(ed - st) < 30:
        return
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
            return

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
            pass
    else:
        print('Unknown strand info')

print('\n=== %s start... ===' % sys.argv[0])

nargv = 1
misoPath = sys.argv[nargv];nargv += 1
bndPath = sys.argv[nargv];nargv += 1
profPath = sys.argv[nargv];nargv += 1
namePath = sys.argv[nargv];nargv += 1
psiPath = sys.argv[nargv];nargv += 1

misoFile = open(misoPath, 'r')
bndFile = open(bndPath, 'r')
profFile = open(profPath, 'r')
nameFile = open(namePath, 'r')
psiFile = open(psiPath, 'w')

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
for line in misoFile:
    substr = line.strip().split('\t')
    if len(substr) < 11:
        continue
    if not substr[1].replace('.', '', 1).isdigit():
        continue
    psi = float(substr[1])
    info = substr[0].split('@')
    if len(info) == 2:
        if '|' in info[0]:
            #continue
            #print('A5SS')
            var = info[0].split(':')
            chrName = var[0]
            idx = var[2].split('|')
            st = int(idx[0])
            ed = int(idx[1])
            strand = var[3]
            add_exon(chrName, strand, min(st, ed), max(st, ed), psi)
        elif '|' in info[1]:
            #continue
            #print('A3SS')
            #print(info)
            var = info[1].split(':')
            chrName = var[0]
            idx = var[1].split('|')
            st = int(idx[0])
            ed = int(idx[1])
            strand = var[3]
            #print(st, ed)
            add_exon(chrName, strand, min(st, ed), max(st, ed), psi)
            #input()
        elif '-' in substr[0]:
            #continue
            #print('RI')
            var = info[0].split(':')
            chrName = var[0]
            strand = var[2]
            var = info[0].split(':')[1].split('-')
            idx = []
            idx.append(int(var[0]))
            idx.append(int(var[1]))
            var = info[1].split(':')[1].split('-')
            idx.append(int(var[0]))
            idx.append(int(var[1]))
            idx = sorted(idx)
            add_exon(chrName, strand, idx[1], idx[2], psi)
    elif len(info) == 3:
        #continue
        #print('SE')
        var = info[1].split(':')
        chrName = var[0]
        st = int(var[1])
        ed = int(var[2])
        strand = var[3]
        add_exon(chrName, strand, min(st, ed), max(st, ed), psi)
    elif len(info) == 4:
        #continue
        #print('MXE')
        var = info[1].split(':')
        chrName = var[0]
        strand = var[3]
        st = int(var[1])
        ed = int(var[2])
        add_exon(chrName, strand, min(st, ed), max(st, ed), psi)
        var = info[2].split(':')
        st = int(var[1])
        ed = int(var[2])
        add_exon(chrName, strand, min(st, ed), max(st, ed), 1 - psi)
    else:
        print("Error event name!")

for g in range(NG):
    for e in range(NE[g]):
        if estMask[g][e] > 0:
            estPsi[g][e] /= estMask[g][e]
    
json.dump(estPsi, psiFile, indent = 4)

print('=== %s finished! ===\n' % sys.argv[0])
