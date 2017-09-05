import random
import json
import pandas 
import scipy.stats as stats
import sys

print('\n=== %s start... ===' % sys.argv[0])

nargv = 1
gtfPath = sys.argv[nargv];nargv += 1
salmonPath = sys.argv[nargv];nargv += 1
bndPath = sys.argv[nargv];nargv += 1
profPath = sys.argv[nargv];nargv += 1
namePath = sys.argv[nargv];nargv += 1
psiPath = sys.argv[nargv];nargv += 1
tpmPath = sys.argv[nargv];nargv += 1

gtfFile = open(gtfPath, 'r')
salmonFile = open(salmonPath, 'r')
bndFile = open(bndPath, 'r')
profFile = open(profPath, 'r')
nameFile = open(namePath, 'r')
psiFile = open(psiPath, 'w')
tpmFile = open(tpmPath, 'w')

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

transidTPM = {}
for line in salmonFile:
    substr = line.strip().split('\t')
    if not 'Name' in substr[0]:
        transidTPM[substr[0]] = float(substr[3])

estTPM = {}
#estFPKMTot = 0.0
for line in gtfFile:
    substr = line.strip().split('\t')
    if substr[2] == 'transcript':
        subInfo = substr[8].split(';')
        for info in subInfo:
            if 'transcript_id' in info.split('"')[0]:
                isofId = info.split('"')[1]
            #if 'FPKM' in info.split('"')[0]:
            #    fpkm = info.split('"')[1]
        #estFPKMTot += float(fpkm)
        estTPM[isofId] = transidTPM[isofId]
#for isofId in estTPM:
#    estTPM[isofId] /= estFPKMTot

estPsi = [[0.0 for e in range(NE[g])] for g in range(NG)]
estTot = [0.0 for g in range(NG)]
refGeneId = -1
refGeneInfo = None
spamGene = {}
lossExon = []
gtfFile = open(gtfPath, 'r')
for line in gtfFile:
    substr = line.strip().split('\t')
    chrName = substr[0]
    region = substr[2]
    st = int(substr[3]) - 1
    ed = int(substr[4])
    strand = substr[6]
    subInfo = substr[8].split(';')
    for info in subInfo:
        if 'transcript_id' in info.split('"')[0]:
            isofId = info.split('"')[1]

    #print((chrName, region, st, ed, strand, isofId, estTPM[isofId]))
    if region == 'transcript':
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
            #print('GeneID = %s, tpm = %f' % (refGeneId, estTPM[isofId]))
            #print(refGeneInfo)
            #print('Match rate:' + str(maxRate))
            if refGeneId == -1:
                spamGene[isofId] = 'Unmapped gene'
            else:
                estTot[refGeneId] += estTPM[isofId]
        else:
            spamGene[isofId] = 'Lack strand info'
    elif region == 'exon':
        if not isofId in spamGene:
            mapped = False
            for e in range(refGeneInfo[2]):
                exonSt = refGeneInfo[3][e][0]
                exonEd = refGeneInfo[3][e][1]
                overlap = float(min(ed, exonEd) - max(st, exonSt))
                rate = overlap / (exonEd - exonSt)
                if rate > 0.8:
                    mapped = True
                    #print('Exon Id = %d, overlap = %f, len = %f, rate = %f' % (e, overlap, exonEd - exonSt, rate))
                    estPsi[refGeneId][e] += estTPM[isofId]
            if not mapped:
                for info in subInfo:
                    if 'exon_number' in info.split('"')[0]:
                        exonId = info.split('"')[1]
                lossExon.append((isofId, exonId))

print(len(spamGene))
empty = 0
for g in range(NG):
    if estTot[g] > 0.0:
        empty += 1
    for e in range(NE[g]):
        if estTot[g] > 0.0:
            estPsi[g][e] /= estTot[g]
        else:
            estPsi[g][e] = 0.0

print(NG)
print(empty)

json.dump(estPsi, psiFile, indent = 4)
json.dump(estTot, tpmFile, indent = 4)

print('=== %s finished! ===\n' % sys.argv[0])
