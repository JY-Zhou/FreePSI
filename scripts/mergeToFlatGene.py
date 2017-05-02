import sys
import os
import json

print('\n=== %s start... ===' % sys.argv[0])

argId = 1
isoformPath = sys.argv[argId]
argId += 1
boundaryPath = sys.argv[argId]
argId += 1
isoformConfigPath = sys.argv[argId]
argId += 1
nameMapPath = sys.argv[argId]

isoformFile = open(isoformPath, 'r')
boundaryFile = open(boundaryPath, 'w')
isoformConfigFile = open(isoformConfigPath, 'w') 
nameMapFile = open(nameMapPath, 'w')

isoform = []
clusterId = []
chrSt = {}
chrEd = {}
k = 0
for line in isoformFile:
    substr = line[:-1].split('\t')
    chrName = substr[0]
    isoform.append(tuple(substr[:-1]))
    clusterId.append(int(substr[-1]))
    if not chrName in chrSt:
        chrSt[chrName] = k
    chrEd[chrName] = k
    k += 1

curCluster = 0
rawGene = []
for i in range(len(isoform)):
    if curCluster != clusterId[i]:
        curCluster = clusterId[i]
        rawGene.append([])
    rawGene[-1].append(i)

gene = []
geneBounds = []
minOverlap = 0.3
for g in range(len(rawGene)):
    j = 0
    gl = gr = 0
    gsymbol = []
    gene.append([])
    geneBounds.append((0, 0))

    while j < len(rawGene[g]):
        i = rawGene[g][j]
        l = int(isoform[i][1])
        r = int(isoform[i][2])
        symbol = isoform[i][4]
        if len(gene[-1]) == 0:
            gl = l
            gr = r
            gsymbol.append(symbol)
            gene[-1].append(i)
            geneBounds[-1] = (gl, gr)
            j += 1
        else:
            if symbol in gsymbol:
                gl = min(gl, l)
                gr = max(gr, r)
                gene[-1].append(i)
                geneBounds[-1] = (gl, gr)
                j += 1
            else:
                inter = min(gr, r) - max(gl, l)
                #union = max(gr, r) - min(gl, l)
                rate = float(inter) / min(gr-gl, r-l)
                #print(symbol + ' -- rate =' + str(rate))
                if rate > minOverlap:
                    gl = min(gl, l)
                    gr = max(gr, r)
                    gsymbol.append(symbol)
                    gene[-1].append(i)
                    geneBounds[-1] = (gl, gr)
                    j += 1
                else:
                    gl = gr = 0
                    gsymbol = []
                    gene.append([])
                    geneBounds.append((0, 0))

print('# of cluster genes = ' + str(len(rawGene)))
print('# of refined cluster genes = ' + str(len(gene)))
#for g in range(len(gene)):
#    for i in gene[g]:
#        print(isoform[i][3], end=',' )
#    print('')

minLen = 30
NG = 0
NE = []
L = []
ISO = []
ISOProf = {}
allNewExon = []
nameMap = {}

for g in range(len(gene)):
    corIsoform = []
    for i in gene[g]:
        corIsoform.append(isoform[i])

    (st, ed) = geneBounds[g]
    chrName = corIsoform[0][0]
    strand = corIsoform[0][5]

    exons = []
    for iso in corIsoform:
        isost = int(iso[1])
        exnum = int(iso[-3])
        exst_str = iso[-2].split(',')
        exed_str = iso[-1].split(',')
        exst = []
        exlen = []
        for i in range(exnum):
            exst.append(str(int(exst_str[i]) - isost))
            exlen.append(str(int(exed_str[i]) - int(exst_str[i])))
        for i in range(exnum):
            exons.append(( isost + int(exst[i]) , isost + int(exst[i]) + int(exlen[i]) ))

    exonbndset = {}
    for e in exons:
        exonbndset[e[0]] = 1
        exonbndset[e[1]] = 1

    exonbnd = exonbndset.keys()
    exonbnd = sorted(exonbnd)

    newExon = []
    for i in range(1, len(exonbnd)):
        for e in exons:
            if exonbnd[i-1] >= e[0] and exonbnd[i] <= e[1]:
                newExon.append((exonbnd[i-1], exonbnd[i]))
                break

    checker = True
    while checker:
        checker = False
        for i in range(len(newExon)):
            if newExon[i][1] - newExon[i][0] <= minLen:
                checker = True
                if i > 0 and newExon[i-1][1] == newExon[i][0]:
                    newExon[i-1] = (newExon[i-1][0], newExon[i][1])
                elif i < len(newExon) - 1 and newExon[i+1][0] == newExon[i][1]:
                    newExon[i+1] = (newExon[i][0], newExon[i+1][1])
                del newExon[i]
                break

    allNewExon.append(newExon)

    if len(newExon) > 0:
        newTxSt = newExon[0][0]
        newTxEd = newExon[0][1]
        for x in newExon:
            newTxSt = min(newTxSt, x[0])
            newTxEd = max(newTxEd, x[1])

        row = str(chrName) + '\t' + str(newTxSt) + '\t' + str(newTxEd) + '\t'
        row = row + 'Gene' + str(NG) + '\t0\t' + str(strand) + '\t'
        row = row + str(newTxSt) + '\t' + str(newTxEd) + '\t0\t' + str(len(newExon)) + '\t'
        for x in newExon:
            row = row + str(x[1] - x[0]) + ','
        row = row + '\t'
        for x in newExon:
            row = row + str(x[0] - newTxSt) + ','
        print(row, file=boundaryFile)
        

        ISO.append([])
        ISOProf['Gene' + str(NG)] = {}
        for i in range(len(corIsoform)):
            isost = int(corIsoform[i][1])
            exnum = int(corIsoform[i][-3])
            exst_str = corIsoform[i][-2].split(',')
            exed_str = corIsoform[i][-1].split(',')
            exst = []
            exlen = []
            for j in range(exnum):
                exst.append(str(int(exst_str[j]) - isost))
                exlen.append(str(int(exed_str[j]) - int(exst_str[j])))

            nameMap['Gene' + str(NG) + 'Iso' + str(i)] = corIsoform[i][3]
            nameMap[corIsoform[i][3]] = 'Gene' + str(NG) + 'Iso' + str(i)

            ISO[-1].append([])
            ISOProf['Gene' + str(NG)]['Gene' + str(NG) + 'Iso' + str(i)] = []

            for k in range(len(newExon)):
                el = newExon[k][0]
                er = newExon[k][1]
                add = False
                for j in range(exnum):
                    if int(exlen[j]) > minLen:
                        l = isost + int(exst[j])
                        r = isost + int(exst[j]) + int(exlen[j])
                        overlap = float(min(er, r) - max(el, l))
                        rate = max(overlap / (er-el), overlap / (r - l))
                        if rate > 0.8:
                            add = True
                if add:
                    ISO[-1][-1].append(k)
                    ISOProf['Gene' + str(NG)]['Gene' + str(NG) + 'Iso' + str(i)].append(k)

            if len(ISO[-1][-1]) == 0:
                print(chrName + ' EMPTY ISO..')
                print(corIsoform[i])
                print(rate)
                for x in corIsoform:
                    print(x)
                print('Gene %d, #exon = %d, %s' % (g, len(newExon), corIsoform[0][4] ))
                continue

        NE.append(len(newExon))
        NG += 1
        L.append([])
        for x in newExon:
            L[-1].append(x[1] - x[0])

    else:
        print('Gene includes no exons ...')
        print('Gene %d, #exon = %d, %s' % (g, len(newExon), corIsoform[0][3] ))
        pass

config = {
        'NG' : NG,
        'NE' : NE,
        'Iso': ISOProf,
        'L' : L
        }

tot = 0
for g in range(NG):
    coverage = [0 for i in range(NE[g])]
    for i in range(len(ISO[g])):
        for e in ISO[g][i]:
            coverage[e] += 1
    for e in range(len(coverage)):
        if coverage[e] == 0:
            tot += 1
            print(g)
            print(e)
print('# of uncovered exons = ' + str(tot))

json.dump(config, isoformConfigFile, indent=4)
json.dump(nameMap, nameMapFile, indent=4)

print('=== %s finish! ===\n' % sys.argv[0])
