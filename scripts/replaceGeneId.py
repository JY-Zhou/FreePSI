import re
import sys

print('\n=== %s start... ===' % sys.argv[0])

argId = 1
mapPath = sys.argv[argId]
argId += 1
gtfInPath = sys.argv[argId]
argId += 1
gtfOutPath = sys.argv[argId]

mapFile = open(mapPath, 'r')
gtfInFile = open(gtfInPath, 'r')
gtfOutFile = open(gtfOutPath, 'w')

nameMap = {}
for line in mapFile:
    if line[0] != '#':
        substr = line.strip().split('\t')
        nameMap[substr[1]] = substr[12]

txMap = {}
k = 0
for line in gtfInFile:
    k += 1
    substr = line.strip().split('\t')[-1].replace('"','').replace(';','').split(' ')
    if substr[1] != substr[3]:
        continue
    
    chrId = line.split('\t')[0]
    if substr[1] in txMap:
        if chrId != txMap[substr[1]]:
            continue
    else:
        txMap[substr[1]] = chrId

    line = re.sub(r'gene_id "' + substr[1] + r'";', 
            r'gene_id "' + nameMap[substr[1]] + r'";',
            line)
    print(line.strip(), file = gtfOutFile)

print('\n=== %s finished!  ===' % sys.argv[0])
