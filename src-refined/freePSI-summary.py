import json
import sys

bedFile = open(sys.argv[1], 'r')
psiFile = open(sys.argv[2], 'r')
summaryFile = open(sys.argv[3], 'w')

header = "Chromosome\tStrand\tExon starting site\tExon ending site\tPSI"
print(header, file=summaryFile)
psi = json.load(psiFile)
k = 0
for line in bedFile:
    substr = line.strip().split('\t')
    if int(substr[-3]) != len(psi[k]):
        print('Error: the output doesn\'t match the input!')
        exit()
    exSt = substr[-1].split(',')
    exLen = substr[-2].split(',')
    for i in range(int(substr[-3])):
        summary = substr[0] + '\t'
        summary += substr[5] + '\t'
        summary += str(int(substr[1]) + int(exSt[i])) + '\t'
        summary += str(int(substr[1]) + int(exSt[i]) + int(exLen[i])) + '\t'
        summary += str(psi[k][i])
        print(summary, file=summaryFile)
    k += 1
