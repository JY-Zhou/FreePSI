import sys
import random

nargv = 1
gtfPath = sys.argv[nargv]; nargv += 1
sampledGtfPath = sys.argv[nargv]; nargv += 1
percent = sys.argv[nargv]; nargv += 1

gtfFile = open(gtfPath, 'r')
sampledGtfFile = open(sampledGtfPath, 'w')

print('=== ' + sys.argv[0] + ' start...')

percent = float(percent)/ 100.0

gtf = {}
for line in gtfFile:
    refname = line.split('\t')[-1].split(';')[-2]
    if refname not in gtf:
        gtf[refname] = [line.strip()]
    else:
        gtf[refname].append(line.strip())

sampled = 0
for iso in gtf:
    if random.random() < percent:
        sampled += 1
        for x in gtf[iso]:
            print(x, file=sampledGtfFile)

print('Sample rate is %f, retained %d isoforms from %d (%f)' % (percent, sampled, len(gtf.keys()), float(sampled)/len(gtf.keys())))

print('=== ' + sys.argv[0] + ' finish!')
