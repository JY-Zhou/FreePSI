#!/bin/bash
set -e
set -x

rm -rf exon_boundary/*

ANNOTATION=hg38_refGene
TIME=exon_boundary/time_usage.csv

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    awk '!($4 in a) {a[$4] = 1;print}' \
    transcriptome_cleaned/${ANNOTATION}_cleaned.shortbed > \
    exon_boundary/${ANNOTATION}_unique.shortbed

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    sort -k1,1 -k2,2n \
    exon_boundary/${ANNOTATION}_unique.shortbed > \
    exon_boundary/${ANNOTATION}_sorted.shortbed

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    bedtools cluster -s -i \
    exon_boundary/${ANNOTATION}_sorted.shortbed > \
    exon_boundary/${ANNOTATION}_cluster.shortbed

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    python3 ../../scripts/mergeToFlatGene.py \
    exon_boundary/${ANNOTATION}_cluster.shortbed \
    exon_boundary/${ANNOTATION}_exonBoundary.bed \
    exon_boundary/${ANNOTATION}_exonProfile.json \
    exon_boundary/${ANNOTATION}_nameMap.json
