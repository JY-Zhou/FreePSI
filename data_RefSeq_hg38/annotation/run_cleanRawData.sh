#!/bin/bash
set -e
set -x

ANNOTATION=hg38_refGene

python3 ../../scripts/replaceGeneId.py \
    transcriptome_raw/${ANNOTATION}.hgtable \
    transcriptome_raw/${ANNOTATION}.gtf \
    transcriptome_cleaned/${ANNOTATION}.gtf 

awk '$1 !~ /_/ && $3 ~ /exon/ && $0 !~ /MIR/' \
    transcriptome_cleaned/${ANNOTATION}.gtf > \
    transcriptome_cleaned/${ANNOTATION}_cleaned.gtf

awk '$0 !~ /\#/ && $3 !~ /_/ && $13 !~ /MIR/ {printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$3,$5,$6,$2,$13,$4,$9,$10,$11)}' \
    transcriptome_raw/${ANNOTATION}.hgtable > \
    transcriptome_cleaned/${ANNOTATION}_cleaned.shortbed
