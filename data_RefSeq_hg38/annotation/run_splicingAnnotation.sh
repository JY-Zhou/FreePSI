#!/bin/bash

set -x
set -e

rm -rf splicing_event/*

awk '$3 !~ /_/' \
    transcriptome_raw/hg38_refGene.hgtable > \
    splicing_event/refGene.txt

TIME=splicing_event/time_usage.csv

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    python ~/utility/bioinfo/rnaseqlib-clip/rnaseqlib/gff/gff_make_annotation.py \
    splicing_event/ \
    splicing_event/gff \
    --flanking-rule commonshortest \
    --genome-label hg38
