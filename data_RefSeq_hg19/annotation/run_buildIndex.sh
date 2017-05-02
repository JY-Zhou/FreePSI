#!/bin/bash

set -x
set -e

rm -rf hisat_index/*

TIME=hisat_index/time_usage.csv

REF_FILES=$(ls -1 ../genome/*.fa | tr '\n' ',')
REF_FILES=${REF_FILES%,}

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    hisat2-build ${REF_FILES} hisat_index/index
/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    hisat2-inspect hisat_index/index > hisat_index/index.fa
