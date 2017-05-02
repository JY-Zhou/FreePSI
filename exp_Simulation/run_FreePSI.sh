#!/bin/bash

set -x
set -e

if [ $# != 4 ] ; then 
    echo "Require reads directory and thread number" 
    exit 1; 
fi 

K=${3}
R=${4}
THREAD=${2}
READS=${1}
OUTPUT=output/${1}_FreePSI_${K}_${R}
RESULT=result/${1}_FreePSI_${K}_${R}
TIME=${RESULT}/time_usage.csv

BND_FILE=../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_exonBoundary.bed
GENOME_DIR=../data_RefSeq_hg38/genome


rm -rf ${OUTPUT} ${RESULT}
mkdir ${OUTPUT} ${RESULT}

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    jellyfish count \
    -m $K -s 2G -t ${THREAD} -Q 5\
    ${READS}/reads_final.1.fastq \
    -o ${OUTPUT}/reads.1.jf \

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    jellyfish dump ${OUTPUT}/reads.1.jf -o ${OUTPUT}/reads.1.fa

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    jellyfish count \
    -m $K -s 2G -t ${THREAD} -Q 5\
    ${READS}/reads_final.2.fastq \
    -o ${OUTPUT}/reads.2.jf \

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    jellyfish dump ${OUTPUT}/reads.2.jf -o ${OUTPUT}/reads.2.fa

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    ../tools/freePSI -k $K -l $R -p ${THREAD}\
    -g ${GENOME_DIR} \
    -1 ${OUTPUT}/reads.1.fa \
    -2 ${OUTPUT}/reads.2.fa \
    -a ${BND_FILE} \
    -o ${OUTPUT}

cp ${OUTPUT}/tpm_freePSI_raw.json ${RESULT}/tpm_FreePSI_${K}_${R}.json

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    python3 ../scripts/computePsiFromFreePSI.py \
    ${OUTPUT}/psi_freePSI_raw.json \
    ${RESULT}/psi_FreePSI_${K}_${R}.json
