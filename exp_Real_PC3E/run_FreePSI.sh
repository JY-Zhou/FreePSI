#!/bin/bash

set -x
set -e

if [ $# != 1 ] ; then 
    echo "Require reads directory!" 
    exit 1; 
fi 

READS=PC3E/RNAseq
OUTPUT=output/${1}_FreePSI
RESULT=result/${1}_FreePSI
TIME=${RESULT}/time_usage.csv

BND_FILE=../data_RefSeq_hg19/annotation/exon_boundary/hg19_refGene_exonBoundary.bed
GENOME_DIR=../data_RefSeq_hg19/genome

K=25

rm -rf ${OUTPUT} ${RESULT}
mkdir ${OUTPUT} ${RESULT}


/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    jellyfish count -m ${K} -s 2G -t 16 ${READS}/${1}_1.fastq -o ${OUTPUT}/reads.1.jf -Q A -L 10

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    jellyfish dump ${OUTPUT}/reads.1.jf -o ${OUTPUT}/reads.1.fa

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    jellyfish count -m ${K} -s 2G -t 16 ${READS}/${1}_2.fastq -o ${OUTPUT}/reads.2.jf -Q A -L 10

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    jellyfish dump ${OUTPUT}/reads.2.jf -o ${OUTPUT}/reads.2.fa

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    ../tools/freePSI -k ${K} -l 30 \
    -g ${GENOME_DIR} \
    -a ${BND_FILE} \
    -1 ${OUTPUT}/reads.1.fa \
    -2 ${OUTPUT}/reads.2.fa \
    -o ${OUTPUT}

cp ${OUTPUT}/tpm_freePSI_raw.json ${RESULT}/tpm_FreePSI.json

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    python3 ../scripts/computePsiFromFreePSI.py \
    ${OUTPUT}/psi_freePSI_raw.json \
    ${RESULT}/psi_FreePSI.json
