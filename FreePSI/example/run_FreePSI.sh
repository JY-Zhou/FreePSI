#!/bin/bash

set -e

# Provide the directory containing Jellyfish (usually named as 'bin')  
Jellyfish=
# Provide the directory containing FreePSI
# E.g. 
#FreePSI=../bin/
FreePSI=

if [ -z ${Jellyfish} ]; then
    echo "Error: Please modify me to provide the directory containing Jellyfish"
    exit;
fi

if [ -z ${FreePSI} ]; then
    echo "Error: Please modify me to provide the directory containing FreePSI"
    exit;
fi

GENOME_DIR=./genome
BND_FILE=./annotation/hg38_refGene_exonBoundary_chr21.bed
READS=RNA-seq
K=27
THREAD=4

set -x
# Count k-mers in RNA-seq reads using jellyfish
${Jellyfish}/jellyfish count -m ${K} -s 2G -t ${THREAD} -Q 5 ${READS}/reads_final.1.fastq -o ${READS}/reads.1.jf
${Jellyfish}/jellyfish dump ${READS}/reads.1.jf -o ${READS}/reads.1.fa
${Jellyfish}/jellyfish count -m ${K} -s 2G -t ${THREAD} -Q 5 ${READS}/reads_final.2.fastq -o ${READS}/reads.2.jf
${Jellyfish}/jellyfish dump ${READS}/reads.2.jf -o ${READS}/reads.2.fa

# Produce raw estimates of PSI values using FreePSI
${FreePSI}/freePSI -k $K -p ${THREAD} \
    -g ${GENOME_DIR} \
    -a ${BND_FILE} \
    -1 ${READS}/reads.1.fa \
    -2 ${READS}/reads.2.fa \
    -o .

# Post-process the raw estimates of PSI values
python3 ${FreePSI}/postProc.py \
    ./psi_freePSI_raw.json \
    ./psi_freePSI.json

# Summarize the PSI values into a readable file
python3 ${FreePSI}/summary.py \
    ./annotation/hg38_refGene_exonBoundary_chr21.bed \
    ./psi_freePSI.json \
    ./psi_freePSI.summary 
