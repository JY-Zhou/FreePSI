#!/bin/bash

set -e

GENOME_DIR=./genome
BND_FILE=./annotation/hg38_refGene_exonBoundary_chr21.bed
READS=RNA-seq
K=27
THREAD=4

set -x
# Count k-mers in RNA-seq reads using jellyfish
jellyfish count -m ${K} -s 100M -t ${THREAD} -Q 5 ${READS}/reads_final.1.fastq -o ${READS}/reads.1.jf
jellyfish dump ${READS}/reads.1.jf -o ${READS}/reads.1.fa
jellyfish count -m ${K} -s 100M -t ${THREAD} -Q 5 ${READS}/reads_final.2.fastq -o ${READS}/reads.2.jf
jellyfish dump ${READS}/reads.2.jf -o ${READS}/reads.2.fa

# Produce raw estimates of PSI values using FreePSI
#Build
freePSI build\
    -k $K -p ${THREAD} \
    -g ${GENOME_DIR} \
    -a ${BND_FILE} \
    -1 ${READS}/reads.1.fa \
    -2 ${READS}/reads.2.fa \
    -o ./hashtable.json

#Quant
freePSI quant\
    -k $K -p ${THREAD} \
    -i ./hashtable.json \
    -o .

# Post-process the raw estimates of PSI values
freePSI-postProc.py \
    ./psi_freePSI_raw.json \
    ./psi_freePSI.json

# Summarize the PSI values into a readable file
freePSI-summary.py \
    ./annotation/hg38_refGene_exonBoundary_chr21.bed \
    ./psi_freePSI.json \
    ./psi_freePSI.summary
