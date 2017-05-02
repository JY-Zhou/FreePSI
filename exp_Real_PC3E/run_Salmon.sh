#!/bin/bash

set -x
set -e

if [ $# != 1 ] ; then 
    echo "Require reads directory!" 
    exit 1; 
fi 

READS=PC3E/RNAseq
OUTPUT=output/${1}_Salmon
RESULT=result/${1}_Salmon
TIME=${RESULT}/time_usage.csv

GTF_FILE=../data_RefSeq_hg19/annotation/transcriptome_cleaned/hg19_refGene_cleaned.gtf
GEN_DIR=../data_RefSeq_hg19/genome

rm -rf ${OUTPUT} ${RESULT}
mkdir ${OUTPUT} ${RESULT}

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    rsem-prepare-reference \
    --gtf ${GTF_FILE} \
    ${GEN_DIR} \
    ${OUTPUT}/salmon

# Now create the Salmon transcript index.
/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    salmon index \
    -t ${OUTPUT}/salmon.transcripts.fa \
    -i ${OUTPUT}/index

# Use Salmon to calculate per-transcript TPMs.
mkdir ${OUTPUT}/quant
/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    salmon quant \
    -p 16 \
    -l ISF \
    -i ${OUTPUT}/index \
    -1 ${READS}/${1}_1.fastq \
    -2 ${READS}/${1}_2.fastq \
    -o ${OUTPUT}/quant

# Compute PSI from output
/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    python3 ../scripts/computePsiFromReferred.py \
    ${OUTPUT}/quant/quant.sf \
    ../data_RefSeq_hg19/annotation/exon_boundary/hg19_refGene_exonProfile.json \
    ../data_RefSeq_hg19/annotation/exon_boundary/hg19_refGene_nameMap.json \
    ${RESULT}/psi_Salmon.json \
    ${RESULT}/tpm_Salmon.json
