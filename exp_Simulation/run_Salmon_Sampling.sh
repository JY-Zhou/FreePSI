#!/bin/bash

set -x
set -e

if [ $# != 2 ] ; then 
    echo "Require reads directory!" 
    exit 1; 
fi 

READS=${1}
OUTPUT=output/${1}_Salmon_Sampling_${2}
RESULT=result/${1}_Salmon_Sampling_${2}
TIME=${RESULT}/time_usage.csv

GTF_FILE=../data_RefSeq_hg38/annotation/transcriptome_cleaned/hg38_refGene_cleaned.gtf
GEN_DIR=../data_RefSeq_hg38/genome

rm -rf ${OUTPUT} ${RESULT}
mkdir ${OUTPUT} ${RESULT}

python3 ../scripts/sampleGtf.py \
    ${GTF_FILE} \
    ${OUTPUT}/hg38_refGene_subset.gtf \
    ${2}

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    rsem-prepare-reference \
    --gtf ${OUTPUT}/hg38_refGene_subset.gtf \
    ${GEN_DIR} \
    ${OUTPUT}/salmon

# Now create the Salmon_Sampling transcript index.
/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    salmon index \
    -t ${OUTPUT}/salmon.transcripts.fa \
    -i ${OUTPUT}/index

# Use Salmon_Sampling to calculate per-transcript TPMs.
mkdir ${OUTPUT}/quant
/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    salmon quant \
    -p 16 \
    -l ISF \
    -i ${OUTPUT}/index \
    -1 ${READS}/reads_final.1.fastq \
    -2 ${READS}/reads_final.2.fastq \
    -o ${OUTPUT}/quant

# Compute PSI from output
/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    python3 ../scripts/computePsiFromReferred.py \
    ${OUTPUT}/quant/quant.sf \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_exonProfile.json \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_nameMap.json \
    ${RESULT}/psi_Salmon_Sampling_${2}.json \
    ${RESULT}/tpm_Salmon_Sampling_${2}.json
