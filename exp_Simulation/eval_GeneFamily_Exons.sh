#!/bin/bash

set -x
set -e

if [ $# != 2 ] ; then 
    echo "Require experiment name and tool name!" 
    exit 1; 
fi 

TRUTH=result/${1}_SimTruth
ESTIMATION=result/${1}_${2}
EVALUATION=evaluation/${1}_geneFamily

if [ ! -d "${EVALUATION}" ]; then
    mkdir ${EVALUATION}
fi

rm -rf ${EVALUATION}/eval_${2}.out

python3 ../scripts/evaluateGeneFamilyPsi.py \
    ../data_RefSeq_hg38/annotation/gene_family/HGNC_Gene_Family.txt \
    multi-reads/multi.out \
    ../data_RefSeq_hg38/annotation/transcriptome_cleaned/hg38_refGene_cleaned.shortbed \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_nameMap.json \
    ${TRUTH}/psi_SimTruth.json \
    ${TRUTH}/tpm_SimTruth.json \
    ${ESTIMATION}/psi_${2}.json \
    ${ESTIMATION}/tpm_${2}.json \
    ${EVALUATION}/eval_${2}.out
