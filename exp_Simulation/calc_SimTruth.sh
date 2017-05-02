#!/bin/bash

set -x
set -e

if [ $# != 1 ]; then
    echo "Require reads directory!"
    exit 1;
fi

READS=${1}
RESULT=result/${1}_SimTruth

rm -rf ${RESULT}
mkdir ${RESULT}

python3 ../scripts/computePsiFromSimTruth_TPMPSI.py \
    ${READS}/flux_simulator_main_expression.pro \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_exonProfile.json \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_nameMap.json \
    ${RESULT}/psi_SimTruth.json \
    ${RESULT}/tpm_SimTruth.json 

python3 ../scripts/computePsiFromSimTruth_MASK.py \
    ../data_RefSeq_hg38/annotation/splicing_event/gff/commonshortest/ \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_exonBoundary.bed \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_exonProfile.json \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_nameMap.json \
    ${RESULT}/mask_SimTruth.json
