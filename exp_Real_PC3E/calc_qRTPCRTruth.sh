#!/bin/bash

set -x
set -e

if [ $# != 1 ]; then
    echo "Require reads directory!"
    exit 1;
fi

QRTPCR=PC3E/qRTPCR/PC3E.txt
RESULT=result/${1}_qRTPCRTruth

rm -rf ${RESULT}
mkdir ${RESULT}

python3 ../scripts/computePsiFromQRTPCRTruth.py \
    ${QRTPCR} \
    ../data_RefSeq_hg19/annotation/exon_boundary/hg19_refGene_exonBoundary.bed \
    ../data_RefSeq_hg19/annotation/exon_boundary/hg19_refGene_exonProfile.json \
    ../data_RefSeq_hg19/annotation/exon_boundary/hg19_refGene_nameMap.json \
    ${RESULT}/psi_qRTPCRTruth.json \
    ${RESULT}/mask_qRTPCRTruth.json 
