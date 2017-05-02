#!/bin/bash

set -x
set -e

if [ $# != 2 ] ; then 
    echo "Require experiment name and tool name!" 
    exit 1; 
fi 

TRUTH=result/${1}_qRTPCRTruth
ESTIMATION=result/${1}_${2}
EVALUATION=evaluation/${1}

if [ ! -d "${EVALUATION}" ]; then
    mkdir ${EVALUATION}
fi

rm -rf ${EVALUATION}/eval_${2}.out

python3 ../scripts/evaluateqRTPCRPsi.py \
    ${TRUTH}/psi_qRTPCRTruth.json \
    ${ESTIMATION}/psi_${2}.json \
    ${TRUTH}/mask_qRTPCRTruth.json \
    ${EVALUATION}/eval_${2}.out
