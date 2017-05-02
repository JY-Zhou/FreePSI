#!/bin/bash

set -x
set -e

if [ $# != 2 ] ; then 
    echo "Require experiment name and tool name!" 
    exit 1; 
fi 

TRUTH=result/${1}_SimTruth
ESTIMATION=result/${1}_${2}
EVALUATION=evaluation/${1}_spliced

if [ ! -d "${EVALUATION}" ]; then
    mkdir ${EVALUATION}
fi

rm -rf ${EVALUATION}/eval_${2}.out

python3 ../scripts/evaluateSplicedPsi.py \
    ${TRUTH}/psi_SimTruth.json \
    ${TRUTH}/tpm_SimTruth.json \
    ${ESTIMATION}/psi_${2}.json \
    ${TRUTH}/mask_SimTruth.json \
    ${EVALUATION}/eval_${2}.out
