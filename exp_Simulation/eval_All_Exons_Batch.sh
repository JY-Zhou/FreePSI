#!/bin/bash

set -x
set -e

if [ $# != 1 ] ; then 
    echo "Require experiment name!" 
    exit 1; 
fi 

EVALUATION=evaluation/${1}_all

rm -rf ${EVALUATION}
mkdir ${EVALUATION}

./eval_All_Exons.sh ${1} Salmon
./eval_All_Exons.sh ${1} Hisat_Cufflinks
./eval_All_Exons.sh ${1} FreePSI
./eval_All_Exons.sh ${1} Hisat_Cufflinks_Annot
./eval_All_Exons.sh ${1} Hisat_Cufflinks_Annot_Sampling_90
./eval_All_Exons.sh ${1} Hisat_Cufflinks_Annot_Sampling_80
./eval_All_Exons.sh ${1} Hisat_Cufflinks_Annot_Sampling_70
./eval_All_Exons.sh ${1} Salmon_Sampling_90
./eval_All_Exons.sh ${1} Salmon_Sampling_80
./eval_All_Exons.sh ${1} Salmon_Sampling_70
./eval_All_Exons.sh ${1} FreePSI_21
./eval_All_Exons.sh ${1} FreePSI_23
./eval_All_Exons.sh ${1} FreePSI_25
./eval_All_Exons.sh ${1} FreePSI_27
./eval_All_Exons.sh ${1} FreePSI_29
./eval_All_Exons.sh ${1} FreePSI_27_27
./eval_All_Exons.sh ${1} FreePSI_27_28
./eval_All_Exons.sh ${1} FreePSI_27_29
./eval_All_Exons.sh ${1} FreePSI_27_30

cat ${EVALUATION}/* > ${EVALUATION}/all.out

