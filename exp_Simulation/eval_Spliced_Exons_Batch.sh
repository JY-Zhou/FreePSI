#!/bin/bash

set -x
set -e

if [ $# != 1 ] ; then 
    echo "Require experiment name!" 
    exit 1; 
fi 

EVALUATION=evaluation/${1}_spliced

rm -rf ${EVALUATION}
mkdir ${EVALUATION}

./eval_Spliced_Exons.sh ${1} Salmon
./eval_Spliced_Exons.sh ${1} Hisat_Cufflinks
./eval_Spliced_Exons.sh ${1} FreePSI
./eval_Spliced_Exons.sh ${1} Hisat_MISO
./eval_Spliced_Exons.sh ${1} Hisat_Cufflinks_Annot
./eval_Spliced_Exons.sh ${1} Hisat_Cufflinks_Annot_Sampling_90
./eval_Spliced_Exons.sh ${1} Hisat_Cufflinks_Annot_Sampling_80
./eval_Spliced_Exons.sh ${1} Hisat_Cufflinks_Annot_Sampling_70
./eval_Spliced_Exons.sh ${1} Salmon_Sampling_90
./eval_Spliced_Exons.sh ${1} Salmon_Sampling_80
./eval_Spliced_Exons.sh ${1} Salmon_Sampling_70
./eval_Spliced_Exons.sh ${1} FreePSI_21
./eval_Spliced_Exons.sh ${1} FreePSI_23
./eval_Spliced_Exons.sh ${1} FreePSI_25
./eval_Spliced_Exons.sh ${1} FreePSI_27
./eval_Spliced_Exons.sh ${1} FreePSI_29
./eval_Spliced_Exons.sh ${1} FreePSI_27_27
./eval_Spliced_Exons.sh ${1} FreePSI_27_28
./eval_Spliced_Exons.sh ${1} FreePSI_27_29
./eval_Spliced_Exons.sh ${1} FreePSI_27_30

cat ${EVALUATION}/* > ${EVALUATION}/all.out
