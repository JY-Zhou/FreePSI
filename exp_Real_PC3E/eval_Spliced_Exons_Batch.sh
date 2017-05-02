#!/bin/bash

set -x
set -e

if [ $# != 1 ] ; then 
    echo "Require experiment name!" 
    exit 1; 
fi 

EVALUATION=evaluation/${1}

rm -rf ${EVALUATION}
mkdir ${EVALUATION}

./eval_Spliced_Exons.sh ${1} Salmon
./eval_Spliced_Exons.sh ${1} Hisat_Cufflinks
./eval_Spliced_Exons.sh ${1} Hisat_Cufflinks_Annot
./eval_Spliced_Exons.sh ${1} FreePSI
./eval_Spliced_Exons.sh ${1} Hisat_MISO

cat ${EVALUATION}/* > ${EVALUATION}/all.out
