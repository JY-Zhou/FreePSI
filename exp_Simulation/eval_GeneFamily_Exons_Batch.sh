#!/bin/bash

set -x
set -e

if [ $# != 1 ] ; then 
    echo "Require experiment name!" 
    exit 1; 
fi 

EVALUATION=evaluation/${1}_geneFamily

rm -rf ${EVALUATION}
mkdir ${EVALUATION}

./eval_GeneFamily_Exons.sh ${1} Hisat_Cufflinks_Annot
./eval_GeneFamily_Exons.sh ${1} Hisat_Cufflinks_Annot_withoutmul
./eval_GeneFamily_Exons.sh ${1} FreePSI
./eval_GeneFamily_Exons.sh ${1} Hisat_Cufflinks
./eval_GeneFamily_Exons.sh ${1} Salmon

cat ${EVALUATION}/* > ${EVALUATION}/all.out

