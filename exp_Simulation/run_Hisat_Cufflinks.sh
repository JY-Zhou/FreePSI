#!/bin/bash

set -x
set -e

if [ $# != 2 ] ; then 
    echo "Require reads directory and thread number" 
    exit 1; 
fi 

THREAD=${2}
READS=${1}
OUTPUT=output/${1}_Hisat_Cufflinks
RESULT=result/${1}_Hisat_Cufflinks
TIME=${RESULT}/time_usage.csv

INDEX_DIR=../data_RefSeq_hg38/annotation/hisat_index

rm -rf ${OUTPUT} ${RESULT}
mkdir ${OUTPUT} ${RESULT}

# Use Hisat_Cufflinks to calculate per-transcript TPMs.
/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    hisat2 --fr --dta-cufflinks -p ${THREAD} \
    -x ${INDEX_DIR}/index \
    -1 ${READS}/reads_final.1.fastq \
    -2 ${READS}/reads_final.2.fastq \
    | samtools view --threads ${THREAD} \
    -Sbo ${OUTPUT}/mapped_reads.bam

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    samtools sort \
    --threads ${THREAD} \
    -o ${OUTPUT}/mapped_reads_sorted.bam \
    ${OUTPUT}/mapped_reads.bam

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    cufflinks\
    -o ${OUTPUT}/transcriptome \
    -u -b ${INDEX_DIR}/index.fa \
    -p ${THREAD} --library-type fr-secondstrand \
    ${OUTPUT}/mapped_reads_sorted.bam

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    python3 ../scripts/computePsiFromAssembly.py \
    ${OUTPUT}/transcriptome/transcripts.gtf \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_exonBoundary.bed \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_exonProfile.json \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_nameMap.json \
    ${RESULT}/psi_Hisat_Cufflinks.json \
    ${RESULT}/tpm_Hisat_Cufflinks.json 
