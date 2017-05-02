#!/bin/bash

set -x
set -e

if [ $# != 1 ] ; then 
    echo "Require reads directory!" 
    exit 1; 
fi 

READS=${1}
OUTPUT=output/${1}_Hisat_Cufflinks_Annot
RESULT=result/${1}_Hisat_Cufflinks_Annot
TIME=${RESULT}/time_usage.csv

INDEX_DIR=../data_RefSeq_hg38/annotation/hisat_index
GTF_FILE=../data_RefSeq_hg38/annotation/transcriptome_cleaned/hg38_refGene_cleaned.gtf

rm -rf ${OUTPUT} ${RESULT}
mkdir ${OUTPUT} ${RESULT}

# Use Hisat_Cufflinks to calculate per-transcript TPMs.
/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    hisat2 --fr --dta-cufflinks -p 16 \
    -x ${INDEX_DIR}/index \
    -1 ${READS}/reads_final.1.fastq \
    -2 ${READS}/reads_final.2.fastq \
    | samtools view --threads 16 \
    -Sbo ${OUTPUT}/mapped_reads.bam

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    samtools sort \
    --threads 16 \
    -o ${OUTPUT}/mapped_reads_sorted.bam \
    ${OUTPUT}/mapped_reads.bam

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    cufflinks -u \
    -G ${GTF_FILE} \
    -u -b ${INDEX_DIR}/index.fa \
    -o ${OUTPUT}/transcriptome \
    -p 16 --library-type fr-secondstrand \
    ${OUTPUT}/mapped_reads_sorted.bam

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    python3 ../scripts/computePsiFromAssembly.py \
    ${OUTPUT}/transcriptome/transcripts.gtf \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_exonBoundary.bed \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_exonProfile.json \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_nameMap.json \
    ${RESULT}/psi_Hisat_Cufflinks_Annot.json \
    ${RESULT}/tpm_Hisat_Cufflinks_Annot.json 
