#!/bin/bash

set -x
set -e

if [ $# != 1 ] ; then 
    echo "Require reads directory!" 
    exit 1; 
fi 

READS=PC3E/RNAseq
OUTPUT=output/${1}_Hisat_Cufflinks_Annot
RESULT=result/${1}_Hisat_Cufflinks_Annot
TIME=${RESULT}/time_usage.csv

INDEX_DIR=../data_RefSeq_hg19/annotation/hisat_index
GTF_FILE=../data_RefSeq_hg19/annotation/transcriptome_cleaned/hg19_refGene_cleaned.gtf

rm -rf ${OUTPUT} ${RESULT}
mkdir ${OUTPUT} ${RESULT}

# Use Hisat_Cufflinks_Annot to calculate per-transcript TPMs.
/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    hisat2 --dta-cufflinks -p 16 \
    -x ${INDEX_DIR}/index \
    -U ${READS}/${1}_1.fastq,${READS}/${1}_2.fastq \
    | samtools view --threads 16 \
    -Sbo ${OUTPUT}/mapped_reads.bam

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    samtools sort \
    --threads 16 \
    -o ${OUTPUT}/mapped_reads_sorted.bam \
    ${OUTPUT}/mapped_reads.bam

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    cufflinks -p 16 \
    -G ${GTF_FILE} \
    -u -b ${INDEX_DIR}/index.fa \
    -o ${OUTPUT}/transcriptome \
    ${OUTPUT}/mapped_reads_sorted.bam

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    python3 ../scripts/computePsiFromAssembly.py \
    ${OUTPUT}/transcriptome/transcripts.gtf \
    ../data_RefSeq_hg19/annotation/exon_boundary/hg19_refGene_exonBoundary.bed \
    ../data_RefSeq_hg19/annotation/exon_boundary/hg19_refGene_exonProfile.json \
    ../data_RefSeq_hg19/annotation/exon_boundary/hg19_refGene_nameMap.json \
    ${RESULT}/psi_Hisat_Cufflinks_Annot.json \
    ${RESULT}/tpm_Hisat_Cufflinks_Annot.json 
