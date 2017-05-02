#!/bin/bash

set -x
set -e

if [ $# != 2 ] ; then 
    echo "Require reads directory and thread number" 
    exit 1; 
fi 

THREAD=${2}
READS=${1}
OUTPUT=output/${1}_Hisat_MISO
RESULT=result/${1}_Hisat_MISO
TIME=${RESULT}/time_usage.csv

INDEX_DIR=../data_RefSeq_hg38/annotation/hisat_index
SPLICE_DIR=../data_RefSeq_hg38/annotation/splicing_event/gff/commonshortest
MISO_SETTING=${OUTPUT}/miso_setting.txt

rm -rf ${OUTPUT} ${RESULT}
mkdir ${OUTPUT} ${RESULT}

echo "[data]" > ${MISO_SETTING}
echo "filter_results = True" >> ${MISO_SETTING}
echo "min_event_reads = 20" >> ${MISO_SETTING}
echo "[sampler]" >> ${MISO_SETTING}
echo "burn_in = 500" >> ${MISO_SETTING}
echo "lag = 10" >> ${MISO_SETTING}
echo "num_iters = 5000" >> ${MISO_SETTING}
echo "num_processors = ${THREAD}" >> ${MISO_SETTING}

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    hisat2 --fr -p ${THREAD} \
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
    samtools index \
    ${OUTPUT}/mapped_reads_sorted.bam \
    ${OUTPUT}/mapped_reads_sorted.bam.bai

splicing_event=("A3SS" "A5SS" "MXE" "RI" "SE")
echo "" > ${OUTPUT}/miso_summary.out
for name in ${splicing_event[@]}
do
    /usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
        index_gff --index \
        ${SPLICE_DIR}/${name}.hg38.gff3 \
        ${OUTPUT}/gff_index_${name}

    /usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
        miso --run \
        ${OUTPUT}/gff_index_${name} \
        ${OUTPUT}/mapped_reads_sorted.bam \
        --settings-filename=${MISO_SETTING} \
        -p ${THREAD} \
        --read-len 76 \
        --output-dir ${OUTPUT}/miso_output_${name}

    /usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
        summarize_miso \
        --summarize-samples ${OUTPUT}/miso_output_${name} \
        ${OUTPUT}/miso_summary_${name}

    cat ${OUTPUT}/miso_summary_${name}/summary/miso_output_${name}.miso_summary \
    >> ${OUTPUT}/miso_summary.out
done

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    python3 ../scripts/computePsiFromSplicing.py \
    ${OUTPUT}/miso_summary.out \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_exonBoundary.bed \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_exonProfile.json \
    ../data_RefSeq_hg38/annotation/exon_boundary/hg38_refGene_nameMap.json \
    ${RESULT}/psi_Hisat_MISO.json
