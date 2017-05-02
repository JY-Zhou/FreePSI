#!/bin/bash

set -x
set -e

if [ $# != 1 ] ; then 
    echo "Require reads directory!" 
    exit 1; 
fi 

READS=PC3E/RNAseq
OUTPUT=output/${1}_Hisat_MISO
RESULT=result/${1}_Hisat_MISO
TIME=${RESULT}/time_usage.csv

INDEX_DIR=../data_RefSeq_hg19/annotation/hisat_index
SPLICE_DIR=../data_RefSeq_hg19/annotation/splicing_event/official/hg19
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
echo "num_processors = 16" >> ${MISO_SETTING}

/usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
    hisat2 -p 16 \
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
    samtools index \
    ${OUTPUT}/mapped_reads_sorted.bam \
    ${OUTPUT}/mapped_reads_sorted.bam.bai

splicing_event=("A3SS" "A5SS" "MXE" "RI" "SE")
echo "" > ${OUTPUT}/miso_summary.out
for name in ${splicing_event[@]}
do
    /usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
        index_gff --index \
        ${SPLICE_DIR}/${name}.hg19.gff3 \
        ${OUTPUT}/gff_index_${name}

    /usr/bin/time -f "%C\nNatural time:\t%e\nCPU time:\t%U+%S\nMemory peak:\t%M\n" -o ${TIME} -a \
        miso --run \
        ${OUTPUT}/gff_index_${name} \
        ${OUTPUT}/mapped_reads_sorted.bam \
        --settings-filename=${MISO_SETTING} \
        -p 16 \
        --read-len 101 \
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
    ../data_RefSeq_hg19/annotation/exon_boundary/hg19_refGene_exonBoundary.bed \
    ../data_RefSeq_hg19/annotation/exon_boundary/hg19_refGene_exonProfile.json \
    ../data_RefSeq_hg19/annotation/exon_boundary/hg19_refGene_nameMap.json \
    ${RESULT}/psi_Hisat_MISO.json
