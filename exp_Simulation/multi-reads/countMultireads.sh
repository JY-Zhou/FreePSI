#!/bin/bash

cat ./region.out | while read line
do
    tot=`samtools view -F 0x0004 -c mapped_reads_sorted.bam "$line"`
    multi=`samtools view -f 0x0100 -c mapped_reads_sorted.bam "$line"`
    printf "$line:\t $tot \t $multi\n"
done
