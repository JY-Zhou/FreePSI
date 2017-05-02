#!/bin/bash

set -o nounset
set -o errexit

# Create temporary directory for FluxSimulator
rm -rf flux_simulator_tmp
mkdir flux_simulator_tmp

# Now use Flux Simulator to simulate reads.
flux-simulator -p flux_simulator_main_simulation.par
rm -rf flux_simulator_tmp

mv main_reads.fastq reads.fastq

# Some isoform quantifiers require reads to be presented in a random order,
# hence we shuffle the reads output by Flux Simulator.
paste - - - - - - - - < reads.fastq | shuf | tr '\t' '\n' > reads.tmp
mv reads.tmp reads.fastq

# We've produced paired-end reads - split the Flux Simulator output into
# files containing left and right reads.
paste - - - - < reads.fastq | awk -F '\t' '$1~/\/1/ {print $0 > "lr.tmp"} $1~/\/2/ {print $0 > "rr.tmp"}'
rm reads.fastq
tr '\t' '\n' < lr.tmp > reads_final.1.fastq
rm lr.tmp
tr '\t' '\n' < rr.tmp > reads_final.2.fastq
rm rr.tmp

# Remove intermediate files not necessary for quantification.
#rm *.lib
#rm *.bed
