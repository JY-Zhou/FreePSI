# FreePSI
An alignment-free approach to estimating exon-inclusion ratios without a reference transcriptome

## Brief introduction
-

## Requirements
FreePSI requires [Jellyfish](https://github.com/gmarcais/Jellyfish/releases) (version 2.2.6 or higher) for counting k-mers in RNA-seq reads and [Python3](https://www.python.org/downloads/) (version 3.5.2 or higher) with [Scipy](https://www.scipy.org/) (version 0.18.1 or higher) and [Numpy](http://www.numpy.org/) (version 1.11.1 or higher) for post-processing and summarizing output.

## Installation
Currently, we provide a binary release [here](https://github.com/JY-Zhou/FreePSI/releases) (compiled with GCC version 5.4.0) for 64-bit Linux system.

Building FreePSI from source requires installing libraries [Eigen](http://eigen.tuxfamily.org) (version 3.3 or higher) and [Boost](http://www.boost.org) (version 1.61.0 or higher).
We will provide an installable version soon.

## Example
We provide an example for testing FreePSI.
The example contains the reference genome with exon boundary annotation of human chromosome 21 and corresponding simulated RNA-seq reads.
Please run the shell script `run_FreePSI.sh` in the folder `example` for a try.

## Usage
### Main procedure
#### K-mer hash table construction
Description:
  The `build` command builds the k-mer hash table from a reference genome and loads the k-mer counts of RNA-seq reads.

Usage:
  `freePSI build -g <GENOME_DIR> -a <EXON_BND_ANNOT> [-r <KMER_COUNT> | -1 <KMER_COUNT> -2 <KMER_COUNT>] -o <HASHTABLE> [options..]`

Options:
  `-g <GENOME_DIR>`                  The directory containing the reference genome of each chromosome (.fasta format)
  
  `-a <EXON_BND_ANNOT>`              The annotation of exon boundary (.bed format)
  
  `-r <KMER_COUNT>`                  The k-mer count of single-end reads produced by Jellyfish (.fasta format)
  
  `-1 <KMER_COUNT> -2 <KMER_COUNT>`  The k-mer count of paired-end reads produced by Jellyfish (.fasta format)
  
  `-o <HASHTABLE>`                   The k-mer hash table (.json format)
  
  `-k [Integer]`                     The length of k-mer (Default: 27)
  
  `-p [Thread number]`               The thread numbers. (Default: 1)
  
  `-h`                               Help information

#### PSI estimation
Description:
  The `quant` command quantifies the PSI values.

Usage:
  `freePSI quant -i <INDEX> -o <OUTPUT>`

Options:
  `-i <HASHTABLE>`                   The k-mer hash table (.json format)
  
  `-o <OUTPUT>`                      The result of PSI values (.json format)
  
  `-k [Integer]`                     The length of k-mer (Default: 27)
  
  `-p [Thread number]`               The thread numbers. (Default: 1)
  
  `-h`                               Help information


### Post-processing
`python3 postProc.py <RAW_OUTPUT> <REFINED_OUTPUT>`

### Summary output
`python3 summary.py <EXON_BND_ANNOT> <REFINED_OUTPUT> <SUMMARY>`
