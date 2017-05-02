# FreePSI
An alignment-free approach to estimating exon-inclusion ratios without a reference transcriptome

## Brief introduction
FreePSI is a new method for genome-wide percent spliced in (PSI) estimation that requires neither a reference transcriptome (hence, transcriptome-free) nor the mapping of RNA-seq reads (hence, alignment-free).
The first freedom allows FreePSI to work effectively when a high quality reference transcriptome is unavailable and the second freedom not only helps make FreePSI more efficient, it also eliminates the necessity of dealing with multi-reads, which is a challenging problem by itself.
Note that this is the first alignment-free method in RNA-seq data analysis that does not require a reference transcriptome.
An outline of the method is sketched below.

FreePSI takes as the input a reference genome with exon boundary annotation and a set of RNA-seq reads.
Since a reference transcriptome is not assumed, it uses a weighted directed bipartite graph (called an abundance flow graph) to represent all possible isoforms of a gene and their expression levels.
In such a graph, each vertex represents an exon boundary and each edge represents either an exon or an exon junction.
The weight of an edge represents the total relative abundance of all isoforms covering the corresponding exon or junction.
Obviously, to estimate the PSI value of each exon, it suffices to infer the edge weights in every abundance flow graph.
By regarding each edge as a sequence of k-mers, FreePSI constructs a novel probabilistic model for generating all observed k-mers in the input RNA-seq reads based on the abundance flow graphs for all genes.
It then employs the expectation-maximization framework and a divide-and-conquer strategy to decompose genome-wide maximum likelihood estimation into independent optimization subproblems for each gene, which is then solved by an ultrafast algorithm, conjugate gradient projection descent.
Finally, it uses a post-processing procedure based on straightforward correlation analysis to 'smooth out' the PSI values in each gene.
FreePSI could have important applications in alternative splicing analysis when a high quality reference transcriptome is unavailable.

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
`freePSI [options] -g <GENOME_DIR> -a <EXON_BND_ANNOT> [-r <KMER_COUNT> | -1 <KMER_COUNT> -2 <KMER_COUNT>] -o <OUTPUT>`

Options:

  `-g <GENOME_DIR>`                  The directory containing the reference genome of each chromosome (.fasta format)
  
  `-a <EXON_BND_ANNOT>`              The annotation of exon boundary (.bed format)
  
  `-r <KMER_COUNT>`                  K-mer counts produced by Jellyfish in single-end reads (.fasta format)
  
  `-1 <KMER_COUNT> -2 <KMER_COUNT>`  K-mer counts produced by Jellyfish of paired-end reads  (.fasta format)
  
  `-o <OUTPUT> `                     The result of PSI values
  
  `-k [Integer] `                    The length of k-mer (Default: 27)
  
  `-p [Thread number]`               The thread numbers (Default: 1)
  
  `-h `                              Help information

### Post-processing
`python3 postProc.py <RAW_OUTPUT> <REFINED_OUTPUT>`

### Summary output
`python3 summary.py <EXON_BND_ANNOT> <REFINED_OUTPUT> <SUMMARY>`
