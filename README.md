Genomic and Transcriptomic Analysis for Primary and Multi-metastatic Tumor Samples 
============================

We anaylzed deep whole genome sequencing data and RNA sequencing data of 30 tumors 
from a metastatic breast cancer patient. Here are the scripts that we used for analyses. 
Results are published in manuscript "Novel temporal and spatial patterns of metastatic 
colonization from rapid-autopsy tumor biopsies" (https://doi.org/10.1101/2021.06.03.446803). 

We divided the analyses to the following 5 sections:

  1. CNV analysis from WGS
  2. Phylogeny inference
  3. Subclone reconstruction
  4. Mutational Signature analysis
  5. RNAseq analysis

In each section, we provide necessary input files (except controlled access bam files) 
and scripts to perform the analysis.
  
CNV analysis from WGS
---------------------
See the `README.md` file in the CNV analysis folder

Phylogeny inference
-------------------
We used two methods to study the phylogenetic relationships among tumor
samples, which were described in the manuscript.  We provide `phylogeny.R`
script and an input vcf file to generate an allele frequency heatmap for all
samples, as well as to perform UPGMA clustering based on hamming distance.  We
also provide `seederSeeker.py` script and an input txt file (contains allele
frequency values) to study the phylogenetic relationships among tumor samples
based on the perfect phylogeny algorithm.  There are two parameters in the
script: one is a `threshold` for allele frequency in a sample to determine
whether a variant is present or absent; the other one is a `cutoff` for the
minimun number of variants that have the same sample distribution. 

The script can be invoked in the following manner

`$ cat strictSomatic_correctAF.txt | python3 seederSeeker.py > output.dot`

The output file can be visualized by graphvis

Subclone Reconstruction
-----------------------
We used `SubcloneSeeker V2` to generate subclone structures. We refer you to
the [SubcloneSeeker](https://github.com/yiq/SubcloneSeeker/tree/v2) repository
for more information. We provide input files, and a `Makefile` to perform all
analysis once SubcloneSeeker is installed

Mutational Signature
--------------------
We provide `ra_mutationpatterns.R` script to perform mutational signature
analysis. We provide input vcf files for the above R script. 

RNAseq analysis
---------------
We provide `CNV_RNAseq.R` script to perform copy number inference from RNAseq
data. In this R script, we provide the function for PCA and clustering
analysis.

We provide `edgeR_analysis.R` for differentially expressed gene analysis. 

The input files required for above analysis are included. 
