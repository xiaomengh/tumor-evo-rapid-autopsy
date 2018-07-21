#!/bin/bash

module load bedtools/2.24.0
echo "start intersect job"
INPUT_FILE=$1
SNP_MAP= <PATH-TO-FILE>/InfiniumOmni2-5-8v1-3_A1.bed
OUTFILE=${INPUT_FILE##*/}
OUTFILE=${OUTFILE%.bed}.SNP_match.txt
bedtools intersect -wa -a $INPUT_FILE -b $SNP_MAP > $OUTFILE
echo "finished intersect job"
