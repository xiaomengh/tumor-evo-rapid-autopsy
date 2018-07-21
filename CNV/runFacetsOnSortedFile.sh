#!/bin/bash


INPUT=$1
SAMPLE=${INPUT%_*}
mkdir $SAMPLE
cd $SAMPLE

module load R
RUNFACETS=runFacets.R
echo pt=\"$SAMPLE\" > temp.r
cat temp.r $RUNFACETS > temp1.r
R CMD BATCH --no-save temp1.r
rm temp.r
rm temp1.r
