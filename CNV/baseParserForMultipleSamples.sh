#!/bin/bash

module load python/3.5.2

BASE_PARSER=baseParserForMultipleSamples.py

INPUT_FILE=$1
OUTFILE=${INPUT_FILE##*/}
OUTFILE=${OUTFILE%%.txt}.parsed.bed

cat $INPUT_FILE | python3 $BASE_PARSER > $OUTFILE

