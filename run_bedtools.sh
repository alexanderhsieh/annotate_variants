#!/bin/bash
## USAGE: sh run_bedtools.sh <variants file>
## Purpose: Preprocessing filter to handle the following:
##			mappability, segmental duplication, LCR
## Dependencies: LCR, mapapbility, and segmental duplication tracks from UCSC genome browser

########################################
## Parse arguments and check files
########################################
# usage message
usage="
	## Usage:
	-v (required) - .txt file containing unfiltered variants
	-b (required) - path to bedtools
	-m (required) - path to mappability BED file
	-s (required) - path to segdup BED file
	-l (required) - path to LCR bed file
	-p (required) - path to parse_bedtools_isec.py script
	-H (flag) - echo this message and exit
"

# parse arguments
while getopts v:b:m:s:l:p:H opt; do 
	case "$opt" in 
		v) VAR="$OPTARG";;
		b) BED="$OPTARG";;
		m) MAP="$OPTARG";;
		s) SEG="$OPTARG";;
		l) LCR="$OPTARG";;
		p) PAR="$OPTARG";;
		H) echo "$usage"; exit;;
		\?) echo "$usage"; exit;;
		:) echo "$uage"; exit;;
	esac
done 

# check for arguments
if [[ -z $VAR ]] | [[ -z $BED ]] | [[ -z $MAP ]] | [[ -z $SEG ]] | [[ -z $LCR ]] | [[ -z $PAR ]]; then
	echo "$usage"; exit
fi

########################################
## Get output file prefix - everything before .txt
########################################
OUTPREFIX=`basename $VAR | sed s/.txt//`

########################################
## Format input variants as BED file and sort
## NOTE: assumes 1st col is id, 2nd = chr, 3rd = pos
########################################
echo "Formatting input variants BED file"
awk -F '\t' '{if($1!="id") print "chr"$2"\t"$3"\t"$3}' $VAR | sort -k1,1 -k2,2n > tmp.bed
echo ""

########################################
## Run bedtools intersection
########################################
echo "Running BEDtools intersect"
"$BED"/intersectBed -wa -wb -a tmp.bed -b $LCR $MAP $SEG -filenames > bed.isec.out.txt
echo ""
########################################
## Parse intersection file and append columns to input variants file
########################################
echo "Parsing BEDtools results and adding new columns in variants file"
python $PAR $VAR bed.isec.out.txt > $OUTPREFIX".RR.txt"
echo ""

echo "DONE"