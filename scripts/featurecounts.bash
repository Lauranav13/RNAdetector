#!/bin/bash

##############################################################################
# Options:
# 	-a FILE GTF
#		-b INPUT BAM (BAM_PATH or ALIGNMENT_PATH)
# 	-t NUMBER OF THREADS
#		-o OUTPUT
##############################################################################
while getopts ":a:b:t:o:" opt; do
	case $opt in
		a ) GTF_FILE=$OPTARG ;;
		b ) INPUT_BAM=$OPTARG ;;
		t ) THREADS=$OPTARG ;;
		o ) OUTPUT=$OPTARG ;;
		\?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
		: ) echo "Option -$OPTARG requires an argument." >&2; exit 2;;
	esac
done

#### Check parameters ####
# Check GTF annotation files
if [ -z $GTF_FILE ] || [ ! -f $GTF_FILE ]; then
	echo "Annotation file does not exist!"
	exit 3
fi

# Check input files
if [ -z $INPUT_BAM ] || [ ! -f $INPUT_BAM ]; then
	echo "Input file does not exist!" >&2
	exit 4
fi

# Check number of threads and set 1 as default value
if [ -z $THREADS ]; then
	THREADS=1
fi

# Check output
if [ -z $OUTPUT ]; then
	echo "Output file must be specified!" >&2
	exit 5
fi

# Check if output directory is writable
if [ ! -w $(dirname $OUTPUT) ]; then
	echo "Output directory is not writable!" >&2
	exit 6
fi

#### Counting ####
SAMPLE_NAME=$(basename $INPUT_BAM ".bam")
SUFF="_fc.txt"
OUTPUT_NAME=$SAMPLE_NAME$SUFF
/mnt/d/tool/subread-1.6.4-Linux-x86_64/bin/featureCounts -T $THREADS -a $GTF_FILE -o $OUTPUT/$OUTPUT_NAME $INPUT_BAM