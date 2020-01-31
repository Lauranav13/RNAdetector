#!/bin/bash

##############################################################################
# Options:
#   -a input GTF/GFF formatted annotation file name (optional)
#   -i input SAM file name (required; generated by BWA-MEM)
#   -t number of threads (default: 1)
#   -o output (required)
#   -f FASTA reference genome.
#   -m max spanning distance (default: 200,000)
# FLAGS: -p for paired sequencing -1 to use CIRI1 instead of CIRI2
#   -h HARMONIZED output file
#   -b BED file for annotation
##############################################################################
PAIRED=false
VERSION=v2
while getopts "p1a:i:t:o:f:m:h:" opt; do
  case $opt in
  p) PAIRED=true ;;
  1) VERSION=v1 ;;
  a) GTF_FILE=$OPTARG ;;
  i) INPUT_SAM_FILE=$OPTARG ;;
  t) THREADS=$OPTARG ;;
  o) OUTPUT=$OPTARG ;;
  f) FASTA_FILE=$OPTARG ;;
  m) SPANNING=$OPTARG ;;
  h) HARMONIZED=$OPTARG ;;
  b) BED_FILE=$OPTARG ;;
  \?)
    echo "Invalid option: -$OPTARG"
    exit 1
    ;;
  :)
    echo "Option -$OPTARG requires an argument."
    exit 2
    ;;
  esac
done

#### Check parameters ####
# Check GTF annotation files
if [ -z "$GTF_FILE" ] || [ ! -f "$GTF_FILE" ]; then
  echo "Annotation file does not exist!"
  exit 3
fi

# Check input files
if [ -z "$INPUT_SAM_FILE" ] || [ ! -f "$INPUT_SAM_FILE" ]; then
  echo "Input file does not exist!"
  exit 4
fi

# Check number of threads and set 1 as default value
if [ -z "$THREADS" ]; then
  THREADS=1
fi

# Check output
if [ -z "$OUTPUT" ]; then
  echo "Output file must be specified!"
  exit 6
fi

# Check if output directory is writable
if [ ! -w "$(dirname "$OUTPUT")" ]; then
  echo "Output directory is not writable!"
  exit 7
fi

# Check fasta file
if [ -z "$FASTA_FILE" ] || [ ! -f "$FASTA_FILE" ]; then
  echo "FASTA file does not exist!"
  exit 8
fi

# Check max spanning distance (default: 200,000)
if [ -z "$SPANNING" ]; then
  SPANNING=200000
fi

#### circRNA identification ####

if [ "$VERSION" = "v1" ]; then
  if [ "$PAIRED" = "true" ]; then
    if ! perl /usr/local/bin/CIRI1.pl -P -I "$INPUT_SAM_FILE" -A "$GTF_FILE" -F "$FASTA_FILE" -M "$SPANNING" -O "$OUTPUT"; then
      echo "CIRI returned non zero exit code!"
      exit 5
    fi
  else
    if ! perl /usr/local/bin/CIRI1.pl -S -I "$INPUT_SAM_FILE" -A "$GTF_FILE" -F "$FASTA_FILE" -M "$SPANNING" -O "$OUTPUT"; then
      echo "CIRI returned non zero exit code!"
      exit 5
    fi

  fi
else
  if ! perl /usr/local/bin/CIRI2.pl -I "$INPUT_SAM_FILE" -A "$GTF_FILE" -F "$FASTA_FILE" -S "$SPANNING" -O "$OUTPUT" -T $THREADS; then
    echo "CIRI returned non zero exit code!"
    exit 5
  fi
fi

# Check SAM file
if [ ! -f "$OUTPUT" ]; then
  echo "Unable to find CIRI output file!"
  exit 9
fi

chmod 777 "$OUTPUT"

if [ ! -z "$HARMONIZED" ]; then
  CURR_DIR=$(pwd)
  SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
  cd $CURR_DIR
  if [ ! -z "$BED_FILE" ]; then
    if ! Rscript "${SCRIPT_PATH}/harmonize.R" -i "$OUTPUT" -a "ciri" -o "$HARMONIZED" -a "$BED_FILE"; then
      echo "Unable to harmonize output file"
      exit 10
    fi
  else
    if ! Rscript "${SCRIPT_PATH}/harmonize.R" -i "$OUTPUT" -a "ciri" -o "$HARMONIZED"; then
      echo "Unable to harmonize output file"
      exit 10
    fi
  fi
  chmod 777 "$HARMONIZED"
fi
