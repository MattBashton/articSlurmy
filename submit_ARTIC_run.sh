#!/usr/bin/bash -eu

# Matt Bashton 2020
# Submit script for ARTIC network analysis

tput bold
echo "Matt Bashton 2020"
echo "Submit script for ARTIC network analysis, submits an array job for an ONT run,"
echo "each barcode is a task in the array, requires barcode - sample ID mapping file."
echo ""
tput sgr0

TIME=$(date)
HOST=$(hostname)

# We need a sample sheet to run
[ $# -ne 5 ] && { echo -en "Error, usage: $(basename $0) <tab delimited barcode dir - sample ID mapping> <input run dir> <sequencing_summary.txt> <output dir for analysis> <global run name>\n\n" ; exit 1; }

# Get our input arguments
BARCODE_MAPPING="$1"
INPUT_DIR=${2%/}
SEQ_SUM="$3"
OUTPUT_DIR=${4%/}
RUN_NAME="$5"

# Check files exist
if [[ ! -f ${BARCODE_MAPPING} ]] ; then
    echo "Error barcode - sample ID mapping file: ${BARCODE_MAPPING} does not exist."
    exit
fi

if [[ ! -f ${SEQ_SUM} ]] ; then
    echo "Error sequencing_summary.txt file: ${SEQ_SUM} does not exist."
    exit
fi

if [[ ! -d ${INPUT_DIR} ]] ; then
    echo "Error input run dir: ${INPUT_DIR} does not exist."
    exit
fi

if [[ -d ${OUTPUT_DIR} ]] ; then
    echo "Error output analysis run dir: ${OUTPUT_DIR} exists! - Will not overwrite"
    exit
fi

# Create output dir
mkdir ${OUTPUT_DIR}

# Size run from barcode mapping file
N=$(wc -l ${BARCODE_MAPPING} | awk '{print $1}')

# List variables
echo "** Variables **"
echo " - Host: ${HOST}"
echo " - Current directory: ${PWD}"
echo " - Time: ${TIME}"
echo " - Barecode - sample ID mapping file: ${BARCODE_MAPPING}"
echo " - Sequencing summary file: ${SEQ_SUM}"
echo " - Input run dir: ${INPUT_DIR}"
echo " - Ouput dir: ${OUTPUT_DIR}"
echo " - Global run name: ${RUN_NAME}"
echo " - Number of barcodes to process: ${N}"
echo ""

# Copy script to output dir
echo "Copying artic.sh to ${OUTPUT_DIR}"
cp -v artic.sh ${OUTPUT_DIR}
echo "Copying getStats.sh to ${OUTPUT_DIR}, (run after array job is finished)"
cp -v getStats.sh ${OUTPUT_DIR}

# Submit array
echo "Submitting array job of ${N} tasks,"
JOB_ID1=$(sbatch --array 1-${N} --job-name artic --parsable -D ${OUTPUT_DIR} ${OUTPUT_DIR}/artic.sh ${PWD}/${BARCODE_MAPPING} ${PWD}/${SEQ_SUM} ${PWD}/${INPUT_DIR} ${PWD}/${OUTPUT_DIR} ${RUN_NAME})
echo "Array job ID is: ${JOB_ID1}"
