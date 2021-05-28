#!/usr/bin/bash
#SBATCH --export=NONE
#SBATCH --mem=16G
#SBATCH --cpus-per-task 4
#SBATCH --time=48:00:00
# Generate SoG/SGE like .o and .e files
#SBATCH -o %x.o%A.%a
#SBATCH -e %x.e%A.%a
# Partition or queue to use
#SBATCH -p compute

# Needed for conda on our set-up
source ~/.bashrc
# Prevent Qt plotting issues with unattached X11 display
# No longer needed with artic v.1.2.0+, (uncomment for older versions of artic)
#DISPLAY=""

TIME=$(date)
HOST=$(hostname)
START_T=$(date +%s)

echo "Matt Bashton 2020"
echo "SLURM Bash script to run the Artic Network nCoV2 pipeline on ONT data"
echo ""

set -e
set -u
set -o pipefail

# Get our input arguments
BARCODE_MAPPING="${1}"
SEQ_SUM="${2}"
INPUT_DIR="${3}"
OUTPUT_DIR="${4}"
RUN_NAME="${5}"

# Parse barcode mapping
LINE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" ${BARCODE_MAPPING})
set ${LINE}
BARCODE="${1}"
SAMPLE="${2}"

# Some defaults for ONT data
# See https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
BASE_DIR=${PWD}
CONDA_ENV="artic-ncov2019"
CPU="4"
MIN_LEN="400"
MAX_LEN="700"
NORMALISE="200"
# This will be down to where you have artic-ncov2019 conda env installed
SCHEME_DIR="${HOME}/artic-ncov2019/primer_schemes"
PRIMERS="nCoV-2019/V3"
ANNOVAR="${HOME}/annovar/table_annovar.pl"
ANNOVAR_DB="${HOME}/annovar/sarscov2db"
AMPLICON_BED="${HOME}/articSlurmy/ampliconsV3.bed"

# Make tmpdir for run
# The /scratch /tmp /local prefix here might be different on your system
TMPDIR=$(mktemp -d -p /local)

# List variables
echo "** Variables **"
echo " - Host: ${HOST}"
echo " - Current working directory: ${PWD}"
echo " - Time: ${TIME}"
echo " - User home dir: ${HOME}"
echo " - Job name: ${SLURM_JOB_NAME}"
echo " - Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo " - Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo " - Tmp dir: ${TMPDIR}"
echo " - Global run name: ${RUN_NAME}"
echo " - Barecode - sample ID mapping file: ${BARCODE_MAPPING}"
echo " - Sequencing summary file: ${SEQ_SUM}"
echo " - Input run dir: ${INPUT_DIR}"
echo " - Ouput dir: ${OUTPUT_DIR}"
echo " - Barcodes for this run: ${BARCODE}"
echo " - Sample ID for output: ${SAMPLE}"
echo " - Threads to run artic minion: ${CPU}"
echo " - Min length for read filtering: ${MIN_LEN}"
echo " - Max length for read filtering: ${MAX_LEN}"
echo " - Guppyplex normalisation setting: ${NORMALISE}"
echo " - Primer scheme dir: ${SCHEME_DIR}"
echo " - Primer type/version: ${PRIMERS}"
echo " - Annovar: ${ANNOVAR}"
echo " - Annovar DB: ${ANNOVAR_DB}"
echo " - Amplicons bed file: ${AMPLICON_BED}"
echo ""

# Activate conda enviroment
echo "Activating Conda enviroment: ${CONDA_ENV}"
conda activate ${CONDA_ENV}
conda env list
echo -ne "Python interpreter located at: "
which python
python --version
echo ""

# cd to tmpdir dir, as output is streamed, keep it off Luster
echo "cd ing to ${TMPDIR}"
cd ${TMPDIR}

# guppyplex filtering output fastq is written to $TMPDIR, single file per barcode, this script works per barcode
echo "Running guppyplex to filter reads, min length = ${MIN_LEN}, max length = ${MAX_LEN}, "
echo "input is: ${INPUT_DIR}/fastq_pass/${BARCODE} prefixing output fastq with: ${RUN_NAME}"
READS_PRE=$(zcat ${INPUT_DIR}/fastq_pass/${BARCODE}/*.gz | awk '{s++}END{print s/4}')
echo "Total input reads pre guppyplex: ${READS_PRE}"
artic guppyplex --skip-quality-check --min-length ${MIN_LEN} --max-length ${MAX_LEN} --directory ${INPUT_DIR}/fastq_pass/${BARCODE} --prefix ${RUN_NAME}
READS_POST=$(awk '{s++}END{print s/4}' ${TMPDIR}/${RUN_NAME}_${BARCODE}.fastq)
PC_PLEX=$(echo "scale=8; (${READS_POST}/${READS_PRE})*100" | bc | xargs printf "%.2f\n")
echo "Reads post guppyplex: ${READS_POST} (${PC_PLEX}%)"

# artic minion final stage, many output / temp files written to ${TMPDIR}
echo "Running artic minion --normalise ${NORMALISE} on ${CPU} threads, fast5 dir is: ${INPUT_DIR}/fast5_pass"
artic minion --normalise ${NORMALISE} --threads ${CPU} --scheme-directory ${SCHEME_DIR} --read-file ${TMPDIR}/${RUN_NAME}_${BARCODE}.fastq --fast5-directory ${INPUT_DIR}/fast5_pass --sequencing-summary ${SEQ_SUM} ${PRIMERS} ${SAMPLE}

# Check output at this stage
echo ""
echo "Current output is:"
ls -lh

# Find number of aligned reads in ${SAMPLE}.sorted.bam
echo ""
ALN_READS=$(samtools view -c -F 2308 ${TMPDIR}/${SAMPLE}.sorted.bam)
echo "Aligned reads in .sorted.bam: ${ALN_READS}"

# Annotate with annovar
# Work out No variants called
VAR=$(zgrep -cv '^#' ${SAMPLE}.pass.vcf.gz || true)
echo ""
echo "Number of variants called: ${VAR}"

if (( ${VAR} > 0 )); then
    echo "Annotating with annovar located at: ${ANNOVAR} using annovar DB ${ANNOVAR_DB}"
    # Fix ref for annovar with sed
    zcat ${SAMPLE}.pass.vcf.gz | sed 's/MN908947.3/NC_045512v2/' > ${SAMPLE}.av.pass.vcf
    ${ANNOVAR} -buildver NC_045512v2 ${SAMPLE}.av.pass.vcf ${ANNOVAR_DB} -protocol avGene -operation g -vcfinput

    echo "Variants called:"
    cut -f 1-10 ${SAMPLE}.av.pass.vcf.NC_045512v2_multianno.txt | column -t

elif (( ${VAR} == 0 )); then
    echo "No variants called, skipping annovar"
fi

# Final run stats

# Work out number of Ns
TOTAL_BP=$(grep -v ">" ${TMPDIR}/${SAMPLE}.consensus.fasta | tr -d "\n\r" | wc -c)
TOTAL_N=$(grep -v ">" ${TMPDIR}/${SAMPLE}.consensus.fasta | grep -o 'N' | wc -l)
PC_N=$(echo "scale=8; (${TOTAL_N}/${TOTAL_BP})*100" | bc | xargs printf "%.2f\n")

# Per amplicon depth stats
mosdepth -t ${CPU} -b $AMPLICON_BED ${SAMPLE}.inputAmp ${TMPDIR}/${SAMPLE}.sorted.bam
mosdepth -t ${CPU} -b $AMPLICON_BED ${SAMPLE}.outputAmp ${TMPDIR}/${SAMPLE}.primertrimmed.rg.sorted.bam
echo ""
echo "Minimap2 per amplicon mean coverage:"
zcat ${SAMPLE}.inputAmp.regions.bed.gz | awk '{print NR "\t" $5}'
zcat ${SAMPLE}.inputAmp.regions.bed.gz | awk '{print NR "\t" $5}' > ${SAMPLE}.inputAmp.depth.tsv
echo ""
echo "Post nanopolish per amplicon mean coverage:"
zcat ${SAMPLE}.outputAmp.regions.bed.gz | awk '{print NR "\t" $5}'
zcat ${SAMPLE}.outputAmp.regions.bed.gz | awk '{print NR "\t" $5}' > ${SAMPLE}.outputAmp.depth.tsv

# Per amplicon N counts
# Fix header for bedtools
sed s'/>.*/>MN908947.3/' ${TMPDIR}/${SAMPLE}.consensus.fasta > ${TMPDIR}/temp.${SAMPLE}.MN908947.3.fasta
echo ""
echo "Per amplicon N counts:"
bedtools getfasta -fi ${TMPDIR}/temp.${SAMPLE}.MN908947.3.fasta -tab -bed ${AMPLICON_BED} | cut -f 2 | tr -d -c 'N\n' | awk '{ print NR "\t" length; }'
bedtools getfasta -fi ${TMPDIR}/temp.${SAMPLE}.MN908947.3.fasta -tab -bed ${AMPLICON_BED} | cut -f 2 | tr -d -c 'N\n' | awk '{ print NR "\t" length; }' > ${SAMPLE}.AmpNs.tsv

# Global depth stats
mosdepth -t ${CPU} ${SAMPLE}.input ${TMPDIR}/${SAMPLE}.sorted.bam
mosdepth -t ${CPU} ${SAMPLE}.output ${TMPDIR}/${SAMPLE}.primertrimmed.rg.sorted.bam
echo ""
echo "Minimap2 mean coverage stats:"
column -t ${SAMPLE}.input.mosdepth.summary.txt | head -n 2
IN_DEPTH=$(sed -sn 2p ${SAMPLE}.input.mosdepth.summary.txt | cut -f 4)
echo ""
echo "Post nanopolish mean coverage stats:"
column -t ${SAMPLE}.output.mosdepth.summary.txt | head -n 2
OUT_DEPTH=$(sed -sn 2p ${SAMPLE}.output.mosdepth.summary.txt | cut -f 4)
echo ""

# Copy back and rename output
# Not all of the output is needed, and files to upload will have to be renamed see:
# https://docs.covid19.climb.ac.uk/upload-instructions

# Dir for data to upload
# Needs alignment.bam and consensus.fa
echo "Creating upload dir at ${OUTPUT_DIR}/upload/${RUN_NAME}/${SAMPLE}"
mkdir -p ${OUTPUT_DIR}/upload/${RUN_NAME}/${SAMPLE}
echo "Copying upload data..."
cp ${TMPDIR}/${SAMPLE}.consensus.fasta ${OUTPUT_DIR}/upload/${RUN_NAME}/${SAMPLE}/consensus.fa
cp ${TMPDIR}/${SAMPLE}.sorted.bam ${OUTPUT_DIR}/upload/${RUN_NAME}/${SAMPLE}/alignment.bam
# Dir for all output
echo "Creating output dir at ${OUTPUT_DIR}/processed/${RUN_NAME}"
mkdir -p ${OUTPUT_DIR}/processed/${RUN_NAME}
echo "Copying ${SAMPLE}.* ..."
cp ${TMPDIR}/${SAMPLE}.* ${OUTPUT_DIR}/processed/${RUN_NAME}
# Needed for older artic v.1.1.3 which produced .pngs with a different output
# convention, commented out now as not produced in v.1.2.0 or higher
#echo "Copying ${SAMPLE}-* ..."
#cp ${TMPDIR}/${SAMPLE}-* ${OUTPUT_DIR}/processed/${RUN_NAME}

echo "Gzipping ${RUN_NAME}_${BARCODE}.fastq ..."
# Needs pigz, conda install pigz, into your artic-ncov2019 env
# Or replace with single threaded gzip, omitting -p ${CPU}
pigz -p ${CPU} ${TMPDIR}/${RUN_NAME}_${BARCODE}.fastq
echo "Copying ${RUN_NAME}_${BARCODE}.fastq.gz ..."
cp ${TMPDIR}/${RUN_NAME}_${BARCODE}.fastq.gz ${OUTPUT_DIR}/processed/${RUN_NAME}

# Run pangolin
echo ""
echo "Running Pangolin: Phylogenetic Assignment of Named Global Outbreak LINeages"
echo "Switching conda envs"
echo "Pangolin stdout written to: ${OUTPUT_DIR}/processed/${RUN_NAME}/${SAMPLE}.pangolin.out"
echo "Snakemake spam written to: ${OUTPUT_DIR}/processed/${RUN_NAME}/${SAMPLE}.pangolin.err"
echo ""
conda deactivate
conda activate pangolin
pangolin -t ${CPU} -o ${TMPDIR}/${SAMPLE}.pangolin --tempdir ${TMPDIR} ${TMPDIR}/${SAMPLE}.consensus.fasta > ${OUTPUT_DIR}/processed/${RUN_NAME}/${SAMPLE}.pangolin.out 2>${OUTPUT_DIR}/processed/${RUN_NAME}/${SAMPLE}.pangolin.err
cp ${TMPDIR}/${SAMPLE}.pangolin/lineage_report.csv ${OUTPUT_DIR}/processed/${RUN_NAME}/${SAMPLE}.pangolin.lineage_report.csv
echo "Lineage report:"
column -t -s ',' ${TMPDIR}/${SAMPLE}.pangolin/lineage_report.csv
echo ""

# Clean-up temp dir
echo "Removing temp dir: ${TMPDIR}"
cd ${BASE_DIR}
rm -rf ${TMPDIR}

echo ""
echo ""
echo "Final stats"
echo "==========="
echo ""
echo "Reads pre guppyplex: ${READS_PRE}"
echo "Reads post guppyplex: ${READS_POST} (${PC_PLEX}%)"
echo "Aligned reads in sorted.bam: ${ALN_READS}"
echo "Number of Ns consensus: ${TOTAL_N} (${PC_N}%)"
echo "Mean depth of input alignment: ${IN_DEPTH}"
echo "Mean depth of final alignment: ${OUT_DEPTH}"
echo "Number of variants called: ${VAR}"
echo ""

# Write these stats to a file for each sample in output dir
echo "Reads pre guppyplex: ${READS_PRE}" > ${OUTPUT_DIR}/processed/${RUN_NAME}/${SAMPLE}.stats.txt
echo "Reads post guppyplex: ${READS_POST}" >> ${OUTPUT_DIR}/processed/${RUN_NAME}/${SAMPLE}.stats.txt
echo "Reads post guppyples %: ${PC_PLEX}" >> ${OUTPUT_DIR}/processed/${RUN_NAME}/${SAMPLE}.stats.txt
echo "Aligned reads in sorted.bam: ${ALN_READS}" >> ${OUTPUT_DIR}/processed/${RUN_NAME}/${SAMPLE}.stats.txt
echo "Number of Ns consensus: ${TOTAL_N}" >> ${OUTPUT_DIR}/processed/${RUN_NAME}/${SAMPLE}.stats.txt
echo "Number of Ns consensus %: ${PC_N}" >> ${OUTPUT_DIR}/processed/${RUN_NAME}/${SAMPLE}.stats.txt
echo "Mean depth of input alignment: ${IN_DEPTH}" >> ${OUTPUT_DIR}/processed/${RUN_NAME}/${SAMPLE}.stats.txt
echo "Mean depth of final alignment: ${OUT_DEPTH}" >> ${OUTPUT_DIR}/processed/${RUN_NAME}/${SAMPLE}.stats.txt
echo "Number of variants called: ${VAR}" >> ${OUTPUT_DIR}/processed/${RUN_NAME}/${SAMPLE}.stats.txt

date
END_T=$(date +%s)
RUN_T=$((${END_T}-${START_T}))
echo -ne "Run time was: "
printf '%dh:%dm:%ds\n' $((${RUN_T}/3600)) $((${RUN_T}%3600/60)) $(($RUN_T%60))
echo ""
echo "END"
