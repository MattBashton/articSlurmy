#!/usr/bin/bash

# Matt Bashton 2020
# Generating .tsv stats for ARTIC network analysis run
# Run this script on *.o* standard out files produced in your output/ dir

tput bold
echo "Matt Bashton 2020"
echo "Generate .tsv stats for ARTIC network analysis run"
echo ""
tput sgr0

date
hostname

# Get all STD out files
shopt -s nullglob
OUTFILES=(*.o*)
echo ""
echo "Found: ${#OUTFILES[@]} run standard out files to parse"

function get_stats {

    local FILE=${1}

    SAMPLE_ID=$(grep -oP '^\s-\sSample\sID\sfor\soutput:\s\K\S+' ${FILE})
    BARCODE=$(grep -oP '^\s-\sBarcodes\sfor\sthis\srun:\s\K\S+' ${FILE})
    IN_READS=$(grep -oP '^Total\sinput\sreads\spre\sguppyplex:\s\K\S+' ${FILE})
    OUT_READS=$(grep -oP '^Reads\spost\sguppyplex:\s\K\S+' ${FILE} | head -n 1)
    PC_R=$(grep -oP '^Reads\spost\sguppyplex:\s\d+\s\(\K\d+\.\d+' ${FILE} | head -n 1)
    IN_DEPTH=$(grep -oP '^Mean\sdepth\sof\sinput\salignment:\s\K\S+' ${FILE})
    OUT_DEPTH=$(grep -oP '^Mean\sdepth\sof\sfinal\salignment:\s\K\S+' ${FILE})
    N=$(grep -oP '^Number\sof\sNs\sconsensus:\s\K\d+' ${FILE})
    PC_N=$(grep -oP '^Number\sof\sNs\sconsensus:\s\d+\s\(\K\d+\.\d+' ${FILE})
    VARS=$(grep -oP '^Number\sof\svariants\scalled:\s\K\d+' ${FILE} | head -n 1)
    LINEAGE=$(grep -A 1 '^taxon' ${FILE} | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)

    echo -e "${SAMPLE_ID}\t${BARCODE}\t${IN_READS}\t${OUT_READS}\t${PC_R}\t${IN_DEPTH}\t${OUT_DEPTH}\t${N}\t${PC_N}\t${VARS}\t${LINEAGE}"

}

for FILE in ${OUTFILES[@]}
do
    get_stats ${FILE} >> run_stats.tmp
done

# Write out header line to file
echo -e "central_sample_id\tbarcode\ttotal_reads\tpost_guppyplex_reads\treads_carried_into_analysis\tinput_depth\touput_depth\tnumber_consensus_Ns\t%N\tvariants_called\tlineage" > run_stats.tsv

# Sort our ouptut by barcode
sort -k 2 run_stats.tmp >> run_stats.tsv
# rm tmp file
rm run_stats.tmp

echo ""
echo "Stats for run:"
echo ""
column -t run_stats.tsv

echo ""
echo "Done: output written to run_stats.tsv"
