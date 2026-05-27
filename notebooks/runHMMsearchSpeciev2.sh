#!/bin/bash
ml hmmer
SPECIE=$1
PROTEIN_DIR=$2
HMM_LIST=$3
P_THRESHOLD=$4
OUT_DIR=$5
INPUTDIR=${PROTEIN_DIR}/${SPECIE}/
mkdir -p $OUT_DIR
for HMM_FILE in $( cat $HMM_LIST );
do
HMM_BASE_PATH=$( basename ${HMM_FILE})
echo "hmmsearch --tblout ${OUT_DIR}/${HMM_BASE_PATH}.out --notextw --noali -E ${P_THRESHOLD} ${HMM_FILE} ${INPUTDIR}/protein.faa"
hmmsearch --tblout ${OUT_DIR}/${HMM_BASE_PATH}.out --notextw --noali -E ${P_THRESHOLD} ${HMM_FILE} ${INPUTDIR}/protein.faa 
done
