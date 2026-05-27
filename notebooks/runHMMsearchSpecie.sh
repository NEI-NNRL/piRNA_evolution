#!/bin/bash
ml hmmer
SPECIE=$1
PROTEIN_DIR=$2
HMM_DIR=$3
P_THRESHOLD=$4
OUT_DIR=$5
INPUTDIR=${PROTEIN_DIR}/${SPECIE}/
for HMM_FILE in $( ls $HMM_DIR/ENSG*.faa.MSA.faa.hmm);
do
HMM_BASE_PATH=$( basename ${HMM_FILE})
hmmsearch --tblout ${OUT_DIR}/${HMM_BASE_PATH}.out --notextw --noali -E ${P_THRESHOLD} ${HMM_FILE} ${INPUTDIR}/protein.faa 
done