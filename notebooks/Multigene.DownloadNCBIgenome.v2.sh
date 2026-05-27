#!/bin/bash
sleep $3 # to avoid Error 439 too many requests
gAcc=$1
echo $gAcc
BaseDir=$2
FinalDir=${BaseDir}/${gAcc}/
echo $FinalDir
mkdir -p $FinalDir
cd /lscratch/${SLURM_JOB_ID}/ 
~/scripts/datasets download genome accession $gAcc --include gff3,rna,cds,protein,genome,seq-report
unzip ncbi_dataset.zip
for f in $(find ./ncbi_dataset/ -type f )
do
echo $f
cp $f $FinalDir
done
