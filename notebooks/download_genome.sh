#!/bin/bash
sleep $3
inName=$1
genomesDIR=$2
outDIR=$genomesDIR/${inName}/
mkdir -p $outDIR
cd $outDIR
~/scripts/datasets download genome accession ${inName} --include gff3,rna,cds,protein,genome,seq-report --filename ${outDIR}${inName}.zip
unzip ${inName}.zip
for f in $(find  ./ -type f | grep ncbi | grep data); do mv $f $(basename $f)  ; done;
rm -r ncbi_dataset/
