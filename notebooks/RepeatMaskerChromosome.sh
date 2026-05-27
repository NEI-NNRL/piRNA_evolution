#!/bin/bash
outDIR=$1
fastaFile=$2
taxonID=$3
#cd $outDIR
cd /lscratch/${SLURM_JOB_ID}/ 
ml repeatmasker/4.2.1
echo "1) We are in $(pwd)"
echo $(ls)
echo "HOME=/lscratch/${SLURM_JOB_ID}/ RepeatMasker -dir /lscratch/${SLURM_JOB_ID}/ -engine rmblast -pa 4  -align -gff  -species $taxonID  $fastaFile"
HOME=/lscratch/${SLURM_JOB_ID}/ RepeatMasker -dir /lscratch/${SLURM_JOB_ID}/ -engine rmblast -pa 4 -align -gff  -species $taxonID  $fastaFile
echo "2) We are in $(pwd)"
echo $(ls)
#cd /lscratch/${SLURM_JOB_ID}/ 
RMputFile=$(ls | grep -e "out$")
mkdir -p $outDIR
cp $RMputFile $outDIR
echo "RMputFile = $RMputFile"
TEfileBED=${fastaFile}.repeatMasker.bed
TEfileBED=$(basename $TEfileBED)
TEfileBED=${outDIR}/${TEfileBED}
echo "TEfileBED =  $TEfileBED"
cat $RMputFile | grep  -e "CMC" -e "Chapaev" -e "Crypton" -e "Dada" -e "Ginger" -e "Harbinger" -e "hAT" -e "Blackjack" -e "Charlie" -e "IS3EU" -e "Kolobok" -e "Maverick" -e "Merlin" -e "MULE" -e "PiggyBac" -e "Sola" -e "LINE" -e "L1" -e "L2" -e "Penelope" -e "RTE" -e "LTR" -e "Copia" -e "DIRS" -e "ERV" -e "Gypsy" -e "PLE" -e "Helitron" -e "MIR" -e "SINE" -e "TcMar" -e "Tigger" -e "MuDR" -e "Zisupton" -e "Jockey" -e "CR1" -e "Dong" | grep -v ARTEFACT | awk '{ s = ($9 == "+" ? "+" : "-"); if ($9 != "+" && $9 !="C"){ exit 1}; print $5"\t"$6"\t"$7"\t"$11"\t0\t"s}' >  ${TEfileBED}
