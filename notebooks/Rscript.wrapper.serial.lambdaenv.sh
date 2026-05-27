#!/bin/bash
# READS file or R commands (args in the file)
source myconda
eval "$(mamba shell hook --shell bash)"
mamba activate lambdaenv
commandfile=$1
lowerlim=$2
upperlim=$3
let "numlines=$upperlim - $lowerlim +1"
echo "$lowerlim - $upperlim"
echo "numlines=$numlines"

cat $commandfile | head -n $upperlim | tail -n $numlines | while IFS= read -r line; do
    echo "Rscript $line"
    Rscript $line
done
