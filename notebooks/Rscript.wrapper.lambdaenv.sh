#!/bin/bash
source myconda
eval "$(mamba shell hook --shell bash)"
mamba activate lambdaenv
echo Rscript $1 $2 $3 $4 $5 $6 $7 $8 $9
Rscript $1 $2 $3 $4 $5 $6 $7 $8 $9
