#!/bin/bash
inFname=$1
workdir=$(dirname $)/
ml muscle
ml hmmer
cd $workdir
muscle -align ${inFname} -output ${inFname}.MSA.faa
hmmbuild ${inFname}.MSA.faa.hmm ${inFname}.MSA.faa