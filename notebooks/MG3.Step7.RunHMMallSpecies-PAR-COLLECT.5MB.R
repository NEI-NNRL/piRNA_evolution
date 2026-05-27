argsT = commandArgs(trailingOnly=FALSE)
message(length(argsT))
library(stringr)
for (arg in argsT){
    if (str_count(arg, "--file=")>0){
        tmpStr=str_remove(arg, "--file=")
        sourceDir=dirname(tmpStr)
        break
    }
}
source(paste0(sourceDir,'/MG3.functions1.R'))
args = commandArgs(trailingOnly=TRUE)
HMMid=args[1]

message("HMMid: ", HMMid)
HMMidDir=paste0(OutputDir,"step7.everyHMMsearch5MB/PARALLEL/",HMMid,"/")
accessions=list.files(HMMidDir, pattern = ".*txt", full.names = FALSE)
accessions=gsub(pattern = ".txt", replacement = "", x = accessions, fixed = TRUE) 
bestGenesList=list()

for (accession in accessions){
    message(accession," ",getTimeStamp())
    bestGene=readLines(paste0(HMMidDir,accession,".txt"))
    if (nchar(bestGene)>0){
        bestGenesList[[accession]]=bestGene
    }
    

}
outFname=paste0(OutputDir,"step7.everyHMMsearch5MB/", HMMid,".RDS")
saveRDS(bestGenesList, outFname)
