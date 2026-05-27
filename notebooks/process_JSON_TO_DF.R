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
fname=args[1]
outdir=args[2]

library(tidyjson)

process_a_sample=function(inrecord){
    outlist=list()
    for (fn in names(inrecord)){
        if (typeof(inrecord[[fn]]) != 'list'){
            outlist[[fn]]=inrecord[[fn]]
        }
    }
    return(as.data.frame(outlist))
}

json_data <- read_json(fname)

DF_list=lapply(X = json_data$`..JSON`[[1]][['d']]$results, FUN = process_a_sample)
DF = bind_rows(DF_list)

outfname=paste0(outdir, '/',basename(fname),'.csv')
write.csv(DF, outfname)
