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
options(scipen = 999)
args = commandArgs(trailingOnly=TRUE)
lowerRuleNumLimit=as.integer(args[1])
upperRuleNumLimit=as.integer(args[2])

AllPassedRules=NULL
tmpDIR=paste0(OutputDir, "/RulesCompliantInKeySpeciesMatrix/")
for (sp in c("hsapiens", relevantSpecies)){
    message(sp)
    fname=paste0(tmpDIR,sp,".",lowerRuleNumLimit,"-",upperRuleNumLimit,".txt")
    message(fname)
    tmpDF=read.table(fname, col.names="ruleID")
    tmpDF$specie=sp
    AllPassedRules=addToDF(AllPassedRules, tmpDF)
}
AllPassedRulesSum=AllPassedRules %>% group_by(ruleID) %>% summarize(Nspecies=length(unique(specie)), N=n()) %>% as.data.table()
saveRDS(AllPassedRulesSum, paste0(tmpDIR, "AllPassedRulesSum.",lowerRuleNumLimit,"-",upperRuleNumLimit,".RDS"))