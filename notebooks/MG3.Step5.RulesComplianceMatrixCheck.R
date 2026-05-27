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
allTransactionsDFspecieFname=args[1]
apriori_rulesFname=args[2]
lowerRuleNumLimit=as.integer(args[3])
upperRuleNumLimit=as.integer(args[4])
specie=args[5]

library(Matrix)
library(mltools)
library(sparseMatrixStats)

apriori_rules=readRDS(apriori_rulesFname)
allTransactionsDFspecie=readRDS(allTransactionsDFspecieFname)


LHSrulesMatr=as(apriori_rules@lhs@data, "TsparseMatrix")
RHSrulesMatr=as(apriori_rules@rhs@data, "TsparseMatrix")
LHSrulesMatrInt <- spMatrix(nrow = LHSrulesMatr@Dim[2], ncol = LHSrulesMatr@Dim[1], j=LHSrulesMatr@i +1, i=LHSrulesMatr@j +1, x=rep(1,length(LHSrulesMatr@i)))
RHSrulesMatrInt <- spMatrix(nrow = RHSrulesMatr@Dim[2], ncol = RHSrulesMatr@Dim[1], j=RHSrulesMatr@i +1, i=RHSrulesMatr@j +1, x=rep(1,length(RHSrulesMatr@i)))
RulesMatrFull=LHSrulesMatrInt+RHSrulesMatrInt
message(paste("RulesMatrFull dim", paste0(RulesMatrFull@Dim, collapse=", ")))

tmpTransactionsDFspecie=allTransactionsDFspecie[ allTransactionsDFspecie$specie==specie, ]
tmpTransactionsDFspecie$geneFactor=as.integer(factor(tmpTransactionsDFspecie$ensembl_gene_id, 
                                                         levels = apriori_rules@lhs@itemInfo$labels))
transactionIDfactors=factor(tmpTransactionsDFspecie$transactionID)
TmatrixSpecie=spMatrix(nrow = length(apriori_rules@lhs@itemInfo$labels), ncol = max(as.integer(transactionIDfactors)), 
                      i=tmpTransactionsDFspecie$geneFactor, j=as.integer(transactionIDfactors),
                       x=rep(1, length(tmpTransactionsDFspecie$geneFactor))
                      )
colnames(TmatrixSpecie)=levels(transactionIDfactors)

message(paste("TmatrixSpecie dim", paste0(TmatrixSpecie@Dim, collapse=", ")))


RulesMatr=RulesMatrFull[lowerRuleNumLimit:upperRuleNumLimit,]
rownames(RulesMatr)<-lowerRuleNumLimit:upperRuleNumLimit

RulesNmembers=rowSums(RulesMatr)

proRulesTransactionsMatrxSpecie=RulesMatr %*% TmatrixSpecie
maxPerRule=rowMaxs(proRulesTransactionsMatrxSpecie)
GoodRulesIDsTMP=which(maxPerRule>=RulesNmembers)
speciePassedRulesIDs= names(GoodRulesIDsTMP)

writeLines(speciePassedRulesIDs, paste0(dirname(allTransactionsDFspecieFname),"/",specie,".",lowerRuleNumLimit,"-",upperRuleNumLimit,".txt"))