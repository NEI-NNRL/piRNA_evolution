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
HMM=args[1]

message("HMMid: ", HMM)

simpleHitCheck=function(specie, HMM){
    #return values -1 no hit at all, 0 maybe, 1 positive
        SpecieDir=paste0(OutputDir,"/HMMsKeySpeciesResults/", specie, "/")
        HMMfname=paste0(SpecieDir,HMM,".faa.MSA.faa.hmm.out")
        HMMsearchDF=readHMMout(HMMfname)
        if (is.null(HMMsearchDF)){
            #essage("\tNo hits. Next genome please") 
            return(list(-1, HMMsearchDF))
        }
        HMMsearchDF=HMMsearchDF[ HMMsearchDF$score>minHMMscoreRatio*max(HMMsearchDF$score), ]
        if (nrow(HMMsearchDF)==0){
            #message("\tNo hits. Next genome please")
            return(list(-1, HMMsearchDF))
        }
        if (nrow(HMMsearchDF)==1){
            if (HMMsearchDF$targetname[1]==HMM){
                return(list(1, HMMsearchDF))
            }else{
                return(list(-1, HMMsearchDF)) # not even try to find it using the association rules
            }
        }
    return(list(0, HMMsearchDF))

}
requestKeyGTFbyHMM=function(inspecie, inHMM){
    queryTxt=sprintf("SELECT * FROM keySpeciesGenomicCoords WHERE 
                          specie = '%s' AND  type='CDS' AND (
                        specie_gene_id IN (SELECT DISTINCT homolog FROM HomologyBijectiveWellDefined WHERE ensembl_gene_id = '%s' AND specie = '%s')
                         OR specie_gene_id='%s') ", 
        inspecie, inHMM, inspecie, inHMM)
    GTFsDF1name=dbGetQuery(mydbRO, queryTxt)
    return(GTFsDF1name)
}

numcompliant_species=0
for (specie in c('hsapiens', relevantSpecies)){
    message(specie)
    tmplist=simpleHitCheck(specie, HMM)
    scheval=tmplist[[1]]
    HMMsearchDF=tmplist[[2]]
    if (scheval==-1) next
    if (scheval==1){
        numcompliant_species=numcompliant_species+1
        message("\t Simple compliance")
        next
    }
    # checking the rules
    ## removing the RHS
    HMMsearchDF=HMMsearchDF[HMMsearchDF$targetname!=HMM,]
    ## get potential neighbors
    queryTxt=sprintf('SELECT * FROM PassedAssociationRules5MBNamed WHERE RHS_ensembl_gene_id="%s"', HMM) # Already compliant
    HMMrulesDF=dbGetQuery(mydbRO, queryTxt)
    if (nrow(HMMrulesDF)==0){
        message("\tNo association rules found. Next")
        next
    }
    possibleLHSneightbours=unique(HMMrulesDF$LHS_ensembl_gene_id)
    # check gene IDs necessary for LHS
    geneIDvsNeighProximityDF=NULL
    for (targetHMMid in unique(HMMsearchDF$targetname)){
        tmpDF=data.frame(geneID=targetHMMid, neighEnsemblID=possibleLHSneightbours, foundnearby=FALSE)
        geneIDvsNeighProximityDF=addToDF(geneIDvsNeighProximityDF, tmpDF)
    }
    # list of candidate genes
    message("\tDEBUG list of candidate genes")
    CandidateGenesGRlist=list()
    for (targetHMMid in HMMsearchDF$targetname){
        # obtaining GRanges for the potential hits
        CandidateGenesGRlist[[targetHMMid]]=GRanges(requestKeyGTFbyHMM(specie, targetHMMid))
    }
            # genomic locations of the neighbours
    message("\tDEBUG NeighboursDF. Num possible neighs: ", length(possibleLHSneightbours))
    NeighboursDF=NULL
    specieChrNames=dbGetQuery(mydbRO, sprintf("SELECT DISTINCT( seqnames ) FROM keySpeciesGenomicCoords WHERE specie='%s'", specie))[[1]]
    for (neighHMMid in possibleLHSneightbours){
        message("\t\tneighHMMid: ",neighHMMid," ",getTimeStamp())
        SpecieDir=paste0(OutputDir,"/HMMsKeySpeciesResults/", specie, "/")
        HMMDFtmp=readHMMout(paste0(SpecieDir,neighHMMid,".faa.MSA.faa.hmm.out"))
        HMMDFtmp=HMMDFtmp[ HMMDFtmp$score>minHMMscoreRatio*max(HMMDFtmp$score), ]
        # get coordinates of the neighbour
        tmpNeighGTFdf=NULL
        for (targetname in HMMDFtmp$targetname){
            tmpNeighGTFdf=addToDF(tmpNeighGTFdf, requestKeyGTFbyHMM(specie, targetname))
        }
        if (is.null(tmpNeighGTFdf) || nrow(tmpNeighGTFdf)==0){ # neighbouring HMM not found
            next
        }
        tmpNeighGTFgr=GRanges(tmpNeighGTFdf)
        if (all(seqnames(tmpNeighGTFgr) %in% specieChrNames)){
            # preventing warnings
            seqlevels(tmpNeighGTFgr)=specieChrNames
        }
        # for each candidate gene check overlaps with each neightbour
        for (targetHMMid in HMMsearchDF$targetname){
            message("\t\t\t targetHMMid = ", targetHMMid)
            message("\t\t\t length(CandidateGenesGRlist[[targetHMMid]]) = ", length(CandidateGenesGRlist[[targetHMMid]]))
            if (length(CandidateGenesGRlist[[targetHMMid]])>0){
                distGeneNeight=getMinDistMultiChr(CandidateGenesGRlist[[targetHMMid]], tmpNeighGTFgr)
                geneIDvsNeighProximityDF$foundnearby[ geneIDvsNeighProximityDF$targetHMMid==targetHMMid & 
                                                 geneIDvsNeighProximityDF$neighEnsemblID == neighHMMid ] = (distGeneNeight <= lookaroundDist)
            }
        }
    }
    # presence of the neighbours near the candidate genes
    rulesDFbyGene=NULL
    message("\tDEBUG rulesDFbyGene")
    for (geneID in unique(geneIDvsNeighProximityDF$geneID)){
        subsetProvymityDF=geneIDvsNeighProximityDF[ geneIDvsNeighProximityDF$geneID == geneID ,]
        tmpDT=HMMrulesDF %>% group_by(ruleID) %>% 
            summarize(geneID=geneID, allLHSfound=all( LHS_ensembl_gene_id %in% 
                                      subsetProvymityDF$neighEnsemblID[ subsetProvymityDF$foundnearby ]
                                     ), numMembers=as.integer(mean(numMembers)) # do not understand why "sample" fails
                     ) %>% as.data.table()
        rulesDFbyGene=addToDF(rulesDFbyGene, tmpDT)
    }
    message("\tDEBUG genesRulesSummaryDF")
    suppressWarnings({ # -Inf warning for no neighbours found
        genesRulesSummaryDF=rulesDFbyGene %>% group_by(geneID) %>% 
            summarize(numCompliantRules=sum(allLHSfound), maxCompliantMembers=max(numMembers[allLHSfound]))
        })
    message("\tDEBUG selection")
    if (max(genesRulesSummaryDF$numCompliantRules)<10){
        message("\t Now enough false positive supporting rules")
        numcompliant_species=numcompliant_species+1
    }
    
}
message("Number of compliant species: ", numcompliant_species)
comToBash=sprintf("echo %d > %s/RulesCompliantInKeySpeciesMatrix/keySpeciesVerification/%s.txt", numcompliant_species, OutputDir, HMM)
#message(comToBash)
system(comToBash, intern = FALSE, wait = TRUE)
