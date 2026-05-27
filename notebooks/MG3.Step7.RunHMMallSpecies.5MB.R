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

queryTxt="SELECT DISTINCT accession FROM GTFs"
accessions=dbGetQuery(mydbRO, queryTxt)[[1]]
bestGenesList=list()

for (accession in accessions){
    message(accession," ",getTimeStamp())
    SpecieDir=paste0(OutputDir,"HMMsAllSpeciesResults/", accession,"/")
    queryTxt=sprintf('SELECT seqnames FROM seqinfo WHERE accession="%s"', 
                     accession)
    accessionChrNames=dbGetQuery(mydbRO, queryTxt)[[1]]

    HMMout=paste0(SpecieDir,HMMid,".faa.MSA.faa.hmm.out")
    HMMsearchDF=readHMMout(HMMout)
    if (is.null(HMMsearchDF)){
        message("\tNo hits. Next genome please") 
        next
    }
    HMMsearchDF=HMMsearchDF[ HMMsearchDF$score>minHMMscoreRatio*max(HMMsearchDF$score), ]
    if (nrow(HMMsearchDF)==0){
        message("\tNo hits. Next genome please")
        next
    }
    geneIDs=NULL
    for (targetname in HMMsearchDF$targetname){
        tmpGTF=requestGTFbyProteinName(accession, targetname)
        geneName=unique(tmpGTF$gene)
        if (any(geneName=="?")){
            geneName=unique(tmpGTF$Name)
        }
        if (length(geneName)!=1){ 
            stop(length(geneName)," genes for ",accession, " " ,targetname)
        }
        geneIDs=unique(addToVec(geneIDs, geneName))
    }
    if (is.null(geneIDs)){ 
        message("\tNo gene names for ", accession, " " ,targetname)
        next
    }
    geneIDsDF=data.frame(geneID=geneIDs, ensembl_gene_id=NA, exact_hit=FALSE)
    for (i in 1:nrow(geneIDsDF)){
        geneName=geneIDsDF$geneID[i]
        queryTxt=sprintf('SELECT DISTINCT ensembl_gene_id FROM BijectiveRelevantSpeciesView WHERE
                         homolog IN (SELECT DISTINCT specieEnsemblID FROM GeneNameSynonims WHERE synonym="%s")', geneName) 
        nonHumanEnseblIDs=dbGetQuery(mydbRO, queryTxt)[[1]]
        queryTxt=sprintf('SELECT DISTINCT specieEnsemblID FROM GeneNameSynonims WHERE synonym="%s" 
                         AND specie="hsapiens" AND ensembl_gene_id IN (SELECT DISTINCT ensembl_gene_id FROM BijectiveRelevantSpeciesView )', 
                         geneName) 
        HumanEnseblIDs=dbGetQuery(mydbRO, queryTxt)[[1]]
        EnseblIDs=unique(c(nonHumanEnseblIDs, HumanEnseblIDs))

        for (EnseblID in EnseblIDs){
           geneIDsDF$ensembl_gene_id[i]=EnseblIDs # check for NA returns
            if (geneIDsDF$ensembl_gene_id[i] == HMMid){
                geneIDsDF$exact_hit[i]=TRUE
                #message("\tHMM ID exact correspondence")
            } 
        }
    }
    if (sum(geneIDsDF$exact_hit)==1 & length(unique(geneIDsDF$ensembl_gene_id))==1){
        message("\tFound. No further search")
        BestGene=unique(geneIDsDF$geneID[ geneIDsDF$exact_hit])
        message("\tBest gene: ", BestGene)
        bestGenesList[[accession]]=BestGene
    }else{
        message("\tNeed more analysis")

        
        # association rules and possible neighbours
        queryTxt=sprintf('SELECT * FROM PassedAssociationRules5MBNamed WHERE RHS_ensembl_gene_id="%s"', HMMid) # Already compliant
        HMMrulesDF=dbGetQuery(mydbRO, queryTxt)
        message("\tDEBUG HMMrulesDF")
        if (nrow(HMMrulesDF)==0){
            message("No association rules found. Next")
            next
        }
        possibleLHSneightbours=unique(HMMrulesDF$LHS_ensembl_gene_id)
        geneIDvsNeighProximityDF=NULL
        for (geneID in geneIDsDF$geneID){
            tmpDF=data.frame(geneID=geneID, neighEnsemblID=possibleLHSneightbours, foundnearby=FALSE)
            geneIDvsNeighProximityDF=addToDF(geneIDvsNeighProximityDF, tmpDF)
        }
        # list of candidate genes
        message("\tDEBUG list of candidate genes")
        CandidateGenesGRlist=list()
        for (geneID in geneIDsDF$geneID){
            CandidateGenesGRlist[[geneID]]=GRanges(requestGTFbyProteinName(accession, geneID))
        }
        # genomic locations of the neighbours
        message("\tDEBUG NeighboursDF. Num possible neighs: ", length(possibleLHSneightbours))
        NeighboursDF=NULL
        for (neighHMMid in possibleLHSneightbours){
            message("\t\tneighHMMid: ",neighHMMid," ",getTimeStamp())
            HMMDFtmp=readHMMout(paste0(SpecieDir,neighHMMid,".faa.MSA.faa.hmm.out"))
            HMMDFtmp=HMMDFtmp[ HMMDFtmp$score>minHMMscoreRatio*max(HMMDFtmp$score), ]
            # get coordinates of the neighbour
            tmpNeighGTFdf=NULL
            for (targetname in HMMDFtmp$targetname){
                tmpNeighGTFdf=addToDF(tmpNeighGTFdf, requestGTFbyProteinName(accession, targetname))
            }
            if (is.null(tmpNeighGTFdf)){ # neighbouring HMM not found
                next
            }
            tmpNeighGTFgr=GRanges(tmpNeighGTFdf)
            if (all(seqnames(tmpNeighGTFgr) %in% accessionChrNames)){
                # preventing warnings
                seqlevels(tmpNeighGTFgr)=accessionChrNames
            }
            # for each candidate gene check overlaps with each neightbour
            for (geneID in geneIDsDF$geneID){
                distGeneNeight=getMinDistMultiChr(CandidateGenesGRlist[[geneID]], tmpNeighGTFgr)
                geneIDvsNeighProximityDF$foundnearby[ geneIDvsNeighProximityDF$geneID==geneID & 
                                                     geneIDvsNeighProximityDF$neighEnsemblID == neighHMMid ] = (distGeneNeight <= lookaroundDist)
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
            message("\tNo compliant genes")
            next
        }else{
            message("\tDEBUG BestGeneBoolMask")
            numRulesBoolMask= ( genesRulesSummaryDF$numCompliantRules>0.1*max(genesRulesSummaryDF$numCompliantRules) )
            numMembersBoolMask= ( genesRulesSummaryDF$maxCompliantMembers==max(genesRulesSummaryDF$maxCompliantMembers) )
            BestGeneBoolMask=( numRulesBoolMask & numMembersBoolMask)
            if (sum(BestGeneBoolMask)!=1){
                message("\tNo clear best gene")
                next
            }else{
                BestGene=genesRulesSummaryDF$geneID[BestGeneBoolMask]
                message("\tBest gene: ", BestGene)
                bestGenesList[[accession]]=BestGene
            }
        }
    }
}
outFname=paste0(OutputDir,"step7.everyHMMsearch5MB/", HMMid,".RDS")
saveRDS(bestGenesList, outFname)
