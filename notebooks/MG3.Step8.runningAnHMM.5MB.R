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
outFname=paste0(OutputDir,"step8.HMMgenomeDescriptions/",HMMid,".RDS")
if (file.exists(outFname)){
    message("\t Found result. Exiting...")
    quit(status = 0)
}

checkGeneNameIntegrity=function(inlist){
    # get genomic coordinates of the gene in genome.
    # check if they make sense
    coordsDFtmp=requestGTFbyProteinName(inlist[2], inlist[3])
    if (length(unique(coordsDFtmp$seqnames))!=1) return(FALSE)
    if (length(unique(coordsDFtmp$strand))!=1) return(FALSE)
    coordsDFtmp$end=as.integer(coordsDFtmp$end) # max value is 2e9 should be enough for any mammalian chromosome
    coordsDFtmp$start=as.integer(coordsDFtmp$start)
    # check if all exons are within 5 MB range. Longest human gene is 2.4 MB (dystrophin)
    if (max(coordsDFtmp$end)-min(coordsDFtmp$start) > 5e6) return(FALSE)
    return(TRUE)
}

queryTxt=sprintf('SELECT * FROM AllSpeciesFoundGeneNames5MB WHERE accession IN (SELECT accession FROM EutherianGenomes) AND ensembl_gene_id="%s"', 
                     HMMid)
accessionsDF=dbGetQuery(mydbRO, queryTxt)
accessionsDF$singleGene=apply(X = accessionsDF, MARGIN = 1, FUN = checkGeneNameIntegrity)
#filter 1
accessionsDFf1=accessionsDF[ accessionsDF$singleGene, ]
if ( nrow( accessionsDFf1 ) < minDistinctGenomesForStats){
    message("Not enough genomes")
    quit(status = 0)
}
# adding suborder info
queryTxt='SELECT name AS familyName, id AS familyID, accession  FROM taxonomy WHERE rank="family"'
taxonomyDF=dbGetQuery(mydbRO, queryTxt)
accessionsDFf1=merge(accessionsDFf1, taxonomyDF, by="accession", all.x = TRUE, all.y = FALSE)
familiesDF=accessionsDFf1 %>% group_by(familyID) %>% summarize(Ngemomes=n())
numFamiliesForStas=nrow(familiesDF[ familiesDF$Ngemomes>=minGenomesInFamiliesForStats,])
if ( numFamiliesForStas < minDistinctFamiliesForStats){
    message("Not enough families with at least ", minGenomesInFamiliesForStats, " Genomes")
    quit(status = 0)
}

#going over the genomes
accessionsDFf1$passedIndividualAnalysis=FALSE
accessionsDFf1$DSlen=NA
accessionsDFf1$GapLen=NA
accessionsDFf1$sRMSKlen=NA
accessionsDFf1$asRMSKlen=NA
accessionsDFf1$sRMSKlen10K=NA
accessionsDFf1$asRMSKlen10K=NA
accessionsDFf1$sRMSKlen50K=NA
accessionsDFf1$asRMSKlen50K=NA
accessionsDFf1$sRMSKlen100K=NA
accessionsDFf1$asRMSKlen100K=NA
message("going over the genomes")
for (i in 1:nrow(accessionsDFf1)){
    accession=accessionsDFf1$accession[i]
    message("\t", accession)
    geneName=accessionsDFf1$geneName[i]
    tmpHMMDF=requestGTFbyProteinName(inaccession = accession, inname = geneName)
    lastCDS=getLastCDS(tmpHMMDF)
    nextDS.CDS=requestDoqnstreamCDSSQL(lastCDS, 'GTFs', 'suspectedTEgenes')
    DSintervalGR=getCDSdowstream(lastCDS, nextDS.CDS)
    if (is.null(DSintervalGR)) next
    # get sequence of the SD interval
    DSseq=getSeqIndexedFA(getFApath(accession), as.vector(seqnames(DSintervalGR)), 
                          as.vector(start(DSintervalGR)), as.vector(width(DSintervalGR)))
    if (is.null(DSseq)) next
    gapLen=str_count(DSseq, "[^ATCGatcg]")
    if (gapLen/str_length(DSseq)>maxGapFraction){
        message("GAP TOO LARGE ",lastCDS$accession )
        next
    }
    # checking RMSK track
    chrLen=getChrLenSQL(accession, lastCDS$seqnames)
    if (is.na(chrLen)) next
    queryTxt=sprintf('SELECT * FROM RMSK WHERE accession="%s" AND seqnames="%s"',
                             accession, lastCDS$seqnames)
    RMSKdf=dbGetQuery(mydbRO, queryTxt)
    RMSKgr=GRanges(RMSKdf)
    RMSKred=GenomicRanges::reduce(RMSKgr)
    RMSKchrFraction=sum(width(RMSKred))/chrLen
    if (RMSKchrFraction<minChrRMSKfraction){
        message("No confidence in RMSK", accession)
        next
    }
    
    accessionsDFf1$passedIndividualAnalysis[i]=TRUE
    # passed individual tests
    # adding stats
    accessionsDFf1$DSlen[i]=str_length(DSseq)
    accessionsDFf1$GapLen[i]=gapLen

    
    DSRMSK=getDSranges(lastCDS, str_length(DSseq), "RMSK")
    DSRMSK.RED=GenomicRanges::reduce(DSRMSK)
    sRMSK.RED=GenomicRanges::intersect(DSRMSK.RED, DSintervalGR)
    asRMSK.RED=GenomicRanges::intersect(DSRMSK.RED, GenomicRanges::invertStrand(DSintervalGR))
    accessionsDFf1$sRMSKlen[i]=sum(width(sRMSK.RED))
    accessionsDFf1$asRMSKlen[i]=sum(width(asRMSK.RED))

    availDist=distToChrEnd(lastCDS)
    if (availDist < 10e3) next
    DSintervalToAnalyze=getDSinterval(lastCDS, 10e3)
    DSRMSK=getDSranges(lastCDS, 10e3, "RMSK")
    DSRMSK.RED=GenomicRanges::reduce(DSRMSK)
    sRMSK.RED=GenomicRanges::intersect(DSRMSK.RED, DSintervalToAnalyze)
    asRMSK.RED=GenomicRanges::intersect(DSRMSK.RED, GenomicRanges::invertStrand(DSintervalToAnalyze))
    accessionsDFf1$sRMSKlen10K[i]=sum(width(sRMSK.RED))
    accessionsDFf1$asRMSKlen10K[i]=sum(width(asRMSK.RED))

    if (availDist < 50e3) next
    DSintervalToAnalyze=getDSinterval(lastCDS, 50e3)
    DSRMSK=getDSranges(lastCDS, 50e3, "RMSK")
    DSRMSK.RED=GenomicRanges::reduce(DSRMSK)
    sRMSK.RED=GenomicRanges::intersect(DSRMSK.RED, DSintervalToAnalyze)
    asRMSK.RED=GenomicRanges::intersect(DSRMSK.RED, GenomicRanges::invertStrand(DSintervalToAnalyze))
    accessionsDFf1$sRMSKlen50K[i]=sum(width(sRMSK.RED))
    accessionsDFf1$asRMSKlen50K[i]=sum(width(asRMSK.RED))

    if (availDist < 1e5) next
    DSintervalToAnalyze=getDSinterval(lastCDS, 1e5)
    DSRMSK=getDSranges(lastCDS, 1e5, "RMSK")
    DSRMSK.RED=GenomicRanges::reduce(DSRMSK)
    sRMSK.RED=GenomicRanges::intersect(DSRMSK.RED, DSintervalToAnalyze)
    asRMSK.RED=GenomicRanges::intersect(DSRMSK.RED, GenomicRanges::invertStrand(DSintervalToAnalyze))
    accessionsDFf1$sRMSKlen100K[i]=sum(width(sRMSK.RED))
    accessionsDFf1$asRMSKlen100K[i]=sum(width(asRMSK.RED))
}
accessionsDFf2 = accessionsDFf1[ accessionsDFf1$passedIndividualAnalysis, ]

familiesDF=accessionsDFf2 %>% group_by(familyID) %>% summarize(Ngemomes=n())
numFamiliesForStas=nrow(familiesDF[ familiesDF$Ngemomes>=minGenomesInFamiliesForStats,])
if ( numFamiliesForStas < minDistinctFamiliesForStats){
    message("Not enough families with at least ", minGenomesInFamiliesForStats, " Genomes. (filter 2)")
    quit(status = 0)
}


saveRDS(accessionsDFf2, outFname)
