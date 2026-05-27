# general functions for the MG3 attempt

suppressMessages({
    library("rjson")
    library(stringr)
    library(data.table)
    library(dplyr)
    library(rtracklayer)
    library(GenomicRanges)
    library(parallel)
    library(pbmcapply)
    library(DBI)
    library(biomaRt)
    library(taxize)
    library(ggplot2)
    library(Biostrings)
})

getTimeStamp<-function(){
    timeMark<-Sys.time()
    timeMark=gsub(pattern = " ", replacement = "_", x = timeMark)
    timeMark=gsub(pattern = ":", replacement = ".", x = timeMark)
    return(timeMark)
}
addToVec=function(inVec, inVal){
    if (is.null(inVec)){
        return(inVal)
    }
    return(c(inVec, inVal))
}
addToDF=function(inDF, inVal){
    if (is.null(inDF)){
        return(inVal)
    }
    extracols=setdiff(colnames(inDF), colnames(inVal))
    for (cn in extracols){
        inVal[[cn]]=NA
    }
    extracols=setdiff(colnames(inVal), colnames(inDF))
    for (cn in extracols){
        inDF[[cn]]=NA
    }
    inDF=rbind(inDF, as.data.frame(inVal))
    return(inDF)
}
options(repr.matrix.max.cols=150, repr.matrix.max.rows=150)
apikey=stop("Set apikey variable value for ENTREZ access in the MG3.functions1.R") # remove the stop() function and just set the value
options(ENTREZ_KEY=apikey)
OutputDir<-"~/data/lambda/output/"
StaticInpDir<-"~/data/lambda/data/"
SQLfname="lambda_project_reproduction.sqlite"
EntrezSQLfname="entrez_lambda_project_reproduction.sqlite"
EntrezTSVfnames=c("Compara.115.protein_default.homologies.tsv", "Compara.112.protein_default.homologies.tsv")
LogDir<-paste0(OutputDir,"slurm_logs")

minCHrLim=1500000
minBEDfileSize=5000
maxPvalue="1E-10" # for HMMsearch
maxDelta=1000000
lookaroundDist=5e6 # for transaction
maxGapFraction=0.1
maxDSanalysisLen=maxDelta*10
minChrRMSKfraction=0.20
minDistToCHRend=100000
minDistinctFamiliesForStats=10
minGenomesInFamiliesForStats=3
minDistinctGenomesForStats=50
minHMMscoreRatioKeySpecies=0.75
minHMMscoreRatio=0.85
maxHMMfalsePositiveHits=2 # number of times an HMM can predict wrong protein in key specie before I find it useless
minRMSKinsertionlen=10 # discard TE insertion shorter than
Ensembl_Version=112


# KEY species
KeyCommonName=c("Human", "Mouse", "Dog", "Cat", "Cow", "Camel", "AfricanSavannaElephant", "Pig", "BlueWhale")
TaxIDs=c(9606, 10090, 9615, 9685, 9913, 9838, 9785, 9823, 9771)
ScientificName=c("homo_sapiens", "mus_musculus", "canis_lupus_familiaris" , "felis_catus", "bos_taurus", "camelus_dromedarius","loxodonta_africana", "sus_scrofa", "balaenoptera_musculus")

EntrezHomolog=c("", "mmusculus", "clfamiliaris" , "fcatus", "btaurus", "cdromedarius","lafricana", "sscrofa", "bmusculus")
GenomeAccession=c("GCF_009914755.1", "GCF_000001635.27", "GCF_011100685.1", "GCF_018350175.1", 
"GCF_002263795.3", "GCF_036321535.1", "GCF_030014295.1", "GCF_000003025.6" ,"GCF_009873245.2")
EnsemblGTF=c("Homo_sapiens.GRCh38.115.gtf", "Mus_musculus.GRCm39.115.gtf", "Canis_lupus_familiaris.ROS_Cfam_1.0.115.gtf", 
             "Felis_catus.F.catus_Fca126_mat1.0.115.gtf", 
            "Bos_taurus.ARS-UCD1.3.112.gtf", 
             "Camelus_dromedarius.CamDro2.115.gtf", "Loxodonta_africana.loxAfr3.115.gtf",
             "Sus_scrofa.Sscrofa11.1.115.gtf" ,"Balaenoptera_musculus.mBalMus1.v2.115.gtf")
EnsemblePepFA=c("Homo_sapiens.GRCh38.pep.all.fa", "Mus_musculus.GRCm39.pep.all.fa", "Canis_lupus_familiaris.ROS_Cfam_1.0.pep.all.fa",
               "Felis_catus.F.catus_Fca126_mat1.0.pep.all.fa", "Bos_taurus.ARS-UCD1.3.pep.all.fa", "Camelus_dromedarius.CamDro2.pep.all.fa",
               "Loxodonta_africana.loxAfr3.pep.all.fa", "Sus_scrofa.Sscrofa11.1.pep.all.fa", "Balaenoptera_musculus.mBalMus1.v2.pep.all.fa")
EnsembleFTPpath=c("http://ftp.ensembl.org/pub/current_mysql/homo_sapiens_core_115_38/",
                  "http://ftp.ensembl.org/pub/current_mysql/mus_musculus_core_115_39/",
                  "http://ftp.ensembl.org/pub/current_mysql/canis_lupus_familiaris_core_115_1/",
                  "http://ftp.ensembl.org/pub/current_mysql/felis_catus_core_115_1/",
                  "https://ftp.ensembl.org/pub/release-112/mysql/bos_taurus_core_112_13/",
                  "https://ftp.ensembl.org/pub/current_mysql/camelus_dromedarius_core_115_2/",
                  "https://ftp.ensembl.org/pub/current_mysql/loxodonta_africana_core_115_3/",
                  "https://ftp.ensembl.org/pub/current_mysql/sus_scrofa_core_115_111/",
                  "https://ftp.ensembl.org/pub/current_mysql/balaenoptera_musculus_core_115_2/"
                 )
KeySpDF=data.frame(common=KeyCommonName, taxon=TaxIDs, HomologueTable=EntrezHomolog, genome=GenomeAccession, EnsemblGTF=EnsemblGTF, EnsemblePepFA=EnsemblePepFA,
                  ScientificName=ScientificName, EnsembleFTPpath=EnsembleFTPpath)
relevantSpecies=KeySpDF$HomologueTable
relevantSpecies=relevantSpecies[ nchar(relevantSpecies) > 0]

#cleaning SQL connections
objs=ls()
for (ob in objs){
    if (class(ob)=="SQLiteConnection"){
        dbDisconnect(mydb)
    }
}
if (file.exists(paste0(OutputDir,EntrezSQLfname))){
    EntrezdbRW <- dbConnect(RSQLite::SQLite(), dbname=paste0(OutputDir,EntrezSQLfname), vfs="unix-none")
    EntrezdbRO <- dbConnect(RSQLite::SQLite(), dbname=paste0(OutputDir,EntrezSQLfname), flags=RSQLite::SQLITE_RO, vfs="unix-none")
}
if (file.exists(paste0(OutputDir,SQLfname))){
    mydbRW <- dbConnect(RSQLite::SQLite(), dbname=paste0(OutputDir,SQLfname), vfs="unix-none")
    mydbRO <- dbConnect(RSQLite::SQLite(), dbname=paste0(OutputDir,SQLfname), flags=RSQLite::SQLITE_RO, vfs="unix-none")
    
    #creating SQL views
    if ("homology" %in% dbListTables(mydbRO)){
        queryTxt=sprintf("CREATE TEMP VIEW relevantSpeciesView AS SELECT * FROM homology WHERE specie IN ('%s')", paste0(relevantSpecies, collapse = "', '"))
        dbExecute(mydbRO, queryTxt)
        queryTxt=sprintf("CREATE TEMP VIEW BijectiveRelevantSpeciesView AS 
        WITH bad_homolog_ids AS 
            (SELECT DISTINCT homolog FROM relevantSpeciesView GROUP BY homolog HAVING COUNT(DISTINCT ensembl_gene_id) !=1),
    bad_ensemble_ids AS 
            (SELECT DISTINCT ensembl_gene_id FROM relevantSpeciesView WHERE homolog IN (SELECT homolog FROM bad_homolog_ids)),
    good_ensemble_ids AS
            ( SELECT ensembl_gene_id FROM 
                ( SELECT ensembl_gene_id, COUNT( DISTINCT specie ) AS numspecies, COUNT() AS numrecords 
                    FROM relevantSpeciesView GROUP BY ensembl_gene_id) 
        	  WHERE numspecies=%d AND numrecords=%d )
    SELECT * FROM relevantSpeciesView WHERE ensembl_gene_id IN (SELECT ensembl_gene_id FROM good_ensemble_ids) AND
        ensembl_gene_id NOT IN (SELECT ensembl_gene_id FROM bad_ensemble_ids) AND
        homolog NOT IN (SELECT homolog FROM bad_homolog_ids)", 
        length(relevantSpecies), length(relevantSpecies)) 
        dbExecute(mydbRO, queryTxt)
    }
    if ("HMMkeySpecieCrossHits" %in% dbListTables(mydbRO)){
        dbExecute(mydbRO,"CREATE TEMP VIEW IF NOT EXISTS HMMkeySpecieCrossHitsCount AS SELECT HMM, specie, COUNT() AS C FROM HMMkeySpecieCrossHits GROUP BY HMM, specie ")
        #queryTxt=sprintf('CREATE TEMP VIEW IF NOT EXISTS HMMkeySpecieCrossHitsWellDefined AS SELECT * FROM HMMkeySpecieCrossHits WHERE HMM NOT IN (SELECT HMM FROM HMMkeySpecieCrossHitsCount WHERE C > %d )', maxHMMfalsePositiveHits)
        #dbExecute(mydbRO, queryTxt)
        queryTxt=sprintf('CREATE TEMP VIEW IF NOT EXISTS HomologyBijectiveWellDefined AS SELECT * FROM BijectiveRelevantSpeciesView WHERE ensembl_gene_id NOT IN (SELECT HMM FROM HMMkeySpecieCrossHitsCount WHERE C > %d )', maxHMMfalsePositiveHits)
        dbExecute(mydbRO, queryTxt)
    }
    if ("PassedAssociationRules5MB" %in% dbListTables(mydbRO)){
        dbExecute(mydbRO,"CREATE TEMP VIEW IF NOT EXISTS PassedAssociationRules5MBNamed AS SELECT dict1.ensembl_gene_id AS RHS_ensembl_gene_id, dict2.ensembl_gene_id AS LHS_ensembl_gene_id, numMembers, ruleID FROM PassedAssociationRules5MB INNER JOIN PassedAssociationGeneFactors5MB dict1 ON PassedAssociationRules5MB.RHSgeneFactor = dict1.geneFactor INNER JOIN PassedAssociationGeneFactors5MB dict2 ON PassedAssociationRules5MB.LHSgeneFactor = dict2.geneFactor")
    }
    
    #functions
    
    renameDFcolumn=function(inDF, incol, innewname){
        oldcol=colnames(inDF)
        if (length(oldcol)==0) stop("No columns in DF")
        boolmask=grepl(paste0("^",incol,"$"), oldcol)
        if (sum(as.integer(boolmask))!=1) stop("Not singular column to rename (0 or more than 1)")
        oldcol[boolmask]=innewname
        colnames(inDF)=oldcol
        return(inDF)
    }
    
    getGenomeChrs=function(inAccession){
        inFiname=paste0(StaticInpDir, inAccession, "/sequence_report.jsonl")
        inLines=readLines(inFiname)
        chrlens=as.integer(str_extract(string = inLines, pattern = ',\\s*"length"\\s*:\\s*(\\d+)\\s*,',group = 1 ))
        boolmask=(chrlens>=minCHrLim)
        boolmask[is.na(boolmask)]=FALSE
        goodLines=inLines[boolmask]
        thetext=paste0(goodLines, collapse = ",\n")
    
        thetext=paste0('{"assemblies": [',thetext,"]}")
        result <- fromJSON(thetext)
        outDF=NULL
        for (assembly in result$assemblies){
            if (is.null(outDF)){
                outDF=as.data.frame(assembly)
            }else{
                extracols=setdiff(colnames(outDF), names(assembly))
                for (cn in extracols){
                    assembly[[cn]]=NA
                }
                extracols=setdiff(names(assembly), colnames(outDF))
                for (cn in extracols){
                    outDF[[cn]]=NA
                }
                outDF=rbind(outDF, as.data.frame(assembly))
            }
        }
        return(outDF)
    }
    
    readSpecieGTF=function(inGTFname){
        TMP.GTF=import.gff(inGTFname)
        TMP.DT=as.data.table(TMP.GTF)
        TMP.DT$seqnames=as.character(TMP.DT$seqnames)
        TMP.DT$chromosome=NULL # to not cofuse the GRanges
        #checking if all the fields are present
        keyfields=c("Name", "protein_id", "gene")
        if (! all(keyfields %in% colnames(TMP.DT))){
            return(NULL)
        }
        return(TMP.DT)
    }
    
    retrieveSpecieProteinsPlain=function(inSp, inGenes){
        message(inSp)
        tmp_ensembl = useEnsembl(biomart="ensembl", dataset=paste0(inSp,"_gene_ensembl"), version = Ensembl_Version)
        tmp_transcriptsDF=getBM(mart = tmp_ensembl, 
            attributes=c('ensembl_gene_id', 'ensembl_peptide_id', 'ensembl_transcript_id', 'transcript_is_canonical'),
            filters    = "ensembl_gene_id",
            values     = inGenes
                        )
        message("NA filter")
        tmp_transcriptsDF=tmp_transcriptsDF[!is.na(tmp_transcriptsDF$transcript_is_canonical),]
        message("Sleep")
        Sys.sleep(10)
        message("Peptides query")
        tmp_peptidesDF=as.data.table(getBM(mart = tmp_ensembl, 
            attributes=c('ensembl_peptide_id', 'peptide', 'peptide_version', 'cds_length'),
            filters    = "ensembl_peptide_id",
            values     = tmp_transcriptsDF$ensembl_peptide_id)
                                     )
        message("Key 1")
        setkey(tmp_peptidesDF, ensembl_peptide_id)
        message("Version filter")
        tmp_peptidesDF = tmp_peptidesDF %>% group_by(ensembl_peptide_id) %>% filter(peptide_version == max(peptide_version))
        message("Length filter")
        tmp_peptidesDF = as.data.table(tmp_peptidesDF %>% group_by(ensembl_peptide_id) %>% filter(cds_length == max(cds_length)))
        message("Key 2")
        setkey(tmp_peptidesDF, ensembl_peptide_id)
        tmp_transcriptsDF=as.data.table(tmp_transcriptsDF)
        message("Key 3")
        setkey(tmp_transcriptsDF, ensembl_peptide_id)
        message("Merge in function")
        tmpDF=merge(tmp_transcriptsDF, tmp_peptidesDF, by = "ensembl_peptide_id", all =FALSE)
        tmpDF$specie=inSp
        return(tmpDF)
    }
    retrieveSpecieProteinsLocalSQL=function(inSp, inGenes){
        message(inSp)
        queryTxt=sprintf('SELECT DISTINCT gene_id AS ensembl_gene_id, gene_version, transcript_id AS ensembl_transcript_id, transcript_version FROM gene_transcripts WHERE specie="%s" AND tag="Ensembl_canonical"', inSp)
        tmp_transcriptsDF=dbGetQuery(EntrezdbRO, queryTxt)
        tmp_transcriptsDF=tmp_transcriptsDF %>% group_by(ensembl_gene_id) %>% filter(gene_version==max(gene_version)) %>% ungroup() %>% 
                group_by(ensembl_transcript_id) %>% filter(transcript_version==max(transcript_version)) %>% ungroup()
        message("Peptides query")
        queryTxt=sprintf('SELECT DISTINCT peptide_id AS ensembl_peptide_id, peptide_ver AS peptide_version, transcript_id AS ensembl_transcript_id, transcript_ver, peptide, LENGTH(peptide) AS peptide_length FROM peptides WHERE specie="%s"', inSp)
        tmp_peptidesDF=dbGetQuery(EntrezdbRO, queryTxt)
        tmp_peptidesDF = tmp_peptidesDF %>% group_by(ensembl_peptide_id) %>% filter(peptide_version == max(peptide_version))
        message("Length filter")
        tmp_peptidesDF = as.data.table(tmp_peptidesDF %>% group_by(ensembl_peptide_id) %>% filter(peptide_length == max(peptide_length)))
        message("Key 2")
        setkey(tmp_peptidesDF, ensembl_transcript_id)
        tmp_transcriptsDF=as.data.table(tmp_transcriptsDF)
        message("Key 3")
        setkey(tmp_transcriptsDF, ensembl_transcript_id)
        message("Merge in function")
        tmpDF=merge(tmp_transcriptsDF, tmp_peptidesDF, by = "ensembl_transcript_id", all =FALSE)
        tmpDF$specie=inSp
        tmpDF=tmpDF[ tmpDF$ensembl_gene_id %in% inGenes, ]
        return(tmpDF)
    }

    
    readHMMout=function(inFname){
        HMMlines=readLines(inFname)
        HMMlines=HMMlines[! grepl(pattern = "^#.*", x = HMMlines)]
        if (length(HMMlines)==0){
            return(NULL)
        }
        HMM.DT=as.data.table(t(as.data.table(str_split(string = HMMlines, pattern = " +", n=19))))
        colnames(HMM.DT)=c("targetname", "accession", "queryname", "accession2", "Evalue", "score", "bias", "Evalue.DOM", "score.DOM", "bias.DOM", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "descriptionoftarget")
        rownames(HMM.DT)=NULL
        HMM.DT=HMM.DT[, str_count(string = colnames(HMM.DT), pattern = "DOM")==0, with=FALSE]
        HMM.DT=HMM.DT[, str_count(string = colnames(HMM.DT), pattern = "accession")==0, with=FALSE]
        HMM.DT=as.data.table(HMM.DT)
        for (cn in c("Evalue", "score", "bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc")){
            HMM.DT[[cn]]=as.numeric(HMM.DT[[cn]])
        }
        return(HMM.DT)
    }
    
    #Minimal distance between 2 set of CDSs
    getSpan<-function(inGR){
        chrName=unique(as.vector(seqnames(inGR)))
        if (length(chrName)!=1) stop("Problem with chr names")
        strandName=unique(as.vector(strand(inGR)))
        if (length(strandName)!=1) stop("Problem with strands")
        startCoord=min(as.vector(start(inGR)))
        endCoord=max(as.vector(end(inGR)))
        retGR=GRanges(seqnames = chrName, ranges = IRanges(start = startCoord, end = endCoord), 
                      strand = strandName)
        seqlevels(retGR)=seqlevels(inGR)
        return(retGR)
    }
    
    getSpansManyChrs<-function(inGR){
        chrNames=unique(as.vector(seqnames(inGR)))
        strandNames=unique(as.vector(strand(inGR)))
        outGR=NULL
        for (chrName in chrNames){
            for (strandName in strandNames){
                if (any( seqnames(inGR) == chrName & strand(inGR) == strandName)){
                    tmpGR=getSpan( inGR[ seqnames(inGR) == chrName & strand(inGR) == strandName ])
                    outGR=addToVec(outGR, tmpGR)
                }
            }
        }
        return(outGR)
    }
    
    getMinDist<-function(inGR1, inGR2, checkCDS=TRUE){
        if (checkCDS){
            inGR1<-inGR1[inGR1$type=="CDS"]
            inGR2<-inGR2[inGR2$type=="CDS"]
        }
        span1=getSpan(inGR1)
        span2=getSpan(inGR2)
        #message(span1)
        #message(span2)
        if (length( GenomicRanges::intersect(span1, span2))>0) return(0)
        if (!any(seqnames(span1)==seqnames(span2))) return(Inf) # different chromosomes
        return(min(c(abs(start(span2)-end(span1)),abs(start(span1)-end(span2)))))
        
    }
    
    getMinDistMultiChr<-function(inGR1, inGR2){
        inGR1<-inGR1[inGR1$type=="CDS"]
        inGR2<-inGR2[inGR2$type=="CDS"]
        spans1=getSpansManyChrs(inGR1)
        spans2=getSpansManyChrs(inGR2)
        minDist=Inf
        for (i in 1:length(spans1)){
            for (j in 1:length(spans2)){
                #message("1")
                #message(spans1[i])
                #message("2")
                #message(spans2[j])
            minDist= min(minDist, getMinDist(spans1[i], spans2[j], checkCDS=FALSE))
            }
        }
        return(minDist)    
    }
    
    requestGTFbyProteinName=function(inaccession, inname){
         queryTxt=sprintf('SELECT * FROM GTFs WHERE accession="%s" AND type="CDS" AND Name="%s"', 
                         inaccession, inname)
        GTFsDF1name=dbGetQuery(mydbRO, queryTxt)
       #message("\t1a")
        queryTxt=sprintf('SELECT * FROM GTFs WHERE accession="%s" AND type="CDS" AND gene="%s"', 
                         inaccession, inname)
        GTFsDF1gene=dbGetQuery(mydbRO, queryTxt)
       #message("\t1b")
        queryTxt=sprintf('SELECT * FROM GTFs WHERE accession="%s" AND type="CDS" AND protein_id="%s"', 
                         inaccession, inname)
        GTFsDF1prot=dbGetQuery(mydbRO, queryTxt)
       #message("\t1c")
        GTFsDF1=rbind(GTFsDF1name, GTFsDF1gene)
        GTFsDF1=rbind(GTFsDF1, GTFsDF1prot)
        return(GTFsDF1)
    }
    
    getLastCDS=function(inDF){
        if (nrow(inDF)==0) return(NULL)
        accessionTMP=inDF$accession[1]
        inDF=as.data.frame( GenomicRanges::reduce(GRanges(inDF)) )
        if (inDF$strand[1]=="+"){
            #message("+")
            lastCDS=as.list(inDF[ order(-inDF$end), ][1,])
        }else{
            #message("-")
            lastCDS=as.list(inDF[ order(inDF$start), ][1,])
        }
        lastCDS$accession=accessionTMP
        return(lastCDS)
    }
    
    requestDoqnstreamCDSSQL=function(inexon, intable, excludetable=NULL){
        excludefilter=""
        if (! is.null(excludetable)){
            excludefilter = sprintf('AND NOT ( gene IN ( SELECT gene FROM suspectedTEgenes WHERE accession="%s") OR Name IN ( SELECT gene FROM suspectedTEgenes WHERE accession="%s"))', inexon$accession, inexon$accession)
        }
        if (inexon$strand=="+"){
            queryTxt=sprintf('SELECT * FROM %s WHERE accession="%s" AND seqnames="%s" AND ( start>%d OR end>%d) AND type="CDS" 
                            %s
                            ORDER BY start LIMIT 1',
                        intable, inexon$accession, inexon$seqnames, inexon$end, inexon$end, excludefilter)
        }else{
            queryTxt=sprintf('SELECT * FROM %s WHERE accession="%s" AND seqnames="%s" AND ( start<%d OR end<%d) AND type="CDS" %s ORDER BY end DESC LIMIT 1',
                        intable, inexon$accession, inexon$seqnames, inexon$start, inexon$start, excludefilter)
        }
        #message(queryTxt)
        tmpDF=dbGetQuery(mydbRO, queryTxt)
        return(tmpDF)
    }
    
    getCDSdowstream=function(lastExon, nextCDS){
        if (typeof(lastExon)=="list") lastExon=as.data.table(lastExon)
        if (typeof(nextCDS)=="list") nextCDS=as.data.table(nextCDS)
        chrstart=NA
        chrend=NA
        if (nrow(nextCDS)==0){
            # the last of first CDS on CHR
            if (lastExon$strand[1]=="+"){
                queryTxt=sprintf('SELECT length FROM seqinfo WHERE accession="%s" AND seqnames="%s"',
                                 lastExon$accession, lastExon$seqnames)
                tmpOUT=dbGetQuery(mydbRO, queryTxt)
                if (nrow(tmpOUT)==1){
                    chrend=tmpOUT[1,1]
                }else{
                    message("Cannot get chr length")
                    return(NULL)
                }
            }else{
                chrstart=1
            }
        }
        tmpGR=GenomicRanges::reduce(c(GRanges(lastExon),GRanges(nextCDS) ))
        theGap=GenomicRanges::gaps(tmpGR, start=chrstart, end=chrend, ignore.strand=TRUE)
        strand(theGap)=lastExon$strand[1]
        if (length(theGap)!=1){
            message("No gap found")
            return(NULL)
        }
        return(theGap)
    }
    
    getSeqIndexedFA=function(inFAname, inChr, inStart, inLen){
        CONtmp = file(inFAname,"rb")
        indDFtmp=read.table(paste0(inFAname,".fai"))
        colnames(indDFtmp)=c("seqnames", "seqlen", "chrStartByte", "NTperLine", "BytesPerLine")
        targetRow=indDFtmp[ indDFtmp$seqnames==inChr,]
        if (nrow(targetRow)!=1){
            stop("CHR name not found in index: ", inChr, " ", inFAname)
        }
        if (targetRow$seqlen<(inStart + inLen)){
            message("CHR too short: ", inChr, " ", inFAname, " " , inStart, " " , inLen)
            return(NULL)
        }
        nLinestmp=floor(inStart/targetRow$NTperLine)
        BytesInLines=nLinestmp*targetRow$BytesPerLine
        extraNTs=inStart-nLinestmp*targetRow$NTperLine
        coordinateFrameShift=2
        if (inStart==1) coordinateFrameShift = 1 # we do not want to capture the chr name
        BytesShift = targetRow$chrStartByte + BytesInLines + extraNTs  -coordinateFrameShift #(-1) due to zero based shift
        seek(CONtmp,BytesShift)
        tmpBuf=readBin(CONtmp,"raw",inLen*2)
        if (tmpBuf[2]!=0x0a){
            tmpBuf=tmpBuf[2:length(tmpBuf)]
        }
        #message(tmpBuf)
        outChar=rawToChar(tmpBuf)
        #message(outChar)
        outChar=str_replace_all(string = outChar, pattern = "\n", replacement = "")
        outChar=substring(outChar,1,inLen)
        if (! is.na(str_match(string = outChar, pattern = "[^ATCGNRYKMSWBDHVatcgnrykmswbdhv]")[1,1])){
            message(outChar)
            message(" ")
            message(paste0(str_extract_all(string=outChar, pattern="[^ATCGNRYKMSWBDHVatcgnrykmswbdhv]")[[1]], collapse=", "))
            stop("Non NT symbol: ", inChr, " ", inFAname, " " , inStart, " " , inLen)
        }
        close.connection(CONtmp)
        return(outChar)
    }
    
    getFApath=function(inAccession){
        genomeDir=paste0(StaticInpDir, inAccession,"/")
        #message(genomeDir)
        fs=list.files(genomeDir, pattern = paste0(inAccession,".*_genomic.fna$"), full.names = TRUE)
        if (length(fs)!=1){
            message(inAccession)
            message(fs)
            message(paste(inAccession, "PROBLEM"))
        }
        return(fs)
    }
    
    getChrLenSQL=function(accession, chrname){
        queryTxt=sprintf('SELECT length FROM seqinfo WHERE accession="%s" AND seqnames="%s"',
                                 accession, chrname)
        tmpOUT=dbGetQuery(mydbRO, queryTxt)
        if (nrow(tmpOUT)==1){
            return(tmpOUT[1,1])
        }else{
            message("Cannot get chr length")
            return(NA)
        }
    }
    
    distToChrEnd=function(inExon){
        if (inExon$strand=="+"){
            chrLenTMP=getChrLenSQL(inExon$accession, inExon$seqnames)
            return(chrLenTMP-inExon$end)
        }
        return(inExon$start -1)
    }
    
    getDSranges=function(inExon, DSlen, inTable){
        if (inExon$strand=="+"){
            searchStart=inExon$end
            searchEnd=searchStart+DSlen
        }else{
            searchStart=inExon$start-DSlen
            searchEnd=inExon$start
        }
    
        queryTxt=sprintf('SELECT * FROM %s WHERE ( ( start BETWEEN %d AND %d ) OR ( end BETWEEN %d AND %d ) ) AND accession="%s" AND seqnames="%s"',
                                 inTable, searchStart, searchEnd, searchStart, searchEnd, inExon$accession, inExon$seqnames)
        outGR=GRanges(dbGetQuery(mydbRO, queryTxt))
        return(outGR)
    }
    
    getDSinterval=function(inExon, DSlen){
        if (inExon$strand=="+"){
            searchStart=inExon$end
            searchEnd=searchStart+DSlen
        }else{
            searchStart=inExon$start-DSlen
            searchEnd=inExon$start
        }
    
        outGR=GRanges(seqnames=inExon$seqnames, ranges=IRanges(start=searchStart, end=searchEnd), strand=inExon$strand)
        return(outGR)
    }
}
getGenesTranscripsCharacteristicsFromGTFs<-function(inGTFname){
    # Reads GTF file, returns dataframe 1 of geneID, transcriptID, gene_version, transcript_version, gene_name, 
    # gene_source, gene_biotype, transcript_name, transcript_source, transcript_biotype, tag
        TMP.GTF=import.gff(inGTFname)
        TMP.GTF$type[ is.na(TMP.GTF$type)]=""
        TMP.GTF=TMP.GTF[ TMP.GTF$type=="transcript"]
        TMP.DT=as.data.table(TMP.GTF)
        TMP.DT$seqnames=as.character(TMP.DT$seqnames)
        #checking if all the fields are present
        keyfields=c("gene_id", "transcript_id", "gene_version", "transcript_version", "gene_name", 
                    "gene_source", "gene_biotype", "transcript_name", "transcript_source", 
                    "transcript_biotype", "tag")
        if (! all(keyfields %in% colnames(TMP.DT))){
            return(NULL)
        }
        TMP.DT=TMP.DT[, ..keyfields]
        TMP.DT$tag[is.na(TMP.DT$tag)]=""
        TMP.DT$gene_version=as.numeric(TMP.DT$gene_version)
        TMP.DT$transcript_version=as.numeric(TMP.DT$transcript_version)
        return(TMP.DT)
}
# moved this command last since it is prone to fail.
#human_ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = Ensembl_Version)

get_table_prefix=function(inURL){
    for (t in rev(unlist(str_split(string = inURL, pattern = "/")))){
        if (nchar(t)>0){
            return(t)
        }
    }
}

testRepeatMaskerDB=function(inTaxonID){
    famdbpath="/usr/local/apps/repeatmasker/4.2.1/famdb.py"
    retval=system(paste0(famdbpath, " lineage ", inTaxonID), wait=TRUE, intern = TRUE)
    if (length(retval)>0) return(TRUE)
    return(FALSE)
}

getTaxonID=function(inAccession){
    queryTxt=sprintf('SELECT * FROM taxonomy WHERE accession="%s"',inAccession)
    taxDFtmp=dbGetQuery(mydbRO, queryTxt)
    return(taxDFtmp[ order(-taxDFtmp$orderTax),])
}

getWorkingRMSKTaxonID=function(inAccession, maxRankDelta=5, DEBUG=FALSE){
    tmpTaxDF=getTaxonID(inAccession)  
    workingTaxon=NULL
    for (j in 1:max(maxRankDelta,nrow(tmpTaxDF))){
        taxIDtmp=tmpTaxDF$id[j]
        if (DEBUG) message("Testing: ", taxIDtmp)
        if (testRepeatMaskerDB(taxIDtmp)){
            workingTaxon=taxIDtmp
            return(workingTaxon)
        }
    }
    return(workingTaxon)
}