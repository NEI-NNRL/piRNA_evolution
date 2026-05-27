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

inargs = commandArgs(trailingOnly=TRUE)
accession = inargs[1]
chromosomesToAnalyze = inargs[2]
writeDir= inargs[3]
#GCF_027475565.1 NC_069380.1 32536 ~/data/lambda/multigene2/Genomes/ ~/data/lambda/multigene2/RMSK/'

writeAndRun<-function(accession, chromosomesToAnalyze, writeDir){
    SpecieDir=paste0(StaticInpDir,'/', accession, "/")
    SpFiles<-list.files(path = SpecieDir)
    #message(paste0(SpFiles, collapse = ", \n"))
    GenomeFiles=str_extract(string = SpFiles, pattern = paste0(accession,".*_genomic.fna$"))
    #message(GenomeFiles)

    GenomeFiles = GenomeFiles[! is.na(GenomeFiles) ]
    #message(GenomeFiles)
    TargetGenome<-Rsamtools::FaFile(paste0(SpecieDir, GenomeFiles[1]))
    TargetGenomeSeq<-Rsamtools::getSeq(TargetGenome)
    TargetGenomeNames<-names(TargetGenomeSeq)
    TargetChrName<-TargetGenomeNames[!is.na(str_extract(string = TargetGenomeNames, 
                                                        pattern = paste0("^",chromosomesToAnalyze)))]
    message(TargetChrName)
    if (length(TargetChrName)!=1){
        error(paste(accession), "has few possible chromosomes for export")
    }
    TargetChrSeq=TargetGenomeSeq[TargetChrName]
    SpecieOutputDir=paste0(writeDir,"/", accession, "/")
    system(sprintf("mkdir %s -p",SpecieOutputDir), wait=TRUE)
    TargetChrFileName=paste0(SpecieOutputDir, GenomeFiles, "_", chromosomesToAnalyze, ".fna")
    #message(TargetChrFileName)
    writeXStringSet(TargetChrSeq, TargetChrFileName)
    taxonIdTMP=getWorkingRMSKTaxonID(accession)
    comToBash=paste0(sourceDir,"/RepeatMaskerChromosome.sh ",SpecieOutputDir," ",TargetChrFileName," ", taxonIdTMP)
    message(comToBash)
    retText=system(comToBash, wait=TRUE)
    message(paste0(retText, collapse=" \n"))
}

writeAndRun(accession, chromosomesToAnalyze, writeDir)
