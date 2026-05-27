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
accession=args[1]
message(accession)

SpecieDir=paste0(StaticInpDir, accession,"/")
GTFname=paste0(SpecieDir,"genomic.gff")
columnsToExport=c("seqnames", "start", "end", "width", "strand", "source", "type", "gene", "protein_id", "Name")

GTFtmp=readSpecieGTF(GTFname)
if (!is.null(GTFtmp)){
    GTFtmp$strand=as.character(GTFtmp$strand)
    GTFtmp=GTFtmp[,..columnsToExport]
    GTFtmp$accession=accession
    for (cn in colnames(GTFtmp)){
        if (typeof(GTFtmp[[cn]])==typeof("A")){
            GTFtmp[[cn]][ is.na(GTFtmp[[cn]])] = "?"
        }
    }
    
    write.csv(GTFtmp, file = paste0(OutputDir,"/tmpGTF2CSV/",accession,".csv"))

}else{
    message("NO GTF FILE")
}

