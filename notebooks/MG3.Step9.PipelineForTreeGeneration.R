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
source(paste0(sourceDir,'/MG3.functions9.R'))

args = commandArgs(trailingOnly=TRUE)

HMMid=args[1]

JOBid=Sys.getenv("SLURM_JOB_ID")
tmpDir=paste0('/lscratch/',JOBid,"/")

queryTxt=sprintf('SELECT * FROM HMMstats5MB WHERE ensembl_gene_id="%s"', HMMid)
HMMstatsDT=as.data.frame(dbGetQuery(mydbRO, queryTxt))

queryTxt='SELECT taxonomy.accession, taxonomy.id AS taxon FROM taxonomy INNER JOIN (SELECT accession, MAX(orderTax) AS maxOrd FROM taxonomy GROUP BY accession) AS maxrank 
        ON taxonomy.accession=maxrank.accession WHERE 
            taxonomy.orderTax=maxrank.maxOrd  AND 
            taxonomy.accession IN ( SELECT * FROM EutherianGenomes )'
taxID.DF=dbGetQuery(mydbRO, queryTxt)
HMMstatsDT=merge(HMMstatsDT, taxID.DF, by="accession", all.x = TRUE, all.y = FALSE)

# iterating over taxons and selecting the best assembly
OnePerTax=NULL
for (taxon in unique(HMMstatsDT$taxon)){
    boolmask=HMMstatsDT$taxon ==taxon
    tmpSubset=HMMstatsDT[ boolmask, ]
    OnePerTax=addToDF(OnePerTax, selectAssembly(tmpSubset))
}

queryTxt=sprintf('SELECT DISTINCT * FROM taxonomy WHERE accession IN ("%s")',
                    paste0(OnePerTax$accession, collapse = '", "'))
PresentTaxizeDF=dbGetQuery(mydbRO, queryTxt)

PresentTaxizeDF=merge(PresentTaxizeDF, OnePerTax[, c("accession", "taxon")], by="accession", all.x = TRUE)

PresentTaxizeList=split(PresentTaxizeDF, PresentTaxizeDF$taxon)
PresentTaxizeList=lapply(PresentTaxizeList, FUN = function(x) {return(x[ order(x$orderTax),])})

AlignedTaxaDF=BringTaxaToCommonCoordinateFrame( PresentTaxizeList )
AlignedTaxaDF.id=BringTaxaToCommonCoordinateFrame( PresentTaxizeList, indexesOnly = TRUE )

currentLevel=ncol(AlignedTaxaDF)-1
prevTaxon=AlignedTaxaDF[1,ncol(AlignedTaxaDF)]

AlignedTaxaBranchingLevelsRaw=recursiveTreeBranchFinder(currentLevel, prevTaxon, AlignedTaxaDF, TRUE)
AlignedTaxaBranchingLevels=unique(sort(AlignedTaxaBranchingLevelsRaw))
AlignedTaxaLevelsToCheck=sort(unique(c(1, AlignedTaxaBranchingLevels)))

AlignedTaxaDFSubset=AlignedTaxaDF[, AlignedTaxaLevelsToCheck]
AlignedTaxaDF.idSubset=AlignedTaxaDF.id[, AlignedTaxaLevelsToCheck]

currentLevel=ncol(AlignedTaxaDFSubset)-1
prevTaxon=AlignedTaxaDFSubset[1,ncol(AlignedTaxaDFSubset)]
message(prevTaxon)

returnTreeList=recursiveTreeBuilderV4(currentLevel, prevTaxon, AlignedTaxaDFSubset, AlignedTaxaDF.idSubset, OnePerTax, TRUE, parentID = 1)
returnTreeString=returnTreeList[["text"]]
returnTreeXML=returnTreeList[["data"]]

edgeMarks=str_extract_all(string = returnTreeString, pattern = "\\{[^\\}]+@[0-9]+-[0-9]+\\}")
i=1
for (edgeMark  in unlist(edgeMarks)){
    returnTreeString=gsub(edgeMark, paste0("{",as.character(i),"}"), returnTreeString,fixed = TRUE)
    returnTreeXMLnew=gsub(edgeMark, as.character(i), returnTreeXML,fixed = TRUE)
    if (paste0(returnTreeXML, collapse = " ")==paste0(returnTreeXMLnew, collapse = " ")){
        message(edgeMark)
    }
    returnTreeXML=returnTreeXMLnew
    i=i+1
}
outputJSON=paste0('{\n"tree": "',returnTreeString,'",\n"data" : [')
for (XMLrecord in returnTreeXML){
    outputJSON=paste0(outputJSON,XMLrecord,"\n,")
}
outputJSON=substring(outputJSON,1,nchar(outputJSON)-1)
outputJSON=paste0(outputJSON,'],\n"metadata": {"info": "script by Alex", "data": "',getTimeStamp(),'"}}')
treeFileName=paste0(tmpDir,"tree",Sys.getpid(), ".json")
fileConn<-file(treeFileName)
writeLines(outputJSON, fileConn)
close(fileConn)


jTREEtmp=read.jtree(file=treeFileName)
jTREEtmp@data$TaxonLevel<-as.integer(jTREEtmp@data$TaxonLevel)
jTREEtmp@data[jTREEtmp@data$node, ] = jTREEtmp@data
TEasNormDetectedTillNextCDS=as.numeric(jTREEtmp@data$asRMSKlen)/as.numeric(jTREEtmp@data$DSlen)
TEasNormDetectedTillNextCDS[abs(TEasNormDetectedTillNextCDS)==Inf]=0

isTip=(jTREEtmp@data$isTip==1)
MaskedNodes=jTREEtmp@data$node[jTREEtmp@data$isMasked==1] 


#jTREEtmp@phylo$tip.label=gsub(pattern = "_",replacement = " ",x = jTREEtmp@phylo$tip.label,fixed = TRUE)
#jTREEtmp@phylo$tip.label=paste(jTREEtmp@data$scientific.name[ jTREEtmp@data$node %in% 1:length(jTREEtmp@phylo$tip.label) ], "(",jTREEtmp@data$common.name[ jTREEtmp@data$node %in% 1:length(jTREEtmp@phylo$tip.label) ],")")

Nodelabels=rep("!!!!", max(jTREEtmp@data$node))
Nodelabels=jTREEtmp@data$text
Nodelabels[ 1:length(jTREEtmp@phylo$tip.label) ] = round(TEasNormDetectedTillNextCDS[isTip], 3)
Nodelabels[ MaskedNodes ] = NA


#adding the fixed (10K, 50K etc) dowstream length data
tmpTipData=jTREEtmp@data[ jTREEtmp@data$isTip==1,]
TipColnames<-colnames(tmpTipData)
tmpTipData$ID=jTREEtmp@phylo$tip.label
DSfixedLens=TipColnames[ grepl("asRMSKlen\\d+K", TipColnames)]
DSfixedLens=as.integer(str_remove_all(string = DSfixedLens, pattern = "[^0-9]"))
tipDownstreamTEasData=NULL
for (DSlen in DSfixedLens){
    asVal=as.numeric(tmpTipData[[sprintf("asRMSKlen%dK", DSlen)]])
    asFraction=asVal/(DSlen*1000)
    IDs=tmpTipData$ID
    tmpDF=data.frame(ID=IDs, asFraction=asFraction, DSlenFixed=DSlen)
    tipDownstreamTEasData=addToDF(tipDownstreamTEasData, tmpDF)
}


myTreeGG=ggtree::ggtree(jTREEtmp)+ geom_tiplab(as_ylab=TRUE, color='black')+
geom_label(aes(x=branch, label=Nodelabels, color=as.numeric(TEasNormDetectedTillNextCDS)))+ 
scale_color_gradient2(limits=c(0,max(as.numeric(na.exclude(as.numeric(TEasNormDetectedTillNextCDS))))),
                      oob = scales::squish,
                      low = "green", mid = "green", high = "red")

myTreeGG$data$DSlen[ MaskedNodes ] = NA
distanceToNextCDSLabels=round(as.numeric(myTreeGG$data$DSlen )/1000, digits = 1)
distanceToNextCDSLabels[ myTreeGG$data$isTip == 0] = NA

#adding taxonomy labels
curX=0
for (myColname in rev(colnames(AlignedTaxaDFSubset))){
    if (str_count(myColname, ".plus")>0) myColname="clade"
    myTreeGG=myTreeGG + annotate(geom = "text", x=curX-0.5, y=0, label=myColname)
    curX=curX+1
}

myTreeGG=myTreeGG+ theme(plot.margin = unit(c(1, 1.5, 1, 1), "inches")) + theme(legend.position = "bottom")+
labs(colour = "TE anti-sense")+
geom_fruit(data=tipDownstreamTEasData, geom=geom_tile,
                  mapping=aes(y=ID, x=factor(DSlenFixed), fill=asFraction),
                  offset = 0.04,size = 0.02, pwidth=0.1)+
scale_fill_gradient2(limits=c(0,max(as.numeric(na.exclude(as.numeric(TEasNormDetectedTillNextCDS))))),
                      oob = scales::squish,
                      low = "green", mid = "green", high = "red")+
new_scale_colour() +
geom_point(aes(x = myTreeGG$data$x, y = myTreeGG$data$y, 
               size=as.numeric(myTreeGG$data$DSlen ),
               color=log2(as.numeric(myTreeGG$data$asRMSKlen )/as.numeric(myTreeGG$data$sRMSKlen ))
              ))+scale_color_gradient2(
                      low = "white", mid = "gray", high = "black")+
labs(colour = "Strandedness, log2")+
annotate(geom = "text", x = myTreeGG$data$x + 0.25, 
         y = myTreeGG$data$y, label=distanceToNextCDSLabels) + 
annotate(geom = "text", x = 0, y = max(myTreeGG$data$y), label = HMMid, vjust = 2, hjust = -2)


###
treesdir=paste0(OutputDir,"/step9Trees/")
system(paste0("mkdir -p ", treesdir), intern = FALSE, wait = TRUE)

ggsave(
  paste0(HMMid,".svg"),
  plot = myTreeGG,
  device = "svg",
  path = treesdir,
  scale = 1,
  width = 45,
  height = 40,
  units = "in",
  dpi = 300
)
#saveFigure(myTreeGG, HMMid, treesdir, 45, 40 )