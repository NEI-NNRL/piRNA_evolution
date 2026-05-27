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

library(latex2exp)
library(fpc)


queryTxt=sprintf('SELECT * FROM HMMstats5MB WHERE ensembl_gene_id="%s"', HMMid)
HMMstatsDT=as.data.table(dbGetQuery(mydbRO, queryTxt))


statsDFtmp = HMMstatsDT %>% group_by(ensembl_gene_id, familyID) %>% filter(n()>2)

mycolors=c('#FF0000', '#0000FF', '#FFFF00', '#00FF00', '#FF00FF', '#FF8000', '#00FFFF', '#000000')
myfills=c('#FF0000', '#0000FF', '#FFFF00', '#00FF00', '#FF00FF', '#FF8000', '#00FFFF', '#000000', '#FFFFFF')

length(mycolors)*length(myfills)*5
visDF=expand.grid(color=mycolors, fill=myfills, shape=21:25)
visDF <- visDF[sample(nrow(visDF)),]

visDF$familyName=''
visDF$familyName[1:length(unique(HMMstatsDT$familyName))]=levels(factor(HMMstatsDT$familyName))
visDF=visDF[1:length(unique(HMMstatsDT$familyName)),]

colorList=as.list(setNames(as.character(visDF$color),visDF$familyName))
fillList=as.list(setNames(as.character(visDF$fill),visDF$familyName))

shapeList=as.list(setNames(visDF$shape,visDF$familyName))
shapeListClusters=list()
shapeListClusters[['0']]=NA
shapeListClusters[['1']]=0
shapeListClusters[['2']]=1
shapeListClusters[['3']]=2
shapeListClusters[['4']]=5
shapeListClusters[['5']]=6

library(ggnewscale)

HMMstatsDT$incolX=HMMstatsDT$DSlen/median(HMMstatsDT$DSlen)#(HMMstatsDT$DSlen)^0.35
HMMstatsDT$incolY=2*(HMMstatsDT$asRMSKlen-HMMstatsDT$sRMSKlen)/HMMstatsDT$DSlen
colYthreshold=10
#Ymask=HMMstatsDT$incolX<(4000^0.35)
epsilon=0.2
minOddsRatio=2
minMeanXforTopCluster=(30000^0.35)
#HMMstatsDT$incolY[ Ymask] =0
HMMstatsDT=HMMstatsDT[ !(is.na(HMMstatsDT$incolX) | is.na(HMMstatsDT$incolY)) ,]

Dbscan_cl=dbscan(HMMstatsDT[, c("incolX", "incolY")], eps = epsilon, MinPts = 10) 
HMMstatsDT$cluster=Dbscan_cl$cluster

my_gg_obj= ggplot() + geom_point(data = HMMstatsDT, aes(y=incolY, x=incolX, color=familyName, fill=familyName, shape=familyName)) +
    scale_shape_manual(values = shapeList)+
    scale_colour_manual(values = colorList)+
    scale_fill_manual(values = fillList) + # adding clusters
    new_scale_colour()+new_scale_fill()+new_scale("shape")+
    geom_point(data = HMMstatsDT, aes(y=incolY, x=incolX, shape=as.character(cluster), size=3))+
    scale_shape_manual(values = shapeListClusters) + scale_size(guide = 'none') + ggtitle(HMMid)


###
imagedir=paste0(OutputDir,"/scratches/clustering/")
system(paste0("mkdir -p ", imagedir), intern = FALSE, wait = TRUE)

ggsave(
  paste0(HMMid,".svg"),
  plot = my_gg_obj,
  device = "svg",
  path = imagedir,
  scale = 1,
  width = 20,
  height = 15,
  units = "in",
  dpi = 300
)


