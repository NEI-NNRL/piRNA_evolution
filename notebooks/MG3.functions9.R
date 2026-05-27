suppressMessages({
    library(ggtree)
    library(treeio)
    library(ape)
    library(tidyverse)
    library(ggnewscale)
    library(ggtreeExtra)
    library(taxize)
    library(digest)
})

fixSVG<-function(inFname){
        #fixing SVG
        con=file(inFname,open="r")
        lines=readLines(con) 
        singleString <- paste(lines, collapse="\n")
        singleString<-str_replace_all(singleString, "textLength='[^']*'\\s*", "")
        singleString<-str_replace_all(singleString, 'font-family:[[:space:]]*"[^"]+"', 'font-family: "Helvetica"')
        singleString<-str_replace_all(singleString, 'stroke: #333333', 'stroke: #000000')
        singleString<-str_replace_all(singleString, 'stroke: #666666', 'stroke: #000000')
    singleString<-str_replace_all(singleString, ';fill: #4D4D4D; font-family:', ';fill: #000000; font-family:')
        singleString<-str_replace_all(singleString, 'stroke-width: "[^"]+"', 'stroke-width: 0.5')
        close(con)
        message(nchar(singleString))
        fileConn<-file(inFname)
        writeLines(singleString, fileConn)
        close(fileConn)
}



saveFigurePlain<-function(GGfigure, figName, outputDir, width.in, height.in ){

    ggsave(filename = paste0(figName, ".pdf"),path = outputDir , width = width.in,
    height = height.in,
    units = "in",
    plot = GGfigure,device = cairo_pdf)

 svgname=paste0(outputDir,"/",figName,".svg")
    ggsave(filename = svgname, width = width.in,
      height = height.in,
      units = "in",
    plot = GGfigure,device = "svg", bg="transparent")
    fixSVG(svgname)
 options(repr.plot.width=width.in, repr.plot.height=height.in)

    
}


get_rank2<-function(x, getNames=FALSE){
    # Based on ropensci/taxize/get_rank function by Dr. Vinh Tran 
    BiotaTaxID=131567
    rank_df <- x[, 'rank']
    rank_df[x[,'id']==BiotaTaxID]="Biota"
      
    shiftNumb=0
    prevKnownTaxRank=""
    for (i in 1:length(rank_df)){
        if (rank_df[i] %in% c(taxize:::mainTaxonomyRank(),"Biota")){
            shiftNumb=0
            prevKnownTaxRank=rank_df[i]
        }else{
            shiftNumb=shiftNumb+1
            rank_df[i]=paste0(prevKnownTaxRank,".plus",shiftNumb)
        }
    
    }
    names(rank_df) <- rank_df
    id_df <- x[, 'id']
    nameDf <- x[, 'name']
    if (getNames){
          joined_df <- cbind(
        data.frame(rank_df,stringsAsFactors=FALSE),
        data.frame(nameDf,stringsAsFactors=FALSE)
      )
      joined_df$name <- paste0(joined_df$nameDf, "##", joined_df$rank_df)
    }else{
      joined_df <- cbind(
        data.frame(rank_df, stringsAsFactors = FALSE),
        data.frame(id_df, stringsAsFactors = FALSE)
      )
    }
      
    if (getNames){
         df <- data.frame(t(data.frame(rev(joined_df$name))), stringsAsFactors = FALSE) 
    }else{
        df <- data.frame(
        t(data.frame(rev(joined_df$rank_df))), stringsAsFactors = FALSE)
    }
      outDf <- data.frame(tip = x[nrow(x), "name"], df, stringsAsFactors = FALSE)
    return(outDf)
}

my_get_id2 <- function (x) {
    # Based on ropensci/taxize/get_rank function  by Dr. Vinh Tran
  BiotaTaxID=131567
  rank_df <- x[, 'rank']
    rank_df[x[,'id']==BiotaTaxID]="Biota"

    shiftNumb=0
    prevKnownTaxRank=""
    for (i in 1:length(rank_df)){
        #message(rank_df[i])
        if (rank_df[i] %in% c(taxize:::mainTaxonomyRank(),"Biota")){
            shiftNumb=0
            prevKnownTaxRank=rank_df[i]
        }else{
            shiftNumb=shiftNumb+1
            rank_df[i]=paste0(prevKnownTaxRank,".plus",shiftNumb)
        }

    }
    names(rank_df) <- rank_df

  nameDf <- x[, 'id']
  id_df <- x[, 'id']

  joined_df <- cbind(
    data.frame(rank_df,stringsAsFactors=FALSE),
    data.frame(nameDf,stringsAsFactors=FALSE)
  )
  joined_df$name <- paste0(joined_df$nameDf, "##", joined_df$rank_df)

  df <- data.frame(t(data.frame(rev(joined_df$name))), stringsAsFactors = FALSE)
  outDf <- data.frame(tip = x[nrow(x), "name"], df, stringsAsFactors = FALSE)
  return(outDf)
}

taxonomy_table_creator2 <- function (nameList, rankList) {
    #Based on ropensci/taxize/taxonomy_table_creator originally written by Dr. Vinh Tran 
    colnames(nameList)[1] <- "tip"
    # get indexed rank list
    index2RankDf <- taxize:::rank_indexing(rankList)
    # get ordered rank list
    orderedRank <- factor(index2RankDf$rank, levels = index2RankDf$rank)
    # create a dataframe containing ordered ranks
    i=length(orderedRank)
    orderedRank2=c()
    while (i>0){
        if (str_count(orderedRank[i], ".plus")==0){
            orderedRank2=c(as.character(orderedRank[i]), orderedRank2)
            #finding plus ones
            plusPattern=paste0(orderedRank[i],".plus")
            plusOnes=grep(pattern = paste0("^",plusPattern), x= orderedRank)
            if (length(plusOnes)>0){
                foundSubRanks=as.integer(gsub(pattern = plusPattern, replacement = "", x=orderedRank[plusOnes]))
                foundSubRanks=sort(foundSubRanks)
                for (subRank in foundSubRanks){
                    orderedRank2=c(paste0(plusPattern,subRank), orderedRank2)
                }
            }
        }
        i=i-1
    }
    orderedRank=factor(orderedRank2, levels=orderedRank2)
    
  full_rank_name_df <- data.frame(
    "rank"= matrix(unlist(orderedRank), nrow = length(orderedRank), byrow=TRUE),
    stringsAsFactors = FALSE
  )
    
  full_rank_name_df$index <- as.numeric(rownames(full_rank_name_df))
  
    
  fullRankIDdf <- data.frame(
    rank = matrix(
      unlist(orderedRank), nrow = length(orderedRank), byrow = TRUE
    ), stringsAsFactors = FALSE)
  fullRankIDdf$index <- as.numeric(rownames(fullRankIDdf))

  for (i in 1:nrow(nameList)) {
    ### get list of all IDs for this taxon
    taxonDf <- data.frame(nameList[i,])
    taxonName <- unlist(
      strsplit(as.character(nameList[i,]$tip), "##", fixed = TRUE)
    )

    ### convert into long format
    # mTaxonDf <- suppressWarnings(melt(taxonDf,id = "tip"))
    mTaxonDf <- data.table::setDF(
      suppressWarnings(
        data.table::melt(data.table::as.data.table(taxonDf), id = "tip")
      )
    )

    ### get rank names and corresponding IDs
    splitCol <- data.frame(
      do.call(
        'rbind', strsplit(as.character(mTaxonDf$value), "##", fixed=TRUE)
      )
    )
    mTaxonDf <- cbind(mTaxonDf,splitCol)

    ### remove NA cases
    mTaxonDf <- mTaxonDf[complete.cases(mTaxonDf),]

    ### subselect mTaxonDf to keep only 2 column rank id and rank name
    mTaxonDf <- mTaxonDf[,c("X1","X2")]
    if(mTaxonDf$X2[1] != index2RankDf$rank[1] &&
       !index2RankDf$rank[1] %in% mTaxonDf$X2){
      mTaxonDf <- rbind(
        data.frame("X1"=mTaxonDf$X1[1], "X2"=index2RankDf$rank[1]), mTaxonDf
      )
    }

    ### rename columns
    colnames(mTaxonDf) <- c(taxonName[1],"rank")

    ### merge with index2RankDf (Df contains all available ranks of input data)
    full_rank_name_df <- merge(
      full_rank_name_df,mTaxonDf, by=c("rank"), all.x = TRUE
    )

    ### reorder ranks
    full_rank_name_df <- full_rank_name_df[order(full_rank_name_df$index),]

    ### replace NA id by id of previous rank
    full_rank_name_df <- zoo::na.locf(full_rank_name_df)
  }

  ### remove index column
  full_rank_name_df <- subset(full_rank_name_df, select = -c(index))

  ### transpose into wide format
  t_full_rank_name_df <- data.table::transpose(full_rank_name_df)

  ### set first row to column names
  colnames(t_full_rank_name_df) = as.character(unlist(t_full_rank_name_df[1,]))
  t_full_rank_name_df <- t_full_rank_name_df[-1,]

  ### replace NA values in the dataframe t_full_rank_name_df
  if(nrow(t_full_rank_name_df[is.na(t_full_rank_name_df),]) > 0){
    t_full_rank_name_dfTMP <- 
      t_full_rank_name_df[complete.cases(t_full_rank_name_df),]
    t_full_rank_name_dfEdit <- t_full_rank_name_df[is.na(t_full_rank_name_df),]

    for(i in 1:nrow(t_full_rank_name_dfEdit)){
      for(j in 1:(ncol(t_full_rank_name_df)-1)){
        if(is.na(t_full_rank_name_dfEdit[i,j])){
          t_full_rank_name_dfEdit[i,j] <- t_full_rank_name_dfEdit[i,j+1]
        }
      }
      if(is.na(t_full_rank_name_dfEdit[i,ncol(t_full_rank_name_df)])){
        t_full_rank_name_dfEdit[i,ncol(t_full_rank_name_df)] <-
        t_full_rank_name_dfEdit[i,ncol(t_full_rank_name_df)-1]
      }
    }
    t_full_rank_name_df <- rbind(t_full_rank_name_dfEdit,t_full_rank_name_dfTMP)
  }

  ### add fullName column into t_full_rank_name_df
  fullName <- colnames(full_rank_name_df)[-c(1)]
  full_rank_name_df <- cbind(fullName,t_full_rank_name_df)

  ### return
  return(full_rank_name_df)
}
        
# converting class2tree into BringTaxaToCommonCoordinateFrame
inc2t=paste0(as.character(body(taxize::class2tree)), collapse = "\n")
c2t2change=str_extract(string = inc2t, pattern = regex("(.+df <- df\\[,\\s*-1\\])", dotall = TRUE), group = 1)
if (digest(c2t2change, algo = "md5", serialize = FALSE) != 'c31c51395e38ba18dbb362c3257be0d9'){
    stop("taxize::class2tree source code changed. MD5 sum does not match v0.10.0 . Contact me for BringTaxaToCommonCoordinateFrame source code")
}
t0=str_extract(string = c2t2change, pattern = "rankList <- dt2df[^\n]+")
t=str_replace(string = t0, pattern = "dt2df", replacement = "taxize:::dt2df")
t=str_replace(string = t, pattern = "get_rank", replacement = "get_rank2")
c2t2change=str_replace(string = c2t2change, pattern = fixed(t0), replacement = t)
t0=str_extract(string = c2t2change, pattern = "nameList <- dt2df[^\n]+")
t="if (!indexesOnly){"
t1=str_replace(string = t0, pattern = "dt2df", replacement = "taxize:::dt2df")
t2a=str_replace(string = t1, pattern = "get_name", replacement = "get_rank2, getNames = TRUE")
t2b=str_replace(string = t1, pattern = "get_name", replacement = "my_get_id2")
t=paste(t, t2a, '}else{', t2b, "}", sep = "\n")
c2t2change=str_replace(string = c2t2change, pattern = fixed(t0), replacement = t)
c2t2change=str_replace(string = c2t2change, pattern = fixed('taxonomy_table_creator'), replacement = 'taxonomy_table_creator2')
c2t2change=str_replace(string = c2t2change, pattern = ',\\s*remove_shared[^\\)]+\\)', replacement = ')')
c2t2change=paste0("BringTaxaToCommonCoordinateFrame<-function(input, varstep = FALSE, check = FALSE, indexesOnly=FALSE)\n", c2t2change, '\nreturn(df)\n}', collapse = '\n')
eval(parse(text =c2t2change))
#done modidying
                
recursiveTreeBranchFinder<-function(currentLevel, prevTaxon, inDF, isroot=FALSE){
    childrenTaxons=unique(inDF[inDF[,currentLevel+1]==prevTaxon,currentLevel])
    if (length(childrenTaxons)>1){
        returnBranchingLevels=c(currentLevel+1)
    }else{
        returnBranchingLevels=c()
    }
    
    #returnXMLlist=list()
    if (currentLevel>1){
        for (chName0 in childrenTaxons){
            deeperReturns=recursiveTreeBranchFinder(currentLevel-1, chName0, inDF)
            returnBranchingLevels=c(deeperReturns, returnBranchingLevels)
        }
    }
    
    return(returnBranchingLevels)
}
                
recursiveTreeBuilderV4<-function(currentLevel, prevTaxon, inDF, in.tax.DF, inDataDF, isroot=FALSE, parentID=""){
    childrenTaxons=unique(inDF[inDF[,currentLevel+1]==prevTaxon,currentLevel])
    childrenTaxons.id=unique(in.tax.DF[inDF[,currentLevel+1]==prevTaxon,currentLevel])
    myTaxonID=unique(in.tax.DF[inDF[,currentLevel+1]==prevTaxon,currentLevel+1])
    myTaxonRank=colnames(inDF[currentLevel+1])
    MaskStatus=0
    returnTreeString=""
    MyAverageData=NULL
    if (currentLevel>1){
        returnXMLlist=list()
        for (childNum in 1:length(childrenTaxons)){
            chName0=childrenTaxons[childNum]
            chID0=childrenTaxons.id[childNum]
            #message(paste0(myTaxonRank, collapse = ", "))
            if ((length(childrenTaxons)==1) && (myTaxonID==chID0) && (str_count(myTaxonRank, ".plus")>0)) {
                MaskStatus=1
            }
            chName=sub(pattern = " ",replacement = "_",x = chName0,fixed = TRUE)
            deeperList=recursiveTreeBuilderV4(currentLevel-1, chName0, inDF, in.tax.DF, inDataDF, parentID=childNum)
            deeperText=deeperList[["text"]]
            deeperXML=deeperList[["data"]]
            deeperAverageData=deeperList[["averageData"]]
            chID=paste0("{",chName,"@",currentLevel,'-',parentID,"}")
            returnTreeString=paste0(returnTreeString,",",deeperText,":1",chID)
            returnXMLlist=c(returnXMLlist, deeperXML)
            XMLrecordTMP=paste0('{"edge_num": ',chID,',"isTip": 0,"text": "',chName,'"')
            XMLrecordTMP=paste0(XMLrecordTMP,',"TaxonLevel": ',currentLevel,',"taxID":',chID0,',"isMasked":', deeperList[["isMasked"]])
            #message(chName0)
            #message(deeperAverageData)
            for (cnumb in 1:ncol(deeperAverageData)){
                cname=colnames(deeperAverageData)[cnumb]
                #message(cname)
                cval=deeperAverageData[1,cnumb]
                XMLrecordTMP=paste0(XMLrecordTMP,',"',cname,'": "',cval,'"')
            }
            returnXMLlist[[length(returnXMLlist)+1]]=paste(XMLrecordTMP,"}")
            if (is.null(MyAverageData)){
                MyAverageData=deeperAverageData
            }else{
                MyAverageData=rbind(MyAverageData,deeperAverageData)
            }
        }
    }else{
        returnXMLlist=list()
        for (childNum in 1:length(childrenTaxons)){
            chName0=childrenTaxons[childNum]
            chID0=childrenTaxons.id[childNum]
            #message(typeof(chID0))
            specieData=inDataDF[inDataDF$taxon == chID0,]
            #message(specieData)
            chName=sub(pattern = " ",replacement = "_",x = chName0,fixed = TRUE)
            chID=paste0("{",chName,"@",currentLevel,'-',parentID,"}")
            returnTreeString=paste0(returnTreeString,',',chName,':1',chID)
            XMLrecordTMP=paste0('{"edge_num": ',chID,',"isTip": 1,"text": "',chName,'"')
            for (i in 1:ncol(specieData)){
                TMPcolname=colnames(specieData)[i]
                TMPcolData=specieData[1,i]
                XMLrecordTMP=paste0(XMLrecordTMP,',"',TMPcolname,'": "',TMPcolData,'"')
            }
            XMLrecordTMP=paste0(XMLrecordTMP,',"TaxonLevel": ',currentLevel,',"taxID":',chID0,',"isMasked":', MaskStatus)
            returnXMLlist[[length(returnXMLlist)+1]]=paste(XMLrecordTMP,"}")
            #adding children numeric values
            numbercBooleanMask=c()
            for (colInd in 1:ncol(specieData)){
                numbercBooleanMask[ length(numbercBooleanMask)+1] = is.numeric(specieData[1,colInd]) & (colnames(specieData)[colInd]!="taxonId")
            }
            specieData.numeric=specieData[1,numbercBooleanMask]
            if (is.null(MyAverageData)){
                MyAverageData=specieData.numeric
            }else{
                MyAverageData=rbind(MyAverageData,specieData.numeric)
            }
        }
    }
    if (! isroot){
        returnTreeString=paste0("(",substring(returnTreeString,2,nchar(returnTreeString)),")")
    }else{
        chID=paste0("{",prevTaxon,"@",currentLevel+1,"-0}")
        rootID=unique(in.tax.DF[inDF[,currentLevel+1]==prevTaxon,currentLevel+1])
        returnTreeString=paste0("(",substring(returnTreeString,2,nchar(returnTreeString)),")",":1",chID,";")
        returnXMLlist[[length(returnXMLlist)+1]]=paste0('{"edge_num": ',chID,',"isTip": 0,"text": "',prevTaxon,'","TaxonLevel": ',currentLevel+1,
                                                        ',"taxID":',rootID,', "isMasked": 0}')
    }
    MyAverageDataCollapsed=apply(MyAverageData, MARGIN = 2, FUN = mean)
    
    return(list(text=returnTreeString, data=returnXMLlist, isMasked=MaskStatus, averageData=t(as.data.frame(MyAverageDataCollapsed))))
}                


selectAssembly=function(tmpSubset){
    if (nrow(tmpSubset)==1) return(tmpSubset)
    if (!all(tmpSubset$DSlen>1e6)){ # DSlen selection
        tmpSubset=tmpSubset[ tmpSubset$DSlen == max(tmpSubset$DSlen), ]
    }
    if (nrow(tmpSubset)==1) return(tmpSubset)
    # minimal gap
    tmpSubset=tmpSubset[ tmpSubset$GapLen/tmpSubset$DSlen == min(tmpSubset$GapLen/tmpSubset$DSlen), ]
    if (nrow(tmpSubset)==1) return(tmpSubset)
    # max asRMSKlen
    tmpSubset=tmpSubset[ tmpSubset$asRMSKlen == max(tmpSubset$asRMSKlen), ]
    if (nrow(tmpSubset)==1) return(tmpSubset)
    # min sRMSKlen
    tmpSubset=tmpSubset[ tmpSubset$sRMSKlen == min(tmpSubset$sRMSKlen), ]
    return(tmpSubset[1,])
}
                

                

                