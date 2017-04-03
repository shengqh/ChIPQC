#################################################################
#################################################################
## Utility
mergeMetadata <- function(object,addMetaData,facetBy,colourBy,lineBy){
  metadata <- QCmetadata(object)
  if(!is.null(addMetaData) & class(addMetaData)=="data.frame"){
    metadata <- merge(metadata,addMetaData,by.x=1,by.y=1,all.x=TRUE,all.y=FALSE)
  }
  if(all(facetBy %in% c("Sample",colnames(metadata)))){
    facet <- facet_wrap(
      formula(paste("~",paste(facetBy,collapse="+")))
    )
  }else{
    facet <- NULL  
  }
  if(all(colourBy %in% c(colnames(metadata),"Sample"))){
    colour <- aes_string(colour=colourBy)
  }else{
    colour <- NULL  
  }
  if(all(lineBy %in% c(colnames(metadata),"Sample"))){
    lineType <- aes_string(linetype=lineBy)
  }else{
    lineType <- NULL  
  }
  
  mDataList <- list(metadata,facet,colour,lineType)
  names(mDataList) <- c("metadata","facetBy","colour","lineType")
  return(mDataList)
}

############################
## CrossCoverage Plotting

makeCCplot <- function(CCDataFrame,shiftlength,readlen,excludedBox=F){
  CCDataFrame <- CCDataFrame[order(CCDataFrame$Shift_Size),]
  if(excludedBox){
    P <- ggplot(CCDataFrame,aes(x=Shift_Size,y=CC_Score))+geom_path(alpha = 1,size=1.3)+xlim(0,shiftlength)+ylab("CC_Score")+
      geom_rect(data=CCDataFrame,xmin=-Inf,colour="grey80",xmax=readlen*2,ymin=-Inf,ymax=Inf,fill="grey20",alpha = 0.002)+
      theme(axis.title.y=element_text(angle=0))
  }else{
    P <- ggplot(CCDataFrame,aes(x=Shift_Size,y=CC_Score))+geom_path(alpha = 1,size=1.3)+xlim(0,shiftlength)+ylab("CC_Score")+
      theme(axis.title.y=element_text(angle=0))     
  }
  return(P)
}

#setGeneric("plotCC", function(object="ChIPQCsample",method="Coverage") standardGeneric("plotCC"))

setGeneric("plotCC", function(object="ChIPQCexperiment",method="Coverage",facet=TRUE,
                              facetBy=c("Tissue","Factor"),
                              colourBy="Replicate",
                              lineBy=NULL,
                              addMetaData=NULL,excludedBox=F
)
  standardGeneric("plotCC")
)

setMethod("plotCC", "ChIPQCexperiment", function(object,method="Coverage",facet=TRUE,
                                                 facetBy=c("Tissue","Factor"),                                                 
                                                 colourBy="Replicate",
                                                 lineBy=NULL,
                                                 addMetaData=NULL,excludedBox=FALSE
){
  
  ccvector <- crosscoverage(object)
  readlen <- max(readlength(object))
  shiftlength <- nrow(ccvector)
  
  toMelt <- data.frame("Shift_Size"=seq(1,shiftlength),
                       #"metadataOfInterest"=metadataOfInterest,
                       ccvector)
  CCDataFrame <- melt(toMelt,id.vars=c("Shift_Size"))
  
  
  metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy,lineBy)        
  
  
  CCDataFrameWithMetaData <- merge(CCDataFrame,metadataOpts$metadata,by.x=2,by.y=1,all=FALSE)
  
  colnames(CCDataFrameWithMetaData)[1:3] <- c("Sample","Shift_Size","CC_Score")
  Plot <- makeCCplot(CCDataFrameWithMetaData,shiftlength,readlen,excludedBox)
  if(facet){
    Plot <- Plot + aes(group=Sample) +
      metadataOpts$facetBy +
      metadataOpts$colour +
      metadataOpts$lineType
  }else{
    Plot <- Plot + aes(group=Sample) +
      metadataOpts$colour +
      metadataOpts$lineType     
  }
  
  return(Plot)
})

setMethod("plotCC", "list", function(object,method="Coverage",facet=TRUE,
                                     facetBy=c("Tissue","Factor"),                                                 
                                     colourBy="Replicate",
                                     lineBy=NULL,
                                     addMetaData=NULL,excludedBox=FALSE
){
  
  ccvector <- crosscoverage(object)
  readlen <- max(readlength(object))
  shiftlength <- nrow(ccvector)
  
  toMelt <- data.frame("Shift_Size"=seq(1,shiftlength),
                       #"metadataOfInterest"=metadataOfInterest,
                       ccvector)
  CCDataFrame <- melt(toMelt,id.vars=c("Shift_Size"))
  
  
  metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy,lineBy)        
  
  
  CCDataFrameWithMetaData <- merge(CCDataFrame,metadataOpts$metadata,by.x=2,by.y=1,all=FALSE)
  
  colnames(CCDataFrameWithMetaData)[1:3] <- c("Sample","Shift_Size","CC_Score")
  
  Plot <- makeCCplot(CCDataFrameWithMetaData,shiftlength,readlen,excludedBox)
  if(facet){
    Plot <- Plot + aes(group=Sample) +
      metadataOpts$facetBy +
      metadataOpts$colour +
      metadataOpts$lineType
  }else{
    Plot <- Plot + aes(group=Sample) +
      metadataOpts$colour +
      metadataOpts$lineType      
  }
  return(Plot)
})
setMethod("plotCC", "ChIPQCsample", function(object,method="Coverage",excludedBox=FALSE){
  if(method=="Coverage"){
    ccvector <- crosscoverage(object)
    readlen <- readlength(object)
    shiftlength <- length(ccvector)
    CCDataFrame <- data.frame(cbind(seq(1,shiftlength),as.numeric(ccvector)))
    colnames(CCDataFrame) <- c("Shift_Size","CC_Score")
    Plot <- makeCCplot(CCDataFrame,shiftlength,readlen,excludedBox)
    return(Plot)
  }
})


#################################################################
#################################################################

## Peak Profile Plotting
makePeakProfilePlot <- function(PSDataFrame,Window){
  P <- ggplot(PSDataFrame,aes(x=Distance,y=Signal))+geom_line(size=1.3)+xlim(-Window/2,Window/2)+ylab("Signal")+
    theme(axis.title.y=element_text(angle=0))
  return(P)   
  
  
}

setGeneric("plotPeakProfile", function(object="ChIPQCexperiment",facet=TRUE,
                                       facetBy=c("Tissue","Factor"),
                                       colourBy="Replicate",
                                       lineBy=NULL,
                                       addMetaData=NULL
)                                      
  standardGeneric("plotPeakProfile"))
setMethod("plotPeakProfile", "ChIPQCexperiment", function(object,facet=TRUE,
                                                          facetBy=c("Tissue","Factor"),
                                                          colourBy="Replicate",
                                                          lineBy=NULL,
                                                          addMetaData=NULL                                                         
)
{
  PeakSignal <- averagepeaksignal(object)
  if(all(!is.na(PeakSignal))){
    Width <- seq(-nrow(PeakSignal)/2,nrow(PeakSignal)/2)[-(nrow(PeakSignal)/2+1)]
    PSDataFrame <- data.frame(Width,PeakSignal)
    
    PSDataFrame <- melt(PSDataFrame,id.vars=c("Width"))
    
    metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy,lineBy)        
    
    PSDataFrameWithMetaData <- merge(PSDataFrame,metadataOpts$metadata,by.x=2,by.y=1,all=FALSE)      
    colnames(PSDataFrameWithMetaData)[1:3] <- c("Sample","Distance","Signal")
    
    
    
    Plot <- makePeakProfilePlot(PSDataFrameWithMetaData,Window=length(Width))
    if(facet){
      Plot <- Plot + aes(group=Sample) +
        metadataOpts$facetBy +
        metadataOpts$colour +
        metadataOpts$lineType
    }else{
      Plot <- Plot + aes(group=Sample) +
        metadataOpts$colour +
        metadataOpts$lineType       
    }
    return(Plot)
  }else{
    stop("No average signal available for sample")
  }
})


setMethod("plotPeakProfile", "list", function(object,facet=TRUE,
                                              facetBy=c("Sample"),
                                              colourBy="Sample",
                                              lineBy=NULL,
                                              addMetaData=NULL                                                         
)
{
  PeakSignal <- averagepeaksignal(object)
  if(all(!is.na(PeakSignal))){
    Width <- seq(-nrow(PeakSignal)/2,nrow(PeakSignal)/2)[-(nrow(PeakSignal)/2+1)]
    PSDataFrame <- data.frame(Width,PeakSignal)
    
    PSDataFrame <- melt(PSDataFrame,id.vars=c("Width"))
    
    metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy,lineBy)        
    
    PSDataFrameWithMetaData <- merge(PSDataFrame,metadataOpts$metadata,by.x=2,by.y=1,all=FALSE)      
    colnames(PSDataFrameWithMetaData)[1:3] <- c("Sample","Distance","Signal")
    
    
    
    Plot <- makePeakProfilePlot(PSDataFrameWithMetaData,Window=length(Width))
    
    if(facet){
      Plot <- Plot + aes(group=Sample) +
        metadataOpts$facetBy +
        metadataOpts$colour +
        metadataOpts$lineType
    }else{
      Plot <- Plot + aes(group=Sample) +
        metadataOpts$colour +
        metadataOpts$lineType
    }
    return(Plot)
  }else{
    stop("No average signal available for sample")
  }
})


setMethod("plotPeakProfile", "ChIPQCsample", function(object){
  PeakSignal <- averagepeaksignal(object)
  if(all(!is.na(PeakSignal))){
    Width <- seq(-length(PeakSignal)/2,length(PeakSignal)/2)[-(length(PeakSignal)/2+1)]
    PSDataFrame <- data.frame(Width,PeakSignal)
    colnames(PSDataFrame) <- c("Distance","Signal")
    Plot <- makePeakProfilePlot(PSDataFrame,Window=length(Width))
    return(Plot)
  }else{
    stop("No average signal available for sample")
  }
})

#########################################
#########################################



makeCoveragePlot <- function(covHistFrame,maxDepthToPlot)
{
  covHistFrame <- covHistFrame[order(covHistFrame$Depth),]
  P <- ggplot(covHistFrame,aes(x=Depth,y=log10_bp,group=Sample))
  P <- P+geom_path(na.rm=TRUE,size=1.2)
  P <- P+xlim(0,maxDepthToPlot)
  P <- P+ylab("log10 BP")
  P <- P+theme(axis.title.y=element_text(angle=0))
  return(P)
}



setGeneric("plotCoverageHist", function(object="ChIPQCexperiment",maxDepthToPlot=100,facet=TRUE,
                                        facetBy=c("Tissue","Factor"),
                                        colourBy="Replicate",
                                        lineBy=NULL,
                                        addMetaData=NULL                                                         
                                        
)
  standardGeneric("plotCoverageHist"))
setMethod("plotCoverageHist", "ChIPQCexperiment", function(object,maxDepthToPlot=100,facet=TRUE,
                                                           facetBy=c("Tissue","Factor"),
                                                           colourBy="Replicate",
                                                           lineBy=NULL,
                                                           addMetaData=NULL                                                         
)
{
  
  
  covHist <- melt(coveragehistogram(object))
  
  covHistFrame <- data.frame(
    Sample=covHist$Var2,
    Depth = as.numeric(as.vector(covHist$Var1)),
    log10_bp=log10(as.numeric(as.vector(covHist$value)))
  )   
  
  metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy,lineBy)        
  
  covHistFrameWithMetaData <- merge(covHistFrame,metadataOpts$metadata,by.x=1,by.y=1,all=FALSE)      
  
  colnames(covHistFrameWithMetaData)[1:3] <- c("Sample","Depth","log10_bp")
  
  
  
  Plot <- makeCoveragePlot(covHistFrameWithMetaData,maxDepthToPlot)
  if(facet){
    Plot <- Plot + 
      metadataOpts$facetBy +
      metadataOpts$colour +
      metadataOpts$lineType
  }else{
    Plot <- Plot +
      metadataOpts$colour +
      metadataOpts$lineType     
  }
  
  return(Plot)
}
)

setMethod("plotCoverageHist", "list", function(object,maxDepthToPlot=100,facet=TRUE,
                                               facetBy=c("Sample"),
                                               colourBy="Sample",
                                               lineBy=NULL,
                                               addMetaData=NULL                                                         
)
{
  
  
  covHist <- melt(coveragehistogram(object))
  
  covHistFrame <- data.frame(
    Sample=covHist$Var2,
    Depth = as.numeric(as.vector(covHist$Var1)),
    log10_bp=log10(as.numeric(as.vector(covHist$value)))
  )   
  
  metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy,lineBy)        
  
  covHistFrameWithMetaData <- merge(covHistFrame,metadataOpts$metadata,by.x=1,by.y=1,all=FALSE)      
  
  colnames(covHistFrameWithMetaData)[1:3] <- c("Sample","Depth","log10_bp")
  
  
  
  Plot <- makeCoveragePlot(covHistFrameWithMetaData,maxDepthToPlot)
  if(facet){
    Plot <- Plot + 
      metadataOpts$facetBy +
      metadataOpts$colour +
      metadataOpts$lineType
  }else{
    Plot <- Plot + 
      metadataOpts$colour +
      metadataOpts$lineType
    
  }
  
  return(Plot)
}
)

setMethod("plotCoverageHist", "ChIPQCsample", function(object,maxDepthToPlot=100){
  
  
  
  covHist <- coveragehistogram(object)
  covHistFrame <- data.frame(Sample="Sample",Depth=seq(1,length(covHist)),log10_bp=log10(covHist))
  Plot <- makeCoveragePlot(covHistFrame,maxDepthToPlot)  
  return(Plot)
}
)



########################################################
########################################################
makeFripPlot <- function(fripDataFrame){
  fripDataFrame$Reads <- factor(as.vector(fripDataFrame$Reads),levels = c("OutSide","Inside"))
  P <- ggplot(fripDataFrame, aes(Sample,FRIP, fill=Reads))
  P <- P+geom_bar(stat="identity")
  P <- P+labs(title = "Percentage of Reads In Peaks")
  P <- P+theme(axis.title.y=element_text(angle=0),panel.background = element_blank())+scale_fill_manual(values=c( "#CCCCFF","#000099"))
  
  return(P)
}

setGeneric("plotFrip", function(object="ChIPQCexperiment",type="barstacked",facet=TRUE,
                                facetBy=c("Tissue","Factor"),
                                addMetaData=NULL,AsPercent=TRUE) standardGeneric("plotFrip"))



setMethod("plotFrip", "ChIPQCexperiment", function(object,type="barstacked",facet=TRUE,
                                                   facetBy=c("Tissue","Factor"),
                                                   addMetaData=NULL,AsPercent=TRUE){
  rip <- rip(object)
  mapped <- mapped(object)  
  ripWithPeaks <- rip[!is.na(rip)]
  mappedWithPeaks <- mapped[!is.na(rip)]
  toMelt <- data.frame(Sample=names(ripWithPeaks),Inside=ripWithPeaks,OutSide=mappedWithPeaks-ripWithPeaks)
  PercentInside <- toMelt[,"Inside"]/rowSums(toMelt[,c("Inside","OutSide")])
  PercentOutside <- toMelt[,"OutSide"]/rowSums(toMelt[,c("Inside","OutSide")])
  if(AsPercent){
    toMelt[,"Inside"] <- PercentInside*100
    toMelt[,"OutSide"] <- PercentOutside*100
  }
  
  
  fripDataFrame <- melt(toMelt)
  
  
  metadataOpts <- ChIPQC:::mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)        
  
  
  fripDataFrameWithMetaData <- merge(fripDataFrame,metadataOpts$metadata,by.x=1,by.y=1,all=FALSE)      
  colnames(fripDataFrameWithMetaData)[1:3] <- c("Sample","Reads","FRIP")
  fripDataFrameWithMetaData <- fripDataFrameWithMetaData[order(fripDataFrameWithMetaData[,"Sample"],fripDataFrameWithMetaData[,"Reads"]),]
  Plot <- ChIPQC:::makeFripPlot(fripDataFrameWithMetaData)
  if(facet){
    metadataOpts$facetBy$params$free$x <- TRUE
    Plot <- Plot + 
      metadataOpts$facetBy 
  }
  
  return(Plot)
})


setMethod("plotFrip", "list", function(object,type="barstacked",facet=TRUE,
                                       facetBy=c("Sample"),
                                       addMetaData=NULL,AsPercent=TRUE){
  rip <- rip(object)
  mapped <- mapped(object)  
  ripWithPeaks <- rip[!is.na(rip)]
  mappedWithPeaks <- mapped[!is.na(rip)]
  toMelt <- data.frame(Sample=names(ripWithPeaks),Inside=ripWithPeaks,OutSide=mappedWithPeaks-ripWithPeaks)
  PercentInside <- toMelt[,"Inside"]/rowSums(toMelt[,c("Inside","OutSide")])
  PercentOutside <- toMelt[,"OutSide"]/rowSums(toMelt[,c("Inside","OutSide")])
  if(AsPercent){
    toMelt[,"Inside"] <- PercentInside*100
    toMelt[,"OutSide"] <- PercentOutside*100
  }
  
  
  fripDataFrame <- melt(toMelt)
  
  
  metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)        
  
  
  fripDataFrameWithMetaData <- merge(fripDataFrame,metadataOpts$metadata,by.x=1,by.y=1,all=FALSE)      
  colnames(fripDataFrameWithMetaData)[1:3] <- c("Sample","Reads","FRIP")
  fripDataFrameWithMetaData <- fripDataFrameWithMetaData[order(fripDataFrameWithMetaData[,"Sample"],fripDataFrameWithMetaData[,"Reads"]),]
  Plot <- makeFripPlot(fripDataFrameWithMetaData)
  if(facet){
    metadataOpts$facetBy$params$free$x <- TRUE
    Plot <- Plot + 
      metadataOpts$facetBy 
  }
  return(Plot)
})

setMethod("plotFrip", "ChIPQCsample", function(object,type="barstacked",facet=TRUE,
                                               facetBy=c("Tissue","Factor"),AsPercent=TRUE){
  rip <- rip(object)
  mapped <- mapped(object)  
  fripDataFrame <- data.frame(
    Sample=c("Sample","Sample"),
    Reads=c("Inside","Outside"),
    FRIP=c(rip,mapped-rip))
  
  if(AsPercent){
    fripDataFrame[,"FRIP"] <- (fripDataFrame[,"FRIP"]/sum(fripDataFrame[,"FRIP"]))*100
  }
  
  
  Plot <- makeFripPlot(fripDataFrame)
  Plot <- Plot+xlab("")
  return(Plot)
})

###############
###############
makeFriblPlot <- function(friblDataFrame){
  friblDataFrame$Reads <- factor(as.vector(friblDataFrame$Reads),levels = c("OutSide","Inside"))
  
  P <- ggplot(friblDataFrame, aes(Sample,FRIBL, fill=Reads))
  P <- P+geom_bar(stat="identity")
  P <- P+labs(title = "Percentage Of Reads In Blacklists")
  P <- P+theme(axis.title.y=element_text(angle=0),panel.background = element_blank())+scale_fill_manual(values=c( "#CCCCFF","#000099"))      
  return(P)
}



setGeneric("plotFribl", function(object="ChIPQCexperiment",type="barstacked",facet=TRUE,
                                 facetBy=c("Tissue","Factor"),
                                 addMetaData=NULL,AsPercent=TRUE) standardGeneric("plotFribl"))

setMethod("plotFribl", "ChIPQCexperiment", function(object,type="barstacked",facet=TRUE,
                                                    facetBy=c("Tissue","Factor"),
                                                    addMetaData=NULL,AsPercent=TRUE){
  AsPercent <- TRUE
  ribl <- ribl(object)
  mapped <- mapped(object)  
  riblWithBLs <- ribl[!is.na(ribl)]
  mappedWithBLs <- mapped[!is.na(ribl)]
  toMelt <- data.frame(Sample=names(riblWithBLs),Inside=riblWithBLs,OutSide=mappedWithBLs-riblWithBLs)
  PercentInside <- toMelt[,"Inside"]/rowSums(toMelt[,c("Inside","OutSide")])
  PercentOutside <- toMelt[,"OutSide"]/rowSums(toMelt[,c("Inside","OutSide")])
  if(AsPercent){
    toMelt[,"Inside"] <- PercentInside*100
    toMelt[,"OutSide"] <- PercentOutside*100
  }
  friblDataFrame <- melt(toMelt)
  
  
  metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)        
  
  
  friblDataFrameWithMetaData <- merge(friblDataFrame,metadataOpts$metadata,by.x=1,by.y=1,all=FALSE)      
  colnames(friblDataFrameWithMetaData)[1:3] <- c("Sample","Reads","FRIBL")
  friblDataFrameWithMetaData <- friblDataFrameWithMetaData[order(friblDataFrameWithMetaData[,"Sample"],friblDataFrameWithMetaData[,"Reads"]),]
  Plot <- makeFriblPlot(friblDataFrameWithMetaData)
  if(facet){
    metadataOpts$facetBy$params$free$x <- TRUE
    Plot <- Plot + 
      metadataOpts$facetBy 
  }
  #   if(AsPercent){
  #      Plot+scale_y_continuous(labels = "percent")
  #   }
  return(Plot)
  
})


setMethod("plotFribl", "list", function(object,type="barstacked",facet=TRUE,
                                        facetBy=c("Sample"),
                                        addMetaData=NULL,AsPercent=TRUE){
  AsPercent <- TRUE
  ribl <- ribl(object)
  mapped <- mapped(object)  
  riblWithBLs <- ribl[!is.na(ribl)]
  mappedWithBLs <- mapped[!is.na(ribl)]
  toMelt <- data.frame(Sample=names(riblWithBLs),Inside=riblWithBLs,OutSide=mappedWithBLs-riblWithBLs)
  PercentInside <- toMelt[,"Inside"]/rowSums(toMelt[,c("Inside","OutSide")])
  PercentOutside <- toMelt[,"OutSide"]/rowSums(toMelt[,c("Inside","OutSide")])
  if(AsPercent){
    toMelt[,"Inside"] <- PercentInside*100
    toMelt[,"OutSide"] <- PercentOutside*100
  }
  friblDataFrame <- melt(toMelt)
  
  
  metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)        
  
  
  friblDataFrameWithMetaData <- merge(friblDataFrame,metadataOpts$metadata,by.x=1,by.y=1,all=FALSE)      
  colnames(friblDataFrameWithMetaData)[1:3] <- c("Sample","Reads","FRIBL")
  friblDataFrameWithMetaData <- friblDataFrameWithMetaData[order(friblDataFrameWithMetaData[,"Sample"],friblDataFrameWithMetaData[,"Reads"]),]
  Plot <- makeFriblPlot(friblDataFrameWithMetaData)
  if(facet){
    metadataOpts$facetBy$params$free$x <- TRUE
    Plot <- Plot + 
      metadataOpts$facetBy 
  }
  #   if(AsPercent){
  #      Plot+scale_y_continuous(labels = "percent")
  #   }
  return(Plot)
  
})

setMethod("plotFribl", "ChIPQCsample", function(object,type="barstacked",AsPercent=TRUE){
  #AsPercent <- TRUE
  ribl <- ribl(object)
  mapped <- mapped(object) 
  friblDataFrame <- data.frame(
    Sample=c("Sample","Sample"),
    Reads=c("Inside","Outside"),
    FRIBL=c(ribl,mapped-ribl)   
  )
  if(AsPercent){
    friblDataFrame[,"FRIBL"] <- (friblDataFrame[,"FRIBL"]/sum(friblDataFrame[,"FRIBL"]))*100
  }
  Plot <- makeFriblPlot(friblDataFrame)
  Plot <- Plot+xlab("")
  return(Plot)
  
})

##############################################
##############################################
makeRapPlot <- function(rapDataFrame){
  
  P <- ggplot(rapDataFrame, aes(Sample, CountsInPeaks))+
    geom_boxplot(fill="lightblue")+
    theme(axis.title.y=element_text(angle=0))+
    xlab("")
  return(P)
}

setGeneric("plotRap", function(object="ChIPQCexperiment",facet=TRUE,
                               facetBy=c("Tissue","Factor"),
                               addMetaData=NULL) standardGeneric("plotRap"))


setMethod("plotRap", "ChIPQCsample", function(object){
  
  Peaks <- peaks(object)
  if(length(Peaks) > 0){
    #CountsInPeaks <- (Peaks$Counts/mapped(object))*10^6  
    CountsInPeaks <- Peaks$Counts
    rapDataFrame <- data.frame(Sample="Sample",CountsInPeaks=CountsInPeaks)
    Plot <- makeRapPlot(rapDataFrame)
    return(Plot)
  }
  
})

setMethod("plotRap", "ChIPQCexperiment", function(object,facet=TRUE,
                                                  facetBy=c("Tissue","Factor"),
                                                  addMetaData=NULL){
  
  Peaks <- peaks(object)
  haspeaks = lapply(Peaks,length)>0
  Peaks = Peaks[haspeaks]
  peakCountList <- sapply(Peaks,function(x)elementMetadata(x)$Counts)
  peakCountMatrix <- list2matrix(peakCountList)
  ylimNew <- c(0,max(unlist(lapply(apply(peakCountMatrix,2,function(x) boxplot.stats(x[x!=0])),
                                   function(x) x$stats[5]))))
  
  rapDataFrame <- melt(peakCountMatrix)
  metadata <- QCmetadata(object)[haspeaks,]         
  
  metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)        
  
  
  rapDataFrameWithMetaData <- merge(rapDataFrame,metadataOpts$metadata,by.x=2,by.y=1,all=FALSE)      
  colnames(rapDataFrameWithMetaData)[1:3] <- c("Sample","PeakNumber","CountsInPeaks")
  
  
  Plot <- makeRapPlot(rapDataFrameWithMetaData)
  
  metadataOpts$facetBy$params[2]$free$x <- TRUE
  if(facet){
    Plot <- Plot + metadataOpts$facetBy 
  }
  Plot <- Plot + coord_cartesian(ylim = ylimNew*1.05)
  return(Plot)
  
})


setMethod("plotRap", "list", function(object,facet=TRUE,
                                      facetBy=c("Sample"),
                                      addMetaData=NULL){
  
  Peaks <- peaks(object)
  haspeaks = lapply(Peaks,length)>0
  Peaks = Peaks[haspeaks]
  peakCountList <- sapply(Peaks,function(x)elementMetadata(x)$Counts)
  peakCountMatrix <- list2matrix(peakCountList)
  ylimNew <- c(0,max(unlist(lapply(apply(peakCountMatrix,2,function(x) boxplot.stats(x[x!=0])),
                                   function(x) x$stats[5]))))
  
  rapDataFrame <- melt(peakCountMatrix)
  metadata <- QCmetadata(object)[haspeaks,]         
  
  metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)        
  
  
  rapDataFrameWithMetaData <- merge(rapDataFrame,metadataOpts$metadata,by.x=2,by.y=1,all=FALSE)      
  colnames(rapDataFrameWithMetaData)[1:3] <- c("Sample","PeakNumber","CountsInPeaks")
  
  
  Plot <- makeRapPlot(rapDataFrameWithMetaData)
  
  metadataOpts$facetBy$params[2]$free$x <- TRUE
  if(facet){
    Plot <- Plot + metadataOpts$facetBy 
  }
  Plot <- Plot + coord_cartesian(ylim = ylimNew*1.05)
  return(Plot)
  
})

#########################################
#########################################

makeRegiPlot <- function(regiScoresFrame){
  regiScoresFrame[,"GenomicIntervals"] <- factor(regiScoresFrame[,"GenomicIntervals"],levels=unique(as.vector(regiScoresFrame[,"GenomicIntervals"])))
  Plot <- ggplot(regiScoresFrame, aes(Sample,GenomicIntervals))  
  Plot <- Plot+geom_tile(aes(y=Sample,x=GenomicIntervals,fill = log2_Enrichment)) +
    scale_fill_gradient2(low="blue",high="yellow",mid="black",midpoint=median(regiScoresFrame$log2_Enrichment))
  return(Plot)
}

setGeneric("plotRegi", function(object="ChIPQCexperiment",facet=TRUE,
                                facetBy=c("Tissue","Factor"),
                                addMetaData=NULL) standardGeneric("plotRegi"))

setMethod("plotRegi", "ChIPQCexperiment", function(object,facet=TRUE,
                                                   facetBy=c("Tissue","Factor"),
                                                   addMetaData=NULL){
  regiScores <- regi(object)
  hasna = apply(is.na(regiScores),2,sum)>0
  if(sum(hasna)==length(hasna)) {
    warning('No genomic annotation computed')
    invisible(NULL)
  } else {
    regiScores = regiScores[,!hasna]
  }
  #rownames(regiScores)[5:7] <- c("500_Upstream","2000To500_Upstream","20000To2000_Upstream")   
  meltedDF <- melt(regiScores)
  
  #metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)        
  regiScoresFrame <- data.frame(Sample=meltedDF$Var2,
                                GenomicIntervals=meltedDF$Var1,
                                log2_Enrichment=meltedDF$value)
  metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)
  metadataOpts$metadata = metadataOpts$metadata[!hasna,]
  
  facetGridForm <- as.formula(paste0(paste(names(metadataOpts$facetBy$params$facets),collapse="+"),"~."))
  
  regiScoresFrameWithMetadata <- merge(regiScoresFrame,metadataOpts$metadata,by.x=1,by.y=1,alQl=FALSE,sort=FALSE)      
  
  #regiScoresFrameWithMetadata[,"GenomicIntervals"] <- ordered(regiScoresFrameWithMetadata[,"GenomicIntervals"],c("20000To2000_Upstream","2000To500_Upstream",
  #                                                                                                               "500_Upstream","5UTRs","Transcripts","Introns","3UTRs"))
  
  Plot <- makeRegiPlot(regiScoresFrameWithMetadata)
  if(facet){
    Plot <- Plot+facet_grid(facetGridForm,scales="free_y")#,ncol=1)
  }
  Plot <- Plot + theme(text = element_text(size=20),
                       axis.text.x = element_text(angle=90, vjust=1))+xlab("")+ylab("")
  
  return(Plot)
}
)




setMethod("plotRegi", "list", function(object,facet=TRUE,
                                       facetBy=c("Sample"),
                                       addMetaData=NULL){
  if(all(sapply(object,function(x)class(x) == "ChIPQCsample")) & all(sapply(names(object),function(x)!is.null(x)))){
    regiScores <- regi(object)
    hasna = apply(is.na(regiScores),2,sum)>0
    if(sum(hasna)==length(hasna)) {
      warning('No genomic annotation computed')
      invisible(NULL)
    } else {
      regiScores = regiScores[,!hasna]
    }
    #rownames(regiScores)[5:7] <- c("500_Upstream","2000To500_Upstream","20000To2000_Upstream")   
    meltedDF <- melt(regiScores)
    
    #metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)        
    regiScoresFrame <- data.frame(Sample=meltedDF$Var2,
                                  GenomicIntervals=meltedDF$Var1,
                                  log2_Enrichment=meltedDF$value)
    metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)
    metadataOpts$metadata = metadataOpts$metadata[!hasna,]
    
    facetGridForm <- as.formula(paste0(paste(names(metadataOpts$facetBy$params$facets),collapse="+"),"~."))
    
    regiScoresFrameWithMetadata <- merge(regiScoresFrame,metadataOpts$metadata,by.x=1,by.y=1,all=FALSE,sort=FALSE)      
    
    #regiScoresFrameWithMetadata[,"GenomicIntervals"] <- ordered(regiScoresFrameWithMetadata[,"GenomicIntervals"],c("20000To2000_Upstream","2000To500_Upstream",
    #                                                                                                               "500_Upstream","5UTRs","Transcripts","Introns","3UTRs"))
    
    Plot <- makeRegiPlot(regiScoresFrameWithMetadata)
    if(facet){
      Plot <- Plot+facet_grid(facetGridForm,scales="free_y")#,ncol=1)
    }
    Plot <- Plot + theme(text = element_text(size=20),
                         axis.text.x = element_text(angle=90, vjust=1))+xlab("")+ylab("")
    
    return(Plot)
  }else{
    NULL
  }
}
)


setMethod("plotRegi", "ChIPQCsample", function(object){
  ### PlaceHolder!
  regiScores <- regi(object)
  if(sum(is.na(regiScores))>0) {
    warning('No genomic annotation computed')
    invisible(NULL)
  }
  regiScoresFrame <- data.frame(Sample="SampleName",GenomicIntervals=names(regiScores),log2_Enrichment=unname(regiScores))
  Plot <- makeRegiPlot(regiScoresFrame)
  
  return(Plot)
}
)

#################################################################
#################################################################

#################################################################
#################################################################
setGeneric("plotSSD", function(object="ChIPQCexperiment",facet=TRUE,
                               facetBy=c("Tissue","Factor"),addMetaData=NULL)
  standardGeneric("plotSSD"))

makeSSDPlot <- function(SSDDataFrame){
  Plot <- ggplot(SSDDataFrame, aes(y = Sample,x=SSD,colour=Filter)) + 
    geom_point(size=5)+xlim(0.1,4)
  return(Plot)
}


setMethod("plotSSD", "ChIPQCsample", function(object,facet=TRUE,
                                              facetBy=c("Tissue","Factor")){
  ssd <- object@SSD
  ssdBL <- object@SSDBL
  #ssdBL[is.na(ssdBL)] <- ssd[is.na(ssdBL)] 
  SSDDataFrame <- data.frame(
    Sample=c("Sample","Sample"),
    Filter=c("Pre_Blacklist","Post_Blacklist"),
    SSD=c(ssd,ssdBL))
  
  Plot <- makeSSDPlot(SSDDataFrame)
  return(Plot)
})

setMethod("plotSSD", "ChIPQCexperiment", function(object,facet=TRUE,
                                                  facetBy=c("Tissue","Factor"),addMetaData=NULL){
  ssd <- unlist(lapply(QCsample(object),function(x) x@SSD))
  ssdBL <- unlist(lapply(QCsample(object),function(x) x@SSDBL))
  #ssdBL[is.na(ssdBL)] <- ssd[is.na(ssdBL)] 
  SSDDataFrame <- data.frame(
    Sample=names(ssdBL),
    Pre_Blacklist=ssd,
    Post_Blacklist=ssdBL)
  meltedSSDDataFrame <- melt(SSDDataFrame)
  metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)
  
  SSDDataFrameWithMetaData <- merge(meltedSSDDataFrame,metadataOpts$metadata,by.x=1,by.y=1,all=FALSE)
  
  facetGridForm <- as.formula(paste0(paste(names(metadataOpts$facetBy$params$facets),collapse="+"),"~."))
  
  colnames(SSDDataFrameWithMetaData)[1:3] <- c("Sample","Filter","SSD")
  
  Plot <- makeSSDPlot(SSDDataFrameWithMetaData)
  if(facet){
    Plot <- Plot +facet_grid(facetGridForm,scales="free_y")
  }
  return(Plot)
})


setMethod("plotSSD", "list", function(object,facet=TRUE,
                                      facetBy=c("Sample"),addMetaData=NULL){
  ssd <- unlist(lapply(object,function(x) x@SSD))
  ssdBL <- unlist(lapply(object,function(x) x@SSDBL))
  #ssdBL[is.na(ssdBL)] <- ssd[is.na(ssdBL)] 
  SSDDataFrame <- data.frame(
    Sample=names(ssdBL),
    Pre_Blacklist=ssd,
    Post_Blacklist=ssdBL)
  meltedSSDDataFrame <- melt(SSDDataFrame)
  metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)
  
  SSDDataFrameWithMetaData <- merge(meltedSSDDataFrame,metadataOpts$metadata,by.x=1,by.y=1,all=FALSE)
  
  facetGridForm <- as.formula(paste0(paste(names(metadataOpts$facetBy$params$facets),collapse="+"),"~."))
  
  colnames(SSDDataFrameWithMetaData)[1:3] <- c("Sample","Filter","SSD")
  
  Plot <- makeSSDPlot(SSDDataFrameWithMetaData)
  if(facet){
    Plot <- Plot +facet_grid(facetGridForm,scales="free_y")
  }
  return(Plot)
  
})


#################################################################
#################################################################

