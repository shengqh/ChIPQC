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
   if(all(colourBy %in% colnames(metadata))){
      colour <- aes_string(colour=colourBy)
   }else{
      colour <- NULL  
   }
   if(all(lineBy %in% colnames(metadata))){
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

makeCCplot <- function(CCDataFrame,shiftlength,readlen){
   P <- ggplot(CCDataFrame,aes(x=Shift_Size,y=CC_Score))+geom_path(alpha = 1,size=1.3)+xlim(0,shiftlength)+ylab("CC_Score")+
      geom_rect(data=CCDataFrame,xmin=-Inf,colour="red",xmax=readlen*2,ymin=-Inf,ymax=Inf,fill="red",alpha = 0.002)+
      theme(axis.title.y=element_text(angle=0))
   return(P)
}

#setGeneric("plotCC", function(object="ChIPQCsample",method="Coverage") standardGeneric("plotCC"))

setGeneric("plotCC", function(object="ChIPQCexperiment",method="Coverage",facet=TRUE,
                              facetBy=c("Tissue","Factor"),
                              colourBy="Replicate",
                              lineBy=NULL,
                              addMetaData=NULL
)
   standardGeneric("plotCC")
)

setMethod("plotCC", "ChIPQCexperiment", function(object,method="Coverage",facet=TRUE,
                                                 facetBy=c("Tissue","Factor"),                                                 
                                                 colourBy="Replicate",
                                                 lineBy=NULL,
                                                 addMetaData=NULL
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
   
   Plot <- makeCCplot(CCDataFrameWithMetaData,shiftlength,readlen)      
   Plot <- Plot + aes(group=Sample) +
      metadataOpts$facetBy +
      metadataOpts$colour +
      metadataOpts$lineType
   
   return(Plot)
})

setMethod("plotCC", "ChIPQCsample", function(object,method="Coverage"){
   if(method=="Coverage"){
      ccvector <- crosscoverage(object)
      readlen <- readlength(object)
      shiftlength <- length(ccvector)
      CCDataFrame <- data.frame(cbind(seq(1,shiftlength),as.numeric(ccvector)))
      colnames(CCDataFrame) <- c("Shift_Size","CC_Score")
      Plot <- makeCCplot(CCDataFrame,shiftlength,readlen)
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
      Plot <- Plot + aes(group=Sample) +
         metadataOpts$facetBy +
         metadataOpts$colour +
         metadataOpts$lineType
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
   Plot <- Plot + 
      metadataOpts$facetBy +
      metadataOpts$colour +
      metadataOpts$lineType
   
   
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
   
   P <- ggplot(fripDataFrame, aes(Sample,FRIP, fill=Reads))
   P <- P+geom_bar(stat="identity")
   P <- P+labs(title = "Reads In Peaks")
   P <- P+theme(axis.title.y=element_text(angle=0),panel.background = element_blank())+scale_fill_manual(values=c( "#000099","#CCCCFF"))
   
   return(P)
}

setGeneric("plotFrip", function(object="ChIPQCexperiment",type="barstacked",facet=TRUE,
                                facetBy=c("Tissue","Factor"),
                                addMetaData=NULL) standardGeneric("plotFrip"))



setMethod("plotFrip", "ChIPQCexperiment", function(object,type="barstacked",facet=TRUE,
                                                   facetBy=c("Tissue","Factor"),
                                                   addMetaData=NULL){
   rip <- rip(object)
   mapped <- mapped(object)  
   ripWithPeaks <- rip[!is.na(rip)]
   mappedWithPeaks <- mapped[!is.na(rip)]
   toMelt <- data.frame(Sample=names(ripWithPeaks),Inside=ripWithPeaks,OutSide=mappedWithPeaks-ripWithPeaks)
   fripDataFrame <- melt(toMelt)
   
   
   metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)        
   
   
   fripDataFrameWithMetaData <- merge(fripDataFrame,metadataOpts$metadata,by.x=1,by.y=1,all=FALSE)      
   colnames(fripDataFrameWithMetaData)[1:3] <- c("Sample","Reads","FRIP")
   fripDataFrameWithMetaData <- fripDataFrameWithMetaData[order(fripDataFrameWithMetaData[,"Sample"],fripDataFrameWithMetaData[,"Reads"]),]
   Plot <- makeFripPlot(fripDataFrameWithMetaData)
   metadataOpts$facetBy$free$x <- TRUE
   Plot <- Plot + 
      metadataOpts$facetBy 
   
   return(Plot)
})

setMethod("plotFrip", "ChIPQCsample", function(object,type="barstacked",facet=TRUE,
                                               facetBy=c("Tissue","Factor")){
   rip <- rip(object)
   mapped <- mapped(object)  
   fripDataFrame <- data.frame(
      Sample=c("Sample","Sample"),
      Reads=c("Inside Peaks","Outside Peaks"),
      FRIP=c(rip,mapped-rip))
   Plot <- makeFripPlot(fripDataFrame)
   Plot <- Plot+xlab("")
   return(Plot)
})

###############
###############
makeFriblPlot <- function(friblDataFrame){
   
   P <- ggplot(friblDataFrame, aes(Sample,FRIBL, fill=Reads))
   P <- P+geom_bar(stat="identity")
   P <- P+labs(title = "Proportion Of Reads In Blacklists")
   P <- P+theme(axis.title.y=element_text(angle=0),panel.background = element_blank())+scale_fill_manual(values=c( "#000099","#CCCCFF"))      
   return(P)
}



setGeneric("plotFribl", function(object="ChIPQCexperiment",type="barstacked",facet=TRUE,
                                 facetBy=c("Tissue","Factor"),
                                 addMetaData=NULL) standardGeneric("plotFribl"))

setMethod("plotFribl", "ChIPQCexperiment", function(object,type="barstacked",facet=TRUE,
                                                    facetBy=c("Tissue","Factor"),
                                                    addMetaData=NULL){
   ribl <- ribl(object)
   mapped <- mapped(object)  
   riblWithBLs <- ribl[!is.na(ribl)]
   mappedWithBLs <- mapped[!is.na(ribl)]
   toMelt <- data.frame(Sample=names(riblWithBLs),Inside=riblWithBLs,OutSide=mappedWithBLs-riblWithBLs)
   friblDataFrame <- melt(toMelt)
   
   
   metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)        
   
   
   friblDataFrameWithMetaData <- merge(friblDataFrame,metadataOpts$metadata,by.x=1,by.y=1,all=FALSE)      
   colnames(friblDataFrameWithMetaData)[1:3] <- c("Sample","Reads","FRIBL")
   friblDataFrameWithMetaData <- friblDataFrameWithMetaData[order(friblDataFrameWithMetaData[,"Sample"],friblDataFrameWithMetaData[,"Reads"]),]
   Plot <- makeFriblPlot(friblDataFrameWithMetaData)
   metadataOpts$facetBy$free$x <- TRUE
   Plot <- Plot + 
      metadataOpts$facetBy 
   return(Plot)
   
})

setMethod("plotFribl", "ChIPQCsample", function(object,type="barstacked"){
   ribl <- ribl(object)
   mapped <- mapped(object) 
   friblDataFrame <- data.frame(
      Sample=c("Sample","Sample"),
      Reads=c("Inside Blacklist","Outside Blacklist"),
      FRIBL=c(ribl,mapped-ribl)   
   )  
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
   peakCountList <- sapply(peaks(object),function(x)elementMetadata(x)$Counts)
   peakCountMatrix <- list2matrix(peakCountList)
   #   rapDataFrame <- data.frame(
   
   
   #   )
   rapDataFrame <- melt(peakCountMatrix)
   metadata <- QCmetadata(object)         
   
   metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)        
   
   
   rapDataFrameWithMetaData <- merge(rapDataFrame,metadataOpts$metadata,by.x=2,by.y=1,all=FALSE)      
   colnames(rapDataFrameWithMetaData)[1:3] <- c("Sample","PeakNumber","CountsInPeaks")
   
   Plot <- makeRapPlot(rapDataFrameWithMetaData)
   
   metadataOpts$facetBy[2]$free$x <- TRUE
   
   Plot <- Plot + metadataOpts$facetBy 
   return(Plot)
   
})


#########################################
#########################################

makeRegiPlot <- function(regiScoresFrame){
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
   rownames(regiScores)[5:7] <- c("500_Upstream","2000To500_Upstream","20000To2000_Upstream")   
   meltedDF <- melt(regiScores)
   
   metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)        
   regiScoresFrame <- data.frame(Sample=meltedDF$Var2,
                                 GenomicIntervals=meltedDF$Var1,
                                 log2_Enrichment=meltedDF$value)
   metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy=NULL,lineBy=NULL)        
   facetGridForm <- as.formula(paste0(paste(names(metadataOpts$facetBy$facets),collapse="+"),"~."))
   
   regiScoresFrameWithMetadata <- merge(regiScoresFrame,metadataOpts$metadata,by.x=1,by.y=1,all=FALSE)      
   
   regiScoresFrameWithMetadata[,"GenomicIntervals"] <- ordered(regiScoresFrameWithMetadata[,"GenomicIntervals"],c("20000To2000_Upstream","2000To500_Upstream",
                                                                                                                  "500_Upstream","5UTRs","Transcripts","Introns","3UTRs"))
   
   Plot <- makeRegiPlot(regiScoresFrameWithMetadata)
   
   Plot <- Plot+facet_grid(facetGridForm,scales="free_y")#,ncol=1)
   Plot <- Plot + theme(text = element_text(size=20),
                        axis.text.x = element_text(angle=90, vjust=1))+xlab("")+ylab("")
   
   return(Plot)
}
)


setMethod("plotRegi", "ChIPQCsample", function(object){
   ### PlaceHolder!
   regiScores <- regi(object)
   regiScoresFrame <- data.frame(Sample="SampleName",GenomicIntervals=names(regiScores),log2_Enrichment=unname(regiScores))
   Plot <- makeRegiPlot(regiScoresFrame)
   
   return(Plot)
}
)

#################################################################
#################################################################
setGeneric("plotSSD", function(object="ChIPQCexperiment",method="Coverage",facet=TRUE,
                               facetBy=c("Tissue","Factor"),
                               colourBy="Replicate",
                               lineBy=NULL,
                               addMetaData=NULL
)
   standardGeneric("plotSSD")
)

setMethod("plotSSD", "ChIPQCexperiment", function(object,method="Coverage",facet=TRUE,
                                                  facetBy=c("Tissue","Factor"),                                                 
                                                  colourBy="Replicate",
                                                  lineBy=NULL,
                                                  addMetaData=NULL
){
   
   ssdvector <- ssd(object)
   for(sample in names(ssdvector)){
      input <- sample
   }
   
   
   
   toMelt <- data.frame("Shift_Size"=seq(1,shiftlength),
                        #"metadataOfInterest"=metadataOfInterest,
                        ccvector)
   CCDataFrame <- melt(toMelt,id.vars=c("Shift_Size"))
   
   
   metadataOpts <- mergeMetadata(object,addMetaData,facetBy,colourBy,lineBy)        
   
   
   CCDataFrameWithMetaData <- merge(CCDataFrame,metadataOpts$metadata,by.x=2,by.y=1,all=FALSE)
   
   colnames(CCDataFrameWithMetaData)[1:3] <- c("Sample","Shift_Size","CC_Score")
   
   Plot <- makeCCplot(CCDataFrameWithMetaData,shiftlength,readlen)      
   Plot <- Plot + 
      metadataOpts$facetBy +
      metadataOpts$colour +
      metadataOpts$lineType
   
   return(Plot)
})

#################################################################
#################################################################

