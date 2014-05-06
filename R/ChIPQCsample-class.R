
setClass("ChIPQCsample",contains = "GRanges",
         slots=c(AveragePeakSignal="list",
                 CrossCoverage="numeric",CrossCorrelation="numeric",SSD="numeric",SSDBL="numeric",CountsInPeaks="numeric",
                 CountsInBlackList="numeric",CountsInFeatures="list",PropInFeatures="list",
                 CoverageHistogram="numeric",FlagAndTagCounts="numeric",
                 readlength="numeric"
         ))


setMethod("show","ChIPQCsample",
          function (object){
             message("\t\t\t\t\t",class(object),"")
             message("Number of Mapped reads: ",object@FlagAndTagCounts[2],"")
             message("Number of Mapped reads passing MapQ filter: ",object@FlagAndTagCounts[4],"")    
             message("Percentage Of Reads as Non-Duplicates (NRF): ",round((((object@FlagAndTagCounts[2])-(object@FlagAndTagCounts[5]))/object@FlagAndTagCounts[2])*100, digits=2),"(",round(object@FlagAndTagCounts[5]/object@FlagAndTagCounts[2],digits=2),")","")            
             message("Percentage Of Reads in Blacklisted Regions: ",round((object@CountsInBlackList/object@FlagAndTagCounts[2])*100),"")            
             
             message("SSD: ",object@SSD,"")
             message("Fragment Length Cross-Coverage: ",FragmentLengthCrossCoverage(object),"")
             message("Relative Cross-Coverage: ",RelativeCrossCoverage(object),"")
             
             message("Percentage Of Reads in GenomicFeature: ")
             print(data.frame(ProportionOfCounts=unlist(object@CountsInFeatures),row.names=gsub("CountsIn","",names(object@CountsInFeatures)))/object@FlagAndTagCounts[4])
             
             message("Percentage Of Reads in Peaks: ",round((object@CountsInPeaks/object@FlagAndTagCounts[4])*100,digits=2),"")            
             message("Number of Peaks: ",length(object),"")    
             
             print(as(object,"GRanges"))
          }
)

setGeneric("crosscoverage", function(object="ChIPQCsample") standardGeneric("crosscoverage"))
setMethod("crosscoverage", signature(object="ChIPQCsample"), function(object) ((object@CrossCoverage[1]-object@CrossCoverage)/object@CrossCoverage[1]))
setGeneric("ssd", function(object="ChIPQCsample") standardGeneric("ssd"))
setMethod("ssd", "ChIPQCsample", function(object) object@SSD)
setGeneric("fragmentlength", function(object="ChIPQCsample",width) standardGeneric("fragmentlength"))
setMethod("fragmentlength", "ChIPQCsample", function(object,width){
   if(missing(width)) {
      width = readlength(object)
   }
   MaxShift <- which.max(running(crosscoverage(object)[-seq(1:(2*readlength(object)))],width=width,allow.fewer=TRUE))+2*readlength(object)
   return(unname(MaxShift))
})

setGeneric("FragmentLengthCrossCoverage", function(object="ChIPQCsample",width) standardGeneric("FragmentLengthCrossCoverage"))
setMethod("FragmentLengthCrossCoverage", signature(object="ChIPQCsample"), function(object){ 
   FragmentLengthCrossCoverage <- crosscoverage(object)[fragmentlength(object,10)]-crosscoverage(object)[1]
   return(FragmentLengthCrossCoverage)
}
)
setGeneric("ReadLengthCrossCoverage", function(object="ChIPQCsample",width) standardGeneric("ReadLengthCrossCoverage"))
setMethod("ReadLengthCrossCoverage", signature(object="ChIPQCsample"), function(object){ 
   ReadLengthCrossCoverage <-  crosscoverage(object)[readlength(object)]-crosscoverage(object)[1]
   return(ReadLengthCrossCoverage)
}
)

setGeneric("RelativeCrossCoverage", function(object="ChIPQCsample",width) standardGeneric("RelativeCrossCoverage"))
setMethod("RelativeCrossCoverage", signature(object="ChIPQCsample"), function(object){ 
   RelativeCrossCoverage <- FragmentLengthCrossCoverage(object)/ReadLengthCrossCoverage(object)
   return(RelativeCrossCoverage)
}
)


setGeneric("flagtagcounts", function(object="ChIPQCsample") standardGeneric("flagtagcounts"))
setMethod("flagtagcounts", "ChIPQCsample", function(object) object@FlagAndTagCounts)
setGeneric("flagtagcounts", function(object="ChIPQCsample") standardGeneric("flagtagcounts"))
setMethod("flagtagcounts", "ChIPQCsample", function(object) object@FlagAndTagCounts)
setGeneric("coveragehistogram", function(object="ChIPQCsample") standardGeneric("coveragehistogram"))
setMethod("coveragehistogram", "ChIPQCsample", function(object) object@CoverageHistogram)
setGeneric("averagepeaksignal", function(object="ChIPQCsample") standardGeneric("averagepeaksignal"))
setMethod("averagepeaksignal", "ChIPQCsample", function(object) object@AveragePeakSignal[[1]])
setGeneric("Normalisedaveragepeaksignal", function(object="ChIPQCsample") standardGeneric("Normalisedaveragepeaksignal"))
setMethod("Normalisedaveragepeaksignal", "ChIPQCsample", function(object) object@AveragePeakSignal[[2]])
setGeneric("peaks", function(object="ChIPQCsample") standardGeneric("peaks"))
setMethod("peaks", "ChIPQCsample", function(object) as(object,"GRanges"))
setGeneric("readlength", function(object="ChIPQCsample") standardGeneric("readlength"))
setMethod("readlength", "ChIPQCsample", function(object) object@readlength)

setGeneric("PropGenomeInFeature", function(object="ChIPQCsample") standardGeneric("PropGenomeInFeature"))
setMethod("PropGenomeInFeature", "ChIPQCsample", function(object) {
   PropInFeatures <- object@PropInFeatures
   names(PropInFeatures) <- gsub("PropIn","",names(PropInFeatures))
   return(PropInFeatures)
})

setGeneric("CountsInFeatures", function(object="ChIPQCsample") standardGeneric("CountsInFeatures"))
setMethod("CountsInFeatures", "ChIPQCsample", function(object)  {
   CountsInFeatures <- object@CountsInFeatures
   names(CountsInFeatures) <- gsub("CountsIn","",names(CountsInFeatures))
   return(CountsInFeatures)
})

setGeneric("PropCountsInFeatures", function(object="ChIPQCsample") standardGeneric("PropCountsInFeatures"))
setMethod("PropCountsInFeatures", "ChIPQCsample", function(object){
   return(as.list(unlist(CountsInFeatures(object))/mapped(object)))
})


setGeneric("regi", function(object="ChIPQCsample") standardGeneric("regi"))
setMethod("regi", "ChIPQCsample", function(object){
   PropCountInFeatures <- data.frame(PropCountInFeatures=unlist(PropCountsInFeatures(object)),row.names=names(PropCountsInFeatures(object)))
   if(sum(is.na(PropCountInFeatures))>0) {
     #warning('No genomic features computed',call.=FALSE)
     savenames = rownames(PropCountInFeatures)
     PropCountInFeatures = PropCountInFeatures[,1]
     names(PropCountInFeatures) = savenames
     return(PropCountInFeatures)
   }
   PropGenomeInFeatures <- data.frame(PropGenomeInFeature=unlist(PropGenomeInFeature(object)),row.names=names(PropGenomeInFeature(object)))
   regiFrame <- merge(PropCountInFeatures,PropGenomeInFeatures,by=0,all=FALSE,sort=FALSE)
   regi <- log2(regiFrame[,"PropCountInFeatures"]/regiFrame[,"PropGenomeInFeature"])
   names(regi) <- regiFrame[,"Row.names"]
   #regi = regi[c(1,2,3,7,6,5,4)]
   #names(regi)[2] = "5UTRs"
   return(regi)
})


setGeneric("frip", function(object="ChIPQCsample") standardGeneric("frip"))
setMethod("frip", "ChIPQCsample", function(object){
   CountsInPeaks <- object@CountsInPeaks  
   TotalCounts <- object@FlagAndTagCounts["Mapped"]
   FRIP <- unname(CountsInPeaks/TotalCounts)  
   return(FRIP)
}
)
setGeneric("rip", function(object="ChIPQCsample") standardGeneric("rip"))
setMethod("rip", "ChIPQCsample", function(object){
   CountsInPeaks <- unname(object@CountsInPeaks)
   return(CountsInPeaks)
}
)
setGeneric("ribl", function(object="ChIPQCsample") standardGeneric("ribl"))
setMethod("ribl", "ChIPQCsample", function(object){
   CountsInBlackList <- unname(object@CountsInBlackList)
   return(CountsInBlackList)
}
)

setGeneric("mapped", function(object="ChIPQCsample") standardGeneric("mapped"))
setMethod("mapped", "ChIPQCsample", function(object){
   MappedCounts <- unname(object@FlagAndTagCounts[2])  
   return(MappedCounts)
}
)


setGeneric("QCmetrics", function(object="ChIPQCsample") standardGeneric("QCmetrics"))
setMethod("QCmetrics", "ChIPQCsample", function(object){
   res        = c(reads(object,FALSE),
                  signif((mapped(object)/reads(object,FALSE))*100,3),
                  signif((1-reads(object,TRUE)/reads(object,FALSE))*100,3),
                  signif(duplicateRate(object)*100,3),
                  readlength(object),
                  fragmentlength(object,width=readlength(object)),
                  signif(RelativeCrossCoverage(object),3),
                  #signif(FragmentLengthCrossCoverage(object),3),
                  #signif(ReadLengthCrossCoverage(object),3),
                  signif(ssd(object),3),
                  signif(frip(object)*100,3))
   names(res) = c("Reads",
                  "Map%",
                  "Filt%",
                  "Dup%",
                  "ReadL",
                  "FragL",
                  "RelCC",
                  #"FragLenCC",
                  #"ReadLenCC",
                  "SSD",
                  "RiP%")
   blk = ribl(object)
   if(!is.na(blk)) {
      names(blk) <- "RiBL%"
      blk = signif(blk/res[1]*100,3)
      res = c(res,blk)
   }
   return(res)
})

setGeneric("reads", function(object="ChIPQCsample", bFiltered) standardGeneric("reads"))
setMethod("reads", "ChIPQCsample", function(object,bFiltered){
   if(missing(bFiltered)) bFiltered=TRUE
   if(!bFiltered) {
      res = object@FlagAndTagCounts[1] + object@FlagAndTagCounts[2]
   } else {
      res = object@FlagAndTagCounts[4]
   }
   return(res)
})

setGeneric("duplicates", function(object="ChIPQCsample", bFiltered) standardGeneric("duplicates"))
setMethod("duplicates", "ChIPQCsample", function(object,bFiltered){
   if(missing(bFiltered)) bFiltered=TRUE
   if(bFiltered) {
      res = object@FlagAndTagCounts[5]
   } else {
      res = object@FlagAndTagCounts[3]
   }
   return(res)
})

setGeneric("duplicateRate", function(object="ChIPQCsample", bFiltered) standardGeneric("duplicateRate"))
setMethod("duplicateRate", "ChIPQCsample", function(object, bFiltered){
   if(missing(bFiltered)) bFiltered=TRUE
   res = duplicates(object,bFiltered) / reads(object, bFiltered)
   return(res)
})






