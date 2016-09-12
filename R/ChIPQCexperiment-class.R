
setOldClass("DBA")
setClass("ChIPQCexperiment",contains = "list", #ChIPQCsampleList",
         slots=c(Samples="list",#ChIPQCsampleList",
                 DBA="DBA",
                 annotation="list"
         ))


showChIPQCexperiment = function (object){
   meta = data.frame(as.matrix(QCmetadata(object)),stringsAsFactors=FALSE)
   cat("Samples:",length(object@Samples),":",selectSome(names(object@Samples)),"\n")
   meta2 = NULL
   if(!sum(meta$Tissue=="")) {
      meta2 = cbind(meta2,meta$Tissue)
      colnames(meta2)[ncol(meta2)] = "Tissue"
   }
   if(!sum(meta$Factor=="")) {
      meta2 = cbind(meta2,meta$Factor)
      colnames(meta2)[ncol(meta2)] = "Factor"
   }
   if(!sum(meta$Condition=="")){
      meta2 = cbind(meta2,meta$Condition)
      colnames(meta2)[ncol(meta2)] = "Condition"
   }
   if(!sum(meta$Treatment=="")) {
      meta2 = cbind(meta2,meta$Treatment)
      colnames(meta2)[ncol(meta2)] = "Treatment"
   }
   if(!sum(meta$Control==""))  {
      meta2 = cbind(meta2,meta$Control)
      colnames(meta2)[ncol(meta2)] = "Control"
   }
   meta$Replicate[is.na(meta$Replicate)]=""
   if(!sum(meta$Replicate=="")) {
      meta2 = cbind(meta2,meta$Replicate)
      colnames(meta2)[ncol(meta2)] = "Replicate"
   }
   if(!is.null(meta$Peaks)) {
      if(!sum(as.numeric(as.character(meta$Peaks)))==0) {
         meta2 = cbind(meta2,meta$Peaks)
         colnames(meta2)[ncol(meta2)] = "Peaks"
      }
   }
   rownames(meta2) = meta$ID
   print(data.frame(meta2))
   print(QCmetrics(object))
}

setMethod("show","ChIPQCexperiment",showChIPQCexperiment)

setGeneric("QCsample", function(object="ChIPQCexperiment",sampleID) standardGeneric("QCsample"))
setMethod("QCsample", "ChIPQCexperiment", function(object,sampleID){
   if(missing(sampleID)) sampleID=0
   if(is.numeric(sampleID)){
      if(sampleID>0) {
         if(sampleID<= length(object@Samples)) {
            return(object@Samples[[sampleID]])
         } else {
            stop('Invalid sample ID')
         }
      } else return(object@Samples)
   } else if(is.character(sampleID)){
      sampnum = which(names(object@Samples) %in% sampleID)
      if(length(sampnum)==0) {
         stop ('Invalid sample ID')
      } else {
         return(object@Samples[[sampnum]])
      }
   } else stop ('Invalid sample ID')
})

setGeneric("QCcontrol", function(object="ChIPQCexperiment",sampleID) standardGeneric("QCcontrol"))
setMethod("QCcontrol", "ChIPQCexperiment", function(object,sampleID){
   if(missing(sampleID)) sampleID=0
   meta = QCmetadata(object)
   controlbams = unique(meta$bamControl)
   controlbams = controlbams[!is.na(controlbams)]
   controlbams = controlbams[!is.na(controlbams)]
   controlids = which(meta$bamRead %in% controlbams)
   numcontrols = length(controlids)
   numchips = (nrow(meta)) - numcontrols
   controls = object@Samples[controlids]
   
   if(is.numeric(sampleID)) {
      if(sampleID==0) {
         return(controls)
      } else  if(sampleID <= numchips) {
         sampleID = meta$ID[sampleID]
      } else {
         stop('Invalid sample number')
      }
   }
   
   if(is.character(sampleID)){
      sampnum = which(meta$ID %in% sampleID)
      if(length(sampnum)==0) {
         stop ('Invalid sample ID')
      } else {
         controlname = meta$Control[sampnum]
         controlnum =  which(names(controls) %in% controlname)
         if(length(controlnum)==0) {
            stop ('No such control ID')
         } else {
            return(controls[[controlnum]])
         }
      }
   } else stop('invalid control ID')
})

setGeneric("QCdba", function(object="ChIPQCexperiment") standardGeneric("QCdba"))
setMethod("QCdba", "ChIPQCexperiment", function(object){
   return(object@DBA)
})

setGeneric("QCannotation", function(object="ChIPQCexperiment",bRetrieve) standardGeneric("QCannotation"))
setMethod("QCannotation", "ChIPQCexperiment", definition=function(object,bRetrieve){
   if(missing(bRetrieve)) {
      bRetrieve=FALSE
   }
   if(bRetrieve==FALSE){
      return(object@annotation$version)
   } else {
      if(length(object@annotation)>1) {
         return(object@annotation)
      } else {
         warning("Annotation unavailable for retrieval")
         return(object@annotation$version)
      }
   }
})

setGeneric("QCmetadata", function(object) standardGeneric("QCmetadata"))
setMethod("QCmetadata", "ChIPQCexperiment", function(object) {
   meta = t(object@DBA$class[-8,])
   peaks = sapply(object@DBA$peaks,nrow)
   if(sum(peaks>0)>0) {
      meta = cbind(meta,peaks)
      colnames(meta)[ncol(meta)] = "Peaks"
   }
   return(data.frame(meta))   
})

setMethod("QCmetadata", "list", function(object) {
   tmmM <- lapply(object,metadata)
   for(k in 1:length(object)){
      sampleMetadata <- as.data.frame(tmmM[[k]])
      
      if(length(sampleMetadata) ==  0){
         sampleMetadata <- data.frame(NoMetadata=TRUE)
      }
      if(k == 1){
         dfM <- t(data.frame(Sample=names(object)[k],sampleMetadata))
         colnames(dfM) <- names(object)[k]
      }else if(k == 2){
         dfMtemp <- t(data.frame(Sample=names(object)[k],sampleMetadata))
         colnames(dfMtemp) <- names(object)[k]
         dfM <- merge(dfM,dfMtemp,by=0)
      }else{
         dfMtemp <- t(data.frame(Sample=names(object)[k],sampleMetadata))
         colnames(dfMtemp) <- names(object)[k]
         dfM <- merge(dfM,dfMtemp,by.x=1,by.y=0)
      }
   }
   dfM <- t(dfM)
   colnames(dfM) <- dfM[1,]
   dfM <- dfM[-1,]
   dfM <- cbind(dfM[,"Sample",drop=FALSE],dfM[,colnames(dfM) != "Sample",drop=FALSE])
   dfM <- dfM[,colnames(dfM) != "NoMetadata",drop=FALSE]
   return(dfM)
})

setMethod("QCmetrics", "ChIPQCexperiment", function(object="ChIPQCexperiment"){
   res = t(sapply(object@Samples,QCmetrics))
   rownames(res) = names(object@Samples)
   return(res)
})

setMethod("QCmetrics", "list", function(object){
   res = t(sapply(object,QCmetrics))
   rownames(res) = names(object)
   return(res)
})
setMethod("crosscoverage", "ChIPQCexperiment", function(object) {
   res = sapply(object@Samples,crosscoverage)
   rownames(res) = 1:nrow(res)
   return(res)
}) 

setMethod("crosscoverage", "list", function(object) {
   res = sapply(object,crosscoverage)
   rownames(res) = 1:nrow(res)
   return(res)
}) 

setMethod("ssd", "ChIPQCexperiment", function(object){
   res = sapply(object@Samples,ssd)
   names(res)=names(object@Samples)
   return(res)
})

setMethod("ssd", "list", function(object){
   res = sapply(object,ssd)
   names(res)=names(object)
   return(res)
})


setMethod("fragmentlength", "ChIPQCexperiment", function(object){
   res = sapply(object@Samples,fragmentlength)
   names(res)=names(object@Samples)
   return(res)   
})

setMethod("fragmentlength", "list", function(object){
   res = sapply(object,fragmentlength)
   names(res)=names(object)
   return(res)   
})
setMethod("FragmentLengthCrossCoverage","ChIPQCexperiment", function(object){ 
   res = sapply(object@Samples,FragmentLengthCrossCoverage)
   names(res)=names(object@Samples)
   return(res)  
})
setMethod("FragmentLengthCrossCoverage","list", function(object){ 
   res = sapply(object,FragmentLengthCrossCoverage)
   names(res)=names(object)
   return(res)  
})
setMethod("ReadLengthCrossCoverage", "ChIPQCexperiment", function(object){ 
   res = sapply(object@Samples,ReadLengthCrossCoverage)
   names(res)=names(object@Samples)
   return(res)
}
)
setMethod("ReadLengthCrossCoverage", "list", function(object){ 
   res = sapply(object,ReadLengthCrossCoverage)
   names(res)=names(object)
   return(res)
}
)

setMethod("RelativeCrossCoverage", "ChIPQCexperiment", function(object){ 
   res = sapply(object@Samples,RelativeCrossCoverage)
   #names(res)=names(object@Samples)
   return(res)
}
)

setMethod("RelativeCrossCoverage", "list", function(object){ 
   res = sapply(object,RelativeCrossCoverage)
   #names(res)=names(object@Samples)
   return(res)
}
)


setMethod("flagtagcounts", "ChIPQCexperiment", function(object) {
   res = sapply(object@Samples,flagtagcounts)
   #names(res)=names(object@Samples)
   return(res)   
})

setMethod("flagtagcounts", "list", function(object) {
   res = sapply(object,flagtagcounts)
   #names(res)=names(object@Samples)
   return(res)   
})


setMethod("regi", "ChIPQCexperiment", function(object){
   res = sapply(object@Samples,regi)
   #names(res)=names(object@Samples)
   return(signif(res,3))   
})

setMethod("regi", "list", function(object){
   res = sapply(object,regi)
   #names(res)=names(object@Samples)
   return(signif(res,3))   
})

setMethod("coveragehistogram", "ChIPQCexperiment", function(object) {
   res = sapply(object@Samples,coveragehistogram)
   res = list2matrix(res)
   return(res)   
})
setMethod("coveragehistogram", "list", function(object) {
   res = sapply(object,coveragehistogram)
   res = list2matrix(res)
   return(res)   
})


setMethod("averagepeaksignal", "ChIPQCexperiment", function(object) {
   res = sapply(object@Samples,averagepeaksignal)
   
   #return(res)
   
   res = list2matrix(res)
   return(res)   
})

setMethod("averagepeaksignal", "list", function(object) {
   res = sapply(object,averagepeaksignal)
   
   #return(res)
   
   res = list2matrix(res)
   return(res)   
})


setMethod("Normalisedaveragepeaksignal", "ChIPQCexperiment", function(object) {
   res = sapply(object@Samples,Normalisedaveragepeaksignal)
   res = list2matrix(res)
   return(res)   
})
setMethod("Normalisedaveragepeaksignal", "list", function(object) {
   res = sapply(object,Normalisedaveragepeaksignal)
   res = list2matrix(res)
   return(res)   
})
setMethod("peaks", "ChIPQCexperiment", function(object){
   res = sapply(object@Samples,peaks)
   res = GRangesList(res)
   return(res)
})

setMethod("peaks", "list", function(object){
   res = sapply(object,peaks)
   res = GRangesList(res)
   return(res)
})

setMethod("readlength", "ChIPQCexperiment", function(object){
   res = sapply(object@Samples,readlength)
   names(res)=names(object@Samples)
   return(res)
})


setMethod("readlength", "list", function(object){
   res = sapply(object,readlength)
   names(res)=names(object)
   return(res)
})

setMethod("frip", "ChIPQCexperiment", function(object){
   res = sapply(object@Samples,frip)
   names(res)=names(object@Samples)
   return(res)
})

setMethod("frip", "list", function(object){
   res = sapply(object,frip)
   names(res)=names(object)
   return(res)
})

setMethod("rip", "ChIPQCexperiment", function(object){
   res = sapply(object@Samples,rip)
   names(res)=names(object@Samples)
   return(res)
})
setMethod("rip", "list", function(object){
   res = sapply(object,rip)
   names(res)=names(object)
   return(res)
})
setMethod("ribl", "ChIPQCexperiment", function(object){
   res = sapply(object@Samples,ribl)
   names(res)=names(object@Samples)
   return(res)
})
setMethod("ribl", "list", function(object){
   res = sapply(object,ribl)
   names(res)=names(object)
   return(res)
})
setMethod("mapped", "ChIPQCexperiment", function(object){
   res = sapply(object@Samples,mapped)
   names(res)=names(object@Samples)
   return(res)
})
setMethod("mapped", "list", function(object){
   res = sapply(object,mapped)
   names(res)=names(object)
   return(res)
})

setMethod("reads", "ChIPQCexperiment", function(object,bFiltered){
   if(missing(bFiltered)) bFiltered=TRUE
   res = sapply(object@Samples,reads,bFiltered)
   names(res)=names(object@Samples)   
   return(res)
})
setMethod("reads", "list", function(object,bFiltered){
   if(missing(bFiltered)) bFiltered=TRUE
   res = sapply(object,reads,bFiltered)
   names(res)=names(object)   
   return(res)
})
setMethod("duplicates", "ChIPQCexperiment", function(object,bFiltered){
   if(missing(bFiltered)) bFiltered=TRUE
   res = sapply(object@Samples,duplicates,bFiltered)
   names(res)=names(object@Samples)
   return(res)
})
setMethod("duplicates", "list", function(object,bFiltered){
   if(missing(bFiltered)) bFiltered=TRUE
   res = sapply(object,duplicates,bFiltered)
   names(res)=names(object)
   return(res)
})

setMethod("duplicateRate", "ChIPQCexperiment", function(object, bFiltered){
   if(missing(bFiltered)) bFiltered=TRUE
   res = sapply(object@Samples,duplicateRate,bFiltered)
   names(res)=names(object@Samples)
   return(res)
})

setMethod("duplicateRate", "list", function(object, bFiltered){
   if(missing(bFiltered)) bFiltered=TRUE
   res = sapply(object,duplicateRate,bFiltered)
   names(res)=names(object)
   return(res)
})

