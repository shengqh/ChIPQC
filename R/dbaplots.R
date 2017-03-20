do_plotCorHeatmap = function(object,attributes,...) {
   if(sum(as.numeric(as.vector(dba.show(QCdba(object))$Intervals)))==0) {
      message('No peaks to plot.')      
   }else{
      atts=getAtts(attributes=attributes)  
      x = suppressWarnings(dba.plotHeatmap(object,attributes=atts,...))
      invisible(x)
   }
}

do_plotPCA = function(object,attributes,label,...) {
   if(sum(as.numeric(as.vector(dba.show(QCdba(object))$Intervals)))==0) {
      message('No peaks to plot.')      
   }else{
   atts=getAtts(attributes=attributes)
   if(!missing(label)){
      label=getAtts(attributes=label)
   }
   x = suppressWarnings(dba.plotPCA(object,attributes=atts,label=label,...))
   return(x)
   }
}

setGeneric("plotCorHeatmap", function(object="ChIPQCexperiment",attributes,...) standardGeneric("plotCorHeatmap"))
setMethod("plotCorHeatmap", "ChIPQCexperiment", do_plotCorHeatmap)

setGeneric("plotPrincomp", function(object="ChIPQCexperiment",attributes,...) standardGeneric("plotPrincomp"))
setMethod("plotPrincomp", "ChIPQCexperiment", do_plotPCA)

getAtts = function(attributes) {
   atts=NULL
   if(missing(attributes)) {
      atts = DBA_ID
   } else {
      if("ID" %in% attributes) {
         atts = c(atts,DBA_ID)
      }
      if("Tissue" %in% attributes) {
         atts = c(atts,DBA_TISSUE)
      }
      if("Factor" %in% attributes) {
         atts = c(atts,DBA_FACTOR)
      }
      if("Condition" %in% attributes) {
         atts = c(atts,DBA_CONDITION)
      }
      if("Treatment" %in% attributes) {
         atts = c(atts,DBA_TREATMENT)
      }
      if("Replicate" %in% attributes) {
         atts = c(atts,DBA_REPLICATE)
      }
   }
   
   if(is.null(atts)) atts = DBA_ID
   return(atts)
}
