sampleQC <- function(bamFile,bedFile=NULL,blklist=NULL,ChrOfInterest=NULL,GeneAnnotation=NULL,Window=400,FragmentLength=150,
                     shiftWindowStart=1,shiftWindowEnd=300,mapQCutoff=15,runCrossCor=FALSE){
  #    require(Rsamtools)
  #    require(GenomicRanges)
  #    require(GenomicAlignments)
  
  ChrLengths <- scanBamHeader(bamFile)[[1]]$targets
  
  if(length(ChrLengths[ChrLengths < shiftWindowEnd - shiftWindowStart]) > 0){
    message("Removing ",length(ChrLengths[ChrLengths < shiftWindowEnd - shiftWindowStart]),
            " chromosomes with length less than cross-coverage shift")
    ChrLengths <- ChrLengths[!ChrLengths < shiftWindowEnd - shiftWindowStart]
  }
  message("Bam file has ",length(names(ChrLengths))," contigs")
  if(!all(ChrOfInterest %in% names(ChrLengths))){
    stop("Contigs of interest are not all in Bam file!!")
  }
  if(!is.null(ChrOfInterest)){
    ChrLengths <- ChrLengths[names(ChrLengths) %in% ChrOfInterest]
  }
  
  ShiftMat <- NULL
  ShiftMatCor <- NULL   
  FlagTagCounts <- NULL
  CovHist <- NULL
  Cov <- NULL
  SSD <- NULL
  SSDBL <- NULL
  readlength <- as.numeric(NA)
  CoverageMatrix <- NULL
  CoverageMatrixInput <- NULL
  bedRangesTemp <- GetGRanges(GRanges(),AllChr=names(ChrLengths))
  bedRangesSummitsTemp <- vector("numeric")
  Counts <- NULL
  MultiRangesCountsTemp <- NULL
  txdb <- NULL
  GeneAnnotationTotalSizes <- NULL
  if(class(GeneAnnotation)!="list" & !is.null(GeneAnnotation)) {
    GeneAnnotation = getAnnotation(GeneAnnotation,AllChr=names(ChrLengths))
    
  }
  if(class(GeneAnnotation)=="list"){
    if(length(GeneAnnotation)>1) {
      GeneAnnotation <- GeneAnnotation[!names(GeneAnnotation) %in% "version"]
      GeneAnnotation <- lapply(GeneAnnotation,function(x)reduce(GetGRanges(x,names(ChrLengths))))
      GeneAnnotationTotalSizes <- unlist(lapply(GeneAnnotation,function(x)sum(width(x))))
    } else {
      GeneAnnotation=NULL
    }
  }
  
  for(k in 1:length(ChrLengths)){
    
    Param <- ScanBamParam(which=GRanges(seqnames=names(ChrLengths)[k],IRanges(start=1,end=unname(ChrLengths[names(ChrLengths) == names(ChrLengths)[k]])-shiftWindowEnd)),
                          what=c("flag","mapq"))
    temp <- readGAlignmentsFromBam(bamFile,param=Param)
    if(length(temp) < 1){
      
      emptyChr_SSD <- 0
      names(emptyChr_SSD) <- names(ChrLengths)[k]
      SSD <- c(SSD,emptyChr_SSD)
      
      emptyChr_SSDBL <- 0
      names(emptyChr_SSDBL) <- names(ChrLengths)[k]
      SSDBL <- c(SSD,emptyChr_SSDBL)
      
      
      ShiftMattemp <- matrix(rep(0,(shiftWindowEnd-shiftWindowStart)+1),ncol=1)
      colnames(ShiftMattemp) <- names(ChrLengths)[k]
      ShiftMat <- cbind(ShiftMat,ShiftMattemp)
    }else{
      if(k == 1){
        tocheckforreads <- 1000
        readlength=round(mean(width(temp[1:tocheckforreads])))
      }
      
      Sample_GIT <- GIntervalTree(GRanges(seqnames=seqnames(temp),ranges=ranges(temp),strand=strand(temp),elementMetadata(temp)))
      
      flagMapQ <- cbind(bamFlagAsBitMatrix(Sample_GIT$flag)[,c("isUnmappedQuery","isDuplicate")],Sample_GIT$mapq)
      colnames(flagMapQ) <- c("A","B","C")
      flagMapQ[is.na(flagMapQ[,"C"]),"C"] <- Inf
      temp <- as.data.frame(xtabs(~A+B+C, data=flagMapQ))
      UnMapped <- sum(temp[temp[,"A"] == 1,"Freq"])
      Mapped <- sum(temp[temp[,"A"] != 1,"Freq"])
      Duplicates <- sum(temp[temp[,"A"] != 1 & temp[,"B"] == 1,"Freq"])
      
      MapQPass <- sum(temp[temp[,"A"] != 1 & as.numeric(as.vector(temp[,"C"])) >= 15,"Freq"])
      MapQPassAndDup <- sum(temp[temp[,"A"] != 1 & temp[,"B"] == 1 & as.numeric(as.vector(temp[,"C"])) >= mapQCutoff,"Freq"])
      FlagTagCounts <- rbind(cbind(UnMapped,Mapped,Duplicates,MapQPass,MapQPassAndDup),FlagTagCounts)
      
      
      Cov <- coverage(Sample_GIT,width=unname(ChrLengths[k]))
      message("Calculating coverage histogram for ",names(ChrLengths)[k],"\n")
      CovHist <- c(CovHist,list(table(Cov[[which(names(Cov) %in% names(ChrLengths)[k])]])))
      message("Calculating SSD for ",names(ChrLengths)[k],"\n")
      SSD <- c(SSD,sd(Cov)[names(ChrLengths)[k]])
      
      PosCoverage <- coverage(IRanges(start(Sample_GIT[strand(Sample_GIT)=="+"]),start(Sample_GIT[strand(Sample_GIT)=="+"])),width=ChrLengths[k])
      NegCoverage <- coverage(IRanges(end(Sample_GIT[strand(Sample_GIT)=="-"]),end(Sample_GIT[strand(Sample_GIT)=="-"])),width=ChrLengths[k])
      message("Calculating shift for ",names(ChrLengths)[k],"\n")
      #ShiftsTemp <- shiftApply(seq(shiftWindowStart,shiftWindowEnd),PosCoverage,NegCoverage,cor)
      
      ShiftsTemp <- shiftApply(seq(shiftWindowStart,shiftWindowEnd),PosCoverage,NegCoverage,RleSumAny, verbose = TRUE)         
      ShiftMat <- cbind(ShiftMat,ShiftsTemp)
      colnames(ShiftMat)[ncol(ShiftMat)] <- names(ChrLengths)[k]
      if(runCrossCor==TRUE){
        ShiftsCorTemp <- shiftApply(seq(shiftWindowStart,shiftWindowEnd),PosCoverage,NegCoverage,cor, verbose = TRUE)         
        ShiftMatCor <- cbind(ShiftMatCor,ShiftsCorTemp)
        colnames(ShiftMatCor)[ncol(ShiftMatCor)] <- names(ChrLengths)[k]
      }
      
      if(!is.null(bedFile)){
        bedRanges <- GetGRanges(bedFile,as.vector(names(ChrLengths)),names(ChrLengths)[k])
        GRangesOfInterestList <- GRangesList(GRanges(seqnames(bedRanges),ranges(bedRanges)))
        names(GRangesOfInterestList) <- c(names(GRangesOfInterestList),"Peaks")
        seqlevels(GRangesOfInterestList) <- names(ChrLengths)
      }else{
        GRangesOfInterestList <- GRangesList()
        seqlevels(GRangesOfInterestList) <- names(ChrLengths)
      }
      
      if(!is.null(blklist)){
        blkRanges <- GetGRanges(blklist,names(ChrLengths),names(ChrLengths)[k])           
        GRangesOfInterestListBL <- GRangesList(GRanges(seqnames(blkRanges),ranges(blkRanges)))
        names(GRangesOfInterestListBL) <- "BlackList"
        GRangesOfInterestList <- c(GRangesOfInterestList,GRangesOfInterestListBL)
        Cov[[names(ChrLengths)[k]]][ranges(blkRanges)] <- 0
        SSDBL <- c(SSDBL,sd(Cov)[names(ChrLengths)[k]])
        
        #print(names(GRangesOfInterestList))
      }else{
        SSDBL <- NULL
      }
      if(class(GeneAnnotation)=="list"){
        #GRangesOfInterestListGA <- GRangesList()
        noVersionGeneAnnotation <- GeneAnnotation[!names(GeneAnnotation)=="version"]
        GRangesOfInterestListGF <- GRangesList()
        seqlevels(GRangesOfInterestListGF) <- names(ChrLengths)
        for(g in 1:length(noVersionGeneAnnotation)){ 
          GRangesOfInterestListGF <- c(
            GRangesOfInterestListGF,
            GRangesList(GetGRanges(noVersionGeneAnnotation[[g]],names(ChrLengths),names(ChrLengths)[k],simplify=TRUE))
          )
        }
        names(GRangesOfInterestListGF) <- names(noVersionGeneAnnotation)
        GRangesOfInterestList <- c(GRangesOfInterestList,GRangesOfInterestListGF)
        #GRangesOfInterestListGA <- 
        #GRangesList(GetGRanges(GeneAnnotation$All5utrs,names(ChrLengths),names(ChrLengths)[k],simplify=TRUE),
        #            GetGRanges(GeneAnnotation$All3utrs,names(ChrLengths),names(ChrLengths)[k],simplify=TRUE),
        #            GetGRanges(GeneAnnotation$Allcds,names(ChrLengths),names(ChrLengths)[k],simplify=TRUE),
        #            GetGRanges(GeneAnnotation$Allintrons,names(ChrLengths),names(ChrLengths)[k],simplify=TRUE),
        #            GetGRanges(GeneAnnotation$Alltranscripts,names(ChrLengths),names(ChrLengths)[k],simplify=TRUE),
        #            GetGRanges(GeneAnnotation$Promoters500,names(ChrLengths),names(ChrLengths)[k],simplify=TRUE),
        #            GetGRanges(GeneAnnotation$Promoters2000to500,names(ChrLengths),names(ChrLengths)[k],simplify=TRUE),
        #            GetGRanges(GeneAnnotation$LongPromoter20000to2000,names(ChrLengths),names(ChrLengths)[k],simplify=TRUE))
        #names(GRangesOfInterestListGA) <-
        #   c("All5utrs","All3utrs","Allcds","Allintrons","Alltranscripts",
        #     "Promoters500","Promoters2000to500","LongPromoter20000to2000")
        
        
        
        #GRangesOfInterestList <- c(GRangesOfInterestList,GRangesOfInterestListGA)
      }
      GRangesOfInterestList <- GRangesOfInterestList
      MultiRangesCountsTemp <- c(MultiRangesCountsTemp,countOverlaps(GRangesOfInterestList,Sample_GIT))
      message("Counting reads in features for ",names(ChrLengths)[k],"\n")        
      if(!is.null(bedFile)){
        CountsTemp <- countOverlaps(bedRanges,Sample_GIT)
        Counts  <- c(Counts,CountsTemp)
        
        bedRangesTemp <- c(bedRangesTemp,bedRanges)
        message("Signal over peaks for ",names(ChrLengths)[k],"\n")                   
        AllFragRanges <- resize(as(Sample_GIT[Sample_GIT %over% bedRanges],"GRanges"),FragmentLength,"start")
        bedRangesSummits <- findCovMaxPos(AllFragRanges,bedRanges,ChrLengths[k],FragmentLength)
        bedRangesSummitsTemp <- c(bedRangesSummitsTemp,as.numeric(as.vector(start(bedRangesSummits))))
        updatedRanges  <- resize(bedRangesSummits,Window,"center")
        AllCov <- coverage(AllFragRanges,width=unname(ChrLengths[k])+Window)
        CoverageMatrix <- rbind(CoverageMatrix,matrix(as.vector(AllCov[[which(names(AllCov) == names(ChrLengths)[k])]][ranges(updatedRanges[seqnames(updatedRanges) == names(ChrLengths)[k]])]),ncol=Window,byrow=TRUE))
      }
    }
  }
  if(!is.null(FlagTagCounts)){
    FlagTagCounts <- colSums(FlagTagCounts)
  }else{
    FlagTagCounts <- as.numeric(data.frame(UnMapped=0,Mapped=0,Duplicates=0,MapQPass=0,MapQPassAndDup=0))
  }
  Weights <- ChrLengths
  CovHistAll <- NULL
  if(!is.null(CovHist)){
    for(l in 1:length(CovHist)){
      print(length(l))
      tempCovHist <- cbind(as.numeric(names(CovHist[[l]])),as.numeric(CovHist[[l]]))
      tempCovHist <- merge(CovHistAll,tempCovHist,by.x=0,by.y=1,all=TRUE,sort=FALSE)
      CovSums <- data.frame(Depth=rowSums(tempCovHist[,-1,drop=FALSE],na.rm=TRUE),row.names=tempCovHist[,1])
      CovHistAll <- CovSums
    }
    tempCovHistAll <- as.numeric(CovHistAll[,1])
    names(tempCovHistAll) <- rownames(CovHistAll)
    CovHistAll <- tempCovHistAll
  }else{
    CovHistAll <- as.numeric(NA)     
  }
  if(!is.null(ShiftMat)){
    ShiftsAv <- apply(ShiftMat,1,function(x)weighted.mean(x,Weights[colnames(ShiftMat)],na.rm=TRUE))
  }else{
    ShiftsAv <- as.numeric(NA)
  }
  if(!is.null(ShiftMatCor)){
    ShiftsCorAv <- apply(ShiftMatCor,1,function(x)weighted.mean(x,Weights[colnames(ShiftMatCor)],na.rm=TRUE))
  }else{
    ShiftsCorAv <- as.numeric(NA)
  }   
  if(!is.null(SSD)){   
    SSDAv <- unname((weighted.mean(SSD[names(ChrLengths)],ChrLengths)*1000)/sqrt(FlagTagCounts[4]))
  }else{
    SSDAv <- as.numeric(NA)     
  }
  
  if(!is.null(MultiRangesCountsTemp)){
    GFCountsMatrix <- matrix(MultiRangesCountsTemp,ncol=length(GRangesOfInterestList),byrow=TRUE)
    colnames(GFCountsMatrix) <- unique(names(GRangesOfInterestList))
  }else{
    GFCountsMatrix <- NULL
  }
  
  if(!is.null(bedFile)){
    AvProfile <- colMeans(CoverageMatrix)
    NormAvProfile <- (AvProfile/FlagTagCounts[4])*1e6
    elementMetadata(bedRangesTemp) <- data.frame(Counts,bedRangesSummitsTemp)
    #print(length(GRangesOfInterestList))
    ReadsInPeaks <- sum(GFCountsMatrix[,"Peaks"])
  }else{
    AvProfile <- NA
    NormAvProfile <- NA
    bedRangesTemp <- GRanges()
    #print(length(GRangesOfInterestList))
    ReadsInPeaks <- as.numeric(NA)
  }  
  if(!is.null(blklist)){   
    ReadsInBLs <- sum(GFCountsMatrix[,"BlackList"])  
  }else{
    ReadsInBLs <- as.numeric(NA)
  }
  if(!is.null(SSDBL)){   
    SSDBLAv <- unname((weighted.mean(SSDBL[names(ChrLengths)],ChrLengths)*1000)/(sqrt(FlagTagCounts[4]-ReadsInBLs)))
  }else{
    SSDBLAv <- as.numeric(NA)     
  }   
  
  if(class(GeneAnnotation)=="list" & !is.null(GFCountsMatrix)){
    #CountsInFeatures <-vector("list",length=ncol(GFCountsMatrix))
    CountsInFeatures <- as.list(apply(GFCountsMatrix,2,function(x)sum(x,na.rm=TRUE)))
    names(CountsInFeatures) <- colnames(GFCountsMatrix)
    #ReadsIn5utrs<- sum(GFCountsMatrix[,"All5utrs"])
    #ReadsIn3utrs <- sum(GFCountsMatrix[,"All3utrs"])
    #ReadsInCDS <- sum(GFCountsMatrix[,"Allcds"])
    #ReadsInIntrons <- sum(GFCountsMatrix[,"Allintrons"])
    #ReadsInTranscripts <- sum(GFCountsMatrix[,"Alltranscripts"])              
    #ReadsInPromoters500 <- sum(GFCountsMatrix[,"Promoters500"])              
    #ReadsInPromoters2000to500 <- sum(GFCountsMatrix[,"Promoters2000to500"])              
    #ReadsInLongPromoter20000to2000 <- sum(GFCountsMatrix[,"LongPromoter20000to2000"])              
    #CountsInFeatures <- list(CountsIn5utrs=ReadsIn5utrs,CountsIn3UTRs=ReadsIn3utrs,
    #                         CountsinReads=ReadsInCDS,CountsInIntrons=ReadsInIntrons,CountsInTranscripts=ReadsInTranscripts,
    #                         CountsInPromoters500=ReadsInPromoters500,CountsInPromoters2000to500=ReadsInPromoters2000to500,
    #                         CountsInLongPromoter20000to2000=ReadsInLongPromoter20000to2000)
    
    PropInFeatures <- as.list(GeneAnnotationTotalSizes/sum(as.numeric(ChrLengths)))
    #TotalLength5utrs<- GeneAnnotationTotalSizes["All5utrs"]
    #TotalLength3utrs <- GeneAnnotationTotalSizes["All3utrs"]
    #TotalLengthCDS <- GeneAnnotationTotalSizes["Allcds"]
    #TotalLengthIntrons <- GeneAnnotationTotalSizes["Allintrons"]
    #TotalLengthTranscripts <- GeneAnnotationTotalSizes["Alltranscripts"]              
    #TotalLengthPromoters500 <- GeneAnnotationTotalSizes["Promoters500"]              
    #TotalLengthPromoters2000to500 <- GeneAnnotationTotalSizes["Promoters2000to500"]              
    #TotalLengthLongPromoter20000to2000 <- GeneAnnotationTotalSizes["LongPromoter20000to2000"]
    #PropInFeatures <- list(PropIn5utrs=TotalLength5utrs/sum(as.numeric(ChrLengths)),PropIn3UTRs=TotalLength3utrs/sum(as.numeric(ChrLengths)),
    #                       PropInCDS=TotalLengthCDS/sum(as.numeric(ChrLengths)),PropInIntrons=TotalLengthIntrons/sum(as.numeric(ChrLengths)),PropInTranscripts=TotalLengthTranscripts/sum(as.numeric(ChrLengths)),
    #                       PropInPromoters500=TotalLengthPromoters500/sum(as.numeric(ChrLengths)),PropInPromoters2000to500=TotalLengthPromoters2000to500/sum(as.numeric(ChrLengths)),
    #                       PropInLongPromoter20000to2000=TotalLengthLongPromoter20000to2000/sum(as.numeric(ChrLengths))
    
    
    #)
    
    
    
    
  }else{
    ReadsIn5utrs<- NA
    ReadsIn3utrs <- NA
    ReadsInCDS <- NA
    ReadsInIntrons <- NA
    ReadsInTranscripts <- NA
    ReadsInPromoters500 <- NA
    ReadsInPromoters2000to500 <- NA
    ReadsInLongPromoter20000to2000 <- NA
    
    CountsInFeatures <- list(CountsIn5utrs=ReadsIn5utrs,CountsIn3UTRs=ReadsIn3utrs,
                             CountsinReads=ReadsInCDS,CountsInIntrons=ReadsInIntrons,CountsInTranscripts=ReadsInTranscripts,
                             CountsInPromoters500=ReadsInPromoters500,CountsInPromoters2000to500=ReadsInPromoters2000to500,
                             CountsInLongPromoter20000to2000=ReadsInLongPromoter20000to2000) 
    
    
    TotalLength5utrs<- NA
    TotalLength3utrs <- NA
    TotalLengthCDS <- NA
    TotalLengthIntrons <- NA
    TotalLengthTranscripts <- NA
    TotalLengthPromoters500 <- NA
    TotalLengthPromoters2000to500 <- NA
    TotalLengthLongPromoter20000to2000 <- GeneAnnotationTotalSizes["LongPromoter20000to2000"]
    PropInFeatures <- list(PropIn5utrs=TotalLength5utrs,PropIn3UTRs=TotalLength3utrs,
                           PropInReads=TotalLengthCDS,PropInIntrons=TotalLengthIntrons,PropInTranscripts=TotalLengthTranscripts,
                           PropInPromoters500=TotalLengthPromoters500,PropInPromoters2000to500=TotalLengthPromoters2000to500,
                           PropInLongPromoter20000to2000=TotalLengthLongPromoter20000to2000)
    
    
  }
  
  
  
  CHQC <- new("ChIPQCsample", bedRangesTemp,AveragePeakSignal=list(AvProfile,NormAvProfile),
              CrossCoverage=ShiftsAv,CrossCorrelation=ShiftsCorAv,SSD=SSDAv,SSDBL=SSDBLAv,CountsInPeaks=ReadsInPeaks,
              CountsInBlackList=ReadsInBLs,
              CountsInFeatures=CountsInFeatures,PropInFeatures=PropInFeatures,
              CoverageHistogram=CovHistAll,FlagAndTagCounts=FlagTagCounts,readlength=readlength)
  
  return(CHQC)    
  
  
}  

RleSumAny <- function (e1, e2)
{
  len <- length(e1)
  stopifnot(len == length(e2))
  x1 <- runValue(e1); s1 <- cumsum(runLength(e1))
  x2 <- runValue(e2); s2 <- cumsum(runLength(e2))
  .Call("rle_sum_any",
        as.integer(x1), as.integer(s1),
        as.integer(x2), as.integer(s2),
        as.integer(len),
        PACKAGE = "chipseq")
}

GetGRanges <- function(LoadFile,AllChr=NULL,ChrOfInterest=NULL,simple=FALSE,sepr="\t",simplify=FALSE){
  #    require(Rsamtools)
  #    require(GenomicRanges)
  
  if(class(LoadFile) == "GRanges"){
    RegionRanges <- LoadFile
    if(simplify){
      RegionRanges <- GRanges(seqnames(RegionRanges),ranges(RegionRanges))
    }
  }else{
    if(class(LoadFile) == "character"){
      RangesTable <- read.delim(LoadFile,sep=sepr,header=TRUE,comment="#")
    }else if(class(LoadFile) == "matrix"){
      RangesTable <- as.data.frame(LoadFile)
    } else{
      RangesTable <- as.data.frame(LoadFile)
    }
    Chromosomes <- as.vector(RangesTable[,1])
    Start <- as.numeric(as.vector(RangesTable[,2]))
    End <- as.numeric(as.vector(RangesTable[,3]))
    RegionRanges <- GRanges(seqnames=Chromosomes,ranges=IRanges(start=Start,end=End))
    if(simple == FALSE){
      if(ncol(RangesTable) > 4){
        ID <- as.vector(RangesTable[,4])
        Score <- as.vector(RangesTable[,5])
        if(ncol(RangesTable) > 6){
          Strand <- rep("*",nrow(RangesTable))
          RemainderColumn <- as.data.frame(RangesTable[,-c(1:6)])
          elementMetadata(RegionRanges) <- cbind(ID,Score,Strand,RemainderColumn)
        }else{
          elementMetadata(RegionRanges) <- cbind(ID,Score)
        }
      }
    }
  }
  if(!is.null(AllChr)){ 
    RegionRanges <- RegionRanges[seqnames(RegionRanges) %in% AllChr]    
    seqlevels(RegionRanges,force=TRUE) <- AllChr
  }
  if(!is.null(ChrOfInterest)){      
    RegionRanges <- RegionRanges[seqnames(RegionRanges) == ChrOfInterest]      
  }
  
  return(RegionRanges)
}

findCovMaxPos <- function(reads,bedRanges,ChrOfInterest,FragmentLength){
  #    require(GenomicRanges)
  #    require(Rsamtools)
  
  cat("done\n")
  cat("Calculating coverage\n")
  MaxRanges <- GRanges()
  if(length(reads) > 0){
    seqlengths(reads)[names(ChrOfInterest)] <- ChrOfInterest
    AllCov <- coverage(reads) 
    cat("Calculating Summits on ",names(ChrOfInterest)," ..")
    covPerPeak <- Views(AllCov[[which(names(AllCov) %in% names(ChrOfInterest))]],ranges(bedRanges[seqnames(bedRanges) == names(ChrOfInterest)]))
    meanSummitLocations <- viewApply(covPerPeak,function(x)round(mean(which(x==max(x)))))
    Maxes <- (start(bedRanges)+meanSummitLocations)-1
    if(any(is.na(Maxes))){ 
      NoSummitRanges <- bedRanges[is.na(Maxes)]
      Maxes[is.na(Maxes)]  <- (start((ranges(NoSummitRanges[seqnames(NoSummitRanges) == names(ChrOfInterest)])))+end((ranges(NoSummitRanges[seqnames(NoSummitRanges) == names(ChrOfInterest)]))))/2
    }
    MaxRanges <- GRanges(seqnames(bedRanges[seqnames(bedRanges) == names(ChrOfInterest)]),IRanges(start=Maxes,end=Maxes),elementMetadata=elementMetadata(bedRanges[seqnames(bedRanges) == names(ChrOfInterest)]))
    #revAllCov <- rev(coverage(reads))
    #revAllCov <- runmean(revAllCov[names(revAllCov) %in% ChrOfInterest],20)
    #cat("Calculating reverse Summits on ",ChrOfInterest," ..")
    #revMaxes <- which.max(Views(revAllCov[[which(names(revAllCov) %in% ChrOfInterest)]],ranges(bedRanges[seqnames(bedRanges) == ChrOfInterest])))
    #if(any(is.na(revMaxes))){ 
    #  revNoSummitRanges <- bedRanges[is.na(revMaxes)]
    #  revMaxes[is.na(revMaxes)]  <- (start((ranges(revNoSummitRanges[seqnames(revNoSummitRanges) == ChrOfInterest])))+end((ranges(revNoSummitRanges[seqnames(revNoSummitRanges) == ChrOfInterest]))))/2
    #}
    #revMaxRanges <- GRanges(seqnames(bedRanges[seqnames(bedRanges) == ChrOfInterest]),IRanges(start=Maxes,end=Maxes),elementMetadata=elementMetadata(bedRanges[seqnames(bedRanges) == ChrOfInterest]))
    #meanMaxes <- rowMeans(cbind(Maxes,revMaxes))
    #meanMaxRanges <- GRanges(seqnames(bedRanges[seqnames(bedRanges) == ChrOfInterest]),IRanges(start=meanMaxes,end=meanMaxes),elementMetadata=elementMetadata(bedRanges[seqnames(bedRanges) == ChrOfInterest]))
    #cat(".done\n")
  }
  #return(meanMaxRanges)
  return(MaxRanges)
}

getAnnotation = function(GeneAnnotation="hg19",AllChr){
  
  if(!is.null(GeneAnnotation)){
    if(GeneAnnotation == "hg19"){
      require(TxDb.Hsapiens.UCSC.hg19.knownGene)
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    } else if(GeneAnnotation == "hg18"){
      require(TxDb.Hsapiens.UCSC.hg18.knownGene)
      txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene        
    } else if(GeneAnnotation == "mm10"){
      require(TxDb.Mmusculus.UCSC.mm10.knownGene)
      txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene        
    } else if(GeneAnnotation == "mm9"){
      require(TxDb.Mmusculus.UCSC.mm9.knownGene)
      txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene        
    } else if(GeneAnnotation == "rn4"){
      require(TxDb.Rnorvegicus.UCSC.rn4.ensGene)
      txdb <- TxDb.Rnorvegicus.UCSC.rn4.ensGene        
    } else if(GeneAnnotation == "ce6"){
      require(TxDb.Celegans.UCSC.ce6.ensGene)
      txdb <- TxDb.Celegans.UCSC.ce6.ensGene        
    } else if(GeneAnnotation == "dm3"){
      require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
      txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene        
    }else {
      stop('Unsupported annotation:',GeneAnnotation)
    }
    All5utrs <- reduce(unique(unlist(fiveUTRsByTranscript(txdb))))
    All3utrs <- reduce(unique(unlist(threeUTRsByTranscript(txdb))))
    Allcds <- reduce(unique(unlist(cdsBy(txdb,"tx"))))
    Allintrons <- reduce(unique(unlist(intronsByTranscript(txdb))))
    Alltranscripts <- reduce(unique(unlist(transcripts(txdb))))
    
    posAllTranscripts <- Alltranscripts[strand(Alltranscripts) == "+"]
    posAllTranscripts <- posAllTranscripts[!(start(posAllTranscripts)-20000 < 0)]
    negAllTranscripts <- Alltranscripts[strand(Alltranscripts) == "-"]
    chrLimits <- seqlengths(negAllTranscripts)[as.character(seqnames(negAllTranscripts))]      
    negAllTranscripts <- negAllTranscripts[!(end(negAllTranscripts)+20000 > chrLimits)]      
    Alltranscripts <- c(posAllTranscripts,negAllTranscripts)
    Promoters500 <-  reduce(flank(Alltranscripts,500))    
    Promoters2000to500 <-  reduce(flank(Promoters500,1500))
    LongPromoter20000to2000  <- reduce(flank(Promoters2000to500,18000))
    if(!missing(AllChr) & !is.null(AllChr)){
      All5utrs <- GetGRanges(All5utrs,AllChr=AllChr)
      All3utrs <- GetGRanges(All3utrs,AllChr=AllChr)
      Allcds <- GetGRanges(Allcds,AllChr=AllChr)
      Allintrons <- GetGRanges(Allintrons,AllChr=AllChr)
      Alltranscripts <- GetGRanges(Alltranscripts,AllChr=AllChr)
      Promoters500 <- GetGRanges(Promoters500,AllChr=AllChr)
      Promoters2000to500 <-  GetGRanges(Promoters2000to500,AllChr=AllChr)
      LongPromoter20000to2000  <- GetGRanges(LongPromoter20000to2000,AllChr=AllChr)        
    }
  }
  return(list(version=GeneAnnotation,LongPromoter20000to2000=LongPromoter20000to2000,
              Promoters2000to500=Promoters2000to500,Promoters500=Promoters500,
              All5utrs=All5utrs,Alltranscripts=Alltranscripts,Allcds=Allcds,
              Allintrons=Allintrons,All3utrs=All3utrs))
}


