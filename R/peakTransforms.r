#' Plot coverage of points or regions.
#'
#' @param testRanges Named character vector of region locations 
#' @param bamFile Named character vector of bamFile locations
#' @param method Method to select reproducible summits to merge.
#' @param summit Only mean avaialble
#' @param resizepeak Only asw available
#' @param overlap Type of overlap to consider for finding consensus sites
#' @param fragmentLength Predicted fragment length. Set to NULL to auto-calculate
#' @param NonPrimaryPeaks A list of parameters to deal with non primary peaks in consensus regions.
#' @return Consensus A GRanges object of consensus summits. 
#' @export
findconsensusRegions <- function(testRanges,bamFiles=NULL,method="majority",summit="mean",resizepeak="asw",overlap="any",fragmentLength=NULL,
                                 NonPrimaryPeaks=list(withinsample="drop",betweensample="mean")){

  testRanges <- GRangesList(
    bplapply(
      testRanges,
      ChIPQC:::GetGRanges)
    )
  
  ans <- lapply(
    names(bamFiles),
    function(x) summitPipeline(
      unlist(bamFiles[x]),
      unlist(testRanges[x]),
      fragmentLength=NULL,readlength=36
      )
  )
  consensusRanges <- runConsensusRegions(
    testRanges,
    method="majority",
    overlap="any"
    )
  if(unlist(NonPrimaryPeaks["withinsample"])=="drop"){
    consensensusAns <- bplapply(ans,function(x)
      dropNonPrimary(x,consensusRanges)
    )
    ansSummits <- do.call(cbind,bplapply(ans,function(x)
      extractSummits(x,consensusRanges)
    ))
    ansSummitScores <- do.call(cbind,bplapply(ans,function(x)
      extractScores(x,consensusRanges)
    )) 
    if(unlist(NonPrimaryPeaks["betweensample"])=="mean"){
      meanSummits <- rowMeans(ansSummits,na.rm=TRUE)
    }
    if(unlist(NonPrimaryPeaks["betweensample"])=="weightedmean"){
      #meanSummits <- rowMeans(ansSummits,na.omit=TRUE)
      meanSummits <- round(sapply(seq(1,nrow(ansSummits)),function(x)weighted.mean(ansSummits[x,],ansSummitScores[x,])))
    }
  }  
  start(consensusRanges) <- end(consensusRanges) <- meanSummits
  return(consensusRanges)  
}

extractSummits <- function(x,consensusRanges){
  summits <- vector("numeric",length=length(consensusRanges))
  overMat <- as.matrix(findOverlaps(consensusRanges,x))
  summits[overMat[,1]] <- as.vector(start(x[overMat[,2]]))
  return(summits)
}
extractScores <- function(x,consensusRanges,score="summitScores"){
  scores <- vector("numeric",length=length(consensusRanges))
  overMat <- as.matrix(findOverlaps(consensusRanges,x))
  scores[overMat[,1]] <- as.vector(elementMetadata(x[overMat[,2]])[,"summitScores"])
  return(scores)
}

dropNonPrimary <- function(x,consensusRanges,id="elementMetadata.ID",score="summitScores"){
  mat <- as.matrix(findOverlaps(consensusRanges,x))
  tempConsensusRanges <- consensusRanges[mat[,1],]
  elementMetadata(tempConsensusRanges) <- elementMetadata(x[mat[,2]])[,c(id,score)]
  tempConsensusRanges <- tempConsensusRanges[order(elementMetadata(tempConsensusRanges)[,score],decreasing=T),]  
  primaryIDs <- elementMetadata(tempConsensusRanges[match(unique(tempConsensusRanges[,id]),tempConsensusRanges[,id])])[,id]
  x <- x[elementMetadata(x)[,id] %in% primaryIDs]
  return(x)
}
#' Returns summits and summmit scores after optional fragment length prediction and read extension
#'
#' @param peakFile GRanges of region 
#' @param reads Character vector of bamFile location or GAlignments object
#' @param fragmentLength Predicted or calculated fragment length. Set as NULL for auto prediction of fragment length
#' @param readlength Read length of alignments.
#' @return Summits A GRanges object of summits and summmit scores.
#' @export
summitPipeline <- function(reads,peakfile,fragmentLength,readlength){
  message("Reading in peaks..",appendLF=FALSE)
  testRanges <- ChIPQC:::GetGRanges(peakfile)
  message("done")  
  if(class(reads) == "GAlignments"){
    message("Alignments loaded")
    ChrLengths <- seqlengths(reads)
  }
  if(class(reads) == "character"){
    message("Reading in alignments..",appendLF=FALSE)
    ChrLengths <- scanBamHeader(reads)[[1]]$targets
    reads <- readGAlignmentsFromBam(reads)
    message("done") 
  }
  message("Calculating fragmentlength..",appendLF=FALSE)
  ccscores <- getShifts(reads,ChrLengths,shiftWindowStart=1,shiftWindowEnd=400)
  fragmentLength <- getFragmentLength(ccscores,readlength)
  message("done")
  message("Extending reads..",appendLF=FALSE)  
  reads <- resize(as(reads,"GRanges"),fragmentLength,"start")
  message("done")
  message("Finding summit locations..",appendLF=FALSE)
  peaks <- runFindSummit(testRanges,reads,fragmentLength=NULL)
  message("done")
  message("Scoring summits..",appendLF=FALSE)
  peaks <- getSummitScore(reads,peaks,fragmentLength=NULL,score="height")
  message("done")
  return(peaks)
}

runConsensusRegions <- function(testRanges,method="majority",overlap="any"){
    if(class(testRanges) == "GRangesList" & length(testRanges) > 1){
      
      reduced <- reduce(unlist(testRanges))
      consensusIDs <- paste0("consensus_",seq(1,length(reduced)))
      elementMetadata(reduced) <- 
      do.call(cbind,lapply(testRanges,function(x)(reduced %over% x)+0))
      if(method=="majority"){
        reducedConsensus <- reduced[rowSums(as.data.frame(elementMetadata(reduced))) > length(testRanges)/2,]
      }
      if(method=="none"){
        reducedConsensus <- reduced
      }
    consensusIDs <- paste0("consensus_",seq(1,length(reducedConsensus)))
    elementMetadata(reducedConsensus) <- cbind(as.data.frame(elementMetadata(reducedConsensus)),consensusIDs)
    return(reducedConsensus)
    
  }
}
runGetShifts <- function(reads,ChrLengths,ChrOfInterestshift,shiftWindowStart=1,shiftWindowEnd=400){

reads <- reads
ChrLengths <- seqlengths(reads)
PosCoverage <- coverage(IRanges(start(reads[strand(reads)=="+"]),start(reads[strand(reads)=="+"])),width=ChrLengths[names(ChrLengths) %in% ChrOfInterestshift])
NegCoverage <- coverage(IRanges(end(reads[strand(reads)=="-"]),end(reads[strand(reads)=="-"])),width=ChrLengths[names(ChrLengths) %in% ChrOfInterestshift])
message("Calculating shift for ",ChrOfInterestshift,"\n")
ShiftsTemp <- shiftApply(seq(shiftWindowStart,shiftWindowEnd),PosCoverage,NegCoverage,RleSumAny, verbose = TRUE)         
return(ShiftsTemp)
}
getShifts <- function(reads,ChrLengths,
                      shiftWindowStart=1,shiftWindowEnd=400){
  if(is.character(reads)){
    reads <- readGAlignmentsFromBam(reads)
  }  
shiftMat <- do.call(cbind,bplapply(names(ChrLengths),function(x)
runGetShifts(reads[seqnames(reads) %in% x],ChrLengths,x,
          shiftWindowStart=1,shiftWindowEnd=400)))
cc_scores <- (rowSums(shiftMat)[1]-rowSums(shiftMat))/rowSums(shiftMat)[1]
return(cc_scores)
}

getFragmentLength <- function(x,readLength){
  #peaks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 0)+1
  MaxShift <- which.max(caTools:::runmean(x[-seq(1,(2*readLength))],10))+2*readLength
  
}
  
runFindSummit <- function(testRanges,reads,fragmentLength=NULL){
  if(is.character(reads)){
    reads <- readGAlignmentsFromBam(reads)
  }
  if(!is.null(fragmentLength)){
    message("Extending reads to fragmentlength of ",fragmentLength,"..",appendLF=F)
    reads <- resize(as(reads,"GRanges"),fragmentLength,"start")
    message("done")
  }
  test <- do.call(c,
                  bplapply(
    unique(seqnames(reads))[unique(seqnames(reads)) %in% unique(seqnames(testRanges))],
    function(x) 
    ChIPQC:::findCovMaxPos(reads[seqnames(testRanges) %in% x],testRanges[seqnames(testRanges) %in% x],seqlengths(reads)[names(seqlengths(reads)) %in% x],fragmentLength)
    )
  )
  return(test)                                        
}

getSummitScore <- function(reads,summits,fragmentLength=NULL,score="height"){
  if(is.character(reads)){
    reads <- readGAlignmentsFromBam(reads)
  }
  if(!is.null(fragmentLength)){
    message("Extending reads to fragmentlength of ",fragmentLength,"..",appendLF=F)
    reads <- resize(as(reads,"GRanges"),fragmentLength,"start")
    message("done")
  }
  test <- do.call(c,
                  bplapply(
                    as.vector(unique(seqnames(reads))[unique(seqnames(reads)) %in% unique(seqnames(summits))]),
                    function(x) 
                      runGetSummitScore(reads[seqnames(reads) %in% x],summits[seqnames(summits) %in% x],seqlengths(reads)[names(seqlengths(reads)) %in% x])
                  )
  )
  return(test)                                         
}
runGetSummitScore <- function(reads,summits,ChrOfInterestshift,FragmentLength=150,score="height"){
    

  cov <- coverage(reads)
    if(score=="height"){
      summitScores <- as.vector(unlist(cov[summits],use.names=F))
    }
  elementMetadata(summits) <- cbind(as.data.frame(elementMetadata(summits)),summitScores)
    return(summits)
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
library(BiocParallel)
library(GenomicAlignments)
