#' The ChIPprofile class
#' @export
setClass("ChIPprofile",contains = "SummarizedExperiment",
         slots=c(params="list"
         ))

#' Merge ChIPprofile objects
#'
#' A Merge two ChIPprofiles into 1 combined ChIPprofile
#' 
#' @usage
#' \S4method{mergeChIPprofiles}{ChIPprofile,ChIPprofile}(ChIPprofile)
#'
#' @docType methods
#' @name mergeChIPprofiles
#' @rdname mergeChIPprofiles
#' @aliases mergeChIPprofiles mergeChIPprofiles,ChIPprofile-method
#' 
#' @author Thomas Carroll
#'
#' @export
#' @param x A ChIPprofile object 
#' @param y A ChIPprofile object
mergeChIPprofiles <- function(x,y){
  #t1 <- rowData(x)
  selfOverlapst1 <- table(as.matrix(findOverlaps(x,x,type="equal"))[,2])
  filtDupt1 <- as.numeric(names(selfOverlapst1[
    selfOverlapst1 >1]))
  if(length(filtDupt1)){
    t1 <- x[-filtDupt1,]
  }else{
    t1 <- x
  }
  selfOverlapst2 <- table(as.matrix(findOverlaps(y,y,type="equal"))[,2])
  filtDupt2 <- as.numeric(names(selfOverlapst2[
    selfOverlapst2 >1]))
  if(length(filtDupt2)){
    t2 <- y[-filtDupt2,]
  }else{
    t2 <- x
  }
  t3 <- t1[as.matrix(findOverlaps(t1,t2,type="equal"))[,2]]
  t4 <- t2[as.matrix(findOverlaps(t1,t2,type="equal"))[,1]]
  t1t2Laps <- findOverlaps(t1,t2,type="equal")
  bothMeta <- NULL
  if(!is.null(elementMetadata(t3)) & !is.null(elementMetadata(t4))){
    tempId <- seq(1,length(t1t2Laps)) 
    t3Data <- cbind(tempId,as.data.frame(elementMetadata(t3)))
    t4Data <- as.data.frame(elementMetadata(t4))
    # Can fix this to always get some metadata if any grange does
    t4Data <- cbind(tempId,t4Data[,!colnames(t4Data) %in% colnames(t3Data)])
    bothMeta <- merge(t3Data,t4Data,by="tempId",all=T)    
  }
  elementMetadata(t3) <- bothMeta[,!colnames(bothMeta) %in% "tempId"]
  filtIndexX <- as.matrix(t1t2Laps)[,2]
  filtIndexY <- as.matrix(t1t2Laps)[,1]
  t3 <- t3[filtIndexX,]
  t4 <- t4[filtIndexY,]
  
  profileSample <- SummarizedExperiment(c(assays(t3),assays(t4)),rowData=rowData(t3))
  exptData(profileSample) <- list(names=c(exptData(x)$names,exptData(y)$names))
  paramList <- x@params
  return(new("ChIPprofile",profileSample,params=paramList))
  
}

setGeneric("mergeChIPprofiles", function(x="ChIPprofile",y="ChIPprofile") standardGeneric("mergeChIPprofiles"))

