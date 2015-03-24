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
#' \S4method{mergeChIPprofiles}{ChIPprofile}(x,y)
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
mergeChIPprofiles.ChIPprofile <- function(x,y){
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

#' @rdname mergeChIPprofiles
#' @export
setMethod("mergeChIPprofiles", signature(x="ChIPprofile",y="ChIPprofile"), mergeChIPprofiles.ChIPprofile)

setMethod("c", "ChIPprofile",
          function (x,...)
          {
            assayList <- lapply(list(x,...),function(x)assays(x)[[1]])
            subsetProfile <- SummarizedExperiment(assayList,rowData=rowData(x))
            exptData(subsetProfile)$names <- unlist(lapply(list(x,...),function(x)exptData(x)$name))
            exptData(subsetProfile)$AlignedReadsInBam <- unlist(lapply(list(x,...),function(x)exptData(x)$AlignedReadsInBam))
            return(new("ChIPprofile",subsetProfile,params=x@params))
            
          }
)


#' Normalise quantile
#'
#' Quantile normalisation across bins/regions.
#' 
#' @usage
#' \S4method{normaliseQuantiles}{ChIPprofile}(object)
#'
#' @docType methods
#' @name normaliseQuantiles
#' @rdname normaliseQuantiles
#' @aliases normaliseQuantiles normaliseQuantiles,ChIPprofile-method
#' 
#' @author Thomas Carroll
#'
#' @export
#' @param object A ChIPprofile object 
normaliseQuantiles.ChIPprofile <-  function (object)
          {
            
            tempT <- do.call(cbind,lapply(assays(object),function(x)as.vector(x)))
            l <- 1
            
            for(j in 1:ncol(object)){
              print(j)
              tempT[l:(l+nrow(object)-1),] <- normalize.quantiles(tempT[l:(l+nrow(object)-1),])
              l <- l+nrow(object)
            }
            
            qnormAssays <- lapply(1:ncol(tempT),function(x) matrix(tempT[,x],nrow=nrow(object),byrow = F))
            for(c in 1:length(qnormAssays)){
              colnames(qnormAssays[[c]]) <- colnames(assays(object)[[c]])
            }     
            normProfile <- SummarizedExperiment(qnormAssays,rowData=rowData(object))
            exptData(normProfile) <- exptData(object)
            return(new("ChIPprofile",normProfile,params=object@params))
          }  


setGeneric("normaliseQuantiles", function(object="ChIPprofile") standardGeneric("normaliseQuantiles"))

#' @rdname normaliseQuantiles
#' @export
setMethod("normaliseQuantiles", signature(object="ChIPprofile"), normaliseQuantiles.ChIPprofile)


setMethod("rbind", "ChIPprofile",
          function (x,...,deparse.level=1)
          {
            assayList <- vector("list",length=length(x))
            regions <- list(x,...)
            for(a in 1:length(x)){
              listTemp <- vector("list",length=length(regions))
              for(r in 1:length(regions)){
                listTemp[[r]] <- assays(regions[[r]])[[a]]
              }
              assayList[[a]] <- do.call(rbind,listTemp)
            }
            newRowData <- unlist(GRangesList(lapply(regions,function(x)rowData(x))))
            subsetProfile <- SummarizedExperiment(assayList,rowData=newRowData)
            exptData(subsetProfile)$names <- exptData(x)$names
            exptData(subsetProfile)$AlignedReadsInBam <- exptData(x)$AlignedReadsInBam
            return(new("ChIPprofile",subsetProfile,params=x@params))            
          }
)

setMethod("cbind", "ChIPprofile",
          function (x,...,deparse.level=1)
          {
            assayList <- vector("list",length=length(x))
            regions <- list(x,...)
            for(a in 1:length(x)){
              listTemp <- vector("list",length=length(regions))
              for(r in 1:length(regions)){
                listTemp[[r]] <- assays(regions[[r]])[[a]]
              }
              assayList[[a]] <- do.call(cbind,listTemp)
            }
            subsetProfile <- SummarizedExperiment(assayList,rowData=rowData(x))
            exptData(subsetProfile)$names <- exptData(x)$names
            exptData(subsetProfile)$AlignedReadsInBam <- exptData(x)$AlignedReadsInBam
            return(new("ChIPprofile",subsetProfile,params=x@params))            
          }
)

setMethod("[[", c("ChIPprofile", "ANY", "missing"),
          function(x, i, j, ...)
          {
            subsetProfile <- SummarizedExperiment(assays(x)[[i, ...]],rowData=rowData(x))
            exptData(subsetProfile)$names <- exptData(x)$names[i]
            exptData(subsetProfile)$AlignedReadsInBam <- exptData(x)$AlignedReadsInBam[i]
            return(new("ChIPprofile",subsetProfile,params=x@params))                        
          })



setMethod("$", "ChIPprofile",
          function(x, name)
          {
            x[[which(exptData(x)$names %in% name)]]
          })

# setMethod("Ops", signature(e1="ChIPprofile", e2="ChIPprofile"),
#           function(e1, e2) {
#             callGeneric(assays(e1)[[1]],assays(e2)[[1]])
#             assaysList <- vector("list",length=length(assays(e1)))
#             for(i in 1:length(assays(e1))){
#               assaysList[[i]] <- callGeneric(assays(e1)[[i]],assays(e2)[[i]]) 
#             }
#             return(assaysList)
#           }
# )

setMethod("Ops", signature(e1="ChIPprofile", e2="ChIPprofile"),
          function(e1, e2) {
            assayList <- vector("list",length=length(assays(e1)))
            for(i in 1:length(assays(e1))){
              assayList[[i]] <- callGeneric(assays(e1)[[i]],assays(e2)[[i]]) 
            }
            subsetProfile <- SummarizedExperiment(assayList,rowData=rowData(e1))
            exptData(subsetProfile)$names <- exptData(e1)$names[i]
            exptData(subsetProfile)$AlignedReadsInBam <- exptData(e1)$AlignedReadsInBam[i]
            return(new("ChIPprofile",subsetProfile,params=e1@params))                                    
          }
)
setMethod("Ops", signature(e1="ChIPprofile", e2="numeric"),
          function(e1, e2) {
            assayList <- vector("list",length=length(assays(e1)))
            for(i in 1:length(assays(e1))){
              assayList[[i]] <- callGeneric(assays(e1)[[i]],e2) 
            }
            subsetProfile <- SummarizedExperiment(assayList,rowData=rowData(e1))
            exptData(subsetProfile)$names <- exptData(e1)$names[i]
            exptData(subsetProfile)$AlignedReadsInBam <- exptData(e1)$AlignedReadsInBam[i]
            return(new("ChIPprofile",subsetProfile,params=e1@params))                                    
          }
)
setMethod("Ops", signature(e1="numeric", e2="ChIPprofile"),
          function(e1, e2) {
            assayList <- vector("list",length=length(assays(e2)))
            for(i in 1:length(assays(e1))){
              assayList[[i]] <- callGeneric(e1,assays(e2)[[i]]) 
            }
            subsetProfile <- SummarizedExperiment(assayList,rowData=rowData(e1))
            exptData(subsetProfile)$names <- exptData(e2)$names[i]
            exptData(subsetProfile)$AlignedReadsInBam <- exptData(e2)$AlignedReadsInBam[i]
            return(new("ChIPprofile",subsetProfile,params=e2@params))                                    
          }
)

setMethod("mean", "ChIPprofile",
          function(x, ...){
            if(missing(...)){
              assayList <- Reduce("+",assays(x))/length(assays(x))
              subsetProfile <- SummarizedExperiment(assayList,rowData=rowData(x))
              exptData(subsetProfile)$names <- paste0(exptData(x)$name,collapse="&")
              exptData(subsetProfile)$AlignedReadsInBam <- exptData(x)$AlignedReadsInBam
              return(new("ChIPprofile",subsetProfile,params=x@params))                                                  
            }else{
              x <- list(x,...)
              assayList <- vector("list",length=length(assays(x[[1]])))
              for(a in 1:length(x[[1]])){
                listTemp <- vector("list",length=length(x))
                for(r in 1:length(x)){
                  listTemp[[r]] <- assays(x[[r]])[[a]]
                }
                assayList[[a]] <- Reduce("+",listTemp)/length(x)
                subsetProfile <- SummarizedExperiment(assayList,rowData=rowData(x[[1]]))
                exptData(subsetProfile)$names <- exptData(x[[1]])$names
                exptData(subsetProfile)$AlignedReadsInBam <- unlist(lapply(x,function(x)exptData(x)$AlignedReadsInBam))
                return(new("ChIPprofile",subsetProfile,params=x[[1]]@params))                                                  
                
              }
            }
          }
)

setMethod("log2", "ChIPprofile",
          function(x){
            assayList <- lapply(assays(x),log2)
            subsetProfile <- SummarizedExperiment(assayList,rowData=rowData(x))
            exptData(subsetProfile)$names <- exptData(x)$names
            exptData(subsetProfile)$AlignedReadsInBam <- exptData(x)$AlignedReadsInBam
            return(new("ChIPprofile",subsetProfile,params=x[[1]]@params))                                                  
            
          }
)


setMethod("log", "ChIPprofile",
          function(x,base=exp(1)){
            assayList <- lapply(assays(x),function(x)log(x,base))
            subsetProfile <- SummarizedExperiment(assayList,rowData=rowData(x))
            exptData(subsetProfile)$names <- exptData(x)$names
            exptData(subsetProfile)$AlignedReadsInBam <- exptData(x)$AlignedReadsInBam
            return(new("ChIPprofile",subsetProfile,params=x[[1]]@params))                                                  
            
          }
)

zeroToMin <- function(x){
  for(r in 1:nrow(x)){
    print(r)
    
    temp[r,temp[r,] == 0] <- min(temp[r,temp[r,] != 0])
  }
  return(temp)
}

zeroToMin2 <- function(x){
  for(a in 1:length(assays(x))){
    temp <- assays(x)[[a]]
    ZeroRows <- rowSums(temp) == 0
    temp[ZeroRows,] <- min(temp[temp != 0])
    for(r in 1:nrow(temp)){
      print(r)
      temp[r,temp[r,] == 0] <- min(temp[r,temp[r,] != 0])
    }
    for(r in 1:nrow(temp)){
      print(r)
      temp[r,temp[r,] == 0] <- min(temp[r,temp[r,] != 0])
    }
    assays(x)[[a]] <- temp
  }
  return(x)                                                  
}



normalise.ChIPprofile <-  function (object,method="rpm",normFactors=NULL)
{
  
  assaylist <- assays(object)
  if(method=="rpm" & object@params$format == "bam"){
    assayList <- lapply(1:length(assays(object)),
                           function(x) assays(object)[[x]]*exptData(object)$Aligned[x]
                        )
    subsetProfile <- SummarizedExperiment(assayList,rowData=rowData(object))
    exptData(subsetProfile)$names <- exptData(object)$names
    exptData(subsetProfile)$AlignedReadsInBam <- exptData(object)$AlignedReadsInBam
    return(new("ChIPprofile",subsetProfile,params=object[[1]]@params))                                                  
    
    return(subsetProfile)
  }
  if(method=="quantile" ){
    return(normaliseQuantiles(object))
  }
  if(method=="normaliseSamples" & is.numeric(normFactors)){
    return(object*normFactors)
  }
  if(method=="normaliseRegions" & is.numeric(normFactors)){
    assayList <- assays(object)
    for(k in 1:length(assayList)){
      assayList[[k]] <- t(t(assayList[[k]]) * normFactors)
    }
    subsetProfile <- SummarizedExperiment(assayList,rowData=rowData(object))
    exptData(subsetProfile)$names <- exptData(object)$names
    exptData(subsetProfile)$AlignedReadsInBam <- exptData(object)$AlignedReadsInBam
    return(new("ChIPprofile",subsetProfile,params=object[[1]]@params))                                                  
  }
  if(method=="signalInRegion"){
    assayList <- assays(object)
    for(k in 1:length(assayList)){
      assayList[[k]] <- t(t(assayList[[k]])/rowSums(assayList[[k]],na.rm=T))
    }
    subsetProfile <- SummarizedExperiment(assayList,rowData=rowData(object))
    exptData(subsetProfile)$names <- exptData(object)$names
    exptData(subsetProfile)$AlignedReadsInBam <- exptData(object)$AlignedReadsInBam
    return(new("ChIPprofile",subsetProfile,params=object[[1]]@params))                                                  
  }  
  }
  

