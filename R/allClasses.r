#' The soggi function and ChIPprofile object.
#'
#' Manual for soggi and ChIPprofile object
#' 
#' @aliases ChIPprofile ChIPprofile-ChIPprofile soggi
#'
#' @references See \url{http://bioinformatics.csc.mrc.ac.uk} for more details on soGGi workflows
#' @rdname ChIPprofile
#' @docType class
#' @param bamFile Character vector for location of BAM file.
#' @param testRanges GRanges object of regions to plot.
#' @param samplename Character vector of sample name. Default is NULL.
#' @param nOfWindows Number of windows to bin regions into for coverage calculations (Default 100)
#' @param FragmentLength Integer vector Predicted or expected fragment length.
#' @param style Point or region (see details)
#' @param distanceAround Distance around centre of region to be used for plotting
#' @param distanceUp Distance upstream from centre of region to be used for plotting
#' @param distanceDown Distance downstream from centre of region to be used for plotting
#' @param distanceInRegionStart Distance into region start 
#' (5' for Watson/positive strand or notspecified strand Regions,3' for Crick/negatie strand regions) 
#' for plotting.
#' @param distanceOutRegionStart Distance out from region start 
#' (5' for Watson/positive strand or notspecified strand Regions,3' for Crick/negatie strand regions) 
#' for plotting.
#' @param distanceInRegionEnd Distance into region end 
#' (3' for Watson/positive strand or notspecified strand Regions,5' for Crick/negatie strand regions) 
#' for plotting.
#' @param distanceOutRegionEnd Distance out from region end 
#' (3' for Watson/positive strand or notspecified strand Regions,5' for Crick/negatie strand regions) 
#' for plotting.
#' @param paired Is data paired end 
#' @param normalize Calculate coverage as RPM. Presently only RPM available.
#' @param plotBy Score to be used for plotting. Presently only coverage.
#' @param removeDup Remove duplicates before calculating coverage.
#' @param verbose TRUE or FALSE
#' @param format BAM or BigWig
#' @param seqlengths Chromosomes to be used. If missing will report all.
#' @param forceFragment Centre fragment and force consistent fragment width.
#' @param method Character vector of value "bp","bin" or "spline". 
#' The bin method divides a region of interest into equal sized bins of number specified in nOfWindows.
#' Coverage or counts are then summarised within these windows.
#' The spline method creates a spline with the number of spline points as specified in nOfWindows argument.
#' @param downSample Down sample GRanges or bamFile to this proportion of orginal.
#' @param genome BSGenome object to be used when using PWM input.
#' @param cutoff Cut-off for idnetifying motifs when using PWM input.
#' @param minFragmentLength Remove fragments smaller than this.
#' @param maxFragmentLength Remove fragments larger than this. 
#' @return ChIPprofile A ChIPprofile object. 
#' @export
setClass("ChIPprofile",contains = "SummarizedExperiment",
         slots=c(params="list"
         ))


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


#' Join, subset and manipulate ChIPprofile objects
#' @rdname manipulateObjects
#' @export
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

#' @rdname manipulateObjects
#' @export
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

#' @rdname manipulateObjects
#' @export
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

#' @rdname manipulateObjects
#' @param j Should be missing
#' @export
setMethod("[[", c("ChIPprofile", "ANY", "missing"),
          function(x, i, j, ...)
          {
            subsetProfile <- SummarizedExperiment(assays(x)[[i, ...]],rowData=rowData(x))
            exptData(subsetProfile)$names <- exptData(x)$names[i]
            exptData(subsetProfile)$AlignedReadsInBam <- exptData(x)$AlignedReadsInBam[i]
            return(new("ChIPprofile",subsetProfile,params=x@params))                        
          })


#' @rdname manipulateObjects
#' @export
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

#' Arithmetic operations
#' @rdname Ops
#' @param e1 ChIPprofile object
#' @param e2 ChIPprofile object
#' @export
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

#' @rdname Ops
#' @export
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

#' @rdname Ops
#' @export
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

#' @rdname Ops
#' @export
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

#' @rdname Ops
#' @export
setMethod("log2", "ChIPprofile",
          function(x){
            assayList <- lapply(assays(x),log2)
            subsetProfile <- SummarizedExperiment(assayList,rowData=rowData(x))
            exptData(subsetProfile)$names <- exptData(x)$names
            exptData(subsetProfile)$AlignedReadsInBam <- exptData(x)$AlignedReadsInBam
            return(new("ChIPprofile",subsetProfile,params=x[[1]]@params))                                                  
            
          }
)

#' @rdname Ops
#' @export
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


#' Normalise ChIPprofiles
#'
#' Various normalisation methods for ChIPprofile objects
#' 
#' @usage
#' \S4method{normalise}{ChIPprofile}(object)
#'
#' @docType methods
#' @name normalise
#' @rdname normalise
#' @aliases normalise normalise,ChIPprofile-method
#' 
#' @author Thomas Carroll
#'
#' @export
#' @param object A ChIPprofile object 
#' @param method Normalisation method 
#' @param normFactors A numeric vector used to scale columns or rows.
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
  

setGeneric("normalise", function(object="ChIPprofile",method="rpm",normFactors=NULL) standardGeneric("normalise"))

#' @rdname normalise
#' @export
setMethod("normalise", signature(object="ChIPprofile",method="character",normFactors="numeric"), normalise.ChIPprofile)

#' Example ChIPprofiles
#'
#' This dataset contains peaks from an in-house EBF1 ChIP-seq 
#'
#' \itemize{
#' \item ChIPprofiles
#' }
#'
#' @docType data
#' @keywords datasets
#' @name chipExampleBig
#' @usage data(chipExampleBig)
#' @return A ChIPprofile object with two rows
NULL