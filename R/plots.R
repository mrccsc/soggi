#' Plot regions
#'
#' A function to plot regions
#' 
#' @usage
#' \S4method{plotRegion}{ChIPprofile}(object,
#' gts,summariseBy,
#' colourBy,lineBy,groupBy,
#' plotregion,outliers)
#'
#'
#' @docType methods
#' @name plotRegion
#' @rdname plotRegion
#' @aliases plotRegion plotRegion,ChIPprofile-method
#' 
#' @author Thomas Carroll
#'
#' @param object A ChIPprofile object 
#' @param gts A list of character vectors and GRanges
#' @param plotregion region to plot
#' @param summariseBy Column names from GRanges elementmetadata. Formula or character vector of column names to use
#' to collapse genomic ranges to summarised profiles.
#' summariseBy can not be used injustion with groups specified by gts argument.
#' @param colourBy Character vector or formula of either column names from colData(object) containing
#' sample metadata Character vector or formula of either column names from colData(object) containing
#' sample metadata or character vector "group" to colour by groups in gts
#' @param lineBy Character vector or formula of either column names from colData(object) containing
#' sample metadata or character vector "group" to set line type by groups in gts
#' @param groupBy Character vector or formula of either column names from colData(object) containing
#' sample metadata or character "group" to colour by groups in gts
#' @param outliers A numeric vector of length 1 containing proportion to exclude from limits 
plotRegion.ChIPprofile <- function(object,gts=NULL,summariseBy=NULL,colourBy=NULL,lineBy=NULL,groupBy=NULL,plotregion="full",outliers=NULL)
{
  
## This function can work using two main grouping logics
## When gts is supplied, groupBy, colourBy, lineBy work only with columns
## from the sample metadata found in colData(object).
## If the gts argument is specified then summariseBy must be single character vector
##  specifies which metadata column which column to be used for grouping
##
## When summariseBy is specified as formula or character vector longer  
##  than one then maintain as/coerced into a formula and gts is ignored.
## groupBy, colourBy, lineBy must be parts of formula as assessed by terms()  
## 
## gts may be named list of character vectors referring to rownames(object) or a granges
## When 
  #app <- lapply(gsets,function(x){colMeans(assays(object)[[1]][rowData(object)$name %in% x,])})
  nOfWindows <- object@params$nOfWindows
  
  ## When running with gts option
  ## SummariseBy can now only be used to select column for gts to be match to
  ## and colourBy, lineBy and groupBy will refer to sample metadata or be group
  
  if(!is.null(gts)){
    
    
    profileList <- list()
    
    for(p in 1:length(assays(object))){
      profileTemp <- assays(object)[[p]]
      if(!is.null(outliers)){
        profileTempList <- lapply(gts,function(x)
          colMeans(winsorizeMatrix(subsetProfile(profileTemp,x,rowData(object),summariseBy),outliers,1-outliers))
          )         
      }else{
        profileTempList <- lapply(gts,function(x)colMeans(subsetProfile(profileTemp,x,rowData(object),summariseBy))) 
      }
      profileMatTemp <- melt(as.data.frame(do.call(cbind,profileTempList)))
      if(object@params$style=="region" & plotregion=="full"){
        axisIndex=c(seq(1,(object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)),
                    (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+seq(1,object@params$nOfWindows)*100,
                    (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+
                      seq(1,(object@params$distanceInRegionEnd+object@params$distanceOutRegionEnd+1)))
      }
      if(object@params$style=="point"){
        axisIndex=c(seq(1,(object@params$distanceAround+object@params$distanceAround+1)))
      }
      if(object@params$style=="percentOfRegion"){
        axisIndex=c(seq(1,((nOfWindows*((object@params$distanceAround)/100))*2)+nOfWindows))
      }      
      profileFrame <-data.frame("xIndex"=axisIndex,Group=profileMatTemp[,1],Sample=basename(unlist(exptData(object)["names"]))[p],Score=profileMatTemp[,2])
      
      profileList[[p]] <- profileFrame
    }  
    meltedProfileFrame <- do.call(rbind,profileList)
    colnames(meltedProfileFrame) <- c("xIndex","Group","Sample","Score")
  }else{
    if(!is.null(outliers)){
      profileList <- lapply(c(assays(object)),function(x)
        colMeans(winsorizeMatrix(x,outliers,1-outliers))
      )         
    }else{    
      profileList <- lapply(c(assays(object)),colMeans)
    }
    profileFrame <- do.call(cbind,profileList)
    colnames(profileFrame) <- basename(unlist(exptData(object)["names"]))
    if(object@params$style=="region"){    
      axisIndex=c(seq(1,(object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)),
                  (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+seq(1,object@params$nOfWindows)*100,
                  (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+
                    seq(1,(object@params$distanceInRegionEnd+object@params$distanceOutRegionEnd+1)))
    }
    if(object@params$style=="point"){
      axisIndex=c(seq(1,(object@params$distanceAround+object@params$distanceAround+1)))
    }
    if(object@params$style=="percentOfRegion"){
      axisIndex=c(seq(1,((object@params$nOfWindows*((object@params$distanceAround)/100))*2)+object@params$nOfWindows))
    }     
    rownames(profileFrame) <- axisIndex
    meltedProfileFrame <- melt(profileFrame)
    colnames(meltedProfileFrame) <- c("xIndex","Sample","Score")
  }
  #profileList <- lapply(c(assays(object)),function(y)lapply(gsets,function(x){colMeans(y[rowData(object)$name %in% x,])}))
  
  
  P <- ggplot(meltedProfileFrame,
              aes(x=xIndex,y=Score))+geom_path(alpha = 1,size=1.3)+xlim(0,max(axisIndex))+ylab("Score")+theme(axis.title.y=element_text(angle=0))
  if(object@params$style=="region" & plotregion =="full"){
    P <- P + scale_x_continuous(breaks=c(1,object@params$distanceOutRegionStart+1,
                                       (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1),
                                       (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(25/100)*(object@params$nOfWindows*100),
                                       (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(75/100)*(object@params$nOfWindows*100),
                                       (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100),
                                       (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+object@params$distanceInRegionEnd+1,
                                       (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+(object@params$distanceInRegionEnd+object@params$distanceOutRegionEnd+1)),
                              labels=c(paste0("Start-",object@params$distanceOutRegionStart),
                                       "Start",
                                       paste0("Start+",object@params$distanceInRegionStart),
                                       "25%",
                                       "75%",
                                       paste0("End-",object@params$distanceInRegionEnd),
                                       "End",
                                       paste0("End+",object@params$distanceOutRegionEnd)
                                       )
                                )+
      theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))
  }
  if(object@params$style=="point"){
    P <- P + scale_x_continuous(breaks=c(1,object@params$distanceAround+1,object@params$distanceAround+1+object@params$distanceAround),
                              labels=c("Centre-1500","Centre","Centre+1500"))+
      theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))
  }
  if(object@params$style=="percentOfRegion"){
    P <- P + scale_x_continuous(breaks=c(1,(nOfWindows*((object@params$distanceAround)/100)),
                                            (nOfWindows*((object@params$distanceAround)/100))+nOfWindows,
                                             2*(nOfWindows*((object@params$distanceAround)/100))+nOfWindows),
                                labels=c(paste0("Start-",(nOfWindows*((object@params$distanceAround)/100)),"%"),
                                                "Start",
                                                "End",
                                                paste0("End+",(nOfWindows*((object@params$distanceAround)/100)),"%")
                                                ))+
      theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))
  }  
  if(object@params$style=="region" & plotregion =="start"){
    P <- P + scale_x_continuous(breaks=c(1,object@params$distanceOutRegionStart+1,object@params$distanceInRegionStart+1+object@params$distanceOutRegionStart),
                              labels=c("Start-1500","Centre","Start+1500"),
                              limits=c(1,object@params$distanceInRegionStart+1+object@params$distanceOutRegionStart))+
      theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))    
  }
  if(object@params$style=="region" & plotregion =="end"){
    P <- P + scale_x_continuous(breaks=(object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+c(1,object@params$distanceInRegionEnd+1,object@params$distanceInRegionEnd+1+object@params$distanceOutRegionEnd),
                              labels=c("End-1500","Centre","End+1500"),
                              limits=(object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+c(1,object@params$distanceInRegionEnd+1+object@params$distanceOutRegionEnd))+
                                theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))    
  }  
  if(!is.null(gts)){
    P <- P+aes(group=Group)
  }
  return(P)
}

winsorizeMatrix <- function(mat,limitlow,limithigh){
  apply(mat,2,function(x)winsorizeVector(x,limitlow,limithigh))
}
winsorizeVector <- function(vect,limitlow,limithigh){
  qs <- quantile(vect,c(limitlow,limithigh))
  vect[vect < qs[1]] <- qs[1]  
  vect[vect > qs[2]] <- qs[2]  
  vect
}
setGeneric("plotRegion", function(object="ChIPprofile",gts=NULL,summariseBy=NULL,colourBy=NULL,lineBy=NULL,groupBy=NULL,plotregion="character",outliers=NULL) standardGeneric("plotRegion"))

#' @rdname plotRegion
#' @export
setMethod("plotRegion", signature(object="ChIPprofile"), plotRegion.ChIPprofile)

subsetProfile <- function(profile,group,granges,summariseColumn){
  print(class(group))
  print(class(group) == "GRanges")  
  if(class(group) == "character"){
    print("")
    if(is.null(summariseColumn)){
      return(profile[rownames(profile) %in% group,])
    }else{
      return(profile[profile[,summariseColumn] %in% group,])
    }
  }
  if(class(group) == "GRanges"){
    print("MakingGRanges subset")
    print(granges %over% group)      
    print(profile[granges %over% group,])          
    return(profile[granges %over% group,])
  }
}