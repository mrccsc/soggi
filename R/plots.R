#' Plot regions
#'
#' A function to plot regions
#' 
#' @usage
#' \S4method{plotRegion}{ChIPprofile}(object)
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
plotRegion.ChIPprofile <- function(object,gts=NULL,plotregion="full")
{
  #app <- lapply(gsets,function(x){colMeans(assays(object)[[1]][rowData(object)$name %in% x,])})
  if(!is.null(gts)){
    profileList <- list()
    for(p in 1:length(assays(object))){
      profileTemp <- assays(object)[[p]]
      profileTempList <- lapply(gts,function(x)colMeans(profileTemp[rowData(object)$name %in% x,])) 
      
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
    profileList <- lapply(c(assays(object)),colMeans)
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
      axisIndex=c(seq(1,((nOfWindows*((object@params$distanceAround)/100))*2)+nOfWindows))
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

setGeneric("plotRegion", function(object="ChIPprofile",gts=NULL,plotregion="character") standardGeneric("plotRegion"))

#' @rdname plotRegion
#' @export
setMethod("plotRegion", signature(object="ChIPprofile"), plotRegion.ChIPprofile)
