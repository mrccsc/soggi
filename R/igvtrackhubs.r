#' Make sample metadata file for use with IGV.
#'
#' Creates sample metadata file for IGV from data.frames of user supplied metadata and 
#' file locations
#' 
#'
#'
#' @docType methods
#' @name MakeIGVSampleMetadata
#' @rdname MakeIGVSampleMetadata
#' 
#' @author Thomas Carroll
#'
#' @param sampleMetadata A data.frame object with metadata information for samples.
#'  First column must contain unique sample ids. 
#' @param SampleSheet A data.frame of file locations. First column must contain the unique sample ids.
#' @param Directory where igv xml files are located
#' @export
MakeIGVSampleMetadata <- function(sampleMetadata,SampleSheet,igvdirectory){
  write.table("#sampleTable",file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,sep="\t")
  colnames(sampleMetadata)[1] <- "Linking_id"
  print(colnames(sampleMetadata))
  sampleMetadata <- as.matrix(sampleMetadata)
  SampleSheet <- as.matrix(SampleSheet)
  write.table(sampleMetadata,file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=T,quote=F,append=T,sep="\t")
  BamMappings <- cbind(paste(SampleSheet[!is.na(SampleSheet[,"bam"]),"SampleName"],"Bam",sep="_"),SampleSheet[,"SampleName"])
  BigWigMappings <- cbind(paste(SampleSheet[!is.na(SampleSheet[,"bigwig"]),"SampleName"],"Bigwig",sep="_"),SampleSheet[,"SampleName"])
  IntervalMappings <- cbind(paste(SampleSheet[!is.na(SampleSheet[,"interval"]),"SampleName"],"Interval",sep="_"),SampleSheet[,"SampleName"])
  write.table("\n#sampleMapping",file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table("#Bams",file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(BamMappings,file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table("\n#BigWigs",file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(BigWigMappings,file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table("\n#Intervals",file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(IntervalMappings,file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,append=T,sep="\t")
}

#' Make session xml.
#'
#' Creates session for IGV from ChIPQC object
#' 
#'
#'
#' @docType methods
#' @name makeTrackHub
#' @rdname makeTrackHub
#' 
#' @author Thomas Carroll
#'
#' @param QCobject A ChIPQC object 
#' @param peaksDir Directory containing peaks
#' @param bigwigDir Directory containing bigwigs
#' @param IGVdirectory Directory to write IGV xml file
#' @param genome genome for IGV
#' @export
makeTrackHub <- function(QCobject,peaksDir,bigwigDir,IGVdirectory,genome){
  dir.create(IGVdirectory,showWarnings=F)  
  ss <- QCmetadata(QCobject)
  SampleSheet <- ss[,c("ID","Tissue","Factor")]
  colnames(SampleSheet)[1] <- "SampleName"
  fileSheet <- cbind(ss[,c("ID"),drop=F],rep(NA,nrow(ss)))
  peaks <- gsub("/+","/",dir(peaksDir,pattern="*peaks.bed",full.names=T))
  names(peaks) <- gsub("_WithInput.*","",basename(peaks))
  bigwigs <- gsub("/+","/",dir(bigwigDir,pattern="*.Dup.*\\.bw",full.names=T))
  names(bigwigs) <- gsub("DupMarkedNormalised.bw","",basename(bigwigs))
  fileSheet <- merge(fileSheet,cbind(names(peaks),peaks),all=T,by=1,all.x=T,all.y=F)
  fileSheet <- merge(fileSheet,cbind(names(bigwigs),bigwigs),all=T,by=1,all.x=T,all.y=F)
  colnames(SampleSheet) <- c("SampleName",'Tissue',"Factor")
  colnames(fileSheet) <- c("SampleName","bam","bigwig","interval") 
  MakeIGVSampleMetadata(SampleSheet,fileSheet,IGVdirectory)
  return(MakeIGVSessionXML(fileSheet,IGVdirectory,"IGVfull",genome,locusName="All"))
}

#' Make session xml per file
#'
#' Creates session for IGV per file
#' 
#'
#'
#' @docType methods
#' @name MakeIGVSessionXML
#' @rdname MakeIGVSessionXML
#' 
#' @author Thomas Carroll
#'
#' @param SampleSheet A samplesheet containing file locations 
#' @param igvdirectory Directory for IGV
#' @param XMLname Name for IGV session xml
#' @param IGVdirectory Directory to write IGV xml file
#' @param genomeName genome for IGV
#' @param locusName locus to display in igv on loading
#' @export
MakeIGVSessionXML <- function(SampleSheet,igvdirectory,XMLname,genomeName,locusName="All"){
  i <- 1
  require(XML)
  SampleSheet <- as.matrix(SampleSheet)
  Output <- file.path(igvdirectory,paste(XMLname,".xml",sep=""))
  GlobalNode <- newXMLNode("Global",attrs=c(genome.value=genomeName,groupTracksBy="Linking_id",locus=locusName,version=3))
  ResourcesNode <- newXMLNode("Resources",parent=GlobalNode)
  MetaDataNode <- newXMLNode("Resource",parent=ResourcesNode,attrs=c(name="SampleMetadata",path=relativePath(file.path(igvdirectory,"SampleMetadata.txt"),Output),relativePath=T))
  PanelDataNode <-  newXMLNode("Panel",attrs=c(height="350",name="DataPanel",width="1115"),parent=GlobalNode)
  #bamFiles <- SampleSheet[!is.na(SampleSheet[,"bam"]),"bam"]
  #bigwigFiles <- SampleSheet[!is.na(SampleSheet[,"bigwig"]),"bigwig"]
  #intervals <- SampleSheet[!is.na(SampleSheet[,"interval"]),"interval"]
  bamFiles <- SampleSheet[,"bam"]
  bigwigFiles <- SampleSheet[,"bigwig"]
  intervalFiles <- SampleSheet[,"interval"]    
  resources <- vector("list")
  #print(Output)
  for(i in 1:nrow(SampleSheet)){
    print(i)
    if(!is.na(SampleSheet[i,"bam"])){
      NewName <- paste(SampleSheet[i,"SampleName"],"_Bam",sep="")
      resources <-  c(resources,list(newXMLNode("Resource",parent=ResourcesNode,attrs=c(label=NewName,name=NewName,path=relativePath(bamFiles[i],Output),relativePath=T))))
      TrackNode <-  newXMLNode("Track",attrs=c(altColor="0,0,178",color="0,0,178",colorOption="UNEXPECTED_PAIR",displayMode="EXPANDED",featureVisibilityWindow="-1",fontSize="10",id=relativePath(bamFiles[i],Output),name=NewName,showDataRange="true",sortByTag="",visible="true"),parent=PanelDataNode)
    }
    if(!is.na(SampleSheet[i,"interval"])){
      NewName <- paste(SampleSheet[i,"SampleName"],"_Interval",sep="")
      resources <-  c(resources,list(newXMLNode("Resource",parent=ResourcesNode,attrs=c(label=NewName,name=NewName,path=relativePath(intervalFiles[i],Output),relativePath=T))))
      TrackNode <-  newXMLNode("Track",attrs=c(altColor="0,0,178",color="0,0,178",displayMode="COLLAPSED",featureVisibilityWindow="-1",fontSize="10",height="45",id=relativePath(intervalFiles[i],Output),name=NewName,renderer="BASIC_FEATURE",showDataRange="true",sortable="false",visible="true",windowFunction="count"),parent=PanelDataNode)
    }
    if(!is.na(SampleSheet[i,"bigwig"])){
      NewName <- paste(SampleSheet[i,"SampleName"],"_Bigwig",sep="")
      print(relativePath(bigwigFiles[i],Output))
      resources <-  c(resources,list(newXMLNode("Resource",parent=ResourcesNode,attrs=c(label=NewName,name=NewName,path=relativePath(bigwigFiles[i],Output),relativePath=T))))
      TrackNode <-  newXMLNode("Track",attrs=c(altColor="0,0,178",autoscale="true",color="0,0,178",displayMode="COLLAPSED",featureVisibilityWindow="-1",fontSize="10",id=relativePath(bigwigFiles[i],Output),name=NewName,renderer="BAR_CHART",showDataRange="true",visible="true",windowFunction="mean"),parent=PanelDataNode)
      DisplayRangeNode <-  newXMLNode("DataRange",attrs=c(baseline="0.0",drawBaseline="true",flipAxis="false",maximum="50",minimum="5",type="LINEAR"),parent=TrackNode)
    }
  }  
  saveXML(GlobalNode,file=Output)
  
  return(Output)
}

#' Make HTML page for IGV sessions
#'
#' Creates HTML page for IGV sessions
#' 
#'
#'
#' @docType methods
#' @name maketracktable
#' @rdname maketracktable
#' 
#' @author Thomas Carroll
#'
#' @param fileSheet A data frame, matrix or character containing sample information locations
#' @param SampleSheet A data frame, matrix or character containing sample information locations 
#' @param igvdirectory Directory for IGV
#' @param filename Name for IGV session xml
#' @param basedirectory Directory to write IGV xml file
#' @param genome genome for IGV
#' @export
maketracktable <- function(fileSheet,SampleSheet,filename,basedirectory,genome){
  
  basedirectory <- gsub("/$","",basedirectory)
  MakeIGVSampleMetadata(SampleSheet,fileSheet,basedirectory)
  
  #genome <- "mm9"
  xmlFiles <- unlist(lapply(seq(1,nrow(fileSheet)),function(x)
    MakeIGVSessionXML(fileSheet[x,,drop=F],
                      basedirectory,
                      paste0(fileSheet[x,1],"igv"),
                      genome,
                      locusName="All")
  ))
  
  dataTableJS <- readLines(system.file(package="tracktables","js","datatables.js"))
  jqueryJS <- readLines(system.file(package="tracktables","js","jquery.min.js"))
  dataTableCSS <- readLines(system.file(package="tracktables","js","jquery.datatables.css"))
  dataTableScroller <- readLines(system.file(package="tracktables","js","dataTables.scroller.min.js"))
  tracktablesCSS <- readLines(system.file(package="tracktables","js","tracktables.css"))
  
  giHTMLs <- vector("character",nrow(fileSheet))
  giHTMLLinks <- vector("character",nrow(fileSheet))
  for(l in 1:nrow(fileSheet)){
    if(!is.na(fileSheet[l,"interval"])){
       giHTMLs[l] <- makebedtable(ChIPQC:::GetGRanges(fileSheet[l,"interval"]),paste0(fileSheet[l,"SampleName"],"GI.html"),basedirectory)  
       giHTMLLinks[l] <- paste0("\"<a href=\\\"",file.path(basedirectory,basename(giHTMLs[l])),"\\\">Intervals</a>\"")
    }else{
      giHTMLLinks[l] <- shQuote("No Intervals")
      
    }
  }
  
  library(RJSONIO)
  files <- unlist(lapply(xmlFiles,function(x)relativePath(x,
                                                          gsub("//","/",file.path(basedirectory,filename))
  )))
  t3mp <- "\"<a href=\\\"http://localhost:60151/load?file=\".concat(dir.concat(\"/"
  t4mp <- "\\\"\".concat(\""
  t5mp <- "</a>\")))"
  jsMat <- cbind(
    matrix(paste0("\"",as.vector(SampleSheet),"\""),ncol=ncol(SampleSheet),byrow=F),
    paste0(t3mp,files,"&merge=true",t4mp,">",SampleSheet[,1],t5mp),
    giHTMLLinks
  )
  setigv <- paste0("var igvtable = [",paste0(
    "[",apply(jsMat,1,function(x)paste0(
      x,collapse=","))
    ,"]\n",collapse=",")
    ,"];",sep="")
  
  jspart1 <- paste0("var loc = window.location.pathname;\n",
                    "var dir = loc.substring(0, loc.lastIndexOf('/'));\n",setigv,"\n")
  jspart2 <- paste0(
    "$(document).ready(function() {
    $('#demo').html( '<table cellpadding=\"0\" cellspacing=\"0\" border=\"0\" class=\"display\" id=\"example\"></table>' );
    $('#example').dataTable( {
    \"data\": igvtable,\ncolumns:",
    paste0("[",paste0(
      unlist(lapply(c(colnames(SampleSheet),"IGV","Intervals"),function(x)paste0(
        c("{\"title\"",paste0(
          "\"",
          x,"\"}")
        ),collapse=":")
      )),collapse=",\n")
      ,"]")
    ,"\n","} );\n","} );\n")
  jspart1.2 <- paste0(jspart1,jspart2)
  doc <- newXMLDoc(isHTML = T)
  html <- newXMLNode("html",parent=doc)
  head <- newXMLNode("head",parent = html)
  title <- newXMLNode("h2",
                      "Tracktables Report",
                      parent=head)
  css <- newXMLNode("style",
                    attrs=c("style type"="text/css","class"="init"),
                    paste0(dataTableCSS,collapse=""),
                    parent=head)
  tracktablescss <- newXMLNode("style",
                    attrs=c("style type"="text/css","class"="init"),
                    paste0(tracktablesCSS,collapse=""),
                    parent=head)  
  jqueryjs <- newXMLNode("script",
                         attrs=c(type="text/javascript",language="javascript"),
                         paste0(jqueryJS,collapse=""),
                         parent=head)
  datatablejs <- newXMLNode("script",
                            attrs=c(type="text/javascript",language="javascript"),
                            paste0(dataTableJS,collapse=""),
                            parent=head)
  jspart1.2js <- newXMLNode("script",
                            attrs=c(type="text/javascript",language="javascript"),
                            jspart1.2,
                            parent=head)
  body <- newXMLNode("body",
                     attrs=c(class="dt-example"),
                     parent=html)
  div <- newXMLNode("div",
                    attrs=c(class="container"),
                    parent=body)
  section <- newXMLNode("section",
                        parent=div)
  divtttext <- newXMLNode("div",
                          attrs=c(id="tttext"),
                          parent=section)
  h1 <- newXMLNode("h1","The Tracktables Sample Report",
                   parent=divtttext)
  p1 <- newXMLNode("p","
                   This report contains sample information and dynamic links to display and control Broad's Integrative Genome Browser (IGV). This report aims to speed up the organisation and visualisation of genomics data by allowing for the passing of metadata and sample information to IGV and the rapid selection of samples and points of interest using HTML tables.
                   ",
                   parent=divtttext)
  p2 <- newXMLNode("p","Getting started:",
                   parent=divtttext)
  ul1 <- newXMLNode("ul","",
                   parent=divtttext)
  li1 <- newXMLNode("li","To take advantage of the integration with IGV, <b>IGV must be already running </b>on your machine or can be launched now from this webstart.",
                    parent=ul1,cdata=TRUE)
  li2 <- newXMLNode("li","To load coverage, BAM and/or interval files (bed, narrow peak format etc) simply click the respective sample link in the IGV column.",
                    parent=ul1)
  li3 <- newXMLNode("li","To open a new tracktable containing information on Sample interval files click the link in that sample's repsective Intervals column.",
                    parent=ul1)
  p3 <- newXMLNode("p","For further information on the use of tracktables, please see our Github page or Bioconductor site.",
                   parent=divtttext)
  div2 <- newXMLNode("div",
                     attrs=c(id="demo"),
                     parent=section)
  saveXML(doc,file=file.path(basedirectory,filename),doctype="html")
  return(doc)
}

#' Make HTML page for interval files or GRanges
#'
#' Creates HTML page for interval files or GRanges
#' 
#'
#'
#' @docType methods
#' @name makebedtable
#' @rdname makebedtable
#' 
#' @author Thomas Carroll
#'
#' @param grangesObject A GRanges or object which may be parsed by ChIPQC GetGRanges method.
#' @param name Name for IGV session xml
#' @param basedirectory Directory to write IGV xml file
#' @export
makebedtable <- function(grangesObject,name,basedirectory){
  
  dataTableJS <- readLines(system.file(package="tracktables","js","datatables.js"))
  jqueryJS <- readLines(system.file(package="tracktables","js","jquery.min.js"))
  dataTableCSS <- readLines(system.file(package="tracktables","js","jquery.datatables.css"))
  dataTableScroller <- readLines(system.file(package="tracktables","js","dataTables.scroller.min.js"))
  tracktablesCSS <- readLines(system.file(package="tracktables","js","tracktables.css"))
  
  grangesFrame <- as.matrix(as.data.frame(grangesObject))
  grangesFrame <- str_trim(grangesFrame)
  jsarray <- paste("[",paste0("[",apply(grangesFrame,1,function(x)paste0(c(shQuote(c(paste0("<a href=\"http://localhost:60151/goto?locus=",x[1],":",x[2],"-",x[3],"\">IGV</a>"))),shQuote(x)),collapse=",")),"]",collapse=",\n"),"]")
  jsArrayForIGV <- paste0("var igvtable =",jsarray,";\n")
  jspart2 <- paste0(
    "$(document).ready(function() {
    $('#demo').html( '<table cellpadding=\"0\" cellspacing=\"0\" border=\"0\" class=\"display\" id=\"example\"></table>' );
    $('#example').dataTable( {
    deferRender:    true,
    dom:            \"frtiS\",
    scrollY:        200,
    scrollCollapse: true,
    
    \"data\": igvtable,\ncolumns:",
    paste0("[",paste0(
      unlist(lapply(c("IGV_Link",colnames(as.data.frame(grangesObject))),function(x)paste0(
        c("{\"title\"",paste0(
          "\"",
          x,"\"}")
        ),collapse=":")
      )),collapse=",\n")
      ,"]")
    ,"\n","} );\n","} );\n")
  
  jspart1.2 <- paste0(jsArrayForIGV,jspart2)
  doc <- newXMLDoc(isHTML = T)
  html <- newXMLNode("html",parent=doc)
  head <- newXMLNode("head",parent = html)
  title <- newXMLNode("h2",
                      "Tracktables Report",
                      parent=head)
  css <- newXMLNode("style",
                    attrs=c("style type"="text/css","class"="init"),
                    paste0(dataTableCSS,collapse=""),
                    parent=head)
  tracktablescss <- newXMLNode("style",
                               attrs=c("style type"="text/css","class"="init"),
                               paste0(tracktablesCSS,collapse=""),
                               parent=head)  
  jqueryjs <- newXMLNode("script",
                         attrs=c(type="text/javascript",language="javascript"),
                         paste0(jqueryJS,collapse=""),
                         parent=head)
  datatablejs <- newXMLNode("script",
                            attrs=c(type="text/javascript",language="javascript"),
                            paste0(dataTableJS,collapse=""),
                            parent=head)
  datatableScroller <- newXMLNode("script",
                                  attrs=c(type="text/javascript",language="javascript"),
                                  paste0(dataTableScroller,collapse=""),
                                  parent=head)  
  jspart1.2js <- newXMLNode("script",
                            attrs=c(type="text/javascript",language="javascript"),
                            jspart1.2,
                            parent=head)
  body <- newXMLNode("body",
                     attrs=c(class="dt-example"),
                     parent=html)
  div <- newXMLNode("div",
                    attrs=c(class="container"),
                    parent=body)
  section <- newXMLNode("section",
                        parent=div)
  divtttext <- newXMLNode("div",
                     attrs=c(id="tttext"),
                     parent=section)
  h1 <- newXMLNode("h1","The Tracktables Interval Report",
                   parent=divtttext)
  p1 <- newXMLNode("p","
                   This report contains genomic interval coordinates,  metadata and dynamic links to control the region displayed within Broad's Integrative Genome Browser (IGV). This alows fort rapid visualisation and interrogation of points of interest within the Genome Browser using HTML tables.
                   ",
                   parent=divtttext)
  p2 <- newXMLNode("p","Getting started:",
                   parent=divtttext)
  ul1 <- newXMLNode("ul","",
                    parent=divtttext)
  li1 <- newXMLNode("li",paste0("To take advantage of the integration with IGV, ","IGV must be already running","on your machine or can be launched now from this webstart."),
                    parent=ul1)
  li2 <- newXMLNode("li","To change IGV display to the region of interest, simply click the respective Interval link in the IGV column.",
                    parent=ul1)
  p3 <- newXMLNode("p","For further information on the use of tracktables, please see our Github page or Bioconductor site.",
                   parent=divtttext)  
  div2 <- newXMLNode("div",
                     attrs=c(id="demo"),
                     parent=section)
  saveXML(doc,file=file.path(basedirectory,name),doctype="html")
  
  
}

