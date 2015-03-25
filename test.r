
ctcf_all_Phastcons <- regionPlot("/Users/tcarroll/Downloads/DP_thymocyte_CTCF_Shih_sorted.bam",allSE,nOfWindows=100,style="percentOfRegion",format="bam",FragmentLength=130,distanceAround = 100,method = "bin")


gts <- c(as.list(unique(names(allSE))),list(unique(names(allSE))[-c(1,2)]))

names(gts) <- c(unique(names(allSE)),"All_SE")
nOfWindows <- 100
pold2 <- plotRegion(ctcf_all_Phastcons,gts=gts)
pold2 <- pold2+facet_wrap(~Group)+aes(colour=Group)

ott <- mm9PC[width(mm9PC) > 1000]



#library(soGGi)
library(GenomicRanges)
mm9Genes <- read.delim("/Users/tcarroll/Downloads/mm9Genes_May2012.txt",sep="\t",h=T)

mm9GeneRanges <- GRanges(seqnames=paste0("chr",mm9Genes[,3]),ranges=IRanges(start=mm9Genes[,1],end=mm9Genes[,2]),strand=mm9Genes[,4],name=mm9Genes[,5],biotype=mm9Genes[,6])
JustChrOfInterest <- unique(as.vector(seqnames(mm9GeneRanges)))[grep("chr\\d.|chr\\d|chrX|chrY|chrM",unique(as.vector(seqnames(mm9GeneRanges))))]
mm9PC <- mm9GeneRanges[mm9GeneRanges$biotype == "protein_coding"]
mm9PC <- mm9PC[order(width(mm9PC),decreasing=T)]
mm9PC <- mm9PC[match(unique(mm9PC$name),mm9PC$name)]
mm9PC <- mm9PC[!mm9PC$name == ""]
mm9PC <- mm9PC[seqnames(mm9PC) %in% JustChrOfInterest]
mm9PC$Feature <- rep("Gene",length(mm9PC))

RNAPII_tam1 <- regionPlot("/Users/tcarroll//Desktop/ziwei//RNAPII-Exp3-4OHTDupMarked.bam",mm9PC,style="region",format="bam")
RNAPII_tam2 <- regionPlot("/Users/tcarroll//Desktop/ziwei//RNAPII-Exp4-4OHTDupMarked.bam",mm9PC,style="region",format="bam")
#RNAPII_tam3 <- regionPlot("/Users/tcarroll//Desktop/ziwei//pol6h_Rep1DupMarked.bam",mm9PC,style="region",format="bam")
#RNAPII_tam4 <- regionPlot("/Users/tcarroll//Desktop/ziwei//pol6h_Rep2DupMarked.bam",mm9PC,style="region",format="bam")

RNAPII_etoh1 <- regionPlot("/Users/tcarroll//Desktop/ziwei//RNAPII-Exp3-ETOHDupMarked.bam",mm9PC,style="region",format="bam")
RNAPII_etoh2 <- regionPlot("/Users/tcarroll//Desktop/ziwei//RNAPII-Exp4-ETOHDupMarked.bam",mm9PC,style="region",format="bam")
#RNAPII_etoh3 <- regionPlot("/Users/tcarroll//Desktop/ziwei//pol0h_Rep1DupMarked.bam",mm9PC,style="region",format="bam")
#RNAPII_etoh4 <- regionPlot("/Users/tcarroll//Desktop/ziwei//pol0h_Rep2DupMarked.bam",mm9PC,style="region",format="bam")

mm9PCTest <- mm9PC[width(mm9PC) > 1000]


RNAPII_tam1 <- regionPlot("/Users/tcarroll//Desktop/ziwei//RNAPII-Exp3-4OHTDupMarked.bam",mm9PCTest,style="percentOfRegion",format="bam",distanceAround = 100)
RNAPII_tam2 <- regionPlot("/Users/tcarroll//Desktop/ziwei//RNAPII-Exp4-4OHTDupMarked.bam",mm9PCTest,style="percentOfRegion",format="bam",distanceAround = 100)
#RNAPII_tam3 <- regionPlot("/Users/tcarroll//Desktop/ziwei//pol6h_Rep1DupMarked.bam",mm9PC,style="region",format="bam")
#RNAPII_tam4 <- regionPlot("/Users/tcarroll//Desktop/ziwei//pol6h_Rep2DupMarked.bam",mm9PC,style="region",format="bam")

RNAPII_etoh1 <- regionPlot("/Users/tcarroll//Desktop/ziwei//RNAPII-Exp3-ETOHDupMarked.bam",mm9PCTest,style="percentOfRegion",format="bam",distanceAround = 100)
RNAPII_etoh2 <- regionPlot("/Users/tcarroll//Desktop/ziwei//RNAPII-Exp4-ETOHDupMarked.bam",mm9PCTest,style="percentOfRegion",format="bam",distanceAround = 100)
#RNAPII_etoh3 <- regionPlot("/Users/tcarroll//Desktop/ziwei//pol0h_Rep1DupMarked.bam",mm9PC,style="region",format="bam")
#RNAPII_etoh4 <- regionPlot("/Users/tcarroll//Desktop/ziwei//pol0h_Rep2DupMarked.bam",mm9PC,style="region",format="bam")
polExp <- c(RNAPII_etoh1,RNAPII_etoh2,RNAPII_tam1,RNAPII_tam2)


gtsDown <- rownames(namedwholeGeneTable[namedwholeGeneTable$log2FoldChange < -1 & namedwholeGeneTable$padj < 0.05 & !is.na(namedwholeGeneTable$padj),])
gtsUp <- rownames(namedwholeGeneTable[namedwholeGeneTable$log2FoldChange > 0.5 & namedwholeGeneTable$padj < 0.05 & !is.na(namedwholeGeneTable$padj),])

setMethod("rbind", "ChIPprofile",
          function (x,...)
          {
            assayList <- lapply(list(x,...),function(x)assays(x)[[1]])
            subsetProfile <- SummarizedExperiment(assayList,rowData=rowData(x))
            exptData(subsetProfile)$names <- unlist(lapply(list(x,...),function(x)exptData(x)$name))
            exptData(subsetProfile)$AlignedReadsInBam <- unlist(lapply(list(x,...),function(x)exptData(x)$AlignedReadsInBam))
            return(new("ChIPprofile",subsetProfile,params=x@params))
            
          }
)

assyCom <- log2(assays(normPolExp)[[3]]+assays(normPolExp)[[4]]) - log2(assays(normPolExp)[[1]]+assays(normPolExp)[[2]])
subsetProfile <- SummarizedExperiment(assyCom,rowData=rowData(normPolExp))
exptData(subsetProfile)$names <- c("Diff")
exptData(subsetProfile)$AlignedReadsInBam <- "Loads"
fid <- (new("ChIPprofile",subsetProfile,params=normPolExp@params))


class1 <- correctedTravelling[correctedTravelling$TSS_log2FoldChange > 0 & 
                                correctedTravelling$Gene_not_TSSlog2FoldChange > 0 & 
                                correctedTravelling$Travellinglog2FoldChange > 0,]

class2 <- correctedTravelling[correctedTravelling$TSS_log2FoldChange > 0 & 
                                correctedTravelling$Gene_not_TSSlog2FoldChange > 0 & 
                                correctedTravelling$Travellinglog2FoldChange < 0,]

class3 <- correctedTravelling[correctedTravelling$TSS_log2FoldChange > 0 & 
                                correctedTravelling$Gene_not_TSSlog2FoldChange < 0 & 
                                correctedTravelling$Travellinglog2FoldChange < 0,]

class4 <- correctedTravelling[correctedTravelling$TSS_log2FoldChange < 0 & 
                                correctedTravelling$Gene_not_TSSlog2FoldChange > 0 & 
                                correctedTravelling$Travellinglog2FoldChange > 0,]

class5 <- correctedTravelling[
  correctedTravelling$Gene_not_TSSlog2FoldChange > 0 & 
    correctedTravelling$Travellinglog2FoldChange > 0,]


class6 <- correctedTravelling[
  correctedTravelling$Gene_not_TSSlog2FoldChange < 0 & 
    correctedTravelling$Travellinglog2FoldChange < 0,]

class7 <- correctedTravelling[
  correctedTravelling$TSS_log2FoldChange < 0 & 
    correctedTravelling$Travellinglog2FoldChange > 0,]

class8 <- correctedTravelling[
  correctedTravelling$TSS_log2FoldChange > 0 & 
    correctedTravelling$Travellinglog2FoldChange < 0,]

class1 <- correctedTravelling[correctedTravelling$TSS_log2FoldChange > 0 & 
                                correctedTravelling$Gene_not_TSSlog2FoldChange > 0 & 
                                correctedTravelling$Travellinglog2FoldChange > 0,]

gts=list(class2=as.vector(class2$Row.names),class3=as.vector(class3$Row.names),class4=as.vector(class4$Row.names),class5=as.vector(class5$Row.names),
         class6=as.vector(class6$Row.names),class7=as.vector(class7$Row.names))


HAranges <- ChIPQC:::GetGRanges("/Users/tcarroll//Downloads//HA_WithInput_Input_peaks.bed")
Endoranges <- ChIPQC:::GetGRanges("/Users/tcarroll//Downloads//Endogenous_WithInput_Input_peaks.bed")

tom <- list(HAranges,Endoranges)
names(tom) <- c("ha","endo")

dnase <- regionPlot("/Users/tcarroll//Downloads/DNAseDupMarkedNormalised.bw",mm9PCTest,style="percentOfRegion",format="bigwig",distanceAround = 100)
pol2 <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/RNAPol2DupMarkedNormalised.bw",mm9PCTest,style="percentOfRegion",format="bigwig",distanceAround = 100)
geneToUse <- rowData(pol2[order(rowSums(assays(pol2)[[1]][,75:125]),decreasing=T)])

pol2ser <- regionPlot("/Users/tcarroll/Downloads/RNAPol2ser2DupMarkedNormalised.bw",mm9PCTest,style="percentOfRegion",format="bigwig",distanceAround = 100)
h3k4me3 <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/H3K4me3_IkNegDupMarkedNormalised.bw",mm9PCTest,style="percentOfRegion",format="bigwig",distanceAround = 100)
h3k9ac <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/H3K9ac_IkNegDupMarkedNormalised.bw",mm9PCTest,style="percentOfRegion",format="bigwig",distanceAround = 100)

geneToUse
pol2_H <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/RNAPol2DupMarkedNormalised.bw",geneToUse[1000:1200],samplename = "highPol",style="percentOfRegion",format="bigwig",distanceAround = 100)
pol2_M <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/RNAPol2DupMarkedNormalised.bw",geneToUse[5000:5200],samplename = "midPol",style="percentOfRegion",format="bigwig",distanceAround = 100)
pol2s_H <- regionPlot("/Users/tcarroll/Downloads/RNAPol2ser2DupMarkedNormalised.bw",geneToUse[1000:1200],samplename = "highPols",style="percentOfRegion",format="bigwig",distanceAround = 100)
pol2s_M <- regionPlot("/Users/tcarroll/Downloads/RNAPol2ser2DupMarkedNormalised.bw",geneToUse[5000:5200],samplename = "midPols",style="percentOfRegion",format="bigwig",distanceAround = 100)
h3k4me3_H <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/H3K4me3_IkNegDupMarkedNormalised.bw",geneToUse[1000:1200],samplename = "highk4me3",style="percentOfRegion",format="bigwig",distanceAround = 100)
h3k4me3_M <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/H3K4me3_IkNegDupMarkedNormalised.bw",geneToUse[5000:5200],samplename = "midk4me3",style="percentOfRegion",format="bigwig",distanceAround = 100)
h3k9ac_H <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/H3K9ac_IkNegDupMarkedNormalised.bw",geneToUse[1000:1200],samplename = "highk9ac",style="percentOfRegion",format="bigwig",distanceAround = 100)
h3k9ac_M <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/H3K9ac_IkNegDupMarkedNormalised.bw",geneToUse[5000:5200],samplename = "midk9ac",style="percentOfRegion",format="bigwig",distanceAround = 100)
dnase_H <- regionPlot("/Users/tcarroll/Downloads/DNAseDupMarkedNormalised.bw",geneToUse[1000:1200],samplename = "highdnase",style="percentOfRegion",format="bigwig",distanceAround = 100)
dnase_M <- regionPlot("/Users/tcarroll/Downloads/DNAseDupMarkedNormalised.bw",geneToUse[5000:5200],samplename = "middnase",style="percentOfRegion",format="bigwig",distanceAround = 100)

pol2s_H <- regionPlot("/Users/tcarroll/Downloads/RNAPol2ser2DupMarkedNormalised.bw",geneToUse[1000:1200],samplename = "highPolser",style="percentOfRegion",format="bigwig",distanceAround = 100)
pol2s_M <- regionPlot("/Users/tcarroll/Downloads/RNAPol2ser2DupMarkedNormalised.bw",geneToUse[5000:5200],samplename = "midPolser",style="percentOfRegion",format="bigwig",distanceAround = 100)
dnase_H <- regionPlot("/Users/tcarroll/Downloads/DNAseDupMarkedNormalised.bw",geneToUse[1000:1200],samplename = "highdnase",style="percentOfRegion",format="bigwig",distanceAround = 100)
dnase_M <- regionPlot("/Users/tcarroll/Downloads/DNAseDupMarkedNormalised.bw",geneToUse[5000:5200],samplename = "middnase",style="percentOfRegion",format="bigwig",distanceAround = 100)

pol <- c(pol2_H,pol2_M,h3k4me3_H,h3k4me3_M,h3k9ac_H,h3k9ac_M)
chipExampleBig <- c(dnase_H,dnase_M,pol2s_H,pol2s_M,pol2_H,pol2_M,h3k4me3_H,h3k4me3_M,h3k9ac_H,h3k9ac_M)
pol2ser <- regionPlot("/Users/tcarroll/Downloads/RNAPol2ser2DupMarkedNormalised.bw",mm9PCTest,style="percentOfRegion",format="bigwig",distanceAround = 100)
h3k4me3 <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/H3K4me3_IkNegDupMarkedNormalised.bw",mm9PCTest,style="percentOfRegion",format="bigwig",distanceAround = 100)
h3k9ac <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/H3K9ac_IkNegDupMarkedNormalised.bw",mm9PCTest,style="percentOfRegion",format="bigwig",distanceAround = 100)

#library(soGGi)
library(GenomicRanges)
mm9Genes <- read.delim("/Users/tcarroll/Downloads/mm9Genes_May2012.txt",sep="\t",h=T)

mm9GeneRanges <- GRanges(seqnames=paste0("chr",mm9Genes[,3]),ranges=IRanges(start=mm9Genes[,1],end=mm9Genes[,2]),strand=mm9Genes[,4],name=mm9Genes[,5],biotype=mm9Genes[,6])
JustChrOfInterest <- unique(as.vector(seqnames(mm9GeneRanges)))[grep("chr\\d.|chr\\d|chrX|chrY|chrM",unique(as.vector(seqnames(mm9GeneRanges))))]
mm9PC <- mm9GeneRanges[mm9GeneRanges$biotype == "protein_coding"]
mm9PC <- mm9PC[order(width(mm9PC),decreasing=T)]
mm9PC <- mm9PC[match(unique(mm9PC$name),mm9PC$name)]
mm9PC <- mm9PC[!mm9PC$name == ""]
mm9PC <- mm9PC[seqnames(mm9PC) %in% JustChrOfInterest]
mm9PC$Feature <- rep("Gene",length(mm9PC))

pol2 <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/RNAPol2DupMarkedNormalised.bw",mm9PCTest,style="percentOfRegion",format="bigwig",distanceAround = 100)
geneToUse <- rowData(pol2[order(rowSums(assays(pol2)[[1]][,75:125]),decreasing=T)])
pol2_H <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/RNAPol2DupMarkedNormalised.bw",geneToUse[1000:1200],samplename = "highPol",style="percentOfRegion",format="bigwig",distanceAround = 100)
pol2_M <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/RNAPol2DupMarkedNormalised.bw",geneToUse[5000:5200],samplename = "midPol",style="percentOfRegion",format="bigwig",distanceAround = 100)
pol2s_H <- regionPlot("/Users/tcarroll/Downloads/RNAPol2ser2DupMarkedNormalised.bw",geneToUse[1000:1200],samplename = "highPols",style="percentOfRegion",format="bigwig",distanceAround = 100)
pol2s_M <- regionPlot("/Users/tcarroll/Downloads/RNAPol2ser2DupMarkedNormalised.bw",geneToUse[5000:5200],samplename = "midPols",style="percentOfRegion",format="bigwig",distanceAround = 100)
h3k4me3_H <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/H3K4me3_IkNegDupMarkedNormalised.bw",geneToUse[1000:1200],samplename = "highk4me3",style="percentOfRegion",format="bigwig",distanceAround = 100)
h3k4me3_M <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/H3K4me3_IkNegDupMarkedNormalised.bw",geneToUse[5000:5200],samplename = "midk4me3",style="percentOfRegion",format="bigwig",distanceAround = 100)
h3k9ac_H <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/H3K9ac_IkNegDupMarkedNormalised.bw",geneToUse[1000:1200],samplename = "highk9ac",style="percentOfRegion",format="bigwig",distanceAround = 100)
h3k9ac_M <- regionPlot("/Users/tcarroll/Downloads/randomTracks-2/H3K9ac_IkNegDupMarkedNormalised.bw",geneToUse[5000:5200],samplename = "midk9ac",style="percentOfRegion",format="bigwig",distanceAround = 100)
dnase_H <- regionPlot("/Users/tcarroll/Downloads/DNAseDupMarkedNormalised.bw",geneToUse[1000:1200],samplename = "highdnase",style="percentOfRegion",format="bigwig",distanceAround = 100)
dnase_M <- regionPlot("/Users/tcarroll/Downloads/DNAseDupMarkedNormalised.bw",geneToUse[5000:5200],samplename = "middnase",style="percentOfRegion",format="bigwig",distanceAround = 100)

tom <- c(dnase,pol2,pol2ser,h3k4me3,h3k9ac)
tom <- log2(zeroToMin2(tom))
p <- plotRegion(tom,outliers=0.01,colourBy="Sample")
