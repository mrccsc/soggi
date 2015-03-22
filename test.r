
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
