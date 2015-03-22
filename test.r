
ctcf_all_Phastcons <- regionPlot("/Users/tcarroll/Downloads/DP_thymocyte_CTCF_Shih_sorted.bam",allSE,nOfWindows=100,style="percentOfRegion",format="bam",FragmentLength=130,distanceAround = 100,method = "bin")


gts <- c(as.list(unique(names(allSE))),list(unique(names(allSE))[-c(1,2)]))

names(gts) <- c(unique(names(allSE)),"All_SE")
nOfWindows <- 100
pold2 <- plotRegion(ctcf_all_Phastcons,gts=gts)
pold2 <- pold2+facet_wrap(~Group)+aes(colour=Group)

ott <- mm9PC[width(mm9PC) > 1000]
