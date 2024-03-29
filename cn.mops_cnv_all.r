# Author: Julie BOGOIN
# Modified by: Jinu Han

library(cn.mops)

segments_auto <- read.table(file="/media/hanjinu/PM883/db/refs/interval_list/UKBB/autosomes.bed",
                    header=FALSE, sep="\t", as.is=TRUE)


all <- list.files(path=".", pattern=".dedup.bam$")


#######################################################################################
print('Working on all...')

setwd(dir=".")

gr <- GRanges(segments_auto[,1], IRanges(segments_auto[,2],segments_auto[,3]))

X <- getSegmentReadCountsFromBAM(all, GR=gr, parallel = 64)

resCNMOPS <- exomecn.mops(X)
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)

#plot(resCNMOPS, which=5)

setwd(dir="./cn.mops_output")

segm <- as.data.frame(segmentation(resCNMOPS))
write.csv(segm, file="segmentation_all.csv")

CNVs <- as.data.frame(cnvs(resCNMOPS))
write.csv(CNVs, file="cnvs_all.csv")

CNVRegions <- as.data.frame(cnvr(resCNMOPS))
write.csv(CNVRegions, file="cnvr_all.csv")
