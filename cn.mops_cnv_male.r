# Author: Julie BOGOIN
# Modified: Jinu Han

library(cn.mops)

segments_XY <- read.table(file="/media/hanjinu/PM883/db/refs/interval_list/UKBB/XY.bed",
                    header=FALSE, sep="\t", as.is=TRUE)


male <- read.table(file="male_list.txt", header=FALSE, sep=" ", as.is=TRUE)
for (i in (1:length(male[,1]))){
    male[,1][i] = paste(male[,1][i],'.dedup.bam',sep='')
}


########################################################################################
print('Working on male...')

setwd(dir=".")

gr <- GRanges(segments_XY[,1], IRanges(segments_XY[,2],segments_XY[,3]))

X <- getSegmentReadCountsFromBAM(male[,1], GR=gr, parallel=64)

resCNMOPS <- exomecn.mops(X, lowerThreshold = -0.9)

resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)

#plot(resCNMOPS, which=5)

setwd(dir="./cn.mops_output")

segm <- as.data.frame(segmentation(resCNMOPS))
write.csv(segm, file="segmentation_male.csv")

CNVs <- as.data.frame(cnvs(resCNMOPS))
write.csv(CNVs, file="cnvs_male.csv")

CNVRegions <- as.data.frame(cnvr(resCNMOPS))
write.csv(CNVRegions, file="cnvr_male.csv")

setwd(dir="../")
