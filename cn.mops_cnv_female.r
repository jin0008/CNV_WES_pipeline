# Author: Julie BOGOIN
# Modified by: Jinu Han

library(cn.mops)

segments_XY <- read.table(file="/media/hanjinu/PM883/db/refs/interval_list/UKBB/XY.bed",
                    header=FALSE, sep="\t", as.is=TRUE)

female <- read.table(file="female_list.txt", header=FALSE, sep=" ", as.is=TRUE)
for (i in (1:length(female[,1]))){
    female[,1][i] = paste(female[,1][i],'.dedup.bam',sep='')
}


#######################################################################################
print('Working on female...')

setwd(dir=".")

gr <- GRanges(segments_XY[,1], IRanges(segments_XY[,2],segments_XY[,3]))

X <- getSegmentReadCountsFromBAM(female[,1], GR=gr, parallel=64)

resCNMOPS <- exomecn.mops(X, lowerThreshold = -0.9)

resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)

#plot(resCNMOPS, which=5)

setwd(dir="./cn.mops_output")

segm <- as.data.frame(segmentation(resCNMOPS))
print('phase5')
write.csv(segm, file="segmentation_female.csv")

CNVs <- as.data.frame(cnvs(resCNMOPS))
write.csv(CNVs, file="cnvs_female.csv")

CNVRegions <- as.data.frame(cnvr(resCNMOPS))
write.csv(CNVRegions, file="cnvr_female.csv")

setwd(dir="../")
