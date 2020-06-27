# Author: Julie BOGOIN

library(cn.mops)

segments_auto <- read.table(file="/media/jbogoin/Data1/jbogoin/ref/gencode/v34_hg38/autosomes/gencode.v34.basic.annotation.autosome.bed",
                    header=FALSE, sep="\t", as.is=TRUE)

segments_XY <- read.table(file="/media/jbogoin/Data1/jbogoin/ref/gencode/v34_hg38/XY/gencode.v34.basic.annotation.XY.bed",
                    header=FALSE, sep="\t", as.is=TRUE)

all <- list.files(path=".", pattern=".dedup.bam$")

female <- read.table(file="female_list.txt", header=FALSE, sep=" ", as.is=TRUE)
for (i in length(female)){
    female[i]=paste(female[i],'.dedup.bam')
}

male <- read.table(file="male_list.txt", header=FALSE, sep=" ", as.is=TRUE)
for (i in length(male)){
    male[i]=paste(female[i],'.dedup.bam')
}

########################################################################################
print('Working on female...')

gr <- GRanges(segments_XY[,1], IRanges(segments_XY[,2],segments_XY[,3]))

X <- getSegmentReadCountsFromBAM(female, GR=gr)

resCNMOPS <- exomecn.mops(X)
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)

#plot(resCNMOPS, which=5)

setwd(dir="./cn.mops_output")

segm <- as.data.frame(segmentation(resCNMOPS))
write.csv(segm, file="segmentation_female.csv")

CNVs <- as.data.frame(cnvs(resCNMOPS))
write.csv(CNVs, file="cnvs_female.csv")

CNVRegions <- as.data.frame(cnvr(resCNMOPS))
write.csv(CNVRegions, file="cnvr_female.csv")

########################################################################################
print('Working on male...')

gr <- GRanges(segments_XY[,1], IRanges(segments_XY[,2],segments_XY[,3]))

X <- getSegmentReadCountsFromBAM(male, GR=gr)

resCNMOPS <- exomecn.mops(X)
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)

#plot(resCNMOPS, which=5)

segm <- as.data.frame(segmentation(resCNMOPS))
write.csv(segm, file="segmentation_male.csv")

CNVs <- as.data.frame(cnvs(resCNMOPS))
write.csv(CNVs, file="cnvs_female.csv")

CNVRegions <- as.data.frame(cnvr(resCNMOPS))
write.csv(CNVRegions, file="cnvr_male.csv")

########################################################################################
print('Working on all...')

gr <- GRanges(segments_auto[,1], IRanges(segments_auto[,2],segments_auto[,3]))

X <- getSegmentReadCountsFromBAM(all, GR=gr)

resCNMOPS <- exomecn.mops(X)
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)

#plot(resCNMOPS, which=5)

segm <- as.data.frame(segmentation(resCNMOPS))
write.csv(segm, file="segmentation_all.csv")

CNVs <- as.data.frame(cnvs(resCNMOPS))
write.csv(CNVs, file="cnvs_all.csv")

CNVRegions <- as.data.frame(cnvr(resCNMOPS))
write.csv(CNVRegions, file="cnvr_all.csv")
