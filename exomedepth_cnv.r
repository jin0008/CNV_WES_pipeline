# Author: Julie BOGOIN

library(ExomeDepth)
library(seqinr)

targets_auto <- read.table(file="/media/Data1/jbogoin/ref/gencode/v34_hg38/autosomes/gencode.v34.basic.annotation.autosome.bed",
                    header=FALSE, sep=" ", as.is=TRUE)

targets_XY <- read.table(file="/media//Data1/jbogoin/ref/gencode/v34_hg38/XY/gencode.v34.basic.annotation.XY.scratch.bed",
                    header=FALSE, sep="\t", as.is=TRUE)

female <- read.table(file="female_list.txt", header=FALSE, sep=" ", as.is=TRUE)
for (i in (1:length(female[,1]))){
    female[,1][i] = paste(female[,1][i],'.dedup.bam',sep='')
}
female_vec <-  unlist(female, recursive = TRUE, use.names = TRUE)

male <- read.table(file="male_list.txt", header=FALSE, sep=" ", as.is=TRUE)
for (i in (1:length(male[,1]))){
    male[,1][i] = paste(male[,1][i],'.dedup.bam',sep='')
}
male_vec <-  unlist(male, recursive = TRUE, use.names = TRUE)

bams_list <- list.files(path=".", pattern=".dedup.bam$")
bams_vec <-  unlist(bams_list, recursive = TRUE, use.names = TRUE)



########################################################################################
print('Working on female...')

# Generate read count data
ExomeCount <- getBamCounts(bed.frame = targets_XY,
                bam.files = female_vec,
                include.chr = FALSE,
                referenceFasta = NULL) #format de sortie: GRanges class

## Remove the annoying X
names(ExomeCount) <- sub("X", "", names(ExomeCount), fixed=TRUE)

# Conversion en dataframe
ExomeCount.dafr <- as(ExomeCount, 'data.frame')

## Remove the annoying X
names(ExomeCount.dafr) <- sub("X", "", names(ExomeCount.dafr), fixed=TRUE)

## Remove the annoying chr letters
ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$chromosome), pattern = 'chr', replacement = '')
ExomeCount.dafr$names <- rep('CDS', length(ExomeCount.dafr$chromosome))

### Prepare the main matrix of read count data
sample_names = grep(names(ExomeCount.dafr), pattern = '*.dedup.bam')
ExomeCount.mat <- as.matrix(ExomeCount.dafr[, sample_names])

nsamples <- ncol(ExomeCount.mat)

setwd(dir="./exomedepth_output/female")

### start looping over each sample
for (i in 1:nsamples) {

	#### Create the aggregate reference set for this sample
	my.choice <- select.reference.set (test.counts = ExomeCount.mat[,i],
					   reference.counts = ExomeCount.mat[,-i],
					   bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
					   n.bins.reduced = 10000)

	my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE], MAR = 1, FUN = sum)

	message('Now creating the ExomeDepth object')
	all.exons <- new('ExomeDepth', test = ExomeCount.mat[,i], reference = my.reference.selected, formula = 'cbind(test, reference) ~ 1')

	################ Now call the CNVs
	all.exons <- CallCNVs(x = all.exons, transition.probability = 10^-4,
			      chromosome = ExomeCount.dafr$chromosome,
			      start = ExomeCount.dafr$start,
			      end = ExomeCount.dafr$end,
			      name =  ExomeCount.dafr$names)

	# Creating CSV files
    col_names <- colnames(ExomeCount.mat)
	col_names <- sub("dedup.bam", "", col_names, fixed=TRUE)
	sample_name <-  col_names[i]
	output.file <- paste(sample_name, 'csv', sep = '')
	write.csv(file = output.file, x = all.exons@CNV.calls, row.names = FALSE)

}


# #######################################################################################
# print('Working on male...')

setwd(dir="../..")

# Generate read count data
ExomeCount <- getBamCounts(bed.frame = targets_XY,
                bam.files = male_vec,
                include.chr = FALSE,
                referenceFasta = NULL) #format de sortie: GRanges class

## Remove the annoying X
names(ExomeCount) <- sub("X", "", names(ExomeCount), fixed=TRUE)

# Conversion en dataframe
ExomeCount.dafr <- as(ExomeCount, 'data.frame')

## Remove the annoying X
names(ExomeCount.dafr) <- sub("X", "", names(ExomeCount.dafr), fixed=TRUE)

## Remove the annoying chr letters
ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$chromosome), pattern = 'chr', replacement = '')
ExomeCount.dafr$names <- rep('CDS', length(ExomeCount.dafr$chromosome))

### Prepare the main matrix of read count data
sample_names = grep(names(ExomeCount.dafr), pattern = '*.dedup.bam')
ExomeCount.mat <- as.matrix(ExomeCount.dafr[, sample_names])

### start looping over each sample
nsamples <- ncol(ExomeCount.mat)

### start looping over each sample
setwd(dir="./exomedepth_output/male")

for (i in 1:nsamples) {

	#### Create the aggregate reference set for this sample
	my.choice <- select.reference.set (test.counts = ExomeCount.mat[,i],
					   reference.counts = ExomeCount.mat[,-i],
					   bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
					   n.bins.reduced = 10000)

	my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE], MAR = 1, FUN = sum)

	message('Now creating the ExomeDepth object')
	all.exons <- new('ExomeDepth', test = ExomeCount.mat[,i], reference = my.reference.selected, formula = 'cbind(test, reference) ~ 1')

	################ Now call the CNVs
	all.exons <- CallCNVs(x = all.exons, transition.probability = 10^-4,
			      chromosome = ExomeCount.dafr$chromosome,
			      start = ExomeCount.dafr$start,
			      end = ExomeCount.dafr$end,
			      name =  ExomeCount.dafr$names)

	# Creating CSV files
    col_names <- colnames(ExomeCount.mat)
	col_names <- sub("dedup.bam", "", col_names, fixed=TRUE)
	sample_name <-  col_names[i]
	output.file <- paste(sample_name, 'csv', sep = '')
	write.csv(file = output.file, x = all.exons@CNV.calls, row.names = FALSE)

}


########################################################################################
print('Working on all...')
#### get the annotation datasets to be used later
#targets.GRanges <- GRanges(seqnames = targets[,1], IRanges(start=targets[,2], end=targets[,3]))

setwd(dir="../..")

# Generate read count data
ExomeCount <- getBamCounts(bed.frame = targets_auto,
                bam.files = bams_vec,
                include.chr = FALSE,
                referenceFasta = NULL) #format de sortie: GRanges class

## Remove the annoying X
names(ExomeCount) <- sub("X", "", names(ExomeCount), fixed=TRUE)

# Conversion en dataframe
ExomeCount.dafr <- as(ExomeCount, 'data.frame')

## Remove the annoying X
names(ExomeCount.dafr) <- sub("X", "", names(ExomeCount.dafr), fixed=TRUE)

## Remove the annoying chr letters
ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$chromosome), pattern = 'chr', replacement = '')
ExomeCount.dafr$names <- rep('CDS', length(ExomeCount.dafr$chromosome))

### Prepare the main matrix of read count data
sample_names = grep(names(ExomeCount.dafr), pattern = '*.dedup.bam')
ExomeCount.mat <- as.matrix(ExomeCount.dafr[, sample_names])

nsamples <- ncol(ExomeCount.mat)

### start looping over each sample
setwd(dir="./exomedepth_output/all")

for (i in 1:nsamples) {

	#### Create the aggregate reference set for this sample
	my.choice <- select.reference.set (test.counts = ExomeCount.mat[,i],
					   reference.counts = ExomeCount.mat[,-i],
					   bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
					   n.bins.reduced = 10000)

	my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE], MAR = 1, FUN = sum)

	message('Now creating the ExomeDepth object')
	all.exons <- new('ExomeDepth', test = ExomeCount.mat[,i], reference = my.reference.selected, formula = 'cbind(test, reference) ~ 1')

	################ Now call the CNVs
	all.exons <- CallCNVs(x = all.exons, transition.probability = 10^-4,
			      chromosome = ExomeCount.dafr$chromosome,
			      start = ExomeCount.dafr$start,
			      end = ExomeCount.dafr$end,
			      name =  ExomeCount.dafr$names)

	# Creating CSV files
    col_names <- colnames(ExomeCount.mat)
	col_names <- sub("dedup.bam", "", col_names, fixed=TRUE)
	sample_name <-  col_names[i]
	output.file <- paste(sample_name, 'csv', sep = '')
	write.csv(file = output.file, x = all.exons@CNV.calls, row.names = FALSE)

}
	#### Create the aggregate reference set for this sample
	my.choice <- select.reference.set (test.counts = ExomeCount.mat[,i],
					   reference.counts = ExomeCount.mat[,-i],
					   bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
					   n.bins.reduced = 10000)

	my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE], MAR = 1, FUN = sum)

	message('Now creating the ExomeDepth object')
	all.exons <- new('ExomeDepth', test = ExomeCount.mat[,i], reference = my.reference.selected, formula = 'cbind(test, reference) ~ 1')

	################ Now call the CNVs
	all.exons <- CallCNVs(x = all.exons, transition.probability = 10^-4,
			      chromosome = ExomeCount.dafr$chromosome,
			      start = ExomeCount.dafr$start,
			      end = ExomeCount.dafr$end,
			      name =  ExomeCount.dafr$names)

	# Creating CSV files
    col_names <- colnames(ExomeCount.mat)
	col_names <- sub("dedup.bam", "", col_names, fixed=TRUE)
	sample_name <-  col_names[i]
	output.file <- paste(sample_name, 'csv', sep = '')
	write.csv(file = output.file, x = all.exons@CNV.calls, row.names = FALSE)

}


########################################################################################
print('Working on all...')
#### get the annotation datasets to be used later
#targets.GRanges <- GRanges(seqnames = targets[,1], IRanges(start=targets[,2], end=targets[,3]))

setwd(dir="../..")

# Generate read count data
ExomeCount <- getBamCounts(bed.frame = targets_auto,
                bam.files = bams_vec,
                include.chr = FALSE,
                referenceFasta = NULL) #format de sortie: GRanges class

## Remove the annoying X
names(ExomeCount) <- sub("X", "", names(ExomeCount), fixed=TRUE)

# Conversion en dataframe
ExomeCount.dafr <- as(ExomeCount, 'data.frame')

## Remove the annoying X
names(ExomeCount.dafr) <- sub("X", "", names(ExomeCount.dafr), fixed=TRUE)

## Remove the annoying chr letters
ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$chromosome), pattern = 'chr', replacement = '')
ExomeCount.dafr$names <- rep('CDS', length(ExomeCount.dafr$chromosome))

print(head(ExomeCount.dafr))

### Prepare the main matrix of read count data
sample_names = grep(names(ExomeCount.dafr), pattern = '*.dedup.bam')
ExomeCount.mat <- as.matrix(ExomeCount.dafr[, sample_names])

nsamples <- ncol(ExomeCount.mat)

### start looping over each sample
setwd(dir="./exomedepth_output/all")

for (i in 1:nsamples) {

	#### Create the aggregate reference set for this sample
	my.choice <- select.reference.set (test.counts = ExomeCount.mat[,i],
					   reference.counts = ExomeCount.mat[,-i],
					   bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
					   n.bins.reduced = 10000)

	my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE], MAR = 1, FUN = sum)

	message('Now creating the ExomeDepth object')
	all.exons <- new('ExomeDepth', test = ExomeCount.mat[,i], reference = my.reference.selected, formula = 'cbind(test, reference) ~ 1')

	################ Now call the CNVs
	all.exons <- CallCNVs(x = all.exons, transition.probability = 10^-4,
			      chromosome = ExomeCount.dafr$chromosome,
			      start = ExomeCount.dafr$start,
			      end = ExomeCount.dafr$end,
			      name =  ExomeCount.dafr$names)

	# Creating CSV files
    col_names <- colnames(ExomeCount.mat)
	col_names <- sub("dedup.bam", "", col_names, fixed=TRUE)
	sample_name <-  col_names[i]
	output.file <- paste(sample_name, 'csv', sep = '')
	write.csv(file = output.file, x = all.exons@CNV.calls, row.names = FALSE)

}
