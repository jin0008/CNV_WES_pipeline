# Author: Julie BOGOIN
# Modified by: Jinu Han

source ~/miniconda3/etc/profile.d/conda.sh
conda activate gatk_env

REF="/media/hanjinu/PM883/db/refs/hg38_broad/Homo_sapiens_assembly38.fasta"
DIC="/media/hanjinu/PM883/db/refs/hg38_broad/Homo_sapiens_assembly38.dict"

CENTROMETIC_AUTO="/media/hanjinu/PM883/db/refs/hg38_centromeric/centromeric_regions_autosomes.bed"
CENTROMETIC_XY="/media/hanjinu/PM883/db/refs/hg38_centromeric/centromeric_regions_XY.bed"

TARGET_AUTO="/media/hanjinu/PM883/db/refs/interval_list/UKBB/autosomes.interval_list"
TARGET_XY="/media/hanjinu/PM883/db/refs/interval_list/UKBB/XY.interval_list"

PLOIDY_AUTO="/media/hanjinu/PM883/db/refs/contig_ploidy_priors/contig_ploidy_prior_autosomes.tsv"
PLOIDY_XY="/media/hanjinu/PM883/db/refs/contig_ploidy_priors/contig_ploidy_prior_XY.tsv"

MAPPABILITY="/media/hanjinu/PM883/db/refs/mappability/k100.umap.bed.gz"
SEGMENTAL_DUPLICATE="/media/hanjinu/PM883/db/refs/segmental_duplicate/segmental_duplicate.bed.gz"

XY_Preprocessed_Interval="/media/hanjinu/PM883/db/refs/interval_list/UKBB/targets.XY.preprocessed.interval_list"
AUTOSOME_Preprocessed_Interval="/media/hanjinu/PM883/db/refs/interval_list/UKBB/targets.autosomes.preprocessed.interval_list"

echo ""
echo "*************************"
echo "GATK4 CNV DETECTION start"
echo "*************************"
echo ""

rm -rf gatkcnv_output
mkdir gatkcnv_output
cd gatkcnv_output
mkdir all
mkdir female
mkdir male
cd ..

########################################################################################
echo ""
echo "Working on female..."
echo ""

FEMALE=""
while read line
do
FEMALE+="$line.dedup.bam ";
done < female_list.txt

echo "list female:"
echo $FEMALE
echo ""

# Define the resolution of the analysis with a genomic intervals list
#gatk PreprocessIntervals \
#    -R $REF \
#    -L $TARGET_XY \
#    -XL $CENTROMETIC_XY \
#    --bin-length 0 \
#    --padding 50 \
#    --interval-merging-rule OVERLAPPING_ONLY \
#    -O gatkcnv_output/female/targets.preprocessed.interval_list \
#    --verbosity ERROR

for sample_id in $FEMALE;
do SAMPLE=${sample_id%%.dedup.bam}; \

# Collect raw integer counts data
gatk CollectReadCounts \
        -L $XY_Preprocessed_Interval \
        -XL $CENTROMETIC_XY \
        -R $REF \
	--tmp-dir $PWD \
        --interval-merging-rule OVERLAPPING_ONLY \
        -I $SAMPLE.dedup.bam \
        --format TSV \
        -O gatkcnv_output/female/$SAMPLE.tsv \
        --verbosity ERROR;

done

# AnnotateIntervals with GC content
gatk AnnotateIntervals \
    -L $XY_Preprocessed_Interval \
    -XL $CENTROMETIC_XY \
    -R $REF \
    --tmp-dir $PWD \
    -imr OVERLAPPING_ONLY \
    -O gatkcnv_output/female/targets.annotated.tsv \
    --mappability-track $MAPPABILITY \
    --segmental-duplication-track $SEGMENTAL_DUPLICATE \
    --verbosity ERROR


cd gatkcnv_output/female

COUNTS_LIST=""
for counts_file in *.tsv; do
        if [ $(grep -v targets.annotated.tsv <<<$counts_file) ]; then
        FILE=${counts_file%%.tsv};
        COUNTS_LIST+="-I $FILE.tsv ";
        fi
done

# FilterIntervals based on GC-content and cohort extreme counts
gatk FilterIntervals \
        -L $XY_Preprocessed_Interval \
        -XL $CENTROMETIC_XY \
	--tmp-dir $PWD \
        --annotated-intervals targets.annotated.tsv \
        $COUNTS_LIST \
        -imr OVERLAPPING_ONLY \
        -O targets.cohort.gc.filtered.interval_list \
        --verbosity ERROR


# DetermineGermlineContigPloidy in COHORT MODE
gatk DetermineGermlineContigPloidy \
        -L targets.cohort.gc.filtered.interval_list \
        --interval-merging-rule OVERLAPPING_ONLY \
        $COUNTS_LIST \
	--tmp-dir $PWD \
        --contig-ploidy-priors $PLOIDY_XY \
        --output . \
        --output-prefix ploidy \
        --verbosity ERROR


# GermlineCNVCaller in COHORT MODE
gatk GermlineCNVCaller \
        --run-mode COHORT \
        -L targets.cohort.gc.filtered.interval_list \
        $COUNTS_LIST \
        --contig-ploidy-calls ploidy-calls \
        --annotated-intervals targets.annotated.tsv \
        --interval-merging-rule OVERLAPPING_ONLY \
        --output . \
        --output-prefix cohort \
	--tmp-dir $PWD \
        --verbosity ERROR

cd ..
cd ..

index=0

for sample_id in $FEMALE;
do SAMPLE=${sample_id%%.dedup.bam};

# PostprocessGermlineCNVCalls COHORT MODE
gatk PostprocessGermlineCNVCalls \
        --model-shard-path gatkcnv_output/female/cohort-model \
        --calls-shard-path gatkcnv_output/female/cohort-calls \
        --allosomal-contig chrX --allosomal-contig chrY \
        --contig-ploidy-calls gatkcnv_output/female/ploidy-calls \
        --sample-index $index \
        --output-genotyped-intervals gatkcnv_output/female/genotyped-intervals.$SAMPLE.vcf.gz \
        --output-genotyped-segments gatkcnv_output/female/genotyped-segments.$SAMPLE.vcf.gz \
        --output-denoised-copy-ratios gatkcnv_output/female/denoised-copy-ratios.$SAMPLE.tsv \
        --sequence-dictionary $DIC \
	--tmp-dir $PWD \
        --verbosity ERROR;

let "index+=1";

done

# ######################################################################################
echo ""
echo "Working on male..."
echo ""

MALE=""
while read line
do
MALE+="$line.dedup.bam "
done < male_list.txt

echo "Liste des hommes:"
echo $MALE
echo ""

# Define the resolution of the analysis with a genomic intervals list
#gatk PreprocessIntervals \
#    -R  $REF \
#    -L  $TARGET_XY \
#    -XL $CENTROMETIC_XY \
#    --bin-length 0 \
#    --padding 50 \
#    --interval-merging-rule OVERLAPPING_ONLY \
#    -O gatkcnv_output/male/targets.preprocessed.interval_list \
#    --verbosity ERROR

for sample_id in $MALE;
do SAMPLE=${sample_id%%.dedup.bam}; \

# Collect raw integer counts data
gatk CollectReadCounts \
        -L $XY_Preprocessed_Interval \
        -XL $CENTROMETIC_XY \
        -R $REF \
	--tmp-dir $PWD \
        --interval-merging-rule OVERLAPPING_ONLY \
        -I $SAMPLE.dedup.bam \
        --format TSV \
        -O gatkcnv_output/male/$SAMPLE.tsv \
        --verbosity ERROR;

done

# AnnotateIntervals with GC content
gatk AnnotateIntervals \
    -L $XY_Preprocessed_Interval \
    -XL $CENTROMETIC_XY \
    -R $REF \
    --tmp-dir $PWD \
    -imr OVERLAPPING_ONLY \
    --mappability-track $MAPPABILITY \
    --segmental-duplication-track $SEGMENTAL_DUPLICATE \
    -O gatkcnv_output/male/targets.annotated.tsv \
    --verbosity ERROR


cd gatkcnv_output/male

COUNTS_LIST=""
for counts_file in *.tsv; do
        if [ $(grep -v targets.annotated.tsv <<<$counts_file) ]; then
        FILE=${counts_file%%.tsv};
        COUNTS_LIST+="-I $FILE.tsv ";
        fi
done

# FilterIntervals based on GC-content and cohort extreme counts
gatk FilterIntervals \
        -L $XY_Preprocessed_Interval \
        -XL $CENTROMETIC_XY \
	--tmp-dir $PWD \
        --annotated-intervals targets.annotated.tsv \
        $COUNTS_LIST \
        -imr OVERLAPPING_ONLY \
        -O targets.cohort.gc.filtered.interval_list \
        --verbosity ERROR


# DetermineGermlineContigPloidy in COHORT MODE
gatk DetermineGermlineContigPloidy \
        -L targets.cohort.gc.filtered.interval_list \
        --interval-merging-rule OVERLAPPING_ONLY \
        $COUNTS_LIST \
	--tmp-dir $PWD \
        --contig-ploidy-priors $PLOIDY_XY \
        --output . \
        --output-prefix ploidy \
        --verbosity ERROR


# GermlineCNVCaller in COHORT MODE
gatk GermlineCNVCaller \
        --run-mode COHORT \
        -L targets.cohort.gc.filtered.interval_list\
        $COUNTS_LIST \
        --contig-ploidy-calls ploidy-calls \
        --annotated-intervals targets.annotated.tsv \
        --interval-merging-rule OVERLAPPING_ONLY \
        --output . \
        --output-prefix cohort \
	--tmp-dir $PWD \
        --verbosity ERROR

cd ..
cd ..

index=0

for sample_id in $MALE;
do SAMPLE=${sample_id%%.dedup.bam};

# PostprocessGermlineCNVCalls COHORT MODE
gatk PostprocessGermlineCNVCalls \
        --model-shard-path gatkcnv_output/male/cohort-model \
        --calls-shard-path gatkcnv_output/male/cohort-calls \
        --allosomal-contig chrX --allosomal-contig chrY \
        --contig-ploidy-calls gatkcnv_output/male/ploidy-calls \
        --sample-index $index \
        --output-genotyped-intervals gatkcnv_output/male/genotyped-intervals.$SAMPLE.vcf.gz \
        --output-genotyped-segments gatkcnv_output/male/genotyped-segments.$SAMPLE.vcf.gz \
        --output-denoised-copy-ratios gatkcnv_output/male/denoised-copy-ratios.$SAMPLE.tsv \
        --sequence-dictionary $DIC \
	--tmp-dir $PWD \
        --verbosity ERROR;

let "index+=1";

done

# ########################################################################################
echo ""
echo "Working on all..."
echo ""

Define the resolution of the analysis with a genomic intervals list
#gatk PreprocessIntervals \
#   -R $REF \
#   -L $TARGET_AUTO \
#   -XL $CENTROMETIC_AUTO \
#   --bin-length 0 \
#   --padding 50 \
#   --interval-merging-rule OVERLAPPING_ONLY \
#   -O gatkcnv_output/all/targets.preprocessed.interval_list \
#   --verbosity ERROR

for sample_id in *.dedup.bam;
do SAMPLE=${sample_id%%.dedup.bam}; \

#Collect raw integer counts data
gatk CollectReadCounts \
   -L $AUTOSOME_Preprocessed_Interval \
   -XL $CENTROMETIC_AUTO \
   -R $REF \
   --tmp-dir $PWD \
   --interval-merging-rule OVERLAPPING_ONLY \
   -I $SAMPLE.dedup.bam \
   --format TSV \
   -O gatkcnv_output/all/$SAMPLE.tsv \
   --verbosity ERROR;

done

#AnnotateIntervals with GC content
gatk AnnotateIntervals \
   -L $AUTOSOME_Preprocessed_Interval \
   -XL $CENTROMETIC_AUTO \
   -R $REF \
   --tmp-dir $PWD \
   --mappability-track $MAPPABILITY \
   --segmental-duplication-track $SEGMENTAL_DUPLICATE \
   -imr OVERLAPPING_ONLY \
   -O gatkcnv_output/all/targets.annotated.tsv \
   --verbosity ERROR


cd gatkcnv_output/all

COUNTS_LIST=""
for counts_file in *.tsv; do
        if [ $(grep -v targets.annotated.tsv <<<$counts_file) ]; then
        FILE=${counts_file%%.tsv};
        COUNTS_LIST+="-I $FILE.tsv ";
        fi
done

#FilterIntervals based on GC-content and cohort extreme counts
gatk FilterIntervals \
       -L $AUTOSOME_Preprocessed_Interval \
       -XL $CENTROMETIC_AUTO \
       --tmp-dir $PWD \
       --annotated-intervals targets.annotated.tsv \
       $COUNTS_LIST \
       -imr OVERLAPPING_ONLY \
       -O targets.cohort.gc.filtered.interval_list \
       --verbosity ERROR


# DetermineGermlineContigPloidy in COHORT MODE
gatk DetermineGermlineContigPloidy \
        -L targets.cohort.gc.filtered.interval_list \
        --interval-merging-rule OVERLAPPING_ONLY \
        $COUNTS_LIST \
	--tmp-dir $PWD \
        --contig-ploidy-priors $PLOIDY_AUTO \
        --output . \
        --output-prefix ploidy \
        --verbosity ERROR

#making directory for scatter interval
mkdir -p scatter

#scatter interval
gatk --java-options "-Xmx8G" IntervalListTools \
--INPUT targets.cohort.gc.filtered.interval_list \
--SUBDIVISION_MODE INTERVAL_COUNT \
--SCATTER_CONTENT 5000 \
--OUTPUT scatter
 	
# GermlineCNVCaller in COHORT MODE
gatk GermlineCNVCaller \
        --run-mode COHORT \
        -L targets.cohort.gc.filtered.interval_list \
	$COUNTS_LIST \
        --contig-ploidy-calls ploidy-calls \
        --annotated-intervals targets.annotated.tsv \
        --interval-merging-rule OVERLAPPING_ONLY \
        --output . \
        --output-prefix cohort \
	--tmp-dir $PWD \
        --verbosity ERROR

cd ..
cd ..

index=0

for sample_id in *.dedup.bam;
do SAMPLE=${sample_id%%.dedup.bam};

# PostprocessGermlineCNVCalls COHORT MODE
gatk PostprocessGermlineCNVCalls \
        --model-shard-path gatkcnv_output/all/cohort-model \
	--calls-shard-path gatkcnv_output/all/cohort-calls \
        --allosomal-contig chrX --allosomal-contig chrY \
        --contig-ploidy-calls gatkcnv_output/all/ploidy-calls \
        --sample-index $index \
        --output-genotyped-intervals gatkcnv_output/all/genotyped-intervals.$SAMPLE.vcf.gz \
        --output-genotyped-segments gatkcnv_output/all/genotyped-segments.$SAMPLE.vcf.gz \
	--output-denoised-copy-ratios gatkcnv_output/all/denoised-copy-ratios.$SAMPLE.tsv \
        --sequence-dictionary $DIC \
	--tmp-dir $PWD \
        --verbosity ERROR;

let "index+=1";

done

echo ""
echo "GATK4 CNV DETECTION job done!"
echo ""
