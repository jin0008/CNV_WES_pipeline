# Author: Julie BOGOIN

source ~/miniconda3/etc/profile.d/conda.sh
conda activate gatk_env

REF="/media/Data1/jbogoin/ref/fa_hg38/hg38_GenDev/hg38_GenDev.fa"
DIC="/media/Data1/jbogoin/ref/fa_hg38/hg38_GenDev/hg38_GenDev.dict"

CENTROMETIC_AUTO="/media/Data1/jbogoin/ref/fa_hg38/hg38_centromeric/centromeric_regions_autosomes.bed"
CENTROMETIC_XY="/media/Data1/jbogoin/ref/fa_hg38/hg38_centromeric/centromeric_regions_XY.bed"

TARGET_AUTO="/media/Data1/jbogoin/ref/gencode/v34_hg38/autosomes/gencode.v34.basic.annotation.autosome.interval_list"
TARGET_XY="/media/Data1/jbogoin/ref/gencode/v34_hg38/XY/gencode.v34.basic.annotation.XY.interval_list"

PLOIDY_AUTO="/media/Data1/jbogoin/ref/contig_ploidy_priors/ploidy_priors_table_autosome.tsv"
PLOIDY_XY="/media/Data1/jbogoin/ref/contig_ploidy_priors/ploidy_priors_table_XY.tsv"

echo ""
echo "GATK4 CNV DETECTION start"
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

echo "Liste des femmes:"
echo $FEMALE
echo ""

# Define the resolution of the analysis with a genomic intervals list
gatk PreprocessIntervals \
    -R $REF \
    -L $TARGET_XY \
    -XL $CENTROMETIC_XY \
    --bin-length 0 \
    --padding 50 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O gatkcnv_output/female/targets.preprocessed.interval_list \
    --verbosity ERROR

for sample_id in $FEMALE;
do SAMPLE=${sample_id%%.dedup.bam}; \

# Collect raw integer counts data
gatk CollectReadCounts \
        -L gatkcnv_output/female/targets.preprocessed.interval_list \
        -XL $CENTROMETIC_XY \
        -R $REF \
        --interval-merging-rule OVERLAPPING_ONLY \
        -I $SAMPLE.dedup.bam \
        --format TSV \
        -O gatkcnv_output/female/$SAMPLE.tsv \
        --verbosity ERROR;

done

# AnnotateIntervals with GC content
gatk AnnotateIntervals \
    -L gatkcnv_output/female/targets.preprocessed.interval_list \
    -XL $CENTROMETIC_XY \
    -R $REF \
    -imr OVERLAPPING_ONLY \
    -O gatkcnv_output/female/targets.annotated.tsv \
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
        -L targets.preprocessed.interval_list \
        -XL $CENTROMETIC_XY \
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
        --verbosity ERROR;

let "index+=1";

done

#######################################################################################
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
gatk PreprocessIntervals \
    -R  $REF \
    -L  $TARGET_XY \
    -XL $CENTROMETIC_XY \
    --bin-length 0 \
    --padding 50 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O gatkcnv_output/male/targets.preprocessed.interval_list \
    --verbosity ERROR

for sample_id in $MALE;
do SAMPLE=${sample_id%%.dedup.bam}; \

# Collect raw integer counts data
gatk CollectReadCounts \
        -L gatkcnv_output/male/targets.preprocessed.interval_list \
        -XL $CENTROMETIC_XY \
        -R $REF \
        --interval-merging-rule OVERLAPPING_ONLY \
        -I $SAMPLE.dedup.bam \
        --format TSV \
        -O gatkcnv_output/male/$SAMPLE.tsv \
        --verbosity ERROR;

done

# AnnotateIntervals with GC content
gatk AnnotateIntervals \
    -L gatkcnv_output/male/targets.preprocessed.interval_list \
    -XL $CENTROMETIC_XY \
    -R $REF \
    -imr OVERLAPPING_ONLY \
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
        -L targets.preprocessed.interval_list \
        -XL $CENTROMETIC_XY \
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
        --verbosity ERROR;

let "index+=1";

done

########################################################################################
echo ""
echo "Working on all..."
echo ""

Define the resolution of the analysis with a genomic intervals list
gatk PreprocessIntervals \
    -R  $REF \
    -L  $TARGET_AUTO \
    -XL $CENTROMETIC_AUTO \
    --bin-length 0 \
    --padding 50 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O gatkcnv_output/all/targets.preprocessed.interval_list \
    --verbosity ERROR

for sample_id in *.dedup.bam;
do SAMPLE=${sample_id%%.dedup.bam}; \

# Collect raw integer counts data
gatk CollectReadCounts \
-L gatkcnv_output/all/targets.preprocessed.interval_list \
    -XL $CENTROMETIC_AUTO \
    -R $REF \
    --interval-merging-rule OVERLAPPING_ONLY \
    -I $SAMPLE.dedup.bam \
    --format TSV \
    -O gatkcnv_output/all/$SAMPLE.tsv \
    --verbosity ERROR;

done

# AnnotateIntervals with GC content
gatk AnnotateIntervals \
    -L gatkcnv_output/all/targets.preprocessed.interval_list \
    -XL $CENTROMETIC_AUTO \
    -R $REF \
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

# FilterIntervals based on GC-content and cohort extreme counts
gatk FilterIntervals \
        -L targets.preprocessed.interval_list \
        -XL $CENTROMETIC_AUTO \
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
        --contig-ploidy-priors $PLOIDY_AUTO \
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
        --verbosity ERROR;

let "index+=1";

done

echo ""
echo "GATK4 CNV DETECTION job done!"
echo ""
