# Author: Julie BOGOIN
# Modified by: Jinu Han

source ~/miniconda2/etc/profile.d/conda.sh
conda activate gatk

REF="/media/hanjinu/SS200/db/refs/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
DIC="/media/hanjinu/SS200/db/refs/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dict"

CENTROMETIC_AUTO="/media/hanjinu/SS200/db/refs/gencode/centromeric_regions_autosomes.bed"
CENTROMETIC_XY="/media/hanjinu/SS200/db/refs/gencode/centromeric_regions_XY.bed"

TARGET_AUTO="/media/hanjinu/SS200/db/refs/gencode/autosomes/gencode.v34.basic.annotation.autosome.interval_list"
TARGET_XY="/media/hanjinu/SS200/db/refs/gencode/gencode.v34.basic.annotation.XY.scratch.interval_list"

PLOIDY_AUTO="/media/hanjinu/SS200/db/refs/contig_ploidy_priors/ploidy_priors_table_autosome.tsv"
PLOIDY_XY="/media/hanjinu/SS200/db/refs/contig_ploidy_priors/ploidy_priors_table_XY.tsv"

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
FEMALE+="$line.CNV.bam ";
done < female_list.txt

echo "Lisf of women:"
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
    -O gatkcnv_output/female/targets.preprocessed.interval_list 

for sample_id in $FEMALE;
do SAMPLE=${sample_id%%.CNV.bam}; \

# Collect raw integer counts data
gatk CollectReadCounts \
        -L gatkcnv_output/female/targets.preprocessed.interval_list \
        -XL $CENTROMETIC_XY \
        -R $REF \
        --interval-merging-rule OVERLAPPING_ONLY \
        -I $SAMPLE.CNV.bam \
        -O gatkcnv_output/female/$SAMPLE.tsv ;

done

# AnnotateIntervals with GC content
gatk AnnotateIntervals \
    -L gatkcnv_output/female/targets.preprocessed.interval_list \
    -XL $CENTROMETIC_XY \
    --mappability-track /media/hanjinu/SS200/db/refs/hg38/k100.umap.bed \
    -R $REF \
    -imr OVERLAPPING_ONLY \
    -O gatkcnv_output/female/targets.annotated.tsv 


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
        --interval-set-rule UNION \
        --annotated-intervals targets.annotated.tsv \
        $COUNTS_LIST \
        -imr OVERLAPPING_ONLY \
        -O targets.cohort.gc.filtered.interval_list 

# DetermineGermlineContigPloidy in COHORT MODE
gatk DetermineGermlineContigPloidy \
        -L targets.cohort.gc.filtered.interval_list \
        --interval-merging-rule OVERLAPPING_ONLY \
        $COUNTS_LIST \
        --contig-ploidy-priors $PLOIDY_XY \
        --output . \
        --output-prefix ploidy \
        --verbosity DEBUG


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
        --verbosity DEBUG

cd ..
cd ..

index=0

for sample_id in $FEMALE;
do SAMPLE=${sample_id%%.CNV.bam};

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
MALE+="$line.CNV.bam "
done < male_list.txt

echo "List of men:"
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
do SAMPLE=${sample_id%%.CNV.bam}; \

# Collect raw integer counts data
gatk CollectReadCounts \
        -L gatkcnv_output/male/targets.preprocessed.interval_list \
        -XL $CENTROMETIC_XY \
        -R $REF \
        --interval-merging-rule OVERLAPPING_ONLY \
        -I $SAMPLE.CNV.bam \
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
do SAMPLE=${sample_id%%.CNV.bam};

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

for sample_id in *.CNV.bam;
do SAMPLE=${sample_id%%.CNV.bam}; \

# Collect raw integer counts data
gatk CollectReadCounts \
-L gatkcnv_output/all/targets.preprocessed.interval_list \
    -XL $CENTROMETIC_AUTO \
    -R $REF \
    --interval-merging-rule OVERLAPPING_ONLY \
    -I $SAMPLE.CNV.bam \
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

for sample_id in *.CNV.bam;
do SAMPLE=${sample_id%%.CNV.bam};

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
