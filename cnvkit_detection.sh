# Author: Julie BOGOIN

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cnvkit_env

REF="/media/Data1/jbogoin/ref/fa_hg38/hg38_GenDev/hg38_GenDev.fa"

TARGET_AUTO="/media/Data1/jbogoin/ref/gencode/v34_hg38/autosomes/gencode.v34.basic.annotation.autosome.interval_list"
TARGET_XY="/media/Data1/jbogoin/ref/gencode/v34_hg38/XY/gencode.v34.basic.annotation.XY.interval_list"

DATA=$PWD

echo ""
echo "CNVKIT DETECTION start"
echo ""

rm -Rf cnvkit_output
mkdir cnvkit_output
cd cnvkit_output
mkdir all
mkdir female
mkdir male
cd ..

####################################################################################
echo ""
echo "Working on female..."
echo ""

FEMALE=""
while read line
do
FEMALE+="$line.dedup.bam "
done < female_list.txt

echo "Liste des femmes:"
echo $FEMALE
echo ""

#create a pooled reference by running batch command specifying only normal samples
~/cnvkit/cnvkit.py batch -n $FEMALE \
	-m hybrid \
    -f $REF \
	--targets $TARGET_XY \
	--output-reference $DATA/cnvkit_output/female/pooled-reference.cnn \
	--short-names \
	-d $DATA/cnvkit_output/female \
	-p 12

# Run WGS batch pipeline using the pooled reference file
~/cnvkit/cnvkit.py batch $FEMALE \
	-m hybrid \
	-r $DATA/cnvkit_output/female/pooled-reference.cnn \
	-d $DATA/cnvkit_output/female \
	-p 12 \
	--scatter \
	--diagram

rm $DATA/cnvkit_output/female/*.target.bed;


####################################################################################
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

#create a pooled reference by running batch command specifying only normal samples
~/cnvkit/cnvkit.py batch -n $MALE \
	-m hybrid \
    -f $REF \
	--targets $TARGET_XY \
	--output-reference $DATA/cnvkit_output/male/pooled-reference.cnn \
	--short-names \
	-d $DATA/cnvkit_output/male \
	-p 12

# Run WGS batch pipeline using the pooled reference file
~/cnvkit/cnvkit.py batch $MALE \
	-m hybrid \
	-r $DATA/cnvkit_output/male/pooled-reference.cnn \
	-d $DATA/cnvkit_output/male \
	-p 12 \
	--scatter \
	--diagram

rm $DATA/cnvkit_output/male/*.target.bed;


####################################################################################
echo ""
echo "Working on all..."
echo ""

SAMPLES=$(ls *.dedup.bam)

#create a pooled reference by running batch command specifying only normal samples
~/cnvkit/cnvkit.py batch -n $SAMPLES \
	-m hybrid \
    -f $REF \
	--targets $TARGET_AUTO \
	--output-reference $DATA/cnvkit_output/all/pooled-reference.cnn \
	--short-names \
	-d $DATA/cnvkit_output/all \
	-p 12

# Run WGS batch pipeline using the pooled reference file
~/cnvkit/cnvkit.py batch $SAMPLES \
	-m hybrid \
	-r $DATA/cnvkit_output/all/pooled-reference.cnn \
	-d $DATA/cnvkit_output/all \
	-p 12 \
	--scatter \
	--diagram

rm $DATA/cnvkit_output/all/*.target.bed;

echo ""
echo "CNVKIT DETECTION job done!"
echo ""
