# Author: Julie BOGOIN

source ~/miniconda3/etc/profile.d/conda.sh


REF="/media/hanjinu/SS200/db/refs/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"

TARGET_AUTO="/media/hanjinu/SS200/db/refs/gencode/gencode.v34.basic.annotation.autosome.interval_list"
TARGET_XY="/media/hanjinu/SS200/db/refs/gencode/gencode.v34.basic.annotation.XY.scratch.interval_list"

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
FEMALE+="$line.CNV.bam "
done < female_list.txt

echo "Liste des femmes:"
echo $FEMALE
echo ""

#create a pooled reference by running batch command specifying only normal samples
cnvkit.py batch -n $FEMALE \
	-m hybrid \
	-f $REF \
	--targets $TARGET_XY \
	--output-reference $DATA/cnvkit_output/female/pooled-reference.cnn \
	--short-names \
	-d $DATA/cnvkit_output/female \
	-p 12

# Run WGS batch pipeline using the pooled reference file
cnvkit.py batch $FEMALE \
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
MALE+="$line.CNV.bam "
done < male_list.txt

echo "Liste des hommes:"
echo $MALE
echo ""

#create a pooled reference by running batch command specifying only normal samples
cnvkit.py batch -n $MALE \
	-m hybrid \
	-f $REF \
	--targets $TARGET_XY \
	--output-reference $DATA/cnvkit_output/male/pooled-reference.cnn \
	--short-names \
	-d $DATA/cnvkit_output/male \
	-p 12

# Run WGS batch pipeline using the pooled reference file
cnvkit.py batch $MALE \
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

SAMPLES=$(ls *.CNV.bam)

#create a pooled reference by running batch command specifying only normal samples
cnvkit.py batch -n $SAMPLES \
	-m hybrid \
	-f $REF \
	--targets $TARGET_AUTO \
	--output-reference $DATA/cnvkit_output/all/pooled-reference.cnn \
	--short-names \
	-d $DATA/cnvkit_output/all \
	-p 12

# Run WGS batch pipeline using the pooled reference file
cnvkit.py batch $SAMPLES \
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
