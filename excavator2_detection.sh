# Author: Julie BOGOIN

source ~/miniconda3/etc/profile.d/conda.sh
conda activate excavator2_env

TARGETS_XY="/media/jbogoin/Data1/jbogoin/ref/gencode/v34_hg38/autosomes/gencode.v34.basic.annotation.auto.merged.bed"
TARGETS_AUTO="/media/Data1/jbogoin/ref/gencode/v34_hg38/XY/gencode.v34.basic.annotation.XY.merged.bed"
TARGETS_ALL="/media/jbogoin/Data1/jbogoin/ref/gencode/v34_hg38/all/gencode.v34.basic.annotation.CDS.merged.bed"

DATA=$PWD

echo ""
echo "EXCAVATOR2 CNV DETECTION start"
echo ""

rm -Rf excavator2_output/
mkdir excavator2_output
cd excavator2_output
mkdir all
mkdir female
mkdir male
cd ..

rm -Rf ~/EXCAVATOR2_Package_v1.1.2/data/targets/hg38/autosomes_w10K/
rm -Rf ~/EXCAVATOR2_Package_v1.1.2/data/targets/hg38/XY_w10K/

# BAM list file
ls *.dedup.bam > $DATA/excavator2_output/all/bam_list.txt

# female list file
FEMALE=""
while read line
do
FEMALE+="$line.dedup.bam\n"
done < female_list.txt

echo -e $FEMALE > $DATA/excavator2_output/female/bam_list.txt

# male list file
MALE=""
while read line
do
MALE+="$line.dedup.bam\n"
done < male_list.txt

echo -e $MALE > $DATA/excavator2_output/male/bam_list.txt


####################################################################################
echo ""
echo "Working on all..."
echo ""

# Target initialization
cd ~/EXCAVATOR2_Package_v1.1.2
perl TargetPerla.pl SourceTarget.txt $TARGETS_ALL all_w10K 10000 hg38

# # Data Prepare Module
# cd $DATA

# for sample_id in *.dedup.bam; \
# do SAMPLE=${sample_id%%.dedup.bam}

# echo $DATA/$SAMPLE.dedup.bam $DATA/excavator2_output/all/$SAMPLE $SAMPLE >> $DATA/excavator2_output/all/ExperimentalFilePrepare.w10K.txt;

# done

# # RC calculations
# cd ~/EXCAVATOR2_Package_v1.1.2
# perl EXCAVATORDataPrepare.pl $DATA/excavator2_output/all/ExperimentalFilePrepare.w10K.txt --processors 12 --target MyTarget_w10K --assembly hg38

# # Experimental analysis file
# labels=("C1" "C2" "C3" "C4" "C5" "C6" "C7" "C8" "C9" "C10" "C11")

# rm $DATA/excavator2_output/all/ExperimentalFileAnalysis.w10K.* 

# cd $DATA

# for sample_id in *.dedup.bam; do 
# 	SAMPLE=${sample_id%%.dedup.bam};	
	
# 	echo 'T1' $DATA/excavator2_output/all/$SAMPLE $SAMPLE > $DATA/excavator2_output/all/ExperimentalFileAnalysis.w10K.$SAMPLE.txt;
	
# 	others=$(ls *.dedup.bam | grep -v $SAMPLE);
	
# 	# TRansformer la liste others en tableau indicable table
# 	table=( ${others// / } )
# 	table_clean=${table[@]/.dedup.bam/}
# 	NORMAL=( ${table_clean// / } )

# 	for i in `seq 0 10`; do
		
# 		echo ${labels[$i]} $DATA/excavator2_output/all/${NORMAL[$i]} ${NORMAL[$i]}\
# 			>> $DATA/excavator2_output/all/ExperimentalFileAnalysis.w10K.$SAMPLE.txt;
	
# 	done

# 	# Segmentation of the WMRC
# 	cd ~/EXCAVATOR2_Package_v1.1.2;
	
# 	perl EXCAVATORDataAnalysis.pl $DATA/excavator2_output/all/ExperimentalFileAnalysis.w10K.$SAMPLE.txt\
# 		--processors 12 --target autosomes_w10K\
# 	    --assembly hg38\
# 		--output $DATA/excavator2_output/all/w10K_results.$SAMPLE\
# 	    --mode pooling;
	
# 	cd $DATA	

# done

# echo ""
# echo "EXCAVATOR2 CNV DETECTION job done!"
# echo ""
