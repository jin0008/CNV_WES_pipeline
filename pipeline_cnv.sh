#!/bin/bash

# Author: Julie BOGOIN
# Modified by: Jinu Han

exec &> cnv_sortie.log

echo ""
echo "pipeline_cnv.sh start"
echo ""

DATA=$PWD

source ~/miniconda3/etc/profile.d/conda.sh

# sex determination
conda activate

python3 /home/hanjinu/src/CNV_WES_pipeline/sex_determination.py

conda deactivate

# generation of men’s and women’s lists
python /home/hanjinu/src/CNV_WES_pipeline/female_male_lists.py


# tools cnv detection
bash /home/hanjinu/src/CNV_WES_pipeline/cn.mops_launch.sh
bash /home/hanjinu/src/CNV_WES_pipeline/exomedepth_launch.sh
bash /home/hanjinu/src/CNV_WES_pipeline/gatk_detection.sh

# tools results generation
conda activate results_cnv

cd cn.mops_output
python /home/hanjinu/src/CNV_WES_pipeline/cn.mops_results.py

cd ../exomedepth_output
python /home/hanjinu/src/CNV_WES_pipeline/exomedepth_results.py

cd ../gatkcnv_output
python /home/hanjinu/src/CNV_WES_pipeline/gatk_results.py

cd $DATA

conda deactivate

echo ""
echo "pipeline_cnv.sh job done!"
echo ""
