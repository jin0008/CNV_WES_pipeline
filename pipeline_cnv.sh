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
conda activate results_cnv

python /home/hanjinu/src/CNV_WES_pipeline/female_male_lists.py

conda deactivate

# tools cnv detection
bash /home/hanjinu/src/CNV_WES_pipeline/cn.mops_launch.sh
bash /home/hanjinu/src/CNV_WES_pipeline/cnvkit_detection.sh
bash /home/hanjinu/src/CNV_WES_pipeline/excavator2_detection.sh
bash /home/hanjinu/src/CNV_WES_pipeline/exomedepth_launch.sh
bash /home/hanjinu/src/CNV_WES_pipeline/gatk_detection.sh

# tools results generation
conda activate results_cnv

cd cn.mops_output
python /home/hanjinu/src/CNV_WES_pipeline/cn.mops_results.py

cd ../cnvkit_output
python /home/hanjinu/src/CNV_WES_pipeline/cnvkit_results.py

cd ../excavator2_output
python /home/hanjinu/src/CNV_WES_pipeline/excavator2_results.py

cd ../exomedepth_output
python /home/hanjinu/src/CNV_WES_pipeline/exomedepth_results.py

cd ../gatkcnv_output
python /home/hanjinu/src/CNV_WES_pipeline/gatk_results.py

cd $DATA

# results summary
python /home/hanjinu/src/CNV_WES_pipeline/cnv_results.py
python /home/hanjinu/src/CNV_WES_pipeline/cnv_interval_objet_sample.py
python /home/hanjinu/src/CNV_WES_pipeline/cnv_interval_objet_run.py

conda deactivate

# annotations

# annovar
conda activate annot_env

bash /home/hanjinu/src/CNV_WES_pipeline/annovar.sh

conda deactivate

conda activate results_cnv

cd annovar_output
python /home/hanjinu/src/CNV_WES_pipeline/annovar_results.py

# ClinVar
# In_gene
# DGV_count

cd $DATA
python /home/hanjinu/src/CNV_WES_pipeline/combine_annot.py

conda deactivate

echo ""
echo "pipeline_cnv.sh job done!"
echo ""
