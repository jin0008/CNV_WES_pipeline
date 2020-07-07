#!/bin/bash

# Author: Julie BOGOIN

echo ""
echo "pipeline_csv.sh start"
echo ""

source ~/miniconda3/etc/profile.d/conda.sh

# sex determination
conda activate sex_env

python3 ~/CNV_WES_pipeline/sex_determination.py

conda deactivate

# generation des listes hommes et femmes
conda activate results_cnv

python ~/CNV_WES_pipeline/female_male_lists.py

conda deactivate

# tools cnv detection
bash ~/CNV_WES_pipeline/cn.mops_launch.sh
bash ~/CNV_WES_pipeline/cnvkit_detection.sh
bash ~/CNV_WES_pipeline/excavator2_detection.sh
bash ~/CNV_WES_pipeline/exomedepth_launch.sh
bash ~/CNV_WES_pipeline/gatk_detection.sh

# tools results generation
conda activate results_cnv

cd cn.mops_output
python ~/CNV_WES_pipeline/cn.mops_results.py

cd ../cnvkit_output
python ~/CNV_WES_pipeline/cnvkit_results.py

cd ../excavator2_output
python ~/CNV_WES_pipeline/excavator2_results.py

cd ../exomedepth_output
python ~/CNV_WES_pipeline/exomedepth_results.py

cd ../gatk_output
python ~/CNV_WES_pipeline/gatk_results.py

cd ..

# results summary
python ~/CNV_WES_pipeline/cnv_results.py
python ~/CNV_WES_pipeline/cnv_interval_object_sample.py
python ~/CNV_WES_pipeline/cnv_interval_object_run.py

conda deactivate

# annotations

conda activate annot_env

bash ~/CNV_WES_pipeline/annovar.sh

conda deactivate

conda activate results_cnv

cd annotations/annovar_output
python ~/CNV_WES_pipeline/annovar_results.py
cd ..

python ~/CNV_WES_pipeline/clinvar.py
python ~/CNV_WES_pipeline/dvg.py
python ~/CNV_WES_pipeline/in_gene.py

conda deactivate

echo ""
echo "pipeline_csv.sh job done!"
echo ""

