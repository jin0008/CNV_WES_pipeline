# Author: Julie BOGOIN

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cn.mops_env

echo ""
echo "cn.mops CNV DETECTION start"
echo ""

rm -Rf cn.mops_output
mkdir cn.mops_output

R -q --vanilla < ~/CNV_WES_pipeline/cn.mops_cnv_female.r
R -q --vanilla < ~/CNV_WES_pipeline/cn.mops_cnv_male.r
R -q --vanilla < ~/CNV_WES_pipeline/cn.mops_cnv_all.r

echo ""
echo "cn.mops CNV DETECTION job done!"
echo ""