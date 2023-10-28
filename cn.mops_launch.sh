# Author: Julie BOGOIN


echo ""
echo "cn.mops CNV DETECTION start"
echo ""

rm -Rf cn.mops_output
mkdir cn.mops_output

R -q --vanilla < /home/hanjinu/src/CNV_WES_pipeline/cn.mops_cnv_female.r
R -q --vanilla < /home/hanjinu/src/CNV_WES_pipeline/cn.mops_cnv_male.r
R -q --vanilla < /home/hanjinu/src/CNV_WES_pipeline/cn.mops_cnv_all.r

echo ""
echo "cn.mops CNV DETECTION job done!"
echo ""
