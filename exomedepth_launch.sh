# Author: Julie BOGOIN

source ~/miniconda3/etc/profile.d/conda.sh
conda activate exomedepth_env

echo ""
echo "*******************************"
echo "ExomeDepth CNV DETECTION start."
echo "*******************************"
echo ""


# rm -Rf exomedepth_output
# mkdir exomedepth_output
# cd exomedepth_output
# mkdir all
# mkdir female
# mkdir male
# cd ..

R -q --vanilla < ~/SCRIPTS/CNV_WES_pipeline/exomedepth_cnv.r

echo ""
echo "ExomeDepth CNV DETECTION job done!"
echo ""
