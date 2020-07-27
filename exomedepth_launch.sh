# Author: Julie BOGOIN

source ~/miniconda3/etc/profile.d/conda.sh
conda activate

echo ""
echo "ExomeDepth CNV DETECTION start"
echo ""

rm -Rf exomedepth_output
mkdir exomedepth_output
cd exomedepth_output
mkdir all
mkdir female
mkdir male
cd ..

R -q --vanilla < /home/hanjinu/src/CNV_WES_pipeline/exomedepth_cnv.r

echo ""
echo "ExomeDepth CNV DETECTION job done!"
echo ""
