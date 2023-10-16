# Author: Julie BOGOIN

source ~/miniconda3/etc/profile.d/conda.sh
conda activate annot_env

echo ""
echo "*******************"
echo "annotation.sh start"
echo "*******************"
echo ""

rm -Rf annovar_output
mkdir -p annovar_output

## Download annotation databases from ANNOVAR or UCSC and save to humandb/ directory
#perl ./annotate_variation.pl -downdb avdblist
#perl annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene humandb

## Retirer les noms de colonnes au fichier interval_results.txt > /annovar_output/ex1.avinput
sed '1d' cnv_with_frequences.txt > ex1.avinput

## Remplacer XF par X et XM par X
sed -i -e "s/XF/X/g" ex1.avinput
sed -i -e "s/XM/X/g" ex1.avinput

perl /home/hanjinu/src/annovar/table_annovar.pl \
ex1.avinput \
/media/hanjinu/PM883/AnnovarDB/humandb/hg38 \
--buildver hg38 \
--out annotation \
--remove \
--otherinfo \
--protocol refGene \
--operation g

mv ex1.avinput annovar_output/
mv annotation.hg38_multianno.txt annovar_output/

#first output file: .variant_function
#contains annotation for all variants, by adding two columns to the beginning of each input line

#second output file: .exonic_variant_function
#contains the amino acid changes as a result of the exonic variant.

echo ""
echo "annotation.sh job done!"
echo ""
