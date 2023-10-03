# Author: Julie BOGOIN

echo ""
echo "REFERENCE INDEXATION start"
echo ""

echo ""
echo "Please indicate the name of the .fa.gz archive to process (without extension):"
read reference 

gunzip $reference.fa.gz

bwa index -a bwtsw $reference.fa

gatk CreateSequenceDictionary -R $reference.fa -O $reference.dict

samtools faidx $reference.fa

echo ""
echo "REFERENCE INDEXATION job done!"
echo ""
