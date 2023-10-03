# Author: Julien Buratti

source ~/miniconda3/etc/profile.d/conda.sh
conda activate exome_targets_env

echo ""
echo "targets_preparation start"
echo ""

#### Preparing bed files

DIC="/media/hanjinu/PM883/db/refs/hg38_broad/Homo_sapiens_assembly38.dict"

# Download the gencode.v44.basic.annotation.gff3 file from the GENECODE website
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gff3.gz

# Retrieve CDS
zgrep "CDS" gencode.v44.basic.annotation.gff3.gz > gencode.v44.basic.annotation.CDS.gff3

# Add a 4th column and sort
awk -F "\t" '{print $1"\t"$4-1"\t"$5}' gencode.v44.basic.annotation.CDS.gff3 \
    | sort -Vu > gencode.v44.basic.annotation.CDS.bed

# Remove overlapping intervals
bedtools merge -i gencode.v44.basic.annotation.CDS.bed > gencode.v44.basic.annotation.CDS.merged.bed 

# Add column 4
awk -F "\t" '{print $1, $2, $3, "CDS"}' gencode.v44.basic.annotation.CDS.merged.bed > gencode.v44.basic.annotation.CDS.merged.4fields.bed

# Transformer bed en interval.list
gatk BedToIntervalList \
    -I gencode.v44.basic.annotation.CDS.merged.4fields.bed \
    -O gencode.v44.basic.annotation.CDS.merged.4fields.interval_list \
    -SD /media/hanjinu/PM883/db/refs/hg38_broad/Homo_sapiens_assembly38.dict

# Make one target file per chromosome
for i in {1..22} X Y
do
grep "^chr${i}" gencode.v44.basic.annotation.CDS.merged.bed > gencode.v44.basic.annotation.CDS.chr${i}.bed &
done

echo ""
echo "targets_preparation job done!"
echo ""




