# Author: Julien Buratti

source ~/miniconda3/etc/profile.d/conda.sh
conda activate exome_targets_env

echo ""
echo "targets_preparation start"
echo ""

#### Préparation fichiers bed

DIC="/media/hanjinu/SS200/db/refs/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dict"

# Télécharger le fichier gencode.v34.basic.annotation.gff3 sur le site de GENECODE
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.basic.annotation.gff3.gz

# Recupérer les CDS
zgrep "CDS" gencode.v34.basic.annotation.gff3.gz > gencode.v34.basic.annotation.CDS.gff3

# Ajouter une 4ème colonne et trier
awk -F "\t" '{print $1"\t"$4-1"\t"$5}' gencode.v34.basic.annotation.CDS.gff3 \
    | sort -Vu > gencode.v34.basic.annotation.CDS.bed

# Supprimer les intervalles chevauchants
bedtools merge -i gencode.v34.basic.annotation.CDS.bed > gencode.v34.basic.annotation.CDS.merged.bed 

# Ajouter une colonne 4
awk -F "\t" '{print $1, $2, $3, "CDS"}' gencode.v34.basic.annotation.CDS.merged.bed > gencode.v34.basic.annotation.CDS.merged.4fields.bed

# Transformer bed en interval.list
gatk BedToIntervalList \
    -I gencode.v34.basic.annotation.CDS.merged.4fields.bed \
    -O gencode.v34.basic.annotation.CDS.merged.4fields.interval_list \
    -SD /media/hanjinu/SS200/db/refs/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dict

# Faire un fichier cible par chromosome
for i in {1..22} X Y; do ;
grep "^chr${i}" gencode.v34.basic.annotation.CDS.merged.bed > gencode.v34.basic.annotation.CDS.chr${i}.bed; 
done

echo ""
echo "targets_preparation job done!"
echo ""




