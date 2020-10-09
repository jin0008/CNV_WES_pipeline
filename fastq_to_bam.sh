#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate fastq_bam_env

# hg38
REF="/media/hanjinu/SS200/db/refs/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"

# hg19
#REF=/media/hanjinu/SS200/db/refs/b37/human_g1k_v37_decoy.fasta

echo ""
echo "bwa_fastq_to_bam.sh start"
echo ""

#rm *.bam

## MAPPING BWA SEQUENTIAL ###
for R1 in *_R1.fastq.gz; 

    do 
    R2=${R1/_R1/_R2}; 

    SAMPLE=${R1%%_*}; 

    FLOWCELL="$(zcat $R1 | head -1 | awk '{print $1}' | cut -d ":" -f 3)"; 

    DEVICE="$(zcat $R1 | head -1 | awk '{print $1}' | cut -d ":" -f 1 | cut -d "@" -f 2)"; 

    BARCODE="$(zcat $R1 | head -1 | awk '{print $2}' | cut -d ":" -f 4)"; 

    RG=$(echo "\"@RG\tID:${DEVICE}.${FLOWCELL}.${SAMPLE}\tPU:${FLOWCELL}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}-${BARCODE}\"")
    
    MAPPING_CMD=$(echo "bwa mem -M -t 12 -R $RG $REF $R1 $R2 \
    | samtools sort -@ 6 -o ${SAMPLE}.bam -");

    eval $MAPPING_CMD;

 done;

### MARK DUPLICATES
for i in *.bam; 
    
    do 
    SAMPLE=${i%.*};
    
    sambamba markdup -t 12 ${SAMPLE}.bam ${SAMPLE}.dedup.bam;
    
done

