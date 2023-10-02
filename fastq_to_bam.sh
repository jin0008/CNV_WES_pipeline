#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate fastq_bam_env

# hg38
REF="/media/hanjinu/PM883/db/refs/hg38_broad/Homo_sapiens_assembly38.fasta"

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
    
    MAPPING_CMD=$(echo "dragen-os -r /media/hanjinu/PM883/db/refs/hg38_broad/ \
    --num-threads $1 --RGID $RG --RGSM ${SAMPLE} -1 $R1 -2 $R2 | \
    samtools view -@ $1 -Sb | samtools sort -n -@ $1 -o ${SAMPLE}.mapped.bam -");
    eval $MAPPING_CMD;

 done;

### MARK DUPLICATES
for i in *.bam; 
    
    do 
    SAMPLE=${i%.*};
    
    gatk --java-options "${java_opt}" MarkDuplicatesSpark \
    -I ${SAMPLE}.mapped.bam \
    -O ${SAMPLE}.dedup.bam \
    -M ${SAMPLE}.mark_dup_metrics.txt \
    --spark-master local[$1] \
    --optical-duplicate-pixel-distance 2500 \
    --tmp-dir $PWD \
    --verbosity ERROR;

    rm -rf ${SAMPLE}.mapped.bam
    
done

