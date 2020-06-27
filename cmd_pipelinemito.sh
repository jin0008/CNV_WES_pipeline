#!/bin/bash

# Author: Julie BOGOIN

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pipmitoDijon_env

DATA=$pwd
REF='/media/jbogoin/Data1/jbogoin/ref/hg38_Mlast/'

echo ""
echo "PIPELINE MITO start"
echo ""

rm -Rf pipelinemito_output/
mkdir pipelinemito_output/

# noms des fastq en minuscle obligatoire!
echo "noms des fastq: conversion des majuscules en minuscules:"
echo ""
for file in *.fastq.gz;
do
    new=`echo $file |tr '[:upper:]' '[:lower:]'`;
    echo "transformation $file => $new";
    mv -i "$file" "$new"
done

echo ""
echo "Liste des fastq à traiter:"
echo ""
if ! ls *.fastq.gz;
then
    echo "Il n'y a pas de fastq dans ce répertoire!";
    echo "Utiliser bam_to_fastq.sh pour les générer.";
    echo ""
    echo "PIPELINE MITO job not done!"
    echo ""
 
else
    cd ~/pipelinemito;
    echo ""
    sudo docker load -i pipelinemitov1.tar;
    echo ""
    sudo docker run -v $DATA:/data:rw\
        -v $REF:/mitopipeline:ro\
        --env THREAD=8\
        --env REFNAME=hg38_GenDev.fa\
        pipelinemitov1;
    
    echo "Déplacement des fichiers de résultats vers /pipelinemito_output"
    echo ""
    sudo mv *.tsv pipelinemito_output;
    sudo mv *.vcf pipelinemito_output;
    sudo mv *.failed pipelinemito_output;
    sudo mv *.log pipelinemito_output;

    echo ""
    echo "PIPELINE MITO job done!"
    echo ""
fi
