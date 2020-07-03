# -*- coding: utf-8 -*-
# Author: Julie BOGOIN

import pandas
import os

print("\nClinVar program openning.\n")

annot = pandas.read_csv('annovar_output/annotation_results.csv',index_col=None, header=[0])
annot['genes'] = annot['genes'].astype('str')
annot['genes'] = annot['genes'].str.split(pat=';', expand=False)

#os.system('cd /media/Data1/jbogoin/ref/clinvar')
#os.system('wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id')

clinvar = pandas.read_csv('/media/Data1/jbogoin/ref/clinvar/gene_condition_source_id', sep= '\t', index_col=None, header=[0])

del clinvar['ConceptID']
del clinvar['SourceID']
del clinvar['DiseaseMIM']
del clinvar['LastUpdated']
del clinvar['#GeneID']

clinvar['RelatedGenes'] = clinvar['RelatedGenes'].astype('str')
clinvar['RelatedGenes'].replace(['nan'], '', inplace=True)

clinvar['AssociatedGenes'] = clinvar['AssociatedGenes'].astype('str')
clinvar['AssociatedGenes'].replace(['nan'], '', inplace=True)

clinvar['clinvar_gene'] = clinvar['RelatedGenes'] + clinvar['AssociatedGenes']
clinvar['clinvar_gene'] = clinvar['clinvar_gene'].astype('str')

del clinvar['RelatedGenes']
del clinvar['AssociatedGenes']

# Recherche des genes annotes dans ClinVar
match_all = []
for gene_list in annot['genes']:
    match_line = []
    for gene in gene_list:
        for index, row in clinvar.iterrows():
            if gene == row['clinvar_gene']:
                match_line.append(row['DiseaseName'])
                continue
        match_line.append(".")
    match_all.append(match_line)

annot['ClinVar_match'] = pandas.Series(match_all)

annot['ClinVar_match'] = annot['ClinVar_match'].str.join('/')
annot['genes'] = annot['genes'].str.join('/')

if os.path.isfile('clinavar_results.csv'):
    os.remove('clinavar_results.csv')
    print('Previous results file removed.')
 
annot.to_csv('clinavar_results.csv', index=False)                                                                                          
print("clinavar_results.csv generated.\n")

print("ClinVar program job done!\n")
