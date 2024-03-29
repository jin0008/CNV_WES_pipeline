# -*- coding: utf-8 -*-
#!/usr/bin/env python3.9

import pandas
import os
from subprocess import Popen, PIPE

def tabix_query(filename, chrom, start, end):
	"""Call tabix and generate an array of strings for each line it returns."""
	query = '{}:{}-{}'.format(chrom, start, end)
	process = Popen(['tabix', '-f', filename, query], stdout=PIPE)
	for line in process.stdout:
		yield line.strip().split()


print("\nAnnots program openning.\n")

DGV = "/media/hanjinu/PM883/db/refs/dgv/DGV_hg38_variant.bed.gz"

###########
# CLINVAR #
###########

#os.system('cd /media/hanjinu/PM883/db/refs/clinvar')
#os.system('wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id')

clinvar = pandas.read_csv('/media/hanjinu/PM883/db/refs/clinvar/gene_condition_source_id', sep= '\t', index_col=None, header=[0])

del clinvar['ConceptID']
del clinvar['SourceID']
del clinvar['DiseaseMIM']
del clinvar['LastUpdated']
del clinvar['#GeneID']
del clinvar['SourceName']

clinvar['RelatedGenes'] = clinvar['RelatedGenes'].astype('str')
clinvar['RelatedGenes'].replace(['nan'], '', inplace=True)

clinvar['AssociatedGenes'] = clinvar['AssociatedGenes'].astype('str')
clinvar['AssociatedGenes'].replace(['nan'], '', inplace=True)

clinvar['clinvar_gene'] = clinvar['RelatedGenes'] + clinvar['AssociatedGenes']
clinvar['clinvar_gene'] = clinvar['clinvar_gene'].astype('str')

clinvar['DiseaseName'] = clinvar['DiseaseName'].astype('str')

del clinvar['RelatedGenes']
del clinvar['AssociatedGenes']

###########
# IN GENE #
###########

uniq = pandas.read_csv('/media/hanjinu/PM883/db/refs/genes_uniq/hg38_genes_uniq.bed', index_col=None, header=None, sep='\t')

uniq.columns = ['contig', 'start', 'stop','gene']

# Remove the 'chr' from the contig column
uniq['contig'] = uniq['contig'].str[3:]

uniq.rename(columns={'gene': 'hg38_gene'}, inplace=True)


###########
# Annovar #
###########

annovar = pandas.read_csv('annovar_output/annotation_results.csv',index_col=None, header=[0])
annovar['genes'] = annovar['genes'].astype('str')
annovar['genes'] = annovar['genes'].str.split(pat=';', expand=False)
annovar.rename(columns={'genes': 'annovar_genes'}, inplace=True)


# Browse the df annovar
clinvar_all = []
ingene_list = []
dgv_details_s = []
count_all = []

for index, row in annovar.iterrows():

#ClinVar
###############################################################
	clinvar_line = []

	for gene in row['annovar_genes']:

		dfc = clinvar[clinvar['clinvar_gene'] == gene]

		if len(dfc.index) == 0:
			clinvar_line.append('.')

		if len(dfc.index) >= 1:
			disease_list = []

			for index_dfc, row_dfc in dfc.iterrows():
				disease_list.append(row_dfc['DiseaseName'])

			clinvar_line.append(" + ".join(disease_list))

#Ingene
###############################################################
	if len(row['annovar_genes']) == 1 :

		dfg = uniq[ (uniq['contig'] == row['contig']) & \
		(uniq['start'] > row['start']) & (uniq['stop'] > int(row['end'])) ]

		if len(dfg.index) == 0:
			ingene_list.append(0)

		if len(dfg.index) >= 1:
			ingene_list.append(1)

	else:
		ingene_list.append('-')

#DGV
##############################################################
	count_r = 0

	dgv_details_l = []

	records = tabix_query(DGV,('chr' + row['contig']), row['start'], row['end'])

	for record in records:

		if len(record) == 6:
			pct = 0.0
			chrom_r = record[0].decode("utf-8")
			start_r = int(record[1].decode("utf-8"))
			end_r = int(record[2].decode("utf-8"))
			obsG_r = int(record[4].decode("utf-8"))
			obsL_r = int(record[5].decode("utf-8"))
			size_r = int(end_r) - int(start_r)

			# identity cnv - dgv OR cnv included in dgv
			if (row['start'] >= start_r) and (row['end'] <= end_r):
				pct = 100.0

				# dgv included in cnv
			elif (row['start'] < start_r) and (row['end'] > end_r):
				pct = (float(size_r) / float(row['size'])) * 100.0

			# overlap
			elif (row['start'] < end_r) and (row['end'] > start_r):
				max_start = max(row['start'], int(start_r))
				min_end = min(row['end'], int(end_r))
				overlap_size = int(min_end) - int(max_start)
				pct = (overlap_size) / row['size'] * 100.0

			if pct >= 80.0 and obsG_r != '' and obsL_r != '':

				if (row['effect'] == "deletion" and int(obsL_r) > 0) or \
					 (row['effect'] == "duplication" and int(obsG_r) > 0):

					count_r += 1

					dgv_details_l.append(chrom_r + ":" + str(start_r) + ":" + str(end_r) + ":" + str(obsG_r) + ":" + str(obsL_r))


	dgv_details_s.append("/".join(dgv_details_l))

	count_all.append(count_r)

	clinvar_all.append("/".join(clinvar_line))


annovar['DiseaseName'] = pandas.Series(clinvar_all)
annovar['in_gene'] = pandas.Series(ingene_list, dtype=str)
annovar['DGV_details'] = pandas.Series(dgv_details_s, dtype=str)
annovar['DGV_count'] = pandas.Series(count_all, dtype=str)

annovar['annovar_genes'] = annovar['annovar_genes'].str.join(' / ')


if os.path.isfile('final_cnv_tab.csv'):
    os.remove('final_cnv_tab.csv')
    print('Previous results file removed.')

annovar.to_csv('final_cnv_tab.csv', sep='\t', index=False)
print("final_cnv_tab.csv generated.\n")

print("Annots program job done!\n")
