# -*- coding: utf-8 -*-
# Author: Julie Bogoin

import pandas
import os

print("\nDVG program openning.\n")

# dgv = pandas.read_csv('dgv_results.csv',index_col=None, header=[0])
# dgv.sort_values(by=['contig', 'start'])

uniq = pandas.read_csv("/media/jbogoin/Data1/jbogoin/ref/fa_hg38/genes_uniq/hg38_genes_uniq.bed.gz", \
    index_col=None, header=None, compression='gzip', sep='\t')

uniq.columns = ['contig', 'start', 'stop','gene']

# Supprimer le 'chr' de la colonne contig
uniq['contig'] = uniq['contig'].str[3:]

inGene = 0

for index_d, row_d in dgv.iterrows():

    nb_genes = len(row_d['genes'].str.split("/"))

    if nb_genes == 1:

        for index_u, row_u in uniq.iterrows():

            if row_d['contig'] == row_d['contig'] \
                row_d['start'] > row_u['start'] and row_d['end'] < row_u['end']:
                inGene = 1
                break




