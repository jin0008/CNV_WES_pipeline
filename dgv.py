# -*- coding: utf-8 -*-
# Author: Julie Bogoin

import pandas
import os

print("\nDVG program openning.\n")

clinvar = pandas.read_csv('clinavar_results.csv',index_col=None, header=[0])

dgv = pandas.read_csv('/media/jbogoin/Data1/jbogoin/ref/dgv/DGV_hg38_variants_2020-02-25_light.bed.gz', \
    index_col=None, header=[0], compression='gzip', sep='\t')

# Supprimer le 'chr' de la colonne contig
dgv['#chr'] = dgv['#chr'].str[3:]

dgv.sort_values(by=['#chr', 'start'])

count_r = 0
pct = 0
dgv_details_list = []

for index_c, row_c in clinvar.iterrows():

    details_list = []

    for index_g, row_g in dgv.iterrows():

        if row_c['contig'] == row_g['#chr']: 
        
            # identité cnv - dgv OU cnv inclus dans dgv
            if (row_c['start'] >= row_g['start']) and (row_c['end'] <= row_g['end']):
                pct = 100.0

            # dgv inclus dans cnv
            elif (row_c['start'] < row_g['start']) and (row_c['end'] > row_g['end']):
                pct = (row_g['samplesize'] / row_c['size']) * 100.0

            # overlap
            elif (row_c['start'] < row_g['end']) and (row_c['end'] > row_g['start']):
                max_start = max(row_c['start'], row_g['start'])
                min_end = min(row_c['end'], row_g['end'])
                overlap_size = int(min_end) - int(max_start)
                pct = (overlap_size / row_c['size']) * 100.0

            if (pct >= 80.0 and row_g['observedgains'] != '' and row_g['observedlosses'] != ''):
                
                if ((row_c['effect'] == 'deletion' and row_g['observedlosses'] > 0) \
                or (row_c['effect'] == 'duplication' and row_g['observedgains'] > 0)):
                    
                    count_r =  count_r + 1
                    details_list.append(max_start) 
                    details_list.append(min_end) 
                    details_list.append(row_g['observedgains']) 
                    details_list.append(row_g['observedlosses'])

            details_list.append(count_r)

        dgv_details_list.append(details_list)

clinvar['DGV_details'] = pandas.Series(dgv_details)
clinvar['DGV_details'] = clinvar['DGV_details'].str.join('/')

if os.path.isfile('dgv_results.csv'):
    os.remove('dgv_results.csv')
    print('Previous results file removed.')
 
clinvar.to_csv('dgv_results.csv', index=False)                                                                                          
print("dgv_results.csv generated.\n")

print("DGV program job done!\n")
