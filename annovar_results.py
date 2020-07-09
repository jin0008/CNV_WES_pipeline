# Author: Julie BOGOIN

import os
import pandas

print('\nAnnotation results program openning.\n')

df = pandas.read_csv('annotation.hg38_multianno.txt', sep='\t', index_col=None)

df.columns = ['contig', 'start', 'end', 'ref', 'alt',\
        'function', 'genes', 'gene_details', 'exonic_function', 'aa_change', \
        'sample', 'sex', 'size', 'effect', 'log2copy_ratio', 'CN', 'cnv_tool', \
        'targets_number', 'empty']

if os.path.isfile('annotation_results.csv'):
    os.remove('annotation_results.csv')
    print('Previous results file removed.')

del df['ref']
del df['alt']
del df['aa_change']
del df['empty']

df.to_csv('annotation_results.csv', index=False)
print("annotation_results.csv generated.\n")
print("Annotation CNV results job done!\n")
