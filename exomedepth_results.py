# Author: Julie BOGOIN
# Modified by: Jinu Han


import os
import pandas
import numpy

print("\n**************************************")
print("ExomeDepth CNV results program openning.")
print("**************************************\n")

if os.path.isfile('exomedepth_results.csv'):
    os.remove('exomedepth_results.csv')
    print('Previous results file removed.')

li = []

path = '.'
folders = os.listdir(path)

for folder in folders:
    
    files = os.listdir(folder)

    for name in files:
    
        if ".csv" in name:
            df = pandas.read_csv((folder + '/' + name), sep=',',\
                index_col=None, header=[0])
            df.dropna(how='all')
            sample_name = name.split('.')
            df['sample'] = sample_name[0]
            li.append(df)


concat  = pandas.concat(li, axis=0, ignore_index=True)

df_sex = pandas.read_csv('../samples.txt', header = [0], sep="\t", index_col=None, dtype={'sample':object})

frame = pandas.merge(concat, df_sex, left_on='sample', right_on='sample')

frame['contig'] = 'chr' + frame['chromosome'].astype('str')

del frame['chromosome']
del frame['reads.expected']
del frame['reads.observed']
del frame['end.p']
del frame['start.p']
del frame['id']
del frame['BF']

frame.rename(columns={'reads.ratio':'cnv_ratio'}, inplace=True)

frame.loc[frame.cnv_ratio==0, 'cnv_ratio'] = 0.000001
frame['log2copy_ratio'] = numpy.log2(frame['cnv_ratio'])

frame['CN'] = frame['cnv_ratio'].astype('int')

frame.rename(columns={'type':'effect'}, inplace=True)
frame.rename(columns={'nexons':'targets_number'}, inplace=True)

total = frame.shape[0] - 12

frame.query('log2copy_ratio>0.4 or log2copy_ratio<-0.7', inplace=True)

cols = ['sample', 'sex', 'contig', 'start', 'end', 'cnv_ratio','log2copy_ratio', 'CN', 'effect', 'targets_number']
frame = frame[cols]

print('{0} CNV lines filtred among {1} lines found by ExomeDepth.'.format(frame.shape[0], total))

frame.to_csv('exomedepth_results.csv', index=False)                                                                                           
print("exomedepth_results.csv generated.\n")
print("ExomeDepth CNV results job done!\n")
