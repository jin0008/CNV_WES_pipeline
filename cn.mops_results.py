# Author: Julie BOGOIN

import os
import pandas

print("\ncn.mops CNV results program openning.\n")

li = []

path = '.'
folders = os.listdir(path)

for folder in folders:
    
    files = os.listdir(folder)

    for name in files:
    
        if "cnvs.csv" in name:

            df = pandas.read_csv('cnvs.csv', sep=',',index_col=None, header=[0])
            df.dropna(how='all')

            sn = df['sampleName'].str.split(pat='.', n=0, expand=True)
            df['sample'] = sn[0]
            li.append(df)

concat  = pandas.concat(li, axis=0, ignore_index=True)

cnv =  concat['CN'].str.split(pat='CN', n=0, expand=True)
concat['CN'] = cnv[1]
concat['CN'] = df['CN'].astype('str').astype('int')

del concat['sampleName']
del concat['Unnamed: 0']
del concat['strand']
del concat['width']
del concat['median']

concat.rename(columns={'seqnames': 'contig'}, inplace=True)
concat.rename(columns={'mean': 'log2copy_ratio'}, inplace=True)

total = concat.shape[0]

concat.query('log2copy_ratio>1 or log2copy_ratio<-1', inplace=True)
concat['cnv_ratio'] = df['log2copy_ratio']**2

concat['effect'] = 'i'

concat.loc[df.log2copy_ratio>1, 'effect'] = "duplication"
concat.loc[df.log2copy_ratio<-1, 'effect'] = "deletion"

df_sex = pandas.read_csv('../samples.txt', header = [0], sep="\t", index_col=None)

frame = pandas.merge(concat, df_sex, left_on='sample', right_on='sample')

cols = ['sample', 'sex', 'contig', 'start', 'end', 'cnv_ratio','log2copy_ratio', 'CN', 'effect']
frame = frame[cols]

print('{0} CNV lines filtred among {1} lines found by cn.mops.'.format(frame.shape[0], total))

if os.path.isfile('cn.mops_results.csv'):
    os.remove('cn.mops_results.csv')
    print('Previous results file removed.')

frame.to_csv('cn.mops_results.csv', index=False)                                                                                           
print("cn.mops_results.csv generated.\n")
print("cn.mops CNV results job done!\n")
