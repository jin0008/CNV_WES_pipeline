# Author: Julie BOGOIN

import os
import pandas

print("\n************************************")
print("cn.mops CNV results program openning.")
print("************************************\n")

if os.path.isfile('cn.mops_results.csv'):
    os.remove('cn.mops_results.csv')
    print('Previous results file removed.')

li = []

path = '.'
files = os.listdir(path)

for name in files:

    if "cnvs" in name:
        
        taille = os.path.getsize(name)

        if taille != 0:

            df = pandas.read_csv(name, sep=',', index_col=None, header=[0])
            df.dropna(how='all')
            sn = df['sampleName'].str.split(pat='.', n=0, expand=True)
            df['sample'] = sn[0]
            li.append(df)

concat  = pandas.concat(li, axis=0, ignore_index=True)

cnv =  concat['CN'].str.split(pat='CN', n=0, expand=True)
concat['CN'] = cnv[1]

del concat['sampleName']
del concat['Unnamed: 0']
del concat['strand']
del concat['width']
del concat['median']

concat.rename(columns={'seqnames': 'contig'}, inplace=True)
concat.rename(columns={'mean': 'log2copy_ratio'}, inplace=True)

total = concat.shape[0]

concat['log2copy_ratio'] = concat['log2copy_ratio'].astype('str').astype('float')

concat.query('log2copy_ratio>0.4 or log2copy_ratio<-0.7', inplace=True)

concat['cnv_ratio'] = concat['log2copy_ratio']**2

concat['effect'] = 'i'
concat.loc[concat.log2copy_ratio>0.4, 'effect'] = "duplication"
concat.loc[concat.log2copy_ratio<-0.7, 'effect'] = "deletion"

df_sex = pandas.read_csv('../samples.txt', header = [0], sep="\t", index_col=None)

frame = pandas.merge(concat, df_sex, left_on='sample', right_on='sample')

cols = ['sample', 'sex', 'contig', 'start', 'end', 'cnv_ratio','log2copy_ratio', 'CN', 'effect']
frame = frame[cols]

#trier par samples puis contigs
frame.sort_values(by=['sample','contig'])

print('{0} CNV lines filtred among {1} lines found by cn.mops.'.format(frame.shape[0], total))

frame.to_csv('cn.mops_results.csv', index=False)                                                                                           
print("cn.mops_results.csv generated.\n")
print("cn.mops CNV results job done!\n")
