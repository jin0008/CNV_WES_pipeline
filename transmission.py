# Author: Julie BOGOIN

import os
import pandas

#### FUNCTION ####
##################
def family (df, pedigree):
    
    related1 = []
    related2 = []
    lo = []

    for index, row in df.iterrows():

        for index_p, row_p in pedigree.iterrows():

            if row['sample'] == row_p['Index'] :

                related1.append(row_p['Mother'])
                related2.append(row_p['Father'])


    df['related1'] = pandas.Series(related1)
    df['related2'] = pandas.Series(related2)
    df['in_related1'] = False
    df['in_related2'] = False
    df['in_sample'] = True

    for index_p, row_p in pedigree.iterrows():

        index = row_p['Index']
        mother = row_p['Mother']
        father = row_p['Father']

        switch = df[(df['sample'] == index)]
        
        lo.append(switch)
    
    return lo


def overlap (families, df):

    family_list = []
  
    for family in families: 

        family.sort_values(by=['effect','contig', 'start', 'sample'], inplace=True)
        family.reset_index(inplace=True)
        del family['index']

        df_1 = family.shift(periods=1)
        
        family['in_related1'] = (family['effect'] == df_1['effect'])\
                                        & (family['contig'] == df_1['contig'])\
                                        & (family['start'] <= df_1['end']) \
                                        & (family['end'] >= df_1['start']) \
                                        & (family['related1'] == df_1['sample'])

        family['in_related2'] = (family['effect'] == df_1['effect'])\
                                        & (family['contig'] == df_1['contig'])\
                                        & (family['start'] <= df_1['end']) \
                                        & (family['end'] >= df_1['start']) \
                                        & (family['related2'] == df_1['sample'])

        family_list.append(family)

    concat = pandas.concat(family_list, axis=0, ignore_index=True)

    return concat 


def final (concat):
    
    # Creation df final
    final = pandas.DataFrame(columns=['contig','effect', \
    'start', 'end', 'in_sample', 'in_related1', 'in_related2'])
    final['start'] = concat['start']
    final['end'] = concat['end']
    final['contig'] = concat['contig']
    final['sample'] = concat['sample']
    final['effect'] = concat['effect']
    final['in_sample'] = concat['in_sample']
    final['in_related1'] = concat['in_related1']
    final['in_related2'] = concat['in_related2']
    final['score'] = 1
    final.loc[((final.in_related1==True) \
        & (final.in_related2==True)), 'score'] = 3
    final.loc[((final.in_related1==False) \
        & (final.in_related2==True)), 'score'] = 2
    final.loc[((final.in_related1==True) \
        & (final.in_related2==False)), 'score'] = 2

    cols = ['sample', 'contig', 'start', 'end', 'effect', \
    'in_sample', 'in_related1', 'in_related2', 'score']
    final = final[cols]

    return final


#### MAIN ####
##############
print("\nTransmission program openning.\n")

path = "."
dirs = os.listdir(path)

# Creation df final

if os.path.isfile('pedigree.txt'):

    pedigree = pandas.read_csv('pedigree.txt', sep='\t',index_col=None, header=[0])
    del pedigree['Comment']
    del pedigree['Sex']
    del pedigree['Fam']

    hm_path = './homemade_results.csv'
    if os.path.isfile(hm_path):
        df_hm = pandas.read_csv(hm_path, sep='\t', index_col=None)
        df_hm.rename(columns = {'chrom': 'contig', 'cnv_type': 'effect', \
        'log2ratio': 'log2copy_ratio'}, inplace=True)
        df_hm.loc[df_hm.effect=='dup', 'effect'] = "duplication"
        df_hm.loc[df_hm.effect=='del', 'effect'] = "deletion"
        df_hm = df_hm[df_hm.algo == 'WCC']
        del df_hm['algo']
        del df_hm['size']
        del df_hm['function']
        del df_hm['inGene']
        del df_hm['nb_genes']
        del df_hm['genes']
        del df_hm['gene_details']
        del df_hm['freq_run']
        del df_hm['clinvar']
        del df_hm['nb_DGV_hits']
        del df_hm['DGV_details']
        df_hm['start'] = df_hm['start'].astype('str').astype('int')
        df_hm['end'] = df_hm['end'].astype('str').astype('int')
        df_hm['effect'] = df_hm['effect'].astype('str')
        df_hm['contig'] = df_hm['contig'].astype('str')
        df_hm['sample'] = df_hm['sample'].astype('str')
        cols = ['sample', 'contig', 'start', 'end', 'log2copy_ratio', 'effect']
        df_hm = df_hm[cols]        

    for directory in dirs:
    
        if directory == 'gatkcnv_output':
            file_path = './' + directory + '/gatk_results.csv'
            if os.path.isfile(file_path):
                df_gatk = pandas.read_csv(file_path,index_col=None)
    
        if directory == 'cnvkit_output':
            file_path = './' + directory + '/cnvkit_results.csv'
            if os.path.isfile(file_path):
                df_cnvkit = pandas.read_csv(file_path,index_col=None)
        
        if directory == 'exomedepth_output':
            file_path = './' + directory + '/exomedepth_results.csv'
            if os.path.isfile(file_path):
                df_exomedepth = pandas.read_csv(file_path,index_col=None)
        
        if directory == 'cn.mops_output':
            file_path = './' + directory + '/cn.mops_results.csv'
            if os.path.isfile(file_path):
                df_cnmops = pandas.read_csv(file_path,index_col=None)
        
        if directory == 'excavator2_output':
            file_path = './' + directory + '/excavator2_results.csv'
            if os.path.isfile(file_path):
                df_excavator2 = pandas.read_csv(file_path,index_col=None)

    if not os.path.exists('./transmission'):
        os.makedirs('./transmission')

### cn.mops
    family_cn = family(df_cnmops, pedigree)
    ol_cn = overlap(family_cn, df_cnmops)
    ol_cn.drop_duplicates(subset=['start', 'end'],\
        keep='first', inplace=True)
    final_cn = final(ol_cn)

    if os.path.isfile('./transmission/cn_transmission.csv'):
        os.remove('./transmission/cn_transmission.csv')
        print('\nPrevious cn_transmission.csv file removed.')
    
    final_cn.to_csv('transmission/cn_transmission.csv', index=False) 
    scores_cn = final_cn['score'].tolist()
    denovo_cn = scores_cn.count(1)
    print('CN.MOPS: Nombre de CNV de novo = {}'.format(denovo_cn))                                                                                           
    print("cn_transmission.csv generated.\n")

### CNVKIT    
    family_cnvkit = family(df_cnvkit, pedigree)
    ol_cnvkit = overlap(family_cnvkit, df_cnvkit)
    ol_cnvkit.drop_duplicates(subset=['start', 'end'],\
        keep='first', inplace=True)
    final_cnvkit = final(ol_cnvkit)

    if os.path.isfile('./transmission/cnvkit_transmission.csv'):
        os.remove('./transmission/cnvkit_transmission.csv')
        print('Previous cnvkit_transmission.csv file removed.')
    
    final_cnvkit.to_csv('transmission/cnvkit_transmission.csv', index=False) 
    scores_cnvkit = final_cnvkit['score'].tolist()
    denovo_cnvkit = scores_cnvkit.count(1)
    print('CNVKIT: Nombre de CNV de novo = {}'.format(denovo_cnvkit))                                                                                          
    print("cnvkit_transmission.csv generated.\n")

### exomedepth
    family_ed = family(df_exomedepth, pedigree)
    ol_ed = overlap(family_ed, df_exomedepth)
    ol_ed.drop_duplicates(subset=['start', 'end'],\
        keep='first', inplace=True)
    final_ed = final(ol_ed)

    if os.path.isfile('./transmission/ed_transmission.csv'):
        os.remove('./transmission/ed_transmission.csv')
        print('Previous ed_transmission.csv file removed.')
    
    final_ed.to_csv('transmission/ed_transmission.csv', index=False)  
    scores_ed = final_ed['score'].tolist()
    denovo_ed = scores_ed.count(1)
    print('EXOMEDEPTH: Nombre de CNV de novo = {}'.format(denovo_ed))                                                                                           
    print("ed_transmission.csv generated.\n")

### excavator2
    family_ex = family(df_excavator2, pedigree)
    ol_ex = overlap(family_ex, df_excavator2)
    ol_ex.drop_duplicates(subset=['start', 'end'],\
        keep='first', inplace=True)
    final_ex = final(ol_ex)

    if os.path.isfile('./transmission/ex_transmission.csv'):
        os.remove('./transmission/ex_transmission.csv')
        print('Previous ex_transmission.csv file removed.')
    
    final_ex.to_csv('transmission/ex_transmission.csv', index=False)
    scores_ex = final_ex['score'].tolist()
    denovo_ex = scores_ex.count(1)
    print('EXCAVATOR2: Nombre de CNV de novo = {}'.format(denovo_ex))                                                                                            
    print("ex_transmission.csv generated.\n")

### GATK
    family_gatk = family(df_gatk, pedigree)
    ol_gatk = overlap(family_gatk, df_gatk)
    ol_gatk.drop_duplicates(subset=['start', 'end'],\
        keep='first', inplace=True)
    final_gatk = final(ol_gatk)

    if os.path.isfile('./transmission/gatk_transmission.csv'):
        os.remove('./transmission/gatk_transmission.csv')
        print('Previous gatk_transmission.csv file removed.')
    
    final_gatk.to_csv('transmission/gatk_transmission.csv', index=False)  
    scores_gatk = final_gatk['score'].tolist()
    denovo_gatk = scores_gatk.count(1)
    print('GATK: Nombre de CNV de novo = {}'.format(denovo_gatk))                                                                                         
    print("gatk_transmission.csv generated.\n")


### HomeMade
    family_hm = family(df_hm, pedigree)
    ol_hm = overlap(family_hm, df_hm)
    ol_hm.drop_duplicates(subset=['start', 'end'],\
        keep='first', inplace=True)
    final_hm = final(ol_hm)

    if os.path.isfile('./transmission/hm_transmission.csv'):
        os.remove('./transmission/hm_transmission.csv')
        print('Previous hm_transmission.csv file removed.')
    
    final_hm.to_csv('transmission/hm_transmission.csv', index=False)  
    scores_hm = final_hm['score'].tolist()
    denovo_hm = scores_hm.count(1)
    print('HomeMade: Nombre de CNV de novo = {}'.format(denovo_hm))                                                                                         
    print("hm_transmission.csv generated.\n")
    
else:
    print('Le fichier pedigree.txt est absent. Calcul impossible.\n')

print("Transmission program job done!\n")