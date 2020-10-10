# Author: Julie BOGOIN

###### IMPORT ######

import pandas
import os

# """
# Classes
# """

class Cnv:
    # """ cnv class containing one cnv's following properties:
    # - start
    # - end
    # - size
    # - effect
    # - log2copy_ratio
    # - targets
    # - contig
    # - copynumber
    # - sample
    # - sex
    # - cnv_tool
    # """

    def __init__(self, sample, sex, contig, start, end, size, log2copy_ratio, copynumber, effect,  cnv_tool, targets):
        self.sample = sample
        self.sex = sex
        self.contig = contig
        self.start = start
        self.end = end
        self.size = size
        self.log2copy_ratio = log2copy_ratio
        self.copynumber = copynumber
        self.effect = effect
        self.cnv_tool = cnv_tool
        self.targets = targets
    
    @staticmethod
    def from_df_row(row):
        return Cnv(
            row['sample'],
            row['sex'],
            row['contig'],
            row['start'],
            row['end'],
            row['size'],
            row['log2copy_ratio'],
            row['CN'],
            row['effect'],
            row['cnv_tool'],
            row['targets_number'])
        
    @staticmethod
    def read_results(cnv_file):
        """Parse result file and return the list of Cnv objects
        """
        try:
            f = open(cnv_file, 'r')
        except IOError:
            print('WARNING! cnv_results.csv does not exist!')
            print('CNV interval checking program aborted.\n')
        else:
            f.close()

        df = pandas.read_csv(cnv_file, header=[0])

        # Supprimer les lignes identiques
        df.drop_duplicates(keep='first', inplace=True)

        # Trier par effect, sample, contig et position start.
        df.sort_values(by=['effect','sample', 'contig', 'start', 'end'], inplace=True, \
            ascending = [True, True, True, True, False], kind = 'mergesort')

        # Ajouter une colonne 'size'.
        df['size'] = df['end'] - df['start']
        
        # Supprimer le 'chr' de la colonne contig
        df['contig'] = df['contig'].str[3:]

        cnv_list = []
        for index, row in df.iterrows():
            cnv_list.append(Cnv.from_df_row(row))
     
        return cnv_list

"""
Functions
"""

def chevauchement_sample(sorted_sample_cnv):
    
    sample_cnv_list = []

    prev_cnv = sorted_sample_cnv[0]
    curr_cnv = sorted_sample_cnv[1]
    
    for i in range(1,len(sorted_sample_cnv)-1):

        if (curr_cnv.contig == prev_cnv.contig)\
            and (curr_cnv.effect == prev_cnv.effect)\
            and (curr_cnv.sample == prev_cnv.sample)\
            and (int(curr_cnv.start) < int(prev_cnv.end)) \
            and (int(curr_cnv.end) > int(prev_cnv.start)):

            min_start = min(int(prev_cnv.start), int(curr_cnv.start))
            max_end = max(int(prev_cnv.end), int(curr_cnv.end)) 

            contig = prev_cnv.contig
            start = min_start
            end = max_end
            size = max_end - min_start
            sample = curr_cnv.sample
            sex =  curr_cnv.sex
            effect = curr_cnv.effect
            cnv_tool = "/".join([prev_cnv.cnv_tool, curr_cnv.cnv_tool])
            copynumber = "/".join([str(prev_cnv.copynumber), str(curr_cnv.copynumber)])
            log2copy_ratio = "/".join([str(prev_cnv.log2copy_ratio), str(curr_cnv.log2copy_ratio)])
            targets = "/".join([str(prev_cnv.targets), str(curr_cnv.targets)])

            sample_cnv = Cnv(sample, sex, contig, start, end, size, \
                log2copy_ratio, copynumber, effect,  cnv_tool, targets)
            prev_cnv = sample_cnv
            curr_cnv = sorted_sample_cnv[i+1]
                    
        else:
            sample_cnv_list.append(prev_cnv)
            prev_cnv = sorted_sample_cnv[i]
            curr_cnv = sorted_sample_cnv[i+1]

    return sample_cnv_list


###### PROGRAMME PRINCIPAL ######

print("\n************************************")
print("Sample intervals checking program openning.")
print("************************************\n")

print('-START Generation des cnv...')
cnv = Cnv.read_results('cnv_results.csv')
cnv_count = len(cnv)
print('-END Generation des cnv: {} cnv generes'.format(cnv_count))

sample_cnv_list = chevauchement_sample(cnv)

print('{} chevauchement(s) trouves.'.format(len(sample_cnv_list) - len(cnv)))
print('-{} cnv(s) au final.'.format(len(sample_cnv_list)))

if os.path.isfile('interval_sample_results.txt'):
    os.remove('interval_sample_results.txt')
    print('\nPrevious results file removed.')

with open('interval_sample_results.txt', 'w') as results_file:
    
    results_file.write('sample\tcontig\tstart\tend\tsex\tsize\teffect\tlog2copy_ratio\tcopynumber\tcnv_tool\ttargets_number\n')
    
    for run_cnv in sample_cnv_list:

        if run_cnv.contig == 'XF' or run_cnv.contig == 'XM':
            run_cnv.contig = 'X'

        results_file.write(str(run_cnv.sample) + '\t'
                           + str(run_cnv.contig) + '\t'
                           + str(run_cnv.start) + '\t'
                           + str(run_cnv.end) + '\t'
                           + str(run_cnv.sex) + '\t'
                           + str(run_cnv.size) + '\t'
                           + str(run_cnv.effect) + '\t'
                           + str(run_cnv.log2copy_ratio) + '\t'
                           + str(run_cnv.copynumber) + '\t'
                           + str(run_cnv.cnv_tool) + '\t'
                           + str(run_cnv.targets) + '\t'
                           + '\n')   

print("interval_sample_results.txt generated.\n")

print("Sample intervals checking: job done!\n")
