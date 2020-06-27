# coding: utf-8

# Author: Julie BOGOIN

###### IMPORT ######

import pandas
import os
import random
from operator import attrgetter

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
    # - cnv_ratio
    # - copynumber
    # - sample
    # - sex
    # - cnv_tool
    # """

    def __init__(self, start, end, size, effect, log2copy_ratio, targets, contig, copynumber, sample, sex, cnv_tool):
        self.start = start
        self.end = end
        self.size = size
        self.effect = effect
        self.log2copy_ratio = log2copy_ratio
        self.targets = targets
        self.contig = contig
        self.copynumber = copynumber
        self.sample = sample
        self.sex = sex
        self.cnv_tool = cnv_tool
    
    @staticmethod
    def from_df_row(row):
        return Cnv(
            row['start'],
            row['end'],
            row['cnv_size'],
            row['effect'],
            row['log2copy_ratio'],
            row['targets_number'],
            row['contig'],
            row['CN'],
            row['sample'],
            row['sex'],
            row['cnv_tool'])
        
    @staticmethod
    def read_results(cnv_file):
        """Parse result file and return the list of Cnv objects
        """
        try:
            f = open("cnv_results.csv", 'r')
        except IOError:
            print('WARNING! cnv_results.csv does not exist!')
            print('CNV interval checking program aborted.\n')
        else:
            f.close()

        df = pandas.read_csv("cnv_results.csv", index_col=None, header=[0])

        # Supprimer les lignes identiques
        df.drop_duplicates(keep='first', inplace=True)

        # Trier par effect, sample, contig et position start.
        df.sort_values(by=['effect','sample', 'contig', 'start'], inplace=True)

        # Ajouter une colonne 'cnv_size'.
        df['cnv_size'] = df['end'] - df['start']
        
        # Supprimer le 'chr' de la colonne contig
        df['contig'] = df['contig'].str[3:]
 
        cnv_list = []
        for index, row in df.iterrows():
            cnv_list.append(Cnv.from_df_row(row))
     
        return cnv_list

"""
Functions
"""

def sort_in_run(sample_cnv_list):
    random.shuffle(sample_cnv_list)
    sample_cnv_list.sort(key=lambda cnv:cnv.effect)
    sample_cnv_list = sorted(sample_cnv_list, key=attrgetter('contig'))
    sample_cnv_list = sorted(sample_cnv_list, key=attrgetter('start'))
    
    return sample_cnv_list

def chevauchement_sample(sorted_sample_cnv):
    
    sample_cnv_list = []
    
    for i in range(1,len(sorted_sample_cnv)):
        prev_cnv = sorted_sample_cnv[i-1]
        curr_cnv = sorted_sample_cnv[i]

        if (curr_cnv.contig == prev_cnv.contig)\
            and (curr_cnv.effect == prev_cnv.effect)\
            and (int(curr_cnv.start) < int(prev_cnv.end)) and (int(curr_cnv.end) > int(prev_cnv.start)):

            min_start = min(int(prev_cnv.start), int(curr_cnv.start))
            max_end = max(int(prev_cnv.end), int(curr_cnv.end)) 

            contig = prev_cnv.contig
            start = min_start
            end = max_end
            size = max_end - min_start
            sample = curr_cnv.sample
            sex =  curr_cnv.sex
            effect = curr_cnv.effect
            cnv_tool = ",".join([prev_cnv.cnv_tool, curr_cnv.cnv_tool])
            copynumber = ",".join([str(prev_cnv.copynumber), str(curr_cnv.copynumber)])
            log2copy_ratio = ",".join([str(prev_cnv.log2copy_ratio), str(curr_cnv.log2copy_ratio)])
            targets = ",".join([str(prev_cnv.targets), str(curr_cnv.targets)])

            sample_cnv = Cnv(start, end, size, effect, log2copy_ratio, targets, contig, copynumber, sample, sex, cnv_tool)
            sample_cnv_list.append(sample_cnv)
                    
        else:
            sample_cnv_list.append(prev_cnv)

    return sample_cnv_list

def chevauchement_run(sorted_run_cnv):

    run_cnv_list = []

    for i in range(1, len(sample_cnv_list)):
        prev_cnv = sample_cnv_list[i-1]
        curr_cnv = sample_cnv_list[i]

        if (curr_cnv.contig == prev_cnv.contig) \
            and (int(curr_cnv.effect == prev_cnv.effect)) \
            and (int(curr_cnv.start) < int(prev_cnv.end)) \
            and (int(curr_cnv.end) > int(prev_cnv.start)):

            max_start = max(int(prev_cnv.start), int(curr_cnv.start))
            min_end = min(int(prev_cnv.end), int(curr_cnv.end)) 

            common_part = min_end - max_start
            overlap_on_prev_cnv = (float(common_part) / float(prev_cnv.size))
            overlap_on_curr_cnv = (float(common_part) / float(curr_cnv.size))
			
            if (overlap_on_prev_cnv > 0.8 and overlap_on_curr_cnv) > 0.8:
                    
                contig = prev_cnv.contig
                start = max_start
                end = min_end
                size = min_end - max_start
                sample = curr_cnv.sample
                sex =  curr_cnv.sex
                effect = curr_cnv.effect

                cnv_tool = ",".join([prev_cnv.cnv_tool, curr_cnv.cnv_tool])
                copynumber = ",".join([str(prev_cnv.copynumber), str(curr_cnv.copynumber)])
                log2copy_ratio = ",".join([str(prev_cnv.log2copy_ratio), str(curr_cnv.log2copy_ratio)])
                targets = ",".join([str(prev_cnv.targets), str(curr_cnv.targets)])

                run_cnv = Cnv(start, end, size, effect, log2copy_ratio, targets, contig, copynumber, sample, sex, cnv_tool)

                run_cnv_list.append(run_cnv)

        else:
            run_cnv_list.append(prev_cnv)        

    return run_cnv_list
    

###### PROGRAMME PRINCIPAL ######

print("\nCNV interval checking program openning.\n")

print('\n-START Generation des cnv...')
cnv = Cnv.read_results('cnv_results.csv')
cnv_count = len(cnv)
print('-END Generation des cnv: {} cnv generes'.format(cnv_count))

sample_cnv_list = chevauchement_sample(cnv)

sorted_run_cnv = sort_in_run(sample_cnv_list)

run_cnv_list = chevauchement_run(sorted_run_cnv)

print('{} CNV intervals found among {} CNV lines.'.format(len(run_cnv_list), len(sorted_run_cnv)))

if os.path.isfile('interval_results.txt'):
    os.remove('interval_results.txt')
    print('Previous results file removed.')


with open('interval_results.txt', 'w') as results_file:
    
    results_file.write('contig\tstart\tend\treference\tobserved\tsample\tsex\tsize\teffect\tlog2copy_ratio\tcopynumber\tcnv_tool\ttargets_number\n')
    
    for run_cnv in run_cnv_list:

        if run_cnv.contig == 'XF' or run_cnv.contig == 'XM':
            run_cnv.contig = 'X'

        results_file.write(str(run_cnv.contig) + '\t'
                           + str(run_cnv.start) + '\t'
                           + str(run_cnv.end) + '\t'
                           + str('0') + '\t'
                           + str('0') + '\t'
                           + str(run_cnv.sample) + '\t'
                           + str(run_cnv.sex) + '\t'
                           + str(run_cnv.size) + '\t'
                           + str(run_cnv.effect) + '\t'
                           + str(run_cnv.log2copy_ratio) + '\t'
                           + str(run_cnv.copynumber) + '\t'
                           + str(run_cnv.cnv_tool) + '\t'
                           + str(run_cnv.targets) + '\t'
                           + '\n')   

print("interval_results.txt generated.\n")

print("CNV interval checking: job done!\n")
