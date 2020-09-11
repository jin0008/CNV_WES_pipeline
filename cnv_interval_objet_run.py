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

    def __init__(self, sample, contig, start, end, sex, size, effect, log2copy_ratio, copynumber, cnv_tool, targets_number):
        self.sample = sample
        self.contig = contig
        self.start = start
        self.end = end
        self.sex = sex
        self.size = size
        self.effect = effect
        self.log2copy_ratio = log2copy_ratio
        self.copynumber = copynumber  
        self.cnv_tool = cnv_tool
        self.targets = targets_number

    @staticmethod
    def from_df_row(row):
        return Cnv(
            row['sample'],
            row['contig'],
            row['start'],
            row['end'],
            row['sex'],
            row['size'],
            row['effect'],
            row['log2copy_ratio'],
            row['copynumber'],
            row['cnv_tool'],
            row['targets_number'])
        
    @staticmethod
    def read_results(cnv_file):
        """Parse result file and return the list of Cnv objects
        """
        try:
            f = open(cnv_file, 'r')
        except IOError:
            print('WARNING! interval_sample_results.csv does not exist!')
            print('CNV interval checking program aborted.\n')
        else:
            f.close()

        df = pandas.read_csv(cnv_file, sep='\t', index_col=False)
        
        # Trier par effect, contig et position start.
        df.sort_values(by=['effect', 'contig', 'start'], inplace=True)

        cnv_list = []
        for index, row in df.iterrows():
            cnv_list.append(Cnv.from_df_row(row))
     
        return cnv_list

"""
Functions
"""

def chevauchement_run(sorted_sample_cnv):
    
    run_cnv_list = []

    prev_cnv = sorted_sample_cnv[0]
    curr_cnv = sorted_sample_cnv[1]

    sample_list = []
    contig_list = []
    start_list = []
    end_list = []
    sex_list = []
    size_list = []
    effect_list = []
    log2copy_list = []
    copynumber_list = []
    tool_list = []
    targets_list = []

    for i in range(0, len(sorted_sample_cnv)-1):
        
        if (curr_cnv.effect == prev_cnv.effect) \
        and (curr_cnv.contig == prev_cnv.contig) \
        and (int(curr_cnv.start) <= int(prev_cnv.end)) \
        and (int(curr_cnv.end) >= int(prev_cnv.start)):
            
            max_start = max(int(prev_cnv.start), int(curr_cnv.start))
            min_end = min(int(prev_cnv.end), int(curr_cnv.end)) 
            min_start = min(int(prev_cnv.start), int(curr_cnv.start))
            max_end = max(int(prev_cnv.end), int(curr_cnv.end))
            print('\nmax start = {}'.format(max_start))
            print('min end = {}'.format(min_end))
            
            common_part = min_end - max_start
            print('common part = {}'.format(common_part))
            
            overlap_on_prev_cnv = (float(common_part) / float(prev_cnv.size))
            overlap_on_curr_cnv = (float(common_part) / float(curr_cnv.size))
	    print('overlap on prev cnv = {}'.format(overlap_on_prev_cnv))
            print('overlap on curr cnv = {}'.format(overlap_on_curr_cnv))
            
            if (overlap_on_prev_cnv > 0.8) and (overlap_on_curr_cnv > 0.8):

                sample_list.append(prev_cnv.sample)
                contig_list.append(prev_cnv.contig)
                start_list.append(min_start)
                end_list.append(max_end)
                sex_list.append(prev_cnv.sex)
                size_list.append(max_end - min_start)
                effect_list.append(prev_cnv.effect)
                log2copy_list.append(prev_cnv.log2copy_ratio)
                copynumber_list.append(prev_cnv.copynumber)
                tool_list.append(prev_cnv.cnv_tool)
                targets_list.append(prev_cnv.targets)
    
                prev_cnv = sorted_sample_cnv[i]
                curr_cnv = sorted_sample_cnv[i+1]
                print('... fusion des intervalles ... sample = {}'.format(prev_cnv.sample))

            else:
                run_cnv_list.append(prev_cnv)    
                prev_cnv = sorted_sample_cnv[i]
                curr_cnv = sorted_sample_cnv[i+1]
                print('\nchevauchement sans overlap ligne = {}'.format(i))
        
        else:
            run_cnv_list.append(prev_cnv)
            prev_cnv = sorted_sample_cnv[i]
            curr_cnv = sorted_sample_cnv[i+1]
            
            print('\nligne solo = {}'.format(i))
            
            for j in range(0, len(sample_list)):
                
                run_cnv = Cnv(sample_list[j], contig_list[j], start_list[j], end_list[j], \
                sex_list[j], size_list[j], effect_list[j], log2copy_list[j], copynumber_list[j], \
                tool_list[j], targets_list[j])

                run_cnv_list.append(run_cnv)

            sample_list = []
            contig_list = []
            start_list = []
            end_list = []
            sex_list = []
            size_list = []
            effect_list = []
            log2copy_list = []
            copynumber_list = []
            tool_list = []
            targets_list = []

    return run_cnv_list
    

###### PROGRAMME PRINCIPAL ######

print("\nRun intervals checking program openning.\n")

print('-START Generation des cnv...')
cnv = Cnv.read_results('interval_sample_results.txt')
cnv_count = len(cnv)
print('-END Generation des cnv: {} cnv generes'.format(cnv_count))

run_cnv_list = chevauchement_run(cnv)

if os.path.isfile('interval_run_results.txt'):
    os.remove('interval_run_results.txt')
    print('\nPrevious results file removed.')

with open('interval_run_results.txt', 'w') as results_file:
    
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

print("interval_run_results.txt generated.\n")

print("Run intervals checking: job done!\n")