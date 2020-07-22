# coding: utf-8

# Author: Julie bogoin

import os

print("\nfemale_male_lists.py program openning.\n")

if os.path.isfile("female_list.txt"):
    print("Lists already generated. No change.\n")
else:
    males=[]
    females=[]
    if os.path.isfile('samples.txt'):
        with open ('samples.txt', 'r') as file:
            for line in file:
                line = file.readline()
                info = line[0].split(sep='\t')
                if info[1]=='M':
                    males.append(info[0])
                if info[1]=='F':
                    females.append(info[0])
        with open ('female_list.txt', 'w') as female_file:
            for female in females:
                female_file.write(female + '\n')
        with open ('male_list.txt', 'w') as male_file:
            for male in males:
                male_file.write(male + '\n')
    else:
        print("\nATTENTION! The samples.txt file is not present.\n \
        Launch sex_determination.py and start over!\n")
    print("female_male_lists.py job done!\n")
   
