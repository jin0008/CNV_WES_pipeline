# coding: utf-8

# Author: Julie bogoin

import os

print("\nfemale_male_lists.py program openning.\n")

if os.path.isfile("female_list.txt"):
    print("Listes deja generees. Pas de changement.\n")
else:

    males=[]
    females=[]

    if os.path.isfile('samples.txt'):
        with open ('samples.txt', 'r') as file:
            for line in file:
                line = line.split(sep='\n')
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
        print("\nATTENTION! Le fichier samples.txt n'est pas present.\n \
        Lancez sex_determination.py et recommencez!\n")

    print("female_male_lists.py job done!\n")
    


