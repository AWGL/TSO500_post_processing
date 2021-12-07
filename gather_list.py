"""
Prints a specifically formatted list of samples for input to the Illumina app.
Must start with the Demultiplex_Output folder and all folders must have the full file path

"""
import os
import sys

# load in system arguments
in_file = sys.argv[1]
wd = sys.argv[2]

# list must start with Demultiplex_Output folder
gather_list = [f'{wd}/Demultiplex_Output']

# add each sample analysis folder to the end of the list
with open(in_file, 'r') as file:
    for line in file:
        gather_list.append(f'{wd}/analysis/{line}')

# print the whole list seperated by spaces
print(' '.join(gather_list))
