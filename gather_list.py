import os
import sys

in_file = sys.argv[1]
wd = sys.argv[2]
#wd = os.getcwd()

gather_list = [f'{wd}/Demultiplex_Output']

with open(in_file, 'r') as file:
    for line in file:
        gather_list.append(f'{wd}/analysis/{line}')

print(' '.join(gather_list))
