#!/usr/bin/env python3

from Bio import SeqIO
import re
from optparse import OptionParser
import subprocess
import time

# defining the arguments which can be passed to the script
arguments = OptionParser()

arguments.add_option('-q', '--query', dest='query_seq', help='input query sequence file')# the query is the sequence which needs to be annotated
arguments.add_option('-i', '--input', dest='tabular_in', help='input tabular file for parsing')# the tabular input is the tabular output of the lastz program
arguments.add_option('-c', '--coding', dest='coding_table', help='input the coding table for the reference database')
arguments.add_option('-o', '--out', dest='output_name', help='name of file to which the bed output will be appended')
(options, args) = arguments.parse_args()
if (options.query_seq is None or options.tabular_in is None or options.coding_table is None or options.output_name is None): # if one of the arguments is not provided
        print('\n----------> A mandatory option is missing !\n')   # raise an error
        arguments.print_help()# and provide the help
        exit(-1) # exit the script

# defining the query
# the query will be a list of lists (for multiple query files)
query = list(SeqIO.parse(options.query_seq,'fasta'))

# everything needs to be done as many times as there are sequences in the query file

# the bed file will be created out of list of lists
# this list of lists needs to be outside the for loop because if must contain all hits to both queries


bed_file = open('./%s.bed' % (options.output_name), 'a') # append to file



# processing the coding table input
# the coding table file is read line by line
d_primary_key_upper = {}
d_primary_key_lower = {}
lst = []
with open(options.coding_table) as c:
        for line in c:
                items = line.split()
                if items[1] == 'F':
                        d1 = { items[0] : items[2] + '__' + items[0] }
                        d_primary_key_upper.update(d1)
                        lst.append(items[2] + '__' + items[0])
                        
                else:
                        d1 = { items[0] : items[2] + '__' + items[0] }
                        d_primary_key_lower.update(d1)
                        lst.append(items[2] + '__' + items[0])

for i in query: # looping over queries
        bed_data = [] # each query needs to have it's own list 
        d = {}
        for name in lst:
                d[name] = list(map(int, '0'*len(query[0])))

# processing the lastz tabular output line by line
        with open(options.tabular_in) as f:
                for line in f: # iterating over lines in file                        
                        if line[0] == '#': # if the lines start with a #, do not process this line
                                continue
                        items = line.split() # spliting the line into a list of strings
                        if i.id == items[0]: # if the first element of line match the name of the query
                                items[5] = re.sub('__.*', '', items[5]) # rename the hit name
                                start = int(items[2]) - 1 # calculating start of alignment on query
                                end = start + int(items[3]) # end of alignment on query
                                if items[9] == '+': # if the target element is '+' oriented
                                        d[d_primary_key_upper[items[5]]][start:end] = [i+1 for i in d[d_primary_key_upper[items[5]]][start:end]]
                                else: # if the target element is '-' oriented
                                  # use the primary key dictionary where pseudocoded names begin with a non capital letter
                                        d[d_primary_key_lower[items[5]]][start:end] = [i+1 for i in d[d_primary_key_lower[items[5]]][start:end]]
        sting = list(map(int, '0'*len(query[0]))) # strng of 0 which will contain the final scores of incrementation
        for key in d.keys(): # looping over keys in dictionary d
                for idx, val in enumerate(d[key]): # iterating over indexes and values for each key
                        if val == 0: # if the value is 0 there is no need to change anything in the final list
                                continue
                        else: # if the value is not 0
                                if sting[idx] == 0: # check if the value in the final list is 0
                                        sting[idx] = re.sub('__.*', '', key) # if yes assign the position in the final list as a corresponding alphabetical value
                                else:
                                        sting[idx] = 'X' # if not assign the conflict as a X
        name = i.id
        for idx, val in enumerate(sting): # zero based
                if val == 0:
                        continue
                if len(bed_data) == 0:
                        bed_row = [name, idx, idx, val]
                        bed_data.append(bed_row)
                if val == bed_data[len(bed_data)-1][3]:
                        if idx-bed_data[len(bed_data)-1][2] > 1:
                                bed_row = [name, idx, idx, val]
                                bed_data.append(bed_row)
                        else:
                                bed_data[len(bed_data)-1][2] = idx
                                
                if val != bed_data[len(bed_data)-1][3]:
                        bed_row = [name, idx, idx, val]
                        bed_data.append(bed_row)
        for row in bed_data:
                row[1] = str(row[1])
                row[2] = str(row[2])
                row = '\t'.join(row)
                bed_file.write(row + '\n')


# now because bed_data is a list of lists I will convert the lists into strings and write into file:

bed_file.close()# close the file

