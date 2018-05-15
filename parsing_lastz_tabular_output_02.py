#!/usr/bin/env python3

from optparse import OptionParser
import re

# defining the arguments that need to be passed to the script
arguments = OptionParser()

arguments.add_option('-i', '--tabular', dest='lastz_tabular_input', help='input tabular file for parsing')
arguments.add_option('-c', '--coding', dest='coding_table', help='input the coding table')
arguments.add_option('-o', '--output', dest='output_name', help='output name')
(options, args) = arguments.parse_args()
if (options.lastz_tabular_input is None or options.coding_table is None or options.output_name is None): # if one of the arguments is missing
    print('\n----------> A mandatory option is missing !\n')  # raise an error
    arguments.print_help()  # and provide the help
    exit(-1)  # exit the script
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# parsing the coding table file
d_primary_key_upper = {} # dictionary which contains names of repeats as keys and uppercase pseudocodes as values
d_primary_key_lower = {} # dictionary which contains names of repeats as keys and lowercase pseudocodes as values

with open(options.coding_table) as f: # opening the file of the coding table
    for line in f: # iterating over all lines
        items = line.split() # splitting the lines into items
        if items[1] == 'F': # if the repeat is forward oriented assign it to d_primary_key_upper
            d_primary_key_upper[items[0]] = items[2]
        else: # if the key is reverse oriented assign it to d_primary_key_lower
            d_primary_key_lower[items[0]] = items[2]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# creating functions which will be called in executing the main function of the script
# function create_new_dictionary creates the initial dictionary whith strings of zeos and increments based on the first line in the
# lastz tabular output

# q_length <- query length, a_length <- length of alignment, start <- position of the start of the alignment
# repeat_name <- name of the repeat in the line that is being processed, orientation <- orientation of the repeat
# upper_coding_table <- coding table with the upper case pseudocodes, lower_coding_table <- lower case pseudocodes
def create_new_dictionary(q_length, start, a_length, repeat_name, orientation, upper_coding_table, lower_coding_table): 
    d = {} # creating the dictionary d which has the pseudocodes as keys and list of zero as value
    for key in upper_coding_table: # iterating over keys
        pseudocode_upper = upper_coding_table[key] 
        pseudocode_lower = lower_coding_table[key]
        d[pseudocode_upper] = list(map(int, '0'* int(q_length)))
        d[pseudocode_lower] = list(map(int, '0'* int(q_length)))
    start_increment = int(start) - 1 # calculating start of th alignment
    end_increment = start_increment + int(a_length) # calculating the position of the end

    if orientation == '+': # if the repeat is forward oriented
        key_for_d = upper_coding_table[repeat_name] # d_primary_key_upper is used to get the pseudocode which is a key in dictionary d
        d[key_for_d][start_increment:end_increment] = [i+1 for i in d[key_for_d][start_increment:end_increment]] # incrementing 
    else: # if the repeat is reverse orientted
        key_for_d = lower_coding_table[repeat_name] # # d_primary_key_lower is used to get the pseudocode which is a key in dictionary d
        d[key_for_d][start_increment:end_increment] = [i+1 for i in d[key_for_d][start_increment:end_increment]] # incrementing
    return (d)  # output of this function is the newly created d dictionary


new_dictionary = create_new_dictionary(20,1,5,'SSR_TAA','+',d_primary_key_upper,d_primary_key_lower)

    
# this function is used when the name of the read in the line that is being processed is the same as the read name od the line
# that has been previously processed
# a_length <- alignment length, start <- position of the start of the alignment,
# repeat_name <- name of the repeat in the line that is being processed, orientation <- orientation of the repeat
# upper_coding_table <- coding table with the upper case pseudocodes, lower_coding_table <- lower case pseudocodes
def continue_incrementing(start, a_length, repeat_name, orientation, d, upper_coding_table, lower_coding_table):

    start_increment = int(start)-1
    end_increment = start_increment + int(a_length)

    if orientation == '+':
        key_for_d = upper_coding_table[repeat_name]
        d[key_for_d][start_increment:end_increment] = [i+1 for i in d[key_for_d][start_increment:end_increment]]
    else:
        key_for_d = lower_coding_table[repeat_name]
        d[key_for_d][start_increment:end_increment] = [i+1 for i in d[key_for_d][start_increment:end_increment]]
        
    return(d)

# in the case where the the name of the read in the line of the lastz output does not match the one that has been previously processed
# the existing dictionay d must the processed into a continuous pseudocoded string
# q_length <- query  length, d <- existing d dictionary
def create_pseudocode_string(q_length, d):
    sting = list(map(str, '0'* int(q_length))) # creating a list of 0 the length of the query
    for key in d.keys(): # iterating over keys in d
        for idx, val in enumerate(d[key]): # iterating over all values and indices of the values in each key of d
            if val == 0: # if the value is zero withing the dictionary, skip over it
                continue
            elif sting[idx] == '0': # if it is not, and the value of the string is 0, assign a new value that is
                # the same as the key that is currently being processed
                sting[idx] = key
            else:
                # if the value in the string is not 0, then it is a conflict, therefore assign value of X
                sting[idx] = 'X'
    return(sting) # the output of this function is is the pseudocoded string

# this function creates the bed file based on the pseudocoded string
# read_name <- name of the read of the existing pseudocoded string, pseudocode_string <- existing continuous pseudocoded string
def create_bed_file(read_name, pseudocode_string):
    bed_data = []
    for idx, val in enumerate(pseudocode_string): # iterating over all values and indices of the string
        if val == '0': # if the value is equal to 0, skip over it
            continue
        if len(bed_data) == 0: # if the bed list is empty add first value that is not 0
            bed_row = [read_name, idx, idx, val]
            bed_data.append(bed_row)
        if bed_data[-1][3] == val: # if the value is the same pseudocode as the last value in the bed list
        # the check for a gap
        # if there is no gap assign a new value for the end of the alignment in the bed list
            if idx - bed_data[-1][2] == 1 or idx - bed_data[-1][2] == 0:
                bed_data[-1][2] = idx
            else:
                bed_row = [read_name, idx, idx, val]
                bed_data.append(bed_row)
        else:
            bed_row = [read_name, idx, idx, val]
            bed_data.append(bed_row)

    for row in bed_data:
        row[1] = str(row[1])
        row[2] = str(row[2])
        row = '\t'.join(row)
        bed_file.write(row + '\n')

    return(bed_data) # the output is the bed list
    

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# creating the main function of the script which will call the rest of the functions
bed_file = open('./%s' % (options.output_name), 'a') # file in which the bed data will bewritten
r_name = '' # variabe that stores the name of the read

with open(options.lastz_tabular_input) as t: # opening lastz output
    for line in t: # iterating over lines
        items = line.split() # splitting columns of lines
        items[5] = re.sub('__.*', '', items[5]) # modifying names of reads
        if len(r_name) == 0:  # if the r_name variable is empty, process the first line on the tabular intput
            r_name = items[0] # assign a new value to r_name
            dictionary = create_new_dictionary(items[1], items[2], items[3], items[5], items[9], d_primary_key_upper, d_primary_key_lower) # create the dictionary
            # which has pseudocodes as keys and lists of 0  with incrementation as values
        if r_name == items[0]: # if the name of the read that is currently being processed
            # and the previously processed one, continue with the incrementation
            dictionary = continue_incrementing(items[2], items[3], items[5], items[9], dictionary, d_primary_key_upper, d_primary_key_lower)
        if r_name != items[0]: # if they are not the same
            string = create_pseudocode_string(items[1], dictionary) # create the pseudocode string
            bed = create_bed_file(r_name, string) # creat the bed file from the string
            r_name = items[0] # assign a new value to r_name
            dictionary = create_new_dictionary(items[1], items[2], items[3], items[5], items[9], d_primary_key_upper, d_primary_key_lower) # create a new dictionary


# these two lines were added because when the loop arrives to the final read
# it ends the loop without creating the pseudocoded string and bed file because
# there is no more sequence to which to compare the last read 
string = create_pseudocode_string(items[1], dictionary) # create string
bed  = create_bed_file(r_name, string) # create bed and write to file 


bed_file.close()
        
