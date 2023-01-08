#!/usr/bin/python

"""
sipros_prepare_protein_database.py

sipros_prepare_protein_database.py is used created a protein database with decoy sequence
appended

Created by Xuan Guo on 02/20/2017.
Copyright (c) 2017 Xuan Guo (ORNL). Allrights reserved.
"""

## Import Python package modules
import sys, getopt, os, random
from datetime import datetime, date, time

## Import Sipros package modules
import sipros_post_module
import parseconfig

## Check file exist
check_file_exist = sipros_post_module.check_file_exist

## Returns the current time in a nice format
curr_time = sipros_post_module.curr_time

## Format time as a pretty string
format_time = sipros_post_module.format_time

def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], "hi:o:c:")

    output_filename = ""
    input_filename  = ""
    config_filename = ''

    # Basic options
    for option, value in opts:
        if option in ("-h"):
            print('sipros_prepare_protein_database.py -i input-file -o output-file -c config-file')
            sys.exit(1)
        if option in ("-i"):
            input_filename = value
        elif option in ("-o"):
            output_filename = value
        elif option in ('-c'):
            config_filename = value
        else:
            print('sipros_prepare_protein_database.py -i input-file -o output-file -c config-file')
            sys.exit(1) 

    if input_filename == '' or output_filename == '' or config_filename == '' :
        print('sipros_prepare_protein_database.py -i input-file -o output-file -c config-file')
        sys.exit(1)
    '''
    if (output_filename == "") :
        (inputFileNameRoot, inputFileNameExt) = os.path.splitext(input_filename)
        output_filename = inputFileNameRoot + "_CFR" + inputFileNameExt
    '''
        
    print('Command:')
    print('sipros_prepare_protein_database.py -i {} -o {} -c {}'.format(input_filename, output_filename, config_filename))
        
    return input_filename, output_filename, config_filename


def reverse_protein_database(input_file_str, output_file_str, all_config_dict) :
    
    probability_1 = 0.5
    probability_2 = 1
    
    training_prefix_str = 'Rev1_'
    testing_prefix_str = 'TestRev_'
    reserved_prefix_str = 'Rev2_'
    
    if '[Protein_Identification]Reserved_Decoy_Prefix' in all_config_dict:
        probability_1 = 0.333333
        probability_2 = 0.666666
        reserved_prefix_str = all_config_dict['[Protein_Identification]Reserved_Decoy_Prefix']
    else:
        probability_2 = 1
    
    if reserved_prefix_str == 'SIP': # no training is needed for SIP search, so only test decoy is generated
        probability_1 = -1
        
    if '[Protein_Identification]Training_Decoy_Prefix' in all_config_dict:
        training_prefix_str = all_config_dict['[Protein_Identification]Training_Decoy_Prefix']
    else:
        print('Value for [Protein_Identification]Training_Decoy_Prefix is missing.')
        exit(1)
        
    if '[Protein_Identification]Testing_Decoy_Prefix' in all_config_dict:
        testing_prefix_str = all_config_dict['[Protein_Identification]Testing_Decoy_Prefix']
    else:
        print('Value for [Protein_Identification]Testing_Decoy_Prefix is missing.')
        exit(1)
    
    output_file = open(output_file_str, "w")
    id_str = ""
    seq_str = ""
    seq_new_str = ""
    input_file = open(input_file_str, "r")
    line_str = ""
    random_float = 0.0
    # avoid different result when running again
    random.seed(9527)
    for line_str in input_file:
        if line_str[0] == '>':
            if seq_str != "":
                seq_new_str = (seq_str[::-1])
                output_file.write(id_str)
                output_file.write(seq_str)
                output_file.write('\n')
                random_float = random.random()
                if random_float <= probability_1:
                    output_file.write('>')
                    output_file.write(training_prefix_str)
                elif random_float <= probability_2:
                    output_file.write('>')
                    output_file.write(testing_prefix_str)
                else:
                    output_file.write('>')
                    output_file.write(reserved_prefix_str)                    
                output_file.write(id_str[1:])
                output_file.write(seq_new_str)
                output_file.write("\n")
            id_str = line_str
            seq_str = ""
        else:
            seq_str += line_str.strip()

    if seq_str != "":
        seq_new_str = (seq_str[::-1])
        output_file.write(id_str)
        output_file.write(seq_str)
        output_file.write('\n')
        random_float = random.random()
        if random_float <= probability_1:
            output_file.write('>')
            output_file.write(training_prefix_str)
        elif random_float <= probability_2:
            output_file.write('>')
            output_file.write(testing_prefix_str)
        else:
            output_file.write('>')
            output_file.write(reserved_prefix_str)                    
        output_file.write(id_str[1:])
        output_file.write(seq_new_str)
        output_file.write("\n")
    id_str = line_str
    seq_str = ""
        
    input_file.close()
    output_file.close()

## Parse config file
def parse_config(config_filename):

    # Save all config values to dictionary
    all_config_dict = {}    # initialize dictionay
    # Save config values to dictionary
    config_dict = {}    # initialize dictionay

    # Call Yinfeng's parseconfig.py module
    check_file_exist(config_filename)
    all_config_dict = parseconfig.parseConfigKeyValues(config_filename)
    
    return all_config_dict

def main(argv=None):

    if argv is None:
        argv = sys.argv
        
    # parse options
    input_file_str, output_file_str, config_file_str = parse_options(argv)
    
    # Display work start and time record
    start_time = datetime.now()
    sys.stderr.write('[%s] Beginning sipros_prepare_protein_database.py\n' % (curr_time()))
    sys.stderr.write('------------------------------------------------------------------------------\n')    

    # parse config file
    sys.stderr.write('[Step 1] Parse options:                           Running -> ')
    all_config_dict = parse_config(config_file_str)
    sys.stderr.write('Done!\n')
    
    # reverse sequence and save file
    sys.stderr.write('[Step 2] Generate new database:                   Running -> ')
    reverse_protein_database(input_file_str, output_file_str, all_config_dict)
    sys.stderr.write('Done!\n')
    
    # Time record, calculate elapsed time, and display work end
    finish_time = datetime.now()
    duration = finish_time - start_time
    sys.stderr.write('------------------------------------------------------------------------------\n')
    sys.stderr.write('[%s] Ending sipros_prepare_protein_database.py \n' % curr_time())
    sys.stderr.write('Run complete [%s elapsed]\n' %  format_time(duration))
    
## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    sys.exit(main())
