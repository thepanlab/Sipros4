#!/usr/bin/python

"""
sipros_peptides_assembling.py

sipros_peptides_assembling is the second post-processing program
after running of sipros_peptides_filtering.py for assembling 
peptides to identify proteins.

Created by Tae-Hyuk (Ted) Ahn on 10/10/2012.
Copyright (c) 2012 Tae-Hyuk Ahn (ORNL). Allrights reserved.
"""

## Import standard modules
from __future__ import division
import sys, getopt, warnings, os, re
from datetime import datetime, date, time
from collections import defaultdict
import csv 

## Import Sipros package modules
import sipros_post_module
import parseconfig

## increase CSV field size limit
csv.field_size_limit(1000000000)

## Version control
def get_version():
    return "4.0.1 (Alpha)"
"""
1. for Sipros4.0 (added two columns in sip file)
   new program should handle both types of sip files
"""



## Import classes
## Class Usage
Usage = sipros_post_module.Usage

## Class for ignoring comments '#' in sipros file
CommentedFile = sipros_post_module.CommentedFile

## Exit system with error message
die = sipros_post_module.die

## Returns the current time in a nice format
curr_time = sipros_post_module.curr_time

## Format time as a pretty string
format_time = sipros_post_module.format_time

## Check file exist
check_file_exist = sipros_post_module.check_file_exist

## Get file(s) list in working dir with specific file extension
get_file_list_with_ext = sipros_post_module.get_file_list_with_ext

## Get base_out filename
get_base_out = sipros_post_module.get_base_out

## Class for PepOutFields object
PepOutFields = sipros_post_module.PepOutFields

## Class for PsmOutFields object
PsmOutFields = sipros_post_module.PsmOutFields
Psm4OutFields = sipros_post_module.Psm4OutFields

## list_to_string
list_to_string = sipros_post_module.list_to_string

## Find string between two substrings
find_between = sipros_post_module.find_between

## Division error handling
divide = sipros_post_module.divide

## PrettyFloat (%0.5f)
PrettyFloat = sipros_post_module.PrettyFloat

## Get item list
check_sub_list = sipros_post_module.check_sub_list

## Get item list
get_item_list = sipros_post_module.get_item_list

## set_float_digit
set_float_digit = sipros_post_module.set_float_digit


## Help message
help_message = '''

Usage:
    python sipros_peptides_assembling.py [options]

Inputs:
    sipros_peptides_filtering output: [*.pep.txt] and [*.psm.txt] files
        (multiple peptide filtering results can be processed)
        (search automatically in current directory)
    sipros config file    (search automatically in current directory)

Options:
    -h/--help
    -v/--version
    -w/--working-dir ./path/    # Directory path containing SIPROS output files
                                # (default = current directory) 
    -c/--config-file SiprosConfig.cfg    # SIPROS config file 
                                         # (default = SiprosConfig.cfg) 

Outputs:
    BaseFilename.pro.txt
    BaseFilename.pro2pep.txt
    BaseFilenameepro2psm.txt
        - where BaseFilename = common prefix of inputs
        - if len(BaseFilename) < 5, then BaseFilename = Sipros_searches
'''

## Glboal variables
pep_file_ext = '.pep.txt'
psm_file_ext = '.psm.txt'

pep_iden_str = '[Peptide_Identification]'
fasta_database_str = 'FASTA_Database'
pro_iden_str = '[Protein_Identification]'
decoy_prefix_str = 'Decoy_Prefix'
min_peptide_per_protein_str = 'Min_Peptide_Per_Protein'
min_unique_peptide_per_protein_str = 'Min_Unique_Peptide_Per_Protein'
remove_decoy_identification_str = 'Remove_Decoy_Identification'

sipros4_psmout_column_length = 20
sipros4_input = None

# defaul value
decoy_prefix = 'Rev_'
min_peptide_per_protein = 2
min_unique_peptide_per_protein = 1
remove_decoy_identification = 'No'


## Parse options
def parse_options(argv):

    try:
        opts, args = getopt.getopt(argv[1:], "hvVw:c:",
                                    ["help",
                                     "version",
                                     "working-dir",
                                     "config-file"])

    # Error handling of options
    except getopt.error, msg:
        raise Usage(msg)

    # Default working dir and config file
    working_dir = "./"
    config_file = "SiprosConfig.cfg"

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            raise Usage(help_message)
        if option in ("-v", "-V", "--version"):
            print "sipros_peptides_assembling.py V%s" % (get_version())
            sys.exit(0)
        if option in ("-w", "--working-dir"):
            working_dir = value
            if working_dir[-1] != '/':
                working_dir = working_dir + '/'
        if option in ("-c", "--config-file"):
            config_file = value

    # only -w is provided
    if working_dir != "./" and config_file == "SiprosConfig.cfg":
        config_file = working_dir + config_file

    return (working_dir, config_file)


## Parse config file
def parse_config(config_filename):

    # Save all config values to dictionary
    all_config_dict = {}    # initialize dictionay
    # Save config values to dictionary
    config_dict = {}    # initialize dictionay

    # Call Yinfeng's parseconfig.py module
    check_file_exist(config_filename)
    all_config_dict = parseconfig.parseConfigKeyValues(config_filename)

    # only save protein_identification config info to config_dict
    config_dict[decoy_prefix_str] = decoy_prefix
    config_dict[min_peptide_per_protein_str] = min_peptide_per_protein
    config_dict[min_unique_peptide_per_protein_str] = min_unique_peptide_per_protein
    config_dict[remove_decoy_identification_str] = remove_decoy_identification
    for key, value in all_config_dict.items():
        if key == (pep_iden_str + fasta_database_str):
            config_dict[fasta_database_str] = value
        elif key == (pro_iden_str + decoy_prefix_str):
            config_dict[decoy_prefix_str] = value
        elif key == (pro_iden_str + min_peptide_per_protein_str):
            config_dict[min_peptide_per_protein_str] = value
        elif key == (pro_iden_str + min_unique_peptide_per_protein_str):
            config_dict[min_unique_peptide_per_protein_str] = value
        elif key == (pro_iden_str + remove_decoy_identification_str):
            config_dict[remove_decoy_identification_str] = value
        else:
            continue

    # return config dictionary
    return config_dict


## Read fasta file and save the description
def read_fasta_file(working_dir, config_dict):

    # get fasta filename from config file
    fasta_filename = config_dict[fasta_database_str].strip()
    fasta_filename_only = fasta_filename.split("/")[-1]

    # get working_dir
    if working_dir[-1] != '/':
        working_dir = working_dir + '/'

    fasta_filename_dir = working_dir + fasta_filename_only

    # check file exist and open the file
    try:
        with open(fasta_filename) as f: pass
        fasta_file = open(fasta_filename, 'r')
    except:
        try:
            with open(fasta_filename_only) as f: pass
            fasta_file = open(fasta_filename_only, 'r')
        except:
            try:
                with open(fasta_filename_dir) as f: pass
                fasta_file = open(fasta_filename_dir, 'r')
            except:
                print >> sys.stderr, '\nCannot open', fasta_filename
                print >> sys.stderr, 'Check your config file!'
                die("Program exit!")

    # save the fasta ID and description to the fasta_ID_dict
    fasta_ID_dict = {}    # initialize dictionary

    # FASTA ID is space delimited file
    fasta_ID_del = ' '

    # read lines
    for line in fasta_file.readlines():
        line = line.strip()
        # check line start with '>'
        if line.startswith('>'):
            # protein ID is the first word without '>'
            protein_ID = line.split(fasta_ID_del)[0][1:]
            # protein description is the whole line without '>'
            protein_desc = line[1:]
            # replace "#" to "$"
            protein_desc = re.sub('#', '$', protein_desc)
            # replace tab to " "
            protein_desc = re.sub('\t', ' ', protein_desc)

            # save the fasta ID and description to the fasta_ID_dict
            fasta_ID_dict[protein_ID] = protein_desc    # initialize dictionary

    return fasta_ID_dict


## check pep and psm files pair set, and save run#
def get_run_num(pep_file_list, psm_file_list):

    # dictionary of Run# for each pep_file
    run_num_dict = {}
    psm_run_num_dict = {}

    # read multiple pep files
    for pep_file_idx, pep_file in enumerate(pep_file_list):
    
        # check file exist and open
        check_file_exist(pep_file)

        # If base common prefix ends with '.pep.txt', then remove '.pep.txt'
        base_pep_file = pep_file.replace(pep_file_ext, "")

        # make a psm filename using pep filename
        psm_file = base_pep_file + psm_file_ext

        # check psm file exist
        check_file_exist(psm_file)

        # for Run#
        run_tag = 'Run' + str(pep_file_idx + 1)

        # run_num_dict
        run_num_dict[pep_file] = run_tag
        psm_run_num_dict[psm_file] = run_tag

    return (run_num_dict, psm_run_num_dict)

        
# Read and load pep and psm files
def read_run_files(run_num_dict):

    # save the pep file data to the defaultdict
    pep_data_dict = defaultdict(list)
    # save the psm file data to the defaultdict
    psm_data_dict = defaultdict(list)
    # save peptides list (value) to unique protein (key)
    pro_pep_dict = defaultdict(list)
    # save protein list (value) to unique peptide (key)
    pep_pro_dict = defaultdict(list)

    # key = pep_file, val = run_num , sorted by Run# index
    for pep_file, run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):

        # read line with csv
        pep_reader = csv.reader(CommentedFile(open(pep_file, 'rb')),
                                   delimiter='\t')
        # skip header
        headline = pep_reader.next()

        # get data
        for pep_line in pep_reader:
            
            # adapt class PepOutFields
            pep_obj = PepOutFields._make(pep_line)
            identified_peptide = pep_obj.IdentifiedPeptide
            identified_peptide = identified_peptide.strip()
         
            # new unique ID for pep_data_dict is identified_peptide + run_num
            pep_run_id = identified_peptide + "_" + run_num

            pep_data_dict[pep_run_id].append(pep_line)

            # get protein item list
            protein_names = pep_obj.ProteinNames
            pro_item_list = get_item_list(protein_names.strip())

            # for loop of pro_item_list
            for pro_item in pro_item_list:

                pro_item = pro_item.strip()
                
                # save IdentifiedPeptide list to unique protein dict
                if identified_peptide not in pro_pep_dict[pro_item]:
                    pro_pep_dict[pro_item].append(identified_peptide)

                # save protein list to unique peptide dict
                if pro_item not in pep_pro_dict[identified_peptide]:
                    pep_pro_dict[identified_peptide].append(pro_item)

    # try to find indistinguishable set
    indistin_pro_dict = defaultdict(list)
    for pro_key, pep_list in pro_pep_dict.items():
        sorted_pep_list = sorted(set(pep_list))
        sorted_pep_list_join = '_'.join(sorted_pep_list)
        indistin_pro_dict[sorted_pep_list_join].append(pro_key)

    # indistin_key = str(sorted(pep_list)), indistin_value=pro_key list
    for indistin_key, indistin_value in indistin_pro_dict.items():
        # if proteins have a same set of peptides

        if len(indistin_value) > 1:
            # get new protein name
            new_pro_key = list_to_string(sorted(indistin_value))
            new_pro_val = pro_pep_dict[indistin_value[0]]
            # append with new key and value to pro_pep_dict
            # delete provious pep_pro_dict key,value, and append new one
            for new_pro_val_one in new_pro_val:
                pro_pep_dict[new_pro_key].append(new_pro_val_one)
                pep_pro_dict[new_pro_val_one].append(new_pro_key)

                # delete pep_pro_dict indistin protein list of for the peptide
                org_pep_pro_dict_list = pep_pro_dict[new_pro_val_one]
                for indistin_value_one in indistin_value:
                    org_pep_pro_dict_list.remove(indistin_value_one)

                if len(org_pep_pro_dict_list) == 0:
                    del pep_pro_dict[new_pro_val_one]
                else:
                    pep_pro_dict[new_pro_val_one] = org_pep_pro_dict_list

            # delete previous indistinguishable proteins in pro_pep_dict
            for indistin_pro in indistin_value:
                del pro_pep_dict[indistin_pro]

    # key = pep_file, val = run_num , sorted by Run# index
    for pep_file, run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):

        # If base common prefix ends with '.pep.txt', then remove '.pep.txt'
        base_pep_file = pep_file.replace(pep_file_ext, "")

        # make a psm filename using pep filename
        psm_file = base_pep_file + psm_file_ext

        # check psm file exist
        check_file_exist(psm_file)

        # read line with csv
        psm_reader = csv.reader(CommentedFile(open(psm_file, 'rb')),
                                   delimiter='\t')
        # skip header
        headline = psm_reader.next()

        # get data
        for psm_line in psm_reader:

            # check sipros3 or sipros4 input
            psm_obj = None 
            global sipros4_input
            if len(psm_line) == sipros4_psmout_column_length:
                psm_obj = Psm4OutFields._make(psm_line)
                sipros4_input = True
            else:
                psm_obj = PsmOutFields._make(psm_line)
                sipros4_input = False

            # get protein item list
            protein_names = psm_obj.ProteinNames
            pro_item_list = get_item_list(protein_names.strip())

            # save key=protein_id, val=line
            for pro_item_one in pro_item_list:

                # new unique ID for psm_data_dict is protein_name + run_num
                pro_item_id = pro_item_one + "_" + run_num

                # save to dictionary
                psm_data_dict[pro_item_id].append(psm_line)

    return (pep_data_dict, psm_data_dict, pro_pep_dict, pep_pro_dict)


## Greedy algorithm to extract proteins
def greedy_alg(config_dict, pro_pep_dict, pep_pro_dict):

    # get config value
    min_peptide_per_protein = int(config_dict[min_peptide_per_protein_str])
    min_unique_peptide_per_protein = int(config_dict[min_unique_peptide_per_protein_str])

    # return config dictionary

    # First, extract proteins that have >= min_peptide_per_proteins 
    # and at least min_unique_peptide_per_protein of those is unique

    # save the extracted proteins to the list
    pro_greedy_list = []

    # copy pro_pep_dict to pro_pep_dict_red for reduction in greedy steps
    pro_pep_dict_red = defaultdict(list)
    for pro_key, pep_list in pro_pep_dict.items():
        for pep_list_one in pep_list:
            pro_pep_dict_red[pro_key].append(pep_list_one)

    # copy pep_pro_dict to pep_pro_dict_red for reduction in greedy steps
    pep_pro_dict_red = defaultdict(list)
    for pep_key, pro_list in pep_pro_dict.items():
        for pro_list_one in pro_list:
            pep_pro_dict_red[pep_key].append(pro_list_one)

    # for loop of pro_pep_dict
    for pro_key, pep_list in pro_pep_dict.items():

        # if proteins that have >= min_peptide_per_protein
        if len(pep_list) >= min_peptide_per_protein:

            # if at least min_unique_peptide_per_protein is unique
            unique_pep_pro_num = 0
            for pep_list_one in pep_list:
                pep_pro_num = len(pep_pro_dict[pep_list_one])
                if pep_pro_num == 1:
                    unique_pep_pro_num += 1
            # unique peptides num should be >= min_unique_peptide_per_protein
            # if min_unique_peptide_per_protein = 0, then set 1
            if unique_pep_pro_num >= max(min_unique_peptide_per_protein,1):

                # append to the pro_greedy_list
                if pro_key not in pro_greedy_list:
                    pro_greedy_list.append(pro_key)

                # remove the protein
                try:
                    del pro_pep_dict_red[pro_key]
                except:
                    pass
            
                # remove all peptides that are covered by the protein
                for pep_list_item in pep_list:
                    # delete by key
                    try:
                        del pep_pro_dict_red[pep_list_item]
                    except:
                        pass

                    # get new peptide list for the proteins from pep_pro_dict
                    for pro_list_item in pep_pro_dict[pep_list_item]:
                        new_pep_list = pro_pep_dict[pro_list_item]

                        # if new_pep_list is sub list of pep_list, then remove
                        if check_sub_list(new_pep_list, pep_list):
                            try:
                                del pro_pep_dict_red[pro_list_item]
                            except:
                                pass

    # Second, iteratively extract a protein at a time that covers the most peptides
    if len(pro_pep_dict_red.keys()) > 0:
        # Run greedy iterations until it converges
        converge = False

        # Iterate greedy algorithm until it converges
        greedy_step = 0
        while (converge == False):

            greedy_step += 1

            # find a protein that covers the most peptides
            ppdr_idx = 0
            for key_ppdr, val_ppdr in pro_pep_dict_red.items():
                if ppdr_idx == 0:
                    max_key_ppdr = key_ppdr
                    max_len_val_ppdr = len(val_ppdr)
                else:
                    # get current one
                    cur_len_val_ppdr = len(val_ppdr)
                    cur_key_ppdr = key_ppdr
                    # get max one
                    if cur_len_val_ppdr > max_len_val_ppdr:
                        max_len_val_ppdr = cur_len_val_ppdr 
                        max_key_ppdr = cur_key_ppdr
                ppdr_idx += 1

            max_pro_one = max_key_ppdr

            # if proteins that have >= min_peptide_per_protein
            if len(pro_pep_dict_red[max_pro_one]) >= min_peptide_per_protein:

                # if at least min_unique_peptide_per_protein is unique
                unique_pep_pro_num = 0
                for pep_list_one in pro_pep_dict_red[max_pro_one]:
                    pep_pro_num = len(pep_pro_dict[pep_list_one])
                    if pep_pro_num == 1:
                        unique_pep_pro_num += 1
                # if at least min_unique_peptide_per_protein is unique
                if unique_pep_pro_num >= min_unique_peptide_per_protein:

                    # append the protein to the pro_greedy_list
                    pro_greedy_list.append(max_pro_one)

                    # loop for pep set
                    for pep_list_sub_one in pro_pep_dict[max_pro_one]:
                        # delete by key
                        try:
                            del pep_pro_dict_red[pep_list_sub_one]
                        except:
                            pass

                        # loop for pro set 
                        for pro_list_sub_one in pep_pro_dict[pep_list_sub_one]:
                            try:
                                del pro_pep_dict_red[pro_list_sub_one]
                            except:
                                pass

                    # remove the protein
                    try:
                        del pro_pep_dict_red[max_pro_one]
                    except:
                        pass
            
                    # remove all peptides that are covered by the protein
                    for pep_list_item in pro_pep_dict[max_pro_one]:
                        # delete by key
                        try:
                            del pep_pro_dict_red[pep_list_item]
                        except:
                            pass

                        # get new peptide list for the proteins from pep_pro_dict
                        for pro_list_item in pep_pro_dict[pep_list_item]:
                            new_pep_list = pro_pep_dict[pro_list_item]
        
                            # if new_pep_list is sub list of pep_list, then remove
                            if check_sub_list(new_pep_list, pep_list):
                                try:
                                    del pro_pep_dict_red[pro_list_item]
                                except:
                                    pass
                else: 
                    del pro_pep_dict_red[max_pro_one]
            # the max peptides number for protein < min_peptide_per_protein
            else:
                converge = True

            # if there is no protein, then converge
            if len(pro_pep_dict_red.keys()) == 0:
                converge = True

        # greedy algorithm done

    pro_greedy_list = sorted(pro_greedy_list)

    return pro_greedy_list


## Get protein description to handle multiple protein IDs
def get_protein_description(protein_ID, fasta_ID_dict):

    # initialize
    protein_description_list = []

    # if multiple IDs
    if (protein_ID.startswith('{')) and (protein_ID.endswith('}')):

        # multiple proteins exist
        protein_ID_list = get_item_list(protein_ID.strip())
        for protein_ID_one in protein_ID_list:
            
            # check protein ID exist
            if protein_ID_one in fasta_ID_dict:
                protein_description_one = fasta_ID_dict[protein_ID_one]
                protein_description_list.append(protein_description_one)
            else: 
                protein_description_one = "N/A"
                protein_description_list.append(protein_description_one)

    # single ProteinID
    else:
        # check protein ID exist
        if protein_ID in fasta_ID_dict:
            protein_description = fasta_ID_dict[protein_ID]
            protein_description_list.append(protein_description)
        else: 
            protein_description = "N/A"
            protein_description_list.append(protein_description)

    # convert list to string
    protein_description = list_to_string(protein_description_list)

    return protein_description


## check decoy match
def check_decoy_match(ProteinNames, decoy_prefix):

    # match type (correct or decoy) and strip
    match_type = ProteinNames
    match_type_list = []

    if match_type.startswith("{"):
        # if starts with {}, then delete parenthesis { }
        match_type = match_type[1:-1]
        # sometimes multiple proteins
        match_type_list = re.split(r"\s*[,]\s*", match_type.strip())
    else:
        match_type_list.append(match_type)
    # TF -> True or False(decoy match)
    TF = False

    # for loop of proteins
    for match_item in match_type_list:
        # if at least one of matches is True, then match is True
        if not match_item.startswith(decoy_prefix):
            TF = True 
            break

    return TF


## Report output files
def report_output(config_dict,
                  run_num_dict,
                  psm_run_num_dict,
                  pep_data_dict,
                  psm_data_dict,
                  pro_pep_dict,
                  pep_pro_dict,
                  pro_greedy_list,
                  fasta_ID_dict):


    # save .pro.txt data
    pro_out_data = []
    # save .pro2pep.txt data
    pro2pep_out_data = []
    # save .pro2psm.txt data
    pro2psm_out_data = []

    # get config value
    min_peptide_per_protein = int(config_dict[min_peptide_per_protein_str])
    min_unique_peptide_per_protein = int(config_dict[min_unique_peptide_per_protein_str])
    remove_decoy_identification = config_dict[remove_decoy_identification_str]

    # to get decoy_prefix
    decoy_prefix = config_dict[decoy_prefix_str]

    # total number of proteins
    total_proteins_before_filtering = len(pro_pep_dict.keys())
    decoy_proteins_before_filtering = 0
    for key, val in pro_pep_dict.items():
        check_decoy_match_val = check_decoy_match(key, decoy_prefix)
        if check_decoy_match_val is False:
            decoy_proteins_before_filtering += 1
    target_proteins_before_filtering = int(total_proteins_before_filtering) - int(decoy_proteins_before_filtering)

    # total number of identified proteins
    total_proteins_after_filtering = len(pro_greedy_list)
    decoy_proteins_after_filtering = 0
    for pro_one in pro_greedy_list:
        check_decoy_match_val = check_decoy_match(pro_one, decoy_prefix)
        if check_decoy_match_val is False:
            decoy_proteins_after_filtering += 1
    target_proteins_after_filtering = int(total_proteins_after_filtering) - int(decoy_proteins_after_filtering)

    # protein FDR
    protein_fdr = 0.0
    if float(total_proteins_after_filtering) != 0:
        protein_fdr = divide(float(decoy_proteins_after_filtering), float(total_proteins_after_filtering))
    protein_fdr = PrettyFloat(protein_fdr)

    # output header
    def_para_msg = ""
    def_para_msg += "#\t########################################\n"
    def_para_msg += "#\t##### Peptide Assembling by Sipros #####\n"
    def_para_msg += "#\t########################################\n"
    def_para_msg += "#\t\n"
    def_para_msg += "#\t###############\n"
    def_para_msg += "#\t# Input Files #\n"
    def_para_msg += "#\t###############\n"
    def_para_msg += "#\t\n"
    def_para_msg += "#\t[Input_Files]\n"
    def_para_msg += "#\t\n"
    # key = psm_file, val = run_num , sorted by Run# index
    for psm_file, run_num in sorted(psm_run_num_dict.items(), key=lambda x: x[1][-1]):
        def_para_msg += "#\tpsm{" + run_num + "} = " + str(psm_file) + "\n"
    def_para_msg += "#\t\n"
    # key = pep_file, val = run_num , sorted by Run# index
    for pep_file, run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):
        def_para_msg += "#\tpep{" + run_num + "} = " + str(pep_file) + "\n"
    def_para_msg += "#\t\n"
    def_para_msg += "#\t########################\n"
    def_para_msg += "#\t# Filtering Parameters #\n"
    def_para_msg += "#\t########################\n"
    def_para_msg += "#\t\n"
    def_para_msg += "#\t" + pro_iden_str + "\n"
    def_para_msg += "#\t\n"
    def_para_msg += "#\t# Minimum number of peptides per protein\n"
    def_para_msg += "#\t" + min_peptide_per_protein_str + " = " + str(min_peptide_per_protein) + "\n"
    def_para_msg += "#\t\n"
    def_para_msg += "#\t# Minimum number of unique peptides per protein\n"
    def_para_msg += "#\t" + min_unique_peptide_per_protein_str + " = " + str(min_unique_peptide_per_protein) + "\n"
    def_para_msg += "#\t\n"
    def_para_msg += "#\t#######################\n"
    def_para_msg += "#\t# Statistical Results #\n"
    def_para_msg += "#\t#######################\n"
    def_para_msg += "#\t\n"
    def_para_msg += "#\t[Statistical_Results]\n"
    def_para_msg += "#\t\n"
    def_para_msg += "#\t# Numbers of proteins before filtering\n"
    def_para_msg += "#\tDecoy_Proteins_Before_Filtering = " + str(decoy_proteins_before_filtering) + "\n"
    def_para_msg += "#\tTarget_Proteins_Before_Filtering = " + str(target_proteins_before_filtering) + "\n"
    def_para_msg += "#\tTotal_Proteins_Before_Filtering = " + str(total_proteins_before_filtering) + "\n"
    def_para_msg += "#\t\n"
    def_para_msg += "#\t# Numbers of proteins after filtering\n"
    def_para_msg += "#\tDecoy_Proteins_After_Filtering = " + str(decoy_proteins_after_filtering) + "\n"
    def_para_msg += "#\tTarget_Proteins_After_Filtering = " + str(target_proteins_after_filtering) + "\n"
    def_para_msg += "#\tTotal_Proteins_After_Filtering = " + str(total_proteins_after_filtering) + "\n"
    def_para_msg += "#\t\n"
    def_para_msg += "#\t# Protein FDR = Decoy_Proteins_After_Filtering / Total_Proteins_After_Filtering\n"
    def_para_msg += "#\tProtein_FDR = " + set_float_digit(protein_fdr) + "\n"
    def_para_msg += "#\t\n"

    pro_out_file.write(def_para_msg)
    pro2pep_out_file.write(def_para_msg)
    pro2psm_out_file.write(def_para_msg)

    # pro.txt column message
    pro_column_msg = ""
    pro_column_msg += "#\t################\n"
    pro_column_msg += "#\t# Column Names #\n"
    pro_column_msg += "#\t################\n"
    pro_column_msg += "#\t\n"
    pro_column_msg += "#\t[Column_Names]\n"
    pro_column_msg += "#\tProteinID = Names of the protein\n"
    pro_column_msg += "#\tRun#_UniquePeptideCounts = Number of unique peptides in a run\n"
    pro_column_msg += "#\tRun#_TotalPeptideCounts = Number of all peptides in a run\n"
    pro_column_msg += "#\tRun#_UniqueSpectrumCounts = Number of unique PSM in a run\n"
    pro_column_msg += "#\tRun#_TotalSpectrumCounts = Number of all PSM in a run\n"
    pro_column_msg += "#\tRun#_BalancedSpectrumCounts = Balanced spectrum count in a run\n"
    pro_column_msg += "#\tRun#_NormalizedBalancedSpectrumCounts = Normalized Balanced spectrum count in a run\n"
    pro_column_msg += "#\tProteinDescription = Protein description\n"
    pro_column_msg += "#\tTargetMatch = T for target match and F for decoy match\n"
    pro_column_msg += "#\t\n"

    pro_out_file.write(pro_column_msg)

    # pro2pep.txt column message
    pro2pep_column_msg = ""
    pro2pep_column_msg += "#\t################\n"
    pro2pep_column_msg += "#\t# Column Names #\n"
    pro2pep_column_msg += "#\t################\n"
    pro2pep_column_msg += "#\t\n"
    pro2pep_column_msg += "#\t[Column_Names]\n"
    pro2pep_column_msg += "#\t\n"
    pro2pep_column_msg += "#\t+ = Marker of a protein line\n"
    pro2pep_column_msg += "#\tProteinID = Names of the protein\n"
    pro2pep_column_msg += "#\tRun#_UniquePeptideCounts = Number of unique peptides in a run\n"
    pro2pep_column_msg += "#\tRun#_TotalPeptideCounts = Number of all peptides in a run\n"
    pro2pep_column_msg += "#\tRun#_UniqueSpectrumCounts = Number of unique PSM in a run\n"
    pro2pep_column_msg += "#\tRun#_TotalSpectrumCounts = Number of all PSM in a run\n"
    pro2pep_column_msg += "#\tRun#_BalancedSpectrumCounts = Balanced spectrum count in a run\n"
    pro2pep_column_msg += "#\tRun#_NormalizedBalancedSpectrumCounts = Normalized Balanced spectrum count in a run\n"
    pro2pep_column_msg += "#\tProteinDescription = Protein description\n"
    pro2pep_column_msg += "#\tTargetMatch = T for target match and F for decoy match\n"
    pro2pep_column_msg += "#\t\n"
    pro2pep_column_msg += "#\t* = Marker of a peptide line\n"
    pro2pep_column_msg += "#\tIdentifiedPeptide = Identified peptide sequence with potential PTMs and mutations\n"
    pro2pep_column_msg += "#\tParentCharge = Charge state of identified peptide\n"
    pro2pep_column_msg += "#\tOriginalPeptide = Original peptide sequence in the FASTA file\n"
    pro2pep_column_msg += "#\tProteinNames = Names of proteins of the peptide\n"
    pro2pep_column_msg += "#\tProteinCount = Number of proteins that the peptide can be assigned to\n"
    pro2pep_column_msg += "#\tTargetMatch = T for target match and F for decoy match\n"
    pro2pep_column_msg += "#\tSpectralCount = Number of PSMs in which the peptide is identified\n"
    pro2pep_column_msg += "#\tBestScore = The best score of those PSMs\n"
    pro2pep_column_msg += "#\tPSMs = List of PSMs for the peptide: FT2_Filename[Scan_Number]\n"
    pro2pep_column_msg += "#\tScanType = Scan type of those PSMs\n"
    pro2pep_column_msg += "#\tSearchName = Sipros search name\n"
    pro2pep_column_msg += "#\t\n"

    pro2pep_out_file.write(pro2pep_column_msg)

    # pro2psm.txt column message
    pro2psm_column_msg = ""
    pro2psm_column_msg += "#\t################\n"
    pro2psm_column_msg += "#\t# Column Names #\n"
    pro2psm_column_msg += "#\t################\n"
    pro2psm_column_msg += "#\t\n"
    pro2psm_column_msg += "#\t[Column_Names]\n"
    pro2psm_column_msg += "#\t\n"
    pro2psm_column_msg += "#\t+ = Marker of a protein line\n"
    pro2psm_column_msg += "#\tProteinID = Names of the protein\n"
    pro2psm_column_msg += "#\tRun#_UniquePeptideCounts = Number of unique peptides in a run\n"
    pro2psm_column_msg += "#\tRun#_TotalPeptideCounts = Number of all peptides in a run\n"
    pro2psm_column_msg += "#\tRun#_UniqueSpectrumCounts = Number of unique PSM in a run\n"
    pro2psm_column_msg += "#\tRun#_TotalSpectrumCounts = Number of all PSM in a run\n"
    pro2psm_column_msg += "#\tRun#_BalancedSpectrumCounts = Balanced spectrum count in a run\n"
    pro2pep_column_msg += "#\tRun#_NormalizedBalancedSpectrumCounts = Normalized Balanced spectrum count in a run\n"
    pro2psm_column_msg += "#\tProteinDescription = Protein description\n"
    pro2psm_column_msg += "#\tTargetMatch = T for target match and F for decoy match\n"
    pro2psm_column_msg += "#\t\n"
    pro2psm_column_msg += "#\t* = Marker of a PSM line\n"
    pro2psm_column_msg += "#\tFilename = Filename of input FT2 file\n"
    pro2psm_column_msg += "#\tScanNumber = Scan number of the PSM\n"
    pro2psm_column_msg += "#\tParentCharge = Charge state of the PSM\n"
    pro2psm_column_msg += "#\tMeasuredParentMass = Measured parent mass\n"
    pro2psm_column_msg += "#\tCalculatedParentMass = Calculated parent mass from peptide sequence\n"
    pro2psm_column_msg += "#\tMassErrorDa = Mass error in Da with 1-Da error correction\n"
    pro2psm_column_msg += "#\tMassErrorPPM = Mass error in PPM with 1-Da error correction\n"
    pro2psm_column_msg += "#\tScanType = Scan type of the PSM\n"
    pro2psm_column_msg += "#\tSearchName = Sipros search name\n"
    pro2psm_column_msg += "#\tScoringFunction = Scoring function used in the search\n"
    pro2psm_column_msg += "#\tScore = Score\n"
    pro2psm_column_msg += "#\tDeltaZ = Difference between the best PSM and the next best PSM of this scan\n"
    pro2psm_column_msg += "#\tDeltaP = Difference between the best modified PSM and its PTM isoform\n"
    pro2psm_column_msg += "#\tIdentifiedPeptide = Identified peptide sequence with potential PTMs and mutations\n"
    pro2psm_column_msg += "#\tOriginalPeptide = Original peptide sequence in the FASTA file\n"
    pro2psm_column_msg += "#\tProteinNames = Names of proteins of the peptide\n"
    pro2psm_column_msg += "#\tProteinCount = Number of proteins that the peptide can be assigned to\n"
    pro2psm_column_msg += "#\tTargetMatch = T for target match and F for decoy match\n"
    if sipros4_input:
        pro2psm_column_msg += "#\tAveAtom% = Average of atom% estimate\n"
        pro2psm_column_msg += "#\tStdAtom% = Standard deviation of atom% estimate\n"
    pro2psm_column_msg += "#\t\n"

    pro2psm_out_file.write(pro2psm_column_msg)

    # pro_out list
    pro_out_list = ['ProteinID']
    for pep_file, run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):
        pro_out_iter_list = [run_num + '_UniquePeptideCounts',
                             run_num + '_TotalPeptideCounts',
                             run_num + '_UniqueSpectrumCounts',
                             run_num + '_TotalSpectrumCounts',
                             run_num + '_BalancedSpectrumCounts',
                             run_num + '_NormalizedBalancedSpectrumCounts']
        pro_out_list = pro_out_list + pro_out_iter_list
    pro_out_list.append('ProteinDescription')
    pro_out_list.append('TargetMatch')

    pro_out_list_2 = ['+'] + pro_out_list

    # write to pro output
    pro_out_file.write('\t'.join(pro_out_list) + '\n')
    pro2psm_out_file.write('\t'.join(pro_out_list_2) + '\n')
    pro2pep_out_file.write('\t'.join(pro_out_list_2) + '\n')


    # if pro_greedy_list length > 0
    if len(pro_greedy_list) > 0:

        # loop pro_greedy_list
        for protein_one in pro_greedy_list:

            # Column1: Protein ID
            ProteinID = protein_one

            # Column2: Run#_UniquePeptideCounts
            unique_peptide_counts_dict = {}
            # Column3: Run#_TotalPeptideCounts
            total_peptide_counts_dict = {}
            # Column4: Run#_UniqueSpectrumCounts
            unique_spectrum_counts_dict = {}
            # Column5: Run#_TotalSpectrumCounts
            total_spectrum_counts_dict = {}
            # Column6: Run#_BalancedSpectrumCounts
            balanced_spectrum_counts_dict = {}
            # Column7: Run#_NormalizedBalancedSpectrumCounts
            normalized_balanced_spectrum_counts_dict = {}

            # for loop Run#
            for pep_file, run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):

                # get pep_list from pro_pep_dict
                pep_list = pro_pep_dict[protein_one]

                # initialize
                unique_peptide_counts = 0
                total_peptide_counts = 0
                unique_spectrum_counts = 0
                total_spectrum_counts = 0
                balanced_spectrum_counts = 0.0
                share_protein_number = 1
                share_protein_list = []

                # for total_peptide_counts and unique_peptide_counts
                for pep_list_one in pep_list:

                    # get pep_run_id for checking exist of pep_data_dict
                    pep_run_id = pep_list_one + "_" + run_num

                    # for total_peptide_count
                    if pep_run_id in pep_data_dict:
                        total_peptide_counts += len(pep_data_dict[pep_run_id])

                        # for total_spectrum_count and balanced spectrum count
                        for pep_data_dict_one in pep_data_dict[pep_run_id]:
                            pep_data_dict_one_obj = PepOutFields._make(pep_data_dict_one)
                            pep_spectral_count = int(pep_data_dict_one_obj.SpectralCount)
                            total_spectrum_counts += pep_spectral_count

                            # initialize for balanced spectrum count
                            each_balanced_pro_num = 0 
                            each_balanced_spectrum_counts = 0.0 

                            # count each_balanced_pro_num for each peptide
                            for pep_pro_dict_item in pep_pro_dict[pep_list_one]:
                                # if protein is included int he protein greedy list, then count it
                                if pep_pro_dict_item in pro_greedy_list:
                                    each_balanced_pro_num += 1

                            # get each_balanced_spectrum_counts by division
                            each_balanced_spectrum_counts = divide(float(pep_spectral_count), float(each_balanced_pro_num))

                            # add to balanced_spectrum_counts
                            balanced_spectrum_counts += each_balanced_spectrum_counts

                        # if unique peptide
                        if len(pep_pro_dict[pep_list_one]) == 1:
                            unique_peptide_counts += len(pep_data_dict[pep_run_id])

                            # for uniq_spectrum_count
                            for pep_data_dict_one in pep_data_dict[pep_run_id]:
                                pep_data_dict_one_obj = PepOutFields._make(pep_data_dict_one)
                                pep_spectral_count = int(pep_data_dict_one_obj.SpectralCount)
                                unique_spectrum_counts += pep_spectral_count
                            
                # save for each run_num
                unique_peptide_counts_dict[run_num] = unique_peptide_counts
                total_peptide_counts_dict[run_num] = total_peptide_counts
                unique_spectrum_counts_dict[run_num] = unique_spectrum_counts
                total_spectrum_counts_dict[run_num] = total_spectrum_counts
                balanced_spectrum_counts_dict[run_num] = balanced_spectrum_counts
                normalized_balanced_spectrum_counts_dict[run_num] = balanced_spectrum_counts

            # Column : ProteinDescription
            ProteinDescription = get_protein_description(ProteinID, fasta_ID_dict)

            # Column : TargetMatch
            target_match_val = check_decoy_match(ProteinID, decoy_prefix)
            TargetMatch = 'T' if target_match_val is True else 'F'

            # pro_out_list
            pro_out_list= []
            pro_out_list.append(ProteinID)
            first_bsc_idx = 0    # the first column index of BalancedSpectrumCounts
            first_nbsc_idx = 0    # the first column index of NormalizedBalancedSpectrumCounts
            inc_int = 0    # the increment integer for each run
            for pep_file, run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):
                run_inc_index = 0    # for column[0]:ProteinID
                pro_out_list.append(str(unique_peptide_counts_dict[run_num]))
                run_inc_index += 1
                pro_out_list.append(str(total_peptide_counts_dict[run_num]))
                run_inc_index += 1
                pro_out_list.append(str(unique_spectrum_counts_dict[run_num]))
                run_inc_index += 1
                pro_out_list.append(str(total_spectrum_counts_dict[run_num]))
                run_inc_index += 1
                pro_out_list.append(str(balanced_spectrum_counts_dict[run_num]))
                run_inc_index += 1
                first_bsc_idx = run_inc_index    # the first column index of BalancedSpectrumCounts
                pro_out_list.append(str(normalized_balanced_spectrum_counts_dict[run_num]))
                run_inc_index += 1
                first_nbsc_idx = run_inc_index    # the first column index of NormalizedBalancedSpectrumCounts
                inc_int = run_inc_index    # the increment integer for each run
            pro_out_list.append(ProteinDescription)
            pro_out_list.append(TargetMatch)

            #pro_out_data.append(pro_out_list)
            if (remove_decoy_identification in ["No", "NO","no"]) or (TargetMatch == 'T'):
                pro_out_data.append(pro_out_list)

        # for calculating NormalizedBalancedSpectrumCounts
        sum_bsc_list = []    # save sum_bsc for each run into list

        # for loop each run to calculate avg_sum_bsc
        for pep_file, run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):
            run_int = int(run_num[3:])    # run integer (Run1 -> 1)
            bsc_idx = first_bsc_idx + inc_int*(run_int-1)    # calculate BalancedSpectrumCounts index for each run
            sum_bsc = sum(float(item[bsc_idx]) for item in pro_out_data)    # sum of BSC
            sum_bsc_list.append(sum_bsc)    # append to list
        avg_sum_bsc = divide(sum(sum_bsc_list),len(sum_bsc_list))
        
        # for loop of pro_out_data
        for pro_out_data_idx, pro_out_data_item in enumerate(pro_out_data):
            # for loop each run
            for pep_file, run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):
                run_int = int(run_num[3:])    # run integer (Run1 -> 1)
                bsc_idx = first_bsc_idx + inc_int*(run_int-1)    # calculate BalancedSpectrumCounts column index
                nbsc_idx = first_nbsc_idx + inc_int*(run_int-1)    # calculate NormalizedBalancedSpectrumCounts column index
                nbsc_factor = divide(avg_sum_bsc, sum_bsc_list[(run_int - 1)])

                bsc_val = float(pro_out_data_item[bsc_idx])
                nbsc_val = bsc_val * nbsc_factor
                nbsc_val = round(nbsc_val, 2)

                # update pro_out_data_item
                pro_out_data_item[nbsc_idx] = str(nbsc_val)

            # update pro_out_data and write to file
            pro_out_data[pro_out_data_idx] = pro_out_data_item
            pro_out_file.write('\t'.join(pro_out_data_item) + '\n') 

    # for pep out
    pep_out_list = ['*',
                    'IdentifiedPeptide',    #0  
                    'ParentCharge',         #1  
                    'OriginalPeptide',      #2  
                    'ProteinNames',         #3  
                    'ProteinCount',         #4  
                    'TargetMatch',          #5  
                    'SpectralCount',        #6
                    'BestScore',            #7
                    'PSMs',                 #8
                    'ScanType',             #9  
                    'SearchName']           #10  
    pro2pep_out_file.write('\t'.join(pep_out_list) + '\n')

    # if pro_greedy_list length > 0
    if len(pro_greedy_list) > 0:

        # for pro2pep_out_file
        for pro_out_data_item in pro_out_data:
            pro2pep_out_file.write('+\t' + '\t'.join(pro_out_data_item) + '\n')
            pro_ID = pro_out_data_item[0]
            pep_list = pro_pep_dict[pro_ID]

            # for loop pep_list
            for pep_list_item in pep_list:                                                                          

                # for loop Run#
                for pep_file, run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):

                    pep_run_id = pep_list_item + "_" + run_num
                    pep_data_one = pep_data_dict[pep_run_id]                                                         

                    for pep_data_one_each in pep_data_one:                                                              
                        pro2pep_out_file.write('*\t' + '\t'.join(pep_data_one_each) + '\n') 

    # for psm out
    psm_out_list = ['*',
                    'Filename',             #0  
                    'ScanNumber',           #1  
                    'ParentCharge',         #2  
                    'MeasuredParentMass',   #3  
                    'CalculatedParentMass', #4
                    'MassErrorDa',          #5  
                    'MassErrorPPM',         #6  
                    'ScanType',             #7  
                    'SearchName',           #8  
                    'ScoringFunction',      #9  
                    'Score',                #10 
                    'DeltaZ',               #11 
                    'DeltaP',               #12 
                    'IdentifiedPeptide',    #13 
                    'OriginalPeptide',      #14 
                    'ProteinNames',         #15 
                    'ProteinCount',         #16 
                    'TargetMatch']          #17 
    if sipros4_input:
        psm_out_list.append('AveAtom%')     #18
        psm_out_list.append('StdAtom%')     #19

    pro2psm_out_file.write('\t'.join(psm_out_list) + '\n')

    # if pro_greedy_list length > 0
    if len(pro_greedy_list) > 0:

        # for pro2psm_out_file
        for pro_out_data_item in pro_out_data:
            pro2psm_out_file.write('+\t' + '\t'.join(pro_out_data_item) + '\n')
            pro_ID = pro_out_data_item[0]

            # to handle indistinguishable proteins
            if pro_ID.startswith('{'):
                pro_ID = pro_ID[1:-1]
                pro_ID_list = re.split(r"\s*[,]\s*", pro_ID.strip())

                # all target psm list
                all_psm_list = []

                for pro_ID_one in pro_ID_list:

                    # for loop Run#
                    for pep_file, run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):

                        psm_run_id = pro_ID_one + "_" + run_num
                        psm_list = psm_data_dict[psm_run_id]

                        for psm_list_item in psm_list:
                            all_psm_list.append(psm_list_item)

                unique_psm_list = [list(x) for x in set(tuple(x) for x in all_psm_list)]

                # print 
                for psm_list_item in unique_psm_list:
                    pro2psm_out_file.write('*\t' + '\t'.join(psm_list_item) + '\n')
                        
                        
            else:
                # for loop Run#
                for pep_file, run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):

                    psm_run_id = pro_ID + "_" + run_num
                    psm_list = psm_data_dict[psm_run_id]

                    # print 
                    for psm_list_item in psm_list:
                        pro2psm_out_file.write('*\t' + '\t'.join(psm_list_item) + '\n')


## +------+
## | Main |
## +------+
def main(argv=None):

    # try to get arguments and error handling
    try:
        if argv is None:
            argv = sys.argv
        try:
            # parse options
            (working_dir, config_filename) = parse_options(argv)

            # Display work start and time record
            start_time = datetime.now()
            sys.stderr.write('[%s] Beginning sipros_peptides_assembling.py run (V%s)\n' % (curr_time(), get_version()))
            sys.stderr.write('------------------------------------------------------------------------------\n')

            # Parse options and get config file
            sys.stderr.write('[Step 1] Parse options and read fasta file:                   Running -> ')
            # Call parse_config to open and read config file
            config_dict = parse_config(config_filename)
            # Read fasta file and retrieve protein ID and description
            fasta_ID_dict = read_fasta_file(working_dir, config_dict)
            # Get .pep.txt output file(s) in working directory
            pep_file_list = get_file_list_with_ext(working_dir, pep_file_ext)
            # Get .psm.txt output file(s) in working directory
            psm_file_list = get_file_list_with_ext(working_dir, psm_file_ext)
            # check pep and psm files pair set, and save run#
            (run_num_dict, psm_run_num_dict) = get_run_num(pep_file_list, psm_file_list)
            # Get base_out for output
            base_out_default = 'Sipros_searches'
            base_out = get_base_out(pep_file_list, base_out_default, working_dir)
            sys.stderr.write('Done!\n')

            # Read and load pep and psm files
            sys.stderr.write('[Step 2] Load %s file(s):                               Running -> ' % (pep_file_ext))
            (pep_data_dict, psm_data_dict, pro_pep_dict, pep_pro_dict) = read_run_files(run_num_dict)
            sys.stderr.write('Done!\n')

            # Merge indistinguishable proteins that have an identical set of peptides
            sys.stderr.write('[Step 3] Merge indistinguishable proteins:                    Running -> Done!\n')

            # extract proteins that have >2 peptides and at least one of those is unique
            # then iteratively extract a protein at a time that covers the most peptides
            sys.stderr.write('[Step 4] Greedy algorithm for identifying a list of proteins: Running -> ')
            (pro_greedy_list) = greedy_alg(config_dict, pro_pep_dict, pep_pro_dict)
            sys.stderr.write('Done!\n')

            # Report output
            sys.stderr.write('[Step 5] Report output:                                       Running -> ')
            # Open output files
            global pro_out_file
            pro_out_file = open(base_out + ".pro.txt", 'wb')
            global pro2pep_out_file
            pro2pep_out_file = open(base_out + ".pro2pep.txt", 'wb')
            global pro2psm_out_file
            pro2psm_out_file = open(base_out + ".pro2psm.txt", 'wb')
            # Report output files
            report_output(config_dict,
                          run_num_dict,
                          psm_run_num_dict,
                          pep_data_dict,
                          psm_data_dict,
                          pro_pep_dict,
                          pep_pro_dict,
                          pro_greedy_list,
                          fasta_ID_dict)
            sys.stderr.write('Done!\n')

            # Time record, calculate elapsed time, and display work end
            finish_time = datetime.now()
            duration = finish_time - start_time
            sys.stderr.write('------------------------------------------------------------------------------\n')
            sys.stderr.write('[%s] Ending sipros_peptides_assembling.py run\n' % curr_time())
            sys.stderr.write('Run complete [%s elapsed]\n' %  format_time(duration))

        # Error handling
        except Usage, err:
            print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
            return 2

    # Error handling
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "for help use -h/--help"
        return 2


## If this program runs as standalone, then exit.
if __name__ == "__main__":
    sys.exit(main())

