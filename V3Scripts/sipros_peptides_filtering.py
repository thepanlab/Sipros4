#!/usr/bin/python

"""
sipros_peptides_filtering.py

sipros_peptides_filtering.py is the first post-processing program
after running of sipros to find statistically significant peptides 
and report them.

Created by Tae-Hyuk (Ted) Ahn on 09/23/2012.
Copyright (c) 2012 Tae-Hyuk Ahn (ORNL). Allrights reserved.
"""

## Import Python package modules
import sys, getopt, warnings, os, re
import time
from datetime import datetime
from collections import namedtuple
from collections import defaultdict
import csv
import math

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

## Import classes and definitions in the Sipros module
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

## Find string between two substrings
find_between = sipros_post_module.find_between

## Get file(s) list in working dir with specific file extension
get_file_list_with_ext = sipros_post_module.get_file_list_with_ext

## Get base_out filename
get_base_out = sipros_post_module.get_base_out

## Division error handling
divide = sipros_post_module.divide

## list_to_string
list_to_string = sipros_post_module.list_to_string

## list_to_bracket
list_to_bracket = sipros_post_module.list_to_bracket

## Class for sipros fields object
SiprosFields = sipros_post_module.SiprosFields
Sipros4Fields = sipros_post_module.Sipros4Fields

## Class for PsmOutFields object
PsmOutFields = sipros_post_module.PsmOutFields
Psm4OutFields = sipros_post_module.Psm4OutFields

## Class for PepSubFields object
PepSubFields = sipros_post_module.PepSubFields

## Class for PepDataFields object
PepDataFields = sipros_post_module.PepDataFields

## Class for PepOutFields object
PepOutFields = sipros_post_module.PepOutFields

## list_to_bracket
frange = sipros_post_module.frange

## check_file_exist
check_file_exist = sipros_post_module.check_file_exist

## get_protein_count
get_protein_count = sipros_post_module.get_protein_count

## set_float_digit
set_float_digit = sipros_post_module.set_float_digit

## peptide delete residues
peptide_delete_residues = sipros_post_module.peptide_delete_residues

## merge protein names
merge_protein_names = sipros_post_module.merge_protein_names




## Help message
help_message = '''

Usage:
    python sipros_peptides_filtering.py [options]

Inputs:
    sipros output file(s) (search automatically in current directory)
    sipros config file    (search automatically in current directory)

Options:
    -h/--help
    -v/--version
    -w/--working-dir ./path/    # Directory path containing SIPROS output files
                                # (default = current directory) 
    -c/--config-file SiprosConfig.cfg    # SIPROS config file 
                                         # (default = SiprosConfig.cfg) 

Outputs:
    BaseFilename.psm.txt
    BaseFilename.pep.txt 
        - where BaseFilename = common prefix of inputs
        - if len(BaseFilename) < 5, then BaseFilename = Sipros_searches
'''


## Global variables
rank_set = 1
charge_set = 3
FDR_parameter = 1.0
sipros_file_ext = '.sip'
carbon_mass_diff = 1.003355
pro_iden_str      = '[Protein_Identification]'
decoy_prefix_str  = 'Decoy_Prefix'
FDR_filtering_str = 'FDR_Filtering'
FDR_threshold_str = 'FDR_Threshold'
sipros4_column_length = 15
sipros4_input = None


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
            print "sipros_peptides_filtering.py V%s" % (get_version())
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

    # valiables were defined in global
    # pro_iden_str      = '[Protein_Identification]'
    # decoy_prefix_str  = 'Decoy_Prefix'
    # FDR_filtering_str = 'FDR_Filtering'
    # FDR_threshold_str = 'FDR_Threshold'

    # only save protein_identification config info to config_dict
    for key, value in all_config_dict.items():
        if key == (pro_iden_str + decoy_prefix_str):
            config_dict[decoy_prefix_str] = value
        elif key == (pro_iden_str + FDR_filtering_str):
            config_dict[FDR_filtering_str] = value
        elif key == (pro_iden_str + FDR_threshold_str):
            config_dict[FDR_threshold_str] = value
        else:
            continue

    # return config dictionary
    return config_dict


## check decoy match
def check_decoy_match(ProteinNames, decoy_prefix):

    # match type (correct or decoy) and strip
    match_type = ProteinNames
    # delete parenthesis { }
    match_type = match_type[1:-1]
    # sometimes multiple proteins
    match_type_list = re.split(r"\s*[,]\s*", match_type.strip())
    # TF -> True or False(decoy match)
    TF = False

    # for loop of proteins
    for match_item in match_type_list:
        # if at least one of matches is True, then match is True
        if not match_item.startswith(decoy_prefix):
            TF = True
            break

    return TF


## Read sipros output files
def read_sipros_files(sipros_file_list, config_dict, base_out):

    # save data for searching cut-off threshold
    sipros_psm_data = defaultdict(list)    # initialize dict of list

    # pep_sub_dict for preparing pep_out
    pep_sub_dict = defaultdict(list)    # initialize dict of list

    # decoy_prefix
    decoy_prefix  = str(config_dict[decoy_prefix_str])

    # for psm_data
    psm_data_w = open(psm_data_filename, 'wb')

    # for pep_data
    pep_data_w = open(pep_data_filename, 'wb')

    # read multiple sipros files
    for file_idx, sipros_file in enumerate(sipros_file_list):
        # read line with csv
        sipros_reader = csv.reader(CommentedFile(open(sipros_file, 'rb')),
                                   delimiter='\t')
        # skip header
        headline = sipros_reader.next()

        # get data
        for line_idx, line in enumerate(sipros_reader):

            # check sipros3 or sipros4 input
            sipros_obj = None
            global sipros4_input
            if len(line) == sipros4_column_length:
                sipros_obj = Sipros4Fields._make(line)
                sipros4_input = True
            else:
                sipros_obj = SiprosFields._make(line)
                sipros4_input = False

            # get from max score
            psm_Filename = sipros_obj.Filename.strip()
            psm_ScanNumber = sipros_obj.ScanNumber.strip()
            psm_Score = sipros_obj.Score.strip()
            psm_IdentifiedPeptide  = peptide_delete_residues(sipros_obj.IdentifiedPeptide.strip())
            psm_MeasuredParentMass = sipros_obj.MeasuredParentMass
            psm_CalculatedParentMass  = sipros_obj.CalculatedParentMass
            psm_mass_diff = math.fabs(float(psm_MeasuredParentMass) - float(psm_CalculatedParentMass))
            psm_ProteinNames = sipros_obj.ProteinNames.strip()

            # unique_ID for hash
            psm_ID = psm_Filename.split("/")[-1].strip() + '_+_' + psm_ScanNumber.strip()
            psm_idx = str(file_idx) + '_' + str(line_idx)

            # Case 1:
            # if the psm does not exist in the dataset (Filename, ScanNumber), then save this one
            if psm_ID not in sipros_psm_data:

                # new psm_DeltaZ (actually, sam as the score)
                psm_best_score = float(psm_Score)
                psm_second_best_score = 'X'
                psm_lowest_score = float(psm_Score)

                # save psm_out_list to the sipros_psm_data list
                psm_data_list = [psm_idx,
                                 psm_best_score,
                                 psm_second_best_score,
                                 psm_mass_diff,
                                 psm_IdentifiedPeptide,
                                 psm_ProteinNames,
                                 psm_lowest_score,
                                 psm_CalculatedParentMass]

                # save sipros_psm_data
                sipros_psm_data[psm_ID] = psm_data_list

            # Case2:
            # if already exist, then compare the score
            else:
                prev_psm_idx = sipros_psm_data[psm_ID][0]
                prev_psm_best_score = sipros_psm_data[psm_ID][1]
                prev_psm_second_best_score = sipros_psm_data[psm_ID][2]
                prev_psm_mass_diff = sipros_psm_data[psm_ID][3]
                prev_psm_IdentifiedPeptide = sipros_psm_data[psm_ID][4]
                prev_psm_ProteinNames = sipros_psm_data[psm_ID][5]
                prev_psm_lowest_score = sipros_psm_data[psm_ID][6]
                prev_psm_CalculatedParentMass = sipros_psm_data[psm_ID][7]

                # Case 2.1: if new score > previous best score, then replace
                if float(psm_Score) > float(prev_psm_best_score):
                    # second best score shoud be previous best score
                    psm_best_score = float(psm_Score)
                    psm_second_best_score = float(prev_psm_best_score)
                    # replace
                    psm_data_list = [psm_idx,
                                     psm_best_score,
                                     psm_second_best_score,
                                     psm_mass_diff,
                                     psm_IdentifiedPeptide,
                                     psm_ProteinNames,
                                     prev_psm_lowest_score,
                                     psm_CalculatedParentMass]

                    # replace sipros_psm_data
                    sipros_psm_data[psm_ID] = psm_data_list

                # Case 2.2: if new score < previous score, just update the second best score
                elif float(psm_Score) < float(prev_psm_best_score):

                    # best score should be kept
                    # second best score shoud be considered
                    # Case2.2.1: prev_psm_second_best_score doesn not exist (X), then psm_second_best_score = psmScore
                    if str(prev_psm_second_best_score) == 'X':
                        psm_second_best_score = float(psm_Score)

                    # Case2.2.2: prev_psm_second_best_score != "X", then compare the values
                    else:
                        # Case 2.2.2.1: psm_Score > prev_psm_second_best_score, then update the second best score to the current score
                        if float(psm_Score) > float(prev_psm_second_best_score):
                            psm_second_best_score = float(psm_Score)
                        # Case 2.2.2.2: psm_Score <= prev_psm_best_score, then do NOT update the second best score
                        else:
                            psm_second_best_score = float(prev_psm_second_best_score)

                    # only update psm_second_best_score
                    sipros_psm_data[psm_ID][2] = psm_second_best_score

                    # update lowest score
                    if float(psm_Score) < float(prev_psm_lowest_score):
                        sipros_psm_data[psm_ID][6] = float(psm_Score)

                # Case 2.3: if new score == previous score, 
                # find min of difference between MeasuredParentMass and CalculatedParentMass
                else:

                    # Case 2.3.1:
                    # if psm mass diff < previous mass diff, then update
                    if float(psm_mass_diff) < float(prev_psm_mass_diff) and abs(float(psm_mass_diff) - float(prev_psm_mass_diff)) > 0.00005:

                        # second best score shoud be previous best score
                        psm_best_score = float(psm_Score)
                        psm_second_best_score = float(prev_psm_best_score)
                        # replace
                        psm_data_list = [psm_idx,
                                         psm_best_score,
                                         psm_second_best_score,
                                         psm_mass_diff,
                                         psm_IdentifiedPeptide,
                                         psm_ProteinNames,
                                         prev_psm_lowest_score,
                                         psm_CalculatedParentMass]
                        sipros_psm_data[psm_ID] = psm_data_list

                    # Case 2.3.2:
                    # if psm mass diff > previous mass diff, then update the second best score only
                    elif float(psm_mass_diff) > float(prev_psm_mass_diff) and abs(float(psm_mass_diff) - float(prev_psm_mass_diff)) > 0.00005:

                        # only update psm_second_best_score as psm_Score
                        psm_second_best_score = float(psm_Score)
                        sipros_psm_data[psm_ID][2] = psm_second_best_score

                    # Case 2.3.3:
                    # if prev mass diff == current mass diff, then compare PTM score
                    else:

                        # calcualte PTM scores
                        s1 = ''.join([char if char.isalnum() else '$' for char in psm_IdentifiedPeptide ])
                        s2 = ''.join([char if char.isalnum() else '$' for char in prev_psm_IdentifiedPeptide ])
                        s1_special_char_count = s1.count('$') - 2
                        s2_special_char_count = s2.count('$') - 2

                        # Case 2.3.3.1
                        # if new PTM score is < previous PTM score, then update
                        if s1_special_char_count < s2_special_char_count:
                  
                            # update 
                            psm_best_score = float(psm_Score)
                            psm_second_best_score = float(prev_psm_best_score)
                            # replace
                            psm_data_list = [psm_idx,
                                             psm_best_score,
                                             psm_second_best_score,
                                             psm_mass_diff,
                                             psm_IdentifiedPeptide,
                                             psm_ProteinNames,
                                             prev_psm_lowest_score,
                                             psm_CalculatedParentMass]
                            sipros_psm_data[psm_ID] = psm_data_list

                        # Case 2.3.3.2
                        # if new PTM score is > previous PTM score, then update the second best score only
                        elif s1_special_char_count > s2_special_char_count:
                  
                            # only update psm_second_best_score as psm_Score
                            psm_second_best_score = float(psm_Score)
                            sipros_psm_data[psm_ID][2] = psm_second_best_score

                        # Case 2.3.3.3:
                        # if new PTM score -= previous PTM score, then consider peptide alphabetical order
                        else:

                            # Case 2.3.3.3.1:
                            # if new peptide alphabetical orderis less than previous one, then update
                            if psm_IdentifiedPeptide.upper() < prev_psm_IdentifiedPeptide.upper():

                                # update
                                psm_best_score = float(psm_Score)
                                psm_second_best_score = float(prev_psm_best_score)
                                # replace
                                psm_data_list = [psm_idx,
                                                 psm_best_score,
                                                 psm_second_best_score,
                                                 psm_mass_diff,
                                                 psm_IdentifiedPeptide,
                                                 psm_ProteinNames,
                                                 prev_psm_lowest_score,
                                                 psm_CalculatedParentMass]
                                sipros_psm_data[psm_ID] = psm_data_list

                            # Case 2.3.3.3.2:
                            # if new peptide alphabetical orderis > previous one, then update the second best score only
                            elif psm_IdentifiedPeptide.upper() > prev_psm_IdentifiedPeptide.upper():

                                # only update psm_second_best_score as psm_Score
                                psm_second_best_score = float(psm_Score)
                                sipros_psm_data[psm_ID][2] = psm_second_best_score

                            # Case 2.3.3.3.3:
                            # if IdentifiedPeptide == Previous IdentifiedPeptide, then merge proteins
                            # but do not update the second best score
                            else:
                                merge_ProteinNames = merge_protein_names(psm_ProteinNames, prev_psm_ProteinNames)
                                sipros_psm_data[psm_ID][5] = merge_ProteinNames

    # finally, if the psm second best score is "X", then save it to the best score
    for psm_ID, psm_data_list in sipros_psm_data.iteritems():
        psm_best_score = str(psm_data_list[1])
        psm_second_best_score = str(psm_data_list[2])
        if psm_second_best_score == 'X':
            sipros_psm_data[psm_ID][2] = float(psm_best_score)
            

    # read multiple sipros files
    for file_idx, sipros_file in enumerate(sipros_file_list):
        # read line with csv
        sipros_reader = csv.reader(CommentedFile(open(sipros_file, 'rb')),
                                   delimiter='\t')
        # skip header
        headline = sipros_reader.next()

        # get data
        for line_idx, line in enumerate(sipros_reader):
            
            # check sipros3 or sipros4 input
            sipros_obj = None
            if sipros4_input:
                sipros_obj = Sipros4Fields._make(line)
            else:
                sipros_obj = SiprosFields._make(line)


            # get from max score
            psm_Filename = sipros_obj.Filename.strip()
            psm_ScanNumber = sipros_obj.ScanNumber.strip()
            psm_Score = sipros_obj.Score.strip()
            psm_IdentifiedPeptide  = peptide_delete_residues(sipros_obj.IdentifiedPeptide.strip())
            psm_MeasuredParentMass = sipros_obj.MeasuredParentMass
            psm_CalculatedParentMass  = sipros_obj.CalculatedParentMass
            psm_mass_diff = math.fabs(float(psm_MeasuredParentMass) - float(psm_CalculatedParentMass))
            psm_ProteinNames = sipros_obj.ProteinNames.strip()

            # unique_ID for hash
            psm_ID = psm_Filename.split("/")[-1].strip() + '_+_' + psm_ScanNumber.strip()
            psm_idx = str(file_idx) + '_' + str(line_idx)

            # get data from spros_psm_data
            psm_data_list = sipros_psm_data[psm_ID]

            #psm_data_list = [psm_idx,                  # psm_data_list[0]
            #                 psm_best_score,           # psm_data_list[1]
            #                 psm_second_best_score,    # psm_data_list[2]
            #                 psm_mass_diff,            # psm_data_list[3]
            #                 psm_IdentifiedPeptide,    # psm_data_list[4]
            #                 psm_ProteinNames,         # psm_data_list[5]
            #                 prev_psm_lowest_score,    # psm_data_list[6]
            #                 psm_CalculatedParentMass] # psm_data_list[7]

            # calcualte PTM scores
            s1 = ''.join([char if char.isalnum() else '$' for char in psm_IdentifiedPeptide ]) # change all special characters to $
            s2 = ''.join([char if char.isalnum() else '$' for char in psm_data_list[4]])
            s1_alnum = ''.join(char for char in s1 if char.isalnum())   # special characters are removed
            s2_alnum = ''.join(char for char in s2 if char.isalnum())
            s1_special_char_count = s1.count('$')-2                     # subtract 2 for [ ]
            s2_special_char_count = s2.count('$')-2 

            # Case 1:
            # if the IdentifiedPeptide doesn't have a PTM, DeltaP = "NA"
            # This will be printed in the next step

            # Case 2:
            # if the highest score IdentifiedPeptide has a PTM
            # if abs(the current psm_CalculatedParentMass - the CalculatedParentMass of the IdentifiedPeptide ) < 0.00005
            # if abs(the current psm_MassDiff - the MassDiff of the IdentifiedPeptide) < 0.00005
            # if the current peptide string without PTM == peptide string without PTM of the IdentifiedPeptide
            # if the number of PTMs of the current peptide == the number of PTMs of the IdentifiedPeptide
            # if the current entire peptide string != the entire string of the IdentifiedPeptide
            # if the current score > lowest score, then update the lowest score to the current score
            if s2_special_char_count > 0 and \
               abs(float(psm_CalculatedParentMass) - float(psm_data_list[7])) < 0.00005 and \
               abs(float(psm_mass_diff) - float(psm_data_list[3])) < 0.00005 and \
               s1_alnum == s2_alnum and \
               s1_special_char_count == s2_special_char_count and \
               psm_IdentifiedPeptide != psm_data_list[4] and \
               float(psm_Score) > psm_data_list[6]:

                sipros_psm_data[psm_ID][6] = float(psm_Score)


    # read multiple sipros files
    for file_idx, sipros_file in enumerate(sipros_file_list):
        # read line with csv
        sipros_reader = csv.reader(CommentedFile(open(sipros_file, 'rb')),
                                   delimiter='\t')
        # skip header
        headline = sipros_reader.next()

        # get data
        for line_idx, line in enumerate(sipros_reader):
            
            # check sipros3 or sipros4 input
            sipros_obj = None
            if sipros4_input:
                sipros_obj = Sipros4Fields._make(line)
            else:
                sipros_obj = SiprosFields._make(line)


            # get from max score
            psm_Filename   = sipros_obj.Filename.strip()
            psm_ScanNumber = sipros_obj.ScanNumber.strip()
            psm_ParentCharge  = sipros_obj.ParentCharge.strip()
            psm_MeasuredParentMass = sipros_obj.MeasuredParentMass.strip()
            psm_CalculatedParentMass  = sipros_obj.CalculatedParentMass.strip()
            psm_MassErrorDa_org  = float(psm_CalculatedParentMass) - float(psm_MeasuredParentMass)
            psm_MassErrorDa_mod = math.fmod(psm_MassErrorDa_org,carbon_mass_diff)
            if psm_MassErrorDa_mod > 0:
                psm_MassErrorDa_one = math.fabs(psm_MassErrorDa_mod - carbon_mass_diff)
                if psm_MassErrorDa_mod < psm_MassErrorDa_one:
                    psm_MassErrorDa = str(psm_MassErrorDa_mod)
                else:
                    psm_MassErrorDa = str(psm_MassErrorDa_one)
            elif psm_MassErrorDa_mod < 0:
                psm_MassErrorDa_one = math.fabs(psm_MassErrorDa_mod + carbon_mass_diff)
                if math.fabs(psm_MassErrorDa_mod) < psm_MassErrorDa_one:
                    psm_MassErrorDa = str(psm_MassErrorDa_mod)
                else:
                    psm_MassErrorDa = str(psm_MassErrorDa_one)
            else:
                psm_MassErrorDa = str(0.0)
            psm_MassErrorPPM_org = divide(float(psm_MassErrorDa), float(psm_CalculatedParentMass))
            psm_MassErrorPPM = str(psm_MassErrorPPM_org * 1000000)
            psm_ScanType  = sipros_obj.ScanType.strip()
            psm_SearchName  = sipros_obj.SearchName.strip()
            psm_ScoringFunction  = sipros_obj.ScoringFunction.strip()
            psm_Score  = sipros_obj.Score.strip()
            psm_IdentifiedPeptide  = peptide_delete_residues(sipros_obj.IdentifiedPeptide.strip())
            psm_OriginalPeptide = peptide_delete_residues(sipros_obj.OriginalPeptide.strip())

            # unique_ID for hash
            psm_ID = psm_Filename.split("/")[-1].strip() + '_+_' + psm_ScanNumber.strip()
            psm_idx = str(file_idx) + '_' + str(line_idx)

            # get data from spros_psm_data
            psm_data_list = sipros_psm_data[psm_ID]
            psm_idx_in_dict = psm_data_list[0]
            psm_DeltaZ = float(psm_data_list[1]) - float(psm_data_list[2])
            psm_ProteinNames = psm_data_list[5]
            psm_ProteinCount = str(get_protein_count(psm_ProteinNames))
            psm_TF = check_decoy_match(psm_ProteinNames, decoy_prefix)
            psm_TargetMatch = 'T' if psm_TF is True else 'F' 
            # calcualte PTM scores
            psm_PTM_diff = float(psm_data_list[1]) - float(psm_data_list[6])
            s2 = ''.join([char if char.isalnum() else '$' for char in psm_data_list[4]])    # for the highest score IdentifiedPeptide PTM
            s2_special_char_count = s2.count('$') - 2 
            psm_DeltaP = 'NA'
            if s2_special_char_count != 0:
                psm_DeltaP = str(psm_PTM_diff)

            # for sipros4
            psm_AveAtom = sipros_obj.AveAtom.strip() if sipros4_input else 'N'
            psm_StdAtom = sipros_obj.StdAtom.strip() if sipros4_input else 'N'

            # check psm_idx, if exist, then save to file
            if psm_idx == psm_idx_in_dict:
                
                # save psm_out_list to the sipros_psm_data list
                psm_out_list = [str(psm_Filename),
                                str(psm_ScanNumber),
                                str(psm_ParentCharge),
                                str(psm_MeasuredParentMass),
                                str(psm_CalculatedParentMass),
                                str(psm_MassErrorDa),
                                str(psm_MassErrorPPM),
                                str(psm_ScanType),
                                str(psm_SearchName),
                                str(psm_ScoringFunction),
                                str(psm_Score),
                                str(psm_DeltaZ),
                                str(psm_DeltaP),
                                str(psm_IdentifiedPeptide),
                                str(psm_OriginalPeptide),
                                str(psm_ProteinNames),
                                str(psm_ProteinCount),
                                str(psm_TargetMatch)]

                if sipros4_input:
                    psm_out_list.append(str(psm_AveAtom))
                    psm_out_list.append(str(psm_StdAtom))

                psm_data_w.write('\t'.join(psm_out_list) + '\n')

    # release memory and close file writing
    sipros_psm_data = defaultdict(list)    # initialize dict of list
    sipros_psm_data.clear()
    psm_data_w.close()

    # open to read psm_data_file
    psm_data_f = open(psm_data_filename, 'rb')
    psm_data_r = csv.reader(psm_data_f, delimiter='\t')

    # get data
    for line_idx, line in enumerate(psm_data_r):

        # check sipros3 or sipros4 input
        psm_obj = None
        if sipros4_input:
            psm_obj = Psm4OutFields._make(line)
        else:
            psm_obj = PsmOutFields._make(line)

        # pep ID is unique with IdentifiedPeptide and ParentCharge
        pep_ID = psm_obj.IdentifiedPeptide.strip() + '_+_' + psm_obj.ParentCharge.strip()

        # pep_sub_list to save pep info into dictionary
        pep_sub_list = [psm_obj.IdentifiedPeptide,  # IdentifiedPeptide
                        psm_obj.ParentCharge,       # ParentCharge
                        psm_obj.OriginalPeptide,    # OriginalPeptide
                        psm_obj.ProteinNames,       # ProteinNames
                        psm_obj.Score,              # Score
                        psm_obj.Filename,           # Filename
                        psm_obj.ScanNumber,         # ScanNumber
                        psm_obj.ScanType,           # ScanType
                        psm_obj.SearchName]         # SearchName

        # pep_sub_list to save pep info into dictionary
        pep_sub_dict[pep_ID].append(pep_sub_list)

    # close psm_data_f
    psm_data_f.close()

    # prepare sipros_pep_data
    for pep_ID, pep_sub_list in sorted(pep_sub_dict.iteritems(), key=lambda (k,v): (v[0][0], int(v[0][1]))):

        # initialize variable
        pep_ID_best_score = 0.0
        pep_protein_names_list = []
        pep_protein_names_str = ''
        pep_data_list = []
        pep_identified_peptide = pep_ID.strip().split("_+_")[0]
        pep_parent_charge = pep_ID.strip().split("_+_")[1]

        # for loop of pep_sub_list for each pep_ID
        for pep_sub_list_one in pep_sub_list:

            # get PepSubFields object
            pep_sub_obj = PepSubFields._make(pep_sub_list_one)

            # update best score for the pep_ID
            if float(pep_sub_obj.Score) > float(pep_ID_best_score):
                pep_ID_best_score = pep_sub_obj.Score

            # get protein names list
            protein_names_wo_bracket = pep_sub_obj.ProteinNames.strip()[1:-1]
            pep_protein_names_list.append(protein_names_wo_bracket)

        # pep_protein_names_list
        pep_protein_names_list = list(set(pep_protein_names_list))
        pep_protein_names_str = list_to_bracket(pep_protein_names_list)

        # pep_data_list
        pep_data_list = [str(pep_identified_peptide),
                         str(pep_parent_charge),
                         str(pep_ID_best_score),
                         str(pep_protein_names_str)]

        # sipros_pep_data dictionary, key=pep_ID, val=pep_ID_best_score
        pep_data_w.write('\t'.join(pep_data_list) + '\n')

    # close pep_data_w
    pep_data_w.close()

    return (pep_sub_dict)


## get the # of list if the element > threshold
def get_listnum_threshold(input_list, thres_val):
    try:
        new_list = [i for i in input_list if i >= thres_val]
        num_new_list = len(new_list)
        return num_new_list
    except:
        return 0


## FDR calculator
def FDR_calculator(FP, TP):
    FDR_numerator = float(FP)*float(FDR_parameter)
    FDR_denoninator = float(FP) + float(TP)
    FDR_accept = True

    if  FDR_denoninator == 0:
        FDR_value = 0.0
        FDR_accept = False
    else:
        FDR_value = divide(FDR_numerator, FDR_denoninator)
        FDR_accept = True

    return (FDR_accept, float(FDR_value))


## Find cut-off score and other data using given FDR threshold and lists of TP and FT
def get_cutoff_data(FDR_threshold, F_list, T_list):

    # save cutoff score and other data to cutoff_data
    cutoff_data = []

    # get info of lists
    FT_list = F_list + T_list
    max_score = max(FT_list) if len(FT_list) > 0 else 0.0
    min_score = min(FT_list) if len(FT_list) > 0 else 0.0
    avg_score = sum(FT_list) / len(FT_list) if len(FT_list) else 0.0
    step_size = float(avg_score/10000)

    # initialize
    final_cutoff_score = 0.0
    final_FDR_value = 0.0
    final_F_num = 0
    final_T_num = 0
    final_TF_num = 0
    final_accept = False
    prev_TF_num = 0
    FDR_accept = False

    # it max_score=min_score=0.0
    if (max_score == 0.0) and (min_score == 0.0):
        # cutoff_data list
        cutoff_data = [0.0, False]
    else:
        # get list from max to min with decrement stepsize
        cutoff_score_range = frange(max_score, min_score, -step_size)

        # for loop of cutoff_score_range
        for cutoff_score in cutoff_score_range:
            F_num = int(get_listnum_threshold(F_list, cutoff_score))
            T_num = int(get_listnum_threshold(T_list, cutoff_score))
            TF_num = F_num + T_num
            (FDR_accept, FDR_value) = FDR_calculator(F_num, T_num)

            # update final values if conditions satisfies
            # 1) FDR_accept is True
            # 2) FDR_value should be less than or equal to FDR_threshold
            # 3) TF_num is greater than to previous TF_num
            if (FDR_accept is True) and (FDR_value <= FDR_threshold) and (TF_num > prev_TF_num) :
                final_cutoff_score = cutoff_score
                final_FDR_value = FDR_value
                final_F_num = F_num
                final_T_num = T_num
                final_TF_num = final_T_num + final_F_num
                final_accept = FDR_accept

            # previous TF_num
            prev_TF_num = TF_num

        # cutoff_data list
        cutoff_data = [final_cutoff_score, final_accept]

    return cutoff_data


## Calculate cutoff score by given FDR value
def get_cutoff_score(config_dict):

    # cutoff score data dictionary for PSM or Peptides
    cutoff_dict = defaultdict(list)

    # dictionary for TF
    TF_dict = defaultdict(list)

    # get config variables
    decoy_prefix  = str(config_dict[decoy_prefix_str])
    FDR_filtering = str(config_dict[FDR_filtering_str])
    FDR_threshold = float(config_dict[FDR_threshold_str])

    # if FDR_filtering == 'PSM'
    if FDR_filtering == 'PSM':

        # open to read psm_data_file
        psm_data_f = open(psm_data_filename, 'rb')
        psm_data_r = csv.reader(psm_data_f, delimiter='\t')

        # get data
        for line_idx, line in enumerate(psm_data_r):

            data_one_obj = PsmOutFields._make(line)

            # data calss to variables
            protein_names = data_one_obj.ProteinNames
            charge_val = int(data_one_obj.ParentCharge)
            if charge_val >= charge_set:
                charge_val = charge_set
            charge_idx = charge_val - 1
            score = float(data_one_obj.Score)

            # if one of matches is True, then True
            match_TF = check_decoy_match(protein_names, decoy_prefix)
            TF = 1 if match_TF is True else 0

            # [charge_idx][TF] -> (charge_idx*2) + TF -> unique idx
            # this structure is flexible for increasing charge set than fixed list
            # [0][0] -> 0
            # [0][1] -> 1
            # [1][0] -> 2
            # [1][1] -> 3
            # [2][0] -> 4
            # [2][1] -> 5
            TF_dict_idx = (charge_idx*2) + TF
            # save the score to the TF_dict
            TF_dict[TF_dict_idx].append(score)

        # close psm_data_f
        psm_data_f.close()

        # loop for each charge
        for idx in range(charge_set):
            F_idx = (idx*2) + 0
            T_idx = (idx*2) + 1

            # save cutoff score info to cutoff_dict[charge_idx]
            cutoff_dict[idx] = get_cutoff_data(FDR_threshold, TF_dict[F_idx], TF_dict[T_idx])

    # if FDR_filtering == 'Peptide'
    elif FDR_filtering == 'Peptide':

        # open to read psm_data_file
        pep_data_f = open(pep_data_filename, 'rb')
        pep_data_r = csv.reader(pep_data_f, delimiter='\t')

        # get data
        for line_idx, line in enumerate(pep_data_r):

            pep_data_obj = PepDataFields._make(line)

            # data calss to variables
            identified_peptide = pep_data_obj.IdentifiedPeptide
            charge_val = int(pep_data_obj.ParentCharge)
            if charge_val >= charge_set:
                charge_val = charge_set
            charge_idx = charge_val - 1
            protein_names = pep_data_obj.ProteinNames
            score = float(pep_data_obj.BestScore)

            # if one of matches is True, then True
            match_TF = check_decoy_match(protein_names, decoy_prefix)
            TF = 1 if match_TF is True else 0
            # [charge_idx][TF] -> (charge_idx*2) + TF -> unique idx
            TF_dict_idx = (charge_idx*2) + TF
            # save the score to the TF_dict
            TF_dict[TF_dict_idx].append(score)

        # close psm_data_f
        pep_data_f.close()

        # loop for each charge
        for idx in range(charge_set):
            F_idx = (idx*2) + 0
            T_idx = (idx*2) + 1

            # save cutoff score info to cutoff_dict[charge_idx]
            cutoff_dict[idx] = get_cutoff_data(FDR_threshold, TF_dict[F_idx], TF_dict[T_idx])

    return cutoff_dict


## Report output files
def report_output(sipros_file_list, 
                  pep_sub_dict, 
                  cutoff_dict, 
                  config_dict):

    # dictionary for matching T and F(decoy match)
    psm_T_dict = defaultdict(list)
    psm_F_dict = defaultdict(list)
    pep_T_dict = defaultdict(list)
    pep_F_dict = defaultdict(list)

    # filtering parameters
    decoy_prefix  = str(config_dict[decoy_prefix_str])
    FDR_filtering = str(config_dict[FDR_filtering_str])
    FDR_threshold = float(config_dict[FDR_threshold_str])
    
    # header message
    psm_header_msg = ""
    psm_header_msg += "#\t############################################################\n"
    psm_header_msg += "#\t##### PSM (Peptide-Spectrum Match) Filtering by Sipros #####\n"
    psm_header_msg += "#\t############################################################\n"
    psm_header_msg += "#\t\n"

    pep_header_msg = ""
    pep_header_msg += "#\t#######################################\n"
    pep_header_msg += "#\t##### Peptide Filtering by Sipros #####\n"
    pep_header_msg += "#\t#######################################\n"
    pep_header_msg += "#\t\n"

    psm_out_file.write(psm_header_msg)
    pep_out_file.write(pep_header_msg)

    # input file message
    input_file_msg = ""
    input_file_msg += "#\t###############\n"
    input_file_msg += "#\t# Input Files #\n"
    input_file_msg += "#\t###############\n"
    input_file_msg += "#\t\n"
    input_file_msg += "#\t[Input_Files]\n"
    input_file_msg += "#\t\n"
    for idx, file_name in enumerate(sipros_file_list):
        input_file_msg += "#\tsipros_file_" + str(idx + 1) + " = " + file_name + "\n"
    input_file_msg += "#\t\n"

    psm_out_file.write(input_file_msg)
    pep_out_file.write(input_file_msg)

    # protein_identification message
    pro_iden_msg = ""
    pro_iden_msg += "#\t########################\n"
    pro_iden_msg += "#\t# Filtering Parameters #\n"
    pro_iden_msg += "#\t########################\n"
    pro_iden_msg += "#\t\n"
    pro_iden_msg += "#\t[Protein_Identification]\n"
    pro_iden_msg += "#\t\n"
    pro_iden_msg += "#\t# The prefix of decoy sequences' locus IDs in the database\n"
    pro_iden_msg += "#\t" + decoy_prefix_str + " = " + str(decoy_prefix) + "\n"
    pro_iden_msg += "#\t\n"
    pro_iden_msg += "#\t# Level of FDR filtering. Options: \"PSM\" and \"Peptide\"\n"
    pro_iden_msg += "#\t" + FDR_filtering_str + " = " + str(FDR_filtering) + "\n"
    pro_iden_msg += "#\t\n"
    pro_iden_msg += "#\t# FDR threshold for filtering peptide identifications\n"
    pro_iden_msg += "#\t" + FDR_threshold_str + " = " + str(FDR_threshold) + "\n"
    pro_iden_msg += "#\t\n"

    psm_out_file.write(pro_iden_msg)
    pep_out_file.write(pro_iden_msg)

    # cutoff_score_list
    cutoff_score_list = []

    # for total and decoy match
    psm_decoy_match_num = 0
    pep_decoy_match_num = 0

    # open to read psm_data_file
    psm_data_f = open(psm_data_filename, 'rb')
    psm_data_r = csv.reader(psm_data_f, delimiter='\t')

    # get data
    for line_idx, line in enumerate(psm_data_r):

        # check sipros3 or sipros4 input
        data_one_obj = None
        if sipros4_input:
            data_one_obj = Psm4OutFields._make(line)
        else:
            data_one_obj = PsmOutFields._make(line)

        # data calss to variables
        protein_names = data_one_obj.ProteinNames
        charge_val = int(data_one_obj.ParentCharge)
        if charge_val >= charge_set:
            charge_val = charge_set
        charge_idx = charge_val - 1
        score = float(data_one_obj.Score)

        # for decoy_match_num
        psm_match_TF = check_decoy_match(protein_names, decoy_prefix)
        if psm_match_TF is False:
            psm_decoy_match_num += 1

        # get final_cutoff_score and final_accept with charge index
        final_cutoff_score = cutoff_dict[charge_idx][0]
        final_accept       = cutoff_dict[charge_idx][1]

        # if score exist
        if final_accept is True:
            # if score >= cutoff_score
            if score >= final_cutoff_score:
                # if one of matches is True, then True
                match_TF = check_decoy_match(protein_names, decoy_prefix)
                if match_TF is True:
                    psm_T_dict[charge_idx].append(score) 
                else:
                    psm_F_dict[charge_idx].append(score)

    # close psm_data_f
    psm_data_f.close()
    total_psm_before_filtering = line_idx + 1

    # open to read pep_data_file
    pep_data_f = open(pep_data_filename, 'rb')
    pep_data_r = csv.reader(pep_data_f, delimiter='\t')

    # get data
    for line_idx, line in enumerate(pep_data_r):

        # call PsmOutFields class
        pep_data_obj = PepDataFields._make(line)

        # get pep_ID
        pep_ID = pep_data_obj.IdentifiedPeptide.strip() + '_+_' + pep_data_obj.ParentCharge.strip()

        # data calss to variables
        protein_names = pep_data_obj.ProteinNames
        charge_val = int(pep_data_obj.ParentCharge)
        if charge_val >= charge_set:
            charge_val = charge_set
        charge_idx = charge_val - 1
        best_score = float(pep_data_obj.BestScore)

        # for total_match_num and decoy_match_num
        pep_match_TF = check_decoy_match(protein_names, decoy_prefix)
        if pep_match_TF is False:
            pep_decoy_match_num += 1

        # get final_cutoff_score and final_accept with charge index
        final_cutoff_score = cutoff_dict[charge_idx][0]
        final_accept       = cutoff_dict[charge_idx][1]

        # if score exist
        if final_accept is True:

            # if score >= cutoff_score
            if best_score >= final_cutoff_score:

                # if one of matches is True, then True
                match_TF = check_decoy_match(protein_names, decoy_prefix)
                if match_TF is True:
                    pep_T_dict[charge_idx].append(score) 
                else:
                    pep_F_dict[charge_idx].append(score)

    # close psm_data_f
    pep_data_f.close()
    total_peptides_before_filtering = line_idx + 1

    # variables for psm_out
    psm_cutoff_num_list = []
    psm_cutoff_per_list = []
    psm_cutoff_per = 0.0
    psm_cutoff_F_sum = 0
    psm_cutoff_TF_sum = 0
    psm_cutoff_num_sum = ''

    # PSM decoy_match_num
    decoy_psm_before_filtering = psm_decoy_match_num
    target_psm_before_filtering = total_psm_before_filtering - decoy_psm_before_filtering

    # peptide decoy_match_num
    decoy_peptides_before_filtering = pep_decoy_match_num
    target_peptides_before_filtering = total_peptides_before_filtering - decoy_peptides_before_filtering

    # for PSM statistics
    for idx in range(charge_set):
        # cutoff_dict
        psm_final_cutoff_score = cutoff_dict[idx][0]
        psm_final_accept       = cutoff_dict[idx][1]
        psm_final_F_num        = len(psm_F_dict[idx])
        psm_final_T_num        = len(psm_T_dict[idx])
        psm_final_TF_num       = psm_final_F_num + psm_final_T_num
 

        # if final_accept is True
        if psm_final_accept is True:
            cutoff_score_list.append(cutoff_dict[idx][0])
            psm_cutoff_num = str(psm_final_F_num) + '/' + str(psm_final_TF_num)
            psm_cutoff_per = divide(float(psm_final_F_num), float(psm_final_TF_num))
            psm_cutoff_F_sum += psm_final_F_num
            psm_cutoff_TF_sum += psm_final_TF_num
            psm_cutoff_num_list.append(psm_cutoff_num)
            psm_cutoff_per_list.append(psm_cutoff_per)
        # else if cannot find score cutoff for given FDR
        else: 
            cutoff_score_list.append('NA')
            psm_cutoff_num = str('NA/NA')
            psm_cutoff_num = str('NA')
            psm_cutoff_num_list.append(psm_cutoff_num)
            psm_cutoff_per_list.append(psm_cutoff_per)

    # PSM decoy_match_num after filtering
    decoy_psm_after_filtering = psm_cutoff_F_sum
    total_psm_after_filtering = psm_cutoff_TF_sum
    target_psm_after_filtering = total_psm_after_filtering - decoy_psm_after_filtering

    # variables for pep_out
    pep_cutoff_num_list = []
    pep_cutoff_per_list = []
    pep_cutoff_per = 0.0
    pep_cutoff_F_sum = 0
    pep_cutoff_TF_sum = 0
    pep_cutoff_num_sum = ''

    # for pep_out
    for idx in range(charge_set):
        # cutoff_dict
        pep_final_cutoff_score = cutoff_dict[idx][0]
        pep_final_accept       = cutoff_dict[idx][1]
        pep_final_F_num        = len(pep_F_dict[idx])
        pep_final_T_num        = len(pep_T_dict[idx])
        pep_final_TF_num       = pep_final_F_num + pep_final_T_num
 

        # if final_accept is True
        if pep_final_accept is True:
            pep_cutoff_per = divide(float(pep_final_F_num) , float(pep_final_TF_num))
            pep_cutoff_num = str(pep_final_F_num) + '/' + str(pep_final_TF_num)
            pep_cutoff_F_sum += pep_final_F_num
            pep_cutoff_TF_sum += pep_final_TF_num
            pep_cutoff_num_list.append(pep_cutoff_num)
            pep_cutoff_per_list.append(pep_cutoff_per)
        # else if cannot find score cutoff for given FDR
        else: 
            pep_cutoff_num = str('NA/NA')
            pep_cutoff_per = str('NA')
            pep_cutoff_num_list.append(pep_cutoff_num)
            pep_cutoff_per_list.append(pep_cutoff_per)

    # peptide decoy_match_num after filtering
    decoy_peptides_after_filtering = pep_cutoff_F_sum
    total_peptides_after_filtering = pep_cutoff_TF_sum
    target_peptides_after_filtering = total_peptides_after_filtering - decoy_peptides_after_filtering

    # stat message
    psm_stat_msg = ""
    psm_stat_msg += "#\t#######################\n"
    psm_stat_msg += "#\t# Statistical_Results #\n"
    psm_stat_msg += "#\t#######################\n"
    psm_stat_msg += "#\t\n"
    psm_stat_msg += "#\t[Statistical_Results]\n"
    psm_stat_msg += "#\t\n"
    psm_stat_msg += "#\t# Numbers of PSM before filtering\n"
    psm_stat_msg += "#\tDecoy_PSM_Before_Filtering = " + str(decoy_psm_before_filtering) + "\n"
    psm_stat_msg += "#\tTarget_PSM_Before_Filtering = " + str(target_psm_before_filtering) + "\n"
    psm_stat_msg += "#\tTotal_PSM_Before_Filtering = " + str(total_psm_before_filtering) + "\n"
    psm_stat_msg += "#\t\n"
    psm_stat_msg += "#\t# Score cutoffs for PSM at charge state 1, 2, and 3 or higher\n"
    psm_stat_msg += "#\tScore_Cutoff{Z: +1} = " + set_float_digit(cutoff_score_list[0]) + "\n"
    psm_stat_msg += "#\tScore_Cutoff{Z: +2} = " + set_float_digit(cutoff_score_list[1]) + "\n"
    psm_stat_msg += "#\tScore_Cutoff{Z: +3 or higher} = " + set_float_digit(cutoff_score_list[2]) + "\n"
    psm_stat_msg += "#\t\n"
    psm_stat_msg += "#\t# FDRs and numbers of peptides at charge state 1, 2, and 3 or higher\n"
    psm_stat_msg += "#\t# FDR = Decoy/Total\n"
    psm_stat_msg += "#\tPSM_FDR{Z: +1} = " + set_float_digit(psm_cutoff_per_list[0]) + " "
    psm_stat_msg += "(Decoy/Total: " + str(psm_cutoff_num_list[0]) + ")\n"
    psm_stat_msg += "#\tPSM_FDR{Z: +2} = " + set_float_digit(psm_cutoff_per_list[1]) + " "
    psm_stat_msg += "(Decoy/Total: " + str(psm_cutoff_num_list[1]) + ")\n"
    psm_stat_msg += "#\tPSM_FDR{Z: +3 or higher} = " + set_float_digit(psm_cutoff_per_list[2]) + " "
    psm_stat_msg += "(Decoy/Total: " + str(psm_cutoff_num_list[2]) + ")\n"
    psm_stat_msg += "#\t\n"
    psm_stat_msg += "#\t# Numbers of PSM after filtering\n"
    psm_stat_msg += "#\tDecoy_PSM_After_Filtering = " + str(decoy_psm_after_filtering) + "\n"
    psm_stat_msg += "#\tTarget_PSM_After_Filtering = " + str(target_psm_after_filtering) + "\n"
    psm_stat_msg += "#\tTotal_PSM_After_Filtering = " + str(total_psm_after_filtering) + "\n"
    psm_stat_msg += "#\t\n"

    psm_out_file.write(psm_stat_msg)

    # PSM column message
    psm_column_msg = ""
    psm_column_msg += "#\t################\n"
    psm_column_msg += "#\t# Column Names #\n"
    psm_column_msg += "#\t################\n"
    psm_column_msg += "#\t\n"
    psm_column_msg += "#\t[Column Names]\n"
    psm_column_msg += "#\t\n"
    psm_column_msg += "#\tFilename = Filename of input FT2 file\n"
    psm_column_msg += "#\tScanNumber = Scan number of the PSM\n"
    psm_column_msg += "#\tParentCharge = Charge state of the PSM\n"
    psm_column_msg += "#\tMeasuredParentMass = Measured parent mass\n"
    psm_column_msg += "#\tCalculatedParentMass = Calculated parent mass from peptide sequence\n"
    psm_column_msg += "#\tMassErrorDa = Mass error in Da with 1-Da error correction\n"
    psm_column_msg += "#\tMassErrorPPM = Mass error in PPM with 1-Da error correction\n"
    psm_column_msg += "#\tScanType = Scan type of the PSM\n"
    psm_column_msg += "#\tSearchName = Sipros search name\n"
    psm_column_msg += "#\tScoringFunction = Scoring function used in the search\n"
    psm_column_msg += "#\tScore = Score\n"
    psm_column_msg += "#\tDeltaZ = Difference between the best PSM score and the next best PSM of this scan\n"
    psm_column_msg += "#\tDeltaP = Difference between the best modified PSM and its PTM isoform\n"
    psm_column_msg += "#\tIdentifiedPeptide = Identified peptide sequence with potential PTMs and mutations\n"
    psm_column_msg += "#\tOriginalPeptide = Original peptide sequence in the FASTA file\n"
    psm_column_msg += "#\tProteinNames = Names of proteins of the peptide\n"
    psm_column_msg += "#\tProteinCount = Number of proteins that the peptide can be assigned to\n"
    psm_column_msg += "#\tTargetMatch = T for target match and F for decoy match\n"
    if sipros4_input:
        psm_column_msg += "#\tAveAtom% = Average of atom% estimate\n"
        psm_column_msg += "#\tStdAtom% = Standard deviation of atom% estimate\n"
    psm_column_msg += "#\t\n"

    psm_out_file.write(psm_column_msg)

    # stat message
    pep_stat_msg = ""
    pep_stat_msg += "#\t#######################\n"
    pep_stat_msg += "#\t# Statistical_Results #\n"
    pep_stat_msg += "#\t#######################\n"
    pep_stat_msg += "#\t\n"
    pep_stat_msg += "#\t[Statistical_Results]\n"
    pep_stat_msg += "#\t\n"
    pep_stat_msg += "#\t# Numbers of peptides before filtering\n"
    pep_stat_msg += "#\tDecoy_Peptides_Before_Filtering = " + str(decoy_peptides_before_filtering) + "\n"
    pep_stat_msg += "#\tTarget_Peptides_Before_Filtering = " + str(target_peptides_before_filtering) + "\n"
    pep_stat_msg += "#\tTotal_Peptides_Before_Filtering = " + str(total_peptides_before_filtering) + "\n"
    pep_stat_msg += "#\t\n"
    pep_stat_msg += "#\t# Score cutoffs for peptides at charge state 1, 2, and 3 or higher\n"
    pep_stat_msg += "#\tScore_Cutoff{Z: +1} = " + set_float_digit(cutoff_score_list[0]) + "\n"
    pep_stat_msg += "#\tScore_Cutoff{Z: +2} = " + set_float_digit(cutoff_score_list[1]) + "\n"
    pep_stat_msg += "#\tScore_Cutoff{Z: +3 or higher} = " + set_float_digit(cutoff_score_list[2]) + "\n"
    pep_stat_msg += "#\t\n"
    pep_stat_msg += "#\t# FDRs and numbers of peptides at charge state 1, 2, and 3 or higher\n"
    pep_stat_msg += "#\t# FDR = Decoy/Total\n"
    pep_stat_msg += "#\tPeptide_FDR{Z: +1} = " + set_float_digit(pep_cutoff_per_list[0]) + " "
    pep_stat_msg += "(Decoy/Total: " + str(pep_cutoff_num_list[0]) + ")\n"
    pep_stat_msg += "#\tPeptide_FDR{Z: +2} = " + set_float_digit(pep_cutoff_per_list[1]) + " "
    pep_stat_msg += "(Decoy/Total: " + str(pep_cutoff_num_list[1]) + ")\n"
    pep_stat_msg += "#\tPeptide_FDR{Z: +3 or higher} = " + set_float_digit(pep_cutoff_per_list[2]) + " "
    pep_stat_msg += "(Decoy/Total: " + str(pep_cutoff_num_list[2]) + ")\n"
    pep_stat_msg += "#\t\n"
    pep_stat_msg += "#\t# Numbers of peptides after filtering\n"
    pep_stat_msg += "#\tDecoy_Peptides_After_Filtering = " + str(decoy_peptides_after_filtering) + "\n"
    pep_stat_msg += "#\tTarget_Peptides_After_Filtering = " + str(target_peptides_after_filtering) + "\n"
    pep_stat_msg += "#\tTotal_Peptides_After_Filtering = " + str(total_peptides_after_filtering) + "\n"
    pep_stat_msg += "#\t\n"

    pep_out_file.write(pep_stat_msg)

    # pep column message
    pep_column_msg = ""
    pep_column_msg += "#\t################\n"
    pep_column_msg += "#\t# Column Names #\n"
    pep_column_msg += "#\t################\n"
    pep_column_msg += "#\t\n"
    pep_column_msg += "#\t[Column Names]\n"
    pep_column_msg += "#\t\n"
    pep_column_msg += "#\tIdentifiedPeptide = Identified peptide sequence with potential PTMs and mutations\n"
    pep_column_msg += "#\tParentCharge = Charge state of identified peptide\n"
    pep_column_msg += "#\tOriginalPeptide = Original peptide sequence in the FASTA file\n"
    pep_column_msg += "#\tProteinNames = Names of proteins of the peptide\n"
    pep_column_msg += "#\tProteinCount = Number of proteins that the peptide can be assigned to\n"
    pep_column_msg += "#\tTargetMatch = T for target match and F for decoy match\n"
    pep_column_msg += "#\tSpectralCount = Number of PSMs in which the peptide is identified\n"
    pep_column_msg += "#\tBestScore = The best score of those PSMs\n"
    pep_column_msg += "#\tPSMs = List of PSMs for the peptide: FT2_Filename[Scan_Number]\n"
    pep_column_msg += "#\tScanType = Scan type of those PSMs\n"
    pep_column_msg += "#\tSearchName = Sipros search name\n"
    pep_column_msg += "#\t\n"

    pep_out_file.write(pep_column_msg)

    # for psm out
    psm_out_list = ['Filename',             #0
                    'ScanNumber',           #1
                    'ParentCharge',         #2
                    'MeasuredParentMass',   #3
                    'CalculatedParentMass', #4
                    'MassErrorDa',          #5 CalculatedParentMass - MeasuredParentMass
                    'MassErrorPPM',         #6 MassErrorDa / CalculatedParentMass
                    'ScanType',             #7
                    'SearchName',           #8
                    'ScoringFunction',      #9
                    'Score',                #10
                    'DeltaZ',               #11 the difference score between the rank 1 and 2
                    'DeltaP',               #12
                    'IdentifiedPeptide',    #13
                    'OriginalPeptide',      #14
                    'ProteinNames',         #15
                    'ProteinCount',         #16
                    'TargetMatch']          #17
    if sipros4_input:
        psm_out_list.append('AveAtom%')
        psm_out_list.append('StdAtom%')

    # for pep out
    pep_out_list = ['IdentifiedPeptide',    #0
                    'ParentCharge',         #1
                    'OriginalPeptide',      #2
                    'ProteinNames',         #3
                    'ProteinCount',         #4
                    'TargetMatch',          #5
                    'SpectralCount',        #6 number of PSMs matched to this peptide
                    'BestScore',            #7 the highest score of those PSMs
                    'PSMs',                 #8 a list of PSMs matched to this peptide. Use{Filename[ScanNumber],Filename[ScanNumber]} format
                    'ScanType',             #9
                    'SearchName']           #10

    # write list to tab-del fil with join 
    psm_out_file.write('\t'.join(psm_out_list) + '\n')
    pep_out_file.write('\t'.join(pep_out_list) + '\n')

    # open to read psm_data_file
    psm_data_f = open(psm_data_filename, 'rb')
    psm_data_r = csv.reader(psm_data_f, delimiter='\t')

    # get data
    for line_idx, line in enumerate(psm_data_r):

        # check sipros3 or sipros4 input
        data_one_obj = None
        if sipros4_input:
            data_one_obj = Psm4OutFields._make(line)
        else:
            data_one_obj = PsmOutFields._make(line)

        # data calss to variables
        protein_names = data_one_obj.ProteinNames
        charge_val = int(data_one_obj.ParentCharge)
        if charge_val >= charge_set:
            charge_val = charge_set
        charge_idx = charge_val - 1
        score = float(data_one_obj.Score)

        # for decoy_match_num
        psm_match_TF = check_decoy_match(protein_names, decoy_prefix)
        if psm_match_TF is False:
            psm_decoy_match_num += 1

        # get final_cutoff_score and final_accept with charge index
        final_cutoff_score = cutoff_dict[charge_idx][0]
        final_accept       = cutoff_dict[charge_idx][1]

        # if score exist
        if final_accept is True:
            # if score >= cutoff_score
            if score >= final_cutoff_score:
                # replace
                psm_out_list = [data_one_obj.Filename,
                                data_one_obj.ScanNumber,
                                data_one_obj.ParentCharge,
                                data_one_obj.MeasuredParentMass,
                                data_one_obj.CalculatedParentMass,
                                data_one_obj.MassErrorDa,
                                data_one_obj.MassErrorPPM,
                                data_one_obj.ScanType,
                                data_one_obj.SearchName,
                                data_one_obj.ScoringFunction,
                                data_one_obj.Score,
                                data_one_obj.DeltaZ,
                                data_one_obj.DeltaP,
                                data_one_obj.IdentifiedPeptide,
                                data_one_obj.OriginalPeptide,
                                data_one_obj.ProteinNames,
                                data_one_obj.ProteinCount,
                                data_one_obj.TargetMatch]
                if sipros4_input:
                    psm_out_list.append(data_one_obj.AveAtom)
                    psm_out_list.append(data_one_obj.StdAtom)

                psm_out_file.write('\t'.join(psm_out_list) + '\n')

    # close psm_data_f
    psm_data_f.close()

    # open to read pep_data_file
    pep_data_f = open(pep_data_filename, 'rb')
    pep_data_r = csv.reader(pep_data_f, delimiter='\t')

    # get data
    for line_idx, line in enumerate(pep_data_r):

        # call PsmOutFields class
        pep_data_obj = PepDataFields._make(line)

        # get pep_ID
        pep_ID = pep_data_obj.IdentifiedPeptide.strip() + '_+_' + pep_data_obj.ParentCharge.strip()

        # data calss to variables
        protein_names = pep_data_obj.ProteinNames
        charge_val = int(pep_data_obj.ParentCharge)
        if charge_val >= charge_set:
            charge_val = charge_set
        charge_idx = charge_val - 1
        best_score = float(pep_data_obj.BestScore)

        # for total_match_num and decoy_match_num
        pep_match_TF = check_decoy_match(protein_names, decoy_prefix)
        if pep_match_TF is False:
            pep_decoy_match_num += 1

        # get final_cutoff_score and final_accept with charge index
        final_cutoff_score = cutoff_dict[charge_idx][0]
        final_accept       = cutoff_dict[charge_idx][1]

        # if score exist
        if final_accept is True:

            # if score >= cutoff_score
            if best_score >= final_cutoff_score:

                # for pep_out_list
                pep_out_identified_peptide = pep_data_obj.IdentifiedPeptide
                pep_out_parent_charge = pep_data_obj.ParentCharge

                # initialize
                pep_out_original_peptide_list = []
                pep_out_protein_names_list = []
                pep_out_spectral_count = 0
                pep_out_psms_list = []
                pep_out_scan_type_list = []
                pep_out_search_name_list = []
                pep_out_scan_type_dict = {}
                pep_out_search_name_dict = {}

                # for loop of pep_sub_dict[pep_ID]
                for pep_sub_list in pep_sub_dict[pep_ID]:
 
                    # get peb_sub_obj
                    pep_sub_obj = PepSubFields._make(pep_sub_list)

                    # check score >=final_cutoff_score
                    if float(pep_sub_obj.Score) >= float(final_cutoff_score):
                        pep_out_original_peptide_list.append(pep_sub_obj.OriginalPeptide)
                        protein_names_wo_bracket = pep_sub_obj.ProteinNames.strip()[1:-1]
                        pep_out_protein_names_list.append(protein_names_wo_bracket)
                        pep_out_spectral_count += 1
                        pep_out_psm = pep_sub_obj.Filename + '[' + pep_sub_obj.ScanNumber + ']'
                        pep_out_psms_list.append(pep_out_psm)
                        pep_out_scan_type_dict[pep_out_psm] = pep_sub_obj.ScanType
                        pep_out_search_name_dict[pep_out_psm] = pep_sub_obj.SearchName

                # get unique list of list
                pep_out_original_peptide_list = list(set(pep_out_original_peptide_list))
                pep_out_protein_names_list = list(set(pep_out_protein_names_list))
                pep_out_psms_list = list(set(pep_out_psms_list))
                for pep_out_psms_list_item in pep_out_psms_list:
                    pep_out_scan_type_list.append(pep_out_scan_type_dict[pep_out_psms_list_item])
                    pep_out_search_name_list.append(pep_out_search_name_dict[pep_out_psms_list_item])

                # convert to string
                pep_out_original_peptide_str = list_to_string(pep_out_original_peptide_list)
                pep_out_protein_names_str = list_to_bracket(pep_out_protein_names_list)
                pep_out_protein_count_str = str(get_protein_count(pep_out_protein_names_str))
                pep_out_target_match_val = check_decoy_match(pep_out_protein_names_str, decoy_prefix)
                pep_out_target_match_str = 'T' if pep_out_target_match_val is True else 'F' 
                pep_out_psms_str = list_to_bracket(pep_out_psms_list)
                pep_out_scan_type_str = list_to_bracket(pep_out_scan_type_list)
                pep_out_search_name_str = list_to_bracket(pep_out_search_name_list)

                pep_out_list = [str(pep_out_identified_peptide),
                                str(pep_out_parent_charge),
                                str(pep_out_original_peptide_str),
                                str(pep_out_protein_names_str),
                                str(pep_out_protein_count_str),
                                str(pep_out_target_match_str),
                                str(pep_out_spectral_count),
                                str(best_score),
                                str(pep_out_psms_str),
                                str(pep_out_scan_type_str),
                                str(pep_out_search_name_str)]
           
                # append to the pep_out_data list
                pep_out_file.write('\t'.join(pep_out_list) + '\n')

    # close psm_data_f
    pep_data_f.close()
    total_peptides_before_filtering = line_idx + 1

    # delete psm_data_file and pep_data_file
    os.remove(psm_data_filename)
    os.remove(pep_data_filename)


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
            sys.stderr.write('[%s] Beginning sipros_peptides_filtering.py run (V%s)\n' % (curr_time(), get_version()))
            sys.stderr.write('------------------------------------------------------------\n')

            # Parse options and get config file
            sys.stderr.write('[Step 1] Parse options and get config file: Running -> ')
            # Call parse_config to open and read config file
            config_dict = parse_config(config_filename)
            # Get sipros output file(s) in working directory
            sipros_file_list = get_file_list_with_ext(working_dir, sipros_file_ext)
            # Get base_out for output
            base_out_default = 'Sipros_searches'
            base_out = get_base_out(sipros_file_list, base_out_default, working_dir)
            # Open output files
            try:
                global psm_out_file
                global pep_out_file
                psm_out_file = open(base_out + ".psm.txt", 'wb')
                pep_out_file = open(base_out + ".pep.txt", 'wb')
            except:
                print >> sys.stderr, '\nCannot write output file! Check permission of working directory.'
                die("Program exit!")
            # for tmp file
            global psm_data_filename
            global pep_data_filename
            psm_data_filename = base_out + ".psm_data.tmp"
            pep_data_filename = base_out + ".pep_data.tmp"

            sys.stderr.write('Done!\n')

            # Read sipros output files and get data
            sys.stderr.write('[Step 2] Load sipros output file(s):        Running -> ')
            pep_sub_dict = read_sipros_files(sipros_file_list, config_dict, base_out)
            sys.stderr.write('Done!\n')
 
            # Calculate cuff-off score by FDR threshold
            sys.stderr.write('[Step 3] Search cut-off score by given FDR: Running -> ')
            cutoff_dict = get_cutoff_score(config_dict)
            sys.stderr.write('Done!\n')

            # Report output
            sys.stderr.write('[Step 4] Report .psm.txt and .pep.txt:      Running -> ')
            # Report output files
            report_output(sipros_file_list, pep_sub_dict, cutoff_dict, config_dict)
            sys.stderr.write('Done!\n')

            # Time record, calculate elapsed time, and display work end
            finish_time = datetime.now()
            duration = finish_time - start_time
            sys.stderr.write('------------------------------------------------------------\n')
            sys.stderr.write('[%s] Ending sipros_peptides_filtering.py run\n' % curr_time())
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
