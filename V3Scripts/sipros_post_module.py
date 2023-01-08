#!/usr/bin/python

"""
sipros_post_module.py

sipros_post_module.py module includes common classes, definitions,
and funcitons for sipros post-processing programs.

Created by Tae-Hyuk (Ted) Ahn on 10/10/2012.
Copyright (c) 2012 Tae-Hyuk Ahn (ORNL). Allrights reserved.
"""

## Import standard modules
import sys, warnings, os, re
from datetime import datetime, date, time
from collections import namedtuple
import math


## Class Usage
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


## Class for ignoring comments '#' in sipros file
class CommentedFile:
    def __init__(self, f, comment_string = "#"):
        self.f = f
        self.comment_string = comment_string
    def next(self):
        line = self.f.next()
        while line.startswith(self.comment_string):
            line = self.f.next()
        return line
    def __iter__(self):
        return self


## Exit system with error message
def die(msg=None):
    if msg is not None:
        print >> sys.stderr, msg
        sys.exit(1)


## Returns the current time in a nice format
def curr_time():
    curr_time = datetime.now()
    return curr_time.strftime("%c")


## Format time as a pretty string
def format_time(td):
    hours = td.seconds // 3600
    minutes = (td.seconds % 3600) // 60
    seconds = td.seconds % 60
    return '%02d:%02d:%02d' % (hours, minutes, seconds)


## Find string between two substrings
def find_between(s, first, last):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""


## Division error handling
def divide(x, y):
    try:
        result = x / y
    except ZeroDivisionError as detail:
        print >> sys.stderr, 'Handling run-time error:', detail
        die('Program exit!')
    else:
        return result


## Check file exist
def check_file_exist(filename):

    try:
        with open(filename) as f: pass
    except IOError as e:
        print >> sys.stderr, '\nCannot open', filename
        die("Program exit!")


## Get file(s) list in working dir with specific file extension
def get_file_list_with_ext(working_dir, file_ext):

    # define sipros file extension 
    file_list = []

    # working directory
    if os.path.exists(working_dir):
        for file_name in os.listdir(working_dir):

            # check the file extension
            if file_name.endswith(file_ext):
                file_path_name = working_dir + file_name
                file_list.append(file_path_name)

        if len(file_list) == 0:
            print >> sys.stderr, "\nCannot open %s file(s)." % (file_ext)
            die("Program exit!")
        file_list = sorted(file_list)

    else:
        print >> sys.stderr, "\nCannot open working directory", working_dir
        die("Program exit!")

    return file_list

## Get base output filename with input file list and base_out_default
def get_base_out(file_list, base_out_default, working_dir):

    # Get base output with common prefix 
    base_out = os.path.commonprefix(file_list)
    base_out_filename = base_out.split('/')[-1]

    # If base common prefix ends with '.pep.txt', then remove '.pep.txt'
    base_out = base_out.replace(".pep.txt", "_")

    # If base_out file name is less than 5, then use default baseout
    if len(base_out_filename) < 5:
        base_out = working_dir + base_out_default

    # If base common prefix ends with '_' or '.', then remove
    base_out = base_out[:-1] if (base_out[-1] in ('_', '.')) else base_out

    return base_out


## list_to_string 
## if single element, then just convert to string
## if multiple elements, then bracket {A,B}
def list_to_string(input_list):

    if len(input_list) > 1:
        converted_str = '{' + ','.join(input_list) + '}'
    else:
        converted_str = ''.join(input_list)

    return converted_str


## list_to_bracket 
## bracket the list
def list_to_bracket(input_list):

    converted_str = '{' + ','.join(input_list) + '}'

    return converted_str


## Class for sipros fields object
class SiprosFields(namedtuple('SiprosFields', 
        ['Filename', 
        'ScanNumber',
        'ParentCharge',
        'MeasuredParentMass',
        'CalculatedParentMass',
        'ScanType',
        'SearchName',
        'ScoringFunction',
        'Rank',
        'Score',
        'IdentifiedPeptide',
        'OriginalPeptide',
        'ProteinNames'])):
    def __init__(self):
        self.data = self


## Class for sipros4 fields object
class Sipros4Fields(namedtuple('SiprosFields', 
        ['Filename', 
        'ScanNumber',
        'ParentCharge',
        'MeasuredParentMass',
        'CalculatedParentMass',
        'ScanType',
        'SearchName',
        'ScoringFunction',
        'Rank',
        'Score',
        'IdentifiedPeptide',
        'OriginalPeptide',
        'ProteinNames',
        'AveAtom',
        'StdAtom'])):

    def __init__(self):
        self.data = self


## Class for PsmOutFields object
class PsmOutFields(namedtuple('PsmOutFields', 
        ['Filename',
         'ScanNumber',
         'ParentCharge',
         'MeasuredParentMass',
         'CalculatedParentMass',
         'MassErrorDa',     # CalculatedParentMass - MeasuredParentMass
         'MassErrorPPM',    # MassErrorDa / CalculatedParentMass
         'ScanType',
         'SearchName',
         'ScoringFunction',
         'Score',
         'DeltaZ',          # The difference score between the rank 1 and 2
         'DeltaP',          # The difference score between isoform
         'IdentifiedPeptide',
         'OriginalPeptide',
         'ProteinNames',
         'ProteinCount',
         'TargetMatch'])):
    def __init__(self):
        self.data = self


## Class for PsmOutFields object (for sipro4)
class Psm4OutFields(namedtuple('PsmOutFields', 
        ['Filename',
         'ScanNumber',
         'ParentCharge',
         'MeasuredParentMass',
         'CalculatedParentMass',
         'MassErrorDa',     # CalculatedParentMass - MeasuredParentMass
         'MassErrorPPM',    # MassErrorDa / CalculatedParentMass
         'ScanType',
         'SearchName',
         'ScoringFunction',
         'Score',
         'DeltaZ',          # The difference score between the rank 1 and 2
         'DeltaP',          # The difference score between isoform
         'IdentifiedPeptide',
         'OriginalPeptide',
         'ProteinNames',
         'ProteinCount',
         'TargetMatch',
         'AveAtom',
         'StdAtom'])):
    def __init__(self):
        self.data = self


## Class for PepSubFields object
class PepSubFields(namedtuple('PepSubFields', 
        ['IdentifiedPeptide',
         'ParentCharge',
         'OriginalPeptide',
         'ProteinNames',
         'Score',
         'Filename',
         'ScanNumber',
         'ScanType',
         'SearchName'])):

    def __init__(self):
        self.data = self


## Class for PepSubFields object
class PepDataFields(namedtuple('PepDataFields', 
        ['IdentifiedPeptide',
         'ParentCharge',
         'BestScore',
         'ProteinNames'])):

    def __init__(self):
        self.data = self


## Class for PepOutFields object
class PepOutFields(namedtuple('PepOutFields', 
        ['IdentifiedPeptide',    #0
         'ParentCharge',         #1
         'OriginalPeptide',      #2
         'ProteinNames',         #3
         'ProteinCount',         #4
         'TargetMatch',          #5
         'SpectralCount',        #6 number of PSMs matched to this peptide
         'BestScore',            #7 the highest score of those PSMs
         'PSMs',                 #8 a list of PSMs matched to this peptide. Use {Filename[ScanNumber],Filename[ScanNumber]} format
         'ScanType',             #9 ScanType
         'SearchName'])):        #10 SearchName
 
    def __init__(self):
        self.data = self


## Class for pretty float
class PrettyFloat(float):
    def __repr__(self):
        return "%0.5f" % self


## A range function, that does accept float increments
def frange(start, end=None, inc=None):

    if end == None:
        end = start + 0.0
        start = 0.0

    if inc == None:
        inc = 1.0

    L = []
    while 1:
        next = start + len(L) * inc
        if inc > 0 and next >= end:
            break
        elif inc < 0 and next <= end:
            break
        L.append(next)
        
    return L


## check sub list
def check_sub_list(list_A, list_B):

    check_status = True

    for list_A_item in list_A:
        if list_A_item not in list_B:
            check_status = False
        else:
            continue

    return check_status


## get item list from parenthesis string as {AA,BB}
def get_item_list(input_string):

    input_string = input_string[1:-1]
    item_list = re.split(r"\s*[,]\s*", input_string.strip())

    return item_list


## get_protein_count
def get_protein_count(protein_names):
 
    input_string = protein_names[1:-1]
    item_list = re.split(r"\s*[,]\s*", input_string.strip())
    protein_count = len(item_list)

    return protein_count


## set float digit
def set_float_digit(input_val):

    if input_val is float:
        output_val = str("{0:.5f}".format(round(input_val,5)))
    else:
        output_val = str(input_val)

    return output_val

## peptide delete residues
def peptide_delete_residues(peptide_string):
    
    try: 
        left_braket_index = peptide_string.index('[')
        right_braket_index = peptide_string.index(']')
        if len(peptide_string) > right_braket_index + 1:
            if peptide_string[right_braket_index + 1].isalpha():
                peptide_output = peptide_string[left_braket_index:right_braket_index+1]
            else:
                peptide_output = peptide_string[left_braket_index:right_braket_index+2]
        else:
           peptide_output = peptide_string[left_braket_index:right_braket_index+1]

        return peptide_output
    except AttributeError:
        print >> sys.stderr, '\nCannot parse peptide correctly.\n'
        die("Program exit!")


## merge protein names
def merge_protein_names(first_protein_names, second_protein_names):

    first_protein_list = get_item_list(first_protein_names)
    second_protein_list = get_item_list(second_protein_names)

    merge_protein_list = list(set(first_protein_list + second_protein_list))

    merge_protein_names = list_to_bracket(merge_protein_list)

    return merge_protein_names
