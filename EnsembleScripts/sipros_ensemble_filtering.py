'''
Created on Sep 7, 2016

@author: xgo
'''

import getopt, sys, os
import numpy as np
import csv
import math
import re

try:
    from sets import Set
except ImportError:
    pass

from datetime import datetime, date, time
from collections import namedtuple
from sklearn import linear_model
from sklearn import preprocessing
from subprocess import call
from multiprocessing import Process
from multiprocessing import Queue, cpu_count

## Import Sipros package modules
import sipros_post_module
import sipros_peptides_assembling
import parseconfig

## Returns the current time in a nice format
curr_time = sipros_post_module.curr_time
## Format time as a pretty string
format_time = sipros_post_module.format_time
## get the file extension
get_file_list_with_ext =  sipros_post_module.get_file_list_with_ext
## Class for ignoring comments '#' in sipros file
CommentedFile = sipros_post_module.CommentedFile

## global variables
## training prefix
train_str = '' 
## testing prefix
test_str = '' 
## reserved prefix
reserve_str = ''
## ratio of testing decoy vs forward
Test_Fwd_Ratio = 1
## maximum precurer mass windows
mass_window_max_int = 0
## feature list
feature_selection_list = []
## fwd psm value
LabelFwd = 1
## training psm value
LabelTrain = 2
## testing psm value
LabelTest = 3
## reserved psm value
LabelReserve = 4
## sip mode training psm value
LabelSipTrainFwd = 1
## forward psms
num_forward_psms_before_filtering = 0
## protein database size
num_proteins = 0

## Class for PepOutFields object
class PsmFields4(namedtuple('PsmFields',
        ['FileName',  # 0
         'ScanNumber',  # 1
         'ParentCharge',  # 2
         'MeasuredParentMass',  # 3
         'ScanType',  # 4
         'SearchName',  # 5
         'IdentifiedPeptide',  # 6
         'OriginalPeptide',  # 7
         'CalculatedParentMass',  # 8
         'MVH',  # 9
         'Xcorr',  # 10
         'WDP',  # 11
         'ProteinNames',  # 12
         'DiffMVH',
         'DiffXcorr',
         'DiffWDP',
         'RetentionTime',
         'DeltaP'])): # 33 ,
    def __init__(self):
        self.data = self

## save filename in this list, for saving memory
filename_list = []
## get the filename index
def get_set_filename(filename):
    if filename in filename_list:
        return filename_list.index(filename)
    else:
        filename_list.append(filename)
        return filename_list.index(filename)

## save scantype in this list, for saving memory
scantype_list = []
## get the scan type index
def get_set_scantype(scantype):
    if scantype in scantype_list:
        return scantype_list.index(scantype)
    else:
        scantype_list.append(scantype)
        return scantype_list.index(scantype)

## save search name in this list, for saving memory
searchname_list = []
## get the index for given search name
def get_set_searchname(searchname):
    if searchname in searchname_list:
        return searchname_list.index(searchname)
    else:
        searchname_list.append(searchname)
        return searchname_list.index(searchname)    

## get the percentage
def get_percentage_back(psm_list):
    l = []
    for e in searchname_list:
        begin_i = e.index('_')
        end_i = e.index('Pct')
        pct_s = e[begin_i+1:end_i]
        l.append(pct_s)
    for psm_o in psm_list:
        psm_o.pct_s = l[psm_o.SearchName]

## class defining the psm
class PSM:
    
    # number of scores
    iNumScores = 3
    # Neutron mass
    fNeutronMass = 1.00867108694132
    # pattern for getting original peptides
    pattern = re.compile('[^\w\[\]]')

    def __init__(self, psm_field):
        self.FileName = get_set_filename(psm_field.FileName)
        self.ScanNumber = int(psm_field.ScanNumber)
        self.ParentCharge = int(psm_field.ParentCharge)
        self.ScanType = get_set_scantype(psm_field.ScanType)
        self.SearchName = get_set_searchname(psm_field.SearchName)
        self.lfScores = [float(psm_field.MVH), float(psm_field.Xcorr), float(psm_field.WDP)]
        self.ProteinNames = psm_field.ProteinNames.strip()
        self.IdentifiedPeptide = psm_field.IdentifiedPeptide
        s1 = ''.join([char if char.isalnum() else '$' for char in self.IdentifiedPeptide ])
        self.PTMscore = s1.count('$') - 2
        self.OriginalPeptide = psm_field.OriginalPeptide
        self.OriginalPeptide = PSM.pattern.sub('', self.IdentifiedPeptide)
        self.protein_list = []
        self.RealLabel = get_protein_type(self.ProteinNames, self.protein_list)
        self.fPredictProbability = 0.0
        self.fMassDiff = 0.0
        self.dM = 0.0
        self.MeasuredParentMass = float(psm_field.MeasuredParentMass)
        self.CalculatedParentMass = float(psm_field.CalculatedParentMass)
        self.set_mass_diff(self.MeasuredParentMass, self.CalculatedParentMass)
        
        self.score_differential_list = []
        self.iLocalRank = 0
        self.DeltaP = 'NA'
        if type(psm_field).__name__ == 'PsmFields4':
            self.score_differential_list = [float(psm_field.DiffMVH), float(psm_field.DiffXcorr), float(psm_field.DiffWDP)]
            self.DeltaP = psm_field.DeltaP
        else:
            print('not support input format.')
            sys.exit(1)
        
        self.NMC = 0
        self.IPSC = 0
        self.OPSC = 0
        self.UPSC = 0 # unique peptide
        self.SPSC = 0 # shared peptide
        self.NRS = 0
        self.PPC = 0
        self.UPPC = 0
        self.SPPC = 0
        self.feature_list = []
        
        self.TrainingLabel = 0
        
        self.pct_s = ''
    
    # extract 10 features
    def get_feature_final_list(self):
        del self.feature_list[:]
        self.feature_list.extend(self.lfScores) # 2, 3, 4: 1, 2, 3
        self.feature_list.append(abs(self.fMassDiff)) # 6: 5
        self.feature_list.extend(self.score_differential_list) # 7 - 24: 6 - 23
        self.feature_list.append(self.NMC) # 25: 24
        self.feature_list.append((self.OPSC)) # 27: 26
        self.feature_list.append((self.SPSC)) # 29: 28
    '''
    def get_feature_list(self):
        del self.feature_list[:]
        self.feature_list.append(self.ParentCharge) # 1: 0
        self.feature_list.extend(self.lfScores) # 2, 3, 4: 1, 2, 3
        self.feature_list.append(self.ScoreAgreement) # 5: 4
        self.feature_list.append(abs(self.fMassDiff)) # 6: 5
        self.feature_list.extend(self.score_differential_list) # 7 - 24: 6 - 23
        self.feature_list.append(self.NMC) # 25: 24
        self.feature_list.append((self.IPSC)) # 26: 25
        self.feature_list.append((self.OPSC)) # 27: 26
        self.feature_list.append((self.UPSC)) # 28: 27
        self.feature_list.append((self.SPSC)) # 29: 28
        self.feature_list.append(abs(self.iMassWindow)) # 30: 29
        self.feature_list.append((self.PPC)) # 31: 30
        
        for c in ptm_selection_list:
            self.feature_list.append(self.IdentifiedPeptide.count(ptm_str[c])) # 32: 31
    '''
    
    # put proteins inside {}
    def set_protein_names(self):
        self.ProteinNames = '{' + ','.join(self.protein_list) + '}'
    
    # add protein to psm, in case some protein missing
    def add_protein(self, protein_l):
        add_bool = False
        for p in protein_l:
            if p not in self.protein_list:
                add_bool = True
                self.protein_list.append(p)
        
        if add_bool:
            self.set_protein_names()
    
    # get the mass difference, considering mass windows 
    def set_mass_diff(self, measured_mass, calculated_mass):
        fDiff = calculated_mass - measured_mass
        fTemp = fDiff
        fCeil = 0
        fDown = 0
        if fDiff >= 0:
            fDiff = fTemp
            fCeil = math.ceil(fTemp)*PSM.fNeutronMass
            fFloor = math.floor(fTemp)*PSM.fNeutronMass
            if fFloor > fTemp:
                fFloor -= PSM.fNeutronMass
            if fCeil - PSM.fNeutronMass > fTemp:
                fCeil -= PSM.fNeutronMass
            if fTemp > fCeil - fTemp:
                fTemp = fCeil - fTemp
            if fDiff > fDiff - fFloor:
                fDiff = abs(fDiff - fFloor)
            if abs(fTemp) < abs(fDiff):
                fDiff = fTemp
                self.dM = -fTemp
            else:
                self.dM = fDiff
        else:
            fCeil = math.ceil(fDiff)*PSM.fNeutronMass
            if fCeil < fDiff:
                fCeil += PSM.fNeutronMass
            fFloor = math.floor(fDiff)*PSM.fNeutronMass
            if fFloor + PSM.fNeutronMass < fDiff:
                fFloor += PSM.fNeutronMass
            fDiff = fTemp
            if abs(fTemp) > fCeil - fTemp:
                fTemp = fCeil - fTemp
            if abs(fDiff) > fDiff - fFloor:
                fDiff = fDiff - fFloor
            fTemp = abs(fTemp)
            fDiff = abs(fDiff)
            if fTemp < fDiff:
                fDiff = fTemp
                self.dM = -fTemp
            else:
                self.dM = fDiff
        self.fMassDiff = fDiff
    
    # remove training proteins and reserved proteins
    def clean_protein_name(self):
        self.ProteinNames = ""
        l = []
        if not reserve_str == "":
            for sProtein in self.protein_list:
                sProtein.strip()
                if not (sProtein.startswith(reserve_str)):
                    l.append(sProtein)
            self.protein_list = l
            l = []
        
        for sProtein in self.protein_list:
            sProtein.strip()
            if train_str == "":
                if sProtein not in l:
                    l.append(sProtein)
            elif not (sProtein.startswith(train_str)):
                if sProtein not in l:
                    l.append(sProtein)
        self.ProteinNames = '{'+','.join(l) + '}'
        self.protein_list = l
    
    # sip mode
    def set_real_label(self):
        self.RealLabel = get_protein_type(self.ProteinNames, self.protein_list)


# # Version control
def get_version():
    return "Sipros Ensemble 1.0.1 (Alpha)"

# # Help message
help_message = '''
Usage:
    python sipros_ensemble_filtering.py [options]

Inputs:
    -i PSM.tab
    -c Sipros Ensemble configuration file

Options:
    -h show help info
    -v show version info

Outputs:
    -o output directory
'''

# # Parse options
def parse_options(argv):
    
    try:
        opts, _args = getopt.getopt(argv[1:], "hvi:c:o:x:")
    except getopt.GetoptError:
        print("illigal option(s)")
        print(help_message)
        sys.exit(0)

    # Default working dir and config file
    input_file = ""
    output_folder = ""
    config_file = ""
    debug_code = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print(help_message)
            sys.exit(0)
        elif option in ("-v", "-V", "--version"):
            print("{} version {}".format(__file__, get_version()))
            sys.exit(0)
        elif option in ("-i"):
            input_file = value
        elif option in ("-o"):
            output_folder = value
        elif option in ("-c"):
            config_file = value
        elif option in ("-x"):
            debug_code = value
            
    if input_file == "" or output_folder == "" or config_file == '':
        print(help_message)
        sys.exit(0)

    output_folder = os.path.join(output_folder, '')

    return (input_file, config_file, output_folder, debug_code)

# # Decoy Reverse Forward protein
def protein_type(protein_sequence, lProtein=None):
    sProteins = protein_sequence.replace('{', '')
    sProteins = sProteins.replace('}', '')
    asProteins = sProteins.split(',')
    if lProtein != None:
        del lProtein[:]
        for sProtein in asProteins:
            sProtein = sProtein.strip()
            if sProtein not in lProtein:
                lProtein.append(sProtein)
    
    if reserve_str != '':
        if train_str != '':
            for sProtein in asProteins:
                if not (sProtein.startswith(train_str) or sProtein.startswith(test_str) or sProtein.startswith(reserve_str)):
                    return LabelFwd
        else:
            for sProtein in asProteins:
                if not (sProtein.startswith(test_str) or sProtein.startswith(reserve_str)):
                    return LabelFwd            
    else:
        if train_str != '':
            for sProtein in asProteins:
                if not (sProtein.startswith(train_str) or sProtein.startswith(test_str)):
                    return LabelFwd
        else:
            for sProtein in asProteins:
                if not (sProtein.startswith(test_str)):
                    return LabelFwd
    
    if test_str != '':
        for sProtein in asProteins:
            if sProtein.startswith(test_str):
                return LabelTest
    
    if reserve_str != '':
        for sProtein in asProteins:
            if sProtein.startswith(reserve_str):
                return LabelReserve
    
    return LabelTrain

def get_protein_type(protein_sequence, lProtein=None):
    """
    get the protein type
    if all reserved type, return LabelReserve
    if all testing type, return LabelTest
    if all training type, return LabelTrain
    otherwise, it is forward protein, return LabelFwd
    """
    sProteins = protein_sequence.replace('{', '')
    sProteins = sProteins.replace('}', '')
    asProteins = sProteins.split(',')
    if lProtein != None:
        del lProtein[:]
        for sProtein in asProteins:
            sProtein = sProtein.strip()
            if sProtein not in lProtein:
                lProtein.append(sProtein)
                
    protein_list_tmp_1 = []
    protein_list_tmp_2 = []
    
    reserve_type = True
    if reserve_str != '':
        for sProtein in asProteins:
            if not sProtein.startswith(reserve_str):
                protein_list_tmp_1.append(sProtein)
                reserve_type = False
        if reserve_type:
            return LabelReserve
    else:
        protein_list_tmp_1.extend(asProteins)
    
    training_type = True
    if train_str != '':
        for sProtein in protein_list_tmp_1:
            if not sProtein.startswith(train_str):
                protein_list_tmp_2.append(sProtein)
                training_type = False
        if training_type:
            return LabelTrain
    else:
        protein_list_tmp_2.extend(protein_list_tmp_1)
    
    testing_type = True
    if test_str != '':
        for sProtein in protein_list_tmp_2:
            if not sProtein.startswith(test_str):
                testing_type = False
        if testing_type:
            return LabelTest
    
    return LabelFwd

# # read the psm table
def read_psm_table(input_file):
    
    sip_files_list = []
    
    # check if it is a folder
    if os.path.isdir(input_file) == True:
        lFileList = get_file_list_with_ext(input_file, '.tab')
        for sFileName in lFileList:
            sip_files_list.append(sFileName)
    elif ',' in input_file:
        file_split_list = input_file.split(',')
        for file_str in file_split_list:
            sip_files_list.append(file_str)
    else:
        sip_files_list.append(input_file)

    # get the base name from sip file list
    base_out = os.path.commonprefix(sip_files_list)
    base_out_filename = os.path.split(base_out)[-1]
    base_out_filename = base_out_filename.replace(".tab", "_")
    
    # If base_out file name is less than 5, then use default baseout
    if len(base_out_filename) < 5:
        base_out = 'sipros_results'
    else:
        # If base common prefix ends with '_' or '.', then remove
        base_out = base_out_filename[:-1] if (base_out_filename[-1] in ('_', '.')) else base_out_filename
    
    psm_list = []
    
    # read line with csv
    for file_str in sip_files_list:
        f = csv.reader(CommentedFile(open(file_str, 'rb')), delimiter='\t')
        # skip header
        _sHeader = f.next()
        # get data
        for sLine in f:
            PsmFields_obj = PsmFields4._make(sLine)
            psm_obj = PSM(PsmFields_obj)
            psm_list.append(psm_obj)
        
    # sorting all PSM, first on file name, than scan number
    psm_list = sorted(psm_list, key=lambda psm: (psm.FileName, psm.ScanNumber))
        
    return (psm_list, base_out)

# # Division error handling
divide = sipros_post_module.divide
FDR_parameter = 1.0

# # FDR calculator
def FDR_calculator(FP, TP):
    FDR_numerator = float(FP) * float(FDR_parameter)
    FDR_denominator = float(TP)
    FDR_accept = True

    if  FDR_denominator == 0:
        FDR_value = 1.0
        FDR_accept = False
    else:
        FDR_value = divide(FDR_numerator, FDR_denominator)
        FDR_accept = True

    return (FDR_accept, float(FDR_value))

## psm level filtering
def show_Fdr(psm_list, fdr_float, charge_left_given = -1, charge_right_given = -1):
    
    # list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    list_sorted = sorted(psm_list, key=lambda x: (-x.fPredictProbability, -x.fMassDiff, -x.PTMscore, x.IdentifiedPeptide))
    num_forward = 0
    num_training = 0
    num_testing = 0
    # forward, training, testing
    best_nums = [0, 0, 0]

    
    psm_filtered_list = []
    cutoff_probability = 1000.0
    # without considering training label
    for oPsm in list_sorted:
        if charge_left_given != -1 and (oPsm.ParentCharge < charge_left_given or oPsm.ParentCharge > charge_right_given):
            continue
        if oPsm.RealLabel == LabelFwd:
            num_forward += 1
        elif oPsm.RealLabel == LabelTrain:
            num_training += 1
        elif oPsm.RealLabel == LabelTest:
            num_testing += 1
        else:
            sys.stderr.write('error 768.\n')
        (FDR_accept, FDR_value) = FDR_calculator(num_testing, num_forward)
        if (FDR_accept is True) and (FDR_value <= fdr_float):
            if (best_nums[0] + best_nums[2]) < (num_forward + num_testing):
                best_nums = [num_forward, num_training, num_testing]
                cutoff_probability = oPsm.fPredictProbability
            
    for oPsm in list_sorted:
        if charge_left_given != -1 and (oPsm.ParentCharge < charge_left_given or oPsm.ParentCharge > charge_right_given):
            continue
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
     
    return psm_filtered_list

## peptide level filtering
def show_Fdr_Pep(psm_list, fdr_float, charge_left_given = -1, charge_right_given = -1):
    
    list_sorted = sorted(psm_list, key=lambda x: (-x.fPredictProbability, -x.fMassDiff, -x.PTMscore, x.IdentifiedPeptide))
    
    peptide_set = Set()
    num_forward_pep = 0
    num_training_pep = 0
    num_testing_pep = 0
    best_forward_pep = 0
    best_testing_pep = 0
    
    psm_filtered_list = []
    cutoff_probability = 1000.0
    # without considering training label
    for oPsm in list_sorted:
        if charge_left_given != -1 and (oPsm.ParentCharge < charge_left_given or oPsm.ParentCharge > charge_right_given):
            continue
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if pep_str not in peptide_set:
            if oPsm.RealLabel == LabelFwd:
                num_forward_pep += 1
                peptide_set.add(pep_str)
            elif oPsm.RealLabel == LabelTest:
                num_testing_pep += 1
                peptide_set.add(pep_str)
            elif oPsm.RealLabel == LabelTrain:
                num_training_pep += 1
                peptide_set.add(pep_str)
            else:
                sys.stderr.write('Error 71341\n')

        (FDR_accept, FDR_value) = FDR_calculator(num_testing_pep, num_forward_pep)
        if (FDR_accept is True) and (FDR_value <= fdr_float):
            if (num_testing_pep + num_forward_pep) > (best_testing_pep + best_forward_pep):
                cutoff_probability = oPsm.fPredictProbability
                best_forward_pep = num_forward_pep
                best_testing_pep = num_testing_pep
            

    for oPsm in list_sorted:
        if charge_left_given != -1 and (oPsm.ParentCharge < charge_left_given or oPsm.ParentCharge > charge_right_given):
            continue
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
     
    return psm_filtered_list

## remove redundant psm, only one unique spectrum kept
def re_rank(psm_list, consider_charge_bool = False):
    psm_new_list = []
    psm_dict = {}
    if consider_charge_bool :
        for oPsm in psm_list:
            sId = filename_list[oPsm.FileName] + '_' + str(oPsm.ScanNumber) + '_' + str(oPsm.ParentCharge)
            if sId in psm_dict:
                if oPsm.fPredictProbability > psm_dict[sId].fPredictProbability:
                    psm_dict[sId] = oPsm
                elif oPsm.fPredictProbability == psm_dict[sId].fPredictProbability:
                    if abs(oPsm.fMassDiff) < abs(psm_dict[sId].fMassDiff): 
                        psm_dict[sId] = oPsm
                    elif abs(oPsm.fMassDiff) == abs(psm_dict[sId].fMassDiff):
                        # calculate PTM scores
                        if oPsm.PTMscore < psm_dict[sId].PTMscore:
                            psm_dict[sId] = oPsm
                        elif oPsm.PTMscore == psm_dict[sId].PTMscore:
                            if oPsm.IdentifiedPeptide.upper() < psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId] = oPsm
                            elif oPsm.IdentifiedPeptide.upper() == psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId].add_protein(oPsm.protein_list)
                    
            else:
                psm_dict[sId] = oPsm
    else :
        for oPsm in psm_list:
            sId = filename_list[oPsm.FileName] + '_' + str(oPsm.ScanNumber)
            if sId in psm_dict:
                if oPsm.fPredictProbability > psm_dict[sId].fPredictProbability:
                    psm_dict[sId] = oPsm
                elif oPsm.fPredictProbability == psm_dict[sId].fPredictProbability:
                    if abs(oPsm.fMassDiff) < abs(psm_dict[sId].fMassDiff): 
                        psm_dict[sId] = oPsm
                    elif abs(oPsm.fMassDiff) == abs(psm_dict[sId].fMassDiff):
                        # calculate PTM scores
                        if oPsm.PTMscore < psm_dict[sId].PTMscore:
                            psm_dict[sId] = oPsm
                        elif oPsm.PTMscore == psm_dict[sId].PTMscore:
                            if oPsm.IdentifiedPeptide.upper() < psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId] = oPsm
                            elif oPsm.IdentifiedPeptide.upper() == psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId].add_protein(oPsm.protein_list)
                    
            else:
                psm_dict[sId] = oPsm

    for _key, value in psm_dict.iteritems():
        psm_new_list.append(value)
    
    return psm_new_list

## sip mode, first round fitering using score cutoff
def cutoff_filtering(psm_list, config_dict=None, fdr_given=None):
    if fdr_given == None:
        fdr_given = float(config_dict[pro_iden_str + FDR_Threshold_str])
    for oPsm in psm_list:
        oPsm.fPredictProbability = oPsm.lfScores[0]
    psm_new_list = re_rank(psm_list, False)
    psm_return_list = []
    charge_set = [[1, 1], [2, 2], [3, 10000]]
    if config_dict[pro_iden_str + FDR_Filtering_str] == 'PSM':
        for e in charge_set:
            psm_filtered_list_local = show_Fdr(psm_new_list, fdr_given, e[0], e[1])
            psm_return_list.extend(psm_filtered_list_local)
    else:
        for e in charge_set:  
            psm_filtered_list_local = show_Fdr_Pep(psm_new_list, fdr_given, e[0], e[1])
            psm_return_list.extend(psm_filtered_list_local)
            
    # psm_return_list = re_rank(psm_return_list, False)
    
    return psm_return_list

## sip mode, logistic regression for sip mode
def logistic_regression_sip(psm_list, config_dict=None):
    fdr_given = float(config_dict[pro_iden_str + FDR_Threshold_str])*(Test_Fwd_Ratio)
    # machine learning
    # # construct training data
    train_data_list = []
    train_label_list = []
    test_data_list = []
    psm_rank_list = []
    positive_int = 1
    negative_int = 0
    num_pos = 0
    num_neg = 0
    for oPsm in psm_list:
        if oPsm.RealLabel == LabelReserve:
            continue
        test_data_list.append(oPsm.feature_list)
        psm_rank_list.append(oPsm)
        if oPsm.TrainingLabel == LabelSipTrainFwd:
            train_data_list.append(oPsm.feature_list)
            train_label_list.append(positive_int)
            num_pos += 1
        elif oPsm.RealLabel == LabelTrain:
            train_data_list.append(oPsm.feature_list)
            train_label_list.append(negative_int)
            num_neg += 1
    # print(str(num_pos))
    # print(str(num_neg))
        
    train_data_np = np.array(train_data_list)[:, feature_selection_list]
    train_label_np = np.array(train_label_list)
    
    # only forward left
    unique_np = np.unique(train_label_np)
    if unique_np.shape[0] == 1:
        psm_filtered_list_local = show_Fdr(psm_list, fdr_given)
        return psm_filtered_list_local
    
    # # training
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False)
                                             # n_jobs=-1)
    logreg.fit(train_data_np, train_label_np)

    # # test
    test_unknown_np = np.array(test_data_list)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    
    for i in range(len(predict_np)):
        psm_rank_list[i].fPredictProbability = predict_np[i, 1]
   
    psm_new_list = re_rank(psm_rank_list)
    
    if config_dict[pro_iden_str + FDR_Filtering_str] == 'PSM':
        psm_filtered_list_local = show_Fdr(psm_new_list, fdr_given)
    else:
        psm_filtered_list_local = show_Fdr_Pep(psm_new_list, fdr_given)        
    return psm_filtered_list_local

## regular search mode
def logistic_regression_no_category(psm_list, config_dict=None):
    fdr_given = float(config_dict[pro_iden_str + FDR_Threshold_str])*(Test_Fwd_Ratio)
    # machine learning
    # # construct training data
    #psm_list_selected = []
    train_data_list = []
    train_label_list = []
    test_data_list = []
    psm_rank_list = []
    positive_int = 1
    negative_int = 0
    bDisableLocalRank = False
    if len(psm_list) < 800000:
        bDisableLocalRank = True
    for oPsm in psm_list:
        if oPsm.RealLabel == LabelReserve:
            continue
        test_data_list.append(oPsm.feature_list)
        psm_rank_list.append(oPsm)
        if oPsm.RealLabel != LabelTrain and (oPsm.iLocalRank == 0 or bDisableLocalRank ):
            train_data_list.append(oPsm.feature_list)
            train_label_list.append(positive_int)
        elif oPsm.RealLabel == LabelTrain:
            train_data_list.append(oPsm.feature_list)
            train_label_list.append(negative_int)
        
    train_data_np = np.array(train_data_list)[:, feature_selection_list]
    train_label_np = np.array(train_label_list)
    
     # only forward left
    unique_np = np.unique(train_label_np)
    if unique_np.shape[0] == 1:
        print('no decoy in the training set')
        sys.exit(1)
        '''
        psm_filtered_list_local = show_Fdr(psm_list, fdr_given)
        return psm_filtered_list_local
        '''
    # # training
    # num_positive = float((train_label_np==LabelPositive).sum())
    # num_negative = float((train_label_np==LabelNegative).sum())
    # num_positive = 100.0
    # num_negative = 1.0
    # class_weight_dict = {0: (num_positive/(num_negative+num_positive)), 1:(num_negative/(num_negative+num_positive))}
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False)
                                             # n_jobs=-1)
    logreg.fit(train_data_np, train_label_np)

    # # test
    test_unknown_np = np.array(test_data_list)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    
    for i in range(len(predict_np)):
        psm_rank_list[i].fPredictProbability = predict_np[i, 1]
   
    psm_new_list = re_rank(psm_rank_list)
    
    if config_dict[pro_iden_str + FDR_Filtering_str] == 'PSM':
        psm_filtered_list_local = show_Fdr(psm_new_list, fdr_given)
    else:
        psm_filtered_list_local = show_Fdr_Pep(psm_new_list, fdr_given)        
    return psm_filtered_list_local

## get the number of missed cleavage sites
def get_num_missed_cleavage_sites(sIdentifiedSeq, sResiduesBeforeCleavage, sResiduesAfterCleavage):
    count_int = 0
    for i in range(len(sIdentifiedSeq) - 1):
        if sIdentifiedSeq[i] in sResiduesBeforeCleavage and sIdentifiedSeq[i + 1] in sResiduesAfterCleavage:
            count_int += 1
    return count_int

## sip mode, if the percentage difference larger than 10, do not use it for the count
def generate_Prophet_features_sip_diff(lPsm, config_dict):
    # get some statistics
    global num_forward_psms_before_filtering
    num_forward_psms_before_filtering = 0
    for one_psm in lPsm:
        if one_psm.RealLabel == LabelFwd:
            num_forward_psms_before_filtering += 1
    simple_feature_bool = False
    if float(num_forward_psms_before_filtering)/float(len(lPsm)) > 0.7 and len(lPsm) > 1500000:
        simple_feature_bool = True
    # peptide with PTM dictionary is for IPSC, identified peptide
    peptide_with_modification_dict = {}
    # peptide without PTM dictionary is for OPSC, original peptide
    peptide_dict = {}
    peptide_protein_dict = {}
    # psm_set = Set()
    for oPsm in lPsm:
        oPsm.NMC = get_num_missed_cleavage_sites(oPsm.OriginalPeptide, 
                                                 config_dict[pep_iden_str + cleave_after_residues_str],
                                                 config_dict[pep_iden_str + cleave_before_residues_str])
           
        if oPsm.IdentifiedPeptide in peptide_with_modification_dict:
            peptide_with_modification_dict[oPsm.IdentifiedPeptide].append(oPsm.pct_s)
        else:
            peptide_with_modification_dict[oPsm.IdentifiedPeptide] = [oPsm.pct_s]
        if oPsm.OriginalPeptide in peptide_protein_dict:
            pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
            for protein in oPsm.protein_list:
                if protein not in pro_list:
                    pro_list.append(protein)
        else:
            pro_list = []
            for pro in oPsm.protein_list:
                if pro not in pro_list:
                    pro_list.append(pro)
            peptide_protein_dict[oPsm.OriginalPeptide] = pro_list
            

    pattern = re.compile('[^\w\[\]]')
    for key, _value in peptide_with_modification_dict.iteritems():
        peptide_str = pattern.sub('', key)
        if peptide_str in peptide_dict:
            # peptide_dict[peptide_str] += 1
            peptide_dict[peptide_str] += _value
        else:
            # peptide_dict[peptide_str] = 1
            peptide_dict[peptide_str] = _value
    
    # # sibling peptides
    pro_unique_dict = {} # number of unique peptide of a protain
    pro_shared_dict = {} # number of shared peptide of a protain
    # debug
    pro_balanced_shared_dict = {}
    # debug
    num_changed = 0
    changed_flag = False
    # psm_set.clear()
    pro_pep_dict = {}
    pro_unique_pep_dict = {}
    pro_shared_pep_dict = {}
    for oPsm in lPsm:
        pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
        changed_flag = False
        print_flag = True
        for protein in pro_list:
            if not protein in oPsm.protein_list:
                changed_flag = True
                if not print_flag:
                    print(pro_list)
                    print(oPsm.protein_list)
                    print_flag = True
                oPsm.protein_list.append(protein)
        if len(oPsm.protein_list) != len(pro_list):
            print('check 4')
            print(oPsm.protein_list)
            print(pro_list)
            exit(1)
            
        if changed_flag:
            oPsm.set_protein_names()
            oPsm.RealLabel = get_protein_type(oPsm.ProteinNames)
            # print(oPsm.OriginalPeptide)
            num_changed += 1
            
        if len(pro_list) > 1:
            num_pro_float = float(len(pro_list))
            for protein in pro_list:
                if protein in pro_shared_dict: 
                    pro_shared_dict[protein] += 1
                else:
                    pro_shared_dict[protein] = 1
                # debug
                if protein in pro_balanced_shared_dict:
                    pro_balanced_shared_dict[protein] += 1.0/num_pro_float
                else:
                    pro_balanced_shared_dict[protein] = 1.0/num_pro_float
                # debug
        else:
            if pro_list[0] in pro_unique_dict:
                pro_unique_dict[pro_list[0]] += 1
            else:
                pro_unique_dict[pro_list[0]] = 1
        
        for pro in pro_list:
            if pro in pro_pep_dict:
                l = pro_pep_dict[pro]
                if oPsm.OriginalPeptide not in l:
                    l.append(oPsm.OriginalPeptide)
            else:
                l = []
                l.append(oPsm.OriginalPeptide)
                pro_pep_dict[pro] = l
        if len(pro_list) == 1:
            if pro_list[0] in pro_unique_pep_dict:
                l = pro_unique_pep_dict[pro_list[0]]
                if oPsm.OriginalPeptide not in l:
                    l.append(oPsm.OriginalPeptide)
            else:
                l = []
                l.append(oPsm.OriginalPeptide)
                pro_unique_pep_dict[pro_list[0]] = l
        else:
            for pro in pro_list:
                if pro in pro_shared_pep_dict:
                    l = pro_shared_pep_dict[pro]
                    if oPsm.OriginalPeptide not in l:
                        l.append(oPsm.OriginalPeptide)
                else:
                    l = []
                    l.append(oPsm.OriginalPeptide)
                    pro_shared_pep_dict[pro] = l
    if num_changed != 0:
        if not print_flag:
            print("num changed %d" % num_changed)
    
    # collect features
    num_unique_per_pro = 0
    num_shared_per_pro = 0
    num_balanced_shared_per_pro = 0
    
    max_unique_per_psm = 0
    max_shared_per_psm = 0
    max_balanced_shared_per_psm = 0
    
    num_unique_per_psm = 0
    num_shared_per_psm = 0
    num_balanced_shared_per_psm = 0
    
    max_per_psm = 0 # linked to protein
    max_balanced_per_psm = 0 # linked to a protein
    
    max_linked_unique_per_psm = 0
    max_linked_shared_per_psm = 0
    max_linked_balanced_unique_per_psm = 0
    max_linked_balanced_shared_per_psm = 0
    
    for oPsm in lPsm:
        oPsm.IPSC = peptide_with_modification_dict[oPsm.IdentifiedPeptide]
        oPsm.OPSC = peptide_dict[oPsm.OriginalPeptide]
        
        max_unique_per_psm = 0
        max_shared_per_psm = 0
        max_balanced_shared_per_psm = 0
    
        num_unique_per_psm = 0
        num_shared_per_psm = 0
        num_balanced_shared_per_psm = 0
        
        max_per_psm = 0 # linked to protein
        max_balanced_per_psm = 0 # linked to a protein
    
        max_linked_unique_per_psm = 0
        max_linked_shared_per_psm = 0
        max_linked_balanced_unique_per_psm = 0
        max_linked_balanced_shared_per_psm = 0
        
        pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
        
        for protein in pro_list:
            if len(pro_pep_dict[protein]) > oPsm.PPC:
                oPsm.PPC = len(pro_pep_dict[protein])
            if len(pro_list) == 1 and len(pro_unique_pep_dict[protein]) > oPsm.UPPC:
                oPsm.UPPC = len(pro_unique_pep_dict[protein])
            if len(pro_list) > 1 and len(pro_shared_pep_dict[protein]) > oPsm.SPPC:
                oPsm.SPPC = len(pro_shared_pep_dict[protein])
            num_unique_per_pro = 0
            num_shared_per_pro = 0
            num_balanced_shared_per_pro = 0
            
            if protein in pro_unique_dict:
                num_unique_per_pro = pro_unique_dict[protein]
                if len(pro_list) == 1:
                    # pass
                    num_unique_per_pro -= 1
                
            if protein in pro_shared_dict:
                num_shared_per_pro = pro_shared_dict[protein]
                if len(pro_list) > 1:
                    # pass
                    num_shared_per_pro -= 1
            
            if protein in pro_balanced_shared_dict:
                num_balanced_shared_per_pro = pro_balanced_shared_dict[protein]
                if len(pro_list) > 1:
                    num_balanced_shared_per_pro -= 1.0/ float(len(pro_list))
            
            if num_unique_per_pro > max_unique_per_psm:
                max_unique_per_psm = num_unique_per_pro
            if num_shared_per_pro > max_shared_per_psm:
                max_shared_per_psm = num_shared_per_pro
            if num_unique_per_pro + num_shared_per_pro > max_per_psm:
                max_per_psm = num_unique_per_pro + num_shared_per_pro
                max_linked_unique_per_psm = num_unique_per_pro
                max_linked_shared_per_psm = num_shared_per_pro
            num_unique_per_psm += num_unique_per_pro
            num_shared_per_psm += num_shared_per_pro
            
            if num_balanced_shared_per_pro > max_balanced_shared_per_psm:
                max_balanced_shared_per_psm = num_balanced_shared_per_pro
            if num_unique_per_pro + num_balanced_shared_per_pro > max_balanced_per_psm:
                max_balanced_per_psm = num_unique_per_pro + num_balanced_shared_per_pro
                max_linked_balanced_unique_per_psm = num_unique_per_pro
                max_linked_balanced_shared_per_psm = num_balanced_shared_per_pro
            num_balanced_shared_per_psm += num_balanced_shared_per_pro
        
        oPsm.UPSC = num_unique_per_psm
        oPsm.SPSC = num_shared_per_psm
        oPsm.SPSC = num_unique_per_psm + num_shared_per_psm
        oPsm.UPSC = max_unique_per_psm
        oPsm.SPSC = max_shared_per_psm
        if len(oPsm.protein_list) == 1:
            oPsm.UPSC = 1
        else:
            oPsm.UPSC = 0
        if simple_feature_bool:
            oPsm.SPSC = max_linked_unique_per_psm
        else:
            oPsm.SPSC = max_linked_unique_per_psm + max_linked_shared_per_psm
        
## collect spectrum and protein level features
def generate_Prophet_features_test(lPsm, config_dict):
    # get some statistics
    global num_forward_psms_before_filtering
    num_forward_psms_before_filtering = 0
    for one_psm in lPsm:
        if one_psm.RealLabel == LabelFwd:
            num_forward_psms_before_filtering += 1
    simple_feature_bool = False
    if float(num_forward_psms_before_filtering)/float(len(lPsm)) > 0.7 and len(lPsm) > 1500000:
        simple_feature_bool = True
    if config_dict[pep_iden_str + search_type_str] == 'SIP':
        simple_feature_bool = False
    # peptide with PTM dictionary is for IPSC
    peptide_with_modification_dict = {}
    # peptide without PTM dictionary is for OPSC
    peptide_dict = {}
    peptide_protein_dict = {}
    # psm_set = Set()
    for oPsm in lPsm:
        oPsm.NMC = get_num_missed_cleavage_sites(oPsm.OriginalPeptide, 
                                                 config_dict[pep_iden_str + cleave_after_residues_str],
                                                 config_dict[pep_iden_str + cleave_before_residues_str])
           
        if oPsm.IdentifiedPeptide in peptide_with_modification_dict:
            peptide_with_modification_dict[oPsm.IdentifiedPeptide] += 1
        else:
            peptide_with_modification_dict[oPsm.IdentifiedPeptide] = 1
        if oPsm.OriginalPeptide in peptide_protein_dict:
            pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
            for protein in oPsm.protein_list:
                if protein not in pro_list:
                    pro_list.append(protein)
        else:
            pro_list = []
            for pro in oPsm.protein_list:
                if pro not in pro_list:
                    pro_list.append(pro)
            peptide_protein_dict[oPsm.OriginalPeptide] = pro_list
            

    pattern = re.compile('[^\w\[\]]')
    for key, _value in peptide_with_modification_dict.iteritems():
        peptide_str = pattern.sub('', key)
        if peptide_str in peptide_dict:
            # peptide_dict[peptide_str] += 1
            peptide_dict[peptide_str] += _value
        else:
            # peptide_dict[peptide_str] = 1
            peptide_dict[peptide_str] = _value
    
    # # sibling peptides
    pro_unique_dict = {} # number of unique peptide of a protain
    pro_shared_dict = {} # number of shared peptide of a protain
    # debug
    pro_balanced_shared_dict = {}
    # debug
    num_changed = 0
    changed_flag = False
    # psm_set.clear()
    pro_pep_dict = {}
    pro_unique_pep_dict = {}
    pro_shared_pep_dict = {}
    for oPsm in lPsm:
        pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
        '''
        if len(oPsm.protein_list) != len(pro_list):
            print('check 3')
            print(oPsm.protein_list)
            print(pro_list)
        '''    
        changed_flag = False
        print_flag = True
        for protein in pro_list:
            if not protein in oPsm.protein_list:
                changed_flag = True
                if not print_flag:
                    print(pro_list)
                    print(oPsm.protein_list)
                    print_flag = True
                oPsm.protein_list.append(protein)
        if len(oPsm.protein_list) != len(pro_list):
            print('check 4')
            print(oPsm.protein_list)
            print(pro_list)
            exit(1)
            
        if changed_flag:
            oPsm.set_protein_names()
            oPsm.RealLabel = get_protein_type(oPsm.ProteinNames)
            # print(oPsm.OriginalPeptide)
            num_changed += 1
        '''
        unique_id_str = oPsm.FileName + '_' + str(oPsm.ScanNumber) + '_' + oPsm.IdentifiedPeptide
        if unique_id_str in psm_set:
            continue
        else:
            psm_set.add(unique_id_str)
        '''
        if len(pro_list) > 1:
            num_pro_float = float(len(pro_list))
            for protein in pro_list:
                if protein in pro_shared_dict: 
                    pro_shared_dict[protein] += 1
                else:
                    pro_shared_dict[protein] = 1
                # debug
                if protein in pro_balanced_shared_dict:
                    pro_balanced_shared_dict[protein] += 1.0/num_pro_float
                else:
                    pro_balanced_shared_dict[protein] = 1.0/num_pro_float
                # debug
        else:
            if pro_list[0] in pro_unique_dict:
                pro_unique_dict[pro_list[0]] += 1
            else:
                pro_unique_dict[pro_list[0]] = 1
        
        for pro in pro_list:
            if pro in pro_pep_dict:
                l = pro_pep_dict[pro]
                if oPsm.OriginalPeptide not in l:
                    l.append(oPsm.OriginalPeptide)
            else:
                l = []
                l.append(oPsm.OriginalPeptide)
                pro_pep_dict[pro] = l
        if len(pro_list) == 1:
            if pro_list[0] in pro_unique_pep_dict:
                l = pro_unique_pep_dict[pro_list[0]]
                if oPsm.OriginalPeptide not in l:
                    l.append(oPsm.OriginalPeptide)
            else:
                l = []
                l.append(oPsm.OriginalPeptide)
                pro_unique_pep_dict[pro_list[0]] = l
        else:
            for pro in pro_list:
                if pro in pro_shared_pep_dict:
                    l = pro_shared_pep_dict[pro]
                    if oPsm.OriginalPeptide not in l:
                        l.append(oPsm.OriginalPeptide)
                else:
                    l = []
                    l.append(oPsm.OriginalPeptide)
                    pro_shared_pep_dict[pro] = l
    if num_changed != 0:
        if not print_flag:
            print("num changed %d" % num_changed)
    
    # collect features
    num_unique_per_pro = 0
    num_shared_per_pro = 0
    num_balanced_shared_per_pro = 0
    
    max_unique_per_psm = 0
    max_shared_per_psm = 0
    max_balanced_shared_per_psm = 0
    
    num_unique_per_psm = 0
    num_shared_per_psm = 0
    num_balanced_shared_per_psm = 0
    
    max_per_psm = 0 # linked to protein
    max_balanced_per_psm = 0 # linked to a protein
    
    max_linked_unique_per_psm = 0
    max_linked_shared_per_psm = 0
    max_linked_balanced_unique_per_psm = 0
    max_linked_balanced_shared_per_psm = 0
    
    for oPsm in lPsm:
        oPsm.IPSC = peptide_with_modification_dict[oPsm.IdentifiedPeptide]
        oPsm.OPSC = peptide_dict[oPsm.OriginalPeptide]
        
        max_unique_per_psm = 0
        max_shared_per_psm = 0
        max_balanced_shared_per_psm = 0
    
        num_unique_per_psm = 0
        num_shared_per_psm = 0
        num_balanced_shared_per_psm = 0
        
        max_per_psm = 0 # linked to protein
        max_balanced_per_psm = 0 # linked to a protein
    
        max_linked_unique_per_psm = 0
        max_linked_shared_per_psm = 0
        max_linked_balanced_unique_per_psm = 0
        max_linked_balanced_shared_per_psm = 0
        
        pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
        
        for protein in pro_list:
            if len(pro_pep_dict[protein]) > oPsm.PPC:
                oPsm.PPC = len(pro_pep_dict[protein])
            if len(pro_list) == 1 and len(pro_unique_pep_dict[protein]) > oPsm.UPPC:
                oPsm.UPPC = len(pro_unique_pep_dict[protein])
            if len(pro_list) > 1 and len(pro_shared_pep_dict[protein]) > oPsm.SPPC:
                oPsm.SPPC = len(pro_shared_pep_dict[protein])
            num_unique_per_pro = 0
            num_shared_per_pro = 0
            num_balanced_shared_per_pro = 0
            
            if protein in pro_unique_dict:
                num_unique_per_pro = pro_unique_dict[protein]
                if len(pro_list) == 1:
                    # pass
                    num_unique_per_pro -= 1
                
            if protein in pro_shared_dict:
                num_shared_per_pro = pro_shared_dict[protein]
                if len(pro_list) > 1:
                    # pass
                    num_shared_per_pro -= 1
            
            if protein in pro_balanced_shared_dict:
                num_balanced_shared_per_pro = pro_balanced_shared_dict[protein]
                if len(pro_list) > 1:
                    num_balanced_shared_per_pro -= 1.0/ float(len(pro_list))
            
            if num_unique_per_pro > max_unique_per_psm:
                max_unique_per_psm = num_unique_per_pro
            if num_shared_per_pro > max_shared_per_psm:
                max_shared_per_psm = num_shared_per_pro
            if num_unique_per_pro + num_shared_per_pro > max_per_psm:
                max_per_psm = num_unique_per_pro + num_shared_per_pro
                max_linked_unique_per_psm = num_unique_per_pro
                max_linked_shared_per_psm = num_shared_per_pro
            num_unique_per_psm += num_unique_per_pro
            num_shared_per_psm += num_shared_per_pro
            
            if num_balanced_shared_per_pro > max_balanced_shared_per_psm:
                max_balanced_shared_per_psm = num_balanced_shared_per_pro
            if num_unique_per_pro + num_balanced_shared_per_pro > max_balanced_per_psm:
                max_balanced_per_psm = num_unique_per_pro + num_balanced_shared_per_pro
                max_linked_balanced_unique_per_psm = num_unique_per_pro
                max_linked_balanced_shared_per_psm = num_balanced_shared_per_pro
            num_balanced_shared_per_psm += num_balanced_shared_per_pro
        
        oPsm.UPSC = num_unique_per_psm
        oPsm.SPSC = num_shared_per_psm
        oPsm.SPSC = num_unique_per_psm + num_shared_per_psm
        oPsm.UPSC = max_unique_per_psm
        oPsm.SPSC = max_shared_per_psm
        # oPsm.SPSC = max_shared_per_psm + max_unique_per_psm
        # oPsm.PPC = oPsm.UPPC + oPsm.SPPC
        # oPsm.UPSC = max_linked_unique_per_psm
        # oPsm.SPSC = max_linked_shared_per_psm
        # oPsm.SPSC = max_linked_balanced_shared_per_psm
        if len(oPsm.protein_list) == 1:
            oPsm.UPSC = 1
        else:
            oPsm.UPSC = 0
        if simple_feature_bool:
            oPsm.SPSC = max_linked_unique_per_psm
        else:
            oPsm.SPSC = max_linked_unique_per_psm + max_linked_shared_per_psm
        # debug
        # oPsm.SPSC = max_linked_unique_per_psm
        # debug
        # oPsm.SPSC = float(max_linked_unique_per_psm)/float(len(oPsm.protein_list)) + float(max_linked_shared_per_psm)/float(len(oPsm.protein_list))
        # oPsm.SPSC = max_per_psm
        
        # oPsm.UPSC = max_unique_per_psm
        # oPsm.SPSC = max_shared_per_psm
        
## Exit system with error message
def die(msg=None):
    if msg is not None:
        print >> sys.stderr, msg
        sys.exit(1)

## Check file exist
def check_file_exist(filename):

    try:
        with open(filename) as _f: pass
    except IOError as _e:
        print >> sys.stderr, '\nCannot open', filename
        die("Program exit!")

## parameters used to read configuration file
pep_iden_str = '[Peptide_Identification]'
search_type_str = 'Search_Type'
fasta_database_str = 'FASTA_Database'
pro_iden_str = '[Protein_Identification]'
FDR_Threshold_str = 'FDR_Threshold'
training_decoy_prefix_str = 'Training_Decoy_Prefix'
testing_decoy_prefix_str = 'Testing_Decoy_Prefix'
reserved_decoy_prefix_str = 'Reserved_Decoy_Prefix'
FDR_Filtering_str = 'FDR_Filtering'
cleave_after_residues_str = 'Cleave_After_Residues'
cleave_before_residues_str = 'Cleave_Before_Residues'
Mass_Tolerance_Parent_Ion_str = 'Filter_Mass_Tolerance_Parent_Ion'
search_Mass_Tolerance_Parent_Ion_str = 'Search_Mass_Tolerance_Parent_Ion'
Parent_Mass_Windows_str = 'Parent_Mass_Windows'
Filter_Mass_Tolerance_Parent_Ion_Unit = 'Filter_Mass_Tolerance_Parent_Ion_Unit'

## Parse config file
def parse_config(config_filename):

    # Save config values to dictionary
    config_dict = {}    # initialize dictionary

    # Call Yinfeng's parseconfig.py module
    check_file_exist(config_filename)
    
    # Save all config values to dictionary
    all_config_dict = parseconfig.parseConfigKeyValues(config_filename)
    
    # get the mass window
    line_str = all_config_dict[pep_iden_str + Parent_Mass_Windows_str]
    e_list = line_str.strip().split(',')
    global mass_window_max_int
    mass_window_max_int = 0
    for e in e_list:
        v = abs(int(e.strip()))
        if v > mass_window_max_int:
            mass_window_max_int = v
    line_str = all_config_dict[pep_iden_str + search_Mass_Tolerance_Parent_Ion_str]
    if float(line_str.strip()) > 0.5:
        mass_window_max_int += 1
    
    global train_str, test_str, reserve_str
    
    # do not care training or testing
    if all_config_dict[pep_iden_str + search_type_str] == 'SIP':
        # test_str = os.path.commonprefix([all_config_dict[pro_iden_str + training_decoy_prefix_str], all_config_dict[pro_iden_str + testing_decoy_prefix_str]])
        test_str = all_config_dict[pro_iden_str + testing_decoy_prefix_str]
        train_str = test_str + '_fake'
        reserve_str = ''
        return all_config_dict
    
    train_str = all_config_dict[pro_iden_str + training_decoy_prefix_str]
    test_str = all_config_dict[pro_iden_str + testing_decoy_prefix_str]
    if pro_iden_str + reserved_decoy_prefix_str in all_config_dict:
        reserve_str = all_config_dict[pro_iden_str + reserved_decoy_prefix_str]
    else:
        reserve_str = ''
    
    # return config dictionary
    return all_config_dict

## reset the configurations
def config_reset(all_config_dict):
    global train_str, test_str, reserve_str
    
    train_str = all_config_dict[pro_iden_str + training_decoy_prefix_str]
    test_str = all_config_dict[pro_iden_str + testing_decoy_prefix_str]
    if pro_iden_str + reserved_decoy_prefix_str in all_config_dict:
        reserve_str = all_config_dict[pro_iden_str + reserved_decoy_prefix_str]
    else:
        reserve_str = ''

## peptide class for generating pep.txt file
class Peptide:
    
    def __init__(self):
        self.IdentifiedPeptide = ''
        self.ParentCharge = ''
        self.OriginalPeptide = ''
        self.ProteinNames = []
        self.ProteinCount = 0
        self.SpectralCount = 0
        self.BestScore = 0.0
        self.PSMs = []
        self.ScanType = []
        self.SearchName = []
        
    def add(self, oPsm):
        self.SpectralCount += 1
        if self.BestScore < oPsm.fPredictProbability:
            self.BestScore = oPsm.fPredictProbability
        self.PSMs.append(filename_list[oPsm.FileName] + '['+str(oPsm.ScanNumber) + ']')
        self.ScanType.append(scantype_list[oPsm.ScanType])
        self.SearchName.append(searchname_list[oPsm.SearchName])
        if oPsm.RealLabel == LabelFwd:
            self.TargetMatch = 'T'
        
    def set(self, oPsm):
        self.IdentifiedPeptide = oPsm.IdentifiedPeptide
        self.ParentCharge = oPsm.ParentCharge
        self.OriginalPeptide = oPsm.OriginalPeptide
        self.ProteinNames = oPsm.ProteinNames
        self.ProteinCount = len(oPsm.protein_list)
        self.SpectralCount = 1
        self.BestScore = oPsm.fPredictProbability
        self.PSMs.append(filename_list[oPsm.FileName] + '['+str(oPsm.ScanNumber) + ']')
        self.ScanType.append(scantype_list[oPsm.ScanType])
        self.SearchName.append(searchname_list[oPsm.SearchName])
        if oPsm.RealLabel == LabelFwd:
            self.TargetMatch = 'T'
        else:
            self.TargetMatch = 'F'
        
    def __repr__(self):
        l = [self.IdentifiedPeptide,
             str(self.ParentCharge),
             self.OriginalPeptide,
             self.ProteinNames,
             str(self.ProteinCount),
             self.TargetMatch,
             str(self.SpectralCount),
             str(self.BestScore),
             ('{'+','.join(self.PSMs)+'}'),
             ('{'+','.join(self.ScanType)+'}'),
             ('{'+','.join(self.SearchName)+'}')]
        
        return '\t'.join(l) 

## generate psm.txt and pep.txt files
def generate_psm_pep_txt(base_out, out_folder, psm_filtered_list, config_dict):
    
    # protein_identification message
    pro_iden_msg = ""
    pro_iden_msg += "#\t########################\n"
    pro_iden_msg += "#\t# Filtering Parameters #\n"
    pro_iden_msg += "#\t########################\n"
    pro_iden_msg += "#\t\n"
    pro_iden_msg += "#\t[Protein_Identification]\n"
    pro_iden_msg += "#\t\n"
    pro_iden_msg += "#\t# the prefix of training decoy sequences' locus IDs in the database\n"
    pro_iden_msg += "#\t" + training_decoy_prefix_str + " = " + str(train_str) + "\n"
    pro_iden_msg += "#\t\n"
    pro_iden_msg += "#\t# the prefix of test decoy sequences' locus IDs in the database\n"
    pro_iden_msg += "#\t" + testing_decoy_prefix_str + " = " + str(test_str) + "\n"
    pro_iden_msg += "#\t\n"
    pro_iden_msg += "#\t# Level of FDR filtering. Options: \"PSM\" and \"Peptide\"\n"
    pro_iden_msg += "#\t" + FDR_Filtering_str + " = " + config_dict[pro_iden_str + FDR_Filtering_str] + "\n"
    pro_iden_msg += "#\t\n"
    pro_iden_msg += "#\t# FDR threshold for filtering peptide identifications\n"
    pro_iden_msg += "#\t" + FDR_Threshold_str + " = " + config_dict[pro_iden_str + FDR_Threshold_str] + "\n"
    pro_iden_msg += "#\t\n"

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
    psm_column_msg += "#\tScore = Predicted Probability of being true PSM\n"
    psm_column_msg += "#\tDeltaZ = Difference between the best PSM score and the next best PSM of this scan\n"
    psm_column_msg += "#\tDeltaP = Difference between the best modified PSM and its PTM isoform\n"
    psm_column_msg += "#\tIdentifiedPeptide = Identified peptide sequence with potential PTMs and mutations\n"
    psm_column_msg += "#\tOriginalPeptide = Original peptide sequence in the FASTA file\n"
    psm_column_msg += "#\tProteinNames = Names of proteins of the peptide\n"
    psm_column_msg += "#\tProteinCount = Number of proteins that the peptide can be assigned to\n"
    psm_column_msg += "#\tTargetMatch = T for target match and F for decoy match\n"
    psm_column_msg += "#\t\n"

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

    # get the FDR, # target, # decoy
    psm_target_int = 0
    psm_decoy_int = 0
    pep_target_int = 0
    pep_decoy_int = 0
    pep_set = Set()
    for oPsm in psm_filtered_list:
        if oPsm.RealLabel == LabelTrain or oPsm.RealLabel == LabelReserve:
            continue
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if pep_str not in pep_set:
            if oPsm.RealLabel == LabelFwd:
                pep_target_int += 1
            else:
                pep_decoy_int += 1
            pep_set.add(pep_str)
        
        if oPsm.RealLabel == LabelFwd:
            psm_target_int += 1
        else:
            psm_decoy_int += 1
    
    # write out into files
    base_out = os.path.join(out_folder, base_out)
    psm_txt_file_str = base_out + ".psm.txt"
    with open(psm_txt_file_str, 'w') as fw:
        fw.write(pro_iden_msg)
        fw.write('#\t########################################\n')
        fw.write('#\t####### PSM Filtering by Sipros ########\n')
        fw.write('#\t########################################\n')
        fw.write('#\t\n')
        fw.write('#\t#######################\n')
        fw.write('#\t# Statistical Results #\n')
        fw.write('#\t#######################\n')
        fw.write('#\t\n')
        fw.write('#\t[Statistical_Results]\n')
        fw.write('#\t\n')
        fw.write('#\t# Numbers of psm after filtering\n')
        fw.write('#\tDecoy_PSMs_After_Filtering = %d\n' % psm_decoy_int)
        fw.write('#\tTarget_PSMs_After_Filtering = %d\n' % psm_target_int)
        fw.write('#\t# PSM FDR = Decoy_PSMs_After_Filtering / Target_PSMs_After_Filtering\n')
        fw.write('#\tPSM_FDR = %.2f%%\n' % ((1/Test_Fwd_Ratio) * 100.0 * psm_decoy_int/psm_target_int))
        fw.write('#\t\n')
        fw.write(psm_column_msg)
        # for psm out
        psm_out_list = ['Filename',  # 0
                    'ScanNumber',  # 1
                    'ParentCharge',  # 2
                    'MeasuredParentMass',  # 3
                    'CalculatedParentMass',  # 4
                    'MassErrorDa',  # 5 CalculatedParentMass - MeasuredParentMass
                    'MassErrorPPM',  # 6 MassErrorDa / CalculatedParentMass
                    'ScanType',  # 7
                    'SearchName',  # 8
                    'ScoringFunction',  # 9
                    'Score',  # 10
                    'DeltaZ',  # 11 the difference score between the rank 1 and 2
                    'DeltaP',  # 12
                    'IdentifiedPeptide',  # 13
                    'OriginalPeptide',  # 14
                    'ProteinNames',  # 15
                    'ProteinCount',  # 16
                    'TargetMatch']  # 17
        fw.write('\t'.join(psm_out_list) + '\n')
        for oPsm in psm_filtered_list:
            if oPsm.RealLabel == LabelTrain:
                continue 
            oPsm.clean_protein_name()
            fw.write(filename_list[oPsm.FileName])
            fw.write('\t')
            fw.write(str(oPsm.ScanNumber))
            fw.write('\t')
            fw.write(str(oPsm.ParentCharge))
            fw.write('\t')
            fw.write('%.3f' % oPsm.MeasuredParentMass)
            fw.write('\t')
            fw.write('%.3f' % oPsm.CalculatedParentMass)
            fw.write('\t')
            fw.write('%.3f' % (oPsm.fMassDiff))
            fw.write('\t')
            fw.write('%.3f' % (1000000*(oPsm.fMassDiff)/oPsm.CalculatedParentMass))
            fw.write('\t')
            fw.write(scantype_list[oPsm.ScanType])
            fw.write('\t')
            fw.write(searchname_list[oPsm.SearchName])
            fw.write('\t')
            fw.write('SiprosEnsemble')
            fw.write('\t')
            fw.write(str(oPsm.fPredictProbability))
            fw.write('\t')
            fw.write('NA')
            fw.write('\t')
            fw.write(oPsm.DeltaP)
            fw.write('\t')
            fw.write(oPsm.IdentifiedPeptide)
            fw.write('\t')
            fw.write(oPsm.OriginalPeptide)
            fw.write('\t')
            fw.write(oPsm.ProteinNames)
            fw.write('\t')
            fw.write(str(len(oPsm.protein_list)))
            fw.write('\t')
            if oPsm.RealLabel == LabelFwd:
                fw.write('T')
            else:
                fw.write('F')
            fw.write('\n')
            
    # pep_sub_dict for preparing pep_out
    pep_sub_dict = {}    # initialize dict of list
    for oPsm in psm_filtered_list:
        if oPsm.RealLabel == LabelTrain or oPsm.RealLabel == LabelReserve:
            continue
        pep_ID = oPsm.IdentifiedPeptide + '_+_' + str(oPsm.ParentCharge)
        if pep_ID in pep_sub_dict:
            pep_sub_dict[pep_ID].add(oPsm)
        else:
            oPeptide = Peptide()
            oPeptide.set(oPsm)
            pep_sub_dict[pep_ID] = oPeptide
    pep_txt_file_str = base_out + ".pep.txt"
    with open(pep_txt_file_str, 'w') as fw:
        fw.write(pro_iden_msg)
        # statistic results
        fw.write('#\t########################################\n')
        fw.write('#\t####### PSM Filtering by Sipros ########\n')
        fw.write('#\t########################################\n')
        fw.write('#\t\n')
        fw.write('#\t#######################\n')
        fw.write('#\t# Statistical Results #\n')
        fw.write('#\t#######################\n')
        fw.write('#\t\n')
        fw.write('#\t[Statistical_Results]\n')
        fw.write('#\t\n')
        fw.write('#\t# Numbers of peptide after filtering\n')
        fw.write('#\tDecoy_peptides_After_Filtering = %d\n' % pep_decoy_int)
        fw.write('#\tTarget_peptides_After_Filtering = %d\n' % pep_target_int)
        fw.write('#\t# peptide FDR = Decoy_peptides_After_Filtering / Target_peptides_After_Filtering\n')
        fw.write('#\tPeptide_FDR = %.2f%%\n' % ((1/Test_Fwd_Ratio) * 100.0*pep_decoy_int/pep_target_int))
        fw.write('#\t\n')
        fw.write(pep_column_msg)
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
        fw.write('\t'.join(pep_out_list) + '\n')
        for _pep_id, oPeptide in pep_sub_dict.iteritems():
            fw.write(repr(oPeptide))
            fw.write('\n')
    
    return (psm_txt_file_str, pep_txt_file_str)

## collect the score agreement, save in iLocalRank
def remark_concensus(psm_list):
    psm_dict = {}
    for oPsm in psm_list:
        unique_id_str = str(oPsm.FileName) + '_' + str(oPsm.ScanNumber) + '_' + oPsm.IdentifiedPeptide
        if unique_id_str in psm_dict:
            psm_dict[unique_id_str] += 1
        else:
            psm_dict[unique_id_str] = 1
    
    psm_set = Set()
    for oPsm in psm_list:
        unique_id_str = str(oPsm.FileName) + '_' + str(oPsm.ScanNumber) + '_' + oPsm.IdentifiedPeptide
        count_int = psm_dict[unique_id_str]
        if count_int >= 3:
            oPsm.iLocalRank = 0
        else:
            oPsm.iLocalRank = 3 - count_int
        if count_int == 2:
            psm_set.add(str(oPsm.FileName) + '_' + str(oPsm.ScanNumber))
    
    for oPsm in psm_list:
        unique_id_str = str(oPsm.FileName) + '_' + str(oPsm.ScanNumber)
        if unique_id_str in psm_set:
            if oPsm.iLocalRank == 2:
                oPsm.iLocalRank = 3 # Mi

## mass filtering
def mass_filter(psm_list, config_dict):
    mass_tolerance = float(config_dict[pro_iden_str + Mass_Tolerance_Parent_Ion_str])
    da_filter = True
    if pro_iden_str + Filter_Mass_Tolerance_Parent_Ion_Unit in config_dict:
        if config_dict[pro_iden_str + Filter_Mass_Tolerance_Parent_Ion_Unit] == "PPM":
            da_filter = False
    psm_new_list = []
    if da_filter:
        for oPsm in psm_list:
            if abs(oPsm.fMassDiff) > mass_tolerance:
                continue
            psm_new_list.append(oPsm)
    else:
        for oPsm in psm_list:
            ppm = abs(1000000*(oPsm.MeasuredParentMass - oPsm.CalculatedParentMass)/oPsm.CalculatedParentMass)
            if ppm > mass_tolerance:
                continue
            psm_new_list.append(oPsm)
    return psm_new_list

# find the testing decoy verse forward 
def find_train_test_ratio(config_dict):
    database_str = config_dict['[Peptide_Identification]FASTA_Database']
    num_train_int = 0
    num_test_int = 0
    num_reserved_int = 0
    num_fwd_int = 0
    less_train_str = '>' + train_str
    less_test_str = '>' + test_str
    less_reserve_str = '>' + reserve_str
    with open(database_str, 'r') as f:
        for line_str in f:
            if line_str.startswith('>'):
                if line_str.startswith(less_train_str):
                    num_train_int += 1
                elif line_str.startswith(less_test_str):
                    num_test_int += 1
                elif reserve_str!= '' and line_str.startswith(less_reserve_str):
                    num_reserved_int += 1
                else:
                    num_fwd_int += 1
    
    global Test_Fwd_Ratio
    Test_Fwd_Ratio = float(num_test_int) / float(num_fwd_int)
    global num_proteins
    num_proteins = num_fwd_int
    return Test_Fwd_Ratio

## debug mode, generate pin format for percolator
def generatePINPercolator(psm_list, output_str):
    with open(output_str, 'w') as fw:
        fw.write("SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tPEP\tPRO\tMVH\tdeltMVH\tXcorr\tdeltXcorr\tWDP\tDdeltWDP\tPepLen\tdM\tabsdM\tCharge1\tCharge2\tCharge3\tCharge4\tCharge5\tenzInt\tPeptide\tProtein")
        fw.write('\n')
        row_list = []
        for psm in psm_list:
            del row_list[:]
            fileID = filename_list[psm.FileName][filename_list[psm.FileName].rfind('_')+1:filename_list[psm.FileName].rfind('.')]
            row_list.append(str(psm.ScanNumber)+fileID+str(psm.ParentCharge)) # SpecId
            if psm.RealLabel == LabelFwd: # Label
                row_list.append("1")
            elif psm.RealLabel == LabelReserve:
                row_list.append("1")
            else:
                row_list.append("-1")
            row_list.append(str(psm.ScanNumber)+fileID) # ScanNr
            row_list.append(str(psm.MeasuredParentMass)) # ExpMass
            row_list.append(str(psm.CalculatedParentMass)) # CalcMass
            row_list.append(str(psm.OPSC)) # PEP
            row_list.append(str(psm.SPSC)) # PRO
            row_list.append(str(psm.lfScores[0])) # MVH
            row_list.append(str(psm.score_differential_list[9])) # deltMVH
            row_list.append(str(psm.lfScores[1])) # Xcorr
            row_list.append(str(psm.score_differential_list[10])) # deltXcorr
            row_list.append(str(psm.lfScores[2])) # WDP
            row_list.append(str(psm.score_differential_list[11])) # deltWDP
            row_list.append(str(len(psm.OriginalPeptide) - 2)) # PepLen
            row_list.append(str(psm.dM)) # dM
            if abs(psm.dM) > PSM.fNeutronMass:
                print("error mass difference")
            row_list.append(str(psm.fMassDiff)) # absdM
            for i in range(1, 4): # Charge 1-4
                if psm.ParentCharge > i:
                    row_list.append("0")
                else:
                    row_list.append("1")
            if psm.ParentCharge > 4: # Charge5
                row_list.append("1")
            else:
                row_list.append("0")
            row_list.append(str(psm.NMC)) # enzInt
            row_list.append(psm.IdentifiedPeptide.replace("[", "-.").replace("]", ".-")) # Peptide
            row_list.append("\t".join(psm.protein_list)) # Protein
            fw.write("\t".join(row_list))
            fw.write("\n")

def debug_mode(psm_list, output_folder):
    generatePINPercolator(psm_list, output_folder[:-1])

## script entrance
def main(argv=None):
    if argv is None:
        argv = sys.argv
    
    # Display work start and time record
    start_time = datetime.now()
    sys.stderr.write('[{}] Beginning {} ({})\n'.format(curr_time(), os.path.basename(__file__), get_version()))
    sys.stderr.write('------------------------------------------------------------------------------\n')
    
    # parse options
    (input_file, config_file, output_folder, debug_code) = parse_options(argv)
    
    # get the configuration parameters
    config_dict = parse_config(config_file)
    
    if config_dict[pep_iden_str + search_type_str] == 'SIP':
        # sip_filtering(input_file, config_dict, output_folder, start_time)
        sip_filtering_LR(input_file, config_dict, output_folder, start_time)
        return
    
    # read the big psm table
    sys.stderr.write('[Step 1] Parse options and read PSM file:                   Running -> ')
    # find out the train and testing ratio
    find_train_test_ratio(config_dict)
    (psm_list, base_out) = read_psm_table(input_file)
    sys.stderr.write('Done!\n')
    
    # get score agreement info
    sys.stderr.write('[Step 2] Feature extraction:                                Running -> ')
    remark_concensus(psm_list)
    # mass filtering
    psm_list = mass_filter(psm_list, config_dict)
    # generate features
    generate_Prophet_features_test(psm_list, config_dict)
    # set feature all PSMs
    for oPsm in psm_list:
        oPsm.get_feature_final_list()
    sys.stderr.write('Done!\n')
    
    # machine learning
    if debug_code != "":
        debug_mode(psm_list, output_folder)
        return
    
    sys.stderr.write('[Step 3] Train ML and re-rank PSMs:                         Running -> ')
    del feature_selection_list[:]
    # feature_selection_list.extend([1, 2, 3, 5, 15, 16, 17, 24, 26, 28])
    feature_selection_list.extend([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    psm_filtered_list = logistic_regression_no_category(psm_list, config_dict)
    sys.stderr.write('Done!\n')

    # write output
    sys.stderr.write('[Step 4] Report output:                                     Running -> ')
    (psm_txt_file_str, pep_txt_file_str) = generate_psm_pep_txt(base_out, output_folder, psm_filtered_list, config_dict)
    sys.stderr.write('Done!\n')
    
    # Time record, calculate elapsed time, and display work end
    finish_time = datetime.now()
    duration = finish_time - start_time
    sys.stderr.write('------------------------------------------------------------------------------\n')
    sys.stderr.write('[{}] Ending {} \n'.format(curr_time(), os.path.basename(__file__)))
    sys.stderr.write('Run complete [%s elapsed]\n' %  format_time(duration))

'''
split psm table to 2 parts
1. below 3%
2. above 3%
'''
def sip_filtering_LR_try_2(input_file, config_dict, output_folder, start_time):
    # read the big psm table
    sys.stderr.write('[Step 1] Parse options and read PSM file:                   Running -> ')
    # reading data from PSM table
    (psm_list, base_out) = read_psm_table(input_file)
    # mass filtering
    psm_list = mass_filter(psm_list, config_dict)
    sys.stderr.write('Done!\n')
    
    # find score-cutoff
    sys.stderr.write('[Step 2] Re-rank PSMs and find score cutoff:                Running -> ')
    psm_filtered_list = cutoff_filtering(psm_list, config_dict)
    sys.stderr.write('Done!\n')
    
    # train a machine learning model: logistic regression
    sys.stderr.write('[Step 3] Feature extraction:                                Running -> ')
    # generate features
    get_percentage_back(psm_list)
    
    psm_list_below_3pct = []
    psm_list_above_3pct = []
    
    for psm_obj in psm_list:
        if int(psm_obj.pct_s) <= 3:
            psm_list_below_3pct.append(psm_obj)
        else:
            psm_list_above_3pct.append(psm_obj)
    
    generate_Prophet_features_test(psm_list_below_3pct, config_dict)
    generate_Prophet_features_test(psm_list_above_3pct, config_dict)
    # generate_Prophet_features_group(psm_list, config_dict)
    # set feature all PSMs
    for oPsm in psm_list:
        oPsm.get_feature_final_list()
    # reset the config
    config_reset(config_dict)
    # find out the train and testing ratio
    find_train_test_ratio(config_dict)
    # mark the positive train data
    for oPsm in psm_list:
        oPsm.set_real_label()
    for oPsm in psm_filtered_list:
        if oPsm.RealLabel == LabelFwd:
            oPsm.TrainingLabel = LabelSipTrainFwd
    sys.stderr.write('Done!\n')
    
    # machine learning
    sys.stderr.write('[Step 4] Train ML and re-rank PSMs:                         Running -> ')
    del feature_selection_list[:]
    # feature_selection_list.extend([1, 2, 3, 5, 15, 16, 17, 24, 26, 28])
    # feature_selection_list.extend([1, 2, 15, 24, 26, 28])
    feature_selection_list.extend([0, 3, 4, 7, 8, 9])
    psm_filtered_below_3pct_list = logistic_regression_sip(psm_list_below_3pct, config_dict)
    psm_filtered_above_3pct_list = logistic_regression_sip(psm_list_above_3pct, config_dict)
    psm_filtered_list = []
    psm_filtered_list.append(psm_filtered_below_3pct_list)
    psm_filtered_list.append(psm_filtered_above_3pct_list)
    sys.stderr.write('Done!\n')
    
    # write output
    sys.stderr.write('[Step 5] Report output:                                     Running -> ')
    (psm_txt_file_str, pep_txt_file_str) = generate_psm_pep_txt(base_out, output_folder, psm_filtered_list, config_dict)
    sys.stderr.write('Done!\n')
    
    # Time record, calculate elapsed time, and display work end
    finish_time = datetime.now()
    duration = finish_time - start_time
    sys.stderr.write('------------------------------------------------------------------------------\n')
    sys.stderr.write('[{}] Ending {}\n'.format(curr_time(), os.path.basename(__file__)))
    sys.stderr.write('Run complete [%s elapsed]\n' %  format_time(duration))

'''
if sip mode
use WDP to find the positive training data
train the logistic model with 3 scores
re-rank PSM based on the predicted probabilities
'''
def sip_filtering_LR(input_file, config_dict, output_folder, start_time):
    # read the big psm table
    sys.stderr.write('[Step 1] Parse options and read PSM file:                   Running -> ')
    # reading data from PSM table
    (psm_list, base_out) = read_psm_table(input_file)
    # mass filtering
    psm_list = mass_filter(psm_list, config_dict)
    sys.stderr.write('Done!\n')
    
    # find score-cutoff
    sys.stderr.write('[Step 2] Re-rank PSMs and find score cutoff:                Running -> ')
    psm_filtered_list = cutoff_filtering(psm_list, config_dict)
    sys.stderr.write('Done!\n')
    
    # write output
    sys.stderr.write('[Step 3] Report output:                                     Running -> ')
    (psm_txt_file_str, pep_txt_file_str) = generate_psm_pep_txt(base_out, output_folder, psm_filtered_list, config_dict)
    sys.stderr.write('Done!\n')
    
    # Time record, calculate elapsed time, and display work end
    finish_time = datetime.now()
    duration = finish_time - start_time
    sys.stderr.write('------------------------------------------------------------------------------\n')
    sys.stderr.write('[{}] Ending {}\n'.format(curr_time(), os.path.basename(__file__)))
    sys.stderr.write('Run complete [%s elapsed]\n' %  format_time(duration))

    
if __name__ == '__main__':
    sys.exit(main())
