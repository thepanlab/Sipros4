#!/usr/bin/python

"""
sipros_post_module.py

sipros_post_module.py module includes common classes, definitions,
and funcitons for sipros post-processing programs.

Created by Tae-Hyuk (Ted) Ahn on 10/10/2012.
Copyright (c) 2012 Tae-Hyuk Ahn (ORNL). Allrights reserved.

Modified by Xuan Guo on 07/13/2016
"""

# # Import standard modules
import sys, os, re, math
from datetime import datetime
from collections import namedtuple
from multiprocessing import Process
try:
    from sets import Set
except ImportError:
    pass

SIP_WDP_score_idx = 0

# # Class Usage
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


# # Class for ignoring comments '#' in sipros file
class CommentedFile:
    def __init__(self, f, comment_string="#"):
        self.f = f
        self.comment_string = comment_string
    def next(self):
        line = self.f.next()
        while line.startswith(self.comment_string):
            line = self.f.next()
        return line
    def __iter__(self):
        return self


# # Exit system with error message
def die(msg=None):
    if msg is not None:
        print >> sys.stderr, msg
        sys.exit(1)


# # Returns the current time in a nice format
def curr_time():
    curr_time = datetime.now()
    return curr_time.strftime("%c")


# # Format time as a pretty string
def format_time(td):
    hours = td.seconds // 3600
    minutes = (td.seconds % 3600) // 60
    seconds = td.seconds % 60
    return '%02d:%02d:%02d' % (hours, minutes, seconds)


# # Find string between two substrings
def find_between(s, first, last):
    try:
        start = s.index(first) + len(first)
        end = s.index(last, start)
        return s[start:end]
    except ValueError:
        return ""


# # Division error handling
def divide(x, y):
    try:
        result = x / y
    except ZeroDivisionError as detail:
        print >> sys.stderr, 'Handling run-time error:', detail
        die('Program exit!')
    else:
        return result


# # Check file exist
def check_file_exist(filename):

    try:
        with open(filename) as _f: pass
    except IOError as _e:
        print >> sys.stderr, '\nCannot open', filename
        die("Program exit!")


# # Get file(s) list in working dir with specific file extension
def get_file_list_with_ext(working_dir, file_ext):

    # define sipros file extension
    file_list = []

    # working directory
    if os.path.exists(working_dir):
        for file_name in os.listdir(working_dir):

            # check the file extension
            if file_name.endswith(file_ext):
                file_path_name = os.path.join(working_dir, file_name)
                # file_path_name = working_dir + file_name
                file_list.append(file_path_name)

        if len(file_list) == 0:
            print >> sys.stderr, "\nCannot open %s file(s)." % (file_ext)
            die("Program exit!")
        file_list = sorted(file_list)

    else:
        print >> sys.stderr, "\nCannot open working directory", working_dir
        die("Program exit!")

    return file_list

# # Get base output filename with input file list and base_out_default
def get_base_out(file_list, base_out_default, working_dir):

    # Get base output with common prefix
    base_out = os.path.commonprefix(file_list)
    base_out_filename = base_out.split('/')[-1]

    # If base common prefix ends with '.pep.txt', then remove '.pep.txt'
    base_out = base_out.replace(".pep.txt", "_")
    
    # If base common prefix ends with '.tab', then remove '.tab'
    base_out = base_out.replace(".tab", "_")

    # If base_out file name is less than 5, then use default baseout
    if len(base_out_filename) < 5:
        base_out = os.path.join(working_dir, base_out_default)

    # If base common prefix ends with '_' or '.', then remove
    base_out = base_out[:-1] if (base_out[-1] in ('_', '.')) else base_out

    return base_out


# # list_to_string
# # if single element, then just convert to string
# # if multiple elements, then bracket {A,B}
def list_to_string(input_list):

    if len(input_list) > 1:
        converted_str = '{' + ','.join(input_list) + '}'
    else:
        converted_str = ''.join(input_list)

    return converted_str


# # list_to_bracket
# # bracket the list
def list_to_bracket(input_list):

    converted_str = '{' + ','.join(input_list) + '}'

    return converted_str


# # Class for sipros fields object
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

# # Class for spectrum fields object
class SpectrumFields(namedtuple('SpectrumFields',
        ['x',
        'Filename',
        'ScanNumber',
        'ParentCharge',
        'MeasuredParentMass',
        'ScanType',
        'SearchName',
        'TotalIntensity',
        'MaxIntensity'])):
    def __init__(self):
        self.data = self

# # Class for psm fields object
class PsmFields(namedtuple('PsmFields',
        ['x',
        'OriginalPeptide',
        'IdentifiedPeptide',
        'CalculatedParentMass',
        'MVH',
        'Xcorr',
        'WDP',
        'ProteinNames'])):
    def __init__(self):
        self.data = self


# # Class for sipros4 fields object
class Sipros4Fields(namedtuple('SiprosFields',
        ['Filename',
        'ScanNumber',
        'ParentCharge',
        'MeasuredParentMass',
        'CalculatedParentMass',
        'ScanType',
        'SearchName',
        'Rank',
        'MVH',
        'Xcorr',
        'WeightDotSum',
        'IdentifiedPeptide',
        'OriginalPeptide',
        'ProteinNames',
        'AveAtom',
        'StdAtom'])):

    def __init__(self):
        self.data = self


# # Class for PsmOutFields object
class PsmOutFields(namedtuple('PsmOutFields',
        ['Filename',
         'ScanNumber',
         'ParentCharge',
         'MeasuredParentMass',
         'CalculatedParentMass',
         'MassErrorDa',  # CalculatedParentMass - MeasuredParentMass
         'MassErrorPPM',  # MassErrorDa / CalculatedParentMass
         'ScanType',
         'SearchName',
         'ScoringFunction',
         'Score',
         'DeltaZ',  # The difference score between the rank 1 and 2
         'DeltaP',  # The difference score between isoform
         'IdentifiedPeptide',
         'OriginalPeptide',
         'ProteinNames',
         'ProteinCount',
         'TargetMatch'])):
    def __init__(self):
        self.data = self


# # Class for PsmOutFields object (for sipro4)
class Psm4OutFields(namedtuple('PsmOutFields',
        ['Filename',
         'ScanNumber',
         'ParentCharge',
         'MeasuredParentMass',
         'CalculatedParentMass',
         'MassErrorDa',  # CalculatedParentMass - MeasuredParentMass
         'MassErrorPPM',  # MassErrorDa / CalculatedParentMass
         'ScanType',
         'SearchName',
         'ScoringFunction',
         'Score',
         'DeltaZ',  # The difference score between the rank 1 and 2
         'DeltaP',  # The difference score between isoform
         'IdentifiedPeptide',
         'OriginalPeptide',
         'ProteinNames',
         'ProteinCount',
         'TargetMatch',
         'AveAtom',
         'StdAtom'])):
    def __init__(self):
        self.data = self


# # Class for PepSubFields object
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


# # Class for PepSubFields object
class PepDataFields(namedtuple('PepDataFields',
        ['IdentifiedPeptide',
         'ParentCharge',
         'BestScore',
         'ProteinNames'])):

    def __init__(self):
        self.data = self


# # Class for PepOutFields object
class PepOutFields(namedtuple('PepOutFields',
        ['IdentifiedPeptide',  # 0
         'ParentCharge',  # 1
         'OriginalPeptide',  # 2
         'ProteinNames',  # 3
         'ProteinCount',  # 4
         'TargetMatch',  # 5
         'SpectralCount',  # 6 number of PSMs matched to this peptide
         'BestScore',  # 7 the highest score of those PSMs
         'PSMs',  # 8 a list of PSMs matched to this peptide. Use {Filename[ScanNumber],Filename[ScanNumber]} format
         'ScanType',  # 9 ScanType
         'SearchName'])):  # 10 SearchName

    def __init__(self):
        self.data = self


# # Class for pretty float
class PrettyFloat(float):
    def __repr__(self):
        return "%0.5f" % self


# # A range function, that does accept float increments
def frange(start, end=None, inc=None):

    if end == None:
        end = start + 0.0
        start = 0.0

    if inc == None:
        inc = 1.0

    L = []
    while 1:
        fnext = start + len(L) * inc
        if inc > 0 and fnext >= end:
            break
        elif inc < 0 and fnext <= end:
            break
        L.append(fnext)

    return L


# # check sub list
def check_sub_list(list_A, list_B):

    check_status = True

    for list_A_item in list_A:
        if list_A_item not in list_B:
            check_status = False
        else:
            continue

    return check_status


# # get item list from parenthesis string as {AA,BB}
def get_item_list(input_string):

    input_string = input_string[1:-1]
    item_list = re.split(r"\s*[,]\s*", input_string.strip())
    item_list_new = []
    for item_one in item_list:
        if item_one not in item_list_new:
            item_list_new.append(item_one)
    return item_list_new


# # get_protein_count
def get_protein_count(protein_names):

    input_string = protein_names[1:-1]
    item_list = re.split(r"\s*[,]\s*", input_string.strip())
    protein_count = len(item_list)

    return protein_count


# # set float digit
def set_float_digit(input_val):

    if input_val is float:
        output_val = str("{0:.5f}".format(round(input_val, 5)))
    else:
        output_val = str(input_val)

    return output_val

# # peptide delete residues
def peptide_delete_residues(peptide_string):

    try:
        left_braket_index = peptide_string.index('[')
        right_braket_index = peptide_string.index(']')
        if len(peptide_string) > right_braket_index + 1:
            if peptide_string[right_braket_index + 1].isalpha():
                peptide_output = peptide_string[left_braket_index:right_braket_index + 1]
            else:
                peptide_output = peptide_string[left_braket_index:right_braket_index + 2]
        else:
            peptide_output = peptide_string[left_braket_index:right_braket_index + 1]

        return peptide_output
    except AttributeError:
        print >> sys.stderr, '\nCannot parse peptide correctly.\n'
        die("Program exit!")


# # merge protein names
def merge_protein_names(first_protein_names, second_protein_names):

    first_protein_list = get_item_list(first_protein_names)
    second_protein_list = get_item_list(second_protein_names)

    merge_protein_list = list(set(first_protein_list + second_protein_list))

    merge_protein_names = list_to_bracket(merge_protein_list)

    return merge_protein_names


# # Task wrapper
class PsmPack:

    def __init__(self, _iSize=1000, _iStartScanNumber=0):
        self.iSize = _iSize
        self.lPsm = []
        for _i in range(_iSize):
            self.lPsm.append([])
        self.iStartScanNumber = _iStartScanNumber
        self.bEmpty = True
        self.current = 0

    def add(self, lOnePsm, iScanNumber):
        self.lPsm[iScanNumber - self.iStartScanNumber].append(lOnePsm)
        self.bEmpty = False

    def empty(self):
        return self.bEmpty

    def __iter__(self):
        self.current = 0
        return self

    def next(self):
        if self.current >= self.iSize:
            raise StopIteration
        else:
            while not self.lPsm[self.current]:
                self.current += 1
                if self.current >= self.iSize:
                    raise StopIteration
            self.current += 1
            return self.lPsm[self.current - 1]

fNeutronMass = 1.00867108694132 # it is Neutron mass
# # the mass difference, inverted, larger better.
def MassDiff(oPepScores):
    fDiff = oPepScores.fCalculatedParentMass - oPepScores.fMeasuredParentMass
    fTemp = fDiff
    fCeil = 0
    fDown = 0
    if fDiff >= 0:
        fDiff = fTemp
        fCeil = math.ceil(fTemp)*fNeutronMass
        fFloor = math.floor(fTemp)*fNeutronMass
        if fFloor > fTemp:
            fFloor -= fNeutronMass
        if fCeil - fNeutronMass > fTemp:
            fCeil -= fNeutronMass
        if fTemp > fCeil - fTemp:
            fTemp = fCeil - fTemp
        if fDiff > fDiff - fFloor:
            fDiff = abs(fDiff - fFloor)
        if abs(fTemp) < abs(fDiff):
            fDiff = fTemp
    else:
        fCeil = math.ceil(fDiff)*fNeutronMass
        if fCeil < fDiff:
            fCeil += fNeutronMass
        fFloor = math.floor(fDiff)*fNeutronMass
        if fFloor + fNeutronMass < fDiff:
            fFloor += fNeutronMass
        fDiff = fTemp
        if abs(fTemp) > fCeil - fTemp:
            fTemp = fCeil - fTemp
        if abs(fDiff) > fDiff - fFloor:
            fDiff = fDiff - fFloor
        fTemp = abs(fTemp)
        fDiff = abs(fDiff)
        if fTemp < fDiff:
            fDiff = fTemp
    return -fDiff

# # the PTM score, the count of PTM, larger better
def PtmScore(oPepScores):
    s1 = ''.join([char if char.isalnum() else '$' for char in oPepScores.sIdentifiedPeptide ])
    return -(s1.count('$') - 2)

# # Rank Product
def RankProductInvert(liRank):
    fProduct = 1.0
    for i in liRank:
        fProduct *= i
    return 1 / fProduct

# # Pep info in the Spe2Pep file
# # non SIP mode
class PepScores:

    def __init__(self, _fMeasuredParentMass, _iCharge, _sSearchName, sPeptideLine, isSIP=False):
        self.fMeasuredParentMass = _fMeasuredParentMass
        self.iCharge = _iCharge
        asWords = sPeptideLine.split('\t')
        self.sIdentifiedPeptide = peptide_delete_residues(asWords[1])
        self.sOriginalPeptide = peptide_delete_residues(asWords[2])
        # self.sIdentifiedPeptide = peptide_delete_residues(asWords[2])
        # self.sOriginalPeptide = peptide_delete_residues(asWords[1])
        self.fCalculatedParentMass = float(asWords[3])
        self.lfScores = []
        self.liRanks = []
        self.lfScoreDiff = []  # difference between the current one with the second best one
        for e in asWords[4:-1]:
            self.lfScores.append(float(e))
        # remove the {}
        self.sProteinNames = (asWords[-1])[1:-1]
        self.fRankProduct = 0.0
        self.iRank = 0
        self.sSearchName = _sSearchName
        self.fPct = 0
        if isSIP:
            self.fPct = float(_sSearchName[_sSearchName.find('_') + 1:_sSearchName.find("Pct")])
        # delta -> comet way
        # diff -> difference between current one to the next best one
        # diffnor -> difference between current one to the next best one normalized by the current one
        self.lfDeltaRankProduct = []
        self.lfDeltaRankScore = []
        self.lfDiffRankProduct = []
        self.lfDiffRankScore = []
        self.lfDiffNorRankProduct = []
        self.lfDiffNorRankScore = []
        self.DeltaP = 'NA'
        if len(self.sOriginalPeptide) != len(self.sIdentifiedPeptide):
            self.DeltaP = 1

def numberTopRanks(liRanks):
    iCount = 0
    for i in liRanks:
        if i == 1:
            iCount += 1
    return iCount

# # count pep ranked as top by how many scores
def get_number_Top_Ranks(pep, all_top_peps):
    iCount = 0
    for i in range(len(pep.lfScores)):
        if pep.lfScores[i] >= all_top_peps[i].lfScores[i]:
            pep.liRanks[i] = 1
            iCount += 1
    return iCount

def zero_divide(a, b):
    if b != 0:
        return a / b
    else:
        return 0

# # PSM info in the Spe2Pep file
class PepSpectrumMatch:

    iPurgeTrigger = 100
    iRankSave = 20

    def __init__(self, sSpectrumLine, SIP=False):
        asWords = sSpectrumLine.split('\t')
        self.sFileName = asWords[1]
        self.iScanNumber = int(asWords[2])
        self.iParentCharge = int(asWords[3])
        self.fMeasuredParentMass = float(asWords[4])
        self.sScanType = asWords[5]
        self.sSearchName = asWords[6]
        # self.fTotalIntensity = float(asWords[7])
        # self.fMaxIntensity = float(asWords[8])
        self.sRTime = asWords[7]
        self.lPepScores = []
        self.oBestPep = None
        self.oSecondBestPep = None
        self.oRestPep = None
        self.lTopPep = []
        self.pep_rank_list = []
        self.SIPmode = SIP

    def addPepScores(self, pep):
        for e in self.lPepScores:
            if e.sIdentifiedPeptide == pep.sIdentifiedPeptide:
                if e.lfScores[0] < pep.lfScores[0]:
                    e.lfScores = pep.lfScores
                    e.fCalculatedParentMass = pep.fCalculatedParentMass
                    e.sSearchName = pep.sSearchName
                    e.fPct = pep.fPct
                elif e.lfScores[0] == pep.lfScores[0] and self.SIPmode:
                    if e.fPct > pep.fPct:
                        e.fCalculatedParentMass = pep.fCalculatedParentMass
                        e.sSearchName = pep.sSearchName
                        e.fPct = pep.fPct
                words = pep.sProteinNames.split(',')
                for sProtein in words:
                    if e.sProteinNames.find(sProtein) == -1:
                        e.sProteinNames += ','
                        e.sProteinNames += sProtein
                return
        self.lPepScores.append(pep)
        if len(self.lPepScores) > self.iPurgeTrigger:
            if self.SIPmode:
                self.purgeSIP()
            else:
                self.purge()
            
    def purgeSIP(self):
        # sort pep accoring to SIP_WDP_score
        lPep = sorted(self.lPepScores, key=lambda pep: (-pep.lfScores[SIP_WDP_score_idx], pep.fPct))
        pep_new_list = []
        # keep all the pep with the highest score
        pep_new_list.extend(lPep[0:self.iPurgeTrigger/2])

    def purge(self):
        iNumScores = len(self.lPepScores[0].lfScores)
        for j in self.lPepScores:
            del j.liRanks[:]
        for i in range(iNumScores):
            lPep = sorted(self.lPepScores, key=lambda pep: (-pep.lfScores[i], -MassDiff(pep), -PtmScore(pep), pep.sIdentifiedPeptide, pep.fPct))
            iRank = 1
            for j in lPep:
                j.liRanks.append(iRank)
                iRank += 1
        liRanksNew = []
        for j in self.lPepScores:
            if any(i <= self.iRankSave for i in j.liRanks):
                liRanksNew.append(j)
        self.lPepScores = liRanksNew
     
    def ranking_sip(self):
        # only kep the top score psm, or the least mass difference psm, or lowest ptm score psm, or the psm with the smallest percent
        lPep = sorted(self.lPepScores, \
                      key=lambda pep: (-pep.lfScores[SIP_WDP_score_idx], \
                                       -MassDiff(pep), -PtmScore(pep), \
                                       pep.fPct, pep.sIdentifiedPeptide))
        self.pep_rank_list.append(lPep)
        iRank = 1
        for j in lPep:
            j.liRanks.append(iRank)
            iRank += 1
        # score rank -> score differential
        for j in range(0, len(lPep) - 1):
            lPep[j].lfDiffRankScore.append(lPep[j].lfScores[SIP_WDP_score_idx] - \
                                           lPep[j + 1].lfScores[SIP_WDP_score_idx])
            lPep[j].lfDiffRankScore.append(0) # for the format sake, SIP does not use the other two score
            lPep[j].lfDiffRankScore.append(0) # for the format sake, SIP does not use the other two score

        lPep[len(lPep) - 1].lfDiffRankScore.append(0)
        lPep[len(lPep) - 1].lfDiffRankScore.append(0)
        lPep[len(lPep) - 1].lfDiffRankScore.append(0)
        del self.lTopPep[:]
        self.lTopPep.append(lPep[0])
        
    def ranking(self):
        del self.lTopPep[:]
        for j in self.lPepScores:
            del j.liRanks[:]
            del j.lfScoreDiff[:]
        iNumScores = len(self.lPepScores[0].lfScores)
        
        for i in range(iNumScores):
            lPep = sorted(self.lPepScores, key=lambda pep: (-pep.lfScores[i], -MassDiff(pep), -PtmScore(pep), pep.sIdentifiedPeptide, pep.fPct))
            self.pep_rank_list.append(lPep)
            iRank = 1
            for j in lPep:
                j.liRanks.append(iRank)
                iRank += 1
            # score rank -> score differential
            for j in range(0, len(lPep) - 1):
                lPep[j].lfDiffRankScore.append(lPep[j].lfScores[i] - lPep[j + 1].lfScores[i])

            lPep[len(lPep) - 1].lfDiffRankScore.append(0)
            
            self.lTopPep.append(lPep[0])
        
        # if SA == 1, calculate DeltaP individually
        for s in range(iNumScores):
            for lPep_local in self.pep_rank_list:
                # contain PTM
                if len(lPep_local[0].sIdentifiedPeptide) != len(lPep_local[0].sOriginalPeptide):
                    #pep_sorted_str = ''.join(sorted(lPep_local[0].sIdentifiedPeptide))
                    for pep in lPep_local:
                        #pep_sorted_compared_str = ''.join(sorted(pep.sIdentifiedPeptide))
                        #if pep_sorted_str == pep_sorted_compared_str and \
                        #pep.sIdentifiedPeptide != lPep_local[0].sIdentifiedPeptide :
                        
                        if pep.sIdentifiedPeptide != lPep_local[0].sIdentifiedPeptide and \
                        abs(pep.fCalculatedParentMass - lPep_local[0].fCalculatedParentMass) < 0.00005 and \
                        abs(math.fabs(pep.fMeasuredParentMass - pep.fCalculatedParentMass) - math.fabs(lPep_local[0].fMeasuredParentMass - lPep[0].fCalculatedParentMass)) < 0.00005 and \
                        pep.sOriginalPeptide == lPep_local[0].sOriginalPeptide and \
                        len(pep.sIdentifiedPeptide) == len(lPep_local[0].sIdentifiedPeptide):
                        
                            avg_deltaP = 0.0
                            for si in range(iNumScores):
                                # if self.lTopPep[si].lfScores[si] != 0 and lPep_local[0].liRanks[si] == 1:
                                if self.lTopPep[si].lfScores[si] != 0:
                                    avg_deltaP += (lPep_local[0].lfScores[si] - pep.lfScores[si])/self.lTopPep[si].lfScores[si]
                            avg_deltaP /= float(iNumScores)
                            lPep_local[0].DeltaP = avg_deltaP
                            break
    
    def all_top_ranked_psm(self):
        str_list = []
        for pep in self.lTopPep:
            feature_list = []
            feature_list.append(self.sFileName)
            feature_list.append(str(self.iScanNumber))
            feature_list.append(str(pep.iCharge))
            feature_list.append(str(self.fMeasuredParentMass))
            feature_list.append(self.sScanType)
            feature_list.append(pep.sSearchName)
            feature_list.append(pep.sIdentifiedPeptide)
            feature_list.append(pep.sOriginalPeptide)
            feature_list.append(str(pep.fCalculatedParentMass))
            feature_list.extend((str(x) for x in pep.lfScores))
            feature_list.append('{'+pep.sProteinNames+'}')
            # feature_list.append(str(num_agreement(pep.liRanks)))
            # feature_list.extend((str(x) for x in pep.lfDeltaRankProduct))
            # feature_list.extend((str(x) for x in pep.lfDeltaRankScore))
            # feature_list.extend((str(x) for x in pep.lfDiffRankProduct))
            feature_list.extend((str(x) for x in pep.lfDiffRankScore))
            # feature_list.extend((str(x) for x in pep.lfDiffNorRankProduct))
            # feature_list.extend((str(x) for x in pep.lfDiffNorRankScore))
            feature_list.append(self.sRTime)
            # feature_list.append((str(pep.iRank)))
            feature_list.append((str(pep.DeltaP)))
            str_list.append('\t'.join(feature_list))
        return '\n'.join(str_list)
    
    def all_top_5_ranked_psm(self, top_n = 5):
        str_list = []
        pep_set = Set()
        for l in self.pep_rank_list:
            n = len(l)
            if n > top_n:
                n = top_n
            for i in range(n):
                pep_set.add(l[i])
        
        for pep in pep_set:
            feature_list = []
            feature_list.append(self.sFileName)
            feature_list.append(str(self.iScanNumber))
            feature_list.append(str(self.iParentCharge))
            feature_list.append(str(self.fMeasuredParentMass))
            feature_list.append(self.sScanType)
            feature_list.append(pep.sSearchName)
            feature_list.append(pep.sIdentifiedPeptide)
            feature_list.append(pep.sOriginalPeptide)
            feature_list.append(str(pep.fCalculatedParentMass))
            feature_list.extend((str(x) for x in pep.lfScores))
            feature_list.append('{'+pep.sProteinNames+'}')
            # feature_list.append(str(num_agreement(pep.liRanks)))
            # feature_list.extend((str(x) for x in pep.lfDeltaRankProduct))
            # feature_list.extend((str(x) for x in pep.lfDeltaRankScore))
            # feature_list.extend((str(x) for x in pep.lfDiffRankProduct))
            feature_list.extend((str(x) for x in pep.lfDiffRankScore))
            # feature_list.extend((str(x) for x in pep.lfDiffNorRankProduct))
            # feature_list.extend((str(x) for x in pep.lfDiffNorRankScore))
            feature_list.append(self.sRTime)
            # feature_list.append((str(pep.iRank)))
            feature_list.append((str(pep.DeltaP)))
            str_list.append('\t'.join(feature_list))
        return '\n'.join(str_list)

    def removeReverse(self, lProteins):
        newlist = []
        for sProtein in lProteins:
            if not sProtein.startswith('Rev_'):
                newlist.append(sProtein)
        return newlist

# # lOnePsm: 
# # + spectrum 
# # * peptide
# # * peptide 
def SelectTopRankedPsm(lOnePsm, isSIP=False):
    psm_dict = {}
    psm = PepSpectrumMatch(lOnePsm[0][0], isSIP)
    psm_dict[psm.iParentCharge] = psm
    '''
    if psm.iScanNumber == 4299:
        print 'check'
    '''
    iCharge = 0
    sSearchName = ''
    for PsmInOneFile in lOnePsm:
        for sline in PsmInOneFile:
            if sline[0] == '+':
                # a spectrum line
                iCharge = get_charge(sline)
                sSearchName = get_search_name(sline)
                if iCharge not in psm_dict:
                    psm = PepSpectrumMatch(sline, isSIP)
                    psm_dict[iCharge] = psm
            else:
                # a peptide line
                pep = PepScores(psm_dict[iCharge].fMeasuredParentMass, iCharge, sSearchName, sline, isSIP)
                psm_dict[iCharge].addPepScores(pep)
    # sorting and then ranking
    psm_list = []
    for k,v in psm_dict.items():
        if isSIP:
            v.ranking_sip()
            if v.lTopPep[0].lfScores[SIP_WDP_score_idx] <= 0:
                continue
        else:
            v.ranking()
        psm_list.append(v)
    return psm_list

# # peak a line from the file
def peek_line(f):
    pos = f.tell()
    sline = f.readline()
    f.seek(pos)
    return sline

# # get the scan number from a line
def get_scan_number(sLine, sDelimiter='\t', iFirstDelimiter=2):
    iPos = -1
    while iFirstDelimiter > 0:
        iPos = sLine.find(sDelimiter, iPos + 1)
        iFirstDelimiter -= 1
    iBegin = iPos + 1
    iPos = sLine.find(sDelimiter, iBegin)
    iScanNumber = int(sLine[iBegin:iPos])
    return iScanNumber

# # get the charge from a line
def get_charge(sLine, sDelimiter='\t', iFirstDelimiter=3):
    iPos = -1
    while iFirstDelimiter > 0:
        iPos = sLine.find(sDelimiter, iPos + 1)
        iFirstDelimiter -= 1
    iBegin = iPos + 1
    iPos = sLine.find(sDelimiter, iBegin)
    iCharge = int(sLine[iBegin:iPos])
    return iCharge

# # get the search name from a line
def get_search_name(sLine, sDelimiter='\t', iFirstDelimiter=6):
    iPos = -1
    while iFirstDelimiter > 0:
        iPos = sLine.find(sDelimiter, iPos + 1)
        iFirstDelimiter -= 1
    iBegin = iPos + 1
    iPos = sLine.find(sDelimiter, iBegin)
    sSearchName = sLine[iBegin:iPos]
    return sSearchName

# # get the PSM with scan number less than the upper scan number bound
def get_psm(f, lPsm, sSpectrum='+', sPeptide='*', iUpperScanNumber=0):
    bEof = True
    lOnePsm = []
    while True:
        pos = f.tell()
        sline = f.readline().strip()
        # # end of file
        if not sline:
            break
        iScanNumber = 0
        if sline[0] == sSpectrum:
            iScanNumber = get_scan_number(sline)
            if iScanNumber < iUpperScanNumber:
                bEof = False
                lOnePsm = []
                lOnePsm.append(sline)
                lPsm.add(lOnePsm, iScanNumber)
            else:
                # roll back the previous position
                f.seek(pos)
                bEof = False
                break
        else:
            lOnePsm.append(sline)
    return bEof

# # skip the comment area and the header
def skip_comment(f, sComment='#', iLineHeader=0):
    pos = f.tell()
    sline = f.readline()
    while(sline[0] == sComment):
        pos = f.tell()
        sline = f.readline()
    f.seek(pos)
    for _i in range(iLineHeader):
        f.readline()
    # print f.readline()

# # Spe2Pep reader
class Spe2PepReader(Process):

    def __init__(self, queue=None, name=None, searchname=None, inputFolder=None):
        super(Spe2PepReader, self).__init__()
        self.name = name
        self.qPsmUnprocessed = queue;
        self.iNumScanProcessed = 0
        self.sSearchName = searchname
        self.FileList = []
        self.iScanInterval = 1000
        self.sInputFolder = inputFolder

    # # list files with 'Spe2Pep.txt' extensions
    # # put the search results for the same FT2 file into a list
    def categorizeSpe2PepFile(self, sWorkingDirectory):
        lFileList = get_file_list_with_ext(sWorkingDirectory, 'Spe2Pep.txt')
        sFt2Name = ''
        lFt2Name = []
        iIndexFt2 = 0
        for sFileName in lFileList:
            iPos = sFileName.rfind(self.sSearchName)
            if iPos != -1:
                # a '.' is before the search name, so iPos-1
                sFt2Name = sFileName[0:iPos - 1]
                if sFt2Name in lFt2Name:
                    iIndexFt2 = lFt2Name.index(sFt2Name)
                    self.FileList[iIndexFt2].append(sFileName)
                else:
                    lFt2Name.append(sFt2Name)
                    lNewFt2 = []
                    lNewFt2.append(sFileName)
                    self.FileList.append(lNewFt2)

    def readData(self):
        for _id, lFiles in enumerate(self.FileList):
            lFileReader = []
            for sFiles in lFiles:
                oFile = open(sFiles, 'r')
                skip_comment(oFile, iLineHeader=2)
                lFileReader.append(oFile)
            # # peek the first scan number
            iSmallestScanNumer = sys.maxint
            for f in lFileReader:
                sLine = peek_line(f)
                iScanNumber = get_scan_number(sLine)
                if iScanNumber < iSmallestScanNumer:
                    iSmallestScanNumer = iScanNumber
            # #
            iLastScanExcluded = iSmallestScanNumer
            bReachEof = False
            while(not bReachEof):
                psmPack = PsmPack(_iSize=self.iScanInterval, _iStartScanNumber=iLastScanExcluded)
                iLastScanExcluded = iLastScanExcluded + self.iScanInterval
                bReachEof = True
                for f in lFileReader:
                    bReachEof = bReachEof & (get_psm(f, psmPack, iUpperScanNumber=iLastScanExcluded))
                if not psmPack.empty():
                    # # add to the job queue
                    for psm in psmPack:
                        self.qPsmUnprocessed.put(psm, True)
                        self.iNumScanProcessed += 1
                        if self.iNumScanProcessed % 100 == 0:
                            pass
                            # print 'Read # scans %d' % self.iNumScanProcessed
                            # print 'Queue size %f' % self.qPsmUnprocessed.qsize()

            # # close the file reader
            for f in lFileReader:
                f.close()

    def run(self):
        if not self.qPsmUnprocessed:
            sys.stderr.write("Job Queue error in PSM reading.")
        self.categorizeSpe2PepFile(self.sInputFolder)
        self.readData()
        self.qPsmUnprocessed.put(None)

# # thread class for ranking the PSM
class RankPsm(Process):

    def __init__(self, qPsmUnprocessed, qPsmProcessed, name=None, isSIP=False):
        super(RankPsm, self).__init__()
        self.name = name
        self.qPsmUnprocessed = qPsmUnprocessed
        self.qPsmProcessed = qPsmProcessed
        self.iCount = 0
        self.SIP = isSIP
        return

    def run(self):
        if not self.qPsmUnprocessed:
            sys.stderr.write("Job Queue error in PSM ranking.")
        if not self.qPsmProcessed:
            sys.stderr.write("Job Queue error in PSM ranking.")
        while True:
            psm = self.qPsmUnprocessed.get(True)
            if psm is None:
                break
            Psm_list = SelectTopRankedPsm(psm, self.SIP)
            del psm
            self.iCount += 1
            if self.iCount % 10 == 0:
                pass
                # print "Rank # scans %i" % self.iCount
            for oPsm in Psm_list:
                self.qPsmProcessed.put(oPsm, True)
        self.qPsmProcessed.put(None)
        self.qPsmUnprocessed.put(None)
        return

# # score agreement recode
def agreement(liRank):
    iIndex = 0
    for i in liRank:
        if i == 1:
            iIndex = (iIndex << 1) + 1;
        else:
            iIndex = (iIndex << 1)
    return iIndex

# # number of score agreement
def num_agreement(liRanks):
    iNum = 0
    for i in liRanks:
        if i == 1:
            iNum += 1
    return iNum

bAdditionalFeatures = True
bExtraNegativePsm = False
bRTime = True
bAdditionPepScoreAgreementNotThree = True

# # write the PSM table
def writePsm(sOutputFile, qPsmProcessed, iNumRankers, pepxml_bool = False):
    iNumTarget = 0
    iNumReverse = 0
    iNumShuffle = 0
    iNumProcessedScans = 0
    
    with open(sOutputFile, 'w') as f:
        # header
        f.write('FileName\t')
        f.write('ScanNumber\t')
        f.write('ParentCharge\t')
        f.write('MeasuredParentMass\t')
        f.write('ScanType\t')
        f.write('SearchName\t')
        f.write('IdentifiedPeptide\t')
        f.write('OriginalPeptide\t')
        f.write('CalculatedParentMass\t')
        f.write('MVH\t')
        f.write('Xcorr\t')
        f.write('WDP\t')
        f.write('ProteinNames\t')
        f.write('Diff_MVH\tDiff_Xcorr\tDiff_WDP\t')
        f.write('RetentionTime\t')
        f.write('DeltaP\n')
        
        # PSM
        while True:
            psm = qPsmProcessed.get(True)
            if psm is None:
                iNumRankers -= 1
                if iNumRankers == 0:
                    break
                else:
                    continue
            if pepxml_bool:
                f.write(psm.all_top_5_ranked_psm())
            else:
                f.write(psm.all_top_ranked_psm())

            f.write('\n')
                
            del psm
            iNumProcessedScans += 1
            if iNumProcessedScans % 100 == 0:
                pass
                # print " Processed & Saved %i Scans\r" % iNumProcessedScans


pep_iden_str = '[Peptide_Identification]'
search_name_str = 'Search_Name'
FASTA_Database_str = 'FASTA_Database'
Maximum_Missed_Cleavages_str = 'Maximum_Missed_Cleavages'
Cleave_After_Residues_str = 'Cleave_After_Residues'
