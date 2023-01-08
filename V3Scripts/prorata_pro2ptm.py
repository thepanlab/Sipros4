#!/usr/bin/python

## Import Python package modules
import sys, getopt, warnings, os, re
from datetime import datetime, date, time
from collections import namedtuple
from collections import defaultdict
import csv
import math
import numpy as np


def parse_options(argv):

    
    opts, args = getopt.getopt(argv[1:], "hw:",
                                    ["help",
                                     "working-dir"])


    # Default working dir
    working_dir = "./"

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-w workingdirectory"
            sys.exit(1)
        if option in ("-w", "--working-dir"):
            working_dir = value
            if working_dir[-1] != '/':
                working_dir = working_dir + '/'

    pro2ptm_filename_list = get_file_list_with_ext(working_dir, ".pro2ptm.txt")
    pro2ptm_filename = pro2ptm_filename_list[0]
    if (len(pro2ptm_filename_list) > 1) :
	print "too many pro2ptm.txt files\n"
	sys.exit(1)	
    
    prorata_filename_list =  get_file_list_with_ext(working_dir, ".ProRata_LabelFree_Peptide.txt")
    prorata_filename = prorata_filename_list[0]    

    if (len(prorata_filename_list) > 1) :
	print "too many pro2ptm.txt files\n"
	sys.exit(1)
	                  

    return [pro2ptm_filename, prorata_filename]



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

       # if len(file_list) == 0:
        #    print >> sys.stderr, "\nCannot open %s file(s)." % (file_ext)
            # die("Program exit!")
	 #   sys.exit(0)
        file_list = sorted(file_list)

    else:
        print >> sys.stderr, "\nCannot open working directory", working_dir
        die("Program exit!")

    return file_list

def getColumnId(sColumnNameLine, ColumnName) :
	ColumnName_list = sColumnNameLine.split("\t")
	try:
		iColumnId = ColumnName_list.index(ColumnName)
	except ValueError:
		print "can't find column "+ColumnName
		sys.exit(0)
	#print iColumnId
	return iColumnId


def GetOutputFileName(pro2ptm_filename) :
	(pro2ptmFileNameRoot, pro2ptmFileNameExt) = os.path.splitext(pro2ptm_filename) #Ext is txt
	outputFileName = pro2ptmFileNameRoot + ".prorata.txt"
	(pro2ptmFileNameRoot, pro2ptmFileNameExt) = os.path.splitext(pro2ptmFileNameRoot) #Ext is pro2ptm
	(pro2ptmFilePath, pro2ptmFileNamePureRoot)= os.path.split(pro2ptmFileNameRoot)
	return (outputFileName,pro2ptmFileNamePureRoot)

def ReadProRataFile(prorata_filename, pro2ptmFileNamePureRoot) :
	
	ProRataFile = open(prorata_filename)	
	ProRataFile_list  = ProRataFile.readlines()
	ProRataFileLength = len(ProRataFile_list)
	sTitleLine = ProRataFile_list[1]
	sTitleLine = sTitleLine.strip()
	iXIC_filename_ColumnId = getColumnId(sTitleLine, "XIC_filename")
	iidentifier_ColumnId   = getColumnId(sTitleLine, "identifier")
	ipeakHeight_ColumnId   = getColumnId(sTitleLine, "peakHeight")
	iIDs_ColumnId          = getColumnId(sTitleLine, "IDs")
	
	XIC_identifier_List = []
	peakHeight_List     = []
	allIDs_List         = []

	for i in range(2, ProRataFileLength) : 
		currentLine = ProRataFile_list[i]
		currentLine = currentLine.strip()
		if (currentLine == "") :
			continue
		currentInfo_List = currentLine.split("\t")
		sXIC_filename    = currentInfo_List[iXIC_filename_ColumnId]
		sidentifier      = currentInfo_List[iidentifier_ColumnId]
		speakHeight      = currentInfo_List[ipeakHeight_ColumnId]
		sIDs             = currentInfo_List[iIDs_ColumnId]
		sIDs_List        = sIDs.split(",")
		iIDs_List_Length = len(sIDs_List)
		currentIDList = []
		for j in range(iIDs_List_Length) :
			eachID = sIDs_List[j]
			#print pro2ptmFileNamePureRoot
			iPosNamePureRoot = eachID.find(pro2ptmFileNamePureRoot)
			if (iPosNamePureRoot != 0) :
				eachID = eachID[iPosNamePureRoot:]
			iPosRightMostDot = eachID.rfind(".")
			iPos2ndRightMostDot = eachID.rfind(".", 0, iPosRightMostDot)
			if (iPosRightMostDot == -1) or (iPos2ndRightMostDot == -1) :
				print "wrong ID: "+eachID
				sys.exit(1)
			eachID = eachID[:iPos2ndRightMostDot] + "[" + eachID[iPos2ndRightMostDot+1:iPosRightMostDot] + "]"
			currentIDList.append(eachID)	
		XIC_identifier_List.append(sXIC_filename+"@"+sidentifier)
		peakHeight_List.append(speakHeight)
		allIDs_List.append(currentIDList)
	ProRataFile.close()
	return [XIC_identifier_List, peakHeight_List, allIDs_List]


def GetPeakHeight(sPSMs, XIC_identifier_List, peakHeight_List, allIDs_List) :
    Used_List  = []
    dSumPeakHeight = 0
    iProrataListLength = len(XIC_identifier_List)
    sPSMs_List = (sPSMs[1:-1]).split(",")
    iPSMsListLength    = len(sPSMs_List) 
    for i in range(iPSMsListLength):
        currentPSM = sPSMs_List[i]
        for j in range(iProrataListLength):
            if currentPSM in allIDs_List[j]:
                currentXIC_identifier = XIC_identifier_List[j]
                if currentXIC_identifier not in Used_List :
                    Used_List.append(currentXIC_identifier)
                    dSumPeakHeight += float(peakHeight_List[j])
    if  dSumPeakHeight == 0 :
        dSumPeakHeight = 1
    return dSumPeakHeight
    




def handleAllFiles(pro2ptm_filename, prorata_filename, outputFileName, pro2ptmFileNamePureRoot) :
	
    [XIC_identifier_List, peakHeight_List, allIDs_List] = ReadProRataFile(prorata_filename, pro2ptmFileNamePureRoot)
#	print len(XIC_identifier_List)
#	for i in range(len(XIC_identifier_List)) :
#		print XIC_identifier_List[i], peakHeight_List[i], allIDs_List[i]
    pro2ptmFile = open(pro2ptm_filename)
    outputFile  = open(outputFileName, "w")

    bProteinTitleLine = True
    bPeptideTitleLine = True
    for eachLine in pro2ptmFile :
            eachLine = eachLine.strip()
            if eachLine == "" :
                continue
            if eachLine.startswith("#") :
                outputFile.write(eachLine+"\n")
            if eachLine.startswith("+") :
                if (bProteinTitleLine) :
                    bProteinTitleLine = False
                    outputFile.write("#\tModifiedPeakHeight = Total peak height of chromatographic peaks of modified peptides\n") 
                    outputFile.write("#\tUnmodifiedPeakHeight = Total peak height of chromatographic peaks of unmodified peptides\n")
                    outputFile.write("#\n")
                outputFile.write(eachLine+"\n")
            if eachLine.startswith("*") :
                if (bPeptideTitleLine) :
                    bPeptideTitleLine = False
                    iModifiedPSMs_ColumnId   = getColumnId(eachLine, "ModifiedPSMs")
                    iUnmodifiedPSMs_ColumnId = getColumnId(eachLine, "UnmodifiedPSMs")
                    outputFile.write(eachLine+"\tModifiedPeakHeight\tUnmodifiedPeakHeight\n")
                else :
                    eachLineInfo_List = eachLine.split("\t")
                    sModifiedPSMs     = eachLineInfo_List[iModifiedPSMs_ColumnId]
                    sUnmodifiedPSMs   = eachLineInfo_List[iUnmodifiedPSMs_ColumnId]
                    ModifiedPeakHeight= GetPeakHeight(sModifiedPSMs, XIC_identifier_List, peakHeight_List, allIDs_List)
                    UnmodifiedPeakHeight= GetPeakHeight(sUnmodifiedPSMs,XIC_identifier_List,peakHeight_List,allIDs_List)
                    outputFile.write(eachLine+"\t"+str(ModifiedPeakHeight)+"\t"+str(UnmodifiedPeakHeight)+"\n")
    pro2ptmFile.close()
    outputFile.close()

## +------+
## | Main |
## +------+
def main(argv=None):

    # try to get arguments and error handling
        if argv is None:
		argv = sys.argv
       		 # parse options
		[pro2ptm_filename, prorata_filename] = parse_options(argv)
	(outputFileName,pro2ptmFileNamePureRoot) = GetOutputFileName(pro2ptm_filename)
	handleAllFiles(pro2ptm_filename, prorata_filename, outputFileName, pro2ptmFileNamePureRoot)


## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()
