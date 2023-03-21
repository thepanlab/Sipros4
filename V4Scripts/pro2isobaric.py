#!/usr/bin/python


## Import Python package modules
import sys, getopt, warnings, os, re
from datetime import datetime, date, time
from collections import namedtuple
from collections import defaultdict
import csv
import math
import numpy as np


## Import Sipros package modules
#import sipros_post_module
import parseconfig
#import HierarchicalClustering

## Parse options
def parse_options(argv):

    
    opts, args = getopt.getopt(argv[1:], "hw:c:",
                                    ["help",
                                     "working-dir",
                                     "config-file"])


    # Default working dir and config file
    working_dir = "./"
    config_file = "SiprosConfig.cfg"

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-c configurefile -w workingdirectory"
            sys.exit(1)
        if option in ("-w", "--working-dir"):
            working_dir = value
            if working_dir[-1] != '/':
                working_dir = working_dir + '/'
        if option in ("-c", "--config-file"):
            config_file = value

    # only -w is provided
    if working_dir != "./" and config_file == "SiprosConfig.cfg":
        config_file = working_dir + config_file

    psm_filename_list = get_file_list_with_ext(working_dir, ".psm.isobaric.txt")
    psm_filename = psm_filename_list[0]
    
    pro2psm_filename_list =  get_file_list_with_ext(working_dir, ".pro2psm.txt")
    pro2psm_filename = pro2psm_filename_list[0]                      

    return [config_file, psm_filename, pro2psm_filename]



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

def ParseFileName(sFile_Name) : 
# return real file name
	(sHeadName, sTailName) = os.path.split(sFile_Name)
	return sTailName


def ReadPsmFile(psm_filename, iReporter_Ion_Number) :
# extract reporter ion value from the .psm.isobaric.txt file
	Reporter_Ion_Info_dict = {}
	sReporter_Ion_Names_list = []
	bFirstLine = True
	PsmFile = open(psm_filename)
	for eachLine in PsmFile :
		eachLine = eachLine.strip()
		if (eachLine.startswith("#") or (eachLine == "")) :
			continue
		lineInfo = eachLine.split("\t")
		if (bFirstLine) :
			iPsm_FileName_ColumnId = getColumnId(eachLine, "Filename")
			iPsm_ScanId_ColumnId = getColumnId(eachLine, "ScanNumber")
			bFirstLine = False
			for i in range((-1)*iReporter_Ion_Number, 0) :
				sReporter_Ion_Names_list.append(lineInfo[i])
		else :
			sFull_Scan_Filename = lineInfo[iPsm_FileName_ColumnId]
			sTail_Scan_Filename = ParseFileName(sFull_Scan_Filename)
			ScanId = lineInfo[iPsm_ScanId_ColumnId]
			dReporter_Ion_Values_list = []
			for i in range((-1)*iReporter_Ion_Number, 0) :
				 #print i, iReporter_Ion_Number, lineInfo[i]
				 dReporter_Ion_Values_list.append(float(lineInfo[i]))
			Reporter_Ion_Info_dict[ScanId+"@"+sTail_Scan_Filename] = dReporter_Ion_Values_list
			#print ScanId+"@"+sTail_Scan_Filename, dReporter_Ion_Values_list

	PsmFile.close()
	return [Reporter_Ion_Info_dict, sReporter_Ion_Names_list]


def ProteinNameListCompare(stringA, stringB):
	stringA = stringA.lstrip("{")
	stringA = stringA.rstrip("}")
	stringA_list = stringA.split(",")
	
	stringB = stringB.lstrip("{")
	stringB = stringB.rstrip("}")
	stringB_list = stringB.split(",")
	
	bComp = (set(stringA_list) == set(stringB_list))
	return bComp

def HandlePro2psmFile(pro2psm_filename, ScansInfo_dict, sReporter_Ion_Names_list, sOnly_Use_Unique_Peptides, sPeptide_Modification_Symbol) :
	
	sExtended_Title_line = ""
	for sIon_Name in sReporter_Ion_Names_list :
		sExtended_Title_line = sExtended_Title_line + "\t" + sIon_Name

	Pro2psmFile = open(pro2psm_filename)	
	sCoreFilename = pro2psm_filename[:-12]
	sPro2psm_Isobaric_Filename = sCoreFilename + ".pro2psm.isobaric.txt"
	sPro_Isobaric_Filename = sCoreFilename + ".pro.isobaric.txt"
	Pro2psm_Isobaric_File = file(sPro2psm_Isobaric_Filename, "w")
	Pro_Isobaric_File     = file(sPro_Isobaric_Filename, "w")
	bProTitleLine = True
	bPsmTitleLine = True
	sPsmUnit_list = [] # first is protein, others are scans
	dSummation_ion_values_list = [0 for x in range(len(sReporter_Ion_Names_list))]
	for eachLine in Pro2psmFile :
		eachLine = eachLine.strip()
		if (eachLine == "") :
			continue
		if eachLine.startswith("#") :
			Pro2psm_Isobaric_File.write(eachLine+"\n")
			Pro_Isobaric_File.write(eachLine+"\n")
		if eachLine.startswith("+") : # protein line
			if (bProTitleLine) : 
				bProTitleLine = False 
				iProteinID_ColumnId  = getColumnId(eachLine, "ProteinID")
				Pro2psm_Isobaric_File.write(eachLine + sExtended_Title_line  + "\n")
				Pro_Isobaric_File.write((eachLine + sExtended_Title_line)[2:]  + "\n")
			else :
				if (len(sPsmUnit_list) > 0  ) :
					sProteinLine = sPsmUnit_list[0]
					for eachValue in dSummation_ion_values_list :
						sProteinLine = sProteinLine + "\t" + str(eachValue)
					Pro2psm_Isobaric_File.write(sProteinLine + "\n")
					Pro_Isobaric_File.write(sProteinLine[2:] + "\n")
					for i in range(1, len(sPsmUnit_list)):
						Pro2psm_Isobaric_File.write(sPsmUnit_list[i] + "\n")
					sPsmUnit_list = []
					dSummation_ion_values_list = [0 for x in range(len(sReporter_Ion_Names_list))]
				sProteinInfo_list = eachLine.split("\t")
				sProteinID   =  sProteinInfo_list[iProteinID_ColumnId]
				sPsmUnit_list.append(eachLine)

		if eachLine.startswith("*") : # psm line
			if (bPsmTitleLine) :
				iProPsm_Scan_File_Fullname_ColumnId = getColumnId(eachLine, "Filename")
				iProPsm_ScanName_ColumnId = getColumnId(eachLine, "ScanNumber")
				iProPsm_ProteinCount_ColumnId = getColumnId(eachLine, "ProteinCount")
				iProPsm_IdentifiedPeptide_ColumnId = getColumnId(eachLine, "IdentifiedPeptide")
				iProPsm_ProteinNames_ColumnId = getColumnId(eachLine, "ProteinNames")
				bPsmTitleLine = False
				Pro2psm_Isobaric_File.write(eachLine + sExtended_Title_line  + "\n")
				#Pro_Isobaric_File.write(eachLine + sExtended_Title_line  + "\n")
			else :
				sPsmInfo_list = eachLine.split("\t")
				sScan_File_Fullname = sPsmInfo_list[iProPsm_Scan_File_Fullname_ColumnId]
				sScan_File_Tailname = ParseFileName(sScan_File_Fullname)
				sScanName = sPsmInfo_list[iProPsm_ScanName_ColumnId]
				sProteinCount = sPsmInfo_list[iProPsm_ProteinCount_ColumnId]
				sIdentifiedPeptide = sPsmInfo_list[iProPsm_IdentifiedPeptide_ColumnId]
				dReporter_Ion_Values_list = ScansInfo_dict.get(sScanName+"@"+sScan_File_Tailname)
				sCurrentProteinName = sPsmInfo_list[iProPsm_ProteinNames_ColumnId]
				if (dReporter_Ion_Values_list == None) :
					print "Can't find "+ sScanName+"@"+sScan_File_Tailname
					sys.exit(0)
				for eachValue in dReporter_Ion_Values_list :
					eachLine = eachLine + "\t" + str(eachValue)
				sPsmUnit_list.append(eachLine)
				if ((sOnly_Use_Unique_Peptides == "True") and (sProteinCount != "1")) :
					if not (ProteinNameListCompare(sCurrentProteinName, sProteinID)) :
						continue
				bPass_Peptide_Modification_Symbol = True
				if (sPeptide_Modification_Symbol != "None"):
					for i in range(len(sPeptide_Modification_Symbol)) :
						if not( sPeptide_Modification_Symbol[i:i+1] in sIdentifiedPeptide ) :
							bPass_Peptide_Modification_Symbol = False
							#print sIdentifiedPeptide
							break
				if (bPass_Peptide_Modification_Symbol ) :
					for i in range(len(dSummation_ion_values_list)) :
						dSummation_ion_values_list[i] += dReporter_Ion_Values_list[i]
	if (len(sPsmUnit_list) > 0  ) :
		sProteinLine = sPsmUnit_list[0]
		for eachValue in dSummation_ion_values_list :
			sProteinLine = sProteinLine + "\t" + str(eachValue)
		Pro2psm_Isobaric_File.write(sProteinLine + "\n")
		Pro_Isobaric_File.write((sProteinLine)[2:] + "\n")
		for i in range(1, len(sPsmUnit_list)):
			Pro2psm_Isobaric_File.write(sPsmUnit_list[i] + "\n")

	Pro2psm_Isobaric_File.close()
	Pro_Isobaric_File.close()
	Pro2psmFile.close()



## +------+
## | Main |
## +------+
def main(argv=None):

    # try to get arguments and error handling
        if argv is None:
		argv = sys.argv
       		 # parse options
		[config_filename, psm_filename, pro2psm_filename] = parse_options(argv)
	wholeDict = parseconfig.parseConfigKeyValues(config_filename)
	sOnly_Use_Unique_Peptides = wholeDict.get("[Isobaric_Chemical_Labeling]Only_Use_Unique_Peptides")
	sReporter_Ion_Key_dict = parseconfig.getConfigMasterKeyValue ("[Isobaric_Chemical_Labeling]Reporter_Ion", wholeDict) 
	sPeptide_Modification_Symbol = wholeDict.get("[Isobaric_Chemical_Labeling]Peptide_Modification_Symbol")

	[ScansInfo_dict, sReporter_Ion_Names_list] = ReadPsmFile(psm_filename, len(sReporter_Ion_Key_dict))
	HandlePro2psmFile(pro2psm_filename, ScansInfo_dict, sReporter_Ion_Names_list, sOnly_Use_Unique_Peptides, sPeptide_Modification_Symbol)


## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()
