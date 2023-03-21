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
    psm_filename_list = get_file_list_with_ext(working_dir, ".psm.txt")
    psm_filename = psm_filename_list[0]

    FT2_filename_list = get_file_list_with_ext(working_dir, ".FT2")
    MS2_filename_list = get_file_list_with_ext(working_dir, ".MS2")
    Scans_filename_list = FT2_filename_list + MS2_filename_list 
    
    #pro2psm_filename_list =  get_file_list_with_ext(working_dir, ".pro2psm.txt")
    #pro2psm_filename = pro2psm_filename_list[0]                      

    return [config_file, psm_filename, Scans_filename_list]



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

def ParseTypeIdentification(sIdentification_type) :
	#index of @
	iAt_index = sIdentification_type.find("@")
#	print sIdentification_type
	if (iAt_index < 0) :
		print "no @ at "+sIdentification_type
		sys.exit(0)
	else :
		sSeekInfo = sIdentification_type[iAt_index+1 :]
		sSeekInfo = sSeekInfo.strip() # remove white spaces
		sSeekInfo = sSeekInfo.strip("]") #remove "]"
#		print sSeekInfo
		sSeekInfo_list = sSeekInfo.split("[")
		SeekInfo_list = []
		if (len(sSeekInfo_list) != 2) :
			print "incorrect Identification_Scan "+sIdentification_type
			sys.exit(0)
		else :
#			print sSeekInfo_list
			SeekInfo_list.append(sSeekInfo_list[0])
			SeekInfo_list.append(int(sSeekInfo_list[1]))
	return SeekInfo_list

def OrganizeScanInfo(FT2Scans_list, HCD_scan_seek_list, CID_scan_seek_list, Reporter_Ion_Key_list) :
	Scan_dict = {}
	Zero_Intensities_list = [ 0 for x in range(len(Reporter_Ion_Key_list))]
	
	for i in range(len(FT2Scans_list)):
		sScan_Type = FT2Scans_list[i][1]
		if (sScan_Type == "HCD") :
			sTarget_Scan_Type = HCD_scan_seek_list[0]
			iTarget_Scan_Index= i + HCD_scan_seek_list[1]
		elif (sScan_Type == "CID") :
			sTarget_Scan_Type = CID_scan_seek_list[0]
			iTarget_Scan_Index= i + CID_scan_seek_list[1]
		else :
			print "incorrect scan type: "+sScan_Type
			sys.exit(0)
		if ((iTarget_Scan_Index >=0) and (iTarget_Scan_Index < len(FT2Scans_list)) and (sTarget_Scan_Type == FT2Scans_list[iTarget_Scan_Index][1])) :
			Scan_dict[FT2Scans_list[i][0]] = [FT2Scans_list[iTarget_Scan_Index][2], FT2Scans_list[iTarget_Scan_Index][0]]
#			print i, iTarget_Scan_Index, FT2Scans_list[i][0]
		else :
#			print i, iTarget_Scan_Index, FT2Scans_list[i][0], sTarget_Scan_Type
			Scan_dict[FT2Scans_list[i][0]] = [Zero_Intensities_list, "Zero Case"]
#	print CID_scan_seek_list, HCD_scan_seek_list
	return Scan_dict


def ParseFileName(sFile_Name) : 
# return real file name
	(sHeadName, sTailName) = os.path.split(sFile_Name)
	return sTailName

    
def ReadFT2File(FT2_filename, HCD_scan_seek_list, CID_scan_seek_list, Reporter_Ion_Key_list) :
# save key information of the FT2 file in a dictionary
# this function put key information to a list first, then save organized information in a dictionary
	FT2File = open(FT2_filename) 
	bFirstScan = True
	FT2Scans_list = []	
	sFT2File_TailName = ParseFileName(FT2_filename)

	for eachline in FT2File.readlines() :
		eachline = eachline.strip()
		if (eachline == "") :
			continue
		if eachline.startswith("S") :
			if (bFirstScan) :
				bFirstScan = False
			else :
				FT2Scans_list.append([sScanId+"@"+sFT2File_TailName, sScan_Type ,dHighest_Intensities_list])
			lineInfo = eachline.split("\t")
			sScanId = lineInfo[1] # scanid is a string
			# initialize the list of highest intensities 
			dHighest_Intensities_list = [ 0 for x in range(len(Reporter_Ion_Key_list))]
			sScan_Type = ""
		elif eachline.startswith("I") :
			lineInfo = eachline.split("\t")
			if (lineInfo[1] == "ScanType"):
				Scan_Type_Info = lineInfo[2].split("@")
				sScan_Type = Scan_Type_Info[1].strip()
		elif eachline[0:1].isdigit() :  # starts with digits
			lineInfo = eachline.split("\t")
			if (len(lineInfo) < 2) :
				continue
			dMassofZ = float(lineInfo[0])
			dIntensity = float(lineInfo[1])  
			for i in range(len(Reporter_Ion_Key_list)) :
				if ((Reporter_Ion_Key_list[i][1] <= dMassofZ) and (dMassofZ <= Reporter_Ion_Key_list[i][2]) and (dIntensity > dHighest_Intensities_list[i])):
					dHighest_Intensities_list[i] = dIntensity
	#last scan
	if not(bFirstScan):
		FT2Scans_list.append([sScanId+"@"+sFT2File_TailName, sScan_Type, dHighest_Intensities_list])
	FT2File.close()
	FT2Scans_dict = OrganizeScanInfo(FT2Scans_list, HCD_scan_seek_list, CID_scan_seek_list, Reporter_Ion_Key_list)
	return FT2Scans_dict


def ReadAllFT2Files(FT2_filename_list, sIdentification_HCD, sIdentification_CID, sMass_Tolerance_Reporter_Ions, sReporter_Ion_Key_dict) :  
	print "Reading FT2 Files ..." 
	# how to seek the suitable HCD scan, 0th item is about scan type, 1st item is about index offset
	HCD_scan_seek_list = ParseTypeIdentification(sIdentification_HCD)
	# how to seek the suitable CID scan, 0th item is about scan type, 1st item is about index offset
	if (sIdentification_CID != ""):
		CID_scan_seek_list = ParseTypeIdentification(sIdentification_CID)
	else :
		CID_scan_seek_list = ["", 0]
	dMass_Tolerance_Reporter_Ions = float (sMass_Tolerance_Reporter_Ions)
	Reporter_Ion_Key_list = []
	for k, v in sReporter_Ion_Key_dict.iteritems():
		# k is Reporter_Ion, v is MassofZ
		Reporter_Ion_Key_list.append([k, float(v)-dMass_Tolerance_Reporter_Ions, float(v)+dMass_Tolerance_Reporter_Ions])
	Reporter_Ion_Key_list.sort(key = lambda x:x[0])
	All_FT2Scans_dict = {}	
	for each_FT2_filename in FT2_filename_list :
		Current_FT2Scans_dict =  ReadFT2File(each_FT2_filename, HCD_scan_seek_list, CID_scan_seek_list, Reporter_Ion_Key_list)
		All_FT2Scans_dict.update(Current_FT2Scans_dict)

#	print Reporter_Ion_Key_list
	#for k, v in All_FT2Scans_dict.iteritems():
		#k is scanid and filename, v is values
	#	print k,v

	return [All_FT2Scans_dict, Reporter_Ion_Key_list]

def NormalizaReporterIonValue(Summation_list, Normalization_list) :
	#calculate Normalization factors
	#print "Normalizing Reporter Ion Values ..."
	Summation_of_All_ReporterIonValue = sum(Summation_list)
	for i in range(len(Normalization_list)) :
		if (Summation_list[i] > 0) :
			Normalization_list[i] = Summation_of_All_ReporterIonValue/(len(Normalization_list) * Summation_list[i])
	


def HandlePsmFile(sPsm_filename, FT2Info_dict, Reporter_Ion_list, sTotal_Intensity_Normalization, sPeptide_Modification_Symbol) :
	print "handling psm file ..."
	PsmFile = open(sPsm_filename)
	sPrefix_Filename = sPsm_filename[:-8]
	sPsm_Isobaric_Filename = sPrefix_Filename + ".psm.isobaric.txt"
	Psm_Isobaric_File = file(sPsm_Isobaric_Filename, "w")
	
	bPsmFirstLine = True # True: hasn't read the title line of the psm file
	Zero_Intensities_list = [ 0 for x in range(len(Reporter_Ion_list))]
	All_Original_Lines_list = []
	All_Used_ScanInfo_list  = []
	Summation_list = [0 for x in range(len(Reporter_Ion_list))]
	Normalization_list = [1 for x in range(len(Reporter_Ion_list))] #factors for normalizing Reporter_Ion_Value
	All_Reporter_Ion_Value_list = []

	for eachLine in PsmFile :
		eachLine = eachLine.strip()
		if eachLine.startswith("#") :
			Psm_Isobaric_File.write(eachLine+"\n")
			continue
		if (eachLine == "") :
			continue
		elif (bPsmFirstLine) :
			iFT2_File_Fullname_ColumnId = getColumnId(eachLine, "Filename")
			iScanNum_ColumnId = getColumnId(eachLine, "ScanNumber")
			iTargetMatch_ColumnId = getColumnId(eachLine, "TargetMatch")
			iIdentifiedPeptide_ColumnId = getColumnId(eachLine, "IdentifiedPeptide")
			bPsmFirstLine = False
			for each_Key_Element in Reporter_Ion_list :
				eachLine = eachLine + "\tReporter_Ion{" + each_Key_Element[0] +"}"
			Psm_Isobaric_File.write(eachLine+"\n")
		else :
			Protein_Info_list = eachLine.split("\t")
			sProtein_Fullname = Protein_Info_list[iFT2_File_Fullname_ColumnId] #sProtein_Fullname is the FT2 file full name
			sProtein_Tailname = ParseFileName(sProtein_Fullname)
			sScanNumber       = Protein_Info_list[iScanNum_ColumnId]
			sTargetMatch      = Protein_Info_list[iTargetMatch_ColumnId]
			sIdentifiedPeptide= Protein_Info_list[iIdentifiedPeptide_ColumnId]
			Reporter_Ion_AllInfo_list = FT2Info_dict.get(sScanNumber+"@"+sProtein_Tailname)
			if (Reporter_Ion_AllInfo_list == None) :
				print "can't find "+sScanNumber+"@"+sProtein_Tailname
				sys.exit(0)
			Reporter_Ion_Value_list = Reporter_Ion_AllInfo_list[0]
			Reporter_Ion_Original_ScanInfo = Reporter_Ion_AllInfo_list[1]
			if (Reporter_Ion_Original_ScanInfo in All_Used_ScanInfo_list) :
				#print "used scan: "+Reporter_Ion_Original_ScanInfo
				All_Used_ScanInfo_list.append(Reporter_Ion_Original_ScanInfo)
				All_Original_Lines_list.append(eachLine)
				All_Reporter_Ion_Value_list.append(Zero_Intensities_list)
				#for i in range(len(Reporter_Ion_Value_list)) :
				#	Summation_list[i] += Zero_Intensities_list[i]
			elif ((sTargetMatch == "T") and (sPeptide_Modification_Symbol in sIdentifiedPeptide)):
				All_Used_ScanInfo_list.append(Reporter_Ion_Original_ScanInfo)
				All_Original_Lines_list.append(eachLine)
				All_Reporter_Ion_Value_list.append(Reporter_Ion_Value_list)
				for i in range(len(Reporter_Ion_Value_list)) :
					Summation_list[i] += Reporter_Ion_Value_list[i]
			else :
				All_Original_Lines_list.append(eachLine)
				All_Reporter_Ion_Value_list.append(Zero_Intensities_list)
	if (sTotal_Intensity_Normalization == "True") :
		NormalizaReporterIonValue(Summation_list, Normalization_list)
	

	for i in range(len(All_Original_Lines_list)) :
		eachLine = All_Original_Lines_list[i]
		for j in range(len(Reporter_Ion_Value_list)) :
			eachValue = All_Reporter_Ion_Value_list[i][j] * Normalization_list[j]
			eachLine = eachLine + "\t" + str(eachValue)
		Psm_Isobaric_File.write(eachLine+"\n")

	
	PsmFile.close()
	Psm_Isobaric_File.close()


def OrganizeReporterIon(Reporter_Ion_Key_list, sReporter_Ion_Key_dict) :
	Reporter_Ion_list = []
	for each_Key_Element in  Reporter_Ion_Key_list :
		current_Key = each_Key_Element[0]
		current_Key_Value = sReporter_Ion_Key_dict.get(current_Key)
		Reporter_Ion_list.append([current_Key, int(round(float(current_Key_Value)))])

	#Reporter_Ion_list.sort(key = lambda x:x[0])
	return Reporter_Ion_list

def IsotopicImpurityCorrect(All_FT2Scans_dict, Reporter_Ion_list, sReporter_Ion_Isotopic_Impurity_Distribution_Key_dict) :
	print "Correcting Isotopic Impurity ..."
# equation set Ax=B
	a  = []
	Zero_list = [0 for x in range(len(Reporter_Ion_list))]
	selfIndex = 2 # index for the current Isotopic_Impurity_Distribution in the Reporter_Ion_Isotopic_Impurity_Distribution_list
	for i in range(len(Reporter_Ion_list)) :
		current_Reporter_Ion = Reporter_Ion_list[i][0]
		sReporter_Ion_Isotopic_Impurity_Distribution = sReporter_Ion_Isotopic_Impurity_Distribution_Key_dict.get(current_Reporter_Ion)
		if (sReporter_Ion_Isotopic_Impurity_Distribution == None) :
			print "cant find reporter ion isotopic impurity distribution of " + sReporter_Ion_Isotopic_Impurity_Distribution
		else :
			sReporter_Ion_Isotopic_Impurity_Distribution = sReporter_Ion_Isotopic_Impurity_Distribution.strip()
			sReporter_Ion_Isotopic_Impurity_Distribution_list =  sReporter_Ion_Isotopic_Impurity_Distribution.split(",")
			a_row = []
			for j in range(len(Reporter_Ion_list)) :
				if (j==i) :
					currentIndex = selfIndex
				else :
					currentIndex = selfIndex + Reporter_Ion_list[j][1] - Reporter_Ion_list[i][1]
				if ((currentIndex >=0) and (currentIndex <=4)) :
					a_row.append(float(sReporter_Ion_Isotopic_Impurity_Distribution_list[currentIndex])/100.0)
				else :
					a_row.append(0)
		a.append(a_row)
#	print a
	A = np.transpose(np.array(a))
	for k, v in All_FT2Scans_dict.iteritems():
		#k is scanid and filename, v is Reporter_Ion_Isotopic_Impurity
		if (v[0] != Zero_list) :
			B= np.array(v[0])
			try :
				X= np.linalg.solve(A, B)
			except np.linalg.linalg.LinAlgError as err:
				if 'Singular matrix' in err.message:
					print "Singular matrix! No isotopic impurity correction is conducted."
					return

			All_FT2Scans_dict[k] = [X, v[1]]
		





## +------+
## | Main |
## +------+
def main(argv=None):

    # try to get arguments and error handling
        if argv is None:
		argv = sys.argv
       		 # parse options
		[config_filename, psm_filename, FT2_filename_list] = parse_options(argv)
	wholeDict = parseconfig.parseConfigKeyValues(config_filename)

	#sDecoy_Prefix = wholeDict.get("[Protein_Identification]Decoy_Prefix")
	sIdentification_HCD  = wholeDict.get("[Isobaric_Chemical_Labeling]Identification_Scan{HCD}")
	if ("[Isobaric_Chemical_Labeling]Identification_Scan{CID}" in wholeDict) :
		sIdentification_CID  = wholeDict.get("[Isobaric_Chemical_Labeling]Identification_Scan{CID}")
	sIdentification_CID  = ""
	sMass_Tolerance_Reporter_Ions = wholeDict.get("[Isobaric_Chemical_Labeling]Mass_Tolerance_Reporter_Ions")
	sReporter_Ion_Key_dict = parseconfig.getConfigMasterKeyValue ("[Isobaric_Chemical_Labeling]Reporter_Ion", wholeDict) 
	sIsotopic_Impurity_Correction = wholeDict.get("[Isobaric_Chemical_Labeling]Isotopic_Impurity_Correction")
	sReporter_Ion_Isotopic_Impurity_Distribution_Key_dict = {}
	if (sIsotopic_Impurity_Correction == "True") :
		sReporter_Ion_Isotopic_Impurity_Distribution_Key_dict = parseconfig.getConfigMasterKeyValue ("[Isobaric_Chemical_Labeling]Reporter_Ion_Isotopic_Impurity_Distribution", wholeDict)
#		print sReporter_Ion_Isotopic_Impurity_Distribution_Key_dict
	sTotal_Intensity_Normalization = wholeDict.get("[Isobaric_Chemical_Labeling]Total_Intensity_Normalization")
	sOnly_Use_Unique_Peptides = wholeDict.get("[Isobaric_Chemical_Labeling]Only_Use_Unique_Peptides")
	sPeptide_Modification_Symbol = wholeDict.get("[Isobaric_Chemical_Labeling]Peptide_Modification_Symbol")
	if (sPeptide_Modification_Symbol == "None") :
		sPeptide_Modification_Symbol = ""
	
	[All_FT2Scans_dict, Reporter_Ion_Key_list]  =  ReadAllFT2Files(FT2_filename_list, sIdentification_HCD, sIdentification_CID, sMass_Tolerance_Reporter_Ions, sReporter_Ion_Key_dict)
	Reporter_Ion_list = OrganizeReporterIon(Reporter_Ion_Key_list, sReporter_Ion_Key_dict)
	if (sIsotopic_Impurity_Correction == "True") :
		IsotopicImpurityCorrect(All_FT2Scans_dict, Reporter_Ion_list, sReporter_Ion_Isotopic_Impurity_Distribution_Key_dict)

	HandlePsmFile(psm_filename, All_FT2Scans_dict, Reporter_Ion_list, sTotal_Intensity_Normalization, sPeptide_Modification_Symbol)
	


		


## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()

