#!/usr/bin/python


## Import Python package modules
import sys, getopt, warnings, os, re
from datetime import datetime, date, time
from collections import namedtuple
from collections import defaultdict
import csv
import math



## Import Sipros package modules
#import sipros_post_module
import parseconfig
import HierarchicalClustering

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
    psm_filename_list = get_file_list_with_ext(working_dir, ".pro2psm.txt")
    psm_filename = psm_filename_list[0]                        

    return (config_file, psm_filename)








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

def getColumnId(sColumnNameLine, ColumnName) :
	ColumnName_list = sColumnNameLine.split("\t")
	try:
		iColumnId = ColumnName_list.index(ColumnName)
	except ValueError:
		print "can't find column "+ColumnName
		sys.exit(0)
	#print iColumnId
	return iColumnId
    
def clusterSpectrum(proteinData, iProteinId_ColumnId, iProteinDescription_ColumnId, iSearchName_ColumnId, iProteinName_ColumId) :
	proteinLine = proteinData[0]
	proteinLineFeatures = proteinLine.split("\t")
	sProteinID = proteinLineFeatures[iProteinId_ColumnId]
	sProteinDescription = proteinLineFeatures[iProteinDescription_ColumnId]
	spectrumLineSet = {}
	allEnrichments = []
	allProteinLines = []
	allSpectrumLines = []
	for i in range(1, len(proteinData)) :
		spectrumFeatures = proteinData[i].split("\t")
		sSearchName = spectrumFeatures[iSearchName_ColumnId]
		SearchNameInfo = sSearchName.split("_")
		currentEnrichment = int(SearchNameInfo[1][:-3])
		newLineList = [proteinData[i]]
		if (currentEnrichment in spectrumLineSet) :
			oldLineList = spectrumLineSet[currentEnrichment]
			spectrumLineSet[currentEnrichment] = oldLineList + newLineList
		else :			
			spectrumLineSet[currentEnrichment] = newLineList
		allEnrichments.append(currentEnrichment)
	#print allEnrichments
	#print dClusteringThreshold
	#sys.exit(1)
	spectrumClusters = HierarchicalClustering.hierarchicalClustering(allEnrichments, dClusteringThreshold)
	for eachSpectrumCluster in spectrumClusters :
		if (len(eachSpectrumCluster) < iMinPSMPerIsotopicCluster) :
			continue
		iTotalSpectrumCounts = len(eachSpectrumCluster)
		iUniqueSpectrumCounts = 0
		sEnrichmentLevels = "{"
		dEnrichmentSummation = 0
		dEnrichmentSequreSummation = 0
		SpectrumClusterLines = []
		SpectrumUsed         = []
		for eachSpectrum in eachSpectrumCluster :
			if (sEnrichmentLevels == "{" ) :
				sEnrichmentLevels = sEnrichmentLevels + str(eachSpectrum)
			else :
				sEnrichmentLevels = sEnrichmentLevels + "," + str(eachSpectrum)
			dEnrichmentSummation += eachSpectrum
			dEnrichmentSequreSummation += eachSpectrum*eachSpectrum
			if (eachSpectrum not in SpectrumUsed) :
				SpectrumUsed.append(eachSpectrum)
				spectrumLine  = spectrumLineSet.get(eachSpectrum) # spectrumLine is a list now
				SpectrumClusterLines = SpectrumClusterLines + spectrumLine
				for eachEle in spectrumLine :
					currentrumInfo = eachEle.split("\t")
					sProteinNames  = currentrumInfo[iProteinName_ColumId]
					if not(","  in  sProteinNames) :
						iUniqueSpectrumCounts += 1
		sEnrichmentLevels = sEnrichmentLevels + "}"
		dAverageEnrichmentLevel = float(dEnrichmentSummation) / len(eachSpectrumCluster)
		iClusterSize = len(eachSpectrumCluster)
		if (iClusterSize == 1):
			dStandardDeviation = 0
		else :
			dStandardDeviation = math.sqrt(float(float(dEnrichmentSequreSummation)/iClusterSize - dAverageEnrichmentLevel*dAverageEnrichmentLevel )*iClusterSize/(iClusterSize-1))
		outputProteinLine = "+\t"+sProteinID+"\t"+str(dAverageEnrichmentLevel)+"\t"+str(dStandardDeviation)+"\t"+str(iUniqueSpectrumCounts)
		outputProteinLine = outputProteinLine + "\t" + str(iTotalSpectrumCounts) + "\t" +sEnrichmentLevels+"\t"+sProteinDescription
		allProteinLines.append(outputProteinLine)
		allSpectrumLines.append(SpectrumClusterLines)
	return (allProteinLines, allSpectrumLines)

def RemoveOutsideBrackets(wholeName) :
	currentName = wholeName
	currentName = currentName.strip("{}")
	return currentName

def EqualProteinName(sName1, sName2, separateStr) :
	bEqual = False
	sCurrentName1 = RemoveOutsideBrackets(sName1)
	sCurrentName2 = RemoveOutsideBrackets(sName2)
	sCurrentName1Info = sCurrentName1.split(separateStr)
	sCurrentName1Info.sort()
	sCurrentName2Info = sCurrentName2.split(separateStr)
	sCurrentName2Info.sort()
	allNameStr1 = ""
	for eachName in sCurrentName1Info:
		allNameStr1 = allNameStr1 + "!" + eachName
	allNameStr2 = ""
	for eachName in sCurrentName2Info:
		allNameStr2 = allNameStr2 + "!" + eachName
	if (allNameStr1 == allNameStr2) :
		bEqual = True
	return bEqual

	
def MergeData(allProteinLineData, allSpectrumLineData, iProteinID_Output_ColumnId, iUniqueSpectrumCounts_Output_ColumnId,   iProteinDescription_Output_ColumnId, iProteinName_ColumId, iFilename_ColumnId, iScanNumber_ColumnId) :	
	# preprocessing:
	allSpectrumStringLists = []
	for eachProtein in allSpectrumLineData :
		allSpectrumStringSubList = []
		for eachProteinCluster in eachProtein :
			allCurrentSpectrums = []
			for eachSpectrumLine in eachProteinCluster :
				eachSpectrumInfo = eachSpectrumLine.split("\t")
				allCurrentSpectrums.append(eachSpectrumInfo[iScanNumber_ColumnId]+"@"+eachSpectrumInfo[iFilename_ColumnId])
			allCurrentSpectrums.sort()
			sCurrentCombineString = ""
			for eachCombine in allCurrentSpectrums :
				sCurrentCombineString = sCurrentCombineString + "!" + eachCombine
			allSpectrumStringSubList.append(sCurrentCombineString)
		allSpectrumStringLists.append(allSpectrumStringSubList)
	#print len(allProteinLineData)
	#print len(allSpectrumStringLists)
	for i in range(len(allSpectrumStringLists)) :
		#print len(allProteinLineData[i])
		for j in range(len(allSpectrumStringLists[i])):
			#print i, j
			sProteinLine = allProteinLineData[i][j]
			allCurrentProteinInfo = sProteinLine.split("\t")
			sProteinID   = RemoveOutsideBrackets(allCurrentProteinInfo[iProteinID_Output_ColumnId])
			sProteinDescription = RemoveOutsideBrackets(allCurrentProteinInfo[iProteinDescription_Output_ColumnId])
			sThisSpectrumCombineString = allSpectrumStringLists[i][j]
			for k in range(i+1, len(allSpectrumStringLists)) :
				iLengthSubList = len(allSpectrumStringLists[k])
				for r in reversed(range(iLengthSubList)) :
					if (allSpectrumStringLists[k][r] == sThisSpectrumCombineString) :
						sThatProteinLine = allProteinLineData[k][r]
						allThatProteinInfo = sThatProteinLine.split("\t")
						sThatProteinID   = RemoveOutsideBrackets(allThatProteinInfo[iProteinID_Output_ColumnId])
						sThatProteinDescription = RemoveOutsideBrackets(allThatProteinInfo[iProteinDescription_Output_ColumnId])
						sProteinID = sProteinID + "," + sThatProteinID
						sProteinDescription = sProteinDescription + "," + sThatProteinDescription
						del allSpectrumStringLists[k][r]
						del allProteinLineData[k][r]
						del allSpectrumLineData[k][r]
			allCurrentProteinInfo[iProteinID_Output_ColumnId] = "{"+sProteinID+"}"
			allCurrentProteinInfo[iProteinDescription_Output_ColumnId] = "{"+sProteinDescription+"}"
			iUniqueSpectrumCounts = 0
			for eachSpectrumLine in  allSpectrumLineData[i][j] :
				eachSpectrumInfo = eachSpectrumLine.split("\t")
				sProteinNameOfSpectrum = eachSpectrumInfo[iProteinName_ColumId]
				if EqualProteinName(allCurrentProteinInfo[iProteinID_Output_ColumnId], sProteinNameOfSpectrum, ",")  :
					iUniqueSpectrumCounts = iUniqueSpectrumCounts + 1
			allCurrentProteinInfo[iUniqueSpectrumCounts_Output_ColumnId] = str(iUniqueSpectrumCounts)
			sUpdatedProteinLine = ""
			for k in range(len(allCurrentProteinInfo)-1) :
				sUpdatedProteinLine = sUpdatedProteinLine + allCurrentProteinInfo[k] + "\t"
			sUpdatedProteinLine = sUpdatedProteinLine + allCurrentProteinInfo[-1]
			allProteinLineData[i][j] = sUpdatedProteinLine
	return (allProteinLineData, allSpectrumLineData)
	
def HandleFiles(psmFileName) :
	clusterFileName = psmFileName[:-12] + ".pro.cluster.txt"
	pro2psmClusterFileName = psmFileName[:-12] + ".pro2psm.cluster.txt"
	inputFile       = open(psmFileName)
	clusterFile     = file(clusterFileName, "w")
	pro2psmClusterFile = file(pro2psmClusterFileName, "w")
	bPassProteinFormatLine  = False # if true, the format line of protein has been gone through
	bPassSpectrumFormatLine = False # if true, the format line of spectrum has been gone through
	sproteinFormatLine = "+\tProteinID\tAverageEnrichmentLevel\tStandardDeviation\tUniqueSpectrumCounts\tTotalSpectrumCounts\tEnrichmentLevels\tProteinDescription"
	iProteinID_Output_ColumnId = getColumnId(sproteinFormatLine, "ProteinID")
	iUniqueSpectrumCounts_Output_ColumnId = getColumnId(sproteinFormatLine, "UniqueSpectrumCounts")
	iProteinDescription_Output_ColumnId = getColumnId(sproteinFormatLine, "ProteinDescription")
	proteinData = []
	allProteinLineData = []
	allSpectrumLineData= []
	for currentLine in inputFile:
		currentLine = currentLine.strip()
		if currentLine.startswith("#") :
			clusterFile.write(currentLine+"\n")
			pro2psmClusterFile.write(currentLine+"\n")
			continue
		if (currentLine.startswith("*") and (not (bPassSpectrumFormatLine))) :
			pro2psmClusterFile.write(currentLine+"\n")
			iSearchName_ColumnId = getColumnId(currentLine, "SearchName")
			iProteinName_ColumId = getColumnId(currentLine, "ProteinNames")
			iFilename_ColumnId   = getColumnId(currentLine, "Filename")
			iScanNumber_ColumnId = getColumnId(currentLine, "ScanNumber")
			bPassSpectrumFormatLine = True
			continue
		if (currentLine.startswith("+") and (not (bPassProteinFormatLine))) :
			clusterFile.write(sproteinFormatLine+"\n")
			pro2psmClusterFile.write(sproteinFormatLine+"\n")
			iProteinId_ColumnId = getColumnId(currentLine, "ProteinID")
			iProteinDescription_ColumnId = getColumnId(currentLine, "ProteinDescription")
			bPassProteinFormatLine = True
			continue
		if (currentLine.startswith("+") and (bPassProteinFormatLine) ) :
			if (len(proteinData) > 0) :
				#clusterSpectrum(proteinData)
			        (allProteinLines, allSpectrumLines) = clusterSpectrum(proteinData, iProteinId_ColumnId, iProteinDescription_ColumnId, iSearchName_ColumnId, iProteinName_ColumId)
			        if len(allProteinLines) > 0:
					allProteinLineData.append(allProteinLines)
					allSpectrumLineData.append(allSpectrumLines)
				#for i in range(len(allProteinLines)) :
               				#clusterFile.write(allProteinLines[i]+"\n")
               				#pro2psmClusterFile.write(allProteinLines[i]+"\n")
                			#for eachSpectrumLine in allSpectrumLines[i] :
                        			#pro2psmClusterFile.write(eachSpectrumLine+"\n")
			proteinData = []
			proteinData.append(currentLine)
		if (currentLine.startswith("*") and (bPassSpectrumFormatLine)) :
			proteinData.append(currentLine)
	(allProteinLines, allSpectrumLines) = clusterSpectrum(proteinData, iProteinId_ColumnId, iProteinDescription_ColumnId, iSearchName_ColumnId, iProteinName_ColumId)
	if len(allProteinLines) > 0:
		allProteinLineData.append(allProteinLines)
		allSpectrumLineData.append(allSpectrumLines)
	#for i in range(len(allProteinLines)) :
	#	clusterFile.write(allProteinLines[i]+"\n")
	#	pro2psmClusterFile.write(allProteinLines[i]+"\n")
	#	for eachSpectrumLine in allSpectrumLines[i] :
	#		pro2psmClusterFile.write(eachSpectrumLine+"\n")
	(allUpdatedProteinLineData, allUpdatedSpectrumLineData) = MergeData(allProteinLineData, allSpectrumLineData, iProteinID_Output_ColumnId, iUniqueSpectrumCounts_Output_ColumnId, iProteinDescription_Output_ColumnId, iProteinName_ColumId, iFilename_ColumnId, iScanNumber_ColumnId)
	
	for i in range(len(allUpdatedProteinLineData)) :
		for j in range(len(allUpdatedProteinLineData[i])):
			clusterFile.write(allUpdatedProteinLineData[i][j]+"\n")
			pro2psmClusterFile.write(allUpdatedProteinLineData[i][j]+"\n")
			for eachUpdatedSpectrum in allUpdatedSpectrumLineData[i][j] :
				pro2psmClusterFile.write(eachUpdatedSpectrum+"\n")
	inputFile.close()
	clusterFile.close()
	pro2psmClusterFile.close()
	
	
## +------+
## | Main |
## +------+
def main(argv=None):

    # try to get arguments and error handling
        if argv is None:
		argv = sys.argv
       		 # parse options
		(config_filename, psm_filename) = parse_options(argv)
	wholeDict = parseconfig.parseConfigKeyValues(config_filename)
	sClusteringThreshold = wholeDict.get("[Stable_Isotope_Probing]Clustering_Threshold")
	sMinPSMPerIsotopicCluster =  wholeDict.get("[Stable_Isotope_Probing]Min_PSM_Per_Isotopic_Cluster")
	global dClusteringThreshold
	global iMinPSMPerIsotopicCluster
	dClusteringThreshold = float(sClusteringThreshold[:-1])
	iMinPSMPerIsotopicCluster = int(sMinPSMPerIsotopicCluster)
	HandleFiles(psm_filename)
		


## If this program runs as standalone, then exit.
if __name__ == "__main__":
    sys.exit(main())

