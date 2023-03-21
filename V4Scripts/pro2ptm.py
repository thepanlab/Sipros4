#!/usr/bin/python


## Import Python package modules
import sys, getopt, warnings, os, re
from datetime import datetime, date, time
from collections import namedtuple
from collections import defaultdict
import csv
from Bio import SeqIO
from Bio.Seq import Seq


import parseconfig


def parse_options(argv):

    
    opts, args = getopt.getopt(argv[1:], "hw:c:o:",
                                    ["help",
                                     "working-dir",
                                     "config-file",
				     "output-file",])


    # Default working dir and config file
    working_dir = "./"
    config_file = "SiprosConfig.cfg"
    outputFileName = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-c configurefile -w workingdirectory -o outputfile"
            sys.exit(1)
        if option in ("-w", "--working-dir"):
            working_dir = value
            if working_dir[-1] != '/':
                working_dir = working_dir + '/'
        if option in ("-c", "--config-file"):
            config_file = value
        if option in ("-o", "--output-file"):
            outputFileName = value

    pro2pep_filename_list = get_file_list_with_ext(working_dir, ".pro2pep.txt")
    pro2pepFileName = pro2pep_filename_list[0]
    
    wholeDict = parseconfig.parseConfigKeyValues(config_file)
    databaseFileName = wholeDict.get("[Peptide_Identification]FASTA_Database")

    if (outputFileName == "") :
	(pro2pepFileNameRoot, pro2pepFileNameExt) = os.path.splitext(pro2pepFileName) # ext is txt 
	expandFileName = pro2pepFileNameRoot + ".expand.txt"
	(pro2pepFileNameRoot, pro2pepFileNameExt) = os.path.splitext(pro2pepFileNameRoot) # ext is pro2pep
	outputFileName = pro2pepFileNameRoot + ".pro2ptm.txt"
	#print outputFileName
    return [databaseFileName, pro2pepFileName, outputFileName, expandFileName]

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

    
def outputComments(outputFileName) :
	outputFile = open(outputFileName, "w")
	outputFile.write("#\t[Column_Names]\n")
	outputFile.write("#\n")
	outputFile.write("#\t+ = Marker of a protein line\n")
	outputFile.write("#\tProteinID = Names of the protein\n")	
	outputFile.write("#\tCoveredLength = Number of residues in covered protein sequence\n")
	outputFile.write("#\tSequenceCoverage = Percentage of protein sequence covered\n")	
	outputFile.write("#\tTotalPTMCounts = Number of all PTMs\n")
	outputFile.write("#\tTotalPeptideCounts = Number of all peptides\n")
	outputFile.write("#\tTotalSpectralCounts = Number of all PSMs\n")
	outputFile.write("#\tProteinDescription = Protein description\n")
	outputFile.write("#\n")
	outputFile.write("#\t* = Marker of a PTM line\n")
	outputFile.write("#\tPosition = Position of the modified residue in the protein\n")
	outputFile.write("#\tResidue = Amino acid type of the modified residue\n")	
	outputFile.write("#\tPTM = PTM type of the modified residue\n")	
	outputFile.write("#\tNsequence = N-terminal flanking sequence\n")
	outputFile.write("#\tCsequence = C-terminal flanking sequence\n")	
	outputFile.write("#\tModifiedPeptideCounts = Number of modified peptides\n")
	outputFile.write("#\tUnmodifiedPeptideCounts = Number of unmodified peptides\n")	
	outputFile.write("#\tUniqueModifiedPeptideCount = Number of unique modified peptides\n")	
	outputFile.write("#\tUniqueUnmodifiedPeptideCount = Number of unique unmodified peptides\n")	
	outputFile.write("#\tModifiedSpectralCount = Number of modified spectra\n")
	outputFile.write("#\tUnmodifiedSpectralCount = Number of unmodified spectra\n")
	outputFile.write("#\tUniqueModifiedSpectralCount = Number of unique modified spectra\n")	
	outputFile.write("#\tUniqueUnmodifiedSpectralCount = Number of unique unmodified spectra\n")
	outputFile.write("#\tBestScore = The best score of all modified peptides\n")
	outputFile.write("#\tPeptides = List of modified peptides\n")
	outputFile.write("#\tModifiedPSMs = List of modified PSMs: FT2_Filename[Scan_Number]\n")
	outputFile.write("#\tUnmodifiedPSMs = List of unmodified PSMs: FT2_Filename[Scan_Number]\n")
	outputFile.write("#\n")
	outputFile.write("+\tProteinID\tCoveredLength\tSequenceCoverage\tTotalPTMCounts\tTotalPeptideCounts\t")
	outputFile.write("TotalSpectralCounts\tProteinDescription\n")
	outputFile.write("*\tPosition\tResidue\tPTM\tNsequence\tCsequence\tModifiedPeptideCounts\tUnmodifiedPeptideCounts\t")
	outputFile.write("UniqueModifiedPeptideCount\tUniqueUnmodifiedPeptideCount\tModifiedSpectralCount\t")
	outputFile.write("UnmodifiedSpectralCount\tUniqueModifiedSpectralCount\tUniqueUnmodifiedSpectralCount\t")
	outputFile.write("BestScore\tPeptides\tModifiedPSMs\tUnmodifiedPSMs\n")	
	outputFile.close()

def getColumnId(sColumnNameLine, ColumnName) :
	ColumnName_list = sColumnNameLine.split("\t")
	try:
		iColumnId = ColumnName_list.index(ColumnName)
	except ValueError:
		print "can't find column "+ColumnName
		sys.exit(0)
	#print iColumnId
	return iColumnId	
	
	
def findSeq(seqId, databaseSeqList):
	for eachRecord in databaseSeqList :
	      if eachRecord.id == seqId :
		    foundRecord = eachRecord
		    break
#	print seqId
	return foundRecord

def IdentifyOriginalPeptide(sProteinSeq, sOriginalPeptide) :
  # identify peptides on the protein
	iNtermPosition = sOriginalPeptide.find("[")
	iCtermPosition = sOriginalPeptide.find("]")
	if ((iNtermPosition == -1) or (iCtermPosition == -1)) :
	      print "illegal original peptide :"+sOriginalPeptide 
	      sys.exit(0)
	sOriginalPeptideRoot       = sOriginalPeptide[iNtermPosition+1 : iCtermPosition]
	iOriginalPeptideRootLength = len(sOriginalPeptideRoot)
	iStartPosition = 0
	iPeptidePositionList = []
	iPeptidePosition = sProteinSeq.find(sOriginalPeptideRoot, iStartPosition)
	while (iPeptidePosition != -1) :
	      iPeptidePositionList.append([iPeptidePosition, iPeptidePosition + iOriginalPeptideRootLength - 1])
	      iStartPosition = iPeptidePosition + 1
	      iPeptidePosition = sProteinSeq.find(sOriginalPeptideRoot, iStartPosition )
	return iPeptidePositionList

def IdentifyPTMPosition(sIdentifiedPeptide_ori, sOriginalPeptide_ori) :
    # identify ptms 
	PTM2PepPositionList = []	
	iNtermPosition = sOriginalPeptide_ori.find("[")
	iCtermPosition = sOriginalPeptide_ori.find("]")
	sOriginalPeptide   = sOriginalPeptide_ori[iNtermPosition : iCtermPosition+1]
	if ((sIdentifiedPeptide_ori[-1] != "]") and (sIdentifiedPeptide_ori[-1].isalpha())) :
	      sIdentifiedPeptide = sIdentifiedPeptide_ori[iNtermPosition : -1]
	else :
	      sIdentifiedPeptide = sIdentifiedPeptide_ori[iNtermPosition : ]
	if (sIdentifiedPeptide != sOriginalPeptide) :
	      sCurrentIdentifiedPeptide = sIdentifiedPeptide
	      for i in range(1, len(sOriginalPeptide)) :
		    if (sCurrentIdentifiedPeptide[i] != sOriginalPeptide[i]) :  
			  #get a ptm
			  PTM2PepPositionList.append([i-1, sCurrentIdentifiedPeptide[i]])
			  #[position on peptide, ptmtype]
			  sCurrentIdentifiedPeptide = sCurrentIdentifiedPeptide[:i] + sCurrentIdentifiedPeptide[i+1:]
	      if (sIdentifiedPeptide[-1] != "]") :
		  #end character is a ptmtype
		    PTM2PepPositionList.append([len(sOriginalPeptide) -1 , sIdentifiedPeptide[-1]])
	#print PTM2PepPositionList
	return PTM2PepPositionList
	

def IdentifyPTM2Pro(PTM2PepPositionList, peptidePositionsList, iProteinLength) :
	PTM2ProteinList = []
	NonPTM2ProteinList = []
	PTM2ProteinArray = [0] * (iProteinLength +2) # 1: unmodified 2: modified 0: uncovered
	for proteinPosition in peptidePositionsList :
	      proteinPositionBegin = proteinPosition[0]
	      proteinPositionEnd   = proteinPosition[1]
	      peptideLength        = proteinPositionEnd - proteinPositionBegin + 1 + 2
	      if (proteinPositionBegin == 0) and (PTM2ProteinArray[0] == 0) :
		    PTM2ProteinArray[0] = 1
	      if (proteinPositionEnd == (iProteinLength-1)) and (PTM2ProteinArray[iProteinLength+1] == 0) :
		    PTM2ProteinArray[iProteinLength+1] = 1
	      for i in range(proteinPositionBegin, proteinPositionEnd+1) :
		    if (PTM2ProteinArray[i+1] == 0) :
			  PTM2ProteinArray[i+1] = 1
	      for PTM2PeptidePositionList  in PTM2PepPositionList :
		    PTM2PeptidePosition = PTM2PeptidePositionList[0]
		    PTM2PeptideType     = PTM2PeptidePositionList[1]
		    PTM2ProteinPosition = proteinPositionBegin + PTM2PeptidePositionList[0] - 1
		    if (PTM2PeptidePosition == 0) and (PTM2ProteinPosition != -1) :
			  continue
		    if (PTM2PeptidePosition == peptideLength -1) and (PTM2ProteinPosition != iProteinLength) :
			  continue
		    # -1 and iProteinLength are [ and ] for protein, respectively
		    PTM2ProteinList.append([PTM2ProteinPosition, PTM2PeptideType])
		    PTM2ProteinArray[PTM2ProteinPosition + 1] = 2
	for i in range(iProteinLength +2) :
	      if (PTM2ProteinArray[i] == 1) :
		    NonPTM2ProteinList.append([i-1, ""])
#	print PTM2ProteinList
#	print NonPTM2ProteinList
	return [PTM2ProteinList,  NonPTM2ProteinList]

def OriginizePTMInfo(sPeptideCoverList, iTotalSpectrumCounts, iTotalPeptideCounts, 
		     sProteinSeq, iProteinSeqCoverList, sProteinId, sProteinDescription) :
	sWholeProteinSeq = "[" + sProteinSeq + "]" 
	sProteinLine     = "+\t" + sProteinId + "\t"
	sPTMLineList     = []
	iRealProteinLength  = len(sProteinSeq)
	iCoveredLength   = sum(iProteinSeqCoverList)
	dSequenceCoverage= iCoveredLength/float(iRealProteinLength)
	#iTotalPTMCounts = 0
	#CoveredLength
	sProteinLine = sProteinLine + str(iCoveredLength) + "\t"
	#SequenceCoverage
	sProteinLine = sProteinLine + str(dSequenceCoverage) + "\t"
#	print sPeptideCoverList[0]
	for i in range (len(sPeptideCoverList)) :
	      if (len(sPeptideCoverList[i]) > 0) :
		    iProteinPosition = i
		    sResidue         = sWholeProteinSeq[i]
		    if (i == 0) :
			  sNsequence = "["
		    elif (i <= 11) :
			  sNsequence = sWholeProteinSeq[:i]
		    else :
			  sNsequence = "[" + sWholeProteinSeq[i-10:i]
		    if (i == len(sWholeProteinSeq) -1) :
			  sCsequence = "]"
		    elif (i >= len(sWholeProteinSeq) -12) :
			  sCsequence = sWholeProteinSeq[i+1:]
		    else :
			  sCsequence = sWholeProteinSeq[i+1: i+11] + "]"
		    iUnmodifiedPeptideCounts = 0
		    iUniqueUnmodifiedPeptideCount = 0
		    iUnmodifiedSpectralCount = 0
		    iUniqueUnmodifiedSpectralCount = 0
		    sUnmodifiedPSMs = ""
		    PTMTypeList      = []
		    PTMOrganizedInfo = []
		    for j in range(len(sPeptideCoverList[i])) :
			  sCurrentPTMType    = sPeptideCoverList[i][j][0]
			  iSpectralCount     = sPeptideCoverList[i][j][1]
			  sProteinCount      = sPeptideCoverList[i][j][2]
			  dBestScore         = sPeptideCoverList[i][j][3]
			  sPSMs              = sPeptideCoverList[i][j][4]
			  sIdentifiedPeptide = sPeptideCoverList[i][j][5]
			  if (sCurrentPTMType != "") :
				if (sProteinCount == "1") :
				      iCurrentUniqueModifiedPeptideCounts = 1
				      iCurrentUniqueModifiedSpectralCount = iSpectralCount
				else :
				      iCurrentUniqueModifiedPeptideCounts = 0
				      iCurrentUniqueModifiedSpectralCount = 0
				if sCurrentPTMType in PTMTypeList :
				      iPTMTypeId = PTMTypeList.index(sCurrentPTMType)
				      sPSMs = ","+ sPSMs[1 : -1] #ModifiedPSMs
				      sIdentifiedPeptide = ","+sIdentifiedPeptide
				      #ModifiedPeptideCounts
				      PTMOrganizedInfo[iPTMTypeId][0] = PTMOrganizedInfo[iPTMTypeId][0] + 1
				      #UniqueModifiedPeptideCount
				      PTMOrganizedInfo[iPTMTypeId][1] = PTMOrganizedInfo[iPTMTypeId][1] + iCurrentUniqueModifiedPeptideCounts
				      #ModifiedSpectralCount
				      PTMOrganizedInfo[iPTMTypeId][2] = PTMOrganizedInfo[iPTMTypeId][2] + iSpectralCount
				      # UniqueModifiedSpectralCount 
				      PTMOrganizedInfo[iPTMTypeId][3] = PTMOrganizedInfo[iPTMTypeId][3] + iCurrentUniqueModifiedSpectralCount
				      # BestScore
				      if (PTMOrganizedInfo[iPTMTypeId][4] < dBestScore) :
					      PTMOrganizedInfo[iPTMTypeId][4] = dBestScore
				      # Peptides
				      PTMOrganizedInfo[iPTMTypeId][5] = PTMOrganizedInfo[iPTMTypeId][5] + sIdentifiedPeptide
				      # ModifiedPSMs
				      PTMOrganizedInfo[iPTMTypeId][6] = PTMOrganizedInfo[iPTMTypeId][6] + sPSMs
				      
				else :
				      iPTMTypeId = len(PTMTypeList)
				      PTMTypeList.append(sCurrentPTMType)
				      sPSMs = sPSMs[1 : -1]
				      # ModifiedPeptideCounts, UniqueModifiedPeptideCount, ModifiedSpectralCount,
				      # UniqueModifiedSpectralCount, BestScore, Peptides, PSMs
				      PTMOrganizedInfo.append([1, iCurrentUniqueModifiedPeptideCounts, iSpectralCount, 
							      iCurrentUniqueModifiedSpectralCount, dBestScore,
							      sIdentifiedPeptide, sPSMs])
			  else :
				iUnmodifiedPeptideCounts = iUnmodifiedPeptideCounts + 1
				iUnmodifiedSpectralCount = iUnmodifiedSpectralCount + iSpectralCount
				if (sProteinCount == "1") :
				      iUniqueUnmodifiedPeptideCount = iUniqueUnmodifiedPeptideCount + 1
				      iUniqueUnmodifiedSpectralCount= iUniqueUnmodifiedSpectralCount+ iSpectralCount
				if (sUnmodifiedPSMs != "") :
				      sUnmodifiedPSMs = sUnmodifiedPSMs + ","
				sUnmodifiedPSMs = sUnmodifiedPSMs + sPSMs[1 : -1]
		    for eachPTMId in range(len(PTMOrganizedInfo)) :
			  sPTMLine = "*\t"
			  #Position, Residue, PTM, Nsequence, Csequence
			  sPTMLine = sPTMLine+str(i)+"\t"+sWholeProteinSeq[i]+"\t"+PTMTypeList[eachPTMId]+"\t"+sNsequence+"\t"+sCsequence+"\t"
			  #ModifiedPeptideCounts, UnmodifiedPeptideCounts
			  sPTMLine = sPTMLine+str(PTMOrganizedInfo[eachPTMId][0])+"\t"+str(iUnmodifiedPeptideCounts)+"\t"
			  #UniqueModifiedPeptideCount, UniqueUnmodifiedPeptideCount
			  sPTMLine = sPTMLine+str(PTMOrganizedInfo[eachPTMId][1])+"\t"+str(iUniqueUnmodifiedPeptideCount)+"\t"
			  #ModifiedSpectralCount, UnmodifiedSpectralCount
			  sPTMLine = sPTMLine+str(PTMOrganizedInfo[eachPTMId][2])+"\t"+str(iUnmodifiedSpectralCount)+"\t"
			  #UniqueModifiedSpectralCount, UniqueUnmodifiedSpectralCount
			  sPTMLine = sPTMLine+str(PTMOrganizedInfo[eachPTMId][3])+"\t"+str(iUniqueUnmodifiedSpectralCount)+"\t"
			  #BestScore, Peptides
			  sPTMLine = sPTMLine+str(PTMOrganizedInfo[eachPTMId][4])+"\t{"+str(PTMOrganizedInfo[eachPTMId][5])+"}\t"
			  #ModifiedPSMs
			  sPTMLine = sPTMLine+"{"+str(PTMOrganizedInfo[eachPTMId][6])+"}\t"
			  #UnmodifiedPSMs
			  sPTMLine = sPTMLine+"{"+sUnmodifiedPSMs+"}"
			  
			  sPTMLineList.append(sPTMLine)
			  
			  
	# TotalPTMCounts
	sProteinLine = sProteinLine + str(len(sPTMLineList)) + "\t"
	# TotalPeptideCounts
	sProteinLine = sProteinLine + str(iTotalPeptideCounts) + "\t"
	# TotalSpectralCounts
	sProteinLine = sProteinLine + str(iTotalSpectrumCounts) + "\t"
	# ProteinDescription
	sProteinLine = sProteinLine + sProteinDescription
	
	return[sProteinLine, sPTMLineList]
	
	
def handlePro2PTM(databaseFileName, pro2pepFileName, outputFileName) :	
	bProTitleLine = True
	bPepTitleLine = True
	bFirstProNonTitleLine = True
	pro2pepFile = open(pro2pepFileName, "r")
	outputFile  = open(outputFileName, "a")
	databaseSeqList = list(SeqIO.parse( databaseFileName, "fasta" ))
	for eachLine in pro2pepFile :
	      eachLine = eachLine.strip()
	      if (eachLine == "") :
		    continue
	      if eachLine.startswith("#") :
		    continue
	      if eachLine.startswith("+") : # protein line
		    if (bProTitleLine) :
			  bProTitleLine = False
			  iProteinId = getColumnId(eachLine, "ProteinID")
			  iProteinDescriptionId = getColumnId(eachLine, "ProteinDescription")
		    else :
			  if (bFirstProNonTitleLine) :
				bFirstProNonTitleLine = False
			  else:
			        [sProteinLine, sPTMLineList] =OriginizePTMInfo(sPeptideCoverList, iTotalSpectrumCounts, iTotalPeptideCounts, 
						 sProteinSeq, iProteinSeqCover, sProteinId, sProteinDescription)
				outputFile.write(sProteinLine+"\n")
				for eachPTMLine in sPTMLineList :
				      outputFile.write(eachPTMLine+"\n")
			  sProteinInfoList = eachLine.split("\t")
			  sProteinId       = sProteinInfoList[iProteinId]
			  sProteinDescription = sProteinInfoList[iProteinDescriptionId]
			  proteinRecord    = findSeq(sProteinId, databaseSeqList)
			  sProteinSeq      = proteinRecord.seq.tostring()
			  iProteinLength   = len(sProteinSeq)
			  iProteinSeqCover = [0]* iProteinLength
			  iTotalSpectrumCounts = 0
			  #iTotalPTMCounts  = 0
			  iTotalPeptideCounts = 0
			  sPeptideCoverList = [[]  for x in range(iProteinLength + 2)] # 0 for Nterm, the last for Cterm
			  #print sPeptideCoverList
			  
			  
	      if eachLine.startswith("*") : # peptide line
		    if (bPepTitleLine) :
			  bPepTitleLine = False
			  iIdentifiedPeptideId = getColumnId(eachLine, "IdentifiedPeptide")
			  iOriginalPeptideId   = getColumnId(eachLine, "OriginalPeptide")
			  iSpectralCountId     = getColumnId(eachLine, "SpectralCount")
			  iBestScoreId         = getColumnId(eachLine, "BestScore")
			  iProteinCountId      = getColumnId(eachLine, "ProteinCount")
			  iPSMsId              = getColumnId(eachLine, "PSMs")
		    else :
			  sPeptideInfoList     = eachLine.split("\t")
			  sIdentifiedPeptide   = sPeptideInfoList[iIdentifiedPeptideId]
			  sOriginalPeptide     = sPeptideInfoList[iOriginalPeptideId]			  
			  iSpectralCount       = int(sPeptideInfoList[iSpectralCountId])
			  sProteinCount        = sPeptideInfoList[iProteinCountId] # "1" : unique
			  dBestScore           = float(sPeptideInfoList[iBestScoreId])
			  sPSMs                = sPeptideInfoList[iPSMsId]
			  peptidePositionsList = IdentifyOriginalPeptide(sProteinSeq, sOriginalPeptide)
			  if (len(peptidePositionsList) == 0) :
			  # no peptide identified due to illegal characters in protein sequence
				continue
			  iTotalSpectrumCounts = iTotalSpectrumCounts + iSpectralCount
			  iTotalPeptideCounts  = iTotalPeptideCounts + 1
			  PTM2PepPositionList  = IdentifyPTMPosition(sIdentifiedPeptide, sOriginalPeptide)
			  for i in range (len(peptidePositionsList)) :
				for j in range(peptidePositionsList[i][0], peptidePositionsList[i][1] + 1) :
				      iProteinSeqCover[j] = 1 
			  [PTM2ProMapList, NonPTM2ProMapList] = IdentifyPTM2Pro(PTM2PepPositionList, peptidePositionsList, iProteinLength)
			  #print sPeptideCoverList[2]
			  for PTM2ProMap in PTM2ProMapList :
				currentPTMInfo = []
				currentPTMInfo.append(PTM2ProMap[1])
				currentPTMInfo.append(iSpectralCount)
				currentPTMInfo.append(sProteinCount)
				currentPTMInfo.append(dBestScore)
				currentPTMInfo.append(sPSMs)
				currentPTMInfo.append(sIdentifiedPeptide)
				#print PTM2ProMap[0]+1
				sPeptideCoverList[PTM2ProMap[0]+1].append(currentPTMInfo)
				#print sPeptideCoverList[2]
			  #print "!"
			  for NonPTM2ProMap in NonPTM2ProMapList :
				currentNonPTMInfo = []
				currentNonPTMInfo.append(NonPTM2ProMap[1])
				currentNonPTMInfo.append(iSpectralCount)
				currentNonPTMInfo.append(sProteinCount)
				currentNonPTMInfo.append(dBestScore)
				currentNonPTMInfo.append(sPSMs)
				currentNonPTMInfo.append(sIdentifiedPeptide)
				#print NonPTM2ProMap[0]+1
				sPeptideCoverList[NonPTM2ProMap[0]+1].append(currentNonPTMInfo)
			  #print sPeptideCoverList
			  #sys.exit(0)
	[sProteinLine, sPTMLineList] = OriginizePTMInfo(sPeptideCoverList, iTotalSpectrumCounts, iTotalPeptideCounts, 
							sProteinSeq, iProteinSeqCover, sProteinId, sProteinDescription)
							
	outputFile.write(sProteinLine+"\n")
	for eachPTMLine in sPTMLineList :
	      outputFile.write(eachPTMLine+"\n")							
	pro2pepFile.close()
	outputFile.close()

def writeoutputFile(sProteinLine, sPeptideLineList, iProteinIdColumnId, iProteindescriptionColumnId, outputFile) :
	sProteinLineList = []
	sProteinInfoList = sProteinLine.split("\t")
	sProteinId       = sProteinInfoList[iProteinIdColumnId]
	sProteinDescription = sProteinInfoList[iProteindescriptionColumnId]
	if (sProteinId.startswith("{") and sProteinId.endswith("}")) :
		sProteinIdList = sProteinId[1:-1].split(",")
		sProteinDescriptionList = sProteinDescription[1:-1].split(",")
		for i in range(len(sProteinIdList)) :
			sCurrentProteinId = sProteinIdList[i].strip()
			sProteinDescription = sProteinDescriptionList[i].strip()
			sCurrentProteinLine = ""
			for j in range(len(sProteinInfoList)) :
				if (sCurrentProteinLine == "" ) :
					sHeadTab = ""
				else :
					sHeadTab = "\t"
				if (j == iProteinIdColumnId) :
					sCurrentProteinLine = sCurrentProteinLine + sHeadTab + sCurrentProteinId
				elif (j == iProteindescriptionColumnId) :
					sCurrentProteinLine = sCurrentProteinLine + sHeadTab + sProteinDescription
				else:
					sCurrentProteinLine = sCurrentProteinLine + sHeadTab + sProteinInfoList[j]
			sProteinLineList.append(sCurrentProteinLine)

	else :
		sProteinLineList = [sProteinLine]
	#print len(sProteinLineList)
	for eachProteinLine in sProteinLineList :
		outputFile.write(eachProteinLine+"\n")
		for eachPeptideLine in sPeptideLineList :
			outputFile.write(eachPeptideLine+"\n")

def ExpandPro2Pep(pro2pepFileName, outputFileName) :

	pro2pepFile = open(pro2pepFileName, "r")
	outputFile  = open(outputFileName,  "w")
	bProteinTitleLine = True
	bPeptideTitleLine = True
	bFirstRealProteinLine = True
	for eachLine in pro2pepFile :
		eachLine = eachLine.strip()
		if (eachLine == "") :
			continue
		if (eachLine.startswith("#")) :
			outputFile.write(eachLine+"\n")
		if (eachLine.startswith("+")) :
			if (bProteinTitleLine) :
				bProteinTitleLine = False
				outputFile.write(eachLine+"\n")
				iProteinIdColumnId = getColumnId(eachLine, "ProteinID")
				iProteindescriptionColumnId = getColumnId(eachLine, "ProteinDescription")
			else:
				if (bFirstRealProteinLine) :
					bFirstRealProteinLine = False
				else :
					writeoutputFile(sProteinLine, sPeptideLineList, iProteinIdColumnId, iProteindescriptionColumnId, outputFile)
				sProteinLine = eachLine
				sPeptideLineList = []
		if (eachLine.startswith("*")) :
			if (bPeptideTitleLine) :
				bPeptideTitleLine = False
				outputFile.write(eachLine+"\n")
			else :
				sPeptideLineList.append(eachLine) 
	writeoutputFile(sProteinLine, sPeptideLineList, iProteinIdColumnId, iProteindescriptionColumnId, outputFile)
	pro2pepFile.close()
	outputFile.close()
	
	
def main(argv=None):

    # try to get arguments and error handling
        if argv is None:
		argv = sys.argv
       		 # parse options
		[databaseFileName, pro2pepFileName, outputFileName, expandFileName] = parse_options(argv)
	ExpandPro2Pep(pro2pepFileName, expandFileName)
	outputComments(outputFileName)
	handlePro2PTM(databaseFileName, expandFileName, outputFileName)
	os.remove(expandFileName)



## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()
