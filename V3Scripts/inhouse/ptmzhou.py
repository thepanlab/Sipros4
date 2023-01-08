#!/usr/bin/python


## Import Python package modules
import sys, getopt, warnings, os, re
from datetime import datetime, date, time
from collections import namedtuple
from collections import defaultdict
import csv




def parse_options(argv):

    
    opts, args = getopt.getopt(argv[1:], "hw:c:o:",
                                    ["help",
                                     "working-dir",
				     "output-file",])


    # Default working dir and config file
    working_dir = "./"
    outputFileName = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-w workingdirectory -o outputfile"
            sys.exit(1)
        if option in ("-w", "--working-dir"):
            working_dir = value
            if working_dir[-1] != '/':
                working_dir = working_dir + '/'
        if option in ("-o", "--output-file"):
            outputFileName = value

    pro2ptm_filename_list = get_file_list_with_ext(working_dir, ".pro2ptm.txt")
    pro2ptmFileName = pro2ptm_filename_list[0]
    
    pro2psm_filename_list = get_file_list_with_ext(working_dir, ".pro2psm.txt")
    pro2psmFileName = pro2psm_filename_list[0]

    if (outputFileName == "") :
        (pro2ptmFileNameRoot, pro2ptmFileNameExt) = os.path.splitext(pro2ptmFileName) # ext is txt 
        (pro2ptmFileNameRoot, pro2ptmFileNameExt) = os.path.splitext(pro2ptmFileNameRoot) # ext is pro2pep
        outputFileName = pro2ptmFileNameRoot + ".pro2ptm_zhou.txt"
	#print outputFileName
    return [pro2psmFileName, pro2ptmFileName, outputFileName]

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
    outputFile.write("#\tDeltaZ = Highest difference between the best PSM and the next best PSM of this scan\n")
    outputFile.write("#\tDeltaZ/Score = Highest DeltaZ/Score\n")
    outputFile.write("#\n")
	#outputFile.write("+\tProteinID\tCoveredLength\tSequenceCoverage\tTotalPTMCounts\tTotalPeptideCounts\t")
	#outputFile.write("TotalSpectralCounts\tProteinDescription\n")
	#outputFile.write("*\tPosition\tResidue\tPTM\tNsequence\tCsequence\tModifiedPeptideCounts\tUnmodifiedPeptideCounts\t")
	#outputFile.write("UniqueModifiedPeptideCount\tUniqueUnmodifiedPeptideCount\tModifiedSpectralCount\t")
	#outputFile.write("UnmodifiedSpectralCount\tUniqueModifiedSpectralCount\tUniqueUnmodifiedSpectralCount\t")
	#outputFile.write("BestScore\tPeptides\tModifiedPSMs\tUnmodifiedPSMs\tDeltaZ\tDeltaZ/Score\n")	
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
	
def pro2psmFileInfo(pro2psm_file) :

    bPro2PsmFirstProteinLine = True
    bPro2PsmFirstSpectrumLine = True

    all_pro2psm = []
    current_pro2psm = []
    sCurrentProteinID = ""
    sNextProteinID = ""
    for sCurrentPro2PsmLine in pro2psm_file :
        sCurrentPro2PsmLine = sCurrentPro2PsmLine.strip()
        if ((sCurrentPro2PsmLine.startswith("#")) or (sCurrentPro2PsmLine == "")) :
            continue
        if (sCurrentPro2PsmLine.startswith("+")):
            if (bPro2PsmFirstProteinLine) : 
                bPro2PsmFirstProteinLine  = False
                Pro2Psm_ProteinID_ColumnID = getColumnId(sCurrentPro2PsmLine, "ProteinID")
            else:
                all_protein_info_list = sCurrentPro2PsmLine.split("\t")
                if (sCurrentProteinID == "") :
                    sCurrentProteinID = all_protein_info_list[Pro2Psm_ProteinID_ColumnID]
                    current_pro2psm = []
                    current_pro2psm.append(sCurrentProteinID)
                else:
                    sNextProteinID = all_protein_info_list[Pro2Psm_ProteinID_ColumnID]
                    all_pro2psm.append(current_pro2psm)
                    current_pro2psm = []
                    sCurrentProtienID = ""
                    current_pro2psm.append(sNextProteinID)

        if (sCurrentPro2PsmLine.startswith("*")):
            if (bPro2PsmFirstSpectrumLine) :
                bPro2PsmFirstSpectrumLine = False
                Pro2Psm_Filename_ColumnID = getColumnId(sCurrentPro2PsmLine, "Filename")
                Pro2Psm_ScanNumber_ColumnID = getColumnId(sCurrentPro2PsmLine, "ScanNumber")
                Pro2Psm_Score_ColumnID = getColumnId(sCurrentPro2PsmLine, "Score")
                Pro2Psm_DeltaZ_ColumnID = getColumnId(sCurrentPro2PsmLine, "DeltaZ")
            else :
                all_peptide_info_list = sCurrentPro2PsmLine.split("\t")
                if (sNextProteinID != "") :
                    sCurentProteinID = sNextProteinID
                    #if (sCurrentProteinID == "{CGL2_08970G0002,UBAL2_7931G0123}") :
                    #    print len(all_pro2psm)
                    sNextProteinID = ""
                    #current_pro2psm = []
                    #current_pro2psm.append(sCurentProteinID)
                psm_id = all_peptide_info_list[Pro2Psm_Filename_ColumnID]+"["+all_peptide_info_list[Pro2Psm_ScanNumber_ColumnID]+"]"
                dDeltaZ = float(all_peptide_info_list[Pro2Psm_DeltaZ_ColumnID])
                dDeltaZScore = dDeltaZ/float(all_peptide_info_list[Pro2Psm_Score_ColumnID])
                current_pro2psm.append((psm_id, [dDeltaZ, dDeltaZScore]))
                
    all_pro2psm.append(current_pro2psm)
    #print all_pro2psm[4039]
    #sys.exit(0)
    return all_pro2psm

def addColumns(pro2psmFileName, pro2ptmFileName, outputFileName) :
    pro2psm_file = open(pro2psmFileName)
    pro2ptm_file = open(pro2ptmFileName)
    output_file  = open(outputFileName, "a")

    all_pro2psm  = pro2psmFileInfo(pro2psm_file)
    
    bTitleProtein = True
    bTitlePtm     = True
    pro2psm_list_id = 0
    bCompressMode = False # pro2psm {p1,p2,p3}, pro2ptm p1, p2, p3
    for sCurrentPro2ptmLine in pro2ptm_file :
        sCurrentPro2ptmLine = sCurrentPro2ptmLine.strip()
        if (sCurrentPro2ptmLine.startswith("#") or (sCurrentPro2ptmLine == "")) :
            continue
        if sCurrentPro2ptmLine.startswith("+") :
            output_file.write(sCurrentPro2ptmLine+"\n")
            if bTitleProtein :
                bTitleProtein = False
                Pro2Ptm_ProteinID_ColumnID = getColumnId(sCurrentPro2ptmLine, "ProteinID")
            else:
                line_info_list = sCurrentPro2ptmLine.split("\t")
                sCurrentProteinID = line_info_list[Pro2Ptm_ProteinID_ColumnID]
    #            print pro2psm_list_id
                currentPsmInfo_list = all_pro2psm[pro2psm_list_id]
                #print currentPsmInfo_list 
                if (currentPsmInfo_list[0] != sCurrentProteinID) :
                    if  (currentPsmInfo_list[0].find(sCurrentProteinID) >= 0 ) :
                        print sCurrentProteinID, "in", currentPsmInfo_list[0]
                        if not(bCompressMode) :
                            bCompressMode = True
                            iNumProteinTogether = len(currentPsmInfo_list[0].split(",")) 
                            iNumProteinTogether -= 1
                        else :
                            iNumProteinTogether -= 1
                            if (iNumProteinTogether == 0) :
                                bCompressMode = False
                                pro2psm_list_id += 1
     #                   print iNumProteinTogether, pro2psm_list_id                        

                    else :
                        print "protein id don't match:", currentPsmInfo_list[0], sCurrentProteinID
                        sys.exit(1)
                else :
                    pro2psm_list_id += 1
                #print sCurrentProteinID
                currentTuple = currentPsmInfo_list[1:]
        if sCurrentPro2ptmLine.startswith("*") :
            output_file.write(sCurrentPro2ptmLine)
            if bTitlePtm :
                bTitlePtm = False
                output_file.write("\tDeltaZ\tDeltaZ/Score\n")
                Pro2Ptm_ModifiedPSMs_ColumnID = getColumnId(sCurrentPro2ptmLine, "ModifiedPSMs")
            else :
                line_info_list = sCurrentPro2ptmLine.split("\t")
                sModifiedPSMs  = line_info_list[Pro2Ptm_ModifiedPSMs_ColumnID]
                sModifiedPSMs  = sModifiedPSMs.strip("{}")
                dDeltaZ = -1
                dDeltaZScore = -1
                if (sModifiedPSMs != "") :
                    psm_list = sModifiedPSMs.split(",")
                    for each_psm in psm_list :
                        psm_info_list = dict(currentTuple).get(each_psm)
                        if (psm_info_list != None) :
                            current_DeltaZ= psm_info_list[0]
                            current_DeltaZScore = psm_info_list[1]
                            if (dDeltaZ < current_DeltaZ) :
                                dDeltaZ = current_DeltaZ
                            if (dDeltaZScore < current_DeltaZScore) :
                                dDeltaZScore = current_DeltaZScore
                        else :
                            print "Can't find", each_psm, "in", sCurrentProteinID, "!"
                if (dDeltaZ < -0.5) :
                    output_file.write("\tNA")
                else :
                    output_file.write("\t"+str(dDeltaZ)) 
                if (dDeltaZScore < -0.5) :
                    output_file.write("\tNA\n")
                else:
                    output_file.write("\t"+str(dDeltaZScore)+"\n")


    pro2psm_file.close()
    pro2ptm_file.close()
    output_file.close()
	
def main(argv=None):

    # try to get arguments and error handling
        if argv is None:
		argv = sys.argv
       		 # parse options
		[pro2psmFileName, pro2ptmFileName, outputFileName] = parse_options(argv)
	outputComments(outputFileName)
	addColumns(pro2psmFileName, pro2ptmFileName, outputFileName)
	



## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()
