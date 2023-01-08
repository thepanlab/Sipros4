
#!/usr/bin/python

import sys, getopt, warnings, os, re
import parseconfig


def parse_options(argv):

    try:
        opts, args = getopt.getopt(argv[1:], "hvw:c:",
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
            print "-c configure file -w working directory"
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
    
    

def GenerateOneSipConfig(lsElementPercent, dCurrentEnrichmentLevel, ConfigList, working_dir, dElementSum, iSipElementIndex, sSIPElement, iElementPercentLineId, iSearchNameLineId, iSIPElementIsotope, iSearchTypeLineId) :

#	print dCurrentEnrichmentLevel

	sSearchName = sSIPElement+str(iSIPElementIsotope)+"_"+str(int(round(dCurrentEnrichmentLevel*100)))+"Pct"
	dFirstPercent = dElementSum-dCurrentEnrichmentLevel
	if (dFirstPercent > 1) :
		dFirstPercent = 1
	elif (dFirstPercent < 0) :
		dFirstPercent = 0
 #       print dCurrentEnrichmentLevel
	sElementPercent = "Element_Percent{"+sSIPElement+"}\t=\t" + str( dFirstPercent  )
	
	for i in range(1, len(lsElementPercent)) :
		if (i == iSipElementIndex) :
			sElementPercent = sElementPercent + ",\t" + str(dCurrentEnrichmentLevel)
		else :
			sElementPercent = sElementPercent + ",\t" + lsElementPercent[i]	
	if not(working_dir.endswith("/")) :
		working_dir = working_dir + "/"
#	print sSearchName
	currentConfigFile = file(working_dir+"SiproConfig."+sSearchName+".cfg", "w")
	for i in range (len(ConfigList)) :
		if (i == iSearchNameLineId) :
			currentConfigFile.write("Search_Name\t=\t"+sSearchName+"\n")
		elif (i == iElementPercentLineId) :
			currentConfigFile.write(sElementPercent+"\n")
		elif (i == iSearchTypeLineId) :
			#currentConfigFile.write("Search_Type = Regular\n")
			currentConfigFile.write(ConfigList[i]) # keep the Search_Type = SIP
		else :
			currentConfigFile.write(ConfigList[i])
	currentConfigFile.close()
  
    
    
def GenerateSipConfig(config_filename, working_dir, lsElementPercent, dMaxEnrichmentLevel, dMinEnrichmentLevel, dEnrichmentLevelIncrement, iSipElementIndex, sSIPElement, iSIPElementIsotope) :
	originalConfig= open (config_filename)
	ConfigList = originalConfig.readlines()
	iSearchNameLineId = -1
	iElementPercentLineId = -1
	iSearchTypeLineId = -1
	for i in range (len(ConfigList)):
		currentLine = ConfigList[i].lstrip()
		if currentLine.startswith("#") :
			continue
		if currentLine.startswith("Element_Percent{"+sSIPElement+"}") :
			iElementPercentLineId = i
		if currentLine.startswith("Search_Name") : 
			iSearchNameLineId = i
		if currentLine.startswith("Search_Type") :
			iSearchTypeLineId = i
	dCurrentEnrichmentLevel = dMinEnrichmentLevel
	dElementSum = float(lsElementPercent[0])+float(lsElementPercent[iSipElementIndex])
	while ((dCurrentEnrichmentLevel < (dMaxEnrichmentLevel+0.0000001)) and (dCurrentEnrichmentLevel <= (dElementSum+0.0000001))) :
#		print dCurrentEnrichmentLevel
		GenerateOneSipConfig(lsElementPercent, dCurrentEnrichmentLevel, ConfigList, working_dir, dElementSum, 
			iSipElementIndex, sSIPElement, iElementPercentLineId, iSearchNameLineId, iSIPElementIsotope, iSearchTypeLineId)
		
		
		dCurrentEnrichmentLevel += dEnrichmentLevelIncrement

	originalConfig.close()
    
    
## +------+
## | Main |
## +------+
def main(argv=None):
# try to get arguments and error handling
	if argv is None:
		argv = sys.argv
		(working_dir, config_filename) = parse_options(argv)
	wholeDict = parseconfig.parseConfigKeyValues( config_filename)	
	sMaxEnrichmentLevel = wholeDict.get("[Stable_Isotope_Probing]Maximum_Enrichment_Level")
	sMinEnrichmentLevel = wholeDict.get("[Stable_Isotope_Probing]Minimum_Enrichment_Level")
	sEnrichmentLevelIncrement = wholeDict.get("[Stable_Isotope_Probing]Enrichment_Level_Increment")

	if ((sMaxEnrichmentLevel == None) or (sMinEnrichmentLevel == None) or (sEnrichmentLevelIncrement == None)) :
		print "Enrichment level information is incomplete."
		sys.exit(0)
	dMaxEnrichmentLevel  =  float(sMaxEnrichmentLevel[0:-1])/100
	dMinEnrichmentLevel  =  float(sMinEnrichmentLevel[0:-1])/100
	dEnrichmentLevelIncrement = float(sEnrichmentLevelIncrement[0:-1])/100
	#print dMaxEnrichmentLevel, dMinEnrichmentLevel, dEnrichmentLevelIncrement
	sSIPElement = wholeDict.get("[Stable_Isotope_Probing]SIP_Element")
	sSIPElementIsotope = wholeDict.get("[Stable_Isotope_Probing]SIP_Element_Isotope")
	if ((sSIPElement == None) or (sSIPElementIsotope == None)) :
		print "SIP information is incomplete."
		sys.exit(0)
	iSIPElementIsotope = int(sSIPElementIsotope)
	#print sSIPElement, iSIPElementIsotope
	sElementMasses = wholeDict.get("[Peptide_Identification]Element_Masses{"+sSIPElement+"}")
	sElementPercent = wholeDict.get("[Peptide_Identification]Element_Percent{"+sSIPElement+"}")
	if ((sElementMasses == None) or (sElementPercent == None)) :
		print "[Peptide_Identification]Element {"+sSIPElement+"}" + "is not available"
		sys.exit(0)
	sElementMasses   = sElementMasses.strip(",")
	sElementPercent  = sElementPercent.strip(",")
	lsElementMasses  = sElementMasses.split(",")
	lsElementPercent = sElementPercent.split(",")
	if (len(lsElementMasses) != len(lsElementPercent)) or (len(lsElementMasses) == 1):
		print "The number of elements is wrong"
		sys.exit(0)
	ldElementPercent = []
	iSipElementIndex = -1
	for i in range(len(lsElementPercent)) :
		ldElementPercent.append(float(lsElementPercent[i]))
		if (int(round(float(lsElementMasses[i])))  ==  iSIPElementIsotope) :
			iSipElementIndex = i
	if (iSipElementIndex == -1) :
		print "can't find the target element parcent."
		sys.exit(0)
	if (iSipElementIndex == 0) :#################################
		print "The first element can't be the target."
		sys.exit(0)
	GenerateSipConfig(config_filename, working_dir, lsElementPercent, dMaxEnrichmentLevel, 
		dMinEnrichmentLevel, dEnrichmentLevelIncrement, iSipElementIndex, sSIPElement, iSIPElementIsotope)





if __name__ == "__main__":
    sys.exit(main())




