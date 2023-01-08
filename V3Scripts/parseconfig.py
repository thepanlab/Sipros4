#!/usr/bin/python

import sys, getopt, warnings, os, re


def getConfigMasterKeyValue (sMasterKey, dictConfigKeyValues) :
	dictKeyValueSet = {}	
	for currentKey, currentValue in dictConfigKeyValues.items() :
		if( currentKey.startswith(sMasterKey+"{") and currentKey.endswith("}") and ((len(sMasterKey) + 2) < len(currentKey))) :
			coreKey = currentKey[len(sMasterKey)+1: len(currentKey)-1]
			dictKeyValueSet [coreKey] = currentValue
	return dictKeyValueSet



def parseConfigLine (sLine, sSectionName) :
# sSectionName is a list, but only the first value is used
	currentKey = ""
	currentValue = ""
	if ((sLine[0] == "[") and (sLine[len(sLine)-1] == "]" )) :
		sSectionName[0] = sLine
	elif (sSectionName[0] == "")   : 
		print "no section name"
	else :
		twoParts = sLine.split("=")
		if (len(twoParts) != 2) :
			print "wrong line: "+sLine
		else :
			currentKey   = sSectionName[0] + twoParts[0].rstrip()
			currentValue = twoParts[1].lstrip()
			#print currentKey+"=>"+currentValue
	return [currentKey, currentValue]


def parseConfigKeyValues (filepath) :
	configFile = open(filepath)
	
	sSectionName = [""]
	dictConfigKeyValues = {}
	for sLine in configFile.readlines() :
		poundPos = sLine.find("#")
		if (poundPos > -1) :
			sLine = sLine[0:poundPos]
		sLine = sLine.strip()
		if (sLine == "") :
			continue
		else :
			#print "!!!"+sLine+"!!!"
			currentKey, currentValue = parseConfigLine (sLine, sSectionName)
			if (currentKey != "") and (currentValue != "" ):
			#	print currentKey+"=>"+currentValue
				if (dictConfigKeyValues.get(currentKey) == None ):
					dictConfigKeyValues[currentKey] = currentValue
				else :
					print currentKey + " has existed"	
	configFile.close()
#	for currentKey, currentValue in dictConfigKeyValues.items() :
#		print currentKey, currentValue 
	return dictConfigKeyValues
