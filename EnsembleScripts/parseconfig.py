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
		print("no section name")
	else :
		twoParts = sLine.split("=")
		if (len(twoParts) != 2) :
			print("wrong line: "+sLine)
		else :
			currentKey   = sSectionName[0] + twoParts[0].rstrip()
			currentValue = twoParts[1].lstrip()
			#print(currentKey+"=>"+currentValue)
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
			#print("!!!"+sLine+"!!!")
			currentKey, currentValue = parseConfigLine (sLine, sSectionName)
			if (currentKey != "") and (currentValue != "" ):
			#	print(currentKey+"=>"+currentValue)
				if (dictConfigKeyValues.get(currentKey) == None ):
					dictConfigKeyValues[currentKey] = currentValue
				else :
					print(currentKey + " has existed")
	configFile.close()
#	for currentKey, currentValue in dictConfigKeyValues.items() :
#		print(currentKey, currentValue)
	return dictConfigKeyValues

def getModificationDictionary(dictConfigKeyValues):
	p = re.compile('{.+}')
	pa = re.compile('[\d]+')
	pat = re.compile('[a-zA-Z]+')
	patt = re.compile('[-\.\d]+')
	patte = re.compile('[A-Z]')
	# get the element mass
	element_mass_list_dict = {}
	element_pert_list_dict = {}
	element_mass_dict = {}
	element_dict = {}
	compound_list_dict = {}
	element_modification_list_dict = {}
	for e_key, e_value in dictConfigKeyValues.iteritems():
		if e_key.startswith('[Peptide_Identification]Element_Masses'):
			m = p.search(e_key)
			e_str = m.group(0)[1:-1]
			element_mass_list_dict[e_str] = e_value
		if e_key.startswith('[Peptide_Identification]Element_Percent'):
			m = p.search(e_key)
			e_str = m.group(0)[1:-1]
			element_pert_list_dict[e_str] = e_value
		if e_key.startswith('[Peptide_Identification]Element_List'):
			element_list = pat.findall(e_value)
			idx = 0
			for e in element_list:
				element_dict[idx] = e
				idx += 1
		if e_key.startswith('[Peptide_Identification]Residue'):
			m = p.search(e_key)
			e_str = m.group(0)[1:-1]
			compound_list_dict[e_str] = map(int, patt.findall(e_value))
		if e_key.startswith('[Peptide_Identification]PTM'):
			m = p.search(e_key)
			e_str = m.group(0)[1:-1]
			if len(e_str) > 1:
				e_str = e_str[0]
			element_list = patte.findall(e_value)
			for e_element in element_list:
				if e_element in element_modification_list_dict:
					element_modification_list_dict[e_element].append(e_str)
				else:
					element_modification_list_dict[e_element] = [e_str]
			
			
	
	for e_key, e_value in element_mass_list_dict.iteritems():
		mass_list = map(float, patt.findall(e_value))
		pert_list = map(float, patt.findall(element_pert_list_dict[e_key]))
		max_value = max(pert_list)
		max_index = pert_list.index(max_value)
		element_mass_dict[e_key] = mass_list[max_index]
	
	# get the mass for all compound
	compound_mass_dict = {}
	for e_key, e_value in compound_list_dict.iteritems():
		'''
		if e_key.isalpha():
			continue
		'''
		mass_float = 0
		for e_idx, e_num in enumerate(e_value):
			e = element_dict[e_idx]
			mass_float += element_mass_dict[e] * float(e_num)
		compound_mass_dict[e_key] = mass_float
	
	return (compound_mass_dict, element_modification_list_dict)
