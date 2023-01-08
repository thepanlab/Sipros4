#ifndef PRONOVOCONFIG_H_
#define PRONOVOCONFIG_H_

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "isotopologue.h"
#include "omp.h"

using namespace std;

class Isotopologue;

class ProNovoConfig
{
public:
	/*
	 * Sets up sessionwide configuration
	 * the configurations are loaded in to memory as static variables
	 */

	static bool setFilename(const string &sConfigFileName);

	static bool setWorkingDirectory(const string &sDirectoryName);

	static vector<pair<string, string>> getNeutralLossList()
	{
		return vpNeutralLossList;
	}

	static string getWorkingDirectory()
	{
		return sWorkingDirectory;
	}

	static string getSearchName()
	{
		return sSearchName;
	}

	static string getSearchType()
	{
		return sSearchType;
	}

	static char getSeparator();

	// retrieve <FASTA_Database>
	static string getFASTAfilename()
	{
		return sFASTAFilename;
	}

	// retrieve <Fragmentation_Method>
	static string getFragmentationMethod()
	{
		return sFragmentationMethod;
	}

	// retrieve the Minimum length of a peptide
	static int getMinPeptideLength()
	{
		return iMinPeptideLength;
	}

	// retrieve the max length of a peptide
	static int getMaxPeptideLength()
	{
		return iMaxPeptideLength;
	}

	// retrieve <Mass_Accuracy>	<Parent_Ion>
	static double getMassAccuracyParentIon()
	{
		return dMassAccuracyParentIon;
	}

	// retrieve <Mass_Accuracy>	<Fragment_Ions>
	static double getMassAccuracyFragmentIon()
	{
		return dMassAccuracyFragmentIon;
	}

	// retrieve <Parent_Mass_Windows>
	static vector<int> getParentMassWindows()
	{
		return viParentMassWindows;
	}

	static bool getPeptideMassWindows(double dPeptideMass, vector<pair<double, double>> &vpPeptideMassWindows);

	// retrieve <Max_PTM_Count>
	static int getMaxPTMcount()
	{
		return iMaxPTMcount;
	}

	// retrieve <Cleavage_Rules>
	static string getCleavageAfterResidues()
	{
		return sCleavageAfterResidues;
	}
	static string getCleavageBeforeResidues()
	{
		return sCleavageBeforeResidues;
	}
	static int getMaxMissedCleavages()
	{
		return iMaxMissedCleavages;
	}
	static bool getTestStartRemoval()
	{
		return bTestStartRemoval;
	}

	static bool getPTMinfo(map<string, string> &mPTMinfo);

	// retrieve <ATOM_ISOTOPIC_COMPOSITION>
	// the input character is the atom name CHONPS
	static bool getAtomIsotopicComposition(
		char cAtom,
		vector<double> &vdAtomicMass,
		vector<double> &vdComposition);

	static Isotopologue configIsotopologue;
	static vector<string> vsSingleResidueNames;
	static vector<double> vdSingleResidueMasses;

	static double getResidueMass(string sResidue);

	static double getTerminusMassN()
	{
		return dTerminusMassN;
	}
	static double getTerminusMassC()
	{
		return dTerminusMassC;
	}

	static double getProtonMass()
	{
		return 1.007276466;
	}

	static double getNeutronMass()
	{
		return neutronMass;
	}

	static double dnorm(double mean, double sd, double x)
	{
		double SQRT2PI = 2.506628;
		double a = (x - mean) / sd;
		return exp(-a * a / 2) / (SQRT2PI * sd);
	}

	static double pnorm(double dMean, double dStandardDeviation, double dRandomVariable)
	{
		double dZScore = (dRandomVariable - dMean) / dStandardDeviation;
		double dProbability = 0.5 * erfc(-dZScore / sqrt(2.0));
		return dProbability;
	}

	static double scoreError(double dMassError)
	{

		//	pnorm function
		return (1.0 - pnorm(0, (getMassAccuracyFragmentIon() / 2), fabs(dMassError))) * 2.0;

		//	dnorm function
		//	return  ( dnorm( 0, (getMassAccuracyFragmentIon() / 2.0), fabs(dMassError) ) ) /
		//			( dnorm( 0, (getMassAccuracyFragmentIon() / 2.0), 0 ) )	;

		//  sigmoid function
		//	return ( 1/(1+exp(dMassError*600-3)));
	}

	static string &getSetSIPelement() { return SIPelement; }
	static double &getSetMinValue() { return minValue; }
	static double &getSetFold() { return fold; }
	// compute deduction coefficient in score function
	// only suitbale for carbon and nitrogen SIP now
	static void setDeductionCoefficient();
	// get deduction coefficient in score function
	static double getDeductionCoefficient() { return deductionCoefficient; }

protected:
	ProNovoConfig();

private:
	static ProNovoConfig *ProNovoConfigSingleton;

	// the filename of the configuration file
	static string sFilename;

	string sSectionName;

	// the working directory
	static string sWorkingDirectory;

	// replace delimitor in a line
	static void replaceDelimitor(string &sLine, char cOldDelimitor, char cNewDelimitor);

	// variables from the PEPTIDE_IDENTIFICATION element
	static string sFASTAFilename;
	static string sFragmentationMethod;
	static string sSearchType;
	static string sSearchName;

	static int iMaxPTMcount;

	static int iMinPeptideLength;
	static int iMaxPeptideLength;

	static string sCleavageAfterResidues;
	static string sCleavageBeforeResidues;
	static int iMaxMissedCleavages;
	static bool bTestStartRemoval;

	static double dMassAccuracyParentIon;
	static double dMassAccuracyFragmentIon;
	static vector<int> viParentMassWindows;

	static vector<pair<double, double>> vpPeptideMassWindowOffset;

	static vector<pair<string, string>> vpNeutralLossList;

	static double dTerminusMassN;
	static double dTerminusMassC;

	static string sElementList;

	static string SIPelement;
	// for deductionCoefficient compute in SIP search
	static double neutronMass, deductionCoefficient, minValue, fold;

	// this is used to setup configIsotopologue
	// retrieve Elemental composition of amino acid residues
	static bool getResidueElementalComposition(string &sResidueElementalComposition);

	static bool calculatePeptideMassWindowOffset();

	/*
	 * New functions for cfg config files
	 */

	// all parameters in key-value pairs
	static map<string, string> mapConfigKeyValues;

	// parse the cfg file to populate mapConfigKeyValues
	bool parseConfigKeyValues();

	bool parseConfigLine(const string &sLine);

	// get the value of a key;
	// return false, if can't find the key in the mapConfigKeyValues
	static bool getConfigValue(string sConfigKey, string &sConfigValue);

	// get a set of key-value pairs, given a master key
	// return false, if can't find the key in the mapConfigKeyValues
	static bool getConfigMasterKeyValue(string sMasterKey, map<string, string> &mapKeyValueSet);

	// new version based on cfg config files
	bool getParameters();

	// handle neural loss
	static void NeutralLoss();
};

#endif /*PRONOVOCONFIG_H_*/
