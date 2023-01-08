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
#include <limits>
#include "omp.h"

using namespace std;

typedef long long INT64;

// To keep time information of functions.
#define CLOCKSTART INT64 mem_start = checkMemoryUsage(); double begin = omp_get_wtime(); cout<<"Currently in file: " << __FILE__ << " Function: "<< __FUNCTION__ << "()" << endl;
#define CLOCKSTOP INT64 mem_end = checkMemoryUsage(); double end = omp_get_wtime(); cout << "Function " << __FUNCTION__ << "() finished in " <<(end - begin)<< " Seconds." << endl << "Memory used: " << mem_end << " - " <<  mem_start << " = "<< mem_end - mem_start << " MB."<< endl << endl;

// Get the memory usage with a Linux kernel.
inline unsigned int checkMemoryUsage()
{
    // get KB memory into count
    unsigned int count=0;

    #if defined(__linux__)
    ifstream f("/proc/self/status"); // read the linux file
    while(!f.eof()){
        string key;
        f>>key;
        if(key=="VmData:"){     // size of data
            f>>count;
        break;
        }
    }
    f.close();
    #endif

    // return MBs memory (size of data)
    return (count/1024);
};


class Isotopologue;

//--------------Comet Begin------------
#define PROTON_MASS 1.00727646688
#define NUM_ION_SERIES 9
#define NUM_SP_IONS 200 // num ions for preliminary scoring

#define ION_SERIES_A                0
#define ION_SERIES_B                1
#define ION_SERIES_C                2
#define ION_SERIES_X                3
#define ION_SERIES_Y                4
#define ION_SERIES_Z                5

#define SPARSE_MATRIX_SIZE          100
#define FLOAT_ZERO                  1e-6     // 0.000001

#define MAX_FRAGMENT_CHARGE         5
#define MAX_PEPTIDE_LEN             150       // max # of AA for a peptide

struct Options             // output parameters
{
	int iNumStored;               // # of search results to store for xcorr analysis
	int iStartCharge;
	int iEndCharge;
	int iMaxFragmentCharge;
	int iMaxPrecursorCharge;
	int iRemovePrecursor;         // 0=no, 1=yes, 2=ETD precursors
	double dMinIntensity;
	double dRemovePrecursorTol;

	Options() {
		iNumStored = 0;               // # of search results to store for xcorr analysis
		iStartCharge = 0;
		iEndCharge = 0;
		iMaxFragmentCharge = 0;
		iMaxPrecursorCharge = 0;
		iRemovePrecursor = 0;         // 0=no, 1=yes, 2=ETD precursors
		dMinIntensity = 0;
		dRemovePrecursorTol = 0;
	}
};

struct IonInfo {
	int iNumIonSeriesUsed;
	int piSelectedIonSeries[NUM_ION_SERIES];
	int bUseNeutralLoss;
	int iIonVal[NUM_ION_SERIES];
	IonInfo() {
		bUseNeutralLoss = 0;
		iIonVal[ION_SERIES_A] = 0;
		iIonVal[ION_SERIES_B] = 1;
		iIonVal[ION_SERIES_C] = 0;
		iIonVal[ION_SERIES_X] = 0;
		iIonVal[ION_SERIES_Y] = 1;
		iIonVal[ION_SERIES_Z] = 0;
		iNumIonSeriesUsed = 2;
		piSelectedIonSeries[0] = 1;
		piSelectedIonSeries[1] = 4;
	}
};

struct PrecalcMasses {
	double dNtermProton;          // dAddNterminusPeptide + PROTON_MASS
	double dCtermOH2Proton;       // dAddCterminusPeptide + dOH2fragment + PROTON_MASS
	double dCtermOH2;       // dAddCterminusPeptide + dOHfragment
	double dOH2ProtonCtermNterm;  // dOH2parent + PROTON_MASS + dAddCterminusPeptide + dAddNterminusPeptide
	int iMinus17HighRes; // BIN'd value of mass(NH3)
	int iMinus17LowRes;
	int iMinus18HighRes; // BIN'd value of mass(H2O)
	int iMinus18LowRes;
	double dCO;
	double dNH3;
	double dNH2;
	double dCOminusH2;
};

#define AminoAcidMassesSize 256
// store the mass for different amino acids
class AminoAcidMasses{
public:
	static double dNULL;
	static double dERROR;
	vector<double> vdMasses;
	// double vdMasses[AminoAcidMassesSize];

	// construct function
	AminoAcidMasses();
	// clear vdMasses
	void clear();
	// reach an empty spot
	double end();
	// return the mass for the given amino acid
	double find(char _cAminoAcid);

	double operator[](char _cAminoAcid) const;

	double & operator[](char _cAminoAcid);

};
//--------------Comet End------------
class ProNovoConfig {
public:
	/*
	 * Sets up sessionwide configuration
	 * the configurations are loaded in to memory as static variables
	 */

	static bool setFilename(const string & sConfigFileName);

	static bool setWorkingDirectory(const string & sDirectoryName);

	static vector<pair<string, string> > getNeutralLossList() {
		return vpNeutralLossList;
	}

	static string getWorkingDirectory() {
		return sWorkingDirectory;
	}

	static string getSearchName() {
		return sSearchName;
	}

	static string getSearchType() {
		return sSearchType;
	}

	static char getSeparator();

	// retrieve <FASTA_Database>
	static string getFASTAfilename() {
		return sFASTAFilename;
	}

	// retrieve <Fragmentation_Method>
	static string getFragmentationMethod() {
		return sFragmentationMethod;
	}

	// retrieve the Minimum length of a peptide
	static int getMinPeptideLength() {
		return iMinPeptideLength;
	}

	// retrieve the max length of a peptide
	static int getMaxPeptideLength() {
		return iMaxPeptideLength;
	}

	// retrieve <Mass_Accuracy>	<Parent_Ion>
	static double getMassAccuracyParentIon() {
		return dMassAccuracyParentIon;
	}

	// retrieve <Mass_Accuracy>	<Fragment_Ions>
	static double getMassAccuracyFragmentIon() {
		return dMassAccuracyFragmentIon;
	}

	// retrieve <Parent_Mass_Windows>
	static vector<int> getParentMassWindows() {
		return viParentMassWindows;
	}

	static bool getPeptideMassWindows(double dPeptideMass,
			vector<pair<double, double> > & vpPeptideMassWindows);

	// retrieve <Max_PTM_Count>
	static int getMaxPTMcount() {
		return iMaxPTMcount;
	}

	// retrieve <Cleavage_Rules>
	static string getCleavageAfterResidues() {
		return sCleavageAfterResidues;
	}
	static string getCleavageBeforeResidues() {
		return sCleavageBeforeResidues;
	}
	static int getMaxMissedCleavages() {
		return iMaxMissedCleavages;
	}
	static bool getTestStartRemoval() {
		return bTestStartRemoval;
	}

	static bool getPTMinfo(map<string, string> & mPTMinfo);

	// retrieve <ATOM_ISOTOPIC_COMPOSITION>
	// the input character is the atom name CHONPS
	static bool getAtomIsotopicComposition(char cAtom,
			vector<double> & vdAtomicMass, vector<double> & vdComposition);

	static Isotopologue configIsotopologue;
	static vector<string> vsSingleResidueNames;
	static vector<double> vdSingleResidueMasses;

	static double getResidueMass(string sResidue);

	static double getTerminusMassN() {
		return dTerminusMassN;
	}
	static double getTerminusMassC() {
		return dTerminusMassC;
	}

	static double getProtonMass() {
		return 1.007276466;
	}

	static double getNeutronMass() {
		return 1.003355;
	}

	static double dnorm(double mean, double sd, double x) {
		double SQRT2PI = 2.506628;
		double a = (x - mean) / sd;
		return exp(-a * a / 2) / (SQRT2PI * sd);
	}

	static double pnorm(double dMean, double dStandardDeviation,
			double dRandomVariable) {
		double dZScore = (dRandomVariable - dMean) / dStandardDeviation;
		double dProbability = 0.5 * erfc(-dZScore / sqrt(2.0));
		return dProbability;
	}

	static double scoreError(double dMassError) {

		//	pnorm function
		return (1.0
				- pnorm(0, (getMassAccuracyFragmentIon() / 2), fabs(dMassError)))
				* 2.0;

		//	dnorm function
		//	return  ( dnorm( 0, (getMassAccuracyFragmentIon() / 2.0), fabs(dMassError) ) ) /
		//			( dnorm( 0, (getMassAccuracyFragmentIon() / 2.0), 0 ) )	;

		//  sigmoid function
		//	return ( 1/(1+exp(dMassError*600-3)));

	}

	//---------------Comet Begin---------------------
	static bool bXcorrEnable;
	static Options options;
	static double dInverseBinWidth; // this is used in BIN() many times so use inverse binWidth to do multiply vs. divide
	static double dOneMinusBinOffset; // this is used in BIN() many times so calculate once
	static IonInfo ionInformation;
	static int iXcorrProcessingOffset;
	static PrecalcMasses precalcMasses;
	static double dMaxMS2ScanMass;
	static double dMaxPeptideMass;
	// static map<char, double> pdAAMassFragment;
	static AminoAcidMasses pdAAMassFragment;
	static double dHighResFragmentBinSize;
	static double dHighResFragmentBinStartOffset;
	static double dLowResFragmentBinSize;
	static double dLowResFragmentBinStartOffset;
	static double dHighResInverseBinWidth;
	static double dLowResInverseBinWidth;
	static double dHighResOneMinusBinOffset;
	static double dLowResOneMinusBinOffset;
	static int iMaxPercusorCharge;
	//---------------Comet End-----------------------

	//---------------Myrimatch Begin-----------------
	static bool bMvhEnable;
	static double ClassSizeMultiplier;
	static int NumIntensityClasses;
	static int minIntensityClassCount;
	static double ticCutoffPercentage;
	static int MaxPeakCount;
	static int MinMatchedFragments;
	static double minObservedMz;
	static double maxObservedMz;
	//---------------Myrimatch End-------------------

	//---------------Sipros Score Begin--------------
	static bool bWeightDotSumEnable;
	static bool bLessIsotopicDistribution;
	static bool bMultiScores;
	static string sDecoyPrefix;
	static int INTTOPKEEP; // the top n PSM for calculation of other two scores
	static int iRank;
	//---------------Sipros Score End----------------
	static string sCleavageAfterResidues;
	static string sCleavageBeforeResidues;
	static int num_threads;

protected:
	ProNovoConfig();

private:

	static ProNovoConfig* ProNovoConfigSingleton;

	// the filename of the configuration file
	static string sFilename;

	string sSectionName;

	// the working directory
	static string sWorkingDirectory;

	// replace delimitor in a line
	static void replaceDelimitor(string & sLine, char cOldDelimitor,
			char cNewDelimitor);

	// variables from the PEPTIDE_IDENTIFICATION element
	static string sFASTAFilename;
	static string sFragmentationMethod;
	static string sSearchType;
	static string sSearchName;

	static int iMaxPTMcount;

	static int iMinPeptideLength;
	static int iMaxPeptideLength;

	static int iMaxMissedCleavages;
	static bool bTestStartRemoval;

	static double dMassAccuracyParentIon;
	static double dMassAccuracyFragmentIon;
	static vector<int> viParentMassWindows;

	static vector<pair<double, double> > vpPeptideMassWindowOffset;

	static vector<pair<string, string> > vpNeutralLossList;

	static double dTerminusMassN;
	static double dTerminusMassC;

	static string sElementList;

	// this is used to setup configIsotopologue
	// retrieve Elemental composition of amino acid residues
	static bool getResidueElementalComposition(
			string & sResidueElementalComposition);

	static bool calculatePeptideMassWindowOffset();

	/*
	 * New functions for cfg config files
	 */

	// all parameters in key-value pairs
	static map<string, string> mapConfigKeyValues;

	// parse the cfg file to populate mapConfigKeyValues
	bool parseConfigKeyValues();

	bool parseConfigLine(const string & sLine);

	// get the value of a key;
	// return false, if can't find the key in the mapConfigKeyValues
	static bool getConfigValue(string sConfigKey, string & sConfigValue);

	// get a set of key-value pairs, given a master key
	// return false, if can't find the key in the mapConfigKeyValues
	static bool getConfigMasterKeyValue(string sMasterKey,
			map<string, string> & mapKeyValueSet);

	// new version based on cfg config files
	bool getParameters();

	//handle neural loss
	static void NeutralLoss();
};

#endif /*PRONOVOCONFIG_H_*/
