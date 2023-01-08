#ifndef MS2SCANVECTOR_H
#define MS2SCANVECTOR_H

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "ms2scan.h"
#include "tokenvector.h"
#include "proNovoConfig.h"
#include "proteindatabase.h"
#include "MVH.h"
#include "CometSearchMod.h"

#define ZERO            0.00000001
#define PEPTIDE_ARRAY_SIZE  2000000
#define PEPTIDE_ARRAY_SIP_SIZE  200000
#define TASKWAIT_SIZE 64
#define TASKRESUME_SIZE 32
#define TASKPEPTIDE_ARRAY_SIZE  20000

using namespace std;

class MS2ScanVector {
	// All MS2 scans to be scored
	// the MS2 scans are sorted by their precursor masses
	vector<MS2Scan *> vpAllMS2Scans;
	// the precursor mass of these MS2 scans
	// this is used for quickly inserting a peptide into corrent MS2 scans
	// vpAllMS2Scans and vpPrecursorMasses are in the same order
	vector<double> vpPrecursorMasses;
	//vector <int> mass_w; // mass window
	string sFT2Filename;    // the FT2 filename
	string sOutputFile; // the output file name
	string sConfigFile; // the configure file name

	// this should be moved to the Peptide class or make it a static member in the Config class
	map<char, double> mapResidueMass; // mass except N and C termini;

	void saveScan(MS2Scan * pMS2Scan);
	bool isMS1HighRes(const string & target);
	bool ChargeDetermination(const vector<double> & vdAllmz, double pmz); //return true if parent_charge should be 1
	bool bScreenOutput; // if true, allows standard output
	static bool mygreater(double i, double j);
	static bool myless(MS2Scan * pMS2Scan1, MS2Scan * pMS2Scan2);
	static bool mylessScanId(MS2Scan * pMS2Scan1, MS2Scan * pMS2Scan2);

	// find every MS2 scan whose precursor mass matches peptide mass
	bool assignPeptides2Scans(Peptide * currentPeptide);
	pair<int, int> GetRangeFromMass(double lb, double ub);

	// regular search functions
	void preProcessAllMs2Mvh(); // pre-process all MS2 scans before mvh
	void searchDatabaseMvh(); // search all ms2 scans against the protein list using mvh
	void processPeptideArrayMvh(vector<Peptide*>& vpPeptideArray); // process peptide array using mvh score
	void searchDatabaseMvhTask(); // search all ms2 scans against the protein list using mvh, task version
	void processPeptideArrayMvhTask(vector<Peptide*>& vpPeptideArray, omp_lock_t * pLck); // process peptide array using mvh score, task version
	void postProcessAllMs2WdpXcorr(); // post processing after mvh scoring
	void postProcessAllMs2Wdp(); // post processing using wdp scoring
	void postProcessAllMs2Xcorr(); // post processing using xcorr scoring

	// sip search functions
	void preProcessAllMs2WdpSip();  // pre-process all MS2 scans before wdp for sip search
	void searchDatabaseWdpSip(); // search all ms2 scans against the protein list using wdp for sip search
	void processPeptideArrayWdpSip(vector<Peptide*>& vpPeptideArray); // process peptide array using wdp score for sip search
	void postProcessAllMs2MvhXcorr(); // post processing after wdp scoring
	void postProcessAllMs2XcorrSip(); // post processing using xcorr scoring for sip search
	void postProcessAllMs2MvhSip(); // post processing using mvh scoring for sip search

	// write results to a SIP file
	// the SIP file will the same base file name as sFT2Filename
	// change the extension filename to ".SIP"
	void setOutputFile(const string & sFT2FilenameInput, const string & sOutputDirectory);
	void writeOutputEnsemble();
	void GetAllRangeFromMass(double dPeptideMass, vector<pair<int, int> > & vpPeptideMassRanges);
	string ParsePath(string sPath);

public:
	MS2ScanVector(const string & sFT2FilenameInput, const string & sOutputDirectory, const string & sConfigFilename,
			bool bScreenOutput);
	~MS2ScanVector();

	// Populate vpAllMS2Scans from the input FT2 file
	// Determine the charge state of every scan by calling the function MS2Scan::isSinglyCharged().
	// Creat a +2 scan and a +3 scan for an unknown multipe charged scan.
	// Return false if there is a problem with the file
	bool loadMassData();
	bool ReadFT2File();    //Read FT2 files
	bool ReadMzmlFile();   //Read mzML files

	// begin the database searching
	void startProcessingMvh(); // start functions to process the loaded FT2 file using mvh as the prime score
	void startProcessingMvhTask(); // start functions to process the loaded FT2 file using mvh as the prime score, task mode
	void startProcessingWdpSip(); // start functions to process the loaded FT2 file with WDP as prime score without tasking

	// variables for the MVH thread
	vector<vector<double> *> _ppdAAforward;
	vector<vector<double> *> _ppdAAreverse;
	vector<vector<double> *> psequenceIonMasses;
	vector<vector<char> *> pSeqs;
	int num_max_threads;
	void preMvh();
	void postMvh();
	// variables for the Xcorr thread
	vector<bool *> vpbDuplFragmentGlobal;
	vector<double *> v_pdAAforwardGlobal;
	vector<double *> v_pdAAreverseGlobal;
	vector<unsigned int ***> v_uiBinnedIonMassesGlobal;
	vector<vector<bool> > vvpbDuplFragmentGlobal;
	vector<vector<double> > vvdBinnedIonMassesGlobal;
	vector<vector<int> > vvdBinGlobal;
	void preXcorr();
	void postXcorr();
	int iOpenMPTaskNum;
	omp_lock_t lckOpenMpTaskNum;
	omp_lock_t lckOpenMpTaskNumHalfed;

};

#endif // MS2SCANVECTOR_H
