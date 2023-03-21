
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

#ifdef Ticktock
#include "ticktock.h"
#endif

#define ZERO 0.00000001
#define PEPTIDE_ARRAY_SIZE 200000

using namespace std;

class MS2ScanVector
{
    // All MS2 scans to be scored
    // the MS2 scans are sorted by their precursor masses
    // the precursor mass of these MS2 scans
    // this is used for quickly inserting a peptide into corrent MS2 scans
    // vpAllMS2Scans and vpPrecursorMasses are in the same order
    vector<double> vpPrecursorMasses;
    // vector <int> mass_w; // mass window
    string sFT2Filename; // the FT2 filename
    string sOutputFile;  // the output file name
    string sConfigFile;  // the configure file name

    // this should be moved to the Peptide class or make it a static member in the Config class
    map<char, double> mapResidueMass; // mass except N and C termini;

    void saveScan(MS2Scan &pMS2Scan);
    bool isMS1HighRes(const string &target);
    bool ChargeDetermination(const vector<double> &vdAllmz, double pmz); // return true if parent_charge should be 1
    bool bScreenOutput;                                                  // if true, allows standard output
    static bool mygreater(double i, double j);
    static bool myless(const MS2Scan &pMS2Scan1, const MS2Scan &pMS2Scan2);
    static bool mylessScanId(MS2Scan *pMS2Scan1, MS2Scan *pMS2Scan2);
    void preProcessAllMS2(); // Preprocessing all MS2 scans by multi-threading
    void searchDatabase();   // Search all MS2 scans against the protein list by multi-threading
    // find every MS2 scan whose precursor mass matches peptide mass
    void assignPeptides2Scans(Peptide *currentPeptide);
    void processPeptideArray(vector<Peptide *> &vpPeptideArray);
    pair<int, int> GetRangeFromMass(double lb, double ub);
    // Postprocessing all MS2 scans' results by multi-threading
    void postProcessAllMS2();
    // write results to a SIP file
    // the SIP file will the same base file name as sFT2Filename
    // change the extension filename to ".SIP"
    void setOutputFile(const string &sFT2FilenameInput, const string &sOutputDirectory);
    void writeOutput();
    // calculate mean and standard deviation of scores of a ms2 scan
    void calculateMeanAndDeviation(int inumberScore, double dScoreSum, double dScoreSquareSum,
                                   double &dMean, double &dDeviation);
    void GetAllRangeFromMass(double dPeptideMass, vector<pair<int, int>> &vpPeptideMassRanges);
    string ParsePath(string sPath);

public:
    vector<MS2Scan *> vpAllMS2ScanPtrs;
    vector<MS2Scan> vpAllMS2Scans;
    MS2ScanVector(const string &sFT2FilenameInput, const string &sOutputDirectory,
                  const string &sConfigFilename, bool bScreenOutput);
    ~MS2ScanVector();

    // Populate vpAllMS2Scans from the input FT2 file
    // Determine the charge state of every scan by calling the function MS2Scan::isSinglyCharged().
    // Creat a +2 scan and a +3 scan for an unknown multipe charged scan.
    // Return false if there is a problem with the file
    bool loadFT2file();
    bool ReadFT2File();     // Read FT2 files
    void startProcessing(); // start functions to process the loaded FT2 file
};

#endif // MS2SCANVECTOR_H
