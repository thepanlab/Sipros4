#ifndef MS2SCAN_H
#define MS2SCAN_H

#include <vector>
#include <string>

#include "peptide.h"
//#include "TandemMassSpectrum.h"

#define BIN_RES 1000
#define LOW_BIN_RES 10
#define TOP_N 5
#define TOPPICKNUM 20
#define ZERO 0.00000001
#define SMALLINCREMENT 0.00000001 // to solve the incorrect cut off

using namespace std;

class ProductIon // for each ion
{
    char cIonType;              // y or b
    int iIonNumber;             // 1: y1 or b1; 2: y2 or b2; ...
    int iCharge;                // charge state
    int iMostAbundantPeakIndex; // related to the item with highest intensity
    double dMostAbundantMass;   // related to the item with highest intensity
    double dMostAbundantMZ;     // related to the item with highest intensity
    double dMZError;
    double dMassError;                   // calculate based on dMZError
    double dScoreWeight;                 // related to mass error and intensity
    bool bComplementaryFragmentObserved; // if y_k and b_n-k co-exist, where n is # of residues
public:
    ProductIon();
    ~ProductIon();
    // set the basic info
    void setProductIon(char cIonTypeInput, int iIonNumberInput, int iChargeInput);
    // set observed info
    void setObservedInfo(double dMZErrorInput, double dWeightInput,
                         double dMostAbundantMZInput, int iMostAbundantPeakIndexInput);
    void setComplementaryFragmentObserved(bool bComplementaryFragmentObservedInput);
    char getIonType() { return cIonType; };
    int getIonNumber() { return iIonNumber; };
    int getCharge() { return iCharge; };
    double getMZError() { return dMZError; };
    double getMassError() { return dMassError; };
    double getScoreWeight() { return dScoreWeight; };
    double getMostAbundantMass() { return dMostAbundantMass; };
    double getMostAbundantMZ() { return dMostAbundantMZ; };
    int getMostAbundantPeakIndex() { return iMostAbundantPeakIndex; };
    bool getComplementaryFragmentObserved() { return bComplementaryFragmentObserved; };
};

/*
class ScanUnit
{
public:
    double intensity;
    bool match;
    ScanUnit(double inten, bool mat){intensity = inten; match=mat;};
};*/

class PeptideUnit
{
public:
    double dCalculatedParentMass;
    double dScore;
    string sIdentifiedPeptide;
    string sOriginalPeptide;
    string sProteinNames;
    string sScoringFunction;
    char cIdentifyPrefix;
    char cIdentifySuffix;
    char cOriginalPrefix;
    char cOriginalSuffix;

    void setPeptideUnitInfo(const Peptide *currentPeptide, const double &dScore, string sScoringFunction);
};

class alignas(64) MS2Scan
{
    int bin_res, iMaxMZ, iMinMZ;

    double dMassTolerance; //<FRAGMENT_IONS>
    double dProtonMass;    // proton mass

    string sScanType; // format: FT-MS1/FT-MS2@CID

    vector<double> vdpreprocessedMZ;
    vector<double> vdpreprocessedIntensity;
    vector<int> vipreprocessedCharge; // the value of vipreprocessedCharge is zero for low-resolution MS2Scan

    // vdMaxMzIntensity[i] is the maximum intensity at M/Z window of i plus and minus iMzRange
    vector<double> vdMaxMzIntensity;
    vector<double> vdMzIntensity;
    vector<double> vdHighIntensity; // thresholds

    vector<int> vbPeakPresenceBins;
    vector<pair<int, int>> vbPeakPresenceBins2D; // first is lower bounder, second is upper bound
    vector<int> viIntensityRank;                 // ranks of preprocessed intensities staring with 0

    void preprocessLowMS2();
    void preprocessHighMS2();
    void initialPreprocess();
    void sortPeakList(); // bubble sort the peak list by MZ in ascending order
    int getMaxValueIndex(const vector<double> &vdData);
    void normalizeMS2scan();
    void setIntensityThreshold();
    void filterMS2scan();
    static bool mygreater(double i, double j);
    void binCalculation();
    void binCalculation2D(); // replace binCalculation();
    void saveScore(const double &dScore, const Peptide *currentPeptide, vector<PeptideUnit *> &vpTopPeptides,
                   vector<double> &vdAllScores, string sScoreFunction = "WeightSum");
    static bool GreaterScore(PeptideUnit *p1, PeptideUnit *p2);
    void WeightCompare(const string &sPeptide, vector<bool> &vbFragmentZ2);
    bool searchMZ(const double &dTarget, int &iIndex4Found);   // corresponds to binCalculation()
    bool searchMZ2D(const double &dTarget, int &iIndex4Found); // corresponds to binCalculation2D()
    bool searchMZ2D(const double &dTarget, const double &dErrRange, int &iIndex4Found);
    // merge same peptide. If no same peptide return false, otherwise, return true
    bool mergePeptide(vector<PeptideUnit *> &vpTopPeptides, const string &sPeptide, const string &sProteinName);
    void sortPreprocessedIntensity(); // sort preprocessed Intensity;

    void cleanup();

    // static bool mySUGreater (ScanUnit s1, ScanUnit s2);

    bool binarySearch(const double &dTarget, const vector<double> &vdList, const double &dTolerance,
                      vector<int> &viIndex4Found);

    // build the map between y or b ions for observed intensity and related mass (only for one ion)
    bool findProductIon(const vector<double> &vdIonMass, // expected mass
                        const vector<double> &vdIonProb,
                        // expected intensity based on the summation of all related intensities
                        const int &iCharge,
                        double &dScoreWeight,
                        double &dAverageMZError,
                        double &dMostAbundantObservedMZ, // with highest intensity
                        int &iMostAbundantPeakIndex);    // start with 0

    bool findProductIonSIP(const vector<double> &vdIonMass, // expected mass
                           const vector<double> &vdIonProb,
                           // expected intensity based on the summation of all related intensities
                           const int &iCharge,
                           double &dScoreWeight,
                           double &dAverageMZError,
                           double &dMostAbundantObservedMZ, // with highest intensity
                           int &iMostAbundantPeakIndex);    // start with 0
public:
    MS2Scan();
    //    MS2Scan(const MS2Scan *& cMS2Scan);
    ~MS2Scan();

    string sFT2Filename;       // FT2file name
    int iParentChargeState;    // Parent ion charge state
    double dParentMZ;          // Parent ion M/Z
    double dParentNeutralMass; // Parent neutral mass
    int iScanId;               // Product ions in the scan

    vector<PeptideUnit *> vpWeightSumTopPeptides;

    // the scores of all scored peptides
    vector<double> vdWeightSumAllScores;

    int inumberofWeightSumScore;
    double dsumofWeightScore;
    double dsumofSquareWeightSumScore;

    vector<double> vdMZ;
    vector<double> vdIntensity;
    vector<int> viCharge; // the value of viCharge is zero for low-resolution MS2

    bool isMS1HighRes;            // Is the MS1 scan a high-resolution scan?
    bool isMS2HighRes;            // Is this MS2 scan a high-resolution scan?
    vector<Peptide *> vpPeptides; // current set of peptides to be scored
    bool bSetMS2Flag;             // false when preprocess fail on bad data

    // preprocess this scan, including
    // (1) remove noise peaks
    // (2) normalize intensity
    void preprocess();
    // only for test
    void print();
    // score all peptides matched to this scan
    void scorePeptides();
    void scorePeptidesLowMS2();
    void scorePeptidesHighMS2();

    // scoring functions
    void scoreWeightSum(Peptide *currentPeptide);
    void scoreRankSum(Peptide *currentPeptide);
    double scoreIntensity(const bool observed, const double realIntensity, const double expectedIntensity);
    // void scoreRankSum_test(Peptide * currentPeptide);
    void scoreWeightSumHighMS2(Peptide *currentPeptide);
    void scoreRankSumHighMS2(Peptide *currentPeptide);
    // normalize raw scores
    void postprocess();
    double CalculateRankSum(double r1, double n1, double n2);

    void setScanType(string sScanType) { this->sScanType = sScanType; };
    string getScanType() { return this->sScanType; };
};

#endif // MS2SCAN_H
