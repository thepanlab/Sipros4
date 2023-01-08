#ifndef PROTEINDATABASE_H
#define PROTEINDATABASE_H

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "peptide.h"
#include "proNovoConfig.h"
#include "ptm.h"
//#include "alphabet.h"

//#define ORDERSTRING "ACDEFGHIJKLMNPQRSTVWY[]"

using namespace std;

class ProteinDatabase
{
protected:
    int imaxPTM;
    int iProteinId; // current protein id
    string scurrentProtein, scurrentProteinName, snextProteinName,
        slegalChar, sOriginalPeptide;
    bool bstayCurrentProtein;         // false if need to read the next protein
    bool bstayCurrentOriginalPeptide; // false if need the next original peptide
    ifstream db_stream;
    PTM_List ptmlist;                             // from configure file
    vector<vector<pair<string, double>>> ptm_map; // reorgnized based on ptmlist
    vector<int> vicleavageSite;
    int iclCheck;          // number of cleavage sites missed;
    int ibeginCleavagePos; // position in vicleavageSite
    double dOriginalPeptideMass;
    string orderstring;
    // all potential ptm position for the current peptide
    vector<pair<int, vector<pair<string, double>>>> ptm_position_all;
    int icurrentMaxPtm;             // maximum ptms allowed to happen in the current peptide
    int iPtmCount;                  // exact number of ptms in the current peptide
    vector<int> comb_order;         // info for combination;
    vector<int> ptm_order, ele_num; // based on given comb_order
    // Amino aa;
    vector<double> mass_map;    // Masses for each position in protein
    vector<double> msum_map;    // Mass sum ending at each position in protein
    int imutationPos;           // mutation position: -1 means no mutation
    int imutationOrder;         // indicate which mutation status it reaches. -1 no mutation
    int imutationCleavageCount; // for mutate one cleavage site
                                // bool bIsAfterCleavage; //true for after_cleavage site, false for before_cleavage site
    string sMutationOrder;      // string for potential mutations on the current position
    bool bLeftSubpeptide;       // true: left subpeptide exists, which is based on the new cleavage site
    bool bRightSubpeptide;      // true: right subpeptide exists, which is based on the new cleavage site
    string sLeftSubpeptide;
    string sRightSubpeptide;
    int iLeftSubPeptideBeginPos;  // beginning position of left subpeptide on the protein
    int iLeftSubPeptideEndPos;    // ending position of left subpeptide on the protein
    int iRightSubPeptideBeginPos; // beginning position of right subpeptide on the protein
    int iRightSubPeptideEndPos;   // ending position of right subpeptide on the protein
    char cLeftSubPeptideSuffix;
    char cRightSubPeptidePrefix;
    string sNonAfterCleavage;
    bool bScreenOutput; // if true, allows standard output
    // return false if there is no more peptide from this protein
    bool getNextPeptidePTM(Peptide *currentPeptide);
    // return false if no new protein
    bool getNextProtein();

    // make decoy protein
    bool bGetDecoy;
    bool getDecoyOrNextProtein();

    void RemoveIllegalResidue(string &seq);
    void Initial_PTM_Map();
    bool getNextOriginalPeptide(Peptide *originalPeptide);
    // In comparison of getNextOriginalPeptide, getNextOriginalPeptideForMutation
    // generates peptides with one more than allowed missed cleavage, those peptide can exist
    // only when the one cleavage was mutated into another non-mutated amino acid
    bool getNextOriginalPeptideForMutation(Peptide *originalPeptide);
    bool getNextPtmPeptide(Peptide *ptmPeptide);
    bool getNextMutationPeptide(Peptide *mutationPeptide);

    bool GenerateNextComb(vector<int> &comb_order, const int &total_num);
    bool GenerateNextPTM(vector<int> &ptm_order, const vector<int> &ele_num);
    void setProtein(); // set cleavage site;
    bool isCleavageSite(char c1, char c2);
    bool isPeptideLengthGood(int ipepLength);
    void initializePtmInfo(); // ptm info for the current original peptide
    void initializeCombOrderPTM();
    void initializePermutationPTM();
    bool getNextPeptideMutation(Peptide *currentPeptide); // get next peptide under mutation mode
    // return characters contained by sStringA, but not by sStringB
    string stringDifference(const string &sStringA, const string &sStringB);
    string mutationString(const string &sWholeString, const char &chCurrentRedidue);
    void mutatePeptide(const string &sOriginalPeptideContent, string &sMutatedPeptide);
    bool mutatePeptideLessCleavage(Peptide *mutationPeptide, const string &sOriginalPeptideContent, int iendCleavagePos, string &sMutatedPeptide);
    bool mutatePeptideEqualCleavage(Peptide *mutationPeptide, const string &sOriginalPeptideContent, int iendCleavagePos, string &sMutatedPeptide);
    bool mutatePeptideMoreCleavage(Peptide *mutationPeptide, const string &sOriginalPeptideContent, int iendCleavagePos, string &sMutatedPeptide);
    // check cleavage info for the new mutated peptide;
    bool verifyPeptide(const string &sMutatedPeptide);
    void setPeptideInfo(Peptide *thePeptide, const string &sIdentifyPeptide, const double &dMass);
    void setLeftSubPeptideInfo(Peptide *thePeptide);
    void setRightSubPeptideInfo(Peptide *thePeptide);
    void subPeptide(const string &sMutatedPeptide);

public:
    ProteinDatabase(bool bScreenOutput);
    //    ProteinDatabase(const ProteinDatabase& other);
    ~ProteinDatabase();

    // load the database from file ProNovoConfig::getFASTAfilename()
    // this function just opens the file. The file is read protein by protein
    void loadDatabase();
    // perform in silico digestion of this protein and return peptides one at a time
    // return true, if a peptide is generated from this protein and its information is returned by reference
    // return false, if there is no more peptide from databasefile
    // currentPeptide: a Peptide object that has been set to be this peptide by calling setPeptide()
    bool getNextPeptide(Peptide *currentPeptide);
    bool getFirstProtein();
};

#endif // PROTEINDATABASE_H
