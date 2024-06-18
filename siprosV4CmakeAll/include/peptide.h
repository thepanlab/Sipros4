#ifndef PEPTIDE_H
#define PEPTIDE_H

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "proNovoConfig.h"

using namespace std;

class alignas(64) Peptide
{
protected:
    // the sequence of the peptide
    string sPeptide;
    // the sequence of the original peptide
    string sOriginalPeptide;
    // the protein name
    string sProteinName;
    // the begining position of the original peptide on the the protein
    int ibeginPos;
    // the mass of the whole peptide
    double dPeptideMass;
    // score if it is in weightsum vector, the score is weightsum score, so do Xorr, Ranksome, and others.
    double dscore;
    // the length of the peptide
    int iPeptideLength;
    // prefix residue of identified peptide,  or "-"
    char cIdentifyPrefix;
    // suffix residue of identified peptide, or "-"
    char cIdentifySuffix;
    // prefix residue of original peptide,  or "-"
    char cOriginalPrefix;
    // suffix residue of original peptide,  or "-"
    char cOriginalSuffix;

public:
    // Low-res MS2 scan just need
    // the most abundant isotopic masses of predicted Y and B ions of the peptide
    vector<double> vdYionMasses;
    vector<double> vdBionMasses;

    // High-res MS2 scan need
    // the full isotopic distribution of predicted Y and B ions of the peptide
    vector<vector<double>> vvdYionMass;
    vector<vector<double>> vvdYionProb;
    vector<vector<double>> vvdBionMass;
    vector<vector<double>> vvdBionProb;

    Peptide();
    ~Peptide();

    void setPeptide(const string &sPeptide, const string &sOriginalPeptide,
                    const string &sProteinName, const int &ibeginPos,
                    const double &dPeptideMass, const char &cIdentifyPrefix, const char &cIdentifySuffix,
                    const char &cOriginalPrefix, const char &cOriginalSuffix);
    string getPeptideSeq() const { return sPeptide; };
    string getOriginalPeptideSeq() const { return sOriginalPeptide; };
    string getProteinName() const { return sProteinName; };
    int getBeginPosProtein() const { return ibeginPos; };
    double getPeptideMass() const { return dPeptideMass; };
    double getPeptideScore() const { return dscore; };
    int getPeptideLength() const { return iPeptideLength; };
    char getIdentifyPrefix() const { return cIdentifyPrefix; };
    char getIdentifySuffix() const { return cIdentifySuffix; };
    char getOriginalPrefix() const { return cOriginalPrefix; };
    char getOriginalSuffix() const { return cOriginalSuffix; };

    void setPeptideScore(double dscore) { this->dscore = dscore; };

    // calculate all expected y and b ions from sPeptide
    void calculateExpectedFragments(const string &sNewPeptide, const map<char, double> &mapResidueMass);

    void calculateExpectedFragmentsSIP(const string &sNewPeptide);

    void preprocessing(bool isMS2HighRes, const map<char, double> &mapResidueMass);

    // this function doesn't need mapResidueMass as input
    void calculateIsotope(const string &sNewPeptide, const map<char, double> &mapResidueMass);

    string neutralLossProcess(const string &sCurrentPeptide);
};

#endif // PEPTIDE_H
