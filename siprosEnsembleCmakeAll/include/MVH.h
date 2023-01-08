/*
 * MVH.h
 *
 *  Created on: May 23, 2016
 *      Author: xgo
 */

#ifndef SCORES_MVH_H_
#define SCORES_MVH_H_

class MS2Scan;

#include "ms2scan.h"
#include "proNovoConfig.h"
#include <bitset>
#include <numeric>
#include <limits>

using namespace std;

///Initialize the elemental masses
const double HYDROGEN_MONO = 01.00783;
const double OXYGEN_MONO = 15.99491;

///Initialize the water and ammonia masses
const double WATER_MONO = 2 * HYDROGEN_MONO + OXYGEN_MONO;
const double Proton = 1.00727646688;

enum FragmentTypes {
	FragmentType_A, FragmentType_B, FragmentType_C, FragmentType_X, FragmentType_Y, FragmentType_Z, FragmentType_Z_Radical, FragmentTypes_Size
};

class lnFactorialTable {
public:
	lnFactorialTable() {
		m_table.push_back(0);
		m_table.push_back(0);
	}

	~lnFactorialTable() {
		m_table.clear();
	}

	double operator[](int index) {
		// Is the table big enough?
		int maxIndex = m_table.size() - 1;
		if (index > maxIndex) {
			cerr << "error lnFactorialTable " << endl;
			exit(1);
		}
		return m_table.at(index);
	}

	void resize(int index) {
		int maxIndex = ((int)m_table.size()) - 1;
		if (index > maxIndex) {
			while (index > maxIndex) {
				m_table.push_back(m_table.at(maxIndex) + log((float) m_table.size()));
				++maxIndex;
			}
		}
	}

private:
	std::vector<double> m_table;
};

class MVH {
public:
	MVH();
	~MVH();

	static bitset<FragmentTypes_Size> fragmentTypes;
	static bool bUseSmartPlusThreeModel;
	static lnFactorialTable * lnTable;
	// for SIP mode, only larger than this cutoff, peak will be considered
	static double ProbabilityCutOff;

	static bool CalculateSequenceIons(string & currentPeptide, int maxIonCharge, bool useSmartPlusThreeModel, vector<double>* sequenceIonMasses,
			vector<double> * _pdAAforward, vector<double> * _pdAAreverse, vector<char> * seq);
	static bool CalculateSequenceIonsSIP(string & currentPeptide, int maxIonCharge, bool useSmartPlusThreeModel, vector<double>* sequenceIonMasses,
			vector<vector<double> > & vvdYionMass, vector<vector<double> > & vvdYionProb, vector<vector<double> > & vvdBionMass,
			vector<vector<double> > & vvdBionProb, vector<char> * seq);
	static bool destroyLnTable();
	static multimap<double, char>::iterator findNear(map<double, char> * peakData, double mz, double tolerance);
	static bool initialLnTable(int maxPeakBins);
	static double lnCombin(int n, int k);
	static bool Preprocess(MS2Scan * Spectrum, multimap<double, double> * IntenSortedPeakPreData);
	static bool ScoreSequenceVsSpectrum(string & currentPeptide, MS2Scan * Spectrum, vector<double>* sequenceIonMasses, vector<double> * _pdAAforward,
			vector<double> * _pdAAreverse, double & dMvh, vector<char> * seq);
	static bool ScoreSequenceVsSpectrumSIP(string & currentPeptide, MS2Scan * Spectrum, vector<double>* sequenceIonMasses,
			vector<vector<double> > & vvdYionMass, vector<vector<double> > & vvdYionProb, vector<vector<double> > & vvdBionMass,
			vector<vector<double> > & vvdBionProb, double & dMvh, vector<char> * seq);
};

#endif /* SCORES_MVH_H_ */
