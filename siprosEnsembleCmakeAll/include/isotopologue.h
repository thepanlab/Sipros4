#ifndef ISOTOPOLOGUE_H
#define ISOTOPOLOGUE_H

#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <iomanip>
#include <algorithm>
#include <iostream>
#include "proNovoConfig.h"

using namespace std;

/**/
class IsotopeDistribution {
public:
	IsotopeDistribution();
	IsotopeDistribution(vector<double> vItsMass, vector<double> vItsComposition);
	~IsotopeDistribution();

	vector<double> vMass;
	vector<double> vProb;

	double getMostAbundantMass();
	double getAverageMass();
	double getLowestMass();

	void filterProbCutoff(double dProCutoff);

	// print out the isotoptic distribution
	// this is mainly used for debuging
	void print();

};

class Isotopologue {
public:
	Isotopologue();
	~Isotopologue();

	// setup all variables from configuration
	bool setupIsotopologue(const string & sTable, const string & AtomNameInput);

	// get the MostAbundant masses of  residues
	bool getSingleResidueMostAbundantMasses(vector<string> & vsResidues, vector<double> & vdMostAbundantMasses, double & dTerminusMassN,
			double & dTerminusMassC);

	// compute the mass of the most abundant isotopologue
	double computeMostAbundantMass(string sSequence);
	double computeAverageMass(string sSequence);
	double computeMonoisotopicMass(string sSequence);

	// compute isotoptic distributions for all product ions of a sequence
	// this isotopologue class is modified by only adding this function
	// The first dimension of vvdYionMass and vvdYionProb is from y1, y2, ...
	// The first dimension of vvdBionMass and vvdBionProb is from b1, b2, ...
	// the mass is calculated assuming cleavage of the peptide bond
	bool computeProductIon(string sSequence, vector<vector<double> > & vvdYionMass, vector<vector<double> > & vvdYionProb,
			vector<vector<double> > & vvdBionMass, vector<vector<double> > & vvdBionProb);

	// compute isotoptic distribution for an amino acid sequence
	bool computeIsotopicDistribution(string sSequence, IsotopeDistribution & myIsotopeDistribution);

	// compute isotoptic distribution for a given atomic composition,
	// which can be that of a residue's or a amino acid sequence's
	bool computeIsotopicDistribution(vector<int> AtomicComposition, IsotopeDistribution & myIsotopeDistribution);

	// compute the atomic composition for an amino acid sequence
	bool computeAtomicComposition(string sSequence, vector<int> & myAtomicComposition);

private:
	// emass functions for IsotopeDistribution's arithmetic
	IsotopeDistribution sum(const IsotopeDistribution & distribution0, const IsotopeDistribution & distribution1);
	// functions for IsotopeDistribution's arithmetic
	IsotopeDistribution sum_backup(const IsotopeDistribution & distribution0, const IsotopeDistribution & distribution1);
	IsotopeDistribution multiply(const IsotopeDistribution & distribution0, int count);
	void shiftMass(IsotopeDistribution & distribution0, double dMass);

	// implementation of max and min
	double maximum(double a, double b) {
		return (a > b) ? a : b;
	}
	double minimum(double a, double b) {
		return (a < b) ? a : b;
	}

	// when two peaks have a mass difference less than the MassPrecision
	// they will be merged into one peak with their average mass and sum intensity
	const double MassPrecision; // 0.01

	// when a peak have a probability less than the ProbabilityCutoff
	// this peak will be ingnored, which makes the total probability space less than 1
	const double ProbabilityCutoff; // 1*10E-9

	// the name of atoms
	string AtomName;

	// the number of natural CHONPS and enriched CHONPS
	unsigned int AtomNumber;

	// variables for this isotopologue
	map<string, vector<int> > mResidueAtomicComposition;
	vector<IsotopeDistribution> vAtomIsotopicDistribution;
	map<string, IsotopeDistribution> vResidueIsotopicDistribution;

	// Sipros Ensemble
	// emass needs the mass to be one nucleus difference
	bool CheckMass(vector<double> & vdMass, vector<double> & vdNaturalCompositionTemp);
};

#endif //ISOTOPOLOGUE_H
