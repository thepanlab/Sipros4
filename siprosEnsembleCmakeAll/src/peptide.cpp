#include "peptide.h"

Peptide::Peptide() {
	iPeptideLength = 0;
	cOriginalPrefix = '-';
	cOriginalSuffix = '-';
	dscore = 0;
	dPeptideMass = 0;
	ibeginPos = 0;
	cIdentifyPrefix = '-';
	cIdentifySuffix = '-';
}

Peptide::~Peptide() {

}

void Peptide::setPeptide(const string& sPeptide, const string& sOriginalPeptide, const string& sProteinName, const int& ibeginPos, const double & dPeptideMass,
		const char & cIdentifyPrefix, const char & cIdentifySuffix, const char & cOriginalPrefix, const char & cOriginalSuffix) {
	this->sPeptide = sPeptide;
	this->sOriginalPeptide = sOriginalPeptide;
	this->sProteinName = sProteinName;
	this->ibeginPos = ibeginPos;
	this->dPeptideMass = dPeptideMass;
	this->cIdentifyPrefix = cIdentifyPrefix;
	this->cIdentifySuffix = cIdentifySuffix;
	this->cOriginalPrefix = cOriginalPrefix;
	this->cOriginalSuffix = cOriginalSuffix;
}

void Peptide::calculateExpectedFragments(const string & sNewPeptide, const map<char, double> & mapResidueMass) {
	int i;
	double dMass = 0;
	map<char, double>::const_iterator iter;
	vdBionMasses.clear();
	vdBionMasses.reserve(sNewPeptide.length());
	vdYionMasses.clear();
	vdYionMasses.reserve(sNewPeptide.length());

	if (sNewPeptide[0] != '[')
		cerr << "ERROR: peptide sequence must start with '[' as N terminus; Invalid sequence = " << sNewPeptide << endl;
	for (i = 1; i < (int) sNewPeptide.length(); ++i) {
		if (sNewPeptide[i] == ']')
			// hit the C terminus and break out of the loop
			break;
		iter = mapResidueMass.find(sNewPeptide[i]);
		if (iter == mapResidueMass.end())
			cerr << "WARNING: Residue " << sNewPeptide[i] << " Peptide " << sNewPeptide << " is not defined in the config." << endl;
		else if (isalpha(sNewPeptide[i])) {
			// this is an amino acid residue
			dMass = dMass + iter->second;
			vdBionMasses.push_back(dMass);
		} else {
			// this is a PTM and its mass should be added to the proceding residue
			//	cout << "PTM " << sNewPeptide[i] << " # " << iter->second << endl;
			dMass = dMass + iter->second;
			if (vdBionMasses.size() > 0)
				vdBionMasses.back() += iter->second;
			// if not, this symbol represents a PTM to the N terminus, whose mass will be added to the next residue.
		}
	}
	for (i = i + 1; i < (int) sNewPeptide.length(); ++i)
		if (!isalpha(sNewPeptide[i])) {
			// this is PTM to the C terminus
			iter = mapResidueMass.find(sNewPeptide[i]);
			if (iter == mapResidueMass.end())
				cerr << "WARNING: Residue " << sNewPeptide[i] << " Peptide " << sNewPeptide << " is not defined in the config." << endl;
			else
				dMass = dMass + iter->second;
		}
	dMass = dMass + ProNovoConfig::getTerminusMassC() + ProNovoConfig::getTerminusMassN();
	vdBionMasses.pop_back();
	for (i = vdBionMasses.size() - 1; i >= 0; --i)
		vdYionMasses.push_back(dMass - vdBionMasses[i]);
}

void Peptide::calculateIsotope(const string & sNewPeptide) {
	ProNovoConfig::configIsotopologue.computeProductIon(sNewPeptide, vvdYionMass, vvdYionProb, vvdBionMass, vvdBionProb);
}

void Peptide::preprocessing(bool isMS2HighRes, const map<char, double>& mapResidueMass) {
	int i;
	string sNewPeptide;
	for (i = 0; i < (int) sPeptide.length(); ++i)
		if (isalpha(sPeptide[i]))
			iPeptideLength = iPeptideLength + 1;

	sNewPeptide = neutralLossProcess(sPeptide);
	sNeutralLossPeptide = sNewPeptide;
	if (isMS2HighRes) {
		calculateIsotope(sNewPeptide); // just for weightsum
		//calculateExpectedFragments(mapResidueMass); // just for ranksum
	} else
		calculateExpectedFragments(sNewPeptide, mapResidueMass);
}

void Peptide::preprocessingMVH() {
	int i;
	iPeptideLength = 0;
	for (i = 0; i < (int) sPeptide.length(); ++i) {
		if (isalpha(sPeptide.at(i))) {
			iPeptideLength = iPeptideLength + 1;
		}
	}

	sNeutralLossPeptide = neutralLossProcess(sPeptide);
}

void Peptide::preprocessingSIP() {
	int i;
	iPeptideLength = 0;
	for (i = 0; i < (int) sPeptide.length(); ++i) {
		if (isalpha(sPeptide.at(i))) {
			iPeptideLength = iPeptideLength + 1;
		}
	}

	sNeutralLossPeptide = neutralLossProcess(sPeptide);
	// get isotope distribution
	calculateIsotope(sNeutralLossPeptide); // just for weightsum
	// re-scale the probability
	int iLen2 = 0, j = 0;
	double dMaxProb = 0;
	// y-ion
	int iLen1 = this->vvdYionProb.size();
	for (i = 0; i < iLen1; ++i) {
		iLen2 = vvdYionProb.at(i).size();
		dMaxProb = 0;
		for (j = 0; j < iLen2; ++j) {
			if (dMaxProb < vvdYionProb.at(i).at(j)) {
				dMaxProb = vvdYionProb.at(i).at(j);
			}
		}
		for (j = 0; j < iLen2; ++j) {
			vvdYionProb.at(i).at(j) /= dMaxProb;
		}
	}
	// b-ion
	iLen1 = this->vvdBionProb.size();
	for (i = 0; i < iLen1; ++i) {
		iLen2 = vvdBionProb.at(i).size();
		dMaxProb = 0;
		for (j = 0; j < iLen2; ++j) {
			if (dMaxProb < vvdBionProb.at(i).at(j)) {
				dMaxProb = vvdBionProb.at(i).at(j);
			}
		}
		for (j = 0; j < iLen2; ++j) {
			vvdBionProb.at(i).at(j) /= dMaxProb;
		}
	}
}

string Peptide::neutralLossProcess(const string& sCurrentPeptide) {
	string sNewPeptide;
	vector<pair<string, string> > vpNeutralLossList;
	size_t i, listLength, pos;
	sNewPeptide = sCurrentPeptide;
	vpNeutralLossList = ProNovoConfig::getNeutralLossList();
	//   cout<<"!!!"<<endl;
	//   cout<< sNewPeptide <<endl;
	listLength = vpNeutralLossList.size();
	for (i = 0; i < listLength; i++) {
		pos = 0;
		pos = sNewPeptide.find(vpNeutralLossList[i].first, pos);
		while (pos != string::npos) {
			sNewPeptide.replace(pos, 1, vpNeutralLossList[i].second);
			pos = sNewPeptide.find(vpNeutralLossList[i].first, pos);
		}
		//cout<<vpNeutralLossList[i].first<<" : "<<vpNeutralLossList[i].second<<endl;
	}
	//cout<< sNewPeptide <<endl;
	return sNewPeptide;
}

void Peptide::preprocessing(string & sPeptide, bool isMS2HighRes, const map<char, double> & mapResidueMass, vector<vector<double> > * vvdYionMass,
		vector<vector<double> > * vvdYionProb, vector<vector<double> > * vvdBionMass, vector<vector<double> > * vvdBionProb, vector<double> * pvdYionMass,
		vector<double> * pvdBionMass) {
	if (isMS2HighRes) {
		ProNovoConfig::configIsotopologue.computeProductIon(sPeptide, (*vvdYionMass), (*vvdYionProb), (*vvdBionMass), (*vvdBionProb));
	} else {
		calculateExpectedFragments(sPeptide, mapResidueMass, pvdYionMass, pvdBionMass);
	}
}

void Peptide::calculateExpectedFragments(const string & sNewPeptide, const map<char, double> & mapResidueMass, vector<double> * pvdYionMass,
		vector<double> * pvdBionMass) {
	int i;
	double dMass = 0;
	map<char, double>::const_iterator iter;
	pvdBionMass->clear();
	pvdBionMass->reserve(sNewPeptide.length());
	pvdYionMass->clear();
	pvdYionMass->reserve(sNewPeptide.length());

	if (sNewPeptide[0] != '[')
		cerr << "ERROR: peptide sequence must start with '[' as N terminus; Invalid sequence = " << sNewPeptide << endl;
	for (i = 1; i < (int) sNewPeptide.length(); ++i) {
		if (sNewPeptide[i] == ']') {	// hit the C terminus and break out of the loop
			break;
		}
		iter = mapResidueMass.find(sNewPeptide[i]);
		if (iter == mapResidueMass.end()) {
			cerr << "WARNING: Residue " << sNewPeptide[i] << " Peptide " << sNewPeptide << " is not defined in the config." << endl;
		} else if (isalpha(sNewPeptide[i])) {
			// this is an amino acid residue
			dMass = dMass + iter->second;
			pvdBionMass->push_back(dMass);
		} else {
			// this is a PTM and its mass should be added to the proceding residue
			//	cout << "PTM " << sNewPeptide[i] << " # " << iter->second << endl;
			dMass = dMass + iter->second;
			if (pvdBionMass->size() > 0) {
				pvdBionMass->back() += iter->second;
			}
			// if not, this symbol represents a PTM to the N terminus, whose mass will be added to the next residue.
		}
	}
	for (i = i + 1; i < (int) sNewPeptide.length(); ++i) {
		if (!isalpha(sNewPeptide[i])) {
			// this is PTM to the C terminus
			iter = mapResidueMass.find(sNewPeptide[i]);
			if (iter == mapResidueMass.end()) {
				cerr << "WARNING: Residue " << sNewPeptide[i] << " Peptide " << sNewPeptide << " is not defined in the config." << endl;
			} else {
				dMass = dMass + iter->second;
			}
		}
	}
	dMass = dMass + ProNovoConfig::getTerminusMassC() + ProNovoConfig::getTerminusMassN();
	pvdBionMass->pop_back();
	for (i = pvdBionMass->size() - 1; i >= 0; --i) {
		pvdYionMass->push_back(dMass - (*pvdBionMass)[i]);
	}
}

