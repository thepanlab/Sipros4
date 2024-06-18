#include "peptide.h"

Peptide::Peptide()
{
	iPeptideLength = 0;
}

Peptide::~Peptide()
{
}

void Peptide::setPeptide(const string &sPeptide, const string &sOriginalPeptide,
						 const string &sProteinName, const int &ibeginPos,
						 const double &dPeptideMass, const char &cIdentifyPrefix, const char &cIdentifySuffix,
						 const char &cOriginalPrefix, const char &cOriginalSuffix)
{
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

void Peptide::calculateExpectedFragmentsSIP(const string &sNewPeptide)
{
	ProNovoConfig::configIsotopologue.computeProductIon(sNewPeptide,
														vvdYionMass, vvdYionProb, vvdBionMass, vvdBionProb);
	vdBionMasses.clear();
	vdBionMasses.reserve(sNewPeptide.length());
	vdYionMasses.clear();
	vdYionMasses.reserve(sNewPeptide.length());

	for (size_t i = 0; i < vvdYionProb.size(); i++)
	{
		// Find the index of the maximum value in vvdYionProb[i]
		auto maxElement = std::max_element(vvdYionProb[i].begin(), vvdYionProb[i].end());
		int maxIndex = std::distance(vvdYionProb[i].begin(), maxElement);
		vdYionMasses.push_back(vvdYionMass[i][maxIndex]);
	}
	for (size_t i = 0; i < vvdBionProb.size(); i++)
	{
		auto maxElement = std::max_element(vvdBionProb[i].begin(), vvdBionProb[i].end());
		int maxIndex = std::distance(vvdBionProb[i].begin(), maxElement);
		vdBionMasses.push_back(vvdBionMass[i][maxIndex]);
	}
}

void Peptide::calculateExpectedFragments(const string &sNewPeptide, const map<char, double> &mapResidueMass)
{
	int i;
	double dMass = 0;
	map<char, double>::const_iterator iter;
	vdBionMasses.clear();
	vdBionMasses.reserve(sNewPeptide.length());
	vdYionMasses.clear();
	vdYionMasses.reserve(sNewPeptide.length());

	if (sNewPeptide[0] != '[')
		cerr << "ERROR: peptide sequence must start with '[' as N terminus; Invalid sequence = " << sNewPeptide << endl;
	for (i = 1; i < (int)sNewPeptide.length(); ++i)
	{
		if (sNewPeptide[i] == ']')
			// hit the C terminus and break out of the loop
			break;
		iter = mapResidueMass.find(sNewPeptide[i]);
		if (iter == mapResidueMass.end())
			cerr << "WARNING: Residue " << sNewPeptide[i] << " Peptide " << sNewPeptide << " is not defined in the config." << endl;
		else if (isalpha(sNewPeptide[i]))
		{
			// this is an amino acid residue
			dMass = dMass + iter->second;
			vdBionMasses.push_back(dMass);
		}
		else
		{
			// this is a PTM and its mass should be added to the proceding residue
			//	cout << "PTM " << sNewPeptide[i] << " # " << iter->second << endl;
			dMass = dMass + iter->second;
			if (vdBionMasses.size() > 0)
				vdBionMasses.back() += iter->second;
			// if not, this symbol represents a PTM to the N terminus, whose mass will be added to the next residue.
		}
	}
	for (i = i + 1; i < (int)sNewPeptide.length(); ++i)
		if (!isalpha(sNewPeptide[i]))
		{
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

void Peptide::calculateIsotope(const string &sNewPeptide, const map<char, double> &mapResidueMass)
{
	ProNovoConfig::configIsotopologue.computeProductIon(sNewPeptide,
														vvdYionMass, vvdYionProb, vvdBionMass, vvdBionProb);
}

void Peptide::preprocessing(bool isMS2HighRes, const map<char, double> &mapResidueMass)
{
	int i;
	string sNewPeptide;
	for (i = 0; i < (int)sPeptide.length(); ++i)
		if (isalpha(sPeptide[i]))
			iPeptideLength = iPeptideLength + 1;

	sNewPeptide = neutralLossProcess(sPeptide);
	if (isMS2HighRes)
	{
		calculateIsotope(sNewPeptide, mapResidueMass); // just for weightsum
													   // calculateExpectedFragments(mapResidueMass); // just for ranksum
	}
	else if (ProNovoConfig::getSearchType() == "SIP")
		calculateExpectedFragmentsSIP(sNewPeptide);
	else
		calculateExpectedFragments(sNewPeptide, mapResidueMass);
}

string Peptide::neutralLossProcess(const string &sCurrentPeptide)
{
	string sNewPeptide;
	vector<pair<string, string>> vpNeutralLossList;
	size_t i, listLength, pos;
	sNewPeptide = sCurrentPeptide;
	vpNeutralLossList = ProNovoConfig::getNeutralLossList();
	//   cout<<"!!!"<<endl;
	//   cout<< sNewPeptide <<endl;
	listLength = vpNeutralLossList.size();
	for (i = 0; i < listLength; i++)
	{
		pos = 0;
		pos = sNewPeptide.find(vpNeutralLossList[i].first, pos);
		while (pos != string::npos)
		{
			sNewPeptide.replace(pos, 1, vpNeutralLossList[i].second);
			pos = sNewPeptide.find(vpNeutralLossList[i].first, pos);
		}
		// cout<<vpNeutralLossList[i].first<<" : "<<vpNeutralLossList[i].second<<endl;
	}
	// cout<< sNewPeptide <<endl;
	return sNewPeptide;
}

// only for test
// void Peptide::print()
// {
// 	cout << "Protein name: " << sProteinName << endl;
// 	cout << "Protein seq: " << sOriginalPeptide << endl;
// }