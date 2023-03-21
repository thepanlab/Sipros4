#include "isotopologue.h"

IsotopeDistribution::IsotopeDistribution()
{
}

IsotopeDistribution::IsotopeDistribution(const vector<double> &vItsMass, const vector<double> &vItsProb)
	: vMass(vItsMass), vProb(vItsProb)
{
}

void IsotopeDistribution::print()
{
	cout << "Mass " << '\t' << "Inten" << endl;
	for (unsigned int i = 0; i < vMass.size(); i++)
	{
		cout << setprecision(8) << vMass[i] << '\t' << vProb[i] << endl;
	}
}

double IsotopeDistribution::getMostAbundantMass()
{
	double dMaxProb = 0;
	double dMass = 0;
	for (unsigned int i = 0; i < vMass.size(); ++i)
	{
		if (dMaxProb < vProb[i])
		{
			dMaxProb = vProb[i];
			dMass = vMass[i];
		}
	}
	return dMass;
}

double IsotopeDistribution::getAverageMass()
{
	double dSumProb = 0;
	double dSumMass = 0;
	for (unsigned int i = 0; i < vMass.size(); ++i)
	{
		dSumProb = dSumProb + vProb[i];
		dSumMass = dSumMass + vProb[i] * vMass[i];
	}

	if (dSumProb <= 0)
		return 1.0;

	return (dSumMass / dSumProb);
}

void IsotopeDistribution::filterProbCutoff(double dProbCutoff)
{
	vector<double> vMassCopy = vMass;
	vector<double> vProbCopy = vProb;
	vMass.clear();
	vProb.clear();
	for (unsigned int i = 0; i < vProbCopy.size(); ++i)
	{
		if (vProbCopy[i] >= dProbCutoff)
		{
			vMass.push_back(vMassCopy[i]);
			vProb.push_back(vProbCopy[i]);
		}
	}
}

double IsotopeDistribution::getLowestMass()
{
	return *min_element(vMass.begin(), vMass.end());
}

Isotopologue::Isotopologue() : MassPrecision(0.01), ProbabilityCutoff(0.000000001)
{
}

Isotopologue::~Isotopologue()
{
	// destructor
}

bool Isotopologue::setupIsotopologue(const string &sTable, const string &AtomNameInput)
{
	// CHONPS
	AtomName = AtomNameInput;
	AtomNumber = AtomName.size();

	istringstream issStream(sTable);
	string sResidue;
	vector<int> viAtomVector;
	int iNumber;
	unsigned int i;

	// parse out the RESIDUE_ATOMIC_COMPOSITION table
	while (!(issStream.eof()))
	{
		// each row is expected to start with the residue name, following by 6 numbers for the natuual CHONPS
		issStream >> sResidue;
		if (sResidue == "")
			continue;
		viAtomVector.clear();
		viAtomVector.reserve(AtomNumber);
		for (i = 0; i < AtomNumber; ++i)
		{
			if (issStream.eof())
			{
				// this row doesn't have 6 fields
				cerr << "ERROR:  the RESIDUE_ATOMIC_COMPOSITION table in ProNovoConfig is not correct!" << endl;
				return false;
			}
			issStream >> iNumber;
			viAtomVector.push_back(iNumber);
		}
		// add this row into the mResidueAtomicComposition table
		mResidueAtomicComposition[sResidue] = viAtomVector;
		sResidue = "";
	}

	// push 6 empty IsotopeDistributions into vAtomIsotopicDistribution
	vAtomIsotopicDistribution.reserve(AtomNumber);
	for (i = 0; i < (AtomNumber); ++i)
	{
		IsotopeDistribution TempDistribution;
		vAtomIsotopicDistribution.push_back(TempDistribution);
	}

	// variables to be passed as reference to ProNovoConfig::getAtomIsotopicComposition
	// to receive its return value
	vector<double> vdMassTemp;
	vector<double> vdNaturalCompositionTemp;

	// the isotopic distribution is pushed into vAtomIsotopicDistribution in the order of
	// natural CHONPS
	for (i = 0; i < AtomName.size(); ++i)
	{
		if (!ProNovoConfig::getAtomIsotopicComposition(
				AtomName[i],
				vdMassTemp,
				vdNaturalCompositionTemp))
		{
			cerr << "ERROR: cannot retrieve isotopic composition for atom " << AtomName[i] << " from ProNovoConfig" << endl;
			return false;
		}
		vAtomIsotopicDistribution[i].vMass = vdMassTemp;
		vAtomIsotopicDistribution[i].vProb = vdNaturalCompositionTemp;
	}

	// calculate Isotopic distribution for all residues
	map<string, vector<int>>::iterator ResidueIter;
	IsotopeDistribution tempIsotopeDistribution;
	for (ResidueIter = mResidueAtomicComposition.begin(); ResidueIter != mResidueAtomicComposition.end(); ResidueIter++)
	{
		if (!computeIsotopicDistribution(ResidueIter->second, tempIsotopeDistribution))
		{
			cerr << "ERROR: cannot calculate the isotopic distribution for residue " << ResidueIter->first << endl;
			return false;
		}

		vResidueIsotopicDistribution[ResidueIter->first] = tempIsotopeDistribution;

		//		cout << "Residue " << ResidueIter->first << endl;
		//		tempIsotopeDistribution.print();
	}

	return true;
}

double Isotopologue::computeMostAbundantMass(string sSequence)
{
	IsotopeDistribution tempIsotopeDistribution;
	if (!computeIsotopicDistribution(sSequence, tempIsotopeDistribution))
		return 0;
	else
		return tempIsotopeDistribution.getMostAbundantMass();
}

double Isotopologue::computeAverageMass(string sSequence)
{
	IsotopeDistribution tempIsotopeDistribution;
	if (!computeIsotopicDistribution(sSequence, tempIsotopeDistribution))
		return 0;
	else
		return tempIsotopeDistribution.getAverageMass();
}

double Isotopologue::computeMonoisotopicMass(string sSequence)
{
	IsotopeDistribution tempIsotopeDistribution;
	if (!computeIsotopicDistribution(sSequence, tempIsotopeDistribution))
		return 0;
	else
		return tempIsotopeDistribution.getLowestMass();
}

bool Isotopologue::getSingleResidueMostAbundantMasses(vector<string> &vsResidues, vector<double> &vdMostAbundantMasses, double &dTerminusMassN, double &dTerminusMassC)
{
	vsResidues.clear();
	vdMostAbundantMasses.clear();

	map<string, IsotopeDistribution>::iterator ResidueIter;
	string sCurrentResidue;
	IsotopeDistribution currentDistribution;
	double dCurrentMostAbundantMasses;

	// for single amino acid

	for (ResidueIter = vResidueIsotopicDistribution.begin(); ResidueIter != vResidueIsotopicDistribution.end(); ResidueIter++)
	{
		sCurrentResidue = ResidueIter->first;
		currentDistribution = ResidueIter->second;
		dCurrentMostAbundantMasses = currentDistribution.getMostAbundantMass();
		if (sCurrentResidue.size() == 1)
		{
			vsResidues.push_back(sCurrentResidue);
			vdMostAbundantMasses.push_back(dCurrentMostAbundantMasses);
			//			cout << sCurrentResidue << "  " << dCurrentMostAbundantMasses << endl;
		}
		else if (sCurrentResidue == "NTerm" || sCurrentResidue == "Nterm")
		{
			dTerminusMassN = dCurrentMostAbundantMasses;
		}
		else if (sCurrentResidue == "CTerm" || sCurrentResidue == "Cterm")
		{
			dTerminusMassC = dCurrentMostAbundantMasses;
		}
		else
		{
			cerr << "ERROR: Cannot recognize the configuration for residue " << sCurrentResidue << endl;
		}
	}

	unsigned int i;

	// bubble sort the list by mass
	unsigned int n = vsResidues.size();
	unsigned int pass;
	double dCurrentMass;
	for (pass = 1; pass < n; pass++)
	{ // count how many times
		// This next loop becomes shorter and shorter
		for (i = 0; i < n - pass; i++)
		{
			if (vdMostAbundantMasses[i] > vdMostAbundantMasses[i + 1])
			{
				// exchange
				dCurrentMass = vdMostAbundantMasses[i];
				sCurrentResidue = vsResidues[i];

				vdMostAbundantMasses[i] = vdMostAbundantMasses[i + 1];
				vsResidues[i] = vsResidues[i + 1];

				vdMostAbundantMasses[i + 1] = dCurrentMass;
				vsResidues[i + 1] = sCurrentResidue;
			}
		}
	}

	return true;
}

bool Isotopologue::computeIsotopicDistribution(string sSequence, IsotopeDistribution &myIsotopeDistribution)
{
	IsotopeDistribution sumDistribution;
	IsotopeDistribution currentDistribution;
	map<string, IsotopeDistribution>::iterator ResidueIter;

	ResidueIter = vResidueIsotopicDistribution.find("Nterm");
	if (ResidueIter != vResidueIsotopicDistribution.end())
	{
		currentDistribution = ResidueIter->second;
		sumDistribution = currentDistribution;
	}
	else
	{
		cerr << "ERROR: can't find the N-terminus" << endl;
		return false;
	}

	ResidueIter = vResidueIsotopicDistribution.find("Cterm");
	if (ResidueIter != vResidueIsotopicDistribution.end())
	{
		currentDistribution = ResidueIter->second;
		sumDistribution = sum(currentDistribution, sumDistribution);
	}
	else
	{
		cerr << "ERROR: can't find the C-terminus" << endl;
		return false;
	}

	// add up all residues's isotopic distribution
	for (unsigned int j = 0; j < sSequence.length(); j++)
	{
		string currentResidue = sSequence.substr(j, 1);
		ResidueIter = vResidueIsotopicDistribution.find(currentResidue);
		if (ResidueIter != vResidueIsotopicDistribution.end())
		{
			currentDistribution = ResidueIter->second;
			sumDistribution = sum(currentDistribution, sumDistribution);
		}
		else
		{
			cerr << "ERROR: can't find the residue " << currentResidue << endl;
			return false;
		}
	}

	myIsotopeDistribution = sumDistribution;

	return true;
}

bool Isotopologue::computeProductIon(string sSequence,
									 vector<vector<double>> &vvdYionMass,
									 vector<vector<double>> &vvdYionProb,
									 vector<vector<double>> &vvdBionMass,
									 vector<vector<double>> &vvdBionProb)
{
	// get the mass for a proton
	double dProtonMass = ProNovoConfig::getProtonMass();
	unsigned int i = 0;

	//	if(!isalpha(sSequence[0]))
	//	{
	//		cout << "ERROR: First character in a peptide sequence can't be a PTM." << endl;
	//		return false;
	//	}

	int iPeptideLength = 0;
	for (i = 0; i < sSequence.length(); ++i)
	{
		if (isalpha(sSequence[i]))
		{
			iPeptideLength = iPeptideLength + 1;
		}
	}

	if (iPeptideLength < ProNovoConfig::getMinPeptideLength())
	{
		cerr << "ERROR: Peptide sequence is too short " << sSequence << endl;
		return false;
	}

	//	cout << "sSequence = " << sSequence << endl;

	vvdYionMass.clear();
	vvdYionProb.clear();
	vvdBionMass.clear();
	vvdBionProb.clear();

	vvdYionMass.reserve(iPeptideLength);
	vvdYionProb.reserve(iPeptideLength);
	vvdBionMass.reserve(iPeptideLength);
	vvdBionProb.reserve(iPeptideLength);

	map<string, IsotopeDistribution>::iterator ResidueIter;

	vector<IsotopeDistribution> vResidueDistribution;
	vResidueDistribution.reserve(iPeptideLength + 2);
	IsotopeDistribution sumDistribution;
	IsotopeDistribution currentDistribution;

	string currentResidue;
	string currentPTM;

	if (sSequence[0] != '[')
	{
		cerr << "ERROR: First character in a peptide sequence must be [." << endl;
		return false;
	}

	unsigned int iStartResidueIndex = 1;
	ResidueIter = vResidueIsotopicDistribution.find("Nterm");
	if (ResidueIter != vResidueIsotopicDistribution.end())
	{
		currentDistribution = ResidueIter->second;

		if (!isalpha(sSequence[1]))
		{
			iStartResidueIndex = 2;
			currentPTM = sSequence[1];
			ResidueIter = vResidueIsotopicDistribution.find(currentPTM);
			if (ResidueIter == vResidueIsotopicDistribution.end())
			{
				cerr << "ERROR: cannot find this PTM in the config file " << currentPTM << endl;
				return false;
			}

			currentDistribution = sum(currentDistribution, ResidueIter->second);
		}

		vResidueDistribution.push_back(currentDistribution);
	}
	else
	{
		cerr << "ERROR: can't find the N-terminus" << endl;
		return false;
	}

	for (i = iStartResidueIndex; i < sSequence.length(); i++)
	{
		if (sSequence[i] == ']')
		{
			break;
		}

		if (!isalpha(sSequence[i]))
		{
			cerr << "ERROR: One residue can only have one PTM (Up to only one symbol after an amino acid)" << endl;
			return false;
		}
		currentResidue = sSequence.substr(i, 1);
		ResidueIter = vResidueIsotopicDistribution.find(currentResidue);
		if (ResidueIter == vResidueIsotopicDistribution.end())
		{
			cerr << "ERROR: cannot find this residue in the config file. " << currentResidue << endl;
			return false;
		}
		currentDistribution = ResidueIter->second;
		if (i + 1 < sSequence.length())
		{
			if (!isalpha(sSequence[i + 1]) && sSequence[i + 1] != ']') // this residue is modified
			{
				currentPTM = sSequence[i + 1];
				ResidueIter = vResidueIsotopicDistribution.find(currentPTM);
				if (ResidueIter == vResidueIsotopicDistribution.end())
				{
					cerr << "ERROR: cannot find this PTM in the config file " << currentPTM << endl;
					return false;
				}

				// this PTM can substract mass from the residue
				// the currentDistribution could even be negative.
				currentDistribution = sum(currentDistribution, ResidueIter->second);

				i = i + 1; // increment index to the next residue
			}
		}

		vResidueDistribution.push_back(currentDistribution);
	}

	ResidueIter = vResidueIsotopicDistribution.find("Cterm");
	if (ResidueIter != vResidueIsotopicDistribution.end())
	{
		currentDistribution = ResidueIter->second;

		if (i + 1 < sSequence.length())
		{
			if (!isalpha(sSequence[i + 1]))
			{
				currentPTM = sSequence[i + 1];
				ResidueIter = vResidueIsotopicDistribution.find(currentPTM);
				if (ResidueIter == vResidueIsotopicDistribution.end())
				{
					cerr << "ERROR: cannot find this PTM in the config file " << currentPTM << endl;
					return false;
				}

				currentDistribution = sum(currentDistribution, ResidueIter->second);
			}
		}

		vResidueDistribution.push_back(currentDistribution);
	}
	else
	{
		cerr << "ERROR: can't find the C-terminus" << endl;
		return false;
	}

	// compute B-ion series

	vector<IsotopeDistribution> vBionDistribution;
	vBionDistribution.reserve(iPeptideLength);

	// start with the N-terminus distibution, which should the first element
	currentDistribution = vResidueDistribution.front();
	sumDistribution = currentDistribution;

	int j;
	for (j = 1; j < iPeptideLength; j++)
	{
		currentDistribution = vResidueDistribution[j];
		sumDistribution = sum(currentDistribution, sumDistribution);
		vBionDistribution.push_back(sumDistribution);
	}

#pragma omp simd
	for (i = 0; i < vBionDistribution.size(); i++)
	{
		vvdBionMass.push_back(move(vBionDistribution[i].vMass));
		vvdBionProb.push_back(move(vBionDistribution[i].vProb));
		//	cout << "b " << i << endl;
		//	vBionDistribution[i].print();
	}

	// compute Y-ion series
	vector<IsotopeDistribution> vYionDistribution;
	vYionDistribution.reserve(iPeptideLength);

	// start with the C-terminus distibution, which should the last element
	currentDistribution = vResidueDistribution.back();
	sumDistribution = currentDistribution;

	for (j = iPeptideLength; j > 1; j--)
	{
		currentDistribution = vResidueDistribution[j];
		sumDistribution = sum(currentDistribution, sumDistribution);
		vYionDistribution.push_back((sumDistribution));
	}

#pragma omp simd
	for (i = 0; i < vYionDistribution.size(); i++)
	{
		vvdYionMass.push_back(move(vYionDistribution[i].vMass));
		vvdYionProb.push_back(move(vYionDistribution[i].vProb));
		//	cout << "y " << i << endl;
		//	vYionDistribution[i].print();
	}

	// change the masses to correct for the proton transfer during peptide bond cleavage
	unsigned int n;
	unsigned int m;
#pragma omp simd
	for (n = 0; n < vvdYionMass.size(); ++n)
	{
		for (m = 0; m < vvdYionMass[n].size(); ++m)
		{
			vvdYionMass[n][m] += dProtonMass;
		}
	}
#pragma omp simd
	for (n = 0; n < vvdBionMass.size(); ++n)
	{
		for (m = 0; m < vvdBionMass[n].size(); ++m)
		{
			vvdBionMass[n][m] -= dProtonMass;
		}
	}

	//	cout << "Y		B" << endl;
	//	for (n = 0; n < vvdBionMass.size(); ++n)
	//	{
	//		cout << vvdYionMass[n][0] << "\t\t" << vvdBionMass[n][0] << endl;
	//	}

	return true;
}

bool Isotopologue::computeIsotopicDistribution(vector<int> AtomicComposition, IsotopeDistribution &myIsotopeDistribution)
{
	IsotopeDistribution sumDistribution;
	IsotopeDistribution currentAtomDistribution;
	currentAtomDistribution = multiply(vAtomIsotopicDistribution[0], AtomicComposition[0]);
	sumDistribution = currentAtomDistribution;

	for (unsigned int i = 1; i < AtomNumber; i++)
	{
		currentAtomDistribution = multiply(vAtomIsotopicDistribution[i], AtomicComposition[i]);
		sumDistribution = sum(currentAtomDistribution, sumDistribution);
	}
	myIsotopeDistribution = sumDistribution;

	return true;
}

bool Isotopologue::computeAtomicComposition(string sSequence, vector<int> &myAtomicComposition)
{
	vector<int> AtomicComposition;
	vector<int> CurrentComposition;
	unsigned int i;
	map<string, vector<int>>::iterator ResidueIter;

	for (i = 0; i < AtomNumber; i++)
		AtomicComposition.push_back(0);

	ResidueIter = mResidueAtomicComposition.find("Nterm");
	if (ResidueIter != mResidueAtomicComposition.end())
	{
		CurrentComposition = ResidueIter->second;
		for (i = 0; i < AtomNumber; i++)
			AtomicComposition[i] = CurrentComposition[i];
	}
	else
	{
		cerr << "ERROR: can't find the atomic composition for the N-terminus" << endl;
		return false;
	}

	ResidueIter = mResidueAtomicComposition.find("Cterm");
	if (ResidueIter != mResidueAtomicComposition.end())
	{
		CurrentComposition = ResidueIter->second;
		for (i = 0; i < AtomNumber; i++)
			AtomicComposition[i] += CurrentComposition[i];
	}
	else
	{
		cerr << "ERROR: can't find the atomic composition for the C-terminus" << endl;
		return false;
	}

	for (unsigned int j = 0; j < sSequence.length(); j++)
	{
		string currentResidue = sSequence.substr(j, 1);
		ResidueIter = mResidueAtomicComposition.find(currentResidue);
		if (ResidueIter != mResidueAtomicComposition.end())
		{
			CurrentComposition = ResidueIter->second;
			for (i = 0; i < AtomNumber; i++)
				AtomicComposition[i] += CurrentComposition[i];
		}
		else
		{
			cerr << "ERROR: can't find the atomic composition for residue/PTM: " << currentResidue << endl;
			return false;
		}
	}

	myAtomicComposition = AtomicComposition;
	return true;
}

IsotopeDistribution Isotopologue::sum(const IsotopeDistribution &distribution0, const IsotopeDistribution &distribution1)
{
	double ProbabilityCutoff_local = 0.000001;

	IsotopeDistribution sumDistribution;
	double currentMass;
	double currentProb;
	int iSizeDistribution0 = distribution0.vMass.size();
	int iSizeDistribution1 = distribution1.vMass.size();
	int iCount = 0;
	double dSum = 0;
	size_t newSize = iSizeDistribution0 + iSizeDistribution1 - 1;
	sumDistribution.vMass.reserve(newSize);
	sumDistribution.vProb.reserve(newSize);
#pragma omp simd
	for (int k = 0; k < newSize; k++)
	{
		double sumweight = 0, summass = 0;
		int start = k < (iSizeDistribution1 - 1) ? 0 : k - iSizeDistribution1 + 1; // max(0, k-f_n+1)
		int end = k < (iSizeDistribution0 - 1) ? k : iSizeDistribution0 - 1;	   // min(g_n - 1, k)
		iCount = 0;
		dSum = 0;
		for (int i = start; i <= end; i++)
		{
			double weight = distribution0.vProb[i] * distribution1.vProb[k - i];
			double mass = distribution0.vMass[i] + distribution1.vMass[k - i];
			sumweight += weight;
			summass += weight * mass;
			iCount += 1;
			dSum += mass;
		}
		if (sumweight == 0)
		{
			currentMass = dSum / ((double)iCount);
		}
		else
		{
			currentMass = summass / sumweight;
		}
		currentProb = sumweight;
		sumDistribution.vMass.push_back(currentMass);
		sumDistribution.vProb.push_back(currentProb);
	}

	//prune small probabilities
	vector<double>::iterator iteProb = sumDistribution.vProb.begin();
	vector<double>::iterator iteMass = sumDistribution.vMass.begin();
	while (iteProb != sumDistribution.vProb.end())
	{
		if ((*iteProb) > ProbabilityCutoff_local)
		{
			break;
		}
		iteProb++;
		iteMass++;
	}
	if (iteProb != sumDistribution.vProb.begin())
	{
		sumDistribution.vProb.erase(sumDistribution.vProb.begin(), iteProb);
		sumDistribution.vMass.erase(sumDistribution.vMass.begin(), iteMass);
	}

	iteProb = sumDistribution.vProb.end() - 1;
	iteMass = sumDistribution.vMass.end() - 1;
	while (iteProb != sumDistribution.vProb.begin())
	{
		if ((*iteProb) > ProbabilityCutoff_local)
		{
			break;
		}
		iteProb--;
		iteMass--;
	}
	if (iteProb != sumDistribution.vProb.end() - 1)
	{
		sumDistribution.vProb.erase(iteProb + 1, sumDistribution.vProb.end());
		sumDistribution.vMass.erase(iteMass + 1, sumDistribution.vMass.end());
	}

	// normalize the probability space to 1
	int iSizeSumDistribution;
	int i;
	double sumProb = 0;
	iSizeSumDistribution = sumDistribution.vMass.size();
	for (i = 0; i < iSizeSumDistribution; ++i)
	{
		sumProb += sumDistribution.vProb[i];
	}

	if (sumProb <= 0)
	{
		cerr << "Error: Sum of distribution is zero" << endl;
		exit(1);
		return sumDistribution;
	}

	for (i = 0; i < iSizeSumDistribution; ++i)
	{
		sumDistribution.vProb[i] = sumDistribution.vProb[i] / sumProb;
	}

	return sumDistribution;
}

IsotopeDistribution Isotopologue::multiply(const IsotopeDistribution &distribution0, int count)
{
	if (count == 1)
		return distribution0;
	IsotopeDistribution productDistribution;
	productDistribution.vMass.push_back(0.0);
	productDistribution.vProb.push_back(1.0);

	if (count < 0)
	{
		IsotopeDistribution negativeDistribution = distribution0;
#pragma omp simd
		for (unsigned int n = 0; n < negativeDistribution.vMass.size(); n++)
		{
			negativeDistribution.vMass[n] = -negativeDistribution.vMass[n];
		}

		for (int i = 0; i < abs(count); ++i)
		{
			productDistribution = sum(productDistribution, negativeDistribution);
		}
	}
	else
	{
		for (int i = 0; i < count; ++i)
		{
			productDistribution = sum(productDistribution, distribution0);
		}
	}
	return productDistribution;
}

void Isotopologue::shiftMass(IsotopeDistribution &distribution0, double dMass)
{
	for (unsigned int i = 0; i < distribution0.vMass.size(); ++i)
	{
		distribution0.vMass[i] = distribution0.vMass[i] + dMass;
	}
}
