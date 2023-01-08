/*
 * MVH.cpp
 *
 *  Created on: May 23, 2016
 *      Author: xgo
 */

#include "MVH.h"

bool MVH::bUseSmartPlusThreeModel = true;
lnFactorialTable * MVH::lnTable = NULL;
bitset<FragmentTypes_Size> MVH::fragmentTypes(string("0010010"));
double MVH::ProbabilityCutOff = 0.25;

MVH::MVH() {
	// TODO Auto-generated constructor stub
}

MVH::~MVH() {
	// TODO Auto-generated destructor stub
}

double round(double f, int precision) {
	if (f == 0.0f)
		return +0.0f;

	double multiplier = pow(10.0, (double) precision); // moves f over <precision> decimal places
	f *= multiplier;
	f = floor(f + 0.5f);
	return f / multiplier;
}

/**
 * variables needed from MS2Scan
 * vdMZ
 * vdIntensity
 * dParentNeutralMass
 * dParentMZ
 */
bool MVH::Preprocess(MS2Scan * Spectrum, multimap<double, double> * IntenSortedPeakPreData) {

	vector<double> *vdMZ = &(Spectrum->vdMZ);
	vector<double> *vdIntensity = &(Spectrum->vdIntensity);
	if (((int) vdMZ->size()) < ProNovoConfig::minIntensityClassCount) {
		Spectrum->bSkip = true;
		return true;
	}
	// Spectrum->mzLowerBound = vdMZ->front();
	// Spectrum->mzUpperBound = vdMZ->back();
	Spectrum->mzLowerBound = ProNovoConfig::minObservedMz;
	Spectrum->mzUpperBound = ProNovoConfig::maxObservedMz;
	double totalPeakSpace = Spectrum->mzUpperBound - Spectrum->mzLowerBound;
	double totalIonCurrent = 0;
	if (IntenSortedPeakPreData == NULL) {
		cerr << "Error 61" << endl;
		exit(1);
		return true;
	}
	IntenSortedPeakPreData->clear();
	multimap<double, double>::iterator ite;
	for (int i = 0; i < (int) vdMZ->size(); i++) {
		if (vdMZ->at(i) <= Spectrum->dParentNeutralMass) {
			ite = IntenSortedPeakPreData->insert(make_pair(vdIntensity->at(i), vdMZ->at(i)));
			totalIonCurrent += ite->first;
		}
	}
	// FilterByTIC
	// Filters out the peaks with the lowest intensities until only <ticCutoffPercentage> of the total ion current remains
	if (IntenSortedPeakPreData->empty()) {
		Spectrum->bSkip = true;
		cerr << "error MVH::Preprocess 1" << endl;
		return true;
	}
	multimap<double, double>::reverse_iterator rite;
	double relativeIntensity = 0.0;
	for (rite = IntenSortedPeakPreData->rbegin(); relativeIntensity < ProNovoConfig::ticCutoffPercentage && rite != IntenSortedPeakPreData->rend(); ++rite) {
		relativeIntensity += (rite->first / totalIonCurrent);
	}
	if (rite == IntenSortedPeakPreData->rend()) {
		--rite;
	}
	ite = IntenSortedPeakPreData->lower_bound(rite->first);
	if (ite != IntenSortedPeakPreData->begin()) {
		IntenSortedPeakPreData->erase(IntenSortedPeakPreData->begin(), ite);
	}

	//FilterByPeakCount
	if (!IntenSortedPeakPreData->empty()) {
		if (((int) IntenSortedPeakPreData->size()) > ProNovoConfig::MaxPeakCount) {
			int peakCount = IntenSortedPeakPreData->size();
			ite = IntenSortedPeakPreData->begin();
			while (peakCount > ProNovoConfig::MaxPeakCount) {
				++ite;
				peakCount--;
			}
			IntenSortedPeakPreData->erase(IntenSortedPeakPreData->begin(), ite);
		}
	}
	//Water Loss
	double dMinOneWater = Spectrum->dParentMZ - WATER_MONO / Spectrum->iParentChargeState - ProNovoConfig::getMassAccuracyParentIon();
	double dMinTwoWater = Spectrum->dParentMZ - 2 * WATER_MONO / Spectrum->iParentChargeState - ProNovoConfig::getMassAccuracyParentIon();
	double dMaxOneWater = Spectrum->dParentMZ - WATER_MONO / Spectrum->iParentChargeState + ProNovoConfig::getMassAccuracyParentIon();
	double dMaxTwoWater = Spectrum->dParentMZ - 2 * WATER_MONO / Spectrum->iParentChargeState + ProNovoConfig::getMassAccuracyParentIon();
	double maxIntenOneWater = 0;
	double dTargeMassOneWater = 0;
	double maxIntenTwoWater = 0;
	double dTargeMassTwoWater = 0;
	for (ite = IntenSortedPeakPreData->begin(); ite != IntenSortedPeakPreData->end(); ++ite) {
		if (ite->second > dMinOneWater && ite->second < dMaxOneWater) {
			if (maxIntenOneWater < ite->first) {
				maxIntenOneWater = ite->first;
				dTargeMassOneWater = ite->second;
			}
		} else if (ite->second > dMinTwoWater && ite->second < dMaxTwoWater) {
			if (maxIntenTwoWater < ite->first) {
				maxIntenTwoWater = ite->first;
				dTargeMassTwoWater = ite->second;
			}
		}
	}
	if (dTargeMassOneWater > 0) {
		pair<multimap<double, double>::iterator, multimap<double, double>::iterator> iterpair = IntenSortedPeakPreData->equal_range(maxIntenOneWater);
		for (ite = iterpair.first; ite != iterpair.second; ++ite) {
			if (ite->second == dTargeMassOneWater) {
				IntenSortedPeakPreData->erase(ite);
				break;
			}
		}
	}
	if (dTargeMassTwoWater > 0) {
		pair<multimap<double, double>::iterator, multimap<double, double>::iterator> iterpair = IntenSortedPeakPreData->equal_range(maxIntenTwoWater);
		for (ite = iterpair.first; ite != iterpair.second; ite++) {
			if (ite->second == dTargeMassTwoWater) {
				IntenSortedPeakPreData->erase(ite);
				break;
			}
		}
	}
	//ClassifyPeakIntensities for mvh
	map<double, char> * peakData = Spectrum->peakData;
	rite = IntenSortedPeakPreData->rbegin();
	for (int i = 0; i < ProNovoConfig::NumIntensityClasses; ++i) {
		int numFragments = (int) round(
				(double) (pow((double) ProNovoConfig::ClassSizeMultiplier, i) * IntenSortedPeakPreData->size())
						/ (double) ProNovoConfig::minIntensityClassCount, 0);
		for (int j = 0; rite != IntenSortedPeakPreData->rend() && j < numFragments; ++j, ++rite) {
			(*peakData)[rite->second] = i + 1;
		}
	}
	vector<int> * intenClassCounts = Spectrum->intenClassCounts;
	intenClassCounts->clear();
	intenClassCounts->resize(ProNovoConfig::NumIntensityClasses + 1, 0);
	map<double, char>::iterator itr;
	for (itr = peakData->begin(); itr != peakData->end(); itr++) {
		if (itr->second > 0) {
			++(intenClassCounts->at(itr->second - 1));
		}
	}
	double fragMassError = ProNovoConfig::getMassAccuracyFragmentIon();
	int totalPeakBins = (int) round(totalPeakSpace / (fragMassError * 2.0f), 0);
	int peakCount = (int) peakData->size();
	intenClassCounts->at(ProNovoConfig::NumIntensityClasses) = totalPeakBins - peakCount;
	Spectrum->totalPeakBins = totalPeakBins;
	Spectrum->bSkip = false;
	IntenSortedPeakPreData->clear();
	return true;
}

bool MVH::CalculateSequenceIons(string & sSequence, int maxIonCharge, bool useSmartPlusThreeModel, vector<double>* sequenceIonMasses,
		vector<double> * _pdAAforward, vector<double> * _pdAAreverse, vector<char> * seq) {
	sequenceIonMasses->clear();
	_pdAAforward->clear();
	_pdAAreverse->clear();
	seq->clear();
	int i = 0, j = 0, k = 0, iPeptideLength, iPos;
	k = ((int) sSequence.length()) - 1;
	char currentPTM = 0;
	iPeptideLength = 0;
	iPos = 0;
	double iterResidueMonoMass;
	for (i = 0; i <= k; ++i) {
		if (isalpha(sSequence.at(i))) {
			iPeptideLength = iPeptideLength + 1;
			seq->push_back(sSequence.at(i));
		}
	}
	int iLenMinus1 = iPeptideLength - 1;
	if (iPeptideLength < ProNovoConfig::getMinPeptideLength()) {
		cerr << "ERROR MVH 1: Peptide sequence is too short " << sSequence << endl;
		exit(1);
		return false;
	}
	j = 0;
	if (sSequence.at(j) != '[') {
		cerr << "ERROR: First character in a peptide sequence must be [." << sSequence << endl;
		return false;
	}
	double dBion = 0;
	double dYion = ProNovoConfig::precalcMasses.dCtermOH2;
	j++;
	if (!isalpha(sSequence.at(j))) {
		currentPTM = sSequence.at(j);
		iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
		if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
			cerr << "ERROR: cannot find this PTM in the config file " << currentPTM << endl;
			return false;
		}
		dBion += iterResidueMonoMass;
		j++;
	}
	if (sSequence.at(k) == ']') {
		k--;
	} else {
		currentPTM = sSequence.at(k);
		iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
		if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
			cerr << "ERROR: cannot find this PTM in the config file " << currentPTM << endl;
			return false;
		}
		dYion += iterResidueMonoMass;
		k--;
		if (sSequence.at(k) != ']') {
			cerr << "ERROR: (second) Last character in a peptide sequence must be ]." << endl;
			return false;
		}
		k--;
	}
	//
	for (i = 0; i <= iLenMinus1; i++) {
		//First forward
		if (!isalpha(sSequence.at(j))) {
			cerr << "ERROR: One residue can only have one PTM (Up to only one symbol after an amino acid)" << endl;
			exit(1);
			return false;
		}
		currentPTM = sSequence.at(j);
		iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
		if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
			cerr << "ERROR: One residue can only have one PTM (Up to only one symbol after an amino acid)" << endl;
			exit(1);
			return false;
		}
		dBion += iterResidueMonoMass;
		j++;
		if (j < (int) sSequence.length() && !isalpha(sSequence.at(j)) && sSequence.at(j) != ']') {
			currentPTM = sSequence.at(j);
			iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
			if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
				cerr << "ERROR: cannot find this PTM in the config file " << currentPTM << endl;
				return false;
			}
			dBion += iterResidueMonoMass;
			j++;
		}
		_pdAAforward->push_back(dBion);
		// now reverse
		if (!isalpha(sSequence.at(k))) {
			currentPTM = sSequence.at(k);
			iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
			if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
				cerr << "ERROR: cannot find this PTM in the config file " << currentPTM << endl;
				exit(1);
				return false;
			}
			dYion += iterResidueMonoMass;
			k--;
		}
		if (!isalpha(sSequence.at(k))) {
			cerr << "ERROR: One residue can only have one PTM (Up to only one symbol after an amino acid)" << endl;
			exit(1);
			return false;
		}
		currentPTM = sSequence.at(k);
		iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
		if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
			cerr << "ERROR: cannot find this PTM in the config file" << currentPTM << endl;
			exit(1);
			return false;
		}
		dYion += iterResidueMonoMass;
		k--;
		_pdAAreverse->push_back(dYion);
		iPos++;
	}

	// calculate y ion MZs
	if (maxIonCharge > 2) {
		if (useSmartPlusThreeModel) {
			int totalStrongBasicCount = 0, totalWeakBasicCount = 0;
			for (i = 0; i < iPeptideLength; i++) {
				if (seq->at(i) == 'R' || seq->at(i) == 'K' || seq->at(i) == 'H') {
					++totalStrongBasicCount;
				} else if (seq->at(i) == 'Q' || seq->at(i) == 'N') {
					++totalWeakBasicCount;
				}
			}
			int totalBasicity = totalStrongBasicCount * 4 + totalWeakBasicCount * 2 + iPeptideLength - 2;
			map<double, int> basicityThresholds;
			basicityThresholds.clear();
			basicityThresholds[0.0] = 1;
			for (int z = 1; z < maxIonCharge - 1; ++z) {
				basicityThresholds[(double) z / (double) (maxIonCharge - 1)] = z + 1;
			}
			for (int c = 0; c <= iPeptideLength; ++c) {
				int bStrongBasicCount = 0, bWeakBasicCount = 0;
				for (i = 0; i < c; ++i) {
					if (seq->at(i) == 'R' || seq->at(i) == 'K' || seq->at(i) == 'H') {
						++bStrongBasicCount;
					} else if (seq->at(i) == 'Q' || seq->at(i) == 'N') {
						++bWeakBasicCount;
					}
				}
				int bScore = bStrongBasicCount * 4 + bWeakBasicCount * 2 + c;
				double basicityRatio = (double) bScore / (double) totalBasicity;
				map<double, int>::iterator itr = basicityThresholds.upper_bound(basicityRatio);
				if (itr == basicityThresholds.begin()) {
					cerr << "error 83" << endl;
				}
				--itr;
				int bZ = itr->second;
				int yZ = maxIonCharge - bZ;
				if (c > 0) {
					if (fragmentTypes[FragmentType_B]) {
						sequenceIonMasses->push_back((_pdAAforward->at(c - 1) + (Proton * bZ)) / bZ);
					}
				}
				if (c < iPeptideLength) {
					if (fragmentTypes[FragmentType_Y]) {
						sequenceIonMasses->push_back((_pdAAreverse->at(iLenMinus1 - c) + (Proton * yZ)) / yZ);
					}
				}
			}
		} else {	// no Smart Plus Three Model
			for (int z = 1; z < maxIonCharge; ++z) {
				for (int c = 0; c <= iLenMinus1; ++c) {
					if (c > 0) {
						if (fragmentTypes[FragmentType_B]) {
							sequenceIonMasses->push_back((_pdAAforward->at(c - 1) + (Proton * z)) / z);
						}
					}
					if (fragmentTypes[FragmentType_Y]) {
						sequenceIonMasses->push_back((_pdAAreverse->at(c) + (Proton * z)) / z);
					}
				}
			}
		}
	} else {
		for (int c = 0; c < iLenMinus1; ++c) {
			if (fragmentTypes[FragmentType_B]) {
				sequenceIonMasses->push_back((_pdAAforward->at(c) + Proton));
			}
			if (fragmentTypes[FragmentType_Y]) {
				sequenceIonMasses->push_back((_pdAAreverse->at(c) + Proton));
			}
		}
		if (fragmentTypes[FragmentType_Y]) {
			sequenceIonMasses->push_back((_pdAAreverse->at(iLenMinus1) + Proton));
		}
	}
	return true;
}

bool MVH::CalculateSequenceIonsSIP(string & sSequence, int maxIonCharge, bool useSmartPlusThreeModel, vector<double>* sequenceIonMasses,
		vector<vector<double> > & vvdYionMass, vector<vector<double> > & vvdYionProb, vector<vector<double> > & vvdBionMass,
		vector<vector<double> > & vvdBionProb, vector<char> * seq) {
	sequenceIonMasses->clear();
	seq->clear();
	int k = 0, iPeptideLength = 0;
	k = ((int) sSequence.length()) - 1;
	for (int i = 0; i <= k; ++i) {
		if (isalpha(sSequence.at(i))) {
			iPeptideLength = iPeptideLength + 1;
			seq->push_back(sSequence.at(i));
		}
	}

	int iIonMassSize = 0, iIsotopicDistSize = 0;
	double dFragmentIonMassZ = 0;
	int iLenMinus1 = iPeptideLength - 1;
	if (((int)vvdYionMass.size()) != iLenMinus1 || ((int)vvdBionMass.size()) != iLenMinus1) {
		cout << "check 111" << endl;
	}
	// calculate y ion MZs
	if (maxIonCharge > 2) {
		if (useSmartPlusThreeModel) {
			int totalStrongBasicCount = 0, totalWeakBasicCount = 0;
			for (int i = 0; i < iPeptideLength; i++) {
				if (seq->at(i) == 'R' || seq->at(i) == 'K' || seq->at(i) == 'H') {
					++totalStrongBasicCount;
				} else if (seq->at(i) == 'Q' || seq->at(i) == 'N') {
					++totalWeakBasicCount;
				}
			}
			int totalBasicity = totalStrongBasicCount * 4 + totalWeakBasicCount * 2 + iPeptideLength - 2;
			map<double, int> basicityThresholds;
			basicityThresholds.clear();
			basicityThresholds[0.0] = 1;
			for (int z = 1; z < maxIonCharge - 1; ++z) {
				basicityThresholds[(double) z / (double) (maxIonCharge - 1)] = z + 1;
			}
			for (int c = 1; c < iPeptideLength; ++c) {
				int bStrongBasicCount = 0, bWeakBasicCount = 0;
				for (int i = 0; i < c; ++i) {
					if (seq->at(i) == 'R' || seq->at(i) == 'K' || seq->at(i) == 'H') {
						++bStrongBasicCount;
					} else if (seq->at(i) == 'Q' || seq->at(i) == 'N') {
						++bWeakBasicCount;
					}
				}
				int bScore = bStrongBasicCount * 4 + bWeakBasicCount * 2 + c;
				double basicityRatio = (double) bScore / (double) totalBasicity;
				map<double, int>::iterator itr = basicityThresholds.upper_bound(basicityRatio);
				if (itr == basicityThresholds.begin()) {
					cerr << "error 83" << endl;
				}
				--itr;
				int bZ = itr->second;
				int yZ = maxIonCharge - bZ;
				if (fragmentTypes[FragmentType_B]) {
					iIsotopicDistSize = vvdBionMass.at(c - 1).size();
					for (int j = 0; j < iIsotopicDistSize; ++j) {
						if (vvdBionProb.at(c - 1).at(j) < ProbabilityCutOff) {
							continue;
						}
						dFragmentIonMassZ = (vvdBionMass.at(c - 1).at(j) / bZ + Proton);
						sequenceIonMasses->push_back(dFragmentIonMassZ);
					}
				}

				if (fragmentTypes[FragmentType_Y]) {
					iIsotopicDistSize = vvdYionMass.at(iLenMinus1 - c).size();
					for (int j = 0; j < iIsotopicDistSize; ++j) {
						if (vvdYionProb.at(iLenMinus1 - c).at(j) < ProbabilityCutOff) {
							continue;
						}
						dFragmentIonMassZ = (vvdYionMass.at(iLenMinus1 - c).at(j) / yZ + Proton);
						sequenceIonMasses->push_back(dFragmentIonMassZ);
					}
				}

			}
		} else {	// no Smart Plus Three Model
			for (int z = 1; z < maxIonCharge; ++z) {
				// b-ion
				if (fragmentTypes[FragmentType_B]) {
					iIonMassSize = vvdBionMass.size();
					for (int i = 0; i < iIonMassSize; ++i) {
						iIsotopicDistSize = vvdBionMass.at(i).size();
						for (int j = 0; j < iIsotopicDistSize; ++j) {
							if (vvdBionProb.at(i).at(j) < ProbabilityCutOff) {
								continue;
							}
							dFragmentIonMassZ = (vvdBionMass.at(i).at(j) / z + Proton);
							sequenceIonMasses->push_back(dFragmentIonMassZ);
						}
					}
				}
				// y-ion
				if (fragmentTypes[FragmentType_Y]) {
					iIonMassSize = vvdYionMass.size();
					for (int i = 0; i < iIonMassSize; ++i) {
						iIsotopicDistSize = vvdYionMass.at(i).size();
						for (int j = 0; j < iIsotopicDistSize; ++j) {
							if (vvdYionProb.at(i).at(j) < ProbabilityCutOff) {
								continue;
							}
							dFragmentIonMassZ = (vvdYionMass.at(i).at(j) / z + Proton);
							sequenceIonMasses->push_back(dFragmentIonMassZ);
						}
					}
				}
			}
		}
	} else {
		// b-ion
		if (fragmentTypes[FragmentType_B]) {
			iIonMassSize = vvdBionMass.size();
			for (int i = 0; i < iIonMassSize; ++i) {
				iIsotopicDistSize = vvdBionMass.at(i).size();
				for (int j = 0; j < iIsotopicDistSize; ++j) {
					if (vvdBionProb.at(i).at(j) < ProbabilityCutOff) {
						continue;
					}
					dFragmentIonMassZ = (vvdBionMass.at(i).at(j) + Proton);
					sequenceIonMasses->push_back(dFragmentIonMassZ);
				}
			}
		}
		// y-ion
		if (fragmentTypes[FragmentType_Y]) {
			iIonMassSize = vvdYionMass.size();
			for (int i = 0; i < iIonMassSize; ++i) {
				iIsotopicDistSize = vvdYionMass.at(i).size();
				for (int j = 0; j < iIsotopicDistSize; ++j) {
					if (vvdYionProb.at(i).at(j) < ProbabilityCutOff) {
						continue;
					}
					dFragmentIonMassZ = (vvdYionMass.at(i).at(j) + Proton);
					sequenceIonMasses->push_back(dFragmentIonMassZ);
				}
			}
		}
	}
	return true;
}

bool MVH::destroyLnTable() {
	if (MVH::lnTable != NULL) {
		delete MVH::lnTable;
		MVH::lnTable = NULL;
	}
	return true;
}

multimap<double, char>::iterator MVH::findNear(map<double, char> * peakData, double mz, double tolerance) {
	multimap<double, char>::iterator cur, min, max, best;
	min = peakData->lower_bound(mz - tolerance);
	max = peakData->upper_bound(mz + tolerance);
	if (min == max) {
		return peakData->end();
	}
	best = min;
	double minDiff = fabs(mz - best->first);
	for (cur = min; cur != max; ++cur) {
		double curDiff = fabs(mz - cur->first);
		if (curDiff < minDiff) {
			minDiff = curDiff;
			best = cur;
		}
	}
	return best;
}

bool MVH::initialLnTable(int maxPeakBins) {
	if (MVH::lnTable != NULL) {
		delete MVH::lnTable;
		MVH::lnTable = NULL;
	}
	MVH::lnTable = new lnFactorialTable();
	MVH::lnTable->resize(maxPeakBins);
	return true;
}

double MVH::lnCombin(int n, int k) {
	if (n < 0 || k < 0 || n < k)
		return -1;

	try {
		return (*lnTable)[n] - (*lnTable)[n - k] - (*lnTable)[k];
	} catch (std::exception& e) {
		cerr << "lnCombin(): caught exception with n=" << n << " and k=" << k << endl;
		throw e;
	}
}

bool MVH::ScoreSequenceVsSpectrum(string & currentPeptide, MS2Scan * Spectrum, vector<double>* seqIons, vector<double> * _pdAAforward,
		vector<double> * _pdAAreverse, double & dMvh, vector<char> * seq) {
	if (!CalculateSequenceIons(currentPeptide, Spectrum->iParentChargeState, bUseSmartPlusThreeModel, seqIons, _pdAAforward, _pdAAreverse, seq)) {
		return false;
	}
	int totalPeaks = (int) seqIons->size();
	char peakItr;
	PeakList * pPeakList = Spectrum->pPeakList;
	vector<int> mvhKey;
	mvhKey.clear();
	mvhKey.resize(ProNovoConfig::NumIntensityClasses + 1, 0);

	for (int j = 0; j < (int) seqIons->size(); ++j) {
		if (seqIons->at(j) < Spectrum->mzLowerBound || seqIons->at(j) > Spectrum->mzUpperBound) {
			--totalPeaks;
			continue;
		}
		peakItr = pPeakList->findNear(seqIons->at(j), ProNovoConfig::getMassAccuracyFragmentIon());
		if (peakItr != pPeakList->end() && peakItr > 0) {
			++(mvhKey.at(peakItr - 1));
		} else {
			++(mvhKey.at(ProNovoConfig::NumIntensityClasses));
		}
	}

	double mvh = 0.0;
	int fragmentsUnmatched = mvhKey.back();

	if (fragmentsUnmatched != totalPeaks) {
		int fragmentsPredicted = accumulate(mvhKey.begin(), mvhKey.end(), 0);
		int fragmentsMatched = fragmentsPredicted - fragmentsUnmatched;

		if (fragmentsMatched >= ProNovoConfig::MinMatchedFragments) {
			int numVoids = Spectrum->intenClassCounts->back();
			int totalPeakBins = numVoids + pPeakList->size();
			for (int i = 0; i < (int) Spectrum->intenClassCounts->size(); ++i) {
				mvh += lnCombin(Spectrum->intenClassCounts->at(i), mvhKey.at(i));
			}
			mvh -= lnCombin(totalPeakBins, fragmentsPredicted);
			dMvh = -mvh;
		} else {
			return false;
		}
	} else {
		return false;
	}
	/*	if (currentPeptide == "[QITNQDSDARTLIVNNTNGNK]" && Spectrum->iScanId == 120) {
	 cout << "check" << endl;
	 }*/
	return true;
}

bool MVH::ScoreSequenceVsSpectrumSIP(string & currentPeptide, MS2Scan * Spectrum, vector<double>* seqIons, vector<vector<double> > & vvdYionMass,
		vector<vector<double> > & vvdYionProb, vector<vector<double> > & vvdBionMass, vector<vector<double> > & vvdBionProb, double & dMvh,
		vector<char> * seq) {
	if (!CalculateSequenceIonsSIP(currentPeptide, Spectrum->iParentChargeState, bUseSmartPlusThreeModel, seqIons, vvdYionMass, vvdYionProb, vvdBionMass,
			vvdBionProb, seq)) {
		return false;
	}
	int totalPeaks = (int) seqIons->size();
	char peakItr;
	PeakList * pPeakList = Spectrum->pPeakList;
	vector<int> mvhKey;
	mvhKey.clear();
	mvhKey.resize(ProNovoConfig::NumIntensityClasses + 1, 0);

	for (int j = 0; j < (int) seqIons->size(); ++j) {
		if (seqIons->at(j) < Spectrum->mzLowerBound || seqIons->at(j) > Spectrum->mzUpperBound) {
			--totalPeaks;
			continue;
		}
		peakItr = pPeakList->findNear(seqIons->at(j), ProNovoConfig::getMassAccuracyFragmentIon());
		if (peakItr != pPeakList->end() && peakItr > 0) {
			++(mvhKey.at(peakItr - 1));
		} else {
			++(mvhKey.at(ProNovoConfig::NumIntensityClasses));
		}
	}

	double mvh = 0.0;
	int fragmentsUnmatched = mvhKey.back();

	if (fragmentsUnmatched != totalPeaks) {
		int fragmentsPredicted = accumulate(mvhKey.begin(), mvhKey.end(), 0);
		int fragmentsMatched = fragmentsPredicted - fragmentsUnmatched;

		if (fragmentsMatched >= ProNovoConfig::MinMatchedFragments) {
			int numVoids = Spectrum->intenClassCounts->back();
			int totalPeakBins = numVoids + pPeakList->size();
			for (int i = 0; i < (int) Spectrum->intenClassCounts->size(); ++i) {
				mvh += lnCombin(Spectrum->intenClassCounts->at(i), mvhKey.at(i));
			}
			mvh -= lnCombin(totalPeakBins, fragmentsPredicted);
			dMvh = -mvh;
		} else {
			return false;
		}
	} else {
		return false;
	}
	return true;
}

