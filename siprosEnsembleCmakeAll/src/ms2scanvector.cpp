#include "ms2scanvector.h"
#include "SiprosReader.h"

MS2ScanVector::MS2ScanVector(const string & sFT2FilenameInput, const string & sOutputDirectory,
		const string & sConfigFilename, bool bScreenOutput) {
	unsigned int n;
	vector<string> vsSingleResidueNames = ProNovoConfig::vsSingleResidueNames;
	vector<double> vdSingleResidueMasses = ProNovoConfig::vdSingleResidueMasses;
	sFT2Filename = sFT2FilenameInput;
	sConfigFile = sConfigFilename;
	//mass_w = ProNovoConfig::getParentMassWindows();
	setOutputFile(sFT2FilenameInput, sOutputDirectory);
	this->bScreenOutput = bScreenOutput;
	for (n = 0; n < vsSingleResidueNames.size(); ++n)
		mapResidueMass[vsSingleResidueNames[n][0]] = vdSingleResidueMasses[n];

	iOpenMPTaskNum = 0;
}

MS2ScanVector::~MS2ScanVector() {
// the destructors will free memory from vpAllMS2Scans
	vector<MS2Scan *>::iterator it;
	for (it = vpAllMS2Scans.begin(); it != vpAllMS2Scans.end(); ++it) {
		delete (*it);
	}
	vpAllMS2Scans.clear();
	vpPrecursorMasses.clear();
}

void MS2ScanVector::setOutputFile(const string& sFT2FilenameInput, const string& sOutputDirectory) {
	string sLeftFileName, sMiddleFileName, sRightFileName = ".sip";
	size_t pos;
	//cout<<ProNovoConfig::getSearchName()<<endl;
	if ((ProNovoConfig::getSearchName() == "") || (ProNovoConfig::getSearchName() == "Null")
			|| (ProNovoConfig::getSearchName() == "NULL") || (ProNovoConfig::getSearchName() == "null"))
		sRightFileName = ".sip";
	else
		sRightFileName = "." + ProNovoConfig::getSearchName() + ".sip";

	sLeftFileName = sFT2FilenameInput.substr(0, sFT2FilenameInput.length() - 4);
	if (sOutputDirectory == "")
		// if output directory is not specified, it is working directory by default
		sOutputFile = sLeftFileName + sRightFileName;
	else {
		pos = sLeftFileName.rfind(ProNovoConfig::getSeparator());
		if (pos == string::npos)
			sMiddleFileName = sLeftFileName;
		else
			sMiddleFileName = sLeftFileName.substr(pos + 1);
		sOutputFile = sOutputDirectory + ProNovoConfig::getSeparator() + sMiddleFileName + sRightFileName;
	}
}

bool MS2ScanVector::ReadFT2File() {
	bool bReVal, flag_1stScan = true; //flag_1stScan true indicates pMS2Scan is empty
	string sline;
	istringstream input;
	MS2Scan * pMS2Scan = NULL;
	ifstream ft2_stream(sFT2Filename.c_str());
	int tmp_charge;
	double tmp_mz, tmp_intensity;

	bReVal = ft2_stream.is_open();
	if (bReVal) {
		while (!ft2_stream.eof()) {
			sline.clear();
//	    string tmp_sline = sline;
			getline(ft2_stream, sline);
//	    cout << sline <<"wyf"<<endl;
			if (sline == "")
				continue;
			if ((sline.at(0) >= '0') && (sline.at(0) <= '9')) {
				TokenVector words(sline, " \t\n\r");
				if (words.size() < 6)
				// The judgement of MS scan resolution here will be overwritten,
				// we still keep the related codes here for future possible usuage

						{
					pMS2Scan->isMS2HighRes = false;
					pMS2Scan->viCharge.push_back(0);
				} else {
					pMS2Scan->isMS2HighRes = true;
					input.clear();
					input.str(words[5]);
					input >> tmp_charge;
					input.clear();
					pMS2Scan->viCharge.push_back(tmp_charge);
				}
				input.clear();
				input.str(words[0]);
				input >> tmp_mz;
				input.clear();
				if (tmp_mz > ProNovoConfig::maxObservedMz) {
					ProNovoConfig::maxObservedMz = tmp_mz;
				}
				if (tmp_mz < ProNovoConfig::minObservedMz) {
					ProNovoConfig::minObservedMz = tmp_mz;
				}
				pMS2Scan->vdMZ.push_back(tmp_mz);
				input.str(words[1]);
				input >> tmp_intensity;
				input.clear();
				pMS2Scan->vdIntensity.push_back(tmp_intensity);
			} else if (sline.at(0) == 'S') {
				if (flag_1stScan)
					flag_1stScan = false;
				else if (pMS2Scan->vdIntensity.empty())
					delete pMS2Scan;
				else {
					if (ProNovoConfig::getMassAccuracyFragmentIon() < 0.1)
						pMS2Scan->isMS2HighRes = true;
					else
						pMS2Scan->isMS2HighRes = false;
					saveScan(pMS2Scan);
				}
				pMS2Scan = new MS2Scan;
				pMS2Scan->sFT2Filename = sFT2Filename;
				TokenVector words(sline, " \r\t\n");
				input.clear();
				input.str(words[1]);
				input >> pMS2Scan->iScanId;
				input.clear();
				input.str(words[2]);
				input >> pMS2Scan->dParentMZ;
				input.clear();
				pMS2Scan->isMS1HighRes = isMS1HighRes(words[2]);
				// in case no Z Line
				pMS2Scan->iParentChargeState = 0;
				pMS2Scan->dParentNeutralMass = 0;
			} else if (sline.at(0) == 'Z') {
				TokenVector words(sline, " \r\t\n");
				input.clear();
				input.str(words[1]);
				input >> pMS2Scan->iParentChargeState;
				input.clear();
				//input.str(words[2]);
				//input >> pMS2Scan->dParentNeutralMass;
				//input.clear();
			} else if (sline.at(0) == 'I') {
				TokenVector words(sline, " \r\t\n");
				if (words[1] == "ScanType")
					pMS2Scan->setScanType(words[2]);
				// retention time
				if (words[1] == "RetentionTime" || words[1] == "RTime") {
					pMS2Scan->setRTime(words.at(2));
				}
			}
		}
		ft2_stream.clear();
		ft2_stream.close();

		//recognition high or low MS2
//	cout<<ProNovoConfig::getMassAccuracyFragmentIon()<<endl;
		if (ProNovoConfig::getMassAccuracyFragmentIon() < 0.1)
			pMS2Scan->isMS2HighRes = true;
		else
			pMS2Scan->isMS2HighRes = false;
		if (!flag_1stScan) // To avoid empty file
			saveScan(pMS2Scan);
	}
	return bReVal;
}

bool MS2ScanVector::ReadMzmlFile() {
	bool bReVal = false;
	MS2Scan * pMS2Scan;
	vector<Spectrum> * _vSpectra = new vector<Spectrum>();
	SiprosReader::MzmlReader(sFT2Filename, _vSpectra);
	int scannumber = _vSpectra->size();
	int peaknumber = 0;
	double tmp_mz = 0;
	for (int i = 0; i < scannumber; ++i) {
		// cout << _vSpectra->at(i).getScanNumber() << endl;
		pMS2Scan = new MS2Scan;
		pMS2Scan->sFT2Filename = sFT2Filename;
		pMS2Scan->iScanId = _vSpectra->at(i).getScanNumber();
		pMS2Scan->dParentMZ = _vSpectra->at(i).getMZ(0);
		// in case no Z Line
		pMS2Scan->iParentChargeState = 0;
		pMS2Scan->dParentNeutralMass = 0;

		if (_vSpectra->at(i).sizeZ() > 0) {
			pMS2Scan->iParentChargeState = _vSpectra->at(i).atZ(0).z;
		}

		pMS2Scan->setRTime(std::to_string((long double) _vSpectra->at(i).getRTime()));

		peaknumber = _vSpectra->at(i).getPeaks()->size();
		for (int j = 0; j < peaknumber; ++j) {
			tmp_mz = _vSpectra->at(i).getPeaks()->at(j).mz;
			if (tmp_mz > ProNovoConfig::maxObservedMz) {
				ProNovoConfig::maxObservedMz = tmp_mz;
			}
			if (tmp_mz < ProNovoConfig::minObservedMz) {
				ProNovoConfig::minObservedMz = tmp_mz;
			}
			pMS2Scan->vdMZ.push_back(_vSpectra->at(i).getPeaks()->at(j).mz);
			pMS2Scan->vdIntensity.push_back(_vSpectra->at(i).getPeaks()->at(j).intensity);
			pMS2Scan->viCharge.push_back(0);
		}

		if (ProNovoConfig::getMassAccuracyFragmentIon() < 0.1)
			pMS2Scan->isMS2HighRes = true;
		else
			pMS2Scan->isMS2HighRes = false;

		if (pMS2Scan->vdIntensity.empty()) {
			delete pMS2Scan;
		} else {
			saveScan(pMS2Scan);
		}
	}
	_vSpectra->clear();
	delete _vSpectra;
	bReVal = true;
	return bReVal;
}

bool MS2ScanVector::loadMassData() {
	// read all MS2 scans from the file and populate vpAllMS2Scans
	// sort all MS2 scans in vpAllProteins by ascending order of their precursor masses
	// save their precursor masses in vpPrecursorMasses to quick look-up in assignPeptides2Scans()

	bool bReVal; //false when the file fails to be opened.
	vector<MS2Scan *>::iterator it;

	// check the suffix of the MS2 data file
	string filename_str = this->sFT2Filename;
	transform(filename_str.begin(), filename_str.end(), filename_str.begin(), (int (*)(int))tolower);

if(	filename_str.substr(filename_str.rfind('.') + 1) == "mzml") {
		bReVal = ReadMzmlFile();
	} else {
		bReVal = ReadFT2File();
	}

	if (bReVal) {
		sort(vpAllMS2Scans.begin(), vpAllMS2Scans.end(), myless);
		for (it = vpAllMS2Scans.begin(); it < vpAllMS2Scans.end(); it++)
			vpPrecursorMasses.push_back((*it)->dParentNeutralMass);
	}

	return bReVal;
}

bool MS2ScanVector::isMS1HighRes(const std::string& target) {
	bool bReVal = true;
	if (target.at(target.size() - 2) == '.')
		bReVal = false;
	return bReVal;
}

bool MS2ScanVector::ChargeDetermination(const std::vector<double>& vdAllmz, double pmz)
// return true, if the charge state is decided to be one

		{
	bool bReVal = false;
	int iPosi;
	vector<double> vdtempMZ = vdAllmz;
	iPosi = max(0, (int) (((int) vdtempMZ.size()) * 0.05 - 1));
	nth_element(vdtempMZ.begin(), vdtempMZ.begin() + iPosi, vdtempMZ.end(), mygreater);
	if (vdtempMZ.at(iPosi) < pmz)
		bReVal = true;
	return bReVal;
}

bool MS2ScanVector::mygreater(double i, double j) {
	return (i > j);
}

bool MS2ScanVector::myless(MS2Scan* pMS2Scan1, MS2Scan* pMS2Scan2) {
	return (pMS2Scan1->dParentNeutralMass < pMS2Scan2->dParentNeutralMass);
}

void MS2ScanVector::saveScan(MS2Scan * pMS2Scan) { //parentChargeState > 0, save, otherwise, try 1, or 2 and 3
	int j;
	bool bchargeOne;
	MS2Scan * pMS2newScan;
	if (pMS2Scan->iParentChargeState == 0) {
		bchargeOne = ChargeDetermination(pMS2Scan->vdMZ, pMS2Scan->dParentMZ);
		if (bchargeOne) {
			j = 1;
			pMS2newScan = new MS2Scan;
			*pMS2newScan = *pMS2Scan;
			pMS2newScan->iParentChargeState = j;
			pMS2newScan->dParentNeutralMass = (pMS2newScan->dParentMZ) * (pMS2newScan->iParentChargeState)
					- pMS2newScan->iParentChargeState * ProNovoConfig::getProtonMass();
			pMS2newScan->dParentMass = (pMS2newScan->dParentMZ) * (pMS2newScan->iParentChargeState);
			vpAllMS2Scans.push_back(pMS2newScan);
		} else
		// if the charge state is not +1, it could be greater than +3,
		// but current we just try +2 and +3.
		{
			for (j = 2; j <= 3; j++) {
				pMS2newScan = new MS2Scan;
				*pMS2newScan = *pMS2Scan;
				pMS2newScan->iParentChargeState = j;
				pMS2newScan->dParentNeutralMass = (pMS2newScan->dParentMZ) * (pMS2newScan->iParentChargeState)
						- pMS2newScan->iParentChargeState * ProNovoConfig::getProtonMass();
				pMS2newScan->dParentMass = (pMS2newScan->dParentMZ) * (pMS2newScan->iParentChargeState);
				vpAllMS2Scans.push_back(pMS2newScan);
			}
		}
	} else {
		pMS2newScan = new MS2Scan;
		*pMS2newScan = *pMS2Scan;
		pMS2newScan->dParentNeutralMass = (pMS2newScan->dParentMZ) * (pMS2newScan->iParentChargeState)
				- (pMS2newScan->iParentChargeState) * ProNovoConfig::getProtonMass();
		pMS2newScan->dParentMass = (pMS2newScan->dParentMZ) * (pMS2newScan->iParentChargeState);
		vpAllMS2Scans.push_back(pMS2newScan);
	}

	if (pMS2newScan->dParentMass > ProNovoConfig::dMaxMS2ScanMass) {
		ProNovoConfig::dMaxMS2ScanMass = pMS2newScan->dParentMass;
	}
	if (pMS2newScan->iParentChargeState > ProNovoConfig::iMaxPercusorCharge) {
		ProNovoConfig::iMaxPercusorCharge = pMS2newScan->iParentChargeState;
	}
	delete pMS2Scan;
}

void MS2ScanVector::preProcessAllMs2Mvh() {
	CLOCKSTART
	;

	int i, iScanSize;
	iScanSize = (int) vpAllMS2Scans.size();
	if (bScreenOutput)
		cout << "Preprocessing " << vpAllMS2Scans.size() << " scans " << endl;

	int num_threads = omp_get_max_threads();
	vector<multimap<double, double> *> vpIntenSortedPeakPreData;
	for (int j = 0; j < num_threads; ++j) {
		vpIntenSortedPeakPreData.push_back(new multimap<double, double>());
	}

#pragma omp parallel for \
    schedule(guided)
	for (i = 0; i < iScanSize; i++) {
		// vpAllMS2Scans.at(i)->preprocess();
		int iThreadId = omp_get_thread_num();
		vpAllMS2Scans.at(i)->preprocessMvh(vpIntenSortedPeakPreData.at(iThreadId));
	}

	for (int j = 0; j < num_threads; ++j) {
		delete vpIntenSortedPeakPreData.at(j);
	}
	vpIntenSortedPeakPreData.clear();
	int maxPeakBins = 0;
	int iNumSkippedScans = 0;
	for (i = 0; i < iScanSize; i++) {
		if (vpAllMS2Scans.at(i)->bSkip) {
			iNumSkippedScans++;
		}
		if (vpAllMS2Scans.at(i)->totalPeakBins > maxPeakBins) {
			maxPeakBins = vpAllMS2Scans.at(i)->totalPeakBins;
		}
	}
	MVH::initialLnTable(maxPeakBins);

	{
#pragma omp parallel for schedule(guided)
		for (i = 0; i < iScanSize; i++) {
			vpAllMS2Scans.at(i)->sumIntensity();
		}
	}

	if (bScreenOutput) {
		cout << "Preprocessing Done." << endl;
	}

	CLOCKSTOP
	;
}

void MS2ScanVector::preProcessAllMs2WdpSip() {
	int i, iScanSize;
	iScanSize = (int) vpAllMS2Scans.size();
	if (bScreenOutput)
		cout << "Preprocessing " << vpAllMS2Scans.size() << " scans " << endl;

#pragma omp parallel for \
    schedule(guided)
	for (i = 0; i < iScanSize; i++)
		vpAllMS2Scans.at(i)->preprocess();
}

void MS2ScanVector::GetAllRangeFromMass(double dPeptideMass, vector<std::pair<int, int> >& vpPeptideMassRanges)
// all ranges of MS2 scans are stored in  vpPeptideMassWindows
		{
	int i;
	pair<int, int> pairMS2Range;
	pair<int, int> lastPairRange(-100, -100);
	vector<pair<double, double> > vpPeptideMassWindows;
	vpPeptideMassWindows.clear();
	vpPeptideMassRanges.clear();
	ProNovoConfig::getPeptideMassWindows(dPeptideMass, vpPeptideMassWindows);
	for (i = 0; i < (int) vpPeptideMassWindows.size(); i++) {
		pairMS2Range = GetRangeFromMass(vpPeptideMassWindows.at(i).first, vpPeptideMassWindows.at(i).second);
		if ((pairMS2Range.first > -1) && (pairMS2Range.second > -1)) {
			if ((lastPairRange.first < 0) || (lastPairRange.second < 0))
				lastPairRange = pairMS2Range;
			else {
				if (lastPairRange.second > pairMS2Range.first)
					lastPairRange.second = pairMS2Range.second;
				else {
					vpPeptideMassRanges.push_back(lastPairRange);
					lastPairRange = pairMS2Range;
				}
			}
		}
	}
	if ((lastPairRange.first > -1) && (lastPairRange.second > -1))
		vpPeptideMassRanges.push_back(lastPairRange);
}

pair<int, int> MS2ScanVector::GetRangeFromMass(double lb, double ub)
//lb and ub are lower and upper bounds of acceptable parent mass values
		{
	pair<int, int> p;
	int low = 0, high = vpPrecursorMasses.size() - 1, mid;
	double target;
	target = (lb + ub) / 2.0;
	//double lb, ub;  // lower and upper bounds on acceptable parent mass values
	//ub = target + error;
	//lb = target - error;
	while ((high - low) > 1) {
		mid = (high + low) / 2;
		if (vpPrecursorMasses[mid] > target)
			high = mid;
		else
			low = mid;
	}

	// Iterate till we get to the first element > than the lower bound
	int ndx = low;
	//cout<<scan_mass_list_.size()<<endl;
	if (vpPrecursorMasses[ndx] >= lb) {
		while (ndx >= 0 && vpPrecursorMasses[ndx] >= lb)
			ndx--;
		ndx++;
	} else
		while (ndx < (int) vpPrecursorMasses.size() && vpPrecursorMasses[ndx] < lb)
			ndx++;
	if (ndx == (int) vpPrecursorMasses.size() || vpPrecursorMasses[ndx] > ub)
		p = make_pair(-1, -1);
	else {
		low = ndx;
		while (ndx < (int) vpPrecursorMasses.size() && vpPrecursorMasses[ndx] <= ub)
			ndx++;
		high = ndx - 1;
		p = make_pair(low, high);
	}

	if (p.first == -1)
		return p;
	if (vpPrecursorMasses[p.first] < lb)
		cerr << "ERROR L " << vpPrecursorMasses[p.first] << " " << lb << endl;
	if (vpPrecursorMasses[p.second] > ub)
		cerr << "ERROR U " << vpPrecursorMasses[p.second] << " " << ub << endl;
	return p;
}

/*
 void MS2ScanVector::assignPeptides2Scans(Peptide * currentPeptide)
 {
 int i;
 pair<int, int> pairMS2Range;
 //   cout<<currentPeptide->getPeptideMass()<<endl;
 pairMS2Range = GetRangeFromMass(currentPeptide->getPeptideMass(),
 ProNovoConfig::getMassAccuracyParentIon());
 if ((pairMS2Range.first > -1) && (pairMS2Range.second > -1))
 for (i= pairMS2Range.first; i<= pairMS2Range.second; i++)
 vpAllMS2Scans.at(i)->vpPeptides.push_back(currentPeptide);
 }*/

bool MS2ScanVector::assignPeptides2Scans(Peptide * currentPeptide) {
	int i, j;
	bool bAssigned = false;
	vector<pair<int, int> > vpPeptideMassRanges;
	pair<int, int> pairMS2Range;

	GetAllRangeFromMass(currentPeptide->getPeptideMass(), vpPeptideMassRanges);

	for (j = 0; j < (int) vpPeptideMassRanges.size(); j++) {
		pairMS2Range = vpPeptideMassRanges.at(j);
		if ((pairMS2Range.first > -1) && (pairMS2Range.second > -1)) {
			for (i = pairMS2Range.first; i <= pairMS2Range.second; i++) {
				vpAllMS2Scans.at(i)->vpPeptides.push_back(currentPeptide);
			}
			bAssigned = true;
		}
	}
	return bAssigned;
}

void MS2ScanVector::processPeptideArrayWdpSip(vector<Peptide*>& vpPeptideArray) {
	int i, iPeptideArraySize, iScanSize;
	iPeptideArraySize = (int) vpPeptideArray.size();

//    for (int i=0; i< (int) vpPeptideArray.size(); i++)
//      cout<<vpPeptideArray.at(i)->getPeptideSeq() <<"\t"
//	  <<vpPeptideArray.at(i)->getOriginalPeptideSeq()<<"\t"
//	  <<vpPeptideArray.at(i)->getProteinName()<<"\t"
//	  <<vpPeptideArray.at(i)->getBeginPosProtein()<<"\t"
//	  <<vpPeptideArray.at(i)->getPeptideMass()<<endl;

	//cout<<"calculating fragments of "<<vpPeptideArray.size()<<"  peptides"<<endl;
#pragma omp parallel for \
    shared(vpPeptideArray) private(i) \
    schedule(guided)

	for (i = 0; i < iPeptideArraySize; i++)
		vpPeptideArray[i]->preprocessing(vpAllMS2Scans.at(0)->isMS2HighRes, mapResidueMass);
	//vpPeptideArray[i]->calculateExpectedFragments(mapResidueMass);

	//cout<<"scoring "<<vpPeptideArray.size()<<"  peptides"<<endl;
	// every MS2 scans scores their matched peptides

	iScanSize = (int) vpAllMS2Scans.size();
#pragma omp parallel for  \
    schedule(guided)

	for (i = 0; i < iScanSize; i++)
		vpAllMS2Scans[i]->scorePeptides();

	// free memory of all peptide objects
	for (i = 0; i < (int) vpPeptideArray.size(); i++)
		delete vpPeptideArray[i];

	// empty peptide array
	vpPeptideArray.clear();

}

void MS2ScanVector::processPeptideArrayMvh(vector<Peptide*>& vpPeptideArray) {
	int i, iPeptideArraySize, iScanSize;
	iPeptideArraySize = (int) vpPeptideArray.size();

#pragma omp parallel for \
    shared(vpPeptideArray) private(i) \
    schedule(guided)

	for (i = 0; i < iPeptideArraySize; i++) {
		vpPeptideArray.at(i)->preprocessingMVH();
	}

	// every MS2 scans scores their matched peptides

	iScanSize = (int) vpAllMS2Scans.size();
#pragma omp parallel for  \
    schedule(guided)

	for (i = 0; i < iScanSize; i++) {
		int iThreadId = omp_get_thread_num();
		vpAllMS2Scans[i]->scorePeptidesMVH(psequenceIonMasses.at(iThreadId), _ppdAAforward.at(iThreadId),
				_ppdAAreverse.at(iThreadId), pSeqs.at(iThreadId));
	}

	// free memory of all peptide objects
	for (i = 0; i < (int) vpPeptideArray.size(); i++)
		delete vpPeptideArray.at(i);

	// empty peptide array
	vpPeptideArray.clear();
}

void MS2ScanVector::processPeptideArrayMvhTask(vector<Peptide*>& vpPeptideArray, omp_lock_t * pLck) {
	int i, j, k, iPeptideArraySize;
	vector<pair<int, int> > vpPeptideMassRanges;
	pair<int, int> pairMS2Range;
	bool bAssigned = false;
	bool bMerged = false, bScored = false, bProcessed = false;
	double dMvh = 0;
	int iThreadId = omp_get_thread_num();

	iPeptideArraySize = vpPeptideArray.size();

	for (i = 0; i < iPeptideArraySize; i++) {
		bAssigned = false;
		bProcessed = false;
		// assign peptide to scan
		GetAllRangeFromMass(vpPeptideArray.at(i)->getPeptideMass(), vpPeptideMassRanges);
		for (j = 0; j < (int) vpPeptideMassRanges.size(); j++) {
			pairMS2Range = vpPeptideMassRanges.at(j);
			if ((pairMS2Range.first > -1) && (pairMS2Range.second > -1)) {
				bAssigned = true;
				if (bAssigned && !bProcessed) {
					vpPeptideArray.at(i)->preprocessingMVH();
					bProcessed = true;
				}
				// calculate the score
				for (k = pairMS2Range.first; k <= pairMS2Range.second; ++k) {
					if (vpAllMS2Scans.at(k)->bSkip) {
						continue;
					}
					omp_set_lock(&(pLck[k]));
					bMerged = vpAllMS2Scans.at(k)->mergePeptide(vpAllMS2Scans.at(k)->vpWeightSumTopPeptides,
							vpPeptideArray.at(i)->getPeptideForScoring(), vpPeptideArray.at(i)->getProteinName());
					omp_unset_lock(&(pLck[k]));
					// not merged then score
					if (!bMerged) {
						bScored = MVH::ScoreSequenceVsSpectrum(vpPeptideArray.at(i)->sNeutralLossPeptide,
								vpAllMS2Scans.at(k), psequenceIonMasses.at(iThreadId), _ppdAAforward.at(iThreadId),
								_ppdAAreverse.at(iThreadId), dMvh, pSeqs.at(iThreadId));
						if (bScored) {
							omp_set_lock(&(pLck[k]));
							vpAllMS2Scans.at(k)->saveScore(dMvh, vpPeptideArray.at(i),
									vpAllMS2Scans.at(k)->vpWeightSumTopPeptides,
									 "MVH");
							omp_unset_lock(&(pLck[k]));
						}
					}
				}
			}
		}
	}

	// free memory of all peptide objects
	for (i = 0; i < iPeptideArraySize; i++) {
		delete vpPeptideArray.at(i);
	}

	// empty peptide array
	vpPeptideArray.clear();

	omp_set_lock(&lckOpenMpTaskNum);
	iOpenMPTaskNum--;
	if (iOpenMPTaskNum < TASKRESUME_SIZE) {
		if (omp_test_lock(&lckOpenMpTaskNumHalfed) == 0) {
			// lock is set'
			omp_unset_lock(&lckOpenMpTaskNumHalfed);
			// cout << endl << "Resume" << endl;
		} else {
			omp_unset_lock(&lckOpenMpTaskNumHalfed);
		}
	}
	omp_unset_lock(&lckOpenMpTaskNum);
}

void MS2ScanVector::searchDatabaseMvh() {
	CLOCKSTART
	;

	ProteinDatabase myProteinDatabase(bScreenOutput);
	vector<Peptide *> vpPeptideArray;
	Peptide * currentPeptide;
	myProteinDatabase.loadDatabase();
	this->preMvh();
	if (myProteinDatabase.getFirstProtein()) {
		currentPeptide = new Peptide;
		// get one peptide from the database at a time, until there is no more peptide
		while (myProteinDatabase.getNextPeptide(currentPeptide)) {
			// assign the pointers of peptides to appropriete MS2Scans
			if (assignPeptides2Scans(currentPeptide)) {
				// save the new peptide to the array
				vpPeptideArray.push_back(currentPeptide);
				if (currentPeptide->getPeptideMass() > ProNovoConfig::dMaxPeptideMass) {
					ProNovoConfig::dMaxPeptideMass = currentPeptide->getPeptideMass();
				}
			} else {
				delete currentPeptide;
			}
			// create a new peptide for the next iteration
			currentPeptide = new Peptide;
			// when the vpPeptideArray is full
			if (vpPeptideArray.size() >= PEPTIDE_ARRAY_SIZE)
				processPeptideArrayMvh(vpPeptideArray);
		}
		// the last peptide object is an empty object and need to be deleted
		delete currentPeptide;
		// there are still unprocessed peptides in the vpPeptideArray
		// need to process them in the same manner
		// the following code is the same as inside if(vpPeptideArray.size() >= PEPTIDE_ARRAY_SIZE )
//    cout<<vpPeptideArray.size()<<endl;
		if (!vpPeptideArray.empty())
			processPeptideArrayMvh(vpPeptideArray);
	}
	CLOCKSTOP
	;

	this->postMvh();
	MVH::destroyLnTable();
	PeptideUnit::iNumScores = 1;
	cout << endl << "Search done." << endl;
}

void MS2ScanVector::searchDatabaseWdpSip() {

	ProteinDatabase myProteinDatabase(bScreenOutput);
	vector<Peptide *> vpPeptideArray;
	Peptide * currentPeptide;
	myProteinDatabase.loadDatabase();

	if (myProteinDatabase.getFirstProtein()) {
		currentPeptide = new Peptide;
		// get one peptide from the database at a time, until there is no more peptide
		while (myProteinDatabase.getNextPeptide(currentPeptide)) {

			// assign the pointers of peptides to appropriete MS2Scans
			if (assignPeptides2Scans(currentPeptide)) {
				// save the new peptide to the array
				vpPeptideArray.push_back(currentPeptide);
			} else {
				// not assigned, delete
				delete currentPeptide;
			}
			// create a new peptide for the next iteration
			currentPeptide = new Peptide;
			// when the vpPeptideArray is full
			if (vpPeptideArray.size() >= PEPTIDE_ARRAY_SIP_SIZE) {
				processPeptideArrayWdpSip(vpPeptideArray);
			}
		}
		// the last peptide object is an empty object and need to be deleted
		delete currentPeptide;
		// there are still unprocessed peptides in the vpPeptideArray
		// need to process them in the same manner
		// the following code is the same as inside if(vpPeptideArray.size() >= PEPTIDE_ARRAY_SIZE )
		if (!vpPeptideArray.empty()) {
			processPeptideArrayWdpSip(vpPeptideArray);
		}

		PeptideUnit::iNumScores = 1;
	}

}

void MS2ScanVector::searchDatabaseMvhTask() {
	ProteinDatabase myProteinDatabase(bScreenOutput);
	vector<Peptide *> vpPeptideArray;
	Peptide * currentPeptide;
	myProteinDatabase.loadDatabase();
	this->preMvh();

	int iScanNum = this->vpAllMS2Scans.size();
	omp_lock_t * pLck = new omp_lock_t[iScanNum];
	iOpenMPTaskNum = 0;
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < iScanNum; ++i) {
		omp_init_lock(&(pLck[i]));
	}
	omp_init_lock(&lckOpenMpTaskNum);
	omp_init_lock(&lckOpenMpTaskNumHalfed);
	// omp_set_lock(&lckOpenMpTaskNumHalfed);
#pragma omp parallel
	{
#pragma omp single
		{
			if (myProteinDatabase.getFirstProtein()) {
				currentPeptide = new Peptide();
				// get one peptide from the database at a time, until there is no more peptide
				while (myProteinDatabase.getNextPeptide(currentPeptide)) {
					if (currentPeptide->getPeptideMass() > ProNovoConfig::dMaxPeptideMass) {
						ProNovoConfig::dMaxPeptideMass = currentPeptide->getPeptideMass();
					}
					vpPeptideArray.push_back(currentPeptide);
					if (vpPeptideArray.size() >= TASKPEPTIDE_ARRAY_SIZE) {
						omp_set_lock(&lckOpenMpTaskNum);
						iOpenMPTaskNum++;
						omp_unset_lock(&lckOpenMpTaskNum);
#pragma omp task firstprivate(vpPeptideArray)
						{
							processPeptideArrayMvhTask(vpPeptideArray, pLck);
						}

						vpPeptideArray.clear();

						if (iOpenMPTaskNum >= TASKWAIT_SIZE) {
							// cout << endl << "Stop" << endl;
							omp_set_lock(&lckOpenMpTaskNumHalfed);
						}
					}
					currentPeptide = new Peptide();
				}
				delete currentPeptide;
			}
			if (!vpPeptideArray.empty()) {
#pragma omp task firstprivate(vpPeptideArray)
				{
					processPeptideArrayMvhTask(vpPeptideArray, pLck);
				}
				vpPeptideArray.clear();
			}
#pragma omp taskwait
		}
	}

	omp_destroy_lock(&lckOpenMpTaskNumHalfed);
	omp_destroy_lock(&lckOpenMpTaskNum);

#pragma omp parallel for schedule(guided)
	for (int i = 0; i < iScanNum; ++i) {
		omp_destroy_lock(&(pLck[i]));
	}

	this->postMvh();
	MVH::destroyLnTable();
	PeptideUnit::iNumScores = 1;
	cout << endl << "Search done." << endl;
}

void MS2ScanVector::startProcessingMvh() {
	// Preprocessing all MS2 scans by mult-threading
	preProcessAllMs2Mvh();

	// Search all MS2 scans against the database by mult-threading
	searchDatabaseMvh();

	// Postprocessing all MS2 scans' results by mult-threading
	postProcessAllMs2WdpXcorr();

	// write results to a SIP file
	writeOutputEnsemble();

}

void MS2ScanVector::startProcessingMvhTask() {
	// Preprocessing all MS2 scans by mult-threading
	preProcessAllMs2Mvh();

	// Search all MS2 scans against the database by mult-threading
	searchDatabaseMvhTask();

	// Postprocessing all MS2 scans' results by mult-threading
	postProcessAllMs2WdpXcorr();

	// write results to a SIP file
	// writeOutput();
	writeOutputEnsemble();
}

void MS2ScanVector::startProcessingWdpSip() {

	// Preprocessing all MS2 scans by multi-threading
	preProcessAllMs2WdpSip();

	// Search all MS2 scans against the database by mult-threading
	searchDatabaseWdpSip();

	// Postprocessing all MS2 scans' results by mult-threading
	postProcessAllMs2MvhXcorr();

	// write results to a SIP file
	writeOutputEnsemble();

}

void MS2ScanVector::postProcessAllMs2WdpXcorr() {

	int i, iScanSize;
	iScanSize = (int) vpAllMS2Scans.size();

	postProcessAllMs2Xcorr();

	postProcessAllMs2Wdp();

#pragma omp parallel for \
    schedule(guided)
	for (i = 0; i < iScanSize; i++) {
		vpAllMS2Scans.at(i)->scoreFeatureCalculation();
	}
}

void MS2ScanVector::postProcessAllMs2MvhXcorr() {
	int i, iScanSize;
	iScanSize = (int) vpAllMS2Scans.size();

	postProcessAllMs2XcorrSip();

	postProcessAllMs2MvhSip();

#pragma omp parallel for \
	    schedule(guided)
	for (i = 0; i < iScanSize; i++) {
		vpAllMS2Scans.at(i)->scoreFeatureCalculationWDPSip();
	}
}

void MS2ScanVector::postProcessAllMs2Wdp() {
	vector<vector<vector<double> > *> vpvvdYionMass;
	vector<vector<vector<double> > *> vpvvdYionProb;
	vector<vector<vector<double> > *> vpvvdBionMass;
	vector<vector<vector<double> > *> vpvvdBionProb;
	vector<vector<double> *> vpvdYionMass;
	vector<vector<double> *> vpvdBionMass;
	num_max_threads = omp_get_max_threads();
	for (int i = 0; i < num_max_threads; ++i) {
		vpvvdYionMass.push_back(new vector<vector<double> >());
		vpvvdYionProb.push_back(new vector<vector<double> >());
		vpvvdBionMass.push_back(new vector<vector<double> >());
		vpvvdBionProb.push_back(new vector<vector<double> >());
		vpvdYionMass.push_back(new vector<double>());
		vpvdBionMass.push_back(new vector<double>());
	}

	int iScanSize;
	iScanSize = (int) vpAllMS2Scans.size();

#pragma omp parallel for \
    schedule(guided)
	for (int i = 0; i < iScanSize; i++) {
		vpAllMS2Scans.at(i)->preprocess();
		bool bHighRes = vpAllMS2Scans.at(i)->isMS2HighRes;
		int iThreadId = omp_get_thread_num();
		for (int j = 0; j < ((int) vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.size()); j++) {
			Peptide::preprocessing(vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sPeptideForScoring, bHighRes,
					mapResidueMass, vpvvdYionMass.at(iThreadId), vpvvdYionProb.at(iThreadId),
					vpvvdBionMass.at(iThreadId), vpvvdBionProb.at(iThreadId), vpvdYionMass.at(iThreadId),
					vpvdBionMass.at(iThreadId));
			double dWeightSum = 0;
			if (bHighRes) {
				dWeightSum = vpAllMS2Scans.at(i)->scoreWeightSumHighMS2(
						&(vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sPeptideForScoring),
						vpvvdYionMass.at(iThreadId), vpvvdYionProb.at(iThreadId), vpvvdBionMass.at(iThreadId),
						vpvvdBionProb.at(iThreadId));
			} else {
				dWeightSum = vpAllMS2Scans.at(i)->scoreWeightSum(
						&(vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sPeptideForScoring),
						vpvdYionMass.at(iThreadId), vpvdBionMass.at(iThreadId));
			}
			vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores.push_back(dWeightSum);
		}

	}

	for (int i = 0; i < num_max_threads; ++i) {
		delete vpvvdYionMass.at(i);
		delete vpvvdYionProb.at(i);
		delete vpvvdBionMass.at(i);
		delete vpvvdBionProb.at(i);
		delete vpvdYionMass.at(i);
		delete vpvdBionMass.at(i);
	}
	++PeptideUnit::iNumScores;
}

void MS2ScanVector::postProcessAllMs2Xcorr() {
	//MH: Must be equal to largest possible array
	int iArraySizePreprocess = (int) ((ProNovoConfig::dMaxMS2ScanMass + 3 + 2.0)
			* ProNovoConfig::dHighResInverseBinWidth);
	//MH: Must be equal to largest possible array
	int iArraySizeScore = (int) ((ProNovoConfig::dMaxPeptideMass + 100) * ProNovoConfig::dHighResInverseBinWidth);
	CometSearchMod::iArraySizePreprocess = iArraySizePreprocess;
	CometSearchMod::iArraySizeScore = iArraySizeScore;
	CometSearchMod::iDimesion2 = 9;
	CometSearchMod::iMAX_PEPTIDE_LEN = MAX_PEPTIDE_LEN;
	CometSearchMod::iMaxPercusorCharge = ProNovoConfig::iMaxPercusorCharge + 1;

	vector<double *> vpdTmpRawData;
	vector<double *> vpdTmpFastXcorrData;
	vector<double *> vpdTmpCorrelationData;
	vector<double *> vpdTmpSmoothedSpectrum;
	vector<double *> vpdTmpPeakExtracted;
	vector<bool *> vpbDuplFragment;
	vector<double *> v_pdAAforward;
	vector<double *> v_pdAAreverse;
	vector<unsigned int ***> v_uiBinnedIonMasses;
	num_max_threads = omp_get_max_threads();
	for (int i = 0; i < num_max_threads; ++i) {
		vpdTmpRawData.push_back(new double[iArraySizePreprocess]());
		vpdTmpFastXcorrData.push_back(new double[iArraySizePreprocess]());
		vpdTmpCorrelationData.push_back(new double[iArraySizePreprocess]());
		vpdTmpSmoothedSpectrum.push_back(new double[iArraySizePreprocess]());
		vpdTmpPeakExtracted.push_back(new double[iArraySizePreprocess]());
		vpbDuplFragment.push_back(new bool[iArraySizeScore]());
		v_pdAAforward.push_back(new double[MAX_PEPTIDE_LEN]());
		v_pdAAreverse.push_back(new double[MAX_PEPTIDE_LEN]());
		unsigned int *** _uiBinnedIonMasses = new unsigned int**[ProNovoConfig::iMaxPercusorCharge + 1]();
		for (int ii = 0; ii < ProNovoConfig::iMaxPercusorCharge + 1; ii++) {
			_uiBinnedIonMasses[ii] = new unsigned int*[9]();
			for (int j = 0; j < 9; j++) {
				_uiBinnedIonMasses[ii][j] = new unsigned int[MAX_PEPTIDE_LEN]();
			}
		}
		v_uiBinnedIonMasses.push_back(_uiBinnedIonMasses);
	}

	int iScanSize;
	iScanSize = (int) vpAllMS2Scans.size();

#pragma omp parallel for \
    schedule(guided)
	for (int i = 0; i < iScanSize; i++) {
		int iThreadId = omp_get_thread_num();
		struct Query * pQuery = new Query();
		if (!CometSearchMod::Preprocess(pQuery, vpAllMS2Scans.at(i), vpdTmpRawData.at(iThreadId),
				vpdTmpFastXcorrData.at(iThreadId), vpdTmpCorrelationData.at(iThreadId),
				vpdTmpSmoothedSpectrum.at(iThreadId), vpdTmpPeakExtracted.at(iThreadId))) {
			cout << "Error Post Xcorr." << endl;
			exit(1);
		} else {
			if (vpAllMS2Scans.at(i)->pQuery != NULL) {
				delete vpAllMS2Scans.at(i)->pQuery;
				vpAllMS2Scans.at(i)->pQuery = NULL;
			}
			vpAllMS2Scans.at(i)->pQuery = pQuery;
			for (int j = 0; j < (int) vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.size(); j++) {
				double dXcorr = 0;
				CometSearchMod::ScorePeptides(&(vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sPeptideForScoring),
						vpbDuplFragment.at(iThreadId), v_pdAAforward.at(iThreadId), v_pdAAreverse.at(iThreadId),
						vpAllMS2Scans.at(i), v_uiBinnedIonMasses.at(iThreadId), dXcorr);
				vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores.push_back(dXcorr);
			}
			delete pQuery;
			vpAllMS2Scans.at(i)->pQuery = NULL;
		}

	}

	for (int i = 0; i < num_max_threads; ++i) {
		delete[] vpdTmpRawData.at(i);
		delete[] vpdTmpFastXcorrData.at(i);
		delete[] vpdTmpCorrelationData.at(i);
		delete[] vpdTmpSmoothedSpectrum.at(i);
		delete[] vpdTmpPeakExtracted.at(i);
		delete[] vpbDuplFragment.at(i);
		delete[] v_pdAAforward.at(i);
		delete[] v_pdAAreverse.at(i);
		for (int ii = 0; ii < ProNovoConfig::iMaxPercusorCharge + 1; ii++) {
			for (int j = 0; j < 9; j++) {
				delete[] v_uiBinnedIonMasses.at(i)[ii][j];
			}
			delete[] v_uiBinnedIonMasses.at(i)[ii];
		}
		delete[] v_uiBinnedIonMasses.at(i);
	}
	++PeptideUnit::iNumScores;
}

void MS2ScanVector::postProcessAllMs2XcorrSip() {
	//MH: Must be equal to largest possible array
	int iArraySizePreprocess = (int) ((ProNovoConfig::dMaxMS2ScanMass + 3 + 2.0)
			* ProNovoConfig::dHighResInverseBinWidth);
	//MH: Must be equal to largest possible array
	int iArraySizeScore = (int) ((ProNovoConfig::dMaxPeptideMass + 100) * ProNovoConfig::dHighResInverseBinWidth);
	CometSearchMod::iArraySizePreprocess = iArraySizePreprocess;
	CometSearchMod::iArraySizeScore = iArraySizeScore;
	CometSearchMod::iDimesion2 = 9;
	CometSearchMod::iMAX_PEPTIDE_LEN = MAX_PEPTIDE_LEN;
	CometSearchMod::iMaxPercusorCharge = ProNovoConfig::iMaxPercusorCharge + 1;

	vector<double *> vpdTmpRawData;
	vector<double *> vpdTmpFastXcorrData;
	vector<double *> vpdTmpCorrelationData;
	vector<double *> vpdTmpSmoothedSpectrum;
	vector<double *> vpdTmpPeakExtracted;
	vector<bool *> vpbDuplFragment;
	vector<double *> v_pdAAforward;
	vector<double *> v_pdAAreverse;
	vector<unsigned int ***> v_uiBinnedIonMasses;
	num_max_threads = omp_get_max_threads();
	for (int i = 0; i < num_max_threads; ++i) {
		vpdTmpRawData.push_back(new double[iArraySizePreprocess]());
		vpdTmpFastXcorrData.push_back(new double[iArraySizePreprocess]());
		vpdTmpCorrelationData.push_back(new double[iArraySizePreprocess]());
		vpdTmpSmoothedSpectrum.push_back(new double[iArraySizePreprocess]());
		vpdTmpPeakExtracted.push_back(new double[iArraySizePreprocess]());
		vpbDuplFragment.push_back(new bool[iArraySizeScore]());
		v_pdAAforward.push_back(new double[MAX_PEPTIDE_LEN]());
		v_pdAAreverse.push_back(new double[MAX_PEPTIDE_LEN]());
		unsigned int *** _uiBinnedIonMasses = new unsigned int**[ProNovoConfig::iMaxPercusorCharge + 1]();
		for (int ii = 0; ii < ProNovoConfig::iMaxPercusorCharge + 1; ii++) {
			_uiBinnedIonMasses[ii] = new unsigned int*[9]();
			for (int j = 0; j < 9; j++) {
				_uiBinnedIonMasses[ii][j] = new unsigned int[MAX_PEPTIDE_LEN]();
			}
		}
		v_uiBinnedIonMasses.push_back(_uiBinnedIonMasses);
	}

	int iScanSize;
	iScanSize = (int) vpAllMS2Scans.size();

	preXcorr();

#pragma omp parallel for \
	    schedule(guided)
	for (int i = 0; i < iScanSize; i++) {
		int iThreadId = omp_get_thread_num();
		struct Query * pQuery = new Query();
		if (!CometSearchMod::Preprocess(pQuery, vpAllMS2Scans.at(i), vpdTmpRawData.at(iThreadId),
				vpdTmpFastXcorrData.at(iThreadId), vpdTmpCorrelationData.at(iThreadId),
				vpdTmpSmoothedSpectrum.at(iThreadId), vpdTmpPeakExtracted.at(iThreadId))) {
			cout << "Error Post Xcorr." << endl;
			exit(1);
		} else {
			if (vpAllMS2Scans.at(i)->pQuery != NULL) {
				delete vpAllMS2Scans.at(i)->pQuery;
				vpAllMS2Scans.at(i)->pQuery = NULL;
			}
			vpAllMS2Scans.at(i)->pQuery = pQuery;
			for (int j = 0; j < (int) vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.size(); j++) {
				double dXcorr = 0;
				PeptideUnit * pepUnit = vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j);
				CometSearchMod::ScorePeptidesSIPNoCancelOut(pepUnit->vvdYionMass, pepUnit->vvdYionProb,
						pepUnit->vvdBionMass, pepUnit->vvdBionProb, vpAllMS2Scans.at(i),
						vvpbDuplFragmentGlobal.at(iThreadId), vvdBinnedIonMassesGlobal.at(iThreadId),
						vvdBinGlobal.at(iThreadId), dXcorr);
				vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores.push_back(dXcorr);
			}
			delete pQuery;
			vpAllMS2Scans.at(i)->pQuery = NULL;
		}

	}

	postXcorr();

	for (int i = 0; i < num_max_threads; ++i) {
		delete[] vpdTmpRawData.at(i);
		delete[] vpdTmpFastXcorrData.at(i);
		delete[] vpdTmpCorrelationData.at(i);
		delete[] vpdTmpSmoothedSpectrum.at(i);
		delete[] vpdTmpPeakExtracted.at(i);
		delete[] vpbDuplFragment.at(i);
		delete[] v_pdAAforward.at(i);
		delete[] v_pdAAreverse.at(i);
		for (int ii = 0; ii < ProNovoConfig::iMaxPercusorCharge + 1; ii++) {
			for (int j = 0; j < 9; j++) {
				delete[] v_uiBinnedIonMasses.at(i)[ii][j];
			}
			delete[] v_uiBinnedIonMasses.at(i)[ii];
		}
		delete[] v_uiBinnedIonMasses.at(i);
	}
	++PeptideUnit::iNumScores;
}

double roundMy(double f, int precision) {
	if (f == 0.0f)
		return +0.0f;

	double multiplier = pow(10.0, (double) precision); // moves f over <precision> decimal places
	f *= multiplier;
	f = floor(f + 0.5f);
	return f / multiplier;
}

void MS2ScanVector::postProcessAllMs2MvhSip() {
	int num_threads = omp_get_max_threads();
	vector<multimap<double, double> *> vpIntenSortedPeakPreData;
	for (int j = 0; j < num_threads; ++j) {
		vpIntenSortedPeakPreData.push_back(new multimap<double, double>());
	}
	this->preMvh();
	int iScanSize;
	iScanSize = (int) vpAllMS2Scans.size();

	double totalPeakSpace = ProNovoConfig::maxObservedMz - ProNovoConfig::minObservedMz;
	int totalPeakBins = (int) roundMy(totalPeakSpace / (ProNovoConfig::getMassAccuracyFragmentIon() * 2.0f), 0);
	MVH::initialLnTable(totalPeakBins);

#pragma omp parallel for \
    schedule(guided)
	for (int i = 0; i < iScanSize; i++) {
		int iThreadId = omp_get_thread_num();
		vpAllMS2Scans.at(i)->preprocessMvh(vpIntenSortedPeakPreData.at(iThreadId));
		if (!vpAllMS2Scans.at(i)->bSkip) {
			for (int j = 0; j < (int) vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.size(); j++) {
				double dMvh = 0;
				PeptideUnit * pepUnit = vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j);
				MVH::ScoreSequenceVsSpectrumSIP(pepUnit->sPeptideForScoring, vpAllMS2Scans.at(i),
						psequenceIonMasses.at(iThreadId), pepUnit->vvdYionMass, pepUnit->vvdYionProb,
						pepUnit->vvdBionMass, pepUnit->vvdBionProb, dMvh, pSeqs.at(iThreadId));
				vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores.push_back(dMvh);
			}
		} else {
			for (int j = 0; j < (int) vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.size(); j++) {
				double dMvh = 0;
				vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores.push_back(dMvh);
			}
		}
		if (vpAllMS2Scans.at(i)->pPeakList != NULL) {
			delete vpAllMS2Scans.at(i)->pPeakList;
			vpAllMS2Scans.at(i)->pPeakList = NULL;
		}
		if (vpAllMS2Scans.at(i)->intenClassCounts != NULL) {
			delete vpAllMS2Scans.at(i)->intenClassCounts;
			vpAllMS2Scans.at(i)->intenClassCounts = NULL;
		}
	}

	for (int j = 0; j < num_threads; ++j) {
		delete vpIntenSortedPeakPreData.at(j);
	}
	vpIntenSortedPeakPreData.clear();
	this->postMvh();
	MVH::destroyLnTable();
	++PeptideUnit::iNumScores;
}

bool MS2ScanVector::mylessScanId(MS2Scan* pMS2Scan1, MS2Scan* pMS2Scan2) {
	return (pMS2Scan1->iScanId < pMS2Scan2->iScanId);
}

string MS2ScanVector::ParsePath(string sPath)
// return the filename without path
		{
	string sTailFileName;
	size_t iPosition;
	iPosition = sPath.rfind(ProNovoConfig::getSeparator());
	if (iPosition == string::npos)
		sTailFileName = sPath;
	else
		sTailFileName = sPath.substr(iPosition + 1);
	return sTailFileName;
}

void MS2ScanVector::writeOutputEnsemble() {

	int i, j;
	int k = 0;
	ofstream outputFile;
	ifstream configStream;
	string sline, sTailFT2FileName;

	sTailFT2FileName = ParsePath(sFT2Filename);
	sort(vpAllMS2Scans.begin(), vpAllMS2Scans.end(), mylessScanId);
	string sOutputFileTxt = sOutputFile.substr(0, sOutputFile.length() - 4);
	sOutputFileTxt += ".Spe2Pep.txt";
	outputFile.open(sOutputFileTxt.c_str());

	configStream.open(sConfigFile.c_str());
	while (!configStream.eof()) {
		sline.clear();
		getline(configStream, sline);
		outputFile << "#\t" << sline << endl;
	}

	//Spectrum level head
	outputFile << "+\tFilename\tScanNumber\tParentCharge\tMeasuredParentMass"
			<< "\tScanType\tSearchName\tRetentionTime" << endl;
	//PSM level head
	outputFile << "*\tIdentifiedPeptide\tOriginalPeptide\tCalculatedParentMass " << "\tMVH\tXcorr\tWDP\tProteinNames"
			<< endl;

	for (i = 0; i < ((int)vpAllMS2Scans.size()); i++) {
		if (!vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.empty()) {
			//Spectrum level information
			outputFile << "+";
			outputFile << "\t" << sTailFT2FileName;
			outputFile << "\t" << vpAllMS2Scans.at(i)->iScanId;
			outputFile << "\t" << vpAllMS2Scans.at(i)->iParentChargeState;
			outputFile << "\t" << setiosflags(ios::fixed) << setprecision(5) << vpAllMS2Scans.at(i)->dParentNeutralMass;
			outputFile << "\t" << vpAllMS2Scans.at(i)->getScanType();
			outputFile << "\t" << ProNovoConfig::getSearchName();
			//sInt sum of intensity
			// outputFile << "\t" << vpAllMS2Scans.at(i)->dSumIntensity;
			//maxInt
			// outputFile << "\t" << vpAllMS2Scans.at(i)->dMaxIntensity;
			// Retention Time
			outputFile << "\t" << vpAllMS2Scans.at(i)->getRTime();
			outputFile << endl;
			//PSM level head
			for (j = 0; j < ((int)vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.size()); j++) {
				if (!vpAllMS2Scans.at(i)->isAnyScoreInTopN(j, ProNovoConfig::INTTOPKEEP)) {
					continue;
				}
				outputFile << "*\t";
				//IdentifiedPeptide
				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifyPrefix != '-') {
					outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifyPrefix;
				}
				outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sIdentifiedPeptide;
				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifySuffix != '-') {
					outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifySuffix;
				}
				outputFile << "\t";
				//OriginalPeptide
				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalPrefix != '-') {
					outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalPrefix;
				}
				outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sOriginalPeptide;
				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalSuffix != '-') {
					outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalSuffix;
				}

				//CalculatedParentMass
				outputFile << "\t";
				outputFile << setiosflags(ios::fixed) << setprecision(4)
						<< vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->dPepNeutralMass << "\t";
				//MVH, Xcorr, WDP
				for (k = 0; k < ((int)vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores.size()); k++) {
					outputFile << setiosflags(ios::fixed) << setprecision(4)
							<< vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores.at(k) << "\t";
				}
				//ProteinNames
				outputFile << "{" << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sProteinNames << "}" << endl;
			}
		}
	}
	configStream.clear();
	configStream.close();
	outputFile.close();
}

void MS2ScanVector::preMvh() {
	num_max_threads = omp_get_max_threads();
	for (int i = 0; i < num_max_threads; ++i) {
		_ppdAAforward.push_back(new vector<double>());
		_ppdAAreverse.push_back(new vector<double>());
		psequenceIonMasses.push_back(new vector<double>());
		pSeqs.push_back(new vector<char>());
	}

}

void MS2ScanVector::postMvh() {
	for (int i = 0; i < num_max_threads; ++i) {
		delete _ppdAAforward.at(i);
		delete _ppdAAreverse.at(i);
		delete psequenceIonMasses.at(i);
		delete pSeqs.at(i);
	}
	_ppdAAforward.clear();
	_ppdAAreverse.clear();
	psequenceIonMasses.clear();
	pSeqs.clear();
}

void MS2ScanVector::preXcorr() {
	//MH: Must be equal to largest possible array
	int iArraySizeScore = (int) ((ProNovoConfig::dMaxMS2ScanMass * 2 + 100) * ProNovoConfig::dHighResInverseBinWidth);
	CometSearchMod::iArraySizeScore = iArraySizeScore;
	num_max_threads = omp_get_max_threads();
	for (int i = 0; i < num_max_threads; ++i) {
		vpbDuplFragmentGlobal.push_back(new bool[iArraySizeScore]());
		v_pdAAforwardGlobal.push_back(new double[MAX_PEPTIDE_LEN]());
		v_pdAAreverseGlobal.push_back(new double[MAX_PEPTIDE_LEN]());
		unsigned int *** _uiBinnedIonMasses = new unsigned int**[ProNovoConfig::iMaxPercusorCharge + 1]();
		for (int ii = 0; ii < ProNovoConfig::iMaxPercusorCharge + 1; ii++) {
			_uiBinnedIonMasses[ii] = new unsigned int*[9]();
			for (int j = 0; j < 9; j++) {
				_uiBinnedIonMasses[ii][j] = new unsigned int[MAX_PEPTIDE_LEN]();
			}
		}
		v_uiBinnedIonMassesGlobal.push_back(_uiBinnedIonMasses);
	}
	// SIP mode variable
	for (int i = 0; i < num_max_threads; ++i) {
		vvpbDuplFragmentGlobal.push_back(vector<bool>());
		vvpbDuplFragmentGlobal.back().clear();
		vvpbDuplFragmentGlobal.back().resize(iArraySizeScore, false);
		vvdBinnedIonMassesGlobal.push_back(vector<double>());
		vvdBinnedIonMassesGlobal.back().clear();
		vvdBinnedIonMassesGlobal.back().resize(iArraySizeScore, 0);
		vvdBinGlobal.push_back(vector<int>());
	}
}

void MS2ScanVector::postXcorr() {
	for (int i = 0; i < num_max_threads; ++i) {
		delete[] vpbDuplFragmentGlobal.at(i);
		delete[] v_pdAAforwardGlobal.at(i);
		delete[] v_pdAAreverseGlobal.at(i);
		for (int ii = 0; ii < ProNovoConfig::iMaxPercusorCharge + 1; ii++) {
			for (int j = 0; j < 9; j++) {
				delete[] v_uiBinnedIonMassesGlobal.at(i)[ii][j];
			}
			delete[] v_uiBinnedIonMassesGlobal.at(i)[ii];
		}
		delete[] v_uiBinnedIonMassesGlobal.at(i);
	}
	vpbDuplFragmentGlobal.clear();
	v_pdAAforwardGlobal.clear();
	v_pdAAreverseGlobal.clear();
	v_uiBinnedIonMassesGlobal.clear();
	// SIP mode variable
	vvpbDuplFragmentGlobal.clear();
	vvdBinnedIonMassesGlobal.clear();
	vvdBinGlobal.clear();
}

