
#include "ms2scanvector.h"
#include "memory.h"

MS2ScanVector::MS2ScanVector(const string &sFT2FilenameInput, const string &sOutputDirectory,
							 const string &sConfigFilename, bool bScreenOutput)
{
	unsigned int n;
	vector<string> vsSingleResidueNames = ProNovoConfig::vsSingleResidueNames;
	vector<double> vdSingleResidueMasses = ProNovoConfig::vdSingleResidueMasses;
	sFT2Filename = sFT2FilenameInput;
	sConfigFile = sConfigFilename;
	// mass_w = ProNovoConfig::getParentMassWindows();
	setOutputFile(sFT2FilenameInput, sOutputDirectory);
	this->bScreenOutput = bScreenOutput;
	for (n = 0; n < vsSingleResidueNames.size(); ++n)
		mapResidueMass[vsSingleResidueNames[n][0]] = vdSingleResidueMasses[n];
}

MS2ScanVector::~MS2ScanVector()
{
	// the destructors will free memory from vpAllMS2Scans
	// vector<MS2Scan *>::iterator it;
	// for (it = vpAllMS2Scans.begin(); it < vpAllMS2Scans.end(); it++)
	// 	delete *it;
	// cout << "MS2ScanVector" << endl;
	// vpAllMS2Scans.clear();
	// vpAllMS2ScanPtrs.clear();
	// vpPrecursorMasses.clear();
}

void MS2ScanVector::setOutputFile(const string &sFT2FilenameInput, const string &sOutputDirectory)
{
	string sLeftFileName, sMiddleFileName, sRightFileName = ".sip";
	size_t pos;
	// cout<<ProNovoConfig::getSearchName()<<endl;
	if ((ProNovoConfig::getSearchName() == "") || (ProNovoConfig::getSearchName() == "Null") || (ProNovoConfig::getSearchName() == "NULL") || (ProNovoConfig::getSearchName() == "null"))
		sRightFileName = ".sip";
	else
		sRightFileName = "." + ProNovoConfig::getSearchName() + ".sip";

	sLeftFileName = sFT2FilenameInput.substr(0, sFT2FilenameInput.length() - 4);
	if (sOutputDirectory == "")
		// if output directory is not specified, it is working directory by default
		sOutputFile = sLeftFileName + sRightFileName;
	else
	{
		pos = sLeftFileName.rfind(ProNovoConfig::getSeparator());
		if (pos == string::npos)
			sMiddleFileName = sLeftFileName;
		else
			sMiddleFileName = sLeftFileName.substr(pos + 1);
		sOutputFile = sOutputDirectory + ProNovoConfig::getSeparator() + sMiddleFileName + sRightFileName;
	}
}

bool MS2ScanVector::ReadFT2File()
{
	bool bReVal, flag_1stScan = true; // flag_1stScan true indicates pMS2Scan is empty
	string sline;
	istringstream input;
	ifstream ft2_stream(sFT2Filename.c_str());
	int tmp_charge;
	double tmp_mz, tmp_intensity;
	MS2Scan pMS2Scan;
	bReVal = ft2_stream.is_open();
	if (bReVal)
	{
		while (!ft2_stream.eof())
		{
			sline.clear();
			//	    string tmp_sline = sline;
			getline(ft2_stream, sline);
			// cout << sline <<"wyf"<<endl;
			if (sline == "")
				continue;
			// read vdMZ, vdIntensity, viCharge from chunk 0, 1, 5 at line started with number
			if ((sline.at(0) >= '0') && (sline.at(0) <= '9'))
			{
				TokenVector words(sline, " \t\n\r");
				if (words.size() < 6)
				// The judgement of MS scan resolution here will be overwritten,
				// we still keep the related codes here for future possible usuage
				{
					pMS2Scan.isMS2HighRes = false;
					pMS2Scan.viCharge.push_back(0);
				}
				else
				{
					pMS2Scan.isMS2HighRes = true;
					input.clear();
					input.str(words[5]);
					input >> tmp_charge;
					input.clear();
					pMS2Scan.viCharge.push_back(tmp_charge);
				}
				input.clear();
				input.str(words[0]);
				input >> tmp_mz;
				input.clear();
				pMS2Scan.vdMZ.push_back(tmp_mz);
				input.str(words[1]);
				input >> tmp_intensity;
				input.clear();
				pMS2Scan.vdIntensity.push_back(tmp_intensity);
			}
			// read iScanId, isMS1HighRes, dParentMZ from chunk 1, 3 at line started with S
			else if (sline.at(0) == 'S')
			{
				if (flag_1stScan)
					flag_1stScan = false;
				// ignore empty scan
				else if (pMS2Scan.vdIntensity.empty())
					;
				else
				{
					if (ProNovoConfig::getMassAccuracyFragmentIon() < 0.1)
						pMS2Scan.isMS2HighRes = true;
					else
						pMS2Scan.isMS2HighRes = false;
					saveScan(pMS2Scan);
				}
				// reset pMS2Scan after saving or ignoring empty scan
				pMS2Scan = MS2Scan();
				pMS2Scan.sFT2Filename = sFT2Filename;
				TokenVector words(sline, " \r\t\n");
				input.clear();
				input.str(words[1]);
				input >> pMS2Scan.iScanId;
				input.clear();
				input.str(words[2]);
				input >> pMS2Scan.dParentMZ;
				input.clear();
				pMS2Scan.isMS1HighRes = isMS1HighRes(words[2]);
				// in case no Z Line
				pMS2Scan.iParentChargeState = 0;
				pMS2Scan.dParentNeutralMass = 0;
			}
			// read iParentChargeState from chunk 1 at line started with Z
			else if (sline.at(0) == 'Z')
			{
				TokenVector words(sline, " \r\t\n");
				input.clear();
				input.str(words[1]);
				input >> pMS2Scan.iParentChargeState;
				input.clear();
				// input.str(words[2]);
				// input >> pMS2Scan->dParentNeutralMass;
				// input.clear();
			}
			// read ScanType from chunk 2 at line started with "I	ScanType"
			else if (sline.at(0) == 'I')
			{
				TokenVector words(sline, " \r\t\n");
				if (words[1] == "ScanType")
					pMS2Scan.setScanType(words[2]);
			}
		}
		ft2_stream.clear();
		ft2_stream.close();

		// recognition high or low MS2
		//	cout<<ProNovoConfig::getMassAccuracyFragmentIon()<<endl;

		// save the last scan
		if (ProNovoConfig::getMassAccuracyFragmentIon() < 0.1)
			pMS2Scan.isMS2HighRes = true;
		else
			pMS2Scan.isMS2HighRes = false;
		if (!flag_1stScan) // To avoid empty file
			saveScan(pMS2Scan);
	}
	return bReVal;
}

bool MS2ScanVector::loadFT2file()
{
#ifdef Ticktock
	TICK(loadFT2file);
#endif
	// read all MS2 scans from the file and populate vpAllMS2Scans
	// sort all MS2 scans in vpAllProteins by ascending order of their precursor masses
	// save their precursor masses in vpPrecursorMasses to quick look-up in assignPeptides2Scans()
	bool bReVal; // false when the file fails to be opened.
	// vector<MS2Scan>::iterator it;
	bReVal = ReadFT2File();
	if (bReVal)
	{
		sort(vpAllMS2Scans.begin(), vpAllMS2Scans.end(), myless);
		vpAllMS2ScanPtrs.reserve(vpAllMS2Scans.size());
		vpPrecursorMasses.reserve(vpAllMS2Scans.size());
		for (size_t i = 0; i < vpAllMS2Scans.size(); i++)
		{
			vpPrecursorMasses.push_back(vpAllMS2Scans[i].dParentNeutralMass);
			vpAllMS2ScanPtrs.push_back(vpAllMS2Scans.data() + i);
		}
	}
#ifdef Ticktock
	TOCK1ST(loadFT2file);
#endif
	return bReVal;
}

bool MS2ScanVector::isMS1HighRes(const std::string &target)
{
	bool bReVal = true;
	if (target.at(target.size() - 2) == '.')
		bReVal = false;
	return bReVal;
}

bool MS2ScanVector::ChargeDetermination(const std::vector<double> &vdAllmz, double pmz)
// return true, if the charge state is decided to be one

{
	bool bReVal = false;
	int iPosi;
	vector<double> vdtempMZ = vdAllmz;
	iPosi = max(0, (int)(((int)vdtempMZ.size()) * 0.05 - 1));
	nth_element(vdtempMZ.begin(), vdtempMZ.begin() + iPosi, vdtempMZ.end(), mygreater);
	if (vdtempMZ.at(iPosi) < pmz)
		bReVal = true;
	return bReVal;
}

bool MS2ScanVector::mygreater(double i, double j)
{
	return (i > j);
}

bool MS2ScanVector::myless(const MS2Scan &pMS2Scan1, const MS2Scan &pMS2Scan2)
{
	return (pMS2Scan1.dParentNeutralMass < pMS2Scan2.dParentNeutralMass);
}

void MS2ScanVector::saveScan(MS2Scan &pMS2Scan)
{ // parentChargeState > 0, save, otherwise, try 1, or 2 and 3
	int j;
	bool bchargeOne;
	if (pMS2Scan.iParentChargeState == 0)
	{
		bchargeOne = ChargeDetermination(pMS2Scan.vdMZ, pMS2Scan.dParentMZ);
		if (bchargeOne)
		{
			j = 1;
			pMS2Scan.iParentChargeState = j;
			pMS2Scan.dParentNeutralMass = pMS2Scan.dParentMZ * pMS2Scan.iParentChargeState - pMS2Scan.iParentChargeState * ProNovoConfig::getProtonMass();
			vpAllMS2Scans.push_back(move(pMS2Scan));
		}
		else
		// if the charge state is not +1, it could be greater than +3,
		// but current we just try +2 and +3.
		{
			MS2Scan pMS2newScan;
			for (j = 2; j <= 3; j++)
			{
				pMS2newScan = pMS2Scan;
				pMS2newScan.iParentChargeState = j;
				pMS2newScan.dParentNeutralMass = pMS2newScan.dParentMZ * pMS2newScan.iParentChargeState - pMS2newScan.iParentChargeState * ProNovoConfig::getProtonMass();
				vpAllMS2Scans.push_back(move(pMS2newScan));
			}
		}
	}
	else
	{
		pMS2Scan.dParentNeutralMass = pMS2Scan.dParentMZ * pMS2Scan.iParentChargeState - pMS2Scan.iParentChargeState * ProNovoConfig::getProtonMass();
		vpAllMS2Scans.push_back(move(pMS2Scan));
	}
}

void MS2ScanVector::preProcessAllMS2()
{
	int i, iScanSize;
	iScanSize = (int)vpAllMS2Scans.size();
	if (bScreenOutput)
		cout << "Preprocessing " << vpAllMS2Scans.size() << " scans " << endl;
#pragma omp parallel for schedule(guided)
	for (i = 0; i < iScanSize; i++)
		vpAllMS2ScanPtrs.at(i)->preprocess();
}

void MS2ScanVector::GetAllRangeFromMass(double dPeptideMass,
										vector<std::pair<int, int>> &vpPeptideMassRanges)
// all ranges of MS2 scans are stored in  vpPeptideMassWindows
{
	int i;
	pair<int, int> pairMS2Range;
	pair<int, int> lastPairRange(-100, -100);
	vector<pair<double, double>> vpPeptideMassWindows;
	vpPeptideMassWindows.clear();
	vpPeptideMassRanges.clear();
	ProNovoConfig::getPeptideMassWindows(dPeptideMass, vpPeptideMassWindows);
	for (i = 0; i < (int)vpPeptideMassWindows.size(); i++)
	{
		pairMS2Range = GetRangeFromMass(vpPeptideMassWindows.at(i).first,
										vpPeptideMassWindows.at(i).second);
		if ((pairMS2Range.first > -1) && (pairMS2Range.second > -1))
		{
			if ((lastPairRange.first < 0) || (lastPairRange.second < 0))
				lastPairRange = pairMS2Range;
			else
			{
				if (lastPairRange.second > pairMS2Range.first)
					lastPairRange.second = pairMS2Range.second;
				else
				{
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
// lb and ub are lower and upper bounds of acceptable parent mass values
{
	pair<int, int> p;
	int low = 0, high = vpPrecursorMasses.size() - 1, mid;
	double target;
	target = (lb + ub) / 2.0;
	// double lb, ub;  // lower and upper bounds on acceptable parent mass values
	// ub = target + error;
	// lb = target - error;
	while ((high - low) > 1)
	{
		mid = (high + low) / 2;
		if (vpPrecursorMasses[mid] > target)
			high = mid;
		else
			low = mid;
	}

	// Iterate till we get to the first element > than the lower bound
	int ndx = low;
	// cout<<scan_mass_list_.size()<<endl;
	if (vpPrecursorMasses[ndx] >= lb)
	{
		while (ndx >= 0 && vpPrecursorMasses[ndx] >= lb)
			ndx--;
		ndx++;
	}
	else
		while (ndx < (int)vpPrecursorMasses.size() && vpPrecursorMasses[ndx] < lb)
			ndx++;
	if (ndx == (int)vpPrecursorMasses.size() || vpPrecursorMasses[ndx] > ub)
		p = make_pair(-1, -1);
	else
	{
		low = ndx;
		while (ndx < (int)vpPrecursorMasses.size() && vpPrecursorMasses[ndx] <= ub)
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

void MS2ScanVector::assignPeptides2Scans(Peptide *currentPeptide)
{
	int i, j;
	vector<pair<int, int>> vpPeptideMassRanges;
	pair<int, int> pairMS2Range;

	GetAllRangeFromMass(currentPeptide->getPeptideMass(), vpPeptideMassRanges);

	for (j = 0; j < (int)vpPeptideMassRanges.size(); j++)
	{
		pairMS2Range = vpPeptideMassRanges.at(j);
		if ((pairMS2Range.first > -1) && (pairMS2Range.second > -1))
			for (i = pairMS2Range.first; i <= pairMS2Range.second; i++)
				vpAllMS2ScanPtrs.at(i)->vpPeptides.push_back(currentPeptide);
	}
}

void MS2ScanVector::processPeptideArray(vector<Peptide *> &vpPeptideArray)
{
	int i, iPeptideArraySize, iScanSize;
	iPeptideArraySize = (int)vpPeptideArray.size();

//    for (int i=0; i< (int) vpPeptideArray.size(); i++)
//      cout<<vpPeptideArray.at(i)->getPeptideSeq() <<"\t"
//	  <<vpPeptideArray.at(i)->getOriginalPeptideSeq()<<"\t"
//	  <<vpPeptideArray.at(i)->getProteinName()<<"\t"
//	  <<vpPeptideArray.at(i)->getBeginPosProtein()<<"\t"
//	  <<vpPeptideArray.at(i)->getPeptideMass()<<endl;

// cout<<"calculating fragments of "<<vpPeptideArray.size()<<"  peptides"<<endl;
#ifdef Ticktock
	TICK(preprocessing);
#endif
#pragma omp parallel for shared(vpPeptideArray) private(i) \
	schedule(guided)

	for (i = 0; i < iPeptideArraySize; i++)
		vpPeptideArray[i]->preprocessing(vpAllMS2ScanPtrs.at(0)->isMS2HighRes, mapResidueMass);
// vpPeptideArray[i]->calculateExpectedFragments(mapResidueMass);

// cout<<"scoring "<<vpPeptideArray.size()<<"  peptides"<<endl;
//  every MS2 scans scores their matched peptides
#ifdef Ticktock
	TOCK(preprocessing);
	TICK(scorePeptides);
#endif
	iScanSize = (int)vpAllMS2Scans.size();
#pragma omp parallel for schedule(guided)

	for (i = 0; i < iScanSize; i++)
		vpAllMS2ScanPtrs[i]->scorePeptides();

	// free memory of all peptide objects's content
	// but not delete the peptide object
	for (i = 0; i < (int)vpPeptideArray.size(); i++)
		vpPeptideArray[i]->~Peptide();
#ifdef Ticktock
	TOCK(scorePeptides);
#endif
	// empty peptide array
	vpPeptideArray.clear();
}

void MS2ScanVector::searchDatabase()
{
	ProteinDatabase myProteinDatabase(bScreenOutput);
	myProteinDatabase.loadDatabase();
	if (myProteinDatabase.getFirstProtein())
	{
#ifdef Ticktock
		TICK(loadDatabase);
#endif
		vector<Peptide *> vpPeptideArray;
		vpPeptideArray.reserve(PEPTIDE_ARRAY_SIZE);
		char *peptidesBuffer = new char[PEPTIDE_ARRAY_SIZE * sizeof(Peptide)];
		// memset(peptidesBuffer, 0, sizeof(Peptide) * PEPTIDE_ARRAY_SIZE);
		Peptide *firstPeptidePtr = (Peptide *)peptidesBuffer;
		Peptide *currentPeptidePtr = new (firstPeptidePtr) Peptide();
		int peptidesCount = 0;
		// get one peptide from the database at a time, until there is no more peptide
		while (myProteinDatabase.getNextPeptide(currentPeptidePtr))
		{
			// save the new peptide to the array
			vpPeptideArray.push_back(currentPeptidePtr);
			// assign the pointers of peptides to appropriete MS2Scans
			assignPeptides2Scans(currentPeptidePtr);
			// create a new peptide for the next iteration
			peptidesCount++;
			// when the vpPeptideArray is full
			if (vpPeptideArray.size() >= PEPTIDE_ARRAY_SIZE)
			{
#ifdef Ticktock
				TOCK(loadDatabase);
#endif
				processPeptideArray(vpPeptideArray);
#ifdef Ticktock
				TICKAGAIN(loadDatabase);
#endif
				currentPeptidePtr = new (firstPeptidePtr) Peptide();
				peptidesCount = 0;
			}
			else
				currentPeptidePtr = new (firstPeptidePtr + peptidesCount) Peptide();
		}
		// the last peptide object is an empty object and need to be deleted
		// delete currentPeptidePtr;
		// there are still unprocessed peptides in the vpPeptideArray
		// need to process them in the same manner
		// the following code is the same as inside if(vpPeptideArray.size() >= PEPTIDE_ARRAY_SIZE )
		//    cout<<vpPeptideArray.size()<<endl;
		if (!vpPeptideArray.empty())
		{
#ifdef Ticktock
			TOCK(loadDatabase);
#endif
			processPeptideArray(vpPeptideArray);
		}
		// delete the peptide objects
		delete[] peptidesBuffer;
	}
}

void MS2ScanVector::startProcessing()
{
// Preprocessing all MS2 scans by mult-threading
#ifdef Ticktock
	TICK(preProcessAllMS2);
#endif
	preProcessAllMS2();
#ifdef Ticktock
	TOCK(preProcessAllMS2);
#endif

// Search all MS2 scans against the database by mult-threading
#ifdef Ticktock
	TICK(searchDatabase);
#endif
	searchDatabase();
#ifdef Ticktock
	TOCK(searchDatabase);
#endif

	// Postprocessing all MS2 scans' results by mult-threading
	postProcessAllMS2();

// write results to a SIP file
#ifdef Ticktock
	TICK(writeOutput);
#endif
	writeOutput();
#ifdef Ticktock
	TOCK(writeOutput);
#endif
}

void MS2ScanVector::postProcessAllMS2()
{
	int i;
	for (i = 0; i < (int)vpAllMS2Scans.size(); i++)
		vpAllMS2ScanPtrs.at(i)->postprocess();
}

bool MS2ScanVector::mylessScanId(MS2Scan *pMS2Scan1, MS2Scan *pMS2Scan2)
{
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

void MS2ScanVector::writeOutput()
{
	int i, j;
	ofstream outputFile;
	ifstream configStream;
	string sline, sTailFT2FileName;
	double dcurrentMeanWeightSum, dcurrentDeviationWeightSum;
	sTailFT2FileName = ParsePath(sFT2Filename);

	// can only sort vpAllMS2ScanPtrs; sort vpAllMS2Scan will corrupt
	sort(vpAllMS2ScanPtrs.begin(), vpAllMS2ScanPtrs.end(), mylessScanId);

	outputFile.open(sOutputFile.c_str());
	configStream.open(sConfigFile.c_str());
	while (!configStream.eof())
	{
		sline.clear();
		getline(configStream, sline);
		outputFile << "#\t" << sline << endl;
	}

	outputFile << "Filename\tScanNumber\tParentCharge\tMeasuredParentMass\t";
	outputFile << "CalculatedParentMass\tScanType\tSearchName\tScoringFunction\tRank\tScore\t";
	outputFile << "IdentifiedPeptide\tOriginalPeptide\tProteinNames" << endl;
	for (i = 0; i < vpAllMS2ScanPtrs.size(); i++)
		if (!vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.empty())
		{
			if (vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.at(0)->sScoringFunction == "WeightSum")
				calculateMeanAndDeviation(vpAllMS2ScanPtrs.at(i)->inumberofWeightSumScore,
										  vpAllMS2ScanPtrs.at(i)->dsumofWeightScore,
										  vpAllMS2ScanPtrs.at(i)->dsumofSquareWeightSumScore,
										  dcurrentMeanWeightSum, dcurrentDeviationWeightSum);
			for (j = 0; j < (int)vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.size(); j++)
			{
				/*		cout<<sFT2Filename<<"\t"<<vpAllMS2Scans.at(i)->iScanId<<"\t";
		cout<<vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->getPeptideScore()<<"\t";
		cout<<vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->getPeptideSeq()<<"\t";
		cout<<vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->getOriginalPeptideSeq()<<"\t";
		cout<<vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->getProteinName()<<"\t";
		cout<<vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->getBeginPosProtein()<<"\t";
		cout<<vpAllMS2Scans.at(i)->iParentChargeState<<"\t";
		cout<<i<<"\t";
		cout<<vpAllMS2Scans.at(i)->vdWeightSumAllScores.size()<<"\t0\thit"<<endl;
		*/
				outputFile << sTailFT2FileName << "\t";
				outputFile << vpAllMS2ScanPtrs.at(i)->iScanId << "\t";
				outputFile << vpAllMS2ScanPtrs.at(i)->iParentChargeState << "\t";
				outputFile << setiosflags(ios::fixed) << setprecision(5) << vpAllMS2ScanPtrs.at(i)->dParentNeutralMass << "\t";
				outputFile << setiosflags(ios::fixed) << setprecision(5) << vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.at(j)->dCalculatedParentMass << "\t";

				outputFile << vpAllMS2ScanPtrs.at(i)->getScanType() << "\t";
				outputFile << ProNovoConfig::getSearchName() << "\t";

				outputFile << vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.at(j)->sScoringFunction << "\t";
				outputFile << (j + 1) << "\t";
				outputFile << setiosflags(ios::fixed) << setprecision(2) << vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.at(j)->dScore << "\t";
				// outputFile<<((vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->dScore-dcurrentMeanWeightSum)/dcurrentDeviationWeightSum)<<"\t";

				if (vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifyPrefix != '-')
					outputFile << vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifyPrefix;
				outputFile << vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.at(j)->sIdentifiedPeptide;
				if (vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifySuffix != '-')
					outputFile << vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifySuffix;
				outputFile << "\t";

				if (vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalPrefix != '-')
					outputFile << vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalPrefix;
				outputFile << vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.at(j)->sOriginalPeptide;
				if (vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalSuffix != '-')
					outputFile << vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalSuffix;
				outputFile << "\t";

				outputFile << "{" << vpAllMS2ScanPtrs.at(i)->vpWeightSumTopPeptides.at(j)->sProteinNames << "}" << endl;
			}
		}
	configStream.clear();
	configStream.close();
	outputFile.close();
}

void MS2ScanVector::calculateMeanAndDeviation(int inumberScore, double dScoreSum,
											  double dScoreSquareSum, double &dMean, double &dDeviation)
{
	dMean = dScoreSum / inumberScore;
	// dDeviation is sample standard deviation

	if (inumberScore > 1)
		dDeviation = sqrt((dScoreSquareSum - dMean * dMean) * inumberScore / (inumberScore - 1));

	if ((inumberScore == 1) || (dDeviation < ZERO))
		dDeviation = ZERO;
}
