#include "ms2scan.h"

/**********switching from weightsum to ranksum needs to turn on some function in ms2scan.cpp and peptide.cpp *******/

MS2Scan::MS2Scan() {
	dMassTolerance = ProNovoConfig::getMassAccuracyFragmentIon();
	dProtonMass = ProNovoConfig::getProtonMass();
	inumberofWeightSumScore = 0;
	dsumofSquareWeightSumScore = 0;
	dsumofWeightScore = 0;
	peakData = NULL;
	pPeakList = NULL;
	intenClassCounts = NULL;
	sRTime = "0.0";
	totalPeakBins = 0;
	pQuery = NULL;
	mzUpperBound = 0;
	dParentMass = 0;
	dMaxIntensity = 0;
	mzLowerBound = 0;
	bSkip = false;
	iNumPeptideAssigned = 0;
	dSumIntensity = 0;
	iScanId = 0;
	iParentChargeState = 0;
	bSetMS2Flag = false;
	isMS1HighRes = true;
	dParentNeutralMass = 0;
	dParentMZ = 0;
	isMS2HighRes = true;
	iMaxMZ = 0;
	bin_res = 0;
	iMinMZ = 0;
	sScanType = "null";
}

MS2Scan::~MS2Scan() {
	int i;
	for (i = 0; i < (int) vpWeightSumTopPeptides.size(); i++)
		delete vpWeightSumTopPeptides.at(i);

	if (peakData != NULL) {
		delete peakData;
		peakData = NULL;
	}

	if (pPeakList != NULL) {
		delete pPeakList;
		pPeakList = NULL;
	}

	if (intenClassCounts != NULL) {
		delete intenClassCounts;
		intenClassCounts = NULL;
	}

	if (pQuery != NULL) {
		delete pQuery;
		pQuery = NULL;
	}
}

void MS2Scan::postprocess() {
	// WDP score

}

bool MS2Scan::mergePeptide(vector<PeptideUnit*>& vpTopPeptides, const string & sPeptide, const string & sProteinName) {
	int i;
	bool bReVal = false, bNewProtein;
	size_t iPosition;
	if (!vpTopPeptides.empty())
		for (i = 0; i < (int) vpTopPeptides.size(); i++)
			if (vpTopPeptides.at(i)->sIdentifiedPeptide == sPeptide) {
				bNewProtein = true; // True: sProteinName doesn't appear in vpTopPeptides.at(i)->sProteinNames
				if (sProteinName == vpTopPeptides.at(i)->sProteinNames)
					bNewProtein = false;
				else {
					iPosition = vpTopPeptides.at(i)->sProteinNames.find("," + sProteinName);
					if (iPosition != string::npos)
						bNewProtein = false;
				}
				if (bNewProtein)
					vpTopPeptides.at(i)->sProteinNames += "," + sProteinName;
				bReVal = true;
				break;
			}
	return bReVal;
}

void MS2Scan::scorePeptidesHighMS2() {
	int i;
	for (i = 0; i < (int) vpPeptides.size(); i++) {

		if (!mergePeptide(vpWeightSumTopPeptides, vpPeptides.at(i)->getPeptideSeq(), vpPeptides.at(i)->getProteinName()))
			scoreWeightSumHighMS2(vpPeptides.at(i));
		//scoreRankSumHighMS2(vpPeptides.at(i));
	}
}

void MS2Scan::scorePeptidesLowMS2() {
	int i;
	for (i = 0; i < (int) vpPeptides.size(); i++) {
		if (!mergePeptide(vpWeightSumTopPeptides, vpPeptides.at(i)->getPeptideSeq(), vpPeptides.at(i)->getProteinName())) {
			scoreWeightSum(vpPeptides.at(i));
//	    scoreRankSum(vpPeptides.at(i));
		}
	}
}

void MS2Scan::cleanup() {
//    cout<<"vdMZ: "<<vdMZ.size()<<endl;
//    cout<<"bin "<<vbPeakPresenceBins.size()<<endl;
	vector<double>().swap(vdMaxMzIntensity);
	vector<double>().swap(vdMzIntensity);
	vector<double>().swap(vdHighIntensity);
	// vector<double>().swap(vdMZ);
	// vector<double>().swap(vdIntensity);
	// vector<int>().swap(viCharge);
	// vdMaxMzIntensity.clear();
	//   vdMzIntensity.clear();
	//   vdHighIntensity.clear();
	//   vdMZ.clear();
	//   vdIntensity.clear();
	//   viCharge.clear();
}

void MS2Scan::scoreWeightSum(Peptide* currentPeptide) {
	double dScore = 0;
	// calculate score
	vector<bool> vbFragmentZ2; //true: y side gets the +2; false: b ion gets +2, only calculated for +3 peptides
	int iNumFragments = 0, n, iIndex4MostAbundunt = 0;
	double dBweight = 0, dYweight = 0, dExpectedMZ, dExYweight = 0, dExBweight = 0;
	int iExYreward = 0, iExBreward = 0, iCurrentYreward, iCurrentBreward;

	if (iParentChargeState >= 3)
		WeightCompare(currentPeptide->getPeptideSeq(), vbFragmentZ2);
	iNumFragments = currentPeptide->vdYionMasses.size();
	for (n = 0; n < iNumFragments; ++n) {
		dYweight = 0;
		dBweight = 0;
		dExpectedMZ = currentPeptide->vdYionMasses[n] + ProNovoConfig::getProtonMass();
		if (searchMZ2D(dExpectedMZ, iIndex4MostAbundunt))
			dYweight = ProNovoConfig::scoreError(fabs(vdpreprocessedMZ[iIndex4MostAbundunt] - dExpectedMZ)) + vdpreprocessedIntensity[iIndex4MostAbundunt];
		dExpectedMZ = currentPeptide->vdBionMasses[iNumFragments - n - 1] + ProNovoConfig::getProtonMass();
		if (searchMZ2D(dExpectedMZ, iIndex4MostAbundunt))
			dBweight = ProNovoConfig::scoreError(fabs(vdpreprocessedMZ[iIndex4MostAbundunt] - dExpectedMZ)) + vdpreprocessedIntensity[iIndex4MostAbundunt];
		if (iParentChargeState >= 3) {
			if (vbFragmentZ2[n]) {
				dExpectedMZ = (currentPeptide->vdYionMasses[n] + ProNovoConfig::getProtonMass() * 2) / 2;
				if (searchMZ2D(dExpectedMZ, iIndex4MostAbundunt))
					dYweight += ProNovoConfig::scoreError(fabs(vdpreprocessedMZ[iIndex4MostAbundunt] - dExpectedMZ))
							+ vdpreprocessedIntensity[iIndex4MostAbundunt];
			} else {
				dExpectedMZ = (currentPeptide->vdBionMasses[iNumFragments - n - 1] + ProNovoConfig::getProtonMass() * 2) / 2;
				if (searchMZ2D(dExpectedMZ, iIndex4MostAbundunt))
					dBweight += ProNovoConfig::scoreError(fabs(vdpreprocessedMZ[iIndex4MostAbundunt] - dExpectedMZ))
							+ vdpreprocessedIntensity[iIndex4MostAbundunt];
			}
		}

		//bouns
		if (dYweight > ZERO && dBweight > ZERO) {
			dScore = dScore + (dYweight + dBweight) * 2;
			if (dExYweight > ZERO) {
				dScore += dYweight;
				iCurrentYreward = 3;
				if (iExYreward < 3)
					dScore += dExYweight;
			} else
				iCurrentYreward = 2;
			if (dExBweight > ZERO) {
				dScore += dBweight;
				iCurrentBreward = 3;
				if (iExBreward < 3)
					dScore += dExBweight;
			} else
				iCurrentBreward = 2;
		} else {
			dScore = dScore + (dYweight + dBweight);
			if (dExYweight > ZERO) {
				dScore += dYweight;
				iCurrentYreward = 2;
				if (iExYreward < 3)
					dScore += dExYweight;
			} else
				iCurrentYreward = 0;
			if (dExBweight > ZERO) {
				dScore += dBweight;
				iCurrentBreward = 2;
				if (iExBreward < 3)
					dScore += dExBweight;
			} else
				iCurrentBreward = 0;
		}
		dExYweight = dYweight;
		dExBweight = dBweight;
		iExYreward = iCurrentYreward;
		iExBreward = iCurrentBreward;
	}
	// save score
	saveScoreSIP(dScore, currentPeptide, vpWeightSumTopPeptides);
}

/*
 void MS2Scan::scoreRankSum(Peptide* currentPeptide)
 {
 double dScore = 0;
 // calculate score
 vector<bool> vbFragmentZ2;//true: y side gets the +2; false: b ion gets +2, only calculated for +3 peptides
 int iNumFragments = 0, iIndex4MostAbundunt=0;
 double dExpectedMZ;
 int iMatchCount = 0, iEmptyMZ;
 int iExMZ, iPeakNum, iUpperBound , iMZSize, ioneSum=0, iallSum=0, iUnOber=0, n_1, n_2;
 int i;
 vector <ScanUnit> vSAllUnits;
 double dUvalue, dDvalue, dMvalue;

 iNumFragments = currentPeptide->vdYionMasses.size();

 iPeakNum = (int)vdpreprocessedMZ.size();
 iUpperBound = iPeakNum -1;
 iExMZ = -100;

 for (i = 0; i<iPeakNum; i++)
 {
 if (((int) vdpreprocessedMZ[i]) > iMaxMZ)
 {
 iUpperBound = i - 1;
 break;
 }
 ScanUnit su(vdpreprocessedIntensity[i], false);
 vSAllUnits.push_back(su);
 if (((int) vdpreprocessedMZ[i]) == iExMZ)
 continue;
 else
 {
 iMatchCount++;
 iExMZ = (int) vdpreprocessedMZ[i];
 }
 }
 iEmptyMZ = iMaxMZ - iMinMZ + 1 - iMatchCount;
 iMZSize = (int) vSAllUnits.size();
 dScore = 0;
 if( iParentChargeState >= 3 )
 {
 WeightCompare(currentPeptide->getPeptideSeq(), vbFragmentZ2);
 for (i = 0; i < iNumFragments; i++)
 {
 if (vbFragmentZ2[i])
 dExpectedMZ = (currentPeptide->vdYionMasses[i]+ ProNovoConfig::getProtonMass()*2)/2;
 else
 dExpectedMZ = currentPeptide->vdYionMasses[i] + ProNovoConfig::getProtonMass();
 if (searchMZ2D( dExpectedMZ, iIndex4MostAbundunt))
 if (iIndex4MostAbundunt <= iUpperBound)
 vSAllUnits[iIndex4MostAbundunt].match = true;
 else
 iUnOber++;
 else
 iUnOber++;
 if (!(vbFragmentZ2[i]))
 dExpectedMZ = (currentPeptide->vdBionMasses[iNumFragments - i - 1] + ProNovoConfig::getProtonMass()*2)/2;
 else
 dExpectedMZ = currentPeptide->vdBionMasses[iNumFragments - i - 1] + ProNovoConfig::getProtonMass();
 if (searchMZ2D( dExpectedMZ, iIndex4MostAbundunt))
 if (iIndex4MostAbundunt <= iUpperBound)
 vSAllUnits[iIndex4MostAbundunt].match = true;
 else
 iUnOber++;
 else
 iUnOber++;
 }
 }else
 for (i = 0; i < iNumFragments; i++)
 {
 dExpectedMZ = currentPeptide->vdYionMasses[i] + ProNovoConfig::getProtonMass();
 if (searchMZ2D( dExpectedMZ, iIndex4MostAbundunt))
 if (iIndex4MostAbundunt <= iUpperBound)
 vSAllUnits[iIndex4MostAbundunt].match = true;
 else
 iUnOber++;
 else
 iUnOber++;
 dExpectedMZ = currentPeptide->vdBionMasses[iNumFragments - i - 1] + ProNovoConfig::getProtonMass();
 if (searchMZ2D( dExpectedMZ, iIndex4MostAbundunt))
 if (iIndex4MostAbundunt <= iUpperBound)
 vSAllUnits[iIndex4MostAbundunt].match = true;
 else
 iUnOber++;
 else
 iUnOber++;
 }
 sort(vSAllUnits.begin(), vSAllUnits.end(), mySUGreater);
 for (i=0; i<iMZSize; i++)
 if (vSAllUnits[i].match)
 iallSum += ioneSum;
 else
 ioneSum++;
 n_1 = iNumFragments*2;
 n_2 = iMZSize - (iNumFragments*2-iUnOber) + iEmptyMZ;
 dUvalue = iallSum + (iMZSize - (iNumFragments*2-iUnOber) + 0.5*iEmptyMZ  )*iUnOber ;
 dMvalue = n_1*n_2/2.0;
 dDvalue = sqrt(n_1*n_2*(n_1+n_2+1)/12.0);
 dScore = (dMvalue-dUvalue)/dDvalue;

 saveScore(dScore, currentPeptide, vpWeightSumTopPeptides, vdWeightSumAllScores, "RankSum");
 }*/

void MS2Scan::scoreRankSumHighMS2(Peptide* currentPeptide) {

	double dScore = 0;

	vector<bool> vbFragmentZ2; //true: y side gets the +2; false: b ion gets +2, only calculated for +3 peptides
	int iNumFragments = 0, iIndex4MostAbundunt = 0;
	double dExpectedMZ;
	int iMatchCount = 0, iEmptyMZ;
	int iLastMZ, iPeakNum, iUpperBound, iMZSize, iCurrentUnMatchPeakNumber = 0, iUnMatchPeakSum = 0, iUnOberserve = 0, iIonNumber, iUnMatchPeakNumber;
	int i;
	//vector <ScanUnit> vSAllUnits;
	vector<bool> vbMatch;
	double dUvalue;

	iNumFragments = currentPeptide->vdYionMasses.size();

	iPeakNum = (int) vdpreprocessedMZ.size();
	iUpperBound = iPeakNum - 1;
	iLastMZ = -100;

	for (i = 0; i < iPeakNum; i++) {
		if (((int) vdpreprocessedMZ[i]) > iMaxMZ) {
			iUpperBound = i - 1;
			break;
		}
		vbMatch.push_back(false);
		if (((int) vdpreprocessedMZ[i]) == iLastMZ)
			continue;
		else {
			iMatchCount++;
			iLastMZ = (int) vdpreprocessedMZ[i];
		}
	}
	iEmptyMZ = iMaxMZ - iMinMZ + 1 - iMatchCount;
	iMZSize = (int) vbMatch.size();
	dScore = 0;
	if (iParentChargeState >= 3) {
		WeightCompare(currentPeptide->getPeptideSeq(), vbFragmentZ2);
		for (i = 0; i < iNumFragments; i++) {
			if (vbFragmentZ2[i])
				dExpectedMZ = (currentPeptide->vdYionMasses[i] + ProNovoConfig::getProtonMass() * 2) / 2;
			else
				dExpectedMZ = currentPeptide->vdYionMasses[i] + ProNovoConfig::getProtonMass();
			if (searchMZ2D(dExpectedMZ, iIndex4MostAbundunt))
				if (iIndex4MostAbundunt <= iUpperBound)
					vbMatch[iIndex4MostAbundunt] = true;
				else
					iUnOberserve++;
			else
				iUnOberserve++;
			if (!(vbFragmentZ2[i]))
				dExpectedMZ = (currentPeptide->vdBionMasses[iNumFragments - i - 1] + ProNovoConfig::getProtonMass() * 2) / 2;
			else
				dExpectedMZ = currentPeptide->vdBionMasses[iNumFragments - i - 1] + ProNovoConfig::getProtonMass();
			if (searchMZ2D(dExpectedMZ, iIndex4MostAbundunt))
				if (iIndex4MostAbundunt <= iUpperBound)
					vbMatch[iIndex4MostAbundunt] = true;
				else
					iUnOberserve++;
			else
				iUnOberserve++;
		}
	} else
		for (i = 0; i < iNumFragments; i++) {
			dExpectedMZ = currentPeptide->vdYionMasses[i] + ProNovoConfig::getProtonMass();
			if (searchMZ2D(dExpectedMZ, iIndex4MostAbundunt))
				if (iIndex4MostAbundunt <= iUpperBound)
					vbMatch[iIndex4MostAbundunt] = true;
				else
					iUnOberserve++;
			else
				iUnOberserve++;
			dExpectedMZ = currentPeptide->vdBionMasses[iNumFragments - i - 1] + ProNovoConfig::getProtonMass();
			if (searchMZ2D(dExpectedMZ, iIndex4MostAbundunt))
				if (iIndex4MostAbundunt <= iUpperBound)
					vbMatch[iIndex4MostAbundunt] = true;
				else
					iUnOberserve++;
			else
				iUnOberserve++;
		}
	for (i = 0; i < iMZSize; i++)
		if (vbMatch[viIntensityRank[i]])
			iUnMatchPeakSum += iCurrentUnMatchPeakNumber;
		else
			iCurrentUnMatchPeakNumber++;
	iIonNumber = iNumFragments * 2;
	iUnMatchPeakNumber = iMZSize - (iNumFragments * 2 - iUnOberserve) + iEmptyMZ;
	dUvalue = iUnMatchPeakSum + (iMZSize - (iNumFragments * 2 - iUnOberserve) + 0.5 * iEmptyMZ) * iUnOberserve;
	if (iMZSize > 0)
		dScore = CalculateRankSum(dUvalue, iIonNumber, iUnMatchPeakNumber);
	saveScore(dScore, currentPeptide, vpWeightSumTopPeptides, "RankSum");
}

void MS2Scan::scoreRankSum(Peptide* currentPeptide) {

	double dScore = 0;

	vector<bool> vbFragmentZ2; //true: y side gets the +2; false: b ion gets +2, only calculated for +3 peptides
	int iNumFragments = 0, iIndex4MostAbundunt = 0;
	double dExpectedMZ;
	int iMatchCount = 0, iEmptyMZ;
	int iLastMZ, iPeakNum, iUpperBound, iMZSize, iCurrentUnMatchPeakNumber = 0, iUnMatchPeakSum = 0, iUnOberserve = 0, iIonNumber, iUnMatchPeakNumber;
	int i;
	//vector <ScanUnit> vSAllUnits;
	vector<bool> vbMatch;
	double dUvalue;

	iNumFragments = currentPeptide->vdYionMasses.size();

	iPeakNum = (int) vdpreprocessedMZ.size();
	iUpperBound = iPeakNum - 1;
	iLastMZ = -100;

	for (i = 0; i < iPeakNum; i++) {
		if (((int) vdpreprocessedMZ[i]) > iMaxMZ) {
			iUpperBound = i - 1;
			break;
		}
		vbMatch.push_back(false);
		if (((int) vdpreprocessedMZ[i]) == iLastMZ)
			continue;
		else {
			iMatchCount++;
			iLastMZ = (int) vdpreprocessedMZ[i];
		}
	}
	iEmptyMZ = iMaxMZ - iMinMZ + 1 - iMatchCount;
	iMZSize = (int) vbMatch.size();
	dScore = 0;
	if (iParentChargeState >= 3) {
		WeightCompare(currentPeptide->getPeptideSeq(), vbFragmentZ2);
		for (i = 0; i < iNumFragments; i++) {
			if (vbFragmentZ2[i])
				dExpectedMZ = (currentPeptide->vdYionMasses[i] + ProNovoConfig::getProtonMass() * 2) / 2;
			else
				dExpectedMZ = currentPeptide->vdYionMasses[i] + ProNovoConfig::getProtonMass();
			if (searchMZ2D(dExpectedMZ, iIndex4MostAbundunt))
				if (iIndex4MostAbundunt <= iUpperBound)
					vbMatch[iIndex4MostAbundunt] = true;
				else
					iUnOberserve++;
			else
				iUnOberserve++;
			if (!(vbFragmentZ2[i]))
				dExpectedMZ = (currentPeptide->vdBionMasses[iNumFragments - i - 1] + ProNovoConfig::getProtonMass() * 2) / 2;
			else
				dExpectedMZ = currentPeptide->vdBionMasses[iNumFragments - i - 1] + ProNovoConfig::getProtonMass();
			if (searchMZ2D(dExpectedMZ, iIndex4MostAbundunt))
				if (iIndex4MostAbundunt <= iUpperBound)
					vbMatch[iIndex4MostAbundunt] = true;
				else
					iUnOberserve++;
			else
				iUnOberserve++;
		}
	} else
		for (i = 0; i < iNumFragments; i++) {
			dExpectedMZ = currentPeptide->vdYionMasses[i] + ProNovoConfig::getProtonMass();
			if (searchMZ2D(dExpectedMZ, iIndex4MostAbundunt))
				if (iIndex4MostAbundunt <= iUpperBound)
					vbMatch[iIndex4MostAbundunt] = true;
				else
					iUnOberserve++;
			else
				iUnOberserve++;
			dExpectedMZ = currentPeptide->vdBionMasses[iNumFragments - i - 1] + ProNovoConfig::getProtonMass();
			if (searchMZ2D(dExpectedMZ, iIndex4MostAbundunt))
				if (iIndex4MostAbundunt <= iUpperBound)
					vbMatch[iIndex4MostAbundunt] = true;
				else
					iUnOberserve++;
			else
				iUnOberserve++;
		}
	for (i = 0; i < iMZSize; i++)
		if (vbMatch[viIntensityRank[i]])
			iUnMatchPeakSum += iCurrentUnMatchPeakNumber;
		else
			iCurrentUnMatchPeakNumber++;
	iIonNumber = iNumFragments * 2;
	iUnMatchPeakNumber = iMZSize - (iNumFragments * 2 - iUnOberserve) + iEmptyMZ;
	dUvalue = iUnMatchPeakSum + (iMZSize - (iNumFragments * 2 - iUnOberserve) + 0.5 * iEmptyMZ) * iUnOberserve;
	if (iMZSize > 0)
		dScore = CalculateRankSum(dUvalue, iIonNumber, iUnMatchPeakNumber);
	saveScore(dScore, currentPeptide, vpWeightSumTopPeptides, "RankSum");
}

/*
 bool MS2Scan::mySUGreater(ScanUnit s1, ScanUnit s2)
 {
 return ((s1.intensity) > (s2.intensity));
 }*/

bool MS2Scan::searchMZ(const double & dTarget, int & iIndex4Found) {
	double dErrRange = ProNovoConfig::getMassAccuracyFragmentIon();

	if (dTarget > (vdpreprocessedMZ.back() + dErrRange))
		return false;
	iIndex4Found = vbPeakPresenceBins[(unsigned long int) (dTarget * bin_res + SMALLINCREMENT)];
//	if (iIndex4Found != -1)
//	    cout<<dTarget<<" "<<vdpreprocessedMZ[iIndex4Found]<<" "<<fabs(vdpreprocessedMZ[iIndex4Found]-dTarget)<<endl;
	if (iIndex4Found == -1)
		return false;
	else
		return true;
}

bool MS2Scan::searchMZ2D(const double& dTarget, int& iIndex4Found) {
	bool bReVal = true;
	double dErrRange = ProNovoConfig::getMassAccuracyFragmentIon();
	double dCurrentIntensity = -1;
	int i, iLowerBound, iUpperBound;
	iIndex4Found = -1;
	if (dTarget > (vdpreprocessedMZ.back() + dErrRange)) {
		bReVal = false;
	} else {
		iLowerBound = vbPeakPresenceBins2D[(unsigned long int) (dTarget + SMALLINCREMENT)].first;
		if (iLowerBound > -1) {
			iUpperBound = vbPeakPresenceBins2D[(unsigned long int) (dTarget + SMALLINCREMENT)].second;
			for (i = iLowerBound; i <= iUpperBound; i++) {
				if ((vdpreprocessedMZ.at(i) >= (dTarget - dErrRange)) && (vdpreprocessedMZ.at(i) <= (dTarget + dErrRange)))
					if (vdpreprocessedIntensity.at(i) > dCurrentIntensity) {
						iIndex4Found = i;
						dCurrentIntensity = vdpreprocessedIntensity.at(i);
					}
			}
		}
	}
	if (iIndex4Found == -1)
		bReVal = false;
//    if (iIndex4Found > -1)
//	cout<<dTarget<<" "<<vdpreprocessedMZ[iIndex4Found]<<" "<<fabs(vdpreprocessedMZ[iIndex4Found]-dTarget)<<endl;
	return bReVal;
}

bool MS2Scan::searchMZ2D(const double& dTarget, const double & dErrRange, int& iIndex4Found) {
	bool bReVal = true;
	double dCurrentIntensity = -1;
	int i, iLowerBound, iUpperBound;
	iIndex4Found = -1;
	if (dTarget > (vdpreprocessedMZ.back() + dErrRange)) {
		bReVal = false;
	} else {
		iLowerBound = vbPeakPresenceBins2D[(unsigned long int) (dTarget + SMALLINCREMENT)].first;
		if (iLowerBound > -1) {
			iUpperBound = vbPeakPresenceBins2D[(unsigned long int) (dTarget + SMALLINCREMENT)].second;
			for (i = iLowerBound; i <= iUpperBound; i++) {
				if ((vdpreprocessedMZ.at(i) >= (dTarget - dErrRange)) && (vdpreprocessedMZ.at(i) <= (dTarget + dErrRange)))
					if (vdpreprocessedIntensity.at(i) > dCurrentIntensity) {
						iIndex4Found = i;
						dCurrentIntensity = vdpreprocessedIntensity.at(i);
					}
			}
		}
	}
	if (iIndex4Found == -1)
		bReVal = false;
//    if (iIndex4Found > -1)
//	cout<<dTarget<<" "<<vdpreprocessedMZ[iIndex4Found]<<" "<<fabs(vdpreprocessedMZ[iIndex4Found]-dTarget)<<endl;
	return bReVal;
}

void MS2Scan::WeightCompare(const string& sPeptide, vector<bool> & vbFragmentZ2) {
	int i, cumscore = 0;
	double threshold;
	vector<int> subscore;

	for (i = 1; i < (int) sPeptide.length(); ++i) {
		if (sPeptide[i] == ']') {
			// hit the C terminus and break out of the loop
			break;
		}

		if (isalpha(sPeptide[i])) {
			if ((sPeptide[i] == 'R') || (sPeptide[i] == 'H') || (sPeptide[i] == 'K'))
				cumscore += 5;
			else if ((sPeptide[i] == 'Q') || (sPeptide[i] == 'N'))
				cumscore += 3;
			else
				cumscore += 1;
			subscore.push_back(cumscore);
		}
	}
	threshold = subscore.back() / 2.0;
	for (i = 0; i < (int) subscore.size(); i++)
		vbFragmentZ2.push_back(subscore.at(i) > threshold);
}

void MS2Scan::saveScore(const double & dScore, const Peptide * currentPeptide, vector<PeptideUnit *> & vpTopPeptides,
		string sScoreFunction) {
	PeptideUnit * copyPeptide;

	if (sScoreFunction == "WeightSum") {
		inumberofWeightSumScore++;
		dsumofWeightScore += dScore;
		dsumofSquareWeightSumScore += dScore * dScore;
	}

	//vdAllScores.push_back(dScore);

	if ((int) vpTopPeptides.size() < TOP_N) {
		copyPeptide = new PeptideUnit;
		copyPeptide->setPeptideUnitInfo(currentPeptide, dScore, sScoreFunction);
		vpTopPeptides.push_back(copyPeptide);
		sort(vpTopPeptides.begin(), vpTopPeptides.end(), GreaterScore);
	} else if (dScore > vpTopPeptides.at(TOP_N - 1)->dScore) {
		copyPeptide = new PeptideUnit;
		copyPeptide->setPeptideUnitInfo(currentPeptide, dScore, sScoreFunction);
		delete vpTopPeptides.at(TOP_N - 1);
		vpTopPeptides.at(TOP_N - 1) = copyPeptide;
		sort(vpTopPeptides.begin(), vpTopPeptides.end(), GreaterScore);
	}
}

void MS2Scan::saveScoreSIP(const double & dScore, const Peptide * currentPeptide, vector<PeptideUnit *> & vpTopPeptides,
		string sScoreFunction) {
	PeptideUnit * copyPeptide;

	if (sScoreFunction == "WeightSum") {
		inumberofWeightSumScore++;
		dsumofWeightScore += dScore;
		dsumofSquareWeightSumScore += dScore * dScore;
	}

	if ((int) vpTopPeptides.size() < TOP_N_SIP) {
		copyPeptide = new PeptideUnit;
		copyPeptide->setPeptideUnitInfo(currentPeptide, dScore, sScoreFunction);
		copyPeptide->setIonMassProb(currentPeptide);
		vpTopPeptides.push_back(copyPeptide);
		sort(vpTopPeptides.begin(), vpTopPeptides.end(), GreaterScore);
	} else if (dScore > vpTopPeptides.at(TOP_N_SIP - 1)->dScore) {
		copyPeptide = new PeptideUnit;
		copyPeptide->setPeptideUnitInfo(currentPeptide, dScore, sScoreFunction);
		copyPeptide->setIonMassProb(currentPeptide);
		delete vpTopPeptides.at(TOP_N_SIP - 1);
		vpTopPeptides.at(TOP_N_SIP - 1) = copyPeptide;
		sort(vpTopPeptides.begin(), vpTopPeptides.end(), GreaterScore);
	}
}

bool MS2Scan::GreaterScore(PeptideUnit * p1, PeptideUnit * p2) {
	return ((p1->dScore) > (p2->dScore));
}

void MS2Scan::scorePeptides() {

	if (bSetMS2Flag) {
		if (isMS2HighRes)
			scorePeptidesHighMS2();
		else

			scorePeptidesLowMS2();

		vpPeptides.clear();
	}
}

void MS2Scan::scorePeptidesMVH(vector<double> * sequenceIonMasses, vector<double> * pdAAforward, vector<double> * pdAAreverse, vector<char> * Seqs) {
	if (!bSkip) {
		int i = 0;
		double dMvh = 0;
		for (i = 0; i < (int) vpPeptides.size(); i++) {
			dMvh = 0;
			if (!mergePeptide(vpWeightSumTopPeptides, vpPeptides.at(i)->getPeptideSeq(), vpPeptides.at(i)->getProteinName())) {
				if (MVH::ScoreSequenceVsSpectrum(vpPeptides.at(i)->sNeutralLossPeptide, this, sequenceIonMasses, pdAAforward, pdAAreverse, dMvh, Seqs)) {
					saveScore(dMvh, vpPeptides.at(i), vpWeightSumTopPeptides);
				}
			}

		}
	} else {
		// cerr << "error scorePeptidesMVH" << endl;
		// exit(1);
	}
	vpPeptides.clear();
}

void MS2Scan::scorePeptidesXcorr(bool *pbDuplFragment, double * _pdAAforward, double * _pdAAreverse, unsigned int *** _uiBinnedIonMasses) {
	int i = 0;
	double dXcorr = 0;
	for (i = 0; i < (int) vpPeptides.size(); i++) {
		dXcorr = 0;
		if (!mergePeptide(vpWeightSumTopPeptides, vpPeptides.at(i)->getPeptideSeq(), vpPeptides.at(i)->getProteinName())) {
			if (CometSearchMod::ScorePeptides(&(vpPeptides.at(i)->sNeutralLossPeptide), pbDuplFragment, _pdAAforward, _pdAAreverse, this, _uiBinnedIonMasses,
					dXcorr)) {
				if (dXcorr > 0) {
					saveScore(dXcorr, vpPeptides.at(i), vpWeightSumTopPeptides);
				}
			}
		}
	}
	vpPeptides.clear();
}

void MS2Scan::preprocessHighMS2() {
	initialPreprocess();
	if (bSetMS2Flag) {
//	normalizeMS2scan(); // for ranksum only
//	setIntensityThreshold(); // for ranksum only
//	filterMS2scan(); // for ranksum only
		vdpreprocessedMZ = vdMZ; // for weightsum only
		vdpreprocessedIntensity = vdIntensity; // for weightsum only
		vipreprocessedCharge = viCharge; // for weightsum only
		binCalculation2D();
		//sortPreprocessedIntensity(); // for ranksum only
		//highMS2scan need the last part for TandemMassSpectrum
	}
}

void MS2Scan::preprocessLowMS2() {
	initialPreprocess();
	if (bSetMS2Flag) {
		normalizeMS2scan();
		setIntensityThreshold();
		filterMS2scan();
		binCalculation2D();
		//sortPreprocessedIntensity(); // for ranksum only
		//highMS2scan need the last part for TandemMassSpectrum
	}
}

void MS2Scan::binCalculation() {
	int i;
	unsigned long int iBinNumber;
	unsigned long int iTarget;
	unsigned long int k;
	unsigned long int iBinRange;

	vbPeakPresenceBins.clear();
// compute the number of vbPeakPresenceBins needed, add 10 to make sure there is a bin for the last peak
	iBinNumber = (unsigned long int) ((vdpreprocessedMZ.back() + 10.0) * bin_res + SMALLINCREMENT);
	// populate all vbPeakPresenceBins with false
//	vbPeakPresenceBins.resize(iBinNumber, false);
	vbPeakPresenceBins.resize(iBinNumber, -1);
	iBinRange = (unsigned long int) (ProNovoConfig::getMassAccuracyFragmentIon() * bin_res + SMALLINCREMENT);	// + 1;
	for (i = 0; i < (int) vdpreprocessedMZ.size(); i++) {
		iTarget = (unsigned long int) (vdpreprocessedMZ[i] * bin_res + SMALLINCREMENT);
//	for(k = iTarget-iBinRange; k <= iTarget+iBinRange; k++)
		for (k = iTarget - iBinRange; k < iTarget + iBinRange; k++)
			if (vbPeakPresenceBins[k] == -1)
				vbPeakPresenceBins[k] = i;
			else if (vdpreprocessedIntensity[vbPeakPresenceBins[k]] < vdpreprocessedIntensity[i])
				vbPeakPresenceBins[k] = i;
	}
}

void MS2Scan::binCalculation2D() {
	unsigned long int iBinNumber, iLowerBound, iUpperBound;
	double dErrRange;
	pair<int, int> initalPair(-1, -1);
	int i;
	unsigned long int j;

	vbPeakPresenceBins2D.clear();
	dErrRange = ProNovoConfig::getMassAccuracyFragmentIon();
	iBinNumber = (unsigned long int) (vdpreprocessedMZ.back() + 10.0 + SMALLINCREMENT);
	vbPeakPresenceBins2D.resize(iBinNumber, initalPair);
	for (i = 0; i < (int) vdpreprocessedMZ.size(); i++) {
		iLowerBound = (unsigned long int) (vdpreprocessedMZ.at(i) - dErrRange + SMALLINCREMENT);
		iLowerBound = ((iLowerBound > 0) ? iLowerBound : 0); // deal with  lower bound less than 0
		iUpperBound = (unsigned long int) (vdpreprocessedMZ.at(i) + dErrRange + SMALLINCREMENT);
		for (j = iLowerBound; j <= iUpperBound; j++) {
			if (vbPeakPresenceBins2D.at(j).first == -1) {
				vbPeakPresenceBins2D.at(j).first = i;
				vbPeakPresenceBins2D.at(j).second = i;
			} else {
				vbPeakPresenceBins2D.at(j).second = i;
			}
		}
	}
}

void MS2Scan::filterMS2scan() {
	int i;
	vector<int> lowpeak;
	lowpeak.clear();
	if (isMS2HighRes) {
		vdpreprocessedMZ = vdMZ;
		vipreprocessedCharge = viCharge;
	}
	for (i = 0; i < (int) vdIntensity.size(); ++i)
		if (isMS2HighRes)
			vdpreprocessedIntensity[i] = vdpreprocessedIntensity[i] / vdMaxMzIntensity[(int) (vdMZ[i] + SMALLINCREMENT)];
		else if (vdpreprocessedIntensity[i] >= vdHighIntensity[(int) (vdMZ[i] + SMALLINCREMENT)]) {
			vdpreprocessedMZ.push_back(vdMZ[i]);
			vdpreprocessedIntensity[i] = vdpreprocessedIntensity[i] / vdMaxMzIntensity[(int) (vdMZ[i] + SMALLINCREMENT)];
			vipreprocessedCharge.push_back(viCharge[i]);
		} else
			lowpeak.push_back(i);
	if (!isMS2HighRes)
		for (i = (int) lowpeak.size() - 1; i >= 0; i--)
			vdpreprocessedIntensity.erase(vdpreprocessedIntensity.begin() + lowpeak.at(i));
	if (vdpreprocessedMZ.front() <= 0) {
		cerr << "ERROR: negative MZ value = " << vdpreprocessedMZ.front() << endl;
		bSetMS2Flag = false;
	}
}

void MS2Scan::setIntensityThreshold() {

	int iMzRange = 50;
	vector<double> vdLocalIntensity;
	int iLowerBound = 0;
	int iUpperBound = 0;
	int iPeakposi = 0;
	int i;

	vdMaxMzIntensity.clear();
	vdHighIntensity.clear();
	vdMaxMzIntensity.resize(5000);
	vdHighIntensity.resize(5000);
	fill(vdMaxMzIntensity.begin(), vdMaxMzIntensity.end(), 1.0);
	for (i = 0; i < (int) vdMaxMzIntensity.size(); ++i) {
		iLowerBound = max(0, i - iMzRange);
		iUpperBound = min(i + iMzRange, (int) vdMzIntensity.size());
		//iPeakposi   = min(iLowerBound+TOPPICKNUM -1, iUpperBound-1);
		vdMaxMzIntensity[i] = *max_element(vdMzIntensity.begin() + iLowerBound, vdMzIntensity.begin() + iUpperBound);
		vdLocalIntensity.clear();
		vdLocalIntensity.resize(2 * iMzRange + 1, 0);
		copy(vdMzIntensity.begin() + iLowerBound, vdMzIntensity.begin() + iUpperBound, vdLocalIntensity.begin());
		iPeakposi = min(TOPPICKNUM - 1, (int) vdLocalIntensity.size() - 1);
		//iPeakposi = max(0, ((int)vdLocalIntensity.size())/2 -1);//wyf : median*2 method
		nth_element(vdLocalIntensity.begin(), vdLocalIntensity.begin() + iPeakposi, vdLocalIntensity.end(), mygreater);
		vdHighIntensity[i] = vdLocalIntensity[iPeakposi];
		//vdHighIntensity[i] = REJECTRATIO*vdLocalIntensity[iPeakposi];//wyf : median*2 method
	}
}

bool MS2Scan::mygreater(double i, double j) {
	return (i > j);
}

void MS2Scan::normalizeMS2scan() {
	int i;
	double dMaxInt;
	// change intensity to relative intensity
	dMaxInt = vdIntensity[getMaxValueIndex(vdIntensity)];
	for (i = 0; i < (int) vdIntensity.size(); ++i)
		vdpreprocessedIntensity.push_back((vdIntensity[i] / dMaxInt) * 100);
//    for (i = 0; i < (int)vdIntensity.size(); ++i)
//	cout<<vdpreprocessedIntensity.at(i)<<endl;
	vdMzIntensity.clear();
	vdMzIntensity.resize(5000);
	fill(vdMzIntensity.begin(), vdMzIntensity.end(), 0.0);
	for (i = 0; i < (int) vdMZ.size(); ++i)
		// if there are multiple peaks in this M/Z bin, use the intensity of the highest peaks
		if (vdpreprocessedIntensity[i] > vdMzIntensity[(int) vdMZ[i]]) //purpose of MzIntensity is for quickly getting the maximum
			vdMzIntensity[(int) vdMZ[i]] = vdpreprocessedIntensity[i];
}

void MS2Scan::initialPreprocess() {
	vdpreprocessedMZ.clear();
	vdpreprocessedIntensity.clear();
	vipreprocessedCharge.clear();

	if (isMS2HighRes)
		bin_res = BIN_RES;
	else
		bin_res = LOW_BIN_RES;
	if ((int) vdMZ.size() == 0 || vdMZ.size() != vdIntensity.size() || vdMZ.size() != viCharge.size()) {
		cerr << "ERROR: Problem with the input mass spectrum" << endl;
		bSetMS2Flag = false;
	} else {
		bSetMS2Flag = true;
		sortPeakList();
		iMaxMZ = (int) min(vdMZ[vdMZ.size() - 1], dParentMZ * iParentChargeState);
		iMinMZ = (int) vdMZ[0];

	}
}

int MS2Scan::getMaxValueIndex(const std::vector<double>& vdData) {
	int iMaxIndex = 0;
	double dMaxValue = 0;
	for (unsigned int i = 0; i < vdData.size(); ++i) {
		if (dMaxValue < vdData[i]) {
			dMaxValue = vdData[i];
			iMaxIndex = (int) i;
		}
	}
	return iMaxIndex;
}

/**
 * sort the MZ in ascending order
 */
void MS2Scan::sortPeakList() {
	int iPeakCount = (int) vdMZ.size();
	int pass, i;
	double dCurrentMZ;
	double dCurrentInt;
	int iCurrentZ;
	for (pass = 1; pass < iPeakCount; pass++) {  // count how many times
												 // This next loop becomes shorter and shorter
		for (i = 0; i < iPeakCount - pass; i++) {
			if (vdMZ[i] > vdMZ[i + 1]) {
				// exchange
				dCurrentMZ = vdMZ[i];
				dCurrentInt = vdIntensity[i];
				iCurrentZ = viCharge[i];

				vdMZ[i] = vdMZ[i + 1];
				vdIntensity[i] = vdIntensity[i + 1];
				viCharge[i] = viCharge[i + 1];

				vdMZ[i + 1] = dCurrentMZ;
				vdIntensity[i + 1] = dCurrentInt;
				viCharge[i + 1] = iCurrentZ;
			}
		}
	}
}

void MS2Scan::sortPreprocessedIntensity() {
	int i, j, ipeakCount, icurrentRank, iUpperBound;
	ipeakCount = (int) vdpreprocessedIntensity.size();
	viIntensityRank.clear();

	iUpperBound = ipeakCount - 1;
	for (i = 0; i < ipeakCount; i++) {
		if (((int) vdpreprocessedMZ[i]) > iMaxMZ) {
			iUpperBound = i - 1;
			break;
		}
		viIntensityRank.push_back(i);
	}
	for (i = iUpperBound; i >= 1; i--)
		for (j = 0; j < i; j++) {												 //move the ith samllest to ith right most position
			if (vdpreprocessedIntensity.at(viIntensityRank.at(j)) < vdpreprocessedIntensity.at(viIntensityRank.at(j + 1))) {
				icurrentRank = viIntensityRank.at(j + 1);
				viIntensityRank.at(j + 1) = viIntensityRank.at(j);
				viIntensityRank.at(j) = icurrentRank;
			}
		}
//    for (i=0; i<ipeakCount; i++)
//      cout<<vdpreprocessedIntensity.at(viIntensityRank.at(i))<<endl;
}

double MS2Scan::CalculateRankSum(double dUvalue, double n1, double n2) {
	double dZscore, dDvalue, dMvalue;

	dMvalue = n1 * n2 / 2.0; // calculating mean value
	dDvalue = sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0); // calculating deviation value
	if (dDvalue > 0)
		dZscore = (dMvalue - dUvalue) / dDvalue;
	else
		dZscore = 0;

	return dZscore;
}

void MS2Scan::preprocess() {
	if (isMS2HighRes)
		preprocessHighMS2();
	else
		preprocessLowMS2();
	cleanup();
}

void MS2Scan::preprocessMvh(multimap<double, double> * pIntenSortedPeakPreData) {
	sortPeakList();
	// MVH
	if (peakData != NULL) {
		delete peakData;
		peakData = NULL;
	}
	if (intenClassCounts != NULL) {
		delete intenClassCounts;
		intenClassCounts = NULL;
	}
	if (pPeakList != NULL) {
		delete pPeakList;
		pPeakList = NULL;
	}
	this->peakData = new map<double, char>();
	this->intenClassCounts = new vector<int>();
	MVH::Preprocess(this, pIntenSortedPeakPreData);
	this->pPeakList = new PeakList(peakData);
	delete peakData;
	peakData = NULL;
}

void MS2Scan::preprocessXcorr() {

}

bool MS2Scan::findProductIon(const vector<double> & vdIonMass, const vector<double> & vdIonProb, const int & iCharge, double & dScoreWeight,
		double & dAverageMZError, double & dMostAbundantObservedMZ, int & iMostAbundantPeakIndex) {
	int iIndex4MaxInt = getMaxValueIndex(vdIonProb);
	double dMaxIntExpectedMZ = (vdIonMass[iIndex4MaxInt] / (double) iCharge) + dProtonMass;
	int iIndex4SelectedFound = 0;
	dScoreWeight = 1.0;
	dAverageMZError = dMassTolerance;

	// search for the most abundant peak
	if (!searchMZ2D(dMaxIntExpectedMZ, iIndex4SelectedFound)) {
		return false;
	}
	// test whether the iCharge is consistant with viZinput
	if (vipreprocessedCharge[iIndex4SelectedFound] != 0) {
		if (vipreprocessedCharge[iIndex4SelectedFound] != iCharge) {
			return false;
		}
	}

	iMostAbundantPeakIndex = iIndex4SelectedFound;
	dMostAbundantObservedMZ = vdpreprocessedMZ[iIndex4SelectedFound];
	double dMostAbundantObservedIntensity = vdpreprocessedIntensity[iIndex4SelectedFound];
	double dMostAbundantMZError = dMostAbundantObservedMZ - dMaxIntExpectedMZ;

	// compute expected MZ and intensity for this product ion
	int iIsotopicEnvelopeSize = vdIonProb.size();
	vector<bool> vbObserved(iIsotopicEnvelopeSize, false);
	vector<double> vdObservedMZ(iIsotopicEnvelopeSize, 0);
	vector<double> vdObservedRelativeInt(iIsotopicEnvelopeSize, 0);
	vector<double> vdMZError(iIsotopicEnvelopeSize, dMassTolerance);
	vector<double> vdExpectedMZ(iIsotopicEnvelopeSize, 0);
	vector<double> vdExpectedRelativeInt(iIsotopicEnvelopeSize, 0);
	// a expected ion have to exceed dMinRelativeExpectedInt to be considered
	int i;
	for (i = 0; i < iIsotopicEnvelopeSize; ++i) {
		vdExpectedMZ[i] = vdIonMass[i] / (double) iCharge + dProtonMass;
		vdExpectedRelativeInt[i] = vdIonProb[i] / vdIonProb[iIndex4MaxInt];
	}

	// set max int value
	vbObserved[iIndex4MaxInt] = true;
	vdObservedMZ[iIndex4MaxInt] = dMostAbundantObservedMZ;
	vdMZError[iIndex4MaxInt] = dMostAbundantMZError;
	vdObservedRelativeInt[iIndex4MaxInt] = 1.0;

	if (iIsotopicEnvelopeSize == 1) {
		// there is only one expected peak in the isotopic distribution
		dScoreWeight = 1.0;
		dAverageMZError = dMostAbundantMZError;
		return true;
	}

	// search for other less abundant peaks ions 
	int iIndex4MostAbundunt;
	double dShiftedExpectedMZ;
	for (i = 0; i < iIsotopicEnvelopeSize; ++i) {
		if (i == iIndex4MaxInt)
			continue;
		// shift expected MZ by the error of the most abundant ion and search for it
		dShiftedExpectedMZ = vdExpectedMZ[i] + dMostAbundantMZError;
		if (searchMZ2D(dShiftedExpectedMZ, dMassTolerance / 2, iIndex4MostAbundunt)) {
			vbObserved[i] = true;
			vdObservedMZ[i] = vdpreprocessedMZ[iIndex4MostAbundunt];
			vdMZError[i] = vdObservedMZ[i] - vdExpectedMZ[i];
			vdObservedRelativeInt[i] = vdpreprocessedIntensity[iIndex4MostAbundunt] / dMostAbundantObservedIntensity;
		}
	}

	// dScoreWeight is a variant of cross correlation 
	// between vdObservedRelativeInt and vdExpectedRelativeInt
	// assume means of both is 0
	double dXY = 0;
	for (i = 0; i < iIsotopicEnvelopeSize; ++i) {
		dXY = dXY + vdExpectedRelativeInt[i] * vdObservedRelativeInt[i];

	}
	double dXX = 0;
	for (i = 0; i < iIsotopicEnvelopeSize; ++i) {
		dXX = dXX + pow(vdExpectedRelativeInt[i], 2);
	}
	double dYY = 0;
	for (i = 0; i < iIsotopicEnvelopeSize; ++i) {
		dYY = dYY + pow(vdObservedRelativeInt[i], 2);
	}
	if (dXX == 0 || dYY == 0)
		dScoreWeight = 0;
	else
		dScoreWeight = dXY / sqrt(dXX * dYY);

	dScoreWeight = dScoreWeight + 1.0;

	// calculate average mass error for this product ion
	double dTotalRelativeIntensity = 0;
	double dTotalMZError = 0;
	for (i = 0; i < (int) vdMZError.size(); ++i) {
		if (vbObserved[i]) {
			dTotalMZError += vdMZError[i] * vdExpectedRelativeInt[i];
			dTotalRelativeIntensity += vdExpectedRelativeInt[i];
		}
	}
	dAverageMZError = dTotalMZError / dTotalRelativeIntensity;
//	dAverageMZError = dMostAbundantMZError;

	return true;
}

////just for SIP
bool MS2Scan::findProductIonSIP(const vector<double> & vdIonMass, const vector<double> & vdIonProb, const int & iCharge, double & dScoreWeight,
		double & dAverageMZError, double & dMostAbundantObservedMZ, int & iMostAbundantPeakIndex) {
	int iIndex4MaxInt = getMaxValueIndex(vdIonProb);
	double dMaxIntExpectedMZ = (vdIonMass[iIndex4MaxInt] / (double) iCharge) + dProtonMass;
	int iIndex4SelectedFound = 0;
	dScoreWeight = 0;
	dAverageMZError = dMassTolerance;

	// search for the most abundant peak
	if (!searchMZ2D(dMaxIntExpectedMZ, iIndex4SelectedFound)) {
		return false;
	}
	iMostAbundantPeakIndex = iIndex4SelectedFound;
	dMostAbundantObservedMZ = vdpreprocessedMZ[iIndex4SelectedFound];
	double dMostAbundantObservedIntensity = vdpreprocessedIntensity[iIndex4SelectedFound];
	double dMostAbundantMZError = dMostAbundantObservedMZ - dMaxIntExpectedMZ;

	// compute expected MZ and intensity for this product ion
	vector<bool> vbObserved(vdIonProb.size(), false);
	vector<double> vdObservedMZ(vdIonProb.size(), 0);
	vector<double> vdObservedRelativeInt(vdIonProb.size(), 0);
	vector<double> vdMZError(vdIonProb.size(), dMassTolerance);
	vector<double> vdExpectedMZ(vdIonProb.size(), 0);
	vector<double> vdExpectedRelativeInt(vdIonProb.size(), 0);
	// a expected ion have to exceed dMinRelativeExpectedInt to be considered
	int i;
	for (i = 0; i < (int) vdIonProb.size(); ++i) {
		vdExpectedMZ[i] = vdIonMass[i] / (double) iCharge + dProtonMass;
		vdExpectedRelativeInt[i] = vdIonProb[i] / vdIonProb[iIndex4MaxInt];
	}

	// set max int value
	vbObserved[iIndex4MaxInt] = true;
	vdObservedMZ[iIndex4MaxInt] = dMostAbundantObservedMZ;
	vdMZError[iIndex4MaxInt] = dMostAbundantMZError;
	vdObservedRelativeInt[iIndex4MaxInt] = 1.0;

	if (vbObserved.size() == 1) {
		// there is only one expected peak in the isotopic distribution
		dScoreWeight = 1.0;
		dAverageMZError = dMostAbundantMZError;
		return true;
	}

	// search for  ions on the left of the most abundant peak
	vector<int> viIndex4Found;
	double dShiftedExpectedMZ;
	unsigned int j;
	double dMinRelativeIntRatio = 0.5;
	double dCurrentRelativeInt = 0;
	for (i = iIndex4MaxInt + 1; i < (int) vdExpectedMZ.size(); ++i) {
		// shift expected MZ by the error of the most abundant ion and search for it
		dShiftedExpectedMZ = vdExpectedMZ[i] + dMostAbundantMZError;
		if (binarySearch(dShiftedExpectedMZ, vdpreprocessedMZ, dMassTolerance / 2, viIndex4Found)) {
			// the observed peaks may be noise peak
			// filter them by requiring their observed relative intensity larger than
			// a minimum ratio of their expected relative intensity
			// and less than a 150% or five times of their expected relative int, whichever smaller
			for (j = 0; j < viIndex4Found.size(); ++j) {
				iIndex4SelectedFound = viIndex4Found[j];
				dCurrentRelativeInt = vdpreprocessedIntensity[iIndex4SelectedFound] / dMostAbundantObservedIntensity;
				if (dCurrentRelativeInt > vdExpectedRelativeInt[i] * dMinRelativeIntRatio && dCurrentRelativeInt < min(vdExpectedRelativeInt[i] * 5, 1.5)) {
					vbObserved[i] = true;
					vdObservedMZ[i] = vdpreprocessedMZ[iIndex4SelectedFound];
					vdMZError[i] = vdObservedMZ[i] - vdExpectedMZ[i];
					vdObservedRelativeInt[i] = dCurrentRelativeInt;
					break;
				}
			}
		} else
			break;
	}
	// search for ions on the left of the most abundant peak
	// identical to the above function except the index
	for (i = iIndex4MaxInt - 1; i >= 0; --i) {
		// shift expected MZ by the error of the most abundant ion and search for it
		dShiftedExpectedMZ = vdExpectedMZ[i] + dMostAbundantMZError;
		if (binarySearch(dShiftedExpectedMZ, vdpreprocessedMZ, dMassTolerance / 2, viIndex4Found))
		//if (searchMZ2D(dShiftedExpectedMZ, viIndex4Found))
				{
			// the observed peaks may be noise peak
			// filter them by requiring their observed relative intensity larger than
			// a minimum ratio of their expected relative intensity
			// and less than a 150% or five times of their expected relative int, whichever smaller
			for (j = 0; j < viIndex4Found.size(); ++j) {
				iIndex4SelectedFound = viIndex4Found[j];
				dCurrentRelativeInt = vdpreprocessedIntensity[iIndex4SelectedFound] / dMostAbundantObservedIntensity;
				if (dCurrentRelativeInt > vdExpectedRelativeInt[i] * dMinRelativeIntRatio && dCurrentRelativeInt < min(vdExpectedRelativeInt[i] * 5, 1.5)) {
					vbObserved[i] = true;
					vdObservedMZ[i] = vdpreprocessedMZ[iIndex4SelectedFound];
					vdMZError[i] = vdObservedMZ[i] - vdExpectedMZ[i];
					vdObservedRelativeInt[i] = dCurrentRelativeInt;
					break;
				}
			}
		} else
			break;
	}

	// calculate score weight for this product ion
	// this formula is still very ad hoc empirical
	dScoreWeight = 0;
	vector<double> vdTempExpectedRelativeInt = vdExpectedRelativeInt;
	vdTempExpectedRelativeInt[iIndex4MaxInt] = 0;
	int iIndex4SecondHighestInt = getMaxValueIndex(vdTempExpectedRelativeInt);
	double dDetectionLimit4RelativeIntensity = 0.5;
	if (vbObserved[iIndex4SecondHighestInt]) {
		// if the second highest peak is found
		// to be implemented
		dScoreWeight = 2.0;
	} else {
		if (vdExpectedRelativeInt[iIndex4SecondHighestInt] > dDetectionLimit4RelativeIntensity) {
			// if the second highest peak is not found and is expected to
			// be found because its relative intensity is higher than the detection limit
			dScoreWeight = 0.5;
		} else {
			// if the second highest peak is not found and is expected to not be found
			dScoreWeight = 1.0;
		}
	}

	// test whether the iCharge is consistant with viZinput
	// if not, lower the dScoreWeight
	if (vipreprocessedCharge[iMostAbundantPeakIndex] != 0) {
		if (vipreprocessedCharge[iMostAbundantPeakIndex] != iCharge) {
			dScoreWeight = dScoreWeight / 2;
		}
	}

	// calculate average mass error for this product ion
	double dTotalRelativeIntensity = 0;
	double dTotalMZError = 0;
	for (i = 0; i < (int) vdMZError.size(); ++i) {
		if (vbObserved[i]) {
			dTotalMZError += vdMZError[i] * vdExpectedRelativeInt[i];
			dTotalRelativeIntensity += vdExpectedRelativeInt[i];
		}
	}
	dAverageMZError = dTotalMZError / dTotalRelativeIntensity;
//      dAverageMZError = dMostAbundantMZError;

	return true;
}

void MS2Scan::scoreWeightSumHighMS2(Peptide* currentPeptide) //it's primaryScorePetide in the previous version
		{

	double dScore = 0;
	int iPeptideLength = currentPeptide->getPeptideLength(), i, j, iMostAbundantPeakIndex = 0;
	int n; // Ion number starting from one
	int z; // charge state
	vector<ProductIon> vFoundIons;
	double dScoreWeight = 0, dMZError = 1, dMostAbundantObservedMZ = 0, dAverageMZError = 0, dBonus4ComplementaryFragmentObserved = 1.0;
	string sPeptide = currentPeptide->getPeptideSeq();
//    cout<<currentPeptide->vvdYionMass.size()<<endl;
	for (n = 0; n < (int) currentPeptide->vvdYionMass.size(); ++n)
		for (z = 1; z <= iParentChargeState; ++z) {
			ProductIon currentIon;
			currentIon.setProductIon('y', n + 1, z);
			// if (ProNovoConfig::getSearchType() == "SIP") {
				if (findProductIonSIP(currentPeptide->vvdYionMass[n], currentPeptide->vvdYionProb[n], z, dScoreWeight, dMZError, dMostAbundantObservedMZ,
						iMostAbundantPeakIndex)) {
					currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
					vFoundIons.push_back(currentIon);
				}
			/* } else {
				if (findProductIon(currentPeptide->vvdYionMass[n], currentPeptide->vvdYionProb[n], z, dScoreWeight, dMZError, dMostAbundantObservedMZ,
						iMostAbundantPeakIndex)) {
					currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
					vFoundIons.push_back(currentIon);
				}
			}*/
		}

	for (n = 0; n < (int) currentPeptide->vvdBionMass.size(); ++n)
		for (z = 1; z <= iParentChargeState; ++z) {
			ProductIon currentIon;
			currentIon.setProductIon('b', n + 1, z);
			// if (ProNovoConfig::getSearchType() == "SIP") {
				if (findProductIonSIP(currentPeptide->vvdBionMass[n], currentPeptide->vvdBionProb[n], z, dScoreWeight, dMZError, dMostAbundantObservedMZ,
						iMostAbundantPeakIndex)) {
					currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
					vFoundIons.push_back(currentIon);
				}
			/*} else {
				if (findProductIon(currentPeptide->vvdBionMass[n], currentPeptide->vvdBionProb[n], z, dScoreWeight, dMZError, dMostAbundantObservedMZ,
						iMostAbundantPeakIndex)) {
					currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
					vFoundIons.push_back(currentIon);
				}
			}*/
		}

	for (i = 0; i < (int) vFoundIons.size(); ++i)
		vFoundIons[i].setComplementaryFragmentObserved(false);

	for (i = 0; i < (int) vFoundIons.size(); ++i)
		for (j = i + 1; j < (int) vFoundIons.size(); ++j)
			if (vFoundIons[i].getIonNumber() + vFoundIons[j].getIonNumber() == iPeptideLength)
				if ((vFoundIons[i].getIonType() == 'y' && vFoundIons[j].getIonType() == 'b')
						|| (vFoundIons[i].getIonType() == 'b' && vFoundIons[j].getIonType() == 'y')) {
					vFoundIons[i].setComplementaryFragmentObserved(true);
					vFoundIons[j].setComplementaryFragmentObserved(true);
				}
	for (i = 0; i < (int) vFoundIons.size(); ++i) {
		dAverageMZError += vFoundIons[i].getMZError();
		//cout<<vFoundIons[i].getMZError()<<endl;
	}
	dAverageMZError = dAverageMZError / (double) vFoundIons.size();

	//   cout<<vFoundIons.size()<<endl;

	for (i = 0; i < (int) vFoundIons.size(); ++i) {
		if (vFoundIons[i].getComplementaryFragmentObserved())
			dBonus4ComplementaryFragmentObserved = 2.0;
		else
			dBonus4ComplementaryFragmentObserved = 1.0;
		// if (ProNovoConfig::getSearchType() == "SIP")
			dScore += ProNovoConfig::scoreError(fabs(vFoundIons[i].getMZError() - dAverageMZError)) * vFoundIons[i].getScoreWeight()
					* dBonus4ComplementaryFragmentObserved;
		/*else
			// no mass error calibration
			dScore += ProNovoConfig::scoreError(fabs(vFoundIons[i].getMZError())) * vFoundIons[i].getScoreWeight() * dBonus4ComplementaryFragmentObserved;
		*/
		//cout<<dScore<<endl;
	}

	// saveScore(dScore, currentPeptide, vpWeightSumTopPeptides, vdWeightSumAllScores);
	saveScoreSIP(dScore, currentPeptide, vpWeightSumTopPeptides);
}

bool MS2Scan::binarySearch(const double& dTarget, const vector<double>& vdList, const double& dTolerance, vector<int>& viIndex4Found) {
	viIndex4Found.clear();
	int low = 0;
	int high = vdList.size() - 1;
	int mid;
	while (low <= high) {
		mid = (low + high) / 2;
		if (vdList[mid] > dTarget + dTolerance) {
			high = mid - 1;
		} else if (vdList[mid] < dTarget - dTolerance) {
			low = mid + 1;
		} else {
			// found vdList[mid] is within the dTolerance range of dTarget
			viIndex4Found.push_back(mid);

			// now find all other values that is within the range
			int upperBound = mid + 1;
			int lowerBound = mid - 1;

			while (upperBound < (int) vdList.size() && vdList[upperBound] < dTarget + dTolerance) {
				viIndex4Found.push_back(upperBound);
				upperBound += 1;
			}

			while (lowerBound >= 0 && vdList[lowerBound] > dTarget - dTolerance) {
				viIndex4Found.push_back(lowerBound);
				lowerBound -= 1;
			}

			sort(viIndex4Found.begin(), viIndex4Found.end());
			return true;
		}
	}
	return false;
}

bool MS2Scan::isAnyScoreInTopN(int _iIndexPeptide, int _iRankThreshold) {
	for (int k = 0; k < ((int)vpWeightSumTopPeptides.at(_iIndexPeptide)->vdScores.size()); k++) {
		if (vpWeightSumTopPeptides.at(_iIndexPeptide)->vdRank.at(k) <= _iRankThreshold) {
			return true;
		}
	}
	return false;
}

bool scoreSort(pair<double, int> _a, pair<double, int> _b) {
	return _a.first > _b.first;
}

void MS2Scan::scoreFeatureCalculationWDPSip(){
	vector<pair<double, int> > vp;
	// double dMax;
	if (vpWeightSumTopPeptides.size() <= 0) {
		return;
	}
	// rank all peptide
	for (int i = 0; i < PeptideUnit::iNumScores; i++) {
		vp.clear();
		for (int j = 0; j < (int) vpWeightSumTopPeptides.size(); j++) {
			vp.push_back(make_pair(vpWeightSumTopPeptides.at(j)->vdScores.at(i), j));
		}
		sort(vp.begin(), vp.end(), scoreSort);
		// dMax = vp.at(0).first;
		if (vp.size() == 1) {
			vpWeightSumTopPeptides.at(0)->vdRank.push_back(1);

		} else {
			for (int j = 0; j < ((int)vpWeightSumTopPeptides.size()); j++) {
				vpWeightSumTopPeptides.at(vp.at(j).second)->vdRank.push_back(j + 1);
			}
		}
	}
}

void MS2Scan::sumIntensity() {
	dSumIntensity = 0;
	dMaxIntensity = 0;
	for (int i = 0; i < (int) vdIntensity.size(); i++) {
		dSumIntensity += vdIntensity.at(i);
		if (dMaxIntensity < vdIntensity.at(i)) {
			dMaxIntensity = vdIntensity.at(i);
		}
	}
}

void MS2Scan::scoreFeatureCalculation() {
	vector<pair<double, int> > vp;
	// double dMax;
	if (vpWeightSumTopPeptides.size() <= 0) {
		return;
	}
	// if any score equals to zero, remove that peptide
	for (int j = 0; j < (int) vpWeightSumTopPeptides.size(); j++) {
		for (int i = 0; i < PeptideUnit::iNumScores; i++) {
			if (vpWeightSumTopPeptides.at(j)->vdScores.at(i) == 0) {
				delete vpWeightSumTopPeptides.at(j);
				vpWeightSumTopPeptides.at(j) = NULL;
				break;
			}
		}
	}
	vector<PeptideUnit *>::iterator i = vpWeightSumTopPeptides.begin();
	while (i != vpWeightSumTopPeptides.end()) {
		if ((*i) == NULL) {
			i = vpWeightSumTopPeptides.erase(i);
		} else {
			++i;
		}
	}
	if (vpWeightSumTopPeptides.size() <= 0) {
		return;
	}
	// rank all peptide
	for (int i = 0; i < PeptideUnit::iNumScores; i++) {
		vp.clear();
		for (int j = 0; j < (int) vpWeightSumTopPeptides.size(); j++) {
			vp.push_back(make_pair(vpWeightSumTopPeptides.at(j)->vdScores.at(i), j));
		}
		sort(vp.begin(), vp.end(), scoreSort);
		// dMax = vp.at(0).first;
		if (vp.size() == 1) {
			vpWeightSumTopPeptides.at(0)->vdRank.push_back(1);

		} else {
			for (int j = 0; j < ((int)vpWeightSumTopPeptides.size()); j++) {
				vpWeightSumTopPeptides.at(vp.at(j).second)->vdRank.push_back(j + 1);
			}
		}
	}
}

double MS2Scan::scoreWeightSum(string * currentPeptide, vector<double> * pvdYionMasses, vector<double> * pvdBionMasses) {
	double dScore = 0;
	// calculate score
	vector<bool> vbFragmentZ2; //true: y side gets the +2; false: b ion gets +2, only calculated for +3 peptides
	int iNumFragments = 0, n, iIndex4MostAbundunt = 0;
	double dBweight = 0, dYweight = 0, dExpectedMZ, dExYweight = 0, dExBweight = 0;
	int iExYreward = 0, iExBreward = 0, iCurrentYreward, iCurrentBreward;

	if (iParentChargeState >= 3) {
		WeightCompare((*currentPeptide), vbFragmentZ2);
	}
	iNumFragments = pvdYionMasses->size();
	for (n = 0; n < iNumFragments; ++n) {
		dYweight = 0;
		dBweight = 0;
		dExpectedMZ = (*pvdYionMasses)[n] + ProNovoConfig::getProtonMass();
		if (searchMZ2D(dExpectedMZ, iIndex4MostAbundunt)) {
			dYweight = ProNovoConfig::scoreError(fabs(vdpreprocessedMZ[iIndex4MostAbundunt] - dExpectedMZ)) + vdpreprocessedIntensity[iIndex4MostAbundunt];
		}
		dExpectedMZ = (*pvdBionMasses)[iNumFragments - n - 1] + ProNovoConfig::getProtonMass();
		if (searchMZ2D(dExpectedMZ, iIndex4MostAbundunt)) {
			dBweight = ProNovoConfig::scoreError(fabs(vdpreprocessedMZ[iIndex4MostAbundunt] - dExpectedMZ)) + vdpreprocessedIntensity[iIndex4MostAbundunt];
		}
		if (iParentChargeState >= 3) {
			if (vbFragmentZ2[n]) {
				dExpectedMZ = ((*pvdYionMasses)[n] + ProNovoConfig::getProtonMass() * 2) / 2;
				if (searchMZ2D(dExpectedMZ, iIndex4MostAbundunt)) {
					dYweight += ProNovoConfig::scoreError(fabs(vdpreprocessedMZ[iIndex4MostAbundunt] - dExpectedMZ))
							+ vdpreprocessedIntensity[iIndex4MostAbundunt];
				}
			} else {
				dExpectedMZ = ((*pvdBionMasses)[iNumFragments - n - 1] + ProNovoConfig::getProtonMass() * 2) / 2;
				if (searchMZ2D(dExpectedMZ, iIndex4MostAbundunt)) {
					dBweight += ProNovoConfig::scoreError(fabs(vdpreprocessedMZ[iIndex4MostAbundunt] - dExpectedMZ))
							+ vdpreprocessedIntensity[iIndex4MostAbundunt];
				}
			}
		}

		//bouns
		if (dYweight > ZERO && dBweight > ZERO) {
			dScore = dScore + (dYweight + dBweight) * 2;
			if (dExYweight > ZERO) {
				dScore += dYweight;
				iCurrentYreward = 3;
				if (iExYreward < 3) {
					dScore += dExYweight;
				}
			} else {
				iCurrentYreward = 2;
			}
			if (dExBweight > ZERO) {
				dScore += dBweight;
				iCurrentBreward = 3;
				if (iExBreward < 3) {
					dScore += dExBweight;
				}
			} else {
				iCurrentBreward = 2;
			}
		} else {
			dScore = dScore + (dYweight + dBweight);
			if (dExYweight > ZERO) {
				dScore += dYweight;
				iCurrentYreward = 2;
				if (iExYreward < 3) {
					dScore += dExYweight;
				}
			} else {
				iCurrentYreward = 0;
			}
			if (dExBweight > ZERO) {
				dScore += dBweight;
				iCurrentBreward = 2;
				if (iExBreward < 3) {
					dScore += dExBweight;
				}
			} else {
				iCurrentBreward = 0;
			}
		}
		dExYweight = dYweight;
		dExBweight = dBweight;
		iExYreward = iCurrentYreward;
		iExBreward = iCurrentBreward;
	}
	return dScore;
}

double MS2Scan::scoreWeightSumHighMS2(string * currentPeptide, vector<vector<double> > * vvdYionMass, vector<vector<double> > * vvdYionProb,
		vector<vector<double> > * vvdBionMass, vector<vector<double> > * vvdBionProb) {

	double dScore = 0;
	int iPeptideLength = 0;
	for (int i = 0; i < (int) currentPeptide->length(); ++i) {
		if (isalpha((*currentPeptide)[i])) {
			iPeptideLength = iPeptideLength + 1;
		}
	}
	int i, j, iMostAbundantPeakIndex = 0;
	int n; // Ion number starting from one
	int z; // charge state
	vector<ProductIon> vFoundIons;
	double dScoreWeight = 0, dMZError = 1, dMostAbundantObservedMZ = 0, dAverageMZError = 0, dBonus4ComplementaryFragmentObserved = 1.0;

	for (n = 0; n < (int) vvdYionMass->size(); ++n) {
		for (z = 1; z <= iParentChargeState; ++z) {
			ProductIon currentIon;
			currentIon.setProductIon('y', n + 1, z);
			if (ProNovoConfig::getSearchType() == "SIP") {
				if (findProductIonSIP((*vvdYionMass)[n], (*vvdYionProb)[n], z, dScoreWeight, dMZError, dMostAbundantObservedMZ, iMostAbundantPeakIndex)) {
					currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
					vFoundIons.push_back(currentIon);
				}
			} else {
				if (findProductIon((*vvdYionMass)[n], (*vvdYionProb)[n], z, dScoreWeight, dMZError, dMostAbundantObservedMZ, iMostAbundantPeakIndex)) {
					currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
					vFoundIons.push_back(currentIon);
				}
			}
		}
	}

	for (n = 0; n < (int) vvdBionMass->size(); ++n) {
		for (z = 1; z <= iParentChargeState; ++z) {
			ProductIon currentIon;
			currentIon.setProductIon('b', n + 1, z);
			if (ProNovoConfig::getSearchType() == "SIP") {
				if (findProductIonSIP((*vvdBionMass)[n], (*vvdBionProb)[n], z, dScoreWeight, dMZError, dMostAbundantObservedMZ, iMostAbundantPeakIndex)) {
					currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
					vFoundIons.push_back(currentIon);
				}
			} else {
				if (findProductIon((*vvdBionMass)[n], (*vvdBionProb)[n], z, dScoreWeight, dMZError, dMostAbundantObservedMZ, iMostAbundantPeakIndex)) {
					currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
					vFoundIons.push_back(currentIon);
				}
			}
		}
	}

	for (i = 0; i < (int) vFoundIons.size(); ++i) {
		vFoundIons[i].setComplementaryFragmentObserved(false);
	}

	for (i = 0; i < (int) vFoundIons.size(); ++i) {
		for (j = i + 1; j < (int) vFoundIons.size(); ++j) {
			if (vFoundIons[i].getIonNumber() + vFoundIons[j].getIonNumber() == iPeptideLength) {
				if ((vFoundIons[i].getIonType() == 'y' && vFoundIons[j].getIonType() == 'b')
						|| (vFoundIons[i].getIonType() == 'b' && vFoundIons[j].getIonType() == 'y')) {
					vFoundIons[i].setComplementaryFragmentObserved(true);
					vFoundIons[j].setComplementaryFragmentObserved(true);
				}
			}
		}
	}
	for (i = 0; i < (int) vFoundIons.size(); ++i) {
		dAverageMZError += vFoundIons[i].getMZError();
	}
	dAverageMZError = dAverageMZError / (double) vFoundIons.size();

	for (i = 0; i < (int) vFoundIons.size(); ++i) {
		if (vFoundIons[i].getComplementaryFragmentObserved()) {
			dBonus4ComplementaryFragmentObserved = 2.0;
		} else {
			dBonus4ComplementaryFragmentObserved = 1.0;
		}
		if (ProNovoConfig::getSearchType() == "SIP") {
			dScore += ProNovoConfig::scoreError(fabs(vFoundIons[i].getMZError() - dAverageMZError)) * vFoundIons[i].getScoreWeight()
					* dBonus4ComplementaryFragmentObserved;
		} else {	// no mass error calibration
			dScore += ProNovoConfig::scoreError(fabs(vFoundIons[i].getMZError())) * vFoundIons[i].getScoreWeight() * dBonus4ComplementaryFragmentObserved;
		}
	}
	return dScore;
}

ProductIon::ProductIon() {
	cIonType = 'x';
	iIonNumber = 0;
	iCharge = 1;
	dMostAbundantMass = 0;
	dMostAbundantMZ = 0;
	dMZError = 0;
	dMassError = 0;
	dScoreWeight = 1.0;
	bComplementaryFragmentObserved = false;
	iMostAbundantPeakIndex = 0;
}

ProductIon::~ProductIon() {

}

void ProductIon::setProductIon(char cIonTypeInput, int iIonNumberInput, int iChargeInput) {
	cIonType = cIonTypeInput;
	iIonNumber = iIonNumberInput;
	iCharge = iChargeInput;
}

void ProductIon::setObservedInfo(double dMZErrorInput, double dWeightInput, double dMostAbundantMZInput, int iMostAbundantPeakIndexInput) {
	dMZError = dMZErrorInput;
	dMassError = dMZError * iCharge;
	dScoreWeight = dWeightInput;
	dMostAbundantMZ = dMostAbundantMZInput;
	double dProtonMass = ProNovoConfig::getProtonMass();
	dMostAbundantMass = dMostAbundantMZ * iCharge - dProtonMass * iCharge;
	iMostAbundantPeakIndex = iMostAbundantPeakIndexInput;
}

void ProductIon::setComplementaryFragmentObserved(bool bComplementaryFragmentObservedInput) {
	bComplementaryFragmentObserved = bComplementaryFragmentObservedInput;
}

int PeptideUnit::iNumScores = 0;

void PeptideUnit::setPeptideUnitInfo(const Peptide* currentPeptide, const double & dScore, string sScoringFunction) {
	dCalculatedParentMass = currentPeptide->getPeptideMass();
	sIdentifiedPeptide = currentPeptide->getPeptideSeq();
	sOriginalPeptide = currentPeptide->getOriginalPeptideSeq();
	sProteinNames = currentPeptide->getProteinName();
	cIdentifyPrefix = currentPeptide->getIdentifyPrefix();
	cIdentifySuffix = currentPeptide->getIdentifySuffix();
	cOriginalPrefix = currentPeptide->getOriginalPrefix();
	cOriginalSuffix = currentPeptide->getOriginalSuffix();

	this->dScore = dScore;
	this->sScoringFunction = sScoringFunction;

	// Sipros Ensemble
	vdScores.push_back(dScore);
	dPepNeutralMass = currentPeptide->getPeptideMass();
	iPepLength = currentPeptide->getPeptideLength();
	sPeptideForScoring = currentPeptide->getPeptideForScoring();
}

void PeptideUnit::setIonMassProb(const Peptide * currentPeptide) {
	this->vvdBionMass = currentPeptide->vvdBionMass;
	this->vvdBionProb = currentPeptide->vvdBionProb;
	this->vvdYionMass = currentPeptide->vvdYionMass;
	this->vvdYionProb = currentPeptide->vvdYionProb;
}

char PeakList::iNULL = 255;

PeakList::PeakList(map<double, char> * _peakData) {
	pMassHub.clear();
	pPeaks.clear();
	pClasses.clear();
	iPeakSize = _peakData->size();
	if (iPeakSize == 0) {
		return;
	}
	pPeaks.resize(iPeakSize, 0);
	pClasses.resize(iPeakSize, 0);
	int j = 0, k = 0;
	for (std::map<double, char>::iterator i = _peakData->begin(); i != _peakData->end(); i++) {
		pPeaks.at(j) = i->first;
		pClasses.at(j) = i->second;
		++j;
	}
	if (j != iPeakSize) {
		cout << "ERRO PeakList 2" << endl;
		exit(1);
	}
	iLowestMass = (int) pPeaks.at(0);
	iHighestMass = (int) pPeaks.at(iPeakSize - 1);
	iMassHubPairSizeMinusOne = (iHighestMass - iLowestMass);
	iMassHubSize = (iHighestMass - iLowestMass + 1) * 2;
	pMassHub.resize(iMassHubSize, ((short)-1));

	int iMass = 0, iMassIdx;
	int iStart, iEnd;
	for (j = 0; j < iPeakSize;) {
		iMass = pPeaks.at(j);
		iStart = j;
		iEnd = j + 1;
		for (k = j + 1; k < iPeakSize; ++k) {
			if (((int) pPeaks.at(k)) > iMass) {
				iEnd = k;
				break;
			}
		}
		if (k == iPeakSize) {
			iEnd = k;
		}
		j = iEnd;
		iMassIdx = 2 * (iMass - iLowestMass);
		if (iMassIdx + 1 >= iMassHubSize) {
			cout << "ERRO PeakList 1" << endl;
			exit(1);
		}
		pMassHub.at(iMassIdx) = iStart;
		pMassHub.at(iMassIdx + 1) = iEnd;
	}
}

PeakList::~PeakList() {
}

char PeakList::end() {
	return iNULL;
}

char PeakList::findNear(double mz, double tolerance) {
	if (iPeakSize == 0) {
		return iNULL;
	}
	double dMin = 1000000, dDiff;
	int iMzU, iMzL, i, j;
	char iClass = iNULL;
	iMzU = (int) (mz + tolerance);
	iMzL = (int) (mz - tolerance);
	if (iMzU < iLowestMass) {
		return iNULL;
	}
	if (iMzL > iHighestMass) {
		return iNULL;
	}
	int iStart = 0, iEnd = iMassHubPairSizeMinusOne;
	if (iMzL >= iLowestMass) {
		iStart = iMzL - iLowestMass;
	}
	if (iMzU <= iHighestMass) {
		iEnd = iMzU - iLowestMass;
	}
	for (; iStart <= iEnd; ++iStart) {
		if (pMassHub.at(iStart * 2) != -1) {
			for (i = pMassHub.at(iStart * 2), j = pMassHub.at(iStart * 2 + 1); i < j; ++i) {
				dDiff = fabs(mz - pPeaks.at(i));
				if (dDiff < dMin) {
					dMin = dDiff;
					iClass = pClasses[i];
				}
			}
		}
	}
	if (dMin < tolerance) {
		return iClass;
	} else {
		return iNULL;
	}
}

int PeakList::size() {
	return iPeakSize;
}

