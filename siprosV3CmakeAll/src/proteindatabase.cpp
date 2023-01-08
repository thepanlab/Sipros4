#include "proteindatabase.h"

ProteinDatabase::ProteinDatabase(bool bScreenOutput)
{
	int i;
	vector<string> vsallNames = ProNovoConfig::vsSingleResidueNames;
	slegalChar = "";
	imaxPTM = ProNovoConfig::getMaxPTMcount();
	bstayCurrentProtein = false;
	bstayCurrentOriginalPeptide = false;
	scurrentProtein = "";
	scurrentProteinName = "";
	snextProteinName = "";
	iclCheck = 0;
	orderstring = "[]";
	for (i = 0; i < (int)ProNovoConfig::vsSingleResidueNames.size(); i++)
		if ((ProNovoConfig::vsSingleResidueNames[i].length() == 1) &&
			(isalpha(ProNovoConfig::vsSingleResidueNames[i][0])))
			orderstring += ProNovoConfig::vsSingleResidueNames[i];
	// cout<<orderstring<<endl;
	// orderstring = ORDERSTRING;
	this->bScreenOutput = bScreenOutput;
	for (i = 0; i < (int)vsallNames.size(); i++)
		if (isalpha(vsallNames.at(i).at(0)))
			slegalChar = slegalChar + vsallNames.at(i);
	if ((ProNovoConfig::getSearchType() == "Regular") || (ProNovoConfig::getSearchType() == "SIP"))
	{
		// Currently, in the peptide generating part, SIP and Regular modes are same.
		if (!ptmlist.populate_from_xml_config())
		// the current configure file is not xml format, but we still use this function
		{
			cerr << "Error in parsing PTM rules from config " << endl;
			exit(1);
		}
		Initial_PTM_Map();
	}
	if (ProNovoConfig::getSearchType() == "Mutation")
		sNonAfterCleavage = stringDifference(slegalChar,
											 ProNovoConfig::getCleavageAfterResidues());
}

// ProteinDatabase::ProteinDatabase(const ProteinDatabase& other)
//{
// }

ProteinDatabase::~ProteinDatabase()
{
	db_stream.clear();
	db_stream.close();
}

void ProteinDatabase::loadDatabase()
{
	// open fasta file containing protein sequences
	string sDBFilename;
	sDBFilename = ProNovoConfig::getFASTAfilename();
	db_stream.open(sDBFilename.c_str());
	if (!(db_stream.is_open()))
	{
		cerr << "fail to open Database file " << sDBFilename << endl;
		exit(1);
	}
}

void ProteinDatabase::initializePtmInfo()
// For a original peptide, this function records all positions where ptm may happen.
// This funcion also records related ptm types to each possible position.
{
	int j, residue_id, max_ptm;
	pair<int, vector<pair<string, double>>> ptm_position;
	ptm_position_all.clear();
	for (j = 0; j < (int)(sOriginalPeptide.length()); j++)
	// identify positions on which ptm may happen
	{
		residue_id = orderstring.find(sOriginalPeptide.at(j));
		if (!(ptm_map[residue_id].empty()))
		{
			// cout<<residue_id<<endl;
			ptm_position = make_pair(j, ptm_map[residue_id]);
			ptm_position_all.push_back(ptm_position);
		}
	}

	max_ptm = ProNovoConfig::getMaxPTMcount(); // max ptm allowed for each peptides
	// If max_ptm is greater than the number of positions allowing ptms,
	// the number of positions is actual upper bounder of ptm
	icurrentMaxPtm = ((max_ptm < ((int)(ptm_position_all.size()))) ? max_ptm
																   : ((int)(ptm_position_all.size())));
	iPtmCount = 1;
	if (iPtmCount <= icurrentMaxPtm)
		initializeCombOrderPTM();
}

void ProteinDatabase::initializeCombOrderPTM()
// In comb_order, put numbers for 0 to iPtmCount-1
{
	int j;
	comb_order.clear();
	for (j = 0; j < iPtmCount; j++)
		comb_order.push_back(j);
	initializePermutationPTM();
}

void ProteinDatabase::initializePermutationPTM()
{
	int i;
	ptm_order.clear();
	// ptm_order is for which ptm happen on positions
	ele_num.clear();
	// ele_num is for numbers of ptm may happen on positions
	for (i = 0; i < (int)comb_order.size(); i++)
	{
		ptm_order.push_back(0);
		ele_num.push_back((int)ptm_position_all.at(comb_order.at(i)).second.size());
	}
}

bool ProteinDatabase::getNextPeptidePTM(Peptide *currentPeptide)
// for original peptides and ptm peptides
{
	bool bReVal;
	if (!bstayCurrentOriginalPeptide)
	{
		// the next peptide must be original
		//  only used for generating the first original Peptide
		bReVal = getNextOriginalPeptide(currentPeptide);
		if (bReVal)
			initializePtmInfo();
	}
	else
	{
		bstayCurrentOriginalPeptide = getNextPtmPeptide(currentPeptide);
		if (!bstayCurrentOriginalPeptide)
		{
			bReVal = getNextOriginalPeptide(currentPeptide);
			if (bReVal)
				initializePtmInfo();
		}
		else
			bReVal = true;
	}
	return bReVal;
}

bool ProteinDatabase::getNextPeptideMutation(Peptide *currentPeptide)
// The residue mutation before or after the cleavage sites can silence or create cleavage sites
// The current version just considers impact of mutation before cleavage sites
{
	bool bReVal;

	if (!bstayCurrentOriginalPeptide)
	{
		// the next peptide must be original
		bReVal = getNextOriginalPeptideForMutation(currentPeptide);
		if ((bReVal) && (iclCheck > ProNovoConfig::getMaxMissedCleavages()))
		{
			// mutate all cleavage for the first peptide, right most for others
			bReVal = getNextMutationPeptide(currentPeptide);
			if (!bReVal)
			{
				bstayCurrentOriginalPeptide = false;
				bReVal = getNextPeptideMutation(currentPeptide);
			}
		}
	}
	else
	{
		bstayCurrentOriginalPeptide = getNextMutationPeptide(currentPeptide);
		if (!bstayCurrentOriginalPeptide)
		{
			bReVal = getNextOriginalPeptideForMutation(currentPeptide);
			if ((bReVal) && (iclCheck > ProNovoConfig::getMaxMissedCleavages()))
			{
				// mutate all cleavage for the first peptide, right most for others
				bReVal = getNextMutationPeptide(currentPeptide);
				if (!bReVal)
				{
					bstayCurrentOriginalPeptide = false;
					bReVal = getNextPeptideMutation(currentPeptide);
				}
			}
		}
		else
			bReVal = true;
	}
	//    if (bReVal)
	//	cout<<currentPeptide->getPeptideSeq()<<endl;
	return bReVal;
}

bool ProteinDatabase::getNextPeptide(Peptide *currentPeptide)
{
	bool bReVal = true;
	if (!(bstayCurrentProtein))
		bReVal = getDecoyOrNextProtein();
	if (bReVal)
	{
		if ((ProNovoConfig::getSearchType() == "Regular") || (ProNovoConfig::getSearchType() == "SIP"))
		{
			bstayCurrentProtein = getNextPeptidePTM(currentPeptide);
			if (!bstayCurrentProtein)
				bReVal = getNextPeptide(currentPeptide);
		}
		if (ProNovoConfig::getSearchType() == "Mutation")
		{
			bstayCurrentProtein = getNextPeptideMutation(currentPeptide);
			if (!bstayCurrentProtein)
				bReVal = getNextPeptide(currentPeptide);
		}
	}
	return bReVal;
}

void ProteinDatabase::RemoveIllegalResidue(string &seq)
{
	string legalstr = slegalChar;
	size_t found;
	//	cout<<legalstr<<"bbb"<<endl;
	found = seq.find_first_not_of(legalstr);
	while (found != string::npos)
	{
		if (bScreenOutput)
			cout << "Remove illegal character " << seq.substr(found, 1) << " of " << seq << endl;
		seq.erase(found, 1);
		// cout<<seq<<endl;
		found = seq.find_first_not_of(legalstr);
	}
}

bool ProteinDatabase::getFirstProtein()
{
	bool bnewProtein = true; // false if fail to retrieve a new protein
	bGetDecoy = true;
	string sline = "";
	iProteinId = 1;
	bstayCurrentOriginalPeptide = false;
	bstayCurrentProtein = true;
	scurrentProteinName = snextProteinName;
	snextProteinName = "";
	scurrentProtein = "";
	iclCheck = 0;
	sline.clear();
	getline(db_stream, sline);
	if (bScreenOutput)
		cout << "Processing protein #" << iProteinId << "\r";
	if (sline.at(0) == '>')
	{
		//	cout<<sline<<endl;
		//	cout<<sline.find_first_of(" \t\f\v\n\r")<<endl;
		scurrentProteinName = sline.substr(1, sline.find_first_of(" \t\f\v\n\r") - 1);
		while (!db_stream.eof())
		{
			sline.clear();
			getline(db_stream, sline);
			if (sline == "")
				continue;
			if (sline.find_first_not_of(" \r\n") != string::npos)
			{
				if (sline.at(0) == '>')
				{
					snextProteinName = sline.substr(1, sline.find_first_of(" \t\f\v\n\r") - 1);
					break;
				}
				else
				{
					RemoveIllegalResidue(sline);
					scurrentProtein = scurrentProtein + sline.substr(0, sline.find_last_not_of("\r\n") + 1);
				}
			}
		}
		setProtein();
	}
	else
		bnewProtein = false;
	return bnewProtein;
}

bool ProteinDatabase::getDecoyOrNextProtein()
{
	if (bGetDecoy)
	{
		scurrentProteinName = "Rev_" + scurrentProteinName;
		reverse(scurrentProtein.begin(), scurrentProtein.end());
		iclCheck = 0;
		bstayCurrentOriginalPeptide = false;
		setProtein();
		bGetDecoy = false;
		return true;
	}
	else
	{
		bGetDecoy = true;
		return getNextProtein();
	}
}

bool ProteinDatabase::getNextProtein()
{
	bool bnewProtein = true; // false if fail to retrieve a new protein
	string sline = "";

	if (snextProteinName == "")
		bnewProtein = false;
	else
	{
		iProteinId++;
		bstayCurrentOriginalPeptide = false;
		scurrentProteinName = snextProteinName;
		snextProteinName = "";
		scurrentProtein = "";
		iclCheck = 0;
		if (bScreenOutput)
			cout << "Processing protein #" << iProteinId << "\r";
		while (!db_stream.eof())
		{
			sline.clear();
			getline(db_stream, sline);
			if (sline == "")
				continue;
			if (sline.find_first_not_of(" \r\n") != string::npos)
			{
				if (sline.at(0) == '>')
				{
					snextProteinName = sline.substr(1, sline.find_first_of(" \t\f\v\n\r") - 1);
					break;
				}
				else
				{
					RemoveIllegalResidue(sline);
					scurrentProtein = scurrentProtein + sline.substr(0, sline.find_last_not_of("\r\n") + 1);
				}
			}
		}
		setProtein();
	}
	return bnewProtein;
}

void ProteinDatabase::Initial_PTM_Map()
// organize ptm information
{
	vector<pair<string, double>> ptm_elm;
	pair<string, double> cur_elm;
	int i, residue_id;

	ptm_elm.clear();
	// cout<< orderstring.length() <<endl;
	for (i = 0; i < (int)orderstring.length(); i++) // initialize the vector
		ptm_map.push_back(ptm_elm);
	for (i = 0; i < ptmlist.size(); i++)
	{
		//		cout<<orderstring<<"aa"<<endl;
		//		cout<<ptmlist.residue(i)<<"bb"<<endl;
		residue_id = orderstring.find(ptmlist.residue(i));
		cur_elm = make_pair<string, double>(ptmlist.symbol(i), ptmlist.mass_shift(i));
		ptm_map[residue_id].push_back(cur_elm);
	}
}

bool ProteinDatabase::GenerateNextComb(std::vector<int> &comb_order, const int &total_num)
{
	int i, j, r_num, ori_val;
	bool re_val = false;
	r_num = comb_order.size();
	for (i = (r_num - 1); i >= 0; i--)
		if (comb_order[i] < total_num - r_num + i)
		{
			ori_val = comb_order[i];
			for (j = i; j < r_num; j++)
				comb_order[j] = ori_val + 1 + j - i;
			re_val = true;
			break;
		}
	//	cout<<re_val<<endl;
	return re_val;
}

bool ProteinDatabase::GenerateNextPTM(std::vector<int> &ptm_order, const std::vector<int> &ele_num)
// based on given comb_order
{
	bool re_val = false;
	int i, j;
	// cout<<ptm_order[0]<<endl;

	for (i = 0; i < (int)ptm_order.size(); i++)
		if (ptm_order[i] < (ele_num[i] - 1))
		{
			ptm_order[i] += 1;
			for (j = 0; j < i; j++)
				ptm_order[j] = 0;
			re_val = true;
			break;
		}

	// cout<<ptm_order[0]<<endl<<endl;
	return re_val;
}

bool ProteinDatabase::isCleavageSite(char c1, char c2)
// Return true, if this is a real cleavage site
{
	return ((ProNovoConfig::getCleavageAfterResidues().find(c1, 0) != string::npos) && (ProNovoConfig::getCleavageBeforeResidues().find(c2, 0) != string::npos));
}

void ProteinDatabase::setProtein()
{
	int i;
	double msum = 0.0;
	//"M" sometimes is the test start
	if ((ProNovoConfig::getTestStartRemoval()) && (scurrentProtein.at(0) == 'M'))
		scurrentProtein = scurrentProtein.substr(1);

	vicleavageSite.clear();
	// Technically, the first residue is considered as a residue after cleavage site,
	// while the last residue is considered as a residue before cleavage site.
	vicleavageSite.push_back(-1);
	for (i = 0; i < (int)(scurrentProtein.length() - 1); i++)
		if (isCleavageSite(scurrentProtein.at(i), scurrentProtein.at(i + 1)))
			vicleavageSite.push_back(i);
	vicleavageSite.push_back(scurrentProtein.length() - 1);

	ibeginCleavagePos = -1;

	mass_map.clear();
	msum_map.clear();
	for (i = 0; i < (int)scurrentProtein.length(); i++)
	{
		// msum += aa.MonoisotopicMass(scurrentProtein.at(i));
		msum += ProNovoConfig::getResidueMass(scurrentProtein.substr(i, 1));
		msum_map.push_back(msum);
		mass_map.push_back(ProNovoConfig::getResidueMass(scurrentProtein.substr(i, 1)));
	}
}

bool ProteinDatabase::isPeptideLengthGood(int ipepLength)
// return true, if the length of given peptide is acceptable
{
	return ((ipepLength >= ProNovoConfig::getMinPeptideLength()) && (ipepLength <= ProNovoConfig::getMaxPeptideLength()));
}

bool ProteinDatabase::getNextOriginalPeptide(Peptide *originalPeptide)
{
	int iendCleavagePos;
	bool bReVal = true;
	bstayCurrentOriginalPeptide = true;
	// get peptides without missed cleavage sites first, then with 1 missed cleavage site,
	// 2 missed cleavage sites, ...
	if (iclCheck <= ProNovoConfig::getMaxMissedCleavages())
	{
		ibeginCleavagePos++;
		iendCleavagePos = ibeginCleavagePos + iclCheck + 1;
		if (iendCleavagePos < (int)vicleavageSite.size())
		{
			if (isPeptideLengthGood(vicleavageSite.at(iendCleavagePos) - vicleavageSite.at(ibeginCleavagePos)))
			{
				sOriginalPeptide = "[" + scurrentProtein.substr(vicleavageSite.at(ibeginCleavagePos) + 1, vicleavageSite.at(iendCleavagePos) - vicleavageSite.at(ibeginCleavagePos)) + "]";
				dOriginalPeptideMass = msum_map[vicleavageSite.at(iendCleavagePos)] + ProNovoConfig::getTerminusMassN() + ProNovoConfig::getTerminusMassC();
				if (vicleavageSite.at(ibeginCleavagePos) >= 0)
					dOriginalPeptideMass = dOriginalPeptideMass - msum_map[vicleavageSite.at(ibeginCleavagePos)];
				setPeptideInfo(originalPeptide, sOriginalPeptide, dOriginalPeptideMass);
				// originalPeptide->setPeptide(sOriginalPeptide, sOriginalPeptide, scurrentProteinName, vicleavageSite.at(ibeginCleavagePos)+1,dOriginalPeptideMass );
			}
			else
				bReVal = getNextOriginalPeptide(originalPeptide);
		}
		else
		{
			iclCheck++;
			ibeginCleavagePos = -1;
			bReVal = getNextOriginalPeptide(originalPeptide);
		}
	}
	else
		bReVal = false;
	return bReVal;
}

bool ProteinDatabase::getNextOriginalPeptideForMutation(Peptide *originalPeptide)
// The current version allows at most ONE mutation in a peptide
{
	int iendCleavagePos;
	bool bReVal = true;
	bstayCurrentOriginalPeptide = true;
	imutationPos = 0;
	imutationOrder = -1;
	imutationCleavageCount = 0;
	bLeftSubpeptide = false;
	bRightSubpeptide = false;
	if (iclCheck <= (ProNovoConfig::getMaxMissedCleavages() + 1)) //"+1" for cleavage site mutation
	// mutation can silence at most one missed cleavage site.
	{
		ibeginCleavagePos++;
		iendCleavagePos = ibeginCleavagePos + iclCheck + 1;
		if (iendCleavagePos < (int)vicleavageSite.size())
		{
			if (isPeptideLengthGood(vicleavageSite.at(iendCleavagePos) - vicleavageSite.at(ibeginCleavagePos)))
			{
				sOriginalPeptide = "[" + scurrentProtein.substr(vicleavageSite.at(ibeginCleavagePos) + 1, vicleavageSite.at(iendCleavagePos) - vicleavageSite.at(ibeginCleavagePos)) + "]";
				dOriginalPeptideMass = msum_map[vicleavageSite.at(iendCleavagePos)] + ProNovoConfig::getTerminusMassN() + ProNovoConfig::getTerminusMassC();
				if (vicleavageSite.at(ibeginCleavagePos) >= 0)
					dOriginalPeptideMass = dOriginalPeptideMass - msum_map[vicleavageSite.at(ibeginCleavagePos)];
				setPeptideInfo(originalPeptide, sOriginalPeptide, dOriginalPeptideMass);
				// originalPeptide->setPeptide(sOriginalPeptide, sOriginalPeptide, scurrentProteinName, vicleavageSite.at(ibeginCleavagePos)+1,dOriginalPeptideMass );
			}
			else
				bReVal = getNextOriginalPeptideForMutation(originalPeptide);
		}
		else
		{
			iclCheck++;
			ibeginCleavagePos = -1;
			bReVal = getNextOriginalPeptideForMutation(originalPeptide);
		}
	}
	else
		bReVal = false;
	return bReVal;
}

bool ProteinDatabase::getNextPtmPeptide(Peptide *ptmPeptide)
// only for ptm peptides
{
	int i, iPtmPos;
	bool bReVal = true;
	string scurrentPeptide, sPtm;
	double dcurrentMass, dPtmMass;

	if (iPtmCount > icurrentMaxPtm)
		bReVal = false;
	else
	{
		scurrentPeptide = sOriginalPeptide;
		dcurrentMass = dOriginalPeptideMass;
		for (i = ((int)comb_order.size() - 1); i >= 0; i--)
		{
			iPtmPos = ptm_position_all.at(comb_order.at(i)).first;
			//	    cout<<"i "<<i<<endl;
			//	    cout<<ptm_position_all.at(comb_order.at(i)).second.size()<<endl;
			sPtm = ptm_position_all.at(comb_order.at(i))
					   .second.at(ptm_order.at(i))
					   .first;

			dPtmMass = ptm_position_all.at(comb_order.at(i))
						   .second.at(ptm_order.at(i))
						   .second;
			dcurrentMass += dPtmMass;
			scurrentPeptide.insert(iPtmPos + 1, sPtm);
		}
		setPeptideInfo(ptmPeptide, scurrentPeptide, dcurrentMass);
		// ptmPeptide->setPeptide(scurrentPeptide, sOriginalPeptide, scurrentProteinName, vicleavageSite.at(ibeginCleavagePos)+1,dcurrentMass);
		if (!GenerateNextPTM(ptm_order, ele_num))
		{
			if (!GenerateNextComb(comb_order, (int)(ptm_position_all.size())))
			{
				iPtmCount++;
				if (iPtmCount <= icurrentMaxPtm)
					initializeCombOrderPTM();
			}
			else
				initializePermutationPTM();
		}
	}
	return bReVal;
}

bool ProteinDatabase::verifyPeptide(const string &sMutatedPeptide)
{
	// one mutation can create at most two cleavage sites, so we have to double check it
	bool bReVal;
	int i, iclNum = 0, iclPos;
	for (i = 1; i <= iclCheck; i++)
	{
		iclPos = vicleavageSite.at(ibeginCleavagePos + i) - vicleavageSite.at(ibeginCleavagePos) - 1;
		if (iclPos == imutationPos)
			continue;
		if (isCleavageSite(sMutatedPeptide.at(iclPos), sMutatedPeptide.at(iclPos + 1)))
			iclNum++;
	}

	if (imutationPos < ((int)sMutatedPeptide.length() - 1))
		if (isCleavageSite(sMutatedPeptide.at(imutationPos), sMutatedPeptide.at(imutationPos + 1)))
			iclNum++;

	if (iclNum > ProNovoConfig::getMaxMissedCleavages())
		bReVal = false;
	else
		bReVal = true;
	return bReVal;
}

bool ProteinDatabase::mutatePeptideLessCleavage(Peptide *mutationPeptide,
												const string &sOriginalPeptideContent, int iendCleavagePos, string &sMutatedPeptide)
{
	// This function is used to mutate original peptides with less than
	// maximum missed cleavage sites, so a new missed cleavage site is acceptable.
	// Once a new cleavage site is created, we may have new peptides.
	// the first mutation must be within "Before_Residues"
	// unless it is the first residue of the whole protein
	// the last mutation must be within "After_Residues"
	// unless it is the last residue of the whole protein
	int bReVal = false;
	double dmutationMass;
	if (imutationOrder == -1)
	{
		if ((imutationPos == 0) && (ibeginCleavagePos != 0))
			sMutationOrder = mutationString(ProNovoConfig::getCleavageBeforeResidues(), sOriginalPeptideContent.at(imutationPos));
		else if ((imutationPos == ((int)sOriginalPeptideContent.length() - 1)) && (iendCleavagePos < (int)vicleavageSite.size() - 1))
			sMutationOrder = mutationString(ProNovoConfig::getCleavageAfterResidues(), sOriginalPeptideContent.at(imutationPos));
		else
			sMutationOrder = mutationString(slegalChar, sOriginalPeptideContent.at(imutationPos));
	}
	imutationOrder++;
	if (imutationOrder >= (int)sMutationOrder.length())
	{
		imutationPos++;
		imutationOrder = -1;
		bReVal = getNextMutationPeptide(mutationPeptide);
	}
	else
	{
		bReVal = true;
		mutatePeptide(sOriginalPeptideContent, sMutatedPeptide);
		// setPeptide;
		dmutationMass = dOriginalPeptideMass - ProNovoConfig::getResidueMass(sOriginalPeptideContent.substr(imutationPos, 1)) + ProNovoConfig::getResidueMass(sMutationOrder.substr(imutationOrder, 1));
		if (imutationPos < ((int)sMutatedPeptide.length() - 1))
			subPeptide(sMutatedPeptide); // call this function before [] are added
		sMutatedPeptide = '[' + sMutatedPeptide + ']';
		setPeptideInfo(mutationPeptide, sMutatedPeptide, dmutationMass);
	}
	return bReVal;
}

bool ProteinDatabase::mutatePeptideEqualCleavage(Peptide *mutationPeptide,
												 const string &sOriginalPeptideContent, int iendCleavagePos, string &sMutatedPeptide)
// This function is used to mutate original peptides with
// maximum missed cleavage sites. If a new cleavage site is created by mutation,
// the current peptide can't exist due to the number of cleavage sites

{
	int bReVal = false;
	double dmutationMass;
	if (imutationOrder == -1)
	{
		if ((imutationPos == 0) && (ibeginCleavagePos != 0))
			sMutationOrder = mutationString(ProNovoConfig::getCleavageBeforeResidues(), sOriginalPeptideContent.at(imutationPos));
		else if ((imutationPos == ((int)sOriginalPeptideContent.length() - 1)) && (iendCleavagePos < (int)vicleavageSite.size() - 1))
			sMutationOrder = mutationString(ProNovoConfig::getCleavageAfterResidues(), sOriginalPeptideContent.at(imutationPos));
		else
			sMutationOrder = mutationString(slegalChar, sOriginalPeptideContent.at(imutationPos));
	}
	imutationOrder++;
	if (imutationOrder >= (int)sMutationOrder.length())
	{
		imutationPos++;
		imutationOrder = -1;
		bReVal = getNextMutationPeptide(mutationPeptide);
	}
	else
	{
		mutatePeptide(sOriginalPeptideContent, sMutatedPeptide);
		if ((imutationPos > 0) && (imutationPos < ((int)sMutatedPeptide.length() - 1)))
			subPeptide(sMutatedPeptide);
		if (verifyPeptide(sMutatedPeptide))
		{
			bReVal = true;
			dmutationMass = dOriginalPeptideMass - ProNovoConfig::getResidueMass(sOriginalPeptideContent.substr(imutationPos, 1)) + ProNovoConfig::getResidueMass(sMutationOrder.substr(imutationOrder, 1));
			sMutatedPeptide = '[' + sMutatedPeptide + ']';
			setPeptideInfo(mutationPeptide, sMutatedPeptide, dmutationMass);
		}
		else if (bLeftSubpeptide)
		{
			bReVal = true;
			setLeftSubPeptideInfo(mutationPeptide); // setLeftPeptide;
		}
		else if (bRightSubpeptide)
		{
			bReVal = true;
			setRightSubPeptideInfo(mutationPeptide); // setRightPeptide
		}
		else
			bReVal = getNextMutationPeptide(mutationPeptide);
	}
	return bReVal;
}

bool ProteinDatabase::mutatePeptideMoreCleavage(Peptide *mutationPeptide,
												const string &sOriginalPeptideContent, int iendCleavagePos, string &sMutatedPeptide)
// This function is used to mutate original peptides with one more
// maximum missed cleavage sites. One missed cleavage site needs to
// be silenced for meeting the requirement of maximum missed cleavage sites
{
	// mutate all cleavage sites in turns
	int bReVal = false;
	double dmutationMass;
	if (imutationOrder == -1)
	{
		imutationCleavageCount++;
		if (imutationCleavageCount > iclCheck)
			bReVal = false;
		else
		{
			imutationPos = vicleavageSite.at(ibeginCleavagePos + imutationCleavageCount) - vicleavageSite.at(ibeginCleavagePos) - 1;
			sMutationOrder = mutationString(sNonAfterCleavage, sOriginalPeptideContent.at(imutationPos));
		}
	}
	imutationOrder++;
	if (imutationCleavageCount <= iclCheck)
	{
		if (imutationOrder >= (int)sMutationOrder.length())
		{
			imutationOrder = -1;
			bReVal = getNextMutationPeptide(mutationPeptide);
		}
		else
		{
			mutatePeptide(sOriginalPeptideContent, sMutatedPeptide);
			if (!(verifyPeptide(sMutatedPeptide)))
				bReVal = getNextMutationPeptide(mutationPeptide);
			else
			{
				bReVal = true;
				dmutationMass = dOriginalPeptideMass - ProNovoConfig::getResidueMass(sOriginalPeptideContent.substr(imutationPos, 1)) + ProNovoConfig::getResidueMass(sMutationOrder.substr(imutationOrder, 1));
				sMutatedPeptide = '[' + sMutatedPeptide + ']';
				setPeptideInfo(mutationPeptide, sMutatedPeptide, dmutationMass);
			}
		}
	}
	return bReVal;
}

bool ProteinDatabase::getNextMutationPeptide(Peptide *mutationPeptide)
// Mutation can create a new cleavage site. So we have to consider new peptides
{
	int iendCleavagePos;
	bool bReVal = false;
	string sOriginalPeptideContent; // without[]
	string sMutatedPeptide;

	iendCleavagePos = ibeginCleavagePos + iclCheck + 1;
	sOriginalPeptideContent = sOriginalPeptide.substr(1, sOriginalPeptide.length() - 2);

	if (bLeftSubpeptide)
	{
		bReVal = true;
		setLeftSubPeptideInfo(mutationPeptide); // setLeftPeptide;
	}
	else if (bRightSubpeptide)
	{
		bReVal = true;
		setRightSubPeptideInfo(mutationPeptide); // setRightPeptide
	}
	else if (imutationPos < (int)sOriginalPeptideContent.length())
	{
		if (iclCheck < ProNovoConfig::getMaxMissedCleavages())
			bReVal = mutatePeptideLessCleavage(mutationPeptide, sOriginalPeptideContent,
											   iendCleavagePos, sMutatedPeptide);
		else if (iclCheck == ProNovoConfig::getMaxMissedCleavages())
			bReVal = mutatePeptideEqualCleavage(mutationPeptide, sOriginalPeptideContent,
												iendCleavagePos, sMutatedPeptide);
		else
			bReVal = mutatePeptideMoreCleavage(mutationPeptide, sOriginalPeptideContent,
											   iendCleavagePos, sMutatedPeptide);
	}
	return bReVal;
}

string ProteinDatabase::stringDifference(const string &sStringA, const string &sStringB)
{
	// return characters contained by sStringA, but not by sStringB
	string sDifferenceString = "";
	int i, iPos;
	for (i = 0; i < (int)sStringA.length(); i++)
	{
		iPos = (int)(sStringB.find(sStringA.substr(i, 1)));
		if (iPos == (int)string::npos)
			sDifferenceString = sDifferenceString + sStringA.substr(i, 1);
	}
	return sDifferenceString;
}

string ProteinDatabase::mutationString(const string &sWholeString, const char &chCurrentRedidue)
{
	// return the chars for mutation
	// the usage of this function is to make sure the cleavage site will be kept or slienced,
	// or no cleavage sites will be created
	int iPos;
	string sReString;
	iPos = (int)sWholeString.find(chCurrentRedidue);
	if (iPos == 0)
		sReString = sWholeString.substr(iPos + 1);
	else
		sReString = sWholeString.substr(0, iPos) + sWholeString.substr(iPos + 1);
	return sReString;
}

void ProteinDatabase::mutatePeptide(const string &sOriginalPeptideContent, string &sMutatedPeptide)
{
	if (imutationPos == 0)
		sMutatedPeptide = sMutationOrder.substr(imutationOrder, 1) + sOriginalPeptideContent.substr(imutationPos + 1);
	else
		sMutatedPeptide = sOriginalPeptideContent.substr(0, imutationPos) + sMutationOrder.substr(imutationOrder, 1) + sOriginalPeptideContent.substr(imutationPos + 1);
}

void ProteinDatabase::setPeptideInfo(Peptide *thePeptide, const string &sIdentifyPeptide, const double &dMass)
{
	char cIdentifyPrefix, cIdentifySuffix, cOriginalPrefix, cOriginalSuffix;
	int iBeginPos, iEndPos; // positions of Peptide on the original protein
	iBeginPos = vicleavageSite.at(ibeginCleavagePos) + 1;
	iEndPos = iBeginPos + sOriginalPeptide.length() - 3;
	if (iBeginPos == 0)
		cOriginalPrefix = '-'; //  since no good empty character, '-' won't be printed finally
	else
		cOriginalPrefix = scurrentProtein.at(iBeginPos - 1);
	cIdentifyPrefix = cOriginalPrefix;
	if (iEndPos == ((int)scurrentProtein.length() - 1))
		cOriginalSuffix = '-'; // since no good empty character, '-' won't be printed finally
	else
		cOriginalSuffix = scurrentProtein.at(iEndPos + 1);
	cIdentifySuffix = cOriginalSuffix;
	thePeptide->setPeptide(sIdentifyPeptide, sOriginalPeptide, scurrentProteinName, vicleavageSite.at(ibeginCleavagePos) + 1,
						   dMass, cIdentifyPrefix, cIdentifySuffix, cOriginalPrefix, cOriginalSuffix);
}

void ProteinDatabase::setLeftSubPeptideInfo(Peptide *thePeptide)
{
	char cIdentifyPrefix, cIdentifySuffix, cOriginalPrefix, cOriginalSuffix;
	string sSubOriginalPeptide, sSubIndentifyPeptide;
	int ilength;
	double dSubMass = 0;

	bLeftSubpeptide = false;
	ilength = iLeftSubPeptideEndPos - iLeftSubPeptideBeginPos + 1;
	sSubIndentifyPeptide = '[' + sLeftSubpeptide + ']';
	sSubOriginalPeptide = '[' + scurrentProtein.substr(iLeftSubPeptideBeginPos, ilength) + ']';
	if (iLeftSubPeptideBeginPos == 0)
		cIdentifyPrefix = '-'; //  since no good empty character, '-' won't be printed finally
	else
		cIdentifyPrefix = scurrentProtein.at(iLeftSubPeptideBeginPos - 1);
	cOriginalPrefix = cIdentifyPrefix;
	cOriginalSuffix = scurrentProtein.at(iLeftSubPeptideEndPos + 1);
	cIdentifySuffix = cLeftSubPeptideSuffix;

	if (iLeftSubPeptideEndPos > 0)
		dSubMass += msum_map[iLeftSubPeptideEndPos - 1];
	if (iLeftSubPeptideBeginPos > 0)
		dSubMass = dSubMass - msum_map[iLeftSubPeptideBeginPos - 1];
	dSubMass += ProNovoConfig::getResidueMass(sLeftSubpeptide.substr(ilength - 1, 1));
	dSubMass += ProNovoConfig::getTerminusMassN() + ProNovoConfig::getTerminusMassC();

	thePeptide->setPeptide(sSubIndentifyPeptide, sSubOriginalPeptide, scurrentProteinName, iLeftSubPeptideBeginPos,
						   dSubMass, cIdentifyPrefix, cIdentifySuffix, cOriginalPrefix, cOriginalSuffix);
}

void ProteinDatabase::setRightSubPeptideInfo(Peptide *thePeptide)
{
	char cIdentifyPrefix, cIdentifySuffix, cOriginalPrefix, cOriginalSuffix;
	string sSubOriginalPeptide, sSubIndentifyPeptide;
	int ilength;
	double dSubMass = 0;

	bRightSubpeptide = false;
	ilength = iRightSubPeptideEndPos - iRightSubPeptideBeginPos;
	sSubIndentifyPeptide = '[' + sRightSubpeptide + ']';
	sSubOriginalPeptide = '[' + scurrentProtein.substr(iRightSubPeptideBeginPos, ilength) + ']';
	cIdentifyPrefix = cRightSubPeptidePrefix;
	cOriginalPrefix = scurrentProtein.at(iRightSubPeptideBeginPos - 1);
	if (iRightSubPeptideEndPos == ((int)scurrentProtein.length() - 1))
		cIdentifySuffix = '-'; //  since no good empty character, '-' won't be printed finally
	else
		cIdentifySuffix = scurrentProtein.at(iRightSubPeptideEndPos + 1);
	cOriginalSuffix = cIdentifySuffix;

	dSubMass += msum_map[iRightSubPeptideEndPos];
	dSubMass = dSubMass - msum_map[iRightSubPeptideBeginPos];
	dSubMass += ProNovoConfig::getResidueMass(sRightSubpeptide.substr(0, 1));
	dSubMass += ProNovoConfig::getTerminusMassN() + ProNovoConfig::getTerminusMassC();

	thePeptide->setPeptide(sSubIndentifyPeptide, sSubOriginalPeptide, scurrentProteinName, iRightSubPeptideBeginPos,
						   dSubMass, cIdentifyPrefix, cIdentifySuffix, cOriginalPrefix, cOriginalSuffix);
}

void ProteinDatabase::subPeptide(const string &sMutatedPeptide)
// this function takes care of potential new subpeptides due to the new missed cleavage site
// created by mutation
{
	int i, ileftMostCleavage, irightMostCleavage;
	bool bnewCleavage = true;
	if (imutationPos < ((int)sMutatedPeptide.length() - 1))
	{
		if (isCleavageSite(sMutatedPeptide.at(imutationPos), sMutatedPeptide.at(imutationPos + 1)))
		{
			for (i = 1; i <= iclCheck; i++)
				if (imutationPos == (vicleavageSite.at(ibeginCleavagePos + i) - vicleavageSite.at(ibeginCleavagePos) - 1))
				{
					bnewCleavage = false;
					break;
				}
		}
		else
			bnewCleavage = false;
	}
	else
		bnewCleavage = false;
	// irightMostCleavage is -1 for no missed cleavage site case, otherwise the relative position of right most missing cleavage site
	// ileftMostCleavage is -1 for no missed cleavage site case, otherwise the relative position of left most missing cleavage site
	if (bnewCleavage)
	{
		irightMostCleavage = vicleavageSite.at(ibeginCleavagePos + iclCheck) - vicleavageSite.at(ibeginCleavagePos) - 1;
		if (imutationPos > irightMostCleavage) // for left peptide
			if (isPeptideLengthGood(imutationPos + 1))
			{
				bLeftSubpeptide = true;
				sLeftSubpeptide = sMutatedPeptide.substr(0, imutationPos + 1);
				cLeftSubPeptideSuffix = sMutatedPeptide.at(imutationPos + 1);
				iLeftSubPeptideBeginPos = vicleavageSite.at(ibeginCleavagePos) + 1;
				iLeftSubPeptideEndPos = vicleavageSite.at(ibeginCleavagePos) + imutationPos + 1;
			}
		ileftMostCleavage = vicleavageSite.at(ibeginCleavagePos + 1) - vicleavageSite.at(ibeginCleavagePos) - 1;
		if ((imutationPos < ileftMostCleavage)) // for right peptide
			if (isPeptideLengthGood(sMutatedPeptide.length() - imutationPos - 1))
			{
				bRightSubpeptide = true;
				sRightSubpeptide = sMutatedPeptide.substr(imutationPos + 1);
				cRightSubPeptidePrefix = sMutatedPeptide.at(imutationPos);
				iRightSubPeptideBeginPos = vicleavageSite.at(ibeginCleavagePos) + imutationPos + 2;
				iRightSubPeptideEndPos = vicleavageSite.at(ibeginCleavagePos) + sMutatedPeptide.length();
			}
	}
}
