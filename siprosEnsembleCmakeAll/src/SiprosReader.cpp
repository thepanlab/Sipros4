/*
 * SiprosReader.cpp
 *
 *  Created on: May 16, 2017
 *      Author: xgo
 */


#include <string>
#include <iostream>
#include "SiprosReader.h"

bool SiprosReader::MzmlReader(string & filename_str, vector<Spectrum> * _vSpectra) {

	// std::string filename = "/home/xgo/Temp/OSU_D2_FASP_Elite_02262014_10.mzML";

	// ICometSearchManager* pCometSearchMgr = GetCometSearchManager();

	// bool bSearchSucceeded = pCometSearchMgr->ReadSpectrumData(filename_str.c_str(), _vSpectra);

	// ReleaseCometSearchManager();

	bool bSearchSucceeded = ReadSpectrumData(filename_str.c_str(), _vSpectra);

	return bSearchSucceeded;
}

static void SetMSLevelFilter(MSReader &mstReader) {
	vector<MSSpectrumType> msLevel;

	msLevel.push_back(MS2);

	mstReader.setFilter(msLevel);
}

bool SiprosReader::ReadSpectrumData(const char* cFile, vector<Spectrum> * _vSpectra) {
	bool bSucceeded = true;

	// For file access using MSToolkit.
	MSReader mstReader;

	// We want to read only MS2 scans.
	SetMSLevelFilter(mstReader);

	bSucceeded = LoadSpectra(mstReader, cFile, _vSpectra);

	return bSucceeded;
}

bool SiprosReader::LoadSpectra(MSReader &mstReader, const char* cFile, vector<Spectrum> *  _vSpectra) {
	Spectrum mstSpectrum;           // For holding spectrum.
	bool _bFirstScan = true;
	int iTotalScans = 0;
	// Load all input spectra.
	while (true) {
		// Loads in MSMS spectrum data.
		if (_bFirstScan) {
			PreloadIonsFile(mstReader, mstSpectrum, cFile, false);
			if(_vSpectra != NULL)_vSpectra->push_back(mstSpectrum);
			_bFirstScan = false;
		} else {
			PreloadIonsFile(mstReader, mstSpectrum, NULL, true);
			if(_vSpectra != NULL)_vSpectra->push_back(mstSpectrum);
		}

		iTotalScans++;

		if (mstSpectrum.getScanNumber() >= mstReader.getLastScan()) {
			// cout << "Done";
			break;
		}
	}

	return true;
}

void SiprosReader::PreloadIonsFile(MSReader &mstReader, Spectrum &spec, const char* cFile, bool bNext, int scNum) {
	if (!bNext) {
		mstReader.readFile(cFile, spec, scNum);
	} else {
		mstReader.readFile(NULL, spec);
	}
}
