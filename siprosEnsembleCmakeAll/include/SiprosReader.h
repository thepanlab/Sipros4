/*
 * SiprosReader.h
 *
 *  Created on: May 17, 2017
 *      Author: xgo
 */

#ifndef SIPROSREADER_H_
#define SIPROSREADER_H_

#include "MSReader.h"
#include "Spectrum.h"

using namespace std;
using namespace MSToolkit;

class SiprosReader {
public:
	bool static MzmlReader(std::string & filename_str, std::vector<Spectrum> * _vSpectra);

	bool static ReadSpectrumData(const char* cFile, vector<Spectrum> * _vSpectra);

	bool static LoadSpectra(MSReader &mstReader, const char* cFile, vector<Spectrum> * _vSpectra);

	void static PreloadIonsFile(MSReader &mstReader, Spectrum &spec, const char* cFile = NULL, bool bNext = false,
			int scNum = 0);
};

#endif /* SIPROSREADER_H_ */
