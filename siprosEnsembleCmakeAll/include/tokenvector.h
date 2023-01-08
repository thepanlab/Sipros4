/********************************************************/
// tokenvector.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Token vector class duplicating various PERL functions,
// including split, join, shift, and pop.  Implemented
// as a vector (i.e. split returns a vector of elements).
/********************************************************/

#ifndef TOKEN_H_
#define TOKEN_H_

#include <string>
#include <vector>
#include <iostream>

using namespace std;

// Token vector is an instance of a regular string vector,
// just with some added useful methods that duplicate
// common PERL parsing functions.

class TokenVector : public vector <string>
{
	
public:
	TokenVector() {}
	~TokenVector() {}
	TokenVector(const string&, const string& delimiters = " ");
	void Split(const string&, const string& delimiters = " ");
	string Join(const string& splice = " ");
	string Shift();
	string Pop();
};

#endif /*TOKENVECTOR_H_*/
