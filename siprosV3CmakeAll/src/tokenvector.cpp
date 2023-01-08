/********************************************************/
// tokenvector.cpp
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Token vector class duplicating various PERL functions,
// including Split, Join, Shift, and Pop.  Implemented
// as a vector (i.e. Split returns a vector of elements).
/********************************************************/

#include "tokenvector.h"

// Constructor.  Given a string and a set of delimiters,
// creates a vector of the substrings.

TokenVector::TokenVector(const string& str, const string& delimiters)
{
  clear();
  Split(str, delimiters);
}

// Removes the first element of the vector and returns it.

string TokenVector::Shift() {
  string str = *begin();
  erase(begin());
  return str;
}

// Removes the last element of the vector and returns it.

string TokenVector::Pop() {
  string str = *(end()-1);
  erase(end()-1);
  return str;
}

// Splits a string into a vector of substrings given
// a set of delimiters.

void TokenVector::Split(const string& str, const string& delimiters)
{
  clear();
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  while (string::npos != pos || string::npos != lastPos) {
    push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiters, pos);
    pos = str.find_first_of(delimiters, lastPos);
  }
}

// Joins a set of substrings together, putting a delimiter
// in between each substring.  Returns the final string.

string TokenVector::Join(const string& splice)
{
  string str = "";
  vector<string>::iterator iter;
  for(iter = begin(); iter < end(); iter++) {
    if(iter != begin()) str.append( splice );
    str.append ( *iter );
  }	
  return str;
}
