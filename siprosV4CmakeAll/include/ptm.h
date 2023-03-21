#ifndef _PTM_H
#define _PTM_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "proNovoConfig.h"

using namespace std;

class PTM_List {

  vector <char> _residue;
  vector <double> _mass_shift;
  vector <string> _symbol;

public:

  PTM_List() {}
  ~PTM_List() {}

  char residue(int n) { return _residue[n]; }
  double mass_shift(int n) { return _mass_shift[n]; }
  string symbol(int n) { return _symbol[n]; }

  int size() { return _residue.size(); }

  bool populate_from_xml_config();
  //The current configure file is not xml format, but this function name is kept
  
  bool add_ptm(char, string, double);
  void dump();

};

#endif
