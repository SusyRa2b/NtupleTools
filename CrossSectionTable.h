// -*- C++ -*-

/*
see .cxx file for some description
*/

#ifndef CROSSSECTIONTABLE_H
#define CROSSSECTIONTABLE_H

enum SUSYProcess {
  ng = 0,
  ns = 1,
  nn = 2,
  ll = 3,
  sb = 4,
  ss = 5,
  tb = 6,
  bb = 7,
  gg = 8,
  sg = 9,
  NotFound = 10
};

#include <map>
#include <iostream>
#include <fstream>

#include <cassert>

#include "TString.h"

class CrossSectionTable {
public:
  CrossSectionTable(const TString & inputFile, bool smsFormat=false);
  ~CrossSectionTable();
  void loadFileToDatabase(const TString & filename); //CMSSM
  void loadFileToDatabaseSMS(const TString & filename); //SMS

  //allow const-safe direct access to the database
  std::map<SUSYProcess,double>  operator[] (const std::pair<int,int> & scanpoint) const;

  //for interactive ROOT
  double getCrossSection(const int m0, const int m12, const SUSYProcess process) const;
  double getSMSCrossSection(const int m0) const;

  std::map<std::pair<int, int>, std::map<SUSYProcess, double> >::iterator begin() { return database_.begin();}
  std::map<std::pair<int, int>, std::map<SUSYProcess, double> >::iterator end()   { return database_.end();}

private:
  std::map<std::pair<int, int>, std::map<SUSYProcess, double> > database_;

};

#endif
