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
  CrossSectionTable(const TString & inputFile, const TString & format="CMSSM" ,const TString & histoname="");
  ~CrossSectionTable();
  void loadFileToDatabase(const TString & filename); //CMSSM
  void loadFileToDatabaseSMS(const TString & filename); //SMS
  void loadFileToDatabaseSMS_simple(const TString & filename); // simple 2 column [mass] [cross-section] format
  void loadFileToDatabaseSMSRoot(const TString & filename,const TString & histoname); //SMS from root file
  void loadFileToDatabasePMSSM(const TString & filename,const bool append); //pMSSM

  void appendFileToDatabasePMSSM(const TString & filename) {loadFileToDatabasePMSSM(filename,true);}

  //allow const-safe direct access to the database
  std::map<SUSYProcess,double>  operator[] (const std::pair<int,int> & scanpoint) const;

  //for interactive ROOT
  double getCrossSection(const int m0, const int m12, const SUSYProcess process) const;
  double getSMSCrossSection(const int m0) const;
  double getPMSSMCrossSection(const int pointindex) const {return getSMSCrossSection(pointindex);}

  std::map<std::pair<int, int>, std::map<SUSYProcess, double> >::iterator begin() { return database_.begin();}
  std::map<std::pair<int, int>, std::map<SUSYProcess, double> >::iterator end()   { return database_.end();}

private:
  std::map<std::pair<int, int>, std::map<SUSYProcess, double> > database_;

};

#endif
