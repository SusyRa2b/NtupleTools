// -*- C++ -*-

#ifndef SYSTINFO_H
#define SYSTINFO_H
//very simple data container used by the SignalEffData class

#include <fstream>
#include <iostream>


//class extracted from drawReducedTrees.h

class SystInfo {
public:
  SystInfo(float p=0, float m=0, int s=0);
  ~SystInfo();
  SystInfo( std::ifstream* input);

  float plus;
  float minus;
  //convention:
  //use 0 for numbers that are completely unset
  //use 1 for numbers that are fixed by SignalEffData
  //use 2 for numbers that are set by some actual event counts
  int status;

  void write(std::ofstream* outfile) const;
};


#endif
