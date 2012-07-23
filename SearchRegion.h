// -*- C++ -*-

#ifndef SEARCHREGION_H
#define SEARCHREGION_H

//class extracted from drawReducedTrees.h

#include "TString.h"

//special container for holding the "search regions of interest"
class SearchRegion {
public:
  SearchRegion(TString btagSel,TString htSel,TString metSel,TString oId,bool isSig=true);
  ~SearchRegion();
  void Print() const;

  double getLowEdgeMET();

  TString htSelection;
  TString metSelection;
  TString btagSelection;

  TString owenId;
  bool isSIG;

  TString id() const;

  bool operator==(const SearchRegion& other);
  bool operator!=(const SearchRegion& other) {return !((*this)==other);}

};
#endif
