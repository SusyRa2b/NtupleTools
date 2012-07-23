#include "SearchRegion.h"

#include "TString.h"
#include "TRegexp.h"
#include <iostream>

SearchRegion::SearchRegion(TString btagSel,TString htSel,TString metSel,TString oId,bool isSig) : 
  htSelection(htSel),metSelection(metSel),btagSelection(btagSel),owenId(oId),isSIG(isSig) {}
SearchRegion::~SearchRegion() {}
void SearchRegion::Print() const {
  std::cout<<" == "<<btagSelection<<" "<<htSelection<<" "<<metSelection<<std::endl;

}

bool SearchRegion::operator==(const SearchRegion& other) {
  //ht, met, btag cuts must match
  if ( htSelection != other.htSelection ) return false;
  if ( metSelection != other.metSelection ) return false;
  if ( btagSelection != other.btagSelection ) return false;
  //also must both be either SIG or SB
  if (isSIG != other.isSIG) return false;

  //don't care about owenId
  return true;
}

double SearchRegion::getLowEdgeMET() {

  //find MET>=xxx
  TRegexp metge("MET>=[0-9]+");   //can't handle white space or > instead of >=
  TRegexp numbersonly("[0-9]+");

  TString metcut=  metSelection(metge);
  TString metcutval=  metcut(numbersonly);

  return metcutval.Atof();
}

TString SearchRegion::id() const {

  TString theid= btagSelection;
  theid += owenId;
  return theid;
}
