// -*- C++ -*-

#ifndef CONFIGURATIONDESCRIPTIONS_H
#define CONFIGURATIONDESCRIPTIONS_H

//class extracted from drawReducedTrees.h
//holds information about the options used to make the reducedTree

#include "TString.h"
#include <map>

class ConfigurationDescriptions {
public:
  ConfigurationDescriptions();
  ~ConfigurationDescriptions();

  //setters
  void setDefault(const TString & description) {default_ = description;}
  void setCorrected(const TString & description) {corrected_=description;}
  void addVariation(const TString & description1, const TString & description2);

  //getters
  TString getDefault() const {return default_;}
  TString getCorrected() const {return corrected_;}
  TString at(const unsigned int i);

  //utilities
  TString getVariedSubstring(const TString &currentVariation);

  //kludge a way to iterate over the elements
  //this fits nicely with the way we were already accessing the old data structure in drawReducedTrees.C
  unsigned int size();

private:
  TString default_;
  TString corrected_;
  std::map< TString, std::pair<TString, TString> > variationPairs;
};

#endif
