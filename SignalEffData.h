// -*- C++ -*-

#ifndef SIGNALEFFDATA_H
#define SIGNALEFFDATA_H

#include "TString.h"
#include <map>

#include "SystInfo.h"

const TString SignalEffDataSuffix_ = "sedf";
const TString sedfPath = "/home/joshmt/30May_T1tttt_output/"; //for signal systematics

class SignalEffData {
  //keep all of the relevant numbers for signal efficiency for =one point=
  // by one point i mean one sample and one selection
public:
  SignalEffData();
  SignalEffData(TString idtoload);
  ~SignalEffData();

  //  TString id;

  float rawYield; //yield in lumiScale_ invpb, or for SMS really the raw yield
  float effCorr; //factor including all corrections to the efficiency

  float yield_JER;         //NOT persisted when object is saved to file
  float yield_JER_PU;      //NOT persisted when object is saved to file
  float yield_JER_PU_HLT;  //NOT persisted when object is saved to file

//   float eff_derivative_b;
//   float eff_derivative_c;
//   float eff_derivative_l;

//  float sigma_btageff; //TODO persist this
  //  float eff_derivative_b_1s; //NOT persisted (cross-check)

  float totalSystematic(); 
  float totalSystematicWithoutB(); 
  float symmetrize(const TString & which);
  float symmetrizeWithSign(const TString & which);

  float value(const TString &which) {return symmetrize(which);}
  float valuePlus(const TString &which);
  float valueMinus(const TString &which);

  void set(const TString & which, float valminus, float valplus);

  void setFixedForScan();

    //pair is for + / -
  //important convention notes:
  //  values must be fractional ( i.e. 2% is 0.02)
  //  values should preserve sign (not after fabs())
  std::map<TString, SystInfo > systematics;
  //motivation for using a map:
  //  want to be able to clear() it and *regain* the memory footprint, without getting rid of the class object itself
  TString translateVariation(const TString & which) ;

  void write(TString id) const; //write the SignalEffData contents to a file

};

#endif
