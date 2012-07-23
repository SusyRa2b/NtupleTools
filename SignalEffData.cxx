#include "SignalEffData.h"

//#include <ifstream>

#include <iostream>
#include <fstream>

#include <cmath>
using namespace std;

void SignalEffData::write(TString id) const {
  //goal: be able to write to an ascii file all of the important data members,
  //such that I can destroy the object, and recreate it later using the content of the file

  //filename constructed from id
  TString filename = "SignalEffData.";
  filename += id;
  filename +=".";
  filename += SignalEffDataSuffix_;

  filename.Prepend(sedfPath);

  ofstream output(filename.Data());
  output<<rawYield<<endl<<effCorr<<endl;

//   output<<eff_derivative_b<<endl;
//   output<<eff_derivative_c<<endl;
//   output<<eff_derivative_l<<endl;

  for (map<TString, SystInfo >::const_iterator isyst=systematics.begin(); isyst!=systematics.end(); ++isyst) {
    output << isyst->first<<" ";
    isyst->second.write( &output );
  }

  output.close();
  cout<<" == done writing output file: "<<filename<<endl;

}

SignalEffData::SignalEffData(TString idtoload) :
  //  id(idtoload),
  rawYield(0),
  effCorr(1),
  yield_JER(0), yield_JER_PU(0),yield_JER_PU_HLT(0)//,eff_derivative_b(0),eff_derivative_c(0),eff_derivative_l(0)
{
  //filename constructed from id
  TString filename = "SignalEffData.";
  filename += idtoload;

  filename.Prepend(sedfPath);

  ifstream input(filename.Data()); //should check that this is good
  input>>rawYield;
  input>>effCorr;

//   input>>eff_derivative_b;
//   input>>eff_derivative_c;
//   input>>eff_derivative_l;

  TString akey;
  while (input>>akey) {
    //load the SystInfo from the file using the special ctor
    cout<<"==loading key "<<akey<<endl; //DEBUG
    systematics[akey] = SystInfo(&input);
  }

  input.close();
}

SignalEffData::SignalEffData() : 
  //  id("noname"),
  rawYield(0),
  effCorr(1),
  yield_JER(0), yield_JER_PU(0),yield_JER_PU_HLT(0)//,eff_derivative_b(0),eff_derivative_c(0),eff_derivative_l(0)
{ 

  systematics["JES"] = SystInfo();             
  systematics["btag"] = SystInfo();            
  systematics["lftag"] = SystInfo();            
  systematics["PDF"] = SystInfo();             
  systematics["MET"] = SystInfo();             
  systematics["PU"] = SystInfo();              
  systematics["JER"] = SystInfo();             
  systematics["kFactor"] = SystInfo();         
  systematics["cleaning"] = SystInfo(1e-2,1e-2,1);
  systematics["LepVeto"] = SystInfo(3e-2,3e-2,1); 
  systematics["trigger"] = SystInfo(); //now going to be done directly in likelihood (at least for MHT leg)
  systematics["lumi"] = SystInfo(2.2e-2 , 2.2e-2,1);   //final 2011 number is 2.2%

  //list of signal systematics:
  //  JES
  // btag efficiency
  // PDFs (acceptance)
  // (RA2 says they don't do PDF uncertainties on the cross section. so i won't either)
  //  unc energy (MET)

  // PU
  //  JER
  //  trigger eff
  //  MET cleaning
  //  lepton veto eff
  // lumi
  // for mSugra, NLO cross section (k factor)

}

SignalEffData::~SignalEffData() 
{
  systematics.clear();
}

void SignalEffData::setFixedForScan() {
  systematics["PU"].plus = 1e-2;
  systematics["PU"].minus = -systematics["PU"].plus;
  systematics["PU"].status = 1;

  systematics["JER"].plus = 1e-2;
  systematics["JER"].minus = -systematics["JER"].plus;
  systematics["JER"].status = 1;

}


TString SignalEffData::translateVariation(const TString & which) {
  //this is a nasty hack that i do not like
  //the variations are known only by their differences in name
  //these are easily human-readable but don't fit in well with the names I chose to use
  //in this class SignalEffData

  //for now I will just "translate" them here in a hard-coded way.
  //maybe i should use .Contains() but that can be dangerous

  // the b-tag ones are not needed anymore
//   if (which == "BTagEff03") return "btag"; //hardcoding this 03 is a really bad idea...
//   else if (which == "BTagEff04") return "btag"; //a bad idea indeed...
  if (which=="JERbias") return "JER";
  else if (which=="JES0") return "JES";
  else if (which=="METunc0") return "MET";
  else if (which =="PUunc0") return "PU";

  return which;

}

void SignalEffData::set(const TString & which, float valminus, float valplus) {

  TString translatedWhich=  translateVariation(which);
  
  map<TString, SystInfo >::iterator it=systematics.find(translatedWhich);
  if (it==systematics.end() ) {
    cout<<"ERROR -- cannot find in systematics list: "<<translatedWhich<<endl;
    return;
  }

  it->second.minus = valminus;
  it->second.plus = valplus;
  it->second.status = 2;
  
}

float SignalEffData::symmetrize(const TString & which) {

  float s=-1;

 map<TString, SystInfo >::iterator it=systematics.find(which);
  if (it==systematics.end() ) {
    cout<<"ERROR -- cannot find in systematics list: "<<which<<endl;
  }
  //realizing now that i symmetrized different results differently!
  //for now maintain strict consistency
  else  if ( it->first=="kFactor" ) { //average the 2 parts of the pair
    s = 0.5* (fabs(it->second.plus) + fabs(it->second.minus));
  }
  else { //return the larger deviation
    //for PDF uncertainties for now just store the pre-processed results
    float var1= fabs( it->second.plus);
    float var2= fabs( it->second.minus);
    s = var2>var1? var2:var1;
  }
  
  return s;
}

float SignalEffData::symmetrizeWithSign(const TString & which) {

  float s=-1;

  map<TString, SystInfo >::iterator it=systematics.find(which);
  if (it==systematics.end() ) {
    cout<<"ERROR -- cannot find in systematics list: "<<which<<endl;
  }
  else {
    s=    (it->second.plus - it->second.minus)*0.5;
  }
  
  return s;
}

float SignalEffData::valuePlus(const TString & which) {

  float s=-1;

  map<TString, SystInfo >::iterator it=systematics.find(which);
  if (it==systematics.end() ) {
    cout<<"ERROR -- cannot find in systematics list: "<<which<<endl;
  }
  else { //return the larger deviation
    s=it->second.plus;
  }  
  return s;
}
float SignalEffData::valueMinus(const TString & which) {

  float s=-1;

  map<TString, SystInfo >::iterator it=systematics.find(which);
  if (it==systematics.end() ) {
    cout<<"ERROR -- cannot find in systematics list: "<<which<<endl;
  }
  else { //return the larger deviation
    s=it->second.minus;
  }  
  return s;
}

float SignalEffData::totalSystematic() {

  float total2=0;
  for (map<TString, SystInfo >::iterator isyst=systematics.begin(); isyst!=systematics.end(); ++isyst) {
    if ( isyst->second.status == 0) { 
      cout<<"WARNING -- systematic is unset! "<<isyst->first<<endl;
    }
    total2 += pow( symmetrize(isyst->first),2);
  }

  return 100*sqrt(total2);
}
float SignalEffData::totalSystematicWithoutB() {

  float total2=0;
  for (map<TString, SystInfo >::iterator isyst=systematics.begin(); isyst!=systematics.end(); ++isyst) {
    if ( isyst->second.status == 0) { 
      cout<<"WARNING -- systematic is unset! "<<isyst->first<<endl;
    }
    if ( isyst->first != "btag" && isyst->first != "lftag")  total2 += pow( symmetrize(isyst->first),2);
  }

  return 100*sqrt(total2);
}
