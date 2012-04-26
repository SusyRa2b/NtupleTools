// -*- C++ -*-

#include "TStopwatch.h"

#include "TRegexp.h"

TString sedfPath = "/home/joshmt/09Feb_mSUGRA_output/"; //for signal systematics

//useful for playing around with plots in interactive ROOT
TH1D* hinteractive=0;
TH2D* h2d=0;

TLatex* text1=0;
TLatex* text2=0;

//holds a list of the *active* samples (to be plotted)
std::vector<TString> samples_;
//hold a list of all samples
std::set<TString> samplesAll_;
//config descriptions defined below
//these maps use the sample names as keys
//std::map<TString, TFile*> files_;
std::map<TString, std::map<TString, TFile*> > files_;
std::map<TString, TH1D*> histos_;
std::map<TString, UInt_t> sampleColor_;
std::map<TString, TString> sampleOwenName_;
std::map<TString, TString> sampleLabel_;
std::map<TString, UInt_t> sampleMarkerStyle_;
std::map<TString, float> sampleScaleFactor_; //1 by default //implemented only for drawPlots!
TChain* dtree=0;
TH1D* hdata=0;
TH2D* hdata2d=0;

TString currentConfig_;
TH2D* scanSMSngen=0;
TH1D* referenceCrossSectionGluino=0;

//default selection
TString selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1";

//intentionally bogus values
float eff_SB_MHT_             = 0.1;
float eff_SB_MHT_err_[2]      = {0.1, 0.1};
float eff_SB_ldp_MHT_         = 0.1;   
float eff_SB_ldp_MHT_err_[2]  = {0.1, 0.1};
float eff_SIG_MHT_            = 0.1; 
float eff_SIG_MHT_err_[2]     = {0.1, 0.1};
float eff_SIG_ldp_MHT_        = 0.1;
float eff_SIG_ldp_MHT_err_[2] = {0.1, 0.1};
float eff_SIG_SL_MHT_         = 0.1; 
float eff_SIG_SL_MHT_err_[2]  = {0.1, 0.1};
float eff_SB_1m_MHT_          = 0.1;
float eff_SB_1m_MHT_err_[2]   = {0.1, 0.1};
float eff_SB_1e_MHT_          = 0.1; 
float eff_SB_1e_MHT_err_[2]   = {0.1, 0.1};

void setTrigEff(const TString which) {


  if (which=="WideSB") {
    cout<<"Using efficiencies for 150-250 GeV SB"<<endl;
    //150-250 GeV SB (with FINAL lumi weighting -April 5) 
    eff_SB_MHT_             = 0.850;
    eff_SB_MHT_err_[0]      = 0.034;
    eff_SB_MHT_err_[1]      = 0.046;
    eff_SB_ldp_MHT_         = 0.912;   
    eff_SB_ldp_MHT_err_[0]  = 0.024;
    eff_SB_ldp_MHT_err_[1]  = 0.057;
    eff_SIG_MHT_            = 0.981; 
    eff_SIG_MHT_err_[0]     = 0.012;
    eff_SIG_MHT_err_[1]     = 0.036;
    eff_SIG_ldp_MHT_        = eff_SIG_MHT_; 
    eff_SIG_ldp_MHT_err_[0] = eff_SIG_MHT_err_[0]; //due to low stats in SIG-LDP, use the SIG numbers
    eff_SIG_ldp_MHT_err_[1] = eff_SIG_MHT_err_[1];
    eff_SIG_SL_MHT_         = 0.999; 
    eff_SIG_SL_MHT_err_[0]  = 0.001;
    eff_SIG_SL_MHT_err_[1]  = 0.003;
    eff_SB_1m_MHT_          = 0.990;
    eff_SB_1m_MHT_err_[0]   = 0.002;
    eff_SB_1m_MHT_err_[1]   = 0.002;
    eff_SB_1e_MHT_          = 0.955; 
    eff_SB_1e_MHT_err_[0]   = 0.004;
    eff_SB_1e_MHT_err_[1]   = 0.004;
  }
  else if (which=="LowSB") {
    //150-200 GeV SB
    // CAREFUL -- if we really move to this region this we have to split e and mu in the SL SB
    eff_SB_MHT_             = 0.832;
    eff_SB_MHT_err_[0]      = 0.042;
    eff_SB_MHT_err_[1]      = 0.054;
    eff_SB_ldp_MHT_         = 0.91;   
    eff_SB_ldp_MHT_err_[0]  = 0.026;
    eff_SB_ldp_MHT_err_[1]  = 0.059;
    eff_SIG_MHT_            = 0.982; 
    eff_SIG_MHT_err_[0]     = 0.012;
    eff_SIG_MHT_err_[1]     = 0.036;
    eff_SIG_ldp_MHT_        = eff_SIG_MHT_; 
    eff_SIG_ldp_MHT_err_[0] = eff_SIG_MHT_err_[0];
    eff_SIG_ldp_MHT_err_[1] = eff_SIG_MHT_err_[1];
    eff_SIG_SL_MHT_         = 0.999; 
    eff_SIG_SL_MHT_err_[0]  = 0.001;
    eff_SIG_SL_MHT_err_[1]  = 0.001;

    //bullshit numbers! (copied from wide SB)
    eff_SB_1m_MHT_          = 0.990;
    eff_SB_1m_MHT_err_[0]   = 0.002;
    eff_SB_1m_MHT_err_[1]   = 0.002;
    eff_SB_1e_MHT_          = 0.955; 
    eff_SB_1e_MHT_err_[0]   = 0.004;
    eff_SB_1e_MHT_err_[1]   = 0.004;

    //old numbers?
//     eff_SB_SL_MHT_          = 0.996;  //not really true but need new numbers
//     eff_SB_SL_MHT_err_[0]   = 0.002;
//     eff_SB_SL_MHT_err_[1]   = 0.003;
  }
  else if (which=="LowSBSpecialTest1") {
    //150-200 GeV SB
    eff_SB_MHT_             = 0.832;
    eff_SB_MHT_err_[0]      = 0.042;
    eff_SB_MHT_err_[1]      = 0.054;
    eff_SB_ldp_MHT_         = 0.91;   
    eff_SB_ldp_MHT_err_[0]  = 0.026;
    eff_SB_ldp_MHT_err_[1]  = 0.059;
    eff_SIG_MHT_            = 0.841;  //use 200-250 SB number
    eff_SIG_MHT_err_[0]     = 0.059;
    eff_SIG_MHT_err_[1]     = 0.090;
    eff_SIG_ldp_MHT_        = 0.936; 
    eff_SIG_ldp_MHT_err_[0] = 0.034;
    eff_SIG_ldp_MHT_err_[1] = 0.118;
    eff_SIG_SL_MHT_         = 0.996; 
    eff_SIG_SL_MHT_err_[0]  = 0.002;
    eff_SIG_SL_MHT_err_[1]  = 0.003;


    //bullshit numbers! (copied from wide SB)
    eff_SB_1m_MHT_          = 0.990;
    eff_SB_1m_MHT_err_[0]   = 0.002;
    eff_SB_1m_MHT_err_[1]   = 0.002;
    eff_SB_1e_MHT_          = 0.955; 
    eff_SB_1e_MHT_err_[0]   = 0.004;
    eff_SB_1e_MHT_err_[1]   = 0.004;
  }
  else if (which=="NominalSB") {
    cout<<"Using efficiencies for 200-250 GeV SB"<<endl;
    //nominal numbers used for frozen AN etc. SB corresponds to 200-250 GeV
    eff_SB_MHT_             = 0.838;
    eff_SB_MHT_err_[0]      = 0.060;
    eff_SB_MHT_err_[1]      = 0.090;
    eff_SB_ldp_MHT_         = 0.937;   
    eff_SB_ldp_MHT_err_[0]  = 0.034;
    eff_SB_ldp_MHT_err_[1]  = 0.120;
    eff_SIG_MHT_            = 0.981; 
    eff_SIG_MHT_err_[0]     = 0.012;
    eff_SIG_MHT_err_[1]     = 0.036;
    eff_SIG_ldp_MHT_        = eff_SIG_MHT_; 
    eff_SIG_ldp_MHT_err_[0] = eff_SIG_MHT_err_[0];
    eff_SIG_ldp_MHT_err_[1] = eff_SIG_MHT_err_[1]; //due to low stats in SIG-LDP, use the SIG numbers for now.
    eff_SIG_SL_MHT_         = 0.999; 
    eff_SIG_SL_MHT_err_[0]  = 0.001;
    eff_SIG_SL_MHT_err_[1]  = 0.003;
    eff_SB_1m_MHT_          = 1.000; 
    eff_SB_1m_MHT_err_[0]   = 0.000;    
    eff_SB_1m_MHT_err_[1]   = 0.004;
    eff_SB_1e_MHT_          = 0.994; 
    eff_SB_1e_MHT_err_[0]   = 0.003;
    eff_SB_1e_MHT_err_[1]   = 0.006;
  }
  else if (which == "NominalSBTypeI") {
    //200-250 GeV type I MET
    cout<<"Using efficiencies for 200-250 GeV SB (Type I MET)"<<endl;
    eff_SB_MHT_             = 0.857;
    eff_SB_MHT_err_[0]      = 0.043; //plus
    eff_SB_MHT_err_[1]      = 0.085; //minus
    eff_SB_ldp_MHT_         = 0.858;   
    eff_SB_ldp_MHT_err_[0]  = 0.053;
    eff_SB_ldp_MHT_err_[1]  = 0.082;
    eff_SIG_MHT_            = 0.982; 
    eff_SIG_MHT_err_[0]     = 0.012;
    eff_SIG_MHT_err_[1]     = 0.036;
    eff_SIG_ldp_MHT_        = eff_SIG_MHT_; 
    eff_SIG_ldp_MHT_err_[0] = eff_SIG_MHT_err_[0];
    eff_SIG_ldp_MHT_err_[1] = eff_SIG_MHT_err_[1]; //due to low stats in SIG-LDP, use the SIG numbers for now.
    //numbers below here are not done yet
    cout<<" WARNING ---- SL trig eff not done yet!"<<endl;
    eff_SIG_SL_MHT_         = 0.999; 
    eff_SIG_SL_MHT_err_[0]  = 0.001;
    eff_SIG_SL_MHT_err_[1]  = 0.001;
    eff_SB_1m_MHT_          = 0.996; 
    eff_SB_1m_MHT_err_[0]   = 0.002;    
    eff_SB_1m_MHT_err_[1]   = 0.003;
    eff_SB_1e_MHT_          = 0.996; 
    eff_SB_1e_MHT_err_[0]   = 0.002;
    eff_SB_1e_MHT_err_[1]   = 0.003;
  }
  else assert(0);

}

void printEff() {

  cout<<"eff_SB_MHT                "<<eff_SB_MHT_<<endl;
  cout<<"eff_SB_MHT_err_plus       "<<eff_SB_MHT_err_[0]<<endl;
  cout<<"eff_SB_MHT_err_minus      "<<eff_SB_MHT_err_[1]<<endl;

  cout<<"eff_SB_ldp_MHT            "<<eff_SB_ldp_MHT_<<endl;
  cout<<"eff_SB_ldp_MHT_err_plus   "<<eff_SB_ldp_MHT_err_[0]<<endl;
  cout<<"eff_SB_ldp_MHT_err_minus  "<<eff_SB_ldp_MHT_err_[1]<<endl;

  //assert(eff_SB_1e_MHT_ == eff_SB_1m_MHT_);
//   cout<<"eff_SB_sl_MHT             "<<eff_SB_1e_MHT_<<endl;
//   cout<<"eff_SB_sl_MHT_err_plus    "<<eff_SB_1e_MHT_err_[0]<<endl;
//   cout<<"eff_SB_sl_MHT_err_minus   "<<eff_SB_1e_MHT_err_[1]<<endl;

  cout<<"eff_SB_1e_MHT             "<<eff_SB_1e_MHT_<<endl;
  cout<<"eff_SB_1e_MHT_err_plus    "<<eff_SB_1e_MHT_err_[0]<<endl;
  cout<<"eff_SB_1e_MHT_err_minus   "<<eff_SB_1e_MHT_err_[1]<<endl;
  cout<<"eff_SB_1m_MHT             "<<eff_SB_1m_MHT_<<endl;
  cout<<"eff_SB_1m_MHT_err_plus    "<<eff_SB_1m_MHT_err_[0]<<endl;
  cout<<"eff_SB_1m_MHT_err_minus   "<<eff_SB_1m_MHT_err_[1]<<endl;

  cout<<"eff_SIG_MHT               "<<eff_SIG_MHT_<<endl;
  cout<<"eff_SIG_MHT_err_plus      "<<eff_SIG_MHT_err_[0]<<endl;
  cout<<"eff_SIG_MHT_err_minus     "<<eff_SIG_MHT_err_[1]<<endl;

  cout<<"eff_SIG_ldp_MHT           "<<eff_SIG_ldp_MHT_<<endl;
  cout<<"eff_SIG_ldp_MHT_err_plus  "<<eff_SIG_ldp_MHT_err_[0]<<endl;
  cout<<"eff_SIG_ldp_MHT_err_minus "<<eff_SIG_ldp_MHT_err_[1]<<endl;

  cout<<"eff_SIG_sl_MHT            "<<eff_SIG_SL_MHT_<<endl;
  cout<<"eff_SIG_sl_MHT_err_plus   "<<eff_SIG_SL_MHT_err_[0]<<endl;
  cout<<"eff_SIG_sl_MHT_err_minus  "<<eff_SIG_SL_MHT_err_[1]<<endl;

}

float leg_x1 , leg_x2, leg_y1, leg_y2;
void resetLegendPosition() {
  leg_x1 = 0.696;
  leg_x2=0.94;
  leg_y1=0.5;
  leg_y2=0.9;
}
void resetLegendPositionR() {
  leg_x1 = 0.2;
  leg_x2=0.45;
  leg_y1=0.65;
  leg_y2=0.9;
}


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
ConfigurationDescriptions::ConfigurationDescriptions() : 
  default_(""),corrected_("") { }
ConfigurationDescriptions::~ConfigurationDescriptions() {}

TString ConfigurationDescriptions::at(const unsigned int i) {
  if (i == 0) return getDefault();
  else  if (i == 1) return getCorrected();
  else {
    const unsigned int j=i-2;
    unsigned int k=0;
    for (std::map<TString, std::pair<TString,TString> >::iterator iconfig=variationPairs.begin(); iconfig!=variationPairs.end(); ++iconfig) {
      if (j == k) return iconfig->second.first;
      else if (j == k+1) return iconfig->second.second;
      k+=2;
    }
  }

  cout<<"WARNING in ConfigurationDescriptions::at() -- asked for element "<<i<<" when there are only "<<this->size()<<" elements!"<<endl;
  return "";
}

unsigned int ConfigurationDescriptions::size() {
  unsigned int s=0;
  if (default_ != "") ++s;
  if (corrected_ != "") ++s;

  s+= 2*variationPairs.size();

  return s;
}

TString ConfigurationDescriptions::getVariedSubstring(const TString & currentVariation) {

  TObjArray* baseline = corrected_.Tokenize("_");

  TObjArray* mine = currentVariation.Tokenize("_");

  TString output="";
  for (int i=0; i<baseline->GetEntries(); i++) {
    TString b=baseline->At(i)->GetName();
    TString m=mine->At(i)->GetName();
    if (b!=m) {
      output+=b;
    }
  }
  return output;
}

void ConfigurationDescriptions::addVariation(const TString & description1, const TString & description2) {

  TString var1=getVariedSubstring(description1);
  TString var2=getVariedSubstring(description2);

  assert(var1 == var2);

  variationPairs[var1] = make_pair(description1,description2);

}

ConfigurationDescriptions configDescriptions_;

//special containers for holding the "search regions of interest"
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
SearchRegion::SearchRegion(TString btagSel,TString htSel,TString metSel,TString oId,bool isSig) : 
  htSelection(htSel),metSelection(metSel),btagSelection(btagSel),owenId(oId),isSIG(isSig) {}
SearchRegion::~SearchRegion() {}
void SearchRegion::Print() const {
  cout<<" == "<<btagSelection<<" "<<htSelection<<" "<<metSelection<<endl;

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

std::vector<SearchRegion > searchRegions_;
std::vector<SearchRegion > sbRegions_;
bool searchRegionsSet_=false;

void setSearchRegions( TString  which="") {
  //note that the search regions can be set exactly once per session. after that they will not be overridden by
  //further calls to this function
  //this is intentional...since this code is ill-designed and there are zillions of calls to setSearchRegions()
  //scattered throughout the code. This way the user is ensured that the first call has primacy
  if (searchRegionsSet_) return;
  
  //nb: some of the code *depends* on the fact that for there are equal numbers of corresponding
  //sbRegions and searchRegions, with the only difference being the MET selection!

  //also, for style reasons the 'owenId' should not contain the number of b tags.
  //everywhere that we use the owenId as an identifier, we combine with the number of b tags

  // i honestly can't remember if the 'owenId' must match between SB and SIG regions

  if (which=="") {cout<<"Setting default 'Moriond Wide SB 2012' search regions"<<endl; which="MoriondWideSB";}
  else    cout<<"Setting search regions to set: "<<which<<endl;


  //27 Jan 2012 -- (preliminary) regions for testing shapes of background....
  if (which=="METbins3B") {
    setTrigEff("WideSB"); //set trigger efficiency
    sbRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=150&&MET<250","METBin1",false));
    searchRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=250 &&MET<300","METBin1"));
    
    sbRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=150&&MET<250","METBin2",false));
    searchRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=300 &&MET<350","METBin2"));
    
    sbRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=150&&MET<250","METBin3",false));
    searchRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=350","METBin3"));
  }
  else if (which=="METfinebins1BL") {
    setTrigEff("WideSB"); //set trigger efficiency
    sbRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=150&&MET<250","METFineBin1",false));
    searchRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=250 &&MET<300","METFineBin1"));
    
    sbRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=150&&MET<250","METFineBin2",false));
    searchRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=300 &&MET<350","METFineBin2"));
    
    sbRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=150&&MET<250","METFineBin3",false));
    searchRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=350 &&MET<400","METFineBin3"));

    sbRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=150&&MET<250","METFineBin4",false));
    searchRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=400 &&MET<450","METFineBin4"));

    sbRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=150&&MET<250","METFineBin5",false));
    searchRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=450 &&MET<500","METFineBin5"));

    sbRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=150&&MET<250","METFineBin6",false));
    searchRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=500","METFineBin6"));
  }
  else if (which=="METfinebins2BL") {
    setTrigEff("WideSB"); //set trigger efficiency
    sbRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=150&&MET<250","METFineBin1",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=250 &&MET<300","METFineBin1"));
    
    sbRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=150&&MET<250","METFineBin2",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=300 &&MET<350","METFineBin2"));
    
    sbRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=150&&MET<250","METFineBin3",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=350 &&MET<400","METFineBin3"));

    sbRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=150&&MET<250","METFineBin4",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=400 &&MET<450","METFineBin4"));

    sbRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=150&&MET<250","METFineBin5",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=450 &&MET<500","METFineBin5"));

    sbRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=150&&MET<250","METFineBin6",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=500","METFineBin6"));
  }
  else if (which=="METfinebins2BT") {
    setTrigEff("WideSB"); //set trigger efficiency
    
    sbRegions_.push_back( SearchRegion( "ge2b","HT>=600","MET>=150&&MET<250","METFineBin1",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=600","MET>=300 &&MET<350","METFineBin1"));
    
    sbRegions_.push_back( SearchRegion( "ge2b","HT>=600","MET>=150&&MET<250","METFineBin2",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=600","MET>=350 &&MET<400","METFineBin2"));

    sbRegions_.push_back( SearchRegion( "ge2b","HT>=600","MET>=150&&MET<250","METFineBin3",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=600","MET>=400 &&MET<450","METFineBin3"));

    sbRegions_.push_back( SearchRegion( "ge2b","HT>=600","MET>=150&&MET<250","METFineBin4",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=600","MET>=450 &&MET<500","METFineBin4"));

    sbRegions_.push_back( SearchRegion( "ge2b","HT>=600","MET>=150&&MET<250","METFineBin5",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=600","MET>=500","METFineBin5"));
  }
  else if (which=="METfinebins3B") {
    setTrigEff("WideSB"); //set trigger efficiency
    sbRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=150&&MET<250","METFineBin1",false));
    searchRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=250 &&MET<300","METFineBin1"));
    
    sbRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=150&&MET<250","METFineBin2",false));
    searchRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=300 &&MET<350","METFineBin2"));
    
    sbRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=150&&MET<250","METFineBin3",false));
    searchRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=350 &&MET<400","METFineBin3"));

    sbRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=150&&MET<250","METFineBin4",false));
    searchRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=400 &&MET<450","METFineBin4"));

    sbRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=150&&MET<250","METFineBin5",false));
    searchRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=450 &&MET<500","METFineBin5"));

    sbRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=150&&MET<250","METFineBin6",false));
    searchRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=500","METFineBin6"));
  }
//These are the nominal search regions for the "Moriond 2012" analysis
  else if (which=="Moriond") {
  //oct25
    setTrigEff("NominalSB"); //set trigger efficiency
    
    sbRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=200&&MET<250","Loose",false));
    searchRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=250","Loose")); //1BL
    
    sbRegions_.push_back( SearchRegion( "ge1b","HT>=500","MET>=200&&MET<250","Tight",false));
    searchRegions_.push_back( SearchRegion( "ge1b","HT>=500","MET>=500","Tight")); //1BT
    
    sbRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=200&&MET<250","Loose",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=250","Loose")); //2BL
  
    sbRegions_.push_back( SearchRegion( "ge2b","HT>=600","MET>=200&&MET<250","Tight",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=600","MET>=300","Tight")); //2BT
    
    sbRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=200&&MET<250","Loose",false));
    searchRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=250","Loose")); //3B
  }
  else if (which=="MoriondLowSB") {
    setTrigEff("LowSB"); //set trigger efficiency
    //move SB down to 150-200
    sbRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=150&&MET<200","LooseLowSB",false));
    searchRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=250","LooseLowSB")); //1BL
    
    sbRegions_.push_back( SearchRegion( "ge1b","HT>=500","MET>=150&&MET<200","TightLowSB",false));
    searchRegions_.push_back( SearchRegion( "ge1b","HT>=500","MET>=500","TightLowSB")); //1BT
    
    sbRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=150&&MET<200","LooseLowSB",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=250","LooseLowSB")); //2BL
    
    sbRegions_.push_back( SearchRegion( "ge2b","HT>=600","MET>=150&&MET<200","TightLowSB",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=600","MET>=300","TightLowSB")); //2BT
    
    sbRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=150&&MET<200","LooseLowSB",false));
    searchRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=250","LooseLowSB")); //3B
  }
  else if (which=="MoriondWideSB") {
    setTrigEff("WideSB"); //set trigger efficiency
    //SB is 150-250
    sbRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=150&&MET<250","LooseWideSB",false));
    searchRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=250","LooseWideSB")); //1BL
        
    sbRegions_.push_back( SearchRegion( "ge1b","HT>=500","MET>=150&&MET<250","TightWideSB",false));
    searchRegions_.push_back( SearchRegion( "ge1b","HT>=500","MET>=500","TightWideSB")); //1BT
    
    sbRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=150&&MET<250","LooseWideSB",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=400","MET>=250","LooseWideSB")); //2BL
    
    sbRegions_.push_back( SearchRegion( "ge2b","HT>=600","MET>=150&&MET<250","TightWideSB",false));
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=600","MET>=300","TightWideSB")); //2BT
    
    sbRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=150&&MET<250","LooseWideSB",false));
    searchRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=250","LooseWideSB")); //3B
    
  }
  else if (which=="eq2BT") {
    //work in progress
    setTrigEff("WideSB"); //set trigger efficiency
    sbRegions_.push_back( SearchRegion( "eq2b","HT>=600","MET>=150&&MET<250","TightHTWideSB",false));
    searchRegions_.push_back( SearchRegion( "eq2b","HT>=600","MET>=300","TightHTWideSB")); //2BT
  }
  else if (which=="SpecialTest1") {
    setTrigEff("LowSBSpecialTest1"); //set trigger efficiency
    sbRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=150&&MET<200","LooseLowSB",false));
    searchRegions_.push_back( SearchRegion( "ge1b","HT>=400","MET>=210 && MET<230","LooseLowSB")); //1BL

  }
  //important note -- i am writing the signal systematics code to _assume_  that the SB region is shared for the shape analysis.
  //so don't try to combine the different sets of HT cuts into one set of regions
  else if (which=="bShapeLoose") {
    setTrigEff("NominalSB"); //set trigger efficiency
    //Tentative search regions for the b-shape analysis
    //Loose exclusive regions
    sbRegions_.push_back( SearchRegion( "eq1b","HT>=400","MET>=200&&MET<250","Loose",false));
    searchRegions_.push_back( SearchRegion( "eq1b","HT>=400","MET>=250","Loose")); //1BL
    
    sbRegions_.push_back( SearchRegion( "eq2b","HT>=400","MET>=200&&MET<250","Loose",false));
    searchRegions_.push_back( SearchRegion( "eq2b","HT>=400","MET>=250","Loose")); //2BL
    
    sbRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=200&&MET<250","Loose",false));
    searchRegions_.push_back( SearchRegion( "ge3b","HT>=400","MET>=250","Loose")); //3B
  }
  else if (which=="bShapeTightHT") {
    setTrigEff("NominalSB"); //set trigger efficiency
    //2BT exclusive regions
    sbRegions_.push_back( SearchRegion( "eq1b","HT>=600","MET>=200&&MET<250","TightHT",false));
    searchRegions_.push_back( SearchRegion( "eq1b","HT>=600","MET>=300","TightHT")); //2BT
    
    sbRegions_.push_back( SearchRegion( "eq2b","HT>=600","MET>=200&&MET<250","TightHT",false));
    searchRegions_.push_back( SearchRegion( "eq2b","HT>=600","MET>=300","TightHT")); //2BT
    
    sbRegions_.push_back( SearchRegion( "ge3b","HT>=600","MET>=200&&MET<250","TightHT",false));
    searchRegions_.push_back( SearchRegion( "ge3b","HT>=600","MET>=300","TightHT")); //2BT
  }
  else if (which=="bShapeTightMET") {
    setTrigEff("NominalSB"); //set trigger efficiency
    //1BT exclusive regions
    sbRegions_.push_back( SearchRegion( "eq1b","HT>=500","MET>=200&&MET<250","TightMET",false));
    searchRegions_.push_back( SearchRegion( "eq1b","HT>=500","MET>=500","TightMET")); //1BT
    
    sbRegions_.push_back( SearchRegion( "eq2b","HT>=500","MET>=200&&MET<250","TightMET",false));
    searchRegions_.push_back( SearchRegion( "eq2b","HT>=500","MET>=500","TightMET")); //1BT
    
    sbRegions_.push_back( SearchRegion( "ge3b","HT>=500","MET>=200&&MET<250","TightMET",false));
    searchRegions_.push_back( SearchRegion( "ge3b","HT>=500","MET>=500","TightMET")); //1BT
  }
  else if (which=="Summer2011") {
    setTrigEff("LowSB"); //set trigger efficiency
    //2011 Summer result
    sbRegions_.push_back( SearchRegion( "ge1b","HT>=350","MET>=150&&MET<200","Loose",false)); //loose SB
    searchRegions_.push_back( SearchRegion( "ge1b","HT>=350","MET>=200","Loose")); //loose Sig
    
    sbRegions_.push_back( SearchRegion( "ge1b","HT>=500","MET>=150&&MET<200","Tight",false)); //tight SB
    searchRegions_.push_back( SearchRegion( "ge1b","HT>=500","MET>=300","Tight")); //tight Sig
    
    sbRegions_.push_back( SearchRegion( "ge2b","HT>=350","MET>=150&&MET<200","Loose",false)); //loose SB
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=350","MET>=200","Loose")); //loose Sig
    
    sbRegions_.push_back( SearchRegion( "ge2b","HT>=500","MET>=150&&MET<200","Tight",false)); //tight SB
    searchRegions_.push_back( SearchRegion( "ge2b","HT>=500","MET>=300","Tight")); //tight Sig
  }
  else assert(0);

  searchRegionsSet_=true;
}

//very simple data container used by the SignalEffData class
class SystInfo {
public:
  SystInfo(float p=0, float m=0, int s=0);
  ~SystInfo();
  SystInfo( ifstream* input);

  float plus;
  float minus;
  //convention:
  //use 0 for numbers that are completely unset
  //use 1 for numbers that are fixed by SignalEffData
  //use 2 for numbers that are set by some actual event counts
  int status;

  void write(ofstream* outfile) const;
};

void SystInfo::write(ofstream* outfile) const {

  (*outfile)<<status<<" "<<minus<<" "<<plus<<endl;

}

SystInfo::SystInfo(ifstream* input) :
  plus(0),minus(0),status(99)
{

  (*input)>>status>>minus>>plus;
  cout<< "  Loaded status,minus,plus = "<<status<<" "<<minus<<" "<<plus<<endl; //DEBUG
}

SystInfo::SystInfo(float p, float m, int s) :
  plus(p),  minus(m),  status(s) {}
SystInfo::~SystInfo() {}

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
  map<TString, SystInfo > systematics;
  //motivation for using a map:
  //  want to be able to clear() it and *regain* the memory footprint, without getting rid of the class object itself
  TString translateVariation(const TString & which) ;

  void write(TString id) const; //write the SignalEffData contents to a file

};

const TString SignalEffDataSuffix_ = "sedf";
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


struct OwenData {
  double Nsig; //number in signal region , data //done
  double Nsb; // number in SB, data             //done
  //  double Nsig_sl; //number in SL SIG, data  //done
  //  double Nsb_sl; // number in SL SB, data   //done

  double Nsig_sl; //number in SL SIG, data  //done
  double Nsb_1e; // number in SL SB, data   //done
  double Nsb_1m; // number in SL SB, data   //done

  double Nsig_ldp; //number in SIG, fail DP //done
  double Nsb_ldp;   // number in SB, fail DP //done

  //owen didn't ask for these
  //  double  Nlsb ;
  //  double  Nlsb_ldp;

  //deprecated
  //  double  Nlsb_0b ;      // done
  //  double  Nlsb_0b_ldp;   // done
  double Rlsb_passfail;
  double Rlsb_passfail_err;


  double Nttbarmc_sig_ldp; //done
  double Nttbarmc_sb_ldp; //done
  //  double lsf_WJmc; //done
  double NWJmc_sig_ldp; //done
  double NWJmc_sb_ldp; //done
  //  double lsf_Znnmc; //done
  double NZnnmc_sig_ldp; //done
  double NZnnmc_sb_ldp; //done
  //  double lsf_Zjmc;
  double NZjmc_sig_ldp;
  double NZjmc_sb_ldp;
  double Nsingletopmc_sig_ldp;
  double Nsingletopmc_sb_ldp;

  //don't need DataLumi...that's just lumiScale_
} ;



struct texData {

  double value;
  double statError;
  double systError;
  double trigErrorPlus;
  double trigErrorMinus;

};

//  btag+owenId  background cat
map<TString, map<TString, texData> > resultsMap_;

//again, this should be C++ style rather than C style, but let's just do it like this for now
//in other words, resultsMap_ should be an instance of a class that has a ::writeToText() function...
void writeResultsMapToText(TString filename) {

  ofstream outfile(filename.Data());

  for (map<TString, map<TString, texData> >::iterator ibin = resultsMap_.begin() ; ibin!=resultsMap_.end() ; ++ibin) {
    map<TString, texData> thisbin = ibin->second;
    for ( map<TString, texData> ::iterator ibackground = thisbin.begin(); ibackground!= thisbin.end(); ++ibackground) {
      outfile<<ibin->first<<" "<<ibackground->first<<" "
	     <<ibackground->second.value<<" "
	     <<ibackground->second.statError<<" "
	     <<ibackground->second.systError<<" "
	     <<ibackground->second.trigErrorPlus<<" "
	     <<ibackground->second.trigErrorMinus<<" "
	     <<endl;
    }
  }

}

void addResults(TString id, TString addcat1, TString addcat2, TString newcat) {

  resultsMap_[id][newcat].value = resultsMap_[id][addcat1].value  + resultsMap_[id][addcat2].value ;
  resultsMap_[id][newcat].statError = sqrt(pow(resultsMap_[id][addcat1].statError,2)  + pow(resultsMap_[id][addcat2].statError,2));
  resultsMap_[id][newcat].systError = sqrt(pow(resultsMap_[id][addcat1].systError,2)  + pow(resultsMap_[id][addcat2].systError,2));

  const float tol=0.001;

  bool hastrig1= resultsMap_[id][addcat1].trigErrorPlus >tol || resultsMap_[id][addcat1].trigErrorMinus >tol;
  bool hastrig2= resultsMap_[id][addcat2].trigErrorPlus >tol || resultsMap_[id][addcat2].trigErrorMinus >tol;
  //it's not ok to add trigger errors because they will be correlated
  assert( !(hastrig1 && hastrig2));

  if (hastrig1) {
    resultsMap_[id][newcat].trigErrorPlus = resultsMap_[id][addcat1].trigErrorPlus;
    resultsMap_[id][newcat].trigErrorMinus = resultsMap_[id][addcat1].trigErrorMinus;
  }
  else if (hastrig2) {
    resultsMap_[id][newcat].trigErrorPlus = resultsMap_[id][addcat2].trigErrorPlus;
    resultsMap_[id][newcat].trigErrorMinus = resultsMap_[id][addcat2].trigErrorMinus;
  }
  else {
    resultsMap_[id][newcat].trigErrorPlus = 0;
    resultsMap_[id][newcat].trigErrorMinus = 0;
  }

}

TString formatLatex(const texData & data) {
  //in principle we should adapt format_nevents() but for now do this

  const float tol=0.001;


  char out[50];
  if (data.trigErrorPlus <tol && data.trigErrorMinus<tol)  sprintf(out,"$%.1f \\pm %.1f \\pm %.1f$",data.value,data.statError,data.systError);
  else {
    sprintf(out,"$%.1f \\pm %.1f^{+%.1f}_{-%.1f}$",data.value,data.statError,
	    sqrt(data.systError*data.systError + data.trigErrorPlus*data.trigErrorPlus),
	    sqrt(data.systError*data.systError + data.trigErrorMinus*data.trigErrorMinus));
  }

  return TString(out);
}



//ok, i would rather implement this as a C++ style class instead of a C style function, but for now I'll do the easiest thing
void writeSignalEffDataMapToFiles(const  map<pair<int,int>, SignalEffData> & thedata, const TString & id ) {

  //goal -- write thedata to a set of files
  //SignalEffData already has a method to write its contents to a file

  for ( map<pair<int,int>, SignalEffData>::const_iterator iscan = thedata.begin(); iscan!= thedata.end(); ++iscan) {

    TString fileid;
    fileid.Form("%s_%d_%d",id.Data(), iscan->first.first, iscan->first.second);

    iscan->second.write(fileid);
  }
}
map<pair<int,int>, SignalEffData> loadSignalEffDataMapFromFiles( const TString & id ) {

  map<pair<int,int>, SignalEffData> datamap;

  //this is more tricky. we have to loop over all files that contain (not begin with) id
  TChain dummychain("dummychain");
  TString searchfor ="*.";
  searchfor += SignalEffDataSuffix_;
  searchfor.Prepend(sedfPath);
  dummychain.Add(searchfor);
  TObjArray* dir1list = dummychain.GetListOfFiles();
  int nfiles1=dir1list->GetEntries();
  
  for (int ifile1=0; ifile1<nfiles1; ifile1++) {
    TString fullname=dir1list->At(ifile1)->GetTitle();
    if (fullname.Contains(id)) {
      TString justMasses = fullname(fullname.Index(id)+id.Length()+1,fullname.Length());
      //now extract the mass values
      int m1=TString(justMasses(0,justMasses.Index("_"))).Atoi();
      int m2=TString(justMasses(justMasses.Index("_")+1,justMasses.Length() )).Atoi();
      justMasses.Prepend("_");
      justMasses.Prepend(id);
      SignalEffData data(justMasses);
      datamap[make_pair(m1,m2)] = data;
    }
  }

  return datamap;
}


std::map<TString, OwenData> owenMap_;

void printOwen(const TString& owenKey) {

  cout<< " === "<<owenKey<<" === "<<endl;

  cout<<"Nsig              "<<  owenMap_[owenKey].Nsig<<endl;
  cout<<"Nsb               "<<  owenMap_[owenKey].Nsb<<endl;

  cout<<"Nsig_sl           "<<  owenMap_[owenKey].Nsig_sl<<endl;
  //  cout<<"Nsb_sl            "<<  owenMap_[owenKey].Nsb_sl<<endl;
  cout<<"Nsb_1e            "<<  owenMap_[owenKey].Nsb_1e<<endl;
  cout<<"Nsb_1m            "<<  owenMap_[owenKey].Nsb_1m<<endl;

  cout<<"Nsig_ldp          "<<  owenMap_[owenKey].Nsig_ldp<<endl;
  cout<<"Nsb_ldp           "<<  owenMap_[owenKey].Nsb_ldp<<endl;

  //  cout<<"Nlsb_0b           "<<  owenMap_[owenKey].Nlsb_0b<<endl;
  //  cout<<"Nlsb_0b_ldp       "<<  owenMap_[owenKey].Nlsb_0b_ldp<<endl;
  cout<<"Rlsb_passfail     "<<owenMap_[owenKey].Rlsb_passfail<<endl;
  cout<<"Rlsb_passfail_err "<<owenMap_[owenKey].Rlsb_passfail_err<<endl;

  cout<<"Nttbarmc_sig_ldp  "<<  owenMap_[owenKey].Nttbarmc_sig_ldp<<endl;
  cout<<"Nttbarmc_sb_ldp   "<<  owenMap_[owenKey].Nttbarmc_sb_ldp<<endl;

  cout<<"Nsingletopmc_sig_ldp  "<<  owenMap_[owenKey].Nsingletopmc_sig_ldp<<endl;
  cout<<"Nsingletopmc_sb_ldp   "<<  owenMap_[owenKey].Nsingletopmc_sb_ldp<<endl;

  //  cout<<"lsf_WJmc          "<<  owenMap_[owenKey].lsf_WJmc<<endl;
  cout<<"NWJmc_sig_ldp     "<<  owenMap_[owenKey].NWJmc_sig_ldp<<endl;
  cout<<"NWJmc_sb_ldp      "<<  owenMap_[owenKey].NWJmc_sb_ldp<<endl;

  //  cout<<"lsf_Znnmc         "<<  owenMap_[owenKey].lsf_Znnmc<<endl;
  cout<<"NZnnmc_sig_ldp    "<<  owenMap_[owenKey].NZnnmc_sig_ldp<<endl;
  cout<<"NZnnmc_sb_ldp     "<<  owenMap_[owenKey].NZnnmc_sb_ldp<<endl;

  //  cout<<"lsf_Zjmc         "<<  owenMap_[owenKey].lsf_Zjmc<<endl;
  cout<<"NZjmc_sig_ldp    "<<  owenMap_[owenKey].NZjmc_sig_ldp<<endl;
  cout<<"NZjmc_sb_ldp     "<<  owenMap_[owenKey].NZjmc_sb_ldp<<endl;

}
void printOwenShape(const TString& owenKey) {

  //  cout<< " === "<<owenKey<<" === "<<endl;

  TString bstr="";
  if (owenKey.BeginsWith("eq1b")) bstr="_1b";
  else  if (owenKey.BeginsWith("eq2b")) bstr="_2b";
  else  if (owenKey.BeginsWith("ge3b")) bstr="_3b";
  else assert(0);

  cout<<"Nsig"<<bstr<<"              "<<  owenMap_[owenKey].Nsig<<endl;
  cout<<"Nsb"<<bstr<<"               "<<  owenMap_[owenKey].Nsb<<endl;

  cout<<"Nsig_sl"<<bstr<<"           "<<  owenMap_[owenKey].Nsig_sl<<endl;
   cout<<"Nsb_1e"<<bstr<<"            "<<  owenMap_[owenKey].Nsb_1e<<endl;
   cout<<"Nsb_1m"<<bstr<<"            "<<  owenMap_[owenKey].Nsb_1m<<endl;

  cout<<"Nsb_sl"<<bstr<<"            "<<  owenMap_[owenKey].Nsb_1e+owenMap_[owenKey].Nsb_1m<<endl;


  cout<<"Nsig_ldp"<<bstr<<"          "<<  owenMap_[owenKey].Nsig_ldp<<endl;
  cout<<"Nsb_ldp"<<bstr<<"           "<<  owenMap_[owenKey].Nsb_ldp<<endl;

  if (owenKey.BeginsWith("eq1b")) {
    cout<<"Rlsb_passfail     "<<owenMap_[owenKey].Rlsb_passfail<<endl;
    cout<<"Rlsb_passfail_err "<<owenMap_[owenKey].Rlsb_passfail_err<<endl;
  }

  cout<<"Nttbarsingletopzjetsmc_sig_ldp"<<bstr<<"  "<<  owenMap_[owenKey].Nttbarmc_sig_ldp+owenMap_[owenKey].Nsingletopmc_sig_ldp+owenMap_[owenKey].NZjmc_sig_ldp<<endl;
  cout<<"Nttbarsingletopzjetsmc_sb_ldp"<<bstr<<"  "<<  owenMap_[owenKey].Nttbarmc_sb_ldp+owenMap_[owenKey].Nsingletopmc_sb_ldp+owenMap_[owenKey].NZjmc_sb_ldp<<endl;

//   cout<<"Nttbarmc_sig_ldp"<<bstr<<"  "<<  owenMap_[owenKey].Nttbarmc_sig_ldp<<endl;
//   cout<<"Nttbarmc_sb_ldp"<<bstr<<"   "<<  owenMap_[owenKey].Nttbarmc_sb_ldp<<endl;

//   cout<<"Nsingletopmc_sig_ldp"<<bstr<<"  "<<  owenMap_[owenKey].Nsingletopmc_sig_ldp<<endl;
//   cout<<"Nsingletopmc_sb_ldp"<<bstr<<"   "<<  owenMap_[owenKey].Nsingletopmc_sb_ldp<<endl;

  cout<<"NWJmc_sig_ldp"<<bstr<<"     "<<  owenMap_[owenKey].NWJmc_sig_ldp<<endl;
  cout<<"NWJmc_sb_ldp"<<bstr<<"      "<<  owenMap_[owenKey].NWJmc_sb_ldp<<endl;

  cout<<"NZnnmc_sig_ldp"<<bstr<<"    "<<  owenMap_[owenKey].NZnnmc_sig_ldp<<endl;
  cout<<"NZnnmc_sb_ldp"<<bstr<<"     "<<  owenMap_[owenKey].NZnnmc_sb_ldp<<endl;

//   cout<<"NZjmc_sig_ldp"<<bstr<<"    "<<  owenMap_[owenKey].NZjmc_sig_ldp<<endl;
//   cout<<"NZjmc_sb_ldp"<<bstr<<"     "<<  owenMap_[owenKey].NZjmc_sb_ldp<<endl;

}


// class for holding inputs to drawSimple (could be expanded)
class toPlot {
  TString varname;
  int nbins;
  double binlow, binhigh;
public:
  toPlot(TString,int,double,double);
  void print();
  TString getVarName() {return varname;}
  int getNBins() {return nbins;}
  double getBinLow() {return binlow;}
  double getBinHigh() {return binhigh;}
};
toPlot::toPlot(TString a, int b, double c, double d) {
  varname = a;
  nbins = b;
  binlow = c;
  binhigh = d;
}
void toPlot::print() {
  cout << varname << ": " << nbins << " bins, from " << binlow << " to " << binhigh << endl;
}


bool quiet_=false;
//bool quiet_=true;
bool doRatio_=false;
bool logy_=false;
bool dostack_=true;
bool doleg_=true;
bool dodata_=true;
bool addOverflow_=true;
//bool doSubtraction_=false;
bool drawMCErrors_=false;
bool renormalizeBins_=false;//no setter function
bool owenColor_ = false;

int m0_=0;
int m12_=0;
TString susyCrossSectionVariation_="";

bool normalized_=false;

bool useFlavorHistoryWeights_=false;//no setter function
float flavorHistoryScaling_=-1;

bool usePUweight_=false;
bool useHTeff_ = false;
bool useMHTeff_ = false;
enum bnnMHTeffMode {kOff=0, kOn, kOnPlus, kOnMinus, kPlot, kPlotPlus, kPlotMinus};
bnnMHTeffMode thebnnMHTeffMode_ = kOff;

TString btagSFweight_="1";

bool savePlots_ = true; //no setter function
bool drawTotalSM_=false; //no setter function
bool drawTotalSMSusy_=false;//no setter function
bool drawSusyOnly_=false;//no setter function
bool drawMarkers_=true;//no setter function

bool doVerticalLine_=false;
double verticalLinePosition_=0;

bool doCustomPlotMax_=false;
double customPlotMax_=0;

bool doCustomPlotMin_=false;
double customPlotMin_=0;

float maxScaleFactor_ = 1.05;

bool splitTTbar_ = false;
bool splitWJets_ = false;
//the next three are automatically configured in slABCD()
bool splitTTbarForClosureTest_ = false;
bool splitWJetsForClosureTest_ = false;
bool splitSingleTopForClosureTest_ = false;

//bool latexMode_=false;
bool latexMode_=true;
const TString pm = latexMode_ ? " \\pm " : " +/- ";

TCanvas* thecanvas=0;
//TCanvas* cratio=0;
TLegend* leg=0;
THStack* thestack=0;
TH1D* totalsm=0;
TH1D* totalsmsusy=0;
TH1D* totalewk=0;
TH1D* totalqcdttbar=0;
TH1D* totalnonttbar=0;
TH1D* totalnonqcd=0;
TH1D* totalqcd=0; //ben - just for ease of doing event counts with drawPlots
TH1D* ratio=0; float ratioMin=0; float ratioMax=2.5;
TLine* ratioLine=0;
TGraphErrors* mcerrors=0;
bool loaded_=false; //bookkeeping
bool loadedSusyHistos_=false;//bookkeeping

// == set configuration options ==
void setQuiet(bool q) {
  quiet_ = q;
}

void doRatioPlot(bool doIt) {
  doRatio_=doIt;
}

void setPlotMaximum(double max) {
  customPlotMax_=max;
  doCustomPlotMax_=true;
}

void resetPlotMaximum() {
  doCustomPlotMax_=false;
}

void setPlotMinimum(double min) {
  customPlotMin_=min;
  doCustomPlotMin_=true;
}

void resetPlotMinimum() {
  doCustomPlotMin_=false;
}

void enableVerticalLine(double position) {
  doVerticalLine_=true;
  verticalLinePosition_ =position;
}

// void showDataMinusMC(bool dosub) {
//   doSubtraction_=dosub;
// }

void resetVerticalLine() {
  doVerticalLine_=false;
}

void doOverflowAddition(bool doOv) {
  addOverflow_ = doOv;
}

void setLogY(bool dolog) {
  logy_=dolog;
  if (logy_) maxScaleFactor_=3;
  else maxScaleFactor_=1.05;
}

void setStackMode(bool dostack, bool normalized=false) {
  dostack_=dostack;
  normalized_=normalized;
}

void doData(bool dodata) {
  dodata_=dodata;
}

void drawLegend(bool doleg) {
  doleg_=doleg;
}

void setLumiScale(double lumiscale){
  lumiScale_ = lumiscale;
}


CrossSectionTable * CrossSectionTable_mSUGRAtanb40_=0;
void loadSusyCrossSections() {
  if (CrossSectionTable_mSUGRAtanb40_==0) {
    CrossSectionTable_mSUGRAtanb40_ = new  CrossSectionTable("NLOxsec_tanb40_10.txt");
  }

}

bool isSampleScan(const TString & name ) {
  if (name.Contains("T1bbbb")) return true;
  if (name.Contains("T2bb")) return true;
  if (name.Contains("T1tttt")) return true;
  if (name.Contains("T2tt")) return true;
  if (name.Contains("SUGRA")) return true;


  return false;
}

bool isSampleSM(const TString & name) {

  if (name.Contains("LM")) return false;

  if (name.Contains("SUGRA")) return false;

  if (name.Contains("T1bbbb")) return false;
  if (name.Contains("T2bb")) return false;
  if (name.Contains("T1tttt")) return false;
  if (name.Contains("T2tt")) return false;

  return true;
}

void loadScanSMSngen(const TString& sampleOfInterest) {
  if (scanSMSngen==0) scanSMSngen = (TH2D*) files_[currentConfig_][sampleOfInterest]->Get("scanSMSngen");
}

void loadReferenceCrossSections() {
  if (referenceCrossSectionGluino==0) {
    cout<<"Loading reference cross section file"<<endl;
    TFile fxs("referenceXSecs.root");
    if (fxs.IsZombie()) {cout<<"Problem with cross section file!"<<endl; assert(0);}
    TH1D* hxs = (TH1D*) fxs.Get("gluino");
    gROOT->cd();
    referenceCrossSectionGluino = (TH1D*) hxs->Clone("referenceCrossSectionGluino");
    fxs.Close();
  }

}

//m0, m12            //pdf set
map<pair<int,int>, map<TString, TH2D*> >  scanProcessTotalsMap;
void loadSusyScanHistograms() {
  if (loadedSusyHistos_) return;


  TFile* susyfile = 0;
  for (unsigned int isample=0; isample<samples_.size(); isample++) {
    if ( !isSampleSM(samples_[isample])) {
      susyfile =  files_[currentConfig_][samples_[isample]];
      cout<<" Loading SUSY scan histograms from file: "<<susyfile->GetName()<<endl;
    }
  }

  if (susyfile==0) {cout<<"did not find SUSY file in loadSusyScanHistograms!"<<endl; return;}
  TString fn=susyfile->GetName();
  bool isSMS = !fn.Contains("SUGRA");

  if (isSMS && scanSMSngen==0) assert(0);

  const   int ntotal = susyfile->GetListOfKeys()->GetEntries();
  cout<<"found this many keys: "<<ntotal<<endl;
  int nloaded=0;
  for (int i=0; i<ntotal; i++) {
    TString objname=  susyfile->GetListOfKeys()->At(i)->GetName();
    if (objname.BeginsWith("scanProcessTotals") ) {
      TString pdfset=   objname.Tokenize("_")->At(0)->GetName();
      pdfset.ReplaceAll("scanProcessTotals","");
      TString m0str=   objname.Tokenize("_")->At(1)->GetName();
      TString m12str=   objname.Tokenize("_")->At(2)->GetName();
      int m0 = m0str.Atoi();
      int m12 = m12str.Atoi();
      if (isSMS && TMath::Nint(scanSMSngen->GetBinContent(scanSMSngen->FindBin(m0,m12))) ==0) continue;
      scanProcessTotalsMap[make_pair(m0,m12)][pdfset] = (TH2D*) susyfile->Get(objname);
      nloaded++;
      if (nloaded%1000==0) cout<<nloaded<<" ; "<<i<<" / "<<ntotal<<endl;
    }
  }

  loadedSusyHistos_=true;
}

void drawPlotHeaderInside() {
  if (text1 != 0 ) delete text1;
  //text1 = new TLatex(3.570061,23.08044,"CMS Preliminary");
  text1 = new TLatex(3.570061,23.08044,"CMS Simulation");
  text1->SetNDC();
  text1->SetTextAlign(13);
  text1->SetX(0.5);
  text1->SetY(.85);
  text1->SetTextFont(42);
  text1->SetTextSizePixels(24);
  text1->Draw();
}

void drawPlotHeader(double xoffset = 0) {
  //  return;
  float ypos = 0.97;
  if(doRatio_) ypos=ypos+0.012;
  // i'm gonna leave this out for now
  if (text1 != 0 ) delete text1;
  //text1 = new TLatex(3.570061,23.08044,"CMS"); //no more preliminary!
  text1 = new TLatex(3.570061,23.08044,"CMS Preliminary"); 
  text1->SetNDC();
  text1->SetTextAlign(13);
  text1->SetX(0.68 + xoffset); //add 0.2 if you get rid of the "Preliminary"
  text1->SetY(ypos+0.007);
  text1->SetTextFont(42);
  text1->SetTextSizePixels(24);
  text1->SetTextSize(0.045); //copied from ben's code. maybe needs a switch so it is only used for AN2011()
  text1->Draw();

  if (normalized_ == false) {
    TString astring;
    //astring.Form("%.0f pb^{-1} at #sqrt{s} = 7 TeV",lumiScale_);
    //astring.Form("%.1f fb^{-1} at #sqrt{s} = 7 TeV",lumiScale_/1000.);
    astring.Form("L_{int} = %.2f fb^{-1}, #sqrt{s} = 7 TeV",lumiScale_/1000.);
    if(lumiScale_>33. && lumiScale_<34.) astring.Form("L_{int} = %.2f fb^{-1}, #sqrt{s} = 7 TeV", 4982.91/1000.);//hardcoded, but don't know what else to do for this...
    if (text2 != 0 ) delete text2;
    text2 = new TLatex(3.570061,23.08044,astring);
    text2->SetNDC();
    text2->SetTextAlign(13);
    text2->SetX(0.17);
    text2->SetY(ypos+0.015);
    text2->SetTextFont(42);
    text2->SetTextSizePixels(24);
    text2->SetTextSize(0.045); //copied from ben's code. maybe needs a switch so it is only used for AN2011()
    text2->Draw();
  }

}


int mainpadWidth; int mainpadHeight;
int ratiopadHeight = 150;
// TPad* mainPad=0;
// TPad* ratioPad=0;
void renewCanvas(const TString opt="") {
  if (thecanvas!=0) delete thecanvas;

  int canvasWidth = mainpadWidth;
  int canvasHeight = opt.Contains("ratio") ? mainpadHeight+ratiopadHeight : mainpadHeight;

  thecanvas= new TCanvas("thecanvas","the canvas",canvasWidth,canvasHeight);
  thecanvas->cd()->SetRightMargin(0.04);
  thecanvas->cd()->SetTopMargin(0.07); //test

  if (opt.Contains("ratio")) {
    //thecanvas->Divide(1,2);
    //const float padding=0.01; const float ydivide=0.2;

    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);

    Float_t small = 1e-5;
    thecanvas->Divide(1,2,small,small);
    const float padding=1e-5; const float ydivide=0.3;
    thecanvas->GetPad(1)->SetPad( padding, ydivide + padding, 1-padding, 1-padding);
    thecanvas->GetPad(2)->SetPad( padding, padding, 1-padding, ydivide-padding);
    thecanvas->GetPad(1)->SetTopMargin(0.06);
    thecanvas->GetPad(1)->SetRightMargin(.05);
    thecanvas->GetPad(2)->SetRightMargin(.05);
    thecanvas->GetPad(2)->SetBottomMargin(.4);
    thecanvas->GetPad(1)->Modified();
    thecanvas->GetPad(2)->Modified();
    thecanvas->cd(1);
    gPad->SetBottomMargin(small);
    gPad->Modified();

    if (!quiet_)  cout<< thecanvas->GetPad(1)->GetXlowNDC() <<"\t"
		      << thecanvas->GetPad(1)->GetWNDC() <<"\t"
		      << thecanvas->GetPad(1)->GetYlowNDC() <<"\t"
		      << thecanvas->GetPad(1)->GetHNDC() <<endl;
    if (logy_) thecanvas->GetPad(1)->SetLogy();
  }
  else { if (logy_) thecanvas->SetLogy(); }


  int cdarg = opt.Contains("ratio") ? 1 : 0;
  thecanvas->cd(cdarg);

}

void resetPadDimensions() {
  mainpadWidth = 600; 
  mainpadHeight=550;
}

void setPadDimensions(int x, int y) {

  mainpadWidth = x; 
  mainpadHeight= y;
}

void renewLegend() {

  if (leg!=0) delete leg;
  leg = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
  leg->SetBorderSize(0);
  leg->SetLineStyle(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);

}

void resetHistos() {
  for ( std::map<TString, TH1D*>::iterator i = histos_.begin(); i!=histos_.end(); ++i) {
    if (i->second != 0) {
      delete  i->second;
      i->second= 0;
    }
  }
}

double findOverallMax(const TH1D* hh) {

  double max=-1e9;

  for (int i=1; i<= hh->GetNbinsX(); i++) {
    double val = hh->GetBinContent(i) + hh->GetBinError(i);
    if (val>max) max=val;
  }
  return max;
}

//code largely lifted from Owen
//returned string is the y title of the renormalized histo
TString renormBins( TH1D* hp, int refbin ) {

  if ( hp==0 ) return "PROBLEM";

  double refbinwid = -99;
  if(refbin==-1) refbinwid = 1.;
  else {
    refbin = hp->GetBinLowEdge( refbin+1 ) - hp->GetBinLowEdge( refbin ) ;
    if (!quiet_)  printf(" reference bin: [%6.1f,%6.1f], width = %6.3f\n",  hp->GetBinLowEdge( refbin ), hp->GetBinLowEdge( refbin+1 ), refbinwid ) ;
  }

  for ( int bi=1; bi<= hp->GetNbinsX(); bi++ ) {
    double binwid = hp->GetBinLowEdge( bi+1 ) - hp->GetBinLowEdge( bi ) ;
    double sf = refbinwid / binwid ;
    if (!quiet_)    printf("  bin %d : width= %6.2f, sf=%7.3f\n", bi, binwid, sf ) ;
    hp->SetBinContent( bi, sf*(hp->GetBinContent( bi )) ) ;
    hp->SetBinError( bi, sf*(hp->GetBinError( bi )) ) ;
  } // bi.

  TString ytitle;
  ytitle.Form("(Events / bin) * (%5.1f / bin width)", refbinwid );

  return ytitle;
}

void fillFlavorHistoryScaling() {

  TTree* tree=0;
  for (unsigned int isample=0; isample<samples_.size(); isample++) {
    if ( samples_[isample].Contains("WJets")) { //will pick out WJets or WJetsZ2
      tree = (TTree*) files_[currentConfig_][samples_[isample]]->Get("reducedTree");
    }
  }
  if (tree==0) {cout<<"Did not find a WJets sample!"<<endl; return;}

  gROOT->cd();
  TH1D dummyU("dummyU","",1,0,1e9);
  TH1D dummyk("dummyk","",1,0,1e9);
  assert(0);//please check that the following cuts have the weights and cuts handled correctly
  tree->Draw("HT>>dummyU","1","goff");
  tree->Draw("HT>>dummyk","flavorHistoryWeight","goff");
  flavorHistoryScaling_ = dummyU.Integral() / dummyk.Integral();
  if (!quiet_) cout<<"flavor history scaling factor = "<<flavorHistoryScaling_<<endl;

}

//try to bring some rationality to the options passed to getcutstring....
enum sampleType {kData, kMC, kmSugraPoint, kmSugraPlane, kSMSPoint, kSMSPlane};
//for more rationality, should make a data struct (or class) to hold the options, so that fewer arguments need to be passed
//TString getCutString(double lumiscale, TString extraWeight="", TString thisSelection="", TString extraSelection="", int pdfWeightIndex=0,TString pdfSet="CTEQ", bool isSusyScan=false, int susySubProcess=-1, const bool isData=false) {
TString getCutString(sampleType type, TString extraWeight="", TString thisSelection="", TString extraSelection="", int pdfWeightIndex=0,TString pdfSet="CTEQ", int susySubProcess=-1,const float scalefactor=1) {

  //an explanation of the options
  /*
kmSugraPoint means put a cut on m0,m12. Good for making plots of one sample point (e.g. MET distribution of 
m0=x, m12=y). But doesn't work for looking at the whole plane at once as we do when we do systematics. So
for that we have kmSugraPlane, where no requirement is placed on m0 and m12, only on the susySubProcess

Similarly, kSMSPoint will select one point in the mGL, mLSP plane, and weight using the reference cross section.
kSMSPlane treats the MC like data -- no weighting at all. Good for making efficiency measurements. (Calling function must 
normalize properly by Ngenerated)

for legacy purposes I am keeping all of the weight and selection TStrings, although I wonder if they are really needed....
  */
  double lumiscale = 0;

  //for now treat sms point just like smsplane

  if (type==kData)    lumiscale=1;
  else if (type==kMC || type==kmSugraPoint || type==kmSugraPlane ||type==kSMSPoint) lumiscale=lumiScale_;
  else if (type==kSMSPlane) lumiscale=1;
  else {assert(0);}

  TString weightedcut="weight"; 
  
  weightedcut += "*(";
  weightedcut +=lumiscale;
  weightedcut+=")";

  //new scale factor
  weightedcut += "*(";
  weightedcut +=scalefactor;
  weightedcut+=")";

  //horrible kludge
  //  if (type==kData) weightedcut+="*(runNumber>=178411)";

  //this flavorHistoryWeight business is too kludgey...someday should fix it
  if (extraWeight=="flavorHistoryWeight" && type!=kData) {
    if (flavorHistoryScaling_ <0) {
      fillFlavorHistoryScaling();
    }
    extraWeight.Form("flavorHistoryWeight*%f",flavorHistoryScaling_);
  }
  if (extraWeight!="") {
    weightedcut += "*(";
    weightedcut +=extraWeight;
    weightedcut+=")";
  }
  if (pdfWeightIndex != 0 && type!=kData) {
    TString pdfString;
    pdfString.Form("*pdfWeights%s[%d]",pdfSet.Data(),pdfWeightIndex);
    weightedcut += pdfString;
  }
  if (usePUweight_ && type!=kData) {
    weightedcut +="*PUweight";
  }
  if (useHTeff_ &&  type!=kData) {
    weightedcut +="*hltHTeff";
  }
  if (useMHTeff_ &&  type!=kData) {
    weightedcut +="*hltMHTeff";
  }
  if (thebnnMHTeffMode_==kOn &&  type==kData) {
    weightedcut +="*(1/hltMHTeffBNN)"; 
  }
  if (thebnnMHTeffMode_==kPlot &&  type!=kData) {
    weightedcut +="*hltMHTeffBNN"; 
  }
  if (thebnnMHTeffMode_==kPlotPlus &&  type!=kData) {
    weightedcut +="*hltMHTeffBNNUp"; 
  }
  if (thebnnMHTeffMode_==kPlotMinus &&  type!=kData) {
    weightedcut +="*hltMHTeffBNNDown"; 
  }
  if (thebnnMHTeffMode_==kOnPlus &&  type==kData) {
    weightedcut +="*(1/hltMHTeffBNNUp)"; 
  }
  if (thebnnMHTeffMode_==kOnMinus &&  type==kData) {
    weightedcut +="*(1/hltMHTeffBNNDown)"; 
  }
  if (btagSFweight_=="") btagSFweight_="1";
  if ( type==kData) {

    TString btagSFweight = btagSFweight_;
    //need to chop off the trailing _HFplus _LFminus etc
    if (btagSFweight.Contains("_")) btagSFweight = btagSFweight(0,btagSFweight.Index("_"));

    if (btagSFweight=="probge1") weightedcut += "*(nbjets>=1)";
    else if (btagSFweight=="probge2") weightedcut += "*(nbjets>=2)";
    else if (btagSFweight=="probge3") weightedcut += "*(nbjets>=3)";
    else if (btagSFweight=="prob2") weightedcut += "*(nbjets==2)";
    else if (btagSFweight=="(1-prob1-prob0-probge3)") weightedcut += "*(nbjets==2)";
    else if (btagSFweight=="prob1") weightedcut += "*(nbjets==1)";
    else if (btagSFweight=="prob0") weightedcut += "*(nbjets==0)";
    else if (btagSFweight=="1") {} //do nothing
    else {cout<<"btagSF problem with data"<<endl; assert(0);}
  }
  else {
    weightedcut += "*";
    weightedcut += btagSFweight_;
  } 

  if (thisSelection!="") {
    weightedcut += "*(";
    weightedcut+=thisSelection;
    if (extraSelection != "") {
      weightedcut += " && ";
      weightedcut +=extraSelection;
    }
    if (type == kmSugraPoint || type==kSMSPoint) {
      weightedcut += " && m0==";
      weightedcut +=m0_;
      weightedcut += " && m12==";
      weightedcut +=m12_;
    }
    weightedcut+=")";
  }
  else if (extraSelection !="") {
    weightedcut += "*(";
    weightedcut +=extraSelection;
    if (type == kmSugraPoint || type==kSMSPoint) {
      weightedcut += " && m0==";
      weightedcut +=m0_;
      weightedcut += " && m12==";
      weightedcut +=m12_;
    }
    weightedcut+=")";
  }

  //the easy part is done above -- keep only events at the current scan point
  //the hard part is to properly weight the events according to NLO cross sections
  if (type == kmSugraPoint) { 
    //the cross section for the subprocess is stored in scanCrossSection (systematics in scanCrossSectionPlus and scanCrossSectionMinus)
    //the normalization for the subprocess is stored in TH1D scanProcessTotals_<m0>_<m12>
    loadSusyScanHistograms();
    TH2D* thishist = scanProcessTotalsMap[make_pair(m0_,m12_)][pdfSet];
    if (thishist==0) cout<<"We've got a problem in getCutString!"<<endl;
    TString susyprocessweight = "(";
    int lowbound=0; int highbound=10;
    if (susySubProcess>=0) { lowbound=susySubProcess; highbound=susySubProcess;}
    for (int i=lowbound; i<=highbound; i++ ) {
      char thisweight[50];
      //thisn does not have to be integer (in the case of pdf weights)
      float thisn = thishist->GetBinContent(i,pdfWeightIndex); //bin 0 has no special PDF weighting
      if (thisn==0) thisn=1; //avoid div by 0. if there are no events anyway then any value is ok
      sprintf(thisweight, "((SUSY_process==%d)*scanCrossSection%s/%f)",i,susyCrossSectionVariation_.Data(),thisn);
      susyprocessweight += thisweight;
      if (i!=10)    susyprocessweight += " + ";
    }
    susyprocessweight += ")";
    weightedcut+="*";
    weightedcut += susyprocessweight;
  }
  else if (type == kmSugraPlane && susySubProcess>=0) {
    char thisweight[50];
    sprintf(thisweight, "*((SUSY_process==%d)*scanCrossSection%s)",susySubProcess,susyCrossSectionVariation_.Data());
    weightedcut += thisweight;
  }
  else if (type == kSMSPoint) {
    //need to weight with the standard sigma / N ; lumi factor is done above
    assert(scanSMSngen!=0);
    int bin=  scanSMSngen->FindBin(m0_,m12_);
    double ngen=scanSMSngen->GetBinContent(bin);
    loadReferenceCrossSections();
    int xsbin=    referenceCrossSectionGluino->FindBin(m0_);
    double xs = referenceCrossSectionGluino->GetBinContent(xsbin);
    char xsweight[100];
    sprintf(xsweight,"*(%f/%f)",xs,ngen);
    weightedcut += xsweight;
  }
  else if (type == kmSugraPlane && susySubProcess<0) { //sanity check
    cout<<"You need to specify the susySubProcess!"<<endl; assert(0);
  }
  //for kSMSPlane we don't need to do anything!

  if (!quiet_)  cout<<weightedcut<<endl;
  return weightedcut;
}

//legacy interface!
//add an interface more like the old one, but with an addition for the data/MC lumi scaling
TString getCutString(bool isData, TString extraSelection="",TString extraWeight="",int pdfWeightIndex=0,TString pdfSet="CTEQ", bool isSusyScan=false, int susySubProcess=-1) {

  sampleType st=kMC;

  if (isData) st=kData;
  else if (isSusyScan) st = kmSugraPoint;
  else if (susySubProcess >=0) st=kmSugraPlane;
  return    getCutString(st, extraWeight, selection_,extraSelection, pdfWeightIndex,pdfSet,susySubProcess) ;
}


void addOverflowBin(TH1D* theHist) {

  if (theHist->GetBinErrorOption() == TH1::kPoisson) { //for new error treatment
    //    cout<<"Using kPoisson overflow addition: ";
    int nbins=    theHist->GetNbinsX();
    //    cout<<theHist->GetBinContent(nbins);
    theHist->SetBinContent(nbins,theHist->GetBinContent(nbins)+theHist->GetBinContent(nbins+1));
    //    cout<<" --> "<<theHist->GetBinContent(nbins)<<endl;
  }
  else { //legacy code
    //this code was written for when there was a customizable plot range (post-histo creation)
    //it could be made a lot simpler now
    
    int lastVisibleBin = theHist->GetNbinsX();
    //  cout<<theHist<<"  "<<lastVisibleBin<<"\t";
    
    //in case there is no custom range, the code should just add the overflow bin to the last bin of the histo
    double lastBinContent = theHist->GetBinContent(lastVisibleBin);
    double lastBinError = pow(theHist->GetBinError(lastVisibleBin),2); //square in prep for addition
    
    //  cout<<"Overflow addition: "<<lastBinContent<<" +/- "<<sqrt(lastBinError)<<" --> ";
    
    //now loop over the bins that aren't being shown at the moment (including the overflow bin)
    for (int ibin = lastVisibleBin+1; ibin <= 1 + theHist->GetNbinsX() ; ++ibin) {
      lastBinContent += theHist->GetBinContent( ibin);
      lastBinError += pow(theHist->GetBinError( ibin),2);
    }
    lastBinError = sqrt(lastBinError);
    
    theHist->SetBinContent(lastVisibleBin,lastBinContent);
    theHist->SetBinError(lastVisibleBin,lastBinError);
  }

}


void addOverflowBin(TH1F* theHist) {
  //copied from TH1D implementation
  if (theHist->GetBinErrorOption() == TH1::kPoisson) { //for new error treatment
    //    cout<<"Using kPoisson overflow addition: ";
    int nbins=    theHist->GetNbinsX();
    //    cout<<theHist->GetBinContent(nbins);
    theHist->SetBinContent(nbins,theHist->GetBinContent(nbins)+theHist->GetBinContent(nbins+1));
    //    cout<<" --> "<<theHist->GetBinContent(nbins)<<endl;
  }
  else { //legacy code
    //this code was written for when there was a customizable plot range (post-histo creation)
    //it could be made a lot simpler now
    
    int lastVisibleBin = theHist->GetNbinsX();
    //  cout<<theHist<<"  "<<lastVisibleBin<<"\t";
    
    //in case there is no custom range, the code should just add the overflow bin to the last bin of the histo
    double lastBinContent = theHist->GetBinContent(lastVisibleBin);
    double lastBinError = pow(theHist->GetBinError(lastVisibleBin),2); //square in prep for addition
    
    //  cout<<"Overflow addition: "<<lastBinContent<<" +/- "<<sqrt(lastBinError)<<" --> ";
    
    //now loop over the bins that aren't being shown at the moment (including the overflow bin)
    for (int ibin = lastVisibleBin+1; ibin <= 1 + theHist->GetNbinsX() ; ++ibin) {
      lastBinContent += theHist->GetBinContent( ibin);
      lastBinError += pow(theHist->GetBinError( ibin),2);
    }
    lastBinError = sqrt(lastBinError);
    
    theHist->SetBinContent(lastVisibleBin,lastBinContent);
    theHist->SetBinError(lastVisibleBin,lastBinError);
  }
}

void drawVerticalLine() {
  if (thecanvas==0) return;

  //this is a fine example of ROOT idiocy
  TVirtualPad* thePad = thecanvas->GetPad(0); //needs fixing for ratio plots
  double xmin,ymin,xmax,ymax;
  thePad->GetRangeAxis(xmin,ymin,xmax,ymax);
  //for academic interest, can get the same numbers using e.g. thePad->GetUymax()
  if (logy_) {
    ymax = pow(10, ymax);
    ymin = pow(10, ymin);
  }
  TLine theLine(verticalLinePosition_,ymin,verticalLinePosition_,ymax);
  theLine.SetLineColor(kBlue);
  theLine.SetLineWidth(3);

  theLine.DrawClone();

}

//add a sample to be plotted to the *end* of the list
void addSample(const TString & newsample) {

  //see if it is already there
  for (std::vector<TString>::iterator it = samples_.begin(); it!=samples_.end(); ++it) {
    if ( *it == newsample) {
      cout<<newsample<<" is already on the list!"<<endl;
      return;
    }
  }
  //if it isn't there, go ahead and add it

  if ( samplesAll_.find(newsample) != samplesAll_.end() ) {
    samples_.push_back(newsample);
  }
  else {
    cout<<"Could not find sample with name "<<newsample<<endl;
  }
}

void resetSampleScaleFactors() {
  for (std::set<TString>::iterator isample=samplesAll_.begin(); isample!=samplesAll_.end(); ++isample) {
    sampleScaleFactor_[*isample] = 1;
  }
}

void setSampleScaleFactor(const TString & sample, const float sf) {
  bool done=false;
  for (std::set<TString>::iterator isample=samplesAll_.begin(); isample!=samplesAll_.end(); ++isample) {
    if (*isample == sample) {sampleScaleFactor_[*isample] = sf; done=true;}
  }
  if (!done) cout<<"Failed to find the sample "<<sample<<endl;
}

void clearSamples() {
  samples_.clear();
}

void removeSample(const TString & sample) {

  for (std::vector<TString>::iterator it = samples_.begin(); it!=samples_.end(); ++it) {
    if ( *it == sample) {
      if (!quiet_) cout<<sample<<" removed from plotting list!"<<endl;
      samples_.erase(it);
      return;
    }
  }

  //if we get down here then we didn't find the sample
  cout<<sample<<" could not be found on the plotting list, so I did not remove it!"<<endl;
}

void setColorScheme(const TString & name) {
  if (name == "stack") {
    sampleColor_["LM13"] = kGray;
    sampleColor_["LM9"] =kGray;
    sampleColor_["mSUGRAtanb40"] =kGray;
    sampleColor_["T1bbbb"] =kGray;
    sampleColor_["T1tttt"] =kGray;
    sampleColor_["T2bb"] =kGray;
    sampleColor_["T2tt"] =kGray;
    sampleColor_["QCD"] = kYellow;
    sampleColor_["PythiaQCD"] = kYellow;
    sampleColor_["PythiaPUQCD"] = kYellow;
    sampleColor_["PythiaPUQCDFlat"] = kYellow;
    sampleColor_["TTbarJets"]=kRed+1;
    sampleColor_["TTbarSingleTopWJetsCombined"]=kRed+1;
    sampleColor_["ttbar"]=kRed+1;
    sampleColor_["TTbarJets-semiMu"]=kViolet;
    sampleColor_["TTbarJets-semiMuGood"]=kViolet;
    sampleColor_["TTbarJets-semiMuFailEta"]=kViolet;
    sampleColor_["TTbarJets-semiMuFailPt"]=kViolet;
    sampleColor_["TTbarJets-semiMuFailRecoIso"]=kViolet;
    sampleColor_["TTbarJets-semiMuFailOther"]=kViolet;
    sampleColor_["TTbarJets-semiEle"]=kViolet-9;
    sampleColor_["TTbarJets-semiEleGood"]=kViolet-9;
    sampleColor_["TTbarJets-semiEleFailEta"]=kViolet-9;
    sampleColor_["TTbarJets-semiEleFailPt"]=kViolet-9;
    sampleColor_["TTbarJets-semiEleFailRecoIso"]=kViolet-9;
    sampleColor_["TTbarJets-semiEleFailOther"]=kViolet-9;
    sampleColor_["TTbarJets-semiTauHad"]=kViolet-7;
    sampleColor_["TTbarJets-dilep"]=kMagenta-10;
    sampleColor_["TTbarJets-had"]=kRed-5;
    sampleColor_["TTbarJets-other"]=kPink-9;
    sampleColor_["SingleTop"] = kMagenta;
    sampleColor_["WJets"] = kGreen-3;
    sampleColor_["WJets-mu"]=kSpring+9;
    sampleColor_["WJets-muGood"]=kSpring+9;
    sampleColor_["WJets-muFailEta"]=kSpring+9;
    sampleColor_["WJets-muFailPt"]=kSpring+9;
    sampleColor_["WJets-muFailRecoIso"]=kSpring+9;
    sampleColor_["WJets-muFailOther"]=kSpring+9;
    sampleColor_["WJets-ele"]=kCyan;
    sampleColor_["WJets-eleGood"]=kCyan;
    sampleColor_["WJets-eleFailEta"]=kCyan;
    sampleColor_["WJets-eleFailPt"]=kCyan;
    sampleColor_["WJets-eleFailRecoIso"]=kCyan;
    sampleColor_["WJets-eleFailOther"]=kCyan;
    sampleColor_["WJets-tauHad"]=kGreen-10;
    sampleColor_["WJetsZ2"] = kGreen-3;
    sampleColor_["ZJets"] = kAzure-2;
    sampleColor_["Zinvisible"] = kOrange-3;
    sampleColor_["SingleTop-sChannel"] = kMagenta+1; //for special cases
    sampleColor_["SingleTop-tChannel"] = kMagenta+2; //for special cases
    sampleColor_["SingleTop-tWChannel"] = kMagenta+3; //for special cases
    sampleColor_["SingleTopBar-sChannel"] = kMagenta+4; //for special cases
    sampleColor_["SingleTopBar-tChannel"] = kMagenta+5; //for special cases
    sampleColor_["SingleTopBar-tWChannel"] = kMagenta+6; //for special cases
    //sampleColor_["SingleTop-sandtCombined"] = kMagenta+1; //for special cases
    //sampleColor_["SingleTop-tWCombined"] = kMagenta+3; //for special cases
    sampleColor_["SingleTop-sandtCombined-muGood"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-muFailEta"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-muFailPt"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-muFailRecoIso"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-muFailOther"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-eleGood"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-eleFailEta"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-eleFailPt"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-eleFailRecoIso"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-eleFailOther"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-tauHad"] = kMagenta+1;
    sampleColor_["SingleTop-tWCombined-semiMuGood"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiMuFailEta"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiMuFailPt"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiMuFailRecoIso"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiMuFailOther"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiEleGood"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiEleFailEta"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiEleFailPt"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiEleFailRecoIso"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiEleFailOther"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiTauHad"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-dilep"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-had"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-other"] = kMagenta+3;
    sampleColor_["TotalSM"] = kBlue+2;
    sampleColor_["Total"] = kGreen+3;
    sampleColor_["VV"] = kCyan+1;
    sampleColor_["HerwigQCDFlat"] = kYellow;
  }
  else if (name == "nostack" || name=="owen") {
    sampleColor_["LM13"] = kBlue+2;
    sampleColor_["LM9"] = kCyan+2;
    sampleColor_["mSUGRAtanb40"] =kCyan+2;
    sampleColor_["T1bbbb"] =kCyan+2;
    sampleColor_["T1tttt"] =kCyan+2;
    sampleColor_["T2bb"] =kCyan+2;
    sampleColor_["T2tt"] =kCyan+2;
    sampleColor_["QCD"] = 2;
    sampleColor_["PythiaQCD"] = 2;
    sampleColor_["PythiaPUQCD"] =2;
    sampleColor_["PythiaPUQCDFlat"] =2;
    sampleColor_["TTbarJets"]=4;
    sampleColor_["TTbarSingleTopWJetsCombined"]=4;
    sampleColor_["TTbarJets-semiMu"]=kGreen-3;
    sampleColor_["TTbarJets-semiMuGood"]=kGreen-3;
    sampleColor_["TTbarJets-semiMuFailEta"]=kGreen-3;
    sampleColor_["TTbarJets-semiMuFailPt"]=kGreen-3;
    sampleColor_["TTbarJets-semiMuFailRecoIso"]=kGreen-3;
    sampleColor_["TTbarJets-semiMuFailOther"]=kGreen-3;
    sampleColor_["TTbarJets-semiEle"]=kAzure-2;
    sampleColor_["TTbarJets-semiEleGood"]=kAzure-2;
    sampleColor_["TTbarJets-semiEleFailEta"]=kAzure-2;
    sampleColor_["TTbarJets-semiEleFailPt"]=kAzure-2;
    sampleColor_["TTbarJets-semiEleFailRecoIso"]=kAzure-2;
    sampleColor_["TTbarJets-semiEleFailOther"]=kAzure-2;
    sampleColor_["TTbarJets-semiTauHad"]=kBlack;
    sampleColor_["TTbarJets-dilep"]=kYellow;
    sampleColor_["TTbarJets-had"]=kRed-5;
    sampleColor_["TTbarJets-other"]=kPink-9;
    sampleColor_["SingleTop"] = kMagenta;
    sampleColor_["WJets"] = kOrange;
    sampleColor_["WJets-mu"]=kSpring+9;
    sampleColor_["WJets-muGood"]=kSpring+9;
    sampleColor_["WJets-muFailEta"]=kSpring+9;
    sampleColor_["WJets-muFailPt"]=kSpring+9;
    sampleColor_["WJets-muFailRecoIso"]=kSpring+9;
    sampleColor_["WJets-muFailOther"]=kSpring+9;
    sampleColor_["WJets-ele"]=kCyan;
    sampleColor_["WJets-eleGood"]=kCyan;
    sampleColor_["WJets-eleFailEta"]=kCyan;
    sampleColor_["WJets-eleFailPt"]=kCyan;
    sampleColor_["WJets-eleFailRecoIso"]=kCyan;
    sampleColor_["WJets-eleFailOther"]=kCyan;
    sampleColor_["WJets-tauHad"]=kGreen-10;
    sampleColor_["WJetsZ2"] = kOrange;
    sampleColor_["ZJets"] = 7;
    sampleColor_["Zinvisible"] = kOrange+7;
    sampleColor_["SingleTop-sChannel"] = kMagenta+1; //for special cases
    sampleColor_["SingleTop-tChannel"] = kMagenta+2; //for special cases
    sampleColor_["SingleTop-tWChannel"] = kMagenta+3; //for special cases
    sampleColor_["SingleTopBar-sChannel"] = kMagenta+4; //for special cases
    sampleColor_["SingleTopBar-tChannel"] = kMagenta+5; //for special cases
    sampleColor_["SingleTopBar-tWChannel"] = kMagenta+6; //for special cases
    //sampleColor_["SingleTop-sandtCombined"] = kMagenta+1; //for special cases
    //sampleColor_["SingleTop-tWCombined"] = kMagenta+3; //for special cases
    sampleColor_["SingleTop-sandtCombined-muGood"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-muFailEta"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-muFailPt"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-muFailRecoIso"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-muFailOther"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-eleGood"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-eleFailEta"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-eleFailPt"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-eleFailRecoIso"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-eleFailOther"] = kMagenta+1;
    sampleColor_["SingleTop-sandtCombined-tauHad"] = kMagenta+1;
    sampleColor_["SingleTop-tWCombined-semiMuGood"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiMuFailEta"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiMuFailPt"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiMuFailRecoIso"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiMuFailOther"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiEleGood"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiEleFailEta"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiEleFailPt"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiEleFailRecoIso"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiEleFailOther"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-semiTauHad"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-dilep"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-had"] = kMagenta+3;
    sampleColor_["SingleTop-tWCombined-other"] = kMagenta+3;
    sampleColor_["TotalSM"] =kGreen+1; //owen requested 3
    sampleColor_["Total"] = 6;
    sampleColor_["VV"] = kOrange-3;
    sampleColor_["HerwigQCDFlat"] = 2;
  }
  else {
    cout<<"Sorry, color scheme "<<name<<" is not known!"<<endl;
  }

}

void resetSamples(bool joinSingleTop=true) {

  samples_.clear();
  //this block controls what samples will enter your plot
  //order of this vector controls order of samples in stack

  //careful -- QCD must have 'QCD' in its name somewhere.
  //samples_.push_back("QCD"); //madgraph
  //samples_.push_back("PythiaQCD");
  samples_.push_back("PythiaPUQCD"); 
  //samples_.push_back("PythiaPUQCDFlat");

  if(splitTTbarForClosureTest_){
    samples_.push_back("TTbarJets-semiMuGood");
    samples_.push_back("TTbarJets-semiMuFailEta");
    samples_.push_back("TTbarJets-semiMuFailPt");
    samples_.push_back("TTbarJets-semiMuFailRecoIso");
    samples_.push_back("TTbarJets-semiMuFailOther");
    samples_.push_back("TTbarJets-semiEleGood");
    samples_.push_back("TTbarJets-semiEleFailEta");
    samples_.push_back("TTbarJets-semiEleFailPt");
    samples_.push_back("TTbarJets-semiEleFailRecoIso");
    samples_.push_back("TTbarJets-semiEleFailOther");
    samples_.push_back("TTbarJets-semiTauHad");
    samples_.push_back("TTbarJets-dilep");
    samples_.push_back("TTbarJets-had");
    samples_.push_back("TTbarJets-other");
  }
  else if(splitTTbar_){// "load" the decay modes separately
    samples_.push_back("TTbarJets-semiMu");
    samples_.push_back("TTbarJets-semiEle");
    samples_.push_back("TTbarJets-semiTauHad");
    samples_.push_back("TTbarJets-dilep");
    samples_.push_back("TTbarJets-had");
    samples_.push_back("TTbarJets-other");
  }
  else{
    samples_.push_back("TTbarJets");
  }
  if(splitWJetsForClosureTest_){// "load" the decay modes separately
    samples_.push_back("WJets-muGood");
    samples_.push_back("WJets-muFailEta");
    samples_.push_back("WJets-muFailPt");
    samples_.push_back("WJets-muFailRecoIso");
    samples_.push_back("WJets-muFailOther");
    samples_.push_back("WJets-eleGood");
    samples_.push_back("WJets-eleFailEta");
    samples_.push_back("WJets-eleFailPt");
    samples_.push_back("WJets-eleFailRecoIso");
    samples_.push_back("WJets-eleFailOther");
    samples_.push_back("WJets-tauHad");
  }
  else if(splitWJets_){// "load" the decay modes separately
    samples_.push_back("WJets-mu");
    samples_.push_back("WJets-ele");
    samples_.push_back("WJets-tauHad");
  }
  else{
    samples_.push_back("WJets");
  }
  //flip this bool to control whether SingleTop is loaded as one piece or 3
  //or two pieces for the tt+W+t enhanced closure test
  if(splitSingleTopForClosureTest_){
    //samples_.push_back("SingleTop-sandtCombined");
    samples_.push_back("SingleTop-sandtCombined-muGood");
    samples_.push_back("SingleTop-sandtCombined-muFailEta");
    samples_.push_back("SingleTop-sandtCombined-muFailPt");
    samples_.push_back("SingleTop-sandtCombined-muFailRecoIso");
    samples_.push_back("SingleTop-sandtCombined-muFailOther");
    samples_.push_back("SingleTop-sandtCombined-eleGood");
    samples_.push_back("SingleTop-sandtCombined-eleFailEta");
    samples_.push_back("SingleTop-sandtCombined-eleFailPt");
    samples_.push_back("SingleTop-sandtCombined-eleFailRecoIso");
    samples_.push_back("SingleTop-sandtCombined-eleFailOther");
    samples_.push_back("SingleTop-sandtCombined-tauHad");

    //samples_.push_back("SingleTop-tWCombined");
    samples_.push_back("SingleTop-tWCombined-semiMuGood");
    samples_.push_back("SingleTop-tWCombined-semiMuFailEta");
    samples_.push_back("SingleTop-tWCombined-semiMuFailPt");
    samples_.push_back("SingleTop-tWCombined-semiMuFailRecoIso");
    samples_.push_back("SingleTop-tWCombined-semiMuFailOther");
    samples_.push_back("SingleTop-tWCombined-semiEleGood");
    samples_.push_back("SingleTop-tWCombined-semiEleFailEta");
    samples_.push_back("SingleTop-tWCombined-semiEleFailPt");
    samples_.push_back("SingleTop-tWCombined-semiEleFailRecoIso");
    samples_.push_back("SingleTop-tWCombined-semiEleFailOther");
    samples_.push_back("SingleTop-tWCombined-semiTauHad");
    samples_.push_back("SingleTop-tWCombined-dilep");
    samples_.push_back("SingleTop-tWCombined-had");
    samples_.push_back("SingleTop-tWCombined-other");
  }
  else if (joinSingleTop) samples_.push_back("SingleTop");
  else {
    samples_.push_back("SingleTop-sChannel");
    samples_.push_back("SingleTop-tChannel");
    samples_.push_back("SingleTop-tWChannel");
    samples_.push_back("SingleTopBar-sChannel");
    samples_.push_back("SingleTopBar-tChannel");
    samples_.push_back("SingleTopBar-tWChannel");
  }
  samples_.push_back("ZJets");
  samples_.push_back("VV");
  samples_.push_back("Zinvisible");
  //samples_.push_back("HerwigQCDFlat");
  samples_.push_back("LM9");

}

TString getSampleLabel(const TString & sample) {
  TString label="Not Found";
  if (  isSampleScan(sample)) {
    char ss[50];
    sprintf(ss, "#splitline{m_{gluino} = %d GeV}{m_{LSP} = %d GeV}",m0_,m12_);
    label=ss;
  }
  else {
    label = sampleLabel_[sample]; //should replace this with find() instead
  }
  return label;

}


void loadSamples(bool joinSingleTop=true, TString signalEffMode="") {
  if (loaded_) return;
  loaded_=true;

  resetPadDimensions();
  resetLegendPosition();

  resetSamples(joinSingleTop);
  //samplesAll_ should have *every* available sample
  //also note that there's no harm in failing to load one of these samples, 
  //as long as you don't actually try to draw it

  //samplesAll_.insert("QCD");
  //samplesAll_.insert("PythiaQCD");
  samplesAll_.insert("PythiaPUQCD");
  //samplesAll_.insert("PythiaPUQCDFlat");
  samplesAll_.insert("TTbarJets");
  samplesAll_.insert("TTbarSingleTopWJetsCombined");
  samplesAll_.insert("TTbarJets-semiMu");
  samplesAll_.insert("TTbarJets-semiMuGood");
  samplesAll_.insert("TTbarJets-semiMuFailEta");
  samplesAll_.insert("TTbarJets-semiMuFailPt");
  samplesAll_.insert("TTbarJets-semiMuFailRecoIso");
  samplesAll_.insert("TTbarJets-semiMuFailOther");
  samplesAll_.insert("TTbarJets-semiEle");
  samplesAll_.insert("TTbarJets-semiEleGood");
  samplesAll_.insert("TTbarJets-semiEleFailEta");
  samplesAll_.insert("TTbarJets-semiEleFailPt");
  samplesAll_.insert("TTbarJets-semiEleFailRecoIso");
  samplesAll_.insert("TTbarJets-semiEleFailOther");
  samplesAll_.insert("TTbarJets-semiTauHad");
  samplesAll_.insert("TTbarJets-dilep");
  samplesAll_.insert("TTbarJets-had");
  samplesAll_.insert("TTbarJets-other");
  samplesAll_.insert("WJets");
  samplesAll_.insert("WJets-muGood");
  samplesAll_.insert("WJets-muFailEta");
  samplesAll_.insert("WJets-muFailPt");
  samplesAll_.insert("WJets-muFailRecoIso");
  samplesAll_.insert("WJets-muFailOther");
  samplesAll_.insert("WJets-eleGood");
  samplesAll_.insert("WJets-eleFailEta");
  samplesAll_.insert("WJets-eleFailPt");
  samplesAll_.insert("WJets-eleFailRecoIso");
  samplesAll_.insert("WJets-eleFailOther");
  samplesAll_.insert("WJets-mu");
  samplesAll_.insert("WJets-ele");
  samplesAll_.insert("WJets-tauHad");
  samplesAll_.insert("ZJets");
  samplesAll_.insert("Zinvisible");
  samplesAll_.insert("SingleTop");
  samplesAll_.insert("SingleTop-sChannel");
  samplesAll_.insert("SingleTop-tChannel");
  samplesAll_.insert("SingleTop-tWChannel");
  samplesAll_.insert("SingleTopBar-sChannel");
  samplesAll_.insert("SingleTopBar-tChannel");
  samplesAll_.insert("SingleTopBar-tWChannel");
  //samplesAll_.insert("SingleTop-sandtCombined");
  //samplesAll_.insert("SingleTop-tWCombined");
  samplesAll_.insert("SingleTop-sandtCombined-muGood");
  samplesAll_.insert("SingleTop-sandtCombined-muFailEta");
  samplesAll_.insert("SingleTop-sandtCombined-muFailPt");
  samplesAll_.insert("SingleTop-sandtCombined-muFailRecoIso");
  samplesAll_.insert("SingleTop-sandtCombined-muFailOther");
  samplesAll_.insert("SingleTop-sandtCombined-eleGood");
  samplesAll_.insert("SingleTop-sandtCombined-eleFailEta");
  samplesAll_.insert("SingleTop-sandtCombined-eleFailPt");
  samplesAll_.insert("SingleTop-sandtCombined-eleFailRecoIso");
  samplesAll_.insert("SingleTop-sandtCombined-eleFailOther");
  samplesAll_.insert("SingleTop-sandtCombined-tauHad");
  samplesAll_.insert("SingleTop-tWCombined-semiMuGood");
  samplesAll_.insert("SingleTop-tWCombined-semiMuFailEta");
  samplesAll_.insert("SingleTop-tWCombined-semiMuFailPt");
  samplesAll_.insert("SingleTop-tWCombined-semiMuFailRecoIso");
  samplesAll_.insert("SingleTop-tWCombined-semiMuFailOther");
  samplesAll_.insert("SingleTop-tWCombined-semiEleGood");
  samplesAll_.insert("SingleTop-tWCombined-semiEleFailEta");
  samplesAll_.insert("SingleTop-tWCombined-semiEleFailPt");
  samplesAll_.insert("SingleTop-tWCombined-semiEleFailRecoIso");
  samplesAll_.insert("SingleTop-tWCombined-semiEleFailOther");
  samplesAll_.insert("SingleTop-tWCombined-semiTauHad");
  samplesAll_.insert("SingleTop-tWCombined-dilep");
  samplesAll_.insert("SingleTop-tWCombined-had");
  samplesAll_.insert("SingleTop-tWCombined-other");
  samplesAll_.insert("HerwigQCDFlat");
  //  samplesAll_.insert("WJetsZ2");
  //  samplesAll_.insert("ZJetsZ2");
  samplesAll_.insert("VV");


  samplesAll_.insert("LM13");
  samplesAll_.insert("LM9");

  samplesAll_.insert("mSUGRAtanb40");
  samplesAll_.insert("T1bbbb");
  samplesAll_.insert("T1tttt");
  samplesAll_.insert("T2bb");
  samplesAll_.insert("T2tt");

 
  
  //FOR PLOTS
  ////////////
  if (signalEffMode=="") {
    configDescriptions_.setDefault("CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0");
    configDescriptions_.setCorrected("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0");
    //configDescriptions_.setDefault("CSVM_PFMETTypeI_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0");
    //configDescriptions_.setCorrected("CSVM_PFMETTypeI_PF2PATjets_JES0_JERbias_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0");

  }
  else if (signalEffMode=="scan" || signalEffMode=="complete") {  //Only for signal systematics
      
    configDescriptions_.setDefault("CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0");
    configDescriptions_.setCorrected("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0");

    //convention: put 'down' variation first and 'up' variation second

    //JES //LM9 and scans
    configDescriptions_.addVariation("CSVM_PF2PATjets_JESdown_JERbias_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0",
				     "CSVM_PF2PATjets_JESup_JERbias_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0");
    //JER //LM9 only
    if (signalEffMode=="complete") {
      configDescriptions_.addVariation("CSVM_PF2PATjets_JES0_JERdown_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0",
				       "CSVM_PF2PATjets_JES0_JERup_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0");
    }
    //unclustered MET //LM9 and scans
    configDescriptions_.addVariation("CSVM_PF2PATjets_JES0_JERbias_PFMET_METuncDown_PUunc0_BTagEff04_HLTEff0", 
				     "CSVM_PF2PATjets_JES0_JERbias_PFMET_METuncUp_PUunc0_BTagEff04_HLTEff0");

    if (signalEffMode=="complete") {
      //PU //LM9 only
      configDescriptions_.addVariation("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUuncDown_BTagEff04_HLTEff0",
				       "CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUuncUp_BTagEff04_HLTEff0");
    }
    
    //UPDATE -- with 'm' series reducedTrees, don't do this one anymore
    //btag eff //LM9 and scans
    //    configDescriptions_.addVariation("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEffdown4_HLTEff0",
    //				     "CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEffup4_HLTEff0");
    
    //HLT eff //never use this one
    //    configDescriptions_.addVariation("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff0_HLTEffdown",
    //"CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff0_HLTEffup");
  }
  else assert(0);
  ///////////////
  //////////////
 

  currentConfig_=configDescriptions_.getDefault();

  //these blocks are just a "dictionary"
  //no need to ever comment these out
  setColorScheme("stack");

  sampleLabel_["mSUGRAtanb40"] = "tan #beta = 40";
  sampleLabel_["T1bbbb"] = "T1bbbb";
  sampleLabel_["T1tttt"] = "T1tttt";
  sampleLabel_["T2bb"] = "T2bb";
  sampleLabel_["T2tt"] = "T2tt";
  sampleLabel_["LM13"] = "LM13";
  sampleLabel_["LM9"] = "LM9";
  sampleLabel_["QCD"] = "QCD";
  sampleLabel_["PythiaQCD"] = "QCD (no PU)";
  sampleLabel_["PythiaPUQCDFlat"] = "QCD"; 
  sampleLabel_["PythiaPUQCD"] = "QCD";
  sampleLabel_["TTbarJets"]="t#bar{t}";
  sampleLabel_["TTbarSingleTopWJetsCombined"]="t#bar{t}+Single-Top+W#rightarrowl#nu";
  sampleLabel_["ttbar"]="t#bar{t}+W+t"; //for DD ttwt
  sampleLabel_["TTbarJets-semiMu"]="t#bar{t}:semi-#mu";
  sampleLabel_["TTbarJets-semiMuGood"]="t#bar{t}:semi-#mu - good";
  sampleLabel_["TTbarJets-semiMuFailEta"]="t#bar{t}:semi-#mu - faileta";
  sampleLabel_["TTbarJets-semiMuFailPt"]="t#bar{t}:semi-#mu - failpt";
  sampleLabel_["TTbarJets-semiMuFailRecoIso"]="t#bar{t}:semi-#mu - failrecoiso";
  sampleLabel_["TTbarJets-semiMuFailOther"]="t#bar{t}:semi-#mu - failother";
  sampleLabel_["TTbarJets-semiEle"]="t#bar{t}:semi-e";
  sampleLabel_["TTbarJets-semiEleGood"]="t#bar{t}:semi-e - good";
  sampleLabel_["TTbarJets-semiEleFailEta"]="t#bar{t}:semi-e -faileta";
  sampleLabel_["TTbarJets-semiEleFailPt"]="t#bar{t}:semi-e - failpt";
  sampleLabel_["TTbarJets-semiEleFailRecoIso"]="t#bar{t}:semi-e - failrecoiso";
  sampleLabel_["TTbarJets-semiEleFailOther"]="t#bar{t}:semi-e - failother";
  sampleLabel_["TTbarJets-semiTauHad"]="t#bar{t}:semi-#tau(#rightarrow had)";
  sampleLabel_["TTbarJets-dilep"]="t#bar{t}:ee,#mu#mu,e#mu";
  sampleLabel_["TTbarJets-had"]="t#bar{t}:fully hadronic";
  sampleLabel_["TTbarJets-other"]="t#bar{t}:other";
  sampleLabel_["SingleTop"] = "Single-Top";
  sampleLabel_["WJets"] = "W#rightarrowl#nu";
  sampleLabel_["WJets-mu"] = "W#rightarrow#mu#nu";
  sampleLabel_["WJets-muGood"] = "W#rightarrow#mu#nu - good";
  sampleLabel_["WJets-muFailEta"] = "W#rightarrow#mu#nu - faileta";
  sampleLabel_["WJets-muFailPt"] = "W#rightarrow#mu#nu - failpt";
  sampleLabel_["WJets-muFailRecoIso"] = "W#rightarrow#mu#nu - failrecoiso";
  sampleLabel_["WJets-muFailOther"] = "W#rightarrow#mu#nu - failother";
  sampleLabel_["WJets-ele"] = "W#rightarrow e#nu";
  sampleLabel_["WJets-eleGood"] = "W#rightarrow e#nu - good";
  sampleLabel_["WJets-eleFailEta"] = "W#rightarrow e#nu - faileta";
  sampleLabel_["WJets-eleFailPt"] = "W#rightarrow e#nu - failpt";
  sampleLabel_["WJets-eleFailRecoIso"] = "W#rightarrow e#nu - failrecoiso";
  sampleLabel_["WJets-eleFailOther"] = "W#rightarrow e#nu - failother";
  sampleLabel_["WJets-tauHad"] = "W#rightarrow#tau(#rightarrow had)#nu";
  sampleLabel_["WJetsZ2"] = "W#rightarrowl#nu (Z2)";
  sampleLabel_["ZJets"] = "Z/#gamma*#rightarrowl^{+}l^{-}";
  sampleLabel_["Zinvisible"] = "Z#rightarrow#nu#nu";
  sampleLabel_["SingleTop-sChannel"] = "Single-Top (s)";
  sampleLabel_["SingleTop-tChannel"] = "Single-Top (t)";
  sampleLabel_["SingleTop-tWChannel"] = "Single-Top (tW)";
  sampleLabel_["SingleTopBar-sChannel"] = "Single-TopBar (s)";
  sampleLabel_["SingleTopBar-tChannel"] = "Single-TopBar (t)";
  sampleLabel_["SingleTopBar-tWChannel"] = "Single-TopBar (tW)";
  //sampleLabel_["SingleTop-sandtCombined"] = "Single-Top (s and t)";
  //sampleLabel_["SingleTop-tWCombined"] = "Single-Top (tW)";
  sampleLabel_["SingleTop-sandtCombined-muGood"] = "Single-Top (s and t)-mu-good";
  sampleLabel_["SingleTop-sandtCombined-muFailEta"] = "Single-Top (s and t)-mu-faileta";
  sampleLabel_["SingleTop-sandtCombined-muFailPt"] = "Single-Top (s and t)-mu-failpt";
  sampleLabel_["SingleTop-sandtCombined-muFailRecoIso"] = "Single-Top (s and t)-mu-failrecoiso";
  sampleLabel_["SingleTop-sandtCombined-muFailOther"] = "Single-Top (s and t)-mu-failother";
  sampleLabel_["SingleTop-sandtCombined-eleGood"] = "Single-Top (s and t)-e-good";
  sampleLabel_["SingleTop-sandtCombined-eleFailEta"] = "Single-Top (s and t)-e-faileta";
  sampleLabel_["SingleTop-sandtCombined-eleFailPt"] = "Single-Top (s and t)-e-failpt";
  sampleLabel_["SingleTop-sandtCombined-eleFailRecoIso"] = "Single-Top (s and t)-e-failrecoiso";
  sampleLabel_["SingleTop-sandtCombined-eleFailOther"] = "Single-Top (s and t)-e-failother";
  sampleLabel_["SingleTop-sandtCombined-tauHad"] = "Single-Top (s and t)-tauhad";
  sampleLabel_["SingleTop-tWCombined-semiMuGood"] = "Single-Top (tW)-mu-good";
  sampleLabel_["SingleTop-tWCombined-semiMuFailEta"] = "Single-Top (tW)-mu-faileta";
  sampleLabel_["SingleTop-tWCombined-semiMuFailPt"] = "Single-Top (tW)-mu-failpt";
  sampleLabel_["SingleTop-tWCombined-semiMuFailRecoIso"] = "Single-Top (tW)-mu-failrecoiso";
  sampleLabel_["SingleTop-tWCombined-semiMuFailOther"] = "Single-Top (tW)-mu-failother";
  sampleLabel_["SingleTop-tWCombined-semiEleGood"] = "Single-Top (tW)-e-good";
  sampleLabel_["SingleTop-tWCombined-semiEleFailEta"] = "Single-Top (tW)-e-faileta";
  sampleLabel_["SingleTop-tWCombined-semiEleFailPt"] = "Single-Top (tW)-e-failpt";
  sampleLabel_["SingleTop-tWCombined-semiEleFailRecoIso"] = "Single-Top (tW)-e-failrecoiso";
  sampleLabel_["SingleTop-tWCombined-semiEleFailOther"] = "Single-Top (tW)-e-failother";
  sampleLabel_["SingleTop-tWCombined-semiTauHad"] = "Single-Top (tW)-tauhad";
  sampleLabel_["SingleTop-tWCombined-dilep"] = "Single-Top (tW)-dilep";
  sampleLabel_["SingleTop-tWCombined-had"] = "Single-Top (tW)-had";
  sampleLabel_["SingleTop-tWCombined-other"] = "Single-Top (tW)-other";
  sampleLabel_["VV"] = "Diboson";
  sampleLabel_["HerwigQCDFlat"] = "Herwig QCD";
  sampleLabel_["TotalSM"] = "SM";
  sampleLabel_["Total"] = "SM + LM13"; //again, this is a hack

  sampleMarkerStyle_["mSUGRAtanb40"] = kFullStar;
  sampleMarkerStyle_["T1bbbb"] = kFullStar;
  sampleMarkerStyle_["T1tttt"] = kFullStar;
  sampleMarkerStyle_["T2bb"] = kFullStar;
  sampleMarkerStyle_["T2tt"] = kFullStar;
  sampleMarkerStyle_["LM13"] = kFullStar;
  sampleMarkerStyle_["LM9"] = kFullStar;
  sampleMarkerStyle_["QCD"] = kFullCircle;
  sampleMarkerStyle_["PythiaQCD"] = kOpenCircle;
  sampleMarkerStyle_["PythiaPUQCDFlat"] = kOpenCircle;  
  sampleMarkerStyle_["PythiaPUQCD"] = kOpenCircle;
  sampleMarkerStyle_["TTbarJets"]= kFullSquare;
  sampleMarkerStyle_["TTbarSingleTopWJetsCombined"]= kFullSquare;
  sampleMarkerStyle_["ttbar"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets-semiMu"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets-semiMuGood"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets-semiMuFailEta"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets-semiMuFailPt"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets-semiMuFailRecoIso"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets-semiMuFailOther"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets-semiEle"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets-semiEleGood"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets-semiEleFailEta"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets-semiEleFailPt"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets-semiEleFailRecoIso"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets-semiEleFailOther"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets-semiTauHad"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets-dilep"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets-had"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets-other"]= kFullSquare;
  sampleMarkerStyle_["SingleTop"] = kOpenSquare;
  sampleMarkerStyle_["WJets"] = kMultiply;
  sampleMarkerStyle_["WJets-mu"] = kMultiply;
  sampleMarkerStyle_["WJets-muGood"] = kMultiply;
  sampleMarkerStyle_["WJets-muFailEta"] = kMultiply;
  sampleMarkerStyle_["WJets-muFailPt"] = kMultiply;
  sampleMarkerStyle_["WJets-muFailRecoIso"] = kMultiply;
  sampleMarkerStyle_["WJets-muFailOther"] = kMultiply;
  sampleMarkerStyle_["WJets-ele"] = kMultiply;
  sampleMarkerStyle_["WJets-eleGood"] = kMultiply;
  sampleMarkerStyle_["WJets-eleFailEta"] = kMultiply;
  sampleMarkerStyle_["WJets-eleFailPt"] = kMultiply;
  sampleMarkerStyle_["WJets-eleFailRecoIso"] = kMultiply;
  sampleMarkerStyle_["WJets-eleFailOther"] = kMultiply;
  sampleMarkerStyle_["WJets-tauHad"] = kMultiply;
  sampleMarkerStyle_["WJetsZ2"] = kMultiply;
  sampleMarkerStyle_["ZJets"] = kFullTriangleUp;
  sampleMarkerStyle_["Zinvisible"] = kFullTriangleDown;
  sampleMarkerStyle_["HerwigQCDFlat"] = kFullCircle;
  sampleMarkerStyle_["VV"] = kOpenCross;
  sampleMarkerStyle_["SingleTop-sChannel"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tChannel"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tWChannel"] = kOpenSquare;
  sampleMarkerStyle_["SingleTopBar-sChannel"] = kOpenSquare;
  sampleMarkerStyle_["SingleTopBar-tChannel"] = kOpenSquare;
  sampleMarkerStyle_["SingleTopBar-tWChannel"] = kOpenSquare;
  //sampleMarkerStyle_["SingleTop-sandtCombined"] = kOpenSquare;
  //sampleMarkerStyle_["SingleTop-tWCombined"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-sandtCombined-muGood"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-sandtCombined-muFailEta"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-sandtCombined-muFailPt"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-sandtCombined-muFailRecoIso"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-sandtCombined-muFailOther"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-sandtCombined-eleGood"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-sandtCombined-eleFailEta"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-sandtCombined-eleFailPt"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-sandtCombined-eleFailRecoIso"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-sandtCombined-eleFailOther"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-sandtCombined-tauHad"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tWCombined-semiMuGood"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tWCombined-semiMuFailEta"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tWCombined-semiMuFailPt"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tWCombined-semiMuFailRecoIso"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tWCombined-semiMuFailOther"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tWCombined-semiEleGood"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tWCombined-semiEleFailEta"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tWCombined-semiEleFailPt"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tWCombined-semiEleFailRecoIso"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tWCombined-semiEleFailOther"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tWCombined-semiTauHad"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tWCombined-dilep"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tWCombined-had"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tWCombined-other"] = kOpenSquare;
  sampleMarkerStyle_["TotalSM"] = kOpenCross; //FIXME?
  sampleMarkerStyle_["Total"] = kDot; //FIXME?

  sampleOwenName_["mSUGRAtanb40"] = "msugra40";
  sampleOwenName_["T1bbbb"] = "t1bbbb";
  sampleOwenName_["T1tttt"] = "t1tttt";
  sampleOwenName_["T2bb"] = "t2bb";
  sampleOwenName_["T2tt"] = "t2tt";
  sampleOwenName_["LM13"] = "lm13";
  sampleOwenName_["LM9"] = "lm9";
  sampleOwenName_["QCD"] = "qcd";
  sampleOwenName_["PythiaQCD"] = "qcd";
  sampleOwenName_["PythiaPUQCDFlat"] = "qcd"; 
  sampleOwenName_["PythiaPUQCD"] = "qcd";
  sampleOwenName_["TTbarJets"]="ttbar";
  sampleOwenName_["TTbarSingleTopWJetsCombined"]="ttbartw";
  sampleOwenName_["ttbar"]="ttbar";
  sampleOwenName_["SingleTop"] = "singletop";
  sampleOwenName_["WJets"] = "wjets";
  sampleOwenName_["WJetsZ2"] = "wjets";
  sampleOwenName_["ZJets"] = "zjets";
  sampleOwenName_["HerwigQCDFlat"] = "herwigqcdflat";
  sampleOwenName_["Zinvisible"] = "zinvis";
  sampleOwenName_["VV"] = "vv";
  sampleOwenName_["SingleTop-sChannel"] = "singletops";
  sampleOwenName_["SingleTop-tChannel"] = "singletopt";
  sampleOwenName_["SingleTop-tWChannel"] = "singletoptw";
  sampleOwenName_["SingleTopBar-sChannel"] = "singletopbars";
  sampleOwenName_["SingleTopBar-tChannel"] = "singletopbart";
  sampleOwenName_["SingleTopBar-tWChannel"] = "singletopbartw";
  //sampleOwenName_["SingleTop-sandtCombined"] = "singletopsandt";
  //sampleOwenName_["SingleTop-tWCombined"] = "singletoptw";
  sampleOwenName_["SingleTop-sandtCombined-muGood"] = "singletopsandt";
  sampleOwenName_["SingleTop-sandtCombined-muFailEta"] = "singletopsandt";
  sampleOwenName_["SingleTop-sandtCombined-muFailPt"] = "singletopsandt";
  sampleOwenName_["SingleTop-sandtCombined-muFailRecoIso"] = "singletopsandt";
  sampleOwenName_["SingleTop-sandtCombined-muFailOther"] = "singletopsandt";
  sampleOwenName_["SingleTop-sandtCombined-eleGood"] = "singletopsandt";
  sampleOwenName_["SingleTop-sandtCombined-eleFailEta"] = "singletopsandt";
  sampleOwenName_["SingleTop-sandtCombined-eleFailPt"] = "singletopsandt";
  sampleOwenName_["SingleTop-sandtCombined-eleFailRecoIso"] = "singletopsandt";
  sampleOwenName_["SingleTop-sandtCombined-eleFailOther"] = "singletopsandt";
  sampleOwenName_["SingleTop-sandtCombined-tauHad"] = "singletopsandt";
  sampleOwenName_["SingleTop-tWCombined-semiMuGood"] = "singletoptw";
  sampleOwenName_["SingleTop-tWCombined-semiMuFailEta"] = "singletoptw";
  sampleOwenName_["SingleTop-tWCombined-semiMuFailPt"] = "singletoptw";
  sampleOwenName_["SingleTop-tWCombined-semiMuFailRecoIso"] = "singletoptw";
  sampleOwenName_["SingleTop-tWCombined-semiMuFailOther"] = "singletoptw";
  sampleOwenName_["SingleTop-tWCombined-semiEleGood"] = "singletoptw";
  sampleOwenName_["SingleTop-tWCombined-semiEleFailEta"] = "singletoptw";
  sampleOwenName_["SingleTop-tWCombined-semiEleFailPt"] = "singletoptw";
  sampleOwenName_["SingleTop-tWCombined-semiEleFailRecoIso"] = "singletoptw";
  sampleOwenName_["SingleTop-tWCombined-semiEleFailOther"] = "singletoptw";
  sampleOwenName_["SingleTop-tWCombined-semiTauHad"] = "singletoptw";
  sampleOwenName_["SingleTop-tWCombined-dilep"] = "singletoptw";
  sampleOwenName_["SingleTop-tWCombined-had"] = "singletoptw";
  sampleOwenName_["SingleTop-tWCombined-other"] = "singletoptw";
  sampleOwenName_["TotalSM"] = "totalsm";
  sampleOwenName_["Total"] = "total";  

  //set all scale factors to 1 to start with.
  resetSampleScaleFactors();

  //  for (std::vector<TString>::iterator iconfig=configDescriptions_.begin(); iconfig!=configDescriptions_.end(); ++iconfig) {
  for (unsigned int iconfig=0; iconfig<configDescriptions_.size(); ++iconfig) {
    for (std::set<TString>::iterator isample=samplesAll_.begin(); isample!=samplesAll_.end(); ++isample) {
      TString fname="reducedTree.";
      TString thisconfig = configDescriptions_.at(iconfig);
      fname += thisconfig;
      fname+=".";

      if((splitTTbarForClosureTest_ || splitTTbar_) && (*isample).Contains("TTbarJets") )
	fname += "TTbarJets";
      else if((splitWJetsForClosureTest_ || splitWJets_) && (*isample).Contains("WJets") )
	fname += "WJets";
      else if((splitSingleTopForClosureTest_) && (*isample).Contains("SingleTop-sandtCombined"))
	fname += "SingleTop-sandtCombined";
      else if((splitSingleTopForClosureTest_) && (*isample).Contains("SingleTop-tWCombined"))
	fname += "SingleTop-tWCombined";
      else
	fname += *isample;

      fname+=".root";

      fname.Prepend(inputPath);

      files_[thisconfig][*isample] = new TFile(fname);
      if (files_[thisconfig][*isample]->IsZombie() ) {cout<<"file error with "<<*isample<<endl; files_[thisconfig][*isample]=0;}
      else { if (!quiet_)    cout<<"Added sample: "<<thisconfig<<"\t"<<*isample<<endl;}
    }
  }

  //load data file too
  TString dname="reducedTree.";
  dname+=currentConfig_;
  //dname+=".data.root";
  //dname+=".ht_run2011a_SUM_promptrecov4only_uptojun24.root";
  //dname+=".data_promptrecoThroughJul1.root";
  dname+=".ht*.root";
  //dname+=".singlemu*.root";
  //dname+="*.root";
  //dname+=".ht_run2011a_SUM_promptrecov4only_uptojul1.root";
  dname.Prepend(dataInputPath);
  dname.ReplaceAll("JERbias","JER0"); //JERbias not relevant for data
  if ( dodata_) {
    dtree = new TChain("reducedTree");
    dtree->Add(dname);
  }

}

//ug...not liking this
sampleType getSampleType(const TString & sample , const TString & planeOrPoint="") {

  assert(planeOrPoint=="" || planeOrPoint=="point" || planeOrPoint=="plane");

  if (sample=="data") return kData;
  else if (sample.Contains("mSUGRA") && planeOrPoint=="point") return kmSugraPoint;
  else if (sample.Contains("mSUGRA") && planeOrPoint=="plane") return kmSugraPlane;
  else if (sample.Contains("T1bbbb") && planeOrPoint=="point") return kSMSPoint;
  else if (sample.Contains("T1bbbb") && planeOrPoint=="plane") return kSMSPlane;
  else if (sample.Contains("T1tttt") && planeOrPoint=="point") return kSMSPoint;
  else if (sample.Contains("T1tttt") && planeOrPoint=="plane") return kSMSPlane;
  else if (sample.Contains("T2bb") && planeOrPoint=="point") return kSMSPoint;
  else if (sample.Contains("T2bb") && planeOrPoint=="plane") return kSMSPlane;
  else if (sample.Contains("T2tt") && planeOrPoint=="point") return kSMSPoint;
  else if (sample.Contains("T2tt") && planeOrPoint=="plane") return kSMSPlane;

  return kMC;

}

//beginning of something i've wanted for a long time. functions that return some standard cuts
TCut getLeptonVetoCut() {
  TCut a = "nElectrons==0 && nMuons==0"; //default cut
  
  //require lost lepton due to pT (ttbar only)
  //TCut a="nElectrons==0 && nMuons==0 && (decayType==201102 ||decayType==101102 ||decayType==101302 ||decayType==201302)";
  return a;
}
TCut getSingleLeptonCut() {
  TCut sl = "(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)";
  return sl;
}

TCut getSingleMuonCut() {
  TCut sl = "(nElectrons==0 && nMuons==1 && MT_Wlep>=0&&MT_Wlep<100)";
  return sl;
}

TCut getSingleElectronCut() {
  TCut sl = "(nElectrons==1 && nMuons==0 && MT_Wlep>=0&&MT_Wlep<100)";
  return sl;
}


//if something is passed to varbins, then low and high will be ignored
float drawSimple(const TString var, const int nbins, const double low, const double high, const TString filename, 
		 const TString histname , const TString samplename, const float* varbins=0) {

  loadSamples();
  gROOT->SetStyle("CMS");
//I would rather implement this functionality via drawPlots(), but I think it will be simpler
//to just write something simple

//no presentation, just fill the histogram and save
  TTree* tree=0;
  if (samplename=="data") {
    tree = dtree;
  }
  else {
    for (unsigned int isample=0; isample<samples_.size(); isample++) {
      if ( samples_[isample] == samplename) {
	if (!quiet_) cout <<samples_[isample]<<endl;
	tree = (TTree*) files_[currentConfig_][samples_[isample]]->Get("reducedTree");
      }
    }
  }
  if (tree==0) {cout<<"Something went wrong finding your sample (" << samplename <<")!"<<endl; return 0;}
  gROOT->cd();

  
  TH1D* hh=0;
  if (varbins==0) {
    hh = new TH1D(histname,histname,nbins,low,high);
  }
  else {
    hh = new TH1D(histname,histname,nbins,varbins);
  }
  hh->Sumw2();
   
  TString optfh= useFlavorHistoryWeights_ && samplename.Contains("WJets") ? "flavorHistoryWeight" : "";
  if (samplename=="data") tree->Project(histname,var,getCutString(true).Data());
  else tree->Project(histname,var,getCutString( getSampleType(samplename,"point"),optfh,selection_,"",0,"").Data());
  float theIntegral = hh->Integral(0,nbins+1);

  if (addOverflow_)  addOverflowBin( hh ); //manipulates the TH1D

  //at this point i've got a histogram. what more could i want?
  hinteractive = (TH1D*)hh->Clone("hinteractive");
  if(savePlots_){
    TFile fout(filename,"UPDATE");
    hh->Write();
    fout.Close();
  }
  delete hh; //deleting ROOT objects can be dangerous...but i've tried carefully to avoid a double deletion here by doing gROOT->cd() before the creation of hh
  return theIntegral;
}

float drawSimple(const TString var, const int nbins, const float* varbins, const TString filename, 
		 const TString histname , const TString samplename) {
  return drawSimple(var, nbins, 0, 1, filename, histname, samplename, varbins);
}


void draw2d(const TString var, const int nbins, const float low, const float high, 
	    const TString vary, const int nbinsy, const float lowy, const float highy,
	    const TString xtitle, TString ytitle, TString filename="",
	    const float* varbins=0,const float* varbinsy=0) {
  //in case varbins is used, the low and high will not be used

  loadSamples();
  if (filename=="") filename=var;
  gROOT->SetStyle("CMS");

  renewCanvas();
  if (h2d!=0) delete h2d;
  const TString hname="h2d";
  if (varbins==0 &&varbinsy==0)  h2d = new TH2D(hname,"",nbins,low,high,nbinsy,lowy,highy);
  else  h2d = new TH2D(hname,"",nbins,varbins,nbinsy,varbinsy);
  h2d->Sumw2();

  TH2D*  h2d_temp =0;
  if (varbins==0 &&varbinsy==0) h2d_temp= new TH2D("h2d_temp","",nbins,low,high,nbinsy,lowy,highy);
  else h2d_temp= new TH2D("h2d_temp","",nbins,varbins,nbinsy,varbinsy);
  h2d_temp->Sumw2();

  if (dodata_) {
    if (hdata2d!=0) delete hdata2d;
    if (varbins==0 &&varbinsy==0)   hdata2d = new TH2D("hdata2d","",nbins,low,high,nbinsy,lowy,highy);
    else   hdata2d = new TH2D("hdata2d","",nbins,varbins,nbinsy,varbinsy);
    hdata2d->Sumw2();
  }

  TString opt="box";
  TString drawstring=vary;
  drawstring+=":";
  drawstring+=var;
  for (unsigned int isample=0; isample<samples_.size() ; isample++) {
    if ( isSampleSM(samples_[isample]) ) {
      gROOT->cd();
      TTree* tree = (TTree*) files_[currentConfig_][samples_[isample]]->Get("reducedTree");
      gROOT->cd();
      TString weightopt= useFlavorHistoryWeights_ && samples_[isample].Contains("WJets") ? "flavorHistoryWeight" : "";
      h2d_temp->Reset(); //just to be safe
      tree->Project("h2d_temp",drawstring,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,"",0,""));
      h2d->Add(h2d_temp);
    }
  }
  h2d->SetXTitle(xtitle);
  h2d->SetYTitle(ytitle);
  h2d->SetLineColor(sampleColor_["TotalSM"]);
  double zmax = h2d->GetMaximum();

  if (dodata_) {
    dtree->Project("hdata2d",drawstring, getCutString(true));
    hdata2d->SetLineColor(kBlack);
    if (hdata2d->GetMaximum() > zmax) zmax= hdata2d->GetMaximum();
  }

  h2d->Draw(opt);

  h2d->SetMaximum(zmax*1.1);
  hdata2d->SetMaximum(zmax*1.1);

  if (dodata_)  hdata2d->Draw("box same");


  TString savename = filename;
  if (logy_) savename += "-logY";

  savename += "-draw2d";

  //amazingly, \includegraphics cannot handle an extra dot in the filename. so avoid it.
  if (savePlots_) {
    thecanvas->SaveAs(savename+".eps"); //for me
    //  thecanvas->Print(savename+".C");    //for formal purposes
    thecanvas->SaveAs(savename+".pdf"); //for pdftex
    thecanvas->SaveAs(savename+".png"); //for twiki
  }

  delete h2d_temp;

}

// // overkill?
// class AxisInfo {
// public:
//   AxisInfo(int low, int high, TString title, TString unit="");
//   ~AxisInfo();

//   int low() const {return low_;}
//   int high() const {return high_;}
//   TString title() const {return title_;}
//   TString unit() const {return unit_;}

// private:
//   int low_;
//   int high_;
//   TString title_;
//   TString unit_;
// };
// AxisInfo::AxisInfo(int low, int high, TString title, TString unit) :
//   low_(low),
//   high_(high),
//   title(title_),
//   unit(unit_)
// {}
// AxisInfo::~AxisInfo() {}
TString appendBinWidth(const TString & ytitle, const double low,const double high,const int nbins,TString unit) {

  if (ytitle=="Arbitrary units") return ytitle; //oh this is so ugly

  double w = (high - low)/nbins;

  if (unit!="") unit.Prepend(" "); //if and only if unit is not null, prepend a space

  TString fulltitle;
  fulltitle.Form("%s/%.0f%s",ytitle.Data(),w,unit.Data());

  return fulltitle;
}

double dataIntegral_, MCIntegral_, MCIntegralErr_;
void drawPlots(const TString var, const int nbins, const float low, const float high, const TString xtitle, TString ytitle, TString filename="", const float* varbins=0, TString unit="") {
  //  cout<<"[drawPlots] var = "<<var<<endl;

  loadSamples();

  if (filename=="") filename=var;

  //  TH1D* thestackH=0;

  gROOT->SetStyle("CMS");
  //gStyle->SetHatchesLineWidth(1);

  TString canvasOpt = doRatio_ ? "ratio" : "";
  const int mainPadIndex = doRatio_ ? 1 : 0;
  renewCanvas(canvasOpt);

  thecanvas->cd(mainPadIndex);
  renewLegend();

  if (dostack_) {
    if (thestack!= 0 ) delete thestack;
    thestack = new THStack("thestack","--");
    if (doRatio_) {
      if (ratio!=0) delete ratio;
      ratio = (varbins==0) ? new TH1D("ratio","data/(SM MC)",nbins,low,high) : new TH1D("ratio","",nbins,varbins);
      ratio->Sumw2();

      ratioLine = new TLine(low, 1, high, 1);

    }
  }
  if (totalsm!=0) delete totalsm;
  totalsm = (varbins==0) ? new TH1D("totalsm","",nbins,low,high) : new TH1D("totalsm","",nbins,varbins);
  totalsm->Sumw2();
  if (totalsmsusy!=0) delete totalsmsusy;
  totalsmsusy = (varbins==0) ? new TH1D("totalsmsusy","",nbins,low,high) : new TH1D("totalsmsusy","",nbins,varbins);
  totalsmsusy->Sumw2();
  if (totalewk!=0) delete totalewk;
  totalewk = (varbins==0) ? new TH1D("totalewk","",nbins,low,high) : new TH1D("totalewk","",nbins,varbins);
  totalewk->Sumw2();
  if (totalqcdttbar!=0) delete totalqcdttbar;
  totalqcdttbar = (varbins==0) ? new TH1D("totalqcdttbar","",nbins,low,high) : new TH1D("totalqcdttbar","",nbins,varbins);
  totalqcdttbar->Sumw2();
  if (totalnonttbar!=0) delete totalnonttbar;
  totalnonttbar = (varbins==0) ? new TH1D("totalnonttbar","",nbins,low,high) : new TH1D("totalnonttbar","",nbins,varbins);
  totalnonttbar->Sumw2();
  if (totalnonqcd!=0) delete totalnonqcd;
  totalnonqcd = (varbins==0) ? new TH1D("totalnonqcd","",nbins,low,high) : new TH1D("totalnonqcd","",nbins,varbins);
  totalnonqcd->Sumw2();
  if (totalqcd!=0) delete totalqcd;
  totalqcd = (varbins==0) ? new TH1D("totalqcd","",nbins,low,high) : new TH1D("totalqcd","",nbins,varbins);
  totalqcd->Sumw2();

  totalsm->SetMarkerColor(sampleColor_["TotalSM"]);
  totalsm->SetLineColor(sampleColor_["TotalSM"]);
  totalsm->SetLineWidth(2);
  totalsm->SetMarkerStyle(sampleMarkerStyle_["TotalSM"]);
  if (!drawMarkers_)  totalsm->SetMarkerSize(0); //no marker for this one

  totalsmsusy->SetMarkerColor(sampleColor_["Total"]);
  totalsmsusy->SetLineColor(sampleColor_["Total"]);
  totalsmsusy->SetLineWidth(2);
  totalsmsusy->SetMarkerStyle(sampleMarkerStyle_["Total"]);
  if (!drawMarkers_)  totalsmsusy->SetMarkerSize(0); //no marker for this one

  if (varbins==0 && !renormalizeBins_) {ytitle = appendBinWidth(ytitle,low,high,nbins,unit);} //add the " / 10 GeV" part to the title

  //here is the part that is really different from the previous implementation
  //need to make new histograms
  resetHistos(); //delete existing histograms
  TString opt="hist e";
  double histMax=-1e9;
  for (unsigned int isample=0; isample<samples_.size(); isample++) {
    if (!quiet_)   cout <<samples_[isample]<<endl;

    gROOT->cd();
    //should each histo have a different name? maybe
    TString hname = jmt::fortranize(var); hname += "_"; hname += samples_[isample];
    histos_[samples_[isample]] = (varbins==0) ? new TH1D(hname,"",nbins,low,high) : new TH1D(hname,"",nbins,varbins);
    histos_[samples_[isample]]->Sumw2();

    TTree* tree = (TTree*) files_[currentConfig_][samples_[isample]]->Get("reducedTree");
    gROOT->cd();
    TString weightopt= useFlavorHistoryWeights_ && samples_[isample].Contains("WJets") ? "flavorHistoryWeight" : "";

    //treat ttbar in a special way 
    if(samples_[isample].Contains("TTbarJets-")){
      if(splitTTbarForClosureTest_){
	//tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,"",0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	if(samples_[isample].Contains("TTbarJets-semiMuGood")){
	  TString semiMuMode = "(((W1decayType==13 && W2decayType==112) || (W1decayType==112 && W2decayType==13) || (W1decayType==1513 && W2decayType==112) || (W1decayType==112 && W2decayType==1513) || (W1decayType==13 && W2decayType==134) || (W1decayType==134 && W2decayType==13) || (W1decayType==1513 && W2decayType==134) || (W1decayType==134 && W2decayType==1513))&&(decayType==101300))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiMuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-semiMuFailEta")){
	  TString semiMuMode = "(((W1decayType==13 && W2decayType==112) || (W1decayType==112 && W2decayType==13) || (W1decayType==1513 && W2decayType==112) || (W1decayType==112 && W2decayType==1513) || (W1decayType==13 && W2decayType==134) || (W1decayType==134 && W2decayType==13) || (W1decayType==1513 && W2decayType==134) || (W1decayType==134 && W2decayType==1513))&&(decayType==101301 || decayType==201301))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiMuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-semiMuFailPt")){
	  TString semiMuMode = "(((W1decayType==13 && W2decayType==112) || (W1decayType==112 && W2decayType==13) || (W1decayType==1513 && W2decayType==112) || (W1decayType==112 && W2decayType==1513) || (W1decayType==13 && W2decayType==134) || (W1decayType==134 && W2decayType==13) || (W1decayType==1513 && W2decayType==134) || (W1decayType==134 && W2decayType==1513))&&(decayType==101302 || decayType==201302))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiMuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-semiMuFailRecoIso")){
	  TString semiMuMode = "(((W1decayType==13 && W2decayType==112) || (W1decayType==112 && W2decayType==13) || (W1decayType==1513 && W2decayType==112) || (W1decayType==112 && W2decayType==1513) || (W1decayType==13 && W2decayType==134) || (W1decayType==134 && W2decayType==13) || (W1decayType==1513 && W2decayType==134) || (W1decayType==134 && W2decayType==1513))&&(decayType==201303))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiMuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-semiMuFailOther")){
	  TString semiMuMode = "(((W1decayType==13 && W2decayType==112) || (W1decayType==112 && W2decayType==13) || (W1decayType==1513 && W2decayType==112) || (W1decayType==112 && W2decayType==1513) || (W1decayType==13 && W2decayType==134) || (W1decayType==134 && W2decayType==13) || (W1decayType==1513 && W2decayType==134) || (W1decayType==134 && W2decayType==1513))&&(decayType==101304))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiMuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-semiEleGood")){
	  TString semiEleMode = "(((W1decayType==11 && W2decayType==112) || (W1decayType==112 && W2decayType==11) || (W1decayType==1511 && W2decayType==112) || (W1decayType==112 && W2decayType==1511) || (W1decayType==11 && W2decayType==134) || (W1decayType==134 && W2decayType==11) || (W1decayType==1511 && W2decayType==134) || (W1decayType==134 && W2decayType==1511))&&(decayType==101100))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiEleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-semiEleFailEta")){
	  TString semiEleMode = "(((W1decayType==11 && W2decayType==112) || (W1decayType==112 && W2decayType==11) || (W1decayType==1511 && W2decayType==112) || (W1decayType==112 && W2decayType==1511) || (W1decayType==11 && W2decayType==134) || (W1decayType==134 && W2decayType==11) || (W1decayType==1511 && W2decayType==134) || (W1decayType==134 && W2decayType==1511))&&(decayType==101101 || decayType==201101))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiEleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-semiEleFailPt")){
	  TString semiEleMode = "(((W1decayType==11 && W2decayType==112) || (W1decayType==112 && W2decayType==11) || (W1decayType==1511 && W2decayType==112) || (W1decayType==112 && W2decayType==1511) || (W1decayType==11 && W2decayType==134) || (W1decayType==134 && W2decayType==11) || (W1decayType==1511 && W2decayType==134) || (W1decayType==134 && W2decayType==1511))&&(decayType==101102 || decayType==201102))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiEleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-semiEleFailRecoIso")){
	  TString semiEleMode = "(((W1decayType==11 && W2decayType==112) || (W1decayType==112 && W2decayType==11) || (W1decayType==1511 && W2decayType==112) || (W1decayType==112 && W2decayType==1511) || (W1decayType==11 && W2decayType==134) || (W1decayType==134 && W2decayType==11) || (W1decayType==1511 && W2decayType==134) || (W1decayType==134 && W2decayType==1511))&&(decayType==201103))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiEleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-semiEleFailOther")){
	  TString semiEleMode = "(((W1decayType==11 && W2decayType==112) || (W1decayType==112 && W2decayType==11) || (W1decayType==1511 && W2decayType==112) || (W1decayType==112 && W2decayType==1511) || (W1decayType==11 && W2decayType==134) || (W1decayType==134 && W2decayType==11) || (W1decayType==1511 && W2decayType==134) || (W1decayType==134 && W2decayType==1511))&&(decayType==101104))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiEleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-semiTauHad")){
	  TString semiTauHadMode = "((W1decayType==15 && W2decayType==112) || (W1decayType==112 && W2decayType==15) || (W1decayType==15 && W2decayType==134) || (W1decayType==134 && W2decayType==15))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiTauHadMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-dilep")){
	  TString dilepMode = "((W1decayType==11 || W1decayType==13 || W1decayType==1511 || W1decayType==1513) && (W2decayType==11 || W2decayType==13 || W2decayType==1511 || W2decayType==1513))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,dilepMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-had")){
	  TString hadMode = "(W1decayType==112 || W1decayType==134) && (W2decayType==112 || W2decayType==134)";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,hadMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-other")){
	  //includes tautau(both had), etau(had), mutau(had), 
	  TString otherMode = "((W1decayType==15 && W2decayType==15)||(W1decayType==11 && W2decayType==15)||(W1decayType==15 && W2decayType==11)||(W1decayType==1511 && W2decayType==15)||(W1decayType==15 && W2decayType==1511)||(W1decayType==13 && W2decayType==15)||(W1decayType==15 && W2decayType==13)||(W1decayType==15 && W2decayType==1513)||(W1decayType==1513 && W2decayType==15))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,otherMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
      }
      else if(splitTTbar_){
	if(samples_[isample].Contains("TTbarJets-semiMu")){
	  TString semiMuMode = "((W1decayType==13 && W2decayType==112) || (W1decayType==112 && W2decayType==13) || (W1decayType==1513 && W2decayType==112) || (W1decayType==112 && W2decayType==1513) || (W1decayType==13 && W2decayType==134) || (W1decayType==134 && W2decayType==13) || (W1decayType==1513 && W2decayType==134) || (W1decayType==134 && W2decayType==1513))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiMuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-semiEle")){
	  TString semiEleMode = "((W1decayType==11 && W2decayType==112) || (W1decayType==112 && W2decayType==11) || (W1decayType==1511 && W2decayType==112) || (W1decayType==112 && W2decayType==1511) || (W1decayType==11 && W2decayType==134) || (W1decayType==134 && W2decayType==11) || (W1decayType==1511 && W2decayType==134) || (W1decayType==134 && W2decayType==1511))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiEleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-semiTauHad")){
	  TString semiTauHadMode = "((W1decayType==15 && W2decayType==112) || (W1decayType==112 && W2decayType==15) || (W1decayType==15 && W2decayType==134) || (W1decayType==134 && W2decayType==15))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiTauHadMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-dilep")){
	  TString dilepMode = "((W1decayType==11 || W1decayType==13 || W1decayType==1511 || W1decayType==1513) && (W2decayType==11 || W2decayType==13 || W2decayType==1511 || W2decayType==1513))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,dilepMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-had")){
	  TString hadMode = "((W1decayType==112 || W1decayType==134) && (W2decayType==112 || W2decayType==134))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,hadMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("TTbarJets-other")){
	  //includes tautau(both had), etau(had), mutau(had), 
	  TString otherMode = "((W1decayType==15 && W2decayType==15)||(W1decayType==11 && W2decayType==15)||(W1decayType==15 && W2decayType==11)||(W1decayType==1511 && W2decayType==15)||(W1decayType==15 && W2decayType==1511)||(W1decayType==13 && W2decayType==15)||(W1decayType==15 && W2decayType==13)||(W1decayType==15 && W2decayType==1513)||(W1decayType==1513 && W2decayType==15))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,otherMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}

      }
      else tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,"",0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
    }
    else if(samples_[isample].Contains("WJets-")){
      if(splitWJetsForClosureTest_){
	if(samples_[isample].Contains("WJets-muGood")){
	  TString MuMode = "((W1decayType==13 || W1decayType==1513) && (decayType==101300))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,MuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("WJets-muFailEta")){
	  TString MuMode = "((W1decayType==13 || W1decayType==1513) &&(decayType==101301 || decayType==201301))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,MuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("WJets-muFailPt")){
	  TString MuMode = "((W1decayType==13 || W1decayType==1513) &&(decayType==101302 || decayType==201302))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,MuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("WJets-muFailRecoIso")){
	  TString MuMode = "((W1decayType==13 || W1decayType==1513) &&(decayType==201303))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,MuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("WJets-muFailOther")){
	  TString MuMode = "((W1decayType==13 || W1decayType==1513) && (decayType==101304))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,MuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("WJets-eleGood")){
	  TString EleMode = "((W1decayType==11 || W1decayType==1511)&& (decayType==101100))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,EleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("WJets-eleFailEta")){
	  TString EleMode = "((W1decayType==11 || W1decayType==1511)&&(decayType==101101 || decayType==201101))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,EleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("WJets-eleFailPt")){
	  TString EleMode = "((W1decayType==11 || W1decayType==1511)&&(decayType==101102 || decayType==201102))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,EleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("WJets-eleFailRecoIso")){
	  TString EleMode = "((W1decayType==11 || W1decayType==1511) &&(decayType==201103))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,EleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("WJets-eleFailOther")){
	  TString EleMode = "((W1decayType==11 || W1decayType==1511) && (decayType==101104))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,EleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("WJets-tauHad")){
	  TString TauHadMode = "(W1decayType==15)";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,TauHadMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
      }
      else if(splitWJets_){
	if(samples_[isample].Contains("WJets-mu")){
	  TString MuMode = "(W1decayType==13 || W1decayType==1513)";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,MuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("WJets-ele")){
	  TString EleMode = "(W1decayType==11 || W1decayType==1511)";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,EleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("WJets-tauHad")){
	  TString TauHadMode = "(W1decayType==15)";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,TauHadMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
      }
      else tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,"",0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
    }
    else if(samples_[isample].Contains("SingleTop-sandtCombined")){
      if(splitSingleTopForClosureTest_){
	if(samples_[isample].Contains("SingleTop-sandtCombined-muGood")){
	  TString MuMode = "((W1decayType==13 || W1decayType==1513) && (decayType==101300))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,MuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-sandtCombined-muFailEta")){
	  TString MuMode = "((W1decayType==13 || W1decayType==1513) &&(decayType==101301 || decayType==201301))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,MuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-sandtCombined-muFailPt")){
	  TString MuMode = "((W1decayType==13 || W1decayType==1513) &&(decayType==101302 || decayType==201302))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,MuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-sandtCombined-muFailRecoIso")){
	  TString MuMode = "((W1decayType==13 || W1decayType==1513) &&(decayType==201303))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,MuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-sandtCombined-muFailOther")){
	  TString MuMode = "((W1decayType==13 || W1decayType==1513) && (decayType==101304))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,MuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-sandtCombined-eleGood")){
	  TString EleMode = "((W1decayType==11 || W1decayType==1511)&& (decayType==101100))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,EleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-sandtCombined-eleFailEta")){
	  TString EleMode = "((W1decayType==11 || W1decayType==1511)&&(decayType==101101 || decayType==201101))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,EleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-sandtCombined-eleFailPt")){
	  TString EleMode = "((W1decayType==11 || W1decayType==1511)&&(decayType==101102 || decayType==201102))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,EleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-sandtCombined-eleFailRecoIso")){
	  TString EleMode = "((W1decayType==11 || W1decayType==1511) &&(decayType==201103))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,EleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-sandtCombined-eleFailOther")){
	  TString EleMode = "((W1decayType==11 || W1decayType==1511) && (decayType==101104))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,EleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-sandtCombined-tauHad")){
	  TString TauHadMode = "(W1decayType==15)";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,TauHadMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
      }
      else tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,"",0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
    }

    else if(samples_[isample].Contains("SingleTop-tWCombined")){
      if(splitSingleTopForClosureTest_){
	if(samples_[isample].Contains("SingleTop-tWCombined-semiMuGood")){
	  TString semiMuMode = "(((W1decayType==13 && W2decayType==112) || (W1decayType==112 && W2decayType==13) || (W1decayType==1513 && W2decayType==112) || (W1decayType==112 && W2decayType==1513) || (W1decayType==13 && W2decayType==134) || (W1decayType==134 && W2decayType==13) || (W1decayType==1513 && W2decayType==134) || (W1decayType==134 && W2decayType==1513))&&(decayType==101300))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiMuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-tWCombined-semiMuFailEta")){
	  TString semiMuMode = "(((W1decayType==13 && W2decayType==112) || (W1decayType==112 && W2decayType==13) || (W1decayType==1513 && W2decayType==112) || (W1decayType==112 && W2decayType==1513) || (W1decayType==13 && W2decayType==134) || (W1decayType==134 && W2decayType==13) || (W1decayType==1513 && W2decayType==134) || (W1decayType==134 && W2decayType==1513))&&(decayType==101301 || decayType==201301))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiMuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-tWCombined-semiMuFailPt")){
	  TString semiMuMode = "(((W1decayType==13 && W2decayType==112) || (W1decayType==112 && W2decayType==13) || (W1decayType==1513 && W2decayType==112) || (W1decayType==112 && W2decayType==1513) || (W1decayType==13 && W2decayType==134) || (W1decayType==134 && W2decayType==13) || (W1decayType==1513 && W2decayType==134) || (W1decayType==134 && W2decayType==1513))&&(decayType==101302 || decayType==201302))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiMuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-tWCombined-semiMuFailRecoIso")){
	  TString semiMuMode = "(((W1decayType==13 && W2decayType==112) || (W1decayType==112 && W2decayType==13) || (W1decayType==1513 && W2decayType==112) || (W1decayType==112 && W2decayType==1513) || (W1decayType==13 && W2decayType==134) || (W1decayType==134 && W2decayType==13) || (W1decayType==1513 && W2decayType==134) || (W1decayType==134 && W2decayType==1513))&&(decayType==201303))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiMuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-tWCombined-semiMuFailOther")){
	  TString semiMuMode = "(((W1decayType==13 && W2decayType==112) || (W1decayType==112 && W2decayType==13) || (W1decayType==1513 && W2decayType==112) || (W1decayType==112 && W2decayType==1513) || (W1decayType==13 && W2decayType==134) || (W1decayType==134 && W2decayType==13) || (W1decayType==1513 && W2decayType==134) || (W1decayType==134 && W2decayType==1513))&&(decayType==101304))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiMuMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-tWCombined-semiEleGood")){
	  TString semiEleMode = "(((W1decayType==11 && W2decayType==112) || (W1decayType==112 && W2decayType==11) || (W1decayType==1511 && W2decayType==112) || (W1decayType==112 && W2decayType==1511) || (W1decayType==11 && W2decayType==134) || (W1decayType==134 && W2decayType==11) || (W1decayType==1511 && W2decayType==134) || (W1decayType==134 && W2decayType==1511))&&(decayType==101100))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiEleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-tWCombined-semiEleFailEta")){
	  TString semiEleMode = "(((W1decayType==11 && W2decayType==112) || (W1decayType==112 && W2decayType==11) || (W1decayType==1511 && W2decayType==112) || (W1decayType==112 && W2decayType==1511) || (W1decayType==11 && W2decayType==134) || (W1decayType==134 && W2decayType==11) || (W1decayType==1511 && W2decayType==134) || (W1decayType==134 && W2decayType==1511))&&(decayType==101101 || decayType==201101))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiEleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-tWCombined-semiEleFailPt")){
	  TString semiEleMode = "(((W1decayType==11 && W2decayType==112) || (W1decayType==112 && W2decayType==11) || (W1decayType==1511 && W2decayType==112) || (W1decayType==112 && W2decayType==1511) || (W1decayType==11 && W2decayType==134) || (W1decayType==134 && W2decayType==11) || (W1decayType==1511 && W2decayType==134) || (W1decayType==134 && W2decayType==1511))&&(decayType==101102 || decayType==201102))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiEleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-tWCombined-semiEleFailRecoIso")){
	  TString semiEleMode = "(((W1decayType==11 && W2decayType==112) || (W1decayType==112 && W2decayType==11) || (W1decayType==1511 && W2decayType==112) || (W1decayType==112 && W2decayType==1511) || (W1decayType==11 && W2decayType==134) || (W1decayType==134 && W2decayType==11) || (W1decayType==1511 && W2decayType==134) || (W1decayType==134 && W2decayType==1511))&&(decayType==201103))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiEleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-tWCombined-semiEleFailOther")){
	  TString semiEleMode = "(((W1decayType==11 && W2decayType==112) || (W1decayType==112 && W2decayType==11) || (W1decayType==1511 && W2decayType==112) || (W1decayType==112 && W2decayType==1511) || (W1decayType==11 && W2decayType==134) || (W1decayType==134 && W2decayType==11) || (W1decayType==1511 && W2decayType==134) || (W1decayType==134 && W2decayType==1511))&&(decayType==101104))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiEleMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-tWCombined-semiTauHad")){
	  TString semiTauHadMode = "((W1decayType==15 && W2decayType==112) || (W1decayType==112 && W2decayType==15) || (W1decayType==15 && W2decayType==134) || (W1decayType==134 && W2decayType==15))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,semiTauHadMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-tWCombined-dilep")){
	  TString dilepMode = "((W1decayType==11 || W1decayType==13 || W1decayType==1511 || W1decayType==1513) && (W2decayType==11 || W2decayType==13 || W2decayType==1511 || W2decayType==1513))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,dilepMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-tWCombined-had")){
	  TString hadMode = "(W1decayType==112 || W1decayType==134) && (W2decayType==112 || W2decayType==134)";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,hadMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
	if(samples_[isample].Contains("SingleTop-tWCombined-other")){
	  //includes tautau(both had), etau(had), mutau(had), 
	  TString otherMode = "((W1decayType==15 && W2decayType==15)||(W1decayType==11 && W2decayType==15)||(W1decayType==15 && W2decayType==11)||(W1decayType==1511 && W2decayType==15)||(W1decayType==15 && W2decayType==1511)||(W1decayType==13 && W2decayType==15)||(W1decayType==15 && W2decayType==13)||(W1decayType==15 && W2decayType==1513)||(W1decayType==1513 && W2decayType==15))";
	  tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,otherMode,0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
	}
      }
      else tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,"",0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    
    }

    else
      tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,"",0,"",-1,sampleScaleFactor_[samples_[isample]]).Data());    

    //now the histo is filled
    
    if (renormalizeBins_) ytitle=renormBins(histos_[samples_[isample]],2 ); //manipulates the TH1D //FIXME hard-coded "2"
    if (addOverflow_)  addOverflowBin( histos_[samples_[isample]] ); //manipulates the TH1D
    histos_[samples_[isample]]->SetXTitle(xtitle);
    histos_[samples_[isample]]->SetYTitle(ytitle);

    hinteractive = histos_[samples_[isample]];// hinteractive will point to the last sample's histo

    if (isSampleSM(samples_[isample])) {
      totalsm->Add(histos_[samples_[isample]]);
      //      if (!quiet_)    cout << "totalsm: " << samples_[isample] << endl;
    }
    if (!samples_[isample].Contains("LM") && !samples_[isample].Contains("QCD") && !samples_[isample].Contains("TTbar") && !samples_[isample].Contains("SUGRA")) {
      totalewk->Add(histos_[samples_[isample]]);
      //      if (!quiet_) cout << "totalewk: " << samples_[isample] << endl;
    }
    if (samples_[isample].Contains("QCD") || samples_[isample].Contains("TTbar")){
      totalqcdttbar->Add(histos_[samples_[isample]]);
      //      if (!quiet_) cout << "totalqcdttbar: " << samples_[isample] << endl;
    }
    if (!samples_[isample].Contains("TTbar") && !samples_[isample].Contains("LM") && !samples_[isample].Contains("SUGRA")) {
      totalnonttbar->Add(histos_[samples_[isample]]);
      //      if (!quiet_) cout << "totalnonttbar: " << samples_[isample] << endl;
    }
    if (!samples_[isample].Contains("QCD") && !samples_[isample].Contains("LM")&& !samples_[isample].Contains("SUGRA")){
       totalnonqcd->Add(histos_[samples_[isample]]);
       //      if (!quiet_) cout << "totalnonqcd: " << samples_[isample] << endl;
    }
    if (samples_[isample].Contains("QCD")){
       totalqcd->Add(histos_[samples_[isample]]);
       //      if (!quiet_) cout << "totalqcd: " << samples_[isample] << endl;
    }
    totalsmsusy->Add(histos_[samples_[isample]]); //add everything!

    //now just do a bunch of histogram formatting
    if (!dostack_) {
      //set line color instead of fill color for this type of plot
      histos_[samples_[isample]]->SetLineColor(sampleColor_[samples_[isample]]);
      histos_[samples_[isample]]->SetMarkerStyle(sampleMarkerStyle_[samples_[isample]]);
      histos_[samples_[isample]]->SetMarkerColor(sampleColor_[samples_[isample]]);
      if (!drawMarkers_) histos_[samples_[isample]]->SetMarkerSize(0);

      //ad hoc additions
      histos_[samples_[isample]]->SetLineWidth(2);
    }
    else {
      histos_[samples_[isample]]->SetFillColor(sampleColor_[samples_[isample]]);
      histos_[samples_[isample]]->SetMarkerSize(0);
    }

    if (dostack_) { //add histo to stack
      thestack->Add(histos_[samples_[isample]] );
    }
    else { //draw non-stacked histo
      //normalize
      if ( normalized_ && histos_[samples_[isample]]->Integral() >0)  histos_[samples_[isample]]->Scale( 1.0 / histos_[samples_[isample]]->Integral() );
      if (!drawSusyOnly_ || samples_[isample].Contains("LM")|| samples_[isample].Contains("SUGRA")) { //drawSusyOnly_ means don't draw SM
	//set max
	if ( findOverallMax( histos_[samples_[isample]]) > histMax) histMax = findOverallMax(histos_[samples_[isample]]);
       
	histos_[samples_[isample]]->Draw(opt);
	if (!opt.Contains("same")) opt+=" same";
      }
    }
  } //loop over samples and fill histograms

  //fill legend
  if (drawTotalSMSusy_) leg->AddEntry(totalsmsusy, sampleLabel_["Total"]);
  if (drawTotalSM_) leg->AddEntry(totalsm, sampleLabel_["TotalSM"]);
  for (int isample=int(samples_.size())-1; isample>=0; isample--) {
    if (dostack_ || (!drawSusyOnly_ || samples_[isample].Contains("LM")|| samples_[isample].Contains("SUGRA"))) { //drawSusyOnly_ means don't draw SM
      leg->AddEntry(histos_[samples_[isample]], sampleLabel_[samples_[isample]]);
    }
  }

  if (!dostack_) {
    //this is all a re-implemenataion of stuff done is HistHolder. Oh well.

    if (drawTotalSM_) histMax = totalsm->GetMaximum();
    for (unsigned int isample=0; isample<samples_.size(); isample++) {
      double pmx= doCustomPlotMax_ ? customPlotMax_ : histMax*maxScaleFactor_;
      histos_[samples_[isample]]->SetMaximum(pmx);
      if (doCustomPlotMin_) histos_[samples_[isample]]->SetMinimum(customPlotMin_);
    }

    if (drawTotalSM_) { 
      totalsm->Draw(opt); 
      if (doCustomPlotMax_) totalsm->SetMaximum(customPlotMax_);
    }
    if (drawTotalSMSusy_) {
      totalsmsusy->Draw(opt); 
      if (doCustomPlotMax_) totalsmsusy->SetMaximum(customPlotMax_);
    }
  }
  else {
    thestack->Draw("hist");
    thestack->GetHistogram()->GetXaxis()->SetTitle(xtitle);
    thestack->GetHistogram()->GetYaxis()->SetTitle(ytitle);

    if (doVerticalLine_) drawVerticalLine(); //i want to draw the data last

    if (drawMCErrors_) {
      if (mcerrors!=0) delete mcerrors;
      mcerrors = new TGraphErrors(totalsm);
      mcerrors->SetFillStyle(3353); //3353 3544
      mcerrors->SetFillColor(1);
      //ack. TGraphs and TH1s use different conventions for numbering.
      for ( int ibin=1; ibin<=totalsm->GetNbinsX(); ibin++) {
	double yerr = mcerrors->GetErrorY(ibin-1);
	double xerr = totalsm->GetBinCenter(ibin) - totalsm->GetBinLowEdge(ibin);
	mcerrors->SetPointError(ibin-1,xerr,yerr);
      }
      mcerrors->Draw("2 same");
    }

    if (doCustomPlotMax_) thestack->SetMaximum(customPlotMax_);
    if (doCustomPlotMin_) thestack->SetMinimum(customPlotMin_);
  } //if doStack_

  if (dodata_) {
    gROOT->cd();
    //    if (!quiet_)     cout<<"Drawing data!"<<endl;
    if (hdata != 0) delete hdata;
    TString hname = jmt::fortranize(var); hname += "_"; hname += "data";
    hdata = (varbins==0) ? new TH1D(hname,"",nbins,low,high) : new TH1D(hname,"",nbins,varbins);
    gROOT->cd();
    dtree->Project(hname,var,getCutString(true));

    hdata->SetBinErrorOption(TH1::kPoisson);
    //now the histo is filled
    
    hdata->SetMarkerColor(kBlack);
    hdata->SetLineWidth(2);
    hdata->SetMarkerStyle(kFullCircle);
    hdata->SetMarkerSize(1);
    if (renormalizeBins_) renormBins(hdata,2 ); //manipulates the histogram //FIXME hard-coded "2"
    if (addOverflow_)     addOverflowBin(hdata); // manipulates the histogram!
    leg->AddEntry(hdata,"Data");

    if (!quiet_)    cout<<"Data underflow: " <<hdata->GetBinContent(0)<<endl;//BEN
    hdata->Draw("SAME E");
    if (!doCustomPlotMax_) {
      double mymax=-1e9;
      if (dostack_) mymax = thestack->GetMaximum();
      if ( findOverallMax(totalsm) >mymax) mymax = findOverallMax(totalsm);
      if (findOverallMax(hdata) > mymax) mymax = findOverallMax(hdata);
      if (dostack_) thestack->SetMaximum( maxScaleFactor_*mymax);
      else { //i don't like repeating this loop constantly; at a minimum it should be abstracted
	for (unsigned int isample=0; isample<samples_.size(); isample++)  histos_[samples_[isample]]->SetMaximum(maxScaleFactor_*mymax);
      }
    }


    if (!quiet_ && !renormalizeBins_) {
      cout<<"Integral of data, EW, total SM: "<<hdata->Integral()<<" ; "<<totalewk->Integral()<<" ; "<<totalsm->Integral()<<endl;
      cout<<"Chi^2 Test results: "<<hdata->Chi2Test(totalsm,"UW P")<<endl;
      cout<<"KS Test results: "<<hdata->KolmogorovTest(totalsm,"N")<<endl;;
      cout<<"KS Test (without \"N\") results: "<<hdata->KolmogorovTest(totalsm)<<endl;

      dataIntegral_ = hdata->Integral();
      MCIntegral_ = totalsm->IntegralAndError(totalsm->GetXaxis()->GetFirst(),totalsm->GetXaxis()->GetLast(),MCIntegralErr_);
    }
    if (doRatio_) {
      thecanvas->cd(2);
      ratio->Divide(hdata,totalsm);
      ratio->SetMinimum(ratioMin);
      ratio->SetMaximum(ratioMax);

      ratio->GetYaxis()->SetNdivisions(200 + int(ratioMax-ratioMin)+1);    //set ticks ; to be seen if this really works
      //ratio->GetYaxis()->SetNdivisions(300 + 5);    //float ratioMin=0.5; float ratioMax=1.5;
      //ratio->GetYaxis()->SetLabelSize(0.2); //make y label bigger
      ratio->GetYaxis()->SetLabelSize(0.1); //make y label bigger
      ratio->GetXaxis()->SetLabelSize(0.1); //make y label bigger
      ratio->GetXaxis()->SetTitleOffset(1.1);
      ratio->GetXaxis()->SetTitle(xtitle); //make y label bigger
      ratio->GetXaxis()->SetLabelSize(0.12);
      ratio->GetXaxis()->SetLabelOffset(0.04);
      ratio->GetXaxis()->SetTitleSize(0.16);
      ratio->GetYaxis()->SetTitle("Data/MC ");
      ratio->GetYaxis()->SetTitleSize(0.16);
      ratio->GetYaxis()->SetTitleOffset(.43);
      gPad->SetTopMargin(1e-5);
      gPad->SetTickx();
      gPad->Modified();      

      ratio->Draw();

      //thecanvas->GetPad(2)->SetTopMargin(0.1);

      ratioLine->Draw();

    }
  }
  
  if (!quiet_ && dostack_ && dodata_ && nbins<11 && doRatio_) {//BEN - 11 is an arbitrary number that isn't too big so we don't print out too much stuff.
    for(int i=1; i<=nbins; i++){
      cout << "data: " << hdata->GetBinContent(i) << " +- " << hdata->GetBinError(i) << ", totalsm: " << totalsm->GetBinContent(i) << " +- " << totalsm->GetBinError(i) << ", ratio: " << ratio->GetBinContent(i) << " +- " << ratio->GetBinError(i) << endl;
    }
  }
  

  thecanvas->cd(mainPadIndex);
  if (doleg_)  leg->Draw();
  drawPlotHeader();

  //  if (doSubtraction_) savename+="-MCSub";
  TString savename = filename;
  if (logy_) savename += "-logY";
  //  savename += scaleAppendToFilename;

  if (!dostack_ && !normalized_)      savename += "-drawPlain";
  else if (!dostack_ && normalized_)  savename += "-drawNorm";
  else savename += "-drawStack";

  //amazingly, \includegraphics cannot handle an extra dot in the filename. so avoid it.
  if (savePlots_) {
    thecanvas->SaveAs(savename+".eps"); //for me
    //  thecanvas->Print(savename+".C");    //for formal purposes
    thecanvas->SaveAs(savename+".pdf"); //for pdftex
    thecanvas->SaveAs(savename+".png"); //for twiki
  }

  //dump some event counts to the screen
  if (!quiet_) {
    for (unsigned int isample=0; isample<samples_.size(); isample++) {
      cout<<samples_[isample]<<" =\t "<<histos_[samples_[isample]]->Integral()<<" +/- "<<jmt::errOnIntegral(histos_[samples_[isample]])<<endl;
    }
    cout<<"total SM =\t "<<totalsm->Integral()<<" +/- "<<jmt::errOnIntegral(totalsm)<<endl;
  }
  
}

void drawPlots(const TString var, const int nbins, const float* varbins, const TString xtitle, const TString ytitle, TString filename="") {
  //provide a more natural modification to the argument list....
  drawPlots( var, nbins, 0, 1, xtitle,ytitle, filename,  varbins);
}

TH1D* getHist(const TString & sample) {
  TH1D* h=0;
  if (sample=="totalsm") h=totalsm;
  else if (sample=="data") h=hdata;
  else h = histos_[sample];

  return h;
}

double getIntegral(const TString & sample) {
  return getHist(sample)->Integral();
}

double getIntegralErr(const TString & sample) {
  return jmt::errOnIntegral(getHist(sample));
}

double getSumOfIntegrals(const std::vector<TString> & samples) {
  double sum=0;
  for (unsigned int j= 0; j<samples.size(); j++) {
    sum += getIntegral(samples.at(j));
  }
  return sum;
}

double getSumOfIntegralsErr(const std::vector<TString> & samples) {
  double sum=0;
  for (unsigned int j= 0; j<samples.size(); j++) {
    sum += pow(getIntegralErr(samples.at(j)),2);
  }
  return sqrt(sum);
}

void drawSignificance(const TString & var, const int nbins, const float low, const float high, const TString & savename) {

  bool oldSaveSetting = savePlots_;
  savePlots_=false;
  drawPlots(var,nbins,low,high,var,"",savename);

  TH1D* hsusy=0;
  for (std::vector<TString>::iterator it = samples_.begin(); it!=samples_.end(); ++it) {
    if ( (*it).Contains("LM") ||(*it).Contains("SUGRA")  ) hsusy = histos_[ *it];
  }
  if (hsusy==0) {
    cout<<"Didn't find a signal sample"<<endl;
    return;
  }

  for (int ibin= 1; ibin<=nbins; ibin++) {
    double B = totalsm->Integral(ibin,nbins);
    double S = hsusy->Integral(ibin,nbins);
    if (B>0)    cout<<totalsm->GetBinLowEdge(ibin)<<"\t"<<S/sqrt(B)<<endl;//"\t"<<totalsm->GetBinContent(ibin)<<endl;
    else     cout<<ibin<<" B is zero"<<endl;
  }



  savePlots_=oldSaveSetting;
}

//could add xtitle and ytitle
void drawR(const TString vary, const float cutVal, const TString var, const int nbins, const float low, const float high, const TString& savename, const float* varbins=0, bool dataOnly=false) {
  const TString ytitle="N pass / N fail";
  TString xtitle = var;
  if(var=="MET") xtitle = "E_{T}^{miss} [GeV]"; 

  cout << "x axis: " << var << endl;
  TString cstring1 = vary, cstring2=vary;
  cstring1 += " >= ";
  cstring2 += " < ";
  cstring1 += cutVal;
  cstring2 += cutVal;

  loadSamples();
  
  gROOT->SetStyle("CMS");
  
  TString opt=doRatio_? "ratio":"";
  renewCanvas(opt);
  
  resetHistos(); //delete existing histograms

  resetLegendPositionR();
  renewLegend();

  if(dodata_){
    if (hdata != 0) delete hdata;
    
    TString hname = var; hname += "_"; hname += "data";
    hdata = (varbins==0) ? new TH1D(hname,"",nbins,low,high) : new TH1D(hname,"",nbins,varbins);
    hdata->Sumw2();
    
    TString hnameP = var; hnameP += "_"; hnameP += "dataPass";
    histos_[hnameP] = (varbins==0) ? new TH1D(hnameP,"",nbins,low,high) : new TH1D(hnameP,"",nbins,varbins);
    histos_[hnameP]->Sumw2();
    
    TString hnameF = var; hnameF += "_"; hnameF += "dataFail";
    histos_[hnameF] = (varbins==0) ? new TH1D(hnameF,"",nbins,low,high) : new TH1D(hnameF,"",nbins,varbins);
    histos_[hnameF]->Sumw2();
    
    gROOT->cd();
    dtree->Project(hnameP,var,getCutString(kData,"",selection_,cstring1).Data());
    dtree->Project(hnameF,var,getCutString(kData,"",selection_,cstring2).Data());
    if (addOverflow_)  addOverflowBin( histos_[hnameP] );
    if (addOverflow_)  addOverflowBin( histos_[hnameF] );
    //compute ratio
    hdata->Divide(histos_[hnameP], histos_[hnameF]);

    hdata->SetMarkerColor(kBlack);
    hdata->SetLineWidth(2);
    hdata->SetMarkerStyle(kFullCircle);
    hdata->SetMarkerSize(1);
    hdata->SetYTitle(ytitle);
    hdata->SetXTitle(xtitle);
    
    leg->AddEntry(hdata,"Data");
  }
  if(dataOnly){
    if(doRatio_ || dodata_) cout << "drawR inconsistent arguments" << endl;
    if (doCustomPlotMax_) {
      hdata->SetMaximum(customPlotMax_);
    }
    if (doCustomPlotMin_) {
      hdata->SetMinimum(customPlotMin_);
    }
    gROOT->cd();//what does this do?
    thecanvas->cd();
    gPad->SetRightMargin(0.1);
    gPad->Modified();
    hinteractive = hdata;
    hdata->Draw();
    drawPlotHeader(-.1);
    if (doleg_)  leg->Draw();
    thecanvas->SaveAs("mindpPassOverFail-"+savename+".eps");
    thecanvas->SaveAs("mindpPassOverFail-"+savename+".pdf");
    thecanvas->SaveAs("mindpPassOverFail-"+savename+".png");
    return;
  }


  TH1D * totalsm_pass = (varbins==0) ? new TH1D("totalsm_pass","",nbins,low,high) : new TH1D("totalsm_pass","",nbins,varbins);
  TH1D * totalsm_fail = (varbins==0) ? new TH1D("totalsm_fail","",nbins,low,high) : new TH1D("totalsm_fail","",nbins,varbins);
  totalsm_pass->Sumw2(); 
  totalsm_fail->Sumw2(); 
  if (totalsm!=0) delete totalsm;
  totalsm =  (varbins==0) ? new TH1D("totalsm","",nbins,low,high) : new TH1D("totalsm","",nbins,varbins);
  totalsm->Sumw2();

  totalsm->SetMarkerColor(sampleColor_["TotalSM"]);
  totalsm->SetLineColor(sampleColor_["TotalSM"]);
  totalsm->SetLineWidth(2);
  totalsm->SetMarkerStyle(0);
  totalsm->SetYTitle(ytitle);

  TString drawopt="hist e";
  float max=-1e9; TString firsthist="";
  for (unsigned int isample=0; isample<samples_.size(); isample++) {

    if (!quiet_) cout <<samples_[isample]<<endl;
    TTree* tree = (TTree*) files_[currentConfig_][samples_[isample]]->Get("reducedTree");

    gROOT->cd();

    //need Pass, Fail, and Ratio for each sample
    TString hnameP = var; hnameP += "_"; hnameP += samples_[isample];
    hnameP += "_Pass";
    histos_[hnameP] = (varbins==0) ? new TH1D(hnameP,"",nbins,low,high) : new TH1D(hnameP,"",nbins,varbins);
    histos_[hnameP]->Sumw2();

    TString hnameF = var; hnameF += "_"; hnameF += samples_[isample];
    hnameF += "_Fail";
    histos_[hnameF] = (varbins==0) ? new TH1D(hnameF,"",nbins,low,high) : new TH1D(hnameF,"",nbins,varbins);
    histos_[hnameF]->Sumw2();

    TString hnameR = var; hnameR += "_"; hnameR += samples_[isample];
    hnameR += "_Ratio";
    histos_[hnameR] = (varbins==0) ? new TH1D(hnameR,"",nbins,low,high) : new TH1D(hnameR,"",nbins,varbins);
    histos_[hnameR]->Sumw2();

    
    //Fill histos
    if (useFlavorHistoryWeights_) assert(0); // this needs to be implemented
    tree->Project(hnameP,var,getCutString(getSampleType(samples_[isample],"point"),"",selection_,cstring1).Data());
    tree->Project(hnameF,var,getCutString(getSampleType(samples_[isample],"point"),"",selection_,cstring2).Data());

    if (addOverflow_)  addOverflowBin( histos_[hnameP] );
    if (addOverflow_)  addOverflowBin( histos_[hnameF] );

    cout<<"Bug check"<<endl;
    for (int ib=1; ib<3; ib++) {
      cout<<histos_[hnameP]->GetBinContent(ib)<<"\t";
      cout<<histos_[hnameF]->GetBinContent(ib)<<"\t";
      cout<<histos_[hnameR]->GetBinContent(ib)<<endl;
    }
    //compute ratio
    histos_[hnameR]->Divide(histos_[hnameP], histos_[hnameF]);

    if (isSampleSM(samples_[isample])) {
      totalsm_pass->Add(histos_[hnameP]);
      totalsm_fail->Add(histos_[hnameF]);
    }

    //   cout<<"content of bin 2: "<<histos_[hnameP]->GetBinContent(2)<<" / "<< histos_[hnameF]->GetBinContent(2)<<" = "<<histos_[hnameR]->GetBinContent(2)<<endl;

    //now format the histograms
    if (!quiet_) cout<<"setting color to: "<<sampleColor_[samples_[isample]]<<endl;
    histos_[hnameR]->SetLineColor(sampleColor_[samples_[isample]]);
    histos_[hnameR]->SetMarkerStyle(sampleMarkerStyle_[samples_[isample]]);
    histos_[hnameR]->SetMarkerColor(sampleColor_[samples_[isample]]);
    histos_[hnameR]->SetYTitle(ytitle);
    histos_[hnameR]->SetXTitle(xtitle);
    histos_[hnameR]->GetYaxis()->SetLabelSize(0.04); //make y label bigger
    histos_[hnameR]->GetXaxis()->SetLabelSize(0.04); //make x label bigger

    //ad hoc additions
    histos_[hnameR]->SetLineWidth(2);
    
    //draw
    thecanvas->cd(1);
    gPad->SetRightMargin(0.1);
    gPad->Modified();
    
    if(dodata_) drawPlotHeader(-.1);
    else  drawPlotHeaderInside(); //for qcd only plots, this works.
    
    if (hnameR.Contains("QCD")) { //HACK draw only qcd
      histos_[hnameR]->Draw(drawopt);
      if (!drawopt.Contains("same")) drawopt+=" same";
      
      if (firsthist=="") firsthist = hnameR;
      if (histos_[hnameR]->GetMaximum() > max) max = histos_[hnameR]->GetMaximum();
      leg->AddEntry(histos_[hnameR], sampleLabel_[samples_[isample]]);
    }
    
  }
  
  histos_[firsthist]->SetMaximum( max*maxScaleFactor_);
  hinteractive =  histos_[firsthist];
  
  totalsm->Divide(totalsm_pass,totalsm_fail);
  if (drawTotalSM_) {
    totalsm->Draw("hist e same");
    leg->AddEntry(totalsm,sampleLabel_["TotalSM"]);
  }
  
  if (dodata_) {
    gROOT->cd();
    if (!quiet_)   cout<<"Drawing data!"<<endl;
    
    thecanvas->cd(1);
    hdata->Draw("SAME");
    
    if (hdata->GetMaximum() > max)  {
      histos_[firsthist]->SetMaximum( maxScaleFactor_*hdata->GetMaximum());
      totalsm->SetMaximum(maxScaleFactor_*hdata->GetMaximum());
    }
    if (doCustomPlotMax_) {
      histos_[firsthist]->SetMaximum( customPlotMax_);
      totalsm->SetMaximum(customPlotMax_);
    }
    if (doCustomPlotMin_) {
      histos_[firsthist]->SetMinimum( customPlotMin_);
      totalsm->SetMinimum(customPlotMin_);
    }
    
    if (doRatio_) {
      thecanvas->cd(2);
      gPad->SetRightMargin(0.1);
      gPad->Modified();
      if (ratio!=0) delete ratio;
      ratio = (varbins==0) ? new TH1D("ratio","data/(SM MC)",nbins,low,high) : new TH1D("ratio","data/(SM MC)",nbins,varbins);
      ratio->Sumw2();
      ratio->Divide(hdata,totalsm); 
      ratio->SetMinimum(ratioMin);
      ratio->SetMaximum(ratioMax);
      ratio->GetYaxis()->SetNdivisions(200 + int(ratioMax-ratioMin)+1);    //set ticks ; to be seen if this really works
      ratio->GetYaxis()->SetLabelSize(0.1); //make y label bigger
      ratio->GetXaxis()->SetLabelSize(0.1); //make x label bigger
      ratio->Draw();
      //cout<<"KS Test results (shape only): "<<hdata->KolmogorovTest(totalsm)<<endl;;
    }
    
  } // end of if (dodata_)
  else {
    //    cout<<"setting custom plot ranges"<<endl;//jmt debug
    cout<<firsthist<<endl;
    cout<<histos_[firsthist]<<endl;
    if (doCustomPlotMax_) {
      histos_[firsthist]->SetMaximum( customPlotMax_);
      totalsm->SetMaximum(customPlotMax_);
    }
    if (doCustomPlotMin_) {
      histos_[firsthist]->SetMinimum( customPlotMin_);
      totalsm->SetMinimum(customPlotMin_);
    }
  }
  
  //  cout<<"about to cd canvas and draw leg"<<endl;//jmt debug
  
  if (doRatio_)  thecanvas->cd(1);
  if (doleg_)  leg->Draw();
  
  thecanvas->SaveAs("mindpPassOverFail-"+savename+".eps");
  thecanvas->SaveAs("mindpPassOverFail-"+savename+".pdf");
  thecanvas->SaveAs("mindpPassOverFail-"+savename+".png");

  delete totalsm_pass;
  delete totalsm_fail;

}

void drawR(const TString vary, const float cutVal, const int nbins, const float low, const float high, const TString& savename) {//for backwards compatibility
  drawR(vary, cutVal, "MET", nbins, low, high, savename, 0); 
}

void drawR(const TString vary, const float cutVal, const TString var, const int nbins, const float* varbins=0, const TString& savename="") {//more natural
  drawR(vary, cutVal, var, nbins, 0, 1, savename, varbins);
}



//typedef map<pair<int,int>, pair<float,float> > susyScanYields;
typedef map<pair<int,int>, float > susyScanYields;
susyScanYields getSusyScanYields(const TString & sampleOfInterest,const int pdfindex=0, const TString & pdfset="CTEQ") {
  cout<<" ~~ begin slow version of getSusyScanYields()"<<endl;

  /*
this function was written to draw mSugra. It can also draw T1bbbb with no problem.
The pdfindex and pdfset can also be specified, but it is not good for checking the whole mSugra plane because it is too slow
  */

  if (sampleOfInterest.Contains("mSUGRA") ) assert( pdfindex==0); //just until i have time to think about it

  TStopwatch timer;

  /*
this function was written for mSugra.
most of it is irrelevant for SMS, but we'll use it anyway
  */
  susyScanYields theYields;

  //sample is an argument
  //other important things are defined by the usual global variables
  if (!quiet_) cout<<sampleOfInterest<<" "<<currentConfig_<<endl;

  //these are for mSugra
  TString varx="m0"; TString xtitle=varx;
  int  nbinsx=310; float lowx=-0.5; float highx=3100-0.5; //UPDATED for new scans

  TString vary="m12"; TString ytitle=vary;
  int  nbinsy=110; float lowy=-0.5; float highy=1100-0.5;

  if (sampleOfInterest== "T1bbbb" || sampleOfInterest== "T2bb" || sampleOfInterest== "T2tt" || sampleOfInterest== "T1tttt") { //change binning
    //m0 -> mGl
    //m12 -> mLSP
    nbinsx=60;//scanSMSngen->GetNbinsX();
    nbinsy=60;//scanSMSngen->GetNbinsY(); //take our histo def'n from the ngen histo
    lowx=0;//scanSMSngen->GetXaxis()->GetBinLowEdge(1);
    lowy=0;//scanSMSngen->GetYaxis()->GetBinLowEdge(1);
    highx=1500;//scanSMSngen->GetXaxis()->GetBinLowEdge(nbinsx+1);
    highy=1500;//scanSMSngen->GetYaxis()->GetBinLowEdge(nbinsy+1);
  }

  TString drawstring = vary+":"+varx;

  TTree* thetree = (TTree*) files_[currentConfig_][sampleOfInterest]->Get("reducedTree");

  const int nsubprocesses = sampleOfInterest.Contains("mSUGRA") ? 10 : 0;
  vector<TH2D*> raw0;
  for (int i=0; i<=nsubprocesses; i++) {
    TString hname="raw0_";
    hname += i;
    raw0.push_back(new TH2D(hname,"raw event counts",nbinsx,lowx,highx,nbinsy,lowy,highy));
    raw0[i]->Sumw2();
    TString thecut = sampleOfInterest.Contains("mSUGRA") ? getCutString(kmSugraPlane,"",selection_,"",pdfindex,pdfset,i) : getCutString(kSMSPlane,"",selection_,"",pdfindex,pdfset);
    thetree->Project(hname,drawstring,thecut.Data());
  }


  cout<<"[getSusyScanYields (slow)] done counting cross-section weighted events"<<endl;
  //[for mSugra]
  //at this point, each bin contains Npass_i * sigma_i for that (m0,m12)
  //need to divide by N_i for each (m0,m12)

  //[for SMS]
  //each bin contains just raw Npass_i for that (mgl,mLSP)

  if (sampleOfInterest.Contains("mSUGRA") ) {
    //loop over i and histo bins
    for (int i=0; i<=nsubprocesses; i++) {
      for (map<pair<int,int>, map<TString,TH2D*> >::iterator iscanpoint = scanProcessTotalsMap.begin(); iscanpoint!=scanProcessTotalsMap.end(); ++iscanpoint) {
	int m0=iscanpoint->first.first;
	int m12=iscanpoint->first.second;
	TH2D* thishist = scanProcessTotalsMap[make_pair(m0,m12)][pdfset];
	//thisn does not have to be integer in the case of pdf weights
	double thisn = thishist->GetBinContent(i,pdfindex); //bin 0 has the raw total events (no pdf weights)
	int bin=  raw0[i]->FindBin(m0,m12);
	double N_i_thispoint = raw0[i]->GetBinContent(bin);
	double err_i_thispoint = raw0[i]->GetBinError(bin);
	if (thisn < 0.0000001) {
	  if (N_i_thispoint > 0.0000001) cout<<"Possible problem: "<<m0<<" "<<m12<<" "<<i<<" "<< N_i_thispoint<<" "<<thisn<<endl;
	  thisn=1; //prevent divide by zero
	  N_i_thispoint = 0; //need to come back to what is going wrong h
	}
	N_i_thispoint /= thisn;
	err_i_thispoint /= thisn;
	raw0[i]->SetBinContent(bin,N_i_thispoint);
	raw0[i]->SetBinError(bin,err_i_thispoint);
      }
    }
    
    cout<<"[getSusyScanYields (slow)] done dividing by Ngen"<<endl;
    //now we have Npass_i * sigma_i * lumi / Ngen_i 
    //all that is left is to make the sum over i
    
    for (map<pair<int,int>, map<TString,TH2D*> >::iterator iscanpoint = scanProcessTotalsMap.begin(); iscanpoint!=scanProcessTotalsMap.end(); ++iscanpoint) {
      double Nraw = 0, errraw=0;
      
      int bin=  raw0[0]->FindBin(iscanpoint->first.first , iscanpoint->first.second);
      for (unsigned int i=0; i<raw0.size(); i++) {
	Nraw += raw0[i]->GetBinContent(bin);
	errraw += pow(raw0[i]->GetBinError(bin),2);
      }
      //cout<<iscanpoint->first.first<<" "<<iscanpoint->first.second<<" "<<Nraw<< " +/- "<<sqrt(errraw)<<endl;
      theYields[iscanpoint->first] = Nraw;//make_pair(Nraw,sqrt(errraw));
    }
    cout<<"[getSusyScanYields (slow)] done summing over subprocesses"<<endl;
  }
  else { //not mSugra
    //just copy the histogram into the data structure
    //kinda stupid to do this instead of just returning a histogram, but such is life
    for (int i=1; i<=nbinsx; i++) {
      for (int j=1; j<=nbinsy; j++) {
	int mgl=TMath::Nint(raw0[0]->GetXaxis()->GetBinLowEdge(i));
	int mlsp=TMath::Nint(raw0[0]->GetYaxis()->GetBinLowEdge(j));
	double nevents = raw0[0]->GetBinContent(raw0[0]->FindBin(mgl,mlsp));
	//throw out points that have ngen too far from 10000
	//note that scanSMSngen could be replaced by scanProcessTotalsMap, but it is easier and safer not to change this for now
	if ( TMath::Nint(scanSMSngen->GetBinContent(scanSMSngen->FindBin(mgl,mlsp))) >= 1000 ) {
	  theYields[make_pair(mgl,mlsp)] = nevents;//make_pair(nevents,0); //skip errors for now
	}
      }
    }
  }

  //try to clean up
  cout<<"[getSusyScanYields (slow)] cleaning up memory"<<endl;
  for (unsigned int i=0; i<raw0.size(); i++) {
    delete raw0[i];
  }

  timer.Stop();
  cout<<"[getSusyScanYields() (slow)] CPU time = "<<timer.CpuTime()<<endl;

  return theYields;
}

//sigh...looks like to use Luke's magic code i have to really go at this from scratch
//this works for T1bbbb. Need to work on it to make it work for mSugra
vector<susyScanYields> getSusyScanYields(const TString & sampleOfInterest,const TString & pdfset) {
  cout<<" ~~ begin fast version of getSusyScanYields()"<<endl;
  TStopwatch timer;

  vector<susyScanYields> theYields;
  if (!quiet_) cout<<sampleOfInterest<<" "<<currentConfig_<<endl;
  //these are for mSugra
  TString varx="m0"; TString xtitle=varx;
  int  nbinsx=310; float lowx=-0.5; float highx=3100-0.5; //UPDATED for new scans

  TString vary="m12"; TString ytitle=vary;
  int  nbinsy=110; float lowy=-0.5; float highy=1100-0.5;
  if (sampleOfInterest== "T1bbbb" || sampleOfInterest== "T2bb" || sampleOfInterest== "T2tt" || sampleOfInterest== "T1tttt") { //change binning
    //m0 -> mGl
    //m12 -> mLSP
    if (scanSMSngen==0) scanSMSngen = (TH2D*) files_[currentConfig_][sampleOfInterest]->Get("scanSMSngen");
    nbinsx=60;//scanSMSngen->GetNbinsX();
    nbinsy=60;//scanSMSngen->GetNbinsY(); //take our histo def'n from the ngen histo
    lowx=0;//scanSMSngen->GetXaxis()->GetBinLowEdge(1);
    lowy=0;//scanSMSngen->GetYaxis()->GetBinLowEdge(1);
    highx=1500;//scanSMSngen->GetXaxis()->GetBinLowEdge(nbinsx+1);
    highy=1500;//scanSMSngen->GetYaxis()->GetBinLowEdge(nbinsy+1);
  }
  TString drawstring = vary+":"+varx;

  TTree* thetree = (TTree*) files_[currentConfig_][sampleOfInterest]->Get("reducedTree");
  //key of the map is the susysubprocess. the vector is over the pdf weight indices
  map<int, vector<TH2F*> > raw0;
  int npdfset=0;
  if (pdfset=="CTEQ") npdfset=44;
  else if (pdfset=="MSTW") npdfset =40;
  else if (pdfset=="NNPDF") {npdfset =99; } //not implemented yet
  else assert(0);

  //work luke's magic
  TSelectorMultiDraw* multiDraw = new TSelectorMultiDraw();

  //for mSugra this is crazy. i need to loop over both pdf sets and subprocesses
  const int nsubprocesses = sampleOfInterest.Contains("mSUGRA") ? 10 : 0;
  for (int isusy=0; isusy<=nsubprocesses; isusy++) {
    for (int ipdf=0; ipdf<=npdfset; ipdf++) { //start at 0 instead of 1...i think that's right
      TString hname="raw0_";
      hname += ipdf;
      hname+="_";
      hname+=isusy;
      raw0[isusy].push_back(new TH2F(hname,"raw event counts",nbinsx,lowx,highx,nbinsy,lowy,highy));
      raw0[isusy][ipdf]->Sumw2();
      TString thecut = sampleOfInterest.Contains("mSUGRA") ? getCutString(kmSugraPlane,"",selection_,"",ipdf,pdfset,isusy) : getCutString(kSMSPlane,"",selection_,"",ipdf,pdfset);
      multiDraw->LoadVariables( TString(drawstring+">>"+hname).Data(), thecut.Data());
    }
  }
  Long64_t numberofentries = thetree->GetEntries();
  thetree->Process(multiDraw,"goff",numberofentries,0);


  //copy the histograms into the data structure
  //kinda stupid to do this instead of just returning the histograms, but such is life
  for (int ipdf=0; ipdf<=npdfset; ipdf++) { //start at 0 instead of 1...i think that's right
    susyScanYields oneSetOfYields;

    if (sampleOfInterest.Contains("mSUGRA")) {
      //[for mSugra]
      //at this point, each bin contains Npass_i * sigma_i * lumi for that (m0,m12)
      //need to divide by N_i for each (m0,m12,ipdf)
      for (int i=0; i<=nsubprocesses; i++) {
	for (map<pair<int,int>, map<TString,TH2D*> >::iterator iscanpoint = scanProcessTotalsMap.begin(); iscanpoint!=scanProcessTotalsMap.end(); ++iscanpoint) {
	  int m0=iscanpoint->first.first;
	  int m12=iscanpoint->first.second;
	  TH2D* thishist = scanProcessTotalsMap[make_pair(m0,m12)][pdfset];
	  double thisn = thishist->GetBinContent(i,ipdf); //ngen is not necessarily an integer! 
	  int bin=  raw0[i][ipdf]->FindBin(m0,m12);
	  double N_i_thispoint = raw0[i][ipdf]->GetBinContent(bin);
	  double err_i_thispoint = raw0[i][ipdf]->GetBinError(bin);
	  if (thisn < 0.0000001) {
	    if (N_i_thispoint > 0.0000001) cout<<"[new getSusyScanYields()] Possible problem: "<<m0<<" "<<m12<<" "<<i<<" "<<ipdf<<" "<< N_i_thispoint<<" "<<thisn<<endl;
	    thisn=1; //prevent divide by zero
	    N_i_thispoint = 0; 
	  }
	  N_i_thispoint /= thisn;
	  err_i_thispoint /= thisn;
	  raw0[i][ipdf]->SetBinContent(bin,N_i_thispoint);
	  raw0[i][ipdf]->SetBinError(bin,err_i_thispoint);
	}
      }
      //now we have Npass_i * sigma_i * lumi / Ngen_i 
      //all that is left is to make the sum over i
      for (map<pair<int,int>, map<TString,TH2D*> >::iterator iscanpoint = scanProcessTotalsMap.begin(); iscanpoint!=scanProcessTotalsMap.end(); ++iscanpoint) {
	double Nraw = 0, errraw=0;
	
	int bin=  raw0[0][ipdf]->FindBin(iscanpoint->first.first , iscanpoint->first.second);
	//since we're using a map, we could do an iterator loop here, but this works too
	for (unsigned int i=0; i<raw0.size(); i++) { //we want the size of the *map* (number of susy subprocesses)
	  Nraw += raw0[i][ipdf]->GetBinContent(bin);
	  errraw += pow(raw0[i][ipdf]->GetBinError(bin),2);
	}
	//cout<<iscanpoint->first.first<<" "<<iscanpoint->first.second<<" "<<Nraw<< " +/- "<<sqrt(errraw)<<endl;
	oneSetOfYields[iscanpoint->first] = Nraw;//make_pair(Nraw,sqrt(errraw));
      }
    }
    else { //for T1bbbb and hopefully other samples that don't have the susy subprocesses 
      for (int i=1; i<=nbinsx; i++) {
	for (int j=1; j<=nbinsy; j++) {
	  int mgl=TMath::Nint(raw0[0][ipdf]->GetXaxis()->GetBinLowEdge(i));
	  int mlsp=TMath::Nint(raw0[0][ipdf]->GetYaxis()->GetBinLowEdge(j));
	  double nevents = raw0[0][ipdf]->GetBinContent(raw0[0][ipdf]->FindBin(mgl,mlsp));
	  if ( TMath::Nint(scanSMSngen->GetBinContent(scanSMSngen->FindBin(mgl,mlsp))) >= 1000 ) {
	    //need to multiply by N_gen_nominal / N_gen_ipdf
	    //N_gen_nominal should match between scanSMSngen and scanProcessTotals!
	    TH2D* thishist = scanProcessTotalsMap[make_pair(mgl,mlsp)][pdfset];
	    assert( thishist->GetBinContent(10,0) == TMath::Nint(scanSMSngen->GetBinContent(scanSMSngen->FindBin(mgl,mlsp))) );
	    //this cout was a useful check, but it's a *lot* of output
	    //	    cout<<"   -- Rescaling PDF weights by "<< thishist->GetBinContent(10,0)<<" / "<< thishist->GetBinContent(10,ipdf)<<endl;
	    nevents *= thishist->GetBinContent(10,0) / thishist->GetBinContent(10,ipdf);

	    oneSetOfYields[make_pair(mgl,mlsp)] = nevents;//make_pair(nevents,0); //skip errors for now
	  }
	}
      }
    }

    theYields.push_back(oneSetOfYields);
  } //end loop over pdf index
  timer.Stop();
  cout<<"[getSusyScanYields() (fast)] CPU time ("<<pdfset<<") = "<<timer.CpuTime()<<endl;

  //try to clean up
  for ( map<int, vector<TH2F*> >::iterator im = raw0.begin(); im!=raw0.end() ; ++im) {
    for (unsigned int i=0; i<im->second.size(); i++) {
      delete im->second[i];
    }
  }
  delete multiDraw;
  cout<<"[getSusyScanYields() (fast)] done with cleanup"<<endl;

  return theYields;

}

//   XplusMax          XminusMax
pair<susyScanYields,susyScanYields> doScanPDFUncertainties(const TString & sample, const TString & pdfset,  susyScanYields & nominal) {
  //can't make nominal const because of the use of [] . should fix it

  int npdfset=0;
  if (pdfset=="CTEQ") npdfset=44;
  else if (pdfset=="MSTW") npdfset =40;
  else if (pdfset=="NNPDF") {npdfset =99; }
  else assert(0);

  vector<susyScanYields> theYields;
  theYields=getSusyScanYields(sample,pdfset);

  //  pair<int,int> apoint=make_pair(375,100);
//   for (unsigned int i=0; i<theYields.size(); i++) {
//     cout<<"DEBUG: "<<nominal[apoint].first<<"\t"<<theYields[i][apoint].first<<endl;
//   }
  
  //in the current implementation, for T1bbbb we do all the calculations above and then
  //just shuffle the results into two vectors below

  //in mSugra we get all the results below

  vector<susyScanYields> cteqXplus;
  vector<susyScanYields> cteqXminus;
  for (int i=1; i<=npdfset; i++) {
    //for the scans, the pdfWeight in the ntuple is already normalized properly
    //i filled element 0 of theYields on purpose so that I could skip it here

    //this is old logic, which assumed that sometimes we would need to do the calculation one pdf index at a time
    //no need to change it though
    susyScanYields cteq= theYields.empty() ?  getSusyScanYields(sample,i,pdfset) : theYields.at(i);
    if (pdfset=="NNPDF") { //for nnpdf, no need to split things up
      cteqXplus.push_back(cteq); 
    }
    else {
      if (i%2==0) {
	//	cout<<" "<<cteq[apoint].first<<endl;
	cteqXminus.push_back(cteq);
      }
      else {
	//	cout<<"DEBUG "<<cteq[apoint].first;
	cteqXplus.push_back(cteq);
      }
    }
  }
  assert((cteqXplus.size() == cteqXminus.size()) || pdfset=="NNPDF");
  susyScanYields XplusMax ;
  susyScanYields XminusMax;
  if (pdfset=="NNPDF" ) {

    int nn=0;
    for (unsigned int i=0; i<cteqXplus.size(); i++) { //only xplus is filled
      for (susyScanYields::iterator ipoint=cteqXplus[i].begin(); ipoint!=cteqXplus[i].end(); ++ipoint) {
	//if this scan point doesn't exist, create it with 0
	if ( XplusMax.count(ipoint->first) ==0) XplusMax[ipoint->first]=0;//make_pair(0,0);
	if ( XminusMax.count(ipoint->first) ==0) XminusMax[ipoint->first]=0;//make_pair(0,0);
	double val = cteqXplus[i][ipoint->first];//.first;
	XplusMax[ipoint->first] += val; //sum up the values for pdf i and this scan point
	XminusMax[ipoint->first] += val*val; //sum up the values for pdf i and this scan point
      }
      nn++;
    }
    for (susyScanYields::iterator ipoint=XplusMax.begin() ; ipoint!=XplusMax.end(); ++ipoint) {
      XplusMax[ipoint->first] = XplusMax[ipoint->first] / double(nn);
      XminusMax[ipoint->first] = XminusMax[ipoint->first] / double(nn);
    }

  }
  else {
    //loop over pdfs and over scan points
    for (unsigned int i=0; i<cteqXplus.size(); i++) {
      for (susyScanYields::iterator ipoint=cteqXplus[i].begin(); ipoint!=cteqXplus[i].end(); ++ipoint) {
	//	if (ipoint->first ==apoint) cout<<ipoint->first.first<<" "<<ipoint->first.second<<" -- "<<nominal[ipoint->first].first<<" " <<cteqXplus[i][ipoint->first].first<<" "<<cteqXminus[i][ipoint->first].first<<endl;
	double diff1 =  cteqXplus[i][ipoint->first] - nominal[ipoint->first] ;
	double diff2 =  cteqXminus[i][ipoint->first] - nominal[ipoint->first] ;
	double larger = diff1>diff2 ? diff1 : diff2;
	if ( 0 > larger ) larger = 0;
	if ( XplusMax.count(ipoint->first) ==0) {
	  XplusMax[ipoint->first]=0;//make_pair(0,0);
	  //  if (ipoint->first ==apoint) cout<<" ***** CREATING XplusMax"<<endl;
	}
	XplusMax[ipoint->first] += larger*larger;
	//	if (ipoint->first ==apoint) cout<<"       XplusMax = "<<XplusMax[ipoint->first].first<<endl;
      }
    }
    for (unsigned int i=0; i<cteqXplus.size(); i++) {
      for (susyScanYields::iterator ipoint=cteqXplus[i].begin(); ipoint!=cteqXplus[i].end(); ++ipoint) {
	double diff1 =  nominal[ipoint->first]  - cteqXplus[i][ipoint->first];
	double diff2 =  nominal[ipoint->first] - cteqXminus[i][ipoint->first] ;
	double larger = diff1>diff2 ? diff1 : diff2;
	if ( 0 > larger ) larger = 0;
	if ( XminusMax.count(ipoint->first) ==0) XminusMax[ipoint->first]=0;//make_pair(0,0);
	XminusMax[ipoint->first] += larger*larger;
      }
    }

    const double scale = (pdfset=="CTEQ") ? 1.645 : 1;
    for (susyScanYields::iterator ipoint=XplusMax.begin() ; ipoint!=XplusMax.end() ; ++ipoint) {
      ipoint->second = sqrt(ipoint->second);
      ipoint->second = ipoint->second / scale;
    }
    for (susyScanYields::iterator ipoint=XminusMax.begin() ; ipoint!=XminusMax.end() ; ++ipoint) {
      ipoint->second = sqrt(ipoint->second);
      ipoint->second = ipoint->second / scale;
    }
  }



  return make_pair(XplusMax,XminusMax);

}

void getCutStringForCutflow(vector<TString> &vectorOfCuts, vector<TString> &stageCut, bool isTightSelection, bool btagSF=false) {

  //careful -- doesn't support mSUGRA right now

  vectorOfCuts.clear(); stageCut.clear();

  //define the HT and MET thresholds
  TString minHT; TString minMET;
  if (!isTightSelection) {minHT = "400"; minMET = "250";}
  else {minHT = "500"; minMET = "300";}

  TString cut;
  TString thisSelection;
  TString selectionPreB;
  TString selectionEq1bLoose;
  TString selectionEq2bLoose;
  TString selectionGe3bLoose;  
  TString selectionGe1bLoose;
  TString selectionGe2bLoose;
  TString selectionEq1bTightHT;
  TString selectionEq2bTightHT;  
  TString selectionGe3bTightHT;  
  TString selectionGe1bTightHT;
  TString selectionGe2bTightHT;  
  TString selectionEq1bTightMET;
  TString selectionEq2bTightMET;  
  TString selectionGe3bTightMET;  
  TString selectionGe1bTightMET;
  TString selectionGe2bTightMET;  

  TString selectionGe1bTight;
  TString selectionGe2bTight;  

  //inclusive
  thisSelection="";
  cut=getCutString(kMC,thisSelection);
  vectorOfCuts.push_back(cut);
  stageCut.push_back("Inclusive");

  //trigger
  if (thisSelection=="") thisSelection += "cutTrigger==1";
  else thisSelection += " && cutTrigger";
  cut=getCutString(kMC,thisSelection);
  vectorOfCuts.push_back(cut);
  stageCut.push_back("Trigger");
  
  //PV
  if (thisSelection=="") thisSelection += "cutPV==1";
  else thisSelection += " && cutPV==1";
  cut=getCutString(kMC,thisSelection);
  vectorOfCuts.push_back(cut);
  stageCut.push_back("PV");

  //HT
  if (thisSelection=="") {thisSelection += "HT>="; thisSelection += minHT;}
  else { thisSelection += " && HT>="; thisSelection += minHT;}
  cut=getCutString(kMC,thisSelection); //not compatible with scans!
  vectorOfCuts.push_back(cut);
  if (latexMode_) stageCut.push_back("HT$\\ge$"+minHT);
  else stageCut.push_back("HT>="+minHT);
  
  //MET
  if (thisSelection=="") {thisSelection += "MET>="; thisSelection += minMET; }
  else {thisSelection += " && MET>="; thisSelection += minMET;}
  cut=getCutString(kMC,thisSelection);
  vectorOfCuts.push_back(cut);
  if (latexMode_) stageCut.push_back("\\MET$\\ge$"+minMET);
  else stageCut.push_back("MET>="+minMET);

  //3 or more jets
  if (thisSelection=="") thisSelection += "cut3Jets==1";
  else thisSelection += " && cut3Jets==1";
  cut=getCutString(kMC,thisSelection);
  vectorOfCuts.push_back(cut);
  if (latexMode_) stageCut.push_back("$\\ge$ 3 jets");
  else stageCut.push_back(">= 3 jets");

  //Ele veto
  if (thisSelection=="") thisSelection += "cutEleVeto==1";
  else thisSelection += " && cutEleVeto==1";
  cut=getCutString(kMC,thisSelection);
  vectorOfCuts.push_back(cut);
  stageCut.push_back("e veto");

  //Mu veto
  if (thisSelection=="") thisSelection += "cutMuVeto==1";
  else thisSelection += " && cutMuVeto==1";
  cut=getCutString(kMC,thisSelection);
  vectorOfCuts.push_back(cut);
  if (latexMode_) stageCut.push_back("$\\mu$ veto");
  else stageCut.push_back("Mu veto");

  //angular cuts

  //  //deltaPhi
  //  if (thisSelection=="") thisSelection += "cutDeltaPhi==1";
  //  else thisSelection += " && cutDeltaPhi==1";
  //  cut=getCutString(lumiScale_,thisSelection);
  //  vectorOfCuts.push_back(cut);
  //  stageCut.push_back("DeltaPhi");

  //deltaPhiN
  //  if (thisSelection=="") thisSelection += "cutDeltaPhiN==1";
  //  else thisSelection += " && cutDeltaPhiN==1";
  if (thisSelection=="") thisSelection += "minDeltaPhiN>4";
  else thisSelection += " && minDeltaPhiN>=4";
  cut=getCutString(kMC,thisSelection);
  vectorOfCuts.push_back(cut);
  if (latexMode_) stageCut.push_back("$\\minDeltaPhiN>$4");
  else stageCut.push_back("minDeltaPhiN>4");

  //  //deltaPhi(MET,taus)
  //  if (thisSelection=="") thisSelection += "cutDeltaPhiTaus==1";
  //  else thisSelection += " && cutDeltaPhiTaus==1";
  //  cut=getCutString(lumiScale_,thisSelection);
  //  vectorOfCuts.push_back(cut);
  //  stageCut.push_back("DeltaPhiTaus");
  
  //Cleaning
  if (thisSelection=="") thisSelection += "passCleaning==1";
  else thisSelection +=" && passCleaning==1";
  cut=getCutString(kMC,thisSelection);
  vectorOfCuts.push_back(cut);
  stageCut.push_back("Cleaning");

  //store selection string pre b cut
  selectionPreB=thisSelection;

  //loose selection
  ////>= 1 b
  selectionGe1bLoose=selectionPreB;
  if (btagSF) btagSFweight_="probge1";
  else selectionGe1bLoose += " && nbjetsCSVM>=1";
  cut=getCutString(kMC,selectionGe1bLoose);
  vectorOfCuts.push_back(cut);
  if (latexMode_) stageCut.push_back("HT$\\ge$400, \\MET$\\ge$250, $\\ge$1 b");
  else stageCut.push_back("HT>=400, MET>=250, >= 1 b");


  //==1b
/* shape
  selectionEq1bLoose=selectionPreB; 
  if (btagSF) btagSFweight_="prob1";
  else selectionEq1bLoose +=" && nbjetsCSVM==1";
  cut=getCutString(kMC,selectionEq1bLoose);
  vectorOfCuts.push_back(cut);
  if (latexMode_) stageCut.push_back("HT$\\ge$400, \\MET$\\ge$250, $==$1 b");
  else stageCut.push_back("HT>=400, MET>=250, == 1 b");
*/

  ////>= 2 b
  selectionGe2bLoose=selectionPreB; 
  if (btagSF) btagSFweight_="probge2";
  else selectionGe2bLoose += " && nbjetsCSVM>=2";
  cut=getCutString(kMC,selectionGe2bLoose);
  vectorOfCuts.push_back(cut);
  if (latexMode_) stageCut.push_back("HT$\\ge$400, \\MET$\\ge$250, $\\ge$2 b");
  else stageCut.push_back("HT>=400, MET>=250, >= 2 b");

  //== 2 b
  /* shape
  selectionEq2bLoose=selectionPreB; 
  if (btagSF) btagSFweight_="(1-prob1-prob0-probge3)";
  else selectionEq2bLoose += " && nbjetsCSVM==2";
  cut=getCutString(kMC,selectionEq2bLoose);
  vectorOfCuts.push_back(cut);
  if (latexMode_) stageCut.push_back("HT$\\ge$400, \\MET$\\ge$250, $==$2 b");
  else stageCut.push_back("HT>=400, MET>=250, == 2 b");
*/

  //>= 3 b
  selectionGe3bLoose=selectionPreB; 
  if (btagSF) btagSFweight_="probge3";
  else selectionGe3bLoose += " && nbjetsCSVM>=3";
  cut=getCutString(kMC,selectionGe3bLoose);
  vectorOfCuts.push_back(cut);
  if (latexMode_) stageCut.push_back("HT$\\ge$400, \\MET$\\ge$250, $\\ge$3 b");
  else stageCut.push_back("HT>=400, MET>=250, >= 3 b");


  //tight selection
  if (!isTightSelection){ //print out the results for the tight selection anyway

    ////OLD 1BT and 2BT cuts
    ////>=1 b
    selectionGe1bTight=selectionPreB; 
    if (btagSF) {btagSFweight_="probge1"; selectionGe1bTight += " && HT>=500 && MET>=500";}
    else selectionGe1bTight += " && nbjetsCSVM>=1 && HT>=500 && MET>=500"; //hard-coded!
    cut=getCutString(kMC,selectionGe1bTight);
    vectorOfCuts.push_back(cut);
    if (latexMode_) stageCut.push_back("HT$\\ge$500, \\MET$\\ge$500, $\\ge$1 b");
    else stageCut.push_back("HT>=500, MET>=500, >= 1 b");
    ////>=2 b
    selectionGe1bTight=selectionPreB;
    if (btagSF) {btagSFweight_="probge2"; selectionGe1bTight += " && HT>=600 && MET>=300";}
    else selectionGe1bTight += " && nbjetsCSVM>=2 && HT>=600 && MET>=300"; //hard-coded!
    cut=getCutString(kMC,selectionGe1bTight);
    vectorOfCuts.push_back(cut);
    if (latexMode_) stageCut.push_back("HT$\\ge$600, \\MET$\\ge$300, $\\ge$2 b");
    else stageCut.push_back("HT>=600, MET>=300, >= 2 b");  

    //TIGHT MET
    //==1 b
/*  shape
   selectionEq1bTightMET=selectionPreB; 
    if (btagSF) {btagSFweight_="prob1"; selectionEq1bTightMET += " && HT>=500 && MET>=500";}
    else selectionEq1bTightMET += " && nbjetsCSVM==1 && HT>=500 && MET>=500"; //hard-coded!
    cut=getCutString(kMC,selectionEq1bTightMET);
    vectorOfCuts.push_back(cut);
    if (latexMode_) stageCut.push_back("HT$\\ge$500, \\MET$\\ge$500, $==$1 b");
    else stageCut.push_back("HT>=500, MET>=500, == 1 b");
    //==2 b
    selectionEq2bTightMET=selectionPreB; 
    if (btagSF) {btagSFweight_="(1-prob1-prob0-probge3)"; selectionEq2bTightMET += " && HT>=500 && MET>=500";}
    else selectionEq2bTightMET += " && nbjetsCSVM==2 && HT>=500 && MET>=500"; //hard-coded!
    cut=getCutString(kMC,selectionEq2bTightMET);
    vectorOfCuts.push_back(cut);
    if (latexMode_) stageCut.push_back("HT$\\ge$500, \\MET$\\ge$500, $==$2 b");
    else stageCut.push_back("HT>=500, MET>=500, == 2 b");
    //>=3 b
    selectionGe3bTightMET=selectionPreB; 
    if (btagSF) {btagSFweight_="probge3"; selectionGe3bTightMET += " && HT>=500 && MET>=500";}
    else selectionGe3bTightMET += " && nbjetsCSVM>=3 && HT>=500 && MET>=500"; //hard-coded!
    cut=getCutString(kMC,selectionGe3bTightMET);
    vectorOfCuts.push_back(cut);
    if (latexMode_) stageCut.push_back("HT$\\ge$500, \\MET$\\ge$500, $\\ge$3 b");
    else stageCut.push_back("HT>=500, MET>=500, >= 3 b");
        
    //TIGHT HT
    //==1 b
    selectionEq1bTightHT=selectionPreB; 
    if (btagSF) {btagSFweight_="prob1"; selectionEq1bTightHT += " && HT>=600 && MET>=300";}
    else selectionEq1bTightHT += " && nbjetsCSVM==1 && HT>=600 && MET>=300"; //hard-coded!
    cut=getCutString(kMC,selectionEq1bTightHT);
    vectorOfCuts.push_back(cut);
    if (latexMode_) stageCut.push_back("HT$\\ge$600, \\MET$\\ge$300, $==$1 b");
    else stageCut.push_back("HT>=600, MET>=300, == 1 b");
    //==2 b
    selectionEq2bTightHT=selectionPreB; 
    if (btagSF) {btagSFweight_="(1-prob1-prob0-probge3)"; selectionEq2bTightHT += " && HT>=600 && MET>=300";}
    else selectionEq2bTightHT += " && nbjetsCSVM==2 && HT>=600 && MET>=300"; //hard-coded!
    cut=getCutString(kMC,selectionEq2bTightHT);
    vectorOfCuts.push_back(cut);
    if (latexMode_) stageCut.push_back("HT$\\ge$600, \\MET$\\ge$300, $==$2 b");
    else stageCut.push_back("HT>=600, MET>=300, == 2 b");
    //>=3 b
    selectionGe3bTightHT=selectionPreB; 
    if (btagSF) {btagSFweight_="probge3"; selectionGe3bTightHT += " && HT>=600 && MET>=300";}
    else selectionGe3bTightHT += " && nbjetsCSVM>=3 && HT>=600 && MET>=300"; //hard-coded!
    cut=getCutString(kMC,selectionGe3bTightHT);
    vectorOfCuts.push_back(cut);
    if (latexMode_) stageCut.push_back("HT$\\ge$600, \\MET$\\ge$300, $\\ge$3 b");
    else stageCut.push_back("HT>=600, MET>=300, >= 3 b");
     */   
  }
  
  //  //check output
  //  for (unsigned int istage=0; istage<vectorOfCuts.size(); istage++){
  //    cout<<stageCut[istage]<<endl;
  //    cout<<vectorOfCuts[istage]<<endl;
  //  }
      
}

void cutflow(bool isTightSelection){

  loadSamples();
  resetHistos();

  //jmt -- turn on weighting
  usePUweight_=true; 
  useHTeff_=true;   
  useMHTeff_=true;
  currentConfig_=configDescriptions_.getCorrected(); //use JERbias

  vector<TString> vectorOfCuts; //each element is a successive cut string for the cutflow table
  vector<TString> stageCut; //name of cut at each stage
  //the 'true' at the end turns on b tag SF
  getCutStringForCutflow(vectorOfCuts, stageCut, isTightSelection,true); //fills vectorOfCuts

  vector< vector<float> > cutflowEntries; //stores events by stage and sample
  vector< vector<float> > cutflowEntriesE; //stores error

  TString var = "HT";
  const float* varbins=0;
  const int nbins=1; const float low=0; const float high=100000000000;

  vector<float> nSumSM;
  vector<float> nSumSME;

  //loop over cuts
  for (unsigned int istage=0; istage<vectorOfCuts.size(); istage++){
    TString thisStageCut = vectorOfCuts[istage];
    resetHistos();

    vector<float> nPass; //number of events passing current cut in each sample
    vector<float> nPassE; //error on events passing cuts    

    //for total background
    if (totalsm!=0) delete totalsm;
    totalsm = (varbins==0) ? new TH1D("totalsm","",nbins,low,high) : new TH1D("totalsm","",nbins,varbins);
    totalsm->Sumw2();

    //loop over samples, copied from drawPlots()
    for (unsigned int isample=0; isample<samples_.size(); isample++) {
      if (!quiet_)   cout <<samples_[isample]<<endl;

      gROOT->cd();
      //should each histo have a different name? maybe
      TString hname = jmt::fortranize(var); hname += "_"; hname += samples_[isample];
      histos_[samples_[isample]] = (varbins==0) ? new TH1D(hname,"",nbins,low,high) : new TH1D(hname,"",nbins,varbins);
      histos_[samples_[isample]]->Sumw2();

      TTree* tree = (TTree*) files_[currentConfig_][samples_[isample]]->Get("reducedTree");
      gROOT->cd();
      TString weightopt= useFlavorHistoryWeights_ && samples_[isample].Contains("WJets") ? "flavorHistoryWeight" : "";
      tree->Project(hname,var,thisStageCut);
      //now the histo is filled

      if (isSampleSM(samples_[isample])) {
	totalsm->Add(histos_[samples_[isample]]);
	if (!quiet_)    cout << "totalsm: " << samples_[isample] << endl;
      }

      nPass.push_back(histos_[samples_[isample]]->GetBinContent(1)); //to get total entries with weighting applied
      nPassE.push_back(histos_[samples_[isample]]->GetBinError(1)); //error
            
    }//end loop over samples

    cutflowEntries.push_back(nPass);
    cutflowEntriesE.push_back(nPassE);
    nSumSM.push_back(totalsm->GetBinContent(1));
    nSumSME.push_back(totalsm->GetBinError(1));

  }//end current cut 

  //now print out the cutflow table
  if (isTightSelection) cout<<"tight selection (HT>=500, MET>=300)"<<endl;
  else cout<<"baseline selection (HT>=400, MET>=250)"<<endl; 

  TString col_start; TString col; TString col_end; TString hline; TString hhline;

  if (latexMode_){
    col_start=""; col=" & "; col_end=" \\\\ "; hline="\\hline"; hhline="\\hline \\hline";
  }
  else{
    col_start=" | "; col=" | "; col_end=" | "; hline=""; hhline="";
  }

  cout<<hhline<<endl;
  //list sample names
  cout<<col_start<<"Cut "<<col;
  for (unsigned int isample=0; isample<samples_.size(); isample++) {
    if ( isSampleSM(samples_[isample])) { //skip LM points to add total background column first 
      if ( samples_[isample].Contains("VV") ) cout<<"Diboson"<<col; 
      else if ( samples_[isample].Contains("QCD") ) cout<<"QCD"<<col;
      else if ( samples_[isample]=="SingleTop" ) cout<<"Single Top"<<col;
      else if ( samples_[isample].Contains("TTbarJets") ){
	if (latexMode_) cout<<"\\ttbar"<<col;
	else cout<<samples_[isample]<<col;
      }
      else if ( samples_[isample].Contains("WJets") ){
	if (latexMode_)cout<<"\\WJets"<<col;
	else cout<<samples_[isample]<<col;
      }
      else if ( samples_[isample].Contains("ZJets") ){
	if (latexMode_)cout<<"\\ZJets"<<col;
        else cout<<samples_[isample]<<col;
      }
      else if ( samples_[isample].Contains("Zinvisible") ){
        if (latexMode_)cout<<"\\Zinvisible"<<col;
        else cout<<samples_[isample]<<col;
      }
      else cout<<samples_[isample]<<col;
    }
  }
  //now add total background, LM
  cout<<"Total SM";
  for (unsigned int isample=0; isample<samples_.size(); isample++) {
    if ( !isSampleSM(samples_[isample])  ) cout<<col<<samples_[isample];
  }
  cout<<col_end<<endl;
  
  //fill table
  for (unsigned int istage=0; istage<vectorOfCuts.size(); istage++){
    cout<<col_start<<stageCut[istage]<<col;
    for (unsigned int isample=0; isample<samples_.size(); isample++){
      if ( isSampleSM(samples_[isample]) ) {
	cout<<jmt::format_nevents(cutflowEntries[istage][isample],cutflowEntriesE[istage][isample])<<col;
      }
    }//end loop over samples
    //now add total background, LM
    cout<<jmt::format_nevents(nSumSM[istage],nSumSME[istage]);
    for (unsigned int isample=0; isample<samples_.size(); isample++) {
      if ( !isSampleSM(samples_[isample]) ) cout<<col<<jmt::format_nevents(cutflowEntries[istage][isample],cutflowEntriesE[istage][isample]);
    }
    cout<<col_end<<endl;
    if (stageCut[istage]=="Cleaning") cout<<hline<<endl;

    if (latexMode_){
      if (stageCut[istage]=="HT$\\ge$400, \\MET$\\ge$250, $\\ge$3 b" && !isTightSelection) cout<<hline<<endl;
    }
    else {
      if (stageCut[istage]=="HT>=400, MET>=250, >= 2 b" && !isTightSelection) cout<<hline<<endl;
    }

  }//end loop over cuts
  cout<<hhline<<endl;
  
}
