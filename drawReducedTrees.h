// -*- C++ -*-

#include "TStopwatch.h"

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
TChain* dtree=0;
TH1D* hdata=0;

TString currentConfig_;

//default selection
TString selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1";

float leg_x1 , leg_x2, leg_y1, leg_y2;
void resetLegendPosition() {
  leg_x1 = 0.696;
  leg_x2=0.94;
  leg_y1=0.5;
  leg_y2=0.92;
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

  TString htSelection;
  TString metSelection;
  TString btagSelection;

  TString owenId;
  bool isSIG;
};
SearchRegion::SearchRegion(TString btagSel,TString htSel,TString metSel,TString oId,bool isSig) : 
  htSelection(htSel),metSelection(metSel),btagSelection(btagSel),owenId(oId),isSIG(isSig) {}
SearchRegion::~SearchRegion() {}
void SearchRegion::Print() const {
  cout<<" == "<<btagSelection<<" "<<htSelection<<" "<<metSelection<<endl;

}

std::vector<SearchRegion > searchRegions_;
std::vector<SearchRegion > sbRegions_;
bool searchRegionsSet_=false;
void setSearchRegions() {
  if (searchRegionsSet_) return;
  
  //nb: some of the code *depends* on the fact that for there are equal numbers of corresponding
  //sbRegions and searchRegions, with the only difference being the MET selection!

  //also, for style reasons the 'owenId' should not contain the number of b tags.
  //everywhere that we use the owenId as an identifier, we combine with the number of b tags

  //oct25
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

  
  /*
  //2011 Summer result
  sbRegions_.push_back( SearchRegion( "ge1b","HT>=350","MET>=150&&MET<200","Loose",false)); //loose SB
  searchRegions_.push_back( SearchRegion( "ge1b","HT>=350","MET>=200","Loose")); //loose Sig

  sbRegions_.push_back( SearchRegion( "ge1b","HT>=500","MET>=150&&MET<200","Tight",false)); //tight SB
  searchRegions_.push_back( SearchRegion( "ge1b","HT>=500","MET>=300","Tight")); //tight Sig

  sbRegions_.push_back( SearchRegion( "ge2b","HT>=350","MET>=150&&MET<200","Loose",false)); //loose SB
  searchRegions_.push_back( SearchRegion( "ge2b","HT>=350","MET>=200","Loose")); //loose Sig

  sbRegions_.push_back( SearchRegion( "ge2b","HT>=500","MET>=150&&MET<200","Tight",false)); //tight SB
  searchRegions_.push_back( SearchRegion( "ge2b","HT>=500","MET>=300","Tight")); //tight Sig
  */

  searchRegionsSet_=true;
}

//very simple data container used by the SignalEffData class
class SystInfo {
public:
  SystInfo(double p=0, double m=0, int s=0);
  ~SystInfo();
  SystInfo( ifstream* input);

  double plus;
  double minus;
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

SystInfo::SystInfo(double p, double m, int s) :
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

  double rawYield; //yield in lumiScale_ invpb, or for SMS really the raw yield

  double effCorr; //factor including all corrections to the efficiency
  double totalSystematic(); 
  double symmetrize(const TString & which);

  double value(const TString &which) {return symmetrize(which);}

  void set(const TString & which, double valminus, double valplus);

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

  ofstream output(filename.Data());
  output<<rawYield<<endl<<effCorr<<endl;

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
  effCorr(1)
{
  //filename constructed from id
  TString filename = "SignalEffData.";
  filename += idtoload;

  ifstream input(filename.Data()); //should check that this is good
  input>>rawYield;
  input>>effCorr;

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
  effCorr(1)
{ 

  systematics["JES"] = SystInfo();             
  systematics["btag"] = SystInfo();            
  systematics["PDF"] = SystInfo();             
  systematics["MET"] = SystInfo();             
  systematics["PU"] = SystInfo();              
  systematics["JER"] = SystInfo();             
  systematics["kFactor"] = SystInfo();         
  systematics["cleaning"] = SystInfo(1e-2,1e-2,1);  
  systematics["LepVeto"] = SystInfo(2e-2,2e-2,1);   
  systematics["trigger"] = SystInfo(2.5e-2,2.5e-2,1);
  systematics["lumi"] = SystInfo(4.5e-2 , 4.5e-2,1);  
  //we don't really want this one this time.
  //keep it for now just for the sake of comparison with old results
  //  systematics["L2L3"] = SystInfo(1e-2 , 1e-2,1);

  //list of signal systematics:
  //  JES
  // btag efficiency
  // PDFs (acceptance)
  // [not done] PDFs (cross-section) -- not considered last time and not relevant for SMS
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
  systematics["PU"].plus = 4e-2;
  systematics["PU"].minus = -systematics["PU"].plus;
  systematics["PU"].status = 1;

  systematics["JER"].plus = 2e-2;
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

  if (which == "BTagEff02") return "btag";
  else if (which=="JERbias") return "JER";
  else if (which=="JES0") return "JES";
  else if (which=="METunc0") return "MET";
  else if (which =="PUunc0") return "PU";

  return which;

}

void SignalEffData::set(const TString & which, double valminus, double valplus) {

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

double SignalEffData::symmetrize(const TString & which) {

  double s=-1;

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
    double var1= fabs( it->second.plus);
    double var2= fabs( it->second.minus);
    s = var2>var1? var2:var1;
  }
  
  return s;
}

double SignalEffData::totalSystematic() {

  double total2=0;
  for (map<TString, SystInfo >::iterator isyst=systematics.begin(); isyst!=systematics.end(); ++isyst) {
    if ( isyst->second.status == 0) { 
      cout<<"WARNING -- systematic is unset! "<<isyst->first<<endl;
    }
    total2 += pow( symmetrize(isyst->first),2);
  }

  return 100*sqrt(total2);
}

struct OwenData {
  double Nsig; //number in signal region , data //done
  double Nsb; // number in SB, data             //done
  double Nsig_sl; //number in SL SIG, data  //done
  double Nsb_sl; // number in SL SB, data   //done
  double Nsig_ldp; //number in SIG, fail DP //done
  double Nsb_ldp;   // number in SB, fail DP //done

  //owen didn't ask for these
  //  double  Nlsb ;
  //  double  Nlsb_ldp;
  double  Nlsb_0b ;      // done
  double  Nlsb_0b_ldp;   // done

  double Nttbarmc_sig_ldp; //done
  double Nttbarmc_sb_ldp; //done

  double lsf_WJmc; //done
  double NWJmc_sig_ldp; //done
  double NWJmc_sb_ldp; //done

  double lsf_Znnmc; //done
  double NZnnmc_sig_ldp; //done
  double NZnnmc_sb_ldp; //done

  double lsf_Zjmc;
  double NZjmc_sig_ldp;
  double NZjmc_sb_ldp;

  double Nsingletopmc_sig_ldp;
  double Nsingletopmc_sb_ldp;

  //don't need DataLumi...that's just lumiScale_
} ;


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
  cout<<"Nsb_sl            "<<  owenMap_[owenKey].Nsb_sl<<endl;

  cout<<"Nsig_ldp          "<<  owenMap_[owenKey].Nsig_ldp<<endl;
  cout<<"Nsb_ldp           "<<  owenMap_[owenKey].Nsb_ldp<<endl;

  cout<<"Nlsb_0b           "<<  owenMap_[owenKey].Nlsb_0b<<endl;
  cout<<"Nlsb_0b_ldp       "<<  owenMap_[owenKey].Nlsb_0b_ldp<<endl;

  cout<<"Nttbarmc_sig_ldp  "<<  owenMap_[owenKey].Nttbarmc_sig_ldp<<endl;
  cout<<"Nttbarmc_sb_ldp   "<<  owenMap_[owenKey].Nttbarmc_sb_ldp<<endl;

  cout<<"Nsingletopmc_sig_ldp  "<<  owenMap_[owenKey].Nsingletopmc_sig_ldp<<endl;
  cout<<"Nsingletopmc_sb_ldp   "<<  owenMap_[owenKey].Nsingletopmc_sb_ldp<<endl;

  cout<<"lsf_WJmc          "<<  owenMap_[owenKey].lsf_WJmc<<endl;
  cout<<"NWJmc_sig_ldp     "<<  owenMap_[owenKey].NWJmc_sig_ldp<<endl;
  cout<<"NWJmc_sb_ldp      "<<  owenMap_[owenKey].NWJmc_sb_ldp<<endl;

  cout<<"lsf_Znnmc         "<<  owenMap_[owenKey].lsf_Znnmc<<endl;
  cout<<"NZnnmc_sig_ldp    "<<  owenMap_[owenKey].NZnnmc_sig_ldp<<endl;
  cout<<"NZnnmc_sb_ldp     "<<  owenMap_[owenKey].NZnnmc_sb_ldp<<endl;

  cout<<"lsf_Zjmc         "<<  owenMap_[owenKey].lsf_Zjmc<<endl;
  cout<<"NZjmc_sig_ldp    "<<  owenMap_[owenKey].NZjmc_sig_ldp<<endl;
  cout<<"NZjmc_sb_ldp     "<<  owenMap_[owenKey].NZjmc_sb_ldp<<endl;

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
bool useHLTeff_ = false;
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
TH1D* ratio=0; float ratioMin=0; float ratioMax=2;
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

map<pair<int,int>, TH1D* >  scanProcessTotalsMap;
void loadSusyScanHistograms() {
  if (loadedSusyHistos_) return;

  TFile* susyfile = 0;
  for (unsigned int isample=0; isample<samples_.size(); isample++) {
    if ( samples_[isample].Contains("mSUGRA")) {
      susyfile =  files_[currentConfig_][samples_[isample]];
      cout<<" Loading SUSY scan histograms from file: "<<susyfile->GetName()<<endl;
    }
  }

  if (susyfile==0) {cout<<"did not find mSugra in loadSusyScanHistograms!"<<endl; return;}

  for (int i=0; i<susyfile->GetListOfKeys()->GetEntries(); i++) {
    TString objname=  susyfile->GetListOfKeys()->At(i)->GetName();
    if (objname.BeginsWith("scanProcessTotals") ) {
      TString m0str=   objname.Tokenize("_")->At(1)->GetName();
      TString m12str=   objname.Tokenize("_")->At(2)->GetName();
      int m0 = m0str.Atoi();
      int m12 = m12str.Atoi();
      scanProcessTotalsMap[make_pair(m0,m12)] = (TH1D*) susyfile->Get(objname);
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

void drawPlotHeader() {
  //  return;
  float ypos = 0.97;
  if(doRatio_) ypos=ypos+0.012;
  // i'm gonna leave this out for now
  if (text1 != 0 ) delete text1;
  text1 = new TLatex(3.570061,23.08044,"CMS"); //no more preliminary!
  text1->SetNDC();
  text1->SetTextAlign(13);
  text1->SetX(0.68 + 0.2); //this 0.2 is because we got rid of the "Preliminary"
  text1->SetY(ypos+0.007);
  text1->SetTextFont(42);
  text1->SetTextSizePixels(24);
  text1->SetTextSize(0.045); //copied from ben's code. maybe needs a switch so it is only used for AN2011()
  text1->Draw();

  if (normalized_ == false) {
    TString astring;
    //astring.Form("%.0f pb^{-1} at #sqrt{s} = 7 TeV",lumiScale_);
    //astring.Form("%.1f fb^{-1} at #sqrt{s} = 7 TeV",lumiScale_/1000.);
    astring.Form("L_{int} = %.1f fb^{-1}, #sqrt{s} = 7 TeV",lumiScale_/1000.);
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
int ratiopadHeight = 250;
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
    thecanvas->Divide(1,2);
    const float padding=0.01; const float ydivide=0.2;
    thecanvas->GetPad(1)->SetPad( padding, ydivide + padding, 1-padding, 1-padding);
    thecanvas->GetPad(2)->SetPad( padding, padding, 1-padding, ydivide-padding);
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
TString getCutString(sampleType type, TString extraWeight="", TString thisSelection="", TString extraSelection="", int pdfWeightIndex=0,TString pdfSet="CTEQ", int susySubProcess=-1) {

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
  else if (type==kMC || type==kmSugraPoint || type==kmSugraPlane ) lumiscale=lumiScale_;
  else if (type==kSMSPlane||type==kSMSPoint) lumiscale=1;
  else {assert(0);}

  TString weightedcut="weight"; 
  
  weightedcut += "*(";
  weightedcut +=lumiscale;
  weightedcut+=")";
  
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
  if (useHLTeff_ &&  type!=kData) {
    weightedcut +="*hltHTeff";
    //    weightedcut +="*hltMHTeff"; //needs to be there now, but not compatible with Summer reducedTrees
  }

  if (btagSFweight_=="") btagSFweight_="1";
  if ( type==kData) {
    if (btagSFweight_=="probge1") weightedcut += "*(nbjets>=1)";
    else if (btagSFweight_=="probge2") weightedcut += "*(nbjets>=2)";
    else if (btagSFweight_=="probge3") weightedcut += "*(nbjets>=3)";
    else if (btagSFweight_=="prob1") weightedcut += "*(nbjets==1)";
    else if (btagSFweight_=="prob0") weightedcut += "*(nbjets==0)";
    else if (btagSFweight_=="1") {} //do nothing
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
    TH1D* thishist = scanProcessTotalsMap[make_pair(m0_,m12_)];
    if (thishist==0) cout<<"We've got a problem in getCutString!"<<endl;
    TString susyprocessweight = "(";
    int lowbound=0; int highbound=10;
    if (susySubProcess>=0) { lowbound=susySubProcess; highbound=susySubProcess;}
    for (int i=lowbound; i<=highbound; i++ ) {
      char thisweight[50];
      int thisn = TMath::Nint(thishist->GetBinContent(i));
      if (thisn==0) thisn=1; //avoid div by 0. if there are no events anyway then any value is ok
      sprintf(thisweight, "((SUSY_process==%d)*scanCrossSection%s/%d)",i,susyCrossSectionVariation_.Data(),thisn);
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
    //need to weight with the standard sigma / N
    //sigma should be stored in the scanCrossSection field. N needs to be pulled out of the histogram
    //not too important so implement later

    //forget cross section and just return without weight

    //    cout<<"not implemented yet!"<<endl;
    //    assert(0);
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


void addOverflowBin(TH1F* theHist) {
  //this code was written for when there was a customizable plot range (post-histo creation)
  //it could be made a lot simpler now
  //this one is copied from the function of the same name that takes a TH1D

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
    sampleColor_["T2bb"] =kGray;
    sampleColor_["T2tt"] =kGray;
    sampleColor_["QCD"] = kYellow;
    sampleColor_["PythiaQCD"] = kYellow;
    sampleColor_["PythiaPUQCD"] = kYellow;
    sampleColor_["PythiaPUQCDFlat"] = kYellow;
    sampleColor_["TTbarJets"]=kRed+1;
    sampleColor_["SingleTop"] = kMagenta;
    sampleColor_["WJets"] = kGreen-3;
    sampleColor_["WJetsZ2"] = kGreen-3;
    sampleColor_["ZJets"] = kAzure-2;
    sampleColor_["Zinvisible"] = kOrange-3;
    sampleColor_["SingleTop-sChannel"] = kMagenta+1; //for special cases
    sampleColor_["SingleTop-tChannel"] = kMagenta+2; //for special cases
    sampleColor_["SingleTop-tWChannel"] = kMagenta+3; //for special cases
    sampleColor_["SingleTopBar-sChannel"] = kMagenta+4; //for special cases
    sampleColor_["SingleTopBar-tChannel"] = kMagenta+5; //for special cases
    sampleColor_["SingleTopBar-tWChannel"] = kMagenta+6; //for special cases
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
    sampleColor_["T2bb"] =kCyan+2;
    sampleColor_["T2tt"] =kCyan+2;
    sampleColor_["QCD"] = 2;
    sampleColor_["PythiaQCD"] = 2;
    sampleColor_["PythiaPUQCD"] =2;
    sampleColor_["PythiaPUQCDFlat"] =2;
    sampleColor_["TTbarJets"]=4;
    sampleColor_["SingleTop"] = kMagenta;
    sampleColor_["WJets"] = kOrange;
    sampleColor_["WJetsZ2"] = kOrange;
    sampleColor_["ZJets"] = 7;
    sampleColor_["Zinvisible"] = kOrange+7;
    sampleColor_["SingleTop-sChannel"] = kMagenta+1; //for special cases
    sampleColor_["SingleTop-tChannel"] = kMagenta+2; //for special cases
    sampleColor_["SingleTop-tWChannel"] = kMagenta+3; //for special cases
    sampleColor_["TotalSM"] =kGreen+2; //owen requested 3
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

  samples_.push_back("TTbarJets");
  //flip this bool to control whether SingleTop is loaded as one piece or 3
  if (joinSingleTop) samples_.push_back("SingleTop");
  else {
    samples_.push_back("SingleTop-sChannel");
    samples_.push_back("SingleTop-tChannel");
    samples_.push_back("SingleTop-tWChannel");
  }
  samples_.push_back("WJets");
  samples_.push_back("ZJets");
  samples_.push_back("VV");
  samples_.push_back("Zinvisible");
  //samples_.push_back("HerwigQCDFlat");
  samples_.push_back("LM9");

}


void loadSamples(bool joinSingleTop=true) {
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
  samplesAll_.insert("WJets");
  samplesAll_.insert("ZJets");
  samplesAll_.insert("Zinvisible");
  samplesAll_.insert("SingleTop");
  samplesAll_.insert("SingleTop-sChannel");
  samplesAll_.insert("SingleTop-tChannel");
  samplesAll_.insert("SingleTop-tWChannel");
  samplesAll_.insert("SingleTopBar-sChannel");
  samplesAll_.insert("SingleTopBar-tChannel");
  samplesAll_.insert("SingleTopBar-tWChannel");
  samplesAll_.insert("HerwigQCDFlat");
  //  samplesAll_.insert("WJetsZ2");
  //  samplesAll_.insert("ZJetsZ2");
  samplesAll_.insert("VV");


  samplesAll_.insert("LM13");
  samplesAll_.insert("LM9");

  samplesAll_.insert("mSUGRAtanb40");
  samplesAll_.insert("T1bbbb");
  samplesAll_.insert("T2bb");
  samplesAll_.insert("T2tt");

 
  
  //FOR PLOTS
  ////////////

  configDescriptions_.setDefault("CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0");
  configDescriptions_.setCorrected("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0");

  
  //JES
  //configDescriptions_.addVariation("CSVM_PF2PATjets_JESdown_JERbias_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0",
  //				   "CSVM_PF2PATjets_JESup_JERbias_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0");
  //JER
  //  configDescriptions_.addVariation("CSVM_PF2PATjets_JES0_JERdown_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0",
  //"SSVHPT_PF2PATjets_JES0_JERup_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0");

  //unclustered MET
  //configDescriptions_.addVariation("CSVM_PF2PATjets_JES0_JERbias_PFMET_METuncDown_PUunc0_BTagEff0_HLTEff0", 
  //				   "CSVM_PF2PATjets_JES0_JERbias_PFMET_METuncUp_PUunc0_BTagEff0_HLTEff0");

  //PU
  //    configDescriptions_.addVariation("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUuncDown_BTagEff0_HLTEff0",
  //"CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUuncUp_BTagEff0_HLTEff0");

  //btag eff
  //configDescriptions_.addVariation("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEffdown_HLTEff0",
  //				   "CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEffup_HLTEff0");

  //HLT eff
  //    configDescriptions_.addVariation("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff0_HLTEffdown",
  //"CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff0_HLTEffup");

  ///////////////
  //////////////
 



/*
  //new btag eff prescription
  configDescriptions_.push_back("CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff02_HLTEff0");
  configDescriptions_.push_back("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff02_HLTEff0");
  //comment here to save time
  
  //JES
  configDescriptions_.push_back("CSVM_PF2PATjets_JESdown_JERbias_PFMET_METunc0_PUunc0_BTagEff02_HLTEff0");
  configDescriptions_.push_back("CSVM_PF2PATjets_JESup_JERbias_PFMET_METunc0_PUunc0_BTagEff02_HLTEff0");
  //JER
    configDescriptions_.push_back("CSVM_PF2PATjets_JES0_JERdown_PFMET_METunc0_PUunc0_BTagEff02_HLTEff0");
    configDescriptions_.push_back("CSVM_PF2PATjets_JES0_JERup_PFMET_METunc0_PUunc0_BTagEff02_HLTEff0");

  //unclustered MET
  configDescriptions_.push_back("CSVM_PF2PATjets_JES0_JERbias_PFMET_METuncDown_PUunc0_BTagEff02_HLTEff0");
  configDescriptions_.push_back("CSVM_PF2PATjets_JES0_JERbias_PFMET_METuncUp_PUunc0_BTagEff02_HLTEff0");

  //PU
      configDescriptions_.push_back("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUuncDown_BTagEff02_HLTEff0");
      configDescriptions_.push_back("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUuncUp_BTagEff02_HLTEff0");

  //btag eff
  configDescriptions_.push_back("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEffdown2_HLTEff0");
  configDescriptions_.push_back("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEffup2_HLTEff0");

  //HLT eff
  //    configDescriptions_.push_back("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff02_HLTEffdown");
  //    configDescriptions_.push_back("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff02_HLTEffup");
*/

/*  
  //new btag eff prescription
 //for summer11 signal systematics!
  configDescriptions_.setDefault("SSVHPT_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff02_HLTEff0");
  configDescriptions_.setCorrected("SSVHPT_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff02_HLTEff0");
  //comment here to save time
  
  //JES
  configDescriptions_.addVariation("SSVHPT_PF2PATjets_JESdown_JERbias_PFMET_METunc0_PUunc0_BTagEff02_HLTEff0",
				   "SSVHPT_PF2PATjets_JESup_JERbias_PFMET_METunc0_PUunc0_BTagEff02_HLTEff0");
  //JER - out for scans
  //  configDescriptions_.addVariation("SSVHPT_PF2PATjets_JES0_JERdown_PFMET_METunc0_PUunc0_BTagEff02_HLTEff0",
				    //  				   "SSVHPT_PF2PATjets_JES0_JERup_PFMET_METunc0_PUunc0_BTagEff02_HLTEff0");

  //unclustered MET
  configDescriptions_.addVariation("SSVHPT_PF2PATjets_JES0_JERbias_PFMET_METuncDown_PUunc0_BTagEff02_HLTEff0",
				   "SSVHPT_PF2PATjets_JES0_JERbias_PFMET_METuncUp_PUunc0_BTagEff02_HLTEff0");

  //PU - out for scans
  //   configDescriptions_.addVariation("SSVHPT_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUuncDown_BTagEff02_HLTEff0",
  //  				   "SSVHPT_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUuncUp_BTagEff02_HLTEff0");

  //btag eff
  configDescriptions_.addVariation("SSVHPT_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEffdown2_HLTEff0",
				   "SSVHPT_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEffup2_HLTEff0");

  //HLT eff
  //    configDescriptions_.push_back("SSVHPT_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff02_HLTEffdown");
  //    configDescriptions_.push_back("SSVHPT_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff02_HLTEffup");
  
*/

  currentConfig_=configDescriptions_.getDefault();

  //these blocks are just a "dictionary"
  //no need to ever comment these out
  setColorScheme("stack");

  sampleLabel_["mSUGRAtanb40"] = "tan #beta = 40";
  sampleLabel_["T1bbbb"] = "T1bbbb";
  sampleLabel_["T2bb"] = "T2bb";
  sampleLabel_["T2tt"] = "T2tt";
  sampleLabel_["LM13"] = "LM13";
  sampleLabel_["LM9"] = "LM9";
  sampleLabel_["QCD"] = "QCD (madgraph)";
  sampleLabel_["PythiaQCD"] = "QCD (no PU)";
  sampleLabel_["PythiaPUQCDFlat"] = "QCD"; 
  sampleLabel_["PythiaPUQCD"] = "QCD";
  sampleLabel_["TTbarJets"]="t#bar{t}";
  sampleLabel_["SingleTop"] = "Single-Top";
  sampleLabel_["WJets"] = "W#rightarrowl#nu";
  sampleLabel_["WJetsZ2"] = "W#rightarrowl#nu (Z2)";
  sampleLabel_["ZJets"] = "Z/#gamma*#rightarrowl^{+}l^{-}";
  sampleLabel_["Zinvisible"] = "Z#rightarrow#nu#nu";
  sampleLabel_["SingleTop-sChannel"] = "Single-Top (s)";
  sampleLabel_["SingleTop-tChannel"] = "Single-Top (t)";
  sampleLabel_["SingleTop-tWChannel"] = "Single-Top (tW)";
  sampleLabel_["SingleTopBar-sChannel"] = "Single-TopBar (s)";
  sampleLabel_["SingleTopBar-tChannel"] = "Single-TopBar (t)";
  sampleLabel_["SingleTopBar-tWChannel"] = "Single-TopBar (tW)";
  sampleLabel_["VV"] = "Diboson";
  sampleLabel_["HerwigQCDFlat"] = "Herwig QCD";
  sampleLabel_["TotalSM"] = "SM";
  sampleLabel_["Total"] = "SM + LM13"; //again, this is a hack

  sampleMarkerStyle_["mSUGRAtanb40"] = kFullStar;
  sampleMarkerStyle_["T1bbbb"] = kFullStar;
  sampleMarkerStyle_["T2bb"] = kFullStar;
  sampleMarkerStyle_["T2tt"] = kFullStar;
  sampleMarkerStyle_["LM13"] = kFullStar;
  sampleMarkerStyle_["LM9"] = kFullStar;
  sampleMarkerStyle_["QCD"] = kFullCircle;
  sampleMarkerStyle_["PythiaQCD"] = kOpenCircle;
  sampleMarkerStyle_["PythiaPUQCDFlat"] = kOpenCircle;  
  sampleMarkerStyle_["PythiaPUQCD"] = kOpenCircle;
  sampleMarkerStyle_["TTbarJets"]= kFullSquare;
  sampleMarkerStyle_["SingleTop"] = kOpenSquare;
  sampleMarkerStyle_["WJets"] = kMultiply;
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
  sampleMarkerStyle_["TotalSM"] = kOpenCross; //FIXME?
  sampleMarkerStyle_["Total"] = kDot; //FIXME?

  sampleOwenName_["mSUGRAtanb40"] = "msugra40";
  sampleOwenName_["T1bbbb"] = "t1bbbb";
  sampleOwenName_["T2bb"] = "t2bb";
  sampleOwenName_["T2tt"] = "t2tt";
  sampleOwenName_["LM13"] = "lm13";
  sampleOwenName_["LM9"] = "lm9";
  sampleOwenName_["QCD"] = "qcd";
  sampleOwenName_["PythiaQCD"] = "qcd";
  sampleOwenName_["PythiaPUQCDFlat"] = "qcd"; 
  sampleOwenName_["PythiaPUQCD"] = "qcd";
  sampleOwenName_["TTbarJets"]="ttbar";
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
  sampleOwenName_["TotalSM"] = "totalsm";
  sampleOwenName_["Total"] = "total";  

  //  for (std::vector<TString>::iterator iconfig=configDescriptions_.begin(); iconfig!=configDescriptions_.end(); ++iconfig) {
  for (unsigned int iconfig=0; iconfig<configDescriptions_.size(); ++iconfig) {
    for (std::set<TString>::iterator isample=samplesAll_.begin(); isample!=samplesAll_.end(); ++isample) {
      TString fname="reducedTree.";
      TString thisconfig = configDescriptions_.at(iconfig);
      fname += thisconfig;
      fname+=".";
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
  else if (sample.Contains("T2bb") && planeOrPoint=="point") return kSMSPoint;
  else if (sample.Contains("T2bb") && planeOrPoint=="plane") return kSMSPlane;
  else if (sample.Contains("T2tt") && planeOrPoint=="point") return kSMSPoint;
  else if (sample.Contains("T2tt") && planeOrPoint=="plane") return kSMSPlane;

  return kMC;

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
	    const TString xtitle, TString ytitle, TString filename="") {

  //for now hard-code this to plot only the first sample and only COLZ!

  loadSamples();
  if (filename=="") filename=var;
  gROOT->SetStyle("CMS");

  renewCanvas();
  if (h2d!=0) delete h2d;
  const TString hname="h2d";
  h2d = new TH2D(hname,"",nbins,low,high,nbinsy,lowy,highy);
  h2d->Sumw2();
  TString opt="colz";
  for (unsigned int isample=0; isample<1 ; isample++) { //plot only the first sample!
    gROOT->cd();
    TTree* tree = (TTree*) files_[currentConfig_][samples_[isample]]->Get("reducedTree");
    gROOT->cd();
    TString weightopt= useFlavorHistoryWeights_ && samples_[isample].Contains("WJets") ? "flavorHistoryWeight" : "";
    TString drawstring=vary;
    drawstring+=":";
    drawstring+=var;
    tree->Project(hname,drawstring,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,"",0,""));

    //now the histo is filled
    
    h2d->SetXTitle(xtitle);
    h2d->SetYTitle(ytitle);
  }
  h2d->Draw(opt);

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

}

bool isSampleSM(const TString & name) {

  if (name.Contains("LM")) return false;

  if (name.Contains("SUGRA")) return false;

  return true;
}


void drawPlots(const TString var, const int nbins, const float low, const float high, const TString xtitle, TString ytitle, TString filename="", const float* varbins=0) {
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
    tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),weightopt,selection_,"",0,"").Data());
    
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
    hdata->Sumw2();
    gROOT->cd();
    dtree->Project(hname,var,getCutString(true));
    //now the histo is filled
    
    hdata->UseCurrentStyle(); //maybe not needed anymore
    hdata->SetMarkerColor(kBlack);
    hdata->SetLineWidth(2);
    hdata->SetMarkerStyle(kFullCircle);
    hdata->SetMarkerSize(1);
    if (renormalizeBins_) renormBins(hdata,2 ); //manipulates the histogram //FIXME hard-coded "2"
    if (addOverflow_)     addOverflowBin(hdata); // manipulates the histogram!
    leg->AddEntry(hdata,"Data");

    if (!quiet_)    cout<<"Data underflow: " <<hdata->GetBinContent(0)<<endl;//BEN
    hdata->Draw("SAME");
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
    }
    if (doRatio_) {
      thecanvas->cd(2);
      ratio->Divide(hdata,totalsm);
      ratio->SetMinimum(ratioMin);
      ratio->SetMaximum(ratioMax);
      ratio->GetYaxis()->SetNdivisions(200 + int(ratioMax-ratioMin)+1);    //set ticks ; to be seen if this really works
      ratio->GetYaxis()->SetLabelSize(0.2); //make y label bigger
      ratio->Draw();
      thecanvas->GetPad(2)->SetTopMargin(0.1);
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
void drawR(const TString vary, const float cutVal, const TString var, const int nbins, const float low, const float high, const TString& savename, const float* varbins=0) {
  const TString ytitle="N pass / N fail";
  TString xtitle = var;
  if(var=="MET") xtitle = "E_{T}^{miss} [GeV]"; 
  bool dataOnly = true;

  //const TString var = "MET"; //hardcoded for now
  cout << "x axis: " << var << endl;
  TString cstring1 = vary, cstring2=vary;
  cstring1 += " >= ";
  cstring2 += " < ";
  cstring1 += cutVal;
  cstring2 += cutVal;

  //terrible hack to decide if bias correction should be calculated
  bool calcBiasCorr = false;
  if(nbins==4 && low>-0.01 && low<0.01 && high>199.99 && high<200.01) calcBiasCorr=true;
  float cb_qcd=0, cb_qcd_err=0, cb_sm=0, cb_sm_err=0, cb_data=0, cb_data_err=0;
  float cp_qcd=0, cp_qcd_err=0, cp_sm=0, cp_sm_err=0, cp_data=0, cp_data_err=0;
  float n_qcd_sb = 0, n_qcd_sb_err = 0, n_qcd_sig = 0, n_qcd_sig_err = 0;
  float n_qcd_a = 0, n_qcd_a_err = 0, n_qcd_d = 0, n_qcd_d_err = 0;

  loadSamples();

  gROOT->SetStyle("CMS");

  TString opt=doRatio_? "ratio":"";
  renewCanvas(opt);

  //in first incarnation, make a separate r(MET) plot for each sample in the list

//   TH1D* qcdPass = new TH1D("qcdPass","",nbins,low,high);
//   TH1D* qcdFail = new TH1D("qcdFail","",nbins,low,high);
//   TH1D* qcdRatio = new TH1D("qcdRatio","",nbins,low,high);

//   qcdPass->Sumw2();
//   qcdFail->Sumw2();
//   qcdRatio->Sumw2();

  resetHistos(); //delete existing histograms

  renewLegend();


  // === begin correlation hack ====
  gROOT->cd();
  TH2D totalsm2d_50("totalsm2d_50","",50,50,100,50,0,TMath::Pi());
  TH2D totalsm2d_SB("totalsm2d_SB","",50,100,150,50,0,TMath::Pi());
  totalsm2d_50.Sumw2();
  totalsm2d_SB.Sumw2();
  TH2D data2d_50("data2d_50","",50,50,100,50,0,TMath::Pi());
  TH2D data2d_SB("data2d_SB","",50,100,150,50,0,TMath::Pi());
  data2d_50.Sumw2();
  data2d_SB.Sumw2();
  // === end correlation hack ===

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

      //comment out filling of these for now to save time
      TH2D this2d_SB("this2d_SB","",50,100,150,50,0,TMath::Pi());
      //      tree->Project("this2d_SB","minDeltaPhi:MET",getCutString().Data());
      TH2D this2d_50("this2d_50","",50,50,100,50,0,TMath::Pi());
      //      tree->Project("this2d_50","minDeltaPhi:MET",getCutString().Data());

      totalsm2d_SB.Add(&this2d_SB);
      totalsm2d_50.Add(&this2d_50);
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

    if(dodata_) drawPlotHeader();

    if (hnameR.Contains("QCD") && !dataOnly) { //HACK draw only qcd
      histos_[hnameR]->Draw(drawopt);
      if (!drawopt.Contains("same")) drawopt+=" same";
      
      //hack to save
      //TH1D* hQCD = (TH1D*)histos_[hnameR]->Clone(savename+"_qcd");
      //hQCD->SaveAs(savename+"_qcd.root");
     
      //drawPlotHeaderInside(); //for qcd only plots, this works.
      
      if (firsthist="") firsthist = hnameR;
      if (histos_[hnameR]->GetMaximum() > max) max = histos_[hnameR]->GetMaximum();
      leg->AddEntry(histos_[hnameR], sampleLabel_[samples_[isample]]);
      
      if(calcBiasCorr){
	cp_qcd = histos_[hnameR]->GetBinContent(3)/ histos_[hnameR]->GetBinContent(2);
	cp_qcd_err = jmt::errAoverB( histos_[hnameR]->GetBinContent(3), histos_[hnameR]->GetBinError(3), histos_[hnameR]->GetBinContent(2), histos_[hnameR]->GetBinError(2)); 
	cb_qcd = histos_[hnameR]->GetBinContent(4)/ histos_[hnameR]->GetBinContent(3);
	cb_qcd_err = jmt::errAoverB( histos_[hnameR]->GetBinContent(4), histos_[hnameR]->GetBinError(4), histos_[hnameR]->GetBinContent(3), histos_[hnameR]->GetBinError(3)); 
	n_qcd_sb = histos_[hnameP]->GetBinContent(3);
	n_qcd_sb_err = histos_[hnameP]->GetBinError(3);
	n_qcd_sig = histos_[hnameP]->GetBinContent(4);
	n_qcd_sig_err = histos_[hnameP]->GetBinError(4);
     	n_qcd_a = histos_[hnameF]->GetBinContent(3);
	n_qcd_a_err = histos_[hnameF]->GetBinError(3);
	n_qcd_d = histos_[hnameF]->GetBinContent(4);
	n_qcd_d_err = histos_[hnameF]->GetBinError(4);
      }
    }
    
  }

  if(!dataOnly){
    histos_[firsthist]->SetMaximum( max*maxScaleFactor_);
    hinteractive =  histos_[firsthist];
  }

  totalsm->Divide(totalsm_pass,totalsm_fail);
  if (drawTotalSM_ && !dataOnly) {
    totalsm->Draw("hist e same");
    //    leg->Clear();
    leg->AddEntry(totalsm,sampleLabel_["TotalSM"]);
  }
  if(calcBiasCorr){
    cp_sm = totalsm->GetBinContent(3)/totalsm->GetBinContent(2);
    cp_sm_err = jmt::errAoverB(totalsm->GetBinContent(3),totalsm->GetBinError(3),totalsm->GetBinContent(2),totalsm->GetBinError(2));
    cb_sm = totalsm->GetBinContent(4)/totalsm->GetBinContent(3);
    cb_sm_err = jmt::errAoverB(totalsm->GetBinContent(4),totalsm->GetBinError(4),totalsm->GetBinContent(3),totalsm->GetBinError(3));
  }


  if (dodata_) {
    gROOT->cd();
    if (!quiet_)   cout<<"Drawing data!"<<endl;
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

    dtree->Project("data2d_SB","minDeltaPhi:MET",getCutString(kData,"",selection_).Data());
    dtree->Project("data2d_50","minDeltaPhi:MET",getCutString(kData,"",selection_).Data());

    //    hdata->UseCurrentStyle(); //maybe not needed anymore
    hdata->SetMarkerColor(kBlack);
    hdata->SetLineWidth(2);
    hdata->SetMarkerStyle(kFullCircle);
    hdata->SetMarkerSize(1);

    thecanvas->cd(1);
    hdata->SetYTitle(ytitle);
    hdata->SetXTitle(xtitle);
    hdata->GetYaxis()->SetLabelSize(0.04); //make y label bigger
    hdata->GetXaxis()->SetLabelSize(0.04); //make x label bigger
    if(dataOnly) hdata->Draw();
    else hdata->Draw("SAME");
    leg->AddEntry(hdata,"Data");

    if (hdata->GetMaximum() > max && !dataOnly)  {
      histos_[firsthist]->SetMaximum( maxScaleFactor_*hdata->GetMaximum());
      totalsm->SetMaximum(maxScaleFactor_*hdata->GetMaximum());
    }
    else if (doCustomPlotMax_ && dataOnly) {
      hdata->SetMaximum(customPlotMax_);
    }
    else if (doCustomPlotMax_) {
      histos_[firsthist]->SetMaximum( customPlotMax_);
      totalsm->SetMaximum(customPlotMax_);
    }
    if (doCustomPlotMin_ && dataOnly) {
      hdata->SetMinimum(customPlotMin_);
    }
    else if (doCustomPlotMin_) {
      histos_[firsthist]->SetMinimum( customPlotMin_);
      totalsm->SetMinimum(customPlotMin_);
    }
    
    //    cratio->cd();
    thecanvas->cd(2);
    gPad->SetRightMargin(0.1);
    gPad->Modified();
    if (ratio!=0) delete ratio;
    if(!dataOnly){
      ratio = (varbins==0) ? new TH1D("ratio","data/(SM MC)",nbins,low,high) : new TH1D("ratio","data/(SM MC)",nbins,varbins);
      ratio->Sumw2();
      ratio->Divide(hdata,totalsm); 
      ratio->SetMinimum(ratioMin);
      ratio->SetMaximum(ratioMax);
      ratio->GetYaxis()->SetNdivisions(200 + int(ratioMax-ratioMin)+1);    //set ticks ; to be seen if this really works
      ratio->GetYaxis()->SetLabelSize(0.1); //make y label bigger
      ratio->GetXaxis()->SetLabelSize(0.1); //make x label bigger
      ratio->Draw();
    }
    cout<<"KS Test results (shape only): "<<hdata->KolmogorovTest(totalsm)<<endl;;

    if(calcBiasCorr){
      cp_data = hdata->GetBinContent(3)/hdata->GetBinContent(2);
      cp_data_err = jmt::errAoverB(hdata->GetBinContent(3),hdata->GetBinError(3),hdata->GetBinContent(2),hdata->GetBinError(2)); 
      cb_data = hdata->GetBinContent(4)/hdata->GetBinContent(3);
      cb_data_err = jmt::errAoverB(hdata->GetBinContent(4),hdata->GetBinError(4),hdata->GetBinContent(3),hdata->GetBinError(3)); 
    }
  }
  else {
    if (doCustomPlotMax_) {
      histos_[firsthist]->SetMaximum( customPlotMax_);
      totalsm->SetMaximum(customPlotMax_);
    }
    if (doCustomPlotMin_) {
      histos_[firsthist]->SetMinimum( customPlotMin_);
      totalsm->SetMinimum(customPlotMin_);
    }
  }
  thecanvas->cd(1);
   if (doleg_)  leg->Draw();

  thecanvas->SaveAs("mindpPassOverFail-"+savename+".eps");
  thecanvas->SaveAs("mindpPassOverFail-"+savename+".pdf");
  thecanvas->SaveAs("mindpPassOverFail-"+savename+".png");

//   TCanvas* c2d=new TCanvas("c2d","2d",800,800);
//   c2d->Divide(2,2);
//   c2d->cd(1);
//   totalsm2d_50.DrawCopy("colz");
//   c2d->cd(2);
//   totalsm2d_SB.DrawCopy("colz");
//   c2d->cd(3);
//   data2d_50.DrawCopy("colz");
//   c2d->cd(4);
//   data2d_SB.DrawCopy("colz");
  cout<<"Total SM MC correlation [50<MET<100]  = "<<totalsm2d_50.GetCorrelationFactor()<<endl;
  cout<<"Total SM MC correlation [100<MET<150] = "<<totalsm2d_SB.GetCorrelationFactor()<<endl;
  cout<<"Data correlation [50<MET<100]         = "<<data2d_50.GetCorrelationFactor()<<endl;
  cout<<"Data correlation [100<MET<150]        = "<<data2d_SB.GetCorrelationFactor()<<endl;
  if(calcBiasCorr){
    cout<<endl;
    cout<<"Pseudo bias correction (using 50<MET<100 and 100<MET<150):" << endl;
    cout<<"QCD MC: "<<cp_qcd<<" +/- "<<cp_qcd_err<<endl;
    cout<<"SM MC: "<<cp_sm<<" +/- "<<cp_sm_err<<endl;
    cout<<"Data: "<<cp_data<<" +/- "<<cp_data_err<<endl;
    cout<<endl;
    cout<<"True bias correction (using 100<MET<150 and 150<MET):" << endl;
    cout<<"QCD MC: "<<cb_qcd<<" \\pm "<<cb_qcd_err<<endl;
    cout<<"SM MC: "<<cb_sm<<" \\pm "<<cb_sm_err<<endl;
    cout<<"Data: "<<cb_data<<" \\pm "<<cb_data_err<<endl;
    cout << endl;
    cout << "QCD Event Counts" << endl;
    cout << "SB: " << n_qcd_sb << " +- " << n_qcd_sb_err << endl;
    cout << "SIG: " << n_qcd_sig << " +- " << n_qcd_sig_err << endl;
    cout << "A: " << n_qcd_a << " +- " << n_qcd_a_err << endl;
    cout << "D: " << n_qcd_d << " +- " << n_qcd_d_err << endl;
    cout << endl;
  }
  cout<<"End of drawR()"<<endl;
}

void drawR(const TString vary, const float cutVal, const int nbins, const float low, const float high, const TString& savename) {//for backwards compatibility
  drawR(vary, cutVal, "MET", nbins, low, high, savename, 0); 
}

void drawR(const TString vary, const float cutVal, const TString var, const int nbins, const float* varbins=0, const TString& savename="") {//more natural
  drawR(vary, cutVal, var, nbins, 0, 1, savename, varbins);
}

//utility function for making output more readable
TString format_nevents(double n,double e) {

  //  const bool moreDigits = false;
  const bool moreDigits = true;
  const int eCutoff = moreDigits ? 10 : 1;
  const int extraDigits = moreDigits ? 1:0;

  TString mathmode = latexMode_ ? "$" : "";
  
  char out[100];
  if (e >= eCutoff || e < 0.00001) { //show whole numbers only
    sprintf(out,"%s%.0f%s%.0f%s",mathmode.Data(),n,pm.Data(),e,mathmode.Data());
  }
  else {
    int nfig = ceil(fabs(log10(e))) + extraDigits;
    TString form="%s%.";
    form+=nfig; form+="f%s%.";
    form+=nfig; form+="f%s";
    sprintf(out,form.Data(),mathmode.Data(),n,pm.Data(),e,mathmode.Data());
  }
  return TString(out);
}

typedef map<pair<int,int>, pair<double,double> > susyScanYields;
TH2D* scanSMSngen=0;
susyScanYields getSusyScanYields(const TString & sampleOfInterest,const int pdfindex=0, const TString & pdfset="") {
  cout<<" ~~ begin slow version of getSusyScanYields()"<<endl;

  /*
this function was written to draw mSugra. It can also draw T1bbbb with no problem.
The pdfindex and pdfset can also be specified, but it is not good for checking the whole mSugra plane because it is too slow
  */

  if (sampleOfInterest.Contains("mSUGRA") ) assert(pdfset=="" && pdfindex==0); //just until i have time to think about it

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
  int  nbinsx=210; float lowx=-0.5; float highx=2100-0.5;

  TString vary="m12"; TString ytitle=vary;
  int  nbinsy=110; float lowy=-0.5; float highy=1100-0.5;

  if (sampleOfInterest== "T1bbbb" || sampleOfInterest== "T2bb" || sampleOfInterest== "T2tt") { //change binning
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

  //[for mSugra]
  //at this point, each bin contains Npass_i * sigma_i for that (m0,m12)
  //need to divide by N_i for each (m0,m12)

  //[for SMS]
  //each bin contains just raw Npass_i for that (mgl,mLSP)

  if (sampleOfInterest.Contains("mSUGRA") ) {
    //loop over i and histo bins
    for (int i=0; i<=nsubprocesses; i++) {
      for (map<pair<int,int>, TH1D* >::iterator iscanpoint = scanProcessTotalsMap.begin(); iscanpoint!=scanProcessTotalsMap.end(); ++iscanpoint) {
	int m0=iscanpoint->first.first;
	int m12=iscanpoint->first.second;
	TH1D* thishist = scanProcessTotalsMap[make_pair(m0,m12)];
	int thisn = TMath::Nint(thishist->GetBinContent(i));
	int bin=  raw0[i]->FindBin(m0,m12);
	double N_i_thispoint = raw0[i]->GetBinContent(bin);
	double err_i_thispoint = raw0[i]->GetBinError(bin);
	if (thisn == 0) {
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
    
    //now we have Npass_i * sigma_i * lumi / Ngen_i 
    //all that is left is to make the sum over i
    
    for (map<pair<int,int>, TH1D* >::iterator iscanpoint = scanProcessTotalsMap.begin(); iscanpoint!=scanProcessTotalsMap.end(); ++iscanpoint) {
      double Nraw = 0, errraw=0;
      
      int bin=  raw0[0]->FindBin(iscanpoint->first.first , iscanpoint->first.second);
      for (unsigned int i=0; i<raw0.size(); i++) {
	Nraw += raw0[i]->GetBinContent(bin);
	errraw += pow(raw0[i]->GetBinError(bin),2);
      }
      //cout<<iscanpoint->first.first<<" "<<iscanpoint->first.second<<" "<<Nraw<< " +/- "<<sqrt(errraw)<<endl;
      theYields[iscanpoint->first] = make_pair(Nraw,sqrt(errraw));
    }
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
	if ( TMath::Nint(scanSMSngen->GetBinContent(scanSMSngen->FindBin(mgl,mlsp))) >= 1000 ) {
	  theYields[make_pair(mgl,mlsp)] = make_pair(nevents,0); //skip errors for now
	}
      }
    }
  }

  //try to clean up
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
  int  nbinsx=210; float lowx=-0.5; float highx=2100-0.5;

  TString vary="m12"; TString ytitle=vary;
  int  nbinsy=110; float lowy=-0.5; float highy=1100-0.5;
  if (sampleOfInterest== "T1bbbb" || sampleOfInterest== "T2bb" || sampleOfInterest== "T2tt") { //change binning
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
  map<int, vector<TH2D*> > raw0;
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
      raw0[isusy].push_back(new TH2D(hname,"raw event counts",nbinsx,lowx,highx,nbinsy,lowy,highy));
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
      //at this point, each bin contains Npass_i * sigma_i for that (m0,m12)
      //need to divide by N_i for each (m0,m12)
      for (int i=0; i<=nsubprocesses; i++) {
	for (map<pair<int,int>, TH1D* >::iterator iscanpoint = scanProcessTotalsMap.begin(); iscanpoint!=scanProcessTotalsMap.end(); ++iscanpoint) {
	  int m0=iscanpoint->first.first;
	  int m12=iscanpoint->first.second;
	  TH1D* thishist = scanProcessTotalsMap[make_pair(m0,m12)];
	  int thisn = TMath::Nint(thishist->GetBinContent(i));
	  int bin=  raw0[i][ipdf]->FindBin(m0,m12);
	  double N_i_thispoint = raw0[i][ipdf]->GetBinContent(bin);
	  double err_i_thispoint = raw0[i][ipdf]->GetBinError(bin);
	  if (thisn == 0) {
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
      for (map<pair<int,int>, TH1D* >::iterator iscanpoint = scanProcessTotalsMap.begin(); iscanpoint!=scanProcessTotalsMap.end(); ++iscanpoint) {
	double Nraw = 0, errraw=0;
	
	int bin=  raw0[0][ipdf]->FindBin(iscanpoint->first.first , iscanpoint->first.second);
	//since we're using a map, we could do an iterator loop here, but this works too
	for (unsigned int i=0; i<raw0.size(); i++) { //we want the size of the *map* (number of susy subprocesses)
	  Nraw += raw0[i][ipdf]->GetBinContent(bin);
	  errraw += pow(raw0[i][ipdf]->GetBinError(bin),2);
	}
	//cout<<iscanpoint->first.first<<" "<<iscanpoint->first.second<<" "<<Nraw<< " +/- "<<sqrt(errraw)<<endl;
	oneSetOfYields[iscanpoint->first] = make_pair(Nraw,sqrt(errraw));
      }
    }
    else { //for T1bbbb and hopefully other samples that don't have the susy subprocesses 
      for (int i=1; i<=nbinsx; i++) {
	for (int j=1; j<=nbinsy; j++) {
	  int mgl=TMath::Nint(raw0[0][ipdf]->GetXaxis()->GetBinLowEdge(i));
	  int mlsp=TMath::Nint(raw0[0][ipdf]->GetYaxis()->GetBinLowEdge(j));
	  double nevents = raw0[0][ipdf]->GetBinContent(raw0[0][ipdf]->FindBin(mgl,mlsp));
	  if ( TMath::Nint(scanSMSngen->GetBinContent(scanSMSngen->FindBin(mgl,mlsp))) >= 1000 ) {
	    oneSetOfYields[make_pair(mgl,mlsp)] = make_pair(nevents,0); //skip errors for now
	  }
	}
      }
    }

    theYields.push_back(oneSetOfYields);
  }
  timer.Stop();
  cout<<"[getSusyScanYields() (fast)] CPU time ("<<pdfset<<") = "<<timer.CpuTime()<<endl;

  //try to clean up
  //need to fix this for the new data structure
  for ( map<int, vector<TH2D*> >::iterator im = raw0.begin(); im!=raw0.end() ; ++im) {
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
	if ( XplusMax.count(ipoint->first) ==0) XplusMax[ipoint->first]=make_pair(0,0);
	if ( XminusMax.count(ipoint->first) ==0) XminusMax[ipoint->first]=make_pair(0,0);
	double val = cteqXplus[i][ipoint->first].first;
	XplusMax[ipoint->first].first += val; //sum up the values for pdf i and this scan point
	XminusMax[ipoint->first].first += val*val; //sum up the values for pdf i and this scan point
      }
      nn++;
    }
    for (susyScanYields::iterator ipoint=XplusMax.begin() ; ipoint!=XplusMax.end(); ++ipoint) {
      XplusMax[ipoint->first].first = XplusMax[ipoint->first].first / double(nn);
      XminusMax[ipoint->first].first = XminusMax[ipoint->first].first / double(nn);
    }

  }
  else {
    //loop over pdfs and over scan points
    for (unsigned int i=0; i<cteqXplus.size(); i++) {
      for (susyScanYields::iterator ipoint=cteqXplus[i].begin(); ipoint!=cteqXplus[i].end(); ++ipoint) {
	//	if (ipoint->first ==apoint) cout<<ipoint->first.first<<" "<<ipoint->first.second<<" -- "<<nominal[ipoint->first].first<<" " <<cteqXplus[i][ipoint->first].first<<" "<<cteqXminus[i][ipoint->first].first<<endl;
	double diff1 =  cteqXplus[i][ipoint->first].first - nominal[ipoint->first].first ;
	double diff2 =  cteqXminus[i][ipoint->first].first - nominal[ipoint->first].first ;
	double larger = diff1>diff2 ? diff1 : diff2;
	if ( 0 > larger ) larger = 0;
	if ( XplusMax.count(ipoint->first) ==0) {
	  XplusMax[ipoint->first]=make_pair(0,0);
	  //  if (ipoint->first ==apoint) cout<<" ***** CREATING XplusMax"<<endl;
	}
	XplusMax[ipoint->first].first += larger*larger;
	//	if (ipoint->first ==apoint) cout<<"       XplusMax = "<<XplusMax[ipoint->first].first<<endl;
      }
    }
    for (unsigned int i=0; i<cteqXplus.size(); i++) {
      for (susyScanYields::iterator ipoint=cteqXplus[i].begin(); ipoint!=cteqXplus[i].end(); ++ipoint) {
	double diff1 =  nominal[ipoint->first].first  - cteqXplus[i][ipoint->first].first ;
	double diff2 =  nominal[ipoint->first].first - cteqXminus[i][ipoint->first].first ;
	double larger = diff1>diff2 ? diff1 : diff2;
	if ( 0 > larger ) larger = 0;
	if ( XminusMax.count(ipoint->first) ==0) XminusMax[ipoint->first]=make_pair(0,0);
	XminusMax[ipoint->first].first += larger*larger;
      }
    }

    const double scale = (pdfset=="CTEQ") ? 1.645 : 1;
    for (susyScanYields::iterator ipoint=XplusMax.begin() ; ipoint!=XplusMax.end() ; ++ipoint) {
      ipoint->second.first = sqrt(ipoint->second.first);
      ipoint->second.first = ipoint->second.first / scale;
    }
    for (susyScanYields::iterator ipoint=XminusMax.begin() ; ipoint!=XminusMax.end() ; ++ipoint) {
      ipoint->second.first = sqrt(ipoint->second.first);
      ipoint->second.first = ipoint->second.first / scale;
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
  TString selectionGe1bLoose;
  TString selectionGe2bLoose;
  TString selectionGe3bLoose;
  TString selectionEq1bTight;
  TString selectionGe1bTight;
  TString selectionGe2bTight;  

  /*
  //inclusive
  thisSelection="";
  cut=getCutString(lumiScale_,thisSelection);
  vectorOfCuts.push_back(cut);
  stageCut.push_back("Inclusive");

  //trigger
  if (thisSelection=="") thisSelection += "cutTrigger==1";
  else thisSelection += " && cutTrigger";
  cut=getCutString(lumiScale_,thisSelection);
  vectorOfCuts.push_back(cut);
  stageCut.push_back("Trigger");

  //PV
  if (thisSelection=="") thisSelection += "cutPV==1";
  else thisSelection += " && cutPV==1";
  cut=getCutString(lumiScale_,thisSelection);
  vectorOfCuts.push_back(cut);
  stageCut.push_back("PV");
  */

  //HT
  if (thisSelection=="") {thisSelection += "HT>="; thisSelection += minHT;}
  else { thisSelection += " && HT>="; thisSelection += minHT;}
  cut=getCutString(kMC,thisSelection); //not compatible with scans!
  vectorOfCuts.push_back(cut);
  if (latexMode_) stageCut.push_back("HT$\\ge$"+minHT);
  else stageCut.push_back("HT>="+minHT);
  
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

  //MET
  if (thisSelection=="") {thisSelection += "MET>="; thisSelection += minMET; }
  else {thisSelection += " && MET>="; thisSelection += minMET;}
  cut=getCutString(kMC,thisSelection);
  vectorOfCuts.push_back(cut);
  if (latexMode_) stageCut.push_back("\\MET$\\ge$"+minMET);
  else stageCut.push_back("MET>="+minMET);

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
//   if (thisSelection=="") thisSelection += "cutCleaning==1";
//   else thisSelection +=" && cutCleaning==1";
//   cut=getCutString(lumiScale_,thisSelection);
//   vectorOfCuts.push_back(cut);
//   stageCut.push_back("Cleaning");

  //store selection string pre b cut
  selectionPreB=thisSelection;

  //loose selection
  //>= 1 b
  selectionGe1bLoose=selectionPreB;
  if (btagSF) btagSFweight_="probge1";
  else selectionGe1bLoose += " && nbjetsCSVM>=1";
  cut=getCutString(kMC,selectionGe1bLoose);
  vectorOfCuts.push_back(cut);
  if (latexMode_) stageCut.push_back("HT$\\ge$400, \\MET$\\ge$250, $\\ge$1 b");
  else stageCut.push_back("HT>=400, MET>=250, >= 1 b");

  //==1b
  //jmt -- don't bother with this
//   selectionEq1bLoose=selectionPreB; selectionEq1bLoose +=" && nbjetsSSVHPT==1";
//   cut=getCutString(lumiScale_,selectionEq1bLoose);
//   vectorOfCuts.push_back(cut);
//   if (latexMode_) stageCut.push_back("HT$\\ge$350, \\MET$\\ge$200, $==$1 b");
//   else stageCut.push_back("HT>=350, MET>=200, == 1 b");

  //>= 2 b
  selectionGe2bLoose=selectionPreB; 
  if (btagSF) btagSFweight_="probge2";
  else selectionGe2bLoose += " && nbjetsCSVM>=2";
  cut=getCutString(kMC,selectionGe2bLoose);
  vectorOfCuts.push_back(cut);
  if (latexMode_) stageCut.push_back("HT$\\ge$400, \\MET$\\ge$250, $\\ge$2 b");
  else stageCut.push_back("HT>=400, MET>=250, >= 2 b");

  //>= 3 b
  selectionGe3bLoose=selectionPreB; 
  if (btagSF) btagSFweight_="probge3";
  else selectionGe2bLoose += " && nbjetsCSVM>=3";
  cut=getCutString(kMC,selectionGe3bLoose);
  vectorOfCuts.push_back(cut);
  if (latexMode_) stageCut.push_back("HT$\\ge$400, \\MET$\\ge$250, $\\ge$3 b");
  else stageCut.push_back("HT>=400, MET>=250, >= 3 b");


  //tight selection
  if (!isTightSelection){ //print out the results for the tight selection anyway
    //>=1 b
    selectionGe1bTight=selectionPreB; 
    if (btagSF) {btagSFweight_="probge1"; selectionGe1bTight += " && HT>=500 && MET>=500";}
    else selectionGe1bTight += " && nbjetsCSVM>=1 && HT>=500 && MET>=500"; //hard-coded!
    cut=getCutString(kMC,selectionGe1bTight);
    vectorOfCuts.push_back(cut);
    if (latexMode_) stageCut.push_back("HT$\\ge$500, \\MET$\\ge$500, $\\ge$1 b");
    else stageCut.push_back("HT>=500, MET>=500, >= 1 b");
    
    //==1 b
//     selectionGe1bTight=selectionPreB; selectionGe1bTight += " && nbjetsSSVHPT==1 && HT>=500 && MET>=300"; //hard-coded!
//     cut=getCutString(lumiScale_,selectionGe1bTight);
//     vectorOfCuts.push_back(cut);
//     if (latexMode_) stageCut.push_back("HT$\\ge$500, \\MET$\\ge$300, $==$1 b");
//     else stageCut.push_back("HT>=500, MET>=300, == 1 b");
    
    //>=2 b
    selectionGe1bTight=selectionPreB;
    if (btagSF) {btagSFweight_="probge2"; selectionGe1bTight += " && HT>=600 && MET>=300";}
    else selectionGe1bTight += " && nbjetsCSVM>=2 && HT>=600 && MET>=300"; //hard-coded!
    cut=getCutString(kMC,selectionGe1bTight);
    vectorOfCuts.push_back(cut);
    if (latexMode_) stageCut.push_back("HT$\\ge$600, \\MET$\\ge$300, $\\ge$2 b");
    else stageCut.push_back("HT>=600, MET>=300, >= 2 b");  
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
  useHLTeff_=true;   
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

      //qcd reweighting not implemented yet

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
	cout<<format_nevents(cutflowEntries[istage][isample],cutflowEntriesE[istage][isample])<<col;
      }
    }//end loop over samples
    //now add total background, LM
    cout<<format_nevents(nSumSM[istage],nSumSME[istage]);
    for (unsigned int isample=0; isample<samples_.size(); isample++) {
      if ( !isSampleSM(samples_[isample]) ) cout<<col<<format_nevents(cutflowEntries[istage][isample],cutflowEntriesE[istage][isample]);
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
