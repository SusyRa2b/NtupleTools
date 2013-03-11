// -*- C++ -*-

#include "TStopwatch.h"

#include "TRegexp.h"

//classes now split into their own files
#include "ConfigurationDescriptions.h"
// #include "SignalEffData.h" //deprecated
//#include "SearchRegion.h" //deprecated

#include "TGraphAsymmErrors.h"

#include <utility>

//CMS style
TStyle *theStyle =0;
void initStyle() {

  //check if the style is already defined
  if (theStyle==0 && gROOT->GetStyle("CMS")==0) {
    theStyle = new TStyle("CMS","Style for P-TDR");


    // For the canvas:
    theStyle->SetCanvasBorderMode(0);
    theStyle->SetCanvasColor(kWhite);
    theStyle->SetCanvasDefH(600); //Height of canvas
    theStyle->SetCanvasDefW(600); //Width of canvas
    theStyle->SetCanvasDefX(0);   //POsition on screen
    theStyle->SetCanvasDefY(0);
    
    // For the Pad:
    theStyle->SetPadBorderMode(0);
    // theStyle->SetPadBorderSize(Width_t size = 1);
    theStyle->SetPadColor(kWhite);
    theStyle->SetPadGridX(false);
    theStyle->SetPadGridY(false);
    theStyle->SetGridColor(0);
    theStyle->SetGridStyle(3);
    theStyle->SetGridWidth(1);
    
    // For the frame:
    theStyle->SetFrameBorderMode(0);
    theStyle->SetFrameBorderSize(1);
    theStyle->SetFrameFillColor(0);
    theStyle->SetFrameFillStyle(0);
    theStyle->SetFrameLineColor(1);
    theStyle->SetFrameLineStyle(1);
    theStyle->SetFrameLineWidth(1);
    
    // For the histo:
    // theStyle->SetHistFillColor(1);
    // theStyle->SetHistFillStyle(0);
    theStyle->SetHistLineColor(1);
    theStyle->SetHistLineStyle(0);
    theStyle->SetHistLineWidth(1);
    // theStyle->SetLegoInnerR(Float_t rad = 0.5);
  // theStyle->SetNumberContours(Int_t number = 20);
    
    theStyle->SetEndErrorSize(2);
    //  theStyle->SetErrorMarker(20);
    theStyle->SetErrorX(0.);
    
    theStyle->SetMarkerStyle(20);
    
    //For the fit/function:
    theStyle->SetOptFit(1);
    theStyle->SetFitFormat("5.4g");
    theStyle->SetFuncColor(2);
    theStyle->SetFuncStyle(1);
    theStyle->SetFuncWidth(1);
    
    //For the date:
    theStyle->SetOptDate(0);
    // theStyle->SetDateX(Float_t x = 0.01);
    // theStyle->SetDateY(Float_t y = 0.01);
    
    // For the statistics box:
    theStyle->SetOptFile(0);
    theStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
    theStyle->SetStatColor(kWhite);
    theStyle->SetStatFont(42);
    theStyle->SetStatFontSize(0.025);
    theStyle->SetStatTextColor(1);
    theStyle->SetStatFormat("6.4g");
    theStyle->SetStatBorderSize(1);
    theStyle->SetStatH(0.1);
    theStyle->SetStatW(0.15);
    // theStyle->SetStatStyle(Style_t style = 1001);
    // theStyle->SetStatX(Float_t x = 0);
    // theStyle->SetStatY(Float_t y = 0);
    
    // Margins:
    theStyle->SetPadTopMargin(0.05);
    theStyle->SetPadBottomMargin(0.13);
    theStyle->SetPadLeftMargin(0.16);
    theStyle->SetPadRightMargin(0.02);
    
    // For the Global title:
    
    theStyle->SetOptTitle(0);
    theStyle->SetTitleFont(42);
    theStyle->SetTitleColor(1);
    theStyle->SetTitleTextColor(1);
    theStyle->SetTitleFillColor(10);
    theStyle->SetTitleFontSize(0.05);
    // theStyle->SetTitleH(0); // Set the height of the title box
    // theStyle->SetTitleW(0); // Set the width of the title box
    // theStyle->SetTitleX(0); // Set the position of the title box
    // theStyle->SetTitleY(0.985); // Set the position of the title box
    // theStyle->SetTitleStyle(Style_t style = 1001);
    // theStyle->SetTitleBorderSize(2);
    
    // For the axis titles:
    
    theStyle->SetTitleColor(1, "XYZ");
    theStyle->SetTitleFont(42, "XYZ");
    theStyle->SetTitleSize(0.06, "XYZ");
    // theStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
    // theStyle->SetTitleYSize(Float_t size = 0.02);
    theStyle->SetTitleXOffset(0.9);
    theStyle->SetTitleYOffset(1.25);
    // theStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset
    
    // For the axis labels:
    
    theStyle->SetLabelColor(1, "XYZ");
    theStyle->SetLabelFont(42, "XYZ");
    theStyle->SetLabelOffset(0.007, "XYZ");
    theStyle->SetLabelSize(0.05, "XYZ");
    
    // For the axis:
    
    theStyle->SetAxisColor(1, "XYZ");
    theStyle->SetStripDecimals(kTRUE);
    theStyle->SetTickLength(0.03, "XYZ");
    theStyle->SetNdivisions(510, "XYZ");
    theStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    theStyle->SetPadTickY(1);
    
    // Change for log plots:
    theStyle->SetOptLogx(0);
    theStyle->SetOptLogy(0);
    theStyle->SetOptLogz(0);
    
    // Postscript options:
    theStyle->SetPaperSize(20.,20.);
    // theStyle->SetLineScalePS(Float_t scale = 3);
    // theStyle->SetLineStyleString(Int_t i, const char* text);
    // theStyle->SetHeaderPS(const char* header);
    // theStyle->SetTitlePS(const char* pstitle);
    
    // theStyle->SetBarOffset(Float_t baroff = 0.5);
    // theStyle->SetBarWidth(Float_t barwidth = 0.5);
    // theStyle->SetPaintTextFormat(const char* format = "g");
    // theStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
    // theStyle->SetTimeOffset(Double_t toffset);
    // theStyle->SetHistMinimumZero(kTRUE);
    
    theStyle->cd(); //what does this do?
    //end CMS style
    
  }

}

//useful for playing around with plots in interactive ROOT
TH1D* hinteractive=0;
TH2D* h2d=0;

TLatex* text1=0;
TLatex* text2=0;
TLatex* text3=0;
TLatex* text4=0;

//holds a list of the *active* samples (to be plotted)
std::vector<TString> samples_;
//hold a list of all samples
std::set<TString> samplesAll_;
//config descriptions defined below
//these maps use the sample names as keys
std::map<TString, std::map<TString, std::pair<TFile*,TChain*> > > files_;
std::map<TString, TH1D*> histos_;
std::map<TString, UInt_t> sampleColor_;
std::map<TString, TString> sampleOwenName_;
std::map<TString, TString> sampleLabel_;
std::map<TString, UInt_t> sampleMarkerStyle_;
std::map<TString, UInt_t> sampleLineStyle_;
std::map<TString, TString> sampleWeightFactor_;//an arbitary string for weighting individual samples
std::map<TString, float> sampleScaleFactor_; //1 by default //implemented only for drawPlots!
std::map<TString, TChain*> primaryDatasets_;
std::vector<TLine*> verticalLines_; //for decorating plots with a set of vertical lines
std::vector<float> verticalLinePositions_; //for decorating plots with a set of vertical lines
TChain* dtree=0;
TH1D* hdata=0;
TH2D* hdata2d=0;

TString currentConfig_;
TH2D* scanSMSngen=0;
//TH1D* referenceCrossSectionGluino=0;
CrossSectionTable * CrossSectionTable_gluino_=0;
CrossSectionTable * CrossSectionTable_stop_=0;


//default selection
TString selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1";

TString nameOfEventWeight_="weight";
TString reducedTreeName_="reducedTree";

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


ConfigurationDescriptions configDescriptions_;






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

bool smsHack_=false; //horrible...

bool quiet_=false;
//bool quiet_=true;
bool doRatio_=false;
bool logy_=false;
bool logx_=false;
bool dostack_=true;
bool doleg_=true;
bool dodata_=true;
bool addOverflow_=true;
//bool doSubtraction_=false;
bool drawMCErrors_=false;
bool renormalizeBins_=false;//no setter function
bool owenColor_ = false;
bool drawFilenameOnPlot_=false;

int m0_=0;
int m12_=0;
TString susyCrossSectionVariation_="";

bool normalized_=false;

bool usePUweight_=false;
bool useTrigEff_ = false;

TString btagSFweight_="1";

bool savePlots_ = true; //no setter function
bool drawTotalSM_=false; //no setter function
bool drawTotalSMSusy_=false;//no setter function
bool drawSusyOnly_=false;//no setter function
bool drawMarkers_=true;//no setter function

bool useMassInLegend_=true;

bool doCustomPlotMax_=false;
double customPlotMax_=0;

bool doCustomPlotMin_=false;
double customPlotMin_=0;

float maxScaleFactor_ = 1.05;
float xlabeloffset_=0.015;

bool splitTTbar_ = false;
bool splitWJets_ = false;
//the next three are automatically configured in slABCD()
bool splitTTbarForClosureTest_ = false;
bool splitWJetsForClosureTest_ = false;
bool splitSingleTopForClosureTest_ = false;

//split the Zinvisible sample according to the neutrino acceptance categories
bool splitZinvisible_ = false;
bool vAccepMuon_ = false; //true(false) = muon(electron) definition of acceptance

//this is for plotting headers and also for dd estimates
bool is2011Data_ = false;
bool isPreliminary_ = (is2011Data_) ? false : true;

bool plotSelectionLabel_ = true;
TString selectionLabel_="";
bool doAppendBinWidth_ = true;

//bool latexMode_=false;
bool latexMode_=true;
const TString pm = latexMode_ ? " \\pm " : " +/- ";

TCanvas* thecanvas=0;
//TCanvas* cratio=0;
TLegend* leg=0;
THStack* thestack=0;
TLatex* extraText=0;
TString extratext_="";
double rLine = -1;
TH1D* totalsm=0;
TH1D* totalsmsusy=0;
TH1D* totalewk=0;
TH1D* totalqcdttbar=0;
TH1D* totalnonttbar=0;
TH1D* totalnonqcd=0;
TH1D* totalqcd=0; //ben - just for ease of doing event counts with drawPlots
TH1D* ratio=0; float ratioMin=0.0; float ratioMax=2.5;
TLine* ratioLine=0;
TGraphErrors* mcerrors=0;
TGraphAsymmErrors* effgraph=0;
bool loaded_=false; //bookkeeping
bool loadedSusyHistos_=false;//bookkeeping

//initialized tchains that point to MC samples
void resetChains() {

  cout<<"[resetting chains]"<<endl;

  //loop over files and set the TChain to point to the reducedTree in that files
  for ( std::map<TString, std::map<TString, std::pair<TFile*,TChain*> > >::iterator ifile = files_.begin();
	ifile != files_.end(); ++ifile) {
    cout<<"\t"<<ifile->first<<endl;

    std::map<TString, std::pair<TFile*,TChain*> > thisguy = ifile->second;

    for ( std::map<TString, std::pair<TFile*,TChain*> >::iterator icfg = thisguy.begin(); icfg!= thisguy.end(); ++icfg) {
      if (icfg->second.first !=0 && (!icfg->second.first->IsZombie())) {
	cout<<"\t"<<icfg->first<<"   "<<icfg->second.first->GetName()<<endl ;
	//now we've got the std::pair
	//	files_[ifile->first][icfg->first].second = (TChain*) icfg->second.first->Get(reducedTreeName_);
	if (files_[ifile->first][icfg->first].second != 0) delete files_[ifile->first][icfg->first].second;
	files_[ifile->first][icfg->first].second = new TChain(reducedTreeName_);
	files_[ifile->first][icfg->first].second->Add( files_[ifile->first][icfg->first].first->GetName() );
      }
      else files_[ifile->first][icfg->first].second =0;
    }
  }

}

TString stripSamplename(TString fullname) {

  //actually, this is ok
  //  if ( fullname.Contains("$") && fullname.Contains(":") )   cout<<"Simultaneous use of $ and : syntax not currently allowed!"<<endl;

  //eg for samplename T1tttt$900$100, return T1tttt
  TString name = fullname.Tokenize("$")->At(0)->GetName();
  //also allow the TTbarJets:cutstring syntax
  name = name.Tokenize(":")->At(0)->GetName();

  return name;
}


TChain* getTree(const TString & samplename) {

  TChain* t = (TChain*) files_[currentConfig_][samplename].second;
  return t;
}

//plot multiple "samples" as one
void chainSamples(const TString & parent, const TString & child) {

  if (parent.Contains("$") || parent.Contains(":") ) {
    cout<<" cannot chain subsamples defined with $ or :"<<endl;
    return;
  }

  //first enforce some sanity
  //the parent must be plotted, and the child must not be plotted
  //the might be a bit overzealous, especially the first part. we'll see

  bool foundit=false;
  for (unsigned int isample = 0; isample<samples_.size(); isample++) {
    if (samples_[isample] == parent || stripSamplename(samples_[isample])==parent) { foundit=true;}
    if (samples_[isample] == child) {
      cout<<child<<" is already on the plotting list! Why do you want to plot it twice? quitting..."<<endl;
      return;
    }
  }

  if (!foundit) {
    cout<<parent<<" is not on the plotting list! Not chaining samples...."<<endl;
    return;
  }


  cout<<"Adding "<<child<<" to "<<parent<<". n entries goes from "<<getTree(parent)->GetEntries()<<" to ";

  getTree(parent)->Add( getTree(child));
  cout<<getTree(parent)->GetEntries()<<endl;

}


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

TLine* getVerticalLine(unsigned int index) {  return verticalLines_[index];}
void addVerticalLine(float position) {
  verticalLinePositions_.push_back(position);
}

// void showDataMinusMC(bool dosub) {
//   doSubtraction_=dosub;
// }

void deleteVerticalLines() {
  for (unsigned int il=0; il<verticalLines_.size(); il++) {
    cout<<" DEBUG deleting line "<<il<<" "<< verticalLines_.at(il)<<endl;
    delete verticalLines_.at(il);
  }
}

void resetVerticalLine() {
  //this can crash the code if run at the wrong time. why?
  deleteVerticalLines();
  verticalLines_.clear();
  verticalLinePositions_.clear();
}

void doOverflowAddition(bool doOv) {
  addOverflow_ = doOv;
}

void setLogY(bool dolog) {
  logy_=dolog;
  if (logy_) maxScaleFactor_=3;
  else maxScaleFactor_=1.05;
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

//CMSSM/mSugra stuff is definitely not up-to-date -- may not work anymore
CrossSectionTable * CrossSectionTable_mSUGRAtanb40_=0;
void loadSusyCrossSections() {
  if (CrossSectionTable_mSUGRAtanb40_==0) {
    CrossSectionTable_mSUGRAtanb40_ = new  CrossSectionTable("NLOxsec_tanb40_10.txt");
  }

}

bool isSampleSMS(const TString & name ) {
  if (name.Contains("T1bbbb")) return true;
  if (name.Contains("T2bb")) return true;
  if (name.Contains("T1tttt")) return true;
  if (name.Contains("T2tt")) return true;

  return false;
}

bool isSampleScan(const TString & name ) {

  if (isSampleSMS(name)) return true;

  if (name.Contains("SUGRA")) return true;


  return false;
}

bool isSampleSM(const TString & name) {

  if (name.Contains("LM")) return false;

  if (name.Contains("sbottom")) return false;

  if (name.Contains("SUGRA")) return false;
  if (isSampleSMS(name)) return false;

  return true;
}

TString extractExtraCut(TString fullname) {
  if ( !fullname.Contains(":") ) return "";

  TString cuts = fullname.Tokenize(":")->At(1)->GetName();
  cuts.Prepend("(");
  cuts += ")";

  return cuts;
}

void setScanPoint(const TString & name) {

  // format is e.g. T1bbbb$m0$m12

  if ((name.Tokenize("$")->At(1) == 0) || (name.Tokenize("$")->At(2)==0) ) assert(0);

  TString m0=name.Tokenize("$")->At(1)->GetName();
  TString m12=name.Tokenize("$")->At(2)->GetName();
  m0_=m0.Atoi();
  m12_=m12.Atoi();
}


int getSMSProdProcess(const TString & sample) {

  if (sample.Contains("T1")) return 8; //gluino-gluino production
  assert(0); //T2 needs implementation

  return 10;
}


void loadScanSMSngen(const TString& sampleOfInterest) {
  if (scanSMSngen==0) scanSMSngen = (TH2D*) files_[currentConfig_][sampleOfInterest].first->Get("scanSMSngen");
  else {
    cout<<"Replacing scanSMSngen"<<endl;
    scanSMSngen = (TH2D*) files_[currentConfig_][sampleOfInterest].first->Get("scanSMSngen");
  }
}

void loadReferenceCrossSections() {
  //bad idea to hard-code to this remote afs?
  if ( CrossSectionTable_gluino_ == 0) {
    cout<<"Loading reference cross section file (g~g~)"<<endl;
    CrossSectionTable_gluino_ = new  CrossSectionTable("/afs/hephy.at/user/w/walten/public/referenceXSecs.root","smsroot","gluino8TeV_NLONLL");
  }

  if (CrossSectionTable_stop_ == 0) {
    cout<<"Loading reference cross section file (t~t~)"<<endl;
    CrossSectionTable_stop_ = new  CrossSectionTable("/afs/hephy.at/user/w/walten/public/referenceXSecs.root","smsroot","stop8TeV_NLONLL");
  }

}

// jmt Dec 2012 -- honestly not sure if this is needed anymore. leave it for now.
//m0, m12            //pdf set
map<pair<int,int>, map<TString, TH2D*> >  scanProcessTotalsMap;
void loadSusyScanHistograms() {
  if (loadedSusyHistos_) return;


  TFile* susyfile = 0;
  for (unsigned int isample=0; isample<samples_.size(); isample++) {
    if ( !isSampleSM(samples_[isample])) {
      susyfile =  files_[currentConfig_][samples_[isample]].first;
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
  text1 = new TLatex(5,23.08044,"CMS Simulation");
  text1->SetNDC();
  text1->SetTextAlign(13);
  text1->SetX(0.6);
  text1->SetY(.9);
  text1->SetTextFont(42);
  text1->SetTextSizePixels(24);
  text1->Draw();
}

void drawPlotHeader(double xoffset = 0) {

  const  bool publicationFormat = true; //to switch on the format suggested by Bill for the paper

  //  return;
  float ypos = 0.97;
  //  if(doRatio_) ypos=ypos+0.012;
  // i'm gonna leave this out for now
  if (text1 != 0 ) delete text1;
  if(isPreliminary_){
    text1 = new TLatex(3.570061,23.08044,"CMS Preliminary"); 
    text1->SetX(0.68 + xoffset ); //add 0.2 if you get rid of the "Preliminary"
  }
  else{ 
    text1 = new TLatex(3.570061,23.08044,"CMS"); //no more preliminary!
    text1->SetX(0.68 + xoffset +0.2); 
  }
  text1->SetNDC();
  text1->SetTextAlign(13);
  text1->SetY(ypos+0.007);
  text1->SetTextFont(42);
  text1->SetTextSizePixels(24);
  text1->SetTextSize(0.045); //copied from ben's code. maybe needs a switch so it is only used for AN2011()
  if (!publicationFormat)  text1->Draw();

  if (normalized_ == false) {
    TString astring;
    if(is2011Data_){
      //astring.Form("%.0f pb^{-1} at #sqrt{s} = 7 TeV",lumiScale_);
      //astring.Form("%.1f fb^{-1} at #sqrt{s} = 7 TeV",lumiScale_/1000.);      
      if (publicationFormat)    astring.Form("CMS, #sqrt{s} = 7 TeV, L_{int} = %.2f fb^{-1}",lumiScale_/1000.);
      else     astring.Form("L_{int} = %.2f fb^{-1}, #sqrt{s} = 7 TeV",lumiScale_/1000.);

      astring.Form("CMS, #sqrt{s} = 7 TeV, L_{int} = %.2f fb^{-1}",lumiScale_/1000.);

    }
    else astring.Form("CMS Preliminary, L_{int} = %.2f fb^{-1}, #sqrt{s} = 8 TeV",lumiScale_/1000.);
    if(lumiScale_>33. && lumiScale_<34.) astring.Form("L_{int} = %.2f fb^{-1}, #sqrt{s} = 7 TeV", 4982.91/1000.);//hardcoded, but don't know what else to do for this...
    if (text2 != 0 ) delete text2;
    text2 = new TLatex(3.570061,23.08044,astring);
    text2->SetNDC();
    text2->SetTextAlign(13);
    text2->SetX(0.17);
    const float extray= publicationFormat ? 0.01: 0;
    text2->SetY(ypos+0.015 + extray);
    text2->SetTextFont(42);
    text2->SetTextSizePixels(24);
    text2->SetTextSize(0.045); //copied from ben's code. maybe needs a switch so it is only used for AN2011()
    text2->Draw();
  }

  if(plotSelectionLabel_){
    if (text4 != 0 ) delete text4;
    text4 = new TLatex(5,23.08044,selectionLabel_);
    text4->SetNDC();
    text4->SetTextAlign(13);
    text4->SetX(0.2);
    text4->SetY(.89);
    text4->SetTextFont(42);
    //text4->SetTextSizePixels(26);
    text4->SetTextSize(0.06);
    text4->Draw();
  }
}

void setSelectionLabel(TString label){
  selectionLabel_ = label;
}



int mainpadWidth; int mainpadHeight;
int ratiopadHeight = 150;
// TPad* mainPad=0;
// TPad* ratioPad=0;
void renewCanvas(const TString opt="") {
  if (thecanvas!=0) delete thecanvas;

  int canvasWidth = mainpadWidth;
  int canvasHeight = opt.Contains("ratio") ? mainpadHeight+ratiopadHeight : mainpadHeight;

  //  if (drawCutsOnPlot_) canvasHeight += 50;

  thecanvas= new TCanvas("thecanvas","the canvas",canvasWidth,canvasHeight);
  thecanvas->cd()->SetRightMargin(0.04);
  thecanvas->cd()->SetTopMargin(0.08); //test

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
    thecanvas->GetPad(1)->SetTopMargin(0.07); //0.06
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
    if (logx_) thecanvas->GetPad(1)->SetLogx();
  }
  else { 
    if (logy_) thecanvas->SetLogy(); 
    if (logx_) thecanvas->SetLogx();
}


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

void renewExtraText() {
  if (extraText!=0) delete extraText;
  extraText = new TLatex();
  extraText->SetNDC();
  extraText->SetTextAlign(13);
  extraText->SetX(0.5);
  extraText->SetY(.85);
  extraText->SetTextFont(42);
  extraText->SetTextSizePixels(24);
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


//jmt Dec 2012 -- clearly this function is stylistically awful, but it does the job and I don't want to mess with it right now
//try to bring some rationality to the options passed to getcutstring....
enum sampleType {kData, kMC, kmSugraPoint, kmSugraPlane, kSMSPointGlGl, kSMSPointStopStop, kSMSPointSbottomSbottom, kSMSPlane};
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
  else if (type==kMC || type==kmSugraPoint || type==kmSugraPlane ||type==kSMSPointGlGl || type==kSMSPointStopStop) lumiscale=lumiScale_;
  else if (type==kSMSPlane) lumiscale=1;
  else {assert(0);}

  TString weightedcut=nameOfEventWeight_;
  
  weightedcut += "*(";
  weightedcut +=lumiscale;
  weightedcut+=")";

  //new scale factor
  weightedcut += "*(";
  weightedcut +=scalefactor;
  weightedcut+=")";

  //horrible kludge
  /*
  if (type!=kData) {
    if (thisSelection.Contains("pass_PFHT350_PFMET100")) thisSelection.ReplaceAll("pass_PFHT350_PFMET100","1");
    else  if (thisSelection.Contains("pass_PFHT350_Mu15_PFMET45")) thisSelection.ReplaceAll("pass_PFHT350_Mu15_PFMET45","1");
    else  if (thisSelection.Contains("pass_IsoMu24_eta2p1")) thisSelection.ReplaceAll("pass_IsoMu24_eta2p1","1");
    else  if (thisSelection.Contains("pass_IsoMu24")) thisSelection.ReplaceAll("pass_IsoMu24","1");
  }
  */

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
  if (useTrigEff_ &&  type!=kData) {
    weightedcut +="*hlteff";
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
    if (type == kmSugraPoint || type==kSMSPointGlGl|| type==kSMSPointStopStop) {
      if (smsHack_) weightedcut += " && mStop=="; //for stop samples where i used the wrong names
      else          weightedcut += " && m0==";
      weightedcut +=m0_;
      if (smsHack_) weightedcut += " && mLSP=="; //for stop samples where i used the wrong names
      else          weightedcut += " && m12==";
      weightedcut +=m12_;
    }
    weightedcut+=")";
  }
  else if (extraSelection !="") {
    weightedcut += "*(";
    weightedcut +=extraSelection;
    if (type == kmSugraPoint || type==kSMSPointGlGl|| type==kSMSPointStopStop) {
      if (smsHack_)  weightedcut += " && mStop=="; //for stop samples where i used the wrong names
      else           weightedcut += " && m0==";
      weightedcut +=m0_;
      if (smsHack_) weightedcut += " && mLSP=="; //for stop samples where i used the wrong names
      else          weightedcut += " && m12=="; 
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
  else if (type == kSMSPointGlGl || type == kSMSPointStopStop ||type == kSMSPointSbottomSbottom ) {
    if (!smsHack_) {
    //need to weight with the standard sigma / N ; lumi factor is done above
    assert(scanSMSngen!=0);
    int bin=  scanSMSngen->FindBin(m0_,m12_);
    double ngen=scanSMSngen->GetBinContent(bin);
    loadReferenceCrossSections();
    double xs =0;
    if (type == kSMSPointGlGl) xs = CrossSectionTable_gluino_->getSMSCrossSection(m0_);
    else if (type == kSMSPointStopStop) xs = CrossSectionTable_stop_->getSMSCrossSection(m0_);
    else assert(0); //have not yet implemented sbottom
    char xsweight[100];
    sprintf(xsweight,"*(%f/%f)",xs,ngen);
    weightedcut += xsweight;
    }
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
  if (thecanvas==0 || verticalLinePositions_.size()==0) return;

  deleteVerticalLines(); //clean up old lines

  thecanvas->Update(); //important!

  int padindex = doRatio_ ? 1 : 0;
  TVirtualPad* thePad = thecanvas->GetPad(padindex);
  thecanvas->cd(padindex);
  //  double xmin,ymin,xmax,ymax;
  //  thePad->GetRangeAxis(xmin,ymin,xmax,ymax);
  double ymax=  thePad->GetUymax();
  double ymin=  thePad->GetUymin();
  //for academic interest, can get the same numbers using e.g. thePad->GetUymax()
  if (logy_) {
    ymax = pow(10, ymax);
    ymin = pow(10, ymin);
  }

  for (unsigned int il = 0; il<verticalLinePositions_.size(); il++) {
    TLine* theLine = new TLine(verticalLinePositions_.at(il),ymin,verticalLinePositions_.at(il),ymax);
    theLine->SetLineColor(kBlack); //hard-coded for now...
    theLine->SetLineWidth(2);
    theLine->Draw();
    verticalLines_.push_back(theLine);
  }


}

void setSampleColor(const TString & sample, UInt_t color)  {  sampleColor_[sample]=color;}
void setSampleLineStyle(const TString & sample, UInt_t style)  {  sampleLineStyle_[sample]=style;}

UInt_t getSampleColor(const TString & sample)  { 
  if ( sampleColor_.count(sample) ) return sampleColor_[sample];
  //if we don't find it, return the default
  return 1;
}

//add a sample to be plotted to the *end* of the list
void addSample(const TString & newsample, const int color=-1, const TString label="") {

  //see if it is already there
  for (std::vector<TString>::iterator it = samples_.begin(); it!=samples_.end(); ++it) {
    if ( *it == newsample) {
      cout<<newsample<<" is already on the list!"<<endl;
      return;
    }
  }
  //if it isn't there, go ahead and add it

  if ( samplesAll_.find(stripSamplename(newsample)) != samplesAll_.end() ) {
    samples_.push_back(newsample);
    if (color>=0) sampleColor_[newsample]=color;
    if (label!="") sampleLabel_[newsample] = label;
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

void resetSampleWeightFactors() {

  for (std::set<TString>::iterator isample=samplesAll_.begin(); isample!=samplesAll_.end(); ++isample) {
    sampleWeightFactor_[*isample] = "";
  }
}

void setSampleScaleFactor(const TString & sample, const float sf) {
  bool done=false;
  for (std::vector<TString>::iterator isample=samples_.begin(); isample!=samples_.end(); ++isample) {
    if (*isample == sample) {sampleScaleFactor_[*isample] = sf; done=true;}
  }
  if (!done) cout<<"Failed to find the sample "<<sample<<endl;
}

float getSampleScaleFactor(const TString & sample) {

  if (  sampleScaleFactor_.count(sample)>0)    return sampleScaleFactor_[sample];
  else if ( sampleScaleFactor_.count( stripSamplename(sample))>0) return sampleScaleFactor_[stripSamplename(sample)];

  return 1;
}

void setSampleWeightFactor(const TString & sample, const TString & factor) {
  bool done=false;
  for (std::vector<TString>::iterator isample=samples_.begin(); isample!=samples_.end(); ++isample) {
    if (*isample == sample) {sampleWeightFactor_[*isample] = factor; done=true;}
  }
  if (!done) cout<<"Failed to find the sample "<<sample<<endl;
}

void printSampleWeightFactors() {
  for (std::map<TString, TString>::iterator thisfactor = sampleWeightFactor_.begin(); thisfactor!=sampleWeightFactor_.end(); ++thisfactor) {
    cout<<thisfactor->first <<" "<<thisfactor->second<<endl;
  }

}

TString getSampleWeightFactor(const TString & sample) {

  if (  sampleWeightFactor_.count(sample)>0)      return sampleWeightFactor_[sample];
  else if ( sampleWeightFactor_.count( stripSamplename(sample))>0)   return sampleWeightFactor_[stripSamplename(sample)];


  return "";
}

void clearSamples() {
  if (!quiet_) cout<<"cleared sample list!"<<endl;
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
    sampleColor_["sbottom-185-250"]=kRed-7;
    sampleColor_["sbottom-189-270"]=kRed-3;
    sampleColor_["sbottom-217-300"]=kRed+2;
    sampleColor_["mSUGRAtanb40"] =kGray;
    sampleColor_["T1bbbb"] =kBlack;
    sampleColor_["T1tttt"] =kBlack;
    sampleColor_["T2bb"] =kGray;
    sampleColor_["T2tt"] =kGray;
    sampleColor_["QCD"] = kYellow;
    sampleColor_["PythiaQCD"] = kYellow;
    sampleColor_["PythiaPUQCD"] = kYellow;
    sampleColor_["PythiaPUQCDFlat"] = kYellow;
    sampleColor_["TTbarJets"]=kAzure-3;
    sampleColor_["TTbarJets0"]=kAzure-3;
    sampleColor_["TTbarJetsPowheg"]=kAzure-3;
    sampleColor_["TTbarJetsMCNLO"]=kAzure-3;
    sampleColor_["TTbarSingleTopWJetsCombined"]=kAzure-3;
    sampleColor_["ttbar"]=kAzure-3;
    sampleColor_["SingleTop"] = kMagenta;
    sampleColor_["TTV"] =kBlue+3 ;
    sampleColor_["WJets"] = kGreen-3;
    sampleColor_["ZJets"] = kViolet-3;
    sampleColor_["Zinvisible"] = kOrange-3;
    sampleColor_["TotalSM"] = kBlue+2;
    sampleColor_["Total"] = kGreen+3;
    sampleColor_["VV"] = kCyan+1;
    sampleColor_["HerwigQCDFlat"] = kYellow;
  }
  else if (name == "nostack" || name=="owen") {
    sampleColor_["LM13"] = kBlue+2;
    sampleColor_["LM9"] = kCyan+2;
    sampleColor_["sbottom-185-250"]=kRed-7;
    sampleColor_["sbottom-189-270"]=kRed-3;
    sampleColor_["sbottom-189-270"]=kRed+2;
    sampleColor_["mSUGRAtanb40"] =kCyan+2;
    sampleColor_["T1bbbb"] =kCyan+2;
    sampleColor_["T1tttt"] =kCyan+2;
    sampleColor_["T2bb"] =kCyan+2;
    sampleColor_["T2tt"] =kCyan+2;
    sampleColor_["QCD"] = 2;
    sampleColor_["PythiaQCD"] = 2;
    sampleColor_["PythiaPUQCD"] =2;
    sampleColor_["PythiaPUQCDFlat"] =2;
    sampleColor_["TTbarJets"]=kBlue;
    sampleColor_["TTV"] =kBlue+4 ;
    sampleColor_["TTbarJets0"]=kCyan;
    sampleColor_["TTbarJetsPowheg"]=kCyan+4;
    sampleColor_["TTbarJetsMCNLO"]=kViolet;
    sampleColor_["SingleTop"] = kMagenta;
    sampleColor_["WJets"] = kOrange;
    sampleColor_["ZJets"] = 7;
    sampleColor_["Zinvisible"] = kOrange+7;
    sampleColor_["TotalSM"] =kGreen+1; //owen requested 3
    sampleColor_["Total"] = 6;
    sampleColor_["VV"] = kOrange-3;
    sampleColor_["HerwigQCDFlat"] = 2;
  }
  else {
    cout<<"Sorry, color scheme "<<name<<" is not known!"<<endl;
  }

}

void setStackMode(bool dostack, bool normalized=false,bool adjustcolors=true) {
  dostack_=dostack;
  normalized_=normalized;

  if (adjustcolors) { //sometimes I don't want the code to be this clever!
  if (dostack) setColorScheme("stack");
  else setColorScheme("nostack");
  }
}

void resetSamples(bool joinSingleTop=true) {

  assert(joinSingleTop); //no longer support any other option

  samples_.clear();
  //this block controls what samples will enter your plot
  //order of this vector controls order of samples in stack

  //careful -- QCD must have 'QCD' in its name somewhere.
  //samples_.push_back("QCD"); //madgraph
  //samples_.push_back("PythiaQCD");
  samples_.push_back("PythiaPUQCD"); 
  //samples_.push_back("PythiaPUQCDFlat");

  samples_.push_back("TTbarJets");

  samples_.push_back("WJets");
  samples_.push_back("SingleTop");
  samples_.push_back("ZJets");
  samples_.push_back("VV");
  samples_.push_back("Zinvisible");
  //samples_.push_back("HerwigQCDFlat");
  //  samples_.push_back("LM9");

}

bool overrideSMSlabels_=false; //another dirty hack
TString getSampleLabel(const TString & sample) {
  TString label="Not Found";
  if (  isSampleScan(sample) && !overrideSMSlabels_){
    char ss[50];
    //more compact format?
    TString shortname = stripSamplename(sample);    
    if( useMassInLegend_) {
      setScanPoint(sample); //set m0 and m12
      //    sprintf(ss, "#splitline{m_{gluino} = %d GeV}{m_{LSP} = %d GeV}",m0_,m12_);
      sprintf(ss, "%s (%d GeV, %d GeV)",shortname.Data(),m0_,m12_);
    }
    else{
      sprintf(ss, "%s ",shortname.Data());
    }
    label=ss;  
  }
  else {
    if (sampleLabel_.count(sample))  label = sampleLabel_[sample];
    else label = sample;
  }
  return label;

}


void addDataToChain(TChain* datachain, const TString & name) {

  TString dname="reducedTree.";
  dname+=currentConfig_;
  dname += ".";
  dname += name;
  dname += "*.root";
  dname.Prepend(dataInputPath);
  dname.ReplaceAll("JERbias","JER0"); //JERbias not relevant for data
  datachain->Add(dname);

}

void setDatasetToDraw(const TString & dataset) {

  if (primaryDatasets_.count(dataset) == 0) {
    cout<<"Dataset not found: "<<dataset<<endl;
    return;
  }
  dtree = primaryDatasets_[dataset];

}

void addToSamplesAll(const TString & name) {
  samplesAll_.insert(name);
}

void loadSamples(bool joinSingleTop=true, TString signalEffMode="") {
  if (loaded_) return;
  loaded_=true; //implies that loadSamples() can only be called once per session

  initStyle();
  resetPadDimensions();
  resetLegendPosition();

  resetSamples(joinSingleTop);
  //samplesAll_ should have *every* available sample
  //however we now have addToSamplesAll() to allow the user to add some stuff before calling loadSamples()
  //also note that there's no harm in failing to load one of these samples, 
  //as long as you don't actually try to draw it

  samplesAll_.insert("QCD1000");
  samplesAll_.insert("QCD120");
  samplesAll_.insert("QCD1400");
  samplesAll_.insert("QCD170");
  samplesAll_.insert("QCD1800");
  samplesAll_.insert("QCD300");
  samplesAll_.insert("QCD470");
  samplesAll_.insert("QCD600");
  samplesAll_.insert("QCD800");
  samplesAll_.insert("WJets_HT250To300");
  samplesAll_.insert("WJets_HT300To400");
  samplesAll_.insert("WJets_HT400ToInf");
  samplesAll_.insert("WW");
  samplesAll_.insert("WZ");
  samplesAll_.insert("ZZ");
  samplesAll_.insert("ZJets_HT200To400");
  samplesAll_.insert("ZJets_HT400ToInf");
  samplesAll_.insert("Zinvisible_HT100To200");
  samplesAll_.insert("Zinvisible_HT200To400");
  samplesAll_.insert("Zinvisible_HT400ToInf");

  //samplesAll_.insert("PythiaQCD");
  samplesAll_.insert("PythiaPUQCD");
  //samplesAll_.insert("PythiaPUQCDFlat");
  samplesAll_.insert("TTbarJets");
  samplesAll_.insert("TTbarJets0");
  samplesAll_.insert("TTbarJetsPowheg");
  samplesAll_.insert("TTbarJetsMCNLO");
  samplesAll_.insert("TTbarEmbedFlipRecoMu");
  samplesAll_.insert("TTbarEmbedFlipRecoEle");
  samplesAll_.insert("TTbarEmbedGenTauHad");
  samplesAll_.insert("TTTo2L2Nu2B");
  samplesAll_.insert("TTbarJets8TeV");
  samplesAll_.insert("TTbarSingleTopWJetsCombined");
  samplesAll_.insert("WJets");
  samplesAll_.insert("WJetsInc");
  samplesAll_.insert("ZJets");
  samplesAll_.insert("ZJetsInc");
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
  samplesAll_.insert("HerwigQCDFlat");
  //  samplesAll_.insert("WJetsZ2");
  //  samplesAll_.insert("ZJetsZ2");
  samplesAll_.insert("VV");
  samplesAll_.insert("TTV");


  samplesAll_.insert("LM13");
  samplesAll_.insert("LM9");

  samplesAll_.insert("mSUGRAtanb40");
  samplesAll_.insert("T1bbbb");
  samplesAll_.insert("T1tttt");
  samplesAll_.insert("T2bb");
  samplesAll_.insert("T2tt");

  samplesAll_.insert("sbottom-185-250"); 
  samplesAll_.insert("sbottom-189-270"); 
  samplesAll_.insert("sbottom-217-300"); 

  //since we don't use JERbias anymore this is basically irrelevant, but let's keep this framework
  ////////////
  if (signalEffMode=="") {
    configDescriptions_.setDefault("CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0");
    configDescriptions_.setCorrected("CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0");
  }
  else if (signalEffMode=="ra2b2012") {
     configDescriptions_.setDefault("CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0");
     configDescriptions_.setCorrected("CSVM_PF2PATjets_JES0_JERbias_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0");
  }
  else if (signalEffMode=="stop") {
    configDescriptions_.setDefault("default");
    configDescriptions_.setCorrected("corrected");
  }
  else if (signalEffMode=="scan" || signalEffMode=="complete") {  //Only for signal systematics [not used anymore in 2012]
      
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
  sampleLabel_["sbottom-185-250"]="sbottom185-250";
  sampleLabel_["sbottom-189-270"]="sbottom189-270";
  sampleLabel_["sbottom-217-300"]="sbottom217-300";
  sampleLabel_["QCD"] = "QCD";
  sampleLabel_["PythiaQCD"] = "QCD (no PU)";
  sampleLabel_["PythiaPUQCDFlat"] = "QCD"; 
  sampleLabel_["PythiaPUQCD"] = "QCD";
  sampleLabel_["TTbarJets"]="t#bar{t}";
  sampleLabel_["TTbarJets0"]="t#bar{t} (small MG)";
  sampleLabel_["TTbarJetsPowheg"]="t#bar{t} (Powheg)";
  sampleLabel_["TTbarJetsMCNLO"]="t#bar{t} (MC@NLO)";
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
  sampleLabel_["WJetsInc"] = "W#rightarrowl#nu (all HT)";
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
  sampleLabel_["ZJetsInc"] = "Z/#gamma*#rightarrowl^{+}l^{-} (all HT)";
  sampleLabel_["Zinvisible"] = "Z#rightarrow#nu#nu";
  sampleLabel_["Zinvisible-2vAcc"] = "Z#rightarrow#nu#nu - 2Acc";
  sampleLabel_["Zinvisible-1vAcc"] = "Z#rightarrow#nu#nu - 1Acc";
  sampleLabel_["Zinvisible-0vAcc"] = "Z#rightarrow#nu#nu - 0Acc";
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

  sampleLineStyle_["T1bbbb"] = 1;
  sampleLineStyle_["T1tttt"] = 7;
  sampleLineStyle_["T2tt"] = 1;
  sampleLineStyle_["T2bb"] = 1;
  sampleLineStyle_["LM9"] = 1;
  sampleLineStyle_["sbottom-185-250"]=1;
  sampleLineStyle_["sbottom-189-270"]=1;
  sampleLineStyle_["sbottom-217-300"]=1;

  sampleMarkerStyle_["mSUGRAtanb40"] = kFullStar;
  sampleMarkerStyle_["T1bbbb"] = kFullStar;
  sampleMarkerStyle_["T1tttt"] = kFullStar;
  sampleMarkerStyle_["T2bb"] = kFullStar;
  sampleMarkerStyle_["T2tt"] = kFullStar;
  sampleMarkerStyle_["LM13"] = kFullStar;
  sampleMarkerStyle_["sbottom-185-250"]=kFullStar;
  sampleMarkerStyle_["sbottom-189-270"]=kFullStar;
  sampleMarkerStyle_["sbottom-217-300"]=kFullStar;
  sampleMarkerStyle_["LM9"] = kFullStar;
  sampleMarkerStyle_["QCD"] = kFullCircle;
  sampleMarkerStyle_["PythiaQCD"] = kOpenCircle;
  sampleMarkerStyle_["PythiaPUQCDFlat"] = kOpenCircle;  
  sampleMarkerStyle_["PythiaPUQCD"] = kOpenCircle;
  sampleMarkerStyle_["TTbarJets"]= kFullSquare;
  sampleMarkerStyle_["TTbarJets0"]= kFullCircle;
  sampleMarkerStyle_["TTbarJetsPowheg"]= kOpenCircle;
  sampleMarkerStyle_["TTbarJetsMCNLO"]= kOpenCircle;
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
  sampleMarkerStyle_["WJetsInc"] = kMultiply;
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
  sampleMarkerStyle_["ZJetsInc"] = kFullTriangleDown;
  sampleMarkerStyle_["Zinvisible"] = kFullTriangleDown;
  sampleMarkerStyle_["Zinvisible-2vAcc"] = kFullTriangleDown;
  sampleMarkerStyle_["Zinvisible-1vAcc"] = kFullTriangleDown;
  sampleMarkerStyle_["Zinvisible-0vAcc"] = kFullTriangleDown;
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

  //jmt Dec 2012 -- do we use these for anything anymore?
  sampleOwenName_["mSUGRAtanb40"] = "msugra40";
  sampleOwenName_["T1bbbb"] = "t1bbbb";
  sampleOwenName_["T1tttt"] = "t1tttt";
  sampleOwenName_["T2bb"] = "t2bb";
  sampleOwenName_["T2tt"] = "t2tt";
  sampleOwenName_["LM13"] = "lm13";
  sampleOwenName_["LM9"] = "lm9";
  sampleOwenName_["sbottom-189-270"]="sbottom189-270";
  sampleOwenName_["sbottom-185-250"]="sbottom185-250";
  sampleOwenName_["sbottom-217-300"]="sbottom217-300";
  sampleOwenName_["QCD"] = "qcd";
  sampleOwenName_["PythiaQCD"] = "qcd";
  sampleOwenName_["PythiaPUQCDFlat"] = "qcd"; 
  sampleOwenName_["PythiaPUQCD"] = "qcd";
  sampleOwenName_["TTbarJets"]="ttbar";
  sampleOwenName_["TTbarJetsPowheg"]="ttbar";
  sampleOwenName_["TTbarSingleTopWJetsCombined"]="ttbartw";
  sampleOwenName_["ttbar"]="ttbar";
  sampleOwenName_["SingleTop"] = "singletop";
  sampleOwenName_["WJets"] = "wjets";
  sampleOwenName_["WJetsInc"] = "wjetsinc";
  sampleOwenName_["WJetsZ2"] = "wjets";
  sampleOwenName_["ZJets"] = "zjets";
  sampleOwenName_["ZJetsInc"] = "zjetsinc";
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
  resetSampleWeightFactors();

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

      TFile* dummyp1=0;
      TChain* dummyp2=0;

      files_[thisconfig][*isample] = make_pair(dummyp1,dummyp2);
      files_[thisconfig][*isample].first = new TFile(fname);
      if (files_[thisconfig][*isample].first->IsZombie() ) {cout<<"file error with "<<*isample<<endl; files_[thisconfig][*isample].first=0;}
      else { if (!quiet_)    cout<<"Added sample: "<<thisconfig<<"\t"<<*isample<<endl;}
    }
  }
  //now set the TChain pointers in files_
  resetChains();

  //load data files too
  //use Run2012B names
  primaryDatasets_["HTMHT"]=0;
  primaryDatasets_["JetHT"]=0;
  //  primaryDatasets_["SingleMu"]=0;
  //  primaryDatasets_["SingleElectron"]=0;
  //  primaryDatasets_["MuHad"]=0;
  //  primaryDatasets_["ElectronHad"]=0;
  //  primaryDatasets_["DoubleElectron"]=0;
  //  //  primaryDatasets_["DoubleMu"]=0;
  //  primaryDatasets_["MuEG"]=0;
  //  primaryDatasets_["Photon"]=0;
  primaryDatasets_["MET"]=0;
 
  primaryDatasets_["2012hybrid"] = 0;
  primaryDatasets_["2012hybridplus"] = 0;

  //there are others, but this is what I have now
  //2011 names  
//   primaryDatasets_["2011HT"]=0;
//   primaryDatasets_["2011SingleMu"]=0;
//   primaryDatasets_["2011SingleElectron"]=0;
//   primaryDatasets_["2011DoubleMu"]=0;
//   primaryDatasets_["2011DoubleElectron"]=0;


  for (map<TString, TChain*>::iterator idataset=primaryDatasets_.begin(); idataset!=primaryDatasets_.end(); ++idataset) {
    if (idataset->first == "HTMHT") {
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"HT_Run2012A");
      addDataToChain(idataset->second,"HTMHT_Run2012B");
      addDataToChain(idataset->second,"HTMHT_Run2012C");
      addDataToChain(idataset->second,"HTMHT_Run2012D");
    }
    else if (idataset->first == "2012hybrid") {
      //HTMHT + MET, to be used with care to avoid duplicate events
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"HT_Run2012A");
      addDataToChain(idataset->second,"MET_Run2012A");
      addDataToChain(idataset->second,"HTMHT_Run2012B");
      addDataToChain(idataset->second,"MET_Run2012B");
      addDataToChain(idataset->second,"HTMHT_Run2012C");
      addDataToChain(idataset->second,"MET_Run2012C");
      addDataToChain(idataset->second,"HTMHT_Run2012D");
      addDataToChain(idataset->second,"MET_Run2012D");
    }
    else if (idataset->first == "2012hybridplus") {
      //HTMHT + MET, to be used with care to avoid duplicate events
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"HT_Run2012A");
      addDataToChain(idataset->second,"MET_Run2012A");
      addDataToChain(idataset->second,"HTMHT_Run2012B");
      addDataToChain(idataset->second,"MET_Run2012B");
      addDataToChain(idataset->second,"HTMHT_Run2012C");
      addDataToChain(idataset->second,"MET_Run2012C");
      addDataToChain(idataset->second,"HTMHT_Run2012D");
      addDataToChain(idataset->second,"MET_Run2012D");
      //for B+C+D, also use JetHT
      addDataToChain(idataset->second,"JetHT_Run2012B");
      addDataToChain(idataset->second,"JetHT_Run2012C");
      addDataToChain(idataset->second,"JetHT_Run2012D");
    }
    else if (idataset->first == "JetHT") {
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"HT_Run2012A");
      addDataToChain(idataset->second,"JetHT_Run2012B");
      addDataToChain(idataset->second,"JetHT_Run2012C");
      addDataToChain(idataset->second,"JetHT_Run2012D");
    }
    else if (idataset->first == "SingleMu") {
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"SingleMu_Run2012A");
      addDataToChain(idataset->second,"SingleMu_Run2012B");
    }
    else if (idataset->first == "SingleElectron") {
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"SingleElectron_Run2012A");
      addDataToChain(idataset->second,"SingleElectron_Run2012B");
    }
    else if (idataset->first == "Photon") {
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"Photon_Run2012A");
      cout<<"Warning -- Photon Run2012B not here"<<endl;
    }
    else if (idataset->first == "MET") {
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"MET_Run2012A");
      addDataToChain(idataset->second,"MET_Run2012B");
      addDataToChain(idataset->second,"MET_Run2012C");
    }
    else if (idataset->first == "MuHad") {
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"MuHad_Run2012A");
      cout<<"Warning -- MuHad Run2012B not here"<<endl;
      //      addDataToChain(idataset->second,"MuHad_Run2012B");
    }
    else if (idataset->first == "ElectronHad") {
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"ElectronHad_Run2012A");
      addDataToChain(idataset->second,"ElectronHad_Run2012B");
    }
    else if (idataset->first == "DoubleElectron") {
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"DoubleElectron_Run2012A");
      addDataToChain(idataset->second,"DoubleElectron_Run2012B");
    }
    else if (idataset->first == "DoubleMu") {
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"DoubleMu_Run2012A");
      addDataToChain(idataset->second,"DoubleMu_Run2012B");
    }
    else if (idataset->first == "MuEG") {
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"MuEG_Run2012A");
      addDataToChain(idataset->second,"MuEG_Run2012B");
    }
    else if (idataset->first == "2011HT") {
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"ht");
    }
    else if (idataset->first == "2011SingleMu") {
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"singlemu");
    }
    else if (idataset->first == "2011SingleElectron") {
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"singleelectron");
    }
    else if (idataset->first == "2011DoubleMu") {
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"doublemu");
    }
    else if (idataset->first == "2011DoubleElectron") {
      idataset->second = new TChain(reducedTreeName_);
      addDataToChain(idataset->second,"doubleelectron");
    }
    else {
      cout<<" Can't find "<<idataset->first<<endl;
      assert(0);
    }
  }

  dtree = primaryDatasets_["HTMHT"]; //set some default

}

//ug...not liking this
sampleType getSampleType(const TString & sample , const TString & planeOrPoint="") {

  assert(planeOrPoint=="" || planeOrPoint=="point" || planeOrPoint=="plane");

  if (sample=="data") return kData;
  else if (sample.Contains("mSUGRA") && planeOrPoint=="point") return kmSugraPoint;
  else if (sample.Contains("mSUGRA") && planeOrPoint=="plane") return kmSugraPlane;
  else if (sample.Contains("T1bbbb") && planeOrPoint=="point") return kSMSPointGlGl;
  else if (sample.Contains("T1bbbb") && planeOrPoint=="plane") return kSMSPlane;
  else if (sample.Contains("T1tttt") && planeOrPoint=="point") return kSMSPointGlGl;
  else if (sample.Contains("T1tttt") && planeOrPoint=="plane") return kSMSPlane;
  else if (sample.Contains("T2bb") && planeOrPoint=="point") return kSMSPointSbottomSbottom;
  else if (sample.Contains("T2bb") && planeOrPoint=="plane") return kSMSPlane;
  else if (sample.Contains("T2tt") && planeOrPoint=="point") return kSMSPointStopStop;
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
		 const TString histname , const TString samplename, const float* varbins=0 ,
		 const TString xtitle="",const TString ytitle="") {

  loadSamples();
  gROOT->SetStyle("CMS");
//I would rather implement this functionality via drawPlots(), but I think it will be simpler
//to just write something simple

//no presentation, just fill the histogram and save
  TChain* tree=0;
  if (samplename=="data") {
    tree = dtree;
  }
  else {
    for (unsigned int isample=0; isample<samples_.size(); isample++) {
      if ( samples_[isample] == samplename) {
	if (!quiet_) cout <<samples_[isample]<<endl;
	//	tree = (TChain*) files_[currentConfig_][stripSamplename(samples_[isample])]->Get(reducedTreeName_);
	tree = getTree(stripSamplename(samples_[isample]));
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
  hh->SetXTitle(xtitle);
  hh->GetYaxis()->SetTitle(ytitle);

  if (samplename=="data") tree->Project(histname,var,getCutString(true).Data());
  else //tree->Project(histname,var,getCutString( getSampleType(samplename,"point"),optfh,selection_,"",0,"").Data());
    tree->Project(histname,var,getCutString( getSampleType(samplename,"point"),getSampleWeightFactor(samplename),selection_,extractExtraCut(samplename),0,"",-1,getSampleScaleFactor(samplename)).Data());

  //  hh->Scale( getSampleScaleFactor(samplename));

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
		 const TString histname, const TString samplename, const TString xtitle="",const TString ytitle="") {
  return drawSimple(var, nbins, 0, 1, filename, histname, samplename, varbins, xtitle, ytitle);
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
      //      TChain* tree = (TChain*) files_[currentConfig_][samples_[isample]]->Get(reducedTreeName_);
      TChain* tree = getTree(samples_[isample]);
      gROOT->cd();
      
      h2d_temp->Reset(); //just to be safe
      tree->Project("h2d_temp",drawstring,getCutString( getSampleType(samples_[isample],"point"),getSampleWeightFactor(samples_[isample]),selection_,extractExtraCut(samples_[isample]),0,"",-1,getSampleScaleFactor(samples_[isample])).Data());

      h2d->Add(h2d_temp);
    }
  }
  h2d->SetXTitle(xtitle);
  h2d->SetYTitle(ytitle);
  h2d->SetLineColor(getSampleColor("TotalSM"));
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
  }
  if (doRatio_) {
    if (ratio!=0) delete ratio;
    ratio = (varbins==0) ? new TH1D("ratio","data/(SM MC)",nbins,low,high) : new TH1D("ratio","",nbins,varbins);
    ratio->Sumw2();
    ratioLine = new TLine(low, 1, high, 1);
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

  totalsm->SetMarkerColor(getSampleColor("TotalSM"));
  totalsm->SetLineColor(getSampleColor("TotalSM"));
  totalsm->SetLineWidth(2);
  totalsm->SetMarkerStyle(sampleMarkerStyle_["TotalSM"]);
  if (!drawMarkers_)  totalsm->SetMarkerSize(0); //no marker for this one

  totalsmsusy->SetMarkerColor(getSampleColor("Total"));
  totalsmsusy->SetLineColor(getSampleColor("Total"));
  totalsmsusy->SetLineWidth(2);
  totalsmsusy->SetMarkerStyle(sampleMarkerStyle_["Total"]);
  if (!drawMarkers_)  totalsmsusy->SetMarkerSize(0); //no marker for this one

  if (doAppendBinWidth_ && varbins==0 && !renormalizeBins_) {ytitle = appendBinWidth(ytitle,low,high,nbins,unit);} //add the " / 10 GeV" part to the title

  //here is the part that is really different from the previous implementation
  //need to make new histograms
  resetHistos(); //delete existing histograms
  TString opt="hist e";
  double histMax=-1e9;
  vector<TString> signals;
  for (unsigned int isample=0; isample<samples_.size(); isample++) {

    //the full sample name can contain "extra" stuff like the mass point and a cut. Get the base sample name
    const    TString samplename= stripSamplename(samples_[isample]);

    if (!quiet_)   cout <<samples_[isample]<<endl;

    gROOT->cd();
    //ive each histo have a different name
    TString hname = jmt::fortranize(var); hname += "_"; hname += jmt::fortranize(samples_[isample]);
    histos_[samples_[isample]] = (varbins==0) ? new TH1D(hname,"",nbins,low,high) : new TH1D(hname,"",nbins,varbins);
    histos_[samples_[isample]]->Sumw2();
    if( !files_[currentConfig_][samplename].first ){
      cout << "ERROR: TFile for " << samples_[isample] << " is null.  Check that the file exists.  Exiting." << endl;
      assert(0);
    }
    TChain* tree = getTree(samplename);
    if (tree==0) {
      cout<<"tree is null! Something must be wrong!"<<endl<<currentConfig_<<endl<<samplename<<endl;
    }
    gROOT->cd();

    //replace the scanSMSngen with the appropriate one for this sample
    if ( isSampleSMS(samplename)) loadScanSMSngen(samplename);
    if (isSampleScan(samplename)) setScanPoint(samples_[isample]);

    //fill the histogram!
    tree->Project(hname,var,getCutString( getSampleType(samples_[isample],"point"),getSampleWeightFactor(samples_[isample]),selection_,extractExtraCut(samples_[isample]),0,"",-1,getSampleScaleFactor(samples_[isample])).Data());

    //now the histo is filled
    if (renormalizeBins_) ytitle=renormBins(histos_[samples_[isample]],2 ); //manipulates the TH1D //FIXME hard-coded "2"
    if (addOverflow_)  addOverflowBin( histos_[samples_[isample]] ); //manipulates the TH1D
    histos_[samples_[isample]]->SetXTitle(xtitle);
    histos_[samples_[isample]]->SetYTitle(ytitle);
    histos_[samples_[isample]]->GetXaxis()->SetLabelOffset(xlabeloffset_);

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
    if (!isSampleSM(samples_[isample]))    signals.push_back(samples_[isample]);
    //add everything! (but only the 1st signal)
    if (signals.size()<=1)  {
      totalsmsusy->Add(histos_[samples_[isample]]);
    }
    //now just do a bunch of histogram formatting
    if (!dostack_) {
      //set line color instead of fill color for this type of plot
      histos_[samples_[isample]]->SetLineColor(getSampleColor(samples_[isample]));
      histos_[samples_[isample]]->SetMarkerStyle(sampleMarkerStyle_[samples_[isample]]);
      histos_[samples_[isample]]->SetMarkerColor(getSampleColor(samples_[isample]));
      if (!drawMarkers_) histos_[samples_[isample]]->SetMarkerSize(0);

      //ad hoc additions
      histos_[samples_[isample]]->SetLineWidth(2);
    }
    else {
      if ( isSampleSM(samples_[isample]))   histos_[samples_[isample]]->SetFillColor(getSampleColor(samples_[isample]));
      if ( !isSampleSM(samples_[isample]))  {
	histos_[samples_[isample]]->SetLineWidth(2);
	histos_[samples_[isample]]->SetLineColor(getSampleColor(samples_[isample]));
	histos_[samples_[isample]]->SetLineStyle(sampleLineStyle_[samples_[isample]]);
      }
      histos_[samples_[isample]]->SetMarkerSize(0);
    }

    if (dostack_) { //add histo to stack
      //but *don't* add more than one signal to the stack
      if (isSampleSM(samples_[isample]) || signals.size()<=1) thestack->Add(histos_[samples_[isample]] );
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
  if (drawTotalSMSusy_) leg->AddEntry(totalsmsusy, sampleLabel_["Total"],"F");
  if (drawTotalSM_) leg->AddEntry(totalsm, sampleLabel_["TotalSM"],"F");
  for (int isample=int(samples_.size())-1; isample>=0; isample--) {
    if (dostack_ || (!drawSusyOnly_ || samples_[isample].Contains("LM")|| samples_[isample].Contains("SUGRA"))) { //drawSusyOnly_ means don't draw SM
      leg->AddEntry(histos_[samples_[isample]], getSampleLabel(samples_[isample]),"F");
    }
  }

  if (!dostack_) {
    //this is all a re-implemenataion of stuff done is HistHolder. Oh well.

    //also unit normalize the totalsm
    if (normalized_ && totalsm->Integral()>0) totalsm->Scale(1.0/totalsm->Integral()); //15 Jan 2013

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
    thestack->GetHistogram()->GetXaxis()->SetLabelOffset(xlabeloffset_);
    thestack->GetHistogram()->GetYaxis()->SetTitle(ytitle);
    thestack->GetHistogram()->GetXaxis()->SetTitleOffset(0.93);

    //drawverticalline used to be called here but i've now moved it lower

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

    //ok, if we've got "extra" signals we need to draw them here (so skip the first one)
    for (unsigned int isig=1; isig<signals.size(); isig++) {
      //first we need a name
      TString extrasignalname="totalsmPlus";
      extrasignalname+= signals.at(isig);
      //now clone the totalsm and add the signal to it
      TH1D* newtotal = (TH1D*) totalsm->Clone(extrasignalname);
      newtotal->Add( histos_[signals.at(isig)]);
      //now plop the newtotal into histos_
      histos_[signals.at(isig)] = newtotal;
      //FIXME need to adjust the plot max in case this is higher than the stack

      //now draw
      cout<<"drawing "<<signals.at(isig)<<endl;
      newtotal->SetLineColor(getSampleColor(signals.at(isig)));
      newtotal->SetLineStyle(sampleLineStyle_[signals.at(isig)]);
      newtotal->Draw("SAME HIST");
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
    hdata->Draw("SAME E0");
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
  } //if dodata_
  if (doRatio_) {
      TString ratioytitle = "Data/MC ";
      thecanvas->cd(2);
      if (dodata_) {
	ratio->Divide(hdata,totalsm);
	if (!quiet_) cout<<"plotting data/MC ratio"<<endl;
      }
      else { //kind of a dirty hack
	unsigned int sampleRatioIndex1 =0;
	unsigned int sampleRatioIndex2 =1;
	if (!quiet_) cout<<"plotting ratio "<<samples_[sampleRatioIndex1] <<" / "<<samples_[sampleRatioIndex2]<<endl;
	ratioytitle = ""; //for now just clear the title
	if (histos_[samples_[sampleRatioIndex1]] != 0 && histos_[samples_[sampleRatioIndex2]] != 0)
	  ratio->Divide(histos_[samples_[sampleRatioIndex1]],histos_[samples_[sampleRatioIndex2]]);
	else
	  cout<<"error: not plotting ratio because one of the histograms doesn't exist"<<endl;
      }
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
      ratio->GetYaxis()->SetTitle(ratioytitle);
      ratio->GetYaxis()->SetTitleSize(0.16);
      ratio->GetYaxis()->SetTitleOffset(.43);
      gPad->SetTopMargin(1e-5);
      gPad->SetTickx();
      gPad->Modified();      

      ratio->Draw();

      //thecanvas->GetPad(2)->SetTopMargin(0.1);

      ratioLine->Draw();

    
  }

  thecanvas->cd(mainPadIndex);
  if (doleg_)  leg->Draw();
  drawPlotHeader();

  if (drawFilenameOnPlot_) {
    if (text3!=0) delete text3;
    text3=new TLatex(5,23,filename.Data());
    text3->SetNDC();
    text3->SetX(0.2);
    text3->SetY(0.89);
    text3->SetTextFont(42);
    text3->SetTextSize(0.03);
    text3->Draw();
  }

  drawVerticalLine();

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
    //    thecanvas->Print(savename+".C");    //for formal purposes
    thecanvas->SaveAs(savename+".pdf"); //for pdftex
    thecanvas->SaveAs(savename+".png"); //for twiki
  }

  //dump some event counts to the screen
  if (!quiet_) {
    int ns=0;
    for (unsigned int isample=0; isample<samples_.size(); isample++) {
      if (!isSampleSM( samples_[isample])) ns++;
      if ( isSampleSM( samples_[isample]) || ns<=1 || !dostack_) 
	cout<<samples_[isample]<<" =\t "<<histos_[samples_[isample]]->Integral()<<" +/- "<<jmt::errOnIntegral(histos_[samples_[isample]])<<endl;
      else  //starting with the 2nd signal, the total sm will be included, so we need to get rid of it (only for dostack_ true)
	cout<<samples_[isample]<<" =\t "<<histos_[samples_[isample]]->Integral()-totalsm->Integral()<<endl;
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
  renewExtraText();

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
    if (extratext_!="") {
      extraText->DrawLatex(.5,.5,extratext_);
      extraText->Draw();
    }
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

  totalsm->SetMarkerColor(getSampleColor("TotalSM"));
  totalsm->SetLineColor(getSampleColor("TotalSM"));
  totalsm->SetLineWidth(2);
  totalsm->SetMarkerStyle(0);
  totalsm->SetYTitle(ytitle);

  TString drawopt="hist e";
  float max=-1e9; TString firsthist="";
  for (unsigned int isample=0; isample<samples_.size(); isample++) {

    if (!quiet_) cout <<samples_[isample]<<endl;
    //    TChain* tree = (TChain*) files_[currentConfig_][samples_[isample]]->Get(reducedTreeName_);
    TChain* tree = getTree(samples_[isample]);

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
    tree->Project(hnameP,var,getCutString(getSampleType(samples_[isample],"point"),getSampleWeightFactor(samples_[isample]),selection_,cstring1).Data());
    tree->Project(hnameF,var,getCutString(getSampleType(samples_[isample],"point"),getSampleWeightFactor(samples_[isample]),selection_,cstring2).Data());

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
    if (!quiet_) cout<<"setting color to: "<<getSampleColor(samples_[isample])<<endl;
    histos_[hnameR]->SetLineColor(getSampleColor(samples_[isample]));
    histos_[hnameR]->SetMarkerStyle(sampleMarkerStyle_[samples_[isample]]);
    histos_[hnameR]->SetMarkerColor(getSampleColor(samples_[isample]));
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
    
    if (hnameR.Contains("QCD")) { //HACK draw only qcd
      histos_[hnameR]->Draw(drawopt);
      if (!drawopt.Contains("same")) drawopt+=" same";
      
      if (firsthist=="") firsthist = hnameR;
      if (histos_[hnameR]->GetMaximum() > max) max = histos_[hnameR]->GetMaximum();
      leg->AddEntry(histos_[hnameR], sampleLabel_[samples_[isample]]);
    }
     
    if(dodata_) drawPlotHeader(-.1);
    else{
      drawPlotHeaderInside(); //for qcd only plots, this works.
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
  if (extratext_!="") {
      extraText->DrawLatex(.5,.85,extratext_);
      extraText->Draw();
  }
  TLine* myRline = 0;
  if (rLine>0) { //this should actually be allowed to go negative
    myRline = new TLine( histos_[firsthist]->GetXaxis()->GetBinLowEdge(1), rLine, histos_[firsthist]->GetXaxis()->GetBinUpEdge(histos_[firsthist]->GetXaxis()->GetLast() ), rLine);
    myRline->SetLineWidth(2);
    myRline->SetLineStyle(7);
    myRline->Draw();
  }
  
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

//code from gala (I think) to make cut flow tables.
//I will not remove it for now, but I think it is quite outdated
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
  useTrigEff_=true;   
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

      //      TChain* tree = (TChain*) files_[currentConfig_][samples_[isample]]->Get(reducedTreeName_);
      TChain* tree = getTree(samples_[isample]);
      gROOT->cd();
      TString weightopt= "";
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

void drawEffVRej(const TString & h1, const TString h2,const TString & xtitle,const TString & ytitle,bool flip=false) {

  //step over the bins in histos_[h1] and histos_[h2]

  //first find h1, h2
  if ( histos_.count(h1) ==0 || (histos_.count(h2) ==0 && h2!="totalsm") ) {
    cout<<" Problem! "<<histos_.count(h1)<<" "<<histos_.count(h2)<<endl;
    return;
  }

  TH1D* hist1=histos_[h1];
  TH1D* hist2= (h2=="totalsm") ? totalsm : histos_[h2];

  const  int nbins=  hist1->GetNbinsX();

  if (effgraph!=0) delete effgraph;
  effgraph = new TGraphAsymmErrors(nbins );

  //silly little code to find best S/root(B)
  double  bestSrootB = 0; //depends on S being on x axis
  double bestCutVal=-99;  double bestX=-99;  double bestY=-99;
  int bestBin=-1;
  for ( int ibin=1; ibin<=nbins; ibin++) {

    double x,y;
    if (flip) {
      x=  hist1->Integral(ibin,nbins) / hist1->Integral();
      y=  hist2->Integral(ibin,nbins) / hist2->Integral();
    }
    else {
      x=  hist1->Integral(1,ibin) / hist1->Integral();
      y=  hist2->Integral(1,ibin) / hist2->Integral();
    }
    //alternate
    effgraph->SetPoint(ibin, x,y);

    if ( y>0 && (x / sqrt(y) > bestSrootB)) {
      bestSrootB = x/sqrt(y);
      bestCutVal = hist1->GetBinLowEdge(ibin);
      bestX=x; bestY=y; bestBin=ibin;
    }

    //    if (!quiet_)   cout<<ibin<<" "<<histos_[h1]->GetBinLowEdge(ibin)<<" \t "<<x<<" "<<y<<endl;
  }

  cout<<"Assuming S is on the x axis, B is on the y axis, best S/root(B) is "<<bestSrootB
      <<" "<<bestBin<<" "<<bestCutVal<<" \t "<<bestX<<" "<<bestY<<endl;

  renewCanvas();
  effgraph->Draw("Apl");
  effgraph->GetHistogram()->SetXTitle(xtitle);
  effgraph->GetHistogram()->SetYTitle(ytitle);
  effgraph->GetHistogram()->SetMaximum(1);
  effgraph->GetHistogram()->SetMinimum(0);

}

void drawTrigEff(const TString & pd, const TCut & tag, const TCut & probe, const TString & var, int nbins, float low,float high) {
  loadSamples(true,"ra2b2012");

  //see ~/scratch0/analysis2/CMSSW_4_2_5/src/NtupleTools/BasicLoopCU/splitByTrigger.C
  renewCanvas();

  setDatasetToDraw(pd);

  TCut thecuts = selection_.Data();

  TH1D hnum("hnum","numerator",nbins,low,high);
  TH1D hden("hden","denominator",nbins,low,high);
  hnum.Sumw2();
  hden.Sumw2();

  dtree->Project("hnum",var,thecuts && probe && tag);
  dtree->Project("hden",var,thecuts && tag);

  if (effgraph != 0) delete effgraph;
  effgraph = new TGraphAsymmErrors();
  effgraph->BayesDivide(&hnum,&hden);

  if (ratioLine!=0) delete ratioLine;
  ratioLine = new TLine(low, 1, high, 1);
  ratioLine->SetLineColor(2);

  effgraph->Draw("AP");
  effgraph->GetHistogram()->SetXTitle(var);
  effgraph->GetHistogram()->SetMaximum(1.1);

  if (doCustomPlotMax_)     effgraph->GetHistogram()->SetMaximum(customPlotMax_);
  if (doCustomPlotMin_)     effgraph->GetHistogram()->SetMinimum(customPlotMin_);

  ratioLine->Draw();
  effgraph->Draw("p");

  TString desc;
  desc.Form("%s,%s",pd.Data(),tag.GetTitle());
  if (text1!=0) delete text1;
  if (text2!=0) delete text2;
  if (text3!=0) delete text3;
  text1 = new TLatex(5,23.08044,desc.Data());
  text1->SetNDC();
  text1->SetX(0.3);
  text1->SetY(0.15);
  text1->SetTextFont(42);
  text1->SetTextSize(0.025);
  text1->Draw();
  text2 = new TLatex(5,23.08044,probe.GetTitle());
  text2->SetNDC();
  text2->SetX(0.3);
  text2->SetY(0.2);
  text2->SetTextFont(42);
  text2->SetTextSize(0.025);
  text2->Draw();

  text3=new TLatex(5,23,selection_.Data());
  text3->SetNDC();
  text3->SetX(0.15);
  text3->SetY(0.875);
  text3->SetTextFont(42);
  text3->SetTextSize(0.022);
  text3->Draw();



  TString filename;
  filename.Form("trigEff_%s-%s_%s_%s",pd.Data(),jmt::fortranize(tag.GetTitle()).Data(),jmt::fortranize(probe.GetTitle()).Data(),var.Data());
  if (savePlots_) {
    thecanvas->SaveAs(filename+".eps");
    thecanvas->SaveAs(filename+".png");
    thecanvas->SaveAs(filename+".pdf");
  }
}
