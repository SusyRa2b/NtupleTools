/*
to compile:

must have symlink to MiscUtil.cxx in the working directory.
TSelectorMultiDraw and CrossSectionTable need to be compiled only when they are changed.

gSystem->Load("TSelectorMultiDraw_C.so");
gSystem->Load("CrossSectionTable_cxx.so");
.L drawReducedTrees.C++

====== this is the nominal code for drawing RA2b data/MC comparison plots =======
The signal efficiency functions are now forked into signalEffSyst.C

2011 updates -- runDataQCD2011() computes all numbers needed for QCD estimate including systematics

-- drawPlots() -- main plotting routine to draw pretty stacks of MC with data on top
-- drawSimple() -- very simple function to draw exactly one sample and put it in a file
-- drawR() -- draws r(MET). a bit more kludgely, but works.
-- drawSignificance -- draws S/sqrt(B) as a function of the variable

The functions that do the heavy lifting, along with various utility functions,
are in drawReducedTrees.h.

This code replaces drawCutflowPlots.C (which had replaced drawBasicPlots.C). 
The input files are the reducedTrees.

Depending on the exact setup in loadSamples(), the QCD and Single-top samples must be added together with 'hadd' 
in order to get one file per sample.

*/

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TCut.h"
#include "TLatex.h"
#include "TChain.h"

//from Luke
#include "TSelectorMultiDraw.h"

//container for mSugra cross-sections
#include "CrossSectionTable.h"

// can be checked out of UserCode/joshmt
#include "MiscUtil.cxx"

#include <fstream>

#include <iostream>
#include <iomanip>
#include <map>
#include <set>
	
//*** AFTER SUMMER
//***************************

  //TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-35g/";
//TString dataInputPath =  "/cu2/ra2b/reducedTrees/V00-02-35h/";
TString inputPath = "/cu2/wteo/reducedTrees/Fall11/V00-02-35k/";
TString dataInputPath =  "/cu2/ra2b/reducedTrees/V00-02-35k/";

//TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-35d/";
//TString inputPath = "/cu1/joshmt/reducedTrees/test/"; 

//-- reducedTrees for Oct 25 SUSY meeting. 3464.581/pb. 
//TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-35a/";
//TString dataInputPath =  "/cu2/ra2b/reducedTrees/V00-02-35a/";


//*** SUMMER RESULT
//***************************
//for rerunning the final background estimates...need standard MC and MET cleaning
//V00-02-24_fullpf2pat was used for 13 Aug update of data results
//stick with this version for now because 25 is missing one ZJets sample and that takes away one event from the LDP
//TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-24_fullpf2pat/";//path for MC
//TString dataInputPath = "/cu3/wteo/reducedTrees/V00-02-05_v3-pickevents/"; //includes MET cleaning but uses a tight skim (not good for plots)

//for making the standard set of plots...need standard MC and all data.
//TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-25_fullpf2pat/";//path for MC	     
//TString dataInputPath = "/cu2/ra2b/reducedTrees/V00-02-25_fullpf2pat/";

//for njet reweighting (needed because not all reducedTrees in nominal V00-02-25 have njets30)
//TString inputPath = "/cu2/ra2b/reducedTrees/benV00-02-25_fullpf2pat/";
//TString dataInputPath = "/cu2/ra2b/reducedTrees/benV00-02-25_fullpf2pat/";

//for signal systematics
  //   TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-25c_fullpf2pat/"; //LM9 with correct pdf weights
//TString inputPath = "/cu2/joshmt/reducedTrees/V00-02-25c_fullpf2pat/"; //with correct pdf weights
//TString inputPath = "/home/joshmt/";//path for MC
//TString dataInputPath = "/cu2/ra2b/reducedTrees/V00-02-24_fullpf2pat/"; //sym links to V00-02-05_v3

//the cutdesc string is now defined in loadSamples()

//double lumiScale_ = 1143; //official summer conf lumi
//double lumiScale_ = 3464.581;//oct25
double lumiScale_ = 4683.719;//nov4

#include "drawReducedTrees.h"

const bool reweightLSBdata_=true; //whether or not LSB data is reweighted based on PV distribution
const bool useScaleFactors_=true; //whether or not to use MC scale factors when doing subtraction for data-driven estimates
const bool useBNNEffCurves_=false; 
const bool btaggedLSB_=true;
const bool use1B_SL_=false; // use ge1b selection for SL sample in ttbar method

void RSL(double SIG, double SB) {

  double rsl = SIG/SB;
  double err = jmt::errAoverB(SIG,sqrt(SIG),SB,sqrt(SB));

  cout<<"R_SL = "<<rsl<< " +/- "<<err<<endl;

}

void countInBoxesBreakdown(const SearchRegion & region) {
  //shows number of events from each background and signal in the 6 boxes
  
  useFlavorHistoryWeights_=false;
  loadSamples();


  //dodata_=false;
  savePlots_ =false;
  setQuiet(true);
  doOverflowAddition(true);  
  //region.Print();
  
  //for susy scan; harmless for non-susy scan
  resetSamples();
  
  TString btagselection = region.btagSelection;
  TCut HTcut=region.htSelection.Data(); 
  TCut SRMET = region.metSelection.Data();
  bool isSIG = region.isSIG;

   
  TCut bcut = "1";
  if (btagselection=="ge1b") {
    bcut="nbjetsCSVM>=1";
  }
  else if (btagselection=="ge2b") {
    bcut="nbjetsCSVM>=2";
  }
  else if (btagselection=="ge3b") {
    bcut="nbjetsCSVM>=3";
  }
  else if (btagselection=="eq1b") {
    bcut="nbjetsCSVM==1";
  }
  else if (btagselection=="eq2b") {
    bcut="nbjetsCSVM==2";
  }
  else {assert(0);}
  
  usePUweight_=false;
  useHTeff_=false;
  useMHTeff_=false;
  thebnnMHTeffMode_ = kOff;
  btagSFweight_="1";
  currentConfig_=configDescriptions_.getDefault(); //completely raw MC


  if (useScaleFactors_) {
    bcut="1";
    usePUweight_=true;
    useHTeff_=true;
    useMHTeff_=true;
    currentConfig_=configDescriptions_.getCorrected(); //add JERbias

    if (btagselection=="ge2b") {
      btagSFweight_="probge2";
    }
    else if (btagselection=="ge1b") {
      btagSFweight_="probge1";
    }
    else if (btagselection=="eq1b") {
      btagSFweight_="prob1";
    }
    else if (btagselection=="eq2b") {
      btagSFweight_="(1-prob1-prob0-probge3)"; //because we don't have prob2 in the reducedTrees
    }
    else if (btagselection=="ge3b") {
      btagSFweight_="probge3";
    }
    else {assert(0);}
  }

  TString thisbox="";
  
  TCut baseline = "cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passCleaning==1";
  TCut baselineSL = "cutPV==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0 && MT_Wlep<100 && passCleaning==1";
  TCut passOther = "minDeltaPhiN>=4";
  TCut failOther="minDeltaPhiN<4";

  TString var="HT"; TString xtitle=var;
  int  nbins=1; float low=0; float high=1e9;
  char output[500];
  //Look at: PythiaPUQCD TTbarJets SingleTop WJets ZJets VV Zinvisible LM9
  
  //double QCD, QCDerr, ttbar, ttbarerr, singletop, singletoperr, wjets, wjetserr, zjets, zjetserr, vv, vverr, zinv, zinverr, lm9, lm9err, data, dataerr, totalsm, totalsmerr;

  // --- nominal (SB or SIG) selection ---
  selection_ = baseline && HTcut && passOther && SRMET && bcut; 
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  thisbox=""; thisbox+= isSIG ? "$N_{SIG}$":"$N_{SB}$";
  sprintf(output,"%s & %s%s & %s & %s & %s & %s & %s & %s & %s & %s& %i & %s \\\\",thisbox.Data(),region.btagSelection.Data(),region.owenId.Data(),
	  jmt::format_nevents(getIntegral("PythiaPUQCD"),getIntegralErr("PythiaPUQCD")).Data(), jmt::format_nevents(getIntegral("TTbarJets"),getIntegralErr("TTbarJets")).Data(),
	  jmt::format_nevents(getIntegral("SingleTop"),getIntegralErr("SingleTop")).Data(), jmt::format_nevents(getIntegral("WJets"),getIntegralErr("WJets")).Data(),
	  jmt::format_nevents(getIntegral("ZJets"),getIntegralErr("ZJets")).Data(), jmt::format_nevents(getIntegral("VV"),getIntegralErr("VV")).Data(),
	  jmt::format_nevents(getIntegral("Zinvisible"),getIntegralErr("Zinvisible")).Data(), 	  jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm")).Data(), 
	  (int)(getIntegral("data")+0.5), jmt::format_nevents(getIntegral("LM9"),getIntegralErr("LM9")).Data());
  cout<<output<<endl;  

  // -- SL selection --
  selection_ = baselineSL && HTcut && passOther && SRMET && bcut;
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  thisbox=""; thisbox+= isSIG ? "$N_{SIG,SL}$":"$N_{SB,SL}$";
  sprintf(output,"%s & %s%s & %s & %s & %s & %s & %s & %s & %s & %s& %i & %s \\\\",thisbox.Data(),region.btagSelection.Data(),region.owenId.Data(),
	  jmt::format_nevents(getIntegral("PythiaPUQCD"),getIntegralErr("PythiaPUQCD")).Data(), jmt::format_nevents(getIntegral("TTbarJets"),getIntegralErr("TTbarJets")).Data(),
	  jmt::format_nevents(getIntegral("SingleTop"),getIntegralErr("SingleTop")).Data(), jmt::format_nevents(getIntegral("WJets"),getIntegralErr("WJets")).Data(),
	  jmt::format_nevents(getIntegral("ZJets"),getIntegralErr("ZJets")).Data(), jmt::format_nevents(getIntegral("VV"),getIntegralErr("VV")).Data(),
	  jmt::format_nevents(getIntegral("Zinvisible"),getIntegralErr("Zinvisible")).Data(), 	  jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm")).Data(), 
	(int)(getIntegral("data")+0.5), jmt::format_nevents(getIntegral("LM9"),getIntegralErr("LM9")).Data());
  cout<<output<<endl; 
  
  
  // -- LDP selection --
  selection_ = baseline && HTcut && failOther && SRMET && bcut;
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  thisbox=""; thisbox+= isSIG ? "$N_{SIG,LDP}$":"$N_{SB,LDP}$";
  sprintf(output,"%s & %s%s & %s & %s & %s & %s & %s & %s & %s & %s& %i & %s \\\\",thisbox.Data(),region.btagSelection.Data(),region.owenId.Data(),
	  jmt::format_nevents(getIntegral("PythiaPUQCD"),getIntegralErr("PythiaPUQCD")).Data(), jmt::format_nevents(getIntegral("TTbarJets"),getIntegralErr("TTbarJets")).Data(),
	  jmt::format_nevents(getIntegral("SingleTop"),getIntegralErr("SingleTop")).Data(), jmt::format_nevents(getIntegral("WJets"),getIntegralErr("WJets")).Data(),
	  jmt::format_nevents(getIntegral("ZJets"),getIntegralErr("ZJets")).Data(), jmt::format_nevents(getIntegral("VV"),getIntegralErr("VV")).Data(),
	  jmt::format_nevents(getIntegral("Zinvisible"),getIntegralErr("Zinvisible")).Data(), 	  jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm")).Data(), 
	  (int)(getIntegral("data")+0.5), jmt::format_nevents(getIntegral("LM9"),getIntegralErr("LM9")).Data());
  cout<<output<<endl; 
}




void runCountInBoxesBreakdown() {

  setSearchRegions();
  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    countInBoxesBreakdown(searchRegions_[i]) ;
    countInBoxesBreakdown(sbRegions_[i]) ;
  }

}


void countInBoxes(const SearchRegion & region,  ofstream * outfile ) {

  /*
owen asked for total MC event counts in the 6 boxes, output in a particular format
  */

  //obviously this event counting is coded elsewhere in many other forms. but let's just redo it

  //now that i know what owen really wants, i could have done this in the signal systematics code, but
  //now this is written....

  useFlavorHistoryWeights_=false;
  loadSamples(); //needs to come first! this is why this should be turned into a class, with this guy in the ctor

  dodata_=false;
  savePlots_ =false;
  //  setQuiet(true);

  region.Print();

  //for susy scan; harmless for non-susy scan
  resetSamples();
  const TString sample="totalsm";

  TString btagselection = region.btagSelection;
  TCut HTcut=region.htSelection.Data(); 
  TCut SRMET = region.metSelection.Data();
  bool isSIG = region.isSIG;
  TCut bcut = "nbjets>=1";
  TString bweightstring="probge1";
  if (btagselection=="ge2b") {
    bcut="nbjets>=2";
    bweightstring="probge2";
  }
  else if (btagselection=="ge1b") {}
  else if (btagselection=="ge3b") {
    bcut="nbjets>=3";
    bweightstring="probge3"; 
  }
  else if (btagselection=="eq1b") {
    bcut="nbjets==1";
    bweightstring="prob1"; 
  }
  else if (btagselection=="eq2b") {
    bcut="nbjets==2";
    bweightstring="(1-prob1-prob0-probge3)"; 
  }
  else {assert(0);}

  TCut baseline = "cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 &&passCleaning==1";
  TCut baselineSL = "cutPV==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) &&passCleaning==1";

  TCut passOther = "minDeltaPhiN>=4";
  TCut failOther="minDeltaPhiN<4";


  //use all of the corrections that are in the AN table
  currentConfig_=configDescriptions_.getCorrected(); //add  JER bias
  usePUweight_=true;
  useHTeff_=true;
  useMHTeff_=true;
  thebnnMHTeffMode_ = kOff;
  btagSFweight_=bweightstring; //use b tag weight instead

  TString var="HT"; TString xtitle=var;
  int  nbins=1; float low=0; float high=1e9;

/*
Nsig              155
Nsb               244
Nsig_sl           103
Nsb_sl            165
Nsig_ldp          89
Nsb_ldp           393
*/

  // --- nominal selection ---
  selection_ = baseline && HTcut && passOther && SRMET; //no b cut because we're using the weighting
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  (*outfile)<<"N";//<<sampleOwenName_[sample]<<"mc"<<"_";
  if (isSIG) (*outfile)<<"sig            "; 
  else       (*outfile)<<"sb             "; 
  (*outfile)<<getIntegral(sample)<<endl;

  //now the error
//   (*outfile)<<"N"<<sampleOwenName_[sample]<<"mc"<<"_";
//   if (isSIG) (*outfile)<<"sig_err        "; 
//   else       (*outfile)<<"sb_err         "; 
//   (*outfile)<<getIntegralErr(sample)<<endl;

  // -- SL selection --
  selection_ = baselineSL && HTcut && passOther && SRMET; //no b cut because we're using the weighting
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  (*outfile)<<"N";//<<sampleOwenName_[sample]<<"mc"<<"_";
  if (isSIG) (*outfile)<<"sig_sl         "; 
  else       (*outfile)<<"sb_sl          "; 
  (*outfile)<<getIntegral(sample)<<endl;

  //now the error
//   (*outfile)<<"N"<<sampleOwenName_[sample]<<"mc"<<"_";
//   if (isSIG) (*outfile)<<"sig_sl_err     "; 
//   else       (*outfile)<<"sb_sl_err      "; 
//   (*outfile)<<getIntegralErr(sample)<<endl;

  // -- LDP selection --
  selection_ = baseline && HTcut && failOther && SRMET; //no b cut because we're using the weighting
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  (*outfile)<<"N";//<<sampleOwenName_[sample]<<"mc"<<"_";
  if (isSIG) (*outfile)<<"sig_ldp        "; 
  else       (*outfile)<<"sb_ldp         "; 
  (*outfile)<<getIntegral(sample)<<endl;

  //now the error
//   (*outfile)<<"N"<<sampleOwenName_[sample]<<"mc"<<"_";
//   if (isSIG) (*outfile)<<"sig_ldp_err    "; 
//   else       (*outfile)<<"sb_ldp_err     "; 
//   (*outfile)<<getIntegralErr(sample)<<endl;

  

}

void runCountInBoxesMC() {

  setSearchRegions();
  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    char effoutput[500];
    sprintf(effoutput,"mcEventCounts.%s%s.%dinvpb.dat",searchRegions_[i].btagSelection.Data(),searchRegions_[i].owenId.Data(),TMath::Nint(lumiScale_));
    ofstream  textfile(effoutput);

    countInBoxes(searchRegions_[i], &textfile ) ;
    countInBoxes(sbRegions_[i], &textfile ) ;
    textfile.close();
  }

}



std::pair<double,double> ABCD_njetRW(TString prescontrol, TString physcontrol, TString Acutjm, TString Bcutjm, TString Ccutjm, TString Dcutjm, TString lsbbtagsf=""){
  cout << "Running closure with MC jet multiplicity reweighted to match data." << endl;
  
  //  bool quietornot = quiet_;
  //  setQuiet(false);

  TString histfilename = "dummy.root";
  TFile fh(histfilename,"RECREATE");//will delete old root file if present 
  fh.Close(); //going to open only when needed 

  double physlumiscale = lumiScale_;
  double preslumiscale = 30.15471;//hardcoded for now
  TString btagSFweight_nom = btagSFweight_;

  //going to use the same binning for the physics and prescaled samples
  //this isn't necessary, but it is simpler
  int jmnbins = 4;
  float jmbins[] = {2.5,3.5,4.5,5.5,6.5};


  //PHYS jet multiplicity scale factors
  ////////////////////////////////////////////////////////////
  physcontrol+= " && minDeltaPhiN<4 && MET<250"; //LDP SB

  lumiScale_ = physlumiscale - preslumiscale;
  if(useScaleFactors_){ 
    btagSFweight_ = "prob0";
  }
  else { physcontrol+= " && nbjetsCSVM==0"; }

  selection_ = physcontrol;
  drawSimple("njets30",jmnbins,jmbins,"dummy.root", "","data");
  TH1D* hJMphysicsData = (TH1D*)hinteractive->Clone("hJMphysicsData");
  
  if(useScaleFactors_){
    float eff_SB_ldp_MHT = eff_SB_ldp_MHT_; 
    //float eff_SB_ldp_MHT_err[2] = {eff_SB_ldp_MHT_err_[0], eff_SB_ldp_MHT_err_[1]};
    hJMphysicsData->Scale(1/eff_SB_ldp_MHT);//for now, ignore the error on the efficiency 
  }

  drawSimple("njets30",jmnbins,jmbins,"dummy.root", "","PythiaPUQCD");
  TH1D* hJMphysicsQCD = (TH1D*)hinteractive->Clone("hJMphysicsQCD");
  drawPlots("njets30",jmnbins,jmbins, "", "", "deleteme");
  TH1D* hJMphysicsNonQCD = (TH1D*)totalnonqcd->Clone("hJMphysicsNonQCD");


  //LSB jet multiplicity scale factors
  ///////////////////////////////////////////
  lumiScale_ = preslumiscale;
  //original design of the analysis uses >=1b for the prescaled sample here. 
  //but in the case of the b-tagged LSB, we want to use the =0b sample to avoid correlations
  if(useScaleFactors_){
    if (btaggedLSB_) btagSFweight_ = "prob0";
    else    btagSFweight_ = "probge1";
  }
  else{ 
    if (btaggedLSB_) prescontrol+= " && nbjetsCSVM==0"; 
    else    prescontrol+= " && nbjetsCSVM>=1"; 
  }

  selection_= prescontrol;
  drawSimple("njets30",jmnbins,jmbins,"dummy.root", "","data");
  TH1D* hJMprescaleData = (TH1D*)hinteractive->Clone("hJMprescaleData");
  drawSimple("njets30",jmnbins,jmbins,"dummy.root", "","PythiaPUQCD");
  TH1D* hJMprescaleQCD = (TH1D*)hinteractive->Clone("hJMprescaleQCD");


  //Now ABCD Regions

  //jmt -- seems like this is the event counting that is the actual MC-based closure test
  //for usebtaggedLSB_, we need to apply b-tagging, I think
  //////////////////////////////////////////
  lumiScale_ = physlumiscale;
  if(useScaleFactors_) {
    if (btaggedLSB_)   btagSFweight_ = (lsbbtagsf!="") ?  lsbbtagsf : btagSFweight_nom;
    else btagSFweight_ = "prob0";
  }

  selection_ = Acutjm;

  drawSimple("njets30",jmnbins,jmbins,"dummy.root", "","PythiaPUQCD");
  TH1D* hJMAqcd = (TH1D*)hinteractive->Clone("hJMAqcd");
  
  selection_ = Bcutjm;
  drawSimple("njets30",jmnbins,jmbins,"dummy.root", "","PythiaPUQCD");
  TH1D* hJMBqcd = (TH1D*)hinteractive->Clone("hJMBqcd");
  

  btagSFweight_ = btagSFweight_nom;

  selection_ = Dcutjm;
  drawSimple("njets30",jmnbins,jmbins,"dummy.root", "","PythiaPUQCD");
  TH1D* hJMDqcd = (TH1D*)hinteractive->Clone("hJMDqcd");

  selection_ = Ccutjm;
  drawSimple("njets30",jmnbins,jmbins,"dummy.root", "","PythiaPUQCD");
  TH1D* hJMCqcd = (TH1D*)hinteractive->Clone("hJMCqcd");
  
  //Now we have all the input we need to do the reweighting..


  //weights////////////////////////////////////////////////
  for(int k=1; k<=hJMphysicsData->GetNbinsX(); k++){
    assert(hJMphysicsData->GetBinContent(k)>0);
    assert(hJMphysicsQCD->GetBinContent(k)>0);
    assert(hJMprescaleData->GetBinContent(k)>0);
    assert(hJMprescaleQCD->GetBinContent(k)>0);
    assert((hJMphysicsData->GetBinContent(k)- hJMphysicsNonQCD->GetBinContent(k))>0);//need to think about this more
  }
  TH1D* hphysW = (TH1D*)hJMphysicsData->Clone("hphysW");
  hphysW->Add(hJMphysicsNonQCD,-1);
  hphysW->Divide(hJMphysicsQCD);
  TH1D* hpresW = (TH1D*)hJMprescaleData->Clone("hpresW"); //assume no contamination for LSB
  hpresW->Divide(hJMprescaleQCD);
      
  TFile fh1(histfilename, "UPDATE");
  hphysW->Write();
  hpresW->Write();
  fh1.Close();

  //reweighted boxes
  TH1D* hJMAqcdRW = (TH1D*)hJMAqcd->Clone("hJMAqcdRW");
  TH1D* hJMBqcdRW = (TH1D*)hJMBqcd->Clone("hJMBqcdRW");
  TH1D* hJMCqcdRW = (TH1D*)hJMCqcd->Clone("hJMCqcdRW"); 
  TH1D* hJMDqcdRW = (TH1D*)hJMDqcd->Clone("hJMDqcdRW");
  hJMAqcdRW->Multiply(hpresW);
  hJMBqcdRW->Multiply(hpresW);
  hJMCqcdRW->Multiply(hphysW);
  hJMDqcdRW->Multiply(hphysW);

  double Arw, Brw, Crw, Drw;
  Arw = hJMAqcdRW->Integral();
  Brw = hJMBqcdRW->Integral();
  Crw = hJMCqcdRW->Integral();
  Drw = hJMDqcdRW->Integral();


  ////////////////////////////////////////////////////////////stat//
  double rwErr2 = 0;
  for(int k=1; k<=hJMphysicsData->GetNbinsX(); k++){
    double QAk = hJMAqcdRW->GetBinContent(k);
    double QBk = hJMBqcdRW->GetBinContent(k);
    double QCk = hJMCqcdRW->GetBinContent(k);
    double QDk = hJMDqcdRW->GetBinContent(k);
    double Qk = hJMphysicsQCD->GetBinContent(k);
    double nk = hJMphysicsNonQCD->GetBinContent(k);
    double Rk = hJMprescaleData->GetBinContent(k);
    double rk = hJMprescaleQCD->GetBinContent(k);
    double Pk = hJMphysicsData->GetBinContent(k);

    double QAkerr = hJMAqcdRW->GetBinError(k);
    double QBkerr = hJMBqcdRW->GetBinError(k);
    double QCkerr = hJMCqcdRW->GetBinError(k);
    double QDkerr = hJMDqcdRW->GetBinError(k);
    double Qkerr = hJMphysicsQCD->GetBinError(k);
    double nkerr = hJMphysicsNonQCD->GetBinError(k);
    double Rkerr = hJMprescaleData->GetBinError(k);
    double rkerr = hJMprescaleQCD->GetBinError(k);
    double Pkerr = hJMphysicsData->GetBinError(k);

    double QAterm=0,QBterm=0,QCterm=0,QDterm=0,Qterm=0, nterm=0, Rterm=0, Pterm=0, rterm=0;
    double X_toppQ=0, X_botpQ=0, X_toppn=0, X_botpn=0, X_toppR=0, X_botpR=0, X_toppP=0, X_botpP=0, X_toppr=0, X_botpr=0;
    double X_top = Crw*Arw;
    double X_bot = Drw*Brw;
       
    QAterm = Crw/Drw/Brw*Rk/rk;
    QBterm = -Crw*Arw/Drw/Brw/Brw*Rk/rk;
    QCterm = Arw/Drw/Brw*(Pk-nk)/Qk;
    QDterm = -Crw*Arw/Brw/Drw/Drw*(Pk-nk)/Qk;

    X_toppQ = -Arw*QCk*(Pk-nk)/Qk/Qk;
    X_botpQ = -Brw*QDk*(Pk-nk)/Qk/Qk;

    X_toppn = -Arw*QCk/Qk;
    X_botpn = -Brw*QDk/Qk;
    
    X_toppR = Crw*QAk/rk;
    X_botpR = Drw*QBk/rk;

    X_toppr = -Crw*QAk*Rk/rk/rk;
    X_botpr = -Drw*QBk*Rk/rk/rk;

    X_toppP = Arw*QCk/Qk;
    X_botpP = Brw*QDk/Qk;

    //quotient rule!
    Qterm = (X_bot*X_toppQ - X_top*X_botpQ)/X_bot/X_bot;
    nterm = (X_bot*X_toppn - X_top*X_botpn)/X_bot/X_bot;
    Rterm = (X_bot*X_toppR - X_top*X_botpR)/X_bot/X_bot;
    rterm = (X_bot*X_toppr - X_top*X_botpr)/X_bot/X_bot;
    Pterm = (X_bot*X_toppP - X_top*X_botpP)/X_bot/X_bot;
    
    cout << "njetrw: QA: " << QAterm*QAterm*QAkerr*QAkerr << endl;
    cout << "njetrw: QB: " << QBterm*QBterm*QBkerr*QBkerr << endl;
    cout << "njetrw: QC: " << QCterm*QCterm*QCkerr*QCkerr << endl;
    cout << "njetrw: QD: " << QDterm*QDterm*QDkerr*QDkerr << endl;
    cout << "njetrw: Q: " << Qterm*Qterm*Qkerr*Qkerr << endl;
    cout << "njetrw: n: " << nterm*nterm*nkerr*nkerr << endl;
    cout << "njetrw: R: " << Rterm*Rterm*Rkerr*Rkerr << endl;
    cout << "njetrw: r: " << rterm*rterm*rkerr*rkerr << endl;
    cout << "njetrw: P: " << Pterm*Pterm*Pkerr*Pkerr << endl;

    rwErr2 += QAterm*QAterm*QAkerr*QAkerr + QBterm*QBterm*QBkerr*QBkerr + QCterm*QCterm*QCkerr*QCkerr + QDterm*QDterm*QDkerr*QDkerr;
    rwErr2 += Qterm*Qterm*Qkerr*Qkerr + nterm*nterm*nkerr*nkerr + Rterm*Rterm*Rkerr*Rkerr + rterm*rterm*rkerr*rkerr + Pterm*Pterm*Pkerr*Pkerr;
  }
  //  setQuiet(quietornot);
  //////////////////////////////////////////////////////////////////
  cout << "njetrw: ABCD: " << Arw << " " << Brw << " " << Crw << " " << Drw << endl; 
  cout << "njetrw: BD/A (rw): " << Brw*Drw/Arw << endl;
  cout << "njetrw: C (rw): " << Crw << endl;
  double closure = 1.0-Crw*Arw/Drw/Brw;
  double closurestat = sqrt(rwErr2);
  cout << "njetrw: closure (rw): " << 100.*closure << " +- " << 100.*closurestat << endl;
  cout << "njetrw: Table: $" << Brw << "$ & $" << Arw << "$ & $" << Drw << "$ & $" << Brw*Drw/Arw << "$ & $" << Crw << "$ & $" <<  100.*closure << " \\pm " << 100.*closurestat << "$ \\\\" << endl;
  
  return make_pair( 100*sqrt(pow(closure,2) + pow(closurestat,2)), 0); //estimate in denominator
}

TString getLSBbsel(const SearchRegion & r) {
  TString lsbbsel = "==0"; //default is use 0b LSB
  if (btaggedLSB_) { //in case we want btagged LSB
    if (r.btagSelection.Contains("3")) lsbbsel = "ge2b"; //use >=2b in place of >=3b
    else lsbbsel = r.btagSelection;
  }

  return lsbbsel;
}

std::pair<double,std::vector<double> > anotherABCD( const SearchRegion & region, bool datamode=false, float subscale=1,float SBshift=0, TString LSBbsel="==0", float PVCorFactor = 0, bool doNjetRW = false, bool forTTbarEstimate=false) {
  //kind of awful, but we'll return the estimate for datamode=true but the closure test bias for datamode=false
  
  doOverflowAddition(true);
  
  TString btagselection = region.btagSelection;
  
  //for user convenience, to allow different forms of input
  if (LSBbsel == "eq0b") LSBbsel = "==0";
  else  if (LSBbsel == "ge1b") LSBbsel = ">=1";
  else  if (LSBbsel == "ge2b") LSBbsel = ">=2";
  else  if (LSBbsel == "ge3b") LSBbsel = ">=3";
  else  if (LSBbsel == "eq1b") LSBbsel = "==1";
  else  if (LSBbsel == "eq2b") LSBbsel = "==2";
  //now just logic checks
  else  if (LSBbsel == "==0") {}
  else  if (LSBbsel == "==1") {}
  else  if (LSBbsel == "==2") {}
  else  if (LSBbsel == ">=1") {}
  else  if (LSBbsel == ">=2") {}
  else  if (LSBbsel == ">=3") {}
  else assert(0);

  TString owenKey = btagselection;
  owenKey += region.owenId;
  OwenData * myOwen = &(owenMap_[owenKey]);
  const  bool isSIG = region.isSIG;
  
  setStackMode(false);
  doData(true);
  setQuiet(true);
  
  useFlavorHistoryWeights_=false;
  
  loadSamples(); //needs to come first! this is why this should be turned into a class, with this guy in the ctor
  setColorScheme("nostack");
  //clearSamples();
  //addSample("PythiaPUQCD");
  
  TString sampleOfInterest = "PythiaPUQCD";
  if (datamode) {
    sampleOfInterest = "data";
    clearSamples();
    addSample("TTbarJets");
    addSample("WJets");
    addSample("ZJets");
    addSample("Zinvisible");
    addSample("VV");
    addSample("SingleTop");
  }
  else {
    resetSamples();
  }
  
  savePlots_=false;
  
  setLogY(false);
  TString var,xtitle;
  int nbins;
  float low,high;
  
  // --  count events
  TCut triggerCutLSB = "1";
  TCut triggerCutLSBInverted = "1";
  TCut triggerCut = "1";
  triggerCutLSB = "pass_utilityHLT==1";
  triggerCutLSBInverted = "(pass_utilityHLT==0 || isRealData==0)";
  triggerCut = "cutTrigger==1";

  //don: hard-code MHT efficiencies
  //offline cuts: HT>400, 0L, mindphin>4(<4), 0b
  float eff_MHT = 1, eff_MHT_err[2];//first entry is + error, second entry is - error
  float eff_ldp_MHT = 1, eff_ldp_MHT_err[2];
  TString metselection = region.metSelection;
  
  TCut ge1b = "nbjetsCSVM>=1";
  btagSFweight_="1";
  TString btagSFweight="1"; //we're going to have to switch this one in and out of the global var
  if (useScaleFactors_) {
    usePUweight_=true;
    useHTeff_=true;
    useMHTeff_=false;
    thebnnMHTeffMode_ = kOff;//dealt with later
    currentConfig_=configDescriptions_.getCorrected(); //add JERbias

    if(metselection=="MET>=200&&MET<250" || metselection=="MET>=150&&MET<200") {
      //      cout<<"Using SB efficiency corrections for QCD"<<endl;
      eff_MHT     = eff_SB_MHT_;      eff_MHT_err[0] = eff_SB_MHT_err_[0];         eff_MHT_err[1] = eff_SB_MHT_err_[1];
      eff_ldp_MHT = eff_SB_ldp_MHT_ ; eff_ldp_MHT_err[0] = eff_SB_ldp_MHT_err_[0]; eff_ldp_MHT_err[1] = eff_SB_ldp_MHT_err_[1];
    }
    else{
      //      cout<<"Using SIG efficiency corrections for QCD"<<endl;
      eff_MHT     = eff_SIG_MHT_;      eff_MHT_err[0] = eff_SIG_MHT_err_[0];         eff_MHT_err[1] = eff_SIG_MHT_err_[1];
      eff_ldp_MHT = eff_SIG_ldp_MHT_ ; eff_ldp_MHT_err[0] = eff_SIG_ldp_MHT_err_[0]; eff_ldp_MHT_err[1] = eff_SIG_ldp_MHT_err_[1];
    }
    
    ge1b="1";
    if (btagselection=="ge2b") {
      btagSFweight="probge2";
    }
    else if (btagselection=="ge1b") {
      btagSFweight="probge1";
    }
    else if (btagselection=="ge3b") {
      btagSFweight="probge3";
    }
    else if (btagselection=="eq1b") {
      btagSFweight="prob1";
    }
    else if (btagselection=="eq2b") {
      btagSFweight="(1-prob1-prob0-probge3)";
    }
    else {assert(0);}
  }
  else {
    usePUweight_=false;
    useHTeff_=false;
    useMHTeff_=false;
    btagSFweight_="1";
    currentConfig_=configDescriptions_.getDefault(); //completely raw MC
    
    if (btagselection=="ge2b") {
      ge1b="nbjetsCSVM>=2";
    }
    else if (btagselection=="ge1b") {}
    else if (btagselection=="ge3b") {
      ge1b="nbjetsCSVM>=3";
    }
    else if (btagselection=="eq1b") {
      ge1b="nbjetsCSVM==1";
    }
    else if (btagselection=="eq2b") {
      ge1b="nbjetsCSVM==2";
    }
    else {assert(0);}
  }
  
  TCut HTcut=region.htSelection.Data(); 
  TCut SRMET = region.metSelection.Data();
  
  SRMET = SRMET && ge1b && triggerCut; //apply physics trigger and physics b cut in high MET region
  
  //  TCut baseline = "cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1";
  TCut baseline = TCut("cutPV==1 && cut3Jets==1") && getLeptonVetoCut();
  baseline = baseline&&HTcut;
  TCut cleaning = "weight<1000 && passCleaning==1";
  
  
  //LSB b-tagging is independent of search region
  TString LSBbtagSFweight="";
  char cutstring0[100];
  sprintf(cutstring0,"nbjetsCSVM%s",LSBbsel.Data());
  TCut LSBbtag = TCut(cutstring0); 
  if(useScaleFactors_){
    LSBbtag = "1";
    if(LSBbsel=="==0"){ LSBbtagSFweight = "prob0"; }
    else if(LSBbsel=="==1"){ LSBbtagSFweight = "prob1"; }
    else if(LSBbsel=="==2"){ LSBbtagSFweight = "(1-prob1-prob0-probge3)"; }
    else if(LSBbsel==">=1"){ LSBbtagSFweight = "probge1"; }
    else if(LSBbsel==">=2"){ LSBbtagSFweight = "probge2"; }
    else if(LSBbsel==">=3"){ LSBbtagSFweight = "probge3"; }
    else{ assert(0);}
  }
  //if  useScaleFactors_, LSBbtag="1" and LSBbtagSFweight is "prob0" or "probge1"
  //if !useScaleFactors_, LSBbtag is the real cut, and LSBbtagSFweight is ""

  char cutstring1[100];
  sprintf(cutstring1,"MET>= %.0f && MET < %.0f", 50.0+SBshift,100.0);
  TCut SBMET = TCut(cutstring1)&&triggerCutLSB&&LSBbtag;


  TCut dpcut = "1";//"minDeltaPhiN>=4";
  //  TCut passOther = "deltaPhiMPTcaloMET<2";
  //  TCut failOther = "deltaPhiMPTcaloMET>=2";
  TCut passOther = "minDeltaPhiN>=4";
  TCut failOther = "minDeltaPhiN<4";

  
  var="HT"; xtitle=var;
  nbins=10; low=0; high=5000;
  double A,B,D,SIG,Aerr,Berr,Derr,SIGerr;
  double D_bnnPlus=0, D_bnnMinus=0;
  //double DerrMinus;
  double RLSB_RW = 0;
  double dRLSB_RW = 0;
  
  //used for jet multiplicity reweighting closure test
  TString Acutjm, Bcutjm, Dcutjm, SIGcutjm;
  
  //baselines used in lsb and njet reweighting -- notice no mdpN or btag cut!
  //  TCut physcontrol =  TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>200") && triggerCut && triggerCutLSBInverted && HTcut && cleaning;
  TCut physcontrol =  TCut("cutPV==1 && cut3Jets==1 && MET>200") && triggerCut && triggerCutLSBInverted && HTcut && cleaning && getLeptonVetoCut();
  TString physcontrolstring; physcontrolstring += physcontrol;
  TCut prescontrol =  TCut("cutPV==1 && cut3Jets==1 && MET>=50 && MET<100") && triggerCutLSB && HTcut && cleaning && getLeptonVetoCut();
  TString prescontrolstring; prescontrolstring += prescontrol;

  //FOR LSB
  btagSFweight_ =  LSBbtagSFweight;
    
  if(datamode &&  reweightLSBdata_) {

    float *pvbins;
    int pvnbins=11;
    float pvbins_1[]= {0.5,2.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,12.5,14.5,16.5};
    float pvbins_2[]={0.5,3.5,5.5,8.5,11.5,16.5};
    float pvbins_3[]={0.5,5.5,10.5,16.5};

    int choice = 1; //this logic is kinda convoluted, but it (barely) works
    if ( LSBbsel.Contains("2") ) choice =2;
    else if (LSBbsel.Contains("3") ) choice =3;
    else if (region.owenId.Contains("TightHT")) choice =2;
    else if (region.owenId.Contains("TightMET")) choice=3;

    //adjust to coarser binning for higher numbers of b tags
    if (choice==2 ) pvnbins=5;
    else if (choice==3 )   pvnbins=3;
    pvbins = new float[pvnbins+1];
    for (int jj=0; jj<pvnbins+1; jj++)  {
      if (choice==2 ) pvbins[jj] = pvbins_2[jj];
      else if (choice==3 ) pvbins[jj] = pvbins_3[jj];
      else if (choice==1) pvbins[jj] = pvbins_1[jj];
      else assert(0);
    }

    //physics triggers control sample -- physics triggers, 0b
    TCut antibcut = "nbjetsCSVM==0"; 
    btagSFweight_ = "1"; //this will always be 0b data
    selection_ = physcontrol && antibcut && passOther; //add mdpN cut to keep sample indepdenent of njet rw
    drawSimple("nGoodPV",pvnbins,pvbins,"dummy", "","data");
    TH1D* hPVphysics = (TH1D*)hinteractive->Clone("hPVphysics");
    btagSFweight_ =  LSBbtagSFweight; //restore nominal LSB b tag cut

    //LSB unweighted
    selection_ = baseline && SBMET && cleaning;
    drawSimple("nGoodPV",pvnbins,pvbins,"dummy", "","data");
    TH1D* hPVprescale = (TH1D*)hinteractive->Clone("hPVprescale");
    
    //calculate weights (physics shape with prescale integral divided by prescale)
    TH1D* hPVprescale_RW = (TH1D*)hPVphysics->Clone("hPVprescale_RW");
    hPVprescale_RW->Scale(hPVprescale->Integral()/hPVphysics->Integral());
    TH1D* hPV_W = (TH1D*)hPVphysics->Clone("hPV_W");
    hPV_W->Reset();
    hPV_W->Divide(hPVprescale_RW,hPVprescale); 
    
    //LSB pass mdp unweighted
    selection_ = baseline && SBMET && passOther && cleaning;
    drawSimple("nGoodPV",pvnbins,pvbins,"dummy", "","data");
    TH1D* hPVprescalePass = (TH1D*)hinteractive->Clone("hPVprescalePass");
    
    //LSB fail mdp unweighted
    selection_ = baseline && SBMET && failOther && cleaning;
    drawSimple("nGoodPV",pvnbins,pvbins,"dummy", "","data");
    TH1D* hPVprescaleFail = (TH1D*)hinteractive->Clone("hPVprescaleFail");
    
    //check for zero entries
    for(int j=1; j<=hPVphysics->GetNbinsX();j++){
      assert(hPVphysics->GetBinContent(j)>0);
      //      cout<<hPVprescalePass->GetBinLowEdge(j)<<" "<<hPVprescalePass->GetBinContent(j)<<endl;
      assert(hPVprescalePass->GetBinContent(j)>0);
      assert(hPVprescaleFail->GetBinContent(j)>0);
    }
    
    //weighted LSB pass and fail mdp
    TH1D* hPVprescalePass_RW = (TH1D*)hPVprescalePass->Clone("hPVprescalePass_RW");
    TH1D* hPVprescaleFail_RW = (TH1D*)hPVprescaleFail->Clone("hPVprescaleFail_RW");
    hPVprescalePass_RW->Multiply(hPV_W);
    hPVprescaleFail_RW->Multiply(hPV_W);

    delete [] pvbins;
    
    ///////////////////////////////////////////////////////////////////////////////stat///////
    double dR2 = 0;
    double myA=0, myB=0;
    for(int k=1; k<=hPVphysics->GetNbinsX(); k++){
      double dNk = hPVphysics->GetBinError(k);
      double dPk = hPVprescalePass->GetBinError(k);
      double dFk = hPVprescaleFail->GetBinError(k);
      double Nk = hPVphysics->GetBinContent(k);
      double Pk = hPVprescalePass->GetBinContent(k);
      double Fk = hPVprescaleFail->GetBinContent(k);
      double Ntot = hPVphysics->Integral();
      double PFtot = hPVprescalePass->Integral()+hPVprescaleFail->Integral();
      
      //pX means derivative w.r.t. X
      double ApN=0, BpN=0;
      double ApP=0, BpP=0;
      double ApF=0, BpF=0;
      double Nterm=0, Pterm=0, Fterm=0;
      for(int i=1; i<=hPVphysics->GetNbinsX(); i++){
	double Ni = hPVphysics->GetBinContent(i);
	double Pi = hPVprescalePass->GetBinContent(i);
	double Fi = hPVprescaleFail->GetBinContent(i);
	
	myB += Ni/Ntot*PFtot/(Pi+Fi)*Pi;
	myA += Ni/Ntot*PFtot/(Pi+Fi)*Fi;

	//part without k
	if(i==k) continue;
	BpN += (-1.)*PFtot/(Pi+Fi)*Pi*Ni/(Ntot*Ntot);
	ApN += (-1.)*PFtot/(Pi+Fi)*Fi*Ni/(Ntot*Ntot);
	BpP += Ni/Ntot*Pi/(Pi+Fi);
	ApP += Ni/Ntot*Fi/(Pi+Fi);
	BpF += Ni/Ntot*Pi/(Pi+Fi);
	ApF += Ni/Ntot*Fi/(Pi+Fi);
      }//i loop
      
      //part with k
      BpN += PFtot/(Pk+Fk)*Pk*(Ntot-Nk)/(Ntot*Ntot);
      ApN += PFtot/(Pk+Fk)*Fk*(Ntot-Nk)/(Ntot*Ntot);
      BpP += Nk/Ntot/(Pk+Fk)/(Pk+Fk)*((Pk+Fk)*(PFtot+Pk)-PFtot*Pk);
      ApP += Nk/Ntot*Fk/(Pk+Fk)/(Pk+Fk)*(Pk+Fk-PFtot);
      BpF += Nk/Ntot*Pk/(Pk+Fk)/(Pk+Fk)*(Pk+Fk-PFtot);
      ApF += Nk/Ntot/(Pk+Fk)/(Pk+Fk)*((Pk+Fk)*(PFtot+Fk)-PFtot*Fk);      
      
      //quotient rule!
      Nterm = (myA*BpN-myB*ApN)/(myA*myA);
      Pterm = (myA*BpP-myB*ApP)/(myA*myA);
      Fterm = (myA*BpF-myB*ApF)/(myA*myA);
      
      dR2 += Nterm*Nterm*dNk*dNk + Pterm*Pterm*dPk*dPk + Fterm*Fterm*dFk*dFk;
      
    }//k loop
    RLSB_RW = myB/myA; //matches hPVprescalePass_RW->Integral()/hPVprescaleFail_RW->Integral()
    dRLSB_RW = sqrt(dR2);
    ///////////////////////////////////////////////////////////////////////////////////////////
    //unweighted
    B = hPVprescalePass->Integral();
    Berr = jmt::errOnIntegral(hPVprescalePass); 
    A = hPVprescaleFail->Integral();
    Aerr = jmt::errOnIntegral(hPVprescaleFail);
  }
  else{
    //A   -- aka 50 - 100 and high MPT,MET
    selection_ = baseline && cleaning && dpcut  && SBMET && failOther; //auto cast to TString seems to work
    Acutjm = selection_;
    drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_A");
    A=getIntegral(sampleOfInterest);
    Aerr=getIntegralErr(sampleOfInterest);
    //B
    selection_ = baseline && cleaning && dpcut  && SBMET && passOther; //auto cast to TString seems to work
    Bcutjm = selection_;
    drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_B");
    B=getIntegral(sampleOfInterest);
    Berr=getIntegralErr(sampleOfInterest);
  }
  
  //Now that we've moved on from the LSB, use search region b-tag requirement
  btagSFweight_ = btagSFweight;
  
  //D
  selection_ = baseline && cleaning && dpcut  && SRMET && failOther; //auto cast to TString seems to work
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_D");
  Dcutjm = selection_;
  D=getIntegral(sampleOfInterest);
  Derr=getIntegralErr(sampleOfInterest);

  double D_temp = D, Derr_temp = Derr;
  if(useScaleFactors_ && datamode){

    if(useBNNEffCurves_){
      thebnnMHTeffMode_ = kOn;
      drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_D");
      D=getIntegral(sampleOfInterest);
      Derr=getIntegralErr(sampleOfInterest);      

      thebnnMHTeffMode_ = kOnPlus;
      drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_D");
      D_bnnPlus=getIntegral(sampleOfInterest);
      thebnnMHTeffMode_ = kOnMinus;
      drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_D");
      D_bnnMinus=getIntegral(sampleOfInterest);

      thebnnMHTeffMode_ = kOff;//turn it back off
    }
    else{
      D = D_temp/eff_ldp_MHT;
      //Derr = jmt::errAoverB(D_temp,Derr_temp,eff_ldp_MHT,eff_ldp_MHT_err[0]);
      //DerrMinus = jmt::errAoverB(D_temp,Derr_temp,eff_ldp_MHT,eff_ldp_MHT_err[1]);
      
      //Derr = Derr > DerrMinus ? Derr: DerrMinus;
    }
  }
  
  double Dsub = 0,Dsuberr=0;
  if (datamode) {
    Dsub = subscale*getIntegral("totalsm");
    Dsuberr = subscale*getIntegralErr("totalsm");
    
    //special stuff for owen
    //    myOwen->Nlsb_0b_ldp = A;
    //    myOwen->Nlsb_0b     = B;
    if (isSIG) myOwen->Nsig_ldp = D_temp;
    else       myOwen->Nsb_ldp = D_temp;
    
    if (isSIG) {
      myOwen->Nttbarmc_sig_ldp = getIntegral("TTbarJets");
      myOwen->Nsingletopmc_sig_ldp = getIntegral("SingleTop");
      //      double lsfw = getIntegral("WJets") > 0 ?  pow(getIntegralErr("WJets"),2)/getIntegral("WJets"): -1;
      //      double lsfzj = getIntegral("ZJets") > 0 ?  pow(getIntegralErr("ZJets"),2)/getIntegral("ZJets"): -1;
      //      double lsfz = getIntegral("Zinvisible") > 0 ?  pow(getIntegralErr("Zinvisible"),2)/getIntegral("Zinvisible"): -1;
      
      //      if (lsfw>0) myOwen->lsf_WJmc=lsfw;
      //      if (lsfz>0) myOwen->lsf_Znnmc=lsfz;
      //      if (lsfzj>0) myOwen->lsf_Zjmc=lsfzj;
      
      myOwen->NWJmc_sig_ldp = getIntegral("WJets");// /lsfw;
      myOwen->NZnnmc_sig_ldp = getIntegral("Zinvisible");// /lsfz;
      myOwen->NZjmc_sig_ldp = getIntegral("ZJets");// /lsfzj;
      
    }
    else {
      myOwen->Nttbarmc_sb_ldp = getIntegral("TTbarJets");
      myOwen->Nsingletopmc_sb_ldp = getIntegral("SingleTop");
      
      //      double lsfw = getIntegral("WJets") > 0 ?  pow(getIntegralErr("WJets"),2)/getIntegral("WJets"): -1;
      //      double lsfz = getIntegral("Zinvisible") > 0 ?  pow(getIntegralErr("Zinvisible"),2)/getIntegral("Zinvisible"): -1;
      //      double lsfzj = getIntegral("ZJets") > 0 ?  pow(getIntegralErr("ZJets"),2)/getIntegral("ZJets"): -1;
      
      //      if (lsfw>0) myOwen->lsf_WJmc=lsfw;
      //      if (lsfz>0) myOwen->lsf_Znnmc=lsfz;
      //      if (lsfzj>0) myOwen->lsf_Zjmc=lsfzj;
      
      myOwen->NWJmc_sb_ldp = getIntegral("WJets");// /lsfw;
      myOwen->NZnnmc_sb_ldp = getIntegral("Zinvisible");// /lsfz;
      myOwen->NZjmc_sb_ldp = getIntegral("ZJets");// /lsfzj;
    }
    //end of special stuff for owen
  }
  if (Dsub>D) {Dsub=D-0.00001; cout<<"Subtraction in D is as big as D!"<<endl;}
  
  //SIG
  selection_ = baseline && cleaning && dpcut  && SRMET && passOther; //auto cast to TString seems to work
  SIGcutjm = selection_;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_C");
  SIG=getIntegral(sampleOfInterest);
  SIGerr=getIntegralErr(sampleOfInterest);
  
  if (datamode) { //for owen 
    if (isSIG)  myOwen->Nsig = SIG;
    else        myOwen->Nsb  = SIG;
  }

  
  if(doNjetRW){ return ABCD_njetRW(prescontrolstring, physcontrolstring, Acutjm, Bcutjm, SIGcutjm, Dcutjm,LSBbtagSFweight);}
  
  //now calculate (B/A)*D=R*D
  double RLSB_UW = B/A; //unweighted
  double RLSB_UWerr = jmt::errAoverB(B,Berr,A,Aerr);
  
  double myR, myRerr;
  if(datamode &&  reweightLSBdata_){
    myR = PVCorFactor*(RLSB_RW-RLSB_UW)+RLSB_RW; //PVCorFactor=0 for reweighted, +/-1 for +/-100% correction
    myRerr = dRLSB_RW;
  }
  else{
    myR = RLSB_UW;
    myRerr = RLSB_UWerr; 
  }
  if (datamode) { //jan13 new. is this ok?
    myOwen-> Rlsb_passfail= myR;
    myOwen-> Rlsb_passfail_err= myRerr;
  }

  double estimate = myR*(D-Dsub);
  double estimateerr = jmt::errAtimesB(myR, myRerr, D-Dsub, sqrt(Derr*Derr+Dsuberr*Dsuberr));

  double estimate_temp = estimate;
  //double estimateerr_temp = estimateerr;
  double estimateerrTrigPlus  = 0;
  double estimateerrTrigMinus = 0;

  if(useScaleFactors_ && datamode){
    //do NOT scale by efficiency if this is a SB region and if this is for the ttbar estimate
    if( !(forTTbarEstimate) ){
      //double estimateerrMinus=0;
      estimate = eff_MHT*estimate_temp;
      //estimateerr = jmt::errAtimesB(estimate_temp, estimateerr_temp, eff_MHT, eff_MHT_err[0]);
      //estimateerrMinus = jmt::errAtimesB(estimate_temp, estimateerr_temp, eff_MHT, eff_MHT_err[1]);
      //std::cout << "estimate: "<< estimate << "+" << estimateerr << "- " <<estimateerrMinus << std::endl;
      
      //use the larger of the two as the error in the estimate
      //estimateerr = estimateerr > estimateerrMinus ? estimateerr: estimateerrMinus;
    }

    //vary the trigger efficiencies by their statistical error to get trigger error on estimate
    //first vary eff_ldp_MHT

    double estimate_up, estimate_down;

    if(useBNNEffCurves_){
      if( !(forTTbarEstimate) ) estimate_up = eff_MHT*myR*(D_bnnPlus - Dsub);
      else estimate_up = myR*(D_bnnPlus - Dsub);
      if( !(forTTbarEstimate) ) estimate_down = eff_MHT*myR*(D_bnnMinus - Dsub);
      else estimate_down = myR*(D_bnnMinus - Dsub);
    }
    else{
      float eff_ldp_MHT_up   = eff_ldp_MHT + eff_ldp_MHT_err[0];
      D = D_temp/eff_ldp_MHT_up;
      if( !(forTTbarEstimate) ) estimate_up = eff_MHT*myR*(D-Dsub);
      else estimate_up = myR*(D-Dsub);
      
      float eff_ldp_MHT_down = eff_ldp_MHT - eff_ldp_MHT_err[1];
      D = D_temp/eff_ldp_MHT_down;
      if( !(forTTbarEstimate) ) estimate_down = eff_MHT*myR*(D-Dsub);
      else estimate_down = myR*(D-Dsub);
    }

    double delta_down = (estimate_up - estimate);    
    double delta_up = (estimate_down - estimate);
    //std::cout << "TRIG: + " << delta_up << " - " << delta_down << std::endl;

    //second vary eff_MHT (except when evaluating SB number for ttbar estimate)
    float eff_MHT_up   = eff_MHT + eff_MHT_err[0];
    estimate_up = eff_MHT_up*estimate_temp;
    double delta_up2;
    if( !(forTTbarEstimate) ) delta_up2 = (estimate_up - estimate);
    else delta_up2 = 0;

    float eff_MHT_down = eff_MHT - eff_MHT_err[1];
    estimate_down = eff_MHT_down*estimate_temp;
    double delta_down2;
    if( !(forTTbarEstimate) ) delta_down2 = (estimate_down - estimate);
    else delta_down2 = 0;
    //std::cout << "TRIG2: + " << delta_up2 << " - " << delta_down2 << std::endl;

    estimateerrTrigPlus  = sqrt(delta_up*delta_up + delta_up2*delta_up2);
    estimateerrTrigMinus = sqrt(delta_down*delta_down + delta_down2*delta_down2);
    //    std::cout << "TrigError: trig+ " << estimateerrTrigPlus << " trig- " << estimateerrTrigMinus << std::endl;
  }

  double closureStat2 = datamode? 0: jmt::errAoverB(SIG,SIGerr,estimate,estimateerr);
  double R0 = myR;
  double R0err = myRerr;

  TString name = region.btagSelection;
  name += region.owenId;
  name += isSIG ? ", SIG":", SB";
  char output[500];

  //for the purpose of the print-out only, revert D back to the observed data counts
  if(useScaleFactors_ && datamode){
    D = D_temp;
    Derr = Derr_temp;
  }

  if (!datamode) {
    float fc = (estimate-SIG)/sqrt(estimateerr*estimateerr + SIGerr*SIGerr);
    sprintf(output,"%s & %s & %s & %s & %s & %s & $%f \\pm %f$ \\\\ %% %f",name.Data(),
	    jmt::format_nevents(B,Berr).Data(),jmt::format_nevents(A,Aerr).Data(),
	    jmt::format_nevents(D,Derr).Data(),jmt::format_nevents(estimate,estimateerr).Data(),
	    jmt::format_nevents(SIG,SIGerr).Data(),100*(estimate-SIG)/estimate,closureStat2*100,fc);    
    cout<<output<<endl;
  }
  else {
    if(!(forTTbarEstimate)){
      if(reweightLSBdata_){
	sprintf(output,"%s & $%f \\pm %f$ & %d & %s & %s$^{+%.2f}_{-%.2f}$  \\\\ %% %f +/- %f",name.Data(),
		RLSB_RW,dRLSB_RW,
		TMath::Nint(D), jmt::format_nevents(Dsub,Dsuberr).Data(),
		jmt::format_nevents(estimate,estimateerr).Data(),estimateerrTrigPlus,estimateerrTrigMinus,R0,R0err);
	cout<<"(qcd) DATA\t";
      }
      else{
	sprintf(output,"%s & %d & %d & %d & %s & %s$^{+%.2f}_{-%.2f}$   \\\\ %% %f +/- %f",name.Data(),
		TMath::Nint(B),TMath::Nint(A),
		TMath::Nint(D), jmt::format_nevents(Dsub,Dsuberr).Data(),
		jmt::format_nevents(estimate,estimateerr).Data(),estimateerrTrigPlus,estimateerrTrigMinus,R0,R0err);
	cout<<"(qcd) DATA\t";
      }
      cout<<output<<endl;
    }
  }

  std::vector<double> v_errors;
  v_errors.push_back(estimateerr); 
  v_errors.push_back(estimateerrTrigPlus);
  v_errors.push_back(estimateerrTrigMinus);
  if (datamode)  return make_pair(estimate,v_errors);
  return make_pair( 100*sqrt(pow((SIG-estimate)/estimate,2) + pow(closureStat2,2)), v_errors); //estimate in denominator
}


void runClosureTest2011(std::map<TString, std::vector<double> > & syst, bool addnjetrw = true)  {
  
  setSearchRegions();
  assert(sbRegions_.size() == searchRegions_.size());

  for (unsigned int i=0; i<sbRegions_.size(); i++) {

    double sb =   fabs(anotherABCD(sbRegions_[i],false,1,0,getLSBbsel(sbRegions_[i])).first);
    double sig =  fabs(anotherABCD(searchRegions_[i],false,1,0,getLSBbsel(searchRegions_[i])).first);
    double sb_rw = 0, sig_rw = 0;
    double finalsb = sb;
    double finalsig = sig;

    if(addnjetrw){
      sb_rw = fabs(anotherABCD(sbRegions_[i],false,1,0,getLSBbsel(sbRegions_[i]),0,true).first);
      sig_rw = fabs(anotherABCD(searchRegions_[i],false,1,0,getLSBbsel(searchRegions_[i]),0,true).first);

      finalsb = (sb_rw>sb) ? sb_rw: sb;
      finalsig = (sig_rw>sig) ? sig_rw: sig;

      cout << "njetrw: ************** int " << i << " ****************" << endl;
      cout << "njetrw: sb: " << sb << ", sb_rw: " << sb_rw << endl;
      cout << "njetrw: sig: " << sig << ", sig_rw: " << sig_rw << endl;
      cout << "njetrw: ******************************************" << endl;
    }
    

    syst["Closure"].push_back( finalsb );
    syst["Closure"].push_back( finalsig ); 
  }
  
}

void runClosureTest2011()  {
  std::map<TString, std::vector<double> > dummy;
  
  //cout << "You are starting closure test without considering njet reweighting!" << endl;
  //runClosureTest2011(dummy,false);

  cout << "You are starting closure test that also considers njet reweighting! (this will do both)" << endl;
  runClosureTest2011(dummy,true);

}


void testQCD(unsigned int i) {
  //do it with shifted SB
  setSearchRegions();
  cout<<" SB +10 GeV"<<endl;
  anotherABCD(sbRegions_[i],true,1,10,getLSBbsel(sbRegions_[i]));


}

//to hold systematics
//making this global because I'm lazy...
//the string is the category of systematic
//the vector is list of values for all of the search and sb regions (not a very good approach)
std::map<TString, std::vector<double> > qcdSystErrors;
void runDataQCD2011(const bool forOwen=false) {

  const bool doNJetRWClosure = true; //should usually be true; set it to false only to save time
  if (!doNJetRWClosure) cout<<"Warning! njet reweighting for closure test is turned off!"<<endl;

  setSearchRegions();
  assert(sbRegions_.size() == searchRegions_.size());

  //data structures to hold results
  vector<std::pair<double,std::vector<double> > > n;
  vector<std::pair<double,std::vector<double> > > subp;
  vector<std::pair<double,std::vector<double> > > subm;
			  
  vector<std::pair<double,std::vector<double> > > sbp;
  vector<std::pair<double,std::vector<double> > > sbm;
			  
  vector<std::pair<double,std::vector<double> > > lsbp;
  vector<std::pair<double,std::vector<double> > > lsbm;

  cout<<" ==== Nominal data results === "<<endl;
  for (unsigned int i=0; i<sbRegions_.size(); i++) {
    n.push_back( anotherABCD(sbRegions_[i],true,1,0,getLSBbsel(sbRegions_[i])));
    std::pair<double,std::vector<double> > qcdresult= anotherABCD(searchRegions_[i],true,1,0,getLSBbsel(searchRegions_[i]) );
    n.push_back( qcdresult);
    resultsMap_[ searchRegions_[i].id()]["QCD"].value = qcdresult.first;
    resultsMap_[ searchRegions_[i].id()]["QCD"].statError = qcdresult.second.at(0);
    resultsMap_[ searchRegions_[i].id()]["QCD"].trigErrorPlus = qcdresult.second.at(1);
    resultsMap_[ searchRegions_[i].id()]["QCD"].trigErrorMinus = qcdresult.second.at(2);
  }
  cout<<" =END Nominal data results === "<<endl;

  /*
    the "owen mode" is collecting in global data structures the quantities that owen needs
    If we run the systematics then we will clobber the regular numbers. So quit now with just the nominal
    values.
  */
  if (forOwen) return;

  //now do it again with +40% subtraction
  cout<<" subtraction +40% "<<endl;
  for (unsigned int i=0; i<sbRegions_.size(); i++) {
    subp.push_back( anotherABCD(sbRegions_[i],true,1.4,0,getLSBbsel(sbRegions_[i])));
    subp.push_back( anotherABCD(searchRegions_[i],true,1.4,0,getLSBbsel(searchRegions_[i])));
  }
  
  //now do it again with -40% subtraction
  cout<<" subtraction -40% "<<endl;
  for (unsigned int i=0; i<sbRegions_.size(); i++) {
    subm.push_back( anotherABCD(sbRegions_[i],true,0.6,0,getLSBbsel(sbRegions_[i])));
    subm.push_back( anotherABCD(searchRegions_[i],true,0.6,0,getLSBbsel(searchRegions_[i])));
  }
  
  cout<<" ==== systematics for MC subtraction ==== "<<endl;
  for (unsigned int j=0; j<n.size(); j++) {
    double var1 = 100*(subp[j].first-n[j].first )/n[j].first;
    double var2 = 100*(subm[j].first-n[j].first )/n[j].first;
    //expect these to be equal and opposite
    if ( fabs(var1)-fabs(var2) > 0.1 ) cout<<j<<" found MC subtraction systematic size discrepancy"<<endl;
    cout<<var1<<"\t"<<var2<<endl;

    qcdSystErrors["MCsub"].push_back( fabs(var1)>fabs(var2) ? fabs(var1) : fabs(var2));
  }
  cout<<" =END systematics for MC subtraction ==== "<<endl;

  //do it with shifted SB
  cout<<" SB +10 GeV"<<endl;
  for (unsigned int i=0; i<sbRegions_.size(); i++) {
    sbp.push_back( anotherABCD(sbRegions_[i],true,1,10,getLSBbsel(sbRegions_[i])));
    sbp.push_back( anotherABCD(searchRegions_[i],true,1,10,getLSBbsel(searchRegions_[i])));
  }

  //now do it again with shift in the other direction
  cout<<" SB -10 GeV "<<endl;
  for (unsigned int i=0; i<sbRegions_.size(); i++) {
    TString lsbbsel = "==0"; //default is use 0b LSB
    if (btaggedLSB_) lsbbsel = sbRegions_[i].btagSelection; //in case we want btagged LSB
    sbm.push_back( anotherABCD(sbRegions_[i],true,1,-10,getLSBbsel(sbRegions_[i])));
    sbm.push_back( anotherABCD(searchRegions_[i],true,1,-10,getLSBbsel(searchRegions_[i])));
  }

  cout<<" ==== systematics for SB shift ==== "<<endl;
  for (unsigned int j=0; j<n.size(); j++) {
    double var1 = 100*(n[j].first  -sbp[j].first)/n[j].first;
    double var2 = 100*(n[j].first -sbm[j].first)/n[j].first;
    cout<<var1<<"\t"<<var2<<endl;
    //expect these to be different...no sanity check is warranted
    qcdSystErrors["SBshift"].push_back( fabs(var1)>fabs(var2) ? fabs(var1) : fabs(var2));
  }
  cout<<" =END systematics for SB shift ==== "<<endl;
  
  if(reweightLSBdata_){
    //do it with LSB RW correction +100%
    cout << "LSB RW correction +100%" << endl;
    for (unsigned int i=0; i<sbRegions_.size(); i++) {
      lsbp.push_back( anotherABCD(sbRegions_[i],true,1,0,getLSBbsel(sbRegions_[i]),1));
      lsbp.push_back( anotherABCD(searchRegions_[i],true,1,0,getLSBbsel(searchRegions_[i]),1));
    }
    
    //do it with LSB RW correction -100%
    cout << "LSB RW correction -100%" << endl;
    for (unsigned int i=0; i<sbRegions_.size(); i++) {
      lsbm.push_back( anotherABCD(sbRegions_[i],true,1,0,getLSBbsel(sbRegions_[i]),-1));
      lsbm.push_back( anotherABCD(searchRegions_[i],true,1,0,getLSBbsel(searchRegions_[i]),-1));
    }
    
    cout<<" ==== systematics for LSB reweighting ==== "<<endl;
    for (unsigned int j=0; j<n.size(); j++) {
      double var1 = 100*(n[j].first  -lsbp[j].first)/n[j].first;
      double var2 = 100*(n[j].first -lsbm[j].first)/n[j].first;
      if ( fabs(var1)-fabs(var2) > 0.1 ) cout<<j<<" found LSB reweighting systematic size discrepancy"<<endl;
      cout<<var1<<"\t"<<var2<<endl;
      qcdSystErrors["LSBrw"].push_back( fabs(var1)>fabs(var2) ? fabs(var1) : fabs(var2));
    }
    cout<<" ==== END systematics for LSB reweighting ==== "<<endl;
  }

  if (!btaggedLSB_) {
    cout<<"== Cross check with >=1 b instead of exactly 0 b =="<<endl;
    for (unsigned int i=0; i<sbRegions_.size(); i++) {
      anotherABCD(sbRegions_[i],true,1,0,">=1");
      anotherABCD(searchRegions_[i],true,1,0,">=1");
    }
  }
  
  cout<<" == Running QCD closure test =="<<endl;
  runClosureTest2011(qcdSystErrors,doNJetRWClosure);
  
  cout<<" == QCD systematics summary =="<<endl;
  cout.setf(ios_base::fixed); //want exactly one decimal place....
  cout.precision(1);

  /*
    stuffing the alternating SB and SIG into a vector was ugly, so here we have to be careful if we want to print out something
    human-readable as a label. Also we are hard-coding the fact that the *odd* elements are SIG
  */
  for (unsigned int j=0; j<n.size(); j++) {
    int index = ( j - (j%2))/2;
    SearchRegion thisregion = (j%2==0) ? sbRegions_.at(index) : searchRegions_.at(index);
    if( reweightLSBdata_){
      qcdSystErrors["Total"].push_back( sqrt( pow(qcdSystErrors["MCsub"].at(j),2) +  pow(qcdSystErrors["Closure"].at(j),2) + pow(qcdSystErrors["SBshift"].at(j),2) +  pow(qcdSystErrors["LSBrw"].at(j),2)  ));
      cout<<thisregion.id() <<"\t& $"<<qcdSystErrors["MCsub"].at(j)<<"$ & $"<<qcdSystErrors["Closure"].at(j)<<"$ & $"<<qcdSystErrors["SBshift"].at(j)<<"$ & $"<<qcdSystErrors["LSBrw"].at(j) << "$ & $" <<qcdSystErrors["Total"].at(j)<< "$ \\\\" << endl;
    }
    else{
      qcdSystErrors["Total"].push_back( sqrt( pow(qcdSystErrors["MCsub"].at(j),2) +  pow(qcdSystErrors["Closure"].at(j),2)+ pow(qcdSystErrors["SBshift"].at(j),2)));
      cout<<thisregion.id()<<"\t& $"<<qcdSystErrors["MCsub"].at(j)<<"$ & $"<<qcdSystErrors["Closure"].at(j)<<"$ & $"<<qcdSystErrors["SBshift"].at(j)<<"$ & $"<<qcdSystErrors["Total"].at(j)<< "$ \\\\"<< endl;
    }
    cout.unsetf(ios_base::fixed); //reset the cout precision
    cout.precision(6);

    //we only want to fill this data for the SIG regions
    if (thisregion.isSIG) resultsMap_[ thisregion.id()]["QCD"].systError = 0.01*qcdSystErrors["Total"].at(j)*resultsMap_[ thisregion.id()]["QCD"].value;
    //note that the qcdSystErrors is in %, so I need to convert to a number of events here
  }
  
}

double slABCD(const unsigned int searchRegionIndex, bool datamode=false, const TString & mode="", const TString & closureMode="nominal", bool fillResultsGlobal=false, const bool bShape=false ) {
  //in datamode, return the estimate; in non-datamode, return the Closure Test results (true-pred)/pred

/*
modifications for SHAPE analysis: (bShape = true)
-- new formula mu_i,jb = mu_SB,1b * ( mu_i,SL,jb / mu_SB,SL,1b )

-- for DATA, this means the data-driven subtraction must always be for the 1b case
*/

  assert( closureMode=="nominal" || closureMode=="justttbar" || closureMode=="wtplus" ||closureMode=="wtminus" ||closureMode=="slonlywtplus" ||closureMode=="slonlywtminus" ||closureMode=="0lonlywtplus" ||closureMode=="0lonlywtminus" || closureMode=="justw"|| closureMode=="justsingletop");

  if (closureMode=="justttbar") cout<<"Will run closure test in ttbar only mode!"<<endl;
  else if (closureMode=="nominal") {}
  else cout<<"Running closure test in mode: "<<closureMode<<endl;

  const SearchRegion region = searchRegions_[searchRegionIndex];
  SearchRegion qcdsubregion = sbRegions_[searchRegionIndex];
  if (bShape) { //qcdsubregion should be the eq1b region
    for (int rr = searchRegionIndex; rr>=0; rr--) { //loop backwards; assumes that the vector is filled a certain way, of course
      if (sbRegions_[rr].btagSelection=="eq1b") {
	qcdsubregion = sbRegions_[rr];
	break;
      }
    }
  }

  TString btagselection = region.btagSelection;
  TString owenKey = btagselection;
  owenKey += region.owenId;
  OwenData * myOwen = &(owenMap_[owenKey]);

  //we have to be careful. because of all of the global variables used in this code, we have to call
  //this other guy before we do anything else (should make it a class!)

  double SBsubQCD =0,SBsubQCDerr=0,SBsubZ=0,SBsubZerr=0,SBsubMisc=0,SBsubMiscerr=0;
  double SBsubQCDTrigErrPlus=0, SBsubQCDTrigErrMinus=0;
  std::pair<double,std::vector<double> > SBsubQCDp;
  std::pair<double,std::vector<double> > SIGQCDp;
  double SIGQCD=0, SIGQCDerr=0,SIGQCDTrigErrPlus=0,SIGQCDTrigErrMinus=0;
  if (datamode)  {
    //regions should only differ in the MET selection (except for b tag, in case of bShape analysis)
    if (!bShape) assert ( region.btagSelection == qcdsubregion.btagSelection);
    assert(region.htSelection == qcdsubregion.htSelection);
    assert(region.owenId == qcdsubregion.owenId);
    //get the QCD estimate for the SB region
    SBsubQCDp = anotherABCD(qcdsubregion, datamode, 1 ,0,getLSBbsel(qcdsubregion),0,false,true) ;//set "forTTbarEstimate" flag to true
    SBsubQCD=SBsubQCDp.first;
    SBsubQCDerr=SBsubQCDp.second.at(0);

    //Need the uncertainty on the qcd component due to trigger efficiency
    if(SBsubQCDp.second.size() < 3){ 
      std::cout << "Not all errors stored in qcd SB estimate error vector. Exiting." << std::endl;
      assert(0);
    }
    SBsubQCDTrigErrPlus = SBsubQCDp.second.at(1);
    SBsubQCDTrigErrMinus = SBsubQCDp.second.at(2);

    //for combining QCD and TTbar SIG estimate
    SIGQCDp = anotherABCD(region, datamode, 1 ,0,getLSBbsel(region),0,false,true) ;//set "forTTbarEstimate" flag to true
    if(SIGQCDp.second.size() < 3){ 
      std::cout << "Not all errors stored in qcd SIG estimate error vector. Exiting." << std::endl;
      assert(0);
    }
    SIGQCD=SIGQCDp.first;
    SIGQCDerr=SIGQCDp.second.at(0);
    SIGQCDTrigErrPlus  = SIGQCDp.second.at(1);
    SIGQCDTrigErrMinus = SIGQCDp.second.at(2);


    if (mode.Contains("QCD")) {
       //SBsubQCDerr is the statistical error on the QCD in SB.
      //but the systematic is needed too
      //factor 0.01 to convert from human-readable % to fractional
      double qcdsbsyst=0.01*qcdSystErrors["Total"].at(2*searchRegionIndex);
      cout<<"Found a QCD systematic of "<<100*qcdsbsyst<<endl;
      SBsubQCDerr = sqrt( SBsubQCDerr*SBsubQCDerr + pow(qcdsbsyst*SBsubQCD,2));
      
      if (mode=="QCDup") SBsubQCD += SBsubQCDerr;
      if (mode=="QCDdown") SBsubQCD -= SBsubQCDerr;
    }
    
    //now hard-coded Z->nunu
    double zv[2];
    double ze[2];
    //sometimes need to average mu mu and ee estimates
    double zsbsyst=0.5;
    double zinvscale = 4.68/4.65;
    bool doMean = true;

    if (qcdsubregion.owenId == "Loose" && qcdsubregion.btagSelection=="ge1b") {
      doMean=false;//averaging already done
      zv[0] = 82; ze[0]=20;
      zsbsyst = 0.0;//ze is stat+syst error combined
    }
    else if (qcdsubregion.owenId == "Tight" && qcdsubregion.btagSelection=="ge1b") {
      doMean=false;//averaging already done
      zv[0] = 44; ze[0]=13;
      zsbsyst = 0.0;
    }
    else if (qcdsubregion.owenId == "Loose" && qcdsubregion.btagSelection=="ge2b") {
      doMean=false;//averaging already done
      zv[0] = 14; ze[0]=9;
      zsbsyst = 0.0;
    }
    else if (qcdsubregion.owenId == "Tight" && qcdsubregion.btagSelection=="ge2b") {
      doMean=false;//averaging already done
      zv[0] = 3.8; ze[0]=2.7;
      zsbsyst = 0.0;
    }
    else if ( (qcdsubregion.owenId == "Loose" || qcdsubregion.owenId.Contains( "METBin")|| qcdsubregion.owenId.Contains( "METFineBin"))  && qcdsubregion.btagSelection=="ge3b") {
      doMean=false;//averaging already done
      zv[0] = 1.9; ze[0]=2.8;
      zsbsyst = 0.0;   
    }
    //for 150-200 GeV SB
    else if (qcdsubregion.owenId == "LooseLowSB" && qcdsubregion.btagSelection=="ge1b") {
      doMean=false;   zv[0] = 116.7; ze[0]=26.9;    zsbsyst = 0.0;
    }
    else if (qcdsubregion.owenId == "TightLowSB" && qcdsubregion.btagSelection=="ge1b") {
      doMean=false;   zv[0] = 58.9; ze[0]=15.8;    zsbsyst = 0.0;
    }
    else if (qcdsubregion.owenId == "LooseLowSB" && qcdsubregion.btagSelection=="ge2b") {
      doMean=false;   zv[0] = 19.5; ze[0]=12.9;    zsbsyst = 0.0;
    }
    else if (qcdsubregion.owenId == "TightLowSB" && qcdsubregion.btagSelection=="ge2b") {
      doMean=false;   zv[0] = 5.0; ze[0]=3.5;    zsbsyst = 0.0;
    }
    else if (qcdsubregion.owenId == "LooseLowSB" && qcdsubregion.btagSelection=="ge3b") {
      doMean=false;   zv[0] = 2.7; ze[0]=3.9;    zsbsyst = 0.0;
    }
    //for 150-250 GeV SB -- take 150-200 and add the 200-250 numbers. PRELIMINARY...probably the errors are wrong
    else if (qcdsubregion.owenId == "LooseWideSB" && qcdsubregion.btagSelection=="ge1b") {
      doMean=false;   zv[0] = 116.7 + 82; ze[0]=jmt::addInQuad(26.9,20);    zsbsyst = 0.0;
    }
    else if (qcdsubregion.owenId == "TightWideSB" && qcdsubregion.btagSelection=="ge1b") {
      doMean=false;   zv[0] = 58.9 + 44; ze[0]=jmt::addInQuad(15.8,13);    zsbsyst = 0.0;
    }
    else if (qcdsubregion.owenId == "LooseWideSB" && qcdsubregion.btagSelection=="ge2b") {
      doMean=false;   zv[0] = 19.5 + 14; ze[0]=jmt::addInQuad(12.9,9);    zsbsyst = 0.0;
    }
    else if (qcdsubregion.owenId == "TightWideSB" && qcdsubregion.btagSelection=="ge2b") {
      doMean=false;   zv[0] = 5.0 + 3.8; ze[0]=jmt::addInQuad(3.5,2.7);    zsbsyst = 0.0;
    }
    else if (qcdsubregion.owenId == "LooseWideSB" && qcdsubregion.btagSelection=="ge3b") {
      doMean=false;   zv[0] = 2.7 + 1.9; ze[0]=jmt::addInQuad(3.9,2.8);    zsbsyst = 0.0;
    }
    else if (qcdsubregion.owenId == "Loose" && qcdsubregion.btagSelection=="eq1b") {
      doMean=false;//averaging already done
      //preliminary numbers from http://www.slac.stanford.edu/~gaz/RA2b/ZinvTable.pdf
      zv[0] =68; ze[0]=15;
      zsbsyst = 0.0;    //this is a *percent* error...but we're not putting it in that way anymore
    }
    else if (qcdsubregion.owenId == "Loose" && qcdsubregion.btagSelection=="eq2b") {
      doMean=false;//averaging already done
      zv[0] =12; ze[0]=8;
      zsbsyst = 0.0;   
    }
    else if (qcdsubregion.owenId == "TightMET" && qcdsubregion.btagSelection=="eq1b") {
      doMean=false;//averaging already done
      zv[0] =36; ze[0]=10;
      zsbsyst = 0.0; 
    }
    else if (qcdsubregion.owenId == "TightMET" && qcdsubregion.btagSelection=="eq2b") {
      doMean=false;//averaging already done
      zv[0] =6; ze[0]=4;
      zsbsyst = 0.0; 
    }
    else if (qcdsubregion.owenId == "TightMET" && qcdsubregion.btagSelection=="ge3b") {
      doMean=false;//averaging already done
      zv[0] =1.0; ze[0]=1.4;
      zsbsyst = 0.0; 
    }
    else if (qcdsubregion.owenId == "TightHT" && qcdsubregion.btagSelection=="eq1b") {
      doMean=false;//averaging already done
      zv[0] =19; ze[0]=7;
      zsbsyst = 0.0; 
    }
    else if (qcdsubregion.owenId == "TightHT" && qcdsubregion.btagSelection=="eq2b") {
      doMean=false;//averaging already done
      zv[0] =3.3; ze[0]=2.4;
      zsbsyst = 0.0; 
    }
    else if (qcdsubregion.owenId == "TightHT" && qcdsubregion.btagSelection=="ge3b") {
      doMean=false;//averaging already done
      zv[0] =0.5; ze[0]=0.8;
      zsbsyst = 0.0; 
    }

    else {assert(0);}

    if(doMean){
      SBsubZ = zinvscale*(jmt::weightedMean(2,zv,ze));
      SBsubZerr = zinvscale*(jmt::weightedMean(2,zv,ze,true));
    }
    else{//if weighted averaging was already done
      SBsubZ = zinvscale*zv[0];
      SBsubZerr = zinvscale*ze[0];
    }

    if (mode.Contains("Z")) {
      SBsubZerr = sqrt( SBsubZerr*SBsubZerr + pow(zsbsyst*SBsubZ,2));
      
      if (mode=="Zup") SBsubZ += SBsubZerr;
      if (mode=="Zdown") SBsubZ -= SBsubZerr;
    }
  }

  setStackMode(false);
  doData(true);
  setQuiet(true);

  TString sampleOfInterest = "totalsm";
  useFlavorHistoryWeights_=false;
  loadSamples(); //needs to come first! this is why this should be turned into a class, with this guy in the ctor
  setColorScheme("nostack");
  clearSamples();
  if (datamode) {
    addSample("ZJets");
    addSample("VV");
    sampleOfInterest="data";
  }
  else {
    if (closureMode!="justw" && closureMode!="justsingletop")  addSample("TTbarJets");
    if (closureMode!="justttbar") {
      if (closureMode!="justsingletop")  addSample("WJets");
      if (closureMode!="justw")          addSample("SingleTop");
    }
  }

  savePlots_=false;

  setLogY(false);
  TString var,xtitle;
  int nbins;
  float low,high;

  //don: hard-code MHT efficiencies
  //offline cuts: HT>400, 0L (1L), mindphin>4
  float eff_MHT = 1, eff_MHT_err[2];//first entry is + error, second entry is - error
  float eff_SB_MHT = 1, eff_SB_MHT_err[2];
  float eff_SL_MHT = 1, eff_SL_MHT_err[2];
  float eff_1e_SB_MHT = 1, eff_1e_SB_MHT_err[2];
  float eff_1m_SB_MHT = 1, eff_1m_SB_MHT_err[2];

  
  // --  count events
  TCut btagcut = "nbjetsCSVM>=1";
  TCut btagcut_1b = bShape ? "nbjetsCSVM==1" : "nbjetsCSVM>=1"; //for the new shape analysis (and also for checks with a ge1b SL sample)
  btagSFweight_="1";
  TString btagSFweight=""; //we're going to have to switch this one in and out of the global var
  TString btagSFweight_1b="";
  if (useScaleFactors_) {
    btagcut="1";
    if (bShape)    {btagcut_1b="1"; btagSFweight_1b="prob1";} // shape
    else           {btagcut_1b="1"; btagSFweight_1b="probge1";} // ge1b SL cross-check
    usePUweight_=true;
    useHTeff_=true;
    useMHTeff_=false;
    thebnnMHTeffMode_ = kOff;//dealt with later
    currentConfig_=configDescriptions_.getCorrected(); //add JERbias
    
    eff_MHT       = eff_SIG_MHT_;    eff_MHT_err[0]       = eff_SIG_MHT_err_[0];    eff_MHT_err[1]       = eff_SIG_MHT_err_[1];
    eff_SB_MHT    = eff_SB_MHT_ ;    eff_SB_MHT_err[0]    = eff_SB_MHT_err_[0];     eff_SB_MHT_err[1]    = eff_SB_MHT_err_[1];
    eff_SL_MHT    = eff_SIG_SL_MHT_; eff_SL_MHT_err[0]    = eff_SIG_SL_MHT_err_[0]; eff_SL_MHT_err[1]    = eff_SIG_SL_MHT_err_[1];
    eff_1e_SB_MHT = eff_SB_1e_MHT_;  eff_1e_SB_MHT_err[0] = eff_SB_1e_MHT_err_[0];  eff_1e_SB_MHT_err[1] = eff_SB_1e_MHT_err_[1];
    eff_1m_SB_MHT = eff_SB_1m_MHT_;  eff_1m_SB_MHT_err[0] = eff_SB_1m_MHT_err_[0];  eff_1m_SB_MHT_err[1] = eff_SB_1m_MHT_err_[1];

    if (btagselection=="ge2b") {
      btagSFweight="probge2";
    }
    else if (btagselection=="ge1b") {
      btagSFweight="probge1";
    }
    else if (btagselection=="ge3b") {
      btagSFweight="probge3";
    }
    else if (btagselection=="eq1b") {
      btagSFweight="prob1";
    }
    else if (btagselection=="eq2b") {
      btagSFweight="(1-prob1-prob0-probge3)";
    }
    else {assert(0);}
    btagSFweight_=btagSFweight;
  }
  else {
    usePUweight_=false;
    useHTeff_=false;
    useMHTeff_=false;
    thebnnMHTeffMode_ = kOff;
    btagSFweight_="1";
    btagSFweight_1b="1"; //shape
    currentConfig_=configDescriptions_.getDefault(); //completely raw MC
    
    if (btagselection=="ge2b") {
      btagcut="nbjetsCSVM>=2";
    }
    else if (btagselection=="ge1b") {}
    else if (btagselection=="ge3b") {
      btagcut="nbjetsCSVM>=3";
    }
    else if (btagselection=="eq1b") {
      btagcut="nbjetsCSVM==1";
    }
    else if (btagselection=="eq2b") {
      btagcut="nbjetsCSVM==2";
    }
    else {assert(0);}
  }

  TCut HTcut=region.htSelection.Data(); 
  TCut SRMET = region.metSelection.Data();
  TCut baseline = "cutPV==1 && cut3Jets==1 && cutTrigger==1";
  baseline = baseline&&HTcut&&btagcut;
  TCut cleaning = "weight<1000 && passCleaning==1";
  TCut SBMET = qcdsubregion.metSelection.Data();
  TCut dpcut = "minDeltaPhiN>=4";
  TCut failOther = getSingleLeptonCut();//"(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)";
  TCut failOther_1e = getSingleElectronCut();//"(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)";
  TCut failOther_1m = getSingleMuonCut();//"(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)";
  TCut passOther = getLeptonVetoCut();//"nElectrons==0 && nMuons==0";

  double A_1e,A_1m,B,D,SIG,Aerr_1e,Aerr_1m,Berr,Derr,SIGerr;
  double B_bnnPlus=0, B_bnnMinus=0;
  //double AerrMinus,BerrMinus,DerrMinus;
  //A = SB, SL
  if (!datamode && (closureMode=="wtplus" )) {
    setSampleScaleFactor("WJets",1.38);
    setSampleScaleFactor("SingleTop",2);
  }
  else if (!datamode && (closureMode=="wtminus" )) {
    setSampleScaleFactor("WJets",1-0.38);
    setSampleScaleFactor("SingleTop",0);
  }
  else if (!datamode && ( closureMode=="slonlywtplus")) {
    setSampleScaleFactor("WJets",1.5);
    setSampleScaleFactor("SingleTop",1.5);
  }
  else if (!datamode && (closureMode=="slonlywtminus")) {
    setSampleScaleFactor("WJets",0.5);
    setSampleScaleFactor("SingleTop",0.5);
  }

  selection_ = baseline && cleaning && dpcut  && SBMET && failOther_1e; //auto cast to TString seems to work
  var="HT"; xtitle=var;
  nbins=10; low=0; high=5000;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_A");
  A_1e=getIntegral(sampleOfInterest);
  Aerr_1e=getIntegralErr(sampleOfInterest);

  selection_ = baseline && cleaning && dpcut  && SBMET && failOther_1m;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_A");
  A_1m=getIntegral(sampleOfInterest);
  Aerr_1m=getIntegralErr(sampleOfInterest);

  //  cout<<" SL SB found e,mu = "<<A_1e<<" , "<<A_1m<<endl;

  //shape -- get the SB,SL number for 1b
  //this code block was written for the SHAPE analysis but it turns out to be useful for overriding the b-tag cut to always use ge1b for the SB
  btagSFweight=btagSFweight_;
  selection_ = TCut("cutPV==1 && cut3Jets==1 && cutTrigger==1") &&HTcut && cleaning && dpcut  && SBMET && failOther_1e && btagcut_1b; //auto cast to TString seems to work
  btagSFweight_ = btagSFweight_1b;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_A");
  double mu_SB_1e_1b=getIntegral(sampleOfInterest);
  double mu_SB_1e_1b_err=getIntegralErr(sampleOfInterest);
  selection_ = TCut("cutPV==1 && cut3Jets==1 && cutTrigger==1") &&HTcut && cleaning && dpcut  && SBMET && failOther_1m && btagcut_1b;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_A");
  double mu_SB_1m_1b=getIntegral(sampleOfInterest);
  double mu_SB_1m_1b_err=getIntegralErr(sampleOfInterest);
  //put back the btag sf
  btagSFweight_ = btagSFweight;
  //end shape

  //D = SIG,SL
  //keep the scaling on W, t as set for region A
  selection_ = baseline && cleaning && dpcut  && SRMET && failOther; //auto cast to TString seems to work
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_D");
  D=getIntegral(sampleOfInterest);
  Derr=getIntegralErr(sampleOfInterest);

  //get the SIG,SL number for ge1b
  //this is only needed for a ge1b SL cross-check (not for the nominal method and not for the shape analysis)
  if (use1B_SL_) {
    btagSFweight=btagSFweight_;
    selection_ = TCut("cutPV==1 && cut3Jets==1 && cutTrigger==1") &&HTcut && cleaning && dpcut  && SRMET && failOther && btagcut_1b;
    btagSFweight_ = btagSFweight_1b;
    drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_A");
    D=getIntegral(sampleOfInterest); //clobber the nominal value of D
    Derr=getIntegralErr(sampleOfInterest);
    //put back the btag sf
    btagSFweight_ = btagSFweight;

    //finally, clobber A as well
    A_1e=mu_SB_1e_1b;
    Aerr_1e = mu_SB_1e_1b_err;

    A_1m=mu_SB_1m_1b;
    Aerr_1m = mu_SB_1m_1b_err;
  } //end new block

  if (datamode) { //for owen //SHAPE TODO
    myOwen->Nsig_sl = D;
    myOwen->Nsb_1e = A_1e;
    myOwen->Nsb_1m = A_1m;
  }

  double A_1e_temp = A_1e, Aerr_1e_temp = Aerr_1e;
  double A_1m_temp = A_1m, Aerr_1m_temp = Aerr_1m;
  double mu_SB_1e_1b_temp = mu_SB_1e_1b; //shape
  double mu_SB_1e_1b_err_temp = mu_SB_1e_1b_err; //shape
  double mu_SB_1m_1b_temp = mu_SB_1m_1b; //shape
  double mu_SB_1m_1b_err_temp = mu_SB_1m_1b_err; //shape
  double D_temp = D, Derr_temp = Derr;
  if(useScaleFactors_ && datamode) {
    if(useBNNEffCurves_){ 
      assert(0);
    }
    else {
      A_1e /= eff_1e_SB_MHT;
      A_1m /= eff_1m_SB_MHT;

      Aerr_1e /= eff_1e_SB_MHT;
      Aerr_1m /= eff_1m_SB_MHT;

      //Aerr = jmt::errAoverB(A_temp,Aerr_temp,eff_SL_SB_MHT,eff_SL_SB_MHT_err[0]);
      //AerrMinus = jmt::errAoverB(A_temp,Aerr_temp,eff_SL_SB_MHT,eff_SL_SB_MHT_err[1]);
      //Aerr = Aerr > AerrMinus ? Aerr: AerrMinus;
      
      D    /= eff_SL_MHT;
      Derr /= eff_SL_MHT; //treat trig eff error as 0, but still need to correct poisson error on D for trig eff
      //Derr = jmt::errAoverB(D_temp,Derr_temp,eff_SL_MHT,eff_SL_MHT_err[0]);
      //DerrMinus = jmt::errAoverB(D_temp,Derr_temp,eff_SL_MHT,eff_SL_MHT_err[1]);
      //Derr = Derr > DerrMinus ? Derr: DerrMinus;
      
      //shape
      mu_SB_1e_1b /= eff_1e_SB_MHT;
      mu_SB_1m_1b /= eff_1m_SB_MHT;

      mu_SB_1e_1b_err /= eff_1e_SB_MHT;
      mu_SB_1m_1b_err /= eff_1m_SB_MHT;
    }
  }


  //for wtplus, wtminus the scale factor is already set
  //for slonly modes, reset them to 1
  if (!datamode && ( closureMode=="slonlywtplus" || closureMode=="slonlywtminus")) {
    setSampleScaleFactor("WJets",1.0);
    setSampleScaleFactor("SingleTop",1.0);
  }
  else if (!datamode && ( closureMode=="0lonlywtplus")) {
    setSampleScaleFactor("WJets",1.5);
    setSampleScaleFactor("SingleTop",1.5);
  }
  else if (!datamode && (closureMode=="0lonlywtminus")) {
    setSampleScaleFactor("WJets",0.5);
    setSampleScaleFactor("SingleTop",0.5);
  }

  //B = SB
  btagSFweight=btagSFweight_;
  if (bShape) btagSFweight_ = btagSFweight_1b;
  if (bShape) selection_=TCut("cutPV==1 && cut3Jets==1 && cutTrigger==1") &&HTcut && cleaning && dpcut  && SBMET && passOther && btagcut_1b;
  else   selection_ = baseline && cleaning && dpcut  && SBMET && passOther;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_B");
  B=getIntegral(sampleOfInterest);
  Berr=getIntegralErr(sampleOfInterest);

  if (datamode) {
    SBsubMisc=getIntegral("totalsm");
    SBsubMiscerr=getIntegralErr("totalsm");
    if (mode.Contains("MC")) {
      const double mcsyst = 1; //100% uncertainty
      SBsubMiscerr = sqrt( SBsubMiscerr*SBsubMiscerr + pow(mcsyst*SBsubMisc,2));

      if (mode=="MCup") SBsubMisc += SBsubMiscerr;
      else if (mode=="MCdown") SBsubMisc -= SBsubMiscerr;
    }
  }

  double B_temp = B, Berr_temp = Berr;
  if(useScaleFactors_ && datamode){

    if(useBNNEffCurves_){
      thebnnMHTeffMode_ = kOn;
      drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_B");
      B=getIntegral(sampleOfInterest);
      Berr=getIntegralErr(sampleOfInterest);

      thebnnMHTeffMode_ = kOnPlus;
      drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_B");
      B_bnnPlus=getIntegral(sampleOfInterest);
      thebnnMHTeffMode_ = kOnMinus;
      drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_B");
      B_bnnMinus=getIntegral(sampleOfInterest);

      thebnnMHTeffMode_ = kOff;//turn it back off
    }
    else{
      B = B/eff_SB_MHT;
      Berr/=eff_SB_MHT; //need to scale Poisson errors by trig eff
      //Berr = jmt::errAoverB(B_temp,Berr_temp,eff_SB_MHT,eff_SB_MHT_err[0]);
      //BerrMinus = jmt::errAoverB(B_temp,Berr_temp,eff_SB_MHT,eff_SB_MHT_err[1]);
      //Berr = Berr > BerrMinus ? Berr: BerrMinus;
    }
  }
  //put back the btag sf
  btagSFweight_ = btagSFweight;

  double SBtrue=0,SBtrueErr=0;
  if (!datamode && bShape) {
    //for closure test of bShape, need the true number in SB
    selection_ = baseline && cleaning && dpcut  && SBMET && passOther;
    drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_B");
    SBtrue=getIntegral(sampleOfInterest);
    SBtrueErr=getIntegralErr(sampleOfInterest);
  }

  //SIG
  selection_ = baseline && cleaning && dpcut  && SRMET && passOther; //auto cast to TString seems to work
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_C");
  SIG=getIntegral(sampleOfInterest);
  SIGerr=getIntegralErr(sampleOfInterest);

  // now we've got our counts in all of the boxes; all have been corrected for trigger eff
  // now find the total number of 'perfect trigger' events in the SL SB
  double A = A_1e + A_1m;
  double Aerr = sqrt(Aerr_1e*Aerr_1e + Aerr_1m*Aerr_1m);

  double mu_SB_SL_1b = mu_SB_1e_1b + mu_SB_1m_1b;
  double mu_SB_SL_1b_err = sqrt(mu_SB_1e_1b_err*mu_SB_1e_1b_err + mu_SB_1m_1b_err*mu_SB_1m_1b_err);

  double suberr = sqrt(SBsubMiscerr*SBsubMiscerr + SBsubQCDerr*SBsubQCDerr + SBsubZerr*SBsubZerr);

  double estimate=0,estimateerr=0;
  double estimate_temp = estimate, estimateerr_temp = estimateerr;
  double estimateerrTrigPlus  = 0, estimateerrTrigMinus = 0;
  double estimateerrTrigPlus_unscaled  = 0, estimateerrTrigMinus_unscaled = 0;

  double estimateQCDANDTTbar=0, estimateerrQCDANDTTbar=0, 
    estimateerrQCDANDTTbarTrigPlus=0, estimateerrQCDANDTTbarTrigMinus=0; 

  TString name = btagselection;
  name += region.owenId;

  double SBestimate=0,SBestimate_err=0;
  if (bShape) { //calculate SB and SIG estimates
    double Bprime = B-SBsubMisc-SBsubQCD-SBsubZ;
    double Bprimeerr= sqrt(suberr*suberr + Berr*Berr);

    SBestimate = Bprime * A/ mu_SB_SL_1b; 
    SBestimate_err = jmt::errAoverB(A,Aerr,mu_SB_SL_1b,mu_SB_SL_1b_err);
    SBestimate_err = jmt::errAtimesB(Bprime,Bprimeerr,A / mu_SB_SL_1b,SBestimate_err);
    
    double SIGestimate = Bprime * D / mu_SB_SL_1b;
    double SIGestimate_err = jmt::errAoverB(D,Derr,mu_SB_SL_1b,mu_SB_SL_1b_err);
    SIGestimate_err = jmt::errAtimesB(Bprime,Bprimeerr,D / mu_SB_SL_1b,SIGestimate_err);
    //closure test output
    cout<<"SHAPE "<<name.Data()
	<<" SB: " <<SBestimate<<" +/- "<<SBestimate_err<<" true: "<<SBtrue<<" +/- "<<SBtrueErr
	<<" SIG: "<<SIGestimate<<" +/- "<<SIGestimate_err<<" true: "<<SIG<<" +/- "<<SIGerr<<endl;
    estimate = SIGestimate;
    estimateerr = SIGestimate_err;
  }
  else {  //now calculate B*D/A (old-fashioned ABCD)
    double numerr=jmt::errAtimesB(B-SBsubMisc-SBsubQCD-SBsubZ,sqrt(suberr*suberr + Berr*Berr),
				  D,Derr);
    double num = (B-SBsubMisc-SBsubQCD-SBsubZ)*D;
    estimate = num / A;
    estimateerr= jmt::errAoverB(num,numerr,A,Aerr);
  }
  estimate_temp = estimate; estimateerr_temp = estimateerr;

  if(useScaleFactors_ && datamode){
    
    //double estimateerrMinus = 0;
    estimate = eff_MHT*estimate_temp;
    estimateerr = eff_MHT*estimateerr_temp;
    //estimateerr = jmt::errAtimesB(estimate_temp, estimateerr_temp, eff_MHT, eff_MHT_err[0]);
    // estimateerrMinus = jmt::errAtimesB(estimate_temp, estimateerr_temp, eff_MHT, eff_MHT_err[1]);
    //std::cout << "estimate: "<< estimate << "+" << estimateerr << "- " <<estimateerrMinus << std::endl;
    //use the larger of the two as the error in the estimate
    //estimateerr = estimateerr > estimateerrMinus ? estimateerr: estimateerrMinus;

    //for bShape, divide by SB,SL,1b; for normal, divide by A
    double mu_SB_SL_jb = bShape ? mu_SB_SL_1b : A; 

    //vary the trigger efficiencies by their statistical error to get trigger error on estimate
    //first vary QCD trigger efficiency
    float SBsubQCD_up   = SBsubQCD + SBsubQCDTrigErrPlus;
    double    num = (B-SBsubMisc-SBsubQCD_up-SBsubZ)*D;
    double estimate_up;
    estimate_up = eff_MHT*num/  mu_SB_SL_jb;
    double delta_down = (estimate_up - estimate);    
    estimate_up = num/  mu_SB_SL_jb;
    double delta_down_unscaled = (estimate_up - estimate_temp);    

    float SBsubQCD_down   = SBsubQCD - SBsubQCDTrigErrMinus;
    num = (B-SBsubMisc-SBsubQCD_down-SBsubZ)*D;
    double estimate_down;
    estimate_down= eff_MHT*num/  mu_SB_SL_jb;
    double delta_up = (estimate_down - estimate); 
    estimate_down= num/  mu_SB_SL_jb;
    double delta_up_unscaled = (estimate_down - estimate_temp); 

    //second vary eff_MHT   
    float eff_MHT_up   = eff_MHT + eff_MHT_err[0];
    estimate_up = eff_MHT_up*estimate_temp;
    double delta_up2;
    delta_up2 = (estimate_up - estimate);

    float eff_MHT_down = eff_MHT - eff_MHT_err[1];
    estimate_down = eff_MHT_down*estimate_temp;
    double delta_down2;
    delta_down2 = (estimate_down - estimate);

    //third vary eff_SB_MHT
    double estimate_up_unscaled, estimate_down_unscaled;
    if(useBNNEffCurves_){
      num = (B_bnnPlus-SBsubMisc-SBsubQCD-SBsubZ)*D;
      estimate_up = eff_MHT*num/  mu_SB_SL_jb;
      estimate_up_unscaled = num/  mu_SB_SL_jb;

      num = (B_bnnMinus-SBsubMisc-SBsubQCD-SBsubZ)*D;
      estimate_down = eff_MHT*num/  mu_SB_SL_jb;
      estimate_down_unscaled = num/  mu_SB_SL_jb;
    }
    else{
      float eff_SB_MHT_up   = eff_SB_MHT + eff_SB_MHT_err[0];
      double B_down = B_temp/eff_SB_MHT_up;
      num = (B_down-SBsubMisc-SBsubQCD-SBsubZ)*D;
      estimate_up = eff_MHT*num/  mu_SB_SL_jb;
      estimate_up_unscaled = num/  mu_SB_SL_jb;
      
      float eff_SB_MHT_down   = eff_SB_MHT - eff_SB_MHT_err[1];
      double B_up = B_temp/eff_SB_MHT_down;
      num = (B_up-SBsubMisc-SBsubQCD-SBsubZ)*D;
      estimate_down = eff_MHT*num/  mu_SB_SL_jb;
      estimate_down_unscaled = num/  mu_SB_SL_jb;
    }
    double delta_down3 = (estimate_up - estimate);
    double delta_up3 = (estimate_down - estimate);

    double delta_down3_unscaled = (estimate_up_unscaled - estimate_temp);
    double delta_up3_unscaled = (estimate_down_unscaled - estimate_temp);


    estimateerrTrigPlus  = sqrt(delta_up*delta_up + delta_up2*delta_up2 + delta_up3*delta_up3);
    estimateerrTrigMinus = sqrt(delta_down*delta_down + delta_down2*delta_down2 + delta_down3*delta_down3);

    estimateerrTrigPlus_unscaled  = sqrt(delta_up_unscaled*delta_up_unscaled + delta_up3_unscaled*delta_up3_unscaled);
    estimateerrTrigMinus_unscaled = sqrt(delta_down_unscaled*delta_down_unscaled + delta_down3_unscaled*delta_down3_unscaled);


    /////////////////////////////////
    //combine QCD and ttbar estimates
    /////////////////////////////////
    estimateQCDANDTTbar = eff_MHT*(estimate_temp+SIGQCD);
    estimateerrQCDANDTTbar = eff_MHT*sqrt(estimateerr_temp*estimateerr_temp + SIGQCDerr*SIGQCDerr);

    //vary eff_MHT   
    estimate_up = eff_MHT_up*(estimate_temp+SIGQCD);
    double delta_qcdandttbar_up = (estimate_up - estimateQCDANDTTbar);
    estimate_down = eff_MHT_down*(estimate_temp+SIGQCD);
    double delta_qcdandttbar_down = (estimate_down - estimateQCDANDTTbar);

    //vary qcd component
    estimate_up = eff_MHT*(estimate_temp + (SIGQCD+SIGQCDTrigErrPlus) );
    double delta_qcdandttbar_up2 = (estimate_up - estimateQCDANDTTbar);
    estimate_down = eff_MHT*(estimate_temp + (SIGQCD-SIGQCDTrigErrMinus) );
    double delta_qcdandttbar_down2 = (estimate_down - estimateQCDANDTTbar);

    //vary ttbar component 
    estimate_up = eff_MHT*( (estimate_temp+estimateerrTrigPlus_unscaled) + SIGQCD );
    double delta_qcdandttbar_up3 = (estimate_up - estimateQCDANDTTbar);
    estimate_down = eff_MHT*( (estimate_temp-estimateerrTrigMinus_unscaled) + SIGQCD );
    double delta_qcdandttbar_down3 = (estimate_down - estimateQCDANDTTbar);

    estimateerrQCDANDTTbarTrigPlus = sqrt( delta_qcdandttbar_up*delta_qcdandttbar_up 
					   + delta_qcdandttbar_up2*delta_qcdandttbar_up2 
					   + delta_qcdandttbar_up3*delta_qcdandttbar_up3 );

    estimateerrQCDANDTTbarTrigMinus = sqrt( delta_qcdandttbar_down*delta_qcdandttbar_down 
					   + delta_qcdandttbar_down2*delta_qcdandttbar_down2 
					   + delta_qcdandttbar_down3*delta_qcdandttbar_down3 );

  }
  

  double closureStat= datamode? 0: jmt::errAoverB(SIG,SIGerr,estimate,estimateerr);


//   cout<<" ==== "<<endl
//       <<"Estimate = "<<estimate<<" +/- "<<estimateerr<<endl
//       <<" truth   = "<<SIG     <<" +/- "<<SIGerr<<endl;
  // btagselection += tight ? " Tight " : " Loose ";



  char output[500];
  char output2[500];

  //for the purpose of the print-out only, revert A,B,D back to the observed data counts
  if(useScaleFactors_ && datamode){
    D = D_temp;
    Derr = Derr_temp;
    A = A_1e_temp + A_1m_temp; //FIXME -- should really output both e and mu separately
    Aerr = sqrt(Aerr_1e_temp*Aerr_1e_temp + Aerr_1m_temp*Aerr_1m_temp);
    B = B_temp;
    Berr = Berr_temp;
    mu_SB_SL_1b = mu_SB_1e_1b_temp+mu_SB_1m_1b_temp;
    mu_SB_SL_1b_err = 0;//mu_SB_SL_1b_err_temp; //FIXME
  }

  if (!datamode) {
    if (bShape) {
      sprintf(output,"%s SB  & %s & %s & %s & %s & %s & $%f \\pm %f$ \\\\ ",name.Data(),
	      jmt::format_nevents(B,Berr).Data(),jmt::format_nevents(mu_SB_SL_1b,mu_SB_SL_1b_err).Data(),
	      jmt::format_nevents(A,Aerr).Data(),jmt::format_nevents(SBestimate,SBestimate_err).Data(),
	      jmt::format_nevents(SBtrue,SBtrueErr).Data(),100*(SBestimate-SBtrue)/SBestimate,
	      100*jmt::errAoverB(SBtrue,SBtrueErr,SBestimate,SBestimate_err));
      sprintf(output2,"%s SIG & %s & %s & %s & %s & %s & $%f \\pm %f$ \\\\ ",name.Data(),
	      jmt::format_nevents(B,Berr).Data(),jmt::format_nevents(mu_SB_SL_1b,mu_SB_SL_1b_err).Data(),
	      jmt::format_nevents(D,Derr).Data(),jmt::format_nevents(estimate,estimateerr).Data(),
	      jmt::format_nevents(SIG,SIGerr).Data(),100*(estimate-SIG)/estimate,
	      100*closureStat);
    }
    else {
    float fc = (estimate-SIG)/sqrt(estimateerr*estimateerr + SIGerr*SIGerr);
      sprintf(output,"%s & %s & %s & %s & %s & %s & $%f \\pm %f$ \\\\ %% %f",name.Data(),
	    jmt::format_nevents(D,Derr).Data(),jmt::format_nevents(A,Aerr).Data(),
	    jmt::format_nevents(B,Berr).Data(),jmt::format_nevents(estimate,estimateerr).Data(),
	      jmt::format_nevents(SIG,SIGerr).Data(), 100*(estimate-SIG)/estimate, 100*closureStat, fc);
    }
  }
  else {
    if (!bShape)
      sprintf(output,"ttbar DATA %s & %d & %d & %d & %s & %s$^{+%.2f}_{-%.2f}$  \\\\",name.Data(),
	      TMath::Nint(D),TMath::Nint(A),
	      TMath::Nint(B), jmt::format_nevents(SBsubMisc+SBsubQCD+SBsubZ,suberr).Data(),
	      jmt::format_nevents(estimate,estimateerr).Data(),estimateerrTrigPlus,estimateerrTrigMinus);
    else 
      sprintf(output,"ttbar SHAPE DATA %s & %d & %d & %d & %s & %s$^{+%.2f}_{-%.2f}$  \\\\",name.Data(),
	      TMath::Nint(D),TMath::Nint(mu_SB_SL_1b),
	      TMath::Nint(B), jmt::format_nevents(SBsubMisc+SBsubQCD+SBsubZ,suberr).Data(),
	      jmt::format_nevents(estimate,estimateerr).Data(),estimateerrTrigPlus,estimateerrTrigMinus);

    sprintf(output2,"ttbar+QCD DATA %s & %s$^{+%.2f}_{-%.2f}$  \\\\",name.Data(),
	    jmt::format_nevents(estimateQCDANDTTbar,estimateerrQCDANDTTbar).Data(),
	    estimateerrQCDANDTTbarTrigPlus,estimateerrQCDANDTTbarTrigMinus);

    //stuff results into the resultsMap_
    if (fillResultsGlobal) {
      resultsMap_[ region.id()]["ttbar"].value = estimate;
      resultsMap_[ region.id()]["ttbar"].statError = estimateerr;
      resultsMap_[ region.id()]["ttbar"].trigErrorPlus = estimateerrTrigPlus;
      resultsMap_[ region.id()]["ttbar"].trigErrorMinus = estimateerrTrigMinus;
      
      resultsMap_[ region.id()]["ttbarQCD"].value = estimateQCDANDTTbar;
      resultsMap_[ region.id()]["ttbarQCD"].statError = estimateerrQCDANDTTbar;
      resultsMap_[ region.id()]["ttbarQCD"].trigErrorPlus = estimateerrQCDANDTTbarTrigPlus;
      resultsMap_[ region.id()]["ttbarQCD"].trigErrorMinus = estimateerrQCDANDTTbarTrigMinus;
    }
  }
  cout<<output<<endl;
  if(datamode ||bShape) cout << output2 << endl;

  resetSampleScaleFactors();

  if (datamode) return estimate;
  return 100*sqrt(pow((SIG-estimate)/estimate,2) + pow(closureStat,2));
}

void runSLClosureTest2011(const bool bShape=false) {

  setSearchRegions();
  cout<<"Note that the following is the joint tt+W+t closure test only!"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++)    slABCD(j,false,"","nominal",false,bShape);

/*
  cout<<"Note that the following is the tt closure test only!"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) slABCD(j,false,"","justttbar");
  for (unsigned int j=0; j<searchRegions_.size();j++) slABCD(j,false,"","justw");
  for (unsigned int j=0; j<searchRegions_.size();j++) slABCD(j,false,"","justsingletop");
*/
  /*
  for (unsigned int j=0; j<searchRegions_.size();j++) slABCD(j,false,"","wtplus");
  for (unsigned int j=0; j<searchRegions_.size();j++) slABCD(j,false,"","wtminus");

  for (unsigned int j=0; j<searchRegions_.size();j++) slABCD(j,false,"","slonlywtplus");
  for (unsigned int j=0; j<searchRegions_.size();j++) slABCD(j,false,"","slonlywtminus");

  for (unsigned int j=0; j<searchRegions_.size();j++) slABCD(j,false,"","0lonlywtplus");
  for (unsigned int j=0; j<searchRegions_.size();j++) slABCD(j,false,"","0lonlywtminus");
  */
}

vector<double> ttbarClosureSyst;
void runTtbarEstimate2011(const bool forOwen=false, const bool bShape=false) {
  setSearchRegions();

  std::vector<double> ttw;
  std::vector<double> p;
  std::vector<double> m;
  
  std::vector<double> qcd;
  std::vector<double> znn;
  std::vector<double> mc;
  //  int i=0;
  cout<<"nominal ttbar"<<endl;
  //this bool set to true at the end tells it to fill the global data for tex output
  //this is not a pretty solution
  for (unsigned int j=0; j<searchRegions_.size();j++) ttw.push_back(slABCD(j,true,"","nominal",true,bShape));
 
  if (forOwen) return;

  // for systematics due to QCD subtraction
  //must do   runDataQCD2011(); first (in order to fill syst errors for qcd)
  if (  qcdSystErrors.size()==0) {cout<<"Need to do runDataQCD2011() first!"<<endl; return;}

  cout<<"vary QCD subtraction"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) {
    p.push_back(slABCD(j,true,"QCDup","nominal",false,bShape));
    m.push_back(slABCD(j,true,"QCDdown","nominal",false,bShape));
    if (fabs(p[j]-ttw[j])/ttw[j] - fabs(m[j]-ttw[j])/ttw[j] > 0.01) cout<<"** Difference in size of QCD sub syst!"<<endl;
    qcd.push_back(100*fabs(p[j]-ttw[j])/ttw[j]);
  }

  p.clear(); m.clear();
  //for Znunu subtraction systematics
  cout<<"vary Znn subtraction"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) {
    p.push_back(slABCD(j,true,"Zup","nominal",false,bShape));
    m.push_back(slABCD(j,true,"Zdown","nominal",false,bShape));
    if (fabs(p[j]-ttw[j])/ttw[j] - fabs(m[j]-ttw[j])/ttw[j] > 0.01) cout<<"** Difference in size of Znn sub syst!"<<endl;
    znn.push_back(100*fabs(p[j]-ttw[j])/ttw[j]);
  }

  p.clear(); m.clear();
  //for MC subtraction systematics
  cout<<"vary MC subtraction"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) {
    p.push_back(slABCD(j,true,"MCup","nominal",false,bShape));
    m.push_back(slABCD(j,true,"MCdown","nominal",false,bShape));
    if (fabs(p[j]-ttw[j])/ttw[j] - fabs(m[j]-ttw[j])/ttw[j] > 0.01) cout<<"** Difference in size of MC sub syst!"<<endl;
    mc.push_back(100*fabs(p[j]-ttw[j])/ttw[j]);
  }

  //finally, run the closure test
  std::vector<double> closure;
  cout<<"Running closure tests"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) {
    double allsamples=fabs(slABCD(j,false,"","nominal",false,bShape)); //mix of samples
    //    double ttbaronly=fabs(slABCD(j,false,"","justttbar")); //just ttbar
    double wtup=fabs(slABCD(j,false,"","wtplus",false,bShape)); //vary W+t
    double wtdown=fabs(slABCD(j,false,"","wtminus",false,bShape)); //vary W+t
    //find the largest of the 3 and store it in allsamples
    if (wtup > wtdown) wtdown = wtup;
    if (wtdown > allsamples) allsamples=wtdown;
    closure.push_back(allsamples);
  }

  ttbarClosureSyst.clear();
  cout<<" == summary (%) == "<<endl;
  cout<<"\tClosure\tQCD\tZ\tMC\tTotal"<<endl;
  cout.setf(ios_base::fixed); //want exactly one decimal place....
  cout.precision(1);

  for (unsigned int j=0; j<searchRegions_.size();j++) {
    double totalsyst = sqrt(closure[j]*closure[j] + qcd[j]*qcd[j] + znn[j]*znn[j] + mc[j]*mc[j]);
    cout<<searchRegions_[j].id() <<"\t & $"<<closure[j]<<"$ & $"<<qcd[j]<<"$ & $"<<znn[j]<<"$ & $"<<mc[j]<<"$ & $"<<totalsyst << "$ \\\\ " << endl;
    ttbarClosureSyst.push_back(closure[j]);

    //i have gone into slABCD and filled the global with central values and stat errors directly inside
    //but syst errors must be done here (note that they are in % so we must convert to events)
    resultsMap_[ searchRegions_[j].id()]["ttbar"].systError = 0.01*totalsyst*resultsMap_[ searchRegions_[j].id()]["ttbar"].value;
    //again, all errors except stat have been filled in slABCD
    //these are uncorrelated so let's just add them in quadrature
    resultsMap_[ searchRegions_[j].id()]["ttbarQCD"].systError = sqrt(
								      pow(resultsMap_[ searchRegions_[j].id()]["ttbar"].systError,2) 
								      + pow(resultsMap_[searchRegions_[j].id()]["QCD"].systError,2));
  }
  cout.unsetf(ios_base::fixed); //want exactly one decimal place....
  cout.precision(6);

 

}

//these 2 functions run the bshape version of the analysis
void runDDbShape() {

  runDataQCD2011(true);
  runTtbarEstimate2011(true,true); //turn on shape analysis

}
void runDDwithSystShape() {

  runDataQCD2011(false);
  runTtbarEstimate2011(false,true); //turn on shape analysis

}


void printOwenAll() { //no b shape

  runDataQCD2011(true);
  runTtbarEstimate2011(true,false); //turn off shape analysis

  for (std::map<TString,OwenData>::iterator i=owenMap_.begin(); i!=owenMap_.end(); ++i) {
    printOwen( i->first);
  }
  cout<<"luminosity = "<<lumiScale_<<endl;

}

void printOwenSyst(TString regions) {
  setSearchRegions(regions);
  //note that because of the way the code is written, i am not at all sure that
  //it is safe to run printOwenAll() and printOwenSyst() in the same session. better to be safe!

  runDataQCD2011();
  runTtbarEstimate2011(); //no shape analysis

  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    char effoutput[500];
    sprintf(effoutput,"backgroundSyst.%s%s.%dinvpb.dat",searchRegions_[i].btagSelection.Data(),searchRegions_[i].owenId.Data(),TMath::Nint(lumiScale_));

    ofstream textfile(effoutput);

    textfile<<"sf_mc            "<<1<<endl;
    textfile<<"sf_mc_err        "<<0.4<<endl;
    textfile<<"sf_qcd_sb        "<<1<<endl;
    textfile<<"sf_qcd_sb_err    "<<0.01*sqrt(pow(qcdSystErrors["Closure"].at(2*i),2)+pow(qcdSystErrors["SBshift"].at(2*i),2) + pow(qcdSystErrors["LSBrw"].at(2*i),2))<<endl;
    textfile<<"sf_qcd_sig       "<<1<<endl;
    textfile<<"sf_qcd_sig_err   "<<0.01*sqrt(pow(qcdSystErrors["Closure"].at(2*i+1),2)+pow(qcdSystErrors["SBshift"].at(2*i+1),2)+ pow(qcdSystErrors["LSBrw"].at(2*i+1),2))<<endl;
    textfile<<"sf_ttwj_sig      "<<1<<endl;
    textfile<<"sf_ttwj_sig_err  "<<0.01*ttbarClosureSyst[i]<<endl;
    textfile.close();
  }

  if (regions=="Moriond") {
    cout<<" == Last but not least, output for Results.tex =="<<endl;
    resultsMap_["ge1bLoose"]["ttbarCC"].value = 322;
    resultsMap_["ge1bLoose"]["ttbarCC"].statError = 30;
    resultsMap_["ge1bLoose"]["ttbarCC"].systError = 37;
    resultsMap_["ge1bLoose"]["ttbarCC"].trigErrorPlus = 0;
    resultsMap_["ge1bLoose"]["ttbarCC"].trigErrorMinus = 0;
    
    resultsMap_["ge1bTight"]["ttbarCC"].value = 5.0;
    resultsMap_["ge1bTight"]["ttbarCC"].statError = 2.3;
    resultsMap_["ge1bTight"]["ttbarCC"].systError = 2.2;
    resultsMap_["ge1bTight"]["ttbarCC"].trigErrorPlus = 0;
    resultsMap_["ge1bTight"]["ttbarCC"].trigErrorMinus = 0;

  resultsMap_["ge2bLoose"]["ttbarCC"].value = 123;
  resultsMap_["ge2bLoose"]["ttbarCC"].statError = 16;
  resultsMap_["ge2bLoose"]["ttbarCC"].systError = 16;
  resultsMap_["ge2bLoose"]["ttbarCC"].trigErrorPlus = 0;
  resultsMap_["ge2bLoose"]["ttbarCC"].trigErrorMinus = 0;

  resultsMap_["ge2bTight"]["ttbarCC"].value = 21.0;
  resultsMap_["ge2bTight"]["ttbarCC"].statError = 6.3;
  resultsMap_["ge2bTight"]["ttbarCC"].systError = 4.3;
  resultsMap_["ge2bTight"]["ttbarCC"].trigErrorPlus = 0;
  resultsMap_["ge2bTight"]["ttbarCC"].trigErrorMinus = 0;

  resultsMap_["ge3bLoose"]["ttbarCC"].value = 13.7;
  resultsMap_["ge3bLoose"]["ttbarCC"].statError = 6.1;
  resultsMap_["ge3bLoose"]["ttbarCC"].systError = 3.1;
  resultsMap_["ge3bLoose"]["ttbarCC"].trigErrorPlus = 0;
  resultsMap_["ge3bLoose"]["ttbarCC"].trigErrorMinus = 0;

  resultsMap_["ge1bLoose"]["Znn"].value = 132;
  resultsMap_["ge1bLoose"]["Znn"].statError = 17;
  resultsMap_["ge1bLoose"]["Znn"].systError = 23;
  resultsMap_["ge1bLoose"]["Znn"].trigErrorPlus = 0;
  resultsMap_["ge1bLoose"]["Znn"].trigErrorMinus = 0;

  resultsMap_["ge1bTight"]["Znn"].value = 2.5;
  resultsMap_["ge1bTight"]["Znn"].statError = 1.9;
  resultsMap_["ge1bTight"]["Znn"].systError = 0.5;
  resultsMap_["ge1bTight"]["Znn"].trigErrorPlus = 0;
  resultsMap_["ge1bTight"]["Znn"].trigErrorMinus = 0;

  resultsMap_["ge2bLoose"]["Znn"].value = 23;
  resultsMap_["ge2bLoose"]["Znn"].statError = 3.5;
  resultsMap_["ge2bLoose"]["Znn"].systError = 11.5;
  resultsMap_["ge2bLoose"]["Znn"].trigErrorPlus = 0;
  resultsMap_["ge2bLoose"]["Znn"].trigErrorMinus = 0;

  resultsMap_["ge2bTight"]["Znn"].value = 4.6;
  resultsMap_["ge2bTight"]["Znn"].statError = 1.4;
  resultsMap_["ge2bTight"]["Znn"].systError = 2.4;
  resultsMap_["ge2bTight"]["Znn"].trigErrorPlus = 0;
  resultsMap_["ge2bTight"]["Znn"].trigErrorMinus = 0;

  resultsMap_["ge3bLoose"]["Znn"].value = 3;
  resultsMap_["ge3bLoose"]["Znn"].statError = 0.9;
  resultsMap_["ge3bLoose"]["Znn"].systError = 3.3;
  resultsMap_["ge3bLoose"]["Znn"].trigErrorPlus = 0;
  resultsMap_["ge3bLoose"]["Znn"].trigErrorMinus = 0;

  cout<<"    & 1BL & 1BT & 2BL & 2BT & 3B \\"<<endl;
  //we'll do this in a C way instead of a C++ way
  cout<<"QCD & "
      <<formatLatex( resultsMap_["ge1bLoose"]["QCD"])<<" & "
      <<formatLatex( resultsMap_["ge1bTight"]["QCD"])<<" & "
      <<formatLatex( resultsMap_["ge2bLoose"]["QCD"])<<" & "
      <<formatLatex( resultsMap_["ge2bTight"]["QCD"])<<" & "
      <<formatLatex( resultsMap_["ge3bLoose"]["QCD"])<<" \\ "<<endl;

  cout<<"\\ttbar/W+jets & "
      <<formatLatex( resultsMap_["ge1bLoose"]["ttbar"])<<" & "
      <<formatLatex( resultsMap_["ge1bTight"]["ttbar"])<<" & "
      <<formatLatex( resultsMap_["ge2bLoose"]["ttbar"])<<" & "
      <<formatLatex( resultsMap_["ge2bTight"]["ttbar"])<<" & "
      <<formatLatex( resultsMap_["ge3bLoose"]["ttbar"])<<" \\ "<<endl;
  
  cout<<"\\ttbar/W+jets CC & "
      <<formatLatex( resultsMap_["ge1bLoose"]["ttbarCC"])<<" & "
      <<formatLatex( resultsMap_["ge1bTight"]["ttbarCC"])<<" & "
      <<formatLatex( resultsMap_["ge2bLoose"]["ttbarCC"])<<" & "
      <<formatLatex( resultsMap_["ge2bTight"]["ttbarCC"])<<" & "
      <<formatLatex( resultsMap_["ge3bLoose"]["ttbarCC"])<<" \\ "<<endl;
  
  cout<<"\\Zinvisible & "
      <<formatLatex( resultsMap_["ge1bLoose"]["Znn"])<<" & "
      <<formatLatex( resultsMap_["ge1bTight"]["Znn"])<<" & "
      <<formatLatex( resultsMap_["ge2bLoose"]["Znn"])<<" & "
      <<formatLatex( resultsMap_["ge2bTight"]["Znn"])<<" & "
      <<formatLatex( resultsMap_["ge3bLoose"]["Znn"])<<" \\ "<<endl;

  cout<<" \\hline"<<endl;

  addResults("ge1bLoose","ttbarQCD","Znn","total");
  addResults("ge1bTight","ttbarQCD","Znn","total");
  addResults("ge2bLoose","ttbarQCD","Znn","total");
  addResults("ge2bTight","ttbarQCD","Znn","total");
  addResults("ge3bLoose","ttbarQCD","Znn","total");

  cout<<"Total SM & "
      <<formatLatex( resultsMap_["ge1bLoose"]["total"])<<" & "
      <<formatLatex( resultsMap_["ge1bTight"]["total"])<<" & "
      <<formatLatex( resultsMap_["ge2bLoose"]["total"])<<" & "
      <<formatLatex( resultsMap_["ge2bTight"]["total"])<<" & "
      <<formatLatex( resultsMap_["ge3bLoose"]["total"])<<" \\ "<<endl;
  }

}

void runBinsOfMET() {
  setSearchRegions("METbins");

  //run estimates with full syst
  runDataQCD2011();
  runTtbarEstimate2011(); 

  //now hard-code Keith's results
  //https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=175106
  //page 2
  resultsMap_["ge3bMETBin1"]["Znn"].value = 1.4;
  resultsMap_["ge3bMETBin1"]["Znn"].statError = 0.6;
  resultsMap_["ge3bMETBin1"]["Znn"].systError = 1.5;
  resultsMap_["ge3bMETBin1"]["Znn"].trigErrorPlus = 0;
  resultsMap_["ge3bMETBin1"]["Znn"].trigErrorMinus = 0;

  resultsMap_["ge3bMETBin2"]["Znn"].value = 0.9;
  resultsMap_["ge3bMETBin2"]["Znn"].statError = 0.4;
  resultsMap_["ge3bMETBin2"]["Znn"].systError = 1.0;
  resultsMap_["ge3bMETBin2"]["Znn"].trigErrorPlus = 0;
  resultsMap_["ge3bMETBin2"]["Znn"].trigErrorMinus = 0;

  resultsMap_["ge3bMETBin3"]["Znn"].value = 0.4 + 0.2 + 0.1 + 0.1;
  resultsMap_["ge3bMETBin3"]["Znn"].statError = sqrt( 0.2*0.2 + 0.2*0.2 + 0.1*0.1 + 0.1*0.1);
  resultsMap_["ge3bMETBin3"]["Znn"].systError = 0.4 + 0.2 + 0.1 + 0.1; //the sum in quadrature of 100% errors is not a 100% error, so i'll put in 100% error here by hand
  resultsMap_["ge3bMETBin3"]["Znn"].trigErrorPlus = 0;
  resultsMap_["ge3bMETBin3"]["Znn"].trigErrorMinus = 0;

  writeResultsMapToText("DDresults_METbins3B_1BLSL.dat");

}


void runFineBinsOfMET() {
  setSearchRegions("METfinebins");

  assert(use1B_SL_==true);

  //run estimates with full syst
  runDataQCD2011();
  runTtbarEstimate2011(); 

  //now hard-code Keith's results
  //https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=175106
  //page 2
  resultsMap_["ge3bMETFineBin1"]["Znn"].value = 1.4;
  resultsMap_["ge3bMETFineBin1"]["Znn"].statError = 0.6;
  resultsMap_["ge3bMETFineBin1"]["Znn"].systError = 1.5;
  resultsMap_["ge3bMETFineBin1"]["Znn"].trigErrorPlus = 0;
  resultsMap_["ge3bMETFineBin1"]["Znn"].trigErrorMinus = 0;

  resultsMap_["ge3bMETFineBin2"]["Znn"].value = 0.9;
  resultsMap_["ge3bMETFineBin2"]["Znn"].statError = 0.4;
  resultsMap_["ge3bMETFineBin2"]["Znn"].systError = 1.0;
  resultsMap_["ge3bMETFineBin2"]["Znn"].trigErrorPlus = 0;
  resultsMap_["ge3bMETFineBin2"]["Znn"].trigErrorMinus = 0;

  resultsMap_["ge3bMETFineBin3"]["Znn"].value = 0.4 ;
  resultsMap_["ge3bMETFineBin3"]["Znn"].statError = 0.2;
  resultsMap_["ge3bMETFineBin3"]["Znn"].systError = 0.4 ;
  resultsMap_["ge3bMETFineBin3"]["Znn"].trigErrorPlus = 0;
  resultsMap_["ge3bMETFineBin3"]["Znn"].trigErrorMinus = 0;

  resultsMap_["ge3bMETFineBin4"]["Znn"].value = 0.2;
  resultsMap_["ge3bMETFineBin4"]["Znn"].statError =  0.2;
  resultsMap_["ge3bMETFineBin4"]["Znn"].systError = 0.2;
  resultsMap_["ge3bMETFineBin4"]["Znn"].trigErrorPlus = 0;
  resultsMap_["ge3bMETFineBin4"]["Znn"].trigErrorMinus = 0;

  resultsMap_["ge3bMETFineBin5"]["Znn"].value =  0.1;
  resultsMap_["ge3bMETFineBin5"]["Znn"].statError = 0.1;
  resultsMap_["ge3bMETFineBin5"]["Znn"].systError = 0.1; 
  resultsMap_["ge3bMETFineBin5"]["Znn"].trigErrorPlus = 0;
  resultsMap_["ge3bMETFineBin5"]["Znn"].trigErrorMinus = 0;

  resultsMap_["ge3bMETFineBin6"]["Znn"].value = 0.1;
  resultsMap_["ge3bMETFineBin6"]["Znn"].statError = 0.1;
  resultsMap_["ge3bMETFineBin6"]["Znn"].systError = 0.1; 
  resultsMap_["ge3bMETFineBin6"]["Znn"].trigErrorPlus = 0;
  resultsMap_["ge3bMETFineBin6"]["Znn"].trigErrorMinus = 0;

  writeResultsMapToText("DDresults_METfinebins3B_1BLSL.dat");

}


//read back the file generated above and make a plot
void drawDD() 
{
  gROOT->SetStyle("CMS");
  loadSamples();

  setSearchRegions("METbins");
  ifstream file("DDresults_METbins3B.dat");

  //  setSearchRegions("METfinebins");
  //  ifstream file("DDresults_METfinebins3B_1BLSL.dat");

  int nbins=(int) searchRegions_.size();
  float *metbins = new float[nbins+1];
  for (int ibin=0; ibin<nbins; ibin++) {
    metbins[ibin]=  searchRegions_.at(ibin).getLowEdgeMET();
  }
  metbins[nbins] = 550; //hard-coded because there is no upper bound on MET

  TString xtitle = "E_{T}^{miss}";

  resetHistos();
  setColorScheme("stack");

  if (totalsm!=0) delete totalsm;
  totalsm =  new TH1D("totalsm","",nbins,metbins);
  totalsm->Sumw2();

  if (thestack!= 0 ) delete thestack;
  thestack = new THStack("thestack","--");

  //histos are now filled
  renewCanvas();
  renewLegend();

  TString binName;
  TString backgroundName;
  texData vals;
  while (file>>binName) {
    file>>backgroundName;
    file>>vals.value;
    file>>vals.statError;
    file>>vals.systError;
    file>>vals.trigErrorPlus;
    file>>vals.trigErrorMinus;

    if (!histos_.count(backgroundName)) {
      histos_[backgroundName] = new TH1D(backgroundName,backgroundName,nbins,metbins);
      histos_[backgroundName]->Sumw2();
    }
    //need to get the bin number
    unsigned int index=0;
    for ( ; index<searchRegions_.size(); index++) {
      if  (binName==searchRegions_[index].id()) break;
    }
    histos_[backgroundName]->SetBinContent(index+1, vals.value);
    double avtrigerr = 0.5*(fabs(vals.trigErrorPlus) + fabs(vals.trigErrorMinus));
    double totalerr = sqrt( vals.systError*vals.systError + vals.statError*vals.statError + avtrigerr*avtrigerr);
    histos_[backgroundName]->SetBinError(index+1, totalerr);

    histos_[backgroundName]->SetFillColor(sampleColor_[backgroundName]);
    histos_[backgroundName]->SetMarkerSize(0);
    histos_[backgroundName]->SetXTitle(xtitle);
    //    histos_[backgroundName]->SetYTitle(ytitle);

  }

  for (map<TString, TH1D*>::iterator ih=histos_.begin(); ih!=histos_.end(); ++ih) {

    backgroundName= ih->first;

    if (backgroundName != "ttbarQCD") {
      thestack->Add(histos_[backgroundName]);
    }
      
    if (backgroundName=="ttbarQCD" || backgroundName=="Zinvisible") {
      totalsm ->Add(histos_[backgroundName]);
    }
  }
  //get the legend filled in the opposite order
  for (map<TString, TH1D*>::reverse_iterator ih=histos_.rbegin(); ih!=histos_.rend(); ++ih) {
    backgroundName= ih->first;
    if (backgroundName != "ttbarQCD") leg->AddEntry(histos_[backgroundName], sampleLabel_[backgroundName]);
  }

  totalsm->SetMarkerStyle(sampleMarkerStyle_["totalsm"]);

  thecanvas->cd();
  thestack->Draw("hist");
    thestack->GetHistogram()->GetXaxis()->SetTitle(xtitle);

  //draw errors on DD
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
  //  totalsm->Draw("same");
  leg->Draw();

  gROOT->cd();
  if (hdata!=0) delete hdata;
  hdata =  new TH1D("hdata","",nbins,metbins);
  hdata->Sumw2();
  //to draw data we have to define all of the cuts....
  
  TCut baseline = "cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passCleaning==1 && minDeltaPhiN>=4";
  TCut htcut = searchRegions_[0].htSelection.Data();
  TCut btagging="nbjets>=3";
  TCut thecut = baseline&&htcut&&btagging;
  gROOT->cd();
  dtree->Project("hdata","MET",thecut.GetTitle()); //oy...we should use getCutString
  hdata->SetMarkerColor(kBlack);
  hdata->SetLineWidth(2);
  hdata->SetMarkerStyle(kFullCircle);
  hdata->SetMarkerSize(1);
  hdata->Draw("SAME");

  double mymax = thestack->GetMaximum();
  if ( findOverallMax(totalsm) >mymax) mymax = findOverallMax(totalsm);
  if (findOverallMax(hdata) > mymax) mymax = findOverallMax(hdata);
  thestack->SetMaximum( maxScaleFactor_*mymax);


}

void AN2011_prescale( TString btagselection="ge1b",const int mode=1 ) {

  loadSamples();

  //this mode thing is kludgey
  TString modestring="";
  if (mode==1) {
    usePUweight_=false;
    useHTeff_=false;
    useMHTeff_=false;
    thebnnMHTeffMode_ = kOff;
    btagSFweight_="1";
    currentConfig_=configDescriptions_.getDefault(); //completely raw MC
  }
  else if (mode==2 || mode==3) {
    usePUweight_=true;
    useHTeff_=true;
    useMHTeff_=false;//SPECIAL FOR PRESCALE DATA!
    thebnnMHTeffMode_ = kOff;
    currentConfig_=configDescriptions_.getCorrected(); //JER bias
    if (mode==2) modestring="-JER-PU-HLT";
    else if(mode==3) modestring="-JER-PU-HLT-bSF";
  }
  else assert(0);

  TCut btagcut = "nbjetsCSVM>=1";
  if (mode==1 || mode==2) {
    if ( btagselection=="ge1b") {} //do nothing
    else  if ( btagselection=="ge2b" ) {
      btagcut = "nbjetsCSVM>=2";
    }
    else if ( btagselection=="eq1b" ) {
      btagcut = "nbjetsCSVM==1";
    }
    else if ( btagselection=="eq2b" ) {
      btagcut = "nbjetsCSVM==2";
    }
    else if ( btagselection=="eq0b" ) {
      btagcut = "nbjetsCSVM==0";
    }
    else if (btagselection=="ge3b" ){
      btagcut = "nbjetsCSVM>=3";
    }
    else {
      assert(0);
    }
  }
  else if (mode==3) {
    if ( btagselection=="ge1b") {
      btagcut="1";
      btagSFweight_="probge1";
    }
    else  if ( btagselection=="ge2b" ) {
      btagcut = "1";
      btagSFweight_="probge2";
    }
    else if ( btagselection=="eq0b" ) {
      btagcut = "1";
      btagSFweight_="prob0";
    }
    else if ( btagselection=="eq1b" ) {
      btagcut = "1";
      btagSFweight_="prob1";
    }
    else if ( btagselection=="eq2b" ) {
      btagcut = "1";
      btagSFweight_="(1-prob1-prob0-probge3)";
    }
    else {
      assert(0);
    }
  }



  int nbins;
  float low,high;
  TString var,xtitle;
  
  doOverflowAddition(true);

  setQuiet(false);

  // ==========================

  TCut util = "pass_utilityHLT==1 && weight<1000";
  //lumiScale_ = 19.23288; //customize lumiScale_
  lumiScale_ = 30.15471; 

  // ========= regular N-1 plots

  resetSamples(); //use all samples
  setStackMode(true); //regular stack
  setColorScheme("stack");
  doData(true);
  drawMCErrors_=true;

  
  //MET
  selection_ =TCut("cutHT==1 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4 && passCleaning==1")&&util&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 20; low=0; high=200;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_MET_"+btagselection+modestring,0,"GeV");
  setLogY(true);  setPlotMinimum(1e-1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_MET_"+btagselection+modestring,0,"GeV");
  

  
  //njets in 50 < MET < 100 region
  selection_ =TCut("cutHT==1 && cutPV==1 && MET>=50 && MET<100 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4 && passCleaning==1")&&util&&btagcut;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 7; low=2; high=9;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_njets_LSB_"+btagselection+modestring);
  

  selection_ =TCut("cutHT==1 && cutPV==1 && MET>=50 && MET<100 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passCleaning==1")&&util&&btagcut;
  var="minDeltaPhiN"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 20; low=0; high=40;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_minDeltaPhiN_LSB_"+btagselection+modestring);
  setLogY(true);  setPlotMinimum(1e-1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_minDeltaPhiN_LSB_"+btagselection+modestring);
  return;

  //look again at MET
  selection_ =TCut("MET>=50 && MET<100 && cutHT==1 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4 && passCleaning==1")&&util&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 10; low=50; high=100;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_METlimited_"+btagselection+modestring,0,"GeV");


  //sanity check
/*
  btagSFweight_="1";
  doRatioPlot(true);
  ratioMin=0; ratioMax=3;
  selection_ =TCut("MET>150 && cutHT==1 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN < 4")&&util;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 10; low=150; high=250;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_METhigh_lowDP");

  lumiScale_ =  1143;
  ratioMin=1; ratioMax=1.5;
  selection_ ="cutTrigger==1 && MET>150 && cutHT==1 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN < 4";
  drawPlots(var,nbins,low,high,xtitle,"Events", "unprescaled_METhigh_lowDP");
  */

 
  //// ben
/*
  doRatioPlot(true);
  selection_ =TCut("HT>=400 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1&&passCleaning==1 && MET<100")&&util&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 10; low=0; high=100;
  ratioMin=0; ratioMax=2;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_MET_"+btagselection+modestring,0,"GeV");
*/
/*  
  doRatioPlot(true);
  selection_ =TCut("HT>=400 && cutPV==1 && MET>=50 && MET<100 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&btagcut;
  var="minDeltaPhiN"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 20; low=0; high=40;
  setLogY(false); resetPlotMinimum();
  //drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_minDeltaPhiN_LSB_"+btagselection+modestring);
  setLogY(true);  setPlotMinimum(1e-1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_minDeltaPhiN_LSB_"+btagselection+modestring);
  */

/*
  doRatioPlot(true);
  selection_ =TCut("HT>=400 && cutPV==1 && MET>=50 && MET<100 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN < 4 &&passCleaning==1")&&btagcut&&util;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 7; low=2; high=9;
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_njets_"+btagselection+modestring);
*/
 
}


void AN2011_PUQCD( TString btagselection="ge1b",const int mode=1 ) {
  
  loadSamples();

  
  //this mode thing is kludgey
  TString modestring="";
  if (mode==1) {
    usePUweight_=false;
    useHTeff_=false;
    useMHTeff_=false;
    thebnnMHTeffMode_ = kOff;
    btagSFweight_="1";
    currentConfig_=configDescriptions_.getDefault(); //completely raw MC
  }
  else if (mode==2 || mode==3) {
    usePUweight_=true;
    useHTeff_=true;
    useMHTeff_=false;//SPECIAL for prescaled data!
    thebnnMHTeffMode_ = kOff;
    currentConfig_=configDescriptions_.getCorrected(); //JER bias
    if (mode==2) modestring="-JER-PU-HLT";
    else if(mode==3) modestring="-JER-PU-HLT-bSF";
  }
  else assert(0);
  
  TCut btagcut = "nbjetsCSVM>=1";
  if (mode==1 || mode==2) {
    if ( btagselection=="ge1b") {} //do nothing
    else  if ( btagselection=="ge2b" ) {
      btagcut = "nbjetsCSVM>=2";
    }
    else if ( btagselection=="eq1b" ) {
      btagcut = "nbjetsCSVM==1";
    }
    else if ( btagselection=="eq0b" ) {
      btagcut = "nbjetsCSVM==0";
    }
    else if ( btagselection=="eq2b" ) {
      btagcut = "nbjetsCSVM==2";
    }
    else if (btagselection=="ge3b" ){
      btagcut = "nbjetsCSVM>=3";
    }
    else {
      assert(0);
    }
  }
  else if (mode==3) {
    if ( btagselection=="ge1b") {
      btagcut="1";
      btagSFweight_="probge1";
    }
    else  if ( btagselection=="ge2b" ) {
      btagcut = "1";
      btagSFweight_="probge2";
    }
    else  if ( btagselection=="ge3b" ) {
      btagcut = "1";
      btagSFweight_="probge3";
    }
    else if ( btagselection=="eq1b" ) {
      btagcut = "1";
      btagSFweight_="prob1";
    }
    else if ( btagselection=="eq0b" ) {
      btagcut = "1";
      btagSFweight_="prob0";
    }
    else if ( btagselection=="eq2b" ) {
      btagcut = "1";
      btagSFweight_="(1-prob1-prob0-probge3)";
    }
    else {
      assert(0);
    }
  }
  
  
  int nbins;
  float low,high;
  TString var,xtitle;
  
  doOverflowAddition(true);
  setQuiet(false);

  // ==========================

  TCut util = "pass_utilityHLT==1 && weight<1000";
  double holdLumiScale = lumiScale_;
  lumiScale_ = 30.15471; 

  resetSamples(); //use all samples


  doData(true);
  drawMCErrors_=true;

  setStackMode(false); //regular stack
  setColorScheme("nostack");
  drawTotalSM_=true;  

  setPlotMaximum(0.5); setPlotMinimum(0);
  selection_ =TCut("HT>=400 && MET >=50 && MET <100 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 &&passCleaning==1")&&util&&btagcut;
  drawR("minDeltaPhiN", 4, "nGoodPV", 14, 0.5, 14.5, "nGoodPV_"+btagselection+modestring);
  return;

    
  drawTotalSM_=false;  
  doOverflowAddition(false);
  setPlotMaximum(0.5); setPlotMinimum(0);
  selection_ =TCut("HT>=400 && MET >=50 && MET <100 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passCleaning==1")&&util&&btagcut;
  drawR("minDeltaPhiN", 4, "runNumber", 10, 160300, 180300, "runNumber_"+btagselection+modestring,0,true);
  
  
  lumiScale_ = holdLumiScale;
  selection_ =TCut("HT>=400 && MET >=50 && MET <100 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passCleaning==1")&&util&&btagcut;
  drawSimple("nGoodPV",14,0.5,14.5,"dummy", "","data");
  TH1D* PVpres = (TH1D*)hinteractive->Clone("PVpres");
  selection_ =TCut("HT>=400 && MET >=200 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passCleaning==1 && cutTrigger==1")&&btagcut;
  drawSimple("nGoodPV",14,0.5,14.5,"dummy", "","data");
  TH1D* PVphys = (TH1D*)hinteractive->Clone("PVphys");
  TCanvas * Cpv = new TCanvas("myC", "myC", mainpadWidth, mainpadHeight);
  Cpv->cd();
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.07);
  PVpres->SetLineColor(kRed);
  PVpres->SetMarkerColor(kRed);
  PVpres->GetXaxis()->SetTitle("nGoodPV");
  PVpres->GetYaxis()->SetTitle("Arbitrary Units");
  PVpres->DrawNormalized();
  PVphys->DrawNormalized("SAME");
  TLegend * myLegend = new TLegend(0.53,0.9,0.81,0.83); 
  myLegend->SetFillColor(0);
  myLegend->SetBorderSize(0);
  myLegend->SetLineStyle(0);
  myLegend->SetTextFont(42);
  myLegend->SetFillStyle(0);
  myLegend->SetTextSize(0.035);
  myLegend->AddEntry(PVpres,"Data, prescaled triggers", "p");
  myLegend->AddEntry(PVphys, "Data, physics triggers", "p");
  myLegend->Draw();
  drawPlotHeader();
  Cpv->Print("PV_phys_pres_"+btagselection + modestring +".pdf");
  Cpv->Print("PV_phys_pres_"+btagselection + modestring +".png");
  

}

void AN2011_ttbarw( double& r_sl, double& r_sl_err, double& r_nom, double& r_nom_err, 
		    TString btagselection="ge1b", TString HTselection="Loose" , TString samplename="TTbarJets", TString ttbarcomponent="All") {

  /*
goal:

show MET distributions for w+jets, ttbar, singletop
in both the SL and nominal samples

this is going to require some serious kludgey stuff because this code has
never been used to plot the same sample with different cuts on top of each 
other.
  */

  /*
.L drawReducedTrees.C++
  */

  TCut btagcut = "nbjetsCSVM>=1";
  if ( btagselection=="ge1b") {} //do nothing
  else  if ( btagselection=="ge2b" ) {
    btagcut = "nbjetsCSVM>=2";
  }
  else if ( btagselection=="eq1b" ) {
    btagcut = "nbjetsCSVM==1";
  }
  else if ( btagselection=="eq2b" ) {
    btagcut = "nbjetsCSVM==2";
  }
  else if ( btagselection=="ge0b" ) {
    btagcut = "nbjetsCSVM>=0";
  }
  else if ( btagselection=="ge3b" ) {
    btagcut = "nbjetsCSVM>=3";
  }
  else {
    assert(0);
  }


  TCut HTcut="HT>=400";
  if (btagselection=="ge1b" && HTselection=="Tight")  HTcut="HT>=500";
  if (btagselection=="ge2b" && HTselection=="Tight")  HTcut="HT>=600";

  loadSamples();

  if (useScaleFactors_) {
    usePUweight_=true;
    useHTeff_=true;
    useMHTeff_=false;
    thebnnMHTeffMode_ = kOff;
    currentConfig_=configDescriptions_.getCorrected(); //add JERbias

    ////The MET comparison plots for the AN have the b-tag SF stuff below commented out!
    //btagcut="1";   
    //if (btagselection=="ge2b") {
    //  btagSFweight_="probge2";
    //}
    //else if (btagselection=="ge1b") {
    //  btagSFweight_="probge1";
    //}
    //else if (btagselection=="ge3b") {
    //  btagSFweight_="probge3";
    //}
    //else {assert(0);}
  }


  int nbins;
  float low,high;
  TString var,xtitle;
  bool vb = true;

  if(samplename=="datamode"){

    if(useScaleFactors_){ 
      useMHTeff_=true;
      thebnnMHTeffMode_ = kOff;

      btagcut="1";   
      if (btagselection=="ge2b") {
	btagSFweight_="probge2";
      }
      else if (btagselection=="ge1b") {
	btagSFweight_="probge1";
      }
      else if (btagselection=="eq1b") {
	btagSFweight_="prob1";
      }
      else if (btagselection=="eq2b") {
	btagSFweight_="(1-prob1-prob0-probge3)";
      }
      else if (btagselection=="ge3b") {
	btagSFweight_="probge3";
      }
      else {assert(0);}
    }

    resetSamples(); //use all samples
    setStackMode(true); //regular stack
    setColorScheme("stack");
    doData(true);
    drawMCErrors_=true;
    doOverflowAddition(false);

    //ST
    //selection_ =TCut("MET>=200 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && minDeltaPhiN >= 4 &&cutEleVeto==1 && cutMuVeto==1 && nTaus==1 &&passCleaning==1")&&btagcut&&HTcut;;
    //SL
    selection_ =TCut("MET>=200 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 &&passCleaning==1")&&btagcut&&HTcut;;
    //SE
    //selection_ =TCut("MET>=200 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && (nElectrons==1 && nMuons==0) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 &&passCleaning==1")&&btagcut&&HTcut;;
    //SMU
    //selection_ =TCut("MET>=200 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && (nElectrons==0 && nMuons==1) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 &&passCleaning==1")&&btagcut&&HTcut;;
    var="MET"; xtitle="E_{T}^{miss} [GeV]";
    nbins = 1; low=200; high=250;
    //nbins = 15; low=200; high=500;
    drawPlots(var,nbins,low,high,xtitle,"Events", "test_"+btagselection);

    double dataSB = dataIntegral_; 
    double MCSB = MCIntegral_, MCSBErr = MCIntegralErr_;

    doOverflowAddition(true);
    nbins = 1; low=250; high=600;
    if (HTselection=="Tight" && btagselection=="ge1b")  low=500;
    else if (HTselection=="Tight" && btagselection=="ge2b")  low=300;
    drawPlots(var,nbins,low,high,xtitle,"Events", "test_"+btagselection);

    double dataSIG = dataIntegral_; 
    double MCSIG = MCIntegral_, MCSIGErr = MCIntegralErr_;

    std::cout << "dataSB = " << dataSB << ", MCSB = " << MCSB << std::endl;
    std::cout << "dataSIG = " << dataSIG << ", MCSIG = " << MCSIG << std::endl;

    double dataRatioSL = dataSIG/dataSB;
    double dataRatioSL_err = jmt::errAoverB(dataSIG,sqrt(dataSIG),dataSB,sqrt(dataSB));

    double MCRatioSL = MCSIG/MCSB;
    double MCRatioSL_err = jmt::errAoverB(MCSIG, MCSIGErr ,MCSB, MCSBErr);

    std::cout << btagselection << "," << HTselection << ":" << "R_SL(data) = " << dataRatioSL << " +/- " << dataRatioSL_err 
	      << ", R_SL(MC) = " << MCRatioSL << " +/- " << MCRatioSL_err << std::endl;

  }
  else{
    //ge1b, Loose
    //const int nvarbins=18;
    //const float varbins[]={100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 550}; //AN and PAS

    //ge1b, Tight
    //const int nvarbins=15;
    //const float varbins[]={100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 350, 400, 450, 500, 550}; 
  
    //ge2b, Loose and Tight
    // ===== This is the one Don has in CVS!
    //    const int nvarbins=15;
    //    const float varbins[]={100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, 600, 650}; //AN and PAS 
    const int nvarbins=11;
    const float varbins[]={100, 125, 150, 175, 200, 225, 250,  300, 350, 400,  500,  650}; //AN and PAS 
    //const int nvarbins=12;
    //const float varbins[]={200, 225, 250, 275, 300, 350, 400, 450, 550,650,750, 850,950}; 

    doOverflowAddition(true);

    setStackMode(false,true); //normalized
    setColorScheme("nostack");
    clearSamples();
    addSample(samplename);
    drawLegend(false);
    doRatioPlot(false);

    doData(false);

    TString sample;
    if(samplename=="TTbarJets") sample = "ttbar";
    else if(samplename=="WJets") sample = "wjets";
    else if(samplename=="SingleTop") sample = "singletop";

    const bool PVhack=false;

    if (PVhack)  selection_ =TCut("cutPV==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==0)||(nElectrons==0 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=-2 && MT_Wlep<100 && nGoodPV>=10 && passCleaning==1")&&btagcut&&HTcut;
    else         selection_ =TCut("cutPV==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 && passCleaning==1")&&btagcut&&HTcut;
    //selection_ =TCut("MET>=150 && cutPV==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4")&&btagcut&&HTcut; //AN v5
    var="MET"; xtitle="E_{T}^{miss} [GeV]";
    //nbins = 40; low=100; high=550; 
    nbins = 55; low=100; high=650; 
    //nbins = 40; low=200; high=600; 
    //nbins = 40; low=150; high=550; //AN v5
    //nbins = 20; low=150; high=550;
    //nbins = 16; low=150; high=550;
    if(vb) drawPlots(var,nvarbins, varbins, xtitle,"Arbitrary units", "SBandSIG_MET_SL_"+sample+"_"+btagselection+"_"+HTselection);
    else drawPlots(var,nbins,low,high,xtitle,"Arbitrary units", "SBandSIG_MET_SL_"+sample+"_"+btagselection+"_"+HTselection);

    ////////////////////////////////////////////////
    //This part sets up the SIG/SB ratio computation
    ////////////////////////////////////////////////

    //the lower boundary on the SB should be set to 200 GeV
    int lowbin = hinteractive->FindBin(200);
    int boundarybin_sb = 0, boundarybin_sig=0;
    //the upper boundary on the SB should be set to 250 GeV
    boundarybin_sb = hinteractive->FindBin(250);

    //the lower boundary on the SR depends on the selection
    if (HTselection=="Tight" && btagselection=="ge1b")  boundarybin_sig = hinteractive->FindBin(500);
    else if (HTselection=="Tight" && btagselection=="ge2b")  boundarybin_sig = hinteractive->FindBin(300);
    else boundarybin_sig = hinteractive->FindBin(250);

    //the upper boundary on the SR region should be the last bin on the plot
    int highbin= hinteractive->GetNbinsX();


    ////////////////////////////////////////
    //This part computes the SL SIG/SB ratio
    ////////////////////////////////////////

    double sl_sb_err=0, sl_sig_err=0;
    double sl_sb = hinteractive->IntegralAndError(lowbin,boundarybin_sb-1,sl_sb_err);//integral from 200 to <250
    double sl_sig = hinteractive->IntegralAndError(boundarybin_sig,highbin,sl_sig_err);//integral from e.g. 300 to last bin
    double r_sl_sigoversb = sl_sig/sl_sb;
    double r_sl_sigoversb_err = jmt::errAoverB(sl_sig,sl_sig_err,sl_sb,sl_sb_err);

    std::cout << "SL: SB content= " << sl_sb << " +/- " << sl_sb_err << std::endl;
    std::cout << "SL: SIG content= " << sl_sig << " +/- " << sl_sig_err << std::endl;
    std::cout << "SL: SIG/SB ratio = " << r_sl_sigoversb << " +/- " << r_sl_sigoversb_err << std::endl;

    r_sl = r_sl_sigoversb;
    r_sl_err = r_sl_sigoversb_err;
  
    TH1D* SLplot = (TH1D*)hinteractive->Clone("SLplot");


    //now switch to the normal selection
    //CAREFUL OF THIS HACK!
    if (PVhack)  selection_ =TCut("cutPV==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==0)||(nElectrons==0 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=-2 && MT_Wlep<100 && passCleaning==1 && nGoodPV<10")&&btagcut&&HTcut;
    else{    

      if(ttbarcomponent=="All"){
	selection_ =TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4 && passCleaning==1 ")&&btagcut&&HTcut;
      }
      else{
	if(!splitTTbar_){std::cout << "ERROR: For mode " << ttbarcomponent << " the splitTTbar_ flag must be true!" << std::endl; assert(0);}
	clearSamples();
	if(ttbarcomponent=="TauHad"){
	  addSample("TTbarJets-semiTauHad");
	  selection_ =TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4 && passCleaning==1 ")&&btagcut&&HTcut;
	}
         
	else if(ttbarcomponent=="SemiEleFailEta"){
	  addSample("TTbarJets-semiEle");
	  selection_ =TCut(" cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN>=4 && passCleaning==1 && (decayType==101101 || decayType==201101)")&&btagcut&&HTcut;
	}
	else if(ttbarcomponent=="SemiEleFailPt"){
	  addSample("TTbarJets-semiEle");
	  selection_ =TCut(" cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN>=4 && passCleaning==1 && (decayType==101102 || decayType==201102)")&&btagcut&&HTcut;
	}
	else if(ttbarcomponent=="SemiEleFailOther"){
	  addSample("TTbarJets-semiEle");
	  selection_ =TCut(" cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN>=4 && passCleaning==1 && decayType==101104")&&btagcut&&HTcut;
	}
	
	else if(ttbarcomponent=="SemiMuFailEta"){
	  addSample("TTbarJets-semiMu");
	  selection_ =TCut(" cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN>=4 && passCleaning==1 && (decayType==101301 || decayType==201301)")&&btagcut&&HTcut;
	}
	else if(ttbarcomponent=="SemiMuFailPt"){
	  addSample("TTbarJets-semiMu");
	  selection_ =TCut(" cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN>=4 && passCleaning==1 && (decayType==101302 || decayType==201302)")&&btagcut&&HTcut;
	}
	else if(ttbarcomponent=="SemiMuFailOther"){
	  addSample("TTbarJets-semiMu");
	  selection_ =TCut(" cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN>=4 && passCleaning==1 && decayType==101304")&&btagcut&&HTcut;
	}
	else{std::cout << "ttbarcomponent not recognized.  Quitting" << std::endl; assert(0);}

      }
    }
    
    //selection_ =TCut("MET>=150 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&btagcut&&HTcut;// AN v5
    if(vb)drawPlots(var,nvarbins, varbins ,xtitle,"Arbitrary units", "SBandSIG_MET_normal_"+sample+"_"+btagselection+"_"+HTselection);
    else drawPlots(var,nbins,low,high,xtitle,"Arbitrary units", "SBandSIG_MET_normal_"+sample+"_"+btagselection+"_"+HTselection);
    TH1D* SIGplot = (TH1D*)hinteractive->Clone("SIGplot");
    SLplot->SetLineColor(kRed);
    SLplot->SetMarkerColor(kRed);
    SLplot->Draw("SAME");
    
    if (PVhack)  thecanvas->SaveAs("METshape_"+sample+"_0LinPVbins_"+btagselection+"_"+HTselection+".pdf");
    else         thecanvas->SaveAs("METshape_"+sample+"_SLandStandard_"+btagselection+"_"+HTselection+".pdf");

  
    ////////////////////////////////////////
    //This part computes the nominal SIG/SB ratio
    ////////////////////////////////////////

    double lv_sb_err=0, lv_sig_err=0;
    double lv_sb = hinteractive->IntegralAndError(lowbin,boundarybin_sb-1,lv_sb_err);
    double lv_sig = hinteractive->IntegralAndError(boundarybin_sig,highbin,lv_sig_err);
    double r_lv_sigoversb = lv_sig/lv_sb;
    double r_lv_sigoversb_err = jmt::errAoverB(lv_sig,lv_sig_err,lv_sb,lv_sb_err);

    std::cout << "LV: SB content= " << lv_sb << " +/- " << lv_sb_err << std::endl;
    std::cout << "LV: SIG content= " << lv_sig << " +/- " << lv_sig_err << std::endl;
    std::cout << "LV: SIG/SB ratio = " << r_lv_sigoversb << " +/- " << r_lv_sigoversb_err << std::endl;

    r_nom = r_lv_sigoversb;
    r_nom_err = r_lv_sigoversb_err;

    //Hack to get it plotted with ratio plot
    TCanvas* myC = 0;
    myC = new TCanvas("myC", "myC", 600,700);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    Float_t small = 1e-5;
    myC->Divide(1,2,small,small);
    //myC->Divide(1,2);
    //const float padding=0.01; const float ydivide=0.2;
    const float padding=1e-5; const float ydivide=0.3;
    myC->GetPad(1)->SetPad( padding, ydivide + padding, 1-padding, 1-padding);
    myC->GetPad(2)->SetPad( padding, padding, 1-padding, ydivide-padding);
    myC->GetPad(1)->SetLogy(1);
    myC->GetPad(1)->SetRightMargin(.05);
    myC->GetPad(2)->SetRightMargin(.05);
    myC->GetPad(2)->SetBottomMargin(.4);
    myC->GetPad(1)->Modified();
    myC->GetPad(2)->Modified();
    myC->cd(1);
    gPad->SetBottomMargin(small);
    gPad->Modified();

    renormBins(SIGplot,-1);
    renormBins(SLplot,-1);
    SIGplot->GetYaxis()->SetTitleSize(0.07);
    SIGplot->GetYaxis()->SetTitleOffset(1.);
    SIGplot->Draw("HIST E");
    SLplot->SetMarkerStyle(kFullTriangleUp);
    SLplot->Draw("HIST E SAME");
    TH1D* myRatio;
    if(vb) myRatio = new TH1D("ratio", "ratio", nvarbins, varbins);
    else  myRatio = new TH1D("ratio", "ratio", nbins,low,high);
    myRatio->Sumw2();
    myRatio->SetLineColor(kBlack);
    myRatio->SetMarkerColor(kBlack);
    myRatio->Divide(SIGplot,SLplot);
    myRatio->SetMinimum(0);
    myRatio->SetMaximum(2.5);
    myRatio->GetYaxis()->SetNdivisions(200 + int(ratioMax-ratioMin)+1);    //set ticks ; to be seen if this really works
    //myRatio->GetYaxis()->SetNdivisions(400);
    myRatio->GetYaxis()->SetLabelSize(0.1); //make y label bigger
    myRatio->GetXaxis()->SetLabelSize(0.1); //make y label bigger
    myRatio->GetXaxis()->SetTitleOffset(1.1);
    myRatio->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]"); //make y label bigger
    myRatio->GetXaxis()->SetLabelSize(0.12);
    myRatio->GetXaxis()->SetLabelOffset(0.04);
    myRatio->GetXaxis()->SetTitleSize(0.16);
    myRatio->GetYaxis()->SetTitle("Ratio ");
    myRatio->GetYaxis()->SetTitleSize(0.16);
    myRatio->GetYaxis()->SetTitleOffset(.43);
    myC->cd(2);
    gPad->SetTopMargin(small);
    gPad->SetTickx();
    gPad->Modified();

    myRatio->Draw();
    TLine* myLine;
    myLine = new TLine(low, 1, high, 1);
    myLine->Draw();
    myC->cd(1);
    //TLatex* mytext = new TLatex(3.570061,23.08044,"CMS Preliminary");
    TLatex* mytext = new TLatex(3.570061,23.08044,"CMS Simulation");
    mytext->SetNDC();
    mytext->SetTextAlign(13);
    mytext->SetX(0.6);
    mytext->SetY(0.9);
    mytext->SetTextFont(42);
    mytext->SetTextSizePixels(24);
    mytext->Draw();
    TLegend* myLegend = new TLegend(0.615,0.84,0.81,0.7);
    myLegend->SetFillColor(0);
    myLegend->SetBorderSize(0);
    myLegend->SetLineStyle(0);
    myLegend->SetTextFont(42);
    myLegend->SetFillStyle(0);
    myLegend->SetTextSize(0.045);
    if(sample=="ttbar"){
      if (PVhack) {
	myLegend->AddEntry(SIGplot,"t#bar{t}, 0 lepton (n_{PV} #leq 9)", "lp");
	myLegend->AddEntry(SLplot, "t#bar{t}, 0 lepton (n_{PV} #geq 10)", "lp");
      }
      else {
	
	if(ttbarcomponent=="All"){
	  myLegend->AddEntry(SIGplot,"t#bar{t}, 0 leptons", "lp");
	}
	else{
	  if(ttbarcomponent=="TauHad"){
	    myLegend->AddEntry(SIGplot,"t#bar{t}, 0 leptons - semi-#tau-had", "lp");
	  }
	  else if(ttbarcomponent=="SemiEleFailEta" || ttbarcomponent=="SemiEleFailPt" || ttbarcomponent=="SemiEleFailOther"){
	    myLegend->AddEntry(SIGplot,"t#bar{t}, 0 leptons - semi-e", "lp");
	  }

	  else if(ttbarcomponent=="SemiMuFailEta" || ttbarcomponent=="SemiMuFailPt" || ttbarcomponent=="SemiMuFailOther"){
	    myLegend->AddEntry(SIGplot,"t#bar{t}, 0 leptons - semi-#mu", "lp");
	  }
	  
	}
	
	myLegend->AddEntry(SLplot, "t#bar{t}, 1 lepton", "lp");
      }
    }
    else if(sample=="wjets"){
      myLegend->AddEntry(SIGplot,"wjets, 0 leptons", "lp");
      myLegend->AddEntry(SLplot, "wjets, 1 lepton", "lp");
    }
    else if(sample=="singletop"){
      myLegend->AddEntry(SIGplot,"single-top, 0 leptons", "lp");
      myLegend->AddEntry(SLplot, "single-top, 1 lepton", "lp");
    }
    myLegend->Draw();
    if (PVhack)    myC->Print("METshape_logAndRatio_"+sample+"_0LinPVbins_"+btagselection+"_"+HTselection+".pdf");
    //else     myC->Print("METshape_logAndRatio_"+sample+"_SLandStandard_"+btagselection+"_"+HTselection+".pdf");
    else{
      if(ttbarcomponent=="All"){
	myC->Print("METshape_logAndRatio_"+sample+"_SLandStandard_"+btagselection+"_"+HTselection+".pdf");	
      }
      else{
	myC->Print("METshape_Fall11_"+ttbarcomponent+"_logAndRatio_"+sample+"_SLandStandard_"+btagselection+"_"+HTselection+".pdf");	
      }
    }
  }
}

//kludgy , but it suffices for now
void makeTTbarWMETPlots(bool forttbarbreakdown=false) {
  
  double ratio_sl = 0, ratio_sl_err = 0, ratio_nom = 0, ratio_nom_err = 0;
  std::vector<double> ratios_sl, ratios_sl_err, ratios_nom, ratios_nom_err;
  std::vector<string> selections;
  

  if(forttbarbreakdown){


    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Loose" , "TTbarJets", "TauHad");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Tight" , "TTbarJets", "TauHad");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Loose" , "TTbarJets", "TauHad");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Tight" , "TTbarJets", "TauHad");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge3b", "Loose" , "TTbarJets", "TauHad");
    
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Loose" , "TTbarJets", "SemiEleFailEta");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Tight" , "TTbarJets", "SemiEleFailEta");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Loose" , "TTbarJets", "SemiEleFailEta");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Tight" , "TTbarJets", "SemiEleFailEta");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge3b", "Loose" , "TTbarJets", "SemiEleFailEta");
    
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Loose" , "TTbarJets", "SemiEleFailPt");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Tight" , "TTbarJets", "SemiEleFailPt");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Loose" , "TTbarJets", "SemiEleFailPt");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Tight" , "TTbarJets", "SemiEleFailPt");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge3b", "Loose" , "TTbarJets", "SemiEleFailPt");
    
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Loose" , "TTbarJets", "SemiEleFailOther");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Tight" , "TTbarJets", "SemiEleFailOther");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Loose" , "TTbarJets", "SemiEleFailOther");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Tight" , "TTbarJets", "SemiEleFailOther");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge3b", "Loose" , "TTbarJets", "SemiEleFailOther");
    
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Loose" , "TTbarJets", "SemiMuFailEta");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Tight" , "TTbarJets", "SemiMuFailEta");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Loose" , "TTbarJets", "SemiMuFailEta");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Tight" , "TTbarJets", "SemiMuFailEta");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge3b", "Loose" , "TTbarJets", "SemiMuFailEta");
    
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Loose" , "TTbarJets", "SemiMuFailPt");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Tight" , "TTbarJets", "SemiMuFailPt");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Loose" , "TTbarJets", "SemiMuFailPt");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Tight" , "TTbarJets", "SemiMuFailPt");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge3b", "Loose" , "TTbarJets", "SemiMuFailPt");
    
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Loose" , "TTbarJets", "SemiMuFailOther");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Tight" , "TTbarJets", "SemiMuFailOther");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Loose" , "TTbarJets", "SemiMuFailOther");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Tight" , "TTbarJets", "SemiMuFailOther");
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge3b", "Loose" , "TTbarJets", "SemiMuFailOther");

  }
  else{
    std::cout << "TTbarJets" << std::endl;

    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Loose" , "TTbarJets");
    selections.push_back("1BL");
    ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Tight" , "TTbarJets");
    selections.push_back("1BT");
    ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Loose" , "TTbarJets");
    selections.push_back("2BL");
    ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Tight" , "TTbarJets");
    selections.push_back("2BT");
    ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge3b", "Loose" , "TTbarJets");
    selections.push_back("3B");
    ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);
  


    std::cout << "WJets" << std::endl;

    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Loose" , "WJets");
    selections.push_back("1BL");
    ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Tight" , "WJets");
    selections.push_back("1BT");
    ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Loose" , "WJets");
    selections.push_back("2BL");
    ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Tight" , "WJets");
    selections.push_back("2BT");
    ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge3b", "Loose" , "WJets");
    selections.push_back("3B");
    ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

    std::cout << "SingleTop" << std::endl;

    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Loose" , "SingleTop");
    selections.push_back("1BL");
    ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Tight" , "SingleTop");
    selections.push_back("1BT");
    ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Loose" , "SingleTop");
    selections.push_back("2BL");
    ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Tight" , "SingleTop");
    selections.push_back("2BT");
    ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge3b", "Loose" , "SingleTop");
    selections.push_back("3B");
    ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);


    for(uint i = 0; i<selections.size(); ++i){
      if(i==0) std::cout << "TTbarJets" << std::endl;
      if(i==5) std::cout << "WJets" << std::endl;
      if(i==10) std::cout << "SingleTop" << std::endl;
      std::cout << " & "<< selections.at(i) <<" & " << setprecision(3) << ratios_nom.at(i) << "$\\pm$" << ratios_nom_err.at(i) << "\t & " << ratios_sl.at(i) << "$\\pm$" << ratios_sl_err.at(i) << " \\\\" << std::endl; 
    }
  
    //print out the SIG/SB ratio in the SL region for data and full MC
    savePlots_ = false;
    TString mode = "datamode";
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Loose" , mode);
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Tight" , mode);
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Loose" , mode);
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Tight" , mode);
    AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge3b", "Loose" , mode);
  
  }

}



void AN2011_r(TString btagselection="ge1b", const int mode=1) { //for paper, data+MC plots at low MET

  loadSamples();

  setColorScheme("nostack");
  setStackMode(false);
  doOverflowAddition(true);
  drawTotalSM_=true; 
  doRatioPlot(false);
  drawLegend(true);
  doData(true);
  
  TString modestring="";
  if (mode==1) {
    usePUweight_=false;
    useHTeff_=false;
    useMHTeff_=false;
    thebnnMHTeffMode_ = kOff;
    btagSFweight_="1";
    currentConfig_=configDescriptions_.getDefault(); //completely raw MC
  }
  else if (mode==2 || mode==3) {
    usePUweight_=true;
    useHTeff_=true;
    useMHTeff_=false; //SPECIAL for prescaled data!
    thebnnMHTeffMode_ = kOff;
    currentConfig_=configDescriptions_.getCorrected(); //JER bias
    if (mode==2) modestring="-JER-PU-HLT";
    else if(mode==3) modestring="-JER-PU-HLT-bSF";
  }
  else assert(0);
  
  TCut btagcut = "nbjetsCSVM>=1";
  if (mode==1 || mode==2) {
    if ( btagselection=="ge1b") {} //do nothing
    else  if ( btagselection=="ge2b" ) {
      btagcut = "nbjetsCSVM>=2";
    }
    else if ( btagselection=="eq1b" ) {
      btagcut = "nbjetsCSVM==1";
    }  
    else if ( btagselection=="eq0b" ) {
      btagcut = "nbjetsCSVM==0";
    }
    else if ( btagselection=="ge3b" ) {
      btagcut = "nbjetsCSVM>=3";
    }
    else {
      assert(0);
    }
  }
  else if (mode==3) {
    if ( btagselection=="eq0b") {
      btagcut="1";
      btagSFweight_="prob0";
    }
    else if ( btagselection=="ge1b") {
      btagcut="1";
      btagSFweight_="probge1";
    }
    else  if ( btagselection=="ge2b" ) {
      btagcut = "1";
      btagSFweight_="probge2";
    }
    else if ( btagselection=="ge3b" ) {
      btagcut = "1";
      btagSFweight_="probge3";
    }
    else {
      assert(0);
    }
  }
  
  lumiScale_ = 30.15471; //customize lumiScale_
  resetSamples();
  setPlotMinimum(0); setPlotMaximum(0.3);
  leg_x1 = 0.2380496; leg_x2=0.481167; leg_y1=0.5029892; leg_y2=0.9231433;
  
  TCut util = "pass_utilityHLT==1 && weight<2"; //don't use weight<1 when plotting data!

  selection_ =TCut("HT>=400 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passCleaning==1")&&util&&btagcut;
  drawR("minDeltaPhiN",4,15,0,150,"Data_HT400_"+btagselection+modestring);
  
  selection_ =TCut("HT>=500 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passCleaning==1")&&util&&btagcut;
  drawR("minDeltaPhiN",4,15,0,150,"Data_HT500_"+btagselection+modestring); 

  selection_ =TCut("HT>=600 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passCleaning==1")&&util&&btagcut;
  drawR("minDeltaPhiN",4,15,0,150,"Data_HT600_"+btagselection+modestring); 

}


void AN2011_r_SLreq() { //for paper, MC only plots all the way to high MET
  //used for PAS
  loadSamples();
  
  // == draw r(MET)
  clearSamples();
  addSample("PythiaPUQCD");
  setColorScheme("nostack");
  
  usePUweight_=true;

  drawLegend(false);
  doRatioPlot(false);
  doData(false);
  setStackMode(false);
  doOverflowAddition(true);
  
  setPlotMaximum(1); setPlotMinimum(0);
 
  const int nvarbins=15;
  // const float varbins[]={0,15,30,45,80,110,140,170,200,230,260,290,320,350}; //PAS binning //13 bins
   //const float varbins[]={0,15,30,45,60,90,120,150,180,210,240,270,300,400,500,600}; //15 bins
  //  const float varbins[]={0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,450,500,600}; //17 bins
  const float varbins[]={0,20,40,60,80,100,120,140,160,180,200,250,300,400,500,600}; //15 bins
  //  const float varbins[]={0,10,20,30,40,50,60,70,100,130,160,200,240,280,320,360,400,450,500,550,600}; //for pub //20 bins

  //  const float varbins[]={0,20,40,60,80,100,120,140,160,180,200,225,250,275,300,350,450,550,650}; //18 bins

  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjets==0&&weight<1 &&passCleaning==1";
  //   selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 &&weight<1000 &&passCleaning==1";
  //  btagSFweight_="prob0";
  drawR("minDeltaPhiN",4, "MET", nvarbins, varbins,"eq0b");
  
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjets>=1&&weight<1 &&passCleaning==1";
  drawR("minDeltaPhiN",4, "MET", nvarbins, varbins,"ge1b");
  
  setPlotMaximum(3); setPlotMinimum(0);
  drawR("minDeltaPhi",0.3, "MET", nvarbins, varbins,"old_ge1b");

  setPlotMaximum(1); setPlotMinimum(0);
    selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjets>=2&&weight<1 &&passCleaning==1";
  drawR("minDeltaPhiN",4, "MET", nvarbins, varbins,"ge2b");
 
  /*
  //other plots looked at, but not in PAS, for example.
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjets==1&&weight<1000";
  drawR("minDeltaPhiN",4, "MET", nvarbins, varbins,"eq1b");
  
  const int nvarbins2=8;
  const float varbins2[]={0,15,30,45,75,125,175,225, 350};
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjets>=2&&weight<1000";
  drawR("minDeltaPhiN",4, "MET", nvarbins2, varbins2,"ge2b");
  */
}


void AN2011_opt() {

  TCut ge1b = "nbjets>=1";
  TCut ge2b = "nbjets>=2";
  TCut ge3b = "nbjets>=3";

  lumiScale_=1800;

  int nbins;
  float low,high;
  TString var,xtitle;
  
  resetSamples(); //use all samples
  setStackMode(true); //regular stack
  setColorScheme("stack");
  doData(false);
  drawMCErrors_=true;

  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=400 && minDeltaPhiN >= 4")&&ge1b;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 10; low=400; high=1400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "opt_HT_MET400_HT400_ge1b");


  selection_ =TCut("HT>=500 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=300 && minDeltaPhiN >= 4")&&ge2b;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 10; low=400; high=1400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "opt_HT_MET300_HT500_ge2b");


  selection_ =TCut("HT>=900 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN >= 4")&&ge1b;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 10; low=400; high=1400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "opt_HT_MET200_HT900_ge1b");



  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN >= 4")&&ge3b;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 10; low=400; high=1400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "opt_HT_MET200_HT400_ge3b");


  //


  selection_ =TCut("HT>=1000 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN >= 4")&&ge1b;
  var="MET"; xtitle="MET (GeV)";
  nbins = 10; low=200; high=600;
  drawPlots(var,nbins,low,high,xtitle,"Events", "opt_HT_MET200_HT1000_ge1b");

}

void drawMSugraTest() {
  /*
.L drawReducedTrees.C++
  */

  loadSamples();
  clearSamples();
  addSample("TTbarJets");
  addSample("mSUGRAtanb40");
  m0_=800;
  m12_=200;
  doOverflowAddition(true);
  setStackMode(true); //regular stack
  setColorScheme("stack");
  doData(false);
  int nbins; float low,high;
  TString var,xtitle;

  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN >= 4 && nbjets>=1");
  var="njets"; xtitle="Jet multiplicity";
  nbins = 8; low=1; high=9;
  drawPlots(var,nbins,low,high,xtitle,"Events", "testmsugra");

}

void AN2011_excess() {

  //these plots certainly duplicate what we already have, but what the hell

  loadSamples();
  doRatioPlot(true);

  //set all corrections to MC
  usePUweight_=true;
  useHTeff_=true;
  useMHTeff_=true;
  thebnnMHTeffMode_=kOff;
  currentConfig_=configDescriptions_.getCorrected(); //JER bias
  TString modestring="-JER-PU-HLT-bSF";

  setStackMode(true); //regular stack
  setColorScheme("stack");
  doData(true);


  int nbins;
  float low,high;
  TString var,xtitle;
  
  doOverflowAddition(true);
  drawMCErrors_=true;
  
/*
  btagSFweight_="probge2";
  selection_ =TCut("HT>=600 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 &&passCleaning==1 && minDeltaPhiN>=4");
  var="MET"; xtitle="E_{T}^{miss}";
  low=200; high = 500; nbins=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_2BT_MET"+modestring,0,"GeV");
  setLogY(true); setPlotMinimum(0.7);
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_2BT_MET"+modestring,0,"GeV");
  setLogY(false); resetPlotMinimum();

  btagSFweight_="probge3";
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 &&passCleaning==1 && minDeltaPhiN>=4");
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_3B_MET"+modestring,0,"GeV");
  setLogY(true); setPlotMinimum(0.7);
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_3B_MET"+modestring,0,"GeV");
  setLogY(false); resetPlotMinimum();


  btagSFweight_="probge2";
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=300 &&passCleaning==1 && minDeltaPhiN>=4");
  var="HT"; xtitle="H_{T}";
  low=400; high = 1500; nbins=11;
  drawPlots(var,nbins,low,high,xtitle,"Events", "HT_2BT"+modestring,0,"GeV");
  setLogY(true); setPlotMinimum(0.7);
  drawPlots(var,nbins,low,high,xtitle,"Events", "HT_2BT"+modestring,0,"GeV");
  setLogY(false); resetPlotMinimum();
*/
  btagSFweight_="probge3";
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=250 &&passCleaning==1 && minDeltaPhiN>=4");

/*
  var="HT"; xtitle="H_{T}";
  low=400; high = 1400; nbins=5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "HT_3B"+modestring,0,"GeV");
  setLogY(true); setPlotMinimum(0.7);
  drawPlots(var,nbins,low,high,xtitle,"Events", "HT_3B"+modestring,0,"GeV");
  setLogY(false); resetPlotMinimum();
*/

  float binsHT[]={400,600,800,1400};
  float binsMET[]={250,300,350,500};
  setColorScheme("nostack");
  //  draw2d("MET", 5, 250, 500,"HT", 5, 400, 1400, "MET", "HT", "3B_METvHT");
  //  draw2d("MET", 3, 250, 500,"HT", 3, 400, 1400, "MET", "HT", "3B_METvHT_varbins",binsMET,binsHT);

  draw2d("HT", 5,400, 1400,"MET", 5, 250, 500, "H_{T}", "E_{T}^{miss}", "3B_METvHT");
  draw2d("HT", 3, 400, 1400,"MET", 3, 250, 500, "H_{T}", "E_{T}^{miss}", "3B_METvHT_varbins",binsHT,binsMET);

  
  {for (int ibinx=1 ; ibinx<=h2d->GetNbinsX(); ibinx++) {
    for (int ibiny=1 ; ibiny<=h2d->GetNbinsY(); ibiny++) {
      cout<<"HT="<<h2d->GetXaxis()->GetBinLowEdge(ibinx)<<"-"<<h2d->GetXaxis()->GetBinLowEdge(ibinx)+h2d->GetXaxis()->GetBinWidth(ibinx)
	  <<" MET="<<h2d->GetYaxis()->GetBinLowEdge(ibiny)<<"-"<<h2d->GetYaxis()->GetBinLowEdge(ibiny)+h2d->GetYaxis()->GetBinWidth(ibiny)
	  <<" SM MC="<<h2d->GetBinContent(ibinx,ibiny)<<" +/- "<< h2d->GetBinError(ibinx,ibiny)<<" data="<<hdata2d->GetBinContent(ibinx,ibiny)<<endl;
    }
  }  }

  TH2D* h2d_sigma = (TH2D*)h2d->Clone("h2d_sigma");
  {
  for (int ibinx=1 ; ibinx<=h2d->GetNbinsX(); ibinx++) {
    for (int ibiny=1 ; ibiny<=h2d->GetNbinsY(); ibiny++) {
      double nMC = h2d->GetBinContent(ibinx,ibiny);
      double nData = hdata2d->GetBinContent(ibinx,ibiny);
      double nMC_err = h2d->GetBinError(ibinx,ibiny);
      double nData_err = hdata2d->GetBinError(ibinx,ibiny);

      double nsigma = (nData-nMC) / sqrt( nData_err*nData_err + nMC_err*nMC_err);
      h2d_sigma->SetBinContent(ibinx,ibiny,nsigma);
    }
  }}
  renewCanvas();
  h2d_sigma->Draw("COLZ");
  h2d_sigma->Draw("text same");


  // Now take a look at the 3B SL sample
  btagSFweight_="probge3";
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0))&& MT_Wlep>=0 && MT_Wlep<100 && MET>=200 &&passCleaning==1 && minDeltaPhiN>=4");

  var="MET"; xtitle="E_{T}^{miss}";
  low=200; high = 500; nbins=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_SL_3B_MET"+modestring,0,"GeV");
//   setLogY(true); setPlotMinimum(0.7);
//   drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_SL_3B_MET"+modestring,0,"GeV");
//   setLogY(false); resetPlotMinimum();

  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0))&& MT_Wlep>=0 && MT_Wlep<100 && MET>=250 &&passCleaning==1 && minDeltaPhiN>=4");
  var="HT"; xtitle="H_{T}";
  low=400; high = 1400; nbins=5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "HT_SL_3B"+modestring,0,"GeV");

}

//for the actual AN and paper, use mode 3, logy and doRatio false
void AN2011( TString btagselection="ge1b",const int mode=1, bool logy=false, bool doRatio=false ) {
  doRatioPlot(doRatio);
  setLogY(logy);
  if(logy)  setPlotMinimum(0.5);
  loadSamples();

  //this mode thing is a bit kludgey
  TString modestring="";
  if (mode==1) {
    usePUweight_=false;
    useHTeff_=false;
    useMHTeff_=false;
    thebnnMHTeffMode_ = kOff;
    btagSFweight_="1";
    currentConfig_=configDescriptions_.getDefault(); //completely raw MC
  }
  else if (mode==2 || mode==3) {
    usePUweight_=true;
    useHTeff_=true;
    useMHTeff_=true;    
    thebnnMHTeffMode_ = kOff;
    currentConfig_=configDescriptions_.getCorrected(); //JER bias
    if (mode==2) modestring="-JER-PU-HLT";
    else if(mode==3) modestring="-JER-PU-HLT-bSF";
  }
  else if (mode==4) {
    usePUweight_=true;
    useHTeff_=true;
    useMHTeff_=false;    
    thebnnMHTeffMode_ = kPlot;
    currentConfig_=configDescriptions_.getCorrected(); //JER bias
    modestring="-JER-PU-HLT-BNN-bSF";
  }
  else assert(0);


  TCut btagcut = "nbjetsCSVM>=1";
  if (mode==1 || mode==2) {
    if ( btagselection=="ge1b") {} //do nothing
    else  if ( btagselection=="ge2b" ) {
      btagcut = "nbjetsCSVM>=2";
    }
    else if ( btagselection=="eq1b" ) {
      btagcut = "nbjetsCSVM==1";
    }
    else if ( btagselection=="ge3b" ) {
      btagcut = "nbjetsCSVM>=3";
    }
    else {
      assert(0);
    }
  }
  else if (mode==3 ||mode==4) {
    if ( btagselection=="eq0b") {
      btagcut="1";
      btagSFweight_="prob0";
    }
    else if ( btagselection=="ge1b") {
      btagcut="1";
      btagSFweight_="probge1";
    }
    else  if ( btagselection=="ge2b" ) {
      btagcut = "1";
      btagSFweight_="probge2";
    }
    else if ( btagselection=="ge3b" ) {
      btagcut = "1";
      btagSFweight_="probge3";
    }
    else {
      assert(0);
    }
  }

  loadSamples();
  
  int nbins;
  float low,high;
  TString var,xtitle;
  
  doOverflowAddition(true);
  
  setStackMode(false,true); //normalized
  setColorScheme("nostack");
  clearSamples();
  addSample("PythiaPUQCD");
  addSample("TTbarJets");
  addSample("LM9");

  doData(false);

  /*
  var="minDeltaPhiN"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 40; low=0; high=40;
  //no delta phi cut, loose MET window (only good for MC)
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100")&&btagcut;
  drawPlots(var,nbins,low,high,xtitle,"Arbitrary units", "minDeltaPhiN_looseMET_MConly_"+btagselection+modestring);
  */
  
  //  return;
  
  // ========= regular N-1 plots

  resetSamples(); //use all samples
  setStackMode(true); //regular stack
  setColorScheme("stack");
  doData(true);
  drawMCErrors_=true;
  
  //ben's experience with the Summer11 analysis is that this mod was necessary for style reasons
  //since these variables are global, we can just set them here.
  //they are reset down at the bottom with resetLegendPosition()
  //  leg_x1 = 0.6; leg_x2=0.96; leg_y1=0.4; leg_y2=0.88;//bensep28 - for data/mc stack comparison
  leg_x2=0.96; leg_y1=0.4; //leg_y2=0.88;//bensep28 - for data/mc stack comparison
  //leg_x1 = 0.6; leg_x2=0.96; leg_y1=0.42; leg_y2=0.9;//bensep28 -for data/mc stack comparison NJETS
  
  var="minDeltaPhiN"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 20; low=0; high=40;
  if (btagselection=="ge3b") nbins=10;
  //no delta phi cut
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 &&passCleaning==1")&&btagcut;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_minDeltaPhiN_"+btagselection+modestring);

  // 2BT but with loose MET
  nbins = 10; low=0; high=40;
  selection_ =TCut("HT>=600 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 &&passCleaning==1")&&btagcut;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_minDeltaPhiN_HT600_"+btagselection+modestring);    

  //signal region N-1 minDeltaPhiN
  nbins = 20; low=0; high=40;
  if (btagselection=="ge1b") {
    selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=250 &&passCleaning==1")&&btagcut;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_minDeltaPhiN_1BL_"+btagselection+modestring);

    nbins=10;
    selection_ =TCut("HT>=500 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=500 &&passCleaning==1")&&btagcut;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_minDeltaPhiN_1BT_"+btagselection+modestring);
  }
  else if (btagselection=="ge2b") {
    selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=250 &&passCleaning==1")&&btagcut;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_minDeltaPhiN_2BL_"+btagselection+modestring);

    nbins=10;
    selection_ =TCut("HT>=600 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=300 &&passCleaning==1")&&btagcut;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_minDeltaPhiN_2BT_"+btagselection+modestring);
  }
  else if (btagselection=="ge3b") {
    nbins=10;
    selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=250 &&passCleaning==1")&&btagcut;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_minDeltaPhiN_3B_"+btagselection+modestring);

  }

  // n Jets
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 6; low=2; high=8;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_njets_"+btagselection+modestring);
  
  // n Nets - QCD dominated (LDP)
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN < 4 &&passCleaning==1")&&btagcut;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 7; low=2; high=9;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_njetsLDP_"+btagselection+modestring);
  
  //HT
  selection_ =TCut("cutHT==1&&cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 20; low=400; high=1100;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_HT_"+btagselection+modestring,0,"GeV");
  
  //MET
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 30; low=200; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_"+btagselection+modestring,0,"GeV");

    
  //MET distribution with tighter HT cut
  selection_ =TCut("HT>500 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  TString thisytitle = "";
  if(btagselection=="ge2b") { nbins = 6; low=200; high=500;} 
  else{ nbins = 17; low=200; high=500;}
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_HT500_"+btagselection+modestring,0,"GeV");
  
  selection_ =TCut("HT>600 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  if(btagselection=="ge2b") { nbins = 6; low=200; high=500; } 
  else{ nbins = 17; low=200; high=500; }
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_HT600_"+btagselection+modestring,0,"GeV");

    // == finally, draw the signal region only!
  if (btagselection=="ge1b") {
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=250 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 20; low=400; high=1100;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_HT_1BL_"+btagselection+modestring,0,"GeV");

  var="MET+HT"; xtitle="M_{eff} [GeV]";
  nbins = 20; low=650; high=1550;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_Meff_1BL_"+btagselection+modestring,0,"GeV");

  selection_ =TCut("HT>=500 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=500 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 20; low=500; high=1050;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_HT_1BT_"+btagselection+modestring,0,"GeV");

  var="MET+HT"; xtitle="M_{eff} [GeV]";
  nbins = 5; low=1000; high=2000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_Meff_1BT_"+btagselection+modestring,0,"GeV");
  }
  else if (btagselection=="ge2b") {
    selection_ =TCut("HT>400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=250 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
    var="HT"; xtitle="H_{T} (GeV)";
    nbins = 20; low=400; high=1100;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_HT_2BL_"+btagselection+modestring,0,"GeV");

    var="MET+HT"; xtitle="M_{eff} [GeV]";
    nbins = 20; low=650; high=1550;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_Meff_2BL_"+btagselection+modestring,0,"GeV");

    selection_ =TCut("HT>600 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=300 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
    var="HT"; xtitle="H_{T} (GeV)";
    nbins = 10; low=600; high=1600;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_HT_2BT_"+btagselection+modestring,0,"GeV");

    var="MET+HT"; xtitle="M_{eff} [GeV]";
    nbins = 6; low=900; high=2000;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_Meff_2BT_"+btagselection+modestring,0,"GeV");

  }
  else if (btagselection=="ge3b") {
    selection_ =TCut("HT>400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=250 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
    var="HT"; xtitle="H_{T} (GeV)";
    nbins = 5; low=400; high=1200;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_HT_3B_"+btagselection+modestring,0,"GeV");

    var="MET+HT"; xtitle="M_{eff} [GeV]";
    nbins = 5; low=650; high=1650;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_Meff_3B_"+btagselection+modestring,0,"GeV");

    selection_ =TCut("HT>400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=350 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
    nbins = 5; low=750; high=1750;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_Meff_3B_TighterMET_"+btagselection+modestring,0,"GeV");

  }
  // ===== single lepton selection
  
  // == MET for the electron sample
  selection_ =TCut("MET>=200 && HT>=400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && nElectrons==1 && nMuons==0 && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=200; high=600;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_1e0mu_"+btagselection+modestring,0,"GeV");
   
  //pT?
  var="eleet1"; xtitle="electron p_{T} [GeV]";
  nbins = 20; low=0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_eleET_1e0mu_"+btagselection+modestring,0,"GeV");
  
  
  // == MET for the muon sample
  selection_ =TCut("MET>=200 && HT>=400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && nElectrons==0 && nMuons==1 && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=200; high=600;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_0e1mu_"+btagselection+modestring,0,"GeV");
  
  var="muonpt1"; xtitle="muon p_{T} [GeV]";
  nbins = 20; low=0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_muonpT_0e1mu_"+btagselection+modestring,0,"GeV");
   
  
  // == MET for the combined sample
  selection_ =TCut("MET>=200 &&cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=200; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_SL_"+btagselection+modestring,0,"GeV");
  
  // == HT for the combined sample
  selection_ =TCut("MET>=200 &&cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100")&&btagcut;
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 20; low=350; high=1050;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_HT_SL_"+btagselection+modestring,0,"GeV");

  /*
  // == njets for the combined sample
  selection_ =TCut("MET>=200 && HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 &&  ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 && MET>=200 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 7; low=2; high=9;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_njets_SL_"+btagselection+modestring);
  
  // == njets for the combined sample (HT>500)
  selection_ =TCut("MET>=200 && HT>=500 && cutPV==1 && cutTrigger==1 && cut3Jets==1 &&  ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 && MET>=200 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 7; low=2; high=9;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_njets_SL_HT500_"+btagselection+modestring);

  // == njets for the combined sample (HT>600)
  selection_ =TCut("MET>=200 && HT>=600 && cutPV==1 && cutTrigger==1 && cut3Jets==1 &&  ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 && MET>=200 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 7; low=2; high=9;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_njets_SL_HT600_"+btagselection+modestring);
  */
  
  // == MET for the combined sample (HT>500)
  selection_ =TCut("MET>=200 && HT>500 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=200; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_SL_HT500_"+btagselection+modestring,0,"GeV");
  
  // == MET for the combined sample (HT>600)
  selection_ =TCut("MET>=200 && HT>600 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=200; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_SL_HT600_"+btagselection+modestring,0,"GeV");
  

  //different scale to compare to Kristen
//   selection_ =TCut("MET>=150 && HT>500 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100")&&btagcut;
//   var="MET"; xtitle="E_{T}^{miss} [GeV]";
//   nbins = 35; low=150; high=500;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_SL_HT500_morebins_"+btagselection+modestring);

  
/*
  // == take a look at WJets
  //doesn't work with the fancy b tag SF code. And we're not using it anyway
  selection_ ="MET>=200 && cutHT==1 && cutPV==1 && cutTrigger==1 && njets==2 && (nElectrons==0 && nMuons==1) && minDeltaPhiN >= 4 && nbjetsCSVM==1";
  var="MT_Wlep"; xtitle="M_{T} [GeV]";
   nbins = 20; low=0; high=200;
   drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MT_0e1mu_eq1b");
*/
/*
  //MET
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN < 4")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 30; low=150; high=450;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_ldp_"+btagselection+modestring);
 
  */

  /*
  //SB njet
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && MET<250 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 7; low=2; high=9;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SB_njets_"+btagselection+modestring);
  //SB HT
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && MET<250 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 20; low=350; high=1050;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SB_HT_"+btagselection+modestring,0,"GeV");
  //SB mindphin
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && MET<250 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="minDeltaPhiN"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 20; low=0; high=40;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SB_mindphin_"+btagselection+modestring);
  */


  /*
  //==single-tau selection
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nTaus==1 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  //nbins = 14; low=150; high=500;
  const int nvarbins=3;
  const float varbins[]={150, 200, 250, 500}; //AN and PAS 
  //drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_ST_"+btagselection+modestring,0,"GeV");
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "SBandSIG_MET_ST_"+btagselection+modestring);
  */

  resetLegendPosition();
}


void morerstuff() {
  /*
.L drawReducedTrees.C++
  */
  setStackMode(false);
  doData(false);

  useFlavorHistoryWeights_=false;

  loadSamples(); //needs to come first! this is why this should be turned into a class, with this guy in the ctor
  clearSamples();
  addSample("PythiaPUQCD");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
  setLogY(true);
  drawR("minDeltaPhi",0.3,13,70,230,"ge1b");
  setLogY(false);
  drawR("minDeltaPhi",0.3,13,70,230,"ge1b");

  //  TF1 *fa = new TF1("fa","[0]*exp(-1*[1]*x)",-3,3); 

}


void biasStudy() {
  /*
.L drawReducedTrees.C++
  */
  useFlavorHistoryWeights_=false;

  setStackMode(false);
  doData(true);

  //no met, no mindeltaphi
  setPlotMinimum(1e-2); //setPlotMaximum(1);

 drawTotalSM_=true;
 leg_y1=0.6;
 setLogY(true);
 
 setQuiet(true);

 double cbias,err;

 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 &&njets==3";
 drawR("minDeltaPhi",0.3,4,0,200,"ge1b-4bins");
 cbias= hinteractive->GetBinContent(4)/hinteractive->GetBinContent(3);
 err = jmt::errAoverB(hinteractive->GetBinContent(4),hinteractive->GetBinError(4),hinteractive->GetBinContent(3),hinteractive->GetBinError(3));
 cout<<"no extra cut cbias = "<<cbias <<" +/- "<<err<<endl;

 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 && (MET/HT<0.2) &&njets==3";
 drawR("minDeltaPhi",0.3,4,0,200,"ge1b-METoverHTcut-4bins");
 cbias= hinteractive->GetBinContent(4)/hinteractive->GetBinContent(3);
 err = jmt::errAoverB(hinteractive->GetBinContent(4),hinteractive->GetBinError(4),hinteractive->GetBinContent(3),hinteractive->GetBinError(3));
 cout<<"MET/HT<0.2   cbias = "<<cbias <<" +/- "<<err<<endl;

 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 && (MET/HT>=0.2) &&njets==3";
 drawR("minDeltaPhi",0.3,4,0,200,"ge1b-METoverHTcut-4bins");
 cbias= hinteractive->GetBinContent(4)/hinteractive->GetBinContent(3);
 err = jmt::errAoverB(hinteractive->GetBinContent(4),hinteractive->GetBinError(4),hinteractive->GetBinContent(3),hinteractive->GetBinError(3));
 cout<<"MET/HT>=0.2  cbias = "<<cbias <<" +/- "<<err<<endl;


}

void drawrplots() {
  /*
.L drawReducedTrees.C++
  */

  useFlavorHistoryWeights_=false;

  setStackMode(false);
  doData(true);

  //no met, no mindeltaphi
  setPlotMinimum(1e-2); //setPlotMaximum(1);

 drawTotalSM_=true;
 leg_y1=0.6;
 setLogY(true);

 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,20,0,200,"ge1b");

 selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,20,0,200,"eq1b");

 selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,20,0,200,"ge2b");

 //coarse
 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,4,0,200,"ge1b-4bins");

 selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,4,0,200,"eq1b-4bins");

 selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,4,0,200,"ge2b-4bins");
//  const double epsilon=-1;
//  double cbias= hinteractive->GetBinContent(4)/hinteractive->GetBinContent(3);
//  double cbias_syst = (1-cbias)/2;
//  double cbias_up = cbias+cbias_syst;
//  double cbias_down = cbias-cbias_syst;

//  double r_up = hinteractive->GetBinContent(3) * cbias_up;
//  double r_down = hinteractive->GetBinContent(3) * cbias_down;
//  TLine lup(hinteractive->GetBinCenter(4)+epsilon, r_down, hinteractive->GetBinCenter(4)+epsilon, r_up);
//  lup.SetLineColor(kMagenta);
//  lup.SetLineWidth(2);
//  lup.Draw("SAME");

//not so useful
 selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 && bestTopMass>300";
 drawR("minDeltaPhi",0.3,4,0,200,"ge2b-m3j300");

 //more useful? maybe
 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 && bestTopMass>450";
 drawR("minDeltaPhi",0.3,4,0,200,"ge1b-m3j450");

 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 && bestTopMass>450 && MET>=100&&MET<150";
 drawR("minDeltaPhi",0.3,5,100,150,"ge1b-SBonly-m3j450");


 //njets == 3
// selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && cutCleaning==1 && njets==3";
// setLogY(true);
//  drawR("minDeltaPhi",0.3,50,0,250);

//  //4 jets
// selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && cutCleaning==1 && njets==4";
// setLogY(true);
//  drawR("minDeltaPhi",0.3,50,0,250);

//  //5 jets
// selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && cutCleaning==1 && njets==5";
// setLogY(true);
//  drawR("minDeltaPhi",0.3,50,0,250);

}

void drawSomething() {
  /* for reference, here is all cuts (>=0 b)
selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && cutCleaning==1";

//here is all cuts but with ECAL cleaning removed
selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";

.L drawReducedTrees.C++

  */

  //here we imitate exactly the plots we made in plotWithScale() of the old code

  setStackMode(true);
  doData(true);
  doRatioPlot(true);

  int nbins;
  float low,high;
  TString var,xtitle;

  //  setQCDErrorMode(true);
  // PROBLEMS -- CMS style ruins the hash marks. central value on the TGraphErrors is wrong!
  //might be easiest to accumulate a total MC uncertainty and plot that!

  //the strength and weakness of this tree-based method is that I need to apply the cuts now!
  //most straightforward way is to manually set the cut string before each plot

  // TO DO add weight<1000 cut everywhere!


  // ==== MET plots ====
  setLogY(true);   setPlotMinimum(1e-1);
  nbins=25; low= 0; high=260;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;

  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET");
  //  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_JER0"); //in case we're doing the non-JERbias version

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_ge2b");

  setLogY(true);   setPlotMinimum(1e-1);
  nbins=50; low= 0; high=500;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0.5; ratioMax = 1.5;
  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_METw");

  // == MET plot in limited range (SB and up) ==
  setLogY(false);   setPlotMinimum(0);
  doRatioPlot(false);
  nbins=10; low= 100; high=200;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_SB_ge1b");

  // == MET plot in fail minDeltaPhi
  nbins=10; low= 100; high=200;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_SBfail_ge1b");

  doRatioPlot(true);
  // === MET plot in ILV/T2 region ===
  setLogY(false);   setPlotMinimum(0);
  nbins=30; low= 0; high=300;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1) || (nElectrons==1 && nMuons==0)) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_ILV_ge1b");

  nbins=20; low= 0; high=300;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_IEV_ge1b");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_IMV_ge1b");

  // == replacements for Don's ILV plots (but without fancy ttbar categories)
  //these are in the MET>150 region
  doRatioPlot(false);
  setLogY(false);   setPlotMinimum(0);
  nbins=6; low=1; high=7;
  var="njets"; xtitle="Jet multiplicity";
  //ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_imv_1btag");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_imv_eq1btag");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_imv_2btag");

  doRatioPlot(true);
  setLogY(false);   setPlotMinimum(0);
  nbins=10; low= 0; high=200;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_imv_1btag");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_imv_eq1btag");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_imv_2btag");

  doRatioPlot(false);
  setLogY(false);   setPlotMinimum(0);
  nbins=6; low=1; high=7;
  var="njets"; xtitle="Jet multiplicity";
  //ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_iev_1btag");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_iev_eq1btag");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_iev_2btag");

  doRatioPlot(true);
  setLogY(false);   setPlotMinimum(0);
  nbins=10; low= 0; high=200;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_iev_1btag");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_iev_eq1btag");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_iev_2btag");


  // ==== MHT plot === (no MET cut)
  setLogY(true);  setPlotMinimum(1e-1);
  nbins=50; low= 0; high=500;
  ratioMin = 0.5; ratioMax = 1.5;
  var="MHT"; xtitle="H_{T}^{miss} [GeV]";

  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MHT");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MHT_ge1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MHT_ge2b");

  //to compare to RA2
  setLogY(true);  setPlotMinimum(1);
  nbins=16; low= 0; high=160;
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MHTveryloose");


  // ==== minDeltaPhi plots ====
  setLogY(false);
  resetPlotMinimum();
  nbins=10; low=0; high = TMath::Pi() + 0.001;
  var="minDeltaPhi"; xtitle="min(#Delta#phi[ jets 1..3, E_{T}^{miss} ] )";

  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj");
  //drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_JER0");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_ge2b");

  // ==== HT plots ==== //in this case this is not an N-1 plot, because we've already cut on HT
  ratioMin = 0; ratioMax = 3;
  resetPlotMinimum();
  setLogY(false);
  nbins=10; low=300; high = 1000;
  var="HT"; xtitle="H_{T} (GeV)";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_HT_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_HT_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_HT_ge2b");

  // ==== n jets ====
  nbins=8; low=0; high = 8;
  resetPlotMinimum();
  setLogY(false);
  var="njets"; xtitle="Jet multiplicity";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hnjets_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hnjets_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 &&  MET>=100 && MET<150&& cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hnjets_ge2b");

  // ==== n b jets ====
  ratioMin = 0; ratioMax = 2;
  nbins=4; low=0; high = 4;
  resetPlotMinimum();
  setLogY(false);
  var="nbjets"; xtitle="Number of b-tagged jets";
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && cut3Jets==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hnbjets");

  // ==== n Electrons ====
  nbins=3; low=0; high= 3;
  resetPlotMinimum();
  setLogY(false);
  var="nElectrons"; xtitle="Number of electrons";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnElectrons_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnElectrons_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnElectrons_ge2b");

  // ==== n Muons ====
  nbins=3; low=0; high= 3;
  resetPlotMinimum();
  setLogY(false);
  var="nMuons"; xtitle="Number of muons";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnMuons_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnMuons_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnMuons_ge2b");

  // ==== lead jet pT ====
  nbins=10; low=50; high= 450;
  resetPlotMinimum();
  setLogY(false);
  var="jetpt1"; xtitle="p_{T} of lead jet (GeV)";

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_ge2b");

  // ======= Owen's variables =====
  ratioMin = 0; ratioMax = 3;
  nbins=10; low=50; high= 550;
  resetPlotMinimum();
  setLogY(false);
  var="bestTopMass"; xtitle="3-jet mass (GeV)";

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge2b");

  //signal region
  nbins=20; low=0; high= 800;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge1b_SIG");
  //SB region
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge1b_wide");
  //LSB region
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET<50";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge1b_LSB");

  // ======= study cleaning  =======
  //  bool passBadPFMuon, passInconsistentMuon, passEcalCleaning;
  //no MET cut -- look at BadPFMuon events
  setLogY(true);  resetPlotMinimum();
  nbins=20; low= 0; high=260;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;

  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passBadPFMuon==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "MET_BadPFMuon");

  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "MET_InconsistentMuon");

  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passEcalCleaning==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "MET_FailBE");

  //plot HT (signal region), no b tagging, with and without cleaning
  resetPlotMinimum();
  setLogY(false);
  nbins=10; low=300; high = 1000;
  var="HT"; xtitle="H_{T} (GeV)";
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passEcalCleaning==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HT_FailBE");


   // ==== lead b jet pT ====
   nbins=10; low=30; high= 200;
   resetPlotMinimum();
   setLogY(false);
   var="bjetpt1"; xtitle="p_{T} of lead b jet (GeV)";

   selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
   drawPlots(var,nbins,low,high,xtitle,"Events", "Hbjetpt1_ge1b");

//   selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
//   drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_eq1b");

//   selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
//   drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_ge2b");
   
// ==== number of primary vertices ====
   nbins =9; low=1; high=9;
   resetPlotMinimum();
   setLogY(false);
   var="nGoodPV"; xtitle="# of good PVs";
   selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "nGoodPV_SB_ge1b");
   selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "nGoodPV_SBfail_ge1b");

  //check bad PV events
   nbins =6; low=0; high=5;
   resetPlotMinimum();
   setLogY(false);
   var="nGoodPV"; xtitle="# of good PVs";
   selection_ ="nbjets>=1 && cutHT==1 && cutPV==0 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&MET>50";
   drawPlots(var,nbins,low,high,xtitle,"Events", "nGoodPV_PVfail_ge1b");
   selection_ ="cutHT==1 && cutPV==0 && cutTrigger==1 && cutDeltaPhi==1";
   drawPlots(var,nbins,low,high,xtitle,"Events", "nGoodPV_PVfail_veryloose");


   //exploration
   setLogY(false);
   resetPlotMinimum();
 
   nbins=20; low=-TMath::Pi(); high = TMath::Pi();
   var="METphi-jetphi1"; xtitle="(METphi-jetphi1)";
   selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>100 && passInconsistentMuon==1 && passBadPFMuon==1";
   drawPlots(var,nbins,low,high,xtitle,"Events", "HdeltaPhiMETjet1_nomindp");

   nbins=20; low=-TMath::Pi(); high = TMath::Pi();
   var="METphi-jetphi2"; xtitle="METphi-jetphi2";
   selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>100 && passInconsistentMuon==1 && passBadPFMuon==1 && (abs(METphi-jetphi1)>1)";
   drawPlots(var,nbins,low,high,xtitle,"Events", "HdeltaPhiMETjet2_metphi1cut");

   var="METphi-jetphi2"; xtitle="METphi-jetphi2";
   selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>150 && passInconsistentMuon==1 && passBadPFMuon==1 && (abs(METphi-jetphi1)>1)";
   drawPlots(var,nbins,low,high,xtitle,"Events", "HdeltaPhiMETjet2_metphi1cut_sig");

   var="METphi-jetphi3"; xtitle="METphi-jetphi3";
   selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>150 && passInconsistentMuon==1 && passBadPFMuon==1 && (abs(METphi-jetphi1)>1) && (abs(METphi-jetphi2)>0.3)";
   drawPlots(var,nbins,low,high,xtitle,"Events", "HdeltaPhiMETjet3_metphi12cut_sig");

  // ==== MET/HT ====
   //addSample("LM13");
  doRatioPlot(false);
  setLogY(true);
  nbins=40;
  low=0;
  high=1;
  doOverflowAddition(true);
  var="MET/HT"; xtitle="MET/HT";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverHt_ge1b_noMetCut");

  setLogY(false);
  resetPlotMinimum();
  nbins=10;

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 && MET>100 &&MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverHt_ge1b_SB");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&MET>=150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverHt_ge1b_SIG");

  // ==== MET/sqrt(HT) ====
   //addSample("LM13");
  doRatioPlot(false);
  nbins=20;  low=0;  high=20;
  doOverflowAddition(true);
  var="MET/sqrt(HT)"; xtitle="MET/#sqrt{HT}";

   setLogY(false);
  resetPlotMinimum();

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 && MET>100 &&MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverRootHt_ge1b_SB");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&MET>=150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverRootHt_ge1b_SIG");


  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&MET>=100";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverRootHt_ge1b_metOver100");

}

void studySignificance() {
/*
.L drawReducedTrees.C++
  */

  setStackMode(true);
  doData(true);

  int nbins;
  float low,high;
  TString var,xtitle;
  doOverflowAddition(true);
  doRatioPlot(false);


  setLogY(false);
  resetPlotMinimum();

  //this compares the effectiveness of MET and MET/sqrt(HT); MET is better
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";

  nbins=20;  low=0;  high=20;
  var="MET/sqrt(HT)"; xtitle="MET/#sqrt{HT}";
  drawSignificance(var,nbins,low,high,"sigfile");

  nbins=30;  low=0;  high=300;
  var="MET"; xtitle="MET";
  drawSignificance(var,nbins,low,high,"sigfile");


  //after MET cut, what can we do?
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=20;  low=0;  high=20;
  var="MET/sqrt(HT)"; xtitle="MET/#sqrt{HT}";
  drawSignificance(var,nbins,low,high,"sigfile"); //nothing good here

  //reoptimize minDeltaPhi
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=31;  low=0;  high=TMath::Pi();
  var="minDeltaPhi"; xtitle="minDeltaPhi";
  drawSignificance(var,nbins,low,high,"sigfile");

  var="minDeltaPhiAll30"; xtitle="minDeltaPhi All 30"; //not as good!
  drawSignificance(var,nbins,low,high,"sigfile");


  //remove 3Jet cut
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=7; low=0; high=7;
  var="njets"; xtitle="jet multiplicity";
  drawSignificance(var,nbins,low,high,"sigfile");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=40; low=200; high=1000;
  var="HT"; xtitle="HT";
  drawSignificance(var,nbins,low,high,"sigfile");

  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=4; low=0; high=4;
  var="nbjets"; xtitle="n b tags";
  drawSignificance(var,nbins,low,high,"sigfile");

  //nominal cuts
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=40;  low=0;  high=400;
  var="jetpt1"; xtitle="lead jet pt";
  drawSignificance(var,nbins,low,high,"sigfile");

  nbins=40;  low=0;  high=400;
  var="bjetpt1"; xtitle="lead b jet pt";
  drawSignificance(var,nbins,low,high,"sigfile");

}

void drawMinDeltaPhiMETslices(){
  useFlavorHistoryWeights_=false;
  doOverflowAddition(true);
  doData(true);
  doRatioPlot(true);
  setQuiet(false);
  ratioMin=0; ratioMax=3;
  
  bool doVB = false;
  bool doSlices = true;
  if(!doSlices){
    leg_y1=0.6; //and also remember to use correct color scheme
    owenColor_=true;
  }  

  int nbins;
  float low,high;
  TString var,xtitle;

  nbins=10; low=0; high = TMath::Pi() + 0.001;
  var="minDeltaPhi"; xtitle="min(#Delta#phi[ jets 1..3, E_{T}^{miss} ] )";
  const int nvarbins=2;
  const float varbins[]={0.,0.3,high};
  
  TCut baseSelection =""; 
  TCut theseCuts = ""; 
  TString extraName = "";

  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjets == 0";
  for (int ibtag = 0; ibtag<5; ibtag++) { 
    TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
    if (ibtag==0) { //nothing to do
    }
    else if (ibtag==1) {
      theBTaggingCut = eq1b; 
      btagstring = "eq1b";
    }
    else if (ibtag==2) {
      theBTaggingCut = ge2b; 
      btagstring = "ge2b";
    }
    else if (ibtag==3) {
      theBTaggingCut = pretag;
      btagstring = "pre";
    }
    else if (ibtag==4) {
      theBTaggingCut = antitag;
      btagstring = "antib";
    }
    else assert(0);
    

    //+++++++++++NOMINAL SELECTION+++++++++++//
    baseSelection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1";//nominal
    extraName = "_nominal";
     
    if(doSlices){
      //MinDeltaPhi in MET slices
      setStackMode(true); setLogY(false); resetPlotMinimum();
      //MET < 50
      theseCuts = "MET<50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50_vb");
      //50 < MET < 100 
      theseCuts = "MET<100. && MET>=50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100_vb");
      //100 < MET < 150
      theseCuts = "MET<150. && MET>=100.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150_vb");
      //MET < 150
      theseCuts = "MET>=150.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf"); 
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf_vb");
    }
    else{
      //Ratio plot
      selection_ = baseSelection && theBTaggingCut;
      setStackMode(false); setPlotMinimum(1e-2); drawTotalSM_=true; setLogY(true); 
      drawR("minDeltaPhi",0.3,20,0,200,btagstring+extraName);
      drawR("minDeltaPhi",0.3,4,0,200,btagstring+extraName+"_4bins");
    }
    
    /*
    //[+++++++++++++njets==2++++++++]//
    baseSelection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1 && njets==2"; //remove cut3Jets and add njets
    extraName = "_njetsEQ2";

    if(doSlices){
      //MinDeltaPhi in MET slices
      setStackMode(true); setLogY(false); resetPlotMinimum();
      //MET < 50
      theseCuts = "MET<50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50_vb");
      //50 < MET < 100 
      theseCuts = "MET<100. && MET>=50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100_vb");
      //100 < MET < 150
      theseCuts = "MET<150. && MET>=100.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150_vb");
      //MET < 150
      theseCuts = "MET>=150.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf_vb");
    }
    else{
      //Ratio plot
      selection_ = baseSelection && theBTaggingCut;
      setStackMode(false); setPlotMinimum(1e-2); drawTotalSM_=true; setLogY(true); 
      drawR("minDeltaPhi",0.3,20,0,200,btagstring+extraName);
      drawR("minDeltaPhi",0.3,4,0,200,btagstring+extraName+"_4bins");
    }
    */
    
    /*
    //[+++++++++++m3j>350++++++++++++]//
    baseSelection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1 && bestTopMass>350";//nominal plus m3j cut
    extraName = "_m3jGT350";
    
    if(doSlices){
      //MinDeltaPhi in MET slices
      setStackMode(true); setLogY(false); resetPlotMinimum();
      //MET < 50
      theseCuts = "MET<50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50_vb");
      //50 < MET < 100 
      theseCuts = "MET<100. && MET>=50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100_vb");
      //100 < MET < 150
      theseCuts = "MET<150. && MET>=100.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150_vb");
      //MET < 150
      theseCuts = "MET>=150.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf_vb");
    }
    else{
      //Ratio plot
      selection_ = baseSelection && theBTaggingCut;
      setStackMode(false); setPlotMinimum(1e-2); drawTotalSM_=true; setLogY(true); 
      drawR("minDeltaPhi",0.3,20,0,200,btagstring+extraName);
      drawR("minDeltaPhi",0.3,4,0,200,btagstring+extraName+"_4bins");
    }
    */

  }//b-tag cut loop 
}


void drawMETPlots() {
  useFlavorHistoryWeights_=false;
  setStackMode(true);
  doData(true);
  doRatioPlot(true);
  setQuiet(false);
  doOverflowAddition(true);
  loadSamples();

  int nbins;
  float low,high;
  TString var,xtitle;

  // ==== MET plots ====
  setLogY(true);   setPlotMinimum(1e-1);
  nbins=10; low= 0.; high=150.;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 3;

  /*
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_ge2b");

  selection_ ="nbjets==0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_antib");
  */ 


  nbins=20; low=0; high=400;
  //double tempLumiScale =  lumiScale_;
  //double thisLumi = 94.432+166.732+216.430+154.830+84.698+57.994+2.101;
  //cout << "thisLumi: " << thisLumi << endl;
  //lumiScale_ = thisLumi/190.5;
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_test");
  //lumiScale_= tempLumiScale;
}


void drawMETPlots_utility(){
  useFlavorHistoryWeights_=false;
  setStackMode(true);
  doData(true);
  doRatioPlot(true);
  setQuiet(false);
  doOverflowAddition(true);
  loadSamples();
  ratioMin = 0; ratioMax = 3;

  int nbins;
  float low,high;
  TString var,xtitle;
  bool doVB=false;

  TCut baseSelection = "cutHT==1 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1"; //cuts not present: trigger, MET, mdp, bjets
  TCut utilSelection = "pass_utilityHLT_HT300>=1";
  TCut nominalTrigger = "cutTrigger==1";
  TCut mdpCut = "cutDeltaPhi==1";
  TCut invertMdpCut = "cutDeltaPhi==0";

  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjets == 0";
  for (int ibtag = 0; ibtag<5; ibtag++) { 
    TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
    if (ibtag==0) { //nothing to do
    }
    else if (ibtag==1) {
      theBTaggingCut = eq1b; 
      btagstring = "eq1b";
    }
    else if (ibtag==2) {
      theBTaggingCut = ge2b; 
      btagstring = "ge2b";
    }
    else if (ibtag==3) {
      theBTaggingCut = pretag;
      btagstring = "pre";
    }
    else if (ibtag==4) {
      theBTaggingCut = antitag;
      btagstring = "antib";
    }
    else assert(0);
    
    TString extraName = "_nominal";
    
    // ==== MET plots ====
    setLogY(true);   setPlotMinimum(1e-1);
    nbins=10; low= 0.; high=150.;
    var="MET"; xtitle="E_{T}^{miss} [GeV]";
    selection_ = baseSelection && utilSelection && mdpCut && theBTaggingCut; //no MET cut
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_METstack");
    setLogY(false); resetPlotMinimum();
    //high MET
    nbins=1; low=150.; high=200.;
    selection_ = baseSelection && utilSelection && invertMdpCut && theBTaggingCut; //no MET cut, inverted mdp
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_imdp_METstack");
    //cross check hack
    double tempLumiScale =  lumiScale_;
    double thisLumi = 94.432+166.732+216.430+154.830+84.698+57.994+2.101;
    cout << "thisLumi: " << thisLumi << endl;
    lumiScale_ = thisLumi/190.5;
    selection_ = baseSelection && nominalTrigger && invertMdpCut && theBTaggingCut; //no MET cut, inverted mdp
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_imdp_METstack_full");
    lumiScale_= tempLumiScale;
    

    if(0){
      // ==== MinDeltaPhi in MET slices ====
      setLogY(false); resetPlotMinimum();
      nbins=7;  low = 0; high = TMath::Pi() + 0.001;
      var="minDeltaPhi"; xtitle="min(#Delta#phi[ jets 1..3, E_{T}^{miss} ] )";
      const int nvarbins=2;
      const float varbins[]={0.,0.3,high};
      TCut theseCuts;
      //MET < 50
      theseCuts = "MET<50.";
      selection_ = baseSelection && utilSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50_vb");
      //50 < MET < 100 
      theseCuts = "MET<100. && MET>=50.";
      selection_ = baseSelection && utilSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100_vb");
      //100 < MET < 150
      theseCuts = "MET<150. && MET>=100.";
      selection_ = baseSelection && utilSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150_vb");
      //MET < 150
      theseCuts = "MET>=150.";
      selection_ = baseSelection && utilSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf"); 
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf_vb");
    }
  }//end loop over b-tag selection
}

void qcdFit(bool loose = true){
  loadSamples(true);
  savePlots_=false; //don't save eps,png,pdf files
  doOverflowAddition(false);

  //output file 
  TString histfilename = "fitHists.root";
  TFile fh(histfilename,"RECREATE");//will delete old root file if present 
  fh.Close(); //going to open only when needed 
  
  //-- define some cuts to use 
  //--------------------------
  const TCut phys = "cutTrigger==1";
  const TCut util = "pass_utilityHLT_HT300>=1";
  const TCut LSB = "MET>=50 && MET<100";
  const TCut SB = "MET>=150 && MET<200";
  TCut SIG;
  if(loose) SIG = "MET>=200";
  else SIG =  "MET>=300";
  TCut HTcut;
  if(loose) HTcut = "HT>=350";
  else HTcut = "HT>=500";
  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjets == 0";
  const TCut passmdpn = "minDeltaPhiN>=4";
  const TCut failmdpn = "minDeltaPhiN<4";
  TCut base = HTcut&&TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && weight<1000");// no trigger, met, minDeltaPhiN, btagger

  TString var = "MET";
  int nbins; double low; double high;
  
  selection_ = base&&util&&LSB&&antitag&& passmdpn;
  nbins=1; low=50; high=100;
  drawSimple(var,nbins,low,high,histfilename, "h_LSB_fDP_antib_MET_data", "data");
  drawSimple(var,nbins,low,high,histfilename, "h_LSB_fDP_antib_MET_qcd", "PythiaPUQCD");
  drawPlots( var,nbins,low,high, "", "",      "deleteme");  
  TFile fh1(histfilename, "UPDATE");
  totalnonqcd->SetName("h_LSB_fDP_antib_MET_nonqcd");
  totalnonqcd->Write();
  fh1.Close();

}

void drawQCDreweight(bool loose = true, bool njets50 = false){
  loadSamples(true);
  savePlots_=false; //don't save eps,png,pdf files
  doOverflowAddition(true);

  //for Summer11 result, run drawQCDreweight(true,false) and drawQCDreweight(false,false) with b12same=true.
  bool b12same=true;

  TString var = "";
  if(njets50) var = "njets";
  else var = "njets30";

  //output file 
  TString histfilename = "";
  histfilename += "qcdReweight";
  if(loose) histfilename += "_loose_";
  else histfilename += "_tight_";
  if(njets50) histfilename += "50";
  else histfilename += "30";
  histfilename += ".root";
  TFile fh(histfilename,"RECREATE");//will delete old root file if present 
  fh.Close(); //going to open only when needed 
  
  //-- define some cuts to use 
  //--------------------------
  const TCut phys = "cutTrigger==1";
  const TCut util = "pass_utilityHLT_HT300>=1";
  const TCut LSB = "MET>=50 && MET<100";
  const TCut SB = "MET>=150 && MET<200";
  TCut SIG;
  if(loose) SIG = "MET>=200";
  else SIG =  "MET>=300";
  TCut HTcut;
  if(loose) HTcut = "HT>=350";
  else HTcut = "HT>=500";
  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjets == 0";
  const TCut passmdpn = "minDeltaPhiN>=4";
  const TCut failmdpn = "minDeltaPhiN<4";
  TCut base = HTcut&&TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && weight<1000");// no trigger, met, minDeltaPhiN, btagger
  
  //binning
  int nbins;
  double low;
  double high;
  
  //-- LSB plots
  //------------
  lumiScale_ = 19.23288;
  if(njets50) {nbins=4; low=2.5; high=6.5;}
  else{        nbins=5; low=2.5; high=7.5;}
  selection_ = base && LSB && util && failmdpn && antitag;
  drawSimple(var,nbins,low,high,histfilename, "h_LSB_fDP_antib_njets_data", "data");
  drawSimple(var,nbins,low,high,histfilename, "h_LSB_fDP_antib_njets_qcd", "PythiaPUQCD");
  drawPlots( var,nbins,low,high, "", "",      "deleteme");  
  TFile fh1(histfilename, "UPDATE");
  totalnonqcd->SetName("h_LSB_fDP_antib_njets_nonqcd");
  totalnonqcd->Write();
  fh1.Close();
  selection_ = base && LSB && util && passmdpn && antitag;
  drawSimple(var,nbins,low,high,histfilename, "h_LSB_pDP_antib_njets_qcd", "PythiaPUQCD");
  
  //-- SB plots
  //-----------
  lumiScale_ = 1143.;
  
  if(njets50) { nbins=4; low=2.5; high=6.5;} 
  else{         nbins=4; low=2.5; high=6.5;}
  selection_ =  base && SB && phys && failmdpn && ge2b;
  drawSimple(var,nbins,low,high,histfilename, "h_SB_fDP_ge2b_njets_data", "data");
  drawSimple(var,nbins,low,high,histfilename, "h_SB_fDP_ge2b_njets_qcd", "PythiaPUQCD");
  drawPlots( var,nbins,low,high, "", "",      "deleteme");
  TFile fh3(histfilename, "UPDATE");
  totalnonqcd->SetName("h_SB_fDP_ge2b_njets_nonqcd");
  totalnonqcd->Write();
  fh3.Close();
  selection_ = base && SB && phys && passmdpn && ge2b;
  drawSimple(var,nbins,low,high,histfilename, "h_SB_pDP_ge2b_njets_qcd", "PythiaPUQCD");
  
  if(!b12same){
    if(njets50) { nbins=4; low=2.5; high=6.5;} 
    else{         nbins=5; low=2.5; high=7.5;}
  }
  selection_ =  base && SB && phys && failmdpn && ge1b;
  drawSimple(var,nbins,low,high,histfilename, "h_SB_fDP_ge1b_njets_data", "data");
  drawSimple(var,nbins,low,high,histfilename, "h_SB_fDP_ge1b_njets_qcd", "PythiaPUQCD");
  drawPlots( var,nbins,low,high, "", "",      "deleteme");   
  TFile fh2(histfilename, "UPDATE");
  totalnonqcd->SetName("h_SB_fDP_ge1b_njets_nonqcd");
  totalnonqcd->Write();
  fh2.Close(); 
  selection_ = base && SB && phys && passmdpn && ge1b;
  drawSimple(var,nbins,low,high,histfilename, "h_SB_pDP_ge1b_njets_qcd", "PythiaPUQCD");

  //-- SIG plots
  //------------
  if(njets50){
    if(loose){ nbins=3; low=2.5; high=5.5; }
    else{      nbins=2; low=2.5; high=4.5; }
  }
  else{
    if(loose){ nbins=3; low=2.5; high=5.5; }
    else{      nbins=2; low=2.5; high=4.5; }
  }
  selection_ =  base && SIG && phys && failmdpn && ge2b;
  drawSimple(var,nbins,low,high,histfilename, "h_SIG_fDP_ge2b_njets_data", "data");
  drawSimple(var,nbins,low,high,histfilename, "h_SIG_fDP_ge2b_njets_qcd", "PythiaPUQCD");
  drawPlots(var,nbins,low,high, "", "",      "deleteme");
  TFile fh5(histfilename, "UPDATE");
  totalnonqcd->SetName("h_SIG_fDP_ge2b_njets_nonqcd");
  totalnonqcd->Write();
  fh5.Close();
  selection_ = base && SIG && phys && passmdpn && ge2b;
  drawSimple(var,nbins,low,high,histfilename, "h_SIG_pDP_ge2b_njets_qcd", "PythiaPUQCD");

  if(!b12same){
    if(njets50){
      if(loose){ nbins=4; low=2.5; high=6.5; }
      else{      nbins=2; low=2.5; high=4.5; }
    }
    else{
      if(loose){ nbins=5; low=2.5; high=7.5; }
      else{      nbins=4; low=2.5; high=6.5; }
    }
  }
  selection_ =  base && SIG && phys && failmdpn && ge1b;
  drawSimple(var,nbins,low,high,histfilename, "h_SIG_fDP_ge1b_njets_data", "data");
  drawSimple(var,nbins,low,high,histfilename, "h_SIG_fDP_ge1b_njets_qcd", "PythiaPUQCD");
  drawPlots(var,nbins,low,high, "", "",      "deleteme");
  TFile fh4(histfilename, "UPDATE");
  totalnonqcd->SetName("h_SIG_fDP_ge1b_njets_nonqcd");
  totalnonqcd->Write();
  fh4.Close();
  selection_ = base && SIG && phys && passmdpn && ge1b;
  drawSimple(var,nbins,low,high,histfilename, "h_SIG_pDP_ge1b_njets_qcd", "PythiaPUQCD");
  
}



void ptresABCD(bool loose = true){
  doOverflowAddition(true);
  loadSamples();
  clearSamples();
  addSample("PythiaPUQCD");
  
  const TCut LSB = "MET>=50 && MET<100";
  const TCut SB = "MET>=150 && MET<200";
  TCut SIG;
  if(loose) SIG = "MET>=200";
  else SIG =  "MET>=300";
  TCut HTcut;
  if(loose) HTcut = "HT>=350";
  else HTcut = "HT>=500";
  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjets == 0";
  const TCut passmdpn = "minDeltaPhiN>=4";
  const TCut failmdpn = "minDeltaPhiN<4";//.1412212
  //const TCut passmdpn = "minDeltaPhiN15>=2.7";
  //const TCut failmdpn = "minDeltaPhiN15<2.7";//.142238
  //const TCut passmdpn = "minDeltaPhiN5>=7.9";
  //const TCut failmdpn = "minDeltaPhiN5<7.9";
  TCut base = HTcut&&TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && weight<1000");// no trigger, met, minDeltaPhiN, btagger
  
  double A, B, C, D, Aerr, Berr, Cerr, Derr;

  TString sampleOfInterest = ("PythiaPUQCD");
  TString var="HT"; TString xtitle=var;
  int nbins=1; double low=0; double high=50000; //ht overflow is on
  

  //A
  selection_ = base&&LSB&&failmdpn&&antitag;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_A");
  A=getIntegral(sampleOfInterest);
  Aerr=getIntegralErr(sampleOfInterest);

  //B
  selection_ = base&&LSB&&passmdpn&&antitag;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_A");
  B=getIntegral(sampleOfInterest);
  Berr=getIntegralErr(sampleOfInterest);

  cout << "**************" << endl;
  cout << "B/A: " << B/A << endl;
  cout << "**************" << endl;
 

  //C
  selection_ = base&&SIG&&passmdpn&&ge1b;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_A");
  C=getIntegral(sampleOfInterest);
  Cerr=getIntegralErr(sampleOfInterest);

  //D
  selection_ = base&&SIG&&failmdpn&&ge1b;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_A");
  D=getIntegral(sampleOfInterest);
  Derr=getIntegralErr(sampleOfInterest);

  double numerr=jmt::errAtimesB(B,Berr,D,Derr);
  double num = B*D;
  double estimate = num / A;
  double estimateerr= jmt::errAoverB(num,numerr,A,Aerr);

  cout << "prediction: " << estimate << " +- " << estimateerr << endl;
  cout << "MC truth  : " << C << " +- " << Cerr << endl;
  
  double num1 = C-estimate;
  double perc = num1/C;
  double percerr = jmt::errAoverB(estimate, estimateerr,C,Cerr);
  cout << "percent   : " << perc*100. << " +- " << percerr*100. << endl;
  cout << "combined  : " << sqrt(perc*perc + percerr*percerr)*100. << endl;

}


void studyPrescale_r(int ibtag = 4) {
  loadSamples();
  drawTotalSM_=true;
  drawLegend(true);
  setColorScheme("nostack");
  setStackMode(false);
  doOverflowAddition(true);
  doRatioPlot(true);
  doData(true);

  lumiScale_ = 30.15471; 

  //TString btagselection="antib";
  //TCut btagcut = "nbjets==0";
  
  TCut util = "pass_utilityHLT==1 && weight<1000";
  TCut LSB = "MET>=50 && MET<100";
  TCut runN = "runNumber>=175832";

  const  TCut ge1b =  "nbjetsCSVM >= 1";
  const  TCut ge2b =  "nbjetsCSVM >= 2";
  const  TCut ge3b =  "nbjetsCSVM >= 3";
  const  TCut eq1b =  "nbjetsCSVM == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjetsCSVM == 0";
 
  //  for (int ibtag = 4; ibtag<5; ibtag++) { 
  TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
  if (ibtag==0) { //nothing to do
  }
  else if (ibtag==1) {
    theBTaggingCut = eq1b; 
    btagstring = "eq1b";
  }
  else if (ibtag==2) {
    theBTaggingCut = ge2b; 
    btagstring = "ge2b";
  }
  else if (ibtag==3) {
    theBTaggingCut = pretag;
    btagstring = "pre";
  }
  else if (ibtag==4) {
    theBTaggingCut = antitag;
    btagstring = "antib";
  }
  else if (ibtag==5) {
    theBTaggingCut = ge3b;
    btagstring = "ge3b";
  }

  else assert(0);
  
  doOverflowAddition(false);
  //cross check with note
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN",4,10,50,150,"Data_eq0b");
  //drawR("minDeltaPhiN",4,1,50,100,"Data_eq0b");
  
  //trigger version
  //set dataOnly and might want to change hdata->SetXTitle()
  doOverflowAddition(false);
  //setPlotMaximum(0.5); setPlotMinimum(0);
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "pass_utilityHLT_HT300", 7, 1.5, 8.5, "triggerVersion_"+btagstring);
  //unexpected value for MC in last bin when adding overflow.
  
  //nGoodPV
  doOverflowAddition(true);
  //setPlotMaximum(0.5); setPlotMinimum(0);
  selection_ =TCut("HT>=400 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "nGoodPV", 20, 0.5, 20.5, "nGoodPV_"+btagstring);
  //const int nvarbins=9;
  //const float varbins[]={0.5,1.5,2.5,3.5,4.5,5.5,6.5,8.5,11.5,16.5};
  //drawR("minDeltaPhiN", 4, "nGoodPV", nvarbins, varbins, "nGoodPV_"+btagstring);
  //dependence is seen here
  //
  
  
  //prescale vs physics check
  selection_ =TCut("HT>=400 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut&&runN;
  drawSimple("nGoodPV",20,0.5,20.5,"nGoodPV.root", "nGoodPV_"+btagstring+"tag_prescale_data","data");
  selection_ =TCut("HT>=400 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=200")&&theBTaggingCut&&runN;
  drawSimple("nGoodPV",20,0.5,20.5,"nGoodPV.root", "nGoodPV_"+btagstring+"tag_physics_data","data");
  


  //njets
  doOverflowAddition(true);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "njets", 8, 2.5, 10.5, "njets_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=150 && weight<1000")&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "njets",8, 2.5, 10.5, "njets_physics_"+btagstring);
  
  //nbjets -- does not depend on theBTaggingCut!!!
  doOverflowAddition(false);
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB;
  //drawR("minDeltaPhiN", 4, "nbjets", 4, -0.5, 3.5, "nbjets_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=150");
  //drawR("minDeltaPhiN", 4, "nbjets", 4, -0.5, 3.5, "nbjets_physics_"+btagstring);
   
  //runnumber
  //set dataOnly to true
  //should hack after first instance of thecanvas->cd(1); to add:  gPad->SetRightMargin(.1); gPad->Modified();
  doOverflowAddition(false);
  setPlotMaximum(0.5); setPlotMinimum(0);
  selection_ =TCut("HT>=400 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "runNumber", 30, 160403.5, 180000.5, "runNumber_"+btagstring);
  //const int nvarbins2 = 7;
  //const float varbins2[] = {160403.5, 161000, 162000, 163500, 165000,166000,167000,168000};
  //drawR("minDeltaPhiN", 4, "runNumber", nvarbins2, varbins2, "runNumber_"+btagstring);
  
  
  //lowMET
  doOverflowAddition(false);
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&theBTaggingCut;
  const int nvarbins=10;
  const float varbins[]={0,10,20,30,40,50,60,70,80,90,100};
  //  drawR("minDeltaPhiN", 4, "MET", nvarbins, varbins, "MET_"+btagstring);
  //}//end btag loop

  //MET fit?
  doOverflowAddition(false);
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "MET", 1, 50, 100, "MET_"+btagstring);

}
 

void lowMETcount(){
  loadSamples();
  doOverflowAddition(false);
  setQuiet(false);
  setStackMode(true); //regular stack
  setColorScheme("stack");
  doData(true);
  drawMCErrors_=true;
  
  lumiScale_ = 19.23288; //customize lumiScale_
  
  TString btagselection="antib";
  TCut theBTaggingCut = "nbjets==0";
  TCut util = "pass_utilityHLT_HT300>=1 && weight<1000";
  
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&theBTaggingCut;
  const int nvarbins=3;
  const float varbins[]={40,50,100,110};
  drawPlots("MET", nvarbins, varbins, "","","deleteme");
  cout << "BEN: Data integral1:::::::: " << hdata->GetBinContent(1) << endl;
  cout << "BEN: Data integral2:::::::: " << hdata->GetBinContent(2) << endl;
  cout << "BEN: Data integral3:::::::: " << hdata->GetBinContent(3) << endl;
}//end of function


void sigbpt(){
  loadSamples();
  doOverflowAddition(true);
  
  TString btagstring="ge1b";
  TCut theBTaggingCut = "nbjets>=1";

  //loose
  selection_ =TCut("(HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN>=4 && MET>=200)")&&theBTaggingCut;
  cout << "LM9: " << drawSimple("bjetpt1",20,30,500,"bjetpt1.root", "bjetpt1_"+btagstring+"tagLooseSignalRegion_lm9", "LM9") << endl;;
}

void studyNjet(int ibtag=3){
  loadSamples();
  drawTotalSM_=true;
  drawLegend(true);
  //setColorScheme("nostack");
  //setStackMode(false);
  doOverflowAddition(true);
  doRatioPlot(true);
  doData(true);

  //lumiScale_ = 19.23288; //customize lumiScale_
  //TString btagselection="antib"; 
  //TCut btagcut = "nbjets==0"; 

  TCut util = "pass_utilityHLT_HT300>=1 && weight<1000";
  TCut LSB = "MET>=50 && MET<100";

  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjets == 0";

  //  for (int ibtag = 4; ibtag<5; ibtag++) {
  TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
  if (ibtag==0) { //nothing to do 
  }
  else if (ibtag==1) {
    theBTaggingCut = eq1b;
    btagstring = "eq1b";
  }
  else if (ibtag==2) {
    theBTaggingCut = ge2b;
    btagstring = "ge2b";
  }
  else if (ibtag==3) {
    theBTaggingCut = pretag;
    btagstring = "pre";
  }
  else if (ibtag==4) {
    theBTaggingCut = antitag;
    btagstring = "antib";
  }
  else assert(0);


  doOverflowAddition(true);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "njets", 7, 2.5, 10.5, "njets_"+btagstring); 
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=150")&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "njets", 7, 2.5, 10.5, "njets_physics_"+btagstring); 

  //njet distributions -- note - NO minDeltahiN cut!
  //selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=150 && weight<1000")&&theBTaggingCut;
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=150 && weight<1000")&&ge2b;
  //drawSimple("njets",10,0,10,"njets.root", "njets_"+btagstring+"tag_sig_qcd","PythiaPUQCD");
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=50 && MET<100 && weight<1000")&&theBTaggingCut;
  //drawSimple("njets",10,0,10,"njets.root", "njets_"+btagstring+"tag_lsb_qcd","PythiaPUQCD");


  //N JET DATA/MC
  setStackMode(true); setLogY(true); resetPlotMinimum();    
  TCut failmdp = "minDeltaPhiN<4.";

  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=150 && weight<1000")&&theBTaggingCut && failmdp;
  drawPlots("njets",7,2.5,10,"njets","Events", "H_njets_physics_"+btagstring);
  
  lumiScale_ = 19.23288; //customize lumiScale_
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut&&failmdp;
  drawPlots("njets",7,2.5,10,"njets","Events", "H_njets_prescaled_"+btagstring);
}

void studyRqcd(int ibtag = 4){
  loadSamples();
  drawTotalSM_=false;
  drawLegend(true);
  setColorScheme("nostack");
  setStackMode(false);
  doOverflowAddition(true);
  doRatioPlot(true);
  doData(false);

  lumiScale_ = 19.23288; //customize lumiScale_

  //TString btagselection="antib";
  //TCut btagcut = "nbjets==0";
  
  TCut util = "pass_utilityHLT_HT300>=1 && weight<1000";
  TCut LSB = "MET>=50 && MET<100";
  TCut SB = "MET>=150 && MET<200";
  TCut SIG = "MET>=200";
  TCut highMET = "MET>=150";

  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjets == 0";
 
  //  for (int ibtag = 4; ibtag<5; ibtag++) { 
  TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
  if (ibtag==0) { //nothing to do
  }
  else if (ibtag==1) {
    theBTaggingCut = eq1b; 
    btagstring = "eq1b";
  }
  else if (ibtag==2) {
    theBTaggingCut = ge2b; 
    btagstring = "ge2b";
  }
  else if (ibtag==3) {
    theBTaggingCut = pretag;
    btagstring = "pre";
  }
  else if (ibtag==4) {
    theBTaggingCut = antitag;
    btagstring = "antib";
  }
  else assert(0);
  
  //clearSamples();
  //addSample("PythiaPUQCD");
  
  /*
    // NJETS
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  drawR("minDeltaPhiN", 4, "njets", 8, 2.5, 10.5, "njets_lsb_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&SB&&theBTaggingCut;
  drawR("minDeltaPhiN", 4, "njets", 8, 2.5, 10.5, "njets_sb_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&SIG&&theBTaggingCut;
  drawR("minDeltaPhiN", 4, "njets", 8, 2.5, 10.5, "njets_sig_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&highMET&&theBTaggingCut;
  drawR("minDeltaPhiN", 4, "njets", 8, 2.5, 10.5, "njets_highMET_"+btagstring);
  */

  /*
  //NJETS FOR JOSH
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  drawR("minDeltaPhiN", 4, "njets", 4, 2.5, 6.5, "njets_joshBin_lsb_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&highMET&&theBTaggingCut;
  drawR("minDeltaPhiN", 4, "njets", 4, 2.5, 6.5, "njets_joshBin_highMET_"+btagstring);
  */  

  //setPlotMaximum(0.35); setPlotMinimum(0);
  //selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cut3Jets==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "MET", 1, 50, 100, "test_"+btagstring);
  
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cut3Jets==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "MET", 15, 50, 350, "test_"+btagstring);

  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && njets==3")&&util&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "MET", 6, 50, 350, "MET_njets3_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && njets==4")&&util&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "MET", 6, 50, 350, "MET_njets4_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && njets>=5")&&util&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "MET", 6, 50, 350, "MET_njetsge5_"+btagstring);

}


void studyRtt_0lep() {

  const  float exval=  gStyle->GetErrorX();

   loadSamples();

   usePUweight_=true;
   useHTeff_=true;
   useMHTeff_=true; // is this what I want????
   thebnnMHTeffMode_ = kOff;
   currentConfig_=configDescriptions_.getCorrected(); //add JERbias

   //   setSearchRegions();

//   double nsig = 118;
//   double nsb = 133;

//   double nsigerr=17;
//   double nsberr=16;

//   double r=nsig/nsb;
//   double rerr = jmt::errAoverB(nsig,nsigerr,nsb,nsberr);
//   cout<<"R_0l = "<<r<<" +/- "<<rerr<<endl;

  //plot R_SL as a function of run number in data

  // ??? plot R_0L as a function of run number in data

  // plot R_SL and R_0L as functions of nPV

   gStyle->SetOptStat(0);

  TCut baselineSL = "cutTrigger==1 && minDeltaPhiN>4 && passCleaning==1 && cutPV==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0 && MT_Wlep<100";
  TCut baselineEle = "cutTrigger==1 && minDeltaPhiN>4 && passCleaning==1 && cutPV==1 && cut3Jets==1 && (nElectrons==1 && nMuons==0) && MT_Wlep>=0 && MT_Wlep<100";
  TCut baselineMu = "cutTrigger==1 && minDeltaPhiN>4 && passCleaning==1 && cutPV==1 && cut3Jets==1 && (nElectrons==0 && nMuons==1) && MT_Wlep>=0 && MT_Wlep<100";
  TCut baseline = "cutTrigger==1 && minDeltaPhiN>4 && passCleaning==1 && cutPV==1 && cut3Jets==1 && nElectrons==0 && nMuons==0";
  TCut HTcut = "HT>400";
  TCut SIG="MET>250";
  TCut SB = "MET>200 && MET<=250";
  //  TCut bcut = "nbjets>=1";
  TCut bcut="1";      btagSFweight_="probge1"; //use b tag SF
  //TCut bcut="1";//no b tag

    float rmin=0.6; float rmax=1.1; //ge1b
  //float rmin=0.5; float rmax=1.5; //ge2b

  int nbins = 5;
  double runs[]={160000,166000,169000,175000,177500,181000};
/* this did not prove to show any compelling
  TH1D SL_SIG_vTime("SL_SIG_vTime","SL SIG versus run number",nbins,runs);
  TH1D SL_SB_vTime("SL_SB_vTime","SL SB versus run number",nbins,runs);
  TH1D SL_R_vTime("SL_R_vTime","R_SL versus run number",nbins,runs);

  SL_SIG_vTime.Sumw2();
  SL_SB_vTime.Sumw2();
  SL_R_vTime.Sumw2();

  selection_ = baselineSL && HTcut && SIG && bcut; //SIG
  dtree->Project("SL_SIG_vTime","runNumber",getCutString(true).Data());

  selection_ = baselineSL && HTcut && SB && bcut; //SB
  dtree->Project("SL_SB_vTime","runNumber",getCutString(true).Data());

  SL_R_vTime.Divide(&SL_SIG_vTime,&SL_SB_vTime);

  renewCanvas();
  SL_R_vTime.Draw();

  thecanvas->SaveAs("SL_R_vTime.eps");
  thecanvas->SaveAs("SL_R_vTime.png");
  thecanvas->SaveAs("SL_R_vTime.pdf");
*/

/*
  TH1D zL_SB_vTime("zL_SB_vTime","zL SB versus run number",nbins,runs);
  TH1D SL_SB_vTime("SL_SB_vTime","SL SB versus run number",nbins,runs);
  TH1D zSL_R_vTime("zSL_R_vTime","SB r_(0L/SL) versus run number",nbins,runs);

  zL_SB_vTime.Sumw2();
  SL_SB_vTime.Sumw2();
  zSL_R_vTime.Sumw2();

  selection_ = baseline && HTcut && SB && bcut; //SIG
  dtree->Project("zL_SB_vTime","runNumber",getCutString(true).Data());

  selection_ = baselineMu && HTcut && SB && bcut; //SB
  dtree->Project("SL_SB_vTime","runNumber",getCutString(true).Data());

  zSL_R_vTime.Divide(&zL_SB_vTime,&SL_SB_vTime);

  renewCanvas();
  zSL_R_vTime.Draw();

  thecanvas->SaveAs("zSL_r_vTime.eps");
  thecanvas->SaveAs("zSL_r_vTime.png");
  thecanvas->SaveAs("zSL_r_vTime.pdf");
*/

  //now do it by nPV
  TString var ="nGoodPV";
  nbins = 3;
  float npv[]={0,4,10,30};
  TH1D SL_SIG_nPV("SL_SIG_nPV","SL SIG versus nPV",nbins,npv);
  TH1D SL_SB_nPV("SL_SB_nPV","SL SB versus nPV",nbins,npv);
  TH1D SL_R_nPV("SL_R_nPV","R_SL versus nPV",nbins,npv);

  SL_SIG_nPV.Sumw2();
  SL_SB_nPV.Sumw2();
  SL_R_nPV.Sumw2();
  SL_R_nPV.SetMarkerStyle(4);

  savePlots_=false;
  drawMarkers_=false;

  selection_ = baselineSL && HTcut && SIG && bcut; //SIG
  dtree->Project("SL_SIG_nPV",var,getCutString(true).Data());
  drawPlots(var,nbins,0,30,var,"R_SL","dummy",npv);
  TH1D* SL_SIG_nPV_MC = (TH1D*) totalsm->Clone("SL_SIG_nPV_MC");

  selection_ = baselineSL && HTcut && SB && bcut; //SB
  dtree->Project("SL_SB_nPV",var,getCutString(true).Data());
  drawPlots(var,nbins,0,30,var,"R_SL","dummy",npv);
  TH1D* SL_SB_nPV_MC = (TH1D*) totalsm->Clone("SL_SB_nPV_MC");

  SL_R_nPV.Divide(&SL_SIG_nPV,&SL_SB_nPV);
  TH1D*  SL_R_nPV_MC= (TH1D*) SL_SB_nPV_MC->Clone("SL_R_nPV_MC");
  SL_R_nPV_MC->Reset(); //prob not needed
  SL_R_nPV_MC->Divide(SL_SIG_nPV_MC,SL_SB_nPV_MC);

  //ok, now do it all again for the SL sample WITH ONLY tt+W+t
  clearSamples();
  addSample("TTbarJets");
  addSample("WJets");
  addSample("SingleTop");

  selection_ = baselineSL && HTcut && SIG && bcut; //SIG
  drawPlots(var,nbins,0,30,var,"R_SL","dummy",npv);
  TH1D* SL_ttWt_SIG_nPV_MC = (TH1D*) totalsm->Clone("SL_ttWt_SIG_nPV_MC");

  selection_ = baselineSL && HTcut && SB && bcut; //SB
  drawPlots(var,nbins,0,30,var,"R_SL","dummy",npv);
  TH1D* SL_ttWt_SB_nPV_MC = (TH1D*) totalsm->Clone("SL_ttWt_SB_nPV_MC");

  TH1D*  SL_ttWt_R_nPV_MC= (TH1D*) SL_ttWt_SB_nPV_MC->Clone("SL_ttWt_R_nPV_MC");
  SL_ttWt_R_nPV_MC->Reset(); //prob not needed
  SL_ttWt_R_nPV_MC->Divide(SL_ttWt_SIG_nPV_MC,SL_ttWt_SB_nPV_MC);


  renewCanvas();
  SL_R_nPV.SetLineColor(kBlack);
  SL_R_nPV_MC->SetLineColor(sampleColor_["TotalSM"]);
  SL_ttWt_R_nPV_MC->SetLineColor(sampleColor_["TTbarJets"]);
  gStyle->SetErrorX(exval);

  SL_R_nPV_MC->SetMaximum(rmax);
  SL_R_nPV_MC->SetMinimum(rmin);
  SL_R_nPV_MC->Draw();
  SL_ttWt_R_nPV_MC->Draw("SAME");
  SL_R_nPV.Draw("same");

 
  thecanvas->SaveAs("SL_R_nPV.eps");
  thecanvas->SaveAs("SL_R_nPV.png");
  thecanvas->SaveAs("SL_R_nPV.pdf");
 return;

  //ok, now do it all again for the 0 lepton sample
  resetSamples(); //all samples
  TH1D zL_SIG_nPV("zL_SIG_nPV","0 lep SIG versus nPV",nbins,npv);
  TH1D zL_SB_nPV("zL_SB_nPV","0 lep SB versus nPV",nbins,npv);
  TH1D zL_R_nPV("zL_R_nPV","R_zl versus nPV",nbins,npv);

  zL_SIG_nPV.Sumw2();
  zL_SB_nPV.Sumw2();
  zL_R_nPV.Sumw2();
  //  zL_R_nPV.SetMarkerStyle(4);

  selection_ = baseline && HTcut && SIG && bcut; //SIG
  dtree->Project("zL_SIG_nPV",var,getCutString(true).Data());
  drawPlots(var,nbins,0,30,var,"R_zL","dummy",npv);
  TH1D* zL_SIG_nPV_MC = (TH1D*) totalsm->Clone("zL_SIG_nPV_MC");

  selection_ = baseline && HTcut && SB && bcut; //SB
  dtree->Project("zL_SB_nPV",var,getCutString(true).Data());
  drawPlots(var,nbins,0,30,var,"R_zL","dummy",npv);
  TH1D* zL_SB_nPV_MC = (TH1D*) totalsm->Clone("zL_SB_nPV_MC");

  zL_R_nPV.Divide(&zL_SIG_nPV,&zL_SB_nPV);
  TH1D*  zL_R_nPV_MC= (TH1D*) zL_SB_nPV_MC->Clone("zL_R_nPV_MC");
  zL_R_nPV_MC->Reset(); //prob not needed
  zL_R_nPV_MC->Divide(zL_SIG_nPV_MC,zL_SB_nPV_MC);

  renewCanvas();
  zL_R_nPV.SetLineColor(kBlack);
  zL_R_nPV_MC->SetLineColor(sampleColor_["TotalSM"]);
  //  zL_R_nPV_MC->SetMarkerColor(sampleColor_["TotalSM"]);

  //ok, now do it all again for the 0 lepton sample WITH ONLY tt+W+t
  clearSamples();
  addSample("TTbarJets");
  addSample("WJets");
  addSample("SingleTop");
//   TH1D zL_ttWt_SIG_nPV("zL_ttWt_SIG_nPV","0 lep SIG (ttWt) versus nPV",nbins,npv);
//   TH1D zL_ttWt_SB_nPV("zL_ttWt_SB_nPV","0 lep SB (ttWt) versus nPV",nbins,npv);
//   TH1D zL_ttWt_R_nPV("zL_ttWt_R_nPV","R_zl (ttWt) versus nPV",nbins,npv);

//   zL_ttWt_SIG_nPV.Sumw2();
//   zL_ttWt_SB_nPV.Sumw2();
//   zL_ttWt_R_nPV.Sumw2();
  //  zL_ttWt_R_nPV.SetMarkerStyle(4);

  selection_ = baseline && HTcut && SIG && bcut; //SIG
  drawPlots(var,nbins,0,30,var,"R_zL","dummy",npv);
  TH1D* zL_ttWt_SIG_nPV_MC = (TH1D*) totalsm->Clone("zL_ttWt_SIG_nPV_MC");

  selection_ = baseline && HTcut && SB && bcut; //SB
  drawPlots(var,nbins,0,30,var,"R_zL","dummy",npv);
  TH1D* zL_ttWt_SB_nPV_MC = (TH1D*) totalsm->Clone("zL_ttWt_SB_nPV_MC");

  TH1D*  zL_ttWt_R_nPV_MC= (TH1D*) zL_ttWt_SB_nPV_MC->Clone("zL_ttWt_R_nPV_MC");
  zL_ttWt_R_nPV_MC->Reset(); //prob not needed
  zL_ttWt_R_nPV_MC->Divide(zL_ttWt_SIG_nPV_MC,zL_ttWt_SB_nPV_MC);

  //now draw

  renewCanvas();
  zL_ttWt_R_nPV_MC->SetLineColor(sampleColor_["TTbarJets"]);

  gStyle->SetErrorX(exval);

  zL_R_nPV_MC->SetMaximum(rmax);
  zL_R_nPV_MC->SetMinimum(rmin);
  zL_R_nPV_MC->Draw();
  zL_ttWt_R_nPV_MC->Draw("same");
  zL_R_nPV.Draw("same");

  thecanvas->SaveAs("zL_R_nPV.eps");
  thecanvas->SaveAs("zL_R_nPV.png");
  thecanvas->SaveAs("zL_R_nPV.pdf");

  resetSamples();
  //now draw data/data comparison of nPv for 0L and SL samples
  //add SL SB and SIG together to get pv distribution in whole SL sample
  TH1D* SL_SIGSB_nPV = (TH1D*) SL_SIG_nPV.Clone("SL_SIGSB_nPV");
  SL_SIGSB_nPV->Add(&SL_SB_nPV); 

  selection_ = baseline && HTcut && bcut && TCut("MET>200"); //SIG+SB
  TH1D* zL_SIGSB_nPV = (TH1D*) SL_SIGSB_nPV->Clone("zL_SIGSB_nPV");
  zL_SIGSB_nPV->Reset();
  dtree->Project("zL_SIGSB_nPV",var,getCutString(true).Data());

  SL_SIGSB_nPV->Scale(1.0 / SL_SIGSB_nPV->Integral());
  zL_SIGSB_nPV->Scale(1.0 / zL_SIGSB_nPV->Integral());

  renewCanvas();
  zL_SIGSB_nPV->SetLineColor(kBlack);
  //  zL_SIGSB_nPV->SetMarkerColor(kBlack);
  SL_SIGSB_nPV->SetLineColor(kBlue);
  //  SL_SIGSB_nPV->SetMarkerColor(kBlue);
  //  SL_SIGSB_nPV->SetMarkerStyle(25);
  //  zL_SIGSB_nPV->SetMarkerStyle(4);

  gStyle->SetErrorX(exval);
  zL_SIGSB_nPV->SetMinimum(0.1);
  zL_SIGSB_nPV->SetMaximum(0.7);

  zL_SIGSB_nPV->Draw();
  SL_SIGSB_nPV->Draw("SAME");

  thecanvas->SaveAs("nPV_SL0L_data.eps");
  thecanvas->SaveAs("nPV_SL0L_data.png");
  thecanvas->SaveAs("nPV_SL0L_data.pdf");

//   TH1D SL_SIG_nPV("SL_SIG_nPV","SL SIG versus nPV",nbins,npv);
//   TH1D SL_SB_nPV("SL_SB_nPV","SL SB versus nPV",nbins,npv);
//   TH1D SL_R_nPV("SL_R_nPV","R_SL versus nPV",nbins,npv);
  //now do a quick rescaling of R_SL

  double n_sig_prime=0;
  double n_sb_prime=0;

  //for x check
  double n_sig=0;
  double n_sb=0;

  for (int ii=0; ii<3; ii++) {
    n_sig_prime += SL_SIG_nPV.GetBinContent(ii+1) * zL_SIGSB_nPV->GetBinContent(ii+1)/SL_SIGSB_nPV->GetBinContent(ii+1);
    n_sb_prime += SL_SB_nPV.GetBinContent(ii+1) * zL_SIGSB_nPV->GetBinContent(ii+1)/SL_SIGSB_nPV->GetBinContent(ii+1);

    n_sig += SL_SIG_nPV.GetBinContent(ii+1);
    n_sb += SL_SB_nPV.GetBinContent(ii+1);
  }

  cout<<"R_SL nominal  = "<<n_sig/n_sb<<endl;
  cout<<"R_SL reweight = "<<n_sig_prime/n_sb_prime<<endl;

}


void drawTrigEff() {

  //loadSamples();

  int nbins;
  float low,high;
  TString var,xtitle;
  
  doOverflowAddition(true);

  setStackMode(false,false); //normalized
  //setColorScheme("nostack");
  clearSamples();
  //addSample("PythiaPUQCD");
  //addSample("TTbarJets");
  //addSample("LM9");
  drawLegend(false);
  doData(true);


  //output file 
  //TString histfilename = "mhteff_dec4.root";
  TString histfilename = "mhteff_dec8_eq0b.root";
  //TString histfilename = "mhteff_dec8_eq0b_HT500.root";
  //TString histfilename = "mhteff_dec8_eq0b_HT600.root";
  //TString histfilename = "mhteff_dec8_ge1b.root";
  TFile fh(histfilename,"RECREATE");//will delete old root file if present 
  fh.Close(); //going to open only when needed 

  //TCut HTcut = "HT>=400";
  //TCut HTcut = "HT>=400";  TCut njetCut = "njets >=3";    TCut dpcut = "minDeltaPhiN>=4";
  //TCut HTcut = "HT>=400";  TCut njetCut = "nbjetsCSVM>=1";    TCut dpcut = "";
  TCut HTcut = "HT>=600";  TCut njetCut = "nbjetsCSVM==0";    TCut dpcut = "";
  //TCut dpcut = "";
  
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 100; low=0; high=1000;
  setLogY(false); resetPlotMinimum();


  /////////////////
  // HT260_MHT60_v2
  /////////////////
  //all
  selection_ =TCut("pass_utilityHLT==1   && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2", "data");
  //all - pass cleaning
  selection_ =TCut("pass_utilityHLT==1   && passCleaning==1 && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_clean_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && passCleaning==1 && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_clean", "data");
  //SL
  selection_ =TCut("pass_utilityHLT==1   && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_SL_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_SL", "data");
  //NoLep
  selection_ =TCut("pass_utilityHLT==1   && cutEleVeto==1 && cutMuVeto==1 && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_NoLep_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && cutEleVeto==1 && cutMuVeto==1 && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_NoLep", "data");
  //highdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_HDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_HDP", "data");
  //ldp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN < 4 && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_LDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN < 4 && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_LDP", "data");
  //SL - hdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_SLHDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_SLHDP", "data");
  //NoLep - hdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_NoLepHDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_NoLepHDP", "data");
  //NoLep - ldp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN < 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_NoLepLDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN < 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=160431 && runNumber<=161204")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht260mht60v2_NoLepLDP", "data");

  /////////////////
  // HT250_MHT60_v*
  /////////////////
  //all
  selection_ =TCut("pass_utilityHLT==1   && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3", "data");
  //all - pass cleaning
  selection_ =TCut("pass_utilityHLT==1   && passCleaning==1 && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_clean_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && passCleaning==1 && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_clean", "data");
  //SL
  selection_ =TCut("pass_utilityHLT==1   && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_SL_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_SL", "data");
  //NoLep
  selection_ =TCut("pass_utilityHLT==1   && cutEleVeto==1 && cutMuVeto==1 && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_NoLep_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && cutEleVeto==1 && cutMuVeto==1 && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_NoLep", "data");
  //highdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_HDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_HDP", "data");
  //ldp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN < 4 && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_LDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN < 4 && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_LDP", "data");
  //SL - hdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_SLHDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_SLHDP", "data");
  //NoLep - hdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_NoLepHDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_NoLepHDP", "data");
  //NoLep - ldp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN < 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_NoLepLDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN < 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=161205 && runNumber<=164923")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht60v2and3_NoLepLDP", "data");


  /////////////////
  // HT250_MHT70_v1
  /////////////////
  //all
  selection_ =TCut("pass_utilityHLT==1   && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1", "data");
  //all - pass cleaning
  selection_ =TCut("pass_utilityHLT==1   && passCleaning==1 && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_clean_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && passCleaning==1 && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_clean", "data");
  //SL
  selection_ =TCut("pass_utilityHLT==1   && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_SL_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_SL", "data");
  //NoLep
  selection_ =TCut("pass_utilityHLT==1   && cutEleVeto==1 && cutMuVeto==1 && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_NoLep_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && cutEleVeto==1 && cutMuVeto==1 && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_NoLep", "data");
  //highdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_HDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_HDP", "data");
  //ldp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN < 4 && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_LDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN < 4 && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_LDP", "data");
  //SL - hdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_SLHDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_SLHDP", "data");
  //NoLep - hdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_NoLepHDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_NoLepHDP", "data");
  //NoLep - ldp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN < 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_NoLepLDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN < 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=164924 && runNumber<=165921")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht250mht70v1_NoLepLDP", "data");


  /////////////////
  // HT300_MHT75_v*
  /////////////////
  //all
  selection_ =TCut("pass_utilityHLT==1   && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8", "data");
  //all - pass cleaning
  selection_ =TCut("pass_utilityHLT==1   && passCleaning==1 && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_clean_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && passCleaning==1 && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_clean", "data");
  //SL
  selection_ =TCut("pass_utilityHLT==1   && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_SL_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_SL", "data");
  //NoLep
  selection_ =TCut("pass_utilityHLT==1   && cutEleVeto==1 && cutMuVeto==1 && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_NoLep_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && cutEleVeto==1 && cutMuVeto==1 && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_NoLep", "data");
  //highdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_HDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_HDP", "data");
  //ldp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN < 4 && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_LDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN < 4 && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_LDP", "data");
  //SL - hdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_SLHDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_SLHDP", "data");
  //NoLep - hdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_NoLepHDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_NoLepHDP", "data");
  //NoLep - ldp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN < 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_NoLepLDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN < 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=165922 && runNumber<=166978")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht75v7and8_NoLepLDP", "data");


  /////////////////
  // HT300_MHT80_v*
  /////////////////
  //all
  selection_ =TCut("pass_utilityHLT==1   && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2", "data");
  //all - pass cleaning
  selection_ =TCut("pass_utilityHLT==1   && passCleaning==1 && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_clean_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && passCleaning==1 && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_clean", "data");
  //SL
  selection_ =TCut("pass_utilityHLT==1   && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_SL_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_SL", "data");
  //NoLep
  selection_ =TCut("pass_utilityHLT==1   && cutEleVeto==1 && cutMuVeto==1 && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_NoLep_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && cutEleVeto==1 && cutMuVeto==1 && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_NoLep", "data");
  //highdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_HDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_HDP", "data");
  //ldp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN < 4 && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_LDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN < 4 && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_LDP", "data");
  //SL - hdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_SLHDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_SLHDP", "data");
  //NoLep - hdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_NoLepHDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_NoLepHDP", "data");
  //NoLep - LDP
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN < 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_NoLepLDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN < 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=166979 && runNumber<=173211")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht80v1and2_NoLepLDP", "data");



  /////////////////
  // HT300_MHT90_v2
  /////////////////
  //all
  selection_ =TCut("pass_utilityHLT==1   && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_den", "data");
  //drawPlots(var,nbins,low,high,xtitle,"Events", "MET_denom");
  //TH1D* DENOMplot = (TH1D*)hinteractive->Clone("DENOMplot");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2", "data");
  //drawPlots(var,nbins,low,high,xtitle,"Events", "MET_num");
  //all -passcleaning
  selection_ =TCut("pass_utilityHLT==1   && passCleaning==1 && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_clean_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && passCleaning==1 && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_clean", "data");
  //SL
  selection_ =TCut("pass_utilityHLT==1   && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_SL_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_SL", "data");
  //NoLep
  selection_ =TCut("pass_utilityHLT==1   && cutEleVeto==1 && cutMuVeto==1 && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_NoLep_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && cutEleVeto==1 && cutMuVeto==1 && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_NoLep", "data");
  //highdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_HDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_HDP", "data");
  //ldp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN < 4 && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_LDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN < 4 && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_LDP", "data");
  //SL-HDP
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_SLHDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_SLHDP", "data");
  //NoLep-HDP
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_NoLepHDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_NoLepHDP", "data");
  //NoLep-LDP
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN < 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_NoLepLDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN < 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=173212 && runNumber<=176544")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht300mht90v2_NoLepLDP", "data");


  /////////////////
  // HT350_MHT90_v1
  /////////////////
  //all
  selection_ =TCut("pass_utilityHLT==1   && runNumber>=176545 && runNumber<=178410")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && runNumber>=176545 && runNumber<=178410")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1", "data");
  //all -pass cleaning
  selection_ =TCut("pass_utilityHLT==1   && passCleaning==1 && runNumber>=176545 && runNumber<=178410")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_clean_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && passCleaning==1 && runNumber>=176545 && runNumber<=178410")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_clean", "data");
  //SL
  selection_ =TCut("pass_utilityHLT==1   && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=176545 && runNumber<=178410")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_SL_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=176545 && runNumber<=178410")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_SL", "data");
  //NoLep
  selection_ =TCut("pass_utilityHLT==1   && cutEleVeto==1 && cutMuVeto==1 && runNumber>=176545 && runNumber<=178410")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_NoLep_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && cutEleVeto==1 && cutMuVeto==1 && runNumber>=176545 && runNumber<=178410")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_NoLep", "data");
  //highdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && runNumber>=176545 && runNumber<=178410 ")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_HDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && runNumber>=176545 && runNumber<=178410")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_HDP", "data");
  //ldp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN < 4 && runNumber>=176545 && runNumber<=178410 ")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_LDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN < 4 && runNumber>=176545 && runNumber<=178410")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_LDP", "data");
  //SL -HDP
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=176545 && runNumber<=178410")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_SLHDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=176545 && runNumber<=178410")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_SLHDP", "data");
  //NoLep -HDP
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=176545 && runNumber<=178410")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_NoLepHDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=176545 && runNumber<=178410")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_NoLepHDP", "data");
  //NoLep -LDP
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN < 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=176545 && runNumber<=178410")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_NoLepLDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN < 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=176545 && runNumber<=178410")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht90v1_NoLepLDP", "data");


  //////////////////
  // HT350_MHT110_v3
  //////////////////
  //all
  selection_ =TCut("pass_utilityHLT==1   && runNumber>=178411")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && runNumber>=178411")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3", "data");
  //all -passcleaning
  selection_ =TCut("pass_utilityHLT==1   && passCleaning==1 && runNumber>=178411")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_clean_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && passCleaning==1 && runNumber>=178411")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_clean", "data");
  //SL
  selection_ =TCut("pass_utilityHLT==1   && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=178411")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_SL_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=178411")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_SL", "data");
  //NoLep
  selection_ =TCut("pass_utilityHLT==1   && cutEleVeto==1 && cutMuVeto==1 && runNumber>=178411")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_NoLep_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && cutEleVeto==1 && cutMuVeto==1 && runNumber>=178411")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_NoLep", "data");
  //highdp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && runNumber>=178411 ")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_HDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && runNumber>=178411")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_HDP", "data");
  //ldp
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN < 4 && runNumber>=178411 ")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_LDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN < 4 && runNumber>=178411")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_LDP", "data");
  //SL -HDP
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=178411")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_SLHDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && runNumber>=178411")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_SLHDP", "data");
  //NoLep -HDP
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN >= 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=178411")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_NoLepHDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN >= 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=178411")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_NoLepHDP", "data");
  //NoLep -LDP
  selection_ =TCut("pass_utilityHLT==1   && minDeltaPhiN < 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=178411")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_NoLepLDP_den", "data");
  selection_ =TCut("cutTrigger==1 && pass_utilityHLT==1  && minDeltaPhiN < 4 && cutEleVeto==1 && cutMuVeto==1 && runNumber>=178411")&&njetCut&&dpcut&&HTcut;
  drawSimple(var,nbins,low,high,histfilename, "MET_ht350mht110v3_NoLepLDP", "data");



  //TH1D* NUMplot = (TH1D*)hinteractive->Clone("NUMplot");

  ////TGraphAsymmErrors * h_eff = new TGraphAsymmErrors(NUMplot,DENOMplot,"cl=0.683 b(1,1) mode");
  //
  //TH1D* h_eff = (TH1D*)NUMplot->Clone("h_num");
  //h_eff->Sumw2();
  //h_eff->Divide(h_eff,DENOMplot,1,1,"B"); 
  //
  //TCanvas* myC = 0;
  //myC = new TCanvas("myC", "myC", 600,700);
  //gStyle->SetPadBorderMode(0);
  //gStyle->SetFrameBorderMode(0);
  //
  //h_eff->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]"); //make y label bigger
  //
  ////h_eff->Draw("AP");
  //h_eff->Draw();
  

  return;
}
