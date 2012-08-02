/*
to compile:

must have symlink to MiscUtil.cxx in the working directory.
TSelectorMultiDraw and CrossSectionTable need to be compiled only when they are changed.

gSystem->Load("TSelectorMultiDraw_C.so");
gSystem->Load("CrossSectionTable_cxx.so");
gSystem->Load("ConfigurationDescriptions_cxx.so");
gSystem->Load("SearchRegion_cxx.so");
gSystem->Load("SystInfo_cxx.so");
gSystem->Load("SignalEffData_cxx.so");
.L drawReducedTrees2012.C++

====== this is the nominal code for drawing RA2b data/MC comparison plots =======
The signal efficiency functions are now forked into signalEffSyst.C

2011 updates -- runDataQCD2011() computes all numbers needed for QCD estimate including systematics

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
	
TString inputPath = "/cu2/ra2b/reducedTrees/v63_0/"; //2012
TString dataInputPath =  "/cu2/ra2b/reducedTrees/v63_0/ORIGINALS/";
//TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-35u/Fall11/";//7 TeV

double lumiScale_ = 4301;//2012 lumi ????
double preLumiScale_ = 30;//god only knows for 2012

//make a symlink that point from this name to drawReducedTree.h
//this is to make the ROOT dictionary generation work correctly
#include "drawReducedTrees2012.h"

//2012 -- remove whatever of these are now crap
const bool reweightLSBdata_=true; //whether or not LSB data is reweighted based on PV distribution
const bool useScaleFactors_=true; //whether or not to use MC scale factors when doing subtraction for data-driven estimates
const bool useBNNEffCurves_=false; 
const bool btaggedLSB_=false;
const bool use1B_SL_=false; // use ge1b selection for SL sample in ttbar method
const bool doNJetRWClosure_ = true; //should usually be true; set it to false only to save time

void RA2bControlPlots2012() {

  nameOfEventWeight_="weight2"; //for 2012 cfA skimmed ntuples

  loadSamples();
  usePUweight_=true;

  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getDefault();

  resetSamples();
  removeSample("LM9");

  //removeSample("VV");

  int nbins;
  float low,high;
  TString var,xtitle;
  
  doOverflowAddition(true);

  dodata_=true;

  // == ttbar control sample --> 1BL
  setDatasetToDraw("HTMHT");
  lumiScale_=5100;
  selection_ =TCut("HT>=400 && cutPV==1 && pass_PFHT350_PFMET100==1  && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1")&&getSingleLeptonCut();
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 30; low=100; high=400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_1BL_SL",0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_1BL_SL",0,"GeV");

  //HT
  setDatasetToDraw("HTMHT");
  lumiScale_= 5100;
  selection_ =TCut("HT>=300 && cutPV==1 && pass_PFHT350_PFMET100==1  && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1 && MET>=175")&&getSingleLeptonCut();
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 22; low=300; high=1400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_1BL_SL",0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_1BL_SL",0,"GeV");

  // == W control sample --> nb == 0, 1 lepton
  setDatasetToDraw("HTMHT");
  lumiScale_=5100;
  selection_ =TCut("HT>=400 && cutPV==1 && pass_PFHT350_PFMET100==1  && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0")&&getSingleLeptonCut();
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 30; low=100; high=400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_eq0b_SL",0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_eq0b_SL",0,"GeV");

  //HT
  setDatasetToDraw("HTMHT");
  lumiScale_=5100;
  selection_ =TCut("HT>=300 && MET>=175 && cutPV==1 && pass_PFHT350_PFMET100==1  && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0")&&getSingleLeptonCut();
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 22; low=300; high=1400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_eq0b_SL",0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_eq0b_SL",0,"GeV");

  //nPV -- just as a check
  setDatasetToDraw("HTMHT");
  lumiScale_=5100;
  selection_ =TCut("HT>=400 && cutPV==1 && pass_PFHT350_PFMET100==1 && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && MET>=175 &&nbjets>=1")&&getSingleLeptonCut();
  var="nGoodPV"; xtitle="num of good PV";
  nbins = 20; low=0; high=40;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_npv_SL",0);

  //check non-ttbar
  setDatasetToDraw("HTMHT");
  lumiScale_=5100;
  selection_ =TCut("HT>=400 && cutPV==1 && pass_PFHT350_PFMET100==1 && njets>=2 && passCleaning==1 && MET>=160 &&nbjets==0")&&getSingleLeptonCut();
  var="nGoodPV"; xtitle="num of good PV";
  nbins = 20; low=0; high=40;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_npv_SL_eq0b",0);

  //njets
  setDatasetToDraw("HTMHT");
  lumiScale_=5100;
  selection_ =TCut("HT>=400 && cutPV==1 && pass_PFHT350_PFMET100==1 && minDeltaPhiN >= 4 && passCleaning==1 && MET>=175 &&nbjets>=1")&&getSingleLeptonCut();
  var="njets"; xtitle="n jets (pT>50 GeV)";
  nbins = 10; low=0; high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets_met175_ht400_SL",0);

  //nbjets
  setDatasetToDraw("HTMHT");
  lumiScale_=5100;
  selection_ =TCut("HT>=400 && cutPV==1 && pass_PFHT350_PFMET100==1 && minDeltaPhiN >= 4 && passCleaning==1 && MET>=175 &&njets>=3")&&getSingleLeptonCut();
  var="nbjets"; xtitle="n CSVM b jets (pT>30 GeV)";
  nbins = 5; low=0; high=5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_nbjets_met175_ht400_SL",0);

  //minDeltaPhiN
  setDatasetToDraw("HTMHT");
  lumiScale_=5100;
  selection_ =TCut("HT>=400 && cutPV==1 && pass_PFHT350_PFMET100==1 && passCleaning==1 && MET>=175 &&nbjets >=1 &&njets>=3")&&getSingleLeptonCut();
  var="minDeltaPhiN"; xtitle="minDeltaPhiN";
  nbins = 20; low=0; high=40;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhiN_met175_ht400_SL",0);

  // same thing but integrated over b-tags
  setDatasetToDraw("HTMHT");
  lumiScale_=5100;
  selection_ =TCut("HT>=400 && cutPV==1 && pass_PFHT350_PFMET100==1 && passCleaning==1 && MET>=175 &&njets>=3")&&getSingleLeptonCut();
  var="minDeltaPhiN"; xtitle="minDeltaPhiN";
  nbins = 20; low=0; high=40;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhiN_met175_ht400_SL_ge0b",0);

  //muon pT
  setDatasetToDraw("HTMHT");
  lumiScale_=5100;
  selection_ =TCut("HT>=400 && cutPV==1 && pass_PFHT350_PFMET100==1 && passCleaning==1 && MET>=175 &&nbjets >=1 &&njets>=3")&&getSingleMuonCut();
  var="muonpt1"; xtitle="muon pT";
  nbins = 19; low=10; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_muonpt_met175_ht400_1m",0);

  // == MuHad sample
  //pass_PFHT350_Mu15_PFMET45
//   setDatasetToDraw("MuHad");
//   lumiScale_=570+155+929+824+1274;
//   selection_ =TCut("HT>=400 && cutPV==1 && pass_PFHT350_Mu15_PFMET45==1  && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1")&&getSingleMuonCut()&&TCut("muonpt1>=20");
//   var="MET"; xtitle="E_{T}^{miss} [GeV]";
//   nbins = 40; low=0; high=400;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_ge1b_1mu20",0,"GeV");

//   //DY should not have so much MET, but we can try
//   setDatasetToDraw("MuHad");
//   lumiScale_=570+155+929+824+1274;
//   selection_ =TCut("HT>=400 && cutPV==1 && pass_PFHT350_Mu15_PFMET45==1  && njets>=2 && passCleaning==1 && MET >=80")&&TCut("muonpt1>=20 && nMuons5==2");
//   var="bestZmass"; xtitle="l+l- mass";
//   nbins = 20; low=0; high=220;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_Zmass_1mu20",0,"GeV");

  // == LDP sample!
  setDatasetToDraw("HTMHT");
  lumiScale_=5100;
  selection_ =TCut("HT>=400 && cutPV==1 && pass_PFHT350_PFMET100==1  && njets>=3  && minDeltaPhiN < 4 && passCleaning==1 && nbjets>=1")&&getLeptonVetoCut();
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 25; low=150; high=400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_1BL_LDP",0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_1BL_LDP",0,"GeV");

  // == 0b sample
  setDatasetToDraw("HTMHT");
  lumiScale_=5100;
  selection_ =TCut("HT>=400 && cutPV==1 && pass_PFHT350_PFMET100==1  && njets>=3  && MET>=175 && passCleaning==1 && nbjets==0")&&getLeptonVetoCut();
  var="minDeltaPhiN"; xtitle="minDeltaPhiN";
  nbins = 20; low=0; high=40;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_mindp_met175_eq0b",0);

  //LDP -- njets b-veto
  setDatasetToDraw("HTMHT");
  lumiScale_=5100;
  selection_ =TCut("HT>=400 && cutPV==1 && pass_PFHT350_PFMET100==1 && MET>=175 && passCleaning==1 && minDeltaPhiN<4 && nbjets==0")&&getLeptonVetoCut();
  var="njets"; xtitle="n jets (pT>50 GeV)";
  nbins = 10; low=0; high=10;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets_met175_LDP_eq0b",0);

  //same thing but ge1b
  setDatasetToDraw("HTMHT");
  lumiScale_=5100;
  selection_ =TCut("HT>=400 && cutPV==1 && pass_PFHT350_PFMET100==1 && MET>=175 && passCleaning==1 && minDeltaPhiN<4 && nbjets>=1")&&getLeptonVetoCut();
  var="njets"; xtitle="n jets (pT>50 GeV)";
  nbins = 10; low=0; high=10;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets_met175_LDP_ge1b",0);

  //LDP nb
  setDatasetToDraw("HTMHT");
  lumiScale_=5100;
  selection_ =TCut("HT>=400 && cutPV==1 && pass_PFHT350_PFMET100==1 && MET>=175 && passCleaning==1 && minDeltaPhiN<4 && njets>=2")&&getLeptonVetoCut();
  var="nbjets"; xtitle="n b jets (pT>30 GeV)";
  nbins = 5; low=0; high=5;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_nbjets_met175_LDP",0);


  //Single Mu -- can look at MET and HT
  setDatasetToDraw("SingleMu");
  lumiScale_ = 5100;

  //HLT_IsoMu24 without eta cut does not exist in 2012A
  selection_ =TCut("MET>=100 && HT>=200 && njets>=2 && cutPV==1 && pass_IsoMu24_eta2p1==1 && passCleaning==1")&&TCut("nMuons==1 && nElectrons==0 && muonpt1>=28 && abs(muoneta1)<2.1");
  var="MET"; xtitle="MET";
  nbins = 20; low=100; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_singlemu_ge0b",0);

 //  setDatasetToDraw("SingleMu");
//   selection_ =TCut("HT>=300 && cutPV==1 && pass_IsoMu24_eta2p1==1 && MET>150&&  passCleaning==1")&&TCut("nMuons==1 && nElectrons==0 && abs(muoneta1)<2.1");
//   var="muonpt1"; xtitle="muon pt";
//   nbins = 40; low=10; high=200;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_muonpt_singlemu_highmet_ge0b",0);

//   setDatasetToDraw("SingleMu");
//   selection_ =TCut("HT>=300 && cutPV==1 && pass_IsoMu24_eta2p1==1 && MET>150&&  passCleaning==1 && nbjets>=1")&&TCut("nMuons==1 && nElectrons==0 && abs(muoneta1)<2.1");
//   var="muonpt1"; xtitle="muon pt";
//   nbins = 40; low=10; high=200;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_muonpt_singlemu_highmet_ge1b",0);


  //inverted MET cut
//   setDatasetToDraw("SingleMu");
//   selection_ =TCut("HT>=300 && cutPV==1 && pass_IsoMu24_eta2p1==1 && MET<150&&  passCleaning==1 && nbjets>=1")&&TCut("nMuons==1 && nElectrons==0 && abs(muoneta1)<2.1");
//   var="muonpt1"; xtitle="muon pt";
//   nbins = 40; low=10; high=200;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_muonpt_singlemu_lowmet_ge1b",0);

//   //try the loosest I can get
//   setDatasetToDraw("SingleMu");
//   selection_ =TCut("STeff>=300 && cutPV==1 && pass_IsoMu24_eta2p1==1 &&passCleaning==1")&&TCut("nMuons==1 && muonpt1 >29 && abs(muoneta1)<2.1");
//   var="njets30"; xtitle="njets30";
//   nbins = 8; low=0; high=8;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets30_singlemu",0);

//   selection_ =TCut("STeff>=300 && cutPV==1 && pass_IsoMu24_eta2p1==1 &&passCleaning==1 && nbjets>=1")&&TCut("nMuons==1 && muonpt1 >29 && abs(muoneta1)<2.1");
//   var="njets30"; xtitle="njets30";
//   nbins = 8; low=0; high=8;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets30_singlemu_ge1b",0);


  //try requiring 2 leptons
  setDatasetToDraw("DoubleMu");
  lumiScale_ = 5100;

  //I am using the DY samples with an HT cut, so need an HT cut here
  selection_ =TCut("cutPV==1 && pass_Mu17_Mu8==1 && passCleaning==1 && HT>=250")&&TCut("muonpt1>=20 && nMuons==2 && muonpt2>=10");
  var="bestZmass"; xtitle="l+l- mass";
  nbins = 50; low=0; high=200;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_Zmass_2mu",0,"GeV");
  setLogY(true);
  setPlotMinimum(1);
  nbins = 100; low=0; high=220;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_Zmass_2mu",0,"GeV");
  resetPlotMinimum();

  //now add b tag
  selection_ =TCut("cutPV==1 && pass_Mu17_Mu8==1 && passCleaning==1 && HT>=250 && nbjets>=1")&&TCut("muonpt1>=20 && nMuons==2 && muonpt2>=10");
  var="bestZmass"; xtitle="l+l- mass";
  nbins = 50; low=0; high=200;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_Zmass_2mu_ge1b",0,"GeV");
  setLogY(true);
  setPlotMinimum(1);
  nbins = 100; low=0; high=220;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_Zmass_2mu_ge1b",0,"GeV");
  resetPlotMinimum();

  //FUN!
  doOverflowAddition(false);
  selection_ =TCut("cutPV==1 && pass_Mu17_Mu8==1")&&TCut("muonpt1>=20 && nMuons==2 && muonpt2>=10");
  var="bestZmass"; xtitle="l+l- mass";
  nbins = 1000; low=2; high=12;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_Zmass_2mu_low",0,"GeV");
 
}


void compareZMC() {

  //would like to compare the shapes of the Z MCs

  /*
ok, after pondering this, I think that I cannot use regular PF MET when making these comparisons
(and I must cut on MET because of the skim)

I need to use the "modified" MET that comes from recalculating MET after dropping the leptons from Z->ll
  */

  nameOfEventWeight_="weight2"; //for 2012 cfA skimmed ntuples

  loadSamples();
  usePUweight_=false;

  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getDefault();

  clearSamples();
  addSample("ZJets");
  addSample("Zinvisible");

  int nbins;
  float low,high;
  TString var,xtitle;
  
  doOverflowAddition(true);

  dodata_=false;
  setStackMode(false,true);

  // == check HT turn-on
  clearSamples();
  addSample("ZJets");
  addSample("ZJetsInc");

  selection_ = "MET>=120 && njets>=2 && HT>=200";
  nbins = 50; low = 200; high= 1000;
  var = "HT"; xtitle="HT";
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"", "ra2b_DYcomp_ht",0,"GeV");

  // == use inclusive sample
  clearSamples();
  addSample("Zinvisible");
  addSample("ZJetsInc");

  selection_ = "MET>=120 && njets>=2 && HT>=200";
  nbins = 50; low = 200; high= 1000;
  var = "HT"; xtitle="HT";
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"", "ra2b_DYvZinv_ht",0,"GeV");

  // compare HT for a tight cut for all 3?
  clearSamples();
  addSample("Zinvisible");
  addSample("ZJetsInc");
  addSample("ZJets");
  selection_ = "MET>=120 && njets>=2 && HT>=400";
  nbins = 50; low = 400; high= 1400;
  var = "HT"; xtitle="HT";
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"", "ra2b_DYvZinv_ht",0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"", "ra2b_DYvZinv_ht",0,"GeV");

  //OK, this isn't very interesting
//   setStackMode(false,false);
//   selection_ = "MET>=120 && njets>=2 && HT>=200";
//   nbins = 50; low = 200; high= 1200;
//   var = "HT"; xtitle="HT";
//   setLogY(false);
//   drawPlots(var,nbins,low,high,xtitle,"", "ra2b_DYvZinv_ht_nonorm",0,"GeV");

//use higher stats sample
  setStackMode(false,true);
  clearSamples();
  addSample("Zinvisible");
  addSample("ZJets");
  selection_ = "MET>=150 && njets>=2 && HT>=400";
  nbins = 9; low = 1; high= 10;
  var = "njets"; xtitle="njets (pT > 50 GeV)";
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"", "ra2b_DYvZinv_njets",0,"GeV");

  // also W
  clearSamples();
  addSample("WJetsInc");
  addSample("WJets");
  setStackMode(false,false);
  selection_ = "MET>=150 && njets>=2 && HT>=200";
  nbins = 50; low = 200; high= 1000;
  var = "HT"; xtitle="HT";
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"", "ra2b_WvWinc_ht_nonorm",0,"GeV");


}

void trg() {

  setPlotMinimum(0);

  //main trigger, MET leg
  selection_ = "cutPV==1 && HT>=400 && nElectrons==0 && nMuons==0 && njets>=2 && minDeltaPhiN>=4";
  drawTrigEff("JetHT", "pass_PFHT350==1","pass_PFHT350_PFMET100==1", "MET", 10, 0, 300);
  
  //alphaT, "MET" leg
  selection_ = "cutPV==1&&HT>=400 ";
  drawTrigEff("JetHT", "pass_HT300==1","pass_HT300_AlphaT0p53==1", "MET", 10, 0, 300);

  //HT leg with SingleMu
  selection_ = "cutPV==1  && pass_utilityPrescaleModuleHLT==1";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_PFHT350==1", "HT", 20, 300, 1100);


  //SL eff with MuHad -- beautiful!
  selection_ = "cutPV==1 && HT>=400 && nElectrons==0 && nMuons==1 && muonpt1>=20 && njets>=2";
  drawTrigEff("MuHad", "pass_PFHT350_Mu15_PFMET45==1","pass_PFHT350_PFMET100==1", "MET", 50, 0, 500);

  //PFMET150 with singleMu -- interesting inefficiency! real?
  selection_ = "cutPV==1 && HT>=300";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_PFMET150==1", "MET", 10, 0, 400);

  //DiCentralPFJet50
  selection_ = "cutPV==1 && HT>=400 && jetpt1>=70 && jetpt2>=70";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_DiCentralPFJet50_PFMET80==1", "MET", 20, 0, 500);

  selection_ = "cutPV==1 && HT>=400 && jetpt1>=70 && jetpt2>=70";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_DiCentralPFJet50_PFMET80==1 || pass_PFHT350_PFMET100==1", "MET", 20, 0, 500);

  //HT/MET trigger in single mu
  selection_ = "cutPV==1 && HT>=400 && jetpt1>=70 && jetpt2>=70";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_PFHT350_PFMET100==1", "MET", 20, 0, 500);

  //can we look at the jet leg? -- seems good
  selection_ = "cutPV==1 && HT>=300 && MET>=200 && njets>=2";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_DiCentralPFJet50_PFMET80==1", "jetpt2", 20, 50, 250);

  //also in SingleMu -- alphaT trigger
  selection_ = "cutPV==1 && HT>=400";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_HT250_AlphaT0p55==1", "MET", 10, 0, 400);

  selection_ = "cutPV==1 && HT>=300 && MET>=250";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_HT300_AlphaT0p53==1", "HT", 10, 300, 800);

  //the other alphaT trigger
  //to do


  //DiCentralPFJet30_PFMET80_BTagCSV07

  //add nominal trigger as an OR -- looks pretty good
  selection_ = "cutPV==1 && HT>=300";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_PFMET150==1 || pass_PFHT350_PFMET100==1", "MET", 10, 0, 500);

  //also try OR with caloMET
  selection_ = "cutPV==1 && HT>=300";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_PFMET150==1 || pass_MET120_HBHENoiseCleaned==1 || pass_MET200==1", "MET", 10, 0, 500);

  //caloMET only
  selection_ = "cutPV==1 && HT>=300";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_MET120_HBHENoiseCleaned==1", "MET", 10, 0, 500);
  selection_ = "cutPV==1 && HT>=300";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_MET200==1", "MET", 10, 0, 500);

  //PFMET150 in MuHad
  selection_ = "cutPV==1 && HT>=500 && nElectrons==0 && nMuons==1 && muonpt1>=20 && njets>=2";
  drawTrigEff("MuHad", "pass_PFHT350_Mu15_PFMET45==1","pass_PFMET150==1", "MET", 50, 0, 500);
  //or in JetHT ?????
  selection_ = "cutPV==1 && HT>=450";
  drawTrigEff("JetHT", "pass_PFHT350==1","pass_PFMET150==1", "MET", 10, 0, 400);


  //PFHT350_PFMET100
  //PFMET150

  //PFHT350


  selection_ = "cutPV==1 && MET>=150 && jetpt1>=70 && jetpt2>=70";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_PFHT350_PFMET100==1", "HT", 20, 300, 1100);

  selection_ = "cutPV==1 && MET>=150 && jetpt1>=70 && jetpt2>=70";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_DiCentralPFJet50_PFMET80==1", "HT", 20, 300, 1100);

  selection_ = "cutPV==1 && MET>=150 && jetpt1>=70 && jetpt2>=70";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_DiCentralPFJet50_PFMET80==1||pass_PFHT350_PFMET100==1", "HT", 20, 300, 1100);


  //how about Photon dataset?
  selection_ = "HT>=250 && njets>=2 && jetpt1>=70 && jetpt2>=70 && nMuons==1 &&passCleaning==1";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_DiCentralPFJet50_PFMET80==1", "MET", 20, 0, 500);

  selection_ = "HT>=250 && njets>=2 && jetpt1>=70 && jetpt2>=70 && nMuons==1 &&passCleaning==1";
  drawTrigEff("Photon", "pass_Photon135==1||pass_Photon150==1","pass_DiCentralPFJet50_PFMET80==1", "MET", 20, 0, 500);

  selection_ = "HT>=250 && njets>=2 && jetpt1>=70 && jetpt2>=70 && (nElectrons+nMuons==0) &&passCleaning==1";
  drawTrigEff("Photon", "pass_Photon135==1||pass_Photon150==1","pass_DiCentralPFJet50_PFMET80==1", "MET", 20, 0, 500);

  //now tighten HT cut
  selection_ = "HT>=400 && njets>=2 && jetpt1>=70 && jetpt2>=70 && (nElectrons+nMuons==0) &&passCleaning==1";
  drawTrigEff("Photon", "pass_Photon135==1||pass_Photon150==1","pass_DiCentralPFJet50_PFMET80==1", "MET", 20, 0, 500);

  selection_ = "HT>=400 && njets>=2 && jetpt1>=70 && jetpt2>=70 && (nElectrons+nMuons==0) &&passCleaning==1";
  drawTrigEff("JetHT", "pass_PFHT350==1","pass_DiCentralPFJet50_PFMET80==1", "MET", 20, 0, 500);

  selection_ = "HT>=700 && njets>=2 && jetpt1>=70 && jetpt2>=70 && (nElectrons+nMuons==0) &&passCleaning==1";
  drawTrigEff("JetHT", "pass_PFHT650==1","pass_DiCentralPFJet50_PFMET80==1", "MET", 20, 0, 500);

  //
  //main trigger, MET leg
  selection_ = "HT>=400 && njets>=2 && jetpt1>=70 && jetpt2>=70 && (nElectrons+nMuons==0) &&passCleaning==1";
  drawTrigEff("JetHT", "pass_PFHT350==1","pass_PFHT350_PFMET100==1", "MET", 20, 0, 500);
  

}

void rpv_basic_plots() {

  lumiScale_ = 20000;

  loadSamples();


  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getDefault();

  removeSample("LM9");

  removeSample("VV");

  addSample("sbottom-189-270");

  setSampleScaleFactor("sbottom-189-270",0.04545424); //was 4/9

  dodata_=false;

  setStackMode(false,true);
  //define some cuts
  TCut dilepton = "(nElectrons + nMuons == 2) && (eleet1>=20 || muonpt1>=20)";
  TCut zveto = "bestZmass <81 || bestZmass >101";

  TCut ST = "STeff>=400";
  TCut MET = "MET>=35";

  TCut jets = "njets30>=4";

  TCut rcuts = "rMET <0.15 && rl <0.15";

  TCut resonance = "mjjdiff<10";
  TCut resonance5 = "mjjdiff_5<10";

  removeSample("PythiaPUQCD");
  removeSample("WJets");
  removeSample("Zinvisible");

  //now plots
  int nbins;
  float low,high;
  TString var,xtitle;
  
  //Z mass
  nbins = 50;
  low = 0; high = 300;
  var = "bestZmass"; xtitle="l+l- mass";

  selection_ = dilepton;
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonSelection_llmass",0,"GeV");

  //add n jets
  selection_ = dilepton && jets;
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjetsSelection_llmass",0,"GeV");

  //add z veto; plot ST
  selection_ = dilepton&& jets && zveto;
  nbins = 50;
  low = 300; high = 1800;
  var = "STeff"; xtitle="HT + lepton pT + MET";
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjetsNoZSelection_STeff",0,"GeV");

  //same selection; plot MET
  selection_ = dilepton&& jets && zveto;
  nbins = 50;
  low = 0; high = 400;
  var = "MET"; xtitle="MET";
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjetsNoZSelection_MET",0,"GeV");

  //same selection; plot HT
  selection_ = dilepton&& jets && zveto;
  nbins = 50;
  low = 0; high = 1200;
  var = "HT"; xtitle="HT";
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjetsNoZSelection_HT",0,"GeV");

  //add STeff and MET cuts
  selection_ = dilepton&& jets && zveto &&ST&&MET;
  //plot r quantities
  nbins = 30;
  low = 0; high = 0.5;
  var = "rMET"; xtitle="MET / STeff";
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjetsNoZ-ST-MET-Selection_rMET",0,"GeV");

  nbins = 30;
  low = 0; high = 0.5;
  var = "rl"; xtitle="pT l1 / STeff";
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjetsNoZ-ST-MET-Selection_rl",0,"GeV");

  //add cuts on r quantities
  selection_ = dilepton&& jets && zveto &&ST&&MET&&rcuts;

  nbins = 20;
  low = 0; high = 100;
  var = "mjjdiff"; xtitle="reco stop mass difference";
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjetsNoZ-ST-MET-r-Selection_recoMassDifference",0,"GeV");

  nbins = 20;
  selection_ = dilepton&& jets && zveto &&ST&&MET&&rcuts &&TCut("njets30>=5");
  low = 0; high = 100;
  var = "mjjdiff_5"; xtitle="reco stop mass difference (2+3)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjetsNoZ-ST-MET-r-Selection_recoMassDifference5j",0,"GeV");



  //finally plot the peak
  removeSample("SingleTop");
  removeSample("ZJets");

  // this is the money plot (with no b tag)
  selection_ = dilepton&& jets && zveto &&ST&&MET&&rcuts&&resonance;
  nbins = 40;
  low = 50; high = 450;
  var = "mjj1"; xtitle="mjj";
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjetsNoZ-ST-MET-r-jjdiff-Selection_mjj",0,"GeV");

  setStackMode(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjetsNoZ-ST-MET-r-jjdiff-Selection_mjj",0,"GeV");

  //plot without the 'resonance' cut
  setStackMode(false,true);
  selection_ = dilepton&& jets && zveto &&ST&&MET&&rcuts;
  nbins = 40;
  low = 50; high = 450;
  var = "mjj1"; xtitle="mjj";
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjetsNoZ-ST-MET-r-Selection_mjj",0,"GeV");
  //conclusion -- the Delta mjj cut tightens up the signal peak

  //but what about taking an average>
  setStackMode(false,true);
  selection_ = dilepton&& jets && zveto &&ST&&MET&&rcuts;
  nbins = 40;
  low = 50; high = 450;
  var = "0.5*(mjj1+mjj2)"; xtitle="mjj average";
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjetsNoZ-ST-MET-r-Selection_mjjav",0,"GeV");
  //still not as nice as after the 'resonance

  //try plotting the peak using the 5 jet variable
  setStackMode(false,true);
  selection_ = dilepton&& jets && zveto &&ST&&MET&&rcuts&&resonance5 &&TCut("njets30>=5");
  nbins = 40;
  low = 50; high = 450;
  var = "mjj1_5"; xtitle="mjj (2+3)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjetsNoZ-ST-MET-r-jjdiff-Selection_mjj5j",0,"GeV");
  //result -- not nice at all!

  //only >=5 jets events but using the leading 4
  //also remove the "resonance" cut
  selection_ = dilepton&& jets && zveto &&ST&&MET&&rcuts &&TCut("njets30>=5");
  nbins = 40;
  low = 50; high = 450;
  var = "mjj1"; xtitle="mjj";
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjets5NoZ-ST-MET-r-jjdiff-Selection_mjj",0,"GeV");
  //not as bad as above...

  setStackMode(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjets5NoZ-ST-MET-r-jjdiff-Selection_mjj",0,"GeV");
  //Conclusion -- don't see the benefit of this 2+3 variable!


  // try add b-jet veto
  setStackMode(false,true);
  selection_ = dilepton&& jets && zveto &&ST&&MET&&rcuts&&resonance && TCut("nbjets == 0");
  nbins = 40;
  low = 50; high = 450;
  var = "mjj1"; xtitle="mjj";
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjetsNoZ-ST-MET-r-jjdiff-bveto-Selection_mjj",0,"GeV");

  setStackMode(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjetsNoZ-ST-MET-r-jjdiff-bveto-Selection_mjj",0,"GeV");

  //now plot signal only
  clearSamples();
  dodata_=false;
  //ideally i want a way to plot with a sample split into mutliple pieces, treated as different subsamples by the code
  addSample("sbottom-189-270:nCorrectRecoStop==0");
  addSample("sbottom-189-270:nCorrectRecoStop==1");
  addSample("sbottom-189-270:nCorrectRecoStop==2");
  setSampleColor("sbottom-189-270:nCorrectRecoStop==0",1);
  setSampleColor("sbottom-189-270:nCorrectRecoStop==1",2);
  setSampleColor("sbottom-189-270:nCorrectRecoStop==2",4);
  setStackMode(false,true);
  selection_ = dilepton&& jets && zveto &&ST&&MET&&rcuts&&resonance && TCut("nbjets == 0");
  nbins = 40;
  low = 50; high = 450;
  var = "mjj1"; xtitle="mjj";
  drawPlots(var,nbins,low,high,xtitle,"Events", "rpv_dileptonjetsNoZ-ST-MET-r-jjdiff-bveto-Selection_mjj_signalSplit",0,"GeV");


}

void profmattstrassler() {


    lumiScale_ = 5000;

    loadSamples();


    useHTeff_=false;
    useMHTeff_=false;
    thebnnMHTeffMode_ = kOff;
    currentConfig_=configDescriptions_.getDefault();

    removeSample("LM9");

    //    removeSample("VV");
    addSample("sbottom-189-270");
    setSampleScaleFactor("sbottom-189-270",0.04545424); //was 4/9                                                                                                   
    dodata_=false;

    TCut lepton = "(nElectrons15 >=1 || nMuons15>=1) && MT_Wlep15>20";
    TCut jets = "njets30>=6";
    //now plots                                                                                                                                                     
    int nbins;
    float low,high;
    TString var,xtitle;

    nbins = 40;
    low = 300; high = 1500;
    var = "STeff"; xtitle="ST";

    selection_ = lepton&&jets;
    drawPlots(var,nbins,low,high,xtitle,"Events", "pms_lepton6jets30_ST",0,"GeV");


}


void T1shape() {

  loadSamples();

  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getCorrected();
  dodata_=false;

  setStackMode(false,true);

  clearSamples();
  addSample("T1bbbb$800$100");
  addSample("T1bbbb$800$400");
  addSample("T1bbbb$400$100");
  useMassInLegend_=true;

   setSampleColor("T1bbbb$800$100",kRed+1);
   setSampleColor("T1bbbb$800$500",kBlue);
   setSampleColor("T1bbbb$400$100",kGreen+2);

   leg_x1=0.5;
  //   setSampleColor("T1bbbb$925$100",kRed+1);
  
  //   setSampleLineStyle("T1bbbb$925$100",1);
  
  int nbins;
  float low,high;
  TString var,xtitle;
  
  //MET
  nbins = 50;
  low = 0; high = 800;
  var = "MET"; xtitle="MET";
  selection_="cutPV==1 && HT>=400 && njets>=2";
  drawPlots(var,nbins,low,high,xtitle,"Events", "SMSexamples_T1bbbb_MET",0,"GeV");

  //b jet pT
  nbins = 50;
  low = 0; high = 400;
  var = "bjetpt3"; xtitle="3rd leading b jet pT";
  selection_="cutPV==1 && HT>=400 && njets>=2 && nbjets>=3";
  drawPlots(var,nbins,low,high,xtitle,"Events", "SMSexamples_T1bbbb_bjetpt",0,"GeV");

}


void stop_misc() {


  //skim applied was ST>=150 && MET>=40
  inputPath = "/cu3/joshmt/reducedTrees/V00-02-35aa/";
  loadSamples();
  clearSamples();
  addSample("TTbarJets");
  addSample("T2tt$250$50",2);
  addSample("T2tt$350$50",3);
  addSample("T2tt$450$50",6);
  addSample("T2tt$900$50",kOrange);

  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getDefault();

  dodata_=false;

  setStackMode(false,true);
 
  TCut basic = "njets30>=6 && MET>=100 && nbjets>=1 && nElectrons5==0 && nMuons5==0";


  int nbins;
  float low,high;
  TString var,xtitle;
  nbins = 40;
  low = 0; high = 1200;
  var = "mjjb1"; xtitle="m_{jjb} [higher]";

  selection_ = basic;
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_6jetsEtc_mjjb1",0);

  nbins = 40;
  low = 0; high = 400;
  var = "mjjb2"; xtitle="m_{jjb} [lower]";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_6jetsEtc_mjjb2",0);


  nbins = 20;
  low = 0; high = 800;
  var = "mjjb1-mjjb2"; xtitle="m_{jjb} difference";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_6jetsEtc_mjjbdiff",0);

  nbins = 20;
  low = 0; high = 600;
  var = "topPT1"; xtitle="top pt 1";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_6jetsEtc_topPt1",0);
  var = "topPT2"; xtitle="top pt 2";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_6jetsEtc_topPt2",0);

  doOverflowAddition(false);
  nbins = 20;
  low = -500; high = 500;
  var = "topPT1-topPT2"; xtitle="top pt diff";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_6jetsEtc_topPtDiff",0);


  doOverflowAddition(true);
  nbins = 30;
  low = 0; high = 800;
  var = "abs(topPT1-topPT2)"; xtitle="abs top pt diff";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_6jetsEtc_topPtDiffAbs",0);

  basic = "njets30>=6 && MET>=175 && nbjets>=1 && nElectrons5==0 && nMuons5==0";
  selection_ = basic;
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_6jetsEtcMET175_topPtDiffAbs",0);


}
