/*
to compile:

must have symlink to MiscUtil.cxx in the working directory.

gSystem->Load("TSelectorMultiDraw_C.so");
gSystem->Load("CrossSectionTable_cxx.so");
gSystem->Load("ConfigurationDescriptions_cxx.so");
gSystem->Load("SearchRegion_cxx.so");
gSystem->Load("SystInfo_cxx.so");
gSystem->Load("SignalEffData_cxx.so");
.L drawReducedTrees2012.C++

====== this is the nominal code for drawing RA2b data/MC comparison plots =======
The signal efficiency functions are now forked into signalEffSyst.C

The functions that do the heavy lifting, along with various utility functions,
are in drawReducedTrees.h.

This code replaces drawCutflowPlots.C (which had replaced drawBasicPlots.C). 
The input files are the reducedTrees.

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
 
  // == test of the HTMHT + MET business ==
  setDatasetToDraw("2012hybrid");
  lumiScale_ = 700;

  selection_ =TCut("HT>=200 && MET>=175 && cutPV==1 && cutTrigger==1 && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0")&&getSingleLeptonCut();
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 22; low=200; high=1300;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_eq0b_SL_hybrid",0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_eq0b_SL_hybrid",0,"GeV");


}

void RA2b2012_mtspectrum() {


  nameOfEventWeight_="weight2"; //for 2012 cfA skimmed ntuples

  loadSamples();
  usePUweight_=false;

  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getDefault();

  clearSamples();
  addSample("TTbarJets:nbjets==1",kRed,"==1 b");
  addSample("TTbarJets:nbjets==2",kGreen,"==2 b");
  addSample("TTbarJets:nbjets>=3",kBlue,">=3 b");

  int nbins;
  float low,high;
  TString var,xtitle;
  
  doOverflowAddition(true);

  dodata_=false;
  setStackMode(false,true);//normalized

  selection_ ="HT>=400 && MET>=150 && cutPV==1 && cutTrigger==1 && njets>=3 && minDeltaPhiN >= 4 && passCleaning==1&& (nElectrons+nMuons==1)";

  var="MT_Wlep"; xtitle="m_{T} [GeV]";
  nbins = 10; low=0; high=200; 
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_mt",0,"GeV");


}

void RA2b2012_ttbarw() {

  nameOfEventWeight_="weight2"; //for 2012 cfA skimmed ntuples

  loadSamples();
  usePUweight_=true;

  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getDefault();

  clearSamples();
  addSample("TTbarJets");
  addSample("WJets");
  addSample("SingleTop");


  int nbins;
  float low,high;
  TString var,xtitle;
  
  doOverflowAddition(false);
  drawMCErrors_=true;

  dodata_=false;


  setStackMode(true);//regular stack
  
  const int nvarbins=11;
  const float varbins[]={100, 125, 150, 175, 200, 225, 250,  300, 350, 400,  500,  650}; //AN and PAS 
  
  selection_ =TCut("HT>=400 && rawPFMET>=150 && cutPV==1 && cutTrigger==1  && njets>=3 && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1")&&getSingleLeptonCut();
  var="rawPFMET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 55; low=100; high=650; 
  drawPlots(var,nvarbins, varbins, xtitle,"Arbitrary units", "ra2b_rawPFMET_SL_ttwt_ge1b_ht400");
  
  TH1D* SLplot = (TH1D*) totalsm->Clone("SLplot");

  SLplot->Scale(1/SLplot->Integral());

  selection_ =TCut("HT>=400 && rawPFMET>=150 && cutPV==1 && cutTrigger==1  && njets>=3 && minDeltaPhiN >= 4  &&passCleaning==1 && nbjets>=1")&&getLeptonVetoCut();
  var="rawPFMET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 55; low=100; high=650; 
  drawPlots(var,nvarbins, varbins, xtitle,"Arbitrary units", "ra2b_rawPFMET_0L_ttwt_ge1b_ht400");
  
  TH1D* ZLplot = (TH1D*) totalsm->Clone("ZLplot");
  ZLplot->Scale(1/ZLplot->Integral());

  renewCanvas("ratio");

  thecanvas->cd(1);
  ZLplot->SetLineColor(kBlue);
  ZLplot->SetMarkerColor(kBlue);
  ZLplot->Draw();
 
  SLplot->SetLineColor(kRed);
  SLplot->SetMarkerColor(kRed);
  SLplot->Draw("SAME");
   
  thecanvas->cd(2);
  ratio = (TH1D*) ZLplot->Clone("ratio");
  ratio->SetLineColor(kBlack);
  ratio->SetMarkerColor(kBlack);
  ratio->Divide(ZLplot,SLplot);
  ratio->Draw();

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


