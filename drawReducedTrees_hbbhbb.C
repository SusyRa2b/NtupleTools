/*
to compile:

must have symlink to MiscUtil.cxx in the working directory.

gSystem->Load("TSelectorMultiDraw_C.so");
gSystem->Load("CrossSectionTable_cxx.so");
gSystem->Load("ConfigurationDescriptions_cxx.so");
.L drawReducedTrees_hbbhbb.C+

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
#include "TH3.h"
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
	
TString inputPath = "/cu4/ra2b/reducedTrees/v66_10/"; //2012
TString dataInputPath =  "/cu4/ra2b/reducedTrees/v66_10/ORIGINALS/";
//TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-35u/Fall11/";//7 TeV

double lumiScale_ = 19399; //Run 2012 ABC+D (update from Keith)

//make a symlink that point from this name to drawReducedTree.h
//this is to make the ROOT dictionary generation work correctly
#include "drawReducedTrees_hbbhbb.h"


void initHiggsSamples(bool useSkim=true,const TString samplelist="") {

  const  TString skimstring = useSkim?"-skim":"";

  addToSamplesAll(TString("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66")+skimstring);
  addToSamplesAll(TString("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1613_v66")+skimstring);
  addToSamplesAll(TString("TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1606_v66")+skimstring);

  addToSamplesAll(TString("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66")+skimstring);
  addToSamplesAll(TString("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1607_v66")+skimstring);
  addToSamplesAll(TString("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1524_v66")+skimstring);
  addToSamplesAll(TString("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1594_v66")+skimstring);
  addToSamplesAll(TString("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1523_v66")+skimstring);
  addToSamplesAll(TString("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1602_v66")+skimstring);

  addToSamplesAll(TString("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66")+skimstring);
  addToSamplesAll(TString("QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1513_v66")+skimstring);
  addToSamplesAll(TString("QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1578_v66")+skimstring);
  addToSamplesAll(TString("QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1561_v66")+skimstring);
  addToSamplesAll(TString("QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1603_v66")+skimstring);
  addToSamplesAll(TString("QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1585_v66")+skimstring);
  addToSamplesAll(TString("QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1560_v66")+skimstring);
  addToSamplesAll(TString("QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_v3_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1609_v66")+skimstring);
  addToSamplesAll(TString("QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1515_v66")+skimstring);
  addToSamplesAll(TString("QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1516_v66")+skimstring);
  addToSamplesAll(TString("QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1559_v66")+skimstring);

  addToSamplesAll("TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1605_v66");
  addToSamplesAll("TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1604_v66");

  addToSamplesAll("TChihh_150_v68");
  addToSamplesAll("TChihh_200_v68");
  addToSamplesAll("TChihh_250_v68");
  addToSamplesAll("TChihh_400_v68");

  addToSamplesAll(TString("WbbJetsToLNu_Massive_TuneZ2star_8TeV-madgraph-pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1677_v67")+skimstring);

  nameOfEventWeight_="weight3";

  inputPath="/cu6/joshmt/reducedTrees/v68_4k/";
  //careful -- cannot mix skimmed and unskimmed data files in one directory because of the way wildcards are used
  dataInputPath="/cu6/joshmt/reducedTrees/v68_5/"; //for skimmed data
  if (!useSkim)   dataInputPath="/cu6/joshmt/reducedTrees/v68_5/unskimmedData/";

  loadSamples(true,"ra2b2012");
  usePUweight_=true; //helps bring signal and background into alignment; not perfect but better than nothing

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;
  clearSamples();

  if (samplelist=="" || samplelist.Contains("qcd")) {
    addSample(TString("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66")+skimstring,kGreen,"QCD");
    chainSamples(TString("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66")+skimstring,TString("QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1513_v66")+skimstring);
    chainSamples(TString("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66")+skimstring,TString("QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1578_v66")+skimstring);
    chainSamples(TString("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66")+skimstring,TString("QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1561_v66")+skimstring);
    chainSamples(TString("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66")+skimstring,TString("QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1603_v66")+skimstring);
    chainSamples(TString("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66")+skimstring,TString("QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1585_v66")+skimstring);
    chainSamples(TString("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66")+skimstring,TString("QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1560_v66")+skimstring);
    chainSamples(TString("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66")+skimstring,TString("QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_v3_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1609_v66")+skimstring);
    chainSamples(TString("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66")+skimstring,TString("QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1515_v66")+skimstring);
    chainSamples(TString("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66")+skimstring,TString("QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1516_v66")+skimstring);
    chainSamples(TString("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66")+skimstring,TString("QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1559_v66")+skimstring);
    //include had ttbar in QCD
    chainSamples(TString("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66")+skimstring,TString("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1613_v66")+skimstring);
  }

  if (samplelist=="" || samplelist.Contains("znunu")) {
    addSample(TString("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66")+skimstring,kOrange,"Z #rightarrow #nu#nu");
    chainSamples(TString("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66")+skimstring,TString("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1607_v66")+skimstring);
    chainSamples(TString("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66")+skimstring,TString("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1524_v66")+skimstring);
    chainSamples(TString("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66")+skimstring,TString("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1594_v66")+skimstring);
    chainSamples(TString("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66")+skimstring,TString("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1523_v66")+skimstring);
    chainSamples(TString("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66")+skimstring,TString("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1602_v66")+skimstring);
  }

  if (samplelist=="" || samplelist.Contains("ttv")) {
    addSample("TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1605_v66",kViolet,"t#bar{t}V");
    chainSamples("TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1605_v66","TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1604_v66");
  }

  if (samplelist=="" || samplelist.Contains("wbb")) {
    addSample(TString("WbbJetsToLNu_Massive_TuneZ2star_8TeV-madgraph-pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1677_v67")+skimstring,kMagenta,"Wbb");
  }

  if (samplelist=="ttbarjoint") {
    addSample(TString("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66")+skimstring,kBlue-10,"t#bar{t} (#geq 1l)");
    chainSamples(TString("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66") +skimstring ,TString("TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1606_v66") +skimstring );

  setSampleWeightFactor(TString("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66")+skimstring,"topPtWeight");
  }
  else if (samplelist=="" || samplelist.Contains("ttbar")) {
    addSample(TString("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66")+skimstring,kBlue-10,"t#bar{t} (2l)");
    addSample(TString("TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1606_v66") +skimstring,kAzure-3,"t#bar{t} (1l)");

    setSampleWeightFactor(TString("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66")+skimstring,"topPtWeight");
    setSampleWeightFactor(TString("TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1606_v66") +skimstring,"topPtWeight");
  }



  //  if (samplelist=="" || samplelist.Contains("hh150")) addSample("TChihh_150_v68",kRed-10,"hh 150"); //remove 150 point
  if (samplelist=="" || samplelist.Contains("hh200")) addSample("TChihh_200_v68",kRed-7,"hh 200");
  if (samplelist=="" || samplelist.Contains("hh250")) addSample("TChihh_250_v68",kRed,"hh 250");
  if (samplelist=="" || samplelist.Contains("hh400")) addSample("TChihh_400_v68",kRed+2,"hh 400");

  //official source, I think https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections8TeVhiggsinoCMS
  const float hbb = 0.577;  //https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR2 -- i'm assuming m=125
  //  if (samplelist=="" || samplelist.Contains("hh150"))setSampleScaleFactor("TChihh_150_v68",2.141*hbb*hbb/9999.0); //sigma x BF / ngen //remove 150 point
  if (samplelist=="" || samplelist.Contains("hh200"))setSampleScaleFactor("TChihh_200_v68",0.6975*hbb*hbb/9999.0); //sigma x BF / ngen
  if (samplelist=="" || samplelist.Contains("hh250"))setSampleScaleFactor("TChihh_250_v68",0.271*hbb*hbb/9999.0); //sigma x BF / ngen
  if (samplelist=="" || samplelist.Contains("hh400"))setSampleScaleFactor("TChihh_400_v68",0.03*hbb*hbb/9999.0); //sigma x BF / ngen

  setDatasetToDraw("MET"); //for now we're using the MET dataset only

}



void higgs_genlevel_signal() {


  initHiggsSamples(false,"hh200 hh250 hh400");


  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";
  //cannot use this with MET PD only
  //  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;

  TCut metsigveryloose="METsig>30";

  TCut drmax = "deltaRmax_hh<2.4";
  TCut drmin = "deltaRmin_hh<1.9";

  //  setSampleScaleFactor("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim",0.8);

  setStackMode(false,true,false); //stack,norm,label override
  setLogY(false);
  usePUweight_=false; //helps bring signal and background into alignment; not perfect but better than nothing

  selection_="";
  //start with no selection

  //pt
  nbins=20; low=0; high=200;
  setLogY(false);
  var="higgsb1pt"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_genlevel_"+var,0);

  var="higgsb2pt"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_genlevel_"+var,0);

  var="higgsb3pt"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_genlevel_"+var,0);

  var="higgsb4pt"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_genlevel_"+var,0);

  //eta
  nbins=20; low=-5; high=5;
  setLogY(false);
  var="higgs1b1eta"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_genlevel_"+var,0);

  var="higgs1b2eta"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_genlevel_"+var,0);

  var="higgs2b1eta"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_genlevel_"+var,0);

  var="higgs2b2eta"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_genlevel_"+var,0);




  clearSamples();
  addSample(TString("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66"),kBlue-10,"t#bar{t} (#geq 1l)");
  chainSamples(TString("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66")  ,TString("TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1606_v66") );

  var="jetpt4"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy"+var,0);

  var="jetpt2"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy"+var,0);

}

void higgs_signal_scult() {

  initHiggsSamples(false,"hh250 hh400");


  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";
  //cannot use this with MET PD only
  //  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;

  TCut metsigveryloose="METsig>30";

  TCut drmax = "deltaRmax_hh<2.4";
  TCut drmin = "deltaRmin_hh<1.9";

  //  setSampleScaleFactor("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim",0.8);

  setStackMode(false,true,false); //stack,norm,label override
  setLogY(false);

  clearSamples();
  addSample("TChihh_250_v68:(1)",kBlue,"hh250");
  addSample("TChihh_250_v68:deltaRmax_hh<2.4",kRed,"hh250 (#Delta R_{max} < 2.4)");

  selection_ = baseline && trigger && zl&&isotk && (njets4||njets5) && btag2&&btag3&&btag4 && higgsSR_d && metsigveryloose;

  //higgs mass
  nbins=20; low=0; high=200;
  setLogY(false);
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{bb}>";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_sculpt_signal",0);


}

void higgs_ttbar_sculpt() {

  initHiggsSamples(true,"ttbarjoint");


  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";
  //cannot use this with MET PD only
  //  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;

  TCut metsigveryloose="METsig>30";

  TCut drmax = "deltaRmax_hh<2.4";
  TCut drmin = "deltaRmin_hh<1.9";

  //  setSampleScaleFactor("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim",0.8);

  setStackMode(false,true,false); //stack,norm,label override
  setLogY(false);

  clearSamples();
  addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim:(1)",kBlue,"t#bar{t}");
  addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim:deltaRmax_hh<2.4",kRed,"t#bar{t} (#Delta R_{max} < 2.4)");

   //  addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim:deltaRmax_hh<2.2",kRed-7,"t#bar{t} (#Delta R_{max} < 2.2)");
  //  addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim:deltaRmax_hh<2.2&&deltaRmin_hh<1.9",kMagenta,"t#bar{t} (#Delta R_{max} < 2.4, #Delta R_{min} < 1.9)");

  selection_ = baseline && trigger && zl&&isotk && (njets4||njets5) && btag2&&btag3&&btag4 && higgsSR_d && metsigveryloose;

  //higgs mass
  nbins=20; low=0; high=200;
  setLogY(false);
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{bb}>";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_sculpt_ttbar",0);


  clearSamples();
  addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim:CSVbest3<0.244",kBlue,"t#bar{t} (2b)");
  addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim:CSVbest3>0.679&&CSVbest4<0.244",kRed,"t#bar{t} (3b)");
  addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim:CSVbest3>0.679&&CSVbest4>0.244",kBlack,"t#bar{t} (4b)");

  selection_ = baseline && trigger && zl&&isotk && (njets4||njets5) && btag2 && higgsSR_d && TCut("METsig>50")&&drmax;

  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_sculpt_ttbar_btag",0);

  //higgs mass
  nbins=30; low=0; high=300;
  setLogY(false);
  var="jetpt1"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_sculpt_ttbar_btag_"+var,0);

  var="jetpt2"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_sculpt_ttbar_btag_"+var,0);

  nbins=20; low=0; high=200;
  var="jetpt3"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_sculpt_ttbar_btag_"+var,0);

  nbins=10; low=0; high=100;
  var="jetpt4"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_sculpt_ttbar_btag_"+var,0);


}

void higgs_dataMC_control() {

  //goal do data / MC comparisons; use skim to save time. means I can't make a useful njets plot
  initHiggsSamples(true,"ttbar wbb ttv znunu qcd"); //use all samples except signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";
  //cannot use this with MET PD only
  //  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;

  TCut metsigveryloose="METsig>30";

  TCut drmax = "deltaRmax_hh<2.4";
  //  TCut drmin = "deltaRmin_hh<1.9";

  //  setSampleScaleFactor("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim",0.8);

  setStackMode(true,false,false); //stack,norm,label override
  setLogY(false);

  //2b pre-selection ; basically the minimum required by the trigger, plus SL
  selection_ = baseline && trigger && sl && (njets4||njets5) && btag2;

  //Owen's METsig bins
  //float binning[]={0,10,20,30,50,100,150,200};
  nbins=30; low=0; high=300;
  setLogY(false);
  var="METsig"; xtitle="MET significance";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig",0);

  //MET, for completeness
  nbins=30; low=0; high=400;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_MET",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_MET",0);

  //CSV3
  nbins=20; low=0; high=1;
  setLogY(false);
  var="CSVbest3"; xtitle="3rd CSV value";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_CSVbest3",0);

  //CSV4
  nbins=20; low=0; high=1;
  setLogY(false);
  var="CSVbest4"; xtitle="4th CSV value";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_CSVbest4",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_CSVbest4",0);

  //apply 3rd b tag
  selection_ = baseline && trigger && sl && (njets4||njets5) && btag2 &&btag3;

  //CSV4
  nbins=10; low=0; high=1;
  setLogY(false);
  var="CSVbest4"; xtitle="4th CSV value";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_3bPreselection_CSVbest4",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_3bPreselection_CSVbest4",0);

  //now apply full b-tag selection
  selection_ = baseline && trigger && sl && (njets4||njets5) && btag2 &&btag3&&btag4;

  //MET sig
  nbins=15; low=0; high=300;
  setLogY(false);
  var="METsig"; xtitle="MET significance";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bPreselection_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bPreselection_METsig",0);
  //MET
  nbins=15; low=0; high=400;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bPreselection_MET",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bPreselection_MET",0);

 //now apply full b-tag selection and METpreselection
  selection_ = baseline && trigger && sl && (njets4||njets5) && btag2 &&btag3&&btag4 &&metsigveryloose;
  nbins=10; low=0; high=100;
  setLogY(false);
  var="abs(higgsMbb1MassDiff-higgsMbb2MassDiff)"; xtitle="|#Delta m_{bb}|";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bMETsig30_DeltaMbb",0);

  //higgs mass
  nbins=20; low=0; high=200;
  setLogY(false);
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{bb}>";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bMETsig30_AvMbb",0);


  //where is this stupid QCD event?
  nbins=10; low=0; high=4;
  setLogY(false);
  var="deltaPhi2"; xtitle="#Delta #phi (jet2,MET)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bMETsig30_DeltaPhi2",0);

  //now add the higgs SR
  selection_ = baseline && trigger && sl && (njets4||njets5) && btag2 &&btag3&&btag4 &&metsigveryloose &&higgsSR;
  nbins=20; low=0; high=5;
  setLogY(false);
  var="deltaRmax_hh"; xtitle="#Delta R_{max}(b,b)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_hhSR-METsig30_DRmax",0);

  nbins=20; low=0; high=5;
  setLogY(false);
  var="deltaRmin_hh"; xtitle="#Delta R_{min}(b,b)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_hhSR-METsig30_DRmin",0);

  //true N - 1 plots
  selection_ = baseline && trigger && sl && (njets4||njets5) && btag2 &&btag3&&btag4 &&metsigveryloose &&higgsSR_d &&drmax;
  //higgs mass
  nbins=20; low=0; high=200;
  setLogY(false);
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{bb}>";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_Nm1METsig30_AvMbb",0);

  selection_ = baseline && trigger && sl && (njets4||njets5) && btag2 &&btag3&&btag4  &&higgsSR &&drmax;
  nbins=30; low=0; high=300;
  setLogY(true);
  var="METsig"; xtitle="MET significance";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_Nm1_METsig",0);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_Nm1_METsig",0);



}

void higgs_SLratios2() {
  //owen asked me to calculate  [ N(4bSL,data) / N(4bSL,MC) ]   /   [ N(2bSL,data) / N(2bSL,MC)]
  initHiggsSamples(true,"znunu ttbar wbb ttv");

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";
  //  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";
  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";
  TCut notbtag3="CSVbest3<0.244";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut metsigloose="METsig>30";

  TCut drmax = "deltaRmax_hh<2.4";
  //  TCut drmin = "deltaRmin_hh<1.9";


  //Owen's METsig bins
  float binning2[]={0,10,20,30,50,100,200};
  float * binning=0;
  binning=binning2;
  nbins=7; low=0; high=200;
  nbins=6; //use binning2
  var="METsig"; xtitle=var;

  savePlots_=false;

  setStackMode(true,false,false); //stack,norm,labels
  setLogY(false);

  selection_=baseline&&trigger&&sl&&(njets4||njets5)&& btag2&&btag3&&btag4 &&higgsSR&&drmax;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_shapes",binning);
  TH1D* ratio_4b = (TH1D*) ratio->Clone("ratio_4b");  

  selection_=baseline&&trigger&&sl&&(njets4||njets5)&& btag2&&notbtag3 &&higgsSR&&drmax;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_shapes",binning);
  TH1D* ratio_2b = (TH1D*) ratio->Clone("ratio_2b");  

  TH1D* doubleratio=(TH1D*) hdata->Clone("doubleratio");
  doubleratio->Reset();
  doubleratio->Divide(ratio_4b,ratio_2b);
  renewCanvas();
  doubleratio->Draw();

  //looks ok, very limited stats

}

void higgs_ttbarjetmult_sbsig() {
  // goals:
  // reproduce owen's SR/SB numbers in bins of b multiplicity and METsig
  // see if there is a dependence on njets
  // can also compare SR and SB for things like jet pT spectra (which would correlate with b-tag efficiency)
  initHiggsSamples(true,"znunu ttbar wbb ttv");
  //  initHiggsSamples(true,"ttbarjoint");

  const  TString sample = "sl"; //or zl
  assert(sample=="sl" || sample=="zl");

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";
  //  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";
  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";
  TCut notbtag3="CSVbest3<0.244";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;

  TCut metsigloose="METsig>50";

  TCut drmax = "deltaRmax_hh<2.4";
  //  TCut drmin = "deltaRmin_hh<1.9";


  //Owen's METsig bins
  float binning1[]={0,10,20,30,50,100,150,200};
  float binning2[]={0,10,20,30,50,100,200};
  float * binning=0;
  if (sample=="zl") binning = binning1;
  else if (sample=="sl") binning=binning2;
  nbins=7; low=0; high=200;
  if (sample=="sl") nbins=6; //use binning2
  var="METsig"; xtitle=var;

  //samples have already been chained. trick is to not reset the chains, then clear the samples and re-add them
  /* nevermind
     clearSamples();
     addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66:CSVbest2>0.898",kBlue,"t#bar{t} (2b)");
     addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66:CSVbest2>0.898&&CSVbest3>0.679",kRed,"t#bar{t} (3b)");
     addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66:CSVbest2>0.898&&CSVbest3>0.679&&CSVbest4>0.244",kBlack,"t#bar{t} (3b)");
  */

  savePlots_=false;

  setStackMode(false,false,false); //no stack
  setLogY(false);

  //try in bins of njets
  TCut njetsv[3];  TString njetsdesc[3];
  njetsv[0]=njets4||njets5; njetsdesc[0]="njets45";
  njetsv[1]=njets4; njetsdesc[1]="njets4";
  njetsv[2]=njets5; njetsdesc[2]="njets5";

  TCut leptoncut="";
  if (sample=="zl")     leptoncut = zl&&isotk; //lepton veto for signal region
  else if (sample=="sl") {    leptoncut = sl; dodata_=true;} //for data studies use SL

  //    for (int ijets=0;ijets<1 /*3*/;ijets++) { //3 is for jet multiplicity studies
  const int ijets=0;
    TCut njets=njetsv[ijets];

    //SR cuts ; 2 btags
    selection_ = baseline&&trigger&&leptoncut&&njets&&drmax && higgsSR && btag2 && notbtag3; // was !btag3
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SRSBshapes",binning);
    TH1D* h_SR_2b = (TH1D*) totalsm->Clone("h_SR_2b");
    TH1D* h_SR_2b_data=0; if (dodata_) h_SR_2b_data = (TH1D*) hdata->Clone("h_SR_2b_data");
    
    //SR cuts ; 3 btags
    selection_ = baseline&&trigger&&leptoncut&&njets&&drmax && higgsSR && btag2&&btag3 &&!btag4;
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SRSBshapes",binning);
    TH1D* h_SR_3b = (TH1D*) totalsm->Clone("h_SR_3b");
    TH1D* h_SR_3b_data=0; if (dodata_) h_SR_3b_data = (TH1D*) hdata->Clone("h_SR_3b_data");
    
    //SR cuts ; 4 btags
    selection_ = baseline&&trigger&&leptoncut&&njets&&drmax && higgsSR && btag2&&btag3&&btag4;
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SRSBshapes",binning);
    TH1D* h_SR_4b = (TH1D*) totalsm->Clone("h_SR_4b");
    TH1D* h_SR_4b_data=0; if (dodata_) h_SR_4b_data = (TH1D*) hdata->Clone("h_SR_4b_data");

    // -- now SB --
    
    //SB cuts ; 2 btags
    selection_ = baseline&&trigger&&leptoncut&&njets&&drmax && higgsSB && btag2 &&notbtag3;
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SRSBshapes",binning);
    TH1D* h_SB_2b = (TH1D*) totalsm->Clone("h_SB_2b");
    TH1D* h_SB_2b_data=0; if (dodata_) h_SB_2b_data = (TH1D*) hdata->Clone("h_SB_2b_data");
    
    //SB cuts ; 3 btags
    selection_ = baseline&&trigger&&leptoncut&&njets&&drmax && higgsSB && btag2&&btag3 &&!btag4;
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SRSBshapes",binning);
    TH1D* h_SB_3b = (TH1D*) totalsm->Clone("h_SB_3b");
    TH1D* h_SB_3b_data=0; if (dodata_) h_SB_3b_data = (TH1D*) hdata->Clone("h_SB_3b_data");

    //SB cuts ; 4 btags
    selection_ = baseline&&trigger&&leptoncut&&njets&&drmax && higgsSB && btag2&&btag3&&btag4;
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SRSBshapes",binning);
    TH1D* h_SB_4b = (TH1D*) totalsm->Clone("h_SB_4b");
    TH1D* h_SB_4b_data=0; if (dodata_) h_SB_4b_data = (TH1D*) hdata->Clone("h_SB_4b_data");

    TH1D* ratio_2b = (TH1D*)  h_SR_2b->Clone("ratio_2b");
    TH1D* ratio_3b = (TH1D*)  h_SR_2b->Clone("ratio_3b");
    TH1D* ratio_4b = (TH1D*)  h_SR_2b->Clone("ratio_4b");

    ratio_2b->Reset();
    ratio_3b->Reset();
    ratio_4b->Reset();

    ratio_2b->Divide(h_SR_2b,h_SB_2b);
    ratio_3b->Divide(h_SR_3b,h_SB_3b);
    ratio_4b->Divide(h_SR_4b,h_SB_4b);
    
    ratio_2b->SetLineColor(kBlue);
    ratio_3b->SetLineColor(kRed);
    ratio_4b->SetLineColor(kMagenta);
    
    ratio_2b->SetMarkerColor(kBlue);
    ratio_3b->SetMarkerColor(kRed);
    ratio_4b->SetMarkerColor(kMagenta);
    
    ratio_2b->SetMarkerSize(1);
    ratio_3b->SetMarkerSize(1);
    ratio_4b->SetMarkerSize(1);

    ratio_2b->SetMarkerStyle(25);
    ratio_3b->SetMarkerStyle(25);
    ratio_4b->SetMarkerStyle(25);

    
    TH1D* ratio_2b_data=0;
    TH1D* ratio_3b_data=0;
    TH1D* ratio_4b_data=0;
    if (dodata_) {
      ratio_2b_data = (TH1D*)  h_SR_2b->Clone("ratio_2b_data");
      ratio_3b_data = (TH1D*)  h_SR_2b->Clone("ratio_3b_data");
      ratio_4b_data = (TH1D*)  h_SR_2b->Clone("ratio_4b_data");
      ratio_2b_data->Reset();
      ratio_3b_data->Reset();
      ratio_4b_data->Reset();
      ratio_2b_data->Divide(h_SR_2b_data,h_SB_2b_data);
      ratio_3b_data->Divide(h_SR_3b_data,h_SB_3b_data);
      ratio_4b_data->Divide(h_SR_4b_data,h_SB_4b_data);
      
      ratio_2b_data->SetLineColor(kBlue-9);
      ratio_3b_data->SetLineColor(kRed-9);
      ratio_4b_data->SetLineColor(kMagenta-9);
      
      ratio_2b_data->SetMarkerColor(kBlue-9);
      ratio_3b_data->SetMarkerColor(kRed-9);
      ratio_4b_data->SetMarkerColor(kMagenta-9);
      
      ratio_2b_data->SetMarkerSize(1);
      ratio_3b_data->SetMarkerSize(1);
      ratio_4b_data->SetMarkerSize(1);

      ratio_2b_data->SetMarkerStyle(31);
      ratio_3b_data->SetMarkerStyle(31);
      ratio_4b_data->SetMarkerStyle(31);

    }

    //need to split onto multiple canvases so that I can see the various curves and their agreement.
    //or improve the legibility somehow

    TString drawopt="hist e";

    if (dodata_)   setPadDimensions(1000,600);
    renewCanvas("",dodata_? 3:1);
    if (dodata_) thecanvas->cd(3);
    ratio_4b->Draw(drawopt);
    ratio_4b->SetMinimum(0);
    ratio_4b->SetMaximum(0.35);// 0.65

    if (!dodata_) drawopt+=" same";
    if (dodata_) thecanvas->cd(2);
    ratio_3b->Draw(drawopt);
    if (dodata_) thecanvas->cd(1);
    ratio_2b->Draw(drawopt);
    
    if (dodata_) {
      ratio_4b->SetMinimum(0);
      ratio_3b->SetMinimum(0);
      ratio_2b->SetMinimum(0);
      ratio_4b->SetMaximum(0.4);
      ratio_3b->SetMaximum(0.4);
      ratio_2b->SetMaximum(0.4);
      
      
      
      thecanvas->cd(1);
      ratio_2b_data->Draw("hist e same");
      thecanvas->cd(2);
      ratio_3b_data->Draw("hist e same");
      thecanvas->cd(3);
      ratio_4b_data->Draw("hist e same");
    }
    thecanvas->SaveAs(TString("higgs_SRSB_shapes_")+njetsdesc[ijets]+".eps");
    
    //    } //end of loop over jet multiplicity

    //now, plot the ratio of ratios: 4b/2b and 4b/3b
    //this is the 'kappa' factor
    TH1D* doubleRatio_4b2b=0;
    TH1D* doubleRatio_4b3b=0;
    TH1D* doubleRatio_3b2b=0;
    TH1D* doubleRatio_4b2b_data=0;
    TH1D* doubleRatio_4b3b_data=0;
    TH1D* doubleRatio_3b2b_data=0;
    //let's do 3b/2b as well

    doubleRatio_4b2b = (TH1D*) ratio_2b->Clone("doubleRatio_4b2b");
    doubleRatio_4b3b = (TH1D*) ratio_3b->Clone("doubleRatio_4b3b");
    doubleRatio_3b2b = (TH1D*) ratio_3b->Clone("doubleRatio_3b2b");

    doubleRatio_4b2b->Reset();
    doubleRatio_4b3b->Reset();
    doubleRatio_3b2b->Reset();

    doubleRatio_4b2b->Divide(ratio_4b,ratio_2b);
    doubleRatio_4b3b->Divide(ratio_4b,ratio_3b);
    doubleRatio_3b2b->Divide(ratio_3b,ratio_2b);

    doubleRatio_3b2b->SetLineColor(kGreen+3);
    doubleRatio_3b2b->SetMarkerColor(kGreen+3);

    setPadDimensions(1000,600);
    renewCanvas("",3);
    thecanvas->cd(1);
    doubleRatio_4b2b->Draw("hist e");
    doubleRatio_4b2b->SetMinimum(0);    doubleRatio_4b2b->SetMaximum(2);
    thecanvas->cd(2);
    doubleRatio_4b3b->Draw("hist e");
    doubleRatio_4b3b->SetMinimum(0);    doubleRatio_4b3b->SetMaximum(2);
    thecanvas->cd(3);
    doubleRatio_3b2b->Draw("hist e");
    doubleRatio_3b2b->SetMinimum(0);    doubleRatio_3b2b->SetMaximum(2);
    if (dodata_) {
      doubleRatio_4b2b_data = (TH1D*) ratio_2b_data->Clone("doubleRatio_4b2b_data");
      doubleRatio_4b3b_data = (TH1D*) ratio_3b_data->Clone("doubleRatio_4b3b_data");
      doubleRatio_3b2b_data = (TH1D*) ratio_2b_data->Clone("doubleRatio_3b2b_data");

      doubleRatio_4b2b_data->Reset();
      doubleRatio_4b3b_data->Reset();
      doubleRatio_3b2b_data->Reset();
      
      doubleRatio_4b2b_data->Divide(ratio_4b_data,ratio_2b_data);
      doubleRatio_4b3b_data->Divide(ratio_4b_data,ratio_3b_data);
      doubleRatio_3b2b_data->Divide(ratio_3b_data,ratio_2b_data);

    doubleRatio_3b2b_data->SetLineColor(kGreen-6);
    doubleRatio_3b2b_data->SetMarkerColor(kGreen-6);

      thecanvas->cd(1);
      doubleRatio_4b2b_data->Draw("hist e same");
      thecanvas->cd(2);
      doubleRatio_4b3b_data->Draw("hist e same");
      thecanvas->cd(3);
      doubleRatio_3b2b_data->Draw("hist e same");
    }

    return;


  //now plot the jetpt1,2,3,4 spectra after all cuts for:
  //    'low SB'
  //    'medium SB'
  //    'high SB'
  //    'SR'
  // split by nb==2 and nb==4
  setStackMode(false,true,false); 
  savePlots_=true;

  for (int kk=0;kk<4;kk++) {
    clearSamples();
    TString filedesc;
    if (kk==0) {
      addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim:CSVbest3<0.679",kBlue,"t#bar{t} (low SB, 2b)");
      addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim:CSVbest3>0.679&&CSVbest4>0.244",kBlack,"t#bar{t} (low SB, 4b)");
      filedesc="lowSB";
      selection_ =  baseline&&trigger&&zl&&isotk&&(njets4||njets5)&&drmax&&metsigloose && higgsSB_avl;
    }
    else if (kk==1) {
      addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim:CSVbest3<0.679",kBlue,"t#bar{t} (med SB, 2b)");
      addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim:CSVbest3>0.679&&CSVbest4>0.244",kBlack,"t#bar{t} (med SB, 4b)");
      filedesc="medSB";
      selection_ =  baseline&&trigger&&zl&&isotk&&(njets4||njets5)&&drmax&&metsigloose && higgsSR_av &&higgsSB_d;

    }
    else if (kk==2) {
      addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim:CSVbest3<0.679",kBlue,"t#bar{t} (high SB, 2b)");
      addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim:CSVbest3>0.679&&CSVbest4>0.244",kBlack,"t#bar{t} (high SB, 4b)");
      filedesc="highSB";
      selection_ =  baseline&&trigger&&zl&&isotk&&(njets4||njets5)&&drmax&&metsigloose && higgsSB_avh;

    }
    else if (kk==3) {
      addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim:CSVbest3<0.679",kBlue,"t#bar{t} (SR, 2b)");
      addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim:CSVbest3>0.679&&CSVbest4>0.244",kBlack,"t#bar{t} (SR, 4b)");
      filedesc="SR";
      selection_ =  baseline&&trigger&&zl&&isotk&&(njets4||njets5)&&drmax&&metsigloose && higgsSR;
    }


    nbins=15; low=0; high=300;
    var="jetpt1"; xtitle=var;
    drawPlots(var,nbins,low,high,xtitle,"Events", TString("higgs_SRSBshapes_jetpt1_")+filedesc);
    var="jetpt2"; xtitle=var;
    drawPlots(var,nbins,low,high,xtitle,"Events", TString("higgs_SRSBshapes_jetpt2_")+filedesc);
    var="jetpt3"; xtitle=var;
    drawPlots(var,nbins,low,high,xtitle,"Events", TString("higgs_SRSBshapes_jetpt3_")+filedesc);
    var="jetpt4"; xtitle=var;
    drawPlots(var,nbins,low,high,xtitle,"Events", TString("higgs_SRSBshapes_jetpt4_")+filedesc);
  }

}

void higgs_ttbarjetmult1() {

  initHiggsSamples(false,"ttbar hh250 hh400"); //skim has triggers -- can't use it!

  int nbins;
  float low,high;
  TString var,xtitle;


  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0"; 
  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";
  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0&&nIsoTracks5_005_03<2";

  TCut njets="njets20==4 || njets20==5";
  TCut btagpre="CSVbest2>0.898";
  TCut btag = "CSVbest3>0.679 && CSVbest4>0.244";

  TCut higgsmass = "((0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)) && abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";

  TCut drmax = "deltaRmax_hh<2.4";
  TCut drmin = "deltaRmin_hh<1.9";

  TCut metsigloose = "METsig>50";

  setSampleWeightFactor("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66","ttbarDecayCode>=3&&ttbarDecayCode<=4");

  selection_ = baseline&&zl&&isotk;

  var="njets20"; xtitle="njets20";
  nbins=13; low= 0; high=13;
  setStackMode(false,true,false); //no stack
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_njetsStudy_baselineOnly_1e1m",0);


  selection_ = baseline;

  var="njets20"; xtitle="njets20";
  nbins=13; low= 0; high=13;
  setStackMode(false,true,false); //no stack
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_njetsStudy_superBaselineOnly_1e1m",0);

  // -- neither of these plots turned out to be very informative
  selection_ = baseline&&zl&&isotk&&TCut("jetpt2>=70");
  var="njets20"; xtitle="njets20";
  nbins=13; low= 0; high=13;
  setStackMode(false,true,false); //no stack
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_njetsStudy_jetpt2cut_1e1m",0);
  //this had the opposite effect from what i expected. mostly pushed ttbar higher

  // -- neither of these plots turned out to be very informative
  selection_ = baseline&&zl&&isotk&&TCut("jetpt2<70");
  var="njets20"; xtitle="njets20";
  nbins=13; low= 0; high=13;
  setStackMode(false,true,false); //no stack
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_njetsStudy_jetpt2cut2_1e1m",0);

  // -- now add trigger and partial btag (only 2 CSVT)
  // NOTE THAT THIS PLOT IS BUGGED -- can't use all 3 triggers in an AND (or you can, but you don't want to)
  selection_ = baseline&&zl&&isotk&&triggerJetMET&&triggerJets&&triggerMET && btagpre;
  var="njets20"; xtitle="njets20";
  nbins=13; low= 0; high=13;
  setStackMode(false,true,false); //no stack
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_njetsStudy_triggerAndCSVT_1e1m",0);


}

void higgs_Wmass() {
  initHiggsSamples(true,"ttbarjoint hh250 hh400"); //skim has triggers -- can't use it!

  int nbins;
  float low,high;
  TString var,xtitle;

  savePlots_=false;
  setQuiet(false);
  //dummy plot for counting events


  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0"; 
  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";
  //  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger=triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0&&nIsoTracks5_005_03<2";

  TCut njets="njets20==4 || njets20==5";
  TCut btagpre="CSVbest2>0.898";
  TCut btag = "CSVbest3>0.679 && CSVbest4>0.244";

  TCut higgsmass = "((0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)) && abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";

  TCut drmax = "deltaRmax_hh<2.4";
  //  TCut drmin = "deltaRmin_hh<1.9"; //not used

  TCut metsigloose = "METsig>50";


  selection_ = baseline &&trigger&& zl&&isotk&&njets&&btagpre&&btag&&higgsmass&&drmax&&metsigloose;

  var="higgsWCandMass"; xtitle=var;
  nbins=20; low= 0; high=200;
  setStackMode(false,false,false); //no stack
  setLogY(false);

  drawPlots(var,nbins,low,high,xtitle,"Events", "wmass",0);


}

void higgs_triggerImportance() {
  initHiggsSamples(false,"znunu ttbar hh150 hh200 hh250 hh400 wbb ttv"); //skim has triggers -- can't use it!

  int nbins;
  float low,high;
  TString var,xtitle;

  savePlots_=false;
  setQuiet(true);
  //dummy plot for counting events
  var="HT"; xtitle="HT";
  nbins=20; low= 0; high=1e9;
  setStackMode(false,false,false); //no stack
  setLogY(false);


  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0"; 
  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";
  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut triggerRA2b = "passMC_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80==1";

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";

  TCut njets="njets20==4 || njets20==5";
  TCut btagpre="CSVbest2>0.898";
  TCut btag = "CSVbest3>0.679 && CSVbest4>0.244";

  TCut higgsmass = "((0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)) && abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";

  TCut drmax = "deltaRmax_hh<2.4";
  //  TCut drmin = "deltaRmin_hh<1.9"; //not used

  TCut metsigloose = "METsig>50";


  cout<<"--- number passing selection, no trigger req'd ---"<<endl;
  selection_ = baseline && zl&&isotk&&njets&&btagpre&&btag&&higgsmass&&drmax&&metsigloose;
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
  cout<<"totalsm  = "<<jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm"),true,false)<<endl;
  cout<<"hh 150   = "<<jmt::format_nevents(getIntegral("TChihh_150_v68"),getIntegralErr("TChihh_150_v68"),true,false)<<endl;
  cout<<"hh 200   = "<<jmt::format_nevents(getIntegral("TChihh_200_v68"),getIntegralErr("TChihh_200_v68"),true,false)<<endl;
  cout<<"hh 250   = "<<jmt::format_nevents(getIntegral("TChihh_250_v68"),getIntegralErr("TChihh_250_v68"),true,false)<<endl;
  cout<<"hh 400   = "<<jmt::format_nevents(getIntegral("TChihh_400_v68"),getIntegralErr("TChihh_400_v68"),true,false)<<endl;

  double ntotalsm_0 = getIntegral("totalsm");
  double n150_0 = getIntegral("TChihh_150_v68");
  double n200_0 = getIntegral("TChihh_200_v68");
  double n250_0 = getIntegral("TChihh_250_v68");
  double n400_0 = getIntegral("TChihh_400_v68");

  cout<<"--- number passing selection and any (JetMETb||JetsOnly||METonly) trigger ---"<<endl;
  selection_ = baseline && zl&&isotk&&njets&&btagpre&&btag&&higgsmass&&drmax&&metsigloose && (triggerJetMET ||triggerJets||triggerMET);
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
  cout<<"totalsm  = "<<jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm"),true,false)<<"  "<<getIntegral("totalsm")/ntotalsm_0 <<endl;
  cout<<"hh 150   = "<<jmt::format_nevents(getIntegral("TChihh_150_v68"),getIntegralErr("TChihh_150_v68"),true,false)<<"  "<<getIntegral("TChihh_150_v68")/n150_0 <<endl;
  cout<<"hh 200   = "<<jmt::format_nevents(getIntegral("TChihh_200_v68"),getIntegralErr("TChihh_200_v68"),true,false)<<"  "<<getIntegral("TChihh_200_v68")/n200_0 <<endl;
  cout<<"hh 250   = "<<jmt::format_nevents(getIntegral("TChihh_250_v68"),getIntegralErr("TChihh_250_v68"),true,false)<<"  "<<getIntegral("TChihh_250_v68")/n250_0<<endl;
  cout<<"hh 400   = "<<jmt::format_nevents(getIntegral("TChihh_400_v68"),getIntegralErr("TChihh_400_v68"),true,false)<<"  "<<getIntegral("TChihh_400_v68")/n400_0<<endl;

  cout<<"--- number passing selection and only the Jet trigger ---"<<endl;
  selection_ = baseline && zl&&isotk&&njets&&btagpre&&btag&&higgsmass&&drmax&&metsigloose && ((!triggerJetMET) &&triggerJets&&(!triggerMET));
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
  cout<<"totalsm  = "<<jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm"),true,false)<<"  "<<getIntegral("totalsm")/ntotalsm_0 <<endl;
  cout<<"hh 150   = "<<jmt::format_nevents(getIntegral("TChihh_150_v68"),getIntegralErr("TChihh_150_v68"),true,false)<<"  "<<getIntegral("TChihh_150_v68")/n150_0 <<endl;
  cout<<"hh 200   = "<<jmt::format_nevents(getIntegral("TChihh_200_v68"),getIntegralErr("TChihh_200_v68"),true,false)<<"  "<<getIntegral("TChihh_200_v68")/n200_0 <<endl;
  cout<<"hh 250   = "<<jmt::format_nevents(getIntegral("TChihh_250_v68"),getIntegralErr("TChihh_250_v68"),true,false)<<"  "<<getIntegral("TChihh_250_v68")/n250_0<<endl;
  cout<<"hh 400   = "<<jmt::format_nevents(getIntegral("TChihh_400_v68"),getIntegralErr("TChihh_400_v68"),true,false)<<"  "<<getIntegral("TChihh_400_v68")/n400_0<<endl;

  cout<<"--- number passing selection and the Jet trigger ---"<<endl;
  selection_ = baseline && zl&&isotk&&njets&&btagpre&&btag&&higgsmass&&drmax&&metsigloose && triggerJets;
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
  cout<<"totalsm  = "<<jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm"),true,false)<<"  "<<getIntegral("totalsm")/ntotalsm_0 <<endl;
  cout<<"hh 150   = "<<jmt::format_nevents(getIntegral("TChihh_150_v68"),getIntegralErr("TChihh_150_v68"),true,false)<<"  "<<getIntegral("TChihh_150_v68")/n150_0 <<endl;
  cout<<"hh 200   = "<<jmt::format_nevents(getIntegral("TChihh_200_v68"),getIntegralErr("TChihh_200_v68"),true,false)<<"  "<<getIntegral("TChihh_200_v68")/n200_0 <<endl;
  cout<<"hh 250   = "<<jmt::format_nevents(getIntegral("TChihh_250_v68"),getIntegralErr("TChihh_250_v68"),true,false)<<"  "<<getIntegral("TChihh_250_v68")/n250_0<<endl;
  cout<<"hh 400   = "<<jmt::format_nevents(getIntegral("TChihh_400_v68"),getIntegralErr("TChihh_400_v68"),true,false)<<"  "<<getIntegral("TChihh_400_v68")/n400_0<<endl;

  cout<<"--- number passing selection and the JetMET trigger ---"<<endl;
  selection_ = baseline && zl&&isotk&&njets&&btagpre&&btag&&higgsmass&&drmax&&metsigloose && triggerJetMET;
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
  cout<<"totalsm  = "<<jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm"),true,false)<<"  "<<getIntegral("totalsm")/ntotalsm_0 <<endl;
  cout<<"hh 150   = "<<jmt::format_nevents(getIntegral("TChihh_150_v68"),getIntegralErr("TChihh_150_v68"),true,false)<<"  "<<getIntegral("TChihh_150_v68")/n150_0 <<endl;
  cout<<"hh 200   = "<<jmt::format_nevents(getIntegral("TChihh_200_v68"),getIntegralErr("TChihh_200_v68"),true,false)<<"  "<<getIntegral("TChihh_200_v68")/n200_0 <<endl;
  cout<<"hh 250   = "<<jmt::format_nevents(getIntegral("TChihh_250_v68"),getIntegralErr("TChihh_250_v68"),true,false)<<"  "<<getIntegral("TChihh_250_v68")/n250_0<<endl;
  cout<<"hh 400   = "<<jmt::format_nevents(getIntegral("TChihh_400_v68"),getIntegralErr("TChihh_400_v68"),true,false)<<"  "<<getIntegral("TChihh_400_v68")/n400_0<<endl;

  cout<<"--- number passing selection and the MET trigger ---"<<endl;
  selection_ = baseline && zl&&isotk&&njets&&btagpre&&btag&&higgsmass&&drmax&&metsigloose && triggerMET;
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
  cout<<"totalsm  = "<<jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm"),true,false)<<"  "<<getIntegral("totalsm")/ntotalsm_0 <<endl;
  cout<<"hh 150   = "<<jmt::format_nevents(getIntegral("TChihh_150_v68"),getIntegralErr("TChihh_150_v68"),true,false)<<"  "<<getIntegral("TChihh_150_v68")/n150_0 <<endl;
  cout<<"hh 200   = "<<jmt::format_nevents(getIntegral("TChihh_200_v68"),getIntegralErr("TChihh_200_v68"),true,false)<<"  "<<getIntegral("TChihh_200_v68")/n200_0 <<endl;
  cout<<"hh 250   = "<<jmt::format_nevents(getIntegral("TChihh_250_v68"),getIntegralErr("TChihh_250_v68"),true,false)<<"  "<<getIntegral("TChihh_250_v68")/n250_0<<endl;
  cout<<"hh 400   = "<<jmt::format_nevents(getIntegral("TChihh_400_v68"),getIntegralErr("TChihh_400_v68"),true,false)<<"  "<<getIntegral("TChihh_400_v68")/n400_0<<endl;

  cout<<"--- number passing selection and either of the MET triggers ---"<<endl;
  selection_ = baseline && zl&&isotk&&njets&&btagpre&&btag&&higgsmass&&drmax&&metsigloose && (triggerMET||triggerJetMET);
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
  cout<<"totalsm  = "<<jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm"),true,false)<<"  "<<getIntegral("totalsm")/ntotalsm_0 <<endl;
  cout<<"hh 150   = "<<jmt::format_nevents(getIntegral("TChihh_150_v68"),getIntegralErr("TChihh_150_v68"),true,false)<<"  "<<getIntegral("TChihh_150_v68")/n150_0 <<endl;
  cout<<"hh 200   = "<<jmt::format_nevents(getIntegral("TChihh_200_v68"),getIntegralErr("TChihh_200_v68"),true,false)<<"  "<<getIntegral("TChihh_200_v68")/n200_0 <<endl;
  cout<<"hh 250   = "<<jmt::format_nevents(getIntegral("TChihh_250_v68"),getIntegralErr("TChihh_250_v68"),true,false)<<"  "<<getIntegral("TChihh_250_v68")/n250_0<<endl;
  cout<<"hh 400   = "<<jmt::format_nevents(getIntegral("TChihh_400_v68"),getIntegralErr("TChihh_400_v68"),true,false)<<"  "<<getIntegral("TChihh_400_v68")/n400_0<<endl;

  cout<<"--- number passing selection and the RA2b JetMET trigger||PFMET150; METsig>50 ---"<<endl;

  selection_ = baseline && zl&&isotk&&njets&&btagpre&&btag&&higgsmass&&drmax&&metsigloose && (triggerMET||triggerRA2b);
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
  cout<<"totalsm  = "<<jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm"),true,false)<<"  "<<getIntegral("totalsm")/ntotalsm_0 <<endl;
  cout<<"hh 150   = "<<jmt::format_nevents(getIntegral("TChihh_150_v68"),getIntegralErr("TChihh_150_v68"),true,false)<<"  "<<getIntegral("TChihh_150_v68")/n150_0 <<endl;
  cout<<"hh 200   = "<<jmt::format_nevents(getIntegral("TChihh_200_v68"),getIntegralErr("TChihh_200_v68"),true,false)<<"  "<<getIntegral("TChihh_200_v68")/n200_0 <<endl;
  cout<<"hh 250   = "<<jmt::format_nevents(getIntegral("TChihh_250_v68"),getIntegralErr("TChihh_250_v68"),true,false)<<"  "<<getIntegral("TChihh_250_v68")/n250_0<<endl;
  cout<<"hh 400   = "<<jmt::format_nevents(getIntegral("TChihh_400_v68"),getIntegralErr("TChihh_400_v68"),true,false)<<"  "<<getIntegral("TChihh_400_v68")/n400_0<<endl;

  cout<<"--- number passing selection and the MET||JetMETb||Ra2b trigger ---"<<endl;
  selection_ = baseline && zl&&isotk&&njets&&btagpre&&btag&&higgsmass&&drmax&&metsigloose && (triggerMET||triggerJetMET||triggerRA2b);
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
  cout<<"totalsm  = "<<jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm"),true,false)<<"  "<<getIntegral("totalsm")/ntotalsm_0 <<endl;
  cout<<"hh 150   = "<<jmt::format_nevents(getIntegral("TChihh_150_v68"),getIntegralErr("TChihh_150_v68"),true,false)<<"  "<<getIntegral("TChihh_150_v68")/n150_0 <<endl;
  cout<<"hh 200   = "<<jmt::format_nevents(getIntegral("TChihh_200_v68"),getIntegralErr("TChihh_200_v68"),true,false)<<"  "<<getIntegral("TChihh_200_v68")/n200_0 <<endl;
  cout<<"hh 250   = "<<jmt::format_nevents(getIntegral("TChihh_250_v68"),getIntegralErr("TChihh_250_v68"),true,false)<<"  "<<getIntegral("TChihh_250_v68")/n250_0<<endl;
  cout<<"hh 400   = "<<jmt::format_nevents(getIntegral("TChihh_400_v68"),getIntegralErr("TChihh_400_v68"),true,false)<<"  "<<getIntegral("TChihh_400_v68")/n400_0<<endl;


  return;

  cout<<"--- number passing selection, no trigger req'd; METsig>30 ---"<<endl;
  selection_ = baseline && zl&&isotk&&njets&&btagpre&&btag&&higgsmass&&drmax&&TCut("METsig>30");
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
  cout<<"totalsm  = "<<jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm"),true,false)<<endl;
  cout<<"hh 150   = "<<jmt::format_nevents(getIntegral("TChihh_150_v68"),getIntegralErr("TChihh_150_v68"),true,false)<<endl;
  cout<<"hh 200   = "<<jmt::format_nevents(getIntegral("TChihh_200_v68"),getIntegralErr("TChihh_200_v68"),true,false)<<endl;
  cout<<"hh 250   = "<<jmt::format_nevents(getIntegral("TChihh_250_v68"),getIntegralErr("TChihh_250_v68"),true,false)<<endl;
  cout<<"hh 400   = "<<jmt::format_nevents(getIntegral("TChihh_400_v68"),getIntegralErr("TChihh_400_v68"),true,false)<<endl;
  ntotalsm_0 = getIntegral("totalsm");
  n150_0 = getIntegral("TChihh_150_v68");
  n200_0 = getIntegral("TChihh_200_v68");
  n250_0 = getIntegral("TChihh_250_v68");
  n400_0 = getIntegral("TChihh_400_v68");

  cout<<"--- number passing selection and either of the MET triggers; METsig>30 ---"<<endl;
  selection_ = baseline && zl&&isotk&&njets&&btagpre&&btag&&higgsmass&&drmax&&TCut("METsig>30") && (triggerMET||triggerJetMET);
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
  cout<<"totalsm  = "<<jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm"),true,false)<<"  "<<getIntegral("totalsm")/ntotalsm_0 <<endl;
  cout<<"hh 150   = "<<jmt::format_nevents(getIntegral("TChihh_150_v68"),getIntegralErr("TChihh_150_v68"),true,false)<<"  "<<getIntegral("TChihh_150_v68")/n150_0 <<endl;
  cout<<"hh 200   = "<<jmt::format_nevents(getIntegral("TChihh_200_v68"),getIntegralErr("TChihh_200_v68"),true,false)<<"  "<<getIntegral("TChihh_200_v68")/n200_0 <<endl;
  cout<<"hh 250   = "<<jmt::format_nevents(getIntegral("TChihh_250_v68"),getIntegralErr("TChihh_250_v68"),true,false)<<"  "<<getIntegral("TChihh_250_v68")/n250_0<<endl;
  cout<<"hh 400   = "<<jmt::format_nevents(getIntegral("TChihh_400_v68"),getIntegralErr("TChihh_400_v68"),true,false)<<"  "<<getIntegral("TChihh_400_v68")/n400_0<<endl;



}

void higgs_Nminus1() {
  initHiggsSamples(false);

  int nbins;
  float low,high;
  TString var,xtitle;

  savePlots_=true;
  setQuiet(false);

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  //20 June -- update
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0";
  //20 june -- remove the all jet trigger: ||passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_PFMET150==1";
  TCut zl = "nMuons==0&&nElectrons==0";
  TCut tauveto="nTausLoose==0";
  TCut isotk1="nIsoTracks15_005_03==0";
  //  TCut isotk2="nIsoTracks5_005_03<2"; //20 June -- don't use this

  TCut njets="njets20==4 || njets20==5";
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut higgsmass = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsdiff = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";

  TCut drmax = "deltaRmax_hh<2.4";
  TCut drmin = "deltaRmin_hh<1.9";

  TCut metsig = "METsig>30";

  setStackMode(true,false,false); //stack, norm, labels

  stackSignal_=false; //let's unstack signal

  //full selection (updated to not use drmin, isotk2, or jet trigger)
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets && btag2 && btag3 && btag4 && higgsmass && higgsdiff && drmax && metsig;
  var="bjetpt1"; xtitle="lead b jet pT";
  nbins=20; low= 0; high=400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_bjetpt1",0);

  //still full selection
  var="jetpt1"; xtitle="lead jet pT";
  nbins=20; low= 0; high=200; //need 10 GeV Bins
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_jetpt1",0);

  var="jetpt2"; xtitle="2nd jet pT";
  nbins=20; low= 0; high=200; //need 10 GeV Bins
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_jetpt2",0);

  var="jetpt3"; xtitle="3rd jet pT";
  nbins=20; low= 0; high=200; //need 10 GeV Bins
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_jetpt3",0);

  var="jetpt4"; xtitle="4th jet pT";
  nbins=20; low= 0; high=200; //need 10 GeV Bins
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_jetpt4",0);

  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets && btag2 && btag3 && btag4 && higgsmass && higgsdiff && drmax && metsig;
  var="njets20_5p0-njets20"; xtitle="number of forward jets (pT>20)";
  nbins=4; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_njets20forward",0);

  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets && btag2 && btag3 && btag4 && higgsmass && higgsdiff && drmax && metsig;
  var="deltaPhi_hh"; xtitle="#Delta #phi (h,h)";
  nbins=10; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_hhDeltaPhi",0);

  //n leptons
  selection_=baseline&&triggers &&  tauveto && isotk1 && njets && btag2 && btag3 && btag4 && higgsmass && higgsdiff && drmax && metsig;
  var="nMuons+nElectrons"; xtitle="e + #mu";
  nbins=4; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_ePlusMu",0);

  //POG taus
  selection_=baseline&&triggers && zl&&  isotk1 && njets && btag2 && btag3 && btag4 && higgsmass && higgsdiff && drmax && metsig;
  var="nTausLoose"; xtitle="Loose taus";
  nbins=3; low= 0; high=3;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_looseTaus",0);

  //iso tracks
  selection_=baseline&&triggers && zl&& tauveto && njets && btag2 && btag3 && btag4 && higgsmass && higgsdiff && drmax && metsig;
  var="nIsoTracks15_005_03"; xtitle="15 GeV iso tracks";
  nbins=4; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_isoTk15",0);

  //iso tracks
  selection_=baseline&&triggers && zl&& tauveto  && isotk1 && njets && btag2 && btag3 && btag4 && higgsmass && higgsdiff && drmax && metsig;
  var="nIsoTracks5_005_03"; xtitle="5 GeV iso tracks";
  nbins=4; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_isoTk5",0);

  //jet multiplicity
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && btag2 && btag3 && btag4 && higgsmass && higgsdiff && drmax  && metsig;
  var="njets20"; xtitle="jet multiplicity (20 GeV)";
  nbins=9; low= 0; high=9;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_njets20",0);

  //3rd b tag
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets && btag2  && higgsmass && higgsdiff && drmax && metsig;
  var="CSVbest3"; xtitle="3rd CSV value";
  nbins=20; low= 0; high=1;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_CSVbest3",0);

  //4th btag
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets && btag2 && btag3 && higgsmass && higgsdiff && drmax && metsig;
  var="CSVbest4"; xtitle="4th CSV value";
  nbins=20; low= 0; high=1;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_CSVbest4",0);

  //higgs mass
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets && btag2 && btag3 && btag4  && higgsdiff && drmax && metsig;
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="average higgs mass";
  nbins=20; low= 0; high=200;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_higgsmass",0);

  //higgs mass difference
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets && btag2 && btag3 && btag4 && higgsmass  && drmax && metsig;
  var="abs(higgsMbb1MassDiff-higgsMbb2MassDiff)"; xtitle="#Delta higgs mass";
  nbins=20; low= 0; high=80;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_higgsmassdiff",0);

  //drmax
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets && btag2 && btag3 && btag4 && higgsmass && higgsdiff && metsig;
  var="deltaRmax_hh"; xtitle="max #Delta R";
  nbins=20; low= 0; high=6;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_maxDR",0);

  //drmin
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets && btag2 && btag3 && btag4 && higgsmass && higgsdiff  && drmax && metsig;
  var="deltaRmin_hh"; xtitle="min #Delta R";
  nbins=20; low= 0; high=6;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_minDR",0);

  //full selection except METsig
  selection_=baseline&&triggers && zl&& tauveto && isotk1&& njets && btag2 && btag3 && btag4 && higgsmass && higgsdiff && drmax;
  var="METsig"; xtitle="MET significance";
  nbins=20; low= 0; high=200;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_METsig",0);

  //full selection except METsig
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets && btag2 && btag3 && btag4 && higgsmass && higgsdiff && drmax;
  var="MET"; xtitle="MET";
  nbins=20; low= 0; high=300;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_MET",0);

}

void higgs_printCutflowTable( const  bool unweighted=false) {
  if (unweighted)  initHiggsSamples(false,"ttbar hh250 hh400"); //reduced list of samples
  else   initHiggsSamples(false);

  int nbins;
  float low,high;
  TString var,xtitle;

  if (unweighted) {
    nameOfEventWeight_="(1)";
    usePUweight_=false;
    lumiScale_=1;
    resetSampleWeightFactors();
    resetSampleScaleFactors();
  }

  savePlots_=false;
  setQuiet(true);
  //dummy plot
  var="HT"; xtitle="HT";
  nbins=20; low= 0; high=1e9;
  setStackMode(false,false,false); //no stack
  setLogY(false);


  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0"; 
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_PFMET150==1";
  //remove this trigger for now passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1
  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2"; //remove this

  TCut njets="njets30>=4&&njets30<=5";
  TCut btagpre="CSVbest2>0.898";
  TCut btag = "CSVbest3>0.679 && CSVbest4>0.244";

  TCut higgsmass = "((0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)) && abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";

  TCut drmax = "deltaRmax_hh<2.4";
  TCut drmin = "deltaRmin_hh<1.9";

  TCut metsigloose = "METsig>30";
  TCut metsig1 = "METsig>50";
  TCut metsig2 = "METsig>100";
  TCut metsig3 = "METsig>150";

  //after pre-selection, compare keith's variable to mine
  std::vector<pair<TString, TCut> > cuts;
  if (true||unweighted) cuts.push_back(make_pair("PV",TCut("cutPV==1")));  
  else  cuts.push_back(make_pair( "Cleanup",baseline));

  cuts.push_back(make_pair( "Trigger",triggers));
  cuts.push_back(make_pair( "njets 4-5 (30)",njets));
  //cuts.push_back(make_pair( "njets 4-5 (30 GeV)",TCut("njets30==4||njets30==5")));

  cuts.push_back(make_pair( "2 CSVT",btagpre));

  //  cuts.push_back(make_pair("jet pT 2 >50, pT 4>30",TCut("jetpt2>50&&jetpt4>30")));

  cuts.push_back(make_pair("MET sig>30",metsigloose));

  cuts.push_back(make_pair( "Lepton vetoes",zl));
  cuts.push_back(make_pair( "IsoTk vetoes",isotk));
  cuts.push_back(make_pair( "b-tagging",btag));
  cuts.push_back(make_pair( "Higgs masses",higgsmass));
  cuts.push_back(make_pair( "max #Delta R",drmax));
  //  cuts.push_back(make_pair( "min #Delta R",drmin));
  
  cuts.push_back(make_pair( "MET sig > 50",metsig1));
  cuts.push_back(make_pair( "MET sig > 100",metsig2));
  cuts.push_back(make_pair( "MET sig > 150",metsig3));
  /*
  cuts.push_back(make_pair( "MET > 125",TCut("MET>125")));
  cuts.push_back(make_pair( "MET > 175",TCut("MET>175")));
  cuts.push_back(make_pair( "MET > 225",TCut("MET>225")));
  */
  //first print a header row
  cout<<"|* Cut \t\t*|* ";
  bool foundsignals=false;
  for (unsigned int isample = 0; isample<samples_.size(); isample++) {
    if ( !isSampleSM(samples_.at(isample)) && !foundsignals) {
      cout<<" Total SM *|* ";
      foundsignals=true;
    }
    cout<<sampleLabel_[samples_.at(isample)];
    cout<<" *|* ";
  }
  cout<<endl;

  TCut totalcut = cuts.at(0).second;
  for (unsigned int icut=0; icut<cuts.size(); icut++) {
    if (icut != 0) totalcut = totalcut&&cuts.at(icut).second;
    selection_=totalcut;

    drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
    //now we counted events for all of the samples for this set of cuts.
    cout<<"| ";
    cout<<cuts.at(icut).first<<" \t\t | ";
    foundsignals=false;
    for (unsigned int isample = 0; isample<samples_.size(); isample++) {

      TString thissamplename = samples_.at(isample);
      if ( !isSampleSM(thissamplename) && !foundsignals) {
	if (unweighted)	cout<<setprecision(10)<<getIntegral("totalsm")<<" |";
	else	cout<<jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm"),false,false)<<" |";
	foundsignals=true;
      }
      double n = getIntegral(thissamplename);
      double e = getIntegralErr(thissamplename); 
      if (unweighted)     cout<<setprecision(10)<<n<<" |";
      else     cout<<jmt::format_nevents(n,e,false,false)<<" |";
    }
    cout<<endl;
  }

}

void higgmb_cutflowE_bbangles() {
  // -- goal -- check keith's maxDeltaPhi_bb_bb against the deltaRmax_hh variable

  addToSamplesAll("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim");
  addToSamplesAll("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1613_v66-skim");
  addToSamplesAll("TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1606_v66-skim");

  addToSamplesAll("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66");
  addToSamplesAll("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1607_v66");
  addToSamplesAll("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1524_v66");
  addToSamplesAll("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1594_v66");
  addToSamplesAll("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1523_v66");
  addToSamplesAll("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1602_v66");

  addToSamplesAll("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66");
  addToSamplesAll("QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1513_v66");
  addToSamplesAll("QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1578_v66");
  addToSamplesAll("QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1561_v66");
  addToSamplesAll("QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1603_v66");
  addToSamplesAll("QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1585_v66");
  addToSamplesAll("QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1560_v66");
  addToSamplesAll("QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_v3_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1609_v66");
  addToSamplesAll("QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1515_v66");
  addToSamplesAll("QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1516_v66");
  addToSamplesAll("QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1559_v66");


  addToSamplesAll("TChihh_250_v68");
  addToSamplesAll("TChihh_400_v68");

  nameOfEventWeight_="weight3";//don't want an event weight

  inputPath="/cu6/joshmt/reducedTrees/v68_4k/";

  loadSamples(true,"ra2b2012");
  usePUweight_=true; //helps bring signal and background into alignment; not perfect but better than nothing

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0"; 
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1||passMC_PFMET150==1";
  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0&&nIsoTracks5_005_03<2";

  TCut njets="njets20==4 || njets20==5";
  TCut btag = "CSVbest2>0.898 && CSVbest3>0.679 && CSVbest4>0.244";

  TCut higgsmass = "((0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)) && abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";

  clearSamples();


  addSample("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66",kGreen,"QCD");
  chainSamples("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66","QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1513_v66");
  chainSamples("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66","QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1578_v66");
  chainSamples("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66","QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1561_v66");
  chainSamples("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66","QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1603_v66");
  chainSamples("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66","QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1585_v66");
  chainSamples("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66","QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1560_v66");
  chainSamples("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66","QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_v3_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1609_v66");
  chainSamples("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66","QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1515_v66");
  chainSamples("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66","QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1516_v66");
  chainSamples("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577_v66","QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1559_v66");


  addSample("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66",kOrange,"Z #rightarrow #nu#nu");
  chainSamples("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66","ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1607_v66");
  chainSamples("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66","ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1524_v66");
  chainSamples("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66","ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1594_v66");
  chainSamples("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66","ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1523_v66");
  chainSamples("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66","ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1602_v66");

  addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim",kAzure-3,"t#bar{t} (#geq 1l)");
  chainSamples("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim","TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1606_v66-skim");
  addSample("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1613_v66-skim",kYellow,"hadronic t#bar{t}");
  addSample("TChihh_250_v68",kRed,"hh 250");
  addSample("TChihh_400_v68",kRed-7,"hh 400");

  setSampleScaleFactor("TChihh_250_v68",0.09/9999.0); //sigma x BF / ngen
  setSampleScaleFactor("TChihh_400_v68",0.01/9999.0); //sigma x BF / ngen

  TCut metsig1 = "METsig>30";
  TCut metsig2 = "METsig>50";
  TCut metsig3 = "METsig>100";
  TCut metsig4 = "METsig>130";

  //after pre-selection, compare keith's variable to mine
  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmass&&metsig2;

  var="maxDeltaPhi_bb_bb"; xtitle="maxDeltaPhi_bb_bb";
  nbins=20; low= 0; high=6;
  setStackMode(false,false,false); //no stack
  setLogY(false);

  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_bbangles_metsig50_dpmax",0);
  drawEffVRej("TChihh_400_v68","totalsm","hh 400","total SM",false);
  TGraphAsymmErrors* dpmax400  =(TGraphAsymmErrors*)  effgraph->Clone("dpmax400");
  drawEffVRej("TChihh_250_v68","totalsm","hh 250","total SM",false);
  TGraphAsymmErrors* dpmax250  =(TGraphAsymmErrors*)  effgraph->Clone("dpmax250");

  var="deltaRmax_hh"; xtitle=var;
  nbins=20; low= 0; high=6;
  setStackMode(false,false,false); //no stack
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_bbangles_metsig50_drmax",0);
  drawEffVRej("TChihh_400_v68","totalsm","hh 400","total SM",false);
  TGraphAsymmErrors* drmax400  =(TGraphAsymmErrors*)  effgraph->Clone("drmax400");
  drawEffVRej("TChihh_250_v68","totalsm","hh 250","total SM",false);
  TGraphAsymmErrors* drmax250  =(TGraphAsymmErrors*)  effgraph->Clone("drmax250");

  renewCanvas();
  dpmax400->SetLineColor(kBlue);
  dpmax400->SetMarkerColor(kBlue);
  dpmax400->Draw("APL");
  drmax400->SetLineColor(kRed);
  drmax400->SetMarkerColor(kRed);
  drmax400->Draw("PL");

  renewCanvas();
  dpmax250->SetLineColor(kBlue);
  dpmax250->SetMarkerColor(kBlue);
  dpmax250->Draw("aPL");
  drmax250->SetLineColor(kRed);
  drmax250->SetMarkerColor(kRed);
  drmax250->Draw("PL");

    //confirm this for tigher metsig
  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmass&&metsig3;

  var="maxDeltaPhi_bb_bb"; xtitle="maxDeltaPhi_bb_bb";
  nbins=20; low= 0; high=6;
  setStackMode(false,false,false); //no stack
  setLogY(false);

  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_bbangles_metsig100_dpmax",0);
  drawEffVRej("TChihh_400_v68","totalsm","hh 400","total SM",false);
  TGraphAsymmErrors* dpmax400_2  =(TGraphAsymmErrors*)  effgraph->Clone("dpmax400_2");
  drawEffVRej("TChihh_250_v68","totalsm","hh 250","total SM",false);
  TGraphAsymmErrors* dpmax250_2  =(TGraphAsymmErrors*)  effgraph->Clone("dpmax250_2");

  var="deltaRmax_hh"; xtitle=var;
  nbins=20; low= 0; high=6;
  setStackMode(false,false,false); //no stack
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_bbangles_metsig100_drmax",0);
  drawEffVRej("TChihh_400_v68","totalsm","hh 400","total SM",false);
  TGraphAsymmErrors* drmax400_2  =(TGraphAsymmErrors*)  effgraph->Clone("drmax400_2");
  drawEffVRej("TChihh_250_v68","totalsm","hh 250","total SM",false);
  TGraphAsymmErrors* drmax250_2  =(TGraphAsymmErrors*)  effgraph->Clone("drmax250_2");

  renewCanvas();
  dpmax400_2->SetLineColor(kBlue);
  dpmax400_2->SetMarkerColor(kBlue);
  dpmax400_2->Draw("APL");
  drmax400_2->SetLineColor(kRed);
  drmax400_2->SetMarkerColor(kRed);
  drmax400_2->Draw("PL");

  renewCanvas();
  dpmax250_2->SetLineColor(kBlue);
  dpmax250_2->SetMarkerColor(kBlue);
  dpmax250_2->Draw("aPL");
  drmax250_2->SetLineColor(kRed);
  drmax250_2->SetMarkerColor(kRed);
  drmax250_2->Draw("PL");

  //change gears
  //MET versus METsig
  TCut drmax = "deltaRmax_hh<2.2";
  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmass&&drmax;

  var="METsig"; xtitle="MET significance";
  nbins=20; low= 0; high=200;
  setStackMode(false,false,false); //no stack
  setLogY(true);
  setPlotMinimum(0.1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_optcuts_metsig",0);
  setLogY(false);
  drawEffVRej("TChihh_400_v68","totalsm","hh 400","total SM",true);
  TGraphAsymmErrors* metsig400  =(TGraphAsymmErrors*)  effgraph->Clone("metsig400");
  drawEffVRej("TChihh_250_v68","totalsm","hh 250","total SM",true);
  TGraphAsymmErrors* metsig250  =(TGraphAsymmErrors*)  effgraph->Clone("metsig250");


  var="MET"; xtitle="MET";
  nbins=20; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  setLogY(true);
  setPlotMinimum(0.1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_optcuts_met",0);
  setLogY(false);
  drawEffVRej("TChihh_400_v68","totalsm","hh 400","total SM",true);
  TGraphAsymmErrors* met400  =(TGraphAsymmErrors*)  effgraph->Clone("met400");
  drawEffVRej("TChihh_250_v68","totalsm","hh 250","total SM",true);
  TGraphAsymmErrors* met250  =(TGraphAsymmErrors*)  effgraph->Clone("met250");

  renewCanvas();
  metsig250->SetMaximum(0.02);
  metsig250->SetLineColor(kRed);
  metsig250->SetMarkerColor(kRed);
  metsig250->Draw("APL");

  met250->SetLineColor(kBlue);
  met250->SetMarkerColor(kBlue);
  met250->Draw("PL");


  renewCanvas();
  metsig400->SetMaximum(0.005);
  metsig400->SetLineColor(kRed);
  metsig400->SetMarkerColor(kRed);
  metsig400->Draw("APL");

  met400->SetLineColor(kBlue);
  met400->SetMarkerColor(kBlue);
  met400->Draw("PL");


}

void higgsmb_cutflowE() {

  //24 May -- update to keith's latest cuts and trees

  //no qcd for now
  /*
  addToSamplesAll("QCD_Pt-1000to1400");
  addToSamplesAll("QCD_Pt-120to170");
  addToSamplesAll("QCD_Pt-1400to1800");
  addToSamplesAll("QCD_Pt-170to300-1");
  addToSamplesAll("QCD_Pt-170to300-2");
  addToSamplesAll("QCD_Pt-1800");
  addToSamplesAll("QCD_Pt-300to470-1");
  addToSamplesAll("QCD_Pt-300to470-2");
  addToSamplesAll("QCD_Pt-470to600");
  addToSamplesAll("QCD_Pt-600to800");
  addToSamplesAll("QCD_Pt-800to1000");
  */

  addToSamplesAll("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim");
  addToSamplesAll("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1613_v66-skim");
  addToSamplesAll("TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1606_v66-skim");

  addToSamplesAll("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66");
  addToSamplesAll("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1607_v66");
  addToSamplesAll("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1524_v66");
  addToSamplesAll("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1594_v66");
  addToSamplesAll("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1523_v66");
  addToSamplesAll("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1602_v66");

  addToSamplesAll("TChihh_250_v68");
  addToSamplesAll("TChihh_400_v68");

  nameOfEventWeight_="weight3";//don't want an event weight

  inputPath="/cu6/joshmt/reducedTrees/v68_4k/";

  loadSamples(true,"ra2b2012");
  usePUweight_=true; //helps bring signal and background into alignment; not perfect but better than nothing

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0"; 
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1||passMC_PFMET150==1";
  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0&&nIsoTracks5_005_03<2";

  TCut njets="njets20==4 || njets20==5";
  TCut btag = "CSVbest2>0.898 && CSVbest3>0.679 && CSVbest4>0.244";

  TCut higgsmassNew = "((0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)) && abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsmassOld = "higgsMbb2MassDiff>=100&&higgsMbb2MassDiff<=140&&higgsMbb1MassDiff>=100&&higgsMbb1MassDiff<=140";

  TCut bbdp = "maxDeltaPhi_bb_bb<2.0";
  TCut bbdr = "deltaRmax_hh<2.0";

  clearSamples();

  addSample("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66",kOrange,"Z #rightarrow #nu#nu");
  chainSamples("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66","ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1607_v66");
  chainSamples("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66","ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1524_v66");
  chainSamples("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66","ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1594_v66");
  chainSamples("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66","ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1523_v66");
  chainSamples("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66","ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1602_v66");

  addSample("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim",kGreen,"2l");
  addSample("TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1606_v66-skim",kBlue,"1l");
  addSample("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1613_v66-skim",kYellow,"0l");
  addSample("TChihh_250_v68",kRed,"hh 250");
  addSample("TChihh_400_v68",kRed-7,"hh 400");

  setSampleScaleFactor("TChihh_250_v68",0.09/9999.0); //sigma x BF / ngen
  setSampleScaleFactor("TChihh_400_v68",0.01/9999.0); //sigma x BF / ngen

  TCut metsig1 = "METsig>30";
  TCut metsig2 = "METsig>50";
  TCut metsig3 = "METsig>100";
  TCut metsig4 = "METsig>130";


  var="MET"; xtitle="MET";
  nbins=15; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  setPlotMinimum(0.1);
  setLogY(true);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassNew&&bbdp;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentest0",0);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassNew&&bbdp&&metsig1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentest1",0);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassNew&&bbdp&&metsig2;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentest2",0);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassNew&&bbdp&&metsig3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentest3",0);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassNew&&bbdp&&metsig4;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentest4",0);

  //
  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassOld&&bbdp;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentestOld0",0);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassOld&&bbdp&&metsig1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentestOld1",0);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassOld&&bbdp&&metsig2;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentestOld2",0);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassOld&&bbdp&&metsig3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentestOld3",0);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassOld&&bbdp&&metsig4;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentestOld4",0);



  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassNew&&bbdr;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentestDR0",0);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassNew&&bbdr&&metsig1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentestDR1",0);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassNew&&bbdr&&metsig2;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentestDR2",0);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassNew&&bbdr&&metsig3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentestDR3",0);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassNew&&bbdr&&metsig4;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentestDR4",0);

  //
  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassOld&&bbdr;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentestDROld0",0);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassOld&&bbdr&&metsig1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentestDROld1",0);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassOld&&bbdr&&metsig2;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentestDROld2",0);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassOld&&bbdr&&metsig3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentestDROld3",0);

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmassOld&&bbdr&&metsig4;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow_owentestDROld4",0);
  


}

void higgsmbb_quickmbbPlots() {

 nameOfEventWeight_="weight3";//don't want an event weight

  inputPath="/cu4/ra2b/reducedTrees/v68_4/";

  addToSamplesAll("TChihh250");
  addToSamplesAll("TChihh400");

  loadSamples(true,"ra2b2012");
  usePUweight_=true; //helps bring signal and background into alignment; not perfect but better than nothing

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  clearSamples();
  //  addSample("TTbarJets0",kGreen,"t#bar{t}");
  addSample("TChihh250:higgsMbb1MassDiff_correct==2",kRed,"hh 250 (hh correct)");
  addSample("TChihh400:higgsMbb1MassDiff_correct==2",kBlue,"hh 400 (hh correct)");

  setSampleScaleFactor("TChihh250",0.09/9999.0); //sigma x BF / ngen
  setSampleScaleFactor("TChihh400",0.01/9999.0); //sigma x BF / ngen

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  setStackMode(false,false,false); //no stack
  TCut baseline = "cutPV==1 &&passCleaning==1"; 
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05||passMC_PFMET150";
  TCut zl = "nMuons==0&&nElectrons==0";
  TCut isotk="nIsoTracks15_005_03==0";

  TCut njets="njets20==4 || njets20==5";
  TCut btag = "CSVbest2>=0.898 && CSVbest3>=0.679 && CSVbest4>=0.244";

  TCut metsig = "METsig>=60";
  //no mdp cut
  selection_=baseline && triggers &&zl &&isotk&&njets&&btag &&metsig;


  var="higgsMbb1MassDiff"; xtitle="mbb (1)";
  nbins=20; low= 0; high=200;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_massshape_mbb1",0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_massshape_mbb1",0);

  var="higgsMbb2MassDiff"; xtitle="mbb (2)";
  nbins=20; low= 0; high=200;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_massshape_mbb2",0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_massshape_mbb2",0);


}

void higgsmbb_backgroundStudies2() { //this time split the off-mbb peak into a couple regions

  addToSamplesAll("TTJets_FullLeptMGDecays");
  addToSamplesAll("TTJets_SemiLeptMGDecays");

  //use keith's cut flow.

  nameOfEventWeight_="weight3";//don't want an event weight

  inputPath="/cu4/ra2b/reducedTrees/v68_4/";

  loadSamples(true,"ra2b2012");
  usePUweight_=true; //helps bring signal and background into alignment; not perfect but better than nothing

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 &&passCleaning==1"; 
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05||passMC_PFMET150";
  TCut zl = "nMuons==0&&nElectrons==0";
  TCut isotk="nIsoTracks15_005_03==0";

  TCut njets="njets20==4 || njets20==5";
  TCut btag = "CSVbest2>=0.898 && CSVbest3>=0.679 && CSVbest4>=0.244";

  TCut metsigVLoose = "METsig>=40";
  TCut metVLoose = "MET>40";

  TCut mdp="minDeltaPhiN_asin>=4";

  //this is an inelegant trick to manipulate the chains
  //it takes advantage of the fact that clearSamples() does not reset the chains.
  //it's a dirty hack but it works
  clearSamples();

  addSample("TTJets_FullLeptMGDecays");
  chainSamples("TTJets_FullLeptMGDecays","TTJets_SemiLeptMGDecays");

  //  for (int iii=0;iii<1; iii++) {
    
  int iii=0;
    clearSamples();

    map<TString,TString> subsamples;
    subsamples["SR"] = "TTJets_FullLeptMGDecays:higgsMbb2MassDiff>=100&&higgsMbb2MassDiff<=140&&higgsMbb1MassDiff>=100&&higgsMbb1MassDiff<=140";
    subsamples["SBhigh"]=   "TTJets_FullLeptMGDecays:higgsMbb2MassDiff>=140&&higgsMbb1MassDiff>=140"; 
    subsamples["SBlow"] = "TTJets_FullLeptMGDecays:higgsMbb2MassDiff<100&&higgsMbb1MassDiff<100";
    subsamples["SBother"]="TTJets_FullLeptMGDecays:(higgsMbb2MassDiff>=100&&higgsMbb2MassDiff<=140&&(higgsMbb1MassDiff<100||higgsMbb1MassDiff>140))||(higgsMbb1MassDiff>=100&&higgsMbb1MassDiff<=140&&(higgsMbb2MassDiff<100||higgsMbb2MassDiff>140))";
    TString thissample="";
    if (iii==0){
      addSample(subsamples["SR"],kGreen,"t#bar{t} (SR)");
      addSample(subsamples["SBhigh"],kBlue,"t#bar{t} (high SB)");
      addSample(subsamples["SBlow"],kRed,"t#bar{t} (low SB)");
      addSample(subsamples["SBother"],kMagenta,"t#bar{t} (other SB)");
      thissample="ttbar";
    } else return;
    //////
    
    //  addSample("TTbarJets0:higgsMbb2MassDiff>=100&&higgsMbb2MassDiff<=140&&higgsMbb1MassDiff>=100&&higgsMbb1MassDiff<=140",kGreen,"t#bar{t} (SR)");
    //  addSample("TTbarJets0:higgsMbb2MassDiff<100||higgsMbb2MassDiff>140||higgsMbb1MassDiff<100||higgsMbb1MassDiff>140",kBlue,"t#bar{t} (off SR)");
    
    //    fitRatio_="pol1";
    selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&mdp &&metsigVLoose;
    var="METsig"; xtitle="MET significance";
    nbins=15; low= 0; high=300;
    setStackMode(false,false,false); //no stack
    setLogY(true);
    doRatio_=false; ratioMin = 0; ratioMax = 0.5;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_metsig-shapes2-"+thissample,0);
    
    setLogY(true);
    setStackMode(false,true,false); //norm
    doRatio_=false; ratioMin = 0; ratioMax = 2.2;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_metsig-shapes2-"+thissample,0);
    

    //now make ratio plots
    TH1D* SRhist = getHist(subsamples["SR"]);
    setLogY(false);
    renewCanvas();
    TString drawopt="hist e";
    for (map<TString,TString>::iterator isubsample=subsamples.begin(); isubsample!=subsamples.end(); ++isubsample) {

      if (isubsample->first=="SR") continue;

      TH1D* numeratorhist = getHist( subsamples[ isubsample->first]);
      TString ratiohistname = numeratorhist->GetName();
      ratiohistname.Append("_ratio");

      TH1D* thisratio= (TH1D*) numeratorhist->Clone(ratiohistname);
      thisratio->Divide(SRhist);
      histos_[ratiohistname] = thisratio;
      thisratio->SetYTitle("ratio");
      thisratio->SetMaximum(2);
      thisratio->SetMinimum(0);
      thisratio->Draw(drawopt);
      if (!drawopt.Contains("same")) drawopt+=" same";
    }
    thecanvas->SaveAs("ewkhh_metsig-shapes2-"+thissample+"_ratios.eps");

    if (false) {
    //should look at MET too
    selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&mdp &&metVLoose;
    var="MET"; xtitle="MET";
    nbins=15; low= 0; high=300;
    setStackMode(false,false,false); //no stack
    setLogY(true);
    doRatio_=false; ratioMin = 0; ratioMax = 0.5;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_met-shapes2-"+thissample,0);
    
    setLogY(true);
    setStackMode(false,true,false); //norm
    doRatio_=false; ratioMin = 0; ratioMax = 2.2;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_met-shapes2-"+thissample,0);
    //}

    //now ratio plots
    SRhist = getHist(subsamples["SR"]);
    setLogY(false);
    renewCanvas();
    drawopt="hist e";
    for (map<TString,TString>::iterator isubsample=subsamples.begin(); isubsample!=subsamples.end(); ++isubsample) {

      if (isubsample->first=="SR") continue;

      TH1D* numeratorhist = getHist( subsamples[ isubsample->first]);
      TString ratiohistname = numeratorhist->GetName();
      ratiohistname.Append("_ratio");

      TH1D* thisratio= (TH1D*) numeratorhist->Clone(ratiohistname);
      thisratio->Divide(SRhist);
      histos_[ratiohistname] = thisratio;
      thisratio->SetYTitle("ratio");
      thisratio->SetMaximum(2);
      thisratio->SetMinimum(0);
      thisratio->Draw(drawopt);
      if (!drawopt.Contains("same")) drawopt+=" same";
    }
    thecanvas->SaveAs("ewkhh_met-shapes2-"+thissample+"_ratios.eps");
    }
    
    //now try MET with a METsig cut
    selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&mdp &&metsigVLoose;
    var="MET"; xtitle="MET";
    nbins=15; low= 0; high=300;
    setStackMode(false,false,false); //no stack
    setLogY(true);
    doRatio_=false; ratioMin = 0; ratioMax = 0.5;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_met-shapes2-metsig40-"+thissample,0);
    
    setLogY(true);
    setStackMode(false,true,false); //norm
    doRatio_=false; ratioMin = 0; ratioMax = 2.2;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_met-shapes2-metsig40"+thissample,0);
    //}

    //now ratio plots
    SRhist = getHist(subsamples["SR"]);
    setLogY(false);
    renewCanvas();
    drawopt="hist e";
    for (map<TString,TString>::iterator isubsample=subsamples.begin(); isubsample!=subsamples.end(); ++isubsample) {

      if (isubsample->first=="SR") continue;

      TH1D* numeratorhist = getHist( subsamples[ isubsample->first]);
      TString ratiohistname = numeratorhist->GetName();
      ratiohistname.Append("_ratio");

      TH1D* thisratio= (TH1D*) numeratorhist->Clone(ratiohistname);
      thisratio->Divide(SRhist);
      histos_[ratiohistname] = thisratio;
      thisratio->SetYTitle("ratio");
      thisratio->SetMaximum(2);
      thisratio->SetMinimum(0);
      thisratio->Draw(drawopt);
      if (!drawopt.Contains("same")) drawopt+=" same";
    }
    thecanvas->SaveAs("ewkhh_met-shapes2-metsig40-"+thissample+"_ratios.eps");
    

  return;


}

void higgsmbb_backgroundStudies1() {
  addToSamplesAll("QCD_Pt-1000to1400");
  addToSamplesAll("QCD_Pt-120to170");
  addToSamplesAll("QCD_Pt-1400to1800");
  addToSamplesAll("QCD_Pt-170to300-1");
  addToSamplesAll("QCD_Pt-170to300-2");
  addToSamplesAll("QCD_Pt-1800");
  addToSamplesAll("QCD_Pt-300to470-1");
  addToSamplesAll("QCD_Pt-300to470-2");
  addToSamplesAll("QCD_Pt-470to600");
  addToSamplesAll("QCD_Pt-600to800");
  addToSamplesAll("QCD_Pt-800to1000");

  addToSamplesAll("TTJets_FullLeptMGDecays");
  addToSamplesAll("TTJets_SemiLeptMGDecays");



  //use keith's cut flow.
  // Goal 1: check if METsig shape is independent of mbb

  nameOfEventWeight_="weight3";//don't want an event weight

  inputPath="/cu4/ra2b/reducedTrees/v68_4/";

  loadSamples(true,"ra2b2012");
  usePUweight_=true; //helps bring signal and background into alignment; not perfect but better than nothing

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 &&passCleaning==1"; 
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05||passMC_PFMET150";
  TCut zl = "nMuons==0&&nElectrons==0";
  TCut isotk="nIsoTracks15_005_03==0";

  TCut njets="njets20==4 || njets20==5";
  TCut btag = "CSVbest2>=0.898 && CSVbest3>=0.679 && CSVbest4>=0.244";

  TCut metsigVLoose = "METsig>=40";

  TCut mdp="minDeltaPhiN_asin>=4";
  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&mdp &&metsigVLoose;

  //this is an inelegant trick to manipulate the chains
  //it takes advantage of the fact that clearSamples() does not reset the chains.
  //it's a dirty hack but it works
  clearSamples();

  addSample("TTJets_FullLeptMGDecays");
  chainSamples("TTJets_FullLeptMGDecays","TTJets_SemiLeptMGDecays");

  addSample("QCD_Pt-1000to1400");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-120to170");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-1400to1800");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-170to300-1");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-170to300-2");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-1800");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-300to470-1");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-300to470-2");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-470to600");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-600to800");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-800to1000");

  for (int iii=0;iii<1; iii++) {
    
    clearSamples();
    
    TString thissample="";
    if (iii==0){
      addSample("TTJets_FullLeptMGDecays:higgsMbb2MassDiff>=100&&higgsMbb2MassDiff<=140&&higgsMbb1MassDiff>=100&&higgsMbb1MassDiff<=140",kGreen,"t#bar{t} (SR)");
      addSample("TTJets_FullLeptMGDecays:higgsMbb2MassDiff<100||higgsMbb2MassDiff>140||higgsMbb1MassDiff<100||higgsMbb1MassDiff>140",kBlue,"t#bar{t} (off SR)");
      thissample="ttbar";
    } else {
      addSample("QCD_Pt-1000to1400:higgsMbb2MassDiff>=100&&higgsMbb2MassDiff<=140&&higgsMbb1MassDiff>=100&&higgsMbb1MassDiff<=140",kGreen,"QCD (SR)");
      addSample("QCD_Pt-1000to1400:higgsMbb2MassDiff<100||higgsMbb2MassDiff>140||higgsMbb1MassDiff<100||higgsMbb1MassDiff>140",kBlue,"QCD (off SR)");
      thissample="qcd";
    }
    //////
    
    //  addSample("TTbarJets0:higgsMbb2MassDiff>=100&&higgsMbb2MassDiff<=140&&higgsMbb1MassDiff>=100&&higgsMbb1MassDiff<=140",kGreen,"t#bar{t} (SR)");
    //  addSample("TTbarJets0:higgsMbb2MassDiff<100||higgsMbb2MassDiff>140||higgsMbb1MassDiff<100||higgsMbb1MassDiff>140",kBlue,"t#bar{t} (off SR)");
    
    fitRatio_="pol1";
    
    var="METsig"; xtitle="MET significance";
    nbins=30; low= 0; high=300;
    setStackMode(false,false,false); //no stack
    setLogY(true);
    doRatio_=true; ratioMin = 0; ratioMax = 0.5;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_metsig-shapes1-"+thissample,0);
    
    setLogY(true);
    setStackMode(false,true,false); //norm
    doRatio_=true; ratioMin = 0; ratioMax = 2.2;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_metsig-shapes1-"+thissample,0);
    
    
    //should look at MET too
    TCut metVLoose = "MET>40";
    selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&mdp &&metVLoose;
    var="MET"; xtitle="MET";
    nbins=30; low= 0; high=300;
    setStackMode(false,false,false); //no stack
    setLogY(true);
    doRatio_=true; ratioMin = 0; ratioMax = 0.5;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_met-shapes1-"+thissample,0);
    
    setLogY(true);
    setStackMode(false,true,false); //norm
    doRatio_=true; ratioMin = 0; ratioMax = 2.2;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_met-shapes1-"+thissample,0);
  }
  return;

  //better SB definition?
  //what does the ttbar shape look like?

  //
  clearSamples();
  addSample("TTbarJets0");
  var="higgsMbb1MassDiff"; xtitle="mass 1";
  nbins=30; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  setLogY(false);
  doRatio_=false; ratioMin = 0; ratioMax = 1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_mbb1",0);

  var="higgsMbb2MassDiff"; xtitle="mass 2";
  nbins=30; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  setLogY(false);
  doRatio_=false; ratioMin = 0; ratioMax = 1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_mbb2",0);

}

void higgsmbb_cutflowDtry2_illustration() {

  //plots for my cornell talk


 nameOfEventWeight_="weight3";//don't want an event weight

  inputPath="/cu4/ra2b/reducedTrees/v68_4/";

  addToSamplesAll("TChihh250");
  addToSamplesAll("TChihh400");

  loadSamples(true,"ra2b2012");
  usePUweight_=true; //helps bring signal and background into alignment; not perfect but better than nothing

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  clearSamples();
  addSample("TTbarJets0",kGreen,"t#bar{t}");
  addSample("TChihh250",kRed,"hh 250");
  addSample("TChihh400",kBlue,"hh 400");

  setSampleScaleFactor("TChihh250",0.09/9999.0); //sigma x BF / ngen
  setSampleScaleFactor("TChihh400",0.01/9999.0); //sigma x BF / ngen

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 &&passCleaning==1"; 
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05||passMC_PFMET150";
  TCut zl = "nMuons==0&&nElectrons==0";
  TCut isotk="nIsoTracks15_005_03==0";

  TCut njets="njets20==4 || njets20==5";
  TCut btag = "CSVbest2>=0.898 && CSVbest3>=0.679 && CSVbest4>=0.244";


  TCut hbbhbb = "higgsMbb2MassDiff>=100&&higgsMbb2MassDiff<=140&&higgsMbb1MassDiff>=100&&higgsMbb1MassDiff<=140";

  TCut mdp="minDeltaPhiN_asin>=4";
  selection_=baseline && triggers &&zl &&isotk;

  var="njets20"; xtitle="jet multiplicity (pT>20 GeV)";
  nbins=10; low= 0; high=10;
  setStackMode(false,true,false); //no stack
  setPlotMinimum(0);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDb0_"+var,0);

  selection_=baseline && triggers &&zl &&isotk &&njets;


  var="CSVbest2"; xtitle=var;
  nbins=10; low= 0; high=1;
  setStackMode(false,true,false); //no stack
  setPlotMinimum(0);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDb1_"+var,0);
  var="CSVbest3"; xtitle=var;
  nbins=10; low= 0; high=1;
  setStackMode(false,true,false); //no stack
  setPlotMinimum(0);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDb1_"+var,0);
  var="CSVbest4"; xtitle=var;
  nbins=10; low= 0; high=1;
  setStackMode(false,true,false); //no stack
  setPlotMinimum(0);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDb1_"+var,0);

  selection_=baseline && triggers &&zl &&isotk &&njets&&btag;
  var="higgsMbb1MassDiff"; xtitle="higgs mass 1";
  nbins=20; low= 0; high=200;
  setStackMode(false,false,false); //no stack
  setPlotMinimum(0);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDb2_"+var,0);

  selection_=baseline && triggers &&zl &&isotk &&njets&&btag;
  var="higgsMbb2MassDiff"; xtitle="higgs mass 2";
  nbins=20; low= 0; high=200;
  setStackMode(false,false,false); //no stack
  setPlotMinimum(0);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDb2_"+var,0);
  setStackMode(false,true,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDb2_"+var,0);

  selection_=baseline && triggers &&zl &&isotk &&njets&&btag &&hbbhbb;

  var="MET"; xtitle="MET";
  nbins=15; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  setPlotMinimum(0.1);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDb3_MET",0);
  setLogY(false);
  drawEffVRej("TChihh400","TTbarJets0","hh400","ttbar",true);
  drawEffVRej("TChihh250","TTbarJets0","hh250","ttbar",true);


  var="METsig"; xtitle="MET significance";
  nbins=15; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  setPlotMinimum(0.1);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDb3_METsig",0);
  setLogY(false);
  drawEffVRej("TChihh400","TTbarJets0","hh400","ttbar",true);
  drawEffVRej("TChihh250","TTbarJets0","hh250","ttbar",true);

}

void higgsmbb_cutflowDtry2() {

  //same as cutflow below, but really follow keith exactly

/*
22 May -- rerun this. Only change is to add QCD and large ttbar.
*/


 nameOfEventWeight_="weight3";//don't want an event weight

  inputPath="/cu4/ra2b/reducedTrees/v68_4/";

  addToSamplesAll("TChihh250");
  addToSamplesAll("TChihh400");

  addToSamplesAll("QCD_Pt-1000to1400");
  addToSamplesAll("QCD_Pt-120to170");
  addToSamplesAll("QCD_Pt-1400to1800");
  addToSamplesAll("QCD_Pt-170to300-1");
  addToSamplesAll("QCD_Pt-170to300-2");
  addToSamplesAll("QCD_Pt-1800");
  addToSamplesAll("QCD_Pt-300to470-1");
  addToSamplesAll("QCD_Pt-300to470-2");
  addToSamplesAll("QCD_Pt-470to600");
  addToSamplesAll("QCD_Pt-600to800");
  addToSamplesAll("QCD_Pt-800to1000");

  addToSamplesAll("TTJets_FullLeptMGDecays");
  addToSamplesAll("TTJets_SemiLeptMGDecays");



  loadSamples(true,"ra2b2012");
  usePUweight_=true; //helps bring signal and background into alignment; not perfect but better than nothing

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  clearSamples();
  addSample("TTJets_FullLeptMGDecays",kGreen,"t#bar{t}");
  chainSamples("TTJets_FullLeptMGDecays","TTJets_SemiLeptMGDecays");

  addSample("QCD_Pt-1000to1400",kMagenta,"QCD");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-120to170");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-1400to1800");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-170to300-1");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-170to300-2");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-1800");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-300to470-1");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-300to470-2");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-470to600");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-600to800");
  chainSamples("QCD_Pt-1000to1400","QCD_Pt-800to1000");

  addSample("TChihh250",kRed,"hh 250");
  addSample("TChihh400",kBlue,"hh 400");

  setSampleScaleFactor("TChihh250",0.09/9999.0); //sigma x BF / ngen
  setSampleScaleFactor("TChihh400",0.01/9999.0); //sigma x BF / ngen

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 &&passCleaning==1"; 
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05||passMC_PFMET150";
  TCut zl = "nMuons==0&&nElectrons==0";
  TCut isotk="nIsoTracks15_005_03==0";

  TCut njets="njets20==4 || njets20==5";
  TCut btag = "CSVbest2>=0.898 && CSVbest3>=0.679 && CSVbest4>=0.244";

  TCut hbbhbb = "higgsMbb2MassDiff>=100&&higgsMbb2MassDiff<=140&&higgsMbb1MassDiff>=100&&higgsMbb1MassDiff<=140";

  TCut mdp="minDeltaPhiN_asin>=4";
  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&hbbhbb&&mdp;

  setStackMode(false,false,false); //no stack

  var="METsig"; xtitle="MET significance";
  nbins=30; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  setPlotMinimum(0.1);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDa0_metsig",0);
//   setLogY(false);
//   setStackMode(false,true,false); //norm
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDa0_metsig",0);

  TCut metsigLoose = "METsig>75";
  TCut metsigTight = "METsig>150";
  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&hbbhbb&&mdp &&metsigLoose;

  //very nice!
  var="deltaRmax_hh"; xtitle="deltaR max h(b,b)";
  nbins=10; low= 0; high=7;
  setLogY(false);
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDa1_dRmax-metsig75",0);

  var="deltaPhi_hh"; xtitle="deltaPhi(h,h)";
  nbins=10; low= 0; high=5;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDa1_dphh-metsig75",0);

  var="sumPt_hh-higgsMbb2MassDiff-higgsMbb1MassDiff"; xtitle="sumPt(b1..4)-mbb1-mbb2";
  nbins=10; low= -300; high=500;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDa1_sumptdiff-metsig75",0);

  var="sumPt_hh"; xtitle="sumPt(b1..4)";
  nbins=10; low= 0; high=800;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDa1_sumpt1234-metsig75",0);

  var="deltaRmin_hh"; xtitle="deltaR min h(b,b)";
  nbins=10; low= 0; high=6;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDa1_dRmin-metsig75",0);

  //i want to see if minDeltaPhiN is sensible or not. Need qcd samples for that though
   selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&hbbhbb &&metsigLoose;
  var="minDeltaPhiN_asin"; xtitle="minDeltaPhiN_asin";
  nbins=10; low= 0; high=40;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDa1_mindpN-metsig75",0);

  var="minDeltaPhi"; xtitle="minDeltaPhi";
  nbins=10; low= 0; high=4;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDa1_mindp-metsig75",0);

  var="deltaPhi1"; xtitle="deltaPhi1";
  nbins=10; low= 0; high=4;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDa1_deltaPhi1-metsig75",0);

  var="deltaPhi2"; xtitle="deltaPhi2";
  nbins=10; low= 0; high=4;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDa1_deltaPhi2-metsig75",0);

  //update (22 may) to get rid of mdp cut; not needed even for QCD (I think)
  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&hbbhbb &&metsigTight;
  //not useful in this case
  var="deltaRmax_hh"; xtitle="deltaR max h(b,b)";
  nbins=10; low= 0; high=7;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDa1_dRmax-metsig150",0);

  var="deltaPhi_hh"; xtitle="deltaPhi(h,h)";
  nbins=10; low= 0; high=5;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDa1_dphh-metsig150",0);

  var="deltaRmin_hh"; xtitle="deltaR min h(b,b)";
  nbins=10; low= 0; high=6;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDa1_dRmin-metsig150",0);

  //a little test -- update to remove mdp cut
  TCut drmin = "deltaRmin_hh<2";
  TCut metsigMedium = "METsig>100";
  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&hbbhbb &&metsigMedium&&drmin;
  var="deltaRmax_hh"; xtitle="deltaR max h(b,b)";
  nbins=7; low= 0; high=7;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDa2_dRmax",0);

  TCut drmax = "deltaRmax_hh<2";
  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&hbbhbb &&metsigMedium&&drmin&&drmax;
  var="deltaPhi_hh"; xtitle="deltaPhi(h,h)";
  nbins=10; low= 0; high=5;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDa3_dphh",0);

}

void higgsmbb_cutflowD() {
/* from Keith, april 29
 Attached is the best I've managed to do so far. The cuts are
-Good PV
-Pass MET tail cleaning
-Pass trigger emulation (including PFMET150 in addition to the two CSV triggers)
-RA2b lepton and iso track vetos
-minDelPhiN>4
-Exactly 4 or 5 jets with pt>20
- >=2 CSVT, >=3 CSVM, >=4 CSVL
-Select best pairing of 4 jets with highest CSV values in event based on smallest mass difference (jets with pt>20)
-100 < m(bb) < 140 for both m(bb) pairs

The plot is of MET significance. I don't actually know how this variable is computed, but I think Josh was using it ("METsig" in the reducedTrees) and it does a nice job discriminating signal from background. The y axis is in events with the samples weighted by luminosity.

For the 250 GeV signal point I can get to S/B~1 with S~20 events cutting METsig>75 GeV.
For the 400 GeV signal point I can get to S/B~2 with S~4 events cutting METsig>150 GeV.
 */

 nameOfEventWeight_="weight3";//don't want an event weight

  inputPath="/cu4/ra2b/reducedTrees/v68_4/";

  addToSamplesAll("TChihh250");
  addToSamplesAll("TChihh400");

  loadSamples(true,"ra2b2012");
  usePUweight_=true; //helps bring signal and background into alignment; not perfect but better than nothing

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  clearSamples();
  addSample("TTbarJets0",kGreen,"t#bar{t}");
  addSample("TChihh250:higgsMbb1MassDiff_correct==0",kRed-10,"hh 250 0right");
  addSample("TChihh250:higgsMbb1MassDiff_correct==1",kRed-7,"hh 250 1right");
  addSample("TChihh250:higgsMbb1MassDiff_correct==2",kRed,"hh 250 2right");
  addSample("TChihh400:higgsMbb1MassDiff_correct==0",kBlue-10,"hh 400 0right");
  addSample("TChihh400:higgsMbb1MassDiff_correct==1",kBlue-7,"hh 400 1right");
  addSample("TChihh400:higgsMbb1MassDiff_correct==2",kBlue,"hh 400 2right");


  setSampleScaleFactor("TChihh250",0.09/9999.0); //sigma x BF / ngen
  setSampleScaleFactor("TChihh400",0.01/9999.0); //sigma x BF / ngen

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1"; //leave out cleaning etc for now
  setStackMode(false,false,false); //no stack
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05||passMC_PFMET150";
  TCut zl = "nMuons==0&&nElectrons==0";
  TCut isotk="nIsoTracks15_005_03==0";

  TCut njets="njets20==4 || njets20==5";
  TCut btag = "CSVbest2>=0.898 && CSVbest3>=0.679 && CSVbest4>=0.244";

  TCut metsig = "METsig>75";

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&metsig;
  //no minDeltaPhiN for now

  var="higgsMbb1MassDiff"; xtitle="mbb (1)";
  nbins=20; low= 0; high=400;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_mbb1",0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_mbb1",0);

  var="higgsMbb2MassDiff"; xtitle="mbb (2)";
  nbins=20; low= 0; high=400;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_mbb2",0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_mbb2",0);

  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<mbb>";
  nbins=20; low= 0; high=400;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_mbb",0);
 

  //mass difference
  var="abs(higgsMbb1MassDiff-higgsMbb2MassDiff)"; xtitle="|mbb1-mbb2|";
  nbins=10; low= 0; high=100;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_mbbdiff",0);
 

  var="deltaPhi_hh"; xtitle="deltaPhi(h,h)";
  nbins=10; low= 0; high=5;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_dphh",0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_dphh",0);

  var="deltaRmax_hh"; xtitle="deltaR max (b,b)";
  nbins=10; low= 0; high=7;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_dRmax",0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_dRmax",0);

  var="deltaRmin_hh"; xtitle="deltaR min (b,b)";
  nbins=10; low= 0; high=7;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_dRmin",0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_dRmin",0);

  var="sumPt_hh"; xtitle="sumPt(b1..4)";
  nbins=10; low= 0; high=600;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_sumpthh",0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_sumpthh",0);


  var="sumPt_hh-higgsMbb2MassDiff-higgsMbb1MassDiff"; xtitle="sumPt(b1..4)-mbb1-mbb2";
  nbins=10; low= -250; high=400;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_sumptdiff",0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_sumptdiff",0);

  var="3.14159-deltaPhi_hh+deltaRmin_hh+deltaRmax_hh"; xtitle="poor mva 1";
  nbins=10; low= 0; high=12;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_mva1",0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_mva1",0);

  var="3.14159-deltaPhi_hh+deltaRmin_hh"; xtitle="poor mva 2";
  nbins=10; low= 0; high=10;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_mva2",0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_mva2",0);

  //maybe the key is to look at these variables for events that are outside the mass window?

  clearSamples();
  addSample("TChihh250:higgsMbb1MassDiff_correct==0",kRed-10,"hh 250 0right");
  addSample("TChihh250:higgsMbb1MassDiff_correct==1",kRed-7,"hh 250 1right");
  addSample("TChihh250:higgsMbb1MassDiff_correct==2",kRed,"hh 250 2right");
  addSample("TChihh400:higgsMbb1MassDiff_correct==0",kBlue-10,"hh 400 0right");
  addSample("TChihh400:higgsMbb1MassDiff_correct==1",kBlue-7,"hh 400 1right");
  addSample("TChihh400:higgsMbb1MassDiff_correct==2",kBlue,"hh 400 2right");

  TCut badmbb = "higgsMbb2MassDiff<100||higgsMbb2MassDiff>140||higgsMbb1MassDiff<100||higgsMbb1MassDiff>140";
  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&metsig&&badmbb;

  var="deltaRmax_hh"; xtitle="deltaR max (b,b)";
  nbins=10; low= 0; high=7;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_dRmax_badmbb",0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_dRmax_badmbb",0);

  //also, i should look at how often we get it right when we require that it is possible to get it right

  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&metsig;
  clearSamples();
  //addSample("TTbarJets0",kAzure-3,"t#bar{t}");
  addSample("TChihh250:ncompleteHiggsReco20==0",kRed-10,"hh 250 0possible");
  addSample("TChihh250:ncompleteHiggsReco20==1",kRed-7,"hh 250 1possible");
  addSample("TChihh250:ncompleteHiggsReco20==2",kRed,"hh 250 2possible");
  addSample("TChihh400:ncompleteHiggsReco20==0",kBlue-10,"hh 400 0possible");
  addSample("TChihh400:ncompleteHiggsReco20==1",kBlue-7,"hh 400 1possible");
  addSample("TChihh400:ncompleteHiggsReco20==2",kBlue,"hh 400 2possible");

  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<mbb>";
  nbins=20; low= 0; high=400;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowD0_mbb-bypossible",0);

  //maybe there's no point in trying to improve the higgs reco?
  //we're in keith's 40 GeV window most of the time, even if we don't get the higgs pairing correct!
  

}


void higgsmbb_cutflowC() {
/* 26 Apr */

  nameOfEventWeight_="weight3";//don't want an event weight

  inputPath="/cu4/ra2b/reducedTrees/v68_3/";

  addToSamplesAll("TChihh250");
  addToSamplesAll("TChihh400");

  loadSamples(true,"ra2b2012");
  usePUweight_=true; //helps bring signal and background into alignment; not perfect but better than nothing

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  clearSamples();
  addSample("TTbarJets0",kAzure-3,"t#bar{t}");
  addSample("TChihh250",kRed,"Chi hh 250");
  addSample("TChihh400",kViolet,"Chi hh 400");

  //keith's numbers
  //150 --> 0.7 
  //200 --> 0.25
  //in fact keith has 0.09 for this one
  setSampleScaleFactor("TChihh250",0.09/9999.0); //sigma x BF / ngen
  setSampleScaleFactor("TChihh400",0.01/9999.0); //sigma x BF / ngen

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1"; //leave out cleaning etc for now
  setStackMode(false,false,false); //no stack
  TCut METtrigger = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";

  selection_=baseline && METtrigger;

  var="nMuons+nElectrons"; xtitle="n e+mu";
  nbins=5; low= 0; high=5;
  //  setStackMode(false,false,false); //no stack
  //  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB1_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC1_nLeptons",0);


  var="nIsoTracks15_005_03"; xtitle="n iso tracks";
  nbins=4; low= 0; high=4;
  //  setStackMode(false,false,false); //no stack
  //  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB1_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC1_nIsoTk",0);

  var="nTausVLoose"; xtitle="n VLoose Taus";
  nbins=4; low= 0; high=4;
  //  setStackMode(false,false,false); //no stack
  //  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB1_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC1_nTausVL",0);

  TCut zl = "nMuons==0&&nElectrons==0";
  selection_=baseline && METtrigger &&zl;

  var="nIsoTracks15_005_03"; xtitle="n iso tracks";
  nbins=4; low= 0; high=4;
  //  setStackMode(false,false,false); //no stack
  //  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB1_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC2_nIsoTk",0);

  var="nTausVLoose"; xtitle="n VLoose Taus";
  nbins=4; low= 0; high=4;
  //  setStackMode(false,false,false); //no stack
  //  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB1_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC2_nTausVL",0);

  TCut isotk="nIsoTracks15_005_03==0";
  selection_=baseline && METtrigger &&zl&&isotk;

  var="nTausVLoose"; xtitle="n VLoose Taus";
  nbins=4; low= 0; high=4;
  //  setStackMode(false,false,false); //no stack
  //  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB1_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC3_nTausVL",0);

  var="nTausLoose"; xtitle="n Loose Taus";
  nbins=3; low= 0; high=3;
  //  setStackMode(false,false,false); //no stack
  //  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB1_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC3_nTausL",0);

  TCut tauveto="nTausLoose==0";
  selection_=baseline && METtrigger &&zl&&isotk &&tauveto;

  var="nbjets20"; xtitle="CSV b tags (20 GeV)";
  nbins=6; low= 0; high=6;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC4_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC4_"+var,0);

  TCut threeb = "nbjets20>=3";
  selection_=baseline && METtrigger &&zl&&isotk &&tauveto &&threeb;

  var="mjj_closestB20"; xtitle="mbb (closest DR, 20 GeV)";
  nbins=30; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);

  //other observables

  var="mjj_closestB20"; xtitle="mbb (closest DR, 20 GeV)";
  nbins=30; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);

  clearSamples();
  addSample("TTbarJets0",kAzure-3,"t#bar{t}");
  addSample("TChihh400:mjj_closestB20_correct==1",kViolet,"hh 400 (correct)");
  addSample("TChihh400:mjj_closestB20_correct==0",kPink,"hh 400 (wrong)");

  var="mjj_closestB20"; xtitle="mbb (closest DR, 20 GeV)";
  nbins=30; low= 0; high=300;
//   setStackMode(false,false,false); //no stack
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_mjj_closestB20_correctOrNot",0);

  clearSamples();
  addSample("TTbarJets0",kAzure-3,"t#bar{t}");
  addSample("TChihh250:mjj_closestB20_correct==1",kViolet,"hh 250 (correct)");
  addSample("TChihh250:mjj_closestB20_correct==0",kPink,"hh 250 (wrong)");

  var="mjj_closestB20"; xtitle="mbb (closest DR, 20 GeV)";
  nbins=30; low= 0; high=300;
//   setStackMode(false,false,false); //no stack
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_mjj_closestB20_correctOrNot250",0);

  clearSamples();
  addSample("TTbarJets0",kAzure-3,"t#bar{t}");
  addSample("TChihh250",kRed,"Chi hh 250");
  addSample("TChihh400",kViolet,"Chi hh 400");

  var="deltaR_closestB20"; xtitle="DR(bb) (closest DR, 20 GeV)";
  nbins=20; low= 0; high=6;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);

  //not useful at all
  var="cosHel_closestB20"; xtitle="cosHel (closest DR, 20 GeV)";
  nbins=20; low= -1; high=1;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);

  //some playing around with h(bb)+h(bb) with =3 b tags in mind
  var="deltaPhi_hb20"; xtitle="DeltaPhi(h,b)";
  nbins=20; low= 0; high=5;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);

  var="sumPtMjjDiff_closestB20"; xtitle="Sum pT jj - mjj";
  nbins=20; low= -200; high=400;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);


  var="MT_bestCSV"; xtitle=var;
  nbins=20; low= 0; high=400;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);


  var="METsig"; xtitle=var;
  nbins=20; low= 0; high=400;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC5_"+var,0);

  TCut metsig = "METsig>60";
  selection_=baseline && METtrigger &&zl&&isotk &&tauveto &&threeb &&metsig;

  var="mjj_closestB20"; xtitle="mbb (closest DR, 20 GeV)";
  nbins=30; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC6_"+var,0);

  clearSamples();
  addSample("TTbarJets0",kAzure-3,"t#bar{t}");
  addSample("TChihh400:mjj_closestB20_correct==1",kViolet,"hh 400 (correct)");
  addSample("TChihh250:mjj_closestB20_correct==1",kPink,"hh 250 (correct)");

  var="mjj_closestB20"; xtitle="mbb (closest DR, 20 GeV)";
  nbins=30; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC6_corr_"+var,0);


  var="njets20"; xtitle=var;
  nbins=10; low= 0; high=10;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC6_corr_"+var,0);

  var="njets30"; xtitle=var;
  nbins=10; low= 0; high=10;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC6_corr_"+var,0);

  TCut njetscut = "njets30>=3 && njets30<=5";
  selection_=baseline && METtrigger &&zl&&isotk &&tauveto &&threeb &&metsig &&njetscut;

  var="mjj_closestB20"; xtitle="mbb (closest DR, 20 GeV)";
  nbins=30; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC7_corr_"+var,0);

  var="nbjetsCSVT"; xtitle=var;
  nbins=5; low= 0; high=5;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC7_corr_"+var,0);

  TCut csvt="nbjetsCSVT>=1";
  selection_=baseline && METtrigger &&zl&&isotk &&tauveto &&threeb &&metsig &&njetscut && csvt;

  var="mjj_closestB20"; xtitle="mbb (closest DR, 20 GeV)";
  nbins=30; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC8_corr_"+var,0);


  var="deltaPhi_hb20"; xtitle=var;
  nbins=20; low= 0; high=5;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC8_corr_"+var,0);
  setStackMode(false,true,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC8_corr_"+var,0);

  //b jet pT

  var="bjetpt1"; xtitle="lead b jet pT";
  nbins=30; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC8_corr_"+var,0);
  setStackMode(false,true,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC8_corr_"+var,0);

  var="bjetpt2"; xtitle="2nd b jet pT";
  nbins=30; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC8_corr_"+var,0);
  setStackMode(false,true,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC8_corr_"+var,0);

  TCut bpt1="bjetpt1>80"; 
  selection_=baseline && METtrigger &&zl&&isotk &&tauveto &&threeb &&metsig &&njetscut && csvt &&bpt1;

  var="bjetpt2"; xtitle="2nd b jet pT";
  nbins=30; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC9_corr_"+var,0);
  setStackMode(false,true,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC9_corr_"+var,0);

/*
TTbarJets0 =     934.439 +/- 26.8795
TChihh400:mjj_closestB20_correct==1 =    22.3715 +/- 1.15763
TChihh250:mjj_closestB20_correct==1 =    64.8097 +/- 6.19091

*/

  TCut bpt2="bjetpt2>50"; 
  selection_=baseline && METtrigger &&zl&&isotk &&tauveto &&threeb &&metsig &&njetscut && csvt &&bpt1&&bpt2;

  var="mjj_closestB20"; xtitle="mbb (closest DR, 20 GeV)";
  nbins=30; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowC10_corr_"+var,0);
/*
TTbarJets0 =     813.605 +/- 25.0547
TChihh400:mjj_closestB20_correct==1 =    21.9365 +/- 1.15229
TChihh250:mjj_closestB20_correct==1 =    60.2526 +/- 5.92588


*/

}

void higgsmbb_cutflowB() {

  /*
24 April
0th order by-hand cutflow -- 2nd try
  */
  nameOfEventWeight_="weight3";//don't want an event weight

  inputPath="/cu4/ra2b/reducedTrees/v68_3/";

  addToSamplesAll("TChihh250");
  addToSamplesAll("TChihh400");

  loadSamples(true,"ra2b2012");
  usePUweight_=true; //helps bring signal and background into alignment; not perfect but better than nothing

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  clearSamples();
  addSample("TTbarJets0",kAzure-3,"t#bar{t}");
  addSample("TChihh250",kRed,"Chi hh 250");
  addSample("TChihh400",kViolet,"Chi hh 400");

  //keith's numbers
  //150 --> 0.7 
  //200 --> 0.25
  //in fact keith has 0.09 for this one
  setSampleScaleFactor("TChihh250",0.09/9999.0); //sigma x BF / ngen
  setSampleScaleFactor("TChihh400",0.01/9999.0); //sigma x BF / ngen

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 && cutTrigger2==1 && buggyEvent==0";
  setStackMode(false,false,false); //no stack

  TCut METtrigger = "MET>125 && jetpt2>=50";
  TCut btag= "nbjets30>=2"; //loose for now
  TCut zl = "nMuons==0&&nElectrons==0&&nIsoTracks15_005_03==0";
  //throw the kitchen sink at it
  TCut jettrigger = "jetpt2>=100 && jetpt4>=50";
  TCut tauveto = "nTausVLoose==0";

  selection_=baseline && (METtrigger||jettrigger) && btag &&zl&&tauveto;

  var="njets30"; xtitle="njets 30";
  nbins=10; low= 0; high=10;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB1_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB1_"+var,0);

  TCut njet4 = "njets30==4";
  selection_=baseline && (METtrigger||jettrigger) && btag &&zl&&tauveto &&njet4;

 
  //let's plot mjj related stuff right away
  var="mjjdiff"; xtitle=var;
  nbins=20; low= 0; high=100;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB2_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB2_"+var,0);

  TCut mjjdiffcut="mjjdiff<50";
  selection_=baseline && (METtrigger||jettrigger) && btag &&zl&&tauveto &&njet4 &&mjjdiffcut;

  var="0.5*(mjj1+mjj2)"; xtitle=var;
  nbins=30; low= 0; high=200;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB3_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB3_"+var,0);

  TCut mjjcut="0.5*(mjj1+mjj2)<140";
  selection_=baseline && (METtrigger||jettrigger) && btag &&zl&&tauveto &&njet4 &&mjjdiffcut&&mjjcut;

 //MET
  var="MET"; xtitle="MET (GeV)";
  nbins=30; low= 100; high=400;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB4_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB4_"+var,0);

 //jet/MET angles
  setPlotMinimum(0);
  var="deltaPhib1"; xtitle="DeltaPhi(bjet1,MET)";
  nbins=30; low= 0; high=TMath::Pi();
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB4_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB4_"+var,0);

  var="MT_bestCSV"; xtitle=var;
  nbins=20; low= 0; high=400;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB4_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB4_"+var,0);

  TCut mtcut = "MT_bestCSV>200";
  selection_=baseline && (METtrigger||jettrigger) && btag &&zl&&tauveto &&njet4 &&mjjdiffcut&&mjjcut&&mtcut;

  var="nbjetsCSVT"; xtitle=var;
  nbins=5; low= 0; high=5;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB5_"+var,0);

  var="nbjetsCSVL"; xtitle=var;
  nbins=5; low= 0; high=5;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB5_"+var,0);

  var="nbjets30"; xtitle=var;
  nbins=5; low= 0; high=5;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB5_"+var,0);

  TCut CSVTcut = "nbjetsCSVT>=1";
  selection_=baseline && (METtrigger||jettrigger) && btag &&zl&&tauveto &&njet4 &&mjjdiffcut&&mjjcut&&mtcut&&CSVTcut;

  var="bjetpt1"; xtitle=var;
  nbins=20; low= 0; high=500;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB6_"+var,0);

  var="nbjets30"; xtitle=var;
  nbins=5; low= 0; high=5;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB6_"+var,0);

  btag="nbjets30>=3";
  selection_=baseline && (METtrigger||jettrigger) && btag &&zl&&tauveto &&njet4 &&mjjdiffcut&&mjjcut&&mtcut&&CSVTcut;

  var="mjj_closestB"; xtitle=var;
  nbins=30; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB7_"+var,0);

  TCut mjjclose="mjj_closestB>105 && mjj_closestB<145";
  selection_=baseline && (METtrigger||jettrigger) && btag &&zl&&tauveto &&njet4 &&mjjdiffcut&&mjjcut&&mtcut&&CSVTcut &&mjjclose;
  var="MET"; xtitle="MET (GeV)";
  nbins=20; low= 100; high=400;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB8_"+var,0);

/*
TTbarJets0 =     11.0822 +/- 2.90624
TChihh250 =      7.32115 +/- 2.11558
TChihh400 =      4.31845 +/- 0.492785

*/

  //relax njets30 cut
  njet4="njets30<=4";
  selection_=baseline && (METtrigger||jettrigger) && btag &&zl&&tauveto &&njet4 &&mjjdiffcut&&mjjcut&&mtcut&&CSVTcut &&mjjclose;
   var="MET"; xtitle="MET (GeV)";
  nbins=20; low= 100; high=400;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB9_"+var,0);

/*
TTbarJets0 =     22.1508 +/- 4.0162
TChihh250 =      11.9557 +/- 2.31865
TChihh400 =      7.81513 +/- 0.692707
*/
  TCut tightermet="MET>200";
  selection_=baseline && (METtrigger||jettrigger) && btag &&zl&&tauveto &&njet4 &&mjjdiffcut&&mjjcut&&mtcut&&CSVTcut &&mjjclose &&tightermet;
  var="HT30"; xtitle=var;
  nbins=20; low= 200; high=600;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowB10_"+var,0);
/*
TTbarJets0 =     7.00484 +/- 2.28059
TChihh250 =      3.90085 +/- 1.47193
TChihh400 =      6.04992 +/- 0.625537
total SM =       16.9556 +/- 2.7855

*/

//need to be more clever!

}


void higgsmbb_cutflow() {

  /*
23 April
0th order by-hand cutflow
  */

  nameOfEventWeight_="weight3";//don't want an event weight

  inputPath="/cu4/ra2b/reducedTrees/v68_3/";

  addToSamplesAll("TChihh250");
  addToSamplesAll("TChihh400");

  loadSamples(true,"ra2b2012");
  usePUweight_=true; //helps bring signal and background into alignment; not perfect but better than nothing

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  clearSamples();
  addSample("TTbarJets0",kAzure-3,"t#bar{t}");
  addSample("TChihh250",kRed,"Chi hh 250");
  addSample("TChihh400",kViolet,"Chi hh 400");

  //keith's numbers
  //150 --> 0.7 
  //200 --> 0.25
  //in fact keith has 0.09 for this one
  setSampleScaleFactor("TChihh250",0.09/9999.0); //sigma x BF / ngen
  setSampleScaleFactor("TChihh400",0.01/9999.0); //sigma x BF / ngen

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 && cutTrigger2==1 && buggyEvent==0";
  setStackMode(false,false,false); //no stack

  TCut METtrigger = "MET>125 && jetpt2>=50";
  TCut btag= "nbjets>=3";
  TCut zl = "nMuons==0&&nElectrons==0&&nIsoTracks15_005_03==0";

  //don't use this for now
  TCut jettrigger = "jetpt2>=100 && jetpt4>=50";

  //bjetpt1
  selection_=baseline && METtrigger && btag &&zl; //no cuts

  //MET
  var="MET"; xtitle="MET (GeV)";
  nbins=40; low= 0; high=400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow1_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow1_"+var,0);

  var="njets20"; xtitle="njets 20";
  nbins=10; low= 0; high=10;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow1_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow1_"+var,0);

  var="nTausVLoose"; xtitle="VLoose Taus";
  nbins=3; low= 0; high=3;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow1_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow1_"+var,0);

  TCut tauveto = "nTausVLoose==0";
  selection_=baseline && METtrigger && btag &&zl&&tauveto; //no cuts

  var="njets20"; xtitle="njets 20";
  nbins=10; low= 0; high=10;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow2_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow2_"+var,0);


  var="njets30"; xtitle="njets 30";
  nbins=10; low= 0; high=10;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow2_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow2_"+var,0);

  TCut njetuppercut = "njets30<=5";
  selection_=baseline && METtrigger && btag &&zl&&tauveto  && njetuppercut; //no cuts

  //jet/MET angles
  var="deltaPhib1"; xtitle="DeltaPhi(bjet1,MET)";
  nbins=30; low= 0; high=TMath::Pi();
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow3_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow3_"+var,0);

  var="MT_bestCSV"; xtitle=var;
  nbins=20; low= 0; high=400;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow3_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow3_"+var,0);

  TCut mtcut = "MT_bestCSV>200";
  selection_=baseline && METtrigger && btag &&zl&&tauveto  && njetuppercut&&mtcut; 


  var="mjj_closestB"; xtitle=var;
  nbins=30; low= 0; high=400;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow4_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow4_"+var,0);

 //MET
  var="MET"; xtitle="MET (GeV)";
  nbins=30; low= 100; high=400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow4_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow4_"+var,0);

  TCut tightmet = "MET>250";
  selection_=baseline && METtrigger && btag &&zl&&tauveto  && njetuppercut&&mtcut&&tightmet; 

  var="mjj_closestB"; xtitle=var;
  nbins=20; low= 0; high=400;
  setStackMode(false,false,false); //no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow5_"+var,0);
  setStackMode(false,true,false); //norm
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflow5_"+var,0);

  //maybe i need a different strategy -- focus on the part of the signal that
  //has ==4 jets and try to get both higgses?

}

void higgsmbb_firstlook() {

/*
22 April 
plot a bunch of signal+ttbar distributions, no cuts.
main goal is to view the sculpting of the invariant mass variables in background
*/

  initHiggsSamples(false,"hh150 hh200 hh250 hh400"); //update to new loading framework

  usePUweight_=true; //helps bring signal and background into alignment; not perfect but better than nothing

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40";
  TCut ht = "HT>=400";
  TCut met = "MET>=125";
  TCut mdp = "minDeltaPhiN_asin>=4";
  TCut ldp = "minDeltaPhiN_asin<4";
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");
  TCut zlwoisotk = getLeptonVetoCut();


  setStackMode(false,true,false); //normalized to unit area

  selection_=""; //no cuts

  //npv (sanity)

  nbins=25; low= 0; high=50;
  var="nGoodPV"; xtitle="n good PV";
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_npv",0);

  //njets20
  var="njets20"; xtitle="n jets pT>20 GeV";
  nbins=10; low= 0; high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  //njets30
  var="njets30"; xtitle="n jets pT>30 GeV";
  nbins=10; low= 0; high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  //nbjets30
  var="nbjets30"; xtitle="n CSVM jets pT>30 GeV";
  nbins=5; low= 0; high=5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  //jetpt1
  var="jetpt1"; xtitle="jet pt 1 (GeV)";
  nbins=40; low= 0; high=400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  //jetpt2
  var="jetpt2"; xtitle="jet pt 2 (GeV)";
  nbins=40; low= 0; high=400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "0ewkhh_"+var,0);

  //jetpt3
  var="jetpt3"; xtitle="jet pt 3 (GeV)";
  nbins=40; low= 0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  //jetpt4
  var="jetpt4"; xtitle="jet pt 4 (GeV)";
  nbins=40; low= 0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);
  //bjetpt1
  var="bjetpt1"; xtitle="b jet pt 1 (GeV)";
  nbins=40; low= 0; high=400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);


  //MET
  var="MET"; xtitle="MET (GeV)";
  nbins=40; low= 0; high=400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "0ewkhh_"+var,0);

  //HT
  var="HT"; xtitle="HT (GeV)";
  nbins=40; low= 0; high=800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  var="HT30"; xtitle="HT (30 GeV jets) (GeV)";
  nbins=40; low= 0; high=800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  //jet/MET angles
  setPlotMinimum(0);
  var="deltaPhi1"; xtitle="DeltaPhi(jet1,MET)";
  nbins=30; low= 0; high=TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  var="deltaPhi2"; xtitle="DeltaPhi(jet2,MET)";
  nbins=30; low= 0; high=TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  var="deltaPhi3"; xtitle="DeltaPhi(jet3,MET)";
  nbins=30; low= 0; high=TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  clearSamples();
  addSample("TTbarJets0:ttbarDecayCode==2",kOrange+10,"t#bar{t} (had)");
  addSample("TTbarJets0:ttbarDecayCode!=2",kAzure-3,"t#bar{t} (#geq 1l )");
  addSample("TChihh250",kRed,"Chi hh 250");
  addSample("TChihh400",kViolet,"Chi hh 400");

  //jet/MET angles
  var="deltaPhib1"; xtitle="DeltaPhi(bjet1,MET)";
  nbins=30; low= 0; high=TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  var="nElectrons+nMuons"; xtitle="n mu+e";
  nbins=3; low= 0; high=3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  //this wasn't invented for this search but let's look at it
  var="mjj1"; xtitle="mjj1 (pair produced jj resonance)";
  nbins=50; low= 0; high=350;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  //
  var="mjjb1"; xtitle="m_jjb";
  nbins=50; low= 0; high=350;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  //
  var="mjjb2"; xtitle="m_jjb 2";
  nbins=50; low= 0; high=350;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  var="rMET"; xtitle="MET/total energy";
  nbins=40; low= 0; high=1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);


  var="METsig"; xtitle="MET significance";
  nbins=40; low= 0; high=100;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);
  //MT variables

  var="MT_bestCSV"; xtitle=var;
  nbins=40; low= 0; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  var="minMT_jetMET"; xtitle=var;
  nbins=40; low= 0; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  //sphericity

  var="transverseSphericity_jets30"; xtitle=var;
  nbins=40; low= 0; high=4;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  //Delta R(bb)

  var="deltaR_bestTwoCSV"; xtitle=var;
  nbins=30; low= 0; high=5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  var="deltaPhi_bestTwoCSV"; xtitle=var;
  nbins=30; low= 0; high=TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);



  //invariant masses!
  var="mjj_bestTwoCSV"; xtitle=var;
  nbins=50; low= 0; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  //invariant masses!
  var="mjj_closestB"; xtitle=var;
  nbins=50; low= 0; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  var="mjj_h125"; xtitle=var;
  nbins=40; low= 0; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  selection_="njets20>=4";
  var="higgsMbb1"; xtitle=var;
  nbins=40; low= 0; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);
  var="higgsMbb2"; xtitle=var;
  nbins=40; low= 0; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  var="higgsMjj1"; xtitle=var;
  nbins=40; low= 0; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);
  var="higgsMjj2"; xtitle=var;
  nbins=40; low= 0; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);

  //EXO magic
  selection_="";
  var="higgsMbb1delta"; xtitle=var;
  nbins=40; low= 0; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);
  var="higgsMbb2delta"; xtitle=var;
  nbins=40; low= 0; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_"+var,0);



}
