/*
to compile:

must have symlink to MiscUtil.cxx in the working directory.

gSystem->Load("TSelectorMultiDraw_C.so");
gSystem->Load("CrossSectionTable_cxx.so");
gSystem->Load("ConfigurationDescriptions_cxx.so");
//gSystem->Load("SearchRegion_cxx.so"); //no longer needed
//gSystem->Load("SystInfo_cxx.so");     //no longer needed
//gSystem->Load("SignalEffData_cxx.so");//no longer needed
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
	
TString inputPath = "/cu4/ra2b/reducedTrees/v66_7/"; //2012
TString dataInputPath =  "/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/";
//TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-35u/Fall11/";//7 TeV

double lumiScale_ = 4301;//2012 lumi ????
double preLumiScale_ = 30;//god only knows for 2012

//make a symlink that point from this name to drawReducedTree.h
//this is to make the ROOT dictionary generation work correctly
#include "drawReducedTrees2012.h"

void dataABCD() {
  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  //must be done before loadSamples()
  addToSamplesAll("HTRun2012A");
  addToSamplesAll("HTMHTRun2012B");
  addToSamplesAll("HTMHTRun2012C1");
  addToSamplesAll("HTMHTRun2012C2");
  addToSamplesAll("HTMHTRun2012D");

  addToSamplesAll("JetHTRun2012B");
  addToSamplesAll("JetHTRun2012C1");
  addToSamplesAll("JetHTRun2012C2");
  addToSamplesAll("JetHTRun2012D");

  addToSamplesAll("METRun2012A");
  addToSamplesAll("METRun2012B");
  addToSamplesAll("METRun2012C1");
  addToSamplesAll("METRun2012C2");
  addToSamplesAll("METRun2012D");

  loadSamples(true,"ra2b2012");

  usePUweight_=false;
  dodata_=false;

  clearSamples();
  addSample("HTRun2012A",kBlue,"2012 A+B+C");
  chainSamples("HTRun2012A","HTMHTRun2012B");
  chainSamples("HTRun2012A","HTMHTRun2012C1");
  chainSamples("HTRun2012A","HTMHTRun2012C2");

  chainSamples("HTRun2012A","JetHTRun2012B");
  chainSamples("HTRun2012A","JetHTRun2012C1");
  chainSamples("HTRun2012A","JetHTRun2012C2");

  chainSamples("HTRun2012A","METRun2012A");
  chainSamples("HTRun2012A","METRun2012B");
  chainSamples("HTRun2012A","METRun2012C1");
  chainSamples("HTRun2012A","METRun2012C2");

  setSampleScaleFactor("HTRun2012A",5580./12034.14);

  addSample("METRun2012D",kRed,"2012 D");
  chainSamples("METRun2012D","HTMHTRun2012D");
  chainSamples("METRun2012D","JetHTRun2012D");

  drawFilenameOnPlot_=true;

  doOverflowAddition(true);
  setStackMode(true);//no stack

  lumiScale_ = 1;
  savePlots_=true;

  doRatio_=true; ratioMin=0; ratioMax=2;

  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40";
  TCut htloose = "HT>=400";
  TCut httight = "HT>=800";
  TCut metloose = "MET>=125";
  TCut metmedium = "MET>=200";
  TCut mdp = "minDeltaPhiN_asin>=4";
  TCut ldp = "minDeltaPhiN_asin<4";
  TCut btag = "nbjets>=1";
  TCut bveto = "nbjets==0";
  TCut sl = getSingleLeptonCut();
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");
  TCut withisotk = TCut("nIsoTracks15_005_03>=1");

  TString desc,var,xtitle;
  int nbins; float low; float high;
  setStackMode(false);

  TCut met = metmedium;
  TCut ht = htloose;
  desc = "_met200_ht400";
  TString sampledesc;
  TCut thisselection;
  for (int ii = 0; ii<2; ii++) {
    if (ii==0) {
      sampledesc = "_sl";
      thisselection = sl&&mdp;
    }
    else if (ii==1) {
      sampledesc="_ldp";
      thisselection = zl && ldp;
    }

  selection_=baseline && cleaning && ht && btag && thisselection && met;
  var="njets"; xtitle="jet multiplicity";
  nbins=10; low=0;high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ABCvD_njets"+sampledesc+desc,0);

  selection_=baseline && cleaning && ht  && btag && thisselection && met;
  var="nGoodPV"; xtitle="n good PV";
  nbins=25; low=0;high=50;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ABCvD_npv"+sampledesc+desc,0);

  selection_=baseline && cleaning && ht  && btag && thisselection && met;
  var = "jetchargedhadronmult1"; xtitle="lead jet ch had mult";
  nbins = 20; low = 0; high = 100;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ABCvD_chhadmult1"+sampledesc+desc,0);

  selection_=baseline && cleaning && ht  && btag && thisselection && met;
  var = "jetpt1"; xtitle="lead jet pT [GeV]";
  nbins = 30; low = 0; high = 1000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ABCvD_jetpt1"+sampledesc+desc,0);

  selection_=baseline && cleaning && ht  && btag && thisselection && met;
  var = "MET/rawPFMET"; xtitle="Type I PF MET / raw PF MET";
  nbins = 25; low = 0.6; high = 1.6;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ABCvD_METoverRawMET"+sampledesc+desc,0);

  selection_=baseline && cleaning && ht  && btag && thisselection ;
  var = "MET"; xtitle="MET";
  setLogY(true);
  nbins = 40; low = 100; high = 500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ABCvD_met"+sampledesc+desc,0);
  setLogY(false);

  selection_=baseline && cleaning  && btag && thisselection &&met ;
  var = "HT"; xtitle="HT";
  setLogY(true);
  nbins = 40; low = 400; high = 1400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ABCvD_ht"+sampledesc+desc,0);
  setLogY(false);
  }

}

void checkttwt() {
  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples
 loadSamples(true,"ra2b2012");

  usePUweight_=true;
  dodata_=false;

  clearSamples();
  //  addSample("TTbarJetsPowheg");
  addSample("TTbarJets0");
  //  addSample("WJets");
  //  addSample("SingleTop");

  doOverflowAddition(true);
  setStackMode(true);//no stack

  lumiScale_ = 12000;
  savePlots_=false;


  int nbins;
  float low,high;
  TString var,xtitle;
  
  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2";//&&maxTOBTECjetDeltaMult<40";

  //  TCut ht = "HT>=500 && HT<800";
  TCut ht = "HT>=400";
  //  TCut met = "MET>=350";
  TCut met = "MET>=125";
  TCut mdp = "minDeltaPhiN_asin>=4";
  //  TCut mdp = "minDeltaPhiN>=4";
  TCut ldp = "minDeltaPhiN_asin<4";

  TCut btag = "nbjets==1";

  //  TCut sl = getSingleLeptonCut();
  TCut zl = getLeptonVetoCut();//&&TCut("nIsoTracks15_005_03==0");
  //  TCut withisotk = TCut("nIsoTracks15_005_03>=1");


  nbins=5;
  low = 3; high=8;
  var="njets"; xtitle=var;

  selection_ =baseline&&cleaning&&ht&&met&&mdp&&btag&&zl;
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy1",0,"GeV");

  selection_ =baseline&&cleaning&&ht&&met&&mdp&&zl;
  btagSFweight_="prob1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy1",0,"GeV");
 


}

void v66_comparettbar() {
  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples
 loadSamples(true,"ra2b2012");

  usePUweight_=false;
  dodata_=false;

  clearSamples();
  addSample("TTbarJets");
  addSample("TTbarJets0");
  addSample("TTbarJetsPowheg");
  addSample("TTbarJetsMCNLO");
    
  doOverflowAddition(true);
  setStackMode(false);//no stack

  lumiScale_ = 10000;

  int nbins;
  float low,high;
  TString var,xtitle;
  

  nbins=30;
  low = 125; high=425;
  var="MET"; xtitle=var;
  selection_ =TCut("HT>=400 && cutPV==1 && njets>=3 && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1 && MET>=125")&&getLeptonVetoCut();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbarcomp_met",0,"GeV");
 
  nbins=10; low=0;high=10;
  var="njets"; xtitle=var;
  selection_ =TCut("HT>=400 && cutPV==1 && minDeltaPhiN >= 4 && passCleaning==1 && MET>=125")&&getLeptonVetoCut();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbarcomp_njets",0);
  setStackMode(false,true);//no stack
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbarcomp_njets",0);
  

}

void v65_chaincheck() {
  nameOfEventWeight_="weight2"; //for 2012 cfA skimmed ntuples
  usePUweight_=true;

 inputPath = "/cu4/ra2b/reducedTrees/v65_0/ORIGINALS/";
 dataInputPath =  "/cu4/ra2b/reducedTrees/v65_0/ORIGINALS/";

 loadSamples(true,"ra2b2012");
  currentConfig_=configDescriptions_.getDefault();

  doOverflowAddition(true);
  dodata_=true;
  setStackMode(true);//stack

 clearSamples();
 addSample("QCD1000",kYellow,"QCD");
 chainSamples("QCD1000","QCD120");
 chainSamples("QCD1000","QCD1400");
 chainSamples("QCD1000","QCD170");
 chainSamples("QCD1000","QCD1800");
 chainSamples("QCD1000","QCD300");
 chainSamples("QCD1000","QCD470");
 chainSamples("QCD1000","QCD600");
 chainSamples("QCD1000","QCD800");

 addSample("TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1403ra2b_v65",kAzure-3,"t#bar{t}");

 addSample("T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1421ra2b_v65",kMagenta,"Single Top");
 chainSamples("T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1421ra2b_v65","T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1470ra2b_v65");
 chainSamples("T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1421ra2b_v65","T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1455ra2b_v65");
 chainSamples("T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1421ra2b_v65","Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1422ra2b_v65");
 chainSamples("T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1421ra2b_v65","Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1461ra2b_v65");
 chainSamples("T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1421ra2b_v65","Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1451ra2b_v65");

 addSample("WJets_HT250To300",kGreen-3,"W+jets");
 chainSamples("WJets_HT250To300","WJets_HT300To400");
 chainSamples("WJets_HT250To300","WJets_HT400ToInf");

 addSample("WW",kCyan+1,"VV");
 chainSamples("WW","WZ");
 chainSamples("WW","ZZ");

 addSample("ZJets_HT200To400",kViolet-3,"Z+jets");
 chainSamples("ZJets_HT200To400","ZJets_HT400ToInf");

 addSample("Zinvisible_HT100To200",kOrange-3,"Z #to #nu#nu");
 chainSamples("Zinvisible_HT100To200","Zinvisible_HT200To400");
 chainSamples("Zinvisible_HT100To200","Zinvisible_HT400ToInf");

 setDatasetToDraw("2012hybrid");  lumiScale_=12034.14; //Run 2012 A+B+C ?

 int nbins;
 float low,high;
 TString var,xtitle;
 
 nbins=10;
 low = 0; high=10;
 var="njets"; xtitle=var;
 
 selection_="HT>=400 && MET>=150 && minDeltaPhiN>=4 && cutPV==1 &&cutTrigger==1 &&passCleaning==1";
 drawPlots(var,nbins,low,high,xtitle,"Events", "chaintest_v65",0);
  

}

void chaincheck() {
  
  nameOfEventWeight_="weight2"; //for 2012 cfA skimmed ntuples

  loadSamples(true,"ra2b2012");
  usePUweight_=false;

  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getDefault();

  int nbins;
  float low,high;
  TString var,xtitle;
  
  doOverflowAddition(true);
  dodata_=false;
  setStackMode(true);//stack

  clearSamples();
  addSample("TTbarJets");
  addSample("WJets");
  addSample("SingleTop");

  nbins=10;
  low = 0; high=10;
  var="njets"; xtitle=var;

  selection_="HT>=400 && MET>=150 && minDeltaPhiN>=4 && cutPV==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "chaintest_ttwt_normal",0,"GeV");

  clearSamples();
  addSample("TTbarSingleTopWJetsCombined");
  drawPlots(var,nbins,low,high,xtitle,"Events", "chaintest_ttwt_hadd",0,"GeV");

  clearSamples();
  addSample("TTbarJets");
  chainSamples("TTbarJets","WJets");
  chainSamples("TTbarJets","SingleTop");
  drawPlots(var,nbins,low,high,xtitle,"Events", "chaintest_ttwt_chain",0,"GeV");

  resetChains();
  clearSamples();
  addSample("TTbarJets");
  drawPlots(var,nbins,low,high,xtitle,"Events", "chaintest_ttwt_normal_check",0,"GeV");


}

void ra2b2011checks() {

  inputPath = "/cu2/ra2b/reducedTrees/V00-02-35u/Fall11/"; //2011 signal scans
  lumiScale_=4982; //i think
  loadSamples(true,"");

  usePUweight_=true;

  dodata_=false;

  int nbins;
  float low,high;
  TString var,xtitle;


  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getCorrected();

  setStackMode(false,true);
  doRatio_=false;

  clearSamples();
  addSample("TTbarJets:((W1decayType==13 && W2decayType==112) || (W1decayType==112 && W2decayType==13) || (W1decayType==1513 && W2decayType==112) || (W1decayType==112 && W2decayType==1513) || (W1decayType==13 && W2decayType==134) || (W1decayType==134 && W2decayType==13) || (W1decayType==1513 && W2decayType==134) || (W1decayType==134 && W2decayType==1513))",kRed,"#mu");
  addSample("TTbarJets:((W1decayType==11 && W2decayType==112) || (W1decayType==112 && W2decayType==11) || (W1decayType==1511 && W2decayType==112) || (W1decayType==112 && W2decayType==1511) || (W1decayType==11 && W2decayType==134) || (W1decayType==134 && W2decayType==11) || (W1decayType==1511 && W2decayType==134) || (W1decayType==134 && W2decayType==1511))",kBlue,"e");
  addSample("TTbarJets:((W1decayType==15 && W2decayType==112) || (W1decayType==112 && W2decayType==15) || (W1decayType==15 && W2decayType==134) || (W1decayType==134 && W2decayType==15))",kMagenta,"#tau #rightarrow h");

  selection_ =TCut("HT>=400 && cutPV==1 && njets>=3 && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1 && MET>=150")&&getLeptonVetoCut();
  nbins=45; low=150; high=600;
  var="MET"; xtitle="MET";
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b2011_metshapes1",0,"GeV");

  clearSamples();
  addSample("TTbarJets:nbjets==1",kRed,"1");
  addSample("TTbarJets:nbjets==2",kGreen,"2");
  addSample("TTbarJets:nbjets>=3",kBlue,"3");
  selection_ =TCut("HT>=400 && cutPV==1 && njets>=3 && minDeltaPhiN >= 4 && passCleaning==1 && MET>=150")&&getLeptonVetoCut();
  nbins=45; low=150; high=600;
  var="MET"; xtitle="MET";
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b2011_metshapes2",0,"GeV");

  nbins=40; low=400; high=1200;
  var="HT"; xtitle="HT";
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b2011_htshapes2",0,"GeV");

}


void RA2bJESshape2012() {
  //2012 concept...but use 2011 samples

  inputPath = "/cu4/ra2b/reducedTrees/v66_0/"; //2011 signal scans
  lumiScale_=12034.14; //Run 2012 A+B+C 
  loadSamples(true,"scan");

  usePUweight_=true;


  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getDefault();

  clearSamples();
  addSample("T1bbbb$900$600$");

  dodata_=false;

  int nbins;
  float low,high;
  TString var,xtitle;
  
  //unfortunately the tools don't lend themselves to auto-plotting these as separate samples.
  //hack it for now

  currentConfig_=configDescriptions_.getCorrected();
  selection_ =TCut("HT>=400 && cutPV==1 && njets>=3 && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1 && MET>=150 && jetpt2>=70")&&getLeptonVetoCut();
  nbins=6; low=150; high=400;
  var="MET"; xtitle="MET";
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb8TeV_met_jes0",0,"GeV");
  TH1D* hjes0 = (TH1D*)  hinteractive->Clone("hjes0");

  currentConfig_=configDescriptions_.at(2);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb8TeV_met_jesdown",0,"GeV");
  TH1D* hjesdown = (TH1D*)  hinteractive->Clone("hjesdown");
  
  currentConfig_=configDescriptions_.at(3);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb8TeV_met_jesup",0,"GeV");
  TH1D* hjesup = (TH1D*)  hinteractive->Clone("hjesup");
  
  renewCanvas();
  hjes0->Draw("hist e");
  hjesdown->Draw("hist e SAME");
  hjesup->Draw("hist e SAME");

  hjesdown->SetLineColor(kRed);
  hjesup->SetLineColor(kBlue);

  renewLegend();
  leg->AddEntry(hjes0,"nominal");
  leg->AddEntry(hjesdown," JES-");
  leg->AddEntry(hjesup," JES+");
  leg->Draw();

}

void RA2b_sl_byMt() {

  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples
  loadSamples(true,"ra2b2012");
  usePUweight_=true;

  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getDefault();

  int nbins;
  float low,high;
  TString var,xtitle;
  
  doOverflowAddition(true);

  dodata_=false;

  setStackMode(false,true);//normalized

  clearSamples();
  selection_ =TCut("MET>=125 && cutPV==1 && njets>=3  && minDeltaPhiN_asin >= 4 && passCleaning==1 && nbjets>=1 && HT>=400");
  addSample("TTbarJets:(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)",kBlue,"MT<100 GeV");
  addSample("TTbarJets:(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=100)",kRed,"MT>100 GeV");

  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 30; low=125; high=425;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_metshapes_bymt",0,"GeV");

  selection_ =TCut("MET>=125 && cutPV==1 && njets>=3  && minDeltaPhiN_asin >= 4 && passCleaning==1  && HT>=400 && jetpt2>=70");
  var="nbjets"; xtitle=var;
  nbins = 4; low=0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_bshapes_bymt",0,"GeV");

  selection_ =TCut("MET>=125 && cutPV==1 && njets>=3  && minDeltaPhiN_asin >= 4 && passCleaning==1  && HT>=500 &&HT<800&& jetpt2>=70 &&nbjets>=1");
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 30; low=125; high=425;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_metshapes_ht500_bymt",0,"GeV");


  //nIsoTracks15_005_03_lepcleaned
  setStackMode(true);
  clearSamples();
  selection_ =TCut("MET>=125 && cutPV==1 && njets>=3&&jetpt2>=70  && minDeltaPhiN_asin >= 4 && passCleaning==1 && nbjets>=1 && HT>=400 && (nElectrons+nMuons)==1 && MET/caloMET<2");
  addSample("TTbarJets:ttbarDecayCode==7",kRed,"#tau #rightarrow h");
  //addSample("TTbarJets:nElectrons==0&&nMuons==0&&ttbarDecayCode==2",kYellow,"all had");
  addSample("TTbarJets:ttbarDecayCode==1||ttbarDecayCode==9||ttbarDecayCode==8",kOrange,"dilepton (incl #tau)");
  addSample("TTbarJets:ttbarDecayCode==3||ttbarDecayCode==4||ttbarDecayCode==5||ttbarDecayCode==6",kMagenta,"e / #mu (incl #tau)");
  //  addSample("TTbarJets:ttbarDecayCode==5||ttbarDecayCode==6",kBlue,"#tau #rightarrow e / #mu");
  var = "nIsoTracks15_005_03_lepcleaned"; xtitle ="number of isolated tracks after lepton";
  nbins = 3; low=0; high=3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_isotracksl",0,"GeV");


}

void RA2b_ttbarstudies() {

  //goal: compare HT (and MET and nb) shapes for different 0-lepton sample versus 1 lepton

  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  loadSamples(true,"ra2b2012");
  usePUweight_=true;

  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getDefault();

  int nbins;
  float low,high;
  TString var,xtitle;
  
  doOverflowAddition(true);

  dodata_=false;
  doRatioPlot(true);

  setStackMode(false,true);//normalized

  clearSamples();
  addSample("TTbarJets:nElectrons==0&&nMuons==0&&nIsoTracks15_005_03==0",kRed,"0 lep");
  addSample("TTbarJets:(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)",kBlue,"1 lep");
   
  chainSamples("TTbarJets","WJets");
  chainSamples("TTbarJets","SingleTop");

  selection_ =TCut("HT>=400 &&MET>=125&& cutPV==1 && njets>=3 &&jetpt2>=70 && minDeltaPhiN_asin >= 4 && passCleaning==1 && nbjets>=1 &&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40");
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 30; low=125; high=425;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_metshapes_0l_1l",0,"GeV");

  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 40; low=400; high=1200;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_htshapes_0l_1l",0,"GeV");

  //njets
  var="nbjets-0.5"; xtitle="n b jets (p_{T}>50 GeV)";
  nbins = 3; low=0.5; high=3.5;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_nbjetsshapes_0l_1l",0);

  //try ST, out of curiosity

  selection_ =TCut("MET>=150 && cutPV==1 && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1");
  var="ST"; xtitle="S_{T} [GeV]";
  nbins = 35; low=300; high=1000;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_stshapes_0l_1l",0,"GeV");
  //ok, so that doesn't work

  //try meff
  selection_ =TCut("MET>=125 && HT>=400 && cutPV==1 && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1");
  var="HT+MET"; xtitle="m_{eff} [GeV]";
  nbins = 35; low=525; high=2000;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_meffshapes_0l_1l",0,"GeV");
  //not fantastically terrible, but worse than MET and HT alone

  clearSamples();
  addSample("TTbarJets:nElectrons==0&&nMuons==0&&ttbarDecayCode==7",kRed,"#tau #rightarrow h");
  addSample("TTbarJets:nElectrons==0&&nMuons==0&&(ttbarDecayCode==3||ttbarDecayCode==4)",kMagenta,"lost lepton");
  addSample("TTbarJets:(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)",kBlue,"1 lep");

  selection_ =TCut("MET>=150 && cutPV==1 && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1");
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 20; low=300; high=1000;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_htshapes_0l_1l_split",0,"GeV");

  //njets
  var="njets"; xtitle="n jets (p_{T}>50 GeV)";
  nbins = 10; low=0; high=10;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_njetsshapes_0l_1l_split",0);

  //jet pT
  var="jetpt1"; xtitle="lead jet pt";
  nbins = 20; low=0; high=500;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_jetptshapes_0l_1l_split",0);

  //"pseudo-HT"
  var="jetpt1+jetpt2+jetpt3"; xtitle="pT sum of lead 3 jets";
  nbins = 20; low=100; high=900;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_leadjetptsumshapes_0l_1l_split",0);

  //make sure it isn't the MT cut

  clearSamples();
  addSample("TTbarJets:(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)",kBlue-5,"1 lep");
  // addSample("TTbarJets:((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0))",kBlue,"1 lep (all MT)");
  addSample("TTbarJets:((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0))&& MT_Wlep>=100",kAzure,"1 lep (high MT)");

  selection_ =TCut("MET>=150 && cutPV==1 && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1");
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 20; low=300; high=1000;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_htshapes_1l_1lhighMT",0,"GeV");
  //inconclusive! no smoking gun

  //split HT in 0-lep by all categories

  clearSamples();
  addSample("TTbarJets:nElectrons==0&&nMuons==0&&ttbarDecayCode==7",kRed,"#tau #rightarrow h");
  //addSample("TTbarJets:nElectrons==0&&nMuons==0&&ttbarDecayCode==2",kYellow,"all had");
  addSample("TTbarJets:nElectrons==0&&nMuons==0&&ttbarDecayCode==1",kOrange,"di e / #mu");
  addSample("TTbarJets:nElectrons==0&&nMuons==0&&(ttbarDecayCode==3||ttbarDecayCode==4)",kMagenta,"e / #mu");
  addSample("TTbarJets:nElectrons==0&&nMuons==0&&(ttbarDecayCode==5||ttbarDecayCode==6)",kBlue,"#tau #rightarrow e / #mu");
  addSample("TTbarJets:nElectrons==0&&nMuons==0&&(ttbarDecayCode==8||ttbarDecayCode==9)",kRed-5,"misc #tau");

  selection_ =TCut("MET>=150 && cutPV==1 && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1");
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 20; low=300; high=1000;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_htshapes_1l_split",0,"GeV");

  //njets
  var="njets"; xtitle="n jets (p_{T}>50 GeV)";
  nbins = 10; low=0; high=10;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_njetsshapes_1l_split",0);

  //HT in PV bins

  clearSamples();
  addSample("TTbarJets:nElectrons==0&&nMuons==0&&nGoodPV<=7",kRed,"PV<=7");
  addSample("TTbarJets:nElectrons==0&&nMuons==0&&nGoodPV>=8&&nGoodPV<=11",kGreen,"PV 8-11");
  addSample("TTbarJets:nElectrons==0&&nMuons==0&&nGoodPV>=12",kBlue,"PV>=12");

  selection_ =TCut("MET>=150 && cutPV==1 && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1");
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 20; low=300; high=1000;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_htshapes_pv",0,"GeV");

  //HT distribution for only 3 jet events

  clearSamples();
  addSample("TTbarJets:nElectrons==0&&nMuons==0",kRed,"0 lep");
  addSample("TTbarJets:(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)",kBlue,"1 lep");

  selection_ =TCut("MET>=150 && cutPV==1 && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1 &&njets==3");
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 35; low=300; high=1000;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_htshapes_0l_1l_3jets",0,"GeV");

  selection_ =TCut("MET>=150 && cutPV==1 && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1 &&njets==4");
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_htshapes_0l_1l_4jets",0,"GeV");

  //HT for e and mu
  setStackMode(false,true);//normalized

  selection_ =TCut("MET>=150 && cutPV==1 && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1");
  clearSamples();
  addSample("TTbarJets:nElectrons==0&&nMuons==0",kRed,"0 lep");
  addSample("TTbarJets:(((nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)",kBlue,"1 e");
  addSample("TTbarJets:(((nElectrons==0 && nMuons==1)) && MT_Wlep>=0&&MT_Wlep<100)",kMagenta,"1 #mu");

  selection_ =TCut("MET>=150 && cutPV==1 && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1");
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 35; low=300; high=1000;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_htshapes_0l_1em",0,"GeV");

  //
  setStackMode(false,true);//normalized
  selection_ =TCut("MET>=150 && cutPV==1 && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets>=1");
  clearSamples();
  addSample("TTbarJets:(((nElectrons==0 && nMuons==1)) && MT_Wlep>=0&&MT_Wlep<100)&&muonpt1<25",kRed,"#mu pT<25");
  addSample("TTbarJets:(((nElectrons==0 && nMuons==1)) && MT_Wlep>=0&&MT_Wlep<100)&&muonpt1>=25",kBlue,"#mu pT>25");

  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 35; low=300; high=1000;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_htshapes_bypt_1mu",0,"GeV");

  var="jetpt1"; xtitle="lead jet pt";
  nbins = 20; low=0; high=500;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_jetptshapes_bypt_1mu",0);


  // new idea -- compare =0b with b-tagged
  clearSamples();
  addSample("TTbarJets:nElectrons==0&&nMuons==0&&nbjets==0",kBlue,"==0 b");
  addSample("TTbarJets:nElectrons==0&&nMuons==0&&nbjets==1",kRed,"==1 b");
  addSample("TTbarJets:nElectrons==0&&nMuons==0&&nbjets>=2",kMagenta,">=2 b");
  //  addSample("TTbarJets:nElectrons==0&&nMuons==0&&nbjets>=3",kCyan,">=3 b");

  selection_ =TCut("HT>=400 && cutPV==1 && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1");
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 30; low=150; high=450;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_metshapes_eq0b_ge1b",0,"GeV");


  selection_ =TCut("MET>=125 && cutPV==1 && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1");
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 35; low=400; high=1200;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_htshapes_eq0b_ge1b",0,"GeV");



}

void RA2b_tauveto() {


  nameOfEventWeight_="weight2"; //for 2012 cfA skimmed ntuples
  loadSamples(true,"ra2b2012");
  usePUweight_=true;
  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getDefault();

  int nbins;
  float low,high;
  TString var,xtitle;
  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2;
  dodata_=false;

  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2";
  TCut ht = "HT>=400";
  TCut met = "MET>=125";
  TCut mdp = "minDeltaPhiN>=4";
  //  TCut sl = getSingleLeptonCut();
  TCut zl = getLeptonVetoCut();

  setDatasetToDraw("2012hybridplus");
  lumiScale_=12034.14; //Run 2012 A+B+C, i think

  btagSFweight_="probge1";
  if (false) {
    clearSamples();
    setStackMode(false,true);
    addSample("TTbarJets0:ttbarDecayCode==7",kBlue,"#tau #rightarrow h");
    addSample("TTbarJets0:ttbarDecayCode==5||ttbarDecayCode==6||ttbarDecayCode==8||ttbarDecayCode==9",kGreen,"#tau (all other)");
    addSample("TTbarJets0:ttbarDecayCode==1||ttbarDecayCode==3||ttbarDecayCode==4",kRed,"e and/or #mu");
    addSample("TTbarJets0:ttbarDecayCode==2",kMagenta,"hadronic");

    addSample("T1bbbb$1300$725",kBlack,"T1bbbb (1300,725)");

    selection_ = baseline&& cleaning &&met && zl && ht && mdp;
    
    nbins=4; low=0;high=4;
    var="nIsoTracks15_005_03"; xtitle=var;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_isotrack15_005_03",0);

    nbins=40; low=0; high=1;
    var="minTrackIso15"; xtitle=var;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_mintrackiso15",0);
    setStackMode(false,false);
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_mintrackiso15",0);
  }
  savePlots_=false;
  clearSamples();
  addSample("TTbarJets0");
  TString signalpoint = "T1bbbb$1300$725";
  addSample(signalpoint);

  //baseline cuts
  TGraph* gr_pog = new TGraph();
  TGraph* gr_ind = new TGraph();
  TGraph* gr_tk20_3 = new TGraph();
  TGraph* gr_tk15_3 = new TGraph();
  TGraph* gr_tk10_3 = new TGraph();
  TGraph* gr_tk5_3 = new TGraph();
  TGraph* gr_tk20_5 = new TGraph();
  TGraph* gr_tk15_5 = new TGraph();
  TGraph* gr_tk10_5 = new TGraph();
  TGraph* gr_tk5_5 = new TGraph();
  TGraph* gr_tk20_7 = new TGraph();
  TGraph* gr_tk15_7 = new TGraph();
  TGraph* gr_tk10_7 = new TGraph();
  TGraph* gr_tk5_7 = new TGraph();
//   TGraph* gr_pf10 = new TGraph();
//   TGraph* gr_pf5 = new TGraph();

  gr_pog->SetName("pog");
  gr_ind->SetName("ind");
  gr_tk20_3->SetName("tk20_3");
  gr_tk15_3->SetName("tk15_3");
  gr_tk10_3->SetName("tk10_3");
  gr_tk5_3->SetName("tk5_3");
  gr_tk20_5->SetName("tk20_5");
  gr_tk15_5->SetName("tk15_5");
  gr_tk10_5->SetName("tk10_5");
  gr_tk5_5->SetName("tk5_5");
  gr_tk20_7->SetName("tk20_7");
  gr_tk15_7->SetName("tk15_7");
  gr_tk10_7->SetName("tk10_7");
  gr_tk5_7->SetName("tk5_7");
  //  gr_pf10->SetName("pf10");
  //  gr_pf5->SetName("pf5");

  const int nvariations = 4 + 4 +4*4*3 +1; //57
  double baselinett,baselinesig;
  for (int igr=1; igr<= nvariations; igr++) {

    //crude but it works

    TCut tauvetocut = "(1)";
    if (igr==1) {}
    else if (igr==2) tauvetocut="nTausVLoose==0";
    else if (igr==3) tauvetocut="nTausLoose==0";
    else if (igr==4) tauvetocut="nTausMedium==0";
    else if (igr==5) tauvetocut="nTausTight==0";
    else if (igr==6) tauvetocut="nIndirectTaus2==0";
    else if (igr==7) tauvetocut="nIndirectTaus3==0";
    else if (igr==8) tauvetocut="nIndirectTaus4==0";
    else if (igr==9) tauvetocut="nIndirectTaus5==0";

    else if (igr==10) tauvetocut="nIsoTracks20_020_03==0";
    else if (igr==11) tauvetocut="nIsoTracks20_015_03==0";
    else if (igr==12) tauvetocut="nIsoTracks20_010_03==0";
    else if (igr==13) tauvetocut="nIsoTracks20_005_03==0";
    else if (igr==14) tauvetocut="nIsoTracks15_020_03==0";
    else if (igr==15) tauvetocut="nIsoTracks15_015_03==0";
    else if (igr==16) tauvetocut="nIsoTracks15_010_03==0";
    else if (igr==17) tauvetocut="nIsoTracks15_005_03==0";
    else if (igr==18) tauvetocut="nIsoTracks10_020_03==0";
    else if (igr==19) tauvetocut="nIsoTracks10_015_03==0";
    else if (igr==20) tauvetocut="nIsoTracks10_010_03==0";
    else if (igr==21) tauvetocut="nIsoTracks10_005_03==0";
    else if (igr==22) tauvetocut="nIsoTracks5_020_03==0";
    else if (igr==23) tauvetocut="nIsoTracks5_015_03==0";
    else if (igr==24) tauvetocut="nIsoTracks5_010_03==0";
    else if (igr==25) tauvetocut="nIsoTracks5_005_03==0";

    else if (igr==26) tauvetocut="nIsoTracks20_020_05==0";
    else if (igr==27) tauvetocut="nIsoTracks20_015_05==0";
    else if (igr==28) tauvetocut="nIsoTracks20_010_05==0";
    else if (igr==29) tauvetocut="nIsoTracks20_005_05==0";
    else if (igr==30) tauvetocut="nIsoTracks15_020_05==0";
    else if (igr==31) tauvetocut="nIsoTracks15_015_05==0";
    else if (igr==32) tauvetocut="nIsoTracks15_010_05==0";
    else if (igr==33) tauvetocut="nIsoTracks15_005_05==0";
    else if (igr==34) tauvetocut="nIsoTracks10_020_05==0";
    else if (igr==35) tauvetocut="nIsoTracks10_015_05==0";
    else if (igr==36) tauvetocut="nIsoTracks10_010_05==0";
    else if (igr==37) tauvetocut="nIsoTracks10_005_05==0";
    else if (igr==38) tauvetocut="nIsoTracks5_020_05==0";
    else if (igr==39) tauvetocut="nIsoTracks5_015_05==0";
    else if (igr==40) tauvetocut="nIsoTracks5_010_05==0";
    else if (igr==41) tauvetocut="nIsoTracks5_005_05==0";


    else if (igr==42) tauvetocut="nIsoTracks20_020_07==0";
    else if (igr==43) tauvetocut="nIsoTracks20_015_07==0";
    else if (igr==44) tauvetocut="nIsoTracks20_010_07==0";
    else if (igr==45) tauvetocut="nIsoTracks20_005_07==0";
    else if (igr==46) tauvetocut="nIsoTracks15_020_07==0";
    else if (igr==47) tauvetocut="nIsoTracks15_015_07==0";
    else if (igr==48) tauvetocut="nIsoTracks15_010_07==0";
    else if (igr==49) tauvetocut="nIsoTracks15_005_07==0";
    else if (igr==50) tauvetocut="nIsoTracks10_020_07==0";
    else if (igr==51) tauvetocut="nIsoTracks10_015_07==0";
    else if (igr==52) tauvetocut="nIsoTracks10_010_07==0";
    else if (igr==53) tauvetocut="nIsoTracks10_005_07==0";
    else if (igr==54) tauvetocut="nIsoTracks5_020_07==0";
    else if (igr==55) tauvetocut="nIsoTracks5_015_07==0";
    else if (igr==56) tauvetocut="nIsoTracks5_010_07==0";
    else if (igr==57) tauvetocut="nIsoTracks5_005_07==0";


    else assert(0);

    selection_ = baseline&& cleaning && met && zl && ht && mdp&&tauvetocut;
    var = "HT"; xtitle="HT";
    nbins = 20; low = 400; high = 4000;
    drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
    
    double ntt=  getIntegral("TTbarJets0");
    double nsig=  getIntegral(signalpoint);

    if (igr==1) { baselinett=ntt; baselinesig=nsig;}
    else {
      TGraph * gr=0;
      if (igr>=2 && igr<=5) gr=gr_pog;
      else if (igr>=6 && igr<=9) gr=gr_ind;
      else if (igr>=10 && igr<=13) gr=gr_tk20_3;
      else if (igr>=14 && igr<=17) gr=gr_tk15_3;
      else if (igr>=18 && igr<=21) gr=gr_tk10_3;
      else if (igr>=22 && igr<=25) gr=gr_tk5_3;

      else if (igr>=26 && igr<=29) gr=gr_tk20_5;
      else if (igr>=30 && igr<=33) gr=gr_tk15_5;
      else if (igr>=34 && igr<=37) gr=gr_tk10_5;
      else if (igr>=38 && igr<=41) gr=gr_tk5_5;

      else if (igr>=42 && igr<=45) gr=gr_tk20_7;
      else if (igr>=46 && igr<=49) gr=gr_tk15_7;
      else if (igr>=50 && igr<=53) gr=gr_tk10_7;
      else if (igr>=54 && igr<=57) gr=gr_tk5_7;
      //      else if (igr>=18 && igr<=21) gr=gr_pf10;
      //      else if (igr>=22 && igr<=25) gr=gr_pf5;
      gr->SetPoint(gr->GetN(),nsig/baselinesig,ntt/baselinett);
    }
  }

  TFile fout("tauVetoGraph.root","recreate");

   gr_pog->Write();
   gr_ind->Write();
   gr_tk20_3->Write();
   gr_tk15_3->Write();
   gr_tk10_3->Write();
   gr_tk5_3->Write();

   gr_tk20_5->Write();
   gr_tk15_5->Write();
   gr_tk10_5->Write();
   gr_tk5_5->Write();

   gr_tk20_7->Write();
   gr_tk15_7->Write();
   gr_tk10_7->Write();
   gr_tk5_7->Write();
 //   gr_pf10->Write();
//    gr_pf5->Write();

  fout.Write();
  fout.Close();

}


void RA2b_investigation_tracking() {

  nameOfEventWeight_="weight2"; //for 2012 cfA skimmed ntuples
  loadSamples(true,"ra2b2012");
  usePUweight_=true;


  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getDefault();

  resetSamples();
  removeSample("LM9");
  //  addSample("T1bbbb$1200$300");
  removeSample("ZJets");

  int nbins;
  float low,high;
  TString var,xtitle;
  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2;
  dodata_=true;

  setSampleScaleFactor("TTbarJets",0.9);
  setSampleScaleFactor("WJets",0.83);
  setSampleScaleFactor("PythiaPUQCD",1.8);

  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2";
  TCut invertcleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET>2";
  TCut ht = "HT>=1000";
  TCut met = "MET>=250";
  TCut mdp = "minDeltaPhiN>=4";
  TCut ldp = "minDeltaPhiN<4";
  TCut sl = getSingleLeptonCut();
  TCut zl = getLeptonVetoCut();



  setDatasetToDraw("2012hybridplus");
  lumiScale_=12034.14; //Run 2012 A+B+C, i think

  btagSFweight_="probge1";

  //sample of bad events
  selection_ = baseline&& invertcleaning && met && zl && ht && ldp;
  var = "jetchargedhadronmult1"; xtitle="lead jet ch had mult";
  nbins = 20; low = 0; high = 100;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_TOBTEC_chhadmult1_ldp_highPFoverCalo",0);

  selection_ = baseline&& invertcleaning && met && zl && ht && ldp;
  var = "jetchargedhadronmult2"; xtitle="2nd jet ch had mult";
  nbins = 20; low = 0; high = 100;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_TOBTEC_chhadmult2_ldp_highPFoverCalo",0);

  selection_ = baseline&& invertcleaning && met && zl && ht && ldp;
  var = "jetchargedhadronmult3"; xtitle="3rd jet ch had mult";
  nbins = 20; low = 0; high = 100;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_TOBTEC_chhadmult3_ldp_highPFoverCalo",0);

  selection_ = baseline&& invertcleaning && met && zl && ht && ldp;
  var = "jetneutralhadronmult1"; xtitle="lead jet netral had mult";
  nbins = 20; low = 0; high = 100;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_TOBTEC_neuhadmult1_ldp_highPFoverCalo",0);


}


void RA2b_investigation_noise() {

  //use the trees from keith. i don't trust mine for data
  inputPath = "/cu4/ra2b/reducedTrees/v66_2/"; //2012
  dataInputPath =  "/cu4/ra2b/reducedTrees/v66_2/ORIGINALS/";

  nameOfEventWeight_="weight2"; //for 2012 cfA skimmed ntuples
  loadSamples(true,"ra2b2012");
  usePUweight_=true;


  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getDefault();

  resetSamples();
  removeSample("LM9");

  int nbins;
  float low,high;
  TString var,xtitle;
  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2;
  dodata_=true;

  //4 Nov -- updating these a bit
  setSampleScaleFactor("TTbarJets",0.9);
  setSampleScaleFactor("WJets",0.9);
  setSampleScaleFactor("PythiaPUQCD",1.8);

  //new plot for david
  lumiScale_=12034.14; //Run 2012 A+B+C
  setDatasetToDraw("2012hybridplus");


  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut invertedcleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET>2";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2";
  TCut invertedcleaning2="passCleaning==0&&buggyEvent==0&&MET/caloMET<2";
  //  TCut htloose = "HT>=400";
  TCut httight = "HT>=800";
  TCut metloose = "MET>=125";
  TCut mettight = "MET>=250";
//   TCut mdp = "minDeltaPhiN>=4";
//   TCut ldp = "minDeltaPhiN<4";
//   TCut btag = "nbjets>=1";
//   TCut bveto = "nbjets==0";
  TCut zl = getLeptonVetoCut();

  //no btag, no minDp cut, tightish HT
  selection_=baseline&&invertedcleaning&&httight&&mettight&&zl;
  var="minDeltaPhiN"; xtitle="minDeltaPhiN";
  nbins = 20; low=0; high=20;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_anomalousMET_minDeltaPhiN_ht800_met250_ge0b",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_anomalousMET_minDeltaPhiN_ht800_met250_ge0b",0);

  //also for david
  selection_=baseline&&invertedcleaning2&&httight&&mettight&&zl;
  var="minDeltaPhiN"; xtitle="minDeltaPhiN";
  nbins = 40; low=0; high=40;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_failFilters_minDeltaPhiN_ht800_met250_ge0b",0);
 
  nbins = 10; low=0; high=40;
  //now dig in -- no events here
  selection_=baseline&&TCut("csctighthaloFilter==0&&passCleaning==0&&MET/caloMET<2&&buggyEvent==0")&&httight&&mettight&&zl;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_failCSC_minDeltaPhiN_ht800_met250_ge0b",0);

  selection_=baseline&&TCut("greedymuonFilter==0&&passCleaning==0&&MET/caloMET<2&&buggyEvent==0")&&httight&&mettight&&zl;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_greedymuon_minDeltaPhiN_ht800_met250_ge0b",0);

  //here is the culprit
  selection_=baseline&&TCut("hbhenoiseFilter==0&&passCleaning==0&&MET/caloMET<2&&buggyEvent==0")&&httight&&mettight&&zl;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_HBHEnoise_minDeltaPhiN_ht800_met250_ge0b",0);

  selection_=baseline&&TCut("inconsistentmuonFilter==0&&passCleaning==0&&MET/caloMET<2&&buggyEvent==0")&&httight&&mettight&&zl;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_inconsistentMuon_minDeltaPhiN_ht800_met250_ge0b",0);
 
  selection_=baseline&&TCut("ra2ecaltpFilter==0&&passCleaning==0&&MET/caloMET<2&&buggyEvent==0")&&httight&&mettight&&zl;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_TPfilter_minDeltaPhiN_ht800_met250_ge0b",0);

  selection_=baseline&&TCut("scrapingvetoFilter==0&&passCleaning==0&&MET/caloMET<2&&buggyEvent==0")&&httight&&mettight&&zl;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_scrapingVeto_minDeltaPhiN_ht800_met250_ge0b",0);

  selection_=baseline&&TCut("trackingfailureFilter==0&&passCleaning==0&&MET/caloMET<2&&buggyEvent==0")&&httight&&mettight&&zl;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_trackingFailure_minDeltaPhiN_ht800_met250_ge0b",0);

  selection_=baseline&&TCut("hbhenoiseFilter==1&&trackingfailureFilter==0&&passCleaning==0&&MET/caloMET<2&&buggyEvent==0")&&httight&&mettight&&zl;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_trackingFailurePassHBHE_minDeltaPhiN_ht800_met250_ge0b",0);

  selection_=baseline&&TCut("badjetFilter==0&&passCleaning==0&&MET/caloMET<2&&buggyEvent==0")&&httight&&mettight&&zl;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_badjet_minDeltaPhiN_ht800_met250_ge0b",0);

  selection_=baseline&&TCut("badjetFilter==0&&hbhenoiseFilter==1&&passCleaning==0&&MET/caloMET<2&&buggyEvent==0")&&httight&&mettight&&zl;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_badjetPassHBHE_minDeltaPhiN_ht800_met250_ge0b",0);

  selection_=baseline&&TCut("PBNRcode<=0&&passCleaning==0&&MET/caloMET<2&&buggyEvent==0")&&httight&&mettight&&zl;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_failPBNR_minDeltaPhiN_ht800_met250_ge0b",0);

  selection_=baseline&&TCut("PBNRcode<=0&&hbhenoiseFilter==1&&passCleaning==0&&MET/caloMET<2&&buggyEvent==0")&&httight&&mettight&&zl;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_failPBNRPassHBHE_minDeltaPhiN_ht800_met250_ge0b",0);


  selection_=baseline&&TCut("hbhenoiseFilter==1&&passCleaning==0&&MET/caloMET<2&&buggyEvent==0")&&httight&&mettight&&zl;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_failCleaningPassHBHE_minDeltaPhiN_ht800_met250_ge0b",0);


  addSample("T1bbbb$1100$700");
  btagSFweight_="probge1";
  selection_=baseline&&cleaning&&TCut("HT>400&&HT<500&&MET>350")&&zl;
  var="minDeltaPhiN"; xtitle="minDeltaPhiN";
  nbins = 20; low=0; high=60;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhiN_lowhtHighmet_ge1b",0);


  //original studies are down here

  //first confirm owen's findings
  setDatasetToDraw("JetHT");
  lumiScale_=12034.14; //Run 2012 A+B+C, i think

  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 40; low=100; high=500;
  setLogY(true);

  //sl sample should look normal
  selection_ =TCut("HT>=1000 && cutPV==1 && (pass_PFNoPUHT650==1||pass_PFHT650==1) && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0")&&getSingleLeptonCut();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jethtPD_ht1000_pfht650trig_met_sl",0,"GeV");

  //ldp should look normal
  selection_ =TCut("HT>=1000 && cutPV==1 && (pass_PFNoPUHT650==1||pass_PFHT650==1) && njets>=3  && minDeltaPhiN < 4 && passCleaning==1 && nbjets==0")&&getLeptonVetoCut();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jethtPD_ht1000_pfht650trig_met_ldp",0,"GeV");

  //zl weird?
  selection_ =TCut("HT>=1000 && cutPV==1 && (pass_PFNoPUHT650==1||pass_PFHT650==1) && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0")&&getLeptonVetoCut();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jethtPD_ht1000_pfht650trig_met_0l",0,"GeV");

  savePlots_=false;
  selection_ =TCut("HT>=1000 && cutPV==1 && (pass_PFNoPUHT650==1||pass_PFHT650==1) && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0")&&getLeptonVetoCut();
  var="nGoodPV"; xtitle="num of good PV";
  nbins = 20; low=0; high=40;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jethtPD_ht1000_pfht650trig_npv",0);

  selection_ =TCut("HT>=1000 && cutPV==1 && (pass_PFNoPUHT650==1||pass_PFHT650==1) && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0 && MET>=250")&&getLeptonVetoCut();
  var="njets"; xtitle="n jets (pT>50 GeV)";
  nbins = 10; low=0; high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jethtPD_ht1000_pfht650trig_njets",0);

 //jet properties
  maxScaleFactor_=1.3;
  selection_ =TCut("HT>=1000 && cutPV==1 && (pass_PFNoPUHT650==1||pass_PFHT650==1) && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0 && MET>=250 && njets>=3")&&getLeptonVetoCut();
  var = "jetchargedhadronfrac1"; xtitle="lead jet ch had frac";
  nbins = 20; low = 0; high = 1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jethtPD_ht1000_pfht650trig_jet1chhadfrac",0);

  //now try SL sample
  selection_ =TCut("HT>=1000 && cutPV==1 && (pass_PFNoPUHT650==1||pass_PFHT650==1) && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0 && MET>=250")&&getSingleLeptonCut();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jethtPD_ht1000_pfht650trig_jet1chhadfrac_sl",0);
  //compare to LDP
  selection_ =TCut("HT>=1000 && cutPV==1 && (pass_PFNoPUHT650==1||pass_PFHT650==1) && njets>=3  && minDeltaPhiN < 4 && passCleaning==1 && nbjets==0 && MET>=250")&&getLeptonVetoCut();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jethtPD_ht1000_pfht650trig_jet1chhadfrac_ldp",0);

  //jet charged hadron mult
  selection_ =TCut("HT>=1000 && cutPV==1 && (pass_PFNoPUHT650==1||pass_PFHT650==1) && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0 && MET>=250 && njets>=3")&&getLeptonVetoCut();
  var = "jetchargedhadronmult1"; xtitle="lead jet ch had mult";
  nbins = 25; low = 0; high = 50;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jethtPD_ht1000_pfht650trig_jet1chhadmult",0);

  //delta Phi*
  setLogY(true);
  selection_ =TCut("HT>=1000 && cutPV==1 && (pass_PFNoPUHT650==1||pass_PFHT650==1) && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0 && MET>=250 && njets>=3")&&getLeptonVetoCut();
  var = "deltaPhiStar"; xtitle="#Delta#phi *";
  nbins = 30; low = 0; high = 2;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jethtPD_ht1000_pfht650trig_deltaphistar",0);

  selection_ =TCut("HT>=1000 && cutPV==1 && (pass_PFNoPUHT650==1||pass_PFHT650==1) && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0 && MET>=250 && njets>=3")&&getLeptonVetoCut();
  var = "minDeltaPhi"; xtitle="minDeltaPhi";
  nbins = 30; low = 0; high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jethtPD_ht1000_pfht650trig_minDeltaPhi",0);

  //lead jet pt
  setLogY(true);
  selection_ =TCut("HT>=1000 && cutPV==1 && (pass_PFNoPUHT650==1||pass_PFHT650==1) && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0 && MET>=250 && njets>=3")&&getLeptonVetoCut();
  var = "jetpt1"; xtitle="lead jet pT";
  nbins = 40; low = 0; high = 2000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jethtPD_ht1000_pfht650trig_jetpt1",0);

  selection_ =TCut("HT>=1000 && cutPV==1 && (pass_PFNoPUHT650==1||pass_PFHT650==1) && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0 && MET>=250 && njets>=3")&&getLeptonVetoCut();
  var = "jetpt2"; xtitle="2nd jet pT";
  nbins = 40; low = 0; high = 1000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jethtPD_ht1000_pfht650trig_jetpt2",0);

  selection_ =TCut("HT>=1000 && cutPV==1 && (pass_PFNoPUHT650==1||pass_PFHT650==1) && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0 && MET>=250 && njets>=3")&&getLeptonVetoCut();
  var = "deltaPhi1"; xtitle="deltaPhi(lead jet, MET)";
  nbins = 30; low = 0; high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jethtPD_ht1000_pfht650trig_deltaPhi1",0);

  //run 2012C only
  lumiScale_=6806;
  maxScaleFactor_=1.3;
  selection_ =TCut("HT>=1000 && cutPV==1 && pass_PFNoPUHT650==1 && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0 && MET>=250 && njets>=3")&&getLeptonVetoCut();
  var = "jetchargedhadronfrac1"; xtitle="lead jet ch had frac";
  nbins = 20; low = 0; high = 1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jethtPD_ht1000_pfht650trig_jet1chhadfrac_Run2012C",0);

  lumiScale_=12034.14 - 6806; //Run 2012 A+B
  maxScaleFactor_=1.3;
  selection_ =TCut("HT>=1000 && cutPV==1 && pass_PFHT650==1 && minDeltaPhiN >= 4 && passCleaning==1 && nbjets==0 && MET>=250 && njets>=3")&&getLeptonVetoCut();
  var = "jetchargedhadronfrac1"; xtitle="lead jet ch had frac";
  nbins = 20; low = 0; high = 1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jethtPD_ht1000_pfht650trig_jet1chhadfrac_Run2012AB",0);

}

void RA2bSignalRegion() {

  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  loadSamples(true,"ra2b2012");
  usePUweight_=true;

  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getDefault();

  resetSamples();
  removeSample("LM9");

  removeSample("VV");
  removeSample("ZJets");

  addSample("VV",kCyan+1,"rare");
  chainSamples("VV","TTV");
  chainSamples("VV","ZJets");

  addSample("T1bbbb$1300$100");
  useMassInLegend_=false;

  //use the 'owen factors' verbatim
  setSampleScaleFactor("TTbarJets",0.9);
  setSampleScaleFactor("WJets",0.9);
  setSampleScaleFactor("PythiaPUQCD",1.8);

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2;
  dodata_=true;

  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40";
  TCut ht = "HT>=1000";
  TCut met = "MET>=150";
  TCut mdp = "minDeltaPhiN_asin>=4";
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");

  setDatasetToDraw("2012hybridplus");
  lumiScale_=12034.14 + 5580; //Run 2012 ABC+D

  //  const  TString btagselection = "probge3";
  const  TString btagselection = "probge3";

  btagSFweight_=btagselection; //use btag sf

  TString desc = "zl_";
  if (btagselection =="probge3") desc+="3b";
  if (btagselection =="probge2") desc+="ge2b";
  if (btagselection =="prob1") desc+="eq1b";
  desc+="_met150_ht1000";

  //MET
  selection_ = baseline&&cleaning && mdp && zl &&ht;
  var="MET"; xtitle="MET [GeV]";
  nbins = 7; low=150; high=500;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_"+desc,0,"GeV");
//   setLogY(true);
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_"+desc,0,"GeV");

  //met significance (after MET cut)
  selection_ = baseline&&cleaning  && mdp && zl &&ht &&met;
  var="METsig"; xtitle="E_{T}^{miss} significance";
  nbins = 5; low=0; high=400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metsig_"+desc,0,"GeV");
  //  setLogY(true);
  //  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metsig_"+desc,0,"GeV");

  //HT
  selection_ = baseline&&cleaning && met && mdp && zl ;
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 12; low=800; high=2000;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_"+desc,0,"GeV");
//   setLogY(true);
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_"+desc,0,"GeV");

  //meff
  selection_ = baseline&&cleaning && met && mdp && zl && ht;
  var="HT+MET"; xtitle="m_{eff} [GeV]";
  nbins = 6; low=1150; high=2350;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_meff_"+desc,0,"GeV");
//   setLogY(true);
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_meff_"+desc,0,"GeV");


  //nPV -- just as a check
  selection_ = baseline&&cleaning && met && mdp && zl && ht;
  var="nGoodPV"; xtitle="num of good PV";
  nbins = 5; low=0; high=40;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_npv_"+desc,0);

  //njets
  selection_ = TCut("cutPV==1 && cutTrigger2==1 && jetpt2>=70") && met && mdp && zl && ht &&cleaning;
  var="njets"; xtitle="n jets (p_{T} > 50 GeV)";
  nbins = 6; low=3; high=9;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets_"+desc,0);

  //njets30
  selection_ = TCut("cutPV==1 && cutTrigger2==1  && jetpt2>=70") && met && mdp && zl && ht &&cleaning;
  var="njets30"; xtitle="n jets (p_{T} > 30 GeV)";
  nbins = 6; low=0; high=12;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets30_"+desc,0);

  //nbjets BLIND
//   btagSFweight_="";
//   selection_ = baseline&&cleaning && met && mdp && zl && ht;
//   var="nbjets"; xtitle="n CSVM b jets (pT>50 GeV)";
//   nbins = 5; low=0; high=5;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_nbjets_"+desc,0);
//   setLogY(true);
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_nbjets_"+desc,0);
//   setLogY(false);

//   btagSFweight_=btagselection;

  //minDeltaPhiN
  selection_ = baseline&&cleaning  && zl &&ht &&met;
  var="minDeltaPhiN_asin"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 10; low=0; high=40;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhiN_"+desc,0);

  //now the fun stuff

  //jet pt1
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "jetpt1"; xtitle="lead jet pT [GeV]";
  nbins = 10; low = 50; high = 1050;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetpt1_"+desc,0,"GeV");

  maxScaleFactor_=1.3;
  //jet phi1
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "jetphi1"; xtitle="lead jet phi";
  nbins = 10; low = -TMath::Pi(); high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetphi1_"+desc,0);
  //jet eta1
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "jeteta1"; xtitle="lead jet eta";
  nbins = 10; low = -3; high = 3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jeteta1_"+desc,0);
  maxScaleFactor_=1;

  //jet pt2
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "jetpt2"; xtitle="2nd jet pT [GeV]";
  nbins = 7; low = 50; high = 750;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetpt2_"+desc,0,"GeV");
  //jet eta2
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "jeteta2"; xtitle="2nd jet eta";
  nbins = 10; low = -3; high = 3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jeteta2_"+desc,0);
  maxScaleFactor_=1;

 
  //jet pt3
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "jetpt3"; xtitle="3rd jet pT [GeV]";
  nbins = 5; low = 50; high = 550;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetpt3_"+desc,0,"GeV");

  //met phi
  maxScaleFactor_=1.3;
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "METphi"; xtitle="phi of MET";
  nbins = 10; low =-TMath::Pi(); high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metphi_"+desc,0);
  maxScaleFactor_=1.05;

  //calo met check
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "MET/caloMET"; xtitle="TypeIPFMET/caloMET";
  nbins = 10; low =0; high = 2;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metOverCalomet_"+desc,0);

  maxScaleFactor_=2;
  //DeltaPhi(j1,MET)
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "deltaPhi1"; xtitle="DeltaPhi(j1,MET)";
  nbins = 10; low =0; high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_deltaphi1_"+desc,0);

  //DeltaPhi(j2,MET)
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "deltaPhi2"; xtitle="DeltaPhi(j2,MET)";
  nbins = 10; low =0; high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_deltaphi2_"+desc,0);

  maxScaleFactor_=1.05;

  //raw met check
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "MET/rawPFMET"; xtitle="TypeIPFMET/rawPFMET";
  nbins = 10; low =0.5; high = 1.5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metOverRawmet_"+desc,0);

  //jet properties
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "jetchargedhadronfrac1"; xtitle="lead jet ch had frac";
  nbins = 10; low = 0; high = 1;
  maxScaleFactor_=1.3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jet1chhadfrac_"+desc,0);
  maxScaleFactor_=1.05;

  selection_ = baseline&&cleaning && met && zl && ht && mdp;
  var = "jetchargedhadronmult1"; xtitle="lead jet ch had mult";
  nbins = 10; low = 0; high = 60;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jet1chhadmult_"+desc,0);

 //jet properties
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "jetchargedhadronfrac2"; xtitle="2nd jet ch had frac";
  nbins = 10; low = 0; high = 1;
  maxScaleFactor_=1.3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jet2chhadfrac_"+desc,0);
  maxScaleFactor_=1.05;

  selection_ = baseline&&cleaning && met && zl && ht && mdp;
  var = "jetchargedhadronmult2"; xtitle="2nd jet ch had mult";
  nbins = 10; low = 0; high = 90;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jet2chhadmult_"+desc,0);

   //minDeltaPhi
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "minDeltaPhi"; xtitle="minDeltaPhi";
  nbins = 10; low = 0; high = TMath::Pi();
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhi_"+desc,0);
  setLogY(false);

  selection_ = baseline&&cleaning && met && zl && ht ;
  var = "minDeltaPhi"; xtitle="minDeltaPhi (no mdpN cut)";
  nbins = 10; low = 0; high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhi_nominDeltaPhiNcut_"+desc,0);

  //deltaphi*
  selection_ = baseline&&cleaning && met && zl && ht && mdp;
  var = "deltaPhiStar"; xtitle="#Delta#phi *";
  nbins = 10; low = 0; high =1;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_deltaPhiStar_"+desc,0);
  setLogY(false);

}

void ra2bControlttbarcomp() {

  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  loadSamples(true,"ra2b2012");
  usePUweight_=true;

  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getDefault();

  resetSamples();
  removeSample("LM9");

  addSample("TTV",kGreen-9,"TTV");

  setSampleScaleFactor("WJets",0.9); //use the owen factor instead of the LO factor (which would be 0.84)
  setSampleScaleFactor("PythiaPUQCD",1.8);

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0.7; ratioMax = 1.3;
  dodata_=true;

  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40";
  TCut ht = "HT>=400";
  TCut met = "MET>=175";
  TCut mdp = "minDeltaPhiN>=4";
  TCut sl = getSingleLeptonCut();
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");
  //  TCut withisotk = TCut("nIsoTracks15_005_03>=1");

  TString desc;

  setDatasetToDraw("2012hybridplus");
  lumiScale_=12034.14; //Run 2012 A+B+C, i think

  //njets
 
  //do this so that ttbar will be last in the stack
  desc = "madgraph";
  removeSample("TTbarJets");
  addSample("TTbarJets");
  setSampleScaleFactor("TTbarJets",0.9); //use the 'owen factor'

  selection_ = TCut("cutPV==1 && cutTrigger2==1 && jetpt2>=70") && met && mdp && sl && ht && cleaning;
  btagSFweight_="probge1"; //use btag sf
  var="njets"; xtitle="n jets (pT>50 GeV)";
  nbins = 10; low=0; high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets_"+desc,0);

  selection_ = baseline && met && mdp && sl && ht && cleaning &&TCut("nMuons==1");
  btagSFweight_="probge1"; //use btag sf
  var="muonpt1"; xtitle="muon pT";
  setLogY(true);
  nbins = 20; low=0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_mupt_"+desc,0);
  setLogY(false);

  // -- now powheg --
  desc = "powheg";
  removeSample("TTbarJets");
  addSample("TTbarJetsPowheg");
  //  setSampleScaleFactor("TTbarJets",0.9); //use the 'owen factor'

  selection_ = TCut("cutPV==1 && cutTrigger2==1 && jetpt2>=70") && met && mdp && sl && ht && cleaning;
  btagSFweight_="probge1"; //use btag sf
  var="njets"; xtitle="n jets (pT>50 GeV)";
  nbins = 10; low=0; high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets_"+desc,0);


  selection_ = baseline && met && mdp && sl && ht && cleaning &&TCut("nMuons==1");
  btagSFweight_="probge1"; //use btag sf
  var="muonpt1"; xtitle="muon pT";
  setLogY(true);
  nbins = 20; low=0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_mupt_"+desc,0);
  setLogY(false);


  // -- now mc@nlo --
  desc = "mcnlo";
  removeSample("TTbarJetsPowheg");
  addSample("TTbarJetsMCNLO");
  //  setSampleScaleFactor("TTbarJets",0.9); //use the 'owen factor'

  selection_ = TCut("cutPV==1 && cutTrigger2==1 && jetpt2>=70") && met && mdp && sl && ht && cleaning;
  btagSFweight_="probge1"; //use btag sf
  var="njets"; xtitle="n jets (pT>50 GeV)";
  nbins = 10; low=0; high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets_"+desc,0);


  selection_ = baseline && met && mdp && sl && ht && cleaning &&TCut("nMuons==1");
  btagSFweight_="probge1"; //use btag sf
  var="muonpt1"; xtitle="muon pT";
  setLogY(true);
  nbins = 20; low=0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_mupt_"+desc,0);
  setLogY(false);


}

void RA2bControlPlots2012(bool tightmet, bool tightht) {


  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  loadSamples(true,"ra2b2012");
  usePUweight_=true;

  useHTeff_=false;
  useMHTeff_=false;    
  thebnnMHTeffMode_ = kOff;
  currentConfig_=configDescriptions_.getDefault();

  resetSamples();
  removeSample("LM9");

  //  removeSample("VV");
  //  removeSample("ZJets");
  addSample("TTV",kGreen-9,"TTV");

  setSampleScaleFactor("TTbarJets",0.85); //made-up factor
  setSampleScaleFactor("WJets",0.84); //remove the k-factor
  setSampleScaleFactor("PythiaPUQCD",1.8);

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=true;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2;
  dodata_=true;

  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40";
  TCut htloose = "HT>=400";
  TCut httight = "HT>=800";
  TCut metloose = "MET>=125";
  TCut mettight = "MET>=250";
  TCut mdp = "minDeltaPhiN_asin>=4";
  TCut ldp = "minDeltaPhiN_asin<4";
  TCut btag = "nbjets>=1";
  TCut bveto = "nbjets==0";
  TCut sl = getSingleLeptonCut();
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");
  TCut withisotk = TCut("nIsoTracks15_005_03>=1");

  TCut met; TString desc2;
  TCut ht;

  if (!tightmet) {
    met = metloose;
    desc2="";
  }
  else {
    met = mettight;
    desc2="_met250";
  }

  if (!tightht) {
    ht = htloose;
    //    desc2+="";
  }
  else {
    ht = httight;
    desc2+="_ht800";
  }

  TString desc = "ge1b_SL"+desc2;

  setDatasetToDraw("2012hybridplus");
  lumiScale_=12034.14 + 5580; //Run 2012 ABC+D

  //compare triggers
//   setDatasetToDraw("HTMHT");
//   selection_ =TCut("HT>=400 && cutPV==1 && (pass_PFHT350_PFMET100==1||pass_PFNoPUHT350_PFMET100==1)  && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1&&buggyEvent==0 && nbjets>=1")&&getSingleLeptonCut();
//   var="MET"; xtitle="E_{T}^{miss} [GeV]";
//   nbins = 30; low=100; high=400;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_"+desc+"_trigHTMHT",0,"GeV");

//   setDatasetToDraw("2012hybrid");
//   selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1  && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1&&buggyEvent==0 && nbjets>=1")&&getSingleLeptonCut();
//   var="MET"; xtitle="E_{T}^{miss} [GeV]";
//   nbins = 30; low=100; high=400;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_"+desc+"_trigHybrid",0,"GeV");

//   //now tight HT
//   selection_ =TCut("HT>=800 && cutPV==1 && cutTrigger==1  && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1&&buggyEvent==0 && nbjets>=1")&&getSingleLeptonCut();
//   var="MET"; xtitle="E_{T}^{miss} [GeV]";
//   nbins = 30; low=100; high=400;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_"+desc+"_trigHybrid_ht800",0,"GeV");

//   setDatasetToDraw("2012hybridplus");
//   selection_ =TCut("HT>=800 && cutPV==1 && cutTrigger2==1  && njets>=3  && minDeltaPhiN >= 4 && passCleaning==1&&buggyEvent==0 && nbjets>=1")&&getSingleLeptonCut();
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_"+desc+"_trigHybridPlus_ht800",0,"GeV");

//  btagSFweight_=""; //reset


  // == ttbar control sample --

  btagSFweight_="probge1"; //use btag sf

  //met
  selection_ = baseline&&cleaning  && mdp && sl &&ht;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 30; low=100; high=400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_"+desc,0,"GeV");

  //met significance (after MET cut)
  selection_ = baseline&&cleaning  && mdp && sl &&ht &&met;
  var="METsig"; xtitle="E_{T}^{miss} significance";
  nbins = 40; low=0; high=400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metsig_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metsig_"+desc,0,"GeV");


  //HT
  selection_ = baseline&&cleaning && met && mdp && sl ;
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 22; low=300; high=1400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_"+desc,0,"GeV");

  //meff
  selection_ = baseline&&cleaning&& met && mdp && sl && ht;
  var="HT+MET"; xtitle="m_{eff} [GeV]";
  nbins = 22; low=400; high=2000;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_meff_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_meff_"+desc,0,"GeV");


  //nPV -- just as a check
  selection_ = baseline&&cleaning && met && mdp && sl && ht;
  var="nGoodPV"; xtitle="num of good PV";
  nbins = 20; low=0; high=40;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_npv_"+desc,0);

  //njets
  selection_ = TCut("cutPV==1 && cutTrigger2==1 && jetpt2>=70") && met && mdp && sl && ht && cleaning;
  var="njets"; xtitle="n jets (pT>50 GeV)";
  nbins = 10; low=0; high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets_"+desc,0);

  //njets30
  selection_ = TCut("cutPV==1 && cutTrigger2==1 && jetpt2>=70") && met && mdp && sl && ht&&cleaning;
  var="njets30"; xtitle="n jets (pT>30 GeV)";
  nbins = 10; low=0; high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets30_"+desc,0);

  //nbjets
  btagSFweight_="";
  selection_ = baseline&&cleaning && met && mdp && sl && ht;
  var="nbjets"; xtitle="n CSVM b jets (pT>50 GeV)";
  nbins = 5; low=0; high=5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_nbjets_"+desc,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_nbjets_"+desc,0);
  setLogY(false);

  btagSFweight_="probge1";

  //minDeltaPhiN
  selection_ = baseline&&cleaning && met && sl && ht;
  var="minDeltaPhiN_asin"; xtitle="minDeltaPhiN (asin)";
  nbins = 20; low=0; high=40;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhiN_"+desc,0);

  //muon pT
  selection_ = baseline&&cleaning && met && getSingleMuonCut() && ht &&mdp;
  var="muonpt1"; xtitle="muon pT";
  nbins = 19; low=10; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_muonpt_"+desc+"_1m",0);

  //muon eta
  selection_ = baseline&&cleaning && met && getSingleMuonCut() && ht &&mdp;
  var="muoneta1"; xtitle="muon #eta";
  nbins = 20; low=-2.5; high=2.5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_muoneta_"+desc+"_1m",0);

  //muon charge
  selection_ = baseline&&cleaning && met && getSingleMuonCut() && ht &&mdp;
  var="muoncharge1"; xtitle="muon charge";
  nbins = 2; low=-1; high=1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_muoncharge_"+desc+"_1m",0);

  //electron charge
  selection_ = baseline&&cleaning && met && getSingleElectronCut() && ht &&mdp;
  var="elecharge1"; xtitle="electron charge";
  nbins = 2; low=-1; high=1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_electroncharge_"+desc+"_1e",0);


  //electron pT
  selection_ = baseline&&cleaning && met && getSingleElectronCut() && ht &&mdp;
  var="eleet1"; xtitle="electron pT";
  nbins = 19; low=10; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_electronpt_"+desc+"_1e",0);

  //electron eta
  selection_ = baseline&&cleaning && met && getSingleElectronCut() && ht &&mdp;
  var="eleeta1"; xtitle="electron #eta";
  nbins = 20; low=-2.5; high=2.5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_electroneta_"+desc+"_1e",0);

  //now the fun stuff

  //jet pt1
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "jetpt1"; xtitle="lead jet pT [GeV]";
  nbins = 20; low = 50; high = 800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetpt1_"+desc,0,"GeV");

  maxScaleFactor_=1.3;
  //jet phi1
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "jetphi1"; xtitle="lead jet phi";
  nbins = 25; low = -TMath::Pi(); high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetphi1_"+desc,0,"GeV");
  //jet eta1

  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "jeteta1"; xtitle="lead jet eta";
  nbins = 30; low = -3; high = 3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jeteta1_"+desc,0,"GeV");
  maxScaleFactor_=1;
  //jet pt2
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "jetpt2"; xtitle="2nd jet pT [GeV]";
  nbins = 20; low = 50; high = 500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetpt2_"+desc,0,"GeV");

  //jet pt3
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "jetpt3"; xtitle="3rd jet pT [GeV]";
  nbins = 20; low = 50; high = 500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetpt3_"+desc,0,"GeV");

  //met phi
  maxScaleFactor_=1.3;
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "METphi"; xtitle="phi of MET";
  nbins = 30; low =-TMath::Pi(); high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metphi_"+desc,0);
  maxScaleFactor_=1.05;

  //unclustered met check
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "unclusteredMET/MET"; xtitle="unclustered MET fraction";
  nbins = 20; low =0; high = 0.4;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_unclusteredMETfrac_"+desc,0);

  maxScaleFactor_=1.3;
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "unclusteredMETphi"; xtitle="phi of unclustered MET";
  nbins = 30; low =-TMath::Pi(); high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_unclusteredmetphi_"+desc,0);
  maxScaleFactor_=1.05;


  //calo met check
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "MET/caloMET"; xtitle="TypeIPFMET/caloMET";
  nbins = 20; low =0; high = 2;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metOverCalomet_"+desc,0);

  //raw met check
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "MET/rawPFMET"; xtitle="TypeIPFMET/rawPFMET";
  nbins = 20; low =0.5; high = 1.5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metOverRawmet_"+desc,0);

  //b jet pt1
  btagSFweight_="";
  selection_ = baseline&&cleaning && met && sl && ht && mdp && btag;
  var = "bjetpt1"; xtitle="lead b jet pT [GeV]";
  nbins = 20; low = 30; high = 500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_bjetpt1_"+desc,0,"GeV");

  btagSFweight_="probge1";

  //jet properties
  maxScaleFactor_=1.3;
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "jetchargedhadronfrac1"; xtitle="lead jet ch had frac";
  nbins = 20; low = 0; high = 1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jet1chhadfrac_"+desc,0);
  maxScaleFactor_=1.05;

  selection_ = baseline&&cleaning && met && sl && ht && mdp;
  var = "jetchargedhadronmult1"; xtitle="lead jet ch had mult";
  nbins = 25; low = 0; high = 50;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jet1chhadmult_"+desc,0);

  selection_ = baseline&& TCut("passCleaning==1&&buggyEvent==0&&MET/caloMET<2") && met && sl && ht && mdp;
  var = "maxTOBTECjetDeltaMult"; xtitle="max TOB/TEC #Delta Multiplicity";
  nbins = 60; low = -100; high = 50;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_maxtobtecjetdeltamult_"+desc,0);


  btagSFweight_="";

  selection_ = baseline&&cleaning && met && sl && ht && mdp && btag;
  var = "bjetchargedhadronfrac1"; xtitle="lead b jet ch had frac";
  nbins = 20; low = 0; high = 1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_bjet1chhadfrac_"+desc,0);

  selection_ = baseline&&cleaning && met && sl && ht && mdp && btag;
  var = "bjetchargedhadronmult1"; xtitle="lead b jet ch had mult";
  nbins = 25; low = 0; high = 50;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_bjet1chhadmult_"+desc,0);

  //MT_b
  selection_ = baseline&&cleaning && met && sl && ht && mdp && btag;
  var = "MT_b"; xtitle="b-jet transverse mass";
  nbins = 40; low = 0; high = 800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_MTb_"+desc,0);

  //regular MT
  btagSFweight_="probge1";
  selection_ = baseline&&cleaning && met && TCut("nMuons+nElectrons==1") && ht && mdp ;
  var = "MT_Wlep"; xtitle="transverse mass";
  nbins = 20; low = 0; high = 200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_MT_"+desc,0);

  //minDeltaPhi
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var="minDeltaPhiN_asin"; xtitle="minDeltaPhiN (asin)";
  nbins = 20; low = 0; high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhi_"+desc,0);

  selection_ = baseline&&cleaning && met && sl && ht ;
  var = "minDeltaPhi"; xtitle="minDeltaPhi (no mdpN cut)";
  nbins = 20; low = 0; high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhi_nominDeltaPhiNcut_"+desc,0);

  //deltaphi*
  selection_ = baseline&&cleaning && met && sl && ht && mdp;
  var = "deltaPhiStar"; xtitle="#Delta#phi *";
  nbins = 30; low = 0; high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_deltaPhiStar_"+desc,0);

  btagSFweight_="";

  selection_ = baseline&&cleaning && met && sl && ht && mdp && btag;
  var = "minDeltaPhiAllb30"; xtitle="minDeltaPhiAllb30";
  nbins = 30; low = 0; high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhiAllb30_"+desc,0);

  btagSFweight_="";//reset

  //iso track sample
  desc="ge0b_ge1isotk"+desc2;
  // pT
  selection_ = baseline&&cleaning && met && withisotk && ht &&mdp;
  var="isotrackpt1"; xtitle="iso track pT";
  nbins = 25; low=0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_isotkpt_"+desc,0,"GeV");
  // eta
  selection_ = baseline&&cleaning && met && withisotk && ht &&mdp;
  var="isotracketa1"; xtitle="iso track eta";
  nbins = 25; low=-2.5; high=2.5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_isotketa_"+desc,0);
  // phi
  selection_ = baseline&&cleaning && met && withisotk && ht &&mdp;
  var="isotrackphi1"; xtitle="iso track phi";
  nbins = 30; low=-TMath::Pi(); high=TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_isotkphi_"+desc,0);


  // == W control sample --> nb == 0, 1 lepton

  btagSFweight_="prob0"; //use btag sf

  desc = "eq0b_SL"+desc2;

  //MET
  selection_ = baseline&&cleaning && met && mdp && sl ;
  var="MET"; xtitle="MET [GeV]";
  nbins = 20; low=100; high=300;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_"+desc,0,"GeV");

  //met significance (after MET cut)
  selection_ = baseline&&cleaning  && mdp && sl &&ht &&met;
  var="METsig"; xtitle="E_{T}^{miss} significance";
  nbins = 40; low=0; high=400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metsig_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metsig_"+desc,0,"GeV");

  //HT
  selection_ = baseline&&cleaning && met && mdp && sl ;
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 22; low=300; high=1400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_"+desc,0,"GeV");

  //meff
  selection_ = baseline&&cleaning && met && mdp && sl && ht;
  var="HT+MET"; xtitle="m_{eff} [GeV]";
  nbins = 22; low=400; high=2000;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_meff_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_meff_"+desc,0,"GeV");


  //nPV -- just as a check
  selection_ = baseline&&cleaning && met && mdp && sl && ht;
  var="nGoodPV"; xtitle="num of good PV";
  nbins = 20; low=0; high=40;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_npv_"+desc,0);

  //njets
  selection_ = TCut("cutPV==1 && cutTrigger2==1  && jetpt2>=70") && met && mdp && sl && ht&&cleaning;
  var="njets"; xtitle="n jets (pT>50 GeV)";
  nbins = 10; low=0; high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets_"+desc,0);

  //njets30
  selection_ = TCut("cutPV==1 && cutTrigger2==1  && jetpt2>=70") && met && mdp && sl && ht&&cleaning;
  var="njets30"; xtitle="n jets (pT>30 GeV)";
  nbins = 10; low=0; high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets30_"+desc,0);


  //minDeltaPhiN
  selection_ = baseline&&cleaning && met && sl && ht;
  var="minDeltaPhiN_asin"; xtitle="minDeltaPhiN (asin)";
  nbins = 20; low=0; high=40;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhiN_"+desc,0);

  //muon pT
  selection_ = baseline&&cleaning && met && getSingleMuonCut() && ht &&mdp;
  var="muonpt1"; xtitle="muon pT";
  nbins = 19; low=10; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_muonpt_"+desc+"_1m",0);

  //muon eta
  selection_ = baseline&&cleaning && met && getSingleMuonCut() && ht &&mdp;
  var="muoneta1"; xtitle="muon #eta";
  nbins = 20; low=-2.5; high=2.5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_muoneta_"+desc+"_1m",0);

  //electron pT
  selection_ = baseline&&cleaning && met && getSingleElectronCut() && ht &&mdp;
  var="eleet1"; xtitle="electron pT";
  nbins = 19; low=10; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_electronpt_"+desc+"_1e",0);

  selection_ = baseline&&cleaning && met && getSingleElectronCut() && ht &&mdp;
  var="eleeta1"; xtitle="electron #eta";
  nbins = 20; low=-2.5; high=2.5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_electroneta_"+desc+"_1e",0);


  //muon charge
  selection_ = baseline&&cleaning && met && getSingleMuonCut() && ht &&mdp;
  var="muoncharge1"; xtitle="muon charge";
  nbins = 2; low=-1; high=1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_muoncharge_"+desc+"_1m",0);

  //electron charge
  selection_ = baseline&&cleaning && met && getSingleElectronCut() && ht &&mdp;
  var="elecharge1"; xtitle="electron charge";
  nbins = 2; low=-1; high=1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_electroncharge_"+desc+"_1e",0);



  //now the fun stuff

  //jet pt1
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "jetpt1"; xtitle="lead jet pT [GeV]";
  nbins = 20; low = 50; high = 800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetpt1_"+desc,0,"GeV");

  maxScaleFactor_=1.3;
  //jet phi1
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "jetphi1"; xtitle="lead jet phi";
  nbins = 25; low = -TMath::Pi(); high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetphi1_"+desc,0,"GeV");
  //jet eta1


  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "jeteta1"; xtitle="lead jet eta";
  nbins = 30; low = -3; high = 3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jeteta1_"+desc,0,"GeV");
  maxScaleFactor_=1;

  //jet pt2
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "jetpt2"; xtitle="2nd jet pT [GeV]";
  nbins = 20; low = 50; high = 500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetpt2_"+desc,0,"GeV");

  //jet pt3
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "jetpt3"; xtitle="3rd jet pT [GeV]";
  nbins = 20; low = 50; high = 400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetpt3_"+desc,0,"GeV");

  //met phi
  maxScaleFactor_=1.3;
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "METphi"; xtitle="phi of MET";
  nbins = 30; low =-TMath::Pi(); high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metphi_"+desc,0);
  maxScaleFactor_=1.05;

  //calo met check
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "MET/caloMET"; xtitle="TypeIPFMET/caloMET";
  nbins = 20; low =0; high = 2;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metOverCalomet_"+desc,0);

  //raw met check
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "MET/rawPFMET"; xtitle="TypeIPFMET/rawPFMET";
  nbins = 20; low =0.5; high = 1.5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metOverRawmet_"+desc,0);

  //b jet pt1
//   btagSFweight_="";
//   selection_ = baseline&&cleaning && met && sl && ht && mdp && btag;
//   var = "bjetpt1"; xtitle="lead b jet pT [GeV]";
//   nbins = 20; low = 30; high = 500;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_bjetpt1_"+desc,0,"GeV");
//   btagSFweight_="probge1";

  //jet properties
  maxScaleFactor_=1.3;
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "jetchargedhadronfrac1"; xtitle="lead jet ch had frac";
  nbins = 20; low = 0; high = 1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jet1chhadfrac_"+desc,0);
  maxScaleFactor_=1.05;

  selection_ = baseline&&cleaning && met && sl && ht && mdp;
  var = "jetchargedhadronmult1"; xtitle="lead jet ch had mult";
  nbins = 25; low = 0; high = 50;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jet1chhadmult_"+desc,0);

//   btagSFweight_="";
//   selection_ = baseline&&cleaning && met && sl && ht && mdp && btag;
//   var = "bjetchargedhadronfrac1"; xtitle="lead b jet ch had frac";
//   nbins = 20; low = 0; high = 1;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_bjet1chhadfrac_"+desc,0);

//   selection_ = baseline&&cleaning && met && sl && ht && mdp && btag;
//   var = "bjetchargedhadronmult1"; xtitle="lead b jet ch had mult";
//   nbins = 20; low = 0; high = 20;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_bjet1chhadmult_"+desc,0);

//   //MT_b
//   selection_ = baseline&&cleaning && met && sl && ht && mdp && btag;
//   var = "MT_b"; xtitle="b-jet transverse mass";
//   nbins = 40; low = 0; high = 400;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_MTb_"+desc,0);

  //regular MT
  //  btagSFweight_="probge1";
  selection_ = baseline&&cleaning && met && TCut("nMuons+nElectrons==1") && ht && mdp ;
  var = "MT_Wlep"; xtitle="transverse mass";
  nbins = 20; low = 0; high = 200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_MT_"+desc,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_MT_"+desc,0);
  setLogY(false);

  //minDeltaPhi
  maxScaleFactor_=1.3;
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "minDeltaPhi"; xtitle="minDeltaPhi";
  nbins = 20; low = 0; high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhi_"+desc,0);
  maxScaleFactor_=1.05;

  maxScaleFactor_=1.3;
  selection_ = baseline&&cleaning && met && sl && ht ;
  var = "minDeltaPhi"; xtitle="minDeltaPhi (no mdpN cut)";
  nbins = 20; low = 0; high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhi_nominDeltaPhiNcut_"+desc,0);
  maxScaleFactor_=1.05;

  //deltaphi*
  selection_ = baseline&&cleaning && met && sl && ht && mdp;
  var = "deltaPhiStar"; xtitle="#Delta#phi *";
  nbins = 30; low = 0; high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_deltaPhiStar_"+desc,0);

//   btagSFweight_="";
//   selection_ = baseline&&cleaning && met && sl && ht && mdp && btag;
//   var = "minDeltaPhiAllb30"; xtitle="minDeltaPhiAllb30";
//   nbins = 30; low = 0; high = TMath::Pi();
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhiAllb30_"+desc,0);
   btagSFweight_="";//reset


   // W+bb region? njets == 2, ge1b

  btagSFweight_="probge1"; //use btag sf
  desc = "ge1b_SL_dijet"+desc2;

  baseline = "cutPV==1 && cutTrigger2==1 && njets==2 && jetpt2>=70";

  //MET
  selection_ = baseline&&cleaning && met && mdp && sl &&ht;
  var="MET"; xtitle="MET [GeV]";
  nbins = 30; low=100; high=400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_"+desc,0,"GeV");

  //met significance (after MET cut)
  selection_ = baseline&&cleaning  && mdp && sl &&ht &&met;
  var="METsig"; xtitle="E_{T}^{miss} significance";
  nbins = 30; low=0; high=400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metsig_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metsig_"+desc,0,"GeV");

  //HT
  selection_ = baseline&&cleaning && met && mdp && sl ;
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 17; low=300; high=2000;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_"+desc,0,"GeV");

  //meff
  selection_ = baseline&&cleaning && met && mdp && sl && ht;
  var="HT+MET"; xtitle="m_{eff} [GeV]";
  nbins = 20; low=400; high=2400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_meff_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_meff_"+desc,0,"GeV");


  //jet pt1
  selection_ = baseline&&cleaning && met && mdp && sl && ht;
  var = "jetpt1"; xtitle="lead jet pT [GeV]";
  nbins = 20; low = 50; high = 800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetpt1_"+desc,0,"GeV");

  maxScaleFactor_=1.3;
  //jet phi1
  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "jetphi1"; xtitle="lead jet phi";
  nbins = 25; low = -TMath::Pi(); high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetphi1_"+desc,0,"GeV");
  //jet eta1

  selection_ = baseline&&cleaning && met && sl && ht && mdp ;
  var = "jeteta1"; xtitle="lead jet eta";
  nbins = 30; low = -3; high = 3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jeteta1_"+desc,0,"GeV");
  maxScaleFactor_=1;


  //muon charge
  selection_ = baseline&&cleaning && met && getSingleMuonCut() && ht &&mdp;
  var="muoncharge1"; xtitle="muon charge";
  nbins = 2; low=-1; high=1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_muoncharge_"+desc+"_1m",0);

  //electron charge
  selection_ = baseline&&cleaning && met && getSingleElectronCut() && ht &&mdp;
  var="elecharge1"; xtitle="electron charge";
  nbins = 2; low=-1; high=1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_electroncharge_"+desc+"_1e",0);

  baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70"; //reset baseline


  // DY sample but using nominal trigger?

//   desc = "dilepton"+desc2;

//   btagSFweight_="prob0"; //b veto
//   selection_ = baseline&&trg && met && mdp && ht&&TCut("bestZmass>10");
//   var="bestZmass"; xtitle="m_{l+l-} [GeV]";
//   nbins = 19; low=10; high=200;
//   setLogY(false);
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_zmass_"+desc,0,"GeV");
//   setLogY(true);

//   desc="dilepton_ge1b"+desc2;
//   btagSFweight_="probge1"; //b veto
//   selection_ = baseline&&trg && met && mdp && ht&&TCut("bestZmass>10");
//   var="bestZmass"; xtitle="m_{l+l-} [GeV]";
//   nbins = 19; low=10; high=200;
//   setLogY(false);
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_zmass_"+desc,0,"GeV");
//   setLogY(true);

  //QCD control sample (with b)

  btagSFweight_="probge1"; //use btag sf
  desc="ge1b_LDP"+desc2;

  //MET
  selection_ = baseline&&cleaning && ldp && zl &&ht;
  var="MET"; xtitle="MET [GeV]";
  nbins = 30; low=100; high=400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_"+desc,0,"GeV");

  //met significance (after MET cut)
  selection_ = baseline&&cleaning  && ldp && zl &&ht &&met;
  var="METsig"; xtitle="E_{T}^{miss} significance";
  nbins = 40; low=0; high=400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metsig_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metsig_"+desc,0,"GeV");

  //HT
  selection_ = baseline&&cleaning && met && ldp && zl ;
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 17; low=300; high=2000;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_"+desc,0,"GeV");

  //meff
  selection_ = baseline&&cleaning && met && ldp && zl && ht;
  var="HT+MET"; xtitle="m_{eff} [GeV]";
  nbins = 20; low=400; high=2400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_meff_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_meff_"+desc,0,"GeV");


  //nPV -- just as a check
  selection_ = baseline&&cleaning && met && ldp && zl && ht;
  var="nGoodPV"; xtitle="num of good PV";
  nbins = 20; low=0; high=40;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_npv_"+desc,0);

  //njets
  selection_ = TCut("cutPV==1 && cutTrigger2==1 && passCleaning==1&&buggyEvent==0 && jetpt2>=70") && met && ldp && zl && ht;
  var="njets"; xtitle="n jets (pT>50 GeV)";
  nbins = 10; low=0; high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets_"+desc,0);

  //njets30
  selection_ = TCut("cutPV==1 && cutTrigger2==1 && passCleaning==1&&buggyEvent==0 && jetpt2>=70") && met && ldp && zl && ht;
  var="njets30"; xtitle="n jets (pT>30 GeV)";
  nbins = 10; low=0; high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets30_"+desc,0);

  //nbjets
  btagSFweight_="";
  selection_ = baseline&&cleaning && met && ldp && zl && ht;
  var="nbjets"; xtitle="n CSVM b jets (pT>50 GeV)";
  nbins = 5; low=0; high=5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_nbjets_"+desc,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_nbjets_"+desc,0);
  setLogY(false);

  btagSFweight_="probge1";

  //minDeltaPhiN
//   selection_ = baseline&&cleaning && met && zl && ht;
//   var="minDeltaPhiN"; xtitle="minDeltaPhiN";
//   nbins = 20; low=0; high=40;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhiN_"+desc,0);

  //now the fun stuff

  //jet pt1
  selection_ = baseline&&cleaning && met && zl && ht && ldp ;
  var = "jetpt1"; xtitle="lead jet pT [GeV]";
  nbins = 20; low = 50; high = 800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetpt1_"+desc,0,"GeV");

  maxScaleFactor_=1.3;
  //jet phi1
  selection_ = baseline&&cleaning && met && zl && ht && ldp ;
  var = "jetphi1"; xtitle="lead jet phi";
  nbins = 25; low = -TMath::Pi(); high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetphi1_"+desc,0,"GeV");
  //jet eta1
  selection_ = baseline&&cleaning && met && zl && ht && ldp ;
  var = "jeteta1"; xtitle="lead jet eta";
  nbins = 30; low = -3; high = 3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jeteta1_"+desc,0,"GeV");
  maxScaleFactor_=1;

  //jet pt2
  selection_ = baseline&&cleaning && met && zl && ht && ldp ;
  var = "jetpt2"; xtitle="2nd jet pT [GeV]";
  nbins = 20; low = 50; high = 600;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetpt2_"+desc,0,"GeV");

  //jet pt3
  selection_ = baseline&&cleaning && met && zl && ht && ldp ;
  var = "jetpt3"; xtitle="3rd jet pT [GeV]";
  nbins = 20; low = 50; high = 500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetpt3_"+desc,0,"GeV");

  //met phi
  maxScaleFactor_=1.3;
  selection_ = baseline&&cleaning && met && zl && ht && ldp ;
  var = "METphi"; xtitle="phi of MET";
  nbins = 30; low =-TMath::Pi(); high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metphi_"+desc,0);
  maxScaleFactor_=1.05;

  //calo met check
  selection_ = baseline&&cleaning && met && zl && ht && ldp ;
  var = "MET/caloMET"; xtitle="TypeIPFMET/caloMET";
  nbins = 20; low =0; high = 2;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metOverCalomet_"+desc,0);

  //raw met check
  selection_ = baseline&&cleaning && met && zl && ht && ldp ;
  var = "MET/rawPFMET"; xtitle="TypeIPFMET/rawPFMET";
  nbins = 20; low =0.5; high = 1.5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metOverRawmet_"+desc,0);

  //b jet pt1
  btagSFweight_="";
  selection_ = baseline&&cleaning && met && zl && ht && ldp && btag;
  var = "bjetpt1"; xtitle="lead b jet pT [GeV]";
  nbins = 20; low = 30; high = 700;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_bjetpt1_"+desc,0,"GeV");

  btagSFweight_="probge1";

  //jet properties
  selection_ = baseline&&cleaning && met && zl && ht && ldp ;
  var = "jetchargedhadronfrac1"; xtitle="lead jet ch had frac";
  nbins = 20; low = 0; high = 1;
  maxScaleFactor_=1.3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jet1chhadfrac_"+desc,0);
  maxScaleFactor_=1.05;

  selection_ = baseline&&cleaning && met && zl && ht && ldp;
  var = "jetchargedhadronmult1"; xtitle="lead jet ch had mult";
  nbins = 25; low = 0; high = 50;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jet1chhadmult_"+desc,0);

  btagSFweight_="";

  selection_ = baseline&&cleaning && met && zl && ht && ldp && btag;
  var = "bjetchargedhadronfrac1"; xtitle="lead b jet ch had frac";
  nbins = 20; low = 0; high = 1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_bjet1chhadfrac_"+desc,0);

  selection_ = baseline&&cleaning && met && zl && ht && ldp && btag;
  var = "bjetchargedhadronmult1"; xtitle="lead b jet ch had mult";
  nbins = 25; low = 0; high = 50;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_bjet1chhadmult_"+desc,0);

  //MT_b
  selection_ = baseline&&cleaning && met && zl && ht && ldp && btag;
  var = "MT_b"; xtitle="b-jet transverse mass";
  nbins = 40; low = 0; high = 400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_MTb_"+desc,0);

  //minDeltaPhi
  selection_ = baseline&&cleaning && met && zl && ht && ldp ;
  var = "minDeltaPhi"; xtitle="minDeltaPhi";
  nbins = 20; low = 0; high = TMath::Pi();
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhi_"+desc,0);
  setLogY(false);

  //BLIND
//   selection_ = baseline&&cleaning && met && sl && ht ;
//   var = "minDeltaPhi"; xtitle="minDeltaPhi (no mdpN cut)";
//   nbins = 20; low = 0; high = TMath::Pi();
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhi_nominDeltaPhiNcut_"+desc,0);

  //deltaphi*
  selection_ = baseline&&cleaning && met && zl && ht && ldp;
  var = "deltaPhiStar"; xtitle="#Delta#phi *";
  nbins = 30; low = 0; high = TMath::Pi();
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_deltaPhiStar_"+desc,0);
  setLogY(false);

  btagSFweight_="";

  selection_ = baseline&&cleaning && met && zl && ht && ldp && btag;
  var = "minDeltaPhiAllb30"; xtitle="minDeltaPhiAllb30";
  nbins = 30; low = 0; high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhiAllb30_"+desc,0);

  btagSFweight_="";//reset


  //signal sample, but with b-veto

  btagSFweight_="prob0"; //use btag sf
  desc="eq0b_ZL"+desc2;

  //MET
  selection_ = baseline&&cleaning && mdp && zl &&ht;
  var="MET"; xtitle="MET [GeV]";
  nbins = 30; low=100; high=400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_met_"+desc,0,"GeV");

  //met significance (after MET cut)
  selection_ = baseline&&cleaning  && mdp && zl &&ht &&met;
  var="METsig"; xtitle="E_{T}^{miss} significance";
  nbins = 40; low=0; high=400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metsig_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metsig_"+desc,0,"GeV");

  //HT
  selection_ = baseline&&cleaning && met && mdp && zl ;
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 17; low=300; high=2000;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht_"+desc,0,"GeV");

  //meff
  selection_ = baseline&&cleaning && met && mdp && zl && ht;
  var="HT+MET"; xtitle="m_{eff} [GeV]";
  nbins = 20; low=400; high=2400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_meff_"+desc,0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_meff_"+desc,0,"GeV");


  //nPV -- just as a check
  selection_ = baseline&&cleaning && met && mdp && zl && ht;
  var="nGoodPV"; xtitle="num of good PV";
  nbins = 20; low=0; high=40;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_npv_"+desc,0);

  //njets
  selection_ = TCut("cutPV==1 && cutTrigger2==1 && jetpt2>=70") && met && mdp && zl && ht &&cleaning;
  var="njets"; xtitle="n jets (pT>50 GeV)";
  nbins = 10; low=0; high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets_"+desc,0);

  //njets30
  selection_ = TCut("cutPV==1 && cutTrigger2==1  && jetpt2>=70") && met && mdp && zl && ht &&cleaning;
  var="njets30"; xtitle="n jets (pT>30 GeV)";
  nbins = 10; low=0; high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_njets30_"+desc,0);

  //nbjets BLIND
//   btagSFweight_="";
//   selection_ = baseline&&cleaning && met && mdp && zl && ht;
//   var="nbjets"; xtitle="n CSVM b jets (pT>50 GeV)";
//   nbins = 5; low=0; high=5;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_nbjets_"+desc,0);
//   setLogY(true);
//   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_nbjets_"+desc,0);
//   setLogY(false);

  btagSFweight_="prob0";

  //minDeltaPhiN
  selection_ = baseline&&cleaning && met && zl && ht;
  var="minDeltaPhiN_asin"; xtitle="minDeltaPhiN (asin)";
  nbins = 20; low=0; high=40;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhiN_"+desc,0);

  //now the fun stuff

  //jet pt1
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "jetpt1"; xtitle="lead jet pT [GeV]";
  nbins = 20; low = 50; high = 800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetpt1_"+desc,0,"GeV");

  maxScaleFactor_=1.3;
  //jet phi1
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "jetphi1"; xtitle="lead jet phi";
  nbins = 25; low = -TMath::Pi(); high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetphi1_"+desc,0);
  //jet eta1
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "jeteta1"; xtitle="lead jet eta";
  nbins = 30; low = -3; high = 3;
  if (tightht) nbins=15;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jeteta1_"+desc,0);
  maxScaleFactor_=1;

  //jet pt2
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "jetpt2"; xtitle="2nd jet pT [GeV]";
  nbins = 20; low = 50; high = 600;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetpt2_"+desc,0,"GeV");

  //jet pt3
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "jetpt3"; xtitle="3rd jet pT [GeV]";
  nbins = 20; low = 50; high = 500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jetpt3_"+desc,0,"GeV");

  //met phi
  maxScaleFactor_=1.3;
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "METphi"; xtitle="phi of MET";
  nbins = 30; low =-TMath::Pi(); high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metphi_"+desc,0);
  maxScaleFactor_=1.05;

  //calo met check
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "MET/caloMET"; xtitle="TypeIPFMET/caloMET";
  nbins = 20; low =0; high = 2;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metOverCalomet_"+desc,0);

  //calo met check (wider range)
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "MET/caloMET"; xtitle="TypeIPFMET/caloMET";
  nbins = 50; low =0; high = 20;
  setLogY(true); setPlotMinimum(0.1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metOverCalometWide_"+desc,0);

  maxScaleFactor_=2;
  //DeltaPhi(j1,MET)
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "deltaPhi1"; xtitle="DeltaPhi(j1,MET)";
  nbins = 25; low =0; high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_deltaphi1_"+desc,0);
  maxScaleFactor_=1.05;

  //this was a special plot for debugging
//    maxScaleFactor_=2;
//   //DeltaPhi(j1,MET) with MET/calotMET cut
//    selection_ = baseline && cleaning && met && zl && ht && mdp&&TCut("MET/caloMET>2") ;
//    var = "deltaPhi1"; xtitle="DeltaPhi(j1,MET)";
//    nbins = 25; low =0; high = TMath::Pi();
//    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_deltaphi1_highMEToverCaloMET_"+desc,0);
//    maxScaleFactor_=1.05;

  //raw met check
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "MET/rawPFMET"; xtitle="TypeIPFMET/rawPFMET";
  nbins = 20; low =0.5; high = 1.5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metOverRawmet_"+desc,0);

  //jet properties
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "jetchargedhadronfrac1"; xtitle="lead jet ch had frac";
  nbins = 20; low = 0; high = 1;
  maxScaleFactor_=1.3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jet1chhadfrac_"+desc,0);
  maxScaleFactor_=1.05;

  selection_ = baseline&&cleaning && met && zl && ht && mdp;
  var = "jetchargedhadronmult1"; xtitle="lead jet ch had mult";
  nbins = 25; low = 0; high = 50;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jet1chhadmult_"+desc,0);

  if (false) {//special plots for investigation...

    selection_ = baseline&&TCut("passCleaning==1&&buggyEvent==0") && met && zl && ht && mdp &&TCut("MET/caloMET>2");
    var = "jetchargedhadronmult1"; xtitle="lead jet ch had mult";
    nbins = 30; low = 0; high = 200;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jet1chhadmult_highMEToverCaloMET_"+desc,0);
    
    selection_ = baseline&&TCut("passCleaning==1&&buggyEvent==0") && met && zl && ht && mdp &&TCut("MET/caloMET>2");
    maxScaleFactor_=1.3;
    var = "METphi"; xtitle="phi of MET";
    nbins = 30; low =-TMath::Pi(); high = TMath::Pi();
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metphi_highMEToverCaloMET_"+desc,0);
    maxScaleFactor_=1.05;

    selection_ = baseline&&TCut("passCleaning==1&&buggyEvent==0") && met && zl && ht && mdp &&TCut("MET/caloMET>2");
    maxScaleFactor_=1.3;
    var = "jetchargedhadronfrac1"; xtitle="lead jet ch had frac";
    nbins = 20; low = 0; high = 1;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jet1chhadfrac_highMEToverCaloMET_"+desc,0);
    maxScaleFactor_=1.05;

    setDatasetToDraw("2012hybrid");
    //calo met check
    //baseline but with old 'cutTrigger' variable
    selection_ = TCut("cutPV==1 && cutTrigger==1 && njets>=3 && passCleaning==1&&buggyEvent==0 && jetpt2>=70") && met && zl && ht && mdp ;
    var = "MET/caloMET"; xtitle="TypeIPFMET/caloMET";
    nbins = 20; low =0; high = 2;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metOverCalomet_HTMHTandMETonly_"+desc,0);


    setDatasetToDraw("2012hybridplus");
    //try after selecting "bad" events
    //met significance (after MET cut)
    selection_ = baseline&&TCut("passCleaning==1&&buggyEvent==0")  && mdp && zl &&ht &&met && TCut("MET/caloMET>2");
    var="METsig"; xtitle="E_{T}^{miss} significance";
    nbins = 40; low=0; high=400;
    setLogY(false);
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metsig_highMEToverCaloMET_"+desc,0);
    setLogY(true);
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metsig_highMEToverCaloMET"+desc,0);

    //flip that around
    selection_ = baseline&&TCut("passCleaning==1&&buggyEvent==0")  && mdp && zl &&ht &&met && TCut("METsig>20");
    //calo met check (wider range)
    var = "MET/caloMET"; xtitle="TypeIPFMET/caloMET";
    nbins = 20; low =0; high = 4;
    setLogY(false);// setPlotMinimum(0.1);
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metOverCalomet_highMETsig_"+desc,0);

    selection_ = baseline&&TCut("passCleaning==1&&buggyEvent==0")  && mdp && zl &&ht &&met && TCut("METsig<20");
    //calo met check (wider range)
    var = "MET/caloMET"; xtitle="TypeIPFMET/caloMET";
    nbins = 20; low =0; high = 4;
    setLogY(false);// setPlotMinimum(0.1);
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_metOverCalomet_lowMETsig_"+desc,0);

    //jet eta1
    selection_ = baseline&&TCut("passCleaning==1&&buggyEvent==0") && met && zl && ht && mdp ;
    var = "jeteta1"; xtitle="lead jet eta";
    nbins = 40; low = -5; high = 5;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jeteta1_NoMetOverCaloMetCut_"+desc,0);
    //maxScaleFactor_=1;

    //with new cleanup cut
    selection_ = baseline&&TCut("passCleaning==1&&buggyEvent==0") && met && zl && ht && mdp && TCut("MET/caloMET<2") ;
    var = "jeteta1"; xtitle="lead jet eta";
    nbins = 40; low = -5; high = 5;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jeteta1_wide_lowMEToverCaloMET_"+desc,0);

    //invert new cleanup cut
    selection_ = baseline&&TCut("passCleaning==1&&buggyEvent==0") && met && zl && ht && mdp && TCut("MET/caloMET>=2") ;
    var = "jeteta1"; xtitle="lead jet eta";
    nbins = 40; low = -5; high = 5;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_jeteta1_wide_highMEToverCaloMET_"+desc,0);

  }

  //minDeltaPhi
  selection_ = baseline&&cleaning && met && zl && ht && mdp ;
  var = "minDeltaPhi"; xtitle="minDeltaPhi";
  nbins = 20; low = 0; high = TMath::Pi();
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhi_"+desc,0);
  setLogY(false);

  selection_ = baseline&&cleaning && met && sl && ht ;
  var = "minDeltaPhi"; xtitle="minDeltaPhi (no mdpN cut)";
  nbins = 20; low = 0; high = TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_minDeltaPhi_nominDeltaPhiNcut_"+desc,0);

  //deltaphi*
  selection_ = baseline&&cleaning && met && zl && ht && mdp;
  var = "deltaPhiStar"; xtitle="#Delta#phi *";
  nbins = 30; low = 0; high = TMath::Pi();
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_deltaPhiStar_"+desc,0);
  setLogY(false);

 
  btagSFweight_="";//reset


  return;

  //old stuff

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

/*
  //try requiring 2 leptons
  setDatasetToDraw("DoubleMu");

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
*/

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
  

  //etm
  selection_ = "HT>=400 && njets>=2 && jetpt1>=70 && jetpt2>=70 && (nElectrons+nMuons==0) &&passCleaning==1";
  drawTrigEff("MET", "pass_L1ETM40==1","pass_PFHT350_PFMET100==1", "MET", 20, 0, 500);

  selection_ = "HT>=400 && njets>=2 && jetpt1>=70 && jetpt2>=70 && (nElectrons+nMuons==0) &&passCleaning==1";
  drawTrigEff("MET", "pass_L1ETM40==1","pass_PFHT350_PFMET100==1||pass_DiCentralPFJet50_PFMET80==1", "MET", 20, 0, 500);


  selection_ = "HT>=400 && njets>=2 && jetpt1>=70 && jetpt2>=70 && nMuons==1 &&passCleaning==1";
  drawTrigEff("MET", "pass_L1ETM40==1","pass_PFHT350_PFMET100==1||pass_DiCentralPFJet50_PFMET80==1", "MET", 20, 0, 500);


  //HT leg (Run2012A+B)
  selection_ = "njets>=2 && jetpt1>=70 && jetpt2>=70 && nMuons==1 &&passCleaning==1 && MET>=125 && runNumber<198021";
  drawTrigEff("SingleMu", "pass_IsoMu24_eta2p1==1","pass_PFHT650==1", "HT", 50, 0, 1000);


}


