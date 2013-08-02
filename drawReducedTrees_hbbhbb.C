/*
to compile:

must have symlink to MiscUtil.cxx in the working directory.

gSystem->Load("TSelectorMultiDraw_C.so");
gSystem->Load("CrossSectionTable_cxx.so");
gSystem->Load("ConfigurationDescriptions_cxx.so");

.L drawReducedTrees_hbbhbb.C+

or

gSystem->Load("drawReducedTrees_hbbhbb_C.so");


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

void initHiggsSamples69(const bool useSkim=true,const TString samplelist="") { 

  //  if (useSkim==false) {cout<<"Must use skim!"<<endl; assert(0);}
  const  TString skimstring = useSkim?"-skim":"";

  addToSamplesAll("TTJets_ResolutionTest_PU_S10_CU0001_v69"); //special test
  addToSamplesAll("TTJets_ResolutionTest_PU_S10_CU0001_JERbias_v69"); // special test
  addToSamplesAll("HiggsinoNLSP_chargino130_to_500_bino1_TChihh_v69");

  //MG ttbar  
  addToSamplesAll(TString("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1800_v69")+skimstring);
  addToSamplesAll(TString("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1799_v69")+skimstring);
  addToSamplesAll(TString("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1848_v69")+skimstring);
  //missing the hadronic ttbar for the moment

  //sherpa
  addToSamplesAll(TString("TTJets_SemiLeptDecays_8TeV-sherpa_Summer12_DR53X-PU_S10_START53_V19-v1_AODSIM_UCSB1806_v69")+skimstring);
  addToSamplesAll(TString("TTJets_DileptDecays_8TeV-sherpa_Summer12_DR53X-PU_S10_START53_V19-v1_AODSIM_UCSB1804_v69")+skimstring);

  //QCD
  addToSamplesAll(TString("QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1814_v69")+skimstring);
  addToSamplesAll(TString("QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1815_v69")+skimstring);
  addToSamplesAll(TString("QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1816_v69")+skimstring);
  addToSamplesAll(TString("QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1817_v69")+skimstring);
  addToSamplesAll(TString("QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1818_v69")+skimstring);
  addToSamplesAll(TString("QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1819_v69")+skimstring);
  addToSamplesAll(TString("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1820_v69")+skimstring);
  addToSamplesAll(TString("QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1821_v69")+skimstring);
  addToSamplesAll(TString("QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1822_v69")+skimstring);

  //Znunu
  addToSamplesAll(TString("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1832_v69")+skimstring);
  addToSamplesAll(TString("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1833_v69")+skimstring);
  addToSamplesAll(TString("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1834_v69")+skimstring);
  addToSamplesAll(TString("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1831_v69")+skimstring);
  addToSamplesAll(TString("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1835_v69")+skimstring);


  //signal (no skim)
  addToSamplesAll("TChihh-HbbHbb_mHiggsino-200_mLSP-1_8TeV-Pythia6Z_jgsmith-SMS-HbbHbb_mHiggsino-200_mLSP-1_1_UCSB1807_v69");
  addToSamplesAll("TChihh-HbbHbb_mHiggsino-350_mLSP-1_8TeV-Pythia6Z_jgsmith-SMS-HbbHbb_mHiggsino-350_mLSP-1_1_UCSB1811_v69");

  addToSamplesAll("SMS-HbbHbb_mHiggsino-200_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1807_v69");
  addToSamplesAll("SMS-HbbHbb_mHiggsino-250_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1809_v69");
  addToSamplesAll("SMS-HbbHbb_mHiggsino-300_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1810_v69");
  addToSamplesAll("SMS-HbbHbb_mHiggsino-350_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1811_v69");
  addToSamplesAll("SMS-HbbHbb_mHiggsino-400_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1812_v69");
  addToSamplesAll("SMS-HbbHbb_mHiggsino-450_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1808_v69");

  nameOfEventWeight_="weight3"; 

  if (samplelist=="metsigtest1") inputPath = "/cu5/joshmt/reducedTrees/v69_2_pt20/";//special test
  else  inputPath="/cu6/joshmt/reducedTrees/v69_2_pt20/"; 
  dataInputPath="/cu6/joshmt/reducedTrees/v69_2_pt20/data/"; //for skimmed data
  if (!useSkim)   dataInputPath+="unskimmed/";

  loadSamples(true,"ra2b2012");
  usePUweight_=true; //helps bring signal and background into alignment

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;
  clearSamples();
  resetChains(); //important in the case that we call initHiggsSamples() more than once

  if (samplelist=="metsigtest1") {
    nameOfEventWeight_="(1)";
    addSample("TTJets_ResolutionTest_PU_S10_CU0001_v69",kAzure-3,"t#bar{t}");

    //manually set scale factor
    setSampleScaleFactor("TTJets_ResolutionTest_PU_S10_CU0001_v69",234.0 / getTree("TTJets_ResolutionTest_PU_S10_CU0001_v69")->GetEntries());
  }

  if (samplelist=="" || samplelist.Contains("qcd")) {
    TString baseqcdname = TString("QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1814_v69")+skimstring;
    addSample(baseqcdname,kGreen,"QCD");
    chainSamples(baseqcdname,TString("QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1815_v69")+skimstring);
    chainSamples(baseqcdname,TString("QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1816_v69")+skimstring);
    chainSamples(baseqcdname,TString("QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1817_v69")+skimstring);
    chainSamples(baseqcdname,TString("QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1818_v69")+skimstring);
    chainSamples(baseqcdname,TString("QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1819_v69")+skimstring);
    chainSamples(baseqcdname,TString("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1820_v69")+skimstring);
    chainSamples(baseqcdname,TString("QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1821_v69")+skimstring);
    chainSamples(baseqcdname,TString("QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1822_v69")+skimstring);

  }

  if (samplelist==""||samplelist.Contains("ttbar")) {
    addSample(TString("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1848_v69")+skimstring,kGreen-9,"t#bar{t} (0l)");
  }

  if (samplelist==""||samplelist.Contains("znunu")) {
    TString basezname = TString("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1832_v69")+skimstring;
    addSample(basezname,kOrange,"Z #rightarrow #nu#nu");
    chainSamples(basezname,TString("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1833_v69")+skimstring);
    chainSamples(basezname,TString("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1834_v69")+skimstring);
    chainSamples(basezname,TString("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1831_v69")+skimstring);
    chainSamples(basezname,TString("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1835_v69")+skimstring);
  }

  if (samplelist==""||samplelist.Contains("ttbar")) {

    addSample(TString("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1799_v69")+skimstring,kBlue-10,"t#bar{t} (2l)");
    addSample(TString("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1800_v69")+skimstring,kAzure-3,"t#bar{t} (1l)");

    //    setSampleWeightFactor(TString("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1800_v69")+skimstring,"topPtWeight");
    //    if (useSkim)    setSampleWeightFactor(TString("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1799_v69")+skimstring,"topPtWeight");
  }


  if (samplelist==""||samplelist.Contains("hh200")) 
    //    addSample("TChihh-HbbHbb_mHiggsino-200_mLSP-1_8TeV-Pythia6Z_jgsmith-SMS-HbbHbb_mHiggsino-200_mLSP-1_1_UCSB1807_v69",kRed-7,"hh 200");
    addSample("SMS-HbbHbb_mHiggsino-200_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1807_v69",kRed-9,"hh 200");
  if (samplelist==""||samplelist.Contains("hh250"))
    addSample("SMS-HbbHbb_mHiggsino-250_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1809_v69",kRed-7,"hh 250");
  if (samplelist==""||samplelist.Contains("hh300"))
    addSample("SMS-HbbHbb_mHiggsino-300_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1810_v69",kRed-4,"hh 300");
  if (samplelist==""||samplelist.Contains("hh350"))
    //    addSample("TChihh-HbbHbb_mHiggsino-350_mLSP-1_8TeV-Pythia6Z_jgsmith-SMS-HbbHbb_mHiggsino-350_mLSP-1_1_UCSB1811_v69",kRed+2,"hh 350");
    addSample("SMS-HbbHbb_mHiggsino-350_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1811_v69",kRed+1,"hh 350");
  if (samplelist==""||samplelist.Contains("hh400"))
    addSample("SMS-HbbHbb_mHiggsino-400_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1812_v69",kRed+2,"hh 400");
  if (samplelist==""||samplelist.Contains("hh450"))
    addSample("SMS-HbbHbb_mHiggsino-450_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1808_v69",kRed+3,"hh 450");


  const float hbbbb=0.561*0.561;
  //https://twiki.cern.ch/twiki/bin/view/CMS/Electrohiggs#Prospino_NLO_Cross_Sections
   //sigma x BF / ngen
  if (samplelist==""||samplelist.Contains("hh200"))
    setSampleScaleFactor("SMS-HbbHbb_mHiggsino-200_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1807_v69", 0.608 *hbbbb/ getTree("SMS-HbbHbb_mHiggsino-200_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1807_v69")->GetEntries());
  if (samplelist==""||samplelist.Contains("hh250"))
    setSampleScaleFactor("SMS-HbbHbb_mHiggsino-250_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1809_v69", 0.244 *hbbbb/ getTree("SMS-HbbHbb_mHiggsino-250_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1809_v69")->GetEntries());
  if (samplelist==""||samplelist.Contains("hh300"))
    setSampleScaleFactor("SMS-HbbHbb_mHiggsino-300_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1810_v69", 0.111  *hbbbb/ getTree("SMS-HbbHbb_mHiggsino-300_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1810_v69")->GetEntries());
  if (samplelist==""||samplelist.Contains("hh350"))
    setSampleScaleFactor("SMS-HbbHbb_mHiggsino-350_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1811_v69", 0.0553 *hbbbb/getTree("SMS-HbbHbb_mHiggsino-350_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1811_v69")->GetEntries());
  if (samplelist==""||samplelist.Contains("hh400"))
    setSampleScaleFactor("SMS-HbbHbb_mHiggsino-400_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1812_v69",  0.0294 *hbbbb/getTree("SMS-HbbHbb_mHiggsino-400_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1812_v69")->GetEntries());
  if (samplelist==""||samplelist.Contains("hh450"))
    setSampleScaleFactor("SMS-HbbHbb_mHiggsino-450_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1808_v69",   0.0163 *hbbbb/getTree("SMS-HbbHbb_mHiggsino-450_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1808_v69")->GetEntries());
  


  if (samplelist.Contains("hhmg200")) {
    addSample("HiggsinoNLSP_chargino130_to_500_bino1_TChihh_v69:m0==200",kBlue-7,"hh MG 200");
    setSampleScaleFactor("HiggsinoNLSP_chargino130_to_500_bino1_TChihh_v69:m0==200", 0.608 *hbbbb/408456.0); //sigma x BF / ngen
  }

  if (samplelist.Contains("hhmg350")) {
    addSample("HiggsinoNLSP_chargino130_to_500_bino1_TChihh_v69:m0==350",kBlue+2,"hh MG 350");
    setSampleScaleFactor("HiggsinoNLSP_chargino130_to_500_bino1_TChihh_v69:m0==350", 0.0553 *hbbbb/72112.0); //sigma x BF / ngen
  }

  setDatasetToDraw("MET"); //when plotting data, use the MET PD


}

void mgSignalCheck() {
  initHiggsSamples69(true,"hh200 hh350 hhmg200 hhmg350");
  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  usePUweight_=false;  //broken 

  setStackMode(false,false,false); //stack,norm,label override

  selection_ = "(1)"; // no cuts


  nbins=30; low=0; high=300;
  var="MET"; xtitle="MET";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_MET",0);


  nbins=9; low=0; high=9;
  var="njets20"; xtitle="njets20";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_njets20",0);


  nbins=20; low=0; high=100;
  var="jetpt4"; xtitle="jetpt4";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_jetpt4",0);
  setStackMode(false,true,false); //stack,norm,label override
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_jetpt4",0);

  setStackMode(false,false,false); //stack,norm,label override
  nbins=40; low=0; high=200;
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="< m_{jj} >";
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_avmh",0);

   //lower level stuff for debugging
   selection_="jetpt1>1";

   nbins=40; low=-5; high=5;
   var="jeteta1"; xtitle=var;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);

   //lower level stuff for debugging
   nbins=40; low=-4; high=4;
   var="jetphi1"; xtitle=var;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);
  setStackMode(false,true,false); //stack,norm,label override
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);


   selection_="jetpt2>1";
   nbins=40; low=-5; high=5;
   var="jeteta2"; xtitle=var;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);

   nbins=40; low=-4; high=4;
   var="jetphi2"; xtitle=var;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);
  

   selection_="jetpt3>1";
   nbins=40; low=-5; high=5;
   var="jeteta3"; xtitle=var;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);

   nbins=40; low=-4; high=4;
   var="jetphi3"; xtitle=var;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);
  
   selection_="jetpt4>1";
   nbins=40; low=-5; high=5;
   var="jeteta4"; xtitle=var;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);

   nbins=40; low=-4; high=4;
   var="jetphi4"; xtitle=var;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);

   selection_="(1)";
   //gen level
   nbins=40; low=-5; high=5;
   var="higgs1b1eta"; xtitle=var;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);

  nbins=40; low=-5; high=5;
   var="higgs1b2eta"; xtitle=var;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);

    nbins=40; low=-5; high=5;
   var="higgs2b1eta"; xtitle=var;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);

  nbins=40; low=-5; high=5;
   var="higgs2b2eta"; xtitle=var;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);
   //phi
  nbins=40; low=-4; high=4;
   var="higgs1b1phi"; xtitle=var;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);

  nbins=40; low=-4; high=4;
   var="higgs1b2phi"; xtitle=var;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);

    nbins=40; low=-4; high=4;
   var="higgs2b1phi"; xtitle=var;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);

  nbins=40; low=-4; high=4;
   var="higgs2b2phi"; xtitle=var;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);

   //
  nbins=5; low=0; high=5;
  var="nPUjets20"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_"+var,0);

  //now try a cut flow
//define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";
  //cannot use this with MET PD only
  //  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets30==4";
  TCut njets5="njets30==5";
  TCut jet2="jetpt2>50";
  TCut mdp = "minDeltaPhi30>0.3";
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

  TCut metsigveryloose="METsig_2012>25";

  TCut drmax = "deltaRmax_hh<2.2";
  //  TCut drmin = "deltaRmin_hh<1.9";
  selection_ = baseline && trigger && zl&&isotk&& (njets4||njets5)&&jet2&& btag2&&btag3&&btag4 && higgsSR && drmax && TCut("MET>110");

  nbins=20; low=0; high=5;
  var="deltaRmin_hh"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_signalLikeSelection_"+var,0);



  //now, another check
  selection_="(1)";
  clearSamples();
    addSample("HiggsinoNLSP_chargino130_to_500_bino1_TChihh_v69:m0==200",kGreen,"hh MG 200");
    addSample("HiggsinoNLSP_chargino130_to_500_bino1_TChihh_v69:m0==250",kBlue,"hh MG 250");
    addSample("HiggsinoNLSP_chargino130_to_500_bino1_TChihh_v69:m0==300",kMagenta,"hh MG 300");
    addSample("HiggsinoNLSP_chargino130_to_500_bino1_TChihh_v69:m0==350",kRed,"hh MG 350");
  setStackMode(false,true,false); //stack,norm,label override
 
  nbins=40; low=0; high=400;
  var="MET"; xtitle="MET";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidationAll4_"+var,0);

  nbins=40; low=0; high=300;
  var="METsig"; xtitle="METsig";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidationAll4_"+var,0);


  nbins=40; low=0; high=300;
  var="METsig_2012"; xtitle="METsig_2012";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidationAll4_"+var,0);


}

void metsigPUcheck() {

  //show shape of ttbar METsig distribution in bins of PU
  initHiggsSamples69(true,"ttbar");
  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";
  //cannot use this with MET PD only
  //  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets30==4";
  TCut njets5="njets30==5";
  TCut jet2="jetpt2>50";
  TCut mdp = "minDeltaPhi30>0.3";
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

  TCut metsigveryloose="METsig_2012>25";

  TCut drmax = "deltaRmax_hh<2.2";
  //  TCut drmin = "deltaRmin_hh<1.9";

  //  setSampleScaleFactor("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim",0.8);

  setStackMode(false,true,false); //stack,norm,label override
  setLogY(false);
  dodata_=false;
  doRatio_=false;

  //2b pre-selection ; basically the minimum required by the trigger, plus SL
  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btag2;

  clearSamples();
  addSample("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1800_v69-skim:nGoodPV<=10",kViolet,"nPV #leq 10");
  addSample("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1800_v69-skim:nGoodPV<=15&&nGoodPV>10",kBlue,"10 < nPV #leq 15");
  addSample("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1800_v69-skim:nGoodPV<=20&&nGoodPV>15",kGreen,"15 < nPV #leq 20");
  addSample("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1800_v69-skim:nGoodPV>20",kRed,"nPV > 20");

  nbins=30; low=0; high=300;
  var="METsig_2012"; xtitle="MET significance (new)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_BinsOfPU_METsig2012",0);

  nbins=30; low=0; high=300;
  var="METsig"; xtitle="MET significance (old)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_BinsOfPU_METsig",0);

  nbins=30; low=0; high=400;
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_BinsOfPU_MET",0);



}

void higgs_genlevel_signal() {


  initHiggsSamples69(false,"hh200 hh250 hh400");


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

  initHiggsSamples69(false,"hh250 hh400");


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

  initHiggsSamples69(true,"ttbarjoint");


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

void higgs_dataMC_control_QCD_noskim() {
  //plots that are looser than the skim


  initHiggsSamples69(false,"ttbar znunu qcd");

  setOutputDirectory("plots_Control_QCD");

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 3.2;
  dodata_=true;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 
  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut ra2btrigger = "passMC_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80==1 || passMC_DiCentralPFJet50_PFMET80==1";

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut jet2high="jetpt2>70";
  TCut jets=(njets4||njets5) &&jet2;
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut mdp = "minDeltaPhi20>0.3";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;

  TCut metsigbound="METsig<60";

  TCut drmax = "deltaRmax_hh<2.2";

  TCut metsigloose = "METsig>30";

  //plot:
  //DeltaPhi(jets,MET)

  setStackMode(true,false,false); //stack,norm,label override
  setLogY(false);


  //we want to study loose (or no) b-tagging. QCD dominated sample

  //first look at MET-only trigger
  selection_ = baseline && triggerMET && zl&&isotk && jets&&jet2 && !btag3 &&!mdp;

  nbins=40; low=0; high=400;
  setLogY(false);
  var="METsig"; xtitle="MET significance";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigMET150_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigMET150_METsig",0);


  nbins=50; low=0; high=500;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigMET150_MET",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigMET150_MET",0);

  //might be better to raise the MET cut above the turn-on curve
  selection_ = baseline && triggerMET && zl&&isotk && jets&&jet2 && !btag3 &&!mdp &&TCut("MET>230");

  nbins=40; low=0; high=400;
  setLogY(false);
  var="METsig"; xtitle="MET significance";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigMET150-MET230_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigMET150-MET230_METsig",0);


  //for RA2b trigger, raise jetpt2 threshold
  selection_ = baseline && ra2btrigger && zl&&isotk && jets&&jet2high && !btag3 &&!mdp;

  nbins=40; low=0; high=400;
  setLogY(false);
  var="METsig"; xtitle="MET significance";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigJetMET_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigJetMET_METsig",0);


  nbins=50; low=0; high=500;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigJetMET_MET",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigJetMET_MET",0);


  selection_ = baseline && ra2btrigger && zl&&isotk && jets&&jet2high && !btag3 &&!mdp&&TCut("MET>160");

  nbins=40; low=0; high=400;
  setLogY(false);
  var="METsig"; xtitle="MET significance";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigJetMET-MET160_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigJetMET-MET160_METsig",0);

  // == look at anomalous events
  TCut weird = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& (MET/caloMET>=2 || maxTOBTECjetDeltaMult>=40)"; 

  selection_ = weird && ra2btrigger && zl&&isotk && jets&&jet2high && !btag3 &&TCut("MET>160");

  nbins=50; low=0; high=500;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_weird_trigJetMET_MET",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_weird_trigJetMET_MET",0);

  nbins=40; low=0; high=400;
  setLogY(false);
  var="METsig"; xtitle="MET significance";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_weird_trigJetMET_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_weird_trigJetMET_METsig",0);

  TCut junk = "cutPV==1 &&passCleaning==0 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  selection_ = junk && ra2btrigger && zl&&isotk && jets&&jet2high && !btag3 &&TCut("MET>160");

  nbins=50; low=0; high=500;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_junk_trigJetMET_MET",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_junk_trigJetMET_MET",0);

  nbins=40; low=0; high=400;
  setLogY(false);
  var="METsig"; xtitle="MET significance";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_junk_trigJetMET_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_junk_trigJetMET_METsig",0);


  // == now other plots (non-LDP)
  selection_ = baseline && trigger && zl && isotk && jet2&& btag2 && !btag3 && mdp && metsigloose;
  nbins=10; low = 0; high=10;
  setLogY(false);
  var="njets20"; xtitle="jet multiplicity (p_{T} > 20 GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bMETsig30mdp_njets20",0);

  selection_ = baseline && trigger && zl && isotk && jet2&& btag2 && !btag3 && mdp && metsigloose;
  nbins=4; low = 0; high=4;
  setLogY(true);
  var="nPUjets20"; xtitle="pileup jets (p_{T} > 20 GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bMETsig30mdp_nPUjets20",0);




}

void higgs_dataMC_control_QCD() {

  //goal: use 2-b sample to get a qcd control region. 

  initHiggsSamples69(true,"ttbar znunu qcd");

  setOutputDirectory("plots_Control_QCD");

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut jets=(njets4||njets5) &&jet2;
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut mdp = "minDeltaPhi20>0.3";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;

  TCut metsigbound="METsig<60";
  TCut metsigloose = "METsig>30";

  TCut drmax = "deltaRmax_hh<2.2";


  //plot:
  //DeltaPhi(jets,MET)

  setStackMode(true,false,false); //stack,norm,label override
  setLogY(false);

   //i think we want: preselection. 2b tags. lepton & isotk vetoes.
  selection_ = baseline && trigger && zl&&isotk && jets && btag2 && !btag3;

  /*
  nbins=30; low=0; high=300;
  setLogY(false);
  var="METsig_2012"; xtitle="MET significance (new)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselection_METsig2012",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselection_METsig2012",0);
  */

  nbins=30; low=0; high=300;
  setLogY(false);
  var="METsig"; xtitle="MET significance";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselection_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselection_METsig",0);

  nbins=30; low=0; high=400;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselection_MET",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselection_MET",0);

  //nPV
  nbins=20; low=0; high=40;
  setLogY(false);
  var="nGoodPV"; xtitle="n Good PV";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselection_nGoodPV",0);

  //same plots but adding mdp cut
  selection_ = baseline && trigger && zl&&isotk && jets && btag2&&!btag3&&mdp;

  /*
  nbins=30; low=0; high=300;
  setLogY(false);
  var="METsig_2012"; xtitle="MET significance (new)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP_METsig2012",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP_METsig2012",0);
  */

  nbins=30; low=0; high=300;
  setLogY(false);
  var="METsig"; xtitle="MET significance";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP_METsig",0);


  nbins=30; low=0; high=400;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP_MET",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP_MET",0);

  //nPV
  nbins=20; low=0; high=40;
  setLogY(false);
  var="nGoodPV"; xtitle="n Good PV";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP_nGoodPV",0);

  //same plots but inverting mdp cut
  selection_ = baseline && trigger && zl&&isotk && jets && btag2&&!btag3&&!mdp;

  /*
  nbins=30; low=0; high=300;
  setLogY(false);
  var="METsig_2012"; xtitle="MET significance (new)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP_METsig2012",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP_METsig2012",0);
  */

  nbins=13; low=0; high=300;
  float varbinsmetsig[14] = {0,10,20,30,40,50,60,70,80,100,120,140,180,300};
  setLogY(false);
  var="METsig"; xtitle="MET significance";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP_METsig",varbinsmetsig);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP_METsig",varbinsmetsig);

  renormalizeBins_=true; renormalizeWidth_=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP_METsig",varbinsmetsig);
  renormalizeBins_=false;

  nbins=30; low=0; high=400;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP_MET",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP_MET",0);

  //nPV
  nbins=20; low=0; high=40;
  setLogY(false);
  var="nGoodPV"; xtitle="n Good PV";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP_nGoodPV",0);

  //now repeat with MET 125 cut
  TCut minmet="MET>125";

  //same plots but adding mdp cut
  selection_ = baseline && trigger && zl&&isotk && jets && btag2&&!btag3&&mdp&&minmet;
  /*
  nbins=30; low=0; high=300;
  setLogY(false);
  var="METsig_2012"; xtitle="MET significance (new)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP-MET125_METsig2012",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP-MET125_METsig2012",0);
  */

  nbins=30; low=0; high=300;
  setLogY(false);
  var="METsig"; xtitle="MET significance";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP-MET125_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP-MET125_METsig",0);

  //nPV
  nbins=20; low=0; high=40;
  setLogY(false);
  var="nGoodPV"; xtitle="n Good PV";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP-MET125_nGoodPV",0);

  //same plots but inverting mdp cut
  selection_ = baseline && trigger && zl&&isotk && jets && btag2&&!btag3&&!mdp&&minmet;

  /*
  nbins=30; low=0; high=300;
  setLogY(false);
  var="METsig_2012"; xtitle="MET significance (new)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP-MET125_METsig2012",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP-MET125_METsig2012",0);
  */

  nbins=13; low=0; high=300;
  setLogY(false);
  var="METsig"; xtitle="MET significance";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP-MET125_METsig",varbinsmetsig);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP-MET125_METsig",varbinsmetsig);

  renormalizeBins_=true;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP-MET125_METsig",varbinsmetsig);
  renormalizeBins_=false;

  //nPV
  nbins=20; low=0; high=40;
  setLogY(false);
  var="nGoodPV"; xtitle="n Good PV";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP-MET125_nGoodPV",0);

  // == now other plots
  selection_ = baseline && trigger && zl && isotk && jets && btag2 && !btag3 && metsigloose;
  nbins=20; low = 0; high=3;
  setLogY(false);
  var="minDeltaPhi20"; xtitle="#Delta #phi_{min}(jet, MET)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bMETsig30_minDeltaPhi20",0);

  //almost an N-1 plot but don't apply DeltaRmax
  selection_ = baseline && trigger && zl && isotk && jets && btag2 && !btag3 && mdp && metsigloose && higgsSR_d;
  nbins=20; low = 0; high=200;
  setLogY(false);
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{jj}> (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bN-2_higgsMass",0);

  //now make it N-1
  selection_ = baseline && trigger && zl && isotk && jets && btag2 && !btag3 && mdp && metsigloose && higgsSR_d &&drmax;
  nbins=20; low = 0; high=200;
  setLogY(false);
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{jj}> (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bN-1_higgsMass",0);

  //N-2 and N-1 of Delta Higgs mass

  //almost an N-1 plot but don't apply DeltaRmax
  selection_ = baseline && trigger && zl && isotk && jets && btag2 && !btag3 && mdp && metsigloose && higgsSR_av;
  nbins=20; low = 0; high=100;
  setLogY(false);
  var="abs(higgsMbb1MassDiff-higgsMbb2MassDiff)"; xtitle="|#Delta m_{jj}| (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bN-2_higgsMassDiff",0);

  //now make it N-1
  selection_ = baseline && trigger && zl && isotk && jets && btag2 && !btag3 && mdp && metsigloose && higgsSR_av &&drmax;
  nbins=20; low = 0; high=100;
  setLogY(false);
  var="abs(higgsMbb1MassDiff-higgsMbb2MassDiff)"; xtitle="|#Delta m_{jj}| (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bN-1_higgsMassDiff",0);

  //N-1 plot of maxDR
  selection_ = baseline && trigger && zl && isotk && jets && btag2 && !btag3 && mdp && metsigloose && higgsSR ;
  nbins=20; low = 0; high=5;
  setLogY(false);
  var="deltaRmax_hh"; xtitle="#Delta R_{max}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bN-1_dRmax",0);


  return;

  selection_=baseline&&trigger&&zl&&isotk&&jets&&btag2&&TCut("CSVbest3<0.244")&&higgsSR&&TCut("METsig_2012>25")&&drmax &&TCut("MET/caloMET<2&&maxTOBTECjetDeltaMult<40");
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_2bT_deltaPhiMin-withMetOverCalo-TOBTEC",0);

  var="abs(min(abs(METphi-caloMETphi),abs(2*3.14159-METphi+caloMETphi))-3.14159)"; xtitle="| |#Delta #phi_{PF,calo}| - #pi |";
  nbins = 30; low=0; high=TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_2bT_RazorNoiseVar-withMetOverCalo-TOBTEC",0);


}

void higgs_dataMC_debug() {

  initHiggsSamples69(false,"ttbar");

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";

  selection_=triggerJetMET;
  nbins=40; low=0; high=40;
  var="nGoodPV"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "debug",0);

  TCut triggerMET = "passMC_PFMET150==1";
  selection_=triggerMET;
  nbins=40; low=0; high=400;
  var="MET"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "debug",0);

  selection_=triggerMET && TCut("nMuons==1");
  drawPlots(var,nbins,low,high,xtitle,"Events", "debug",0);

  selection_=triggerMET && TCut("nMuons==1 && nbjets>=1");
  drawPlots(var,nbins,low,high,xtitle,"Events", "debug",0);

  selection_=triggerMET && TCut("cutPV==1&&nMuons==1 && nbjets>=1");
  drawPlots(var,nbins,low,high,xtitle,"Events", "debug",0);

  selection_=triggerMET && TCut("cutPV==1&&passCleaning==1&&nMuons==1 && nbjets>=1");
  drawPlots(var,nbins,low,high,xtitle,"Events", "debug",0);

  selection_=triggerMET && TCut("cutPV==1&&passCleaning==1&&nMuons==1 && nbjets>=1 && MET/caloMET<2&& maxTOBTECjetDeltaMult<40");
  drawPlots(var,nbins,low,high,xtitle,"Events", "debug",0);

  selection_=triggerMET && TCut("cutPV==1&&passCleaning==1&&nMuons==1 && nbjets>=1 && MET/caloMET<2&& maxTOBTECjetDeltaMult<40 &&jetpt2>50");
  drawPlots(var,nbins,low,high,xtitle,"Events", "debug",0);

  selection_=triggerMET && TCut("cutPV==1&&passCleaning==1&&nMuons+nElectrons==1 && nbjets>=1 && MET/caloMET<2&& maxTOBTECjetDeltaMult<40 &&jetpt2>50");
  drawPlots(var,nbins,low,high,xtitle,"Events", "debug",0);


  selection_=triggerMET && TCut("cutPV==1&&passCleaning==1&&nMuons+nElectrons==1 && nbjets>=1 && MET/caloMET<2&& maxTOBTECjetDeltaMult<40 &&jetpt2>50&&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0");
  drawPlots(var,nbins,low,high,xtitle,"Events", "debug",0);


  selection_=triggerMET && TCut("cutPV==1&&passCleaning==1&&nMuons+nElectrons==1 && CSVbest2>0.898 && MET/caloMET<2&& maxTOBTECjetDeltaMult<40 &&jetpt2>50&&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0");
  drawPlots(var,nbins,low,high,xtitle,"Events", "debug",0);

  selection_=triggerMET && TCut("cutPV==1&&passCleaning==1&&nMuons+nElectrons==1 && nbjets>=1 && MET/caloMET<2&& maxTOBTECjetDeltaMult<40 &&jetpt2>50&&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0 &&MET>180");
  nbins=20; low=0; high=1;
  var="CSVbest2"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "debug",0);

  nbins=3; low=0; high=3;
  var="nPUjets30"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "debug",0);

  nbins=8; low=0; high=8;
  var="njets30"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "debug",0);

}

void higgs_dataMC_control_SL_JERtest() {
  initHiggsSamples69(false,"metsigtest1");

  setOutputDirectory("plots_JERtest");

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";
  //cannot use this with MET PD only
  //  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets20==4"; //switch back to 20 GeV jets
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut mdp = "minDeltaPhi30>0.3"; // does this need to go to 20?
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

  TCut metsigveryloose="METsig_2012>25";

  TCut drmax = "deltaRmax_hh<2.2";
  //  TCut drmin = "deltaRmin_hh<1.9";

  //  setSampleScaleFactor("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim",0.8);

  setStackMode(true,false,false); //stack,norm,label override
  setLogY(false);

  for (int iii=0; iii<2; iii++) {
  //2b pre-selection ; basically the minimum required by the trigger, plus SL

  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btag2;

  if (iii==1) setSampleWeightFactor("TTJets_ResolutionTest_PU_S10_CU0001_v69","topPtWeight");
  else resetSampleWeightFactors();

  TString suffix=getSampleWeightFactor("TTJets_ResolutionTest_PU_S10_CU0001_v69");
  if (suffix!="") suffix.Prepend("-");

  currentConfig_=configDescriptions_.getDefault();
  nbins=30; low=0; high=300;
  var="METsig_2012"; xtitle="MET significance (new)";
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig2012_JER0"+suffix,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig2012_JER0"+suffix,0);

  currentConfig_=configDescriptions_.getCorrected();
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig2012_JERbias"+suffix,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig2012_JERbias"+suffix,0);

  currentConfig_=configDescriptions_.getDefault();
  var="METsig"; xtitle="MET significance (old)";
  nbins=30; low=0; high=300;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig_JER0"+suffix,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig_JER0"+suffix,0);

  currentConfig_=configDescriptions_.getCorrected();
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig_JERbias"+suffix,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig_JERbias"+suffix,0);

  //MET, for completeness
  currentConfig_=configDescriptions_.getDefault();
  nbins=30; low=0; high=400;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_MET_JER0"+suffix,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_MET_JER0"+suffix,0);

  currentConfig_=configDescriptions_.getCorrected();
  nbins=30; low=0; high=400;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_MET_JERbias"+suffix,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_MET_JERbias"+suffix,0);

  //now add a MET cut; remove other extraneous cuts that are fairly specific to our analysis
  //keep trigger, SL, jet2 (for the trigger), CSVT
  selection_ = baseline && trigger &&TCut("MET>125") && sl &&jet2 && btag2;

  //first, no correction
  currentConfig_=configDescriptions_.getDefault();
  nbins=60; low=0; high=300;
  setLogY(false);
  var="jetpt1"; xtitle="lead jet pT (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_loose_jetpt1_JER0"+suffix,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_loose_jetpt1_JER0"+suffix,0);

  nbins=60; low=0; high=300;
  setLogY(false);
  var="jetpt2"; xtitle="2nd jet pT (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_loose_jetpt2_JER0"+suffix,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_loose_jetpt2_JER0"+suffix,0);

  currentConfig_=configDescriptions_.getCorrected();
  nbins=60; low=0; high=300;
  setLogY(false);
  var="jetpt1"; xtitle="lead jet pT (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_loose_jetpt1_JERbias"+suffix,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_loose_jetpt1_JERbias"+suffix,0);

  nbins=60; low=0; high=300;
  setLogY(false);
  var="jetpt2"; xtitle="2nd jet pT (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_loose_jetpt2_JERbias"+suffix,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_loose_jetpt2_JERbias"+suffix,0);

  //dilepton selection
  //i'm assuming that non-ttbar will be very small due to MET and b-tag requirements
  TCut dl = "nMuons+nElectrons==2";
  selection_ = baseline && trigger &&TCut("MET>125") && dl &&jet2 && btag2;

  currentConfig_=configDescriptions_.getDefault();
  nbins=30; low=100; high=400;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_DL_loose_MET_JER0"+suffix,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_DL_loose_MET_JER0"+suffix,0);

  currentConfig_=configDescriptions_.getCorrected();

  nbins=30; low=100; high=400;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_DL_loose_MET_JERbias"+suffix,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_DL_loose_MET_JERbias"+suffix,0);



  //first, no correction
  currentConfig_=configDescriptions_.getDefault();
  nbins=40; low=0; high=300;
  setLogY(false);
  var="jetpt1"; xtitle="lead jet pT (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_DL_loose_jetpt1_JER0"+suffix,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_DL_loose_jetpt1_JER0"+suffix,0);

  nbins=40; low=0; high=300;
  setLogY(false);
  var="jetpt2"; xtitle="2nd jet pT (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_DL_loose_jetpt2_JER0"+suffix,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_DL_loose_jetpt2_JER0"+suffix,0);

  currentConfig_=configDescriptions_.getCorrected();
  nbins=40; low=0; high=300;
  setLogY(false);
  var="jetpt1"; xtitle="lead jet pT (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_DL_loose_jetpt1_JERbias"+suffix,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_DL_loose_jetpt1_JERbias"+suffix,0);

  nbins=40; low=0; high=300;
  setLogY(false);
  var="jetpt2"; xtitle="2nd jet pT (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_DL_loose_jetpt2_JERbias"+suffix,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_DL_loose_jetpt2_JERbias"+suffix,0);

  currentConfig_=configDescriptions_.getDefault();
  }

  // =============== now, MC-MC comparisons ===============
  clearSamples();
  dodata_=false;

  currentConfig_=configDescriptions_.getDefault();

  addSample("TTJets_ResolutionTest_PU_S10_CU0001_v69",kBlue,"Nominal");
  addSample("TTJets_ResolutionTest_PU_S10_CU0001_JERbias_v69",kRed,"JER smearing");

  //this is quite a hack.
  //i don't have a pre-coded way to plot JERbias and non-JERbias at the same time, so as a hack, just make a symlink from the regular JERbias name to 
  //a name with JER0 in the first stub and JERbias in the second. kind of ugly, but ok for now
  //manually set scale factor for JERbias sample
  setSampleScaleFactor("TTJets_ResolutionTest_PU_S10_CU0001_JERbias_v69",234.0 / getTree("TTJets_ResolutionTest_PU_S10_CU0001_JERbias_v69")->GetEntries());
  resetSampleWeightFactors(); //critical!

  //use only gen-level
  TCut ttbar1l = "ttbarDecayCode>=3 && ttbarDecayCode<=7";
  TCut ttbar2l = "ttbarDecayCode==1 || ttbarDecayCode>=8";
  TCut ttbar0l = "ttbarDecayCode==2";
  selection_ = TCut("cutPV==1")&&ttbar1l; //that's it!

  doRatio_=true; ratioMin = 0.85; ratioMax = 1.15;

  setStackMode(false,false,false);
  nbins=60; low=0; high=300;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "MConly_JER0_JERbias_ttbar1l_MET");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "MConly_JER0_JERbias_ttbar1l_MET");

  nbins=60; low=0; high=400;
  setLogY(false);
  var="jetpt1"; xtitle="lead jet pT (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "MConly_JER0_JERbias_ttbar1l_jetpt1");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "MConly_JER0_JERbias_ttbar1l_jetpt1");

  nbins=60; low=0; high=400;
  setLogY(false);
  var="jetpt2"; xtitle="2nd jet pT (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "MConly_JER0_JERbias_ttbar1l_jetpt2");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "MConly_JER0_JERbias_ttbar1l_jetpt2");


  nbins=40; low=0; high=300;
  setLogY(false);
  var="jetpt3"; xtitle="3rd jet pT (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "MConly_JER0_JERbias_ttbar1l_jetpt3");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "MConly_JER0_JERbias_ttbar1l_jetpt3");

  nbins=40; low=0; high=200;
  setLogY(false);
  var="jetpt4"; xtitle="4th jet pT (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "MConly_JER0_JERbias_ttbar1l_jetpt4");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "MConly_JER0_JERbias_ttbar1l_jetpt4");


  // == now try a zero-MET sample ==
  doRatio_=true; ratioMin = 0.7; ratioMax = 1.3;

  selection_ = TCut("cutPV==1")&&ttbar0l; //that's it!
  nbins=50; low=0; high=100;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "MConly_JER0_JERbias_ttbar0l_MET");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "MConly_JER0_JERbias_ttbar0l_MET");


}

void higgs_dataMC_control_SL_noskim() {

  //goal do data / MC comparisons for plots that can't be made with the skim
  initHiggsSamples69(false,"ttbar znunu qcd"); //use all samples except signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  setOutputDirectory("plots_Control_SL");

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";
  //cannot use this with MET PD only
  //  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets20==4"; //switch back to 20 GeV jets
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut mdp = "minDeltaPhi30>0.3"; // does this need to go to 20?
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

  TCut metsigveryloose="METsig>30"; //go back to old METsig

  TCut drmax = "deltaRmax_hh<2.2";
  //  TCut drmin = "deltaRmin_hh<1.9";

  //  setSampleScaleFactor("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim",0.8);

  setStackMode(true,false,false); //stack,norm,label override
  setLogY(false);

  //2b pre-selection: 
  selection_ = baseline && trigger && sl && jet2&&mdp && btag2;

  nbins=10; low=0; high=10;
  setLogY(false);
  var="njets20"; xtitle="Jet multiplicity (p_{T} > 20 GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_njets20",0);

  nbins=4; low=0; high=4;
  setLogY(true);
  var="nPUjets20"; xtitle="Number of pileup jets (p_{T} > 20 GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_nPUjets20",0);
 
 //2b pre-selection + METsig
  selection_ = baseline && trigger && sl && jet2&&mdp && btag2 &&metsigveryloose;

  nbins=10; low=0; high=10;
  setLogY(false);
  var="njets20"; xtitle="Jet multiplicity (p_{T} > 20 GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselectionMETsig30_njets20",0);

  nbins=4; low=0; high=4;
  setLogY(true);
  var="nPUjets20"; xtitle="Number of pileup jets (p_{T} > 20 GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselectionMETsig30_nPUjets20",0);


}

void higgs_dataMC_control_SL() {

  //goal do data / MC comparisons; use skim to save time. means I can't make a useful njets plot
  initHiggsSamples69(true,"ttbar znunu qcd"); //use all samples except signal
  //this is the current default
  //  initHiggsSamples69(true,"qcd ttbar"); //don't bother with negligible znunu and qcd

  //special test
  //  initHiggsSamples69(false,"metsigtest1");

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  setOutputDirectory("plots_Control_SL");


  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1";
  //cannot use this with MET PD only
  //  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets20==4"; //switch back to 20 GeV jets
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut mdp = "minDeltaPhi30>0.3"; // does this need to go to 20?
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

  TCut metsigveryloose="METsig>30"; //go back to old METsig

  TCut drmax = "deltaRmax_hh<2.2";
  //  TCut drmin = "deltaRmin_hh<1.9";

  //  setSampleScaleFactor("TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-skim",0.8);

  setStackMode(true,false,false); //stack,norm,label override
  setLogY(false);

  //2b pre-selection ; basically the minimum required by the trigger, plus SL
  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btag2;

  //Owen's METsig bins
  //float binning[]={0,10,20,30,50,100,150,200};
  /* no need to mess with new METsig for the moment
  nbins=30; low=0; high=300;
  setLogY(false);
  var="METsig_2012"; xtitle="MET significance (new)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig2012",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig2012",0);
  */

  //this was a special test
//   double fullLumi = lumiScale_;

//   // to check for run dependence
//   setDatasetToDraw("META"); lumiScale_ = 807.1;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig2012_DEBUG_2012A",0);

//   setDatasetToDraw("METB"); lumiScale_ = 4421;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig2012_DEBUG_2012B",0);

//   setDatasetToDraw("METC"); lumiScale_ = 495+6402;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig2012_DEBUG_2012C",0);

//   setDatasetToDraw("METD"); lumiScale_ = 5956+1317;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig2012_DEBUG_2012D",0);

//   setDatasetToDraw("MET"); 
//   lumiScale_=fullLumi;

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

  //nPV
  nbins=20; low=0; high=40;
  setLogY(false);
  var="nGoodPV"; xtitle="n Good PV";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_nGoodPV",0);


  // === now in bins of PV
  /* for now, don't make these plots. may want to add them back in eventually
  drawFilenameOnPlot_=true;

  TCut pvcut[4];
  pvcut[0] = "nGoodPV <=10";
  pvcut[1] = "nGoodPV >10 && nGoodPV<=15";
  pvcut[2] = "nGoodPV >15 && nGoodPV<=20";
  pvcut[3] = "nGoodPV >20";

  TH1D* ratiosByPV[4];

  for (int ipv=0; ipv<4;ipv++) {
    selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btag2 && pvcut[ipv];
    TString pvdesc=jmt::fortranize(pvcut[ipv].GetTitle());
    nbins=30; low=0; high=300;
    setLogY(false);
    var="METsig_2012"; xtitle="MET significance (new)";
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig2012_"+pvdesc,0);
    ratiosByPV[ipv] = (TH1D*)ratio->Clone("ratiosByPV_"+pvdesc);
    ratiosByPV[ipv]->SetLineColor(ipv+1);
    ratiosByPV[ipv]->SetMarkerColor(ipv+1);
    setLogY(true);
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig2012_"+pvdesc,0);

    nbins=30; low=0; high=300;
    setLogY(false);
    var="METsig"; xtitle="MET significance (old)";
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig_"+pvdesc,0);
    setLogY(true);
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig_"+pvdesc,0);
    nbins=30; low=0; high=400;
    setLogY(false);
    var="MET"; xtitle="MET (GeV)";
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_MET_"+pvdesc,0);
    setLogY(true);
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_MET_"+pvdesc,0);
  }
  */
/*
root [101] ratiosByPV[0]->SetLineColor(kMagenta+4)
root [103] ratiosByPV[1]->SetLineColor(kMagenta+2)
root [104] ratiosByPV[2]->SetLineColor(kMagenta)
root [105] ratiosByPV[3]->SetLineColor(kMagenta-7)
root [106] ratiosByPV[3]->SetMarkerColor(kMagenta-7)
root [107] ratiosByPV[2]->SetMarkerColor(kMagenta)
root [108] ratiosByPV[1]->SetMarkerColor(kMagenta+2)
root [109] ratiosByPV[0]->SetMarkerColor(kMagenta+4)

*/ 
/*
 renewCanvas();
  ratiosByPV[0]->Draw();
  ratiosByPV[1]->Draw("same");
  ratiosByPV[2]->Draw("same");
  ratiosByPV[3]->Draw("same");
  thecanvas->SaveAs("higgs_dataMC_SL_preselection_METsig2012_PVratios.eps");

  return;
  drawFilenameOnPlot_=false;
  /// === done with bins of PV
  */

  //minDeltaPhi distribution (no MET or METsig cut)
  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2 && btag2;
  nbins=30; low=0; high=3;
  setLogY(false);
  var="minDeltaPhi30"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselectionNoMdpNoMETsig_minDeltaPhi30",0);

  //add the standard METsig cut
  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2 && btag2 &&metsigveryloose;
  nbins=30; low=0; high=3;
  setLogY(false);
  var="minDeltaPhi30"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselectionNoMdp_minDeltaPhi30",0);

  //CSV3
  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btag2 &&metsigveryloose;
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
  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btag2 &&btag3 &&metsigveryloose;

  //CSV4
  nbins=10; low=0; high=1;
  setLogY(false);
  var="CSVbest4"; xtitle="4th CSV value";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_3bPreselection_CSVbest4",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_3bPreselection_CSVbest4",0);

  //now apply full b-tag selection
  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btag2 &&btag3&&btag4;

  /*
  //MET sig
  nbins=15; low=0; high=300;
  setLogY(false);
  var="METsig_2012"; xtitle="MET significance (new)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bPreselection_METsig2012",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bPreselection_METsig2012",0);
  */
  //old MET sig
  nbins=15; low=0; high=300;
  setLogY(false);
  var="METsig"; xtitle="MET significance (old)";
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
  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btag2 &&btag3&&btag4 &&metsigveryloose;
  nbins=10; low=0; high=100;
  setLogY(false);
  var="abs(higgsMbb1MassDiff-higgsMbb2MassDiff)"; xtitle="|#Delta m_{jj}|";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bMETsig30_DeltaMbb",0);

  //higgs mass
  nbins=20; low=0; high=200;
  setLogY(false);
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{jj}>";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bMETsig30_AvMbb",0);

  //now add the higgs SR
  selection_ = baseline && trigger && sl && (njets4||njets5)&&jet2&&mdp && btag2 &&btag3&&btag4 &&metsigveryloose &&higgsSR;
  nbins=20; low=0; high=5;
  setLogY(false);
  var="deltaRmax_hh"; xtitle="#Delta R_{max}(b,b)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_hhSR-METsig30_DRmax",0);

  //true N - 1 plots -- for the AN
  selection_ = baseline && trigger && sl && (njets4||njets5)&&jet2&&mdp && btag2 &&btag3&&btag4 &&metsigveryloose &&higgsSR_d &&drmax;
  //higgs mass
  nbins=20; low=0; high=200;
  setLogY(false);
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{jj}>";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_Nm1METsig30_AvMbb",0);

  selection_ = baseline && trigger && sl && (njets4||njets5)&&jet2&&mdp && btag2 &&btag3&&btag4 &&metsigveryloose &&higgsSR_av &&drmax;
  nbins=10; low=0; high=100;
  setLogY(false);
  var="abs(higgsMbb1MassDiff-higgsMbb2MassDiff)"; xtitle="|#Delta m_{jj}|";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_Nm1METsig30_DeltaMbb",0);

  //
  selection_ = baseline && trigger && sl && (njets4||njets5)&&jet2&&mdp && btag2 &&btag3&&btag4  &&higgsSR &&drmax;
  nbins=30; low=0; high=300;
  setLogY(true);
  var="METsig"; xtitle="MET significance";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_Nm1_METsig",0);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_Nm1_METsig",0);



}

void higgs_SLratios2() {
  //owen asked me to calculate  [ N(4bSL,data) / N(4bSL,MC) ]   /   [ N(2bSL,data) / N(2bSL,MC)]
  initHiggsSamples69(true,"znunu ttbar");

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

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

void higgs_whyLeptonLost() {

  initHiggsSamples69();
  setOutputDirectory("plots_higgs_whyLeptonLost");

  int nbins;
  float low,high;
  TString var,xtitle;

  savePlots_=true;
  setQuiet(false);

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  //20 June -- update
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0 && MET/caloMET<2 && maxTOBTECjetDeltaMult<40";
  //20 june -- remove the all jet trigger: ||passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_PFMET150==1";
  TCut zl = "nMuons==0&&nElectrons==0";
  TCut tauveto="nTausLoose==0";
  TCut isotk1="nIsoTracks15_005_03==0";

  TCut njets="njets20==4 || njets20==5";
  TCut jet2="jetpt2>50";
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut mdp = "minDeltaPhi20>0.3"; //4 july

  TCut higgsmass = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsdiff = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";

  TCut drmax = "deltaRmax_hh<2.2";
  TCut drmin = "deltaRmin_hh<1.9";

  TCut metsig = "METsig>30";

  setStackMode(true,false,false); //stack, norm, labels

  stackSignal_=false; //let's unstack signal

  clearSamples();
  addSample("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1799_v69-skim",kBlue-10,"t#bar{t} (2l)");
  addSample("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1800_v69-skim:ttbarDecayCode==7",kMagenta-7,"t#bar{t} (#tau #rightarrow h)");
  addSample("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1800_v69-skim:ttbarDecayCode==3||ttbarDecayCode==4",kRed-7,"t#bar{t} (W #rightarrow e/#mu)");
  addSample("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1800_v69-skim:ttbarDecayCode==5||ttbarDecayCode==6",kOrange+10,"t#bar{t} (#tau #rightarrow e/#mu)");

  //higgs mass
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2&&mdp && btag2 && btag3 && btag4  && higgsdiff && drmax && metsig;
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{jj}>";
  nbins=20; low= 0; high=200;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_higgsmass_splitTtbar",0,"GeV");

  //all cuts except metsig
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2&&mdp && btag2 && btag3 && btag4  && higgsdiff&&higgsmass && drmax ;
  var="METsig"; xtitle="MET significance";
  nbins=8; low= 0; high=160;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_METsig_splitTtbar",0);

  //all cuts
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2&&mdp && btag2 && btag3 && btag4  && higgsdiff&&higgsmass && drmax &&metsig;
  var="higgs1tauMatch1+higgs1tauMatch2+higgs2tauMatch1+higgs2tauMatch2"; xtitle="Higgs jets matched to gen tau";
  nbins=3; low= 0; high=3;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_HtauGenMatch_splitTtbar",0);


  //full signal selection except CSV4
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2&&mdp && btag2 && btag3 && higgsdiff&&higgsmass && drmax && metsig;

  clearSamples();
  addSample("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1799_v69-skim",kBlue-10,"t#bar{t} (2l)");
  addSample("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1800_v69-skim:hadWcode==11",kBlue+2,"t#bar{t} (W #rightarrow ud)");
  addSample("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1800_v69-skim:hadWcode==1100",kGreen+2,"t#bar{t} (W #rightarrow cs)");

  var="CSVbest4"; xtitle="4th best CSV value";
  nbins=20; low= 0; high=1;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_splitTtbarByHadW_CSV4",0);

  var="nPUjets20"; xtitle="Number of PU jets (pT>20 GeV)";
  nbins=3; low= 0; high=3;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_splitTtbarByHadW_nPU20",0);


}

void higgs_Nminus1() {
  initHiggsSamples69(true,"qcd ttbar znunu hh250 hh400");
  setOutputDirectory("plots_higgs_Nminus1");

  int nbins;
  float low,high;
  TString var,xtitle;

  savePlots_=true;
  setQuiet(false);

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  //20 June -- update
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0 && MET/caloMET<2 && maxTOBTECjetDeltaMult<40";
  //20 june -- remove the all jet trigger: ||passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_PFMET150==1";
  TCut zl = "nMuons==0&&nElectrons==0";
  TCut tauveto="nTausLoose==0";
  TCut isotk1="nIsoTracks15_005_03==0";

  TCut njets="njets20==4 || njets20==5"; //back to 20 GeV
  TCut jet2="jetpt2>50";
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut mdp = "minDeltaPhi20>0.3"; //decided to use this everywhere

  TCut higgsmass = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsdiff = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";

  TCut drmax = "deltaRmax_hh<2.2";

  TCut metsig = "METsig>30";

  setStackMode(true,false,false); //stack, norm, labels

  stackSignal_=false; //let's unstack signal

  //full selection 
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && btag3 && btag4 &&mdp && higgsmass && higgsdiff && drmax && metsig;
  var="bjetpt1"; xtitle="lead b jet pT";
  nbins=20; low= 0; high=400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_bjetpt1",0,"GeV");

  //still full selection
  var="jetpt1"; xtitle="lead jet pT";
  nbins=30; low= 0; high=300; //need 10 GeV Bins
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_jetpt1",0,"GeV");

  var="jetpt2"; xtitle="2nd jet pT";
  nbins=30; low= 0; high=200; //need 10 GeV Bins
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_jetpt2",0,"GeV");

  var="jetpt3"; xtitle="3rd jet pT";
  nbins=20; low= 0; high=200; //need 10 GeV Bins
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_jetpt3",0,"GeV");

  var="jetpt4"; xtitle="4th jet pT";
  nbins=20; low= 0; high=200; //need 10 GeV Bins
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_jetpt4",0,"GeV");

  //not too interesting
  var="higgsWCandMass"; xtitle="W Candidate mass (GeV)";
  nbins=40; low= 0; high=200;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_WcandMass",0,"GeV");


  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && btag3 && btag4&&mdp && higgsmass && higgsdiff && drmax && metsig;
  var="njets20_5p0-njets20"; xtitle="number of forward jets (pT>20)";
  nbins=4; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_njets20forward",0);

  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && btag3 && btag4&&mdp && higgsmass && higgsdiff && drmax && metsig;
  var="deltaPhi_hh"; xtitle="#Delta #phi (h,h)";
  nbins=10; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_hhDeltaPhi",0);

  //full selection except mdp
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && btag3 && btag4 && higgsmass && higgsdiff && drmax && metsig;
  var="deltaPhi1"; xtitle="#Delta #phi (jet1,MET)";
  nbins=10; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_deltaPhi1",0);
  var="deltaPhi2"; xtitle="#Delta #phi (jet2,MET)";
  nbins=10; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_deltaPhi2",0);
  var="deltaPhi3"; xtitle="#Delta #phi (jet3,MET)";
  nbins=10; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_deltaPhi3",0);
  var="minDeltaPhi20"; xtitle="min #Delta #phi (jet1..3,MET)";
  nbins=10; low= 0; high=3;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_allcuts_deltaPhiMin",0);

  //3b sample; veto on 4th b
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && btag3 && (!btag4) && higgsmass && higgsdiff && drmax && metsig;
  var="deltaPhi2"; xtitle="#Delta #phi (jet2,MET)";
  nbins=10; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_3b_deltaPhi2",0);

  //2b sample; veto on 3rd b
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && (!btag3)&& higgsmass && higgsdiff && drmax && metsig;
  var="deltaPhi2"; xtitle="#Delta #phi (jet2,MET)";
  nbins=20; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_2b_deltaPhi2",0);

  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && (!btag3)&& higgsmass && higgsdiff && drmax && metsig;
  var="minDeltaPhi20"; xtitle="min #Delta #phi (jet1..3,MET)";
  nbins=40; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_2b_deltaPhiMin",0);


  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && TCut("CSVbest3<0.244")&& higgsmass && higgsdiff && drmax && metsig;
  var="minDeltaPhi20"; xtitle="min #Delta #phi (jet1..3,MET)";
  nbins=40; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_2bT_deltaPhiMin",0);


  //n leptons
  selection_=baseline&&triggers &&  tauveto && isotk1 && njets&&jet2 && btag2 && btag3 && btag4&&mdp && higgsmass && higgsdiff && drmax && metsig;
  var="nMuons+nElectrons"; xtitle="e + #mu";
  nbins=4; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_ePlusMu",0);

  //POG taus
  selection_=baseline&&triggers && zl&&  isotk1 && njets&&jet2 && btag2 && btag3 && btag4&&mdp && higgsmass && higgsdiff && drmax && metsig;
  var="nTausLoose"; xtitle="Loose taus";
  nbins=3; low= 0; high=3;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_looseTaus",0);

  //iso tracks
  selection_=baseline&&triggers && zl&& tauveto && njets&&jet2 && btag2 && btag3 && btag4&&mdp && higgsmass && higgsdiff && drmax && metsig;
  var="nIsoTracks15_005_03"; xtitle="15 GeV iso tracks";
  nbins=4; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_isoTk15",0);

  //iso tracks
  selection_=baseline&&triggers && zl&& tauveto  && isotk1 && njets&&jet2 && btag2 && btag3 && btag4&&mdp && higgsmass && higgsdiff && drmax && metsig;
  var="nIsoTracks5_005_03"; xtitle="5 GeV iso tracks";
  nbins=4; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_isoTk5",0);

  //jet multiplicity
  selection_=baseline&&triggers && zl&& tauveto && isotk1&&jet2 && btag2 && btag3 && btag4&&mdp && higgsmass && higgsdiff && drmax  && metsig;
  var="njets20"; xtitle="jet multiplicity (20 GeV)";
  nbins=9; low= 0; high=9;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_njets20",0);

  //3rd b tag
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2&&mdp  && higgsmass && higgsdiff && drmax && metsig;
  var="CSVbest3"; xtitle="3rd CSV value";
  nbins=20; low= 0; high=1;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_CSVbest3",0);

  //4th btag
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && btag3&&mdp && higgsmass && higgsdiff && drmax && metsig;
  var="CSVbest4"; xtitle="4th CSV value";
  nbins=20; low= 0; high=1;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_CSVbest4",0);

  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && btag3&&mdp && higgsmass && higgsdiff && drmax && TCut("METsig>100");
  var="CSVbest4"; xtitle="4th CSV value";
  nbins=20; low= 0; high=1;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_CSVbest4-METsig100",0);


  //higgs mass
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && btag3 && btag4 &&mdp && higgsdiff && drmax && metsig;
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{jj}>";
  nbins=20; low= 0; high=200;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_higgsmass",0,"GeV");

  //higgs mass difference
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && btag3 && btag4&&mdp && higgsmass  && drmax && metsig;
  var="abs(higgsMbb1MassDiff-higgsMbb2MassDiff)"; xtitle="#Delta higgs mass";
  nbins=20; low= 0; high=80;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_higgsmassdiff",0);

  //drmax
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && btag3 && btag4&&mdp && higgsmass && higgsdiff && metsig;
  var="deltaRmax_hh"; xtitle="max #Delta R";
  nbins=30; low= 0; high=6;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_maxDR",0);

  //2d plot
  //draw2d(var,15,low,high,"deltaRmin_hh",15,0,4,xtitle,"min #Delta R","higgs_Nm1_maxDRminDR_hh200",0,0,"SMS-HbbHbb_mHiggsino-200_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1807_v69");
  //draw2d(var,15,low,high,"deltaRmin_hh",15,0,4,xtitle,"min #Delta R","higgs_Nm1_maxDRminDR_hh350",0,0,"SMS-HbbHbb_mHiggsino-350_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1811_v69");

  //linear combination
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && btag3 && btag4&&mdp && higgsmass && higgsdiff && metsig;
  var="3.318-0.139*deltaRmax_hh-0.622*deltaRmin_hh"; xtitle="#DeltaR Fisher";
  nbins=30; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_DRF",0);


  //drmin
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && btag3 && btag4&&mdp && higgsmass && higgsdiff  && drmax && metsig;
  var="deltaRmin_hh"; xtitle="min #Delta R";
  nbins=40; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_minDR",0);

  //full selection except METsig
  setPlotMinimum(0.1);
  selection_=baseline&&triggers && zl&& tauveto && isotk1&& njets&&jet2 && btag2 && btag3 && btag4 &&mdp&& higgsmass && higgsdiff && drmax;
  var="METsig"; xtitle="MET significance";
  nbins=20; low= 0; high=200;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_METsig",0);
  resetPlotMinimum();

  //full selection except METsig
  setPlotMinimum(0.1);
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && btag3 && btag4&&mdp && higgsmass && higgsdiff && drmax;
  var="MET"; xtitle="MET";
  nbins=20; low= 0; high=300;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_MET",0);
  resetPlotMinimum();

}


void higgs_printCutflowTable(const bool latexMode=false,TString region="4bSR", const  bool unweighted=false, const bool makePlot=false) {

  if (unweighted)  initHiggsSamples69(true,"ttbar hh200 hh350"); //reduced list of samples
  else initHiggsSamples69(true);

  if (!(region.Contains("4bSR")||region.Contains("2bSR")||region.Contains("4bSB")||region.Contains("2bSB"))) {cout<<"check region!"<<endl;return;}

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
  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";

  TCut njets="njets20>=4&&njets20<=5";
  //TCut njets="njets20==4";
  //  TCut secondjet="jetpt2>50";
  TCut secondjet="jetpt2>50 && bjetpt1>100";

  TCut btagpre="CSVbest2>0.898";
  TCut btag3 = "CSVbest3>0.679";
  TCut btag4 = "CSVbest4>0.244"; //TCut btag4 = "CSVbest4>0.679";

  TCut higgsmass = "((0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)) && abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;

  TCut drmax = "deltaRmax_hh<2.2"; //owen's been using 2.2

  TCut metsigloose = "METsig>30";

  TCut skim = baseline&&triggers&&njets&&btagpre;

  std::vector<pair<TString, TCut> > cuts;
  if (unweighted) {
    cuts.push_back(make_pair("PV",TCut("cutPV==1")));  

    cuts.push_back(make_pair("MET cleaning",TCut("passCleaning==1 && buggyEvent==0")));  

    cuts.push_back(make_pair( "Trigger",triggers));
    cuts.push_back(make_pair( "njets 4-5",njets));
    
    cuts.push_back(make_pair( "2 CSVT",btagpre));
  }
  else     cuts.push_back(make_pair( "Preselection",skim));


  cuts.push_back(make_pair("jet2 >50",secondjet));
  cuts.push_back(make_pair("minDeltaPhi20",TCut("minDeltaPhi20>0.3")));
  //
  cuts.push_back(make_pair( "Lepton/isotk vetoes",zl&&isotk)); 

  if (region.Contains("4b")) {
    cuts.push_back(make_pair( "3rd btag",btag3));
    cuts.push_back(make_pair( "4th btag",btag4));
  }
  else {
    cuts.push_back(make_pair( "not 3rd btag",!btag3));
  }

  if (region.Contains("SR"))  cuts.push_back(make_pair( "Higgs masses",higgsmass));
  else  cuts.push_back(make_pair( "Higgs mass SB",higgsSB));

  cuts.push_back(make_pair( "#Delta R_{max}",drmax));

  //  cuts.push_back(make_pair("#Delta R_{min}",drmin));

  int max= region.Contains("MET")? 300 : 150; //MET : METsig
  int step = region.Contains("MET")? 10 : 5;
  for (int icut=10;icut<=max;icut+=step) {
    TString label;
    //    if (region.Contains("oldsig")) label = "METsig>";
     if (region.Contains("MET"))  label = "MET>";
    else                              label = "METsig>";
    label+=icut;
    cuts.push_back(make_pair( label,label.Data()));
  }

  //first print a header row

  const TString dividerHead = latexMode ? " & ":" *|* ";
  const TString divider = latexMode ? " & ":" | ";
  const TString lineEnd=latexMode ? "\\\\":" |";
  if (!latexMode) cout<<"|* ";
  cout<<"Cut \t\t"<<dividerHead;
  bool foundsignals=false;
  for (unsigned int isample = 0; isample<samples_.size(); isample++) {
    if ( !isSampleSM(samples_.at(isample)) && !foundsignals) {
      cout<<" Total SM"<<dividerHead;
      foundsignals=true;
    }
    cout<<sampleLabel_[samples_.at(isample)];
    if (isample==samples_.size()-1) cout<<lineEnd;
    else   cout<<dividerHead;
  }
  cout<<endl;

  TCut totalcut = cuts.at(0).second;
  for (unsigned int icut=0; icut<cuts.size(); icut++) {

    if (icut != 0) totalcut = totalcut&&cuts.at(icut).second;
    selection_=totalcut;

    drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
    //now we counted events for all of the samples for this set of cuts.
    if (!latexMode)    cout<<"| ";
    cout<<cuts.at(icut).first<<" \t\t"<<divider;
    foundsignals=false;
    double totalsmvalue=0;
    for (unsigned int isample = 0; isample<samples_.size(); isample++) {

      TString thissamplename = samples_.at(isample);
      if ( !isSampleSM(thissamplename) && !foundsignals) {
	totalsmvalue = getIntegral("totalsm");
	if (unweighted)	cout<<setprecision(10)<<totalsmvalue<<divider;
	else	cout<<jmt::format_nevents(totalsmvalue,getIntegralErr("totalsm"),false,latexMode)<<divider;
	foundsignals=true;
      }
      double n = getIntegral(thissamplename);
      double e = getIntegralErr(thissamplename); 
      if (unweighted)     cout<<setprecision(10)<<n;
      else     cout<<jmt::format_nevents(n,e,false,latexMode);

      if (foundsignals) cout<<" ("<<setprecision(2)<<n/sqrt(totalsmvalue)<<") ";
      

      if (isample==samples_.size()-1) cout<<lineEnd<<endl;
      else   cout<<divider;
    }
  }

}

void higgmb_cutflowE_bbangles() {
  // check keith's maxDeltaPhi_bb_bb against the deltaRmax_hh variable
  //also check MET versus METsig

  // 27 June -- updating to comment out the maxDeltaPhi study and repeat the METsig study with MET versus old METsig versus new METsig
  //3 july -- update again

  initHiggsSamples69();

    TString name200="TChihh-HbbHbb_mHiggsino-200_mLSP-1_8TeV-Pythia6Z_jgsmith-SMS-HbbHbb_mHiggsino-200_mLSP-1_1_UCSB1807_v69";
    //TString name250="TChihh_250_2";
  TString name350="TChihh-HbbHbb_mHiggsino-350_mLSP-1_8TeV-Pythia6Z_jgsmith-SMS-HbbHbb_mHiggsino-350_mLSP-1_1_UCSB1811_v69";
  //  TString name450="TChihh_450_2";

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0"; 
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_PFMET150==1";
  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";

  TCut njets="(njets30==4 || njets30==5)&&jetpt2>50";
  TCut btag = "CSVbest2>0.898 && CSVbest3>0.679 && CSVbest4>0.244";

  TCut higgsmass = "((0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)) && abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";


/* comment out the DeltaPhi study
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
*/

  //change gears
  //MET versus METsig
  TCut drmax = "deltaRmax_hh<2.2";
  selection_=baseline && triggers &&zl &&isotk&&njets&&btag&&higgsmass&&drmax;

  var="METsig"; xtitle="Old MET significance";
  nbins=20; low= 0; high=200;
  setStackMode(false,false,false); //no stack
  setLogY(true);
  setPlotMinimum(0.1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_optcuts_metsig",0);
  setLogY(false);
  drawEffVRej(name200,"totalsm","hh 200","total SM",true);
  TGraphAsymmErrors* metsig200  =(TGraphAsymmErrors*)  effgraph->Clone("metsig200");
  drawEffVRej(name350,"totalsm","hh 350","total SM",true);
  TGraphAsymmErrors* metsig350  =(TGraphAsymmErrors*)  effgraph->Clone("metsig350");
//   drawEffVRej(name450,"totalsm","hh 450","total SM",true);
//   TGraphAsymmErrors* metsig450  =(TGraphAsymmErrors*)  effgraph->Clone("metsig450");

  var="METsig_2012"; xtitle="New MET significance";
  nbins=20; low= 0; high=200;
  setStackMode(false,false,false); //no stack
  setLogY(true);
  setPlotMinimum(0.1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_optcuts_metsignew",0);
  setLogY(false);
  drawEffVRej(name200,"totalsm","hh 200","total SM",true);
  TGraphAsymmErrors* metsig200new  =(TGraphAsymmErrors*)  effgraph->Clone("metsig200new");
  drawEffVRej(name350,"totalsm","hh 350","total SM",true);
  TGraphAsymmErrors* metsig350new  =(TGraphAsymmErrors*)  effgraph->Clone("metsig350new");
  //  drawEffVRej(name450,"totalsm","hh 450","total SM",true);
  //  TGraphAsymmErrors* metsig450new  =(TGraphAsymmErrors*)  effgraph->Clone("metsig450new");


  var="MET"; xtitle="MET";
  nbins=20; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  setLogY(true);
  setPlotMinimum(0.1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_optcuts_met",0);
  setLogY(false);
  drawEffVRej(name200,"totalsm","hh 200","total SM",true);
  TGraphAsymmErrors* met200  =(TGraphAsymmErrors*)  effgraph->Clone("met200");
 drawEffVRej(name350,"totalsm","hh 350","total SM",true);
  TGraphAsymmErrors* met350  =(TGraphAsymmErrors*)  effgraph->Clone("met350");
  // drawEffVRej(name450,"totalsm","hh 450","total SM",true);
  //  TGraphAsymmErrors* met450  =(TGraphAsymmErrors*)  effgraph->Clone("met450");

  renewCanvas();
  //  metsig350->SetMaximum(0.02);
  metsig350->SetLineColor(kRed);
  metsig350->SetMarkerColor(kRed);
  metsig350->Draw("APL");

  metsig350new->SetLineColor(kMagenta);
  metsig350new->SetMarkerColor(kMagenta);
  metsig350new->Draw("PL");

  met350->SetLineColor(kBlue);
  met350->SetMarkerColor(kBlue);
  met350->Draw("PL");

  met350->SetFillStyle(0);
  met350->SetFillColor(kWhite);
  metsig350->SetFillStyle(0);
  metsig350->SetFillColor(kWhite);
  metsig350new->SetFillStyle(0);
  metsig350new->SetFillColor(kWhite);

  resetLegendPositionR();
  renewLegend();
  leg->AddEntry(metsig350,"Old METsig");
  leg->AddEntry(metsig350new,"New METsig");
  leg->AddEntry(met350,"MET");
  leg->Draw();

  thecanvas->SaveAs("higgs_optim_METsigsCompared_350_2.eps");
  thecanvas->SaveAs("higgs_optim_METsigsCompared_350_2.png");
  thecanvas->SaveAs("higgs_optim_METsigsCompared_350_2.pdf");


 renewCanvas();
  //  metsig200->SetMaximum(0.02);
  metsig200->SetLineColor(kRed);
  metsig200->SetMarkerColor(kRed);
  metsig200->Draw("APL");

  metsig200new->SetLineColor(kMagenta);
  metsig200new->SetMarkerColor(kMagenta);
  metsig200new->Draw("PL");

  met200->SetLineColor(kBlue);
  met200->SetMarkerColor(kBlue);
  met200->Draw("PL");

  met200->SetFillStyle(0);
  met200->SetFillColor(kWhite);
  metsig200->SetFillStyle(0);
  metsig200->SetFillColor(kWhite);
  metsig200new->SetFillStyle(0);
  metsig200new->SetFillColor(kWhite);

  resetLegendPositionR();
  renewLegend();
  leg->AddEntry(metsig200,"Old METsig");
  leg->AddEntry(metsig200new,"New METsig");
  leg->AddEntry(met200,"MET");
  leg->Draw();

  thecanvas->SaveAs("higgs_optim_METsigsCompared_200_2.eps");
  thecanvas->SaveAs("higgs_optim_METsigsCompared_200_2.png");
  thecanvas->SaveAs("higgs_optim_METsigsCompared_200_2.pdf");

  /*
 renewCanvas();
  //  metsig450->SetMaximum(0.02);
  metsig450->SetLineColor(kRed);
  metsig450->SetMarkerColor(kRed);
  metsig450->Draw("APL");

  metsig450new->SetLineColor(kMagenta);
  metsig450new->SetMarkerColor(kMagenta);
  metsig450new->Draw("PL");

  met450->SetLineColor(kBlue);
  met450->SetMarkerColor(kBlue);
  met450->Draw("PL");

  met450->SetFillStyle(0);
  met450->SetFillColor(kWhite);
  metsig450->SetFillStyle(0);
  metsig450->SetFillColor(kWhite);
  metsig450new->SetFillStyle(0);
  metsig450new->SetFillColor(kWhite);

  resetLegendPositionR();
  renewLegend();
  leg->AddEntry(metsig450,"Old METsig");
  leg->AddEntry(metsig450new,"New METsig");
  leg->AddEntry(met450,"MET");
  leg->Draw();

  thecanvas->SaveAs("higgs_optim_METsigsCompared_450.eps");
  thecanvas->SaveAs("higgs_optim_METsigsCompared_450.png");
  thecanvas->SaveAs("higgs_optim_METsigsCompared_450.pdf");
  */

}


void higgsmbb_quickmbbPlots() {

  initHiggsSamples69(true,"");
  setOutputDirectory("quickmbbPlots");

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  setStackMode(true,false,false); // stack,norm,labels
  //20 June -- update
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0 && MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; //4july add full cleanup
  //20 june -- remove the all jet trigger: ||passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_PFMET150==1";
  TCut zl = "nMuons==0&&nElectrons==0";
  TCut tauveto="nTausLoose==0";
  TCut isotk1="nIsoTracks15_005_03==0";

  TCut njets="njets20==4 || njets20==5";
  TCut jet2="jetpt2>50";
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut mdp = "minDeltaPhi20>0.3"; //4 july

  TCut higgsmass = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsdiff = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";

  TCut drmax = "deltaRmax_hh<2.2";
  TCut drmin = "deltaRmin_hh<1.9";

  TCut metsig = "METsig>30";

  //no higgs mass cut
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2 && btag3 && btag4 && higgsdiff && drmax && metsig;

  clearSamples(); //SMS-HbbHbb_mHiggsino-350_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1811_v69
  addSample("SMS-HbbHbb_mHiggsino-350_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1811_v69:higgsMbb1MassDiff_correct==0",kRed,"hh 350 (0)");
  addSample("SMS-HbbHbb_mHiggsino-350_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1811_v69:higgsMbb1MassDiff_correct==1",kBlue,"hh 350 (1)");
  addSample("SMS-HbbHbb_mHiggsino-350_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1811_v69:higgsMbb1MassDiff_correct==2",kGreen+4,"hh 350 (2)");

  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{jj}> (GeV)";
  nbins=20; low= 0; high=200;
  setStackMode(false,false,false); 
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_massshape_mbb_350",0,"GeV");

  clearSamples(); //SMS-HbbHbb_mHiggsino-200_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1807_v69
  addSample("SMS-HbbHbb_mHiggsino-200_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1807_v69:higgsMbb1MassDiff_correct==0",kRed,"hh 200 (0)");
  addSample("SMS-HbbHbb_mHiggsino-200_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1807_v69:higgsMbb1MassDiff_correct==1",kBlue,"hh 200 (1)");
  addSample("SMS-HbbHbb_mHiggsino-200_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1807_v69:higgsMbb1MassDiff_correct==2",kGreen+4,"hh 200 (2)");
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_massshape_mbb_200",0,"GeV");



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
