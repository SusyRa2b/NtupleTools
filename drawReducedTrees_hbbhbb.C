/*
====== this is the nominal code for drawing RA2b data/MC comparison plots =======

The functions that do the heavy lifting, along with various utility functions,
are in drawReducedTrees.h.

Suggested version of ROOT is ROOT 5.34/04
/afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.04/x86_64-slc5-gcc43-opt/root/bin/thisroot.csh
ROOT 5.32 will not work.

--- setup --

You must have symlink to or copy of MiscUtil.cxx in the working directory.
This is available at:
https://github.com/joshmt1/UserCode/blob/master/MiscUtil.cxx

You must also make a symlink as follows:
ln -s drawReducedTrees.h drawReducedTrees_hbbhbb.h

--- to compile ---

.L TSelectorMultiDraw.C+
.L CrossSectionTable.cxx+
.L ConfigurationDescriptions.cxx+

.L drawReducedTrees_hbbhbb.C+

-- input file locations --

The paths to the reducedTrees are defined at the top of initHiggsSamples69()

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

TString inputPath = "";
TString dataInputPath =  "";

double lumiScale_ = 19306.647; //updated with pixel lumi (oct 2013)

//make a symlink that point from this name to drawReducedTree.h
//this is to make the ROOT dictionary generation work correctly
#include "drawReducedTrees_hbbhbb.h"

TString addEnding(const TString & samplename, bool useSkim) {

  const int version =71;
  const  TString skimType="owen"; // are we using the owen skim or the david skim/slim?

  TString cfAskim="";
  if (skimType=="david" && useSkim) cfAskim="ra2b";

  TString tailpiece="";
  if (useSkim && skimType=="owen") tailpiece = "-skim";
  else if (useSkim && skimType=="david") tailpiece = "s";

  TString out;
  out.Form("%s%s_v%d%s",samplename.Data(),cfAskim.Data(),version,tailpiece.Data());
  return out;
}

void initHiggsSamples69(const bool useSkim=true,const TString samplelist="") { 
  
  //  inputPath="/cu6/joshmt/reducedTrees/v71_5b/"; 
  //  dataInputPath="/cu6/joshmt/reducedTrees/v71_5b/data/"; //for skimmed data

  inputPath="/cu2/ra2b/reducedTrees/v71_5_pj/";
  dataInputPath=inputPath;

  if (useSkim)   {
    inputPath+="skimmed/";
    dataInputPath=inputPath;
  }
    

  samplesAll_.clear(); // hack to suppress the amount of output we get. not needed for proper functioning.

  //go back to using the unskimmed cfA+Owen's skim
  //i'm going to try to make this more modular
  addToSamplesAll(addEnding("BJets_HT-1000ToInf_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1895",useSkim));
  addToSamplesAll(addEnding("BJets_HT-250To500_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1893",useSkim));
  addToSamplesAll(addEnding("BJets_HT-500To1000_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1894",useSkim));
  
  addToSamplesAll(addEnding("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1903",useSkim));
  addToSamplesAll(addEnding("QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1897",useSkim));
  addToSamplesAll(addEnding("QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1904",useSkim));
  addToSamplesAll(addEnding("QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1898",useSkim));
  addToSamplesAll(addEnding("QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1905",useSkim));
  addToSamplesAll(addEnding("QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1899",useSkim));
  addToSamplesAll(addEnding("QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1900",useSkim));
  addToSamplesAll(addEnding("QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1901",useSkim));
  addToSamplesAll(addEnding("QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1902",useSkim));
  addToSamplesAll(addEnding("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883",useSkim));
  addToSamplesAll(addEnding("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880",useSkim));
  addToSamplesAll(addEnding("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884",useSkim));
  addToSamplesAll(addEnding("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V19_ext1-v1_AODSIM_UCSB1962",useSkim));
  addToSamplesAll(addEnding("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V19_ext2-v1_AODSIM_UCSB1959",useSkim));
  
  addToSamplesAll(addEnding("W2JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1877",useSkim));
  addToSamplesAll(addEnding("W3JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1878",useSkim));
  addToSamplesAll(addEnding("W4JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1879",useSkim));
  addToSamplesAll(addEnding("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1887",useSkim));
  addToSamplesAll(addEnding("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1888",useSkim));
  addToSamplesAll(addEnding("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1889",useSkim));
  addToSamplesAll(addEnding("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1890",useSkim));
  addToSamplesAll(addEnding("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1891",useSkim));
  addToSamplesAll(addEnding("T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1860",useSkim));
  addToSamplesAll(addEnding("T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1861",useSkim));
  addToSamplesAll(addEnding("T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1862",useSkim));
  addToSamplesAll(addEnding("Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1864",useSkim));
  addToSamplesAll(addEnding("Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1865",useSkim));
  addToSamplesAll(addEnding("Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1866",useSkim));

    //official MG signal
  addToSamplesAll(addEnding("SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872",useSkim));
  addToSamplesAll(addEnding("SMS-TChiHH_2b2b_2J_mChargino-350to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871",useSkim));
		  //  addToSamplesAll("SMS-TChiZH_ZccbbHbb_2J_mChargino-130to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1873",useSkim));

  //fullsim signal
  addToSamplesAll("HbbHbb_FullSim_400_v71");
  addToSamplesAll("HbbHbb_FullSim_250_v71");

  //  addToSamplesAll("TT_8TeV-mcatnlo_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1885ra2b_v71s");

  //sherpa
  //  addToSamplesAll("TTJets_SemiLeptDecays_8TeV-sherpa_Summer12_DR53X-PU_S10_START53_V19-v1_AODSIM_UCSB1892ra2b_v71s");
  //  addToSamplesAll(TString("TTJets_DileptDecays_8TeV-sherpa_Summer12_DR53X-PU_S10_START53_V19-v1_AODSIM_UCSB1804_v69"));

  //ttV
  //  addToSamplesAll("TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1857ra2b_v71s");
  //  addToSamplesAll("TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1856ra2b_v71s");
  //  addToSamplesAll("TTH_HToBB_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1855ra2b_v71s");

  //  addToSamplesAll("WbbJetsToLNu_Massive_TuneZ2star_8TeV-madgraph-pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1859ra2b_v71s");

  //vv
  //  addToSamplesAll("WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1874ra2b_v71s");
  //  addToSamplesAll("WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1875ra2b_v71s");
  //  addToSamplesAll("ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1876ra2b_v71s");
  //  addToSamplesAll("WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1858ra2b_v71s");

  //fastsim ttbar
  //  addToSamplesAll("TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V7C_FSIM-v1_AODSIM_UCSB1882ra2b_v71s");

  nameOfEventWeight_="weight3"; 

  loadSamples("hbbhbb");
  usePUweight_=true;

  useTrigEff_=true;
  //for METsig <30 use the same weight as 30-50.
  trigEffWeight_ = "*(0.804*(METsig<50)+0.897*(METsig>=50&&METsig<100)+0.944*(METsig>=100))"; //must include the leading *
  currentConfig_=configDescriptions_.getDefault();

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;
  clearSamples();
  resetChains(); //important in the case that we call initHiggsSamples() more than once

  std::cout << "clear samples reset chains chainsamples" << std::endl;

  if (samplelist=="" || samplelist.Contains("qcd")) {
      TString baseqcdname = addEnding("QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1897",useSkim);
      addSample(baseqcdname,kGreen,"QCD");
      chainSamples(baseqcdname,addEnding("QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1898",useSkim));
      chainSamples(baseqcdname,addEnding("QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1899",useSkim));
      chainSamples(baseqcdname,addEnding("QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1900",useSkim));
      chainSamples(baseqcdname,addEnding("QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1901",useSkim));
      chainSamples(baseqcdname,addEnding("QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1902",useSkim));
      chainSamples(baseqcdname,addEnding("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1903",useSkim));
      chainSamples(baseqcdname,addEnding("QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1904",useSkim));
      chainSamples(baseqcdname,addEnding("QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1905",useSkim));
  }

  if (samplelist.Contains("bjets")) {
    TString     bname=addEnding("BJets_HT-250To500_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1893",useSkim);
    addSample(bname,kGreen,"QCD"); //used to be called 'BJets'
    chainSamples(bname,addEnding("BJets_HT-1000ToInf_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1895",useSkim));
    chainSamples(bname,addEnding("BJets_HT-500To1000_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1894",useSkim));
    setSampleScaleFactor(bname,1.3);
  }

  if (samplelist==""||samplelist.Contains("ttbar")) {
    addSample(addEnding("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880",useSkim),kGreen-9,"t#bar{t} (0l)");
  }

/*
  //no unskimmed version right now
  if (samplelist==""||samplelist.Contains("ttv")) {
    TString  ttvname="TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1857ra2b_v71s";
    addSample(ttvname,kCyan-7,"TTV");
    chainSamples(ttvname,"TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1856ra2b_v71s");
    chainSamples(ttvname,"TTH_HToBB_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1855ra2b_v71s");
  }
  //no unskimmed version right now
  if (samplelist==""||samplelist.Contains("vv")) {
    TString  vvname="WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1874ra2b_v71s";
    addSample(vvname,kCyan+3,"VV");
    chainSamples(vvname,"WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1875ra2b_v71s");
    chainSamples(vvname,"ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1876ra2b_v71s");
    chainSamples(vvname,"WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1858ra2b_v71s");
  }
*/

  if (samplelist==""||samplelist.Contains("singlet")) {
    TString  stname=addEnding("T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1860",useSkim);
    addSample(stname,kBlue,"Single Top");
    chainSamples(stname,addEnding("T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1861",useSkim));
    chainSamples(stname,addEnding("T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1862",useSkim));
    chainSamples(stname,addEnding("Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1864",useSkim));
    chainSamples(stname,addEnding("Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1865",useSkim));
    chainSamples(stname,addEnding("Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1866",useSkim));
  }


  if (samplelist=="" ||samplelist.Contains("wjets")) {
    TString w2=addEnding("W2JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1877",useSkim);
    addSample(w2,kViolet,"W+jets");
    chainSamples(w2,addEnding("W3JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1878",useSkim));
    chainSamples(w2,addEnding("W4JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1879",useSkim));
  }

/*
  //no unskimmed version right now
  if (samplelist.Contains("wbb")) {
    addSample("WbbJetsToLNu_Massive_TuneZ2star_8TeV-madgraph-pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1859ra2b_v71s",kViolet+2,"W+bb");
  }
*/

  if (samplelist==""||samplelist.Contains("znunu")) {
    TString basezname = addEnding("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1887",useSkim);
    addSample(basezname,kOrange,"Z #rightarrow #nu#nu");
    chainSamples(basezname,addEnding("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1889",useSkim));
    chainSamples(basezname,addEnding("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1888",useSkim));
    chainSamples(basezname,addEnding("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1891",useSkim));
    chainSamples(basezname,addEnding("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1890",useSkim));

  }


  if (samplelist==""||samplelist.Contains("ttbar")) {
    addSample(addEnding("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883",useSkim),kBlue-10,"t#bar{t} (2l)");
    TString tslname = addEnding("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884",useSkim);
    addSample(tslname,kAzure-3,"t#bar{t} (1l)");
    chainSamples(tslname,addEnding("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V19_ext1-v1_AODSIM_UCSB1962",useSkim));
    chainSamples(tslname,addEnding("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V19_ext2-v1_AODSIM_UCSB1959",useSkim));
  }


  const float hbbbb=0.561*0.561;

  if (samplelist.Contains("hhmg175") || samplelist=="")
    addSample(addEnding("SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872",useSkim)+"$175$1",kRed-5,"m_{#tilde{H}} = 175 GeV");

  if (samplelist.Contains("hhmg200") || samplelist=="") //used to be kRed-9
    addSample(addEnding("SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872",useSkim)+"$200$1",kRed  ,"m_{#tilde{H}} = 200 GeV");

  if (samplelist.Contains("hhmg250") || samplelist=="")
    addSample(addEnding("SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872",useSkim)+"$250$1",kRed-7,"m_{#tilde{H}} = 250 GeV");

  if (samplelist.Contains("hhmg300") || samplelist=="")
    addSample(addEnding("SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872",useSkim)+"$300$1",kRed-4,"m_{#tilde{H}} = 300 GeV");

  if (samplelist.Contains("hhmg350") || samplelist=="")
    addSample(addEnding("SMS-TChiHH_2b2b_2J_mChargino-350to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871",useSkim)+"$350$1",kRed+1,"m_{#tilde{H}} = 350 GeV");

  if (samplelist.Contains("hhmg400") || samplelist=="")
    addSample(addEnding("SMS-TChiHH_2b2b_2J_mChargino-350to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871",useSkim)+"$400$1",kRed+2,"m_{#tilde{H}} = 400 GeV");

  if (samplelist.Contains("hhmg450") || samplelist=="")
    addSample(addEnding("SMS-TChiHH_2b2b_2J_mChargino-350to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871",useSkim)+"$450$1",kRed+3,"m_{#tilde{H}} = 450 GeV");

  if (samplelist.Contains("hhmg500") || samplelist=="")
    addSample(addEnding("SMS-TChiHH_2b2b_2J_mChargino-350to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871",useSkim)+"$500$1",kRed+4,"m_{#tilde{H}} = 500 GeV");


  if (samplelist.Contains("hhfull250")) {
    loadReferenceCrossSections();
    addSample("HbbHbb_FullSim_250_v71",kCyan,"Full 250");
    double sf250 = hbbbb*CrossSectionTable_higgsino_->getSMSCrossSection(250)/double(getTree("HbbHbb_FullSim_250_v71")->GetEntries());
    cout<<" HH 250 scale factor = "<<hbbbb<<" * "<<CrossSectionTable_higgsino_->getSMSCrossSection(250)<<" / "<<getTree("HbbHbb_FullSim_250_v71")->GetEntries()<<endl;
    setSampleScaleFactor("HbbHbb_FullSim_250_v71",sf250);
  }
  if (samplelist.Contains("hhfull400")) {
    loadReferenceCrossSections();
    addSample("HbbHbb_FullSim_400_v71",kBlue,"Full 400");
    double sf400 = hbbbb*CrossSectionTable_higgsino_->getSMSCrossSection(400)/double(getTree("HbbHbb_FullSim_400_v71")->GetEntries());
    cout<<" HH 400 scale factor = "<<hbbbb<<" * "<<CrossSectionTable_higgsino_->getSMSCrossSection(400)<<" / "<<getTree("HbbHbb_FullSim_400_v71")->GetEntries()<<endl;
    setSampleScaleFactor("HbbHbb_FullSim_400_v71",sf400);
  }


  //scale by h->bb squared factor
  if (samplelist.Contains("hhmg") || samplelist=="") {
    setSampleScaleFactor(addEnding("SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872",useSkim),hbbbb);
    setSampleScaleFactor(addEnding("SMS-TChiHH_2b2b_2J_mChargino-350to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871",useSkim),hbbbb);
  }

  overrideSMSlabels_=true; //don't use the auto-labels for SMS

  setDatasetToDraw("MET"); //when plotting data, use the MET PD
}

void compareWbcMassShape() {
  initHiggsSamples69(true,"ttbar");

  clearSamples();
  addSample("TT_8TeV-mcatnlo_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1885ra2b_v71s:ttbarDecayCode!=2&&hadWcode==11000",kRed,"W #rightarrow bc");
  addSample("TT_8TeV-mcatnlo_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1885ra2b_v71s:ttbarDecayCode!=2&&(hadWcode==11||hadWcode==1100)",kBlue,"Diagonal");

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;
  useTrigEff_=false;
  usePUweight_=true; 
  savePlots_=true;

  setStackMode(false,true,false); //stack,norm,label override
  //define useful cuts

  nbins=20; low=0; high=200;
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="av higgs mass";
  selection_ = "";
  drawPlots(var,nbins,low,high,xtitle,"Events", "massShape_Wbc_CKMdiag",0);
  

}

void SLtoZLbackgroundCheck() {

  const  bool useSkim=true;
  initHiggsSamples69(useSkim,"ttbar znunu bjets wjets singlet"); //use all samples except signal


  //use top pT weights
  setSampleWeightFactor(addEnding("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883",useSkim),"topPtWeightOfficial");
  setSampleWeightFactor(addEnding("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884",useSkim),"topPtWeightOfficial");
  setSampleWeightFactor(addEnding("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880",useSkim),   "topPtWeightOfficial");


  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoPFcands10_010==0";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0";

  TCut njets4="njets20==4"; //switch back to 20 GeV jets
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";

  TCut jets = (njets4||njets5) && jet2;

  //  TCut mdp = "minDeltaPhi20>0.5 || (minDeltaPhi20>0.3&&METsig>50)";
  TCut mdp = "minDeltaPhi20_eta5_noIdAll_nobeta>0.5 || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3 && METsig>50)";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;

  TCut metsig[4];
  metsig[0]="METsig>=30&&METsig<50";
  metsig[1]="METsig>=50&&METsig<100";
  metsig[2]="METsig>=100&&METsig<150";
  metsig[3]="METsig>=150";

  TCut drmax = "deltaRmax_hh<2.2";

  setStackMode(true,false,false); //stack,norm,label override
  setLogY(false);

  TCut btag2sf="nbtag2_nomSF==1";
  TCut btag3sf="nbtag3_nomSF==1";
  TCut btag4sf="nbtag4_nomSF==1";

  setLogY(false);
  savePlots_=false;
  /* 

-- in each METsig bin --

compute in MC: R01 = (4b ZL)/(4b SL) 
count in data: 4b SL
compute: Npred = R01 * (4b SL)_data

compare to Nobs (4b ZL)_data

  */

  double nSL_MC[3][4];
  double nZL_MC[3][4];
  double nSL_MC_err[3][4];
  double nZL_MC_err[3][4];
  double nSL_data[3][4];
  double nZL_data[3][4];

  TCut btagcuts[3];
  btagcuts[0] = btag2sf;
  btagcuts[1] = btag3sf;
  btagcuts[2] = btag4sf;

  for (int ib = 0; ib<3; ib++) {
    for (int imet = 0; imet<4; imet++) {
      
      nbins=10; low=0; high=1e6;
      var="HT"; xtitle="HT";
 
      // SL selection ; 4b SIG
      selection_ = baseline && trigger&& sl && jets && mdp && higgsSR && metsig[imet] && btagcuts[ib] && drmax;
      
      drawPlots(var,nbins,low,high,xtitle,"Events", "blah",0);
      nSL_MC[ib][imet] = getIntegral("totalsm");
      nSL_MC_err[ib][imet] = getIntegralErr("totalsm");
      nSL_data[ib][imet] = getIntegral("data");
      
      // ZL selection ; 4b SIG
      selection_ = baseline && trigger&& zl && isotk && jets && mdp && higgsSR && metsig[imet] && btagcuts[ib] && drmax;
      drawPlots(var,nbins,low,high,xtitle,"Events", "blah",0);
      nZL_MC[ib][imet] = getIntegral("totalsm");
      nZL_MC_err[ib][imet] = getIntegralErr("totalsm");
      nZL_data[ib][imet] = getIntegral("data");
      
    }
  }

  printf("nb ; S bin ; ZL MC / SL MC ; R01  ; SL data ; prediction  ; observed\n");
  for (int ib = 0; ib<3; ib++) {
    cout<< " --------------------------------------------"<<endl;
    for (int imet = 0; imet<4; imet++) {
      //compute R01
    double R01 = nZL_MC[ib][imet] / nSL_MC[ib][imet];
    double R01_err = jmt::errAoverB( nZL_MC[ib][imet],nZL_MC_err[ib][imet],nSL_MC[ib][imet],nSL_MC_err[ib][imet]);

    double prediction = R01 * nSL_data[ib][imet];
    double prediction_err = jmt::errAtimesB( R01,R01_err,nSL_data[ib][imet],sqrt(nSL_data[ib][imet]));

    cout<<ib+2<<" "<<imet+1<<" ;  "<<nZL_MC[ib][imet]<< " / "<<nSL_MC[ib][imet]<<" = "<<R01<<" +/- "<<R01_err<<" ; "<<nSL_data[ib][imet]<<" ; "<<jmt::format_nevents(prediction,prediction_err,true)<<" ; "<<nZL_data[ib][imet]<<endl;
    
    }
  }


}

void SimpleABCDcheck(bool usedata=false) {
  initHiggsSamples69(true,"bjets znunu ttbar");
  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0.4; ratioMax = 1.6;
  dodata_=usedata;
  useTrigEff_=false;
  usePUweight_=true; 
  savePlots_=false;

  setStackMode(false,false,false); //stack,norm,label override
  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 
  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoPFcands10_010==0";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut mdp = "minDeltaPhi20_eta5_noIdAll_nobeta>0.5 || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3 && METsig>50)";
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut btag2sf="nbtag2_nomSF==1";
  TCut btag3sf="nbtag3_nomSF==1";
  TCut btag4sf="nbtag4_nomSF==1";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;

  TCut metsigveryloose="METsig>30";

  TCut drmax = "deltaRmax_hh<2.2";

  TCut SIG_4b = baseline&&trigger&&zl&&isotk&&(njets4||njets5)&&jet2&&mdp&&drmax&&btag4sf && higgsSR;
  TCut SB_4b = baseline&&trigger&&zl&&isotk&&(njets4||njets5)&&jet2&&mdp&&drmax&& btag4sf && higgsSB;
  TCut SIG_3b = baseline&&trigger&&zl&&isotk&&(njets4||njets5)&&jet2&&mdp&&drmax&&btag3sf && higgsSR;
  TCut SB_3b = baseline&&trigger&&zl&&isotk&&(njets4||njets5)&&jet2&&mdp&&drmax&& btag3sf && higgsSB;
  TCut SIG_2b = baseline&&trigger&&zl&&isotk&&(njets4||njets5)&&jet2&&mdp&&drmax&&btag2sf  && higgsSR;
  TCut SB_2b = baseline&&trigger&&zl&&isotk&&(njets4||njets5)&&jet2&&mdp&&drmax&&btag2sf  && higgsSB;

  float varbinsmetsig[6] = {0,30,50,100,150,500};
  nbins=5; low=0; high=500;
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";

  TH1D* h1=0;

  selection_ = SIG_4b;
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",varbinsmetsig);
  h1 = usedata ? hdata : totalsm;
  TH1D* h_SIG_4b = (TH1D*) h1->Clone("h_SIG_4b");

  selection_ = SB_4b;
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",varbinsmetsig);
  h1 = usedata ? hdata : totalsm;
  TH1D* h_SB_4b = (TH1D*) h1->Clone("h_SB_4b");

  selection_ = SIG_3b;
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",varbinsmetsig);
  h1 = usedata ? hdata : totalsm;
  TH1D* h_SIG_3b = (TH1D*) h1->Clone("h_SIG_3b");

  selection_ = SB_3b;
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",varbinsmetsig);
  h1 = usedata ? hdata : totalsm;
  TH1D* h_SB_3b = (TH1D*) h1->Clone("h_SB_3b");

  selection_ = SIG_2b;
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",varbinsmetsig);
  h1 = usedata ? hdata : totalsm;
  TH1D* h_SIG_2b = (TH1D*) h1->Clone("h_SIG_2b");

  selection_ = SB_2b;
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",varbinsmetsig);
  h1 = usedata ? hdata : totalsm;
  TH1D* h_SB_2b = (TH1D*) h1->Clone("h_SB_2b");

  TString truthlabel = dodata_ ? " Observed = " : " Truth      = ";

  // pred = SB_4b * SIG_2b / SB_2b
  for (int ibin =2; ibin<=5; ibin++) { //skip first bin
    cout<<" == ABCD for bin starting at "<<h_SIG_4b->GetBinLowEdge(ibin)<<endl;

    cout<<" Prediction = "<<h_SB_4b->GetBinContent(ibin) * h_SIG_2b->GetBinContent(ibin) / h_SB_2b->GetBinContent(ibin);
    cout<<" +/- ";
    double Rerr=jmt::errAoverB(h_SIG_2b->GetBinContent(ibin),h_SIG_2b->GetBinError(ibin), h_SB_2b->GetBinContent(ibin),h_SB_2b->GetBinError(ibin));
    double err = jmt::errAtimesB(h_SB_4b->GetBinContent(ibin),h_SB_4b->GetBinError(ibin),h_SIG_2b->GetBinContent(ibin) / h_SB_2b->GetBinContent(ibin),Rerr);
    cout<<err<<endl;
    
    cout<<truthlabel<<h_SIG_4b->GetBinContent(ibin)<<" +/- "<<h_SIG_4b->GetBinError(ibin)<<endl;
  }

  // pred = SB_3b * SIG_2b / SB_2b
  for (int ibin =2; ibin<=5; ibin++) { //skip first bin
    cout<<" == ABCD for bin starting at "<<h_SIG_3b->GetBinLowEdge(ibin)<<endl;

    cout<<" Prediction = "<<h_SB_3b->GetBinContent(ibin) * h_SIG_2b->GetBinContent(ibin) / h_SB_2b->GetBinContent(ibin);
    cout<<" +/- ";
    double Rerr=jmt::errAoverB(h_SIG_2b->GetBinContent(ibin),h_SIG_2b->GetBinError(ibin), h_SB_2b->GetBinContent(ibin),h_SB_2b->GetBinError(ibin));
    double err = jmt::errAtimesB(h_SB_3b->GetBinContent(ibin),h_SB_3b->GetBinError(ibin),h_SIG_2b->GetBinContent(ibin) / h_SB_2b->GetBinContent(ibin),Rerr);
    cout<<err<<endl;
    
    cout<<truthlabel<<h_SIG_3b->GetBinContent(ibin)<<" +/- "<<h_SIG_3b->GetBinError(ibin)<<endl;
  }

}

void wjetsCompare() {

  initHiggsSamples69(true,"wjets wbb");
  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0.4; ratioMax = 1.6;
  dodata_=false;
  useTrigEff_=false;
  usePUweight_=true; 

  setStackMode(false,false,false); //stack,norm,label override
  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  //cannot use this with MET PD only
  //  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut mdp = "minDeltaPhi20>0.5 || (minDeltaPhi20>0.3&&METsig>50)";
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

  TCut drmax = "deltaRmax_hh<2.2";

  selection_=baseline&&zl&&isotk&&(njets4||njets5)&&jet2&&mdp&&btag2;

  nbins=30; low=0; high=300;
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "wjets_validation_METsig",0);

  selection_=baseline&&jet2&&btag2;
  nbins=8; low=0; high=8;
  var="njets20"; xtitle="jet multiplicity";
  drawPlots(var,nbins,low,high,xtitle,"Events", "wjets_validation_njets20",0);

  selection_=baseline&&jet2&&btag2&&TCut("njets20>=4");
  nbins=40; low=0; high=1;
  var="CSVbest3"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "wjets_validation_CSVbest3",0);


  selection_=baseline&&jet2&&btag2&&TCut("njets20>=4")&&btag3;
  nbins=40; low=0; high=1;
  var="CSVbest4"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "wjets_validation_CSVbest4",0);


}

void fullsimSignalCheck() {

  initHiggsSamples69(true,"hhmg250 hhmg400 hhfull250 hhfull400");
  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  usePUweight_=false; 

  useTrigEff_=false;

  //400 gev
  setSampleWeightFactor("HbbHbb_FullSim_400_v71","0 + (nGoodPV==1)*0.000000 + (nGoodPV==2)*0.000000 + (nGoodPV==3)*0.806114 + (nGoodPV==4)*0.550517 + (nGoodPV==5)*0.784298 + (nGoodPV==6)*0.850867 + (nGoodPV==7)*0.765688 + (nGoodPV==8)*0.724206 + (nGoodPV==9)*0.807895 + (nGoodPV==10)*0.839914 + (nGoodPV==11)*1.055339 + (nGoodPV==12)*1.135399 + (nGoodPV==13)*1.327954 + (nGoodPV==14)*1.712817 + (nGoodPV==15)*1.617481 + (nGoodPV==16)*1.299388 + (nGoodPV==17)*0.995578 + (nGoodPV==18)*0.771311 + (nGoodPV==19)*0.800530 + (nGoodPV==20)*0.798977 + (nGoodPV==21)*0.899122 + (nGoodPV==22)*0.892036 + (nGoodPV==23)*0.901595 + (nGoodPV==24)*0.996259 + (nGoodPV==25)*0.928075 + (nGoodPV==26)*0.766434 + (nGoodPV==27)*0.988833 + (nGoodPV==28)*0.599485 + (nGoodPV==29)*0.983066 + (nGoodPV==30)*0.941055 + (nGoodPV==31)*1.057779 + (nGoodPV==32)*0.885022 + (nGoodPV==33)*1.278620 + (nGoodPV==34)*3.954734 + (nGoodPV==35)*18.560288 + (nGoodPV==36)*0.000000 + (nGoodPV==37)*0.000000 + (nGoodPV==38)*0.000000 + (nGoodPV==39)*0.000000 + (nGoodPV==40)*0.000000 + (nGoodPV==41)*0.000000 + (nGoodPV==42)*0.000000 + (nGoodPV==43)*0.000000 + (nGoodPV==44)*0.000000 + (nGoodPV==45)*0.000000 + (nGoodPV==46)*0.000000 + (nGoodPV==47)*0.000000 + (nGoodPV==48)*0.000000 + (nGoodPV==49)*0.000000 + (nGoodPV==50)*0.000000 + (nGoodPV==51)*0.000000");
  //250 gev
  setSampleWeightFactor("HbbHbb_FullSim_250_v71","0 + (nGoodPV==1)*0.000000 + (nGoodPV==2)*0.000000 + (nGoodPV==3)*0.189076 + (nGoodPV==4)*0.476019 + (nGoodPV==5)*0.774264 + (nGoodPV==6)*0.855738 + (nGoodPV==7)*0.789157 + (nGoodPV==8)*0.701390 + (nGoodPV==9)*0.765677 + (nGoodPV==10)*0.862419 + (nGoodPV==11)*1.043154 + (nGoodPV==12)*1.126524 + (nGoodPV==13)*1.391946 + (nGoodPV==14)*1.666829 + (nGoodPV==15)*1.629645 + (nGoodPV==16)*1.286200 + (nGoodPV==17)*0.984172 + (nGoodPV==18)*0.789754 + (nGoodPV==19)*0.782547 + (nGoodPV==20)*0.795632 + (nGoodPV==21)*0.932901 + (nGoodPV==22)*0.861056 + (nGoodPV==23)*0.836711 + (nGoodPV==24)*0.979746 + (nGoodPV==25)*0.979451 + (nGoodPV==26)*0.830089 + (nGoodPV==27)*1.051759 + (nGoodPV==28)*0.649837 + (nGoodPV==29)*0.892681 + (nGoodPV==30)*1.070287 + (nGoodPV==31)*1.196061 + (nGoodPV==32)*1.043479 + (nGoodPV==33)*1.402871 + (nGoodPV==34)*1.958256 + (nGoodPV==35)*13.388514 + (nGoodPV==36)*19.901845 + (nGoodPV==37)*0.000000 + (nGoodPV==38)*0.000000 + (nGoodPV==39)*0.000000 + (nGoodPV==40)*0.000000 + (nGoodPV==41)*0.000000 + (nGoodPV==42)*0.000000 + (nGoodPV==43)*0.000000 + (nGoodPV==44)*0.000000 + (nGoodPV==45)*0.000000 + (nGoodPV==46)*0.000000 + (nGoodPV==47)*0.000000 + (nGoodPV==48)*0.000000 + (nGoodPV==49)*0.000000 + (nGoodPV==50)*0.000000 + (nGoodPV==51)*0.000000");

  setOutputDirectory("FastFullValidation");

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 
  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoPFcands10_010==0";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0"; //remove veto of 2nd isotk

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut jets=(njets4||njets5) &&jet2;
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut mdp = "minDeltaPhi20_eta5_noIdAll_nobeta>0.5 || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3&&METsig>50)";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;
  TCut drmax = "deltaRmax_hh<2.2";

  TCut metsigloose = "METsig>30";

  float metsigbins[5] = {30,50,100,150,200};

  //honestly don't know if i should normalize or not. do we care about shapes or normalization or both?

  //selection with *no* trigger requirement and *no* b-tagging (and no METsig)
  selection_ = baseline && zl&&isotk&& (njets4||njets5)&&jet2&& higgsSR && drmax ;

  //MET spectrum
  nbins=40; low=0; high=400;
  var="MET"; xtitle="MET";

  setStackMode(false,false,false); //stack,norm,label override
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_FullFast_MET",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_FullFast_MET",0);

  setStackMode(false,true,false); //stack,norm,label override
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_FullFast_MET",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_FullFast_MET",0);


  //METsig spectrum
  nbins=30; low=0; high=300;
  var="METsig"; xtitle="METsig";

  setStackMode(false,false,false); //stack,norm,label override
  setLogY(false);setPlotMinimum(0);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_FullFast_METsig",0);
  drawPlots(var,4,low,high,xtitle,"Events","higgs_FullFast_METsig_bins",metsigbins);
  resetPlotMinimum();
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_FullFast_METsig",0);

  setStackMode(false,true,false); //stack,norm,label override
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_FullFast_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_FullFast_METsig",0);

  resetSampleWeightFactors();

  setStackMode(false,false,false); //stack,norm,label override
  setLogY(false);setPlotMinimum(0);
  drawPlots(var,4,low,high,xtitle,"Events","higgs_FullFast_METsig_bins-noPUreweight",metsigbins);
 

  //now try in bins of PV

  //nb -- this doesn't make sense unless the plots are normalized to unit area
  resetSampleWeightFactors();

  selection_ = baseline && zl&&isotk&& (njets4||njets5)&&jet2&& higgsSR && drmax &&TCut("nGoodPV>=10&&nGoodPV<=12") &&TCut("METsig>30");
  setStackMode(false,true,false); //stack,norm,label override
  setLogY(false);
  setPlotMinimum(0);
  drawPlots(var,4,low,high,xtitle,"Events","higgs_FullFast_METsig_bins_lowPU",metsigbins);
  resetPlotMinimum();

  selection_ = baseline && zl&&isotk&& (njets4||njets5)&&jet2&& higgsSR && drmax &&TCut("nGoodPV>=19&&nGoodPV<=21");
  setStackMode(false,true,false); //stack,norm,label override
  setLogY(false);  setPlotMinimum(0);

  drawPlots(var,4,low,high,xtitle,"Events","higgs_FullFast_METsig_bins_highPU",metsigbins);
  resetPlotMinimum();


}

void mgSignalCheck() {
  initHiggsSamples69(true,"hh200 hh350 hhmg200 hhmg350");
  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  usePUweight_=true; 

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

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  //cannot use this with MET PD only
  //  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut mdp = "minDeltaPhi20>0.5 || (minDeltaPhi20>0.3&&METsig>50)";
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

  TCut metsig="METsig>30";

  TCut drmax = "deltaRmax_hh<2.2";
  //  TCut drmin = "deltaRmin_hh<1.9";
  selection_ = baseline && trigger && zl&&isotk&& (njets4||njets5)&&jet2&& btag2&&btag3&&btag4 && higgsSR && drmax &&metsig;

  nbins=20; low=0; high=300;
  var="bjetpt1"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidation_signalLikeSelection_"+var,0);



  //now, another check
  selection_="(1)";
  clearSamples();
    addSample("HiggsinoNLSP_chargino130_to_500_bino1-PU_S10-TChihh_v69$200$1",kGreen,"hh MG 200");
    addSample("HiggsinoNLSP_chargino130_to_500_bino1-PU_S10-TChihh_v69$250$1",kBlue,"hh MG 250");
    addSample("HiggsinoNLSP_chargino130_to_500_bino1-PU_S10-TChihh_v69$300$1",kMagenta,"hh MG 300");
    addSample("HiggsinoNLSP_chargino130_to_500_bino1-PU_S10-TChihh_v69$350$1",kRed,"hh MG 350");
    addSample("HiggsinoNLSP_chargino130_to_500_bino1-PU_S10-TChihh_v69$500$1",kOrange,"hh MG 500");

  setStackMode(false,true,false); //stack,norm,label override
 
  nbins=40; low=0; high=400;
  var="MET"; xtitle="MET";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SignalValidationAll4_"+var,0);

  nbins=40; low=0; high=300;
  var="METsig"; xtitle="METsig";
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

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
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

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
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

void higgs_signal_sculpt() {

  initHiggsSamples69(false,"hh250 hh400");


  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
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

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
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

void higgs_dataMC_control_QCD_noskim(TString options="") {
  //plots that are looser than the skim

  initHiggsSamples69(false,"ttbar qcd wjets znunu singlet"); //must remove znunu for now

  //use top pT weights
  setSampleWeightFactor("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71","topPtWeight");
  setSampleWeightFactor("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71","topPtWeight");
  setSampleWeightFactor("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71",   "topPtWeight");

  TString outputdir = "plots_Control_QCD";

  if (options!="") {
    outputdir+="_";
    outputdir+=options;

    if (options.Contains("qcd1p5")) {//needs update to use this option
      setSampleScaleFactor("QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1814_v69",1.5);
    }
  }

  setOutputDirectory(outputdir);

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 3.2;
  dodata_=true;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 
  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut ra2btrigger = "passMC_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80==1 || passMC_DiCentralPFJet50_PFMET80==1";

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoPFcands10_010==0";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0"; //remove veto of 2nd isotk

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut jet2high="jetpt2>70";
  TCut jets=(njets4||njets5) &&jet2;
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  //TCut mdp = "minDeltaPhi20>0.5 || (minDeltaPhi20>0.3&&METsig>50)";
  TCut mdp = "minDeltaPhi20_eta5_noIdAll_nobeta>0.5 || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3&&METsig>50)";

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

  addVerticalLine(30);
  nbins=40; low=0; high=400;
  setLogY(false);
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigMET150_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigMET150_METsig",0);
  resetVerticalLine();


  nbins=40; low=0; high=700;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigMET150_MET",0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigMET150_MET",0,"GeV");

  //might be better to raise the MET cut above the turn-on curve
  selection_ = baseline && triggerMET && zl&&isotk && jets&&jet2 && !btag3 &&!mdp &&TCut("MET>230");

   addVerticalLine(30);
 nbins=40; low=0; high=400;
  setLogY(false);
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigMET150-MET230_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigMET150-MET230_METsig",0);
  resetVerticalLine();


    addVerticalLine(30);
 //for RA2b trigger, raise jetpt2 threshold
  selection_ = baseline && ra2btrigger && zl&&isotk && jets&&jet2high && !btag3 &&!mdp;

  nbins=40; low=0; high=400;
  setLogY(false);
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigJetMET_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigJetMET_METsig",0);
  resetVerticalLine();


  nbins=40; low=0; high=700;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigJetMET_MET",0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigJetMET_MET",0,"GeV");


  selection_ = baseline && ra2btrigger && zl&&isotk && jets&&jet2high && !btag3 &&!mdp&&TCut("MET>160");
    addVerticalLine(30);

   nbins=40; low=0; high=400;
  setLogY(false);
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigJetMET-MET160_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_lt3bLDP_trigJetMET-MET160_METsig",0);
  resetVerticalLine();

  // == look at anomalous events
  TCut weird = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& (MET/caloMET>=2 || maxTOBTECjetDeltaMult>=40)"; 

  selection_ = weird && ra2btrigger && zl&&isotk && jets&&jet2high && !btag3 &&TCut("MET>160");

  nbins=40; low=0; high=700;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_weird_trigJetMET_MET",0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_weird_trigJetMET_MET",0,"GeV");

    addVerticalLine(30);

  nbins=40; low=0; high=400;
  setLogY(false);
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_weird_trigJetMET_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_weird_trigJetMET_METsig",0);
  resetVerticalLine();

  TCut junk = "cutPV==1 &&passCleaning==0 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  selection_ = junk && ra2btrigger && zl&&isotk && jets&&jet2high && !btag3 &&TCut("MET>160");

  nbins=40; low=0; high=700;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_junk_trigJetMET_MET",0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_junk_trigJetMET_MET",0,"GeV");

    addVerticalLine(30);

  nbins=40; low=0; high=400;
  setLogY(false);
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_junk_trigJetMET_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_junk_trigJetMET_METsig",0);
  resetVerticalLine();


  // == now other plots (non-LDP)
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;

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

void higgs_dataMC_control_PUtests() {
  initHiggsSamples69(true,"ttbar znunu bjets wjets ttv vv singlet");

  setOutputDirectory("plots_Control_BJets-PUtest");


  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoPFcands10_010==0";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut jets=(njets4||njets5) &&jet2;
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut mdp = "minDeltaPhi20>0.5 || (minDeltaPhi20>0.3&&METsig>50)";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;

  TCut metsigloose = "METsig>30";

  TCut drmax = "deltaRmax_hh<2.2";

  //plot:
  //DeltaPhi(jets,MET)

  useTrigEff_=false;


  setSampleScaleFactor("BJets_HT-250To500_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1893ra2b_v71s",1.3);

  setStackMode(true,false,false); //stack,norm,label override
  setLogY(false);

  //a plot to compare to owen's selection
  selection_ = baseline && trigger && zl && isotk && jets && btag2 && !btag3 && mdp && drmax;
  //no higgs mass cuts, but do use DRmax
  nbins=50; low=0; high=50;
  var="nGoodPV"; xtitle=var;

  usePUweight_=false;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bWithDrMax_nGoodPV_noPUweight",0);

  usePUweight_=true;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bWithDrMax_nGoodPV_withPUweight",0);

  //now enrich it in qcd and minimize trigger problems
  selection_ = baseline && trigger && zl && isotk && jets && btag2 && !btag3 && !mdp && TCut("MET>150");
  usePUweight_=false;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bLDP-MET150_nGoodPV_noPUweight",0);
  usePUweight_=true;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bLDP-MET150_nGoodPV_withPUweight",0);


}

void datacounts() {

  initHiggsSamples69(true,"ttbar");
  //  setOutputDirectory("plots_Signal");

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoPFcands10_010==0";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut jets=(njets4||njets5) &&jet2;
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  //  TCut mdp = "minDeltaPhi20>0.5 || (minDeltaPhi20>0.3&&METsig>50)";
  TCut mdp = "minDeltaPhi20_eta5_noIdAll_nobeta>0.5 || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3 && METsig>50)";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;

  setQuiet(true);
  TCut drmax = "deltaRmax_hh<2.2";

  //plot:
  //DeltaPhi(jets,MET)

  setStackMode(true,false,false); //stack,norm,label override
  setLogY(false);

  savePlots_=false;

  TCut met[4] = {"METsig>=30 && METsig<50","METsig>=50 && METsig<100","METsig>=100 && METsig<150","METsig>=150"};
  TCut masswindows[2]={higgsSB,higgsSR}; //is this legal?
  TCut bcut[3]={btag2&&!btag3,btag2&&btag3&&!btag4,btag2&&btag3&&btag4};

  for (int ib=0;ib<3; ib++) {
    for (int imw=0;imw<2; imw++) {
      for (int imet=0;imet<4; imet++) {	
	selection_ = baseline && trigger && zl && isotk && jets && bcut[ib] && mdp && masswindows[imw] && drmax && met[imet];
	nbins=1; low=0; high=20000;
	var="HT"; xtitle = "e";
	drawPlots(var,nbins,low,high,xtitle,"Events","",0);
	cout<<"MET bin "<<imet+1<<" nb = "<<ib+2<<" isSIG "<<imw<<" counts = ";
	cout<<	hdata->Integral()<<endl;
      }
    }
  }

}

void higgs_dataMC_SR() {

  initHiggsSamples69(true,"ttbar znunu bjets wjets ttv vv singlet hhmg250 hhmg300 hhmg350 hhmg400");
  setOutputDirectory("plots_Signal");

  //use top pT weights
  setSampleWeightFactor("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883ra2b_v71s","topPtWeight");
  setSampleWeightFactor("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884ra2b_v71s","topPtWeight");
  setSampleWeightFactor("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880ra2b_v71s",   "topPtWeight");


  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoPFcands10_010==0";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut jets=(njets4||njets5) &&jet2;
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  //  TCut mdp = "minDeltaPhi20>0.5 || (minDeltaPhi20>0.3&&METsig>50)";
  TCut mdp = "minDeltaPhi20_eta5_noIdAll_nobeta>0.5 || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3 && METsig>50)";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;


  TCut drmax = "deltaRmax_hh<2.2";

  //plot:
  //DeltaPhi(jets,MET)

  setStackMode(true,false,false); //stack,norm,label override
  setLogY(false);

  TCut metloose = "METsig>=30";
  TCut met0 = "METsig>=15 && METsig<30";
  TCut met1 = "METsig>=30 && METsig<50";
  TCut met2 = "METsig>=50 && METsig<100";
  TCut met3 = "METsig>=100 && METsig<150";
  TCut met4 = "METsig>=150";

  selection_ = baseline && trigger && zl && isotk && jets && btag2 && btag3 && btag4 && mdp && higgsSR && drmax && metloose;
  nbins=20; low=0; high=200;
  var="METsig"; xtitle = "E^{miss}_{T} significance #it{S}";
  drawPlots(var,nbins,low,high,xtitle,"Events","SR4b_allcuts_METsig",0);


}

void higgs_dataMC_control_QCD() {

  //goal: use 2-b sample to get a qcd control region. 

  //  initHiggsSamples69(true,"ttbar znunu bjets wjets ttv vv singlet"); //all backgrounds
  //for now i only have unskimmed and i don't have every sample
  initHiggsSamples69(false,"ttbar znunu bjets wjets singlet"); //all backgrounds
  //  initHiggsSamples69(true,"ttbar znunu bjets wjets ttv vv singlet"); //all backgrounds

  setOutputDirectory("plots_Control_QCD");


  //use top pT weights
  setSampleWeightFactor("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883ra2b_v71s","topPtWeightOfficial");
  setSampleWeightFactor("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884ra2b_v71s","topPtWeightOfficial");
  setSampleWeightFactor("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880ra2b_v71s",   "topPtWeightOfficial");

  //use top pT weights
  setSampleWeightFactor("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71","topPtWeightOfficial");
  setSampleWeightFactor("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71","topPtWeightOfficial");
  setSampleWeightFactor("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71",   "topPtWeightOfficial");


  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoPFcands10_010==0";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut jets=(njets4||njets5) &&jet2;
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut btag2sf="nbtag2_nomSF==1";
  TCut btag3sf="nbtag3_nomSF==1";
  TCut btag4sf="nbtag4_nomSF==1";

  //  TCut mdp = "minDeltaPhi20>0.5 || (minDeltaPhi20>0.3&&METsig>50)";
  TCut mdp = "minDeltaPhi20_eta5_noIdAll_nobeta>0.5 || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3 && METsig>50)";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;

  TCut metsigloose = "METsig>30";

  TCut drmax = "deltaRmax_hh<2.2";

  //plot:
  //DeltaPhi(jets,MET)

  setStackMode(true,false,false); //stack,norm,label override
  setLogY(false);

  // === convert to using b-tag SF

  //a plot to compare to owen's selection
  selection_ = baseline && trigger && zl && isotk && jets && btag2sf && mdp && drmax;
  //no higgs mass cuts, but do use DRmax
  addVerticalLine(30);
  nbins=80; low=0; high=400;
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bWithDrMax_METsig",0);
  resetVerticalLine();

   //i think we want: preselection. 2b tags. lepton & isotk vetoes.
  selection_ = baseline && trigger && zl&&isotk && jets && btag2sf;

  addVerticalLine(30);
  nbins=30; low=0; high=300;
  setLogY(false);
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselection_METsig",0);
  setLogY(true); setPlotMinimum(0.5);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselection_METsig",0);
  resetVerticalLine();
  resetPlotMinimum();

  nbins=30; low=0; high=400;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselection_MET",0);
  setLogY(true); setPlotMinimum(0.5);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselection_MET",0);
  resetPlotMinimum();

  //let's see what these look like with a loose METsig cut
  //  selection_ = baseline && trigger && zl&&isotk && jets && btag2 && !btag3 && metsigloose;
  //use b-tag SF
  selection_ = baseline && trigger && zl&&isotk && jets && btag2sf && metsigloose;

  addVerticalLine(30);
  nbins=30; low=0; high=300;
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  setLogY(true); setPlotMinimum(0.5);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMETsig30_METsig",0);
  resetVerticalLine();
  resetPlotMinimum();

  nbins=30; low=0; high=400;
  var="MET"; xtitle="E^{miss}_{T} (GeV)";
  setLogY(true); setPlotMinimum(0.5);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMETsig30_MET",0,"GeV");
  resetPlotMinimum();

  
  //nPV
  nbins=20; low=0; high=40;
  setLogY(false);
  var="nGoodPV"; xtitle="n Good PV";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselection_nGoodPV",0);

  //same plots but adding mdp cut
  selection_ = baseline && trigger && zl&&isotk && jets && btag2sf &&mdp;

  addVerticalLine(30);
  nbins=30; low=0; high=300;
  setLogY(false);
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP_METsig",0);
  resetVerticalLine();

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



  /*
     Bill Gary:
     It might also be good to make distributions of MET for events
     in each of our four METsig bins
  */
  //keep the nominal very loose analysis selection but with 3rd b veto
  nbins=50; low=0; high=500;
  var="MET"; xtitle="MET (GeV)";


  selection_ = baseline && trigger && zl&&isotk && jets && btag2sf &&mdp &&TCut("METsig>=30 && METsig<50");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP-METsig30to50_MET",0);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP-METsig30to50_MET",0);

  selection_ = baseline && trigger && zl&&isotk && jets && btag2sf &&mdp &&TCut("METsig>=50 && METsig<100");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP-METsig50to100_MET",0);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP-METsig50to100_MET",0);

  selection_ = baseline && trigger && zl&&isotk && jets && btag2sf &&mdp &&TCut("METsig>=100 && METsig<150");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP-METsig100to150_MET",0);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP-METsig100to150_MET",0);

  selection_ = baseline && trigger && zl&&isotk && jets && btag2sf &&mdp &&TCut("METsig>=150");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP-METsig150toInf_MET",0);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP-METsig150toInf_MET",0);


  //same plots but inverting mdp cut
  selection_ = baseline && trigger && zl&&isotk && jets && btag2sf &&!mdp;

  addVerticalLine(30);

  nbins=13; low=0; high=300;
  float varbinsmetsig[14] = {0,10,20,30,40,50,60,70,80,100,120,140,180,300};
  setLogY(false);
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP_METsig",varbinsmetsig);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP_METsig",varbinsmetsig);

  renormalizeBins_=true; renormalizeWidth_=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP_METsig",varbinsmetsig);
  renormalizeBins_=false;

  resetVerticalLine();

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
  selection_ = baseline && trigger && zl&&isotk && jets && btag2sf &&mdp&&minmet;
  addVerticalLine(30);
  nbins=30; low=0; high=300;
  setLogY(false);
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP-MET125_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP-MET125_METsig",0);

  resetVerticalLine();

  //nPV
  nbins=20; low=0; high=40;
  setLogY(false);
  var="nGoodPV"; xtitle="n Good PV";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionMDP-MET125_nGoodPV",0);

  //same plots but inverting mdp cut
  selection_ = baseline && trigger && zl&&isotk && jets && btag2sf &&!mdp&&minmet;

  addVerticalLine(30);

  nbins=13; low=0; high=300;
  setLogY(false);
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP-MET125_METsig",varbinsmetsig);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP-MET125_METsig",varbinsmetsig);

  renormalizeBins_=true;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP-MET125_METsig",varbinsmetsig);
  renormalizeBins_=false;

  resetVerticalLine();

  //nPV
  nbins=20; low=0; high=40;
  setLogY(false);
  var="nGoodPV"; xtitle="n Good PV";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_QCDcontrol_preselectionLDP-MET125_nGoodPV",0);

  // == now other plots
  selection_ = baseline && trigger && zl && isotk && jets && btag2sf && metsigloose;
  nbins=30; low = 0; high=3;
  setLogY(false);
  var="minDeltaPhi20_eta5_noIdAll_nobeta"; xtitle="#Delta #phi_{min}(jet, MET)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bMETsig30_minDeltaPhi20uber",0);
  //this plot is worrisome in that there are extra events outside of the veto region

  selection_ = baseline && trigger && zl && isotk && jets && btag2sf && TCut("METsig>50 && METsig<=100");
  nbins=30; low = 0; high=3;
  setLogY(false);
  var="minDeltaPhi20_eta5_noIdAll_nobeta"; xtitle="#Delta #phi_{min}(jet, MET)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bMETsig50to100_minDeltaPhi20uber",0);

  selection_ = baseline && trigger && zl && isotk && jets && btag2sf && TCut("METsig>100 && METsig<=150");
  nbins=30; low = 0; high=3;
  setLogY(false);
  var="minDeltaPhi20_eta5_noIdAll_nobeta"; xtitle="#Delta #phi_{min}(jet, MET)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bMETsig100to150_minDeltaPhi20uber",0);

/* these were some studies...not really needed anymore
  nbins=30; low = 0; high=3;
  setLogY(false);
  var="deltaPhi1"; xtitle="#Delta #phi(jet 1 (50), MET)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bMETsig30_deltaPhi1",0);
  nbins=30; low = 0; high=3;
  setLogY(false);
  var="deltaPhi2"; xtitle="#Delta #phi(jet 2 (50), MET)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bMETsig30_deltaPhi2",0);
  nbins=30; low = 0; high=3;
  setLogY(false);
  var="deltaPhi3"; xtitle="#Delta #phi(jet 3 (50), MET)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bMETsig30_deltaPhi3",0);

  //these plots are a test
  selection_ = baseline && trigger && zl && isotk && jets && btag2 && !btag3 && metsigloose &&TCut("deltaPhi2>0.3");
  nbins=30; low = 0; high=3;
  setLogY(false);
  var="minDeltaPhi20"; xtitle="#Delta #phi_{min}(jet, MET)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bMETsig30_minDeltaPhi20-Dp2gt0p3",0);

  selection_ = baseline && trigger && zl && isotk && jets && btag2 && !btag3 && metsigloose &&TCut("deltaPhi2>0.5");
  nbins=30; low = 0; high=3;
  setLogY(false);
  var="minDeltaPhi20"; xtitle="#Delta #phi_{min}(jet, MET)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bMETsig30_minDeltaPhi20-Dp2gt0p5",0);
*/

  //now, we want: Data/MC minDeltaPhi distributions of 2b SB, 3b SB, 2b SIG; focus on 1st METsig bin
  nbins=30; low = 0; high=3;
  setLogY(false);
  var="minDeltaPhi20_eta5_noIdAll_nobeta"; xtitle="#Delta #phi_{min}(jet, MET)";

  selection_=baseline&&trigger && zl && isotk && jets && btag2sf && TCut("METsig>30 && METsig<=50") && drmax && higgsSB;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_2bSB_minDeltaPhi20uber",0);
  selection_=baseline&&trigger && zl && isotk && jets && btag3sf && TCut("METsig>30 && METsig<=50") && drmax && higgsSB;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_3bSB_minDeltaPhi20uber",0);
  selection_=baseline&&trigger && zl && isotk && jets && btag2sf && TCut("METsig>30 && METsig<=50") && drmax && higgsSR;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_2bSIG_minDeltaPhi20uber",0);
  //can't add 3b SIG, 4b SB, or 4b SIG until we're unblind
  selection_=baseline&&trigger && zl && isotk && jets && btag4sf && TCut("METsig>30 && METsig<=50") && drmax && higgsSB;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_4bSB_minDeltaPhi20uber",0);
//   selection_=baseline&&trigger && zl && isotk && jets && btag2 && btag3&&!btag4 && TCut("METsig>30 && METsig<=50") && drmax && higgsSR;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_3bSIG_minDeltaPhi20uber",0);
//   selection_=baseline&&trigger && zl && isotk && jets && btag2 && btag3&&btag4 && TCut("METsig>30 && METsig<=50") && drmax && higgsSR;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_4bSIG_minDeltaPhi20uber",0);

  selection_=baseline&&trigger && zl && isotk && jets && btag4sf && TCut("METsig>30 && METsig<=50") && drmax && higgsSB && mdp;
   drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_4bSB_minDeltaPhi20uber-mdpcut",0);

  //LDP plot of maxDR (no higgs mass window)
  selection_ = baseline && trigger && zl && isotk && jets && btag2sf && !mdp && metsigloose  ;
  nbins=20; low = 0; high=5;
  setLogY(false);
  var="deltaRmax_hh"; xtitle="#Delta R_{max}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bMETsig30_dRmax",0);

  //LDP plot of maxDR (no higgs mass window)
  //3b
  selection_ = baseline && trigger && zl && isotk && jets && btag3sf && !mdp && metsigloose  ;
  nbins=20; low = 0; high=5;
  setLogY(false);
  var="deltaRmax_hh"; xtitle="#Delta R_{max}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l3bMETsig30_dRmax",0);
  //4b
  selection_ = baseline && trigger && zl && isotk && jets &&btag4sf && !mdp && metsigloose  ;
  nbins=20; low = 0; high=5;
  setLogY(false);
  var="deltaRmax_hh"; xtitle="#Delta R_{max}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l4bMETsig30_dRmax",0);


  //almost an N-1 plot but don't apply DeltaRmax
  selection_ = baseline && trigger && zl && isotk && jets && btag2sf && mdp && metsigloose && higgsSR_d;
  nbins=20; low = 0; high=200;
  setLogY(false);
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{jj}> (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bN-2_higgsMass",0);

  //now make it N-1
  selection_ = baseline && trigger && zl && isotk && jets && btag2sf && mdp && metsigloose && higgsSR_d &&drmax;
  nbins=20; low = 0; high=200;
  setLogY(false);
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{jj}> (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bN-1_higgsMass",0);

  //N-2 and N-1 of Delta Higgs mass

  //almost an N-1 plot but don't apply DeltaRmax
  selection_ = baseline && trigger && zl && isotk && jets && btag2sf && mdp && metsigloose && higgsSR_av;
  nbins=20; low = 0; high=100;
  setLogY(false);
  var="abs(higgsMbb1MassDiff-higgsMbb2MassDiff)"; xtitle="|#Deltam_{jj}| (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bN-2_higgsMassDiff",0);

  //now make it N-1
  selection_ = baseline && trigger && zl && isotk && jets && btag2sf && mdp && metsigloose && higgsSR_av &&drmax;
  nbins=20; low = 0; high=100;
  setLogY(false);
  var="abs(higgsMbb1MassDiff-higgsMbb2MassDiff)"; xtitle="|#Deltam_{jj}| (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bN-1_higgsMassDiff",0);

  //N-1 plot of maxDR
  selection_ = baseline && trigger && zl && isotk && jets && btag2sf && mdp && metsigloose && higgsSR ;
  nbins=20; low = 0; high=5;
  setLogY(false);
  var="deltaRmax_hh"; xtitle="#DeltaR_{max}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bN-1_dRmax",0);

  // don't bother with b-tag SF for these diagnostic plots

  //apply full selection with 2b && !3b
  selection_=baseline&&trigger&&zl&&isotk&&jets&&btag2&&!btag3 &&higgsSR && metsigloose &&drmax &&mdp;
  setLogY(true);
  var="abs(min(abs(METphi-caloMETphi),abs(2*3.14159-METphi+caloMETphi))-3.14159)"; xtitle="| |#Delta #phi_{PF,calo}| - #pi |";
  nbins = 30; low=0; high=TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bAllCuts_RazorNoiseVar",0);

  //apply full selection (2b, !3b) except no higgs window cut
  selection_=baseline&&trigger&&zl&&isotk&&jets&&btag2&&!btag3  && metsigloose &&drmax &&mdp;
  var="abs(min(abs(METphi-caloMETphi),abs(2*3.14159-METphi+caloMETphi))-3.14159)"; xtitle="| |#Delta #phi_{PF,calo}| - #pi |";
  setLogY(true);
  nbins = 30; low=0; high=TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bNoHiggsMass_RazorNoiseVar",0);

  //apply full selection with 2b && !3b ; ldp region
  selection_=baseline&&trigger&&zl&&isotk&&jets&&btag2&&!btag3 &&higgsSR && metsigloose &&drmax &&!mdp;
  setLogY(true);
  var="abs(min(abs(METphi-caloMETphi),abs(2*3.14159-METphi+caloMETphi))-3.14159)"; xtitle="| |#Delta #phi_{PF,calo}| - #pi |";
  nbins = 30; low=0; high=TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2bLDP_2bT_RazorNoiseVar",0);

/*
  //2b && !3b, no mdp cut
  //=== begin using this selection
  selection_=baseline&&trigger&&zl&&isotk&&jets&&btag2&&!btag3 &&higgsSR && metsigloose &&drmax;
   setLogY(false);

  var="jetpt1"; xtitle=var;
  nbins = 30; low=0; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_jetpt1",0);

  var="jetpt2"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_jetpt2",0);

  var="jetpt3"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_jetpt3",0);

  var="jetpt4"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_jetpt4",0);

  var="jeteta1"; xtitle=var;
  nbins = 30; low=-3; high=3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_jeteta1",0);

  var="jeteta2"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_jeteta2",0);

  var="jeteta3"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_jeteta3",0);

  var="jeteta4"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_jeteta4",0);

  var="jetphi1"; xtitle=var;
  nbins = 30; low=-TMath::Pi(); high=TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_jetphi1",0);

  var="jetphi2"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_jetphi2",0);

  var="jetphi3"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_jetphi3",0);

  var="jetphi4"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_jetphi4",0);

  //now the variables for "all eta" for the lead jet
  var="alletajetpt1"; xtitle=var;
  nbins = 30; low=0; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_alletajetpt1",0);

  var="alletajetphi1"; xtitle=var;
  nbins = 30; low=-TMath::Pi(); high=TMath::Pi();
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_alletajetphi1",0);

  var="alletajeteta1"; xtitle=var;
  nbins = 30; low=-5; high=5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_alletajeteta1",0);

  var="alletajetneutralhadronfrac1"; xtitle=var;
  nbins = 30; low=0; high=1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_alletajetneutralhadronfrac1",0);

  var="alletajetneutralemfrac1"; xtitle=var;
  nbins = 30; low=0; high=1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_alletajetneutralemfrac1",0);

  var="alletajetneutralphotonfrac1"; xtitle=var;
  nbins = 30; low=0; high=1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_alletajetneutralphotonfrac1",0);

  //some b-jet variables
  var="bjetchargedhadronmult1"; xtitle=var;
  nbins = 30; low=0; high=100;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_bjetchargedhadronmult1",0);

  var="bjetchargedhadronmult2"; xtitle=var;
  nbins = 30; low=0; high=100;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_0l2b_bjetchargedhadronmult2",0);
*/

}

void higgs_dataMC_debug() {

  initHiggsSamples69(false,"ttbar");

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";

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

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
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
  initHiggsSamples69(false,"ttbar wjets qcd znunu singlet"); 

  //use top pT weights
  setSampleWeightFactor("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71","topPtWeight");
  setSampleWeightFactor("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71","topPtWeight");
  setSampleWeightFactor("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71",   "topPtWeight");

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  setOutputDirectory("plots_Control_SL");

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  //cannot use this with MET PD only
  //  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoPFcands10_010==0";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0";//remove veto on 2nd iso tk

  TCut njets4="njets20==4"; //switch back to 20 GeV jets
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  //TCut mdp = "minDeltaPhi20>0.5 || (minDeltaPhi20>0.3&&METsig>50)";
  TCut mdp = "minDeltaPhi20_eta5_noIdAll_nobeta>0.5 || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3&&METsig>50)";
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

void higgs_dataMC_weightFactorTest_SL(TString options) {
  initHiggsSamples69(true,"ttbar znunu bjets wjets ttv vv singlet"); 

  //top pt weights
  TString weightfactor="";
  if (options.Contains("toppt")&&options.Contains("btagsf") ) weightfactor = "topPtWeight*BTagWeightLMT";
  else  if ( options.Contains("btagsf") ) weightfactor = "BTagWeightLMT";
  else  if ( options.Contains("toppt")  )  weightfactor = "topPtWeight";

  //set weights
  setSampleWeightFactor("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883ra2b_v71s",weightfactor);
  setSampleWeightFactor("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884ra2b_v71s",weightfactor);
  //this sample is baically irrelevant
  //  setSampleWeightFactor("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880ra2b_v71s",weightfactor);

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  setOutputDirectory("plots_Control_SL_weights");


  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoPFcands10_010==0";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0";
  TCut sle = "nMuons==0&&nElectrons==1 &&nTausLoose==0";
  TCut slm = "nMuons==1&&nElectrons==0 &&nTausLoose==0";


  TCut njets4="njets20==4"; //switch back to 20 GeV jets
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut jets = (njets4||njets5) && jet2;
  TCut mdp = "minDeltaPhi20>0.5 || (minDeltaPhi20>0.3&&METsig>50)";
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

  TCut metsig[4];
  metsig[0]="METsig>=30&&METsig<50";
  metsig[1]="METsig>=50&&METsig<100";
  metsig[2]="METsig>=100&&METsig<150";
  metsig[3]="METsig>=150";

  TCut metsigloose = "METsig>=30";

  TCut drmax = "deltaRmax_hh<2.2";

  setStackMode(true,false,false); //stack,norm,label override
  setLogY(false);

  map<int,TCut> nbcut;
  nbcut[2] =  btag2 && !btag3;
  nbcut[3] =  btag2 && btag3 && !btag4;
  nbcut[4] =  btag2 && btag3 && btag4;


  //plot METsig distribution in b-tag bins (SB only)

  for (int ib = 2;ib<=4; ib++) {
    TString filename;

    selection_ = baseline && trigger && sl && jets && mdp && higgsSB && drmax && metsigloose && nbcut[ib]; //nb SB SL

    float varbinsmetsig[6] = {0,30,50,100,150,500};
    nbins=5; low=0; high=500;
    var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
    filename.Form("SL_SB%db_METsig_%s",ib,jmt::fortranize(weightfactor).Data());
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,varbinsmetsig);
    
    //may want to try log scale here too

    //plot lepton pt spectra -- for all METsig bins
    nbins=10; low=0; high=150;
    var="muonpt1"; xtitle=var;
    selection_=baseline && trigger &&slm &&jets&&mdp&&higgsSB && drmax && metsigloose && nbcut[ib]; //2b SB SLmuon
    filename.Form("SL_SB%db_muPt_%s",ib,jmt::fortranize(weightfactor).Data());
    drawPlots(var,nbins,low,high,xtitle,"Events", filename);
    
    nbins=10; low=0; high=150;
    var="eleet1"; xtitle=var;
    selection_=baseline && trigger &&sle &&jets&&mdp&&higgsSB && drmax && metsigloose && nbcut[ib]; //2b SB SLelectron
    filename.Form("SL_SB%db_elEt_%s",ib,jmt::fortranize(weightfactor).Data());
    drawPlots(var,nbins,low,high,xtitle,"Events", filename);

  }

}


void higgs_dataMC_control_SL_btagSF_test() {

  //goal do data / MC comparisons -- use btag SF: necessitates strictly plotting
  //  initHiggsSamples69(true,"ttbar znunu bjets wjets ttv vv singlet"); //use all samples except signal
  initHiggsSamples69(false,"wjets ttbar singlet"); //use only ttbar; also these samples are unskimmed

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0.4; ratioMax = 1.6;
  dodata_=true;

  setOutputDirectory("plots_Control_SL_btagSF");

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoPFcands10_010==0";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0";
  //we used to veto "lepton cleaned" isolated tracks from the SL sample (but not in RA2b)
  //i didn't add this lepton-cleaned equivalent for the new variable
  //could use the old one, or just forget it. I think I'll just keep things as simple as possible
  // &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets20==4"; //switch back to 20 GeV jets
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  //  TCut mdp = "minDeltaPhi20>0.5 || (minDeltaPhi20>0.3&&METsig>50)";
  TCut mdp = "minDeltaPhi20_eta5_noIdAll_nobeta>0.5 || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3 && METsig>50)";
  TCut btag2="nbtag2_rawMC==1";
  TCut btag3="nbtag3_rawMC==1";
  TCut btag4="nbtag4_rawMC==1";

  TCut btag2sf="nbtag2_nomSF==1";
  TCut btag2sfu="nbtag2_SFp1sig==1";
  TCut btag2sfd="nbtag2_SFm1sig==1";

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

  setStackMode(true,false,false); //stack,norm,label override
  //only needed in the case where we had reset the samples above


  nbins=20; low=0; high=200;
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{jj}>";

  //truly raw MC
  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btag2;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_amh_rawMC",0);

  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btag2sf;
  drawBTagErrors_=true;
  drawTopPtErrors_=false;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_amh_BtagSFOnlyWithErrors",0);

  //use "official" top pT weights
  setSampleWeightFactor("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71","topPtWeightOfficial");
  setSampleWeightFactor("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71","topPtWeightOfficial");
  setSampleWeightFactor("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71",   "topPtWeightOfficial");

  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btag2;
  drawBTagErrors_=false;
  drawTopPtErrors_=true;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_amh_topPtOnlyWithErrors",0);

  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btag2sf;
  drawBTagErrors_=true;
  drawTopPtErrors_=false;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_amh_topPtAndBtagSF_BTagErrors",0);

  drawBTagErrors_=false;
  drawTopPtErrors_=true;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_amh_topPtAndBtagSF_TopPtErrors",0);

  drawBTagErrors_=true;
  drawTopPtErrors_=true;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_amh_topPtAndBtagSF_BtagAndTopPtErrors",0);

}

void higgs_dataMC_control_SL() {

  //goal do data / MC comparisons; use skim to save time. restricts what plots can be made.
  //  initHiggsSamples69(true,"ttbar znunu bjets wjets ttv vv singlet"); //use all samples except signal
  //for now we have to use no skim
  initHiggsSamples69(false,"ttbar znunu bjets wjets singlet"); //use all samples except signal

  //use top pT weights
  setSampleWeightFactor("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883ra2b_v71s","topPtWeightOfficial");
  setSampleWeightFactor("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884ra2b_v71s","topPtWeightOfficial");
  setSampleWeightFactor("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880ra2b_v71s",   "topPtWeightOfficial");

  //use top pT weights
  setSampleWeightFactor("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71","topPtWeightOfficial");
  setSampleWeightFactor("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71","topPtWeightOfficial");
  setSampleWeightFactor("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71",   "topPtWeightOfficial");


  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  setOutputDirectory("plots_Control_SL");

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoPFcands10_010==0";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0";
  //we used to veto "lepton cleaned" isolated tracks from the SL sample (but not in RA2b)
  //i didn't add this lepton-cleaned equivalent for the new variable
  //could use the old one, or just forget it. I think I'll just keep things as simple as possible
  // &&nIsoTracks15_005_03_lepcleaned==0";

  TCut njets4="njets20==4"; //switch back to 20 GeV jets
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  //  TCut mdp = "minDeltaPhi20>0.5 || (minDeltaPhi20>0.3&&METsig>50)";
  TCut mdp = "minDeltaPhi20_eta5_noIdAll_nobeta>0.5 || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3 && METsig>50)";
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

  setStackMode(true,false,false); //stack,norm,label override
  setLogY(false);

  TCut btagge2="nbtag2_nomSF==1||nbtag3_nomSF==1||nbtag4_nomSF==1";
  TCut btagge3="nbtag3_nomSF==1||nbtag4_nomSF==1";
  TCut btag4sf="nbtag4_nomSF==1";
  //2b pre-selection ; basically the minimum required by the trigger, plus SL
  //  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btag2;
  //switch to using btag SF
  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btagge2;
  addVerticalLine(30);
  nbins=30; low=0; high=300;
  setLogY(false);
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig",0);
  setLogY(true); setPlotMinimum(0.5);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig",0);
  resetPlotMinimum();
  resetVerticalLine();
  //MET, for completeness
  nbins=30; low=0; high=500;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_MET",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_MET",0);

  //check data/MC agreement in bins of time
  //one trick here is that the PU reweighting might really get this wrong -- these plots might be deceptive in that way
  nbins=30; low=0; high=300;
  setLogY(false);
  addVerticalLine(30);
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";

  double fullLumi = lumiScale_;

  // to check for run dependence
  setDatasetToDraw("META"); lumiScale_ = 807.1;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig_2012A",0);

  setDatasetToDraw("METB"); lumiScale_ = 4421;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig_2012B",0);

  setDatasetToDraw("METC"); lumiScale_ = 495+6402;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig_2012C",0);

  setDatasetToDraw("METD"); lumiScale_ = 5956+1317;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig_2012D",0);

  //reset to full dataset
  setDatasetToDraw("MET"); 
  lumiScale_=fullLumi;

  resetVerticalLine();
  //nPV
  nbins=20; low=0; high=40;
  setLogY(false);
  var="nGoodPV"; xtitle="n Good PV";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_nGoodPV",0);


  // === now plot in bins of PV; this might be better than in bins of time

  drawFilenameOnPlot_=true;

  TCut pvcut[3];
  pvcut[0] = "nGoodPV <=12";
  pvcut[1] = "nGoodPV >12 && nGoodPV<=21";
  pvcut[2] = "nGoodPV >21";

  //  TH1D* ratiosByPV[4];
  addVerticalLine(30);
  for (int ipv=0; ipv<3;ipv++) {
    selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btagge2 && pvcut[ipv];
    TString pvdesc=jmt::fortranize(pvcut[ipv].GetTitle());

    nbins=30; low=0; high=300;
    var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
    setLogY(true);
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselection_METsig_"+pvdesc,0);

  }

  resetVerticalLine();

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

  /// === done with bins of PV
  */
  drawFilenameOnPlot_=false;

  //minDeltaPhi distribution (no MET or METsig cut)
  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2 && btagge2;
  nbins=30; low=0; high=3;
  setLogY(false);
  var="minDeltaPhi20_eta5_noIdAll_nobeta"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselectionNoMdpNoMETsig_minDeltaPhi20",0);

  //add the standard METsig cut
  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2 && btagge2 &&metsigveryloose;
  nbins=30; low=0; high=3;
  setLogY(false);
  var="minDeltaPhi20_eta5_noIdAll_nobeta"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_preselectionNoMdp_minDeltaPhi20",0);

  //CSV3
  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btagge2 &&metsigveryloose;
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
  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp &&btagge3 &&metsigveryloose;

  //CSV4
  nbins=10; low=0; high=1;
  setLogY(false);
  var="CSVbest4"; xtitle="4th CSV value";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_3bPreselection_CSVbest4",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_3bPreselection_CSVbest4",0);

  //now apply full b-tag selection
  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btag4sf;

  //old MET sig
  addVerticalLine(30);
  nbins=15; low=0; high=300;
  setLogY(false);
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bPreselection_METsig",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bPreselection_METsig",0);
  resetVerticalLine();
  //MET
  nbins=15; low=0; high=400;
  setLogY(false);
  var="MET"; xtitle="MET (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bPreselection_MET",0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bPreselection_MET",0);

 //now apply full b-tag selection and METpreselection
  selection_ = baseline && trigger && sl && (njets4||njets5) &&jet2&&mdp && btag4sf &&metsigveryloose;
  nbins=10; low=0; high=100;
  setLogY(false);
  var="abs(higgsMbb1MassDiff-higgsMbb2MassDiff)"; xtitle="|#Deltam_{jj}| (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bMETsig30_DeltaMbb",0,"GeV");

  //higgs mass
  nbins=20; low=0; high=200;
  setLogY(false);
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{jj}> (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_4bMETsig30_AvMbb",0,"GeV");

  //now add the higgs SR
  selection_ = baseline && trigger && sl && (njets4||njets5)&&jet2&&mdp && btag4sf &&metsigveryloose &&higgsSR;
  nbins=20; low=0; high=5;
  setLogY(false);
  var="deltaRmax_hh"; xtitle="#Delta R_{max}(b,b)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_hhSR-METsig30_DRmax",0);

  //true N - 1 plots -- for the AN
  selection_ = baseline && trigger && sl && (njets4||njets5)&&jet2&&mdp && btag4sf &&metsigveryloose &&higgsSR_d &&drmax;
  //higgs mass
  nbins=20; low=0; high=200;
  setLogY(false);
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{jj}>";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_Nm1METsig30_AvMbb",0);

  selection_ = baseline && trigger && sl && (njets4||njets5)&&jet2&&mdp && btag4sf &&metsigveryloose &&higgsSR_av &&drmax;
  nbins=10; low=0; high=100;
  setLogY(false);
  var="abs(higgsMbb1MassDiff-higgsMbb2MassDiff)"; xtitle="|#Deltam_{jj}|";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_Nm1METsig30_DeltaMbb",0);

  //
  addVerticalLine(30);
  selection_ = baseline && trigger && sl && (njets4||njets5)&&jet2&&mdp && btag4sf  &&higgsSR &&drmax;
  nbins=30; low=0; high=300;
  setLogY(true);
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_Nm1_METsig",0);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_dataMC_SL_Nm1_METsig",0);
  resetVerticalLine();


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

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
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

  TCut drmax = "deltaRmax_hh<2.2";
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

  selection_=baseline&&trigger&&sl&&(njets4||njets5)&& btag2&&!btag3 &&higgsSR&&drmax;
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

  //this code also generates the plots used for the ABCD cross-check in SL data
  //there is a "sl" or "zl" setting below. The SL setting is used for the data-driven cross-check

  initHiggsSamples69(true,"znunu ttbar");
  //  initHiggsSamples(true,"ttbarjoint");

  const  TString sample = "sl"; //or zl
  assert(sample=="sl" || sample=="zl");

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  //  TCut triggerJets = "passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoPFcands10_010==0";
  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut mdp = "minDeltaPhi20_eta5_noIdAll_nobeta>0.5 || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3 && METsig>50)";
  TCut jet2="jetpt2>50";
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

  TCut drmax = "deltaRmax_hh<2.2";
  //  TCut drmin = "deltaRmin_hh<1.9";


  //Owen's METsig bins
  float binning1[]={0,10,20,30,50,100,150,200};
  float binning2[]={0,20,30,50,100,200};
  float * binning=0;
  if (sample=="zl") binning = binning1;
  else if (sample=="sl") binning=binning2;
  nbins=7; low=0; high=200;
  if (sample=="sl") nbins=5; //use binning2
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
    selection_ = baseline&&trigger&&leptoncut&&njets&&drmax&&mdp&&jet2 && higgsSR && btag2 && !btag3; // was notbtag3
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SRSBshapes",binning);
    TH1D* h_SR_2b = (TH1D*) totalsm->Clone("h_SR_2b");
    TH1D* h_SR_2b_data=0; if (dodata_) h_SR_2b_data = (TH1D*) hdata->Clone("h_SR_2b_data");
    
    //SR cuts ; 3 btags
    selection_ = baseline&&trigger&&leptoncut&&njets&&drmax&&mdp&&jet2 && higgsSR && btag2&&btag3 &&!btag4;
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SRSBshapes",binning);
    TH1D* h_SR_3b = (TH1D*) totalsm->Clone("h_SR_3b");
    TH1D* h_SR_3b_data=0; if (dodata_) h_SR_3b_data = (TH1D*) hdata->Clone("h_SR_3b_data");
    
    //SR cuts ; 4 btags
    selection_ = baseline&&trigger&&leptoncut&&njets&&drmax&&mdp&&jet2 && higgsSR && btag2&&btag3&&btag4;
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SRSBshapes",binning);
    TH1D* h_SR_4b = (TH1D*) totalsm->Clone("h_SR_4b");
    TH1D* h_SR_4b_data=0; if (dodata_) h_SR_4b_data = (TH1D*) hdata->Clone("h_SR_4b_data");

    // -- now SB --
    
    //SB cuts ; 2 btags
    selection_ = baseline&&trigger&&leptoncut&&njets&&drmax&&mdp&&jet2 && higgsSB && btag2 &&!btag3;
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SRSBshapes",binning);
    TH1D* h_SB_2b = (TH1D*) totalsm->Clone("h_SB_2b");
    TH1D* h_SB_2b_data=0; if (dodata_) h_SB_2b_data = (TH1D*) hdata->Clone("h_SB_2b_data");
    
    //SB cuts ; 3 btags
    selection_ = baseline&&trigger&&leptoncut&&njets&&drmax&&mdp&&jet2 && higgsSB && btag2&&btag3 &&!btag4;
    drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_SRSBshapes",binning);
    TH1D* h_SB_3b = (TH1D*) totalsm->Clone("h_SB_3b");
    TH1D* h_SB_3b_data=0; if (dodata_) h_SB_3b_data = (TH1D*) hdata->Clone("h_SB_3b_data");

    //SB cuts ; 4 btags
    selection_ = baseline&&trigger&&leptoncut&&njets&&drmax&&mdp&&jet2 && higgsSB && btag2&&btag3&&btag4;
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

    if (dodata_) {//just want a simple plot of SR/SB ratio in bins of btag

      ratio_2b_data->SetLineColor(kBlue);
      ratio_3b_data->SetLineColor(kRed);
      ratio_4b_data->SetLineColor(kMagenta);
      
      ratio_2b_data->SetMarkerColor(kBlue);
      ratio_3b_data->SetMarkerColor(kRed);
      ratio_4b_data->SetMarkerColor(kMagenta);

      ratio_2b_data->SetMarkerStyle(21);
      ratio_3b_data->SetMarkerStyle(24);
      ratio_4b_data->SetMarkerStyle(25);
      
      ratio_2b_data->SetYTitle("Signal region / sideband");
      ratio_2b_data->SetXTitle("MET significance");
      ratio_3b_data->SetYTitle("Signal region / sideband");
      ratio_3b_data->SetXTitle("MET significance");
      ratio_4b_data->SetYTitle("Signal region / sideband");
      ratio_4b_data->SetXTitle("MET significance");

      renewCanvas();
      renewLegend();
      leg->AddEntry(ratio_2b_data,"2b");
      leg->AddEntry(ratio_3b_data,"3b");
      leg->AddEntry(ratio_4b_data,"4b");
      ratio_3b_data->Draw("hist e");
      ratio_4b_data->Draw("hist e same");
      ratio_2b_data->Draw("hist e same");
      ratio_3b_data->SetMinimum(0);
      ratio_3b_data->SetMaximum(0.5);
      leg->Draw();
      thecanvas->SaveAs(TString("higgs_SL_closureInData_")+njetsdesc[ijets]+".eps");
      thecanvas->SaveAs(TString("higgs_SL_closureInData_")+njetsdesc[ijets]+".pdf");
      thecanvas->SaveAs(TString("higgs_SL_closureInData_")+njetsdesc[ijets]+".png");

      double chi2_23=0;
      double chi2_24=0;
      for (int ibin=3; ibin<=5 ; ++ibin) {
	chi2_23 += pow(ratio_2b_data->GetBinContent(ibin)-ratio_3b_data->GetBinContent(ibin),2) / pow(jmt::addInQuad(ratio_2b_data->GetBinError(ibin),ratio_3b_data->GetBinError(ibin)),2);
	chi2_24 += pow(ratio_2b_data->GetBinContent(ibin)-ratio_4b_data->GetBinContent(ibin),2) / pow(jmt::addInQuad(ratio_2b_data->GetBinError(ibin),ratio_4b_data->GetBinError(ibin)),2);
      }
      cout<<"2b3b chi2 (/dof) = "<<chi2_23<<" "<<chi2_23 / 3.0<<endl;
      cout<<"2b4b chi2 (/dof) = "<<chi2_24<<" "<<chi2_24 / 3.0<<endl;

      //ok, in fact this is too hard to read. so do one at a time:
      renewCanvas();
      renewLegend();
      leg->AddEntry(ratio_2b_data,"2b");
      leg->AddEntry(ratio_3b_data,"3b");
      ratio_3b_data->Draw("hist e");
      ratio_2b_data->Draw("hist e same");
      ratio_3b_data->SetMinimum(0);
      ratio_3b_data->SetMaximum(0.5);
      leg->Draw();
      thecanvas->SaveAs(TString("higgs_SL_closureInData_2b3b_")+njetsdesc[ijets]+".eps");
      thecanvas->SaveAs(TString("higgs_SL_closureInData_2b3b_")+njetsdesc[ijets]+".pdf");
      thecanvas->SaveAs(TString("higgs_SL_closureInData_2b3b_")+njetsdesc[ijets]+".png");

      renewCanvas();
      renewLegend();
      leg->AddEntry(ratio_2b_data,"2b");
      leg->AddEntry(ratio_4b_data,"4b");
      ratio_4b_data->Draw("hist e");
      ratio_2b_data->Draw("hist e same");
      ratio_4b_data->SetMinimum(0);
      ratio_4b_data->SetMaximum(0.5);
      leg->Draw();
      thecanvas->SaveAs(TString("higgs_SL_closureInData_2b4b_")+njetsdesc[ijets]+".eps");
      thecanvas->SaveAs(TString("higgs_SL_closureInData_2b4b_")+njetsdesc[ijets]+".pdf");
      thecanvas->SaveAs(TString("higgs_SL_closureInData_2b4b_")+njetsdesc[ijets]+".png");

    }

    return; //screw the kappa factor stuff. simple plots are better

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

  initHiggsSamples69(false,"ttbar hh250 hh400"); //skim has triggers -- can't use it!

  int nbins;
  float low,high;
  TString var,xtitle;


  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0"; 
  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
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
  initHiggsSamples69(false,"znunu ttbar hh150 hh200 hh250 hh400 wbb ttv"); //skim has triggers -- can't use it!

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
  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
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
  useTrigEff_=false ; //doesn't seem needed here

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0 && MET/caloMET<2 && maxTOBTECjetDeltaMult<40";
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1||passMC_PFMET150==1";
  TCut zl = "nMuons==0&&nElectrons==0";
  TCut tauveto="nTausLoose==0";
  TCut isotk1="nIsoPFcands10_010==0";

  TCut njets="njets20==4 || njets20==5";
  TCut jet2="jetpt2>50";
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut mdp = "minDeltaPhi20_eta5_noIdAll_nobeta>0.5 || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3 && METsig>50)";

  TCut higgsmass = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsdiff = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";

  TCut drmax = "deltaRmax_hh<2.2";

  TCut metsig = "METsig>30";

  setStackMode(true,false,false); //stack, norm, labels

  stackSignal_=false; //let's unstack signal

  clearSamples();
  addSample("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883ra2b_v71s",kBlue-10,"t#bar{t} (2l)");
  addSample("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884ra2b_v71s:ttbarDecayCode==7",kMagenta-7,"t#bar{t} (#tau #rightarrow h)");
  addSample("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884ra2b_v71s:ttbarDecayCode==3||ttbarDecayCode==4",kRed-7,"t#bar{t} (W #rightarrow e/#mu)");
  addSample("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884ra2b_v71s:ttbarDecayCode==5||ttbarDecayCode==6",kOrange+10,"t#bar{t} (#tau #rightarrow e/#mu)");

  //higgs mass
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2&&mdp && btag2 && btag3 && btag4  && higgsdiff && drmax && metsig;
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{jj}>";
  nbins=20; low= 0; high=200;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_higgsmass_splitTtbar",0,"GeV");

  //all cuts except metsig
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2&&mdp && btag2 && btag3 && btag4  && higgsdiff&&higgsmass && drmax ;
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
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
  addSample("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883ra2b_v71s",kBlue-10,"t#bar{t} (2l)");
  addSample("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884ra2b_v71s:hadWcode==11",kBlue+2,"t#bar{t} (W #rightarrow ud)");
  addSample("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884ra2b_v71s:hadWcode==1100",kGreen+2,"t#bar{t} (W #rightarrow cs)");

  var="CSVbest4"; xtitle="4th best CSV value";
  nbins=20; low= 0; high=1;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_splitTtbarByHadW_CSV4",0);

  var="nPUjets20"; xtitle="Number of PU jets (pT>20 GeV)";
  nbins=3; low= 0; high=3;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_splitTtbarByHadW_nPU20",0);


}


void higgs_Nminus1(bool plotdata=false,TString options="4b") {
  if (!plotdata)
    initHiggsSamples69(true,"bjets wjets ttv vv singlet ttbar znunu hhmg200 hhmg400");
  else //can't use skim right now, nor can i use ttv vv
    initHiggsSamples69(false,"bjets wjets singlet ttbar znunu hhmg200 hhmg400");

  //use "official" top pT weights
  setSampleWeightFactor("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71","topPtWeightOfficial");
  setSampleWeightFactor("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71","topPtWeightOfficial");
  setSampleWeightFactor("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71",   "topPtWeightOfficial");

  useTrigEff_=plotdata; //use trig eff correction only if we're plotting data
  usePUweight_=true; 

  if (plotdata)  setOutputDirectory("plots_higgs_Nminus1_unblind");
  else   setOutputDirectory("plots_higgs_Nminus1");

  int nbins;
  float low,high;
  TString var,xtitle;

  savePlots_=true;
  setQuiet(false);

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=plotdata;

  //20 June -- update
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0 && MET/caloMET<2 && maxTOBTECjetDeltaMult<40";
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1||passMC_PFMET150==1||passMC_DiCentralPFJet30_PFMHT80==1";
  TCut zl = "nMuons==0&&nElectrons==0";
  TCut tauveto="nTausLoose==0";
  //  TCut isotk1="nIsoTracks15_005_03==0";
  TCut isotk1="nIsoPFcands10_010==0";

  TCut njets="njets20==4 || njets20==5"; //back to 20 GeV
  TCut jet2="jetpt2>50";
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut btag2sf="nbtag2_nomSF==1";
  TCut btag3sf="nbtag3_nomSF==1";
  TCut btag4sf="nbtag4_nomSF==1";


  //  TCut mdp = "minDeltaPhi20>0.5 || (minDeltaPhi20>0.3&&METsig>50)";
  TCut mdp = "minDeltaPhi20_eta5_noIdAll_nobeta>0.5 || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3 && METsig>50)";

  TCut higgsmass = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsdiff = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";

  TCut drmax = "deltaRmax_hh<2.2";

  TCut metsig = "METsig>30";

  /*
  TCut btagcut = btag2;
  if ( options.Contains("3b")) btagcut = btag2&&btag3&&!btag4;
  else if (options.Contains("4b")) btagcut = btag2&&btag3&&btag4;
  */
  TCut btagcut = btag2sf;
  if ( options.Contains("3b")) btagcut = btag3sf;
  else if (options.Contains("4b")) btagcut = btag4sf;

  setStackMode(true,false,false); //stack, norm, labels

  stackSignal_=false; //let's unstack signal

  TString nm1label="higgs_Nm1";
  if (options.Contains("3b")) nm1label+= "-3b";
  nm1label+="_";

  //full selection 
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 &&btagcut &&mdp && higgsmass && higgsdiff && drmax && metsig;
  var="bjetpt1"; xtitle="lead b jet pT";
  nbins=20; low= 0; high=400;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"allcuts_bjetpt1",0,"GeV");

  //still full selection
  var="jetpt1"; xtitle="lead jet pT";
  nbins=30; low= 0; high=300; //need 10 GeV Bins
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"allcuts_jetpt1",0,"GeV");

  var="jetpt2"; xtitle="2nd jet pT";
  nbins=30; low= 0; high=200; //need 10 GeV Bins
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"allcuts_jetpt2",0,"GeV");

  var="jetpt3"; xtitle="3rd jet pT";
  nbins=20; low= 0; high=200; //need 10 GeV Bins
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"allcuts_jetpt3",0,"GeV");

  var="jetpt4"; xtitle="4th jet pT";
  nbins=20; low= 0; high=200; //need 10 GeV Bins
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"allcuts_jetpt4",0,"GeV");

  //not too interesting
  var="higgsWCandMass"; xtitle="W Candidate mass (GeV)";
  nbins=40; low= 0; high=200;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"allcuts_WcandMass",0,"GeV");


  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 &&btagcut&&mdp && higgsmass && higgsdiff && drmax && metsig;
  var="njets20_5p0-njets20"; xtitle="number of forward jets (pT>20)";
  nbins=4; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"allcuts_njets20forward",0);

  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 &&btagcut&&mdp && higgsmass && higgsdiff && drmax && metsig;
  var="deltaPhi_hh"; xtitle="#Delta #phi (h,h)";
  nbins=10; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"allcuts_hhDeltaPhi",0);

  //full selection except mdp
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 &&btagcut && higgsmass && higgsdiff && drmax && metsig;
/* not very interesting anymore
  var="deltaPhi1"; xtitle="#Delta #phi (jet1,MET)";
  nbins=10; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"allcuts_deltaPhi1",0);
  var="deltaPhi2"; xtitle="#Delta #phi (jet2,MET)";
  nbins=10; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"allcuts_deltaPhi2",0);
  var="deltaPhi3"; xtitle="#Delta #phi (jet3,MET)";
  nbins=10; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"allcuts_deltaPhi3",0);
*/

//these could be interesting -- shows how signal changes with the different versions
  var="minDeltaPhi20"; xtitle="min #Delta #phi (jet1..3,MET)";
  nbins=10; low= 0; high=3;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"allcuts_deltaPhiMin",0);

  var="minDeltaPhi20_eta5_noIdAll"; xtitle="min #Delta #phi (all jets, MET)";
  nbins=10; low= 0; high=3;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"allcuts_deltaPhiMinAll",0);

  var="minDeltaPhi20_eta5_noIdAll_nobeta"; xtitle="min #Delta #phi (all jets, MET)";
  nbins=10; low= 0; high=3;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"allcuts_deltaPhiMinUber",0);

/* no longer interesting
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
*/

  //n leptons
  selection_=baseline&&triggers &&  tauveto && isotk1 && njets&&jet2 &&btagcut &&mdp && higgsmass && higgsdiff && drmax && metsig;
  var="nMuons+nElectrons"; xtitle="e + #mu";
  nbins=4; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"ePlusMu",0);

  //POG taus
  selection_=baseline&&triggers && zl&&  isotk1 && njets&&jet2 && btagcut &&mdp && higgsmass && higgsdiff && drmax && metsig;
  var="nTausLoose"; xtitle="Loose taus";
  nbins=3; low= 0; high=3;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"looseTaus",0);

  //iso tracks
  selection_=baseline&&triggers && zl&& tauveto && njets&&jet2 && btagcut &&mdp && higgsmass && higgsdiff && drmax && metsig;
  var="nIsoPFcands10_010"; xtitle="10 GeV iso tracks (PF)";
  nbins=4; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"isoPFtk10",0);

  //jet multiplicity
  selection_=baseline&&triggers && zl&& tauveto && isotk1&&jet2 && btagcut &&mdp && higgsmass && higgsdiff && drmax  && metsig;
  var="njets20"; xtitle="jet multiplicity (20 GeV)";
  nbins=9; low= 0; high=9;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"njets20",0);

  //3rd b tag -- override btagcut (and thus do not use b-tag SFs)
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btag2&&mdp  && higgsmass && higgsdiff && drmax && metsig;
  var="CSVbest3"; xtitle="3rd CSV value";
  nbins=20; low= 0; high=1;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "higgs_Nm1_CSVbest3",0);

  //4th btag -- override btagcut
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
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btagcut &&mdp && higgsdiff && drmax && metsig;
  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{jj}> (GeV)";
  nbins=20; low= 0; high=200;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"higgsmass",0,"GeV");
  //try with fewer bins too
  nbins=10; low= 0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"higgsmass_coarse",0,"GeV");

  //higgs mass difference
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btagcut&&mdp && higgsmass  && drmax && metsig;
  var="abs(higgsMbb1MassDiff-higgsMbb2MassDiff)"; xtitle="|#Deltam_{jj}| (GeV)";
  nbins=8; low= 0; high=80;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"higgsmassdiff",0,"GeV");

  //drmax
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btagcut &&mdp && higgsmass && higgsdiff && metsig;
  var="deltaRmax_hh"; xtitle="#DeltaR_{max}";
  nbins=30; low= 0; high=6;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"maxDR",0);

  //2d plot
  //draw2d(var,15,low,high,"deltaRmin_hh",15,0,4,xtitle,"min #Delta R",nm1label+"maxDRminDR_hh200",0,0,"SMS-HbbHbb_mHiggsino-200_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1807_v69");
  //draw2d(var,15,low,high,"deltaRmin_hh",15,0,4,xtitle,"min #Delta R",nm1label+"maxDRminDR_hh350",0,0,"SMS-HbbHbb_mHiggsino-350_mLSP-1_8TeV-Pythia6Z_jgsmith_UCSB1811_v69");

  //linear combination
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btagcut &&mdp && higgsmass && higgsdiff && metsig;
  var="3.318-0.139*deltaRmax_hh-0.622*deltaRmin_hh"; xtitle="#DeltaR Fisher";
  nbins=30; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"DRF",0);


  //drmin
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btagcut &&mdp && higgsmass && higgsdiff  && drmax && metsig;
  var="deltaRmin_hh"; xtitle="min #Delta R";
  nbins=40; low= 0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"minDR",0);

  //full selection except METsig
  selection_=baseline&&triggers && zl&& tauveto && isotk1&& njets&&jet2 && btagcut &&mdp&& higgsmass && higgsdiff && drmax;
  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  nbins=20; low= 0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"METsig",0);

  //just a different number of bins
  nbins=12; low= 0; high=180;
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"METsig_coarse",0);


  //full selection except METsig
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && btagcut &&mdp && higgsmass && higgsdiff && drmax;
  var="MET"; xtitle="MET";
  nbins=20; low= 0; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", nm1label+"MET",0);

}



void higgs_printCutflowTable(const bool latexMode=false,TString region="4bSR", const  bool unweighted=false, const bool makePlot=false,const bool doCutEff=false) {

/*
  if (unweighted)  initHiggsSamples69(true,"ttbar wbb"); //reduced list of samples
  else initHiggsSamples69(true);
*/
  initHiggsSamples69(true,"hhmg250 hhfull250");
  //   initHiggsSamples69(true,"hhmg400 hhfull400");

  if (!(region.Contains("4bSR")||region.Contains("2bSR")||region.Contains("4bSB")||region.Contains("2bSB")||region.Contains("special"))) {cout<<"check region!"<<endl;return;}

  int nbins;
  float low,high;
  TString var,xtitle;

  if (unweighted) {
    nameOfEventWeight_="(1)";
    usePUweight_=false;
    useTrigEff_=false;
    lumiScale_=1;
    resetSampleWeightFactors();
    resetSampleScaleFactors();
  }

  if (region.Contains("special")) {
    usePUweight_=false;
    useTrigEff_=false;
  setSampleWeightFactor("HbbHbb_FullSim_400_v71","0 + (nGoodPV==1)*0.000000 + (nGoodPV==2)*0.000000 + (nGoodPV==3)*0.806114 + (nGoodPV==4)*0.550517 + (nGoodPV==5)*0.784298 + (nGoodPV==6)*0.850867 + (nGoodPV==7)*0.765688 + (nGoodPV==8)*0.724206 + (nGoodPV==9)*0.807895 + (nGoodPV==10)*0.839914 + (nGoodPV==11)*1.055339 + (nGoodPV==12)*1.135399 + (nGoodPV==13)*1.327954 + (nGoodPV==14)*1.712817 + (nGoodPV==15)*1.617481 + (nGoodPV==16)*1.299388 + (nGoodPV==17)*0.995578 + (nGoodPV==18)*0.771311 + (nGoodPV==19)*0.800530 + (nGoodPV==20)*0.798977 + (nGoodPV==21)*0.899122 + (nGoodPV==22)*0.892036 + (nGoodPV==23)*0.901595 + (nGoodPV==24)*0.996259 + (nGoodPV==25)*0.928075 + (nGoodPV==26)*0.766434 + (nGoodPV==27)*0.988833 + (nGoodPV==28)*0.599485 + (nGoodPV==29)*0.983066 + (nGoodPV==30)*0.941055 + (nGoodPV==31)*1.057779 + (nGoodPV==32)*0.885022 + (nGoodPV==33)*1.278620 + (nGoodPV==34)*3.954734 + (nGoodPV==35)*18.560288 + (nGoodPV==36)*0.000000 + (nGoodPV==37)*0.000000 + (nGoodPV==38)*0.000000 + (nGoodPV==39)*0.000000 + (nGoodPV==40)*0.000000 + (nGoodPV==41)*0.000000 + (nGoodPV==42)*0.000000 + (nGoodPV==43)*0.000000 + (nGoodPV==44)*0.000000 + (nGoodPV==45)*0.000000 + (nGoodPV==46)*0.000000 + (nGoodPV==47)*0.000000 + (nGoodPV==48)*0.000000 + (nGoodPV==49)*0.000000 + (nGoodPV==50)*0.000000 + (nGoodPV==51)*0.000000");
  //250 gev
  setSampleWeightFactor("HbbHbb_FullSim_250_v71","0 + (nGoodPV==1)*0.000000 + (nGoodPV==2)*0.000000 + (nGoodPV==3)*0.189076 + (nGoodPV==4)*0.476019 + (nGoodPV==5)*0.774264 + (nGoodPV==6)*0.855738 + (nGoodPV==7)*0.789157 + (nGoodPV==8)*0.701390 + (nGoodPV==9)*0.765677 + (nGoodPV==10)*0.862419 + (nGoodPV==11)*1.043154 + (nGoodPV==12)*1.126524 + (nGoodPV==13)*1.391946 + (nGoodPV==14)*1.666829 + (nGoodPV==15)*1.629645 + (nGoodPV==16)*1.286200 + (nGoodPV==17)*0.984172 + (nGoodPV==18)*0.789754 + (nGoodPV==19)*0.782547 + (nGoodPV==20)*0.795632 + (nGoodPV==21)*0.932901 + (nGoodPV==22)*0.861056 + (nGoodPV==23)*0.836711 + (nGoodPV==24)*0.979746 + (nGoodPV==25)*0.979451 + (nGoodPV==26)*0.830089 + (nGoodPV==27)*1.051759 + (nGoodPV==28)*0.649837 + (nGoodPV==29)*0.892681 + (nGoodPV==30)*1.070287 + (nGoodPV==31)*1.196061 + (nGoodPV==32)*1.043479 + (nGoodPV==33)*1.402871 + (nGoodPV==34)*1.958256 + (nGoodPV==35)*13.388514 + (nGoodPV==36)*19.901845 + (nGoodPV==37)*0.000000 + (nGoodPV==38)*0.000000 + (nGoodPV==39)*0.000000 + (nGoodPV==40)*0.000000 + (nGoodPV==41)*0.000000 + (nGoodPV==42)*0.000000 + (nGoodPV==43)*0.000000 + (nGoodPV==44)*0.000000 + (nGoodPV==45)*0.000000 + (nGoodPV==46)*0.000000 + (nGoodPV==47)*0.000000 + (nGoodPV==48)*0.000000 + (nGoodPV==49)*0.000000 + (nGoodPV==50)*0.000000 + (nGoodPV==51)*0.000000");

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

  TCut nocut = "(1)";

  TCut baseline = "cutPV==1"; 
  TCut cleaning = "passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40";
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1||passMC_PFMET150==1";
  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoPFcands10_010==0";

  TCut njets="njets20>=4&&njets20<=5";
  //TCut njets="njets20==4";
  TCut secondjet="jetpt2>50";

  TCut btagpre="CSVbest2>0.898";
  TCut btag3 = "CSVbest3>0.679";
  TCut btag4 = "CSVbest4>0.244"; 
  //  TCut btag4 = "CSVbest4>0.679";

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
    cuts.push_back(make_pair( "2 CSVT",btagpre));

    cuts.push_back(make_pair("PV",baseline));  

    cuts.push_back(make_pair( "Trigger",triggers));
    cuts.push_back(make_pair( "njets 4-5",njets));
    
  }
  else if (region.Contains("special")) {
    cuts.push_back(make_pair("no cuts",nocut));  
    cuts.push_back(make_pair("PV",baseline));  
    cuts.push_back(make_pair( "njets 4-5",njets));

  }
  else     cuts.push_back(make_pair( "Preselection",skim));


  cuts.push_back(make_pair("jet2 >50",secondjet));
  cuts.push_back(make_pair("minDeltaPhi",  TCut("minDeltaPhi20_eta5_noIdAll_nobeta>0.5 || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3 && METsig>50)")));

  cuts.push_back(make_pair( "Lepton vetoes",zl)); 
  cuts.push_back(make_pair( "iso tk veto",isotk)); 

  if (region.Contains("4b")) {
    cuts.push_back(make_pair( "3rd btag",btag3));
    cuts.push_back(make_pair( "4th btag",btag4));
  }
  else if (region.Contains("special")) {}//do nothing
  else {
    cuts.push_back(make_pair( "not 3rd btag",!btag3));
  }

  if (region.Contains("SR"))  cuts.push_back(make_pair( "Higgs masses",higgsmass));
  else  cuts.push_back(make_pair( "Higgs mass SB",higgsSB));

  cuts.push_back(make_pair( "#Delta R_{max}",drmax));
  if (!region.Contains("special")) cuts.push_back(make_pair("MET cleaning",cleaning));  


  int max= region.Contains("MET")? 300 : 150; //MET : METsig
  int step = region.Contains("MET")? 10 : 10; //was 10: 5
  for (int icut=10;icut<=max;icut+=step) {
    TString label;
    //    if (region.Contains("oldsig")) label = "METsig>";
     if (region.Contains("MET"))  label = "MET>";
    else                              label = "METsig>";
    label+=icut;
    cuts.push_back(make_pair( label,label.Data()));
  }

  //first print a header row

  //also make TGraphs of S/root(B)
  map<TString, TGraph*> SRootBGraphs;
  map<TString, TGraph*> SGraphs;

  const TString dividerHead = latexMode ? " & ":" *|* ";
  const TString divider = latexMode ? " & ":" | ";
  const TString lineEnd=latexMode ? "\\\\":" |";
  if (!latexMode) cout<<"|* ";
  cout<<"Cut \t\t"<<dividerHead;
  bool foundsignals=false;
  int nsmsamples=0;
  for (unsigned int isample = 0; isample<samples_.size(); isample++) {

    //when we get to the first signal, add a Total SM column
    if ( !isSampleSM(samples_.at(isample)) ) { 
      if (!foundsignals && nsmsamples>0)      cout<<" Total SM"<<dividerHead;
      foundsignals=true;
    }
    else nsmsamples++;

    if (!isSampleSM(samples_.at(isample))) { //signal sample: create sensitivity TGraph
      TGraph* SRootBGraph = new TGraph();
      SRootBGraph->SetName( jmt::fortranize(TString("SOverRootB_")+samples_.at(isample),"dash"));
      SRootBGraph->SetLineColor( getSampleColor( samples_.at(isample)));
      SRootBGraph->SetMarkerColor( getSampleColor( samples_.at(isample)));
      SRootBGraphs[samples_.at(isample)] = SRootBGraph;

      TGraph* SGraph = new TGraph();
      SGraph->SetName( jmt::fortranize(TString("S_")+samples_.at(isample),"dash"));
      SGraph->SetLineColor( getSampleColor( samples_.at(isample)));
      SGraph->SetMarkerColor( getSampleColor( samples_.at(isample)));
      SGraphs[samples_.at(isample)] = SGraph;
    }

    cout<<sampleLabel_[samples_.at(isample)];

    //extra column for cut efficiency
    if (doCutEff)    cout<<dividerHead<<"   ";
    
    if (isample==samples_.size()-1) cout<<lineEnd;
    else   cout<<dividerHead;
  }
  cout<<endl;


  TCut totalcut = cuts.at(0).second;
  double * lastrow = new double[samples_.size()];
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
      if ( !isSampleSM(thissamplename) ) {
	if (!foundsignals && nsmsamples>0) {
	  totalsmvalue = getIntegral("totalsm");
	  if (unweighted)	cout<<setprecision(10)<<totalsmvalue<<divider;
	  else	cout<<jmt::format_nevents(totalsmvalue,getIntegralErr("totalsm"),false,latexMode)<<divider;
	}
	foundsignals=true;
      }
      double n = getIntegral(thissamplename);
      double e = getIntegralErr(thissamplename); 
      if (unweighted)     cout<<setprecision(10)<<n;
      else     cout<<jmt::format_nevents(n,e,false,latexMode);

      //for signals, calculate S/sqrt(B)
      if (foundsignals && nsmsamples>0) {
	double soverrootb = n/sqrt(totalsmvalue);
	cout<<" ("<<setprecision(2)<<soverrootb <<") ";
	//using icut as the x-value is a bit crude but that's the best we can do for now
	SRootBGraphs[ thissamplename ]->SetPoint(SRootBGraphs[ thissamplename ]->GetN(),icut,soverrootb);
	SGraphs[ thissamplename ]->SetPoint(SGraphs[ thissamplename ]->GetN(),icut,n);
      }

      if (doCutEff) {
	cout<<divider;
	if (icut==0) cout<<"  ";
	else {
	  double eff = (lastrow[isample>0]) ? n/lastrow[isample] : 0;
	  printf(" %.2f ",eff);
	}
	lastrow[isample] = n;
      }

      if (isample==samples_.size()-1) cout<<lineEnd<<endl;
      else   cout<<divider;
    }
  }

  delete [] lastrow;

  if ( SRootBGraphs.size() >0 ) {
    TFile fcutflowoutput("higgs_cutflow_graphs.root","recreate"); //hardcode filename for now
    for (  map<TString, TGraph*>::iterator ig = SRootBGraphs.begin(); ig!= SRootBGraphs.end() ; ++ig)    ig->second->Write();
    for (  map<TString, TGraph*>::iterator ig = SGraphs.begin(); ig!= SGraphs.end() ; ++ig)    ig->second->Write();
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
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1||passMC_PFMET150==1";
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

  initHiggsSamples69(true,"hhmg200 hhmg400");
  setOutputDirectory("quickmbbPlots");

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  useTrigEff_=false;

  setStackMode(true,false,false); // stack,norm,labels
  //20 June -- update
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0 && MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; //4july add full cleanup
  //20 june -- remove the all jet trigger: ||passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1||passMC_PFMET150==1";
  TCut zl = "nMuons==0&&nElectrons==0";
  TCut tauveto="nTausLoose==0";
  TCut isotk1="nIsoTracks15_005_03==0";

  TCut njets="njets20==4 || njets20==5";
  TCut jet2="jetpt2>50";
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut mdp = "minDeltaPhi20_eta5_noIdAll_nobeta>0.5 || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3 && METsig>50)";

  TCut higgsmass = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsdiff = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";

  TCut drmax = "deltaRmax_hh<2.2";

  TCut metsig = "METsig>30";

  //no higgs mass cut
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && mdp && btag2 && btag3 && btag4 && higgsdiff && drmax && metsig;

  clearSamples(); 
  addSample("SMS-TChiHH_2b2b_2J_mChargino-350to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871_v71$400$1:higgsMbb1MassDiff_correct==0",kRed,"0 correct");
  addSample("SMS-TChiHH_2b2b_2J_mChargino-350to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871_v71$400$1:higgsMbb1MassDiff_correct==1",kBlue,"1 correct");
  addSample("SMS-TChiHH_2b2b_2J_mChargino-350to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871_v71$400$1:higgsMbb1MassDiff_correct==2",kGreen+4,"2 correct");
  selectionLabel_ = "m_{#tilde{H}} = 400 GeV";

  var="0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)"; xtitle="<m_{jj}> (GeV)";
  nbins=20; low= 0; high=200;
  setStackMode(false,false,false); 
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_massshape_mbb_400",0,"GeV");

 //make the same plots but requiring that it is possible to reconstruct the higgs (all 4 b from higgs were reconstructed)
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && mdp && btag2 && btag3 && btag4 && higgsdiff && drmax && metsig &&TCut("njetsHiggsMatch20==4");
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_massshape_mbb_400_4possible",0,"GeV");
  //another iteration but this time requiring that all 4 b from higgs are reco'd with eta<5
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && mdp && btag2 && btag3 && btag4 && higgsdiff && drmax && metsig &&TCut("njetsHiggsMatch20_eta5==4");
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_massshape_mbb_400_4possible_eta5",0,"GeV");

  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && mdp && btag2 && btag3 && btag4 && higgsdiff && drmax && metsig;
  clearSamples(); 
  addSample("SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872_v71$200$1:higgsMbb1MassDiff_correct==0",kRed,"0 correct");
  addSample("SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872_v71$200$1:higgsMbb1MassDiff_correct==1",kBlue,"1 correct");
  addSample("SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872_v71$200$1:higgsMbb1MassDiff_correct==2",kGreen+4,"2 correct");
  selectionLabel_ = "m_{#tilde{H}} = 200 GeV";

  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_massshape_mbb_200",0,"GeV");

  //make the same plots but requiring that it is possible to reconstruct the higgs (all 4 b from higgs were reconstructed)
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && mdp && btag2 && btag3 && btag4 && higgsdiff && drmax && metsig &&TCut("njetsHiggsMatch20==4");
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_massshape_mbb_200_4possible",0,"GeV");
  //now with eta5
  selection_=baseline&&triggers && zl&& tauveto && isotk1 && njets&&jet2 && mdp && btag2 && btag3 && btag4 && higgsdiff && drmax && metsig &&TCut("njetsHiggsMatch20_eta5==4");
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_massshape_mbb_200_4possible_eta5",0,"GeV");



}

void higgsmbb_backgroundStudies2() { //this time split the off-mbb peak into a couple regions

  addToSamplesAll("TTJets_FullLeptMGDecays");
  addToSamplesAll("TTJets_SemiLeptMGDecays");

  //use keith's cut flow.

  nameOfEventWeight_="weight3";//don't want an event weight

  inputPath="/cu4/ra2b/reducedTrees/v68_4/";

  loadSamples("ra2b2012");
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
  TCut triggers = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1||passMC_PFMET150";
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
    var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
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

void higgsmbb_cutflowDtry2_illustration() {

  //plots for my cornell talk


 nameOfEventWeight_="weight3";//don't want an event weight

  inputPath="/cu4/ra2b/reducedTrees/v68_4/";

  addToSamplesAll("TChihh250");
  addToSamplesAll("TChihh400");

  loadSamples("ra2b2012");
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


  var="METsig"; xtitle="E^{miss}_{T} significance #it{S}";
  nbins=15; low= 0; high=300;
  setStackMode(false,false,false); //no stack
  setPlotMinimum(0.1);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ewkhh_cutflowDb3_METsig",0);
  setLogY(false);
  drawEffVRej("TChihh400","TTbarJets0","hh400","ttbar",true);
  drawEffVRej("TChihh250","TTbarJets0","hh250","ttbar",true);

}
