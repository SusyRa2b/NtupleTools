/*
--- ROOT version ---
source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/x86_64-slc5-gcc46-opt/root/bin/thisroot.csh

--- setup --

You must have symlink to or copy of MiscUtil.cxx in the working directory.
This is available at:
https://github.com/joshmt1/UserCode/blob/master/MiscUtil.cxx

For interactive mode, you must also make a symlink as follows:
ln -s drawReducedTrees.h drawDelphes.h


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

double lumiScale_ = 300e3;//300 fb-1

#include "drawReducedTrees.h"

void initSamples(TString option="all") {

  if (option.Contains("all")) option = "tt bj signal";

  inputPath = "/cu3/joshmt/Upgrade/PhaseII_Configuration4v2_140PileUp/v02/"; //delcmscornell
    //inputPath = "/cu3/joshmt/Upgrade/PhaseI_Configuration0_140PileUp/v01/";
  //inputPath = "root://eoscms.cern.ch//eos/cms/store/user/joshmt/upgrade/v01/PhaseI/Configuration0/140PileUp/"; //direct from EOS
  dodata_=false; //no data for upgrade studies!
  cmEnergy_=14; //TeV
  reducedTreeName_ = "simpleTree";


  addToSamplesAll("tt-4p-0-600-v1510_14TEV");
  addToSamplesAll("tt-4p-600-1100-v1510_14TEV");
  addToSamplesAll("tt-4p-1100-1700-v1510_14TEV");
  addToSamplesAll("tt-4p-1700-2500-v1510_14TEV");
  addToSamplesAll("tt-4p-2500-100000-v1510_14TEV");

  addToSamplesAll("tj-4p-0-500-v1510_14TEV");
  addToSamplesAll("tj-4p-500-1000-v1510_14TEV");
  addToSamplesAll("tj-4p-1000-1600-v1510_14TEV");
  addToSamplesAll("tj-4p-1600-2400-v1510_14TEV");
  addToSamplesAll("tj-4p-2400-100000-v1510_14TEV");

  addToSamplesAll("tB-4p-0-500-v1510_14TEV");
  addToSamplesAll("tB-4p-500-900-v1510_14TEV");
  addToSamplesAll("tB-4p-900-1500-v1510_14TEV");
  addToSamplesAll("tB-4p-1500-2200-v1510_14TEV");
  addToSamplesAll("tB-4p-2200-100000-v1510_14TEV");
 
  addToSamplesAll("LL-4p-0-100-v1510_14TEV");
  addToSamplesAll("LL-4p-100-200-v1510_14TEV");
  addToSamplesAll("LL-4p-200-500-v1510_14TEV");
  addToSamplesAll("LL-4p-500-900-v1510_14TEV");
  addToSamplesAll("LL-4p-900-1400-v1510_14TEV");
  addToSamplesAll("LL-4p-1400-100000-v1510_14TEV");

  addToSamplesAll("Bj-4p-0-300-v1510_14TEV");
  addToSamplesAll("Bj-4p-1100-1800-v1510_14TEV");
  addToSamplesAll("Bj-4p-1800-2700-v1510_14TEV");
  addToSamplesAll("Bj-4p-2700-3700-v1510_14TEV");
  addToSamplesAll("Bj-4p-300-600-v1510_14TEV");
  addToSamplesAll("Bj-4p-3700-100000-v1510_14TEV");
  addToSamplesAll("Bj-4p-600-1100-v1510_14TEV");

  addToSamplesAll("B-4p-0-1-v1510_14TEV");

  TString signal1 = "susyhit_Scenario1_v02";
  addToSamplesAll(signal1);


  loadSamples("delphes");
  usePUweight_=false;
  useTrigEff_=false;


  currentConfig_=configDescriptions_.getDefault();
 clearSamples();
  resetChains(); //important in the case that we call initHiggsSamples() more than once

  if (option.Contains("signal")) {
    addSample(signal1,kRed,"Scenario 1");
    double Nsig = getTree(signal1)->GetEntries();
    setSampleScaleFactor(signal1,0.1547/Nsig ); //from Isabell (to be updated?)
    setSampleWeightFactor(signal1,"1.0/weight"); //cancel out weighting stored in tree
  }

  if (option.Contains("combinesm")) {
    TString smname = "tt-4p-0-600-v1510_14TEV";
    addSample(smname,kAzure-3,"Total SM");
    chainSamples(smname,"tt-4p-600-1100-v1510_14TEV");
    chainSamples(smname,"tt-4p-1100-1700-v1510_14TEV");
    chainSamples(smname,"tt-4p-1700-2500-v1510_14TEV");
    chainSamples(smname,"tt-4p-2500-100000-v1510_14TEV");
    
    TString tj_name = "tj-4p-0-500-v1510_14TEV"; 
    chainSamples(smname,tj_name);
    chainSamples(smname,"tj-4p-500-1000-v1510_14TEV");
    chainSamples(smname,"tj-4p-1000-1600-v1510_14TEV");
    chainSamples(smname,"tj-4p-1600-2400-v1510_14TEV");
    chainSamples(smname,"tj-4p-2400-100000-v1510_14TEV");

    /*    
    TString tb_name = "tB-4p-0-500-v1510_14TEV"; 
    chainSamples(smname,tb_name);
    chainSamples(smname,"tB-4p-500-900-v1510_14TEV");
    chainSamples(smname,"tB-4p-900-1500-v1510_14TEV");
    chainSamples(smname,"tB-4p-1500-2200-v1510_14TEV");
    chainSamples(smname,"tB-4p-2200-100000-v1510_14TEV");
    */

    TString ll_name = "LL-4p-0-100-v1510_14TEV";
    chainSamples(smname,ll_name);
    chainSamples(smname,"LL-4p-100-200-v1510_14TEV");
    chainSamples(smname,"LL-4p-200-500-v1510_14TEV");
    chainSamples(smname,"LL-4p-500-900-v1510_14TEV");
    chainSamples(smname,"LL-4p-900-1400-v1510_14TEV");
    chainSamples(smname,"LL-4p-1400-100000-v1510_14TEV");

    TString bj_name = "Bj-4p-0-300-v1510_14TEV";
    chainSamples(smname,bj_name);
    chainSamples(smname,"Bj-4p-1100-1800-v1510_14TEV");
    chainSamples(smname,"Bj-4p-1800-2700-v1510_14TEV");
    chainSamples(smname,"Bj-4p-2700-3700-v1510_14TEV");
    chainSamples(smname,"Bj-4p-300-600-v1510_14TEV");
    chainSamples(smname,"Bj-4p-3700-100000-v1510_14TEV");
    chainSamples(smname,"Bj-4p-600-1100-v1510_14TEV");

    chainSamples(smname,"B-4p-0-1-v1510_14TEV");

  }
  else {
    if (option.Contains("tt") ) {
      TString tt_name = "tt-4p-0-600-v1510_14TEV";
      addSample(tt_name,kAzure-3,"tt");
      chainSamples(tt_name,"tt-4p-600-1100-v1510_14TEV");
      chainSamples(tt_name,"tt-4p-1100-1700-v1510_14TEV");
      chainSamples(tt_name,"tt-4p-1700-2500-v1510_14TEV");
      chainSamples(tt_name,"tt-4p-2500-100000-v1510_14TEV");
    
      TString tj_name = "tj-4p-0-500-v1510_14TEV"; 
      addSample(tj_name,kMagenta,"tj");
      chainSamples(tj_name,"tj-4p-500-1000-v1510_14TEV");
      chainSamples(tj_name,"tj-4p-1000-1600-v1510_14TEV");
      chainSamples(tj_name,"tj-4p-1600-2400-v1510_14TEV");
      chainSamples(tj_name,"tj-4p-2400-100000-v1510_14TEV");
      /*      
      TString tb_name = "tB-4p-0-500-v1510_14TEV"; 
      addSample(tb_name,kBlue-3,"tB");
      chainSamples(tb_name,"tB-4p-500-900-v1510_14TEV");
      chainSamples(tb_name,"tB-4p-900-1500-v1510_14TEV");
      chainSamples(tb_name,"tB-4p-1500-2200-v1510_14TEV");
      chainSamples(tb_name,"tB-4p-2200-100000-v1510_14TEV");
      */
    }
    if (option.Contains("bj")) {
      TString ll_name = "LL-4p-0-100-v1510_14TEV";
      addSample(ll_name,kGreen,"LL");
      chainSamples(ll_name,"LL-4p-100-200-v1510_14TEV");
      chainSamples(ll_name,"LL-4p-200-500-v1510_14TEV");
      chainSamples(ll_name,"LL-4p-500-900-v1510_14TEV");
      chainSamples(ll_name,"LL-4p-900-1400-v1510_14TEV");
      chainSamples(ll_name,"LL-4p-1400-100000-v1510_14TEV");
      
      TString bj_name = "Bj-4p-0-300-v1510_14TEV";
      addSample(bj_name,kOrange,"Bj");
      chainSamples(bj_name,"Bj-4p-1100-1800-v1510_14TEV");
      chainSamples(bj_name,"Bj-4p-1800-2700-v1510_14TEV");
      chainSamples(bj_name,"Bj-4p-2700-3700-v1510_14TEV");
      chainSamples(bj_name,"Bj-4p-300-600-v1510_14TEV");
      chainSamples(bj_name,"Bj-4p-3700-100000-v1510_14TEV");
      chainSamples(bj_name,"Bj-4p-600-1100-v1510_14TEV");

      addSample("B-4p-0-1-v1510_14TEV",kOrange+3,"B");
    }
  }
}

