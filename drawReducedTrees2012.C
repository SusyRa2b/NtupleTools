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
#include "drawReducedTrees2012.h"

void SLsample_PLBref() {


  nameOfEventWeight_="weight3";//don't want an event weight

  inputPath="/cu4/ra2b/reducedTrees/v66_10/";

  addToSamplesAll("TTJets_FullLeptMGDecays");
  addToSamplesAll("TTJets_SemiLeptMGDecays");

  loadSamples(true,"ra2b2012");
  usePUweight_=true;

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  resetChains();
  clearSamples();

  addSample("TTJets_FullLeptMGDecays",kAzure,"tt 2l");
  addSample("TTJets_SemiLeptMGDecays",kAzure-3,"tt 1l");
  addSample("SingleTop");
  addSample("TTV");
  addSample("WJets");
  addSample("VV");
  addSample("ZJets");

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
  TCut sl = getSingleLeptonCut();

  selection_ = baseline&&cleaning&&ht&&met&&mdp&&sl;

  setStackMode(true,false,false);

  var="nIsoTracks15_005_03_lepcleaned"; xtitle="nIsoTracks15_005_03_lepcleaned";
  nbins = 3; low=0; high=3;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_SL_baseline_nadditionalisotracks",0,"GeV");


}

void Sijk_investigation() {

  nameOfEventWeight_="weight3";//don't want an event weight

  inputPath="./";

  addToSamplesAll("TTJets");

  loadSamples(true,"ra2b2012");
  usePUweight_=false;

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  //use a different sample order than default
  resetChains();
  clearSamples();

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


  selection_ = baseline&&cleaning && mdp  &&ht&&met; //baseline selection w/out any lepton cuts

  TString samplebase="TTJets";
  TString lTotal[2];
  lTotal[0] =samplebase + ":ttbarDecayCode==3";
  lTotal[1] =samplebase + ":ttbarDecayCode==4";
  addSample(lTotal[0]);
  addSample(lTotal[1]);
  TString lLost[2][5];
  {
    for (int ilep=0;ilep<2;ilep++) {
      for (int icode=0;icode<5;icode++) {
	lLost[ilep][icode]   =samplebase;
	TString therest; therest.Form(":ttbarDecayCode==%d&&lostLeptonCode==%d",(ilep==0) ? 3 : 4, icode+1);
	lLost[ilep][icode].Append(therest);
	addSample(lLost[ilep][icode]);
      }
    }
  }

  savePlots_=false;

  setStackMode(false,false,false);

  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 8; low=400; high=1200;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "null",0,"GeV");

  //now make ratio plots (and draw them)
  int someColors[5]={kGreen+1,kCyan+1,kViolet,kRed,kGray};
  TGraphAsymmErrors * rateLoss[2][5];
  {
    for (int ilep=0;ilep<2;ilep++) {
      renewCanvas();
      for (int icode=0;icode<5;icode++) { 
	rateLoss[ilep][icode]=new TGraphAsymmErrors();
	rateLoss[ilep][icode]->BayesDivide(getHist(lLost[ilep][icode]),getHist(lTotal[ilep]));
	rateLoss[ilep][icode]->SetMarkerColor(someColors[icode]); 
	rateLoss[ilep][icode]->SetLineColor(someColors[icode]); 
	TString optdraw = (icode==0)? "ap" : "p";
	rateLoss[ilep][icode]->Draw(optdraw);
      }
      rateLoss[ilep][0]->SetMinimum(0);
      rateLoss[ilep][0]->SetMaximum(0.5);
      TString fname;  fname.Form("lostLeptons_%s.eps",ilep==0? "e":"mu");
      thecanvas->SaveAs(fname);
    }
  }

  //now make the same plots but with a minDeltaRLeptonJet>1.5 cut
  selection_ = baseline&&cleaning && mdp  &&ht&&met &&TCut("minDeltaRLeptonJet>1.2"); //baseline selection w/out any lepton cuts
  drawPlots(var,nbins,low,high,xtitle,"Events", "null",0,"GeV");
  TGraphAsymmErrors * rateLoss2[2][5];
  {
    for (int ilep=0;ilep<2;ilep++) {
      renewCanvas();
      for (int icode=0;icode<5;icode++) { 
	rateLoss2[ilep][icode]=new TGraphAsymmErrors();
	rateLoss2[ilep][icode]->BayesDivide(getHist(lLost[ilep][icode]),getHist(lTotal[ilep]));
	rateLoss2[ilep][icode]->SetMarkerColor(someColors[icode]); 
	rateLoss2[ilep][icode]->SetLineColor(someColors[icode]); 
	TString optdraw = (icode==0)? "ap" : "p";
	rateLoss2[ilep][icode]->Draw(optdraw);
      }
      rateLoss2[ilep][0]->SetMinimum(0);
      rateLoss2[ilep][0]->SetMaximum(0.5);
      TString fname;  fname.Form("lostLeptons_highMinDRLepJet_%s.eps",ilep==0? "e":"mu");
      thecanvas->SaveAs(fname);
    }
  }


  //now study the minDeltaRLeptonJet as a function of HT

  clearSamples();
  TCut genSL="ttbarDecayCode==3||ttbarDecayCode==4";
  addSample("TTJets:HT>=400&&HT<500",kViolet,"400<HT<500 GeV");
  addSample("TTJets:HT>=500&&HT<800",kBlue,"500<HT<800 GeV");
  addSample("TTJets:HT>=800&&HT<1000",kOrange+8,"800<HT<1000 GeV");
  addSample("TTJets:HT>=1000",kRed,"HT>1000 GeV");

  selection_ = baseline&&cleaning && mdp  &&genSL &&met; //baseline selection w/out any lepton cuts

  savePlots_=true;//now we want to use the built-in plot saving

  var="minDeltaRLeptonJet"; xtitle="gen min #Delta R(lepton,j)";
  nbins = 20; low=0; high=5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_Sijk_minDeltaRlj_HTbins",0,"");
  setStackMode(false,true,false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_Sijk_minDeltaRlj_HTbins",0,"");

  //use the same splitting and plot the flavor of the nearest parton
  var="minDeltaRJetFlavor"; xtitle="nearest parton flavor";
  nbins = 10; low=0; high=10; //gluon will be in overflow bin
  setStackMode(false,false,false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_Sijk_nearestParton_HTbins",0,"");
  setStackMode(false,true,false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_Sijk_nearestParton_HTbins",0,"");

  //try to  split between reco and non-reco leptons
  clearSamples();
  addSample("TTJets:(HT>=400&&HT<500)&&(nMuons+nElectrons==1)",kViolet,"400<HT<500 GeV (SL)");
  addSample("TTJets:(HT>=400&&HT<500)&&(nMuons+nElectrons==0)",kViolet-9,"400<HT<500 GeV (ZL)");

  addSample("TTJets:(HT>=500&&HT<800)&&(nMuons+nElectrons==1)",kBlue,"500<HT<800 GeV(SL)");
  addSample("TTJets:(HT>=500&&HT<800)&&(nMuons+nElectrons==0)",kBlue-9,"500<HT<800 GeV (ZL)");

  addSample("TTJets:(HT>=800&&HT<1000)&&(nMuons+nElectrons==1)",kOrange+9,"800<HT<1000 GeV(SL)");
  addSample("TTJets:(HT>=800&&HT<1000)&&(nMuons+nElectrons==0)",kOrange-4,"800<HT<1000 GeV (ZL)");

  addSample("TTJets:(HT>=1000)&&(nMuons+nElectrons==1)",kRed,"HT>1000 GeV(SL)");
  addSample("TTJets:(HT>=1000)&&(nMuons+nElectrons==0)",kRed-9,"HT>1000 GeV (ZL)");

  var="minDeltaRLeptonJet"; xtitle="gen min #Delta R(lepton,j)";
  nbins = 10; low=0; high=5;
  setStackMode(false,false,false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_Sijk_minDeltaRlj_HTbins_0l1l",0,"");
  setStackMode(false,true,false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_Sijk_minDeltaRlj_HTbins_0l1l",0,"");

  //use the same splitting and plot the flavor of the nearest parton
  var="minDeltaRJetFlavor"; xtitle="nearest parton flavor";
  nbins = 10; low=0; high=10; //gluon will be in overflow bin
  setStackMode(false,false,false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_Sijk_nearestParton_HTbins_0l1l",0,"");
  setStackMode(false,true,false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_Sijk_nearestParton_HTbins_0l1l",0,"");


}

void compT1ttttMGVpythia_check() {

  //loop over scanSMSngen for the MG sample

  // for each point !=0, see if there is a corresponding Pythia point

  //if yes, plot a few things
  //save each plot


  inputPath = "/cu4/ra2b/reducedTrees/v68_2/";
  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  addToSamplesAll("T1tttt-UCSB1741");

  loadSamples(true,"ra2b2012");

  clearSamples();

  addSample("T1tttt-UCSB1741",kRed,"MG");
  addSample("T1tttt",kBlue,"Pythia");

  dodata_=false;
  int nbins;
  float low,high;
  TString var,xtitle;

  overrideSMSlabels_=true;
  setStackMode(false,true,false);

  drawFilenameOnPlot_=true;

  selection_="(1)";

  loadScanSMSngen("T1tttt-UCSB1741");
  int nbinsx=scanSMSngen->GetNbinsX();
  int nbinsy=scanSMSngen->GetNbinsY();
  for (int ix=1;ix<=nbinsx; ix++) {
    for (int iy=1;iy<=nbinsy; iy++) {

      loadScanSMSngen("T1tttt-UCSB1741");
      if (scanSMSngen->GetBinContent(ix,iy)>0) {

	loadScanSMSngen("T1tttt");
	if (scanSMSngen->GetBinContent(ix,iy)>0) {

	  m0_=scanSMSngen->GetXaxis()->GetBinLowEdge(ix);
	  m12_=scanSMSngen->GetYaxis()->GetBinLowEdge(iy);

	  TString outfilename;
	  outfilename.Form("T1tttt_MG_Pythia_HT-%d-%d",m0_,m12_);

	  var="HT"; xtitle="HT (GeV)";
	  nbins = 40; low=0; high=2400;
	  //drawPlots(var,nbins,low,high,xtitle,"Events", outfilename,0);
	}
	else {
	  cout<<" no pythia for "<<ix<<" "<<iy<<endl;
	}
      }
    }
  }


}

void t1t1t() { //compare t1t1t to t1tttt

  inputPath = "/cu4/ra2b/reducedTrees/signalComparisons/";
  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  addToSamplesAll("SMS-MadGraph_T1t1t-1000-200to750-100to650_UCSB1753reshuf_v68");
  addToSamplesAll("SMS-MadGraph_T1tttt-775to1075-25to500_UCSB1741reshuf_v68.1025");
  addToSamplesAll("SMS-MadGraph_T1tttt-1100to1400-25to500_UCSB1732reshuf_v68.1150");
  addToSamplesAll("SMS-MadGraph_T5tttt-1075to1175-225to1025-50_UCSB1750reshuf_v68.1125");

  //how about 

  loadSamples(true,"ra2b2012");

  clearSamples();

  //T1tttt -- 1024, 325
  addSample("SMS-MadGraph_T1tttt-775to1075-25to500_UCSB1741reshuf_v68.1025$1025$325",kBlue,"1025, 325");
  addSample("SMS-MadGraph_T1t1t-1000-200to750-100to650_UCSB1753reshuf_v68$1000$500$325",kRed,"1000, 500, 325");

  drawMCErrors_=true;
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;

  dodata_=false;
  int nbins;
  float low,high;
  TString var,xtitle;

  overrideSMSlabels_=true;
  setStackMode(false,false,false);

  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40";
  TCut ht = "HT>=400";
  TCut met = "MET>=125";
  TCut mdp = "minDeltaPhiN_asin>=4";
  TCut ldp = "minDeltaPhiN_asin<4";
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");
  TCut zlwoisotk = getLeptonVetoCut();

  selection_=baseline&&cleaning&&ht&&met&&mdp&&zl;

  var="njets"; xtitle="n jets pT>50 GeV";
  nbins = 13; low=0; high=13;
  drawPlots(var,nbins,low,high,xtitle,"Events", "T1ttttT1t1t-comp_"+var,0);

  var="HT"; xtitle="HT (GeV)";
  nbins = 40; low=400; high=2200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "T1ttttT1t1t-comp_"+var,0);

  var="MET"; xtitle="MET (GeV)";
  nbins = 40; low=100; high=800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "T1ttttT1t1t-comp_"+var,0);

  var="jetpt1"; xtitle="lead jet pT (GeV)";
  nbins = 50; low=0; high=1000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "T1ttttT1t1t-comp_"+var,0);

  clearSamples();
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  addSample("SMS-MadGraph_T1tttt-1100to1400-25to500_UCSB1732reshuf_v68.1150$1150$25",kBlue,"1150, 25");
  addSample("SMS-MadGraph_T5tttt-1075to1175-225to1025-50_UCSB1750reshuf_v68.1125$1125$225$50",kRed+3,"1125, 225, 50");
  addSample("SMS-MadGraph_T5tttt-1075to1175-225to1025-50_UCSB1750reshuf_v68.1125$1125$700$50",kRed-7,"1125, 700, 50");


  var="njets"; xtitle="n jets pT>50 GeV";
  nbins = 14; low=0; high=14;
  drawPlots(var,nbins,low,high,xtitle,"Events", "T1ttttT5tttt-comp_"+var,0);

  var="HT"; xtitle="HT (GeV)";
  nbins = 40; low=400; high=2200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "T1ttttT5tttt-comp_"+var,0);

  var="MET"; xtitle="MET (GeV)";
  nbins = 40; low=100; high=800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "T1ttttT5tttt-comp_"+var,0);

  selection_="(1)";
  var="MET"; xtitle="MET (GeV)";
  nbins = 40; low=0; high=800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "T1ttttT5tttt-nocuts_"+var,0);

  var="HT"; xtitle="HT (GeV)";
  nbins = 40; low=0; high=3000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "T1ttttT5tttt-nocuts_"+var,0);

  setStackMode(false,true,false);

  selection_="(1)";
  var="MET"; xtitle="MET (GeV)";
  nbins = 40; low=0; high=800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "T1ttttT5tttt-nocuts_"+var,0);

  var="HT"; xtitle="HT (GeV)";
  nbins = 40; low=0; high=3000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "T1ttttT5tttt-nocuts_"+var,0);


}

void ARC2_MTcut() {

  //question: what is eff of MT cut on backgrounds and signal?
  //i
  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  loadSamples(true,"ra2b2012");
  usePUweight_=true;

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  //use only some samples
  resetChains();
  clearSamples();

  addSample("TTbarJets");
  addSample("SingleTop");
  addSample("WJets");

  TString t4t1="T1tttt$1100$100";
  TString t4t2="T1tttt$850$350";

  addSample(t4t1,kRed);
  addSample(t4t2,kRed+5);
  useMassInLegend_=true;


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
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");
  TCut sl = getSingleLeptonCut();
  TCut slNoMT = "nElectrons+nMuons==1";

 
  //define analysis bins
  float htedges[]={400,500,800,1000,100000};
  float metedges[]={125,150,250,350,100000};

  TH1D* hmteff_sm[3];
  TH1D* hmteff_t4t1[3];
  TH1D* hmteff_t4t2[3];
  savePlots_=false;
  setStackMode(false,false,false); //important not to have a stack
  //loop over analysis bins
  for (int ib=0;ib<3;ib++) {
    TString histoname = "mteffhist_";
    hmteff_sm[ib] = new TH1D(histoname + TString(ib+1),histoname,16,0,16);
    hmteff_t4t1[ib] = new TH1D(histoname+TString("T1tttt1_") + TString(ib+1),histoname,16,0,16);
    hmteff_t4t2[ib] = new TH1D(histoname+TString("T1tttt2_") + TString(ib+1),histoname,16,0,16);

    int htmetbin=0;
    for (int iht=0;iht<4;iht++) {
      for (int imet=0;imet<4;imet++) {

	//in each bin: count events in SL region with and without MT cut

	var="njets"; xtitle="n jets pT>50 GeV";
	nbins = 20; low=0; high=20;

	float minht=htedges[iht];
	float maxht=htedges[iht+1];
	float minmet=metedges[imet];
	float maxmet=metedges[imet+1];
	TString htcut,metcut;
	htcut.Form("HT>=%f && HT<%f",minht,maxht);
	metcut.Form("MET>=%f && MET<%f",minmet,maxmet);

	//normal SL selection with MT cut
	selection_ = baseline&&cleaning&&mdp&&sl&&TCut(htcut.Data())&&TCut(metcut.Data());
	if (ib+1 ==1) btagSFweight_="prob1";
	else if (ib+1 ==2) btagSFweight_="prob2";
	else if (ib+1 ==3) btagSFweight_="probge3";
	else assert(0);

	drawPlots(var,nbins,low,high,xtitle,"Events", "d",0);
	float totalsmv = getIntegral("totalsm");
	float t4t1v=getIntegral(t4t1);
	float t4t2v=getIntegral(t4t2);

	//without MT cut
	selection_ = baseline&&cleaning&&mdp&&slNoMT&&TCut(htcut.Data())&&TCut(metcut.Data());
	drawPlots(var,nbins,low,high,xtitle,"Events", "d",0);
	float totalsmv_nocut = getIntegral("totalsm");
	float t4t1v_nocut=getIntegral(t4t1);
	float t4t2v_nocut=getIntegral(t4t2);
	//store withCut/Total
	htmetbin++;
	hmteff_sm[ib]->SetBinContent(htmetbin, totalsmv_nocut>0 ? totalsmv/totalsmv_nocut : 0);
	hmteff_t4t1[ib]->SetBinContent(htmetbin, t4t1v_nocut>0 ? t4t1v/t4t1v_nocut : 0);
	hmteff_t4t2[ib]->SetBinContent(htmetbin, t4t2v_nocut>0 ? t4t2v/t4t2v_nocut : 0);


      }
    }
  }

  TFile fout("mtcuthistos.root","RECREATE");
  for (int ib=0;ib<3;ib++) {
    hmteff_sm[ib]->Write();
    hmteff_t4t1[ib]->Write();
    hmteff_t4t2[ib]->Write();
  }
  fout.Close();

}

void T1bbbbcomp() {

  TString isrweight0="((SUSY_recoilPt<=120) + (SUSY_recoilPt>120&&SUSY_recoilPt<=150)*0.95 + (SUSY_recoilPt>150&&SUSY_recoilPt<=250)*0.9+ (SUSY_recoilPt>250)*0.8)";
  TString isrweightD="((SUSY_recoilPt<=120) + (SUSY_recoilPt>120&&SUSY_recoilPt<=150)*0.90 + (SUSY_recoilPt>150&&SUSY_recoilPt<=250)*0.8+ (SUSY_recoilPt>250)*0.6)";
  inputPath = "/cu4/ra2b/reducedTrees/v68_1/";

  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples
  addToSamplesAll("SMS-MadGraph_T1bbbb-1125to1400-0to500_UCSB1712reshuf_v68");
  addToSamplesAll("SMS-MadGraph_T1bbbb-775to1100-525to850_UCSB1703reshuf_v68.850");

  loadSamples(true,"ra2b2012");
  usePUweight_=true;

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  //use a different sample order than default
  resetChains();
  clearSamples();

  //points -- look at 850,600
  // also 1325,200

  //showing: Pythia, ISR0, ISR+, ISR-

  addSample("SMS-MadGraph_T1bbbb-1125to1400-0to500_UCSB1712reshuf_v68$1325$200",kBlue,"MG high #Delta");
  addSample("T1bbbb$1325$200",kRed,"Pythia high #Delta");
  //  addSample("SMS-MadGraph_T1bbbb-775to1100-525to850_UCSB1703reshuf_v68.850$850$600",kBlue,"MG low #Delta");
  //  addSample("T1bbbb$850$600",kRed,"Pythia low #Delta");
  useMassInLegend_=false;
  overrideSMSlabels_=true;

  //  setSampleScaleFactor("SMS-MadGraph_T1bbbb-1125to1400-0to500_UCSB1712reshuf_v68$1325$200",10);
  //  setSampleScaleFactor("T1bbbb$1325$200",10);

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0.3; ratioMax = 1.7;
  dodata_=false;

  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40";
  TCut ht = "HT>=400";
  TCut met = "MET>=125";
  TCut mdp = "minDeltaPhiN_asin>=4";
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");
  TCut sl = getSingleLeptonCut();

  TCut btag="nbjets>=1";

  setStackMode(false,false,false);
  if (false) {
  selection_="(1)";
  var="njets"; xtitle="n jets pT>50 GeV";
  nbins = 14; low=0; high=14;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_njets",0);

  var="njets30"; xtitle="n jets pT>30 GeV";
  nbins = 18; low=0; high=18;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_njets30",0);

  var="MET"; xtitle="Type I PF E_{T}^{miss}(GeV)";
  nbins = 100; low=0; high=1000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_met",0,"GeV");

  var="HT"; xtitle="H_{T}";
  nbins = 100; low=0; high=2000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_ht",0,"GeV");

  var="SUSY_recoilPt"; xtitle="ISR boost pT (GeV)";
  setLogY(true);
  nbins = 100; low=0; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_recoilpt",0,"GeV");
  setLogY(false);

  var="nbjets"; xtitle="n (CSVM) b jets pT>50 GeV";
  nbins = 6; low=0; high=6;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_nbjets",0);

  setLogY(true);
  var="jetpt1"; xtitle="lead jet pT (GeV)";
  nbins = 100; low=0; high=2000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_jetpt1",0,"GeV");

  var="jetpt2"; xtitle="2nd jet pT (GeV)";
  nbins = 100; low=0; high=1000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_jetpt2",0,"GeV");

  var="jetpt3"; xtitle="3rd jet pT (GeV)";
  nbins = 100; low=0; high=1000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_jetpt3",0,"GeV");

  setLogY(false);
  var="minDeltaPhiN_asin"; xtitle="normalized min#Delta#phi";
  nbins = 40; low=0; high=40;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_minDeltaPhiN",0);

  var="nElectrons"; xtitle="n electrons pT>10 GeV";
  nbins = 5; low=0; high=5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_nelectrons",0);

  var="nMuons"; xtitle="n muons pT>10 GeV";
  nbins = 5; low=0; high=5;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_nmuons",0);

  var="nGoodPV"; xtitle="good PVs";
  nbins = 40; low=0; high=40;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_ngoodPV",0);

  // selection_=baseline&&cleaning&&ht&&met&&mdp&&zl&&btag;
  //could plot lepton momenta also
  selection_="nMuons==1";

  setLogY(true);
  var="muonpt1"; xtitle="lead muon pT (GeV)";
  nbins = 50; low=0; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_muonpt1",0,"GeV");

  selection_="nElectrons==1";
  var="eleet1"; xtitle="lead electron pT (GeV)";
  nbins = 50; low=0; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_eleet1",0,"GeV");

  selection_="nMuons==1";
  setLogY(false);
  var="muonpt1+MET"; xtitle="MET + muon pT (GeV)";
  nbins = 100; low=0; high=1000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_muMET",0,"GeV");

  selection_="nElectrons==1";
  setLogY(false);
  var="eleet1+MET"; xtitle="MET + electron pT (GeV)";
  nbins = 100; low=0; high=1000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_eleMET",0,"GeV");

  selection_="nMuons+nElectrons==1";
  setLogY(false);
  var="MT_Wlep"; xtitle="leptonic transverse mass (GeV)";
  nbins = 90; low=0; high=900;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_MT",0,"GeV");
  }

  //now with cuts
  selection_ = baseline&&cleaning&&ht&&met&&mdp&&zl;

  btagSFweight_="probge1";
  setSampleWeightFactor("SMS-MadGraph_T1bbbb-1125to1400-0to500_UCSB1712reshuf_v68$1325$200",isrweight0);
  //  setSampleWeightFactor("SMS-MadGraph_T1bbbb-775to1100-525to850_UCSB1703reshuf_v68.850$850$600",isrweight0);

  if (false) {
  var="njets"; xtitle="n jets pT>50 GeV";
  nbins = 14; low=0; high=14;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp-baseline_njets-isrweight0",0);

  var="MET"; xtitle="Type I PF E_{T}^{miss}(GeV)";
  nbins = 40; low=125; high=525;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp-baseline_met",0,"GeV");

  var="HT"; xtitle="H_{T}";
  nbins = 100; low=400; high=1400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp-baseline_ht",0,"GeV");

  }

  usePUweight_=false; // just for clarity
  var="SUSY_recoilPt"; xtitle="ISR boost pT (GeV)";
  setLogY(true);
  nbins = 18; low=0; high=540;
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp-baseline_recoilpt-highdm",0,"GeV");

    btagSFweight_="(1)";
    selection_="(1)";
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp-nocuts_recoilpt-highdm",0,"GeV");

    clearSamples();
    addSample("SMS-MadGraph_T1bbbb-775to1100-525to850_UCSB1703reshuf_v68.850$850$600",kBlue,"MG low #Delta");
    addSample("T1bbbb$850$600",kRed,"Pythia low #Delta");
setSampleWeightFactor("SMS-MadGraph_T1bbbb-775to1100-525to850_UCSB1703reshuf_v68.850$850$600",isrweight0);


  btagSFweight_="probge1";
  selection_ = baseline&&cleaning&&ht&&met&&mdp&&zl;

  setSampleWeightFactor("SMS-MadGraph_T1bbbb-1125to1400-0to500_UCSB1712reshuf_v68$1325$200",isrweight0);
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp-baseline_recoilpt-lowdm",0,"GeV");

    btagSFweight_="(1)";
    selection_="(1)";
    drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp-nocuts_recoilpt-lowdm",0,"GeV");


  return;
  //
  var="SUSY_recoilPt"; xtitle="ISR boost pT (GeV)";
  setLogY(true);
  nbins = 50; low=0; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp-baseline_recoilpt",0,"GeV");

  setLogY(false);

  var="nbjets"; xtitle="n (CSVM) b jets pT>50 GeV";
  nbins = 6; low=0; high=6;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp-baseline_nbjets",0);


/*
  //some more gen-level stuff
  setLogY(false);
  selection_="(1)";
  var="SUSY_gluino_pt[0]"; xtitle="gluino1 pT (GeV)";
  nbins = 100; low=0; high=1200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1bbbb-MGcomp_gluinoPt1",0);
  var="SUSY_gluino_pt[1]"; xtitle="gluino2 pT (GeV)";
  nbins = 100; low=0; high=1200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1tttt-MGcomp_gluinoPt2",0);

  var="SUSY_top_pt[0]"; xtitle="top1 pT (GeV)";
  nbins = 100; low=0; high=1200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1tttt-MGcomp_topPt1",0);
  var="SUSY_top_pt[1]"; xtitle="top2 pT (GeV)";
  nbins = 100; low=0; high=1200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1tttt-MGcomp_topPt2",0);

  var="SUSY_topbar_pt[0]"; xtitle="topbar1 pT (GeV)";
  nbins = 100; low=0; high=1200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1tttt-MGcomp_topbarPt1",0);
  var="SUSY_topbar_pt[1]"; xtitle="topbar2 pT (GeV)";
  nbins = 100; low=0; high=1200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1tttt-MGcomp_topbarPt2",0);


  var="SUSY_chi0_pt[0]"; xtitle="neutralino 1 pT (GeV)";
  nbins = 100; low=0; high=1200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1tttt-MGcomp_chi0Pt1",0);
  var="SUSY_chi0_pt[1]"; xtitle="neutralino 2 pT (GeV)";
  nbins = 100; low=0; high=1200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1tttt-MGcomp_chi0Pt2",0);

  //1d dalitz projections
  var="sqrt(SUSY_msq12[0])"; xtitle="m_{12} (GeV)";
  nbins=100; low=0; high=1200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1tttt-MGcomp_m12",0);
  var="sqrt(SUSY_msq12[1])"; xtitle="m_{12} [2] (GeV)";
  nbins=100; low=0; high=1200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1tttt-MGcomp_m12_2",0);

  var="sqrt(SUSY_msq23[0])"; xtitle="m_{23} (GeV)";
  nbins=100; low=0; high=1200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1tttt-MGcomp_m23",0);
  var="sqrt(SUSY_msq23[1])"; xtitle="m_{23} [2] (GeV)";
  nbins=100; low=0; high=1200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_T1tttt-MGcomp_m23_2",0);

  // draw2d("SUSY_msq12[0]",50,0,700000,"SUSY_msq23[0]",50,0,700000,"m^{2}_{top-LSP} (GeV^{2})","m^{2}_{topbar-LSP} (GeV^{2})","ra2b_T1tttt-MGcomp_
  renewCanvas();
  getTree("T1ttttPythia")->Draw("SUSY_msq23[0]:SUSY_msq12[0]>>hpythia2d","m0==1150 && m12==50","colz");
  renewCanvas();
  getTree("T1ttttMG")->Draw("SUSY_msq23[0]:SUSY_msq12[0]>>hpythia2d","m0==1150 && m12==50","colz");
*/

}

void quickmbb() { //study high mass mbb!


  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  loadSamples(true,"ra2b2012");
  usePUweight_=true;

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  //use a different sample order than default
  resetChains();
  clearSamples();

  addSample("VV");
  addSample("TTV");
  addSample("ZJets");

  addSample("Zinvisible");
  addSample("PythiaPUQCD");
  addSample("TTbarJets:ttbarDecayCode==2",kYellow+3,"t#bar{t} (had)");
  addSample("SingleTop");
  addSample("WJets");
  addSample("TTbarJets:ttbarDecayCode!=2",kAzure-3,"t#bar{t} (#geq 1l)");

  //    addSample("T1bbbb$1200$250",kRed);
  //  useMassInLegend_=true;

  //use the 'owen factors' verbatim
  // for Fig 3 of the PAS we would need a factor from the fit
  //  setSampleScaleFactor("TTbarJets:ttbarDecayCode!=2",0.9);
  //  setSampleScaleFactor("TTbarJets:ttbarDecayCode==2",0.9);
  //  setSampleScaleFactor("WJets",0.9);
  //  setSampleScaleFactor("PythiaPUQCD",1.8);

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
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");
  TCut sl = getSingleLeptonCut();

  setDatasetToDraw("2012hybridplus");
  TCut ge2b = "nbjets>=2";
  TCut ge3b = "nbjets>=3";


  //for MET plot, add vertical lines
  //  addVerticalLine(150);
  //  addVerticalLine(250);
  //  addVerticalLine(350);

  maxScaleFactor_=3;

  selection_=baseline&&cleaning&&ht&&met&&mdp&&zl&&ge2b;
  var="mjj_bestTwoCSV"; xtitle="m_{bb} (best CSV)";
  //  btagSFweight_="probge1";
  nbins = 50; low=0; high=1200;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b-mbb_bestTwoCSV_ge2b",0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b-mbb_bestTwoCSV_ge2b",0,"GeV");

  //  btagSFweight_="probge3";
  selection_=baseline&&cleaning&&ht&&met&&mdp&&zl&&ge3b;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b-mbb_bestTwoCSV_ge3b",0,"GeV");

  var="mjj_closestB"; xtitle="m_{bb} (closest DR)";
  //  btagSFweight_="probge1";
  selection_=baseline&&cleaning&&ht&&met&&mdp&&zl&&ge2b;
  nbins = 50; low=0; high=1200;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b-mbb_closestDR_ge2b",0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b-mbb_closestDR_ge2b",0,"GeV");

  //  btagSFweight_="probge3";
  selection_=baseline&&cleaning&&ht&&met&&mdp&&zl&&ge3b;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b-mbb_closestDR_ge3b",0,"GeV");


}

void prescaled() { //turned out to be a fool's errand. there are not enough stats in the data.

  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  loadSamples(true,"ra2b2012");
  usePUweight_=true;

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  //use a different sample order than default
  resetChains();
  clearSamples();

  addSample("VV");
  addSample("TTV");
  addSample("ZJets");

  addSample("Zinivisible");
  addSample("PythiaPUQCD");
  addSample("TTbarJets");
  addSample("SingleTop");
  addSample("WJets");

  //use the 'owen factors' verbatim
  // for Fig 3 of the PAS we would need a factor from the fit
  setSampleScaleFactor("TTbarJets",0.9);
  setSampleScaleFactor("WJets",0.9);
  setSampleScaleFactor("PythiaPUQCD",1.8);

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 1.5;
  dodata_=true;

  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40";
  TCut ht = "HT>=400";
  TCut met = "MET>=125";
  TCut mdp = "minDeltaPhiN_asin>=4";
  TCut ldp = "minDeltaPhiN_asin<4";
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");
  TCut sl = getSingleLeptonCut();

  setDatasetToDraw("2012hybridplus");

  selection_=baseline&&cleaning&&ht&&met&&ldp&&zl;
  btagSFweight_="prob1";
  var="HT"; xtitle="HT";
  nbins=32; low=400;high=2000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ldp_HT_eq1b-notrig",0);
  useTrigEff_=true;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ldp_HT_eq1b-oldtrig",0);

  setDatasetToDraw("JetHT");
  lumiScale_=20000*0.001;
  selection_=TCut("cutPV==1 && (pass_PFHT350==1||pass_PFNoPUHT350==1) && njets>=3 && jetpt2>=70") &&cleaning&&ht&&met&&ldp&&zl;
  useTrigEff_=false;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ldp_HT_eq1b-PFHT350",0);

}

void MTbdataMC() {

  inputPath = "/cu4/ra2b/reducedTrees/v66_8/"; //switch to v66_8
  dataInputPath =  "/cu4/ra2b/reducedTrees/v66_8/ORIGINALS/";
  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  loadSamples(true,"ra2b2012");
  usePUweight_=true;

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  clearSamples();

  addSample("TTbarJets");
  addSample("WJets");
  addSample("SingleTop");
  addSample("Zinvisible");
  addSample("VV",kCyan+1,"other");
  chainSamples("VV","TTV");
  chainSamples("VV","ZJets");

  setSampleScaleFactor("TTbarJets",0.9);//add the owen factors
  setSampleScaleFactor("WJets",0.9);//add the owen factors
  setSampleScaleFactor("SingleTop",0.9);//add the owen factors

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;
  setDatasetToDraw("2012hybridplus");

  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40";
  TCut ht = "HT>=400";
  TCut met = "MET>=125";
  TCut mdp = "minDeltaPhiN_asin>=4";
  TCut ldp = "minDeltaPhiN_asin<4";
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");
  TCut sl = getSingleLeptonCut();
  TCut btag = "nbjets>=1";

  setStackMode(true);

  //MT_bestCSV
  selection_=baseline && cleaning && ht && btag && met && mdp && sl;
  var="MT_bestCSV"; xtitle="MT(b-like,MET)";
  nbins=30; low=0;high=600;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_mtcsv_SL_inclusive",0);

  selection_=baseline && cleaning && ht && btag && TCut("MET>=200") && mdp && sl;
  var="MT_bestCSV"; xtitle="MT(b-like,MET)";
  nbins=30; low=0;high=600;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_mtcsv_SL_met200",0);

    addSample("PythiaPUQCD"); //save time
    setSampleScaleFactor("PythiaPUQCD",1.8);

  //LDP
  selection_=baseline && cleaning && ht && btag && TCut("MET>=200") && ldp && zl;
  var="MT_bestCSV"; xtitle="MT(b-like,MET)";
  nbins=30; low=0;high=600;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_mtcsv_LDP_met200",0);

  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_mtcsv_LDP_met200",0);
  //

}

void T2ttstudies() {

  inputPath = "/cu4/ra2b/reducedTrees/v66_8/"; //switch to v66_8

  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  loadSamples(true,"ra2b2012");
  usePUweight_=false;

  useTrigEff_=false;

  
  currentConfig_=configDescriptions_.getDefault();


  //use a different sample order than default
  //resetChains();
  clearSamples();

  addSample("TTbarJets",kBlue,"t#bar{t}+W+t");
  chainSamples("TTbarJets","WJets");
  chainSamples("TTbarJets","SingleTop");
  addSample("PythiaPUQCD",kGreen);
  addSample("Zinvisible");

  addSample("T2tt$600$0",kRed);
  addSample("T2tt$500$150",kRed-3);
  addSample("T2tt$400$125",kRed-7);
  useMassInLegend_=true;

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=false;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40";
  TCut ht = "HT>=400";
  TCut met = "MET>=125";
  TCut mdp = "minDeltaPhiN_asin>=4";
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");
  TCut btag = "nbjets>=1";

  setStackMode(false,true,false);

  selection_=TCut("cutPV==1 && cutTrigger2==1 && jetpt2>=70") && cleaning && ht && btag && met && mdp && zl;
  var="njets"; xtitle="jet multiplicity";
  nbins=8; low=2;high=10;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_njets",0);

  selection_=baseline && cleaning && btag && met && mdp && zl;
  var="HT"; xtitle="HT";
  nbins=20; low=200;high=800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_ht",0);


  selection_=baseline && cleaning && btag && ht && mdp && zl;
  var="MET"; xtitle="MET";
  nbins=30; low=125;high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_met",0);


  selection_=baseline && cleaning && ht && mdp && zl;
  var="nbjets"; xtitle="n b tags";
  nbins=4; low=0;high=4;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_nb",0);

  resetChains();
  clearSamples();
  addSample("TTbarJets",kBlue,"t#bar{t}");
  addSample("WJets");
  //addSample("SingleTop");
  addSample("PythiaPUQCD",kGreen);
  addSample("Zinvisible");

  addSample("T2tt$600$0",kRed);

  selection_=baseline && cleaning && ht && mdp && zl;
  var="nbjets"; xtitle="n b tags";
  nbins=4; low=0;high=4;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_nb-split",0);


  resetChains();
  clearSamples();
  addSample("TTbarJets",kBlue,"t#bar{t}+W+t");
  chainSamples("TTbarJets","WJets");
  chainSamples("TTbarJets","SingleTop");
  addSample("PythiaPUQCD",kGreen);
  addSample("Zinvisible");
  addSample("T2tt$600$0",kRed);
  addSample("T2tt$500$150",kRed-3);
  addSample("T2tt$400$125",kRed-7);


  selection_=baseline && cleaning && ht && btag && met && mdp && zl;
  var="MT_b"; xtitle="MT(b,MET)";
  nbins=60; low=0;high=800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_mtb",0);


  //
  drawEffVRej("T2tt$600$0","totalsm","T2tt (600,0)","total SM",true);
  TGraphAsymmErrors* g_t2tt_600_0 = (TGraphAsymmErrors*) effgraph->Clone("g_t2tt_600_0");

  drawEffVRej("T2tt$400$125","totalsm","T2tt (400,125)","total SM",true);
  TGraphAsymmErrors* g_t2tt_400_125 = (TGraphAsymmErrors*) effgraph->Clone("g_t2tt_400_125");

  //MT_bestCSV
  selection_=baseline && cleaning && ht && btag && met && mdp && zl;
  var="MT_bestCSV"; xtitle="MT(b-like,MET)";
  nbins=60; low=0;high=800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_mtcsv",0);
  //
  drawEffVRej("T2tt$600$0","totalsm","T2tt (600,0)","total SM",true);
  TGraphAsymmErrors* g_t2tt_mtcsv_600_0 = (TGraphAsymmErrors*) effgraph->Clone("g_t2tt_mtcsv_600_0");

  drawEffVRej("T2tt$400$125","totalsm","T2tt (400,125)","total SM",true);
  TGraphAsymmErrors* g_t2tt_mtscv_400_125 = (TGraphAsymmErrors*) effgraph->Clone("g_t2tt_mtcsv_400_125");

  //new plots trying to get at jim's request

  //lumi normalized, log scale
  resetChains();
  clearSamples();
  setStackMode(false,false,false);
  addSample("TTbarJets",kBlue,"t#bar{t}+W+t + QCD + Z #nu #nu");
  chainSamples("TTbarJets","WJets");
  chainSamples("TTbarJets","SingleTop");
  chainSamples("TTbarJets","PythiaPUQCD");
  chainSamples("TTbarJets","Zinvisible");
  addSample("T2tt$600$0",kRed);
  addSample("T2tt$500$150",kRed-3);
  addSample("T2tt$400$125",kRed-7);


  selection_=baseline && cleaning && ht && btag && met && mdp && zl;
  var="MT_b"; xtitle="MT(b,MET)";
  setLogY(true);
  nbins=60; low=0;high=800;
  setPlotMinimum(0.1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_mtb-joinedsm",0);
  resetPlotMinimum();

  //special met cut, no log scale
  selection_=baseline && cleaning && ht && btag && TCut("MET>200") && mdp && zl;
  var="MT_b"; xtitle="MT(b,MET)";
  setLogY(false);
  nbins=60; low=0;high=800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_mtb_met200",0);
  setPlotMaximum(150);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_mtb_met200_max",0);
  //why can't I see
  drawEffVRej("T2tt$600$0","TTbarJets","T2tt (600,0)","total SM",true);

  //now try unit normalized linear scale
  resetPlotMaximum();
  selection_=baseline && cleaning && ht && btag && met && mdp && zl;
  setStackMode(false,true,false);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_mtb_totalsm",0);


  //same thing but   var="MT_bestCSV"; xtitle="MT(b-like,MET)";
  selection_=baseline && cleaning && ht && btag && met && mdp && zl;
  var="MT_bestCSV"; xtitle="MT(b-like,MET)";
  setLogY(true);
  nbins=60; low=0;high=800;
  setPlotMinimum(0.1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_mtcsv-joinedsm",0);
  resetPlotMinimum();

  //special met cut, no log scale
  selection_=baseline && cleaning && ht && btag && TCut("MET>200") && mdp && zl;
  var="MT_bestCSV"; xtitle="MT(b-like,MET)";
  setLogY(false);
  nbins=60; low=0;high=800;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_mtcsv_met200",0);
  setPlotMaximum(150);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_mtcsv_met200_max",0);
  //why can't I see
  drawEffVRej("T2tt$600$0","TTbarJets","T2tt (600,0)","total SM",true);

  //now try unit normalized linear scale
  resetPlotMaximum();
  selection_=baseline && cleaning && ht && btag && met && mdp && zl;
  setStackMode(false,true,false);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_mtcsv_totalsm",0);

  drawEffVRej("T2tt$500$150","TTbarJets","T2tt (500,150)","total SM",true);


  //old plots
  clearSamples();
  addSample("TTbarJets:ttbarDecayCode!=2&&nElectrons==0&&nMuons==0&&nIsoTracks15_005_03==0",kRed,"ttbar 0l");
  addSample("TTbarJets:ttbarDecayCode!=2&&(nElectrons+nMuons==1)&&MT_Wlep>=0&&MT_Wlep<=100",kBlue,"ttbar 1l");
  selection_=baseline && cleaning && ht && btag && met && mdp;
  var="MT_b"; xtitle="MT(b,MET)";
  nbins=50; low=0;high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ttbar_mtb-0l-1l",0);


  selection_=baseline && cleaning && ht && btag  && mdp&&TCut("MT_b>=200");
  var="MET"; xtitle="MET";
  setLogY(true);
  nbins=30; low=125;high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ttbar_met_mtb200-0l-1l",0);


  selection_=baseline && cleaning && met && btag  && mdp&&TCut("MT_b>=200");
  var="HT"; xtitle="HT";
  setLogY(true);
  nbins=30; low=400;high=1200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ttbar_ht_mtb200-0l-1l",0);


  //resetChains();
  setLogY(false);
  clearSamples();
  addSample("TTbarJets");
  addSample("T2tt$600$0",kRed);
 
  selection_=baseline && cleaning && btag && ht && mdp && zl&&TCut("MT_b>=200");
  var="MET"; xtitle="MET";
  nbins=30; low=125;high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_met_mtb200",0);


}


void rocstudyT2tt() {


  inputPath = "/cu4/ra2b/reducedTrees/v66_8/"; //change to 8 for 'new' MT vari

  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  loadSamples(true,"ra2b2012");
  usePUweight_=false;

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();


  //use a different sample order than default
  //resetChains();
  clearSamples();

  addSample("TTbarJets");

  addSample("T2tt$600$0",kRed);
  useMassInLegend_=true;

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=false;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40";
  TCut ht = "HT>=400";
  TCut met = "MET>=125";
  TCut mdp = "minDeltaPhiN_asin>=4";
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");
  TCut btag = "nbjets>=1";

  setStackMode(false,true);

  selection_=baseline && cleaning && ht && btag && met && mdp && zl;
  var="MT_b"; xtitle="MT(closest b,MET)";
  nbins=80; low=0;high=1000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_mtb",0);

  drawEffVRej("T2tt$600$0","TTbarJets","T2tt (600,0)","ttbar",true);
  TGraphAsymmErrors* g_mtb = (TGraphAsymmErrors*) effgraph->Clone("g_mtb");

  var="MT_bestCSV"; xtitle="MT(b-like,MET)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_mtcsv",0);

  drawEffVRej("T2tt$600$0","TTbarJets","T2tt (600,0)","ttbar",true);
  TGraphAsymmErrors* g_mtcsv = (TGraphAsymmErrors*) effgraph->Clone("g_mtcsv");

  var="minMT_jetMET"; xtitle="min MT(jets1..4,MET)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_minmtjetmet",0);
  drawEffVRej("T2tt$600$0","TTbarJets","T2tt (600,0)","ttbar",true);
  TGraphAsymmErrors* g_mtjetmet = (TGraphAsymmErrors*) effgraph->Clone("g_mtjetmet");

  var="MT_jim"; xtitle="hybrid MT";
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_t2tt_mtjim",0);

  drawEffVRej("T2tt$600$0","TTbarJets","T2tt (600,0)","ttbar",true);
  TGraphAsymmErrors* g_jim = (TGraphAsymmErrors*) effgraph->Clone("g_jim");

  g_mtb->SetLineColor(kRed);
  g_mtb->SetMarkerColor(kRed);

  g_mtcsv->SetLineColor(kBlue);
  g_mtcsv->SetMarkerColor(kBlue);

  g_mtjetmet->SetLineColor(kMagenta);
  g_mtjetmet->SetMarkerColor(kMagenta);

  g_jim->SetLineColor(kGreen);
  g_jim->SetMarkerColor(kGreen);

  g_mtjetmet->Draw("pl");
  g_mtcsv->Draw("PL");

  g_mtb->Draw("PL");
  g_jim->Draw("PL");


}

void ttbar_MTb_tail() {

  inputPath = "/cu4/ra2b/reducedTrees/v66_8/"; //change to 8 for 'new' MT vari

  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  loadSamples(true,"ra2b2012");
  usePUweight_=true;

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();


  //use a different sample order than default
  //resetChains();



  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=false;

  doOverflowAddition(true);
  doRatio_=false; ratioMin = 0; ratioMax = 2.2;
  dodata_=false;

  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40";
  TCut ht = "HT>=400";
  TCut met = "MET>=125";
  TCut mdp = "minDeltaPhiN_asin>=4";
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");
  TCut btag = "nbjets>=1";

  // -- as a function of my usual W decay categories

  clearSamples();
  addSample("TTbarJets0:ttbarDecayCode==2",kYellow,"all had");
  addSample("TTbarJets0:ttbarDecayCode==7||ttbarDecayCode==5||ttbarDecayCode==6",kRed,"W #rightarrow #tau");
  addSample("TTbarJets0:ttbarDecayCode==3||ttbarDecayCode==4",kMagenta,"W #rightarrow e / #mu");
  addSample("TTbarJets0:ttbarDecayCode==1||ttbarDecayCode==9||ttbarDecayCode==8",kOrange+7,"dilepton (incl #tau)");

  //baseline ZL selection
  selection_=baseline && cleaning && ht && btag && met && mdp && zl;
  var="MT_bestCSV"; xtitle="MT(b-like,MET)";
  nbins=50; low=0;high=500;

  //normalized
  setStackMode(false,true,false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ttbar_mtb_leptonicsplit",0);
  //stack
  setStackMode(true,false,false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ttbar_mtb_leptonicsplit",0);

  // -- as a function of number of b-tags
  clearSamples();
  addSample("TTbarJets0:nbjets==1",kBlue+4,"1 b");
  addSample("TTbarJets0:nbjets==2",kBlue,"2 b");
  addSample("TTbarJets0:nbjets>=3",kBlue-9,">=3 b");

 nbins=40;

  //normalized
  setStackMode(false,true,false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ttbar_mtb_bsplit",0);
  //stack
  setStackMode(true,false,false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ttbar_mtb_bsplit",0);


  // -- as a function of the new categories
  clearSamples();
  addSample("TTbarJets0:ttbarDecayCode==2",kYellow,"all had");
  addSample("TTbarJets0:ttbarDecayCode>=3&&ttbarDecayCode<=7&&MT_bestCSV_gencode==1",kGreen+3,"b in right top");
  addSample("TTbarJets0:ttbarDecayCode>=3&&ttbarDecayCode<=7&&MT_bestCSV_gencode==2",kRed,"b in wrong top");
  addSample("TTbarJets0:ttbarDecayCode>=3&&ttbarDecayCode<=7&&MT_bestCSV_gencode>=3",kRed-6,"other");
  addSample("TTbarJets0:ttbarDecayCode==1||ttbarDecayCode==9||ttbarDecayCode==8",kOrange+7,"dileptonic");

  //normalized
  setStackMode(false,true,false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ttbar_mtb_matchsplit",0);
  //stack
  setStackMode(true,false,false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ttbar_mtb_matchsplit",0);


}


void dataABCD() {

  //updated to compare final 1.7/fb to the rest
  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  //must be done before loadSamples()
  addToSamplesAll("HTRun2012A");
  addToSamplesAll("HTMHTRun2012B");
  addToSamplesAll("HTMHTRun2012C1");
  addToSamplesAll("HTMHTRun2012C2");
  addToSamplesAll("HTMHTRun2012D");
  addToSamplesAll("HTMHTRun2012D2");

  addToSamplesAll("JetHTRun2012B");
  addToSamplesAll("JetHTRun2012C1");
  addToSamplesAll("JetHTRun2012C2");
  addToSamplesAll("JetHTRun2012D");
  addToSamplesAll("JetHTRun2012D2");

  addToSamplesAll("METRun2012A");
  addToSamplesAll("METRun2012B");
  addToSamplesAll("METRun2012C1");
  addToSamplesAll("METRun2012C2");
  addToSamplesAll("METRun2012D");
  addToSamplesAll("METRun2012D2");

  loadSamples(true,"ra2b2012");

  usePUweight_=false;
  dodata_=false;

  clearSamples();
  addSample("HTRun2012A",kBlue,"2012 A+B+C+D1");
  chainSamples("HTRun2012A","HTMHTRun2012B");
  chainSamples("HTRun2012A","HTMHTRun2012C1");
  chainSamples("HTRun2012A","HTMHTRun2012C2");
  chainSamples("HTRun2012A","HTMHTRun2012D");

  chainSamples("HTRun2012A","JetHTRun2012B");
  chainSamples("HTRun2012A","JetHTRun2012C1");
  chainSamples("HTRun2012A","JetHTRun2012C2");
  chainSamples("HTRun2012A","JetHTRun2012D");

  chainSamples("HTRun2012A","METRun2012A");
  chainSamples("HTRun2012A","METRun2012B");
  chainSamples("HTRun2012A","METRun2012C1");
  chainSamples("HTRun2012A","METRun2012C2");
  chainSamples("HTRun2012A","METRun2012D");

  //  setSampleScaleFactor("HTRun2012A",5580./12034.14);
  setSampleScaleFactor("HTRun2012A",1700./(5580.+12034.14)); //1700 is approx

//   addSample("METRun2012D",kRed,"2012 D");
//   chainSamples("METRun2012D","HTMHTRun2012D");
//   chainSamples("METRun2012D","JetHTRun2012D");

  addSample("METRun2012D2",kRed,"2012 D2");
  chainSamples("METRun2012D2","HTMHTRun2012D2");
  chainSamples("METRun2012D2","JetHTRun2012D2");

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
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ABCDvD2_njets"+sampledesc+desc,0);

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

  useTrigEff_=false;
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


  useTrigEff_=false;
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


  useTrigEff_=false;
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

  useTrigEff_=false;
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

void RA2b_ARC_round1() {


  inputPath = "/cu4/ra2b/reducedTrees/v66_8/"; //switch to v66_8

  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  loadSamples(true,"ra2b2012");
  usePUweight_=true;

  useTrigEff_=false;
  currentConfig_=configDescriptions_.getDefault();

  resetChains();
  clearSamples();
  
  addSample("TTbarJets:ttbarDecayCode==2",kYellow,"all hadronic");
  //  addSample("TTbarJets:ttbarDecayCode==7",kRed,"#tau #rightarrow h + W #rightarrow qq");
  addSample("TTbarJets:ttbarDecayCode==3||ttbarDecayCode==4||ttbarDecayCode==5||ttbarDecayCode==6||ttbarDecayCode==7",kMagenta,"1 lepton (including #tau)");
  addSample("TTbarJets:ttbarDecayCode==1||ttbarDecayCode==9||ttbarDecayCode==8",kOrange,"dilepton");


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
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");

  btagSFweight_="probge1"; //no btag req

  setStackMode(true);

  selection_ = baseline&&cleaning && met &&ht&& mdp && zl;
  var="MET"; xtitle="E^{T}_{miss} [GeV]";
  nbins = 40; low=125; high=325;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ARC_ttbar_met",0);


  clearSamples();
  setStackMode(false,true);
  addSample("TTbarJets:ttbarDecayCode==2",kRed,"hadronic t#bar{t}");
  addSample("TTbarJets:ttbarDecayCode==3||ttbarDecayCode==4||ttbarDecayCode==5||ttbarDecayCode==6||ttbarDecayCode==7",kMagenta,"1 lepton t#bar{t}");
  addSample("PythiaPUQCD",kBlue,"QCD");

  selection_ = baseline&&cleaning && met &&ht && zl;
  var="minDeltaPhiN_asin"; xtitle="#Delta#hat{#phi}_{min}";
  nbins = 40; low=0; high=40;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ARC_ttbar_mindpn",0);

  //check contribution of "other samples"
  clearSamples();
  setStackMode(true,false);
  addSample("VV");
  addSample("ZJets");
  selection_ = baseline&&cleaning && met &&ht&& mdp && zl;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ARC_misc_met",0);
  btagSFweight_ = "probge2";
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ARC_misc_met_ge2b",0);

}

void RA2b_ttbar_njets() {

  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  loadSamples(true,"ra2b2012");
  usePUweight_=true;

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

  setStackMode(false,false);

  TCut base = "HT>=400 &&HT<500&& MET>=125&& cutPV==1 && njets>=3 &&jetpt2>=70 && minDeltaPhiN_asin >= 4 && passCleaning==1 &&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40 && ttbarDecayCode!=2";

  selection_ =base&&TCut("nbjets==2");
  ratioMin=0.4;
  ratioMax=1.6;

  var="njets"; xtitle="jet multiplicity";
  nbins = 4; low=3; high=7;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_njetsshapes_0l_1l_nb2");

  TH1D* ratio_2b_mg = (TH1D*)ratio->Clone("ratio_2b_mg");

  selection_ =base&&TCut("nbjets>=3");
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_njetsshapes_0l_1l_nb3");

  TH1D* ratio_3b_mg = (TH1D*)ratio->Clone("ratio_3b_mg");

  //-- powheg --

  clearSamples();
  addSample("TTbarJetsPowheg:nElectrons==0&&nMuons==0&&nIsoTracks15_005_03==0",kRed,"0 lep");
  addSample("TTbarJetsPowheg:(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)",kBlue,"1 lep");
   
  chainSamples("TTbarJetsPowheg","WJets");
  chainSamples("TTbarJetsPowheg","SingleTop");

  selection_ =base&&TCut("nbjets==2");
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_njetsshapes_0l_1l_nb2_powheg");
  TH1D* ratio_2b_ph = (TH1D*)ratio->Clone("ratio_2b_ph");

  selection_ =base&&TCut("nbjets>=3");
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_njetsshapes_0l_1l_nb3_powheg");
  TH1D* ratio_3b_ph = (TH1D*)ratio->Clone("ratio_3b_ph");

  //-mc@nlo

  clearSamples();
  addSample("TTbarJetsMCNLO:nElectrons==0&&nMuons==0&&nIsoTracks15_005_03==0",kRed,"0 lep");
  addSample("TTbarJetsMCNLO:(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)",kBlue,"1 lep");
   
  chainSamples("TTbarJetsMCNLO","WJets");
  chainSamples("TTbarJetsMCNLO","SingleTop");
  selection_ =base&&TCut("nbjets==2");
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_njetsshapes_0l_1l_nb2_mcnlo");
  TH1D* ratio_2b_nlo = (TH1D*)ratio->Clone("ratio_2b_nlo");

  selection_ =base&&TCut("nbjets>=3");
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_njetshapes_0l_1l_nb3_mcnlo");
  TH1D* ratio_3b_nlo = (TH1D*)ratio->Clone("ratio_3b_nlo");


  ratio_2b_mg->SetMarkerColor(kRed);
  ratio_3b_mg->SetMarkerColor(kRed-7);

  ratio_2b_ph->SetMarkerColor(kBlue);
  ratio_3b_ph->SetMarkerColor(kBlue-7);

  ratio_2b_nlo->SetMarkerColor(kGreen+4);
  ratio_3b_nlo->SetMarkerColor(kGreen-6);

  ratio_2b_mg->SetLineColor(kRed);
  ratio_3b_mg->SetLineColor(kRed-7);

  ratio_2b_ph->SetLineColor(kBlue);
  ratio_3b_ph->SetLineColor(kBlue-7);

  ratio_2b_nlo->SetLineColor(kGreen+4);
  ratio_3b_nlo->SetLineColor(kGreen-6);

  setLogY(false);
  renewCanvas();

  ratio_2b_mg->GetYaxis()->SetLabelSize(0.05); //make y label bigger
  ratio_2b_mg->GetXaxis()->SetLabelSize(0.05); //make x label bigger

  ratio_2b_mg->SetMaximum(2.5);
  ratio_2b_mg->SetMinimum(0.5);

  ratio_2b_mg->GetYaxis()->SetTitle("0l / 1l");
  ratio_2b_mg->GetYaxis()->SetTitleSize(0.07);

  ratio_2b_mg->Draw();
  ratio_3b_mg->Draw("same");

  ratio_2b_ph->Draw("same");
  ratio_3b_ph->Draw("same");

  ratio_2b_nlo->Draw("same");
  ratio_3b_nlo->Draw("same");

  leg_x1 = 0.2; leg_x2=0.5;
  leg_y1=0.5; leg_y2=0.8;
  renewLegend();
  leg->AddEntry(ratio_2b_mg,"2b MG");
  leg->AddEntry(ratio_3b_mg,"3b MG");
  leg->AddEntry(ratio_2b_ph,"2b Powheg");
  leg->AddEntry(ratio_3b_ph,"3b Powheg");
  leg->AddEntry(ratio_2b_nlo,"2b MC@NLO");
  leg->AddEntry(ratio_3b_nlo,"3b MC@NLO");
  leg->Draw();

}

void RA2b_ttbarstudies() {

  //goal: compare HT (and MET and nb) shapes for different 0-lepton sample versus 1 lepton

  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  addToSamplesAll("TTbarJetsMatchingDown");
  addToSamplesAll("TTbarJetsMatchingUp");
  addToSamplesAll("TTbarJetsScaleUp");

  loadSamples(true,"ra2b2012");
  usePUweight_=true;

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

  /*
  resetChains();
  clearSamples();
  addSample("TTbarJetsMatchingUp:nElectrons==0&&nMuons==0&&nIsoTracks15_005_03==0",kRed,"0 lep");
  addSample("TTbarJetsMatchingUp:(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)",kBlue,"1 lep");
  chainSamples("TTbarJetsMatchingUp","WJets");
  chainSamples("TTbarJetsMatchingUp","SingleTop");

  resetChains();
  clearSamples();
  addSample("TTbarJetsMatchingDown:nElectrons==0&&nMuons==0&&nIsoTracks15_005_03==0",kRed,"0 lep");
  addSample("TTbarJetsMatchingDown:(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)",kBlue,"1 lep");
  chainSamples("TTbarJetsMatchingDown","WJets");
  chainSamples("TTbarJetsMatchingDown","SingleTop");

  resetChains();
  clearSamples();
  addSample("TTbarJetsScaleUp:nElectrons==0&&nMuons==0&&nIsoTracks15_005_03==0",kRed,"0 lep");
  addSample("TTbarJetsScaleUp:(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)",kBlue,"1 lep");
  chainSamples("TTbarJetsScaleUp","WJets");
  chainSamples("TTbarJetsScaleUp","SingleTop");

  */

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



  //--compare 2b and 3b

  /* large block of code that was specific to a study of 2b versus 3b with different generators

  setStackMode(false,false);
  selection_ =TCut("HT>=400 &&MET>=125&& cutPV==1 && njets>=3 &&jetpt2>=70 && minDeltaPhiN_asin >= 4 && passCleaning==1 && nbjets==2 &&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40 && ttbarDecayCode!=2");
  ratioMin=0.4;
  ratioMax=1.6;

  float nhtbins=4;
  float htbins[] = {400,500,800,1000,1200};
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = nhtbins; low=400; high=1200;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_htshapes_0l_1l_nb2",htbins,"GeV");

  TH1D* ratio_2b_mg = (TH1D*)ratio->Clone("ratio_2b_mg");

  selection_ =TCut("HT>=400 &&MET>=125&& cutPV==1 && njets>=3 &&jetpt2>=70 && minDeltaPhiN_asin >= 4 && passCleaning==1 && nbjets>=3 &&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40 && ttbarDecayCode!=2");
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_htshapes_0l_1l_nb3",htbins,"GeV");

  TH1D* ratio_3b_mg = (TH1D*)ratio->Clone("ratio_3b_mg");

  //-- powheg --

  clearSamples();
  addSample("TTbarJetsPowheg:nElectrons==0&&nMuons==0&&nIsoTracks15_005_03==0",kRed,"0 lep");
  addSample("TTbarJetsPowheg:(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)",kBlue,"1 lep");
   
  chainSamples("TTbarJetsPowheg","WJets");
  chainSamples("TTbarJetsPowheg","SingleTop");
  selection_ =TCut("HT>=400 &&MET>=125&& cutPV==1 && njets>=3 &&jetpt2>=70 && minDeltaPhiN_asin >= 4 && passCleaning==1 && nbjets==2 &&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40 && ttbarDecayCode!=2");
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_htshapes_0l_1l_nb2_powheg",htbins,"GeV");
  TH1D* ratio_2b_ph = (TH1D*)ratio->Clone("ratio_2b_ph");

  selection_ =TCut("HT>=400 &&MET>=125&& cutPV==1 && njets>=3 &&jetpt2>=70 && minDeltaPhiN_asin >= 4 && passCleaning==1 && nbjets>=3 &&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40 && ttbarDecayCode!=2");
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_htshapes_0l_1l_nb3_powheg",htbins,"GeV");
  TH1D* ratio_3b_ph = (TH1D*)ratio->Clone("ratio_3b_ph");

  //-mc@nlo

  clearSamples();
  addSample("TTbarJetsMCNLO:nElectrons==0&&nMuons==0&&nIsoTracks15_005_03==0",kRed,"0 lep");
  addSample("TTbarJetsMCNLO:(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)",kBlue,"1 lep");
   
  chainSamples("TTbarJetsMCNLO","WJets");
  chainSamples("TTbarJetsMCNLO","SingleTop");
  selection_ =TCut("HT>=400 &&MET>=125&& cutPV==1 && njets>=3 &&jetpt2>=70 && minDeltaPhiN_asin >= 4 && passCleaning==1 && nbjets==2 &&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40 && ttbarDecayCode!=2");
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_htshapes_0l_1l_nb2_mcnlo",htbins,"GeV");
  TH1D* ratio_2b_nlo = (TH1D*)ratio->Clone("ratio_2b_nlo");

  selection_ =TCut("HT>=400 &&MET>=125&& cutPV==1 && njets>=3 &&jetpt2>=70 && minDeltaPhiN_asin >= 4 && passCleaning==1 && nbjets>=3 &&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40 && ttbarDecayCode!=2");
  drawPlots(var,nbins,low,high,xtitle,"au", "ra2b_htshapes_0l_1l_nb3_mcnlo",htbins,"GeV");
  TH1D* ratio_3b_nlo = (TH1D*)ratio->Clone("ratio_3b_nlo");


  ratio_2b_mg->SetMarkerColor(kRed);
  ratio_3b_mg->SetMarkerColor(kRed-7);

  ratio_2b_ph->SetMarkerColor(kBlue);
  ratio_3b_ph->SetMarkerColor(kBlue-7);

  ratio_2b_nlo->SetMarkerColor(kGreen+4);
  ratio_3b_nlo->SetMarkerColor(kGreen-6);

  ratio_2b_mg->SetLineColor(kRed);
  ratio_3b_mg->SetLineColor(kRed-7);

  ratio_2b_ph->SetLineColor(kBlue);
  ratio_3b_ph->SetLineColor(kBlue-7);

  ratio_2b_nlo->SetLineColor(kGreen+4);
  ratio_3b_nlo->SetLineColor(kGreen-6);

  setLogY(false);
  renewCanvas();

  ratio_2b_mg->GetYaxis()->SetLabelSize(0.05); //make y label bigger
  ratio_2b_mg->GetXaxis()->SetLabelSize(0.05); //make x label bigger

  ratio_2b_mg->SetMaximum(2.1);
  ratio_2b_mg->SetMinimum(0.6);

  ratio_2b_mg->GetYaxis()->SetTitle("0l / 1l");
  ratio_2b_mg->GetYaxis()->SetTitleSize(0.07);

  ratio_2b_mg->Draw();
  ratio_3b_mg->Draw("same");

  ratio_2b_ph->Draw("same");
  ratio_3b_ph->Draw("same");

  ratio_2b_nlo->Draw("same");
  ratio_3b_nlo->Draw("same");

  leg_x1 = 0.2; leg_x2=0.5;
  leg_y1=0.5; leg_y2=0.8;
  renewLegend();
  leg->AddEntry(ratio_2b_mg,"2b MG");
  leg->AddEntry(ratio_3b_mg,"3b MG");
  leg->AddEntry(ratio_2b_ph,"2b Powheg");
  leg->AddEntry(ratio_3b_ph,"3b Powheg");
  leg->AddEntry(ratio_2b_nlo,"2b MC@NLO");
  leg->AddEntry(ratio_3b_nlo,"3b MC@NLO");
  leg->Draw();

  */

  //exercise where i was comparing 0l and 1l shapes for misc stuff; not relevant here
  /*
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
  */

  //here is the original study; now updating with isotk veto IN PROGRESS

  clearSamples();
  addSample("TTbarJets:nElectrons==0&&nMuons==0&&nIsoTracks15_005_03==0&&ttbarDecayCode==7",kRed,"#tau #rightarrow h");
  addSample("TTbarJets:nElectrons==0&&nMuons==0&&nIsoTracks15_005_03==0&&(ttbarDecayCode==3||ttbarDecayCode==4)",kMagenta,"lost lepton");
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

  //Dec 2012 -- update to modern cuts

  clearSamples();
  addSample("TTbarJets:nElectrons==0&&nMuons==0",kRed,"0 lep");
  addSample("TTbarJets:(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)",kBlue,"1 lep");

  selection_ =TCut("HT>=400 &&MET>=125&& cutPV==1 && njets>=3 &&jetpt2>=70 && minDeltaPhiN_asin >= 4 && passCleaning==1 && nbjets>=1 &&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40 && njets==3");
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 40; low=400; high=1200;
  setLogY(true);
  doRatioPlot(true);
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
  useTrigEff_=false;
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


  useTrigEff_=false;
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

void RA2b_investigation_RazorNoise() {

  //jan 2013 -- look into razor noise variable
  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples
  loadSamples(true,"ra2b2012");
  usePUweight_=true;

  // "owen factors"
  setSampleScaleFactor("TTbarJets",0.9);
  setSampleScaleFactor("WJets",0.9);
  setSampleScaleFactor("PythiaPUQCD",1.8);
  int nbins;
  float low,high;
  TString var,xtitle;
  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2;
  dodata_=true;
  setDatasetToDraw("2012hybridplus");

  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40";
  TCut ht = "HT>=400";
  TCut met = "MET>=125";
  TCut mdp = "minDeltaPhiN_asin>=4";
  TCut ldp = "minDeltaPhiN_asin<4";
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");
  TCut sl = getSingleLeptonCut();

  //special
  TCut tightht = "HT>=800";
  TCut lightcleaning = "passCleaning==1&&buggyEvent==0";
  TCut tobtec = "maxTOBTECjetDeltaMult<40";
  TCut invertedcleaning = "MET/caloMET>=2 || maxTOBTECjetDeltaMult>=40";

  selection_= baseline && cleaning && ht && met && zl &&mdp;
  //NO b-tagging
  var="min(abs(METphi-caloMETphi),abs(2*3.14159-METphi+caloMETphi))"; xtitle="|#Delta #phi_{PF,calo}|";
  nbins = 30; low=0; high=2*TMath::Pi();
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_baseline_ge0b_deltaPhiMetPfCalo",0);
 
  selection_= baseline && lightcleaning && ht && met && zl &&mdp;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_baseline-noMetRatio-noTobtec_ge0b_deltaPhiMetPfCalo",0);

  selection_= baseline && lightcleaning && ht && met && zl&&mdp;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_baseline-noMetRatio-noTobtec_ge0b_deltaPhiMetPfCalo",0);

  selection_= baseline && lightcleaning&&invertedcleaning && ht && met && zl&&mdp;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_baseline-invertedMetRatio-invertedTobtec_ge0b_deltaPhiMetPfCalo",0);

  selection_= TCut("cutPV==1 && cutTrigger2==1 && njets==2 && jetpt2>=70") && lightcleaning && ht && met && zl&&mdp;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_dijet-noMetRatio-noTobtec_ge0b_deltaPhiMetPfCalo",0);

  //hmmm...not sure what to conclude
  TCut badrazor = "abs(min(abs(METphi-caloMETphi),abs(2*3.14159-METphi+caloMETphi))-3.14159)<1";
  TCut goodrazor = "abs(min(abs(METphi-caloMETphi),abs(2*3.14159-METphi+caloMETphi))-3.14159)>=1";

  var="MET/caloMET"; xtitle="Type I PFMET / caloMET";
  nbins = 40; low=0; high=4;
  setLogY(true);

  selection_= baseline && lightcleaning && ht && met && zl &&badrazor&&mdp;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_baseline_badrazor_ge0b_",0);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_baseline_badrazor_ge0b_",0);

  selection_= baseline && lightcleaning && ht && met && zl &&goodrazor&&mdp;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_baseline_goodrazor_ge0b_",0);
 
  //remake these plots with the tobtec veto

  setLogY(true);
  selection_= baseline && lightcleaning && ht && met && zl &&badrazor &&tobtec&&mdp;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_baseline-withtobtec_badrazor_ge0b_",0);
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_baseline-withtobtec_badrazor_ge0b_",0);

  selection_= baseline && lightcleaning && ht && met && zl &&goodrazor &&tobtec&&mdp;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_baseline-withtobtec_goodrazor_ge0b_",0);
 
  //our actual ZL selection
  selection_= baseline && cleaning && ht && met && zl&&mdp;
  btagSFweight_="probge1";
  var="abs(min(abs(METphi-caloMETphi),abs(2*3.14159-METphi+caloMETphi))-3.14159)"; xtitle="| |#Delta #phi_{PF,calo}| - #pi |";
  nbins = 30; low=0; high=TMath::Pi();
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_baseline_ge1b_deltaPhiMetPfCalo-pi",0);
 
  //try the events we reject with our special cuts
  selection_= baseline && lightcleaning && ht && met && zl &&invertedcleaning&&mdp;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_baseline-invertedMetRatio-invertedTobtec_ge1b_deltaPhiMetPfCalo-pi",0);

  //a more sensitive ZL selection
  selection_= baseline && cleaning && TCut("HT>800") && TCut("MET>200") && zl&&mdp;
  btagSFweight_="probge1";
  var="abs(min(abs(METphi-caloMETphi),abs(2*3.14159-METphi+caloMETphi))-3.14159)"; xtitle="| |#Delta #phi_{PF,calo}| - #pi |";
  nbins = 30; low=0; high=TMath::Pi();
  ratioMax=3;
  setPlotMinimum(0.1);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht800met200_ge1b_deltaPhiMetPfCalo-pi",0);
 
  //same selection except look at the LDP
  selection_= baseline && cleaning && TCut("HT>800") && TCut("MET>200") && zl&&ldp;
  btagSFweight_="probge1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ht800met200_ge1b_LDP_deltaPhiMetPfCalo-pi",0);

}


void RA2b_investigation_noise() {

  //use the trees from keith. i don't trust mine for data
  inputPath = "/cu4/ra2b/reducedTrees/v66_2/"; //2012
  dataInputPath =  "/cu4/ra2b/reducedTrees/v66_2/ORIGINALS/";

  nameOfEventWeight_="weight2"; //for 2012 cfA skimmed ntuples
  loadSamples(true,"ra2b2012");
  usePUweight_=true;


  useTrigEff_=false;
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

  //updating using full A+B+C+D

  nameOfEventWeight_="weight3"; //for 2012 cfA skimmed ntuples

  loadSamples(true,"ra2b2012");
  usePUweight_=true;

  useTrigEff_=true;
  currentConfig_=configDescriptions_.getDefault();

  //use a different sample order than default
  resetChains();
  clearSamples();

  addSample("VV");
  addSample("TTV",kPink,"TTV");
  addSample("ZJets");

  addSample("Zinvisible");
  addSample("PythiaPUQCD");
  addSample("SingleTop");
  addSample("WJets");
  addSample("TTbarJets");

  //  addSample("T1bbbb$1200$200",kRed);
  //  useMassInLegend_=true;

  //use the 'owen factors' verbatim
  // for Fig 3 of the PAS we would need a factor from the fit
  //  setSampleScaleFactor("TTbarJets",0.9);
  //  setSampleScaleFactor("WJets",0.9);
  //  setSampleScaleFactor("PythiaPUQCD",1.8);

  int nbins;
  float low,high;
  TString var,xtitle;

  drawFilenameOnPlot_=false;
  drawMCErrors_=true;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  TCut baseline = "cutPV==1 && cutTrigger2==1 && njets>=3 && jetpt2>=70";
  TCut cleaning="passCleaning==1&&buggyEvent==0&&MET/caloMET<2&&maxTOBTECjetDeltaMult<40";
  TCut ht = "HT>=400";
  TCut met = "MET>=125";
  TCut mdp = "minDeltaPhiN_asin>=4";
  TCut ldp = "minDeltaPhiN_asin<4";
  TCut zl = getLeptonVetoCut()&&TCut("nIsoTracks15_005_03==0");

  setDatasetToDraw("2012hybridplus");

  //ht 2,3,4 and 2-4
  btagSFweight_="prob2";
  TCut ht1 = "HT>=400&&HT<500";
  TCut ht2 = "HT>=500&&HT<800";
  TCut ht3 = "HT>=800&&HT<1000";
  TCut ht4 = "HT>=1000";
  TCut ht234= "HT>500";
  TCut ht1234= "HT>=400";
  TCut met2="MET>=150&&MET<250";
  TCut met3="MET>=250&&MET<350";
  TCut met4="MET>=350";
  TCut met234="MET>=150";

  //  savePlots_=false;
  //  quiet_=true;
  var="MET"; xtitle="E^{T}_{miss} [GeV]";
  nbins = 30; low=125; high=500;

  if (false) {
  //MET 4, HT2
  selection_ = baseline&&cleaning && met4 && mdp && zl && ht2;
  drawPlots(var,nbins,low,high,xtitle,"Events", "blah",0,"GeV");

  selection_ = baseline&&cleaning && met4 && mdp && zl && ht3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "blah",0,"GeV");

  selection_ = baseline&&cleaning && met4 && mdp && zl && ht4;
  drawPlots(var,nbins,low,high,xtitle,"Events", "blah",0,"GeV");

  selection_ = baseline&&cleaning && met4 && mdp && zl && ht234;
  drawPlots(var,nbins,low,high,xtitle,"Events", "blah",0,"GeV");

  btagSFweight_="probge3";
  for (int iht = 0;iht<5; iht++) {
    for (int imet=0;imet<4;imet++) {

      TCut thisht,thismet;
      if (iht==0) thisht=ht1;
      if (iht==1) thisht=ht2;
      if (iht==2) thisht=ht3;
      if (iht==3) thisht=ht4;
      if (iht==4) thisht=ht1234;

      if (imet==0) thismet=met2;
      if (imet==1) thismet=met3;
      if (imet==2) thismet=met4;
      if (imet==3) thismet=met234;

      selection_ = baseline&&cleaning && thismet && mdp && zl && thisht;
      drawPlots(var,nbins,low,high,xtitle,"Events", "blah",0,"GeV");

    }
  }
  }
  btagSFweight_="probge3";
  setLogY(true);

  selection_ = baseline&&cleaning && met && ldp && zl && ht;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ldp_ge3b",0,"GeV");

  btagSFweight_="";
  selection_ = baseline&&cleaning && met && ldp && zl && ht&&TCut("nbjets>=3");
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ldp_ge3b-nobtagsf",0,"GeV");

  btagSFweight_="";
  selection_ = baseline&&cleaning && met && ldp && zl && ht;
  var="nbjets"; xtitle="n bjets";
  nbins = 6; low=0; high=6;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_ldp_nbjets",0);

  savePlots_=false;
  var="HT"; xtitle="ht";
  nbins = 1; low=0; high=9999;

  selection_ = baseline&&cleaning && met && ldp && zl && ht;

  btagSFweight_="prob0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0,"GeV");
  TH1D* bin0=(TH1D*) totalsm->Clone("bin0");
  TH1D* bin0data=(TH1D*) hdata->Clone("bin0data");


  btagSFweight_="prob1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0,"GeV");
  TH1D* bin1=(TH1D*) totalsm->Clone("bin1");
  TH1D* bin1data=(TH1D*) hdata->Clone("bin1data");

  btagSFweight_="prob2";
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0,"GeV");
  TH1D* bin2=(TH1D*) totalsm->Clone("bin2");
  TH1D* bin2data=(TH1D*) hdata->Clone("bin2data");

  btagSFweight_="probge3";
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0,"GeV");
  TH1D* bin3=(TH1D*) totalsm->Clone("bin3");
  TH1D* bin3data=(TH1D*) hdata->Clone("bin3data");



  return;


  //  const  TString btagselection = "probge3";
  const  TString btagselection = "probge3";

  //updating to imitate PAS Fig 3
  btagSFweight_=""; //no btag req

  selection_ = baseline&&cleaning && met && mdp && zl && ht;
  var="nbjets"; xtitle="Number of b-tagged jets";
  nbins = 5; low=0.5; high=5.5;
  //   drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_nbjets_"+desc,0);
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_PAS_nbjets",0);
//   setLogY(false);

  //add vertical lines
  addVerticalLine(150);
  addVerticalLine(250);
  addVerticalLine(350);

  btagSFweight_=btagselection;
  selection_ = baseline&&cleaning && mdp && zl &&ht&&met;
  var="MET"; xtitle="E^{T}_{miss} [GeV]";
  nbins = 30; low=125; high=500;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_PAS_met_ge3b",0,"GeV");
  useTrigEff_=true;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ra2b_PAS_met_ge3b-trigeff",0,"GeV");


//old stuff down here

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

  useTrigEff_=false;
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

  useTrigEff_=false;
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

  useTrigEff_=false;
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

  useTrigEff_=false;
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

  useTrigEff_=false;
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


