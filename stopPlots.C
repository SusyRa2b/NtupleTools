/*
to compile:

must have symlink to MiscUtil.cxx in the working directory.

gSystem->Load("TSelectorMultiDraw_C.so");
gSystem->Load("CrossSectionTable_cxx.so");
gSystem->Load("ConfigurationDescriptions_cxx.so");
//gSystem->Load("SearchRegion_cxx.so");
//gSystem->Load("SystInfo_cxx.so");
//gSystem->Load("SignalEffData_cxx.so");
.L stopPlots.C+

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
	
//these are redefined in the functions below
TString inputPath = "/cu2/ra2b/reducedTrees/v63_0/"; //2012
TString dataInputPath =  "/cu2/ra2b/reducedTrees/v63_0/ORIGINALS/";

double lumiScale_ = 5000;//dummy number
double preLumiScale_ = 30;//god only knows for 2012

//make a symlink that point from this name to drawReducedTree.h
//this is to make the ROOT dictionary generation work correctly
#include "stopPlots.h"

void prettyplots() { //not important except maybe not sanity checks

  inputPath = "/cu3/joshmt/stop/embedTrees/v2/";
  reducedTreeName_ = "maketree/reducedTree"; 
  loadSamples(true,"stop");
  clearSamples();
  addSample("TTbarJetsPowheg:ttbarDecayCode==4||ttbarDecayCode==6",kAzure-3,"#mu, #tau #rightarrow #mu");
  addSample("TTbarJetsPowheg:ttbarDecayCode==3||ttbarDecayCode==5",kRed,"e, #tau #rightarrow e");
  addSample("TTbarJetsPowheg:ttbarDecayCode==1||ttbarDecayCode==8||ttbarDecayCode==9",kBlue,"dileptonic");
  addSample("TTbarJetsPowheg:ttbarDecayCode==2",kYellow,"hadronic");
  addSample("TTbarJetsPowheg:ttbarDecayCode==7",kOrange,"#tau #rightarrow h");
  currentConfig_=configDescriptions_.getDefault();
  dodata_=false;

  TCut baseline = "MET>175 && njets70>=2 && njets50>=4 &&njets30>=5";
  TCut lepveto = "electronpt1<5 && muonpt1<5";
  TCut tkveto="isotrackpt1<10";

  //these trees have no weight!
  nameOfEventWeight_="(1)";
  lumiScale_=20000;
  setSampleScaleFactor("TTbarJetsPowheg",234.0/double(getTree("TTbarJetsPowheg")->GetEntries())); //assumes no skim

  int nbins;
  float low,high;
  TString var,xtitle;

  selection_=baseline;

  var="MET"; xtitle=var;
  nbins = 8; low=175; high=350;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbar_met_baseline",0);

  selection_=baseline&&lepveto;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbar_met_lepveto",0);

  selection_=baseline&&lepveto&&tkveto;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbar_met_tkveto",0);

  TCut dp="DeltaPhiMetJet1>0.5 && DeltaPhiMetJet2>0.5 && DeltaPhiMetJet3>0.3";

  selection_=baseline&&lepveto&&tkveto&&dp;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbar_met_dp",0);



}

void acceptance() { 

  inputPath = "/cu3/joshmt/stop/embedTrees/v2/";
  reducedTreeName_ = "maketree/reducedTree"; 
  loadSamples(true,"stop");
  clearSamples();
  addSample("TTbarJetsPowheg:ttbarDecayCode==4",kAzure-3,"MC W #rightarrow #mu");
  currentConfig_=configDescriptions_.getDefault();
  dodata_=false;

  //these trees have no weight!
  nameOfEventWeight_="(1)";
  lumiScale_=20000;
  setSampleScaleFactor("TTbarJetsPowheg",234.0/double(getTree("TTbarJetsPowheg")->GetEntries())); //assumes no skim

  int nbins;
  float low,high;
  TString var,xtitle;

  selection_ ="MET>150 && njets70>=2 && njets50>=4 &&njets30>=5"; //baseline (with looser MET)
  setLogY(false);
  var="MET"; xtitle="MET";
  nbins = 6; low=150; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbar_met_loosecuts",0);
  TH1D* nocuts = (TH1D*) getHist("TTbarJetsPowheg:ttbarDecayCode==4")->Clone("all");

  selection_ = "MET>150 && njets70>=2 && njets50>=4 &&njets30>=5 && genLeptonPt1>=5 && genLeptonEta1>-2.4 && genLeptonEta1<2.4";
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbar_met_inac_loosecuts",0);
  TH1D* inac = (TH1D*) getHist("TTbarJetsPowheg:ttbarDecayCode==4")->Clone("inac2");
  renewCanvas();
  inac->Divide(nocuts);
  inac->SetFillColor(0);
  inac->Draw("he");
}

void efficiency_tauhad() { //18 Dec 2012

  /*
    goal -- make the following histograms:
    TH1 -- eff_track(eta)
    TH1 -- eff_mu|track(eta)
    TH1 -- eff_id,is|mu(pT)
    
    output is in the form of graphs (for human viewing) and TH1D for use in code
    
    These histograms are used for closure tests in MC
  */

  inputPath = "/cu3/joshmt/stop/embedTrees/v4/";
  reducedTreeName_ = "maketree/reducedTree"; 
  loadSamples(true,"stop");
  clearSamples();
  addSample("TTbarJetsPowheg:ttbarDecayCode==4",kAzure-3,"MC W #rightarrow #mu");
  currentConfig_=configDescriptions_.getDefault();
  dodata_=false;

  //these trees have no weight!
  nameOfEventWeight_="(1)";
  lumiScale_=1;
  //setSampleScaleFactor("TTbarJetsPowheg",234.0/double(getTree("TTbarJetsPowheg")->GetEntries())); //assumes no skim

  int nbins;
  float low,high;
  TString var,xtitle;

  savePlots_=false;

  // ------ as a function of eta ------

  //denominator is all W->mu events pT>5, |eta|<2.4
  selection_ = "MET>175 && njets70>=2 && njets50>=4 &&njets30>=5 && genLeptonEta1>=-2.4 && genLeptonEta1<=2.4 && genLeptonPt1>=5"; //baseline
  var="genLeptonEta1"; xtitle="#eta";
  nbins = 40; low=-2.4; high=2.4;
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbar_muoneta_all",0);
  TH1D* muoneta_all = (TH1D*)  hinteractive->Clone("muoneeta_all");

  //require a matched track
  selection_ = "MET>175 && njets70>=2 && njets50>=4 &&njets30>=5 && genLeptonEta1>=-2.4 && genLeptonEta1<=2.4 && genLeptonPt1>=5 && genmatchedTrackPt1>=5";
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbar_muoneta_track",0);
  TH1D* muoneta_track = (TH1D*)  hinteractive->Clone("muoneeta_track");

  //require a track and a muon
  selection_ = "MET>175 && njets70>=2 && njets50>=4 &&njets30>=5 && genLeptonEta1>=-2.4 && genLeptonEta1<=2.4 && genLeptonPt1>=5 && genmatchedTrackPt1>=5 && genmatchedMuPt1>=5";
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbar_muoneta_mu",0);
  TH1D* muoneta_mu = (TH1D*)  hinteractive->Clone("muoneeta_mu");

  //now the reco matched sample
  selection_ = "MET>175 && njets70>=2 && njets50>=4 &&njets30>=5 && genLeptonEta1>=-2.4 && genLeptonEta1<=2.4 && genLeptonPt1>=5 && genmatchedTrackPt1>=5 && genmatchedMuPt1>=5 && muonpt1>=5 && sqrt(acos(cos(genLeptonPhi1-muonphi1))*acos(cos(genLeptonPhi1-muonphi1)) + (genLeptonEta1-muoneta1)*(genLeptonEta1-muoneta1))<0.3";
  //first plot in the same binning and variable as above
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbar_muoneta_recod",0);
  TH1D* muoneta_recod = (TH1D*)  hinteractive->Clone("muoneeta_recod");

  if (effgraph != 0) delete effgraph;
  effgraph = new TGraphAsymmErrors();
  effgraph->BayesDivide(muoneta_track,muoneta_all);
  effgraph->GetHistogram()->SetYTitle("tracking eff");
  renewCanvas();
  effgraph->Draw("AP");
  thecanvas->SaveAs("stop_eff_MC_tracking.eps");
  thecanvas->SaveAs("stop_eff_MC_tracking.png");
  thecanvas->SaveAs("stop_eff_MC_tracking.pdf");

  if (effgraph != 0) delete effgraph;
  effgraph = new TGraphAsymmErrors();
  effgraph->BayesDivide(muoneta_mu,muoneta_track);
  effgraph->GetHistogram()->SetYTitle("muon reco eff");
  renewCanvas();
  effgraph->Draw("AP");
  thecanvas->SaveAs("stop_eff_MC_muon.eps");
  thecanvas->SaveAs("stop_eff_MC_muon.png");
  thecanvas->SaveAs("stop_eff_MC_muon.pdf");

  if (effgraph != 0) delete effgraph;
  effgraph = new TGraphAsymmErrors();
  effgraph->BayesDivide(muoneta_recod  ,muoneta_mu);
  effgraph->GetHistogram()->SetYTitle("id/iso eff");
  renewCanvas();
  effgraph->Draw("AP");
  thecanvas->SaveAs("stop_eff_MC_idiso.eps");
  thecanvas->SaveAs("stop_eff_MC_idiso.png");
  thecanvas->SaveAs("stop_eff_MC_idiso.pdf");


  //in 1d -- deltaR(mu,jet)
  //require a track and a muon
  var="minDRmuonJet"; xtitle="min #Delta R(#mu, jet)";
  nbins = 10; low=0.3; high=6;

  selection_ = "MET>175 && njets70>=2 && njets50>=4 &&njets30>=5 && genLeptonEta1>=-2.4 && genLeptonEta1<=2.4 && genLeptonPt1>=5 && genmatchedTrackPt1>=5 && genmatchedMuPt1>=5";
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbar_muonDeltaR_mu",0);
  TH1D* muonDeltaR_mu = (TH1D*)  hinteractive->Clone("muonDeltaR_mu");

  //now the reco matched sample
  selection_ = "MET>175 && njets70>=2 && njets50>=4 &&njets30>=5 && genLeptonEta1>=-2.4 && genLeptonEta1<=2.4 && genLeptonPt1>=5 && genmatchedTrackPt1>=5 && genmatchedMuPt1>=5 && muonpt1>=5 && sqrt(acos(cos(genLeptonPhi1-muonphi1))*acos(cos(genLeptonPhi1-muonphi1)) + (genLeptonEta1-muoneta1)*(genLeptonEta1-muoneta1))<0.3";
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbar_muonDeltaR_recod",0);
  TH1D* muonDeltaR_recod = (TH1D*)  hinteractive->Clone("muonDeltaR_recod");

  if (effgraph != 0) delete effgraph;
  effgraph = new TGraphAsymmErrors();
  effgraph->BayesDivide(muonDeltaR_recod  ,muonDeltaR_mu);
  effgraph->GetHistogram()->SetYTitle("id/iso eff");
  renewCanvas();
  effgraph->Draw("AP");
  thecanvas->SaveAs("stop_eff_MC_idiso_muonDR.eps");
  thecanvas->SaveAs("stop_eff_MC_idiso_muonDR.png");
  thecanvas->SaveAs("stop_eff_MC_idiso_muonDR.pdf");


  //in 1d -- pT

  var="genLeptonPt1"; xtitle="muon p_{T} [GeV]";
  nbins = 5; low=5; high=200;
  float ptbins[]={5,10,20,40,90,200};
  //require a track and a muon
  selection_ = "MET>175 && njets70>=2 && njets50>=4 &&njets30>=5 && genLeptonEta1>=-2.4 && genLeptonEta1<=2.4 && genLeptonPt1>=5 && genmatchedTrackPt1>=5 && genmatchedMuPt1>=5";
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbar_muonpt_mu",ptbins);
  TH1D* muonpt_mu = (TH1D*)  hinteractive->Clone("muonpt_mu");

  //now the reco matched sample
  selection_ = "MET>175 && njets70>=2 && njets50>=4 &&njets30>=5 && genLeptonEta1>=-2.4 && genLeptonEta1<=2.4 && genLeptonPt1>=5 && genmatchedTrackPt1>=5 && genmatchedMuPt1>=5 && muonpt1>=5 && sqrt(acos(cos(genLeptonPhi1-muonphi1))*acos(cos(genLeptonPhi1-muonphi1)) + (genLeptonEta1-muoneta1)*(genLeptonEta1-muoneta1))<0.3";
  drawPlots(var,nbins,low,high,xtitle,"Events", "ttbar_muonpt_recod",ptbins);
  TH1D* muonpt_recod = (TH1D*)  hinteractive->Clone("muonpt_recod");

  if (effgraph != 0) delete effgraph;
  effgraph = new TGraphAsymmErrors();
  effgraph->BayesDivide(muonpt_recod  ,muonpt_mu);
  effgraph->GetHistogram()->SetYTitle("id/iso eff");
  renewCanvas();
  effgraph->Draw("AP");
  thecanvas->SaveAs("stop_eff_MC_idiso_muonpt.eps");
  thecanvas->SaveAs("stop_eff_MC_idiso_muonpt.png");
  thecanvas->SaveAs("stop_eff_MC_idiso_muonpt.pdf");


  // now save histograms
  TFile fout("stop_eff_MC.root","recreate");

  if (ratio!=0) delete ratio;
  ratio = (TH1D*) muoneta_track->Clone("eff_tracking");
  ratio->Reset();
  ratio->Divide(muoneta_track,muoneta_all);
  ratio->Write();


  if (ratio!=0) delete ratio;
  ratio = (TH1D*) muoneta_track->Clone("eff_mu");
  ratio->Reset();
  ratio->Divide(muoneta_mu,muoneta_track);
  ratio->Write();


  if (ratio!=0) delete ratio;
  ratio = (TH1D*) muonpt_recod->Clone("eff_idiso");
  ratio->Reset();
  ratio->Divide(muonpt_recod,muonpt_mu);
  ratio->Write();

  fout.Close();

  /*
  TString  vary="minDRmuonJet";
  TString varx="genLeptonPt1";
  logx_=true; //log scale x axis
  setLogY(true);
  draw2d(varx, 5,5, 200, vary, 5,0.3,6,"muon pT","min #Delta R(#mu , jet)","stop_eff2d_MC_muon",0,0);
  */

}



void tauembedtest2_luke() { //5 Dec 2012
  /*
idea from luke to embed _all_ gen-level W->mu events.

this gave strange results. Does not seem to have been useful
  */

  inputPath = "/cu3/joshmt/stop/embedTrees/v3/";
  reducedTreeName_ = "maketree/reducedTree"; 

  //add special test samples
  addToSamplesAll("TTbarEmbedGenAllTauHad");

  loadSamples(true,"stop");
  clearSamples();
  addSample("TTbarJetsPowheg:ttbarDecayCode==7",kAzure-3,"MC #tau #rightarrow h");
  //accidentally embedded dilepton events as well. Enforce decay code 4 to get W->mu sample
  //decay code was calculated with the original PAT genParticles (before embedding)
  addSample("TTbarEmbedGenAllTauHad:ttbarDecayCode==4",kMagenta,"MC embedded #tau #rightarrow h");
  currentConfig_=configDescriptions_.getDefault();
  dodata_=false;
  doRatioPlot(true); ratioMin = 0; ratioMax = 2;

  //these trees have no weight!
  nameOfEventWeight_="(1)";
  lumiScale_=1;

  //now compensate for branching fractions
  setSampleScaleFactor("TTbarEmbedGenAllTauHad:ttbarDecayCode==4",11.25*0.6479 / 10.57); //assumes no skim

  int nbins;
  float low,high;
  TString var,xtitle;
  
  setStackMode(false,false);

  TCut jets="njets70>=2 && njets50>=4 && njets30>=5";
  TCut met = "MET>=175";
  TCut lepveto = "electronpt1 <5 && muonpt1<5";

  selection_ =jets&&met&&lepveto;
  //
  setLogY(false);
  var="njets30"; xtitle="njets (pT>30 GeV)";
  nbins = 8; low=4; high=12;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_njets30",0);

  TCut isotkveto = "isotrackpt1<10";
  selection_ =jets&&met&&lepveto&&isotkveto;
  var="njets30"; xtitle="njets (pT>30 GeV)";
  nbins = 8; low=4; high=12;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_njets30",0);

}

void compareseeds() { //5 Dec 2012


  inputPath = "/cu3/joshmt/stop/embedTrees/v3/";
  reducedTreeName_ = "maketree/reducedTree"; 
  loadSamples(true,"stop");
  clearSamples();
  addSample("TTbarJetsPowheg:ttbarDecayCode==7",kAzure-3,"MC #tau #rightarrow h");
  addSample("TTbarEmbedGenTauHad",kMagenta,"MC embedded #tau #rightarrow h");
  currentConfig_=configDescriptions_.getDefault();
  dodata_=false;
  doRatioPlot(true); ratioMin = 0; ratioMax = 2;

  //these trees have no weight!
  nameOfEventWeight_="(1)";
  lumiScale_=1;


  //now compensate for branching fractions
  setSampleScaleFactor("TTbarEmbedGenTauHad",11.25*0.6479 / 10.57); //assumes no skim

  int nbins;
  float low,high;
  TString var,xtitle;
  
  //pT spectrum of (a) gen-level taus (b) seed sample

  setStackMode(false,false);

  TCut jets="njets70>=2 && njets50>=4 && njets30>=5";
  TCut met = "MET>=175";
  TCut lepveto = "electronpt1 <5 && muonpt1<5";

  // -- pt plot --

  selection_ =jets&&met;//  && TCut("minDRmuonJet>0.5");
  setLogY(false);
  //(a) gen-level taus
  var="genLeptonPt1"; xtitle="gen tau pT";
  nbins = 50; low=0; high=200;
  //  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_njets30",0);
  drawSimple(var, nbins, low, high,"filename","genleveltaupt" ,"TTbarJetsPowheg:ttbarDecayCode==7",0 ,xtitle,"Events");
  TH1D* htaupt=  (TH1D*)hinteractive->Clone("genleveltaupt");
  var="seedmuonpt"; xtitle="seed muon pT";
  drawSimple(var, nbins, low, high,"filename","seedmuonpt" ,"TTbarEmbedGenTauHad",0 ,xtitle,"Events");
  TH1D* hseedpt = (TH1D*)hinteractive->Clone("seedmuonpt");

  renewCanvas();
  htaupt->SetLineColor(kRed); 
  htaupt->SetMarkerSize(0);
  hseedpt->SetMarkerSize(0);
  htaupt->Draw("he");
  hseedpt->SetLineColor(kBlue);
  hseedpt->Draw("He same");

  // -- eta plot --

  selection_ =jets&&met;
  setLogY(false);
  //(a) gen-level taus
  var="genLeptonEta1"; xtitle="gen tau eta";
  nbins = 50; low=-5; high=5;
  drawSimple(var, nbins, low, high,"filename","genleveltaueta" ,"TTbarJetsPowheg:ttbarDecayCode==7",0 ,xtitle,"Events");
  TH1D* htaueta=  (TH1D*)hinteractive->Clone("genleveltaueta");
  var="seedmuoneta"; xtitle="seed muon eta";
  drawSimple(var, nbins, low, high,"filename","seedmuoneta" ,"TTbarEmbedGenTauHad",0 ,xtitle,"Events");
  TH1D* hseedeta = (TH1D*)hinteractive->Clone("seedmuoneta");

  renewCanvas();
  htaueta->SetLineColor(kRed); 
  htaueta->SetMarkerSize(0);
  hseedeta->SetMarkerSize(0);
  htaueta->Draw("he");
  hseedeta->SetLineColor(kBlue);
  hseedeta->Draw("He same");


}

void tauembedtest() { //29 Nov 2012 ; updated 19 Dec
  /*
first real embedding closure test on Steven's agenda
mixed results...sort of promising but certaintly some weird stuff going on

updated for v4

9 Jan -- update for v5

16 Jan -- update for v5b (just switch to raw jets)
  */

  inputPath = "/cu3/joshmt/stop/embedTrees/v5b/";
  reducedTreeName_ = "maketree/reducedTree"; 
  //  addToSamplesAll("TTbarEmbedGenTauHad");
  loadSamples(true,"stop");
  clearSamples();
  //enforce that the MC truth is only the part inside of acceptance
  addSample("TTbarJetsPowheg:ttbarDecayCode==7&&abs(genLeptonEta1)<2.4",kAzure-3,"MC #tau #rightarrow h");
  addSample("TTbarEmbedGenTauHad",kMagenta,"MC embedded #tau #rightarrow h");
    //addSample("TTbarEmbedGenTauHadSwap",kMagenta+2,"MC embedded #tau #rightarrow h (swap)");
  currentConfig_=configDescriptions_.getDefault();
  dodata_=false;
  doRatioPlot(true); ratioMin = 0; ratioMax = 2;

  //these trees have no weight!
  nameOfEventWeight_="(1)";
  lumiScale_=1;


  //now compensate for branching fractions
  setSampleScaleFactor("TTbarEmbedGenTauHad",11.25*0.6479 / 10.57);
  setSampleWeightFactor("TTbarEmbedGenTauHad","1.0/(eff_tr*eff_mu*eff_iso)"); //weight embedded events by event-by-event CS eff correction

  int nbins;
  float low,high;
  TString var,xtitle;
  

  setStackMode(false,false);

  TCut jets="njets70raw>=2 && njets50raw>=4 && njets30raw>=5";
  TCut met = "MET>=175";
  TCut lepveto = "electronpt1 <5 && muonpt1<5";
  TCut isotkveto = "isotrackpt1<10";

  //do not expect anything with lep veto to give good results

  selection_ =jets&&met&&lepveto;
  //
  setLogY(false);
  var="njets30"; xtitle="njets (pT>30 GeV)";
  nbins = 8; low=4; high=12;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_njets30",0);

  selection_ =jets&&met&&lepveto;
  //
  setLogY(false);
  var="njets50"; xtitle="njets (pT>50 GeV)";
  nbins = 8; low=3; high=11;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_njets50",0);

  selection_ = TCut("njets50>=4 && njets30>=5")&&met&&lepveto;
  setLogY(false);
  var="njets70"; xtitle="njets (pT>70 GeV)";
  nbins = 8; low=0; high=8;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_njets70",0);

  selection_ = met&&lepveto;
  setLogY(false);
  var="njets30"; xtitle="njets (pT>30 GeV)";
  nbins = 12; low=0; high=12;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_njets30_nojetcuts",0);


//   clearSamples();
//   addSample("TTbarJetsPowheg:ttbarDecayCode==7",kAzure-3,"MC #tau #rightarrow h");
//   addSample("TTbarEmbedGenTauHad",kMagenta,"MC embedded #tau #rightarrow h");
//   addSample("TTbarJetsPowheg:ttbarDecayCode==4",kBlue,"MC W #rightarrow #mu (rescaled)");
//   setSampleScaleFactor("TTbarJetsPowheg:ttbarDecayCode==4",11.25*0.6479 / 10.57); //assumes no skim

  selection_ = met;
  setLogY(false);
  var="njets30"; xtitle="njets (pT>30 GeV)";
  nbins = 12; low=0; high=12;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_njets30_metcutonly",0);


  selection_ = jets && lepveto;
  var = "MET"; xtitle="E^{T}_{miss} (GeV)";
  nbins = 20; low =0; high=400;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_met",0,"GeV");

  //add isotk veto
  setLogY(false);
  selection_ =jets&&met&&lepveto&&isotkveto;
  var="njets30"; xtitle="njets (pT>30 GeV)";
  nbins = 8; low=4; high=12;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_njets30",0);


  //look at e/mu/isotk distributions.
  selection_ = met &&jets;
  var = "muonpt1"; xtitle="muon pT (GeV)";
  nbins = 20; low =0; high=200;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_muonpt",0,"GeV");


  //look at e/mu/isotk distributions.
  selection_ = met &&jets;
  var = "electronpt1"; xtitle="electron pT (GeV)";
  nbins = 20; low =0; high=100;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_electronpt",0,"GeV");


  //look at e/mu/isotk distributions.
  selection_ = met &&jets;
  var = "isotrackpt1"; xtitle="iso track pT (GeV)";
  nbins = 40; low =-100; high=100;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_isotrackpt",0,"GeV");

  //no lepton vetoes but isotrack veto
  selection_ =jets&&met&&isotkveto;
  var = "MET"; xtitle="E^{T}_{miss} (GeV)";
  nbins = 20; low =175; high=500;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_met_jets-met-isotk",0,"GeV");

  setLogY(false);
  var="njets30raw"; xtitle="njets (raw pT>30 GeV)";
  nbins = 8; low=4; high=12;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_njets30_jets-met-isotk",0);

  var="njets50"; xtitle="njets (pT>50 GeV)";
  nbins = 8; low=4; high=12;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_njets50_jets-met-isotk",0);

  var="DeltaPhiMetJet1"; xtitle="#Delta #phi(jet1, MET)";
  nbins = 24; low=0; high=3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_dp1_jets-met-isotk",0);

  var="DeltaPhiMetJet2"; xtitle="#Delta #phi(jet2, MET)";
  nbins = 24; low=0; high=3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_dp2_jets-met-isotk",0);

  var="DeltaPhiMetJet3"; xtitle="#Delta #phi(jet3, MET)";
  nbins = 24; low=0; high=3;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_dp3_jets-met-isotk",0);

  var="jeteta1"; xtitle="lead jet eta";
  nbins = 30; low=-2.4; high=2.4;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_jeteta1_jets-met-isotk",0);

  //
  setLogY(false);
  var="nbjets30raw"; xtitle="nbjets (raw pT>30 GeV)";
  nbins = 4; low=0; high=4;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_nbjets30_jets-met-isotk",0);

  //some mess down here

  setLogY(false);
  var="njets30"; xtitle="njets (pT>30 GeV)";
  nbins = 8; low=4; high=12;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_njets30_jets-met",0);

  var="njets50"; xtitle="njets (pT>50 GeV)";
  nbins = 8; low=4; high=12;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_njets50_jets-met",0);

  //top tagging
  TCut toptag = "bestTopJetMass>80 && bestTopJetMass<270";
  TCut loosetoptag = "bestTopJetIdx>=0"; //redundant if toptag is applied
  TCut otherb = "remainPassCSVS==1";
  TCut mt2 = "MT2>300";
  TCut mtcuts = "mTbJet+0.5*mTbestTopJet > 500";
  //what is this?   
  // if( !taggingMode_ && type3TopTaggerPtr->pickedRemainingCombfatJetIdx == -1 && oriJetsVec.size()>=6 ) return false; 

  selection_ =jets&&met&&isotkveto && toptag&&otherb&&mt2&&mtcuts;
  var = "MET"; xtitle="E^{T}_{miss} (GeV)";
  nbins = 20; low =175; high=500;
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_met_kitchensink",0,"GeV");


}

void embedtest() {

  inputPath = "/cu3/joshmt/stop/embedTrees/v2/";
  reducedTreeName_ = "maketree/reducedTree"; 
  loadSamples(true,"stop");
  clearSamples();
  addSample("TTbarJetsPowheg:ttbarDecayCode==4",kAzure-3,"MC W #rightarrow #mu");
  addSample("TTbarEmbedFlipRecoMu",kMagenta,"MC embedded #mu");
  currentConfig_=configDescriptions_.getDefault();
  dodata_=false;

  //these trees have no weight!
  nameOfEventWeight_="(1)";
  lumiScale_=1;
  //setSampleScaleFactor("TTbarJetsPowheg",234.0/double(getTree("TTbarJetsPowheg")->GetEntries())); //assumes no skim

  int nbins;
  float low,high;
  TString var,xtitle;
  

  setStackMode(false,false);

  TCut jets="njets70>=2 && njets50>=4 && njets30>=5";
  TCut met = "MET>=175";
  TCut lepveto = "electronpt1 <5 && muonpt1<5";

  selection_ =jets&&met&&lepveto;
  //
  setLogY(false);
  var="njets30"; xtitle="njets (pT>30 GeV)";
  nbins = 8; low=4; high=12;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embed_njets30",0);


  nbins=10; low=175; high=300;
  var="MET";xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embedMu_met",0);

  //now try electrons
  clearSamples();
  addSample("TTbarJetsPowheg:ttbarDecayCode==3",kAzure-3,"MC W #rightarrow e");
  addSample("TTbarEmbedFlipRecoEle",kMagenta,"MC embedded e");
  drawPlots(var,nbins,low,high,xtitle,"Events", "embedEle_njets30",0);

  nbins=10; low=175; high=300;
  var="MET";xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "embedEle_met",0);



}

void stopPAT() {

  inputPath = "/cu3/joshmt/stop/reducedTrees/v0_1/";
  reducedTreeName_ = "MakeTauTree/reducedTree"; //for the TFileService-style output instead of the hand-rolled output

  loadSamples(true,"stop");
  clearSamples();
  addSample("TTbarJets");

  currentConfig_=configDescriptions_.getDefault();
  dodata_=false;

  //these trees have no weight!
  nameOfEventWeight_="(1)";

  lumiScale_=5000;
  setSampleScaleFactor("TTbarJets",165.0/double(getTree("TTbarJets")->GetEntries())); // manually add sigma / Ngen weight factor (7 TeV) //assumes no skim

  is2011Data_ = true;

  //define variables

  int nbins;
  float low,high;
  TString var,xtitle;
  
  //don't stack and don't normalize
  setStackMode(false,false);

  //sanity checks -- raw MET spectrum
  selection_ ="(1)"; //no cuts
  var="MET"; xtitle="MET";
  nbins = 30; low=0; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met",0,"GeV");
  //looks ok i think

  //try splitting by topDecayCategory
  clearSamples();
  addSample("TTbarJets:ttbarDecayCode==4||ttbarDecayCode==6",kMagenta,"#mu, #tau #rightarrow #mu");
  addSample("TTbarJets:ttbarDecayCode==5||ttbarDecayCode==3",kRed,"e, #tau #rightarrow e");
  addSample("TTbarJets:ttbarDecayCode==1||ttbarDecayCode==9",kGreen,"misc e,#mu");
  addSample("TTbarJets:ttbarDecayCode==2",kOrange,"hadronic");
  addSample("TTbarJets:ttbarDecayCode==7||ttbarDecayCode==8",kBlue,"#tau #rightarrow h");
  //  addSample("TTbarJets:ttbarDecayCode==0",1,"?");

  selection_ ="(1)"; //no cuts
  var="MET"; xtitle="MET";
  nbins = 30; low=0; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_bycat",0,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_bycat",0,"GeV");

  //n jets
  var="njets30"; xtitle="njets (pT>30 GeV)";
  nbins = 12; low=0; high=12;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_njets",0,"GeV");


  selection_ ="MET>70"; //get rid of irrelevant hadronic part
  var="nVetoElectrons"; xtitle="n electrons";
  nbins = 3; low=0; high=3;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_nelectrons_met70",0);


  selection_ ="MET>70"; //get rid of irrelevant hadronic part
  var="nVetoMuons"; xtitle="n muons";
  nbins = 3; low=0; high=3;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_nmuons_met70",0);


  //veto electrons + muons
  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0";
  var="nVetoTaus"; xtitle="n taus";
  nbins = 8; low=0; high=8;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_ntaus_met70_eVeto_muVeto",0);

  // == ok, now try to make something more interesting
  clearSamples();
  addSample("TTbarJets:nVetoTaus==0",kRed,"0 reco'd #tau");
  addSample("TTbarJets:nVetoTaus>0",kBlue,">0 reco'd #tau");
  //require tau->had
  selection_ ="nVetoElectrons==0 && nVetoMuons==0 && MET>70 && (ttbarDecayCode==7||ttbarDecayCode==8)";

  //lordy...this is how i have to get the lead tau pt
  var="(tauGenPt[0]>tauGenPt[1])*tauGenPt[0]+(tauGenPt[0]<tauGenPt[1])*tauGenPt[1]"; xtitle="lead gen tau pT";
  nbins = 20; low=0; high=200;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_taugenpt_met70_eVeto_muVeto",0,"GeV");

  //now make the same plot but using the 'visible' gen pT
  var="(tauGenVisPt[0]>tauGenVisPt[1])*tauGenVisPt[0]+(tauGenVisPt[0]<tauGenVisPt[1])*tauGenVisPt[1]"; xtitle="lead gen vis tau pT";
  nbins = 20; low=0; high=100;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_tauvisgenpt_met70_eVeto_muVeto",0,"GeV");


  // try to draw tau rejection eff
  //use cat 7 [one W->tau->h + W->qq' ] only
  clearSamples();
  addSample("TTbarJets");
  var="(tauGenVisPt[0]>tauGenVisPt[1])*tauGenVisPt[0]+(tauGenVisPt[0]<tauGenVisPt[1])*tauGenVisPt[1]"; xtitle="lead gen vis tau pT";
  selection_ = "ttbarDecayCode==7 && nVetoTaus>=1"; //events rejected by veto
  setLogY(false);
  nbins = 20; low=0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_tauvisgenpt_cat7_rejected",0,"GeV");

  TH1D* hrejected = (TH1D*)  hinteractive->Clone("hrejected");

  //now look at the total
  selection_ = "ttbarDecayCode==7"; //all events in this category
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_tauvisgenpt_cat7",0,"GeV");
  TH1D* hcat7 = (TH1D*)  hinteractive->Clone("hcat7");

  if (effgraph != 0) delete effgraph;
  effgraph = new TGraphAsymmErrors();
  effgraph->BayesDivide(hrejected,hcat7);
  renewCanvas();
  effgraph->Draw("AP");

  //OK, so that's great.
  //now, can i check signal efficiency (using hadronic ttbar as a signal proxy)?
  clearSamples();
  addSample("TTbarJets");
  var = "nVetoTaus"; xtitle="n taus";
  setLogY(false);
  nbins = 4; low=0; high=4;
  selection_ = "ttbarDecayCode==2"; //all hadronic ttbar decays
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_nvetotaus_hadronic",0,"GeV");


  // == try looking at jets matched to gen taus?
  //do not look at nVetoTaus. instead look at the componenent variables
  clearSamples();
  addSample("TTbarJets");
  selection_="tauGenPt>0"; // selects entries with a gen tau
  setLogY(false);
  //plot taujetPt spectrum?
  var = "tauGenVisPt"; xtitle="tau gen vis pT"; //xtitle="reco jet pT (matched to gen tau)";
  nbins = 20; low=0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_tauvisgenpt_all",0,"GeV");

  //split the sample into categories
  clearSamples();
  selection_="tauGenPt>0 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0"; // selects entries with a gen tau *and* a tau->h decay
  addSample("TTbarJets");
  //  addSample("TTbarJets:taujetPt>15",kRed,"reco pT>15");
  addSample("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4",kRed-7,"pT>15 & |#eta|<2.4");
  addSample("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898",kRed-9,"pT,|#eta|,!b");
  addSample("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898&&taujetMT<90",kPink-9,"pT,|#eta|,!b,m_{T}<90 GeV");
  addSample("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898&&taujetMT<90&&taujetChm<=4",kGreen,"pT,|#eta|,!b,m_{T},CHM");
  TString finalsamplename = "TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898&&taujetMT<90&&taujetChm<=4&&((taujetPt>=40 && taujetLrm<=0.05)||(taujetPt<40&&(taujetLrm<=(-0.000875*taujetPt+0.0875))))";
  addSample(finalsamplename,kGreen+4,"pT,|#eta|,!b,m_{T},CHM,LRM");



  setLogY(false);
  //plot taujetPt spectrum?
  var = "tauGenVisPt"; xtitle="tau gen vis pT"; //xtitle="reco jet pT (matched to gen tau)";
  nbins = 10; low=0; high=200;
  float varbins[11] = {0,7.5,15,22.5,30,40,50,70,100,150,200};
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_tauvisgenpt_all_split",varbins,"GeV");
  setLogY(true);
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_tauvisgenpt_all_split",varbins,"GeV");
  setLogY(false);
  
  //now make some efficiency plots
  TGraphAsymmErrors*  effgraph1 = new TGraphAsymmErrors();
  effgraph1->BayesDivide(getHist("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898"),getHist("TTbarJets"));

  TGraphAsymmErrors*  effgraph2 = new TGraphAsymmErrors();
  effgraph2->BayesDivide(getHist("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898&&taujetMT<90"),getHist("TTbarJets"));

  TGraphAsymmErrors*  effgraph3 = new TGraphAsymmErrors();
  effgraph3->BayesDivide(getHist("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898&&taujetMT<90&&taujetChm<=4"),getHist("TTbarJets"));

  TGraphAsymmErrors*  effgraph4 = new TGraphAsymmErrors();
  effgraph4->BayesDivide(getHist(finalsamplename),getHist("TTbarJets"));

  renewCanvas();
  effgraph1->SetMarkerColor(kRed-9);
  effgraph2->SetMarkerColor(kPink-9);
  effgraph3->SetMarkerColor(kGreen-9);
  effgraph4->SetMarkerColor(kGreen+4);
  
  effgraph1->Draw("AP"); 
  effgraph2->Draw("P"); 
  effgraph3->Draw("P");
  effgraph4->Draw("P");

  // remake everything with no mT cut
  //need to clean up by hand
  delete effgraph1; delete effgraph2; delete effgraph3; delete effgraph4;
  clearSamples();
  selection_="tauGenPt>0 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0"; // selects entries with a gen tau
  addSample("TTbarJets");
  addSample("TTbarJets:taujetPt>15",kRed,"reco pT>15");
  addSample("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4",kRed-7,"pT>15 & |#eta|<2.4");
  addSample("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898",kRed-9,"pT,|#eta|,!b");
  addSample("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898&&taujetChm<=4",kGreen,"pT,|#eta|,!b,CHM");
  TString finalsamplename2 = "TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898&&taujetChm<=4&&((taujetPt>=40 && taujetLrm<=0.05)||(taujetPt<40&&(taujetLrm<=(-0.000875*taujetPt+0.0875))))";
  addSample(finalsamplename2,kGreen+4,"pT,|#eta|,!b,CHM,LRM");

  setLogY(false);
  //plot taujetPt spectrum?
  var = "tauGenVisPt"; xtitle="tau gen vis pT"; //xtitle="reco jet pT (matched to gen tau)";
  nbins = 10; low=0; high=200;
  float varbins2[11] = {0,7.5,15,22.5,30,40,50,70,100,150,200};
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_tauvisgenpt_all_split_nomt",varbins2,"GeV");
 
  //now make some efficiency plots
  TGraphAsymmErrors*  effgraph1b = new TGraphAsymmErrors();
  effgraph1b->BayesDivide(getHist("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898"),getHist("TTbarJets"));

  TGraphAsymmErrors*  effgraph3b = new TGraphAsymmErrors();
  effgraph3b->BayesDivide(getHist("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898&&taujetChm<=4"),getHist("TTbarJets"));

  TGraphAsymmErrors*  effgraph4b = new TGraphAsymmErrors();
  effgraph4b->BayesDivide(getHist(finalsamplename2),getHist("TTbarJets"));

  renewCanvas();
  effgraph1b->SetMarkerColor(kRed-9);
  effgraph3b->SetMarkerColor(kGreen-9);
  effgraph4b->SetMarkerColor(kGreen+4);
  
  effgraph1b->Draw("AP"); 
  effgraph3b->Draw("P");
  effgraph4b->Draw("P");

  //clean up
  delete effgraph1b;
  delete effgraph3b;
  delete effgraph4b;

  // -- try non cum eff plots --
  //split the sample into categories
  clearSamples();
  selection_="tauGenPt>0 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0"; // selects entries with a gen tau *and* a tau->h decay // && topDecayCode==4
  addSample("TTbarJets");
  //  addSample("TTbarJets:taujetPt>15",kRed,"reco pT>15");
  addSample("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898",kRed-9,"!b");
  addSample("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetMT<90",kPink-9,"m_{T}<90 GeV");
  addSample("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetChm<=4",kGreen,"CHM");
  TString finalsamplename4 = "TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&((taujetPt>=40 && taujetLrm<=0.05)||(taujetPt<40&&(taujetLrm<=(-0.000875*taujetPt+0.0875))))";
  addSample(finalsamplename4,kGreen+4,"LRM");

  setLogY(false);
  var = "tauGenVisPt"; xtitle="tau gen vis pT"; //xtitle="reco jet pT (matched to gen tau)";
  nbins = 10; low=0; high=200;
  float varbins4[11] = {0,7.5,15,22.5,30,40,50,70,100,150,200};
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_tauvisgenpt_all_split",varbins4,"GeV");
  
  //now make some efficiency plots
  TGraphAsymmErrors*  effgraph1c = new TGraphAsymmErrors();
  effgraph1c->BayesDivide(getHist("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898"),getHist("TTbarJets"));

  TGraphAsymmErrors*  effgraph2c = new TGraphAsymmErrors();
  effgraph2c->BayesDivide(getHist("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetMT<90"),getHist("TTbarJets"));

  TGraphAsymmErrors*  effgraph3c = new TGraphAsymmErrors();
  effgraph3c->BayesDivide(getHist("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetChm<=4"),getHist("TTbarJets"));

  TGraphAsymmErrors*  effgraph4c = new TGraphAsymmErrors();
  effgraph4c->BayesDivide(getHist(finalsamplename4),getHist("TTbarJets"));

  renewCanvas();
  effgraph1c->SetMarkerColor(kRed-9);
  effgraph2c->SetMarkerColor(kPink-9);
  effgraph3c->SetMarkerColor(kGreen-9);
  effgraph4c->SetMarkerColor(kGreen+4);
  
  effgraph1c->Draw("AP"); 
  effgraph2c->Draw("P"); 
  effgraph3c->Draw("P");
  effgraph4c->Draw("P");



  // -- make  stack plot --
  setStackMode(false,false);
  clearSamples();
  selection_="tauGenPt>0 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0"; // selects entries with a gen tau
  TString finalsamplename3 = "TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898&&taujetMT<90&&taujetChm<=4&&((taujetPt>=40 && taujetLrm<=0.05)||(taujetPt<40&&(taujetLrm<=(-0.000875*taujetPt+0.0875))))";
  addSample(finalsamplename3,kBlue,"vetoed");
  addSample("TTbarJets:taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898",kRed,"pT,|#eta|,!b");
  addSample("TTbarJets:taujetPt>15",kGreen,"pT");
  //  addSample("TTbarJets:tauGenVisPt>15",kAzure,"vis genTau>15");
  addSample("TTbarJets",kYellow,"all");

  var = "tauGenPt"; xtitle="tau gen pT";
  nbins = 50; low=0; high=300;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_taugenpt_split_fauxstack",0,"GeV");

  // == look at *non*tau jets and tau jets
  setStackMode(false,false);
  clearSamples();
  TString taujets = "TTbarJets:jetMatchGenTau==1&&jetPt>15&&jetEta<2.4&&jetEta>-2.4&&jetCSV<0.898&&jetMT<90&&jetChm<=4&&((jetPt>=40 && jetLrm<=0.05)||(jetPt<40&&(jetLrm<=(-0.000875*jetPt+0.0875))))";
  TString alltaujets = "TTbarJets:jetMatchGenTau==1";
  TString nottaujets = "TTbarJets:jetMatchGenTau==0&&jetPt>15&&jetEta<2.4&&jetEta>-2.4&&jetCSV<0.898&&jetMT<90&&jetChm<=4&&((jetPt>=40 && jetLrm<=0.05)||(jetPt<40&&(jetLrm<=(-0.000875*jetPt+0.0875))))";
  TString allnottaujets = "TTbarJets:jetMatchGenTau==0";
  addSample(taujets,kRed,"rejected taus");
  addSample(alltaujets,kBlue,"all tau jets");
  addSample(nottaujets,kGreen,"rejected q/g jets");
  addSample(allnottaujets,kMagenta,"all q/g jets");

  selection_ ="(1)";
  var = "jetPt"; xtitle="reco jet pT";
  nbins = 40; low=0; high=200;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_recojets",0,"GeV");


  //now make some efficiency plots
  TGraphAsymmErrors*  effgraph1d = new TGraphAsymmErrors();
  effgraph1d->BayesDivide(getHist(taujets),getHist(alltaujets));

  TGraphAsymmErrors*  effgraph2d = new TGraphAsymmErrors();
  effgraph2d->BayesDivide(getHist(nottaujets),getHist(allnottaujets));
  renewCanvas();
  effgraph1d->SetMarkerColor(kRed);
  effgraph2d->SetMarkerColor(kBlue);

  effgraph1d->Draw("AP"); 
  effgraph2d->Draw("P"); 

  

}

void stopPAT_PU() {

  inputPath = "/cu3/joshmt/stop/reducedTrees/v0_1/";
  reducedTreeName_ = "MakeTauTree/reducedTree"; //for the TFileService-style output instead of the hand-rolled output

  loadSamples(true,"stop");
  clearSamples();
  addSample("TTbarJets");
  is2011Data_ = true; //use true for 7TeV and false for 8TeV

  currentConfig_=configDescriptions_.getDefault();
  dodata_=false;

  //these trees have no weight!
  nameOfEventWeight_="(1)";

  lumiScale_=5000;
  setSampleScaleFactor("TTbarJets",165.0/double(getTree("TTbarJets")->GetEntries())); // manually add sigma / Ngen weight factor (7 TeV) //assumes no skim


  //define variables

  int nbins;
  float low,high;
  TString var,xtitle;
  
  //don't stack and don't normalize
  setStackMode(false,false);

  clearSamples();
  selection_="tauGenPt>0 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0"; // selects entries with a gen tau

  TString allLow = "TTbarJets:nGoodPV<=6";
  TString vetoLow = "TTbarJets:nGoodPV<=6&&taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898&&taujetMT<90&&taujetChm<=4&&((taujetPt>=40 && taujetLrm<=0.05)||(taujetPt<40&&(taujetLrm<=(-0.000875*taujetPt+0.0875))))";

  TString allMed = "TTbarJets:nGoodPV>6&&nGoodPV<=11";
  TString vetoMed = "TTbarJets:nGoodPV>6&&nGoodPV<=11&&taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898&&taujetMT<90&&taujetChm<=4&&((taujetPt>=40 && taujetLrm<=0.05)||(taujetPt<40&&(taujetLrm<=(-0.000875*taujetPt+0.0875))))";

  TString allHigh = "TTbarJets:nGoodPV>11";
  TString vetoHigh = "TTbarJets:nGoodPV>11&&taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898&&taujetMT<90&&taujetChm<=4&&((taujetPt>=40 && taujetLrm<=0.05)||(taujetPt<40&&(taujetLrm<=(-0.000875*taujetPt+0.0875))))";

  addSample(allLow);
  addSample(vetoLow);
  addSample(allMed);
  addSample(vetoMed);
  addSample(allHigh);
  addSample(vetoHigh);

  setLogY(false);
  var = "tauGenVisPt"; xtitle="tau gen vis pT"; //xtitle="reco jet pT (matched to gen tau)";
  nbins = 10; low=0; high=200;
  float varbins[11] = {0,7.5,15,22.5,30,40,50,70,100,150,200};
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_tauvisgenpt_all_split",varbins,"GeV");

  //now make eff plots
  TGraphAsymmErrors*  effgraph1a = new TGraphAsymmErrors();
  effgraph1a->BayesDivide(getHist(vetoLow),getHist(allLow));

  TGraphAsymmErrors*  effgraph1b = new TGraphAsymmErrors();
  effgraph1b->BayesDivide(getHist(vetoMed),getHist(allMed));

  TGraphAsymmErrors*  effgraph1c = new TGraphAsymmErrors();
  effgraph1c->BayesDivide(getHist(vetoHigh),getHist(allHigh));

  renewCanvas();
  effgraph1a->SetLineColor(kRed);
  effgraph1b->SetLineColor(kGreen);
  effgraph1c->SetLineColor(kBlue);
  effgraph1a->SetMarkerColor(kRed);
  effgraph1b->SetMarkerColor(kGreen);
  effgraph1c->SetMarkerColor(kBlue);

  effgraph1a->Draw("AP");
  effgraph1b->Draw("P");
  effgraph1c->Draw("P");

  delete effgraph1a; delete effgraph1b; delete effgraph1c;


  //try splitting by topDecayCategory
  clearSamples();

  addSample("TTbarJets:nGoodPV<=6",kRed,"<7");
  addSample("TTbarJets:nGoodPV>=7&&nGoodPV<=11",kGreen,"7-11");
  addSample("TTbarJets:nGoodPV>=12",kBlue,">11");
  setStackMode(false,true);
  //veto electrons + muons, require tau->h
  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0 && (ttbarDecayCode==7||ttbarDecayCode==8)";
  var="nVetoTaus"; xtitle="n taus";
  nbins = 4; low=0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_ntaus_fullselection_bynpv",0);


}

void stopPAT_PU_8TeV() {

  inputPath = "/cu3/joshmt/stop/reducedTrees/v0_1/";
  reducedTreeName_ = "MakeTauTree/reducedTree"; //for the TFileService-style output instead of the hand-rolled output

  loadSamples(true,"stop");
  clearSamples();
  addSample("TTbarJets8TeV");
  is2011Data_ = false; //use true for 7TeV and false for 8TeV

  currentConfig_=configDescriptions_.getDefault();
  dodata_=false;

  //these trees have no weight!
  nameOfEventWeight_="(1)";

  lumiScale_=5000;
  setSampleScaleFactor("TTbarJets8TeV",234.0/double(getTree("TTbarJets8TeV")->GetEntries())); // manually add sigma / Ngen weight factor  //assumes no skim


  //define variables

  int nbins;
  float low,high;
  TString var,xtitle;
  
  //don't stack and don't normalize
  setStackMode(false,false);

  clearSamples();
  selection_="tauGenPt>0 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0"; // selects entries with a gen tau

  TString allLow = "TTbarJets8TeV:nGoodPV<=18";
  TString vetoLow = "TTbarJets8TeV:nGoodPV<=18&&taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898&&taujetMT<90&&taujetChm<=4&&((taujetPt>=40 && taujetLrm<=0.05)||(taujetPt<40&&(taujetLrm<=(-0.000875*taujetPt+0.0875))))";

  TString allMed = "TTbarJets8TeV:nGoodPV>18&&nGoodPV<=28";
  TString vetoMed = "TTbarJets8TeV:nGoodPV>18&&nGoodPV<=28&&taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898&&taujetMT<90&&taujetChm<=4&&((taujetPt>=40 && taujetLrm<=0.05)||(taujetPt<40&&(taujetLrm<=(-0.000875*taujetPt+0.0875))))";

  TString allHigh = "TTbarJets8TeV:nGoodPV>28";
  TString vetoHigh = "TTbarJets8TeV:nGoodPV>28&&taujetPt>15&&taujetEta<2.4&&taujetEta>-2.4&&taujetCSV<0.898&&taujetMT<90&&taujetChm<=4&&((taujetPt>=40 && taujetLrm<=0.05)||(taujetPt<40&&(taujetLrm<=(-0.000875*taujetPt+0.0875))))";

  addSample(allLow);
  addSample(vetoLow);
  addSample(allMed);
  addSample(vetoMed);
  addSample(allHigh);
  addSample(vetoHigh);

  setLogY(false);
  var = "tauGenVisPt"; xtitle="tau gen vis pT"; //xtitle="reco jet pT (matched to gen tau)";
  nbins = 10; low=0; high=200;
  float varbins[11] = {0,7.5,15,22.5,30,40,50,70,100,150,200};
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop8TeV_tauvisgenpt_all_split",varbins,"GeV");

  //now make eff plots
  TGraphAsymmErrors*  effgraph1a = new TGraphAsymmErrors();
  effgraph1a->BayesDivide(getHist(vetoLow),getHist(allLow));

  TGraphAsymmErrors*  effgraph1b = new TGraphAsymmErrors();
  effgraph1b->BayesDivide(getHist(vetoMed),getHist(allMed));

  TGraphAsymmErrors*  effgraph1c = new TGraphAsymmErrors();
  effgraph1c->BayesDivide(getHist(vetoHigh),getHist(allHigh));

  renewCanvas();
  effgraph1a->SetLineColor(kRed);
  effgraph1b->SetLineColor(kGreen);
  effgraph1c->SetLineColor(kBlue);
  effgraph1a->SetMarkerColor(kRed);
  effgraph1b->SetMarkerColor(kGreen);
  effgraph1c->SetMarkerColor(kBlue);

  effgraph1a->Draw("AP");
  effgraph1b->Draw("P");
  effgraph1c->Draw("P");

  delete effgraph1a; delete effgraph1b; delete effgraph1c;


  //try splitting by topDecayCategory
  clearSamples();

  addSample("TTbarJets",kRed,"7 TeV");
  addSample("TTbarJets8TeV",kBlue,"8 TeV");
  setStackMode(false,true);
  //veto electrons + muons, require tau->h
  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0 && (ttbarDecayCode==7||ttbarDecayCode==8)";
  var="nVetoTaus"; xtitle="n taus";
  nbins = 4; low=0; high=4;
  setLogY(false);
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_ntaus_fullselection_byEnergy",0);


}

void distributions() {


  inputPath = "/cu3/joshmt/stop/reducedTrees/v0_1/";
  reducedTreeName_ = "MakeTauTree/reducedTree"; //for the TFileService-style output instead of the hand-rolled output

  loadSamples(true,"stop");
  clearSamples();
  addSample("TTbarJets");

  currentConfig_=configDescriptions_.getDefault();
  dodata_=false;

  //these trees have no weight!
  nameOfEventWeight_="(1)";

  lumiScale_=5000;
  setSampleScaleFactor("TTbarJets",165.0/double(getTree("TTbarJets")->GetEntries())); // manually add sigma / Ngen weight factor (7 TeV) //assumes no skim

  is2011Data_ = true;

  //define variables

  int nbins;
  float low,high;
  TString var,xtitle;
  
  //don't stack
  setStackMode(false,true);

  clearSamples();
  selection_="jetPt>20 && jetEta<2.4 && jetEta>-2.4 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0"; // selects entries with a gen tau


  clearSamples();
  addSample("TTbarJets:jetMatchGenTau==1",kRed,"jets matched to gen tau");
  addSample("TTbarJets:jetMatchGenTau==0",kBlue,"jets not matched");


  nbins=30; low=-1.05; high=1.2;
  var = "jetCSV"; xtitle="CSV value";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_csv");


  nbins=50; low=0; high=500;
  var = "jetMT"; xtitle="jet m_{T}";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_mt");


  nbins=10; low=0; high=10;
  var = "jetChm"; xtitle="jet chm";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_chm");


  nbins=30; low=0; high=0.2;
  var = "jetLrm"; xtitle="jet lrm";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_lrm");
  //new tool to look at eff v rejection curve. depends on drawPlots!
  drawEffVRej("TTbarJets:jetMatchGenTau==1","TTbarJets:jetMatchGenTau==0","tau matched","not tau");

  // try looking at these as a function of PU

  clearSamples();
  addSample("TTbarJets:jetMatchGenTau==1&&nGoodPV<10",kRed,"tau-match (PV<10)");
  addSample("TTbarJets:jetMatchGenTau==1&&nGoodPV>=10",kRed-7,"tau-match (PV>=10)");
  addSample("TTbarJets:jetMatchGenTau==0&&nGoodPV<10",kBlue,"non-tau (PV<10)");
  addSample("TTbarJets:jetMatchGenTau==0&&nGoodPV>=10",kBlue-7,"non-tau (PV>=10)");

  nbins=30; low=0; high=0.2;
  var = "jetLrm"; xtitle="jet lrm";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_lrm_pv");


  nbins=10; low=0; high=10;
  var = "jetChm"; xtitle="jet chm";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_chm_pv");



  nbins=30; low=-1.05; high=1.2;
  var = "jetCSV"; xtitle="CSV value";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_csv_pv");


  nbins=50; low=0; high=500;
  var = "jetMT"; xtitle="jet m_{T}";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_mt_pv");

  //try looking only at low pT jets
  selection_="jetPt>20 && jetPt<60 && jetEta<2.4 && jetEta>-2.4 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0";

  nbins=30; low=0; high=0.2;
  var = "jetLrm"; xtitle="jet lrm";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_lrm_pv_lowpt");



}



void stopPAT_vetoComp() {

  inputPath = "/cu3/joshmt/stop/reducedTrees/v1_1/";
  reducedTreeName_ = "MakeTauTree/reducedTree"; //for the TFileService-style output instead of the hand-rolled output

  loadSamples(true,"stop");
  clearSamples();
  addSample("TTbarJets8TeV");
  is2011Data_ = false; //use true for 7TeV and false for 8TeV

  currentConfig_=configDescriptions_.getDefault();
  dodata_=false;

  //these trees have no weight!
  nameOfEventWeight_="(1)";

  lumiScale_=5000;
  setSampleScaleFactor("TTbarJets8TeV",234.0/double(getTree("TTbarJets8TeV")->GetEntries())); // manually add sigma / Ngen weight factor (7 TeV) //assumes no skim


  //define variables

  int nbins;
  float low,high;
  TString var,xtitle;
  
  //don't stack and don't normalize
  setStackMode(false,false);

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0 && ttbarDecayCode==7"; //select tau->h

    selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0 && (ttbarDecayCode==5 || ttbarDecayCode==6 || ttbarDecayCode==8|| ttbarDecayCode==9)"; //other tau decays

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0 && (ttbarDecayCode==1 || ttbarDecayCode==3 || ttbarDecayCode==4)"; //not tau; ttbar -> lep

  selection_="nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0 && ttbarDecayCode==2"; //ttbar->hadronic

  nbins = 4; low=0; high=4;
  var="nVetoTaus"; xtitle="n Indirect Tau Veto";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_tauh_indirecttau");
  TH1D* indirecttau=(TH1D*)  hinteractive->Clone("indirecttau");

  var="nIsolatedTracks"; xtitle="n Isolated Tracks";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_tauh_isolatedtrack");
  TH1D* isolatedtrack=(TH1D*)  hinteractive->Clone("isolatedtrack");

  var="nVLooseTaus15"; xtitle="n VL POG Taus 15";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_tauh_vloose15");
  TH1D* vloose15 = (TH1D*) hinteractive->Clone("vloose15");

  renewCanvas();
  renewLegend();
  indirecttau->SetLineColor(kRed);
  isolatedtrack->SetLineColor(kBlue);
  vloose15->SetLineColor(kMagenta);

  indirecttau->Draw("hist e");
  isolatedtrack->Draw("hist e same");
  vloose15->Draw("hist e same");

  indirecttau->SetXTitle("n ID'd taus");
  //  indirecttau->SetMaximum(500); //kludge for tau->h
  //indirecttau->SetMaximum(20);
  indirecttau->SetMaximum(160e3); 
  leg->AddEntry(indirecttau,"Indirect");
  leg->AddEntry(isolatedtrack,"IsoTrack");
  leg->AddEntry(vloose15,"VLoose pT>15");
  leg->Draw();

}

void stopPAT_vetoComp2() {

  inputPath = "/cu3/joshmt/stop/reducedTrees/v1_1/";
  reducedTreeName_ = "MakeTauTree/reducedTree"; //for the TFileService-style output instead of the hand-rolled output

  loadSamples(true,"stop");
  clearSamples();
  addSample("TTbarJets8TeV");
  is2011Data_ = false; //use true for 7TeV and false for 8TeV

  currentConfig_=configDescriptions_.getDefault();
  dodata_=false;

  //these trees have no weight!
  nameOfEventWeight_="(1)";

  lumiScale_=5000;
  setSampleScaleFactor("TTbarJets8TeV",234.0/double(getTree("TTbarJets8TeV")->GetEntries())); // manually add sigma / Ngen weight factor (7 TeV) //assumes no skim


  //define variables

  int nbins;
  float low,high;
  TString var,xtitle;
  
  setStackMode(true);  //stacked
  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0";

  nbins = 20; low=175; high=500;
  var="MET"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_preselection");

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0 &&nVetoTaus==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_indirectveto");

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0 &&nIsolatedTracks==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_isotracks");

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0 &&nVLooseTaus15==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_isotracks");

  //try splitting by topDecayCategory
  clearSamples();
  addSample("TTbarJets8TeV:ttbarDecayCode==1",kMagenta-5,"e/mu+e/mu");
  addSample("TTbarJets8TeV:ttbarDecayCode==2",kOrange,"all hadronic");
  addSample("TTbarJets8TeV:ttbarDecayCode==3",kAzure,"e");
  addSample("TTbarJets8TeV:ttbarDecayCode==4",kMagenta,"#mu");
  addSample("TTbarJets8TeV:ttbarDecayCode==5",kRed-5,"#tau #rightarrow e");
  addSample("TTbarJets8TeV:ttbarDecayCode==6",kRed,"#tau #rightarrow #mu");
  addSample("TTbarJets8TeV:ttbarDecayCode==7",kBlue,"#tau #rightarrow h");
  addSample("TTbarJets8TeV:ttbarDecayCode==8",kGreen,"#tau + #tau (all)");
  addSample("TTbarJets8TeV:ttbarDecayCode==9",kYellow,"#tau + e/#mu");
  addSample("T2tt$425$50:ttbarDecayCode==2",kOrange,"T2tt (had)");

  setSampleScaleFactor("T2tt",0.243755 / 50000.0); //manually set xsec / N factor                                                         
  smsHack_=true; //essential for using the T2tt files tau trees
  overrideSMSlabels_=true; //want to use the labels i set above instead of the auto-SMS labels

  setStackMode(true);
  nbins = 20; low=175; high=500;
  var="MET"; xtitle=var;

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_prepreselection");

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_prepreselection_eveto");

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoMuons==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_prepreselection_muveto");

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_preselection");

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nIsolatedTracks==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_prepreselection_isotracks");

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0 &&nVetoTaus==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_indirectveto");

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0 &&nIsolatedTracks==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_isotracks");

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 &&nVetoTaus==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_onlyindirectveto");

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 &&nVLooseTaus15==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_onlyPOGtau");

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 &&nVetoElectrons==0 && nVetoMuons==0 &&nVLooseTaus15==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_emuPOGtau");

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nIsolatedTracks==0 && nVLooseTaus15==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_isoTkPOGtau");

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0 &&nIsolatedTracks==0&&nVetoTaus==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_emutauisotrack");

  //request from rick: pre + ind tau + isotrack
  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoTaus==0 &&nIsolatedTracks==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_met_tauisotrack");



  //more properties
  clearSamples();
  addSample("TTbarJets8TeV:nVetoTaus==1&&nIsolatedTracks==0",kRed,"Indirect only");
  addSample("TTbarJets8TeV:nVetoTaus==0&&nIsolatedTracks==1",kBlue,"isoTrack only");
  addSample("TTbarJets8TeV:nVetoTaus==1&&nIsolatedTracks==1",kBlack,"both");
  addSample("TTbarJets8TeV:nVetoTaus==0&&nIsolatedTracks==0",kGreen,"neither");

  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && ttbarDecayCode==7";

  setStackMode(false,false);
  var="(tauGenVisPt[0]>tauGenVisPt[1])*tauGenVisPt[0]+(tauGenVisPt[0]<tauGenVisPt[1])*tauGenVisPt[1]"; xtitle="lead gen vis tau pT";
  nbins = 20; low=0; high=100;
  drawPlots(var,nbins,low,high,xtitle,"Events", "stop_taugenvispt_byveto");

  //take a close look at the "isotrack only" guys
  selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && ttbarDecayCode==7 && jetMatchGenTau==1";

  clearSamples();
  addSample("TTbarJets8TeV:nVetoTaus==1&&nIsolatedTracks==0",kRed,"Indirect only");
  addSample("TTbarJets8TeV:nVetoTaus==0&&nIsolatedTracks==1",kBlue,"isoTrack only");
  addSample("TTbarJets8TeV:nVetoTaus==1&&nIsolatedTracks==1",kBlack,"both");
  addSample("TTbarJets8TeV:nVetoTaus==0&&nIsolatedTracks==0",kGreen,"neither");


  nbins=30; low=0; high=0.2;
  var = "jetLrm"; xtitle="jet lrm";
  drawPlots(var,nbins,low,high,xtitle,"", "stop_jets_lrm_byid");


  nbins=10; low=0; high=10;
  var = "jetChm"; xtitle="jet chm";
  drawPlots(var,nbins,low,high,xtitle,"", "stop_jets_chm_byid");



  nbins=30; low=-1.05; high=1.2;
  var = "jetCSV"; xtitle="CSV value";
  drawPlots(var,nbins,low,high,xtitle,"", "stop_jets_csv_byid");


  nbins=50; low=0; high=500;
  var = "jetMT"; xtitle="jet m_{T}";
  drawPlots(var,nbins,low,high,xtitle,"", "stop_jets_mt_byid");


  //selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && ttbarDecayCode==7";
  //  getTree("TTbarJets8TeV")->Draw("nVetoTaus:nIsolatedTracks",getCutString(false,"","4.00570590552320108e-05"),"colz")


}

void effVrej_7TeV() {

  //first 7 TeV
  inputPath = "/cu3/joshmt/stop/reducedTrees/v1_1/";
  reducedTreeName_ = "MakeTauTree/reducedTree"; //for the TFileService-style output instead of the hand-rolled output

  loadSamples(true,"stop");
  clearSamples();
  addSample("TTbarJets8TeV");

  currentConfig_=configDescriptions_.getDefault();
  dodata_=false;

  //these trees have no weight!
  nameOfEventWeight_="(1)";

  lumiScale_=5000;
  setSampleScaleFactor("TTbarJets8TeV",234.0/double(getTree("TTbarJets8TeV")->GetEntries())); // manually add sigma / Ngen weight factor (7 TeV) //assumes no skim

  is2011Data_ = false;

  //define variables

  int nbins;
  float low,high;
  TString var,xtitle;
  
  //don't stack
  setStackMode(false,true);

  clearSamples();
  addSample("TTbarJets8TeV:jetMatchGenTau==1",kRed,"jets matched to gen tau");
  addSample("TTbarJets8TeV:jetMatchGenTau==0",kBlue,"jets not matched");


  nbins=30; low=0; high=0.2;
  var = "jetLrm"; xtitle="jet lrm";

  //split into bins of 20 GeV in pT
  selection_="jetPt>10&&jetPt<30 && jetEta<2.4 && jetEta>-2.4 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0"; // jet pt 10-30
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_lrm_pt10-30");

  selection_="jetPt>30&&jetPt<50 && jetEta<2.4 && jetEta>-2.4 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_lrm_pt30-50");


  selection_="jetPt>50&&jetPt<70 && jetEta<2.4 && jetEta>-2.4 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_lrm_pt50-70");

  selection_="jetPt>70&&jetPt<90 && jetEta<2.4 && jetEta>-2.4 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_lrm_pt70-90");

  selection_="jetPt>90 && jetEta<2.4 && jetEta>-2.4 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_lrm_pt90");

  // -- now do chm
  //split into bins of 20 GeV in pT

  nbins=10; low=0; high=10;
  var = "jetChm"; xtitle="jet chm";

  selection_="jetPt>10&&jetPt<30 && jetEta<2.4 && jetEta>-2.4 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0"; // jet pt 10-30
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_chm_pt10-30");

  selection_="jetPt>30&&jetPt<50 && jetEta<2.4 && jetEta>-2.4 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_chm_pt30-50");


  selection_="jetPt>50&&jetPt<70 && jetEta<2.4 && jetEta>-2.4 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_chm_pt50-70");

  selection_="jetPt>70&&jetPt<90 && jetEta<2.4 && jetEta>-2.4 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_chm_pt70-90");

  selection_="jetPt>90 && jetEta<2.4 && jetEta>-2.4 && MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0";
  drawPlots(var,nbins,low,high,xtitle,"au", "stop_jets_chm_pt90");

  //new tool to look at eff v rejection curve. depends on drawPlots!
  //  drawEffVRej("TTbarJets:jetMatchGenTau==1","TTbarJets:jetMatchGenTau==0","tau matched","not tau");


}

void T2ttsanity() {

  smsHack_=true; //avoid use of scansmsngen and all of the other ra2b infrastructure

    inputPath = "/cu3/joshmt/stop/reducedTrees/v1_1/";
    reducedTreeName_ = "MakeTauTree/reducedTree";

    loadSamples(true,"stop");
    clearSamples();
    addSample("T2tt$425$50",kMagenta);

    currentConfig_=configDescriptions_.getDefault();
    dodata_=false;

    //these trees have no weight!                                                                                                                     
    nameOfEventWeight_="(1)";

    lumiScale_=5000;

    // the sampleScaleFactor may not work for T2tt in quite the expected way
    setSampleScaleFactor("T2tt",0.243755 / 50000.0); //manually set xsec / N factor
    //also just getting the number of entries in the tree doesn't work either
    //for now just use an arbitrary normalization ....

    is2011Data_ = false;
    
    int nbins;
    float low,high;
    TString var,xtitle;


    //don't stack and don't normalize
    setStackMode(false,false);

    //sanity checks -- raw MET spectrum                                                                                    
    selection_ ="(1)"; //no cuts         

    var="MET"; xtitle="MET";
    nbins = 40; low=0; high=400;
    drawPlots(var,nbins,low,high,xtitle,"Events", "stop_t2tt_met",0,"GeV");

    //seems sane i guess

    selection_="MET>175 && nbjets30>=1 && njets30>=5 && njets50>=3 && nVetoElectrons==0 && nVetoMuons==0";  //preselection

    //try splitting by topDecayCategory                                                                           
    //in tests, this almost still works in T2tt...maybe not quite as well as for ttbar

    clearSamples();
    overrideSMSlabels_=true; //want to use the labels i set below instead of the auto-SMS labels
    addSample("T2tt$425$50:ttbarDecayCode==4||ttbarDecayCode==6",kMagenta,"#mu, #tau #rightarrow #mu");
    addSample("T2tt$425$50:ttbarDecayCode==5||ttbarDecayCode==3",kRed,"e, #tau #rightarrow e");
    addSample("T2tt$425$50:ttbarDecayCode==1||ttbarDecayCode==9",kGreen,"misc e,#mu");
    addSample("T2tt$425$50:ttbarDecayCode==2",kOrange,"hadronic");
    addSample("T2tt$425$50:ttbarDecayCode==7||ttbarDecayCode==8",kBlue,"#tau #rightarrow h");
    addSample("T2tt$425$50:ttbarDecayCode==0",1,"?");               

    var="MET"; xtitle="MET";
    nbins = 40; low=0; high=400;
    drawPlots(var,nbins,low,high,xtitle,"Events", "stop_t2tt_met_bycat",0,"GeV");

}
