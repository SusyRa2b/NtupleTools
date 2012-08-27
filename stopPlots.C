/*
to compile:

must have symlink to MiscUtil.cxx in the working directory.

gSystem->Load("TSelectorMultiDraw_C.so");
gSystem->Load("CrossSectionTable_cxx.so");
gSystem->Load("ConfigurationDescriptions_cxx.so");
gSystem->Load("SearchRegion_cxx.so");
gSystem->Load("SystInfo_cxx.so");
gSystem->Load("SignalEffData_cxx.so");
.L stopPlots.C++

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
#include "stopPlots.h"


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
  setSampleScaleFactor("TTbarJets8TeV",165.0/double(getTree("TTbarJets8TeV")->GetEntries())); // manually add sigma / Ngen weight factor (7 TeV) //assumes no skim


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
