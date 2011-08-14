/*
====== this is the nominal code for drawing RA2b data/MC comparison plots =======
2011 updates -- runDataQCD2011() computes all numbers needed for QCD estimate including systematics


-- drawPlots() -- main plotting routine to draw pretty stacks of MC with data on top
-- drawSimple() -- very simple function to draw exactly one sample and put it in a file
-- drawR() -- draws r(MET). a bit more kludgely, but works.
-- drawSignificance -- draws S/sqrt(B) as a function of the variable

The functions that do the heavy lifting, along with various utility functions,
are in drawReducedTrees.h.

Samples to plot are currently hard-coded in loadSamples().
New functions addSample() and removeSample() add some run-time flexibility.

Various routines here call the above functions. Some of them (e.g. drawSomething() ) I usually
use in a quasi-interactive mode.

-- drawOwen() -- uses drawPlots() and drawSimple() to make a file with the histograms needed for
toys and fits to data.

Details:
This code replaces drawCutflowPlots.C (which had replaced drawBasicPlots.C). 
The input files are the reducedTrees.

Depending on the exact setup in loadSamples(), the QCD and Single-top samples must be added together with 'hadd' 
in order to get one file per sample.

Potential improvements:
 -- there are becoming way too many configuration options. i've stopped adding setter functions for them out of laziness
 -- the calls to renormBins() currently have a kludge that hard-codes the reference bin to 2

 -- The cuts are defined by the user settings the selection_ string
directly. This can be error prone. A better interface would allow the user
to set cuts more intutitively and with less change of e.g. accidentally 
forgotten cuts (while still preserving the current flexibility).
 -- someday I should test whether I can get rid of duplicate
functionality for TH1F and TH1D e.g. the case of addOverflowBin()
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

// can be checked out of UserCode/joshmt
#include "MiscUtil.cxx"

#include <fstream>

#include <iostream>
#include <map>
#include <set>
	

//for rerunning the final background estimates...need standard MC and MET cleaning
//TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-24_fullpf2pat/";//path for MC	     
//TString dataInputPath = "/cu3/wteo/reducedTrees/V00-02-05_v3-pickevents/"; //includes MET cleaning but uses a tight skim (not good for plots)

//for making the standard set of plots...need standard MC and all data
TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-24_fullpf2pat/";//path for MC	     
TString dataInputPath = "/cu2/ra2b/reducedTrees/V00-02-24_fullpf2pat/"; //sym links to V00-02-05_v3

//for special tests and signal systematics
//TString inputPath = "/tmp/joshmt/";
//TString inputPath = "/home/joshmt/";//path for MC



//the cutdesc string is now defined in loadSamples()

double lumiScale_ = 1091.891;

#include "drawReducedTrees.h"

SignalEffData signalSystematics2011(const SearchRegion & region, bool isSL=false, bool isLDP=false, TString sampleOfInterest="LM9", int m0=0,int m12=0) {
  useFlavorHistoryWeights_=false;
  loadSamples(); //needs to come first! this is why this should be turned into a class, with this guy in the ctor

  dodata_=false;
  savePlots_ =false;
  setQuiet(true);

  SignalEffData results;

  region.Print();
  assert ( !( isSL&&isLDP));
  if (isLDP) cout<<"   == LDP "<<endl;
  if (isSL) cout<<"   == SL "<<endl;

  //for susy scan; harmless for non-susy scan
  m0_=m0;
  m12_=m12;

  clearSamples();
  addSample(sampleOfInterest);

  TString btagselection = region.btagSelection;
  TCut HTcut=region.htSelection.Data(); 
  TCut SRMET = region.metSelection.Data();
  TCut bcut = "nbjetsSSVHPT>=1";
  TString bweightstring="probge1";
  if (btagselection=="ge2b") {
    bcut="nbjetsSSVHPT>=2";
    bweightstring="probge2";
  }
  else if (btagselection=="ge1b") {}
  else if (btagselection=="ge3b") {
    bcut="nbjetsSSVHPT>=3";
    bweightstring="probge3"; //this one isn't calculated right now, i think
  }
  else {assert(0);}

  TCut baseline = "cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1";
  if (isSL) baseline = "cutPV==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0))";

  TCut passOther = "minDeltaPhiN>=4";
  if (isLDP) passOther="minDeltaPhiN<4";

  selection_ = baseline && HTcut && passOther && SRMET && bcut;

  //critical to reset these!
  usePUweight_=false;
  useHLTeff_=false;
  btagSFweight_="1";
  currentConfig_=configDescriptions_[0]; //completely raw MC

  //the logic here is very fragile to user error, so try to ensure that this is really completely raw MC
  cout<<currentConfig_<<endl;
  assert( currentConfig_.Contains("JER0") && currentConfig_.Contains("JES0") && currentConfig_.Contains("METunc0")
	  && currentConfig_.Contains("PUunc0")&& currentConfig_.Contains("BTagEff0")&& currentConfig_.Contains("HLTEff0") );

  TString var="HT"; TString xtitle=var;
  int  nbins=1; float low=0; float high=1e9;

  /*
some explanation:
the completely raw yield is needed to give to owen, both for comparison purposes and also
to calculate the net effect of all corrections = corrected yield / raw yield

But for calculating systematics, I want to compare the fractional change to the JERbias, because
each of the systematics variations was done with JERbias turned on.

So the convention is:
configDescriptions_[0] is raw
configDescriptions_[1] is with JERbias
  */

  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  double raw0=  getIntegral(sampleOfInterest);
  double rawerr0=  getIntegralErr(sampleOfInterest);
  cout<<"raw yield = "<<raw0<<" +/- "<<rawerr0<<endl;
  TH1D* totalPdfWeightsCTEQ = (TH1D*) files_[currentConfig_][sampleOfInterest]->Get("pdfWeightSumCTEQ");
  const  float nominalweight=   totalPdfWeightsCTEQ->GetBinContent(1);

  TH1D* totalPdfWeightsMSTW = (TH1D*) files_[currentConfig_][sampleOfInterest]->Get("pdfWeightSumMSTW");
  cout<< nominalweight<<"\t"<<   totalPdfWeightsMSTW->GetBinContent(1)<<endl;

  TH1D* totalPdfWeightsNNPDF = (TH1D*) files_[currentConfig_][sampleOfInterest]->Get("pdfWeightSumNNPDF");
  cout<< nominalweight<<"\t"<<   totalPdfWeightsNNPDF->GetBinContent(1)<<endl;

  currentConfig_=configDescriptions_[1]; //add  JER bias
  //the logic here is very fragile to user error, so try to ensure that this is the correct thing
  assert( currentConfig_.Contains("JERbias") && currentConfig_.Contains("JES0") && currentConfig_.Contains("METunc0")
	  && currentConfig_.Contains("PUunc0")&& currentConfig_.Contains("BTagEff0")&& currentConfig_.Contains("HLTEff0") );

  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  cout<<"add JER = "<<getIntegral(sampleOfInterest)<<" +/- "<<getIntegralErr(sampleOfInterest) <<endl;

  //this will be used for PDF uncertainties
  TTree* thetree = (TTree*) files_[currentConfig_][sampleOfInterest]->Get("reducedTree");

  usePUweight_=true;
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  cout<<"add PU weight yield = "<<getIntegral(sampleOfInterest)<<" +/- "<<getIntegralErr(sampleOfInterest) <<endl;

  useHLTeff_=true;
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  cout<<"add HLT eff = "<<getIntegral(sampleOfInterest)<<" +/- "<<getIntegralErr(sampleOfInterest) <<endl;

  selection_ = baseline && HTcut && passOther && SRMET; //remove b tag
  btagSFweight_=bweightstring;
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  double nominal=  getIntegral(sampleOfInterest);
  double nominalerr=  getIntegralErr(sampleOfInterest);
  cout<<"add b tag SF = "<<nominal<<" +/- "<<nominalerr <<endl;

  results.effCorr = nominal / raw0;
  results.rawYield = raw0;

  // NLO k factor uncertainty
  susyCrossSectionVariation_="Plus";
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  double kfactorplus= getIntegral(sampleOfInterest);

  susyCrossSectionVariation_="Minus";
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  double kfactorminus= getIntegral(sampleOfInterest);

  //do everything in %
  kfactorplus=  100*fabs(kfactorplus-nominal)/nominal;
  kfactorminus= 100*fabs(kfactorminus-nominal)/nominal;
  cout<<"k factor uncertainties = "<<kfactorminus<<" "<<kfactorplus<<endl;
  results.totalSystematic = pow(0.5*(kfactorplus+kfactorminus),2); //square of the error

  susyCrossSectionVariation_=""; //reset to default k factors
  //all other uncertainties

  //we're gonna hard code the fact that these variations come in pairs
  double var1=0,var2=0;
  for (unsigned int j=2; j<configDescriptions_.size(); j++) {
    currentConfig_=configDescriptions_[j];
    drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
    double thisn=getIntegral(sampleOfInterest);
    if (j%2==0) { //this one runs first
      var2=fabs(100*(thisn-nominal)/nominal);
    }
    else { //this one runs second
      var1=fabs(100*(thisn-nominal)/nominal);
      double bigger = var2>var1? var2:var1;
      cout<<getVariedSubstring(currentConfig_)<<"\t"<<bigger<<endl;
      results.totalSystematic += bigger*bigger; //add in quadrature
    }
  }

  //reset to the nominal with JER bias
  currentConfig_=configDescriptions_[1];

  double largest=0;
  double smallest=1e9;

  //debug
//   TH1D*  HyieldsCTEQ = new TH1D("HyieldsCTEQ","CTEQ yields",100,0,200);
//   TH1D*  HyieldsMSTW = new TH1D("HyieldsMSTW","MSTW yields",100,0,200);
//   TH1D*  HyieldsNNPDF= new TH1D("HyieldsNNPDF","NNPDF yields",100,0,200);

  double percentCTEQ=0,percentMSTW=0,percentNNPDF=0;
  if (false) { //don't do pdf uncertainties
  // ~~~~~~~~~~~~ CTEQ ~~~~~~~~~~
  { //because i'm using horrible unstructured code here, give the variables a scope
  double Xplus[22];// Xplus[0]=0;
  double Xminus[22];// Xminus[0]=0;
  for (int i=1; i<=44; i++) { //there are 44+1 pdf weights
    //get the sum of weight[i] for the whole sample with no cuts
    float totalweight=   totalPdfWeightsCTEQ->GetBinContent(i+1);
    TString extraWeight=""; extraWeight+=nominalweight; extraWeight+="/"; extraWeight+=totalweight;
    //then get the yield after cuts with extra weight pdfWeights[i]/sum
    TH1D pdfEventCounter("pdfEventCounter","pdfEventCounter",1,0,1e9);
    pdfEventCounter.Sumw2();
    thetree->Project("pdfEventCounter","HT",getCutString(false, "",extraWeight,i,"CTEQ").Data());
    if (i%2==0) {
      Xminus[(i-2)/2] = pdfEventCounter.Integral();
    }
    else {
      Xplus[(i-1)/2] = pdfEventCounter.Integral();
    }
    //for legacy comparison
    if (pdfEventCounter.Integral() > largest) largest = pdfEventCounter.Integral();
    if (pdfEventCounter.Integral() < smallest) smallest = pdfEventCounter.Integral();
    //    HyieldsCTEQ->Fill(pdfEventCounter.Integral());
    pdfEventCounter.Reset();
  }

  double XplusMax = 0;
  double XminusMax = 0;
  for (int i=0; i<22; i++) {
    double diff1 =  Xplus[i] - nominal ;
    double diff2 =  Xminus[i] - nominal ;
    double larger = diff1>diff2 ? diff1 : diff2;
    if ( 0 > larger ) larger = 0;
    XplusMax += larger*larger;
  }
  for (int i=0; i<22; i++) {
    double diff1 =  nominal - Xplus[i] ;
    double diff2 =  nominal - Xminus[i] ;
    double larger = diff1>diff2 ? diff1 : diff2;
    if ( 0 > larger ) larger = 0;
    XminusMax += larger*larger;
  }
  XplusMax = sqrt(XplusMax);
  XminusMax = sqrt(XminusMax);

  //  double pdfunc= ( fabs((smallest-nominal)/nominal) > fabs((largest-nominal)/nominal) ) ? fabs((smallest-nominal)/nominal):fabs((largest-nominal)/nominal);
  //  cout<<"PDF old method\t"<<100*pdfunc<<endl;
  const double scale = 1.645;

  cout<<"CTEQ old method\t"<<nominal-smallest<<"\t"<<largest-nominal<<endl;

  cout<<"CTEQ new method\t"<<XminusMax/scale<<"\t"<<XplusMax/scale<<endl;
  percentCTEQ = 0.5* ( XminusMax/scale + XplusMax/scale);
  percentCTEQ = 100* percentCTEQ/nominal;
  }

  // ~~~~~~~~~~~~ MSTW ~~~~~~~~~~
  // i really dislike copy/pasting code like this, but i'm tired....
  {
  double Xplus[20];// Xplus[0]=0;
  double Xminus[20];// Xminus[0]=0;
  for (int i=1; i<=40; i++) { //there are 40+1 pdf weights
    //get the sum of weight[i] for the whole sample with no cuts
    float totalweight=   totalPdfWeightsMSTW->GetBinContent(i+1);
    TString extraWeight=""; extraWeight+=nominalweight; extraWeight+="/"; extraWeight+=totalweight;
    //then get the yield after cuts with extra weight pdfWeights[i]/sum
    TH1D pdfEventCounter("pdfEventCounter","pdfEventCounter",1,0,1e9);
    pdfEventCounter.Sumw2();
    thetree->Project("pdfEventCounter","HT",getCutString(false, "",extraWeight,i,"MSTW").Data());
    if (i%2==0) {
      Xminus[(i-2)/2] = pdfEventCounter.Integral();
    }
    else {
      Xplus[(i-1)/2] = pdfEventCounter.Integral();
    }
    //    HyieldsMSTW->Fill(pdfEventCounter.Integral());
  }

  double XplusMax = 0;
  double XminusMax = 0;
  for (int i=0; i<20; i++) {
    double diff1 =  Xplus[i] - nominal ;
    double diff2 =  Xminus[i] - nominal ;
    double larger = diff1>diff2 ? diff1 : diff2;
    if ( 0 > larger ) larger = 0;
    XplusMax += larger*larger;
  }
  for (int i=0; i<20; i++) {
    double diff1 =  nominal - Xplus[i] ;
    double diff2 =  nominal - Xminus[i] ;
    double larger = diff1>diff2 ? diff1 : diff2;
    if ( 0 > larger ) larger = 0;
    XminusMax += larger*larger;
  }
  XplusMax = sqrt(XplusMax);
  XminusMax = sqrt(XminusMax);
  cout<<"MSTW new method\t"<<XminusMax<<"\t"<<XplusMax<<endl;
  percentMSTW = 0.5* ( XminusMax + XplusMax);
  percentMSTW = 100* percentMSTW/nominal;

  }

  // ~~~~~~~~~~~~~~~~ NNPDF ~~~~~~~~~~~(see emails from H&S -- just take the RMS as the error)
  {
    double totalX=0,totalX2=0;
    int n=0;
    for (int i=1; i<=99; i++) { //0 is just weight 1. there are 99 others
      //get the sum of weight[i] for the whole sample with no cuts
      float totalweight=   totalPdfWeightsNNPDF->GetBinContent(i+1); //i really did store the total in bin i+1
      TString extraWeight=""; extraWeight+=nominalweight; extraWeight+="/"; extraWeight+=totalweight;
      //then get the yield after cuts with extra weight pdfWeights[i]/sum
      TH1D pdfEventCounter("pdfEventCounter","pdfEventCounter",1,0,1e9);
      pdfEventCounter.Sumw2();
      thetree->Project("pdfEventCounter","HT",getCutString(false, "",extraWeight,i,"NNPDF").Data());
      totalX += pdfEventCounter.Integral();
      totalX2 += pdfEventCounter.Integral()*pdfEventCounter.Integral();
      n++;
      //      HyieldsNNPDF->Fill(pdfEventCounter.Integral());
      pdfEventCounter.Reset();
    }
    totalX/= double(n);
    totalX2/=double(n);
    cout<<"NNPDF method\t"<<sqrt(totalX2-(totalX*totalX))<<endl;
    percentNNPDF = 100*sqrt(totalX2-(totalX*totalX)) / nominal;
  }

  if (percentCTEQ > percentMSTW) percentMSTW=percentCTEQ;
  if (percentNNPDF > percentMSTW) percentMSTW = percentNNPDF;
  cout<<"PDF uncertainty (%) = "<<percentMSTW<<endl;

  results.totalSystematic += percentMSTW*percentMSTW;
  }
  else { 
    //    results.totalSystematic += 13*13;
    cout<<"Skipping PDF uncertainties!"<<endl;
  }

  results.totalSystematic += pow(4.5,2); //lumi uncertainty

  //ok, we're done. take the square root to get the total systematics
  results.totalSystematic = sqrt(results.totalSystematic);

  return results;
}

void runSystematics2011_LM9() {
  setSearchRegions();
  //  signalSystematics2011(searchRegions_[0]);
  //    return;

  TString sampleOfInterest="LM9";

  vector<ofstream*> textfiles;
  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    char effoutput[500];
    sprintf(effoutput,"signalSyst.%s.%s%s.dat",sampleOfInterest.Data(),searchRegions_[i].btagSelection.Data(),searchRegions_[i].owenId.Data());

    textfiles.push_back( new ofstream(effoutput));
  }

  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    (*textfiles[i])<<m0_<<" "<<m12_<<" ";
    SignalEffData SB =  signalSystematics2011(sbRegions_[i],false,false,sampleOfInterest);
    SignalEffData SIG=  signalSystematics2011(searchRegions_[i],false,false,sampleOfInterest);

    SignalEffData SBSL= signalSystematics2011(sbRegions_[i],true,false,sampleOfInterest);
    SignalEffData SIGSL=signalSystematics2011(searchRegions_[i],true,false,sampleOfInterest);

    SignalEffData SBLDP= signalSystematics2011(sbRegions_[i],false,true,sampleOfInterest);
    SignalEffData SIGLDP=signalSystematics2011(searchRegions_[i],false,true,sampleOfInterest);

    (*textfiles[i])<<SIG.rawYield<<" "<<SB.rawYield<<" "<<SIGSL.rawYield<<" "<<SBSL.rawYield<<" "<<SIGLDP.rawYield<<" "<<SBLDP.rawYield<<" "
		   <<SIG.effCorr<<" "<<SB.effCorr<<" "<<SIGSL.effCorr<<" "<<SBSL.effCorr<<" "<<SIGLDP.effCorr<<" "<<SBLDP.effCorr<<" "
	    <<SIG.totalSystematic<<" "<<SB.totalSystematic<<" "<<SIGSL.totalSystematic<<" "<<SBSL.totalSystematic<<" "<<SIGLDP.totalSystematic<<" "<<SBLDP.totalSystematic<<endl;
  }

  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    textfiles.at(i)->close();
  }

}

void runSystematics2011_mSugra() {
  TString sampleOfInterest="mSUGRAtanb40";
  loadSamples();
  clearSamples();
  addSample(sampleOfInterest);

  setSearchRegions();
  loadSusyScanHistograms();
  //  signalSystematics2011(searchRegions_[0]);
  //    return;


  vector<ofstream*> textfiles;
  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    char effoutput[500];
    sprintf(effoutput,"effCorrSyst.%s.%s%s.dat",sampleOfInterest.Data(),searchRegions_[i].btagSelection.Data(),searchRegions_[i].owenId.Data());

    textfiles.push_back( new ofstream(effoutput));
  }

  //loop over the scan points
  for (map<pair<int,int>, TH1D* >::iterator iscanpoint = scanProcessTotalsMap.begin(); iscanpoint!=scanProcessTotalsMap.end(); ++iscanpoint) {
  


    m0_=iscanpoint->first.first;
    m12_=iscanpoint->first.second;

    //need to check that there are 10k points, in order to throw out points that didn't work
    int nentries=    scanProcessTotalsMap[iscanpoint->first]->GetEntries();
    cout<<m0_<<" "<<m12_<<" NEntries = "<<nentries<<endl;
//     if (nentries != 10000) {
//       cout<<m0_<<" "<<m12_<<" being skipped! NEntries = "<<nentries<<endl;
//       continue;
//     }

    if (nentries==0) continue;

    for (unsigned int i=0; i<searchRegions_.size(); i++) {
      (*textfiles[i])<<m0_<<" "<<m12_<<" ";
//the logic here is screwy...i've got m0 and m12 as globals but then signalSystematics2011() takes them as arguments
//and overwrites the globals with the arguments.
//for now let's just live with it.
  
      SignalEffData SB =  signalSystematics2011(sbRegions_[i],false,false,sampleOfInterest,m0_,m12_);
      SignalEffData SIG=  signalSystematics2011(searchRegions_[i],false,false,sampleOfInterest,m0_,m12_);
      
      SignalEffData SBSL= signalSystematics2011(sbRegions_[i],true,false,sampleOfInterest,m0_,m12_);
      SignalEffData SIGSL=signalSystematics2011(searchRegions_[i],true,false,sampleOfInterest,m0_,m12_);
      
      SignalEffData SBLDP= signalSystematics2011(sbRegions_[i],false,true,sampleOfInterest,m0_,m12_);
      SignalEffData SIGLDP=signalSystematics2011(searchRegions_[i],false,true,sampleOfInterest,m0_,m12_);
      
      (*textfiles[i])<<SIG.effCorr<<" "<<SB.effCorr<<" "<<SIGSL.effCorr<<" "<<SBSL.effCorr<<" "<<SIGLDP.effCorr<<" "<<SBLDP.effCorr<<" "
		     <<SIG.totalSystematic<<" "<<SB.totalSystematic<<" "<<SIGSL.totalSystematic<<" "<<SBSL.totalSystematic<<" "<<SIGLDP.totalSystematic<<" "<<SBLDP.totalSystematic<<endl;
    }
  }

  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    textfiles.at(i)->close();
  }

}

void averageZ() {
  cout<<"Loose SIG ge1b"<<endl;
  double zn[]={22.8,17.9};
  double ze[]={11.6,10.6};
  jmt::weightedMean(2,zn,ze);

  cout<<"Tight SIG ge1b"<<endl;
  double zn2[]={3.4,3.6};
  double ze2[]={2.2,2.3};
  jmt::weightedMean(2,zn2,ze2);



}

std::pair<double,double> anotherABCD( const SearchRegion & region, bool datamode=false, float subscale=1,float SBshift=0, const TString LSBbsel="==0") {
  //kind of awful, but we'll return the estimate for datamode=true but the closure test bias for datamode=false

  /*
.L drawReducedTrees.C++
  */

  TString btagselection = region.btagSelection;

  TString owenKey = btagselection;
  owenKey += region.owenId;
  OwenData * myOwen = &(owenMap_[owenKey]);
  const  bool isSIG = region.isSIG;

  setStackMode(false);
  doData(datamode);
  setQuiet(true);

  useFlavorHistoryWeights_=false;
  loadSamples(); //needs to come first! this is why this should be turned into a class, with this guy in the ctor
  setColorScheme("nostack");
//   clearSamples();
//   addSample("PythiaPUQCD");

  TString sampleOfInterest = "PythiaPUQCD";
  if (datamode) {
    sampleOfInterest = "data";//"PythiaPUQCD";
    clearSamples();
    addSample("TTbarJets");
    addSample("WJets");
    addSample("ZJets");
    addSample("Zinvisible");
    addSample("VV");
    addSample("SingleTop");
  }
  else {
    resetSamples();
  }

  savePlots_=false;

  setLogY(false);
  TString var,xtitle;
  int nbins;
  float low,high;

 // --  count events
  TCut triggerCutLSB = "1";
  TCut triggerCut = "1";
  if (datamode) {
    triggerCutLSB = "pass_utilityHLT_HT300>=1";
    triggerCut = "cutTrigger==1";
  }


  TCut HTcut=region.htSelection.Data(); 
  TCut SRMET = region.metSelection.Data();
  TCut ge1b = "nbjetsSSVHPT>=1";
  if (btagselection=="ge2b") {
    ge1b="nbjetsSSVHPT>=2";
  }
  else if (btagselection=="ge1b") {}
  else if (btagselection=="ge3b") {
    ge1b="nbjetsSSVHPT>=3";
  }
  else {assert(0);}

  SRMET = SRMET && ge1b && triggerCut; //apply physics trigger and physics b cut in high MET region

  TCut baseline = "cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1";
  baseline = baseline&&HTcut;
  TCut cleaning = "weight<1000";

  char cutstring1[100];
  sprintf(cutstring1,"MET>= %.0f && MET < %.0f && nbjetsSSVHPT%s", 50+SBshift,100+SBshift,LSBbsel.Data());
  //  cout<<"*** SB cut is "<<cutstring1<<endl<<endl;
  TCut SBMET = TCut(cutstring1)&&triggerCutLSB;
  TCut dpcut = "1";//"minDeltaPhiN>=4";
  //  TCut passOther = "deltaPhiMPTcaloMET<2";
  //  TCut failOther = "deltaPhiMPTcaloMET>=2";
  TCut passOther = "minDeltaPhiN>=4";
  TCut failOther = "minDeltaPhiN<4";

  double A,B,D,SIG,Aerr,Berr,Derr,SIGerr;
  //50 - 100 and high MPT,MET ("A")
  selection_ = baseline && cleaning && dpcut  && SBMET && failOther; //auto cast to TString seems to work

  var="HT"; xtitle=var;
  nbins=10; low=0; high=5000;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_A");
  A=getIntegral(sampleOfInterest);
  Aerr=getIntegralErr(sampleOfInterest);
  //B
  selection_ = baseline && cleaning && dpcut  && SBMET && passOther; //auto cast to TString seems to work
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_B");
  B=getIntegral(sampleOfInterest);
  Berr=getIntegralErr(sampleOfInterest);
  //D
  selection_ = baseline && cleaning && dpcut  && SRMET && failOther; //auto cast to TString seems to work
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_D");
  D=getIntegral(sampleOfInterest);
  Derr=getIntegralErr(sampleOfInterest);

  double Dsub = 0,Dsuberr=0;
  if (datamode) {
    Dsub = subscale*getIntegral("totalsm");
    Dsuberr = subscale*getIntegralErr("totalsm");

    //special stuff for owen
    myOwen->Nlsb_0b_ldp = A;
    myOwen->Nlsb_0b     = B;
    if (isSIG) myOwen->Nsig_ldp = D;
    else       myOwen->Nsb_ldp = D;

    if (isSIG) {
      myOwen->Nttbarmc_sig_ldp = getIntegral("TTbarJets");
      myOwen->Nsingletopmc_sig_ldp = getIntegral("SingleTop");
      double lsfw = getIntegral("WJets") > 0 ?  pow(getIntegralErr("WJets"),2)/getIntegral("WJets"): -1;
      double lsfzj = getIntegral("ZJets") > 0 ?  pow(getIntegralErr("ZJets"),2)/getIntegral("ZJets"): -1;
      double lsfz = getIntegral("Zinvisible") > 0 ?  pow(getIntegralErr("Zinvisible"),2)/getIntegral("Zinvisible"): -1;

      if (lsfw>0) myOwen->lsf_WJmc=lsfw;
      if (lsfz>0) myOwen->lsf_Znnmc=lsfz;
      if (lsfzj>0) myOwen->lsf_Zjmc=lsfzj;

      myOwen->NWJmc_sig_ldp = getIntegral("WJets")/lsfw;
      myOwen->NZnnmc_sig_ldp = getIntegral("Zinvisible")/lsfz;
      myOwen->NZjmc_sig_ldp = getIntegral("ZJets")/lsfzj;

    }
    else {
      myOwen->Nttbarmc_sb_ldp = getIntegral("TTbarJets");
      myOwen->Nsingletopmc_sb_ldp = getIntegral("SingleTop");

      double lsfw = getIntegral("WJets") > 0 ?  pow(getIntegralErr("WJets"),2)/getIntegral("WJets"): -1;
      double lsfz = getIntegral("Zinvisible") > 0 ?  pow(getIntegralErr("Zinvisible"),2)/getIntegral("Zinvisible"): -1;
      double lsfzj = getIntegral("ZJets") > 0 ?  pow(getIntegralErr("ZJets"),2)/getIntegral("ZJets"): -1;

      if (lsfw>0) myOwen->lsf_WJmc=lsfw;
      if (lsfz>0) myOwen->lsf_Znnmc=lsfz;
      if (lsfzj>0) myOwen->lsf_Zjmc=lsfzj;

      myOwen->NWJmc_sb_ldp = getIntegral("WJets")/lsfw;
      myOwen->NZnnmc_sb_ldp = getIntegral("Zinvisible")/lsfz;
      myOwen->NZjmc_sb_ldp = getIntegral("ZJets")/lsfzj;
    }
    //end of special stuff for owen
  }
  if (Dsub>D) {Dsub=D-0.00001; cout<<"Subtraction in D is as big as D!"<<endl;}

  //SIG
  selection_ = baseline && cleaning && dpcut  && SRMET && passOther; //auto cast to TString seems to work
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_C");
  SIG=getIntegral(sampleOfInterest);
  SIGerr=getIntegralErr(sampleOfInterest);

  if (datamode) { //for owen
    if (isSIG)  myOwen->Nsig = SIG;
    else        myOwen->Nsb  = SIG;
  }

  //now calculate B*D/A
  double numerr=jmt::errAtimesB(B,Berr,D-Dsub,sqrt(Derr*Derr+Dsuberr*Dsuberr));
  double num = B*(D-Dsub);
  double estimate = num / A;
  double estimateerr= jmt::errAoverB(num,numerr,A,Aerr);
//   cout<<" ==== "<<endl
//       <<"Estimate = "<<estimate<<" +/- "<<estimateerr<<endl
//       <<" truth   = "<<SIG     <<" +/- "<<SIGerr<<endl;
  btagselection += region.owenId;
  btagselection += isSIG ? "SIG":"SB";
  char output[500];
  if (!datamode) {
  sprintf(output,"%s & %s & %s & %s & %s & %s \\\\ %% %f",btagselection.Data(),
	  jmt::format_nevents(B,Berr).Data(),jmt::format_nevents(A,Aerr).Data(),
	  jmt::format_nevents(D,Derr).Data(),jmt::format_nevents(estimate,estimateerr).Data(),
	  jmt::format_nevents(SIG,SIGerr).Data(),100*(SIG-estimate)/SIG);}
  else {
    sprintf(output,"%s & %d & %d & %d & %s & %s   \\\\",btagselection.Data(),
	    TMath::Nint(B),TMath::Nint(A),
	    TMath::Nint(D), jmt::format_nevents(Dsub,Dsuberr).Data(),
	    jmt::format_nevents(estimate,estimateerr).Data());
    cout<<"(qcd) DATA\t";
  }
  cout<<output<<endl;
  if (datamode)  return make_pair(estimate,estimateerr);
  return make_pair(100*(SIG-estimate)/SIG, 0);
}

void runClosureTest2011(std::map<TString, std::vector<double> > & syst)  {

  setSearchRegions();
  assert(sbRegions_.size() == searchRegions_.size());

  for (unsigned int i=0; i<sbRegions_.size(); i++) {
    syst["Closure"].push_back( anotherABCD(sbRegions_[i]).first);
    syst["Closure"].push_back( anotherABCD(searchRegions_[i]).first);
  }

}

void runClosureTest2011()  {
  std::map<TString, std::vector<double> > dummy;
  runClosureTest2011(dummy);

}


//to hold systematics
//making this global because I'm lazy...
std::map<TString, std::vector<double> > qcdSystErrors;
void runDataQCD2011(const bool forOwen=false) {

  setSearchRegions();
  assert(sbRegions_.size() == searchRegions_.size());

  //data structures to hold results
  vector<std::pair<double,double> > n;
  vector<std::pair<double,double> > subp;
  vector<std::pair<double,double> > subm;

  vector<std::pair<double,double> > sbp;
  vector<std::pair<double,double> > sbm;

  cout<<" ==== Nominal data results === "<<endl;
  for (unsigned int i=0; i<sbRegions_.size(); i++) {
    n.push_back( anotherABCD(sbRegions_[i],true));
    n.push_back( anotherABCD(searchRegions_[i],true));
  }
  cout<<" =END Nominal data results === "<<endl;

  /*
    the "owen mode" is collecting in global data structures the quantities that owen needs
    If we run the systematics then we will clobber the regular numbers. So quit now with just the nominal
    values.
  */
  if (forOwen) return;
  
  //now do it again with +50% subtraction
  cout<<" subtraction +50% "<<endl;
  for (unsigned int i=0; i<sbRegions_.size(); i++) {
    subp.push_back( anotherABCD(sbRegions_[i],true,1.5));
    subp.push_back( anotherABCD(searchRegions_[i],true,1.5));
  }

  //now do it again with -50% subtraction
  cout<<" subtraction -50% "<<endl;
  for (unsigned int i=0; i<sbRegions_.size(); i++) {
    subm.push_back( anotherABCD(sbRegions_[i],true,0.5));
    subm.push_back( anotherABCD(searchRegions_[i],true,0.5));
  }

  cout<<" ==== systematics for MC subtraction ==== "<<endl;
  for (unsigned int j=0; j<n.size(); j++) {
    double var1 = 100*(n[j].first  -subp[j].first)/n[j].first;
    double var2 = 100*(n[j].first -subm[j].first)/n[j].first;
    //expect these to be equal and opposite
    if ( fabs(var1)-fabs(var2) > 0.1 ) cout<<j<<" found MC subtraction systematic size discrepancy"<<endl;
    cout<<var1<<"\t"<<var2<<endl;

    qcdSystErrors["MCsub"].push_back( fabs(var1)>fabs(var2) ? fabs(var1) : fabs(var2));
  }
  cout<<" =END systematics for MC subtraction ==== "<<endl;


  //do it with shifted SB
  cout<<" SB +10 GeV"<<endl;
  for (unsigned int i=0; i<sbRegions_.size(); i++) {
    sbp.push_back( anotherABCD(sbRegions_[i],true,1,10));
    sbp.push_back( anotherABCD(searchRegions_[i],true,1,10));
  }

  //now do it again with -50% subtraction
  cout<<" SB -10 GeV "<<endl;
  for (unsigned int i=0; i<sbRegions_.size(); i++) {
    sbm.push_back( anotherABCD(sbRegions_[i],true,1,-10));
    sbm.push_back( anotherABCD(searchRegions_[i],true,1,-10));
  }

  cout<<" ==== systematics for SB shift ==== "<<endl;
  for (unsigned int j=0; j<n.size(); j++) {
    double var1 = 100*(n[j].first  -sbp[j].first)/n[j].first;
    double var2 = 100*(n[j].first -sbm[j].first)/n[j].first;
    cout<<var1<<"\t"<<var2<<endl;
    //expect these to be different...no sanity check is warranted
    qcdSystErrors["SBshift"].push_back( fabs(var1)>fabs(var2) ? fabs(var1) : fabs(var2));
  }
  cout<<" =END systematics for SB shift ==== "<<endl;

  cout<<"== Cross check with >=1 b instead of exactly 0 b =="<<endl;
  for (unsigned int i=0; i<sbRegions_.size(); i++) {
    anotherABCD(sbRegions_[i],true,1,0,">=1");
    anotherABCD(searchRegions_[i],true,1,0,">=1");
  }

  cout<<" == Running QCD closure test =="<<endl;
  runClosureTest2011(qcdSystErrors);

  cout<<" == QCD systematics summary =="<<endl;
  for (unsigned int j=0; j<n.size(); j++) {
    qcdSystErrors["Total"].push_back( sqrt( pow(qcdSystErrors["MCsub"].at(j),2) +  pow(qcdSystErrors["Closure"].at(j),2)+ pow(qcdSystErrors["SBshift"].at(j),2)));
    cout<<j<<"\t&"<<qcdSystErrors["MCsub"].at(j)<<" & "<<qcdSystErrors["Closure"].at(j)<<" & "<<qcdSystErrors["SBshift"].at(j)<<" & "<<qcdSystErrors["Total"].at(j)<<endl;
  }
  
}

//i don't like passing the index instead of the region itself, but it makes some things easier.
//this code is all around a big kludge...
double slABCD(const unsigned int searchRegionIndex, bool datamode=false, const TString & mode="" ) {
  //in datamode, return the estimate; in non-datamode, return the Closure Test results (true-pred)/true

  /*
.L drawReducedTrees.C++
  */

  const SearchRegion region = searchRegions_[searchRegionIndex];
  const SearchRegion qcdsubregion = sbRegions_[searchRegionIndex];

  TString btagselection = region.btagSelection;
  TString owenKey = btagselection;
  owenKey += region.owenId;
  OwenData * myOwen = &(owenMap_[owenKey]);

  //we have to be careful. because of all of the global variables used in this code, we have to call
  //this other guy before we do anything else (should make it a class!)

  double SBsubQCD =0,SBsubQCDerr=0,SBsubZ=0,SBsubZerr=0,SBsubMisc=0,SBsubMiscerr=0;
  std::pair<double,double> SBsubQCDp;
  if (datamode)  {
    //regions should only differ in the MET selection
    assert ( region.btagSelection == qcdsubregion.btagSelection);
    assert(region.htSelection == qcdsubregion.htSelection);
    //get the QCD estimate for the SB region
    SBsubQCDp = anotherABCD(qcdsubregion, datamode, 1 ) ;
    SBsubQCD=SBsubQCDp.first;
    SBsubQCDerr=SBsubQCDp.second;

    if (mode.Contains("QCD")) {
       //SBsubQCDerr is the statistical error on the QCD in SB.
      //but the systematic is needed too
      //factor 0.01 to convert from human-readable % to fractional
      double qcdsbsyst=0.01*qcdSystErrors["Total"].at(2*searchRegionIndex);
      cout<<"Found a QCD systematic of "<<100*qcdsbsyst<<endl;
//       if (tight) {
// 	if ( btagselection=="ge1b")     qcdsbsyst=0.37;
// 	else if (btagselection=="ge2b") qcdsbsyst=0.54;
//       }
//       else {
// 	if ( btagselection=="ge1b")     qcdsbsyst=0.16;
// 	else if (btagselection=="ge2b") qcdsbsyst=0.20;
//       }
      SBsubQCDerr = sqrt( SBsubQCDerr*SBsubQCDerr + pow(qcdsbsyst*SBsubQCD,2));
      
      if (mode=="QCDup") SBsubQCD += SBsubQCDerr;
      if (mode=="QCDdown") SBsubQCD -= SBsubQCDerr;
    }

    //now hard-coded Z->nunu
    double zv[2];
    double ze[2];
 //need to average mu mu and ee estimates
    double zsbsyst=0.5;
    if ( qcdsubregion.owenId == "Tight" && btagselection=="ge1b") { //i think this should work for now...
      zv[0] = 7.0; ze[0]=2.9; //Tight SB ge1b mumu
      zv[1] = 7.3; ze[1]=3.4; //Tight SB ge1b ee
      zsbsyst = 0.55;
    }
    else    if ( qcdsubregion.owenId == "Loose"&&btagselection=="ge1b") {//    else {
      zv[0] = 16.0; ze[0]=11.4; //Loose SB ge1b mumu
      zv[1] = 30.2; ze[1]=17.8; //Loose SB ge1b ee
      zsbsyst=0.2;
    }
    else if (qcdsubregion.owenId == "Tight" && btagselection=="ge2b") {
      zv[0] = 1.6; ze[0]=0.8; //Tight SB ge2b mumu
      zv[1] = 0; ze[1]=0.8; //Tight SB ge2b ee (error is bullshit)
      zsbsyst=1.03;
    }
    else if (qcdsubregion.owenId == "Loose" && btagselection=="ge2b") {
      zv[0] = 6.7; ze[0]=3.5; //Loose SB ge2b mumu
      zv[1] = 0; ze[1]=3.5; //Loose SB ge2b ee (error is bullshit)
      zsbsyst = 0.55;
    }
    else {assert(0);}
    SBsubZ = jmt::weightedMean(2,zv,ze);
    SBsubZerr = jmt::weightedMean(2,zv,ze,true);

    if (mode.Contains("Z")) {
      SBsubZerr = sqrt( SBsubZerr*SBsubZerr + pow(zsbsyst*SBsubZ,2));
      
      if (mode=="Zup") SBsubZ += SBsubZerr;
      if (mode=="Zdown") SBsubZ -= SBsubZerr;
    }
  }

  setStackMode(false);
  doData(true);
  setQuiet(true);

  TString sampleOfInterest = "totalsm";//"PythiaPUQCD";
  useFlavorHistoryWeights_=false;
  loadSamples(); //needs to come first! this is why this should be turned into a class, with this guy in the ctor
  setColorScheme("nostack");
  clearSamples();
  if (datamode) {
    //    addSample("SingleTop");
    addSample("ZJets");
    addSample("VV");
    sampleOfInterest="data";
  }
  else {
    addSample("TTbarJets");
//     addSample("WJets");
//     addSample("SingleTop");
  }

  //  TString sampleOfInterest = "PythiaPUQCD";

  savePlots_=false;

  setLogY(false);
  TString var,xtitle;
  int nbins;
  float low,high;

 // --  count events
  TCut ge1b = "nbjetsSSVHPT>=1";
  if (btagselection=="ge2b") {
    ge1b="nbjetsSSVHPT>=2";
  }
  else if (btagselection=="ge1b") {}
  else if (btagselection=="ge3b") {
    ge1b="nbjetsSSVHPT>=3";
  }
  else {assert(0);}

  TCut HTcut=region.htSelection.Data(); 
  TCut SRMET = region.metSelection.Data();
//   if (tight) {
//     HTcut = "HT>=500";
//     SRMET = "MET >= 300";
//   }
//   else {
//     HTcut = "HT>=350";
//     SRMET = "MET >= 200";
//   }

  TCut baseline = "cutPV==1 && cut3Jets==1 && cutTrigger==1";
  baseline = baseline&&HTcut&&ge1b;
  TCut cleaning = "weight<1000";
  TCut SBMET = qcdsubregion.metSelection.Data();//"MET>=150 && MET<200";
  TCut dpcut = "minDeltaPhiN>=4";
  //  TCut passOther = "deltaPhiMPTcaloMET<2";
  //  TCut failOther = "deltaPhiMPTcaloMET>=2";
  TCut failOther = "((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0))";
  TCut passOther = "nElectrons==0 && nMuons==0";

  double A,B,D,SIG,Aerr,Berr,Derr,SIGerr;
  //  double dA,dB,dD,dSIG,dAerr,dBerr,dDerr,dSIGerr;
  //A = SB, SL
  selection_ = baseline && cleaning && dpcut  && SBMET && failOther; //auto cast to TString seems to work
  var="HT"; xtitle=var;
  nbins=10; low=0; high=5000;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_A");
  A=getIntegral(sampleOfInterest);
  Aerr=getIntegralErr(sampleOfInterest);
  //B = SB
  selection_ = baseline && cleaning && dpcut  && SBMET && passOther; //auto cast to TString seems to work
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_B");
  B=getIntegral(sampleOfInterest);
  Berr=getIntegralErr(sampleOfInterest);
  if (datamode) {
    SBsubMisc=getIntegral("totalsm");
    SBsubMiscerr=getIntegralErr("totalsm");

    if (mode.Contains("MC")) {
      const double mcsyst = 1; //100% uncertainty
      SBsubMiscerr = sqrt( SBsubMiscerr*SBsubMiscerr + pow(mcsyst*SBsubMisc,2));

      if (mode=="MCup") SBsubMisc += SBsubMiscerr;
      else if (mode=="MCdown") SBsubMisc -= SBsubMiscerr;
    }
  }
  //D = SIG,SL
  selection_ = baseline && cleaning && dpcut  && SRMET && failOther; //auto cast to TString seems to work
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_D");
  D=getIntegral(sampleOfInterest);
  Derr=getIntegralErr(sampleOfInterest);

  if (datamode) { //for owen
    myOwen->Nsig_sl = D;
    myOwen->Nsb_sl = A;
  }

  //SIG
  selection_ = baseline && cleaning && dpcut  && SRMET && passOther; //auto cast to TString seems to work
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_C");
  SIG=getIntegral(sampleOfInterest);
  SIGerr=getIntegralErr(sampleOfInterest);

  //now calculate B*D/A
  double suberr = sqrt(SBsubMiscerr*SBsubMiscerr + SBsubQCDerr*SBsubQCDerr + SBsubZerr*SBsubZerr);

  double numerr=jmt::errAtimesB(B-SBsubMisc-SBsubQCD-SBsubZ,sqrt(suberr*suberr + Berr*Berr),
				D,Derr);
  double num = (B-SBsubMisc-SBsubQCD-SBsubZ)*D;
  double estimate = num / A;
  double estimateerr= jmt::errAoverB(num,numerr,A,Aerr);


//   cout<<" ==== "<<endl
//       <<"Estimate = "<<estimate<<" +/- "<<estimateerr<<endl
//       <<" truth   = "<<SIG     <<" +/- "<<SIGerr<<endl;
  // btagselection += tight ? " Tight " : " Loose ";
  btagselection += region.owenId;

  char output[500];
  if (!datamode) {
    sprintf(output,"%s & %s & %s & %s & %s & %s \\\\ %% %f",btagselection.Data(),
	    jmt::format_nevents(D,Derr).Data(),jmt::format_nevents(A,Aerr).Data(),
	    jmt::format_nevents(B,Berr).Data(),jmt::format_nevents(estimate,estimateerr).Data(),
	    jmt::format_nevents(SIG,SIGerr).Data(),100*(estimate-SIG)/SIG);
  }
  else {
    sprintf(output,"ttbar DATA %s & %d & %d & %d & %s & %s  \\\\",btagselection.Data(),
	    TMath::Nint(D),TMath::Nint(A),
	    TMath::Nint(B), jmt::format_nevents(SBsubMisc+SBsubQCD+SBsubZ,suberr).Data(),
	    jmt::format_nevents(estimate,estimateerr).Data());
  }
  cout<<output<<endl;

  if (datamode) return estimate;
  return 100*(estimate-SIG)/SIG;
}

void runSLClosureTest2011() {

  setSearchRegions();

  for (unsigned int j=0; j<searchRegions_.size();j++) slABCD(j);

}

void runTtbarEstimate2011(const bool forOwen=false) {
  setSearchRegions();

  double ttw[4];
  double p[4];
  double m[4];

  double qcd[4];
  double znn[4];
  double mc[4];
  //  int i=0;
  cout<<"nominal ttbar"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) ttw[j]=slABCD(j,true);
 
  if (forOwen) return;

  // for systematics due to QCD subtraction
  //must do   runDataQCD2011(); first (in order to fill syst errors for qcd)
  if (  qcdSystErrors.size()==0) {cout<<"Need to do runDataQCD2011() first!"<<endl; return;}

  cout<<"vary QCD subtraction"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) {
    p[j]=slABCD(j,true,"QCDup");
    m[j]=slABCD(j,true,"QCDdown");
    if (fabs(p[j]-ttw[j])/ttw[j] - fabs(m[j]-ttw[j])/ttw[j] > 0.01) cout<<"** Difference in size of QCD sub syst!"<<endl;
    qcd[j] = 100*fabs(p[j]-ttw[j])/ttw[j];
  }

  //for Znunu subtraction systematics
  cout<<"vary Znn subtraction"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) {
    p[j]=slABCD(j,true,"Zup");
    m[j]=slABCD(j,true,"Zdown");
    if (fabs(p[j]-ttw[j])/ttw[j] - fabs(m[j]-ttw[j])/ttw[j] > 0.01) cout<<"** Difference in size of Znn sub syst!"<<endl;
    znn[j] = 100*fabs(p[j]-ttw[j])/ttw[j];
  }

  //for MC subtraction systematics
  cout<<"vary MC subtraction"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) {
    p[j]=slABCD(j,true,"MCup");
    m[j]=slABCD(j,true,"MCdown");
    if (fabs(p[j]-ttw[j])/ttw[j] - fabs(m[j]-ttw[j])/ttw[j] > 0.01) cout<<"** Difference in size of MC sub syst!"<<endl;
    mc[j] = 100*fabs(p[j]-ttw[j])/ttw[j];
  }

  //finally, run the closure test
  double closure[4];
  cout<<"Running ttbar closure test"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++)   closure[j]=fabs(slABCD(j));

  cout<<" == summary (%) == "<<endl;
  cout<<"\tClosure\tQCD\tZ\tMC\tTotal"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) {
    double totalsyst = sqrt(closure[j]*closure[j] + qcd[j]*qcd[j] + znn[j]*znn[j] + mc[j]*mc[j]);
    cout<<j<<"\t"<<closure[j]<<"\t"<<qcd[j]<<"\t"<<znn[j]<<"\t"<<mc[j]<<"\t"<<totalsyst <<endl;
  }


}

void printOwenAll() {

  runDataQCD2011(true);
  runTtbarEstimate2011(true);

  for (std::map<TString,OwenData>::iterator i=owenMap_.begin(); i!=owenMap_.end(); ++i) {
    printOwen( i->first);
  }
  cout<<"luminosity = "<<lumiScale_<<endl;

}

void AN2011_prescale( TString btagselection="ge1b" ) {
  /*
.L drawReducedTrees.C++
  */
  TCut btagcut = "nbjetsSSVHPT>=1";
  if ( btagselection=="ge1b") {} //do nothing
  else  if ( btagselection=="ge2b" ) {
    btagcut = "nbjetsSSVHPT>=2";
  }
  else if ( btagselection=="eq1b" ) {
    btagcut = "nbjetsSSVHPT==1";
  }
  else if ( btagselection=="eq0b" ) {
    btagcut = "nbjetsSSVHPT==0";
  }
  else {
    assert(0);
  }

  loadSamples();

  int nbins;
  float low,high;
  TString var,xtitle;
  
  doOverflowAddition(true);

  setQuiet(false);

  // ==========================

  TCut util = "pass_utilityHLT_HT300>=1 && weight<1000";
  lumiScale_ = 19.23288; //customize lumiScale_

  // ========= regular N-1 plots

  resetSamples(); //use all samples
  setStackMode(true); //regular stack
  setColorScheme("stack");
  doData(true);
  drawMCErrors_=true;

  //MET
  selection_ =TCut("cutHT==1 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&util&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 20; low=0; high=200;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_MET_"+btagselection);
  setLogY(true);  setPlotMinimum(1e-1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_MET_"+btagselection);

  //njets in 50 < MET < 100 region
  selection_ =TCut("cutHT==1 && cutPV==1 && MET>=50 && MET<100 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&util&&btagcut;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 8; low=1; high=9;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_njets_LSB_"+btagselection);

  selection_ =TCut("cutHT==1 && cutPV==1 && MET>=50 && MET<100 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && weight<10000")&&util&&btagcut;
  var="minDeltaPhiN"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 20; low=0; high=40;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_minDeltaPhiN_LSB_"+btagselection);
  setLogY(true);  setPlotMinimum(1e-1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_minDeltaPhiN_LSB_"+btagselection);


  //look again at MET
  selection_ =TCut("MET>=50 && MET<100 && cutHT==1 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&util&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 10; low=50; high=100;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_METlimited_"+btagselection);


  //sanity check
  doRatioPlot(true);
  ratioMin=0; ratioMax=3;
  selection_ =TCut("MET>150 && cutHT==1 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN < 4")&&util;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 10; low=150; high=250;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_METhigh_lowDP");

  lumiScale_ =  1091.891;
  ratioMin=1; ratioMax=1.5;
  selection_ ="cutTrigger==1 && MET>150 && cutHT==1 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN < 4";
  drawPlots(var,nbins,low,high,xtitle,"Events", "unprescaled_METhigh_lowDP");


  
}


void AN2011_ttbarw( TString btagselection="ge1b", TString HTselection="Loose" , TString samplename="TTbarJets") {

  /*
goal:

show MET distributions for w+jets, ttbar
in both the SL and nominal samples

this is going to require some serious kludgey stuff because this code has
never been used to plot the sample sample with different cuts on top of each 
other.
  */

  /*
.L drawReducedTrees.C++
  */

  TCut btagcut = "nbjetsSSVHPT>=1";
  if ( btagselection=="ge1b") {} //do nothing
  else  if ( btagselection=="ge2b" ) {
    btagcut = "nbjetsSSVHPT>=2";
  }
  else if ( btagselection=="eq1b" ) {
    btagcut = "nbjetsSSVHPT==1";
  }
  else if ( btagselection=="ge0b" ) {
    btagcut = "nbjetsSSVHPT>=0";
  }
  else {
    assert(0);
  }


  TCut HTcut="HT>=350";
  if (HTselection=="Tight")  HTcut="HT>=500";

  loadSamples();

  int nbins;
  float low,high;
  TString var,xtitle;
  
  doOverflowAddition(true);

  setStackMode(false,true); //normalized
  setColorScheme("nostack");
  clearSamples();
  addSample(samplename);
  drawLegend(false);
  doRatioPlot(false);

  doData(false);

  TString sample;
  if(samplename=="TTbarJets") sample = "ttbar";
  else if(samplename=="WJets") sample = "wjets";

  selection_ =TCut("MET>=150 && cutPV==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4")&&btagcut&&HTcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 40; low=150; high=550;
  //nbins = 20; low=150; high=550;
  //nbins = 16; low=150; high=550;
  drawPlots(var,nbins,low,high,xtitle,"Arbitrary units", "SBandSIG_MET_SL_"+sample+"_"+btagselection+"_"+HTselection);

  int lowbin = hinteractive->FindBin(low);
  int boundarybin_sb = 0, boundarybin_sig=0;

  boundarybin_sb = hinteractive->FindBin(200);

  //MET Signal region depends on "Loose" or "Tight" selection
  if (HTselection=="Tight")  boundarybin_sig = hinteractive->FindBin(300);
  else boundarybin_sig = hinteractive->FindBin(200);

  int highbin = hinteractive->FindBin(high);

  //get the SIG/SB Ratio
  double sl_sb_err=0, sl_sig_err=0;
  double sl_sb = hinteractive->IntegralAndError(lowbin,boundarybin_sb-1,sl_sb_err);
  double sl_sig = hinteractive->IntegralAndError(boundarybin_sig,highbin-1,sl_sig_err);
  double r_sl_sigoversb = sl_sig/sl_sb;
  double r_sl_sigoversb_err = jmt::errAoverB(sl_sig,sl_sig_err,sl_sb,sl_sb_err);

  std::cout << "SL: SB content= " << sl_sb << " +/- " << sl_sb_err << std::endl;
  std::cout << "SL: SIG content= " << sl_sig << " +/- " << sl_sig_err << std::endl;
  std::cout << "SL: SIG/SB ratio = " << r_sl_sigoversb << " +/- " << r_sl_sigoversb_err << std::endl;


  TH1D* SLplot = (TH1D*)hinteractive->Clone("SLplot");
  //now switch to the normal selection
  selection_ =TCut("MET>=150 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&btagcut&&HTcut;
  drawPlots(var,nbins,low,high,xtitle,"Arbitrary units", "SBandSIG_MET_normal_"+sample+"_"+btagselection+"_"+HTselection);
  TH1D* SIGplot = (TH1D*)hinteractive->Clone("SIGplot");
  SLplot->SetLineColor(kRed);
  SLplot->SetMarkerColor(kRed);
  SLplot->Draw("SAME");
  
  thecanvas->SaveAs("METshape_"+sample+"_SLandStandard_"+btagselection+"_"+HTselection+".pdf");
  //not enough stats in WJets to really show anything

  //Hack to get it plotted with ratio plot
  TCanvas* myC = 0;
  myC = new TCanvas("myC", "myC", 600,700);
  myC->Divide(1,2);
  const float padding=0.01; const float ydivide=0.2;
  myC->GetPad(1)->SetPad( padding, ydivide + padding, 1-padding, 1-padding);
  myC->GetPad(2)->SetPad( padding, padding, 1-padding, ydivide-padding);
  myC->GetPad(1)->SetLogy(1);
  myC->GetPad(1)->SetRightMargin(.05);
  myC->GetPad(2)->SetRightMargin(.05);
  myC->GetPad(1)->Modified();
  myC->GetPad(2)->Modified();
  myC->cd(1);
  SIGplot->Draw();
  SLplot->Draw("SAME");
  TH1D* myRatio = new TH1D("ratio", "ratio", nbins,low,high);
  myRatio->Sumw2();
  myRatio->SetLineColor(kBlack);
  myRatio->SetMarkerColor(kBlack);
  myRatio->Divide(SIGplot,SLplot);
  myRatio->SetMinimum(0);
  myRatio->SetMaximum(2);
  myRatio->GetYaxis()->SetNdivisions(200 + int(ratioMax-ratioMin)+1);    //set ticks ; to be seen if this really works
  myRatio->GetYaxis()->SetLabelSize(0.15); //make y label bigger
  myC->cd(2);
  myRatio->Draw();
  TLine* myLine = new TLine(low, 1, high, 1);
  myLine->Draw();
  myC->cd(1);
  TLatex* mytext = new TLatex(3.570061,23.08044,"CMS Preliminary");
  mytext->SetNDC();
  mytext->SetTextAlign(13);
  mytext->SetX(0.6);
  mytext->SetY(0.9);
  mytext->SetTextFont(42);
  mytext->SetTextSizePixels(24);
  mytext->Draw();
  TLegend* myLegend = new TLegend(0.63,0.84,0.81,0.7);
  myLegend->SetFillColor(0);
  myLegend->SetBorderSize(0);
  myLegend->SetLineStyle(0);
  myLegend->SetTextFont(42);
  myLegend->SetFillStyle(0);
  myLegend->SetTextSize(0.04);
  myLegend->AddEntry(SIGplot,"TTbar, SIG", "l");
  myLegend->AddEntry(SLplot, "TTbar, SL", "l");
  myLegend->Draw();
  myC->Print("METshape_logAndRatio_"+sample+"_SLandStandard_"+btagselection+"_"+HTselection+".pdf");

  //get the SIG/SB Ratio
  double lv_sb_err=0, lv_sig_err=0;
  double lv_sb = hinteractive->IntegralAndError(lowbin,boundarybin_sb-1,lv_sb_err);
  double lv_sig = hinteractive->IntegralAndError(boundarybin_sig,highbin-1,lv_sig_err);
  double r_lv_sigoversb = lv_sig/lv_sb;
  double r_lv_sigoversb_err = jmt::errAoverB(lv_sig,lv_sig_err,lv_sb,lv_sb_err);

  std::cout << "LV: SB content= " << lv_sb << " +/- " << lv_sb_err << std::endl;
  std::cout << "LV: SIG content= " << lv_sig << " +/- " << lv_sig_err << std::endl;
  std::cout << "LV: SIG/SB ratio = " << r_lv_sigoversb << " +/- " << r_lv_sigoversb_err << std::endl;


}

void AN2011_r() {
  TString btagselection="ge1b";
  TCut btagcut = "nbjetsSSVHPT>=1";
  loadSamples();

  // == draw r(MET)
  clearSamples();
  addSample("PythiaPUQCD");
  setColorScheme("nostack");

  drawLegend(false);
  doRatioPlot(false);
  doData(false);
  setStackMode(false);
  doOverflowAddition(true);

  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1&&weight<1000")&&btagcut;
  setPlotMaximum(0.7); setPlotMinimum(0);
  drawR("minDeltaPhi",0.3,15,50,350,"old_"+btagselection);
  drawR("minDeltaPhiN",4,15,50,350,btagselection);

  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjetsSSVHPT==0&&weight<1000";
  drawR("minDeltaPhiN",4,15,50,350,"eq0b");

  //now try drawing it in data
  //careful -- there is a bug in the drawR() implementation
  //it doesn't work right unless you turn on doRatioPlot
  drawTotalSM_=true;
  doRatioPlot(true);
  drawLegend(true);
  setPlotMinimum(0); setPlotMaximum(0.5);
  lumiScale_ = 19.23288; //customize lumiScale_
  resetSamples();
  doData(true);
  TCut util = "pass_utilityHLT_HT300>=1 && weight<1000";
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjetsSSVHPT==0")&&util;
  drawR("minDeltaPhiN",4,10,50,150,"Data_eq0b");


}

void AN2011_opt() {

  TCut ge1b = "nbjetsSSVHPT>=1";
  TCut ge2b = "nbjetsSSVHPT>=2";
  TCut ge3b = "nbjetsSSVHPT>=3";

  lumiScale_=1800;

  int nbins;
  float low,high;
  TString var,xtitle;
  
  resetSamples(); //use all samples
  setStackMode(true); //regular stack
  setColorScheme("stack");
  doData(false);
  drawMCErrors_=true;

  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=400 && minDeltaPhiN >= 4")&&ge1b;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 10; low=400; high=1400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "opt_HT_MET400_HT400_ge1b");


  selection_ =TCut("HT>=500 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=300 && minDeltaPhiN >= 4")&&ge2b;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 10; low=400; high=1400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "opt_HT_MET300_HT500_ge2b");


  selection_ =TCut("HT>=900 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN >= 4")&&ge1b;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 10; low=400; high=1400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "opt_HT_MET200_HT900_ge1b");



  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN >= 4")&&ge3b;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 10; low=400; high=1400;
  drawPlots(var,nbins,low,high,xtitle,"Events", "opt_HT_MET200_HT400_ge3b");


  //


  selection_ =TCut("HT>=1000 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN >= 4")&&ge1b;
  var="MET"; xtitle="MET (GeV)";
  nbins = 10; low=200; high=600;
  drawPlots(var,nbins,low,high,xtitle,"Events", "opt_HT_MET200_HT1000_ge1b");

}

void drawMSugraTest() {
  /*
.L drawReducedTrees.C++
  */

  loadSamples();
  clearSamples();
  addSample("TTbarJets");
  addSample("mSUGRAtanb40");
  m0_=800;
  m12_=200;
  doOverflowAddition(true);
  setStackMode(true); //regular stack
  setColorScheme("stack");
  doData(false);
  int nbins; float low,high;
  TString var,xtitle;

  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN >= 4 && nbjets>=1");
  var="njets"; xtitle="Jet multiplicity";
  nbins = 8; low=1; high=9;
  drawPlots(var,nbins,low,high,xtitle,"Events", "testmsugra");

}

void AN2011( TString btagselection="ge1b" ) {
  /*
.L drawReducedTrees.C++
  */


  TCut btagcut = "nbjetsSSVHPT>=1";
  if ( btagselection=="ge1b") {} //do nothing
  else  if ( btagselection=="ge2b" ) {
    btagcut = "nbjetsSSVHPT>=2";
  }
  else if ( btagselection=="eq1b" ) {
    btagcut = "nbjetsSSVHPT==1";
  }
  else {
    assert(0);
  }

  loadSamples();

  int nbins;
  float low,high;
  TString var,xtitle;
  
  doOverflowAddition(true);

  setStackMode(false,true); //normalized
  setColorScheme("nostack");
  clearSamples();
  addSample("PythiaPUQCD");
  addSample("TTbarJets");
  addSample("LM9");

  doData(false);

  var="minDeltaPhiN"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 40; low=0; high=40;
  //no delta phi cut, loose MET window (only good for MC)
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100")&&btagcut;
  drawPlots(var,nbins,low,high,xtitle,"Arbitrary units", "minDeltaPhiN_looseMET_MConly_"+btagselection);


  // ========= regular N-1 plots

  resetSamples(); //use all samples
  setStackMode(true); //regular stack
  setColorScheme("stack");
  doData(true);
  drawMCErrors_=true;

  var="minDeltaPhiN"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 20; low=0; high=40;
  //no delta phi cut
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=150")&&btagcut;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_minDeltaPhiN_"+btagselection);

  // n Jets
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=150 && minDeltaPhiN >= 4")&&btagcut;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 8; low=1; high=9;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_njets_"+btagselection);

  //HT
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=150 && minDeltaPhiN >= 4")&&btagcut;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 20; low=350; high=1050;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_HT_"+btagselection);

  //MET
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 35; low=150; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_"+btagselection);

  //MET distribution with tighter HT cut
  selection_ =TCut("HT>500 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  if(btagselection=="ge2b") { nbins = 7; low=150; high=500; }//requested by 
  else{ nbins = 17; low=150; high=500;}
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_HT500_"+btagselection);


  // == finally, draw the signal region only!
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN >= 4")&&btagcut;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 20; low=350; high=1050;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_HT_"+btagselection);

  var="MET+HT"; xtitle="M_{eff} [GeV]";
  nbins = 20; low=550; high=1550;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_Meff__"+btagselection);

  selection_ =TCut("HT>=500 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=300 && minDeltaPhiN >= 4")&&btagcut;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 20; low=500; high=1050;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SIGtight_HTtight_"+btagselection);

  var="MET+HT"; xtitle="M_{eff} [GeV]";
  nbins = 10; low=800; high=1600;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SIGtight_Meff__"+btagselection);

  // ===== single lepton selection

  // == MET for the electron sample
  selection_ =TCut("MET>=150 &&cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && nElectrons==1 && nMuons==0 && minDeltaPhiN >= 4")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=150; high=450;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_1e0mu_"+btagselection);

  //pT?
  var="eleet1"; xtitle="electron p_{T} [GeV]";
  nbins = 20; low=0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_eleET_1e0mu_"+btagselection);

  // == MET for the muon sample
  selection_ =TCut("MET>=150 &&cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && nElectrons==0 && nMuons==1 && minDeltaPhiN >= 4")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=150; high=450;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_0e1mu_"+btagselection);

  var="muonpt1"; xtitle="muon p_{T} [GeV]";
  nbins = 20; low=0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_muonpT_0e1mu_"+btagselection);

  // == MET for the combined sample
  selection_ =TCut("MET>=150 &&cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=150; high=450;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_SL_"+btagselection);

  // == HT for the combined sample
  selection_ =TCut("MET>=150 &&cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4")&&btagcut;
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 20; low=350; high=1050;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_HT_SL_"+btagselection);

  // == njets for the combined sample
  selection_ =TCut("MET>=150 &&cutHT==1 && cutPV==1 && cutTrigger==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4")&&btagcut;
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 20; low=350; high=1050;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_njets_SL_"+btagselection);

  // == MET for the combined sample (HT>500)
  selection_ =TCut("MET>=150 && HT>500 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=150; high=450;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_SL_HT500_"+btagselection);

  //different scale to compare to Kristen
  selection_ =TCut("MET>=150 && HT>500 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 35; low=150; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_SL_HT500_morebins_"+btagselection);


  // == take a look at WJets
  selection_ ="MET>=150 && cutHT==1 && cutPV==1 && cutTrigger==1 && njets==2 && (nElectrons==0 && nMuons==1) && minDeltaPhiN >= 4 && nbjetsSSVHPT==0";
  var="MT_Wlep"; xtitle="M_{T} [GeV]";
  nbins = 20; low=0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MT_0e1mu_0b");

  //MET
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN < 4")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 30; low=150; high=450;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_ldp_"+btagselection);


}

void morerstuff() {
  /*
.L drawReducedTrees.C++
  */
  setStackMode(false);
  doData(false);

  useFlavorHistoryWeights_=false;

  loadSamples(); //needs to come first! this is why this should be turned into a class, with this guy in the ctor
  clearSamples();
  addSample("PythiaPUQCD");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
  setLogY(true);
  drawR("minDeltaPhi",0.3,13,70,230,"ge1b");
  setLogY(false);
  drawR("minDeltaPhi",0.3,13,70,230,"ge1b");

  //  TF1 *fa = new TF1("fa","[0]*exp(-1*[1]*x)",-3,3); 

}


void biasStudy() {
  /*
.L drawReducedTrees.C++
  */
  useFlavorHistoryWeights_=false;

  setStackMode(false);
  doData(true);

  //no met, no mindeltaphi
  setPlotMinimum(1e-2); //setPlotMaximum(1);

 drawTotalSM_=true;
 leg_y1=0.6;
 setLogY(true);
 
 setQuiet(true);

 double cbias,err;

 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 &&njets==3";
 drawR("minDeltaPhi",0.3,4,0,200,"ge1b-4bins");
 cbias= hinteractive->GetBinContent(4)/hinteractive->GetBinContent(3);
 err = jmt::errAoverB(hinteractive->GetBinContent(4),hinteractive->GetBinError(4),hinteractive->GetBinContent(3),hinteractive->GetBinError(3));
 cout<<"no extra cut cbias = "<<cbias <<" +/- "<<err<<endl;

 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 && (MET/HT<0.2) &&njets==3";
 drawR("minDeltaPhi",0.3,4,0,200,"ge1b-METoverHTcut-4bins");
 cbias= hinteractive->GetBinContent(4)/hinteractive->GetBinContent(3);
 err = jmt::errAoverB(hinteractive->GetBinContent(4),hinteractive->GetBinError(4),hinteractive->GetBinContent(3),hinteractive->GetBinError(3));
 cout<<"MET/HT<0.2   cbias = "<<cbias <<" +/- "<<err<<endl;

 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 && (MET/HT>=0.2) &&njets==3";
 drawR("minDeltaPhi",0.3,4,0,200,"ge1b-METoverHTcut-4bins");
 cbias= hinteractive->GetBinContent(4)/hinteractive->GetBinContent(3);
 err = jmt::errAoverB(hinteractive->GetBinContent(4),hinteractive->GetBinError(4),hinteractive->GetBinContent(3),hinteractive->GetBinError(3));
 cout<<"MET/HT>=0.2  cbias = "<<cbias <<" +/- "<<err<<endl;


}

void drawrplots() {
  /*
.L drawReducedTrees.C++
  */

  useFlavorHistoryWeights_=false;

  setStackMode(false);
  doData(true);

  //no met, no mindeltaphi
  setPlotMinimum(1e-2); //setPlotMaximum(1);

 drawTotalSM_=true;
 leg_y1=0.6;
 setLogY(true);

 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,20,0,200,"ge1b");

 selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,20,0,200,"eq1b");

 selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,20,0,200,"ge2b");

 //coarse
 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,4,0,200,"ge1b-4bins");

 selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,4,0,200,"eq1b-4bins");

 selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,4,0,200,"ge2b-4bins");
//  const double epsilon=-1;
//  double cbias= hinteractive->GetBinContent(4)/hinteractive->GetBinContent(3);
//  double cbias_syst = (1-cbias)/2;
//  double cbias_up = cbias+cbias_syst;
//  double cbias_down = cbias-cbias_syst;

//  double r_up = hinteractive->GetBinContent(3) * cbias_up;
//  double r_down = hinteractive->GetBinContent(3) * cbias_down;
//  TLine lup(hinteractive->GetBinCenter(4)+epsilon, r_down, hinteractive->GetBinCenter(4)+epsilon, r_up);
//  lup.SetLineColor(kMagenta);
//  lup.SetLineWidth(2);
//  lup.Draw("SAME");

//not so useful
 selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 && bestTopMass>300";
 drawR("minDeltaPhi",0.3,4,0,200,"ge2b-m3j300");

 //more useful? maybe
 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 && bestTopMass>450";
 drawR("minDeltaPhi",0.3,4,0,200,"ge1b-m3j450");

 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 && bestTopMass>450 && MET>=100&&MET<150";
 drawR("minDeltaPhi",0.3,5,100,150,"ge1b-SBonly-m3j450");


 //njets == 3
// selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && cutCleaning==1 && njets==3";
// setLogY(true);
//  drawR("minDeltaPhi",0.3,50,0,250);

//  //4 jets
// selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && cutCleaning==1 && njets==4";
// setLogY(true);
//  drawR("minDeltaPhi",0.3,50,0,250);

//  //5 jets
// selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && cutCleaning==1 && njets==5";
// setLogY(true);
//  drawR("minDeltaPhi",0.3,50,0,250);

}

void drawForPreAp_smonly() {
  /*
comment out LM13

.L drawReducedTrees.C++
  */
  setStackMode(false,false);
  doData(true);
  doRatioPlot(false);
  //  setPadDimensions(700,500);
  setLogY(true);

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  drawTotalSM_=true;
  drawMarkers_=false;

  var="MET"; xtitle="E_{T}^{miss} [GeV]";

  //log scale, no stack, full MET range
  nbins=35; low=0; high=350;
  setPlotMinimum(0.1);
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_ge1b");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_eq1b");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_ge2b");

  //now linear scale with custom y range
  resetPlotMinimum();
  setLogY(false);

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  setPlotMaximum(35);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_ge1b");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  setPlotMaximum(25);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_eq1b");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  setPlotMaximum(11);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_ge2b");

  //now fail minDeltaPhi for high-ish MET
  resetPlotMinimum();
  setLogY(false);
  nbins=35; low=0; high=350;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>0";
  setPlotMaximum(125);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_IDP_ge1b");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>0";
  setPlotMaximum(105);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_IDP_eq1b");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>0";
  setPlotMaximum(23);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_IDP_ge2b");

  setLogY(true);
  setPlotMinimum(0.1);
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>0";
  setPlotMaximum(150);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_IDP_ge1b");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>0";
  setPlotMaximum(130);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_IDP_eq1b");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>0";
  setPlotMaximum(50);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_IDP_ge2b");

}

void drawForPreAp_susy() {

  /*
uncomment LM13

.L drawReducedTrees.C++
  */
  setStackMode(false,false);
  doData(true);
  doRatioPlot(false);
  //  setPadDimensions(700,500);
  setLogY(true);

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  drawTotalSM_=true;
  drawTotalSMSusy_=true;
  drawSusyOnly_=true;
  drawMarkers_=false;

  nbins=60; low=0; high=600;
  setPlotMinimum(0.1);
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_withSusy_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_withSusy_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_withSusy_ge2b");

  setLogY(false);
  resetPlotMinimum();

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  setPlotMaximum(11);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_withSusy_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  setPlotMaximum(7);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_withSusy_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  setPlotMaximum(5);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_withSusy_ge2b");



}

void drawForANandPAS () {
  /*
.L drawReducedTrees.C++
  */
  setStackMode(true);
  doData(true);
  doRatioPlot(false);
  resetPlotMinimum();
  setPadDimensions(700,500);
  setLogY(false);

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(false);

  var="bestTopMass"; xtitle="best 3-jet mass (GeV)";
  // === plots for AN and PAS (aspect ratio adjusted)
  const int nvarbins=5;
  const float varbins[]={0.,160.,180.,260.,400.,800.};
  //SB region
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  renormalizeBins_=false;
  setPlotMaximum(25);
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_100_150_ge1btag_vb");
  renormalizeBins_=true;
  resetPlotMaximum();
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_100_150_ge1btag_vbrn");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  renormalizeBins_=false;
  setPlotMaximum(20);
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_100_150_eq1btag_vb");
  renormalizeBins_=true;
  resetPlotMaximum();
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_100_150_eq1btag_vbrn");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  renormalizeBins_=false;
  //  setPlotMaximum(8);
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_100_150_ge2btag_vb");
  renormalizeBins_=true;
  resetPlotMaximum();
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_100_150_ge2btag_vbrn");

  //-------------- signal region ---------------
  nbins=10; low=0; high=800;

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && cutMET==1";
  renormalizeBins_=false;
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_50_inf_ge1btag_vb");
  drawPlots(var,nbins,low,high,xtitle,"Events", "mcdata_bestM3j_met_50_inf_ge1btag");
  renormalizeBins_=true;
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_50_inf_ge1btag_vbrn");


  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && cutMET==1";
  renormalizeBins_=false;
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_50_inf_eq1btag_vb");
  drawPlots(var,nbins,low,high,xtitle,"Events", "mcdata_bestM3j_met_50_inf_eq1btag");
  renormalizeBins_=true;
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_50_inf_eq1btag_vbrn");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && cutMET==1";
  renormalizeBins_=false;
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_50_inf_ge2btag_vb");
  drawPlots(var,nbins,low,high,xtitle,"Events", "mcdata_bestM3j_met_50_inf_ge2btag");
  renormalizeBins_=true;
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_50_inf_ge2btag_vbrn");

  // -- SB region, compare shapes --
  setStackMode(false);
  renormalizeBins_=false;
  nbins=10; low=0; high=800;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Arbitrary Units", "m3j_shapesInMC_SB_ge1btag");
  drawPlots(var,nvarbins,varbins,xtitle,"Arbitrary Units", "m3j_shapesInMC_SB_ge1btag_vb");


  resetPadDimensions();

}


void drawSoftJets() {
/*
.L drawReducedTrees.C++
*/

//this seems completely uninteresting!

  setStackMode(true);
  doData(true);
  int nbins;
  float low,high;
  TString var,xtitle;
  doRatioPlot(false);
  setLogY(false);

  nbins=10; low=0; high=10;
  var="nLooseJets20_30"; xtitle="Soft Jet multiplicity";
  //ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 &&cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "nsoftjets_ge1b");
  

}

void drawSomething() {
  /* for reference, here is all cuts (>=0 b)
selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && cutCleaning==1";

//here is all cuts but with ECAL cleaning removed
selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";

.L drawReducedTrees.C++

  */

  //here we imitate exactly the plots we made in plotWithScale() of the old code

  setStackMode(true);
  doData(true);
  doRatioPlot(true);

  int nbins;
  float low,high;
  TString var,xtitle;

  //  setQCDErrorMode(true);
  // PROBLEMS -- CMS style ruins the hash marks. central value on the TGraphErrors is wrong!
  //might be easiest to accumulate a total MC uncertainty and plot that!

  //the strength and weakness of this tree-based method is that I need to apply the cuts now!
  //most straightforward way is to manually set the cut string before each plot

  // TO DO add weight<1000 cut everywhere!


  // ==== MET plots ====
  setLogY(true);   setPlotMinimum(1e-1);
  nbins=25; low= 0; high=260;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;

  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET");
  //  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_JER0"); //in case we're doing the non-JERbias version

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_ge2b");

  setLogY(true);   setPlotMinimum(1e-1);
  nbins=50; low= 0; high=500;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0.5; ratioMax = 1.5;
  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_METw");

  // == MET plot in limited range (SB and up) ==
  setLogY(false);   setPlotMinimum(0);
  doRatioPlot(false);
  nbins=10; low= 100; high=200;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_SB_ge1b");

  // == MET plot in fail minDeltaPhi
  nbins=10; low= 100; high=200;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_SBfail_ge1b");

  doRatioPlot(true);
  // === MET plot in ILV/T2 region ===
  setLogY(false);   setPlotMinimum(0);
  nbins=30; low= 0; high=300;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1) || (nElectrons==1 && nMuons==0)) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_ILV_ge1b");

  nbins=20; low= 0; high=300;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_IEV_ge1b");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_IMV_ge1b");

  // == replacements for Don's ILV plots (but without fancy ttbar categories)
  //these are in the MET>150 region
  doRatioPlot(false);
  setLogY(false);   setPlotMinimum(0);
  nbins=6; low=1; high=7;
  var="njets"; xtitle="Jet multiplicity";
  //ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_imv_1btag");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_imv_eq1btag");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_imv_2btag");

  doRatioPlot(true);
  setLogY(false);   setPlotMinimum(0);
  nbins=10; low= 0; high=200;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_imv_1btag");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_imv_eq1btag");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_imv_2btag");

  doRatioPlot(false);
  setLogY(false);   setPlotMinimum(0);
  nbins=6; low=1; high=7;
  var="njets"; xtitle="Jet multiplicity";
  //ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_iev_1btag");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_iev_eq1btag");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_iev_2btag");

  doRatioPlot(true);
  setLogY(false);   setPlotMinimum(0);
  nbins=10; low= 0; high=200;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_iev_1btag");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_iev_eq1btag");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_iev_2btag");


  // ==== MHT plot === (no MET cut)
  setLogY(true);  setPlotMinimum(1e-1);
  nbins=50; low= 0; high=500;
  ratioMin = 0.5; ratioMax = 1.5;
  var="MHT"; xtitle="H_{T}^{miss} [GeV]";

  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MHT");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MHT_ge1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MHT_ge2b");

  //to compare to RA2
  setLogY(true);  setPlotMinimum(1);
  nbins=16; low= 0; high=160;
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MHTveryloose");


  // ==== minDeltaPhi plots ====
  setLogY(false);
  resetPlotMinimum();
  nbins=10; low=0; high = TMath::Pi() + 0.001;
  var="minDeltaPhi"; xtitle="min(#Delta#phi[ jets 1..3, E_{T}^{miss} ] )";

  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj");
  //drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_JER0");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_ge2b");

  // ==== HT plots ==== //in this case this is not an N-1 plot, because we've already cut on HT
  ratioMin = 0; ratioMax = 3;
  resetPlotMinimum();
  setLogY(false);
  nbins=10; low=300; high = 1000;
  var="HT"; xtitle="H_{T} (GeV)";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_HT_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_HT_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_HT_ge2b");

  // ==== n jets ====
  nbins=8; low=0; high = 8;
  resetPlotMinimum();
  setLogY(false);
  var="njets"; xtitle="Jet multiplicity";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hnjets_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hnjets_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 &&  MET>=100 && MET<150&& cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hnjets_ge2b");

  // ==== n b jets ====
  ratioMin = 0; ratioMax = 2;
  nbins=4; low=0; high = 4;
  resetPlotMinimum();
  setLogY(false);
  var="nbjets"; xtitle="Number of b-tagged jets";
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && cut3Jets==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hnbjets");

  // ==== n Electrons ====
  nbins=3; low=0; high= 3;
  resetPlotMinimum();
  setLogY(false);
  var="nElectrons"; xtitle="Number of electrons";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnElectrons_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnElectrons_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnElectrons_ge2b");

  // ==== n Muons ====
  nbins=3; low=0; high= 3;
  resetPlotMinimum();
  setLogY(false);
  var="nMuons"; xtitle="Number of muons";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnMuons_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnMuons_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnMuons_ge2b");

  // ==== lead jet pT ====
  nbins=10; low=50; high= 450;
  resetPlotMinimum();
  setLogY(false);
  var="jetpt1"; xtitle="p_{T} of lead jet (GeV)";

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_ge2b");

  // ======= Owen's variables =====
  ratioMin = 0; ratioMax = 3;
  nbins=10; low=50; high= 550;
  resetPlotMinimum();
  setLogY(false);
  var="bestTopMass"; xtitle="3-jet mass (GeV)";

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge2b");

  //signal region
  nbins=20; low=0; high= 800;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge1b_SIG");
  //SB region
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge1b_wide");
  //LSB region
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET<50";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge1b_LSB");

  // ======= study cleaning  =======
  //  bool passBadPFMuon, passInconsistentMuon, passEcalCleaning;
  //no MET cut -- look at BadPFMuon events
  setLogY(true);  resetPlotMinimum();
  nbins=20; low= 0; high=260;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;

  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passBadPFMuon==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "MET_BadPFMuon");

  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "MET_InconsistentMuon");

  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passEcalCleaning==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "MET_FailBE");

  //plot HT (signal region), no b tagging, with and without cleaning
  resetPlotMinimum();
  setLogY(false);
  nbins=10; low=300; high = 1000;
  var="HT"; xtitle="H_{T} (GeV)";
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passEcalCleaning==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HT_FailBE");


   // ==== lead b jet pT ====
   nbins=10; low=30; high= 200;
   resetPlotMinimum();
   setLogY(false);
   var="bjetpt1"; xtitle="p_{T} of lead b jet (GeV)";

   selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
   drawPlots(var,nbins,low,high,xtitle,"Events", "Hbjetpt1_ge1b");

//   selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
//   drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_eq1b");

//   selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
//   drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_ge2b");
   
// ==== number of primary vertices ====
   nbins =9; low=1; high=9;
   resetPlotMinimum();
   setLogY(false);
   var="nGoodPV"; xtitle="# of good PVs";
   selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "nGoodPV_SB_ge1b");
   selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "nGoodPV_SBfail_ge1b");

  //check bad PV events
   nbins =6; low=0; high=5;
   resetPlotMinimum();
   setLogY(false);
   var="nGoodPV"; xtitle="# of good PVs";
   selection_ ="nbjets>=1 && cutHT==1 && cutPV==0 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&MET>50";
   drawPlots(var,nbins,low,high,xtitle,"Events", "nGoodPV_PVfail_ge1b");
   selection_ ="cutHT==1 && cutPV==0 && cutTrigger==1 && cutDeltaPhi==1";
   drawPlots(var,nbins,low,high,xtitle,"Events", "nGoodPV_PVfail_veryloose");


   //exploration
   setLogY(false);
   resetPlotMinimum();
 
   nbins=20; low=-TMath::Pi(); high = TMath::Pi();
   var="METphi-jetphi1"; xtitle="(METphi-jetphi1)";
   selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>100 && passInconsistentMuon==1 && passBadPFMuon==1";
   drawPlots(var,nbins,low,high,xtitle,"Events", "HdeltaPhiMETjet1_nomindp");

   nbins=20; low=-TMath::Pi(); high = TMath::Pi();
   var="METphi-jetphi2"; xtitle="METphi-jetphi2";
   selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>100 && passInconsistentMuon==1 && passBadPFMuon==1 && (abs(METphi-jetphi1)>1)";
   drawPlots(var,nbins,low,high,xtitle,"Events", "HdeltaPhiMETjet2_metphi1cut");

   var="METphi-jetphi2"; xtitle="METphi-jetphi2";
   selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>150 && passInconsistentMuon==1 && passBadPFMuon==1 && (abs(METphi-jetphi1)>1)";
   drawPlots(var,nbins,low,high,xtitle,"Events", "HdeltaPhiMETjet2_metphi1cut_sig");

   var="METphi-jetphi3"; xtitle="METphi-jetphi3";
   selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>150 && passInconsistentMuon==1 && passBadPFMuon==1 && (abs(METphi-jetphi1)>1) && (abs(METphi-jetphi2)>0.3)";
   drawPlots(var,nbins,low,high,xtitle,"Events", "HdeltaPhiMETjet3_metphi12cut_sig");

  // ==== MET/HT ====
   //addSample("LM13");
  doRatioPlot(false);
  setLogY(true);
  nbins=40;
  low=0;
  high=1;
  doOverflowAddition(true);
  var="MET/HT"; xtitle="MET/HT";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverHt_ge1b_noMetCut");

  setLogY(false);
  resetPlotMinimum();
  nbins=10;

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 && MET>100 &&MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverHt_ge1b_SB");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&MET>=150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverHt_ge1b_SIG");

  // ==== MET/sqrt(HT) ====
   //addSample("LM13");
  doRatioPlot(false);
  nbins=20;  low=0;  high=20;
  doOverflowAddition(true);
  var="MET/sqrt(HT)"; xtitle="MET/#sqrt{HT}";

   setLogY(false);
  resetPlotMinimum();

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 && MET>100 &&MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverRootHt_ge1b_SB");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&MET>=150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverRootHt_ge1b_SIG");


  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&MET>=100";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverRootHt_ge1b_metOver100");

}

void studySignificance() {
/*
.L drawReducedTrees.C++
  */

  setStackMode(true);
  doData(true);

  int nbins;
  float low,high;
  TString var,xtitle;
  doOverflowAddition(true);
  doRatioPlot(false);


  setLogY(false);
  resetPlotMinimum();

  //this compares the effectiveness of MET and MET/sqrt(HT); MET is better
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";

  nbins=20;  low=0;  high=20;
  var="MET/sqrt(HT)"; xtitle="MET/#sqrt{HT}";
  drawSignificance(var,nbins,low,high,"sigfile");

  nbins=30;  low=0;  high=300;
  var="MET"; xtitle="MET";
  drawSignificance(var,nbins,low,high,"sigfile");


  //after MET cut, what can we do?
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=20;  low=0;  high=20;
  var="MET/sqrt(HT)"; xtitle="MET/#sqrt{HT}";
  drawSignificance(var,nbins,low,high,"sigfile"); //nothing good here

  //reoptimize minDeltaPhi
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=31;  low=0;  high=TMath::Pi();
  var="minDeltaPhi"; xtitle="minDeltaPhi";
  drawSignificance(var,nbins,low,high,"sigfile");

  var="minDeltaPhiAll30"; xtitle="minDeltaPhi All 30"; //not as good!
  drawSignificance(var,nbins,low,high,"sigfile");


  //remove 3Jet cut
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=7; low=0; high=7;
  var="njets"; xtitle="jet multiplicity";
  drawSignificance(var,nbins,low,high,"sigfile");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=40; low=200; high=1000;
  var="HT"; xtitle="HT";
  drawSignificance(var,nbins,low,high,"sigfile");

  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=4; low=0; high=4;
  var="nbjets"; xtitle="n b tags";
  drawSignificance(var,nbins,low,high,"sigfile");

  //nominal cuts
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=40;  low=0;  high=400;
  var="jetpt1"; xtitle="lead jet pt";
  drawSignificance(var,nbins,low,high,"sigfile");

  nbins=40;  low=0;  high=400;
  var="bjetpt1"; xtitle="lead b jet pt";
  drawSignificance(var,nbins,low,high,"sigfile");

}

void drawMinDeltaPhiMETslices(){
  useFlavorHistoryWeights_=false;
  doOverflowAddition(true);
  doData(true);
  doRatioPlot(true);
  setQuiet(false);
  ratioMin=0; ratioMax=3;
  
  bool doVB = false;
  bool doSlices = true;
  if(!doSlices){
    leg_y1=0.6; //and also remember to use correct color scheme
    owenColor_=true;
  }  

  int nbins;
  float low,high;
  TString var,xtitle;

  nbins=10; low=0; high = TMath::Pi() + 0.001;
  var="minDeltaPhi"; xtitle="min(#Delta#phi[ jets 1..3, E_{T}^{miss} ] )";
  const int nvarbins=2;
  const float varbins[]={0.,0.3,high};
  
  TCut baseSelection =""; 
  TCut theseCuts = ""; 
  TString extraName = "";

  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjets == 0";
  for (int ibtag = 0; ibtag<5; ibtag++) { 
    TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
    if (ibtag==0) { //nothing to do
    }
    else if (ibtag==1) {
      theBTaggingCut = eq1b; 
      btagstring = "eq1b";
    }
    else if (ibtag==2) {
      theBTaggingCut = ge2b; 
      btagstring = "ge2b";
    }
    else if (ibtag==3) {
      theBTaggingCut = pretag;
      btagstring = "pre";
    }
    else if (ibtag==4) {
      theBTaggingCut = antitag;
      btagstring = "antib";
    }
    else assert(0);
    

    //+++++++++++NOMINAL SELECTION+++++++++++//
    baseSelection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1";//nominal
    extraName = "_nominal";
     
    if(doSlices){
      //MinDeltaPhi in MET slices
      setStackMode(true); setLogY(false); resetPlotMinimum();
      //MET < 50
      theseCuts = "MET<50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50_vb");
      //50 < MET < 100 
      theseCuts = "MET<100. && MET>=50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100_vb");
      //100 < MET < 150
      theseCuts = "MET<150. && MET>=100.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150_vb");
      //MET < 150
      theseCuts = "MET>=150.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf"); 
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf_vb");
    }
    else{
      //Ratio plot
      selection_ = baseSelection && theBTaggingCut;
      setStackMode(false); setPlotMinimum(1e-2); drawTotalSM_=true; setLogY(true); 
      drawR("minDeltaPhi",0.3,20,0,200,btagstring+extraName);
      drawR("minDeltaPhi",0.3,4,0,200,btagstring+extraName+"_4bins");
    }
    
    /*
    //[+++++++++++++njets==2++++++++]//
    baseSelection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1 && njets==2"; //remove cut3Jets and add njets
    extraName = "_njetsEQ2";

    if(doSlices){
      //MinDeltaPhi in MET slices
      setStackMode(true); setLogY(false); resetPlotMinimum();
      //MET < 50
      theseCuts = "MET<50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50_vb");
      //50 < MET < 100 
      theseCuts = "MET<100. && MET>=50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100_vb");
      //100 < MET < 150
      theseCuts = "MET<150. && MET>=100.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150_vb");
      //MET < 150
      theseCuts = "MET>=150.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf_vb");
    }
    else{
      //Ratio plot
      selection_ = baseSelection && theBTaggingCut;
      setStackMode(false); setPlotMinimum(1e-2); drawTotalSM_=true; setLogY(true); 
      drawR("minDeltaPhi",0.3,20,0,200,btagstring+extraName);
      drawR("minDeltaPhi",0.3,4,0,200,btagstring+extraName+"_4bins");
    }
    */
    
    /*
    //[+++++++++++m3j>350++++++++++++]//
    baseSelection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1 && bestTopMass>350";//nominal plus m3j cut
    extraName = "_m3jGT350";
    
    if(doSlices){
      //MinDeltaPhi in MET slices
      setStackMode(true); setLogY(false); resetPlotMinimum();
      //MET < 50
      theseCuts = "MET<50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50_vb");
      //50 < MET < 100 
      theseCuts = "MET<100. && MET>=50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100_vb");
      //100 < MET < 150
      theseCuts = "MET<150. && MET>=100.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150_vb");
      //MET < 150
      theseCuts = "MET>=150.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf_vb");
    }
    else{
      //Ratio plot
      selection_ = baseSelection && theBTaggingCut;
      setStackMode(false); setPlotMinimum(1e-2); drawTotalSM_=true; setLogY(true); 
      drawR("minDeltaPhi",0.3,20,0,200,btagstring+extraName);
      drawR("minDeltaPhi",0.3,4,0,200,btagstring+extraName+"_4bins");
    }
    */

  }//b-tag cut loop 
}


void drawMETPlots() {
  useFlavorHistoryWeights_=false;
  setStackMode(true);
  doData(true);
  doRatioPlot(true);
  setQuiet(false);
  doOverflowAddition(true);
  loadSamples();

  int nbins;
  float low,high;
  TString var,xtitle;

  // ==== MET plots ====
  setLogY(true);   setPlotMinimum(1e-1);
  nbins=10; low= 0.; high=150.;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 3;

  /*
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_ge2b");

  selection_ ="nbjets==0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_antib");
  */ 


  nbins=20; low=0; high=400;
  //double tempLumiScale =  lumiScale_;
  //double thisLumi = 94.432+166.732+216.430+154.830+84.698+57.994+2.101;
  //cout << "thisLumi: " << thisLumi << endl;
  //lumiScale_ = thisLumi/190.5;
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_test");
  //lumiScale_= tempLumiScale;
}


void drawMETPlots_utility(){
  useFlavorHistoryWeights_=false;
  setStackMode(true);
  doData(true);
  doRatioPlot(true);
  setQuiet(false);
  doOverflowAddition(true);
  loadSamples();
  ratioMin = 0; ratioMax = 3;

  int nbins;
  float low,high;
  TString var,xtitle;
  bool doVB=false;

  TCut baseSelection = "cutHT==1 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1"; //cuts not present: trigger, MET, mdp, bjets
  TCut utilSelection = "pass_utilityHLT_HT300>=1";
  TCut nominalTrigger = "cutTrigger==1";
  TCut mdpCut = "cutDeltaPhi==1";
  TCut invertMdpCut = "cutDeltaPhi==0";

  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjets == 0";
  for (int ibtag = 0; ibtag<5; ibtag++) { 
    TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
    if (ibtag==0) { //nothing to do
    }
    else if (ibtag==1) {
      theBTaggingCut = eq1b; 
      btagstring = "eq1b";
    }
    else if (ibtag==2) {
      theBTaggingCut = ge2b; 
      btagstring = "ge2b";
    }
    else if (ibtag==3) {
      theBTaggingCut = pretag;
      btagstring = "pre";
    }
    else if (ibtag==4) {
      theBTaggingCut = antitag;
      btagstring = "antib";
    }
    else assert(0);
    
    TString extraName = "_nominal";
    
    // ==== MET plots ====
    setLogY(true);   setPlotMinimum(1e-1);
    nbins=10; low= 0.; high=150.;
    var="MET"; xtitle="E_{T}^{miss} [GeV]";
    selection_ = baseSelection && utilSelection && mdpCut && theBTaggingCut; //no MET cut
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_METstack");
    setLogY(false); resetPlotMinimum();
    //high MET
    nbins=1; low=150.; high=200.;
    selection_ = baseSelection && utilSelection && invertMdpCut && theBTaggingCut; //no MET cut, inverted mdp
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_imdp_METstack");
    //cross check hack
    double tempLumiScale =  lumiScale_;
    double thisLumi = 94.432+166.732+216.430+154.830+84.698+57.994+2.101;
    cout << "thisLumi: " << thisLumi << endl;
    lumiScale_ = thisLumi/190.5;
    selection_ = baseSelection && nominalTrigger && invertMdpCut && theBTaggingCut; //no MET cut, inverted mdp
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_imdp_METstack_full");
    lumiScale_= tempLumiScale;
    

    if(0){
      // ==== MinDeltaPhi in MET slices ====
      setLogY(false); resetPlotMinimum();
      nbins=7;  low = 0; high = TMath::Pi() + 0.001;
      var="minDeltaPhi"; xtitle="min(#Delta#phi[ jets 1..3, E_{T}^{miss} ] )";
      const int nvarbins=2;
      const float varbins[]={0.,0.3,high};
      TCut theseCuts;
      //MET < 50
      theseCuts = "MET<50.";
      selection_ = baseSelection && utilSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50_vb");
      //50 < MET < 100 
      theseCuts = "MET<100. && MET>=50.";
      selection_ = baseSelection && utilSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100_vb");
      //100 < MET < 150
      theseCuts = "MET<150. && MET>=100.";
      selection_ = baseSelection && utilSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150_vb");
      //MET < 150
      theseCuts = "MET>=150.";
      selection_ = baseSelection && utilSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf"); 
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf_vb");
    }
  }//end loop over b-tag selection
}


void drawQCDreweight(bool loose = true){
  loadSamples(true);
  savePlots_=false; //don't save eps,png,pdf files
  doOverflowAddition(true);

  //output file 
  TString histfilename = "";
  histfilename += "qcdReweight";
  if(loose) histfilename += "_loose";
  else histfilename += "_tight";
  histfilename += ".root";
  TFile fh(histfilename,"RECREATE");//will delete old root file if present 
  fh.Close(); //going to open only when needed 
  
  //-- define some cuts to use 
  //--------------------------
  const TCut phys = "cutTrigger==1";
  const TCut util = "pass_utilityHLT_HT300>=1";
  const TCut LSB = "MET>=50 && MET<100";
  const TCut SB = "MET>=150 && MET<200";
  TCut SIG;
  if(loose) SIG = "MET>=200";
  else SIG =  "MET>=300";
  TCut HTcut;
  if(loose) HTcut = "HT>=350";
  else HTcut = "HT>=500";
  const  TCut ge1b =  "nbjetsSSVHPT >= 1";
  const  TCut ge2b =  "nbjetsSSVHPT >= 2";
  const  TCut eq1b =  "nbjetsSSVHPT == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjetsSSVHPT == 0";
  const TCut passmdpn = "minDeltaPhiN>=4";
  const TCut failmdpn = "minDeltaPhiN<4";
  TCut base = HTcut&&TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && weight<1000");// no trigger, met, minDeltaPhiN, btagger
  
  //binning
  int nbins;
  double low;
  double high;
  
  //-- LSB plots
  //------------
  lumiScale_ = 19.23288;

  nbins=4;
  low=2.5;
  high=6.5;
  selection_ = base && LSB && util && failmdpn && antitag;
  drawSimple("njets",nbins,low,high,histfilename, "h_LSB_fDP_antib_njets_data", "data");
  drawSimple("njets",nbins,low,high,histfilename, "h_LSB_fDP_antib_njets_qcd", "PythiaPUQCD");
  drawPlots( "njets",nbins,low,high, "", "",      "deleteme");   //drawPlots( "njets", 4, 2.5, 6.5, "", "", "test");
  TFile fh1(histfilename, "UPDATE");
  totalnonqcd->SetName("h_LSB_fDP_antib_njets_nonqcd");
  totalnonqcd->Write();
  fh1.Close();

  selection_ = base && LSB && util && passmdpn && antitag;
  drawSimple("njets",nbins,low,high,histfilename, "h_LSB_pDP_antib_njets_qcd", "PythiaPUQCD");
  
  //-- SB plots
  //-----------
  lumiScale_ = 1091.891;
  
  nbins=4;
  low=2.5;
  high=6.5;
  selection_ =  base && SB && phys && failmdpn && ge1b;
  drawSimple("njets",nbins,low,high,histfilename, "h_SB_fDP_ge1b_njets_data", "data");
  drawSimple("njets",nbins,low,high,histfilename, "h_SB_fDP_ge1b_njets_qcd", "PythiaPUQCD");
  drawPlots( "njets",nbins,low,high, "", "",      "deleteme");   
  TFile fh2(histfilename, "UPDATE");
  totalnonqcd->SetName("h_SB_fDP_ge1b_njets_nonqcd");
  totalnonqcd->Write();
  fh2.Close();
  
  selection_ =  base && SB && phys && failmdpn && ge2b;
  drawSimple("njets",nbins,low,high,histfilename, "h_SB_fDP_ge2b_njets_data", "data");
  drawSimple("njets",nbins,low,high,histfilename, "h_SB_fDP_ge2b_njets_qcd", "PythiaPUQCD");
  drawPlots( "njets",nbins,low,high, "", "",      "deleteme");
  TFile fh3(histfilename, "UPDATE");
  totalnonqcd->SetName("h_SB_fDP_ge2b_njets_nonqcd");
  totalnonqcd->Write();
  fh3.Close();

  selection_ = base && SB && phys && passmdpn && ge1b;
  drawSimple("njets",nbins,low,high,histfilename, "h_SB_pDP_ge1b_njets_qcd", "PythiaPUQCD");

  selection_ = base && SB && phys && passmdpn && ge2b;
  drawSimple("njets",nbins,low,high,histfilename, "h_SB_pDP_ge2b_njets_qcd", "PythiaPUQCD");

  //-- SIG plots
  //------------
  if(loose){
    nbins=4;
    low=2.5;
    high=6.5;
  }
  else{
    nbins=3;
    low=2.5;
    high=5.5;
  }

  selection_ =  base && SIG && phys && failmdpn && ge1b;
  drawSimple("njets",nbins,low,high,histfilename, "h_SIG_fDP_ge1b_njets_data", "data");
  drawSimple("njets",nbins,low,high,histfilename, "h_SIG_fDP_ge1b_njets_qcd", "PythiaPUQCD");
  drawPlots( "njets",nbins,low,high, "", "",      "deleteme");
  TFile fh4(histfilename, "UPDATE");
  totalnonqcd->SetName("h_SIG_fDP_ge1b_njets_nonqcd");
  totalnonqcd->Write();
  fh4.Close();

  if(loose){
    nbins=3;
    low=2.5;
    high=5.5;
  }
  else{
    nbins=2;
    low=2.5;
    high=4.5;
  }
  selection_ =  base && SIG && phys && failmdpn && ge2b;
  drawSimple("njets",nbins,low,high,histfilename, "h_SIG_fDP_ge2b_njets_data", "data");
  drawSimple("njets",nbins,low,high,histfilename, "h_SIG_fDP_ge2b_njets_qcd", "PythiaPUQCD");
  drawPlots( "njets",nbins,low,high, "", "",      "deleteme");
  TFile fh5(histfilename, "UPDATE");
  totalnonqcd->SetName("h_SIG_fDP_ge2b_njets_nonqcd");
  totalnonqcd->Write();
  fh5.Close();

  nbins=4;
  low=2.5;
  high=6.5;
  selection_ = base && SIG && phys && passmdpn && ge1b;
  drawSimple("njets",nbins,low,high,histfilename, "h_SIG_pDP_ge1b_njets_qcd", "PythiaPUQCD");
  
  nbins=3;
  low=2.5;
  high=5.5;
  selection_ = base && SIG && phys && passmdpn && ge2b;
  drawSimple("njets",nbins,low,high,histfilename, "h_SIG_pDP_ge2b_njets_qcd", "PythiaPUQCD");

}


void studyPrescale_r(int ibtag = 4) {
  loadSamples();
  drawTotalSM_=true;
  drawLegend(true);
  setColorScheme("nostack");
  setStackMode(false);
  doOverflowAddition(true);
  doRatioPlot(true);
  doData(true);

  lumiScale_ = 19.23288; //customize lumiScale_

  //TString btagselection="antib";
  //TCut btagcut = "nbjetsSSVHPT==0";
  
  TCut util = "pass_utilityHLT_HT300>=1 && weight<1000";
  TCut LSB = "MET>=50 && MET<100";

  const  TCut ge1b =  "nbjetsSSVHPT >= 1";
  const  TCut ge2b =  "nbjetsSSVHPT >= 2";
  const  TCut eq1b =  "nbjetsSSVHPT == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjetsSSVHPT == 0";
 
  //  for (int ibtag = 4; ibtag<5; ibtag++) { 
  TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
  if (ibtag==0) { //nothing to do
  }
  else if (ibtag==1) {
    theBTaggingCut = eq1b; 
    btagstring = "eq1b";
  }
  else if (ibtag==2) {
    theBTaggingCut = ge2b; 
    btagstring = "ge2b";
  }
  else if (ibtag==3) {
    theBTaggingCut = pretag;
    btagstring = "pre";
  }
  else if (ibtag==4) {
    theBTaggingCut = antitag;
    btagstring = "antib";
  }
  else assert(0);
  
  doOverflowAddition(false);
  //cross check with note
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN",4,10,50,150,"Data_eq0b");
  //drawR("minDeltaPhiN",4,1,50,100,"Data_eq0b");
  
  //trigger version
  //set dataOnly and might want to change hdata->SetXTitle()
  doOverflowAddition(false);
  //setPlotMaximum(0.5); setPlotMinimum(0);
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "pass_utilityHLT_HT300", 7, 1.5, 8.5, "triggerVersion_"+btagstring);
  //unexpected value for MC in last bin when adding overflow.
  
  //nGoodPV
  doOverflowAddition(true);
  //setPlotMaximum(0.5); setPlotMinimum(0);
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "nGoodPV", 20, 0.5, 20.5, "nGoodPV_"+btagstring);
  //const int nvarbins=9;
  //const float varbins[]={0.5,1.5,2.5,3.5,4.5,5.5,6.5,8.5,11.5,16.5};
  //drawR("minDeltaPhiN", 4, "nGoodPV", nvarbins, varbins, "nGoodPV_"+btagstring);
  //dependence is seen here
  //
  //prescale vs physics check
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  drawSimple("nGoodPV",20,0,20,"nGoodPV.root", "nGoodPV_"+btagstring+"tag_prescale_data","data");
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=150")&&theBTaggingCut;
  drawSimple("nGoodPV",20,0,20,"nGoodPV.root", "nGoodPV_"+btagstring+"tag_physics_data","data");
  


  //njets
  doOverflowAddition(true);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "njets", 8, 2.5, 10.5, "njets_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=150 && weight<1000")&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "njets",8, 2.5, 10.5, "njets_physics_"+btagstring);
  
  //nbjets -- does not depend on theBTaggingCut!!!
  doOverflowAddition(false);
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB;
  //drawR("minDeltaPhiN", 4, "nbjetsSSVHPT", 4, -0.5, 3.5, "nbjets_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=150");
  //drawR("minDeltaPhiN", 4, "nbjetsSSVHPT", 4, -0.5, 3.5, "nbjets_physics_"+btagstring);
   
  //runnumber
  //set dataOnly to true
  //should hack after first instance of thecanvas->cd(1); to add:  gPad->SetRightMargin(.1); gPad->Modified();
  doOverflowAddition(false);
  //setPlotMaximum(0.5); setPlotMinimum(0);
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "runNumber", 10, 160403.5, 167913.5, "runNumber_"+btagstring);
  //const int nvarbins2 = 7;
  //const float varbins2[] = {160403.5, 161000, 162000, 163500, 165000,166000,167000,168000};
  //drawR("minDeltaPhiN", 4, "runNumber", nvarbins2, varbins2, "runNumber_"+btagstring);
  
  
  //lowMET
  doOverflowAddition(false);
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&theBTaggingCut;
  const int nvarbins=10;
  const float varbins[]={0,10,20,30,40,50,60,70,80,90,100};
  //  drawR("minDeltaPhiN", 4, "MET", nvarbins, varbins, "MET_"+btagstring);
  //}//end btag loop
}
 
void lowMETcount(){
  loadSamples();
  doOverflowAddition(false);
  setQuiet(false);
  setStackMode(true); //regular stack
  setColorScheme("stack");
  doData(true);
  drawMCErrors_=true;
  
  lumiScale_ = 19.23288; //customize lumiScale_
  
  TString btagselection="antib";
  TCut theBTaggingCut = "nbjetsSSVHPT==0";
  TCut util = "pass_utilityHLT_HT300>=1 && weight<1000";
  
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&theBTaggingCut;
  const int nvarbins=3;
  const float varbins[]={40,50,100,110};
  drawPlots("MET", nvarbins, varbins, "","","deleteme");
  cout << "BEN: Data integral1:::::::: " << hdata->GetBinContent(1) << endl;
  cout << "BEN: Data integral2:::::::: " << hdata->GetBinContent(2) << endl;
  cout << "BEN: Data integral3:::::::: " << hdata->GetBinContent(3) << endl;
}//end of function


void sigbpt(){
  loadSamples();
  doOverflowAddition(true);
  
  TString btagstring="ge1b";
  TCut theBTaggingCut = "nbjets>=1";

  //loose
  selection_ =TCut("(HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN>=4 && MET>=200)")&&theBTaggingCut;
  cout << "LM9: " << drawSimple("bjetpt1",20,30,500,"bjetpt1.root", "bjetpt1_"+btagstring+"tagLooseSignalRegion_lm9", "LM9") << endl;;
}

void studyNjet(int ibtag=3){
  loadSamples();
  drawTotalSM_=true;
  drawLegend(true);
  //setColorScheme("nostack");
  //setStackMode(false);
  doOverflowAddition(true);
  doRatioPlot(true);
  doData(true);

  //lumiScale_ = 19.23288; //customize lumiScale_
  //TString btagselection="antib"; 
  //TCut btagcut = "nbjetsSSVHPT==0"; 

  TCut util = "pass_utilityHLT_HT300>=1 && weight<1000";
  TCut LSB = "MET>=50 && MET<100";

  const  TCut ge1b =  "nbjetsSSVHPT >= 1";
  const  TCut ge2b =  "nbjetsSSVHPT >= 2";
  const  TCut eq1b =  "nbjetsSSVHPT == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjetsSSVHPT == 0";

  //  for (int ibtag = 4; ibtag<5; ibtag++) {
  TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
  if (ibtag==0) { //nothing to do 
  }
  else if (ibtag==1) {
    theBTaggingCut = eq1b;
    btagstring = "eq1b";
  }
  else if (ibtag==2) {
    theBTaggingCut = ge2b;
    btagstring = "ge2b";
  }
  else if (ibtag==3) {
    theBTaggingCut = pretag;
    btagstring = "pre";
  }
  else if (ibtag==4) {
    theBTaggingCut = antitag;
    btagstring = "antib";
  }
  else assert(0);


  doOverflowAddition(true);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "njets", 7, 2.5, 10.5, "njets_"+btagstring); 
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=150")&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "njets", 7, 2.5, 10.5, "njets_physics_"+btagstring); 

  //njet distributions -- note - NO minDeltahiN cut!
  //selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=150 && weight<1000")&&theBTaggingCut;
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=150 && weight<1000")&&ge2b;
  //drawSimple("njets",10,0,10,"njets.root", "njets_"+btagstring+"tag_sig_qcd","PythiaPUQCD");
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=50 && MET<100 && weight<1000")&&theBTaggingCut;
  //drawSimple("njets",10,0,10,"njets.root", "njets_"+btagstring+"tag_lsb_qcd","PythiaPUQCD");


  //N JET DATA/MC
  setStackMode(true); setLogY(true); resetPlotMinimum();    
  TCut failmdp = "minDeltaPhiN<4.";

  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=150 && weight<1000")&&theBTaggingCut && failmdp;
  drawPlots("njets",7,2.5,10,"njets","Events", "H_njets_physics_"+btagstring);
  
  lumiScale_ = 19.23288; //customize lumiScale_
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut&&failmdp;
  drawPlots("njets",7,2.5,10,"njets","Events", "H_njets_prescaled_"+btagstring);
}

void studyRqcd(int ibtag = 4){
  loadSamples();
  drawTotalSM_=false;
  drawLegend(true);
  setColorScheme("nostack");
  setStackMode(false);
  doOverflowAddition(true);
  doRatioPlot(true);
  doData(false);

  lumiScale_ = 19.23288; //customize lumiScale_

  //TString btagselection="antib";
  //TCut btagcut = "nbjetsSSVHPT==0";
  
  TCut util = "pass_utilityHLT_HT300>=1 && weight<1000";
  TCut LSB = "MET>=50 && MET<100";
  TCut SB = "MET>=150 && MET<200";
  TCut SIG = "MET>=200";
  TCut highMET = "MET>=150";

  const  TCut ge1b =  "nbjetsSSVHPT >= 1";
  const  TCut ge2b =  "nbjetsSSVHPT >= 2";
  const  TCut eq1b =  "nbjetsSSVHPT == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjetsSSVHPT == 0";
 
  //  for (int ibtag = 4; ibtag<5; ibtag++) { 
  TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
  if (ibtag==0) { //nothing to do
  }
  else if (ibtag==1) {
    theBTaggingCut = eq1b; 
    btagstring = "eq1b";
  }
  else if (ibtag==2) {
    theBTaggingCut = ge2b; 
    btagstring = "ge2b";
  }
  else if (ibtag==3) {
    theBTaggingCut = pretag;
    btagstring = "pre";
  }
  else if (ibtag==4) {
    theBTaggingCut = antitag;
    btagstring = "antib";
  }
  else assert(0);
  
  //clearSamples();
  //addSample("PythiaPUQCD");
  
  /*
    // NJETS
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  drawR("minDeltaPhiN", 4, "njets", 8, 2.5, 10.5, "njets_lsb_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&SB&&theBTaggingCut;
  drawR("minDeltaPhiN", 4, "njets", 8, 2.5, 10.5, "njets_sb_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&SIG&&theBTaggingCut;
  drawR("minDeltaPhiN", 4, "njets", 8, 2.5, 10.5, "njets_sig_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&highMET&&theBTaggingCut;
  drawR("minDeltaPhiN", 4, "njets", 8, 2.5, 10.5, "njets_highMET_"+btagstring);
  */

  //NJETS FOR JOSH
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  drawR("minDeltaPhiN", 4, "njets", 4, 2.5, 6.5, "njets_joshBin_lsb_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&highMET&&theBTaggingCut;
  drawR("minDeltaPhiN", 4, "njets", 4, 2.5, 6.5, "njets_joshBin_highMET_"+btagstring);
  
  setPlotMaximum(0.35); setPlotMinimum(0);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cut3Jets==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "MET", 1, 50, 100, "test_"+btagstring);

  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cut3Jets==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "MET", 15, 50, 350, "test_"+btagstring);

  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && njets==3")&&util&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "MET", 6, 50, 350, "MET_njets3_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && njets==4")&&util&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "MET", 6, 50, 350, "MET_njets4_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && njets>=5")&&util&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "MET", 6, 50, 350, "MET_njetsge5_"+btagstring);

}


void drawVJets() {
  
  useFlavorHistoryWeights_=true;

  doOverflowAddition(true);
  setStackMode(true);
  doData(true);
  doRatioPlot(true);
  setLogY(false);
  
  setQuiet(true);

  int nbins;
  float low,high;
  TString var,xtitle;

  TCut baseSelection ="cutHT==1 && cutPV==1 && cutTrigger==1 && passInconsistentMuon==1 && passBadPFMuon==1 && cutDeltaPhi==1 && (nElectrons==0 && nMuons==1)";// && bestTopMass>350" ;
  TString extraName = ""; //"_m3jGT350";
  TCut njetCut = "njets ==2";  

  const TCut METcut = "MET>40";
  const TCut METcut2 = "MET>60";
  const TCut METcut3 = "MET>100";
  const TCut MTcut = "MT_Wlep>30";

  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretagCut =  "1";
  for (int ibtag = 0; ibtag<4; ibtag++) { //do this an ugly way for now
    TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
    if (ibtag==0) { //nothing to do
    }
    else if (ibtag==1) {
      theBTaggingCut = eq1b; 
      btagstring = "eq1b";
    }
    else if (ibtag==2) {
      theBTaggingCut = ge2b; 
      btagstring = "ge2b";
    }
    else if(ibtag==3) {
      theBTaggingCut = pretagCut;
      btagstring = "pre";
    }
    else assert(0);
   
    resetPlotMinimum();
    nbins=10; low=0; high = 10;
    var="njets"; xtitle="njets";
    

    selection_ = baseSelection && theBTaggingCut;
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_1muon_njets");


    selection_ = baseSelection && njetCut && theBTaggingCut;  

    resetPlotMinimum();
    nbins=10; low=0; high = 200;
    var="MET"; xtitle="MET (GeV)"; 
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_1muon2jet_met");
   
    resetPlotMinimum();
    nbins=20; low=10; high = 300;
    var="muonpt1"; xtitle="muonpt1"; 
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_1muon2jet_muonpt1");
    
    { //give this its own scope
      cout<<"no extra cuts"<<endl;
      double ndata= hdata->Integral();
      double nnonW=0,nW=0;
      double nnonW_err=0,nW_err=0;
      for ( std::map<TString, TH1D*>::iterator ihist = histos_.begin(); ihist!=histos_.end(); ++ihist) {
	if (ihist->second != 0) {
	  if (TString(ihist->second->GetName()).Contains("WJets")) {
	    nW = ihist->second->Integral();
	    nW_err = jmt::errOnIntegral(ihist->second);
	  }
	  else {
	    nnonW += ihist->second->Integral();
	    nnonW_err += pow(jmt::errOnIntegral(ihist->second),2);
	  }
	}
      }
      cout<<"(data - nonW)/W = ("<<ndata<<" - "<<nnonW<<") / "<<nW<<endl;
      double err = jmt::errAoverB( ndata-nnonW, sqrt(ndata+nnonW_err), nW,nW_err);
      cout<<"                = "<<(ndata-nnonW)/nW<<" +/- "<<err<<endl;
    }

    selection_ = baseSelection && njetCut && theBTaggingCut && METcut;  
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_1muon2jet_muonpt1_MET40");
    { //give this its own scope
      cout<<"MET>40"<<endl;
      double ndata= hdata->Integral();
      double nnonW=0,nW=0;
      double nnonW_err=0,nW_err=0;
      for ( std::map<TString, TH1D*>::iterator ihist = histos_.begin(); ihist!=histos_.end(); ++ihist) {
	if (ihist->second != 0) {
	  if (TString(ihist->second->GetName()).Contains("WJets")) {
	    nW = ihist->second->Integral();
	    nW_err = jmt::errOnIntegral(ihist->second);
	  }
	  else {
	    nnonW += ihist->second->Integral();
	    nnonW_err += pow(jmt::errOnIntegral(ihist->second),2);
	  }
	}
      }
      cout<<"(data - nonW)/W = ("<<ndata<<" - "<<nnonW<<") / "<<nW<<endl;
      double err = jmt::errAoverB( ndata-nnonW, sqrt(ndata+nnonW_err), nW,nW_err);
      cout<<"                = "<<(ndata-nnonW)/nW<<" +/- "<<err<<endl;
    }

    selection_ = baseSelection && njetCut && theBTaggingCut && METcut2;  
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_1muon2jet_muonpt1_MET60");
    { //give this its own scope
      cout<<"MET>60"<<endl;
      double ndata= hdata->Integral();
      double nnonW=0,nW=0;
      double nnonW_err=0,nW_err=0;
      for ( std::map<TString, TH1D*>::iterator ihist = histos_.begin(); ihist!=histos_.end(); ++ihist) {
	if (ihist->second != 0) {
	  if (TString(ihist->second->GetName()).Contains("WJets")) {
	    nW = ihist->second->Integral();
	    nW_err = jmt::errOnIntegral(ihist->second);
	  }
	  else {
	    nnonW += ihist->second->Integral();
	    nnonW_err += pow(jmt::errOnIntegral(ihist->second),2);
	  }
	}
      }
      cout<<"(data - nonW)/W = ("<<ndata<<" - "<<nnonW<<") / "<<nW<<endl;
      double err = jmt::errAoverB( ndata-nnonW, sqrt(ndata+nnonW_err), nW,nW_err);
      cout<<"                = "<<(ndata-nnonW)/nW<<" +/- "<<err<<endl;
    }

    selection_ = baseSelection && njetCut && theBTaggingCut && METcut3;  
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_1muon2jet_muonpt1_MET100");
    { //give this its own scope
      cout<<"MET>100"<<endl;
      double ndata= hdata->Integral();
      double nnonW=0,nW=0;
      double nnonW_err=0,nW_err=0;
      for ( std::map<TString, TH1D*>::iterator ihist = histos_.begin(); ihist!=histos_.end(); ++ihist) {
	if (ihist->second != 0) {
	  if (TString(ihist->second->GetName()).Contains("WJets")) {
	    nW = ihist->second->Integral();
	    nW_err = jmt::errOnIntegral(ihist->second);
	  }
	  else {
	    nnonW += ihist->second->Integral();
	    nnonW_err += pow(jmt::errOnIntegral(ihist->second),2);
	  }
	}
      }
      cout<<"(data - nonW)/W = ("<<ndata<<" - "<<nnonW<<") / "<<nW<<endl;
      double err = jmt::errAoverB( ndata-nnonW, sqrt(ndata+nnonW_err), nW,nW_err);
      cout<<"                = "<<(ndata-nnonW)/nW<<" +/- "<<err<<endl;
    }


    selection_ = baseSelection && njetCut && theBTaggingCut && MTcut;  
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_1muon2jet_muonpt1_MT15");
    { //give this its own scope
      cout<<"MT>30"<<endl;
      double ndata= hdata->Integral();
      double nnonW=0,nW=0;
      double nnonW_err=0,nW_err=0;
      for ( std::map<TString, TH1D*>::iterator ihist = histos_.begin(); ihist!=histos_.end(); ++ihist) {
	if (ihist->second != 0) {
	  if (TString(ihist->second->GetName()).Contains("WJets")) {
	    nW = ihist->second->Integral();
	    nW_err = jmt::errOnIntegral(ihist->second);
	  }
	  else {
	    nnonW += ihist->second->Integral();
	    nnonW_err += pow(jmt::errOnIntegral(ihist->second),2);
	  }
	}
      }
      cout<<"(data - nonW)/W = ("<<ndata<<" - "<<nnonW<<") / "<<nW<<endl;
      double err = jmt::errAoverB( ndata-nnonW, sqrt(ndata+nnonW_err), nW,nW_err);
      cout<<"                = "<<(ndata-nnonW)/nW<<" +/- "<<err<<endl;
      cout<<"W / SM = "<<nW/(nW+nnonW)<<endl;
    }


    selection_ = baseSelection && njetCut && theBTaggingCut;  
    nbins=20; low=0; high=300;
    //nbins=25; low=-5; high=10;
    var="MT_Wlep"; xtitle="leptonic W M_{T}";
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_1muon2jet_MT");

 }

}

void countABCD() {

  double fitResult[]={52.6, 40.5, 14.1};
  double fitResultErr[]={15.1, 15.9, 8.0};

  loadSamples();

  const  TCut baseSelection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1"; //no MET, no minDeltaPhi, no cleaning, no b tag
  const  TCut passCleaning ="passInconsistentMuon == 1 && passBadPFMuon==1 && weight<1000"; //apply cleaning but without ECAL dead cells
  const  TCut passMinDeltaPhi = "cutDeltaPhi==1";
  const  TCut failMinDeltaPhi = "cutDeltaPhi==0";
  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";

  for (int ibtag = 0; ibtag<3; ibtag++) { //do this an ugly way for now
    TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
    if (ibtag==0) { //nothing to do
    }
    else if (ibtag==1) {
      theBTaggingCut = eq1b; 
      btagstring = "eq1b";
    }
    else if (ibtag==2) {
      theBTaggingCut = ge2b; 
      btagstring = "ge2b";
    }
    else assert(0);

    TTree* tree= dtree;//(TTree*) fdata->Get("reducedTree");
    TH1D dummyhist("dummyhist","",1,0,1e9); //kludge!
    dummyhist.Sumw2(); //not really needed for data

    cout<<" ---- "<<btagstring<<" ---- "<<endl;

    TCut METC = "MET>=150 && MET<5000";
    TCut theSIGSelection = baseSelection && passCleaning && passMinDeltaPhi && theBTaggingCut && METC;
    selection_ = theSIGSelection.GetTitle();
    tree->Project("dummyhist","HT",getCutString(1.,"",selection_, "",0).Data());
    cout<<"N_SIG = "<<dummyhist.Integral()<<endl;

    dummyhist.Reset();
    if (dummyhist.Integral() != 0) assert(0);

    TCut METD = "MET>=150 && MET<5000";
    TCut theDSelection = baseSelection && passCleaning && failMinDeltaPhi && theBTaggingCut && METD;
    selection_ = theDSelection.GetTitle();
    tree->Project("dummyhist","HT",getCutString(1.,"",selection_,"",0).Data());
    double nd = dummyhist.Integral();
    cout<<"N_D = "<<nd<<endl;

    TCut META = "MET>=100 && MET<150";
    TCut theASelection = baseSelection && passCleaning && failMinDeltaPhi && theBTaggingCut && META;
    selection_ = theASelection.GetTitle();
    tree->Project("dummyhist","HT",getCutString(1.,"",selection_,"",0).Data());
    double na = dummyhist.Integral();
    cout<<"N_A = "<<na<<endl;

    double R = nd/na;
    double Rerr = R*sqrt(1/nd + 1/na);
    cout<<"R   = "<<R<<" +/- "<<Rerr<<endl;

    double est=R*fitResult[ibtag];
    double estErr = est*sqrt( pow(Rerr/R,2) + pow(fitResultErr[ibtag]/fitResult[ibtag],2));
    cout<<"estimate = "<<est<<" +/- "<<estErr<<endl;

  }
}

void flvHistReweighting() {

  //we want to make a way to apply the k-factors 

  //can use the extraSelection part of the getCutString()

  dodata_=false;
  const TString samplename = "WJets"; //if this is changed to data, lumiScale_ should not be applied
  setQuiet(true);
  loadSamples();
  
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";

  TTree* tree=0;
  if (samplename=="data") {
    tree = dtree;//(TTree*) fdata->Get("reducedTree");
  }
  else {
    for (unsigned int isample=0; isample<samples_.size(); isample++) {
      if ( samples_[isample].Contains(samplename)) {
	tree = (TTree*) files_[currentConfig_][samples_[isample]]->Get("reducedTree");
      }
    }
  }
  if (tree==0) {cout<<"Something went wrong finding your sample!"<<endl; return;}

  gROOT->cd();

  TH1D dummyhist("dummyhist","",1,0,1e9); //the typical kludge to count events; 1 bin histogram with large range
  dummyhist.Sumw2();

  //nominal case
  tree->Project("dummyhist","HT",getCutString(lumiScale_,"",selection_,"",0).Data());
  float n_nom = dummyhist.GetBinContent(1);
  float nerr_nom = dummyhist.GetBinError(1);
  cout<<" nominal     = "<<n_nom<<" +/- "<<nerr_nom<<endl;
  dummyhist.Reset();

  //with flv hist
  tree->Project("dummyhist","HT",getCutString(lumiScale_,"flavorHistoryWeight",selection_,"",0).Data());
  float  n_fh = dummyhist.GetBinContent(1);
  float  nerr_fh = dummyhist.GetBinError(1);
  cout<<" full reweight    = "<<n_fh<<" +/- "<<nerr_fh<<endl;

  const int ncat=11;
  const float kfactors[]  = {2,2,1,1,2,1,2,1,2,1,1};
  const float kerr_up[]   = {1,1,1,1,1,1,1,1,1,1,0};
  const float kerr_down[] = {1,1,0.5,0.5,1,0.5,1,0.5,1,0.5,0};

  float total_up=0;
  float total_down=0;

  //loop over the variations and up/down
  for (int iupdown = 0; iupdown<=1; iupdown++ ) {
    for (int ivariations = 0; ivariations<ncat; ivariations++) {
      dummyhist.Reset();
      TString manualreweight="";
      for (int ii=0; ii<ncat; ii++) {
	float k = kfactors[ii];
	//vary only one at a time
	if (ivariations == ii) {
	  if (iupdown==0) k += kerr_up[ii];
	  else            k -= kerr_down[ii];
	}
	TString thisterm;
	thisterm.Form( "((flavorHistory==%d) * %f)", ii+1, k);
	manualreweight += thisterm;
	if (ii+1 != ncat) manualreweight += " + ";
      }
      tree->Project("dummyhist","HT",getCutString(lumiScale_,manualreweight,selection_,"",0).Data()); //this will be normalized *wrong* because it doesn't have the N/sum(k_i) factor
      float  n = dummyhist.GetBinContent(1);
      float  nerr = dummyhist.GetBinError(1);
      dummyhist.Reset();
      assert(0);//please double check that the following few lines of code have the weights and cuts handled correctly
      tree->Project("dummyhist","HT","1"); //weight of 1
      double unweighted = dummyhist.GetBinContent(1);
      dummyhist.Reset();
      tree->Project("dummyhist","HT",manualreweight); //weighted with just the k factors
      double reweighted = dummyhist.GetBinContent(1);

      n *= unweighted/reweighted;
      nerr *= unweighted/reweighted;

      TString uddesc= iupdown==0 ? " up " : "down";
      cout<<"reweight "<< uddesc<<" [" << ivariations+1<<"] = "<<n<<" +/- "<<nerr<<" ; "<<n-n_fh<<endl;
      if (iupdown==0) total_up += pow(n-n_fh,2);
      else total_down += pow(n-n_fh,2);
    }
  }

  cout<<" up  syst = "<<sqrt(total_up)<<endl;
  cout<<"down syst = "<<sqrt(total_down)<<endl;

}

void pdfUncertainties() {
  dodata_=false;
  const TString samplename = "WJetsZ2"; //if this is changed to data, lumiScale_ should not be applied
  setQuiet(true);
  loadSamples();

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";

  TTree* tree=0;
  if (samplename=="data") {
    tree = dtree;
  }
  else {
    for (unsigned int isample=0; isample<samples_.size(); isample++) {
      if ( samples_[isample] == samplename) {
	tree = (TTree*) files_[currentConfig_][samples_[isample]]->Get("reducedTree");
      }
    }
  }
  if (tree==0) {cout<<"Something went wrong finding your sample!"<<endl; return;}

  gROOT->cd();
  TString optfh= useFlavorHistoryWeights_ && samplename.Contains("WJets") ? "flavorHistoryWeight" : "";

  TH1D dummyhist("dummyhist","",1,0,1e9); //the typical kludge to count events; 1 bin histogram with large range
  dummyhist.Sumw2();
  //nominal case
  tree->Project("dummyhist","HT",getCutString(lumiScale_,optfh,selection_,"",0).Data());
  float n = dummyhist.GetBinContent(1);
  float nerr = dummyhist.GetBinError(1);

  TH1D HpdfUnc("HpdfUnc","",100,n - 3*nerr, n+3*nerr);
  for (int i=1; i<=44; i++) { //hard code that there are 44+1 pdf weights
    dummyhist.Reset(); //maybe not needed
    tree->Project("dummyhist","HT",getCutString(lumiScale_,optfh,selection_,"",i).Data());
    HpdfUnc.Fill(dummyhist.GetBinContent(1));
  }

  cout<<" nominal = "<<n<<" +/- "<<nerr<<endl;
  cout<<" w/pdf   = "<<HpdfUnc.GetMean()<<" +/- "<<HpdfUnc.GetRMS()<<endl;

}

void countQCDMC(){
  useFlavorHistoryWeights_=false;
  setStackMode(true);
  doData(true);
  doRatioPlot(true);
  setQuiet(false);
  doOverflowAddition(true);

  loadSamples();
 
  const TCut myBase ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1"; //no MET, no mdp, no b-tag

  const TCut mdpCut = "cutDeltaPhi==1";
  const TCut invertedmdpCut = "cutDeltaPhi==0";

  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjets == 0";

  //now for a flexible MET region
  for (int metCutLow = 125; metCutLow <=150; metCutLow+=25) {
    //int metCutHigh=metCutLow+50;
    for(int metCutHigh=175; metCutHigh <= 225; metCutHigh+=25){
      
      int metCutHighSB = metCutLow+50;

      TString SBmetCutString;
      SBmetCutString.Form("MET >= %d && MET < %d",metCutLow,metCutHighSB);
      TCut SBMETselection(SBmetCutString.Data());
      TString SIGmetCutString; 
      SIGmetCutString.Form("MET >= %d",metCutHigh);
      TCut SIGMETselection(SIGmetCutString.Data());

      
      for (int ibtag = 0; ibtag<3; ibtag++) { //for now, skip pre and antib
	TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
	if (ibtag==0) { //nothing to do
	}
	else if (ibtag==1) {
	  theBTaggingCut = eq1b; 
	  btagstring = "eq1b";
	}
	else if (ibtag==2) {
	  theBTaggingCut = ge2b; 
	  btagstring = "ge2b";
	}
	else if (ibtag==3) {
	  theBTaggingCut = pretag;
	  btagstring = "pre";
	}
	else if (ibtag==4) {
	  theBTaggingCut = antitag;
	  btagstring = "antib";
	}
	else assert(0);
	
	double nA, nD, nSB, nSIG, nABCD, c_bias;
	double nA_err, nD_err, nSB_err, nSIG_err, nABCD_err, c_bias_err;
	TCut Aselection = myBase && theBTaggingCut && invertedmdpCut && SBMETselection;
	TCut Dselection = myBase && theBTaggingCut && invertedmdpCut && SIGMETselection;
	TCut SBselection = myBase && theBTaggingCut && mdpCut && SBMETselection;    
	TCut SIGselection = myBase && theBTaggingCut && mdpCut && SIGMETselection;
	
	selection_=Aselection.GetTitle();
	drawPlots("HT",1,0,10000,"","","deleteme");
	nA = totalqcd->GetBinContent(1);
	nA_err = totalqcd->GetBinError(1);
	
	selection_=Dselection.GetTitle();
	drawPlots("HT",1,0,10000,"","","deleteme");
	nD = totalqcd->GetBinContent(1);
	nD_err = totalqcd->GetBinError(1);
	
	selection_=SBselection.GetTitle();
	drawPlots("HT",1,0,10000,"","","deleteme");
	nSB = totalqcd->GetBinContent(1);
	nSB_err = totalqcd->GetBinError(1);
	
	selection_=SIGselection.GetTitle();
	drawPlots("HT",1,0,10000,"","","deleteme");
	nSIG = totalqcd->GetBinContent(1);
	nSIG_err = totalqcd->GetBinError(1);
	
	
	double up, down, up_err, down_err;
	nABCD = nD*nSB/nA;
	up = nD*nSIG;
	up_err = jmt::errAtimesB(nD, nD_err, nSB, nSB_err);
	down = nA;
	down_err = nA_err;
	nABCD_err = jmt::errAoverB(up, up_err, down, down_err);
	
	c_bias = nSIG*nA/(nSB*nD);
	up = nSIG*nA;
	up_err = jmt::errAtimesB(nSIG, nSIG_err, nA, nA_err);
	down = nSB*nD;
	down_err = jmt::errAtimesB(nSB, nSB_err, nD, nD_err);
	c_bias_err = jmt::errAoverB(up, up_err, down, down_err);
	
	double pSys, mSys;
	pSys = nABCD*(c_bias+fabs(1.0-c_bias)/2.0) - nABCD*c_bias;
	mSys = nABCD*c_bias - nABCD*(c_bias-fabs(1-c_bias)/2.0);
	

	cout << "countQCDMC: ******************************************" << endl;
	cout << "countQCDMC: SB is " << metCutLow << " GeV <= MET < " << metCutHighSB << " GeV, and SIG is MET >= " << metCutHigh << " GeV" << endl;
	cout << "countQCDMC: " << theBTaggingCut.GetTitle() << endl;
	cout << "countQCDMC: " << metCutLow <<"-"<<metCutHighSB<<" >" << metCutHigh << " " << nA << "+-" << nA_err << " " << nD << "+-" << nD_err << " " << nSB << " +-" << nSB_err << " " << nSIG << "+-" << nSIG_err << " " << nABCD << "+-" << nABCD_err << " " << c_bias << "+-" << c_bias_err << " +" << pSys <<"-"<<mSys << endl;
	cout << "countQCDMC: ******************************************" << endl;
	
      }//end b-tag loop
    }//end met high loop
  }//end met low loop
}// end countQCDMC()


void countILV() {
  setQuiet(true);
  loadSamples();
  const  TCut baseSelection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutDeltaPhi==1 && cutMET==1"; //no cleaning, no b tag
  const  TCut passCleaning ="passInconsistentMuon == 1 && passBadPFMuon==1 && weight<1000"; //apply cleaning but without ECAL dead cells

  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";

  const TCut muVeto = "nMuons==0";
  const TCut eVeto = "nElectrons==0";

  const TCut IMV = "nMuons==1";
  const TCut IEV = "nElectrons==1";

  //want to compare SM MC to data
  TCut IMVselection = baseSelection && passCleaning && eVeto && IMV;
  TCut IEVselection = baseSelection && passCleaning && muVeto && IEV;

  TCut IMV_ge1b = IMVselection && ge1b;
  TCut IMV_eq1b = IMVselection && eq1b;
  TCut IMV_ge2b = IMVselection && ge2b;

  TCut IEV_ge1b = IEVselection && ge1b;
  TCut IEV_eq1b = IEVselection && eq1b;
  TCut IEV_ge2b = IEVselection && ge2b;

  doOverflowAddition(false);
  bool oldSaveSetting = savePlots_;
  savePlots_=false;

  selection_=IMV_ge1b.GetTitle();
  drawPlots("HT",1,0,10000,"","","deleteme");
  cout<<"[IMV ge1b] Integral of data, total SM: "<<hdata->GetBinContent(1)<<" ; "<<totalsm->GetBinContent(1)<< " +/- "<<totalsm->GetBinError(1)<<endl;

  selection_=IMV_eq1b.GetTitle();
  drawPlots("HT",1,0,10000,"","","deleteme");
  cout<<"[IMV eq1b] Integral of data, total SM: "<<hdata->GetBinContent(1)<<" ; "<<totalsm->GetBinContent(1)<< " +/- "<<totalsm->GetBinError(1)<<endl;

  selection_=IMV_ge2b.GetTitle();
  drawPlots("HT",1,0,10000,"","","deleteme");
  cout<<"[IMV ge2b] Integral of data, total SM: "<<hdata->GetBinContent(1)<<" ; "<<totalsm->GetBinContent(1)<< " +/- "<<totalsm->GetBinError(1)<<endl;

  selection_=IEV_ge1b.GetTitle();
  drawPlots("HT",1,0,10000,"","","deleteme");
  cout<<"[IEV ge1b] Integral of data, total SM: "<<hdata->GetBinContent(1)<<" ; "<<totalsm->GetBinContent(1)<< " +/- "<<totalsm->GetBinError(1)<<endl;

  selection_=IEV_eq1b.GetTitle();
  drawPlots("HT",1,0,10000,"","","deleteme");
  cout<<"[IEV eq1b] Integral of data, total SM: "<<hdata->GetBinContent(1)<<" ; "<<totalsm->GetBinContent(1)<< " +/- "<<totalsm->GetBinError(1)<<endl;

  selection_=IEV_ge2b.GetTitle();
  drawPlots("HT",1,0,10000,"","","deleteme");
  cout<<"[IEV ge2b] Integral of data, total SM: "<<hdata->GetBinContent(1)<<" ; "<<totalsm->GetBinContent(1)<< " +/- "<<totalsm->GetBinError(1)<<endl;

  savePlots_=oldSaveSetting;
}

void drawOwen(TCut extracut="") {

  /*
.L drawReducedTrees.C++
  */
  
  //can't handle a tighter MET cut this way, because it affects the method differently than, say, a tighter HT cut
  assert( !TString(extracut.GetTitle()).Contains("MET"));

  savePlots_=false; //don't save eps,png,pdf files

  loadSamples(false); //make sure to load single top as 3 pieces
  //this only works if loadSamples() hasn't been called yet in this root session!

  doOverflowAddition(false);

  const  int nbins = 35;
  const  float min=0;
  const  float max=800;

  bool vb = true;//variable binning ON or OFF
  const int nvarbins=5;
  const float varbins[]={0.,160.,180.,260.,400.,800.};
  //use like this drawSimple("bestTopMass",nvarbins,varbins,histfilename, "bestM3j_met_150_10000_"+btagstring+"tag_qcd","PythiaPUQCDFlat");


  const  TCut baseSelection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1"; //no MET, no minDeltaPhi, no cleaning, no b tag
  // "cutCleaning==1"; //apply all cleaning 
  const  TCut passCleaning ="passInconsistentMuon == 1 && passBadPFMuon==1 && weight<1000"; //apply cleaning but without ECAL dead cells.  Notice the weight cut!
  const  TCut passMinDeltaPhi = "cutDeltaPhi==1";
  const  TCut failMinDeltaPhi = "cutDeltaPhi==0";
  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";

  //
  TString histfilename;
  histfilename.Form("bestM3j-%s.%s.root", vb ? "vb" : "bins", TString(extracut.GetTitle())=="" ? "baseline" : jmt::fortranize(extracut.GetTitle()).Data());
  TString textfilename=  "drawOwen.output";

  ofstream ofile(textfilename.Data());

  for (int ibtag = 0; ibtag<3; ibtag++) { //do this an ugly way for now
    TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
    if (ibtag==0) { //nothing to do
    }
    else if (ibtag==1) {
      theBTaggingCut = eq1b; 
      btagstring = "eq1b";
    }
    else if (ibtag==2) {
      theBTaggingCut = ge2b; 
      btagstring = "ge2b";
    }
    else assert(0);
    ofile<<" == "<<btagstring<<" == "<<endl;
 
   //these regions define the PDFs (templates)
    const TCut LSBMET = "MET>=0 && MET<50";
    TCut theLSBSelection = baseSelection && passCleaning && passMinDeltaPhi && theBTaggingCut && LSBMET && extracut;
    selection_ = theLSBSelection.GetTitle();
    if(vb){
     drawSimple("bestTopMass",nvarbins,varbins,histfilename, "bestM3j_met_0_50_"+btagstring+"tag_data","data");
     drawPlots("bestTopMass",nvarbins,varbins,"","","deleteme");
    }
    else{
      drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_0_50_"+btagstring+"tag_data","data");
      drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
    }
    TFile fh(histfilename,"UPDATE");  
    totalsm->SetName("bestM3j_met_0_50_"+btagstring+"tag_allsm");
    totalsm->Write();          
    fh.Close();   
        
    TCut T2MET = "MET>=50 && MET<5000";
    const  TCut baseT2Selection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1) || (nElectrons==1 && nMuons==0))"; //no MET, no minDeltaPhi, no cleaning, no b tag
    TCut theT2Selection = baseT2Selection && passCleaning && passMinDeltaPhi && theBTaggingCut && T2MET && extracut;
    selection_ = theT2Selection.GetTitle();
    if (vb) {
      drawSimple("bestTopMass",nvarbins,varbins,histfilename, "bestM3j_met_50_5000_t2_"+btagstring+"tag_data","data");
      drawPlots("bestTopMass",nvarbins,varbins,"","","deleteme");
    }
    else{
      drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_50_5000_t2_"+btagstring+"tag_data","data");
      drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
    }
    TFile fh2(histfilename,"UPDATE");  
    totalsm->SetName("bestM3j_met_50_5000_t2_"+btagstring+"tag_allsm");
    totalsm->Write();          
    fh2.Close(); 

    //need T2 in data,ttbar in the SR
    const TCut SRMET = "MET>=150 && MET<5000";
    theT2Selection = baseT2Selection && passCleaning && passMinDeltaPhi && theBTaggingCut && SRMET && extracut;
    selection_ = theT2Selection.GetTitle();
    if (vb) {
      drawSimple("bestTopMass",nvarbins,varbins,histfilename, "bestM3j_met_150_5000_t2_"+btagstring+"tag_ttbar","TTbarJets");
      drawSimple("bestTopMass",nvarbins,varbins,histfilename, "bestM3j_met_150_5000_t2_"+btagstring+"tag_data","data");
      drawPlots("bestTopMass",nvarbins,varbins,"","","deleteme");
    }
    else{
      drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_150_5000_t2_"+btagstring+"tag_ttbar","TTbarJets");
      drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_150_5000_t2_"+btagstring+"tag_data","data");
      drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
    }
    TFile fh22(histfilename,"UPDATE");  
    totalnonttbar->SetName("bestM3j_met_150_5000_t2_"+btagstring+"tag_nonttbar");
    totalewk->SetName("bestM3j_met_150_5000_t2_"+btagstring+"tag_allewk");
    totalnonttbar->Write();          
    totalewk->Write();          
    fh22.Close(); 
    
    // == signal region ==
    TCut theSRSelection = baseSelection && passCleaning && passMinDeltaPhi && theBTaggingCut && SRMET && extracut;
    selection_ = theSRSelection.GetTitle();
    if (vb) {
      //plot all samples
      for (unsigned int isample=0; isample<samples_.size(); isample++) {
	TString oname=sampleOwenName_[samples_[isample]];
	drawSimple("bestTopMass",nvarbins,varbins,histfilename, "bestM3j_met_150_5000_"+btagstring+"tag_"+oname,samples_[isample]);
      }
    }
    else {
      //plot all samples
      for (unsigned int isample=0; isample<samples_.size(); isample++) {
	TString oname=sampleOwenName_[samples_[isample]];
	drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_150_5000_"+btagstring+"tag_"+oname,samples_[isample]);
      }
    }

    //need invdphi in data in the SR
    TCut invdpSelection = baseSelection && passCleaning && failMinDeltaPhi && theBTaggingCut && SRMET && extracut;
    TString nameOfIDPhist = "bestM3j_met_150_5000_invdphi_";
    nameOfIDPhist+=btagstring; nameOfIDPhist+="tag_";
    selection_ = invdpSelection.GetTitle();
    if (vb) {
      drawSimple("bestTopMass",nvarbins,varbins,histfilename, nameOfIDPhist+"qcd" ,"PythiaPUQCD");
      drawSimple("bestTopMass",nvarbins,varbins,histfilename, nameOfIDPhist+"data" ,"data");
      drawPlots("bestTopMass",nvarbins,varbins,"","","deleteme");
    }
    else {
      drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfIDPhist+"qcd" ,"PythiaPUQCD");
      drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfIDPhist+"data" ,"data");
      drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
    }
    TFile fh3(histfilename,"UPDATE");
    totalnonqcd->SetName(nameOfIDPhist+"nonqcd");
    totalnonqcd->Write();
    fh3.Close();

    //now for a flexible MET region
    for (int metCutLow = 80; metCutLow <=100; metCutLow+=10) {
      for (int metCutHigh = 150; metCutHigh <=150; metCutHigh+=5) {
	TString metCutString; metCutString.Form("MET >= %d && MET < %d",metCutLow,metCutHigh);
	ofile<<metCutString<<endl;
	TCut METselection(metCutString.Data());

	theT2Selection = baseT2Selection && passCleaning && passMinDeltaPhi && theBTaggingCut && METselection && extracut;
	TString nameOfT2Hist;
	nameOfT2Hist.Form( "bestM3j_met_%d_%d_t2_%stag_",metCutLow,metCutHigh,btagstring.Data());

	TCut theSelection = baseSelection && passCleaning && passMinDeltaPhi && theBTaggingCut && METselection && extracut;
	selection_ = theSelection.GetTitle();
	TString nameOfHist;
	nameOfHist.Form( "bestM3j_met_%d_%d_%stag_",metCutLow,metCutHigh,btagstring.Data());
	
	//*** nominal selection ***//
	if (vb) {
	  //plot all samples
	  for (unsigned int isample=0; isample<samples_.size(); isample++) {
	    TString oname=sampleOwenName_[samples_[isample]];
	    ofile<<oname<<" = "<<drawSimple("bestTopMass",nvarbins,varbins,histfilename, nameOfHist+oname,samples_[isample])<<endl;
	  }
	  //plot data
	  drawSimple("bestTopMass",nvarbins,varbins,histfilename, nameOfHist+"data","data");
	  //fills plots that are combinations of various samples (to be accessed via global pointers)
	  drawPlots("bestTopMass",nvarbins,varbins,"","","deleteme");
	}
	else {
	  //plot all samples
	  for (unsigned int isample=0; isample<samples_.size(); isample++) {
	    TString oname=sampleOwenName_[samples_[isample]];
	    ofile<<oname<<" = "<<drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfHist+oname,samples_[isample])<<endl;
	  }
	  //plot data
	  drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfHist+"data","data");
	  //fills plots that are combinations of various samples (to be accessed via global pointers)
	  drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
	}
	TFile fh4(histfilename,"UPDATE");
	totalnonttbar->SetName(nameOfHist+"nonttbar");
	totalnonqcd->SetName(nameOfHist+"nonqcd");
	totalewk->SetName(nameOfHist+"allewk");
	totalsm->SetName(nameOfHist+"allsm");
	totalnonttbar->Write();
	totalnonqcd->Write();
	totalewk->Write();
	totalsm->Write();
	fh4.Close();
	
	//*** T2 selection ***//
	selection_ = theT2Selection.GetTitle(); //change to t2 selection
	if (vb) {
	  drawSimple("bestTopMass",nvarbins,varbins,histfilename, nameOfT2Hist+"data","data");
	  drawSimple("bestTopMass",nvarbins,varbins,histfilename, nameOfT2Hist+"ttbar","TTbarJets");
	  drawPlots("bestTopMass",nvarbins,varbins,"","","deleteme");
	}
	else {
	  drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfT2Hist+"data","data");
	  drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfT2Hist+"ttbar","TTbarJets");
	  drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
	}
	TFile fh44(histfilename,"UPDATE");  
	totalnonttbar->SetName(nameOfT2Hist+"nonttbar");
	totalewk->SetName(nameOfT2Hist+"allewk");
	totalnonttbar->Write();          
	totalewk->Write();          
	fh44.Close(); 

	//*** fail mdp selection ***//
      	theSelection = baseSelection && passCleaning && failMinDeltaPhi && theBTaggingCut && METselection && extracut;
	selection_ = theSelection.GetTitle();
	nameOfHist.Form( "bestM3j_met_%d_%d_invdphi_%stag_",metCutLow,metCutHigh,btagstring.Data());
	if(vb){
	  drawSimple("bestTopMass",nvarbins,varbins,histfilename, nameOfHist+"data","data");
	  drawSimple("bestTopMass",nvarbins,varbins,histfilename, nameOfHist+"qcd","PythiaPUQCD");
	  drawPlots("bestTopMass",nvarbins,varbins,"","","deleteme");
	}
	else{
	  drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfHist+"data","data");
	  drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfHist+"qcd","PythiaPUQCD");
	  drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
	}
	TFile fh5(histfilename,"UPDATE");
	totalnonqcd->SetName(nameOfHist+"nonqcd");
	totalsm->SetName(nameOfHist+"allsm");
	totalnonqcd->Write();
	totalsm->Write();
	fh5.Close();
      }
    }
     
  }
  ofile.close();

}

