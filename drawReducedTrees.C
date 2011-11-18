/*
to compile:

must have symlink to MiscUtil.cxx in the working directory.
TSelectorMultiDraw and CrossSectionTable need to be compiled only when they are changed.

gSystem->Load("TSelectorMultiDraw_C.so");
gSystem->Load("CrossSectionTable_cxx.so");
.L drawReducedTrees.C++

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

//from Luke
#include "TSelectorMultiDraw.h"

//container for mSugra cross-sections
#include "CrossSectionTable.h"

// can be checked out of UserCode/joshmt
#include "MiscUtil.cxx"

#include <fstream>

#include <iostream>
#include <map>
#include <set>
	

//*** AFTER SUMMER
//***************************

//-- reducedTrees for Oct 25 SUSY meeting. 3464.581/pb. 
TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-35d/";
TString dataInputPath =  "/cu2/ra2b/reducedTrees/V00-02-35c/";



//*** SUMMER RESULT
//***************************
//for rerunning the final background estimates...need standard MC and MET cleaning
//V00-02-24_fullpf2pat was used for 13 Aug update of data results
//stick with this version for now because 25 is missing one ZJets sample and that takes away one event from the LDP
//TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-24_fullpf2pat/";//path for MC
//TString dataInputPath = "/cu3/wteo/reducedTrees/V00-02-05_v3-pickevents/"; //includes MET cleaning but uses a tight skim (not good for plots)

//for making the standard set of plots...need standard MC and all data.
//TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-25_fullpf2pat/";//path for MC	     
//TString dataInputPath = "/cu2/ra2b/reducedTrees/V00-02-25_fullpf2pat/";

//for njet reweighting (needed because not all reducedTrees in nominal V00-02-25 have njets30)
//TString inputPath = "/cu2/ra2b/reducedTrees/benV00-02-25_fullpf2pat/";
//TString dataInputPath = "/cu2/ra2b/reducedTrees/benV00-02-25_fullpf2pat/";

//for signal systematics
  // TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-25c_fullpf2pat/"; //LM9 with correct pdf weights
//TString inputPath = "/cu2/joshmt/reducedTrees/V00-02-25c_fullpf2pat/"; //with correct pdf weights
//TString inputPath = "/home/joshmt/";//path for MC
//TString dataInputPath = "/cu2/ra2b/reducedTrees/V00-02-24_fullpf2pat/"; //sym links to V00-02-05_v3

//the cutdesc string is now defined in loadSamples()

//double lumiScale_ = 1091.891;
  //double lumiScale_ = 1143; //official summer conf lumi
//double lumiScale_ = 3464.581;//oct25
double lumiScale_ = 4683.719;//nov4

#include "drawReducedTrees.h"


//this is the original implementation, using the traditional method of drawing a 1 bin histo to count events
//works fine for small samples, completely unusable for mSugra
SignalEffData signalSystematics2011(const SearchRegion & region, bool isSL=false, bool isLDP=false, TString sampleOfInterest="LM9", int m0=0,int m12=0) {
  useFlavorHistoryWeights_=false;
  loadSamples(); //needs to come first! this is why this should be turned into a class, with this guy in the ctor

  assert(configDescriptions_.size() == 12); //to avoid stupid mistakes on LM9 (no more HLT eff)

  dodata_=false;
  savePlots_ =false;
  //setQuiet(true);

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
  TCut bcut = "nbjets>=1";
  TString bweightstring="probge1";
  if (btagselection=="ge2b") {
    bcut="nbjets>=2";
    bweightstring="probge2";
  }
  else if (btagselection=="ge1b") {}
  else if (btagselection=="ge3b") {
    bcut="nbjets>=3";
    bweightstring="probge3"; //this one isn't calculated right now, i think
  }
  else {assert(0);}

  TCut baseline = "cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1";
  if (isSL) baseline = "cutPV==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0))";

  //  assert(0); //need to decide if SL should have MT cut

  TCut passOther = "minDeltaPhiN>=4";
  if (isLDP) passOther="minDeltaPhiN<4";

  selection_ = baseline && HTcut && passOther && SRMET && bcut;

  //critical to reset these!
  usePUweight_=false;
  useHLTeff_=false;
  btagSFweight_="1";
  currentConfig_=configDescriptions_.getDefault(); //completely raw MC

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

  currentConfig_=configDescriptions_.getCorrected(); //add  JER bias
  //the logic here used to be very fragile, so try to ensure that this is the correct thing
  //(things are improved a bit now, but still a good check to make)
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

  //do NOT do everything in %
  kfactorplus=  (kfactorplus-nominal)/nominal;
  kfactorminus= (kfactorminus-nominal)/nominal;
  cout<<"k factor uncertainties = "<<kfactorminus<<" "<<kfactorplus<<endl;
  results.set("kFactor",kfactorminus,kfactorplus);

  susyCrossSectionVariation_=""; //reset to default k factors
  //all other uncertainties

  //we're gonna hard code the fact that these variations come in pairs
  double var1=0,var2=0;
  for (unsigned int j=2; j<configDescriptions_.size(); j++) {
    currentConfig_=configDescriptions_.at(j);
    drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
    double thisn=getIntegral(sampleOfInterest);
    if (j%2==0) { //this one runs first
      var2=(thisn-nominal)/nominal;
    }
    else { //this one runs second
      var1=(thisn-nominal)/nominal;
      results.set(configDescriptions_.getVariedSubstring(currentConfig_), var1,var2);
    }
  }

  //reset to the nominal with JER bias
  currentConfig_=configDescriptions_.getCorrected();

  double largest=0;
  double smallest=1e9;

  //debug
//   TH1D*  HyieldsCTEQ = new TH1D("HyieldsCTEQ","CTEQ yields",100,0,200);
//   TH1D*  HyieldsMSTW = new TH1D("HyieldsMSTW","MSTW yields",100,0,200);
//   TH1D*  HyieldsNNPDF= new TH1D("HyieldsNNPDF","NNPDF yields",100,0,200);

  double percentCTEQ=0,percentMSTW=0,percentNNPDF=0;

  if (true) { //do do pdf uncertainties!
    
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
    if (sampleOfInterest!="T1bbbb"&&sampleOfInterest!="T2bb"&&sampleOfInterest!="T2tt")    thetree->Project("pdfEventCounter","HT",getCutString(false, "",extraWeight,i,"CTEQ").Data());
    else      thetree->Project("pdfEventCounter","HT",getCutString(kSMSPoint,"",selection_,"",i,"CTEQ").Data());
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

  cout<<"CTEQ\t"<<XminusMax/scale<<"\t"<<XplusMax/scale<<endl;
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
    if (sampleOfInterest!="T1bbbb"&&sampleOfInterest!="T2bb"&&sampleOfInterest!="T2tt")    thetree->Project("pdfEventCounter","HT",getCutString(false, "",extraWeight,i,"MSTW").Data());
    else      thetree->Project("pdfEventCounter","HT",getCutString(kSMSPoint,"",selection_,"",i,"MSTW").Data());
    //    thetree->Project("pdfEventCounter","HT",getCutString(false, "",extraWeight,i,"MSTW").Data());
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
  cout<<"MSTW\t"<<XminusMax<<"\t"<<XplusMax<<endl;
  percentMSTW = 0.5* ( XminusMax + XplusMax);
  percentMSTW = 100* percentMSTW/nominal;

  }
    
  // ~~~~~~~~~~~~~~~~ NNPDF ~~~~~~~~~~~(see emails from H&S -- just take the RMS as the error)
  {
    double totalX=0,totalX2=0;
    int n=0;
    TH1D nnpdfYields("nnpdfYields","nnpdf yields",2,0,1e9); //the limits really don't matter (root finds the RMS in an unbinned way)
    for (int i=1; i<=99; i++) { //0 is just weight 1. there are 99 others
      //get the sum of weight[i] for the whole sample with no cuts
      float totalweight=   totalPdfWeightsNNPDF->GetBinContent(i+1); //i really did store the total in bin i+1
      TString extraWeight=""; extraWeight+=nominalweight; extraWeight+="/"; extraWeight+=totalweight;
      //then get the yield after cuts with extra weight pdfWeights[i]/sum
      TH1D pdfEventCounter("pdfEventCounter","pdfEventCounter",1,0,1e9);
      pdfEventCounter.Sumw2();
      //      thetree->Project("pdfEventCounter","HT",getCutString(false, "",extraWeight,i,"NNPDF").Data());
      if (sampleOfInterest!="T1bbbb"&&sampleOfInterest!="T2bb"&&sampleOfInterest!="T2tt")    thetree->Project("pdfEventCounter","HT",getCutString(false, "",extraWeight,i,"NNPDF").Data());
      else      thetree->Project("pdfEventCounter","HT",getCutString(kSMSPoint,"",selection_,"",i,"NNPDF").Data());
      cout<<"    yield = "<<pdfEventCounter.Integral()<<endl;
      nnpdfYields.Fill(pdfEventCounter.Integral());
      totalX += pdfEventCounter.Integral();
      totalX2 += pdfEventCounter.Integral()*pdfEventCounter.Integral();
      n++;
      //      HyieldsNNPDF->Fill(pdfEventCounter.Integral());
      pdfEventCounter.Reset();
    }
    totalX/= double(n);
    totalX2/=double(n);
    cout<<"NNPDF\t"<<sqrt(totalX2-(totalX*totalX))<<endl;   
    cout<<"NNPDF (new)\t"<<nnpdfYields.GetRMS() <<endl;
    //    percentNNPDF = 100*sqrt(totalX2-(totalX*totalX)) / nominal;
    percentNNPDF = 100*nnpdfYields.GetRMS() / nominal;
  }

  if (percentCTEQ > percentMSTW) percentMSTW=percentCTEQ;
  if (percentNNPDF > percentMSTW) percentMSTW = percentNNPDF;
  cout<<"PDF uncertainty (%) = "<<percentMSTW<<endl;

  results.set("PDF",percentMSTW*0.01, percentMSTW*0.01); //not in %
  }
  else { 
    //    results.totalSystematic += 13*13;
    cout<<"Skipping PDF uncertainties!"<<endl;
  }

  return results;
}

/* example of how to use the new file-writing feature
void runSystematics2011_LM9_test1( TString sampleOfInterest="LM9" ) {
  setSearchRegions();

  for (ULong_t i=0; i<searchRegions_.size(); i++) {
    SignalEffData SB =  signalSystematics2011(sbRegions_[i],false,false,sampleOfInterest);
    SignalEffData SIG=  signalSystematics2011(searchRegions_[i],false,false,sampleOfInterest);

    SignalEffData SBSL= signalSystematics2011(sbRegions_[i],true,false,sampleOfInterest);
    SignalEffData SIGSL=signalSystematics2011(searchRegions_[i],true,false,sampleOfInterest);

    SignalEffData SBLDP= signalSystematics2011(sbRegions_[i],false,true,sampleOfInterest);
    SignalEffData SIGLDP=signalSystematics2011(searchRegions_[i],false,true,sampleOfInterest);

    //using ULong_t above because it aids the compiler in this concatenation
    SB.write(TString("SB")+i);
    SIG.write(TString("SIG")+i);

    SBSL.write(TString("SBSL")+i);
    SIGSL.write(TString("SIGSL")+i);

    SBLDP.write(TString("SBLDP")+i);
    SIGLDP.write(TString("SIGLDP")+i);
  }
}

void runSystematics2011_LM9_test2( TString sampleOfInterest="LM9" ) {
  setSearchRegions();

  for (ULong_t i=0; i<searchRegions_.size(); i++) {
    char effoutput[500];
    sprintf(effoutput,"signalSyst.%s.%s%s.dat",sampleOfInterest.Data(),searchRegions_[i].btagSelection.Data(),searchRegions_[i].owenId.Data());

    ofstream    textfiles(effoutput);
    textfiles<<m0_<<" "<<m12_<<" "<<0<<" ";
    SignalEffData SB(TString("SB")+i);
    SignalEffData SIG(TString("SIG")+i);
    SignalEffData SBSL(TString("SBSL")+i);
    SignalEffData SIGSL(TString("SIGSL")+i);

    SignalEffData SBLDP(TString("SBLDP")+i);
    SignalEffData SIGLDP(TString("SIGLDP")+i);

    textfiles<<SIG.rawYield<<" "<<SB.rawYield<<" "<<SIGSL.rawYield<<" "<<SBSL.rawYield<<" "<<SIGLDP.rawYield<<" "<<SBLDP.rawYield<<" "
	     <<SIG.effCorr<<" "<<SB.effCorr<<" "<<SIGSL.effCorr<<" "<<SBSL.effCorr<<" "<<SIGLDP.effCorr<<" "<<SBLDP.effCorr<<" "
	     <<SIG.totalSystematic()<<" "<<SB.totalSystematic()<<" "<<SIGSL.totalSystematic()<<" "<<SBSL.totalSystematic()<<" "<<SIGLDP.totalSystematic()<<" "<<SBLDP.totalSystematic()<<endl;
    textfiles.close();
  }
  
}
end example */

//not only for LM9; can also be used for other samples e.g. ttbar. Just change sampleOfInterest
void runSystematics2011_LM9(  TString sampleOfInterest="LM9" ) {
  setSearchRegions();

  //  signalSystematics2011(sbRegions_[0],false,false,"T1bbbb",1175,75);
  //     return;

  

  vector<ofstream*> textfiles;
  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    char effoutput[500];
    sprintf(effoutput,"signalSyst.%s.%s%s.dat",sampleOfInterest.Data(),searchRegions_[i].btagSelection.Data(),searchRegions_[i].owenId.Data());

    textfiles.push_back( new ofstream(effoutput));
  }

  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    (*textfiles[i])<<m0_<<" "<<m12_<<" "<<0<<" ";
    SignalEffData SB =  signalSystematics2011(sbRegions_[i],false,false,sampleOfInterest);
    SignalEffData SIG=  signalSystematics2011(searchRegions_[i],false,false,sampleOfInterest);

    SignalEffData SBSL= signalSystematics2011(sbRegions_[i],true,false,sampleOfInterest);
    SignalEffData SIGSL=signalSystematics2011(searchRegions_[i],true,false,sampleOfInterest);

    SignalEffData SBLDP= signalSystematics2011(sbRegions_[i],false,true,sampleOfInterest);
    SignalEffData SIGLDP=signalSystematics2011(searchRegions_[i],false,true,sampleOfInterest);

    (*textfiles[i])<<SIG.rawYield<<" "<<SB.rawYield<<" "<<SIGSL.rawYield<<" "<<SBSL.rawYield<<" "<<SIGLDP.rawYield<<" "<<SBLDP.rawYield<<" "
		   <<SIG.effCorr<<" "<<SB.effCorr<<" "<<SIGSL.effCorr<<" "<<SBSL.effCorr<<" "<<SIGLDP.effCorr<<" "<<SBLDP.effCorr<<" "
		   <<SIG.totalSystematic()<<" "<<SB.totalSystematic()<<" "<<SIGSL.totalSystematic()<<" "<<SBSL.totalSystematic()<<" "<<SIGLDP.totalSystematic()<<" "<<SBLDP.totalSystematic()<<endl;

//     cout<<" == Signed JES systematic in SB and SIG "<<endl;
//     cout<<"["<<100*SB.systematics["JES"].minus<<" "<<100*SB.systematics["JES"].plus<<"]"
// 	<<"\t"
// 	<<"["<<100*SIG.systematics["JES"].minus<<" "<<100*SIG.systematics["JES"].plus<<"]"
// 	<<endl;
//     cout<<" ====================================== "<<endl;
  }

  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    textfiles.at(i)->close();
  }

}


void countInBoxesBreakdown(const SearchRegion & region) {
  //shows number of events from each background and signal in the 6 boxes
  
  useFlavorHistoryWeights_=false;
  loadSamples();


  //dodata_=false;
  savePlots_ =false;
  setQuiet(true);
  doOverflowAddition(true);  
  //region.Print();
  
  //for susy scan; harmless for non-susy scan
  resetSamples();
  
  TString btagselection = region.btagSelection;
  TCut HTcut=region.htSelection.Data(); 
  TCut SRMET = region.metSelection.Data();
  bool isSIG = region.isSIG;

   
  TCut bcut = "1";
  if (btagselection=="ge1b") {
    bcut="nbjetsCSVM>=1";
  }
  else if (btagselection=="ge2b") {
    bcut="nbjetsCSVM>=2";
  }
  else if (btagselection=="ge3b") {
    bcut="nbjetsCSVM>=3";
  }
  else {assert(0);}
  
  usePUweight_=false;
  useHLTeff_=false;
  btagSFweight_="1";
  currentConfig_=configDescriptions_.getDefault(); //completely raw MC
  
  /*
  //block to cross check old cuts
  TString bweightstring="probge1";
  if (btagselection=="ge2b") {
    bweightstring="probge2";
  }
  else if (btagselection=="ge1b") {}
  else if (btagselection=="ge3b") {
    bweightstring="probge3"; //this one isn't calculated right now, i think
  }
  else {assert(0);}
  //use all of the corrections that are in the AN table
  currentConfig_=configDescriptions_.getCorrected(); //add  JER bias
  usePUweight_=true;
  useHLTeff_=true;
  btagSFweight_=bweightstring; //use b tag weight instead
  */
  
  TString thisbox="";
  
  TCut baseline = "cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1";
  TCut baselineSL = "cutPV==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0 && MT_Wlep<100";
  TCut passOther = "minDeltaPhiN>=4";
  TCut failOther="minDeltaPhiN<4";

  TString var="HT"; TString xtitle=var;
  int  nbins=1; float low=0; float high=1e9;
  char output[500];
  //Look at: PythiaPUQCD TTbarJets SingleTop WJets ZJets VV Zinvisible LM9
  
  //double QCD, QCDerr, ttbar, ttbarerr, singletop, singletoperr, wjets, wjetserr, zjets, zjetserr, vv, vverr, zinv, zinverr, lm9, lm9err, data, dataerr, totalsm, totalsmerr;

  // --- nominal (SB or SIG) selection ---
  selection_ = baseline && HTcut && passOther && SRMET && bcut; 
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  thisbox=""; thisbox+= isSIG ? "$N_{SIG}$":"$N_{SB}$";
  sprintf(output,"%s & %s%s & %s & %s & %s & %s & %s & %s & %s & %s& %i & %s \\\\",thisbox.Data(),region.btagSelection.Data(),region.owenId.Data(),
	  jmt::format_nevents(getIntegral("PythiaPUQCD"),getIntegralErr("PythiaPUQCD")).Data(), jmt::format_nevents(getIntegral("TTbarJets"),getIntegralErr("TTbarJets")).Data(),
	  jmt::format_nevents(getIntegral("SingleTop"),getIntegralErr("SingleTop")).Data(), jmt::format_nevents(getIntegral("WJets"),getIntegralErr("WJets")).Data(),
	  jmt::format_nevents(getIntegral("ZJets"),getIntegralErr("ZJets")).Data(), jmt::format_nevents(getIntegral("VV"),getIntegralErr("VV")).Data(),
	  jmt::format_nevents(getIntegral("Zinvisible"),getIntegralErr("Zinvisible")).Data(), 	  jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm")).Data(), 
	  (int)(getIntegral("data")+0.5), jmt::format_nevents(getIntegral("LM9"),getIntegralErr("LM9")).Data());
  cout<<output<<endl;  

  // -- SL selection --
  selection_ = baselineSL && HTcut && passOther && SRMET && bcut;
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  thisbox=""; thisbox+= isSIG ? "$N_{SIG,SL}$":"$N_{SB,SL}$";
  sprintf(output,"%s & %s%s & %s & %s & %s & %s & %s & %s & %s & %s& %i & %s \\\\",thisbox.Data(),region.btagSelection.Data(),region.owenId.Data(),
	  jmt::format_nevents(getIntegral("PythiaPUQCD"),getIntegralErr("PythiaPUQCD")).Data(), jmt::format_nevents(getIntegral("TTbarJets"),getIntegralErr("TTbarJets")).Data(),
	  jmt::format_nevents(getIntegral("SingleTop"),getIntegralErr("SingleTop")).Data(), jmt::format_nevents(getIntegral("WJets"),getIntegralErr("WJets")).Data(),
	  jmt::format_nevents(getIntegral("ZJets"),getIntegralErr("ZJets")).Data(), jmt::format_nevents(getIntegral("VV"),getIntegralErr("VV")).Data(),
	  jmt::format_nevents(getIntegral("Zinvisible"),getIntegralErr("Zinvisible")).Data(), 	  jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm")).Data(), 
	(int)(getIntegral("data")+0.5), jmt::format_nevents(getIntegral("LM9"),getIntegralErr("LM9")).Data());
  cout<<output<<endl; 
  
  
  // -- LDP selection --
  selection_ = baseline && HTcut && failOther && SRMET && bcut;
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  thisbox=""; thisbox+= isSIG ? "$N_{SIG,LDP}$":"$N_{SB,LDP}$";
  sprintf(output,"%s & %s%s & %s & %s & %s & %s & %s & %s & %s & %s& %i & %s \\\\",thisbox.Data(),region.btagSelection.Data(),region.owenId.Data(),
	  jmt::format_nevents(getIntegral("PythiaPUQCD"),getIntegralErr("PythiaPUQCD")).Data(), jmt::format_nevents(getIntegral("TTbarJets"),getIntegralErr("TTbarJets")).Data(),
	  jmt::format_nevents(getIntegral("SingleTop"),getIntegralErr("SingleTop")).Data(), jmt::format_nevents(getIntegral("WJets"),getIntegralErr("WJets")).Data(),
	  jmt::format_nevents(getIntegral("ZJets"),getIntegralErr("ZJets")).Data(), jmt::format_nevents(getIntegral("VV"),getIntegralErr("VV")).Data(),
	  jmt::format_nevents(getIntegral("Zinvisible"),getIntegralErr("Zinvisible")).Data(), 	  jmt::format_nevents(getIntegral("totalsm"),getIntegralErr("totalsm")).Data(), 
	  (int)(getIntegral("data")+0.5), jmt::format_nevents(getIntegral("LM9"),getIntegralErr("LM9")).Data());
  cout<<output<<endl; 
}




void runCountInBoxesBreakdown() {

  setSearchRegions();
  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    countInBoxesBreakdown(searchRegions_[i]) ;
    countInBoxesBreakdown(sbRegions_[i]) ;
  }

}


void countInBoxes(const SearchRegion & region,  ofstream * outfile ) {

  /*
owen asked for total MC event counts in the 6 boxes, output in a particular format
  */

  //obviously this event counting is coded elsewhere in many other forms. but let's just redo it

  //now that i know what owen really wants, i could have done this in the signal systematics code, but
  //now this is written....

  useFlavorHistoryWeights_=false;
  loadSamples(); //needs to come first! this is why this should be turned into a class, with this guy in the ctor

  dodata_=false;
  savePlots_ =false;
  //  setQuiet(true);

  region.Print();

  //for susy scan; harmless for non-susy scan
  resetSamples();
  const TString sample="totalsm";

  TString btagselection = region.btagSelection;
  TCut HTcut=region.htSelection.Data(); 
  TCut SRMET = region.metSelection.Data();
  bool isSIG = region.isSIG;
  TCut bcut = "nbjets>=1";
  TString bweightstring="probge1";
  if (btagselection=="ge2b") {
    bcut="nbjets>=2";
    bweightstring="probge2";
  }
  else if (btagselection=="ge1b") {}
  else if (btagselection=="ge3b") {
    bcut="nbjets>=3";
    bweightstring="probge3"; //this one isn't calculated right now, i think
  }
  else {assert(0);}

  TCut baseline = "cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1";
  TCut baselineSL = "cutPV==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0))";

  TCut passOther = "minDeltaPhiN>=4";
  TCut failOther="minDeltaPhiN<4";


  //use all of the corrections that are in the AN table
  currentConfig_=configDescriptions_.getCorrected(); //add  JER bias
  usePUweight_=true;
  useHLTeff_=true;
  btagSFweight_=bweightstring; //use b tag weight instead

  TString var="HT"; TString xtitle=var;
  int  nbins=1; float low=0; float high=1e9;

/*
Nsig              155
Nsb               244
Nsig_sl           103
Nsb_sl            165
Nsig_ldp          89
Nsb_ldp           393
*/

  // --- nominal selection ---
  selection_ = baseline && HTcut && passOther && SRMET; //no b cut because we're using the weighting
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  (*outfile)<<"N";//<<sampleOwenName_[sample]<<"mc"<<"_";
  if (isSIG) (*outfile)<<"sig            "; 
  else       (*outfile)<<"sb             "; 
  (*outfile)<<getIntegral(sample)<<endl;

  //now the error
//   (*outfile)<<"N"<<sampleOwenName_[sample]<<"mc"<<"_";
//   if (isSIG) (*outfile)<<"sig_err        "; 
//   else       (*outfile)<<"sb_err         "; 
//   (*outfile)<<getIntegralErr(sample)<<endl;

  // -- SL selection --
  selection_ = baselineSL && HTcut && passOther && SRMET; //no b cut because we're using the weighting
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  (*outfile)<<"N";//<<sampleOwenName_[sample]<<"mc"<<"_";
  if (isSIG) (*outfile)<<"sig_sl         "; 
  else       (*outfile)<<"sb_sl          "; 
  (*outfile)<<getIntegral(sample)<<endl;

  //now the error
//   (*outfile)<<"N"<<sampleOwenName_[sample]<<"mc"<<"_";
//   if (isSIG) (*outfile)<<"sig_sl_err     "; 
//   else       (*outfile)<<"sb_sl_err      "; 
//   (*outfile)<<getIntegralErr(sample)<<endl;

  // -- LDP selection --
  selection_ = baseline && HTcut && failOther && SRMET; //no b cut because we're using the weighting
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  (*outfile)<<"N";//<<sampleOwenName_[sample]<<"mc"<<"_";
  if (isSIG) (*outfile)<<"sig_ldp        "; 
  else       (*outfile)<<"sb_ldp         "; 
  (*outfile)<<getIntegral(sample)<<endl;

  //now the error
//   (*outfile)<<"N"<<sampleOwenName_[sample]<<"mc"<<"_";
//   if (isSIG) (*outfile)<<"sig_ldp_err    "; 
//   else       (*outfile)<<"sb_ldp_err     "; 
//   (*outfile)<<getIntegralErr(sample)<<endl;

  

}

void runCountInBoxesMC() {

  setSearchRegions();
  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    char effoutput[500];
    sprintf(effoutput,"mcEventCounts.%s%s.%dinvpb.dat",searchRegions_[i].btagSelection.Data(),searchRegions_[i].owenId.Data(),TMath::Nint(lumiScale_));
    ofstream  textfile(effoutput);

    countInBoxes(searchRegions_[i], &textfile ) ;
    countInBoxes(sbRegions_[i], &textfile ) ;
    textfile.close();
  }

}




//much more efficienct signal efficiency systematics calculation for mSugra and SMS
map<pair<int,int>, SignalEffData>  runTH2Syst2011_mSugra(const SearchRegion & region, 
                 const bool isSL=false, const bool isLDP=false, const TString & sampleOfInterest="mSUGRAtanb40") {

  //some of this is prob not needed, but to be safe
  useFlavorHistoryWeights_=false;
  dodata_=false;
 
  loadSamples();
  clearSamples();
  addSample(sampleOfInterest);

  setSearchRegions();
  const bool ismSugra= sampleOfInterest.Contains("mSUGRA");
  if (ismSugra)  loadSusyScanHistograms();


  region.Print();
  assert ( !( isSL&&isLDP));
  if (isLDP) cout<<"   == LDP "<<endl;
  if (isSL) cout<<"   == SL "<<endl;


  TString btagselection = region.btagSelection;
  TCut HTcut=region.htSelection.Data(); 
  TCut SRMET = region.metSelection.Data();
  TCut bcut = "nbjets>=1";
  TString bweightstring="probge1";
  if (btagselection=="ge2b") {
    bcut="nbjets>=2";
    bweightstring="probge2";
  }
  else if (btagselection=="ge1b") {}
  else if (btagselection=="ge3b") {
    bcut="nbjets>=3";
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
  currentConfig_=configDescriptions_.getDefault(); //completely raw MC

  //the logic here is very fragile to user error, so try to ensure that this is really completely raw MC
  assert( currentConfig_.Contains("JER0") && currentConfig_.Contains("JES0") && currentConfig_.Contains("METunc0")
	  && currentConfig_.Contains("PUunc0")&& currentConfig_.Contains("BTagEff0")&& currentConfig_.Contains("HLTEff0") );

  susyScanYields rawYields=getSusyScanYields(sampleOfInterest);

  currentConfig_=configDescriptions_.getCorrected(); //add  JER bias
  assert( currentConfig_.Contains("JERbias") && currentConfig_.Contains("JES0") && currentConfig_.Contains("METunc0")
	  && currentConfig_.Contains("PUunc0")&& currentConfig_.Contains("BTagEff0")&& currentConfig_.Contains("HLTEff0") );

  usePUweight_=true;
  useHLTeff_=true;
  selection_ = baseline && HTcut && passOther && SRMET; //remove b tag
  btagSFweight_=bweightstring; // add b tag back in the form of the SF
  susyScanYields nominalYields=getSusyScanYields(sampleOfInterest);

  //only applicable to mSugra
  //NLO k factor uncertainty -- vary it up
  susyCrossSectionVariation_="Plus";
  susyScanYields kFactorPlusYields= ismSugra ? getSusyScanYields(sampleOfInterest) : nominalYields;
  
  //NLO k factor uncertainty -- vary it down
  susyCrossSectionVariation_="Minus";
  susyScanYields kFactorMinusYields=ismSugra ? getSusyScanYields(sampleOfInterest): nominalYields;

  susyCrossSectionVariation_=""; //reset to default k factors

  // == do PDF uncertainties ==
  pair<susyScanYields,susyScanYields> cteqXplusminus;
  pair<susyScanYields,susyScanYields> mstwXplusminus; 
  pair<susyScanYields,susyScanYields> nnpdfXX2; 
  //  if (fullPdfUncertainties) {
  //    cout<<"Doing full pdf uncertainties!"<<endl;
    cteqXplusminus = doScanPDFUncertainties(sampleOfInterest, "CTEQ",nominalYields);
    mstwXplusminus = doScanPDFUncertainties(sampleOfInterest, "MSTW",nominalYields);
    nnpdfXX2 = doScanPDFUncertainties(sampleOfInterest,"NNPDF",nominalYields);
    //}

  //again, this is rather fragile. the configDescriptions_ list needs to be set just right
  //where just right is: variations in pairs, starting at the 3rd (index==2) guy in the list

  //make a map from the varied parameter to a pair of varied yields
  map<TString, pair<susyScanYields,susyScanYields> > variedYields;

  for (unsigned int j=2; j<configDescriptions_.size(); j+=2) {
    currentConfig_=configDescriptions_.at(j);
    TString varied1=  configDescriptions_.getVariedSubstring(currentConfig_);
    susyScanYields variation1Yields=getSusyScanYields(sampleOfInterest);

    currentConfig_=configDescriptions_.at(j+1);
    TString varied2=   configDescriptions_.getVariedSubstring(currentConfig_);
    susyScanYields variation2Yields=getSusyScanYields(sampleOfInterest);

    assert(varied1 == varied2);

    variedYields[varied1] = make_pair(variation1Yields,variation2Yields);
  }


  
// at this point we've got yields (raw , eff corrected , kfactor +/-, variations of JES, etc +/-) for every point in mSugra

// need to :
// convert to % wrt nominal

// return it in a data structure
  

  map<pair<int,int>, SignalEffData> allResults;

  for (  susyScanYields::iterator imsugra=rawYields.begin(); imsugra!=rawYields.end(); ++imsugra) {
    SignalEffData theseResults;

    if (rawYields[imsugra->first].first >0) {
      cout<<imsugra->first.first<<" "<<imsugra->first.second<<" " //m0 and m12
 	  <<rawYields[imsugra->first].first<<" "<<nominalYields[imsugra->first].first<<" "
 	  <<kFactorPlusYields[imsugra->first].first<<" "<<kFactorMinusYields[imsugra->first].first<<" ";

      theseResults.rawYield = rawYields[imsugra->first].first;

      double nominalYield = nominalYields[imsugra->first].first;
      theseResults.effCorr = nominalYield/theseResults.rawYield;
      if (nominalYield==0) nominalYield=0.000001; //should make this better

      //remove fabs()
      double kfactorP = ((kFactorPlusYields[imsugra->first].first - nominalYield) / nominalYield);
      double kfactorM = ((kFactorMinusYields[imsugra->first].first - nominalYield) / nominalYield);
      theseResults.set("kFactor",kfactorM,kfactorP);

      //pdf uncertainties
      //      if (fullPdfUncertainties) {
	
	double cteqXminusMax = cteqXplusminus.first[imsugra->first].first;
	double cteqXplusMax = cteqXplusminus.second[imsugra->first].first;

	//cout<<"\nCTEQ nominal xMinusMax xPlusMax "<<nominalYield<<" "<<cteqXminusMax<<" "<<cteqXplusMax<<endl;
	double  cteqFracPDF = 0.5* ( cteqXminusMax + cteqXplusMax);
	cteqFracPDF = fabs( cteqFracPDF ) /nominalYield;
	
	cout<<endl;
	cout<<imsugra->first.first<<" "<<imsugra->first.second<<" CTEQ  % = "<<100*cteqFracPDF<<endl;

	double mstwXminusMax = mstwXplusminus.first[imsugra->first].first;
	double mstwXplusMax = mstwXplusminus.second[imsugra->first].first;
	double  mstwFracPDF = 0.5* ( mstwXminusMax + mstwXplusMax);
	mstwFracPDF = fabs( mstwFracPDF ) /nominalYield;
	cout<<imsugra->first.first<<" "<<imsugra->first.second<<" MSTW  % = "<<100*mstwFracPDF<<endl;
	if (cteqFracPDF > mstwFracPDF) mstwFracPDF = cteqFracPDF;
	
	double nnpdfTotalX  = nnpdfXX2.first[imsugra->first].first;
	double nnpdfTotalX2 = nnpdfXX2.second[imsugra->first].first;
	double nnpdfFracPDF = sqrt(nnpdfTotalX2 - (nnpdfTotalX*nnpdfTotalX)) / nominalYield;
	cout<<imsugra->first.first<<" "<<imsugra->first.second<<" NNPDF % = "<<100*nnpdfFracPDF<<endl;
 	if (nnpdfFracPDF > mstwFracPDF) mstwFracPDF = nnpdfFracPDF;

	theseResults.set("PDF",mstwFracPDF,mstwFracPDF);

//       }
//       else {
//       //constant term for PDF uncertainties
//         theseResults.totalSystematic += 0.13*0.13; //PDF
//       }
	
      for ( map<TString, pair<susyScanYields,susyScanYields> >::iterator ivariation=variedYields.begin(); ivariation!=variedYields.end(); ++ivariation) {
	cout<<ivariation->second.first[imsugra->first].first<<" "<<ivariation->second.second[imsugra->first].first<<" ";
	
	//get rid of fabs()
	double percent1 = ((ivariation->second.first[imsugra->first].first - nominalYield)/nominalYield);
	double percent2 = ((ivariation->second.second[imsugra->first].first - nominalYield)/nominalYield);
	theseResults.set(ivariation->first,percent1,percent2);

      }
      cout<<endl;

      theseResults.setFixedForScan(); //set some uncertainties (PU, JER) that we fix only for the scans
     
    }
    //stuff theseResults into the map
    allResults[imsugra->first] = theseResults;
  }

  //so these are results for all msugra for one particular selection
  return allResults;

}



void runSystematics2011_scan(TString sampleOfInterest, const unsigned int i) {
  //now converting this to work also on mSugra "mSUGRAtanb40" as well as simplified models e.g. T1bbbb
  
  loadSamples();
  clearSamples();
  addSample(sampleOfInterest);
  
  setSearchRegions();

  if (sampleOfInterest.Contains("SUGRA"))  loadSusyScanHistograms();
  
  if (i>=searchRegions_.size() ) { cout<<"There are only "<<searchRegions_.size()<<" search regions!"<<endl; return;}
  
  //open the output files
  //-- these are the text, root files for Owen et al
  char effoutput[500];
  sprintf(effoutput,"signalSyst.%s.%s%s.dat",sampleOfInterest.Data(),searchRegions_[i].btagSelection.Data(),searchRegions_[i].owenId.Data());
  ofstream*  textfiles= new ofstream(effoutput);
  sprintf(effoutput,"RA2b.%s.%s%s.root",sampleOfInterest.Data(),searchRegions_[i].btagSelection.Data(),searchRegions_[i].owenId.Data());
  TFile*  rootfiles  =new TFile(effoutput,"RECREATE");

  //as specified by mariarosaria for T1bbbb
  // -- also works for T2bb and T2tt (see https://twiki.cern.ch/twiki/bin/viewauth/CMS?SUSY42XSUSYScan)
  int    nbinsx=60;
  int    nbinsy=60;
  double    lowx=0;
  double    lowy=0;
  double    highx=1500;
  double    highy=1500;
  if ( sampleOfInterest.Contains("mSUGRA")) {
    nbinsx=200;
    nbinsy=100;
    lowx=0;
    lowy=0;
    highx=2000;
    highy=1000;
  }

  //despite the name, runTH2Syst2011_mSugra also works for SMS
    map<pair<int,int>, SignalEffData> SB =  runTH2Syst2011_mSugra(sbRegions_[i],false,false,sampleOfInterest);
    map<pair<int,int>, SignalEffData> SIG =  runTH2Syst2011_mSugra(searchRegions_[i],false,false,sampleOfInterest);
    
    map<pair<int,int>, SignalEffData> SBSL =  runTH2Syst2011_mSugra(sbRegions_[i],true,false,sampleOfInterest);
    map<pair<int,int>, SignalEffData> SIGSL =  runTH2Syst2011_mSugra(searchRegions_[i],true,false,sampleOfInterest);
    
    map<pair<int,int>, SignalEffData> SBLDP =  runTH2Syst2011_mSugra(sbRegions_[i],false,true,sampleOfInterest);
    map<pair<int,int>, SignalEffData> SIGLDP =  runTH2Syst2011_mSugra(searchRegions_[i],false,true,sampleOfInterest);

    //big question: does this loop over SB work for mSugra? I think it should, but in the mSugra version I used scanProcessTotalsMap

    for (map<pair<int,int>, SignalEffData >::iterator iscanpoint = SB.begin(); iscanpoint!= SB.end(); ++iscanpoint) {
      int nentries=  sampleOfInterest.Contains("mSUGRA") ?
	scanProcessTotalsMap[iscanpoint->first]->GetEntries() :
	TMath::Nint(scanSMSngen->GetBinContent(scanSMSngen->FindBin(iscanpoint->first.first,iscanpoint->first.second ))); 
      
      (*textfiles)<<iscanpoint->first.first <<" "<<iscanpoint->first.second <<" "<<nentries<<" "
		     <<SIG[iscanpoint->first].rawYield<<" "<<SB[iscanpoint->first].rawYield<<" "<<SIGSL[iscanpoint->first].rawYield<<" "<<SBSL[iscanpoint->first].rawYield<<" "<<SIGLDP[iscanpoint->first].rawYield<<" "<<SBLDP[iscanpoint->first].rawYield<<" "
		     <<SIG[iscanpoint->first].effCorr<<" "<<SB[iscanpoint->first].effCorr<<" "<<SIGSL[iscanpoint->first].effCorr<<" "<<SBSL[iscanpoint->first].effCorr<<" "<<SIGLDP[iscanpoint->first].effCorr<<" "<<SBLDP[iscanpoint->first].effCorr<<" "
		     <<SIG[iscanpoint->first].totalSystematic()<<" "<<SB[iscanpoint->first].totalSystematic()<<" "<<SIGSL[iscanpoint->first].totalSystematic()<<" "<<SBSL[iscanpoint->first].totalSystematic()<<" "<<SIGLDP[iscanpoint->first].totalSystematic()<<" "<<SBLDP[iscanpoint->first].totalSystematic()<<endl;
    }

    //now output to root. need to convert back to TH2D
    rootfiles->cd();
    TString hsuffix = sampleOfInterest; hsuffix+="_"; hsuffix+=searchRegions_[i].btagSelection; hsuffix+=searchRegions_[i].owenId;
    TString hname=sampleOfInterest.Contains("mSUGRA") ? "correctedYield_": "efficiency_"; hname+=hsuffix;
    TH2D efficiency(hname,hname,nbinsx,lowx,highx,nbinsy,lowy,highy);
    hname="statError_"; hname+=hsuffix;
    TH2D statError(hname,hname,nbinsx,lowx,highx,nbinsy,lowy,highy);
    hname="btagError_"; hname+=hsuffix;
    TH2D btagError(hname,hname,nbinsx,lowx,highx,nbinsy,lowy,highy);
    hname="jesError_"; hname+=hsuffix;
    TH2D jesError(hname,hname,nbinsx,lowx,highx,nbinsy,lowy,highy);
    hname="metError_"; hname+=hsuffix;
    TH2D metError(hname,hname,nbinsx,lowx,highx,nbinsy,lowy,highy);
    hname="pdfError_"; hname+=hsuffix;
    TH2D pdfError(hname,hname,nbinsx,lowx,highx,nbinsy,lowy,highy);
    //    hname="otherError_"; hname+=hsuffix;
    //    TH2D otherError(hname,hname,nbinsx,lowx,highx,nbinsy,lowy,highy);

    hname="totalSystematicError_"; hname+=hsuffix;
    TH2D totalSystematicError(hname,hname,nbinsx,lowx,highx,nbinsy,lowy,highy);
    //only SIG for now
    for (map<pair<int,int>, SignalEffData>::iterator iscan=SIG.begin(); iscan!=SIG.end(); ++iscan) {
      double N = iscan->second.rawYield * iscan->second.effCorr ;

      double effval=0;
      double staterrval=0;
      if (sampleOfInterest.Contains("mSUGRA") ) {
	effval = N;
      }
      else {
	double ngen=  TMath::Nint(scanSMSngen->GetBinContent(scanSMSngen->FindBin(iscan->first.first,iscan->first.second ))); 
	effval = N/ngen;
	staterrval = sqrt( N*(1-N/ngen));
      }

      efficiency.Fill( iscan->first.first, iscan->first.second, effval );
      statError.Fill(  iscan->first.first, iscan->first.second, staterrval );
      btagError.Fill(  iscan->first.first, iscan->first.second, iscan->second.value("btag") );
      jesError.Fill(  iscan->first.first, iscan->first.second, iscan->second.value("JES") );
      metError.Fill(  iscan->first.first, iscan->first.second, iscan->second.value("MET") );
      pdfError.Fill(  iscan->first.first, iscan->first.second, iscan->second.value("PDF") );
      //      otherError.Fill(  iscan->first.first, iscan->first.second,  );
      totalSystematicError.Fill(iscan->first.first, iscan->first.second, iscan->second.totalSystematic() );
    }
    rootfiles->Write();
    
    textfiles->close();
    rootfiles->Close();
    //}


}

/* hope to deprecate this function
//run mSugra systematics in the efficient way
void runSystematics2011_mSugra(const unsigned int i) {
  //add index i argument to define the search region we're interested in.

  TString sampleOfInterest="mSUGRAtanb40";
  loadSamples();
  clearSamples();
  addSample(sampleOfInterest);

  setSearchRegions();
  loadSusyScanHistograms();
  //  runTH2Syst2011_mSugra(searchRegions_[0]);
  //    return;

  if (i>=searchRegions_.size() ) { cout<<"There are only "<<searchRegions_.size()<<" search regions!"<<endl; return;}


  //open the output files
  vector<ofstream*> textfiles;
  //  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    char effoutput[500];
    sprintf(effoutput,"signalSyst.%s.%s%s.dat",sampleOfInterest.Data(),searchRegions_[i].btagSelection.Data(),searchRegions_[i].owenId.Data());
    textfiles.push_back( new ofstream(effoutput));
    //  }

    //  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    
    map<pair<int,int>, SignalEffData> SB =  runTH2Syst2011_mSugra(sbRegions_[i],false,false,sampleOfInterest);
    map<pair<int,int>, SignalEffData> SIG =  runTH2Syst2011_mSugra(searchRegions_[i],false,false,sampleOfInterest);
    
    map<pair<int,int>, SignalEffData> SBSL =  runTH2Syst2011_mSugra(sbRegions_[i],true,false,sampleOfInterest);
    map<pair<int,int>, SignalEffData> SIGSL =  runTH2Syst2011_mSugra(searchRegions_[i],true,false,sampleOfInterest);
    
    map<pair<int,int>, SignalEffData> SBLDP =  runTH2Syst2011_mSugra(sbRegions_[i],false,true,sampleOfInterest);
    map<pair<int,int>, SignalEffData> SIGLDP =  runTH2Syst2011_mSugra(searchRegions_[i],false,true,sampleOfInterest);
     
  
    for (map<pair<int,int>, TH1D* >::iterator iscanpoint = scanProcessTotalsMap.begin(); iscanpoint!=scanProcessTotalsMap.end(); ++iscanpoint) {
      int nentries=    scanProcessTotalsMap[iscanpoint->first]->GetEntries();
      cout<<iscanpoint->first.first<<" "<<iscanpoint->first.second<<" NEntries = "<<nentries<<endl;
      (*textfiles[0])<<iscanpoint->first.first <<" "<<iscanpoint->first.second <<" "<<nentries<<" "
		     <<SIG[iscanpoint->first].rawYield<<" "<<SB[iscanpoint->first].rawYield<<" "<<SIGSL[iscanpoint->first].rawYield<<" "<<SBSL[iscanpoint->first].rawYield<<" "<<SIGLDP[iscanpoint->first].rawYield<<" "<<SBLDP[iscanpoint->first].rawYield<<" "
		     <<SIG[iscanpoint->first].effCorr<<" "<<SB[iscanpoint->first].effCorr<<" "<<SIGSL[iscanpoint->first].effCorr<<" "<<SBSL[iscanpoint->first].effCorr<<" "<<SIGLDP[iscanpoint->first].effCorr<<" "<<SBLDP[iscanpoint->first].effCorr<<" "
		     <<SIG[iscanpoint->first].totalSystematic()<<" "<<SB[iscanpoint->first].totalSystematic()<<" "<<SIGSL[iscanpoint->first].totalSystematic()<<" "<<SBSL[iscanpoint->first].totalSystematic()<<" "<<SIGLDP[iscanpoint->first].totalSystematic()<<" "<<SBLDP[iscanpoint->first].totalSystematic()<<endl;
    }
    //  }

    //  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    textfiles.at(0)->close();
    //  }

}
*/

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

const bool reweightLSBdata_=false; //whether or not LSB data is reweighted based on PV distribution
const bool useScaleFactors_=true; //whether or not to use MC scale factors when doing subtraction for data-driven estimates

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
    sampleOfInterest = "data";
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
    triggerCutLSB = "(pass_utilityHLT==1 || !isRealData)"; //because bug in MC for pass_utilityHLT
    triggerCut = "cutTrigger==1";
  }

  TCut ge1b = "nbjetsCSVM>=1";
  btagSFweight_="1";
  TString btagSFweight=""; //we're going to have to switch this one in and out of the global var
  if (useScaleFactors_) {
    ge1b="1";
    usePUweight_=true;
    useHLTeff_=true;
    currentConfig_=configDescriptions_.getCorrected(); //add JERbias

    if (btagselection=="ge2b") {
      btagSFweight="probge2";
    }
    else if (btagselection=="ge1b") {
      btagSFweight="probge1";
    }
    else if (btagselection=="ge3b") {
      btagSFweight="probge3";
    }
    else {assert(0);}
    btagSFweight_=btagSFweight;
  }
  else {
    usePUweight_=false;
    useHLTeff_=false;
    btagSFweight_="1";
    currentConfig_=configDescriptions_.getDefault(); //completely raw MC

    if (btagselection=="ge2b") {
      ge1b="nbjetsCSVM>=2";
    }
    else if (btagselection=="ge1b") {}
    else if (btagselection=="ge3b") {
      ge1b="nbjetsCSVM>=3";
    }
    else {assert(0);}
  }

  TCut HTcut=region.htSelection.Data(); 
  TCut SRMET = region.metSelection.Data();

  SRMET = SRMET && ge1b && triggerCut; //apply physics trigger and physics b cut in high MET region

  TCut baseline = "cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1";
  baseline = baseline&&HTcut;
  TCut cleaning = "weight<1000";

  char cutstring1[100];
  sprintf(cutstring1,"MET>= %.0f && MET < %.0f && nbjetsCSVM%s", 50+SBshift,100+SBshift,LSBbsel.Data());
  //  cout<<"*** SB cut is "<<cutstring1<<endl<<endl;
  TCut SBMET = TCut(cutstring1)&&triggerCutLSB;
  TCut dpcut = "1";//"minDeltaPhiN>=4";
  //  TCut passOther = "deltaPhiMPTcaloMET<2";
  //  TCut failOther = "deltaPhiMPTcaloMET>=2";
  TCut passOther = "minDeltaPhiN>=4";
  TCut failOther = "minDeltaPhiN<4";

  var="HT"; xtitle=var;
  nbins=10; low=0; high=5000;
  double A,B,D,SIG,Aerr,Berr,Derr,SIGerr;
  if(datamode &&  reweightLSBdata_){
    assert(0);//error is wrong
    int pvnbins=18;
    float pvbins[]={0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5};
    
    //physics triggers control sample -- physics triggers, 0b
    //--consider varying MET and HT cuts and PV binning!!
    selection_ = TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjetsCSVM==0") && triggerCut;
    drawSimple("nGoodPV",pvnbins,pvbins,"dummy", "","data");
    TH1D* hPVphysics = (TH1D*)hinteractive->Clone("hPVphysics");
    for(int j=1; j<=hPVphysics->GetNbinsX();j++){
      assert(hPVphysics->GetBinContent(j)>0);
    }
    
    //LSB unweighted
    selection_ = baseline && SBMET;
    drawSimple("nGoodPV",pvnbins,pvbins,"dummy", "","data");
    TH1D* hPVprescale = (TH1D*)hinteractive->Clone("hPVprescale");

    //calculate weights (physics shape with prescale integral divided by prescale)
    TH1D* hPVprescale_RW = (TH1D*)hPVphysics->Clone("hPVprescale_RW");
    //Now essentially do "hPVprescale_RW->Scale(hPVprescale->Integral()/hPVphysics->Integral());" but with correct error (small usually)
    double s = hPVprescale->Integral()/hPVphysics->Integral();
    double serr = jmt::errAoverB(hPVprescale->Integral(),jmt::errOnIntegral(hPVprescale),hPVphysics->Integral(),jmt::errOnIntegral(hPVphysics));
    for(int j=1; j<=hPVprescale_RW->GetNbinsX(); j++){
      assert(hPVprescale->GetBinContent(j)>0);
      hPVprescale_RW->SetBinContent(j,s*hPVprescale_RW->GetBinContent(j));
      double fillerror = sqrt(serr*serr+hPVprescale_RW->GetBinError(j)*hPVprescale_RW->GetBinError(j));
      hPVprescale_RW->SetBinError(j,fillerror);
    }
    TH1D* hPV_W = (TH1D*)hPVphysics->Clone("hPV_W");
    hPV_W->Reset();
    hPV_W->Divide(hPVprescale_RW,hPVprescale); 
    
    //LSB pass mdp unweighted
    selection_ = baseline && SBMET && passOther;
    drawSimple("nGoodPV",pvnbins,pvbins,"dummy", "","data");
    TH1D* hPVprescalePass = (TH1D*)hinteractive->Clone("hPVprescalePass");
    
    //LSB fail mdp unweighted
    selection_ = baseline && SBMET && failOther;
    drawSimple("nGoodPV",pvnbins,pvbins,"dummy", "","data");
    TH1D* hPVprescaleFail = (TH1D*)hinteractive->Clone("hPVprescaleFail");

    //weighted LSB pass and fail mdp
    TH1D* hPVprescalePass_RW = (TH1D*)hPVprescalePass->Clone("hPVprescalePass_RW");
    TH1D* hPVprescaleFail_RW = (TH1D*)hPVprescaleFail->Clone("hPVprescaleFail_RW");
    hPVprescalePass_RW->Multiply(hPV_W);
    hPVprescaleFail_RW->Multiply(hPV_W);

    A = hPVprescaleFail_RW->Integral();
    Aerr = jmt::errOnIntegral(hPVprescaleFail_RW);    
    B = hPVprescalePass_RW->Integral();
    Berr = jmt::errOnIntegral(hPVprescalePass_RW);
  }
  else{
    //A   -- aka 50 - 100 and high MPT,MET
    selection_ = baseline && cleaning && dpcut  && SBMET && failOther; //auto cast to TString seems to work
    drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_A");
    A=getIntegral(sampleOfInterest);
    Aerr=getIntegralErr(sampleOfInterest);
    //B
    selection_ = baseline && cleaning && dpcut  && SBMET && passOther; //auto cast to TString seems to work
    drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_B");
    B=getIntegral(sampleOfInterest);
    Berr=getIntegralErr(sampleOfInterest);
  }
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
  //  double closureStat= datamode? 0: jmt::errAoverB(estimate,estimateerr,SIG,SIGerr); //comment out unused var

  double closureStat2 = datamode? 0: jmt::errAoverB(SIG,SIGerr,estimate,estimateerr);

  //for a cross-check
  double R0 = B/A;
  double R0err = jmt::errAoverB(B,Berr,A,Aerr);


//   cout<<" ==== "<<endl
//       <<"Estimate = "<<estimate<<" +/- "<<estimateerr<<endl
//       <<" truth   = "<<SIG     <<" +/- "<<SIGerr<<endl;
  TString name = region.btagSelection;
  name += region.owenId;
  name += isSIG ? ", SIG":", SB";
  char output[500];
  if (!datamode) {
    //sprintf(output,"%s & %s & %s & %s & %s & %s \\\\ %% %f ++ %f %% alternate denominator: %f ++ %f",btagselection.Data(),
    //    jmt::format_nevents(B,Berr).Data(),jmt::format_nevents(A,Aerr).Data(),
    //    jmt::format_nevents(D,Derr).Data(),jmt::format_nevents(estimate,estimateerr).Data(),
    //    jmt::format_nevents(SIG,SIGerr).Data(),100*(SIG-estimate)/SIG,closureStat*100,100*(SIG-estimate)/estimate,closureStat2*100);
      
    sprintf(output,"%s & %s & %s & %s & %s & %s & $%f \\pm %f$ \\\\ ",name.Data(),
	    jmt::format_nevents(B,Berr).Data(),jmt::format_nevents(A,Aerr).Data(),
	    jmt::format_nevents(D,Derr).Data(),jmt::format_nevents(estimate,estimateerr).Data(),
	    jmt::format_nevents(SIG,SIGerr).Data(),100*(estimate-SIG)/estimate,closureStat2*100);
 
  }
  else {
    sprintf(output,"%s & %d & %d & %d & %s & %s   \\\\ %% %f +/- %f",name.Data(),
	    TMath::Nint(B),TMath::Nint(A),
	    TMath::Nint(D), jmt::format_nevents(Dsub,Dsuberr).Data(),
	    jmt::format_nevents(estimate,estimateerr).Data(),R0,R0err);
    cout<<"(qcd) DATA\t";
  }
  cout<<output<<endl;

  if (datamode)  return make_pair(estimate,estimateerr);
  //return make_pair( 100*sqrt(pow((SIG-estimate)/SIG,2) + pow(closureStat,2)), 0); //MC truth in denominator
  return make_pair( 100*sqrt(pow((SIG-estimate)/estimate,2) + pow(closureStat2,2)), 0); //estimate in denominator
}


void runClosureTest2011(std::map<TString, std::vector<double> > & syst, bool addReweightedQCDClosureTest = true)  {
  if(!addReweightedQCDClosureTest) cout << "You are starting closure test without considering njet reweighting!" << endl;
  
  setSearchRegions();
  assert(sbRegions_.size() == searchRegions_.size());
  
  string selString;
  double unweightedClosure, weightedClosure;
  double unweightedClosureA_sb[]={0,0,0,0};
  double weightedClosureA_sb[]={0,0,0,0};
  double unweightedClosureA_sig[]={0,0,0,0};
  double weightedClosureA_sig[]={0,0,0,0};
  
  if(addReweightedQCDClosureTest){
    ifstream inFile("weightedClosure.dat", std::ios::in);//file generated by additionalTools/reweightQCD.C (root -l -b -q reweightQCD.C++)
    if(!inFile.good()){
      cout << "could not open weightedClosure.dat" << endl;
      assert(0);
    }
    
    while(inFile>>selString>>unweightedClosure>>weightedClosure){
      if(selString=="ge1bloosesb")  {unweightedClosureA_sb[0]=unweightedClosure; weightedClosureA_sb[0]=weightedClosure;}
      if(selString=="ge1btightsb")  {unweightedClosureA_sb[1]=unweightedClosure; weightedClosureA_sb[1]=weightedClosure;}
      if(selString=="ge2bloosesb")  {unweightedClosureA_sb[2]=unweightedClosure; weightedClosureA_sb[2]=weightedClosure;}
      if(selString=="ge2btightsb")  {unweightedClosureA_sb[3]=unweightedClosure; weightedClosureA_sb[3]=weightedClosure;}
      if(selString=="ge1bloosesig") {unweightedClosureA_sig[0]=unweightedClosure; weightedClosureA_sig[0]=weightedClosure;}
      if(selString=="ge1btightsig") {unweightedClosureA_sig[1]=unweightedClosure; weightedClosureA_sig[1]=weightedClosure;}
      if(selString=="ge2bloosesig") {unweightedClosureA_sig[2]=unweightedClosure; weightedClosureA_sig[2]=weightedClosure;}
      if(selString=="ge2btightsig") {unweightedClosureA_sig[3]=unweightedClosure; weightedClosureA_sig[3]=weightedClosure;}
    }
  }
  
  for (unsigned int i=0; i<sbRegions_.size(); i++) {
    double jsb =   fabs(anotherABCD(sbRegions_[i]).first);
    double jsig =  fabs(anotherABCD(searchRegions_[i]).first);
    
    double finalsb = jsb;
    double finalsig = jsig;
    if(addReweightedQCDClosureTest){
      double bsb =  fabs(unweightedClosureA_sb[i]);
      double bsig = fabs(unweightedClosureA_sig[i]);
      double wsb =  fabs(weightedClosureA_sb[i]);
      double wsig = fabs(weightedClosureA_sig[i]);
      cout << "SB: " << jsb << " " << bsb << endl;
      cout << "SIG: " << jsig << " " << bsig << endl;
      assert( fabs(jsb-bsb)   <0.001);//should be true if using same ntuples for reweighted study and current running. could take out.
      assert( fabs(jsig-bsig) <0.001);
      finalsb = jsb>wsb ? jsb:wsb ;
      finalsig = jsig>wsig ? jsig:wsig ;
    }
    syst["Closure"].push_back( finalsb );
    syst["Closure"].push_back( finalsig ); 
  }
  
}

void runClosureTest2011()  {
  std::map<TString, std::vector<double> > dummy;
  
  cout << "You are starting closure test without considering njet reweighting!" << endl;
  runClosureTest2011(dummy,false);

}


//to hold systematics
//making this global because I'm lazy...
//the string is the category of systematic
//the vector is list of values for all of the search and sb regions (not a very good approach)
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
  runClosureTest2011(qcdSystErrors,false);
  
  cout<<" == QCD systematics summary =="<<endl;
  for (unsigned int j=0; j<n.size(); j++) {
    qcdSystErrors["Total"].push_back( sqrt( pow(qcdSystErrors["MCsub"].at(j),2) +  pow(qcdSystErrors["Closure"].at(j),2)+ pow(qcdSystErrors["SBshift"].at(j),2)));
    cout<<j<<"\t&"<<qcdSystErrors["MCsub"].at(j)<<" & "<<qcdSystErrors["Closure"].at(j)<<" & "<<qcdSystErrors["SBshift"].at(j)<<" & "<<qcdSystErrors["Total"].at(j)<<endl;
  }
  

}

//i don't like passing the index instead of the region itself, but it makes some things easier.
//this code is all around a big kludge...
double slABCD(const unsigned int searchRegionIndex, bool datamode=false, const TString & mode="", const bool justttbar=false ) {
  //in datamode, return the estimate; in non-datamode, return the Closure Test results (true-pred)/pred

  if (justttbar == true) cout<<"Will run closure test in ttbar only mode!"<<endl;

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
      SBsubQCDerr = sqrt( SBsubQCDerr*SBsubQCDerr + pow(qcdsbsyst*SBsubQCD,2));
      
      if (mode=="QCDup") SBsubQCD += SBsubQCDerr;
      if (mode=="QCDdown") SBsubQCD -= SBsubQCDerr;
    }
    
    //now hard-coded Z->nunu
    double zv[2];
    double ze[2];
    //sometimes need to average mu mu and ee estimates
    double zsbsyst=0.5;
    double zinvscale = 3.5/3.2;
    bool doMean = true;

    if (qcdsubregion.owenId == "Loose" && qcdsubregion.btagSelection=="ge1b") {
      doMean=false;//averaging already done
      zv[0] = 63.9127; ze[0]=10.7967;
      zsbsyst = 0.3;
    }
    else if (qcdsubregion.owenId == "Tight" && qcdsubregion.btagSelection=="ge1b") {
      doMean=false;//averaging already done
      zv[0] = 35.3791; ze[0]=7.95003;
      zsbsyst = 0.3;
    }
    else if (qcdsubregion.owenId == "case3" ) {
      doMean=false;//averaging already done
      zv[0] = 4.91838; ze[0]=2.92228;
      zsbsyst = 0.3;
    }
    else if (qcdsubregion.owenId == "Loose" && qcdsubregion.btagSelection=="ge2b") {
      doMean=false;//averaging already done
      zv[0] = 10.7258; ze[0]=2.18597;
      zsbsyst = 0.5;
    }
    else if (qcdsubregion.owenId == "Tight" && qcdsubregion.btagSelection=="ge2b") {
      doMean=false;//averaging already done
      zv[0] = 3.41633; ze[0]=1.08638;
      zsbsyst = 0.5;
    }
    else if (qcdsubregion.owenId == "Loose" && qcdsubregion.btagSelection=="ge3b") {
      doMean=false;//averaging already done
      zv[0] = 1.48325; ze[0]=0.527006;
      zsbsyst = 0.7;   
    }
    else {assert(0);}
    if(doMean){
      SBsubZ = zinvscale*(jmt::weightedMean(2,zv,ze));
      SBsubZerr = zinvscale*(jmt::weightedMean(2,zv,ze,true));
    }
    else{//if weighted averaging was already done
      SBsubZ = zinvscale*zv[0];
      SBsubZerr = zinvscale*ze[0];
    }

    if (mode.Contains("Z")) {
      SBsubZerr = sqrt( SBsubZerr*SBsubZerr + pow(zsbsyst*SBsubZ,2));
      
      if (mode=="Zup") SBsubZ += SBsubZerr;
      if (mode=="Zdown") SBsubZ -= SBsubZerr;
    }
  }

  setStackMode(false);
  doData(true);
  setQuiet(true);

  TString sampleOfInterest = "totalsm";
  useFlavorHistoryWeights_=false;
  loadSamples(); //needs to come first! this is why this should be turned into a class, with this guy in the ctor
  setColorScheme("nostack");
  clearSamples();
  if (datamode) {
    addSample("ZJets");
    addSample("VV");
    sampleOfInterest="data";
  }
  else {
    addSample("TTbarJets");
    if (!justttbar) {
      addSample("WJets");
      addSample("SingleTop");
    }
  }

  savePlots_=false;

  setLogY(false);
  TString var,xtitle;
  int nbins;
  float low,high;
  
  // --  count events
  TCut ge1b = "nbjetsCSVM>=1";
  btagSFweight_="1";
  TString btagSFweight=""; //we're going to have to switch this one in and out of the global var
  if (useScaleFactors_) {
    ge1b="1";
    usePUweight_=true;
    useHLTeff_=true;
    currentConfig_=configDescriptions_.getCorrected(); //add JERbias
    
    if (btagselection=="ge2b") {
      btagSFweight="probge2";
    }
    else if (btagselection=="ge1b") {
      btagSFweight="probge1";
    }
    else if (btagselection=="ge3b") {
      btagSFweight="probge3";
    }
    else {assert(0);}
    btagSFweight_=btagSFweight;
  }
  else {
    usePUweight_=false;
    useHLTeff_=false;
    btagSFweight_="1";
    currentConfig_=configDescriptions_.getDefault(); //completely raw MC
    
    if (btagselection=="ge2b") {
      ge1b="nbjetsCSVM>=2";
    }
    else if (btagselection=="ge1b") {}
    else if (btagselection=="ge3b") {
      ge1b="nbjetsCSVM>=3";
    }
    else {assert(0);}
  }

  TCut HTcut=region.htSelection.Data(); 
  TCut SRMET = region.metSelection.Data();
  TCut baseline = "cutPV==1 && cut3Jets==1 && cutTrigger==1";
  baseline = baseline&&HTcut&&ge1b;
  TCut cleaning = "weight<1000";
  TCut SBMET = qcdsubregion.metSelection.Data();//"MET>=150 && MET<200";
  TCut dpcut = "minDeltaPhiN>=4";
  TCut failOther = "(((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0&&MT_Wlep<100)";
  TCut passOther = "nElectrons==0 && nMuons==0";

  double A,B,D,SIG,Aerr,Berr,Derr,SIGerr;
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
  double closureStat= datamode? 0: jmt::errAoverB(SIG,SIGerr,estimate,estimateerr);


//   cout<<" ==== "<<endl
//       <<"Estimate = "<<estimate<<" +/- "<<estimateerr<<endl
//       <<" truth   = "<<SIG     <<" +/- "<<SIGerr<<endl;
  // btagselection += tight ? " Tight " : " Loose ";
  TString name = btagselection;
  name += region.owenId;

  char output[500];
  if (!datamode) {
    //sprintf(output,"%s & %s & %s & %s & %s & %s \\\\ %% %f ++ %f",btagselection.Data(),
    //	    jmt::format_nevents(D,Derr).Data(),jmt::format_nevents(A,Aerr).Data(),
    //	    jmt::format_nevents(B,Berr).Data(),jmt::format_nevents(estimate,estimateerr).Data(),
    //	    jmt::format_nevents(SIG,SIGerr).Data(),100*(estimate-SIG)/estimate,100*closureStat);

    sprintf(output,"%s & %s & %s & %s & %s & %s & $%f \\pm %f$ \\\\",name.Data(),
	    jmt::format_nevents(D,Derr).Data(),jmt::format_nevents(A,Aerr).Data(),
	    jmt::format_nevents(B,Berr).Data(),jmt::format_nevents(estimate,estimateerr).Data(),
	    jmt::format_nevents(SIG,SIGerr).Data(), 100*(estimate-SIG)/estimate, 100*closureStat);
  }
  else {
    sprintf(output,"ttbar DATA %s & %d & %d & %d & %s & %s  \\\\",name.Data(),
	    TMath::Nint(D),TMath::Nint(A),
	    TMath::Nint(B), jmt::format_nevents(SBsubMisc+SBsubQCD+SBsubZ,suberr).Data(),
	    jmt::format_nevents(estimate,estimateerr).Data());
  }
  cout<<output<<endl;

  if (datamode) return estimate;
  return 100*sqrt(pow((SIG-estimate)/estimate,2) + pow(closureStat,2));
}

void runSLClosureTest2011() {

  setSearchRegions();

  cout<<"Note that the following is the joint tt+W+t closure test only!"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) slABCD(j);

  cout<<"Note that the following is the tt closure test only!"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) slABCD(j,false,"",true);

}

vector<double> ttbarClosureSyst;
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
  cout<<"Running closure tests"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) {
    double allsamples=fabs(slABCD(j)); //mix of samples
    double ttbaronly=fabs(slABCD(j,false,"",true)); //just ttbar
    closure[j] = (allsamples>ttbaronly) ? allsamples : ttbaronly;
  }

  ttbarClosureSyst.clear();
  cout<<" == summary (%) == "<<endl;
  cout<<"\tClosure\tQCD\tZ\tMC\tTotal"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) {
    double totalsyst = sqrt(closure[j]*closure[j] + qcd[j]*qcd[j] + znn[j]*znn[j] + mc[j]*mc[j]);
    cout<<j<<"\t"<<closure[j]<<"\t"<<qcd[j]<<"\t"<<znn[j]<<"\t"<<mc[j]<<"\t"<<totalsyst <<endl;
    ttbarClosureSyst.push_back(closure[j]);
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

void printOwenSyst() {

  //note that because of the way the code is written, i am not at all sure that
  //it is safe to run printOwenAll() and printOwenSyst() in the same session. better to be safe!

  runDataQCD2011();
  runTtbarEstimate2011();

  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    char effoutput[500];
    sprintf(effoutput,"backgroundSyst.%s%s.%dinvpb.dat",searchRegions_[i].btagSelection.Data(),searchRegions_[i].owenId.Data(),TMath::Nint(lumiScale_));

    ofstream textfile(effoutput);

    textfile<<"sf_mc            "<<1<<endl;
    textfile<<"sf_mc_err        "<<0.5<<endl;
    textfile<<"sf_qcd_sb        "<<1<<endl;
    textfile<<"sf_qcd_sb_err    "<<0.01*sqrt(pow(qcdSystErrors["Closure"].at(i),2)+pow(qcdSystErrors["SBshift"].at(i),2))<<endl;
    textfile<<"sf_qcd_sig       "<<1<<endl;
    textfile<<"sf_qcd_sig_err   "<<0.01*sqrt(pow(qcdSystErrors["Closure"].at(2*i+1),2)+pow(qcdSystErrors["SBshift"].at(2*i+1),2))<<endl;
    textfile<<"sf_ttwj_sig      "<<1<<endl;
    textfile<<"sf_ttwj_sig_err  "<<0.01*ttbarClosureSyst[i]<<endl;
    textfile.close();
  }

}

void AN2011_prescale( TString btagselection="ge1b",const int mode=1 ) {
  /*
.L drawReducedTrees.C++
  */
  loadSamples();

  //this mode thing is kludgey
  TString modestring="";
  if (mode==1) {
    usePUweight_=false;
    useHLTeff_=false;
    btagSFweight_="1";
    currentConfig_=configDescriptions_.getDefault(); //completely raw MC
  }
  else if (mode==2 || mode==3) {
    assert(0); //not ready
    usePUweight_=true;
    useHLTeff_=true;
    currentConfig_=configDescriptions_.getCorrected(); //JER bias
    if (mode==2) modestring="-JER-PU-HLT";
    else if(mode==3) modestring="-JER-PU-HLT-bSF";
  }
  else assert(0);

  TCut btagcut = "nbjetsCSVM>=1";
  if (mode==1 || mode==2) {
    if ( btagselection=="ge1b") {} //do nothing
    else  if ( btagselection=="ge2b" ) {
      btagcut = "nbjetsCSVM>=2";
    }
    else if ( btagselection=="eq1b" ) {
      btagcut = "nbjetsCSVM==1";
    }
    else if ( btagselection=="eq0b" ) {
      btagcut = "nbjetsCSVM==0";
    }
    else if (btagselection=="ge3b" ){
      btagcut = "nbjetsCSVM>=3";
    }
    else {
      assert(0);
    }
  }
  else if (mode==3) {
    if ( btagselection=="ge1b") {
      btagcut="1";
      btagSFweight_="probge1";
    }
    else  if ( btagselection=="ge2b" ) {
      btagcut = "1";
      btagSFweight_="probge2";
    }
    else if ( btagselection=="eq0b" ) {
      btagcut = "1";
      btagSFweight_="prob0";
    }
    else {
      assert(0);
    }
  }



  int nbins;
  float low,high;
  TString var,xtitle;
  
  doOverflowAddition(true);

  setQuiet(false);

  // ==========================

  TCut util = "(pass_utilityHLT==1 || !isRealData) && weight<1000";
  //lumiScale_ = 19.23288; //customize lumiScale_
  lumiScale_ = 30.15471; 

  // ========= regular N-1 plots

  resetSamples(); //use all samples
  setStackMode(true); //regular stack
  setColorScheme("stack");
  doData(true);
  drawMCErrors_=true;

  /*
  //MET
  selection_ =TCut("cutHT==1 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&util&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 20; low=0; high=200;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_MET_"+btagselection+modestring);
  setLogY(true);  setPlotMinimum(1e-1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_MET_"+btagselection+modestring);

  //njets in 50 < MET < 100 region
  selection_ =TCut("cutHT==1 && cutPV==1 && MET>=50 && MET<100 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&util&&btagcut;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 8; low=1; high=9;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_njets_LSB_"+btagselection+modestring);

  selection_ =TCut("cutHT==1 && cutPV==1 && MET>=50 && MET<100 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && weight<10000")&&util&&btagcut;
  var="minDeltaPhiN"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 20; low=0; high=40;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_minDeltaPhiN_LSB_"+btagselection+modestring);
  setLogY(true);  setPlotMinimum(1e-1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_minDeltaPhiN_LSB_"+btagselection+modestring);


  //look again at MET
  selection_ =TCut("MET>=50 && MET<100 && cutHT==1 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&util&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 10; low=50; high=100;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_METlimited_"+btagselection+modestring);


  //sanity check
  btagSFweight_="1";
  doRatioPlot(true);
  ratioMin=0; ratioMax=3;
  selection_ =TCut("MET>150 && cutHT==1 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN < 4")&&util;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 10; low=150; high=250;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_METhigh_lowDP");

  lumiScale_ =  1143;
  ratioMin=1; ratioMax=1.5;
  selection_ ="cutTrigger==1 && MET>150 && cutHT==1 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN < 4";
  drawPlots(var,nbins,low,high,xtitle,"Events", "unprescaled_METhigh_lowDP");
  */

 
  //// ben
  
  doRatioPlot(true);
  selection_ =TCut("HT>=400 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 10; low=0; high=100;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_MET_"+btagselection+modestring);
  
  
  doRatioPlot(true);
  selection_ =TCut("HT>=400 && cutPV==1 && MET>=50 && MET<100 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&btagcut;
  var="minDeltaPhiN"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 20; low=0; high=40;
  setLogY(false); resetPlotMinimum();
  //drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_minDeltaPhiN_LSB_"+btagselection+modestring);
  setLogY(true);  setPlotMinimum(1e-1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_minDeltaPhiN_LSB_"+btagselection+modestring);
  



 
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

  TCut btagcut = "nbjetsCSVM>=1";
  if ( btagselection=="ge1b") {} //do nothing
  else  if ( btagselection=="ge2b" ) {
    btagcut = "nbjetsCSVM>=2";
  }
  else if ( btagselection=="eq1b" ) {
    btagcut = "nbjetsCSVM==1";
  }
  else if ( btagselection=="ge0b" ) {
    btagcut = "nbjetsCSVM>=0";
  }
  else if ( btagselection=="ge3b" ) {
    btagcut = "nbjetsCSVM>=3";
  }
  else {
    assert(0);
  }


  TCut HTcut="HT>=400";
  if (btagselection=="ge1b" && HTselection=="Tight")  HTcut="HT>=500";
  if (btagselection=="ge2b" && HTselection=="Tight")  HTcut="HT>=600";

  loadSamples();

  int nbins;
  float low,high;
  TString var,xtitle;
  bool vb = true;
  
  //ge1b, Loose
  //const int nvarbins=18;
  //const float varbins[]={100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 550}; //AN and PAS

  //ge1b, Tight
  //const int nvarbins=15;
  //const float varbins[]={100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 350, 400, 450, 500, 550}; 
  
  //ge2b, Loose and Tight
  const int nvarbins=13;
  const float varbins[]={100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550}; //AN and PAS 
  //const int nvarbins=12;
  //const float varbins[]={200, 225, 250, 275, 300, 350, 400, 450, 550,650,750, 850,950}; 

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
  else if(samplename=="SingleTop") sample = "singletop";

  selection_ =TCut("cutPV==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100")&&btagcut&&HTcut;
  //selection_ =TCut("MET>=150 && cutPV==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4")&&btagcut&&HTcut; //AN v5
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 40; low=100; high=550; 
  //nbins = 40; low=200; high=600; 
  //nbins = 40; low=150; high=550; //AN v5
  //nbins = 20; low=150; high=550;
  //nbins = 16; low=150; high=550;
  if(vb) drawPlots(var,nvarbins, varbins, xtitle,"Arbitrary units", "SBandSIG_MET_SL_"+sample+"_"+btagselection+"_"+HTselection);
  else drawPlots(var,nbins,low,high,xtitle,"Arbitrary units", "SBandSIG_MET_SL_"+sample+"_"+btagselection+"_"+HTselection);

  int lowbin = hinteractive->FindBin(low);
  int boundarybin_sb = 0, boundarybin_sig=0;

  boundarybin_sb = hinteractive->FindBin(250);

  //MET Signal region depends on "Loose" or "Tight" selection
  if (HTselection=="Tight" && btagselection=="ge1b")  boundarybin_sig = hinteractive->FindBin(500);
  else if (HTselection=="Tight" && btagselection=="ge2b")  boundarybin_sig = hinteractive->FindBin(300);
  else boundarybin_sig = hinteractive->FindBin(250);

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
  selection_ =TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&btagcut&&HTcut;
  //selection_ =TCut("MET>=150 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&btagcut&&HTcut;// AN v5
  if(vb)drawPlots(var,nvarbins, varbins ,xtitle,"Arbitrary units", "SBandSIG_MET_normal_"+sample+"_"+btagselection+"_"+HTselection);
  else drawPlots(var,nbins,low,high,xtitle,"Arbitrary units", "SBandSIG_MET_normal_"+sample+"_"+btagselection+"_"+HTselection);
  TH1D* SIGplot = (TH1D*)hinteractive->Clone("SIGplot");
  SLplot->SetLineColor(kRed);
  SLplot->SetMarkerColor(kRed);
  SLplot->Draw("SAME");
  
  thecanvas->SaveAs("METshape_"+sample+"_SLandStandard_"+btagselection+"_"+HTselection+".pdf");
  //not enough stats in WJets to really show anything

  //Hack to get it plotted with ratio plot
   //Hack to get it plotted with ratio plot
  TCanvas* myC = 0;
  myC = new TCanvas("myC", "myC", 600,700);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  Float_t small = 1e-5;
  myC->Divide(1,2,small,small);
  //myC->Divide(1,2);
  //const float padding=0.01; const float ydivide=0.2;
  const float padding=1e-5; const float ydivide=0.3;
  myC->GetPad(1)->SetPad( padding, ydivide + padding, 1-padding, 1-padding);
  myC->GetPad(2)->SetPad( padding, padding, 1-padding, ydivide-padding);
  myC->GetPad(1)->SetLogy(1);
  myC->GetPad(1)->SetRightMargin(.05);
  myC->GetPad(2)->SetRightMargin(.05);
  myC->GetPad(2)->SetBottomMargin(.4);
  myC->GetPad(1)->Modified();
  myC->GetPad(2)->Modified();
  myC->cd(1);
  gPad->SetBottomMargin(small);
  gPad->Modified();

  renormBins(SIGplot,-1);
  renormBins(SLplot,-1);
  SIGplot->GetYaxis()->SetTitleSize(0.07);
  SIGplot->GetYaxis()->SetTitleOffset(1.);
  SIGplot->Draw("HIST E");
  SLplot->SetMarkerStyle(kFullTriangleUp);
  SLplot->Draw("HIST E SAME");
  TH1D* myRatio;
  if(vb) myRatio = new TH1D("ratio", "ratio", nvarbins, varbins);
  else  myRatio = new TH1D("ratio", "ratio", nbins,low,high);
  myRatio->Sumw2();
  myRatio->SetLineColor(kBlack);
  myRatio->SetMarkerColor(kBlack);
  myRatio->Divide(SIGplot,SLplot);
  myRatio->SetMinimum(0);
  myRatio->SetMaximum(2.5);
  myRatio->GetYaxis()->SetNdivisions(200 + int(ratioMax-ratioMin)+1);    //set ticks ; to be seen if this really works
  //myRatio->GetYaxis()->SetNdivisions(400);
  myRatio->GetYaxis()->SetLabelSize(0.1); //make y label bigger
  myRatio->GetXaxis()->SetLabelSize(0.1); //make y label bigger
  myRatio->GetXaxis()->SetTitleOffset(1.1);
  myRatio->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]"); //make y label bigger
  myRatio->GetXaxis()->SetLabelSize(0.12);
  myRatio->GetXaxis()->SetLabelOffset(0.04);
  myRatio->GetXaxis()->SetTitleSize(0.16);
  myRatio->GetYaxis()->SetTitle("Ratio ");
  myRatio->GetYaxis()->SetTitleSize(0.16);
  myRatio->GetYaxis()->SetTitleOffset(.43);
  myC->cd(2);
  gPad->SetTopMargin(small);
  gPad->SetTickx();
  gPad->Modified();

  myRatio->Draw();
  TLine* myLine;
  myLine = new TLine(low, 1, high, 1);
  myLine->Draw();
  myC->cd(1);
  //TLatex* mytext = new TLatex(3.570061,23.08044,"CMS Preliminary");
  TLatex* mytext = new TLatex(3.570061,23.08044,"CMS Simulation");
  mytext->SetNDC();
  mytext->SetTextAlign(13);
  mytext->SetX(0.6);
  mytext->SetY(0.9);
  mytext->SetTextFont(42);
  mytext->SetTextSizePixels(24);
  mytext->Draw();
  TLegend* myLegend = new TLegend(0.615,0.84,0.81,0.7);
  myLegend->SetFillColor(0);
  myLegend->SetBorderSize(0);
  myLegend->SetLineStyle(0);
  myLegend->SetTextFont(42);
  myLegend->SetFillStyle(0);
  myLegend->SetTextSize(0.045);
  if(sample=="ttbar"){
    myLegend->AddEntry(SIGplot,"t#bar{t}, 0 leptons", "lp");
    myLegend->AddEntry(SLplot, "t#bar{t}, 1 lepton", "lp");
  }
  else if(sample=="wjets"){
    myLegend->AddEntry(SIGplot,"wjets, 0 leptons", "lp");
    myLegend->AddEntry(SLplot, "wjets, 1 lepton", "lp");
  }
  else if(sample=="singletop"){
    myLegend->AddEntry(SIGplot,"single-top, 0 leptons", "lp");
    myLegend->AddEntry(SLplot, "singletop, 1 lepton", "lp");
  }
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
  TCut btagcut = "nbjetsCSVM>=1";
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
  //drawR("minDeltaPhi",0.3,15,50,350,"old_"+btagselection);
  //drawR("minDeltaPhiN",4,15,50,350,btagselection);

  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjets==0&&weight<1000";
  //drawR("minDeltaPhiN",4,15,50,350,"eq0b");

  //now try drawing it in data
  //careful -- there is a bug in the drawR() implementation
  //it doesn't work right unless you turn on doRatioPlot
  drawTotalSM_=true; 
  //usePUweight_=true;
  doRatioPlot(true);
  drawLegend(true);
  setPlotMinimum(0); setPlotMaximum(0.5);
  lumiScale_ = 30.15471; //customize lumiScale_
  resetSamples();
  doData(true);
  TCut util = "(pass_utilityHLT==1 || !isRealData) && weight<1000";
  selection_ =TCut("HT>=400 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjetsCSVM==0")&&util;
  //drawR("minDeltaPhiN",4,10,50,150,"Data_eq0b");
  drawR("minDeltaPhiN",4,15,0,150,"Data_eq0b"); //used in PAS
  
}


void AN2011_r_SLreq() {
  //used for PAS
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
  
  setPlotMaximum(0.5); setPlotMinimum(0);
 
  const int nvarbins=13;
  const float varbins[]={0,15,30,45,80,110,140,170,200,230,260,290,320,350}; //use for PAS
  
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjets==0&&weight<1000";
  drawR("minDeltaPhiN",4, "MET", nvarbins, varbins,"eq0b");
  
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjets>=1&&weight<1000";
  drawR("minDeltaPhiN",4, "MET", nvarbins, varbins,"ge1b");
  
  setPlotMaximum(5); setPlotMinimum(0);
  drawR("minDeltaPhi",0.3, "MET", nvarbins, varbins,"old_ge1b");
  
  
  /*
  //other plots looked at, but not in PAS, for example.
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjets==1&&weight<1000";
  drawR("minDeltaPhiN",4, "MET", nvarbins, varbins,"eq1b");
  
  const int nvarbins2=8;
  const float varbins2[]={0,15,30,45,75,125,175,225, 350};
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjets>=2&&weight<1000";
  drawR("minDeltaPhiN",4, "MET", nvarbins2, varbins2,"ge2b");
  */
}


void AN2011_opt() {

  TCut ge1b = "nbjets>=1";
  TCut ge2b = "nbjets>=2";
  TCut ge3b = "nbjets>=3";

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

void AN2011( TString btagselection="ge1b",const int mode=1, bool logy=false, bool doRatio=false ) {
  /*
.L drawReducedTrees.C++
  */  
  doRatioPlot(doRatio);
  setLogY(logy);
  if(logy)  setPlotMinimum(0.5);
  loadSamples();

  //this mode thing is a bit kludgey
  TString modestring="";
  if (mode==1) {
    usePUweight_=false;
    useHLTeff_=false;
    btagSFweight_="1";
    currentConfig_=configDescriptions_.getDefault(); //completely raw MC
  }
  else if (mode==2 || mode==3) {
    assert(0); //not ready yet
    usePUweight_=true;
    useHLTeff_=true;
    currentConfig_=configDescriptions_.getCorrected(); //JER bias
    if (mode==2) modestring="-JER-PU-HLT";
    else if(mode==3) modestring="-JER-PU-HLT-bSF";
  }
  else assert(0);


  TCut btagcut = "nbjetsCSVM>=1";
  if (mode==1 || mode==2) {
    if ( btagselection=="ge1b") {} //do nothing
    else  if ( btagselection=="ge2b" ) {
      btagcut = "nbjetsCSVM>=2";
    }
    else if ( btagselection=="eq1b" ) {
      btagcut = "nbjetsCSVM==1";
    }
    else if ( btagselection=="ge3b" ) {
      btagcut = "nbjetsCSVM>=3";
    }
    else {
      assert(0);
    }
  }
  else if (mode==3) {
    assert(0);//not ready yet
    if ( btagselection=="ge1b") {
      btagcut="1";
      btagSFweight_="probge1";
    }
    else  if ( btagselection=="ge2b" ) {
      btagcut = "1";
      btagSFweight_="probge2";
    }
    else if ( btagselection=="eq1b" ) {
      assert(0); //not implemented
    }
    else {
      assert(0);
    }
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
  drawPlots(var,nbins,low,high,xtitle,"Arbitrary units", "minDeltaPhiN_looseMET_MConly_"+btagselection+modestring);
  
  return;
  
  // ========= regular N-1 plots

  resetSamples(); //use all samples
  setStackMode(true); //regular stack
  setColorScheme("stack");
  doData(true);
  drawMCErrors_=true;
  
  /*  
  var="minDeltaPhiN"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 20; low=0; high=40;
  //no delta phi cut
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200")&&btagcut;
  drawPlots(var,nbins,low,high,xtitle,"Events/2", "SBandSIG_minDeltaPhiN_"+btagselection+modestring);
  */
  /*
  // n Jets
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN >= 4")&&btagcut;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 8; low=1; high=9;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_njets_"+btagselection+modestring);
  */
  /*
  //HT
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=150 && minDeltaPhiN >= 4")&&btagcut;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 20; low=350; high=1050;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_HT_"+btagselection+modestring);
  */

  /*
  //MET
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 35; low=150; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events/10 GeV", "SBandSIG_MET_"+btagselection+modestring);
  
  //MET distribution with tighter HT cut
  selection_ =TCut("HT>500 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  TString thisytitle = "";
  if(btagselection=="ge2b") { nbins = 7; low=150; high=500; thisytitle = "Events/50 GeV";} 
  else{ nbins = 17; low=150; high=500; thisytitle = "Events/20.6 GeV";}
  drawPlots(var,nbins,low,high,xtitle,thisytitle, "SBandSIG_MET_HT500_"+btagselection+modestring);
  */ 
  
  /*
  //MET distribution with tighter HT cut
  TString thisytitle = "";
  //2011b
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 20; low=200; high=600; thisytitle = "Events";
  drawPlots(var,nbins,low,high,xtitle,thisytitle, "SBandSIG_MET_HT400_"+btagselection+modestring);
  
  //2011b
  selection_ =TCut("HT>=500 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 20; low=200; high=600; thisytitle = "Events";
  drawPlots(var,nbins,low,high,xtitle,thisytitle, "SBandSIG_MET_HT500_"+btagselection+modestring);
  
  //2011b
  selection_ =TCut("HT>=600 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 20; low=200; high=600; thisytitle = "Events";
  drawPlots(var,nbins,low,high,xtitle,thisytitle, "SBandSIG_MET_HT600_"+btagselection+modestring);
  */

  /*
  
  // == finally, draw the signal region only!
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN >= 4")&&btagcut;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 20; low=350; high=1050;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_HT_"+btagselection+modestring);

  var="MET+HT"; xtitle="M_{eff} [GeV]";
  nbins = 20; low=550; high=1550;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_Meff__"+btagselection+modestring);

  selection_ =TCut("HT>=500 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=300 && minDeltaPhiN >= 4")&&btagcut;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 20; low=500; high=1050;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SIGtight_HTtight_"+btagselection+modestring);

  var="MET+HT"; xtitle="M_{eff} [GeV]";
  nbins = 10; low=800; high=1600;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SIGtight_Meff__"+btagselection+modestring);

  // ===== single lepton selection
  */
  /*
  // == MET for the electron sample
  selection_ =TCut("MET>=200 && HT>=400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && nElectrons==1 && nMuons==0 && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=200; high=600;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_1e0mu_"+btagselection+modestring);
  */ 

  /*
  //pT?
  var="eleet1"; xtitle="electron p_{T} [GeV]";
  nbins = 20; low=0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_eleET_1e0mu_"+btagselection+modestring);
  */
  /*
  // == MET for the muon sample
  selection_ =TCut("MET>=200 && HT>=400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && nElectrons==0 && nMuons==1 && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=200; high=600;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_0e1mu_"+btagselection+modestring);
  */
  /*
  var="muonpt1"; xtitle="muon p_{T} [GeV]";
  nbins = 20; low=0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_muonpT_0e1mu_"+btagselection+modestring);
  */
  /*  
  // == MET for the combined sample
  selection_ =TCut("MET>=150 &&cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=150; high=450;
  drawPlots(var,nbins,low,high,xtitle,"Events/20 GeV", "SBandSIG_MET_SL_"+btagselection+modestring);
  
  // == HT for the combined sample
  selection_ =TCut("MET>=150 &&cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100")&&btagcut;
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 20; low=350; high=1050;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_HT_SL_"+btagselection+modestring);

  // == njets for the combined sample
  selection_ =TCut("MET>=150 &&cutHT==1 && cutPV==1 && cutTrigger==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100")&&btagcut;
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 20; low=350; high=1050;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_njets_SL_"+btagselection+modestring);
  
  // == MET for the combined sample (HT>500)
  selection_ =TCut("MET>=150 && HT>500 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=150; high=450;
  drawPlots(var,nbins,low,high,xtitle,"Events/20 GeV", "SBandSIG_MET_SL_HT500_"+btagselection+modestring);
  
  //different scale to compare to Kristen
  selection_ =TCut("MET>=150 && HT>500 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 35; low=150; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_SL_HT500_morebins_"+btagselection+modestring);

  */
  // == take a look at WJets
  //doesn't work with the fancy b tag SF code. And we're not using it anyway
  selection_ ="MET>=200 && cutHT==1 && cutPV==1 && cutTrigger==1 && njets==2 && (nElectrons==0 && nMuons==1) && minDeltaPhiN >= 4 && nbjetsCSVM==1";
  var="MT_Wlep"; xtitle="M_{T} [GeV]";
   nbins = 20; low=0; high=200;
   drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MT_0e1mu_eq1b");
/*
  //MET
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN < 4")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 30; low=150; high=450;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_ldp_"+btagselection+modestring);
 
  */
}

// working but not polished code to draw efficiency over mSugra space
void drawmSugraEfficiency() {
  gROOT->SetStyle("CMS");
  loadSamples();
  clearSamples();
  addSample("mSUGRAtanb40");
  setSearchRegions();

  loadSusyScanHistograms();

  loadSusyCrossSections();

  TFile fout("mSugraEffMap.root","RECREATE");
  gROOT->cd();

  for (unsigned int iregion=0; iregion<searchRegions_.size(); iregion++) {
    currentConfig_ = configDescriptions_.getCorrected();
    selection_ =TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4")
      && TCut(searchRegions_[iregion].metSelection.Data()) && TCut(searchRegions_[iregion].htSelection.Data());
    usePUweight_=true;
    useHLTeff_=true;
    TString bweightstring="probge1";
    TString btagselection = searchRegions_[iregion].btagSelection;
    if (btagselection=="ge2b") {
      bweightstring="probge2";
    }
    else if (btagselection=="ge1b") {}
    else if (btagselection=="ge3b") {
      bweightstring="probge3"; //this one isn't calculated right now, i think
    }
    else {assert(0);}
    btagSFweight_ = bweightstring;

    //these are for mSugra [copied and pasted from getSusyScanYields() ]
    TString varx="m0"; TString xtitle=varx;
    int  nbinsx=210; float lowx=-0.5; float highx=2100-0.5;
    
    TString vary="m12"; TString ytitle=vary;
    int  nbinsy=110; float lowy=-0.5; float highy=1100-0.5;
    TString drawstring = vary+":"+varx;
    
    TTree* thetree = (TTree*) files_[currentConfig_]["mSUGRAtanb40"]->Get("reducedTree");
    
    //get the efficiency = Nreco / Ngen for each subprocess
    //then reweight by the subprocess cross section
    
    const int nsubprocesses =  10;
    vector<TH2D*> raw0; //for Npass*sigma*L
    TH2D totalSigma("totalSigma","",nbinsx,lowx,highx,nbinsy,lowy,highy);
    for (int i=0; i<=nsubprocesses; i++) {
      TString hname="raw0_";
      hname += i;
      raw0.push_back(new TH2D(hname,"raw event counts",nbinsx,lowx,highx,nbinsy,lowy,highy));
      TString thecut = getCutString(kmSugraPlane,"",selection_,"",0,"",i);
      thetree->Project(hname,drawstring,thecut.Data());
    }
    //[for mSugra]
    //at this point, each bin contains Npass_i * sigma_i * L for that (m0,m12)
    //need to divide by Ngen_i for each (m0,m12)
    
    //loop over i and histo bins
    for (int i=0; i<=nsubprocesses; i++) {
      for (map<pair<int,int>, TH1D* >::iterator iscanpoint = scanProcessTotalsMap.begin(); iscanpoint!=scanProcessTotalsMap.end(); ++iscanpoint) {
	int m0=iscanpoint->first.first;
	int m12=iscanpoint->first.second;
	
	TH1D* thishist = scanProcessTotalsMap[make_pair(m0,m12)];
	int thisn = TMath::Nint(thishist->GetBinContent(i)); // n gen
	int bin=  raw0[i]->FindBin(m0,m12);
	double N_i_thispoint = raw0[i]->GetBinContent(bin); // Npass * sigma * L
	double this_sigma = CrossSectionTable_mSUGRAtanb40_->getCrossSection(m0,m12,SUSYProcess(i));
	if (thisn == 0) {
	  if (N_i_thispoint > 0.0000001) cout<<"Possible problem: "<<m0<<" "<<m12<<" "<<i<<" "<< N_i_thispoint<<" "<<thisn<<endl;
	  thisn=1; //prevent divide by zero
	  N_i_thispoint = 0;
	}
	N_i_thispoint /= thisn;
	totalSigma.SetBinContent(bin, totalSigma.GetBinContent(bin) + this_sigma);
	raw0[i]->SetBinContent(bin,N_i_thispoint);
      }
    }
    
    //now we have Npass_i * sigma_i * lumi / Ngen_i 
    //all that is left is to make the sum over i
    susyScanYields theEff;
    for (map<pair<int,int>, TH1D* >::iterator iscanpoint = scanProcessTotalsMap.begin(); iscanpoint!=scanProcessTotalsMap.end(); ++iscanpoint) {
      int nentries=    iscanpoint->second->GetEntries();
      if (nentries >10000 || nentries<9500) continue;
      
      double Nraw = 0;
      
      int bin=  raw0[0]->FindBin(iscanpoint->first.first , iscanpoint->first.second);
      for (unsigned int i=0; i<raw0.size(); i++) {
	Nraw += raw0[i]->GetBinContent(bin);
      }
      //the badly named Nraw is sum_i Npass_i * sigma_i * lumi / Ngen_i = sum_i eff_i * sgima_i * lumi = sum_i Nobs_i = N observed
      
      //total Sigma is already the sum over i
      double crossSectionTimesLumi = totalSigma.GetBinContent(bin) * lumiScale_;
      if (crossSectionTimesLumi >0)    theEff[iscanpoint->first] = make_pair(Nraw / crossSectionTimesLumi,0);
      else if (Nraw < 0.00001 && crossSectionTimesLumi <0.00001) {}//do nothing
      else {cout<<"Something weird? crossSectionTimeLumi is zero!"<<endl;}
    }
    
    renewCanvas();
    thecanvas->cd()->SetRightMargin(0.2);
    if (h2d!=0) delete h2d;

    //redo the binning with what appears to be the actual binning (for better presentation)
    nbinsx=105;  lowx=-0.5; highx=2100-0.5; //m0
    nbinsy=40;  lowy=-0.5;  highy=800-0.5; //m12

    TString hname="mSugraEffMap";
    h2d = new TH2D(hname,hname,nbinsx,lowx,highx,nbinsy,lowy,highy);
    TString opt="colz";
    for (susyScanYields::iterator ieff = theEff.begin(); ieff!=theEff.end(); ++ieff) {
      int bin=    h2d->FindBin( ieff->first.first, ieff->first.second);
      h2d->SetBinContent(bin, ieff->second.first);
    }
    
    h2d->SetXTitle("m_{0} [GeV]");
    h2d->SetYTitle("m_{1/2} [GeV]");
    h2d->Draw(opt);

    thecanvas->SaveAs(hname+"."+btagselection+searchRegions_[iregion].owenId+".eps");

    fout.cd();
    h2d->SetName(hname+"_"+btagselection+searchRegions_[iregion].owenId);
    h2d->Write();
    gROOT->cd();

    //clean stuff up
    for (unsigned int ir=0; ir < raw0.size(); ir++) delete raw0[ir];
  }

  fout.Close();

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

void qcdFit(bool loose = true){
  loadSamples(true);
  savePlots_=false; //don't save eps,png,pdf files
  doOverflowAddition(false);

  //output file 
  TString histfilename = "fitHists.root";
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
  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjets == 0";
  const TCut passmdpn = "minDeltaPhiN>=4";
  const TCut failmdpn = "minDeltaPhiN<4";
  TCut base = HTcut&&TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && weight<1000");// no trigger, met, minDeltaPhiN, btagger

  TString var = "MET";
  int nbins; double low; double high;
  
  selection_ = base&&util&&LSB&&antitag&& passmdpn;
  nbins=1; low=50; high=100;
  drawSimple(var,nbins,low,high,histfilename, "h_LSB_fDP_antib_MET_data", "data");
  drawSimple(var,nbins,low,high,histfilename, "h_LSB_fDP_antib_MET_qcd", "PythiaPUQCD");
  drawPlots( var,nbins,low,high, "", "",      "deleteme");  
  TFile fh1(histfilename, "UPDATE");
  totalnonqcd->SetName("h_LSB_fDP_antib_MET_nonqcd");
  totalnonqcd->Write();
  fh1.Close();

}

void drawQCDreweight(bool loose = true, bool njets50 = false){
  loadSamples(true);
  savePlots_=false; //don't save eps,png,pdf files
  doOverflowAddition(true);

  //for Summer11 result, run drawQCDreweight(true,false) and drawQCDreweight(false,false) with b12same=true.
  bool b12same=true;

  TString var = "";
  if(njets50) var = "njets";
  else var = "njets30";

  //output file 
  TString histfilename = "";
  histfilename += "qcdReweight";
  if(loose) histfilename += "_loose_";
  else histfilename += "_tight_";
  if(njets50) histfilename += "50";
  else histfilename += "30";
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
  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjets == 0";
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
  if(njets50) {nbins=4; low=2.5; high=6.5;}
  else{        nbins=5; low=2.5; high=7.5;}
  selection_ = base && LSB && util && failmdpn && antitag;
  drawSimple(var,nbins,low,high,histfilename, "h_LSB_fDP_antib_njets_data", "data");
  drawSimple(var,nbins,low,high,histfilename, "h_LSB_fDP_antib_njets_qcd", "PythiaPUQCD");
  drawPlots( var,nbins,low,high, "", "",      "deleteme");  
  TFile fh1(histfilename, "UPDATE");
  totalnonqcd->SetName("h_LSB_fDP_antib_njets_nonqcd");
  totalnonqcd->Write();
  fh1.Close();
  selection_ = base && LSB && util && passmdpn && antitag;
  drawSimple(var,nbins,low,high,histfilename, "h_LSB_pDP_antib_njets_qcd", "PythiaPUQCD");
  
  //-- SB plots
  //-----------
  lumiScale_ = 1143.;
  
  if(njets50) { nbins=4; low=2.5; high=6.5;} 
  else{         nbins=4; low=2.5; high=6.5;}
  selection_ =  base && SB && phys && failmdpn && ge2b;
  drawSimple(var,nbins,low,high,histfilename, "h_SB_fDP_ge2b_njets_data", "data");
  drawSimple(var,nbins,low,high,histfilename, "h_SB_fDP_ge2b_njets_qcd", "PythiaPUQCD");
  drawPlots( var,nbins,low,high, "", "",      "deleteme");
  TFile fh3(histfilename, "UPDATE");
  totalnonqcd->SetName("h_SB_fDP_ge2b_njets_nonqcd");
  totalnonqcd->Write();
  fh3.Close();
  selection_ = base && SB && phys && passmdpn && ge2b;
  drawSimple(var,nbins,low,high,histfilename, "h_SB_pDP_ge2b_njets_qcd", "PythiaPUQCD");
  
  if(!b12same){
    if(njets50) { nbins=4; low=2.5; high=6.5;} 
    else{         nbins=5; low=2.5; high=7.5;}
  }
  selection_ =  base && SB && phys && failmdpn && ge1b;
  drawSimple(var,nbins,low,high,histfilename, "h_SB_fDP_ge1b_njets_data", "data");
  drawSimple(var,nbins,low,high,histfilename, "h_SB_fDP_ge1b_njets_qcd", "PythiaPUQCD");
  drawPlots( var,nbins,low,high, "", "",      "deleteme");   
  TFile fh2(histfilename, "UPDATE");
  totalnonqcd->SetName("h_SB_fDP_ge1b_njets_nonqcd");
  totalnonqcd->Write();
  fh2.Close(); 
  selection_ = base && SB && phys && passmdpn && ge1b;
  drawSimple(var,nbins,low,high,histfilename, "h_SB_pDP_ge1b_njets_qcd", "PythiaPUQCD");

  //-- SIG plots
  //------------
  if(njets50){
    if(loose){ nbins=3; low=2.5; high=5.5; }
    else{      nbins=2; low=2.5; high=4.5; }
  }
  else{
    if(loose){ nbins=3; low=2.5; high=5.5; }
    else{      nbins=2; low=2.5; high=4.5; }
  }
  selection_ =  base && SIG && phys && failmdpn && ge2b;
  drawSimple(var,nbins,low,high,histfilename, "h_SIG_fDP_ge2b_njets_data", "data");
  drawSimple(var,nbins,low,high,histfilename, "h_SIG_fDP_ge2b_njets_qcd", "PythiaPUQCD");
  drawPlots(var,nbins,low,high, "", "",      "deleteme");
  TFile fh5(histfilename, "UPDATE");
  totalnonqcd->SetName("h_SIG_fDP_ge2b_njets_nonqcd");
  totalnonqcd->Write();
  fh5.Close();
  selection_ = base && SIG && phys && passmdpn && ge2b;
  drawSimple(var,nbins,low,high,histfilename, "h_SIG_pDP_ge2b_njets_qcd", "PythiaPUQCD");

  if(!b12same){
    if(njets50){
      if(loose){ nbins=4; low=2.5; high=6.5; }
      else{      nbins=2; low=2.5; high=4.5; }
    }
    else{
      if(loose){ nbins=5; low=2.5; high=7.5; }
      else{      nbins=4; low=2.5; high=6.5; }
    }
  }
  selection_ =  base && SIG && phys && failmdpn && ge1b;
  drawSimple(var,nbins,low,high,histfilename, "h_SIG_fDP_ge1b_njets_data", "data");
  drawSimple(var,nbins,low,high,histfilename, "h_SIG_fDP_ge1b_njets_qcd", "PythiaPUQCD");
  drawPlots(var,nbins,low,high, "", "",      "deleteme");
  TFile fh4(histfilename, "UPDATE");
  totalnonqcd->SetName("h_SIG_fDP_ge1b_njets_nonqcd");
  totalnonqcd->Write();
  fh4.Close();
  selection_ = base && SIG && phys && passmdpn && ge1b;
  drawSimple(var,nbins,low,high,histfilename, "h_SIG_pDP_ge1b_njets_qcd", "PythiaPUQCD");
  
}



void ptresABCD(bool loose = true){
  doOverflowAddition(true);
  loadSamples();
  clearSamples();
  addSample("PythiaPUQCD");
  
  const TCut LSB = "MET>=50 && MET<100";
  const TCut SB = "MET>=150 && MET<200";
  TCut SIG;
  if(loose) SIG = "MET>=200";
  else SIG =  "MET>=300";
  TCut HTcut;
  if(loose) HTcut = "HT>=350";
  else HTcut = "HT>=500";
  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjets == 0";
  const TCut passmdpn = "minDeltaPhiN>=4";
  const TCut failmdpn = "minDeltaPhiN<4";//.1412212
  //const TCut passmdpn = "minDeltaPhiN15>=2.7";
  //const TCut failmdpn = "minDeltaPhiN15<2.7";//.142238
  //const TCut passmdpn = "minDeltaPhiN5>=7.9";
  //const TCut failmdpn = "minDeltaPhiN5<7.9";
  TCut base = HTcut&&TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && weight<1000");// no trigger, met, minDeltaPhiN, btagger
  
  double A, B, C, D, Aerr, Berr, Cerr, Derr;

  TString sampleOfInterest = ("PythiaPUQCD");
  TString var="HT"; TString xtitle=var;
  int nbins=1; double low=0; double high=50000; //ht overflow is on
  

  //A
  selection_ = base&&LSB&&failmdpn&&antitag;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_A");
  A=getIntegral(sampleOfInterest);
  Aerr=getIntegralErr(sampleOfInterest);

  //B
  selection_ = base&&LSB&&passmdpn&&antitag;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_A");
  B=getIntegral(sampleOfInterest);
  Berr=getIntegralErr(sampleOfInterest);

  cout << "**************" << endl;
  cout << "B/A: " << B/A << endl;
  cout << "**************" << endl;
 

  //C
  selection_ = base&&SIG&&passmdpn&&ge1b;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_A");
  C=getIntegral(sampleOfInterest);
  Cerr=getIntegralErr(sampleOfInterest);

  //D
  selection_ = base&&SIG&&failmdpn&&ge1b;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_A");
  D=getIntegral(sampleOfInterest);
  Derr=getIntegralErr(sampleOfInterest);

  double numerr=jmt::errAtimesB(B,Berr,D,Derr);
  double num = B*D;
  double estimate = num / A;
  double estimateerr= jmt::errAoverB(num,numerr,A,Aerr);

  cout << "prediction: " << estimate << " +- " << estimateerr << endl;
  cout << "MC truth  : " << C << " +- " << Cerr << endl;
  
  double num1 = C-estimate;
  double perc = num1/C;
  double percerr = jmt::errAoverB(estimate, estimateerr,C,Cerr);
  cout << "percent   : " << perc*100. << " +- " << percerr*100. << endl;
  cout << "combined  : " << sqrt(perc*perc + percerr*percerr)*100. << endl;

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

  lumiScale_ = 30.15471; 

  //TString btagselection="antib";
  //TCut btagcut = "nbjets==0";
  
  TCut util = "(pass_utilityHLT==1 || !isRealData) && weight<1000";
  TCut LSB = "MET>=50 && MET<100";
  TCut runN = "runNumber>=175832";

  const  TCut ge1b =  "nbjetsCSVM >= 1";
  const  TCut ge2b =  "nbjetsCSVM >= 2";
  const  TCut ge3b =  "nbjetsCSVM >= 3";
  const  TCut eq1b =  "nbjetsCSVM == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjetsCSVM == 0";
 
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
  else if (ibtag==5) {
    theBTaggingCut = ge3b;
    btagstring = "ge3b";
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
  selection_ =TCut("HT>=400 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "nGoodPV", 20, 0.5, 20.5, "nGoodPV_"+btagstring);
  //const int nvarbins=9;
  //const float varbins[]={0.5,1.5,2.5,3.5,4.5,5.5,6.5,8.5,11.5,16.5};
  //drawR("minDeltaPhiN", 4, "nGoodPV", nvarbins, varbins, "nGoodPV_"+btagstring);
  //dependence is seen here
  //
  
  
  //prescale vs physics check
  selection_ =TCut("HT>=400 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut&&runN;
  drawSimple("nGoodPV",20,0.5,20.5,"nGoodPV.root", "nGoodPV_"+btagstring+"tag_prescale_data","data");
  selection_ =TCut("HT>=400 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=200")&&theBTaggingCut&&runN;
  drawSimple("nGoodPV",20,0.5,20.5,"nGoodPV.root", "nGoodPV_"+btagstring+"tag_physics_data","data");
  


  //njets
  doOverflowAddition(true);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "njets", 8, 2.5, 10.5, "njets_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=150 && weight<1000")&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "njets",8, 2.5, 10.5, "njets_physics_"+btagstring);
  
  //nbjets -- does not depend on theBTaggingCut!!!
  doOverflowAddition(false);
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB;
  //drawR("minDeltaPhiN", 4, "nbjets", 4, -0.5, 3.5, "nbjets_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=150");
  //drawR("minDeltaPhiN", 4, "nbjets", 4, -0.5, 3.5, "nbjets_physics_"+btagstring);
   
  //runnumber
  //set dataOnly to true
  //should hack after first instance of thecanvas->cd(1); to add:  gPad->SetRightMargin(.1); gPad->Modified();
  doOverflowAddition(false);
  setPlotMaximum(0.5); setPlotMinimum(0);
  selection_ =TCut("HT>=400 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "runNumber", 30, 160403.5, 180000.5, "runNumber_"+btagstring);
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

  //MET fit?
  doOverflowAddition(false);
  selection_ =TCut("HT>=350 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&theBTaggingCut;
  //drawR("minDeltaPhiN", 4, "MET", 1, 50, 100, "MET_"+btagstring);

}
 
void LSB_PVreweight(){
  loadSamples();
  doOverflowAddition(true);
  savePlots_=false;

  TString histfilename = "dummy.root";
  TFile fh(histfilename,"RECREATE");//will delete old root file if present 
  fh.Close(); //going to open only when needed 

  TCut HTcut = "HT>=600";

  //PV binning
  //int pvnbins=12;
  //double pvlow = 0.5;
  //double pvhigh = 12.5;
  
  int pvnbins=10;
  float pvbins[]={0.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,10.5,11.5,12.5};//3.5/fb. used with allPrescale=false and 10 run bins. 10 pv bins.
  //float pvbins[]=  {0.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,10.5,11.5,12.5,13.5,14.5,15.5}; //0.150042 +- 0.00443286. 13 pv bins.
  //float pvbins[]=  {0.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5};//0.150014 +- 0.00443216
  //float pvbins[]=    {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5};//0.149939 +- 0.00443764

  
  //make a plot to form the run number binning
  selection_=TCut("HT>999");
  drawSimple("runNumber",1, 160403.5, 180000.5,histfilename, "h_binningDummy", "data");
  TH1D* hbinning = (TH1D*)hinteractive->Clone("hbinning");
  TH1D* hpf =      (TH1D*)hinteractive->Clone("hpf");
  TH1D* hpf_RW =   (TH1D*)hinteractive->Clone("hpf_RW");
  hbinning->Reset(); //because only interested in binning
  hpf     ->Reset();
  hpf_RW  ->Reset();

  //PV distribution in physics --> used for reweighting
  selection_ =TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutTrigger==1 && MET>=200")&&HTcut;//note no mdp or btag
  //  drawSimple("nGoodPV",pvnbins,pvlow,pvhigh,histfilename, "nGoodPV_physics_data","data");
  drawSimple("nGoodPV",pvnbins,pvbins,histfilename, "nGoodPV_physics_data","data");
  TH1D* hPVphysics = (TH1D*)hinteractive->Clone("hPVphysics");
  for(int j=1; j<=hPVphysics->GetNbinsX();j++){
    assert(hPVphysics->GetBinContent(j)>0);
  }
  bool allPrescale = true;
  
  for(int i=1; i<=hbinning->GetNbinsX(); i++){
    double runlow = hbinning->GetBinLowEdge(i);
    double runhigh = hbinning->GetBinLowEdge(i+1);
    TString runCutString = "runNumber>=";
    runCutString+=runlow;
    runCutString+=" && runNumber<";
    runCutString+=runhigh;
    TCut runCut = TCut(runCutString);
    
    
    if(allPrescale){
      selection_ =TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=50 && MET<100 && nbjetsCSVM==0 && (pass_utilityHLT || !isRealData)")&&HTcut;//note no mdp or run cut
    }
    else{
      selection_ =TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=50 && MET<100 && nbjetsCSVM==0 && (pass_utilityHLT || !isRealData)")&&HTcut&&runCut;//note no mdp
    }
    //drawSimple("nGoodPV",pvnbins,pvlow,pvhigh,histfilename, "nGoodPV_prescale_data","data");
    drawSimple("nGoodPV",pvnbins,pvbins,histfilename, "nGoodPV_prescale_data","data");
    TH1D* hPVprescale = (TH1D*)hinteractive->Clone("hPVprescale");
    if(hPVprescale->Integral()<1.) continue;
    
    TH1D* hPVprescale_RW = (TH1D*)hPVphysics->Clone("hPVprescale_RW");
    //hPVprescale_RW->Scale(hPVprescale->Integral()/hPVphysics->Integral());
    //include error on scaling -- turns out to be negligible but leave it in anyway...
    double s = hPVprescale->Integral()/hPVphysics->Integral();
    double serr = jmt::errAoverB(hPVprescale->Integral(),jmt::errOnIntegral(hPVprescale),hPVphysics->Integral(),jmt::errOnIntegral(hPVphysics));
    for(int j=1; j<=hPVprescale_RW->GetNbinsX(); j++){
      assert(hPVprescale->GetBinContent(j)>0);
      hPVprescale_RW->SetBinContent(j,s*hPVprescale_RW->GetBinContent(j));
      double fillerror = sqrt(serr*serr+hPVprescale_RW->GetBinError(j)*hPVprescale_RW->GetBinError(j));
      cout << j << " " << hPVprescale_RW->GetBinError(j) << " " << fillerror << endl;
      hPVprescale_RW->SetBinError(j,fillerror);
    }
    
    TH1D* hPV_W = (TH1D*)hPVphysics->Clone("hPV_W");
    hPV_W->Reset();
    hPV_W->Divide(hPVprescale_RW,hPVprescale); 
    
    
    selection_ =TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=50 && MET<100 && nbjetsCSVM==0 && (pass_utilityHLT || !isRealData) && minDeltaPhiN>=4")&&HTcut&&runCut;
    //drawSimple("nGoodPV",pvnbins,pvlow,pvhigh,histfilename, "nGoodPVPass_prescale_data","data");
    drawSimple("nGoodPV",pvnbins,pvbins,histfilename, "nGoodPVPass_prescale_data","data");
    TH1D* hPVprescalePass = (TH1D*)hinteractive->Clone("hPVprescalePass");
    
    selection_ =TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=50 && MET<100 && nbjetsCSVM==0 && (pass_utilityHLT || !isRealData) && minDeltaPhiN<4")&&HTcut&&runCut;
    //drawSimple("nGoodPV",pvnbins,pvlow,pvhigh,histfilename, "nGoodPVFail_prescale_data","data");
    drawSimple("nGoodPV",pvnbins,pvbins,histfilename, "nGoodPVFail_prescale_data","data");
    TH1D* hPVprescaleFail = (TH1D*)hinteractive->Clone("hPVprescaleFail");

    TH1D* hPVprescalePass_RW = (TH1D*)hPVprescalePass->Clone("hPVprescalePass_RW");
    TH1D* hPVprescaleFail_RW = (TH1D*)hPVprescaleFail->Clone("hPVprescaleFail_RW");
    hPVprescalePass_RW->Multiply(hPV_W);
    hPVprescaleFail_RW->Multiply(hPV_W);
    
    double p, perr, f, ferr, pf, pferr;
    double p_RW, p_RWerr, f_RW, f_RWerr, pf_RW, pf_RWerr;
    p = hPVprescalePass->Integral();
    perr = jmt::errOnIntegral(hPVprescalePass);
    f = hPVprescaleFail->Integral();
    ferr = jmt::errOnIntegral(hPVprescaleFail);
    pf = p/f;
    pferr = jmt::errAoverB(p,perr,f,ferr);
    if(!(pf>0. && pf<1e9)){
      pf=0;
      pferr=0;
    }
    hpf->SetBinContent(i,pf);
    hpf->SetBinError(i,pferr);
    cout << "PF: " << pf << " +- " << pferr << endl;

    p_RW = hPVprescalePass_RW->Integral();
    p_RWerr = jmt::errOnIntegral(hPVprescalePass_RW);
    f_RW = hPVprescaleFail_RW->Integral();
    f_RWerr = jmt::errOnIntegral(hPVprescaleFail_RW);
    pf_RW = p_RW/f_RW;
    pf_RWerr = jmt::errAoverB(p_RW, p_RWerr, f_RW, f_RWerr);
    if(!(pf_RW>0. && pf_RW<1e9)){
      pf_RW=0;
      pf_RWerr=0;
    }
    hpf_RW->SetBinContent(i,pf_RW);
    hpf_RW->SetBinError(i,pf_RWerr);
    cout << "PF_RW: " << pf_RW << " +- " << pf_RWerr << endl;
  }
  
  hpf_RW->SetLineColor(kRed);
  hpf_RW->SetMarkerColor(kRed);
  hpf->GetXaxis()->SetTitle("runNumber");
  hpf->GetYaxis()->SetTitle("N pass/N fail");
  hpf->Draw("e1");
  hpf_RW->Draw("e1 same");

  TFile fh2(histfilename,"UPDATE");
  hpf_RW->Write();
  hpf->Write();
  fh2.Close(); 

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
  TCut theBTaggingCut = "nbjets==0";
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
  //TCut btagcut = "nbjets==0"; 

  TCut util = "pass_utilityHLT_HT300>=1 && weight<1000";
  TCut LSB = "MET>=50 && MET<100";

  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjets == 0";

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
  //TCut btagcut = "nbjets==0";
  
  TCut util = "pass_utilityHLT_HT300>=1 && weight<1000";
  TCut LSB = "MET>=50 && MET<100";
  TCut SB = "MET>=150 && MET<200";
  TCut SIG = "MET>=200";
  TCut highMET = "MET>=150";

  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjets == 0";
 
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

  /*
  //NJETS FOR JOSH
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&LSB&&theBTaggingCut;
  drawR("minDeltaPhiN", 4, "njets", 4, 2.5, 6.5, "njets_joshBin_lsb_"+btagstring);
  selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&highMET&&theBTaggingCut;
  drawR("minDeltaPhiN", 4, "njets", 4, 2.5, 6.5, "njets_joshBin_highMET_"+btagstring);
  */  

  //setPlotMaximum(0.35); setPlotMinimum(0);
  //selection_ =TCut("HT>=350 && cutPV==1 && cutEleVeto==1 && cutMuVeto==1 && cut3Jets==1")&&util&&LSB&&theBTaggingCut;
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

/* jmt -- deprecated Nov 2011
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
*/

/* deprecated function (jmt -- November 2011)
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
    tree->Project("dummyhist","HT",getCutString(true).Data());
    cout<<"N_SIG = "<<dummyhist.Integral()<<endl;

    dummyhist.Reset();
    if (dummyhist.Integral() != 0) assert(0);

    TCut METD = "MET>=150 && MET<5000";
    TCut theDSelection = baseSelection && passCleaning && failMinDeltaPhi && theBTaggingCut && METD;
    selection_ = theDSelection.GetTitle();
    tree->Project("dummyhist","HT",getCutString(true).Data());
    double nd = dummyhist.Integral();
    cout<<"N_D = "<<nd<<endl;

    TCut META = "MET>=100 && MET<150";
    TCut theASelection = baseSelection && passCleaning && failMinDeltaPhi && theBTaggingCut && META;
    selection_ = theASelection.GetTitle();
    tree->Project("dummyhist","HT",getCutString(true).Data());
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
*/

//jmt -- don't feel like fixing getCutString here and we're not doing flvHistReweighting anyway...
/*
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
*/

/* jmt -- 14 Nov 2011 -- commenting out a huge chunk that is very old (2010 analysis)

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
	
	//  *** nominal selection ***   
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
	
	//  *** T2 selection 
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

	//   *** fail mdp selection ***
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

//end of large commented out region */
