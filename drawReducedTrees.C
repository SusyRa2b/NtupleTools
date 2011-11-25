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
#include <iomanip>
#include <map>
#include <set>
	

//*** AFTER SUMMER
//***************************

  //TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-35d/";
  //TString inputPath = "/cu1/joshmt/reducedTrees/test/";
TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-35e/";
TString dataInputPath =  "/cu2/ra2b/reducedTrees/V00-02-35d/";

//-- reducedTrees for Oct 25 SUSY meeting. 3464.581/pb. 
//TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-35a/";
//TString dataInputPath =  "/cu2/ra2b/reducedTrees/V00-02-35a/";


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
  //   TString inputPath = "/cu2/ra2b/reducedTrees/V00-02-25c_fullpf2pat/"; //LM9 with correct pdf weights
//TString inputPath = "/cu2/joshmt/reducedTrees/V00-02-25c_fullpf2pat/"; //with correct pdf weights
//TString inputPath = "/home/joshmt/";//path for MC
//TString dataInputPath = "/cu2/ra2b/reducedTrees/V00-02-24_fullpf2pat/"; //sym links to V00-02-05_v3

//the cutdesc string is now defined in loadSamples()

//double lumiScale_ = 1091.891;
//double lumiScale_ = 1143; //official summer conf lumi
//double lumiScale_ = 3464.581;//oct25
double lumiScale_ = 4683.719;//nov4

#include "drawReducedTrees.h"

const bool reweightLSBdata_=true; //whether or not LSB data is reweighted based on PV distribution
const bool useScaleFactors_=true; //whether or not to use MC scale factors when doing subtraction for data-driven estimates

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
    bweightstring="probge3"; 
  }
  else {assert(0);}

  TCut baseline = "cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passCleaning==1";
  if (isSL) baseline = "cutPV==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0 && MT_Wlep<100 && passCleaning==1";


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
  TH2D* totalPdfWeightsCTEQ = (TH2D*) files_[currentConfig_][sampleOfInterest]->Get("scanProcessTotalsCTEQ");
  const  float nominalweight=   totalPdfWeightsCTEQ->GetBinContent(10,0);

  TH2D* totalPdfWeightsMSTW = (TH2D*) files_[currentConfig_][sampleOfInterest]->Get("scanProcessTotalsMSTW");
  cout<< nominalweight<<"\t"<<   totalPdfWeightsMSTW->GetBinContent(10,0)<<endl;

  TH2D* totalPdfWeightsNNPDF = (TH2D*) files_[currentConfig_][sampleOfInterest]->Get("scanProcessTotalsNNPDF");
  cout<< nominalweight<<"\t"<<   totalPdfWeightsNNPDF->GetBinContent(10,0)<<endl;

  currentConfig_=configDescriptions_.getCorrected(); //add  JER bias
  //the logic here used to be very fragile, so try to ensure that this is the correct thing
  //(things are improved a bit now, but still a good check to make)
  assert( currentConfig_.Contains("JERbias") && currentConfig_.Contains("JES0") && currentConfig_.Contains("METunc0")
	  && currentConfig_.Contains("PUunc0")&& currentConfig_.Contains("BTagEff0")&& currentConfig_.Contains("HLTEff0") );

  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  cout<<"add JER = "<<getIntegral(sampleOfInterest)<<" +/- "<<getIntegralErr(sampleOfInterest) <<endl;
  results.yield_JER = getIntegral(sampleOfInterest);

  //this will be used for PDF uncertainties
  TTree* thetree = (TTree*) files_[currentConfig_][sampleOfInterest]->Get("reducedTree");

  usePUweight_=true;
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  cout<<"add PU weight yield = "<<getIntegral(sampleOfInterest)<<" +/- "<<getIntegralErr(sampleOfInterest) <<endl;
  results.yield_JER_PU = getIntegral(sampleOfInterest);

  useHLTeff_=true;
  drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
  cout<<"add HLT eff = "<<getIntegral(sampleOfInterest)<<" +/- "<<getIntegralErr(sampleOfInterest) <<endl;
  results.yield_JER_PU_HLT = getIntegral(sampleOfInterest);

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
    float totalweight=   totalPdfWeightsCTEQ->GetBinContent(10,i); //used to be i+1
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
    float totalweight=   totalPdfWeightsMSTW->GetBinContent(10,i);//was i+1
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
      float totalweight=   totalPdfWeightsNNPDF->GetBinContent(10,i); //used to be i+1
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

//not only for LM9; can also be used for other samples e.g. ttbar. Just change sampleOfInterest
void runSystematics2011_LM9(  TString sampleOfInterest="LM9" ) {
  setSearchRegions();

  //  signalSystematics2011(sbRegions_[0],false,false,"T1bbbb",1175,75);
  //     return;

  

  vector<ofstream*> textfiles;
  vector<ofstream*> textfiles_human;
  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    char effoutput[500];
    sprintf(effoutput,"signalSyst.%s.%s%s.dat",sampleOfInterest.Data(),searchRegions_[i].btagSelection.Data(),searchRegions_[i].owenId.Data());

    textfiles.push_back( new ofstream(effoutput));
    sprintf(effoutput,"signalSystDetails.%s.%s%s.dat",sampleOfInterest.Data(),searchRegions_[i].btagSelection.Data(),searchRegions_[i].owenId.Data());
    textfiles_human.push_back( new ofstream(effoutput));
  }

  for (unsigned int i=0; i< searchRegions_.size(); i++) {
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
    (*textfiles_human[i])<<" Raw yield                &"<<TMath::Nint(SIG.rawYield)<<endl;
    (*textfiles_human[i])<<" After JER correction     &"<<TMath::Nint(SIG.yield_JER)<<endl;
    (*textfiles_human[i])<<" After PU reweighting     &"<<TMath::Nint(SIG.yield_JER_PU)<<endl;
    (*textfiles_human[i])<<" After HLT correction     &"<<TMath::Nint(SIG.yield_JER_PU_HLT)<<endl;
    (*textfiles_human[i])<<" After btag SF            &"<<TMath::Nint(SIG.rawYield * SIG.effCorr)<<endl;
    (*textfiles_human[i])<<" \\hline                   &"<<endl;
    (*textfiles_human[i])<<" Total correction         &"<< setprecision(2)<<SIG.effCorr<<endl;
    (*textfiles_human[i])<<" \\hline \\hline            &"<<endl;
    (*textfiles_human[i])<<"   JES                    &"<<setprecision(3)<<100*SIG.value("JES")<<endl;
    (*textfiles_human[i])<<"   JER                    &"<<100*SIG.value("JER")<<endl;
    (*textfiles_human[i])<<"   Unclustered energy     &"<<100*SIG.value("MET")<<endl;
    (*textfiles_human[i])<<"   PU                     &"<<100*SIG.value("PU")<<endl;
    (*textfiles_human[i])<<"  \\b tag efficiency       &"<<100*SIG.value("btag")<<endl;
    (*textfiles_human[i])<<"   PDF                    &"<<100*SIG.value("PDF")<<endl;
    (*textfiles_human[i])<<"   kFactor                &"<<100*SIG.value("kFactor")<<endl;
    (*textfiles_human[i])<<"   Trigger efficiency     &"<<100*SIG.value("trigger")<<endl;
    (*textfiles_human[i])<<"   \\MET cleaning          &"<<100*SIG.value("cleaning")<<endl;
    (*textfiles_human[i])<<"   Lepton Veto            &"<<100*SIG.value("LepVeto")<<endl;
    (*textfiles_human[i])<<"   Luminosity             &"<<100*SIG.value("lumi")<<endl;
    (*textfiles_human[i])<<" \\hline                   &"<<endl;
    (*textfiles_human[i])<<" Total systematic         &"<<setprecision(3)<<SIG.totalSystematic()<<endl;
    (*textfiles_human[i])<<" \\hline                   &"<<endl;

  }

  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    textfiles.at(i)->close();
    textfiles_human.at(i)->close();

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


  if (useScaleFactors_) {
    bcut="1";
    usePUweight_=true;
    useHLTeff_=true;
    currentConfig_=configDescriptions_.getCorrected(); //add JERbias

    if (btagselection=="ge2b") {
      btagSFweight_="probge2";
    }
    else if (btagselection=="ge1b") {
      btagSFweight_="probge1";
    }
    else if (btagselection=="ge3b") {
      btagSFweight_="probge3";
    }
    else {assert(0);}
  }

  
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
  
  TCut baseline = "cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passCleaning==1";
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
  TCut baseline = "cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passCleaning==1";
  if (isSL) baseline = "cutPV==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && MT_Wlep>=0 && MT_Wlep<100 && passCleaning==1";
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

//now try to use the 'modular' version (1 box at a time)
void runSystematics2011_scanBoxByBox(TString sampleOfInterest, const unsigned int iregion, const unsigned int ibox) {

  loadSamples();
  clearSamples();
  addSample(sampleOfInterest);
  
  setSearchRegions();

  if (sampleOfInterest.Contains("SUGRA"))  loadSusyScanHistograms();
  
  if (iregion>=searchRegions_.size() ) { cout<<"There are only "<<searchRegions_.size()<<" search regions!"<<endl; return;}
 
  if (ibox >= 6) { cout<<"There are only 6 boxes!"<<endl; return;}

  SearchRegion region = ( ibox%2 == 0) ? sbRegions_[iregion] : searchRegions_[iregion];

  bool isSL  = (ibox == 2 || ibox==3) ? true : false;
  bool isLDP = (ibox == 4 || ibox==5) ? true : false;

  map<pair<int,int>, SignalEffData> sedmap = runTH2Syst2011_mSugra(region,isSL,isLDP,sampleOfInterest);

  TString id;
  id.Form("%s_%d_%d",sampleOfInterest.Data(),iregion,ibox);

  writeSignalEffDataMapToFiles(sedmap, id);

}

void loadSystematics2011_scanBoxByBox(TString sampleOfInterest, const unsigned int iregion) {
 
  loadSamples();
  clearSamples();
  addSample(sampleOfInterest);
  setSearchRegions();
  if (sampleOfInterest.Contains("SUGRA"))  loadSusyScanHistograms();
  if (iregion>=searchRegions_.size() ) { cout<<"There are only "<<searchRegions_.size()<<" search regions!"<<endl; return;}
  //open the output files
  //-- these are the text, root files for Owen et al
  char effoutput[500];
  sprintf(effoutput,"signalSyst.%s.%s%s.dat",sampleOfInterest.Data(),searchRegions_[iregion].btagSelection.Data(),searchRegions_[iregion].owenId.Data());
  ofstream*  textfiles= new ofstream(effoutput);
  sprintf(effoutput,"RA2b.%s.%s%s.root",sampleOfInterest.Data(),searchRegions_[iregion].btagSelection.Data(),searchRegions_[iregion].owenId.Data());
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

  TString id;
  int ibox=0;
  id.Form("%s_%d_%d",sampleOfInterest.Data(),iregion,ibox++);
  map<pair<int,int>, SignalEffData> SB =  loadSignalEffDataMapFromFiles(id);
  id.Form("%s_%d_%d",sampleOfInterest.Data(),iregion,ibox++);
  map<pair<int,int>, SignalEffData> SIG = loadSignalEffDataMapFromFiles(id);
  
  id.Form("%s_%d_%d",sampleOfInterest.Data(),iregion,ibox++);
  map<pair<int,int>, SignalEffData> SBSL = loadSignalEffDataMapFromFiles(id);
  id.Form("%s_%d_%d",sampleOfInterest.Data(),iregion,ibox++);
  map<pair<int,int>, SignalEffData> SIGSL =loadSignalEffDataMapFromFiles(id);
  
  id.Form("%s_%d_%d",sampleOfInterest.Data(),iregion,ibox++);
  map<pair<int,int>, SignalEffData> SBLDP =loadSignalEffDataMapFromFiles(id);
  id.Form("%s_%d_%d",sampleOfInterest.Data(),iregion,ibox++);
  map<pair<int,int>, SignalEffData> SIGLDP =loadSignalEffDataMapFromFiles(id);
  
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
  TString hsuffix = sampleOfInterest; hsuffix+="_"; hsuffix+=searchRegions_[iregion].btagSelection; hsuffix+=searchRegions_[iregion].owenId;
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

std::pair<double,double> ABCD_njetRW(TString phys0bcontrol, TString Acutjm, TString Bcutjm, TString Ccutjm, TString Dcutjm){
  cout << "Running closure with MC jet multiplicity reweighted to match data." << endl;
  
  TString histfilename = "dummy.root";
  TFile fh(histfilename,"RECREATE");//will delete old root file if present 
  fh.Close(); //going to open only when needed 

  double holdlumiscale = lumiScale_;
  TString btagSFweight_nom = btagSFweight_;

  //going to use the same binning for the physics and prescaled samples
  //this isn't necessary, but it is simpler
  int jmnbins = 4;
  float jmbins[] = {2.5,3.5,4.5,5.5,6.5};

  //0b control physics
  if(useScaleFactors_) btagSFweight_ = "prob0";
  phys0bcontrol+= " && minDeltaPhiN<4";
  selection_ = phys0bcontrol;
  drawSimple("njets30",jmnbins,jmbins,"dummy.root", "","data");
  TH1D* hJMphysicsData = (TH1D*)hinteractive->Clone("hJMphysicsData");
  drawSimple("njets30",jmnbins,jmbins,"dummy.root", "","PythiaPUQCD");
  TH1D* hJMphysicsQCD = (TH1D*)hinteractive->Clone("hJMphysicsQCD");
  drawPlots("njets30",jmnbins,jmbins, "", "", "deleteme");
  TH1D* hJMphysicsNonQCD = (TH1D*)totalnonqcd->Clone("hJMphysicsNonQCD");
  
  

  //LSB-failmdpn
  lumiScale_ = 30.15471; 
  selection_ = Acutjm;
  drawSimple("njets30",jmnbins,jmbins,"dummy.root", "","PythiaPUQCD");
  TH1D* hJMAqcd = (TH1D*)hinteractive->Clone("hJMAqcd");
  drawSimple("njets30",jmnbins,jmbins,"dummy.root", "","data");
  TH1D* hJMAdata = (TH1D*)hinteractive->Clone("hJMAdata");
  
  //LSB-passmdpn
  lumiScale_ = 30.15471; 
  selection_ = Bcutjm;
  drawSimple("njets30",jmnbins,jmbins,"dummy.root", "","PythiaPUQCD");
  TH1D* hJMBqcd = (TH1D*)hinteractive->Clone("hJMBqcd");

  //D and SIG (aka C)
  lumiScale_ = holdlumiscale;
  btagSFweight_ = btagSFweight_nom;

  selection_ = Dcutjm;
  drawSimple("njets30",jmnbins,jmbins,"dummy.root", "","PythiaPUQCD");
  TH1D* hJMDqcd = (TH1D*)hinteractive->Clone("hJMDqcd");

  selection_ = Ccutjm;
  drawSimple("njets30",jmnbins,jmbins,"dummy.root", "","PythiaPUQCD");
  TH1D* hJMCqcd = (TH1D*)hinteractive->Clone("hJMCqcd");

  //Now we have all the input we need to do the reweighting..

  //weights////////////////////////////////////////////////
  for(int k=1; k<=hJMphysicsData->GetNbinsX(); k++){
    assert(hJMphysicsData->GetBinContent(k)>0);
    assert(hJMphysicsQCD->GetBinContent(k)>0);
    assert(hJMAdata->GetBinContent(k)>0);
    assert(hJMAqcd->GetBinContent(k)>0);
    assert((hJMphysicsData->GetBinContent(k)- hJMphysicsNonQCD->GetBinContent(k))>0);//need to think about this more
  }
  TH1D* hphysW = (TH1D*)hJMphysicsData->Clone("hphysW");
  hphysW->Add(hJMphysicsNonQCD,-1);
  hphysW->Divide(hJMphysicsQCD);
  TH1D* hpresW = (TH1D*)hJMAdata->Clone("hpresW"); //assume no contamination for LSB
  hpresW->Divide(hJMAqcd);
      
  TFile fh1(histfilename, "UPDATE");
  hphysW->Write();
  hpresW->Write();
  fh1.Close();

  //reweighted boxes
  TH1D* hJMAqcdRW = (TH1D*)hJMAqcd->Clone("hJMAqcdRW");
  TH1D* hJMBqcdRW = (TH1D*)hJMBqcd->Clone("hJMBqcdRW");
  TH1D* hJMCqcdRW = (TH1D*)hJMCqcd->Clone("hJMCqcdRW"); 
  TH1D* hJMDqcdRW = (TH1D*)hJMDqcd->Clone("hJMDqcdRW");
  hJMAqcdRW->Multiply(hpresW);
  hJMBqcdRW->Multiply(hpresW);
  hJMCqcdRW->Multiply(hphysW);
  hJMDqcdRW->Multiply(hphysW);

  double Arw, Brw, Crw, Drw;
  Arw = hJMAqcdRW->Integral();
  Brw = hJMBqcdRW->Integral();
  Crw = hJMCqcdRW->Integral();
  Drw = hJMDqcdRW->Integral();


  ////////////////////////////////////////////////////////////stat//
  double rwErr2 = 0;
  for(int k=1; k<=hJMphysicsData->GetNbinsX(); k++){
    double QAk = hJMAqcdRW->GetBinContent(k);
    double QBk = hJMBqcdRW->GetBinContent(k);
    double QCk = hJMCqcdRW->GetBinContent(k);
    double QDk = hJMDqcdRW->GetBinContent(k);
    double Qk =  hJMphysicsQCD->GetBinContent(k);
    double nk = hJMphysicsNonQCD->GetBinContent(k);
    double Rk = hJMAdata->GetBinContent(k);
    double Pk = hJMphysicsData->GetBinContent(k);

    double QAkerr = hJMAqcdRW->GetBinError(k);
    double QBkerr = hJMBqcdRW->GetBinError(k);
    double QCkerr = hJMCqcdRW->GetBinError(k);
    double QDkerr = hJMDqcdRW->GetBinError(k);
    double Qkerr =  hJMphysicsQCD->GetBinError(k);
    double nkerr = hJMphysicsNonQCD->GetBinError(k);
    double Rkerr = hJMAdata->GetBinError(k);
    double Pkerr = hJMphysicsData->GetBinError(k);

    double QAterm=0,QBterm=0,QCterm=0,QDterm=0,Qterm=0, nterm=0, Rterm=0, Pterm=0;
    double X_toppQ=0, X_botpQ=0, X_toppn=0, X_botpn=0, X_toppR=0, X_botpR=0, X_toppP=0, X_botpP=0;
    double X_top = Crw*Arw;
    double X_bot = Drw*Brw;
    
    QAterm = Crw*Arw/Drw/Brw/Brw*QBk*Rk/QAk/QAk;
    QBterm = - Crw*Arw/Drw/Brw/Brw*Rk/QAk;
    QCterm = Arw/Drw/Brw*(Pk-nk)/Qk;
    QDterm = - Crw*Arw/Brw/Drw/Drw*(Pk-nk)/Qk; 

    X_toppQ = - Arw*QCk/Qk/Qk*(Pk-nk);
    X_botpQ = - Brw*QDk/Qk/Qk*(Pk-nk);

    X_toppn = - Arw*QCk/Qk;
    X_botpn = - Brw*QDk/Qk;

    X_toppR = Crw;
    X_botpR = Drw*QBk/QAk;

    X_toppP = Arw*QCk/Qk;
    X_botpP = Brw*QDk/Qk;

    //quotient rule!
    Qterm = (X_bot*X_toppQ - X_top*X_botpQ)/X_bot/X_bot;
    nterm = (X_bot*X_toppn - X_top*X_botpn)/X_bot/X_bot;
    Rterm = (X_bot*X_toppR - X_top*X_botpR)/X_bot/X_bot;
    Pterm = (X_bot*X_toppP - X_top*X_botpP)/X_bot/X_bot;
    
    cout << "njetrw: QA: " << QAterm*QAterm*QAkerr*QAkerr << endl;
    cout << "njetrw: QB: " << QBterm*QBterm*QBkerr*QBkerr << endl;
    cout << "njetrw: QC: " << QCterm*QCterm*QCkerr*QCkerr << endl;
    cout << "njetrw: QD: " << QDterm*QDterm*QDkerr*QDkerr << endl;
    cout << "njetrw: Q: " << Qterm*Qterm*Qkerr*Qkerr << endl;
    cout << "njetrw: n: " << nterm*nterm*nkerr*nkerr << endl;
    cout << "njetrw: R: " << Rterm*Rterm*Rkerr*Rkerr << endl;
    cout << "njetrw: P: " << Pterm*Pterm*Pkerr*Pkerr << endl;

    rwErr2 += QAterm*QAterm*QAkerr*QAkerr + QBterm*QBterm*QBkerr*QBkerr + QCterm*QCterm*QCkerr*QCkerr + QDterm*QDterm*QDkerr*QDkerr;
    rwErr2 += Qterm*Qterm*Qkerr*Qkerr + nterm*nterm*nkerr*nkerr + Rterm*Rterm*Rkerr*Rkerr + Pterm*Pterm*Pkerr*Pkerr;
  }
  //////////////////////////////////////////////////////////////////
  cout << "njetrw: ABCD: " << Arw << " " << Brw << " " << Crw << " " << Drw << endl; 
  cout << "njetrw: BD/A (rw): " << Brw*Drw/Arw << endl;
  cout << "njetrw: C (rw): " << Crw << endl;
  double closure = 1.0-Crw*Arw/Drw/Brw;
  double closurestat = sqrt(rwErr2);
  cout << "njetrw: closure (rw): " << 100.*closure << " +- " << 100.*closurestat << endl;
  cout << "njetrw: Table: $" << Brw << "$ & $" << Arw << "$ & $" << Drw << "$ & $" << Brw*Drw/Arw << "$ & $" << Crw << "$ & $" <<  100.*closure << " \\pm " << 100.*closurestat << "$ \\\\" << endl;
  
  return make_pair( 100*sqrt(pow(closure,2) + pow(closurestat,2)), 0); //estimate in denominator
}

std::pair<double,double> anotherABCD( const SearchRegion & region, bool datamode=false, float subscale=1,float SBshift=0, const TString LSBbsel="==0", float PVCorFactor = 0, bool doNjetRW = false) {
  //kind of awful, but we'll return the estimate for datamode=true but the closure test bias for datamode=false
  
  /*
    .L drawReducedTrees.C++
  */
  doOverflowAddition(true);
  
  TString btagselection = region.btagSelection;
  
  TString owenKey = btagselection;
  owenKey += region.owenId;
  OwenData * myOwen = &(owenMap_[owenKey]);
  const  bool isSIG = region.isSIG;
  
  setStackMode(false);
  doData(true);
  setQuiet(true);
  
  useFlavorHistoryWeights_=false;
  
  loadSamples(); //needs to come first! this is why this should be turned into a class, with this guy in the ctor
  setColorScheme("nostack");
  //clearSamples();
  //addSample("PythiaPUQCD");
  
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
  triggerCutLSB = "pass_utilityHLT==1";
  triggerCut = "cutTrigger==1";
  
  TCut ge1b = "nbjetsCSVM>=1";
  btagSFweight_="1";
  TString btagSFweight=""; //we're going to have to switch this one in and out of the global var
  if (useScaleFactors_) {
    usePUweight_=true;
    useHLTeff_=true;
    currentConfig_=configDescriptions_.getCorrected(); //add JERbias
    
    ge1b="1";
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
  TCut cleaning = "weight<1000 && passCleaning==1";
  
  
  //LSB b-tagging is independent of search region
  TString LSBbtagSFweight="";
  char cutstring0[100];
  sprintf(cutstring0,"nbjetsCSVM%s",LSBbsel.Data());
  TCut LSBbtag = TCut(cutstring0);
  if(useScaleFactors_){
    LSBbtag = "1";
    if(LSBbsel=="==0"){ LSBbtagSFweight = "prob0"; }
    else if(LSBbsel==">=1"){ LSBbtagSFweight = "probge1"; }
    else{ assert(0);}
  }
  btagSFweight_ =  LSBbtagSFweight;
  char cutstring1[100];
  sprintf(cutstring1,"MET>= %.0f && MET < %.0f", 50+SBshift,100+SBshift);
  TCut SBMET = TCut(cutstring1)&&triggerCutLSB&&LSBbtag;
  
  TCut dpcut = "1";//"minDeltaPhiN>=4";
  //  TCut passOther = "deltaPhiMPTcaloMET<2";
  //  TCut failOther = "deltaPhiMPTcaloMET>=2";
  TCut passOther = "minDeltaPhiN>=4";
  TCut failOther = "minDeltaPhiN<4";
  
  var="HT"; xtitle=var;
  nbins=10; low=0; high=5000;
  double A,B,D,SIG,Aerr,Berr,Derr,SIGerr;
  double RLSB_RW = 0;
  double dRLSB_RW = 0;
  
  //used for jet multiplicity reweighting closure test
  TString Acutjm, Bcutjm, Dcutjm, SIGcutjm;
  
  //used in lsb and njet reweighting (the latter adds mdpN<4)
  TCut phys0bcontrol =  TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjetsCSVM==0 && MET>200") && triggerCut && HTcut && cleaning;
  TString phys0bcontrolstring;
  phys0bcontrolstring += phys0bcontrol;
  if(datamode &&  reweightLSBdata_){
    int pvnbins=11;
    float pvbins[]={0.5,2.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,12.5,14.5,16.5};
    
    //physics triggers control sample -- physics triggers, 0b
    //--consider varying MET cut and PV binning!!
    if(useScaleFactors_) btagSFweight_ = "prob0";
    selection_ = phys0bcontrol;
    drawSimple("nGoodPV",pvnbins,pvbins,"dummy", "","data");
    TH1D* hPVphysics = (TH1D*)hinteractive->Clone("hPVphysics");
    btagSFweight_ =  LSBbtagSFweight;
    
    //LSB unweighted
    selection_ = baseline && SBMET && cleaning;
    drawSimple("nGoodPV",pvnbins,pvbins,"dummy", "","data");
    TH1D* hPVprescale = (TH1D*)hinteractive->Clone("hPVprescale");
    
    //calculate weights (physics shape with prescale integral divided by prescale)
    TH1D* hPVprescale_RW = (TH1D*)hPVphysics->Clone("hPVprescale_RW");
    hPVprescale_RW->Scale(hPVprescale->Integral()/hPVphysics->Integral());
    TH1D* hPV_W = (TH1D*)hPVphysics->Clone("hPV_W");
    hPV_W->Reset();
    hPV_W->Divide(hPVprescale_RW,hPVprescale); 
    
    //LSB pass mdp unweighted
    selection_ = baseline && SBMET && passOther && cleaning;
    drawSimple("nGoodPV",pvnbins,pvbins,"dummy", "","data");
    TH1D* hPVprescalePass = (TH1D*)hinteractive->Clone("hPVprescalePass");
    
    //LSB fail mdp unweighted
    selection_ = baseline && SBMET && failOther && cleaning;
    drawSimple("nGoodPV",pvnbins,pvbins,"dummy", "","data");
    TH1D* hPVprescaleFail = (TH1D*)hinteractive->Clone("hPVprescaleFail");
    
    //check for zero entries
    for(int j=1; j<=hPVphysics->GetNbinsX();j++){
      assert(hPVphysics->GetBinContent(j)>0);
      assert(hPVprescalePass->GetBinContent(j)>0);
      assert(hPVprescaleFail->GetBinContent(j)>0);
    }
    
    //weighted LSB pass and fail mdp
    TH1D* hPVprescalePass_RW = (TH1D*)hPVprescalePass->Clone("hPVprescalePass_RW");
    TH1D* hPVprescaleFail_RW = (TH1D*)hPVprescaleFail->Clone("hPVprescaleFail_RW");
    hPVprescalePass_RW->Multiply(hPV_W);
    hPVprescaleFail_RW->Multiply(hPV_W);
    
    ///////////////////////////////////////////////////////////////////////////////stat///////
    double dR2 = 0;
    double myA=0, myB=0;
    for(int k=1; k<=hPVphysics->GetNbinsX(); k++){
      double dNk = hPVphysics->GetBinError(k);
      double dPk = hPVprescalePass->GetBinError(k);
      double dFk = hPVprescaleFail->GetBinError(k);
      double Nk = hPVphysics->GetBinContent(k);
      double Pk = hPVprescalePass->GetBinContent(k);
      double Fk = hPVprescaleFail->GetBinContent(k);
      double Ntot = hPVphysics->Integral();
      double PFtot = hPVprescalePass->Integral()+hPVprescaleFail->Integral();
      
      //pX means derivative w.r.t. X
      double ApN=0, BpN=0;
      double ApP=0, BpP=0;
      double ApF=0, BpF=0;
      double Nterm=0, Pterm=0, Fterm=0;
      for(int i=1; i<=hPVphysics->GetNbinsX(); i++){
	double Ni = hPVphysics->GetBinContent(i);
	double Pi = hPVprescalePass->GetBinContent(i);
	double Fi = hPVprescaleFail->GetBinContent(i);
	
	myB += Ni/Ntot*PFtot/(Pi+Fi)*Pi;
	myA += Ni/Ntot*PFtot/(Pi+Fi)*Fi;

	//part without k
	if(i==k) continue;
	BpN += (-1.)*PFtot/(Pi+Fi)*Pi*Ni/(Ntot*Ntot);
	ApN += (-1.)*PFtot/(Pi+Fi)*Fi*Ni/(Ntot*Ntot);
	BpP += Ni/Ntot*Pi/(Pi+Fi);
	ApP += Ni/Ntot*Fi/(Pi+Fi);
	BpF += Ni/Ntot*Pi/(Pi+Fi);
	ApF += Ni/Ntot*Fi/(Pi+Fi);
      }//i loop
      
      //part with k
      BpN += PFtot/(Pk+Fk)*Pk*(Ntot-Nk)/(Ntot*Ntot);
      ApN += PFtot/(Pk+Fk)*Fk*(Ntot-Nk)/(Ntot*Ntot);
      BpP += Nk/Ntot/(Pk+Fk)/(Pk+Fk)*((Pk+Fk)*(PFtot+Pk)-PFtot*Pk);
      ApP += Nk/Ntot*Fk/(Pk+Fk)/(Pk+Fk)*(Pk+Fk-PFtot);
      BpF += Nk/Ntot*Pk/(Pk+Fk)/(Pk+Fk)*(Pk+Fk-PFtot);
      ApF += Nk/Ntot/(Pk+Fk)/(Pk+Fk)*((Pk+Fk)*(PFtot+Fk)-PFtot*Fk);      
      
      //quotient rule!
      Nterm = (myA*BpN-myB*ApN)/(myA*myA);
      Pterm = (myA*BpP-myB*ApP)/(myA*myA);
      Fterm = (myA*BpF-myB*ApF)/(myA*myA);
      
      dR2 += Nterm*Nterm*dNk*dNk + Pterm*Pterm*dPk*dPk + Fterm*Fterm*dFk*dFk;
      
    }//k loop
    RLSB_RW = myB/myA; //matches hPVprescalePass_RW->Integral()/hPVprescaleFail_RW->Integral()
    dRLSB_RW = sqrt(dR2);
    ///////////////////////////////////////////////////////////////////////////////////////////
    //unweighted
    B = hPVprescalePass->Integral();
    Berr = jmt::errOnIntegral(hPVprescalePass); 
    A = hPVprescaleFail->Integral();
    Aerr = jmt::errOnIntegral(hPVprescaleFail);
  }
  else{
    //A   -- aka 50 - 100 and high MPT,MET
    selection_ = baseline && cleaning && dpcut  && SBMET && failOther; //auto cast to TString seems to work
    Acutjm = selection_;
    drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_A");
    A=getIntegral(sampleOfInterest);
    Aerr=getIntegralErr(sampleOfInterest);
    //B
    selection_ = baseline && cleaning && dpcut  && SBMET && passOther; //auto cast to TString seems to work
    Bcutjm = selection_;
    drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_B");
    B=getIntegral(sampleOfInterest);
    Berr=getIntegralErr(sampleOfInterest);
  }
  
  //Now that we've moved on from the LSB, use search region b-tag requirement
  if(useScaleFactors_) btagSFweight_ = btagSFweight;
  
  //D
  selection_ = baseline && cleaning && dpcut  && SRMET && failOther; //auto cast to TString seems to work
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_D");
  Dcutjm = selection_;
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
  SIGcutjm = selection_;
  drawPlots(var,nbins,low,high,xtitle,"events","qcdstudy_ABCDkludge_row1_C");
  SIG=getIntegral(sampleOfInterest);
  SIGerr=getIntegralErr(sampleOfInterest);
  
  if (datamode) { //for owen 
    if (isSIG)  myOwen->Nsig = SIG;
    else        myOwen->Nsb  = SIG;
  }

  
  if(doNjetRW){ return ABCD_njetRW(phys0bcontrolstring, Acutjm, Bcutjm, SIGcutjm, Dcutjm);}
  
  //now calculate (B/A)*D=R*D
  double RLSB_UW = B/A; //unweighted
  double RLSB_UWerr = jmt::errAoverB(B,Berr,A,Aerr);
  
  double myR, myRerr;
  if(datamode &&  reweightLSBdata_){
    myR = PVCorFactor*(RLSB_RW-RLSB_UW)+RLSB_RW; //PVCorFactor=0 for reweighted, +/-1 for +/-100% correction
    myRerr = dRLSB_RW;
  }
  else{
    myR = RLSB_UW;
    myRerr = RLSB_UWerr; 
  }
  double estimate = myR*(D-Dsub);
  double estimateerr = jmt::errAtimesB(myR, myRerr, D-Dsub, sqrt(Derr*Derr+Dsuberr*Dsuberr));
  double closureStat2 = datamode? 0: jmt::errAoverB(SIG,SIGerr,estimate,estimateerr);
  double R0 = myR;
  double R0err = myRerr;

  TString name = region.btagSelection;
  name += region.owenId;
  name += isSIG ? ", SIG":", SB";
  char output[500];
  if (!datamode) {
    sprintf(output,"%s & %s & %s & %s & %s & %s & $%f \\pm %f$ \\\\ ",name.Data(),
	    jmt::format_nevents(B,Berr).Data(),jmt::format_nevents(A,Aerr).Data(),
	    jmt::format_nevents(D,Derr).Data(),jmt::format_nevents(estimate,estimateerr).Data(),
	    jmt::format_nevents(SIG,SIGerr).Data(),100*(estimate-SIG)/estimate,closureStat2*100);
    
  }
  else {
    if(reweightLSBdata_){
      sprintf(output,"%s & $%f \\pm %f$ & %d & %s & %s   \\\\ %% %f +/- %f",name.Data(),
	      RLSB_RW,dRLSB_RW,
	      TMath::Nint(D), jmt::format_nevents(Dsub,Dsuberr).Data(),
	      jmt::format_nevents(estimate,estimateerr).Data(),R0,R0err);
      cout<<"(qcd) DATA\t";
    }
    else{
      sprintf(output,"%s & %d & %d & %d & %s & %s   \\\\ %% %f +/- %f",name.Data(),
	      TMath::Nint(B),TMath::Nint(A),
	      TMath::Nint(D), jmt::format_nevents(Dsub,Dsuberr).Data(),
	      jmt::format_nevents(estimate,estimateerr).Data(),R0,R0err);
      cout<<"(qcd) DATA\t";
    }
  }
  cout<<output<<endl;
  
  if (datamode)  return make_pair(estimate,estimateerr);
  return make_pair( 100*sqrt(pow((SIG-estimate)/estimate,2) + pow(closureStat2,2)), 0); //estimate in denominator
}


void runClosureTest2011(std::map<TString, std::vector<double> > & syst, bool addnjetrw = true)  {
  
  setSearchRegions();
  assert(sbRegions_.size() == searchRegions_.size());

  for (unsigned int i=0; i<sbRegions_.size(); i++) {

    double sb =   fabs(anotherABCD(sbRegions_[i]).first);
    double sig =  fabs(anotherABCD(searchRegions_[i]).first);
    double sb_rw = 0, sig_rw = 0;
    double finalsb = sb;
    double finalsig = sig;

    if(addnjetrw){
      sb_rw = fabs(anotherABCD(sbRegions_[i],false,1,0,"==0",0,true).first);
      sig_rw = fabs(anotherABCD(searchRegions_[i],false,1,0,"==0",0,true).first);

      finalsb = (sb_rw>sb) ? sb_rw: sb;
      finalsig = (sig_rw>sig) ? sig_rw: sig;

      cout << "njetrw: ************** int " << i << " ****************" << endl;
      cout << "njetrw: sb: " << sb << ", sb_rw: " << sb_rw << endl;
      cout << "njetrw: sig: " << sig << ", sig_rw: " << sig_rw << endl;
      cout << "njetrw: ******************************************" << endl;
    }
    

    syst["Closure"].push_back( finalsb );
    syst["Closure"].push_back( finalsig ); 
  }
  
}

void runClosureTest2011()  {
  std::map<TString, std::vector<double> > dummy;
  
  //cout << "You are starting closure test without considering njet reweighting!" << endl;
  //runClosureTest2011(dummy,false);

  cout << "You are starting closure test that also considers njet reweighting! (this will do both)" << endl;
  runClosureTest2011(dummy,true);

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

  vector<std::pair<double,double> > lsbp;
  vector<std::pair<double,double> > lsbm;

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
  
  if(reweightLSBdata_){
    //do it with LSB RW correction +100%
    cout << "LSB RW correction +100%" << endl;
    for (unsigned int i=0; i<sbRegions_.size(); i++) {
      lsbp.push_back( anotherABCD(sbRegions_[i],true,1,0,"==0",1));
      lsbp.push_back( anotherABCD(searchRegions_[i],true,1,0,"==0",1));
    }
    
    //do it with LSB RW correction -100%
    cout << "LSB RW correction -100%" << endl;
    for (unsigned int i=0; i<sbRegions_.size(); i++) {
      lsbm.push_back( anotherABCD(sbRegions_[i],true,1,0,"==0",-1));
      lsbm.push_back( anotherABCD(searchRegions_[i],true,1,0,"==0",-1));
    }
    
    cout<<" ==== systematics for LSB reweighting ==== "<<endl;
    for (unsigned int j=0; j<n.size(); j++) {
      double var1 = 100*(n[j].first  -lsbp[j].first)/n[j].first;
      double var2 = 100*(n[j].first -lsbm[j].first)/n[j].first;
      if ( fabs(var1)-fabs(var2) > 0.1 ) cout<<j<<" found LSB reweighting systematic size discrepancy"<<endl;
      cout<<var1<<"\t"<<var2<<endl;
      qcdSystErrors["LSBrw"].push_back( fabs(var1)>fabs(var2) ? fabs(var1) : fabs(var2));
    }
    cout<<" ==== END systematics for LSB reweighting ==== "<<endl;
  }


  cout<<"== Cross check with >=1 b instead of exactly 0 b =="<<endl;
  for (unsigned int i=0; i<sbRegions_.size(); i++) {
    anotherABCD(sbRegions_[i],true,1,0,">=1");
    anotherABCD(searchRegions_[i],true,1,0,">=1");
  }
  
  cout<<" == Running QCD closure test =="<<endl;
  runClosureTest2011(qcdSystErrors,true);
  
  cout<<" == QCD systematics summary =="<<endl;
  for (unsigned int j=0; j<n.size(); j++) {
    if( reweightLSBdata_){
      qcdSystErrors["Total"].push_back( sqrt( pow(qcdSystErrors["MCsub"].at(j),2) +  pow(qcdSystErrors["Closure"].at(j),2) + pow(qcdSystErrors["SBshift"].at(j),2) +  pow(qcdSystErrors["LSBrw"].at(j),2)  ));
      cout<<j<<"\t& $"<<qcdSystErrors["MCsub"].at(j)<<"$ & $"<<qcdSystErrors["Closure"].at(j)<<"$ & $"<<qcdSystErrors["SBshift"].at(j)<<"$ & $"<<qcdSystErrors["LSBrw"].at(j) << "$ & $" <<qcdSystErrors["Total"].at(j)<< "$ \\\\" << endl;
    }
    else{
      qcdSystErrors["Total"].push_back( sqrt( pow(qcdSystErrors["MCsub"].at(j),2) +  pow(qcdSystErrors["Closure"].at(j),2)+ pow(qcdSystErrors["SBshift"].at(j),2)));
      cout<<j<<"\t& $"<<qcdSystErrors["MCsub"].at(j)<<"$ & $"<<qcdSystErrors["Closure"].at(j)<<"$ & $"<<qcdSystErrors["SBshift"].at(j)<<"$ & $"<<qcdSystErrors["Total"].at(j)<< "$ \\\\"<< endl;
    }
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
    double zinvscale = 4.68/4.65;
    bool doMean = true;

    if (qcdsubregion.owenId == "Loose" && qcdsubregion.btagSelection=="ge1b") {
      doMean=false;//averaging already done
      zv[0] = 82; ze[0]=20;
      zsbsyst = 0.0;//ze is stat+syst error combined
    }
    else if (qcdsubregion.owenId == "Tight" && qcdsubregion.btagSelection=="ge1b") {
      doMean=false;//averaging already done
      zv[0] = 44; ze[0]=13;
      zsbsyst = 0.0;
    }
    else if (qcdsubregion.owenId == "Loose" && qcdsubregion.btagSelection=="ge2b") {
      doMean=false;//averaging already done
      zv[0] = 14; ze[0]=9;
      zsbsyst = 0.0;
    }
    else if (qcdsubregion.owenId == "Tight" && qcdsubregion.btagSelection=="ge2b") {
      doMean=false;//averaging already done
      zv[0] = 3.8; ze[0]=2.7;
      zsbsyst = 0.0;
    }
    else if (qcdsubregion.owenId == "Loose" && qcdsubregion.btagSelection=="ge3b") {
      doMean=false;//averaging already done
      zv[0] = 1.9; ze[0]=2.8;
      zsbsyst = 0.0;   
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
  TCut cleaning = "weight<1000 && passCleaning==1";
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

  std::vector<double> ttw;
  std::vector<double> p;
  std::vector<double> m;
  
  std::vector<double> qcd;
  std::vector<double> znn;
  std::vector<double> mc;
  //  int i=0;
  cout<<"nominal ttbar"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) ttw.push_back(slABCD(j,true));
 
  if (forOwen) return;

  // for systematics due to QCD subtraction
  //must do   runDataQCD2011(); first (in order to fill syst errors for qcd)
  if (  qcdSystErrors.size()==0) {cout<<"Need to do runDataQCD2011() first!"<<endl; return;}

  cout<<"vary QCD subtraction"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) {
    p.push_back(slABCD(j,true,"QCDup"));
    m.push_back(slABCD(j,true,"QCDdown"));
    if (fabs(p[j]-ttw[j])/ttw[j] - fabs(m[j]-ttw[j])/ttw[j] > 0.01) cout<<"** Difference in size of QCD sub syst!"<<endl;
    qcd.push_back(100*fabs(p[j]-ttw[j])/ttw[j]);
  }

  p.clear(); m.clear();
  //for Znunu subtraction systematics
  cout<<"vary Znn subtraction"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) {
    p.push_back(slABCD(j,true,"Zup"));
    m.push_back(slABCD(j,true,"Zdown"));
    if (fabs(p[j]-ttw[j])/ttw[j] - fabs(m[j]-ttw[j])/ttw[j] > 0.01) cout<<"** Difference in size of Znn sub syst!"<<endl;
    znn.push_back(100*fabs(p[j]-ttw[j])/ttw[j]);
  }

  p.clear(); m.clear();
  //for MC subtraction systematics
  cout<<"vary MC subtraction"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) {
    p.push_back(slABCD(j,true,"MCup"));
    m.push_back(slABCD(j,true,"MCdown"));
    if (fabs(p[j]-ttw[j])/ttw[j] - fabs(m[j]-ttw[j])/ttw[j] > 0.01) cout<<"** Difference in size of MC sub syst!"<<endl;
    mc.push_back(100*fabs(p[j]-ttw[j])/ttw[j]);
  }

  //finally, run the closure test
  std::vector<double> closure;
  cout<<"Running closure tests"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) {
    double allsamples=fabs(slABCD(j)); //mix of samples
    double ttbaronly=fabs(slABCD(j,false,"",true)); //just ttbar
    if(allsamples>ttbaronly)
      closure.push_back(allsamples);
    else closure.push_back(ttbaronly);
  }

  ttbarClosureSyst.clear();
  cout<<" == summary (%) == "<<endl;
  cout<<"\tClosure\tQCD\tZ\tMC\tTotal"<<endl;
  for (unsigned int j=0; j<searchRegions_.size();j++) {
    double totalsyst = sqrt(closure[j]*closure[j] + qcd[j]*qcd[j] + znn[j]*znn[j] + mc[j]*mc[j]);
    cout<<j<<"\t & $"<<closure[j]<<"$ & $"<<qcd[j]<<"$ & $"<<znn[j]<<"$ & $"<<mc[j]<<"$ & $"<<totalsyst << "$ \\\\ " << endl;
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

  TCut util = "pass_utilityHLT==1 && weight<1000";
  //lumiScale_ = 19.23288; //customize lumiScale_
  lumiScale_ = 30.15471; 

  // ========= regular N-1 plots

  resetSamples(); //use all samples
  setStackMode(true); //regular stack
  setColorScheme("stack");
  doData(true);
  drawMCErrors_=true;

  //MET
  selection_ =TCut("cutHT==1 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4 && passCleaning==1")&&util&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 20; low=0; high=200;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_MET_"+btagselection+modestring,0,"GeV");
  setLogY(true);  setPlotMinimum(1e-1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_MET_"+btagselection+modestring,0,"GeV");

  //njets in 50 < MET < 100 region
  selection_ =TCut("cutHT==1 && cutPV==1 && MET>=50 && MET<100 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4 && passCleaning==1")&&util&&btagcut;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 7; low=2; high=9;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_njets_LSB_"+btagselection+modestring);

  selection_ =TCut("cutHT==1 && cutPV==1 && MET>=50 && MET<100 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passCleaning==1")&&util&&btagcut;
  var="minDeltaPhiN"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 20; low=0; high=40;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_minDeltaPhiN_LSB_"+btagselection+modestring);
  setLogY(true);  setPlotMinimum(1e-1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_minDeltaPhiN_LSB_"+btagselection+modestring);


  //look again at MET
  selection_ =TCut("MET>=50 && MET<100 && cutHT==1 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4 && passCleaning==1")&&util&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 10; low=50; high=100;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_METlimited_"+btagselection+modestring);


  //sanity check
/*
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
/*
  doRatioPlot(true);
  selection_ =TCut("HT>=400 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1&&passCleaning==1 && MET<100")&&util&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 10; low=0; high=100;
  ratioMin=0; ratioMax=2;
  setLogY(false); resetPlotMinimum();
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_MET_"+btagselection+modestring,0,"GeV");
*/
/*  
  doRatioPlot(true);
  selection_ =TCut("HT>=400 && cutPV==1 && MET>=50 && MET<100 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1")&&util&&btagcut;
  var="minDeltaPhiN"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 20; low=0; high=40;
  setLogY(false); resetPlotMinimum();
  //drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_minDeltaPhiN_LSB_"+btagselection+modestring);
  setLogY(true);  setPlotMinimum(1e-1);
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_minDeltaPhiN_LSB_"+btagselection+modestring);
  */

/*
  doRatioPlot(true);
  selection_ =TCut("HT>=400 && cutPV==1 && MET>=50 && MET<100 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN < 4 &&passCleaning==1")&&btagcut&&util;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 7; low=2; high=9;
  drawPlots(var,nbins,low,high,xtitle,"Events", "prescaled_njets_"+btagselection+modestring);
*/
 
}


void AN2011_ttbarw( double& r_sl, double& r_sl_err, double& r_nom, double& r_nom_err, 
		    TString btagselection="ge1b", TString HTselection="Loose" , TString samplename="TTbarJets") {

  /*
goal:

show MET distributions for w+jets, ttbar, singletop
in both the SL and nominal samples

this is going to require some serious kludgey stuff because this code has
never been used to plot the same sample with different cuts on top of each 
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

  if (useScaleFactors_) {
    usePUweight_=true;
    useHLTeff_=true;
    currentConfig_=configDescriptions_.getCorrected(); //add JERbias

    //when making the MET comparison plots for the AN, disable b-tag SF for now
    //btagcut="1";   
    //if (btagselection=="ge2b") {
    //  btagSFweight_="probge2";
    //}
    //else if (btagselection=="ge1b") {
    //  btagSFweight_="probge1";
    //}
    //else if (btagselection=="ge3b") {
    //  btagSFweight_="probge3";
    //}
    //else {assert(0);}
  }


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
  const int nvarbins=15;
  const float varbins[]={100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, 600, 650}; //AN and PAS 
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
  //nbins = 40; low=100; high=550; 
  nbins = 55; low=100; high=650; 
  //nbins = 40; low=200; high=600; 
  //nbins = 40; low=150; high=550; //AN v5
  //nbins = 20; low=150; high=550;
  //nbins = 16; low=150; high=550;
  if(vb) drawPlots(var,nvarbins, varbins, xtitle,"Arbitrary units", "SBandSIG_MET_SL_"+sample+"_"+btagselection+"_"+HTselection);
  else drawPlots(var,nbins,low,high,xtitle,"Arbitrary units", "SBandSIG_MET_SL_"+sample+"_"+btagselection+"_"+HTselection);

  ////////////////////////////////////////////////
  //This part sets up the SIG/SB ratio computation
  ////////////////////////////////////////////////

  //the lower boundary on the SB should be set to 200 GeV
  int lowbin = hinteractive->FindBin(200);
  int boundarybin_sb = 0, boundarybin_sig=0;
  //the upper boundary on the SB should be set to 250 GeV
  boundarybin_sb = hinteractive->FindBin(250);

  //the lower boundary on the SR depends on the selection
  if (HTselection=="Tight" && btagselection=="ge1b")  boundarybin_sig = hinteractive->FindBin(500);
  else if (HTselection=="Tight" && btagselection=="ge2b")  boundarybin_sig = hinteractive->FindBin(300);
  else boundarybin_sig = hinteractive->FindBin(250);

  //the upper boundary on the SR region should be the last bin on the plot
  int highbin= hinteractive->GetNbinsX();


  ////////////////////////////////////////
  //This part computes the SL SIG/SB ratio
  ////////////////////////////////////////

  double sl_sb_err=0, sl_sig_err=0;
  double sl_sb = hinteractive->IntegralAndError(lowbin,boundarybin_sb-1,sl_sb_err);//integral from 200 to <250
  double sl_sig = hinteractive->IntegralAndError(boundarybin_sig,highbin,sl_sig_err);//integral from e.g. 300 to last bin
  double r_sl_sigoversb = sl_sig/sl_sb;
  double r_sl_sigoversb_err = jmt::errAoverB(sl_sig,sl_sig_err,sl_sb,sl_sb_err);

  std::cout << "SL: SB content= " << sl_sb << " +/- " << sl_sb_err << std::endl;
  std::cout << "SL: SIG content= " << sl_sig << " +/- " << sl_sig_err << std::endl;
  std::cout << "SL: SIG/SB ratio = " << r_sl_sigoversb << " +/- " << r_sl_sigoversb_err << std::endl;

  r_sl = r_sl_sigoversb;
  r_sl_err = r_sl_sigoversb_err;
  
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


  ////////////////////////////////////////
  //This part computes the nominal SIG/SB ratio
  ////////////////////////////////////////

  double lv_sb_err=0, lv_sig_err=0;
  double lv_sb = hinteractive->IntegralAndError(lowbin,boundarybin_sb-1,lv_sb_err);
  double lv_sig = hinteractive->IntegralAndError(boundarybin_sig,highbin,lv_sig_err);
  double r_lv_sigoversb = lv_sig/lv_sb;
  double r_lv_sigoversb_err = jmt::errAoverB(lv_sig,lv_sig_err,lv_sb,lv_sb_err);

  std::cout << "LV: SB content= " << lv_sb << " +/- " << lv_sb_err << std::endl;
  std::cout << "LV: SIG content= " << lv_sig << " +/- " << lv_sig_err << std::endl;
  std::cout << "LV: SIG/SB ratio = " << r_lv_sigoversb << " +/- " << r_lv_sigoversb_err << std::endl;

  r_nom = r_lv_sigoversb;
  r_nom_err = r_lv_sigoversb_err;

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
    myLegend->AddEntry(SLplot, "single-top, 1 lepton", "lp");
  }
  myLegend->Draw();
  myC->Print("METshape_logAndRatio_"+sample+"_SLandStandard_"+btagselection+"_"+HTselection+".pdf");

}

//kludgy , but it suffices for now
void makeTTbarWMETPlots() {

  double ratio_sl = 0, ratio_sl_err = 0, ratio_nom = 0, ratio_nom_err = 0;
  std::vector<double> ratios_sl, ratios_sl_err, ratios_nom, ratios_nom_err;
  std::vector<string> selections;

  std::cout << "TTbarJets" << std::endl;

  AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Loose" , "TTbarJets");
  selections.push_back("1BL");
  ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

  AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Tight" , "TTbarJets");
  selections.push_back("1BT");
  ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

  AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Loose" , "TTbarJets");
  selections.push_back("2BL");
  ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

  AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Tight" , "TTbarJets");
  selections.push_back("2BT");
  ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

  AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge3b", "Loose" , "TTbarJets");
  selections.push_back("3B");
  ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

  std::cout << "WJets" << std::endl;

  AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Loose" , "WJets");
  selections.push_back("1BL");
  ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

  AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Tight" , "WJets");
  selections.push_back("1BT");
  ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

  AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Loose" , "WJets");
  selections.push_back("2BL");
  ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

  AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Tight" , "WJets");
  selections.push_back("2BT");
  ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

  AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge3b", "Loose" , "WJets");
  selections.push_back("3B");
  ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

  std::cout << "SingleTop" << std::endl;

  AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Loose" , "SingleTop");
  selections.push_back("1BL");
  ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

  AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge1b", "Tight" , "SingleTop");
  selections.push_back("1BT");
  ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

  AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Loose" , "SingleTop");
  selections.push_back("2BL");
  ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

  AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge2b", "Tight" , "SingleTop");
  selections.push_back("2BT");
  ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);

  AN2011_ttbarw(ratio_sl, ratio_sl_err, ratio_nom, ratio_nom_err,"ge3b", "Loose" , "SingleTop");
  selections.push_back("3B");
  ratios_sl.push_back(ratio_sl);  ratios_sl_err.push_back(ratio_sl_err);  ratios_nom.push_back(ratio_nom);  ratios_nom_err.push_back(ratio_nom_err);


  for(uint i = 0; i<selections.size(); ++i){
    if(i==0) std::cout << "TTbarJets" << std::endl;
    if(i==5) std::cout << "WJets" << std::endl;
    if(i==10) std::cout << "SingleTop" << std::endl;
    std::cout << " & "<< selections.at(i) <<" & " << ratios_nom.at(i) << "$\\pm$" << ratios_nom_err.at(i) << "\t & " << ratios_sl.at(i) << "$\\pm$" << ratios_sl_err.at(i) << " \\\\" << std::endl; 
  }


}



void AN2011_r() { //for paper, data+MC plots at low MET
  TString btagselection="ge1b";
  TCut btagcut = "nbjetsCSVM>=1";
  loadSamples();

  // == draw r(MET)
  clearSamples();
  addSample("PythiaPUQCD");
  setColorScheme("nostack");

  setStackMode(false);
  doOverflowAddition(true);

  leg_x1 = 0.2380496; leg_x2=0.481167; leg_y1=0.5029892; leg_y2=0.9231433;

  //now try drawing it in data
  //careful -- there is a bug in the drawR() implementation
  //it doesn't work right unless you turn on doRatioPlot
  drawTotalSM_=true; 
  usePUweight_=true;
  doRatioPlot(false);
  drawLegend(true);
  setPlotMinimum(0); setPlotMaximum(0.5);
  lumiScale_ = 30.15471; //customize lumiScale_
  resetSamples();
  doData(true);
  TCut util = "pass_utilityHLT==1 && weight<2"; //don't use weight<1 when plotting data!
  selection_ =TCut("HT>=400 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjetsCSVM==0 &&passCleaning==1")&&util;
  //drawR("minDeltaPhiN",4,10,50,150,"Data_eq0b");
  drawR("minDeltaPhiN",4,15,0,150,"Data_eq0b");
  
  selection_ =TCut("HT>=500 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjetsCSVM==0 &&passCleaning==1")&&util;
  drawR("minDeltaPhiN",4,15,0,150,"Data_HT500_eq0b"); 

  selection_ =TCut("HT>=600 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjetsCSVM==0 &&passCleaning==1")&&util;
  drawR("minDeltaPhiN",4,15,0,150,"Data_HT600_eq0b"); 

}


void AN2011_r_SLreq() { //for paper, MC only plots all the way to high MET
  //used for PAS
  loadSamples();
  
  // == draw r(MET)
  clearSamples();
  addSample("PythiaPUQCD");
  setColorScheme("nostack");
  
  usePUweight_=true;

  drawLegend(false);
  doRatioPlot(false);
  doData(false);
  setStackMode(false);
  doOverflowAddition(true);
  
  setPlotMaximum(1); setPlotMinimum(0);
 
  const int nvarbins=15;
  // const float varbins[]={0,15,30,45,80,110,140,170,200,230,260,290,320,350}; //PAS binning //13 bins
   //const float varbins[]={0,15,30,45,60,90,120,150,180,210,240,270,300,400,500,600}; //15 bins
  //  const float varbins[]={0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,450,500,600}; //17 bins
  const float varbins[]={0,20,40,60,80,100,120,140,160,180,200,250,300,400,500,600}; //15 bins
  //  const float varbins[]={0,10,20,30,40,50,60,70,100,130,160,200,240,280,320,360,400,450,500,550,600}; //for pub //20 bins

  //  const float varbins[]={0,20,40,60,80,100,120,140,160,180,200,225,250,275,300,350,450,550,650}; //18 bins

  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjets==0&&weight<1 &&passCleaning==1";
  //   selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 &&weight<1000 &&passCleaning==1";
  //  btagSFweight_="prob0";
  drawR("minDeltaPhiN",4, "MET", nvarbins, varbins,"eq0b");
  
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && nbjets>=1&&weight<1 &&passCleaning==1";
  drawR("minDeltaPhiN",4, "MET", nvarbins, varbins,"ge1b");
  
  setPlotMaximum(3); setPlotMinimum(0);
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

//for the actual AN and paper, use mode 3, logy and doRatio false
void AN2011( TString btagselection="ge1b",const int mode=1, bool logy=false, bool doRatio=false ) {
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
    if ( btagselection=="eq0b") {
      btagcut="1";
      btagSFweight_="prob0";
    }
    else if ( btagselection=="ge1b") {
      btagcut="1";
      btagSFweight_="probge1";
    }
    else  if ( btagselection=="ge2b" ) {
      btagcut = "1";
      btagSFweight_="probge2";
    }
    else if ( btagselection=="ge3b" ) {
      btagcut = "1";
      btagSFweight_="probge3";
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

  /*
  var="minDeltaPhiN"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 40; low=0; high=40;
  //no delta phi cut, loose MET window (only good for MC)
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100")&&btagcut;
  drawPlots(var,nbins,low,high,xtitle,"Arbitrary units", "minDeltaPhiN_looseMET_MConly_"+btagselection+modestring);
  */
  
  //  return;
  
  // ========= regular N-1 plots

  resetSamples(); //use all samples
  setStackMode(true); //regular stack
  setColorScheme("stack");
  doData(true);
  drawMCErrors_=true;
  
  //ben's experience with the Summer11 analysis is that this mod was necessary for style reasons
  //since these variables are global, we can just set them here.
  //they are reset down at the bottom with resetLegendPosition()
  //  leg_x1 = 0.6; leg_x2=0.96; leg_y1=0.4; leg_y2=0.88;//bensep28 - for data/mc stack comparison
  leg_x2=0.96; leg_y1=0.4; //leg_y2=0.88;//bensep28 - for data/mc stack comparison
  //leg_x1 = 0.6; leg_x2=0.96; leg_y1=0.42; leg_y2=0.9;//bensep28 -for data/mc stack comparison NJETS

  
  var="minDeltaPhiN"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 20; low=0; high=40;
  //no delta phi cut
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 &&passCleaning==1")&&btagcut;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_minDeltaPhiN_"+btagselection+modestring);
  
  // 2BT but with loose MET
  selection_ =TCut("HT>=600 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 &&passCleaning==1")&&btagcut;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_minDeltaPhiN_HT600_"+btagselection+modestring);    

  //signal region N-1 minDeltaPhiN
  if (btagselection=="ge1b") {
    selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=250 &&passCleaning==1")&&btagcut;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_minDeltaPhiN_1BL_"+btagselection+modestring);

    nbins=10;
    selection_ =TCut("HT>=500 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=500 &&passCleaning==1")&&btagcut;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_minDeltaPhiN_1BT_"+btagselection+modestring);
  }
  else if (btagselection=="ge2b") {
    selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=250 &&passCleaning==1")&&btagcut;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_minDeltaPhiN_2BL_"+btagselection+modestring);

    nbins=10;
    selection_ =TCut("HT>=600 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=300 &&passCleaning==1")&&btagcut;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_minDeltaPhiN_2BT_"+btagselection+modestring);
  }
  else if (btagselection=="ge3b") {
    nbins=10;
    selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=250 &&passCleaning==1")&&btagcut;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_minDeltaPhiN_3B_"+btagselection+modestring);

  }

  // n Jets
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 7; low=2; high=9;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_njets_"+btagselection+modestring);
  
  // n Nets - QCD dominated (LDP)
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN < 4 &&passCleaning==1")&&btagcut;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 7; low=2; high=9;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_njetsLDP_"+btagselection+modestring);
  
  //HT
  selection_ =TCut("cutHT==1&&cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 20; low=400; high=1100;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_HT_"+btagselection+modestring,0,"GeV");
  
  //MET
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 35; low=150; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_"+btagselection+modestring,0,"GeV");
  
  //MET distribution with tighter HT cut
  selection_ =TCut("HT>500 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  TString thisytitle = "";
  if(btagselection=="ge2b") { nbins = 7; low=200; high=500;} 
  else{ nbins = 17; low=200; high=500;}
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_HT500_"+btagselection+modestring,0,"GeV");
  
  selection_ =TCut("HT>600 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  if(btagselection=="ge2b") { nbins = 7; low=150; high=500; } 
  else{ nbins = 17; low=200; high=500; }
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_HT600_"+btagselection+modestring,0,"GeV");

    // == finally, draw the signal region only!
  if (btagselection=="ge1b") {
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=250 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 20; low=400; high=1100;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_HT_1BL_"+btagselection+modestring,0,"GeV");

  var="MET+HT"; xtitle="M_{eff} [GeV]";
  nbins = 20; low=650; high=1550;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_Meff_1BL_"+btagselection+modestring,0,"GeV");

  selection_ =TCut("HT>=500 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=500 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="HT"; xtitle="H_{T} (GeV)";
  nbins = 20; low=500; high=1050;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_HT_1BT_"+btagselection+modestring,0,"GeV");

  var="MET+HT"; xtitle="M_{eff} [GeV]";
  nbins = 10; low=1000; high=2000;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_Meff_1BT_"+btagselection+modestring,0,"GeV");
  }
  else if (btagselection=="ge2b") {
    selection_ =TCut("HT>400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=250 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
    var="HT"; xtitle="H_{T} (GeV)";
    nbins = 20; low=400; high=1100;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_HT_2BL_"+btagselection+modestring,0,"GeV");

    var="MET+HT"; xtitle="M_{eff} [GeV]";
    nbins = 20; low=650; high=1550;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_Meff_2BL_"+btagselection+modestring,0,"GeV");

    selection_ =TCut("HT>600 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=300 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
    var="HT"; xtitle="H_{T} (GeV)";
    nbins = 10; low=600; high=1600;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_HT_2BT_"+btagselection+modestring,0,"GeV");

    var="MET+HT"; xtitle="M_{eff} [GeV]";
    nbins = 11; low=900; high=2000;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_Meff_2BT_"+btagselection+modestring,0,"GeV");

  }
  else if (btagselection=="ge3b") {
    selection_ =TCut("HT>400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=250 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
    var="HT"; xtitle="H_{T} (GeV)";
    nbins = 5; low=400; high=1200;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_HT_3B_"+btagselection+modestring,0,"GeV");

    var="MET+HT"; xtitle="M_{eff} [GeV]";
    nbins = 5; low=650; high=1650;
    drawPlots(var,nbins,low,high,xtitle,"Events", "SIG_Meff_3B_"+btagselection+modestring,0,"GeV");
  }
  // ===== single lepton selection
  
  // == MET for the electron sample
  selection_ =TCut("MET>=200 && HT>=400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && nElectrons==1 && nMuons==0 && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=200; high=600;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_1e0mu_"+btagselection+modestring,0,"GeV");
   
  //pT?
  var="eleet1"; xtitle="electron p_{T} [GeV]";
  nbins = 20; low=0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_eleET_1e0mu_"+btagselection+modestring,0,"GeV");
  
  
  // == MET for the muon sample
  selection_ =TCut("MET>=200 && HT>=400 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && nElectrons==0 && nMuons==1 && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=200; high=600;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_0e1mu_"+btagselection+modestring,0,"GeV");
  
  var="muonpt1"; xtitle="muon p_{T} [GeV]";
  nbins = 20; low=0; high=200;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_muonpT_0e1mu_"+btagselection+modestring,0,"GeV");
  
  
  // == MET for the combined sample
  selection_ =TCut("MET>=200 &&cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=200; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events/20 GeV", "SBandSIG_MET_SL_"+btagselection+modestring);
  /*
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
  */
  // == MET for the combined sample (HT>500)
  selection_ =TCut("MET>=200 && HT>500 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=200; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_SL_HT500_"+btagselection+modestring,0,"GeV");
  
  // == MET for the combined sample (HT>600)
  selection_ =TCut("MET>=200 && HT>600 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100 &&passCleaning==1")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 15; low=200; high=500;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_SL_HT600_"+btagselection+modestring,0,"GeV");
  
  //different scale to compare to Kristen
//   selection_ =TCut("MET>=150 && HT>500 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0)) && minDeltaPhiN >= 4 && MT_Wlep>=0 && MT_Wlep<100")&&btagcut;
//   var="MET"; xtitle="E_{T}^{miss} [GeV]";
//   nbins = 35; low=150; high=500;
//   drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_SL_HT500_morebins_"+btagselection+modestring);

  
/*
  // == take a look at WJets
  //doesn't work with the fancy b tag SF code. And we're not using it anyway
  selection_ ="MET>=200 && cutHT==1 && cutPV==1 && cutTrigger==1 && njets==2 && (nElectrons==0 && nMuons==1) && minDeltaPhiN >= 4 && nbjetsCSVM==1";
  var="MT_Wlep"; xtitle="M_{T} [GeV]";
   nbins = 20; low=0; high=200;
   drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MT_0e1mu_eq1b");
*/
/*
  //MET
  selection_ =TCut("cutHT==1 && cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN < 4")&&btagcut;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  nbins = 30; low=150; high=450;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SBandSIG_MET_ldp_"+btagselection+modestring);
 
  */

  /*
  //SB njet
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && MET<250 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="njets"; xtitle="Jet multiplicity";
  nbins = 7; low=2; high=9;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SB_njets_"+btagselection+modestring);
  //SB HT
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && MET<250 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="HT"; xtitle="H_{T} [GeV]";
  nbins = 20; low=350; high=1050;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SB_HT_"+btagselection+modestring);
  //SB mindphin
  selection_ =TCut("HT>=400 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=200 && MET<250 && minDeltaPhiN >= 4 &&passCleaning==1")&&btagcut;
  var="minDeltaPhiN"; xtitle="#Delta #phi_{N}^{min}";
  nbins = 20; low=0; high=40;
  drawPlots(var,nbins,low,high,xtitle,"Events", "SB_mindphin_"+btagselection+modestring);
  */


  resetLegendPosition();
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
  
  TCut util = "pass_utilityHLT==1 && weight<1000";
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
      selection_ =TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=50 && MET<100 && nbjetsCSVM==0 && pass_utilityHLT==1")&&HTcut;//note no mdp or run cut
    }
    else{
      selection_ =TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=50 && MET<100 && nbjetsCSVM==0 && pass_utilityHLT==1")&&HTcut&&runCut;//note no mdp
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
    
    
    selection_ =TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=50 && MET<100 && nbjetsCSVM==0 && pass_utilityHLT==1 && minDeltaPhiN>=4")&&HTcut&&runCut;
    //drawSimple("nGoodPV",pvnbins,pvlow,pvhigh,histfilename, "nGoodPVPass_prescale_data","data");
    drawSimple("nGoodPV",pvnbins,pvbins,histfilename, "nGoodPVPass_prescale_data","data");
    TH1D* hPVprescalePass = (TH1D*)hinteractive->Clone("hPVprescalePass");
    
    selection_ =TCut("cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=50 && MET<100 && nbjetsCSVM==0 && pass_utilityHLT==1 && minDeltaPhiN<4")&&HTcut&&runCut;
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
