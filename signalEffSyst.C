/* 
this is a fork of drawReducedTrees.C
It shares the header drawReducedTrees.h


gSystem->Load("TSelectorMultiDraw_C.so");
gSystem->Load("CrossSectionTable_cxx.so");
.L signalEffSyst.C++

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
TString inputPath = "/cu3/joshmt/reducedTrees/V00-02-35l/"; 
TString dataInputPath =  "dummy"; //no data needed!
TString inputPathTTbar = "dummy";

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

  //take care of trigger efficiency systematics, which differs between SB and SIG, SL, LDP, etc
  if (region.isSIG && !isSL &&!isLDP) results.set("trigger",eff_SIG_MHT_err_[0],eff_SIG_MHT_err_[1]);
  else if (!region.isSIG && !isSL && !isLDP) results.set("trigger",eff_SB_MHT_err_[0],eff_SB_MHT_err_[1]);
  else if (region.isSIG && isSL) results.set("trigger",eff_SIG_SL_MHT_err_[0],eff_SIG_SL_MHT_err_[1]);
  else if (!region.isSIG && isSL) results.set("trigger",eff_SB_1m_MHT_err_[0],eff_SB_1m_MHT_err_[1]);
  else if (region.isSIG && isLDP) results.set("trigger",eff_SIG_ldp_MHT_err_[0],eff_SIG_ldp_MHT_err_[1]);
  else if (!region.isSIG && isLDP) results.set("trigger",eff_SB_ldp_MHT_err_[0],eff_SB_ldp_MHT_err_[1]);
  else assert(0);

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
  else if (btagselection=="eq1b") {
    bcut="nbjets==1";
    bweightstring="prob1";  
  }
  else if (btagselection=="eq2b") {
    bcut="nbjets==2";
    bweightstring="prob2";
  }
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
  useHTeff_=false;
  useMHTeff_=false;
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

  useMHTeff_=true;
  useHTeff_=true;
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

  // this code was for special checks
  //now we can test the effect of the changing the charm and LF tag rates
  //all of this is done with the regular JERbias file
//   TString btagprobbase=btagSFweight_;

//   //b +
//   btagSFweight_ = btagprobbase + "_bplus";
//   drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
//   double bplus=  getIntegral(sampleOfInterest);

//   //b -
//   btagSFweight_ = btagprobbase + "_bminus";
//   drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
//   double bminus=  getIntegral(sampleOfInterest);

//   //charm +
//   btagSFweight_ = btagprobbase + "_cplus";
//   drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
//   double cplus=  getIntegral(sampleOfInterest);

//   //charm -
//   btagSFweight_ = btagprobbase + "_cminus";
//   drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
//   double cminus=  getIntegral(sampleOfInterest);

//   //LF +
//   btagSFweight_ = btagprobbase + "_lplus";
//   drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
//   double lplus=  getIntegral(sampleOfInterest);

//   //LF -
//   btagSFweight_ = btagprobbase + "_lminus";
//   drawPlots(var,nbins,low,high,xtitle,"events","dummyH");
//   double lminus=  getIntegral(sampleOfInterest);

//   results.eff_derivative_b = ( ((1.0/0.1)* (bplus - nominal) / nominal) + ((-1.0/0.1)* (bminus - nominal) / nominal))*0.5;
//   results.eff_derivative_c = ( ((1.0/0.1)* (cplus - nominal) / nominal) + ((-1.0/0.1)* (cminus - nominal) / nominal))*0.5;
//   results.eff_derivative_l = ( ((1.0/0.1)* (lplus - nominal) / nominal) + ((-1.0/0.1)* (lminus - nominal) / nominal))*0.5;

//   //reset the sf weight string
//   btagSFweight_ = bweightstring;

//   //for LM9 the histogram should have exactly 1 bin
//   TH1D* btageff_avg_nominal = (TH1D*) files_[currentConfig_][sampleOfInterest]->Get("btageff_avg");
//   double av_btageff_nominal = btageff_avg_nominal->GetBinContent(1);
//   double av_btageff_var1=0;
//   double av_btageff_var2=0;

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
    if (j%2==0) { //this one runs first; by convention this one should be "down"
      var2=(thisn-nominal)/nominal;
//       if (configDescriptions_.getVariedSubstring(currentConfig_).Contains("BTag")) {
// 	TH1D* btageff_avg_2 = (TH1D*) files_[currentConfig_][sampleOfInterest]->Get("btageff_avg");
// 	av_btageff_var2=btageff_avg_2->GetBinContent(1);
// 	bminus=thisn;
//       }
    }
    else { //this one runs second; by convention this one should be "up"
      var1=(thisn-nominal)/nominal;
      results.set(configDescriptions_.getVariedSubstring(currentConfig_), var2,var1);
//       if (configDescriptions_.getVariedSubstring(currentConfig_).Contains("BTag")) { //capitalization ok?
// 	TH1D* btageff_avg_1 = (TH1D*) files_[currentConfig_][sampleOfInterest]->Get("btageff_avg");
// 	av_btageff_var1=btageff_avg_1->GetBinContent(1);
// 	bplus=thisn;
//       }
    }
  }
  
  //reset to the nominal with JER bias
  currentConfig_=configDescriptions_.getCorrected();

  //var2 is down, var1 is up
//   double delta_btageff_plus =  (av_btageff_var1 - av_btageff_nominal) / av_btageff_nominal;
//   double delta_btageff_minus = (av_btageff_var2 - av_btageff_nominal) / av_btageff_nominal;
//   //this is the error on the <btageff> for +/- 1 sigma variations on the SF
//   results.sigma_btageff = 0.5*(fabs(delta_btageff_plus) + fabs(delta_btageff_minus));
//   results.eff_derivative_b_1s =  0.5* ( (((bplus - nominal) / nominal)/delta_btageff_plus ) 
// 					+ (((bminus - nominal) / nominal)/delta_btageff_minus)
// 					);

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
    if (sampleOfInterest!="T1bbbb"&&sampleOfInterest!="T2bb"&&sampleOfInterest!="T2tt"&&sampleOfInterest!="T1tttt")    thetree->Project("pdfEventCounter","HT",getCutString(false, "",extraWeight,i,"CTEQ").Data());
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
    if (sampleOfInterest!="T1bbbb"&&sampleOfInterest!="T2bb"&&sampleOfInterest!="T2tt"&&sampleOfInterest!="T1tttt")    thetree->Project("pdfEventCounter","HT",getCutString(false, "",extraWeight,i,"MSTW").Data());
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
      if (sampleOfInterest!="T1bbbb"&&sampleOfInterest!="T2bb"&&sampleOfInterest!="T2tt"&&sampleOfInterest!="T1tttt")    thetree->Project("pdfEventCounter","HT",getCutString(false, "",extraWeight,i,"NNPDF").Data());
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
    (*textfiles_human[i])<<" \\hline                   &"<<endl;

//     cout<<" === cross-check "<<endl;
//     cout<<"SB DeltaEff / DeltaBtagEff = "<<SB.eff_derivative_b_1s<<endl;
//     cout<<"SIG DeltaEff / DeltaBtagEff = "<<SIG.eff_derivative_b_1s<<endl;

//     cout<<"SBSL DeltaEff / DeltaBtagEff = "<<SBSL.eff_derivative_b_1s<<endl;
//     cout<<"SIGSL DeltaEff / DeltaBtagEff = "<<SIGSL.eff_derivative_b_1s<<endl;

//     cout<<"SBLDP DeltaEff / DeltaBtagEff = "<<SBLDP.eff_derivative_b_1s<<endl;
//     cout<<"SIGLDP DeltaEff / DeltaBtagEff = "<<SIGLDP.eff_derivative_b_1s<<endl;
//     cout<<" === new "<<endl;
//     cout<<"SB DeltaEff / DeltaBtagEff = "<<SB.eff_derivative_b<<endl;
//     cout<<"SIG DeltaEff / DeltaBtagEff = "<<SIG.eff_derivative_b<<endl;

//     cout<<"SBSL DeltaEff / DeltaBtagEff = "<<SBSL.eff_derivative_b<<endl;
//     cout<<"SIGSL DeltaEff / DeltaBtagEff = "<<SIGSL.eff_derivative_b<<endl;

//     cout<<"SBLDP DeltaEff / DeltaBtagEff = "<<SBLDP.eff_derivative_b<<endl;
//     cout<<"SIGLDP DeltaEff / DeltaBtagEff = "<<SIGLDP.eff_derivative_b<<endl;
//     cout<<" === charm "<<endl;
//     cout<<"SB DeltaEff / DeltaCtagEff = "<<SB.eff_derivative_c<<endl;
//     cout<<"SIG DeltaEff / DeltaCtagEff = "<<SIG.eff_derivative_c<<endl;

//     cout<<"SBSL DeltaEff / DeltaCtagEff = "<<SBSL.eff_derivative_c<<endl;
//     cout<<"SIGSL DeltaEff / DeltaCtagEff = "<<SIGSL.eff_derivative_c<<endl;

//     cout<<"SBLDP DeltaEff / DeltaCtagEff = "<<SBLDP.eff_derivative_c<<endl;
//     cout<<"SIGLDP DeltaEff / DeltaCtagEff = "<<SIGLDP.eff_derivative_c<<endl;
//     cout<<" === LF "<<endl;
//     cout<<"SB DeltaEff / DeltaLFtagEff = "<<SB.eff_derivative_l<<endl;
//     cout<<"SIG DeltaEff / DeltaLFtagEff = "<<SIG.eff_derivative_l<<endl;

//     cout<<"SBSL DeltaEff / DeltaLFtagEff = "<<SBSL.eff_derivative_l<<endl;
//     cout<<"SIGSL DeltaEff / DeltaLFtagEff = "<<SIGSL.eff_derivative_l<<endl;

//     cout<<"SBLDP DeltaEff / DeltaLFtagEff = "<<SBLDP.eff_derivative_l<<endl;
//     cout<<"SIGLDP DeltaEff / DeltaLFtagEff = "<<SIGLDP.eff_derivative_l<<endl;

  }

  for (unsigned int i=0; i<searchRegions_.size(); i++) {
    textfiles.at(i)->close();
    textfiles_human.at(i)->close();

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
  loadSusyScanHistograms();


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
  useMHTeff_=false;
  useHTeff_=false;
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
  useMHTeff_=true;
  useHTeff_=true;
  selection_ = baseline && HTcut && passOther && SRMET; //remove b tag
  btagSFweight_=bweightstring; // add b tag back in the form of the SF
  susyScanYields nominalYields=getSusyScanYields(sampleOfInterest);

  //only applicable to mSugra
  const bool do_kfactor_variation=false;  //for now...*disable* the cross-section uncertainty
  if (ismSugra && !do_kfactor_variation)  cout<<" NLO cross-section uncertainty disabled!"<<endl;
  //NLO k factor uncertainty -- vary it up
  susyCrossSectionVariation_="Plus";
  susyScanYields kFactorPlusYields= ismSugra && do_kfactor_variation ? getSusyScanYields(sampleOfInterest) : nominalYields;
  
  //NLO k factor uncertainty -- vary it down
  susyCrossSectionVariation_="Minus";
  susyScanYields kFactorMinusYields=ismSugra && do_kfactor_variation ? getSusyScanYields(sampleOfInterest): nominalYields;

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

    assert( eff_SB_1e_MHT_err_[0] == eff_SB_1m_MHT_err_[0]); //until i update this code to work with the split e/mu
    assert( eff_SB_1e_MHT_err_[1] == eff_SB_1m_MHT_err_[1]);

    //take care of trigger efficiency systematics, which differs between SB and SIG, SL, LDP, etc
    if (region.isSIG && !isSL &&!isLDP) theseResults.set("trigger",eff_SIG_MHT_err_[0],eff_SIG_MHT_err_[1]);
    else if (!region.isSIG && !isSL && !isLDP) theseResults.set("trigger",eff_SB_MHT_err_[0],eff_SB_MHT_err_[1]);
    else if (region.isSIG && isSL) theseResults.set("trigger",eff_SIG_SL_MHT_err_[0],eff_SIG_SL_MHT_err_[1]);
    else if (!region.isSIG && isSL) theseResults.set("trigger",eff_SB_1m_MHT_err_[0],eff_SB_1m_MHT_err_[1]);
    else if (region.isSIG && isLDP) theseResults.set("trigger",eff_SIG_ldp_MHT_err_[0],eff_SIG_ldp_MHT_err_[1]);
    else if (!region.isSIG && isLDP) theseResults.set("trigger",eff_SB_ldp_MHT_err_[0],eff_SB_ldp_MHT_err_[1]);
    else assert(0);

    if (rawYields[imsugra->first] >0) {
      cout<<imsugra->first.first<<" "<<imsugra->first.second<<" " //m0 and m12
 	  <<rawYields[imsugra->first]<<" "<<nominalYields[imsugra->first]<<" "
 	  <<kFactorPlusYields[imsugra->first]<<" "<<kFactorMinusYields[imsugra->first]<<" ";

      theseResults.rawYield = rawYields[imsugra->first];

      double nominalYield = nominalYields[imsugra->first];
      theseResults.effCorr = nominalYield/theseResults.rawYield;
      if (nominalYield==0) nominalYield=0.000001; //should make this better

      //remove fabs()
      double kfactorP = ((kFactorPlusYields[imsugra->first] - nominalYield) / nominalYield);
      double kfactorM = ((kFactorMinusYields[imsugra->first] - nominalYield) / nominalYield);
      theseResults.set("kFactor",kfactorM,kfactorP);

      //pdf uncertainties
      //      if (fullPdfUncertainties) {
	
	double cteqXminusMax = cteqXplusminus.first[imsugra->first];
	double cteqXplusMax = cteqXplusminus.second[imsugra->first];

	//cout<<"\nCTEQ nominal xMinusMax xPlusMax "<<nominalYield<<" "<<cteqXminusMax<<" "<<cteqXplusMax<<endl;
	double  cteqFracPDF = 0.5* ( cteqXminusMax + cteqXplusMax);
	cteqFracPDF = fabs( cteqFracPDF ) /nominalYield;
	
	cout<<endl;
	cout<<imsugra->first.first<<" "<<imsugra->first.second<<" CTEQ  % = "<<100*cteqFracPDF<<endl;

	double mstwXminusMax = mstwXplusminus.first[imsugra->first];
	double mstwXplusMax = mstwXplusminus.second[imsugra->first];
	double  mstwFracPDF = 0.5* ( mstwXminusMax + mstwXplusMax);
	mstwFracPDF = fabs( mstwFracPDF ) /nominalYield;
	cout<<imsugra->first.first<<" "<<imsugra->first.second<<" MSTW  % = "<<100*mstwFracPDF<<endl;
	if (cteqFracPDF > mstwFracPDF) mstwFracPDF = cteqFracPDF;
	
	double nnpdfTotalX  = nnpdfXX2.first[imsugra->first];
	double nnpdfTotalX2 = nnpdfXX2.second[imsugra->first];
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
	cout<<ivariation->second.first[imsugra->first]<<" "<<ivariation->second.second[imsugra->first]<<" ";
	
	//get rid of fabs()
	double percent1 = ((ivariation->second.first[imsugra->first] - nominalYield)/nominalYield);
	double percent2 = ((ivariation->second.second[imsugra->first] - nominalYield)/nominalYield);
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

  if (!sampleOfInterest.Contains("SUGRA")) loadScanSMSngen(sampleOfInterest);
  loadSusyScanHistograms();
  
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
  if (!sampleOfInterest.Contains("SUGRA")) loadScanSMSngen(sampleOfInterest);
  loadSusyScanHistograms();

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
    nbinsx=300;
    nbinsy=100;
    lowx=0;
    lowy=0;
    highx=3000;
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
      scanProcessTotalsMap[iscanpoint->first]["CTEQ"]->Integral(0,10,0,0) : //this is supposed to integrate over the subprocesses but not the pdf indices
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

/* the version that runs all 6 boxes serially
void runSystematics2011_scan(TString sampleOfInterest, const unsigned int i) {
  //now converting this to work also on mSugra "mSUGRAtanb40" as well as simplified models e.g. T1bbbb
  
  loadSamples();
  clearSamples();
  addSample(sampleOfInterest);
  
  setSearchRegions();

  loadSusyScanHistograms();
  
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
*/
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
    useHTeff_=true;
    useMHTeff_=true;
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
    int  nbinsx=310; float lowx=-0.5; float highx=3100-0.5; //UPDATED for new scans
    
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
      for (map<pair<int,int>, map<TString,TH2D*> >::iterator iscanpoint = scanProcessTotalsMap.begin(); iscanpoint!=scanProcessTotalsMap.end(); ++iscanpoint) {
	int m0=iscanpoint->first.first;
	int m12=iscanpoint->first.second;
	
	TH2D* thishist = scanProcessTotalsMap[make_pair(m0,m12)]["CTEQ"];
	int thisn = TMath::Nint(thishist->GetBinContent(i,0)); // n gen
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
    for (map<pair<int,int>, map<TString,TH2D*> >::iterator iscanpoint = scanProcessTotalsMap.begin(); iscanpoint!=scanProcessTotalsMap.end(); ++iscanpoint) {
      int nentries=    iscanpoint->second["CTEQ"]->Integral(0,10,0,0); //integrate over subprocesses but not pdf indices
      if (nentries >10000 || nentries<9500) continue;
      
      double Nraw = 0;
      
      int bin=  raw0[0]->FindBin(iscanpoint->first.first , iscanpoint->first.second);
      for (unsigned int i=0; i<raw0.size(); i++) {
	Nraw += raw0[i]->GetBinContent(bin);
      }
      //the badly named Nraw is sum_i Npass_i * sigma_i * lumi / Ngen_i = sum_i eff_i * sgima_i * lumi = sum_i Nobs_i = N observed
      
      //total Sigma is already the sum over i
      double crossSectionTimesLumi = totalSigma.GetBinContent(bin) * lumiScale_;
      if (crossSectionTimesLumi >0)    theEff[iscanpoint->first] = Nraw / crossSectionTimesLumi;//make_pair(Nraw / crossSectionTimesLumi,0);
      else if (Nraw < 0.00001 && crossSectionTimesLumi <0.00001) {}//do nothing
      else {cout<<"Something weird? crossSectionTimeLumi is zero!"<<endl;}
    }
    
    renewCanvas();
    thecanvas->cd()->SetRightMargin(0.2);
    if (h2d!=0) delete h2d;

    //redo the binning with what appears to be the actual binning (for better presentation)
    assert(0); //need to update binning for new scans
    nbinsx=105;  lowx=-0.5; highx=2100-0.5; //m0
    nbinsy=40;  lowy=-0.5;  highy=800-0.5; //m12

    TString hname="mSugraEffMap";
    h2d = new TH2D(hname,hname,nbinsx,lowx,highx,nbinsy,lowy,highy);
    TString opt="colz";
    for (susyScanYields::iterator ieff = theEff.begin(); ieff!=theEff.end(); ++ieff) {
      int bin=    h2d->FindBin( ieff->first.first, ieff->first.second);
      h2d->SetBinContent(bin, ieff->second);
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
