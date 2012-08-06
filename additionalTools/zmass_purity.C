//////////////////////////////////////////////////////////////////////////
// 
// Simple fitting routine to extract the purity of the "VL" sample
// used for the 2011 Z->vv estimate
// 
// Global variables used to set configuration parameters.
//
// To make the input histograms/trees needed to do the fits, one needs
// to run:
//  -> drawZllPurity() in drawReducedTrees.C  (for the binned histos)
//  -> zllPuritySkim.C  (for trees)
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooWorkspace.h"
#include "RooAddPdf.h"
#include "RooPlot.h"

#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLatex.h"


using namespace RooFit ;

// *** Set the desired signal and background pdf's **** //
enum signalMode {kBWxCB=0, kCB, kVxCB, kDoubleGauss};
signalMode theSignalMode_ = kBWxCB;
//signalMode theSignalMode_ = kVxCB;
//signalMode theSignalMode_ = kDoubleGauss;
enum bgMode {kPoly1=0, kPoly2, kPoly3, kCheby1, kCheby2, kCheby3};
bgMode theBGMode = kPoly1;


// *** SB or SIG **** //
bool isSIGregion = true;

// *** Mu or Electron **** //
bool mumode = true;


void drawPlotHeader(double xoff=0.0){
  //misc plot header
  TLatex* text1 = new TLatex(3.570061,23.08044,"CMS");
  text1->SetX(0.68 +0.2);
  text1->SetNDC();
  text1->SetTextAlign(13);
  text1->SetY(0.97 + 0.007);
  text1->SetTextFont(42);
  text1->SetTextSizePixels(24);
  text1->SetTextSize(0.045);
  text1->Draw();

                                                                                                                                                  
  double lumiScale_ = 4982;//final pixel-based 2011 lumi
  TString astring; astring.Form("L_{int} = %.2f fb^{-1}, #sqrt{s} = 7 TeV",lumiScale_/1000.);
  TLatex* text2 = new TLatex(3.570061,23.08044,astring);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(xoff+0.17);
  text2->SetY(0.97+0.015);
  text2->SetTextFont(42);
  text2->SetTextSizePixels(24);
  text2->SetTextSize(0.045); //copied from ben's code. maybe needs a switch so it is only used for AN2011()                                                                                                  
  text2->Draw();

}



void zmass_purity()
{

  gStyle->SetOptStat(111111);


  //////////////////////////////////////////
  //put all the pdf stuff in a workspace
  //////////////////////////////////////////

  RooWorkspace* w = new RooWorkspace("w");

  //  //w->factory("BreitWigner::bwshape(m[60,120],bw_mean[91.2,89,93],bw_width[2.49,2,3])");
  w->factory("BreitWigner::bwshape(m[60,120],bw_mean[91.2,89,93],bw_width[2.7,2,5])");
  //  w->factory("CBShape::cbshape(m[60,120],cb_m0[0.0,-5,1.0],cb_sigma[1.8,0.0,15],alpha[0.959,0.3,2],n[1.266,0.5,10])");
  
  //  w->factory("BreitWigner::bwshape(m[60,120],bw_mean[91.1,89,93],bw_width[3.0,2,5])");//testing for SB electron (will this remove the non-pos-def error?)
  w->factory("CBShape::cbshape(m[60,120],cb_m0[0.0,-5,1.0],cb_sigma[1.8,0.0,5],alpha[1.4,0.3,2],n[2.1,0.5,10])");
  w->factory("FCONV::bw_x_cb(m,bwshape,cbshape)");

  w->factory("Voigtian::voigtshape(m[60,120],v_mean[91.2,89,93],v_width[2.5,2,5],v_sigma[1,0.1,10])");
  w->factory("FCONV::v_x_cb(m,voigtshape,cbshape)");

  w->factory("Gaussian::sig1(m[60,120],g1_mean[91,87,93],g1_sigma[2.5,0,5])") ;
  w->factory("Gaussian::sig2(m[60,120],g1_mean,g2_sigma[9,0,20])") ;
  w->factory("SUM::doublegauss(g_fract[0.9,0.5,1]*sig1,sig2)");

  w->factory("Polynomial::polyshape(m[60,120],{poly_c1[0.0,-1,1],poly_c2[0.1,-1,1],poly_c3[0.1,-1,1]})"); // for fitting cubic background

  //w->factory("Chebychev::chebyshape(m[60,120],{cheby_c1[0.1,-1.5,1],cheby_c2[0.1,-1,1],cheby_c3[0.1,-1,1]})");
  w->factory("Chebychev::chebyshape(m[60,120],{cheby_c1[0.0,-1.5,1],cheby_c2[0.0,-1,1],cheby_c3[0.0,-1,1]})");
  //w->factory("Polynomial::polyshape(m[60,120],{poly_c1[0.05,-1,1],poly_c2[-0.0005,-1,1]})");
  //w->factory("Exponential::expshape(m[60,120],lambda[1,0,20])");//doesn't converge well

 
  //for Breit-Wigner (+) Crystal-Ball for signal, polynomial for bg
  if(theBGMode==kPoly1 || theBGMode==kPoly2 || theBGMode==kPoly3){
    if(theSignalMode_==kBWxCB)
      w->factory("SUM::sum(fract[0.85,0.75,1] * bw_x_cb, polyshape)");
    //w->factory("SUM::sum(fract[0.90,0,1] * bw_x_cb, expshape)");
    else if(theSignalMode_ == kCB)
      w->factory("SUM::sum(fract[0.90,0,1] * cbshape2, polyshape)");
    else if(theSignalMode_ == kVxCB)
      w->factory("SUM::sum(fract[0.85,0.75,1] * v_x_cb, polyshape)");
    else if(theSignalMode_ == kDoubleGauss)
      w->factory("SUM::sum(fract[0.85,0.75,1] * doublegauss, polyshape)");
  }
  else if(theBGMode==kCheby1 || theBGMode==kCheby2 || theBGMode==kCheby3){
    if(theSignalMode_==kBWxCB)
      w->factory("SUM::sum(fract[0.85,0.75,1] * bw_x_cb, chebyshape)");
    //w->factory("SUM::sum(fract[0.90,0,1] * bw_x_cb, expshape)");
    else if(theSignalMode_ == kCB)
      w->factory("SUM::sum(fract[0.90,0,1] * cbshape2, chebyshape)");
    else if(theSignalMode_ == kVxCB)
      w->factory("SUM::sum(fract[0.85,0.75,1] * v_x_cb, chebyshape)");
    else if(theSignalMode_ == kDoubleGauss)
      w->factory("SUM::sum(fract[0.85,0.75,1] * doublegauss, chebyshape)");
  }


  // Print workspace contents
  w->Print() ;

    
  //Set the relevant terms to 0
  if(theBGMode!=kPoly3){ w->var("poly_c3")->setVal(0);  w->var("poly_c3")->setConstant();}
  if(theBGMode!=kPoly2){ w->var("poly_c2")->setVal(0);  w->var("poly_c2")->setConstant();}

  if(theBGMode!=kCheby3){ w->var("cheby_c3")->setVal(0);  w->var("cheby_c3")->setConstant();}
  if(theBGMode!=kCheby2){ w->var("cheby_c2")->setVal(0);  w->var("cheby_c2")->setConstant();}


  ////////////////////////////////
  //Import histograms/trees
  ////////////////////////////////

  //HISTOGRAM FORMAT
  TFile * f_zmass;
  if(mumode) f_zmass = new TFile("zmass_purity_mu.root");
  else f_zmass = new TFile("zmass_purity_ele.root");

  TH1* h_zmass;
  if(mumode) h_zmass= (TH1 *)f_zmass->Get("zmumumass_loose_doublemudata");
  else h_zmass = (TH1 *)f_zmass->Get("zeemass_loose_doubleeledata");  

  // all-inclusive sample (to get signal shape parameters)
  RooRealVar* myvar = w->var("m");
  RooDataHist dh("dh","dh",*myvar,Import(*h_zmass));


  // load data points in VL sample in an unbinned way
  // apparently it reduces the error in the mumu case: 20% -> 10% 
  TFile * f_zmass_VL;
  if(mumode) f_zmass_VL = new TFile("puritytree_zmumu.root");
  else f_zmass_VL = new TFile("puritytree_zee.root");
  TTree* ztree_vl;
  if(isSIGregion)
    ztree_vl = (TTree*) gDirectory->Get("purTreeMVL");  
  else
    ztree_vl = (TTree*) gDirectory->Get("purTreeMVL_SB");


  // Make plot of all-incluseve (binned) sample
  // Also fit this sample
  RooPlot* mframe2 = myvar->frame(Title("Imported TH1 with Poisson error bars")) ;
  if(mumode)  mframe2->SetXTitle("M(#mu#mu) [GeV]");
  else mframe2->SetXTitle("M(ee) [GeV]");
  mframe2->GetYaxis()->SetTitle("Events / 1.5 GeV"); //hard-coded for now...
  dh.plotOn(mframe2) ; 

  RooAddPdf* mysumshape = (RooAddPdf*) w->pdf("sum");
  mysumshape->fitTo(dh); 
  //mysumshape->fitTo(dh,PrintLevel(-1),Extended());
  mysumshape->plotOn(mframe2);


  //show the fit parameters on the plot
  //mysumshape->paramOn(mframe2,Layout(0.65,0.99));
  //mframe2->getAttText()->SetTextSize(0.02); 


  //overlay the BG fit
  RooPolynomial* mypolyshape = (RooPolynomial*) w->pdf("polyshape");
  RooChebychev* mychebyshape = (RooChebychev*) w->pdf("chebyshape");
  if(theBGMode==kPoly1 || theBGMode==kPoly2 || theBGMode==kPoly3)
    mysumshape->plotOn(mframe2,Components(*mypolyshape),LineStyle(kDashed)) ;
  else if(theBGMode==kCheby1 || theBGMode==kCheby2 || theBGMode==kCheby3)
    mysumshape->plotOn(mframe2,Components(*mychebyshape),LineStyle(kDashed)) ;
  //RooExponential* myexpshape = (RooExponential*) w->pdf("expshape");
  //mysumshape->plotOn(mframe2,Components(*myexpshape),LineStyle(kDashed)) ;

  RooGaussian* sig1shape = (RooGaussian*) w->pdf("sig1");
  RooGaussian* sig2shape = (RooGaussian*) w->pdf("sig2");
  if(theSignalMode_== kDoubleGauss){
    mysumshape->plotOn(mframe2,Components(*sig1shape),LineStyle(kDashed)) ;
    mysumshape->plotOn(mframe2,Components(*sig2shape),LineStyle(kDashed)) ;
  }
  
  //TCanvas* c = new TCanvas("test","test",600,600);
  TCanvas* c = new TCanvas("thecanvas","the canvas",600,550);
  c->SetRightMargin(0.04);
  c->SetTopMargin(0.07);

  mframe2->Draw();

  drawPlotHeader(0.12);

  if(mumode) c->SaveAs("zmumu_mass_purityfit_inclusive.pdf");
  else c->SaveAs("zee_mass_purityfit_inclusive.pdf");
 
  //return; //for debugging

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //now fit the mass distribution in the VL sample, fixing all the
  //distribution parameters and floating only the background fraction
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  //TH1* h_zmass_FC = (TH1 *)f_zmass->Get("zeemass_FC_stack_doubleeledata");

  TH1* h_zmass_VL;
  if(mumode) h_zmass_VL = (TH1 *)f_zmass->Get("zmumumass_VL_stack_doublemudata");
  //if(mumode) h_zmass_VL = (TH1 *)f_zmass->Get("zmumumass_VL_stack_HT300_doublemudata");
  //if(mumode) h_zmass_VL = (TH1 *)f_zmass->Get("zmumumass_VL_stack_MET50to150_doublemudata");
  //if(mumode) h_zmass_VL = (TH1 *)f_zmass->Get("zmumumass_VL_stack_MET150to250_doublemudata");
  else h_zmass_VL = (TH1 *)f_zmass->Get("zeemass_VL_stack_doubleeledata");
  //else h_zmass_VL = (TH1 *)f_zmass->Get("zeemass_VL_stack_HT300_doubleeledata");
  //else h_zmass_VL = (TH1 *)f_zmass->Get("zeemass_VL_stack_MET50to150_doubleeledata");
  //else h_zmass_VL = (TH1 *)f_zmass->Get("zeemass_VL_stack_MET150to250_doubleeledata");
  
  RooDataHist dh_VL("dh_VL","dh_VL",*myvar,Import(*h_zmass_VL));

  RooDataSet data_VL("data","dataset with m",ztree_vl,*myvar) ;


  
  //RooPlot* mframe3 = myvar->frame(Title("Imported TH1 with Poisson error bars")) ;
  RooPlot* mframe3 = myvar->frame(Title("Imported dataset with Poisson error bars")) ;
  if(mumode) mframe3->SetXTitle("M(#mu#mu) [GeV]");
  else mframe3->SetXTitle("M(ee) [GeV]");

  mframe3->GetYaxis()->SetTitle("Events / 3 GeV"); //hard-coded for now

  //dh_VL.plotOn(mframe3) ; 
  data_VL.plotOn(mframe3,Binning(20)) ;   
  

  // **** fix signal parameters for the fit to the VL region **** //
  w->var("bw_mean")->setConstant();
  w->var("bw_width")->setConstant();
  w->var("v_mean")->setConstant();
  w->var("v_width")->setConstant();
  w->var("v_sigma")->setConstant();
  w->var("g1_mean")->setConstant();
  w->var("g1_sigma")->setConstant();
  w->var("g2_sigma")->setConstant();
  w->var("g_fract")->setConstant();
  w->var("cb_m0")->setConstant();
  w->var("cb_sigma")->setConstant();
  w->var("alpha")->setConstant();
  w->var("n")->setConstant();

  //Do we want to constrain background shape too?  Currently, no.
  //w->var("poly_c1")->setConstant();
  //w->var("poly_c2")->setConstant();
  //w->var("poly_c3")->setConstant();
  //w->var("cheby_c1")->setConstant();
  //w->var("cheby_c2")->setConstant();
  //w->var("cheby_c3")->setConstant();  
  

  //mysumshape->fitTo(dh_VL); //binned fit
  mysumshape->fitTo(data_VL); //unbineed fit
  mysumshape->plotOn(mframe3);
  RooPolynomial* mypolyshape_VL = (RooPolynomial*) w->pdf("polyshape");
  RooChebychev* mychebyshape_VL = (RooChebychev*) w->pdf("chebyshape");

  if(theBGMode==kPoly1 || theBGMode==kPoly2 || theBGMode==kPoly3)
    mysumshape->plotOn(mframe3,Components(*mypolyshape_VL),LineStyle(kDashed)) ;
  else if(theBGMode==kCheby1 || theBGMode==kCheby2 || theBGMode==kCheby3)
    mysumshape->plotOn(mframe3,Components(*mychebyshape_VL),LineStyle(kDashed)) ;

  if(theSignalMode_== kDoubleGauss){
    mysumshape->plotOn(mframe3,Components(*sig1shape),LineStyle(kDashed)) ;
    mysumshape->plotOn(mframe3,Components(*sig2shape),LineStyle(kDashed)) ;
  }  
  

  //TCanvas* c2 = new TCanvas("test2","test2",600,600);  
  TCanvas* c2 = new TCanvas("thecanvas2","the canvas2",600,550);
  c2->SetRightMargin(0.04);
  c2->SetTopMargin(0.07);

  mframe3->Draw();

  drawPlotHeader();

  if(mumode) c2->SaveAs("zmumu_mass_purityfit_VL.pdf");
  else c2->SaveAs("zee_mass_purityfit_VL.pdf");

  //w->var("fract")->Print();


  //////////////////////////////////////////////////////////////////
  //Compute the purity in the region of interest (i.e. |m-mZ|<15 )
  //using the fitted signal and background shapes
  //////////////////////////////////////////////////////////////////

  double mZ = 91.188;
  myvar->setRange("cut", mZ-15, mZ+15); 
  RooAbsReal* sigandbg = mysumshape->createIntegral (*myvar, NormSet (*myvar), Range("cut"));
  double IntegralValue_window = sigandbg->getVal();
  cout << "sigandbg = " << IntegralValue_window << endl; 

  myvar->setRange("cut", mZ-15, mZ+15); 
  RooAbsReal* sig=0;
  if(theSignalMode_==kBWxCB)
    sig = w->pdf("bw_x_cb")->createIntegral (*myvar, NormSet (*myvar), Range("cut"));
  else if(theSignalMode_==kVxCB)
    sig = w->pdf("v_x_cb")->createIntegral (*myvar, NormSet (*myvar), Range("cut"));
  else if(theSignalMode_==kDoubleGauss)
    sig = w->pdf("doublegauss")->createIntegral (*myvar, NormSet (*myvar), Range("cut"));

  double IntegralValue_sig_win = sig->getVal();
  cout << "sig = " << IntegralValue_sig_win << endl; 

  double fit_fract_val = w->var("fract")->getVal();
  double fit_fract_err = w->var("fract")->getError();
  cout << "fract = " << fit_fract_val
       << " +/- " << fit_fract_err << endl;
  //<< " + " << w->var("fract")->getAsymErrorHi()
  // << " - " << w->var("fract")->getAsymErrorLo() << endl;


  double purity = fit_fract_val * ( IntegralValue_sig_win / IntegralValue_window );
  double purity_err = fit_fract_err * ( IntegralValue_sig_win / IntegralValue_window );
  
  cout << "purity = " << purity << " +/- " << purity_err << endl;


  return;
}
