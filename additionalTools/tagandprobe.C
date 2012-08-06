//////////////////////////////////////////////////////////////////////////
// 
// Simple fitting routine to extract the lepton selection efficiency
// used for the 2011 Z->vv estimate
// 
// Global variables used to set configuration parameters.
//
// To make the input trees needed to do the fits, one needs
// to run tagAndProbeSkim.C.
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
//#include "RooGenericPdf.h"

#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TTree.h"
#include "TLatex.h"

using namespace RooFit ;

//////////////////////////////////////////////
// Options 
//////////////////////////////////////////////

const char *signalMode_str[]={ "sigBWxCB","sigCB","sigVxCB","sigDoubleGauss" ,"sigMCxGaus"};
enum signalMode {kBWxCB=0, kCB, kVxCB, kDoubleGauss,kMCxGaus};
//signalMode theSignalMode_ = kBWxCB;
//signalMode theSignalMode_ = kDoubleGauss;
signalMode theSignalMode_ = kMCxGaus; //this is now the best option (simply take the signal shape from MC)

const char *bgMode_str[]={ "bgPoly1","bgPoly2","bgPoly3","bgCheby1","bgCheby2","bgCheby3","bgExp" };
enum bgMode {kPoly1=0, kPoly2, kPoly3, kCheby1, kCheby2, kCheby3, kExp};
bgMode theBGMode_ = kExp;

//**** switch for Z->mumu vs Z->ee ****//
bool mumode = true;

//**** switch between fitting MC and data ****//
//(if this is on, then use the fitted shape from MC and apply to data (with gaussian convolution))
bool mcmode = false; 

//**** switch for fixing signal shape between pass and fail distribution ****// 
//(for data, the fit works best if you constain the signal shape to be
//the same between pass and fail distributions, unless you use kMCxGauss mode for signal) 
//for Zll MC, it works best if you use separate signal shapes
bool samePassFailSigShape = false;

const int nbins = 50;


//////////////////////////////////////////////

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



//////////////////////////////////////////////
// Main function
//////////////////////////////////////////////


void tagandprobe()
{

  //gStyle->SetOptStat(111111);


  //import the MC shape always (though it's only used in kMCxGaus mode)
  TString s_MCfilename = "tagandprobe_slim_ZJetsMC.root";
  TFile * f_zmass_MC = new TFile(s_MCfilename);

  TTree *ztree_MC_pass, *ztree_MC_fail;
  if(mumode){
    ztree_MC_pass = (TTree*) gDirectory->Get("tnpTreeMuPass");  
    ztree_MC_fail = (TTree*) gDirectory->Get("tnpTreeMuFail");  
  }
  else{
    ztree_MC_pass = (TTree*) gDirectory->Get("tnpTreeElePass");  
    ztree_MC_fail = (TTree*) gDirectory->Get("tnpTreeEleFail");  
  }    

  TH1F *histTemplatePass = new TH1F("histTemplatePass","histTemplatePass",nbins,60,120);
  ztree_MC_pass->Draw("m>>histTemplatePass","","goff");
  TH1F *histTemplateFail = new TH1F("histTemplateFail","histTemplateFail",nbins,60,120);
  ztree_MC_fail->Draw("m>>histTemplateFail","","goff");

  //import data points
  TString s_filename;
  if(mcmode) s_filename = "tagandprobe_slim_ZJetsMC.root";
  //if(mcmode) s_filename = "tagandprobe_slim_withHT_with3jetcut_ZJetsMC.root";
  //if(mcmode) s_filename = "tagandprobe_slim_withHT_ZJetsMC.root";
  else s_filename = "zmass_tnptree_noTrignoHT_samevar.root"; //the nominal data file!!
  //else s_filename = "tagandprobe_slim_noHT_geq2jetcut.root"; //the nominal data file - njets>=2
  //else s_filename = "tagandprobe_slim_noHT_singlemu.root"; // testing with singlemu PD!!
  //else s_filename = "tagandprobe_slim_noHT_geq1jetcut_singlemu.root"; // testing with singlemu PD!! >=1 JET CUT
  //else s_filename = "tagandprobe_slim_noHT_geq2jetcut_singlemu.root"; // testing with singlemu PD!! >=2 JET CUT
  //else s_filename = "tagandprobe_slim_noHT_singleelectron.root"; // testing with singleele PD!!
  //else s_filename = "tagandprobe_slim_noHT_nojetcut_singleelectron.root"; // testing with singleele PD - NO JET CUT
  //else s_filename = "tagandprobe_slim_noHT_geq1jetcut_singleelectron.root"; // testing with singleele PD - >=1 JET CUT
  //else s_filename = "tagandprobe_slim_noHT_geq2jetcut_singleelectron.root"; // testing with singleele PD - >=2 JET CUT
  //else s_filename = "tagandprobe_slim_noHT_geq3jetcut_IDbase_singleelectron.root"; // testing with singleele PD - >=3 JET CUT probe has ID cuts
  TFile * f_zmass = new TFile(s_filename);

  TTree* ztree_pass;
  TTree*  ztree_fail;
  if(mumode){
    ztree_pass = (TTree*) gDirectory->Get("tnpTreeMuPass");  
    ztree_fail = (TTree*) gDirectory->Get("tnpTreeMuFail");  
  }
  else{
    ztree_pass = (TTree*) gDirectory->Get("tnpTreeElePass");  
    ztree_fail = (TTree*) gDirectory->Get("tnpTreeEleFail");  
  }

  /*  Old way - now we use trees
  //Import data histogram
  //TFile * f_zmass = new TFile("zmass_tnp.root");
  //TFile * f_zmass = new TFile("zmass_tnp_noTrig.root");
  TFile * f_zmass = new TFile("zmass_tnp_noTrignoHT.root");
  //TFile * f_zmass = new TFile("zmass_tnp_zjetsMC.root");
  TH1D* h_zmass_muprobe1_pass = (TH1D *)f_zmass->Get("muprobe1_mass_pass");
  TH1D* h_zmass_muprobe2_pass = (TH1D *)f_zmass->Get("muprobe2_mass_pass");
  TH1D* h_zmass_muprobe3_pass = (TH1D *)f_zmass->Get("muprobe3_mass_pass");
  TH1D* h_zmass_muprobe1_fail = (TH1D *)f_zmass->Get("muprobe1_mass_fail");
  TH1D* h_zmass_muprobe2_fail = (TH1D *)f_zmass->Get("muprobe2_mass_fail");
  TH1D* h_zmass_muprobe3_fail = (TH1D *)f_zmass->Get("muprobe3_mass_fail");

  TH1D* h_zmass_eleprobe1_pass = (TH1D *)f_zmass->Get("eleprobe1_mass_pass");
  TH1D* h_zmass_eleprobe2_pass = (TH1D *)f_zmass->Get("eleprobe2_mass_pass");
  TH1D* h_zmass_eleprobe3_pass = (TH1D *)f_zmass->Get("eleprobe3_mass_pass");
  TH1D* h_zmass_eleprobe1_fail = (TH1D *)f_zmass->Get("eleprobe1_mass_fail");
  TH1D* h_zmass_eleprobe2_fail = (TH1D *)f_zmass->Get("eleprobe2_mass_fail");
  TH1D* h_zmass_eleprobe3_fail = (TH1D *)f_zmass->Get("eleprobe3_mass_fail");

  //combine all the parts
  TH1D* h_zmass_sum_muprobe_pass = (TH1D*) h_zmass_muprobe1_pass->Clone();
  h_zmass_sum_muprobe_pass->Sumw2();
  h_zmass_sum_muprobe_pass->Add(h_zmass_muprobe2_pass);
  h_zmass_sum_muprobe_pass->Add(h_zmass_muprobe3_pass);
  TH1D* h_zmass_sum_muprobe_fail = (TH1D*) h_zmass_muprobe1_fail->Clone();
  h_zmass_sum_muprobe_fail->Sumw2();
  h_zmass_sum_muprobe_fail->Add(h_zmass_muprobe2_fail);
  h_zmass_sum_muprobe_fail->Add(h_zmass_muprobe3_fail);
  TH1D* h_zmass_sum_eleprobe_pass = (TH1D*) h_zmass_eleprobe1_pass->Clone();
  h_zmass_sum_eleprobe_pass->Sumw2();
  h_zmass_sum_eleprobe_pass->Add(h_zmass_eleprobe2_pass);
  h_zmass_sum_eleprobe_pass->Add(h_zmass_eleprobe3_pass);
  TH1D* h_zmass_sum_eleprobe_fail = (TH1D*) h_zmass_eleprobe1_fail->Clone();
  h_zmass_sum_eleprobe_fail->Sumw2();
  h_zmass_sum_eleprobe_fail->Add(h_zmass_eleprobe2_fail);
  h_zmass_sum_eleprobe_fail->Add(h_zmass_eleprobe3_fail);
  */


  //h_zmass_sum_muprobe_pass->Draw();



  //////////////////////////////////////////
  //
  //put all the pdfs in a workspace
  //
  //////////////////////////////////////////

  RooWorkspace* w = new RooWorkspace("w");

  TString massvar = "m";
  //TString massvar = "muprobe_mass_pass";

  //from PhysicsTools/TagAndProbe/src/TagProbeFitter.cc
  //pdf framework to extract measured efficiency *as a parameter* (i.e. with error evaluation)
  w->factory("EXPR::nSignalPass('efficiency*fSigAll*numTot', efficiency[0,1], fSigAll[.9,0,1],numTot[1,0,1e10])");
  //w->factory("EXPR::nSignalPass('efficiency*fSigAll*numTot', efficiency[0.7,0,1], fSigAll[.9,0,1],numTot[1,0,1e10])");//testing single electron!
  w->factory("EXPR::nSignalFail('(1-efficiency)*fSigAll*numTot', efficiency, fSigAll,numTot)");  
  w->factory("EXPR::nBkgPass('effBkg*(1-fSigAll)*numTot', effBkg[.9,0,1],fSigAll,numTot)");
  w->factory("EXPR::nBkgFail('(1-effBkg)*(1-fSigAll)*numTot', effBkg,fSigAll,numTot)");  


  //import the MC histgoram
  w->factory(massvar+"[60,120]");

  RooRealVar* thevar = w->var(massvar);    
  //thevar->setBins(nbins,"fft");

  RooDataHist dh_Pass("dh_Pass","boogalie",*thevar,Import(*histTemplatePass));
  RooDataHist dh_Fail("dh_Fail","boogalie",*thevar,Import(*histTemplateFail));
  
  w->import(dh_Pass);
  w->import(dh_Fail);



  /////////////////////////////
  //
  //Specify the signal pdf  
  //
  /////////////////////////////

  TString theSignalModel, theSignalModelFail;
  if(theSignalMode_ == kBWxCB || theSignalMode_ == kVxCB){
    
    //w->factory("CBShape::cbshape(m[60,120],cb_m0[0.0,-5,1.0],cb_sigma[1.8,0.0,15],alpha[0.959,0.3,2],n[1.266,0.5,10])");
    //TString s_CB = "CBShape::cbshape("+massvar+"[60,120],cb_m0[0.0,-5,1.0],cb_sigma[1.8,0,15],alpha[0.959,0.3,4],n[1.266,-10.0,10])";
    TString s_CB = "CBShape::cbshape("+massvar+"[60,120],cb_m0[0.0,-5,1.0],cb_sigma[3,0,15],alpha[0.959,0.3,4],n[1.266,-10.0,10])";
    //TString s_CB = "CBShape::cbshape("+massvar+"[60,120],cb_m0[0.0,-5,2.0],cb_sigma[1.8,0.0,15],alpha[0.959,0.3,4],n[1.266,-10.0,10])";
    TString s_CB_fail = "CBShape::cbshapeFail("+massvar+"[60,120],cb_m0Fail[0.0,-5,1.0],cb_sigmaFail[1.8,0.0,15],alpha[0.959,0.3,4],nFail[1.266,-10.0,10])";
    w->factory(s_CB);
    w->factory(s_CB_fail);

    if(theSignalMode_ == kBWxCB){
      //w->factory("BreitWigner::bwshape("+massvar+"[60,120],bw_mean[91.2,89,93],bw_width[2.49,2,3])");
      //w->factory("BreitWigner::bwshape("+massvar+"[60,120],bw_mean[91.2,89,93],bw_width[2.7,2,5])");

      //this is the one we've been using
      TString s_BW =  "BreitWigner::bwshape("+massvar+"[60,120],bw_mean[91.2,89,93],bw_width[2.7,2,9])";

      //TString s_BW =  "BreitWigner::bwshape("+massvar+"[60,120],bw_mean[91.2,89,93],bw_width[2.7,2,15])";
      TString s_BW_fail = "BreitWigner::bwshapeFail("+massvar+"[60,120],bw_meanFail[91.2,89,93],bw_widthFail[2.7,2,9])";

      //w->factory("BreitWigner::bwshape(m[60,120],bw_mean[91.2,89,93],bw_width[2.7,2,5])");
      //w->factory("CBShape::cbshape(m[60,120],cb_m0[0.4,-10.0,3.0],cb_sigma[2.3,0.0,15],alpha[0.959,0.3,2],n[8.5,0.5,10])");
      //w->factory("CBShape::cbshape(m[60,120],cb_m0[0.0,-5,1.0],cb_sigma[2.3,0.0,15],alpha[0.959,0.3,2],n[1,0.5,10])");

      w->factory(s_BW);
      w->factory(s_BW_fail);

      TString s_BWxCB = "FCONV::bw_x_cb("+massvar+",bwshape,cbshape)";
      TString s_BWxCB_fail = "FCONV::bw_x_cbFail("+massvar+",bwshapeFail,cbshapeFail)";

      w->factory(s_BWxCB);
      w->factory(s_BWxCB_fail);

      theSignalModel = "bw_x_cb";
      theSignalModelFail = theSignalModel+"Fail";
    }
    else if(theSignalMode_ ==kVxCB){
      TString s_V = "Voigtian::voigtshape("+massvar+"[60,120],v_mean[91.2,89,93],v_width[2.5,2,5],v_sigma[1,0.1,10])";
      //TString s_V = "Voigtian::voigtshapeFail("+massvar+"[60,120],v_meanFail[91.2,89,93],v_widthFail[2.5,2,5],v_sigmaFail[1,0.1,10])";

      TString s_VxCB = "FCONV::v_x_cb("+massvar+",voigtshape,cbshape)";
      w->factory(s_V);
      w->factory(s_VxCB);

      theSignalModel = "v_x_cb";
    }
  }
  else if (theSignalMode_ == kCB){
    w->factory("CBShape::cbshape2(m[60,120],cb_m02[91,89,93],cb_sigma2[2.3,0.0,15],alpha2[0.959,0.3,2],n2[1.3,0.5,10])");

    theSignalModel = "cbshape2";
  }
  else if(theSignalMode_ == kDoubleGauss){
    TString s_gaus1 = "Gaussian::sig1("+massvar+"[60,120],g1_mean[91,87,93],g1_sigma[2.5,0,5])";
    TString s_gaus2 = "Gaussian::sig2("+massvar+"[60,120],g1_mean,g2_sigma[9,0,20])";
    TString s_doublegaus = "SUM::doublegauss(g_fract[0.9,0.5,1]*sig1,sig2)";
    w->factory(s_gaus1);
    w->factory(s_gaus2);
    w->factory(s_doublegaus);   

    theSignalModel = "doublegauss";   
  }
  else if(theSignalMode_ == kMCxGaus){

    TString s_histPdfPass = "HistPdf::hpPass("+massvar+"[60,120],dh_Pass)";
    TString s_histPdfFail = "HistPdf::hpFail("+massvar+"[60,120],dh_Fail)";
    
    w->factory(s_histPdfPass);
    w->factory(s_histPdfFail);

    TString s_gaussianPass = "GaussModel::gausPass("+massvar+"[60,120],gausPassMean[0,-1,1],gausPassSigma[2.5,0,10])";
    TString s_gaussianFail = "GaussModel::gausFail("+massvar+"[60,120],gausFailMean[0,-1,1],gausFailSigma[2.5,0,5])";
    //TString s_gaussianFail = "GaussModel::gausFail("+massvar+"[60,120],gausFailMean[0,-1,1],gausFailSigma[4,0,50])";//testing it for singleelectron pd

    w->factory(s_gaussianPass);
    w->factory(s_gaussianFail);


    TString s_MCxGaus = "FCONV::mc_x_gaus("+massvar+",hpPass,gausPass)";
    TString s_MCxGausFail = "FCONV::mc_x_gausFail("+massvar+",hpFail,gausFail)";
    //TString s_MCxGaus = "FCONV::mc_x_gaus("+massvar+",gausPass,hpPass)";
    //TString s_MCxGausFail = "FCONV::mc_x_gausFail("+massvar+",gausFail,hpFail)";
    w->factory(s_MCxGaus);
    w->factory(s_MCxGausFail);

    theSignalModel = "mc_x_gaus";
    theSignalModelFail = theSignalModel+"Fail";
    
  }


  ///////////////////////////////////
  //
  //Specify the background pdf  
  //
  ///////////////////////////////////

  TString theBGModelPass, theBGModelFail;
  if(theBGMode_==kPoly1 || theBGMode_==kPoly2 || theBGMode_==kPoly3){  
    //w->factory("Polynomial::polyshape("+massvar+"[60,120],{poly_c1[0.1,-1,1],poly_c2[0.1,-1,1]})");
    //w->factory("Polynomial::polyshape("+massvar+"[60,120],{poly_c1[0,-1,1],poly_c2[0.1,-1,1]})"); // for fitting linear background
    //w->factory("Polynomial::polyshape("+massvar+"[60,120],{poly_c1[0.1,-1,1],poly_c2[0.1,-1,1],poly_c3[0.1,-1,1]})"); // for fitting cubic background
    w->factory("Polynomial::polyshape("+massvar+"[60,120],{poly_c1[0.0,-1,1],poly_c2[0.1,-1,1],poly_c3[0.1,-1,1]})"); // for fitting cubic background
    w->factory("Polynomial::polyshapeFail("+massvar+"[60,120],{poly_c1Fail[0.0,-1,1],poly_c2Fail[0.1,-1,1],poly_c3Fail[0.1,-1,1]})"); 
    //w->factory("Polynomial::polyshapeFail("+massvar+"[60,120],{poly_c1Fail[-0.5,-1,0],poly_c2Fail[-0.1,-1,0],poly_c3Fail[0.1,-1,1]})"); //testing for singleelectron njet=0

    if(theBGMode_!=kPoly3){ 
      w->var("poly_c3")->setVal(0);  w->var("poly_c3")->setConstant();
      w->var("poly_c3Fail")->setVal(0);  w->var("poly_c3Fail")->setConstant();
    }
    if(theBGMode_!=kPoly2){ 
      w->var("poly_c2")->setVal(0);  w->var("poly_c2")->setConstant();
      w->var("poly_c2Fail")->setVal(0);  w->var("poly_c2Fail")->setConstant();
    }


    ////TEMPORARY FOR MC FITS, SET BACKGROUND TO 0
    //w->var("poly_c1")->setVal(0);
    //w->var("poly_c1")->setConstant();
    //w->var("poly_c1Fail")->setVal(0);
    //w->var("poly_c1Fail")->setConstant();

    theBGModelPass = "polyshape"; 
    theBGModelFail = theBGModelPass+"Fail";

  }
  else if(theBGMode_==kCheby1 || theBGMode_==kCheby2 || theBGMode_==kCheby3){
    //w->factory("Chebychev::chebyshape("+massvar+"[60,120],{cheby_c1[0.1,-1.5,1],cheby_c2[0.1,-1,1],cheby_c3[0.1,-1,1]})");
    w->factory("Chebychev::chebyshape("+massvar+"[60,120],{cheby_c1[0.0,-1.5,1],cheby_c2[0.0,-1,1],cheby_c3[0.0,-1,1]})");


    //w->factory("Chebychev::chebyshapeFail("+massvar+"[60,120],{cheby_c1Fail[0.0,-1.5,1],cheby_c2Fail[0.0,-1,1],cheby_c3Fail[0.0,-1,1]})");
    w->factory("Chebychev::chebyshapeFail("+massvar+"[60,120],{cheby_c1Fail[-0.5,-1.5,0],cheby_c2Fail[-0.5,-1,0],cheby_c3Fail[0.0,-1,1]})");//testing for singleelectron njet>=3

    if(theBGMode_!=kCheby3){ 
      w->var("cheby_c3")->setVal(0);  w->var("cheby_c3")->setConstant();
      w->var("cheby_c3Fail")->setVal(0);  w->var("cheby_c3Fail")->setConstant();
    }
    if(theBGMode_!=kCheby2){ 
      w->var("cheby_c2")->setVal(0);  w->var("cheby_c2")->setConstant();
      w->var("cheby_c2Fail")->setVal(0);  w->var("cheby_c2Fail")->setConstant();
    }


    theBGModelPass = "chebyshape"; 
    theBGModelFail = theBGModelPass+"Fail";
  }
  else if(theBGMode_ == kExp){
    w->factory("Exponential::expshape("+massvar+"[60,120],lambda[0,-10,10])");
    //w->factory("Exponential::expshapeFail("+massvar+"[60,120],lambdaFail[0,-10,10])");
    w->factory("Exponential::expshapeFail("+massvar+"[60,120],lambdaFail[-5,-10,0])");//test singleelectron

    theBGModelPass = "expshape"; 
    theBGModelFail = theBGModelPass+"Fail";
  }

  /////////////////////////////////////////////
  //
  // Combine Signal and Background PDFs
  //
  /////////////////////////////////////////////

  w->factory("SUM::pdfPass(nSignalPass * "+theSignalModel+", nBkgPass*"+theBGModelPass+")");

  if(samePassFailSigShape){
    w->factory("SUM::pdfFail(nSignalFail * "+theSignalModel+", nBkgFail*"+theBGModelFail+")");
  }
  else{
    w->factory("SUM::pdfFail(nSignalFail * "+theSignalModelFail+", nBkgFail*"+theBGModelFail+")");
  }

  w->factory("SIMUL::simPdf(simindex[Passed=1,Failed=0], Passed=pdfPass, Failed=pdfFail)");


  // Print workspace contents
  w->Print();  

  //return;

  //TH1D* myhistpointer;
  //if(mumode) myhistpointer = h_zmass_sum_muprobe_pass; 
  //else myhistpointer = h_zmass_sum_eleprobe_pass; 
  
  RooRealVar* myvar = w->var(massvar);
  //RooDataHist dh("dh","dh",*myvar,Import(*myhistpointer));

  //RooDataSet data_pass("data","dataset with m",ztree_pass,*myvar) ;
  //RooDataSet data_fail("data","dataset with m",ztree_fail,*myvar) ;

  //needed for ZLL MC? (too many stats for unbinned fit?)
  TH1F *myhist = new TH1F("myhist","myhist",nbins,60,120);
  ztree_pass->Draw("m>>myhist"); 
  RooDataHist data_pass("data","dataset with m",*myvar,Import(*myhist));
  ztree_fail->Draw("m>>myhist"); 
  RooDataHist data_fail("data","dataset with m",*myvar,Import(*myhist));

  
  // Make joint dataset
  //w->factory("index[Passed,Failed]") ;
  //RooDataSet* joint = new RooDataSet("joint","joint",RooArgSet(w::x,w::wgt),Index(w::index),Import("Passed",data_pass),Import("Failed",data_fail),WeightVar(w::wgt)) ; 
  //RooDataSet jointData("jointData", "Joint Dataset", RooArgSet(*ws.var("nOn"), *ws.var("nOff")), Index(setLabel), Import(fakeData));

  RooCategory* mycat = w->cat("simindex");


  //RooDataSet jointData("jointData", "Joint Dataset", *myvar, Index( *mycat ),Import("Passed",data_pass),Import("Failed",data_fail));

  RooDataHist jointData("jointData", "Joint Dataset", *myvar, Index( *mycat ),Import("Passed",data_pass),Import("Failed",data_fail));


  double totPassing = data_pass.sumEntries();
  double totFailing = data_fail.sumEntries();

  cout << "totPassing = " << totPassing << endl;
  cout << "totFailing = " << totFailing << endl;

  //set initial value for N_tot
  w->var("numTot")->setVal(totPassing+totFailing);


  /////////////////////////////////////////////
  //
  // Fit!
  //
  /////////////////////////////////////////////
  
  //w->pdf("simPdf")->fitTo(jointData,Extended(true));
  RooAbsPdf* mymodel = w->pdf("simPdf");
  mymodel->fitTo(jointData,Extended(true));
  //mymodel->fitTo(jointData, Save(true), Extended(true), NumCPU(2), Strategy(2), PrintLevel(1), PrintEvalErrors(1), Warnings(1));
  
  /////////////////////////////////////////////
  //
  // Plot results
  //
  /////////////////////////////////////////////

  RooAbsData *dataPass = jointData.reduce(Cut("simindex == 1"));
  RooAbsData *dataFail = jointData.reduce(Cut("simindex == 0"));

  TString s_xtitle;

  if(mumode)  s_xtitle = "M(#mu#mu) [GeV]";
  else s_xtitle = "M(ee) [GeV]";

  //TCanvas* c1 = new TCanvas("test","test",1200,550);
  TCanvas* c1 = new TCanvas("test","test",600,550);
  c1->SetRightMargin(0.04);
  c1->SetTopMargin(0.07);  
  //c1->Divide(2,1);
  //c1->cd(1);
  RooPlot *plotPass = myvar->frame(Bins(20), Title("Pass"));
  data_pass.plotOn(plotPass);
  mymodel->plotOn(plotPass, Slice(*mycat, "Passed"), ProjWData(*dataPass), LineColor(kBlue), Components(theBGModelPass), LineStyle(kDashed));
  mymodel->plotOn(plotPass, Slice(*mycat, "Passed"), ProjWData(*dataPass), LineColor(kBlue));
  plotPass->SetXTitle(s_xtitle);
  plotPass->GetYaxis()->SetTitle("Events / 1.2 GeV"); //hard-coded for now
  plotPass->Draw();
  drawPlotHeader();

  
  TCanvas* c2 = new TCanvas("test2","test2",600,550);
  //c1->cd(2);
  c2->SetRightMargin(0.04);
  c2->SetTopMargin(0.07);  
  RooPlot *plotFail = myvar->frame(Bins(20), Title("Fail"));
  data_fail.plotOn(plotFail);
  mymodel->plotOn(plotFail, Slice(*mycat, "Failed"), ProjWData(*dataFail), LineColor(kBlue), Components(theBGModelFail), LineStyle(kDashed));
  //mymodel->plotOn(plotFail, Slice(*mycat, "Failed"), ProjWData(*dataFail), LineColor(kGreen), Components(theSignalModelFail));//signal fit only
  mymodel->plotOn(plotFail, Slice(*mycat, "Failed"), ProjWData(*dataFail), LineColor(kBlue));
  plotFail->SetXTitle(s_xtitle);
  plotFail->GetYaxis()->SetTitle("Events / 1.2 GeV"); //hard-coded for now
  plotFail->Draw();
  drawPlotHeader();
  
  TString s_sigmode = signalMode_str[theSignalMode_];
  TString s_bgmode = bgMode_str[theBGMode_];


  if(mumode){
    c1->SaveAs("zmumu_mass_tnpfit_pass_"+s_sigmode+"_"+s_bgmode+".pdf");
    c2->SaveAs("zmumu_mass_tnpfit_fail_"+s_sigmode+"_"+s_bgmode+".pdf");
  }
  else{
    c1->SaveAs("zee_mass_tnpfit_pass_"+s_sigmode+"_"+s_bgmode+".pdf");
    c2->SaveAs("zee_mass_tnpfit_fail_"+s_sigmode+"_"+s_bgmode+".pdf");
  }
  

  return;
 
}
