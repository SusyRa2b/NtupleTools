#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "../MiscUtil.cxx"
#include <cassert>
#include <iostream>

using namespace std;




void printABCD(double A, double Aerr, double B, double Berr, double C, double Cerr, double D, double Derr) {
  
  double numerr=jmt::errAtimesB(B,Berr,D,Derr);
  double num = B*D;
  double estimate = num / A;
  double estimateerr= jmt::errAoverB(num,numerr,A,Aerr);
  
  cout << "prediction: " << estimate << " +- " << estimateerr << endl;
  cout << "MC truth  : " << C << " +- " << Cerr << endl;
  
  double num1 = C-estimate;
  double num1err = sqrt(Cerr*Cerr + estimateerr*estimateerr);
  double perc = num1/C;
  double percerr = jmt::errAoverB(num1,num1err,C,Cerr);
  cout << "percent   : " << perc*100. << " +- " << percerr*100. << endl;
  cout << "combined  : " << sqrt(perc*perc + percerr*percerr)*100. << endl;
}



void reweightQCD(){

  //TFile* ifl = TFile::Open("../qcdReweight_loose_50.root","READ"); 
  //TFile* ifl = TFile::Open("../qcdReweight_loose_30.root","READ"); 
  //TFile* ifl = TFile::Open("../qcdReweight_tight_50.root","READ"); 
  TFile* ifl = TFile::Open("../qcdReweight_tight_30.root","READ"); 
  
  bool d = true; //true to double size of correction

  TH1D* h_LSB_fDP_antib_njets_qcd    = (TH1D*)ifl->Get("h_LSB_fDP_antib_njets_qcd");
  TH1D* h_LSB_fDP_antib_njets_nonqcd = (TH1D*)ifl->Get("h_LSB_fDP_antib_njets_nonqcd");
  TH1D* h_LSB_fDP_antib_njets_data   = (TH1D*)ifl->Get("h_LSB_fDP_antib_njets_data");
  TH1D* h_LSB_pDP_antib_njets_qcd    = (TH1D*)ifl->Get("h_LSB_pDP_antib_njets_qcd");

  TH1D* h_SB_fDP_ge1b_njets_qcd    = (TH1D*)ifl->Get("h_SB_fDP_ge1b_njets_qcd");
  TH1D* h_SB_fDP_ge1b_njets_nonqcd = (TH1D*)ifl->Get("h_SB_fDP_ge1b_njets_nonqcd");
  TH1D* h_SB_fDP_ge1b_njets_data   = (TH1D*)ifl->Get("h_SB_fDP_ge1b_njets_data");
  TH1D* h_SB_pDP_ge1b_njets_qcd    = (TH1D*)ifl->Get("h_SB_pDP_ge1b_njets_qcd");
  TH1D* h_SB_fDP_ge2b_njets_qcd    = (TH1D*)ifl->Get("h_SB_fDP_ge2b_njets_qcd");
  TH1D* h_SB_fDP_ge2b_njets_nonqcd = (TH1D*)ifl->Get("h_SB_fDP_ge2b_njets_nonqcd");
  TH1D* h_SB_fDP_ge2b_njets_data   = (TH1D*)ifl->Get("h_SB_fDP_ge2b_njets_data");
  TH1D* h_SB_pDP_ge2b_njets_qcd    = (TH1D*)ifl->Get("h_SB_pDP_ge2b_njets_qcd");

  TH1D* h_SIG_fDP_ge1b_njets_qcd    = (TH1D*)ifl->Get("h_SIG_fDP_ge1b_njets_qcd");
  TH1D* h_SIG_fDP_ge1b_njets_nonqcd = (TH1D*)ifl->Get("h_SIG_fDP_ge1b_njets_nonqcd");
  TH1D* h_SIG_fDP_ge1b_njets_data   = (TH1D*)ifl->Get("h_SIG_fDP_ge1b_njets_data");
  TH1D* h_SIG_pDP_ge1b_njets_qcd    = (TH1D*)ifl->Get("h_SIG_pDP_ge1b_njets_qcd");
  TH1D* h_SIG_fDP_ge2b_njets_qcd    = (TH1D*)ifl->Get("h_SIG_fDP_ge2b_njets_qcd");
  TH1D* h_SIG_fDP_ge2b_njets_nonqcd = (TH1D*)ifl->Get("h_SIG_fDP_ge2b_njets_nonqcd");
  TH1D* h_SIG_fDP_ge2b_njets_data   = (TH1D*)ifl->Get("h_SIG_fDP_ge2b_njets_data");
  TH1D* h_SIG_pDP_ge2b_njets_qcd    = (TH1D*)ifl->Get("h_SIG_pDP_ge2b_njets_qcd");


  //clone for weight histograms
  TH1D* wh_LSB_antib = (TH1D*)h_LSB_fDP_antib_njets_qcd->Clone("wh_LSB_antib");
  TH1D* wh_SB_ge1b   = (TH1D*)h_SB_fDP_ge1b_njets_qcd->Clone("wh_SB_ge1b");
  TH1D* wh_SIG_ge1b  = (TH1D*)h_SIG_fDP_ge1b_njets_qcd->Clone("wh_SIG_ge1b");
  TH1D* wh_SB_ge2b   = (TH1D*)h_SB_fDP_ge2b_njets_qcd->Clone("wh_SB_ge2b");
  TH1D* wh_SIG_ge2b  = (TH1D*)h_SIG_fDP_ge2b_njets_qcd->Clone("wh_SIG_ge2b");
  wh_LSB_antib->Reset();
  wh_SB_ge1b->Reset();
  wh_SIG_ge1b->Reset();
  wh_SB_ge2b->Reset();
  wh_SIG_ge2b->Reset();
  wh_LSB_antib->SetTitle("wh_LSB_antib");
  wh_SB_ge1b->SetTitle("wh_SB_ge1b");  
  wh_SIG_ge1b->SetTitle("wh_SIG_ge1b");
  wh_SB_ge2b->SetTitle("wh_SB_ge2b");
  wh_SIG_ge2b->SetTitle("wh_SIG_ge2b");
  TString xtitle = "njets";
  TString ytitle = "scale factor";
  wh_LSB_antib->GetXaxis()->SetTitle(xtitle);
  wh_SB_ge1b  ->GetXaxis()->SetTitle(xtitle);  
  wh_SIG_ge1b ->GetXaxis()->SetTitle(xtitle);
  wh_SB_ge2b  ->GetXaxis()->SetTitle(xtitle);
  wh_SIG_ge2b ->GetXaxis()->SetTitle(xtitle);
  wh_LSB_antib->GetYaxis()->SetTitle(ytitle);
  wh_SB_ge1b  ->GetYaxis()->SetTitle(ytitle);  
  wh_SIG_ge1b ->GetYaxis()->SetTitle(ytitle);
  wh_SB_ge2b  ->GetYaxis()->SetTitle(ytitle);
  wh_SIG_ge2b ->GetYaxis()->SetTitle(ytitle);
  
  //clone for weighted histograms
  TH1D* h_LSB_fDP_antib_njets_qcd_w = (TH1D*)h_LSB_fDP_antib_njets_qcd->Clone("h_LSB_fDP_antib_njets_qcd_w");
  TH1D* h_LSB_pDP_antib_njets_qcd_w = (TH1D*)h_LSB_pDP_antib_njets_qcd->Clone("h_LSB_pDP_antib_njets_qcd_w");
  TH1D* h_SB_fDP_ge1b_njets_qcd_w = (TH1D*)h_SB_fDP_ge1b_njets_qcd->Clone("h_SB_fDP_ge1b_njets_qcd_w");
  TH1D* h_SB_pDP_ge1b_njets_qcd_w = (TH1D*)h_SB_pDP_ge1b_njets_qcd->Clone("h_SB_pDP_ge1b_njets_qcd_w");
  TH1D* h_SIG_fDP_ge1b_njets_qcd_w = (TH1D*)h_SIG_fDP_ge1b_njets_qcd->Clone("h_SIG_fDP_ge1b_njets_qcd_w");
  TH1D* h_SIG_pDP_ge1b_njets_qcd_w = (TH1D*)h_SIG_pDP_ge1b_njets_qcd->Clone("h_SIG_pDP_ge1b_njets_qcd_w");
  TH1D* h_SB_fDP_ge2b_njets_qcd_w = (TH1D*)h_SB_fDP_ge2b_njets_qcd->Clone("h_SB_fDP_ge2b_njets_qcd_w");
  TH1D* h_SB_pDP_ge2b_njets_qcd_w = (TH1D*)h_SB_pDP_ge2b_njets_qcd->Clone("h_SB_pDP_ge2b_njets_qcd_w");
  TH1D* h_SIG_fDP_ge2b_njets_qcd_w = (TH1D*)h_SIG_fDP_ge2b_njets_qcd->Clone("h_SIG_fDP_ge2b_njets_qcd_w");
  TH1D* h_SIG_pDP_ge2b_njets_qcd_w = (TH1D*)h_SIG_pDP_ge2b_njets_qcd->Clone("h_SIG_pDP_ge2b_njets_qcd_w");
  h_LSB_fDP_antib_njets_qcd_w->Reset();
  h_LSB_pDP_antib_njets_qcd_w->Reset();
  h_SB_fDP_ge1b_njets_qcd_w->Reset();
  h_SB_pDP_ge1b_njets_qcd_w->Reset();
  h_SIG_fDP_ge1b_njets_qcd_w->Reset();
  h_SIG_pDP_ge1b_njets_qcd_w->Reset();
  h_SB_fDP_ge2b_njets_qcd_w->Reset();
  h_SB_pDP_ge2b_njets_qcd_w->Reset();
  h_SIG_fDP_ge2b_njets_qcd_w->Reset();
  h_SIG_pDP_ge2b_njets_qcd_w->Reset();
  
  
  //make weights

  
  if(d){
    wh_LSB_antib->Add(h_LSB_fDP_antib_njets_data, 2);
    wh_LSB_antib->Add(h_LSB_fDP_antib_njets_nonqcd, -2);
    wh_LSB_antib->Add( h_LSB_fDP_antib_njets_qcd, -1);
    wh_LSB_antib->Divide(h_LSB_fDP_antib_njets_qcd);
  }
  else{
    wh_LSB_antib->Add(h_LSB_fDP_antib_njets_data);
    wh_LSB_antib->Add(h_LSB_fDP_antib_njets_nonqcd,-1);
    wh_LSB_antib->Divide(h_LSB_fDP_antib_njets_qcd);
  }
  for(int i=1; i<= wh_LSB_antib->GetNbinsX(); i++){
    assert(h_LSB_fDP_antib_njets_data->GetBinContent(i)>0);
    assert(h_LSB_fDP_antib_njets_qcd->GetBinContent(i)>0);
    if(wh_LSB_antib->GetBinContent(i)<0.) wh_LSB_antib->SetBinContent(i,0);
  }
  
  if(d){
    wh_SB_ge1b->Add(h_SB_fDP_ge1b_njets_data, 2);
    wh_SB_ge1b->Add(h_SB_fDP_ge1b_njets_nonqcd, -2);
    wh_SB_ge1b->Add(h_SB_fDP_ge1b_njets_qcd, -1);
    wh_SB_ge1b->Divide(h_SB_fDP_ge1b_njets_qcd);
  }
  else{
    wh_SB_ge1b->Add(h_SB_fDP_ge1b_njets_data);
    wh_SB_ge1b->Add(h_SB_fDP_ge1b_njets_nonqcd, -1);
    wh_SB_ge1b->Divide(h_SB_fDP_ge1b_njets_qcd);
  }
  for(int i=1; i<= wh_SB_ge1b->GetNbinsX(); i++){
    assert(h_SB_fDP_ge1b_njets_data->GetBinContent(i)>0);
    assert(h_SB_fDP_ge1b_njets_qcd->GetBinContent(i)>0);
    if(wh_SB_ge1b->GetBinContent(i)<0) wh_SB_ge1b->SetBinContent(i,0);
  }
  
  
  if(d){
    wh_SIG_ge1b->Add(h_SIG_fDP_ge1b_njets_data, 2);
    wh_SIG_ge1b->Add(h_SIG_fDP_ge1b_njets_nonqcd, -2);
    wh_SIG_ge1b->Add(h_SIG_fDP_ge1b_njets_qcd, -1);
    wh_SIG_ge1b->Divide(h_SIG_fDP_ge1b_njets_qcd);
  }
  else{
    wh_SIG_ge1b->Add(h_SIG_fDP_ge1b_njets_data);
    wh_SIG_ge1b->Add(h_SIG_fDP_ge1b_njets_nonqcd, -1);
    wh_SIG_ge1b->Divide(h_SIG_fDP_ge1b_njets_qcd);
  }
  for(int i=1; i<= wh_SIG_ge1b->GetNbinsX(); i++){
    assert(h_SIG_fDP_ge1b_njets_data->GetBinContent(i)>0);
    assert(h_SIG_fDP_ge1b_njets_qcd->GetBinContent(i)>0);
    if(wh_SIG_ge1b->GetBinContent(i)<0) wh_SIG_ge1b->SetBinContent(i,0);
  }
  
  if(d){
    wh_SB_ge2b->Add(h_SB_fDP_ge2b_njets_data, 2);
    wh_SB_ge2b->Add(h_SB_fDP_ge2b_njets_nonqcd, -2);
    wh_SB_ge2b->Add(h_SB_fDP_ge2b_njets_qcd, -1);
    wh_SB_ge2b->Divide(h_SB_fDP_ge2b_njets_qcd);
  }
  else{
    wh_SB_ge2b->Add(h_SB_fDP_ge2b_njets_data);
    wh_SB_ge2b->Add(h_SB_fDP_ge2b_njets_nonqcd, -1);
    wh_SB_ge2b->Divide(h_SB_fDP_ge2b_njets_qcd);
  }
  for(int i=1; i<= wh_SB_ge1b->GetNbinsX(); i++){
    assert(h_SB_fDP_ge1b_njets_data->GetBinContent(i)>0);
    assert(h_SB_fDP_ge1b_njets_qcd->GetBinContent(i)>0);
    if(wh_SB_ge1b->GetBinContent(i)<0) wh_SB_ge1b->SetBinContent(i,0);
  }
  for(int i=1; i<= wh_SB_ge2b->GetNbinsX(); i++){
    assert(h_SB_fDP_ge2b_njets_data->GetBinContent(i)>0);
    assert(h_SB_fDP_ge2b_njets_qcd->GetBinContent(i)>0);
    if(wh_SB_ge2b->GetBinContent(i)<0) wh_SB_ge2b->SetBinContent(i,0);
  }
  
  if(d){
    wh_SIG_ge2b->Add(h_SIG_fDP_ge2b_njets_data, 2);
    wh_SIG_ge2b->Add(h_SIG_fDP_ge2b_njets_nonqcd, -2);
    wh_SIG_ge2b->Add(h_SIG_fDP_ge2b_njets_qcd, -1);
    wh_SIG_ge2b->Divide(h_SIG_fDP_ge2b_njets_qcd);
  }
  else{
    wh_SIG_ge2b->Add(h_SIG_fDP_ge2b_njets_data);
    wh_SIG_ge2b->Add(h_SIG_fDP_ge2b_njets_nonqcd, -1);
    wh_SIG_ge2b->Divide(h_SIG_fDP_ge2b_njets_qcd);
  }
  for(int i=1; i<= wh_SIG_ge2b->GetNbinsX(); i++){
    assert(h_SIG_fDP_ge2b_njets_data->GetBinContent(i)>0);
    assert(h_SIG_fDP_ge2b_njets_qcd->GetBinContent(i)>0);
    if(wh_SIG_ge2b->GetBinContent(i)<0) wh_SIG_ge2b->SetBinContent(i,0);
  }
  

  //make weighted histograms  
  h_LSB_fDP_antib_njets_qcd_w->Add(h_LSB_fDP_antib_njets_qcd);
  h_LSB_fDP_antib_njets_qcd_w->Multiply(wh_LSB_antib);
  h_LSB_pDP_antib_njets_qcd_w->Add(h_LSB_pDP_antib_njets_qcd);
  h_LSB_pDP_antib_njets_qcd_w->Multiply(wh_LSB_antib);
  
  h_SB_fDP_ge1b_njets_qcd_w->Add(h_SB_fDP_ge1b_njets_qcd);
  h_SB_fDP_ge1b_njets_qcd_w->Multiply(wh_SB_ge1b);
  h_SB_pDP_ge1b_njets_qcd_w->Add(h_SB_pDP_ge1b_njets_qcd);
  h_SB_pDP_ge1b_njets_qcd_w->Multiply(wh_SB_ge1b);
  
  h_SIG_fDP_ge1b_njets_qcd_w->Add(h_SIG_fDP_ge1b_njets_qcd);
  h_SIG_fDP_ge1b_njets_qcd_w->Multiply(wh_SIG_ge1b);
  h_SIG_pDP_ge1b_njets_qcd_w->Add(h_SIG_pDP_ge1b_njets_qcd);
  h_SIG_pDP_ge1b_njets_qcd_w->Multiply(wh_SIG_ge1b);
  
  h_SB_fDP_ge2b_njets_qcd_w->Add(h_SB_fDP_ge2b_njets_qcd);
  h_SB_fDP_ge2b_njets_qcd_w->Multiply(wh_SB_ge2b);
  h_SB_pDP_ge2b_njets_qcd_w->Add(h_SB_pDP_ge2b_njets_qcd);
  h_SB_pDP_ge2b_njets_qcd_w->Multiply(wh_SB_ge2b);
  
  h_SIG_fDP_ge2b_njets_qcd_w->Add(h_SIG_fDP_ge2b_njets_qcd);
  h_SIG_fDP_ge2b_njets_qcd_w->Multiply(wh_SIG_ge2b);
  h_SIG_pDP_ge2b_njets_qcd_w->Add(h_SIG_pDP_ge2b_njets_qcd);
  h_SIG_pDP_ge2b_njets_qcd_w->Multiply(wh_SIG_ge2b);
  

  
  //do ABCD
  double A, Aerr;
  double B, Berr;
  double C, Cerr;
  double D, Derr;
  
  cout << "ge1b, SB: unweighted" << endl;
  A = h_LSB_fDP_antib_njets_qcd->Integral();
  B = h_LSB_pDP_antib_njets_qcd->Integral();
  C = h_SB_pDP_ge1b_njets_qcd->Integral();
  D = h_SB_fDP_ge1b_njets_qcd->Integral();
  Aerr = jmt::errOnIntegral(h_LSB_fDP_antib_njets_qcd);
  Berr = jmt::errOnIntegral(h_LSB_pDP_antib_njets_qcd);
  Cerr = jmt::errOnIntegral(h_SB_pDP_ge1b_njets_qcd);
  Derr = jmt::errOnIntegral(h_SB_fDP_ge1b_njets_qcd);
  printABCD(A, Aerr, B, Berr, C, Cerr, D, Derr);
  cout << endl;
  
  cout << "ge1b, SB: weighted" << endl;
  A = h_LSB_fDP_antib_njets_qcd_w->Integral();
  B = h_LSB_pDP_antib_njets_qcd_w->Integral();
  C = h_SB_pDP_ge1b_njets_qcd_w->Integral();
  D = h_SB_fDP_ge1b_njets_qcd_w->Integral();
  Aerr = jmt::errOnIntegral(h_LSB_fDP_antib_njets_qcd_w);
  Berr = jmt::errOnIntegral(h_LSB_pDP_antib_njets_qcd_w);
  Cerr = jmt::errOnIntegral(h_SB_pDP_ge1b_njets_qcd_w);
  Derr = jmt::errOnIntegral(h_SB_fDP_ge1b_njets_qcd_w);
  printABCD(A, Aerr, B, Berr, C, Cerr, D, Derr);
  cout << endl;

  cout << "ge1b, SIG: unweighted" << endl;
  A = h_LSB_fDP_antib_njets_qcd->Integral();
  B = h_LSB_pDP_antib_njets_qcd->Integral();
  C = h_SIG_pDP_ge1b_njets_qcd->Integral();
  D = h_SIG_fDP_ge1b_njets_qcd->Integral();
  Aerr = jmt::errOnIntegral(h_LSB_fDP_antib_njets_qcd);
  Berr = jmt::errOnIntegral(h_LSB_pDP_antib_njets_qcd);
  Cerr = jmt::errOnIntegral(h_SIG_pDP_ge1b_njets_qcd);
  Derr = jmt::errOnIntegral(h_SIG_fDP_ge1b_njets_qcd);
  printABCD(A, Aerr, B, Berr, C, Cerr, D, Derr);
  cout << endl;

  cout << "ge1b, SIG: weighted" << endl;
  A = h_LSB_fDP_antib_njets_qcd_w->Integral();
  B = h_LSB_pDP_antib_njets_qcd_w->Integral();
  C = h_SIG_pDP_ge1b_njets_qcd_w->Integral();
  D = h_SIG_fDP_ge1b_njets_qcd_w->Integral();
  Aerr = jmt::errOnIntegral(h_LSB_fDP_antib_njets_qcd_w);
  Berr = jmt::errOnIntegral(h_LSB_pDP_antib_njets_qcd_w);
  Cerr = jmt::errOnIntegral(h_SIG_pDP_ge1b_njets_qcd_w);
  Derr = jmt::errOnIntegral(h_SIG_fDP_ge1b_njets_qcd_w);
  printABCD(A, Aerr, B, Berr, C, Cerr, D, Derr);  
  cout << endl;

  cout << "ge2b, SB: unweighted" << endl;
  A = h_LSB_fDP_antib_njets_qcd->Integral();
  B = h_LSB_pDP_antib_njets_qcd->Integral();
  C = h_SB_pDP_ge2b_njets_qcd->Integral();
  D = h_SB_fDP_ge2b_njets_qcd->Integral();
  Aerr = jmt::errOnIntegral(h_LSB_fDP_antib_njets_qcd);
  Berr = jmt::errOnIntegral(h_LSB_pDP_antib_njets_qcd);
  Cerr = jmt::errOnIntegral(h_SB_pDP_ge2b_njets_qcd);
  Derr = jmt::errOnIntegral(h_SB_fDP_ge2b_njets_qcd);
  printABCD(A, Aerr, B, Berr, C, Cerr, D, Derr);
  cout << endl;

  cout << "ge2b, SB: weighted" << endl;
  A = h_LSB_fDP_antib_njets_qcd_w->Integral();
  B = h_LSB_pDP_antib_njets_qcd_w->Integral();
  C = h_SB_pDP_ge2b_njets_qcd_w->Integral();
  D = h_SB_fDP_ge2b_njets_qcd_w->Integral();
  Aerr = jmt::errOnIntegral(h_LSB_fDP_antib_njets_qcd_w);
  Berr = jmt::errOnIntegral(h_LSB_pDP_antib_njets_qcd_w);
  Cerr = jmt::errOnIntegral(h_SB_pDP_ge2b_njets_qcd_w);
  Derr = jmt::errOnIntegral(h_SB_fDP_ge2b_njets_qcd_w);
  printABCD(A, Aerr, B, Berr, C, Cerr, D, Derr);
  cout << endl;

  cout << "ge2b, SIG: unweighted" << endl;
  A = h_LSB_fDP_antib_njets_qcd->Integral();
  B = h_LSB_pDP_antib_njets_qcd->Integral();
  C = h_SIG_pDP_ge2b_njets_qcd->Integral();
  D = h_SIG_fDP_ge2b_njets_qcd->Integral();
  Aerr = jmt::errOnIntegral(h_LSB_fDP_antib_njets_qcd);
  Berr = jmt::errOnIntegral(h_LSB_pDP_antib_njets_qcd);
  Cerr = jmt::errOnIntegral(h_SIG_pDP_ge2b_njets_qcd);
  Derr = jmt::errOnIntegral(h_SIG_fDP_ge2b_njets_qcd);
  printABCD(A, Aerr, B, Berr, C, Cerr, D, Derr);
  cout << endl;
  
  cout << "ge2b, SIG: weighted" << endl;
  A = h_LSB_fDP_antib_njets_qcd_w->Integral();
  B = h_LSB_pDP_antib_njets_qcd_w->Integral();
  C = h_SIG_pDP_ge2b_njets_qcd_w->Integral();
  D = h_SIG_fDP_ge2b_njets_qcd_w->Integral();
  Aerr = jmt::errOnIntegral(h_LSB_fDP_antib_njets_qcd_w);
  Berr = jmt::errOnIntegral(h_LSB_pDP_antib_njets_qcd_w);
  Cerr = jmt::errOnIntegral(h_SIG_pDP_ge2b_njets_qcd_w);
  Derr = jmt::errOnIntegral(h_SIG_fDP_ge2b_njets_qcd_w);
  printABCD(A, Aerr, B, Berr, C, Cerr, D, Derr);
  cout << endl;
  
  //
  TCanvas* myC = new TCanvas("myC", "myC", 600, 900);
  myC->Divide(2,3);
  myC->cd(1);    
  wh_LSB_antib->Draw();
  myC->cd(3);
  wh_SB_ge1b->Draw();
  myC->cd(4);
  wh_SIG_ge1b->Draw();
  myC->cd(5);
  wh_SB_ge2b->Draw();
  myC->cd(6);
  wh_SIG_ge2b->Draw();
  

  cout << "reweightQCD() is done." << endl;
}

