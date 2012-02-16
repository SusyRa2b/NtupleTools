#include "TFile.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TCanvas.h"

#include <iostream>

using namespace std;

TRandom3* r = 0; 

void runlsbpvrw(TString name = "ge1bLoose") {

  TFile* fin = TFile::Open("../pvrw.root","READ");
  if(fin->IsZombie()) cout << "Problem opening input file!" << endl;
  
  TH1D* hPVphysics = (TH1D*)fin->Get("hPVphysics_"+name);
  TH1D* hPVprescalePass = (TH1D*)fin->Get("hPVprescalePass_"+name);
  TH1D* hPVprescaleFail = (TH1D*)fin->Get("hPVprescaleFail_"+name);
  
  TH1D* tPVphysics = (TH1D*)hPVphysics->Clone("tPVphysics");
  TH1D* tPVprescalePass = (TH1D*)hPVprescalePass->Clone("tPVprescalePass");
  TH1D* tPVprescaleFail = (TH1D*)hPVprescaleFail->Clone("tPVprescaleFail");
  tPVphysics->Reset();
  tPVprescalePass->Reset();
  tPVprescaleFail->Reset();
  
  TH1D* hPFresults = new TH1D("hPFresults", "hPFresults", 500, 0.0, 0.3);


  //toy loop
  for(int i=0; i<1000; i++) {
    
    //make toy inputs
    for(int n=1; n<=hPVphysics->GetNbinsX(); n++) {
      //tPVphysics->SetBinContent(n,r->Gaus(hPVphysics->GetBinContent(n), hPVphysics->GetBinError(n)));
      //tPVprescalePass->SetBinContent(n,r->Gaus(hPVprescalePass->GetBinContent(n),hPVprescalePass->GetBinError(n)));
      //tPVprescaleFail->SetBinContent(n,r->Gaus(hPVprescaleFail->GetBinContent(n),hPVprescaleFail->GetBinError(n)));

      tPVphysics->SetBinContent(n,r->Poisson(hPVphysics->GetBinContent(n)));
      tPVprescalePass->SetBinContent(n,r->Poisson(hPVprescalePass->GetBinContent(n)));
      tPVprescaleFail->SetBinContent(n,r->Poisson(hPVprescaleFail->GetBinContent(n)));
    }

    
    //do cross check (no toy)
    if(i==0) {
      TH1D* hPVprescale = (TH1D*)hPVprescalePass->Clone("hPVprescale");
      hPVprescale->Reset();
      hPVprescale->Add(hPVprescalePass,hPVprescaleFail);

      TH1D* hPVprescale_RW = (TH1D*)hPVphysics->Clone("hPVprescale_RW");
      hPVprescale_RW->Scale(hPVprescale->Integral()/hPVphysics->Integral());
      TH1D* hPV_W = (TH1D*)hPVphysics->Clone("hPV_W");
      hPV_W->Reset();
      hPV_W->Divide(hPVprescale_RW,hPVprescale); 
      
      //weighted LSB pass and fail mdp
      TH1D* hPVprescalePass_RW = (TH1D*)hPVprescalePass->Clone("hPVprescalePass_RW");
      TH1D* hPVprescaleFail_RW = (TH1D*)hPVprescaleFail->Clone("hPVprescaleFail_RW");
      hPVprescalePass_RW->Multiply(hPV_W);
      hPVprescaleFail_RW->Multiply(hPV_W);
    
      cout << name << " PFratio: " << hPVprescalePass_RW->Integral()/hPVprescaleFail_RW->Integral() << endl;
    }


    //now find PFratio for this toy
    TH1D* tPVprescale = (TH1D*)tPVprescalePass->Clone("tPVprescale");
    tPVprescale->Reset();
    tPVprescale->Add(tPVprescalePass,tPVprescaleFail);
    
    TH1D* tPVprescale_RW = (TH1D*)tPVphysics->Clone("tPVprescale_RW");
    tPVprescale_RW->Scale(tPVprescale->Integral()/tPVphysics->Integral());
    TH1D* tPV_W = (TH1D*)tPVphysics->Clone("tPV_W");
    tPV_W->Reset();
    tPV_W->Divide(tPVprescale_RW,tPVprescale); 
    
    //weighted LSB pass and fail mdp
    TH1D* tPVprescalePass_RW = (TH1D*)tPVprescalePass->Clone("tPVprescalePass_RW");
    TH1D* tPVprescaleFail_RW = (TH1D*)tPVprescaleFail->Clone("tPVprescaleFail_RW");
    tPVprescalePass_RW->Multiply(tPV_W);
    tPVprescaleFail_RW->Multiply(tPV_W);
    
    double tPFratio = tPVprescalePass_RW->Integral()/tPVprescaleFail_RW->Integral();
      
    hPFresults->Fill(tPFratio);

  }//end of toy loop
  
  cout << name << " result: " << hPFresults->GetMean() << " +- " << hPFresults->GetRMS() << endl;

  TCanvas* c1 = new TCanvas("c1", "c1", 640, 480);
  c1->cd();
  hPFresults->Draw();  
  c1->SaveAs("lsbpvrw_"+name+".png");

  fin->Close();
}

void lsbpvrw(){

  r = new TRandom3(4444);
  
  runlsbpvrw("ge1bLoose");
  runlsbpvrw("ge1bTight");
  runlsbpvrw("ge2bLoose");
  runlsbpvrw("ge2bTight");
  //runlsbpvrw("ge3bLoose");

}
