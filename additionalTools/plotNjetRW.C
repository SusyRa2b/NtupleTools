#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TPad.h"
#include <iostream>

using namespace std;

void label(TH1D* h) {
  h->GetYaxis()->SetTitle("Data/MC Scale Factor");
  h->GetXaxis()->SetTitle("Jet Multiplicity (pT>30GeV)");
  return;
}

void plotNjetRW() {
  gStyle->SetOptStat(0);

  TFile* f1BLsb = new TFile("../njetRW_ge1bLooseWideSB_SB.root", "READ");
  TFile* f1BLsig = new TFile("../njetRW_ge1bLooseWideSB_SIG.root", "READ");

  TFile* f1BTsb = new TFile("../njetRW_ge1bTightWideSB_SB.root", "READ");
  TFile* f1BTsig = new TFile("../njetRW_ge1bTightWideSB_SIG.root", "READ");

  
  TFile* f2BLsb = new TFile("../njetRW_ge2bLooseWideSB_SB.root", "READ");
  TFile* f2BLsig = new TFile("../njetRW_ge2bLooseWideSB_SIG.root", "READ");
  TFile* f2BTsb = new TFile("../njetRW_ge2bTightWideSB_SB.root", "READ");
  TFile* f2BTsig = new TFile("../njetRW_ge2bTightWideSB_SIG.root", "READ");
  
  /*
  TFile* f3BLsb = new TFile("../njetRW_ge3bLooseWideSB_SB.root", "READ");
  TFile* f3BLsig = new TFile("../njetRW_ge3bLooseWideSB_SIG.root", "READ");
  */

  TH1D* h1BLsb_pres = (TH1D*)f1BLsb->Get("hpresW");
  TH1D* h1BLsb_phys = (TH1D*)f1BLsb->Get("hphysW");
  TH1D* h1BLsig_pres = (TH1D*)f1BLsig->Get("hpresW");
  TH1D* h1BLsig_phys = (TH1D*)f1BLsig->Get("hphysW");

  TH1D* h1BTsb_pres = (TH1D*)f1BTsb->Get("hpresW");
  TH1D* h1BTsb_phys = (TH1D*)f1BTsb->Get("hphysW");
  TH1D* h1BTsig_pres = (TH1D*)f1BTsig->Get("hpresW");
  TH1D* h1BTsig_phys = (TH1D*)f1BTsig->Get("hphysW");

  
  TH1D* h2BLsb_pres = (TH1D*)f2BLsb->Get("hpresW");
  TH1D* h2BLsb_phys = (TH1D*)f2BLsb->Get("hphysW");
  TH1D* h2BLsig_pres = (TH1D*)f2BLsig->Get("hpresW");
  TH1D* h2BLsig_phys = (TH1D*)f2BLsig->Get("hphysW");

  TH1D* h2BTsb_pres = (TH1D*)f2BTsb->Get("hpresW");
  TH1D* h2BTsb_phys = (TH1D*)f2BTsb->Get("hphysW");
  TH1D* h2BTsig_pres = (TH1D*)f2BTsig->Get("hpresW");
  TH1D* h2BTsig_phys = (TH1D*)f2BTsig->Get("hphysW");
  
  /*
  TH1D* h3BLsb_pres = (TH1D*)f3BLsb->Get("hpresW");
  TH1D* h3BLsb_phys = (TH1D*)f3BLsb->Get("hphysW");
  TH1D* h3BLsig_pres = (TH1D*)f3BLsig->Get("hpresW");
  TH1D* h3BLsig_phys = (TH1D*)f3BLsig->Get("hphysW");
  */


  //SB and SIG use the SAME factors!
  //HT cut follows selection
  //Factors are independent of Btag

  //so 1BL, 2BL, and 3B use the same. 1BT and 2BT are different.
  //just need to show plots for SB or SIG alone.

  label(h1BLsb_pres);
  label(h1BLsb_phys);
  TCanvas* c1BL = new TCanvas("1BL", "1BL", 2*640,480);
  c1BL->Divide(2,1);
  c1BL->cd(1);
  gPad->SetLeftMargin(.2);
  gPad->SetBottomMargin(.2);
  h1BLsb_pres->Draw();
  c1BL->cd(2);
  gPad->SetLeftMargin(.2);
  gPad->SetBottomMargin(.2);
  h1BLsb_phys->Draw();
  c1BL->Print("c1BL.png");

  label(h1BTsb_pres);
  label(h1BTsb_phys);
  TCanvas* c1BT = new TCanvas("1BT", "1BT", 2*640,480);
  c1BT->Divide(2,1);
  c1BT->cd(1);
  gPad->SetLeftMargin(.2);
  gPad->SetBottomMargin(.2);
  h1BTsb_pres->Draw();
  c1BT->cd(2);
  gPad->SetLeftMargin(.2);
  gPad->SetBottomMargin(.2);
  h1BTsb_phys->Draw();
  c1BT->Print("c1BT.png");
  
  label(h2BTsb_pres);
  label(h2BTsb_phys);
  TCanvas* c2BT = new TCanvas("2BT", "2BT", 2*640,480);
  c2BT->Divide(2,1);
  c2BT->cd(1);
  gPad->SetLeftMargin(.2);
  gPad->SetBottomMargin(.2);
  h2BTsb_pres->Draw();
  c2BT->cd(2);
  gPad->SetLeftMargin(.2);
  gPad->SetBottomMargin(.2);
  h2BTsb_phys->Draw();
  c2BT->Print("c2BT.png");
  
  return;


}
