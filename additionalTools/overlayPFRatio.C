#include "TPad.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLine.h"
#include "TFile.h"
#include "TLegend.h"
#include "TBox.h"
#include "TH1D.h"
#include "TH1D.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include <vector>
#include <iostream>
#include <fstream>


void overlayPFRatio(){
  // input is files created by drawR using the "hack to save" qcd only. it creates a ratio in its own root file, both with the same name.
  // run by doing root -l overlayPFRatio.C++

  gROOT->SetStyle("CMS");

  TString name1 = "minDeltaPhiN_HT350_ge1b_qcd";
  TString name2 = "minDeltaPhiN_DJR_HT350_ge1b_qcd";
  //TString name3 = "MET_15-2.7_"+bstring+"b_qcd";

  TFile *f1=TFile::Open("../"+name1+".root","READ");
  TFile *f2=TFile::Open("../"+name2+".root","READ");
  //TFile *f3=TFile::Open(name3+".root","READ");

  TH1D* h1=(TH1D*)f1->Get(name1);
  TH1D* h2=(TH1D*)f2->Get(name2);
  //TH1D* h3=(TH1D*)f3->Get(name3);
  h1->UseCurrentStyle();
  h2->UseCurrentStyle();
  //h3->UseCurrentStyle();
    
  h2->SetLineColor(kRed);
  h2->SetMarkerColor(kRed);
  //h3->SetLineColor(kBlue);
  //h3->SetMarkerColor(kBlue);

  h1->GetYaxis()->SetTitle("N pass / N fail");
  h1->GetXaxis()->SetTitle("MET");
  //h1->SetMaximum(1.);
  h1->SetMinimum(0.);
  
  h1->Draw("hist e");
  h2->Draw("SAME hist e");  
  //h3->Draw("SAME hist e");

  TLegend *leg = new TLegend(0.2,.7,0.6,0.88);
  leg->AddEntry(h1,name1, "P");
  leg->AddEntry(h2, name2, "P");
  //leg->AddEntry(h3, name3, "P");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetLineStyle(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);
  leg->Draw();

  gPad->SetRightMargin(0.1);
  gPad->Modified();

}
