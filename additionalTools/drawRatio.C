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


void drawRatio(){

  //input from running drawSimple for pass and for fail with both histograms going into "ben.root"
  //run by doing root -l drawRatio.C++

  gROOT->SetStyle("CMS");

  TString HT = "500";
  TString btag = "eq0b";
  TString var = "minDeltaPhiN";

  TString name1 = var+"_HT"+HT+"_"+btag+"_pass";
  TString name2 = var+"_HT"+HT+"_"+btag+"_fail";

  TFile *f1=TFile::Open("../ben.root","READ");

  TCanvas * myC = new TCanvas("myC", "myC", 800, 480); 
  myC->SetLogy();
  TH1D* h1=(TH1D*)f1->Get(name1);
  TH1D* h2=(TH1D*)f1->Get(name2);
  h1->UseCurrentStyle();
  h2->UseCurrentStyle();
  
  TH1D* h3 = (TH1D*)h1->Clone("h3");
  h3->Reset();
  h3->Divide(h1,h2);
  h3->Draw();
  

  h1->SetLineColor(40);
  h1->SetMarkerColor(40);
  
  h2->SetLineColor(30);
  h2->SetMarkerColor(30);
 
  h3->SetLineColor(kBlack);
  h3->SetMarkerColor(kBlack);
  h3->SetLineWidth(2);


  //h2->GetYaxis()->SetTitle("Events (blue, green), Ratio (black)");
  h2->GetXaxis()->SetTitle("MET");
  h3->GetXaxis()->SetTitle("MET");  
  h3->GetYaxis()->SetTitle("N pass / N fail");
  h2->SetMinimum(0.0009);
  
  h2->Draw("hist e");
  h1->Draw("SAME hist e");  
  h3->Draw("SAME hist e");

  TLegend *leg = new TLegend(0.55,.7,0.88,0.88);
  leg->AddEntry(h1,"pass MDP cut [Events]", "P");
  leg->AddEntry(h2,"fail MDP cut [Events]", "P");
  leg->AddEntry(h3,"pass/fail [Ratio]", "P");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetLineStyle(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->Draw();

  gPad->SetRightMargin(0.1);
  gPad->Modified();
  myC->Print(var+"_HT"+HT+"_"+btag+".png");
  
  TCanvas * myC2 = new TCanvas("myC2", "myC2", 800, 480);
  myC2->cd();
  h3->Draw("hist e");
  h3->SetMaximum(1.5);
  gPad->SetRightMargin(0.1);
  gPad->Modified();
  myC2->Print(var+"_HT"+HT+"_"+btag+"_ratio.png");
  
}
