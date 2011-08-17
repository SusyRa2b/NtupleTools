#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TH1D.h"
#include <iostream>

using namespace std;

void qcdWeights(){
  //gROOT->SetStyle("CMS");
  TFile* ifl = TFile::Open("/cu2/ra2b/reducedTrees/V00-02-24_fullpf2pat/reducedTree.SSVHPT_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.PythiaPUQCD.root","READ"); 
  TTree* reducedTree = (TTree*)ifl->Get("reducedTree");

  TCanvas* myC = new TCanvas("myC", "myC", 800,800);
  myC->Divide(2,2);
  myC->GetPad(1)->SetLogy();
  myC->GetPad(1)->SetLogx();
  myC->GetPad(2)->SetLogy();
  myC->GetPad(2)->SetLogx();
  myC->GetPad(3)->SetLogy();
  myC->GetPad(3)->SetLogx();
  myC->GetPad(4)->SetLogy();
  myC->GetPad(4)->SetLogx();

  TCut base = "cutPV==1 && cutTrigger==1  && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && minDeltaPhiN >= 4";
  TCut loose = "HT>=350&&MET>=200";
  TCut tight = "HT>=500&&MET>=300";
  TCut ge1b = "nbjetsSSVHPT>=1";
  TCut ge2b = "nbjetsSSVHPT>=2";
  
  TString selection_ = "";

  /*
  //for weight histograms to look at
  TH1D* myH1 = new TH1D("weights1","loose, ge1b signal region", 1000000,0,.01);
  TH1D* myH2 = new TH1D("weights2","loose, ge2b signal region", 1000000,0,.01);
  TH1D* myH3 = new TH1D("weights3","tight, ge1b signal region", 1000000,0,.01);
  TH1D* myH4 = new TH1D("weights4","tight, ge2b signal region", 1000000,0,.01);
  */

  /*
  //for counting high weight events
  TH1D* myH1 = new TH1D("weights1","loose, ge1b signal region", 1,0.0008,.01);
  TH1D* myH2 = new TH1D("weights2","loose, ge2b signal region", 1,0.0008,.01);
  TH1D* myH3 = new TH1D("weights3","tight, ge1b signal region", 1,0.0008,.01);
  TH1D* myH4 = new TH1D("weights4","tight, ge2b signal region", 1,0.0008,.01);
  */

  //look at high weight events
  TH1D* myH1 = new TH1D("weights1","loose, ge1b signal region", 1000,0.0008,.01);
  TH1D* myH2 = new TH1D("weights2","loose, ge2b signal region", 1000,0.0008,.01);
  TH1D* myH3 = new TH1D("weights3","tight, ge1b signal region", 1000,0.0008,.01);
  TH1D* myH4 = new TH1D("weights4","tight, ge2b signal region", 1000,0.0008,.01);

  myH1->GetYaxis()->SetRangeUser(.1,3000);
  myH2->GetYaxis()->SetRangeUser(.1,3000);
  myH3->GetYaxis()->SetRangeUser(.1,3000);
  myH4->GetYaxis()->SetRangeUser(.1,3000);
  myH1->GetYaxis()->SetTitle("raw number");
  myH2->GetYaxis()->SetTitle("raw number");
  myH3->GetYaxis()->SetTitle("raw number");
  myH4->GetYaxis()->SetTitle("raw number");
  myH1->GetXaxis()->SetTitle("weight for 1/pb");
  myH2->GetXaxis()->SetTitle("weight for 1/pb");
  myH3->GetXaxis()->SetTitle("weight for 1/pb");
  myH4->GetXaxis()->SetTitle("weight for 1/pb");
  myH1->SetLineWidth(2);
  myH2->SetLineWidth(2);
  myH3->SetLineWidth(2);  
  myH4->SetLineWidth(2);
   
  //cross check to check cuts
  TString temp = "";
  temp = base&&loose&&ge1b;
  selection_ = "("+temp+")*weight*1091.891";
  TH1D* hCheck = new TH1D("hCheck","hCheck",1,200,1E9);
  reducedTree->Project("hCheck","MET",selection_);
  cout << "loose, ge1b signal, weighted integral: " << hCheck->Integral() << endl;
  
  myC->cd(1);
  selection_ = base&&loose&&ge1b;
  reducedTree->Project("weights1","weight",selection_);
  myH1->Draw();
  
  myC->cd(2);
  selection_ = base&&loose&&ge2b;
  reducedTree->Project("weights2","weight",selection_);
  myH2->Draw();
  
  myC->cd(3);
  selection_ = base&&tight&&ge1b;
  reducedTree->Project("weights3","weight",selection_);
  myH3->Draw();
  
  myC->cd(4);
  selection_ = base&&tight&&ge2b;
  reducedTree->Project("weights4","weight",selection_);
  myH4->Draw();
  
  myC->Print("qcd_signal_weights.png");
  ifl->Close();
}

