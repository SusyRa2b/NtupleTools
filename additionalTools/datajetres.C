#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"

#include <iostream>

using namespace std;

void datajetres(){

  //Data are based on reading off red line of Fig. 6 and 7 of JME-10-14 by eye (with a little help from Power Point's ruler...).
  //-- Be careful using this. It's not official and there are no uncertainties
  //-- Assumed constant extrapolation for pT<50 GeV and pT>end of curve. 

  TFile fout("datajetres.root", "RECREATE");
  
  const double bins1[] = {1,50,60,70,80,90,100,130,200,300,400,500,600,700,800};
  const double bins2[] = {1,50,60,70,80,90,100,130,200,300,400,500,600};
  const double bins3[] = {1,50,60,70,80,90,100,130,200,300};
  const double bins4[] = {1,50,60,70,80,90,100};
  
  const int nbins1 = 14;
  const int nbins2 = 12;
  const int nbins3 = 9;
  const int nbins4 = 6;

  TH1D* hEta0_0p5 = new TH1D("hEta0_0p5","0 < |eta| <= 0.5", nbins1, bins1);
  TH1D* hEta0p5_1 = new TH1D("hEta0p5_1","0.5 < |eta| <= 1", nbins1, bins1);
  TH1D* hEta1_1p5 = new TH1D("hEta1_1p5","1 < |eta| <= 1.5", nbins1, bins1);
  TH1D* hEta1p5_2 = new TH1D("hEta1p5_2","1.5 < |eta| <= 2", nbins1, bins1);
  TH1D* hEta2_2p5 = new TH1D("hEta2_2p5","2 < |eta| <= 2.5", nbins2, bins2);
  TH1D* hEta2p5_3 = new TH1D("hEta2p5_3","2.5 < |eta| <= 3", nbins3, bins3);
  TH1D* hEta3_5 = new TH1D("hEta3_5","3 < |eta| <= 5", nbins4, bins4);
  
  hEta0_0p5->SetBinContent(1, 0.122);
  hEta0_0p5->SetBinContent(2, 0.118);
  hEta0_0p5->SetBinContent(3, 0.107);
  hEta0_0p5->SetBinContent(4, 0.101);
  hEta0_0p5->SetBinContent(5, 0.0969);
  hEta0_0p5->SetBinContent(6, 0.0938);
  hEta0_0p5->SetBinContent(7, 0.0875);
  hEta0_0p5->SetBinContent(8, 0.0806);
  hEta0_0p5->SetBinContent(9, 0.0744);
  hEta0_0p5->SetBinContent(10, 0.0688);
  hEta0_0p5->SetBinContent(11, 0.0656);
  hEta0_0p5->SetBinContent(12, 0.0593);
  hEta0_0p5->SetBinContent(13, 0.0575);
  hEta0_0p5->SetBinContent(14, 0.0563);
  hEta0_0p5->Fill(1000, 0.0550);//overflow
  
  hEta0p5_1->SetBinContent(1, 0.122);
  hEta0p5_1->SetBinContent(2, 0.116);
  hEta0p5_1->SetBinContent(3, 0.108);
  hEta0p5_1->SetBinContent(4, 0.101);
  hEta0p5_1->SetBinContent(5, 0.0975);
  hEta0p5_1->SetBinContent(6, 0.0938);
  hEta0p5_1->SetBinContent(7, 0.0875);
  hEta0p5_1->SetBinContent(8, 0.0800);
  hEta0p5_1->SetBinContent(9, 0.0725);
  hEta0p5_1->SetBinContent(10, 0.0669);
  hEta0p5_1->SetBinContent(11, 0.0625);
  hEta0p5_1->SetBinContent(12, 0.0606);
  hEta0p5_1->SetBinContent(13, 0.0575);
  hEta0p5_1->SetBinContent(14, 0.0563);
  hEta0p5_1->Fill(1000, 0.0556);//overflow
  
  hEta1_1p5->SetBinContent(1, 0.141);
  hEta1_1p5->SetBinContent(2, 0.134);
  hEta1_1p5->SetBinContent(3, 0.124);
  hEta1_1p5->SetBinContent(4, 0.116);
  hEta1_1p5->SetBinContent(5, 0.111);
  hEta1_1p5->SetBinContent(6, 0.107);
  hEta1_1p5->SetBinContent(7, 0.100);
  hEta1_1p5->SetBinContent(8, 0.0900);
  hEta1_1p5->SetBinContent(9, 0.0831);
  hEta1_1p5->SetBinContent(10, 0.0775);
  hEta1_1p5->SetBinContent(11, 0.0744);
  hEta1_1p5->SetBinContent(12, 0.0713);
  hEta1_1p5->SetBinContent(13, 0.0700);
  hEta1_1p5->SetBinContent(14, 0.0688);
  hEta1_1p5->Fill(1000, 0.0675);//overflow

  hEta1p5_2->SetBinContent(1, 0.135);
  hEta1p5_2->SetBinContent(2, 0.127);
  hEta1p5_2->SetBinContent(3, 0.117);
  hEta1p5_2->SetBinContent(4, 0.110);
  hEta1p5_2->SetBinContent(5, 0.104);
  hEta1p5_2->SetBinContent(6, 0.0994);
  hEta1p5_2->SetBinContent(7, 0.0913);
  hEta1p5_2->SetBinContent(8, 0.0831);
  hEta1p5_2->SetBinContent(9, 0.0750);
  hEta1p5_2->SetBinContent(10, 0.0700);
  hEta1p5_2->SetBinContent(11, 0.0681);
  hEta1p5_2->SetBinContent(12, 0.0650);
  hEta1p5_2->SetBinContent(13, 0.0638);
  hEta1p5_2->SetBinContent(14, 0.0625);
  hEta1p5_2->Fill(1000, 0.0619);//overflow
  
  hEta2_2p5->SetBinContent(1, 0.125);
  hEta2_2p5->SetBinContent(2, 0.108);
  hEta2_2p5->SetBinContent(3, 0.100);
  hEta2_2p5->SetBinContent(4, 0.0956);
  hEta2_2p5->SetBinContent(5, 0.0925);
  hEta2_2p5->SetBinContent(6, 0.0889);
  hEta2_2p5->SetBinContent(7, 0.0838);
  hEta2_2p5->SetBinContent(8, 0.0763);
  hEta2_2p5->SetBinContent(9, 0.0713);
  hEta2_2p5->SetBinContent(10, 0.0675);
  hEta2_2p5->SetBinContent(11, 0.0644);
  hEta2_2p5->SetBinContent(12, 0.0625);
  hEta2_2p5->Fill(1000, 0.0619);//overflow
  
  hEta2p5_3->SetBinContent(1, 0.124);
  hEta2p5_3->SetBinContent(2, 0.122);
  hEta2p5_3->SetBinContent(3, 0.119);
  hEta2p5_3->SetBinContent(4, 0.114);
  hEta2p5_3->SetBinContent(5, 0.112);
  hEta2p5_3->SetBinContent(6, 0.109);
  hEta2p5_3->SetBinContent(7, 0.104);
  hEta2p5_3->SetBinContent(8, 0.0969);
  hEta2p5_3->SetBinContent(9, 0.0900);
  hEta2p5_3->Fill(1000, 0.0863);//overflow

  hEta3_5->SetBinContent(1, 0.136);
  hEta3_5->SetBinContent(2, 0.133);
  hEta3_5->SetBinContent(3, 0.130);
  hEta3_5->SetBinContent(4, 0.126);
  hEta3_5->SetBinContent(5, 0.125);
  hEta3_5->SetBinContent(6, 0.124);
  hEta3_5->Fill(1000, 0.122);//overflow
  
  bool drawdatajetres = true;
  if(drawdatajetres){
    TCanvas * myC = new TCanvas("myC", "myC", 900,500);
    myC->Divide(4,2);
    TString xtitle = "p_{T} [GeV]";
    TString ytitle = "jet p_{T} resolution";
    hEta0_0p5->GetXaxis()->SetTitle(xtitle);
    hEta0p5_1->GetXaxis()->SetTitle(xtitle);
    hEta1_1p5->GetXaxis()->SetTitle(xtitle);
    hEta1p5_2->GetXaxis()->SetTitle(xtitle);
    hEta2_2p5->GetXaxis()->SetTitle(xtitle);
    hEta2p5_3->GetXaxis()->SetTitle(xtitle);
    hEta3_5->GetXaxis()->SetTitle(xtitle);

    hEta0_0p5->GetYaxis()->SetTitle(ytitle);
    hEta0p5_1->GetYaxis()->SetTitle(ytitle);
    hEta1_1p5->GetYaxis()->SetTitle(ytitle);
    hEta1p5_2->GetYaxis()->SetTitle(ytitle);
    hEta2_2p5->GetYaxis()->SetTitle(ytitle);
    hEta2p5_3->GetYaxis()->SetTitle(ytitle);
    hEta3_5->GetYaxis()->SetTitle(ytitle);

    myC->cd(1); //myC->GetPad(1)->SetLogx();
    //hEta0_0p5->GetXaxis()->SetRangeUser(30,hEta0_0p5->GetBinLowEdge(hEta0_0p5->GetNbinsX())); //doesnt work...
    hEta0_0p5->Draw();
    myC->cd(2); //myC->GetPad(2)->SetLogx();
    hEta0p5_1->Draw();
    myC->cd(3); //myC->GetPad(3)->SetLogx();
    hEta1_1p5->Draw();
    myC->cd(4); //myC->GetPad(4)->SetLogx();
    hEta1p5_2->Draw();
    myC->cd(5); //myC->GetPad(5)->SetLogx();
    hEta2_2p5->Draw();
    myC->cd(6); //myC->GetPad(6)->SetLogx();
    hEta2p5_3->Draw();
    myC->cd(7); //myC->GetPad(7)->SetLogx();
    hEta3_5->Draw();
    myC->Print("djr.png");
  }


  hEta0_0p5->Write();
  hEta0p5_1->Write();
  hEta1_1p5->Write();
  hEta1p5_2->Write();
  hEta2_2p5->Write();
  hEta2p5_3->Write();
  hEta3_5->Write();
  fout.Close();
}


