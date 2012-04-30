/*
computeOverallEff()

This code just calculates the 2011 lumi-averaged MHT efficiency for
the MET range of interest.

One must specify the SB range with the sbType enum, and the kinematic
region of interest with the regionType enum.

*/

#include "TF2.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TFrame.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TF1.h"
#include "TList.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLatex.h"

#include "TGraphAsymmErrors.h"

//#include "tdrstyle.C" 

//#include "binomial_interval.C"


void computeOverallEff() {


  //gROOT->ProcessLine(".L tdrstyle.C");
  //setTDRStyle();
  //gStyle->SetHatchesSpacing(1.0);
  //gROOT->ForceStyle(); 
  ////gROOT->Reset();
  //gROOT->SetStyle("Plain");
  ////gStyle->SetOptFit(1011);
  //gStyle->SetOptStat(0);
  //gStyle->SetMarkerStyle(20);
  gStyle->SetPadGridY(true);
  gStyle->SetPadGridX(true);
  gStyle->SetPadRightMargin(0.05);
  


  bool is1BT = false;
  bool is2BT = false;
  bool isSBANDSIG = false;



  //which sideband?
  enum sbType {kNominal=0, kLow, kWide};// kNominal = 200-250, kLow = 150-200, kWide = 150-250
  enum regionType {kZeroLepHighMDP=0, kZeroLepLowMDP, kSingleLepHighMDP, kSingleMu, kSingleEle};
  //sbType theSBType_ = kNominal;
  sbType theSBType_ = kWide;
  //sbType theSBType_ = kLow;

  regionType theRegionType_ = kZeroLepHighMDP;
  //regionType theRegionType_ = kZeroLepLowMDP;
  //regionType theRegionType_ = kSingleLeppHighMDP; //not used -> use Keith's numbers
  //regionType theRegionType_ = kSingleMu;



  int sb_low, sb_high;
  if(theSBType_==kNominal)   {sb_low=200; sb_high=250;}
  else if(theSBType_==kLow)  {sb_low=150; sb_high=200;}
  else if(theSBType_==kWide) {sb_low=150; sb_high=250;}


  //TFile* f_trig1, f_trig2;
  TFile * f_trig1 = new TFile(     "mhteff_dec8_eq0b.root","READ");
  TFile * f_trig2 = new TFile(     "mhteff_dec8_eq0b.root","READ");
  if( theRegionType_ == kZeroLepHighMDP || theRegionType_ == kZeroLepLowMDP || theRegionType_ == kSingleLepHighMDP){
    f_trig1 = new TFile("mhteff_dec8_eq0b.root");
    f_trig2 = new TFile("mhteff_dec8_eq0b.root");
    //f_trig1 = new TFile("mhteff_apr18_ge1b.root");
    //f_trig2 = new TFile("mhteff_apr18_ge1b.root");
    // f_trig1 = new TFile("mhteff_apr18_eq0b_mindphin8.root");
    //f_trig2 = new TFile("mhteff_apr18_eq0b_mindphin8.root");
    //f_trig1 = new TFile("mhteff_apr18_ge1b_mindphin10.root");
    //f_trig2 = new TFile("mhteff_apr18_ge1b_mindphin10.root");
    //f_trig1 = new TFile("mhteff_apr18_ge1b_mindphin12.root");
    //f_trig2 = new TFile("mhteff_apr18_ge1b_mindphin12.root");
    //f_trig1 = new TFile("mhteff_apr18_ge1b_mindphin14.root");
    //f_trig2 = new TFile("mhteff_apr18_ge1b_mindphin14.root");
    //f_trig1 = new TFile("mhteff_apr18_eq0b_mindphin10.root");
    //f_trig2 = new TFile("mhteff_apr18_eq0b_mindphin10.root");
    //f_trig1 = new TFile("mhteff_apr18_eq0b_mindphin12.root");
    //f_trig2 = new TFile("mhteff_apr18_eq0b_mindphin12.root");
    //f_trig1 = new TFile("mhteff_apr18_eq0b_mindphin14.root");
    //f_trig2 = new TFile("mhteff_apr18_eq0b_mindphin14.root");
  }
  else if( theRegionType_ == kSingleMu || theRegionType_ == kSingleEle){
    f_trig1 = new TFile("1LTrig_v2.root");
    f_trig2 = new TFile("1LTrig_earlier.root");
  }

  //TFile * f_trig1 = new TFile(     "mhteff_nov29.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_dec3.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_dec4.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_dec4_ht500.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_dec4_ht600.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_dec8_eq0b.root","READ");
  //TFile * f_trig2 = new TFile(     "mhteff_dec8_eq0b.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_dec8_ge1b.root","READ");
  //TFile * f_trig2 = new TFile(     "mhteff_dec8_ge1b.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_apr18_ge1b.root","READ");
  //TFile * f_trig2 = new TFile(     "mhteff_apr18_ge1b.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_apr18_ge1b_with3jetcut.root","READ");
  //TFile * f_trig2 = new TFile(     "mhteff_apr18_ge1b_with3jetcut.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_apr18_eq0b_with3jetcut.root","READ");
  //TFile * f_trig2 = new TFile(     "mhteff_apr18_eq0b_with3jetcut.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_apr18_eq0b_mindphin6.root","READ");
  //TFile * f_trig2 = new TFile(     "mhteff_apr18_eq0b_mindphin6.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_apr18_ge1b_mindphin6.root","READ");
  //TFile * f_trig2 = new TFile(     "mhteff_apr18_ge1b_mindphin6.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_apr18_eq0b_mindphin8.root","READ");
  //TFile * f_trig2 = new TFile(     "mhteff_apr18_eq0b_mindphin8.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_apr18_ge1b_mindphin8.root","READ");
  //TFile * f_trig2 = new TFile(     "mhteff_apr18_ge1b_mindphin8.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_apr18_eq0b_mindphin10.root","READ");
  //TFile * f_trig2 = new TFile(     "mhteff_apr18_eq0b_mindphin10.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_apr18_ge1b_mindphin10.root","READ");
  //TFile * f_trig2 = new TFile(     "mhteff_apr18_ge1b_mindphin10.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_apr18_ge2b.root","READ");
  //TFile * f_trig2 = new TFile(     "mhteff_apr18_ge2b.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_dec8_eq0b_typeIMET.root","READ");
  //TFile * f_trig2 = new TFile(     "mhteff_dec8_eq0b_typeIMET.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_dec8_eq0b_HT500.root","READ");
  //TFile * f_trig1 = new TFile(     "mhteff_dec8_eq0b_HT600.root","READ");
  //TFile * f_trig1 = new TFile(       "mhteff_nov29_ge1b.root","READ");
  //TFile * f_trig1 = new TFile(     "1LTrig_v2.root","READ");
  //TFile * f_trig2 = new TFile(     "1LTrig_earlier.root","READ");



  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;

  if(is1BT) std::cout << "IMPORTANT: YOU HAVE THE is1BT FLAG SET" << std::endl;
  if(is2BT) std::cout << "IMPORTANT: YOU HAVE THE is2BT FLAG SET" << std::endl;
  if(isSBANDSIG) std::cout << "IMPORTANT: YOU HAVE THE isSBANDSIG FLAG SET" << std::endl;

  if(theSBType_==kNominal)   std::cout << "SB range = 200-250"<< std::endl;
  else if(theSBType_==kLow)  std::cout << "SB range = 150-200"<< std::endl;
  else if(theSBType_==kWide) std::cout << "SB range = 150-250"<< std::endl;

  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;


  double totallumi = 0;
  double totallumi_SB = 0; //the total lumi is computed on the fly (to deal with chunks in which no efficiency measurement is possible)
  double totallumi_SIG = 0; //the total lumi is computed on the fly (to deal with chunks in which no efficiency measurement is possible)
  double totallumi_ignore_SB = 0;
  double totallumi_ignore_SIG = 0;

  std::vector<double> vlumi;
  std::vector<bool> vlumi_ignore_SB,vlumi_ignore_SIG;


  //pixelLumiCalc numbers (April 5)
  totallumi = 4.982; //this is exactly the sum of the below lumi numbers
  vlumi.push_back(0.0063);//ht260_mht60
  vlumi.push_back(0.0407+0.1687);//ht250_mht60
  vlumi.push_back(0.1364);//ht250_mht70
  vlumi.push_back(0.1005+0.0044+0.4366);//ht300_mht75
  vlumi.push_back(0.2773+0.8313);//ht300_mht80
  vlumi.push_back(0.6522);//ht300_mht90
  vlumi.push_back(1.4421);//ht350_mht95
  vlumi.push_back(0.8855);//ht350_mht110

  //pre-pixelLumiCalc
  //if(isSL) totallumi = 2.7574;
  //else 
  //totallumi = 4.6838; //this is exactly the sum of the below lumi numbers
  //if(!isSL){
  //vlumi.push_back(0.0067);//ht260_mht60
  //vlumi.push_back(0.0409+0.1686);//ht250_mht60
  //vlumi.push_back(0.1391);//ht250_mht70
  //vlumi.push_back(0.0982+0.0043+0.4225);//ht300_mht75
  //vlumi.push_back(0.2657+0.7804);//ht300_mht80
  ////}
  //vlumi.push_back(0.6044);//ht300_mht90
  //vlumi.push_back(1.3405);//ht350_mht95
  //vlumi.push_back(0.8125);//ht350_mht110


  std::vector<std::string> vhname_denom, vhname_num;
  
  ////HT>=400
  //vhname_denom.push_back("MET_ht260mht60v2_den");
  //vhname_denom.push_back("MET_ht250mht60v2and3_den");
  //vhname_denom.push_back("MET_ht250mht70v1_den");
  //vhname_denom.push_back("MET_ht300mht75v7and8_den");
  //vhname_denom.push_back("MET_ht300mht80v1and2_den");
  //vhname_denom.push_back("MET_ht300mht90v2_den");
  //vhname_denom.push_back("MET_ht350mht90v1_den");
  //vhname_denom.push_back("MET_ht350mht110v3_den");
  //vhname_num.push_back("MET_ht260mht60v2");
  //vhname_num.push_back("MET_ht250mht60v2and3");
  //vhname_num.push_back("MET_ht250mht70v1");
  //vhname_num.push_back("MET_ht300mht75v7and8");
  //vhname_num.push_back("MET_ht300mht80v1and2");
  //vhname_num.push_back("MET_ht300mht90v2");
  //vhname_num.push_back("MET_ht350mht90v1");
  //vhname_num.push_back("MET_ht350mht110v3");

  //Nominal sample: HT>400, 0L, mindphin>4
  if( theRegionType_ == kZeroLepHighMDP){
    vhname_denom.push_back("MET_ht260mht60v2_NoLepHDP_den");
    vhname_denom.push_back("MET_ht250mht60v2and3_NoLepHDP_den");
    vhname_denom.push_back("MET_ht250mht70v1_NoLepHDP_den");
    vhname_denom.push_back("MET_ht300mht75v7and8_NoLepHDP_den");
    vhname_denom.push_back("MET_ht300mht80v1and2_NoLepHDP_den");
    vhname_denom.push_back("MET_ht300mht90v2_NoLepHDP_den");
    vhname_denom.push_back("MET_ht350mht90v1_NoLepHDP_den");
    vhname_denom.push_back("MET_ht350mht110v3_NoLepHDP_den");
    vhname_num.push_back("MET_ht260mht60v2_NoLepHDP");
    vhname_num.push_back("MET_ht250mht60v2and3_NoLepHDP");
    vhname_num.push_back("MET_ht250mht70v1_NoLepHDP");  
    vhname_num.push_back("MET_ht300mht75v7and8_NoLepHDP");  
    vhname_num.push_back("MET_ht300mht80v1and2_NoLepHDP");  
    vhname_num.push_back("MET_ht300mht90v2_NoLepHDP");
    vhname_num.push_back("MET_ht350mht90v1_NoLepHDP");
    vhname_num.push_back("MET_ht350mht110v3_NoLepHDP");
  }
  //SL sample: HT>400, 1L, mindphin>4  -----------NOT USED!!!
  else if( theRegionType_ == kSingleLepHighMDP){
    vhname_denom.push_back("MET_ht260mht60v2_SLHDP_den");
    vhname_denom.push_back("MET_ht250mht60v2and3_SLHDP_den");
    vhname_denom.push_back("MET_ht250mht70v1_SLHDP_den");
    vhname_denom.push_back("MET_ht300mht75v7and8_SLHDP_den");
    vhname_denom.push_back("MET_ht300mht80v1and2_SLHDP_den");
    vhname_denom.push_back("MET_ht300mht90v2_SLHDP_den");
    vhname_denom.push_back("MET_ht350mht90v1_SLHDP_den");
    vhname_denom.push_back("MET_ht350mht110v3_SLHDP_den");
    vhname_num.push_back("MET_ht260mht60v2_SLHDP");
    vhname_num.push_back("MET_ht250mht60v2and3_SLHDP");
    vhname_num.push_back("MET_ht250mht70v1_SLHDP");  
    vhname_num.push_back("MET_ht300mht75v7and8_SLHDP");  
    vhname_num.push_back("MET_ht300mht80v1and2_SLHDP");
    vhname_num.push_back("MET_ht300mht90v2_SLHDP");
    vhname_num.push_back("MET_ht350mht90v1_SLHDP");
    vhname_num.push_back("MET_ht350mht110v3_SLHDP");
  }
  //LDP sample: HT>400, 0L, mindphin<4
  else if( theRegionType_ == kZeroLepLowMDP){
    vhname_denom.push_back("MET_ht260mht60v2_NoLepLDP_den");
    vhname_denom.push_back("MET_ht250mht60v2and3_NoLepLDP_den");
    vhname_denom.push_back("MET_ht250mht70v1_NoLepLDP_den");
    vhname_denom.push_back("MET_ht300mht75v7and8_NoLepLDP_den");
    vhname_denom.push_back("MET_ht300mht80v1and2_NoLepLDP_den");
    vhname_denom.push_back("MET_ht300mht90v2_NoLepLDP_den");
    vhname_denom.push_back("MET_ht350mht90v1_NoLepLDP_den");
    vhname_denom.push_back("MET_ht350mht110v3_NoLepLDP_den");
    vhname_num.push_back("MET_ht260mht60v2_NoLepLDP");
    vhname_num.push_back("MET_ht250mht60v2and3_NoLepLDP");
    vhname_num.push_back("MET_ht250mht70v1_NoLepLDP");  
    vhname_num.push_back("MET_ht300mht75v7and8_NoLepLDP");  
    vhname_num.push_back("MET_ht300mht80v1and2_NoLepLDP");  
    vhname_num.push_back("MET_ht300mht90v2_NoLepLDP");
    vhname_num.push_back("MET_ht350mht90v1_NoLepLDP");
    vhname_num.push_back("MET_ht350mht110v3_NoLepLDP");
  }
  //SL sample: HT>400, 1L, MuHad and EleHad datasets
  else if( theRegionType_ == kSingleMu){
    vhname_denom.push_back("histoHLT_HT300_MHT75_muTot");//Assume MHT75 efficiency for first 0.3 /fb (MHT60 & MHT70) as well
    vhname_denom.push_back("histoHLT_HT300_MHT75_muTot");//Assume MHT75 efficiency for first 0.3 /fb (MHT60 & MHT70) as well
    vhname_denom.push_back("histoHLT_HT300_MHT75_muTot");//Assume MHT75 efficiency for first 0.3 /fb (MHT60 & MHT70) as well
    vhname_denom.push_back("histoHLT_HT300_MHT75_muTot");
    vhname_denom.push_back("histoHLT_HT300_MHT80_muTot");
    vhname_denom.push_back("forbayesHLT_HT300_Mu15_PFMHT40_v1HT");
    vhname_denom.push_back("forbayesHLT_HT350_Mu5_PFMHT45_v8HT");
    vhname_denom.push_back("forbayesHLT_HT350_Mu5_PFMHT45_v12HT");
    vhname_num.push_back("histoHLT_HT300_MHT75_muPass");//Assume MHT75 efficiency for first 0.3 /fb (MHT60 & MHT70) as well
    vhname_num.push_back("histoHLT_HT300_MHT75_muPass");//Assume MHT75 efficiency for first 0.3 /fb (MHT60 & MHT70) as well
    vhname_num.push_back("histoHLT_HT300_MHT75_muPass");//Assume MHT75 efficiency for first 0.3 /fb (MHT60 & MHT70) as well
    vhname_num.push_back("histoHLT_HT300_MHT75_muPass");
    vhname_num.push_back("histoHLT_HT300_MHT80_muPass");
    vhname_num.push_back("forbayesHLT_HT300_MHT90_v2muHT");
    vhname_num.push_back("forbayesHLT_HT350_MHT90_v1muHT");
    vhname_num.push_back("forbayesHLT_HT350_MHT110_v3muHT");
  }
  else if( theRegionType_ == kSingleEle){
    vhname_denom.push_back("histoHLT_HT300_MHT75_elTot");//Assume MHT75 efficiency for first 0.3 /fb (MHT60 & MHT70) as well
    vhname_denom.push_back("histoHLT_HT300_MHT75_elTot");//Assume MHT75 efficiency for first 0.3 /fb (MHT60 & MHT70) as well
    vhname_denom.push_back("histoHLT_HT300_MHT75_elTot");//Assume MHT75 efficiency for first 0.3 /fb (MHT60 & MHT70) as well
    vhname_denom.push_back("histoHLT_HT300_MHT75_elTot");
    vhname_denom.push_back("histoHLT_HT300_MHT80_elTot"); 
    vhname_denom.push_back("forbayesHLT_HT300_Ele5_PFMHT40_v6HT");
    vhname_denom.push_back("forbayesHLT_HT350_Ele5_PFMHT45_v6HT");
    vhname_denom.push_back("forbayesHLT_HT350_Ele5_PFMHT45_v10HT");
    vhname_num.push_back("histoHLT_HT300_MHT75_elPass");//Assume MHT75 efficiency for first 0.3 /fb (MHT60 & MHT70) as well
    vhname_num.push_back("histoHLT_HT300_MHT75_elPass");//Assume MHT75 efficiency for first 0.3 /fb (MHT60 & MHT70) as well
    vhname_num.push_back("histoHLT_HT300_MHT75_elPass");//Assume MHT75 efficiency for first 0.3 /fb (MHT60 & MHT70) as well
    vhname_num.push_back("histoHLT_HT300_MHT75_elPass"); 
    vhname_num.push_back("histoHLT_HT300_MHT80_elPass");
    vhname_num.push_back("forbayesHLT_HT300_MHT90_v2elHT");
    vhname_num.push_back("forbayesHLT_HT350_MHT90_v1elHT");
    vhname_num.push_back("forbayesHLT_HT350_MHT110_v3elHT");
  }
  else{std::cout << "didn't specify kinematic region of interest" << std::endl; assert(0);}

  std::vector<TH1F*> vh_denom;
  std::vector<TH1F*> vh_num;

  std::vector<double> effSB, effSB_plus, effSB_minus;
  std::vector<double> effSIG, effSIG_plus, effSIG_minus;

  for(uint i = 0; i< vhname_denom.size(); ++i){
    char fhname_num[200],fhname_denom[200];
    sprintf(fhname_num,"%s",vhname_num.at(i).c_str());
    sprintf(fhname_denom,"%s",vhname_denom.at(i).c_str());
  
    if( (theRegionType_==kSingleMu || theRegionType_==kSingleEle) && i<5){//kludge for SL case, first two triggers are in different file
      vh_denom.push_back( (TH1F *)f_trig2->Get(fhname_denom) );
      vh_num.push_back( (TH1F *)f_trig2->Get(fhname_num) );
    }
    else{
      vh_denom.push_back( (TH1F *)f_trig1->Get(fhname_denom) );
      vh_num.push_back( (TH1F *)f_trig1->Get(fhname_num) );
    }

    //vh_denom.at(i)->Sumw2();
    //vh_num.at(i)->Sumw2();

    //std::cout << "nbins = " <<   vh_num.at(i)->GetNbinsX() << std::endl;

    //COMPUTE SB EFFICIENCY
    float SB_den = vh_denom.at(i)->Integral(vh_denom.at(i)->FindBin(sb_low), vh_denom.at(i)->FindBin(sb_high)-1);
    float SB_num = vh_num.at(i)->Integral(vh_num.at(i)->FindBin(sb_low), vh_num.at(i)->FindBin(sb_high)-1);

    //std::cout << "Integral SB den= " << SB_den << std::endl;
    //std::cout << "Integral SB num= " << SB_num << std::endl;

    TH1F * h_num_temp= new TH1F("h_num_temp","temp",1,0,1);
    TH1F * h_den_temp= new TH1F("h_den_temp","temp",1,0,1);
    double x,y;
      
    if(SB_den){
      vlumi_ignore_SB.push_back(0);
      totallumi_SB += vlumi.at(i);
      h_num_temp->Sumw2(); h_den_temp->Sumw2();
      h_num_temp->SetBinContent(1,SB_num);
      h_den_temp->SetBinContent(1,SB_den);
      TGraphAsymmErrors * h_effSB = new TGraphAsymmErrors(h_num_temp,h_den_temp,"cl=0.683 b(1,1) mode");
      //h_eff->Print();
      h_effSB->GetPoint(0,x,y);
      //std::cout << "x = " << x << ",y = " << y << std::endl;
      std::cout << vhname_num.at(i)<< ": Efficiency (SB) = " << y 
		<< " + " << h_effSB->GetErrorYhigh(0) 
		<< " - " << h_effSB->GetErrorYlow(0) << std::endl;

      //if(i<3){
      //	std::cout << "TEMPORARRRYYYYYYYYYYYYYYYYYYY   ASSUME 100%!!!!" << std::endl;
      //	effSB.push_back(1); effSB_plus.push_back(0); effSB_minus.push_back(0);
      //}
      //else 
      effSB.push_back(y); effSB_plus.push_back(h_effSB->GetErrorYhigh(0)); effSB_minus.push_back(h_effSB->GetErrorYlow(0));


      delete h_effSB;
    }
    else{
      vlumi_ignore_SB.push_back(1);
      totallumi_ignore_SB += vlumi.at(i);
      effSB.push_back(-99); effSB_plus.push_back(-99); effSB_minus.push_back(-99);//will be ignored anyway
      std::cout << vhname_num.at(i)<< ": Efficiency (SB) = " << "IGNORED - no denominator stats" << std::endl;
    }
    

    //COMPUTE SIG EFFICIENCY
    float SIG_den = vh_denom.at(i)->Integral(vh_denom.at(i)->FindBin(250), vh_denom.at(i)->GetNbinsX()+1);
    float SIG_num = vh_num.at(i)->Integral(vh_num.at(i)->FindBin(250), vh_num.at(i)->GetNbinsX()+1);

    if(is1BT){
      SIG_den = vh_denom.at(i)->Integral(vh_denom.at(i)->FindBin(500), vh_denom.at(i)->GetNbinsX()+1);
      SIG_num = vh_num.at(i)->Integral(vh_num.at(i)->FindBin(500), vh_num.at(i)->GetNbinsX()+1);
    }
    if(is2BT){
      SIG_den = vh_denom.at(i)->Integral(vh_denom.at(i)->FindBin(300), vh_denom.at(i)->GetNbinsX()+1);
      SIG_num = vh_num.at(i)->Integral(vh_num.at(i)->FindBin(300), vh_num.at(i)->GetNbinsX()+1);
    }

    if(isSBANDSIG){
      SIG_den = vh_denom.at(i)->Integral(vh_denom.at(i)->FindBin(200), vh_denom.at(i)->GetNbinsX()+1);
      SIG_num = vh_num.at(i)->Integral(vh_num.at(i)->FindBin(200), vh_num.at(i)->GetNbinsX()+1);
    }

    if(SIG_den){
      vlumi_ignore_SIG.push_back(0);
      totallumi_SIG += vlumi.at(i);
      h_num_temp->SetBinContent(1,SIG_num);
      h_den_temp->SetBinContent(1,SIG_den);
      TGraphAsymmErrors * h_effSIG = new TGraphAsymmErrors(h_num_temp,h_den_temp,"cl=0.683 b(1,1) mode");
      h_effSIG->GetPoint(0,x,y);
      std::cout << vhname_num.at(i)<< ": Efficiency (SIG) = "<< y 
		<< " + " << h_effSIG->GetErrorYhigh(0) 
		<< " - " << h_effSIG->GetErrorYlow(0) << std::endl;
      effSIG.push_back(y); effSIG_plus.push_back(h_effSIG->GetErrorYhigh(0)); effSIG_minus.push_back(h_effSIG->GetErrorYlow(0));
      
      delete h_effSIG;
    }
    else{
      vlumi_ignore_SIG.push_back(1);
      totallumi_ignore_SIG += vlumi.at(i);
      effSIG.push_back(-99); effSIG_plus.push_back(-99); effSIG_minus.push_back(-99);//will be ignored anyway
      std::cout << vhname_num.at(i)<< ": Efficiency (SIG) = " << "IGNORED - no denominator stats" << std::endl;
    }

    delete h_num_temp;
    delete h_den_temp;


  }


  //compute overall efficiency
  double effSB_overall = 0, effSB_overall_errp = 0, effSB_overall_errm = 0;
  double effSIG_overall = 0, effSIG_overall_errp = 0, effSIG_overall_errm = 0;
  for(uint i =0; i<effSB.size(); ++i){

    if(!vlumi_ignore_SB.at(i)){
      effSB_overall += effSB.at(i)*vlumi.at(i)/totallumi_SB;
      effSB_overall_errp += effSB_plus.at(i)*effSB_plus.at(i)*( vlumi.at(i)*vlumi.at(i) / totallumi_SB / totallumi_SB );
      effSB_overall_errm += effSB_minus.at(i)*effSB_minus.at(i)*( vlumi.at(i)*vlumi.at(i) / totallumi_SB / totallumi_SB );
    }
    
    if(!vlumi_ignore_SIG.at(i)){
      effSIG_overall += effSIG.at(i)*vlumi.at(i)/totallumi_SIG;
      effSIG_overall_errp += effSIG_plus.at(i)*effSIG_plus.at(i)*( vlumi.at(i)*vlumi.at(i) / totallumi_SIG / totallumi_SIG );
      effSIG_overall_errm += effSIG_minus.at(i)*effSIG_minus.at(i)*( vlumi.at(i)*vlumi.at(i) / totallumi_SIG / totallumi_SIG );
    }
  }
  effSB_overall_errp = sqrt(effSB_overall_errp);
  effSB_overall_errm = sqrt(effSB_overall_errm);

  effSIG_overall_errp = sqrt(effSIG_overall_errp);
  effSIG_overall_errm = sqrt(effSIG_overall_errm);

  std::cout << "*********SUMMARY**********" << std::endl;

  std::cout << "effSB_overall = " << effSB_overall 
	    << " + " << effSB_overall_errp << " - " << effSB_overall_errm 
	    << ", fraction of lumi ignored = " << totallumi_ignore_SB/totallumi << std::endl;
  std::cout << "effSIG_overall = " << effSIG_overall 
	    << " + " << effSIG_overall_errp << " - " << effSIG_overall_errm 
	    << ", fraction of lumi ignored = " << totallumi_ignore_SIG/totallumi << std::endl;

  std::cout << "**************************" << std::endl;
}


