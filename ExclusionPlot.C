#include <sstream>
#include <string>
#include <iomanip>
#include "TSpline.h"

#include "ExclusionPlot.hh"

#include "TROOT.h" 
#include "TList.h"

#include "TLegend.h"
#include "TLegendEntry.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TMarker.h"
#include <vector>
#include "TMath.h"
#include "TString.h"

//Produce the limit plot with the command: root -l ExclusionPlot.C+

//some RA2b-specific kludging
//relevant modes for our PAS are 'all' and 'ge1btight'
const TString RA2bmode = "all"; //"all", "allPlusRA1", ge1bloose etc

const int fillstyle = 3002;

void ExclusionPlot(){
  gStyle->SetPalette(1);

  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);

  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");



  Int_t tanBeta = 40;
  Bool_t plotLO = false;
   
  CommandMSUGRA("35pb_expected_11.root",tanBeta, plotLO);
}


void CommandMSUGRA(TString plotName_,Int_t tanBeta_, Bool_t plotLO_){
  gROOT->SetStyle("CMS");//jmt specific
  gROOT->ForceStyle();

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1); 
  gStyle->SetTextFont(42);
  gStyle->SetFrameBorderMode(0);

  //convert tanb value to string
  std::stringstream tmp;
  tmp << tanBeta_;
  TString tanb( tmp.str() );
  
  
  // Output file
  std::cout << " create " << plotName_ << std::endl;
  TFile* output = new TFile( plotName_, "RECREATE" );
  if ( !output || output->IsZombie() ) { std::cout << " zombie alarm output is a zombie " << std::endl; }
  

  //set old exclusion Limits
  TGraph* LEP_ch = set_lep_ch(tanBeta_);
  TGraph* LEP_sl = set_lep_sl(tanBeta_);//slepton curve
  TGraph* TEV_sg_cdf = set_tev_sg_cdf(tanBeta_);//squark gluino cdf
  TGraph* TEV_sg_d0 = set_tev_sg_d0(tanBeta_);//squark gluino d0
  //  TGraph* TEV_tlp_cdf = set_tev_tlp_cdf(tanBeta_);//trilepton cdf
  //  TGraph* TEV_tlp_d0 = set_tev_tlp_d0(tanBeta_);//trilepton d0
  TGraph* stau   = set_tev_stau(tanBeta_);//stau 
  TGraph* NoEWSB = set_NoEWSB(tanBeta_); 

  TGraph* TEV_sn_d0_1 = set_sneutrino_d0_1(tanBeta_);
  TGraph* TEV_sn_d0_2 = set_sneutrino_d0_2(tanBeta_);

  int nPoints = nSusyGridPoints();
  double m0[nPoints],m12[nPoints],squarkMass[nPoints],gluinoMass[nPoints];

  susyGrid(m0,m12,squarkMass,gluinoMass);

  TGraph2D* squarkMasses = new TGraph2D("squarkMasses","",nPoints,m0,m12,squarkMass);
  TGraph2D* gluinoMasses = new TGraph2D("gluinoMasses","",nPoints,m0,m12,gluinoMass);

  TH2D* gluinoMassPlot = gluinoMasses->GetHistogram();
  TH2D* squarkMassPlot = squarkMasses->GetHistogram();

  //constant ssqquark and gluino lines
  TF1* lnsq[15];
  TF1* lngl[15];

  TGraph* lnsq_40[15];
  TGraph* lngl_40[15];
  
  TLatex* sq_text[15];
  TLatex* gl_text[15];

  TLatex* sq_40_text[15];
  TLatex* gl_40_text[15];

  for(int i = 1; i < 15; i++){
    //lnsq[i] = constant_squark(tanBeta_,i);
    //sq_text[i] = constant_squark_text(i,*lnsq[i],tanBeta_);
    //lngl[i] = constant_gluino(tanBeta_,i);
    //gl_text[i] = constant_gluino_text(i,*lngl[i]);
    lnsq_40[i] = constant_mass(i*250,squarkMasses);
    lngl_40[i] = constant_mass(i*250,gluinoMasses);
    sq_40_text[i] = constant_squark_text_tanBeta40(i*250,lnsq_40[i]);
    gl_40_text[i] = constant_gluino_text_tanBeta40(i*250,lngl_40[i]);;
  }


  //Legends
  TLegend* legst  = makeStauLegend(0.05,tanBeta_);
  TLegend* legNoEWSB  = makeNoEWSBLegend(0.05,tanBeta_);
  TLegend* legexp = makeExpLegend( *TEV_sg_cdf,*TEV_sg_d0,*LEP_ch,*LEP_sl,*TEV_sn_d0_1,0.035,tanBeta_);
  
 
  //make Canvas
  TCanvas* cvsSys = new TCanvas("cvsnm","cvsnm",0,0,800,600);
  gStyle->SetOptTitle(0);
  cvsSys->SetFillColor(0);
  cvsSys->GetPad(0)->SetRightMargin(0.07);
  cvsSys->Range(-120.5298,26.16437,736.0927,750);
  //  cvsSys->Range(-50.5298,26.16437,736.0927,500);
  cvsSys->SetFillColor(0);
  cvsSys->SetBorderMode(0);
  cvsSys->GetPad(0)->SetBorderMode(0);
  cvsSys->GetPad(0)->SetBorderSize(2);
  cvsSys->GetPad(0)->SetLeftMargin(0.1407035);
  cvsSys->GetPad(0)->SetTopMargin(0.08);
  cvsSys->GetPad(0)->SetBottomMargin(0.13);

  cvsSys->SetTitle("tan#beta="+tanb);
 
  output->cd();
  
//and now  the exclusion limits


  TGraph* SSdilep;
  TGraphErrors* OSdilep;
  TGraphErrors* RA1;

  TGraphErrors* RA1_old;
  TGraphErrors* RA5_old;
  TGraphErrors* RA6_old;

  TGraph* RA2b_1b_loose;
  TGraph* RA2b_1b_tight;
  TGraph* RA2b_2b_loose;
  TGraph* RA2b_2b_tight;

  TGraph* RA2b_1b_loose_exp;
  TGraph* RA2b_1b_tight_exp;
  TGraph* RA2b_2b_loose_exp;
  TGraph* RA2b_2b_tight_exp;

  TGraph* RA2b_1b_loose_exp_p;
  TGraph* RA2b_1b_tight_exp_p;
  TGraph* RA2b_2b_loose_exp_p;
  TGraph* RA2b_2b_tight_exp_p;

  TGraph* RA2b_1b_loose_exp_m;
  TGraph* RA2b_1b_tight_exp_m;
  TGraph* RA2b_2b_loose_exp_m;
  TGraph* RA2b_2b_tight_exp_m;

  TGraph* RA2b_1b_loose_shade;
  TGraph* RA2b_1b_tight_shade;
  TGraph* RA2b_2b_loose_shade;
  TGraph* RA2b_2b_tight_shade;

  TSpline3* RA1_tb40 =getCLs1080ObsNLOtb40();

  if (tanBeta_ == 10) {
    SSdilep = SSdilep_NLO();
    OSdilep = OSdilep_NLO();
    RA1 = RA1_NLO();

    RA1_old = getRA1Observed_NLO_tanBeta10();
    RA5_old = getRA5Observed_NLO_tanBeta10();
    RA6_old = getRA6Observed_NLO_tanBeta10();

  }
  if(tanBeta_ == 40)
    {
//       RA2b_1b_loose = RA2b_limit("an-scanplot-unblind-tb40-withcontam-ge1b-loose.root", "hsusyscanExcluded");
//       RA2b_2b_loose = RA2b_limit("an-scanplot-unblind-tb40-withcontam-ge2b-loose.root", "hsusyscanExcluded");
//       RA2b_1b_tight = RA2b_limit("an-scanplot-unblind-tb40-withcontam-ge1b-tight.root", "hsusyscanExcluded");
//       RA2b_2b_tight = RA2b_limit("an-scanplot-unblind-tb40-withcontam-ge2b-tight.root", "hsusyscanExcluded");
//       RA2b_1b_loose = RA2b_limit("/afs/cern.ch/user/o/owen/public/RA2b/clsplots-tb40-ge1bloose.root", "hcls");
//       RA2b_2b_loose = RA2b_limit("/afs/cern.ch/user/o/owen/public/RA2b/clsplots-tb40-ge2bloose.root", "hcls");
//       RA2b_1b_tight = RA2b_limit("/afs/cern.ch/user/o/owen/public/RA2b/clsplots-tb40-ge1btight.root", "hcls");
//       RA2b_2b_tight = RA2b_limit("/afs/cern.ch/user/o/owen/public/RA2b/clsplots-tb40-ge2btight.root", "hcls");

/* 
     TString ra2bfile= "RA2b_tb40_exclusion.25Sep.root";
      RA2b_1b_loose = RA2b_limit(ra2bfile,"curve4_ge1bloose");
      RA2b_1b_tight = RA2b_limit(ra2bfile,"curve4_ge1btight");
      RA2b_2b_loose = RA2b_limit(ra2bfile,"curve4_ge2bloose");
      RA2b_2b_tight = RA2b_limit(ra2bfile,"curve4_ge2btight");

      RA2b_1b_loose_exp = RA2b_limit(ra2bfile,"curve4_1bloose_exp");
      RA2b_1b_tight_exp = RA2b_limit(ra2bfile,"curve4_1btight_exp");
      RA2b_2b_loose_exp = RA2b_limit(ra2bfile,"curve4_2bloose_exp");
      RA2b_2b_tight_exp = RA2b_limit(ra2bfile,"curve4_2btight_exp");

      RA2b_1b_loose_exp_p = RA2b_limit(ra2bfile,"curve4_1bloose_exp_plus");
      RA2b_1b_tight_exp_p = RA2b_limit(ra2bfile,"curve4_1btight_exp_plus");
      RA2b_2b_loose_exp_p = RA2b_limit(ra2bfile,"curve4_2bloose_exp_plus");
      RA2b_2b_tight_exp_p = RA2b_limit(ra2bfile,"curve4_2btight_exp_plus");

      RA2b_1b_loose_exp_m = RA2b_limit(ra2bfile,"curve4_1bloose_exp_minus");
      RA2b_1b_tight_exp_m = RA2b_limit(ra2bfile,"curve4_1btight_exp_minus");
      RA2b_2b_loose_exp_m = RA2b_limit(ra2bfile,"curve4_2bloose_exp_minus");
      RA2b_2b_tight_exp_m = RA2b_limit(ra2bfile,"curve4_2btight_exp_minus");
*/
      RA2b_1b_loose = get_RA2b_1bloose();
      RA2b_1b_tight = get_RA2b_1btight();
      RA2b_2b_loose = get_RA2b_2bloose();
      RA2b_2b_tight = get_RA2b_2btight();

      RA2b_1b_loose_exp = get_RA2b_1bloose_exp();
      RA2b_1b_tight_exp = get_RA2b_1btight_exp();
      RA2b_2b_loose_exp = get_RA2b_2bloose_exp();
      RA2b_2b_tight_exp = get_RA2b_2btight_exp();
    
      RA2b_1b_loose_exp_p = get_RA2b_1bloose_exp_p();
      RA2b_1b_tight_exp_p = get_RA2b_1btight_exp_p();
      RA2b_2b_loose_exp_p = get_RA2b_2bloose_exp_p();
      RA2b_2b_tight_exp_p = get_RA2b_2btight_exp_p();

      RA2b_1b_loose_exp_m = get_RA2b_1bloose_exp_m();
      RA2b_1b_tight_exp_m = get_RA2b_1btight_exp_m();
      RA2b_2b_loose_exp_m = get_RA2b_2bloose_exp_m();
      RA2b_2b_tight_exp_m = get_RA2b_2btight_exp_m();
   
      cout<<"Getting the shaded regions"<<endl;
      RA2b_1b_loose_shade = getShadedRegion(RA2b_1b_loose_exp_p,RA2b_1b_loose_exp_m);
      RA2b_1b_tight_shade = getShadedRegion(RA2b_1b_tight_exp_p,RA2b_1b_tight_exp_m);
      RA2b_2b_loose_shade = getShadedRegion(RA2b_2b_loose_exp_p,RA2b_2b_loose_exp_m);
      RA2b_2b_tight_shade = getShadedRegion(RA2b_2b_tight_exp_p,RA2b_2b_tight_exp_m);
      cout<<"DONE Getting the shaded regions"<<endl;
    }



  double m0min = 0;
  if (tanBeta_ == 40) m0min=400; 
  TH2D* hist = new TH2D("h","h",100,m0min,2000,100,120,700);
  hist->Draw();  
  hist->GetXaxis()->SetTitle("m_{0} [GeV]");
  hist->GetYaxis()->SetTitle("m_{1/2} [GeV]");
  hist->GetXaxis()->SetTitleOffset(.9);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleOffset(1.0);
  hist->GetYaxis()->SetTitleSize(0.06);

  hist->GetXaxis()->SetNdivisions(506);
  //  if (tanBeta_ == 50)  hist->GetXaxis()->SetNdivisions(504);
  hist->GetYaxis()->SetNdivisions(506);

  int col[]={2,3,4};

  //SSdilep->SetLineColor(kGreen+2);
  //SSdilep->SetLineStyle(1);
  //SSdilep->SetLineWidth(3);
  //
  //OSdilep->SetLineColor(kCyan+2);
  //OSdilep->SetLineStyle(1);
  //OSdilep->SetLineWidth(3);
  //
  //RA1->SetLineColor(kRed+2);
  //RA1->SetLineStyle(1);
  //RA1->SetLineWidth(3);
  RA1_tb40->SetLineColor(kBlack);
  RA1_tb40->SetLineStyle(7);
  RA1_tb40->SetLineWidth(3);
  RA1_tb40->SetName("RA1_tb40");
  //
  //TSpline3 *sRA1 = new TSpline3("sRA1",RA1_old);
  //sRA1->SetLineColor(kRed+2);
  ////sRA1->SetLineStyle(5);
  //sRA1->SetLineStyle(2);
  //sRA1->SetLineWidth(3);
  //
  //RA5_old->SetLineColor(kGreen+2);
  ////RA5_old->SetLineStyle(5);
  //RA5_old->SetLineStyle(2);
  //RA5_old->SetLineWidth(3);
  //
  //RA6_old->SetLineColor(kCyan+2);
  ////RA6_old->SetLineStyle(1);
  //RA6_old->SetLineStyle(2);
  //RA6_old->SetLineWidth(3);

  if (RA2bmode.Contains("all")) {
  RA2b_1b_loose->SetLineColor(kRed+2);
  RA2b_1b_loose->SetLineStyle(2);
  RA2b_1b_loose->SetLineWidth(3);
         
  RA2b_1b_tight->SetLineColor(kRed+2);
  RA2b_1b_tight->SetLineStyle(1);
  RA2b_1b_tight->SetLineWidth(3);
         
  RA2b_2b_loose->SetLineColor(kGreen+2);
  RA2b_2b_loose->SetLineStyle(2);
  RA2b_2b_loose->SetLineWidth(3);
         
  RA2b_2b_tight->SetLineColor(kGreen+2);
  RA2b_2b_tight->SetLineStyle(1);
  RA2b_2b_tight->SetLineWidth(3);
  }
  else {
  RA2b_1b_loose->SetLineColor(kRed);
  RA2b_1b_loose->SetLineStyle(1);
  RA2b_1b_loose->SetLineWidth(3);
         
  RA2b_1b_tight->SetLineColor(kRed);
  RA2b_1b_tight->SetLineStyle(1);
  RA2b_1b_tight->SetLineWidth(3);
         
  RA2b_2b_loose->SetLineColor(kRed);
  RA2b_2b_loose->SetLineStyle(1);
  RA2b_2b_loose->SetLineWidth(3);
         
  RA2b_2b_tight->SetLineColor(kRed);
  RA2b_2b_tight->SetLineStyle(1);
  RA2b_2b_tight->SetLineWidth(3);

  RA2b_1b_loose_exp->SetLineColor(kBlue);
  RA2b_2b_loose_exp->SetLineColor(kBlue);
  RA2b_1b_tight_exp->SetLineColor(kBlue);
  RA2b_2b_tight_exp->SetLineColor(kBlue);

  RA2b_1b_loose_exp->SetLineStyle(5);
  RA2b_2b_loose_exp->SetLineStyle(5);
  RA2b_1b_tight_exp->SetLineStyle(5);
  RA2b_2b_tight_exp->SetLineStyle(5);

  RA2b_1b_loose_exp->SetLineWidth(3);
  RA2b_2b_loose_exp->SetLineWidth(3);
  RA2b_1b_tight_exp->SetLineWidth(3);
  RA2b_2b_tight_exp->SetLineWidth(3);

  int acolor=kCyan+2;
  RA2b_1b_loose_exp_p->SetLineColor(acolor);
  RA2b_2b_loose_exp_p->SetLineColor(acolor);
  RA2b_1b_tight_exp_p->SetLineColor(acolor);
  RA2b_2b_tight_exp_p->SetLineColor(acolor);

  RA2b_1b_loose_exp_p->SetLineStyle(1);
  RA2b_2b_loose_exp_p->SetLineStyle(1);
  RA2b_1b_tight_exp_p->SetLineStyle(1);
  RA2b_2b_tight_exp_p->SetLineStyle(1);

  RA2b_1b_loose_exp_p->SetLineWidth(3);
  RA2b_2b_loose_exp_p->SetLineWidth(3);
  RA2b_1b_tight_exp_p->SetLineWidth(3);
  RA2b_2b_tight_exp_p->SetLineWidth(3);

  RA2b_1b_loose_exp_m->SetLineColor(acolor);
  RA2b_2b_loose_exp_m->SetLineColor(acolor);
  RA2b_1b_tight_exp_m->SetLineColor(acolor);
  RA2b_2b_tight_exp_m->SetLineColor(acolor);

  RA2b_1b_loose_exp_m->SetLineStyle(1);
  RA2b_2b_loose_exp_m->SetLineStyle(1);
  RA2b_1b_tight_exp_m->SetLineStyle(1);
  RA2b_2b_tight_exp_m->SetLineStyle(1);

  RA2b_1b_loose_exp_m->SetLineWidth(3);
  RA2b_2b_loose_exp_m->SetLineWidth(3);
  RA2b_1b_tight_exp_m->SetLineWidth(3);
  RA2b_2b_tight_exp_m->SetLineWidth(3);

  RA2b_1b_tight_shade->SetFillStyle(fillstyle);
  RA2b_1b_tight_shade->SetFillColor(acolor);

  RA2b_1b_loose_shade->SetFillStyle(fillstyle);
  RA2b_1b_loose_shade->SetFillColor(acolor);

  RA2b_2b_tight_shade->SetFillStyle(fillstyle);
  RA2b_2b_tight_shade->SetFillColor(acolor);

  RA2b_2b_loose_shade->SetFillStyle(fillstyle);
  RA2b_2b_loose_shade->SetFillColor(acolor);

  }
  
  TLegend* myleg;

  float leg_x1=0.39+0.23;
  float leg_y1=0.65+0.05;
  float leg_x2= 0.55+0.25;
  float leg_y2= 0.84+0.05;

  if (RA2bmode.Contains("all")) {
    leg_y1 -= 0.1;
  }

  if( plotLO_ ) myleg = new TLegend(0.3,0.65,0.65,0.8,NULL,"brNDC");
  else          myleg = new TLegend(leg_x1,leg_y1,leg_x2,leg_y2,NULL,"brNDC");


  myleg->SetFillColor(0); 
  myleg->SetShadowColor(0);
  myleg->SetTextSize(0.04);
  myleg->SetBorderSize(0);

  TLegendEntry *entry=0;
//   entry= myleg->AddEntry("ge1bLoose","LEP2 #tilde{#chi}_{1}^{#pm}","f");

//   entry->SetFillColor(3);
//   entry->SetLineColor(3);
//   entry->SetFillStyle(1001);


  if (RA2bmode.Contains("all")) {
    entry= myleg->AddEntry("ge1bLoose","#geq 1b Loose","l");
    entry->SetLineColor(1);
    entry->SetLineStyle(2);
    entry->SetLineWidth(3);
    entry->SetLineColor(kRed+2);
    entry->SetTextColor(kRed+2);
    
    entry=myleg->AddEntry("ge2bLoose","#geq 2b Loose","l");
    entry->SetLineColor(1);
    entry->SetLineStyle(2);
    entry->SetLineWidth(3);
    entry->SetLineColor(kGreen+2);
    entry->SetTextColor(kGreen+2);
    
    entry=myleg->AddEntry("ge1bTight","#geq 1b Tight","l");
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(3);
    entry->SetLineColor(kRed+2);
    entry->SetTextColor(kRed+2);
    
    entry=myleg->AddEntry("ge2bTight","#geq 2b Tight","l");
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(3);
    entry->SetLineColor(kGreen+2);
    entry->SetTextColor(kGreen+2);
  }
  else if (RA2bmode.Contains("ge")) {
    entry=myleg->AddEntry("ge1bTight","Observed Limit","l");
    entry->SetLineStyle(1);
    entry->SetLineWidth(3);
    entry->SetLineColor(kRed);
    entry->SetTextColor(kBlack);

    entry=myleg->AddEntry("ge1bTight_exp","Expected Limit #pm 1#sigma","lf");
    entry->SetFillStyle (fillstyle);
    entry->SetFillColor (kCyan+2);
    entry->SetLineStyle(5);
    entry->SetLineWidth(3);
    entry->SetLineColor(kBlue);
    entry->SetTextColor(kBlack);


  }
  
   if (RA2bmode=="allPlusRA1") {
    entry=myleg->AddEntry("RA1_tb40","CMS #alpha_{T}","l");
    entry->SetLineColor(kBlack);
    entry->SetLineStyle(7);
    entry->SetLineWidth(3);
    entry->SetTextColor(kBlack);
   }  

  //constant squark and gluino mass contours
  for (int it=2;it<9;it++) {  
    if(it<7){
      if(lngl_40[it]!=0)lngl_40[it]->Draw("samec");   
      if(gl_40_text[it]!=0)gl_40_text[it]->Draw();
    }
      if(lnsq_40[it]!=0)lnsq_40[it]->Draw("samec");
    if(it<6){
      if(sq_40_text[it]!=0)sq_40_text[it]->Draw();
    }
  }
  //SSdilep->Draw("samec");
  //OSdilep->Draw("samec");
  //RA1->Draw("samec");
  //
  //sRA1->Draw("same"); 
  //RA5_old->Draw("c same");
  //RA6_old->Draw("c same");
  TString drawopt="samel"; //default choice
  if (RA2bmode.Contains("all"))  RA2b_1b_loose->Draw(drawopt);
  if (RA2bmode.Contains("all")) RA2b_2b_loose->Draw(drawopt);
  if (RA2bmode.Contains("all")) RA2b_1b_tight->Draw(drawopt);
  if (RA2bmode.Contains("all")) RA2b_2b_tight->Draw(drawopt);
  if (RA2bmode=="allPlusRA1") RA1_tb40->Draw(drawopt);

  if (RA2bmode=="ge1btight") {
    RA2b_1b_tight_shade->Draw("f");
    RA2b_1b_tight_exp_p->Draw(drawopt);
    RA2b_1b_tight_exp_m->Draw(drawopt);
    RA2b_1b_tight_exp->Draw(drawopt);
    RA2b_1b_tight->Draw(drawopt);
  }
  else if (RA2bmode=="ge1bloose") {
    RA2b_1b_loose_shade->Draw("f");
    RA2b_1b_loose_exp_p->Draw(drawopt);
    RA2b_1b_loose_exp_m->Draw(drawopt);
    RA2b_1b_loose_exp->Draw(drawopt);
    RA2b_1b_loose->Draw(drawopt);
  }
  else if (RA2bmode=="ge2bloose") {
    RA2b_2b_loose_shade->Draw("f");
    RA2b_2b_loose_exp_p->Draw(drawopt);
    RA2b_2b_loose_exp_m->Draw(drawopt);
    RA2b_2b_loose_exp->Draw(drawopt);
    RA2b_2b_loose->Draw(drawopt);
  }
  else if (RA2bmode=="ge2btight") {
    RA2b_2b_tight_shade->Draw("f");
    RA2b_2b_tight_exp_p->Draw(drawopt);
    RA2b_2b_tight_exp_m->Draw(drawopt);
    RA2b_2b_tight_exp->Draw(drawopt);
    RA2b_2b_tight->Draw(drawopt);
  }
  //
  //
  //TLatex* RA1label = new TLatex(670,430.,"#alpha_{T}");
  ////TLatex* RA1label = new TLatex(80,288.,"#alpha_{T}");
  //RA1label->SetTextFont(42);
  //RA1label->SetTextSize(0.05);
  //RA1label->SetTextColor(kRed+2);
  //RA1label->Draw("same");

  TLatex* RA2blabel_2b=0;
  TLatex* RA2blabel_1b=0;
  if (false) {
    RA2blabel_2b = new TLatex(1150,330.,"#geq 2 b-tags");
    RA2blabel_2b->SetTextFont(42);
    RA2blabel_2b->SetTextSize(0.05);
    RA2blabel_2b->SetTextColor(kGreen+2);
    RA2blabel_2b->Draw("same");
    
    RA2blabel_1b = new TLatex(1150,430.,"#geq 1 b-tags");
    RA2blabel_1b->SetTextFont(42);
    RA2blabel_1b->SetTextSize(0.05);
    RA2blabel_1b->SetTextColor(kRed+2);
    RA2blabel_1b->Draw("same");
  }
  //
  //TLatex* RA5label = new TLatex(400,370.,"SS Dilepton");
  //RA5label->SetTextFont(42);
  ////RA5label->SetTextAngle(20);
  //RA5label->SetTextSize(0.04);
  //RA5label->SetTextColor(kGreen+2);
  //RA5label->Draw("same");
  //
  //TLatex* RA6label = new TLatex(650,215.,"OS Dilepton");
  //RA6label->SetTextFont(42);
  ////RA6label->SetTextAngle(8);
  //RA6label->SetTextSize(0.04);
  //RA6label->SetTextColor(kCyan+2);
  //RA6label->Draw("same");
  
  
  //exclusion limits previous experiments
  if(tanBeta_ == 3){
    TEV_sn_d0_1->Draw("fsame");
    TEV_sn_d0_2->Draw("fsame");
  }
  //  LEP_ch->Draw("fsame");
  if (tanBeta_ != 40) LEP_sl->Draw("fsame");

  //remove CDF/D0 excluded regions
  //  TEV_sg_cdf->Draw("fsame");
  //  TEV_sg_d0->Draw("same");  
  //  TEV_sg_d0->Draw("fsame");


  //other labels
  Double_t xpos = 0;
  Double_t xposi = 0;
  Double_t ypos = 0;
  if(tanBeta_ == 40) xposi = 180+160;
  if(tanBeta_ == 40) xpos = 400;//240;
  if(tanBeta_ == 40) ypos = -10;
  
  //TLatex* lumilabel = new TLatex(750 +xposi + 100,767.-154,"#sqrt{s} = 7 TeV, #scale[0.65]{#int}Ldt = 0.98 fb^{-1}");
  TLatex* lumilabel = new TLatex(925+xpos-50,767.-154+105,"#sqrt{s} = 7 TeV, L_{int} = 1.1 fb^{-1}");
  TLatex* integral_symbol = new TLatex(1287 +xposi + 100-85,767.-145+95,"#int");

  lumilabel->SetTextSize(0.05);
  integral_symbol->SetTextSize(0.03);
  lumilabel->Draw("same");
  //  integral_symbol->Draw("same");

  TLatex* cmslabel = new TLatex(10.+xpos,767.-154+105,"CMS Preliminary");
  cmslabel->SetTextSize(0.05);
  cmslabel->Draw("same");

  TString text_tanBeta;
  text_tanBeta =  "tan#beta = "+tanb+",  A_{0} = -500 GeV,  #mu > 0";
  TLatex* cmssmpars = new TLatex(/*530.+xpos,690.+ypos-130*/150+xpos,660,text_tanBeta);

  cmssmpars->SetTextSize(0.04);
  cmssmpars->Draw("same");

  TLatex* lep_chargino = new TLatex(250,135,"LEP2 #tilde{#chi}_{1}^{#pm}");
  lep_chargino->SetTextSize(0.03);
  lep_chargino->SetTextFont(42);
  //    lep_chargino->Draw("same");

  TLatex* lep_slepton = new TLatex(26,190,"LEP2 #tilde{#font[12]{l}}^{#pm}");
  lep_slepton->SetTextSize(0.03);
  lep_slepton->SetTextAngle(-83);
  lep_slepton->SetTextFont(42);
  //  lep_slepton->Draw("same");



  //LM points
  TMarker* LM0 = new TMarker(200.,160.,20);
  TMarker* LM1 = new TMarker(60.,250.,20);
  TMarker* LM3 = new TMarker(330.,240.,20);
  TMarker* LM6 = new TMarker(80.,400.,20);
    
  LM0->SetMarkerSize(1.2);
  LM1->SetMarkerSize(1.2);
    
  TLatex* tLM0 = new TLatex(205.,160.," LM0");
  tLM0->SetTextSize(0.035);
    
  TLatex* tLM1 = new TLatex(80.,245.,"LM1");
  tLM1->SetTextSize(0.035);
  
  //TLatex* tLM3 = new TLatex(350.,235.,"LM3 (tan#beta=20)");
  TLatex* tLM3 = new TLatex(350.,235.,"LM3");
  tLM3->SetTextSize(0.035);
  
  TLatex* tLM6 = new TLatex(100.,395.,"LM6");
  tLM6->SetTextSize(0.035);
  
  //  if (tanBeta_ != 50){
  //  LM0->Draw("same");   
  //  tLM0->Draw("same");
  //  LM1->Draw("same");   
  //  tLM1->Draw("same");
  // }

  /*
  if (tanBeta_ == 10){ 
    LM1->Draw("same");
    tLM1->Draw("same");
    LM3->Draw("same");
    tLM3->Draw("same");
    LM6->Draw("same");
    tLM6->Draw("same");
  }
  */



  //stau=LSP contour
  stau->Draw("fsame");
  //  NoEWSB->Draw("fsame");
 
  //legends
  //  legexp->Draw();
  //  legst->Draw();
  //legNoEWSB->Draw();
    myleg->Draw();

  hist->Draw("sameaxis");
  cvsSys->RedrawAxis();
  cvsSys->Update();
  cvsSys->Write();
  
  if( plotLO_ ){
    cvsSys->SaveAs("ExclusionLimit_tanb"+tanb+"_LO.pdf");
    cvsSys->SaveAs("ExclusionLimit_tanb"+tanb+"_LO.png");
  }else{
    cvsSys->SaveAs("ExclusionLimit_tanb"+tanb+"_"+RA2bmode+".eps");
    cvsSys->SaveAs("ExclusionLimit_tanb"+tanb+"_"+RA2bmode+".pdf");
    cvsSys->SaveAs("ExclusionLimit_tanb"+tanb+"_"+RA2bmode+".png");
  }
  
  output->Write();
  //output->Close();
  //delete output; 
  
}


void setPlottingStyle(TH1F& hsig){
  
  hsig.SetStats(kFALSE);
  
  hsig.SetAxisRange(80,750,"Y");
  hsig.SetAxisRange(0,520,"X");
  hsig.SetAxisRange(200,520,"X");

  hsig.GetXaxis()->SetTitle("m_{0} (GeV)");
  hsig.GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hsig.GetYaxis()->SetTitleOffset(0.8);
  hsig.GetYaxis()->SetTitleSize(0.06);
  hsig.GetYaxis()->SetLabelSize(0.06);
  hsig.GetXaxis()->SetTitleOffset(0.9);
  hsig.GetXaxis()->SetTitleSize(0.06);
  hsig.GetXaxis()->SetLabelSize(0.06);

  hsig.SetLineWidth(1);  
  hsig.SetLineColor(kBlue);  
  
}


TGraph* set_sneutrino_d0_1(Int_t tanBeta){
  double sn_m0[14]= {0,  0, 48, 55, 80, 90,100,105,109,105,100, 72, 55,0};
  double sn_m12[14]={0,140,210,220,237,241,242,241,230,220,210,170,150,0};

  TGraph* sn_d0_gr = new TGraph(14,sn_m0,sn_m12);

  sn_d0_gr->SetFillColor(kGreen+3);
  sn_d0_gr->SetFillStyle(1001);

  return sn_d0_gr;
}

TGraph* set_sneutrino_d0_2(Int_t tanBeta){
  double sn_m0[9]= {0, 45, 75,115,130,150,163,185,0};
  double sn_m12[9]={0,140,170,213,202,183,168,140,0};

  TGraph* sn_d0_gr_2 = new TGraph(9,sn_m0,sn_m12);

  sn_d0_gr_2->SetFillColor(kGreen+3);
  sn_d0_gr_2->SetFillStyle(1001);

  return sn_d0_gr_2;
}

TGraph* set_lep_ch(Int_t tanBeta){
  if(tanBeta == 10) return set_lep_ch_tanBeta10();
  if(tanBeta == 40) return set_lep_ch_tanBeta40();
}

TGraph* set_lep_ch_tanBeta10(){

double ch_m0[] ={0  ,100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,2000,0.};  
double ch_m12[]={162,162, 161, 160, 160, 159, 158, 157, 156, 155, 154, 154, 153, 152, 150, 149, 148,  147 , 146 , 146 , 146 , 147 , 148 , 149 , 151 , 154 ,159, 0., 0.}; 

  TGraph* ch_gr = new TGraph(29,ch_m0,ch_m12);

  ch_gr->SetFillColor(3);
  ch_gr->SetLineColor(3);
  //  ch_gr->SetLineWidth(3);
  ch_gr->SetFillStyle(1001);

  return ch_gr;

}

TGraph* set_lep_ch_tanBeta40(){

  double ch_m0[] = {0,   240, 400, 500.,700.,800.,1000,1200., 1300., 1400., 1450., 1500., 1550., 1600., 1650., 1700., 1750., 1800., 2150., 2150.,0};
  double ch_m12[]= {140, 140, 139, 138, 137, 136, 136,  137,   138,  139,   140,   142,   143,   145,   147,   150,   154,   158,  134,   0,   0};

  TGraph* ch_gr = new TGraph(21,ch_m0,ch_m12);

  ch_gr->SetFillColor(3);
  ch_gr->SetLineColor(3);
  ch_gr->SetFillStyle(1001);

  return ch_gr;

}



TGraph* set_tev_stau(Int_t tanBeta){

// taken from Frederic Ronga's calculation

    double st_m0_tanBeta10[] =  {0,   10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110 , 130,  147, 0, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 470, 490, 510, 530, 550, 570, 590, 610, 630, 6500};
    double st_m12_tanBeta10[] = {213,220,240,275,312,351,393,435,476,518, 559, 600., 682., 750.,750, 842., 921., 999., 1076, 1152, 1228, 1304, 1378, 1453, 1527, 1600, 1673, 1746, 1818, 1890, 1962, 2034, 2105, 2175, 2246, 2316, 2386, 2456, 2526, 2595}; 


    double st_m0_tanBeta40[] = {240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 0,440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880};
    double st_m12_tanBeta40[] = {186, 256, 329, 400, 470, 537, 603, 666, 727, 787, 787, 845, 902, 958, 1013, 1067, 1121, 1174, 1226, 1278, 1330, 1381, 1431, 1481, 1531, 1581, 1630, 1679, 1728, 1779, 1825, 1874, 1920, 1971};

    TGraph* st_gr_tanBeta10 = new TGraph(15,st_m0_tanBeta10,st_m12_tanBeta10);
    TGraph* st_gr_tanBeta40 = new TGraph(11,st_m0_tanBeta40,st_m12_tanBeta40);

    
    st_gr_tanBeta40->SetFillColor(40);
    st_gr_tanBeta40->SetFillStyle(1001);
    
    st_gr_tanBeta10->SetFillColor(40);
    st_gr_tanBeta10->SetFillStyle(1001);


    if(tanBeta == 10)return st_gr_tanBeta10;
    if(tanBeta == 40)return st_gr_tanBeta40;
}

TGraph* set_NoEWSB(Int_t tanBeta){

    double st_m0_tanBeta10[]  = { 0, 70, 90,  80,  60, 1000, 1090, 1170, 1240, 1320, 1370, 1440, 1500, 1550, 1610, 1660, 1720, 1780, 1830, 1860, 1920, 1970 , 2000, 2000, 0};
    double st_m12_tanBeta10[] = {10.,10.,20., 30., 40., 50.,  60.,  70.,  80.,  90.,  100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220,   0,   0}; 

//// Needs to be modified
    double st_m0_tanBeta40[] = { 0, 70, 90,  80,  60, 1000, 1090, 1170, 1240, 1320, 1370, 1440, 1500, 1550, 1610, 1660, 1720, 1780, 1830, 1860, 1920, 1970 , 2000, 2000, 0};
    double st_m12_tanBeta40[]= {10.,10.,20., 30., 40., 50.,  60.,  70.,  80.,  90.,  100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220,   0,   0};

    TGraph* st_gr_tanBeta10 = new TGraph(25,st_m0_tanBeta10,st_m12_tanBeta10);
    TGraph* st_gr_tanBeta40 = new TGraph(25,st_m0_tanBeta40,st_m12_tanBeta40);

    
    st_gr_tanBeta40->SetFillColor(40);
    st_gr_tanBeta40->SetFillStyle(1001);
    
    st_gr_tanBeta10->SetFillColor(40);
    st_gr_tanBeta10->SetFillStyle(1001);

    if(tanBeta == 10)return st_gr_tanBeta10;
    if(tanBeta == 40)return st_gr_tanBeta40;
}





TGraph* set_lep_sl(Int_t tanBeta){


  //contour from D0 trilepton paper (PLB 680 (2009) 34-43)

  double *sl_m0 = 0;
  double *sl_m12 = 0;
  int n = 0;

  double sl_m0_3[] ={0,  0, 10, 20, 30, 40, 50, 60, 70, 77,88,95};
  double sl_m12_3[]={0,245,242,239,232,222,209,189,165,140,60,0};
  int n_3 = 12;

  double sl_m0_10[]={ 0,  0, 11, 20, 24, 49, 70, 82,88,90};
  double sl_m12_10[]={0,240,237,233,230,200,150,100,50,0};
  int n_10 = 10;

  if (tanBeta==3){
    sl_m0 = sl_m0_3;
    sl_m12 = sl_m12_3;
    n = n_3;
  }
  //CMS PTDR-II
  //* Selectron_R line mass=99, ISASUGRA7.69, A0=0, m_top=175, tan(beta]=10
  if (tanBeta==10 || tanBeta==50){
    sl_m0 = sl_m0_10;
    sl_m12 = sl_m12_10;
    n = n_10;
  }

  TGraph* lep_sl = new TGraph(n,sl_m0,sl_m12);

  lep_sl->SetFillColor(5);
  lep_sl->SetLineColor(5);
  lep_sl->SetFillStyle(1001);
  
  return lep_sl;
}


TGraph* set_tev_sg_cdf(Int_t tanBeta){

  //New CHF from CDF plot in ICHEP2010 talk (E. Halkiadakis)
  double sg_m0[]= {0,  0, 30, 75,150,185,225,310,360,400,430,500,600,600};
  double sg_m12[]={0,162,168,170,160,150,130,120,109,108,100, 96, 95,  0};
  int np=14;

  TGraph* sg_gr = new TGraph(np,sg_m0,sg_m12);

  //  gStyle->SetHatchesLineWidth(3);

  sg_gr->SetFillColor(2);
  sg_gr->SetLineColor(2);
  //  sg_gr->SetLineWidth(3);
  sg_gr->SetFillStyle(1001); 

  return sg_gr;

}

TGraph* set_tev_sg_d0(Int_t tanBeta){

  double sgd_m0[]= {0,  0, 30, 80,150,240,320,400,500,600,600,0};
  double sgd_m12[]={0,167,166,162,156,138,121,109,105,105,  0,0};
  int npd=12;

  TGraph* sgd_gr = new TGraph(npd,sgd_m0,sgd_m12);

  gStyle->SetHatchesLineWidth(3);

  sgd_gr->SetFillColor(kMagenta+3);
  sgd_gr->SetLineColor(kMagenta+3);
  sgd_gr->SetLineWidth(3);
  sgd_gr->SetFillStyle(3335);

  return sgd_gr;

}



//From Sanjay
TF1* constant_squark(int tanBeta,int i){
//---lines of constant gluino/squark.
// Min squark mass from 1st and 2nd generations using fit for tanbeta = 10.

  double coef1[] = {2.67058e+04, 6.39642e+04, 1.16565e+05, 1.95737e+05, 2.86190e+05};
  double coef2[] = {1.98772e-01, 2.11242e-01, 2.17734e-01, 2.39535e-01, 2.39768e-01};
  double coef3[] = {2.67058e+04, 6.39641e+04, 1.16565e+05, 1.95736e+05, 2.86189e+05};
 
  char hname[200];

  sprintf(hname,"lnsq_%i",i);
  TF1* lnsq = new TF1(hname,"sqrt([0]-x*x*[1]+[2])",0,2000);
  lnsq->SetParameter(0,coef1[i-1]);
  lnsq->SetParameter(1,coef2[i-1]);
  lnsq->SetParameter(2,coef3[i-1]);
  lnsq->SetLineWidth(1);
  lnsq->SetLineColor(kGray);

  return lnsq;
}


TF1* constant_gluino(int tanBeta,int i){
//---lines of constant gluino/squark
  char hname[200];
  sprintf(hname,"lngl_%i",i);

  double coef1[] = {201.77, 311.027, 431.582, 553.895, 676.137};
  double coef2[] = {-0.0146608, -0.01677, -0.022244, -0.0271851, -0.0292212};
   
  TF1* lngl = new TF1(hname,"[0]+x*[1]",0,2000);
  lngl->SetParameter(0,coef1[i-1]);
  lngl->SetParameter(1,coef2[i-1]);
  lngl->SetLineWidth(1);
  lngl->SetLineColor(kGray);

  return lngl;
}

TLatex* constant_squark_text(Int_t it,TF1& lnsq,Int_t tanBeta_){
  char legnm[200];
  it--; //For Sanjay's code

  sprintf(legnm,"#font[92]{#tilde{q}(%i)GeV}",500+250*it);
  Double_t place_x = 170;
  if(tanBeta_ == 40)place_x = 280+70;
  TLatex* t3 = new TLatex(place_x+37*it,lnsq.Eval(-50+place_x+37*it)+5,legnm);
  t3->SetTextSize(0.03);
  t3->SetTextAngle(-15+it*2);
  t3->SetTextColor(kGray+2);

  return t3;
}

TLatex* constant_gluino_text(Int_t it,TF1& lngl){ //, Int_t tanBeta_){
  char legnm[200];
  Int_t tanBeta_ = 10;
  it--; //For Sanjay's code

  sprintf(legnm,"#font[12]{#tilde{g}}#font[92]{(%i)GeV}",500+250*it);

  Double_t place_x = 423;
  Double_t place_y = 18;
  if (tanBeta_ == 10 ) {
    place_x = 1525;
    place_y = -15;
  }
  if (tanBeta_ == 40 ) {
    place_x = 543;
    place_y = 13;
  }
  TLatex* t4 = new TLatex(place_x,25+lngl.Eval(place_x),legnm);
  t4->SetTextSize(0.03);
  t4->SetTextAlign(13);
  t4->SetTextColor(kGray+2);

  return t4;
}

TGraph* constant_mass(double mass,TGraph2D* massGrid)
{
  if(!(mass<massGrid->GetZmax() && mass>massGrid->GetZmin())) return 0;
  TGraph* theGraph = (TGraph*)massGrid->GetContourList(mass)->At(0);
  theGraph->SetLineWidth(1);
  theGraph->SetLineColor(kGray);
  return theGraph;
}

TLatex* constant_squark_text_tanBeta40(int mass,TGraph* massLine){
  if(massLine==0) return 0;
  char legnm[200];
  sprintf(legnm,"#font[92]{#tilde{q}(%i)GeV}",mass);
  Double_t place_x = 180 + 100;
  if(place_x+0.16*mass < 400 || massLine->Eval(place_x+0.16*(mass))+10 <100) return 0;
  TLatex* t3 = new TLatex(place_x+0.16*mass,massLine->Eval(place_x+0.16*(mass))+10,legnm);
  t3->SetTextSize(0.03);
  t3->SetTextAngle(-25);
  t3->SetTextColor(kGray+2);
  return t3;
}

TLatex* constant_gluino_text_tanBeta40(int mass,TGraph* massLine){
  if(!massLine) return 0;
  char legnm[200];
  sprintf(legnm,"#font[12]{#tilde{g}}#font[92]{(%i)GeV}",mass);

  Double_t place_x = 1600;
  if(25+massLine->Eval(place_x) <100) return 0;
  TLatex* t4 = new TLatex(place_x,25+massLine->Eval(place_x),legnm);
  t4->SetTextSize(0.03);
  t4->SetTextAlign(13);
  t4->SetTextColor(kGray+2);

  return t4;
}

TLegend* makeStauLegend(Double_t txtsz,Int_t tanBeta_){
  Double_t ypos_1 = 0.78;
  Double_t ypos_2 = 0.80;
  Double_t xpos_1 = 0.16;
  Double_t xpos_2 = 0.17;
  if(tanBeta_ == 40){
    xpos_1 = 0.14;
    xpos_2 = 0.15;
    ypos_1 = 0.76;
    ypos_2 = 0.78;

  }
    
  TLegend* legst = new TLegend(xpos_1+0.025,ypos_1,xpos_2+0.025,ypos_2);
  legst->SetHeader("#tilde{#tau} = LSP");
  legst->SetFillStyle(0);
  legst->SetBorderSize(0);
  legst->SetTextSize(0.03);
  legst->SetTextAngle(80);

  return legst;
}

TLegend* makeNoEWSBLegend(Double_t txtsz,Int_t tanBeta_){
  Double_t ypos_1 = 0.10+0.02;
  Double_t ypos_2 = 0.20+0.02;
  Double_t xpos_1 = 0.82;
  Double_t xpos_2 = 0.92;
  if(tanBeta_ == 40){
    xpos_1 = 0.10;
    xpos_2 = 0.20;
    ypos_1 = 0.85;
    ypos_2 = 0.95;

  }

  TLegend* legst = new TLegend(xpos_1,ypos_1,xpos_2,ypos_2);
  legst->SetHeader("NoEWSB");
  legst->SetFillStyle(0);
  legst->SetBorderSize(0);
  legst->SetTextSize(0.03);
  legst->SetTextAngle(20);

  return legst;
}


TLegend* makeExpLegend(TGraph& sg_gr, TGraph& sgd_gr,TGraph& ch_gr,TGraph& sl_gr,TGraph& tev_sn,Double_t txtsz,Int_t tanbeta){

  //TLegend* legexp = new TLegend(0.61,0.65,0.91,0.9,NULL,"brNDC");
  //TLegend* legexp = new TLegend(0.70,0.70,0.91,0.9,NULL,"brNDC");
  TLegend* legexp = new TLegend(0.61,0.70,0.99,0.9,NULL,"brNDC");


  legexp->SetFillColor(0);
  legexp->SetShadowColor(0);
  legexp->SetTextSize(txtsz);
  legexp->SetBorderSize(0);

  sg_gr.SetLineColor(1);
  
  legexp->AddEntry(&sg_gr,"CDF  #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=5, #mu<0}","f"); 
  
  legexp->AddEntry(&sgd_gr,"D0   #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=3, #mu<0}","f");  

  ch_gr.SetLineColor(1);
  legexp->AddEntry(&ch_gr,"LEP2   #tilde{#chi}_{1}^{#pm}","f");  
  
  sl_gr.SetLineColor(1);
  if(tanbeta != 40) legexp->AddEntry(&sl_gr,"LEP2   #tilde{#font[12]{l}}^{#pm}","f"); 
  if(tanbeta == 3) legexp->AddEntry(&tev_sn,"D0  #chi^{#pm}_{1}, #chi^{0}_{2}","f");  
 

  return legexp;

}

//////////////
TList* getContours(const TH2* hist,double contourLevel,const TString filename)
{
  TH2* h = (TH2*)hist->Clone("_clone");
  double limitValue[1] = {contourLevel};
  h->SetContour(1,limitValue);
  TCanvas* c = new TCanvas("contour_canvas","Contour List",0,0,600,600);
  h->Draw("CONT LIST");
  c->Update();
  TList* contours = (TList*)((TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours"))->At(0);
  TGraph* contour = (TGraph*)contours->First();
  if(filename!="") 
    {
      for(unsigned int j = 0; j < contours->GetSize(); j++)
	{
	  TString newFilename = filename+"_";
	  newFilename+=j;
	  contour->SaveAs(newFilename+".C");
	  contour = (TGraph*)contours->After(contour); // Get Next graph
	}
    }
  delete h;
  delete c;
  return contours;
}

void smoothHistAcross(TH2* hist,int xMin)
{
  int xMax = hist->GetNbinsX();
  int yMax = hist->GetNbinsY();
  
  for(int xBin = 1; xBin <= xMax; xBin++)
    {
      for(int yBin = 1; yBin <= yMax; yBin++)
	{
	  if(hist->GetBinContent(xBin,yBin)<1 && xBin>=xMin)
	    {
	      if(xBin != xMax && yBin != yMax)
		{
		  if(hist->GetBinContent(xBin + 1,yBin + 1)>0 ||
		     hist->GetBinContent(xBin + 1,yBin)>0 ||
		     hist->GetBinContent(xBin,yBin + 1)>0) hist->SetBinContent(xBin,yBin,1);
		}
	      else if(xBin == xMax)
		{
		  if(hist->GetBinContent(xBin,yBin + 1)>0) hist->SetBinContent(xBin,yBin,1);
		}
	      else if(yBin == yMax)
		{
		  if(hist->GetBinContent(xBin+1,yBin)>0) hist->SetBinContent(xBin,yBin,1);
		}
	    }
	}
    } 
}

void interpolateHistAcross(TH2* hist, int xMin) {
  int xMax = hist->GetNbinsX();
  int yMax = hist->GetNbinsY();
  for(int xBin = xMin; xBin <= xMax; xBin++)  {
    for(int yBin = 1; yBin <= yMax; yBin++) {
      if (hist->GetBinContent(xBin,yBin) <=0 ) {
	double av=0;
	int n=0;
	if (hist->GetBinContent(xBin + 1,yBin + 1)>0 && xBin!=xMax && yBin!=yMax) { av+= hist->GetBinContent(xBin + 1,yBin + 1); ++n;}
	if (hist->GetBinContent(xBin + 1,yBin)    >0 && xBin!=xMax) { av+= hist->GetBinContent(xBin + 1,yBin);     ++n;}
	if (hist->GetBinContent(xBin ,yBin + 1)   >0 && yBin!=yMax) { av+= hist->GetBinContent(xBin ,yBin + 1);    ++n;}
	
	if (n>0)  hist->SetBinContent(xBin,yBin,av/double(n));
      }
    }
  }
}

TH2* smoothHoles(const TH2* originalHist)
{
  TH2* hist = (TH2*)originalHist->Clone("_smoothed");
  int xMax = hist->GetNbinsX();
  int yMax = hist->GetNbinsY();

  int xMin = 0;
  
  for(int xBin = 1; xBin <= xMax; xBin++)
    {
      for(int yBin = 1; yBin <= yMax; yBin++)
	{
	  if(hist->GetBinContent(xBin,yBin)>0)
	    {
	      xMin = xBin;
	      yBin = yMax+1;
	      xBin = xMax+1;
	    }
	}
    } 
  //  for(unsigned int i = 0; i< 1000; i++) smoothHistAcross(hist,xMin);
  for(unsigned int i = 0; i< 1000; i++) interpolateHistAcross(hist,xMin);
  return hist;
}

/*
TGraph* RA2b_limit(const TString filename, const TString histname)
{
  TFile* theFile = TFile::Open(filename);


  TH2* inputHist = (TH2*) theFile->Get(histname);
  
  TH2* inputHistInt = smoothHoles(inputHist);

  TList* contours = getContours(inputHistInt,0.05,"");  //not working at the moment

  TGraph* gr= (TGraph*)(contours->First()->Clone("_save"));

  TFile fout("msugra_debug.root","UPDATE");
  TString id=filename(filename.Index("ge"),9);
  inputHist->SetName(histname+"_"+id);
  inputHist->Write();
  inputHistInt->SetName(histname+"_"+id+"_interpolated");
  inputHistInt->Write();
  gr->SetName("contour_"+id);
  gr->Write();
  fout.Close();

  return gr;
}
*/

TSpline3* RA2b_limit(const TString filename, const TString histname)
{
  TFile* theFile = TFile::Open(filename);

  TSpline3* spline = (TSpline3*) theFile->Get(histname);

  return spline;

}

////////////
