#include "TROOT.h"
#include "TFile.h"
#include "TH2.h"
#include <vector>
#include <fstream>
#include "TString.h"
#include <iostream>

#include "TCanvas.h"

bool isNormalBin(const TString & name) {

  if (name.Contains("LDP")) return false;
  if (name.Contains("SL")) return false;//catches SLsig as well

  return true;
}

void signalEff2012_PDF_details_writetxt(const TString sample,const int minnjets=3) {

  //assumes that b tag bins are merged

  TString njetsstring=".";
  if (minnjets==3) { 
    //do nothing 
  }
  else if (minnjets==5) {
    njetsstring = ".minnjets5.";
  }
  else {  
    cout<<" minnjets = "<<minnjets<<" is not ok"<<endl; return; 
  }

  TFile f1("signalsyst_PDF_CTEQMSTW."+sample+njetsstring+"root");
  TFile f2("signalsyst_PDF_CTEQNNPDF."+sample+njetsstring+"root");

  vector<TH2D*> vh1;
  vector<TH2D*> vh2;


  for (int ih = 0; ih<f1.GetListOfKeys()->GetEntries(); ih++) {
    TString histname = f1.GetListOfKeys()->At(ih)->GetName();
    vh1.push_back((TH2D*) f1.Get(histname));
  }
  //horrible copy and paste
  for (int ih = 0; ih<f2.GetListOfKeys()->GetEntries(); ih++) {
    TString histname = f2.GetListOfKeys()->At(ih)->GetName();
    vh2.push_back((TH2D*) f2.Get(histname));
  }
  assert(vh1.size() == vh2.size());

  ofstream txtfile( "sigsystematics."+sample+njetsstring+"PDF.txt");

 
  gROOT->SetStyle("CMS");//plotting


  //now i've got pointers to all of the histograms
  //loop over mass points
  for (int ix=1; ix<= vh1[0]->GetNbinsX(); ix++) {
    for (int iy=1; iy<=vh1[0]->GetNbinsY(); iy++) {

      //some zero suppression -- but need to really look through all of the bins
      bool bad=true;
      for (unsigned int ih=0; ih<vh1.size(); ih++) {
	if (vh1[ih]->GetBinContent(ix,iy) >0.01) bad=false;
      }
      if (bad) continue;

      TCanvas*  thecanvas=new TCanvas("thecanvas"); //plotting
      int nbinsforplot=0;
      for (unsigned int ih=0; ih<vh1.size(); ih++) {
	if ( isNormalBin(vh1.at(ih)->GetName())) nbinsforplot++;
      }
      TH1D*  hpretty = new TH1D("hpretty","PDF systematics", nbinsforplot ,0,nbinsforplot); //plotting


      txtfile<<vh1[0]->GetXaxis()->GetBinLowEdge(ix)<<" "
	     <<vh1[0]->GetYaxis()->GetBinLowEdge(iy)<<" ";
     cout<<" == "<<vh1[0]->GetXaxis()->GetBinLowEdge(ix)<<" "
	  <<vh1[0]->GetYaxis()->GetBinLowEdge(iy)<<" =="<<endl;

      //for this mass point, compute an error somehow for each bin
      for (unsigned int ih=0; ih<vh1.size(); ih++) {

	double syst1 = vh1[ih]->GetBinContent(ix,iy);
	double syst2 = vh2[ih]->GetBinContent(ix,iy);

	double syst=0;

	//check for bogus values
	if (syst1>0.1 && syst2>0.1) {
	  syst1 = syst1-1;
	  syst2 = syst2-1;
	  syst= (syst1+syst2) / 2.0;
	}
	  //	if ( fabs(syst1)>fabs(syst2) ) syst= syst1;
	  //	else syst=syst2;

	//ignore syst2 for the moment
	txtfile<<syst<<" ";
	txtfile<<syst<<" ";
	txtfile<<syst<<" ";

	cout<<syst<<" ";
	//logic here assumes that the *first* bins are not LDP
	if (isNormalBin(vh1.at(ih)->GetName())) {
	  hpretty->SetBinContent(ih+1,syst);
	  TString name1=TString(vh1[ih]->GetName()).Tokenize("_")->At(2)->GetName();
	  TString name2=TString(vh1[ih]->GetName()).Tokenize("_")->At(3)->GetName();
	  hpretty->GetXaxis()->SetBinLabel(ih+1,name1+" "+name2);
	}
      }
      cout<<endl;
      txtfile<<endl;
      hpretty->GetXaxis()->LabelsOption("v");
      hpretty->DrawCopy();
      TString fname = "pdfsyst_"; fname+=sample; fname+=njetsstring; 
      fname+= vh1[0]->GetXaxis()->GetBinLowEdge(ix);
      fname+=vh1[0]->GetYaxis()->GetBinLowEdge(iy);
      fname+=".png";

      //      if (vh1[0]->GetXaxis()->GetBinLowEdge(ix) >800 && vh1[0]->GetXaxis()->GetBinLowEdge(ix)<1400)	thecanvas->SaveAs(fname);
      delete hpretty;
      delete thecanvas;


    }

  }

}
