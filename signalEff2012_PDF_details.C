#include "TFile.h"
#include "TH2.h"
#include "TString.h"

#include <vector>
#include <iostream>

void signalEff2012_PDF_details(const TString what="CTEQMSTW",const TString sample="T1bbbb",const int minnjets=3,const bool useisr=false) {

  TString njetsstring=".";

  if (minnjets==3) {
    //do nothing
  }
  else if (minnjets==5) {
    njetsstring = ".minnjets5.";
  }
  else {
    cout<<" minnjets = "<<minnjets<<" is not ok"<<endl;
    return;
  }

  TString stub1="eventcounts2x2.mergebbins";
  if (sample.Contains("pMSSM") ||sample.Contains("T1ttcc") ||sample.Contains("14TeV")||sample.Contains("TChi")) stub1="eventcounts.mergebbins";

  TString stub2=stub1;
  stub2+=".withpdfs";

  if ( useisr) {
    stub1+=".Isr0";
    stub2+=".Isr0";
  }

  //  TString nominalstub="CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0."; //old
  TString nominalstub="JES0_JER0_PFMETTypeI_METunc0_PUunc0_hpt20."; //new
  TString f0file = stub1+njetsstring+nominalstub+sample+".root";
  TString fpdffile = stub2+njetsstring+nominalstub+sample+".root";

  TFile f0(f0file);
  TFile fpdf(fpdffile);

  vector<TH2D*> effratios;

   for (int ih = 0; ih<fpdf.GetListOfKeys()->GetEntries(); ih++) {
     TString histname = fpdf.GetListOfKeys()->At(ih)->GetName();
     if (! (histname.BeginsWith("events_")&&histname.Contains(what))) continue;

     TString histnametotal = histname;
     histnametotal.ReplaceAll("events_","eventstotal_");

     TString histname0=histname;
     histname0.ReplaceAll("_"+what+"0","");

     TString histname0total = histnametotal;
     histname0total.ReplaceAll("_"+what+"0","");


     TH2D* h0c = (TH2D*) f0.Get(histname0);
     TH2D* h0t = (TH2D*) f0.Get(histname0total);
     
     TH2D* hpdfc = (TH2D*) fpdf.Get(histname);
     TH2D* hpdft = (TH2D*) fpdf.Get(histnametotal);
     
     //now calculate efficiency
     TH2D* h0r = (TH2D*) h0c->Clone("h0r_"+histname);
     h0r->Reset();
     h0r->Divide(h0c,h0t);

     TH2D* hpdfr = (TH2D*) hpdfc->Clone("hpdfr_"+histname);
     hpdfr->Reset();
     hpdfr->Divide(hpdfc,hpdft);

     //and then the ratio of efficiencies
     TString rname = histname;
     rname.ReplaceAll("events_","effratio_");
     TH2D*  heffratio = (TH2D*) h0r->Clone(rname);
     heffratio->Reset();
     heffratio->Divide(hpdfr,h0r);
     effratios.push_back(heffratio);

   }

   TFile fout("signalsyst_PDF_"+what+"."+sample+njetsstring+"root","recreate");
   for (unsigned int ii=0; ii<effratios.size();ii++) {
     effratios.at(ii)->Write();
   }
   fout.Close();


   /* quick and dirty version
  TH2D* h0c = (TH2D*) f0.Get("events_b1_HT1000to100000_MET350to100000_2x2");
  TH2D* h0t = (TH2D*) f0.Get("eventstotal_b1_HT1000to100000_MET350to100000_2x2");

  TH2D* h0r = (TH2D*) h0c->Clone("h0r");
  h0r->Reset();
  h0r->Divide(h0c,h0t);

  TH2D* hpdfc = (TH2D*) fpdf.Get("events_b1_HT1000to100000_MET350to100000_CTEQMSTW0_2x2");
  TH2D* hpdft = (TH2D*) fpdf.Get("eventstotal_b1_HT1000to100000_MET350to100000_CTEQMSTW0_2x2");

  TH2D* hpdfr = (TH2D*) hpdfc->Clone("hpdfr");
  hpdfr->Reset();
  hpdfr->Divide(hpdfc,hpdft);


  TH2D*  heffratio = (TH2D*) h0r->Clone("heffratio");
  heffratio->Reset();
  heffratio->Divide(hpdfr,h0r);
  heffratio->Draw("COLZ");
   */


}
