#include "TROOT.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TH2.h>
#include <TObjArray.h>
#include "TText.h"
#include "TMath.h"

#include "TCanvas.h"

#include "MiscUtil.cxx"

/*
.L signalEff2012_writetxt.C+

combineScanBins(TString) --  takes an eventcounts file and generates an eventcounts2x2 file with blocks of 4 bins combined

writetxt() -- takes eventcounts files and writes out text files

*/

/* concrete example:
root [0] .L signalEff2012_writetxt.C+
Info in <TUnixSystem::ACLiC>: creating shared library /afs/cern.ch/work/j/joshmt/private/cfaAnalysis/CMSSW_5_2_5/src/NtupleTools/BasicLoopCU/./signalEff2012_writetxt_C.so
root [1] combineScanBins("eventcounts.CSVM_PF2PATjets_JES0_JERbias_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.T1bbbb.root")

-- quit root to be safe, then repeat for JESup and JESdown --
-- then do the writetxt() step --
root [0] .L signalEff2012_writetxt.C+
Info in <TUnixSystem::ACLiC>: creating shared library /afs/cern.ch/work/j/joshmt/private/cfaAnalysis/CMSSW_5_2_5/src/NtupleTools/BasicLoopCU/./signalEff2012_writetxt_C.so
root [1] writetxt("JES","eventcounts2x2.")


*/

void combineScanBins(TString infilename) {

  TFile infile(infilename);

  TH2D* scanSMSngen = (TH2D*) infile.Get("scanSMSngen");
  //first we have to load the histograms from the input file
  vector<TH2D*> vh0;
  for (int ih = 0; ih<infile.GetListOfKeys()->GetEntries(); ih++) {
    TString histname = infile.GetListOfKeys()->At(ih)->GetName();
    if (! (histname.BeginsWith("events_") ||histname.BeginsWith("eventstotal_")) ) continue;
    TH2D* h0 = (TH2D*) infile.Get(histname);
    vh0.push_back(h0);
  }

  //now  new histograms for output
  TString outfilename = infilename;
  outfilename.ReplaceAll("eventcounts.","eventcounts2x2.");
  TFile outfile(outfilename,"RECREATE");


  vector<TH2D*> vh2x2;
  for (unsigned int ih=0; ih<vh0.size(); ih++) {
    TString histname = vh0[ih]->GetName();
    TString newname = histname+"_2x2";
    TH2D* hnew = (TH2D*) vh0[ih]->Clone(newname);    //make the new histograms exact clones of the old
    hnew->Reset();
    vh2x2.push_back(hnew);
  }


  //loop from low to high in both x and y, skipping every other bin in the loop
  //for each bin that we do not skip, combine:
  //  (ix,iy) (ix+1, iy) (ix+1,iy+1) (ix,iy+1)

  for (int ix=1; ix<= scanSMSngen->GetNbinsX(); ix+=2) {
    for (int iy=1; iy<=scanSMSngen->GetNbinsY(); iy+=2) {
      double ngentotal = 
	scanSMSngen->GetBinContent(ix,iy) +
	scanSMSngen->GetBinContent(ix+1,iy) +
	scanSMSngen->GetBinContent(ix,iy+1) +
	scanSMSngen->GetBinContent(ix+1,iy+1);
      scanSMSngen->SetBinContent(ix,iy,ngentotal);
      scanSMSngen->SetBinContent(ix+1,iy,0);
      scanSMSngen->SetBinContent(ix,iy+1,0);
      scanSMSngen->SetBinContent(ix+1,iy+1,0);
      //now loop over all of the event count histos
      for (unsigned int ih=0; ih<vh0.size(); ih++) {
	double ntotal = 
	  vh0[ih]->GetBinContent(ix,iy) +
	  vh0[ih]->GetBinContent(ix+1,iy) +
	  vh0[ih]->GetBinContent(ix,iy+1) +
	  vh0[ih]->GetBinContent(ix+1,iy+1);
	vh2x2[ih]->SetBinContent(ix,iy,ntotal);
      }

    }
  }

  scanSMSngen->Write();
  outfile.Write();
  outfile.Close();

}

void writetxt(TString which, const TString sample, const TString prefix="eventcounts.") {

  assert( which=="counts" || which=="JES"||which=="MET" ||which=="JER");

  //this should be the 'unvaried' (eg JES0) stub.
  //this differentiation between JER0 and JERbias here is a stopgap measure
  TString stub0 = "CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0." ;
  //for JER, need to compare variations with JERbias
  if (which=="JER") stub0= "CSVM_PF2PATjets_JES0_JERbias_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.";

  stub0+=sample;

  TString stub_up = stub0;
  TString stub_down = stub0;

  if (which=="JES") {
    stub_up.ReplaceAll("JES0","JESup");
    stub_down.ReplaceAll("JES0","JESdown");
  }
  else if (which=="MET") {
    stub_up.ReplaceAll("METunc0","METuncUp");
    stub_down.ReplaceAll("METunc0","METuncDown");
  }
  else if (which=="JER") {
    stub_up.ReplaceAll("JERbias","JERup");
    stub_down.ReplaceAll("JERbias","JERdown");
  }

  TString outfilename = (which=="counts") ? "sigcounts." : "sigsystematics.";
  outfilename += stub0.Tokenize(".")->At(1)->GetName();

  if (prefix.Contains("minnjets5")) {
    outfilename+=".minnjets5";
  }

  if (which=="counts")   outfilename+=".txt";
  else if (which=="JES") outfilename += ".JES.txt";
  else if (which=="JER") outfilename += ".JER.txt";
  else if (which=="MET") outfilename += ".MET.txt";


  ofstream txtfile( outfilename.Data());

  TString fn0 = prefix;
  TString fnu = prefix;
  TString fnd = prefix;
  fn0 += stub0;
  fnu += stub_up;
  fnd += stub_down;
  fn0+=".root";
  fnu+=".root";
  fnd+=".root";

  vector<TH2D*> vh0;
  vector<TH2D*> vhu;
  vector<TH2D*> vhd;

  TFile * f0=new TFile(fn0);
  TFile* fu=0;
  TFile* fd=0;
  if (which!="counts") {
    fu=new TFile(fnu);
    fd=new TFile(fnd);
  }

  TH2D* scanSMSngen = (TH2D*) f0->Get("scanSMSngen");


  int nhist=0;
  for (int ih = 0; ih<f0->GetListOfKeys()->GetEntries(); ih++) {
    TString histname = f0->GetListOfKeys()->At(ih)->GetName();
    if (!histname.BeginsWith("events_")) continue;
    nhist++;
    TH2D* h0 = (TH2D*) f0->Get(histname);
    vh0.push_back(h0);
    if (fu!=0) {
      TH2D* hu = (TH2D*) fu->Get(histname);
      TH2D* hd = (TH2D*) fd->Get(histname);
      vhu.push_back(hu);
      vhd.push_back(hd);
    }
  }

  //open a file for output
  outfilename.ReplaceAll(".txt",".root");
  TFile fout(outfilename,"RECREATE");

  vector<TH1D*> vjes0;
  vector<TH1D*> vjesU;
  vector<TH1D*> vjesD;

  //now i've got pointers to all of the histograms
  //loop over mass points
  for (int ix=1; ix<= vh0[0]->GetNbinsX(); ix++) {
    for (int iy=1; iy<=vh0[0]->GetNbinsY(); iy++) {
      
      //zero suppression
      //might want a better check here (use scanSMSngen?)
      //      if (vh0[0]->GetBinContent(ix,iy) ==0) continue;
      if (scanSMSngen->GetBinContent(ix,iy) ==0) continue;

      txtfile<<vh0[0]->GetXaxis()->GetBinLowEdge(ix)<<" "
	     <<vh0[0]->GetYaxis()->GetBinLowEdge(iy)<<" ";

      if (which=="counts") txtfile<<scanSMSngen->GetBinContent(ix,iy)<<" ";

      cout<<" == "<<vh0[0]->GetXaxis()->GetBinLowEdge(ix)<<" "
	  <<vh0[0]->GetYaxis()->GetBinLowEdge(iy)<<" =="<<endl;

      //make a histogram for sanity-check purposes
      TString hname; 
      hname.Form("h%s_%d_%d",which.Data(),TMath::Nint(vh0[0]->GetXaxis()->GetBinLowEdge(ix)),TMath::Nint(vh0[0]->GetYaxis()->GetBinLowEdge(iy)));
      TH1D* hjes0 = new TH1D(hname+"_0",hname,nhist,0,nhist);
      hjes0->Sumw2();
      vjes0.push_back(hjes0);
      TH1D* hjesU=0;      TH1D* hjesD=0;
      if (fu!=0) {
	hjesU = new TH1D(hname+"_Up",hname,nhist,0,nhist); 
	hjesD = new TH1D(hname+"_Down",hname,nhist,0,nhist);
	hjesU->Sumw2();
	hjesD->Sumw2();
	vjesU.push_back(hjesU);
	vjesD.push_back(hjesD);
      }

      //loop over histograms
      for (int ih=0; ih<nhist; ih++) {
	//	cout<<ih<<" "<<vh0[ih]->GetName()<<" ";

	if (fu!=0) {
	  //signed symmetrization of the up and down errors
	  const double cutoff=10; //enforce some sanity check to get rid of NaNs and crazy values
	  bool notok= vhu[ih]->GetBinContent(ix,iy)<cutoff || vhd[ih]->GetBinContent(ix,iy)<cutoff || vh0[ih]->GetBinContent(ix,iy)<cutoff;
	  double delta = notok ? 0 : 0.5*(vhu[ih]->GetBinContent(ix,iy) - vhd[ih]->GetBinContent(ix,iy)) / vh0[ih]->GetBinContent(ix,iy);
	  
	  //unsymmetrized
	  //	  double deltau = vhu[ih]->GetBinContent(ix,iy) / vh0[ih]->GetBinContent(ix,iy) - 1.0;	
	  //	  double deltad = vhd[ih]->GetBinContent(ix,iy) / vh0[ih]->GetBinContent(ix,iy) - 1.0;	
	  //	  cout<<vh0[ih]->GetBinContent(ix,iy)<<" "<<deltau<<" "<<deltad<<" "<<delta<<"\n";
	  txtfile<<delta<<" ";
	  //if we're merging the 3 b bins, then repeat the number twice more
	  if (prefix.Contains("mergebbins")) txtfile<<delta<<" "<<delta<<" ";

	  hjesU->SetBinContent(ih+1,vhu[ih]->GetBinContent(ix,iy));
	  hjesD->SetBinContent(ih+1,vhd[ih]->GetBinContent(ix,iy));
	  hjesU->SetBinError(ih+1,vhu[ih]->GetBinError(ix,iy));
	  hjesD->SetBinError(ih+1,vhd[ih]->GetBinError(ix,iy));
	}
	else {
	  txtfile<< vh0[ih]->GetBinContent(ix,iy)<<" ";
	}
	hjes0->SetBinContent(ih+1,vh0[ih]->GetBinContent(ix,iy));
	hjes0->SetBinError(ih+1,vh0[ih]->GetBinError(ix,iy));

      }

      if (which=="counts") { //need to also print the stat errors
	for (int ih=0; ih<nhist; ih++) {
	  txtfile<< vh0[ih]->GetBinError(ix,iy)<<" ";
	}
      }

      txtfile<<endl;
      //      cout<<endl;
    }
  }

  fout.Write();
  fout.Close();

}

TCanvas * thecanvas;
TH1D* hJES_up;
TH1D* hJES_down;
TH1D* hJES_sym;
TH1D* hJES_pretty;
TText* text;
void drawstuff(TString pointstring = "1100_700", TString what="JES",const TString sample="T1bbbb") {

  TString outfilename = what;
  outfilename += "syst_"; outfilename+=sample; outfilename+="_";
  outfilename+=pointstring; outfilename+=".pdf";

  TFile f("sigsystematics."+sample+"."+what+".root");

  TString histonamebase ="h"+what+"_";
  histonamebase+=pointstring;

  TH1D* hJES_0=(TH1D*) f.Get(histonamebase+"_0");
  TH1D* hJES_Up=(TH1D*) f.Get(histonamebase+"_Up");
  TH1D* hJES_Down=(TH1D*) f.Get(histonamebase+"_Down");

//   hJES_0->Draw();

//   hJES_Up->SetLineColor(kRed);
//   hJES_Up->Draw("same");
//   hJES_Down->SetLineColor(kBlue);
//   hJES_Down->Draw("same");

//aack, this is worse than awful
  hJES_up = (TH1D*) hJES_0->Clone("hJES_up");
  hJES_down = (TH1D*) hJES_0->Clone("hJES_down");

  hJES_sym = (TH1D*) hJES_0->Clone("hJES_sym");

  hJES_up->Reset();
  hJES_down->Reset();
  hJES_sym->Reset();

  for (int i=1; i<= hJES_0->GetNbinsX(); i++) {

    if ( hJES_0->GetBinContent(i) == 0 ) {
      hJES_up->SetBinContent(i,0);
      hJES_down->SetBinContent(i,0);
      hJES_up->SetBinError(i,0);
      hJES_down->SetBinError(i,0);
      continue;
    }

    hJES_up->SetBinContent(i, hJES_Up->GetBinContent(i) / hJES_0->GetBinContent(i) - 1);
    hJES_down->SetBinContent(i, hJES_Down->GetBinContent(i) / hJES_0->GetBinContent(i) - 1);

    hJES_up->SetBinError(i, jmt::errAoverB( hJES_Up->GetBinContent(i), hJES_Up->GetBinError(i),
					    hJES_0->GetBinContent(i),hJES_0->GetBinError(i)) );
    hJES_down->SetBinError(i, jmt::errAoverB( hJES_Down->GetBinContent(i), hJES_Down->GetBinError(i),
					    hJES_0->GetBinContent(i),hJES_0->GetBinError(i)) );

    hJES_sym->SetBinContent(i,  (hJES_up->GetBinContent(i) - hJES_down->GetBinContent(i))/2);

  }

  gROOT->SetStyle("CMS");

  thecanvas=new TCanvas("thecanvas");
//   hJES_up->SetLineColor(kRed);
//   hJES_down->SetLineColor(kBlue);

  // hJES_sym->DrawCopy();
  //hJES_down->DrawCopy("SAME");

  //need to make a prettier plot; skip the LDP
  hJES_pretty = new TH1D("hJES_pretty",what+" systematics", hJES_0->GetNbinsX()/2 , 0,hJES_0->GetNbinsX()/2 );
  int extrabin=0;
  for ( int ibin = 1; ibin <= hJES_0->GetNbinsX() / 2; ibin++ ) {

    double val = hJES_sym->GetBinContent(ibin);

    hJES_pretty->SetBinContent(ibin + extrabin, val);
    TString binlabel;
    binlabel.Form("M%d_H%d",(ibin-1)/4 +1,(ibin-1)%4 +1);
    hJES_pretty->GetXaxis()->SetBinLabel(ibin+extrabin,binlabel);
    //    if ((ibin-1)%4==0) extrabin++;
  }
  hJES_pretty->GetXaxis()->LabelsOption("v");

  double max=  hJES_pretty->GetMaximum();
  double min=  hJES_pretty->GetMinimum();
  max = fabs(max)>fabs(min) ? fabs(max) : fabs(min);

  hJES_pretty->SetMaximum(1.4*max);
  hJES_pretty->SetMinimum(-1.4*max);

  hJES_pretty->DrawCopy();
  pointstring.ReplaceAll("_",",");
  text = new TText(1,0.83*1.4*max,what+" for "+sample+" "+pointstring);
  text->Draw();

  thecanvas->SaveAs(outfilename);

}
