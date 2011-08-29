//void combine()
{

  TString stubs[4];
  stubs[0] = "ge1b-loose.root";
  stubs[1] = "ge1b-tight.root";
  stubs[2] = "ge2b-loose.root";
  stubs[3] = "ge2b-tight.root";

  TFile * files[4];

  TH2F * ul[4];
  for (int i=0; i<4; i++) {

    files[i] = new TFile(TString("/afs/cern.ch/user/o/owen/public/RA2b/an-scanplot-unblind-t1bbbb-withcontam-")+stubs[i]);
    ul[i] = (TH2F*) files[i]->Get("hsusyscanXsecul");
  }

  TH2F* bestUL = (TH2F*) ul[0]->Clone("bestUL");
  bestUL->Reset();
  TH2F* bestULlog10 = (TH2F*) ul[0]->Clone("bestULlog10");
  bestULlog10->Reset();

  TH2F* whichIsBest = (TH2F*) ul[0]->Clone("whichIsBest");
  whichIsBest->Reset();
  whichIsBest->SetTitle("Selection with best UL");

  for (int i=1; i<=bestUL->GetXaxis()->GetNbins(); i++) {
    for (int j=1; j<=bestUL->GetYaxis()->GetNbins(); j++) {
      double minul=1e9;
      int selectionIsBest=0;
      for (int ifile=0; ifile<4; ifile++) {
	double thisul = ul[ifile]->GetBinContent(i,j);
	if (thisul<minul && thisul>0) { minul=thisul; selectionIsBest = ifile+1;}
      }
      if (selectionIsBest!=0) {
	bestULlog10->SetBinContent(i,j,log10(minul));
	bestUL->SetBinContent(i,j,minul);
	whichIsBest->SetBinContent(i,j,selectionIsBest);
      }
    }
  }

  TFile fout("ULcombination.T1bbbb.PL.root","RECREATE");
  bestUL->Write();
  bestULlog10->Write();
  whichIsBest->Write();
  fout.Close();

}
