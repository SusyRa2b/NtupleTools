/* Bill wants both efficiency and best UL for each point. */

//work well for SMS.
//the presentation side for mSugra will be much more difficult

//void combine()
{
  gROOT->SetStyle("CMS");
  gROOT->ForceStyle();

  const  TString pathUl = "/afs/cern.ch/user/o/owen/public/RA2b/an-scanplot-unblind-t1bbbb-withcontam-";
  const TString pathEff = "/afs/cern.ch/user/j/joshmt/public/RA2b/RA2b.T1bbbb.";

  TString stubs[4];
  stubs[0] = "ge1b-loose.root";
  stubs[1] = "ge1b-tight.root";
  stubs[2] = "ge2b-loose.root";
  stubs[3] = "ge2b-tight.root";

  TString stubsEff[4];
  stubsEff[0] = "ge1bLoose.root";
  stubsEff[1] = "ge1bTight.root";
  stubsEff[2] = "ge2bLoose.root";
  stubsEff[3] = "ge2bTight.root";

  TFile * files[4];
  TFile * filesEff[4];

  TH2F * ul[4];
  TH2D * eff[4];
  for (int i=0; i<4; i++) {

    files[i] = new TFile(pathUl+stubs[i]);
    ul[i] = (TH2F*) files[i]->Get("hsusyscanXsecul");

    filesEff[i] = new TFile(pathEff+stubsEff[i]);
    TString hname="efficiency_T1bbbb_";
    TString selection=stubsEff[i];
    selection.ReplaceAll(".root","");
    eff[i] = (TH2D*) filesEff[i]->Get(hname+selection);
  }

  TH2F* bestUL = (TH2F*) ul[0]->Clone("bestUL");
  bestUL->Reset();
  TH2F* bestULlog10 = (TH2F*) ul[0]->Clone("bestULlog10");
  bestULlog10->Reset();

  TH2F* effAtBestUL = (TH2F*) ul[0]->Clone("effAtBestUL");
  effAtBestUL->Reset();
  effAtBestUL->SetTitle("signal efficiency for selection with best UL");

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
	TH2D* bestEff = eff[selectionIsBest-1];
	//need to look up the efficiency without assuming that the histograms have the same binning
	//(i checked, and they do not!)
	//translate the bin numbers into actual coordinates
	double xcoord=	bestUL->GetXaxis()->GetBinCenter(i);
	double ycoord=	bestUL->GetYaxis()->GetBinCenter(j);
	int thisbin=	bestEff->FindBin(xcoord,ycoord); //now look up the bin number on eff histogram
	double thiseff = bestEff->GetBinContent(thisbin);
	//this guy has the same binning by definition
	effAtBestUL->SetBinContent(i,j,thiseff);
      }
    }
  }

  TFile fout("ULcombination.T1bbbb.PL.root","RECREATE");
  bestUL->Write();
  bestULlog10->Write();
  whichIsBest->Write();
  effAtBestUL->Write();
  fout.Close();

  //now presentation ... always the biggest pain
  const float LSPmax=1200;
  const float GLmax=1210;

  const float ypos = 0.92;
  const float xpos = 0.18;

  // ~~~~~~~ eff ~~~~~~
  TCanvas effCanvas("effCanvas","eff canvas",700,500);
  effCanvas.cd()->SetRightMargin(0.14);
  effAtBestUL->SetXTitle("m_{gluino} [GeV]");
  effAtBestUL->SetYTitle("m_{LSP} [GeV]");
  effAtBestUL->GetYaxis()->SetRangeUser(0,LSPmax);
  effAtBestUL->GetXaxis()->SetRangeUser(0,GLmax);
  effAtBestUL->Draw("COLZ");

  TLatex* text1= new TLatex(3.570061,23.08044,"CMS Preliminary");
  text1->SetNDC();
  text1->SetTextAlign(13);
  text1->SetX(xpos);
  text1->SetY(ypos-0.05);
  text1->SetTextFont(42);
  text1->SetTextSizePixels(24);
  text1->Draw();

  TString astring;
  astring.Form("%.0f pb^{-1} at #sqrt{s} = 7 TeV",1143.0);
  text2 = new TLatex(3.570061,23.08044,astring);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(xpos);
  text2->SetY(ypos+0.005);//0.88); //0.005 is a kluge to make things look right. i don't get it
  text2->SetTextFont(42);
  text2->SetTextSizePixels(24);
  text2->Draw();

  effCanvas.SaveAs("effAtBestUL_T1bbbb_PL.eps");
  effCanvas.SaveAs("effAtBestUL_T1bbbb_PL.pdf");

  // ~~~~~~~ best UL ~~~~~~
  TCanvas ulCanvas("ulCanvas","ul canvas",700,500);
  ulCanvas.cd()->SetRightMargin(0.14);
  bestUL->SetXTitle("m_{gluino} [GeV]");
  bestUL->SetYTitle("m_{LSP} [GeV]");
  bestUL->GetYaxis()->SetRangeUser(0,LSPmax);
  bestUL->GetXaxis()->SetRangeUser(0,GLmax);
  bestUL->Draw("COLZ");
  ulCanvas.SetLogz();
  bestUL->SetMinimum(1e-2); //hard coded for T1bbbb PL results

  text1->Draw();
  text2->Draw();
  ulCanvas.SaveAs("bestUL_T1bbbb_PL.eps");
  ulCanvas.SaveAs("bestUL_T1bbbb_PL.pdf");

  // ~~~~~~~ which selection ~~~~~~
  TCanvas whichCanvas("whichCanvas","which canvas",700,500);
  whichCanvas.cd()->SetRightMargin(0.14);
  whichIsBest->SetXTitle("m_{gluino} [GeV]");
  whichIsBest->SetYTitle("m_{LSP} [GeV]");
  whichIsBest->GetYaxis()->SetRangeUser(0,LSPmax);
  whichIsBest->GetXaxis()->SetRangeUser(0,GLmax);
  whichIsBest->Draw("COLZ");

  text1->Draw();
  text2->Draw();

  TLatex * key[4];
  TLatex* key[0]= new TLatex(3.570061,23.08044,"1 = #geq 1b Loose");
  TLatex* key[1]= new TLatex(3.570061,23.08044,"2 = #geq 1b Tight");
  TLatex* key[2]= new TLatex(3.570061,23.08044,"3 = #geq 2b Loose");
  TLatex* key[3]= new TLatex(3.570061,23.08044,"4 = #geq 2b Tight");
  for (int i=0; i<4; i++) {
    key[i]->SetNDC();
    key[i]->SetTextAlign(13);
    key[i]->SetX(xpos);
    key[i]->SetY(ypos-(i+3)*0.05);
    key[i]->SetTextFont(42);
    key[i]->SetTextSizePixels(24);
    key[i]->Draw();
  }

  whichCanvas.SaveAs("bestSelection_T1bbbb_PL.eps");
  whichCanvas.SaveAs("bestSelection_T1bbbb_PL.pdf");

}
