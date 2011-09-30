// use: .x combineScans.C
// setting for public result is aleCLsCustom

//void combine()
{

  const TString which = "aleCLsCustom";
  TString outfileid="";
  if (which=="old") {
    outfileid=   "oldPL";
  }
  else if (which=="owenPL") {
    outfileid=   "alternativePL";
  }
  else if (which=="aleCLs") {
    outfileid=   "officialCLs";
  }
  else if (which=="alePL") {
    outfileid=   "officialPL";
  }
  else if (which=="aleCLsExpected") {
    outfileid=   "officialCLsExpected";
  }
  else if (which=="aleCLsObsExp") { //show observed, based on expected
    outfileid=   "officialCLsObsExp";
  }
  else if (which=="aleCLsCustom") {
    outfileid=   "officialCLsCustom";
  }
  else {assert(0);}

  gROOT->SetStyle("CMS");
  gROOT->ForceStyle();

  TFile fcurve("htCurveT1bbbb.root");
  TH2D* htfrac =  (TH2D*) fcurve.Get("genbHTfrac");
  TGraph* htcurve = (TGraph*) fcurve.Get("htCurveT1bbbb");
  TLine* diagonal = new TLine(300,300,1200,1200);
  diagonal->SetLineStyle(2);

  TFile fmaria("/afs/cern.ch/user/d/dalfonso/public/myFileGraph_RA2b.root");
  TGraph* gr10 = (TGraph*) fmaria.Get("graph_h_RA2b_best_T1bbbb_cloned_1.000000");
  TGraph* gr03 = (TGraph*) fmaria.Get("graph_h_RA2b_best_T1bbbb_cloned_0.333333");
  TGraph* gr30 = (TGraph*) fmaria.Get("graph_h_RA2b_best_T1bbbb_cloned_3.000000");
  gr10->SetMarkerSize(0);
  gr03->SetMarkerSize(0);
  gr30->SetMarkerSize(0);

  TString pathUl="";
  TString pathUlExp="";
  if (which=="old")      pathUl = "/afs/cern.ch/user/o/owen/public/RA2b/an-scanplot-unblind-t1bbbb-withcontam-"; //old
  else if (which=="owenPL") pathUl = "/afs/cern.ch/user/o/owen/public/RA2b/t1bbbb-all-plots.root"; //new
  else if (which=="aleCLs") pathUl = "/afs/cern.ch/user/g/gaz/public/T1bbbb_xsULs.root";
  else if (which=="aleCLsCustom") pathUl = "/afs/cern.ch/user/g/gaz/public/T1bbbb_xsULs.root";
  else if (which=="aleCLsObsExp") pathUl = "/afs/cern.ch/user/g/gaz/public/T1bbbb_xsULs.root";
  else if (which=="aleCLsExpected") pathUl = "/afs/cern.ch/user/g/gaz/public/T1bbbb_xsULs_expected.root";
  else if (which=="alePL") pathUl = "/afs/cern.ch/user/g/gaz/public/T1bbbb_PLxsULs.root";
  const TString pathEff = "/afs/cern.ch/user/j/joshmt/public/RA2b/RA2b.T1bbbb.";

  if (which=="aleCLsObsExp") pathUlExp = "/afs/cern.ch/user/g/gaz/public/T1bbbb_xsULs_expected.root";

  TString stubs[4];
  if (which=="old") {
    stubs[0] = "ge1b-loose.root";
    stubs[1] = "ge1b-tight.root";
    stubs[2] = "ge2b-loose.root";
    stubs[3] = "ge2b-tight.root";
  }
  else if (which=="owenPL") {
    //new
    stubs[0] = "ge1bloose";
    stubs[1] = "ge1btight";
    stubs[2] = "ge2bloose";
    stubs[3] = "ge2btight";
  }
  else if (which=="aleCLs" || which=="alePL" ||which=="aleCLsExpected" ||which=="aleCLsCustom" || which=="aleCLsObsExp") {
    //new
    stubs[0] = "1bloose";
    stubs[1] = "1btight";
    stubs[2] = "2bloose";
    stubs[3] = "2btight";
  }

  TString stubsEff[4];
  stubsEff[0] = "ge1bLoose.root";
  stubsEff[1] = "ge1bTight.root";
  stubsEff[2] = "ge2bLoose.root";
  stubsEff[3] = "ge2bTight.root";

  TFile * files[4];
  TFile * filesEff[4];
  TFile * filesExp[4];

  //new
  if (which!="old")  files[0] = new TFile(pathUl);

  TH2F * ul[4];
  TH2F * ulExp[4];
  TH2D * eff[4];
  for (int i=0; i<4; i++) {

    if (which=="old")    files[i] = new TFile(pathUl+stubs[i]); //old
    else        files[i] = files[0];    //new -- only 1 file

    if (which=="old")    ul[i] = (TH2F*) files[i]->Get("hsusyscanXsecul");
    else        {
      TString histoname="";
      if (which=="owenPL") {
	histoname = stubs[i];
	histoname += "_hplxsecul";
      }
      else if (which=="aleCLs" || which=="aleCLsCustom"|| which=="aleCLsObsExp") {
	histoname = "h2d_";
	histoname+=stubs[i];
	histoname+="_xsUL";
      }
      else if (which=="aleCLsExpected") {
	histoname = "h2d_";
	histoname+=stubs[i];
	histoname+="_xsUL_exp";
      }
      else if (which=="alePL") {
	histoname = "h2d_";
	histoname+=stubs[i];
	histoname+="_PLxsUL";
      }
      ul[i] = (TH2F*) files[i]->Get(histoname);

      if ( which=="aleCLsObsExp") {
	histoname = "h2d_";
	histoname+=stubs[i];
	histoname+="_xsUL_exp";
	if (i==0) filesExp[i]= new TFile(pathUlExp);
	else filesExp[i]=filesExp[0];

	ulExp[i] = (TH2F*) filesExp[i]->Get(histoname);
      }
      else {
	filesExp[i]=0;
	ulExp[i]=0;
      }
    }

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
  bestUL->SetTitle("best upper limit");
  bestULlog10->SetTitle("log10 of best upper limit");

  TH2F* effAtBestUL = (TH2F*) ul[0]->Clone("effAtBestUL");
  effAtBestUL->Reset();
  effAtBestUL->SetTitle("signal efficiency for selection with best UL");

  TH2F* whichIsBest = (TH2F*) ul[0]->Clone("whichIsBest");
  whichIsBest->Reset();
  whichIsBest->SetTitle("Selection with best UL");

  TFile * whichLimitFile = 0;
  if (which=="aleCLsCustom") {
    //replace the whichIsBest histo with Ale's histo
    whichLimitFile = new TFile("/afs/cern.ch/user/g/gaz/public/best_selection_expected.root");
    whichIsBest = (TH2F*) whichLimitFile->Get("h2d_best_sel_two");
    //cout<<"wib = "<<whichIsBest<<endl;
  }

  for (int i=1; i<=bestUL->GetXaxis()->GetNbins(); i++) {
    for (int j=1; j<=bestUL->GetYaxis()->GetNbins(); j++) {
      double minul=1e9;
      int selectionIsBest=0;
      if (whichLimitFile==0) {
	for (int ifile=0; ifile<4; ifile++) {
	  //	  double thisul = ul[ifile]->GetBinContent(i,j);
	  double thisul = (ulExp[ifile]==0) ? ul[ifile]->GetBinContent(i,j) : ulExp[ifile]->GetBinContent(i,j);
	  if (thisul<minul && thisul>0) { minul=thisul; selectionIsBest = ifile+1;}
	}
      }
      else { //don't find the best, because the external file is already telling us which one is best
	selectionIsBest = whichIsBest->GetBinContent(i,j);
	if (selectionIsBest>0)	minul = ul[selectionIsBest-1]->GetBinContent(i,j);
      }
      if (selectionIsBest!=0) {
	minul = ul[selectionIsBest-1]->GetBinContent(i,j);
	bestULlog10->SetBinContent(i,j,log10(minul));
	bestUL->SetBinContent(i,j,minul);
	//	cout<<"selectionIsBest = "<<selectionIsBest<<endl;
	whichIsBest->SetBinContent(i,j,selectionIsBest);
	TH2D* bestEff = eff[selectionIsBest-1];
	//need to look up the efficiency without assuming that the histograms have the same binning
	//(i checked, and they do not!)
	//translate the bin numbers into actual coordinates
	double xcoord=	bestUL->GetXaxis()->GetBinCenter(i);
	double ycoord=	bestUL->GetYaxis()->GetBinCenter(j);
	int thisbin=	bestEff->FindBin(xcoord,ycoord); //now look up the bin number on eff histogram
	double thiseff = 100*bestEff->GetBinContent(thisbin);
	//this guy has the same binning by definition
	effAtBestUL->SetBinContent(i,j,thiseff);
      }
    }
  }


  //kill off the 'forbidden' bins according to Gala's histo

  for (int i=1; i<=bestUL->GetXaxis()->GetNbins(); i++) {
    for (int j=1; j<=bestUL->GetYaxis()->GetNbins(); j++) {
      double x1,x2,y1,y2;
      //find coordinates of all 4 corners
      x1=  bestUL->GetXaxis()->GetBinLowEdge(i);
      x2=  bestUL->GetXaxis()->GetBinLowEdge(i+1);
      y1=  bestUL->GetXaxis()->GetBinLowEdge(j);
      y2=  bestUL->GetXaxis()->GetBinLowEdge(j+1);

      //check ht frac in all for corners
      //i'm going to all this work because i think the histos have different binning
      double  f1=  htfrac->GetBinContent( htfrac->FindBin(x1,y1));
      double  f2=  htfrac->GetBinContent( htfrac->FindBin(x1,y2));
      double  f3=  htfrac->GetBinContent( htfrac->FindBin(x2,y1));
      double  f4=  htfrac->GetBinContent( htfrac->FindBin(x2,y2));

      //      cout<<f1<<"\t"<<f2<<"\t"<<f3<<"\t"<<f4<<endl;
      const double cutval = 0.2;
      //      if (f1<cutval || f2<cutval*/||f3<cutval||f4<cutval) {
      if (f3<cutval) { //seems to work best to reproduce Mariarosaria's cut for MT2b
	bestUL->SetBinContent(i,j,0);
	whichIsBest->SetBinContent(i,j,0);
	effAtBestUL->SetBinContent(i,j,0);
      }
      else {
	if (fabs(x1-y1)<250) {
	  cout<<"NOT Throwing out point with DeltaM = "<<x1-y1<<endl;
	  cout<<x1<<"\t"<<y1<<endl;
	}
      }
    }
  }

  //max is still a holdover from the cloned histo
  //this trick seems effective at finding the real max
  whichIsBest->SetMaximum(4);
  whichIsBest->SetMinimum(0);
  effAtBestUL->SetMaximum( effAtBestUL->GetBinContent(effAtBestUL->GetMaximumBin()));
  effAtBestUL->SetMinimum(0);
  bestUL->SetMinimum(1e-2); //hard coded for T1bbbb PL results
  bestUL->SetMaximum( bestUL->GetBinContent(bestUL->GetMaximumBin()));

  TFile fout("ULcombination.T1bbbb."+outfileid+".root","RECREATE");
  bestUL->Write();
  bestULlog10->Write();
  whichIsBest->Write();
  effAtBestUL->Write();
  fout.Close();

  //now presentation ... always the biggest pain
  const float LSPmax=1200;
  const float GLmax=1210;
  const float GLmin=300;

  //we were using  a "wide" aspect ratio at the request of Bill et al
  //ARC asked for standard square format
  //here s a switch to go between the two
  bool wide = false;

  const float ypos = wide ? 0.92 : 0.88;
  const float xpos = 0.18;

  const int width= wide? 900 : 700;
  const int height=wide? 500 : 650;

  const float rightmarginZ = wide ? 0.18 : 0.24;

  TString legT1bbbb="pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow 2b + LSP; m(#tilde{q})>>m(#tilde{g})"; 
  // ~~~~~~~ eff ~~~~~~
  gStyle->SetPaintTextFormat("2.0f");
  TCanvas effCanvas("effCanvas","eff canvas",width,height);
  if (!wide)  effCanvas.cd()->SetTopMargin(0.1);
  //  effCanvas.cd()->SetRightMargin(0.14);
  effCanvas.cd()->SetRightMargin(rightmarginZ); //14
  effAtBestUL->SetXTitle("m_{#tilde{g}} [GeV]");
  effAtBestUL->SetYTitle("m_{LSP} [GeV]");
  effAtBestUL->SetZTitle("A #times #varepsilon [%]");
  effAtBestUL->GetYaxis()->SetRangeUser(0,LSPmax);
  effAtBestUL->GetXaxis()->SetRangeUser(GLmin,GLmax);
  if (!wide) effAtBestUL->GetXaxis()->SetNdivisions(505);
  effAtBestUL->Draw("COLZ");
  effAtBestUL->Draw("text same");
  diagonal->Draw("L");

  TLatex* text1= new TLatex(3.570061,23.08044,"CMS Preliminary");
  text1->SetNDC();
  text1->SetTextAlign(13);
  text1->SetX(xpos);
  text1->SetY(ypos-0.06);
  text1->SetTextFont(42);
  text1->SetTextSizePixels(23);
  text1->Draw();

  TString astring;
  astring.Form("L_{int} = %.1f fb^{-1}, #sqrt{s} = 7 TeV",1.1);
  text2 = new TLatex(3.570061,23.08044,astring);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(xpos);
  text2->SetY(ypos+0.005);//0.88); //0.005 is a kluge to make things look right. i don't get it
  text2->SetTextFont(42);
  text2->SetTextSizePixels(23);
  text2->Draw();

  TLatex* text3= new TLatex(3.570061,23.08044,legT1bbbb);
  text3->SetNDC();
  text3->SetTextAlign(13);
  text3->SetX(xpos);
  text3->SetY(ypos+0.07);
  text3->SetTextFont(42);
  text3->SetTextSize(4e-2);
  text3->Draw();


  effCanvas.SaveAs("effAtBestUL_T1bbbb_"+outfileid+".eps");
  effCanvas.SaveAs("effAtBestUL_T1bbbb_"+outfileid+".pdf");

  // ~~~~~~~ best UL ~~~~~~
  gStyle->SetPaintTextFormat("3.2f");
  TCanvas ulCanvas("ulCanvas","ul canvas",width,height);
  if (!wide)  ulCanvas.cd()->SetTopMargin(0.1);
  ulCanvas.cd()->SetRightMargin(rightmarginZ); //14
  bestUL->SetXTitle("m_{#tilde{g}} [GeV]");
  bestUL->SetYTitle("m_{LSP} [GeV]");
  bestUL->SetZTitle("Cross section UL at 95% CL [pb]");
  bestUL->GetYaxis()->SetRangeUser(0,LSPmax);
  bestUL->GetXaxis()->SetRangeUser(GLmin,GLmax);
  if (!wide)  bestUL->GetZaxis()->SetTitleOffset(1.3);
  if (!wide) bestUL->GetXaxis()->SetNdivisions(505);
  bestUL->SetMaximum(40);
  bestUL->Draw("COLZ");
  bestUL->SetMarkerSize(0.8);
  bestUL->Draw("text same");
  diagonal->Draw("L");

  ulCanvas.SetLogz();

  text1->Draw();
  text2->Draw();
  text3->Draw();

  gr10->SetLineWidth(3);
  gr30->SetLineWidth(3);
  gr03->SetLineWidth(3);

  gr30->SetLineStyle(3);
  gr03->SetLineStyle(2);

  gr10->Draw("L");
  gr03->Draw("L");
  gr30->Draw("L");

  TLegend smsleg(xpos,ypos-0.29,xpos+0.2,ypos-0.1);
  smsleg.AddEntry(gr10,"#sigma^{prod} = #sigma^{NLO-QCD}");
  smsleg.AddEntry(gr30,"#sigma^{prod} = 3 #times #sigma^{NLO-QCD}");
  smsleg.AddEntry(gr03,"#sigma^{prod} = 1/3 #times #sigma^{NLO-QCD}");
  smsleg.SetFillColor(0); 
  smsleg.SetShadowColor(0);
  smsleg.SetTextSize(0.03);
  smsleg.SetBorderSize(0);
  smsleg.Draw();

  ulCanvas.SaveAs("bestUL_T1bbbb_"+outfileid+".eps");
  ulCanvas.SaveAs("bestUL_T1bbbb_"+outfileid+".pdf");

  // ~~~~~~~ which selection ~~~~~~
  gStyle->SetPaintTextFormat("1.0f");
  TCanvas whichCanvas("whichCanvas","which canvas",width,height);
  if (!wide)  whichCanvas.cd()->SetTopMargin(0.1);
  whichCanvas.cd()->SetRightMargin(0.14);
  whichIsBest->SetXTitle("m_{#tilde{g}} [GeV]");
  whichIsBest->SetYTitle("m_{LSP} [GeV]");
  whichIsBest->GetYaxis()->SetRangeUser(0,LSPmax);
  whichIsBest->GetXaxis()->SetRangeUser(GLmin,GLmax);
  if (!wide) whichIsBest->GetXaxis()->SetNdivisions(505);
  whichIsBest->Draw("COL");
  if (which!="aleCLsCustom")  whichIsBest->Draw("text same");
  diagonal->Draw("L");

  //draw the "David Stuart" code on the plot.
  // this is only relevant for aleCLsCustom mode

  TPaveText * pave=0;
  if (which=="aleCLsCustom") {

    for (int i= 1; i<=whichIsBest->GetXaxis()->GetNbins(); i++) {
      for (int j= 1; j<=whichIsBest->GetYaxis()->GetNbins(); j++) {
	if (whichIsBest->GetBinContent(i,j) >=1 && whichIsBest->GetBinContent(i,j)<=4) {
	  if (whichIsBest->GetXaxis()->GetBinLowEdge(i) >=GLmin) {
	    pave= new TPaveText(whichIsBest->GetXaxis()->GetBinLowEdge(i),whichIsBest->GetYaxis()->GetBinLowEdge(j),
				whichIsBest->GetXaxis()->GetBinLowEdge(i+1),whichIsBest->GetYaxis()->GetBinLowEdge(j+1),"");
	    TString code = ( whichIsBest->GetBinContent(i,j) == 1 || whichIsBest->GetBinContent(i,j) == 2) ? "1" : "2";
	    code += ( whichIsBest->GetBinContent(i,j) == 1 || whichIsBest->GetBinContent(i,j) == 3) ? "L" : "T";
	    pave->AddText(code);
	    pave->SetTextAlign(22); //centered in both vert and horz
	    pave->SetFillColor(kWhite);
	    pave->SetShadowColor(0);
	    pave->SetTextSize(0.025);
	    pave->SetBorderSize(0);
	    pave->SetFillStyle(0);
	    pave->Draw();
	  }
	}
      }
    }
  }

  text1->Draw();
  text2->Draw();
  text3->Draw();

  TLatex * key[4];
  if (which!="aleCLsCustom") {
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
      key[i]->SetTextSizePixels(23);
      key[i]->Draw();
    }
  }

  whichCanvas.SaveAs("bestSelection_T1bbbb_"+outfileid+".eps");
  whichCanvas.SaveAs("bestSelection_T1bbbb_"+outfileid+".pdf");

}
