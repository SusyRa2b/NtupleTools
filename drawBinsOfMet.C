//void plotSlices( const TCut btag = "nbjets>=1")
{
  gROOT->SetStyle("CMS");
  
  const TCut btag = "nbjetsSSVHPT>=1";
  //const TCut btag = "nbjetsSSVHPT>=2";
  //const TCut btag = "nbjetsSSVHPT==0";
  //const TCut btag = "nbjetsSSV0>=1";

  const TString var ="minDeltaPhiN";
  //  const TString var ="minDeltaPhi";
  //  const TString var ="bestTopMass";
  //  const TString var ="MET/(MET+HT)";
  //  const TString var ="deltaPhiMPTcaloMET";

  const bool drawLSB = true;
  const bool drawMSB = true;
  const bool drawSB = true;
  const bool drawSIG = true;

  const float customMax = 1;
  //const float customMax = -1;
  const bool doRatio=false; //hack for a specific ratio plot

  const  int nbins=10;
  float min=0;
  float max=0;
  if (var=="bestTopMass") max=800;
  else   if (var=="minDeltaPhi") max=TMath::Pi();
  else   if (var=="minDeltaPhiN") max=20;
  else   if (var=="deltaPhiMPTcaloMET") max=TMath::Pi();
  else max = 0.6; 

  //updated to use reducedTrees

  bool addOverflow = true;
  if (var=="minDeltaPhi") addOverflow=false;
  else  if (var=="bestTopMass") addOverflow=false;
  else  if (var=="deltaPhiMPTcaloMET") addOverflow=false;

  //  gROOT->SetStyle("CMS");
  gStyle->SetOptStat(0);

  TCut baseline ="HT>=350 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && weight<1000"; //no mindp, no MET
  TCut minDPcut = "cutDeltaPhi==1";
  //TCut minDPcut = "minDeltaPhi > 0.2"; cout<<"warning -- using special minDeltaPhi cut! "<<minDPcut.GetTitle()<<endl;
  if (var!="minDeltaPhi" && var!="minDeltaPhiN") baseline = baseline&&minDPcut;
  TCut LSB = "MET<50";
  TCut MSB = "MET>=50 && MET<100";
  TCut SB = "MET>=100 && MET <150";
  TCut SIG = "MET>150";

  TCut cut1=baseline && LSB && btag;
  TCut cut2=baseline && MSB && btag;
  TCut cut3=baseline && SB && btag;
  TCut cut4=baseline && SIG && btag;

  TString selection1 = TString("weight*(")+cut1.GetTitle()+")";
  TString selection2 = TString("weight*(")+cut2.GetTitle()+")";
  TString selection3 = TString("weight*(")+cut3.GetTitle()+")";
  TString selection4 = TString("weight*(")+cut4.GetTitle()+")";

//   TChain madgraph("reducedTree");
//   madgraph.Add("/cu2/joshmt/V00-03-01_6/reducedTree.Baseline0_PF_JERbias6_pfMEThigh_PFLepRA20e0mu_minDP_MuonCleaning.QCD.root");

  TChain pypu("reducedTree");
  //pypu.Add("/cu2/ra2b/reducedTrees/V00-02-05_v2/reducedTree.SSVHPT.PythiaPUQCD.root");
  pypu.Add("/cu2/ra2b/reducedTrees/V00-02-25_fullpf2pat/reducedTree.SSVHPT_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.PythiaPUQCD.root");
  

//   TChain py("reducedTree");
//   py.Add("/cu2/joshmt/V00-03-01_6/reducedTree.Baseline0_PF_JERbias6_pfMEThigh_PFLepRA20e0mu_minDP_MuonCleaning.PythiaQCD.root");

  TH1::SetDefaultSumw2(); //trick to turn on Sumw2 for all histos

  //const int nbins=6;
  //  const double varbins[]={0.,160.,180.,260.,400.,800.,2000};

  int height= doRatio ? 800 : 600;
  TCanvas * thecanvas= new TCanvas("thecanvas","the canvas",700,height);
  if (doRatio) {
    thecanvas->Divide(1,2);
    const float padding=0.01; const float ydivide=0.2;
    thecanvas->GetPad(1)->SetPad( padding, ydivide + padding, 1-padding, 1-padding);
    thecanvas->GetPad(2)->SetPad( padding, padding, 1-padding, ydivide-padding);
    thecanvas->cd(1);
  }
  
  TH1D Hlow("Hlow","",nbins,min,max);
  TH1D Hmed("Hmed","",nbins,min,max);
  TH1D Hmh("Hmh","",nbins,min,max);
  TH1D Hhigh("Hhigh","",nbins,min,max);
  TH1D Hratio("Hratio","",nbins,min,max);
  /*
  TH1D Hlow("Hlow","",nbins,varbins);
  TH1D Hmed("Hmed","",nbins,varbins);
  TH1D Hmh("Hmh","",nbins,varbins);
  TH1D Hhigh("Hhigh","",nbins,varbins);
  */
  Hlow.SetLineColor(kBlue);
  Hmed.SetLineColor(1);
  Hmh.SetLineColor(kRed);
  Hhigh.SetLineColor(6);
  Hlow.SetMarkerColor(kBlue);
  Hmed.SetMarkerColor(1);
  Hmh.SetMarkerColor(kRed);
  Hhigh.SetMarkerColor(6);

  float width=2.5;
  Hlow.SetLineWidth(width);
  Hmed.SetLineWidth(width);
  Hmh.SetLineWidth(width);
  Hhigh.SetLineWidth(width);

  //was minDeltaPhi, bestTopMass
  if (drawLSB) pypu.Draw(var+">>Hlow",selection1);
  if (drawMSB) pypu.Draw(var+">>Hmed",selection2);
  if (drawSB) pypu.Draw(var+">>Hmh",selection3);
  if (drawSIG) pypu.Draw(var+">>Hhigh",selection4);
  gPad->SetRightMargin(0.1);
  gPad->Modified();

  if (addOverflow) {
    Hlow.SetBinContent(nbins, Hlow.GetBinContent(nbins)+Hlow.GetBinContent(nbins+1));
    Hmed.SetBinContent(nbins, Hmed.GetBinContent(nbins)+Hmed.GetBinContent(nbins+1));
    Hmh.SetBinContent(nbins, Hmh.GetBinContent(nbins)+Hmh.GetBinContent(nbins+1));
    Hhigh.SetBinContent(nbins, Hhigh.GetBinContent(nbins)+Hhigh.GetBinContent(nbins+1));
  }

  //  cout<<"SB histogram integral = "<<Hmh.Integral()<<endl;

  if (Hlow.Integral()>0) Hlow.Scale( 1.0 / Hlow.Integral());
  if (Hmed.Integral()>0) Hmed.Scale( 1.0 / Hmed.Integral());
  if (Hmh.Integral()>0) Hmh.Scale( 1.0 / Hmh.Integral());
  if (Hhigh.Integral()>0) Hhigh.Scale( 1.0 / Hhigh.Integral());
  
  if (var=="minDeltaPhi")  Hhigh.SetXTitle("#Delta #phi_{min}");
  else if (var=="bestTopMass")  {
    TString title="best 3-jet mass (GeV)";
    Hhigh.SetXTitle(title);
    Hmh.SetXTitle(title);
    Hmed.SetXTitle(title);
    Hlow.SetXTitle(title);
  }
  else if (var=="minDeltaPhiN") {
    Hhigh.SetXTitle("#Delta #phi_{N}^{min}");
  }
  TString ytitle="Arbitrary units";
  Hhigh.SetYTitle(ytitle);
  Hmh.SetYTitle(ytitle);
  Hmed.SetYTitle(ytitle);
  Hlow.SetYTitle(ytitle);
  //TString thetitle=btag.GetTitle();
  TString thetitle="";
  //  thetitle += " (fail minDeltaPhi)";
  Hhigh.SetTitle(thetitle);
  Hmh.SetTitle(thetitle);
  Hmed.SetTitle(thetitle);
  Hlow.SetTitle(thetitle);

  TString drawopt="hist e";
  if (drawSIG)  {Hhigh.Draw(drawopt); drawopt="hist e SAME"; if (customMax>0) Hhigh.SetMaximum(customMax);}
  if (drawSB)   {Hmh.Draw(drawopt); drawopt="hist e SAME";   if (customMax>0) Hmh.SetMaximum(customMax);}
  if (drawMSB)  {Hmed.Draw(drawopt); drawopt="hist e SAME";  if (customMax>0) Hmed.SetMaximum(customMax);}
  if (drawLSB)  {Hlow.Draw(drawopt); drawopt="hist e SAME";  if (customMax>0) Hlow.SetMaximum(customMax);}

//   if (var=="minDeltaPhi") {
//     Hhigh.SetMinimum(0);
//     Hhigh.GetXaxis()->SetRangeUser(0,2);
//   }

  TLegend leg(0.45,0.55,0.85,0.85);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  leg.SetLineStyle(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(.04);
  if (drawLSB)  leg.AddEntry(&Hlow,"E_{T}^{miss} < 50 GeV");
  if (drawMSB)  leg.AddEntry(&Hmed,"50 < E_{T}^{miss} < 100 GeV");
  if (drawSB)   leg.AddEntry(&Hmh,"100 < E_{T}^{miss} < 150 GeV");
  if (drawSIG)  leg.AddEntry(&Hhigh,"E_{T}^{miss} > 150 GeV");
  leg.Draw();

  if (doRatio) {
    Hratio.Divide(&Hmh,&Hlow);
    Hratio.SetLineWidth(2);
    thecanvas->cd(2);
    Hratio.Draw();
    Hratio.SetMinimum(0);
    Hratio.SetMaximum(3);
    TLine* l1 = new TLine(min,1,max,1);
    l1->SetLineColor(kMagenta);
    l1->SetLineWidth(2);
    l1->Draw();
    cout<<"KS test results = "<<Hmh.KolmogorovTest(&Hlow)<<endl;
  }

//   double chi2=0;
//   for (int i=1; i<=nbins; i++) {
//     double denom=Hlow.GetBinError(i)*Hlow.GetBinError(i) + Hmh.GetBinError(i)*Hmh.GetBinError(i);
//     double c2= denom>0 ? pow( Hlow.GetBinContent(i) - Hmh.GetBinContent(i) ,2) / denom : 0;
//     chi2+=c2;
//   }

//   cout<<"Hand chi^2 = "<<chi2<<endl;
//   Hlow.Chi2Test(&Hmh,"WW p");


  thecanvas->Print(var+".pdf");
}

// plotAll() {

//   TCanvas cplots("cplots","plots",900,600);
//   cplots.Divide(1,3);

//   cplots.cd(1);
//   plotSlices("nbjets==1");

//   cplots.cd(2);
//   plotSlices("nbjets>=1");

//   cplots.cd(3);
//   plotSlices("nbjets>=2");

//   cplots.SaveAs("compM3jsbvslsb.eps");

// }
