#include "drawDelphesBase.C"
/*
--- to compile ---

.L TSelectorMultiDraw.C+
.L CrossSectionTable.cxx+
.L ConfigurationDescriptions.cxx+

.L drawDelphes_Edge.C+

*/


void  drawDelphes_Edge(TString plotsToMake="all") {

  initSamples("tt bj nm1 skimmed");
  setOutputDirectory("DelphesEdge"); 

  treatAllAsSM_=true; //don't treat signal as signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  //  stackSignal_=false;
  leg_x1 = 0.68;

  stackSignal_=true;
  setStackMode(true,false,false); //stack,norm,label override
  removeSample("naturalModel1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1!=1",kRed-9,"Misc Signal");
  addSample("naturalModel1:leptonsMatchChi2ToChi1==1",kRed,"True Edge Signal");

  //edge selection cuts
  TCut dileptons="mll>20";
  TCut sf = "isSF==1";
  TCut noZ = "mll<85|| mll>95";

  TCut dileptons_loose="mll_loose>20";
  TCut sf_loose = "isSF_loose==1";
  TCut noZ_loose = "mll_loose<85|| mll_loose>95";

  TCut jetsloose = "njets40eta3p0>=2";
  TCut jetstight = "njets40eta3p0>=3";
  TCut metloose = "MET>100";
  TCut mettight = "MET>150";


  //first the 8 TeV Edge, with Z veto added
  selection_ = dileptons && sf && jetsloose && mettight && noZ; //sel1

  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_sel1_mll",0,"GeV");

  //first the 8 TeV Edge, with Z veto added
  selection_ = dileptons && sf && jetstight && metloose && noZ; //sel2

  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_sel2_mll",0,"GeV");

  //then tighten
  selection_ = dileptons && sf && TCut("njets40eta3p0>=5") && TCut("HT>1000") && TCut("MET>200") && noZ && TCut("nbjets40tight>=2");
  var="mll"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new1_mll",0,"GeV");

  //try without b-tags
  selection_ = dileptons && sf && TCut("njets40eta3p0>=5") && TCut("HT>1000") && TCut("MET>200") && noZ;
  var="mll"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new1nob_mll",0,"GeV");

  //plot n-btag distribution
  selection_ = dileptons && sf && TCut("njets40eta3p0>=5") && TCut("HT>1000") && TCut("MET>200") && noZ;
  nbins=8; low=0; high=8;
  var="nbjets40tight"; xtitle="tight b tags";
  stackSignal_=false;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new1nob_nbT",0,"GeV");

  var="nbjets40loose"; xtitle="loose b tags";
  stackSignal_=false;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new1nob_nbL",0,"GeV");

  //add the Z pieces back in
  selection_ = dileptons && sf && TCut("njets40eta3p0>=5") && TCut("HT>1000") && TCut("MET>200");

  var="nbjets40tight"; xtitle="tight b tags";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new1nob-trueEdge_nbT",0,"GeV");

  var="nbjets40loose"; xtitle="loose b tags";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new1nob-trueEdge_nbL",0,"GeV");

  //
  nbins=40; low=0; high=800;
  var="MET"; xtitle="MET";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("004")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new1nob-trueEdge_MET",0,"GeV");

  nbins=40; low=0; high=3000;
  var="HT"; xtitle="HT";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("004")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new1nob-trueEdge_HT",0,"GeV");

  //new2
  selection_ = dileptons && sf && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1");
  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-} (GeV)";
  setPlotMaximum(600); //instead of Z veto and upper cut on mll
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2_mll",0,"GeV");

  selection_ = dileptons && sf && TCut("njets40eta3p0>=4") && TCut("MET>400") && TCut("nbjets40tight>=1") && TCut("mll_maxEta<1.4");
  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-} (GeV)";
  setPlotMaximum(600); //instead of Z veto and upper cut on mll
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2-central_mll",0,"GeV");

  removeSample("naturalModel1:leptonsMatchChi2ToChi1!=1");
  removeSample("naturalModel1:leptonsMatchChi2ToChi1==1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1",kRed-9,"Misc Signal");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1",kRed,"True Edge Signal");
  selection_ = dileptons_loose && sf_loose && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1");
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  setPlotMaximum(800); //instead of Z veto and upper cut on mll
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2_mllLoose",0,"GeV");

  //now the mll signal region
  TCut mll_loose_signal = "mll_loose>=20 && mll_loose<70";
  selection_ = mll_loose_signal && sf_loose && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1");
  var = "HT"; xtitle=var;
  nbins=30; low = 0; high = 3000;
  resetPlotMaximum();
  if (plotsToMake.Contains("all")||plotsToMake.Contains("006")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2mllLooseLowMass_HT",0,"GeV");

  removeSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1");
  removeSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1!=1",kRed-9,"Misc Signal");
  addSample("naturalModel1:leptonsMatchChi2ToChi1==1",kRed,"True Edge Signal");

  TCut mll_signal = "mll>=20 && mll<70";
  selection_ = mll_signal && sf && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1");
  var = "HT"; xtitle=var;
  resetPlotMaximum();
  if (plotsToMake.Contains("all")||plotsToMake.Contains("006")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2mllLowMass_HT",0,"GeV");

  selection_ = mll_signal && sf && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&& TCut("mll_maxEta<1.4");
  var = "HT"; xtitle=var;
  resetPlotMaximum();
  if (plotsToMake.Contains("all")||plotsToMake.Contains("007")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2mllLowMass-central_HT",0,"GeV");

  selection_ =dileptons && sf && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&& TCut("mll_maxEta<1.4")&&TCut("HT>1500");
  var = "mll"; xtitle=var;
  nbins=40; low=0; high=200;
  setPlotMaximum(200); //instead of Z veto and upper cut on mll
  if (plotsToMake.Contains("all")||plotsToMake.Contains("008")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500-central_mll",0,"GeV");

  removeSample("naturalModel1:leptonsMatchChi2ToChi1!=1");
  removeSample("naturalModel1:leptonsMatchChi2ToChi1==1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1",kRed-9,"Misc Signal");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1",kRed,"True Edge Signal");
  selection_ = dileptons_loose && sf_loose && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500"); //old NOMINAL SELECTION, used at pre-approval time
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  resetPlotMaximum();
  if (plotsToMake.Contains("all")||plotsToMake.Contains("009")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_mllLoose",0,"GeV");


  //now plot OF (but otherwise the same cuts)
  selection_ = dileptons_loose && !sf_loose && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500");
  resetPlotMaximum();
  if (plotsToMake.Contains("all")||plotsToMake.Contains("009")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_mllLoose_OF",0,"GeV");

  selection_ = dileptons_loose && sf_loose && TCut("njets40>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500"); //NOMINAL SELECTION
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  resetPlotMaximum();
  if (plotsToMake.Contains("all")||plotsToMake.Contains("009")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_mllLoose",0,"GeV");


  //now plot OF (but otherwise the same cuts)
  selection_ = dileptons_loose && !sf_loose && TCut("njets40>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500");
  resetPlotMaximum();
  if (plotsToMake.Contains("all")||plotsToMake.Contains("009")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_mllLoose_OF",0,"GeV");


  //repeat previous plots but with finer binning
  nbins=100; low=0; high=200;

  //for SF-OF
  TH1D* total_sf=0;
  TH1D* total_of = 0;
  TH1D* true_edge = 0;

  selection_ = dileptons_loose && sf_loose && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500");
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  setPlotMaximum(300); //instead of Z veto and upper cut on mll
  if (plotsToMake.Contains("all")||plotsToMake.Contains("010")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_mllLoose-finebins",0,"GeV");
    total_sf = (TH1D*)  totalsmsusy->Clone("total_sf");
    true_edge =  (TH1D*)getHist("naturalModel1:leptonsMatchChi2ToChi1_loose==1")->Clone("true_edge");
  }

  //now plot OF (but otherwise the same cuts)
  selection_ = dileptons_loose && !sf_loose && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500");
  resetPlotMaximum();
  if (plotsToMake.Contains("all")||plotsToMake.Contains("010")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_mllLoose_OF-finebins",0,"GeV");
    total_of = (TH1D*)  totalsmsusy->Clone("total_of");
  }

  if (total_sf!=0 && total_of!=0) { //do OF subtraction
    TH1D* sf_minus_of = (TH1D*)total_sf->Clone("sf_minus_of");
    sf_minus_of->Add(total_of,-1); //subtraction

    renewCanvas();
    true_edge->Draw("hist");

    double max = findOverallMax(true_edge);
    if ( findOverallMax(sf_minus_of)>max) max = findOverallMax(sf_minus_of);
    double min = findOverallMin(true_edge);
    if ( findOverallMin(sf_minus_of)<min) min = findOverallMin(sf_minus_of);
    true_edge->SetMinimum(min);
    true_edge->SetMaximum(max);
    sf_minus_of->SetMarkerStyle(4);
    sf_minus_of->Draw("same");

    TFile fout("DelphesEdge/sfMinusOf-fine.root","recreate");
    true_edge->Write();
    total_sf->Write();
    total_of->Write();
    sf_minus_of->Write();
    thecanvas->Write();
    fout.Close();
  }

  resetPlotMaximum();
  //with tight selection, want to look at
  // max eta
  // lepton isolation
  selection_ = dileptons_loose && sf_loose && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500");
  nbins=24; low=0; high=2.4;
  var="mll_maxEta_loose"; xtitle="maximum lepton eta";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("011")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_maxEtaLoose",0,"");

  selection_ = dileptons_loose && sf_loose && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500");
  nbins=30; low=-1.5; high=1.5;
  var="leptonIso1_loose"; xtitle="RelIso lepton 1";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("011")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_leptonIso1Loose",0,"");
  var="leptonIso2_loose"; xtitle="RelIso lepton 2";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("011")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_leptonIso2Loose",0,"");

  setStackMode(false,true,false); //stack,norm,label override
  var="leptonIso1_loose"; xtitle="RelIso lepton 1";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("011")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_leptonIso1Loose",0,"");
  var="leptonIso2_loose"; xtitle="RelIso lepton 2";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("011")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_leptonIso2Loose",0,"");

}

void print_zbi() {     //let's try to spit out Zbi for 10% systematic
  //assumes histos are already drawn, and assumes the name of the edge signal

    double ntotal = totalsm->Integral();//contains SUSY because we see treatAllAsSM_ flag
    double ntrueedge = getIntegral("naturalModel1:leptonsMatchChi2ToChi1_loose==1");
    double nbackground = ntotal - ntrueedge;
    double zbi = jmt::zbi( ntotal, nbackground, 0.1*nbackground);
    cout<<"Zbi = "<<zbi<<endl;


}

void drawDelphes_Edge_CutAndCount(TString plotsToMake="all") {

  initSamples("tt bj nm1 skimmed");
  setOutputDirectory("DelphesEdge"); 

  treatAllAsSM_=true; //don't treat signal as signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  stackSignal_=false;
  setStackMode(true,false,false); //stack,norm,label override
  removeSample("naturalModel1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1",kRed-9,"Misc Signal");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1",kRed,"True Edge Signal");

  //edge selection cuts
  TCut dileptons_loose="mll_loose>20";
  TCut sf_loose = "isSF_loose==1";
  TCut cutncount = "mll_loose<70";

  //NOMINAL (pre-appoval) CUT AND COUNT SELECTION
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_cutncount1",0,"");
    print_zbi();
  }

  //Same as above but using medium b-tags
 selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40medium>=1")&&TCut("HT>1500");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_cutncount1b",0,"");
    print_zbi();

  }

 //Same as above but using central jets
 selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=4") && TCut("MET>400")  && TCut("nbjets40medium>=1")&&TCut("HT>1500");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_cutncount1c",0,"");
    print_zbi();
  }


  //NOMINAL (pre-appoval) CUT AND COUNT SELECTION but with central jets
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_cutncount1d",0,"");
    print_zbi();
  }


  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500") &&TCut("mll_maxEta_loose<1.4");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002")) 
    {  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_cutncount2",0,""); print_zbi();}

  //plot HT 
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>800");
  nbins=48; low=800; high=2000;
  var="HT"; xtitle="HT";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_cutncount_test1",0,"");

  //plot MET
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=4") && TCut("MET>200")  && TCut("nbjets40tight>=1")&&TCut("HT>1400");
  nbins=30; low=200; high=800;
  var="MET"; xtitle="MET";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_cutncount_test2",0,"");

  //plot njets
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=2") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500");
  nbins=8; low=2; high=10;
  var="njets40"; xtitle="njets40";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("004")) 
  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_cutncount_test3",0,"");

  //nbjets after tighter njets
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>300")  &&TCut("HT>1200");
  nbins=8; low=0; high=8;
  var="nbjets40medium"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("008")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_cutncount_test6",0,"");

  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>300")  &&TCut("HT>1200");
  nbins=8; low=0; high=8;
  var="nbjets40tight"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("008")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_cutncount_test6b",0,"");

  //redo with tighter njets
  //plot HT 
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>300")  && TCut("nbjets40tight>=1")&&TCut("HT>800");
  nbins=48; low=800; high=2000;
  var="HT"; xtitle="HT";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005")) 
  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_cutncount_test4",0,"");

  //plot MET
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>200")  && TCut("nbjets40tight>=1")&&TCut("HT>1200");
  nbins=30; low=200; high=800;
  var="MET"; xtitle="MET";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005")) 
  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_cutncount_test5",0,"");

  // -- new cut and count region -- 
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>450")  && TCut("nbjets40tight>=1")&&TCut("HT>1250");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("006")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new3_cutncount1",0,"");
    print_zbi();} //NEW SAMPLES -- this is good! Zbi of 4!

  // try very tight MET -- also promising
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>650")  && TCut("nbjets40tight>=1")&&TCut("HT>1200");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("008")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new5_cutncount1",0,"");
    print_zbi();} 
  
  // try very tight btags
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>300")  && TCut("nbjets40medium>=3")&&TCut("HT>1200");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("009")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new6_cutncount1",0,"");
    print_zbi();} 

  if (plotsToMake.Contains("cutflowtest")) {
    savePlots_=false;
    for (int htcut = 1000; htcut<2250; htcut+=250) {
      for (int metcut = 250; metcut<750; metcut+=50) {
	TString metcutstring,htcutstring;
	metcutstring.Form("MET> %d",metcut);
	htcutstring.Form("HT> %d",htcut);
	selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut(metcutstring.Data())  && TCut("nbjets40medium>=2")&&TCut(htcutstring.Data());
	drawPlots(var,nbins,low,high,xtitle,"Events", "blah",0,"");
	cout<<"Cutflow "<<htcut<<" "<<metcut<<" "; print_zbi();
      }
    }
    savePlots_=true;
  }

  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>450")  && TCut("nbjets40tight>=1")&&TCut("HT>1250")&&TCut("mll_maxEta_loose<1.4");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("006")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new3_cutncount1-central",0,"");
    print_zbi();}

  // -- new cut and count region 2 --
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=7") && TCut("MET>350")  && TCut("nbjets40tight>=1")&&TCut("HT>1200");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("007")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new4_cutncount1",0,"");
    print_zbi();
  }

}

void ewkino(TString plotsToMake="all") {

  initSamples("skimmed tt bj nm1");
  setOutputDirectory("DelphesEdgeEwkino");

  treatAllAsSM_=true; //don't treat signal as signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  stackSignal_=true;
  setStackMode(true,false,false); //stack,norm,label override
  removeSample("naturalModel1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1",kRed-9,"Misc Signal");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1",kRed,"True Edge Signal");

  //edge selection cuts
  TCut dileptons_loose="mll_loose>10";
  TCut sf_loose = "isSF_loose==1";
  TCut lowmass = "mll_loose<80"; //gotta get rid of the Z

  TCut bveto = "nbjets40loose==0";
  //  TCut tau="nTaus>0";
  TCut l3="MT_l3MET_loose>=0";
  //  TCut ht="HT<500";
  TCut loosemet = "MET>50";

  //basic trilepton selection 
  selection_ =dileptons_loose && sf_loose && lowmass && l3 &&loosemet;

  nbins=60; low=0; high=120;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l_mll",0,"GeV");

  //plot nb
  setStackMode(false,true,false); //stack,norm,label override
  nbins=5; low=0; high=5;
  var="nbjets40loose"; xtitle="nb loose 40";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l_nb40loose",0,"");

  //plot met
  setStackMode(false,true,false); //stack,norm,label override
  nbins=50; low=0; high=500;
  var="MET"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l_met",0,"GeV");

  //plot mt
  setStackMode(false,true,false); //stack,norm,label override
  nbins=50; low=0; high=500;
  var="MT_l3MET_loose"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l_mtl3",0,"GeV");

  // lessons:
  //b veto, high met
  selection_ =dileptons_loose && sf_loose && lowmass && l3 && bveto && TCut("MET>200") && TCut("MT_l3MET_loose>150");
  setStackMode(true,false,false); //stack,norm,label override
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-2_mll",0,"GeV");

  //let in more background for some studies
  selection_ =dileptons_loose && sf_loose && lowmass && l3 && bveto && TCut("MET>150") && TCut("MT_l3MET_loose>130");
  setStackMode(false,true,false); //stack,norm,label override

  nbins=50; low=0; high=1000;
  var="HT"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-3_ht",0,"GeV");

  nbins=50; low=0; high=300;
  var="leptonPt1_loose"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-3_leppt1",0,"GeV");

  nbins=50; low=0; high=200;
  var="leptonPt2_loose"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-3_leppt2",0,"GeV");

  nbins=8; low=0; high=8;
  var="njets40eta3p0"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-3_njets",0,"GeV");

  //try complete jet veto
  selection_ =dileptons_loose && sf_loose && lowmass && l3 && TCut("njets40eta3p0==0") && TCut("MET>150") && TCut("MT_l3MET_loose>140");
  setStackMode(true,false,false); //stack,norm,label override
  nbins=60; low=0; high=120;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("004")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-4_mll",0,"GeV");

  setStackMode(false,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-4_mll",0,"GeV");

  setStackMode(true,false,false); //stack,norm,label override
  removeSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1");
  removeSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1");
  addSample("naturalModel1:isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&(Chi2ToChi1Code==1||Chi2ToChi1Code==2)",kRed-9,"SF (edge misreco)");
  addSample("naturalModel1:isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&Chi2ToChi1Code>=10",kMagenta,"SF (N2 #rightarrow N1 with tau)");
  addSample("naturalModel1:isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&(!(Chi2ToChi1Code==1||Chi2ToChi1Code==2||Chi2ToChi1Code>=10))",kTeal,"SF (no edge)");
  addSample("naturalModel1:isSF_loose==1&&leptonsMatchChi2ToChi1_loose==1",kRed,"SF (edge reco ok)");
  if (plotsToMake.Contains("all")||plotsToMake.Contains("006")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-4-split_mll",0,"GeV");

  nbins=7; low=0; high=7;
  var="nTrueElMu"; xtitle=var;
 if (plotsToMake.Contains("all")||plotsToMake.Contains("006")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-4-split_ntrueElMu",0,"");


 //another tack. no jet cuts except b veto, but tight in MT and MET
 // lessons:
  //b veto, high met
  selection_ =dileptons_loose && sf_loose && lowmass && l3 && bveto && TCut("MET>250") && TCut("MT_l3MET_loose>200");
  setStackMode(true,false,false); //stack,norm,label override
  nbins=60; low=0; high=120;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("007")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-5_mll",0,"GeV");

}


void ewkino_signalonly(TString plotsToMake="all") {

  TString signalName="naturalModel1";
  initSamples(signalName);
  setOutputDirectory("DelphesSignal");

  treatAllAsSM_=true; //don't treat signal as signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  stackSignal_=true;
  setStackMode(true,false,false); //stack,norm,label override
  /*
  removeSample("naturalModel1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1",kRed-9,"Misc Signal");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1",kRed,"True Edge Signal");
  */

  //edge selection cuts
  TCut dileptons_loose="mll_loose>0";
  TCut sf_loose = "isSF_loose==1";
  TCut lowmass = "mll_loose<200";

  TCut bveto = "nbjets40loose==0";
  TCut l3="MT_l3MET_loose>=0";
  TCut loosemet = "MET>50";
  selection_ =dileptons_loose && sf_loose && lowmass && l3 && TCut("njets40eta3p0==0") && TCut("MET>150") && TCut("MT_l3MET_loose>140");

  clearSamples();
  addSample(signalName+":SusyProductionMode==2000000",kYellow+1,"slepton");
  addSample(signalName+":SusyProductionMode==200000",kOrange+1,"EWKino");
  addSample(signalName+":SusyProductionMode==2000",kBlue,"stop");
  addSample(signalName+":SusyProductionMode==200",kGreen,"sbottom");
  addSample(signalName+":SusyProductionMode==20000",kMagenta,"gluino");
  addSample(signalName+":SusyProductionMode==10010||SusyProductionMode==20",kRed,"other");
  setStackMode(true,false,false); //stack,norm,label override
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-4_prod_mll",0,"GeV");
  setStackMode(false,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-4_prod_mll",0,"GeV");


  clearSamples();
  addSample(signalName+":isSF_loose==0",kRed,"OF");
  addSample(signalName+":isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&(Chi2ToChi1Code==1||Chi2ToChi1Code==2)",kOrange,"SF (edge misreco)");
  addSample(signalName+":isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&Chi2ToChi1Code>=10",kMagenta,"SF (N2 #rightarrow N1 with tau)");
  addSample(signalName+":isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&(!(Chi2ToChi1Code==1||Chi2ToChi1Code==2||Chi2ToChi1Code>=10))",kTeal,"SF (no edge)");
  addSample(signalName+":isSF_loose==1&&leptonsMatchChi2ToChi1_loose==1",kBlue,"SF (edge reco ok)");
  setStackMode(false,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-4_truth_mll",0,"GeV");

  //no OF for stack
  clearSamples();
  addSample(signalName+":isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&(Chi2ToChi1Code==1||Chi2ToChi1Code==2)",kOrange,"SF (edge misreco)");
  addSample(signalName+":isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&Chi2ToChi1Code>=10",kMagenta,"SF (N2 #rightarrow N1 with tau)");
  addSample(signalName+":isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&(!(Chi2ToChi1Code==1||Chi2ToChi1Code==2||Chi2ToChi1Code>=10))",kTeal,"SF (no edge)");
  addSample(signalName+":isSF_loose==1&&leptonsMatchChi2ToChi1_loose==1",kBlue,"SF (edge reco ok)");
  setStackMode(true,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-4_truth_mll",0,"GeV");

}

void preselection(TString plotsToMake="all") {
  setOutputDirectory("DelphesEdge");

  initSamples("combinesm");
  setOutputDirectory("DelphesEdge");

  treatAllAsSM_=true; //don't treat signal as signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  //  stackSignal_=false;

  setStackMode(false,true,false); //stack,norm,label override
  removeSample("naturalModel1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1!=1",kRed-9,"Misc Signal");
  addSample("naturalModel1:leptonsMatchChi2ToChi1==1",kRed,"True Edge Signal");

  //edge selection cuts
  TCut dileptons="mll>20";
  TCut dileptons_loose="mll_loose>20";
  TCut sf = "isSF==1";
  TCut sf_loose = "isSF_loose==1";
  TCut jetsloose = "njets40eta3p0>=2";
  TCut jetstight = "njets40eta3p0>=3";
  TCut metloose = "MET>100";
  TCut mettight = "MET>150";
  //  TCut noZ = "mll<85|| mll>95";

  //use 8 TeV edge as a preselection

  selection_ = dileptons && sf && jetstight && metloose ;

  nbins=35; low=100; high=800;
  var="MET"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_presel_met",0,"GeV");
 
  nbins=100; low=0; high=3000;
  var="HT"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_presel_ht",0,"GeV");

  nbins=12; low=0; high=12;
  var="njets40eta3p0"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_presel_njets40",0,"GeV");

  nbins=6; low=0; high=6;
  var="nbjets40loose"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_presel_nbjets40loose",0,"GeV");

  nbins=6; low=0; high=6;
  var="nbjets40tight"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_presel_nbjets40tight",0,"GeV");

  nbins=40; low=0; high=200;
  var="mll"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_presel_mll",0,"GeV");

  selection_ = dileptons_loose && sf_loose && jetstight && metloose ;
  nbins=40; low=0; high=200;
  var="mll_loose"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_presel_mll_loose",0,"GeV");


}


void compareOFSF(TString option="combinesm") { //
  //compare SF/OF
  setOutputDirectory("DelphesEdge");

  if (option.Contains("combinesm")) {
    initSamples("combinesm skimmed");
    
    clearSamples();
    addSample("tt-4p-0-600-v1510_14TEV:isSF_loose==1",kRed,"SF (tt)");
    addSample("tt-4p-0-600-v1510_14TEV:isSF_loose==0",kBlue,"OF (tt)");

  }
  else if (option.Contains("all")) {
    initSamples("combinesm skimmed");
    chainSamples("tt-4p-0-600-v1510_14TEV","naturalModel1");//combine SM and SUSY into one sample
    clearSamples();
    addSample("tt-4p-0-600-v1510_14TEV:isSF_loose==1",kRed ,"SF");
    addSample("tt-4p-0-600-v1510_14TEV:isSF_loose==0",kBlue,"OF");
  }
  else if (option.Contains("signal")) {
    initSamples("signal");
    
    clearSamples();
    addSample("naturalModel1:isSF_loose==1&&leptonsMatchChi2ToChi1_loose!=1",kRed,"SF (signal; no edge)");
    addSample("naturalModel1:isSF_loose==0",kBlue,"OF (signal)");
  }
  //  else if (option=="total") {
  //initSamples();
    
  //clearSamples();
    //my tricks don't work here
    //    addSample("naturalModel1:isSF==1",kRed,"SF (signal)");
    //    addSample("naturalModel1:isSF==0",kBlue,"OF (signal)");
  //}
  else assert(0);

  doRatio_=true;  
  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);

  stackSignal_=false;

  setStackMode(false,false,false); //stack,norm,label override

  //edge selection cuts
  TCut dileptons="mll_loose>20";
  TCut btag="nbjets40tight>=1";
  TCut rejectedge="leptonsMatchChi2ToChi1_loose!=1";
  //tightened selection without b-tags (remove SF cut)
  selection_ = dileptons && TCut("njets40eta3p0>=4") && TCut("HT>1500") && TCut("MET>400") &&rejectedge&&btag;
  nbins=60; low=20; high=200;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  doRatio_=true;  
  if (option.Contains("001")) drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_mllloose_OFSFr_"+jmt::fortranize(option),0,"GeV");
  doRatio_=false;  
  if (option.Contains("001")) drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_mllloose_OFSF_"+jmt::fortranize(option),0,"GeV");

  //plot only dilepton ttbar events
  selection_ = dileptons && TCut("njets40eta3p0>=4") && TCut("HT>1500") && TCut("MET>400") &&TCut("ttbarDecayCode==1")&&btag;
  nbins=60; low=20; high=200;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  doRatio_=true;  ratioMin=0.5; ratioMax = 1.5;
  if (option.Contains("002")) drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500-genTt2l_mllloose_OFSFr_"+jmt::fortranize(option),0,"GeV");

}


void drawSignalByProductionMode() {
 lumiScale_=3000e3;

  initSamples("nm1 skimmed");
  setOutputDirectory("DelphesEdge");
  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(false);
  doRatio_=false;  
  treatAllAsSM_=true; //don't treat signal as signal

  setStackMode(true,false,false); //stack,norm,label override

  //edge selection cuts
  TCut dileptons="mll_loose>20";
  TCut sf = "isSF_loose==1";
  TCut jets = "njets40eta3p0>=4";
  TCut met = "MET>400";
  TCut ht = "HT>1500";  

  TCut btags = "nbjets40tight>=1";
  TCut bveto = "nbjets40tight==0";

  // does not include sf cut
  TCut analysis_selection = dileptons && jets && met && ht && btags && sf;
  nbins = 40; low = 20; high=100;
  var = "mll_loose"; xtitle="mll";

  clearSamples();
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1",kRed-9,"Non-edge SUSY");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1&&(SusyProductionMode==2000000||SusyProductionMode==200000)",kYellow+1,"l~+EWKino edge");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1&&(SusyProductionMode==2000||SusyProductionMode==200)",kBlue,"t~t~+b~b~ edge");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1&&SusyProductionMode==20000",kMagenta,"g~g~ edge");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1&&SusyProductionMode==10010",kCyan,"g~q~ edge");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1&&SusyProductionMode==20",kGreen,"q~q~ edge");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1&&!(SusyProductionMode==2000000||SusyProductionMode==200000||SusyProductionMode==2000||SusyProductionMode==200||SusyProductionMode==20000||SusyProductionMode==10010||SusyProductionMode==20)",kRed,"other edge");
  

  selection_ = analysis_selection;
  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_signal_by_prod",0,"GeV");
 
}

void generateHistosForFit() {
  lumiScale_=3000e3;

  initSamples("tt bj nm1 skimmed");
  setOutputDirectory("DelphesEdge");
  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  


  setStackMode(true,false,false); //stack,norm,label override
  stackSignal_=false; //the logic of this is a bit broken at the moment, so just don't stack it

  //edge selection cuts
  TCut dileptons="mll_loose>20";
  TCut sf = "isSF_loose==1";
  TCut jets = "njets40eta3p0>=4";
  TCut met = "MET>400";
  TCut ht = "HT>1500";  

  TCut btags = "nbjets40tight>=1";
  TCut bveto = "nbjets40tight==0";

  TCut ee = "abs(leptonFlavor1_loose)==11";
  TCut mm = "abs(leptonFlavor1_loose)==13";

  // does not include sf cut
  TCut analysis_selection = dileptons && jets && met && ht && btags;
  //could split by e and mu
  TCut dy_selection = dileptons && sf && jets && ht && TCut("MET<40") && bveto; 

  nbins = 280; low = 20; high=300;
  var = "mll_loose"; xtitle="mll";

  selection_ = analysis_selection && sf && ee;
  drawPlots(var,nbins,low,high,xtitle,"Events", "templates_new2HT1500_mll_ee",0,"GeV");
  TH1D* smsusy_mll_ee = (TH1D*)  totalsmsusy->Clone("smsusy_mll_ee");
  TH1D* sm_mll_ee = (TH1D*)  totalsm->Clone("sm_mll_ee");

  selection_ = analysis_selection && sf && mm;
  drawPlots(var,nbins,low,high,xtitle,"Events", "templates_new2HT1500_mll_mm",0,"GeV");
  TH1D* smsusy_mll_mm = (TH1D*)  totalsmsusy->Clone("smsusy_mll_mm");
  TH1D* sm_mll_mm = (TH1D*)  totalsm->Clone("sm_mll_mm");

  selection_ = analysis_selection && !sf;
  drawPlots(var,nbins,low,high,xtitle,"Events", "templates_new2HT1500_mll_OF",0,"GeV");
  TH1D* smsusy_mll_OF = (TH1D*)  totalsmsusy->Clone("smsusy_mll_OF");
  TH1D* sm_mll_OF = (TH1D*)  totalsm->Clone("sm_mll_OF");

  selection_ = dy_selection;
  nbins = 280*5; low = 20; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "templates_new2HT1500_mll_DY",0,"GeV");
  TH1D* smsusy_mll_DY = (TH1D*)  totalsmsusy->Clone("smsusy_mll_DY");
  TH1D* sm_mll_DY = (TH1D*)  totalsm->Clone("sm_mll_DY");

  // new histograms: include sm and non-edge susy
  //i.e. the "background only" hypothesis when fitting for edge signal
  removeSample("naturalModel1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1",kRed-9,"Misc Signal");

  nbins = 280; low = 20; high=300;
  selection_ = analysis_selection && sf && ee;
  drawPlots(var,nbins,low,high,xtitle,"Events", "templates_new2HT1500_mll_ee_noedge",0,"GeV");
  TH1D* smsusynoedge_mll_ee = (TH1D*)  totalsmsusy->Clone("smsusynoedge_mll_ee");

  selection_ = analysis_selection && sf && mm;
  drawPlots(var,nbins,low,high,xtitle,"Events", "templates_new2HT1500_mll_mm_noedge",0,"GeV");
  TH1D* smsusynoedge_mll_mm = (TH1D*)  totalsmsusy->Clone("smsusynoedge_mll_mm");

  selection_ = analysis_selection && !sf ;
  drawPlots(var,nbins,low,high,xtitle,"Events", "templates_new2HT1500_mll_OF_noedge",0,"GeV");
  TH1D* smsusynoedge_mll_OF = (TH1D*)  totalsmsusy->Clone("smsusynoedge_mll_OF");
 
  TFile fout("DelphesEdge/templates.root","recreate");
  smsusy_mll_ee->Write();
  smsusy_mll_mm->Write();
  smsusy_mll_OF->Write();
  smsusy_mll_DY->Write();
  smsusynoedge_mll_ee->Write();
  smsusynoedge_mll_mm->Write();
  smsusynoedge_mll_OF->Write();
  sm_mll_ee->Write();
  sm_mll_mm->Write();
  sm_mll_OF->Write();
  sm_mll_DY->Write();
  fout.Close();

}
