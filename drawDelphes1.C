#include "drawDelphesBase.C"
/*
--- to compile ---

.L TSelectorMultiDraw.C+
.L CrossSectionTable.cxx+
.L ConfigurationDescriptions.cxx+

.L drawDelphes1.C+

*/

void  drawDelphes1(TString plotsToMake="all") {


  initSamples();
  setOutputDirectory("Delphes1");

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  //normalize to unit area
  setStackMode(false,true,false); //stack,norm,label override
  stackSignal_=false;

  selection_ = "(1)"; //no cuts!
  nbins=40; low=0; high=1000;
  var="MET"; xtitle="MET (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001"))  drawPlots(var,nbins,low,high,xtitle,"Events", "MET",0,"GeV");

  selection_ = "(1)"; //no cuts!
  nbins=50; low=0; high=3000;
  var="HT"; xtitle="HT (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002"))   drawPlots(var,nbins,low,high,xtitle,"Events", "HT",0,"GeV");

  selection_ = "(1)"; //no cuts!
  nbins=14; low=0; high=14;
  var="njets40"; xtitle="njets (40 GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003"))  drawPlots(var,nbins,low,high,xtitle,"Events", "njets40",0,"");

  selection_ = "(1)"; //no cuts!
  nbins=6; low=0; high=6;
  var="nbjets40tight"; xtitle="n tight b jets (40 GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("004"))   drawPlots(var,nbins,low,high,xtitle,"Events", "nTbjets40",0,"");

  selection_ = "(1)"; //no cuts!
  nbins=50; low=0; high=2000;
  var="MT2"; xtitle="MT2 (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005"))   drawPlots(var,nbins,low,high,xtitle,"Events", "MT2",0,"GeV");

  //edge selection cuts
  TCut dileptons="mll>20";
  TCut sf = "isSF==1";
  TCut jetsloose = "njets30eta3p0>=2";
  TCut jetstight = "njets30eta3p0>=3";
  TCut metloose = "MET>100";
  TCut mettight = "MET>150";

  selection_ = dileptons&&sf&&jetsloose&&mettight; //Edge selection 1
  nbins=50; low=0; high=300;
  var="mll"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("006"))   drawPlots(var,nbins,low,high,xtitle,"Events", "edge_sel1_mll",0,"GeV");

  //same plot, zoom on low mass region
  nbins=100; low=0; high=100;
  var="mll"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("012"))   drawPlots(var,nbins,low,high,xtitle,"Events", "edge_sel1_mll_zoom",0,"GeV");

  if (plotsToMake=="013") {//not compatible with other plots
    clearSamples();
    addSample("susyhit_slhaScenario1_v02",kRed,"Scenario 1");
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_sel1_mll_zoom_sigonly",0,"GeV");

    nbins=100; low=0; high=500;
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_sel1_mll_wide_sigonly",0,"GeV");

    selection_ = dileptons&&!sf&&jetsloose&&mettight; //sanity check -- veto SF
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_sel1_mllOF_zoom_sigonly",0,"GeV");


  }

  //also show with 'real' normalization
  setStackMode(true,false,false); //stack,norm,label override
  stackSignal_=false;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("007"))   drawPlots(var,nbins,low,high,xtitle,"Events", "edge_sel1_mll",0,"GeV");

  setStackMode(false,true,false); //stack,norm,label override
  selection_ = dileptons&&sf&&jetstight&&metloose; //Edge selection 2
  nbins=50; low=0; high=300;
  var="mll"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("008"))   drawPlots(var,nbins,low,high,xtitle,"Events", "edge_sel2_mll",0,"GeV");

  //also show with 'real' normalization
  setStackMode(true,false,false); //stack,norm,label override
  stackSignal_=false;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("009"))   drawPlots(var,nbins,low,high,xtitle,"Events", "edge_sel2_mll",0,"GeV");

  //now try a hadronic analysis selection
  TCut noleptons = "nElectrons==0 && nMuons==0";
  TCut ht600 = "HT>=600";
  TCut jets="njets40>=3";
  selection_=noleptons && ht600 && metloose && jets;

  nbins=40; low=0; high=800;
  var="MT2"; xtitle="MT2 (GeV)";
  setStackMode(false,true,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("010"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_presel_MT2",0,"GeV");
  setStackMode(true,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("011"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_presel_MT2",0,"GeV");

  //try an MT2 cut instead?
  TCut mt2baseline = "MT2>=200";
  selection_=noleptons&&ht600&&jets&&mt2baseline;
  nbins=40; low=0; high=2000;
  var="HT"; xtitle="HT (GeV)";
  setStackMode(false,true,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("014"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2sel_HT",0,"GeV");
  setStackMode(true,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("014"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2sel_HT",0,"GeV");


  //cut much tighter on HT and plot MT2
  TCut tighterht = "HT>1500";
  selection_ = noleptons && tighterht && jets && mt2baseline;
  //plot MT2 and njets
  nbins=40; low=0; high=2000;
  var="MT2"; xtitle="MT2 (GeV)";
  setStackMode(false,true,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("015"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2sel2_MT2",0,"GeV");
  setStackMode(true,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("015"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2sel2_MT2",0,"GeV");

  nbins=12; low=0; high=12;
  var="njets40"; xtitle="njets40 (GeV)";
  setStackMode(false,true,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("015"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2sel2_njets40",0,"GeV");
  setStackMode(true,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("015"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2sel2_njets40",0,"GeV");

  //tighter further!
  TCut tightermt2 = "MT2>350";
  TCut tighterjets = "njets40>=8";
  selection_ = noleptons&&TCut("HT>1000")&&tightermt2&&tighterjets;
  //plot HT again
  nbins=40; low=1000; high=3000;
  var="HT"; xtitle="HT (GeV)";

  setStackMode(true,false,false); //stack,norm,label override
  customDivisionsVal_ = 505;
  customDivisions_=true;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("016"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2sel3_HT",0,"GeV");
  customDivisions_=false;
  //ok, now we're getting somewhere!

  TCut ht2000="HT>2000";
  selection_ = noleptons && ht2000 && tightermt2 && tighterjets;
  //nbjets
  nbins=9; low=0; high=9;
  var="nbjets40tight"; xtitle="n tight b tags";
  setStackMode(true,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("017"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2sel4_nbjets",0,"");
  //ntaus
  nbins=3; low=0; high=3;
  var="nTaus"; xtitle="n tau tags";
  setStackMode(true,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("017"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2sel4_ntaus",0,"");
  setStackMode(false,true,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("017"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2sel4_ntaus",0,"");
  //and MT2
  setStackMode(true,false,false); //stack,norm,label override
  selection_ = noleptons && ht2000 && TCut("MT2>250") && tighterjets;
  nbins=65; low=250; high=1000;
  var="MT2"; xtitle="MT2 (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("017"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2sel4_mt2",0,"GeV");

  //try a 2d plot
  selection_ = noleptons && TCut("HT>800") && TCut("MT2>200") && tighterjets;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("018")) 
    draw2d("MT2",15,200,500,"HT",15,800,3000,"MT2","HT", "hadronic_mt2sel5_mt2Vht",0,0,"susyhit_slhaScenario1_v02");

}
