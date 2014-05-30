#include "drawDelphesBase.C"
/*
--- to compile ---

.L TSelectorMultiDraw.C+
.L CrossSectionTable.cxx+
.L ConfigurationDescriptions.cxx+

.L drawDelphes_signalOnly.C+

*/

void  drawDelphes_signalOnly(TString plotsToMake="all") {

  initSamples("signal");
  setOutputDirectory("DelphesSignal");

  treatAllAsSM_=true; //don't treat signal as signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  stackSignal_=false;

  setStackMode(false,false,false); //stack,norm,label override
  selection_ = "(1)"; //no cuts!
  nbins=8; low=1; high=9;
  //slepton, ewkino, stop, sbottom, gluino, gluino+squark, sq/sq
  var="(SusyProductionMode==2000000)*1+(SusyProductionMode==200000)*2+(SusyProductionMode==2000)*3+(SusyProductionMode==200)*4+(SusyProductionMode==20000)*5+(SusyProductionMode==10010)*6+(SusyProductionMode==20)*7"; xtitle="Susy Production Mode";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_prod",0,"");

  setStackMode(false,true,false); //stack,norm,label override
  //split by SUSY mode
  clearSamples();
  addSample("susyhit_Scenario1_v02:SusyProductionMode==2000000",kYellow+1,"slepton");
  addSample("susyhit_Scenario1_v02:SusyProductionMode==200000",kOrange+1,"EWKino");
  addSample("susyhit_Scenario1_v02:SusyProductionMode==2000",kBlue,"stop");
  addSample("susyhit_Scenario1_v02:SusyProductionMode==200",kGreen,"sbottom");
  addSample("susyhit_Scenario1_v02:SusyProductionMode==20000",kMagenta,"gluino");
  addSample("susyhit_Scenario1_v02:SusyProductionMode==10010||SusyProductionMode==20",kRed,"other");

  selection_ = "(1)"; //no cuts!
  nbins=50; low=0; high=800;
  var="MET"; xtitle="MET";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_byMode_MET",0,"");

  nbins=50; low=0; high=3000;
  var="HT"; xtitle="HT";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_byMode_HT",0,"");

  nbins=12; low=0; high=12;
  var="njets30eta3p0"; xtitle="njets30eta3p0";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_byMode_njets30eta3p0",0,"");

  nbins=6; low=0; high=6;
  var="nbjets40tight"; xtitle="nbjets40tight";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_byMode_nbjets40tight",0,"");

  selection_ = "mll>0 && isSF==1"; //loosest possible dilepton selection
  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-}";

  setStackMode(true,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_byMode_mll",0,"");
  setStackMode(false,true,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_byMode_mll",0,"");


  //split instead by Chi2ToChi1Code
 setStackMode(false,true,false); //stack,norm,label override
  //split by SUSY mode
  clearSamples();
  addSample("susyhit_Scenario1_v02:Chi2ToChi1Code==0",kYellow+1,"No N2 #rightarrow N1");
  addSample("susyhit_Scenario1_v02:Chi2ToChi1Code==1",kRed,"N2 #rightarrow N1 via e/#mu");
  addSample("susyhit_Scenario1_v02:Chi2ToChi1Code==2",kOrange,"N2 #rightarrow N1 via e/#mu x 2");
  addSample("susyhit_Scenario1_v02:Chi2ToChi1Code>=10",kBlue,"N2 #rightarrow N1 via #tau");

  selection_ = "(1)"; //no cuts!
  nbins=8; low=1; high=9;
  //slepton, ewkino, stop, sbottom, gluino, gluino+squark, sq/sq
  var="(SusyProductionMode==2000000)*1+(SusyProductionMode==200000)*2+(SusyProductionMode==2000)*3+(SusyProductionMode==200)*4+(SusyProductionMode==20000)*5+(SusyProductionMode==10010)*6+(SusyProductionMode==20)*7"; xtitle="Susy Production Mode";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_prod_byN2toN1",0,"");

  setStackMode(false,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_prod_byN2toN1",0,"");

  selection_ = "mll>0 && isSF==1"; //loosest possible dilepton selection
  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-}";

  setStackMode(false,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_byN2toN1_mll",0,"");
  setStackMode(false,true,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_byN2toN1_mll",0,"");

  //now, forget gen-level (or not completely)
  clearSamples();
  addSample("susyhit_Scenario1_v02:isSF==0",kRed,"OF");
  addSample("susyhit_Scenario1_v02:isSF==1&&leptonsMatchChi2ToChi1==0&&(Chi2ToChi1Code==1||Chi2ToChi1Code==2)",kOrange,"SF (edge misreco)");
  addSample("susyhit_Scenario1_v02:isSF==1&&leptonsMatchChi2ToChi1==0&&Chi2ToChi1Code>=10",kMagenta,"SF (N2 #rightarrow N1 with tau)");
  addSample("susyhit_Scenario1_v02:isSF==1&&leptonsMatchChi2ToChi1==0&&(!(Chi2ToChi1Code==1||Chi2ToChi1Code==2||Chi2ToChi1Code>=10))",kTeal,"SF (no edge)");
  //  addSample("susyhit_Scenario1_v02:isSF==1&&leptonsMatchChi2ToChi1==0&&nZFromSusy==0",kTeal,"SF (no reco match; no Z)");
  //  addSample("susyhit_Scenario1_v02:isSF==1&&leptonsMatchChi2ToChi1==0&&nZFromSusy>0",kOrange,"SF (no reco match; with Z)");
  addSample("susyhit_Scenario1_v02:isSF==1&&leptonsMatchChi2ToChi1==1",kBlue,"SF (edge reco ok)");

  //leptonsMatchChi2ToChi1

  setStackMode(false,false,false); //stack,norm,label override 
  selection_ = "mll>0"; //loosest possible dilepton selection
  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-}";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("004"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_bySF_mll",0,"");


  clearSamples(); //now just plot SF part for stack
  addSample("susyhit_Scenario1_v02:isSF==1&&leptonsMatchChi2ToChi1==0&&(Chi2ToChi1Code==1||Chi2ToChi1Code==2)",kOrange,"SF (edge misreco)");
  addSample("susyhit_Scenario1_v02:isSF==1&&leptonsMatchChi2ToChi1==0&&Chi2ToChi1Code>=10",kMagenta,"SF (N2 #rightarrow N1 with tau)");
  addSample("susyhit_Scenario1_v02:isSF==1&&leptonsMatchChi2ToChi1==0&&(!(Chi2ToChi1Code==1||Chi2ToChi1Code==2||Chi2ToChi1Code>=10))",kTeal,"SF (no edge)");
  //  addSample("susyhit_Scenario1_v02:isSF==1&&leptonsMatchChi2ToChi1==0&&nZFromSusy==0",kTeal,"SF (no reco match; no Z)");
  //  addSample("susyhit_Scenario1_v02:isSF==1&&leptonsMatchChi2ToChi1==0&&nZFromSusy>0",kOrange,"SF (no reco match; with Z)");
  addSample("susyhit_Scenario1_v02:isSF==1&&leptonsMatchChi2ToChi1==1",kBlue,"SF (edge reco ok)");
  setStackMode(true,false,false); //stack,norm,label override 
  if (plotsToMake.Contains("all")||plotsToMake.Contains("004"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_bySF_mll",0,"");

  //finally, let's go for the actual cut of the Edge analysis
  clearSamples();
  addSample("susyhit_Scenario1_v02:mll_maxEta<1.4",kRed,"Scenario 1 (central)");
  addSample("susyhit_Scenario1_v02:mll_maxEta>1.4",kBlue,"Scenario 1 (forward)");

  setStackMode(false,false,false); //stack,norm,label override 
  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-}";

  selection_="mll>20&&isSF==1 && njets40eta3p0>=2 && MET>150"; //one search region
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_sel1_mll",0,"");

  selection_="mll>20&&isSF==1 && njets40eta3p0>=3 && MET>100"; //the other search region
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_sel2_mll",0,"");

  //do it again, but split by SUSY production mode and stack ; integrate over maxEta
  clearSamples();
  addSample("susyhit_Scenario1_v02:SusyProductionMode==2000000",kYellow+1,"slepton");
  addSample("susyhit_Scenario1_v02:SusyProductionMode==200000",kOrange+1,"EWKino");
  addSample("susyhit_Scenario1_v02:SusyProductionMode==2000",kBlue,"stop");
  addSample("susyhit_Scenario1_v02:SusyProductionMode==200",kGreen,"sbottom");
  addSample("susyhit_Scenario1_v02:SusyProductionMode==20000",kMagenta,"gluino");
  addSample("susyhit_Scenario1_v02:SusyProductionMode==10010||SusyProductionMode==20",kRed,"other");

  setStackMode(true,false,false); //stack,norm,label override 
  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-}";

  selection_="mll>20&&isSF==1 && njets40eta3p0>=2 && MET>150"; //one search region
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_sel1_byMode_mll",0,"");

  selection_="mll>20&&isSF==1 && njets40eta3p0>=3 && MET>100"; //the other search region
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005"))  drawPlots(var,nbins,low,high,xtitle,"Events", "Sc1_sel2_byMode_mll",0,"");

}
