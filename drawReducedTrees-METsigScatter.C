void go() {

  /*
                                                                                                                                                                 
.L TSelectorMultiDraw.C+          
.L CrossSectionTable.cxx+        
 .L ConfigurationDescriptions.cxx+                                                                                                                                                                                    
.L drawReducedTrees_hbbhbb.C+      
                  
  */

  initHiggsSamples69(false,"qcd ttbar znunu");

  //just want a 2d plot of data for MET:METsig

  useTrigEff_=false;

  //define useful cuts                                                                                                                               
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40";
  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1"; //||passMC_DiCentralPFJet30_PFMHT80==1
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;
  
  TCut ra2btrigger = "passMC_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80==1 || passMC_DiCentralPFJet50_PFMET80==1";
  
  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoTracks15_005_03==0";//&&nIsoTracks5_005_03<2";                                                 
  
  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0 &&nIsoTracks15_005_03_lepcleaned==0";
  
  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut jet2high="jetpt2>70";
  TCut jets=(njets4||njets5) &&jet2;
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";
  
  TCut mdp = "minDeltaPhi20>0.3";
  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;

  TCut metsigbound="METsig<60";

  TCut drmax = "deltaRmax_hh<2.2";

  TCut metsigloose = "METsig>30";


  //for RA2b trigger, raise jetpt2 threshold                                                                                                                                       
  selection_ = baseline && ra2btrigger && zl&&isotk && jets&&jet2high && !btag3 &&!mdp;

 //  renewCanvas();
//   if (hdata2d!=0)  delete hdata2d;
//   hdata2d=new TH2D("hdata2d","hdata",100,0,500,100,0,500);
//   dtree->Draw("MET:METsig>>hdata2d",selection_,"colz");
//   hdata2d->SetXTitle("MET significance");
//   hdata2d->SetYTitle("MET (GeV)");

//   thecanvas->SaveAs("plots_Control_QCD/METversusMETsig_LDP_ra2bTrigger.pdf");

  logz_=true;
 renewCanvas();
  if (hdata2d!=0)  delete hdata2d;
  hdata2d=new TH2D("hdata2d","hdata",100,0,500,100,0,500);
  dtree->Draw("MET:METsig>>hdata2d",selection_,"colz");
  hdata2d->SetXTitle("MET significance");
  hdata2d->SetYTitle("MET (GeV)");

  //thecanvas->SaveAs("plots_Control_QCD/METversusMETsig_LDP_ra2bTrigger_logz.pdf");
  
  TProfile * hproj =   hdata2d->ProfileX("hproj");
  hproj->SetLineWidth(2);
  hproj->Draw("SAME");

  thecanvas->SaveAs("plots_Control_QCD/METversusMETsig_LDP_ra2bTrigger_logz.pdf");

  logz_=false;

  //then make a cut on low-ish MET

  selection_ = baseline && ra2btrigger && zl&&isotk && jets&&jet2high && !btag3 &&!mdp &&TCut("MET<100");
  setLogY(true);
  renewCanvas();
  if (hdata!=0) delete hdata;
  hdata=new TH1D("hdata","MET<100 GeV data",50,0,100);
  dtree->Draw("METsig>>hdata",selection_);
  hdata->SetXTitle("MET significance");

  thecanvas->SaveAs("plots_Control_QCD/METsig_LDP_ra2bTrigger_METlt100.pdf");

  setLogY(false);
  // now the DeltaPhi>0.3 sample with b-tagging

  selection_ = baseline && trigger && zl&&isotk && jets&&jet2 && btag2&&!btag3 &&mdp;

  logz_=true;
  renewCanvas();
  delete hdata2d;
  hdata2d=new TH2D("hdata2d","hdata",100,0,500,100,0,500);
  dtree->Draw("MET:METsig>>hdata2d",selection_,"colz");
  hdata2d->SetXTitle("MET significance");
  hdata2d->SetYTitle("MET (GeV)");

  //thecanvas->SaveAs("plots_Control_QCD/METversusMETsig_MDP.pdf");

 TProfile * hproj2 =   hdata2d->ProfileX("hproj2");
  hproj2->SetLineWidth(2);
  hproj2->Draw("SAME");
  thecanvas->SaveAs("plots_Control_QCD/METversusMETsig_MDP_logz_wproj.pdf");

  //now look at MC

  TChain* qcd=  getTree("QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1814_v69");
  //TChain* tt0=  getTree("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1848_v69");
TChain* tt=  getTree("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1799_v69");
 TChain* znn= getTree("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1832_v69");

 selection_ = baseline && trigger && zl&&isotk && jets&&jet2 && btag2&&!btag3 &&mdp;
 
 logz_=true;
  renewCanvas();
  delete hdata2d;
  hdata2d=new TH2D("hdata2d","ttbar",100,0,500,100,0,500);
  tt->Draw("MET:METsig>>hdata2d",getCutString(false),"colz");
  hdata2d->SetXTitle("MET significance");
  hdata2d->SetYTitle("MET (GeV)");

 TProfile * hprojtt =   hdata2d->ProfileX("hprojtt");
  hprojtt->SetLineWidth(2);
  hprojtt->Draw("SAME");
  thecanvas->SaveAs("plots_Control_QCD/METversusMETsig_MDP_ttbarMC.pdf");


 logz_=true;
  renewCanvas();
  delete hdata2d;
  hdata2d=new TH2D("hdata2d","znn",100,0,500,100,0,500);
  znn->Draw("MET:METsig>>hdata2d",getCutString(false),"colz");
  hdata2d->SetXTitle("MET significance");
  hdata2d->SetYTitle("MET (GeV)");
  //HERE

 TProfile * hprojznn =   hdata2d->ProfileX("hprojznn");
  hprojznn->SetLineWidth(2);
  hprojznn->Draw("SAME");
  thecanvas->SaveAs("plots_Control_QCD/METversusMETsig_MDP_znnMC.pdf");

 logz_=true;
  renewCanvas();
  delete hdata2d;
  hdata2d=new TH2D("hdata2d","qcd (mdp)",100,0,500,100,0,500);
  qcd->Draw("MET:METsig>>hdata2d",getCutString(false),"colz");
  hdata2d->SetXTitle("MET significance");
  hdata2d->SetYTitle("MET (GeV)");
//HERE
 TProfile * hprojqcd =   hdata2d->ProfileX("hprojqcd");
  hprojqcd->SetLineWidth(2);
  hprojqcd->Draw("SAME");
  thecanvas->SaveAs("plots_Control_QCD/METversusMETsig_MDP_qcdMC.pdf");

  selection_ = baseline && ra2btrigger && zl&&isotk && jets&&jet2high  &&!mdp;

 logz_=true;
  renewCanvas();
  delete hdata2d;
  hdata2d=new TH2D("hdata2d","qcd (ldp)",100,0,500,100,0,500);
  qcd->Draw("MET:METsig>>hdata2d",getCutString(false),"colz");
  hdata2d->SetXTitle("MET significance");
  hdata2d->SetYTitle("MET (GeV)");

  hdata2d->SetTitle("hqcd");

 TProfile * hprojqcd2 =   hdata2d->ProfileX("hprojqcd2");
  hprojqcd2->SetLineWidth(2);
  hprojqcd2->Draw("SAME");
  thecanvas->SaveAs("plots_Control_QCD/METversusMETsig_LDP_qcdMC.pdf");



  //  logz_=false;
  //  dodata_=true;
  //  draw2d("METsig",100,0,500,"MET",100,0,500,"MET significance","MET","plots_Control_QCD/METversusMETsig_MDP_dataMC");
  //  doRatio_=true;
  //  draw2d("METsig",50,0,500,"MET",50,0,500,"MET significance","MET","plots_Control_QCD/METversusMETsig_MDP_dataMC");

}
