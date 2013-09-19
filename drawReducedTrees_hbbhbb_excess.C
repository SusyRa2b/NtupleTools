/*

.L TSelectorMultiDraw.C+
.L CrossSectionTable.cxx+
.L ConfigurationDescriptions.cxx+


.L drawReducedTrees_hbbhbb_excess.C+

*/

#include "drawReducedTrees_hbbhbb.C"

void makeplots1() {

  //initHiggsSamples69(true,"ttbar znunu qcd wjets ttv vv singlet"); //all backgrounds
  initHiggsSamples69(true,"ttbar znunu bjets wjets ttv vv singlet"); //all backgrounds

  setOutputDirectory("plots_unblind");

  //use top pT weights
  setSampleWeightFactor("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883ra2b_v71s","topPtWeight");
  setSampleWeightFactor("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884ra2b_v71s","topPtWeight");
  setSampleWeightFactor("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880ra2b_v71s",   "topPtWeight");

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoPFcands10_010==0";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut jets=(njets4||njets5) &&jet2;
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut mdp = "minDeltaPhi20>0.5 || (minDeltaPhi20>0.3&&METsig>50)";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;

  TCut metsigloose = "METsig>30";

  TCut drmax = "deltaRmax_hh<2.2";

  setStackMode(true,false,false); //stack,norm,label override
  setLogY(false);

  //define 2b, 3b, 4b SB
  TCut sb2b = baseline && trigger && zl&& isotk && jets && btag2 && !btag3 && mdp && higgsSB && drmax;
  TCut sb3b = baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && higgsSB && drmax;
  TCut sb4b = baseline && trigger && zl&& isotk && jets && btag2 && btag3  &&btag4 && mdp && higgsSB && drmax;

  TCut sb2bsl = baseline && trigger && sl && jets && btag2 && !btag3 && mdp && higgsSB && drmax;
  TCut sb3bsl = baseline && trigger && sl && jets && btag2 && btag3  && !btag4&&mdp && higgsSB && drmax;
  TCut sb4bsl = baseline && trigger && sl && jets && btag2 && btag3  &&btag4 && mdp && higgsSB && drmax;


  TCut met1 = "METsig>=30 && METsig<50";
  TCut met2 = "METsig>=50 && METsig<100";
      TString filename;

  TString dm = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)";
  TString am = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)";

//        - 3b, METsig bin 1, |delta_mjj| with no cut on ave_mjj
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax && met1;
      nbins=40; low=0; high=100;
      var=dm; xtitle=var; filename = "3b_METsig1_dm";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

//        - 3b, METsig bin 2, |delta_mjj| with no cut on ave_mjj
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax && met2;
      nbins=40; low=0; high=100;
      var=dm; xtitle=var; filename = "3b_METsig2_dm";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

//        - 3b, METsig bin 1, ave_mjj with no cut on |delta_mjj|
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax && met1;
      nbins=40; low=0; high=200;
      var=am; xtitle=var; filename = "3b_METsig1_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 3b, METsig bin 2, ave_mjj with no cut on |delta_mjj|
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax && met2;
      nbins=40; low=0; high=200;
      var=am; xtitle=var; filename = "3b_METsig2_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 3b, METsig bin 1, |delta_mjj| excluding 90 < ave_mjj < 150 (for comparison with same for 4b)
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax && met1 && (higgsSB_avl||higgsSB_avh);
      nbins=40; low=0; high=100;
      var=dm; xtitle=var; filename = "3b_METsig1_amSB_dm";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 3b, METsig bin 2, |delta_mjj| excluding 90 < ave_mjj < 150 (for comparison with same for 4b)
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax && met2 && (higgsSB_avl||higgsSB_avh);
      nbins=40; low=0; high=100;
      var=dm; xtitle=var; filename = "3b_METsig2_amSB_dm";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 3b, METsig bin 1, |delta_mjj| excluding |delta_mjj|<30 (for comparison with same for 4b)
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax && met1 && higgsSB_d;
      nbins=40; low=0; high=100;
      var=dm; xtitle=var; filename = "3b_METsig1_dmSB_dm";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 3b, METsig bin 2, |delta_mjj| excluding |delta_mjj|<30 (for comparison with same for 4b)
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax && met2 && higgsSB_d;
      nbins=40; low=0; high=100;
      var=dm; xtitle=var; filename = "3b_METsig2_dmSB_dm";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 3b, METsig bin 1, ave_mjj excluding 90 < ave_mjj < 150 (for comparison with same for 4b)
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax && met1 && (higgsSB_avl||higgsSB_avh);
      nbins=40; low=0; high=100;
      var=am; xtitle=var; filename = "3b_METsig1_amSB_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 3b, METsig bin 2, ave_mjj excluding 90 < ave_mjj < 150 (for comparison with same for 4b)
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax && met2 && (higgsSB_avl||higgsSB_avh);
      nbins=40; low=0; high=100;
      var=am; xtitle=var; filename = "3b_METsig2_amSB_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 3b, METsig bin 1, ave_mjj excluding |delta_mjj|<30 (for comparison with same for 4b)
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax && met1 && higgsSB_d;
      nbins=40; low=0; high=100;
      var=am; xtitle=var; filename = "3b_METsig1_dmSB_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 3b, METsig bin 2, ave_mjj excluding |delta_mjj|<30 (for comparison with same for 4b)
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax && met2 && higgsSB_d;
      nbins=40; low=0; high=100;
      var=am; xtitle=var; filename = "3b_METsig2_dmSB_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

//        - 4b, METsig bin 1, |delta_mjj| excluding 90 < ave_mjj < 150
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && btag4&&mdp && drmax && met1 && (higgsSB_avl||higgsSB_avh);
      nbins=40; low=0; high=100;
      var=dm; xtitle=var; filename = "4b_METsig1_amSB_dm";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 4b, METsig bin 2, |delta_mjj| excluding 90 < ave_mjj < 150
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && btag4&&mdp && drmax && met2 && (higgsSB_avl||higgsSB_avh);
      nbins=40; low=0; high=100;
      var=dm; xtitle=var; filename = "4b_METsig2_amSB_dm";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 4b, METsig bin 1, |delta_mjj| excluding |delta_mjj|<30
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && btag4&&mdp && drmax && met1 && higgsSB_d;
      nbins=40; low=0; high=100;
      var=dm; xtitle=var; filename = "4b_METsig1_dmSB_dm";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 4b, METsig bin 2, |delta_mjj| excluding |delta_mjj|<30
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && btag4&&mdp && drmax && met2 && higgsSB_d;
      nbins=40; low=0; high=100;
      var=dm; xtitle=var; filename = "4b_METsig2_dmSB_dm";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 4b, METsig bin 1, ave_mjj excluding 90 < ave_mjj < 150
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && btag4&&mdp && drmax && met1 && (higgsSB_avl||higgsSB_avh);
      nbins=40; low=0; high=100;
      var=am; xtitle=var; filename = "4b_METsig1_amSB_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 4b, METsig bin 2, ave_mjj excluding 90 < ave_mjj < 150
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && btag4&&mdp && drmax && met2 && (higgsSB_avl||higgsSB_avh);
      nbins=40; low=0; high=100;
      var=am; xtitle=var; filename = "4b_METsig2_amSB_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 4b, METsig bin 1, ave_mjj excluding |delta_mjj|<30
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && btag4&&mdp && drmax && met1 && higgsSB_d;
      nbins=40; low=0; high=100;
      var=am; xtitle=var; filename = "4b_METsig1_dmSB_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 4b, METsig bin 2, ave_mjj excluding |delta_mjj|<30
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && btag4&&mdp && drmax && met2 && higgsSB_d;
      nbins=40; low=0; high=100;
      var=am; xtitle=var; filename = "4b_METsig2_dmSB_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      const bool redoplots=true;

  for ( int i=2; i<=4; i++) {
    TCut sel1;
    if (i==2) sel1=sb2b;
    else if (i==3) sel1=sb3b;
    else if (i==4) sel1=sb4b;

    for ( int j=1; j<=2; j++) {
      TCut sel2;
      if (j==1) sel2=met1;
      else if (j==2) sel2=met2;

      selection_ = sel1&&sel2;

      //CSV output
      setLogY(false);

      if (redoplots) { //1st round of plots
      nbins=10; low=0; high=1;
      var="CSVbest1"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      var="CSVbest2"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      var="CSVbest3"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      setLogY(true);
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      setLogY(false);

      var="CSVbest4"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      setLogY(true);
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      setLogY(false);

      //jet pt

      nbins=20; low=0; high=400;
      var="jetpt1"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=400;
      var="jetpt2"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=200;
      var="jetpt3"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=200;
      var="jetpt4"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //njets
      nbins=9; low=0; high=9;
      var="njets20"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //MET spectrum
      nbins=40; low=0; high=400;
      var="MET"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);


      //max deltaPhi
      nbins=60; low=0; high=6;
      var="maxDeltaPhiAll30"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      } //1st round of plots

      // == 2nd round of plots ==
      if (redoplots) {
      //deltaPhiStar
      nbins=30; low=0; high=3;
      var="deltaPhiStar"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=40; low=0; high=400;
      var="deltaPhiStar_badjet_pt"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=30; low=-TMath::Pi(); high=TMath::Pi();
      var="deltaPhiStar_badjet_phi"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=40; low=-4; high=4;
      var="deltaPhiStar_badjet_eta"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //leptons without isolation
      nbins=4; low=0; high=4;
      var="nMuonsNoRelIso"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=4; low=0; high=4;
      var="nElectronsNoRelIso"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      

      //3rd round
      //look at jets with maxDeltaPhiAll30>2.9
      selection_ = sel1&&sel2 && TCut("maxDeltaPhiAll30>2.9");

      //log scale?
      nbins=20; low=0; high=400;
      var="jetpt1"; xtitle=var;
      filename.Form("SB%db_METsig%d_maxDeltaPhiAll30gt2p9_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=400;
      var="jetpt2"; xtitle=var;
      filename.Form("SB%db_METsig%d_maxDeltaPhiAll30gt2p9_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=200;
      var="jetpt3"; xtitle=var;
      filename.Form("SB%db_METsig%d_maxDeltaPhiAll30gt2p9_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=200;
      var="jetpt4"; xtitle=var;
      filename.Form("SB%db_METsig%d_maxDeltaPhiAll30gt2p9_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //eta
      nbins=20; low=-4; high=4;
      var="jeteta1"; xtitle=var;
      filename.Form("SB%db_METsig%d_maxDeltaPhiAll30gt2p9_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      var="jeteta2"; xtitle=var;
      filename.Form("SB%db_METsig%d_maxDeltaPhiAll30gt2p9_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      var="jeteta3"; xtitle=var;
      filename.Form("SB%db_METsig%d_maxDeltaPhiAll30gt2p9_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      var="jeteta4"; xtitle=var;
      filename.Form("SB%db_METsig%d_maxDeltaPhiAll30gt2p9_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //phi
      nbins=20; low=-TMath::Pi(); high=TMath::Pi();
      var="jetphi1"; xtitle=var;
      filename.Form("SB%db_METsig%d_maxDeltaPhiAll30gt2p9_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      var="jetphi2"; xtitle=var;
      filename.Form("SB%db_METsig%d_maxDeltaPhiAll30gt2p9_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      var="jetphi3"; xtitle=var;
      filename.Form("SB%db_METsig%d_maxDeltaPhiAll30gt2p9_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      var="jetphi4"; xtitle=var;
      filename.Form("SB%db_METsig%d_maxDeltaPhiAll30gt2p9_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //4th round
      //look at deltaPhiStar<0.2
      selection_ = sel1&&sel2 && TCut("deltaPhiStar<0.2");

      nbins=40; low=0; high=400;
      var="deltaPhiStar_badjet_pt"; xtitle=var;
      filename.Form("SB%db_METsig%d_deltaPhiStarlt0.2_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=30; low=-TMath::Pi(); high=TMath::Pi();
      var="deltaPhiStar_badjet_phi"; xtitle=var;
      filename.Form("SB%db_METsig%d_deltaPhiStarlt0.2_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=40; low=-4; high=4;
      var="deltaPhiStar_badjet_eta"; xtitle=var;
      filename.Form("SB%db_METsig%d_deltaPhiStarlt0.2_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      }

      //5th round 
      //look at SL sample
      if (i==2) sel1=sb2bsl;
      else if (i==3) sel1=sb3bsl;
      else if (i==4) sel1=sb4bsl;

      selection_ = sel1&&sel2;
      //deltaPhiStar
      nbins=30; low=0; high=3;
      var="deltaPhiStar"; xtitle=var;
      filename.Form("SL_SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //max deltaPhi
      nbins=60; low=0; high=6;
      var="maxDeltaPhiAll30"; xtitle=var;
      filename.Form("SL_SB%db_METsig%d_%s",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);


    }
  }

}

void makeplots2() {

  //initHiggsSamples69(true,"ttbar znunu qcd wjets ttv vv singlet"); //all backgrounds
  initHiggsSamples69(true,"ttbar znunu bjets wjets ttv vv singlet"); //all backgrounds

  setOutputDirectory("plots_unblind");

  //use top pT weights
  setSampleWeightFactor("TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883ra2b_v71s","topPtWeight");
  setSampleWeightFactor("TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884ra2b_v71s","topPtWeight");
  setSampleWeightFactor("TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880ra2b_v71s",   "topPtWeight");

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true; ratioMin = 0; ratioMax = 2.2;
  dodata_=true;

  //define useful cuts
  TCut baseline = "cutPV==1 &&passCleaning==1 &&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40"; 

  TCut triggerJetMET = "passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiCentralPFJet30_PFMHT80==1";
  TCut triggerMET = "passMC_PFMET150==1";
  TCut trigger = triggerJetMET||triggerMET;

  TCut zl = "nMuons==0&&nElectrons==0&&nTausLoose==0";
  TCut isotk="nIsoPFcands10_010==0";

  TCut sl = "nMuons+nElectrons==1 &&nTausLoose==0";

  TCut njets4="njets20==4";
  TCut njets5="njets20==5";
  TCut jet2="jetpt2>50";
  TCut jets=(njets4||njets5) &&jet2;
  TCut btag2="CSVbest2>0.898";
  TCut btag3="CSVbest3>0.679";
  TCut btag4="CSVbest4>0.244";

  TCut mdp = "minDeltaPhi20>0.5 || (minDeltaPhi20>0.3&&METsig>50)";

  TCut higgsSR_av = "(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)";
  TCut higgsSR_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20";
  TCut higgsSR = higgsSR_av && higgsSR_d;

  TCut higgsSB_avl = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90";
  TCut higgsSB_avh = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150";
  TCut higgsSB_d = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30";
  TCut higgsSB = higgsSB_avl||higgsSB_avh||higgsSB_d;

  TCut metsigloose = "METsig>30";

  TCut drmax = "deltaRmax_hh<2.2";

  setStackMode(true,false,false); //stack,norm,label override
  setLogY(false);

  //define 2b, 3b, 4b SB
  TCut sb2b = baseline && trigger && zl&& isotk && jets && btag2 && !btag3 && mdp && higgsSB && drmax;
  TCut sb3b = baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && higgsSB && drmax;
  TCut sb4b = baseline && trigger && zl&& isotk && jets && btag2 && btag3  &&btag4 && mdp && higgsSB && drmax;

  TCut all2b = baseline && trigger && zl&& isotk && jets && btag2 && !btag3 && mdp && drmax;
  TCut all3b = baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax;
  TCut all4b = baseline && trigger && zl&& isotk && jets && btag2 && btag3  &&btag4 && mdp && drmax;

  TCut sb2bsl = baseline && trigger && sl && jets && btag2 && !btag3 && mdp && higgsSB && drmax;
  TCut sb3bsl = baseline && trigger && sl && jets && btag2 && btag3  && !btag4&&mdp && higgsSB && drmax;
  TCut sb4bsl = baseline && trigger && sl && jets && btag2 && btag3  &&btag4 && mdp && higgsSB && drmax;

  TCut met1 = "METsig>=30 && METsig<50";
  TCut met2 = "METsig>=50 && METsig<100";
  TString filename;

  TString dm = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)";
  TString am = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)";

 doRatio_=false; 
//        - 3b, METsig bin 1, |delta_mjj| vs ave_mjj
  selection_=all3b && met1;
  draw2d(dm,50,0,100,am,50,0,200,dm,am,"3b_METsig1_dm_am");
//        - 3b, METsig bin 2, |delta_mjj| vs ave_mjj
  selection_=all3b && met2;
  draw2d(dm,50,0,100,am,50,0,200,dm,am,"3b_METsig2_dm_am");
//        - 4b, METsig bin 1, |delta_mjj| vs ave_mjj, only for the SB (sig box + buffer removed).
  selection_=sb4b && met1;
  draw2d(dm,50,0,100,am,50,0,200,dm,am,"4b_METsig1_dm_am");
//        - 4b, METsig bin 2, |delta_mjj| vs ave_mjj, only for the SB (sig box + buffer removed).
  selection_=sb4b && met2;
  draw2d(dm,50,0,100,am,50,0,200,dm,am,"4b_METsig2_dm_am");


}
