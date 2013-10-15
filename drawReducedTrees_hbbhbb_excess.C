/*

.L TSelectorMultiDraw.C+
.L CrossSectionTable.cxx+
.L ConfigurationDescriptions.cxx+


.L drawReducedTrees_hbbhbb_excess.C+

*/

#include "drawReducedTrees_hbbhbb.C"

//to draw everything use default argument
void makeplots1(TString todraw="table pv mass 2b 3b 4b met0 met1 met2 met3 met4 sb sig") {

  //initHiggsSamples69(true,"ttbar znunu qcd wjets ttv vv singlet"); //all backgrounds
  initHiggsSamples69(true,"ttbar znunu bjets wjets ttv vv singlet hhmg300"); //all backgrounds

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

  //  TCut mdp = "minDeltaPhi20>0.5 || (minDeltaPhi20>0.3&&METsig>50)";
  TCut mdp = "minDeltaPhi20_eta5_noIdAll_nobeta>0.5 || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3 && METsig>50)";

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
  TCut nb2b = baseline && trigger && zl&& isotk && jets && btag2 && !btag3  && mdp && drmax;
  TCut nb3b = baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax;
  TCut nb4b = baseline && trigger && zl&& isotk && jets && btag2 && btag3  &&btag4 && mdp  && drmax;

  TCut sb2b = nb2b && higgsSB;
  TCut sb3b = nb3b && higgsSB;
  TCut sb4b = nb4b && higgsSB;


  TCut nb2bnodr = baseline && trigger && zl&& isotk && jets && btag2 && !btag3 && mdp ;
  TCut nb3bnodr = baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp ;
  TCut nb4bnodr = baseline && trigger && zl&& isotk && jets && btag2 && btag3  &&btag4 && mdp ;

  TCut nb2bldp = baseline && trigger && zl&& isotk && jets && btag2 && !btag3 && !mdp  && drmax;
  TCut nb3bldp = baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&!mdp  && drmax;
  TCut nb4bldp = baseline && trigger && zl&& isotk && jets && btag2 && btag3  &&btag4 && !mdp  && drmax;

  TCut nb2bsl = baseline && trigger && sl && jets && btag2 && !btag3 && mdp && drmax;
  TCut nb3bsl = baseline && trigger && sl && jets && btag2 && btag3  && !btag4&&mdp && drmax;
  TCut nb4bsl = baseline && trigger && sl && jets && btag2 && btag3  &&btag4 && mdp && drmax;


  TCut met0 = "METsig>=15 && METsig<30";
  TCut met1 = "METsig>=30 && METsig<50";
  TCut met2 = "METsig>=50 && METsig<100";
  TCut met3 = "METsig>=100 && METsig<150";
  TCut met4 = "METsig>=150";
  TString filename;

  TString dm = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)";
  TString am = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)";

  if (todraw.Contains("mass")) {

    // --- unblinded plots ---
    // average mass, in METsig bins, after making all other cuts
    //4b METsig1
    selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && btag4&&mdp && drmax && met1 && higgsSR_d;
    nbins=40; low=0; high=200;
    var=am; xtitle=var; filename = "4b_METsig1_dmSR_am";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    //4b METsig2
    selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && btag4&&mdp && drmax && met2 && higgsSR_d;
    nbins=40; low=0; high=200;
    var=am; xtitle=var; filename = "4b_METsig2_dmSR_am";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    //4b METsig3
    selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && btag4&&mdp && drmax && met3 && higgsSR_d;
    nbins=40; low=0; high=200;
    var=am; xtitle=var; filename = "4b_METsig3_dmSR_am";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    //4b METsig4
    selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && btag4&&mdp && drmax && met4 && higgsSR_d;
    nbins=40; low=0; high=200;
    var=am; xtitle=var; filename = "4b_METsig4_dmSR_am";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    
   //2b METsig1
    selection_= baseline && trigger && zl&& isotk && jets && btag2 && !btag3  &&mdp && drmax && met1 && higgsSR_d;
    nbins=40; low=0; high=200;
    var=am; xtitle=var; filename = "2b_METsig1_dmSR_am";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    //2b METsig2
    selection_= baseline && trigger && zl&& isotk && jets && btag2 && !btag3  &&mdp && drmax && met2 && higgsSR_d;
    nbins=40; low=0; high=200;
    var=am; xtitle=var; filename = "2b_METsig2_dmSR_am";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    //2b METsig3
    selection_= baseline && trigger && zl&& isotk && jets && btag2 && !btag3  &&mdp && drmax && met3 && higgsSR_d;
    nbins=40; low=0; high=200;
    var=am; xtitle=var; filename = "2b_METsig3_dmSR_am";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    //2b METsig4
    selection_= baseline && trigger && zl&& isotk && jets && btag2 && !btag3  &&mdp && drmax && met4 && higgsSR_d;
    nbins=40; low=0; high=200;
    var=am; xtitle=var; filename = "2b_METsig4_dmSR_am";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    

    // -- old blind plots ---
    
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
      nbins=40; low=0; high=200;
      var=am; xtitle=var; filename = "3b_METsig1_amSB_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 3b, METsig bin 2, ave_mjj excluding 90 < ave_mjj < 150 (for comparison with same for 4b)
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax && met2 && (higgsSB_avl||higgsSB_avh);
      nbins=40; low=0; high=200;
      var=am; xtitle=var; filename = "3b_METsig2_amSB_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 3b, METsig bin 1, ave_mjj excluding |delta_mjj|<30 (for comparison with same for 4b)
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax && met1 && higgsSB_d;
      nbins=40; low=0; high=200;
      var=am; xtitle=var; filename = "3b_METsig1_dmSB_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 3b, METsig bin 2, ave_mjj excluding |delta_mjj|<30 (for comparison with same for 4b)
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && drmax && met2 && higgsSB_d;
      nbins=40; low=0; high=200;
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
      nbins=40; low=0; high=200;
      var=am; xtitle=var; filename = "4b_METsig1_amSB_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 4b, METsig bin 2, ave_mjj excluding 90 < ave_mjj < 150
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && btag4&&mdp && drmax && met2 && (higgsSB_avl||higgsSB_avh);
      nbins=40; low=0; high=200;
      var=am; xtitle=var; filename = "4b_METsig2_amSB_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 4b, METsig bin 1, ave_mjj excluding |delta_mjj|<30
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && btag4&&mdp && drmax && met1 && higgsSB_d;
      nbins=40; low=0; high=200;
      var=am; xtitle=var; filename = "4b_METsig1_dmSB_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
//        - 4b, METsig bin 2, ave_mjj excluding |delta_mjj|<30
      selection_= baseline && trigger && zl&& isotk && jets && btag2 && btag3  && btag4&&mdp && drmax && met2 && higgsSB_d;
      nbins=40; low=0; high=200;
      var=am; xtitle=var; filename = "4b_METsig2_dmSB_am";
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

  }

  //pv studies
  if (todraw.Contains("pv")) {

    //can i add the mean and rms to these plots?

    setLogY(true);
    setPlotMinimum(0.1);

    // PV0 - beamspot in transverse plane
    nbins=60; low=0; high=0.03;
    var = "sqrt((PV0_x-BS0_x)*(PV0_x-BS0_x) + (PV0_y-BS0_y)*(PV0_y-BS0_y))"; xtitle="#sqrt{(PVx-BSx)^2 + (PVy-BSy)^2}";

    // 4b SB METsig1
    selection_ = sb4b && met1;
    filename = "SB4b_METsig1_transversePVdisplacement";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 problematic events
    selection_ = sb4b && met1 && TCut("deltaPhiStar<0.2");
    filename = "SB4b_METsig1_deltaPhiStarlt0p2_transversePVdisplacement";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 SL
    selection_ = nb4bsl&&higgsSB && met1;
    filename = "SL_SB4b_METsig1_transversePVdisplacement";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig2
    selection_ = sb4b && met2;
    filename = "SB4b_METsig2_transversePVdisplacement";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    

    // PV uncertainty
    nbins=40; low=0; high=0.01;
    var = "sqrt(PV0_xErr*PV0_xErr+PV0_yErr*PV0_yErr)"; xtitle="#sqrt{PVerrX^2 + PVerrY^2}";
     // 4b SB METsig1
    selection_ = sb4b && met1;
    filename = "SB4b_METsig1_transversePVerr";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 problematic events
    selection_ = sb4b && met1 && TCut("deltaPhiStar<0.2");
    filename = "SB4b_METsig1_deltaPhiStarlt0p2_transversePVerr";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 SL
    selection_ = nb4bsl&&higgsSB && met1;
    filename = "SL_SB4b_METsig1_transversePVerr";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig2
    selection_ = sb4b && met2;
    filename = "SB4b_METsig2_transversePVerr";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

    setLogY(false);    resetPlotMinimum();

    //distance from PV0 to the nearest PV in z
    nbins=25; low=0; high=1;
    var = "abs(PV0_z-PV_z[PV_nearestZindex])"; xtitle="| PV0_z - PVnearest_z |";
     // 4b SB METsig1
    selection_ = sb4b && met1;
    filename = "SB4b_METsig1_PV01dz";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 problematic events
    selection_ = sb4b && met1 && TCut("deltaPhiStar<0.2");
    filename = "SB4b_METsig1_deltaPhiStarlt0p2_PV01dz";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 SL
    selection_ = nb4bsl&&higgsSB && met1;
    filename = "SL_SB4b_METsig1_PV01dz";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig2
    selection_ = sb4b && met2;
    filename = "SB4b_METsig2_PV01dz";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

    //nGoodPV
    nbins=10; low=0; high=30; //purposefully use a coarse binning
    var = "nGoodPV"; xtitle=var;
     // 4b SB METsig1
    selection_ = sb4b && met1;
    filename = "SB4b_METsig1_nGoodPV";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 problematic events
    selection_ = sb4b && met1 && TCut("deltaPhiStar<0.2");
    filename = "SB4b_METsig1_deltaPhiStarlt0p2_nGoodPV";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 SL
    selection_ = nb4bsl&&higgsSB && met1;
    filename = "SL_SB4b_METsig1_nGoodPV";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig2
    selection_ = sb4b && met2;
    filename = "SB4b_METsig2_nGoodPV";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
     // 2b SB METsig1
    selection_ = sb2b && met1;
    filename = "SB2b_METsig1_nGoodPV";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

    //not sure about log scale or not
    setLogY(false);    resetPlotMinimum();

    //n tracks in PV0
    nbins=40; low=0; high=200;
    var="PV_tracksSize[0]"; xtitle=var;
     // 4b SB METsig1
    selection_ = sb4b && met1;
    filename = "SB4b_METsig1_PV0ntracks";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 problematic events
    selection_ = sb4b && met1 && TCut("deltaPhiStar<0.2");
    filename = "SB4b_METsig1_deltaPhiStarlt0p2_PV0ntracks";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 SL
    selection_ = nb4bsl&&higgsSB && met1;
    filename = "SL_SB4b_METsig1_PV0ntracks";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig2
    selection_ = sb4b && met2;
    filename = "SB4b_METsig2_PV0ntracks";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
     // 2b SB METsig1
    selection_ = sb2b && met1;
    filename = "SB2b_METsig1_PV0ntracks";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);


    //n tracks in PV1
    nbins=40; low=0; high=200;
    var="PV_tracksSize[1]"; xtitle=var; // i think no special cut is needed. underflow will just not appear on plot
     // 4b SB METsig1
    selection_ = sb4b && met1 ;
    filename = "SB4b_METsig1_PV1ntracks";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 problematic events
    selection_ = sb4b && met1 && TCut("deltaPhiStar<0.2");
    filename = "SB4b_METsig1_deltaPhiStarlt0p2_PV1ntracks";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 SL
    selection_ = nb4bsl&&higgsSB && met1;
    filename = "SL_SB4b_METsig1_PV1ntracks";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig2
    selection_ = sb4b && met2;
    filename = "SB4b_METsig2_PV1ntracks";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
     // 2b SB METsig1
    selection_ = sb2b && met1;
    filename = "SB2b_METsig1_PV1ntracks";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

    //signed difference between the two
    nbins=40; low=-100; high=100;
    var="PV_tracksSize[0]-PV_tracksSize[1]"; xtitle=var;
    TCut pv1exists = "PV_tracksSize[0]>=1 && PV_tracksSize[1]>=1"; //avoid weird results from dummy values
     // 4b SB METsig1
    selection_ = sb4b && met1 &&pv1exists;
    filename = "SB4b_METsig1_PV0minus1ntracks";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 problematic events
    selection_ = sb4b && met1 && TCut("deltaPhiStar<0.2")&&pv1exists;
    filename = "SB4b_METsig1_deltaPhiStarlt0p2_PV0minus1ntracks";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 SL
    selection_ = nb4bsl&&higgsSB && met1&&pv1exists;
    filename = "SL_SB4b_METsig1_PV0minus1ntracks";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig2
    selection_ = sb4b && met2&&pv1exists;
    filename = "SB4b_METsig2_PV0minus1ntracks";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
     // 2b SB METsig1
    selection_ = sb2b && met1&&pv1exists;
    filename = "SB2b_METsig1_PV0minus1ntracks";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

    //PV0 chi2

    nbins=50; low=0; high=300;
    var="PV_chi2[0]"; xtitle=var;
     // 4b SB METsig1
    selection_ = sb4b && met1;
    filename = "SB4b_METsig1_PV0chi2";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 problematic events
    selection_ = sb4b && met1 && TCut("deltaPhiStar<0.2");
    filename = "SB4b_METsig1_deltaPhiStarlt0p2_PV0chi2";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 SL
    selection_ = nb4bsl&&higgsSB && met1;
    filename = "SL_SB4b_METsig1_PV0chi2";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig2
    selection_ = sb4b && met2;
    filename = "SB4b_METsig2_PV0chi2";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
     // 2b SB METsig1
    selection_ = sb2b && met1;
    filename = "SB2b_METsig1_PV0chi2";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

    //PV0 ndof
    nbins=50; low=0; high=300;
    var="PV_ndof[0]"; xtitle=var;
     // 4b SB METsig1
    selection_ = sb4b && met1;
    filename = "SB4b_METsig1_PV0ndof";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 problematic events
    selection_ = sb4b && met1 && TCut("deltaPhiStar<0.2");
    filename = "SB4b_METsig1_deltaPhiStarlt0p2_PV0ndof";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 SL
    selection_ = nb4bsl&&higgsSB && met1;
    filename = "SL_SB4b_METsig1_PV0ndof";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig2
    selection_ = sb4b && met2;
    filename = "SB4b_METsig2_PV0ndof";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
     // 2b SB METsig1
    selection_ = sb2b && met1;
    filename = "SB2b_METsig1_PV0ndof";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    //PV0 chi2/ndof
    nbins=40; low=0; high=2;
    var="PV_chi2[0]/PV_ndof[0]"; xtitle=var;
    TCut div0protect="PV_ndof[0]>0";
     // 4b SB METsig1
    selection_ = sb4b && met1&&div0protect;
    filename = "SB4b_METsig1_PV0chi2PerNdof";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 problematic events
    selection_ = sb4b && met1&&div0protect && TCut("deltaPhiStar<0.2");
    filename = "SB4b_METsig1_deltaPhiStarlt0p2_PV0chi2PerNdof";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig1 SL
    selection_ = nb4bsl&&higgsSB && met1&&div0protect;
    filename = "SL_SB4b_METsig1_PV0chi2PerNdof";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    // 4b SB METsig2
    selection_ = sb4b && met2&&div0protect;
    filename = "SB4b_METsig2_PV0chi2PerNdof";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
     // 2b SB METsig1
    selection_ = sb2b && met1&&div0protect;
    filename = "SB2b_METsig1_PV0chi2PerNdof";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    //2b SB METsig1 LDP
    selection_ = nb2bldp&&higgsSB && met1&&div0protect;
    filename = "LDP_SB2b_METsig1_PV0chi2PerNdof";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
    //4b SB METsig1 LDP
    selection_ = nb4bldp&&higgsSB && met1&&div0protect;
    filename = "LDP_SB4b_METsig1_PV0chi2PerNdof";
    drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

  }

  if (todraw.Contains("table")) {
    //quantitative study of the "excess" events in the SBs
    savePlots_=false;
    nbins=10; low=0; high=1000; var="MET"; //var choice not very important
    xtitle="MET";
 
    //count data events in SB with deltaPhiStar<0.2
    //count non-QCD MC expectation
    selection_=sb2b && met1 && TCut("deltaPhiStar<0.2");
    drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
    double data_sb_2b = getIntegral("data");
    double data_sb_2b_err = sqrt(data_sb_2b);
    double nonqcd_sb_2b = totalnonqcd->Integral();
    double nonqcd_sb_2b_err = jmt::errOnIntegral(totalnonqcd);

    selection_=sb3b && met1 && TCut("deltaPhiStar<0.2");
    drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
    double data_sb_3b = getIntegral("data");
    double data_sb_3b_err = sqrt(data_sb_3b);
    double nonqcd_sb_3b = totalnonqcd->Integral();
    double nonqcd_sb_3b_err = jmt::errOnIntegral(totalnonqcd);

    selection_=sb4b && met1 && TCut("deltaPhiStar<0.2");
    drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0);
    double data_sb_4b = getIntegral("data");
    double data_sb_4b_err = sqrt(data_sb_4b);
    double nonqcd_sb_4b = totalnonqcd->Integral();
    double nonqcd_sb_4b_err = jmt::errOnIntegral(totalnonqcd);

    //quantify "excess" from QCD (with error)

    double excess_sb_2b = data_sb_2b - nonqcd_sb_2b;
    double excess_sb_3b = data_sb_3b - nonqcd_sb_3b;
    double excess_sb_4b = data_sb_4b - nonqcd_sb_4b;

    //calculate excess/nonqcd with error
    double ratio_sb_2b = excess_sb_2b / nonqcd_sb_2b;
    double ratio_sb_3b = excess_sb_3b / nonqcd_sb_3b;
    double ratio_sb_4b = excess_sb_4b / nonqcd_sb_4b;

    // (data-nonqcd)/nonqcd = data/nonqcd - 1
    double ratio_sb_2b_err = jmt::errAoverB(data_sb_2b,data_sb_2b_err,nonqcd_sb_2b,nonqcd_sb_2b_err);
    double ratio_sb_3b_err = jmt::errAoverB(data_sb_3b,data_sb_3b_err,nonqcd_sb_3b,nonqcd_sb_3b_err);
    double ratio_sb_4b_err = jmt::errAoverB(data_sb_4b,data_sb_4b_err,nonqcd_sb_4b,nonqcd_sb_4b_err);

    //use printf (god forbid)
    printf("2b  %d   %.1f   %.1f    %.2f +/- %.2f",int(data_sb_2b),nonqcd_sb_2b,excess_sb_2b,ratio_sb_2b,ratio_sb_2b_err);
    printf("3b  %d   %.1f   %.1f    %.2f +/- %.2f",int(data_sb_3b),nonqcd_sb_3b,excess_sb_3b,ratio_sb_3b,ratio_sb_3b_err);
    printf("4b  %d   %.1f   %.1f    %.2f +/- %.2f",int(data_sb_4b),nonqcd_sb_4b,excess_sb_4b,ratio_sb_4b,ratio_sb_4b_err);

    savePlots_=true;
  }


    TCut sel1;
    TCut sel2;
    TCut sel3;

    for (int k=0;k<2;k++) { //sb or sig

      for ( int i=2; i<=4; i++) {

	for ( int j=0; j<=4; j++) {

      if (k==0) sel3 = higgsSR;
      else if (k==1) sel3 = higgsSB;

      if (!todraw.Contains("sig") && k==0) continue;
      if (!todraw.Contains("sb") && k==1) continue;
      TString sel3name = k==0 ? "SIG" : "SB";

      if (i==2) sel1=nb2b;
      else if (i==3) sel1=nb3b;
      else if (i==4) sel1=nb4b;
      
      if (!todraw.Contains("2b") && i==2) continue;
      if (!todraw.Contains("3b") && i==3) continue;
      if (!todraw.Contains("4b") && i==4) continue;

      if (j==1) sel2=met1;
      else if (j==2) sel2=met2;
      else if (j==0) sel2=met0;
      else if (j==3) sel2=met3;
      else if (j==4) sel2=met4;

      if (!todraw.Contains("met0") && j==0) continue;
      if (!todraw.Contains("met1") && j==1) continue;
      if (!todraw.Contains("met2") && j==2) continue;
      if (!todraw.Contains("met3") && j==3) continue;
      if (!todraw.Contains("met4") && j==4) continue;

      selection_ = sel1&&sel2&&sel3;

      //CSV output
      setLogY(false);

      nbins=10; low=0; high=1;
      var="CSVbest1"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      var="CSVbest2"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      var="CSVbest3"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      setLogY(true);
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      setLogY(false);

      var="CSVbest4"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      setLogY(true);
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      setLogY(false);

      //jet pt

      nbins=20; low=0; high=400;
      var="jetpt1"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=400;
      var="jetpt2"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=200;
      var="jetpt3"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=200;
      var="jetpt4"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //njets
      nbins=9; low=0; high=9;
      var="njets20"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //n b jets
      nbins=7; low=0; high=7;
      var="nbjetsCSVL"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //out to eta of 5
      nbins=4; low=0; high=4;
      var="njets20_5p0-njets20"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,"nForwardJets");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //npujets
      nbins=9; low=0; high=9;
      var="nPUjets20"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //MET spectrum
      nbins=40; low=0; high=400;
      var="MET"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //METsig spectrum
      nbins=40; low=0; high=100;
      var="METsig"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //HT
      nbins=60; low=0; high=600;
      var="HT20"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //HT+MET
      nbins=60; low=0; high=800;
      var="HT20+MET"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,"HTplusMET");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //|MET-MHT|
      nbins=20; low=0; high=200;
      var="abs(MET-MHT)"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,"METminusMHT");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //rho_kt6PFJetsForIsolation
      nbins=20; low=0; high=40;
      var="rho_kt6PFJetsForIsolation"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);


      //max deltaPhi
      nbins=60; low=0; high=6;
      var="maxDeltaPhi20_eta5_noIdAll"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);



      // == 2nd round of plots ==
      //minDeltaPhiAll20
      nbins=30; low=0; high=3;
      var="minDeltaPhiAll20"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //minDeltaPhiAll20; no eta cut, etc
      nbins=30; low=0; high=3;
      var="minDeltaPhi20_eta5_noIdAll"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //deltaPhiStar
      nbins=30; low=0; high=3;
      var="deltaPhiStar"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=40; low=0; high=400;
      var="deltaPhiStar_badjet_pt"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=40; low=0; high=400;
      var="deltaPhiStar_badjet_estimatedPt"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);


      nbins=60; low=-300; high=300;
      var="deltaPhiStar_badjet_pt-deltaPhiStar_badjet_estimatedPt"; xtitle="dphi* badjet pt-estPt";
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,"dphistar_badjet_ptMinusEstPt");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);


      nbins=30; low=-TMath::Pi(); high=TMath::Pi();
      var="deltaPhiStar_badjet_phi"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //angle between this phi and MET
      nbins=30; low=0; high=TMath::Pi();
      var="acos(cos(deltaPhiStar_badjet_phi-METphi))"; xtitle="Delta Phi(badjet,MET)";
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,"deltaPhiBadjetMET");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);


      nbins=40; low=-4; high=4;
      var="deltaPhiStar_badjet_eta"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //some more properties
      nbins=30; low=0; high=30;
      var="deltaPhiStar_badjet_n60"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=50; low=-100; high=100;
      var="deltaPhiStar_badjet_chMult-deltaPhiStar_badjet_neuMult"; xtitle="charged - neutral mult";
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,"chMinusNeuMult");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=1;
      var="deltaPhiStar_badjet_PUbeta"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      // charged energy fraction
      nbins=20; low=0; high=1;
      var="(deltaPhiStar_badjet_chHadE+deltaPhiStar_badjet_chEmE)/deltaPhiStar_badjet_energy"; xtitle="chHad+chEM energy fraction";
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,"chargedEnergyFraction");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //neutral energy fraction
      nbins=20; low=0; high=1;
      var="(deltaPhiStar_badjet_neuHadE+deltaPhiStar_badjet_neuEmE)/deltaPhiStar_badjet_energy"; xtitle="neuHad+neuEM energy fraction";
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,"neutralEnergyFraction");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);


      // muon energy fraction
      setLogY(true);
      nbins=20; low=0; high=1;
      var="deltaPhiStar_badjet_chMuE/deltaPhiStar_badjet_energy"; xtitle="muon energy fraction";
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,"muonEnergyFraction");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      setLogY(false);

      //slice it another way: HCAL only fractions
      // charged 
      nbins=20; low=0; high=1;
      var="(deltaPhiStar_badjet_chHadE)/deltaPhiStar_badjet_energy"; xtitle="chHad energy fraction";
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,"chargedHadEnergyFraction");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //neutral
      nbins=20; low=0; high=1;
      var="(deltaPhiStar_badjet_neuHadE)/deltaPhiStar_badjet_energy"; xtitle="neuHad energy fraction";
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,"neutralHadEnergyFraction");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);


      //leptons without isolation
      nbins=4; low=0; high=4;
      var="nMuonsNoRelIso"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=4; low=0; high=4;
      var="nElectronsNoRelIso"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      

      //remove these for now
/*
      //3rd round
      //look at jets with maxDeltaPhi20_eta5_noIdAll>2.9
      selection_ = sel1&&sel2 && TCut("maxDeltaPhi20_eta5_noIdAll>2.9");

      //log scale?
      nbins=20; low=0; high=400;
      var="jetpt1"; xtitle=var;
      filename.Form("%s%db_METsig%d_maxDeltaPhi20_eta5_noIdAllgt2p9_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=400;
      var="jetpt2"; xtitle=var;
      filename.Form("%s%db_METsig%d_maxDeltaPhi20_eta5_noIdAllgt2p9_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=200;
      var="jetpt3"; xtitle=var;
      filename.Form("%s%db_METsig%d_maxDeltaPhi20_eta5_noIdAllgt2p9_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=200;
      var="jetpt4"; xtitle=var;
      filename.Form("%s%db_METsig%d_maxDeltaPhi20_eta5_noIdAllgt2p9_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //eta
      nbins=20; low=-4; high=4;
      var="jeteta1"; xtitle=var;
      filename.Form("%s%db_METsig%d_maxDeltaPhi20_eta5_noIdAllgt2p9_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      var="jeteta2"; xtitle=var;
      filename.Form("%s%db_METsig%d_maxDeltaPhi20_eta5_noIdAllgt2p9_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      var="jeteta3"; xtitle=var;
      filename.Form("%s%db_METsig%d_maxDeltaPhi20_eta5_noIdAllgt2p9_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      var="jeteta4"; xtitle=var;
      filename.Form("%s%db_METsig%d_maxDeltaPhi20_eta5_noIdAllgt2p9_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //phi
      nbins=20; low=-TMath::Pi(); high=TMath::Pi();
      var="jetphi1"; xtitle=var;
      filename.Form("%s%db_METsig%d_maxDeltaPhi20_eta5_noIdAllgt2p9_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      var="jetphi2"; xtitle=var;
      filename.Form("%s%db_METsig%d_maxDeltaPhi20_eta5_noIdAllgt2p9_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      var="jetphi3"; xtitle=var;
      filename.Form("%s%db_METsig%d_maxDeltaPhi20_eta5_noIdAllgt2p9_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      var="jetphi4"; xtitle=var;
      filename.Form("%s%db_METsig%d_maxDeltaPhi20_eta5_noIdAllgt2p9_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
*/

      //N-1 plots of drmax
     
      if (i==2) selection_ = nb2bnodr && sel2 && sel3;
      else if (i==3) selection_ = nb3bnodr && sel2 &&sel3;
      else if (i==4) selection_ = nb4bnodr && sel2 &&sel3;
 
      //n b jets
      nbins=40; low=0; high=6;
      var="deltaRmax_hh"; xtitle=var;
      filename.Form("%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);


      //4th round
      //look at deltaPhiStar<0.2
      selection_ = sel1&&sel2&&sel3 && TCut("deltaPhiStar<0.2");

      //n b jets
      nbins=7; low=0; high=7;
      var="nbjetsCSVL"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=30; low=0; high=3;
      var="minDeltaPhiAll20"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

    //rho_kt6PFJetsForIsolation
      nbins=20; low=0; high=40;
      var="rho_kt6PFJetsForIsolation"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //minDeltaPhiAll20; no eta cut, etc
      nbins=30; low=0; high=3;
      var="minDeltaPhi20_eta5_noIdAll"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=40; low=0; high=400;
      var="deltaPhiStar_badjet_pt"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=30; low=-TMath::Pi(); high=TMath::Pi();
      var="deltaPhiStar_badjet_phi"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //angle between this phi and MET
      nbins=30; low=0; high=TMath::Pi();
      var="acos(cos(deltaPhiStar_badjet_phi-METphi))"; xtitle="Delta Phi(badjet,MET)";
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,"deltaPhiBadjetMET");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=60; low=-300; high=300;
      var="deltaPhiStar_badjet_pt-deltaPhiStar_badjet_estimatedPt"; xtitle="dphi* badjet pt-estPt";
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,"dphistar_badjet_ptMinusEstPt");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);


      nbins=40; low=-4; high=4;
      var="deltaPhiStar_badjet_eta"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=1;
      var="deltaPhiStar_badjet_CSV"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=1;
      var="deltaPhiStar_badjet_PUbeta"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=80;
      var="deltaPhiStar_badjet_chMult"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=80;
      var="deltaPhiStar_badjet_neuMult"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=50; low=-100; high=100;
      var="deltaPhiStar_badjet_chMult-deltaPhiStar_badjet_neuMult"; xtitle="charged - neutral mult";
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,"chMinusNeuMult");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=30; low=0; high=30;
      var="deltaPhiStar_badjet_n60"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=30; low=0; high=60;
      var="deltaPhiStar_badjet_n90"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);


      // charged energy fraction
      nbins=20; low=0; high=1;
      var="(deltaPhiStar_badjet_chHadE+deltaPhiStar_badjet_chEmE)/deltaPhiStar_badjet_energy"; xtitle="chHad+chEM energy fraction";
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,"chargedEnergyFraction");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //neutral energy fraction
      nbins=20; low=0; high=1;
      var="(deltaPhiStar_badjet_neuHadE+deltaPhiStar_badjet_neuEmE)/deltaPhiStar_badjet_energy"; xtitle="neuHad+neuEM energy fraction";
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,"neutralEnergyFraction");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);


      // muon energy fraction
      setLogY(true);
      nbins=20; low=0; high=1;
      var="deltaPhiStar_badjet_chMuE/deltaPhiStar_badjet_energy"; xtitle="muon energy fraction";
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,"muonEnergyFraction");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
      setLogY(false);

      //slice it another way: HCAL only fractions
      // charged 
      nbins=20; low=0; high=1;
      var="(deltaPhiStar_badjet_chHadE)/deltaPhiStar_badjet_energy"; xtitle="chHad energy fraction";
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,"chargedHadEnergyFraction");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //neutral
      nbins=20; low=0; high=1;
      var="(deltaPhiStar_badjet_neuHadE)/deltaPhiStar_badjet_energy"; xtitle="neuHad energy fraction";
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,"neutralHadEnergyFraction");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //max deltaPhi
      nbins=60; low=0; high=6;
      var="maxDeltaPhi20_eta5_noIdAll"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //npujets
      nbins=9; low=0; high=9;
      var="nPUjets20"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //out to eta of 5
      nbins=4; low=0; high=4;
      var="njets20_5p0-njets20"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,"nForwardJets");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //METsig spectrum
      nbins=40; low=0; high=100;
      var="METsig"; xtitle=var;
      filename.Form("%s%db_METsig%d_deltaPhiStarlt0p2_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);


      //5th round 
      //look at SL sample
      if (i==2) sel1=nb2bsl;
      else if (i==3) sel1=nb3bsl;
      else if (i==4) sel1=nb4bsl;

      selection_ = sel1&&sel2&&sel3;
      //deltaPhiStar
      nbins=30; low=0; high=3;
      var="deltaPhiStar"; xtitle=var;
      filename.Form("SL_%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=30; low=0; high=30;
      var="deltaPhiStar_badjet_n60"; xtitle=var;
      filename.Form("SL_%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=50; low=-100; high=100;
      var="deltaPhiStar_badjet_chMult-deltaPhiStar_badjet_neuMult"; xtitle="charged - neutral mult";
      filename.Form("SL_%s%db_METsig%d_%s",sel3name.Data(),i,j,"chMinusNeuMult");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=20; low=0; high=1;
      var="deltaPhiStar_badjet_PUbeta"; xtitle=var;
      filename.Form("SL_%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);


      //max deltaPhi
      nbins=60; low=0; high=6;
      var="maxDeltaPhi20_eta5_noIdAll"; xtitle=var;
      filename.Form("SL_%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //min deltaPhi
      nbins=30; low=0; high=3;
      var="minDeltaPhi20_eta5_noIdAll"; xtitle=var;
      filename.Form("SL_%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //min deltaPhi
      nbins=30; low=0; high=3;
      var="minDeltaPhiAll20"; xtitle=var;
      filename.Form("SL_%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);


      //npujets
      nbins=9; low=0; high=9;
      var="nPUjets20"; xtitle=var;
      filename.Form("SL_%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //out to eta of 5
      nbins=4; low=0; high=4;
      var="njets20_5p0-njets20"; xtitle=var;
      filename.Form("SL_%s%db_METsig%d_%s",sel3name.Data(),i,j,"nForwardJets");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //METsig spectrum
      nbins=40; low=0; high=100;
      var="METsig"; xtitle=var;
      filename.Form("SL_%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //n b jets
      nbins=7; low=0; high=7;
      var="nbjetsCSVL"; xtitle=var;
      filename.Form("SL_%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

    //rho_kt6PFJetsForIsolation
      nbins=20; low=0; high=40;
      var="rho_kt6PFJetsForIsolation"; xtitle=var;
      filename.Form("SL_%s%db_METsig%d_%s",sel3name.Data(),i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      nbins=60; low=-300; high=300;
      var="deltaPhiStar_badjet_pt-deltaPhiStar_badjet_estimatedPt"; xtitle="dphi* badjet pt-estPt";
      filename.Form("SL_%s%db_METsig%d_%s",sel3name.Data(),i,j,"dphistar_badjet_ptMinusEstPt");
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

	}
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

//duplicate a bit of the above -- special version that adds some signal curves
void makeplots1signal() {

  //initHiggsSamples69(true,"ttbar znunu qcd wjets ttv vv singlet"); //all backgrounds
  initHiggsSamples69(true,"ttbar znunu bjets wjets ttv vv singlet hhmg175"); //all backgrounds

  setOutputDirectory("plots_unblind");
  stackSignal_=true;

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

  TCut sb2bnodr = baseline && trigger && zl&& isotk && jets && btag2 && !btag3 && mdp && higgsSB;
  TCut sb3bnodr = baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&mdp && higgsSB;
  TCut sb4bnodr = baseline && trigger && zl&& isotk && jets && btag2 && btag3  &&btag4 && mdp && higgsSB;

  TCut sb2bldp = baseline && trigger && zl&& isotk && jets && btag2 && !btag3 && !mdp && higgsSB && drmax;
  TCut sb3bldp = baseline && trigger && zl&& isotk && jets && btag2 && btag3  && !btag4&&!mdp && higgsSB && drmax;
  TCut sb4bldp = baseline && trigger && zl&& isotk && jets && btag2 && btag3  &&btag4 && !mdp && higgsSB && drmax;

  TCut sb2bsl = baseline && trigger && sl && jets && btag2 && !btag3 && mdp && higgsSB && drmax;
  TCut sb3bsl = baseline && trigger && sl && jets && btag2 && btag3  && !btag4&&mdp && higgsSB && drmax;
  TCut sb4bsl = baseline && trigger && sl && jets && btag2 && btag3  &&btag4 && mdp && higgsSB && drmax;

  TCut met0 = "METsig>=15 && METsig<30";
  TCut met1 = "METsig>=30 && METsig<50";
  TCut met2 = "METsig>=50 && METsig<100";
      TString filename;

  TString dm = "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)";
  TString am = "0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)";


    TCut sel1;
      TCut sel2;


  for ( int i=2; i<=4; i++) {

    for ( int j=1; j<=2; j++) {
      if (i==2) sel1=sb2b;
      else if (i==3) sel1=sb3b;
      else if (i==4) sel1=sb4b;

/*      
      if (!todraw.Contains("2b") && i==2) continue;
      if (!todraw.Contains("3b") && i==3) continue;
      if (!todraw.Contains("4b") && i==4) continue;
*/
      if (j==1) sel2=met1;
      else if (j==2) sel2=met2;
      else if (j==0) sel2=met0;

/*
      if (!todraw.Contains("met0") && j==0) continue;
      if (!todraw.Contains("met1") && j==1) continue;
      if (!todraw.Contains("met2") && j==2) continue;
*/

      selection_ = sel1&&sel2;

      //CSV output
      setLogY(false);

      nbins=10; low=0; high=1;
      var="CSVbest4"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s-withSignal",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

/* need new reducedTrees from Pawandeep
      //HT
      nbins=60; low=0; high=600;
      var="HT20"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s-withSignal",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

      //minDeltaPhiAll20
      nbins=30; low=0; high=3;
      var="minDeltaPhiAll20"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s-withSignal",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);
*/
      //deltaPhiStar
      nbins=50; low=0; high=2.5;
      var="deltaPhiStar"; xtitle=var;
      filename.Form("SB%db_METsig%d_%s-withSignal",i,j,var.Data());
      drawPlots(var,nbins,low,high,xtitle,"Events", filename,0);

    }
  }

}
