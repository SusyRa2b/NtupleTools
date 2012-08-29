//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 14 10:52:54 2012 by ROOT version 5.33/02
// from TTree reducedTree/tree with minimal cuts
// found on file: /cu3/joshmt/reducedTrees/V00-02-35v/reducedTree.CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.T1bbbb.root
//////////////////////////////////////////////////////////

#ifndef ISRweights_h
#define ISRweights_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TH1F.h>
#include <TH2D.h>
#include <map>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class ISRweights {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        weight;
   Double_t        scanCrossSection;
   Double_t        scanCrossSectionPlus;
   Double_t        scanCrossSectionMinus;
   ULong64_t       runNumber;
   ULong64_t       lumiSection;
   ULong64_t       eventNumber;
   Int_t           m0;
   Int_t           m12;
   Int_t           W1decayType;
   Int_t           W2decayType;
   Int_t           decayType;
   Float_t         btagIPweight;
   Float_t         PUweight;
   Float_t         hltHTeff;
   Float_t         hltMHTeff;
   Float_t         hltMHTeffBNN;
   Float_t         hltMHTeffBNNUp;
   Float_t         hltMHTeffBNNDown;
   Float_t         pdfWeightsCTEQ[45];
   Float_t         pdfWeightsMSTW[41];
   Float_t         pdfWeightsNNPDF[100];
   Float_t         prob0;
   Float_t         probge1;
   Float_t         prob1;
   Float_t         probge2;
   Float_t         prob2;
   Float_t         probge3;
   Float_t         prob0_HFplus;
   Float_t         probge1_HFplus;
   Float_t         prob1_HFplus;
   Float_t         probge2_HFplus;
   Float_t         prob2_HFplus;
   Float_t         probge3_HFplus;
   Float_t         prob0_HFminus;
   Float_t         probge1_HFminus;
   Float_t         prob1_HFminus;
   Float_t         probge2_HFminus;
   Float_t         prob2_HFminus;
   Float_t         probge3_HFminus;
   Float_t         prob0_LFplus;
   Float_t         probge1_LFplus;
   Float_t         prob1_LFplus;
   Float_t         probge2_LFplus;
   Float_t         prob2_LFplus;
   Float_t         probge3_LFplus;
   Float_t         prob0_LFminus;
   Float_t         probge1_LFminus;
   Float_t         prob1_LFminus;
   Float_t         probge2_LFminus;
   Float_t         prob2_LFminus;
   Float_t         probge3_LFminus;
   Bool_t          cutHT;
   Bool_t          cutPV;
   Bool_t          cutTrigger;
   Bool_t          cut3Jets;
   Bool_t          cutEleVeto;
   Bool_t          cutMuVeto;
   Bool_t          cutMET;
   Bool_t          cutDeltaPhi;
   Bool_t          csctighthaloFilter;
   Bool_t          eenoiseFilter;
   Bool_t          greedymuonFilter;
   Bool_t          hbhenoiseFilter;
   Bool_t          inconsistentmuonFilter;
   Bool_t          ra2ecaltpFilter;
   Bool_t          scrapingvetoFilter;
   Bool_t          trackingfailureFilter;
   Bool_t          trackingfailureFilterPFLOW;
   Bool_t          recovrechitFilter;
   Bool_t          passCleaning;
   Int_t           PBNRcode;
   Bool_t          buggyEvent;
   Int_t           nGoodPV;
   Int_t           SUSY_nb;
   Int_t           SUSY_process;
   Float_t         SUSY_recoilPt;
   Int_t           njets;
   Int_t           njets30;
   Int_t           nbjets;
   Int_t           ntruebjets;
   Int_t           nElectrons;
   Int_t           nMuons;
   Int_t           nTaus;
   Int_t           nElectrons5;
   Int_t           nMuons5;
   Int_t           nElectrons15;
   Int_t           nMuons15;
   Int_t           nbjetsSSVM;
   Int_t           nbjetsSSVHPT;
   Int_t           nbjetsTCHPT;
   Int_t           nbjetsTCHET;
   Int_t           nbjetsTCHPM;
   Int_t           nbjetsCSVM;
   Int_t           nbjetsCSVL;
   Bool_t          isRealData;
   Bool_t          pass_utilityHLT;
   UInt_t          prescaleUtilityHLT;
   UInt_t          versionUtilityHLT;
   UInt_t          pass_utilityHLT_HT300_CentralJet30_BTagIP;
   UInt_t          prescale_utilityHLT_HT300_CentralJet30_BTagIP;
   Float_t         HT;
   Float_t         ST;
   Float_t         MET;
   Float_t         METphi;
   Float_t         MHT;
   Float_t         caloMET;
   Float_t         METPFType1;
   Float_t         METPFType1phi;
   Float_t         bestWMass;
   Float_t         bestTopMass;
   Float_t         topCosHel;
   Float_t         WCosHel;
   Float_t         MT_Wlep;
   Float_t         MT_Wlep5;
   Float_t         MT_Wlep15;
   Float_t         minDeltaPhi;
   Float_t         minDeltaPhiAll;
   Float_t         minDeltaPhiAll30;
   Float_t         minDeltaPhi30_eta5_noIdAll;
   Float_t         minDeltaPhiMetTau;
   Float_t         deltaPhi1;
   Float_t         deltaPhi2;
   Float_t         deltaPhi3;
   Float_t         maxDeltaPhi;
   Float_t         maxDeltaPhiAll;
   Float_t         maxDeltaPhiAll30;
   Float_t         maxDeltaPhi30_eta5_noIdAll;
   Float_t         sumDeltaPhi;
   Float_t         diffDeltaPhi;
   Float_t         deltaPhiStar;
   Float_t         deltaPhiStar_badjet_pt;
   Float_t         deltaPhiStar_badjet_phi;
   Float_t         deltaPhiStar_badjet_eta;
   Float_t         minDeltaPhiN;
   Float_t         deltaPhiN1;
   Float_t         deltaPhiN2;
   Float_t         deltaPhiN3;
   UInt_t          minDeltaPhi_chosenJet;
   UInt_t          minDeltaPhiN_chosenJet;
   UInt_t          maxJetMis_chosenJet;
   UInt_t          maxJetFracMis_chosenJet;
   Int_t           jetMisCategoryT1;
   Int_t           jetMisCategoryT1_noMETfraction;
   Int_t           jetMisCategoryT2;
   Int_t           jetMisCategoryT2_noMETfraction;
   Float_t         minDeltaPhiN_deltaT;
   Float_t         deltaT1;
   Float_t         deltaT2;
   Float_t         deltaT3;
   Float_t         minDeltaPhiN_MC;
   Float_t         deltaPhiN_MC1;
   Float_t         deltaPhiN_MC2;
   Float_t         deltaPhiN_MC3;
   Float_t         deltaT_MC1;
   Float_t         deltaT_MC2;
   Float_t         deltaT_MC3;
   Float_t         minDeltaPhiN_MCq;
   Float_t         deltaPhiN_MCq1;
   Float_t         deltaPhiN_MCq2;
   Float_t         deltaPhiN_MCq3;
   Float_t         deltaT_MCq1;
   Float_t         deltaT_MCq2;
   Float_t         deltaT_MCq3;
   Float_t         minDeltaPhiN_MCk;
   Float_t         deltaPhiN_MCk1;
   Float_t         deltaPhiN_MCk2;
   Float_t         deltaPhiN_MCk3;
   Float_t         deltaT_MCk1;
   Float_t         deltaT_MCk2;
   Float_t         deltaT_MCk3;
   Float_t         deltaT1_otherPt40;
   Float_t         deltaT1_otherPt50;
   Float_t         deltaT1_otherPt70;
   Float_t         deltaT2_otherPt40;
   Float_t         deltaT2_otherPt50;
   Float_t         deltaT2_otherPt70;
   Float_t         deltaT3_otherPt40;
   Float_t         deltaT3_otherPt50;
   Float_t         deltaT3_otherPt70;
   Float_t         minDeltaPhiNAll;
   Float_t         minDeltaPhiNAll30;
   Float_t         minDeltaPhiN_otherEta5;
   Float_t         minDeltaPhiN_otherEta5idNo;
   Float_t         minDeltaPhiN_mainPt30_otherEta5idNo;
   Float_t         minDeltaPhiN_mainPt30Eta5_otherEta5idNo;
   Float_t         minDeltaPhiN_mainPt30Eta5idNo_otherEta5idNo;
   Float_t         minDeltaPhiK;
   Float_t         minDeltaPhiKw;
   Float_t         minDeltaPhiK_otherEta5;
   Float_t         minDeltaPhiK_otherEta5idNo;
   Float_t         minDeltaPhiK_mainPt30_otherEta5idNo;
   Float_t         minDeltaPhiK_mainPt30Eta5_otherEta5idNo;
   Float_t         minDeltaPhiK_mainPt30Eta5idNo_otherEta5idNo;
   Float_t         minDeltaPhiN_DJR;
   Float_t         minDeltaPhiN_DJR_otherEta5;
   Float_t         minDeltaPhiK_DJR;
   Float_t         minDeltaPhiK_DJR_otherEta5;
   Float_t         minDeltaPhiN_otherPt10;
   Float_t         minDeltaPhiN_otherPt20;
   Float_t         minDeltaPhiN_otherPt30;
   Float_t         minDeltaPhiN_otherPt40;
   Float_t         minDeltaPhiN_otherPt50;
   Float_t         minDeltaPhiN_DJR_otherPt10;
   Float_t         minDeltaPhiN_DJR_otherPt20;
   Float_t         minDeltaPhiN_DJR_otherPt30;
   Float_t         minDeltaPhiN_DJR_otherPt40;
   Float_t         minDeltaPhiN_DJR_otherPt50;
   Float_t         minDeltaPhiN_withLepton;
   Float_t         maxJetMis;
   Float_t         max2JetMis;
   Float_t         maxJetMisAll30;
   Float_t         max2JetMisAll30;
   Float_t         maxJetMis_signed;
   Float_t         maxJetMis30_signed;
   Float_t         maxJetFracMis;
   Float_t         max2JetFracMis;
   Float_t         maxJetFracMisAll30;
   Float_t         max2JetFracMisAll30;
   Float_t         deltaPhiMETJetMaxMis;
   Float_t         deltaPhiMETJetMaxMis30;
   UInt_t          rankMaxJetMis;
   UInt_t          rankMaxJetMis30;
   Float_t         CSVout1;
   Float_t         CSVout2;
   Float_t         CSVout3;
   Float_t         minDeltaPhiAllb30;
   Float_t         deltaPhib1;
   Float_t         deltaPhib2;
   Float_t         deltaPhib3;
   Float_t         minDeltaPhiMETMuonsAll;
   Float_t         minDeltaPhiN_lostJet;
   Float_t         deltaPhiN1_lostJet;
   Float_t         deltaPhiN2_lostJet;
   Float_t         deltaPhiN3_lostJet;
   Bool_t          passBadECAL_METphi53_n10_s12;
   Float_t         worstMisA_badECAL_METphi53_n10_s12;
   Float_t         worstMisF_badECAL_METphi53_n10_s12;
   Bool_t          passBadECAL_METphi33_n10_s12;
   Float_t         worstMisA_badECAL_METphi33_n10_s12;
   Float_t         worstMisF_badECAL_METphi33_n10_s12;
   Bool_t          passBadECAL_METphi52_n10_s12;
   Float_t         worstMisA_badECAL_METphi52_n10_s12;
   Float_t         worstMisF_badECAL_METphi52_n10_s12;
   Bool_t          passBadECAL_METphi53_n5_s12;
   Bool_t          passBadECAL_METphi53_crack;
   Float_t         worstMisA_badECAL_METphi53_crack;
   Float_t         worstMisF_badECAL_METphi53_crack;
   Bool_t          passBadECAL_METphi52_crack;
   Float_t         worstMisA_badECAL_METphi52_crack;
   Float_t         worstMisF_badECAL_METphi52_crack;
   Bool_t          passBadECAL_METphi53_crackF;
   Float_t         worstMisA_badECAL_METphi53_crackF;
   Float_t         worstMisF_badECAL_METphi53_crackF;
   Bool_t          passBadECAL_METphi52_crackF;
   Float_t         worstMisA_badECAL_METphi52_crackF;
   Float_t         worstMisF_badECAL_METphi52_crackF;
   Float_t         jetpt1;
   Float_t         jetgenpt1;
   Float_t         jeteta1;
   Float_t         jetgeneta1;
   Float_t         jetphi1;
   Float_t         jetgenphi1;
   Float_t         jetenergy1;
   Int_t           jetflavor1;
   Float_t         jetchargedhadronfrac1;
   Int_t           jetchargedhadronmult1;
   Float_t         jetpt2;
   Float_t         jetgenpt2;
   Float_t         jeteta2;
   Float_t         jetgeneta2;
   Float_t         jetphi2;
   Float_t         jetgenphi2;
   Float_t         jetenergy2;
   Int_t           jetflavor2;
   Float_t         jetchargedhadronfrac2;
   Int_t           jetchargedhadronmult2;
   Float_t         jetpt3;
   Float_t         jetgenpt3;
   Float_t         jeteta3;
   Float_t         jetgeneta3;
   Float_t         jetphi3;
   Float_t         jetgenphi3;
   Float_t         jetenergy3;
   Int_t           jetflavor3;
   Float_t         jetchargedhadronfrac3;
   Int_t           jetchargedhadronmult3;
   Float_t         bjetpt1;
   Float_t         bjeteta1;
   Float_t         bjetphi1;
   Float_t         bjetenergy1;
   Int_t           bjetflavor1;
   Float_t         bjetchargedhadronfrac1;
   Int_t           bjetchargedhadronmult1;
   Float_t         bjetpt2;
   Float_t         bjeteta2;
   Float_t         bjetphi2;
   Float_t         bjetenergy2;
   Int_t           bjetflavor2;
   Float_t         bjetchargedhadronfrac2;
   Int_t           bjetchargedhadronmult2;
   Float_t         bjetpt3;
   Float_t         bjeteta3;
   Float_t         bjetphi3;
   Float_t         bjetenergy3;
   Int_t           bjetflavor3;
   Float_t         bjetchargedhadronfrac3;
   Int_t           bjetchargedhadronmult3;
   Float_t         eleet1;
   Int_t           elecharge1;
   Float_t         elephi1;
   Float_t         eleeta1;
   Float_t         muonpt1;
   Int_t           muoncharge1;
   Float_t         muonphi1;
   Float_t         muoneta1;
   Float_t         muoniso1;
   Float_t         muonchhadiso1;
   Float_t         muonphotoniso1;
   Float_t         muonneutralhadiso1;
   Float_t         eleet2;
   Float_t         elephi2;
   Float_t         eleeta2;
   Float_t         muonpt2;
   Float_t         muonphi2;
   Float_t         muoneta2;
   Float_t         taupt1;
   Float_t         taueta1;
   Float_t         eleRelIso;
   Float_t         muonRelIso;
   Float_t         recomuonpt1;
   Float_t         recomuonphi1;
   Float_t         recomuoneta1;
   Float_t         recomuoniso1;
   Float_t         recomuonmindphijet1;
   Float_t         recomuonpt2;
   Float_t         recomuonphi2;
   Float_t         recomuoneta2;
   Float_t         recomuoniso2;
   Float_t         recomuonmindphijet2;
   Float_t         minDeltaPhiN_Luke;
   Float_t         maxDeltaPhiN_Luke;
   Float_t         deltaPhiN1_Luke;
   Float_t         deltaPhiN2_Luke;
   Float_t         deltaPhiN3_Luke;
   Float_t         minTransverseMETSignificance;
   Float_t         maxTransverseMETSignificance;
   Float_t         transverseMETSignificance1;
   Float_t         transverseMETSignificance2;
   Float_t         transverseMETSignificance3;
   Int_t           njets_lostJet;
   Int_t           nbjets_lostJet;
   Float_t         minDeltaPhiN_Luke_lostJet;
   Float_t         maxDeltaPhiN_Luke_lostJet;
   Float_t         deltaPhiN1_Luke_lostJet;
   Float_t         deltaPhiN2_Luke_lostJet;
   Float_t         deltaPhiN3_Luke_lostJet;
   Float_t         minTransverseMETSignificance_lostJet;
   Float_t         maxTransverseMETSignificance_lostJet;
   Float_t         transverseMETSignificance1_lostJet;
   Float_t         transverseMETSignificance2_lostJet;
   Float_t         transverseMETSignificance3_lostJet;
   Float_t         nLostJet;

   // List of branches
   TBranch        *b_weight;   //!
   TBranch        *b_scanCrossSection;   //!
   TBranch        *b_scanCrossSectionPlus;   //!
   TBranch        *b_scanCrossSectionMinus;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_lumiSection;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_m0;   //!
   TBranch        *b_m12;   //!
   TBranch        *b_W1decayType;   //!
   TBranch        *b_W2decayType;   //!
   TBranch        *b_decayType;   //!
   TBranch        *b_btagIPweight;   //!
   TBranch        *b_PUweight;   //!
   TBranch        *b_hltHTeff;   //!
   TBranch        *b_hltMHTeff;   //!
   TBranch        *b_hltMHTeffBNN;   //!
   TBranch        *b_hltMHTeffBNNUp;   //!
   TBranch        *b_hltMHTeffBNNDown;   //!
   TBranch        *b_pdfWeightsCTEQ;   //!
   TBranch        *b_pdfWeightsMSTW;   //!
   TBranch        *b_pdfWeightsNNPDF;   //!
   TBranch        *b_prob0;   //!
   TBranch        *b_probge1;   //!
   TBranch        *b_prob1;   //!
   TBranch        *b_probge2;   //!
   TBranch        *b_prob2;   //!
   TBranch        *b_probge3;   //!
   TBranch        *b_prob0_HFplus;   //!
   TBranch        *b_probge1_HFplus;   //!
   TBranch        *b_prob1_HFplus;   //!
   TBranch        *b_probge2_HFplus;   //!
   TBranch        *b_prob2_HFplus;   //!
   TBranch        *b_probge3_HFplus;   //!
   TBranch        *b_prob0_HFminus;   //!
   TBranch        *b_probge1_HFminus;   //!
   TBranch        *b_prob1_HFminus;   //!
   TBranch        *b_probge2_HFminus;   //!
   TBranch        *b_prob2_HFminus;   //!
   TBranch        *b_probge3_HFminus;   //!
   TBranch        *b_prob0_LFplus;   //!
   TBranch        *b_probge1_LFplus;   //!
   TBranch        *b_prob1_LFplus;   //!
   TBranch        *b_probge2_LFplus;   //!
   TBranch        *b_prob2_LFplus;   //!
   TBranch        *b_probge3_LFplus;   //!
   TBranch        *b_prob0_LFminus;   //!
   TBranch        *b_probge1_LFminus;   //!
   TBranch        *b_prob1_LFminus;   //!
   TBranch        *b_probge2_LFminus;   //!
   TBranch        *b_prob2_LFminus;   //!
   TBranch        *b_probge3_LFminus;   //!
   TBranch        *b_cutHT;   //!
   TBranch        *b_cutPV;   //!
   TBranch        *b_cutTrigger;   //!
   TBranch        *b_cut3Jets;   //!
   TBranch        *b_cutEleVeto;   //!
   TBranch        *b_cutMuVeto;   //!
   TBranch        *b_cutMET;   //!
   TBranch        *b_cutDeltaPhi;   //!
   TBranch        *b_csctighthaloFilter;   //!
   TBranch        *b_eenoiseFilter;   //!
   TBranch        *b_greedymuonFilter;   //!
   TBranch        *b_hbhenoiseFilter;   //!
   TBranch        *b_inconsistentmuonFilter;   //!
   TBranch        *b_ra2ecaltpFilter;   //!
   TBranch        *b_scrapingvetoFilter;   //!
   TBranch        *b_trackingfailureFilter;   //!
   TBranch        *b_trackingfailureFilterPFLOW;   //!
   TBranch        *b_recovrechitFilter;   //!
   TBranch        *b_passCleaning;   //!
   TBranch        *b_PBNRcode;   //!
   TBranch        *b_buggyEvent;   //!
   TBranch        *b_nGoodPV;   //!
   TBranch        *b_SUSY_nb;   //!
   TBranch        *b_SUSY_process;   //!
   TBranch        *b_SUSY_recoilPt;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_njets30;   //!
   TBranch        *b_nbjets;   //!
   TBranch        *b_ntruebjets;   //!
   TBranch        *b_nElectrons;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_nTaus;   //!
   TBranch        *b_nElectrons5;   //!
   TBranch        *b_nMuons5;   //!
   TBranch        *b_nElectrons15;   //!
   TBranch        *b_nMuons15;   //!
   TBranch        *b_nbjetsSSVM;   //!
   TBranch        *b_nbjetsSSVHPT;   //!
   TBranch        *b_nbjetsTCHPT;   //!
   TBranch        *b_nbjetsTCHET;   //!
   TBranch        *b_nbjetsTCHPM;   //!
   TBranch        *b_nbjetsCSVM;   //!
   TBranch        *b_nbjetsCSVL;   //!
   TBranch        *b_isRealData;   //!
   TBranch        *b_pass_utilityHLT;   //!
   TBranch        *b_prescaleUtilityHLT;   //!
   TBranch        *b_versionUtilityHLT;   //!
   TBranch        *b_pass_utilityHLT_HT300_CentralJet30_BTagIP;   //!
   TBranch        *b_prescale_utilityHLT_HT300_CentralJet30_BTagIP;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_ST;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METphi;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_caloMET;   //!
   TBranch        *b_METPFType1;   //!
   TBranch        *b_METPFType1phi;   //!
   TBranch        *b_bestWMass;   //!
   TBranch        *b_bestTopMass;   //!
   TBranch        *b_topCosHel;   //!
   TBranch        *b_WCosHel;   //!
   TBranch        *b_MT_Wlep;   //!
   TBranch        *b_MT_Wlep5;   //!
   TBranch        *b_MT_Wlep15;   //!
   TBranch        *b_minDeltaPhi;   //!
   TBranch        *b_minDeltaPhiAll;   //!
   TBranch        *b_minDeltaPhiAll30;   //!
   TBranch        *b_minDeltaPhi30_eta5_noIdAll;   //!
   TBranch        *b_minDeltaPhiMetTau;   //!
   TBranch        *b_deltaPhi1;   //!
   TBranch        *b_deltaPhi2;   //!
   TBranch        *b_deltaPhi3;   //!
   TBranch        *b_maxDeltaPhi;   //!
   TBranch        *b_maxDeltaPhiAll;   //!
   TBranch        *b_maxDeltaPhiAll30;   //!
   TBranch        *b_maxDeltaPhi30_eta5_noIdAll;   //!
   TBranch        *b_sumDeltaPhi;   //!
   TBranch        *b_diffDeltaPhi;   //!
   TBranch        *b_deltaPhiStar;   //!
   TBranch        *b_deltaPhiStar_badjet_pt;   //!
   TBranch        *b_deltaPhiStar_badjet_phi;   //!
   TBranch        *b_deltaPhiStar_badjet_eta;   //!
   TBranch        *b_minDeltaPhiN;   //!
   TBranch        *b_deltaPhiN1;   //!
   TBranch        *b_deltaPhiN2;   //!
   TBranch        *b_deltaPhiN3;   //!
   TBranch        *b_minDeltaPhi_chosenJet;   //!
   TBranch        *b_minDeltaPhiN_chosenJet;   //!
   TBranch        *b_maxJetMis_chosenJet;   //!
   TBranch        *b_maxJetFracMis_chosenJet;   //!
   TBranch        *b_jetMisCategoryT1;   //!
   TBranch        *b_jetMisCategoryT1_noMETfraction;   //!
   TBranch        *b_jetMisCategoryT2;   //!
   TBranch        *b_jetMisCategoryT2_noMETfraction;   //!
   TBranch        *b_minDeltaPhiN_deltaT;   //!
   TBranch        *b_deltaT1;   //!
   TBranch        *b_deltaT2;   //!
   TBranch        *b_deltaT3;   //!
   TBranch        *b_minDeltaPhiN_MC;   //!
   TBranch        *b_deltaPhiN_MC1;   //!
   TBranch        *b_deltaPhiN_MC2;   //!
   TBranch        *b_deltaPhiN_MC3;   //!
   TBranch        *b_deltaT_MC1;   //!
   TBranch        *b_deltaT_MC2;   //!
   TBranch        *b_deltaT_MC3;   //!
   TBranch        *b_minDeltaPhiN_MCq;   //!
   TBranch        *b_deltaPhiN_MCq1;   //!
   TBranch        *b_deltaPhiN_MCq2;   //!
   TBranch        *b_deltaPhiN_MCq3;   //!
   TBranch        *b_deltaT_MCq1;   //!
   TBranch        *b_deltaT_MCq2;   //!
   TBranch        *b_deltaT_MCq3;   //!
   TBranch        *b_minDeltaPhiN_MCk;   //!
   TBranch        *b_deltaPhiN_MCk1;   //!
   TBranch        *b_deltaPhiN_MCk2;   //!
   TBranch        *b_deltaPhiN_MCk3;   //!
   TBranch        *b_deltaT_MCk1;   //!
   TBranch        *b_deltaT_MCk2;   //!
   TBranch        *b_deltaT_MCk3;   //!
   TBranch        *b_deltaT1_otherPt40;   //!
   TBranch        *b_deltaT1_otherPt50;   //!
   TBranch        *b_deltaT1_otherPt70;   //!
   TBranch        *b_deltaT2_otherPt40;   //!
   TBranch        *b_deltaT2_otherPt50;   //!
   TBranch        *b_deltaT2_otherPt70;   //!
   TBranch        *b_deltaT3_otherPt40;   //!
   TBranch        *b_deltaT3_otherPt50;   //!
   TBranch        *b_deltaT3_otherPt70;   //!
   TBranch        *b_minDeltaPhiNAll;   //!
   TBranch        *b_minDeltaPhiNAll30;   //!
   TBranch        *b_minDeltaPhiN_otherEta5;   //!
   TBranch        *b_minDeltaPhiN_otherEta5idNo;   //!
   TBranch        *b_minDeltaPhiN_mainPt30_otherEta5idNo;   //!
   TBranch        *b_minDeltaPhiN_mainPt30Eta5_otherEta5idNo;   //!
   TBranch        *b_minDeltaPhiN_mainPt30Eta5idNo_otherEta5idNo;   //!
   TBranch        *b_minDeltaPhiK;   //!
   TBranch        *b_minDeltaPhiKw;   //!
   TBranch        *b_minDeltaPhiK_otherEta5;   //!
   TBranch        *b_minDeltaPhiK_otherEta5idNo;   //!
   TBranch        *b_minDeltaPhiK_mainPt30_otherEta5idNo;   //!
   TBranch        *b_minDeltaPhiK_mainPt30Eta5_otherEta5idNo;   //!
   TBranch        *b_minDeltaPhiK_mainPt30Eta5idNo_otherEta5idNo;   //!
   TBranch        *b_minDeltaPhiN_DJR;   //!
   TBranch        *b_minDeltaPhiN_DJR_otherEta5;   //!
   TBranch        *b_minDeltaPhiK_DJR;   //!
   TBranch        *b_minDeltaPhiK_DJR_otherEta5;   //!
   TBranch        *b_minDeltaPhiN_otherPt10;   //!
   TBranch        *b_minDeltaPhiN_otherPt20;   //!
   TBranch        *b_minDeltaPhiN_otherPt30;   //!
   TBranch        *b_minDeltaPhiN_otherPt40;   //!
   TBranch        *b_minDeltaPhiN_otherPt50;   //!
   TBranch        *b_minDeltaPhiN_DJR_otherPt10;   //!
   TBranch        *b_minDeltaPhiN_DJR_otherPt20;   //!
   TBranch        *b_minDeltaPhiN_DJR_otherPt30;   //!
   TBranch        *b_minDeltaPhiN_DJR_otherPt40;   //!
   TBranch        *b_minDeltaPhiN_DJR_otherPt50;   //!
   TBranch        *b_minDeltaPhiN_withLepton;   //!
   TBranch        *b_maxJetMis;   //!
   TBranch        *b_max2JetMis;   //!
   TBranch        *b_maxJetMisAll30;   //!
   TBranch        *b_max2JetMisAll30;   //!
   TBranch        *b_maxJetMis_signed;   //!
   TBranch        *b_maxJetMis30_signed;   //!
   TBranch        *b_maxJetFracMis;   //!
   TBranch        *b_max2JetFracMis;   //!
   TBranch        *b_maxJetFracMisAll30;   //!
   TBranch        *b_max2JetFracMisAll30;   //!
   TBranch        *b_deltaPhiMETJetMaxMis;   //!
   TBranch        *b_deltaPhiMETJetMaxMis30;   //!
   TBranch        *b_rankMaxJetMis;   //!
   TBranch        *b_rankMaxJetMis30;   //!
   TBranch        *b_CSVout1;   //!
   TBranch        *b_CSVout2;   //!
   TBranch        *b_CSVout3;   //!
   TBranch        *b_minDeltaPhiAllb30;   //!
   TBranch        *b_deltaPhib1;   //!
   TBranch        *b_deltaPhib2;   //!
   TBranch        *b_deltaPhib3;   //!
   TBranch        *b_minDeltaPhiMETMuonsAll;   //!
   TBranch        *b_minDeltaPhiN_lostJet;   //!
   TBranch        *b_deltaPhiN1_lostJet;   //!
   TBranch        *b_deltaPhiN2_lostJet;   //!
   TBranch        *b_deltaPhiN3_lostJet;   //!
   TBranch        *b_passBadECAL_METphi53_n10_s12;   //!
   TBranch        *b_worstMisA_badECAL_METphi53_n10_s12;   //!
   TBranch        *b_worstMisF_badECAL_METphi53_n10_s12;   //!
   TBranch        *b_passBadECAL_METphi33_n10_s12;   //!
   TBranch        *b_worstMisA_badECAL_METphi33_n10_s12;   //!
   TBranch        *b_worstMisF_badECAL_METphi33_n10_s12;   //!
   TBranch        *b_passBadECAL_METphi52_n10_s12;   //!
   TBranch        *b_worstMisA_badECAL_METphi52_n10_s12;   //!
   TBranch        *b_worstMisF_badECAL_METphi52_n10_s12;   //!
   TBranch        *b_passBadECAL_METphi53_n5_s12;   //!
   TBranch        *b_passBadECAL_METphi53_crack;   //!
   TBranch        *b_worstMisA_badECAL_METphi53_crack;   //!
   TBranch        *b_worstMisF_badECAL_METphi53_crack;   //!
   TBranch        *b_passBadECAL_METphi52_crack;   //!
   TBranch        *b_worstMisA_badECAL_METphi52_crack;   //!
   TBranch        *b_worstMisF_badECAL_METphi52_crack;   //!
   TBranch        *b_passBadECAL_METphi53_crackF;   //!
   TBranch        *b_worstMisA_badECAL_METphi53_crackF;   //!
   TBranch        *b_worstMisF_badECAL_METphi53_crackF;   //!
   TBranch        *b_passBadECAL_METphi52_crackF;   //!
   TBranch        *b_worstMisA_badECAL_METphi52_crackF;   //!
   TBranch        *b_worstMisF_badECAL_METphi52_crackF;   //!
   TBranch        *b_jetpt1;   //!
   TBranch        *b_jetgenpt1;   //!
   TBranch        *b_jeteta1;   //!
   TBranch        *b_jetgeneta1;   //!
   TBranch        *b_jetphi1;   //!
   TBranch        *b_jetgenphi1;   //!
   TBranch        *b_jetenergy1;   //!
   TBranch        *b_jetflavor1;   //!
   TBranch        *b_jetchargedhadronfrac1;   //!
   TBranch        *b_jetchargedhadronmult1;   //!
   TBranch        *b_jetpt2;   //!
   TBranch        *b_jetgenpt2;   //!
   TBranch        *b_jeteta2;   //!
   TBranch        *b_jetgeneta2;   //!
   TBranch        *b_jetphi2;   //!
   TBranch        *b_jetgenphi2;   //!
   TBranch        *b_jetenergy2;   //!
   TBranch        *b_jetflavor2;   //!
   TBranch        *b_jetchargedhadronfrac2;   //!
   TBranch        *b_jetchargedhadronmult2;   //!
   TBranch        *b_jetpt3;   //!
   TBranch        *b_jetgenpt3;   //!
   TBranch        *b_jeteta3;   //!
   TBranch        *b_jetgeneta3;   //!
   TBranch        *b_jetphi3;   //!
   TBranch        *b_jetgenphi3;   //!
   TBranch        *b_jetenergy3;   //!
   TBranch        *b_jetflavor3;   //!
   TBranch        *b_jetchargedhadronfrac3;   //!
   TBranch        *b_jetchargedhadronmult3;   //!
   TBranch        *b_bjetpt1;   //!
   TBranch        *b_bjeteta1;   //!
   TBranch        *b_bjetphi1;   //!
   TBranch        *b_bjetenergy1;   //!
   TBranch        *b_bjetflavor1;   //!
   TBranch        *b_bjetchargedhadronfrac1;   //!
   TBranch        *b_bjetchargedhadronmult1;   //!
   TBranch        *b_bjetpt2;   //!
   TBranch        *b_bjeteta2;   //!
   TBranch        *b_bjetphi2;   //!
   TBranch        *b_bjetenergy2;   //!
   TBranch        *b_bjetflavor2;   //!
   TBranch        *b_bjetchargedhadronfrac2;   //!
   TBranch        *b_bjetchargedhadronmult2;   //!
   TBranch        *b_bjetpt3;   //!
   TBranch        *b_bjeteta3;   //!
   TBranch        *b_bjetphi3;   //!
   TBranch        *b_bjetenergy3;   //!
   TBranch        *b_bjetflavor3;   //!
   TBranch        *b_bjetchargedhadronfrac3;   //!
   TBranch        *b_bjetchargedhadronmult3;   //!
   TBranch        *b_eleet1;   //!
   TBranch        *b_elecharge1;   //!
   TBranch        *b_elephi1;   //!
   TBranch        *b_eleeta1;   //!
   TBranch        *b_muonpt1;   //!
   TBranch        *b_muoncharge1;   //!
   TBranch        *b_muonphi1;   //!
   TBranch        *b_muoneta1;   //!
   TBranch        *b_muoniso1;   //!
   TBranch        *b_muonchhadiso1;   //!
   TBranch        *b_muonphotoniso1;   //!
   TBranch        *b_muonneutralhadiso1;   //!
   TBranch        *b_eleet2;   //!
   TBranch        *b_elephi2;   //!
   TBranch        *b_eleeta2;   //!
   TBranch        *b_muonpt2;   //!
   TBranch        *b_muonphi2;   //!
   TBranch        *b_muoneta2;   //!
   TBranch        *b_taupt1;   //!
   TBranch        *b_taueta1;   //!
   TBranch        *b_eleRelIso;   //!
   TBranch        *b_muonRelIso;   //!
   TBranch        *b_recomuonpt1;   //!
   TBranch        *b_recomuonphi1;   //!
   TBranch        *b_recomuoneta1;   //!
   TBranch        *b_recomuoniso1;   //!
   TBranch        *b_recomuonmindphijet1;   //!
   TBranch        *b_recomuonpt2;   //!
   TBranch        *b_recomuonphi2;   //!
   TBranch        *b_recomuoneta2;   //!
   TBranch        *b_recomuoniso2;   //!
   TBranch        *b_recomuonmindphijet2;   //!
   TBranch        *b_minDeltaPhiN_Luke;   //!
   TBranch        *b_maxDeltaPhiN_Luke;   //!
   TBranch        *b_deltaPhiN1_Luke;   //!
   TBranch        *b_deltaPhiN2_Luke;   //!
   TBranch        *b_deltaPhiN3_Luke;   //!
   TBranch        *b_minTransverseMETSignificance;   //!
   TBranch        *b_maxTransverseMETSignificance;   //!
   TBranch        *b_transverseMETSignificance1;   //!
   TBranch        *b_transverseMETSignificance2;   //!
   TBranch        *b_transverseMETSignificance3;   //!
   TBranch        *b_njets_lostJet;   //!
   TBranch        *b_nbjets_lostJet;   //!
   TBranch        *b_minDeltaPhiN_Luke_lostJet;   //!
   TBranch        *b_maxDeltaPhiN_Luke_lostJet;   //!
   TBranch        *b_deltaPhiN1_Luke_lostJet;   //!
   TBranch        *b_deltaPhiN2_Luke_lostJet;   //!
   TBranch        *b_deltaPhiN3_Luke_lostJet;   //!
   TBranch        *b_minTransverseMETSignificance_lostJet;   //!
   TBranch        *b_maxTransverseMETSignificance_lostJet;   //!
   TBranch        *b_transverseMETSignificance1_lostJet;   //!
   TBranch        *b_transverseMETSignificance2_lostJet;   //!
   TBranch        *b_transverseMETSignificance3_lostJet;   //!
   TBranch        *b_nLostJet;   //!

   ISRweights(TString sample);
   virtual ~ISRweights();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString selection);
   void cleanup3B();
   void cleanup2BT();
   void cleanup1BT();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void loadIsrHistos(TString sample) ;
   double getIsrWeight(int m0,int m12, double isrpt);

   std::map< pair<int,int>, TH1F*> isrHistos_;
   TFile* fISRweights_;

   TH2D* scanSMSngen_;
   TString sample_;
};

#endif

#ifdef ISRweights_cxx
ISRweights::ISRweights(TString sample) : fChain(0), fISRweights_(0) , scanSMSngen_(0), sample_(sample)
{
  TTree* tree=0;

  TString infile = "/cu3/joshmt/reducedTrees/V00-02-35v/reducedTree.CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.";
  infile += sample;
  infile +=".root";

// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(infile);
      if (!f || !f->IsOpen()) {
         f = new TFile(infile);
      }
      f->GetObject("reducedTree",tree);
      scanSMSngen_ = (TH2D*) f->Get("scanSMSngen");
   }
   Init(tree);

   loadIsrHistos(sample);
}

ISRweights::~ISRweights()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ISRweights::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ISRweights::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ISRweights::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("scanCrossSection", &scanCrossSection, &b_scanCrossSection);
   fChain->SetBranchAddress("scanCrossSectionPlus", &scanCrossSectionPlus, &b_scanCrossSectionPlus);
   fChain->SetBranchAddress("scanCrossSectionMinus", &scanCrossSectionMinus, &b_scanCrossSectionMinus);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("m0", &m0, &b_m0);
   fChain->SetBranchAddress("m12", &m12, &b_m12);
   fChain->SetBranchAddress("W1decayType", &W1decayType, &b_W1decayType);
   fChain->SetBranchAddress("W2decayType", &W2decayType, &b_W2decayType);
   fChain->SetBranchAddress("decayType", &decayType, &b_decayType);
   fChain->SetBranchAddress("btagIPweight", &btagIPweight, &b_btagIPweight);
   fChain->SetBranchAddress("PUweight", &PUweight, &b_PUweight);
   fChain->SetBranchAddress("hltHTeff", &hltHTeff, &b_hltHTeff);
   fChain->SetBranchAddress("hltMHTeff", &hltMHTeff, &b_hltMHTeff);
   fChain->SetBranchAddress("hltMHTeffBNN", &hltMHTeffBNN, &b_hltMHTeffBNN);
   fChain->SetBranchAddress("hltMHTeffBNNUp", &hltMHTeffBNNUp, &b_hltMHTeffBNNUp);
   fChain->SetBranchAddress("hltMHTeffBNNDown", &hltMHTeffBNNDown, &b_hltMHTeffBNNDown);
   fChain->SetBranchAddress("pdfWeightsCTEQ", pdfWeightsCTEQ, &b_pdfWeightsCTEQ);
   fChain->SetBranchAddress("pdfWeightsMSTW", pdfWeightsMSTW, &b_pdfWeightsMSTW);
   fChain->SetBranchAddress("pdfWeightsNNPDF", pdfWeightsNNPDF, &b_pdfWeightsNNPDF);
   fChain->SetBranchAddress("prob0", &prob0, &b_prob0);
   fChain->SetBranchAddress("probge1", &probge1, &b_probge1);
   fChain->SetBranchAddress("prob1", &prob1, &b_prob1);
   fChain->SetBranchAddress("probge2", &probge2, &b_probge2);
   fChain->SetBranchAddress("prob2", &prob2, &b_prob2);
   fChain->SetBranchAddress("probge3", &probge3, &b_probge3);
   fChain->SetBranchAddress("prob0_HFplus", &prob0_HFplus, &b_prob0_HFplus);
   fChain->SetBranchAddress("probge1_HFplus", &probge1_HFplus, &b_probge1_HFplus);
   fChain->SetBranchAddress("prob1_HFplus", &prob1_HFplus, &b_prob1_HFplus);
   fChain->SetBranchAddress("probge2_HFplus", &probge2_HFplus, &b_probge2_HFplus);
   fChain->SetBranchAddress("prob2_HFplus", &prob2_HFplus, &b_prob2_HFplus);
   fChain->SetBranchAddress("probge3_HFplus", &probge3_HFplus, &b_probge3_HFplus);
   fChain->SetBranchAddress("prob0_HFminus", &prob0_HFminus, &b_prob0_HFminus);
   fChain->SetBranchAddress("probge1_HFminus", &probge1_HFminus, &b_probge1_HFminus);
   fChain->SetBranchAddress("prob1_HFminus", &prob1_HFminus, &b_prob1_HFminus);
   fChain->SetBranchAddress("probge2_HFminus", &probge2_HFminus, &b_probge2_HFminus);
   fChain->SetBranchAddress("prob2_HFminus", &prob2_HFminus, &b_prob2_HFminus);
   fChain->SetBranchAddress("probge3_HFminus", &probge3_HFminus, &b_probge3_HFminus);
   fChain->SetBranchAddress("prob0_LFplus", &prob0_LFplus, &b_prob0_LFplus);
   fChain->SetBranchAddress("probge1_LFplus", &probge1_LFplus, &b_probge1_LFplus);
   fChain->SetBranchAddress("prob1_LFplus", &prob1_LFplus, &b_prob1_LFplus);
   fChain->SetBranchAddress("probge2_LFplus", &probge2_LFplus, &b_probge2_LFplus);
   fChain->SetBranchAddress("prob2_LFplus", &prob2_LFplus, &b_prob2_LFplus);
   fChain->SetBranchAddress("probge3_LFplus", &probge3_LFplus, &b_probge3_LFplus);
   fChain->SetBranchAddress("prob0_LFminus", &prob0_LFminus, &b_prob0_LFminus);
   fChain->SetBranchAddress("probge1_LFminus", &probge1_LFminus, &b_probge1_LFminus);
   fChain->SetBranchAddress("prob1_LFminus", &prob1_LFminus, &b_prob1_LFminus);
   fChain->SetBranchAddress("probge2_LFminus", &probge2_LFminus, &b_probge2_LFminus);
   fChain->SetBranchAddress("prob2_LFminus", &prob2_LFminus, &b_prob2_LFminus);
   fChain->SetBranchAddress("probge3_LFminus", &probge3_LFminus, &b_probge3_LFminus);
   fChain->SetBranchAddress("cutHT", &cutHT, &b_cutHT);
   fChain->SetBranchAddress("cutPV", &cutPV, &b_cutPV);
   fChain->SetBranchAddress("cutTrigger", &cutTrigger, &b_cutTrigger);
   fChain->SetBranchAddress("cut3Jets", &cut3Jets, &b_cut3Jets);
   fChain->SetBranchAddress("cutEleVeto", &cutEleVeto, &b_cutEleVeto);
   fChain->SetBranchAddress("cutMuVeto", &cutMuVeto, &b_cutMuVeto);
   fChain->SetBranchAddress("cutMET", &cutMET, &b_cutMET);
   fChain->SetBranchAddress("cutDeltaPhi", &cutDeltaPhi, &b_cutDeltaPhi);
   fChain->SetBranchAddress("csctighthaloFilter", &csctighthaloFilter, &b_csctighthaloFilter);
   fChain->SetBranchAddress("eenoiseFilter", &eenoiseFilter, &b_eenoiseFilter);
   fChain->SetBranchAddress("greedymuonFilter", &greedymuonFilter, &b_greedymuonFilter);
   fChain->SetBranchAddress("hbhenoiseFilter", &hbhenoiseFilter, &b_hbhenoiseFilter);
   fChain->SetBranchAddress("inconsistentmuonFilter", &inconsistentmuonFilter, &b_inconsistentmuonFilter);
   fChain->SetBranchAddress("ra2ecaltpFilter", &ra2ecaltpFilter, &b_ra2ecaltpFilter);
   fChain->SetBranchAddress("scrapingvetoFilter", &scrapingvetoFilter, &b_scrapingvetoFilter);
   fChain->SetBranchAddress("trackingfailureFilter", &trackingfailureFilter, &b_trackingfailureFilter);
   fChain->SetBranchAddress("trackingfailureFilterPFLOW", &trackingfailureFilterPFLOW, &b_trackingfailureFilterPFLOW);
   fChain->SetBranchAddress("recovrechitFilter", &recovrechitFilter, &b_recovrechitFilter);
   fChain->SetBranchAddress("passCleaning", &passCleaning, &b_passCleaning);
   fChain->SetBranchAddress("PBNRcode", &PBNRcode, &b_PBNRcode);
   fChain->SetBranchAddress("buggyEvent", &buggyEvent, &b_buggyEvent);
   fChain->SetBranchAddress("nGoodPV", &nGoodPV, &b_nGoodPV);
   fChain->SetBranchAddress("SUSY_nb", &SUSY_nb, &b_SUSY_nb);
   fChain->SetBranchAddress("SUSY_process", &SUSY_process, &b_SUSY_process);
   fChain->SetBranchAddress("SUSY_recoilPt", &SUSY_recoilPt, &b_SUSY_recoilPt);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("njets30", &njets30, &b_njets30);
   fChain->SetBranchAddress("nbjets", &nbjets, &b_nbjets);
   fChain->SetBranchAddress("ntruebjets", &ntruebjets, &b_ntruebjets);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   fChain->SetBranchAddress("nTaus", &nTaus, &b_nTaus);
   fChain->SetBranchAddress("nElectrons5", &nElectrons5, &b_nElectrons5);
   fChain->SetBranchAddress("nMuons5", &nMuons5, &b_nMuons5);
   fChain->SetBranchAddress("nElectrons15", &nElectrons15, &b_nElectrons15);
   fChain->SetBranchAddress("nMuons15", &nMuons15, &b_nMuons15);
   fChain->SetBranchAddress("nbjetsSSVM", &nbjetsSSVM, &b_nbjetsSSVM);
   fChain->SetBranchAddress("nbjetsSSVHPT", &nbjetsSSVHPT, &b_nbjetsSSVHPT);
   fChain->SetBranchAddress("nbjetsTCHPT", &nbjetsTCHPT, &b_nbjetsTCHPT);
   fChain->SetBranchAddress("nbjetsTCHET", &nbjetsTCHET, &b_nbjetsTCHET);
   fChain->SetBranchAddress("nbjetsTCHPM", &nbjetsTCHPM, &b_nbjetsTCHPM);
   fChain->SetBranchAddress("nbjetsCSVM", &nbjetsCSVM, &b_nbjetsCSVM);
   fChain->SetBranchAddress("nbjetsCSVL", &nbjetsCSVL, &b_nbjetsCSVL);
   fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
   fChain->SetBranchAddress("pass_utilityHLT", &pass_utilityHLT, &b_pass_utilityHLT);
   fChain->SetBranchAddress("prescaleUtilityHLT", &prescaleUtilityHLT, &b_prescaleUtilityHLT);
   fChain->SetBranchAddress("versionUtilityHLT", &versionUtilityHLT, &b_versionUtilityHLT);
   fChain->SetBranchAddress("pass_utilityHLT_HT300_CentralJet30_BTagIP", &pass_utilityHLT_HT300_CentralJet30_BTagIP, &b_pass_utilityHLT_HT300_CentralJet30_BTagIP);
   fChain->SetBranchAddress("prescale_utilityHLT_HT300_CentralJet30_BTagIP", &prescale_utilityHLT_HT300_CentralJet30_BTagIP, &b_prescale_utilityHLT_HT300_CentralJet30_BTagIP);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("ST", &ST, &b_ST);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("METphi", &METphi, &b_METphi);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("caloMET", &caloMET, &b_caloMET);
   fChain->SetBranchAddress("METPFType1", &METPFType1, &b_METPFType1);
   fChain->SetBranchAddress("METPFType1phi", &METPFType1phi, &b_METPFType1phi);
   fChain->SetBranchAddress("bestWMass", &bestWMass, &b_bestWMass);
   fChain->SetBranchAddress("bestTopMass", &bestTopMass, &b_bestTopMass);
   fChain->SetBranchAddress("topCosHel", &topCosHel, &b_topCosHel);
   fChain->SetBranchAddress("WCosHel", &WCosHel, &b_WCosHel);
   fChain->SetBranchAddress("MT_Wlep", &MT_Wlep, &b_MT_Wlep);
   fChain->SetBranchAddress("MT_Wlep5", &MT_Wlep5, &b_MT_Wlep5);
   fChain->SetBranchAddress("MT_Wlep15", &MT_Wlep15, &b_MT_Wlep15);
   fChain->SetBranchAddress("minDeltaPhi", &minDeltaPhi, &b_minDeltaPhi);
   fChain->SetBranchAddress("minDeltaPhiAll", &minDeltaPhiAll, &b_minDeltaPhiAll);
   fChain->SetBranchAddress("minDeltaPhiAll30", &minDeltaPhiAll30, &b_minDeltaPhiAll30);
   fChain->SetBranchAddress("minDeltaPhi30_eta5_noIdAll", &minDeltaPhi30_eta5_noIdAll, &b_minDeltaPhi30_eta5_noIdAll);
   fChain->SetBranchAddress("minDeltaPhiMetTau", &minDeltaPhiMetTau, &b_minDeltaPhiMetTau);
   fChain->SetBranchAddress("deltaPhi1", &deltaPhi1, &b_deltaPhi1);
   fChain->SetBranchAddress("deltaPhi2", &deltaPhi2, &b_deltaPhi2);
   fChain->SetBranchAddress("deltaPhi3", &deltaPhi3, &b_deltaPhi3);
   fChain->SetBranchAddress("maxDeltaPhi", &maxDeltaPhi, &b_maxDeltaPhi);
   fChain->SetBranchAddress("maxDeltaPhiAll", &maxDeltaPhiAll, &b_maxDeltaPhiAll);
   fChain->SetBranchAddress("maxDeltaPhiAll30", &maxDeltaPhiAll30, &b_maxDeltaPhiAll30);
   fChain->SetBranchAddress("maxDeltaPhi30_eta5_noIdAll", &maxDeltaPhi30_eta5_noIdAll, &b_maxDeltaPhi30_eta5_noIdAll);
   fChain->SetBranchAddress("sumDeltaPhi", &sumDeltaPhi, &b_sumDeltaPhi);
   fChain->SetBranchAddress("diffDeltaPhi", &diffDeltaPhi, &b_diffDeltaPhi);
   fChain->SetBranchAddress("deltaPhiStar", &deltaPhiStar, &b_deltaPhiStar);
   fChain->SetBranchAddress("deltaPhiStar_badjet_pt", &deltaPhiStar_badjet_pt, &b_deltaPhiStar_badjet_pt);
   fChain->SetBranchAddress("deltaPhiStar_badjet_phi", &deltaPhiStar_badjet_phi, &b_deltaPhiStar_badjet_phi);
   fChain->SetBranchAddress("deltaPhiStar_badjet_eta", &deltaPhiStar_badjet_eta, &b_deltaPhiStar_badjet_eta);
   fChain->SetBranchAddress("minDeltaPhiN", &minDeltaPhiN, &b_minDeltaPhiN);
   fChain->SetBranchAddress("deltaPhiN1", &deltaPhiN1, &b_deltaPhiN1);
   fChain->SetBranchAddress("deltaPhiN2", &deltaPhiN2, &b_deltaPhiN2);
   fChain->SetBranchAddress("deltaPhiN3", &deltaPhiN3, &b_deltaPhiN3);
   fChain->SetBranchAddress("minDeltaPhi_chosenJet", &minDeltaPhi_chosenJet, &b_minDeltaPhi_chosenJet);
   fChain->SetBranchAddress("minDeltaPhiN_chosenJet", &minDeltaPhiN_chosenJet, &b_minDeltaPhiN_chosenJet);
   fChain->SetBranchAddress("maxJetMis_chosenJet", &maxJetMis_chosenJet, &b_maxJetMis_chosenJet);
   fChain->SetBranchAddress("maxJetFracMis_chosenJet", &maxJetFracMis_chosenJet, &b_maxJetFracMis_chosenJet);
   fChain->SetBranchAddress("jetMisCategoryT1", &jetMisCategoryT1, &b_jetMisCategoryT1);
   fChain->SetBranchAddress("jetMisCategoryT1_noMETfraction", &jetMisCategoryT1_noMETfraction, &b_jetMisCategoryT1_noMETfraction);
   fChain->SetBranchAddress("jetMisCategoryT2", &jetMisCategoryT2, &b_jetMisCategoryT2);
   fChain->SetBranchAddress("jetMisCategoryT2_noMETfraction", &jetMisCategoryT2_noMETfraction, &b_jetMisCategoryT2_noMETfraction);
   fChain->SetBranchAddress("minDeltaPhiN_deltaT", &minDeltaPhiN_deltaT, &b_minDeltaPhiN_deltaT);
   fChain->SetBranchAddress("deltaT1", &deltaT1, &b_deltaT1);
   fChain->SetBranchAddress("deltaT2", &deltaT2, &b_deltaT2);
   fChain->SetBranchAddress("deltaT3", &deltaT3, &b_deltaT3);
   fChain->SetBranchAddress("minDeltaPhiN_MC", &minDeltaPhiN_MC, &b_minDeltaPhiN_MC);
   fChain->SetBranchAddress("deltaPhiN_MC1", &deltaPhiN_MC1, &b_deltaPhiN_MC1);
   fChain->SetBranchAddress("deltaPhiN_MC2", &deltaPhiN_MC2, &b_deltaPhiN_MC2);
   fChain->SetBranchAddress("deltaPhiN_MC3", &deltaPhiN_MC3, &b_deltaPhiN_MC3);
   fChain->SetBranchAddress("deltaT_MC1", &deltaT_MC1, &b_deltaT_MC1);
   fChain->SetBranchAddress("deltaT_MC2", &deltaT_MC2, &b_deltaT_MC2);
   fChain->SetBranchAddress("deltaT_MC3", &deltaT_MC3, &b_deltaT_MC3);
   fChain->SetBranchAddress("minDeltaPhiN_MCq", &minDeltaPhiN_MCq, &b_minDeltaPhiN_MCq);
   fChain->SetBranchAddress("deltaPhiN_MCq1", &deltaPhiN_MCq1, &b_deltaPhiN_MCq1);
   fChain->SetBranchAddress("deltaPhiN_MCq2", &deltaPhiN_MCq2, &b_deltaPhiN_MCq2);
   fChain->SetBranchAddress("deltaPhiN_MCq3", &deltaPhiN_MCq3, &b_deltaPhiN_MCq3);
   fChain->SetBranchAddress("deltaT_MCq1", &deltaT_MCq1, &b_deltaT_MCq1);
   fChain->SetBranchAddress("deltaT_MCq2", &deltaT_MCq2, &b_deltaT_MCq2);
   fChain->SetBranchAddress("deltaT_MCq3", &deltaT_MCq3, &b_deltaT_MCq3);
   fChain->SetBranchAddress("minDeltaPhiN_MCk", &minDeltaPhiN_MCk, &b_minDeltaPhiN_MCk);
   fChain->SetBranchAddress("deltaPhiN_MCk1", &deltaPhiN_MCk1, &b_deltaPhiN_MCk1);
   fChain->SetBranchAddress("deltaPhiN_MCk2", &deltaPhiN_MCk2, &b_deltaPhiN_MCk2);
   fChain->SetBranchAddress("deltaPhiN_MCk3", &deltaPhiN_MCk3, &b_deltaPhiN_MCk3);
   fChain->SetBranchAddress("deltaT_MCk1", &deltaT_MCk1, &b_deltaT_MCk1);
   fChain->SetBranchAddress("deltaT_MCk2", &deltaT_MCk2, &b_deltaT_MCk2);
   fChain->SetBranchAddress("deltaT_MCk3", &deltaT_MCk3, &b_deltaT_MCk3);
   fChain->SetBranchAddress("deltaT1_otherPt40", &deltaT1_otherPt40, &b_deltaT1_otherPt40);
   fChain->SetBranchAddress("deltaT1_otherPt50", &deltaT1_otherPt50, &b_deltaT1_otherPt50);
   fChain->SetBranchAddress("deltaT1_otherPt70", &deltaT1_otherPt70, &b_deltaT1_otherPt70);
   fChain->SetBranchAddress("deltaT2_otherPt40", &deltaT2_otherPt40, &b_deltaT2_otherPt40);
   fChain->SetBranchAddress("deltaT2_otherPt50", &deltaT2_otherPt50, &b_deltaT2_otherPt50);
   fChain->SetBranchAddress("deltaT2_otherPt70", &deltaT2_otherPt70, &b_deltaT2_otherPt70);
   fChain->SetBranchAddress("deltaT3_otherPt40", &deltaT3_otherPt40, &b_deltaT3_otherPt40);
   fChain->SetBranchAddress("deltaT3_otherPt50", &deltaT3_otherPt50, &b_deltaT3_otherPt50);
   fChain->SetBranchAddress("deltaT3_otherPt70", &deltaT3_otherPt70, &b_deltaT3_otherPt70);
   fChain->SetBranchAddress("minDeltaPhiNAll", &minDeltaPhiNAll, &b_minDeltaPhiNAll);
   fChain->SetBranchAddress("minDeltaPhiNAll30", &minDeltaPhiNAll30, &b_minDeltaPhiNAll30);
   fChain->SetBranchAddress("minDeltaPhiN_otherEta5", &minDeltaPhiN_otherEta5, &b_minDeltaPhiN_otherEta5);
   fChain->SetBranchAddress("minDeltaPhiN_otherEta5idNo", &minDeltaPhiN_otherEta5idNo, &b_minDeltaPhiN_otherEta5idNo);
   fChain->SetBranchAddress("minDeltaPhiN_mainPt30_otherEta5idNo", &minDeltaPhiN_mainPt30_otherEta5idNo, &b_minDeltaPhiN_mainPt30_otherEta5idNo);
   fChain->SetBranchAddress("minDeltaPhiN_mainPt30Eta5_otherEta5idNo", &minDeltaPhiN_mainPt30Eta5_otherEta5idNo, &b_minDeltaPhiN_mainPt30Eta5_otherEta5idNo);
   fChain->SetBranchAddress("minDeltaPhiN_mainPt30Eta5idNo_otherEta5idNo", &minDeltaPhiN_mainPt30Eta5idNo_otherEta5idNo, &b_minDeltaPhiN_mainPt30Eta5idNo_otherEta5idNo);
   fChain->SetBranchAddress("minDeltaPhiK", &minDeltaPhiK, &b_minDeltaPhiK);
   fChain->SetBranchAddress("minDeltaPhiKw", &minDeltaPhiKw, &b_minDeltaPhiKw);
   fChain->SetBranchAddress("minDeltaPhiK_otherEta5", &minDeltaPhiK_otherEta5, &b_minDeltaPhiK_otherEta5);
   fChain->SetBranchAddress("minDeltaPhiK_otherEta5idNo", &minDeltaPhiK_otherEta5idNo, &b_minDeltaPhiK_otherEta5idNo);
   fChain->SetBranchAddress("minDeltaPhiK_mainPt30_otherEta5idNo", &minDeltaPhiK_mainPt30_otherEta5idNo, &b_minDeltaPhiK_mainPt30_otherEta5idNo);
   fChain->SetBranchAddress("minDeltaPhiK_mainPt30Eta5_otherEta5idNo", &minDeltaPhiK_mainPt30Eta5_otherEta5idNo, &b_minDeltaPhiK_mainPt30Eta5_otherEta5idNo);
   fChain->SetBranchAddress("minDeltaPhiK_mainPt30Eta5idNo_otherEta5idNo", &minDeltaPhiK_mainPt30Eta5idNo_otherEta5idNo, &b_minDeltaPhiK_mainPt30Eta5idNo_otherEta5idNo);
   fChain->SetBranchAddress("minDeltaPhiN_DJR", &minDeltaPhiN_DJR, &b_minDeltaPhiN_DJR);
   fChain->SetBranchAddress("minDeltaPhiN_DJR_otherEta5", &minDeltaPhiN_DJR_otherEta5, &b_minDeltaPhiN_DJR_otherEta5);
   fChain->SetBranchAddress("minDeltaPhiK_DJR", &minDeltaPhiK_DJR, &b_minDeltaPhiK_DJR);
   fChain->SetBranchAddress("minDeltaPhiK_DJR_otherEta5", &minDeltaPhiK_DJR_otherEta5, &b_minDeltaPhiK_DJR_otherEta5);
   fChain->SetBranchAddress("minDeltaPhiN_otherPt10", &minDeltaPhiN_otherPt10, &b_minDeltaPhiN_otherPt10);
   fChain->SetBranchAddress("minDeltaPhiN_otherPt20", &minDeltaPhiN_otherPt20, &b_minDeltaPhiN_otherPt20);
   fChain->SetBranchAddress("minDeltaPhiN_otherPt30", &minDeltaPhiN_otherPt30, &b_minDeltaPhiN_otherPt30);
   fChain->SetBranchAddress("minDeltaPhiN_otherPt40", &minDeltaPhiN_otherPt40, &b_minDeltaPhiN_otherPt40);
   fChain->SetBranchAddress("minDeltaPhiN_otherPt50", &minDeltaPhiN_otherPt50, &b_minDeltaPhiN_otherPt50);
   fChain->SetBranchAddress("minDeltaPhiN_DJR_otherPt10", &minDeltaPhiN_DJR_otherPt10, &b_minDeltaPhiN_DJR_otherPt10);
   fChain->SetBranchAddress("minDeltaPhiN_DJR_otherPt20", &minDeltaPhiN_DJR_otherPt20, &b_minDeltaPhiN_DJR_otherPt20);
   fChain->SetBranchAddress("minDeltaPhiN_DJR_otherPt30", &minDeltaPhiN_DJR_otherPt30, &b_minDeltaPhiN_DJR_otherPt30);
   fChain->SetBranchAddress("minDeltaPhiN_DJR_otherPt40", &minDeltaPhiN_DJR_otherPt40, &b_minDeltaPhiN_DJR_otherPt40);
   fChain->SetBranchAddress("minDeltaPhiN_DJR_otherPt50", &minDeltaPhiN_DJR_otherPt50, &b_minDeltaPhiN_DJR_otherPt50);
   fChain->SetBranchAddress("minDeltaPhiN_withLepton", &minDeltaPhiN_withLepton, &b_minDeltaPhiN_withLepton);
   fChain->SetBranchAddress("maxJetMis", &maxJetMis, &b_maxJetMis);
   fChain->SetBranchAddress("max2JetMis", &max2JetMis, &b_max2JetMis);
   fChain->SetBranchAddress("maxJetMisAll30", &maxJetMisAll30, &b_maxJetMisAll30);
   fChain->SetBranchAddress("max2JetMisAll30", &max2JetMisAll30, &b_max2JetMisAll30);
   fChain->SetBranchAddress("maxJetMis_signed", &maxJetMis_signed, &b_maxJetMis_signed);
   fChain->SetBranchAddress("maxJetMis30_signed", &maxJetMis30_signed, &b_maxJetMis30_signed);
   fChain->SetBranchAddress("maxJetFracMis", &maxJetFracMis, &b_maxJetFracMis);
   fChain->SetBranchAddress("max2JetFracMis", &max2JetFracMis, &b_max2JetFracMis);
   fChain->SetBranchAddress("maxJetFracMisAll30", &maxJetFracMisAll30, &b_maxJetFracMisAll30);
   fChain->SetBranchAddress("max2JetFracMisAll30", &max2JetFracMisAll30, &b_max2JetFracMisAll30);
   fChain->SetBranchAddress("deltaPhiMETJetMaxMis", &deltaPhiMETJetMaxMis, &b_deltaPhiMETJetMaxMis);
   fChain->SetBranchAddress("deltaPhiMETJetMaxMis30", &deltaPhiMETJetMaxMis30, &b_deltaPhiMETJetMaxMis30);
   fChain->SetBranchAddress("rankMaxJetMis", &rankMaxJetMis, &b_rankMaxJetMis);
   fChain->SetBranchAddress("rankMaxJetMis30", &rankMaxJetMis30, &b_rankMaxJetMis30);
   fChain->SetBranchAddress("CSVout1", &CSVout1, &b_CSVout1);
   fChain->SetBranchAddress("CSVout2", &CSVout2, &b_CSVout2);
   fChain->SetBranchAddress("CSVout3", &CSVout3, &b_CSVout3);
   fChain->SetBranchAddress("minDeltaPhiAllb30", &minDeltaPhiAllb30, &b_minDeltaPhiAllb30);
   fChain->SetBranchAddress("deltaPhib1", &deltaPhib1, &b_deltaPhib1);
   fChain->SetBranchAddress("deltaPhib2", &deltaPhib2, &b_deltaPhib2);
   fChain->SetBranchAddress("deltaPhib3", &deltaPhib3, &b_deltaPhib3);
   fChain->SetBranchAddress("minDeltaPhiMETMuonsAll", &minDeltaPhiMETMuonsAll, &b_minDeltaPhiMETMuonsAll);
   fChain->SetBranchAddress("minDeltaPhiN_lostJet", &minDeltaPhiN_lostJet, &b_minDeltaPhiN_lostJet);
   fChain->SetBranchAddress("deltaPhiN1_lostJet", &deltaPhiN1_lostJet, &b_deltaPhiN1_lostJet);
   fChain->SetBranchAddress("deltaPhiN2_lostJet", &deltaPhiN2_lostJet, &b_deltaPhiN2_lostJet);
   fChain->SetBranchAddress("deltaPhiN3_lostJet", &deltaPhiN3_lostJet, &b_deltaPhiN3_lostJet);
   fChain->SetBranchAddress("passBadECAL_METphi53_n10_s12", &passBadECAL_METphi53_n10_s12, &b_passBadECAL_METphi53_n10_s12);
   fChain->SetBranchAddress("worstMisA_badECAL_METphi53_n10_s12", &worstMisA_badECAL_METphi53_n10_s12, &b_worstMisA_badECAL_METphi53_n10_s12);
   fChain->SetBranchAddress("worstMisF_badECAL_METphi53_n10_s12", &worstMisF_badECAL_METphi53_n10_s12, &b_worstMisF_badECAL_METphi53_n10_s12);
   fChain->SetBranchAddress("passBadECAL_METphi33_n10_s12", &passBadECAL_METphi33_n10_s12, &b_passBadECAL_METphi33_n10_s12);
   fChain->SetBranchAddress("worstMisA_badECAL_METphi33_n10_s12", &worstMisA_badECAL_METphi33_n10_s12, &b_worstMisA_badECAL_METphi33_n10_s12);
   fChain->SetBranchAddress("worstMisF_badECAL_METphi33_n10_s12", &worstMisF_badECAL_METphi33_n10_s12, &b_worstMisF_badECAL_METphi33_n10_s12);
   fChain->SetBranchAddress("passBadECAL_METphi52_n10_s12", &passBadECAL_METphi52_n10_s12, &b_passBadECAL_METphi52_n10_s12);
   fChain->SetBranchAddress("worstMisA_badECAL_METphi52_n10_s12", &worstMisA_badECAL_METphi52_n10_s12, &b_worstMisA_badECAL_METphi52_n10_s12);
   fChain->SetBranchAddress("worstMisF_badECAL_METphi52_n10_s12", &worstMisF_badECAL_METphi52_n10_s12, &b_worstMisF_badECAL_METphi52_n10_s12);
   fChain->SetBranchAddress("passBadECAL_METphi53_n5_s12", &passBadECAL_METphi53_n5_s12, &b_passBadECAL_METphi53_n5_s12);
   fChain->SetBranchAddress("passBadECAL_METphi53_crack", &passBadECAL_METphi53_crack, &b_passBadECAL_METphi53_crack);
   fChain->SetBranchAddress("worstMisA_badECAL_METphi53_crack", &worstMisA_badECAL_METphi53_crack, &b_worstMisA_badECAL_METphi53_crack);
   fChain->SetBranchAddress("worstMisF_badECAL_METphi53_crack", &worstMisF_badECAL_METphi53_crack, &b_worstMisF_badECAL_METphi53_crack);
   fChain->SetBranchAddress("passBadECAL_METphi52_crack", &passBadECAL_METphi52_crack, &b_passBadECAL_METphi52_crack);
   fChain->SetBranchAddress("worstMisA_badECAL_METphi52_crack", &worstMisA_badECAL_METphi52_crack, &b_worstMisA_badECAL_METphi52_crack);
   fChain->SetBranchAddress("worstMisF_badECAL_METphi52_crack", &worstMisF_badECAL_METphi52_crack, &b_worstMisF_badECAL_METphi52_crack);
   fChain->SetBranchAddress("passBadECAL_METphi53_crackF", &passBadECAL_METphi53_crackF, &b_passBadECAL_METphi53_crackF);
   fChain->SetBranchAddress("worstMisA_badECAL_METphi53_crackF", &worstMisA_badECAL_METphi53_crackF, &b_worstMisA_badECAL_METphi53_crackF);
   fChain->SetBranchAddress("worstMisF_badECAL_METphi53_crackF", &worstMisF_badECAL_METphi53_crackF, &b_worstMisF_badECAL_METphi53_crackF);
   fChain->SetBranchAddress("passBadECAL_METphi52_crackF", &passBadECAL_METphi52_crackF, &b_passBadECAL_METphi52_crackF);
   fChain->SetBranchAddress("worstMisA_badECAL_METphi52_crackF", &worstMisA_badECAL_METphi52_crackF, &b_worstMisA_badECAL_METphi52_crackF);
   fChain->SetBranchAddress("worstMisF_badECAL_METphi52_crackF", &worstMisF_badECAL_METphi52_crackF, &b_worstMisF_badECAL_METphi52_crackF);
   fChain->SetBranchAddress("jetpt1", &jetpt1, &b_jetpt1);
   fChain->SetBranchAddress("jetgenpt1", &jetgenpt1, &b_jetgenpt1);
   fChain->SetBranchAddress("jeteta1", &jeteta1, &b_jeteta1);
   fChain->SetBranchAddress("jetgeneta1", &jetgeneta1, &b_jetgeneta1);
   fChain->SetBranchAddress("jetphi1", &jetphi1, &b_jetphi1);
   fChain->SetBranchAddress("jetgenphi1", &jetgenphi1, &b_jetgenphi1);
   fChain->SetBranchAddress("jetenergy1", &jetenergy1, &b_jetenergy1);
   fChain->SetBranchAddress("jetflavor1", &jetflavor1, &b_jetflavor1);
   fChain->SetBranchAddress("jetchargedhadronfrac1", &jetchargedhadronfrac1, &b_jetchargedhadronfrac1);
   fChain->SetBranchAddress("jetchargedhadronmult1", &jetchargedhadronmult1, &b_jetchargedhadronmult1);
   fChain->SetBranchAddress("jetpt2", &jetpt2, &b_jetpt2);
   fChain->SetBranchAddress("jetgenpt2", &jetgenpt2, &b_jetgenpt2);
   fChain->SetBranchAddress("jeteta2", &jeteta2, &b_jeteta2);
   fChain->SetBranchAddress("jetgeneta2", &jetgeneta2, &b_jetgeneta2);
   fChain->SetBranchAddress("jetphi2", &jetphi2, &b_jetphi2);
   fChain->SetBranchAddress("jetgenphi2", &jetgenphi2, &b_jetgenphi2);
   fChain->SetBranchAddress("jetenergy2", &jetenergy2, &b_jetenergy2);
   fChain->SetBranchAddress("jetflavor2", &jetflavor2, &b_jetflavor2);
   fChain->SetBranchAddress("jetchargedhadronfrac2", &jetchargedhadronfrac2, &b_jetchargedhadronfrac2);
   fChain->SetBranchAddress("jetchargedhadronmult2", &jetchargedhadronmult2, &b_jetchargedhadronmult2);
   fChain->SetBranchAddress("jetpt3", &jetpt3, &b_jetpt3);
   fChain->SetBranchAddress("jetgenpt3", &jetgenpt3, &b_jetgenpt3);
   fChain->SetBranchAddress("jeteta3", &jeteta3, &b_jeteta3);
   fChain->SetBranchAddress("jetgeneta3", &jetgeneta3, &b_jetgeneta3);
   fChain->SetBranchAddress("jetphi3", &jetphi3, &b_jetphi3);
   fChain->SetBranchAddress("jetgenphi3", &jetgenphi3, &b_jetgenphi3);
   fChain->SetBranchAddress("jetenergy3", &jetenergy3, &b_jetenergy3);
   fChain->SetBranchAddress("jetflavor3", &jetflavor3, &b_jetflavor3);
   fChain->SetBranchAddress("jetchargedhadronfrac3", &jetchargedhadronfrac3, &b_jetchargedhadronfrac3);
   fChain->SetBranchAddress("jetchargedhadronmult3", &jetchargedhadronmult3, &b_jetchargedhadronmult3);
   fChain->SetBranchAddress("bjetpt1", &bjetpt1, &b_bjetpt1);
   fChain->SetBranchAddress("bjeteta1", &bjeteta1, &b_bjeteta1);
   fChain->SetBranchAddress("bjetphi1", &bjetphi1, &b_bjetphi1);
   fChain->SetBranchAddress("bjetenergy1", &bjetenergy1, &b_bjetenergy1);
   fChain->SetBranchAddress("bjetflavor1", &bjetflavor1, &b_bjetflavor1);
   fChain->SetBranchAddress("bjetchargedhadronfrac1", &bjetchargedhadronfrac1, &b_bjetchargedhadronfrac1);
   fChain->SetBranchAddress("bjetchargedhadronmult1", &bjetchargedhadronmult1, &b_bjetchargedhadronmult1);
   fChain->SetBranchAddress("bjetpt2", &bjetpt2, &b_bjetpt2);
   fChain->SetBranchAddress("bjeteta2", &bjeteta2, &b_bjeteta2);
   fChain->SetBranchAddress("bjetphi2", &bjetphi2, &b_bjetphi2);
   fChain->SetBranchAddress("bjetenergy2", &bjetenergy2, &b_bjetenergy2);
   fChain->SetBranchAddress("bjetflavor2", &bjetflavor2, &b_bjetflavor2);
   fChain->SetBranchAddress("bjetchargedhadronfrac2", &bjetchargedhadronfrac2, &b_bjetchargedhadronfrac2);
   fChain->SetBranchAddress("bjetchargedhadronmult2", &bjetchargedhadronmult2, &b_bjetchargedhadronmult2);
   fChain->SetBranchAddress("bjetpt3", &bjetpt3, &b_bjetpt3);
   fChain->SetBranchAddress("bjeteta3", &bjeteta3, &b_bjeteta3);
   fChain->SetBranchAddress("bjetphi3", &bjetphi3, &b_bjetphi3);
   fChain->SetBranchAddress("bjetenergy3", &bjetenergy3, &b_bjetenergy3);
   fChain->SetBranchAddress("bjetflavor3", &bjetflavor3, &b_bjetflavor3);
   fChain->SetBranchAddress("bjetchargedhadronfrac3", &bjetchargedhadronfrac3, &b_bjetchargedhadronfrac3);
   fChain->SetBranchAddress("bjetchargedhadronmult3", &bjetchargedhadronmult3, &b_bjetchargedhadronmult3);
   fChain->SetBranchAddress("eleet1", &eleet1, &b_eleet1);
   fChain->SetBranchAddress("elecharge1", &elecharge1, &b_elecharge1);
   fChain->SetBranchAddress("elephi1", &elephi1, &b_elephi1);
   fChain->SetBranchAddress("eleeta1", &eleeta1, &b_eleeta1);
   fChain->SetBranchAddress("muonpt1", &muonpt1, &b_muonpt1);
   fChain->SetBranchAddress("muoncharge1", &muoncharge1, &b_muoncharge1);
   fChain->SetBranchAddress("muonphi1", &muonphi1, &b_muonphi1);
   fChain->SetBranchAddress("muoneta1", &muoneta1, &b_muoneta1);
   fChain->SetBranchAddress("muoniso1", &muoniso1, &b_muoniso1);
   fChain->SetBranchAddress("muonchhadiso1", &muonchhadiso1, &b_muonchhadiso1);
   fChain->SetBranchAddress("muonphotoniso1", &muonphotoniso1, &b_muonphotoniso1);
   fChain->SetBranchAddress("muonneutralhadiso1", &muonneutralhadiso1, &b_muonneutralhadiso1);
   fChain->SetBranchAddress("eleet2", &eleet2, &b_eleet2);
   fChain->SetBranchAddress("elephi2", &elephi2, &b_elephi2);
   fChain->SetBranchAddress("eleeta2", &eleeta2, &b_eleeta2);
   fChain->SetBranchAddress("muonpt2", &muonpt2, &b_muonpt2);
   fChain->SetBranchAddress("muonphi2", &muonphi2, &b_muonphi2);
   fChain->SetBranchAddress("muoneta2", &muoneta2, &b_muoneta2);
   fChain->SetBranchAddress("taupt1", &taupt1, &b_taupt1);
   fChain->SetBranchAddress("taueta1", &taueta1, &b_taueta1);
   fChain->SetBranchAddress("eleRelIso", &eleRelIso, &b_eleRelIso);
   fChain->SetBranchAddress("muonRelIso", &muonRelIso, &b_muonRelIso);
   fChain->SetBranchAddress("recomuonpt1", &recomuonpt1, &b_recomuonpt1);
   fChain->SetBranchAddress("recomuonphi1", &recomuonphi1, &b_recomuonphi1);
   fChain->SetBranchAddress("recomuoneta1", &recomuoneta1, &b_recomuoneta1);
   fChain->SetBranchAddress("recomuoniso1", &recomuoniso1, &b_recomuoniso1);
   fChain->SetBranchAddress("recomuonmindphijet1", &recomuonmindphijet1, &b_recomuonmindphijet1);
   fChain->SetBranchAddress("recomuonpt2", &recomuonpt2, &b_recomuonpt2);
   fChain->SetBranchAddress("recomuonphi2", &recomuonphi2, &b_recomuonphi2);
   fChain->SetBranchAddress("recomuoneta2", &recomuoneta2, &b_recomuoneta2);
   fChain->SetBranchAddress("recomuoniso2", &recomuoniso2, &b_recomuoniso2);
   fChain->SetBranchAddress("recomuonmindphijet2", &recomuonmindphijet2, &b_recomuonmindphijet2);
   fChain->SetBranchAddress("minDeltaPhiN_Luke", &minDeltaPhiN_Luke, &b_minDeltaPhiN_Luke);
   fChain->SetBranchAddress("maxDeltaPhiN_Luke", &maxDeltaPhiN_Luke, &b_maxDeltaPhiN_Luke);
   fChain->SetBranchAddress("deltaPhiN1_Luke", &deltaPhiN1_Luke, &b_deltaPhiN1_Luke);
   fChain->SetBranchAddress("deltaPhiN2_Luke", &deltaPhiN2_Luke, &b_deltaPhiN2_Luke);
   fChain->SetBranchAddress("deltaPhiN3_Luke", &deltaPhiN3_Luke, &b_deltaPhiN3_Luke);
   fChain->SetBranchAddress("minTransverseMETSignificance", &minTransverseMETSignificance, &b_minTransverseMETSignificance);
   fChain->SetBranchAddress("maxTransverseMETSignificance", &maxTransverseMETSignificance, &b_maxTransverseMETSignificance);
   fChain->SetBranchAddress("transverseMETSignificance1", &transverseMETSignificance1, &b_transverseMETSignificance1);
   fChain->SetBranchAddress("transverseMETSignificance2", &transverseMETSignificance2, &b_transverseMETSignificance2);
   fChain->SetBranchAddress("transverseMETSignificance3", &transverseMETSignificance3, &b_transverseMETSignificance3);
   fChain->SetBranchAddress("njets_lostJet", &njets_lostJet, &b_njets_lostJet);
   fChain->SetBranchAddress("nbjets_lostJet", &nbjets_lostJet, &b_nbjets_lostJet);
   fChain->SetBranchAddress("minDeltaPhiN_Luke_lostJet", &minDeltaPhiN_Luke_lostJet, &b_minDeltaPhiN_Luke_lostJet);
   fChain->SetBranchAddress("maxDeltaPhiN_Luke_lostJet", &maxDeltaPhiN_Luke_lostJet, &b_maxDeltaPhiN_Luke_lostJet);
   fChain->SetBranchAddress("deltaPhiN1_Luke_lostJet", &deltaPhiN1_Luke_lostJet, &b_deltaPhiN1_Luke_lostJet);
   fChain->SetBranchAddress("deltaPhiN2_Luke_lostJet", &deltaPhiN2_Luke_lostJet, &b_deltaPhiN2_Luke_lostJet);
   fChain->SetBranchAddress("deltaPhiN3_Luke_lostJet", &deltaPhiN3_Luke_lostJet, &b_deltaPhiN3_Luke_lostJet);
   fChain->SetBranchAddress("minTransverseMETSignificance_lostJet", &minTransverseMETSignificance_lostJet, &b_minTransverseMETSignificance_lostJet);
   fChain->SetBranchAddress("maxTransverseMETSignificance_lostJet", &maxTransverseMETSignificance_lostJet, &b_maxTransverseMETSignificance_lostJet);
   fChain->SetBranchAddress("transverseMETSignificance1_lostJet", &transverseMETSignificance1_lostJet, &b_transverseMETSignificance1_lostJet);
   fChain->SetBranchAddress("transverseMETSignificance2_lostJet", &transverseMETSignificance2_lostJet, &b_transverseMETSignificance2_lostJet);
   fChain->SetBranchAddress("transverseMETSignificance3_lostJet", &transverseMETSignificance3_lostJet, &b_transverseMETSignificance3_lostJet);
   fChain->SetBranchAddress("nLostJet", &nLostJet, &b_nLostJet);
   Notify();
}

Bool_t ISRweights::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ISRweights::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ISRweights::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ISRweights_cxx
