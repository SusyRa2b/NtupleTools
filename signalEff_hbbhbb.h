//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 18 02:38:27 2013 by ROOT version 5.32/00
// from TTree reducedTree/tree with minimal cuts
// found on file: reducedTrees/reducedTree.JES0_JER0_PFMETTypeI_METunc0_PUunc0_hpt20.TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71.root
//////////////////////////////////////////////////////////

#ifndef signalEff_hbbhbb_h
#define signalEff_hbbhbb_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include "TH3.h"

#include <iostream>


enum IsrMode { kNoIsrWeight, kIsr0, kIsrUp,kIsrDown };
/* TO DO
enum TriggerMode { kNoTrigWeight, kTrig0, kTrigUp,kTrigDown };
*/
enum BTagMode {kBTag0,kBTagHfUp,kBTagHfDown,kBTagLfUp,kBTagLfDown};


// Header file for the classes stored in the TTree if any.
class SearchRegion { //a DIFFERENT class than my other SearchRegion class...oh well
public:
  SearchRegion(int nb, float minMETs, float maxMETs, bool isSB, TString pdfset, int pdfindex);
  ~SearchRegion();
  //  void Print() const;

  int nb_; //no min or max: use the mutually exclusive bins defined for the hbbhbb analysis
  float minMETs_;
  float maxMETs_;
  bool isSB_;

  TString pdfset_;
  int pdfindex_;

  TString id() const;

  //  bool operator==(const SearchRegion& other);
  //  bool operator!=(const SearchRegion& other) {return !((*this)==other);}

};


SearchRegion::SearchRegion(int nb, float minMETs, float maxMETs, bool isSB, TString pdfset, int pdfindex) :
  nb_(nb),minMETs_(minMETs),maxMETs_(maxMETs),isSB_(isSB),pdfset_(pdfset),pdfindex_(pdfindex)
{
}

SearchRegion::~SearchRegion() {}

TString SearchRegion::id() const {

  TString sbstring = isSB_ ? "_SB" : "" ;

  TString theid;
  theid.Form("b%d_MET%.0fto%.0f%s",nb_,minMETs_,maxMETs_,sbstring.Data());

  if (pdfset_ != "none") {
    TString suffix;
    suffix.Form("_%s%d",pdfset_.Data(),pdfindex_);
    theid+=suffix;
  }

  return theid;
}

// Fixed size dimensions of array or collections stored in the TTree if any.

class signalEff_hbbhbb {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        weight;
   Double_t        weight2;
   Double_t        weight3;
   Double_t        scanCrossSection;
   ULong64_t       runNumber;
   ULong64_t       lumiSection;
   ULong64_t       eventNumber;
   Int_t           m0;
   Int_t           m12;
   Int_t           mIntermediate;
   Int_t           ttbarDecayCode;
   Int_t           lostLeptonCode;
   Int_t           hadWcode;
   Float_t         genTopPtLep;
   Float_t         genTopPtHad;
   Float_t         genWPtLep;
   Float_t         genWPtHad;
   Float_t         minDeltaRLeptonJet;
   Float_t         minDeltaRLeptonJetReco;
   Int_t           minDeltaRJetFlavor;
   Float_t         PUweight;
   Float_t         PUweightSystVar;
   Float_t         hlteff;
   Float_t         pdfWeightsCTEQ[45];
   Float_t         pdfWeightsMSTW[41];
   Float_t         pdfWeightsNNPDF[100];
   Float_t         PIDweight;
   Float_t         BTagWeight;
   Float_t         BTagWeightLMT;
   Float_t         BTagWeightLMT_HFplus;
   Float_t         BTagWeightLMT_HFminus;
   Float_t         BTagWeightLMT_LFplus;
   Float_t         BTagWeightLMT_LFminus;
   Float_t         SUSY_msq12[2];
   Float_t         SUSY_msq23[2];
   Float_t         SUSY_gluino_pt[2];
   Float_t         SUSY_top_pt[2];
   Float_t         SUSY_topbar_pt[2];
   Float_t         SUSY_chi0_pt[2];
   Float_t         topPtWeight;
   Float_t         topPt;
   Float_t         topPtWeightOfficial;
   Float_t         higgsMass[6];
   Float_t         higgsSumpt[6];
   Float_t         hhDeltaR;
   Float_t         hhDeltaEta;
   Float_t         hhDeltaPhi;
   Float_t         hhDeltaPhiWrong1;
   Float_t         hhDeltaPhiWrong2;
   Float_t         hhDeltaEtaWrong1;
   Float_t         hhDeltaEtaWrong2;
   Float_t         higgs1CosHelCorrect;
   Float_t         higgs1CosHelWrong1;
   Float_t         higgs1CosHelWrong2;
   Float_t         higgs1b1pt;
   Float_t         higgs1b1phi;
   Float_t         higgs1b1eta;
   Float_t         higgs1b2pt;
   Float_t         higgs1b2phi;
   Float_t         higgs1b2eta;
   Float_t         higgs2b1pt;
   Float_t         higgs2b1phi;
   Float_t         higgs2b1eta;
   Float_t         higgs2b2pt;
   Float_t         higgs2b2phi;
   Float_t         higgs2b2eta;
   Float_t         higgsb1pt;
   Float_t         higgsb2pt;
   Float_t         higgsb3pt;
   Float_t         higgsb4pt;
   Float_t         higgsDeltaRminCorrect;
   Float_t         higgsDeltaRminWrong1;
   Float_t         higgsDeltaRminWrong2;
   Float_t         higgsDeltaRmaxCorrect;
   Float_t         higgsDeltaRmaxWrong1;
   Float_t         higgsDeltaRmaxWrong2;
   Float_t         higgsDeltaPhiRminCorrect;
   Float_t         higgsDeltaPhiRminWrong1;
   Float_t         higgsDeltaPhiRminWrong2;
   Float_t         higgsDeltaPhiRmaxCorrect;
   Float_t         higgsDeltaPhiRmaxWrong1;
   Float_t         higgsDeltaPhiRmaxWrong2;
   Float_t         prob2b;
   Float_t         prob3b;
   Float_t         prob4b;
   Float_t         prob0;
   Float_t         probge1;
   Float_t         prob1;
   Float_t         probge2;
   Float_t         prob2;
   Float_t         probge3;
   Float_t         prob3;
   Float_t         probge4;
   Float_t         prob0_HFplus;
   Float_t         probge1_HFplus;
   Float_t         prob1_HFplus;
   Float_t         probge2_HFplus;
   Float_t         prob2_HFplus;
   Float_t         probge3_HFplus;
   Float_t         prob3_HFplus;
   Float_t         probge4_HFplus;
   Float_t         prob0_HFminus;
   Float_t         probge1_HFminus;
   Float_t         prob1_HFminus;
   Float_t         probge2_HFminus;
   Float_t         prob2_HFminus;
   Float_t         probge3_HFminus;
   Float_t         prob3_HFminus;
   Float_t         probge4_HFminus;
   Float_t         prob0_LFplus;
   Float_t         probge1_LFplus;
   Float_t         prob1_LFplus;
   Float_t         probge2_LFplus;
   Float_t         prob2_LFplus;
   Float_t         probge3_LFplus;
   Float_t         prob3_LFplus;
   Float_t         probge4_LFplus;
   Float_t         prob0_LFminus;
   Float_t         probge1_LFminus;
   Float_t         prob1_LFminus;
   Float_t         probge2_LFminus;
   Float_t         prob2_LFminus;
   Float_t         probge3_LFminus;
   Float_t         prob3_LFminus;
   Float_t         probge4_LFminus;
   Bool_t          cutPV;
   Bool_t          cutTrigger;
   Bool_t          cutTrigger2;
   Bool_t          csctighthaloFilter;
   Bool_t          eenoiseFilter;
   Bool_t          greedymuonFilter;
   Bool_t          hbhenoiseFilter;
   Bool_t          inconsistentmuonFilter;
   Bool_t          ra2ecaltpFilter;
   Bool_t          scrapingvetoFilter;
   Bool_t          trackingfailureFilter;
   Bool_t          badjetFilter;
   Bool_t          passCleaning;
   Int_t           PBNRcode;
   Bool_t          buggyEvent;
   Int_t           maxTOBTECjetDeltaMult;
   Int_t           TOBTECjetChMult;
   Int_t           nGoodPV;
   Float_t         PV0_x;
   Float_t         PV0_xErr;
   Float_t         PV0_y;
   Float_t         PV0_yErr;
   Float_t         PV0_z;
   Float_t         PV0_zErr;
   Float_t         PV_x[60];
   Float_t         PV_xErr[60];
   Float_t         PV_y[60];
   Float_t         PV_yErr[60];
   Float_t         PV_z[60];
   Float_t         PV_zErr[60];
   Float_t         PV_chi2[60];
   Float_t         PV_ndof[60];
   Float_t         PV_isFake[60];
   Float_t         PV_isGood[60];
   Float_t         PV_isValid[60];
   Float_t         PV_tracksSize[60];
   Int_t           PV_nearestZindex;
   Float_t         rho_kt6PFJetsForIsolation;
   Float_t         BS0_x;
   Float_t         BS0_xErr;
   Float_t         BS0_y;
   Float_t         BS0_yErr;
   Float_t         BS0_z;
   Float_t         BS0_zErr;
   Int_t           SUSY_nb;
   Int_t           SUSY_process;
   Float_t         SUSY_recoilPt;
   Float_t         SUSY_ISRweight;
   Float_t         SUSY_ISRweightSystDown;
   Int_t           njets;
   Int_t           njets30;
   Int_t           njets20;
   Int_t           nPUjets30;
   Int_t           nPUjets20;
   Int_t           njets20_5p0;
   Int_t           njetsHiggsMatch20;
   Int_t           ncompleteHiggsReco20;
   Int_t           njetsHiggsMatch20_eta5;
   Int_t           nbjets;
   Int_t           nbjets30;
   Int_t           nbjetsTweaked;
   Int_t           ntruebjets;
   Int_t           nGluonsSplitToLF;
   Int_t           nGluonsSplitToC;
   Int_t           nGluonsSplitToB;
   Int_t           nlfjetsFromGluons;
   Int_t           ncjetsFromGluons;
   Int_t           nbjetsFromGluons;
   Int_t           nElectrons;
   Int_t           nElectronsBug;
   Int_t           nMuons;
   Int_t           nElectrons5;
   Int_t           nMuons5;
   Int_t           nElectrons15;
   Int_t           nMuons15;
   Int_t           nTightMuons;
   Int_t           nElectrons20;
   Int_t           nMuons20;
   Int_t           nElectronsNoRelIso;
   Int_t           nMuonsNoRelIso;
   Int_t           nTausVLoose;
   Int_t           nTausLoose;
   Int_t           nTausMedium;
   Int_t           nTausTight;
   Int_t           nCorrectRecoStop;
   Float_t         bestZmass;
   Float_t         mjj1;
   Float_t         mjj2;
   Float_t         mjjdiff;
   Float_t         higgsMbb1;
   Float_t         higgsMbb2;
   Float_t         higgsMjj1;
   Float_t         higgsMjj2;
   Float_t         higgsMbb1delta;
   Float_t         higgsMbb2delta;
   Float_t         higgsMbb1MassDiff;
   Float_t         higgsMbb2MassDiff;
   Bool_t          higgsMbbMassDiff_tiebreak;
   Float_t         higgsMbb1MassDiffAll;
   Float_t         higgsMbb2MassDiffAll;
   Int_t           higgsMbb1MassDiff_correct;
   Int_t           higgsMbb1MassDiffAll_correct;
   Float_t         deltaPhi_hh;
   Float_t         deltaRmax_hh;
   Float_t         deltaRmin_hh;
   Float_t         deltaEta_hh;
   Float_t         sumPt_hh;
   Float_t         higgsPt1;
   Float_t         higgsPt2;
   Float_t         deltaRmax_hhAll;
   Float_t         deltaRmin_hhAll;
   Float_t         higgs1jetpt1;
   Float_t         higgs1jetpt2;
   Float_t         higgs2jetpt1;
   Float_t         higgs2jetpt2;
   Float_t         higgs1jeteta1;
   Float_t         higgs1jeteta2;
   Float_t         higgs2jeteta1;
   Float_t         higgs2jeteta2;
   Float_t         higgs1CSV1;
   Float_t         higgs1CSV2;
   Float_t         higgs2CSV1;
   Float_t         higgs2CSV2;
   Float_t         higgs1partonId1;
   Float_t         higgs1partonId2;
   Float_t         higgs2partonId1;
   Float_t         higgs2partonId2;
   Float_t         higgs1partonFlavor1;
   Float_t         higgs1partonFlavor2;
   Float_t         higgs2partonFlavor1;
   Float_t         higgs2partonFlavor2;
   Float_t         higgs1partonMomId1;
   Float_t         higgs1partonMomId2;
   Float_t         higgs2partonMomId1;
   Float_t         higgs2partonMomId2;
   Bool_t          higgs1tauMatch1;
   Bool_t          higgs1tauMatch2;
   Bool_t          higgs2tauMatch1;
   Bool_t          higgs2tauMatch2;
   Float_t         higgsWCandMass;
   Float_t         higgsWCandMassAll;
   Float_t         maxDelPhiThrustJet;
   Float_t         minDelPhiThrustMET;
   Float_t         mjjb1;
   Float_t         mjjb2;
   Float_t         topPT1;
   Float_t         topPT2;
   Int_t           nbjetsSSVM;
   Int_t           nbjetsSSVHPT;
   Int_t           nbjetsTCHPT;
   Int_t           nbjetsTCHET;
   Int_t           nbjetsTCHPM;
   Int_t           nbjetsCSVT;
   Int_t           nbjetsCSVM;
   Int_t           nbjetsCSVL;
   Bool_t          isRealData;
   Bool_t          pass_utilityHLT;
   UInt_t          prescaleUtilityHLT;
   UInt_t          versionUtilityHLT;
   Bool_t          pass_utilityPrescaleModuleHLT;
   Bool_t          pass_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45;
   Bool_t          pass_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45;
   Bool_t          pass_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned;
   Bool_t          pass_DiCentralPFJet30_PFMET80;
   Bool_t          pass_DiCentralPFJet30_PFMET80_BTagCSV07;
   Bool_t          pass_DiCentralPFJet30_PFMHT80;
   Bool_t          pass_DiCentralPFJet50_PFMET80;
   Bool_t          pass_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80;
   Bool_t          pass_DiPFJet80_DiPFJet30_BTagCSVd07d05;
   Bool_t          pass_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Bool_t          pass_HT200;
   Bool_t          pass_HT250;
   Bool_t          pass_HT250_AlphaT0p55;
   Bool_t          pass_HT300;
   Bool_t          pass_HT300_AlphaT0p53;
   Bool_t          pass_IsoMu24;
   Bool_t          pass_IsoMu24_eta2p1;
   Bool_t          pass_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV;
   Bool_t          pass_L1ETM40;
   Bool_t          pass_MET120_HBHENoiseCleaned;
   Bool_t          pass_MET200;
   Bool_t          pass_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Bool_t          pass_Mu17_Mu8;
   Bool_t          pass_Mu8_DiJet30;
   Bool_t          pass_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Bool_t          pass_PFHT350;
   Bool_t          pass_PFHT350_Mu15_PFMET45;
   Bool_t          pass_PFHT350_PFMET100;
   Bool_t          pass_PFHT650;
   Bool_t          pass_PFJet80;
   Bool_t          pass_PFMET150;
   Bool_t          pass_PFNoPUHT350;
   Bool_t          pass_PFNoPUHT350_Mu15_PFMET45;
   Bool_t          pass_PFNoPUHT350_PFMET100;
   Bool_t          pass_PFNoPUHT650;
   Bool_t          pass_Photon135;
   Bool_t          pass_Photon150;
   Bool_t          pass_QuadJet80;
   Bool_t          pass_SixJet45;
   Bool_t          passMC_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45;
   Bool_t          passMC_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45;
   Bool_t          passMC_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned;
   Bool_t          passMC_DiCentralPFJet30_PFMET80;
   Bool_t          passMC_DiCentralPFJet30_PFMET80_BTagCSV07;
   Bool_t          passMC_DiCentralPFJet30_PFMHT80;
   Bool_t          passMC_DiCentralPFJet50_PFMET80;
   Bool_t          passMC_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80;
   Bool_t          passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05;
   Bool_t          passMC_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Bool_t          passMC_HT200;
   Bool_t          passMC_HT250;
   Bool_t          passMC_HT250_AlphaT0p55;
   Bool_t          passMC_HT300;
   Bool_t          passMC_HT300_AlphaT0p53;
   Bool_t          passMC_IsoMu24;
   Bool_t          passMC_IsoMu24_eta2p1;
   Bool_t          passMC_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV;
   Bool_t          passMC_L1ETM40;
   Bool_t          passMC_MET120_HBHENoiseCleaned;
   Bool_t          passMC_MET200;
   Bool_t          passMC_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Bool_t          passMC_Mu17_Mu8;
   Bool_t          passMC_Mu8_DiJet30;
   Bool_t          passMC_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Bool_t          passMC_PFHT350;
   Bool_t          passMC_PFHT350_Mu15_PFMET45;
   Bool_t          passMC_PFHT350_PFMET100;
   Bool_t          passMC_PFHT650;
   Bool_t          passMC_PFJet80;
   Bool_t          passMC_PFMET150;
   Bool_t          passMC_PFNoPUHT350;
   Bool_t          passMC_PFNoPUHT350_Mu15_PFMET45;
   Bool_t          passMC_PFNoPUHT350_PFMET100;
   Bool_t          passMC_PFNoPUHT650;
   Bool_t          passMC_Photon135;
   Bool_t          passMC_Photon150;
   Bool_t          passMC_QuadJet80;
   Bool_t          passMC_SixJet45;
   Float_t         HT;
   Float_t         HT30;
   Float_t         HT20;
   Float_t         ST;
   Float_t         STeff;
   Float_t         MET;
   Float_t         METphi;
   Float_t         MHT;
   Float_t         METsig;
   Float_t         METsig00;
   Float_t         METsig10;
   Float_t         METsig11;
   Float_t         METsig_2012;
   Float_t         METsig00_2012;
   Float_t         METsig10_2012;
   Float_t         METsig11_2012;
   Float_t         caloMET;
   Float_t         caloMETphi;
   Float_t         rawPFMET;
   Float_t         rawPFMETphi;
   Float_t         caloMETwithHO;
   Float_t         caloMETwithHOphi;
   Float_t         unclusteredMET;
   Float_t         unclusteredMETphi;
   Float_t         bestWMass;
   Float_t         bestTopMass;
   Float_t         topCosHel;
   Float_t         WCosHel;
   Float_t         MT_b;
   Float_t         MT_Hbb1;
   Float_t         MT_Hbb2;
   Float_t         MT_Wlep;
   Float_t         MT_Wlep5;
   Float_t         MT_Wlep15;
   Float_t         maxDeltaPhi_bbbb;
   Float_t         maxDeltaPhi_bb_bb;
   Int_t           MT_bestCSV_gencode;
   Float_t         MT_bestCSV;
   Float_t         MT_jim;
   Float_t         minMT_jetMET;
   Float_t         minMT_bMET;
   Float_t         maxMT_bMET;
   Float_t         deltaThetaT;
   Float_t         deltaR_bestTwoCSV;
   Float_t         deltaPhi_bestTwoCSV;
   Float_t         mjj_bestTwoCSV;
   Float_t         deltaR_closestB;
   Float_t         mjj_closestB;
   Float_t         mjj_closestB20;
   Bool_t          mjj_closestB20_correct;
   Float_t         deltaR_closestB20;
   Float_t         cosHel_closestB20;
   Float_t         deltaPhi_hb20;
   Float_t         sumPtMjjDiff_closestB20;
   Float_t         mjj_h125;
   Float_t         mjj_highptCSVM;
   Float_t         mjj_highptCSVT;
   Float_t         minDeltaPhi;
   Float_t         minDeltaPhiAll;
   Float_t         minDeltaPhiAll20;
   Float_t         minDeltaPhi20;
   Float_t         minDeltaPhi30;
   Float_t         minDeltaPhi20_eta5_noIdAll;
   Float_t         minDeltaPhi20_eta5_noIdAll_nobeta;
   Float_t         deltaPhi1;
   Float_t         deltaPhi2;
   Float_t         deltaPhi3;
   Float_t         deltaPhi1_20;
   Float_t         deltaPhi2_20;
   Float_t         deltaPhi3_20;
   Float_t         maxDeltaPhi;
   Float_t         maxDeltaPhiAll;
   Float_t         maxDeltaPhiAll30;
   Float_t         maxDeltaPhi20_eta5_noIdAll;
   Float_t         sumDeltaPhi;
   Float_t         diffDeltaPhi;
   Float_t         deltaPhiStar;
   Float_t         deltaPhiStar_badjet_pt;
   Float_t         deltaPhiStar_badjet_estimatedPt;
   Float_t         deltaPhiStar_badjet_genPt;
   Float_t         deltaPhiStar_badjet_energy;
   Float_t         deltaPhiStar_badjet_phi;
   Float_t         deltaPhiStar_badjet_eta;
   Float_t         deltaPhiStar_badjet_CSV;
   Float_t         deltaPhiStar_badjet_chEmE;
   Float_t         deltaPhiStar_badjet_chHadE;
   Float_t         deltaPhiStar_badjet_photonE;
   Float_t         deltaPhiStar_badjet_neuEmE;
   Float_t         deltaPhiStar_badjet_neuHadE;
   Float_t         deltaPhiStar_badjet_chMuE;
   Float_t         deltaPhiStar_badjet_etaetaMoment;
   Float_t         deltaPhiStar_badjet_etaphiMoment;
   Float_t         deltaPhiStar_badjet_phiphiMoment;
   Float_t         deltaPhiStar_badjet_rawPt;
   Float_t         deltaPhiStar_badjet_mass;
   Float_t         deltaPhiStar_badjet_PUbeta;
   Float_t         deltaPhiStar_badjet_PUbetaStar;
   Float_t         deltaPhiStar_badjet_PUbetaStarClassic;
   Int_t           deltaPhiStar_badjet_chMult;
   Int_t           deltaPhiStar_badjet_neuMult;
   Int_t           deltaPhiStar_badjet_muMult;
   Int_t           deltaPhiStar_badjet_n60;
   Int_t           deltaPhiStar_badjet_n90;
   Float_t         hh_mostPUlike_beta;
   Float_t         hh_mostPUlike_betaStar;
   Float_t         minDeltaPhiN;
   Float_t         deltaPhiN1;
   Float_t         deltaPhiN2;
   Float_t         deltaPhiN3;
   UInt_t          minDeltaPhi_chosenJet;
   UInt_t          minDeltaPhiN_chosenJet;
   UInt_t          maxJetMis_chosenJet;
   UInt_t          maxJetFracMis_chosenJet;
   Float_t         minDeltaPhiN_deltaT;
   Float_t         deltaT1;
   Float_t         deltaT2;
   Float_t         deltaT3;
   Float_t         minDeltaPhiN_asin;
   Float_t         deltaPhiN_asin1;
   Float_t         deltaPhiN_asin2;
   Float_t         deltaPhiN_asin3;
   Float_t         CSVout1;
   Float_t         CSVout2;
   Float_t         CSVout3;
   Float_t         CSVbest1;
   Float_t         CSVbest2;
   Float_t         CSVbest3;
   Float_t         CSVbest4;
   Float_t         minDeltaPhiAllb30;
   Float_t         deltaPhib1;
   Float_t         deltaPhib2;
   Float_t         deltaPhib3;
   Float_t         minDeltaPhiMETMuonsAll;
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
   Int_t           jetneutralhadronmult1;
   Int_t           jetmumult1;
   Float_t         alletajetpt1;
   Float_t         alletajetphi1;
   Float_t         alletajeteta1;
   Float_t         alletajetneutralhadronfrac1;
   Float_t         alletajetneutralemfrac1;
   Float_t         alletajetphotonfrac1;
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
   Int_t           jetneutralhadronmult2;
   Int_t           jetmumult2;
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
   Int_t           jetneutralhadronmult3;
   Int_t           jetmumult3;
   Float_t         jetpt4;
   Float_t         jetgenpt4;
   Float_t         jeteta4;
   Float_t         jetgeneta4;
   Float_t         jetphi4;
   Float_t         jetgenphi4;
   Float_t         jetenergy4;
   Int_t           jetflavor4;
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
   Float_t         bjetpt4;
   Float_t         bjeteta4;
   Float_t         bjetphi4;
   Float_t         bjetenergy4;
   Int_t           bjetflavor4;
   Float_t         tauMatch_jetpt;
   Float_t         tauMatch_chhadmult;
   Float_t         tauMatch_jetcsv;
   Float_t         tauMatch_MT;
   Float_t         eleet1;
   Float_t         elephi1;
   Float_t         eleeta1;
   Int_t           elecharge1;
   Float_t         muonpt1;
   Float_t         muonphi1;
   Float_t         muoneta1;
   Int_t           muoncharge1;
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
   Float_t         isotrackpt1;
   Float_t         isotracketa1;
   Float_t         isotrackphi1;
   Float_t         eleRelIso;
   Float_t         muonRelIso;
   Float_t         rl;
   Float_t         rMET;
   Float_t         transverseSphericity_jets;
   Float_t         transverseSphericity_jetsMet;
   Float_t         transverseSphericity_jetsMetLeptons;
   Float_t         transverseSphericity_jetsLeptons;
   Float_t         transverseSphericity_jets30;
   Float_t         transverseSphericity_jets30Met;
   Float_t         transverseSphericity_jets30MetLeptons;
   Float_t         transverseSphericity_jets30Leptons;
   Float_t         minTransverseMETSignificance;
   Float_t         maxTransverseMETSignificance;
   Float_t         transverseMETSignificance1;
   Float_t         transverseMETSignificance2;
   Float_t         transverseMETSignificance3;
   Int_t           nIsoTracks20_005_03;
   Int_t           nIsoTracks15_005_03;
   Int_t           nIsoTracks10_005_03;
   Int_t           nIsoTracks5_005_03;
   Int_t           nIsoTracks15_005_03_lepcleaned;
   Int_t           nIsoPFcands5_010;
   Int_t           nIsoPFcands10_010;
   Int_t           nIsoPFcands15_010;
   Float_t         minTrackIso10;
   Float_t         minTrackIso5;
   Float_t         minTrackIso15;
   Float_t         minTrackIso20;
   Float_t         trackd0_5;
   Float_t         trackd0_10;
   Float_t         trackd0_15;
   Float_t         trackd0_20;
   Float_t         trackpt_5;
   Float_t         trackpt_10;
   Float_t         trackpt_15;
   Float_t         trackpt_20;

   // List of branches
   TBranch        *b_weight;   //!
   TBranch        *b_weight2;   //!
   TBranch        *b_weight3;   //!
   TBranch        *b_scanCrossSection;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_lumiSection;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_m0;   //!
   TBranch        *b_m12;   //!
   TBranch        *b_mIntermediate;   //!
   TBranch        *b_ttbarDecayCode;   //!
   TBranch        *b_lostLeptonCode;   //!
   TBranch        *b_hadWcode;   //!
   TBranch        *b_genTopPtLep;   //!
   TBranch        *b_genTopPtHad;   //!
   TBranch        *b_genWPtLep;   //!
   TBranch        *b_genWPtHad;   //!
   TBranch        *b_minDeltaRLeptonJet;   //!
   TBranch        *b_minDeltaRLeptonJetReco;   //!
   TBranch        *b_minDeltaRJetFlavor;   //!
   TBranch        *b_PUweight;   //!
   TBranch        *b_PUweightSystVar;   //!
   TBranch        *b_hlteff;   //!
   TBranch        *b_pdfWeightsCTEQ;   //!
   TBranch        *b_pdfWeightsMSTW;   //!
   TBranch        *b_pdfWeightsNNPDF;   //!
   TBranch        *b_PIDweight;   //!
   TBranch        *b_BTagWeight;   //!
   TBranch        *b_BTagWeightLMT;   //!
   TBranch        *b_BTagWeightLMT_HFplus;   //!
   TBranch        *b_BTagWeightLMT_HFminus;   //!
   TBranch        *b_BTagWeightLMT_LFplus;   //!
   TBranch        *b_BTagWeightLMT_LFminus;   //!
   TBranch        *b_SUSY_msq12;   //!
   TBranch        *b_SUSY_msq23;   //!
   TBranch        *b_SUSY_gluino_pt;   //!
   TBranch        *b_SUSY_top_pt;   //!
   TBranch        *b_SUSY_topbar_pt;   //!
   TBranch        *b_SUSY_chi0_pt;   //!
   TBranch        *b_topPtWeight;   //!
   TBranch        *b_topPt;   //!
   TBranch        *b_topPtWeightOfficial;   //!
   TBranch        *b_higgsMass;   //!
   TBranch        *b_higgsSumpt;   //!
   TBranch        *b_hhDeltaR;   //!
   TBranch        *b_hhDeltaEta;   //!
   TBranch        *b_hhDeltaPhi;   //!
   TBranch        *b_hhDeltaPhiWrong1;   //!
   TBranch        *b_hhDeltaPhiWrong2;   //!
   TBranch        *b_hhDeltaEtaWrong1;   //!
   TBranch        *b_hhDeltaEtaWrong2;   //!
   TBranch        *b_higgs1CosHelCorrect;   //!
   TBranch        *b_higgs1CosHelWrong1;   //!
   TBranch        *b_higgs1CosHelWrong2;   //!
   TBranch        *b_higgs1b1pt;   //!
   TBranch        *b_higgs1b1phi;   //!
   TBranch        *b_higgs1b1eta;   //!
   TBranch        *b_higgs1b2pt;   //!
   TBranch        *b_higgs1b2phi;   //!
   TBranch        *b_higgs1b2eta;   //!
   TBranch        *b_higgs2b1pt;   //!
   TBranch        *b_higgs2b1phi;   //!
   TBranch        *b_higgs2b1eta;   //!
   TBranch        *b_higgs2b2pt;   //!
   TBranch        *b_higgs2b2phi;   //!
   TBranch        *b_higgs2b2eta;   //!
   TBranch        *b_higgsb1pt;   //!
   TBranch        *b_higgsb2pt;   //!
   TBranch        *b_higgsb3pt;   //!
   TBranch        *b_higgsb4pt;   //!
   TBranch        *b_higgsDeltaRminCorrect;   //!
   TBranch        *b_higgsDeltaRminWrong1;   //!
   TBranch        *b_higgsDeltaRminWrong2;   //!
   TBranch        *b_higgsDeltaRmaxCorrect;   //!
   TBranch        *b_higgsDeltaRmaxWrong1;   //!
   TBranch        *b_higgsDeltaRmaxWrong2;   //!
   TBranch        *b_higgsDeltaPhiRminCorrect;   //!
   TBranch        *b_higgsDeltaPhiRminWrong1;   //!
   TBranch        *b_higgsDeltaPhiRminWrong2;   //!
   TBranch        *b_higgsDeltaPhiRmaxCorrect;   //!
   TBranch        *b_higgsDeltaPhiRmaxWrong1;   //!
   TBranch        *b_higgsDeltaPhiRmaxWrong2;   //!
   TBranch        *b_prob2b;   //!
   TBranch        *b_prob3b;   //!
   TBranch        *b_prob4b;   //!
   TBranch        *b_prob0;   //!
   TBranch        *b_probge1;   //!
   TBranch        *b_prob1;   //!
   TBranch        *b_probge2;   //!
   TBranch        *b_prob2;   //!
   TBranch        *b_probge3;   //!
   TBranch        *b_prob3;   //!
   TBranch        *b_probge4;   //!
   TBranch        *b_prob0_HFplus;   //!
   TBranch        *b_probge1_HFplus;   //!
   TBranch        *b_prob1_HFplus;   //!
   TBranch        *b_probge2_HFplus;   //!
   TBranch        *b_prob2_HFplus;   //!
   TBranch        *b_probge3_HFplus;   //!
   TBranch        *b_prob3_HFplus;   //!
   TBranch        *b_probge4_HFplus;   //!
   TBranch        *b_prob0_HFminus;   //!
   TBranch        *b_probge1_HFminus;   //!
   TBranch        *b_prob1_HFminus;   //!
   TBranch        *b_probge2_HFminus;   //!
   TBranch        *b_prob2_HFminus;   //!
   TBranch        *b_probge3_HFminus;   //!
   TBranch        *b_prob3_HFminus;   //!
   TBranch        *b_probge4_HFminus;   //!
   TBranch        *b_prob0_LFplus;   //!
   TBranch        *b_probge1_LFplus;   //!
   TBranch        *b_prob1_LFplus;   //!
   TBranch        *b_probge2_LFplus;   //!
   TBranch        *b_prob2_LFplus;   //!
   TBranch        *b_probge3_LFplus;   //!
   TBranch        *b_prob3_LFplus;   //!
   TBranch        *b_probge4_LFplus;   //!
   TBranch        *b_prob0_LFminus;   //!
   TBranch        *b_probge1_LFminus;   //!
   TBranch        *b_prob1_LFminus;   //!
   TBranch        *b_probge2_LFminus;   //!
   TBranch        *b_prob2_LFminus;   //!
   TBranch        *b_probge3_LFminus;   //!
   TBranch        *b_prob3_LFminus;   //!
   TBranch        *b_probge4_LFminus;   //!
   TBranch        *b_cutPV;   //!
   TBranch        *b_cutTrigger;   //!
   TBranch        *b_cutTrigger2;   //!
   TBranch        *b_csctighthaloFilter;   //!
   TBranch        *b_eenoiseFilter;   //!
   TBranch        *b_greedymuonFilter;   //!
   TBranch        *b_hbhenoiseFilter;   //!
   TBranch        *b_inconsistentmuonFilter;   //!
   TBranch        *b_ra2ecaltpFilter;   //!
   TBranch        *b_scrapingvetoFilter;   //!
   TBranch        *b_trackingfailureFilter;   //!
   TBranch        *b_badjetFilter;   //!
   TBranch        *b_passCleaning;   //!
   TBranch        *b_PBNRcode;   //!
   TBranch        *b_buggyEvent;   //!
   TBranch        *b_maxTOBTECjetDeltaMult;   //!
   TBranch        *b_TOBTECjetChMult;   //!
   TBranch        *b_nGoodPV;   //!
   TBranch        *b_PV0_x;   //!
   TBranch        *b_PV0_xErr;   //!
   TBranch        *b_PV0_y;   //!
   TBranch        *b_PV0_yErr;   //!
   TBranch        *b_PV0_z;   //!
   TBranch        *b_PV0_zErr;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_xErr;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_yErr;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_zErr;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_isFake;   //!
   TBranch        *b_PV_isGood;   //!
   TBranch        *b_PV_isValid;   //!
   TBranch        *b_PV_tracksSize;   //!
   TBranch        *b_PV_nearestZindex;   //!
   TBranch        *b_rho_kt6PFJetsForIsolation;   //!
   TBranch        *b_BS0_x;   //!
   TBranch        *b_BS0_xErr;   //!
   TBranch        *b_BS0_y;   //!
   TBranch        *b_BS0_yErr;   //!
   TBranch        *b_BS0_z;   //!
   TBranch        *b_BS0_zErr;   //!
   TBranch        *b_SUSY_nb;   //!
   TBranch        *b_SUSY_process;   //!
   TBranch        *b_SUSY_recoilPt;   //!
   TBranch        *b_SUSY_ISRweight;   //!
   TBranch        *b_SUSY_ISRweightSystDown;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_njets30;   //!
   TBranch        *b_njets20;   //!
   TBranch        *b_nPUjets30;   //!
   TBranch        *b_nPUjets20;   //!
   TBranch        *b_njets20_5p0;   //!
   TBranch        *b_njetsHiggsMatch20;   //!
   TBranch        *b_ncompleteHiggsReco20;   //!
   TBranch        *b_njetsHiggsMatch20_eta5;   //!
   TBranch        *b_nbjets;   //!
   TBranch        *b_nbjets30;   //!
   TBranch        *b_nbjetsTweaked;   //!
   TBranch        *b_ntruebjets;   //!
   TBranch        *b_nGluonsSplitToLF;   //!
   TBranch        *b_nGluonsSplitToC;   //!
   TBranch        *b_nGluonsSplitToB;   //!
   TBranch        *b_nlfjetsFromGluons;   //!
   TBranch        *b_ncjetsFromGluons;   //!
   TBranch        *b_nbjetsFromGluons;   //!
   TBranch        *b_nElectrons;   //!
   TBranch        *b_nElectronsBug;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_nElectrons5;   //!
   TBranch        *b_nMuons5;   //!
   TBranch        *b_nElectrons15;   //!
   TBranch        *b_nMuons15;   //!
   TBranch        *b_nTightMuons;   //!
   TBranch        *b_nElectrons20;   //!
   TBranch        *b_nMuons20;   //!
   TBranch        *b_nElectronsNoRelIso;   //!
   TBranch        *b_nMuonsNoRelIso;   //!
   TBranch        *b_nTausVLoose;   //!
   TBranch        *b_nTausLoose;   //!
   TBranch        *b_nTausMedium;   //!
   TBranch        *b_nTausTight;   //!
   TBranch        *b_nCorrectRecoStop;   //!
   TBranch        *b_bestZmass;   //!
   TBranch        *b_mjj1;   //!
   TBranch        *b_mjj2;   //!
   TBranch        *b_mjjdiff;   //!
   TBranch        *b_higgsMbb1;   //!
   TBranch        *b_higgsMbb2;   //!
   TBranch        *b_higgsMjj1;   //!
   TBranch        *b_higgsMjj2;   //!
   TBranch        *b_higgsMbb1delta;   //!
   TBranch        *b_higgsMbb2delta;   //!
   TBranch        *b_higgsMbb1MassDiff;   //!
   TBranch        *b_higgsMbb2MassDiff;   //!
   TBranch        *b_higgsMbbMassDiff_tiebreak;   //!
   TBranch        *b_higgsMbb1MassDiffAll;   //!
   TBranch        *b_higgsMbb2MassDiffAll;   //!
   TBranch        *b_higgsMbb1MassDiff_correct;   //!
   TBranch        *b_higgsMbb1MassDiffAll_correct;   //!
   TBranch        *b_deltaPhi_hh;   //!
   TBranch        *b_deltaRmax_hh;   //!
   TBranch        *b_deltaRmin_hh;   //!
   TBranch        *b_deltaEta_hh;   //!
   TBranch        *b_sumPt_hh;   //!
   TBranch        *b_higgsPt1;   //!
   TBranch        *b_higgsPt2;   //!
   TBranch        *b_deltaRmax_hhAll;   //!
   TBranch        *b_deltaRmin_hhAll;   //!
   TBranch        *b_higgs1jetpt1;   //!
   TBranch        *b_higgs1jetpt2;   //!
   TBranch        *b_higgs2jetpt1;   //!
   TBranch        *b_higgs2jetpt2;   //!
   TBranch        *b_higgs1jeteta1;   //!
   TBranch        *b_higgs1jeteta2;   //!
   TBranch        *b_higgs2jeteta1;   //!
   TBranch        *b_higgs2jeteta2;   //!
   TBranch        *b_higgs1CSV1;   //!
   TBranch        *b_higgs1CSV2;   //!
   TBranch        *b_higgs2CSV1;   //!
   TBranch        *b_higgs2CSV2;   //!
   TBranch        *b_higgs1partonId1;   //!
   TBranch        *b_higgs1partonId2;   //!
   TBranch        *b_higgs2partonId1;   //!
   TBranch        *b_higgs2partonId2;   //!
   TBranch        *b_higgs1partonFlavor1;   //!
   TBranch        *b_higgs1partonFlavor2;   //!
   TBranch        *b_higgs2partonFlavor1;   //!
   TBranch        *b_higgs2partonFlavor2;   //!
   TBranch        *b_higgs1partonMomId1;   //!
   TBranch        *b_higgs1partonMomId2;   //!
   TBranch        *b_higgs2partonMomId1;   //!
   TBranch        *b_higgs2partonMomId2;   //!
   TBranch        *b_higgs1tauMatch1;   //!
   TBranch        *b_higgs1tauMatch2;   //!
   TBranch        *b_higgs2tauMatch1;   //!
   TBranch        *b_higgs2tauMatch2;   //!
   TBranch        *b_higgsWCandMass;   //!
   TBranch        *b_higgsWCandMassAll;   //!
   TBranch        *b_maxDelPhiThrustJet;   //!
   TBranch        *b_minDelPhiThrustMET;   //!
   TBranch        *b_mjjb1;   //!
   TBranch        *b_mjjb2;   //!
   TBranch        *b_topPT1;   //!
   TBranch        *b_topPT2;   //!
   TBranch        *b_nbjetsSSVM;   //!
   TBranch        *b_nbjetsSSVHPT;   //!
   TBranch        *b_nbjetsTCHPT;   //!
   TBranch        *b_nbjetsTCHET;   //!
   TBranch        *b_nbjetsTCHPM;   //!
   TBranch        *b_nbjetsCSVT;   //!
   TBranch        *b_nbjetsCSVM;   //!
   TBranch        *b_nbjetsCSVL;   //!
   TBranch        *b_isRealData;   //!
   TBranch        *b_pass_utilityHLT;   //!
   TBranch        *b_prescaleUtilityHLT;   //!
   TBranch        *b_versionUtilityHLT;   //!
   TBranch        *b_pass_utilityPrescaleModuleHLT;   //!
   TBranch        *b_pass_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45;   //!
   TBranch        *b_pass_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45;   //!
   TBranch        *b_pass_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned;   //!
   TBranch        *b_pass_DiCentralPFJet30_PFMET80;   //!
   TBranch        *b_pass_DiCentralPFJet30_PFMET80_BTagCSV07;   //!
   TBranch        *b_pass_DiCentralPFJet30_PFMHT80;   //!
   TBranch        *b_pass_DiCentralPFJet50_PFMET80;   //!
   TBranch        *b_pass_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80;   //!
   TBranch        *b_pass_DiPFJet80_DiPFJet30_BTagCSVd07d05;   //!
   TBranch        *b_pass_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_pass_HT200;   //!
   TBranch        *b_pass_HT250;   //!
   TBranch        *b_pass_HT250_AlphaT0p55;   //!
   TBranch        *b_pass_HT300;   //!
   TBranch        *b_pass_HT300_AlphaT0p53;   //!
   TBranch        *b_pass_IsoMu24;   //!
   TBranch        *b_pass_IsoMu24_eta2p1;   //!
   TBranch        *b_pass_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV;   //!
   TBranch        *b_pass_L1ETM40;   //!
   TBranch        *b_pass_MET120_HBHENoiseCleaned;   //!
   TBranch        *b_pass_MET200;   //!
   TBranch        *b_pass_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_pass_Mu17_Mu8;   //!
   TBranch        *b_pass_Mu8_DiJet30;   //!
   TBranch        *b_pass_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_pass_PFHT350;   //!
   TBranch        *b_pass_PFHT350_Mu15_PFMET45;   //!
   TBranch        *b_pass_PFHT350_PFMET100;   //!
   TBranch        *b_pass_PFHT650;   //!
   TBranch        *b_pass_PFJet80;   //!
   TBranch        *b_pass_PFMET150;   //!
   TBranch        *b_pass_PFNoPUHT350;   //!
   TBranch        *b_pass_PFNoPUHT350_Mu15_PFMET45;   //!
   TBranch        *b_pass_PFNoPUHT350_PFMET100;   //!
   TBranch        *b_pass_PFNoPUHT650;   //!
   TBranch        *b_pass_Photon135;   //!
   TBranch        *b_pass_Photon150;   //!
   TBranch        *b_pass_QuadJet80;   //!
   TBranch        *b_pass_SixJet45;   //!
   TBranch        *b_passMC_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45;   //!
   TBranch        *b_passMC_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45;   //!
   TBranch        *b_passMC_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned;   //!
   TBranch        *b_passMC_DiCentralPFJet30_PFMET80;   //!
   TBranch        *b_passMC_DiCentralPFJet30_PFMET80_BTagCSV07;   //!
   TBranch        *b_passMC_DiCentralPFJet30_PFMHT80;   //!
   TBranch        *b_passMC_DiCentralPFJet50_PFMET80;   //!
   TBranch        *b_passMC_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80;   //!
   TBranch        *b_passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05;   //!
   TBranch        *b_passMC_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_passMC_HT200;   //!
   TBranch        *b_passMC_HT250;   //!
   TBranch        *b_passMC_HT250_AlphaT0p55;   //!
   TBranch        *b_passMC_HT300;   //!
   TBranch        *b_passMC_HT300_AlphaT0p53;   //!
   TBranch        *b_passMC_IsoMu24;   //!
   TBranch        *b_passMC_IsoMu24_eta2p1;   //!
   TBranch        *b_passMC_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV;   //!
   TBranch        *b_passMC_L1ETM40;   //!
   TBranch        *b_passMC_MET120_HBHENoiseCleaned;   //!
   TBranch        *b_passMC_MET200;   //!
   TBranch        *b_passMC_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_passMC_Mu17_Mu8;   //!
   TBranch        *b_passMC_Mu8_DiJet30;   //!
   TBranch        *b_passMC_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_passMC_PFHT350;   //!
   TBranch        *b_passMC_PFHT350_Mu15_PFMET45;   //!
   TBranch        *b_passMC_PFHT350_PFMET100;   //!
   TBranch        *b_passMC_PFHT650;   //!
   TBranch        *b_passMC_PFJet80;   //!
   TBranch        *b_passMC_PFMET150;   //!
   TBranch        *b_passMC_PFNoPUHT350;   //!
   TBranch        *b_passMC_PFNoPUHT350_Mu15_PFMET45;   //!
   TBranch        *b_passMC_PFNoPUHT350_PFMET100;   //!
   TBranch        *b_passMC_PFNoPUHT650;   //!
   TBranch        *b_passMC_Photon135;   //!
   TBranch        *b_passMC_Photon150;   //!
   TBranch        *b_passMC_QuadJet80;   //!
   TBranch        *b_passMC_SixJet45;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_HT30;   //!
   TBranch        *b_HT20;   //!
   TBranch        *b_ST;   //!
   TBranch        *b_STeff;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METphi;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_METsig;   //!
   TBranch        *b_METsig00;   //!
   TBranch        *b_METsig10;   //!
   TBranch        *b_METsig11;   //!
   TBranch        *b_METsig_2012;   //!
   TBranch        *b_METsig00_2012;   //!
   TBranch        *b_METsig10_2012;   //!
   TBranch        *b_METsig11_2012;   //!
   TBranch        *b_caloMET;   //!
   TBranch        *b_caloMETphi;   //!
   TBranch        *b_rawPFMET;   //!
   TBranch        *b_rawPFMETphi;   //!
   TBranch        *b_caloMETwithHO;   //!
   TBranch        *b_caloMETwithHOphi;   //!
   TBranch        *b_unclusteredMET;   //!
   TBranch        *b_unclusteredMETphi;   //!
   TBranch        *b_bestWMass;   //!
   TBranch        *b_bestTopMass;   //!
   TBranch        *b_topCosHel;   //!
   TBranch        *b_WCosHel;   //!
   TBranch        *b_MT_b;   //!
   TBranch        *b_MT_Hbb1;   //!
   TBranch        *b_MT_Hbb2;   //!
   TBranch        *b_MT_Wlep;   //!
   TBranch        *b_MT_Wlep5;   //!
   TBranch        *b_MT_Wlep15;   //!
   TBranch        *b_maxDeltaPhi_bbbb;   //!
   TBranch        *b_maxDeltaPhi_bb_bb;   //!
   TBranch        *b_MT_bestCSV_gencode;   //!
   TBranch        *b_MT_bestCSV;   //!
   TBranch        *b_MT_jim;   //!
   TBranch        *b_minMT_jetMET;   //!
   TBranch        *b_minMT_bMET;   //!
   TBranch        *b_maxMT_bMET;   //!
   TBranch        *b_deltaThetaT;   //!
   TBranch        *b_deltaR_bestTwoCSV;   //!
   TBranch        *b_deltaPhi_bestTwoCSV;   //!
   TBranch        *b_mjj_bestTwoCSV;   //!
   TBranch        *b_deltaR_closestB;   //!
   TBranch        *b_mjj_closestB;   //!
   TBranch        *b_mjj_closestB20;   //!
   TBranch        *b_mjj_closestB20_correct;   //!
   TBranch        *b_deltaR_closestB20;   //!
   TBranch        *b_cosHel_closestB20;   //!
   TBranch        *b_deltaPhi_hb20;   //!
   TBranch        *b_sumPtMjjDiff_closestB20;   //!
   TBranch        *b_mjj_h125;   //!
   TBranch        *b_mjj_highptCSVM;   //!
   TBranch        *b_mjj_highptCSVT;   //!
   TBranch        *b_minDeltaPhi;   //!
   TBranch        *b_minDeltaPhiAll;   //!
   TBranch        *b_minDeltaPhiAll20;   //!
   TBranch        *b_minDeltaPhi20;   //!
   TBranch        *b_minDeltaPhi30;   //!
   TBranch        *b_minDeltaPhi20_eta5_noIdAll;   //!
   TBranch        *b_minDeltaPhi20_eta5_noIdAll_nobeta;   //!
   TBranch        *b_deltaPhi1;   //!
   TBranch        *b_deltaPhi2;   //!
   TBranch        *b_deltaPhi3;   //!
   TBranch        *b_deltaPhi1_20;   //!
   TBranch        *b_deltaPhi2_20;   //!
   TBranch        *b_deltaPhi3_20;   //!
   TBranch        *b_maxDeltaPhi;   //!
   TBranch        *b_maxDeltaPhiAll;   //!
   TBranch        *b_maxDeltaPhiAll30;   //!
   TBranch        *b_maxDeltaPhi20_eta5_noIdAll;   //!
   TBranch        *b_sumDeltaPhi;   //!
   TBranch        *b_diffDeltaPhi;   //!
   TBranch        *b_deltaPhiStar;   //!
   TBranch        *b_deltaPhiStar_badjet_pt;   //!
   TBranch        *b_deltaPhiStar_badjet_estimatedPt;   //!
   TBranch        *b_deltaPhiStar_badjet_genPt;   //!
   TBranch        *b_deltaPhiStar_badjet_energy;   //!
   TBranch        *b_deltaPhiStar_badjet_phi;   //!
   TBranch        *b_deltaPhiStar_badjet_eta;   //!
   TBranch        *b_deltaPhiStar_badjet_CSV;   //!
   TBranch        *b_deltaPhiStar_badjet_chEmE;   //!
   TBranch        *b_deltaPhiStar_badjet_chHadE;   //!
   TBranch        *b_deltaPhiStar_badjet_photonE;   //!
   TBranch        *b_deltaPhiStar_badjet_neuEmE;   //!
   TBranch        *b_deltaPhiStar_badjet_neuHadE;   //!
   TBranch        *b_deltaPhiStar_badjet_chMuE;   //!
   TBranch        *b_deltaPhiStar_badjet_etaetaMoment;   //!
   TBranch        *b_deltaPhiStar_badjet_etaphiMoment;   //!
   TBranch        *b_deltaPhiStar_badjet_phiphiMoment;   //!
   TBranch        *b_deltaPhiStar_badjet_rawPt;   //!
   TBranch        *b_deltaPhiStar_badjet_mass;   //!
   TBranch        *b_deltaPhiStar_badjet_PUbeta;   //!
   TBranch        *b_deltaPhiStar_badjet_PUbetaStar;   //!
   TBranch        *b_deltaPhiStar_badjet_PUbetaStarClassic;   //!
   TBranch        *b_deltaPhiStar_badjet_chMult;   //!
   TBranch        *b_deltaPhiStar_badjet_neuMult;   //!
   TBranch        *b_deltaPhiStar_badjet_muMult;   //!
   TBranch        *b_deltaPhiStar_badjet_n60;   //!
   TBranch        *b_deltaPhiStar_badjet_n90;   //!
   TBranch        *b_hh_mostPUlike_beta;   //!
   TBranch        *b_hh_mostPUlike_betaStar;   //!
   TBranch        *b_minDeltaPhiN;   //!
   TBranch        *b_deltaPhiN1;   //!
   TBranch        *b_deltaPhiN2;   //!
   TBranch        *b_deltaPhiN3;   //!
   TBranch        *b_minDeltaPhi_chosenJet;   //!
   TBranch        *b_minDeltaPhiN_chosenJet;   //!
   TBranch        *b_maxJetMis_chosenJet;   //!
   TBranch        *b_maxJetFracMis_chosenJet;   //!
   TBranch        *b_minDeltaPhiN_deltaT;   //!
   TBranch        *b_deltaT1;   //!
   TBranch        *b_deltaT2;   //!
   TBranch        *b_deltaT3;   //!
   TBranch        *b_minDeltaPhiN_asin;   //!
   TBranch        *b_deltaPhiN_asin1;   //!
   TBranch        *b_deltaPhiN_asin2;   //!
   TBranch        *b_deltaPhiN_asin3;   //!
   TBranch        *b_CSVout1;   //!
   TBranch        *b_CSVout2;   //!
   TBranch        *b_CSVout3;   //!
   TBranch        *b_CSVbest1;   //!
   TBranch        *b_CSVbest2;   //!
   TBranch        *b_CSVbest3;   //!
   TBranch        *b_CSVbest4;   //!
   TBranch        *b_minDeltaPhiAllb30;   //!
   TBranch        *b_deltaPhib1;   //!
   TBranch        *b_deltaPhib2;   //!
   TBranch        *b_deltaPhib3;   //!
   TBranch        *b_minDeltaPhiMETMuonsAll;   //!
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
   TBranch        *b_jetneutralhadronmult1;   //!
   TBranch        *b_jetmumult1;   //!
   TBranch        *b_alletajetpt1;   //!
   TBranch        *b_alletajetphi1;   //!
   TBranch        *b_alletajeteta1;   //!
   TBranch        *b_alletajetneutralhadronfrac1;   //!
   TBranch        *b_alletajetneutralemfrac1;   //!
   TBranch        *b_alletajetphotonfrac1;   //!
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
   TBranch        *b_jetneutralhadronmult2;   //!
   TBranch        *b_jetmumult2;   //!
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
   TBranch        *b_jetneutralhadronmult3;   //!
   TBranch        *b_jetmumult3;   //!
   TBranch        *b_jetpt4;   //!
   TBranch        *b_jetgenpt4;   //!
   TBranch        *b_jeteta4;   //!
   TBranch        *b_jetgeneta4;   //!
   TBranch        *b_jetphi4;   //!
   TBranch        *b_jetgenphi4;   //!
   TBranch        *b_jetenergy4;   //!
   TBranch        *b_jetflavor4;   //!
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
   TBranch        *b_bjetpt4;   //!
   TBranch        *b_bjeteta4;   //!
   TBranch        *b_bjetphi4;   //!
   TBranch        *b_bjetenergy4;   //!
   TBranch        *b_bjetflavor4;   //!
   TBranch        *b_tauMatch_jetpt;   //!
   TBranch        *b_tauMatch_chhadmult;   //!
   TBranch        *b_tauMatch_jetcsv;   //!
   TBranch        *b_tauMatch_MT;   //!
   TBranch        *b_eleet1;   //!
   TBranch        *b_elephi1;   //!
   TBranch        *b_eleeta1;   //!
   TBranch        *b_elecharge1;   //!
   TBranch        *b_muonpt1;   //!
   TBranch        *b_muonphi1;   //!
   TBranch        *b_muoneta1;   //!
   TBranch        *b_muoncharge1;   //!
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
   TBranch        *b_isotrackpt1;   //!
   TBranch        *b_isotracketa1;   //!
   TBranch        *b_isotrackphi1;   //!
   TBranch        *b_eleRelIso;   //!
   TBranch        *b_muonRelIso;   //!
   TBranch        *b_rl;   //!
   TBranch        *b_rMET;   //!
   TBranch        *b_transverseSphericity_jets;   //!
   TBranch        *b_transverseSphericity_jetsMet;   //!
   TBranch        *b_transverseSphericity_jetsMetLeptons;   //!
   TBranch        *b_transverseSphericity_jetsLeptons;   //!
   TBranch        *b_transverseSphericity_jets30;   //!
   TBranch        *b_transverseSphericity_jets30Met;   //!
   TBranch        *b_transverseSphericity_jets30MetLeptons;   //!
   TBranch        *b_transverseSphericity_jets30Leptons;   //!
   TBranch        *b_minTransverseMETSignificance;   //!
   TBranch        *b_maxTransverseMETSignificance;   //!
   TBranch        *b_transverseMETSignificance1;   //!
   TBranch        *b_transverseMETSignificance2;   //!
   TBranch        *b_transverseMETSignificance3;   //!
   TBranch        *b_nIsoTracks20_005_03;   //!
   TBranch        *b_nIsoTracks15_005_03;   //!
   TBranch        *b_nIsoTracks10_005_03;   //!
   TBranch        *b_nIsoTracks5_005_03;   //!
   TBranch        *b_nIsoTracks15_005_03_lepcleaned;   //!
   TBranch        *b_nIsoPFcands5_010;   //!
   TBranch        *b_nIsoPFcands10_010;   //!
   TBranch        *b_nIsoPFcands15_010;   //!
   TBranch        *b_minTrackIso10;   //!
   TBranch        *b_minTrackIso5;   //!
   TBranch        *b_minTrackIso15;   //!
   TBranch        *b_minTrackIso20;   //!
   TBranch        *b_trackd0_5;   //!
   TBranch        *b_trackd0_10;   //!
   TBranch        *b_trackd0_15;   //!
   TBranch        *b_trackd0_20;   //!
   TBranch        *b_trackpt_5;   //!
   TBranch        *b_trackpt_10;   //!
   TBranch        *b_trackpt_15;   //!
   TBranch        *b_trackpt_20;   //!

   signalEff_hbbhbb(TString path, TString filestub, bool joinbtagbins, bool usebtagsf, bool dopdfs, bool pusyst, int isrmode=99); //BEGIN END ra2b-jmt
   virtual ~signalEff_hbbhbb();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   //BEGIN ra2b-jmt
   TH2D* scanSMSngen_;
   std::vector<TH2D*> scanProcessTotals_;
   TString filestub_;
   bool joinbtagbins_;
   bool usebtagsf_;
   bool dopdfs_;
   bool pusyst_;
   IsrMode theIsrMode_;
   BTagMode theBTagMode_;
   //END ra2b-jmt

};

#endif

#ifdef signalEff_hbbhbb_cxx
//BEGIN ra2b-jmt mod
signalEff_hbbhbb::signalEff_hbbhbb(TString path, TString filestub,bool joinbtagbins, bool usebtagsf, bool dopdfs, bool pusyst, int isrmode) : fChain(0) 
																	    ,scanSMSngen_(0)
																	    ,filestub_(filestub)
																	    ,joinbtagbins_(joinbtagbins)
																	    ,usebtagsf_(usebtagsf)
																	    ,dopdfs_(dopdfs)
																	    ,pusyst_(pusyst)
																	    ,theIsrMode_(kNoIsrWeight)
																	    ,theBTagMode_(kBTag0)
				  //END ra2b-jmt mod
{
 //BEGIN ra2b-jmt
  if (isrmode == -1)      theIsrMode_=kIsrDown;
  else if (isrmode ==  1) theIsrMode_=kIsrUp;
  else if (isrmode ==  0) theIsrMode_=kIsr0;
  else theIsrMode_=kNoIsrWeight;

  TTree* tree=0;
  TString filename = path;
  filename+="reducedTree.";
  filename+=filestub;
  filename+=".root";

  std::cout<<filename<<std::endl;
  std::cout<<" join b = "<<joinbtagbins<<std::endl;
  std::cout<<" b tag sf = "<<usebtagsf<<std::endl;
  std::cout<<" pdfs = "<<dopdfs<<std::endl;
  std::cout<<" isrmode = "<<int(theIsrMode_)<<std::endl;
  //END ra2b-jmt


// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename); //BEGIN END ra2b-jmt
      if (!f || !f->IsOpen()) {
	f = new TFile(filename); //BEGIN END ra2b-jmt
      }
      f->GetObject("reducedTree",tree);
      scanSMSngen_ = (TH2D*) f->Get("scanSMSngen"); //BEGIN  ra2b-jmt mod

      assert(scanSMSngen_!=0); // something went wrong in production of the reducedTree

      //jmt -- i don't understand what point this code has
      //Probably should clone scanSMSngen into member variable
      /*
      for ( int ix = 1; ix<=scanSMSngen_->GetXaxis()->GetNbins(); ix++) {
        for ( int iy = 1; iy<=scanSMSngen_->GetYaxis()->GetNbins() ; iy++) {
          double val=scanSMSngen_->GetBinContent(ix,iy);

          scanSMSngen_->SetBinContent(ix,iy, val);
	  
	}
      }
      */
      //need to loop over the keys in the file and load all scan process totals histos
      //(this is now deprecated...but leave it in place just in case)
      for (int ih = 0; ih<f->GetListOfKeys()->GetEntries(); ih++) {
	TString histname = f->GetListOfKeys()->At(ih)->GetName();
	if (!histname.BeginsWith("scanProcessTotals")) continue;
	scanProcessTotals_.push_back( (TH2D*) f->Get(histname));
      }

      //END ra2b-jmt mod

   }
   Init(tree);
}

signalEff_hbbhbb::~signalEff_hbbhbb()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t signalEff_hbbhbb::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t signalEff_hbbhbb::LoadTree(Long64_t entry)
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

void signalEff_hbbhbb::Init(TTree *tree)
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
   fChain->SetBranchAddress("weight2", &weight2, &b_weight2);
   fChain->SetBranchAddress("weight3", &weight3, &b_weight3);
   fChain->SetBranchAddress("scanCrossSection", &scanCrossSection, &b_scanCrossSection);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("m0", &m0, &b_m0);
   fChain->SetBranchAddress("m12", &m12, &b_m12);
   fChain->SetBranchAddress("mIntermediate", &mIntermediate, &b_mIntermediate);
   fChain->SetBranchAddress("ttbarDecayCode", &ttbarDecayCode, &b_ttbarDecayCode);
   fChain->SetBranchAddress("lostLeptonCode", &lostLeptonCode, &b_lostLeptonCode);
   fChain->SetBranchAddress("hadWcode", &hadWcode, &b_hadWcode);
   fChain->SetBranchAddress("genTopPtLep", &genTopPtLep, &b_genTopPtLep);
   fChain->SetBranchAddress("genTopPtHad", &genTopPtHad, &b_genTopPtHad);
   fChain->SetBranchAddress("genWPtLep", &genWPtLep, &b_genWPtLep);
   fChain->SetBranchAddress("genWPtHad", &genWPtHad, &b_genWPtHad);
   fChain->SetBranchAddress("minDeltaRLeptonJet", &minDeltaRLeptonJet, &b_minDeltaRLeptonJet);
   fChain->SetBranchAddress("minDeltaRLeptonJetReco", &minDeltaRLeptonJetReco, &b_minDeltaRLeptonJetReco);
   fChain->SetBranchAddress("minDeltaRJetFlavor", &minDeltaRJetFlavor, &b_minDeltaRJetFlavor);
   fChain->SetBranchAddress("PUweight", &PUweight, &b_PUweight);
   fChain->SetBranchAddress("PUweightSystVar", &PUweightSystVar, &b_PUweightSystVar);
   fChain->SetBranchAddress("hlteff", &hlteff, &b_hlteff);
   fChain->SetBranchAddress("pdfWeightsCTEQ", pdfWeightsCTEQ, &b_pdfWeightsCTEQ);
   fChain->SetBranchAddress("pdfWeightsMSTW", pdfWeightsMSTW, &b_pdfWeightsMSTW);
   fChain->SetBranchAddress("pdfWeightsNNPDF", pdfWeightsNNPDF, &b_pdfWeightsNNPDF);
   fChain->SetBranchAddress("PIDweight", &PIDweight, &b_PIDweight);
   fChain->SetBranchAddress("BTagWeight", &BTagWeight, &b_BTagWeight);
   fChain->SetBranchAddress("BTagWeightLMT", &BTagWeightLMT, &b_BTagWeightLMT);
   fChain->SetBranchAddress("BTagWeightLMT_HFplus", &BTagWeightLMT_HFplus, &b_BTagWeightLMT_HFplus);
   fChain->SetBranchAddress("BTagWeightLMT_HFminus", &BTagWeightLMT_HFminus, &b_BTagWeightLMT_HFminus);
   fChain->SetBranchAddress("BTagWeightLMT_LFplus", &BTagWeightLMT_LFplus, &b_BTagWeightLMT_LFplus);
   fChain->SetBranchAddress("BTagWeightLMT_LFminus", &BTagWeightLMT_LFminus, &b_BTagWeightLMT_LFminus);
   fChain->SetBranchAddress("SUSY_msq12", SUSY_msq12, &b_SUSY_msq12);
   fChain->SetBranchAddress("SUSY_msq23", SUSY_msq23, &b_SUSY_msq23);
   fChain->SetBranchAddress("SUSY_gluino_pt", SUSY_gluino_pt, &b_SUSY_gluino_pt);
   fChain->SetBranchAddress("SUSY_top_pt", SUSY_top_pt, &b_SUSY_top_pt);
   fChain->SetBranchAddress("SUSY_topbar_pt", SUSY_topbar_pt, &b_SUSY_topbar_pt);
   fChain->SetBranchAddress("SUSY_chi0_pt", SUSY_chi0_pt, &b_SUSY_chi0_pt);
   fChain->SetBranchAddress("topPtWeight", &topPtWeight, &b_topPtWeight);
   fChain->SetBranchAddress("topPt", &topPt, &b_topPt);
   fChain->SetBranchAddress("topPtWeightOfficial", &topPtWeightOfficial, &b_topPtWeightOfficial);
   fChain->SetBranchAddress("higgsMass", higgsMass, &b_higgsMass);
   fChain->SetBranchAddress("higgsSumpt", higgsSumpt, &b_higgsSumpt);
   fChain->SetBranchAddress("hhDeltaR", &hhDeltaR, &b_hhDeltaR);
   fChain->SetBranchAddress("hhDeltaEta", &hhDeltaEta, &b_hhDeltaEta);
   fChain->SetBranchAddress("hhDeltaPhi", &hhDeltaPhi, &b_hhDeltaPhi);
   fChain->SetBranchAddress("hhDeltaPhiWrong1", &hhDeltaPhiWrong1, &b_hhDeltaPhiWrong1);
   fChain->SetBranchAddress("hhDeltaPhiWrong2", &hhDeltaPhiWrong2, &b_hhDeltaPhiWrong2);
   fChain->SetBranchAddress("hhDeltaEtaWrong1", &hhDeltaEtaWrong1, &b_hhDeltaEtaWrong1);
   fChain->SetBranchAddress("hhDeltaEtaWrong2", &hhDeltaEtaWrong2, &b_hhDeltaEtaWrong2);
   fChain->SetBranchAddress("higgs1CosHelCorrect", &higgs1CosHelCorrect, &b_higgs1CosHelCorrect);
   fChain->SetBranchAddress("higgs1CosHelWrong1", &higgs1CosHelWrong1, &b_higgs1CosHelWrong1);
   fChain->SetBranchAddress("higgs1CosHelWrong2", &higgs1CosHelWrong2, &b_higgs1CosHelWrong2);
   fChain->SetBranchAddress("higgs1b1pt", &higgs1b1pt, &b_higgs1b1pt);
   fChain->SetBranchAddress("higgs1b1phi", &higgs1b1phi, &b_higgs1b1phi);
   fChain->SetBranchAddress("higgs1b1eta", &higgs1b1eta, &b_higgs1b1eta);
   fChain->SetBranchAddress("higgs1b2pt", &higgs1b2pt, &b_higgs1b2pt);
   fChain->SetBranchAddress("higgs1b2phi", &higgs1b2phi, &b_higgs1b2phi);
   fChain->SetBranchAddress("higgs1b2eta", &higgs1b2eta, &b_higgs1b2eta);
   fChain->SetBranchAddress("higgs2b1pt", &higgs2b1pt, &b_higgs2b1pt);
   fChain->SetBranchAddress("higgs2b1phi", &higgs2b1phi, &b_higgs2b1phi);
   fChain->SetBranchAddress("higgs2b1eta", &higgs2b1eta, &b_higgs2b1eta);
   fChain->SetBranchAddress("higgs2b2pt", &higgs2b2pt, &b_higgs2b2pt);
   fChain->SetBranchAddress("higgs2b2phi", &higgs2b2phi, &b_higgs2b2phi);
   fChain->SetBranchAddress("higgs2b2eta", &higgs2b2eta, &b_higgs2b2eta);
   fChain->SetBranchAddress("higgsb1pt", &higgsb1pt, &b_higgsb1pt);
   fChain->SetBranchAddress("higgsb2pt", &higgsb2pt, &b_higgsb2pt);
   fChain->SetBranchAddress("higgsb3pt", &higgsb3pt, &b_higgsb3pt);
   fChain->SetBranchAddress("higgsb4pt", &higgsb4pt, &b_higgsb4pt);
   fChain->SetBranchAddress("higgsDeltaRminCorrect", &higgsDeltaRminCorrect, &b_higgsDeltaRminCorrect);
   fChain->SetBranchAddress("higgsDeltaRminWrong1", &higgsDeltaRminWrong1, &b_higgsDeltaRminWrong1);
   fChain->SetBranchAddress("higgsDeltaRminWrong2", &higgsDeltaRminWrong2, &b_higgsDeltaRminWrong2);
   fChain->SetBranchAddress("higgsDeltaRmaxCorrect", &higgsDeltaRmaxCorrect, &b_higgsDeltaRmaxCorrect);
   fChain->SetBranchAddress("higgsDeltaRmaxWrong1", &higgsDeltaRmaxWrong1, &b_higgsDeltaRmaxWrong1);
   fChain->SetBranchAddress("higgsDeltaRmaxWrong2", &higgsDeltaRmaxWrong2, &b_higgsDeltaRmaxWrong2);
   fChain->SetBranchAddress("higgsDeltaPhiRminCorrect", &higgsDeltaPhiRminCorrect, &b_higgsDeltaPhiRminCorrect);
   fChain->SetBranchAddress("higgsDeltaPhiRminWrong1", &higgsDeltaPhiRminWrong1, &b_higgsDeltaPhiRminWrong1);
   fChain->SetBranchAddress("higgsDeltaPhiRminWrong2", &higgsDeltaPhiRminWrong2, &b_higgsDeltaPhiRminWrong2);
   fChain->SetBranchAddress("higgsDeltaPhiRmaxCorrect", &higgsDeltaPhiRmaxCorrect, &b_higgsDeltaPhiRmaxCorrect);
   fChain->SetBranchAddress("higgsDeltaPhiRmaxWrong1", &higgsDeltaPhiRmaxWrong1, &b_higgsDeltaPhiRmaxWrong1);
   fChain->SetBranchAddress("higgsDeltaPhiRmaxWrong2", &higgsDeltaPhiRmaxWrong2, &b_higgsDeltaPhiRmaxWrong2);
   fChain->SetBranchAddress("prob2b", &prob2b, &b_prob2b);
   fChain->SetBranchAddress("prob3b", &prob3b, &b_prob3b);
   fChain->SetBranchAddress("prob4b", &prob4b, &b_prob4b);
   fChain->SetBranchAddress("prob0", &prob0, &b_prob0);
   fChain->SetBranchAddress("probge1", &probge1, &b_probge1);
   fChain->SetBranchAddress("prob1", &prob1, &b_prob1);
   fChain->SetBranchAddress("probge2", &probge2, &b_probge2);
   fChain->SetBranchAddress("prob2", &prob2, &b_prob2);
   fChain->SetBranchAddress("probge3", &probge3, &b_probge3);
   fChain->SetBranchAddress("prob3", &prob3, &b_prob3);
   fChain->SetBranchAddress("probge4", &probge4, &b_probge4);
   fChain->SetBranchAddress("prob0_HFplus", &prob0_HFplus, &b_prob0_HFplus);
   fChain->SetBranchAddress("probge1_HFplus", &probge1_HFplus, &b_probge1_HFplus);
   fChain->SetBranchAddress("prob1_HFplus", &prob1_HFplus, &b_prob1_HFplus);
   fChain->SetBranchAddress("probge2_HFplus", &probge2_HFplus, &b_probge2_HFplus);
   fChain->SetBranchAddress("prob2_HFplus", &prob2_HFplus, &b_prob2_HFplus);
   fChain->SetBranchAddress("probge3_HFplus", &probge3_HFplus, &b_probge3_HFplus);
   fChain->SetBranchAddress("prob3_HFplus", &prob3_HFplus, &b_prob3_HFplus);
   fChain->SetBranchAddress("probge4_HFplus", &probge4_HFplus, &b_probge4_HFplus);
   fChain->SetBranchAddress("prob0_HFminus", &prob0_HFminus, &b_prob0_HFminus);
   fChain->SetBranchAddress("probge1_HFminus", &probge1_HFminus, &b_probge1_HFminus);
   fChain->SetBranchAddress("prob1_HFminus", &prob1_HFminus, &b_prob1_HFminus);
   fChain->SetBranchAddress("probge2_HFminus", &probge2_HFminus, &b_probge2_HFminus);
   fChain->SetBranchAddress("prob2_HFminus", &prob2_HFminus, &b_prob2_HFminus);
   fChain->SetBranchAddress("probge3_HFminus", &probge3_HFminus, &b_probge3_HFminus);
   fChain->SetBranchAddress("prob3_HFminus", &prob3_HFminus, &b_prob3_HFminus);
   fChain->SetBranchAddress("probge4_HFminus", &probge4_HFminus, &b_probge4_HFminus);
   fChain->SetBranchAddress("prob0_LFplus", &prob0_LFplus, &b_prob0_LFplus);
   fChain->SetBranchAddress("probge1_LFplus", &probge1_LFplus, &b_probge1_LFplus);
   fChain->SetBranchAddress("prob1_LFplus", &prob1_LFplus, &b_prob1_LFplus);
   fChain->SetBranchAddress("probge2_LFplus", &probge2_LFplus, &b_probge2_LFplus);
   fChain->SetBranchAddress("prob2_LFplus", &prob2_LFplus, &b_prob2_LFplus);
   fChain->SetBranchAddress("probge3_LFplus", &probge3_LFplus, &b_probge3_LFplus);
   fChain->SetBranchAddress("prob3_LFplus", &prob3_LFplus, &b_prob3_LFplus);
   fChain->SetBranchAddress("probge4_LFplus", &probge4_LFplus, &b_probge4_LFplus);
   fChain->SetBranchAddress("prob0_LFminus", &prob0_LFminus, &b_prob0_LFminus);
   fChain->SetBranchAddress("probge1_LFminus", &probge1_LFminus, &b_probge1_LFminus);
   fChain->SetBranchAddress("prob1_LFminus", &prob1_LFminus, &b_prob1_LFminus);
   fChain->SetBranchAddress("probge2_LFminus", &probge2_LFminus, &b_probge2_LFminus);
   fChain->SetBranchAddress("prob2_LFminus", &prob2_LFminus, &b_prob2_LFminus);
   fChain->SetBranchAddress("probge3_LFminus", &probge3_LFminus, &b_probge3_LFminus);
   fChain->SetBranchAddress("prob3_LFminus", &prob3_LFminus, &b_prob3_LFminus);
   fChain->SetBranchAddress("probge4_LFminus", &probge4_LFminus, &b_probge4_LFminus);
   fChain->SetBranchAddress("cutPV", &cutPV, &b_cutPV);
   fChain->SetBranchAddress("cutTrigger", &cutTrigger, &b_cutTrigger);
   fChain->SetBranchAddress("cutTrigger2", &cutTrigger2, &b_cutTrigger2);
   fChain->SetBranchAddress("csctighthaloFilter", &csctighthaloFilter, &b_csctighthaloFilter);
   fChain->SetBranchAddress("eenoiseFilter", &eenoiseFilter, &b_eenoiseFilter);
   fChain->SetBranchAddress("greedymuonFilter", &greedymuonFilter, &b_greedymuonFilter);
   fChain->SetBranchAddress("hbhenoiseFilter", &hbhenoiseFilter, &b_hbhenoiseFilter);
   fChain->SetBranchAddress("inconsistentmuonFilter", &inconsistentmuonFilter, &b_inconsistentmuonFilter);
   fChain->SetBranchAddress("ra2ecaltpFilter", &ra2ecaltpFilter, &b_ra2ecaltpFilter);
   fChain->SetBranchAddress("scrapingvetoFilter", &scrapingvetoFilter, &b_scrapingvetoFilter);
   fChain->SetBranchAddress("trackingfailureFilter", &trackingfailureFilter, &b_trackingfailureFilter);
   fChain->SetBranchAddress("badjetFilter", &badjetFilter, &b_badjetFilter);
   fChain->SetBranchAddress("passCleaning", &passCleaning, &b_passCleaning);
   fChain->SetBranchAddress("PBNRcode", &PBNRcode, &b_PBNRcode);
   fChain->SetBranchAddress("buggyEvent", &buggyEvent, &b_buggyEvent);
   fChain->SetBranchAddress("maxTOBTECjetDeltaMult", &maxTOBTECjetDeltaMult, &b_maxTOBTECjetDeltaMult);
   fChain->SetBranchAddress("TOBTECjetChMult", &TOBTECjetChMult, &b_TOBTECjetChMult);
   fChain->SetBranchAddress("nGoodPV", &nGoodPV, &b_nGoodPV);
   fChain->SetBranchAddress("PV0_x", &PV0_x, &b_PV0_x);
   fChain->SetBranchAddress("PV0_xErr", &PV0_xErr, &b_PV0_xErr);
   fChain->SetBranchAddress("PV0_y", &PV0_y, &b_PV0_y);
   fChain->SetBranchAddress("PV0_yErr", &PV0_yErr, &b_PV0_yErr);
   fChain->SetBranchAddress("PV0_z", &PV0_z, &b_PV0_z);
   fChain->SetBranchAddress("PV0_zErr", &PV0_zErr, &b_PV0_zErr);
   fChain->SetBranchAddress("PV_x", PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_xErr", PV_xErr, &b_PV_xErr);
   fChain->SetBranchAddress("PV_y", PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_yErr", PV_yErr, &b_PV_yErr);
   fChain->SetBranchAddress("PV_z", PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_zErr", PV_zErr, &b_PV_zErr);
   fChain->SetBranchAddress("PV_chi2", PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_ndof", PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_isFake", PV_isFake, &b_PV_isFake);
   fChain->SetBranchAddress("PV_isGood", PV_isGood, &b_PV_isGood);
   fChain->SetBranchAddress("PV_isValid", PV_isValid, &b_PV_isValid);
   fChain->SetBranchAddress("PV_tracksSize", PV_tracksSize, &b_PV_tracksSize);
   fChain->SetBranchAddress("PV_nearestZindex", &PV_nearestZindex, &b_PV_nearestZindex);
   fChain->SetBranchAddress("rho_kt6PFJetsForIsolation", &rho_kt6PFJetsForIsolation, &b_rho_kt6PFJetsForIsolation);
   fChain->SetBranchAddress("BS0_x", &BS0_x, &b_BS0_x);
   fChain->SetBranchAddress("BS0_xErr", &BS0_xErr, &b_BS0_xErr);
   fChain->SetBranchAddress("BS0_y", &BS0_y, &b_BS0_y);
   fChain->SetBranchAddress("BS0_yErr", &BS0_yErr, &b_BS0_yErr);
   fChain->SetBranchAddress("BS0_z", &BS0_z, &b_BS0_z);
   fChain->SetBranchAddress("BS0_zErr", &BS0_zErr, &b_BS0_zErr);
   fChain->SetBranchAddress("SUSY_nb", &SUSY_nb, &b_SUSY_nb);
   fChain->SetBranchAddress("SUSY_process", &SUSY_process, &b_SUSY_process);
   fChain->SetBranchAddress("SUSY_recoilPt", &SUSY_recoilPt, &b_SUSY_recoilPt);
   fChain->SetBranchAddress("SUSY_ISRweight", &SUSY_ISRweight, &b_SUSY_ISRweight);
   fChain->SetBranchAddress("SUSY_ISRweightSystDown", &SUSY_ISRweightSystDown, &b_SUSY_ISRweightSystDown);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("njets30", &njets30, &b_njets30);
   fChain->SetBranchAddress("njets20", &njets20, &b_njets20);
   fChain->SetBranchAddress("nPUjets30", &nPUjets30, &b_nPUjets30);
   fChain->SetBranchAddress("nPUjets20", &nPUjets20, &b_nPUjets20);
   fChain->SetBranchAddress("njets20_5p0", &njets20_5p0, &b_njets20_5p0);
   fChain->SetBranchAddress("njetsHiggsMatch20", &njetsHiggsMatch20, &b_njetsHiggsMatch20);
   fChain->SetBranchAddress("ncompleteHiggsReco20", &ncompleteHiggsReco20, &b_ncompleteHiggsReco20);
   fChain->SetBranchAddress("njetsHiggsMatch20_eta5", &njetsHiggsMatch20_eta5, &b_njetsHiggsMatch20_eta5);
   fChain->SetBranchAddress("nbjets", &nbjets, &b_nbjets);
   fChain->SetBranchAddress("nbjets30", &nbjets30, &b_nbjets30);
   fChain->SetBranchAddress("nbjetsTweaked", &nbjetsTweaked, &b_nbjetsTweaked);
   fChain->SetBranchAddress("ntruebjets", &ntruebjets, &b_ntruebjets);
   fChain->SetBranchAddress("nGluonsSplitToLF", &nGluonsSplitToLF, &b_nGluonsSplitToLF);
   fChain->SetBranchAddress("nGluonsSplitToC", &nGluonsSplitToC, &b_nGluonsSplitToC);
   fChain->SetBranchAddress("nGluonsSplitToB", &nGluonsSplitToB, &b_nGluonsSplitToB);
   fChain->SetBranchAddress("nlfjetsFromGluons", &nlfjetsFromGluons, &b_nlfjetsFromGluons);
   fChain->SetBranchAddress("ncjetsFromGluons", &ncjetsFromGluons, &b_ncjetsFromGluons);
   fChain->SetBranchAddress("nbjetsFromGluons", &nbjetsFromGluons, &b_nbjetsFromGluons);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
   fChain->SetBranchAddress("nElectronsBug", &nElectronsBug, &b_nElectronsBug);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   fChain->SetBranchAddress("nElectrons5", &nElectrons5, &b_nElectrons5);
   fChain->SetBranchAddress("nMuons5", &nMuons5, &b_nMuons5);
   fChain->SetBranchAddress("nElectrons15", &nElectrons15, &b_nElectrons15);
   fChain->SetBranchAddress("nMuons15", &nMuons15, &b_nMuons15);
   fChain->SetBranchAddress("nTightMuons", &nTightMuons, &b_nTightMuons);
   fChain->SetBranchAddress("nElectrons20", &nElectrons20, &b_nElectrons20);
   fChain->SetBranchAddress("nMuons20", &nMuons20, &b_nMuons20);
   fChain->SetBranchAddress("nElectronsNoRelIso", &nElectronsNoRelIso, &b_nElectronsNoRelIso);
   fChain->SetBranchAddress("nMuonsNoRelIso", &nMuonsNoRelIso, &b_nMuonsNoRelIso);
   fChain->SetBranchAddress("nTausVLoose", &nTausVLoose, &b_nTausVLoose);
   fChain->SetBranchAddress("nTausLoose", &nTausLoose, &b_nTausLoose);
   fChain->SetBranchAddress("nTausMedium", &nTausMedium, &b_nTausMedium);
   fChain->SetBranchAddress("nTausTight", &nTausTight, &b_nTausTight);
   fChain->SetBranchAddress("nCorrectRecoStop", &nCorrectRecoStop, &b_nCorrectRecoStop);
   fChain->SetBranchAddress("bestZmass", &bestZmass, &b_bestZmass);
   fChain->SetBranchAddress("mjj1", &mjj1, &b_mjj1);
   fChain->SetBranchAddress("mjj2", &mjj2, &b_mjj2);
   fChain->SetBranchAddress("mjjdiff", &mjjdiff, &b_mjjdiff);
   fChain->SetBranchAddress("higgsMbb1", &higgsMbb1, &b_higgsMbb1);
   fChain->SetBranchAddress("higgsMbb2", &higgsMbb2, &b_higgsMbb2);
   fChain->SetBranchAddress("higgsMjj1", &higgsMjj1, &b_higgsMjj1);
   fChain->SetBranchAddress("higgsMjj2", &higgsMjj2, &b_higgsMjj2);
   fChain->SetBranchAddress("higgsMbb1delta", &higgsMbb1delta, &b_higgsMbb1delta);
   fChain->SetBranchAddress("higgsMbb2delta", &higgsMbb2delta, &b_higgsMbb2delta);
   fChain->SetBranchAddress("higgsMbb1MassDiff", &higgsMbb1MassDiff, &b_higgsMbb1MassDiff);
   fChain->SetBranchAddress("higgsMbb2MassDiff", &higgsMbb2MassDiff, &b_higgsMbb2MassDiff);
   fChain->SetBranchAddress("higgsMbbMassDiff_tiebreak", &higgsMbbMassDiff_tiebreak, &b_higgsMbbMassDiff_tiebreak);
   fChain->SetBranchAddress("higgsMbb1MassDiffAll", &higgsMbb1MassDiffAll, &b_higgsMbb1MassDiffAll);
   fChain->SetBranchAddress("higgsMbb2MassDiffAll", &higgsMbb2MassDiffAll, &b_higgsMbb2MassDiffAll);
   fChain->SetBranchAddress("higgsMbb1MassDiff_correct", &higgsMbb1MassDiff_correct, &b_higgsMbb1MassDiff_correct);
   fChain->SetBranchAddress("higgsMbb1MassDiffAll_correct", &higgsMbb1MassDiffAll_correct, &b_higgsMbb1MassDiffAll_correct);
   fChain->SetBranchAddress("deltaPhi_hh", &deltaPhi_hh, &b_deltaPhi_hh);
   fChain->SetBranchAddress("deltaRmax_hh", &deltaRmax_hh, &b_deltaRmax_hh);
   fChain->SetBranchAddress("deltaRmin_hh", &deltaRmin_hh, &b_deltaRmin_hh);
   fChain->SetBranchAddress("deltaEta_hh", &deltaEta_hh, &b_deltaEta_hh);
   fChain->SetBranchAddress("sumPt_hh", &sumPt_hh, &b_sumPt_hh);
   fChain->SetBranchAddress("higgsPt1", &higgsPt1, &b_higgsPt1);
   fChain->SetBranchAddress("higgsPt2", &higgsPt2, &b_higgsPt2);
   fChain->SetBranchAddress("deltaRmax_hhAll", &deltaRmax_hhAll, &b_deltaRmax_hhAll);
   fChain->SetBranchAddress("deltaRmin_hhAll", &deltaRmin_hhAll, &b_deltaRmin_hhAll);
   fChain->SetBranchAddress("higgs1jetpt1", &higgs1jetpt1, &b_higgs1jetpt1);
   fChain->SetBranchAddress("higgs1jetpt2", &higgs1jetpt2, &b_higgs1jetpt2);
   fChain->SetBranchAddress("higgs2jetpt1", &higgs2jetpt1, &b_higgs2jetpt1);
   fChain->SetBranchAddress("higgs2jetpt2", &higgs2jetpt2, &b_higgs2jetpt2);
   fChain->SetBranchAddress("higgs1jeteta1", &higgs1jeteta1, &b_higgs1jeteta1);
   fChain->SetBranchAddress("higgs1jeteta2", &higgs1jeteta2, &b_higgs1jeteta2);
   fChain->SetBranchAddress("higgs2jeteta1", &higgs2jeteta1, &b_higgs2jeteta1);
   fChain->SetBranchAddress("higgs2jeteta2", &higgs2jeteta2, &b_higgs2jeteta2);
   fChain->SetBranchAddress("higgs1CSV1", &higgs1CSV1, &b_higgs1CSV1);
   fChain->SetBranchAddress("higgs1CSV2", &higgs1CSV2, &b_higgs1CSV2);
   fChain->SetBranchAddress("higgs2CSV1", &higgs2CSV1, &b_higgs2CSV1);
   fChain->SetBranchAddress("higgs2CSV2", &higgs2CSV2, &b_higgs2CSV2);
   fChain->SetBranchAddress("higgs1partonId1", &higgs1partonId1, &b_higgs1partonId1);
   fChain->SetBranchAddress("higgs1partonId2", &higgs1partonId2, &b_higgs1partonId2);
   fChain->SetBranchAddress("higgs2partonId1", &higgs2partonId1, &b_higgs2partonId1);
   fChain->SetBranchAddress("higgs2partonId2", &higgs2partonId2, &b_higgs2partonId2);
   fChain->SetBranchAddress("higgs1partonFlavor1", &higgs1partonFlavor1, &b_higgs1partonFlavor1);
   fChain->SetBranchAddress("higgs1partonFlavor2", &higgs1partonFlavor2, &b_higgs1partonFlavor2);
   fChain->SetBranchAddress("higgs2partonFlavor1", &higgs2partonFlavor1, &b_higgs2partonFlavor1);
   fChain->SetBranchAddress("higgs2partonFlavor2", &higgs2partonFlavor2, &b_higgs2partonFlavor2);
   fChain->SetBranchAddress("higgs1partonMomId1", &higgs1partonMomId1, &b_higgs1partonMomId1);
   fChain->SetBranchAddress("higgs1partonMomId2", &higgs1partonMomId2, &b_higgs1partonMomId2);
   fChain->SetBranchAddress("higgs2partonMomId1", &higgs2partonMomId1, &b_higgs2partonMomId1);
   fChain->SetBranchAddress("higgs2partonMomId2", &higgs2partonMomId2, &b_higgs2partonMomId2);
   fChain->SetBranchAddress("higgs1tauMatch1", &higgs1tauMatch1, &b_higgs1tauMatch1);
   fChain->SetBranchAddress("higgs1tauMatch2", &higgs1tauMatch2, &b_higgs1tauMatch2);
   fChain->SetBranchAddress("higgs2tauMatch1", &higgs2tauMatch1, &b_higgs2tauMatch1);
   fChain->SetBranchAddress("higgs2tauMatch2", &higgs2tauMatch2, &b_higgs2tauMatch2);
   fChain->SetBranchAddress("higgsWCandMass", &higgsWCandMass, &b_higgsWCandMass);
   fChain->SetBranchAddress("higgsWCandMassAll", &higgsWCandMassAll, &b_higgsWCandMassAll);
   fChain->SetBranchAddress("maxDelPhiThrustJet", &maxDelPhiThrustJet, &b_maxDelPhiThrustJet);
   fChain->SetBranchAddress("minDelPhiThrustMET", &minDelPhiThrustMET, &b_minDelPhiThrustMET);
   fChain->SetBranchAddress("mjjb1", &mjjb1, &b_mjjb1);
   fChain->SetBranchAddress("mjjb2", &mjjb2, &b_mjjb2);
   fChain->SetBranchAddress("topPT1", &topPT1, &b_topPT1);
   fChain->SetBranchAddress("topPT2", &topPT2, &b_topPT2);
   fChain->SetBranchAddress("nbjetsSSVM", &nbjetsSSVM, &b_nbjetsSSVM);
   fChain->SetBranchAddress("nbjetsSSVHPT", &nbjetsSSVHPT, &b_nbjetsSSVHPT);
   fChain->SetBranchAddress("nbjetsTCHPT", &nbjetsTCHPT, &b_nbjetsTCHPT);
   fChain->SetBranchAddress("nbjetsTCHET", &nbjetsTCHET, &b_nbjetsTCHET);
   fChain->SetBranchAddress("nbjetsTCHPM", &nbjetsTCHPM, &b_nbjetsTCHPM);
   fChain->SetBranchAddress("nbjetsCSVT", &nbjetsCSVT, &b_nbjetsCSVT);
   fChain->SetBranchAddress("nbjetsCSVM", &nbjetsCSVM, &b_nbjetsCSVM);
   fChain->SetBranchAddress("nbjetsCSVL", &nbjetsCSVL, &b_nbjetsCSVL);
   fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
   fChain->SetBranchAddress("pass_utilityHLT", &pass_utilityHLT, &b_pass_utilityHLT);
   fChain->SetBranchAddress("prescaleUtilityHLT", &prescaleUtilityHLT, &b_prescaleUtilityHLT);
   fChain->SetBranchAddress("versionUtilityHLT", &versionUtilityHLT, &b_versionUtilityHLT);
   fChain->SetBranchAddress("pass_utilityPrescaleModuleHLT", &pass_utilityPrescaleModuleHLT, &b_pass_utilityPrescaleModuleHLT);
   fChain->SetBranchAddress("pass_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45", &pass_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45, &b_pass_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45);
   fChain->SetBranchAddress("pass_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45", &pass_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45, &b_pass_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45);
   fChain->SetBranchAddress("pass_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned", &pass_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned, &b_pass_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned);
   fChain->SetBranchAddress("pass_DiCentralPFJet30_PFMET80", &pass_DiCentralPFJet30_PFMET80, &b_pass_DiCentralPFJet30_PFMET80);
   fChain->SetBranchAddress("pass_DiCentralPFJet30_PFMET80_BTagCSV07", &pass_DiCentralPFJet30_PFMET80_BTagCSV07, &b_pass_DiCentralPFJet30_PFMET80_BTagCSV07);
   fChain->SetBranchAddress("pass_DiCentralPFJet30_PFMHT80", &pass_DiCentralPFJet30_PFMHT80, &b_pass_DiCentralPFJet30_PFMHT80);
   fChain->SetBranchAddress("pass_DiCentralPFJet50_PFMET80", &pass_DiCentralPFJet50_PFMET80, &b_pass_DiCentralPFJet50_PFMET80);
   fChain->SetBranchAddress("pass_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80", &pass_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80, &b_pass_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80);
   fChain->SetBranchAddress("pass_DiPFJet80_DiPFJet30_BTagCSVd07d05", &pass_DiPFJet80_DiPFJet30_BTagCSVd07d05, &b_pass_DiPFJet80_DiPFJet30_BTagCSVd07d05);
   fChain->SetBranchAddress("pass_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &pass_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_pass_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("pass_HT200", &pass_HT200, &b_pass_HT200);
   fChain->SetBranchAddress("pass_HT250", &pass_HT250, &b_pass_HT250);
   fChain->SetBranchAddress("pass_HT250_AlphaT0p55", &pass_HT250_AlphaT0p55, &b_pass_HT250_AlphaT0p55);
   fChain->SetBranchAddress("pass_HT300", &pass_HT300, &b_pass_HT300);
   fChain->SetBranchAddress("pass_HT300_AlphaT0p53", &pass_HT300_AlphaT0p53, &b_pass_HT300_AlphaT0p53);
   fChain->SetBranchAddress("pass_IsoMu24", &pass_IsoMu24, &b_pass_IsoMu24);
   fChain->SetBranchAddress("pass_IsoMu24_eta2p1", &pass_IsoMu24_eta2p1, &b_pass_IsoMu24_eta2p1);
   fChain->SetBranchAddress("pass_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV", &pass_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV, &b_pass_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV);
   fChain->SetBranchAddress("pass_L1ETM40", &pass_L1ETM40, &b_pass_L1ETM40);
   fChain->SetBranchAddress("pass_MET120_HBHENoiseCleaned", &pass_MET120_HBHENoiseCleaned, &b_pass_MET120_HBHENoiseCleaned);
   fChain->SetBranchAddress("pass_MET200", &pass_MET200, &b_pass_MET200);
   fChain->SetBranchAddress("pass_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &pass_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_pass_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("pass_Mu17_Mu8", &pass_Mu17_Mu8, &b_pass_Mu17_Mu8);
   fChain->SetBranchAddress("pass_Mu8_DiJet30", &pass_Mu8_DiJet30, &b_pass_Mu8_DiJet30);
   fChain->SetBranchAddress("pass_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &pass_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_pass_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("pass_PFHT350", &pass_PFHT350, &b_pass_PFHT350);
   fChain->SetBranchAddress("pass_PFHT350_Mu15_PFMET45", &pass_PFHT350_Mu15_PFMET45, &b_pass_PFHT350_Mu15_PFMET45);
   fChain->SetBranchAddress("pass_PFHT350_PFMET100", &pass_PFHT350_PFMET100, &b_pass_PFHT350_PFMET100);
   fChain->SetBranchAddress("pass_PFHT650", &pass_PFHT650, &b_pass_PFHT650);
   fChain->SetBranchAddress("pass_PFJet80", &pass_PFJet80, &b_pass_PFJet80);
   fChain->SetBranchAddress("pass_PFMET150", &pass_PFMET150, &b_pass_PFMET150);
   fChain->SetBranchAddress("pass_PFNoPUHT350", &pass_PFNoPUHT350, &b_pass_PFNoPUHT350);
   fChain->SetBranchAddress("pass_PFNoPUHT350_Mu15_PFMET45", &pass_PFNoPUHT350_Mu15_PFMET45, &b_pass_PFNoPUHT350_Mu15_PFMET45);
   fChain->SetBranchAddress("pass_PFNoPUHT350_PFMET100", &pass_PFNoPUHT350_PFMET100, &b_pass_PFNoPUHT350_PFMET100);
   fChain->SetBranchAddress("pass_PFNoPUHT650", &pass_PFNoPUHT650, &b_pass_PFNoPUHT650);
   fChain->SetBranchAddress("pass_Photon135", &pass_Photon135, &b_pass_Photon135);
   fChain->SetBranchAddress("pass_Photon150", &pass_Photon150, &b_pass_Photon150);
   fChain->SetBranchAddress("pass_QuadJet80", &pass_QuadJet80, &b_pass_QuadJet80);
   fChain->SetBranchAddress("pass_SixJet45", &pass_SixJet45, &b_pass_SixJet45);
   fChain->SetBranchAddress("passMC_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45", &passMC_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45, &b_passMC_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45);
   fChain->SetBranchAddress("passMC_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45", &passMC_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45, &b_passMC_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45);
   fChain->SetBranchAddress("passMC_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned", &passMC_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned, &b_passMC_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned);
   fChain->SetBranchAddress("passMC_DiCentralPFJet30_PFMET80", &passMC_DiCentralPFJet30_PFMET80, &b_passMC_DiCentralPFJet30_PFMET80);
   fChain->SetBranchAddress("passMC_DiCentralPFJet30_PFMET80_BTagCSV07", &passMC_DiCentralPFJet30_PFMET80_BTagCSV07, &b_passMC_DiCentralPFJet30_PFMET80_BTagCSV07);
   fChain->SetBranchAddress("passMC_DiCentralPFJet30_PFMHT80", &passMC_DiCentralPFJet30_PFMHT80, &b_passMC_DiCentralPFJet30_PFMHT80);
   fChain->SetBranchAddress("passMC_DiCentralPFJet50_PFMET80", &passMC_DiCentralPFJet50_PFMET80, &b_passMC_DiCentralPFJet50_PFMET80);
   fChain->SetBranchAddress("passMC_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80", &passMC_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80, &b_passMC_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80);
   fChain->SetBranchAddress("passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05", &passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05, &b_passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05);
   fChain->SetBranchAddress("passMC_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &passMC_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_passMC_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("passMC_HT200", &passMC_HT200, &b_passMC_HT200);
   fChain->SetBranchAddress("passMC_HT250", &passMC_HT250, &b_passMC_HT250);
   fChain->SetBranchAddress("passMC_HT250_AlphaT0p55", &passMC_HT250_AlphaT0p55, &b_passMC_HT250_AlphaT0p55);
   fChain->SetBranchAddress("passMC_HT300", &passMC_HT300, &b_passMC_HT300);
   fChain->SetBranchAddress("passMC_HT300_AlphaT0p53", &passMC_HT300_AlphaT0p53, &b_passMC_HT300_AlphaT0p53);
   fChain->SetBranchAddress("passMC_IsoMu24", &passMC_IsoMu24, &b_passMC_IsoMu24);
   fChain->SetBranchAddress("passMC_IsoMu24_eta2p1", &passMC_IsoMu24_eta2p1, &b_passMC_IsoMu24_eta2p1);
   fChain->SetBranchAddress("passMC_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV", &passMC_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV, &b_passMC_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV);
   fChain->SetBranchAddress("passMC_L1ETM40", &passMC_L1ETM40, &b_passMC_L1ETM40);
   fChain->SetBranchAddress("passMC_MET120_HBHENoiseCleaned", &passMC_MET120_HBHENoiseCleaned, &b_passMC_MET120_HBHENoiseCleaned);
   fChain->SetBranchAddress("passMC_MET200", &passMC_MET200, &b_passMC_MET200);
   fChain->SetBranchAddress("passMC_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &passMC_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_passMC_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("passMC_Mu17_Mu8", &passMC_Mu17_Mu8, &b_passMC_Mu17_Mu8);
   fChain->SetBranchAddress("passMC_Mu8_DiJet30", &passMC_Mu8_DiJet30, &b_passMC_Mu8_DiJet30);
   fChain->SetBranchAddress("passMC_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &passMC_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_passMC_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("passMC_PFHT350", &passMC_PFHT350, &b_passMC_PFHT350);
   fChain->SetBranchAddress("passMC_PFHT350_Mu15_PFMET45", &passMC_PFHT350_Mu15_PFMET45, &b_passMC_PFHT350_Mu15_PFMET45);
   fChain->SetBranchAddress("passMC_PFHT350_PFMET100", &passMC_PFHT350_PFMET100, &b_passMC_PFHT350_PFMET100);
   fChain->SetBranchAddress("passMC_PFHT650", &passMC_PFHT650, &b_passMC_PFHT650);
   fChain->SetBranchAddress("passMC_PFJet80", &passMC_PFJet80, &b_passMC_PFJet80);
   fChain->SetBranchAddress("passMC_PFMET150", &passMC_PFMET150, &b_passMC_PFMET150);
   fChain->SetBranchAddress("passMC_PFNoPUHT350", &passMC_PFNoPUHT350, &b_passMC_PFNoPUHT350);
   fChain->SetBranchAddress("passMC_PFNoPUHT350_Mu15_PFMET45", &passMC_PFNoPUHT350_Mu15_PFMET45, &b_passMC_PFNoPUHT350_Mu15_PFMET45);
   fChain->SetBranchAddress("passMC_PFNoPUHT350_PFMET100", &passMC_PFNoPUHT350_PFMET100, &b_passMC_PFNoPUHT350_PFMET100);
   fChain->SetBranchAddress("passMC_PFNoPUHT650", &passMC_PFNoPUHT650, &b_passMC_PFNoPUHT650);
   fChain->SetBranchAddress("passMC_Photon135", &passMC_Photon135, &b_passMC_Photon135);
   fChain->SetBranchAddress("passMC_Photon150", &passMC_Photon150, &b_passMC_Photon150);
   fChain->SetBranchAddress("passMC_QuadJet80", &passMC_QuadJet80, &b_passMC_QuadJet80);
   fChain->SetBranchAddress("passMC_SixJet45", &passMC_SixJet45, &b_passMC_SixJet45);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("HT30", &HT30, &b_HT30);
   fChain->SetBranchAddress("HT20", &HT20, &b_HT20);
   fChain->SetBranchAddress("ST", &ST, &b_ST);
   fChain->SetBranchAddress("STeff", &STeff, &b_STeff);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("METphi", &METphi, &b_METphi);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("METsig", &METsig, &b_METsig);
   fChain->SetBranchAddress("METsig00", &METsig00, &b_METsig00);
   fChain->SetBranchAddress("METsig10", &METsig10, &b_METsig10);
   fChain->SetBranchAddress("METsig11", &METsig11, &b_METsig11);
   fChain->SetBranchAddress("METsig_2012", &METsig_2012, &b_METsig_2012);
   fChain->SetBranchAddress("METsig00_2012", &METsig00_2012, &b_METsig00_2012);
   fChain->SetBranchAddress("METsig10_2012", &METsig10_2012, &b_METsig10_2012);
   fChain->SetBranchAddress("METsig11_2012", &METsig11_2012, &b_METsig11_2012);
   fChain->SetBranchAddress("caloMET", &caloMET, &b_caloMET);
   fChain->SetBranchAddress("caloMETphi", &caloMETphi, &b_caloMETphi);
   fChain->SetBranchAddress("rawPFMET", &rawPFMET, &b_rawPFMET);
   fChain->SetBranchAddress("rawPFMETphi", &rawPFMETphi, &b_rawPFMETphi);
   fChain->SetBranchAddress("caloMETwithHO", &caloMETwithHO, &b_caloMETwithHO);
   fChain->SetBranchAddress("caloMETwithHOphi", &caloMETwithHOphi, &b_caloMETwithHOphi);
   fChain->SetBranchAddress("unclusteredMET", &unclusteredMET, &b_unclusteredMET);
   fChain->SetBranchAddress("unclusteredMETphi", &unclusteredMETphi, &b_unclusteredMETphi);
   fChain->SetBranchAddress("bestWMass", &bestWMass, &b_bestWMass);
   fChain->SetBranchAddress("bestTopMass", &bestTopMass, &b_bestTopMass);
   fChain->SetBranchAddress("topCosHel", &topCosHel, &b_topCosHel);
   fChain->SetBranchAddress("WCosHel", &WCosHel, &b_WCosHel);
   fChain->SetBranchAddress("MT_b", &MT_b, &b_MT_b);
   fChain->SetBranchAddress("MT_Hbb1", &MT_Hbb1, &b_MT_Hbb1);
   fChain->SetBranchAddress("MT_Hbb2", &MT_Hbb2, &b_MT_Hbb2);
   fChain->SetBranchAddress("MT_Wlep", &MT_Wlep, &b_MT_Wlep);
   fChain->SetBranchAddress("MT_Wlep5", &MT_Wlep5, &b_MT_Wlep5);
   fChain->SetBranchAddress("MT_Wlep15", &MT_Wlep15, &b_MT_Wlep15);
   fChain->SetBranchAddress("maxDeltaPhi_bbbb", &maxDeltaPhi_bbbb, &b_maxDeltaPhi_bbbb);
   fChain->SetBranchAddress("maxDeltaPhi_bb_bb", &maxDeltaPhi_bb_bb, &b_maxDeltaPhi_bb_bb);
   fChain->SetBranchAddress("MT_bestCSV_gencode", &MT_bestCSV_gencode, &b_MT_bestCSV_gencode);
   fChain->SetBranchAddress("MT_bestCSV", &MT_bestCSV, &b_MT_bestCSV);
   fChain->SetBranchAddress("MT_jim", &MT_jim, &b_MT_jim);
   fChain->SetBranchAddress("minMT_jetMET", &minMT_jetMET, &b_minMT_jetMET);
   fChain->SetBranchAddress("minMT_bMET", &minMT_bMET, &b_minMT_bMET);
   fChain->SetBranchAddress("maxMT_bMET", &maxMT_bMET, &b_maxMT_bMET);
   fChain->SetBranchAddress("deltaThetaT", &deltaThetaT, &b_deltaThetaT);
   fChain->SetBranchAddress("deltaR_bestTwoCSV", &deltaR_bestTwoCSV, &b_deltaR_bestTwoCSV);
   fChain->SetBranchAddress("deltaPhi_bestTwoCSV", &deltaPhi_bestTwoCSV, &b_deltaPhi_bestTwoCSV);
   fChain->SetBranchAddress("mjj_bestTwoCSV", &mjj_bestTwoCSV, &b_mjj_bestTwoCSV);
   fChain->SetBranchAddress("deltaR_closestB", &deltaR_closestB, &b_deltaR_closestB);
   fChain->SetBranchAddress("mjj_closestB", &mjj_closestB, &b_mjj_closestB);
   fChain->SetBranchAddress("mjj_closestB20", &mjj_closestB20, &b_mjj_closestB20);
   fChain->SetBranchAddress("mjj_closestB20_correct", &mjj_closestB20_correct, &b_mjj_closestB20_correct);
   fChain->SetBranchAddress("deltaR_closestB20", &deltaR_closestB20, &b_deltaR_closestB20);
   fChain->SetBranchAddress("cosHel_closestB20", &cosHel_closestB20, &b_cosHel_closestB20);
   fChain->SetBranchAddress("deltaPhi_hb20", &deltaPhi_hb20, &b_deltaPhi_hb20);
   fChain->SetBranchAddress("sumPtMjjDiff_closestB20", &sumPtMjjDiff_closestB20, &b_sumPtMjjDiff_closestB20);
   fChain->SetBranchAddress("mjj_h125", &mjj_h125, &b_mjj_h125);
   fChain->SetBranchAddress("mjj_highptCSVM", &mjj_highptCSVM, &b_mjj_highptCSVM);
   fChain->SetBranchAddress("mjj_highptCSVT", &mjj_highptCSVT, &b_mjj_highptCSVT);
   fChain->SetBranchAddress("minDeltaPhi", &minDeltaPhi, &b_minDeltaPhi);
   fChain->SetBranchAddress("minDeltaPhiAll", &minDeltaPhiAll, &b_minDeltaPhiAll);
   fChain->SetBranchAddress("minDeltaPhiAll20", &minDeltaPhiAll20, &b_minDeltaPhiAll20);
   fChain->SetBranchAddress("minDeltaPhi20", &minDeltaPhi20, &b_minDeltaPhi20);
   fChain->SetBranchAddress("minDeltaPhi30", &minDeltaPhi30, &b_minDeltaPhi30);
   fChain->SetBranchAddress("minDeltaPhi20_eta5_noIdAll", &minDeltaPhi20_eta5_noIdAll, &b_minDeltaPhi20_eta5_noIdAll);
   fChain->SetBranchAddress("minDeltaPhi20_eta5_noIdAll_nobeta", &minDeltaPhi20_eta5_noIdAll_nobeta, &b_minDeltaPhi20_eta5_noIdAll_nobeta);
   fChain->SetBranchAddress("deltaPhi1", &deltaPhi1, &b_deltaPhi1);
   fChain->SetBranchAddress("deltaPhi2", &deltaPhi2, &b_deltaPhi2);
   fChain->SetBranchAddress("deltaPhi3", &deltaPhi3, &b_deltaPhi3);
   fChain->SetBranchAddress("deltaPhi1_20", &deltaPhi1_20, &b_deltaPhi1_20);
   fChain->SetBranchAddress("deltaPhi2_20", &deltaPhi2_20, &b_deltaPhi2_20);
   fChain->SetBranchAddress("deltaPhi3_20", &deltaPhi3_20, &b_deltaPhi3_20);
   fChain->SetBranchAddress("maxDeltaPhi", &maxDeltaPhi, &b_maxDeltaPhi);
   fChain->SetBranchAddress("maxDeltaPhiAll", &maxDeltaPhiAll, &b_maxDeltaPhiAll);
   fChain->SetBranchAddress("maxDeltaPhiAll30", &maxDeltaPhiAll30, &b_maxDeltaPhiAll30);
   fChain->SetBranchAddress("maxDeltaPhi20_eta5_noIdAll", &maxDeltaPhi20_eta5_noIdAll, &b_maxDeltaPhi20_eta5_noIdAll);
   fChain->SetBranchAddress("sumDeltaPhi", &sumDeltaPhi, &b_sumDeltaPhi);
   fChain->SetBranchAddress("diffDeltaPhi", &diffDeltaPhi, &b_diffDeltaPhi);
   fChain->SetBranchAddress("deltaPhiStar", &deltaPhiStar, &b_deltaPhiStar);
   fChain->SetBranchAddress("deltaPhiStar_badjet_pt", &deltaPhiStar_badjet_pt, &b_deltaPhiStar_badjet_pt);
   fChain->SetBranchAddress("deltaPhiStar_badjet_estimatedPt", &deltaPhiStar_badjet_estimatedPt, &b_deltaPhiStar_badjet_estimatedPt);
   fChain->SetBranchAddress("deltaPhiStar_badjet_genPt", &deltaPhiStar_badjet_genPt, &b_deltaPhiStar_badjet_genPt);
   fChain->SetBranchAddress("deltaPhiStar_badjet_energy", &deltaPhiStar_badjet_energy, &b_deltaPhiStar_badjet_energy);
   fChain->SetBranchAddress("deltaPhiStar_badjet_phi", &deltaPhiStar_badjet_phi, &b_deltaPhiStar_badjet_phi);
   fChain->SetBranchAddress("deltaPhiStar_badjet_eta", &deltaPhiStar_badjet_eta, &b_deltaPhiStar_badjet_eta);
   fChain->SetBranchAddress("deltaPhiStar_badjet_CSV", &deltaPhiStar_badjet_CSV, &b_deltaPhiStar_badjet_CSV);
   fChain->SetBranchAddress("deltaPhiStar_badjet_chEmE", &deltaPhiStar_badjet_chEmE, &b_deltaPhiStar_badjet_chEmE);
   fChain->SetBranchAddress("deltaPhiStar_badjet_chHadE", &deltaPhiStar_badjet_chHadE, &b_deltaPhiStar_badjet_chHadE);
   fChain->SetBranchAddress("deltaPhiStar_badjet_photonE", &deltaPhiStar_badjet_photonE, &b_deltaPhiStar_badjet_photonE);
   fChain->SetBranchAddress("deltaPhiStar_badjet_neuEmE", &deltaPhiStar_badjet_neuEmE, &b_deltaPhiStar_badjet_neuEmE);
   fChain->SetBranchAddress("deltaPhiStar_badjet_neuHadE", &deltaPhiStar_badjet_neuHadE, &b_deltaPhiStar_badjet_neuHadE);
   fChain->SetBranchAddress("deltaPhiStar_badjet_chMuE", &deltaPhiStar_badjet_chMuE, &b_deltaPhiStar_badjet_chMuE);
   fChain->SetBranchAddress("deltaPhiStar_badjet_etaetaMoment", &deltaPhiStar_badjet_etaetaMoment, &b_deltaPhiStar_badjet_etaetaMoment);
   fChain->SetBranchAddress("deltaPhiStar_badjet_etaphiMoment", &deltaPhiStar_badjet_etaphiMoment, &b_deltaPhiStar_badjet_etaphiMoment);
   fChain->SetBranchAddress("deltaPhiStar_badjet_phiphiMoment", &deltaPhiStar_badjet_phiphiMoment, &b_deltaPhiStar_badjet_phiphiMoment);
   fChain->SetBranchAddress("deltaPhiStar_badjet_rawPt", &deltaPhiStar_badjet_rawPt, &b_deltaPhiStar_badjet_rawPt);
   fChain->SetBranchAddress("deltaPhiStar_badjet_mass", &deltaPhiStar_badjet_mass, &b_deltaPhiStar_badjet_mass);
   fChain->SetBranchAddress("deltaPhiStar_badjet_PUbeta", &deltaPhiStar_badjet_PUbeta, &b_deltaPhiStar_badjet_PUbeta);
   fChain->SetBranchAddress("deltaPhiStar_badjet_PUbetaStar", &deltaPhiStar_badjet_PUbetaStar, &b_deltaPhiStar_badjet_PUbetaStar);
   fChain->SetBranchAddress("deltaPhiStar_badjet_PUbetaStarClassic", &deltaPhiStar_badjet_PUbetaStarClassic, &b_deltaPhiStar_badjet_PUbetaStarClassic);
   fChain->SetBranchAddress("deltaPhiStar_badjet_chMult", &deltaPhiStar_badjet_chMult, &b_deltaPhiStar_badjet_chMult);
   fChain->SetBranchAddress("deltaPhiStar_badjet_neuMult", &deltaPhiStar_badjet_neuMult, &b_deltaPhiStar_badjet_neuMult);
   fChain->SetBranchAddress("deltaPhiStar_badjet_muMult", &deltaPhiStar_badjet_muMult, &b_deltaPhiStar_badjet_muMult);
   fChain->SetBranchAddress("deltaPhiStar_badjet_n60", &deltaPhiStar_badjet_n60, &b_deltaPhiStar_badjet_n60);
   fChain->SetBranchAddress("deltaPhiStar_badjet_n90", &deltaPhiStar_badjet_n90, &b_deltaPhiStar_badjet_n90);
   fChain->SetBranchAddress("hh_mostPUlike_beta", &hh_mostPUlike_beta, &b_hh_mostPUlike_beta);
   fChain->SetBranchAddress("hh_mostPUlike_betaStar", &hh_mostPUlike_betaStar, &b_hh_mostPUlike_betaStar);
   fChain->SetBranchAddress("minDeltaPhiN", &minDeltaPhiN, &b_minDeltaPhiN);
   fChain->SetBranchAddress("deltaPhiN1", &deltaPhiN1, &b_deltaPhiN1);
   fChain->SetBranchAddress("deltaPhiN2", &deltaPhiN2, &b_deltaPhiN2);
   fChain->SetBranchAddress("deltaPhiN3", &deltaPhiN3, &b_deltaPhiN3);
   fChain->SetBranchAddress("minDeltaPhi_chosenJet", &minDeltaPhi_chosenJet, &b_minDeltaPhi_chosenJet);
   fChain->SetBranchAddress("minDeltaPhiN_chosenJet", &minDeltaPhiN_chosenJet, &b_minDeltaPhiN_chosenJet);
   fChain->SetBranchAddress("maxJetMis_chosenJet", &maxJetMis_chosenJet, &b_maxJetMis_chosenJet);
   fChain->SetBranchAddress("maxJetFracMis_chosenJet", &maxJetFracMis_chosenJet, &b_maxJetFracMis_chosenJet);
   fChain->SetBranchAddress("minDeltaPhiN_deltaT", &minDeltaPhiN_deltaT, &b_minDeltaPhiN_deltaT);
   fChain->SetBranchAddress("deltaT1", &deltaT1, &b_deltaT1);
   fChain->SetBranchAddress("deltaT2", &deltaT2, &b_deltaT2);
   fChain->SetBranchAddress("deltaT3", &deltaT3, &b_deltaT3);
   fChain->SetBranchAddress("minDeltaPhiN_asin", &minDeltaPhiN_asin, &b_minDeltaPhiN_asin);
   fChain->SetBranchAddress("deltaPhiN_asin1", &deltaPhiN_asin1, &b_deltaPhiN_asin1);
   fChain->SetBranchAddress("deltaPhiN_asin2", &deltaPhiN_asin2, &b_deltaPhiN_asin2);
   fChain->SetBranchAddress("deltaPhiN_asin3", &deltaPhiN_asin3, &b_deltaPhiN_asin3);
   fChain->SetBranchAddress("CSVout1", &CSVout1, &b_CSVout1);
   fChain->SetBranchAddress("CSVout2", &CSVout2, &b_CSVout2);
   fChain->SetBranchAddress("CSVout3", &CSVout3, &b_CSVout3);
   fChain->SetBranchAddress("CSVbest1", &CSVbest1, &b_CSVbest1);
   fChain->SetBranchAddress("CSVbest2", &CSVbest2, &b_CSVbest2);
   fChain->SetBranchAddress("CSVbest3", &CSVbest3, &b_CSVbest3);
   fChain->SetBranchAddress("CSVbest4", &CSVbest4, &b_CSVbest4);
   fChain->SetBranchAddress("minDeltaPhiAllb30", &minDeltaPhiAllb30, &b_minDeltaPhiAllb30);
   fChain->SetBranchAddress("deltaPhib1", &deltaPhib1, &b_deltaPhib1);
   fChain->SetBranchAddress("deltaPhib2", &deltaPhib2, &b_deltaPhib2);
   fChain->SetBranchAddress("deltaPhib3", &deltaPhib3, &b_deltaPhib3);
   fChain->SetBranchAddress("minDeltaPhiMETMuonsAll", &minDeltaPhiMETMuonsAll, &b_minDeltaPhiMETMuonsAll);
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
   fChain->SetBranchAddress("jetneutralhadronmult1", &jetneutralhadronmult1, &b_jetneutralhadronmult1);
   fChain->SetBranchAddress("jetmumult1", &jetmumult1, &b_jetmumult1);
   fChain->SetBranchAddress("alletajetpt1", &alletajetpt1, &b_alletajetpt1);
   fChain->SetBranchAddress("alletajetphi1", &alletajetphi1, &b_alletajetphi1);
   fChain->SetBranchAddress("alletajeteta1", &alletajeteta1, &b_alletajeteta1);
   fChain->SetBranchAddress("alletajetneutralhadronfrac1", &alletajetneutralhadronfrac1, &b_alletajetneutralhadronfrac1);
   fChain->SetBranchAddress("alletajetneutralemfrac1", &alletajetneutralemfrac1, &b_alletajetneutralemfrac1);
   fChain->SetBranchAddress("alletajetphotonfrac1", &alletajetphotonfrac1, &b_alletajetphotonfrac1);
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
   fChain->SetBranchAddress("jetneutralhadronmult2", &jetneutralhadronmult2, &b_jetneutralhadronmult2);
   fChain->SetBranchAddress("jetmumult2", &jetmumult2, &b_jetmumult2);
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
   fChain->SetBranchAddress("jetneutralhadronmult3", &jetneutralhadronmult3, &b_jetneutralhadronmult3);
   fChain->SetBranchAddress("jetmumult3", &jetmumult3, &b_jetmumult3);
   fChain->SetBranchAddress("jetpt4", &jetpt4, &b_jetpt4);
   fChain->SetBranchAddress("jetgenpt4", &jetgenpt4, &b_jetgenpt4);
   fChain->SetBranchAddress("jeteta4", &jeteta4, &b_jeteta4);
   fChain->SetBranchAddress("jetgeneta4", &jetgeneta4, &b_jetgeneta4);
   fChain->SetBranchAddress("jetphi4", &jetphi4, &b_jetphi4);
   fChain->SetBranchAddress("jetgenphi4", &jetgenphi4, &b_jetgenphi4);
   fChain->SetBranchAddress("jetenergy4", &jetenergy4, &b_jetenergy4);
   fChain->SetBranchAddress("jetflavor4", &jetflavor4, &b_jetflavor4);
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
   fChain->SetBranchAddress("bjetpt4", &bjetpt4, &b_bjetpt4);
   fChain->SetBranchAddress("bjeteta4", &bjeteta4, &b_bjeteta4);
   fChain->SetBranchAddress("bjetphi4", &bjetphi4, &b_bjetphi4);
   fChain->SetBranchAddress("bjetenergy4", &bjetenergy4, &b_bjetenergy4);
   fChain->SetBranchAddress("bjetflavor4", &bjetflavor4, &b_bjetflavor4);
   fChain->SetBranchAddress("tauMatch_jetpt", &tauMatch_jetpt, &b_tauMatch_jetpt);
   fChain->SetBranchAddress("tauMatch_chhadmult", &tauMatch_chhadmult, &b_tauMatch_chhadmult);
   fChain->SetBranchAddress("tauMatch_jetcsv", &tauMatch_jetcsv, &b_tauMatch_jetcsv);
   fChain->SetBranchAddress("tauMatch_MT", &tauMatch_MT, &b_tauMatch_MT);
   fChain->SetBranchAddress("eleet1", &eleet1, &b_eleet1);
   fChain->SetBranchAddress("elephi1", &elephi1, &b_elephi1);
   fChain->SetBranchAddress("eleeta1", &eleeta1, &b_eleeta1);
   fChain->SetBranchAddress("elecharge1", &elecharge1, &b_elecharge1);
   fChain->SetBranchAddress("muonpt1", &muonpt1, &b_muonpt1);
   fChain->SetBranchAddress("muonphi1", &muonphi1, &b_muonphi1);
   fChain->SetBranchAddress("muoneta1", &muoneta1, &b_muoneta1);
   fChain->SetBranchAddress("muoncharge1", &muoncharge1, &b_muoncharge1);
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
   fChain->SetBranchAddress("isotrackpt1", &isotrackpt1, &b_isotrackpt1);
   fChain->SetBranchAddress("isotracketa1", &isotracketa1, &b_isotracketa1);
   fChain->SetBranchAddress("isotrackphi1", &isotrackphi1, &b_isotrackphi1);
   fChain->SetBranchAddress("eleRelIso", &eleRelIso, &b_eleRelIso);
   fChain->SetBranchAddress("muonRelIso", &muonRelIso, &b_muonRelIso);
   fChain->SetBranchAddress("rl", &rl, &b_rl);
   fChain->SetBranchAddress("rMET", &rMET, &b_rMET);
   fChain->SetBranchAddress("transverseSphericity_jets", &transverseSphericity_jets, &b_transverseSphericity_jets);
   fChain->SetBranchAddress("transverseSphericity_jetsMet", &transverseSphericity_jetsMet, &b_transverseSphericity_jetsMet);
   fChain->SetBranchAddress("transverseSphericity_jetsMetLeptons", &transverseSphericity_jetsMetLeptons, &b_transverseSphericity_jetsMetLeptons);
   fChain->SetBranchAddress("transverseSphericity_jetsLeptons", &transverseSphericity_jetsLeptons, &b_transverseSphericity_jetsLeptons);
   fChain->SetBranchAddress("transverseSphericity_jets30", &transverseSphericity_jets30, &b_transverseSphericity_jets30);
   fChain->SetBranchAddress("transverseSphericity_jets30Met", &transverseSphericity_jets30Met, &b_transverseSphericity_jets30Met);
   fChain->SetBranchAddress("transverseSphericity_jets30MetLeptons", &transverseSphericity_jets30MetLeptons, &b_transverseSphericity_jets30MetLeptons);
   fChain->SetBranchAddress("transverseSphericity_jets30Leptons", &transverseSphericity_jets30Leptons, &b_transverseSphericity_jets30Leptons);
   fChain->SetBranchAddress("minTransverseMETSignificance", &minTransverseMETSignificance, &b_minTransverseMETSignificance);
   fChain->SetBranchAddress("maxTransverseMETSignificance", &maxTransverseMETSignificance, &b_maxTransverseMETSignificance);
   fChain->SetBranchAddress("transverseMETSignificance1", &transverseMETSignificance1, &b_transverseMETSignificance1);
   fChain->SetBranchAddress("transverseMETSignificance2", &transverseMETSignificance2, &b_transverseMETSignificance2);
   fChain->SetBranchAddress("transverseMETSignificance3", &transverseMETSignificance3, &b_transverseMETSignificance3);
   fChain->SetBranchAddress("nIsoTracks20_005_03", &nIsoTracks20_005_03, &b_nIsoTracks20_005_03);
   fChain->SetBranchAddress("nIsoTracks15_005_03", &nIsoTracks15_005_03, &b_nIsoTracks15_005_03);
   fChain->SetBranchAddress("nIsoTracks10_005_03", &nIsoTracks10_005_03, &b_nIsoTracks10_005_03);
   fChain->SetBranchAddress("nIsoTracks5_005_03", &nIsoTracks5_005_03, &b_nIsoTracks5_005_03);
   fChain->SetBranchAddress("nIsoTracks15_005_03_lepcleaned", &nIsoTracks15_005_03_lepcleaned, &b_nIsoTracks15_005_03_lepcleaned);
   fChain->SetBranchAddress("nIsoPFcands5_010", &nIsoPFcands5_010, &b_nIsoPFcands5_010);
   fChain->SetBranchAddress("nIsoPFcands10_010", &nIsoPFcands10_010, &b_nIsoPFcands10_010);
   fChain->SetBranchAddress("nIsoPFcands15_010", &nIsoPFcands15_010, &b_nIsoPFcands15_010);
   fChain->SetBranchAddress("minTrackIso10", &minTrackIso10, &b_minTrackIso10);
   fChain->SetBranchAddress("minTrackIso5", &minTrackIso5, &b_minTrackIso5);
   fChain->SetBranchAddress("minTrackIso15", &minTrackIso15, &b_minTrackIso15);
   fChain->SetBranchAddress("minTrackIso20", &minTrackIso20, &b_minTrackIso20);
   fChain->SetBranchAddress("trackd0_5", &trackd0_5, &b_trackd0_5);
   fChain->SetBranchAddress("trackd0_10", &trackd0_10, &b_trackd0_10);
   fChain->SetBranchAddress("trackd0_15", &trackd0_15, &b_trackd0_15);
   fChain->SetBranchAddress("trackd0_20", &trackd0_20, &b_trackd0_20);
   fChain->SetBranchAddress("trackpt_5", &trackpt_5, &b_trackpt_5);
   fChain->SetBranchAddress("trackpt_10", &trackpt_10, &b_trackpt_10);
   fChain->SetBranchAddress("trackpt_15", &trackpt_15, &b_trackpt_15);
   fChain->SetBranchAddress("trackpt_20", &trackpt_20, &b_trackpt_20);

   Notify();
}

Bool_t signalEff_hbbhbb::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void signalEff_hbbhbb::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t signalEff_hbbhbb::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef signalEff_hbbhbb_cxx
