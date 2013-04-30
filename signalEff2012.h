//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct  2 16:20:01 2012 by ROOT version 5.32/00
// from TTree reducedTree/tree with minimal cuts
// found on file: reducedTree.CSVM_PF2PATjets_JES0_JERbias_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.SMS-T1bbbb-v65.root
//////////////////////////////////////////////////////////

#ifndef signalEff2012_h
#define signalEff2012_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include "TH3.h"

#include <iostream>

enum IsrMode { kNoIsrWeight, kIsr0, kIsrUp,kIsrDown };

// Header file for the classes stored in the TTree if any.
class SearchRegion { //a DIFFERENT class than my other SearchRegion class...oh well
public:
  SearchRegion(int minb, int maxb,float minHT, float maxHT,float minMET,float maxMET,bool isLDP,bool isSL,bool isSLsig,int minnjets,TString pdfset,int pdfindex);
  ~SearchRegion();
  //  void Print() const;

  int minb_;
  int maxb_;
  float minHT_;
  float maxHT_;
  float minMET_;
  float maxMET_;
  bool isLDP_;
  bool isSL_;
  bool isSLsig_;
  int minnjets_;

  TString pdfset_;
  int pdfindex_;

  TString id() const;

  //  bool operator==(const SearchRegion& other);
  //  bool operator!=(const SearchRegion& other) {return !((*this)==other);}

};


SearchRegion::SearchRegion(int minb, int maxb,float minHT, float maxHT,float minMET,float maxMET,bool isLDP,bool isSL,bool isSLsig,int minnjets,
			   TString pdfset,int pdfindex) :
  minb_(minb), maxb_(maxb),minHT_(minHT),maxHT_(maxHT),minMET_(minMET),maxMET_(maxMET),isLDP_(isLDP),isSL_(isSL),isSLsig_(isSLsig),minnjets_(minnjets)
  ,pdfset_(pdfset),pdfindex_(pdfindex)
{

  assert( !(isSL_ && isSLsig));

}

SearchRegion::~SearchRegion() {}

//void SearchRegion::Print() const {
  //  std::cout<<this->id()<<std::endl;
  //  std::cout<<"hello world"<<std::endl;

//}

//bool SearchRegion::operator==(const SearchRegion& other) {
//  //ht, met, btag cuts must match
//  if ( htSelection != other.htSelection ) return false;
//  if ( metSelection != other.metSelection ) return false;
//  if ( btagSelection != other.btagSelection ) return false;
//  //also must both be either SIG or SB
//  if (isSIG != other.isSIG) return false;

//  //don't care about owenId
//  return true;
//}

TString SearchRegion::id() const {

  TString ldpstring = isLDP_ ? "_LDP" : "";
  TString slstring = isSL_ ? "_SL" : "";
  TString slsigstring = isSLsig_ ? "_SLsig" : "";
  TString njetsstring = (minnjets_==3) ? "" : "_minnjets"+minnjets_;

  TString theid;
  if ((minb_==maxb_) || (maxb_>10))
    theid.Form("b%d_HT%.0fto%.0f_MET%.0fto%.0f%s%s%s%s",minb_,minHT_,maxHT_,minMET_,maxMET_,njetsstring.Data(),ldpstring.Data(),slstring.Data(),slsigstring.Data());
  else
    theid.Form("b%dto%d_HT%.0fto%.0f_MET%.0fto%.0f%s%s%s%s",minb_,maxb_,minHT_,maxHT_,minMET_,maxMET_,njetsstring.Data(),ldpstring.Data(),slstring.Data(),slsigstring.Data());

  if (pdfset_ != "none") {
    TString suffix;
    suffix.Form("_%s%d",pdfset_.Data(),pdfindex_);
    theid+=suffix;
  }

  return theid;
}


// Fixed size dimensions of array or collections stored in the TTree if any.

class signalEff2012 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        weight;
   Double_t        weight2;
   Double_t        weight3;
   Double_t        scanCrossSection;
   Double_t        scanCrossSectionPlus;
   Double_t        scanCrossSectionMinus;
   ULong64_t       runNumber;
   ULong64_t       lumiSection;
   ULong64_t       eventNumber;
   Int_t           m0;
   Int_t           m12;
   Int_t           mIntermediate;
   Int_t           ttbarDecayCode;
   Float_t         PUweight;
   Float_t PUweightSystVar;
   Int_t maxTOBTECjetDeltaMult;
   Int_t nIsoTracks15_005_03;
   Float_t         hltHTeff;
   Float_t         hltMHTeff;
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
   Bool_t          cutPV;
   Bool_t          cutTrigger;
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
   Int_t           nGoodPV;
   Int_t           SUSY_nb;
   Int_t           SUSY_process;
   Float_t         SUSY_recoilPt;
   Float_t SUSY_ISRweight;
   Float_t SUSY_ISRweightSystUp;
   Float_t SUSY_ISRweightSystDown;
   Int_t           nCorrectRecoStop;
   Int_t           njets;
   Int_t           njets30;
   Int_t           nbjets;
   Int_t           nbjets30;
   Int_t           ntruebjets;
   Int_t           nElectrons;
   Int_t           nMuons;
   Int_t           nTaus;
   Int_t           nElectrons5;
   Int_t           nMuons5;
   Int_t           nElectrons15;
   Int_t           nMuons15;
   Int_t           nElectrons20;
   Int_t           nMuons20;
   Float_t         bestZmass;
   Float_t         mjj1;
   Float_t         mjj2;
   Float_t         mjjdiff;
   Float_t         mjj1_5;
   Float_t         mjj2_5;
   Float_t         mjjdiff_5;
   Float_t         mjjb1;
   Float_t         mjjb2;
   Float_t         topPT1;
   Float_t         topPT2;
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
   Bool_t          pass_utilityPrescaleModuleHLT;
   Bool_t          pass_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45;
   Bool_t          pass_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45;
   Bool_t          pass_DiCentralPFJet30_PFMET80;
   Bool_t          pass_DiCentralPFJet30_PFMET80_BTagCSV07;
   Bool_t          pass_DiCentralPFJet50_PFMET80;
   Bool_t          pass_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80;
   Bool_t          pass_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Bool_t          pass_HT200;
   Bool_t          pass_HT250;
   Bool_t          pass_HT250_AlphaT0p55;
   Bool_t          pass_HT300;
   Bool_t          pass_HT300_AlphaT0p53;
   Bool_t          pass_IsoMu24;
   Bool_t          pass_IsoMu24_eta2p1;
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
   Bool_t          pass_PFMET150;
   Bool_t          pass_PFNoPUHT350;
   Bool_t          pass_PFNoPUHT350_Mu15_PFMET45;
   Bool_t          pass_PFNoPUHT350_PFMET100;
   Bool_t          pass_PFNoPUHT650;
   Bool_t          pass_Photon135;
   Bool_t          pass_Photon150;
   Bool_t          pass_QuadJet80;
   Bool_t          pass_SixJet45;
   Float_t         HT;
   Float_t         HT30;
   Float_t         ST;
   Float_t         STeff;
   Float_t         MET;
   Float_t         METphi;
   Float_t         MHT;
   Float_t         METsig;
   Float_t         METsig00;
   Float_t         METsig10;
   Float_t         METsig11;
   Float_t         caloMET;
   Float_t         rawPFMET;
   Float_t         rawPFMETphi;
   Float_t         bestWMass;
   Float_t         bestTopMass;
   Float_t         topCosHel;
   Float_t         WCosHel;
   Float_t         MT_b;
   Float_t MT_bestCSV;
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
   Float_t         minDeltaPhiAllb30;
   Float_t         deltaPhib1;
   Float_t         deltaPhib2;
   Float_t         deltaPhib3;
   Float_t         minDeltaPhiMETMuonsAll;
   Float_t         minDeltaPhiN_lostJet;
   Float_t         deltaPhiN1_lostJet;
   Float_t         deltaPhiN2_lostJet;
   Float_t         deltaPhiN3_lostJet;
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
   Float_t         elephi1;
   Float_t         eleeta1;
   Float_t         elecharge1;
   Float_t         muonpt1;
   Float_t         muonphi1;
   Float_t         muoneta1;
   Float_t         muoncharge1;
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
   TBranch        *b_weight2;   //!
   TBranch        *b_weight3;   //!
   TBranch        *b_scanCrossSection;   //!
   TBranch        *b_scanCrossSectionPlus;   //!
   TBranch        *b_scanCrossSectionMinus;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_lumiSection;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_m0;   //!
   TBranch        *b_m12;   //!
   TBranch        *b_mIntermediate;   //!
   TBranch* b_maxTOBTECjetDeltaMult;
   TBranch *b_nIsoTracks15_005_03;
   TBranch        *b_ttbarDecayCode;   //!
   TBranch        *b_PUweight;   //!
   TBranch *b_PUweightSystVar;
   TBranch        *b_hltHTeff;   //!
   TBranch        *b_hltMHTeff;   //!
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
   TBranch        *b_cutPV;   //!
   TBranch        *b_cutTrigger;   //!
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
   TBranch        *b_nGoodPV;   //!
   TBranch        *b_SUSY_nb;   //!
   TBranch        *b_SUSY_process;   //!
   TBranch        *b_SUSY_recoilPt;   //!
   TBranch        *b_SUSY_ISRweight;   //!
   TBranch        *b_nCorrectRecoStop;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_njets30;   //!
   TBranch        *b_nbjets;   //!
   TBranch        *b_nbjets30;   //!
   TBranch        *b_ntruebjets;   //!
   TBranch        *b_nElectrons;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_nTaus;   //!
   TBranch        *b_nElectrons5;   //!
   TBranch        *b_nMuons5;   //!
   TBranch        *b_nElectrons15;   //!
   TBranch        *b_nMuons15;   //!
   TBranch        *b_nElectrons20;   //!
   TBranch        *b_nMuons20;   //!
   TBranch        *b_bestZmass;   //!
   TBranch        *b_mjj1;   //!
   TBranch        *b_mjj2;   //!
   TBranch        *b_mjjdiff;   //!
   TBranch        *b_mjj1_5;   //!
   TBranch        *b_mjj2_5;   //!
   TBranch        *b_mjjdiff_5;   //!
   TBranch        *b_mjjb1;   //!
   TBranch        *b_mjjb2;   //!
   TBranch        *b_topPT1;   //!
   TBranch        *b_topPT2;   //!
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
   TBranch        *b_pass_utilityPrescaleModuleHLT;   //!
   TBranch        *b_pass_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45;   //!
   TBranch        *b_pass_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45;   //!
   TBranch        *b_pass_DiCentralPFJet30_PFMET80;   //!
   TBranch        *b_pass_DiCentralPFJet30_PFMET80_BTagCSV07;   //!
   TBranch        *b_pass_DiCentralPFJet50_PFMET80;   //!
   TBranch        *b_pass_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80;   //!
   TBranch        *b_pass_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_pass_HT200;   //!
   TBranch        *b_pass_HT250;   //!
   TBranch        *b_pass_HT250_AlphaT0p55;   //!
   TBranch        *b_pass_HT300;   //!
   TBranch        *b_pass_HT300_AlphaT0p53;   //!
   TBranch        *b_pass_IsoMu24;   //!
   TBranch        *b_pass_IsoMu24_eta2p1;   //!
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
   TBranch        *b_pass_PFMET150;   //!
   TBranch        *b_pass_PFNoPUHT350;   //!
   TBranch        *b_pass_PFNoPUHT350_Mu15_PFMET45;   //!
   TBranch        *b_pass_PFNoPUHT350_PFMET100;   //!
   TBranch        *b_pass_PFNoPUHT650;   //!
   TBranch        *b_pass_Photon135;   //!
   TBranch        *b_pass_Photon150;   //!
   TBranch        *b_pass_QuadJet80;   //!
   TBranch        *b_pass_SixJet45;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_HT30;   //!
   TBranch        *b_ST;   //!
   TBranch        *b_STeff;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METphi;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_METsig;   //!
   TBranch        *b_METsig00;   //!
   TBranch        *b_METsig10;   //!
   TBranch        *b_METsig11;   //!
   TBranch        *b_caloMET;   //!
   TBranch        *b_rawPFMET;   //!
   TBranch        *b_rawPFMETphi;   //!
   TBranch        *b_bestWMass;   //!
   TBranch        *b_bestTopMass;   //!
   TBranch        *b_topCosHel;   //!
   TBranch        *b_WCosHel;   //!
   TBranch        *b_MT_b;   //!
   TBranch *b_MT_bestCSV;
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
   TBranch        *b_minDeltaPhiAllb30;   //!
   TBranch        *b_deltaPhib1;   //!
   TBranch        *b_deltaPhib2;   //!
   TBranch        *b_deltaPhib3;   //!
   TBranch        *b_minDeltaPhiMETMuonsAll;   //!
   TBranch        *b_minDeltaPhiN_lostJet;   //!
   TBranch        *b_deltaPhiN1_lostJet;   //!
   TBranch        *b_deltaPhiN2_lostJet;   //!
   TBranch        *b_deltaPhiN3_lostJet;   //!
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

   signalEff2012(TString path, TString filestub, bool joinbtagbins, bool usebtagsf, bool dopdfs, bool pusyst, int minnjets, int isrmode=99); //BEGIN END jmt
   virtual ~signalEff2012();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   //BEGIN jmt
   bool dimensionIsFilled(TString dim,TH3D* h3,int & filledbin);

   TH2D* scanSMSngen_;
   bool isPMSSM_;
   bool is3D_;
   TString unfilledDimension_;
   vector<TH2D*> scanProcessTotals_;
   TString filestub_;
   bool joinbtagbins_;
   bool usebtagsf_;
   bool dopdfs_;
   bool pusyst_;
   int minnjets_;
   IsrMode theIsrMode_;
   //END jmt
};

#endif

#ifdef signalEff2012_cxx
//BEGIN jmt mod
signalEff2012::signalEff2012(TString path, TString filestub,bool joinbtagbins, bool usebtagsf, bool dopdfs, bool pusyst, int minnjets, int isrmode) : fChain(0) 
													    , scanSMSngen_(0)
													    ,isPMSSM_(false)
													    ,is3D_(false)
												            ,unfilledDimension_("")
													    ,filestub_(filestub)
													    ,joinbtagbins_(joinbtagbins)
													    ,usebtagsf_(usebtagsf)
													    ,dopdfs_(dopdfs)
													    ,pusyst_(pusyst)
													    ,minnjets_(minnjets)
													    ,theIsrMode_(kNoIsrWeight)
			    //END jmt mod
{
 //BEGIN jmt
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
  std::cout<<" minnjets = "<<minnjets<<std::endl;
  std::cout<<" isrmode = "<<int(theIsrMode_)<<std::endl;
  //END jmt


// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename); //BEGIN END jmt
      if (!f || !f->IsOpen()) {
	f = new TFile(filename); //BEGIN END jmt
      }
      f->GetObject("reducedTree",tree);
      scanSMSngen_ = (TH2D*) f->Get("scanSMSngen"); //BEGIN  jmt mod
      TH1D*   scanpMSSMngen = (TH1D*) f->Get("scanpMSSMngen");
      TH3D*   scanSMSngen3D = (TH3D*) f->Get("scanSMSngen3D");

      //dimensions are:
      //                x == mParent
      //                y == mLSP
      //                z == mIntermediate
      if ( filestub.Contains("T5tttt") ) {
	is3D_=true;
	unfilledDimension_="y"; //T5tttt fixes the LSP mass
      }
      else if ( filestub.Contains("T1t1t") ) {
	is3D_=true;
	unfilledDimension_="x"; //T5tttt fixes the LSP mass
      }

      if (scanpMSSMngen!=0) {
	assert(scanSMSngen_==0); //only one of these should exist.
	//the coding is more generic if we treat pMSSM as a 1d SMS.
	scanSMSngen_ = new TH2D("scanSMSngen","scan sms ngen",scanpMSSMngen->GetNbinsX(),scanpMSSMngen->GetXaxis()->GetXmin(),scanpMSSMngen->GetXaxis()->GetXmax(),1,-1,1);
	//now copy the contents of the pMSSM histo into this new histo
	for ( int ix = 1; ix<= scanpMSSMngen->GetNbinsX(); ix++) {
	  scanSMSngen_->SetBinContent(ix, 1, scanpMSSMngen->GetBinContent(ix));
	}
	isPMSSM_=true;
      }
      else if (is3D_ ) {
	//all of our scans are actually 2D. Map the 3d scan onto 2d.
	//advantage -- code stays unchanged after this point
	//disadvantage, need to manually keep track of 3rd mass (probably not a big deal)

	//i was going to try to automatically find the unfilled dimension
	//but since i process the scans split into strips of gluino mass, this is not possible

	int filledbin;
	bool isfilled = dimensionIsFilled(unfilledDimension_,scanSMSngen3D,filledbin);
	assert( !isfilled);

	TAxis * xaxis=0;
	TAxis * yaxis=0;
	if ( unfilledDimension_=="z") {
	  xaxis = scanSMSngen3D->GetXaxis();
	  yaxis = scanSMSngen3D->GetYaxis();
	}
	else if (unfilledDimension_=="y") {
	  xaxis = scanSMSngen3D->GetXaxis();
	  yaxis = scanSMSngen3D->GetZaxis();
	}
	else if (unfilledDimension_=="x") {
	  xaxis = scanSMSngen3D->GetZaxis();
	  yaxis = scanSMSngen3D->GetYaxis();
	}
	else assert(0);

	scanSMSngen_=new TH2D("scanSMSngen","scan sms ngen", xaxis->GetNbins(),xaxis->GetXmin(),xaxis->GetXmax() ,yaxis->GetNbins(), yaxis->GetXmin(),yaxis->GetXmax()  );

	for ( int ix = 1; ix<=xaxis->GetNbins(); ix++) {
	  for ( int iy = 1; iy<=yaxis->GetNbins() ; iy++) {
	    double val=0;
	    if (unfilledDimension_=="z") val =scanSMSngen3D->GetBinContent( ix,iy,filledbin);
	    else  if (unfilledDimension_=="y") val= scanSMSngen3D->GetBinContent( ix,filledbin,iy);
	    else  if (unfilledDimension_=="x") val=scanSMSngen3D->GetBinContent( filledbin,iy,ix);

	    scanSMSngen_->SetBinContent(ix,iy, val);
	  }
	}
      }
      if (scanSMSngen_==0) {
	std::cout<<" WARNING -- scanSMSngen not found (only ok for non-SMS)"<<std::endl;
	scanSMSngen_ = new TH2D("scanSMSngen","scan sms ngen",1,-1,1,1,-1,1);
      }
      //need to loop over the keys in the file and load all scan process totals histos
      //(this is now deprecated...but leave it in place just in case)
      for (int ih = 0; ih<f->GetListOfKeys()->GetEntries(); ih++) {
	TString histname = f->GetListOfKeys()->At(ih)->GetName();
	if (!histname.BeginsWith("scanProcessTotals")) continue;
	scanProcessTotals_.push_back( (TH2D*) f->Get(histname));
      }

      //END jmt mod

   }
   Init(tree);
}

signalEff2012::~signalEff2012()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t signalEff2012::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t signalEff2012::LoadTree(Long64_t entry)
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

void signalEff2012::Init(TTree *tree)
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
   fChain->SetBranchAddress("scanCrossSectionPlus", &scanCrossSectionPlus, &b_scanCrossSectionPlus);
   fChain->SetBranchAddress("scanCrossSectionMinus", &scanCrossSectionMinus, &b_scanCrossSectionMinus);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("m0", &m0, &b_m0);
   fChain->SetBranchAddress("m12", &m12, &b_m12);
   fChain->SetBranchAddress("mIntermediate", &mIntermediate, &b_mIntermediate);
   fChain->SetBranchAddress("nIsoTracks15_005_03", &nIsoTracks15_005_03, &b_nIsoTracks15_005_03);
   fChain->SetBranchAddress("maxTOBTECjetDeltaMult", &maxTOBTECjetDeltaMult, &b_maxTOBTECjetDeltaMult);
   fChain->SetBranchAddress("ttbarDecayCode", &ttbarDecayCode, &b_ttbarDecayCode);
   fChain->SetBranchAddress("PUweight", &PUweight, &b_PUweight);
   fChain->SetBranchAddress("PUweightSystVar", &PUweightSystVar, &b_PUweightSystVar);
   fChain->SetBranchAddress("hltHTeff", &hltHTeff, &b_hltHTeff);
   fChain->SetBranchAddress("hltMHTeff", &hltMHTeff, &b_hltMHTeff);
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
   fChain->SetBranchAddress("cutPV", &cutPV, &b_cutPV);
   fChain->SetBranchAddress("cutTrigger", &cutTrigger, &b_cutTrigger);
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
   fChain->SetBranchAddress("nGoodPV", &nGoodPV, &b_nGoodPV);
   fChain->SetBranchAddress("SUSY_nb", &SUSY_nb, &b_SUSY_nb);
   fChain->SetBranchAddress("SUSY_process", &SUSY_process, &b_SUSY_process);
   fChain->SetBranchAddress("SUSY_recoilPt", &SUSY_recoilPt, &b_SUSY_recoilPt);
   fChain->SetBranchAddress("SUSY_ISRweight", &SUSY_ISRweight, &b_SUSY_ISRweight);
   fChain->SetBranchAddress("nCorrectRecoStop", &nCorrectRecoStop, &b_nCorrectRecoStop);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("njets30", &njets30, &b_njets30);
   fChain->SetBranchAddress("nbjets", &nbjets, &b_nbjets);
   fChain->SetBranchAddress("nbjets30", &nbjets30, &b_nbjets30);
   fChain->SetBranchAddress("ntruebjets", &ntruebjets, &b_ntruebjets);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   //   fChain->SetBranchAddress("nTaus", &nTaus, &b_nTaus);
   fChain->SetBranchAddress("nElectrons5", &nElectrons5, &b_nElectrons5);
   fChain->SetBranchAddress("nMuons5", &nMuons5, &b_nMuons5);
   fChain->SetBranchAddress("nElectrons15", &nElectrons15, &b_nElectrons15);
   fChain->SetBranchAddress("nMuons15", &nMuons15, &b_nMuons15);
   fChain->SetBranchAddress("nElectrons20", &nElectrons20, &b_nElectrons20);
   fChain->SetBranchAddress("nMuons20", &nMuons20, &b_nMuons20);
   fChain->SetBranchAddress("bestZmass", &bestZmass, &b_bestZmass);
   fChain->SetBranchAddress("mjj1", &mjj1, &b_mjj1);
   fChain->SetBranchAddress("mjj2", &mjj2, &b_mjj2);
   fChain->SetBranchAddress("mjjdiff", &mjjdiff, &b_mjjdiff);
   fChain->SetBranchAddress("mjj1_5", &mjj1_5, &b_mjj1_5);
   fChain->SetBranchAddress("mjj2_5", &mjj2_5, &b_mjj2_5);
   fChain->SetBranchAddress("mjjdiff_5", &mjjdiff_5, &b_mjjdiff_5);
   fChain->SetBranchAddress("mjjb1", &mjjb1, &b_mjjb1);
   fChain->SetBranchAddress("mjjb2", &mjjb2, &b_mjjb2);
   fChain->SetBranchAddress("topPT1", &topPT1, &b_topPT1);
   fChain->SetBranchAddress("topPT2", &topPT2, &b_topPT2);
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
   fChain->SetBranchAddress("pass_utilityPrescaleModuleHLT", &pass_utilityPrescaleModuleHLT, &b_pass_utilityPrescaleModuleHLT);
   fChain->SetBranchAddress("pass_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45", &pass_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45, &b_pass_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45);
   fChain->SetBranchAddress("pass_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45", &pass_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45, &b_pass_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45);
   fChain->SetBranchAddress("pass_DiCentralPFJet30_PFMET80", &pass_DiCentralPFJet30_PFMET80, &b_pass_DiCentralPFJet30_PFMET80);
   fChain->SetBranchAddress("pass_DiCentralPFJet30_PFMET80_BTagCSV07", &pass_DiCentralPFJet30_PFMET80_BTagCSV07, &b_pass_DiCentralPFJet30_PFMET80_BTagCSV07);
   fChain->SetBranchAddress("pass_DiCentralPFJet50_PFMET80", &pass_DiCentralPFJet50_PFMET80, &b_pass_DiCentralPFJet50_PFMET80);
   fChain->SetBranchAddress("pass_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80", &pass_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80, &b_pass_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80);
   fChain->SetBranchAddress("pass_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &pass_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_pass_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("pass_HT200", &pass_HT200, &b_pass_HT200);
   fChain->SetBranchAddress("pass_HT250", &pass_HT250, &b_pass_HT250);
   fChain->SetBranchAddress("pass_HT250_AlphaT0p55", &pass_HT250_AlphaT0p55, &b_pass_HT250_AlphaT0p55);
   fChain->SetBranchAddress("pass_HT300", &pass_HT300, &b_pass_HT300);
   fChain->SetBranchAddress("pass_HT300_AlphaT0p53", &pass_HT300_AlphaT0p53, &b_pass_HT300_AlphaT0p53);
   fChain->SetBranchAddress("pass_IsoMu24", &pass_IsoMu24, &b_pass_IsoMu24);
   fChain->SetBranchAddress("pass_IsoMu24_eta2p1", &pass_IsoMu24_eta2p1, &b_pass_IsoMu24_eta2p1);
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
   fChain->SetBranchAddress("pass_PFMET150", &pass_PFMET150, &b_pass_PFMET150);
   fChain->SetBranchAddress("pass_PFNoPUHT350", &pass_PFNoPUHT350, &b_pass_PFNoPUHT350);
   fChain->SetBranchAddress("pass_PFNoPUHT350_Mu15_PFMET45", &pass_PFNoPUHT350_Mu15_PFMET45, &b_pass_PFNoPUHT350_Mu15_PFMET45);
   fChain->SetBranchAddress("pass_PFNoPUHT350_PFMET100", &pass_PFNoPUHT350_PFMET100, &b_pass_PFNoPUHT350_PFMET100);
   fChain->SetBranchAddress("pass_PFNoPUHT650", &pass_PFNoPUHT650, &b_pass_PFNoPUHT650);
   fChain->SetBranchAddress("pass_Photon135", &pass_Photon135, &b_pass_Photon135);
   fChain->SetBranchAddress("pass_Photon150", &pass_Photon150, &b_pass_Photon150);
   fChain->SetBranchAddress("pass_QuadJet80", &pass_QuadJet80, &b_pass_QuadJet80);
   fChain->SetBranchAddress("pass_SixJet45", &pass_SixJet45, &b_pass_SixJet45);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("HT30", &HT30, &b_HT30);
   fChain->SetBranchAddress("ST", &ST, &b_ST);
   fChain->SetBranchAddress("STeff", &STeff, &b_STeff);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("METphi", &METphi, &b_METphi);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("METsig", &METsig, &b_METsig);
   fChain->SetBranchAddress("METsig00", &METsig00, &b_METsig00);
   fChain->SetBranchAddress("METsig10", &METsig10, &b_METsig10);
   fChain->SetBranchAddress("METsig11", &METsig11, &b_METsig11);
   fChain->SetBranchAddress("caloMET", &caloMET, &b_caloMET);
   fChain->SetBranchAddress("rawPFMET", &rawPFMET, &b_rawPFMET);
   fChain->SetBranchAddress("rawPFMETphi", &rawPFMETphi, &b_rawPFMETphi);
   fChain->SetBranchAddress("bestWMass", &bestWMass, &b_bestWMass);
   fChain->SetBranchAddress("bestTopMass", &bestTopMass, &b_bestTopMass);
   fChain->SetBranchAddress("topCosHel", &topCosHel, &b_topCosHel);
   fChain->SetBranchAddress("WCosHel", &WCosHel, &b_WCosHel);
   fChain->SetBranchAddress("MT_b", &MT_b, &b_MT_b);
   fChain->SetBranchAddress("MT_bestCSV", &MT_bestCSV, &b_MT_bestCSV);
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
   fChain->SetBranchAddress("minDeltaPhiAllb30", &minDeltaPhiAllb30, &b_minDeltaPhiAllb30);
   fChain->SetBranchAddress("deltaPhib1", &deltaPhib1, &b_deltaPhib1);
   fChain->SetBranchAddress("deltaPhib2", &deltaPhib2, &b_deltaPhib2);
   fChain->SetBranchAddress("deltaPhib3", &deltaPhib3, &b_deltaPhib3);
   fChain->SetBranchAddress("minDeltaPhiMETMuonsAll", &minDeltaPhiMETMuonsAll, &b_minDeltaPhiMETMuonsAll);
   fChain->SetBranchAddress("minDeltaPhiN_lostJet", &minDeltaPhiN_lostJet, &b_minDeltaPhiN_lostJet);
   fChain->SetBranchAddress("deltaPhiN1_lostJet", &deltaPhiN1_lostJet, &b_deltaPhiN1_lostJet);
   fChain->SetBranchAddress("deltaPhiN2_lostJet", &deltaPhiN2_lostJet, &b_deltaPhiN2_lostJet);
   fChain->SetBranchAddress("deltaPhiN3_lostJet", &deltaPhiN3_lostJet, &b_deltaPhiN3_lostJet);
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
   fChain->SetBranchAddress("eleRelIso", &eleRelIso, &b_eleRelIso);
   fChain->SetBranchAddress("muonRelIso", &muonRelIso, &b_muonRelIso);
   /*
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
   */
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

Bool_t signalEff2012::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void signalEff2012::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t signalEff2012::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef signalEff2012_cxx
