//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 11 19:49:45 2011 by ROOT version 5.22/00d
// from TChain reducedTree/
//////////////////////////////////////////////////////////

#ifndef myCutflow_h
#define myCutflow_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
using namespace std;

class myCutflow {
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
   Float_t         btagIPweight;
   Float_t         PUweight;
   Float_t         hltHTeff;
   Float_t         pdfWeightsCTEQ[45];
   Float_t         pdfWeightsMSTW[41];
   Float_t         pdfWeightsNNPDF[100];
   Float_t         prob0;
   Float_t         probge1;
   Float_t         prob1;
   Float_t         probge2;
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
   Bool_t          passCleaning;
   Int_t           PBNRcode;
   Int_t           nGoodPV;
   Int_t           SUSY_nb;
   Int_t           SUSY_process;
   Int_t           njets;
   Int_t           njets30;
   Int_t           nbjets;
   Int_t           nElectrons;
   Int_t           nMuons;
   Int_t           nbjetsSSVM;
   Int_t           nbjetsSSVHPT;
   Int_t           nbjetsTCHPT;
   Int_t           nbjetsTCHET;
   Int_t           nbjetsTCHPM;
   Int_t           nbjetsCSVM;
   Bool_t          isRealData;
   Bool_t          pass_utilityHLT;
   UInt_t          prescaleUtilityHLT;
   UInt_t          versionUtilityHLT;
   UInt_t          pass_utilityHLT_HT300_CentralJet30_BTagIP;
   UInt_t          prescale_utilityHLT_HT300_CentralJet30_BTagIP;
   Float_t         HT;
   Float_t         MET;
   Float_t         METsig;
   Float_t         METphi;
   Float_t         MHT;
   Float_t         correctedMET;
   Float_t         correctedMETphi;
   Float_t         bestWMass;
   Float_t         bestTopMass;
   Float_t         topCosHel;
   Float_t         WCosHel;
   Float_t         MT_Wlep;
   Float_t         minDeltaPhi;
   Float_t         minDeltaPhiAll;
   Float_t         minDeltaPhiAll30;
   Float_t         minDeltaPhi30_eta5_noIdAll;
   Float_t         minDeltaPhiMetTau;
   Float_t         maxDeltaPhi;
   Float_t         maxDeltaPhiAll;
   Float_t         maxDeltaPhiAll30;
   Float_t         maxDeltaPhi30_eta5_noIdAll;
   Float_t         sumDeltaPhi;
   Float_t         diffDeltaPhi;
   Float_t         minDeltaPhiN;
   Float_t         deltaPhiN1;
   Float_t         deltaPhiN2;
   Float_t         deltaPhiN3;
   Float_t         minDeltaPhiN_otherEta5;
   Float_t         minDeltaPhiN_otherEta5idNo;
   Float_t         minDeltaPhiN_mainPt30_otherEta5idNo;
   Float_t         minDeltaPhiN_mainPt30Eta5_otherEta5idNo;
   Float_t         minDeltaPhiN_mainPt30Eta5idNo_otherEta5idNo;
   Float_t         minDeltaPhiK;
   Float_t         minDeltaPhiK_otherEta5;
   Float_t         minDeltaPhiK_otherEta5idNo;
   Float_t         minDeltaPhiK_mainPt30_otherEta5idNo;
   Float_t         minDeltaPhiK_mainPt30Eta5_otherEta5idNo;
   Float_t         minDeltaPhiK_mainPt30Eta5idNo_otherEta5idNo;
   Float_t         minDeltaPhiN_DJR;
   Float_t         minDeltaPhiN_DJR_otherEta5;
   Float_t         minDeltaPhiK_DJR;
   Float_t         minDeltaPhiK_DJR_otherEta5;
   Float_t         minDeltaPhiN_lostJet;
   Float_t         deltaPhiN1_lostJet;
   Float_t         deltaPhiN2_lostJet;
   Float_t         deltaPhiN3_lostJet;
   Float_t         jetpt1;
   Float_t         jeteta1;
   Float_t         jetphi1;
   Float_t         jetenergy1;
   Float_t         jetpt2;
   Float_t         jeteta2;
   Float_t         jetphi2;
   Float_t         jetenergy2;
   Float_t         jetpt3;
   Float_t         jeteta3;
   Float_t         jetphi3;
   Float_t         jetenergy3;
   Float_t         bjetpt1;
   Float_t         bjeteta1;
   Float_t         bjetphi1;
   Float_t         bjetenergy1;
   Float_t         bjetpt2;
   Float_t         bjeteta2;
   Float_t         bjetphi2;
   Float_t         bjetenergy2;
   Float_t         bjetpt3;
   Float_t         bjeteta3;
   Float_t         bjetphi3;
   Float_t         bjetenergy3;
   Float_t         eleet1;
   Float_t         muonpt1;
   Float_t         lambda1_allJets;
   Float_t         lambda2_allJets;
   Float_t         determinant_allJets;
   Float_t         lambda1_allJetsPlusMET;
   Float_t         lambda2_allJetsPlusMET;
   Float_t         determinant_allJetsPlusMET;
   Float_t         lambda1_topThreeJets;
   Float_t         lambda2_topThreeJets;
   Float_t         determinant_topThreeJets;
   Float_t         lambda1_topThreeJetsPlusMET;
   Float_t         lambda2_topThreeJetsPlusMET;
   Float_t         determinant_topThreeJetsPlusMET;
   Float_t         transverseThrust;
   Float_t         transverseThrustPhi;
   Float_t         transverseThrustWithMET;
   Float_t         transverseThrustWithMETPhi;
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
   TBranch        *b_btagIPweight;   //!
   TBranch        *b_PUweight;   //!
   TBranch        *b_hltHTeff;   //!
   TBranch        *b_pdfWeightsCTEQ;   //!
   TBranch        *b_pdfWeightsMSTW;   //!
   TBranch        *b_pdfWeightsNNPDF;   //!
   TBranch        *b_prob0;   //!
   TBranch        *b_probge1;   //!
   TBranch        *b_prob1;   //!
   TBranch        *b_probge2;   //!
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
   TBranch        *b_passCleaning;   //!
   TBranch        *b_PBNRcode;   //!
   TBranch        *b_nGoodPV;   //!
   TBranch        *b_SUSY_nb;   //!
   TBranch        *b_SUSY_process;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_njets30;   //!
   TBranch        *b_nbjets;   //!
   TBranch        *b_nElectrons;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_nbjetsSSVM;   //!
   TBranch        *b_nbjetsSSVHPT;   //!
   TBranch        *b_nbjetsTCHPT;   //!
   TBranch        *b_nbjetsTCHET;   //!
   TBranch        *b_nbjetsTCHPM;   //!
   TBranch        *b_nbjetsCSVM;   //!
   TBranch        *b_isRealData;   //!
   TBranch        *b_pass_utilityHLT;   //!
   TBranch        *b_prescaleUtilityHLT;   //!
   TBranch        *b_versionUtilityHLT;   //!
   TBranch        *b_pass_utilityHLT_HT300_CentralJet30_BTagIP;   //!
   TBranch        *b_prescale_utilityHLT_HT300_CentralJet30_BTagIP;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METsig;   //!
   TBranch        *b_METphi;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_correctedMET;   //!
   TBranch        *b_correctedMETphi;   //!
   TBranch        *b_bestWMass;   //!
   TBranch        *b_bestTopMass;   //!
   TBranch        *b_topCosHel;   //!
   TBranch        *b_WCosHel;   //!
   TBranch        *b_MT_Wlep;   //!
   TBranch        *b_minDeltaPhi;   //!
   TBranch        *b_minDeltaPhiAll;   //!
   TBranch        *b_minDeltaPhiAll30;   //!
   TBranch        *b_minDeltaPhi30_eta5_noIdAll;   //!
   TBranch        *b_minDeltaPhiMetTau;   //!
   TBranch        *b_maxDeltaPhi;   //!
   TBranch        *b_maxDeltaPhiAll;   //!
   TBranch        *b_maxDeltaPhiAll30;   //!
   TBranch        *b_maxDeltaPhi30_eta5_noIdAll;   //!
   TBranch        *b_sumDeltaPhi;   //!
   TBranch        *b_diffDeltaPhi;   //!
   TBranch        *b_minDeltaPhiN;   //!
   TBranch        *b_deltaPhiN1;   //!
   TBranch        *b_deltaPhiN2;   //!
   TBranch        *b_deltaPhiN3;   //!
   TBranch        *b_minDeltaPhiN_otherEta5;   //!
   TBranch        *b_minDeltaPhiN_otherEta5idNo;   //!
   TBranch        *b_minDeltaPhiN_mainPt30_otherEta5idNo;   //!
   TBranch        *b_minDeltaPhiN_mainPt30Eta5_otherEta5idNo;   //!
   TBranch        *b_minDeltaPhiN_mainPt30Eta5idNo_otherEta5idNo;   //!
   TBranch        *b_minDeltaPhiK;   //!
   TBranch        *b_minDeltaPhiK_otherEta5;   //!
   TBranch        *b_minDeltaPhiK_otherEta5idNo;   //!
   TBranch        *b_minDeltaPhiK_mainPt30_otherEta5idNo;   //!
   TBranch        *b_minDeltaPhiK_mainPt30Eta5_otherEta5idNo;   //!
   TBranch        *b_minDeltaPhiK_mainPt30Eta5idNo_otherEta5idNo;   //!
   TBranch        *b_minDeltaPhiN_DJR;   //!
   TBranch        *b_minDeltaPhiN_DJR_otherEta5;   //!
   TBranch        *b_minDeltaPhiK_DJR;   //!
   TBranch        *b_minDeltaPhiK_DJR_otherEta5;   //!
   TBranch        *b_minDeltaPhiN_lostJet;   //!
   TBranch        *b_deltaPhiN1_lostJet;   //!
   TBranch        *b_deltaPhiN2_lostJet;   //!
   TBranch        *b_deltaPhiN3_lostJet;   //!
   TBranch        *b_jetpt1;   //!
   TBranch        *b_jeteta1;   //!
   TBranch        *b_jetphi1;   //!
   TBranch        *b_jetenergy1;   //!
   TBranch        *b_jetpt2;   //!
   TBranch        *b_jeteta2;   //!
   TBranch        *b_jetphi2;   //!
   TBranch        *b_jetenergy2;   //!
   TBranch        *b_jetpt3;   //!
   TBranch        *b_jeteta3;   //!
   TBranch        *b_jetphi3;   //!
   TBranch        *b_jetenergy3;   //!
   TBranch        *b_bjetpt1;   //!
   TBranch        *b_bjeteta1;   //!
   TBranch        *b_bjetphi1;   //!
   TBranch        *b_bjetenergy1;   //!
   TBranch        *b_bjetpt2;   //!
   TBranch        *b_bjeteta2;   //!
   TBranch        *b_bjetphi2;   //!
   TBranch        *b_bjetenergy2;   //!
   TBranch        *b_bjetpt3;   //!
   TBranch        *b_bjeteta3;   //!
   TBranch        *b_bjetphi3;   //!
   TBranch        *b_bjetenergy3;   //!
   TBranch        *b_eleet1;   //!
   TBranch        *b_muonpt1;   //!
   TBranch        *b_lambda1_allJets;   //!
   TBranch        *b_lambda2_allJets;   //!
   TBranch        *b_determinant_allJets;   //!
   TBranch        *b_lambda1_allJetsPlusMET;   //!
   TBranch        *b_lambda2_allJetsPlusMET;   //!
   TBranch        *b_determinant_allJetsPlusMET;   //!
   TBranch        *b_lambda1_topThreeJets;   //!
   TBranch        *b_lambda2_topThreeJets;   //!
   TBranch        *b_determinant_topThreeJets;   //!
   TBranch        *b_lambda1_topThreeJetsPlusMET;   //!
   TBranch        *b_lambda2_topThreeJetsPlusMET;   //!
   TBranch        *b_determinant_topThreeJetsPlusMET;   //!
   TBranch        *b_transverseThrust;   //!
   TBranch        *b_transverseThrustPhi;   //!
   TBranch        *b_transverseThrustWithMET;   //!
   TBranch        *b_transverseThrustWithMETPhi;   //!
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

   TString Event_; //used to control the data sample

   myCutflow(TString sampleName, TTree *tree=0);
   virtual ~myCutflow();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
//   virtual void     Loop(TString sampleName, bool isTightSelection, bool writeFiles);
   virtual void     Loop(bool isTightSelection, bool writeFiles);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   TString getDataSample(){return Event_;}
   vector<TString> getFileNames(TString sampleName, bool isTightSelection);// vector<TString> &fileNames);
   vector<ofstream*> getFiles(TString sampleName, bool isTightSelection);
//   int writeToFile(TString sampleName, bool isTightSelection);


};

#endif

#ifdef myCutflow_cxx
myCutflow::myCutflow(TString sampleName, TTree *tree):Event_(sampleName)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f) {
         f = new TFile("Memory Directory");
         f->cd("Rint:/");
      }
      tree = (TTree*)gDirectory->Get("reducedTree");

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("reducedTree","");

      if (sampleName=="TTbarJets"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-35a/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.TTbarJets.root/reducedTree");
      }

      else if (sampleName=="SingleTop"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-35a/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.SingleTop.root/reducedTree");
      }

      else if (sampleName=="WJets"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-35a/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.WJets.root/reducedTree");
      }

      else if (sampleName=="ZJets"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-35a/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.ZJets.root/reducedTree");
      }

      else if (sampleName=="Zinvisible"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-35a/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.Zinvisible.root/reducedTree");
      }

      else if (sampleName=="VV"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-35a/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.VV.root/reducedTree");
      }

      else if (sampleName=="PythiaPUQCD"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-35a/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.PythiaPUQCD.root/reducedTree");
      }

      else if (sampleName=="LM9"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-35a/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.LM9.root/reducedTree");
      }

      else if (sampleName=="Data"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-35a/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.ht_run2011a_aug5rereco.root/reducedTree");
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-35a/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.ht_run2011a_may10rereco.root/reducedTree");
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-35a/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.ht_run2011a_promptrecov4_try3.root/reducedTree");
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-35a/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.ht_run2011a_promptrecov6.root/reducedTree");
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-35a/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.ht_run2011b_promptrecov1_oct14.root/reducedTree");
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-35a/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.ht_run2011b_promptrecov1_oct7.root/reducedTree");
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-35a/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.ht_run2011b_promptrecov1.root/reducedTree");

      }

      else cout<<"I couldn't find the sample you asked for"<<endl;
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

myCutflow::~myCutflow()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myCutflow::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t myCutflow::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void myCutflow::Init(TTree *tree)
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
   fChain->SetBranchAddress("btagIPweight", &btagIPweight, &b_btagIPweight);
   fChain->SetBranchAddress("PUweight", &PUweight, &b_PUweight);
   fChain->SetBranchAddress("hltHTeff", &hltHTeff, &b_hltHTeff);
   fChain->SetBranchAddress("pdfWeightsCTEQ", pdfWeightsCTEQ, &b_pdfWeightsCTEQ);
   fChain->SetBranchAddress("pdfWeightsMSTW", pdfWeightsMSTW, &b_pdfWeightsMSTW);
   fChain->SetBranchAddress("pdfWeightsNNPDF", pdfWeightsNNPDF, &b_pdfWeightsNNPDF);
   fChain->SetBranchAddress("prob0", &prob0, &b_prob0);
   fChain->SetBranchAddress("probge1", &probge1, &b_probge1);
   fChain->SetBranchAddress("prob1", &prob1, &b_prob1);
   fChain->SetBranchAddress("probge2", &probge2, &b_probge2);
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
   fChain->SetBranchAddress("passCleaning", &passCleaning, &b_passCleaning);
   fChain->SetBranchAddress("PBNRcode", &PBNRcode, &b_PBNRcode);
   fChain->SetBranchAddress("nGoodPV", &nGoodPV, &b_nGoodPV);
   fChain->SetBranchAddress("SUSY_nb", &SUSY_nb, &b_SUSY_nb);
   fChain->SetBranchAddress("SUSY_process", &SUSY_process, &b_SUSY_process);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("njets30", &njets30, &b_njets30);
   fChain->SetBranchAddress("nbjets", &nbjets, &b_nbjets);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   fChain->SetBranchAddress("nbjetsSSVM", &nbjetsSSVM, &b_nbjetsSSVM);
   fChain->SetBranchAddress("nbjetsSSVHPT", &nbjetsSSVHPT, &b_nbjetsSSVHPT);
   fChain->SetBranchAddress("nbjetsTCHPT", &nbjetsTCHPT, &b_nbjetsTCHPT);
   fChain->SetBranchAddress("nbjetsTCHET", &nbjetsTCHET, &b_nbjetsTCHET);
   fChain->SetBranchAddress("nbjetsTCHPM", &nbjetsTCHPM, &b_nbjetsTCHPM);
   fChain->SetBranchAddress("nbjetsCSVM", &nbjetsCSVM, &b_nbjetsCSVM);
   fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
   fChain->SetBranchAddress("pass_utilityHLT", &pass_utilityHLT, &b_pass_utilityHLT);
   fChain->SetBranchAddress("prescaleUtilityHLT", &prescaleUtilityHLT, &b_prescaleUtilityHLT);
   fChain->SetBranchAddress("versionUtilityHLT", &versionUtilityHLT, &b_versionUtilityHLT);
   fChain->SetBranchAddress("pass_utilityHLT_HT300_CentralJet30_BTagIP", &pass_utilityHLT_HT300_CentralJet30_BTagIP, &b_pass_utilityHLT_HT300_CentralJet30_BTagIP);
   fChain->SetBranchAddress("prescale_utilityHLT_HT300_CentralJet30_BTagIP", &prescale_utilityHLT_HT300_CentralJet30_BTagIP, &b_prescale_utilityHLT_HT300_CentralJet30_BTagIP);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("METsig", &METsig, &b_METsig);
   fChain->SetBranchAddress("METphi", &METphi, &b_METphi);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("correctedMET", &correctedMET, &b_correctedMET);
   fChain->SetBranchAddress("correctedMETphi", &correctedMETphi, &b_correctedMETphi);
   fChain->SetBranchAddress("bestWMass", &bestWMass, &b_bestWMass);
   fChain->SetBranchAddress("bestTopMass", &bestTopMass, &b_bestTopMass);
   fChain->SetBranchAddress("topCosHel", &topCosHel, &b_topCosHel);
   fChain->SetBranchAddress("WCosHel", &WCosHel, &b_WCosHel);
   fChain->SetBranchAddress("MT_Wlep", &MT_Wlep, &b_MT_Wlep);
   fChain->SetBranchAddress("minDeltaPhi", &minDeltaPhi, &b_minDeltaPhi);
   fChain->SetBranchAddress("minDeltaPhiAll", &minDeltaPhiAll, &b_minDeltaPhiAll);
   fChain->SetBranchAddress("minDeltaPhiAll30", &minDeltaPhiAll30, &b_minDeltaPhiAll30);
   fChain->SetBranchAddress("minDeltaPhi30_eta5_noIdAll", &minDeltaPhi30_eta5_noIdAll, &b_minDeltaPhi30_eta5_noIdAll);
   fChain->SetBranchAddress("minDeltaPhiMetTau", &minDeltaPhiMetTau, &b_minDeltaPhiMetTau);
   fChain->SetBranchAddress("maxDeltaPhi", &maxDeltaPhi, &b_maxDeltaPhi);
   fChain->SetBranchAddress("maxDeltaPhiAll", &maxDeltaPhiAll, &b_maxDeltaPhiAll);
   fChain->SetBranchAddress("maxDeltaPhiAll30", &maxDeltaPhiAll30, &b_maxDeltaPhiAll30);
   fChain->SetBranchAddress("maxDeltaPhi30_eta5_noIdAll", &maxDeltaPhi30_eta5_noIdAll, &b_maxDeltaPhi30_eta5_noIdAll);
   fChain->SetBranchAddress("sumDeltaPhi", &sumDeltaPhi, &b_sumDeltaPhi);
   fChain->SetBranchAddress("diffDeltaPhi", &diffDeltaPhi, &b_diffDeltaPhi);
   fChain->SetBranchAddress("minDeltaPhiN", &minDeltaPhiN, &b_minDeltaPhiN);
   fChain->SetBranchAddress("deltaPhiN1", &deltaPhiN1, &b_deltaPhiN1);
   fChain->SetBranchAddress("deltaPhiN2", &deltaPhiN2, &b_deltaPhiN2);
   fChain->SetBranchAddress("deltaPhiN3", &deltaPhiN3, &b_deltaPhiN3);
   fChain->SetBranchAddress("minDeltaPhiN_otherEta5", &minDeltaPhiN_otherEta5, &b_minDeltaPhiN_otherEta5);
   fChain->SetBranchAddress("minDeltaPhiN_otherEta5idNo", &minDeltaPhiN_otherEta5idNo, &b_minDeltaPhiN_otherEta5idNo);
   fChain->SetBranchAddress("minDeltaPhiN_mainPt30_otherEta5idNo", &minDeltaPhiN_mainPt30_otherEta5idNo, &b_minDeltaPhiN_mainPt30_otherEta5idNo);
   fChain->SetBranchAddress("minDeltaPhiN_mainPt30Eta5_otherEta5idNo", &minDeltaPhiN_mainPt30Eta5_otherEta5idNo, &b_minDeltaPhiN_mainPt30Eta5_otherEta5idNo);
   fChain->SetBranchAddress("minDeltaPhiN_mainPt30Eta5idNo_otherEta5idNo", &minDeltaPhiN_mainPt30Eta5idNo_otherEta5idNo, &b_minDeltaPhiN_mainPt30Eta5idNo_otherEta5idNo);
   fChain->SetBranchAddress("minDeltaPhiK", &minDeltaPhiK, &b_minDeltaPhiK);
   fChain->SetBranchAddress("minDeltaPhiK_otherEta5", &minDeltaPhiK_otherEta5, &b_minDeltaPhiK_otherEta5);
   fChain->SetBranchAddress("minDeltaPhiK_otherEta5idNo", &minDeltaPhiK_otherEta5idNo, &b_minDeltaPhiK_otherEta5idNo);
   fChain->SetBranchAddress("minDeltaPhiK_mainPt30_otherEta5idNo", &minDeltaPhiK_mainPt30_otherEta5idNo, &b_minDeltaPhiK_mainPt30_otherEta5idNo);
   fChain->SetBranchAddress("minDeltaPhiK_mainPt30Eta5_otherEta5idNo", &minDeltaPhiK_mainPt30Eta5_otherEta5idNo, &b_minDeltaPhiK_mainPt30Eta5_otherEta5idNo);
   fChain->SetBranchAddress("minDeltaPhiK_mainPt30Eta5idNo_otherEta5idNo", &minDeltaPhiK_mainPt30Eta5idNo_otherEta5idNo, &b_minDeltaPhiK_mainPt30Eta5idNo_otherEta5idNo);
   fChain->SetBranchAddress("minDeltaPhiN_DJR", &minDeltaPhiN_DJR, &b_minDeltaPhiN_DJR);
   fChain->SetBranchAddress("minDeltaPhiN_DJR_otherEta5", &minDeltaPhiN_DJR_otherEta5, &b_minDeltaPhiN_DJR_otherEta5);
   fChain->SetBranchAddress("minDeltaPhiK_DJR", &minDeltaPhiK_DJR, &b_minDeltaPhiK_DJR);
   fChain->SetBranchAddress("minDeltaPhiK_DJR_otherEta5", &minDeltaPhiK_DJR_otherEta5, &b_minDeltaPhiK_DJR_otherEta5);
   fChain->SetBranchAddress("minDeltaPhiN_lostJet", &minDeltaPhiN_lostJet, &b_minDeltaPhiN_lostJet);
   fChain->SetBranchAddress("deltaPhiN1_lostJet", &deltaPhiN1_lostJet, &b_deltaPhiN1_lostJet);
   fChain->SetBranchAddress("deltaPhiN2_lostJet", &deltaPhiN2_lostJet, &b_deltaPhiN2_lostJet);
   fChain->SetBranchAddress("deltaPhiN3_lostJet", &deltaPhiN3_lostJet, &b_deltaPhiN3_lostJet);
   fChain->SetBranchAddress("jetpt1", &jetpt1, &b_jetpt1);
   fChain->SetBranchAddress("jeteta1", &jeteta1, &b_jeteta1);
   fChain->SetBranchAddress("jetphi1", &jetphi1, &b_jetphi1);
   fChain->SetBranchAddress("jetenergy1", &jetenergy1, &b_jetenergy1);
   fChain->SetBranchAddress("jetpt2", &jetpt2, &b_jetpt2);
   fChain->SetBranchAddress("jeteta2", &jeteta2, &b_jeteta2);
   fChain->SetBranchAddress("jetphi2", &jetphi2, &b_jetphi2);
   fChain->SetBranchAddress("jetenergy2", &jetenergy2, &b_jetenergy2);
   fChain->SetBranchAddress("jetpt3", &jetpt3, &b_jetpt3);
   fChain->SetBranchAddress("jeteta3", &jeteta3, &b_jeteta3);
   fChain->SetBranchAddress("jetphi3", &jetphi3, &b_jetphi3);
   fChain->SetBranchAddress("jetenergy3", &jetenergy3, &b_jetenergy3);
   fChain->SetBranchAddress("bjetpt1", &bjetpt1, &b_bjetpt1);
   fChain->SetBranchAddress("bjeteta1", &bjeteta1, &b_bjeteta1);
   fChain->SetBranchAddress("bjetphi1", &bjetphi1, &b_bjetphi1);
   fChain->SetBranchAddress("bjetenergy1", &bjetenergy1, &b_bjetenergy1);
   fChain->SetBranchAddress("bjetpt2", &bjetpt2, &b_bjetpt2);
   fChain->SetBranchAddress("bjeteta2", &bjeteta2, &b_bjeteta2);
   fChain->SetBranchAddress("bjetphi2", &bjetphi2, &b_bjetphi2);
   fChain->SetBranchAddress("bjetenergy2", &bjetenergy2, &b_bjetenergy2);
   fChain->SetBranchAddress("bjetpt3", &bjetpt3, &b_bjetpt3);
   fChain->SetBranchAddress("bjeteta3", &bjeteta3, &b_bjeteta3);
   fChain->SetBranchAddress("bjetphi3", &bjetphi3, &b_bjetphi3);
   fChain->SetBranchAddress("bjetenergy3", &bjetenergy3, &b_bjetenergy3);
   fChain->SetBranchAddress("eleet1", &eleet1, &b_eleet1);
   fChain->SetBranchAddress("muonpt1", &muonpt1, &b_muonpt1);
   fChain->SetBranchAddress("lambda1_allJets", &lambda1_allJets, &b_lambda1_allJets);
   fChain->SetBranchAddress("lambda2_allJets", &lambda2_allJets, &b_lambda2_allJets);
   fChain->SetBranchAddress("determinant_allJets", &determinant_allJets, &b_determinant_allJets);
   fChain->SetBranchAddress("lambda1_allJetsPlusMET", &lambda1_allJetsPlusMET, &b_lambda1_allJetsPlusMET);
   fChain->SetBranchAddress("lambda2_allJetsPlusMET", &lambda2_allJetsPlusMET, &b_lambda2_allJetsPlusMET);
   fChain->SetBranchAddress("determinant_allJetsPlusMET", &determinant_allJetsPlusMET, &b_determinant_allJetsPlusMET);
   fChain->SetBranchAddress("lambda1_topThreeJets", &lambda1_topThreeJets, &b_lambda1_topThreeJets);
   fChain->SetBranchAddress("lambda2_topThreeJets", &lambda2_topThreeJets, &b_lambda2_topThreeJets);
   fChain->SetBranchAddress("determinant_topThreeJets", &determinant_topThreeJets, &b_determinant_topThreeJets);
   fChain->SetBranchAddress("lambda1_topThreeJetsPlusMET", &lambda1_topThreeJetsPlusMET, &b_lambda1_topThreeJetsPlusMET);
   fChain->SetBranchAddress("lambda2_topThreeJetsPlusMET", &lambda2_topThreeJetsPlusMET, &b_lambda2_topThreeJetsPlusMET);
   fChain->SetBranchAddress("determinant_topThreeJetsPlusMET", &determinant_topThreeJetsPlusMET, &b_determinant_topThreeJetsPlusMET);
   fChain->SetBranchAddress("transverseThrust", &transverseThrust, &b_transverseThrust);
   fChain->SetBranchAddress("transverseThrustPhi", &transverseThrustPhi, &b_transverseThrustPhi);
   fChain->SetBranchAddress("transverseThrustWithMET", &transverseThrustWithMET, &b_transverseThrustWithMET);
   fChain->SetBranchAddress("transverseThrustWithMETPhi", &transverseThrustWithMETPhi, &b_transverseThrustWithMETPhi);
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

Bool_t myCutflow::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void myCutflow::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t myCutflow::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef myCutflow_cxx
