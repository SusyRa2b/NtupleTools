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
   ULong64_t       runNumber;
   ULong64_t       lumiSection;
   ULong64_t       eventNumber;
   Bool_t          cutHT;
   Bool_t          cutPV;
   Bool_t          cutTrigger;
   Bool_t          cut3Jets;
   Bool_t          cutEleVeto;
   Bool_t          cutMuVeto;
   Bool_t          cutMET;
   Bool_t          cutDeltaPhi;
   Bool_t          cutCleaning;
   Bool_t          passBadPFMuon;
   Bool_t          passInconsistentMuon;
   Int_t           njets;
   Int_t           nbjets;
   Int_t           nElectrons;
   Int_t           nMuons;
   Float_t         HT;
   Float_t         MET;
   Float_t         METphi;
   Float_t         MHT;
   Float_t         bestWMass;
   Float_t         bestTopMass;
   Float_t         topCosHel;
   Float_t         WCosHel;
   Float_t         MT_Wlep;
   Float_t         minDeltaPhi;
   Float_t         minDeltaPhiAll;
   Float_t         minDeltaPhiAll30;
   Float_t         minDeltaPhi30_eta5_noIdAll;
   Float_t         maxDeltaPhi;
   Float_t         maxDeltaPhiAll;
   Float_t         maxDeltaPhiAll30;
   Float_t         maxDeltaPhi30_eta5_noIdAll;
   Float_t         sumDeltaPhi;
   Float_t         diffDeltaPhi;
   Float_t         minDeltaPhiN;
   Float_t         jetpt1;
   Float_t         jeteta1;
   Float_t         jetphi1;
   Float_t         jetpt2;
   Float_t         jeteta2;
   Float_t         jetphi2;
   Float_t         jetpt3;
   Float_t         jeteta3;
   Float_t         jetphi3;
   Float_t         bjetpt1;
   Float_t         bjeteta1;
   Float_t         bjetphi1;
   Float_t         bjetpt2;
   Float_t         bjeteta2;
   Float_t         bjetphi2;
   Float_t         bjetpt3;
   Float_t         bjeteta3;
   Float_t         bjetphi3;
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

   // List of branches
   TBranch        *b_weight;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_lumiSection;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_cutHT;   //!
   TBranch        *b_cutPV;   //!
   TBranch        *b_cutTrigger;   //!
   TBranch        *b_cut3Jets;   //!
   TBranch        *b_cutEleVeto;   //!
   TBranch        *b_cutMuVeto;   //!
   TBranch        *b_cutMET;   //!
   TBranch        *b_cutDeltaPhi;   //!
   TBranch        *b_cutCleaning;   //!
   TBranch        *b_passBadPFMuon;   //!
   TBranch        *b_passInconsistentMuon;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_nbjets;   //!
   TBranch        *b_nElectrons;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METphi;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_bestWMass;   //!
   TBranch        *b_bestTopMass;   //!
   TBranch        *b_topCosHel;   //!
   TBranch        *b_WCosHel;   //!
   TBranch        *b_MT_Wlep;   //!
   TBranch        *b_minDeltaPhi;   //!
   TBranch        *b_minDeltaPhiAll;   //!
   TBranch        *b_minDeltaPhiAll30;   //!
   TBranch        *b_minDeltaPhi30_eta5_noIdAll;   //!
   TBranch        *b_maxDeltaPhi;   //!
   TBranch        *b_maxDeltaPhiAll;   //!
   TBranch        *b_maxDeltaPhiAll30;   //!
   TBranch        *b_maxDeltaPhi30_eta5_noIdAll;   //!
   TBranch        *b_sumDeltaPhi;   //!
   TBranch        *b_diffDeltaPhi;   //!
   TBranch        *b_minDeltaPhiN;   //!
   TBranch        *b_jetpt1;   //!
   TBranch        *b_jeteta1;   //!
   TBranch        *b_jetphi1;   //!
   TBranch        *b_jetpt2;   //!
   TBranch        *b_jeteta2;   //!
   TBranch        *b_jetphi2;   //!
   TBranch        *b_jetpt3;   //!
   TBranch        *b_jeteta3;   //!
   TBranch        *b_jetphi3;   //!
   TBranch        *b_bjetpt1;   //!
   TBranch        *b_bjeteta1;   //!
   TBranch        *b_bjetphi1;   //!
   TBranch        *b_bjetpt2;   //!
   TBranch        *b_bjeteta2;   //!
   TBranch        *b_bjetphi2;   //!
   TBranch        *b_bjetpt3;   //!
   TBranch        *b_bjeteta3;   //!
   TBranch        *b_bjetphi3;   //!
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
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-05/reducedTree.SSVHPT.TTbarJets.root/reducedTree");
      }

      else if (sampleName=="SingleTop"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-05/reducedTree.SSVHPT.SingleTop.root/reducedTree");
      }

      else if (sampleName=="WJets"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-05/reducedTree.SSVHPT.WJets.root/reducedTree");
      }

      else if (sampleName=="ZJets"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-05/reducedTree.SSVHPT.ZJets.root/reducedTree");
      }

      else if (sampleName=="Zinvisible"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-05/reducedTree.SSVHPT.Zinvisible.root/reducedTree");
      }

      else if (sampleName=="PythiaPUQCD"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-05/reducedTree.SSVHPT.PythiaPUQCD.root/reducedTree");
      }

      else if (sampleName=="LM9"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-05/reducedTree.SSVHPT.LM9.root/reducedTree");
      }

      else if (sampleName=="LM13"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-05/reducedTree.SSVHPT.LM13.root/reducedTree");
      }

      else if (sampleName=="Data"){
	chain->Add("/cu2/ra2b/reducedTrees/V00-02-05/reducedTree.SSVHPT.ht_run2011a_uptojul6.root/reducedTree");
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
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("cutHT", &cutHT, &b_cutHT);
   fChain->SetBranchAddress("cutPV", &cutPV, &b_cutPV);
   fChain->SetBranchAddress("cutTrigger", &cutTrigger, &b_cutTrigger);
   fChain->SetBranchAddress("cut3Jets", &cut3Jets, &b_cut3Jets);
   fChain->SetBranchAddress("cutEleVeto", &cutEleVeto, &b_cutEleVeto);
   fChain->SetBranchAddress("cutMuVeto", &cutMuVeto, &b_cutMuVeto);
   fChain->SetBranchAddress("cutMET", &cutMET, &b_cutMET);
   fChain->SetBranchAddress("cutDeltaPhi", &cutDeltaPhi, &b_cutDeltaPhi);
   fChain->SetBranchAddress("cutCleaning", &cutCleaning, &b_cutCleaning);
   fChain->SetBranchAddress("passBadPFMuon", &passBadPFMuon, &b_passBadPFMuon);
   fChain->SetBranchAddress("passInconsistentMuon", &passInconsistentMuon, &b_passInconsistentMuon);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("nbjets", &nbjets, &b_nbjets);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("METphi", &METphi, &b_METphi);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("bestWMass", &bestWMass, &b_bestWMass);
   fChain->SetBranchAddress("bestTopMass", &bestTopMass, &b_bestTopMass);
   fChain->SetBranchAddress("topCosHel", &topCosHel, &b_topCosHel);
   fChain->SetBranchAddress("WCosHel", &WCosHel, &b_WCosHel);
   fChain->SetBranchAddress("MT_Wlep", &MT_Wlep, &b_MT_Wlep);
   fChain->SetBranchAddress("minDeltaPhi", &minDeltaPhi, &b_minDeltaPhi);
   fChain->SetBranchAddress("minDeltaPhiAll", &minDeltaPhiAll, &b_minDeltaPhiAll);
   fChain->SetBranchAddress("minDeltaPhiAll30", &minDeltaPhiAll30, &b_minDeltaPhiAll30);
   fChain->SetBranchAddress("minDeltaPhi30_eta5_noIdAll", &minDeltaPhi30_eta5_noIdAll, &b_minDeltaPhi30_eta5_noIdAll);
   fChain->SetBranchAddress("maxDeltaPhi", &maxDeltaPhi, &b_maxDeltaPhi);
   fChain->SetBranchAddress("maxDeltaPhiAll", &maxDeltaPhiAll, &b_maxDeltaPhiAll);
   fChain->SetBranchAddress("maxDeltaPhiAll30", &maxDeltaPhiAll30, &b_maxDeltaPhiAll30);
   fChain->SetBranchAddress("maxDeltaPhi30_eta5_noIdAll", &maxDeltaPhi30_eta5_noIdAll, &b_maxDeltaPhi30_eta5_noIdAll);
   fChain->SetBranchAddress("sumDeltaPhi", &sumDeltaPhi, &b_sumDeltaPhi);
   fChain->SetBranchAddress("diffDeltaPhi", &diffDeltaPhi, &b_diffDeltaPhi);
   fChain->SetBranchAddress("minDeltaPhiN", &minDeltaPhiN, &b_minDeltaPhiN);
   fChain->SetBranchAddress("jetpt1", &jetpt1, &b_jetpt1);
   fChain->SetBranchAddress("jeteta1", &jeteta1, &b_jeteta1);
   fChain->SetBranchAddress("jetphi1", &jetphi1, &b_jetphi1);
   fChain->SetBranchAddress("jetpt2", &jetpt2, &b_jetpt2);
   fChain->SetBranchAddress("jeteta2", &jeteta2, &b_jeteta2);
   fChain->SetBranchAddress("jetphi2", &jetphi2, &b_jetphi2);
   fChain->SetBranchAddress("jetpt3", &jetpt3, &b_jetpt3);
   fChain->SetBranchAddress("jeteta3", &jeteta3, &b_jeteta3);
   fChain->SetBranchAddress("jetphi3", &jetphi3, &b_jetphi3);
   fChain->SetBranchAddress("bjetpt1", &bjetpt1, &b_bjetpt1);
   fChain->SetBranchAddress("bjeteta1", &bjeteta1, &b_bjeteta1);
   fChain->SetBranchAddress("bjetphi1", &bjetphi1, &b_bjetphi1);
   fChain->SetBranchAddress("bjetpt2", &bjetpt2, &b_bjetpt2);
   fChain->SetBranchAddress("bjeteta2", &bjeteta2, &b_bjeteta2);
   fChain->SetBranchAddress("bjetphi2", &bjetphi2, &b_bjetphi2);
   fChain->SetBranchAddress("bjetpt3", &bjetpt3, &b_bjetpt3);
   fChain->SetBranchAddress("bjeteta3", &bjeteta3, &b_bjeteta3);
   fChain->SetBranchAddress("bjetphi3", &bjetphi3, &b_bjetphi3);
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
