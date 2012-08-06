#include <iostream>
#include "TChain.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

//root -l zllPuritySkim.C++

void zllPuritySkim(){

  bool mumode = false;  //mu or ele?
  bool datamode = true; //set to false if running over MC
  
  TChain *cData = new TChain("reducedTree");

  if(datamode){
    if(mumode){
      //doublemu data
      cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.doublemu_run2011a_aug5rereco.root");
      cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.doublemu_run2011a_may10rereco.root");
      cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.doublemu_run2011a_promptrecov4.root");
      cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.doublemu_run2011a_promptrecov6.root");
      cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.doublemu_run2011b_promptrecov1.root");
    }
    else{
      //doubleele data
      cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.doubleelectron_run2011a_aug5rereco.root");
      cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.doubleelectron_run2011a_may10rereco.root");
      cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.doubleelectron_run2011a_promptrecov4.root");
      cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.doubleelectron_run2011a_promptrecov6.root");
      cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.doubleelectron_run2011b_promptrecov1.root");
    }
  }
  else cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/Fall11/reducedTree.CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ZJets.root");

  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011a_may10rereco.root");
  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011a_aug5rereco.root");
  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011b_promptrecov1_nov4_try2.root");
  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011b_promptrecov1_oct28_plusLostJobs.root");
  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011a_promptrecov4_try3.root");
  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011a_promptrecov6.root");
  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011b_promptrecov1_oct7.root");
  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011b_promptrecov1.root");
  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011b_promptrecov1_oct21.root");
  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011b_promptrecov1_oct14.root");

  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singlemu_run2011b_promptrecov1_part1.root");
  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singlemu_run2011b_promptrecov1_part2.root");
  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singlemu_run2011b_promptrecov1_part3.root");
  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singlemu_run2011b_promptrecov1_part4.root");
  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singlemu_run2011b_promptrecov1_part5.root");
  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singlemu_run2011b_promptrecov1_part6.root");

  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singleelectron_run2011b_promptrecov1_part1.root");
  //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singleelectron_run2011b_promptrecov1_part2.root");


  //cData->SetBranchStatus("*",0);
  //cData->SetBranchStatus("HT",1);


  Long64_t nentries = cData->GetEntries();

  //these must be floats (i.e. their original data type)! (making them doubles will return nonsense)
  float zmumuMass, zmumuMinDeltaPhiN, zmumuMET;
  float zeeMass, zeeMinDeltaPhiN, zeeMET;

  bool pass_ZmumuHLT, passZmumuCandLoose;
  bool pass_ZeeHLT, passZeeCandLoose;
  bool cut3Jets,cutPV, cutEleVeto, cutMuVeto,cutHT;
  int njets,nMuons,nElectrons, nbjetsCSVL;
  float MET, HT;
  cData->SetBranchAddress("zmumuMass",&zmumuMass);
  cData->SetBranchAddress("zeeMass",&zeeMass);
  cData->SetBranchAddress("pass_ZmumuHLT",&pass_ZmumuHLT);
  cData->SetBranchAddress("pass_ZeeHLT",&pass_ZeeHLT);
  cData->SetBranchAddress("passZmumuCandLoose",&passZmumuCandLoose);
  cData->SetBranchAddress("passZeeCandLoose",&passZeeCandLoose);
  cData->SetBranchAddress("cut3Jets",&cut3Jets);
  cData->SetBranchAddress("cutMuVeto",&cutMuVeto);
  cData->SetBranchAddress("cutEleVeto",&cutEleVeto);
  cData->SetBranchAddress("cutPV",&cutPV);
  cData->SetBranchAddress("HT",&HT);
  cData->SetBranchAddress("njets",&njets);
  cData->SetBranchAddress("nMuons",&nMuons);
  cData->SetBranchAddress("nElectrons",&nElectrons);
  cData->SetBranchAddress("nbjetsCSVL",&nbjetsCSVL);

  cData->SetBranchAddress("zmumuMET",&zmumuMET);
  cData->SetBranchAddress("zeeMET",&zeeMET);
  cData->SetBranchAddress("zmumuMinDeltaPhiN",&zmumuMinDeltaPhiN);
  cData->SetBranchAddress("zeeMinDeltaPhiN",&zeeMinDeltaPhiN);
  //cData->SetBranchAddress("cutMET",&cutMET);
  cData->SetBranchAddress("cutHT",&cutHT);
  //cData->SetBranchAddress("MET",&MET);


  int npass=0;

  //these must be floats!  (making them doubles will return nonsense)
  float mass_Loose;
  float mass_VL;
  float mass_VL_SB;

  //Create a new file + a clone of old tree in new file
  TFile *newfile;
  if(mumode) newfile = new TFile("puritytree_zmumu.root","recreate");
  else newfile = new TFile("puritytree_zee.root","recreate");
  //TFile *newfile = new TFile("tagandprobe_slim_noHT_geq2jetcut_singleelectron.root","recreate");
  //TTree *newtree = cData->CloneTree(0);

  TTree purTreeMVL("purTreeMVL","tree for purity VL");
  purTreeMVL.Branch("m",&mass_VL,"m/F");

  TTree purTreeMVL_SB("purTreeMVL_SB","tree for purity VL");
  purTreeMVL_SB.Branch("m",&mass_VL_SB,"m/F");

  TTree purTreeMLoose("purTreeMLoose","tree for purity Loose");
  purTreeMLoose.Branch("m",&mass_Loose,"m/F");
  

  for (Long64_t i=0;i<nentries; i++) {

    if(!(i%100000)) cout << "processing entry " << i << endl;
    //if(i>300000) break;
    cData->GetEntry(i);
    ////if (event->GetNtrack() > 605) newtree->Fill();

    if(mumode){
      //doublemu
      if( pass_ZmumuHLT==1 && cutPV==1 && passZmumuCandLoose ){
	mass_Loose = zmumuMass; 
	purTreeMLoose.Fill();      
      }
      if( pass_ZmumuHLT==1 && cutHT==1 && cutPV==1 && cut3Jets==1 && passZmumuCandLoose && zmumuMinDeltaPhiN>4 && zmumuMET>250 && nMuons<3 && cutEleVeto==1 && nbjetsCSVL>=1){
	mass_VL = zmumuMass;
	purTreeMVL.Fill();
      }
      if( pass_ZmumuHLT==1 && cutHT==1 && cutPV==1 && cut3Jets==1 && passZmumuCandLoose && zmumuMinDeltaPhiN>4 && zmumuMET>150 && zmumuMET<250 && nMuons<3 && cutEleVeto==1 && nbjetsCSVL>=1){
	mass_VL_SB = zmumuMass;
	purTreeMVL_SB.Fill();
      }
    }
    else{
      //doubleele
      if( pass_ZeeHLT==1 && cutPV==1 && passZeeCandLoose ){
	mass_Loose = zeeMass; 
	purTreeMLoose.Fill();
      }
      if( pass_ZeeHLT==1 && cutHT==1 && cutPV==1 && cut3Jets==1 && passZeeCandLoose && zeeMinDeltaPhiN>4 && zeeMET>250 && nElectrons<3 && cutMuVeto==1 && nbjetsCSVL>=1){
	mass_VL = zeeMass;
	purTreeMVL.Fill();
      }
    
      if( pass_ZeeHLT==1 && cutHT==1 && cutPV==1 && cut3Jets==1 && passZeeCandLoose && zeeMinDeltaPhiN>4 && zeeMET>150 && zeeMET<250 && nElectrons<3 && cutMuVeto==1 && nbjetsCSVL>=1){
	mass_VL_SB = zeeMass;
	purTreeMVL_SB.Fill();
      }
    }

    //event->Clear();
    //if(npass>100) break;
  }


  //newtree->Print();
  //newtree->AutoSave();
  newfile->Write();

  delete cData;
  delete newfile;

}
