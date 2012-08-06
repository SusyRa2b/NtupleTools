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

//root -l tagAndProbeSkim.C++

void tagAndProbeSkim(){

  bool datamode = true; //set to false if running over MC

  TChain *cData = new TChain("reducedTree");

  if(datamode){

    cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011a_may10rereco.root");
    cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011a_aug5rereco.root");
    cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011b_promptrecov1_nov4_try2.root");
    cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011b_promptrecov1_oct28_plusLostJobs.root");
    cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011a_promptrecov4_try3.root");
    cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011a_promptrecov6.root");
    cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011b_promptrecov1_oct7.root");
    cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011b_promptrecov1.root");
    cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011b_promptrecov1_oct21.root");
    cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/OnlySTSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ht_run2011b_promptrecov1_oct14.root");

    //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singlemu_run2011b_promptrecov1_part1.root");
    //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singlemu_run2011b_promptrecov1_part2.root");
    //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singlemu_run2011b_promptrecov1_part3.root");
    //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singlemu_run2011b_promptrecov1_part4.root");
    //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singlemu_run2011b_promptrecov1_part5.root");
    //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singlemu_run2011b_promptrecov1_part6.root");

    //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singleelectron_run2011b_promptrecov1_part1.root");
    //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singleelectron_run2011b_promptrecov1_part2.root");
    //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/v3/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singleelectron_run2011b_promptrecov1_part1.root");
    //cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/NoSkim/v3/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.singleelectron_run2011b_promptrecov1_part2.root");
  }
  else{
    cData->Add("/cu3/wteo/reducedTrees/V00-02-35y_ZllMET_v2/Fall11/reducedTree.CSVM_PF2PATjets_JES0_JERbias_PFMET_METunc0_PUunc0_BTagEff04_HLTEff0.ZJets.root");
  }

  //cData->SetBranchStatus("*",0);
  //cData->SetBranchStatus("HT",1);


  Long64_t nentries = cData->GetEntries();

  //these must be floats (i.e. their original data type)! (making them doubles will return nonsense)
  float muprobe1_mass,muprobe2_mass,muprobe3_mass;
  float eleprobe1_mass,eleprobe2_mass,eleprobe3_mass;

  bool muprobe1_result,muprobe2_result,muprobe3_result;
  bool eleprobe1_result,eleprobe2_result,eleprobe3_result;
  bool cut3Jets,cutPV,cutHT;
  int njets;
  float MET, HT;
  cData->SetBranchAddress("muprobe1_mass",&muprobe1_mass);
  cData->SetBranchAddress("muprobe2_mass",&muprobe2_mass);
  cData->SetBranchAddress("muprobe3_mass",&muprobe3_mass);

  cData->SetBranchAddress("eleprobe1_mass",&eleprobe1_mass);
  cData->SetBranchAddress("eleprobe2_mass",&eleprobe2_mass);
  cData->SetBranchAddress("eleprobe3_mass",&eleprobe3_mass);
  //cData->SetBranchAddress("eleprobeID1_mass",&eleprobe1_mass);
  //cData->SetBranchAddress("eleprobeID2_mass",&eleprobe2_mass);
  //cData->SetBranchAddress("eleprobeID3_mass",&eleprobe3_mass);

  cData->SetBranchAddress("muprobe1_result",&muprobe1_result);
  cData->SetBranchAddress("muprobe2_result",&muprobe2_result);
  cData->SetBranchAddress("muprobe3_result",&muprobe3_result);

  cData->SetBranchAddress("eleprobe1_result",&eleprobe1_result);
  cData->SetBranchAddress("eleprobe2_result",&eleprobe2_result);
  cData->SetBranchAddress("eleprobe3_result",&eleprobe3_result);
  //cData->SetBranchAddress("eleprobeID1_result",&eleprobe1_result);
  //cData->SetBranchAddress("eleprobeID2_result",&eleprobe2_result);
  //cData->SetBranchAddress("eleprobeID3_result",&eleprobe3_result);

  cData->SetBranchAddress("cut3Jets",&cut3Jets);
  cData->SetBranchAddress("cutHT",&cutHT);
  cData->SetBranchAddress("cutPV",&cutPV);
  cData->SetBranchAddress("HT",&HT);
  cData->SetBranchAddress("njets",&njets);
  //cData->SetBranchAddress("cutMET",&cutMET);
  //cData->SetBranchAddress("MET",&MET);


  int npass=0;

  //these must be floats!  (making them doubles will return nonsense)
  float muprobe_mass_pass;
  float muprobe_mass_fail;
  float eleprobe_mass_pass;
  float eleprobe_mass_fail;

  //Create a new file + a clone of old tree in new file
  //TFile *newfile = new TFile("tagandprobe_slim_noHT_geq2jetcut.root","recreate");
  //TFile *newfile = new TFile("tagandprobe_slim_withHT_with3jetcut_ZJetsMC.root","recreate");
  //TFile *newfile = new TFile("tagandprobe_slim_noHT_geq1jetcut_singlemu.root","recreate");
  //TFile *newfile = new TFile("tagandprobe_slim_noHT_geq2jetcut_singleelectron.root","recreate");
  TFile *newfile = new TFile("tagandprobe_slim_noHT_geq3jetcut_IDbase_singleelectron.root","recreate");
  //TTree *newtree = cData->CloneTree(0);

  TTree tnpTreeMuPass("tnpTreeMuPass","tree for tag and probe");
  //tnpTreeMuPass.Branch("muprobe_mass_pass",&muprobe_mass_pass,"muprobe_mass_pass/F");
  tnpTreeMuPass.Branch("m",&muprobe_mass_pass,"m/F");

  TTree tnpTreeMuFail("tnpTreeMuFail","tree for tag and probe");
  //tnpTreeMuFail.Branch("muprobe_mass_fail",&muprobe_mass_fail,"muprobe_mass_fail/F");
  tnpTreeMuFail.Branch("m",&muprobe_mass_fail,"m/F");
  
  TTree tnpTreeElePass("tnpTreeElePass","tree for tag and probe");
  //tnpTreeElePass.Branch("eleprobe_mass_pass",&eleprobe_mass_pass,"eleprobe_mass_pass/F");
  tnpTreeElePass.Branch("m",&eleprobe_mass_pass,"m/F");
  
  TTree tnpTreeEleFail("tnpTreeEleFail","tree for tag and probe");
  //tnpTreeEleFail.Branch("eleprobe_mass_fail",&eleprobe_mass_fail,"eleprobe_mass_fail/F");
  tnpTreeEleFail.Branch("m",&eleprobe_mass_fail,"m/F");
  

  bool alreadyPrinted = false;
  for (Long64_t i=0;i<nentries; i++) {

    if(!(i%100000)) cout << "processing entry " << i << endl;
    //if(i>300000) break;
    cData->GetEntry(i);
    ////if (event->GetNtrack() > 605) newtree->Fill();
    if(  njets>=3 && cutPV==1 && muprobe1_mass>0 && muprobe1_result==1){ muprobe_mass_pass = muprobe1_mass; 
      //cout << "muprobe1_mass (pass) = " << muprobe1_mass << endl; 
      tnpTreeMuPass.Fill();npass++;alreadyPrinted=false;}
    //if(!(npass%10) && !alreadyPrinted) {cout << npass << endl;alreadyPrinted=true;}
    if(  njets>=3 && cutPV==1 && muprobe2_mass>0 && muprobe2_result==1){ muprobe_mass_pass = muprobe2_mass; tnpTreeMuPass.Fill();npass++;alreadyPrinted=false;}
    //if(!(npass%10) && !alreadyPrinted) {cout << npass << endl;alreadyPrinted=true;}
    if(  njets>=3 && cutPV==1 && muprobe3_mass>0 && muprobe3_result==1){ muprobe_mass_pass = muprobe3_mass; tnpTreeMuPass.Fill();npass++;alreadyPrinted=false;}
    //if(!(npass%10) && !alreadyPrinted) {cout << npass << endl;alreadyPrinted=true;}    

    if(  njets>=3 && cutPV==1 && muprobe1_mass>0 && muprobe1_result==0){ muprobe_mass_fail = muprobe1_mass; 
      //cout << "muprobe1_mass (fail) = " << muprobe1_mass << endl; 
      tnpTreeMuFail.Fill();}
    if(  njets>=3 && cutPV==1 && muprobe2_mass>0 && muprobe2_result==0){ muprobe_mass_fail = muprobe2_mass; tnpTreeMuFail.Fill();}
    if(  njets>=3 && cutPV==1 && muprobe3_mass>0 && muprobe3_result==0){ muprobe_mass_fail = muprobe3_mass; tnpTreeMuFail.Fill();}


    if(  njets>=3 && cutPV==1 && eleprobe1_mass>0 && eleprobe1_result==1){ eleprobe_mass_pass = eleprobe1_mass; tnpTreeElePass.Fill();}
    if(  njets>=3 && cutPV==1 && eleprobe2_mass>0 && eleprobe2_result==1){ eleprobe_mass_pass = eleprobe2_mass; tnpTreeElePass.Fill();}
    if(  njets>=3 && cutPV==1 && eleprobe3_mass>0 && eleprobe3_result==1){ eleprobe_mass_pass = eleprobe3_mass; tnpTreeElePass.Fill();}
    
    if(  njets>=3 && cutPV==1 && eleprobe1_mass>0 && eleprobe1_result==0){ eleprobe_mass_fail = eleprobe1_mass; tnpTreeEleFail.Fill();}
    if(  njets>=3 && cutPV==1 && eleprobe2_mass>0 && eleprobe2_result==0){ eleprobe_mass_fail = eleprobe2_mass; tnpTreeEleFail.Fill();}
    if(  njets>=3 && cutPV==1 && eleprobe3_mass>0 && eleprobe3_result==0){ eleprobe_mass_fail = eleprobe3_mass; tnpTreeEleFail.Fill();}

    //event->Clear();

    //if(npass>100) break;
  }




  //newtree->Print();
  //newtree->AutoSave();
  newfile->Write();

  delete cData;
  delete newfile;


}
