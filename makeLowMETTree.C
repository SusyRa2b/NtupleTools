#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

#include <iostream>

using namespace std;

void makeLowMETTree(){

  TFile fout("/cu2/kreis/reducedTrees/V00-02-05/lowMETTree.ht_run2011a_promptrecov4_throughJun24_v2.root","RECREATE");

  TChain* myChain = new TChain("reducedTree");
  myChain->Add("/cu2/kreis/reducedTrees/V00-02-05/reducedTree.SSVHPT.ht_run2011a_promptrecov4_throughJune24.root");

  bool cutHT, cutPV, cut3Jets, cutEleVeto, cutMuVeto, pass_utilityHLT_HT300;
  float HT, MET, minDeltaPhiN, deltaPhiN1, deltaPhiN2, deltaPhiN3;
  int nbjetsSSVHPT;
  double weight;

  myChain->SetBranchAddress("cutHT",&cutHT);
  myChain->SetBranchAddress("cutPV",&cutPV);
  myChain->SetBranchAddress("cut3Jets",&cut3Jets);
  myChain->SetBranchAddress("cutEleVeto",&cutEleVeto);
  myChain->SetBranchAddress("cutMuVeto",&cutMuVeto);
  myChain->SetBranchAddress("pass_utilityHLT_HT300",&pass_utilityHLT_HT300);
  myChain->SetBranchAddress("HT",&HT);
  myChain->SetBranchAddress("MET",&MET);
  myChain->SetBranchAddress("minDeltaPhiN",&minDeltaPhiN);
  myChain->SetBranchAddress("deltaPhiN1",&deltaPhiN1);
  myChain->SetBranchAddress("deltaPhiN2",&deltaPhiN2);
  myChain->SetBranchAddress("deltaPhiN3",&deltaPhiN3);
  myChain->SetBranchAddress("nbjetsSSVHPT",&nbjetsSSVHPT);
  myChain->SetBranchAddress("weight",&weight);


  TTree lowMETTree("lowMETTree","low MET data from prescaled triggers");
  lowMETTree.Branch("HT",&HT,"HT/F");
  lowMETTree.Branch("MET",&MET,"MET/F");
  lowMETTree.Branch("minDeltaPhiN",&minDeltaPhiN,"minDeltaPhiN/F");
  lowMETTree.Branch("deltaPhiN1",&deltaPhiN1,"deltaPhiN1/F");
  lowMETTree.Branch("deltaPhiN2",&deltaPhiN2,"deltaPhiN2/F");
  lowMETTree.Branch("deltaPhiN3",&deltaPhiN3,"deltaPhiN3/F");
  lowMETTree.Branch("nbjetsSSVHPT",&nbjetsSSVHPT,"nbjetsSSVHPT/I"); 
  lowMETTree.Branch("weight",&weight,"weight/D");


  const int numE = myChain->GetEntries();
  cout << "numEntries: " << numE << endl;
  for(int i=0; i<numE; i++){
    if(i%100000==0) cout << "  entry: " << i << ", percent done=" << (int)(i/(double)numE*100.) << endl;

    myChain->GetEvent(i);
    if(! (cutHT==1 && cutPV==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && pass_utilityHLT_HT300==1 ) ) continue;
    lowMETTree.Fill();
  }
  
  fout.Write();
  fout.Close();
}
