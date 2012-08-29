#define ISRweights_cxx
#include "ISRweights.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>

using namespace std;


void ISRweights::cleanup3B() {
  //complete manual!
  TFile f("ISRoutput.3B.T1bbbb.root","update");
  TH2D* isrOverNominalRatioBinary = (TH2D*) f.Get("isrOverNominalRatioBinary");
  TH2D* isrOverNominalRatioBinaryClean = (TH2D*) isrOverNominalRatioBinary->Clone("isrOverNominalRatioBinaryClean");

  isrOverNominalRatioBinaryClean->SetBinContent(19,16,0);
  isrOverNominalRatioBinaryClean->SetBinContent(19,17,0);
  isrOverNominalRatioBinaryClean->SetBinContent(21,18,0);
  isrOverNominalRatioBinaryClean->Write();

  f.Close();

}

void ISRweights::cleanup2BT() {
  //complete manual!
  TFile f("ISRoutput.2BT.T1bbbb.root","update");
  TH2D* isrOverNominalRatioBinary = (TH2D*) f.Get("isrOverNominalRatioBinary");
  TH2D* isrOverNominalRatioBinaryClean = (TH2D*) isrOverNominalRatioBinary->Clone("isrOverNominalRatioBinaryClean");

  isrOverNominalRatioBinaryClean->SetBinContent(14,4,0);
  isrOverNominalRatioBinaryClean->SetBinContent(15,3,0);
  isrOverNominalRatioBinaryClean->SetBinContent(16,3,0);
  isrOverNominalRatioBinaryClean->SetBinContent(16,4,0);
  isrOverNominalRatioBinaryClean->SetBinContent(16,5,0);

  isrOverNominalRatioBinaryClean->SetBinContent(17,3,0);
  isrOverNominalRatioBinaryClean->SetBinContent(17,4,0);
  isrOverNominalRatioBinaryClean->SetBinContent(17,5,0);
  isrOverNominalRatioBinaryClean->SetBinContent(17,6,0);
  isrOverNominalRatioBinaryClean->SetBinContent(17,7,0);

  isrOverNominalRatioBinaryClean->SetBinContent(24,21,0);
  isrOverNominalRatioBinaryClean->SetBinContent(24,19,0);

  isrOverNominalRatioBinaryClean->SetBinContent(27,25,0);

  isrOverNominalRatioBinaryClean->Write();

  f.Close();

}

void ISRweights::cleanup1BT() {
  //complete manual!
  TFile f("ISRoutput.1BT.T1bbbb.root","update");
  TH2D* isrOverNominalRatioBinary = (TH2D*) f.Get("isrOverNominalRatioBinary");
  TH2D* isrOverNominalRatioBinaryClean = (TH2D*) isrOverNominalRatioBinary->Clone("isrOverNominalRatioBinaryClean");

  isrOverNominalRatioBinaryClean->SetBinContent(14,4,0);
  isrOverNominalRatioBinaryClean->SetBinContent(16,5,0);
  isrOverNominalRatioBinaryClean->SetBinContent(17,3,0);

  isrOverNominalRatioBinaryClean->SetBinContent(19,3,0);
  isrOverNominalRatioBinaryClean->SetBinContent(20,3,0);

  isrOverNominalRatioBinaryClean->SetBinContent(19,6,0);
  isrOverNominalRatioBinaryClean->SetBinContent(20,6,0);

  isrOverNominalRatioBinaryClean->SetBinContent(24,13,0);

  isrOverNominalRatioBinaryClean->SetBinContent(32,27,0);
  isrOverNominalRatioBinaryClean->SetBinContent(32,26,0);

  isrOverNominalRatioBinaryClean->SetBinContent(33,26,0);

  isrOverNominalRatioBinaryClean->SetBinContent(33,31,0);

  isrOverNominalRatioBinaryClean->SetBinContent(34,32,0);

  isrOverNominalRatioBinaryClean->SetBinContent(36,32,1);
  isrOverNominalRatioBinaryClean->SetBinContent(35,28,1);
  isrOverNominalRatioBinaryClean->SetBinContent(33,23,1);


  isrOverNominalRatioBinaryClean->Write();

  f.Close();

}
void ISRweights::loadIsrHistos(TString sample) {
  cout<<"Loading Seema's histograms"<<endl;

  if (sample.BeginsWith("T1"))
    fISRweights_ = new TFile("/afs/cern.ch/user/s/seema/public/ISRWeights_14Sep2011/ISRWeights_TopologyT1.root");
  else  if (sample.BeginsWith("T2"))
    fISRweights_ = new TFile("/afs/cern.ch/user/s/seema/public/ISRWeights_14Sep2011/ISRWeights_TopologyT2.root");
  else assert(0);
  int nloaded=0;
  
  if (fISRweights_==0 || fISRweights_->IsZombie()) assert(0);

  const   int ntotal = fISRweights_->GetListOfKeys()->GetEntries();
  for (int ii=0; ii<ntotal; ii++) {
    TString objname=  fISRweights_->GetListOfKeys()->At(ii)->GetName();
    if (objname.BeginsWith("h_ISRWeight") ) {
      TString m0string = objname.Tokenize("_")->At(3)->GetName();
      TString m12string = objname.Tokenize("_")->At(4)->GetName();

      int mym0=m0string.Atoi();
      int mym12=m12string.Atoi();
      isrHistos_[make_pair(mym0,mym12)] = (TH1F*) fISRweights_->Get(objname);
      nloaded++;
    }
  }
  cout<<"Loaded "<<nloaded<<" histograms!"<<endl;

}

double ISRweights::getIsrWeight(int mym0,int mym12, double isrpt) {
  double thisweight=1;
  if ( isrHistos_.count(make_pair(mym0,mym12))>0 ) {
    int bin=  isrHistos_[make_pair(mym0,mym12)]->FindBin(isrpt);
    thisweight= isrHistos_[make_pair(mym0,mym12)]->GetBinContent(bin);
  }
  else {
    cout<<"did not find a histogram for point "<<mym0<<" "<<mym12<<endl;
  }
  return thisweight;
}

void ISRweights::Loop(TString selection)
{
//   In a ROOT session, you can do:
//      Root > .L ISRweights.C++
//      Root > ISRweights t("T1bbbb") (or T1tttt)
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  int    nbinsx=60;
  int    nbinsy=60;
  double    lowx=0;
  double    lowy=0;
  double    highx=1500;
  double    highy=1500;
  TString hname = "efficiency";
  TH2D efficiency(hname,hname,nbinsx,lowx,highx,nbinsy,lowy,highy);
  hname = "isrWeightedEfficiency";
  TH2D isrWeightedEfficiency(hname,hname,nbinsx,lowx,highx,nbinsy,lowy,highy);

  float htcut=400;
  float metcut = 250;

  if (selection=="1BT") {
    htcut = 500;
    metcut = 500;
  }
  else if (selection=="2BT") {
    htcut=600;
    metcut=300;
  }
  else if (selection=="SpecialLoose") {
    htcut=400;
    metcut=150;
  }

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries(); //no FAST!

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     if (jentry%100000 ==0) cout<<100*double(jentry)/double(nentries)<<" %"<<endl;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //apply cuts
      if (passCleaning==1 && njets>=3 && minDeltaPhiN>=4 && nElectrons==0 && nMuons==0) { //baseline
	if (HT>=htcut && MET>=metcut) { //kinematics

	  // ge3b
	  double thisweight = hltHTeff*PUweight;
	  if (selection=="1BL" || selection=="1BT" || selection=="SpecialLoose") thisweight *= probge1;
	  else if (selection=="2BL" || selection=="2BT") thisweight *= probge2;
	  else if (selection=="3B") thisweight *= probge3;
	  else assert(0);
	  efficiency.Fill(m0,m12,thisweight);
	  double isrWeight = getIsrWeight(m0,m12,SUSY_recoilPt);
	  //	  cout<<m0<<" "<<m12<<" "<<SUSY_recoilPt<<" "<<isrWeight<<endl;
	  isrWeightedEfficiency.Fill(m0,m12,thisweight*isrWeight);
	}
      }

   }
   cout<< " ====== "<<endl;
   //now divide by ngen
   for (int ix=1; ix<=nbinsx; ix++) {
     for (int iy=1; iy<=nbinsy; iy++) {
       float x= efficiency.GetXaxis()->GetBinLowEdge(ix);
       float y= efficiency.GetYaxis()->GetBinLowEdge(iy);
       int ngenbin= scanSMSngen_->FindBin(x,y);
       float ngen = scanSMSngen_->GetBinContent(ngenbin);
       //       cout<<x<<" "<<y<<" "<<ngen<<endl;
       if (ngen>0) {
	 efficiency.SetBinContent(ix,iy, efficiency.GetBinContent(ix,iy)/ngen);
	 isrWeightedEfficiency.SetBinContent(ix,iy, isrWeightedEfficiency.GetBinContent(ix,iy)/ngen);
       }
     }
   }

   TH2D* isrOverNominalRatio=(TH2D*) efficiency.Clone("isrOverNominalRatio");
   isrOverNominalRatio->Reset();
   isrOverNominalRatio->Divide(&isrWeightedEfficiency,&efficiency);

   //now prepare a histogram of 0 and 1
   TH2D* isrOverNominalRatioBinary=(TH2D*) isrOverNominalRatio->Clone("isrOverNominalRatioBinary");
   for (int ix=1; ix<=nbinsx; ix++) {
     for (int iy=1; iy<=nbinsy; iy++) {
       float ratioval=    isrOverNominalRatioBinary->GetBinContent(ix,iy);
       if (ratioval > 0.5 ) isrOverNominalRatioBinary->SetBinContent(ix,iy, 1);
       else                 isrOverNominalRatioBinary->SetBinContent(ix,iy, 0);
     }
   }

   TString outfilename="ISRoutput.";
   outfilename+=selection;
   outfilename+=".";
   outfilename+= sample_;
   outfilename+=".root";

   TFile fout(outfilename,"RECREATE");
   efficiency.Write();
   isrWeightedEfficiency.Write();
   isrOverNominalRatio->Write();
   isrOverNominalRatioBinary->Write();
   fout.Close();

}

