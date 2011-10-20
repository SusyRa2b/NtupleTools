#define myCutflow_cxx
#include "myCutflow.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>

//void myCutflow::Loop(TString sampleName, bool isTightSelection, bool writeFiles)
void myCutflow::Loop(bool isTightSelection, bool writeFiles)//:Event_(sampleName)
{
//   In a ROOT session, you can do:
//      Root > .L myCutflow.C
//      Root > myCutflow t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
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
   if (fChain == 0) return;

   const float minHT = isTightSelection ? 500 : 400; 
   const float minMET = isTightSelection ? 300 : 250;
   cout<<"HT >= "<<minHT<<", MET >= "<<minMET<<endl;

   //some counters for the different stages of cuts
   int nPassTrigger = 0;
   int nPassPV = 0;
   int nPassHT = 0;
   int nPass3jets = 0;
   int nPassEleVeto = 0;
   int nPassMuVeto = 0;
   int nPassMET = 0;
   int nPassDPhi = 0;
//   int nPassCleaning = 0;
   int nEq1b = 0;
   int nGe1b = 0;
   int nGe2b = 0;
   int nGe3b = 0;

   vector<double> highestMETevent; 
   vector<double> highestHTevent; 

   vector<TString> textfilenames = getFileNames(Event_,isTightSelection);
   vector<ofstream*> textfiles = getFiles(Event_,isTightSelection);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


      int istage=0;

      /*______________________________________________________ trigger ______________________________________*/

      if (cutTrigger){
	nPassTrigger++;
	if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
	istage++;

      /*______________________________________________________ PV ___________________________________________*/
	
	if (cutPV){
	  nPassPV++;
	  if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
          istage++;

      /*______________________________________________________ HT ____________________________________________*/

//	  if (cutHT){
	  if (HT>=minHT){
	    nPassHT++;
	    if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
	    istage++;

      /*______________________________________________________ >= 3 jets______________________________________*/

	    if (cut3Jets){
	      nPass3jets++;
	      if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
	      istage++;

      /*______________________________________________________ Ele veto ______________________________________*/

	      if (cutEleVeto){
		nPassEleVeto++;
		if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
		istage++;

      /*______________________________________________________ Mu veto _______________________________________*/

		if (cutMuVeto){
		  nPassMuVeto++;
		  if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
		  istage++;

      /*______________________________________________________ MET ___________________________________________*/

//		  if (cutMET){
		  if (MET>=minMET){
		    nPassMET++;
		    if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
		    istage++;

      /*______________________________________________________ min dPhi / min dPhiN __________________________*/

//		    if (cutDeltaPhi){
		    if (minDeltaPhiN>=4){
		      nPassDPhi++;
		      if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
		      istage++;

      /*______________________________________________________ cleaning ______________________________________*/

//		  if (cutCleaning){
//		    nPassCleaning++;
//		    if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
//		    istage++;

      /*______________________________________________________ nB ____________________________________________*/

		      if (nbjetsCSVM>=1) {
			nGe1b++;
			if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
			istage++;
			if (nbjetsCSVM==1) {
			  nEq1b++;
			  if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
			}
			istage++;
			if (nbjetsCSVM>=2){
			  nGe2b++;
			  if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
			  istage++;
			  if (nbjetsCSVM>=3) {
			    nGe3b++;
			    if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
			    istage++;
			  }
			}
		      }//n b jets
//		    }//cleaning
		    }//deltaPhi
		  }//MET
		}//mu veto
	      }//ele veto
	    }//ge 3 jets
	  }//HT
	}//PV
      }//trigger

   }

     
   //   this doesn't weigh the events at all!!!

   //print cutflow table
   cout<<"unweighted cutflow table"<<endl;
   cout<<"Trigger | "    <<nPassTrigger  <<" +/- "<<sqrt(nPassTrigger)<<endl;
   cout<<"PV      | "    <<nPassPV       <<" +/- "<<sqrt(nPassPV)<<endl;
   cout<<"HT | "         <<nPassHT       <<" +/- "<<sqrt(nPassHT)<<endl;
   cout<<">= 3 jets | "  <<nPass3jets    <<" +/- "<<sqrt(nPass3jets)<<endl;
   cout<<"Ele veto | "   <<nPassEleVeto  <<" +/- "<<sqrt(nPassEleVeto)<<endl;
   cout<<"Mu veto | "    <<nPassMuVeto   <<" +/- "<<sqrt(nPassMuVeto)<<endl;
   cout<<"MET | "        <<nPassMET      <<" +/- "<<sqrt(nPassMET)<<endl;
   cout<<"DeltaPhiN | "  <<nPassDPhi     <<" +/- "<<sqrt(nPassDPhi)<<endl;
//   cout<<"Cleaning | "   <<nPassCleaning <<" +/- "<<sqrt(nPassCleaning)<<endl;
   cout<<">= 1 b | "     <<nGe1b         <<" +/- "<<sqrt(nGe1b)<<endl;
   cout<<"== 1 b | "     <<nEq1b         <<" +/- "<<sqrt(nGe1b)<<endl;
   cout<<">= 2 b | "     <<nGe2b         <<" +/- "<<sqrt(nGe2b)<<endl;
   cout<<">= 3 b | "     <<nGe3b         <<" +/- "<<sqrt(nGe3b)<<endl;
   

}


vector<TString> myCutflow::getFileNames(TString sampleName, bool isTightSelection)//, vector<TString> &fileNames)
{
  vector<TString> fileNames;
//  TString file = "/cu3/galank/SUSY/2011/EventCount/"; file += sampleName; file += "/"; //name of directory to dump files in
  TString file = "";
  file +="cutflow.";
  file += isTightSelection ? "TightSelection." : "BaselineSelection.";
  file += sampleName; file += ".";

  fileNames.push_back(file+"cutHT");
  fileNames.push_back(file+"cut3Jets");
  fileNames.push_back(file+"cutEleVeto");
  fileNames.push_back(file+"cutMuVeto");
  fileNames.push_back(file+"cutMET");
  fileNames.push_back(file+"cutDeltaPhiN");
  fileNames.push_back(file+"cutCleaning");
  fileNames.push_back(file+"cut1b");
  fileNames.push_back(file+"cutEq1b");
  fileNames.push_back(file+"cut2b");
  fileNames.push_back(file+"cut3b");

//  for (unsigned int iCutStage=0; iCutStage<fileNames.size(); iCutStage++){
//    cout<<fileNames[iCutStage]<<endl;
//  }

  return fileNames;
}

vector<ofstream*> myCutflow::getFiles(TString sampleName, bool isTightSelection){
  vector<TString> textfilenames = getFileNames(sampleName,isTightSelection);
  vector<ofstream*> textfiles;
  for (unsigned int istage=0; istage<textfilenames.size(); istage++) textfiles.push_back( new ofstream(textfilenames[istage]) );
  return textfiles;
}
