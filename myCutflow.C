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

   float minHT = isTightSelection ? 500 : 400; 
   float minMET = isTightSelection ? 300 : 250;

   bool geq1b = false;
   bool geq2b = true;
   if(isTightSelection && geq1b) {minHT = 500; minMET=500;}
   else if(isTightSelection && geq2b) {minHT = 600; minMET=300;}

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
   int nPassCleaning = 0;

   int nPassScapingVeto = 0;
   int nPassHBHENoise = 0;
   int nPassCSCTightHalo = 0;
   int nPassTrackingFailure = 0;
   int nPassECALTP = 0;
   int nPassEENoise = 0;
   int nPassGreedyMuon = 0;
   int nPassInconsistentMuon = 0;


   int nCleanFirstFive = 0;
   int nCleanLastTwo = 0;


   int nEq1b = 0;
   int nGe1b = 0;
   int nGe2b = 0;
   int nGe3b = 0;

   float nGe1bProb = 0;
   float nGe2bProb = 0;
   float nGe3bProb = 0;


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




      //if(!(runNumber==177319 && eventNumber==369647624 && lumiSection==227)) continue;
      /*
      if(runNumber==162762 && lumiSection==11 && eventNumber==3920149    ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==162762 && lumiSection==35 && eventNumber==13086213   ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==162762 && lumiSection==44 && eventNumber==16726674	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==162762 && lumiSection==52 && eventNumber==19595141	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==162762 && lumiSection==79 && eventNumber==30142641	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==162803 && lumiSection==127 && eventNumber==43788905	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==162803 && lumiSection==135 && eventNumber==48357843	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==162803 && lumiSection==77 && eventNumber==16104566	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==162808 && lumiSection==21 && eventNumber==10650481	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==162808 && lumiSection==28 && eventNumber==14265458	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==162808 && lumiSection==43 && eventNumber==22148878	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==162808 && lumiSection==7 && eventNumber==2883206	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==162811 && lumiSection==107 && eventNumber==49276600	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==162811 && lumiSection==10 && eventNumber==4296297	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==162811 && lumiSection==31 && eventNumber==14269926	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==162811 && lumiSection==48 && eventNumber==21898728	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==162811 && lumiSection==98 && eventNumber==45341778	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163255 && lumiSection==918 && eventNumber==574758535 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163270 && lumiSection==295 && eventNumber==184599363 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163302 && lumiSection==164 && eventNumber==82813919	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163332 && lumiSection==110 && eventNumber==63153304	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163332 && lumiSection==443 && eventNumber==300924904 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163333 && lumiSection==2 && eventNumber==877665	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163340 && lumiSection==207 && eventNumber==104822592 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163340 && lumiSection==330 && eventNumber==165764521 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163385 && lumiSection==269 && eventNumber==156269438 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163402 && lumiSection==98 && eventNumber==59226016	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163588 && lumiSection==127 && eventNumber==86700074	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163657 && lumiSection==92 && eventNumber==76367690	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163659 && lumiSection==517 && eventNumber==379050220 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163738 && lumiSection==278 && eventNumber==227402021 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163738 && lumiSection==56 && eventNumber==43894899	 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163758 && lumiSection==454 && eventNumber==323038090 ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      if(runNumber==163759 && lumiSection==84 && eventNumber==55120943   ) std::cout << "at discrepant event, HT=" << HT << std::endl;
      */

      //std::cout << runNumber << ":"<<lumiSection<<":"<<eventNumber<<", HT = " << HT << ", njets = " << njets << std::endl;

      //if( !passCleaning ) continue;

      /*
      //print out discrepant online calo MHT /  offline pfMET events
      //ht300mht80v1andv2
      if(cutTrigger==0 && pass_utilityHLT==1 && HT>=400 && runNumber>=166979 && runNumber<=173211 && MET>200)
	std::cout << runNumber << ":"<< lumiSection << ":" << eventNumber<<", MET = " << MET << std::endl;
      //ht300mht90v2
      if(cutTrigger==0 && pass_utilityHLT==1 && HT>=400 && runNumber>=173212 && runNumber<=176544 && MET>200)
	std::cout << runNumber << ":"<< lumiSection << ":" << eventNumber<<", MET = " << MET << std::endl;
      //ht350mht90v1
      if(cutTrigger==0 && pass_utilityHLT==1 && HT>=400 && runNumber>=176545 && runNumber<=178410 && MET>200)
	std::cout << runNumber << ":"<< lumiSection << ":" << eventNumber<<", MET = " << MET << std::endl;
      //ht350mht110v1
      if(cutTrigger==0 && pass_utilityHLT==1 && HT>=400 && runNumber>=178411 && MET>200)
	std::cout << runNumber << ":"<< lumiSection << ":" << eventNumber<<", MET = " << MET << std::endl;
      */

      /*
      if( scrapingvetoFilter==1 && hbhenoiseFilter==1 && csctighthaloFilter==1 && trackingfailureFilterPFLOW==1 && ra2ecaltpFilter==1){
	if( greedymuonFilter==1 && inconsistentmuonFilter==1){

      if(cutTrigger==0 && pass_utilityHLT==1 && HT>=400 && MET>200)
	std::cout << runNumber << ":"<< lumiSection << ":" << eventNumber<<", MET = " << MET << std::endl;
	}
      }
      */
      int istage=0;


//      //for met cleaning table
//      if( scrapingvetoFilter==1){
//        nPassScapingVeto++;
//        if( hbhenoiseFilter==1){
//          nPassHBHENoise++;
//          if(csctighthaloFilter==1){
//            nPassCSCTightHalo++;
//            if(trackingfailureFilterPFLOW==1){
//	      nPassTrackingFailure++;
//	      if(ra2ecaltpFilter==1){
//		nPassECALTP++;
//
      /*______________________________________________________ trigger ______________________________________*/

      if (cutTrigger){
      //      if (pass_utilityHLT==1){
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
	      //	      if((nElectrons==0 && nMuons==1)||(nElectrons==1 && nMuons==0) && MT_Wlep>=0 && MT_Wlep<100){
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
//		  if (200<MET && MET<=250){
//		  if (50<=MET && MET<100){
		    nPassMET++;
		    if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
		    istage++;

      /*______________________________________________________ min dPhi / min dPhiN __________________________*/

//		    if (cutDeltaPhi){
		    if (minDeltaPhiN>=4){
//		    if (minDeltaPhiN<4){
		      nPassDPhi++;
		      if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
		      istage++;

      /*______________________________________________________ cleaning ______________________________________*/

//		  if (passCleaning){
//		    nPassCleaning++;
//		    if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
//		    istage++;
//
      /*______________________________________________________ nB ____________________________________________*/
//			if( scrapingvetoFilter==1 && hbhenoiseFilter==1 && csctighthaloFilter==1 && trackingfailureFilterPFLOW==1 && ra2ecaltpFilter==1){
//			  nCleanFirstFive++;
//
//			  if( greedymuonFilter==1 && inconsistentmuonFilter==1){
//			    nCleanLastTwo++;
//


		      //std::cout << runNumber<< " "<< lumiSection<<" " << eventNumber <<" "<< probge1<<" " << probge2 <<" "<< probge3 << std::endl;

		      nGe1bProb += probge1;
		      nGe2bProb += probge2;
		      nGe3bProb += probge3;



		      if (nbjetsCSVM>=1) {
		      //if (nbjetsCSVM==0) {
			nGe1b++;
			//nGe1b += weight;
			if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
			//std::cout <<runNumber<<" "<<eventNumber<<" "<<lumiSection<< std::endl;
			istage++;
			if (nbjetsCSVM==1) {
			  nEq1b++;
			  if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
			}
			istage++;

			//for met cleaning table
			if( scrapingvetoFilter==1){
			  nPassScapingVeto++;
			  if( hbhenoiseFilter==1){
			    nPassHBHENoise++;
			    if(csctighthaloFilter==1){
			      nPassCSCTightHalo++;
			      if(trackingfailureFilterPFLOW==1){
				nPassTrackingFailure++;
				if(ra2ecaltpFilter==1){
				  nPassECALTP++;
				  if(eenoiseFilter==1){
				    nPassEENoise++;
				    if( greedymuonFilter==1){
				      nPassGreedyMuon++; 
				      if(inconsistentmuonFilter==1){
					nPassInconsistentMuon++;
			//	      }
			//	    }
			//	  }
			//	}
			//      }
			//    }
			//  }
			//}



			if (nbjetsCSVM>=2){
			  nGe2b++;
			  std::cout << runNumber<< ":"<< lumiSection<<":" << eventNumber << std::endl;

			  //nGe2b += weight;
			  if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
			  istage++;
			  if (nbjetsCSVM>=3) {
			    nGe3b++;
			    //nGe3b+=weight;
			    if (writeFiles) *textfiles[istage]<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
			    istage++;
			  }
			}
				      }}}}}}}}
			//}//cleaning - last two
			//}//cleaning - first five

		      }//n b jets
		      //}//cleaning
		    }//deltaPhi
		  }//MET
		}//mu veto
	      }//ele veto
	    }//ge 3 jets
	  }//HT
	}//PV
      }//trigger
      //      	}}}}}
   }

     
   //   this doesn't weigh the events at all!!!

   //print cutflow table
   //float W = weight*4684;
   //float W = weight*5000;
   //float W = 158*5000./3701947.;
   float W = 1;
   //float W = 5000;


   cout<<"unweighted cutflow table"<<endl;
   cout<<"Trigger | "    <<nPassTrigger  *W <<" +/- "<<sqrt(nPassTrigger)*W<<endl;
   cout<<"PV      | "    <<nPassPV       *W <<" +/- "<<sqrt(nPassPV)     *W<<endl;
   cout<<"HT | "         <<nPassHT       *W <<" +/- "<<sqrt(nPassHT)     *W<<endl;
   cout<<">= 3 jets | "  <<nPass3jets    *W <<" +/- "<<sqrt(nPass3jets)  *W<<endl;
   cout<<"Ele veto | "   <<nPassEleVeto  *W <<" +/- "<<sqrt(nPassEleVeto)*W<<endl;
   cout<<"Mu veto | "    <<nPassMuVeto   *W <<" +/- "<<sqrt(nPassMuVeto) *W<<endl;
   cout<<"MET | "        <<nPassMET      *W <<" +/- "<<sqrt(nPassMET)    *W<<endl;
   cout<<"DeltaPhiN | "  <<nPassDPhi     *W <<" +/- "<<sqrt(nPassDPhi)   *W<<endl;
//   cout<<"Cleaning | "   <<nPassCleaning <<" +/- "<<sqrt(nPassCleaning)<<endl;
   cout<<">= 1 b | "     <<nGe1b         *W <<" +/- "<<sqrt(nGe1b)*W<<endl;


   cout << "nPassScapingVeto	  | " <<  nPassScapingVeto	 *W << " +/- " << sqrt(nPassScapingVeto	    )*W <<endl;
   cout << "nPassHBHENoise	  | " <<  nPassHBHENoise	 *W << " +/- " << sqrt(nPassHBHENoise	    )*W <<endl;
   cout << "nPassCSCTightHalo	  | " <<  nPassCSCTightHalo	 *W << " +/- " << sqrt(nPassCSCTightHalo    )*W <<endl;	 
   cout << "nPassTrackingFailure  | " <<  nPassTrackingFailure 	 *W << " +/- " << sqrt(nPassTrackingFailure )*W <<endl;	 
   cout << "nPassECALTP		  | " <<  nPassECALTP		 *W << " +/- " << sqrt(nPassECALTP	    )*W <<endl;	 
   cout << "nPassEENoise	  | " <<  nPassEENoise		 *W << " +/- " << sqrt(nPassEENoise	    )*W <<endl;	 
   cout << "nPassGreedyMuon	  | " <<  nPassGreedyMuon	 *W << " +/- " << sqrt(nPassGreedyMuon	    )*W <<endl;
   cout << "nPassInconsistentMuon | " <<  nPassInconsistentMuon  *W << " +/- " << sqrt(nPassInconsistentMuon)*W <<endl;  





   cout<<"Clean(first five) | "    <<nCleanFirstFive*W  <<" +/- "<<sqrt(nCleanFirstFive)*W<<endl;
   cout<<"Clean(last two) | "    <<nCleanLastTwo*W  <<" +/- "<<sqrt(nCleanLastTwo)*W<<endl;
   cout<<"== 1 b | "     <<nEq1b*W         <<" +/- "<<sqrt(nEq1b)*W<<endl;
   cout<<">= 2 b | "     <<nGe2b*W         <<" +/- "<<sqrt(nGe2b)*W<<endl;
   cout<<">= 3 b | "     <<nGe3b*W         <<" +/- "<<sqrt(nGe3b)*W<<endl;
   
   cout<<">= 1 b (Prob)| "     <<nGe1bProb*W         <<" +/- "<<sqrt(nGe1bProb)*W<<endl;
   cout<<">= 2 b (Prob)| "     <<nGe2bProb*W         <<" +/- "<<sqrt(nGe2bProb)*W<<endl;
   cout<<">= 3 b (Prob)| "     <<nGe3bProb*W         <<" +/- "<<sqrt(nGe3bProb)*W<<endl;
   

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
  //fileNames.push_back(file+"cutCleaning");
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
