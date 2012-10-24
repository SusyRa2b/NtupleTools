#define signalEff2012_cxx
#include "signalEff2012.h"
//#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TStopwatch.h"
#include <iostream>

void signalEff2012::Loop()
{

//   In a ROOT session, you can do:
//      Root > .L signalEff2012.C++
//      Root > signalEff2012 t(path, stub, joinbtagbins,usebtagsf);
//everything between reducedTree. and .root is the stub
//      Root > t.Loop();       // Loop on all entries

//define search regions; need to also add SL and LDP

//this will define the order of the columns in the final output! (I think)

  float htedges[]={400,500,800,1000,100000};
  float metedges[]={125,150,250,350,100000};

  vector<SearchRegion> searchregions;

  for (int ildp=0;ildp<2;ildp++) {
    for (int imet = 0;imet<4; imet++) {
      for (int iht = 0;iht<4; iht++) {
	bool theldp = ildp==1 ? true:false;
	if (joinbtagbins_) 
	  searchregions.push_back( SearchRegion(1,99,htedges[iht],htedges[iht+1],metedges[imet],metedges[imet+1],theldp,false));
	else {
	  searchregions.push_back( SearchRegion(1,1,htedges[iht],htedges[iht+1],metedges[imet],metedges[imet+1],theldp,false));
	  searchregions.push_back( SearchRegion(2,2,htedges[iht],htedges[iht+1],metedges[imet],metedges[imet+1],theldp,false));
	  searchregions.push_back( SearchRegion(3,99,htedges[iht],htedges[iht+1],metedges[imet],metedges[imet+1],theldp,false));
	}
      }
    }
  }


  //open an output file
  TString outfilename = "eventcounts.";
  if (joinbtagbins_) outfilename+="mergebbins.";
  outfilename += filestub_;
  outfilename += ".root";
  TFile fout(outfilename,"RECREATE");


  //same indexing as searchregions
  //this ignores that we have to do multiples for things like PDF uncertainties
  //it should work for JES
  vector<TH2D*> eventcounts; 
  for ( size_t ii = 0; ii<searchregions.size(); ii++) {
    TString hname = TString("events_")+searchregions[ii].id();
    //    searchregions[ii].Print();
    eventcounts.push_back( (TH2D*) scanSMSngen_->Clone(hname));
  }

  for ( size_t ii = 0; ii<searchregions.size(); ii++)   {
    eventcounts[ii]->Reset();
    eventcounts[ii]->SetTitle(searchregions[ii].id());
    eventcounts[ii]->Sumw2();
  }
  

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

   Long64_t nentries = fChain->GetEntries(); //BEGIN END jmt mod. no Fast!

   TStopwatch watch;
   watch.Start();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) assert(0); //jmt mod to increase the severity of this (be sure to notice)
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      if (jentry%500000==0) cout<<100*double(jentry)/double(nentries)<<" % done"<<endl;

      // first apply baseline selection
      //njets, jet pT
      if (!( njets>=3 && jetpt2>=70)) continue;

      //cleaning and PV; for MC I will just not even bother to cut on trigger
      if (! ( cutPV && passCleaning)) continue;

      //the buggyEvent cut only matters on data

      //new cleanup cut
      if (! (MET/caloMET<2)) continue;

      //loop over search regions and check whether the event falls into one
      for ( size_t ii = 0; ii<searchregions.size(); ii++) {

	//lepton veto or SL
	if (searchregions[ii].isSL_) {
	  if (! ( ( nMuons+nElectrons==1) && (MT_Wlep<100 && MT_Wlep>0) )) continue;
	}
	else {
	  if (!( nMuons==0 && nElectrons==0 ) ) continue;
	}

	//minDeltaPhiN
	if (searchregions[ii].isLDP_) {
	  if (!( minDeltaPhiN<4)) continue;
	}
	else {
	  if (!( minDeltaPhiN>=4)) continue;
	}

	//note differing use of < versus <=. this is intentional
	
	//passb always true when we're using btag SF
	bool passb = usebtagsf_ || (nbjets >= searchregions[ii].minb_ && nbjets <= searchregions[ii].maxb_);
	bool passmet = MET >= searchregions[ii].minMET_ &&  MET < searchregions[ii].maxMET_;
	bool passht = HT >= searchregions[ii].minHT_ &&  HT < searchregions[ii].maxHT_;

	if (passb && passmet && passht) {
	  //for SMS do not use regular weights
	  float thisweight = PUweight;
	  if (usebtagsf_) { //support only =1,=2,>=3 b bins for now
	    if ( (searchregions[ii].minb_ == searchregions[ii].maxb_) && (searchregions[ii].minb_ == 1) ) 
	      thisweight *= prob1;
	    else  if ( (searchregions[ii].minb_ == searchregions[ii].maxb_) && (searchregions[ii].minb_ == 2) ) 
	      thisweight *= prob2;
	    else  if ( (searchregions[ii].maxb_>10) && (searchregions[ii].minb_ == 3) ) 
	      thisweight *= probge3;
	    else assert(0);
	  }
	  eventcounts[ii]->Fill(m0,m12,thisweight);
	}
      }


   }
   watch.Stop();
   scanSMSngen_->Write();
   fout.Write();
   fout.Close();

   cout<<" done. Wall clock rate = "<<nentries<<" / "<<watch.RealTime()<<" = "<< nentries /watch.RealTime()<<" Hz"<<endl;

}
