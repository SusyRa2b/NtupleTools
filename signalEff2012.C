#define signalEff2012_cxx
#include "signalEff2012.h"
//#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TStopwatch.h"
#include <iostream>

#include <map>

void signalEff2012::Loop()
{

//   In a ROOT session, you can do:
//      Root > .L signalEff2012.C++
//      Root > signalEff2012 t(path, stub, joinbtagbins,usebtagsf,dopdfs,doPUsyst);
//everything between reducedTree. and .root is the stub
//      Root > t.Loop();       // Loop on all entries

//define search regions; need to also add SL and LDP

//this will define the order of the columns in the final output! (I think)

  float htedges[]={400,500,800,1000,100000};
  float metedges[]={125,150,250,350,100000};

  vector<TString> pdfsets;
  map<TString, int> pdfsetsize;
  pdfsets.push_back("none");
  pdfsetsize["none"]=1;
  if ( dopdfs_) {
//     pdfsets.push_back("CTEQ");
//     pdfsets.push_back("MSTW");
//     pdfsets.push_back("NNPDF");
    pdfsets.push_back("CTEQMSTW");
    pdfsets.push_back("CTEQNNPDF");
//     pdfsetsize["CTEQ"]=45;
//     pdfsetsize["MSTW"]=41;
//     pdfsetsize["NNPDF"]=100;
    pdfsetsize["CTEQMSTW"]=1;
    pdfsetsize["CTEQNNPDF"]=1;
  }

  vector<SearchRegion> searchregions;

  //not incredibly elegant, but it will do
  const bool doLeptons =  !(filestub_.Contains("T1bbbb"));
  const int nvariations = doLeptons ? 4 : 2;
  for (int ivariation=0;ivariation<nvariations;ivariation++) {

    bool thesl=false,theslsig=false,theldp=false;

    if (doLeptons ) {
      theslsig = ivariation==1 ? true:false;
      thesl = ivariation==2 ? true:false;
      theldp = ivariation==3 ? true:false;
    }
    else {
      theldp = ivariation==1 ? true:false;
    }

    for (int imet = 0;imet<4; imet++) {
      for (int iht = 0;iht<4; iht++) {

	for (unsigned int ipdfset = 0; ipdfset<pdfsets.size();ipdfset++) {
	  for (int ipdfindex=0; ipdfindex<pdfsetsize[pdfsets.at(ipdfset)]; ipdfindex++) {

	    if (joinbtagbins_) 
	      searchregions.push_back( SearchRegion(1,99,htedges[iht],htedges[iht+1],metedges[imet],metedges[imet+1],theldp,thesl,theslsig,minnjets_,pdfsets[ipdfset],ipdfindex));
	    else {
	      searchregions.push_back( SearchRegion(1,1,htedges[iht],htedges[iht+1],metedges[imet],metedges[imet+1],theldp,thesl,theslsig,minnjets_,pdfsets[ipdfset],ipdfindex));
	      searchregions.push_back( SearchRegion(2,2,htedges[iht],htedges[iht+1],metedges[imet],metedges[imet+1],theldp,thesl,theslsig,minnjets_,pdfsets[ipdfset],ipdfindex));
	      searchregions.push_back( SearchRegion(3,99,htedges[iht],htedges[iht+1],metedges[imet],metedges[imet+1],theldp,thesl,theslsig,minnjets_,pdfsets[ipdfset],ipdfindex));
	    }

	  }
	}
      }
    }
  }


  //open an output file
  TString outfilename = "eventcounts.";
  if (joinbtagbins_) outfilename+="mergebbins.";
  if (dopdfs_) outfilename+="withpdfs.";
  if (pusyst_) outfilename+="pusyst.";
  if (minnjets_ !=3) {
    outfilename += "minnjets";
    outfilename += minnjets_;
    outfilename +=".";
  }
  outfilename += filestub_;
  outfilename += ".root";
  TFile fout(outfilename,"RECREATE");


  //same indexing as searchregions
  vector<TH2D*> eventcounts; 
  vector<TH2D*> eventcountsTotal; 
  for ( size_t ii = 0; ii<searchregions.size(); ii++) {
    TString hname = TString("events_")+searchregions[ii].id();
    TString hnameall = TString("eventstotal_")+searchregions[ii].id();
    //    searchregions[ii].Print();
    eventcounts.push_back( (TH2D*) scanSMSngen_->Clone(hname));
    eventcountsTotal.push_back( (TH2D*) scanSMSngen_->Clone(hnameall));
    cout<<searchregions[ii].id()<<endl;
  }

  for ( size_t ii = 0; ii<searchregions.size(); ii++)   {
    eventcounts[ii]->Reset();
    eventcounts[ii]->SetTitle(searchregions[ii].id());
    eventcounts[ii]->Sumw2();

    eventcountsTotal[ii]->Reset();
    eventcountsTotal[ii]->SetTitle(searchregions[ii].id()); //i think this is ok (name will still be unique)
    eventcountsTotal[ii]->Sumw2();
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
//jmt -- trying to speed things up -- be careful. 
//failure to include a needed variable here leads to bogus results
  if (dopdfs_)  fChain->SetBranchStatus("*",1);  //paranoia
  else {
    fChain->SetBranchStatus("*",0);  // disable all branches
    //activate branches
    fChain->SetBranchStatus("m0",1);
    fChain->SetBranchStatus("m12",1);
    fChain->SetBranchStatus("weight3",1); 
    fChain->SetBranchStatus("PUweight",1); 
    fChain->SetBranchStatus("PUweightSystVar",1);
    fChain->SetBranchStatus("prob1",1);
    fChain->SetBranchStatus("prob2",1);
    fChain->SetBranchStatus("probge3",1);
    //   if (dopdfs_) {
    //     fChain->SetBranchStatus("pdfWeightsCTEQ",1);
    //     fChain->SetBranchStatus("pdfWeightsMSTW",1);
    //     fChain->SetBranchStatus("pdfWeightsNNPDF",1);
    //   }
    fChain->SetBranchStatus("njets",1);
    fChain->SetBranchStatus("jetpt2",1);
    fChain->SetBranchStatus("cutPV",1);
    fChain->SetBranchStatus("passCleaning",1);
    fChain->SetBranchStatus("MET",1);
    fChain->SetBranchStatus("HT",1);
    fChain->SetBranchStatus("caloMET",1);
    fChain->SetBranchStatus("maxTOBTECjetDeltaMult",1);
    fChain->SetBranchStatus("nMuons",1);
    fChain->SetBranchStatus("nElectrons",1);
    fChain->SetBranchStatus("MT_Wlep",1);
    fChain->SetBranchStatus("nIsoTracks15_005_03",1);
    fChain->SetBranchStatus("minDeltaPhiN_asin",1);
    fChain->SetBranchStatus("nbjets",1);
    //    fChain->SetBranchStatus("MT_bestCSV",1);
  }

// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries(); //BEGIN END jmt mod. no Fast!

   TStopwatch watch;
   watch.Start();
   Long64_t nok=0,nbad=0;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) assert(0); //jmt mod to increase the severity of this (be sure to notice)
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      if (jentry%500000==0) {
	cout<<100*double(jentry)/double(nentries)<<" % done"<<endl;
	 if (dopdfs_) cout<<" bad fraction = "<<double(nbad)/ double (nbad+nok)<<endl;
      }

      //loop over search regions and check whether the event falls into one
      for ( size_t ii = 0; ii<searchregions.size(); ii++) {

	// == get event weight ==
	
	//Use PU weight
	double thisweight = pusyst_ ? PUweightSystVar : PUweight;
	thisweight*=weight3; //for SMS this is just 1; for other samples it is 1 pb-1 of data
	double pdfweight=1;
	//use b tag SF if desired
	if (usebtagsf_) { //support only =1,=2,>=3 b bins for now
	  if ( (searchregions[ii].minb_ == searchregions[ii].maxb_) && (searchregions[ii].minb_ == 1) ) 
	    thisweight *= prob1;
	  else  if ( (searchregions[ii].minb_ == searchregions[ii].maxb_) && (searchregions[ii].minb_ == 2) ) 
	    thisweight *= prob2;
	  else  if ( (searchregions[ii].maxb_>10) && (searchregions[ii].minb_ == 3) ) 
	    thisweight *= probge3;
	  else assert(0);
	}
	//if this search region is a pdf variation then use pdf weight
	if ( searchregions[ii].pdfset_ == "CTEQ" ) { //for speed reasons i am scared of doing all of these string comparisons. still outweighed by i/o speed?
	  pdfweight = pdfWeightsCTEQ[searchregions[ii].pdfindex_];
	}
	else if (searchregions[ii].pdfset_ == "MSTW" ) {
	  pdfweight = pdfWeightsMSTW[searchregions[ii].pdfindex_];
	}
	else if (searchregions[ii].pdfset_ == "NNPDF" ) {
	  pdfweight = pdfWeightsNNPDF[searchregions[ii].pdfindex_];
	}
	else if (searchregions[ii].pdfset_ == "CTEQMSTW" ) {
	  assert(searchregions[ii].pdfindex_==0);
	  double w1 = pdfWeightsMSTW[0];
	  double w2 = pdfWeightsCTEQ[0];
	  //	  pdfweight = pdfWeightsMSTW[searchregions[ii].pdfindex_]/pdfWeightsCTEQ[searchregions[ii].pdfindex_];
	  if (w1<=0 || w2<=0) {
	    //	    cout<<w1<<" "<<w2<<endl;
	    nbad++;
	    continue;//skip the event!
	  }
	  else ++nok;
	  pdfweight = w1/w2;
	}
	else if (searchregions[ii].pdfset_ == "CTEQNNPDF" ) {
	  assert(searchregions[ii].pdfindex_==0);
	  double w1 = pdfWeightsNNPDF[0];
	  double w2 = pdfWeightsCTEQ[0];
	  if (w1<=0 || w2<=0) {
	    nbad++;
	    continue;//skip the event!
	  }
	  else ++nok;
	  pdfweight = w1/w2;
	}
	
	
	//Fill the "total number of events" histogram, but ignore the PU and btag weights
	eventcountsTotal[ii]->Fill(m0,m12,pdfweight);
	
	thisweight *= pdfweight;
	
	// == now apply cuts ==
	
	// ONLY FOR T2TT and only for tests!
	//if ( ! (MT_bestCSV>150 )) continue;

	// first apply baseline selection
	//njets, jet pT
	if (!( njets>= minnjets_ && jetpt2>=70)) continue;
	
	//cleaning and PV; for MC I will just not even bother to cut on trigger
	if (!cutPV) continue;	

	//HACK -- there is a bug in some of my trees. Disable cleaning as a work-around, except for the nominal case
	if (!joinbtagbins_ && !passCleaning) continue;

	//the buggyEvent cut only matters on data

	//new cleanup cut
	if (! (MET/caloMET<2)) continue;
	
	//yet another new cleanup cut
	if (maxTOBTECjetDeltaMult>=40) continue;
	
	//lepton veto or SL
	if (searchregions[ii].isSL_) {
	  if (! ( ( nMuons+nElectrons==1) && (MT_Wlep<100 && MT_Wlep>0) )) continue;
	}
	else if (searchregions[ii].isSLsig_) {
	  if (! ( ( nMuons+nElectrons==1) && (MT_Wlep>=100 ) )) continue;
	}
	else { //ldp or signal region
	  if (!( nMuons==0 && nElectrons==0 ) ) continue;
	  if ( nIsoTracks15_005_03>0) continue; //isolated track veto for any region that is not SL
	}

	//minDeltaPhiN
	if (searchregions[ii].isLDP_) {
	  if (!( minDeltaPhiN_asin<4)) continue;
	}
	else {
	  if (!( minDeltaPhiN_asin>=4)) continue;
	}

	//note differing use of < versus <=. this is intentional
	
	//passb always true when we're using btag SF
	bool passb = usebtagsf_ || (nbjets >= searchregions[ii].minb_ && nbjets <= searchregions[ii].maxb_);
	bool passmet = MET >= searchregions[ii].minMET_ &&  MET < searchregions[ii].maxMET_;
	bool passht = HT >= searchregions[ii].minHT_ &&  HT < searchregions[ii].maxHT_;

	if (passb && passmet && passht) eventcounts[ii]->Fill(m0,m12,thisweight);
	
      }


   }
   watch.Stop();
   scanSMSngen_->Write();
   for (unsigned int ih=0;ih<scanProcessTotals_.size();ih++)   scanProcessTotals_[ih]->Write();
   fout.Write();
   fout.Close();
   cout<<" bad / good "<<nbad<<" "<<nok<<endl;

   cout<<" done. Wall clock rate = "<<nentries<<" / "<<watch.RealTime()<<" = "<< nentries /watch.RealTime()<<" Hz"<<endl;

}
