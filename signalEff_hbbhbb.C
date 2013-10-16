#define signalEff_hbbhbb_cxx
#include "signalEff_hbbhbb.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TStopwatch.h"
#include <iostream>

#include <map>


void signalEff_hbbhbb::Loop()
{

//   In a ROOT session, you can do:
//      Root > .L signalEff_hbbhbb.C++
//      Root > signalEff_hbbhbb t(path, stub, joinbtagbins,usebtagsf,dopdfs,doPUsyst);
//everything between reducedTree. and .root is the stub
//      Root > t.Loop();       // Loop on all entries


//this will define the order of the columns in the final output! (I think)

  float metsedges[]={30,50,100,150,9999};
  //float metsedges[]={0,30,50,100,150,500};

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

  bool madgraphIsr=false;
  if (true) { //used to have some logic here based on sample name; let's just let it fly like this for now
    //and will use the isr systematic
    assert(theIsrMode_ != kNoIsrWeight);
    madgraphIsr=true;
    cout<<" madgraph Isr Mode enabled"<<endl;
  }

  const int nvariations = 2;
  for (int ivariation=0;ivariation<nvariations;ivariation++) {

    bool thesb = ivariation==1 ? true:false;//SIG first

    for (int imets = 0;imets<4; imets++) {
      for (unsigned int ipdfset = 0; ipdfset<pdfsets.size();ipdfset++) {
        for (int ipdfindex=0; ipdfindex<pdfsetsize[pdfsets.at(ipdfset)]; ipdfindex++) {

          if (joinbtagbins_) 
            searchregions.push_back( SearchRegion(0,metsedges[imets],metsedges[imets+1],thesb,pdfsets[ipdfset],ipdfindex));
          else {
            searchregions.push_back( SearchRegion(2,metsedges[imets],metsedges[imets+1],thesb,pdfsets[ipdfset],ipdfindex));
            searchregions.push_back( SearchRegion(3,metsedges[imets],metsedges[imets+1],thesb,pdfsets[ipdfset],ipdfindex));
            searchregions.push_back( SearchRegion(4,metsedges[imets],metsedges[imets+1],thesb,pdfsets[ipdfset],ipdfindex));
          }

	}//ipdfindex
      }//ipdfset
    }//metsig
  }//ivariation


  //open an output file
  TString outfilename = "eventcounts.";
  if (joinbtagbins_) outfilename+="mergebbins.";
  if (dopdfs_) outfilename+="withpdfs.";
  if (pusyst_) outfilename+="pusyst.";
  if (theIsrMode_==kIsr0) outfilename+="Isr0.";
  else  if (theIsrMode_==kIsrUp) outfilename+="IsrUp.";
  else  if (theIsrMode_==kIsrDown) outfilename+="IsrDown.";
  outfilename += filestub_;
  outfilename += ".root";
  TFile fout(outfilename,"RECREATE");


  //same indexing as searchregions
  vector<TH2D*> eventcounts; 
  vector<TH2D*> eventcountsTotal;  //actually stupid to make N copies of this....
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

  TString hnameallisr = "eventstotalisr";
  TH2D*  eventcountsTotalISR = (TH2D*) scanSMSngen_->Clone(hnameallisr);
  eventcountsTotalISR->Reset();
  eventcountsTotalISR->SetTitle("eventstotal after isr weight");
  eventcountsTotalISR->Sumw2();  

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//ra2b-jmt -- trying to speed things up -- be careful. 
//failure to include a needed variable here leads to bogus results
//for now do it the slow way
  if (dopdfs_ ||true)  fChain->SetBranchStatus("*",1);  //paranoia
  else {
    fChain->SetBranchStatus("*",0);  // disable all branches
    //activate branches
    fChain->SetBranchStatus("m0",1);
    fChain->SetBranchStatus("m12",1);
    fChain->SetBranchStatus("mIntermediate",1);
    fChain->SetBranchStatus("weight3",1); 
    fChain->SetBranchStatus("PUweight",1); 
    fChain->SetBranchStatus("PUweightSystVar",1);
    fChain->SetBranchStatus("njets20",1);
    fChain->SetBranchStatus("jetpt2",1);
    fChain->SetBranchStatus("cutPV",1);
    fChain->SetBranchStatus("passCleaning",1);
    fChain->SetBranchStatus("METsig",1);
    fChain->SetBranchStatus("caloMET",1);
    fChain->SetBranchStatus("MET",1);
    fChain->SetBranchStatus("maxTOBTECjetDeltaMult",1);
    fChain->SetBranchStatus("nMuons",1);
    fChain->SetBranchStatus("nElectrons",1);
    fChain->SetBranchStatus("nTausLoose",1);
    fChain->SetBranchStatus("nIsoPFcands10_010",1);
    fChain->SetBranchStatus("minDeltaPhi20_eta5_noIdAll_nobeta",1);
    fChain->SetBranchStatus("passMC_DiCentralPFJet30_PFMET80_BTagCSV07",1);
    fChain->SetBranchStatus("passMC_DiCentralPFJet30_PFMHT80",1);
    fChain->SetBranchStatus("passMC_PFMET150",1);
    fChain->SetBranchStatus("nbjetsCSVT",1);
    fChain->SetBranchStatus("nbjetsCSVM",1);
    fChain->SetBranchStatus("nbjetsCSVL",1);
    fChain->SetBranchStatus("higgsMbb1MassDiff",1);
    fChain->SetBranchStatus("higgsMbb2MassDiff",1);
    fChain->SetBranchStatus("deltaRmax_hh",1);
 //   fChain->SetBranchStatus("",1);
  
  }

  if (theIsrMode_ !=kNoIsrWeight ) {
    fChain->SetBranchStatus("SUSY_recoilPt",1);
  }

// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries(); //BEGIN END ra2b-jmt mod. no Fast!

   TStopwatch watch;
   watch.Start();
   Long64_t nok=0,nbad=0;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) assert(0); //ra2b-jmt mod to increase the severity of this (be sure to notice)
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

	if (theIsrMode_!=kNoIsrWeight) {
	  float theisrweight=1;
	  int sigmavar=0;
	  if (theIsrMode_==kIsrDown) sigmavar=-1;
	  else if (theIsrMode_==kIsrUp) sigmavar=1;

	  if (madgraphIsr ) {
	    if (SUSY_recoilPt<=120)      theisrweight = 1;
	    else if (SUSY_recoilPt<=150) theisrweight = 0.95 + 0.05*sigmavar;
	    else if (SUSY_recoilPt<=250) theisrweight = 0.9  + 0.1*sigmavar;
	    else                         theisrweight = 0.8  + 0.2*sigmavar;
	  }
	  // Dominick prescription with conservative central value
	  else {//pythia Isr
	    if (SUSY_recoilPt<=40)      theisrweight = 1;
	    else if (SUSY_recoilPt<=100) theisrweight = 1 + 0.06*sigmavar;
	    else if (SUSY_recoilPt<=150) theisrweight = 1  + 0.1*sigmavar;
	    else                         theisrweight = 1  + 0.15*sigmavar;
	  }

	  thisweight *= theisrweight;
	  if (ii==0)  eventcountsTotalISR->Fill(m0,m12,theisrweight); //only fill once per event!
	}

	//TO DO -- add trigger efficiency

	double pdfweight=1;
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
	
	// First apply baseline selection
	if (!cutPV) continue;	
	bool passtrigger = passMC_DiCentralPFJet30_PFMET80_BTagCSV07||passMC_DiCentralPFJet30_PFMHT80||passMC_PFMET150;
	if (!passtrigger) continue;
	if (!passCleaning) continue;

	//new cleanup cut
	if (! (MET/caloMET<2)) continue;
	
	//yet another new cleanup cut
	if (maxTOBTECjetDeltaMult>=40) continue;
       
        //No leptons
        if ( !(nMuons==0&&nElectrons==0&&nTausLoose==0) ) continue;
        if ( !(nIsoPFcands10_010==0) ) continue;

        //Njets cuts
        if ( !(njets20==4 || njets20==5) ) continue;
        if ( !(jetpt2>50) ) continue;

        //minDeltaPhi cuts
        if ( !((minDeltaPhi20_eta5_noIdAll_nobeta>0.5) || (minDeltaPhi20_eta5_noIdAll_nobeta>0.3&&METsig>50)) ) continue;
       
        //deltaR cut 
        if ( !(deltaRmax_hh<2.2) ) continue;
	
	//SB or SIG
	if (searchregions[ii].isSB_) {
	  bool insb = 
	    (0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<90) ||
	    (0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>150) ||
	    (fabs(higgsMbb1MassDiff-higgsMbb2MassDiff)>30);
	  if (!insb) continue;
	}
	else { //signal region
          if ( !((0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140)) ) continue;
          if ( !(fabs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20) ) continue;
	}

	//TO DO -- change the nbjetsCSMX variables to other names,
	//then load in the nbjetsCSVX values into them up here
	//this will allow flexibility when we have SF variations

	bool passb=false;
	if ( (searchregions[ii].nb_ == 4) //2 CSVT + 1 CSVM + 1 CSVL
	     && ((nbjetsCSVT>=2)&&((nbjetsCSVM)>=3)&&((nbjetsCSVL)>=4))) passb=true;
	else if ( (searchregions[ii].nb_ == 3) //2 CSVT + exactly 1 CSVM and 0 additional CSVL
		  && (((nbjetsCSVT>=2)&&(nbjetsCSVM==3)&&nbjetsCSVL==3))) passb=true;
	else if ( (searchregions[ii].nb_ == 2)
		  && (nbjetsCSVT==2&&nbjetsCSVM==2)) passb=true;
	else if ( (searchregions[ii].nb_ == 0) //special case defined as >=2 Tight b-tags
		  && (nbjetsCSVT>=2)) passb=true;

	bool passb2=false;
	if ( (searchregions[ii].nb_ == 4) //2 CSVT + 1 CSVM + 1 CSVL
	     && (CSVbest2>0.898 && CSVbest3>0.679 && CSVbest4>0.244)) passb2=true;
	else if ( (searchregions[ii].nb_ == 3) //2 CSVT + exactly 1 CSVM and 0 additional CSVL
		  && (CSVbest2>0.898 && CSVbest3>0.679 && CSVbest4<=0.244)) passb2=true;
	else if ( (searchregions[ii].nb_ == 2)
		  && (CSVbest2>0.898&& CSVbest3<=0.679)) passb2=true;
	else if ( (searchregions[ii].nb_ == 0) //special case defined as >=2 Tight b-tags
		  && (CSVbest2>0.898)) passb2=true;
	
	if (passb != passb2) cout<<" problem with btag counting!"<<endl;

	bool passmets = (METsig >= searchregions[ii].minMETs_ &&  METsig < searchregions[ii].maxMETs_);

	if ( passmets && passb ) eventcounts[ii]->Fill(m0,m12,thisweight);
	
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


