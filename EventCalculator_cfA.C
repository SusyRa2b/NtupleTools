#include "EventCalculator_cfA.h"

#include "PUConstants.h"

#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TFile.h"
#include "TMatrixT.h"
#include "TMatrixDEigen.h"
#include "TStopwatch.h"
#include "TRegexp.h"

#include "TVector3.h"

#include "inJSON2012.h"

#include "MiscUtil.cxx"

//#include <RooRealVar.h>
//#include <RooMinuit.h>
//#include <RooTransverseThrustVar.h>

#include <cassert>
#include <fstream>
#include <algorithm>

using namespace std;


EventCalculator::EventCalculator(const TString & sampleName, const vector<string> inputFiles, jetType theJetType, METType theMETType) :
  sampleName_(sampleName),
  sampleIsSignal_(false),
  theScanType_(kNotScan),
  theMETType_(theMETType),
  theJetType_(theJetType),
  theJESType_(kJES0),
  theJERType_(kJER0),
  theMETuncType_(kMETunc0),
  thePUuncType_(kPUunc0),
  theBTagEffType_(kBTagEff04),
  theHLTEffType_(kHLTEff0),
  theBTaggerType_(kCSVM),
  crossSectionTanb40_10_(0),
  crossSectionTanb40_05_(0),
  crossSectionTanb40_20_(0),
  smsCrossSectionFile_(0),
  f_eff_ht300_(0),
  f_eff_ht350_(0),
  htgraph_ht300_(0),
  htgraphPlus_ht300_(0),
  htgraphMinus_ht300_(0),
  htgraph_ht350_(0),
  htgraphPlus_ht350_(0),
  htgraphMinus_ht350_(0),
  f_eff_mht_(0),
  mhtgraph_(0),
  mhtgraphPlus_(0),
  mhtgraphMinus_(0),
  f_tageff_(0),
  fDataJetRes_(0),
  //2011 values, to be used for code validation only
  ResJetPar_(new JetCorrectorParameters("JESfiles/START53_V7F_L2L3Residual_AK5PFchs.txt") ),
  L3JetPar_(new JetCorrectorParameters("JESfiles/START53_V7F_L3Absolute_AK5PFchs.txt") ),
  L2JetPar_(new JetCorrectorParameters("JESfiles/START53_V7F_L2Relative_AK5PFchs.txt") ),
  L1JetPar_(new JetCorrectorParameters("JESfiles/START53_V7F_L1FastJet_AK5PFchs.txt") ),
  JetCorrector_(0 ),
  jecUnc_(new JetCorrectionUncertainty("JESfiles/START53_V7F_Uncertainty_AK5PFchs.txt")),
  starttime_(0),
  recalculatedVariables_(false),
  watch_(0),
  //new cfA stuff
  chainB( new TChain("/configurableAnalysis/eventB")),
  chainV( 0),
  chainA( new TChain("configurableAnalysis/eventA"))
{

  if ( sampleName_.Contains("mSUGRA") ) {
    theScanType_ = kmSugra;
    std::cout<<"\tDetected that I'm running over an mSugra scan!"<<std::endl;
    sampleIsSignal_=true;
  }
  else if (sampleName_.Contains("T1bbbb") || sampleName_.Contains("T2bb") || sampleName_.Contains("T2tt") || sampleName_.Contains("T1tttt")) {
    theScanType_ = kSMS;
    std::cout<<"\tDetected that I'm running over an SMS scan!"<<std::endl;
    sampleIsSignal_=true;
  }
  else if (sampleName_.Contains("LM")) {
    sampleIsSignal_=true;
    std::cout<<"\tDetected that I'm running over a SUSY sample!"<<std::endl;
  }
  else std::cout<<"\tRunning on a SM sample"<<std::endl;

  initEnumNames();

  loadSusyScanCrossSections();

  //  checkConsistency();

  //loadHLTMHTeff();
  loadHLTHTeff();
  loadDataJetRes();

  loadECALStatus();//could also be in reducedTree

  vPar_.clear(); //should be overkill
  vPar_.push_back(*L1JetPar_);
  vPar_.push_back(*L2JetPar_);
  vPar_.push_back(*L3JetPar_);
  vPar_.push_back(*ResJetPar_);
  JetCorrector_ =  new FactorizedJetCorrector(vPar_);

  loadJetTagEffMaps();  

  //finally, initialize the TBraches for cfA
  for (unsigned int ii = 0; ii<inputFiles.size(); ++ii) {
    chainB->Add(inputFiles.at(ii).c_str());
    chainA->Add(inputFiles.at(ii).c_str());
  }

  InitializeB(chainB);
  InitializeA(chainA);

  cout<<sampleName_<<endl;

}

EventCalculator::~EventCalculator() {
  //turns out the program does a segv at terminate if i don't properly delete these
  delete chainA;
  delete chainB;
  cout<<"EventCalculator dead"<<endl;
}


void EventCalculator::initEnumNames() {

  theBTaggerNames_[kSSVM]="SSVHEM";
  theBTaggerNames_[kTCHET]="TCHET";
  theBTaggerNames_[kSSVHPT]="SSVHPT";
  theBTaggerNames_[kTCHPT]="TCHPT";
  theBTaggerNames_[kTCHPM]="TCHPM";
  theBTaggerNames_[kCSVM]="CSVM";

  theMETNames_[kPFMET] = "PFMET";
  theMETNames_[kPFMETTypeI] = "PFMETTypeI";

  theJetNames_[kPF2PAT]="PF2PATjets";
  theJetNames_[kRECOPF]="RecoPFjets";

  theJESNames_[kJES0]="JES0";
  theJESNames_[kJESup]="JESup";
  theJESNames_[kJESdown]="JESdown";
  //  theJESNames_[kJESFLY]="JESFLY"; //disable this

  theMETuncNames_[kMETunc0]="METunc0";
  theMETuncNames_[kMETuncUp]="METuncUp";
  theMETuncNames_[kMETuncDown]="METuncDown";

  thePUuncNames_[kPUunc0]="PUunc0";
  thePUuncNames_[kPUuncUp]="PUuncUp";
  thePUuncNames_[kPUuncDown]="PUuncDown";

  theBTagEffNames_[kBTagEff0]="BTagEff0";
  theBTagEffNames_[kBTagEffup]="BTagEffup";
  theBTagEffNames_[kBTagEffdown]="BTagEffdown";
  theBTagEffNames_[kBTagEff02]="BTagEff02";
  theBTagEffNames_[kBTagEffup2]="BTagEffup2";
  theBTagEffNames_[kBTagEffdown2]="BTagEffdown2";
  theBTagEffNames_[kBTagEff03]="BTagEff03";
  theBTagEffNames_[kBTagEffup3]="BTagEffup3";
  theBTagEffNames_[kBTagEffdown3]="BTagEffdown3";
  theBTagEffNames_[kBTagEff04]="BTagEff04";
  theBTagEffNames_[kBTagEffup4]="BTagEffup4";
  theBTagEffNames_[kBTagEffdown4]="BTagEffdown4";


  theHLTEffNames_[kHLTEff0]="HLTEff0";
  theHLTEffNames_[kHLTEffup]="HLTEffup";
  theHLTEffNames_[kHLTEffdown]="HLTEffdown";

  theJERNames_[kJER0]="JER0";
  theJERNames_[kJERup]="JERup";
  theJERNames_[kJERbias]="JERbias";
  theJERNames_[kJERra2]="JERra2";
  theJERNames_[kJERdown]="JERdown";
 
}



void EventCalculator::stopTimer(const Long64_t ntotal) {
  TDatime stoptime; //default ctor is for current time
  double elapsed= stoptime.Convert() - starttime_->Convert();
  std::cout<<"events / time = "<<ntotal<<" / "<<elapsed<<" = "<<double(ntotal)/double(elapsed)<<" Hz"<<std::endl;
  
  delete starttime_;
  starttime_=0;
}

bool EventCalculator::isSampleRealData() {

  if (sampleName_.BeginsWith("HT_Run2012A")) return true;
  if (sampleName_.BeginsWith("HTMHT_Run2012")) return true;
  if (sampleName_.BeginsWith("JetHT_Run2012")) return true;

  if (sampleName_.BeginsWith("MET_Run2012")) return true;

  if (sampleName_.BeginsWith("MuEG_Run2012")) return true;

  if (sampleName_.BeginsWith("DoubleElectron_Run2012")) return true;

  if (sampleName_.BeginsWith("DoubleMu_Run2012")) return true;

  if (sampleName_.BeginsWith("ElectronHad_Run2012")) return true;
  if (sampleName_.BeginsWith("MuHad_Run2012")) return true;

  if (sampleName_.BeginsWith("SingleMu_Run2012")) return true;
  if (sampleName_.BeginsWith("SingleElectron_Run2012")) return true;

  if (sampleName_.BeginsWith("Photon_Run2012")) return true;

  if (sampleName_.BeginsWith("PhotonHad_Run2012")) return true;

  return false;
}

void  EventCalculator::loadSusyScanCrossSections() {

  if (theScanType_ == kmSugra ) {
    crossSectionTanb40_10_ = new CrossSectionTable("NLOxsec_tanb40_10.txt");
    crossSectionTanb40_05_ = new CrossSectionTable("NLOxsec_tanb40_05.txt");
    crossSectionTanb40_20_ = new CrossSectionTable("NLOxsec_tanb40_20.txt");
  }

}

void EventCalculator::setOptions( const TString & opt) {
  //cannot set the b tagger, or the jet type, or the met type, or the scan type
  //but all other options must be set here

  if (opt=="") return;

  cout<<opt<<endl;

  //i wish i could think of a more clever way to code this
  if ( getOptPiece("JES",opt)== theJESNames_[kJES0]) theJESType_ = kJES0;
  else  if ( getOptPiece("JES",opt)== theJESNames_[kJESup]) theJESType_ = kJESup;
  else  if ( getOptPiece("JES",opt)== theJESNames_[kJESdown]) theJESType_ = kJESdown;
  //  else  if ( getOptPiece("JES",opt)== theJESNames_[kJESFLY]) theJESType_ = kJESFLY;
  else {cout<<"problem with opt in JES"<<endl; assert(0) ;} //enforce a complete set of options!

  if ( getOptPiece("JER",opt)== theJERNames_[kJER0]) theJERType_ = kJER0;
  else  if ( getOptPiece("JER",opt)== theJERNames_[kJERup]) theJERType_ = kJERup;
  else  if ( getOptPiece("JER",opt)== theJERNames_[kJERdown]) theJERType_ = kJERdown;
  else  if ( getOptPiece("JER",opt)== theJERNames_[kJERbias]) theJERType_ = kJERbias;
  else  if ( getOptPiece("JER",opt)== theJERNames_[kJERra2]) theJERType_ = kJERra2;
  else {cout<<"problem with opt in JER"<<endl; assert(0);} //enforce a complete set of options!

  if ( getOptPiece("METunc",opt)== theMETuncNames_[kMETunc0]) theMETuncType_ = kMETunc0;
  else  if ( getOptPiece("METunc",opt)== theMETuncNames_[kMETuncUp]) theMETuncType_ = kMETuncUp;
  else  if ( getOptPiece("METunc",opt)== theMETuncNames_[kMETuncDown]) theMETuncType_ = kMETuncDown;
  else {cout<<"problem with opt in METunc"<<endl; assert(0) ;} //enforce a complete set of options!

  if ( getOptPiece("PUunc",opt)== thePUuncNames_[kPUunc0]) thePUuncType_ = kPUunc0;
  else  if ( getOptPiece("PUunc",opt)== thePUuncNames_[kPUuncUp]) thePUuncType_ = kPUuncUp;
  else  if ( getOptPiece("PUunc",opt)== thePUuncNames_[kPUuncDown]) thePUuncType_ = kPUuncDown;
  else {cout<<"problem with opt in PU"<<endl; assert(0) ;} //enforce a complete set of options!

  if ( getOptPiece("BTag",opt)== theBTagEffNames_[kBTagEff0]) theBTagEffType_ = kBTagEff0;
  else  if ( getOptPiece("BTag",opt)== theBTagEffNames_[kBTagEffup]) theBTagEffType_ = kBTagEffup;
  else  if ( getOptPiece("BTag",opt)== theBTagEffNames_[kBTagEffdown]) theBTagEffType_ = kBTagEffdown;
  else  if ( getOptPiece("BTag",opt)== theBTagEffNames_[kBTagEff02]) theBTagEffType_ = kBTagEff02;
  else  if ( getOptPiece("BTag",opt)== theBTagEffNames_[kBTagEffup2]) theBTagEffType_ = kBTagEffup2;
  else  if ( getOptPiece("BTag",opt)== theBTagEffNames_[kBTagEffdown2]) theBTagEffType_ = kBTagEffdown2;
  else  if ( getOptPiece("BTag",opt)== theBTagEffNames_[kBTagEff03]) theBTagEffType_ = kBTagEff03;
  else  if ( getOptPiece("BTag",opt)== theBTagEffNames_[kBTagEffup3]) theBTagEffType_ = kBTagEffup3;
  else  if ( getOptPiece("BTag",opt)== theBTagEffNames_[kBTagEffdown3]) theBTagEffType_ = kBTagEffdown3;
  else  if ( getOptPiece("BTag",opt)== theBTagEffNames_[kBTagEff04]) theBTagEffType_ = kBTagEff04;
  else  if ( getOptPiece("BTag",opt)== theBTagEffNames_[kBTagEffup4]) theBTagEffType_ = kBTagEffup4;
  else  if ( getOptPiece("BTag",opt)== theBTagEffNames_[kBTagEffdown4]) theBTagEffType_ = kBTagEffdown4;
  else {cout<<"problem with opt in btag"<<endl; assert(0) ;} //enforce a complete set of options!

  if ( getOptPiece("HLT",opt)== theHLTEffNames_[kHLTEff0]) theHLTEffType_ = kHLTEff0;
  else  if ( getOptPiece("HLT",opt)== theHLTEffNames_[kHLTEffup]) theHLTEffType_ = kHLTEffup;
  else  if ( getOptPiece("HLT",opt)== theHLTEffNames_[kHLTEffdown]) theHLTEffType_ = kHLTEffdown;
  else {cout<<"problem with opt in HLT"<<endl;assert(0) ;} //enforce a complete set of options!

cout<<"Got options: "<<endl
    <<theJESNames_[theJESType_]<<endl
    <<theJERNames_[theJERType_]<<endl
    <<theMETuncNames_[theMETuncType_]<<endl
    <<thePUuncNames_[thePUuncType_]<<endl
    <<theBTagEffNames_[theBTagEffType_]<<endl
    <<theHLTEffNames_[theHLTEffType_]<<endl;

}


double EventCalculator::getWeight(Long64_t nentries) {
  if( isSampleRealData() ) return 1;

  if (theScanType_==kmSugra) return 1;//special weighting in effect
  else  if (theScanType_==kSMS) return 1;//special weighting in effect

  double sigma = getCrossSection();
  double w = lumi_ * sigma / double(nentries);

  if (sampleName_.Contains("QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6")) w *= weight;

  return  w;
}



//1d reweighting
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupMCReweightingUtilities
float EventCalculator::getPUWeight( reweight::LumiReWeighting lumiWeights ) {
  if (isSampleRealData() ) return 1;

  float weight=1;
  int npv = -1;
  for ( unsigned int i = 0; i<PU_bunchCrossing->size() ; i++) {
    //consider only in-time PU
    int BX =PU_bunchCrossing->at(i) ;
    if(BX == 0) { 
      npv =PU_TrueNumInteractions->at(i);
    }
  }
  //in-time PU only
  weight = lumiWeights.ITweight( npv );
  return weight;
}

/* need to figure this out for 2012
//3d reweighting
float EventCalculator::getPUWeight( Lumi3DReWeighting lumiWeights ) {
  if (isSampleRealData() ) return 1;

  
  float weight=1;
  int nm1 = -1; int n0 = -1; int np1 = -1;

  for ( unsigned int i = 0; i<pileupsummaryinfo.size() ; i++) {
    //npv = pileupsummaryinfo.at(i).addpileupinfo_getPU_NumInteractions;
    //sum_nvtx += float(npv);

    int BX = pileupsummaryinfo.at(i).addpileupinfo_getBunchCrossing;

    if(BX == -1) { 
      nm1 = pileupsummaryinfo.at(i).addpileupinfo_getPU_NumInteractions;
    }
    else if(BX == 0) {
      n0 = pileupsummaryinfo.at(i).addpileupinfo_getPU_NumInteractions;
    }
    else if(BX == 1) {
      np1 = pileupsummaryinfo.at(i).addpileupinfo_getPU_NumInteractions;
    }

  }


  //3d reweighting
  weight = lumiWeights.weight3D( nm1,n0,np1);

  //if this is ttbar Fall11 MC 
  if (sampleName_.Contains("TTJets_TuneZ2_7TeV-madgraph-tauola_Fall11_v2") ) {
    int nPV =  countGoodPV();
    //From Kristen (via RA4 group)
    double pvweights[25]={0 , 0.549398 , 0.854301 , 1.0274 , 1.19307 ,
			  1.29101 , 1.33746 , 1.25162 , 1.16767 , 1.12114 ,
			  0.954877 , 0.771361 , 0.583225 , 0.46853 , 0.343897 ,
			  0.238567 , 0.160576 , 0.11326 , 0.0544237 , 0.0236446 ,
			  0 , 0.0850062 , 0 , 0 , 0};
    //this is applied to normalize the nGoodPV distribution
    //of ttbar Fall MC to that of data after applying the following cuts
    //{HT>400,MET>250,njets>=3,bjets>=1}
    //(Kristen is applying the same weighting)
    double corr[25] = {0., 0.8919, 0.8289, 1.1152, 1.0013, 0.8529, 1.0224, 0.9010, 1.0397,
		       0.9782, 1.0765, 1.2540, 1.3285, 0.8122, 0.9539, 1.5913, 1.3473, 1.9760,
		       0., 0., 0., 0., 0., 0., 0.};

    if(nPV > 24) weight =0;
    else weight = pvweights[nPV]*0.975*corr[nPV];
    //else weight = pvweights[nPV];  
  }

  return weight;
}
*/

bool EventCalculator::isGoodMuon(const unsigned int k, const bool disableRelIso, const float ptthreshold) {

  // updated for 2012...basically just copying keith's code

  //sanity check
  if ( k >= pf_mus_pt->size() ) return false;

  if (!isGoodPV(0)) return false;

  if (pf_mus_pt->at(k) < 10) return false;
  if (fabs(pf_mus_eta->at(k)) >= 2.4 ) return false; 
  
  //i like to use TMath::Nint for quantities stored as float that are really integers. just to be safe

  if ( TMath::Nint(pf_mus_id_GlobalMuonPromptTight->at(k)) == 0) return false; 
  // GlobalMuonPromptTight includes: isGlobal, globalTrack()->normalizedChi2() < 10, numberOfValidMuonHits() > 0

  if ( TMath::Nint(pf_mus_numberOfMatchedStations->at(k)) <= 1 ) return false;
  // At least one matched station includes requirement of arbitrated tracker muon, so no need for that explicitly

  const float beamx = beamSpot_x->at(0);   
  const float beamy = beamSpot_y->at(0);   
  float d0 = pf_mus_tk_d0dum->at(k) - beamx*sin(pf_mus_tk_phi->at(k)) + beamy*cos(pf_mus_tk_phi->at(k));
  if ( fabs(d0) >= 0.2 ) return false;
  if ( fabs(pf_mus_tk_vz->at(k) - pv_z->at(0) ) >= 0.5 ) return false;
  if ( TMath::Nint(pf_mus_tk_numvalPixelhits->at(k)) == 0 ) return false;

  if ( TMath::Nint(pf_mus_tk_LayersWithMeasurement->at(k)) <= 5 ) return false;
  // now isolation with delta beta corrections are recommended with DelR cone of 0.4
  float isoNeutral = pf_mus_pfIsolationR04_sumNeutralHadronEt->at(k) + pf_mus_pfIsolationR04_sumPhotonEt->at(k) - 0.5*pf_mus_pfIsolationR04_sumPUPt->at(k);
  isoNeutral = ( isoNeutral > 0) ? isoNeutral : 0;
  float  muRelIso = (pf_mus_pfIsolationR04_sumChargedHadronPt->at(k) + isoNeutral) / pf_mus_pt->at(k);
  if (muRelIso >= 0.2) return false; 

  //if we haven't failed anything yet...then it's good
  return true;
}

bool EventCalculator::isGoodRecoMuon(const unsigned int imuon, const bool disableRelIso, const float ptthreshold) {

  assert(disableRelIso); //the iso calculation below is completely worthless. it uses PF quantities on RECO muons.

  return true; //FIXME CFA
/*
  if (myMuonsRECO->at(imuon).pt >= ptthreshold
      && fabs(myMuonsRECO->at(imuon).eta)<2.4
      && myMuonsRECO->at(imuon).GlobalMuonPromptTight == 1
      && myMuonsRECO->at(imuon).isTrackerMuon == 1 
      && myMuonsRECO->at(imuon).innerTrack_numberOfValidHits >=11
      && myMuonsRECO->at(imuon).track_hitPattern_numberOfValidPixelHits >= 1
      && fabs(myMuonsRECOhelper->at(imuon).dxywrtBeamSpot) < 0.02
      && fabs(myMuonsRECO->at(imuon).vz - myVertex->at(0).z ) <1
      && ((myMuonsRECO->at(imuon).chargedHadronIso 
	   + myMuonsRECO->at(imuon).photonIso 
	   + myMuonsRECO->at(imuon).neutralHadronIso)/myMuonsRECO->at(imuon).pt <0.2 ||disableRelIso)
      ) {
    return true;
  }
  */
  return false;
}

bool EventCalculator::isCleanMuon(const unsigned int imuon, const float ptthreshold) {

  if (!isGoodMuon(imuon,false,ptthreshold)) return false;

  //clean muons if using reco-pfjets
  if(theJetType_ == kRECOPF) {    
    bool isNearJet = false;
    for ( unsigned int j = 0; j< jets_AK5PF_pt->size(); j++) {
      if( isGoodJet(j) && jmt::deltaR( jets_AK5PF_eta->at(j), jets_AK5PF_phi->at(j),pf_mus_eta->at(imuon), pf_mus_phi->at(imuon))<0.3 ) {
	isNearJet = true ; break;
      }
    }
    if(isNearJet) return false;
  }

  return true;

}

//this bool could be turned into a more powerful selector for N-1 studies. keep it simple for now
bool EventCalculator::isGoodElectron(const unsigned int k, const bool disableRelIso, const float ptthreshold) {

  //sanity check
  if ( k >= pf_els_pt->size() ) return false;

  if (!isGoodPV(0)) return false;

  const float  beamx = beamSpot_x->at(0);
  const float  beamy = beamSpot_y->at(0); 

  // new electron selection for 2012. Use "cut based" Veto selection:
  /* selection now split into endcap and barrel with different cuts for each.
     barrel  endcap
     defined as:     pf_els_isEB     pf_els_isEE  // no explicit veto on gap any more?
     dEtaIn	      0.007   0.01
     dPhiIn	      0.8     0.7
     sigmaiEtaiEta   0.01    0.03
     H/E	      0.15    --
     d0	      0.04    0.04
     dZ	      0.2     0.2
     PF isolation    0.15    0.15   // now done with rho corrections
  */

  int ndone=0;

  if ( TMath::Nint(pf_els_isEB->at(k)) == 1) {
    ++ndone;
    if (pf_els_dEtaIn->at(k) > 0.007)  return false;
    if (pf_els_dPhiIn->at(k) > 0.8)  return false;
    if (pf_els_sigmaIEtaIEta->at(k) > 0.01) return false;
    if (pf_els_hadOverEm->at(k) > 0.15) return false;
  }
  if (TMath::Nint(pf_els_isEE->at(k)) == 1) {
    ++ndone;
    if (pf_els_dEtaIn->at(k) > 0.01)  return false;
    if (pf_els_dPhiIn->at(k) > 0.7)  return false;
    if (pf_els_sigmaIEtaIEta->at(k) > 0.03) return false;
  }
  
  if (ndone != 1) { //sanity
    cout<<"[isGoodElectron] I went through "<<ndone<<" of the barrel+endcap checks. Weird!"<<endl;
  }

  float d0 = pf_els_d0dum->at(k) - beamx*sin(pf_els_phi->at(k)) + beamy*cos(pf_els_phi->at(k));
  if ( fabs(d0) >= 0.04 ) return false; 
  if ( fabs(pf_els_vz->at(k) - pv_z->at(0) ) >= 0.2 ) return false; 
  
  // keep same pt and eta cuts (except no explicit exclusion of gap region)
  if (pf_els_pt->at(k) < 10) return false;
  if (fabs(pf_els_scEta->at(k)) >= 2.5 ) return false; 
  
  // now recommended isolation is PF relative isolation with cone = 0.3 with rho (effective area) corrections for PU
    
  // isocorr = PFChargedIso (PFNoPU) + max(PFIso(&gamma;+NH) - rho * Aeff(&gamma;+NH), 0.)               
  float rho = rho_kt6PFJetsForIsolation2011;
  // get effective area from delR=0.3 2011 data table for neutral+gamma based on supercluster eta pf_els_scEta->at(k)
  float AE = 0.10; 
  //jmt -- this strikes me as inefficient, but it's not a big enough deal to change it
  if ( fabs( pf_els_scEta->at(k) ) > 1.0 ) AE = 0.12;      
  if ( fabs( pf_els_scEta->at(k) ) > 1.479 ) AE = 0.085;      
  if ( fabs( pf_els_scEta->at(k) ) > 2.0 ) AE = 0.11;      
  if ( fabs( pf_els_scEta->at(k) ) > 2.2 ) AE = 0.12;      
  if ( fabs( pf_els_scEta->at(k) ) > 2.3 ) AE = 0.12;
  if ( fabs( pf_els_scEta->at(k) ) > 2.4 ) AE = 0.13;      
  float eleIso = pf_els_PFphotonIsoR03->at(k) + pf_els_PFneutralHadronIsoR03->at(k) - rho*AE;
  float  elRelIso = ( pf_els_PFchargedHadronIsoR03->at(k) + ( eleIso > 0 ? eleIso : 0.0 ) )/pf_els_pt->at(k);
  if (elRelIso >= 0.15) return false; 
  
  return true;
}


unsigned int EventCalculator::countEle(const float ptthreshold) {

  unsigned int ngoodele=0;
  for (unsigned int i=0; i <pf_els_pt->size() ; i++) {
    if(isGoodElectron(i,false,ptthreshold))      ++ngoodele;
  }
  return ngoodele;
}

unsigned int EventCalculator::countMu(const float ptthreshold) {

  unsigned int ngoodmu=0;
  unsigned int nmu = pf_mus_pt->size();
  for ( unsigned int i = 0; i< nmu; i++) {
    if (isCleanMuon(i,ptthreshold)) {
      //once we reach here we've got a good muon in hand
      ++ngoodmu;   
    }
  }
  
  return ngoodmu;
}

//eventually we might want to get more info, but for now this is good enough
float EventCalculator::getBestZCandidate(const float pt_threshold1,const float pt_threshold2) {

  float mee = getBestZeeCandidate(pt_threshold1,pt_threshold2);
  float mmm = getBestZmmCandidate(pt_threshold1,pt_threshold2);

  float best = ( fabs(mee-mZ_) < fabs(mmm-mZ_)) ? mee : mmm;

  return best;

}

float EventCalculator::getBestZeeCandidate(const float pt_threshold1,const float pt_threshold2) {
  //code copied from Zmm code below (ugly solution)
  //cout<<"Zee"<<endl;
  assert(pt_threshold1>=pt_threshold2);

  //my comment -- it pains me that the structure of this code needs to be duplicated for ee and mm, but i don't feel like being more clever
  float best_mll = -99;
  //give me lists of all the indices of the good, high pt muons
  std::vector<uint> goodeles1;
  std::vector<uint> goodeles2;
  unsigned int nmu = pf_els_pt->size();
  for ( unsigned int i = 0; i< nmu; i++) {
    if (isGoodElectron(i,pt_threshold1)) goodeles1.push_back( i );
    if (isGoodElectron(i,pt_threshold2)) goodeles2.push_back( i );
  }

  //need to have at least 2 leptons, and 1 of them must be above the higher pt threshold
  if( goodeles1.size() >=1 && goodeles2.size()>=2 ) {
    //loop through the lepton pairs, and check if they are
    //opposite-signed and form an invariant mass within the Z-mass
    for ( unsigned int i = 0; i< goodeles1.size(); i++) {
      for ( unsigned int j = 0; j< goodeles2.size(); j++) {
	
	//check opposite charge
	if (TMath::Nint(pf_els_charge->at( goodeles1.at(i) )) != TMath::Nint(pf_els_charge->at( goodeles2.at(j) )) ) {
	  
	  //check invariant mass
	  double z_en = pf_els_energy->at(goodeles1.at(i)) + pf_els_energy->at(goodeles2.at(j));
	  double z_px = pf_els_pt->at(goodeles1.at(i))*cos( pf_els_phi->at(goodeles1.at(i))) + pf_els_pt->at(goodeles2.at(j))*cos( pf_els_phi->at(goodeles2.at(j)));
	  double z_py = pf_els_pt->at(goodeles1.at(i))*sin( pf_els_phi->at(goodeles1.at(i))) + pf_els_pt->at(goodeles2.at(j))*sin( pf_els_phi->at(goodeles2.at(j)));
	  double z_pz = pf_els_pt->at(goodeles1.at(i))*sinh(pf_els_eta->at(goodeles1.at(i))) + pf_els_pt->at(goodeles2.at(j))*sinh(pf_els_eta->at(goodeles2.at(j)));
	  float thismll = sqrt(z_en*z_en - z_px*z_px - z_py*z_py - z_pz*z_pz);
	  
	  if ( fabs(thismll - mZ_) < fabs(best_mll-mZ_)) best_mll=thismll;
	}

      }// j loop
    }//i loop

  } //two high pt leptons
  //cout<<"done"<<endl;
  return best_mll;
}

float EventCalculator::getBestZmmCandidate(const float pt_threshold1,const float pt_threshold2) {
  //code lifted from Don and converted for cfA
  //also adapted to work with two different pt thresholds (i hope)

  //cout<<"Zmm"<<endl;

  assert(pt_threshold1>=pt_threshold2);

  //my comment -- it pains me that the structure of this code needs to be duplicated for ee and mm, but i don't feel like being more clever
  float best_mll = -99;
  //give me lists of all the indices of the good, high pt muons
  std::vector<uint> goodmuons1;
  std::vector<uint> goodmuons2;
  unsigned int nmu = pf_mus_pt->size();
  for ( unsigned int i = 0; i< nmu; i++) {
    if (isCleanMuon(i,pt_threshold1)) goodmuons1.push_back( i );
    if (isCleanMuon(i,pt_threshold2)) goodmuons2.push_back( i );
  }

  //need to have at least 2 leptons, and 1 of them must be above the higher pt threshold
  if( goodmuons1.size() >=1 && goodmuons2.size()>=2 ) {
    //loop through the lepton pairs, and check if they are
    //opposite-signed and form an invariant mass within the Z-mass
    for ( unsigned int i = 0; i< goodmuons1.size(); i++) {
      for ( unsigned int j = 0; j< goodmuons2.size(); j++) {
	
	//check opposite charge
	if (TMath::Nint(pf_mus_charge->at( goodmuons1.at(i) )) != TMath::Nint(pf_mus_charge->at( goodmuons2.at(j) )) ) {
	  
	  //check invariant mass
	  double z_en = pf_mus_energy->at(goodmuons1.at(i)) + pf_mus_energy->at(goodmuons2.at(j));
	  double z_px = pf_mus_pt->at(goodmuons1.at(i))*cos( pf_mus_phi->at(goodmuons1.at(i))) + pf_mus_pt->at(goodmuons2.at(j))*cos( pf_mus_phi->at(goodmuons2.at(j)));
	  double z_py = pf_mus_pt->at(goodmuons1.at(i))*sin( pf_mus_phi->at(goodmuons1.at(i))) + pf_mus_pt->at(goodmuons2.at(j))*sin( pf_mus_phi->at(goodmuons2.at(j)));
	  double z_pz = pf_mus_pt->at(goodmuons1.at(i))*sinh(pf_mus_eta->at(goodmuons1.at(i))) + pf_mus_pt->at(goodmuons2.at(j))*sinh(pf_mus_eta->at(goodmuons2.at(j)));
	  float thismll = sqrt(z_en*z_en - z_px*z_px - z_py*z_py - z_pz*z_pz);
	  
	  if ( fabs(thismll - mZ_) < fabs(best_mll-mZ_)) best_mll=thismll;
	}

      }// j loop
    }//i loop

  } //two high pt leptons
  //cout<<"done"<<endl;
  return best_mll;
}


bool EventCalculator::isGoodTau(const unsigned int itau, const float pTthreshold, const float etaMax) {
  return false; //FIXME CFA
/*
  if ( getTauPt(itau) >= pTthreshold
      && fabs(taus_eta->at(itau)) < etaMax
      && taus_againstElectron->at(itau) > 0
      && taus_againstMuon->at(itau) > 0
      && taus_taNC_half->at(itau) > 0
      ) {
    return true;
  }

  return false;
*/
}

float EventCalculator::getTauPt( unsigned int itau) {
  return 0; //FIXME CFA
  //  return taus_pt->at(itau);
}

unsigned int EventCalculator::countTau() {
  return 0; //FIXME CFA
  //  if (taus_pt==0) return 0; //FIXME CFA
/*
  unsigned int ngoodtau=0;
  for (unsigned int i=0; i <taus_pt->size() ; i++) {
    if(isGoodTau(i))      ++ngoodtau;
  }
  return ngoodtau;
*/
}

TString EventCalculator::stripTriggerVersion(const TString & fullname, int & version) {

  //example input: HLT_HT400_AlphaT0p52_v10
  TRegexp rv("_v[0-9]*$"); //i suppose we could save CPU by making this static or something.....

  Ssiz_t cutoff=  fullname.Index(rv);
  TString stripped=fullname(0,cutoff);
  version =  TString(fullname(cutoff+2,fullname.Length())).Atoi();
  //output: HLT_HT400_AlphaT0p52 , version = 10
  return stripped;
}

bool EventCalculator::passHLT(map<TString,triggerData> & triggerlist, bool alwaysPassMc) {

  if ( !isSampleRealData() && alwaysPassMc) {
    for ( map<TString,triggerData>::iterator ii=triggerlist.begin(); ii!=triggerlist.end(); ++ii) ii->second.pass = true;
    return true;
  }

  bool passanything = false;
  int version;
  //cout<<" --- "<<endl;
  //loop over triggers in ntuple
  for (unsigned int itrig=0; itrig<trigger_name->size(); itrig++) {
    TString thistrigger = TString(trigger_name->at(itrig).c_str());
    if (thistrigger(0,4) != "HLT_") continue; //some triggers don't start with HLT_ !!!
    TString strippedname=stripTriggerVersion(thistrigger,version); //removes the _vXX on the end
    strippedname= strippedname(4,strippedname.Length()); //get rid of HLT_ on the front
    map<TString,triggerData>::iterator ittrig = triggerlist.find(strippedname);
    if ( ittrig != triggerlist.end() ) {//found it
      ittrig->second.pass = ( TMath::Nint(trigger_decision->at(itrig))==1 );
      ittrig->second.prescale =  TMath::Nint(trigger_prescalevalue->at(itrig));
      ittrig->second.version = version;
      passanything = passanything || ittrig->second.pass;
      // cout<<"\t"<<thistrigger<<" "<<ittrig->second.pass<<"("<<ittrig->second.prescale<<" "<<ittrig->second.version <<")"<<" "<<passanything<<endl;
    }

  }
  return passanything;
}




bool EventCalculator::passUtilityPrescaleModuleHLT() {

   if( !isSampleRealData()) return true; 

   return ( passprescalePFHT350filter_decision>0);

}


void EventCalculator::loadHLTMHTeff() {
  if (f_eff_mht_==0) {

    //there are four curves
    //not sure what the best way to combine them are
    f_eff_mht_ = new TFile("mht_eff_uptojul6.root","READ");
    //mhtgraph  = (TGraphAsymmErrors *) f_eff_mht->Get("eff_MHT60");
    //mhtgraph  = (TGraphAsymmErrors *) f_eff_mht->Get("eff_MHT70");
    //mhtgraph  = (TGraphAsymmErrors *) f_eff_mht->Get("eff_MHT75");
    mhtgraph_  = (TGraphAsymmErrors *) f_eff_mht_->Get("eff_MHT80");

    //fill the graphs with plus and minus variations
    mhtgraphPlus_ = new TGraphAsymmErrors(mhtgraph_->GetN());
    mhtgraphMinus_ = new TGraphAsymmErrors(mhtgraph_->GetN());
    for (int i=0; i<mhtgraph_->GetN(); i++) {
      double x,y;
      mhtgraph_->GetPoint(i,x,y);
      double exl= mhtgraph_->GetErrorXlow(i);
      double exh= mhtgraph_->GetErrorXhigh(i);
      double eyl= mhtgraph_->GetErrorYlow(i);
      double eyh= mhtgraph_->GetErrorYhigh(i);
      //shift y
      double yup = y+eyh > 1? 1: y+eyh;
      double ydown = y-eyl < 0 ? 0: y-eyl;
      mhtgraphPlus_->SetPoint(i,x,yup);
      mhtgraphMinus_->SetPoint(i,x,ydown);
      mhtgraphPlus_->SetPointError(i,exl,exh,eyl,eyh); //the errors don't matter
      mhtgraphMinus_->SetPointError(i,exl,exh,eyl,eyh);
    }
  }
}

float EventCalculator::getHLTMHTeff(float offMET, float offHT, uint nElectrons, uint nMuons, double mindphin) {

  float eff=1;

  bool isSingleE = (nElectrons==1 && nMuons==0);
  bool isSingleMu = (nElectrons==0 && nMuons==1);
  bool is0L = (nElectrons==0 && nMuons==0);

  //TGraphAsymmErrors * gr=mhtgraph_;
  //if ( theHLTEffType_ ==kHLTEffup) gr=htgraphPlus;
  //else if ( theHLTEffType_ ==kHLTEffdown) gr=htgraphMinus;

  //for prelim summer result
  //eff = gr->Eval(offMET);

  //for 2011 full result
  //updated with pixelLumiCalc numbers
  if(offHT>400 && offMET>=150 && offMET<250){
    if( is0L && mindphin > 4 ){
      eff =  0.850;
    }
    else if( isSingleE && mindphin > 4){
      eff = 0.955;
    }
    else if( isSingleMu && mindphin > 4){
      eff = 0.990;
    }
    else if ( is0L && mindphin <= 4){
      eff = 0.912;
    }
  }
  //if(offHT>400 && offMET>=200 && offMET<250){
  //  //    eff = 0.94;
  //  //eff = 0.859362; //after HT cut
  //
  //  //errors are wrong!
  //
  //  ////assign a 3% error in the SB region
  //  //if ( theHLTEffType_ ==kHLTEffup){
  //  //  eff = eff*1.03; assert(0);
  //  //}
  //  //else if ( theHLTEffType_ ==kHLTEffdown){
  //  //  eff = eff*0.97; assert(0);
  //  //}
  //
  //  //after HT>400, ==0b
  //  if( is0L && mindphin > 4 ){
  //    eff = 0.841;
  //  }
  //  else if( isSingleE && mindphin > 4){
  //    eff = 0.991;
  //  }
  //  else if( isSingleMu && mindphin > 4){
  //    eff = 1.0;
  //  }
  //  else if ( is0L && mindphin <= 4){
  //    eff = 0.936;
  //  }
  //}
  else if(offHT>400 && offMET>250){
    //    eff = 0.998;
    //eff = 0.975299; //after HT cut

    ////assign a +0.1%(-0.5%) error in the SIG region
    //if ( theHLTEffType_ ==kHLTEffup){
    //  eff = eff*1.001; assert(0);
    //}
    //else if ( theHLTEffType_ ==kHLTEffdown){
    //  eff = eff*0.995; assert(0);
    //}

    //after HT>400, ==0b
    //updated with pixelLumiCalc numbers
    if( is0L && mindphin > 4 ){
      eff = 0.981;
    }
    else if( isSingleE && mindphin > 4){
      eff = 1.0;
    }
    else if( isSingleMu && mindphin > 4){
      eff = 0.999;
    }

    else if ( is0L && mindphin <= 4){
      eff = 0.981; // same as nominal for now
    }


  }


  return eff;
}

void EventCalculator::loadHLTHTeff() {
  if (f_eff_ht300_==0 && f_eff_ht350_==0) {

    /*
there is a horrible hack here where i have manually set the error on the trigger eff to 1.5% at high HT
    */

  //get this file from ~/public/wteo/
    f_eff_ht300_ = new TFile("ht300_out_166864_eff_histogram.root","READ");
    f_eff_ht350_ = new TFile("ht350_out_178866_eff_histogram.root","READ");

    htgraph_ht300_  = (TGraphAsymmErrors *) f_eff_ht300_->Get("myeff");
    htgraph_ht350_  = (TGraphAsymmErrors *) f_eff_ht350_->Get("ht350eff");

    //fill the graphs with plus and minus variations
    htgraphPlus_ht300_ = new TGraphAsymmErrors(htgraph_ht300_->GetN());
    htgraphMinus_ht300_ = new TGraphAsymmErrors(htgraph_ht300_->GetN());
    htgraphPlus_ht350_ = new TGraphAsymmErrors(htgraph_ht350_->GetN());
    htgraphMinus_ht350_ = new TGraphAsymmErrors(htgraph_ht350_->GetN());
    for (int i=0; i<htgraph_ht300_->GetN(); i++) {
      double x,y;
      htgraph_ht300_->GetPoint(i,x,y);
      double exl= htgraph_ht300_->GetErrorXlow(i);
      double exh= htgraph_ht300_->GetErrorXhigh(i);
      double eyl= htgraph_ht300_->GetErrorYlow(i);
      double eyh= htgraph_ht300_->GetErrorYhigh(i);
      //shift y
      double yup = y+eyh > 1? 1: y+eyh;
      double ydown = y-eyl < 0 ? 0: y-eyl;
      //if (x>500) ydown=0.985; //Horrible hack!
      htgraphPlus_ht300_->SetPoint(i,x,yup);
      htgraphMinus_ht300_->SetPoint(i,x,ydown);
      htgraphPlus_ht300_->SetPointError(i,exl,exh,eyl,eyh); //the errors don't matter
      htgraphMinus_ht300_->SetPointError(i,exl,exh,eyl,eyh);

      double x2,y2;
      htgraph_ht350_->GetPoint(i,x2,y2);
      double ex2l= htgraph_ht350_->GetErrorXlow(i);
      double ex2h= htgraph_ht350_->GetErrorXhigh(i);
      double ey2l= htgraph_ht350_->GetErrorYlow(i);
      double ey2h= htgraph_ht350_->GetErrorYhigh(i);
      //shift y
      double y2up = y2+ey2h > 1? 1: y2+ey2h;
      double y2down = y2-ey2l < 0 ? 0: y2-ey2l;
      //if (x>500) ydown=0.985; //Horrible hack!
      htgraphPlus_ht350_->SetPoint(i,x2,y2up);
      htgraphMinus_ht350_->SetPoint(i,x2,y2down);
      htgraphPlus_ht350_->SetPointError(i,ex2l,ex2h,ey2l,ey2h); //the errors don't matter
      htgraphMinus_ht350_->SetPointError(i,ex2l,ex2h,ey2l,ey2h);


    }
  }

}

float EventCalculator::getHLTHTeff(float offHT) {

  float eff=100;

  TGraphAsymmErrors * gr1=htgraph_ht300_;
  TGraphAsymmErrors * gr2=htgraph_ht350_;
  if ( theHLTEffType_ ==kHLTEffup){
    gr1=htgraphPlus_ht300_;
    gr2=htgraphPlus_ht350_;
  }
  else if ( theHLTEffType_ ==kHLTEffdown){
    gr1=htgraphMinus_ht300_;
    gr2=htgraphMinus_ht350_;
  }

  //2011 dataset - weighted average of HT efficiencies
  double frac_HT2XX = 0.071;
  double frac_HT300 = 0.462;
  double frac_HT350 = 0.467;
  eff = frac_HT2XX*1.0 + frac_HT300*gr1->Eval(offHT) + frac_HT350*gr2->Eval(offHT);

  return eff;
}

void EventCalculator::loadDataJetRes(){
  if(fDataJetRes_ ==0){
    fDataJetRes_ = new TFile("additionalTools/datajetres.root", "READ");
    if (fDataJetRes_->IsZombie() ) cout<<" ERROR loading additionalTools/datajetres.root"<<endl;
    else cout<<"Loaded Data Jet Resolution file"<<endl;
    hEta0_0p5_ = (TH1D*) fDataJetRes_->Get("hEta0_0p5");
    hEta0p5_1_ = (TH1D*) fDataJetRes_->Get("hEta0p5_1");
    hEta1_1p5_ = (TH1D*) fDataJetRes_->Get("hEta1_1p5");
    hEta1p5_2_ = (TH1D*) fDataJetRes_->Get("hEta1p5_2");
    hEta2_2p5_ = (TH1D*) fDataJetRes_->Get("hEta2_2p5");
    hEta2p5_3_ = (TH1D*) fDataJetRes_->Get("hEta2p5_3");
    hEta3_5_   = (TH1D*) fDataJetRes_->Get("hEta3_5");

    cout<<    hEta0_0p5_ <<endl
	<<    hEta0p5_1_ <<endl
    << hEta1_1p5_  <<endl
	 <<hEta1p5_2_ <<endl 
	 <<hEta2_2p5_ <<endl 
	 <<hEta2p5_3_ <<endl 
	<<hEta3_5_   <<endl ;

  }
}

float EventCalculator::getDataJetRes(float pt, float eta){
  float abseta = fabs(eta);
  if(abseta<=0.5)                    return hEta0_0p5_->GetBinContent(hEta0_0p5_->FindFixBin(pt));
  else if(abseta>0.5 && abseta<=1.0) return hEta0p5_1_->GetBinContent(hEta0p5_1_->FindFixBin(pt));
  else if(abseta>1.0 && abseta<=1.5) return hEta1_1p5_->GetBinContent(hEta1_1p5_->FindFixBin(pt));
  else if(abseta>1.5 && abseta<=2.0) return hEta1p5_2_->GetBinContent(hEta1p5_2_->FindFixBin(pt));
  else if(abseta>2.0 && abseta<=2.5) return hEta2_2p5_->GetBinContent(hEta2_2p5_->FindFixBin(pt));
  else if(abseta>2.5 && abseta<=3.0) return hEta2p5_3_->GetBinContent(hEta2p5_3_->FindFixBin(pt));
  else if(abseta>3.0 && abseta<=5.0) return hEta3_5_  ->GetBinContent(hEta3_5_  ->FindFixBin(pt));
  else {assert(0);}
}

bool EventCalculator::isGoodPV(unsigned int ipv) {
  
  //sanity check
  if ( ipv >= pv_z->size()) return false;

  bool isgood=false;
  
  if ( ! pv_isFake->at(ipv) ) {
    if ( fabs(pv_z->at(ipv)) < 24 ) {
      if (sqrt(pow(pv_x->at(ipv),2.0)+pow(pv_y->at(ipv),2.0))  < 2 ) {
	if ( pv_ndof->at(ipv) > 4 ) {
	  isgood=true;
	}
      }
    }
  }
  return isgood;
}

int EventCalculator::countGoodPV() {
  
  int ngood=0;

  for (unsigned int ipv = 0; ipv<pv_x->size(); ipv++) {
    if (isGoodPV(ipv)) ++ngood;
  }
  return ngood; 
}

bool EventCalculator::passPV() {
  //consider only the first PV
  return isGoodPV(0);

}

float EventCalculator::getHT(float ptthreshold) {
  float ht=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
    if (isGoodJet( i,ptthreshold ) ) ht+= getJetPt(i);
  }
  return ht;
}

float EventCalculator::getST(float jetthreshold,float leptonthreshold) {
  float st=getHT(jetthreshold);

  for (unsigned int i=0; i <pf_els_pt->size() ; i++) {
    if (isGoodElectron(i,false,leptonthreshold)) st += pf_els_pt->at(i);
  }

  for ( unsigned int i = 0; i< pf_mus_pt->size(); i++) {
    if (isCleanMuon(i,leptonthreshold)) st += pf_mus_pt->at(i);
  }
  
  return st;
}


float EventCalculator::getMET() {
  float myMET=-1;
  float myMETphi =-1000;

  if(theMETType_== kPFMET){ 
    getUncorrectedMET(myMET,myMETphi);

  }
  else if(theMETType_ == kPFMETTypeI){
    getCorrectedMET(myMET,myMETphi);
  }

  //JER and JES are automatically handled inside getCorrectedMET()
  if (theMETuncType_!=kMETunc0) {
    getSmearedUnclusteredMET(myMET,myMETphi);
  }
  return myMET;
}

float EventCalculator::getMETphi() {
  float myMET=-1;
  float myMETphi=-99;

  if(theMETType_== kPFMET){
    getUncorrectedMET(myMET,myMETphi);
  }
  else if(theMETType_ == kPFMETTypeI){
    getCorrectedMET(myMET,myMETphi);
  }

  if (theMETuncType_!=kMETunc0) {
    getSmearedUnclusteredMET(myMET,myMETphi);
  }

  return myMETphi;
}


float EventCalculator::getMHT() {
  std::pair<float,float> mht=getJERAdjustedMHTxy();
  
  return sqrt(mht.first*mht.first + mht.second*mht.second);
}

float EventCalculator::getMHTphi() {
  std::pair<float,float> mht=getJERAdjustedMHTxy();

  return atan2(mht.second,mht.first);
}

float EventCalculator::getMHTphi(int ignoredJet) {
  std::pair<float,float> mht=getJERAdjustedMHTxy(ignoredJet);

  return atan2(mht.second,mht.first);
}



float EventCalculator::getDeltaPhiStar(int & badjet) {

  float dpstar=1e9;
  badjet=-1;

  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    if ( isGoodJetMHT(i) ) {
      float newmhtphi = getMHTphi(i);
      float thisdp = getDeltaPhi( jets_AK5PF_phi->at(i), newmhtphi);
      if ( thisdp < dpstar) {
	dpstar = thisdp;
	badjet = i;
      }
    }
  }
  return dpstar;
}

float EventCalculator::getDeltaPhiMETJetMaxMis(float jetpt = 50){
  if( isSampleRealData() ) return -1;

  const unsigned int BOGUSNUM=999999;
  unsigned int myjet = BOGUSNUM;
  float dp = 99;
  float mymax = 0;  
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    if (isGoodJet(i,jetpt)) {
      float genpt = jets_AK5PF_gen_pt->at(i);
      float pt = getJetPt(i);
      float mis = fabs(genpt - pt);
      if(mis>=mymax){
	mymax = mis;
	myjet = i;
      }
    }    
  }
  if (myjet!=BOGUSNUM) dp =  getDeltaPhi( jets_AK5PF_phi->at(myjet) , getMETphi());
  return dp;
}


float EventCalculator::getMaxJetMis(unsigned int rank=1, unsigned int maxjets=3, float jetpt=50) {
  if( isSampleRealData() ) return -1;
  
  float mymax = 0;  //highest
  float mymax2 = 0; //2nd highest
  unsigned int ngood=0;
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    if (isGoodJet(i,jetpt)) {
      ++ngood;
      float genpt = jets_AK5PF_gen_pt->at(i);
      float pt = getJetPt(i);
      float mis = fabs(genpt - pt);
      if(mis>=mymax){
	mymax2 = mymax;
	mymax = mis;
      }
      else if(mis>mymax2){
	mymax2 = mis;
      }
      if (ngood >= maxjets) break;
    }
  }
  if(rank==1) return mymax;
  else if(rank==2) return mymax2;
  else { assert(0); }
}


float EventCalculator::getMaxJetFracMis(unsigned int rank=1, unsigned int maxjets=3, float jetpt=50){
  if( isSampleRealData() ) return -1;

  float mymax = 0;  //highest
  float mymax2 = 0; //2nd highest
  unsigned int ngood=0;
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    if (isGoodJet(i,jetpt)) {
      ++ngood;
      float genpt = jets_AK5PF_gen_pt->at(i);
      float pt = getJetPt(i);
      float mis = fabs(1.0-pt/genpt);
      if(mis>=mymax){
	mymax2 = mymax;
	mymax = mis;
      }
      else if(mis>mymax2){
	mymax2 = mis;
      }
      if (ngood >= maxjets) break;
    }
  }
  if(rank==1) return mymax;
  else if(rank==2) return mymax2;
  else { assert(0); }
}

double EventCalculator::getMinDeltaPhiMET(unsigned int maxjets) {

  double mindp=99;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    
    if (isGoodJet(i)) {
      ++ngood;
      double dp =  getDeltaPhi( jets_AK5PF_phi->at(i) , getMETphi());
      if (dp<mindp) mindp=dp;
      if (ngood >= maxjets) break;
    }
  }
  
  return mindp;
}


double EventCalculator::getTransverseMETError(unsigned int thisJet) {

  if(!(thisJet<jets_AK5PF_pt->size())) return -99;

  double deltaT=0;

  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    
    if (isGoodJet30(i) && i!=thisJet) {
      double dp = getDeltaPhi( jets_AK5PF_phi->at(i) , jets_AK5PF_phi->at(thisJet) );
      deltaT+=pow(getJetPt(i)*sin(dp),2);
    }
  }
  return 0.1*sqrt(deltaT);
}

double EventCalculator::getDeltaPhiNMET(unsigned int thisJet) {//Luke

  if(!(thisJet<jets_AK5PF_pt->size())) return 1E12;

  double dp =  getDeltaPhi( jets_AK5PF_phi->at(thisJet) , getMETphi() );
  double deltaT = getTransverseMETError(thisJet);
  return dp/atan2(deltaT,getMET());

}

unsigned int EventCalculator::getNthGoodJet(unsigned int goodJetN, float mainpt, float maineta, bool mainid) {
  //return nth good jet starting at 0, or 999999 if it doesn't exist
 
  unsigned int ijet = 999999;
  unsigned int goodJetI=0;
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    if (isGoodJet(i, mainpt, maineta, mainid)) {
      if(goodJetI == goodJetN){
        ijet = i;
        break;
      }
      goodJetI++;
    }
  }
  return ijet;
}

double EventCalculator::getDeltaPhiMETN_deltaT(unsigned int ijet, float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith) {

  if(ijet==999999) return -99;

  double METx = getMET() * cos(getMETphi());
  double METy = getMET() * sin(getMETphi());
  
  //get sum for deltaT
  double sum = 0;
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    if(i==ijet) continue;
    if(isGoodJet(i, otherpt, othereta, otherid)){
      double jetres = dataJetRes ? getDataJetRes(getJetPt(i), jets_AK5PF_eta->at(i)) : 0.1;
      if(keith)  sum += pow( jetres*(METx*getJetPy(i) - METy*getJetPx(i)), 2);
      else sum += pow( jetres*(getJetPx(ijet)*getJetPy(i) - getJetPy(ijet)*getJetPx(i)), 2);
    }//is good jet
  }//i
  
  double deltaT = keith ? sqrt(sum)/getMET() : sqrt(sum)/getJetPt(ijet);
  return deltaT;
}


double EventCalculator::getDeltaPhiMETN_electron_deltaT(unsigned int ielectron, float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith) {

  double METx = getMET() * cos(getMETphi());
  double METy = getMET() * sin(getMETphi());
  
  const double px = pf_els_pt->at(ielectron) * cos(pf_els_phi->at(ielectron));
  const double py = pf_els_pt->at(ielectron) * sin(pf_els_phi->at(ielectron));

  //get sum for deltaT
  double sum = 0;
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    if(isGoodJet(i, otherpt, othereta, otherid)){
      double jetres = dataJetRes ? getDataJetRes(getJetPt(i), jets_AK5PF_eta->at(i)) : 0.1;
      if(keith)  sum += pow( jetres*(METx*getJetPy(i) - METy*getJetPx(i)), 2);
      else sum += pow( jetres*( px*getJetPy(i) - py*getJetPx(i)), 2);
    }//is good jet
  }//i
  
  double deltaT = keith ? sqrt(sum)/getMET() : sqrt(sum)/ (pf_els_pt->at(ielectron));
  return deltaT;
}

double EventCalculator::getDeltaPhiMETN_muon_deltaT(unsigned int imuon, float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith) {

  double METx = getMET() * cos(getMETphi());
  double METy = getMET() * sin(getMETphi());
  
  const double px = pf_mus_pt->at(imuon) * cos(pf_mus_phi->at(imuon));
  const double py = pf_mus_pt->at(imuon) * sin(pf_mus_phi->at(imuon));

  //get sum for deltaT
  double sum = 0;
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    if(isGoodJet(i, otherpt, othereta, otherid)){
      double jetres = dataJetRes ? getDataJetRes(getJetPt(i), jets_AK5PF_eta->at(i)) : 0.1;
      if(keith)  sum += pow( jetres*(METx*getJetPy(i) - METy*getJetPx(i)), 2);
      else sum += pow( jetres*( px*getJetPy(i) - py*getJetPx(i)), 2);
    }//is good jet
  }//i
  
  double deltaT = keith ? sqrt(sum)/getMET() : sqrt(sum)/ (pf_mus_pt->at(imuon));
  return deltaT;
}

//the ielectron is a different logic than is used for the jet index in getDeltaPhiMETN()
//it usable directly on the electron list. we have already checked that it is a good electron
double EventCalculator::getDeltaPhiMETN_electron( const unsigned int ielectron, 
						  float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith, bool useArcsin) {
  double deltaT = getDeltaPhiMETN_electron_deltaT(ielectron, otherpt, othereta, otherid, dataJetRes, keith);

  //calculate deltaPhiMETN
  double dp =  getDeltaPhi(pf_els_phi->at(ielectron), getMETphi());
  double dpN;
  if(useArcsin) {
    if( deltaT/getMET() >= 1.0) dpN = dp / (TMath::Pi()/2.0);
    else dpN = dp / asin(deltaT/getMET());
  }
  else dpN = dp / atan2(deltaT, getMET());

  
  return dpN;
}
//the imuon is a different logic than is used for the jet index in getDeltaPhiMETN()
//it usable directly on the muon list. we have already checked that it is a good muon
double EventCalculator::getDeltaPhiMETN_muon( const unsigned int imuon, 
					      float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith, bool useArcsin) {
  double deltaT = getDeltaPhiMETN_muon_deltaT(imuon, otherpt, othereta, otherid, dataJetRes, keith);

  //calculate deltaPhiMETN
  double dp =  getDeltaPhi(pf_mus_phi->at(imuon), getMETphi());
  double dpN;
  if(useArcsin) {
    if( deltaT/getMET() >= 1.0) dpN = dp / (TMath::Pi()/2.0);
    else dpN = dp / asin(deltaT/getMET());
  }
  else dpN = dp / atan2(deltaT, getMET());

  
  return dpN;
}

double EventCalculator::getDeltaPhiMETN( unsigned int goodJetN, float mainpt, float maineta, bool mainid,
					 float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith, bool useArcsin) {//Ben
  
  //find the goodJetN-th good jet -- this is the jet deltaPhiN will be calculated for
  unsigned int ijet = 999999;
  unsigned int goodJetI=0;
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    if (isGoodJet(i, mainpt, maineta, mainid)) {
      if(goodJetI == goodJetN){
	ijet = i;
	break;
      }
      goodJetI++;
    }
  }
  if(ijet == 999999) return -99;

  double deltaT = getDeltaPhiMETN_deltaT(ijet, otherpt, othereta, otherid, dataJetRes, keith);

  //calculate deltaPhiMETN
  double dp =  getDeltaPhi(jets_AK5PF_phi->at(ijet), getMETphi());
  double dpN;
  if(useArcsin) {
    if( deltaT/getMET() >= 1.0) dpN = dp / (TMath::Pi()/2.0);
    else dpN = dp / asin(deltaT/getMET());
  }
  else dpN = dp / atan2(deltaT, getMET());
  
  return dpN;
}


double EventCalculator::getMinDeltaPhiMETN(unsigned int maxjets, float mainpt, float maineta, bool mainid, 
					   float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith, bool includeLeptons, bool useArcsin) {//Ben
  
  double mdpN=1E12;
  
  for (unsigned int i=0; i<maxjets; i++) {
    if(i>=nGoodJets()) break;
    double dpN =  getDeltaPhiMETN(i, mainpt, maineta, mainid, otherpt, othereta, otherid, dataJetRes, keith, useArcsin);//i is for i'th *good* jet, starting at i=0. returns -99 if bad jet.
    if (dpN>=0 && dpN<mdpN) mdpN=dpN;//checking that dpN>=0 shouldn't be necessary after break statement above, but add it anyway 
  }

  //this option causes leptons to be treated like jets in calculating DeltaPhi(lepton, MET). the resolution term remains jets only
  if (includeLeptons) { 

    double  mdpN_e=2e12;
    double  mdpN_m=2e12;
    //need to get DeltaPhiMETN for any good leptons in the event
    for (unsigned int i=0; i <pf_els_pt->size() ; i++) {
      if (isGoodElectron(i,false,10)) { //pt threshold hard-coded to 10
	double  dPN_e=  getDeltaPhiMETN_electron(i, otherpt, othereta, otherid, dataJetRes, keith, useArcsin);
	if (dPN_e>=0 && dPN_e<mdpN_e) mdpN_e=dPN_e;
      }
    }
    for (unsigned int i=0; i <pf_mus_pt->size() ; i++) {
      if (isCleanMuon(i,10)) { //pt threshold hard-coded to 10
	double  dPN_m=  getDeltaPhiMETN_muon(i, otherpt, othereta, otherid, dataJetRes, keith, useArcsin);
	if (dPN_m>=0 && dPN_m<mdpN_m) mdpN_m=dPN_m;
      }
    }
    
    if (mdpN_e < mdpN) mdpN = mdpN_e;
    if (mdpN_m < mdpN) mdpN = mdpN_m;
  }

  return mdpN;
}

double EventCalculator::getMinDeltaPhiNMET(unsigned int maxjets) {//Luke
  
  double mindpn=1E12;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    
    if (isGoodJet(i)) {
      ++ngood;
      double dpn =  getDeltaPhiNMET( i );
      if (dpn<mindpn) mindpn=dpn;
      if (ngood >= maxjets) break;
    }
  }
  return mindpn;
}



double EventCalculator::getMaxDeltaPhiNMET(unsigned int maxjets){
  
  double maxdpn=-99.;

  unsigned int ngood=0;
  //get the maximum angle between the first n jets and MET
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    
    if (isGoodJet(i)) {
      ++ngood;
      double dpn =  getDeltaPhiNMET( i );
      if (dpn>maxdpn) maxdpn=dpn;
      if (ngood >= maxjets) break;
    }
  }
  return maxdpn;
}



double EventCalculator::getTransverseMETSignificance(unsigned int thisJet){

  if(!(thisJet<jets_AK5PF_pt->size())) return -99;

  double deltaT=getTransverseMETError(thisJet);
  double dp = getDeltaPhi( jets_AK5PF_phi->at(thisJet) , getMETphi() );
  return getMET()*sin(dp)/deltaT;
}



double EventCalculator::getMaxTransverseMETSignificance(unsigned int maxjets) {

  double maxTransverseMETSignificance=-100;

  unsigned int ngood=0;
  //get the maximum Transverse MET Significance
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    
    if (isGoodJet(i)) {
      ++ngood;
      double transverseMETSignificance = getTransverseMETSignificance(i);
      if (transverseMETSignificance>maxTransverseMETSignificance) maxTransverseMETSignificance=transverseMETSignificance;
      if (ngood >= maxjets) break;
    }
  }
  
  return maxTransverseMETSignificance;
}

double EventCalculator::getMinTransverseMETSignificance(unsigned int maxjets) {

  double minTransverseMETSignificance=1E12;

  unsigned int ngood=0;
  //get the minimum Transverse MET Significance
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    
    if (isGoodJet(i)) {
      ++ngood;
      double transverseMETSignificance = getTransverseMETSignificance(i);
      if (transverseMETSignificance<minTransverseMETSignificance) minTransverseMETSignificance=transverseMETSignificance;
      if (ngood >= maxjets) break;
    }
  }
  
  return minTransverseMETSignificance;
}

double EventCalculator::getDeltaPhiMET(unsigned int n, float ptThreshold, bool bjetsonly){
  
  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
    
    if (isGoodJet(i,ptThreshold)) {
      if(bjetsonly && passBTagger(i)==0) continue;
      ngood++;
      if (ngood==n) return  getDeltaPhi( jets_AK5PF_phi->at(i) , getMETphi());
    }
  }
  return 0;
}

double EventCalculator::getMinDeltaPhiMET30(unsigned int maxjets, bool bjetsonly) {

  double mindp=99;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    
    if (isGoodJet30(i)) {
      if(bjetsonly && passBTagger(i)==0) continue;
      ++ngood;
      double dp =  getDeltaPhi( jets_AK5PF_phi->at(i) , getMETphi());
      if (dp<mindp) mindp=dp;
      if (ngood >= maxjets) break;
    }
  }
  
  return mindp;
}

double EventCalculator::getMinDeltaPhiMET30_eta5(unsigned int maxjets) {

  double mindp=99;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    
    if (getJetPt(i) > 30 && fabs(jets_AK5PF_eta->at(i)) < 5 && jetPassLooseID(i)) {
      ++ngood;
      double dp =  getDeltaPhi( jets_AK5PF_phi->at(i) , getMETphi());
      if (dp<mindp) mindp=dp;
      if (ngood >= maxjets) break;
    }
  }
  
  return mindp;
}

double EventCalculator::getMinDeltaPhiMET30_eta5_noId(unsigned int maxjets) {

  double mindp=99;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    
    if (getJetPt(i) > 30 && fabs(jets_AK5PF_eta->at(i)) < 5) {
      ++ngood;
      double dp =  getDeltaPhi( jets_AK5PF_phi->at(i) , getMETphi());
      if (dp<mindp) mindp=dp;
      if (ngood >= maxjets) break;
    }
  }
  
  return mindp;
}

double EventCalculator::getMaxDeltaPhiMET(unsigned int maxjets) {

  double maxdp=0;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    
    if (isGoodJet(i)) {
      ++ngood;
      double dp =  getDeltaPhi( jets_AK5PF_phi->at(i) , getMETphi());
      if (dp>maxdp) maxdp=dp;
      if (ngood >= maxjets) break;
    }
  }
  
  return maxdp;
}
double EventCalculator::getMaxDeltaPhiMET30(unsigned int maxjets) {

  double maxdp=0;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    
    if (isGoodJet30(i)) {
      ++ngood;
      double dp =  getDeltaPhi( jets_AK5PF_phi->at(i) , getMETphi());
      if (dp>maxdp) maxdp=dp;
      if (ngood >= maxjets) break;
    }
  }
  
  return maxdp;
}

double EventCalculator::getMaxDeltaPhiMET30_eta5(unsigned int maxjets) {

  double maxdp=0;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    
    if (getJetPt(i) > 30 && fabs(jets_AK5PF_eta->at(i)) < 5 && jetPassLooseID(i)) {
      ++ngood;
      double dp =  getDeltaPhi( jets_AK5PF_phi->at(i) , getMETphi());
      if (dp>maxdp) maxdp=dp;
      if (ngood >= maxjets) break;
    }
  }
  
  return maxdp;
}
double EventCalculator::getMaxDeltaPhiMET30_eta5_noId(unsigned int maxjets) {

  double maxdp=0;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    
    if (getJetPt(i) > 30 && fabs(jets_AK5PF_eta->at(i)) < 5) {
      ++ngood;
      double dp =  getDeltaPhi( jets_AK5PF_phi->at(i) , getMETphi());
      if (dp>maxdp) maxdp=dp;
      if (ngood >= maxjets) break;
    }
  }
  
  return maxdp;
}


double EventCalculator::getMinDeltaPhiMETTaus() {

  double mindp=99;

/* //FIXME CFA
  for (unsigned int i=0; i< taus_pt->size(); i++) {
    if (taus_pt->at(i) > 30) {
      double dp =  getDeltaPhi( taus_phi->at(i) , getMETphi());
      if (dp<mindp) mindp=dp;
    }
  }
*/
  return mindp;
}

double EventCalculator::getMinDeltaPhiMETMuons(unsigned int maxmuons) {

  double mindp=99;

  unsigned int ngood=0;
  //get the minimum angle between the first n muons and MET
  for (unsigned int i=0; i< mus_pt->size(); i++) {
    
    if (isGoodRecoMuon(i,true,10)) {//disable iso
      ++ngood;
      double dp =  getDeltaPhi( mus_phi->at(i) , getMETphi());
      if (dp<mindp) mindp=dp;
      if (ngood >= maxmuons) break;
    }
  }
  
  return mindp;
}


void EventCalculator::getCorrectedMET(float& correctedMET, float& correctedMETPhi) {

  correctedMET=pfTypeImets_et->at(0);
  correctedMETPhi = pfTypeImets_phi->at(0);

  if (theJESType_==kJES0 && theJERType_==kJER0) return;

  double METx = correctedMET * cos(correctedMETPhi);
  double METy = correctedMET * sin(correctedMETPhi);

  for(unsigned int thisJet = 0; thisJet < jets_AK5PF_pt->size(); thisJet++)   {
    if (isCleanJet(thisJet) ){//this only has an effect for recopfjets
      if ( jets_AK5PF_pt->at(thisJet) >10) {
	METx += jets_AK5PF_pt->at(thisJet) * cos(jets_AK5PF_phi->at(thisJet));
	METx -= getJetPt(thisJet,false) * cos(jets_AK5PF_phi->at(thisJet));
	METy += jets_AK5PF_pt->at(thisJet) * sin(jets_AK5PF_phi->at(thisJet));
	METy -= getJetPt(thisJet,false) * sin(jets_AK5PF_phi->at(thisJet));
      }
    }
  }
  correctedMET = sqrt(METx*METx + METy*METy);
  correctedMETPhi = atan2(METy,METx);
}


void EventCalculator::getUncorrectedMET(float& uncorrectedMET, float& uncorrectedMETPhi) {

  uncorrectedMET=pfmets_et->at(0);
  uncorrectedMETPhi = pfmets_phi->at(0);

  if (theJESType_==kJES0 && theJERType_==kJER0) return;

  double METx = uncorrectedMET * cos(uncorrectedMETPhi);
  double METy = uncorrectedMET * sin(uncorrectedMETPhi);

  //  cout<<endl<< " == computing MET "<<endl; //JMT DEBUG
  for(unsigned int thisJet = 0; thisJet < jets_AK5PF_pt->size(); thisJet++)
    {
      if (isCleanJet(thisJet) ){//this only has an effect for recopfjets
	if ( jets_AK5PF_pt->at(thisJet)>10) {
	  const float uncorrectedpt = jets_AK5PF_pt->at(thisJet) * jets_AK5PF_corrFactorRaw->at(thisJet);
	  METx += uncorrectedpt * cos(jets_AK5PF_phi->at(thisJet)); //used to use uncorrected phi. should be no change!
	  METx -= getUncorrectedJetPt(thisJet,true) * cos(jets_AK5PF_phi->at(thisJet));
	  METy += uncorrectedpt * sin(jets_AK5PF_phi->at(thisJet));
	  METy -= getUncorrectedJetPt(thisJet,true) * sin(jets_AK5PF_phi->at(thisJet));
	}
      }
    }
  uncorrectedMET = sqrt(METx*METx + METy*METy);
  uncorrectedMETPhi = atan2(METy,METx);
}

int EventCalculator::doPBNR() {
  //RA2 particle-based noise rejection

  bool nhBad=false;
  bool phBad=false;
  for (unsigned int it = 0; it<jets_AK5PF_pt->size(); it++) {
    //cfA version from Keith
    double NHF = jets_AK5PF_neutralHadE->at(it)/(jets_AK5PF_energy->at(it)*jets_AK5PF_corrFactorRaw->at(it));
    double PEF = jets_AK5PF_photonEnergy->at(it)/(jets_AK5PF_energy->at(it)*jets_AK5PF_corrFactorRaw->at(it));
    if (NHF > 0.9)  nhBad = true;
    if (PEF > 0.95) phBad = true;
  }
  
  if (nhBad && phBad) return -3;
  else if (phBad) return -2;
  else if (nhBad) return -1;
  return 1;
}


bool EventCalculator::isCleanJet(const unsigned int ijet){
    
  if(theJetType_ == kPF2PAT) return true; //pf2pat jets are clean by construction
  assert(0);//FIXME CFA
  /* 
  //if it's near a good muon, it's not clean
  bool isNearMuon = false;
  for ( unsigned int j = 0; j< pf_mus_pt->size(); j++) {
    if( isGoodMuon(j) && jmt::deltaR( myJetsPF->at(ijet).eta, myJetsPF->at(ijet).phi,myMuonsPF->at(j).eta, myMuonsPF->at(j).phi)<0.1 ){
      isNearMuon = true; break;
    }
  }
  
  if(isNearMuon) return false;
  
  //if it's near a good electron, it's not clean
  bool isNearElectron = false;
  for ( unsigned int j = 0; j< myElectronsPF->size(); j++) {
    if( isGoodElectron(j) && jmt::deltaR( myJetsPF->at(ijet).eta, myJetsPF->at(ijet).phi,myElectronsPF->at(j).eta, myElectronsPF->at(j).phi)<0.3 ){
      isNearElectron = true; break;
    }
  }
    
  if(isNearElectron) return false;
  */
  return true;
}


//i'm tempted to muck with the jet pT threhsold....
float EventCalculator::deltaRBestTwoBjets(int & jetindex1, int & jetindex2 ) {

  jetindex1 = -1;
  jetindex2=-2;

  float bestcsv1=-98;
  float bestcsv2=-99;

  for (unsigned int it = 0; it<jets_AK5PF_pt->size(); it++) {
    if (isGoodJet(it) ) {
      float csv = jets_AK5PF_btag_secVertexCombined->at(it) ;

      if (csv> bestcsv1) {
	bestcsv2 = bestcsv1;
	bestcsv1 = csv;
	jetindex2=jetindex1;
	jetindex1=it;
      }
      else if (csv>bestcsv2) {
	bestcsv2 = csv;
	jetindex2=it;
      }
    }
  }

  float dr=-1;
  if (jetindex1>=0 && jetindex2>=0) dr= jmt::deltaR( jets_AK5PF_eta->at(jetindex1), jets_AK5PF_phi->at(jetindex1),jets_AK5PF_eta->at(jetindex2), jets_AK5PF_phi->at(jetindex2));

  return dr;
}

float EventCalculator::deltaRClosestTwoBjets(int & jetindex1, int & jetindex2 ) {

  float mindr=1e9;
  jetindex1 = -1;
  jetindex2=-1;

  for (unsigned int j1 = 0; j1<jets_AK5PF_pt->size(); j1++) {
    if (isGoodJet(j1) && passBTagger(j1)) {

      //avoid duplicates by starting loop at the jet after j1
      for (unsigned int j2 = j1+1; j2<jets_AK5PF_pt->size(); j2++) {
	
	if (isGoodJet(j2) && passBTagger(j2)) {
	  float   dr= jmt::deltaR( jets_AK5PF_eta->at(j1), jets_AK5PF_phi->at(j1),jets_AK5PF_eta->at(j2), jets_AK5PF_phi->at(j2));
	  if (dr<mindr) {
	    mindr=dr;
	    jetindex1 = j1;
	    jetindex2 = j2;
	  }
	}
      }
    }
  }

  return mindr;
}


bool EventCalculator::isGoodJet(const unsigned int ijet, const float pTthreshold, const float etaMax, const bool jetid) {

  if ( getJetPt(ijet) <pTthreshold) return false;
  if ( fabs(jets_AK5PF_eta->at(ijet)) > etaMax) return false;
  if ( jetid && !jetPassLooseID(ijet) ) return false;

  //do manual cleaning for reco pfjets
  if(theJetType_ == kRECOPF){
    if ( !isCleanJet(ijet) ) return false;
  }

  return true;
}

bool EventCalculator::passBadJetFilter() {

  //look for any jet with pT>30 and failing jet id

  for (unsigned int i=0; i < jets_AK5PF_pt->size(); ++i) {
    //bool argument turns jetID off. 99 is the eta max
    if ( isGoodJet(i,30,99,false) && !isGoodJet(i,30,99,true) ) return false;
  }

  return true;
}

bool EventCalculator::isGoodJetMHT(unsigned int ijet) {

  if ( getJetPt(ijet) <30) return false;
  if ( fabs(jets_AK5PF_eta->at(ijet)) > 5) return false;
  //no jet id for MHT

  return true;
}


bool EventCalculator::passBTagger(int ijet, BTaggerType btagger ) {

  //if we are passed the dummy value, then set to the default type
  if (btagger==Nbtaggers) btagger = theBTaggerType_;

  if (btagger==kSSVM){
    return ( jets_AK5PF_btag_secVertexHighEff->at(ijet) >= 1.74 );
  }
  else if (btagger==kTCHET){
    return ( jets_AK5PF_btag_TC_highEff->at(ijet) >= 10.2);
  }
  else if (btagger==kSSVHPT) return jets_AK5PF_btag_secVertexHighPur->at(ijet) >= 2;
  else if (btagger==kTCHPT ) return jets_AK5PF_btag_TC_highPur->at(ijet) >= 3.41;
  else if (btagger==kTCHPM ) return jets_AK5PF_btag_TC_highPur->at(ijet) >= 1.93;
  else if (btagger==kCSVM  ) return jets_AK5PF_btag_secVertexCombined->at(ijet) >=0.679;
  else if (btagger==kCSVL  ) return jets_AK5PF_btag_secVertexCombined->at(ijet) >=0.244;
  else{
    cout << "Invalid b tagger!" << endl;
    assert(0);
    return false;
  }
}

float EventCalculator::getJetCSV(unsigned int ijet){
  return jets_AK5PF_btag_secVertexCombined->at(ijet);
}

unsigned int EventCalculator::nGoodJets() {
  
  unsigned int njets=0;
  for (unsigned int i=0; i < jets_AK5PF_pt->size(); ++i) {
    if (isGoodJet(i) )   njets++;
  }
  return njets;
}

//this was for testing only...to fill these special histograms
/*
unsigned int EventCalculator::nGoodJets(TH2D* count,TH2D* unc,TH2D* l2l3) {
  
  unsigned int njets=0;
  for (unsigned int i=0; i < myJetsPF->size(); ++i) {
    if (isGoodJet(i) )   {
      njets++;
      
      count->Fill(myJetsPF->at(i).eta, myJetsPF->at(i).pt);
      unc->Fill(myJetsPF->at(i).eta, myJetsPF->at(i).pt,100*0.5*( myJetsPFhelper->at(i).jetUncPlus +myJetsPFhelper->at(i).jetUncMinus ));

      JetCorrector_->setJetEta( myJetsPF->at(i).eta);
      JetCorrector_->setJetPt( myJetsPF->at(i).uncor_pt );
      JetCorrector_->setJetA( myJetsPF->at(i).jetArea );
      JetCorrector_->setRho( sdouble_kt6pfjets_rho_value ); 
      std::vector<float> factors = JetCorrector_->getSubCorrections();
      float L2L3Residualonly = factors[3]/factors[2];
      
      l2l3->Fill(myJetsPF->at(i).eta, myJetsPF->at(i).pt, 100*(L2L3Residualonly-1));
    }
  }
  return njets;
}
*/

unsigned int EventCalculator::nGoodJets30() {
  
  unsigned int njets=0;
  for (unsigned int i=0; i < jets_AK5PF_pt->size(); ++i) {
    if (isGoodJet30(i) )   njets++;
  }
  return njets;
}

unsigned int EventCalculator::nGoodBJets( float ptthreshold,BTaggerType btagger) {
  unsigned int nb=0;
  for (unsigned int i = 0; i < jets_AK5PF_pt->size(); ++i) {
    if (isGoodJet(i,ptthreshold) ) {
      if ( passBTagger(i,btagger) ) nb++;
    }
  }
  return nb;
}

unsigned int EventCalculator::nTrueBJets() {
  unsigned int nb=0;

  if (jets_AK5PF_partonFlavour==0) return nb; //not actually in the tree? FIXME CFA

  for (unsigned int i = 0; i < jets_AK5PF_pt->size(); ++i) {
    if (isGoodJet(i,30) ) {
      if ( abs(jets_AK5PF_partonFlavour->at(i))==5 ) nb++;
    }
  }
  return nb;
}

void EventCalculator::getSmearedUnclusteredMET(float & myMET, float & myMETphi) {
  assert(theJetType_ == kPF2PAT); //enforce that jets and leptons are already disambiguated in the ntuple

  //input values myMET and myMETphi should correspond to the 
  //type of MET being used (corrected or uncorrected)

  float myMETx = myMET * cos(myMETphi);
  float myMETy = myMET * sin(myMETphi);

  //first remove jets from MET
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
    //    if (isCleanJet(i) ){//this only has an effect for recopfjets    
    if ( jets_AK5PF_pt->at(i) >10 ) {
      if (theMETType_ == kPFMET) {    //remove the uncorrected jet pT from the raw MET
	myMETx += getUncorrectedJetPt(i) * cos(jets_AK5PF_phi->at(i));
	myMETy += getUncorrectedJetPt(i) * sin(jets_AK5PF_phi->at(i));
      }
      else if (theMETType_ == kPFMETTypeI) { //remove the corrected jet pT from the corrected MET
	myMETx += getJetPt(i) * cos(jets_AK5PF_phi->at(i));
	myMETy += getJetPt(i) * sin(jets_AK5PF_phi->at(i));
      }
      else assert(0);
    }
  }
  //then muons
  for ( unsigned int i = 0; i<pf_mus_pt->size() ; i++) {
    //    if (isCleanMuon(i)) {
      myMETx += pf_mus_pt->at(i) * cos(pf_mus_phi->at(i));
      myMETy += pf_mus_pt->at(i) * sin(pf_mus_phi->at(i));
      //    }
  }
  //electrons
  for ( unsigned int i = 0; i<pf_els_pt->size() ; i++) {
    //    if (isGoodElectron(i)) {
      myMETx += pf_els_pt->at(i) * cos(pf_els_phi->at(i));
      myMETy += pf_els_pt->at(i) * sin(pf_els_phi->at(i));
      //    }
  }


//   if (sqrt(myMETx*myMETx + myMETy*myMETy) >100) {
//     cout<<myMET<<" "<<sqrt(myMETx*myMETx + myMETy*myMETy)<<"\t"<<atan2(myMETy,myMETx)<<endl; dumpEvent();
//   }

  //now we've got the unclustered component only in myMETxy
  float factor=1;
  if (theMETuncType_ ==kMETuncDown) factor=0.9;
  else if (theMETuncType_ ==kMETuncUp) factor=1.1;
  else {assert(0);}
  myMETx *= factor;
  myMETy *= factor;

  //now repeat everything but change all += to -=

  //jets
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
    if ( jets_AK5PF_pt->at(i) >10 ) {
      if (theMETType_ == kPFMET) {    //add the uncorrected jet pT to the raw MET
	myMETx -= getUncorrectedJetPt(i) * cos(jets_AK5PF_phi->at(i));
	myMETy -= getUncorrectedJetPt(i) * sin(jets_AK5PF_phi->at(i));
      }
      else if (theMETType_ == kPFMETTypeI) { //add the corrected jet pT to the corrected MET
	myMETx -= getJetPt(i) * cos(jets_AK5PF_phi->at(i));
	myMETy -= getJetPt(i) * sin(jets_AK5PF_phi->at(i));
      }
      else assert(0);
    }
  }
  //muons
  for ( unsigned int i = 0; i< pf_mus_pt->size() ; i++) {
    //    if (isCleanMuon(i)) {
      myMETx -= pf_mus_pt->at(i) * cos(pf_mus_phi->at(i));
      myMETy -= pf_mus_pt->at(i) * sin(pf_mus_phi->at(i));
      //    }
  }
  //electrons
  for ( unsigned int i = 0; i< pf_els_pt->size() ; i++) {
    //    if (isGoodElectron(i)) {
      myMETx -= pf_els_pt->at(i) * cos(pf_els_phi->at(i));
      myMETy -= pf_els_pt->at(i) * sin(pf_els_phi->at(i));
      //    }
  }

  //  cout<<sqrt(myMETx*myMETx + myMETy*myMETy)<<endl;
  myMET = sqrt(myMETx*myMETx + myMETy*myMETy);
  myMETphi = atan2(myMETy,myMETx);
}

float EventCalculator::getJERbiasFactor(unsigned int ijet) {
  float abseta = fabs(jets_AK5PF_eta->at(ijet));

  if (theJERType_ == kJERup || theJERType_ ==kJERbias ||theJERType_ ==kJERdown) {
    //old numbers were from JME-10-014-pas
    //now replaced by 2011 numbers from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    float m = 0;
    if (theJERType_ == kJERup) m=1;
    else if (theJERType_==kJERdown) m=-1;

    if (abseta < 0.5)                     return 0.052   + m*0.063;
    else if (abseta <1.1 && abseta>=0.5)  return 0.057   + m*0.057;
    else if (abseta <1.7 && abseta>=1.1)  return 0.096   + m*0.065;
    else if (abseta <2.3 && abseta>=1.7)  return 0.134   + m*0.094;
    else if (abseta<=5 && abseta>= 2.3)   return 0.288   + m*0.200;
    else return 0; //eta>5 not given
  }
  else if (theJERType_ == kJERra2) {
    //from the RA2 AN 2011/247
    float m = 0;

    if (abseta < 0.5) {
      return 0.052 + m*0.065;
    }
    else if (abseta <1.1 && abseta>=0.5) {
      return 0.057 + m*0.06;
    }
    else if (abseta <1.7 && abseta>=1.1) {
      return  0.096 + m*0.07;
    }
    else if (abseta <2.3 && abseta>=1.7) {
      return  0.134 + m*0.1;
    }
    else if (abseta<=5 && abseta>= 2.3) {
      return 0.288 + m*0.222;
    }
    else { //eta>5 not given
      return 0;
    }
  }

  assert(0);
  return 0;
}


float EventCalculator::getJESUncertainty( unsigned int ijet, bool addL2L3toJES ) {
  assert ( theJESType_ != kJES0 );

  float uncertainty = -1000;


/* for timing tests only
  if ( watch_==0 ) watch_ = new TStopwatch(); //ctor starts the timer
  else watch_->Continue(); //kFALSE == don't reset the timer
*/

  bool uncarg=false;

  if ( theJESType_ == kJESup ) {
    //    uncertainty = myJetsPFhelper->at(ijet).jetUncPlus; //cfa ntuples don't appear to have these stored
    uncarg=true;
  }
  else if (theJESType_ == kJESdown) {
    //    uncertainty = myJetsPFhelper->at(ijet).jetUncMinus;  //cfa ntuples don't appear to have these stored
    uncarg=false;
  }
  
  //now get the uncertainty (from the files)
  if ( fabs(jets_AK5PF_eta->at(ijet)) <5 ) {
    jecUnc_->setJetEta( jets_AK5PF_eta->at(ijet));
    jecUnc_->setJetPt( jets_AK5PF_pt->at(ijet)); // here you must use the CORRECTED jet pt
    uncertainty = jecUnc_->getUncertainty(uncarg);
    
    //      cout<<"[EventCalculator::getJESUncertainty] "<<myJetsPFhelper->at(ijet).jetUncPlus<<"\t"<<uncertainty<<endl;
    //      if ( fabs(  myJetsPFhelper->at(ijet).jetUncPlus - uncertainty) >0.01) cout<<"[EventCalculator::getJESUncertainty] something Bad!"<<endl;
    
  }

  //  watch_->Stop(); //for timing tests only

  if (fabs(uncertainty) > 1) uncertainty=0; //important sanity check

  //cout<<"\t JES unc ("<< jets_AK5PF_pt->at(ijet) <<") = "<<uncertainty<<endl;//JMT DEBUG

  if (addL2L3toJES ) { //FIXME tbc that this is really the right recipe. THis is what was used for 2011
    JetCorrector_->setJetEta( jets_AK5PF_eta->at(ijet) );
    JetCorrector_->setJetPt( jets_AK5PF_rawPt->at(ijet) );
    JetCorrector_->setJetA( jets_AK5PF_area->at(ijet) );
    JetCorrector_->setRho( rho_kt6PFJetsForIsolation2012 ); //FIXME CFA this may be right but i'm not sure
    std::vector<float> factors = JetCorrector_->getSubCorrections();
    float L2L3Residual = factors[3]/factors[2] - 1;

    uncertainty = sqrt(uncertainty*uncertainty + pow(L2L3Residual,2) );

    //    cout<<" ++ "<<L2L3Residual<<" = "<<uncertainty<<endl; //JMT DEBUG
  }

  return uncertainty;

}

float EventCalculator::getJetPt( unsigned int ijet, bool addL2L3toJES ) {

  //i am going to write this to allow simultaneous use of JER and JES
  //(hence use if repeated 'if' instead of 'else if')

  float pt = jets_AK5PF_pt->at(ijet);
  //if ( theJESType_ == kJES0 && theJERType_ == kJER0) return pt;

    /* leave this code intact but we don't really want to use this right now
  if (theJESType_ == kJESFLY) {
    assert(0);
    JetCorrector_->setJetEta( myJetsPF->at(ijet).eta);
    JetCorrector_->setJetPt( myJetsPF->at(ijet).uncor_pt );
    JetCorrector_->setJetA( myJetsPF->at(ijet).jetArea );
    JetCorrector_->setRho( sdouble_kt6pfjets_rho_value ); 
    
    std::vector<float> factors = JetCorrector_->getSubCorrections();
    //[0] = L1Fast, [1] = [0]*L2Residual, [2] = [1]*L3Absolute, [3] = [2]*L2L3Residual

    //this will return the L3Absolute corrected pt (i.e. the exact same pt as is normally used for MC)
    //pt = myJetsPF->at(ijet).uncor_pt * factors[2];

    //this will return the pt with "ANTI" L2L3Residual corrections!!
    float L2L3Residualonly = factors[3]/factors[2];
    pt = myJetsPF->at(ijet).uncor_pt * factors[2]/L2L3Residualonly;
    
    //for study of L2L3Res
	    //if (pt>30)  cout<<ijet<<" pT L2L3Residual Unc+ Unc- : "<<pt<<"\t"<<1-L2L3Residualonly<<" "<<myJetsPFhelper->at(ijet).jetUncPlus<<" "<<myJetsPFhelper->at(ijet).jetUncMinus <<" "<< sqrt( pow(1-L2L3Residualonly,2) + pow(myJetsPFhelper->at(ijet).jetUncPlus,2) ) <<endl;

  }
end of kJESFLY block  */
  
 //first JER
  if ( theJERType_ != kJER0 ) {
    float genpt = jets_AK5PF_gen_pt->at(ijet);
    if (genpt > 15) {
      float factor = getJERbiasFactor(ijet);
      float deltapt = (pt - genpt) * factor;
      float frac = (pt+deltapt)/pt;
      float ptscale = frac>0 ? frac : 0;
      pt *= ptscale;
    }
    
  }
  
  //  cout<<endl<< " == computing jet pT "<<endl; //JMT DEBUG
  //then JES
  if ( theJESType_ == kJESup ) {
    //sanity check on unc is now moved inside of getJESUncertainty
    const    float unc = getJESUncertainty(ijet,addL2L3toJES);
    pt *= (1+ unc);
  }
  else if (theJESType_ == kJESdown) {
    const    float unc = getJESUncertainty(ijet,addL2L3toJES);
    pt *= (1- unc);
  }
  
  return pt;
  
}


//code here largely copied and pasted from getJetPt(). could be done better
float EventCalculator::getUncorrectedJetPt( unsigned int ijet, bool addL2L3toJES ) {

  //i am going to write this to allow simultaneous use of JER and JES
  //(hence use if repeated 'if' instead of 'else if')
  //I hope this is sensible!

  //this should be uncorrected pt (cf Paul Geffert)
  float pt = jets_AK5PF_pt->at(ijet) * jets_AK5PF_corrFactorRaw->at(ijet);

  //first JER
  if ( theJERType_ != kJER0 ) {
    float genpt = jets_AK5PF_gen_pt->at(ijet);
    if (genpt > 15) {
      float recopt = jets_AK5PF_pt->at(ijet); //use corrected jet pt here
      float factor = getJERbiasFactor(ijet);
      float deltapt = (recopt - genpt) * factor;
      float frac = (recopt+deltapt)/recopt;
      float ptscale = frac>0 ? frac : 0;
      pt *= ptscale;
    }
  }
  //then JES
  if ( theJESType_ == kJESup ) {
    //sanity check on unc is now moved inside of getJESUncertainty.
    //it will return 0 in case of failed sanity check
    const    float unc = getJESUncertainty(ijet,addL2L3toJES);
    pt *= (1+unc);
  }
  else if (theJESType_ == kJESdown) {
    const    float unc = getJESUncertainty(ijet,addL2L3toJES);
    pt *= (1- unc);
  }

  return pt;

}


float EventCalculator::getJetPx( unsigned int ijet ) {
  return getJetPt(ijet) * cos(jets_AK5PF_phi->at(ijet));
}

float EventCalculator::getJetPy( unsigned int ijet ) {
  return getJetPt(ijet) * sin(jets_AK5PF_phi->at(ijet));
}

float EventCalculator::getJetPz( unsigned int ijet ) {
  return getJetPt(ijet) * sinh(jets_AK5PF_eta->at(ijet));
}

float EventCalculator::getJetEnergy( unsigned int ijet ) {
  //need to move to using jecFactor 
  return jets_AK5PF_energy->at(ijet);
}

bool EventCalculator::jetPassLooseID( unsigned int ijet ) {

  //want the uncorrected energy
  const float jetenergy = jets_AK5PF_energy->at(ijet) * jets_AK5PF_corrFactorRaw->at(ijet);
  const int numConst = jets_AK5PF_mu_Mult->at(ijet)+jets_AK5PF_neutral_Mult->at(ijet)+jets_AK5PF_chg_Mult->at(ijet); //stealing from Keith

  if (jets_AK5PF_neutralHadE->at(ijet) /jetenergy < 0.99
     && jets_AK5PF_neutralEmE->at(ijet) / jetenergy < 0.99
     && numConst > 1
     && ( fabs(jets_AK5PF_eta->at(ijet))>=2.4 
	  || (fabs(jets_AK5PF_eta->at(ijet))<2.4 && jets_AK5PF_chgHadE->at(ijet)/jetenergy>0))
     && ( fabs(jets_AK5PF_eta->at(ijet))>=2.4 
	  || (fabs(jets_AK5PF_eta->at(ijet))<2.4 && jets_AK5PF_chgEmE->at(ijet)/jetenergy<0.99))
     && ( fabs(jets_AK5PF_eta->at(ijet))>=2.4 
	  || (fabs(jets_AK5PF_eta->at(ijet))<2.4 && jets_AK5PF_chg_Mult->at(ijet)>0))
      ) {
    return true;
  }

  return false;
}

float EventCalculator::getRelIsoForIsolationStudyEle() {

  //loop over electrons
  //consider any electron that passes all cuts except reliso
  //return the most isolated reliso of that group

  float bestRelIso=1e9;

  for (unsigned int iele=0; iele < pf_els_pt->size(); ++iele) {
    if ( isGoodElectron(iele,true)) { //true means disable RelIso cut
      float reliso = (pf_els_chargedHadronIso->at(iele)
		      + pf_els_photonIso->at(iele) 
		      + pf_els_neutralHadronIso->at(iele))/pf_els_pt->at(iele);
      if (reliso < bestRelIso) bestRelIso = reliso;
    }
  }

  //if there are no good electrons we'll get 1e9 as the return value
  return bestRelIso;
}

float EventCalculator::getRelIsoForIsolationStudyMuon() {

  //loop over muon
  //consider any muon that passes all cuts except reliso
  //return the most isolated reliso of that group

  float bestRelIso=1e9;
/* FIXME CFA
  for (unsigned int imu=0; imu < mus_pt->size(); ++imu) {
    if ( isGoodRecoMuon(imu,true,10)) { //true means disable RelIso cut
      float reliso = (myMuonsRECO->at(imu).trackIso
		      + myMuonsRECO->at(imu).ecalIso 
		      + myMuonsRECO->at(imu).hcalIso)/myMuonsRECO->at(imu).pt;
      
      if (reliso < bestRelIso) bestRelIso = reliso;
    }
  }
*/
  //if there are no good muons we'll get 1e9 as the return value
  return bestRelIso;
}

float EventCalculator::elePtOfN(unsigned int n, const float ptthreshold) {

  unsigned int ngood=0;
  for (unsigned int i=0; i < pf_els_pt->size(); i++) {
    if(isGoodElectron(i,false,ptthreshold)){
      ngood++;
      if (ngood==n) return pf_els_et->at(i);
    }
  }
  return 0;
}

float EventCalculator::eleEtaOfN(unsigned int n, const float ptthreshold) {

  unsigned int ngood=0;
  for (unsigned int i=0; i < pf_els_pt->size(); i++) {
    if(isGoodElectron(i,false,ptthreshold)){
      ngood++;
      if (ngood==n) return pf_els_eta->at(i);
    }
  }
  return 0;
}

float EventCalculator::elePhiOfN(unsigned int n, const float ptthreshold) {

  unsigned int ngood=0;
  for (unsigned int i=0; i < pf_els_pt->size(); i++) {
    if(isGoodElectron(i,false,ptthreshold)){
      ngood++;
      if (ngood==n) return pf_els_phi->at(i);
    }
  }
  return 0;
}

int EventCalculator::eleChargeOfN(unsigned int n, const float ptthreshold) {

  unsigned int ngood=0;
  for (unsigned int i=0; i < pf_els_pt->size(); i++) {
    if(isGoodElectron(i,false,ptthreshold)){
      ngood++;
      if (ngood==n) return TMath::Nint(pf_els_charge->at(i));
    }
  }
  return 0;
}

int EventCalculator::muonChargeOfN(unsigned int n, const float ptthreshold) {

  unsigned int ngood=0;
  for (unsigned int i=0; i < pf_mus_pt->size(); i++) {
    if (isCleanMuon(i,ptthreshold)) {
      ngood++;
      if (ngood==n) return TMath::Nint(pf_mus_charge->at(i));
    }
  }
  return 0;
}

float EventCalculator::muonPtOfN(unsigned int n, const float ptthreshold) {

  unsigned int ngood=0;
  for (unsigned int i=0; i < pf_mus_pt->size(); i++) {
    if (isCleanMuon(i,ptthreshold)) {
      ngood++;
      if (ngood==n) return pf_mus_pt->at(i);
    }
  }
  return 0;
}

float EventCalculator::muonEtaOfN(unsigned int n, const float ptthreshold) {

  unsigned int ngood=0;
  for (unsigned int i=0; i < pf_mus_pt->size(); i++) {
    if (isCleanMuon(i,ptthreshold)) {
      ngood++;
      if (ngood==n) return pf_mus_eta->at(i);
    }
  }
  return 0;
}


float EventCalculator::muonPhiOfN(unsigned int n, const float ptthreshold) {

  unsigned int ngood=0;
  for (unsigned int i=0; i < pf_mus_phi->size(); i++) {
    if (isCleanMuon(i,ptthreshold)) {
      ngood++;
      if (ngood==n) return pf_mus_phi->at(i);
    }
  }
  return 0;
}

float EventCalculator::muonIsoOfN(unsigned int n, const float ptthreshold) {

  unsigned int ngood=0;
  for (unsigned int i=0; i < pf_mus_phi->size(); i++) {
    if (isCleanMuon(i,ptthreshold)) {
      ngood++;
      if (ngood==n) {
	float reliso = ((pf_mus_chargedHadronIso->at(i) 
			 + pf_mus_photonIso->at(i) 
			 + pf_mus_neutralHadronIso->at(i))/pf_mus_pt->at(i));
	return reliso;
      }
    }
  }
  return 0;
}


float EventCalculator::muonChHadIsoOfN(unsigned int n, const float ptthreshold) {

  unsigned int ngood=0;
  for (unsigned int i=0; i < pf_mus_pt->size(); i++) {
    if (isCleanMuon(i,ptthreshold)) {
      ngood++;
      if (ngood==n) {
	return pf_mus_chargedHadronIso->at(i) /pf_mus_pt->at(i);
      }
    }
  }
  return 0;
}

float EventCalculator::muonPhotonIsoOfN(unsigned int n, const float ptthreshold) {

  unsigned int ngood=0;
  for (unsigned int i=0; i < pf_mus_pt->size(); i++) {
    if (isCleanMuon(i,ptthreshold)) {
      ngood++;
      if (ngood==n) {
	return pf_mus_photonIso->at(i) /pf_mus_pt->at(i);
      }
    }
  }
  return 0;
}

float EventCalculator::muonNeutralHadIsoOfN(unsigned int n, const float ptthreshold) {

  unsigned int ngood=0;
  for (unsigned int i=0; i < pf_mus_pt->size(); i++) {
    if (isCleanMuon(i,ptthreshold)) {
      ngood++;
      if (ngood==n) {
	return pf_mus_neutralHadronIso->at(i)/pf_mus_pt->at(i);
      }
    }
  }
  return 0;
}

//FIXME CFA
float EventCalculator::recoMuonPtOfN(unsigned int n, const float ptthreshold) {
/*
  unsigned int ngood=0;
  for (unsigned int imu=0; imu < myMuonsRECO->size(); ++imu) {
    if ( isGoodRecoMuon(imu,true,ptthreshold)) { //true means disable RelIso cut
      ngood++;
      if (ngood==n) {
	return myMuonsRECO->at(imu).pt;
      }
    }
  }
*/
  return 0;
}

//FIXME CFA
float EventCalculator::recoMuonEtaOfN(unsigned int n, const float ptthreshold) {
/*
  unsigned int ngood=0;
  for (unsigned int imu=0; imu < myMuonsRECO->size(); ++imu) {
    if ( isGoodRecoMuon(imu,true,ptthreshold)) { //true means disable RelIso cut
      ngood++;
      if (ngood==n) {
	return myMuonsRECO->at(imu).eta;;
      }
    }
  }*/
  return 0;
}


//FIXME CFA
float EventCalculator::recoMuonPhiOfN(unsigned int n, const float ptthreshold) {
/*
  unsigned int ngood=0;
  for (unsigned int imu=0; imu < myMuonsRECO->size(); ++imu) {
    if ( isGoodRecoMuon(imu,true,ptthreshold)) { //true means disable RelIso cut
      ngood++;
      if (ngood==n) {
	return myMuonsRECO->at(imu).phi;;
      }
    }
  }*/
  return 0;
}

float EventCalculator::recoMuonIsoOfN(unsigned int n, const float ptthreshold) {
  //FIXME CFA
/*
  unsigned int ngood=0;
  for (unsigned int imu=0; imu < myMuonsRECO->size(); ++imu) {
    if ( isGoodRecoMuon(imu,true,ptthreshold)) { //true means disable RelIso cut
      ngood++;
      if (ngood==n) {

      float reliso = (myMuonsRECO->at(imu).trackIso
		      + myMuonsRECO->at(imu).ecalIso 
		      + myMuonsRECO->at(imu).hcalIso)/myMuonsRECO->at(imu).pt;

	return reliso;
      }
    }
  }*/
  return 0;
}

//i want to know if this reco muon is nearby a jet
float EventCalculator::recoMuonMinDeltaPhiJetOfN(unsigned int n, const float ptthreshold) {
  //FIXME CFA


  float mindp=99;
/*
  unsigned int ngood=0;

  for (unsigned int i=0; i< myMuonsRECO->size(); i++) {
    
    if (isGoodRecoMuon(i,true,ptthreshold)) {//disable iso
      ++ngood;
      if(ngood==n){

	//get the minimum angle between the this muon and all the good jets
	for (unsigned int ijet=0; ijet<jets_AK5PF_pt->size(); ijet++) {
	  if (isGoodJet( ijet ) ){
	    
	    float dp =  getDeltaPhi( myMuonsRECO->at(i).phi , myJetsPF->at(ijet).phi);
	    if (dp<mindp) mindp=dp;

	  }
	}
	return mindp;
      }
    }
  }*/
  
  return mindp;
}


float EventCalculator::tauPtOfN(unsigned int n) {

 return 0; //FIXME CFA
/*
  unsigned int ngood=0;
  for (unsigned int i=0; i < taus_pt->size(); i++) {
    if(isGoodTau(i)){
      ngood++;
      if (ngood==n) return getTauPt(i);
    }
  }
  return 0;
*/
}

float EventCalculator::tauEtaOfN(unsigned int n) {
  return 99; //FIXME CFA
/*
  unsigned int ngood=0;
  for (unsigned int i=0; i < taus_pt->size(); i++) {
    if(isGoodTau(i)){
      ngood++;
      if (ngood==n) return taus_eta->at(i);
    }
  }
  return 0;
*/
}


float  EventCalculator::jetGenPtOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {

    if( isSampleRealData() ) return 0;

    bool pass=false;
    pass = isGoodJet(i);

    if (pass ) {
      ngood++;
      if (ngood==n) return jets_AK5PF_gen_pt->at(i);
    }
  }
  return 0;
}


float  EventCalculator::jetPtOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {

    bool pass=false;
    pass = isGoodJet(i);

    if (pass ) {
      ngood++;
      if (ngood==n) return  getJetPt(i);
    }
  }
  return 0;
}


float  EventCalculator::jetGenPhiOfN(unsigned int n) {
  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {

    if( isSampleRealData() ) return 0;

    bool pass=false;
    pass = isGoodJet(i);

    if (pass ) {
      ngood++;
      if (ngood==n) return jets_AK5PF_gen_phi->at(i);
    }
  }
  return 0;
}


float  EventCalculator::jetPhiOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {

    bool pass=false;
    pass = isGoodJet(i);

    if (pass ) {
      ngood++;
      if (ngood==n) return  jets_AK5PF_phi->at(i);
    }
  }
  return 0;
}

float  EventCalculator::jetGenEtaOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {

    if( isSampleRealData() ) return 0;

    bool pass=false;
    pass = isGoodJet(i);

    if (pass ) {
      ngood++;
      if (ngood==n) return jets_AK5PF_gen_eta->at(i);
    }
  }
  return 0;

}

float  EventCalculator::jetEtaOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {

    bool pass=false;
    pass = isGoodJet(i);

    if (pass ) {
      ngood++;
      if (ngood==n) return jets_AK5PF_eta->at(i);
    }
  }
  return 0;
}

float  EventCalculator::jetEnergyOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {

    bool pass=false;
    pass = isGoodJet(i);

    if (pass ) {
      ngood++;
      if (ngood==n) return  getJetEnergy(i);
    }
  }
  return 0;
}

int  EventCalculator::jetFlavorOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {

    bool pass=false;
    pass = isGoodJet(i);

    if (pass ) {
      ngood++;
      if (ngood==n) return  jets_AK5PF_partonFlavour->at(i);
    }
  }
  return 0;
}

float  EventCalculator::jetChargedHadronFracOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
    bool pass=false;
    pass = isGoodJet(i);

    if (pass ) {
      ++ngood;
      if (ngood==n) return jets_AK5PF_chgHadE->at(i) / (jets_AK5PF_energy->at(i) * jets_AK5PF_corrFactorRaw->at(i) );
    }
  }
  return 0;
}

int  EventCalculator::jetChargedHadronMultOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
    bool pass=false;
    pass = isGoodJet(i);
    if (pass ) {
      ++ngood;
      if (ngood==n) return  TMath::Nint(jets_AK5PF_chg_Mult->at(i));
    }
  }

  return 0;
}

int  EventCalculator::jetMuMultOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
    bool pass=false;
    pass = isGoodJet(i);
    if (pass ) {
      ++ngood;
      if (ngood==n) return  TMath::Nint(jets_AK5PF_mu_Mult->at(i));
    }
  }

  return 0;
}
int  EventCalculator::jetNeutralHadronMultOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
    bool pass=false;
    pass = isGoodJet(i);
    if (pass ) {
      ++ngood;
      if (ngood==n) return  TMath::Nint(jets_AK5PF_neutral_Mult->at(i));
    }
  }

  return 0;
}

float EventCalculator::bjetCSVOfN(unsigned int n) {
  
  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
    
    bool pass=false;
    pass = (isGoodJet30(i) && passBTagger(i));
    
    if (pass ) {
      ngood++;
      if (ngood==n) return  getJetCSV(i);
    }
  }
  return 0;
}


float EventCalculator::bjetPtOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {

    bool pass=false;
    pass = (isGoodJet30(i) && passBTagger(i));

    if (pass ) {
      ngood++;
      if (ngood==n) return  getJetPt(i);
    }
  }
  return 0;
}

float EventCalculator::bjetPhiOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {

    bool pass=false;
    pass = (isGoodJet30(i) && passBTagger(i));

    if (pass ) {
      ngood++;
      if (ngood==n) return  jets_AK5PF_phi->at(i);
    }
  }
  return 0;
}

float EventCalculator::bjetEtaOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {

    bool pass=false;
    pass = (isGoodJet30(i) && passBTagger(i));

    if (pass ) {
      ngood++;
      if (ngood==n) return  jets_AK5PF_eta->at(i);
    }
  }
  return 0;
}

float EventCalculator::bjetEnergyOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {

    bool pass=false;
    pass = (isGoodJet30(i) && passBTagger(i));

    if (pass ) {
      ngood++;
      if (ngood==n) return  getJetEnergy(i);
    }
  }
  return 0;
}

int  EventCalculator::bjetFlavorOfN(unsigned int n) {

  if (jets_AK5PF_partonFlavour==0) return 0; //FIXME CFA

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {

    bool pass=false;
    pass = (isGoodJet30(i) && passBTagger(i));

    if (pass ) {
      ngood++;
      if (ngood==n) return  jets_AK5PF_partonFlavour->at(i);
    }
  }
  return 0;
}


float  EventCalculator::bjetChargedHadronFracOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
    bool pass=false;
    pass = isGoodJet(i) && passBTagger(i);

    if (pass ) {
      ++ngood;
      if (ngood==n) return jets_AK5PF_chgHadE->at(i) / (jets_AK5PF_energy->at(i) * jets_AK5PF_corrFactorRaw->at(i) );
    }
  }
  return 0;
}

int  EventCalculator::bjetChargedHadronMultOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
    bool pass=false;
    pass = isGoodJet(i) && passBTagger(i);

    if (pass ) {
      ++ngood;
      if (ngood==n) return TMath::Nint(jets_AK5PF_chg_Mult->at(i));
    }
  }
  return 0;
}


//-- first two jets are from W.  third is b jet
void EventCalculator::calcCosHel( unsigned int j1i, unsigned int j2i, unsigned int j3i, float & wcoshel,float &tcoshel) {

  //bool verb(true) ;
  bool verb(false) ;
  if ( verb ) { printf( "\n" ) ; }

  if ( j1i == j2i || j1i == j3i || j2i == j3i) {
    tcoshel = -99;
    wcoshel = -99;
    return;
  }

  double   sumPx = getJetPx(j1i) + getJetPx(j2i) + getJetPx(j3i);
  double   sumPy = getJetPy(j1i) + getJetPy(j2i) + getJetPy(j3i);
  double   sumPz = getJetPz(j1i) + getJetPz(j2i) + getJetPz(j3i); 

  //--- ignore quark masses.
  double j1Elab = sqrt( pow( getJetPx(j1i), 2) + pow( getJetPy(j1i), 2) + pow( getJetPz(j1i), 2) ) ;
  double j2Elab = sqrt( pow( getJetPx(j2i), 2) + pow( getJetPy(j2i), 2) + pow( getJetPz(j2i), 2) ) ;
  double j3Elab = sqrt( pow( getJetPx(j3i), 2) + pow( getJetPy(j3i), 2) + pow( getJetPz(j3i), 2) ) ;

  double sumP2 = sumPx*sumPx + sumPy*sumPy + sumPz*sumPz ;
  double sumP = sqrt( sumP2 ) ;

  //--- assume the top mass when computing the energy.

  double sumE = sqrt( mtop_*mtop_ + sumP2 ) ;

  double betatop = sumP / sumE ;
  double gammatop = 1. / sqrt( 1. - betatop*betatop ) ;

  TLorentzVector j1p4lab( getJetPx(j1i), getJetPy(j1i), getJetPz(j1i), j1Elab );
  TLorentzVector j2p4lab( getJetPx(j2i), getJetPy(j2i), getJetPz(j2i), j2Elab );
  TLorentzVector j3p4lab( getJetPx(j3i), getJetPy(j3i), getJetPz(j3i), j3Elab );

  if ( verb ) printf( " calc_wcoshel: betatop = %5.3f,  gammatop = %6.4f\n", betatop, gammatop ) ;

  double betatopx = sumPx / sumE ;
  double betatopy = sumPy / sumE ;
  double betatopz = sumPz / sumE ;

  TVector3 betatopvec( betatopx, betatopy, betatopz ) ;
  TVector3 negbetatopvec = -1.0 * betatopvec ;

  TLorentzVector j1p4trf( j1p4lab ) ;
  TLorentzVector j2p4trf( j2p4lab ) ;
  TLorentzVector j3p4trf( j3p4lab ) ;
  if (verb) { 
    printf(" calc_wcoshel: j1 tlorentzvector before boost:  %6.1f,  %6.1f,  %6.1f,   %6.1f\n",
	   j1p4trf.Px(), j1p4trf.Py(), j1p4trf.Pz(), j1p4trf.E() ) ;
  }
  j1p4trf.Boost( negbetatopvec ) ;
  j2p4trf.Boost( negbetatopvec ) ;
  j3p4trf.Boost( negbetatopvec ) ;
  if (verb) { 
    printf(" calc_wcoshel: j1 tlorentzvector after  boost:  %6.1f,  %6.1f,  %6.1f,   %6.1f\n",
	   j1p4trf.Px(), j1p4trf.Py(), j1p4trf.Pz(), j1p4trf.E() ) ;
  }


  tcoshel = j3p4trf.Vect().Dot( betatopvec ) / ( betatopvec.Mag() * j3p4trf.Vect().Mag() ) ;
  if ( verb ) { printf(" calc_wcoshel: top cos hel = %6.3f\n", tcoshel ) ; }

  //-- now boost j1 and j2 into W frame from top RF.

  TVector3 j1j2p3trf = j1p4trf.Vect() + j2p4trf.Vect() ;

  //-- Use W mass to compute energy
  double j1j2Etrf = sqrt( mW_*mW_ + j1j2p3trf.Mag2() );
  TLorentzVector j1j2p4trf( j1j2p3trf, j1j2Etrf ) ;
  TVector3 betawvec = j1j2p4trf.BoostVector() ;
  TVector3 negbetawvec = -1.0 * betawvec ;

  TLorentzVector j1p4wrf( j1p4trf ) ;
  TLorentzVector j2p4wrf( j2p4trf ) ;
  j1p4wrf.Boost( negbetawvec ) ;
  j2p4wrf.Boost( negbetawvec ) ;
  wcoshel = j1p4wrf.Vect().Dot( betawvec ) / ( betawvec.Mag() * j1p4wrf.Vect().Mag() ) ;
  if (verb) {
    printf(" calc_wcoshel: j1p4 in W rf: %6.1f,  %6.1f,  %6.1f,   %6.1f\n",
	   j1p4wrf.Px(), j1p4wrf.Py(), j1p4wrf.Pz(), j1p4wrf.E() ) ;
    printf(" calc_wcoshel: j2p4 in W rf: %6.1f,  %6.1f,  %6.1f,   %6.1f\n",
	   j2p4wrf.Px(), j2p4wrf.Py(), j2p4wrf.Pz(), j2p4wrf.E() ) ;
    TLorentzVector j1j2wrf = j1p4wrf + j2p4wrf ;
    printf(" calc_wcoshel: j1j2 mass in W rf: %7.2f\n", j1j2wrf.M() ) ;
    printf(" calc_wcoshel: w cos hel = = %6.3f\n", wcoshel ) ;
  }

  //   cout<<"cosHel calculator (t,W): "<<tcoshel<<" "<<wcoshel<<endl;
}


void EventCalculator::calcTopDecayVariables(float & wmass, float & tmass, float & wcoshel, float & tcoshel) {

  // cout<<" == event =="<<endl;

  double bestM2j=-9999;//, bestM2j_j1pt=0, bestM2j_j2pt=0;
  double bestM3j=-9999;//, bestM3j_j3pt=0;
  wcoshel = -99;
  tcoshel = -99;
  if(jets_AK5PF_pt->size()){
    //adopting this code from Owen -- note the loop goes to the second to last jet only
    for (unsigned int j1i = 0; j1i < jets_AK5PF_pt->size() -1; j1i++) {
      if ( isGoodJet30(j1i)) { //owen is using pT>30 cut

	//use exactly the same logic as Owen, to avoid bugs
	if (passBTagger(j1i) ) continue; //veto b jets

	//note how owen does the loop indexing here
	for (unsigned int j2i =j1i+1; j2i<jets_AK5PF_pt->size(); j2i++) {
	  if ( isGoodJet10(j2i)) { //owen is using a pT>10 cut here!

	    if (isGoodJet30(j2i) && passBTagger(j2i)) continue; //veto b jets with >30 gev

	    double m2j = calc_mNj(j1i,j2i);
	    if ( fabs(m2j- mW_) < fabs(bestM2j - mW_) ) {

	      bestM2j = m2j;
	      // bestM2j_j1pt = getLooseJetPt(j1i);
	      //bestM2j_j2pt = getLooseJetPt(j2i);

	      for ( unsigned int j3i=0; j3i<jets_AK5PF_pt->size(); j3i++) {
		if (j3i==j1i || j3i==j2i) continue;

		if ( isGoodJet30(j3i) && passBTagger(j3i)) { //owen uses 30 GeV pT cut

		  double m3j = calc_mNj(j1i,j2i,j3i);

		  if ( fabs(m3j-mtop_) < fabs(bestM3j-mtop_) ) {
		    bestM3j=m3j;
		    calcCosHel(j1i,j2i,j3i,wcoshel,tcoshel); //update helicity angles
		    //bestM3j_j3pt = getLooseJetPt(j3i);
		    //  cout<<"New best!"<<endl;
		  }

		} //is j3 good b jet
	      } //loop over j3

	    } // compare with W mass

	  } //jet 2 is good
	} //j2i

      } //if jet is good

    } //j1i
  }

  wmass = bestM2j;
  tmass = bestM3j;
}

float EventCalculator::getMT_Wlep(const float ptthreshold) {

  float MT = -1.;
  float MT2;
  int nE = countEle(ptthreshold);
  int nM = countMu(ptthreshold);
  
  float myMET = getMET();
  float myMETphi = getMETphi();
  float myMETx = myMET * cos(myMETphi);
  float myMETy = myMET * sin(myMETphi);

  float myP=0., myPphi=0.;
  bool newMT = false;
  if(nE==1 && nM==0){
    newMT = true;
    myP = elePtOfN(1,ptthreshold);
    myPphi = elePhiOfN(1,ptthreshold); 
  }
  else if(nE==0 && nM==1){
    newMT=true;
    myP = muonPtOfN(1,ptthreshold);
    myPphi = muonPhiOfN(1,ptthreshold);
  }
  float px = myP * cos(myPphi);
  float py = myP * sin(myPphi); 
  
  if (newMT) {
    float Et = sqrt( px*px + py*py );
    MT2 = 2*( Et * myMET - (px*myMETx + py*myMETy) ); 
    if(MT2>0.) {MT=sqrt(MT2);}
    else {MT = -2.;}

    //    cout<<" --- "<<endl;
    //    cout <<"Owen MT = "<<MT<<endl;
    //    float mt_ewk  = 2*myP*myMET*(1 - cos( getDeltaPhi(myPphi,myMETphi)));
    //    float mt_ewk2 = 2*myP*myMET*(1 - cos( myPphi-myMETphi));
    //    if (mt_ewk>0) mt_ewk = sqrt(mt_ewk);
    //    if (mt_ewk2>0) mt_ewk2 = sqrt(mt_ewk2);
    //    cout<<" new MT  = "<<mt_ewk<<endl;
    //    cout<<" new MT2 = "<<mt_ewk2<<endl;
  }
  //else cout<<" == "<<endl;

  return MT;
}

//get the transverse mass of the b-jet nearest the MET
float EventCalculator::getMT_bMET() {

  const float ptThreshold = 30;

  float mindp = 1e9;
  int thejet = -99;

  //find b jet closest to MET
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
    
    if (isGoodJet(i,ptThreshold)) {
      if (passBTagger(i)) {
	double dp=  getDeltaPhi( jets_AK5PF_phi->at(i) , getMETphi());
	if (dp<mindp) {
	  mindp=dp;
	  thejet = i;
	}
      }
    }
  }

  float mt = -1;
  if (thejet>=0) {
    mt  = 2*getJetPt(thejet)*getMET()*(1 - cos( getDeltaPhi( jets_AK5PF_phi->at(thejet) , getMETphi())));
    if (mt>0) mt = sqrt(mt);
    else cout<<"MT_b is negative"<<endl;
  }

  return mt;

}

//the ignored jet will be left out of the MHT calculation
//default is -1 (meaning no jet will be ignored)
std::pair<float,float> EventCalculator::getJERAdjustedMHTxy(int ignoredJet) {
  //no difference from MHT calculation -- just return both pieces
  float mhtx=0;
  float mhty=0;

  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
    if (isGoodJetMHT( i ) && (int(i) != ignoredJet) ) {
      mhtx -= getJetPx(i);
      mhty -= getJetPy(i);
    }
  }

  //used to have taus here. 
  //but i think that is not appropriate if the jet collection does not have the taus removed
  
  return make_pair(mhtx,mhty);
}


void EventCalculator::hadronicTopFinder_DeltaR(float & mjjb1, float & mjjb2 , float & topPT1, float & topPT2) {
  mjjb1=-1;
  mjjb2=-1;

  topPT1 = -1;
  topPT2 = -1;

  if (nGoodJets30() < 6) return; 
  //this algorithm forms 2 bjj bjj pairs using the top two bdisc jets as the b's.
  //in principle it might act weird for an event without 2 jets with any CSV value but that doesn't happen much in 6 jet events
  // (there is a commented-out line to "bail out" in case there are not two jets with a CSV value, but I think it is not needed)
  //it is mandatory that the event have 6 jets, hence the cut above.

  //the "light" jets are any jets that are not the ones used as b's

  //the jets are grouped into bqq using a rather crude "Minimum DeltaR" algorithm (that I invented myself)

  //also the pT of each 'top candidate' is returned

  // ===== find the best 2 b jets =====
  set<pair< float, int> > jets_sorted_by_bdisc;
  for (unsigned int kk=0; kk<jets_AK5PF_pt->size(); kk++) {
    if (isGoodJet30(kk) ) { 
      float bdiscval = jets_AK5PF_btag_secVertexCombined->at(kk);
      jets_sorted_by_bdisc.insert(make_pair(bdiscval,kk));
    }
  }

  //copy the best 2 guys to a simpler variable
  int bindices[2];
  int iii=0;
  //  cout<<" ~~~ b jets disc values"<<endl;
  for (set<pair< float, int> >::reverse_iterator ib=jets_sorted_by_bdisc.rbegin(); ib!=jets_sorted_by_bdisc.rend(); ++ib) {
    //    cout<<ib->first<<"\t"<<ib->second<<endl;
    if (iii<2) bindices[iii] = ib->second;
    //should check this choice of 0
    //    if (iii==1 && ib->first < 0) return; //key -- bail out completely if there are not 2 b-like jets!
    ++iii;
  }

  float mindrlight[2]={1e9,1e9};

  int mindrlight_ii[2]={-1,-1};
  int mindrlight_jj[2]={-1,-1};

  int nlightjets=0;

  for ( unsigned int iteration = 0; iteration<2; iteration++) {
    //we need to do the loop twice

    //loop over light jets
    for ( int ii=0; ii<int(jets_AK5PF_pt->size()); ii++) {
      //on the 2nd iteration, skip jets that are already used
      if (iteration==1 && ( (ii == mindrlight_ii[0]) || (ii == mindrlight_jj[0]) )) continue;

      //don't require that these "light" jets fail the btagger, but don't use the 2 best b-jets
      if (isGoodJet30(ii) && ii!=bindices[0] && ii!=bindices[1]) {
	if(iteration==0) nlightjets++;
	
	//loop over other light jets in the event
	for ( int jj=0; jj<int(jets_AK5PF_pt->size()); jj++) {
	  if ( ii==jj) continue;
	  //on the 2nd iteration, skip jets that are already used
	  if (iteration==1 && ( (jj == mindrlight_ii[0]) || (jj == mindrlight_jj[0]) )) continue;
	  
	  if (isGoodJet30(jj) &&  jj!=bindices[0] && jj!=bindices[1]) { //same consideration as above
	    
	    float thisDeltaR = jmt::deltaR( jets_AK5PF_eta->at(ii), jets_AK5PF_phi->at(ii), jets_AK5PF_eta->at(jj), jets_AK5PF_phi->at(jj));
	    if (thisDeltaR < mindrlight[iteration]) {
	      mindrlight[iteration] = thisDeltaR;
	      mindrlight_ii[iteration]=ii;
	      mindrlight_jj[iteration]=jj;
	    }
	  }
	}
      }
    }
  }

  //should never happen
  if (mindrlight_ii[0]== -1 || mindrlight_jj[0]==-1 || mindrlight_ii[1]==-1 ||  mindrlight_jj[1]==-1) {
    cout<<"nlight = "<<nlightjets<<endl;
    return;
  }

  //now we have light jet pairs (mindrlight_ii[0],mindrlight_jj[0])
  // and                        (mindrlight_ii[1],mindrlight_jj[1])

  //make 3-vectors of these jets
  TVector3 jet_ii0(getJetPx(mindrlight_ii[0]),getJetPy(mindrlight_ii[0]),getJetPz(mindrlight_ii[0]));
  TVector3 jet_jj0(getJetPx(mindrlight_jj[0]),getJetPy(mindrlight_jj[0]),getJetPz(mindrlight_jj[0]));
  TVector3 jet_ii1(getJetPx(mindrlight_ii[1]),getJetPy(mindrlight_ii[1]),getJetPz(mindrlight_ii[1]));
  TVector3 jet_jj1(getJetPx(mindrlight_jj[1]),getJetPy(mindrlight_jj[1]),getJetPz(mindrlight_jj[1]));

  //need to add the 3-vectors of these pairs together to form Ws (no mass constraint...just 4 vector addition)
  TVector3 W1 = jet_ii0 + jet_jj0;
  TVector3 W2 = jet_ii1 + jet_jj1;

  //crude (and wrong?) way to see the pairs that are nearest each other
  float drSum1 = jmt::deltaR( jets_AK5PF_eta->at( bindices[0]), jets_AK5PF_phi->at(bindices[0]), W1.Eta(), W1.Phi()) 
    + jmt::deltaR( jets_AK5PF_eta->at( bindices[1]), jets_AK5PF_phi->at(bindices[1]), W2.Eta(), W2.Phi()) ;

  float drSum2 = jmt::deltaR( jets_AK5PF_eta->at( bindices[1]), jets_AK5PF_phi->at(bindices[1]), W1.Eta(), W1.Phi()) 
    + jmt::deltaR( jets_AK5PF_eta->at( bindices[0]), jets_AK5PF_phi->at(bindices[0]), W2.Eta(), W2.Phi()) ;

  TVector3 bjet1(getJetPx(bindices[0]),getJetPy(bindices[0]),getJetPz(bindices[0]));
  TVector3 bjet2(getJetPx(bindices[1]),getJetPy(bindices[1]),getJetPz(bindices[1]));
  if (drSum1<drSum2 ) {
    mjjb1=    calc_mNj(bindices[0], mindrlight_ii[0],mindrlight_jj[0] );
    mjjb2=    calc_mNj(bindices[1], mindrlight_ii[1],mindrlight_jj[1] );
    //assemble top candidates
    TVector3 top1 = bjet1+W1;
    TVector3 top2 = bjet2+W2;
    topPT1 = top1.Perp();
    topPT2 = top2.Perp();
  }
  else {
    mjjb1=    calc_mNj(bindices[1], mindrlight_ii[0],mindrlight_jj[0] );
    mjjb2=    calc_mNj(bindices[0], mindrlight_ii[1],mindrlight_jj[1] );
    TVector3 top1 = bjet2+W1;
    TVector3 top2 = bjet1+W2;
    topPT1 = top1.Perp();
    topPT2 = top2.Perp();
  }

  //could also return the W mass

  //sort by mjjb, keeping the association with topPT
  if (mjjb2>mjjb1) {
    swap(mjjb1,mjjb2);
    swap(topPT1,topPT2);
  }
  return;
}

void EventCalculator::printDecay() {
  cout<<" -- printDecay --"<<endl;
  for (unsigned int k = 0; k< mc_doc_id->size(); k++) {
    // cout<<k<<"\t"<<TMath::Nint(mc_doc_id->at(k))<<" "<<mc_doc_mass->at(k)<<" mom="<<TMath::Nint(mc_doc_mother_id->at(k))<<endl;
    if ( abs(TMath::Nint(mc_doc_mother_id->at(k))) == 1000006) {
      cout<<TMath::Nint(mc_doc_mother_id->at(k))<<"\t"<<mc_doc_eta->at(k)<<" "<<mc_doc_phi->at(k)<<" "<<mc_doc_pt->at(k)<<endl;
    }
  }

}

int EventCalculator::jjResonance_mcTruthCheck(int jj1, int jj2) {

  if (!sampleName_.Contains("sbottom8lnotaus")) return 0; //don't evaluate this for data or any MC except the RPV sbottom signal
  /*
Goal:
check if the two reco'd jets actually came from the resonance decay of the stop

return either 0,1,2

0 == these jets do not match a stop
1 == these jets match stop 1
2 == these jets match stop 2

  */

  //after some checks by eye, i don't trust the partonFlavour variable
  float jet1_eta = jets_AK5PF_eta->at( jj1);
  float jet1_phi = jets_AK5PF_phi->at( jj1);

  float jet2_eta = jets_AK5PF_eta->at( jj2);
  float jet2_phi = jets_AK5PF_phi->at( jj2);

  float stop1_eta[2]={-99,-99};
  float stop1_phi[2]={-99,-99};
  int stop1_i = 0;

  float stop2_eta[2]={-99,-99};
  float stop2_phi[2]={-99,-99};
  int stop2_i = 0;

  //loop over gen particles and find the (anti)stops
  for (unsigned int k = 0; k< mc_doc_id->size(); k++) {
    if (  TMath::Nint(mc_doc_mother_id->at(k)) == 1000006) { // found the stop
      if (stop1_i <2) {
	stop1_eta[stop1_i] = mc_doc_eta->at(k);
	stop1_phi[stop1_i] = mc_doc_phi->at(k);

      }
      else cout<<"Too many stop daughters!"<<endl;
      ++stop1_i;
    }
    else if (  TMath::Nint(mc_doc_mother_id->at(k)) == -1000006) { // found the anti-stop
      if (stop2_i <2) {
	stop2_eta[stop2_i] = mc_doc_eta->at(k);
	stop2_phi[stop2_i] = mc_doc_phi->at(k);
      }
      else cout<<"Too many antistop daughters!"<<endl;
      ++stop2_i;
    }
  }

  const float matchingRadius=0.3; //seems to work ok in 10 test events

  //for stop1, do both jets match?
  bool match1 
    =   jmt::deltaR(jet1_eta,jet1_phi , stop1_eta[0], stop1_phi[0])<matchingRadius 
    &&  jmt::deltaR(jet2_eta,jet2_phi , stop1_eta[1], stop1_phi[1])<matchingRadius;
  //or the other combo
  bool match2 
    =   jmt::deltaR(jet2_eta,jet2_phi , stop1_eta[0], stop1_phi[0])<matchingRadius 
    &&  jmt::deltaR(jet1_eta,jet1_phi , stop1_eta[1], stop1_phi[1])<matchingRadius;
  
  //now try stop2
  bool match3 
    =   jmt::deltaR(jet1_eta,jet1_phi , stop2_eta[0], stop2_phi[0])<matchingRadius 
    &&  jmt::deltaR(jet2_eta,jet2_phi , stop2_eta[1], stop2_phi[1])<matchingRadius;
  //or the other combo
  bool match4 
    =   jmt::deltaR(jet2_eta,jet2_phi , stop2_eta[0], stop2_phi[0])<matchingRadius 
    &&  jmt::deltaR(jet1_eta,jet1_phi , stop2_eta[1], stop2_phi[1])<matchingRadius;
  
  //debug
//   if (match1) cout<<" match 1 DR = "<<jmt::deltaR(jet1_eta,jet1_phi , stop1_eta[0], stop1_phi[0])<<"\t"<<jmt::deltaR(jet2_eta,jet2_phi , stop1_eta[1], stop1_phi[1])<<endl;
//   else if (match2)cout<<" match 2 DR = "<<jmt::deltaR(jet2_eta,jet2_phi , stop1_eta[0], stop1_phi[0])<<"\t"<<jmt::deltaR(jet1_eta,jet1_phi , stop1_eta[1], stop1_phi[1])<<endl;

//   if (match3) cout<<" match 3 DR = "<<jmt::deltaR(jet1_eta,jet1_phi , stop2_eta[0], stop2_phi[0])<<"\t"<<jmt::deltaR(jet2_eta,jet2_phi , stop2_eta[1], stop2_phi[1])<<endl;
//   else if (match4)cout<<" match 4 DR = "<<jmt::deltaR(jet2_eta,jet2_phi , stop2_eta[0], stop2_phi[0])<<"\t"<<jmt::deltaR(jet1_eta,jet1_phi , stop2_eta[1], stop2_phi[1])<<endl;

  if ( (match1 || match2) && (match3||match4)) cout<<"WAIT -- Double match!"<<endl;
  if (match1 ||match2) return 1;
  if (match3 ||match4) return 2;
  return 0;

}

unsigned int EventCalculator::getSUSYnb() {

  unsigned int SUSY_nb=0;

  //how to revive for cfA?
  // answer: a much less fancy algo that just looks at materal ids

  //testing on T1tttt it seems to get the right answer most of the time (not all of the time...)

  for (unsigned int k = 0; k<mc_doc_id->size(); k++) {
    //important to require status 3
    if ( abs(TMath::Nint(mc_doc_id->at(k)))==5 && TMath::Nint(mc_doc_status->at(k))) {
      unsigned int id_mom  = abs(TMath::Nint(mc_doc_mother_id->at(k)));
      unsigned int id_gmom  = abs(TMath::Nint(mc_doc_grandmother_id->at(k)));
      unsigned int id_ggmom  = abs(TMath::Nint(mc_doc_ggrandmother_id->at(k)));

      if ( (id_mom>= 1000001 && id_mom <=1000039) || (id_mom>= 2000001 && id_mom <=2000015)) {SUSY_nb++;}
      else if ( (id_gmom>= 1000001 && id_gmom <=1000039) || (id_gmom>= 2000001 && id_gmom <=2000015)) {SUSY_nb++;}
      else if ( (id_ggmom>= 1000001 && id_ggmom <=1000039) || (id_ggmom>= 2000001 && id_ggmom <=2000015)) {SUSY_nb++;}

    }
  }

  //cout<<"SUSY_nb = "<<SUSY_nb<<endl;
  return SUSY_nb;

}

int EventCalculator::findJetMatchGenTau() {
  const float dr = 0.3;

  //finds the *leading* jet that matches within dr of a gen tau
  //arguably one should take the _closest_ jet instead.....

  set<int> jet_indices;

  //loop over gen taus
  for (size_t kk = 0; kk < mc_taus_id->size(); kk++) {

    float tauphi = mc_taus_phi->at(kk);
    float taueta = mc_taus_eta->at(kk);

    //find the closest jet to this tau
    float mindr = 1e9;
    int index=-1;

    //is there a jet within dr?
    for (size_t jj = 0; jj < jets_AK5PF_pt->size(); jj++) {
      float jetphi = jets_AK5PF_phi->at(jj);
      float jeteta = jets_AK5PF_eta->at(jj);

      float thisdr = jmt::deltaR(taueta,tauphi,jeteta,jetphi);
      if ((thisdr < dr) && (thisdr<mindr)) {
	mindr = thisdr;
	index = jj;
	//	cout<<"\t\t"<<taueta<<" "<<tauphi<<" "<<mc_taus_pt->at(kk)<<"   ;    "<<jeteta<<" "<<jetphi<<" "<<jets_AK5PF_pt->at(jj)<<endl;
      }
    }

    if (index>=0) {
      jet_indices.insert(index);
    }

  }

  //jets should be sorted by pt, so lowest indices are highest pt
  if (jet_indices.empty()) return -1;

  set<int>::iterator the_first_one = jet_indices.begin();
  return *the_first_one;
}

/*
http://svnweb.cern.ch/world/wsvn/icfsusy/trunk/AnalysisV2/SUSYSignalScan/include/NLOTools.hh 
[svnweb.cern.ch]
http://svnweb.cern.ch/world/wsvn/icfsusy/trunk/AnalysisV2/SUSYSignalScan/src/common/NLOTools.cc 
[svnweb.cern.ch]
*/
SUSYProcess EventCalculator::getSUSYProcess(float & pt1, float & phi1, float & pt2, float & phi2) {

  int nfound=0;
  float pt[2],phi[2];

  bool verbose = false;

  int squarks = 0;
  int antisquarks = 0;
  int gluinos = 0;
  int sleptons = 0;
  int neutralinos = 0;
  int charginos = 0;
  int sbottoms = 0;
  int stops = 0;

  int antisbottoms = 0; //jmt addition

  SUSYProcess process = NotFound;

  //  for (std::vector<Event::GenObject>::const_iterator j = ev.GenParticles().begin();  j != ev.GenParticles().end(); ++j) {
  for (unsigned int k = 0; k< mc_doc_id->size(); k++) {

    //can't use this in cfA
    //    int motheridx = TMath::Nint(myGenParticles->at(k).firstMother);
    //    if (motheridx<0) continue;
    int motherId = abs(TMath::Nint( mc_doc_mother_id->at(k) ));

    if (motherId == 21 || motherId <=6  ) {

      bool foundit=false;

      int myid = TMath::Nint( mc_doc_id->at(k) );
      //select squarks
      if (( myid >= 1000001 && myid <= 1000004) ||
	  (myid >= 2000001 && myid <= 2000004) ) {
        squarks++; foundit=true;
      }
      //select antisquarks
      if( ( myid <= -1000001 && myid >= -1000004) ||
	  ( myid <= -2000001 && myid >= -2000004) ){
        antisquarks++; foundit=true;
      }
      if ( abs( myid ) == 1000005 || abs( myid ) == 2000005) {sbottoms++; foundit=true;}
      if ( myid == -1000005 || myid == -2000005) {antisbottoms++;  foundit=true;} //jmt -- tacked on
      if ( abs( myid ) == 1000006 || abs( myid ) == 2000006) {stops++; foundit=true;}

      //select gluinos
      if ( abs( myid) == 1000021 ) {gluinos++; foundit=true;}

      //select sleptons
      if ( (abs( myid) >= 1000011 && fabs(myid) <= 1000016) ||
	   ( abs(myid) >= 2000011 && abs(myid) <= 2000016)) {sleptons++; foundit=true;}
      //select neutralinos
      if ( abs( myid) == 1000022 || abs( myid) == 1000023 ||  abs(myid) == 1000025 || abs(myid) == 1000035 ||  abs(myid) == 1000045  ) {neutralinos++; foundit=true;}

      if ( abs(myid) == 1000024 || abs(myid) == 1000037  ) {charginos++; foundit=true;}

      if (foundit) {
	if (nfound<2) {
	  pt[nfound]=mc_doc_pt->at(k);
	  phi[nfound]=mc_doc_phi->at(k);
	}
	++nfound;
      }

    }
  }

  if((neutralinos + charginos) == 1 && gluinos == 1) process = ng;
  if((neutralinos + charginos) == 1 && (squarks + antisquarks) == 1) process = ns;
  if(neutralinos + charginos == 2) process = nn;
  if(sleptons == 2) process = ll;
  if(squarks ==1 && antisquarks == 1) process = sb;
  if(squarks == 2) process = ss;
  if(stops == 2 ) process = tb;
  if(sbottoms == 2) process = bb;
  if(gluinos == 2) process = gg;
  if(squarks + antisquarks + sbottoms + stops == 1 && gluinos == 1) process = sg;

  if (process == NotFound || nfound>2) verbose = true;
  //jmt -- We seem to have no cross section for neutralino+sbottom
  //also squark + sbotttom

  //adjust sbottom number to just have b~ with no bbar~
  int   n_sbottoms = sbottoms - antisbottoms;
  if ( antisbottoms == 1 && squarks == 1) process=sb; //treat this as squark-antisquark
  if ( n_sbottoms == 1 && squarks == 1) process=ss; //treat this as squark-squark
  if ( n_sbottoms == 1 && antisquarks == 1) process=sb; //treat this as squark-antisquark
  if ( sbottoms==1 && (neutralinos + charginos) == 1) process=ns; //treat this as n+LF

  if(verbose) cout << "************ Process = "  <<  process << endl;
  if(verbose == true)cout << " neutralinos " << neutralinos << "\n charginos " << charginos << "\n gluinos " << gluinos << "\n squarks " << squarks << "\n antisquarks " << antisquarks << "\n sleptons " << sleptons << "\n stops " << stops << "\n sbottoms " << sbottoms 
			  <<" = (b~ + bbar~) = "<<n_sbottoms <<" + "<<antisbottoms
			  << endl;
  if (verbose)  cout<<"\tnfound = "<<nfound<<" ; pt,phi = "<<pt[0] << " "<<phi[0] <<" ; "<<pt[1] << " "<<phi[1]<<endl;

  pt1=pt[0];
  phi1 = phi[0];

  pt2=pt[1];
  phi2 = phi[1];

  return process;
}


std::pair<int,int> EventCalculator::getSMSmasses() {
  assert(theScanType_==kSMS);

  //the string looks like this:
  //# model T1tttt_1100_50  3.782E-12 

  TString modelstring = (*model_params).c_str();
  TString thesubstring=  modelstring.Tokenize(" ")->At(2)->GetName();

  int m0=  TString(thesubstring.Tokenize("_")->At(1)->GetName()).Atoi();
  int m12 = TString(thesubstring.Tokenize("_")->At(2)->GetName()).Atoi();

  //  cout<<m0<<"\t"<<m12<<endl;

  return  make_pair(m0,m12);

}

double EventCalculator::checkPdfWeightSanity( double a) {

  //not sure what to do here now

  if ( std::isnan(a) ) {cout<<"Found PDF NaN!"<<endl; return 0;}
  //  if ( a<0 ) return -1;
  //  if ( a>10) return 1; //i'm not sure if this is a good idea at all, nor do i know what threshold to set

  return a;
}

double EventCalculator::getHLTMHTeffBNN(float offMET, float offHT, uint nElectrons, uint nMuons, double mindphin,
					double& effUp, double& effDown) {

  float eff=1;
  effUp =1; effDown =1;
  return eff; //FIXME CFA

/*
  bool isSingleE = (nElectrons==1 && nMuons==0);
  bool isSingleMu = (nElectrons==0 && nMuons==1);
  bool is0L = (nElectrons==0 && nMuons==0);

  std::vector<double> vlumi;
  vlumi.push_back(0.0063);//ht260_mht60
  vlumi.push_back(0.0407);//ht250_mht60_v2
  vlumi.push_back(0.1687);//ht250_mht60_v3
  vlumi.push_back(0.1364);//ht250_mht70
  vlumi.push_back(0.1005+0.4366);//ht300_mht75_v7
  vlumi.push_back(0.0044);//ht300_mht75_v8
  vlumi.push_back(0.2773);//ht300_mht80_v1
  vlumi.push_back(0.8313);//ht300_mht80_v2
  vlumi.push_back(0.6522);//ht300_mht90
  vlumi.push_back(1.4421);//ht350_mht95
  vlumi.push_back(0.8855);//ht350_mht110  
  double totallumi = 4.982; //this is exactly the sum of the above lumi numbers

  // BNNs:
  std::vector<double> pfmet;
  pfmet.push_back( (double) getMET() );

  double effn[100];

  // Get the efficiency value from the BNNs:

  if( offHT>400 && is0L && mindphin > 4 ){

    if ( isSampleRealData() ) {//if data, just use the curve corresponding to the trigger that was used
      ULong64_t runnumber = getRunNumber();
      double stddeveff = 0;
      if      (runnumber >= 160431 && runnumber <= 161204){
	eff = HT260_MHT60_v2_HT260_v2_160431_161204_sel1( pfmet );	
	//compute the standard deviation of the ensemble of curves
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT260_MHT60_v2_HT260_v2_160431_161204_sel1(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 161205 && runnumber <= 163268){
	eff = HT250_MHT60_v2_HT250_v2_161205_163268_sel1( pfmet ); 
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT250_MHT60_v2_HT250_v2_161205_163268_sel1(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 163269 && runnumber <= 164923){
	eff = HT250_MHT60_v3_HT250_v3_163269_164923_sel1( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT250_MHT60_v3_HT250_v3_163269_164923_sel1(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 164924 && runnumber <= 165921){
	eff = HT250_MHT70_v1_HT250_v4_164924_165921_sel1( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT250_MHT70_v1_HT250_v4_164924_165921_sel1(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 165922 && runnumber <= 166300){
	eff = HT300_MHT75_v7_HT300_v6_165922_166300_166374_166978_sel1( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT300_MHT75_v7_HT300_v6_165922_166300_166374_166978_sel1(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 166301 && runnumber <= 166373){  
	eff = HT300_MHT75_v8_HT300_v7_166301_166373_sel1( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT300_MHT75_v8_HT300_v7_166301_166373_sel1(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 166374 && runnumber <= 166978){
	eff = HT300_MHT75_v7_HT300_v6_165922_166300_166374_166978_sel1( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT300_MHT75_v7_HT300_v6_165922_166300_166374_166978_sel1(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 166979 && runnumber <= 170064){ 
	eff = HT300_MHT80_v1_HT300_v6_166979_170064_sel1( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT300_MHT80_v1_HT300_v6_166979_170064_sel1(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      //end of summer 11 result data			  
      else if (runnumber >= 170065 && runnumber <= 173211){  
	eff = HT300_MHT80_v2_HT300_v9_170065_173211_sel1( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT300_MHT80_v2_HT300_v9_170065_173211_sel1(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 173212 && runnumber <= 176544){
	eff = HT300_MHT90_v2_HT300_v9_173212_176544_sel1( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT300_MHT90_v2_HT300_v9_173212_176544_sel1(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 176545 && runnumber <= 178410){  
	eff = HT350_MHT90_v1_HT350_v8_176545_178410_sel1( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT350_MHT90_v1_HT350_v8_176545_178410_sel1(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 178411){
	eff = HT350_MHT110_v3_HT350_v11_178411_180252_sel1( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT350_MHT110_v3_HT350_v11_178411_180252_sel1(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else {cout<<"No trigger assigned for run = "<<runnumber<<endl; assert(0);}

      stddeveff = sqrt(stddeveff / 99);

      effUp = (eff+stddeveff > 1) ? 1 : eff + stddeveff;
      effDown = eff - stddeveff;
    }
    else{//for MC, use the lumi-weighted average of all trigger curves
      std::vector<double> trigeffs;
      std::vector<double> trigeffserr;
      double stddeveff = 0, efferr = 0;

      trigeffs.push_back(HT260_MHT60_v2_HT260_v2_160431_161204_sel1( pfmet ));
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT260_MHT60_v2_HT260_v2_160431_161204_sel1(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT250_MHT60_v2_HT250_v2_161205_163268_sel1( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT250_MHT60_v2_HT250_v2_161205_163268_sel1(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT250_MHT60_v3_HT250_v3_163269_164923_sel1( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT250_MHT60_v3_HT250_v3_163269_164923_sel1(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT250_MHT70_v1_HT250_v4_164924_165921_sel1( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT250_MHT70_v1_HT250_v4_164924_165921_sel1(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT300_MHT75_v7_HT300_v6_165922_166300_166374_166978_sel1( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT300_MHT75_v7_HT300_v6_165922_166300_166374_166978_sel1(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT300_MHT75_v8_HT300_v7_166301_166373_sel1( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT300_MHT75_v8_HT300_v7_166301_166373_sel1(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT300_MHT80_v1_HT300_v6_166979_170064_sel1( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT300_MHT80_v1_HT300_v6_166979_170064_sel1(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT300_MHT80_v2_HT300_v9_170065_173211_sel1( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT300_MHT80_v2_HT300_v9_170065_173211_sel1(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT300_MHT90_v2_HT300_v9_173212_176544_sel1( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT300_MHT90_v2_HT300_v9_173212_176544_sel1(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT350_MHT90_v1_HT350_v8_176545_178410_sel1( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT350_MHT90_v1_HT350_v8_176545_178410_sel1(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT350_MHT110_v3_HT350_v11_178411_180252_sel1( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT350_MHT110_v3_HT350_v11_178411_180252_sel1(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);


      eff = 0;
      for(uint i = 0; i < trigeffs.size(); ++i){
	//std::cout << "eff " << i << " = " << trigeffs.at(i) << ", err = " << trigeffserr.at(i) << std::endl;
	eff += trigeffs.at(i) * (vlumi.at(i) / totallumi);
	efferr += trigeffserr.at(i) * trigeffserr.at(i) * (vlumi.at(i)*vlumi.at(i) / totallumi /  totallumi );
      }
      efferr = sqrt(efferr);
      //std::cout << "final eff = " << eff << ", err = " << efferr << std::endl;

      effUp = (eff+efferr > 1) ? 1 : eff + efferr;
      effDown = eff - efferr;

    }
  }
  else if ( offHT>400 && is0L && mindphin <= 4){

    if ( isSampleRealData() ) {//if data, just use the curve corresponding to the trigger that was used
      ULong64_t runnumber = getRunNumber();
      double stddeveff = 0;
      if      (runnumber >= 160431 && runnumber <= 161204){
	eff = HT260_MHT60_v2_HT260_v2_160431_161204_sel2( pfmet );
	
	//compute the standard deviation of the ensemble of curves
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT260_MHT60_v2_HT260_v2_160431_161204_sel2(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 161205 && runnumber <= 163268){
	eff = HT250_MHT60_v2_HT250_v2_161205_163268_sel2( pfmet ); 
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT250_MHT60_v2_HT250_v2_161205_163268_sel2(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 163269 && runnumber <= 164923){
	eff = HT250_MHT60_v3_HT250_v3_163269_164923_sel2( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT250_MHT60_v3_HT250_v3_163269_164923_sel2(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 164924 && runnumber <= 165921){
	eff = HT250_MHT70_v1_HT250_v4_164924_165921_sel2( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT250_MHT70_v1_HT250_v4_164924_165921_sel2(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 165922 && runnumber <= 166300){
	eff = HT300_MHT75_v7_HT300_v6_165922_166300_166374_166978_sel2( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT300_MHT75_v7_HT300_v6_165922_166300_166374_166978_sel2(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 166301 && runnumber <= 166373){  
	eff = HT300_MHT75_v8_HT300_v7_166301_166373_sel2( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT300_MHT75_v8_HT300_v7_166301_166373_sel2(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 166374 && runnumber <= 166978){
	eff = HT300_MHT75_v7_HT300_v6_165922_166300_166374_166978_sel2( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT300_MHT75_v7_HT300_v6_165922_166300_166374_166978_sel2(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 166979 && runnumber <= 170064){ 
	eff = HT300_MHT80_v1_HT300_v6_166979_170064_sel2( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT300_MHT80_v1_HT300_v6_166979_170064_sel2(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      //end of summer 11 result data			  
      else if (runnumber >= 170065 && runnumber <= 173211){  
	eff = HT300_MHT80_v2_HT300_v9_170065_173211_sel2( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT300_MHT80_v2_HT300_v9_170065_173211_sel2(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 173212 && runnumber <= 176544){
	eff = HT300_MHT90_v2_HT300_v9_173212_176544_sel2( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT300_MHT90_v2_HT300_v9_173212_176544_sel2(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 176545 && runnumber <= 178410){  
	eff = HT350_MHT90_v1_HT350_v8_176545_178410_sel2( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT350_MHT90_v1_HT350_v8_176545_178410_sel2(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else if (runnumber >= 178411){
	eff = HT350_MHT110_v3_HT350_v11_178411_180252_sel2( pfmet );
	for(uint i = 0; i < 100; ++i){
	  effn[i] = HT350_MHT110_v3_HT350_v11_178411_180252_sel2(pfmet, i, i);
	  stddeveff += (effn[i] - eff)*(effn[i] - eff);
	}      
      }
      else {cout<<"No trigger assigned for run = "<<runnumber<<endl; assert(0);}

      stddeveff = sqrt(stddeveff / 99);

      effUp = (eff+stddeveff > 1) ? 1 : eff + stddeveff;
      effDown = eff - stddeveff;
    }
    else{//for MC, use the lumi-weighted average of all trigger curves
      std::vector<double> trigeffs;
      std::vector<double> trigeffserr;
      double stddeveff = 0, efferr = 0;

      trigeffs.push_back(HT260_MHT60_v2_HT260_v2_160431_161204_sel2( pfmet ));
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT260_MHT60_v2_HT260_v2_160431_161204_sel2(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);


      trigeffs.push_back(HT250_MHT60_v2_HT250_v2_161205_163268_sel2( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT250_MHT60_v2_HT250_v2_161205_163268_sel2(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT250_MHT60_v3_HT250_v3_163269_164923_sel2( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT250_MHT60_v3_HT250_v3_163269_164923_sel2(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT250_MHT70_v1_HT250_v4_164924_165921_sel2( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT250_MHT70_v1_HT250_v4_164924_165921_sel2(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT300_MHT75_v7_HT300_v6_165922_166300_166374_166978_sel2( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT300_MHT75_v7_HT300_v6_165922_166300_166374_166978_sel2(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT300_MHT75_v8_HT300_v7_166301_166373_sel2( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT300_MHT75_v8_HT300_v7_166301_166373_sel2(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT300_MHT80_v1_HT300_v6_166979_170064_sel2( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT300_MHT80_v1_HT300_v6_166979_170064_sel2(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT300_MHT80_v2_HT300_v9_170065_173211_sel2( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT300_MHT80_v2_HT300_v9_170065_173211_sel2(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT300_MHT90_v2_HT300_v9_173212_176544_sel2( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT300_MHT90_v2_HT300_v9_173212_176544_sel2(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT350_MHT90_v1_HT350_v8_176545_178410_sel2( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT350_MHT90_v1_HT350_v8_176545_178410_sel2(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

      trigeffs.push_back(HT350_MHT110_v3_HT350_v11_178411_180252_sel2( pfmet ));
      stddeveff = 0;
      for(uint i = 0; i < 100; ++i){
	effn[i] = HT350_MHT110_v3_HT350_v11_178411_180252_sel2(pfmet, i, i);
	stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
      }      
      stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);


      eff = 0;
      for(uint i = 0; i < trigeffs.size(); ++i){
	//std::cout << "eff " << i << " = " << trigeffs.at(i) << ", err = " << trigeffserr.at(i) << std::endl;
	eff += trigeffs.at(i) * (vlumi.at(i) / totallumi);
	efferr += trigeffserr.at(i) * trigeffserr.at(i) * (vlumi.at(i)*vlumi.at(i) / totallumi /  totallumi );
      }
      efferr = sqrt(efferr);
      //std::cout << "final eff = " << eff << ", err = " << efferr << std::endl;

      effUp = (eff+efferr > 1) ? 1 : eff + efferr;
      effDown = eff - efferr;

    }
  }

  //SingleE and SingleMu samples (separately)
  //for *both* data and MC, use the lumi-weighted average of all trigger curves
  //(we only have the curves for the last three triggers)
  else if( offHT>400 && isSingleE && mindphin > 4){ 

    std::vector<double> vlumi_sl;
    vlumi_sl.push_back(0.0063);//ht260_mht60
    vlumi_sl.push_back(0.0407+0.1687);//ht250_mht60
    vlumi_sl.push_back(0.1364);//ht250_mht70 
    vlumi_sl.push_back(0.1005+0.0044+0.4366);//ht300_mht75
    vlumi_sl.push_back(0.2773+0.8313);//ht300_mht80
    vlumi_sl.push_back(0.6522);//ht300_mht90
    vlumi_sl.push_back(1.4421);//ht350_mht95
    vlumi_sl.push_back(0.8855);//ht350_mht110  
    double totallumi_sl = 4.982; //this is exactly the sum of the above lumi numbers

    std::vector<double> trigeffs;
    std::vector<double> trigeffserr;
    double stddeveff = 0, efferr = 0;

    trigeffs.push_back(ElHad_300_75( pfmet ));
    stddeveff = 0;
    for(uint i = 0; i < 100; ++i){
      effn[i] = ElHad_300_75(pfmet, i, i);
      stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
    }      
    stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

    trigeffs.push_back(ElHad_300_75( pfmet ));
    stddeveff = 0;
    for(uint i = 0; i < 100; ++i){
      effn[i] = ElHad_300_75(pfmet, i, i);
      stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
    }      
    stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

    trigeffs.push_back(ElHad_300_75( pfmet ));
    stddeveff = 0;
    for(uint i = 0; i < 100; ++i){
      effn[i] = ElHad_300_75(pfmet, i, i);
      stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
    }      
    stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

    trigeffs.push_back(ElHad_300_75( pfmet ));
    stddeveff = 0;
    for(uint i = 0; i < 100; ++i){
      effn[i] = ElHad_300_75(pfmet, i, i);
      stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
    }      
    stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

    trigeffs.push_back(ElHad_300_80( pfmet ));
    stddeveff = 0;
    for(uint i = 0; i < 100; ++i){
      effn[i] = ElHad_300_80(pfmet, i, i);
      stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
    }      
    stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);
  
    trigeffs.push_back(HT300MHT90elNew( pfmet ));
    stddeveff = 0;
    for(uint i = 0; i < 100; ++i){
      effn[i] = HT300MHT90elNew(pfmet, i, i);
      stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
    }      
    stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

    trigeffs.push_back(HT350MHT90elNew( pfmet ));
    stddeveff = 0;
    for(uint i = 0; i < 100; ++i){
      effn[i] = HT350MHT90elNew(pfmet, i, i);
      stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
    }      
    stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

    trigeffs.push_back(HT350MHT110elNew( pfmet ));
    stddeveff = 0;
    for(uint i = 0; i < 100; ++i){
      effn[i] = HT350MHT110elNew(pfmet, i, i);
      stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
    }      
    stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

    eff = 0;
    for(uint i = 0; i < trigeffs.size(); ++i){
      //std::cout << "eff " << i << " = " << trigeffs.at(i) << ", err = " << trigeffserr.at(i) << std::endl;
      eff += trigeffs.at(i) * (vlumi_sl.at(i) / totallumi_sl);
      efferr += trigeffserr.at(i) * trigeffserr.at(i) * (vlumi_sl.at(i)*vlumi_sl.at(i) / totallumi_sl /  totallumi_sl );
    }
    efferr = sqrt(efferr);
    //std::cout << "final eff = " << eff << ", err = " << efferr << std::endl;
    
    effUp = (eff+efferr > 1) ? 1 : eff + efferr;
    effDown = eff - efferr;
    
    //if(offMET>=200 && offMET<250){ 
    //  eff = 0.996;
    //  
    //  effUp = eff + 0;//assume error in SL region negligible
    //  effDown = eff - 0;
    //}
    //else if(offMET>250){
    //  eff = 0.999;
    //
    //  effUp = eff + 0;//assume error in SL region negligible
    //  effDown = eff - 0;
    //}
  }

  else if( offHT>400 && isSingleMu && mindphin > 4){ 

    std::vector<double> vlumi_sl;
    vlumi_sl.push_back(0.0063);//ht260_mht60
    vlumi_sl.push_back(0.0407+0.1687);//ht250_mht60
    vlumi_sl.push_back(0.1364);//ht250_mht70 
    vlumi_sl.push_back(0.1005+0.0044+0.4366);//ht300_mht75
    vlumi_sl.push_back(0.2773+0.8313);//ht300_mht80
    vlumi_sl.push_back(0.6522);//ht300_mht90
    vlumi_sl.push_back(1.4421);//ht350_mht95
    vlumi_sl.push_back(0.8855);//ht350_mht110  
    double totallumi_sl = 4.982; //this is exactly the sum of the above lumi numbers

    std::vector<double> trigeffs;
    std::vector<double> trigeffserr;
    double stddeveff = 0, efferr = 0;
  
    trigeffs.push_back(MuHad_300_75( pfmet ));
    stddeveff = 0;
    for(uint i = 0; i < 100; ++i){
      effn[i] = MuHad_300_75(pfmet, i, i);
      stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
    }      
    stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

    trigeffs.push_back(MuHad_300_75( pfmet ));
    stddeveff = 0;
    for(uint i = 0; i < 100; ++i){
      effn[i] = MuHad_300_75(pfmet, i, i);
      stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
    }      
    stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

    trigeffs.push_back(MuHad_300_75( pfmet ));
    stddeveff = 0;
    for(uint i = 0; i < 100; ++i){
      effn[i] = MuHad_300_75(pfmet, i, i);
      stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
    }      
    stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

    trigeffs.push_back(MuHad_300_75( pfmet ));
    stddeveff = 0;
    for(uint i = 0; i < 100; ++i){
      effn[i] = MuHad_300_75(pfmet, i, i);
      stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
    }      
    stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

    trigeffs.push_back(MuHad_300_80( pfmet ));
    stddeveff = 0;
    for(uint i = 0; i < 100; ++i){
      effn[i] = MuHad_300_80(pfmet, i, i);
      stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
    }      
    stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

    trigeffs.push_back(HT300MHT90muNew( pfmet ));
    stddeveff = 0;
    for(uint i = 0; i < 100; ++i){
      effn[i] = HT300MHT90muNew(pfmet, i, i);
      stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
    }      
    stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

    trigeffs.push_back(HT350MHT90muNew( pfmet ));
    stddeveff = 0;
    for(uint i = 0; i < 100; ++i){
      effn[i] = HT350MHT90muNew(pfmet, i, i);
      stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
    }      
    stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

    trigeffs.push_back(HT350MHT110muNew( pfmet ));
    stddeveff = 0;
    for(uint i = 0; i < 100; ++i){
      effn[i] = HT350MHT110muNew(pfmet, i, i);
      stddeveff += (effn[i] - trigeffs.back())*(effn[i] - trigeffs.back());
    }      
    stddeveff = sqrt(stddeveff / 99); trigeffserr.push_back(stddeveff);

    eff = 0;
    for(uint i = 0; i < trigeffs.size(); ++i){
      //std::cout << "eff " << i << " = " << trigeffs.at(i) << ", err = " << trigeffserr.at(i) << std::endl;
      eff += trigeffs.at(i) * (vlumi_sl.at(i) / totallumi_sl);
      efferr += trigeffserr.at(i) * trigeffserr.at(i) * (vlumi_sl.at(i)*vlumi_sl.at(i) / totallumi_sl /  totallumi_sl );
    }
    efferr = sqrt(efferr);
    //std::cout << "final eff = " << eff << ", err = " << efferr << std::endl;
    
    effUp = (eff+efferr > 1) ? 1 : eff + efferr;
    effDown = eff - efferr;
    
    //if(offMET>=200 && offMET<250){ 
    //  eff = 0.996;
    //  
    //  effUp = eff + 0;//assume error in SL region negligible
    //  effDown = eff - 0;
    //}
    //else if(offMET>250){
    //  eff = 0.999;
    //
    //  effUp = eff + 0;//assume error in SL region negligible
    //  effDown = eff - 0;
    //}
  }



  return eff;
*/
}

Long64_t EventCalculator::getNEventsGenerated() {
  
  if (isSampleRealData()) return 1;
  
  //normalization is handled in a special way for scans
  if (theScanType_ == kSMS) return 1;

  /*
    in the past we have auto-calculated event weights based on the number of events in the parent TChain
    
    This scheme is very convenient but has two problems.
    
    (1) we cannot split the reducedTree production into multiple jobs. This is problematic for very large samples (Fall11 ttbar)
    (2) if a skim has already been applied to the parent ntuples, there is no way to know the total number of events before the skim.
    
    Problem 1 could be dealt with by counting events in TChains and storing in a histogram. These could then be combined at the drawReducedTrees
    stage (similar to the machinery already used for signal scans). However, this does not work for problem 2.
    
    The only solution I see for problem 2 is to hard-code the number of generated events for each sample.
    I don't like this solution, but I do not see an alternative.
    This solution then also covers problem 1, and has the advantage that no new machinery is needed in drawReducedTrees.
    
    This function is intended to fulfill this solution by returning the number of events in each sample's original (unskimmed) ntuples.
    This may vary from one round of cfA production to another due to failed jobs, so the =entire= directory string is needed.
    UPDATE -- I think I can use the UCSBxxxx code as a unique identifier...


    The cfA page that catelogs the samples has this information automatically filled in -- just MAKE SURE TO TAKE THE
    NUMBER FROM THE UNSKIMMED SAMPLE!
  */
  
  // ============= v63 samples ============
  // SKIM is applied, so this is now important

  //numbers come from 
  //http://cms2.physics.ucsb.edu/cgi-bin/cfA.pl?Institute=ALL&process=ALL&version=v63

  if (sampleName_.Contains("UCSB1318")) return 3789889 ;
  if (sampleName_.Contains("UCSB1319")) return 1703863 ;
  if (sampleName_.Contains("UCSB1341")) return 30461028 ;

  if (sampleName_.Contains("UCSB1311")) return 1964088 ;
  if (sampleName_.Contains("UCSB1347")) return 5935732 ;
  if (sampleName_.Contains("UCSB1312")) return 2000062 ;
  if (sampleName_.Contains("UCSB1317")) return 1025622 ;
  if (sampleName_.Contains("UCSB1336")) return 2947160 ;
  if (sampleName_.Contains("UCSB1314")) return 977586 ;
  if (sampleName_.Contains("UCSB1292")) return 5927300 ;
  if (sampleName_.Contains("UCSB1332")) return 5480000 ;
  if (sampleName_.Contains("UCSB1297")) return 3994848 ;
  if (sampleName_.Contains("UCSB1298")) return 3992760 ;
  if (sampleName_.Contains("UCSB1310")) return 3998563 ;
  
  if (sampleName_.Contains("UCSB1355")) return 7618593 ;
  if (sampleName_.Contains("UCSB1338")) return 1087272 ;
  if (sampleName_.Contains("UCSB1291")) return 13880079 ;
  if (sampleName_.Contains("UCSB1344")) return 259961 ;
  if (sampleName_.Contains("UCSB1306")) return 23777 ;
  if (sampleName_.Contains("UCSB1307")) return 497658 ;
  if (sampleName_.Contains("UCSB1345")) return 139974 ;
  if (sampleName_.Contains("UCSB1330")) return 1935072 ;
  if (sampleName_.Contains("UCSB1305")) return 493460 ;
  
  if (sampleName_.Contains("UCSB1322")) return 1634582 ;
  if (sampleName_.Contains("UCSB1328")) return 1698744 ;
  if (sampleName_.Contains("UCSB1339")) return 1647807 ;
  if (sampleName_.Contains("UCSB1331")) return 18393090 ;
  
  if (sampleName_.Contains("UCSB1337")) return 6800431 ;
  if (sampleName_.Contains("UCSB1326")) return 9996622 ;
  if (sampleName_.Contains("UCSB1329")) return 4416646 ;
  if (sampleName_.Contains("UCSB1327")) return 5066608 ;
  if (sampleName_.Contains("UCSB1324")) return 1006928 ;
  if (sampleName_.Contains("UCSB1323")) return 4053786 ;
  if (sampleName_.Contains("UCSB1288")) return 9799908 ;


  //  if (sampleName_=="sbottom8lnotaus-189-270") return 50000; //probably deprecated in favor of the names below

  if (sampleName_=="RPVSUSY_sbottom8lnotaus-185-250") return 50000;
  if (sampleName_=="RPVSUSY_sbottom8lnotaus-189-270") return 50000;
  if (sampleName_=="RPVSUSY_sbottom8lnotaus-217-300") return 50000;

  // ============= v65 samples ============
  // SKIM is applied, so this is now important

  //numbers come from 
  //http://cms2.physics.ucsb.edu/cgi-bin/cfA.pl?Institute=ALL&process=ALL&version=v65
  if (sampleName_.Contains("UCSB1403")) return 6923750 ; //TT madgraph
  if (sampleName_.Contains("UCSB1439")) return 21673270; //TT powheg
  if (sampleName_.Contains("UCSB1443")) return 1634582 ; //W 250-300
  if (sampleName_.Contains("UCSB1405")) return 1699486 ; //W 300-400
  if (sampleName_.Contains("UCSB1406")) return 1647807 ; //W 400-
  if (sampleName_.Contains("UCSB1464")) return 1647807 ; //W 400- (52X)
  if (sampleName_.Contains("UCSB1422")) return 139974 ;  //tbar s channel
  if (sampleName_.Contains("UCSB1421")) return 259961 ;  //t s channel
  if (sampleName_.Contains("UCSB1451")) return 493460 ;  //tbar tW channel
  if (sampleName_.Contains("UCSB1455")) return 497658 ;  //t tW channel
  if (sampleName_.Contains("UCSB1461")) return 1923628;	//tbar t channel
  if (sampleName_.Contains("UCSB1470")) return 23777; //t t channel
  if (sampleName_.Contains("UCSB1438")) return 5995944; // QCD 50-80
                                                        //QCD 80-120
  if (sampleName_.Contains("UCSB1477")) return 5755732; //QCD120-170
  if (sampleName_.Contains("UCSB1453")) return 5814398 ;//QCD170-300
  if (sampleName_.Contains("UCSB1425")) return 5927300 ; //QCD 300-470
  if (sampleName_.Contains("UCSB1446")) return 3994848 ; //QCD 470-600
  if (sampleName_.Contains("UCSB1447")) return 3992760 ; //QCD 600-800
  if (sampleName_.Contains("UCSB1459")) return 3998563 ; //QCD 800-1000
  if (sampleName_.Contains("UCSB1448")) return 1964088; //QCD 1000-1400
  if (sampleName_.Contains("UCSB1449")) return 2000062; //QCD 1400-1800
  if (sampleName_.Contains("UCSB1450")) return 977586 ; // QCD 1800-
  if (sampleName_.Contains("UCSB1469")) return 4416646; //Znn100-200
  if (sampleName_.Contains("UCSB1468")) return 5055885; //Znn200-400
  if (sampleName_.Contains("UCSB1391")) return 1006928 ; //Znn 400-
                                           
  if (sampleName_.Contains("UCSB1389")) return 3789889 ; //DY 200-400
  if (sampleName_.Contains("UCSB1390")) return 1703863 ; //DY 400-
  if (sampleName_.Contains("UCSB1388")) return 30459503 ;//DY inclusive
  if (sampleName_.Contains("UCSB1424")) return 9799908 ; //ZZ
  if (sampleName_.Contains("UCSB1458")) return 10000431; //WW
  if (sampleName_.Contains("UCSB1445")) return 10000283 ; //WZ


    // ============= v66 samples ============
  // SKIM is applied, so this is now important

  //numbers come from 
  //`http://cms2.physics.ucsb.edu/cgi-bin/cfA.pl?Institute=ALL&process=ALL&version=v66
  if (sampleName_.Contains("UCSB1489")) return 6923750 ; //TT madgraph
  if (sampleName_.Contains("UCSB1558")) return 21675970; //TT powheg
  if (sampleName_.Contains("UCSB1488")) return 1634582 ; //W 250-300
  if (sampleName_.Contains("UCSB1512")) return 1699486 ; //W 300-400
  if (sampleName_.Contains("UCSB1487")) return 1647807 ; //W 400-
  if (sampleName_.Contains("UCSB1587")) return 4971847 ; //W 400- (not sure what the difference is here)
  if (sampleName_.Contains("UCSB1533")) return 139974 ;  //tbar s channel
  if (sampleName_.Contains("UCSB1530")) return 259961 ;  //t s channel
  if (sampleName_.Contains("UCSB1534")) return 493460 ;  //tbar tW channel
  if (sampleName_.Contains("UCSB1532")) return 497658 ;  //t tW channel
  if (sampleName_.Contains("UCSB1562")) return 1935072;	//tbar t channel
  if (sampleName_.Contains("UCSB1531")) return 23777; //t t channel
//  if (sampleName_.Contains("UCSB1438")) return 5995944; // QCD 50-80
                                                        //QCD 80-120
  if (sampleName_.Contains("UCSB1513")) return 5755732; //QCD120-170
  if (sampleName_.Contains("UCSB1561")) return 5814398 ;//QCD170-300
  if (sampleName_.Contains("UCSB1560")) return 5927300 ; //QCD 300-470
  if (sampleName_.Contains("UCSB1515")) return 3994848 ; //QCD 470-600
  if (sampleName_.Contains("UCSB1516")) return 3992760 ; //QCD 600-800
  if (sampleName_.Contains("UCSB1559")) return 3998563 ; //QCD 800-1000
  if (sampleName_.Contains("UCSB1577")) return 1964088; //QCD 1000-1400
  if (sampleName_.Contains("UCSB1578")) return 2000062; //QCD 1400-1800
  if (sampleName_.Contains("UCSB1517")) return 977586 ; // QCD 1800-
  if (sampleName_.Contains("UCSB1585")) return 977586 ; // QCD 1800- (53X)
  if (sampleName_.Contains("UCSB1579")) return 4040980; //Znn50-100
  if (sampleName_.Contains("UCSB1525")) return 4416646; //Znn100-200
  if (sampleName_.Contains("UCSB1524")) return 5055885; //Znn200-400
  if (sampleName_.Contains("UCSB1523")) return 1006928 ; //Znn 400-
                                           
  if (sampleName_.Contains("UCSB1526")) return 3789889 ; //DY 200-400
  if (sampleName_.Contains("UCSB1535")) return 1703863 ; //DY 400-
//  if (sampleName_.Contains("UCSB1388")) return 30459503 ;//DY inclusive
  if (sampleName_.Contains("UCSB1551")) return 9799908 ; //ZZ
  if (sampleName_.Contains("UCSB1563")) return 10000431; //WW
  if (sampleName_.Contains("UCSB1552")) return 10000283 ; //WZ

  if (sampleName_.Contains("UCSB1571")) return 4246444 ; //TT MG full lep
  if (sampleName_.Contains("UCSB1575")) return 11229902 ; //TT MG semi lep
  if (sampleName_.Contains("UCSB1586")) return 10537444;//TT MG hadronic

  cout<<"[getNEventsGenerated] unknown sample "<<sampleName_<<endl;
  assert(0);
  return 1;

}

double EventCalculator::getCrossSection(){
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV

  //const double bf = 0.32442;

  //Drell Yan
  if (sampleName_.BeginsWith("DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph")) return 3503.71; //NNLO
  if (sampleName_.BeginsWith("DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph")) return 23.43; //from AN //19.73; //LO PREP
  if (sampleName_.BeginsWith("DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph")) return  3.36; //from AN // 2.826; //LO PREP


  //Pythia QCD - all LO PREP
  if (sampleName_.BeginsWith("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6")) return 0.737844;
  if (sampleName_.BeginsWith("QCD_Pt-120to170_TuneZ2star_8TeV_pythia6")) return 156293.3;
  if (sampleName_.BeginsWith("QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6")) return 0.03352235;
  if (sampleName_.BeginsWith("QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6")) return 2.99815997E10; //correct thing to return?
  if (sampleName_.BeginsWith("QCD_Pt-170to300_TuneZ2star_8TeV_pythia6")) return 34138.15;
  if (sampleName_.BeginsWith("QCD_Pt-1800_TuneZ2star_8TeV_pythia6")) return 0.001829005;
  if (sampleName_.BeginsWith("QCD_Pt-300to470_TuneZ2star_8TeV_pythia6")) return 1759.549;
  if (sampleName_.BeginsWith("QCD_Pt-30to50_TuneZ2star_8TeV_pythia6")) return 6.6285328E7;
  if (sampleName_.BeginsWith("QCD_Pt-470to600_TuneZ2star_8TeV_pythia6")) return 113.8791;
  if (sampleName_.BeginsWith("QCD_Pt-600to800_TuneZ2star_8TeV_pythia6")) return 26.9921;
  if (sampleName_.BeginsWith("QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6")) return 3.550036;
  if (sampleName_.BeginsWith("QCD_Pt-80to120_TuneZ2star_8TeV_pythia6")) return 1033680.0;

  //single top
  if (sampleName_.BeginsWith("Tbar_t-channel_TuneZ2star_8TeV-powheg")) return 30.7; //from AN //25;//LO PREP 
  if (sampleName_.BeginsWith("Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg")) return 11.1;//from AN //10.7;//LO PREP
  if (sampleName_.BeginsWith("T_t-channel_TuneZ2star_8TeV-powheg")) return 56.4;// from AN //47.0;//LO PREP
  if (sampleName_.BeginsWith("T_tW-channel-DR_TuneZ2star_8TeV-powheg")) return 11.1; //from AN //10.7;//LO PREP
  if (sampleName_.BeginsWith("T_s-channel_TuneZ2star_8TeV-powheg-tauola")) return 3.79; // from AN //2.82; //LO PREP
  if (sampleName_.BeginsWith("Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola")) return  1.76;//from AN //1.57; //LO PREP

  //W+Jets
  if (sampleName_.BeginsWith("WJetsToLNu_HT-250To300_8TeV-madgraph")) return 57.3; //from AN //48.01; //LO PREP
  if (sampleName_.BeginsWith("WJetsToLNu_HT-300To400_8TeV-madgraph")) return 45.7; //from AN //38.3; //LO PREP
  if (sampleName_.BeginsWith("WJetsToLNu_HT-400ToInf_8TeV-madgraph")) return 30.1; //from AN //25.22; //LO PREP
  if (sampleName_.BeginsWith("WJetsToLNu_TuneZ2Star_8TeV-madgraph")) return 36257.2; //NNLO
  //diboson
  if (sampleName_.BeginsWith("WZ_TuneZ2star_8TeV_pythia6_tauola")) return 32.3;// from AN //12.63; //LO PREP
  if (sampleName_.BeginsWith("WW_TuneZ2star_8TeV_pythia6_tauola")) return 55; //from AN //33.61; //LO PREP
  if (sampleName_.BeginsWith("ZZ_TuneZ2star_8TeV_pythia6_tauola")) return 5.196; //LO PREP

  //Z -> nu nu  
  if (sampleName_.BeginsWith("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph")) return 160.3 *1.19 ; //LO PREP times DY k-factor
  if (sampleName_.BeginsWith("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph")) return 41.49 *1.19; //LO PREP times DY k-factor
  if (sampleName_.BeginsWith("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph")) return 5.274 *1.19; //LO PREP times DY k-factor
  if (sampleName_.BeginsWith("ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph")) return 381.2  *1.19; //LO PREP times DY k-factor

  //LM
  if (sampleName_.BeginsWith("SUSY_LM9_sftsht_8TeV"))          return 9.287; //LO 8TeV from PREP

  //ttbar
  if (sampleName_.BeginsWith("TTJets_TuneZ2star_8TeV-madgraph-tauola") )                return 234;  //approx NNLO
  if (sampleName_.BeginsWith("TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola")) return 234;  //approx NNLO
  if (sampleName_.BeginsWith("TT_CT10_TuneZ2star_8TeV-powheg")) return 234;  //approx NNLO

  if (sampleName_.BeginsWith("TTTo2L2Nu2B_8TeV-powheg-pythia6_Summer12")) return 22.14 *(234.0 / 136.3); //PREP corrected to NNLO

  if (sampleName_.BeginsWith("TTJets_FullLeptMGDecays_8TeV-madgraph")) return 13.43*(234.0 / (13.43+53.4+53.2)); // PREP corrected to NNLO
  if (sampleName_.BeginsWith("TTJets_HadronicMGDecays_8TeV-madgraph")) return 53.4*(234.0 / (13.43+53.4+53.2)); // PREP corrected to NNLO
  if (sampleName_.BeginsWith("TTJets_SemiLeptMGDecays_8TeV-madgraph")) return 53.2*(234.0 / (13.43+53.4+53.2)); //PREP corrected to NNLO

  //RPV 8 TeV -- cross-sections from 1206.2353. cross section is a function of the larger (sbottom) mass
  if (sampleName_.Contains("sbottom8lnotaus-185-250")) return 5.7; //LO madgraph
  if (sampleName_.Contains("sbottom8lnotaus-189-270")) return 3.7; //LO madgraph
  if (sampleName_.Contains("sbottom8lnotaus-217-300")) return 2.0; //LO madgraph


  if (sampleName_.Contains("mSUGRA")) return 1; //NLO cross sections will be specially stored per point
  if (sampleName_.Contains("T1bbbb")) return 1;
  if (sampleName_.Contains("T1tttt")) return 1;
  if (sampleName_.Contains("T2bb")) return 1;
  if (sampleName_.Contains("T2tt")) return 1;

  std::cout<<"Cannot find cross section for this sample!"<<std::endl;
  assert(0); 
  return -1;
}

TString EventCalculator::getSampleNameOutputString(){

  //strategy: as much as possible, give the name that drawReducedTrees expects,
  //and return sampleName for samples that have to be 'hadd'ed afterwards anyway

  //part of me likes this because it makes the filenames clean. but it can be dangerous if there are multiple samples with small differences (like only a difference in PU scenario)

  //8TeV

  //Drell Yan
  if (sampleName_.BeginsWith("DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph")) return "ZJets";
  if (sampleName_.BeginsWith("DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph")) return "ZJets_HT200To400";
  if (sampleName_.BeginsWith("DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph")) return "ZJets_HT400ToInf";


  //Pythia QCD - all LO PREP
  if (sampleName_.BeginsWith("QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6")) return "QCD1000";
  if (sampleName_.BeginsWith("QCD_Pt-120to170_TuneZ2star_8TeV_pythia6")) return "QCD120";
  if (sampleName_.BeginsWith("QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6")) return "QCD1400";
  if (sampleName_.BeginsWith("QCD_Pt-170to300_TuneZ2star_8TeV_pythia6")) return "QCD170";
  if (sampleName_.BeginsWith("QCD_Pt-1800_TuneZ2star_8TeV_pythia6")) return "QCD1800";
  if (sampleName_.BeginsWith("QCD_Pt-300to470_TuneZ2star_8TeV_pythia6")) return "QCD300";
  if (sampleName_.BeginsWith("QCD_Pt-30to50_TuneZ2star_8TeV_pythia6")) return "QCD30";
  if (sampleName_.BeginsWith("QCD_Pt-470to600_TuneZ2star_8TeV_pythia6")) return "QCD470";
  if (sampleName_.BeginsWith("QCD_Pt-600to800_TuneZ2star_8TeV_pythia6")) return "QCD600"; 
  if (sampleName_.BeginsWith("QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6")) return "QCD800";
  if (sampleName_.BeginsWith("QCD_Pt-80to120_TuneZ2star_8TeV_pythia6")) return  "QCD80";

  //single top (don't bother)

  //W+Jets
  if (sampleName_.BeginsWith("WJetsToLNu_HT-400ToInf_8TeV-madgraph")) return "WJets_HT400ToInf";
  if (sampleName_.BeginsWith("WJetsToLNu_HT-250To300_8TeV-madgraph")) return "WJets_HT250To300";
  if (sampleName_.BeginsWith("WJetsToLNu_HT-300To400_8TeV-madgraph")) return "WJets_HT300To400";
  if (sampleName_.BeginsWith("WJetsToLNu_TuneZ2Star_8TeV-madgraph")) return "WJets";
  //diboson
  if (sampleName_.BeginsWith("WZ_TuneZ2star_8TeV_pythia6_tauola")) return "WZ";
  if (sampleName_.BeginsWith("WW_TuneZ2star_8TeV_pythia6_tauola")) return "WW";
  if (sampleName_.BeginsWith("ZZ_TuneZ2star_8TeV_pythia6_tauola")) return "ZZ";

  //Z -> nu nu  
  if (sampleName_.BeginsWith("ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph")) return "Zinvisible_HT100To200";
  if (sampleName_.BeginsWith("ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph")) return "Zinvisible_HT200To400";
  if (sampleName_.BeginsWith("ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph")) return "Zinvisible_HT400ToInf";
  if (sampleName_.BeginsWith("ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph")) return "Zinvisible_HT50To100";

  //LM
  if (sampleName_.BeginsWith("SUSY_LM9_sftsht_8TeV"))          return "LM9";

  //ttbar
  //don't do this -- there are already 2 ttbar samples with very similar names


  //if it isn't found, just use the full name 
  return sampleName_;

}

TString EventCalculator::getCutDescriptionString(){
  
  TString cut = "";
  cut += theBTaggerNames_[theBTaggerType_];
  cut += "_";

  cut += theJetNames_[theJetType_];
  cut += "_";

  cut += theJESNames_[theJESType_];
  cut += "_";

  cut += theJERNames_[theJERType_];
  cut += "_";

  cut += theMETNames_[theMETType_];
  cut += "_";

  cut += theMETuncNames_[theMETuncType_];
  cut += "_";

  cut += thePUuncNames_[thePUuncType_];
  cut += "_";

  cut += theBTagEffNames_[theBTagEffType_];
  cut += "_";

  cut += theHLTEffNames_[theHLTEffType_];

  return cut;
}


double EventCalculator::getScanCrossSection( SUSYProcess p, const TString & variation ) {
  
  //will need to code for tan beta changes too.
  //but for now we only care about tan beta of 40
  if (p==NotFound) return 0;


  //FIXME CFA    
  if (theScanType_==kmSugra) {
    /*
    std::pair<int,int> thispoint = make_pair(TMath::Nint(eventlhehelperextra_m0),TMath::Nint(eventlhehelperextra_m12));
    if (variation=="")   return (*crossSectionTanb40_10_)[thispoint][p];
    else if (variation=="Plus")   return (*crossSectionTanb40_20_)[thispoint][p];
    else if (variation=="Minus")   return (*crossSectionTanb40_05_)[thispoint][p];
    else {assert(0);}
    */
  }
  else if (theScanType_==kSMS) {
    assert(0);
  }


  return 0;
}

void EventCalculator::jjResonanceFinder(float & mjj1, float & mjj2, int & ngoodMC) {//simple first try
  mjj1=0;
  mjj2=0;
  ngoodMC=0;

  //find the good jets
  vector<unsigned int> gjets;
  for (unsigned int i=0; i < jets_AK5PF_pt->size(); ++i) {
    if (isGoodJet30(i) )  gjets.push_back(i);
  }

  if ( gjets.size() <4 ) return;

  //0,1  2,3
  //0,2  1,3
  //0,3  1,2

  float mindiff=1e9;

  int r1_index1=-1;
  int r1_index2=-1;
  int r2_index1=-1;
  int r2_index2=-1;

  //we only want to look at the lead 4 jets
  const unsigned int maxjets=4;
  unsigned int ii=0; //i took out the loop over this. it seems like it is not needed for the strict lead 4 jet case
  for (unsigned int jj=ii+1; jj<maxjets; jj++) {
    
    int otherjet1=-1; int otherjet2=-1;
    for (unsigned int kk=0; kk<maxjets; kk++) {
      if (kk!=ii && kk!=jj) {
	if (otherjet1 == -1) otherjet1=kk;
	else if (otherjet2==-1) otherjet2=kk;
      }
    }
    //    cout<<ii<<" "<<jj<<"\tOther jets "<<otherjet1<<" "<<otherjet2<<endl;      
    
    //get invariant mass of ii,jj and the otherjet pair 
    float this_mjj1 = calc_mNj(gjets.at(ii),gjets.at(jj));
    float this_mjj2 = calc_mNj(gjets.at(otherjet1),gjets.at(otherjet2));

    if ( fabs(this_mjj1-this_mjj2) < mindiff) {
      mindiff = fabs(this_mjj1-this_mjj2);
      mjj1=this_mjj1;
      mjj2=this_mjj2;

      r1_index1 = gjets.at(ii);
      r1_index2 = gjets.at(jj);
      r2_index1 = gjets.at(otherjet1);
      r2_index2 = gjets.at(otherjet2);
    }

  }

//   cout<<"\t1a) "<<jets_AK5PF_eta->at(r1_index1)<<" "<<jets_AK5PF_phi->at(r1_index1)<<" "<<jets_AK5PF_pt->at(r1_index1)<<endl;
//   cout<<"\t1b) "<<jets_AK5PF_eta->at(r1_index2)<<" "<<jets_AK5PF_phi->at(r1_index2)<<" "<<jets_AK5PF_pt->at(r1_index2)<<endl;

//   cout<<"\t2a) "<<jets_AK5PF_eta->at(r2_index1)<<" "<<jets_AK5PF_phi->at(r2_index1)<<" "<<jets_AK5PF_pt->at(r2_index1)<<endl;
//   cout<<"\t2b) "<<jets_AK5PF_eta->at(r2_index2)<<" "<<jets_AK5PF_phi->at(r2_index2)<<" "<<jets_AK5PF_pt->at(r2_index2)<<endl;

  //now check if we found the right jets (in MC)
  int stop1index=  jjResonance_mcTruthCheck(r1_index1,r1_index2);
  int stop2index=  jjResonance_mcTruthCheck(r2_index1,r2_index2);
  //cout<<"[jjResFinder] "<<stop1index<< " "<<stop2index<<endl;
  if ( (stop1index ==stop2index) && stop1index>0) cout<<"WAIT -- we properly reco'd the same stop twice?"<<endl;
  if (stop1index>0) ++ngoodMC;
  if (stop2index>0) ++ngoodMC;

}

std::vector<unsigned int> EventCalculator::jetsetToVector(const vector<unsigned int> & goodjets, const set<unsigned int> & myset) {

  vector<unsigned int> outvec;

  for (  set<unsigned int>::iterator it=myset.begin() ; it != myset.end(); ++it ) {
    outvec.push_back( goodjets.at(*it));
  }

  return outvec;
}

void EventCalculator::jjResonanceFinder5(float & mjj1, float & mjj2) {
  //  cout<<"=="<<endl;
  mjj1=0;
  mjj2=0;

 //find the good jets                                                                                                                                                                               
  vector<unsigned int> gjets;
  for (unsigned int i=0; i < jets_AK5PF_pt->size(); ++i) {
    //for the moment don't bother with the 30 versus 25 GeV cuts. Just use 30
    if (isGoodJet30(i) )  gjets.push_back(i);
  }
  if (gjets.size() < 5) return;

  float mindiff=1e9;
  // " try all possible pairings of 2 and 3 jets"

  //map is overkill here; could just make this a set. but map works so let's stick with it
  map<  pair<set<unsigned int>, set<unsigned int> > , pair<float,float> > combos;

  for (unsigned int i1=0; i1<gjets.size(); i1++) {
    for (unsigned int i2=0; i2<gjets.size(); i2++) {
      for (unsigned int i3=0; i3<gjets.size(); i3++) {
	for (unsigned int i4=0; i4<gjets.size(); i4++) {
	  for (unsigned int i5=0; i5<gjets.size(); i5++) {
	    
	    //can only use a jet once
	    if ( i1==i2  || i1==i3 || i1==i4 || i1==i5) continue;
	    if (  i2==i3 || i2==i4 || i2==i5) continue;
	    if ( i3==i4 || i3==i5) continue;
	    if ( i4==i5) continue;

	    //use of set ensures that the indices are sorted, so that 1,3,5 and 5,3,1 are the same
	    //2-jet combo
	    set<unsigned int> r2;
	    r2.insert(i1);
	    r2.insert(i2);

	    //3-jet combo
	    set<unsigned int> r3;
	    r3.insert(i3);
	    r3.insert(i4);
	    r3.insert(i5);

	    pair<set<unsigned int>, set<unsigned int> > rpair = make_pair(r2,r3);
	    if (combos.count(rpair) == 0) { //have we already calculated this?
	      //get invariant masses
	      float this_mjj1= calc_mNj(jetsetToVector( gjets, r2));//gjets.at(i1),gjets.at(i2));
	      float this_mjj2= calc_mNj(jetsetToVector(gjets, r3));
	      //store inv masses
	      combos[rpair]=make_pair(mjj1,mjj2);
	      //is this the best one?
	      if ( fabs(this_mjj1-this_mjj2) < mindiff) {
		mindiff = fabs(this_mjj1-this_mjj2);
		mjj1=this_mjj1;
		mjj2=this_mjj2;
	      }
	      //	      cout<<"("<<i1<<" "<<i2<<") ("<<i3<<" "<<i4<<" "<<i5<<")\t"<<mjj1<<" "<<mjj2<<endl;
	    }

	  }
	}
      }
    }
  }

}


void EventCalculator::loadJetTagEffMaps() {

  assert(f_tageff_ ==0);

  TString filename = assembleBTagEffFilename(true); //pass true because we want to cut off the .XXXX in scan sample names
  filename.Prepend("btagEffMaps/");
  f_tageff_ = new TFile(filename,"READ");
  if (f_tageff_->IsZombie()) {
    cout<<"Failed to load the b-tag eff map for sample "<<sampleName_<<endl;
    delete f_tageff_;
    f_tageff_=0;
  }

}

void EventCalculator::calculateTagProb(float &Prob0, float &ProbGEQ1, float &Prob1, float &ProbGEQ2, float &Prob2, float &ProbGEQ3,
				       const float extraSFb, const float extraSFc, const float extraSFl, BTagEffModifier modifier) {

  //must initialize correctly
  Prob2 = 0;
  Prob1 = 0; ProbGEQ1 = 1; Prob0 = 1; ProbGEQ2 = 0;

  if (f_tageff_ == 0) { //if the b-tag eff file is not there make sure it will be obvious to the user -- fill all with zero
    Prob0 = 0;
    ProbGEQ1=0;
    Prob1=0;
    ProbGEQ2=0;
    Prob2=0;
    ProbGEQ3=0;
    return;
  }

  char btageffname[200], ctageffname[200], ltageffname[200];
  std::string sbtageff = "h_btageff";  std::string sctageff = "h_ctageff";  std::string sltageff = "h_ltageff";
  sprintf(btageffname,"%s",sbtageff.c_str());   
  sprintf(ctageffname,"%s",sctageff.c_str());   
  sprintf(ltageffname,"%s",sltageff.c_str());   
  TH1F * h_btageff  = (TH1F *)f_tageff_->Get(btageffname);
  TH1F * h_ctageff  = (TH1F *)f_tageff_->Get(ctageffname);
  TH1F * h_ltageff  = (TH1F *)f_tageff_->Get(ltageffname);

  for (unsigned int ijet=0; ijet<jets_AK5PF_pt->size(); ++ijet) {
    float subprob1=0;
    if(isGoodJet(ijet,50)){ //switch to 50 GeV threshold for b jets

      float effi = jetTagEff(ijet, h_btageff, h_ctageff, h_ltageff, extraSFb,extraSFc,extraSFl,modifier);
      //      cout<<"effi = "<<effi<<endl;
      Prob0 = Prob0* ( 1 - effi);
      
      double product = 1;
      for (unsigned int kjet=0; kjet<jets_AK5PF_pt->size(); ++kjet) {
	if(isGoodJet(kjet,50)){
	  float effk = jetTagEff(kjet, h_btageff, h_ctageff, h_ltageff,extraSFb,extraSFc,extraSFl,modifier);
	  if(kjet != ijet){
	    product = product*(1-effk);
	    //      cout<<"effk = "<<effk<<endl;
	  }
	  if(kjet > ijet){
	    double subproduct = 1;
	    for (unsigned int jjet=0; jjet<jets_AK5PF_pt->size(); ++jjet) {
	      //if(jjet > kjet && jjet > ijet){
	      if(jjet != kjet && jjet != ijet){
		if(isGoodJet(jjet,50)){
		  float effj = jetTagEff(jjet, h_btageff, h_ctageff, h_ltageff,extraSFb,extraSFc,extraSFl,modifier);
		  //      cout<<"effj = "<<effj<<endl;
		  subproduct = subproduct*(1-effj);
		}
	      }
	    }//j loop
	    subprob1 += effk*subproduct;
	  }

	}
      }//k loop
      
      Prob1 += effi*product;
      
      Prob2 += effi*subprob1;

    }
  }

  //  std::cout << "prob0 = " << Prob0 << ", prob1 = " << Prob1 << ", prob2 = " << Prob2 << std::endl;  

  ProbGEQ1 = 1 - Prob0;
  ProbGEQ2 = 1 - Prob1 - Prob0;
  ProbGEQ3 = 1 - Prob2 - Prob1 - Prob0;

}

int EventCalculator::getPtBinIndex(const float pt) {
  //binning in AN-12-175

  const int n=15;
  float bins[]={30,40,50,60,70,80,100,120,160,210,260,320,400,500,1e9}; //15 values
 //real max is 670, but use last bin value for jets about 670

  if (pt<bins[0] || pt>bins[n-1]) return -1;

  int jj=0; //highest value of jj in loop is 13
  for ( ; jj<n-1; ++jj)     if (pt >= bins[jj] && pt<bins[jj+1]) break;

  return jj;  
}

float EventCalculator::get_AN_12_175_Table2_Value(const float pt) {

  int bin = getPtBinIndex(pt);

  if (bin<0) return 1;

  float values[]={0.982,0.981,0.992,0.994,0.997,0.9998,1.001,1.000,0.992,0.979,0.947,0.928,0.87,0.84};

  return values[bin];
}

float EventCalculator::get_AN_12_175_Table2_Error(const float pt) {

  int bin = getPtBinIndex(pt);
  if (bin<0) return 1;

  float values[]={0.002,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.002,0.004,0.006,0.009,0.01,0.02};
  return values[bin];
}

float EventCalculator::get_AN_12_175_Table3_Value(const float pt) {

  int bin = getPtBinIndex(pt);
  if (bin<0) return 1;

  float values[]={0.980,0.984,0.992,0.995,0.999,1.002,0.999,0.999,0.994,0.978,0.943,0.92,0.88,0.86};
  return values[bin];
}

float EventCalculator::get_AN_12_175_Table4_Value(const float pt) {
  int bin = getPtBinIndex(pt);
  if (bin<0) return 1;

  //T1bbbb only!
  //(tempted to put in an assert, but i guess that is paranoia)

  float values[]={0.984,0.997,1.002,1.006,1.012,1.0115,1.0159,1.0171,1.0051,0.974,0.943,0.900,0.853,0.78};
  return values[bin];

}

float EventCalculator::get_AN_12_175_Table5_Value(const float pt) {
  int bin = getPtBinIndex(pt);
  if (bin<0) return 1;

  //T1tttt only!
  //(tempted to put in an assert, but i guess that is paranoia)

  float values[]={1.0046,1.0064,1.0136,1.0157,1.0173,1.0189,1.0224,1.0239,1.0087,0.981,0.945,0.899,0.858,0.800};
  return values[bin];

}

float EventCalculator::get_AN_12_175_Table6_Value(const float pt) {
  int bin = getPtBinIndex(pt);
  if (bin<0) return 1;

  float values[]={1.007,1.005,1.012,1.016,1.016,1.015,1.017,1.0191,1.0068,0.982,0.952,0.910,0.878,0.839};
  return values[bin];

}

float EventCalculator::get_AN_12_175_Table8_Value(const float pt) {
  int bin = getPtBinIndex(pt);
  if (bin<0) return 1;

  float values[]={0.9937,0.9945,1.0022,1.002,1.0031,1.0011,1.0006,1.0001,0.9861,0.967,0.946,0.914,0.884,0.840};
  return values[bin];

}


//removed unused code getBTagIPWeight()
float EventCalculator::bJetFastsimSF(const TString & what, int flavor,float pt) {
  assert( what == "value" || what=="stat" || what=="syst");
  //stat will be an absolute error
  //syst will be a fractional error

  float returnVal = 0;
  if (what == "value") returnVal=1;

  //first check if we're in a FASTSIM model
  if (theScanType_ != kNotScan ) {
    if ( abs(flavor) != 5) return returnVal; //no correction (for now) for non-b
    
    if (what == "syst") {//syst errors (fractional)
      const float table2_value = get_AN_12_175_Table2_Value(pt);

      //table 2 versus table 3 // note that this is a FRACTIONAL uncertainty
      float systerr2 = (table2_value - get_AN_12_175_Table3_Value(pt)) / table2_value ;

      //now the model-dependent part
      float model_value = 1;
      if (sampleName_.Contains("T1bbbb")) 	model_value = get_AN_12_175_Table4_Value(pt);
      else if (sampleName_.Contains("T1tttt")) 	model_value = get_AN_12_175_Table5_Value(pt);
      else if (sampleName_.Contains("T2bb")) 	model_value = get_AN_12_175_Table6_Value(pt);
      else if (sampleName_.Contains("T2tt")) 	model_value = get_AN_12_175_Table8_Value(pt);
      else assert(0);
      //Table 2 versus Table 4
      //Tom says syst_sample = (SF_avg - SF_ttbar)/SF_ttbar 
      //syst_sample = ( (SF_T1bbbb + SF_ttbar)/2 - SF_ttbar)/SF_ttbar 
      //note that this is a FRACTIONAL uncertainty

      float systerr1 = ( ( model_value+table2_value)/2 - table2_value) / table2_value ;

      //no need for fabs() on these errors because we're adding in quadrature
      returnVal = sqrt(systerr1*systerr1 + systerr2*systerr2);
    }
    //table 2 applies to any FASTSIM sample
    else if (what=="stat")    returnVal = get_AN_12_175_Table2_Error(pt);
    else if (what=="value")   returnVal = get_AN_12_175_Table2_Value(pt);
    else assert(0);
  }
  //cout<<what<<" "<<flavor<<" "<<pt<<"\t"<<returnVal<<endl;
  return returnVal;
}

//get MC btag efficiency
float EventCalculator::jetTagEff(unsigned int ijet, TH1F* h_btageff, TH1F* h_ctageff, TH1F* h_ltageff,
				 const float extraSFb, const float extraSFc, const float extraSFl,//options that were added but aren't used much
				 BTagEffModifier modifier ) { 

  float tageff=0;
  const float pt = getJetPt(ijet);
  const float x = pt>670 ? 670 : pt;
  const float eta = fabs(jets_AK5PF_eta->at(ijet));
  int flavor = jets_AK5PF_partonFlavour->at(ijet);

  float noEfficiencyThreshold=350; //threshold for the highest SF bin, which we sometimes set to 0 efficiency

  if(isGoodJet(ijet,50)) {

    if (theBTagEffType_ == kBTagEff04 || theBTagEffType_ == kBTagEffup4 || theBTagEffType_ == kBTagEffdown4) { //new  BTV-11-004 prescription 

      if (abs(flavor) ==4 || abs(flavor)==5) { //heavy flavor
	float errFactor = 1;

	// Tagger: CSVM within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
	if (pt >670) { //use SF for 670 with twice the errors
	  errFactor=2;
	}
	if (abs(flavor) == 4)  errFactor=2; //charm has double the errors   "SFc = SFb with twice the quoted uncertainty"

	//calculate sf using 'x' which is copy of pT with upper limit of 670
	float  SFb = 0.6981*((1.+(0.414063*x))/(1.+(0.300155*x)));

	//apply FASTSTIM correction where needed
	SFb *= bJetFastsimSF("value",flavor,pt);

	if (theBTagEffType_ == kBTagEffup4 || theBTagEffType_ == kBTagEffdown4 || modifier==kHFup || modifier == kHFdown) {
	  const int nbins=14;
	  float ptmin[] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500};
	  float ptmax[] = {40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 670};
	  float  SFb_error[] = {
	    0.0295675, //1
	    0.0295095,
	    0.0210867,
	    0.0219349,
	    0.0227033,
	    0.0204062,
	    0.0185857,
	    0.0256242,
	    0.0383341,
	    0.0409675,
	    0.0420284,
	    0.0541299,
	    0.0578761,
	    0.0655432 }; //14
	//need to figure out which bin i'm in
	  float myErr = SFb_error[nbins-1]; //put the high pT guys in the last bin
	  for (int i=0; i<nbins; i++) {
	    if (pt >= ptmin[i] && pt<ptmax[i]) { //found it
	      myErr = SFb_error[i]; 
	      break;
	    }
	  }

	  //fastsim corrections to the error
	  float fastsimerr_stat = bJetFastsimSF("stat",flavor,pt);
	  //syst is a fractional err, so multiply by SFb to get absolute err
	  float fastsimerr_syst = bJetFastsimSF("syst",flavor,pt) * SFb;
	  myErr = sqrt( myErr*myErr + fastsimerr_stat*fastsimerr_stat + fastsimerr_syst*fastsimerr_syst);

	  myErr *= errFactor; //high pT and charm get scaled up

	  // preliminary 2012 correction to the SFb and SFb uncertainty 
	  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC
	  myErr *= 1.5;

	  if ((theBTagEffType_ == kBTagEffup4) || (modifier==kHFup)) SFb += myErr;
	  else if ((theBTagEffType_ == kBTagEffdown4) || (modifier==kHFdown)) SFb -= myErr;
	  else assert(0);
	} //if syst variation

	//	cout<<"jet flavor, pt, SF = "<<abs(flavor)<<" "<<pt<<" "<<SFb<<endl;
	if      (abs(flavor) == 5) tageff = extraSFb * SFb * h_btageff->GetBinContent( h_btageff->FindBin( pt ) );
	else if (abs(flavor) == 4) tageff = extraSFc * SFb * h_ctageff->GetBinContent( h_ctageff->FindBin( pt ) );
	else assert(0);
      } // if heavy flavor
      else { //light flavor [ see https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs.C ]
	//note that these are valid only to 670 GeV, so use 'x' and not 'pt'
	float SF=0;
	float nominalSF=0;

	if ( eta < 0.8 ) {
	  nominalSF =  ((1.06182+(0.000617034*x))+(-1.5732e-06*(x*x)))+(3.02909e-10*(x*(x*x))); // SF without + or - variation in mistag rate
	  if       (theBTagEffType_ == kBTagEff04 && (modifier==kBTagModifier0 || modifier==kHFdown||modifier==kHFup))    SF = nominalSF;
	  else  if (theBTagEffType_ == kBTagEffdown4 || (modifier==kLFdown)) SF = ((0.972455+(7.51396e-06*x))+(4.91857e-07*(x*x)))+(-1.47661e-09*(x*(x*x)));
	  else  if (theBTagEffType_ == kBTagEffup4 || (modifier==kLFup))   SF = ((1.15116+(0.00122657*x))+(-3.63826e-06*(x*x)))+(2.08242e-09*(x*(x*x)));
	}
	else if (eta>=0.8 && eta<1.6) {
	  nominalSF = ((1.111+(-9.64191e-06*x))+(1.80811e-07*(x*x)))+(-5.44868e-10*(x*(x*x)));
	  if       (theBTagEffType_ == kBTagEff04 && (modifier==kBTagModifier0 || modifier==kHFdown||modifier==kHFup))    SF = nominalSF;
	  else  if (theBTagEffType_ == kBTagEffdown4 || (modifier==kLFdown)) SF = ((1.02055+(-0.000378856*x))+(1.49029e-06*(x*x)))+(-1.74966e-09*(x*(x*x)));
	  else  if (theBTagEffType_ == kBTagEffup4 || (modifier==kLFup))   SF = ((1.20146+(0.000359543*x))+(-1.12866e-06*(x*x)))+(6.59918e-10*(x*(x*x)));
	}
	else if (eta>=1.6 && eta<=2.4) {
	  nominalSF = ((1.08498+(-0.000701422*x))+(3.43612e-06*(x*x)))+(-4.11794e-09*(x*(x*x)));
	  if       (theBTagEffType_ == kBTagEff04 && (modifier==kBTagModifier0 || modifier==kHFdown||modifier==kHFup))    SF = nominalSF;
	  else  if (theBTagEffType_ == kBTagEffdown4 || (modifier==kLFdown)) SF = ((0.983476+(-0.000607242*x))+(3.17997e-06*(x*x)))+(-4.01242e-09*(x*(x*x)));
	  else  if (theBTagEffType_ == kBTagEffup4 || (modifier==kLFup))   SF = ((1.18654+(-0.000795808*x))+(3.69226e-06*(x*x)))+(-4.22347e-09*(x*(x*x)));
	}
	//design question -- what do to for jets at eta>2.4? assert? or return tageff=0?
	//i guess tageff 0 makes the most sense, so leave SF = 0

	//2012 -- modify the SF by a 2012/2011 factor given here: https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlighCorrFunc_2012_06_06.txt
	SF *= (1.10422 + (-0.000523856)*x + (1.14251e-06)*x*x);
	nominalSF *=  (1.10422 + (-0.000523856)*x + (1.14251e-06)*x*x);

	//for pT > 670, need to double the uncertainty
	if (pt>670) {
	  float deltaSF = SF - nominalSF; //here is the nominal uncertainty
	  deltaSF *= 2; //double it
	  SF = nominalSF + deltaSF;
	}

	tageff = SF * extraSFl * h_ltageff->GetBinContent( h_ltageff->FindBin( pt )); 
	//	cout<<"jet flavor, pt, SF = "<<abs(flavor)<<" "<<pt<<" "<<SF<<endl;

      } // if light flavor

    } //if BTV-11-004 prescription
    else {

      assert(modifier==kBTagModifier0);// not implementing this for the old options

    //CSVM
    float SF[3] = {0.96,0.96,0.96}; //3 bins of pT (<240, 240 to 350, and >350)
    float SFU[3] = {1,1,1};

    if (theBTagEffType_ == kBTagEffup) {
      SFU[0] = 1.10; // 10% uncertainty on SF from BTV-11-001
      SFU[1] = 1.20; // double the uncertainty for jets>240 GeV (need to check with B-POG)
      SFU[2] = 1.20; // double the uncertainty for jets>240 GeV (need to check with B-POG)
    }
    else if (theBTagEffType_ == kBTagEffdown) {
      SFU[0] = 0.9;
      SFU[1] = 0.8;
      SFU[2] = 0.8;
    }

    //Prescription used for Summer 2011 result (PAS)
    if (theBTagEffType_ == kBTagEff02 || theBTagEffType_ == kBTagEffup2 || theBTagEffType_ == kBTagEffdown2) {
      SF[2]=0; // assume no efficiency for pt>350 GeV
    }
    if (theBTagEffType_ == kBTagEffup2) {
      SFU[0] = 1.10; // 10% uncertainty on SF from BTV-11-001
      SFU[1] = 1.32; // 32% from UCSB's studies
      SFU[2] = 1.32; // this number is irrelevant
    }
    else if (theBTagEffType_ == kBTagEffdown2) {
      SFU[0] = 0.9;
      SFU[1] = 0.68;
      SFU[2] = 0.68;
    }

    //Fall 2011 preliminary prescription
    if (theBTagEffType_ == kBTagEff03 || theBTagEffType_ == kBTagEffup3 || theBTagEffType_ == kBTagEffdown3) {
      noEfficiencyThreshold=500;
      SF[2]=0; // assume no efficiency for pt>500 GeV
    }
    if (theBTagEffType_ == kBTagEffup3) {
      SFU[0] = 1.10; // 10% uncertainty on SF from BTV-11-001
      SFU[1] = 1.36; // 36% choice made on 25 Nov 2011
      SFU[2] = 1.36; // this number is irrelevant
    }
    else if (theBTagEffType_ == kBTagEffdown3) {
      SFU[0] = 0.9;
      SFU[1] = 0.64; 
      SFU[2] = 0.64;
    }

   
    //check the MC efficiencies from Summer11 TTbar
    //this histogram has bin edges {30,50,75,100,150,200,240,500,1000}
    
    if( abs(flavor) == 5){
      tageff = h_btageff->GetBinContent( h_btageff->FindBin( pt ) );
      
      if(pt<240) tageff *= SF[0]*SFU[0]*extraSFb;
      else if (pt>240 && pt<noEfficiencyThreshold) tageff *= SF[1]*SFU[1]*extraSFb;
      else tageff *= SF[2]*SFU[2]*extraSFb;

      //std::cout << "b: tag eff = " << tageff << std::endl;
    }
    else if (abs(flavor) == 4){
      tageff = extraSFc*h_ctageff->GetBinContent( h_ctageff->FindBin( pt ) );
      //std::cout << "c: tag eff = " << tageff << std::endl;
    }
    else if (abs(flavor) == 1 || abs(flavor) == 2
	     || abs(flavor) == 3 || abs(flavor) == 21){
      tageff = extraSFl*h_ltageff->GetBinContent( h_ltageff->FindBin( pt ) );
      //std::cout << "l: tag eff = " << tageff << std::endl;
    }
    }
    
  }

  return tageff;
}

void EventCalculator::dumpEvent () {

  //  if (!calledThisEvent_) {
  //    calledThisEvent_=true;    

    cout<<" == jets"<<endl;
    for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
      //    if (isCleanJet(i) ){//this only has an effect for recopfjets    
      cout<<"\t"<<jets_AK5PF_pt->at(i)<<" "<<jets_AK5PF_eta->at(i)<<" "<<jets_AK5PF_phi->at(i)<<"\t"<<isCleanJet(i)<<endl;
    }
    cout<<" == muons"<<endl;
    for ( unsigned int i = 0; i<pf_mus_et->size() ; i++) {
      cout<<"\t"<<pf_mus_pt->at(i)<<" "<<pf_mus_eta->at(i)<<" "<<pf_mus_phi->at(i)<<"\t"<<isCleanMuon(i)<<endl;
    }
    cout<<" == electrons"<<endl;
    for ( unsigned int i = 0; i<pf_els_et->size() ; i++) {
      cout<<"\t"<<pf_els_et->at(i)<<" "<<pf_els_eta->at(i)<<" "<<pf_els_phi->at(i)<<"\t"<<isGoodElectron(i)<<endl;
    }
    //  }

}


TString EventCalculator::getOptPiece(const TString &key, const TString & opt) {
  TObjArray * pieces = opt.Tokenize("_");

  for (int i=0; i<pieces->GetEntries(); i++) {
    TString thispiece = pieces->At(i)->GetName();
    if (thispiece.Contains(key)) {
      cout<<"Found piece = "<<thispiece<<endl;
      return thispiece;
    }
  }

  return "";
}

/* an idea to really study taus that I don't have time to execute
void EventCalculator::tauStudy(TString outputpath) {


  //open output file
  TString outfilename="tauTree";
  outfilename += ".";
  outfilename += getCutDescriptionString();
  outfilename += ".";

  outfilename += getSampleNameOutputString();
  outfilename+=".root";
  if (outputpath[outputpath.Length()-1] != '/') outputpath += "/";
  outfilename.Prepend(outputpath);
  TFile fout(outfilename,"RECREATE");

  // define the TTree
  TTree tauTree("tauTree","tree with minimal cuts");
  tauTree.SetMaxTreeSize(30000000000LL); //increase maximum tree size to 30GB

  //variables for tree
  double jetpt; //vector, actually


}
*/

void EventCalculator::reducedTree(TString outputpath) {

  //JSON file reading. For now the JSON file must be called GoldenJSON.txt
  //not clear to me if a sym link will be good enough w.r.t. the tarball, etc process
  // from /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-203002_8TeV_PromptReco_Collisions12_JSON_v2.txt
  vector< vector<int> > VRunLumi = MakeVRunLumi("GoldenJSON.txt");
  // from /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt
  vector< vector<int> > VRunLumiJuly13 = MakeVRunLumi("July13JSON.txt");
  // from /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt
  vector< vector<int> > VRunLumiAug24 = MakeVRunLumi("Aug24JSON.txt");


  //open output file
  TString outfilename="reducedTree";
  outfilename += ".";
  outfilename += getCutDescriptionString();
  outfilename += ".";

  outfilename += getSampleNameOutputString();
  outfilename+=".root";
  if (outputpath[outputpath.Length()-1] != '/') outputpath += "/";
  outfilename.Prepend(outputpath);
  cout<<"will save reducedTree to file: "<<outfilename<<endl;
  TFile fout(outfilename,"RECREATE");
  cout<<"will write to file "<<outfilename<<endl;
  if (fout.IsZombie()) assert(0); //save a lot of time in case of file problems

  // define the TTree
  TTree reducedTree("reducedTree","tree with minimal cuts");
  reducedTree.SetMaxTreeSize(30000000000LL); //increase maximum tree size to 30GB

  TRandom3 random(getSeed());

  // ~~~~~~~~~ declare reducedTree variables ~~~~~~~~~~
  //we're making an ntuple, so size matters -- use float not double		     

  //for cfA I am changing this name from 'weight' to 'eventweight' because weight is already a variable in the cfA ntuple.
  //oddly the compiler did not complain about this double-declaration.
  double eventweight; //one exception to the float rule
  double eventweight2; //one exception to the float rule
  //  float btagIPweight;//, pfmhtweight; //
  float PUweight;
  float PUweightSystVar;
  float hltHTeff;
  float hltMHTeff;
  //  float hltMHTeffBNN,hltMHTeffBNNUp,hltMHTeffBNNDown;

  double scanCrossSection,scanCrossSectionPlus,scanCrossSectionMinus;
  int m0,m12;

//   int W1decayType = -1, W2decayType = -1;
//   int decayType = -1;;
  int ttbarDecayCode=0;

  ULong64_t lumiSection, eventNumber, runNumber;
  //  float METsig;
  float ST, STeff, HT, HT30,MHT, MET, METphi, minDeltaPhi, minDeltaPhiAll, minDeltaPhiAll30,minDeltaPhi30_eta5_noIdAll;
  //  float correctedMET, correctedMETphi
  float caloMET,caloMETphi;
  float rawPFMET,rawPFMETphi;
  float maxDeltaPhi, maxDeltaPhiAll, maxDeltaPhiAll30, maxDeltaPhi30_eta5_noIdAll;
  float deltaPhi1, deltaPhi2, deltaPhi3;
  float sumDeltaPhi, diffDeltaPhi;
  float minDeltaPhiN, deltaPhiN1, deltaPhiN2, deltaPhiN3;
  float minDeltaPhiN_deltaT, deltaT1, deltaT2, deltaT3; 
  float minDeltaPhiN_otherEta5, minDeltaPhiN_otherEta5idNo, minDeltaPhiN_mainPt30_otherEta5idNo, minDeltaPhiN_mainPt30Eta5_otherEta5idNo, minDeltaPhiN_mainPt30Eta5idNo_otherEta5idNo;
  float minDeltaPhiK, minDeltaPhiK_otherEta5, minDeltaPhiK_otherEta5idNo, minDeltaPhiK_mainPt30_otherEta5idNo, minDeltaPhiK_mainPt30Eta5_otherEta5idNo, minDeltaPhiK_mainPt30Eta5idNo_otherEta5idNo;
  float minDeltaPhiN_DJR, minDeltaPhiN_DJR_otherEta5, minDeltaPhiK_DJR, minDeltaPhiK_DJR_otherEta5;
  float minDeltaPhiN_otherPt10, minDeltaPhiN_otherPt20, minDeltaPhiN_otherPt30, minDeltaPhiN_otherPt40, minDeltaPhiN_otherPt50;
  float minDeltaPhiN_DJR_otherPt10, minDeltaPhiN_DJR_otherPt20, minDeltaPhiN_DJR_otherPt30, minDeltaPhiN_DJR_otherPt40, minDeltaPhiN_DJR_otherPt50;
  UInt_t minDeltaPhi_chosenJet, minDeltaPhiN_chosenJet;
  UInt_t maxJetMis_chosenJet, maxJetFracMis_chosenJet;
  float minDeltaPhiN_withLepton;
  float minDeltaPhiN_asin, deltaPhiN_asin1, deltaPhiN_asin2, deltaPhiN_asin3;

  float maxJetMis, max2JetMis, maxJetMisAll30, max2JetMisAll30;
  float maxJetFracMis, max2JetFracMis, maxJetFracMisAll30, max2JetFracMisAll30;
  float deltaPhiMETJetMaxMis, deltaPhiMETJetMaxMis30;

  float deltaPhiStar, deltaPhiStar_badjet_pt, deltaPhiStar_badjet_eta, deltaPhiStar_badjet_phi;

  float minDeltaPhiMetTau;

  //items to investigate possible heavy flavor understimate in SIG -- looking at semi leptonic decays, bs, etc
  float CSVout1, CSVout2, CSVout3; //CSV tagger output for the lead three b-tagged jets
  float minDeltaPhiAllb30, deltaPhib1, deltaPhib2, deltaPhib3;
  float minDeltaPhiMETMuonsAll;

  bool cutPV;

  bool	csctighthaloFilter; // + 
  bool	eenoiseFilter; // - for 2012
  bool	greedymuonFilter; // + 
  bool	hbhenoiseFilter; // +
  // hcal laser
  // ee bad sc
  bool	inconsistentmuonFilter; // + 
  bool	ra2ecaltpFilter; // + 
  bool	scrapingvetoFilter; // +
  bool	trackingfailureFilter; // +
  bool badjetFilter; // new for 2012

  bool passCleaning;
  int PBNRcode; //being used in 2012

  bool buggyEvent;

  bool isRealData;

  //triggers
  bool cutTrigger;
  bool cutTrigger2;
  bool pass_utilityHLT;
  UInt_t prescaleUtilityHLT;
  UInt_t versionUtilityHLT;
  bool pass_utilityPrescaleModuleHLT;

  //defines the list of triggers to store info about
  map<TString, triggerData > triggerlist;
  //HTMHT
  triggerlist["PFHT350_PFMET100"]=triggerData();
  triggerlist["PFNoPUHT350_PFMET100"]=triggerData();
  triggerlist["HT250_AlphaT0p55"]=triggerData();
  triggerlist["HT300_AlphaT0p53"]=triggerData();
  //JetHT
  triggerlist["HT200"]=triggerData();
  triggerlist["HT250"]=triggerData();
  triggerlist["HT300"]=triggerData();
  triggerlist["PFHT350"]=triggerData();
  triggerlist["PFHT650"]=triggerData();
  triggerlist["PFNoPUHT350"]=triggerData();
  triggerlist["PFNoPUHT650"]=triggerData();
  //MET
  triggerlist["DiCentralPFJet50_PFMET80"]=triggerData();
  triggerlist["DiCentralPFJet30_PFMET80"]=triggerData();
  triggerlist["DiCentralPFJet30_PFMET80_BTagCSV07"]=triggerData();
  triggerlist["DiCentralPFNoPUJet50_PFMETORPFMETNoMu80"]=triggerData();
  triggerlist["PFMET150"]=triggerData();
  triggerlist["MET200"]=triggerData(); //calo MET
  triggerlist["MET120_HBHENoiseCleaned"]=triggerData(); //cal oMET
  triggerlist["L1ETM40"]; //prescaled L1 pass through
  //MuHad
  triggerlist["PFHT350_Mu15_PFMET45"]=triggerData();
  triggerlist["PFNoPUHT350_Mu15_PFMET45"]=triggerData();
  triggerlist["Mu8_DiJet30"]=triggerData(); //potentially useful for trigger studies
  //ElectronHad
  triggerlist["CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45"]=triggerData();
  triggerlist["CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45"]=triggerData();
  //SingleMuon
  triggerlist["IsoMu24_eta2p1"]=triggerData();
  triggerlist["IsoMu24"]=triggerData();
  //SinglePhoton
  triggerlist["Photon135"]=triggerData();
  triggerlist["Photon150"]=triggerData();
  //DoubleMu
  triggerlist["Mu17_Mu8"]=triggerData();
  //DoubleElectron
  triggerlist["Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL"]=triggerData();
  //MuEG
  triggerlist["Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL"]=triggerData();
  triggerlist["Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL"]=triggerData();
  //MultiJet
  triggerlist["SixJet45"]=triggerData();
  triggerlist["QuadJet80"]=triggerData();

  //want a copy of triggerlist, for storing mc trigger results
  map<TString, triggerData > triggerlist_mc(triggerlist);

  int njets, njets30, nbjets, nbjets30,ntruebjets, nElectrons, nMuons, nTaus;
  int nbjetsSSVM,nbjetsTCHET,nbjetsSSVHPT,nbjetsTCHPT,nbjetsTCHPM,nbjetsCSVM,nbjetsCSVL; 
  int nElectrons5, nElectrons15,nElectrons20;
  int nMuons5, nMuons15,nMuons20;

  float bestZmass;
  float mjj1,mjj2,mjjdiff;
  float mjj1_5,mjj2_5,mjjdiff_5;
  int nCorrectRecoStop;

  float jetpt1, jetgenpt1, jetphi1, jetgenphi1, jeteta1, jetgeneta1, jetenergy1, bjetpt1, bjetphi1, bjeteta1, bjetenergy1;
  float jetpt2, jetgenpt2, jetphi2, jetgenphi2, jeteta2, jetgeneta2, jetenergy2, bjetpt2, bjetphi2, bjeteta2, bjetenergy2;
  float jetpt3, jetgenpt3, jetphi3, jetgenphi3, jeteta3, jetgeneta3, jetenergy3, bjetpt3, bjetphi3, bjeteta3, bjetenergy3;
  int jetflavor1, jetflavor2, jetflavor3, bjetflavor1, bjetflavor2, bjetflavor3;

  //lead jet quality info
  float alletajetpt1,alletajetphi1,alletajeteta1;
  float alletajetneutralhadronfrac1,alletajetneutralemfrac1,alletajetphotonfrac1;


  float jetchargedhadronfrac1, jetchargedhadronfrac2, jetchargedhadronfrac3, bjetchargedhadronfrac1, bjetchargedhadronfrac2, bjetchargedhadronfrac3;
  int jetchargedhadronmult1, jetchargedhadronmult2, jetchargedhadronmult3, bjetchargedhadronmult1, bjetchargedhadronmult2, bjetchargedhadronmult3;
  int jetneutralhadronmult1, jetneutralhadronmult2,jetneutralhadronmult3;
  int jetmumult1, jetmumult2,jetmumult3;


  float mjjb1,mjjb2,topPT1,topPT2;

  float eleet1, elephi1, eleeta1, muonpt1, muonphi1, muoneta1;
  int elecharge1, muoncharge1;
  float eleet2, elephi2, eleeta2, muonpt2, muonphi2, muoneta2;
  float muoniso1,muonchhadiso1,muonphotoniso1,muonneutralhadiso1;
  float taupt1, taueta1;
  float eleRelIso,muonRelIso;

  float recomuonpt1, recomuonphi1, recomuoneta1;
  float recomuonpt2, recomuonphi2, recomuoneta2;
  float recomuoniso1, recomuoniso2;
  float recomuonmindphijet1, recomuonmindphijet2;

  float rl,rMET;

  float METsig,METsig00,METsig10,METsig11;

  float MT_b;
  float MT_Wlep;
  float MT_Wlep5,MT_Wlep15;
  float wMass, topMass, wCosHel, topCosHel;

  int nGoodPV;

  int SUSY_nb;
  int SUSY_process;
  float SUSY_recoilPt;

  //sphericity variables (updated)
  float transverseSphericity_jets;
  float transverseSphericity_jetsMet;
  float transverseSphericity_jetsMetLeptons;
  float transverseSphericity_jetsLeptons;

  float transverseSphericity_jets30;
  float transverseSphericity_jets30Met;
  float transverseSphericity_jets30MetLeptons;
  float transverseSphericity_jets30Leptons;

  float deltaR_bestTwoCSV;
  float mjj_bestTwoCSV;
  float deltaR_closestB;
  float mjj_closestB;

/*
  float transverseThrust,transverseThrustPhi;
  float transverseThrustWithMET,transverseThrustWithMETPhi;
*/
  float minDeltaPhiN_Luke, maxDeltaPhiN_Luke, deltaPhiN1_Luke, deltaPhiN2_Luke, deltaPhiN3_Luke;
  float minTransverseMETSignificance, maxTransverseMETSignificance, transverseMETSignificance1, transverseMETSignificance2, transverseMETSignificance3;

  int nLostJet;
  int njets_lostJet, nbjets_lostJet;
  float minDeltaPhiN_lostJet, deltaPhiN1_lostJet, deltaPhiN2_lostJet, deltaPhiN3_lostJet;
  float minDeltaPhiN_Luke_lostJet, maxDeltaPhiN_Luke_lostJet, deltaPhiN1_Luke_lostJet, deltaPhiN2_Luke_lostJet, deltaPhiN3_Luke_lostJet;
  float minTransverseMETSignificance_lostJet, maxTransverseMETSignificance_lostJet, transverseMETSignificance1_lostJet, transverseMETSignificance2_lostJet, transverseMETSignificance3_lostJet;

  bool passBadECAL_METphi53_n10_s12, passBadECAL_METphi33_n10_s12, passBadECAL_METphi52_n10_s12; 
  float worstMisA_badECAL_METphi53_n10_s12, worstMisF_badECAL_METphi53_n10_s12;
  float worstMisA_badECAL_METphi33_n10_s12, worstMisF_badECAL_METphi33_n10_s12;
  float worstMisA_badECAL_METphi52_n10_s12, worstMisF_badECAL_METphi52_n10_s12;
  bool passBadECAL_METphi53_n5_s12;
  bool passBadECAL_METphi53_crack;
  float worstMisA_badECAL_METphi53_crack, worstMisF_badECAL_METphi53_crack;
  bool passBadECAL_METphi52_crack;
  float worstMisA_badECAL_METphi52_crack, worstMisF_badECAL_METphi52_crack;
  bool passBadECAL_METphi53_crackF;
  float worstMisA_badECAL_METphi53_crackF, worstMisF_badECAL_METphi53_crackF;
  bool passBadECAL_METphi52_crackF;
  float worstMisA_badECAL_METphi52_crackF, worstMisF_badECAL_METphi52_crackF;

  int nTracks1,nTracks10;

  float prob0,probge1,prob1,probge2,probge3,prob2;
  
  float prob0_HFplus,probge1_HFplus,prob1_HFplus,probge2_HFplus,probge3_HFplus,prob2_HFplus;
  float prob0_HFminus,probge1_HFminus,prob1_HFminus,probge2_HFminus,probge3_HFminus,prob2_HFminus;
  float prob0_LFplus,probge1_LFplus,prob1_LFplus,probge2_LFplus,probge3_LFplus,prob2_LFplus;
  float prob0_LFminus,probge1_LFminus,prob1_LFminus,probge2_LFminus,probge3_LFminus,prob2_LFminus;

  std::vector<int> vrun,vlumi,vevent;
  loadEventList(vrun, vlumi, vevent);


  Float_t pdfWeightsCTEQ[45];
  Float_t pdfWeightsMSTW[41];
  Float_t pdfWeightsNNPDF[100];

  // bookkeeping for screen output only
  pair<int,int> lastpoint = make_pair(0,0);

  /*
scanProcessTotalsMap has been extended to handle PDF weight sums.

TH1[process] -> TH2[process, pdf index]

3 scanProcessTotalsMaps (one per pdf set)

should make a scanProcessTotalsMap for any kind of scan (mSugra or SMS)

~~~ things to deal with later ~~~
This scheme makes scanSMSngen redundant. But it is so convenient that i leave it in
Also the pdfWeightSum* histograms that are used for LM9.
  */

  //for scans, I want to keep track of how many of each susy process are generated at each point
  map<pair<int,int>, TH2D* >  scanProcessTotalsMapCTEQ;
  map<pair<int,int>, TH2D* >  scanProcessTotalsMapMSTW;
  map<pair<int,int>, TH2D* >  scanProcessTotalsMapNNPDF;
  if (!sampleIsSignal_) {} //don't create any of these histograms
  else if (theScanType_ == kmSugra) {// should be equiv to below
    //  if (crossSectionTanb40_10_ != 0) {
    for (map<pair<int, int>, map<SUSYProcess, double> >::iterator iscanpoint=crossSectionTanb40_10_->begin(); iscanpoint!=crossSectionTanb40_10_->end(); ++iscanpoint) {
      TString histoname = "scanProcessTotals"; 
      histoname += "_";
      histoname +=iscanpoint->first.first; // m0
      histoname += "_";
      histoname +=iscanpoint->first.second; // m12
      //      scanProcessTotalsMap[iscanpoint->first] = new TH1D(histoname,histoname,int(NotFound),int(ng),int(NotFound));

      //new 2D variations
      histoname.ReplaceAll("scanProcessTotals","scanProcessTotalsCTEQ");
      scanProcessTotalsMapCTEQ[iscanpoint->first] = new TH2D(histoname,histoname,int(NotFound),int(ng),int(NotFound),45,0,45);
      histoname.ReplaceAll("scanProcessTotalsCTEQ","scanProcessTotalsMSTW");
      scanProcessTotalsMapMSTW[iscanpoint->first] = new TH2D(histoname,histoname,int(NotFound),int(ng),int(NotFound),41,0,41);
      histoname.ReplaceAll("scanProcessTotalsMSTW","scanProcessTotalsNNPDF");
      scanProcessTotalsMapNNPDF[iscanpoint->first] = new TH2D(histoname,histoname,int(NotFound),int(ng),int(NotFound),100,0,100);
    }
  }
  else if (theScanType_==kSMS) {
    //will need a loop over scan points. probably i need to do the loop by hand.
    //for each scan point generate a histogram, as above
    for (int im0= 0; im0<= 2000; im0+=25) {
      for (int im12= 0; im12<= 2000; im12+=25) {
	TString histoname = "scanProcessTotals"; 
	histoname += "_";
	histoname +=im0;
	histoname += "_";
	histoname +=im12;

	pair<int,int> iscanpoint = make_pair(im0,im12);

	//new 2D variations
	//for SMS we have a choice -- do we leave the histo with 10 bins in x and only just bin NotFound? or do we just make it 1 bin? TODO
	histoname.ReplaceAll("scanProcessTotals","scanProcessTotalsCTEQ");
	scanProcessTotalsMapCTEQ[iscanpoint] = new TH2D(histoname,histoname,int(NotFound),int(ng),int(NotFound),45,0,45);
	histoname.ReplaceAll("scanProcessTotalsCTEQ","scanProcessTotalsMSTW");
	scanProcessTotalsMapMSTW[iscanpoint] = new TH2D(histoname,histoname,int(NotFound),int(ng),int(NotFound),41,0,41);
	histoname.ReplaceAll("scanProcessTotalsMSTW","scanProcessTotalsNNPDF");
	scanProcessTotalsMapNNPDF[iscanpoint] = new TH2D(histoname,histoname,int(NotFound),int(ng),int(NotFound),100,0,100);
      }
    }
  }
  else if (theScanType_==kNotScan) {
    pair<int,int> iscanpoint = make_pair(0,0);
    TString histoname = "scanProcessTotals"; 
    //same question here -- do we leave the histo with 10 bins in x and only just bin NotFound? or do we just make it 1 bin? TODO
    histoname.ReplaceAll("scanProcessTotals","scanProcessTotalsCTEQ");
    scanProcessTotalsMapCTEQ[iscanpoint] = new TH2D(histoname,histoname,int(NotFound),int(ng),int(NotFound),45,0,45);
    histoname.ReplaceAll("scanProcessTotalsCTEQ","scanProcessTotalsMSTW");
    scanProcessTotalsMapMSTW[iscanpoint] = new TH2D(histoname,histoname,int(NotFound),int(ng),int(NotFound),41,0,41);
    histoname.ReplaceAll("scanProcessTotalsMSTW","scanProcessTotalsNNPDF");
    scanProcessTotalsMapNNPDF[iscanpoint] = new TH2D(histoname,histoname,int(NotFound),int(ng),int(NotFound),100,0,100);
  }
  else assert(0);

  TH2D* scanSMSngen=0; //should be redundant to the scanProcessTotals histograms, but let's keep it for now
  if (theScanType_==kSMS) scanSMSngen = new TH2D("scanSMSngen","number of generated events",80,0,2000,80,0,2000); //mgluino,mLSP

  const  bool puReweightIs1D = true;//((theScanType_!=kNotScan) || sampleName_.Contains("QCD")); //CFA -- no 3D for now

  //initialize PU things
  std::vector< float > DataDist;
  std::vector< float > DataDistSystVar;
  std::vector< float > MCDist;
  for( int i=0; i<60; ++i) {
    //1D reweighting -- all I have done for 2012 (not sure if 3D is used anymore)
    if (puReweightIs1D) {    
      DataDist.push_back(pu::RunsThrough203002[i]);
      DataDistSystVar.push_back(pu::RunsThrough203002systVar[i]);
    }
    //for 3dPU reweighting 
    //    else DataDist2011.push_back(pu::TrueDist2011_f[i]);
    if (sampleName_.Contains("PU_S7"))  MCDist.push_back(pu::Summer2012[i]);
    else if (sampleName_.Contains("PU_S6")) MCDist.push_back(pu::Fall2011[i]);
    else if (sampleName_.Contains("PU_S10")) MCDist.push_back(pu::Summer2012_S10[i]);
    else MCDist.push_back(pu::Summer2012[i]); //just a safety valve
  }

  reweight::LumiReWeighting * LumiWeights=0;
  reweight::LumiReWeighting * LumiWeightsSystVar=0;
  // FIXME CFA
  //Lumi3DReWeighting * LumiWeights3D=0;
  if (puReweightIs1D)  {
    LumiWeights = new reweight::LumiReWeighting( MCDist, DataDist );
    LumiWeightsSystVar = new reweight::LumiReWeighting( MCDist, DataDistSystVar );
  }
  else                 assert(0); //LumiWeights3D = new Lumi3DReWeighting( MCDist2011, DataDist2011);

  //FIXME CFA...this is for PU systematics
//   if (!puReweightIs1D) {
//     if(thePUuncType_ == kPUunc0) {
//       LumiWeights3D->weight3D_init(1);    
//     }
//     //8% uncertainty in total inelastic cross-section (68 mb vs 73.5 mb)
//     //see https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/1479/3/1/1.html
//     else if(thePUuncType_ == kPUuncDown) {
//       LumiWeights3D->weight3D_init(0.92);    
//     }
//     else if(thePUuncType_ == kPUuncUp) {
//       LumiWeights3D->weight3D_init(1.08);    
//     }
//   }

  // ~~~~~~~ define tree branches ~~~~~~~
  reducedTree.Branch("weight",&eventweight,"weight/D"); //special case; 'weight' is already taken here but I want the reducedTree to use it
  reducedTree.Branch("weight2",&eventweight2,"weight2/D");
  reducedTree.Branch("scanCrossSection",&scanCrossSection,"scanCrossSection/D");
  reducedTree.Branch("scanCrossSectionPlus",&scanCrossSectionPlus,"scanCrossSectionPlus/D");
  reducedTree.Branch("scanCrossSectionMinus",&scanCrossSectionMinus,"scanCrossSectionMinus/D");
  reducedTree.Branch("runNumber",&runNumber,"runNumber/l");
  reducedTree.Branch("lumiSection",&lumiSection,"lumiSection/l");
  reducedTree.Branch("eventNumber",&eventNumber,"eventNumber/l");

  reducedTree.Branch("m0",&m0,"m0/I");
  reducedTree.Branch("m12",&m12,"m12/I");


//   reducedTree.Branch("W1decayType",&W1decayType,"W1decayType/I");
//   reducedTree.Branch("W2decayType",&W2decayType,"W2decayType/I");
//   reducedTree.Branch("decayType",&decayType,"decayType/I");
  reducedTree.Branch("ttbarDecayCode",&ttbarDecayCode,"ttbarDecayCode/I");

  //  reducedTree.Branch("btagIPweight",&btagIPweight,"btagIPweight/F");
  //  reducedTree.Branch("pfmhtweight",&pfmhtweight,"pfmhtweight/F");
  reducedTree.Branch("PUweight",&PUweight,"PUweight/F");
  reducedTree.Branch("PUweightSystVar",&PUweightSystVar,"PUweightSystVar/F");
  reducedTree.Branch("hltHTeff",&hltHTeff,"hltHTeff/F");
  reducedTree.Branch("hltMHTeff",&hltMHTeff,"hltMHTeff/F");

//   reducedTree.Branch("hltMHTeffBNN",&hltMHTeffBNN,"hltMHTeffBNN/F");
//   reducedTree.Branch("hltMHTeffBNNUp",&hltMHTeffBNNUp,"hltMHTeffBNNUp/F");
//   reducedTree.Branch("hltMHTeffBNNDown",&hltMHTeffBNNDown,"hltMHTeffBNNDown/F");

  //got to store the whole vector. very big, unfortunately
  reducedTree.Branch("pdfWeightsCTEQ",&pdfWeightsCTEQ,"pdfWeightsCTEQ[45]/F");
  reducedTree.Branch("pdfWeightsMSTW",&pdfWeightsMSTW,"pdfWeightsMSTW[41]/F");
  reducedTree.Branch("pdfWeightsNNPDF",&pdfWeightsNNPDF,"pdfWeightsNNPDF[100]/F");

  reducedTree.Branch("prob0",&prob0,"prob0/F");
  reducedTree.Branch("probge1",&probge1,"probge1/F");
  reducedTree.Branch("prob1",&prob1,"prob1/F");
  reducedTree.Branch("probge2",&probge2,"probge2/F");
  reducedTree.Branch("prob2",&prob2,"prob2/F");
  reducedTree.Branch("probge3",&probge3,"probge3/F");

  //for tag eff variations
  reducedTree.Branch("prob0_HFplus",&prob0_HFplus,"prob0_HFplus/F");
  reducedTree.Branch("probge1_HFplus",&probge1_HFplus,"probge1_HFplus/F");
  reducedTree.Branch("prob1_HFplus",&prob1_HFplus,"prob1_HFplus/F");
  reducedTree.Branch("probge2_HFplus",&probge2_HFplus,"probge2_HFplus/F");
  reducedTree.Branch("prob2_HFplus",&prob2_HFplus,"prob2_HFplus/F");
  reducedTree.Branch("probge3_HFplus",&probge3_HFplus,"probge3_HFplus/F");

  reducedTree.Branch("prob0_HFminus",&prob0_HFminus,"prob0_HFminus/F");
  reducedTree.Branch("probge1_HFminus",&probge1_HFminus,"probge1_HFminus/F");
  reducedTree.Branch("prob1_HFminus",&prob1_HFminus,"prob1_HFminus/F");
  reducedTree.Branch("probge2_HFminus",&probge2_HFminus,"probge2_HFminus/F");
  reducedTree.Branch("prob2_HFminus",&prob2_HFminus,"prob2_HFminus/F");
  reducedTree.Branch("probge3_HFminus",&probge3_HFminus,"probge3_HFminus/F");

  reducedTree.Branch("prob0_LFplus",&prob0_LFplus,"prob0_LFplus/F");
  reducedTree.Branch("probge1_LFplus",&probge1_LFplus,"probge1_LFplus/F");
  reducedTree.Branch("prob1_LFplus",&prob1_LFplus,"prob1_LFplus/F");
  reducedTree.Branch("probge2_LFplus",&probge2_LFplus,"probge2_LFplus/F");
  reducedTree.Branch("prob2_LFplus",&prob2_LFplus,"prob2_LFplus/F");
  reducedTree.Branch("probge3_LFplus",&probge3_LFplus,"probge3_LFplus/F");

  reducedTree.Branch("prob0_LFminus",&prob0_LFminus,"prob0_LFminus/F");
  reducedTree.Branch("probge1_LFminus",&probge1_LFminus,"probge1_LFminus/F");
  reducedTree.Branch("prob1_LFminus",&prob1_LFminus,"prob1_LFminus/F");
  reducedTree.Branch("probge2_LFminus",&probge2_LFminus,"probge2_LFminus/F");
  reducedTree.Branch("prob2_LFminus",&prob2_LFminus,"prob2_LFminus/F");
  reducedTree.Branch("probge3_LFminus",&probge3_LFminus,"probge3_LFminus/F");

  reducedTree.Branch("cutPV",&cutPV,"cutPV/O");
  reducedTree.Branch("cutTrigger",&cutTrigger,"cutTrigger/O");
  reducedTree.Branch("cutTrigger2",&cutTrigger2,"cutTrigger2/O");
  
  reducedTree.Branch("csctighthaloFilter",&csctighthaloFilter,"csctighthaloFilter/O");
  reducedTree.Branch("eenoiseFilter",&eenoiseFilter,"eenoiseFilter/O");
  reducedTree.Branch("greedymuonFilter",&greedymuonFilter,"greedymuonFilter/O");
  reducedTree.Branch("hbhenoiseFilter",&hbhenoiseFilter,"hbhenoiseFilter/O");
  reducedTree.Branch("inconsistentmuonFilter",&inconsistentmuonFilter,"inconsistentmuonFilter/O");
  reducedTree.Branch("ra2ecaltpFilter",&ra2ecaltpFilter,"ra2ecaltpFilter/O");
  reducedTree.Branch("scrapingvetoFilter",&scrapingvetoFilter,"scrapingvetoFilter/O");
  reducedTree.Branch("trackingfailureFilter",&trackingfailureFilter,"trackingfailureFilter/O");
  reducedTree.Branch("badjetFilter",&badjetFilter,"badjetFilter/O");
  reducedTree.Branch("passCleaning",&passCleaning,"passCleaning/O");
  reducedTree.Branch("PBNRcode",&PBNRcode,"PBNRcode/I");

  reducedTree.Branch("buggyEvent",&buggyEvent,"buggyEvent/O");

  reducedTree.Branch("nGoodPV",&nGoodPV,"nGoodPV/I");

  reducedTree.Branch("SUSY_nb",&SUSY_nb,"SUSY_nb/I");
  reducedTree.Branch("SUSY_process",&SUSY_process,"SUSY_process/I");
  reducedTree.Branch("SUSY_recoilPt",&SUSY_recoilPt,"SUSY_recoilPt/F");
  reducedTree.Branch("nCorrectRecoStop",&nCorrectRecoStop,"nCorrectRecoStop/I");

  reducedTree.Branch("njets",&njets,"njets/I");
  reducedTree.Branch("njets30",&njets30,"njets30/I");
  reducedTree.Branch("nbjets",&nbjets,"nbjets/I");
  reducedTree.Branch("nbjets30",&nbjets30,"nbjets30/I");
  reducedTree.Branch("ntruebjets",&ntruebjets,"ntruebjets/I");
  reducedTree.Branch("nElectrons",&nElectrons,"nElectrons/I");
  reducedTree.Branch("nMuons",&nMuons,"nMuons/I");
  reducedTree.Branch("nTaus",&nTaus,"nTaus/I");
  reducedTree.Branch("nElectrons5",&nElectrons5,"nElectrons5/I");
  reducedTree.Branch("nMuons5",&nMuons5,"nMuons5/I");
  reducedTree.Branch("nElectrons15",&nElectrons15,"nElectrons15/I");
  reducedTree.Branch("nMuons15",&nMuons15,"nMuons15/I");

  reducedTree.Branch("nElectrons20",&nElectrons20,"nElectrons20/I"); //jmt
  reducedTree.Branch("nMuons20",&nMuons20,"nMuons20/I");

  reducedTree.Branch("bestZmass",&bestZmass,"bestZmass/F");
  reducedTree.Branch("mjj1",&mjj1,"mjj1/F");
  reducedTree.Branch("mjj2",&mjj2,"mjj2/F");
  reducedTree.Branch("mjjdiff",&mjjdiff,"mjjdiff/F");
  reducedTree.Branch("mjj1_5",&mjj1_5,"mjj1_5/F");
  reducedTree.Branch("mjj2_5",&mjj2_5,"mjj2_5/F");
  reducedTree.Branch("mjjdiff_5",&mjjdiff_5,"mjjdiff_5/F");

  reducedTree.Branch("mjjb1",&mjjb1,"mjjb1/F");
  reducedTree.Branch("mjjb2",&mjjb2,"mjjb2/F");
  reducedTree.Branch("topPT1",&topPT1,"topPT1/F");
  reducedTree.Branch("topPT2",&topPT2,"topPT2/F");

  reducedTree.Branch("nbjetsSSVM",&nbjetsSSVM,"nbjetsSSVM/I");
  reducedTree.Branch("nbjetsSSVHPT",&nbjetsSSVHPT,"nbjetsSSVHPT/I");
  reducedTree.Branch("nbjetsTCHPT",&nbjetsTCHPT,"nbjetsTCHPT/I");
  reducedTree.Branch("nbjetsTCHET",&nbjetsTCHET,"nbjetsTCHET/I");
  reducedTree.Branch("nbjetsTCHPM",&nbjetsTCHPM,"nbjetsTCHPM/I");
  reducedTree.Branch("nbjetsCSVM",&nbjetsCSVM,"nbjetsCSVM/I");
  reducedTree.Branch("nbjetsCSVL",&nbjetsCSVL,"nbjetsCSVL/I");

  reducedTree.Branch("isRealData",&isRealData,"isRealData/O");
  reducedTree.Branch("pass_utilityHLT",&pass_utilityHLT,"pass_utilityHLT/O");
  reducedTree.Branch("prescaleUtilityHLT", &prescaleUtilityHLT, "prescaleUtilityHLT/i");
  reducedTree.Branch("versionUtilityHLT", &versionUtilityHLT, "versionUtilityHLT/i");
  reducedTree.Branch("pass_utilityPrescaleModuleHLT",&pass_utilityPrescaleModuleHLT,"pass_utilityPrescaleModuleHLT/O");

  for ( map<TString,triggerData>::iterator itrigger = triggerlist.begin() ; itrigger!=triggerlist.end() ; ++itrigger) {
    TString treename = "pass_";
    treename+=itrigger->first;
    TString withcode = treename+"/O";
    reducedTree.Branch(treename,&(itrigger->second.pass),withcode);
  }

  //now for the mc trigger results
  for ( map<TString,triggerData>::iterator itrigger = triggerlist_mc.begin() ; itrigger!=triggerlist_mc.end() ; ++itrigger) {
    TString treename = "passMC_";
    treename+=itrigger->first;
    TString withcode = treename+"/O";
    reducedTree.Branch(treename,&(itrigger->second.pass),withcode);
  }

  reducedTree.Branch("HT",&HT,"HT/F");
  reducedTree.Branch("HT30",&HT30,"HT30/F");
  reducedTree.Branch("ST",&ST,"ST/F"); //includes HT + leptons
  reducedTree.Branch("STeff",&STeff,"STeff/F"); //includes HT + leptons + MET
  reducedTree.Branch("MET",&MET,"MET/F");
  reducedTree.Branch("METphi",&METphi,"METphi/F");
  reducedTree.Branch("MHT",&MHT,"MHT/F");

  reducedTree.Branch("METsig",&METsig,"METsig/F");
  reducedTree.Branch("METsig00",&METsig00,"METsig00/F");
  reducedTree.Branch("METsig10",&METsig10,"METsig10/F");
  reducedTree.Branch("METsig11",&METsig11,"METsig11/F");


  reducedTree.Branch("caloMET",&caloMET,"caloMET/F");
  reducedTree.Branch("caloMETphi",&caloMETphi,"caloMETphi/F");
  reducedTree.Branch("rawPFMET",&rawPFMET, "rawPFMET/F");
  reducedTree.Branch("rawPFMETphi",&rawPFMETphi, "rawPFMETphi/F");
  
  reducedTree.Branch("bestWMass",&wMass,"bestWMass/F");
  reducedTree.Branch("bestTopMass",&topMass,"bestTopMass/F");
  reducedTree.Branch("topCosHel",&topCosHel,"topCosHel/F");
  reducedTree.Branch("WCosHel",&wCosHel,"WCosHel/F");
  reducedTree.Branch("MT_b",&MT_b, "MT_b/F");
  reducedTree.Branch("MT_Wlep",&MT_Wlep, "MT_Wlep/F");
  reducedTree.Branch("MT_Wlep5",&MT_Wlep5, "MT_Wlep5/F");
  reducedTree.Branch("MT_Wlep15",&MT_Wlep15, "MT_Wlep15/F");

  reducedTree.Branch("deltaR_bestTwoCSV",&deltaR_bestTwoCSV,"deltaR_bestTwoCSV/F");
  reducedTree.Branch("mjj_bestTwoCSV",&mjj_bestTwoCSV,"mjj_bestTwoCSV/F");
  reducedTree.Branch("deltaR_closestB",&deltaR_closestB,"deltaR_closestB/F");
  reducedTree.Branch("mjj_closestB",&mjj_closestB,"mjj_closestB/F");

  reducedTree.Branch("minDeltaPhi",&minDeltaPhi,"minDeltaPhi/F");
  reducedTree.Branch("minDeltaPhiAll",&minDeltaPhiAll,"minDeltaPhiAll/F");
  reducedTree.Branch("minDeltaPhiAll30",&minDeltaPhiAll30,"minDeltaPhiAll30/F");
  reducedTree.Branch("minDeltaPhi30_eta5_noIdAll",&minDeltaPhi30_eta5_noIdAll,"minDeltaPhi30_eta5_noIdAll/F");

  reducedTree.Branch("minDeltaPhiMetTau",&minDeltaPhiMetTau,"minDeltaPhiMetTau/F");

  reducedTree.Branch("deltaPhi1", &deltaPhi1, "deltaPhi1/F");
  reducedTree.Branch("deltaPhi2", &deltaPhi2, "deltaPhi2/F");
  reducedTree.Branch("deltaPhi3", &deltaPhi3, "deltaPhi3/F");

  reducedTree.Branch("maxDeltaPhi",&maxDeltaPhi,"maxDeltaPhi/F");
  reducedTree.Branch("maxDeltaPhiAll",&maxDeltaPhiAll,"maxDeltaPhiAll/F");
  reducedTree.Branch("maxDeltaPhiAll30",&maxDeltaPhiAll30,"maxDeltaPhiAll30/F");
  reducedTree.Branch("maxDeltaPhi30_eta5_noIdAll",&maxDeltaPhi30_eta5_noIdAll,"maxDeltaPhi30_eta5_noIdAll/F");

  reducedTree.Branch("sumDeltaPhi",&sumDeltaPhi,"sumDeltaPhi/F");
  reducedTree.Branch("diffDeltaPhi",&diffDeltaPhi,"diffDeltaPhi/F");

  reducedTree.Branch("deltaPhiStar",&deltaPhiStar,"deltaPhiStar/F");
  reducedTree.Branch("deltaPhiStar_badjet_pt",&deltaPhiStar_badjet_pt,"deltaPhiStar_badjet_pt/F");
  reducedTree.Branch("deltaPhiStar_badjet_phi",&deltaPhiStar_badjet_phi,"deltaPhiStar_badjet_phi/F");
  reducedTree.Branch("deltaPhiStar_badjet_eta",&deltaPhiStar_badjet_eta,"deltaPhiStar_badjet_eta/F");

  reducedTree.Branch("minDeltaPhiN", &minDeltaPhiN, "minDeltaPhiN/F");
  reducedTree.Branch("deltaPhiN1", &deltaPhiN1, "deltaPhiN1/F");
  reducedTree.Branch("deltaPhiN2", &deltaPhiN2, "deltaPhiN2/F");
  reducedTree.Branch("deltaPhiN3", &deltaPhiN3, "deltaPhiN3/F");
  reducedTree.Branch("minDeltaPhi_chosenJet", &minDeltaPhi_chosenJet, "minDeltaPhi_chosenJet/i");
  reducedTree.Branch("minDeltaPhiN_chosenJet", &minDeltaPhiN_chosenJet, "minDeltaPhiN_chosenJet/i");
  reducedTree.Branch("maxJetMis_chosenJet", &maxJetMis_chosenJet, "maxJetMis_chosenJet/i");
  reducedTree.Branch("maxJetFracMis_chosenJet", &maxJetFracMis_chosenJet, "maxJetFracMis_chosenJet/i");

  reducedTree.Branch("minDeltaPhiN_deltaT", &minDeltaPhiN_deltaT, "minDeltaPhiN_deltaT/F");
  reducedTree.Branch("deltaT1", &deltaT1, "deltaT1/F");
  reducedTree.Branch("deltaT2", &deltaT2, "deltaT2/F");
  reducedTree.Branch("deltaT3", &deltaT3, "deltaT3/F");

  reducedTree.Branch("minDeltaPhiN_asin", &minDeltaPhiN_asin, "minDeltaPhiN_asin/F");
  reducedTree.Branch("deltaPhiN_asin1", &deltaPhiN_asin1, "deltaPhiN_asin1/F");
  reducedTree.Branch("deltaPhiN_asin2", &deltaPhiN_asin2, "deltaPhiN_asin2/F");
  reducedTree.Branch("deltaPhiN_asin3", &deltaPhiN_asin3, "deltaPhiN_asin3/F");

  //variables for qcd studies. don't include if this isn't QCD MC (just to save time and space)
  if (isSampleQCD() ) {
  reducedTree.Branch("minDeltaPhiN_otherEta5", &minDeltaPhiN_otherEta5, "minDeltaPhiN_otherEta5/F");
  reducedTree.Branch("minDeltaPhiN_otherEta5idNo", &minDeltaPhiN_otherEta5idNo, "minDeltaPhiN_otherEta5idNo/F");
  reducedTree.Branch("minDeltaPhiN_mainPt30_otherEta5idNo", &minDeltaPhiN_mainPt30_otherEta5idNo, "minDeltaPhiN_mainPt30_otherEta5idNo/F");
  reducedTree.Branch("minDeltaPhiN_mainPt30Eta5_otherEta5idNo", &minDeltaPhiN_mainPt30Eta5_otherEta5idNo, "minDeltaPhiN_mainPt30Eta5_otherEta5idNo/F");
  reducedTree.Branch("minDeltaPhiN_mainPt30Eta5idNo_otherEta5idNo", &minDeltaPhiN_mainPt30Eta5idNo_otherEta5idNo, "minDeltaPhiN_mainPt30Eta5idNo_otherEta5idNo/F");
  reducedTree.Branch("minDeltaPhiK", &minDeltaPhiK, "minDeltaPhiK/F");
  reducedTree.Branch("minDeltaPhiK_otherEta5", &minDeltaPhiK_otherEta5, "minDeltaPhiK_otherEta5/F");
  reducedTree.Branch("minDeltaPhiK_otherEta5idNo", &minDeltaPhiK_otherEta5idNo, "minDeltaPhiK_otherEta5idNo/F");
  reducedTree.Branch("minDeltaPhiK_mainPt30_otherEta5idNo", &minDeltaPhiK_mainPt30_otherEta5idNo, "minDeltaPhiK_mainPt30_otherEta5idNo/F");
  reducedTree.Branch("minDeltaPhiK_mainPt30Eta5_otherEta5idNo", &minDeltaPhiK_mainPt30Eta5_otherEta5idNo, "minDeltaPhiK_mainPt30Eta5_otherEta5idNo/F");
  reducedTree.Branch("minDeltaPhiK_mainPt30Eta5idNo_otherEta5idNo", &minDeltaPhiK_mainPt30Eta5idNo_otherEta5idNo, "minDeltaPhiK_mainPt30Eta5idNo_otherEta5idNo/F");
  reducedTree.Branch("minDeltaPhiN_DJR", &minDeltaPhiN_DJR, "minDeltaPhiN_DJR/F");
  reducedTree.Branch("minDeltaPhiN_DJR_otherEta5", &minDeltaPhiN_DJR_otherEta5, "minDeltaPhiN_DJR_otherEta5/F");
  reducedTree.Branch("minDeltaPhiK_DJR", &minDeltaPhiK_DJR, "minDeltaPhiK_DJR/F");
  reducedTree.Branch("minDeltaPhiK_DJR_otherEta5", &minDeltaPhiK_DJR_otherEta5, "minDeltaPhiK_DJR_otherEta5/F");
 
  reducedTree.Branch("minDeltaPhiN_otherPt10", &minDeltaPhiN_otherPt10, "minDeltaPhiN_otherPt10/F");
  reducedTree.Branch("minDeltaPhiN_otherPt20", &minDeltaPhiN_otherPt20, "minDeltaPhiN_otherPt20/F");
  reducedTree.Branch("minDeltaPhiN_otherPt30", &minDeltaPhiN_otherPt30, "minDeltaPhiN_otherPt30/F");
  reducedTree.Branch("minDeltaPhiN_otherPt40", &minDeltaPhiN_otherPt40, "minDeltaPhiN_otherPt40/F");
  reducedTree.Branch("minDeltaPhiN_otherPt50", &minDeltaPhiN_otherPt50, "minDeltaPhiN_otherPt50/F");
  reducedTree.Branch("minDeltaPhiN_DJR_otherPt10", &minDeltaPhiN_DJR_otherPt10, "minDeltaPhiN_DJR_otherPt10/F");
  reducedTree.Branch("minDeltaPhiN_DJR_otherPt20", &minDeltaPhiN_DJR_otherPt20, "minDeltaPhiN_DJR_otherPt20/F");
  reducedTree.Branch("minDeltaPhiN_DJR_otherPt30", &minDeltaPhiN_DJR_otherPt30, "minDeltaPhiN_DJR_otherPt30/F");
  reducedTree.Branch("minDeltaPhiN_DJR_otherPt40", &minDeltaPhiN_DJR_otherPt40, "minDeltaPhiN_DJR_otherPt40/F");
  reducedTree.Branch("minDeltaPhiN_DJR_otherPt50", &minDeltaPhiN_DJR_otherPt50, "minDeltaPhiN_DJR_otherPt50/F");
  
  reducedTree.Branch("minDeltaPhiN_withLepton", &minDeltaPhiN_withLepton, "minDeltaPhiN_withLepton/F");

  reducedTree.Branch("maxJetMis",&maxJetMis,"maxJetMis/F");
  reducedTree.Branch("max2JetMis",&max2JetMis,"max2JetMis/F");
  reducedTree.Branch("maxJetMisAll30",&maxJetMisAll30,"maxJetMisAll30/F");
  reducedTree.Branch("max2JetMisAll30",&max2JetMisAll30,"max2JetMisAll30/F");
  reducedTree.Branch("maxJetFracMis",&maxJetFracMis,"maxJetFracMis/F");
  reducedTree.Branch("max2JetFracMis",&max2JetFracMis,"max2JetFracMis/F");
  reducedTree.Branch("maxJetFracMisAll30",&maxJetFracMisAll30,"maxJetFracMisAll30/F");
  reducedTree.Branch("max2JetFracMisAll30",&max2JetFracMisAll30,"max2JetFracMisAll30/F");
  reducedTree.Branch("deltaPhiMETJetMaxMis",&deltaPhiMETJetMaxMis,"deltaPhiMETJetMaxMis/F");
  reducedTree.Branch("deltaPhiMETJetMaxMis30",&deltaPhiMETJetMaxMis30,"deltaPhiMETJetMaxMis30/F");
  } //if sample is qcd
  
  reducedTree.Branch("CSVout1",&CSVout1,"CSVout1/F");
  reducedTree.Branch("CSVout2",&CSVout2,"CSVout2/F");
  reducedTree.Branch("CSVout3",&CSVout3,"CSVout3/F");
  reducedTree.Branch("minDeltaPhiAllb30",&minDeltaPhiAllb30,"minDeltaPhiAllb30/F");
  reducedTree.Branch("deltaPhib1",&deltaPhib1,"deltaPhib1/F");
  reducedTree.Branch("deltaPhib2",&deltaPhib2,"deltaPhib2/F");
  reducedTree.Branch("deltaPhib3",&deltaPhib3,"deltaPhib3/F");
  reducedTree.Branch("minDeltaPhiMETMuonsAll",&minDeltaPhiMETMuonsAll,"minDeltaPhiMETMuonsAll/F");

  reducedTree.Branch("minDeltaPhiN_lostJet", &minDeltaPhiN_lostJet, "minDeltaPhiN_lostJet/F");
  reducedTree.Branch("deltaPhiN1_lostJet", &deltaPhiN1_lostJet, "deltaPhiN1_lostJet/F");
  reducedTree.Branch("deltaPhiN2_lostJet", &deltaPhiN2_lostJet, "deltaPhiN2_lostJet/F");
  reducedTree.Branch("deltaPhiN3_lostJet", &deltaPhiN3_lostJet, "deltaPhiN3_lostJet/F");

  if (isSampleQCD() ) {
  reducedTree.Branch("passBadECAL_METphi53_n10_s12", &passBadECAL_METphi53_n10_s12, "passBadECAL_METphi53_n10_s12/O");
  reducedTree.Branch("worstMisA_badECAL_METphi53_n10_s12",&worstMisA_badECAL_METphi53_n10_s12, "worstMisA_badECAL_METphi53_n10_s12/F");
  reducedTree.Branch("worstMisF_badECAL_METphi53_n10_s12",&worstMisF_badECAL_METphi53_n10_s12, "worstMisF_badECAL_METphi53_n10_s12/F");

  reducedTree.Branch("passBadECAL_METphi33_n10_s12", &passBadECAL_METphi33_n10_s12, "passBadECAL_METphi33_n10_s12/O");
  reducedTree.Branch("worstMisA_badECAL_METphi33_n10_s12",&worstMisA_badECAL_METphi33_n10_s12, "worstMisA_badECAL_METphi33_n10_s12/F");
  reducedTree.Branch("worstMisF_badECAL_METphi33_n10_s12",&worstMisF_badECAL_METphi33_n10_s12, "worstMisF_badECAL_METphi33_n10_s12/F");


  reducedTree.Branch("passBadECAL_METphi52_n10_s12", &passBadECAL_METphi52_n10_s12, "passBadECAL_METphi52_n10_s12/O");
  reducedTree.Branch("worstMisA_badECAL_METphi52_n10_s12",&worstMisA_badECAL_METphi52_n10_s12, "worstMisA_badECAL_METphi52_n10_s12/F");
  reducedTree.Branch("worstMisF_badECAL_METphi52_n10_s12",&worstMisF_badECAL_METphi52_n10_s12, "worstMisF_badECAL_METphi52_n10_s12/F");

  reducedTree.Branch("passBadECAL_METphi53_n5_s12", &passBadECAL_METphi53_n5_s12, "passBadECAL_METphi53_n5_s12/O");
  reducedTree.Branch("passBadECAL_METphi53_crack", &passBadECAL_METphi53_crack, "passBadECAL_METphi53_crack/O");
  reducedTree.Branch("worstMisA_badECAL_METphi53_crack",&worstMisA_badECAL_METphi53_crack, "worstMisA_badECAL_METphi53_crack/F");
  reducedTree.Branch("worstMisF_badECAL_METphi53_crack",&worstMisF_badECAL_METphi53_crack, "worstMisF_badECAL_METphi53_crack/F");
  reducedTree.Branch("passBadECAL_METphi52_crack", &passBadECAL_METphi52_crack, "passBadECAL_METphi52_crack/O");
  reducedTree.Branch("worstMisA_badECAL_METphi52_crack",&worstMisA_badECAL_METphi52_crack, "worstMisA_badECAL_METphi52_crack/F");
  reducedTree.Branch("worstMisF_badECAL_METphi52_crack",&worstMisF_badECAL_METphi52_crack, "worstMisF_badECAL_METphi52_crack/F");

  reducedTree.Branch("passBadECAL_METphi53_crackF", &passBadECAL_METphi53_crackF, "passBadECAL_METphi53_crackF/O");
  reducedTree.Branch("worstMisA_badECAL_METphi53_crackF",&worstMisA_badECAL_METphi53_crackF, "worstMisA_badECAL_METphi53_crackF/F");
  reducedTree.Branch("worstMisF_badECAL_METphi53_crackF",&worstMisF_badECAL_METphi53_crackF, "worstMisF_badECAL_METphi53_crackF/F");
  reducedTree.Branch("passBadECAL_METphi52_crackF", &passBadECAL_METphi52_crackF, "passBadECAL_METphi52_crackF/O");
  reducedTree.Branch("worstMisA_badECAL_METphi52_crackF",&worstMisA_badECAL_METphi52_crackF, "worstMisA_badECAL_METphi52_crackF/F");
  reducedTree.Branch("worstMisF_badECAL_METphi52_crackF",&worstMisF_badECAL_METphi52_crackF, "worstMisF_badECAL_METphi52_crackF/F");
  }

  reducedTree.Branch("jetpt1",&jetpt1,"jetpt1/F");
  reducedTree.Branch("jetgenpt1",&jetgenpt1,"jetgenpt1/F");
  reducedTree.Branch("jeteta1",&jeteta1,"jeteta1/F");
  reducedTree.Branch("jetgeneta1",&jetgeneta1,"jetgeneta1/F");
  reducedTree.Branch("jetphi1",&jetphi1,"jetphi1/F");
  reducedTree.Branch("jetgenphi1",&jetgenphi1,"jetgenphi1/F");
  reducedTree.Branch("jetenergy1",&jetenergy1,"jetenergy1/F");
  reducedTree.Branch("jetflavor1",&jetflavor1,"jetflavor1/I");
  reducedTree.Branch("jetchargedhadronfrac1",&jetchargedhadronfrac1,"jetchargedhadronfrac1/F");
  reducedTree.Branch("jetchargedhadronmult1",&jetchargedhadronmult1,"jetchargedhadronmult1/I");
  reducedTree.Branch("jetneutralhadronmult1",&jetneutralhadronmult1,"jetneutralhadronmult1/I");
  reducedTree.Branch("jetmumult1",&jetmumult1,"jetmumult1/I");

  reducedTree.Branch("alletajetpt1",&alletajetpt1,"alletajetpt1/F");
  reducedTree.Branch("alletajetphi1",&alletajetphi1,"alletajetphi1/F");
  reducedTree.Branch("alletajeteta1",&alletajeteta1,"alletajeteta1/F");
  reducedTree.Branch("alletajetneutralhadronfrac1",&alletajetneutralhadronfrac1,"alletajetneutralhadronfrac1/F");
  reducedTree.Branch("alletajetneutralemfrac1",&alletajetneutralemfrac1,"alletajetneutralemfrac1/F");
  reducedTree.Branch("alletajetphotonfrac1",&alletajetphotonfrac1,"alletajetphotonfrac1/F");

  reducedTree.Branch("jetpt2",&jetpt2,"jetpt2/F");
  reducedTree.Branch("jetgenpt2",&jetgenpt2,"jetgenpt2/F");
  reducedTree.Branch("jeteta2",&jeteta2,"jeteta2/F");
  reducedTree.Branch("jetgeneta2",&jetgeneta2,"jetgeneta2/F");
  reducedTree.Branch("jetphi2",&jetphi2,"jetphi2/F");
  reducedTree.Branch("jetgenphi2",&jetgenphi2,"jetgenphi2/F");
  reducedTree.Branch("jetenergy2",&jetenergy2,"jetenergy2/F");
  reducedTree.Branch("jetflavor2",&jetflavor2,"jetflavor2/I");
  reducedTree.Branch("jetchargedhadronfrac2",&jetchargedhadronfrac2,"jetchargedhadronfrac2/F");
  reducedTree.Branch("jetchargedhadronmult2",&jetchargedhadronmult2,"jetchargedhadronmult2/I");
  reducedTree.Branch("jetneutralhadronmult2",&jetneutralhadronmult2,"jetneutralhadronmult2/I");
  reducedTree.Branch("jetmumult2",&jetmumult2,"jetmumult2/I");

  reducedTree.Branch("jetpt3",&jetpt3,"jetpt3/F");
  reducedTree.Branch("jetgenpt3",&jetgenpt3,"jetgenpt3/F");
  reducedTree.Branch("jeteta3",&jeteta3,"jeteta3/F");
  reducedTree.Branch("jetgeneta3",&jetgeneta3,"jetgeneta3/F");
  reducedTree.Branch("jetphi3",&jetphi3,"jetphi3/F");
  reducedTree.Branch("jetgenphi3",&jetgenphi3,"jetgenphi3/F");
  reducedTree.Branch("jetenergy3",&jetenergy3,"jetenergy3/F");
  reducedTree.Branch("jetflavor3",&jetflavor3,"jetflavor3/I");
  reducedTree.Branch("jetchargedhadronfrac3",&jetchargedhadronfrac3,"jetchargedhadronfrac3/F");
  reducedTree.Branch("jetchargedhadronmult3",&jetchargedhadronmult3,"jetchargedhadronmult3/I");
  reducedTree.Branch("jetneutralhadronmult3",&jetneutralhadronmult3,"jetneutralhadronmult3/I");
  reducedTree.Branch("jetmumult3",&jetmumult3,"jetmumult3/I");

  reducedTree.Branch("bjetpt1",&bjetpt1,"bjetpt1/F");
  reducedTree.Branch("bjeteta1",&bjeteta1,"bjeteta1/F");
  reducedTree.Branch("bjetphi1",&bjetphi1,"bjetphi1/F");
  reducedTree.Branch("bjetenergy1",&bjetenergy1,"bjetenergy1/F");  
  reducedTree.Branch("bjetflavor1",&bjetflavor1,"bjetflavor1/I");
  reducedTree.Branch("bjetchargedhadronfrac1",&bjetchargedhadronfrac1,"bjetchargedhadronfrac1/F");
  reducedTree.Branch("bjetchargedhadronmult1",&bjetchargedhadronmult1,"bjetchargedhadronmult1/I");

  reducedTree.Branch("bjetpt2",&bjetpt2,"bjetpt2/F");
  reducedTree.Branch("bjeteta2",&bjeteta2,"bjeteta2/F");
  reducedTree.Branch("bjetphi2",&bjetphi2,"bjetphi2/F");
  reducedTree.Branch("bjetenergy2",&bjetenergy2,"bjetenergy2/F");
  reducedTree.Branch("bjetflavor2",&bjetflavor2,"bjetflavor2/I");
  reducedTree.Branch("bjetchargedhadronfrac2",&bjetchargedhadronfrac2,"bjetchargedhadronfrac2/F");
  reducedTree.Branch("bjetchargedhadronmult2",&bjetchargedhadronmult2,"bjetchargedhadronmult2/I");
  
  reducedTree.Branch("bjetpt3",&bjetpt3,"bjetpt3/F");
  reducedTree.Branch("bjeteta3",&bjeteta3,"bjeteta3/F");
  reducedTree.Branch("bjetphi3",&bjetphi3,"bjetphi3/F");
  reducedTree.Branch("bjetenergy3",&bjetenergy3,"bjetenergy3/F");  
  reducedTree.Branch("bjetflavor3",&bjetflavor3,"bjetflavor3/I");
  reducedTree.Branch("bjetchargedhadronfrac3",&bjetchargedhadronfrac3,"bjetchargedhadronfrac3/F");
  reducedTree.Branch("bjetchargedhadronmult3",&bjetchargedhadronmult3,"bjetchargedhadronmult3/I");

  float tauMatch_jetpt,tauMatch_chhadmult,tauMatch_jetcsv,tauMatch_MT;
  reducedTree.Branch("tauMatch_jetpt",&tauMatch_jetpt,"tauMatch_jetpt/F");
  reducedTree.Branch("tauMatch_chhadmult",&tauMatch_chhadmult,"tauMatch_chhadmult/F");
  reducedTree.Branch("tauMatch_jetcsv",&tauMatch_jetcsv,"tauMatch_jetcsv/F");
  reducedTree.Branch("tauMatch_MT",&tauMatch_MT,"tauMatch_MT/F");


  reducedTree.Branch("eleet1",&eleet1,"eleet1/F");
  reducedTree.Branch("elephi1",&elephi1,"elephi1/F");
  reducedTree.Branch("eleeta1",&eleeta1,"eleeta1/F");
  reducedTree.Branch("elecharge1",&elecharge1,"elecharge1/I");
  reducedTree.Branch("muonpt1",&muonpt1,"muonpt1/F");
  reducedTree.Branch("muonphi1",&muonphi1,"muonphi1/F");
  reducedTree.Branch("muoneta1",&muoneta1,"muoneta1/F");
  reducedTree.Branch("muoncharge1",&muoncharge1,"muoncharge1/I");
  reducedTree.Branch("muoniso1",&muoniso1,"muoniso1/F");
  reducedTree.Branch("muonchhadiso1",&muonchhadiso1,"muonchhadiso1/F");
  reducedTree.Branch("muonphotoniso1",&muonphotoniso1,"muonphotoniso1/F");
  reducedTree.Branch("muonneutralhadiso1",&muonneutralhadiso1,"muonneutralhadiso1/F");

  reducedTree.Branch("eleet2",&eleet2,"eleet2/F");
  reducedTree.Branch("elephi2",&elephi2,"elephi2/F");
  reducedTree.Branch("eleeta2",&eleeta2,"eleeta2/F");
  reducedTree.Branch("muonpt2",&muonpt2,"muonpt2/F");
  reducedTree.Branch("muonphi2",&muonphi2,"muonphi2/F");
  reducedTree.Branch("muoneta2",&muoneta2,"muoneta2/F");

  reducedTree.Branch("taupt1",&taupt1,"taupt1/F");
  reducedTree.Branch("taueta1",&taueta1,"taueta1/F");

  reducedTree.Branch("eleRelIso",&eleRelIso,"eleRelIso/F");
  reducedTree.Branch("muonRelIso",&muonRelIso,"muonRelIso/F");


  reducedTree.Branch("recomuonpt1",&recomuonpt1,"recomuonpt1/F");
  reducedTree.Branch("recomuonphi1",&recomuonphi1,"recomuonphi1/F");
  reducedTree.Branch("recomuoneta1",&recomuoneta1,"recomuoneta1/F");
  reducedTree.Branch("recomuoniso1",&recomuoniso1,"recomuoniso1/F");
  reducedTree.Branch("recomuonmindphijet1",&recomuonmindphijet1,"recomuonmindphijet1/F");

  reducedTree.Branch("recomuonpt2",&recomuonpt2,"recomuonpt2/F");
  reducedTree.Branch("recomuonphi2",&recomuonphi2,"recomuonphi2/F");
  reducedTree.Branch("recomuoneta2",&recomuoneta2,"recomuoneta2/F");
  reducedTree.Branch("recomuoniso2",&recomuoniso2,"recomuoniso2/F");
  reducedTree.Branch("recomuonmindphijet2",&recomuonmindphijet1,"recomuonmindphijet2/F");

  reducedTree.Branch("rl",&rl,"rl/F");
  reducedTree.Branch("rMET",&rMET,"rMET/F");

  reducedTree.Branch("transverseSphericity_jets",&transverseSphericity_jets,"transverseSphericity_jets/F");
  reducedTree.Branch("transverseSphericity_jetsMet",&transverseSphericity_jetsMet,"transverseSphericity_jetsMet/F");
  reducedTree.Branch("transverseSphericity_jetsMetLeptons",&transverseSphericity_jetsMetLeptons,"transverseSphericity_jetsMetLeptons/F");
  reducedTree.Branch("transverseSphericity_jetsLeptons",&transverseSphericity_jetsLeptons,"transverseSphericity_jetsLeptons/F");

  reducedTree.Branch("transverseSphericity_jets30",&transverseSphericity_jets30,"transverseSphericity_jets30/F");
  reducedTree.Branch("transverseSphericity_jets30Met",&transverseSphericity_jets30Met,"transverseSphericity_jets30Met/F");
  reducedTree.Branch("transverseSphericity_jets30MetLeptons",&transverseSphericity_jets30MetLeptons,"transverseSphericity_jets30MetLeptons/F");
  reducedTree.Branch("transverseSphericity_jets30Leptons",&transverseSphericity_jets30Leptons,"transverseSphericity_jets30Leptons/F");




/*
  reducedTree.Branch("transverseThrust",&transverseThrust,"transverseThrust/F");
  reducedTree.Branch("transverseThrustPhi",&transverseThrustPhi,"transverseThrustPhi/F");

  reducedTree.Branch("transverseThrustWithMET",&transverseThrustWithMET,"transverseThrustWithMET/F");
  reducedTree.Branch("transverseThrustWithMETPhi",&transverseThrustWithMETPhi,"transverseThrustWithMETPhi/F");
*/
  reducedTree.Branch("minDeltaPhiN_Luke", &minDeltaPhiN_Luke, "minDeltaPhiN_Luke/F");
  reducedTree.Branch("maxDeltaPhiN_Luke", &maxDeltaPhiN_Luke, "maxDeltaPhiN_Luke/F");
  reducedTree.Branch("deltaPhiN1_Luke", &deltaPhiN1_Luke, "deltaPhiN1_Luke/F");
  reducedTree.Branch("deltaPhiN2_Luke", &deltaPhiN2_Luke, "deltaPhiN2_Luke/F");
  reducedTree.Branch("deltaPhiN3_Luke", &deltaPhiN3_Luke, "deltaPhiN3_Luke/F");

  reducedTree.Branch("minTransverseMETSignificance", &minTransverseMETSignificance, "minTransverseMETSignificance/F");
  reducedTree.Branch("maxTransverseMETSignificance", &maxTransverseMETSignificance, "maxTransverseMETSignificance/F");
  reducedTree.Branch("transverseMETSignificance1", &transverseMETSignificance1, "transverseMETSignificance1/F");
  reducedTree.Branch("transverseMETSignificance2", &transverseMETSignificance2, "transverseMETSignificance2/F");
  reducedTree.Branch("transverseMETSignificance3", &transverseMETSignificance3, "transverseMETSignificance3/F");

  reducedTree.Branch("njets_lostJet",&njets_lostJet,"njets_lostJet/I");
  reducedTree.Branch("nbjets_lostJet",&nbjets_lostJet,"nbjets_lostJet/I");

  reducedTree.Branch("minDeltaPhiN_Luke_lostJet", &minDeltaPhiN_Luke_lostJet, "minDeltaPhiN_Luke_lostJet/F");
  reducedTree.Branch("maxDeltaPhiN_Luke_lostJet", &maxDeltaPhiN_Luke_lostJet, "maxDeltaPhiN_Luke_lostJet/F");
  reducedTree.Branch("deltaPhiN1_Luke_lostJet", &deltaPhiN1_Luke_lostJet, "deltaPhiN1_Luke_lostJet/F");
  reducedTree.Branch("deltaPhiN2_Luke_lostJet", &deltaPhiN2_Luke_lostJet, "deltaPhiN2_Luke_lostJet/F");
  reducedTree.Branch("deltaPhiN3_Luke_lostJet", &deltaPhiN3_Luke_lostJet, "deltaPhiN3_Luke_lostJet/F");

  reducedTree.Branch("minTransverseMETSignificance_lostJet", &minTransverseMETSignificance_lostJet, "minTransverseMETSignificance_lostJet/F");
  reducedTree.Branch("maxTransverseMETSignificance_lostJet", &maxTransverseMETSignificance_lostJet, "maxTransverseMETSignificance_lostJet/F");
  reducedTree.Branch("transverseMETSignificance1_lostJet", &transverseMETSignificance1_lostJet, "transverseMETSignificance1_lostJet/F");
  reducedTree.Branch("transverseMETSignificance2_lostJet", &transverseMETSignificance2_lostJet, "transverseMETSignificance2_lostJet/F");
  reducedTree.Branch("transverseMETSignificance3_lostJet", &transverseMETSignificance3_lostJet, "transverseMETSignificance3_lostJet/F");

  reducedTree.Branch("nLostJet", &nLostJet, "nLostJet/F");

  reducedTree.Branch("nTracks1",&nTracks1,"nTracks1/I");
  reducedTree.Branch("nTracks10",&nTracks10,"nTracks10/I");

  const int nskimCounter=8;
  TH1D skimCounter("skimCounter","[nevents] ntotal nselected npasstrig nfailtrig npassjson nfailjson npassST nfailST [ntotalgen]",nskimCounter,0,nskimCounter);

  const Long64_t nevents = chainA->GetEntries();
  const Long64_t neventsB = chainB->GetEntries();
  assert(nevents==neventsB);

  //store some extra info just in case....
  skimCounter.SetBinContent(0,nevents);
  skimCounter.SetBinContent(nskimCounter+1,getNEventsGenerated());

  startTimer();
  // ~~~~ now start the real event loop
  for(Long64_t entry=0; entry < nevents; ++entry) {
    chainB->GetEntry(entry);
    chainA->GetEntry(entry);
    //    calledThisEvent_=false;

    skimCounter.Fill(0); //ntotal

    //some output to watch while it's running
    if (entry==0) {
      if ( !isSampleRealData() ) {
	cout << "MC xsec: " << getCrossSection() << endl;
      }
      else{
	cout << "This is data!"<< endl;
      }
      cout << "Running over n out of total generated = "<< nevents<<" / "<<getNEventsGenerated() << endl;
    }
    if (entry%100000==0 ) cout << "  entry: " << entry << ", percent done=" << (int)(entry/(double)nevents*100.)<<  endl;

    //    printDecay();

    pair<int,int> thispoint;
    if (theScanType_==kmSugra) assert(0);// FIXME CFA // thispoint=make_pair(TMath::Nint(eventlhehelperextra_m0),TMath::Nint(eventlhehelperextra_m12)) ;
    else if (theScanType_==kSMS) thispoint=getSMSmasses();
    else thispoint=make_pair(0,0);
    
    if (thispoint != lastpoint) {
      if (theScanType_==kmSugra)    cout<<"At mSugra point m0  = "<<thispoint.first<<" m12 = "<<thispoint.second<<endl;
      else  if (theScanType_==kSMS) cout<<"At SMS point m_gluino = "<<thispoint.first<<" m_LSP = "<<thispoint.second<<endl;
      lastpoint=thispoint;
      if ( theScanType_==kmSugra && scanProcessTotalsMapCTEQ.count(thispoint)==0 )	cout<<"m0 m12 = "<<thispoint.first<<" "<<thispoint.second<<" does not exist in NLO map!"<<endl;
    }

    // ~~~~ special stuff that must be done on all events (whether it passes the skim cuts or not)
    //all of it is for MC. if it was needed for data, we'd want to put the lumi mask cut out here
    float susy_pt1=0,susy_phi1=0,susy_pt2=0,susy_phi2=0;
    SUSYProcess prodprocess= (theScanType_==kmSugra ||theScanType_==kSMS) ? getSUSYProcess(susy_pt1,susy_phi1,susy_pt2,susy_phi2) : NotFound; //don't do LM here, at least not for now
    SUSY_process = int(prodprocess);
    m0 = thispoint.first;
    m12=thispoint.second;

    float susy_px1 = susy_pt1 * cos(susy_phi1);
    float susy_py1 = susy_pt1 * sin(susy_phi1);
    float susy_px2 = susy_pt2 * cos(susy_phi2);
    float susy_py2 = susy_pt2 * sin(susy_phi2);
    float susy_px = susy_px1 + susy_px2;
    float susy_py = susy_py1 + susy_py2;
    SUSY_recoilPt = sqrt(susy_px*susy_px + susy_py*susy_py);

    if (theScanType_==kmSugra) {
      assert(0); //FIXME CFA
      if ( scanProcessTotalsMapCTEQ.count(thispoint) ) {
	//do the new 2D maps
	for (int ipdf=0 ; ipdf<45; ipdf++) {
	  pdfWeightsCTEQ[ipdf] = checkPdfWeightSanity(1 /* geneventinfoproducthelper1.at(ipdf).pdfweight */ ); //FIXME CFA
	  scanProcessTotalsMapCTEQ[thispoint]->SetBinContent( int(prodprocess),  ipdf,
							      scanProcessTotalsMapCTEQ[thispoint]->GetBinContent( int(prodprocess),  ipdf) 
							      + pdfWeightsCTEQ[ipdf] );
	}
	for (int ipdf=0 ; ipdf<41; ipdf++) { 
	  pdfWeightsMSTW[ipdf] = checkPdfWeightSanity(1 /*geneventinfoproducthelper2.at(ipdf).pdfweight*/); //FIXME CFA
	  scanProcessTotalsMapMSTW[thispoint]->SetBinContent( int(prodprocess),  ipdf,
							      scanProcessTotalsMapMSTW[thispoint]->GetBinContent( int(prodprocess),  ipdf) 
							      + pdfWeightsMSTW[ipdf] );
	}
	for (int ipdf=0 ; ipdf<100; ipdf++) { 
	  pdfWeightsNNPDF[ipdf] = checkPdfWeightSanity(1 /*geneventinfoproducthelper.at(ipdf).pdfweight*/); //FIXME CFA
	  scanProcessTotalsMapNNPDF[thispoint]->SetBinContent( int(prodprocess),  ipdf,
							       scanProcessTotalsMapNNPDF[thispoint]->GetBinContent( int(prodprocess),  ipdf) 
							       + pdfWeightsNNPDF[ipdf] );
	}
      }
      else 	continue; // skip this event
    }
    else if (theScanType_==kSMS) {
      //increment a 2d histogram of mGL, mLSP
      //we know 10k were generated everywhere, but what if we have failed jobs?
      // -- this scheme is not compatible with skim, unfortunately
      scanSMSngen->Fill(m0,m12);

      //do the new 2D maps as well
      //the v66 ntuples have a problem -- i didn't store the central pdf weight
      /* 	 to handle this:
	 make counters (ipdf) start from 1 instead of 0
	 access pdfweights arrays at ipdf-1 instead of ipdf
	 calculate an average pdf weight
	 after each loop, set the values for the [0] elements by hand
      */

      /* temporary hack for v65 -- store 1 everywhere instead of a real number */

      double av=0; //v66 kludge
      for (int ipdf=1 ; ipdf<45; ipdf++) {
	pdfWeightsCTEQ[ipdf] = checkPdfWeightSanity(pdfweights_cteq->at(ipdf-1));
	av += pdfWeightsCTEQ[ipdf]; //v66 kludge
	scanProcessTotalsMapCTEQ[thispoint]->SetBinContent( int(prodprocess),  ipdf,
							      scanProcessTotalsMapCTEQ[thispoint]->GetBinContent( int(prodprocess),  ipdf) 
							    + pdfWeightsCTEQ[ipdf] );
      }

      av /= 44; //v66 kludge
      pdfWeightsCTEQ[0] = av; //v66 kludge
      scanProcessTotalsMapCTEQ[thispoint]->SetBinContent( int(prodprocess),  0, //v66 kludge
							  scanProcessTotalsMapCTEQ[thispoint]->GetBinContent( int(prodprocess),  0) 
							  + pdfWeightsCTEQ[0] );
      av=0;   //v66 kludge    
      for (int ipdf=1 ; ipdf<41; ipdf++) { 
	pdfWeightsMSTW[ipdf] = checkPdfWeightSanity(pdfweights_mstw->at(ipdf-1));
	av += pdfWeightsMSTW[ipdf]; //v66 kludge
	scanProcessTotalsMapMSTW[thispoint]->SetBinContent( int(prodprocess),  ipdf,
							    scanProcessTotalsMapMSTW[thispoint]->GetBinContent( int(prodprocess),  ipdf) 
							    + pdfWeightsMSTW[ipdf] );
      }
      av/=40; //v66 kludge
      pdfWeightsMSTW[0] = av; //v66 kludge
      scanProcessTotalsMapMSTW[thispoint]->SetBinContent( int(prodprocess),  0, //v66 kludge
							  scanProcessTotalsMapMSTW[thispoint]->GetBinContent( int(prodprocess),  0) 
							  + pdfWeightsMSTW[0] );
      av=0; //v66 kludge
      for (int ipdf=1 ; ipdf<100; ipdf++) { 
	pdfWeightsNNPDF[ipdf] = checkPdfWeightSanity(pdfweights_nnpdf->at(ipdf-1));
	av += pdfWeightsNNPDF[ipdf]; //v66 kludge
	scanProcessTotalsMapNNPDF[thispoint]->SetBinContent( int(prodprocess),  ipdf,
							     scanProcessTotalsMapNNPDF[thispoint]->GetBinContent( int(prodprocess),  ipdf) 
							     + pdfWeightsNNPDF[ipdf] );
      }
      av/=100; //v66 kludge
      pdfWeightsNNPDF[0] = av; //v66 kludge
      scanProcessTotalsMapNNPDF[thispoint]->SetBinContent( int(prodprocess),  0, //v66 kludge
							  scanProcessTotalsMapNNPDF[thispoint]->GetBinContent( int(prodprocess),  0) 
							  + pdfWeightsNNPDF[0] );

    }
    else if (theScanType_==kNotScan && sampleIsSignal_) {
      //do the new 2D maps as well
      for (int ipdf=0 ; ipdf<45; ipdf++) {
	pdfWeightsCTEQ[ipdf] = checkPdfWeightSanity(1 /*geneventinfoproducthelper1.at(ipdf).pdfweight*/); //FIXME CFA
	  scanProcessTotalsMapCTEQ[thispoint]->SetBinContent( int(prodprocess),  ipdf,
							      scanProcessTotalsMapCTEQ[thispoint]->GetBinContent( int(prodprocess),  ipdf) 
							      + pdfWeightsCTEQ[ipdf] );
      }
      for (int ipdf=0 ; ipdf<41; ipdf++) { 
	pdfWeightsMSTW[ipdf] = checkPdfWeightSanity(1 /*geneventinfoproducthelper2.at(ipdf).pdfweight*/); //FIXME CFA
	scanProcessTotalsMapMSTW[thispoint]->SetBinContent( int(prodprocess),  ipdf,
							    scanProcessTotalsMapMSTW[thispoint]->GetBinContent( int(prodprocess),  ipdf) 
							    + pdfWeightsMSTW[ipdf] );
      }
      for (int ipdf=0 ; ipdf<100; ipdf++) { 
	pdfWeightsNNPDF[ipdf] = checkPdfWeightSanity(1 /*geneventinfoproducthelper.at(ipdf).pdfweight*/); //FIXME CFA
	scanProcessTotalsMapNNPDF[thispoint]->SetBinContent( int(prodprocess),  ipdf,
							     scanProcessTotalsMapNNPDF[thispoint]->GetBinContent( int(prodprocess),  ipdf) 
							     + pdfWeightsNNPDF[ipdf] );
      }
    }


    //very loose skim for reducedTrees
    //now that we have a skim in cfA with a fairly hard MET cut, I will remove all kinematic skims here
    ST = getST(30,10); //ST is always bigger than HT; reducing the jet cut to 30 will only make it larger
    MET=getMET();
    STeff = ST + MET; //obviously STeff is always bigger than ST
    HT=getHT();

    bool passAnyTrigger=passHLT(triggerlist,true); //evaluate all trigger decisions
    passHLT(triggerlist_mc,false); //re-evaluate, storing real results for mc (will be a duplicate for data)

    //idea -- we should do a dataset-aware OR of the two groups of physics triggers here, to allow easy combination of the MET and HTMHT datasets
    bool passMETtriggers = triggerlist["DiCentralPFJet50_PFMET80"].pass || triggerlist["DiCentralPFNoPUJet50_PFMETORPFMETNoMu80"].pass;
    bool passHTMHTtriggers = triggerlist["PFHT350_PFMET100"].pass || triggerlist["PFNoPUHT350_PFMET100"].pass;
    if (sampleName_.BeginsWith("HT_Run2012A") || sampleName_.BeginsWith("HTMHT_Run2012B")|| sampleName_.BeginsWith("HTMHT_Run2012C")) {
      //from MHT dataset, use only events that don't pass the MET triggers
      cutTrigger =   	passHTMHTtriggers && 	(!passMETtriggers);
    }
    else if (sampleName_.BeginsWith("MET_Run2012")) {
      cutTrigger = passMETtriggers;
    }
    else { //MC etc
      cutTrigger = passHTMHTtriggers || passMETtriggers;
    }

    //now we do a super-OR of MET/HTMHT/JetHT
    //i kind of coded this upside down; the logic is:
    //    take all events that pass MET triggers from MET dataset
    //    for those that don't pass MET triggers, look in JetHT dataset
    //    for those that don't pass either of the above, look in HTMHT
    //with an additional complication from the fact that the dataset breakdown is a bit different in Run2012A
    bool passHTtriggers = triggerlist["PFHT650"].pass || triggerlist["PFNoPUHT650"].pass;
    //if (sampleName_.BeginsWith("HT_Run2012A") ) { //contains both HTMHT and HT triggers
    //  cutTrigger2 = (passHTMHTtriggers||passHTtriggers) && (!passMETtriggers);
    //}
    //else if ( sampleName_.BeginsWith("HTMHT_Run2012B")|| sampleName_.BeginsWith("HTMHT_Run2012C")) { //contains HTMHT
    //  cutTrigger2 = passHTMHTtriggers && (!( passMETtriggers||passHTtriggers));
    //}
    //else if ( sampleName_.BeginsWith("JetHT_Run2012B")|| sampleName_.BeginsWith("JetHT_Run2012C")) { //contains HT
    //  //pick up events here that have failed the MET triggers
    //  cutTrigger2 = passHTtriggers && (!passMETtriggers);
    //}
    if (sampleName_.BeginsWith("HT_Run2012A") ) { //contains both HTMHT and HT triggers
      cutTrigger2 = (passHTMHTtriggers||passHTtriggers) && (!passMETtriggers);
    }
    else if ( sampleName_.BeginsWith("JetHT_Run2012B")|| sampleName_.BeginsWith("JetHT_Run2012C")) { //contains HT
      //pick up events here for HT that have failed the MET and HTMHT triggers
      cutTrigger2 = passHTtriggers && (!(passMETtriggers||passHTMHTtriggers));
    }
    else if ( sampleName_.BeginsWith("HTMHT_Run2012B")|| sampleName_.BeginsWith("HTMHT_Run2012C")) { //contains HTMHT
      cutTrigger2 = passHTMHTtriggers && ( !passMETtriggers );
    }
    else if (sampleName_.BeginsWith("MET_Run2012")) { //MET triggers
      //always take events in MET dataset passing the trigger
      cutTrigger2 = passMETtriggers;
    }
    else { //MC etc
      cutTrigger2 = passHTMHTtriggers||passMETtriggers||passHTtriggers;
    }

  //these variables are basically not needed anymore, but leave them in for now
    pass_utilityHLT = triggerlist["PFHT350"].pass; //shortcut name for HT-only trigger
    prescaleUtilityHLT = triggerlist["PFHT350"].prescale;
    versionUtilityHLT = triggerlist["PFHT350"].version;

    isRealData = isSampleRealData();
    runNumber = getRunNumber();
    lumiSection = getLumiSection();
    eventNumber = getEventNumber();
    //stuff that we only do on data....
    bool      passJSON=true;
    //check the JSON file for data
    if (isRealData)   {
      // use appropriate rereco json as needed
      if ( (runNumber>=190456) && (runNumber<=196531) ) passJSON=  inJSON(VRunLumiJuly13,runNumber,lumiSection);
      else if ( (runNumber>=198022) && (runNumber<=198913) ) passJSON=  inJSON(VRunLumiAug24,runNumber,lumiSection);
      // otherwise use prompt reco json
      else passJSON=  inJSON(VRunLumi,runNumber,lumiSection);    
    }

    //keep some tallies of why we failed the skim
    if (passAnyTrigger)  skimCounter.Fill(2); else skimCounter.Fill(3);
    if (passJSON)  skimCounter.Fill(4); else skimCounter.Fill(5);
    if ( true ) skimCounter.Fill(6); else skimCounter.Fill(7); //used to look at the STeff skim cut here

    if (passJSON && passAnyTrigger ) { //very loose skim cut
      //begin main block of filling reducedTree variables
      skimCounter.Fill(1); //nselected

      eventweight = getWeight(nevents);
      eventweight2 = getWeight( getNEventsGenerated());
      
      if (theScanType_!=kSMS) {
	scanCrossSection = getScanCrossSection(prodprocess,"");
	scanCrossSectionPlus = getScanCrossSection(prodprocess,"Plus");
	scanCrossSectionMinus = getScanCrossSection(prodprocess,"Minus");
      }
      else {
	scanCrossSection = 1;
	scanCrossSectionPlus = scanCrossSection;
	scanCrossSectionMinus = scanCrossSection; //there are no cross section errors for SMS
      }

      //if we are running over ttbar, fill info on decay mode
      if (sampleName_.Contains("TTJets")) { //revived for cfA, using the old jmt-ntuples coding scheme
	ttbarDecayCode = getTTbarDecayType();
	  
      }
      //somewhat experimental code; maybe this is a backwards way to do it (look for tau match then store stuff)
      //anyway this will give me a first look
      int tauindex=findJetMatchGenTau();
      if (tauindex>=0) { //crude tau id variables
	tauMatch_jetpt = jets_AK5PF_pt->at(tauindex);
	tauMatch_chhadmult = jets_AK5PF_chg_Mult->at(tauindex);
	tauMatch_jetcsv = jets_AK5PF_btag_secVertexCombined->at(tauindex);
	tauMatch_MT = 2*getJetPt(tauindex)*getMET()*(1 - cos( getDeltaPhi( jets_AK5PF_phi->at(tauindex) , getMETphi())));
      }
      else {
	tauMatch_jetpt=-1;
	tauMatch_chhadmult=-1;
	tauMatch_jetcsv = -1;
	tauMatch_MT = -1;
      }
/*
      //single-top tW-channel has two W's (one from a top)
      if (sampleName_.Contains("T_TuneZ2_tW-channel") || sampleName_.Contains("Tbar_TuneZ2_tW-channel")) {
	int W1, W1daughter, W2, W2daughter;	
	getWDecayType(W2decayType, W2, W2daughter, false); //FIXME CFA
	decayType = getTTbarDecayType(W1decayType, W2decayType, W1, W1daughter, W2, W2daughter,true); //FIXME CFA
      }
      //if we are running of w+jets fill info on decay mode
      //single-top s/t-channel also has one W (from a top)
      if (sampleName_.Contains("WJetsToLNu") 
	  || sampleName_.Contains("T_TuneZ2_s-channel") || sampleName_.Contains("Tbar_TuneZ2_s-channel") 
	  || sampleName_.Contains("T_TuneZ2_t-channel") || sampleName_.Contains("Tbar_TuneZ2_t-channel")){
	int W1, W1daughter;
	decayType = getWDecayType(W1decayType, W1, W1daughter, false); //FIXME CFA
      }
*/

      HT30=getHT(30);
      hltHTeff = getHLTHTeff(HT);
      PUweight =  puReweightIs1D ? getPUWeight(*LumiWeights) : 1;// FIXME CFA getPUWeight(*LumiWeights3D);
      PUweightSystVar =  puReweightIs1D ? getPUWeight(*LumiWeightsSystVar) : 1;// FIXME CFA getPUWeight(*LumiWeights3D);
      //      pfmhtweight = getPFMHTWeight();

      cutPV = passPV();

      calculateTagProb(prob0,probge1,prob1,probge2,prob2,probge3);

      calculateTagProb(prob0_HFplus,probge1_HFplus,prob1_HFplus,probge2_HFplus,prob2_HFplus,probge3_HFplus,1,1,1,kHFup);
      calculateTagProb(prob0_HFminus,probge1_HFminus,prob1_HFminus,probge2_HFminus,prob2_HFminus,probge3_HFminus,1,1,1,kHFdown);

      calculateTagProb(prob0_LFplus,probge1_LFplus,prob1_LFplus,probge2_LFplus,prob2_LFplus,probge3_LFplus,1,1,1,kLFup);
      calculateTagProb(prob0_LFminus,probge1_LFminus,prob1_LFminus,probge2_LFminus,prob2_LFminus,probge3_LFminus,1,1,1,kLFdown);

      pass_utilityPrescaleModuleHLT = passUtilityPrescaleModuleHLT();

      nGoodPV = countGoodPV();

      SUSY_nb = sampleIsSignal_ ? getSUSYnb() : 0;
      //      bjetSumSUSY[thispoint] += SUSY_nb;
      //if(SUSY_process==NotFound) cout<<"SUSY_nb = "<<SUSY_nb<<endl;

      //count jets
      njets = nGoodJets();
      njets30 = nGoodJets30();

      //look for jj resonances
      jjResonanceFinder(mjj1,mjj2, nCorrectRecoStop);
      mjjdiff=fabs(mjj1-mjj2);
      jjResonanceFinder5(mjj1_5,mjj2_5);
      mjjdiff_5=fabs(mjj1_5-mjj2_5);

      //count b jets
      ntruebjets = nTrueBJets();
      nbjets = nGoodBJets();
      nbjets30 = nGoodBJets(30);
      nbjetsSSVM = nGoodBJets( kSSVM);
      nbjetsSSVHPT = nGoodBJets( kSSVHPT);
      nbjetsTCHET = nGoodBJets( kTCHET);
      nbjetsTCHPT = nGoodBJets( kTCHPT);
      nbjetsTCHPM = nGoodBJets( kTCHPM);
      nbjetsCSVM = nGoodBJets( kCSVM);
      nbjetsCSVL = nGoodBJets( kCSVL);

      hadronicTopFinder_DeltaR(mjjb1,mjjb2,topPT1,topPT2);

      CSVout1=bjetCSVOfN(1);
      CSVout2=bjetCSVOfN(2);
      CSVout3=bjetCSVOfN(3);
      minDeltaPhiAllb30=getMinDeltaPhiMET30(99,true);
      deltaPhib1=getDeltaPhiMET(1,30,true);
      deltaPhib2=getDeltaPhiMET(2,30,true); 
      deltaPhib3=getDeltaPhiMET(3,30,true);
      minDeltaPhiMETMuonsAll=getMinDeltaPhiMETMuons(99);

      //count leptons
      nElectrons = countEle();
      nMuons = countMu();
      nElectrons5 = countEle(5);
      nMuons5 = countMu(5);
      nElectrons15 = countEle(15);
      nMuons15 = countMu(15);
      nElectrons20 = countEle(20);
      nMuons20 = countMu(20);

      bestZmass=getBestZCandidate(20,10);

      nTaus = countTau();
      MHT=getMHT();
      METphi = getMETphi();

      METsig = pfmets_fullSignif;
      METsig00 = pfmets_fullSignifCov00;
      METsig10 = pfmets_fullSignifCov10;
      METsig11 = pfmets_fullSignifCov11;

      caloMET = mets_AK5_et->at(0);
      caloMETphi = mets_AK5_phi->at(0);
      getUncorrectedMET(rawPFMET, rawPFMETphi);

      minDeltaPhi = getMinDeltaPhiMET(3);
      minDeltaPhiAll = getMinDeltaPhiMET(99);
      minDeltaPhiAll30 = getMinDeltaPhiMET30(99);
      minDeltaPhi30_eta5_noIdAll = getMinDeltaPhiMET30_eta5_noId(99);
      maxDeltaPhi = getMaxDeltaPhiMET(3);
      maxDeltaPhiAll = getMaxDeltaPhiMET(99);
      maxDeltaPhiAll30 = getMaxDeltaPhiMET30(99);
      maxDeltaPhi30_eta5_noIdAll = getMaxDeltaPhiMET30_eta5_noId(99);
      
      int j_index_1,j_index_2;
      deltaR_bestTwoCSV = deltaRBestTwoBjets(j_index_1,j_index_2);
      mjj_bestTwoCSV = (j_index_1>=0 && j_index_2>=0) ? calc_mNj(j_index_1,j_index_2) : -1;
      deltaR_closestB=deltaRClosestTwoBjets(j_index_1,j_index_2);
      mjj_closestB= (j_index_1>=0 && j_index_2>=0) ? calc_mNj(j_index_1,j_index_2) : -1;

      if (isSampleQCD()) { //don't fill except for QCD MC
      minDeltaPhiN_otherPt10                      = getMinDeltaPhiMETN(3,50,2.4,true, 10,2.4,true, false,false);
      minDeltaPhiN_otherPt20                      = getMinDeltaPhiMETN(3,50,2.4,true, 20,2.4,true, false,false);
      minDeltaPhiN_otherPt30                      = getMinDeltaPhiMETN(3,50,2.4,true, 30,2.4,true, false,false);
      minDeltaPhiN_otherPt40                      = getMinDeltaPhiMETN(3,50,2.4,true, 40,2.4,true, false,false);
      minDeltaPhiN_otherPt50                      = getMinDeltaPhiMETN(3,50,2.4,true, 50,2.4,true, false,false);
      minDeltaPhiN_DJR_otherPt10                  = getMinDeltaPhiMETN(3,50,2.4,true, 10,2.4,true,true,false);
      minDeltaPhiN_DJR_otherPt20                  = getMinDeltaPhiMETN(3,50,2.4,true, 20,2.4,true,true,false);
      minDeltaPhiN_DJR_otherPt30                  = getMinDeltaPhiMETN(3,50,2.4,true, 30,2.4,true,true,false);
      minDeltaPhiN_DJR_otherPt40                  = getMinDeltaPhiMETN(3,50,2.4,true, 40,2.4,true,true,false);
      minDeltaPhiN_DJR_otherPt50                  = getMinDeltaPhiMETN(3,50,2.4,true, 50,2.4,true,true,false);
            
      minDeltaPhiN_withLepton = getMinDeltaPhiMETN(3,50,2.4,true,30,2.4,true,false,false,true);
      }

      int badjet;
      deltaPhiStar = getDeltaPhiStar(badjet);
      deltaPhiStar_badjet_pt = (badjet>=0) ? jets_AK5PF_pt->at(badjet) : -1;
      deltaPhiStar_badjet_phi = (badjet>=0) ? jets_AK5PF_phi->at(badjet) : -1;
      deltaPhiStar_badjet_eta = (badjet>=0) ? jets_AK5PF_eta->at(badjet) : -1;

      if (isSampleQCD()) { //don't fill except for QCD MC
      maxJetMis=getMaxJetMis(1,3,50);
      max2JetMis=getMaxJetMis(2,3,50);
      maxJetMisAll30=getMaxJetMis(1,99,30);
      max2JetMisAll30=getMaxJetMis(2,99,30);
      maxJetFracMis=getMaxJetFracMis(1,3,50);
      max2JetFracMis=getMaxJetFracMis(2,3,50);
      maxJetFracMisAll30=getMaxJetFracMis(1,99,30);
      max2JetFracMisAll30=getMaxJetFracMis(2,99,30);
      deltaPhiMETJetMaxMis = getDeltaPhiMETJetMaxMis(50);
      deltaPhiMETJetMaxMis30 = getDeltaPhiMETJetMaxMis(30);
      }

      maxJetMis_chosenJet=1;
      float maxJetMis2_tmp =getMaxJetMis(1,2,50); 
      if(maxJetMis2_tmp > getMaxJetMis(1,1,50)) maxJetMis_chosenJet=2;
      if(maxJetMis > maxJetMis2_tmp) maxJetMis_chosenJet=3;

      maxJetFracMis_chosenJet=1;
      float maxJetFracMis2_tmp =getMaxJetFracMis(1,2,50); 
      if(maxJetFracMis2_tmp > getMaxJetFracMis(1,1,50)) maxJetFracMis_chosenJet=2;
      if(maxJetFracMis > maxJetFracMis2_tmp) maxJetFracMis_chosenJet=3;

      minDeltaPhiMetTau = getMinDeltaPhiMETTaus();

      sumDeltaPhi = maxDeltaPhi + minDeltaPhi;
      diffDeltaPhi = maxDeltaPhi - minDeltaPhi;
            
      minDeltaPhiN = getMinDeltaPhiMETN(3);


      deltaPhi1 = getDeltaPhiMET(1,50,false);
      deltaPhi2 = getDeltaPhiMET(2,50,false);
      deltaPhi3 = getDeltaPhiMET(3,50,false);

      //a hack that could be put into getMinDeltaPhiMETN with some work
      if(deltaPhiN1<=deltaPhiN2 && deltaPhiN1<=deltaPhiN3) minDeltaPhiN_chosenJet = 1;
      else if(deltaPhiN2<=deltaPhiN1 && deltaPhiN2<=deltaPhiN3) minDeltaPhiN_chosenJet = 2;
      else if(deltaPhiN3<=deltaPhiN1 && deltaPhiN3<=deltaPhiN2) minDeltaPhiN_chosenJet = 3;
      else minDeltaPhiN_chosenJet=0;

      if(deltaPhi1<=deltaPhi2 && deltaPhi1<=deltaPhi3) minDeltaPhi_chosenJet = 1;
      else if(deltaPhi2<=deltaPhi1 && deltaPhi2<=deltaPhi3) minDeltaPhi_chosenJet = 2;
      else if(deltaPhi3<=deltaPhi1 && deltaPhi3<=deltaPhi2) minDeltaPhi_chosenJet = 3;
      else minDeltaPhi_chosenJet=0;

      deltaT1 = getDeltaPhiMETN_deltaT( getNthGoodJet(0,50,2.4,true) );
      deltaT2 = getDeltaPhiMETN_deltaT( getNthGoodJet(1,50,2.4,true) );
      deltaT3 = getDeltaPhiMETN_deltaT( getNthGoodJet(2,50,2.4,true) );
      minDeltaPhiN_deltaT = getDeltaPhiMETN_deltaT( getNthGoodJet((unsigned int)(minDeltaPhiN_chosenJet-1),50,2.4,true) );//

      hltMHTeff = getHLTMHTeff(MET, HT, nElectrons, nMuons, minDeltaPhiN);
//       double effUp, effDown;
//       hltMHTeffBNN = getHLTMHTeffBNN(MET, HT, nElectrons, nMuons, minDeltaPhiN, effUp, effDown);
//       hltMHTeffBNNUp = effUp;
//       hltMHTeffBNNDown = effDown;

      //alternatives for testing of qcd methods
      if (isSampleQCD() ) {
      minDeltaPhiN_otherEta5                      = getMinDeltaPhiMETN(3,50,2.4,true, 30,5,true, false,false);
      minDeltaPhiN_otherEta5idNo                  = getMinDeltaPhiMETN(3,50,2.4,true, 30,5,false,false,false);
      minDeltaPhiN_mainPt30_otherEta5idNo         = getMinDeltaPhiMETN(3,30,2.4,true, 30,5,false,false,false);
      minDeltaPhiN_mainPt30Eta5_otherEta5idNo     = getMinDeltaPhiMETN(3,30,5,  true, 30,5,false,false,false);
      minDeltaPhiN_mainPt30Eta5idNo_otherEta5idNo = getMinDeltaPhiMETN(3,30,5,  false,30,5,false,false,false);

      minDeltaPhiK                                = getMinDeltaPhiMETN(3,50,2.4,true, 30,2.4,true, false,true); //K is for Keith
      minDeltaPhiK_otherEta5                      = getMinDeltaPhiMETN(3,50,2.4,true, 30,5,  true, false,true);
      minDeltaPhiK_otherEta5idNo                  = getMinDeltaPhiMETN(3,50,2.4,true, 30,5,  false,false,true);
      minDeltaPhiK_mainPt30_otherEta5idNo         = getMinDeltaPhiMETN(3,30,2.4,true, 30,5,  false,false,true);
      minDeltaPhiK_mainPt30Eta5_otherEta5idNo     = getMinDeltaPhiMETN(3,30,5,  true, 30,5,  false,false,true);
      minDeltaPhiK_mainPt30Eta5idNo_otherEta5idNo = getMinDeltaPhiMETN(3,30,5,  false,30,5,  false,false,true);
      
      minDeltaPhiN_DJR           = getMinDeltaPhiMETN(3,50,2.4,true,30,2.4,true,true,false); //DJR is for data jet resolution
      minDeltaPhiN_DJR_otherEta5 = getMinDeltaPhiMETN(3,50,2.4,true,30,5,  true,true,false);
      minDeltaPhiK_DJR           = getMinDeltaPhiMETN(3,50,2.4,true,30,2.4,true,true,true );
      minDeltaPhiK_DJR_otherEta5 = getMinDeltaPhiMETN(3,50,2.4,true,30,5,  true,true,true );
      } //don't calculate for non-qcd samples

      minDeltaPhiN_asin = getMinDeltaPhiMETN(3,50,2.4,true,30,2.4,true, false,false,false,true);
      deltaPhiN_asin1 = getDeltaPhiMETN(0,50,2.4,true,30,2.4,true, false,false,true);
      deltaPhiN_asin2 = getDeltaPhiMETN(1,50,2.4,true,30,2.4,true, false,false,true);
      deltaPhiN_asin3 = getDeltaPhiMETN(2,50,2.4,true,30,2.4,true, false,false,true);

      minDeltaPhiN_Luke = getMinDeltaPhiNMET(3);
      maxDeltaPhiN_Luke = getMaxDeltaPhiNMET(3);
      deltaPhiN1_Luke = getDeltaPhiNMET(0);
      deltaPhiN2_Luke = getDeltaPhiNMET(1);
      deltaPhiN3_Luke = getDeltaPhiNMET(2);
      
      if (isSampleQCD() ) { //don't fill for non-qcd (save space)
      passBadECAL_METphi53_n10_s12 = passBadECALFilter("METphi",0.5,"deadCell",0.3,10,12);
      worstMisA_badECAL_METphi53_n10_s12 = passBadECALFilter_worstMis("METphi",0.5,"deadCell",0.3,10,12,"abs");
      worstMisF_badECAL_METphi53_n10_s12 = passBadECALFilter_worstMis("METphi",0.5,"deadCell",0.3,10,12,"frac");
      
      passBadECAL_METphi33_n10_s12 = passBadECALFilter("METphi",0.3,"deadCell",0.3,10,12);
      worstMisA_badECAL_METphi33_n10_s12 = passBadECALFilter_worstMis("METphi",0.3,"deadCell",0.3,10,12,"abs");
      worstMisF_badECAL_METphi33_n10_s12 = passBadECALFilter_worstMis("METphi",0.3,"deadCell",0.3,10,12,"frac");
      
      passBadECAL_METphi52_n10_s12 = passBadECALFilter("METphi",0.5,"deadCell",0.2,10,12);
      worstMisA_badECAL_METphi52_n10_s12 = passBadECALFilter_worstMis("METphi",0.5,"deadCell",0.2,10,12,"abs");
      worstMisF_badECAL_METphi52_n10_s12 = passBadECALFilter_worstMis("METphi",0.5,"deadCell",0.2,10,12,"frac");
      
      passBadECAL_METphi53_n5_s12 = passBadECALFilter("METphi",0.5,"deadCell",0.3,5,12);
      passBadECAL_METphi53_crack = passBadECALFilter("METphi",0.5,"crackEBEE",0.3,0,0);
      worstMisA_badECAL_METphi53_crack = passBadECALFilter_worstMis("METphi",0.5,"crackEBEE",0.3,0,0,"abs");
      worstMisF_badECAL_METphi53_crack = passBadECALFilter_worstMis("METphi",0.5,"crackEBEE",0.3,0,0,"frac");
      passBadECAL_METphi52_crack = passBadECALFilter("METphi",0.5,"crackEBEE",0.2,0,0);
      worstMisA_badECAL_METphi52_crack = passBadECALFilter_worstMis("METphi",0.5,"crackEBEE",0.2,0,0,"abs");
      worstMisF_badECAL_METphi52_crack = passBadECALFilter_worstMis("METphi",0.5,"crackEBEE",0.2,0,0,"frac");

      passBadECAL_METphi53_crackF = passBadECALFilter("METphi",0.5,"crackEEEF",0.3,0,0);
      worstMisA_badECAL_METphi53_crackF = passBadECALFilter_worstMis("METphi",0.5,"crackEEEF",0.3,0,0,"abs");
      worstMisF_badECAL_METphi53_crackF = passBadECALFilter_worstMis("METphi",0.5,"crackEEEF",0.3,0,0,"frac");
      passBadECAL_METphi52_crackF = passBadECALFilter("METphi",0.5,"crackEEEF",0.2,0,0);
      worstMisA_badECAL_METphi52_crackF = passBadECALFilter_worstMis("METphi",0.5,"crackEEEF",0.2,0,0,"abs");
      worstMisF_badECAL_METphi52_crackF = passBadECALFilter_worstMis("METphi",0.5,"crackEEEF",0.2,0,0,"frac");
  }

      minTransverseMETSignificance = getMinTransverseMETSignificance(3);
      maxTransverseMETSignificance = getMaxTransverseMETSignificance(3);
      transverseMETSignificance1 = getTransverseMETSignificance(0);
      transverseMETSignificance2 = getTransverseMETSignificance(1);
      transverseMETSignificance3 = getTransverseMETSignificance(2);
      
      MT_Wlep = getMT_Wlep();
      MT_Wlep5 = getMT_Wlep(5);
      MT_Wlep15 = getMT_Wlep(15);
      
      MT_b = getMT_bMET();

      calcTopDecayVariables(  wMass, topMass, wCosHel, topCosHel);
      

      jetpt1 = jetPtOfN(1);
      jetgenpt1 = jetGenPtOfN(1);
      jetphi1 = jetPhiOfN(1);
      jetgenphi1 = jetGenPhiOfN(1);
      jeteta1 = jetEtaOfN(1);
      jetgeneta1 = jetGenEtaOfN(1);
      jetenergy1 = jetEnergyOfN(1);
      jetflavor1 = jetFlavorOfN(1);
      jetchargedhadronfrac1 = jetChargedHadronFracOfN(1);
      jetchargedhadronmult1 = jetChargedHadronMultOfN(1);
      jetneutralhadronmult1 = jetNeutralHadronMultOfN(1);
      jetmumult1 = jetMuMultOfN(1);

      jetpt2 = jetPtOfN(2);
      jetgenpt2 = jetGenPtOfN(2); 
      jetphi2 = jetPhiOfN(2);
      jetgenphi2 = jetGenPhiOfN(2);
      jeteta2 = jetEtaOfN(2);
      jetgeneta2 = jetGenEtaOfN(2);
      jetenergy2 = jetEnergyOfN(2);
      jetflavor2 = jetFlavorOfN(2);
      jetchargedhadronfrac2 = jetChargedHadronFracOfN(2);
      jetchargedhadronmult2 = jetChargedHadronMultOfN(2);
      jetneutralhadronmult2 = jetNeutralHadronMultOfN(2);
      jetmumult2 = jetMuMultOfN(2);

      jetpt3 = jetPtOfN(3);
      jetgenpt3 = jetGenPtOfN(3); 
      jetphi3 = jetPhiOfN(3);
      jetgenphi3 = jetGenPhiOfN(3);
      jeteta3 = jetEtaOfN(3);
      jetgeneta3 = jetGenEtaOfN(3);
      jetenergy3 = jetEnergyOfN(3);      
      jetflavor3 = jetFlavorOfN(3);
      jetchargedhadronfrac3 = jetChargedHadronFracOfN(3);
      jetchargedhadronmult3 = jetChargedHadronMultOfN(3);
      jetneutralhadronmult3 = jetNeutralHadronMultOfN(3);
      jetmumult3 = jetMuMultOfN(3);

      bjetpt1 = bjetPtOfN(1);
      bjetphi1 = bjetPhiOfN(1);
      bjeteta1 = bjetEtaOfN(1);
      bjetenergy1 = bjetEnergyOfN(1); 
      bjetflavor1 = bjetFlavorOfN(1);
      bjetchargedhadronfrac1 = bjetChargedHadronFracOfN(1);
      bjetchargedhadronmult1 = bjetChargedHadronMultOfN(1);

      bjetpt2 = bjetPtOfN(2);
      bjetphi2 = bjetPhiOfN(2);
      bjeteta2 = bjetEtaOfN(2);
      bjetenergy2 = bjetEnergyOfN(2); 
      bjetflavor2 = bjetFlavorOfN(2);
      bjetchargedhadronfrac2 = bjetChargedHadronFracOfN(2);
      bjetchargedhadronmult2 = bjetChargedHadronMultOfN(2);

      bjetpt3 = bjetPtOfN(3);
      bjetphi3 = bjetPhiOfN(3);
      bjeteta3 = bjetEtaOfN(3);
      bjetenergy3 = bjetEnergyOfN(3);       
      bjetflavor3 = bjetFlavorOfN(3);
      bjetchargedhadronfrac3 = bjetChargedHadronFracOfN(3);
      bjetchargedhadronmult3 = bjetChargedHadronMultOfN(3);

      eleet1 = elePtOfN(1,5);
      elephi1 = elePhiOfN(1,5);
      eleeta1 = eleEtaOfN(1,5);
      elecharge1 = eleChargeOfN(1,5);
      muonpt1 = muonPtOfN(1,5);
      muonphi1 = muonPhiOfN(1,5);
      muoneta1 = muonEtaOfN(1,5);
      muoniso1 = muonIsoOfN(1,5);
      muoncharge1 = muonChargeOfN(1,5);
      muonchhadiso1 = muonChHadIsoOfN(1,5);
      muonphotoniso1 = muonPhotonIsoOfN(1,5);
      muonneutralhadiso1 = muonNeutralHadIsoOfN(1,5);

      eleet2 = elePtOfN(2,5);
      elephi2 = elePhiOfN(2,5);
      eleeta2 = eleEtaOfN(2,5);
      muonpt2 = muonPtOfN(2,5);
      muonphi2 = muonPhiOfN(2,5);
      muoneta2 = muonEtaOfN(2,5);

      taupt1 = tauPtOfN(1);     
      taueta1 = tauEtaOfN(1);     

      //i'm giving these awkward names on purpose, so that they won't be used without understanding what they do
      eleRelIso = getRelIsoForIsolationStudyEle();
      muonRelIso = getRelIsoForIsolationStudyMuon();


      recomuonpt1  = recoMuonPtOfN(1,5);
      recomuonphi1 = recoMuonPhiOfN(1,5);
      recomuoneta1 = recoMuonEtaOfN(1,5);
      recomuoniso1 = recoMuonIsoOfN(1,5);
      recomuonmindphijet1 = recoMuonMinDeltaPhiJetOfN(1,5);

      recomuonpt2  = recoMuonPtOfN(2,5);
      recomuonphi2 = recoMuonPhiOfN(2,5);
      recomuoneta2 = recoMuonEtaOfN(2,5);
      recomuoniso2 = recoMuonIsoOfN(2,5);
      recomuonmindphijet2 = recoMuonMinDeltaPhiJetOfN(2,5);

      rl=0;
      if (eleet1>muonpt1) rl= eleet1/STeff;
      else                rl= muonpt1/STeff;
      rMET = MET/STeff;

      //updated for 2012
      csctighthaloFilter = cschalofilter_decision==1 ? true : false;
      eenoiseFilter = eenoisefilter_decision==1 ? true : false;
      greedymuonFilter = greedymuonfilter_decision==1 ? true : false;
      //shoudl check if 2012 fastsim still cannot use this
      hbhenoiseFilter = (theScanType_!=kNotScan || (hbhefilter_decision==1)) ? true : false;
      inconsistentmuonFilter = inconsistentPFmuonfilter_decision==1 ? true : false;
      ra2ecaltpFilter =ecalTPfilter_decision==1 ? true:false; 
      scrapingvetoFilter = scrapingVeto_decision==1?true:false;
      trackingfailureFilter = trackingfailurefilter_decision ==1 ? true:false;

      bool hcallaser = hcallaserfilter_decision == 1 ? true:false;
      bool eebadsc = eebadscfilter_decision == 1 ? true:false;

      badjetFilter = passBadJetFilter();

      //Keith "RA2 claims eenoise overcleans in 2012; instead use jetID filter"
      //after looking at keith's code it is the new name for the old particle-based noise rejection
      PBNRcode = doPBNR();

      passCleaning = csctighthaloFilter 
	//	&& eenoiseFilter 
	&& greedymuonFilter && hbhenoiseFilter && inconsistentmuonFilter 
	&& ra2ecaltpFilter && scrapingvetoFilter && trackingfailureFilter
	&& hcallaser && eebadsc && (PBNRcode>0) && badjetFilter; //last line are new for 2012

      buggyEvent = inEventList(vrun, vlumi, vevent);
      
      //jet info for the lead jet, no matter where it is
      if (jets_AK5PF_pt->size() >0) {
	alletajetpt1 = jets_AK5PF_pt->at(0);
	alletajeteta1 = jets_AK5PF_eta->at(0);
	alletajetphi1 = jets_AK5PF_phi->at(0);
	//same as used in PBNR, with 0.9 cut
	const float rawjetenergy = jets_AK5PF_energy->at(0) * jets_AK5PF_corrFactorRaw->at(0);
	alletajetneutralhadronfrac1 = jets_AK5PF_neutralHadE->at(0)/rawjetenergy;
	alletajetneutralemfrac1 =  jets_AK5PF_neutralEmE->at(0) / rawjetenergy;
	alletajetphotonfrac1 = jets_AK5PF_photonEnergy->at(0)/rawjetenergy;
      }
      else {
	alletajetpt1 = -99;
	alletajeteta1 =-99;
	alletajetphi1 =-99;
	alletajetneutralhadronfrac1 =-99;
	alletajetneutralemfrac1 =  -99;
	alletajetphotonfrac1 = -99;
      }

      //transverse sphericity
      getSphericity(transverseSphericity_jets,false,false,50);
      getSphericity(transverseSphericity_jetsMet,true,false,50);
      getSphericity(transverseSphericity_jetsMetLeptons,true,true,50);
      getSphericity(transverseSphericity_jetsLeptons,false,true,50);

      getSphericity(transverseSphericity_jets30,false,false,30);
      getSphericity(transverseSphericity_jets30Met,true,false,30);
      getSphericity(transverseSphericity_jets30MetLeptons,true,true,30);
      getSphericity(transverseSphericity_jets30Leptons,false,true,30);

      //tracking...should really have quality cuts but i don't know the appropriate cuts
      nTracks1=0;
      nTracks10=0;
      for (unsigned int iii=0; iii<tracks_chi2->size(); iii++) {
	if (tracks_pt->at(iii) >= 1) {
	  ++nTracks1;
	  if (tracks_pt->at(iii) >= 10) {
	    ++nTracks10;
	  }
	}
      }

      //Uncomment next two lines to do thrust calculations
      //getTransverseThrustVariables(transverseThrust, transverseThrustPhi, false);
      //getTransverseThrustVariables(transverseThrustWithMET, transverseThrustWithMETPhi, true);

      //FIXME CFA
      //      changeVariables(&random,0.05,nLostJet);

      //changeVariablesGenSmearing(&random);
/* FIXME CFA
      njets_lostJet = nGoodJets();
      nbjets_lostJet = nGoodBJets();

      minDeltaPhiN_lostJet = getMinDeltaPhiMETN(3);
      deltaPhiN1_lostJet = getDeltaPhiMETN(0);
      deltaPhiN2_lostJet = getDeltaPhiMETN(1);
      deltaPhiN3_lostJet = getDeltaPhiMETN(2);

      minDeltaPhiN_Luke_lostJet = getMinDeltaPhiNMET(3);
      maxDeltaPhiN_Luke_lostJet = getMaxDeltaPhiNMET(3);
      deltaPhiN1_Luke_lostJet = getDeltaPhiNMET(0);
      deltaPhiN2_Luke_lostJet = getDeltaPhiNMET(1);
      deltaPhiN3_Luke_lostJet = getDeltaPhiNMET(2);

      minTransverseMETSignificance_lostJet = getMinTransverseMETSignificance(3);
      maxTransverseMETSignificance_lostJet = getMaxTransverseMETSignificance(3);
      transverseMETSignificance1_lostJet = getTransverseMETSignificance(0);
      transverseMETSignificance2_lostJet = getTransverseMETSignificance(1);
      transverseMETSignificance3_lostJet = getTransverseMETSignificance(2);
      resetVariables();
      */
      //Fill the tree
      reducedTree.Fill();
      

    } //end of reduced tree skim
  } //end of loop over events
  stopTimer(nevents);

  if (watch_!=0) watch_->Print();

  //now we need to store this the root output
  fout.cd();

  //typically i do this with a TH2(m0, m12), but i don't like this because i need to get the histogram grid correct
  //instead do:
  // TH1 storing m0
  // TH1 storing m12
  // TH1 storing <eff>
  /* probably don't need this anymore
  unsigned int nhistobins = bjetEffSum.size();
  TH1I * btageff_m0 = new TH1I("btageff_m0","m0/mGL coordinate for btag eff",nhistobins, 0,nhistobins);
  TH1I * btageff_m12 = new TH1I("btageff_m12","m12/mLSP coordinate for btag eff",nhistobins, 0,nhistobins);
  TH1D * btageff_avg = new TH1D("btageff_avg","average btag eff",nhistobins, 0,nhistobins);
  unsigned int ibin=1;
  for (map<pair<int,int> , double>::iterator ipoint = bjetEffSum.begin(); ipoint!=bjetEffSum.end(); ++ipoint) {
    double avgval=bjetSumSUSY[ipoint->first]>0 ? bjetEffSum[ipoint->first]/bjetSumSUSY[ipoint->first] : -1;
    btageff_m0->SetBinContent(ibin,ipoint->first.first);
    btageff_m12->SetBinContent(ibin,ipoint->first.second);
    btageff_avg->SetBinContent(ibin,avgval);

    ++ibin;
  }
  */


  fout.Write();
  cout<<"About to close output file"<<endl;
  fout.Close();  

  cout<<"Done closing output file. now terminate"<<endl;


}


unsigned int EventCalculator::getSeed(){

  //V00-02-35

  //if (sampleName_.Contains("ttZ") )                                                  return 4380;
  if (sampleName_.Contains("ttjets_madgraph") )                                      return 4381;
  if (sampleName_.Contains("TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12") )         return 4381;
  if (sampleName_.Contains("zjets") )                                                return 4382;
  if (sampleName_.Contains("ww") )                                                   return 4383;
  if (sampleName_.Contains("wz") )                                                   return 4384;
  if (sampleName_.Contains("zz") )                                                   return 4385;
  //if (sampleName_.Contains("t_s-channel") )                                          return 4386;
  //if (sampleName_.Contains("tbar_s-channel") )                                       return 4387;
  //if (sampleName_.Contains("t_t-channel") )                                          return 4388;
  //if (sampleName_.Contains("tbar_t-channel") )                                       return 4389;
  //if (sampleName_.Contains("t_tW-channel") )                                         return 4390;
  //if (sampleName_.Contains("tbar_tW-channel") )                                      return 4391;
  if (sampleName_.Contains("T_TuneZ2_s-channel_7TeV-powheg-tauola") )                return 4386;
  if (sampleName_.Contains("Tbar_TuneZ2_s-channel_7TeV-powheg-tauola") )             return 4387;
  if (sampleName_.Contains("T_TuneZ2_t-channel_7TeV-powheg-tauola") )                return 4388;
  if (sampleName_.Contains("Tbar_TuneZ2_t-channel_7TeV-powheg-tauola") )             return 4389;
  if (sampleName_.Contains("T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola") )            return 4390;
  if (sampleName_.Contains("Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola") )         return 4391;



  //if (sampleName_.Contains("wjets") )                                                return 4392;
  if (sampleName_.Contains("WJetsToLNu_TuneZ2_7TeV-madgraph-tauola") )               return 4392;
  if (sampleName_.Contains("WJetsToLNu_250_HT_300_TuneZ2_7TeV-madgraph-tauola") )    return 4401;
  if (sampleName_.Contains("WJetsToLNu_300_HT_inf_TuneZ2_7TeV-madgraph-tauola") )    return 4402;
  //if (sampleName_.Contains("DY") )                                                   return 4393;
  if (sampleName_.Contains("DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola") )          return 4393;
    
  if (sampleName_.Contains("QCD_Pt-0to5_TuneZ2_7TeV_pythia6") )                    return 4361;
  if (sampleName_.Contains("QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6") )                return 4362;
  if (sampleName_.Contains("QCD_Pt-120to170_TuneZ2_7TeV_pythia6") )                  return 4363;
  if (sampleName_.Contains("QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6") )                return 4364;
  if (sampleName_.Contains("QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6") )             return 4365;
  if (sampleName_.Contains("QCD_Pt-15to30_TuneZ2_7TeV_pythia6") )                    return 4366;
  if (sampleName_.Contains("QCD_Pt-170to300_TuneZ2_7TeV_pythia6") )                  return 4367;
  if (sampleName_.Contains("QCD_Pt-1800_TuneZ2_7TeV_pythia6") )                      return 4368;
  if (sampleName_.Contains("QCD_Pt-300to470_TuneZ2_7TeV_pythia6") )                  return 4369;
  if (sampleName_.Contains("QCD_Pt-30to50_TuneZ2_7TeV_pythia6") )                    return 4370;
  if (sampleName_.Contains("QCD_Pt-470to600_TuneZ2_7TeV_pythia6") )                  return 4371;
  if (sampleName_.Contains("QCD_Pt-50to80_TuneZ2_7TeV_pythia6") )                    return 4372;
  if (sampleName_.Contains("QCD_Pt-5to15_TuneZ2_7TeV_pythia6") )                     return 4373;
  if (sampleName_.Contains("QCD_Pt-600to800_TuneZ2_7TeV_pythia6") )                  return 4374;
  if (sampleName_.Contains("QCD_Pt-800to1000_TuneZ2_7TeV_pythia6") )                 return 4375;
  if (sampleName_.Contains("QCD_Pt-80to120_TuneZ2_7TeV_pythia6") )                 return 4376;

  if (sampleName_.Contains("LM9_SUSY_sftsht_7TeV-pythia6_v2") )                         return 4400;

  if (sampleName_.Contains("ht_run2011a_may10rereco") )                              return 4500;
  if (sampleName_.Contains("ht_run2011a_aug5rereco") )                               return 4501;
  if (sampleName_.Contains("ht_run2011a_promptrecov6") )                             return 4502;
  if (sampleName_.Contains("ht_run2011b_promptrecov1") )                             return 4503;
  if (sampleName_.Contains("ht_run2011b_promptrecov1_oct7") )                        return 4504;
  if (sampleName_.Contains("ht_run2011b_promptrecov1_oct14") )                       return 4505;
  if (sampleName_.Contains("ht_run2011a_promptrecov4_try3") )                        return 4506;

  if (sampleName_.Contains("T1bbbb") )                        return 4901;
  if (sampleName_.Contains("T1tttt") )                        return 4902;
  if (sampleName_.Contains("T2bb") )                        return 4903;
  if (sampleName_.Contains("T2tt") )                        return 4904;


  /*
  if (sampleName_.Contains("DYJetsToLL_TuneD6T_M-10To50_7TeV-madgraph-tauola") )     return 4357;
  if (sampleName_.Contains("DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola") )         return 4358;
  if (sampleName_.Contains("LM13_SUSY_sftsht_7TeV-pythia6") )                        return 4359;
  if (sampleName_.Contains("LM9_SUSY_sftsht_7TeV-pythia6_2") )                       return 4360;
  
  if (sampleName_.Contains("QCD_Pt_0to5_TuneZ2_7TeV_pythia6_3") )                    return 4361;
  if (sampleName_.Contains("QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6") )                return 4362;
  if (sampleName_.Contains("QCD_Pt_120to170_TuneZ2_7TeV_pythia6") )                  return 4363;
  if (sampleName_.Contains("QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6_3") )              return 4364;
  if (sampleName_.Contains("QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6") )             return 4365;
  if (sampleName_.Contains("QCD_Pt_15to30_TuneZ2_7TeV_pythia6") )                    return 4366;
  if (sampleName_.Contains("QCD_Pt_170to300_TuneZ2_7TeV_pythia6") )                  return 4367;
  if (sampleName_.Contains("QCD_Pt_1800_TuneZ2_7TeV_pythia6") )                      return 4368;
  if (sampleName_.Contains("QCD_Pt_300to470_TuneZ2_7TeV_pythia6") )                  return 4369;
  if (sampleName_.Contains("QCD_Pt_30to50_TuneZ2_7TeV_pythia6") )                    return 4370;
  if (sampleName_.Contains("QCD_Pt_470to600_TuneZ2_7TeV_pythia6") )                  return 4371;
  if (sampleName_.Contains("QCD_Pt_50to80_TuneZ2_7TeV_pythia6") )                    return 4372;
  if (sampleName_.Contains("QCD_Pt_5to15_TuneZ2_7TeV_pythia6") )                     return 4373;
  if (sampleName_.Contains("QCD_Pt_600to800_TuneZ2_7TeV_pythia6") )                  return 4374;
  if (sampleName_.Contains("QCD_Pt_800to1000_TuneZ2_7TeV_pythia6") )                 return 4375;
  if (sampleName_.Contains("QCD_Pt_80to120_TuneZ2_7TeV_pythia6_3") )                 return 4376;
  
  if (sampleName_.Contains("qcd_tunez2_pt0to5_summer11") )                           return 4361;
  if (sampleName_.Contains("qcd_tunez2_pt1000to1400_summer11") )                     return 4362;
  if (sampleName_.Contains("qcd_tunez2_pt120to170_summer11") )                       return 4363;
  if (sampleName_.Contains("qcd_tunez2_pt1400to1800_summer11") )                     return 4364;
  if (sampleName_.Contains("qcd_tunez2_pt15to3000_summer11") )                       return 4365;
  if (sampleName_.Contains("qcd_tunez2_pt15to30_summer11") )                         return 4366;
  if (sampleName_.Contains("qcd_tunez2_pt170to300_summer11") )                       return 4367;
  if (sampleName_.Contains("qcd_tunez2_pt1800_summer11") )                           return 4368;
  if (sampleName_.Contains("qcd_tunez2_pt300to470_summer11") )                       return 4369;
  if (sampleName_.Contains("qcd_tunez2_pt30to50_summer11") )                         return 4370;
  if (sampleName_.Contains("qcd_tunez2_pt470to600_summer11") )                       return 4371;
  if (sampleName_.Contains("qcd_tunez2_pt50to80_summer11") )                         return 4372;
  if (sampleName_.Contains("qcd_tunez2_pt5to15_summer11") )                          return 4373;
  if (sampleName_.Contains("qcd_tunez2_pt600to800_summer11") )                       return 4374;
  if (sampleName_.Contains("qcd_tunez2_pt800to1000_summer11") )                      return 4375;
  if (sampleName_.Contains("qcd_tunez2_pt80to120_summer11") )                        return 4376;
  
  if (sampleName_.Contains("TTJets_TuneD6T_7TeV-madgraph-tauola") )                  return 4377;
  if (sampleName_.Contains("TT_TuneZ2_7TeV-pythia6-tauola") )                  return 4412;
  if (sampleName_.Contains("TToBLNu_TuneZ2_s-channel_7TeV-madgraph") )               return 4378;
  if (sampleName_.Contains("TToBLNu_TuneZ2_t-channel_7TeV-madgraph") )               return 4379;
  if (sampleName_.Contains("TToBLNu_TuneZ2_tW-channel_7TeV-madgraph") )              return 4380;
  if (sampleName_.Contains("WJetsToLNu_TuneZ2_7TeV-madgraph-tauola") )               return 4381;
  if (sampleName_.Contains("WWtoAnything_TuneZ2_7TeV-pythia6-tauola") )              return 4382;
  if (sampleName_.Contains("WZtoAnything_TuneZ2_7TeV-pythia6-tauola") )              return 4383;
  if (sampleName_.Contains("ZinvisibleJets_7TeV-madgraph") )                         return 4384;
  if (sampleName_.Contains("ZZtoAnything_TuneZ2_7TeV-pythia6-tauola") )              return 4385;  
  
  if (sampleName_.Contains("ttjets_tunez2_madgraph_tauola_summer11") )               return 4400;
  
  if (sampleName_.Contains("ht_run2011a_423_jul1.txt") )                             return 4401;
  if (sampleName_.Contains("ht_run2011a_423_jul6.txt") )                             return 4402;
  if (sampleName_.Contains("ht_run2011a_423_jun10.txt") )                            return 4403;
  if (sampleName_.Contains("ht_run2011a_423_jun17.txt") )                            return 4404;
  if (sampleName_.Contains("ht_run2011a_423_jun24.txt") )                            return 4405;
  if (sampleName_.Contains("ht_run2011a_423_jun3.txt") )                             return 4406;
  if (sampleName_.Contains("ht_run2011a_423_may10rereco.txt") )                      return 4407;
  if (sampleName_.Contains("ht_run2011a_423_may10rereco_added.txt") )                return 4408;
  if (sampleName_.Contains("ht_run2011a_423_may10rereco_added_b.txt") )              return 4409;
  if (sampleName_.Contains("ht_run2011a_423_may24.txt") )                            return 4410;
  if (sampleName_.Contains("ht_run2011a_423_may27.txt") )                            return 4411;
  */
  
  cout << "Could not find seed for sample!" << endl;
  //  assert(0);
  return 123455;
}

double EventCalculator::calc_mNj( std::vector<unsigned int> jNi ) {

  //we could use an std::set which enforces that the elements are unique
  for (unsigned int i=0; i<jNi.size()-1; i++) {
    for (unsigned int j=i+1; j<jNi.size(); j++) {
      if (jNi.at(i) == jNi.at(j)) {
	cout<<"Problem in calc_mNj!"<<endl;
	return -1;
      }
    }
  }

  double sumE =0;
  double sumPx=0;
  double sumPy=0;
  double sumPz=0;

  for (unsigned int i=0; i<jNi.size(); i++) {
    unsigned int j1i = jNi.at(i);
    sumE += jets_AK5PF_energy->at( j1i ); 

    sumPx += getJetPx(j1i);
    sumPy += getJetPy(j1i);
    sumPz += getJetPz(j1i);
  }

  double sumP2 = sumPx*sumPx + sumPy*sumPy + sumPz*sumPz ;

  if ( (sumE*sumE) < sumP2 ) return -2. ;

  return sqrt( sumE*sumE - sumP2 );
}

double EventCalculator::calc_mNj( unsigned int j1i, unsigned int j2i) {
  std::vector<unsigned int> v;

  v.push_back(j1i);
  v.push_back(j2i);
  return calc_mNj(v);
}

double EventCalculator::calc_mNj( unsigned int j1i, unsigned int j2i, unsigned int j3i) {
  std::vector<unsigned int> v;

  v.push_back(j1i);
  v.push_back(j2i);
  v.push_back(j3i);
  return calc_mNj(v);
}



// void EventCalculator::changeVariables(TRandom3* random, double jetLossProbability, int& nLostJets)
// {
//   if(recalculatedVariables_) return;
//   recalculatedVariables_=true;
//   myJetsPF_temp                   = myJetsPF;		  
//   myMETPF_temp		          = myMETPF;		  

//   myJetsPF                = new std::vector<jet2_s>;	   //jmt -- switch to PF2PAT
//   myMETPF		  = new std::vector<met1_s>(*myMETPF_temp);	  

//   for(vector<jet2_s>::iterator thisJet = myJetsPF_temp->begin(); thisJet != myJetsPF_temp->end(); thisJet++)
//     {
//       if(random->Rndm() > jetLossProbability)
// 	{
// 	  myJetsPF->push_back(*thisJet);
// 	}
//       else
// 	{
// 	  double jetPx = thisJet->uncor_pt * cos(thisJet->uncor_phi);
// 	  double jetPy = thisJet->uncor_pt * sin(thisJet->uncor_phi);
// 	  double METx = myMETPF->at(0).pt * cos(myMETPF->at(0).phi) - jetPx;
// 	  double METy = myMETPF->at(0).pt * sin(myMETPF->at(0).phi) - jetPy;
// 	  myMETPF->at(0).pt = sqrt(METx*METx + METy*METy);
// 	  myMETPF->at(0).phi = atan2(METy,METx);
// 	}
//     }
//   nLostJets = myJetsPF_temp->size() - jets_AK5PF_pt->size();

// }

// void EventCalculator::resetVariables()
// {
//   if(!recalculatedVariables_) return;
//   recalculatedVariables_ = false;
//   if(myJetsPF != 0) delete myJetsPF;
//   else cout << "you've done something wrong!" << endl;
//   if(myMETPF != 0)delete myMETPF;		
//   else cout << "you've done something wrong!" << endl;

//   myJetsPF                = myJetsPF_temp;		  
//   myMETPF		  = myMETPF_temp;		  
    
// }


void EventCalculator::getSphericity(float & sph, bool addMET, bool addLeptons, const float jetthreshold) {

  float lambda1,lambda2;

  //misnomer in the sense that we're going to use all jets
  TMatrixD Smatrix(2,2);
  
  const  unsigned int njets=  jets_AK5PF_pt->size();
  for (unsigned int i=0; i<njets ; i++) {
    if ( isGoodJet(i, jetthreshold) ) {
      double px = getJetPx(i);
      double py = getJetPy(i);
      Smatrix[0][0] += px*px;
      Smatrix[0][1] += px*py;
      Smatrix[1][0] += px*py;
      Smatrix[1][1] += py*py;
    }
  }

  if (addMET) {
    double phi = getMETphi();
    double eT = getMET();
    double eX = eT*cos(phi);
    double eY = eT*sin(phi);
    Smatrix[0][0] += eX*eX;
    Smatrix[0][1] += eX*eY;
    Smatrix[1][0] += eX*eY;
    Smatrix[1][1] += eY*eY;
  }

  if (addLeptons) {
    for ( unsigned int i = 0; i< pf_mus_pt->size(); i++) {
      if (isCleanMuon(i,10)) {
	double px = pf_mus_pt->at(i) * cos(pf_mus_phi->at(i));
	double py = pf_mus_pt->at(i) * sin(pf_mus_phi->at(i));
	Smatrix[0][0] += px*px;
	Smatrix[0][1] += px*py;
	Smatrix[1][0] += px*py;
	Smatrix[1][1] += py*py;
      }
    }
    for (unsigned int i=0; i <pf_els_pt->size() ; i++) {
      if(isGoodElectron(i,false,10)) {
	double px = pf_els_pt->at(i) * cos(pf_els_phi->at(i));
	double py = pf_els_pt->at(i) * sin(pf_els_phi->at(i));
	Smatrix[0][0] += px*px;
	Smatrix[0][1] += px*py;
	Smatrix[1][0] += px*py;
	Smatrix[1][1] += py*py;
      }
    }
  }

  

  TMatrixDEigen eigenv(Smatrix);
  lambda1 = eigenv.GetEigenValuesRe()[0];
  lambda2 = eigenv.GetEigenValuesRe()[1];

  //for historical reasons, the output is returned via the reference rather than as a return value of the function
  sph = 2*lambda2 / (lambda1 + lambda2);
  
}



/* FIXME CFA
void EventCalculator::getTransverseThrustVariables(float & thrust, float & thrustPhi, bool addMET)
{
    //Begin Thrust Calculations
    //Copy Formula from http://home.thep.lu.se/~torbjorn/pythia/lutp0613man2.pdf 

    unsigned int njets = jets_AK5PF_pt->size();
     
    vector<double> goodJetPx;
    vector<double> goodJetPy;


    for (unsigned int i = 0; i<njets; i++) {
      if ( !isGoodJet(i) ) continue;
      //ngoodjets++;
      double phi = jets_AK5PF_phi->at(i);
      double pT = getJetPt(i);
      double pX = pT*cos(phi);
      double pY = pT*sin(phi);
      goodJetPx.push_back(pX);
      goodJetPy.push_back(pY);
    }

    if(addMET)
      {
	double metPhi = getMETphi();
	double eT = getMET();
	double eX = eT*cos(metPhi);
	double eY = eT*sin(metPhi);
	goodJetPx.push_back(eX);
	goodJetPy.push_back(eY);
      }

    double phiGuess = jets_AK5PF_phi->at(0);
    RooRealVar phi("phi","phi",phiGuess,-1.0*TMath::Pi(),TMath::Pi());
    RooTransverseThrustVar negThrustValue("negThrustValue","negThrustValue",phi,goodJetPx,goodJetPy);
    RooMinuit rm(negThrustValue);
    rm.migrad();
    thrust = -1*negThrustValue.getVal();
    thrustPhi = (phi.getVal()>0) ? phi.getVal()-TMath::Pi() : phi.getVal()+TMath::Pi();
}
*/

/* FIXME CFA
void EventCalculator::cutflow(itreestream& stream, int maxevents=-1){

  std::vector<int> npass;
  std::vector<double> sumw; //sum of weights
  std::vector<double> sumw2; //sum of weights squared
  int npass_eq1b=0;
  double sumw_eq1b=0,sumw2_eq1b=0; //special kludge for eq1b case

  int nevents = stream.size();

  setCutScheme();

  for (unsigned int i=0 ; i<cutTags_.size(); i++) {
    //std::cout << "cutTags_["<< i<< "]=" << cutTags_.at(i) << std::endl;
    npass.push_back(0);
    sumw.push_back(0);
    sumw2.push_back(0);
  }

  cout<<"Running..."<<endl;  



  startTimer();
  for(int entry=0; entry < nevents; ++entry){

    if(maxevents>0 && entry>=maxevents) break;

    // Read event into memory
    stream.read(entry);
    fillObjects();

    //some output to watch while it's running
    if(entry==0){
      if(!isSampleRealData()){
        cout << "MC xsec: " << getCrossSection() << endl;
      }
      else{
        cout << "This is data!"<< endl;
      }
      cout << "weight: "<< getWeight(nevents) << endl;
    }
    if(entry%100000==0) cout << "  entry: " << entry << ", percent done=" << (int)(entry/(double)nevents*100.)<<  endl;
    //if (entry%1000000==0) checkTimer(entry,nevents);
    
    const double weight = getWeight(nevents); //calculate weight
      
    for (unsigned int i=0 ; i<cutTags_.size(); i++) {

      if (cutRequired(cutTags_[i]) && passCut(cutTags_[i]) ) {
	npass.at(i) = npass.at(i) +1;
	sumw.at(i) = sumw.at(i) + weight;
	sumw2.at(i) = sumw2.at(i) + weight*weight;
      }
      else if (cutRequired(cutTags_[i]) && !passCut(cutTags_[i]) ) break;

      //special kludge for eq1b
      //if we reach this point then we have passed the cut. if the cut is >=1b, then we might also
      //pass the ==1b cut
      if ( cutRequired(cutTags_[i]) && cutTags_[i]=="cut1b") {
      	if (passCut("cutEq1b")){
      	  ++npass_eq1b;
      	  sumw_eq1b+=weight;
      	  sumw2_eq1b+= weight*weight;
      	}
      }

    }
  }

  cout<<endl;
  stopTimer(nevents);

  for (unsigned int i=0 ; i<npass.size(); i++) {
    
    if (cutRequired(cutTags_[i])) {
      
      double weighted = sumw.at(i);
      double weighted_error = sqrt(sumw2.at(i));
      
      char ccc[150];
      sprintf(ccc,"%20s %15d | %.2f | Weighted = %f +/- %f",cutNames_[cutTags_[i]].Data(),npass.at(i),100*weighted/sumw.at(0),weighted,weighted_error);
      cout<<ccc<<endl;

      //special kludge for ==1b case
      //this is an ugly duplication of code, to be dealt with later
      if ( cutTags_[i]=="cut1b") {
	//error = sqrt(npass_eq1b);
	weighted = sumw_eq1b;
	weighted_error = sqrt(sumw2_eq1b);
	sprintf(ccc,"%20s %15d | %.2f | Weighted = %f +/- %f","==1b",npass_eq1b,100*sumw_eq1b/sumw.at(0),weighted,weighted_error);
	cout<<ccc<<endl;
	//file <<"==1b"<<"\t"<<setprecision(20) << weighted<<"\t" << weighted_error<<endl;
	//fileU <<"==1b"<<"\t"<<setprecision(20) << npass_eq1b<<endl;

      }
    }
  }

  return;
}

*/


void EventCalculator::loadEventList(std::vector<int> &vrun, std::vector<int> &vlumi, std::vector<int> &vevent){

  std::ifstream inFile;
  //inFile.open("pfmht100_ormht150_eventlist_full.txt");        
  //inFile.open("premht25_pfmht100_ormht150_eventlist_full.txt");        
  //inFile.open("premht50_pfmht100_ormht150_eventlist_full.txt");        
  //inFile.open("premht100_pfmht100_ormht150_eventlist_full.txt");        
  //inFile.open("premht25_pfmht100_ormet150_eventlist_full.txt");        
  //inFile.open("premht25_pfmht100_ormet125_eventlist_full.txt");        
  //inFile.open("premht25_pfmht100_ormet100_eventlist_full.txt");        
  //inFile.open("eventlist_ht350_r178866_SUM.txt");        
  //inFile.open("eventlist_ht350l1fastjet_r178866_SUM.txt");        
  //inFile.open("eventlist_pfht350_r178866_SUM.txt"); 
  if (sampleName_.Contains("HT_Run2012A") )
    inFile.open("eventFilterLists/ecalLaserFilter_HT_Run2012A.txt");
  else  if (sampleName_.Contains("HTMHT_Run2012B") )
    inFile.open("eventFilterLists/ecalLaserFilter_HTMHT_Run2012B.txt");
  else  if (sampleName_.Contains("MET_Run2012A") )
    inFile.open("eventFilterLists/ecalLaserFilter_MET_Run2012A.txt");
  else  if (sampleName_.Contains("MET_Run2012B") )
    inFile.open("eventFilterLists/ecalLaserFilter_MET_Run2012B.txt");
  else  if (sampleName_.Contains("JetHT_Run2012B") )
    inFile.open("eventFilterLists/ecalLaserFilter_JetHT_Run2012B.txt");

  else return;

  if(!inFile || !inFile.good()) {std::cout << "ERROR: can't open event list" << std::endl;  assert(0);}
  while(!inFile.eof()) {
    std::string line;
    getline(inFile,line);
    std::string field;
    std::stringstream ss( line );
    //uint array[3];                                                          

    int value;
    uint i = 0;
    while (getline( ss, field, ':')){
      std::stringstream fs(field);
      value = 0;
      fs >> value;
      //fs >> array[i] ;                                                                       

      if(i==0)vrun.push_back(value);
      else if (i==1) vlumi.push_back(value);
      else vevent.push_back(value);
      i++;
    }
    //std::cout << array[0] << ":" << array[1] << ":" << array[2] << std::endl; 
    //if(event.id().run() == array[0] && event.id().luminosityBlock() == array[1] && event.id().event() == array[2]){ 
    //  eventPassTrigger = true; break;
    //}                                                                        
  }
  inFile.close();


}

bool EventCalculator::inEventList(std::vector<int> &vrun, std::vector<int> &vlumi, std::vector<int> &vevent){
  int lumi =  getLumiSection();
  int run = getRunNumber();
  int event = getEventNumber();
  
  bool onList = false;
  for(uint i = 0; i< vevent.size(); ++i){
    if( event == vevent.at(i)){
      if( lumi == vlumi.at(i)){
	if( run == vrun.at(i)){
	  onList = true;
	  break;
	}
      }
    }
  }
  return onList; 
}

//return the number of tops found (and the indices of the two top quarks)
int EventCalculator::findTop(int& top1, int& top2)
{
  bool foundOneTop = false;
  unsigned int maxParticles = mc_doc_id->size();
  for(unsigned int thisParticle = 0; thisParticle<maxParticles; thisParticle++)
    {
      if(abs(TMath::Nint(mc_doc_id->at(thisParticle))) == 6)
	{
	  if(!foundOneTop){
	    top1 = thisParticle;
	    foundOneTop=true;
	  }
	  else{
	    top2 = thisParticle;
	    return 2;
	  }
	}
    }
  if (foundOneTop) return 1;
  return 0;
}

//returns the decay mode of the W
//if W->enu, returns 11
//if W->munu, returns 13
//if W->taunu->enunu, returns 1511
//if W->taunu->mununu, returns 1513
//if W->taunu->hadrnic decay, return 15
//else if W->ud, returns 112
//else if W->cs, returns 134
//else returns 1 (should never be the case)
int EventCalculator::WDecayType(const int Wparent,int& Wdaughter)
{
  return 0; //FIXME CFA

/*
  int thisParticle = TMath::Nint(myGenParticles->at(Wparent).firstDaughter); 
  //int maxParticle = TMath::Nint(myGenParticles->at(Wparent).lastDaughter)+1;
  int maxParticle = TMath::Nint(myGenParticles->at(Wparent).firstDaughter)+2;//look at first two daughters only
  //for(int i = thisParticle; i < maxParticle ; i++){
  //  std::cout << "\t " << i << ", id = " << TMath::Nint(myGenParticles->at(i).pdgId) << std::endl;
  //}
  if(thisParticle==-1||maxParticle==0) return -1;
  for(; thisParticle < maxParticle ; thisParticle++)
    {
      if( abs(TMath::Nint(myGenParticles->at(thisParticle).pdgId)) == 11) {Wdaughter = thisParticle; return 11;}
      if( abs(TMath::Nint(myGenParticles->at(thisParticle).pdgId)) == 13) {Wdaughter = thisParticle; return 13;} 
      if( abs(TMath::Nint(myGenParticles->at(thisParticle).pdgId)) == 15)
	{
	  int tauDaughter = TMath::Nint(myGenParticles->at(thisParticle).firstDaughter);
	  int lastTauDaughter =  TMath::Nint(myGenParticles->at(thisParticle).lastDaughter);
	  if(tauDaughter == -1 || lastTauDaughter == -1) {Wdaughter = thisParticle; return -15;}
	  if(lastTauDaughter>=(int)myGenParticles->size()) {
	    cout << "tau daughters fall off the end of the generated particle list.  Tau daughter ranges from " 
		 << tauDaughter << " to " << lastTauDaughter << ". Size of the list is " << myGenParticles->size() 
		 << endl; 
	    Wdaughter = thisParticle; return -15;}
	  int currentTauDaughter = 0;
	  int WTauDaughterLastTauDaughterPlaceHolder = 0;
	  while(tauDaughter <= lastTauDaughter)
	    {
	      int tauDaugterPdgId = abs(TMath::Nint(myGenParticles->at(tauDaughter).pdgId));
	      if( tauDaugterPdgId != 11 && tauDaugterPdgId != 12 && tauDaugterPdgId != 13 && tauDaugterPdgId != 14 && tauDaugterPdgId != 15 && tauDaugterPdgId != 16 && tauDaugterPdgId != 22 && tauDaugterPdgId != 24 ){Wdaughter = tauDaughter; return 15;}
	      if(tauDaugterPdgId == 11) {Wdaughter = tauDaughter; return 1511;}
	      if(tauDaugterPdgId == 13) {Wdaughter = tauDaughter; return 1513;}
	      if(tauDaugterPdgId == 15) 
		{
		  lastTauDaughter =  TMath::Nint(myGenParticles->at(tauDaughter).lastDaughter);
		  tauDaughter = TMath::Nint(myGenParticles->at(tauDaughter).firstDaughter);
		  if(tauDaughter == -1 || lastTauDaughter == -1) {Wdaughter = tauDaughter; return -15;}
		  if(lastTauDaughter>(int)myGenParticles->size()) {
		    cout << "tau daughters fall off the end of the generated particle list.  Tau daughter ranges from " 
			 << tauDaughter << " to " << lastTauDaughter << ". Size of the list is " << myGenParticles->size() 
			 << endl; 
		    Wdaughter = thisParticle; return -15;}
		}
	      if(tauDaugterPdgId == 24) 
		{
		  currentTauDaughter =  tauDaughter;			    
		  WTauDaughterLastTauDaughterPlaceHolder = lastTauDaughter;
		  lastTauDaughter= TMath::Nint(myGenParticles->at(tauDaughter).lastDaughter);
		  tauDaughter = TMath::Nint(myGenParticles->at(tauDaughter).firstDaughter);
		}
	      else tauDaughter++;
	      if(tauDaughter > lastTauDaughter && WTauDaughterLastTauDaughterPlaceHolder!=0)
		{			    
		  tauDaughter = currentTauDaughter+1;
		  lastTauDaughter = WTauDaughterLastTauDaughterPlaceHolder;
		}
	    }
	  Wdaughter = thisParticle;
	  return -15;
	}
      if( abs(TMath::Nint(myGenParticles->at(thisParticle).pdgId)) == 1){Wdaughter = Wparent; return 112;}
      if( abs(TMath::Nint(myGenParticles->at(thisParticle).pdgId)) == 3){Wdaughter = Wparent; return 134;}
    }
  {Wdaughter = Wparent;return 1;}
*/
}

int EventCalculator::findW(int& W, int& Wdaughter,int parent=0, bool fromtop=true)
{
/*
   int thisParticle;
   int maxParticle;
   if(parent == 0)
     {
       thisParticle = 0;
       maxParticle = myGenParticles->size();
     }
   else
     {
       thisParticle = TMath::Nint(myGenParticles->at(parent).firstDaughter);
       maxParticle =  min(TMath::Nint(myGenParticles->at(parent).lastDaughter)+1, (int)myGenParticles->size()); 
     }
   for(;thisParticle < maxParticle;thisParticle++)
     {
       if(abs(TMath::Nint(myGenParticles->at(thisParticle).pdgId)) == 24)
	 {
	   //if fromtop is false, this W is required to NOT be from a top
	   //the flag is only relevant for the singletop-tW-channel sample
	   if( !fromtop && abs(TMath::Nint(myGenParticles->at( myGenParticles->at(thisParticle).firstMother ).pdgId))==6) continue;

	   int nextDaughter = TMath::Nint(myGenParticles->at(thisParticle).firstDaughter);
	   if(abs(TMath::Nint(myGenParticles->at(nextDaughter).pdgId)) == 24)
	     {
	       W = nextDaughter;
	       return WDecayType(W,Wdaughter);
	     }
	   else
	     {
	       W = thisParticle;
	       return WDecayType(W,Wdaughter);
	     }
	 }
     }
*/
   return -1;
}

/* FIXME CFA
int EventCalculator::muonMatch(const int trueMuon)
{
  //cout << "muon match" << endl;
  unsigned int maxMuons = myMuonsPF->size(); 
  for(unsigned int thisMuon = 0;thisMuon < maxMuons ; thisMuon++)
    {
      //cout << "deltaR " << deltaR(myMuonsPF->at(thisMuon).eta,myMuonsPF->at(thisMuon).phi,myGenParticles->at(trueMuon).eta,myGenParticles->at(trueMuon).phi) << endl;
      //if(isGoodMuon(thisMuon)) cout << "is good" << endl;
      //else cout << "is bad" << endl;
      if( isGoodMuon(thisMuon) && jmt::deltaR(myMuonsPF->at(thisMuon).eta,myMuonsPF->at(thisMuon).phi,myGenParticles->at(trueMuon).eta,myGenParticles->at(trueMuon).phi) < 0.3 )
	{ 
	  return thisMuon;
	}
    }
  for(unsigned int thisMuon = 0;thisMuon < maxMuons ; thisMuon++)
    {
      //check proximity to any reconstructed muon at deltaR<0.3
      if( jmt::deltaR(myMuonsPF->at(thisMuon).eta,myMuonsPF->at(thisMuon).phi,myGenParticles->at(trueMuon).eta,myGenParticles->at(trueMuon).phi) < 0.3 )
	{ 
	  return thisMuon;
	}
    }
  return -1;
}

int EventCalculator::electronMatch(const int trueElectron)
{
  unsigned int maxElectrons = myElectronsPF->size(); 
  for(unsigned int thisElectron = 0;thisElectron < maxElectrons ; thisElectron++)
    {
      //check proximity to any reconstructed muon at deltaR<0.3
      if( isGoodElectron(thisElectron) && jmt::deltaR(myElectronsPF->at(thisElectron).eta,myElectronsPF->at(thisElectron).phi,myGenParticles->at(trueElectron).eta,myGenParticles->at(trueElectron).phi) < 0.3 )
	{ 
	  return thisElectron;
	}
    }
  for(unsigned int thisElectron = 0;thisElectron < maxElectrons ; thisElectron++)
    {
      //check proximity to any reconstructed muon at deltaR<0.3
      if( jmt::deltaR(myElectronsPF->at(thisElectron).eta,myElectronsPF->at(thisElectron).phi,myGenParticles->at(trueElectron).eta,myGenParticles->at(trueElectron).phi) < 0.3 )
	{ 
	  return thisElectron;
	}
    }
  return -1;
}

int EventCalculator::tauMatch(const int trueTau)
{
  unsigned int maxJets = jets_AK5PF_pt->size(); 
  for(unsigned int thisJet = 0;thisJet < maxJets ; thisJet++)
    {
      //check proximity to any reconstructed muon at deltaR<0.4
      if( jmt::deltaR(myJetsPF->at(thisJet).eta,myJetsPF->at(thisJet).phi,myGenParticles->at(trueTau).eta,myGenParticles->at(trueTau).phi) < 0.4 )
	{ 
	  return thisJet;
	}
    }
  return -1;
}

//return the index of the reco object that the W daughter is matched to
int EventCalculator::daughterMatch(const int Wdaughter, const int WdecayType)
{
  if(WdecayType == 11 || WdecayType == 1511)
    {
      return electronMatch(Wdaughter);
    }
  if(WdecayType == 13 || WdecayType == 1513)
    {
      return muonMatch(Wdaughter);
    }
  if(WdecayType == 15)
    {
      return tauMatch(Wdaughter);
    }
  if(WdecayType == 1 || WdecayType == 112 || WdecayType == 134 )
    {
      return 0;
    }
  return -1;
}
*/

int EventCalculator::getTauDecayType( int tauid) {
  //  cout<<" -- tau decays -- "<<tauid<<endl;
  int founde=0;
  int foundmu=0;

  //kristen says:
  //loop over the "mc_electrons_*" and "mc_mus_*" to find those with "mother_id" equal to the tau, and "grandmother_id" equal to the W.

  //also, allow the gmother to be a tau as well:
  /*
    In the past, there has been occassional glitches with the mother and grandmother ids of the mc_mus and mc_electrons (where the mother id of the particle is entered in twice in the sequence), so to find if it came from a W->tau, I also allow the GRANDmother of the mc_mus or mc_electrons to be a tau
  */

  for ( size_t ii = 0; ii < mc_electrons_id->size(); ii++) {
    //    cout<<"\telec "<<mc_electrons_mother_id->at(ii)<<" "<<mc_electrons_grandmother_id->at(ii)<<" "<<mc_electrons_ggrandmother_id->at(ii)<<endl;
    //this should check the charge at well
    if ((mc_electrons_mother_id->at(ii) == tauid) && (abs(mc_electrons_grandmother_id->at(ii))==24  ) ) founde++;
    else if ( (mc_electrons_mother_id->at(ii) == tauid) && (mc_electrons_grandmother_id->at(ii) == tauid) && (abs(mc_electrons_ggrandmother_id->at(ii)) == 24)) founde++;

  }

  for ( size_t ii = 0; ii < mc_mus_id->size(); ii++) {
    //    cout<<"\tmu "<<mc_mus_mother_id->at(ii)<<" "<<mc_mus_grandmother_id->at(ii)<<" "<<mc_mus_ggrandmother_id->at(ii)<<endl;
    //this should check the charge at well
    if ((mc_mus_mother_id->at(ii) == tauid) && (abs(mc_mus_grandmother_id->at(ii))==24  )) foundmu++;
    else if ( (mc_mus_mother_id->at(ii) == tauid) && (mc_mus_grandmother_id->at(ii) == tauid) && (abs(mc_mus_ggrandmother_id->at(ii)) == 24)) foundmu++;
  }

  //  cout<<" tau decays -- "<<founde<<" "<<foundmu<<endl;


//   cout<<" ~~~~~~~~~ now look at MC taus"<<endl;
//   for (size_t ii = 0; ii<mc_taus_id->size(); ii++) {
//     cout<<ii<<"\t"<<mc_taus_id->at(ii)<<" "<<mc_taus_grandmother_id->at(ii)<<" "<<mc_taus_ggrandmother_id->at(ii)<<endl;
//   }


//seems that there's no way to dig in and actually look for tau->had. just have to trust that all other decays are tau->had

  if (founde>0) return 1;
  if (foundmu>0) return 2;
  return 0;
}

 //trying to rewrite from scratch for cfA
int EventCalculator::getTTbarDecayType() {
   
  //use my old encoding scheme instead of don's

  int founde=0;
  int foundmu=0;
  int foundhad=0;
  int foundtau=0;

  int taudecay[2]={-1,-1};

  //  cout<< "----"<<endl;
  for (size_t jj=0; jj<  mc_doc_id->size(); ++jj) {
    //    cout<<jj <<"\t"<<mc_doc_id->at(jj)<<" "<<mc_doc_mother_id->at(jj)<<" "<<mc_doc_grandmother_id->at(jj)<<endl;

    //go looking for particle who have the top as ~grandmother~
    if ( abs(TMath::Nint(mc_doc_grandmother_id->at(jj))) == 6 ) {

      //these should be the W decay products: q,l,nu
      int wdau=    abs(TMath::Nint(mc_doc_id->at(jj)));
      if (wdau == 11) { //e
	founde++;
      }
      else if (wdau == 13) { // mu
	foundmu++;
      }
      else if (wdau == 15) { //tau
	taudecay[foundtau]=getTauDecayType(TMath::Nint(mc_doc_id->at(jj)) );
	foundtau++;
      }
      else if (wdau == 12 || wdau ==14 || wdau==16) {
	//neutrinos
      }
      else if (wdau>=1 && wdau<=4) { //quarks
	foundhad++;
      }
      else {
	cout<<" Wdau = "<<wdau<<endl;
      }

    }
  }

  //  enum TopDecayCategory {kTTbarUnknown=0,kAllLeptons=1,kAllHadronic=2,kOneElectron=3,kOneMuon=4,kOneTauE=5,kOneTauMu=6,kOneTauHadronic=7,kAllTau=8,kTauPlusLepton=9, nTopCategories=10};

  if (!(foundhad%2==0)) cout<<"PROBLEM: foundhad = "<<foundhad<<endl; //W->qq' should always generate an even number!

  foundhad /=2;

  if ( founde+foundmu == 2) return 1;
  else if ( foundhad == 2) return 2;
  else if ( founde==1 && foundhad==1 ) return 3;
  else if ( foundmu==1 && foundhad==1 ) return 4;
  else if (foundtau==1 && foundhad==1 && taudecay[0]==1) return 5;
  else if (foundtau==1 && foundhad==1 && taudecay[0]==2) return 6;
  else if (foundtau==1 && foundhad==1 && taudecay[0]==0) return 7;
  else if ( foundtau==2) return 8;
  else if (foundtau==1 && founde==1) return 9;
  else if (foundtau==1 && foundmu==1) return 9;

  else {
    cout<<"PROBLEM -- "<<founde<<" "<<foundmu<<" "<<foundtau<<" "<<foundhad<<" tau info -- "<<taudecay[0]<<" "<<taudecay[1]<<endl;
  }

  return 0;
}
 
/* old version pre cfA
int EventCalculator::getTTbarDecayType(int& W1decayType, int& W2decayType, int& W1, int& W1daughter, int& W2, int& W2daughter, bool passW2info=false)
{
  if(myGenParticles == 0 || myGenParticles->size() == 0) return -1;
  int top1=-1;
  int top2=-1;
  int nTops = findTop(top1,top2);
  W1 = 0;
  W1daughter = 0;
  W1decayType = -1;
  if(!passW2info){//for single-top-tW-channel
    W2 = 0;
    W2daughter = 0;
    W2decayType = -1;
  }

  if(nTops>0)
    {
      W1decayType = findW(W1,W1daughter,top1);
      if(nTops>1) W2decayType = findW(W2,W2daughter,top2);
    }

  //std::cout << "W1decayType = " << W1decayType << ", W2decayType = " << W2decayType << std::endl;
  //std::cout << "W1daughter = " << W1daughter << ", W2daughter = " << W2daughter << std::endl;
  
  int W1daughterMatch = -1;
  int W2daughterMatch = -1;
  if( W1decayType > 0 ) W1daughterMatch = daughterMatch(W1daughter,W1decayType);
  if( W2decayType > 0 ) W2daughterMatch = daughterMatch(W2daughter,W2decayType);

  //cout << "done with matching" << endl;

  //for ttbar, considering only (W->had,W->e/mu), (W->tau->had,W->e/mu), (W->had,W->had), (W->had,W->tau->had) 
  if( ((sampleName_.Contains("TTJets_TuneZ2") || sampleName_.Contains("ttjets_madgraph")
	|| sampleName_.Contains("T_TuneZ2_tW") || sampleName_.Contains("Tbar_TuneZ2_tW"))
       && (W1decayType == 112 || W2decayType == 112 || W1decayType == 134 || W2decayType == 134 || W1decayType == 15 || W2decayType == 15) 
       && !(W1decayType == 15 && W2decayType == 15))
      || (sampleName_.Contains("T_TuneZ2_s-channel") || sampleName_.Contains("Tbar_TuneZ2_s-channel")
	  || sampleName_.Contains("T_TuneZ2_t-channel") || sampleName_.Contains("Tbar_TuneZ2_t-channel"))
      )
  //considering only semileptonic decay: (W->had,W->e/mu)
  //if( ((W1decayType==13 && W2decayType==1) || (W1decayType==1 && W2decayType==13) || (W1decayType==1513 && W2decayType==1) || (W1decayType==1 && W2decayType==1513))
  //    || ((W1decayType==11 && W2decayType==1) || (W1decayType==1 && W2decayType==11) || (W1decayType==1511 && W2decayType==1) || (W1decayType==1 && W2decayType==1511)))
    {
      //find the e/mu decay one
      int WdecayType, Wdaughter,WdaughterMatch;
      if(W1decayType!=112 && W1decayType!=134) //if W1 is the non-hadronic one
	{
	  WdecayType =  W1decayType;
	  Wdaughter = W1daughter;
	  WdaughterMatch = W1daughterMatch;
	}
      else if(W2decayType!=112 && W2decayType!=134)
	{
	  WdecayType =  W2decayType;
	  Wdaughter = W2daughter;
	  WdaughterMatch = W2daughterMatch;
	}
      else if(W1decayType!=15)//otherwise, it's a (W->tau->had,W->e/mu) one
	{
	  WdecayType =  W1decayType;
	  Wdaughter = W1daughter;
	  WdaughterMatch = W1daughterMatch;
	}
      else
	{
	  WdecayType =  W2decayType;
	  Wdaughter = W2daughter;
	  WdaughterMatch = W2daughterMatch;
	}
      if(WdecayType == 11 || WdecayType == 1511) //electron
	{
	  //cout << "electron match" << endl;
	  if( WdaughterMatch > -1 ) // electron on electron list
	    {
	      if( isGoodElectron(WdaughterMatch) ) 
		return 101100;//there is match, and it passes all cuts
	      else if( 
		      fabs(myElectronsPF->at(WdaughterMatch).superCluster_eta) > 2.5
		      || fabs(myGenParticles->at(Wdaughter).eta) > 2.5 ) 
		//there is a match, but it is out of eta range
		return 101101;
	      else if( myElectronsPF->at(WdaughterMatch).pt < 10 
		       || myGenParticles->at(Wdaughter).pt < 10 ) 
		//there is a match, it is in eta range, but it is out of pt range
		return 101102;
	      else if ((myElectronsPF->at(WdaughterMatch).chargedHadronIso
			+ myElectronsPF->at(WdaughterMatch).photonIso
			+ myElectronsPF->at(WdaughterMatch).neutralHadronIso)/myElectronsPF->at(WdaughterMatch).pt >=0.2 )
		//there is a match, it is in eta and pt range, but it fails isolation cut
		return 101103;
	      else return 101104;//it fails some other quality cut
	    }
	  else
	    {
	      if(fabs(myGenParticles->at(Wdaughter).eta) > 2.5) return 201101;//
	      else if(myGenParticles->at(Wdaughter).pt < 10) return 201102;//
	      else return 201103;
	    }
	}
      if(WdecayType == 13 || WdecayType == 1513) //muon
	{
	  //cout << "muon match" << endl;
	  if( WdaughterMatch > -1 ) // muon on muon list
	    {
	      //cout << "matched daughter" << endl;
	      //cout << "true index " << Wdaughter << endl;
	      //cout << "match index " << WdaughterMatch << endl;
	      if( isGoodMuon(WdaughterMatch) ) return 101300;
	      else if( fabs(myMuonsPF->at(WdaughterMatch).eta) > 2.4 
		       || fabs(myGenParticles->at(Wdaughter).eta) > 2.4 ) 
		return 101301;
	      else if( myMuonsPF->at(WdaughterMatch).pt < 10 
		       || myGenParticles->at(Wdaughter).pt < 10 ) 
		return 101302;
	      else if ( (myMuonsPF->at(WdaughterMatch).chargedHadronIso
			 + myMuonsPF->at(WdaughterMatch).photonIso
			 + myMuonsPF->at(WdaughterMatch).neutralHadronIso)/myMuonsPF->at(WdaughterMatch).pt >=0.2 )
		return 101303;
	      else return 101304;
	    }
	  else
	    {
	      //cout << "unmatched daughter" << endl;
	      //cout << "true index " << Wdaughter << endl;
	      if(fabs(myGenParticles->at(Wdaughter).eta) > 2.4) return 201301;//
	      else if(myGenParticles->at(Wdaughter).pt < 10) return 201302;//
	      else return 201303;
	    }
	}
      if(WdecayType == 15)//tau
	{
	  //cout << "tau match" << endl;
	  if( WdaughterMatch > -1 ) // tau on jet list
	    {
	      if( isGoodJet(WdaughterMatch) ) return 101500;
	      else return 101501;
	    }
	  else
	    {
	      return 201500;
	    }
	}
    }
  
  return 0;
}


int EventCalculator::getWDecayType(int& WdecayType, int& W, int& Wdaughter, bool fromtop=true)
{
  if(myGenParticles == 0 || myGenParticles->size() == 0) return -1;

  WdecayType = findW(W,Wdaughter, fromtop);
  
  int WdaughterMatch = -1;

  if( WdecayType > 0 ) WdaughterMatch = daughterMatch(Wdaughter,WdecayType);

  //cout << "done with matching" << endl;

  if(WdecayType == 11 || WdecayType == 1511) //electron
    {
      //cout << "electron match" << endl;
      if( WdaughterMatch > -1 ) // electron on electron list
	{
	  if( isGoodElectron(WdaughterMatch) ) 
	    return 101100;//there is match, and it passes all cuts
	  else if( 
		  fabs(myElectronsPF->at(WdaughterMatch).superCluster_eta) > 2.5
		  || fabs(myGenParticles->at(Wdaughter).eta) > 2.5 ) 
	    //there is a match, but it is out of eta range
	    return 101101;
	  else if( myElectronsPF->at(WdaughterMatch).pt < 10 
		   || myGenParticles->at(Wdaughter).pt < 10 ) 
	    //there is a match, it is in eta range, but it is out of pt range
	    return 101102;
	  else if ((myElectronsPF->at(WdaughterMatch).chargedHadronIso
		    + myElectronsPF->at(WdaughterMatch).photonIso
		    + myElectronsPF->at(WdaughterMatch).neutralHadronIso)/myElectronsPF->at(WdaughterMatch).pt >=0.2 )
	    //there is a match, it is in eta and pt range, but it fails isolation cut
	    return 101103;
	  else return 101104;//it fails some other quality cut
	}
      else
	{
	  if(fabs(myGenParticles->at(Wdaughter).eta) > 2.5) return 201101;//
	  else if(myGenParticles->at(Wdaughter).pt < 10) return 201102;//
	  else return 201103;
	}
    }
  if(WdecayType == 13 || WdecayType == 1513) //muon
    {
      //cout << "muon match" << endl;
      if( WdaughterMatch > -1 ) // muon on muon list
	{
	  //cout << "matched daughter" << endl;
	  //cout << "true index " << Wdaughter << endl;
	  //cout << "match index " << WdaughterMatch << endl;
	  if( isGoodMuon(WdaughterMatch) ) return 101300;
	  else if( fabs(myMuonsPF->at(WdaughterMatch).eta) > 2.4 
		   || fabs(myGenParticles->at(Wdaughter).eta) > 2.4 ) 
	    return 101301;
	  else if( myMuonsPF->at(WdaughterMatch).pt < 10 
		   || myGenParticles->at(Wdaughter).pt < 10 ) 
	    return 101302;
	  else if ( (myMuonsPF->at(WdaughterMatch).chargedHadronIso
		     + myMuonsPF->at(WdaughterMatch).photonIso
		     + myMuonsPF->at(WdaughterMatch).neutralHadronIso)/myMuonsPF->at(WdaughterMatch).pt >=0.2 )
	    return 101303;
	  else return 101304;
	}
      else
	{
	  //cout << "unmatched daughter" << endl;
	  //cout << "true index " << Wdaughter << endl;
	  if(fabs(myGenParticles->at(Wdaughter).eta) > 2.4) return 201301;//
	  else if(myGenParticles->at(Wdaughter).pt < 10) return 201302;//
	  else return 201303;
	}
    }
  if(WdecayType == 15)//tau
    {
      //cout << "tau match" << endl;
      if( WdaughterMatch > -1 ) // tau on jet list
	{
	  if( isGoodJet(WdaughterMatch) ) return 101500;
	  else return 101501;
	}
      else
	{
	  return 201500;
	}
    }
  
  
  return 0;
}
*/

void EventCalculator::sampleAnalyzer() {


  TFile fout("histos.root","RECREATE");
  
  //TH1F * h_MCPU = new TH1F("h_MCPU","unweighted PU distribution",35,0.5,34.5);
  //TH1F * h_MCPUr = new TH1F("h_MCPUr","reweighted PU distribution",35,-0.5,34.5);
  //TH1F * h_MCPUBXp1 = new TH1F("h_MCPUBXp1","unweighted PU distribution",35,0.5,34.5);
  //TH1F * h_MCPUrBXp1 = new TH1F("h_MCPUrBxp1","reweighted PU distribution",35,-0.5,34.5);
  //TH1F * h_MCPUBXm1 = new TH1F("h_MCPUBXm1","unweighted PU distribution",35,0.5,34.5);
  //TH1F * h_MCPUrBXm1 = new TH1F("h_MCPUrBXm1","reweighted PU distribution",35,-0.5,34.5);

  //TH1F * h_ht = new TH1F("h_ht","HT",1000,0,1000);
  //TH1F * h_ht_trig = new TH1F("h_ht_trig","HT",1000,0,1000);

  //TH1F * h_met = new TH1F("h_met","HT",1000,0,1000);
  //TH1F * h_met_trig = new TH1F("h_met_trig","HT",1000,0,1000);

  //TH1F * h_met_den_mht75 = new TH1F("h_met_den_mht75","HT",1000,0,1000);
  //TH1F * h_met_trig_mht75 = new TH1F("h_met_trig_mht75","HT",1000,0,1000);
  //TH1F * h_met_den_mht80 = new TH1F("h_met_den_mht80","HT",1000,0,1000);
  //TH1F * h_met_trig_mht80 = new TH1F("h_met_trig_mht80","HT",1000,0,1000);
  //TH1F * h_met_den_mht90_ht300 = new TH1F("h_met_den_mht90_ht300","HT",1000,0,1000);
  //TH1F * h_met_trig_mht90_ht300 = new TH1F("h_met_trig_mht90_ht300","HT",1000,0,1000);
  //TH1F * h_met_den_mht90_ht350 = new TH1F("h_met_den_mht90_ht350","HT",1000,0,1000);
  //TH1F * h_met_trig_mht90_ht350 = new TH1F("h_met_trig_mht90_ht350","HT",1000,0,1000);
  //TH1F * h_met_den_mht110 = new TH1F("h_met_den_mht110","HT",1000,0,1000);
  //TH1F * h_met_trig_mht110 = new TH1F("h_met_trig_mht110","HT",1000,0,1000);
  //
  //TH1F * h_met_0L_den_mht75       = new TH1F("h_met_0L_den_mht75","HT",1000,0,1000);
  //TH1F * h_met_0L_trig_mht75      = new TH1F("h_met_0L_trig_mht75","HT",1000,0,1000);
  //TH1F * h_met_0L_den_mht80       = new TH1F("h_met_0L_den_mht80","HT",1000,0,1000);
  //TH1F * h_met_0L_trig_mht80      = new TH1F("h_met_0L_trig_mht80","HT",1000,0,1000);
  //TH1F * h_met_0L_den_mht90_ht300 = new TH1F("h_met_0L_den_mht90_ht300","HT",1000,0,1000);
  //TH1F * h_met_0L_trig_mht90_ht300= new TH1F("h_met_0L_trig_mht90_ht300","HT",1000,0,1000);
  //TH1F * h_met_0L_den_mht90_ht350 = new TH1F("h_met_0L_den_mht90_ht350","HT",1000,0,1000);
  //TH1F * h_met_0L_trig_mht90_ht350= new TH1F("h_met_0L_trig_mht90_ht350","HT",1000,0,1000);
  //TH1F * h_met_0L_den_mht110      = new TH1F("h_met_0L_den_mht110","HT",1000,0,1000);
  //TH1F * h_met_0L_trig_mht110     = new TH1F("h_met_0L_trig_mht110","HT",1000,0,1000);
  //
  //TH1F * h_met_1L_den_mht75       = new TH1F("h_met_1L_den_mht75","HT",1000,0,1000);
  //TH1F * h_met_1L_trig_mht75      = new TH1F("h_met_1L_trig_mht75","HT",1000,0,1000);
  //TH1F * h_met_1L_den_mht80       = new TH1F("h_met_1L_den_mht80","HT",1000,0,1000);
  //TH1F * h_met_1L_trig_mht80      = new TH1F("h_met_1L_trig_mht80","HT",1000,0,1000);
  //TH1F * h_met_1L_den_mht90_ht300 = new TH1F("h_met_1L_den_mht90_ht300","HT",1000,0,1000);
  //TH1F * h_met_1L_trig_mht90_ht300= new TH1F("h_met_1L_trig_mht90_ht300","HT",1000,0,1000);
  //TH1F * h_met_1L_den_mht90_ht350 = new TH1F("h_met_1L_den_mht90_ht350","HT",1000,0,1000);
  //TH1F * h_met_1L_trig_mht90_ht350= new TH1F("h_met_1L_trig_mht90_ht350","HT",1000,0,1000);
  //TH1F * h_met_1L_den_mht110      = new TH1F("h_met_1L_den_mht110","HT",1000,0,1000);
  //TH1F * h_met_1L_trig_mht110     = new TH1F("h_met_1L_trig_mht110","HT",1000,0,1000);
  //
  //TH1F * h_met_1e0mu_den_mht75       = new TH1F("h_met_1e0mu_den_mht75","HT",1000,0,1000);
  //TH1F * h_met_1e0mu_trig_mht75      = new TH1F("h_met_1e0mu_trig_mht75","HT",1000,0,1000);
  //TH1F * h_met_1e0mu_den_mht80       = new TH1F("h_met_1e0mu_den_mht80","HT",1000,0,1000);
  //TH1F * h_met_1e0mu_trig_mht80      = new TH1F("h_met_1e0mu_trig_mht80","HT",1000,0,1000);
  //TH1F * h_met_1e0mu_den_mht90_ht300 = new TH1F("h_met_1e0mu_den_mht90_ht300","HT",1000,0,1000);
  //TH1F * h_met_1e0mu_trig_mht90_ht300= new TH1F("h_met_1e0mu_trig_mht90_ht300","HT",1000,0,1000);
  //TH1F * h_met_1e0mu_den_mht90_ht350 = new TH1F("h_met_1e0mu_den_mht90_ht350","HT",1000,0,1000);
  //TH1F * h_met_1e0mu_trig_mht90_ht350= new TH1F("h_met_1e0mu_trig_mht90_ht350","HT",1000,0,1000);
  //TH1F * h_met_1e0mu_den_mht110      = new TH1F("h_met_1e0mu_den_mht110","HT",1000,0,1000);
  //TH1F * h_met_1e0mu_trig_mht110     = new TH1F("h_met_1e0mu_trig_mht110","HT",1000,0,1000);
  //
  //TH1F * h_met_1mu0e_den_mht75       = new TH1F("h_met_1mu0e_den_mht75","HT",1000,0,1000);
  //TH1F * h_met_1mu0e_trig_mht75      = new TH1F("h_met_1mu0e_trig_mht75","HT",1000,0,1000);
  //TH1F * h_met_1mu0e_den_mht80       = new TH1F("h_met_1mu0e_den_mht80","HT",1000,0,1000);
  //TH1F * h_met_1mu0e_trig_mht80      = new TH1F("h_met_1mu0e_trig_mht80","HT",1000,0,1000);
  //TH1F * h_met_1mu0e_den_mht90_ht300 = new TH1F("h_met_1mu0e_den_mht90_ht300","HT",1000,0,1000);
  //TH1F * h_met_1mu0e_trig_mht90_ht300= new TH1F("h_met_1mu0e_trig_mht90_ht300","HT",1000,0,1000);
  //TH1F * h_met_1mu0e_den_mht90_ht350 = new TH1F("h_met_1mu0e_den_mht90_ht350","HT",1000,0,1000);
  //TH1F * h_met_1mu0e_trig_mht90_ht350= new TH1F("h_met_1mu0e_trig_mht90_ht350","HT",1000,0,1000);
  //TH1F * h_met_1mu0e_den_mht110      = new TH1F("h_met_1mu0e_den_mht110","HT",1000,0,1000);
  //TH1F * h_met_1mu0e_trig_mht110     = new TH1F("h_met_1mu0e_trig_mht110","HT",1000,0,1000);


  //TH2F * h_Wlnueta = new TH2F("h_Wlnueta","lepton eta vs neutrino eta",100,-5,5,100,-5,5);
  //TH2F * h_Wlnupt = new TH2F("h_Wlnupt","lepton pt vs neutrino pt",100,0,300,100,0,300);
  //
  //TH2F * h_Wtaulnueta = new TH2F("h_Wtaulnueta","tautolepton eta vs neutrino eta",100,-5,5,100,-5,5);
  //TH2F * h_Wtaulnupt = new TH2F("h_Wtaulnupt","tautolepton pt vs neutrino pt",100,0,300,100,0,300);

  //std::vector<int> vrun,vlumi,vevent;
  //loadEventList(vrun, vlumi, vevent);

  const Long64_t nevents = chainA->GetEntries();
  const Long64_t neventsB = chainB->GetEntries();
  assert(nevents==neventsB);

  //float prob0,probge1,prob1,probge2,probge3;

  //float n0b = 0, nge1b = 0, neq1b = 0, nge2b = 0, nge3b = 0;

  ////initialize PU things
  //std::vector< float > DataDist2011;
  //std::vector< float > MCDist2011;
  //for( int i=0; i<35; ++i) {
  //  //for PU_S3 (only in-time)
  //  //DataDist2011.push_back(pu::ObsDist2011_f[i]);
  //  //MCDist2011.push_back(pu::PoissonOneXDist_f[i]);
  //  //for 3dPU reweighting 
  //  DataDist2011.push_back(pu::TrueDist2011_f[i]);
  //  if(i<25)
  //    MCDist2011.push_back(pu::probdistFlat10_f[i]);
  //  //std::cout << i << " " << pu::probdistFlat10_f[i] << std::endl;
  //}  
  ////for( int i=0; i<1000; ++i) {
  ////  DataDist2011.push_back(pu::TrueDist2011_finebin_f[i]);
  ////}
  ////reweight::LumiReWeighting LumiWeights = reweight::LumiReWeighting( MCDist2011, DataDist2011 );
  ////Lumi3DReWeighting LumiWeights = Lumi3DReWeighting("Summer11MC_PUDistTruth.root", "PileupTruth_finebin_data_SUM_upto180252.root", "mcpileup", "pileup");
  //Lumi3DReWeighting LumiWeights = Lumi3DReWeighting( MCDist2011, DataDist2011);
  //
  //std::cout << "size of MC = " << MCDist2011.size() << ", size of Data = " << DataDist2011.size() << std::endl;
  //LumiWeights.weight3D_init(1);
  ////LumiWeights.weight3D_init();
  ////LumiWeights.weight3D_init("Weight3D.root");

  cout<<"Running..."<<endl;  
  int npass = 0;

  //int decayType;
  int W1decayType;
  //int W2decayType;
  int nSemiMu=0,nSemiEle=0,nSemiTauHad=0,nDilep=0,nHad=0,nOther=0;

  startTimer();
  for(Long64_t entry=0; entry < nevents; ++entry){
    chainB->GetEntry(entry);
    chainA->GetEntry(entry);


    if(entry%10000==0) cout << "entry: " << entry << ", percent done=" << (int)(entry/(double)nevents*100.)<<  endl;
    //int event = getEventNumber();
    //if(event !=86684 ) continue;
    
    //std::cout << " HT = " << getHT() << std::endl;
    //
    //for (unsigned int i = 0; i < jets_AK5PF_pt->size(); ++i) {
    //  std::cout << "jet " << i << " pt = " << getJetPt(i)<< ", eta = " << myJetsPF->at(i).eta << ", passID = " << jetPassLooseID(i) <<std::endl;
    //
    //  std::cout << "\tneutralHadEnFrac (should be <0.99) = " << myJetsPF->at(i).neutralHadronEnergyFraction << std::endl;
    //  std::cout << "\tneutralEmEnFrac  (should be <0.99) = " << myJetsPF->at(i).neutralEmEnergyFraction << std::endl;
    //  std::cout << "\tnumDaughters (should be >1)        = " << myJetsPF->at(i).numberOfDaughters << std::endl;
    //  std::cout << "\tchargedHadEnFrac (should be >0)    = " << myJetsPF->at(i).chargedHadronEnergyFraction << std::endl;
    //  std::cout << "\tchargedEmEnFrac  (should be <0.99) = " << myJetsPF->at(i).chargedEmEnergyFraction << std::endl;
    //  std::cout << "\tchargedMult (should be >0)         = " << myJetsPF->at(i).chargedMultiplicity << std::endl;
    // 
    //
    //
    //}
    //for(uint i = 0; i< myGenParticles->size(); ++i){
    //  if(myGenParticles->at(i).status == 3){
    //  std::cout << i << ": ID = " << myGenParticles->at(i).pdgId 
    //		<< ", status = " << myGenParticles->at(i).status 
    //		<< ", pt = " << myGenParticles->at(i).pt 
    //		<< ", eta = " << myGenParticles->at(i).eta 
    //		<< ", phi = " << myGenParticles->at(i).phi 
    //		<< ", mmin = " << myGenParticles->at(i).firstMother 
    //		<< ", mmax = " << myGenParticles->at(i).lastMother
    //		<< ", dmin = " << myGenParticles->at(i).firstDaughter 
    //		<< ", dmax = " << myGenParticles->at(i).lastDaughter
    //		<< std::endl;
    //  }
    //}

    //for (unsigned int k = 0; k<myGenParticles->size(); k++) {
    //  //      for debugging
    //  if(myGenParticles->at(k).status==3)
    //	cout<<k<<"\t"<<TMath::Nint(myGenParticles->at(k).pdgId )<<"\t"<< TMath::Nint(myGenParticles->at(k).firstMother)<<"\t"<<TMath::Nint(myGenParticles->at(k).lastMother)<<"\t"<<TMath::Nint(myGenParticles->at(k).status ) <<endl;
    //}

    int W1, W1daughter, returnvalue;
    //int W2, W2daughter;

    returnvalue = 1; //getWDecayType(W1decayType, W1, W1daughter, false);  //FIXME CFA
    //decayType = getTTbarDecayType(W1decayType, W2decayType, W1, W1daughter, W2, W2daughter,true);


    //std::cout << "decaytype = " << decayType << ", W1decayType = " << W1decayType << ", W2decayType = " << W2decayType 
    //  		<< ", W1 = " << W1 << ", W1 daughter = " << W1daughter << ", W2 = " << W2 << ", W2daughter = " << W2daughter << std::endl;

    //
    ////W1decayType = findW(W1,W1daughter);
    //if(returnvalue==101302){
      std::cout << "W1decayType = " << W1decayType << ", W1 = " << W1 << ", W1daughter = " << W1daughter << std::endl;
      std::cout << "\t returnvalue = " << returnvalue << std::endl;
    //}
    //std::cout << "myGenparticles size = " << myGenParticles->size() << std::endl;
    /*
    decayType = getTTbarDecayType(W1decayType, W2decayType, W1, W1daughter, W2, W2daughter);

    if(W1decayType == 11 || W1decayType == 13 ){
      //|| W1decayType == 1511 || W1decayType == 1513 ){
    //if(W1decayType == 1511 || W1decayType == 1513 ){
       //|| W2decayType == 11 || W2decayType == 13 || W2decayType == 1511 || W2decayType == 1513){

      //std::cout << "decaytype = " << decayType << ", W1decayType = " << W1decayType << ", W2decayType = " << W2decayType 
      //		<< ", W1 = " << W1 << ", W1 daughter = " << W1daughter << ", W2 = " << W2 << ", W2daughter = " << W2daughter << std::endl;
      
      //std::cout << "\t W1 daughter: pdgId = " << myGenParticles->at(W1daughter).pdgId << ", eta = " << myGenParticles->at(W1daughter).eta << ", pt = " << myGenParticles->at(W1daughter).pt  << std::endl;
      //std::cout << "\t W1 daughter+1: pdgId = " << myGenParticles->at(W1daughter+1).pdgId << ", eta = " << myGenParticles->at(W1daughter+1).eta << ", pt = " << myGenParticles->at(W1daughter+1).pt  << std::endl;

       int W1firstdaughter = TMath::Nint(myGenParticles->at(W1).firstDaughter);
       //std::cout << "W1 first daughter = " << W1firstdaughter << std::endl;
       int d1_id = abs(myGenParticles->at(W1firstdaughter).pdgId);
       int d2_id = abs(myGenParticles->at(W1firstdaughter+1).pdgId);
       float d1_eta = myGenParticles->at(W1firstdaughter).eta;
       float d2_eta = myGenParticles->at(W1firstdaughter+1).eta;
       float d1_pt = myGenParticles->at(W1firstdaughter).pt;
       float d2_pt = myGenParticles->at(W1firstdaughter+1).pt;
       //std::cout << "\t W1 daughter 1: pdgId = " << d1_id << ", eta = " << d1_eta << ", pt = " << d1_pt  << std::endl;
       //std::cout << "\t W1 daughter 2: pdgId = " << d2_id << ", eta = " << d2_eta << ", pt = " << d2_pt  << std::endl;

       float l_eta, nu_eta, l_pt, nu_pt;
       if( d1_id == 12 || d1_id == 14 || d1_id == 16){
	 nu_eta = d1_eta; nu_pt = d1_pt;
	 l_eta = d2_eta; l_pt = d2_pt;
       }
       else if( d2_id == 12 || d2_id == 14 || d2_id == 16){
	 nu_eta = d2_eta; nu_pt = d2_pt;
	 l_eta = d1_eta; l_pt = d1_pt;
       }
       else{std::cout << "error - expected a neutrino" << std::endl; assert(0);}
      //std::cout << "\t W2 daughter: pdgID = " << myGenParticles->at(W2daughter).pdgid << ", eta = " << myGenParticles->at(W2daughter).eta << ", pt = " << myGenParticles->at(W2daughter).pt  << std::endl;
      
       h_Wlnueta->Fill(l_eta, nu_eta);
       h_Wlnupt->Fill(l_pt, nu_pt);

    }
    
    if(W2decayType == 11 || W2decayType == 13 ){

       int W2firstdaughter = TMath::Nint(myGenParticles->at(W2).firstDaughter);
       //std::cout << "W1 first daughter = " << W1firstdaughter << std::endl;
       int d1_id = abs(myGenParticles->at(W2firstdaughter).pdgId);
       int d2_id = abs(myGenParticles->at(W2firstdaughter+1).pdgId);
       float d1_eta = myGenParticles->at(W2firstdaughter).eta;
       float d2_eta = myGenParticles->at(W2firstdaughter+1).eta;
       float d1_pt = myGenParticles->at(W2firstdaughter).pt;
       float d2_pt = myGenParticles->at(W2firstdaughter+1).pt;

       float l_eta, nu_eta, l_pt, nu_pt;
       if( d1_id == 12 || d1_id == 14 || d1_id == 16){
	 nu_eta = d1_eta; nu_pt = d1_pt;
	 l_eta = d2_eta; l_pt = d2_pt;
       }
       else if( d2_id == 12 || d2_id == 14 || d2_id == 16){
	 nu_eta = d2_eta; nu_pt = d2_pt;
	 l_eta = d1_eta; l_pt = d1_pt;
       }
       else{std::cout << "error - expected a neutrino" << std::endl; assert(0);}

      
       h_Wlnueta->Fill(l_eta, nu_eta);
       h_Wlnupt->Fill(l_pt, nu_pt);

    }

    //tautolepton
    if(W1decayType == 1511 || W1decayType == 1513 ){

      //W1 daughter is always the e/mu
      int wd_pdgid = abs(myGenParticles->at(W1daughter).pdgId);
      float wd_eta = myGenParticles->at(W1daughter).eta;
      float wd_pt = myGenParticles->at(W1daughter).pt;

      if(wd_pdgid != 11 && wd_pdgid !=13){std::cout << "error- expected an electron or muon" << std::endl; assert(0);}

      //std::cout << "decaytype = " << decayType << ", W1decayType = " << W1decayType << ", W2decayType = " << W2decayType 
      //		<< ", W1 = " << W1 << ", W1 daughter = " << W1daughter << ", W2 = " << W2 << ", W2daughter = " << W2daughter << std::endl;
      //std::cout << "\t W1 daughter: pdgId = " << wd_pdgid << ", eta = " << wd_eta << ", pt = " << wd_pt  << std::endl;
      //std::cout << "\t W1 daughter+1: pdgId = " << myGenParticles->at(W1daughter+1).pdgId << ", eta = " << myGenParticles->at(W1daughter+1).eta << ", pt = " << myGenParticles->at(W1daughter+1).pt  << std::endl;
      
      int W1firstdaughter = TMath::Nint(myGenParticles->at(W1).firstDaughter);
      //std::cout << "W1 first daughter = " << W1firstdaughter << std::endl;
      int d1_id = abs(myGenParticles->at(W1firstdaughter).pdgId);
      int d2_id = abs(myGenParticles->at(W1firstdaughter+1).pdgId);
      float d1_eta = myGenParticles->at(W1firstdaughter).eta;
      float d2_eta = myGenParticles->at(W1firstdaughter+1).eta;
      float d1_pt = myGenParticles->at(W1firstdaughter).pt;
      float d2_pt = myGenParticles->at(W1firstdaughter+1).pt;
      //std::cout << "\t W1 daughter 1: pdgId = " << d1_id << ", eta = " << d1_eta << ", pt = " << d1_pt  << std::endl;
      //std::cout << "\t W1 daughter 2: pdgId = " << d2_id << ", eta = " << d2_eta << ", pt = " << d2_pt  << std::endl;
         
      float l_eta, nu_eta, l_pt, nu_pt;
      l_eta = wd_eta; l_pt = wd_pt;
      if( d1_id == 16){
	nu_eta = d1_eta; nu_pt = d1_pt;
      }
      else if( d2_id == 16){
	nu_eta = d2_eta; nu_pt = d2_pt;
      }
      else{std::cout << "d1_id = " << d1_id << std::endl; std::cout << "error - expected a neutrino" << std::endl; assert(0);}
      ////std::cout << "\t W2 daughter: pdgID = " << myGenParticles->at(W2daughter).pdgid << ", eta = " << myGenParticles->at(W2daughter).eta << ", pt = " << myGenParticles->at(W2daughter).pt  << std::endl;
      
      h_Wtaulnueta->Fill(l_eta, nu_eta);
      h_Wtaulnupt->Fill(l_pt, nu_pt);
    }

    if(W2decayType == 1511 || W2decayType == 1513 ){
      
      //W1 daughter is always the e/mu
      int wd_pdgid = abs(myGenParticles->at(W2daughter).pdgId);
      float wd_eta = myGenParticles->at(W2daughter).eta;
      float wd_pt = myGenParticles->at(W2daughter).pt;

      if(wd_pdgid != 11 && wd_pdgid !=13){std::cout << "error- expected an electron or muon" << std::endl; assert(0);}

      //std::cout << "decaytype = " << decayType << ", W1decayType = " << W1decayType << ", W2decayType = " << W2decayType 
      //		<< ", W1 = " << W1 << ", W1 daughter = " << W1daughter << ", W2 = " << W2 << ", W2daughter = " << W2daughter << std::endl;
      //std::cout << "\t W1 daughter: pdgId = " << wd_pdgid << ", eta = " << wd_eta << ", pt = " << wd_pt  << std::endl;

      int W2firstdaughter = TMath::Nint(myGenParticles->at(W2).firstDaughter);
      //std::cout << "W1 first daughter = " << W1firstdaughter << std::endl;
      int d1_id = abs(myGenParticles->at(W2firstdaughter).pdgId);
      int d2_id = abs(myGenParticles->at(W2firstdaughter+1).pdgId);
      float d1_eta = myGenParticles->at(W2firstdaughter).eta;
      float d2_eta = myGenParticles->at(W2firstdaughter+1).eta;
      float d1_pt = myGenParticles->at(W2firstdaughter).pt;
      float d2_pt = myGenParticles->at(W2firstdaughter+1).pt;
      //std::cout << "\t W2 daughter 1: pdgId = " << d1_id << ", eta = " << d1_eta << ", pt = " << d1_pt  << std::endl;
      //std::cout << "\t W2 daughter 2: pdgId = " << d2_id << ", eta = " << d2_eta << ", pt = " << d2_pt  << std::endl;
   

      float l_eta, nu_eta, l_pt, nu_pt;
      l_eta = wd_eta; l_pt = wd_pt;
      if( d1_id == 16){
	nu_eta = d1_eta; nu_pt = d1_pt;
      }
      else if( d2_id == 16){
	nu_eta = d2_eta; nu_pt = d2_pt;
      }
      else{std::cout << "d1_id = " << d1_id << std::endl; std::cout << "error - expected a neutrino" << std::endl; assert(0);}
      ////std::cout << "\t W2 daughter: pdgID = " << myGenParticles->at(W2daughter).pdgid << ", eta = " << myGenParticles->at(W2daughter).eta << ", pt = " << myGenParticles->at(W2daughter).pt  << std::endl;
      
      h_Wtaulnueta->Fill(l_eta, nu_eta);
      h_Wtaulnupt->Fill(l_pt, nu_pt);
    }
    */

    
    //if((W1decayType==13 && W2decayType==1) || (W1decayType==1 && W2decayType==13) || (W1decayType==1513 && W2decayType==1) || (W1decayType==1 && W2decayType==1513)){
    //  nSemiMu++;
    //}
    //else if((W1decayType==11 && W2decayType==1) || (W1decayType==1 && W2decayType==11) || (W1decayType==1511 && W2decayType==1) || (W1decayType==1 && W2decayType==1511)){
    //  nSemiEle++;
    //}
    //else if((W1decayType==15 && W2decayType==1) || (W1decayType==1 && W2decayType==15)){
    //  nSemiTauHad++;
    //}
    //else if((W1decayType==11 || W1decayType==13 || W1decayType==1511 || W1decayType==1513) && (W2decayType==11 || W2decayType==13 || W2decayType==1511 || W2decayType==1513)){
    //  nDilep++;
    //}
    //else if(W1decayType==1 && W2decayType==1){
    //  nHad++;
    //}
    //else if((W1decayType==15 && W2decayType==15)||(W1decayType==11 && W2decayType==15)||(W1decayType==15 && W2decayType==11)||(W1decayType==1511 && W2decayType==15)||(W1decayType==15 && W2decayType==1511)||(W1decayType==13 && W2decayType==15)||(W1decayType==15 && W2decayType==13)||(W1decayType==15 && W2decayType==1513)||(W1decayType==1513 && W2decayType==15)){
    //  nOther++;
    //}
    //else{
    //  std::cout << "We got a problem boss" << std::endl;
    //  std::cout << "\t W1decayType = " << W1decayType << ", W2decayType = " << W2decayType << std::endl;
    //}

    //if(Cut()==1){


    //std::cout << "ngood jets = " << nGoodJets() << std::endl;
    //for (unsigned int i=0; i < jets_AK5PF_pt->size(); ++i) {
    //  if(getJetPt(i)<30) continue;
    //  if (isGoodJet30(i) ){
    //	std::cout << "found good jet" << std::endl;
    //  }
    //  else{
    //	std::cout << "jet failed id" << std::endl;
    //	std::cout << "\tneutralHadEnFrac (should be <0.99) = " << myJetsPF->at(i).neutralHadronEnergyFraction << std::endl;
    //	std::cout << "\tneutralEmEnFrac  (should be <0.99) = " << myJetsPF->at(i).neutralEmEnergyFraction << std::endl;
    //	std::cout << "\tnumDaughters (should be >1)        = " << myJetsPF->at(i).numberOfDaughters << std::endl;
    //	std::cout << "\tchargedHadEnFrac (should be >0)    = " << myJetsPF->at(i).chargedHadronEnergyFraction << std::endl;
    //	std::cout << "\tchargedEmEnFrac  (should be <0.99) = " << myJetsPF->at(i).chargedEmEnergyFraction << std::endl;
    //	std::cout << "\tchargedMult (should be >0)         = " << myJetsPF->at(i).chargedMultiplicity << std::endl;
    //
    //  }
    //  std::cout << "\tjet pt = " << getJetPt(i) << ", eta = " << myJetsPF->at(i).eta << ", phi = " << myJetsPF->at(i).phi;
    //  std::cout << ", CSV disc = " <<  myJetsPF->at(i).combinedSecondaryVertexBJetTags << std::endl;
    //
    //}


      npass++;

      //double effUp,effDown;     
      //if( (countMu()==1 && countEle()==0) || (countMu()==0 && countEle()==1)){
      //	std::cout << "HT = " << getHT() << ",MET = " << getMET() << ", eff = "
      //		  << getHLTMHTeffBNN(getMET(),getHT(),countEle(),countMu(),getMinDeltaPhiMETN(3),effUp,effDown) << std::endl;
      //	std::cout << ", effUp = " << effUp << ", effDown = " << effDown << std::endl;
      //}
      //std::cout << "MET = " << getMET() << ", eff = " << getHLTMHTeff(getMET()) << std::endl;
      //std::cout << "HT = " << getHT() << ", eff = " << getHLTHTeff(getHT()) << std::endl;

      /*
      bool cutEleVeto,cutMuVeto;
      cutEleVeto = passCut("cutEleVeto");
      cutMuVeto = passCut("cutMuVeto");
      int nElectrons, nMuons;
      nElectrons = countEle();
      nMuons = countMu();
                  
      //int lumi =  getLumiSection();
      int run = getRunNumber();
      //if(run ==170722 && lumi >=110 && lumi<=287) return false;
      if(passLumiMask()){
	int version = 0, prescale = 0;

	/////////////
	//HT300_MHT75
	/////////////
	if(run>=165922 && run<=166978){
	  if(passUtilityHLT(version,prescale)){
	    h_met_den_mht75->Fill(getMET());	
	    if(passHLT())
	      h_met_trig_mht75->Fill(getMET());

	    //0L sample
	    if(cutEleVeto && cutMuVeto){
	      h_met_0L_den_mht75->Fill(getMET());	
	      if(passHLT())
		h_met_0L_trig_mht75->Fill(getMET());
	    }
	    //1L sample
	    if( (nElectrons==1 && nMuons ==0)||(nElectrons==0 && nMuons==1) ){
	      h_met_1L_den_mht75->Fill(getMET());	
	      if(passHLT())
		h_met_1L_trig_mht75->Fill(getMET());
	    }
	    //1e0mu sample
	    if( nElectrons==1 && nMuons ==0 ){
	      h_met_1e0mu_den_mht75->Fill(getMET());	
	      if(passHLT())
		h_met_1e0mu_trig_mht75->Fill(getMET());	      
	    }
	    //1mu0e sample
	    if( nElectrons==0 && nMuons ==1){
	      h_met_1mu0e_den_mht75->Fill(getMET());	
	      if(passHLT())
		h_met_1mu0e_trig_mht75->Fill(getMET());
	    }

	  }
	}

	/////////////
	//HT300_MHT80
	/////////////
	if(run>=166979 && run<=173211){
	  if(passUtilityHLT(version,prescale)){
	    h_met_den_mht80->Fill(getMET());	
	    if(passHLT())
	      h_met_trig_mht80->Fill(getMET());


	    //0L sample
	    if(cutEleVeto && cutMuVeto){
	      h_met_0L_den_mht80->Fill(getMET());	
	      if(passHLT())
		h_met_0L_trig_mht80->Fill(getMET());
	    }
	    //1L sample
	    if( (nElectrons==1 && nMuons ==0)||(nElectrons==0 && nMuons==1) ){
	      h_met_1L_den_mht80->Fill(getMET());	
	      if(passHLT())
		h_met_1L_trig_mht80->Fill(getMET());
	    }
	    //1e0mu sample
	    if( nElectrons==1 && nMuons ==0 ){
	      h_met_1e0mu_den_mht80->Fill(getMET());	
	      if(passHLT())
		h_met_1e0mu_trig_mht80->Fill(getMET());	      
	    }
	    //1mu0e sample
	    if( nElectrons==0 && nMuons ==1){
	      h_met_1mu0e_den_mht80->Fill(getMET());	
	      if(passHLT())
		h_met_1mu0e_trig_mht80->Fill(getMET());
	    }


	  }
	}

	/////////////
	//HT300_MHT90
	/////////////

	if(run>=173212 && run<=176544){
	  if(passUtilityHLT(version,prescale)){
	    h_met_den_mht90_ht300->Fill(getMET());	
	    if(passHLT())
	      h_met_trig_mht90_ht300->Fill(getMET());


	    //0L sample
	    if(cutEleVeto && cutMuVeto){
	      h_met_0L_den_mht90_ht300->Fill(getMET());	
	      if(passHLT())
		h_met_0L_trig_mht90_ht300->Fill(getMET());
	    }
	    //1L sample
	    if( (nElectrons==1 && nMuons ==0)||(nElectrons==0 && nMuons==1) ){
	      h_met_1L_den_mht90_ht300->Fill(getMET());	
	      if(passHLT())
		h_met_1L_trig_mht90_ht300->Fill(getMET());
	    }
	    //1e0mu sample
	    if( nElectrons==1 && nMuons ==0 ){
	      h_met_1e0mu_den_mht90_ht300->Fill(getMET());	
	      if(passHLT())
		h_met_1e0mu_trig_mht90_ht300->Fill(getMET());	      
	    }
	    //1mu0e sample
	    if( nElectrons==0 && nMuons ==1){
	      h_met_1mu0e_den_mht90_ht300->Fill(getMET());	
	      if(passHLT())
		h_met_1mu0e_trig_mht90_ht300->Fill(getMET());
	    }


	  }
	}

	/////////////
	//HT350_MHT90
	/////////////

	if(run>=176545 && run<=178410){
	  if(passUtilityHLT(version,prescale)){
	    h_met_den_mht90_ht350->Fill(getMET());	
	    if(passHLT())
	      h_met_trig_mht90_ht350->Fill(getMET());

	    //0L sample
	    if(cutEleVeto && cutMuVeto){
	      h_met_0L_den_mht90_ht350->Fill(getMET());	
	      if(passHLT())
		h_met_0L_trig_mht90_ht350->Fill(getMET());
	    }
	    //1L sample
	    if( (nElectrons==1 && nMuons ==0)||(nElectrons==0 && nMuons==1) ){
	      h_met_1L_den_mht90_ht350->Fill(getMET());	
	      if(passHLT())
		h_met_1L_trig_mht90_ht350->Fill(getMET());
	    }
	    //1e0mu sample
	    if( nElectrons==1 && nMuons ==0 ){
	      h_met_1e0mu_den_mht90_ht350->Fill(getMET());	
	      if(passHLT())
		h_met_1e0mu_trig_mht90_ht350->Fill(getMET());	      
	    }
	    //1mu0e sample
	    if( nElectrons==0 && nMuons ==1){
	      h_met_1mu0e_den_mht90_ht350->Fill(getMET());	
	      if(passHLT())
		h_met_1mu0e_trig_mht90_ht350->Fill(getMET());
	    }



	  }
	}

	//////////////
	//HT350_MHT110
	//////////////

	if(run>=178411 && run<=180252){
	  if(passUtilityHLT(version,prescale)){
	    h_met_den_mht110->Fill(getMET());	
	    if(passHLT())
	      h_met_trig_mht110->Fill(getMET());

	    //0L sample
	    if(cutEleVeto && cutMuVeto){
	      h_met_0L_den_mht110->Fill(getMET());	
	      if(passHLT())
		h_met_0L_trig_mht110->Fill(getMET());
	    }
	    //1L sample
	    if( (nElectrons==1 && nMuons ==0)||(nElectrons==0 && nMuons==1) ){
	      h_met_1L_den_mht110->Fill(getMET());	
	      if(passHLT())
		h_met_1L_trig_mht110->Fill(getMET());
	    }
	    //1e0mu sample
	    if( nElectrons==1 && nMuons ==0 ){
	      h_met_1e0mu_den_mht110->Fill(getMET());	
	      if(passHLT())
		h_met_1e0mu_trig_mht110->Fill(getMET());	      
	    }
	    //1mu0e sample
	    if( nElectrons==0 && nMuons ==1){
	      h_met_1mu0e_den_mht110->Fill(getMET());	
	      if(passHLT())
		h_met_1mu0e_trig_mht110->Fill(getMET());
	    }

	  }
	}
      }      
      */
      /*          
      //h_ht->Fill(getHT());
      if(getHT()>400){
      h_met->Fill(getMET());
      
      int lumi =  getLumiSection();
      int run = getRunNumber();
      int event = getEventNumber();
      
      bool passTrig = false;
      for(uint i = 0; i< vevent.size(); ++i){
      	if( event == vevent.at(i)){
      	  if( lumi == vlumi.at(i)){
      	    if( run == vrun.at(i)){
      	      passTrig = true;
      	      break;
      	    }
      	  }
      	}
      }
      if(passTrig)
      	//h_ht_trig->Fill(getHT());
	h_met_trig->Fill(getMET());
      }
      */     


      /*               
      int npv = 0, npvbxp1 = 0, npvbxm1 = 0;
      for ( unsigned int i = 0; i<pileupsummaryinfo.size() ; i++) {
	//npv = pileupsummaryinfo.at(i).addpileupinfo_getPU_NumInteractions;
	//sum_nvtx += float(npv);
	
	//consider only in-time PU
	int BX = pileupsummaryinfo.at(i).addpileupinfo_getBunchCrossing;
	if(BX == 0) { 
	  npv = pileupsummaryinfo.at(i).addpileupinfo_getPU_NumInteractions;
	  //continue;
	}
	else if(BX == 1) { 
	  npvbxp1 = pileupsummaryinfo.at(i).addpileupinfo_getPU_NumInteractions;
	  //continue;
	}
	else if(BX == -1) { 
	  npvbxm1 = pileupsummaryinfo.at(i).addpileupinfo_getPU_NumInteractions;
	  //continue;
	}
	
      }
      //std::cout << "npv = " << npv << std::endl;
      h_MCPU->Fill(npv);
      //std::cout << "weight = " << getPUWeight(LumiWeights) << std::endl;
      h_MCPUr->Fill(npv, getPUWeight(LumiWeights) );
	
      h_MCPUBXp1->Fill(npvbxp1);
      h_MCPUrBXp1->Fill(npvbxp1, getPUWeight(LumiWeights) );
      h_MCPUBXm1->Fill(npvbxm1);
      h_MCPUrBXm1->Fill(npvbxm1, getPUWeight(LumiWeights) );
      */

      //calculateTagProb(prob0,probge1,prob1,probge2,probge3);
      //n0b += prob0; 
      //nge1b += probge1; 
      //neq1b += prob1; 
      //nge2b += probge2;      
      //nge3b += probge3;

      //}//end Cut

  }


  cout<<endl;
  stopTimer(nevents);

  //std::cout << "n0b = " << n0b << std::endl;
  //std::cout << "nge1b = " << nge1b << std::endl;
  //std::cout << "neq1b = " << neq1b << std::endl;
  //std::cout << "nge2b = " << nge2b << std::endl;
  //std::cout << "nge3b = " << nge3b << std::endl;
  
  //std::cout << "npass = " << npass << std::endl;
  
  std::cout << "nSemiMu = " << nSemiMu << std::endl;  
  std::cout << "nSemiEle = " << nSemiEle << std::endl;  
  std::cout << "nSemiTauHad = " << nSemiTauHad << std::endl;  
  std::cout << "nDilep = " << nDilep << std::endl;  
  std::cout << "nHad = " << nHad << std::endl;  
  std::cout << "nOther = " << nOther << std::endl;  
  
  
  fout.Write();
  fout.Close();


  return;
}

TString EventCalculator::assembleBTagEffFilename(bool cutnametail) {

  TString outfile = "histos_btageff_csvm_";
  TString basesamplename = sampleName_;
  if (cutnametail) {
    //goal: take e.g. "SMS-T1bbbb_Mgluino-100to2000_mLSP-0to2000_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v1_AODSIM_UCSB1548_v66.9"
    //and chop off the end to produce "SMS-T1bbbb_Mgluino-100to2000_mLSP-0to2000_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v1_AODSIM_UCSB1548_v66"
    int npieces= basesamplename.Tokenize(".")->GetEntries();
    if (npieces>1) {
      //make sure that this last piece is an index. otherwise we might be in a situation we did not expect
      assert( TString(basesamplename.Tokenize(".")->At(npieces-1)->GetName()).IsDigit());
      //we want to chop off the end. how many digits are there?
      int sizetochop = TString(basesamplename.Tokenize(".")->At(npieces-1)->GetName()).Length();
      ++sizetochop; //include the length of the period
      for ( int ii=0;ii<sizetochop;ii++) basesamplename.Chop(); //chop off the end
    }
  }
  outfile+=basesamplename;
  outfile+=".root";
  return outfile;
}

//separate out the jettag efficiency code from sampleAnalyzer for ease of use
void EventCalculator::plotBTagEffMC( ) {

  TString outfile = assembleBTagEffFilename();

  TFile fout(outfile,"RECREATE");
  
  double bins[16] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 670, 2000};
  TH1F * h_bjet = new TH1F("h_bjet","bjet",15,bins);
  TH1F * h_cjet = new TH1F("h_cjet","cjet",15,bins);
  TH1F * h_ljet = new TH1F("h_ljet","ljet",15,bins);
  TH1F * h_btag = new TH1F("h_btag","btag",15,bins);
  TH1F * h_ctag = new TH1F("h_ctag","ctag",15,bins);
  TH1F * h_ltag = new TH1F("h_ltag","ltag",15,bins);
  h_bjet->Sumw2();
  h_cjet->Sumw2();
  h_ljet->Sumw2();
  h_btag->Sumw2();
  h_ctag->Sumw2();
  h_ltag->Sumw2();


  const Long64_t nevents = chainA->GetEntries();
  const Long64_t neventsB = chainB->GetEntries();
  assert(nevents==neventsB);


//   setCutScheme();
//   setIgnoredCut("cut1b");
//   setIgnoredCut("cut2b");
//   setIgnoredCut("cut3b");

  cout<<"Running..."<<endl;  
  int npass = 0;

  int ntaggedjets = 0;
  int ntaggedjets_b = 0;
  int ntaggedjets_c = 0;
  int ntaggedjets_l = 0;

  //this threshold does not need to be synced with the global b-tag pt threshold, because this histogram is binned in pT
  const float jetptthreshold = 30;

  startTimer();
  for(Long64_t entry=0; entry < nevents; ++entry){
    chainB->GetEntry(entry);
    chainA->GetEntry(entry);

    if(entry%10000==0) cout << "entry: " << entry << ", percent done=" << (int)(entry/(double)nevents*100.)<<  endl;
    
    npass++;
    
    for (unsigned int i = 0; i < jets_AK5PF_pt->size(); ++i) {
      int flavor = jets_AK5PF_partonFlavour->at(i);

            
      if(isGoodJet(i,jetptthreshold)){
      
	  
	if(abs(flavor)==5)
	  h_bjet->Fill( getJetPt(i) ) ;
	else if(abs(flavor)==4)
	  h_cjet->Fill( getJetPt(i) ) ;
	else if(abs(flavor)==3 ||abs(flavor)==2 ||abs(flavor)==1 ||abs(flavor)==21 )
	  h_ljet->Fill( getJetPt(i) ) ;
      
	if(passBTagger(i)){
	  if(abs(flavor)==5)
	    h_btag->Fill( getJetPt(i) ) ;
	  else if(abs(flavor)==4)
	    h_ctag->Fill( getJetPt(i) ) ;
	  else if(abs(flavor)==3 ||abs(flavor)==2 ||abs(flavor)==1 ||abs(flavor)==21 )
	    h_ltag->Fill( getJetPt(i) ) ;
	}
      
      }
	
      if (isGoodJet(i,jetptthreshold) && passBTagger(i) ){
	ntaggedjets++;
	  
	if(abs(flavor)==5)
	  ntaggedjets_b++;
	else if(abs(flavor)==4)
	  ntaggedjets_c++;
	else if(abs(flavor)==3 ||abs(flavor)==2 ||abs(flavor)==1 ||abs(flavor)==21 )
	  ntaggedjets_l++;
      }
	
    }
      
  }
  
  cout<<endl;
  stopTimer(nevents);

  std::cout << "npass = " << npass << std::endl;
  std::cout << "ntaggedjets = "   << ntaggedjets << std::endl;
  std::cout << "ntaggedjets_b = " << ntaggedjets_b << std::endl;
  std::cout << "ntaggedjets_c = " << ntaggedjets_c << std::endl;
  std::cout << "ntaggedjets_l = " << ntaggedjets_l << std::endl;
    
  //TH1F *h_btageff = (TH1F*) h_btag->Clone();
  TH1F * h_btageff = new TH1F("h_btageff","btageff",15,bins);
  h_btageff->Divide(h_btag,h_bjet,1,1,"B");
  //TH1F *h_ctageff = (TH1F*) h_ctag->Clone();
  TH1F * h_ctageff = new TH1F("h_ctageff","ctageff",15,bins);
  h_ctageff->Divide(h_ctag,h_cjet,1,1,"B");
  //TH1F *h_ltageff = (TH1F*) h_ltag->Clone();
  TH1F * h_ltageff = new TH1F("h_ltageff","ltageff",15,bins);
  h_ltageff->Divide(h_ltag,h_ljet,1,1,"B");
  
  
  fout.Write();
  fout.Close();


  return;
}


void EventCalculator::loadECALStatus() {
  cout << "loading ECAL Status" << endl;
  assert(badECAL_.size()==0);//don't want to load it twice
  
  TChain *tDetector = new TChain("tree_detector");
  tDetector->Add("ECALStatus_QCD.root");//this file is for QCD MC - only one run.
  //will need to make code more sophisticated if input file contains multiple runs
  
  double td_eta, td_phi;
  int td_status;
  //unsigned int td_run;
  //int td_det;
  //int td_status_simp;
  
  tDetector->SetBranchAddress("eta",&td_eta);
  tDetector->SetBranchAddress("phi",&td_phi);
  tDetector->SetBranchAddress("status",&td_status);
  //tDetector->SetBranchAddress("run",&td_run);
  //tDetector->SetBranchAddress("det",&td_det);
  //tDetector->SetBranchAddress("status_simp",&td_status_simp);

  const int numD = tDetector->GetEntries();
  for(int i=0; i<numD; i++) {
    tDetector->GetEvent(i);
   
    if(td_status>=10) {
      cellECAL thisCell(td_eta,td_phi,td_status);
      badECAL_.push_back(thisCell);
    }
 
  }//end loop over input tree

  return;
}


bool EventCalculator::passBadECALFilter(TString nearmettype, double nearmetphicut, TString nearecaltype, double nearecalRcut, int nbadcellsthr, int badcellstatusthr) {
  
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    if (getJetPt(i) > 30 && fabs(jets_AK5PF_eta->at(i)) < 5) {
      
      if( jetNearMET(i, nearmettype, nearmetphicut) ) {
	if( jetNearBadECALCell(i, nearecaltype, nearecalRcut, nbadcellsthr, badcellstatusthr) ) return 0;
      }
    }
  }
  return 1;
}


float EventCalculator::passBadECALFilter_worstMis(TString nearmettype, double nearmetphicut, TString nearecaltype, double nearecalRcut, int nbadcellsthr, int badcellstatusthr, TString mistype) {
  
  float worst = -9e9;//assume undermeasurement!
  if (isSampleRealData()) return worst;

  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    if (getJetPt(i) > 30 && fabs(jets_AK5PF_eta->at(i)) < 5) {
      
      if( jetNearMET(i, nearmettype, nearmetphicut) ) {
	if( jetNearBadECALCell(i, nearecaltype, nearecalRcut, nbadcellsthr, badcellstatusthr) ) {

	  float genpt = jets_AK5PF_gen_pt->at(i);
	  float pt = getJetPt(i);
	  
	  float thisworst;
	  if(mistype=="abs") {
	    thisworst = genpt - pt;//large for undermeasurement
	  }
	  else if(mistype=="frac") {
	    thisworst = genpt/pt;//large for undermeasurement
	  }
	  else assert(0);
	  
	  if(thisworst>worst) worst=thisworst;
	  
	}
      }
    }
  }

  return worst;
}

  
bool EventCalculator::jetNearMET(unsigned int i, TString type, double nearmetphicut) {
  
  if(type=="METphi") {
    double dp =  getDeltaPhi( jets_AK5PF_phi->at(i) , getMETphi());
    if(dp<=nearmetphicut) return 1;
    else return 0;
  }
  else assert(0);//no other types yet. future: RA1 deltaPhi* method
  
}

bool EventCalculator::jetNearBadECALCell(unsigned int i, TString type, double nearecalRcut, int nbadcellsthr, int badcellstatusthr) {
  double jeteta= jets_AK5PF_eta->at(i);

  if(type=="deadCell") {
    if (!sampleName_.Contains("QCD")) return 0;

    double jetphi= jets_AK5PF_phi->at(i);
    
    int nbad =0;
    for(unsigned int j=0; j<badECAL_.size(); j++) {
      if( (jmt::deltaR(jeteta,jetphi,badECAL_[j].eta,badECAL_[j].phi) <= nearecalRcut) && badECAL_[j].status>= badcellstatusthr) nbad++;
    }
    if(nbad>=nbadcellsthr) return 1;
    else return 0;
    
  }
  else if(type=="crackEBEE") {
    
    if( fabs(fabs(jeteta) - 1.5) <= nearecalRcut ) return 1;
    else return 0;

  }
  else if(type=="crackEEEF") {

    if( fabs(fabs(jeteta) - 3.0) <= nearecalRcut ) return 1;
    else return 0;

  }
  else assert(0);

}


void EventCalculator::InitializeA(TChain *fChain) 
{
   // Set object pointer
   trigger_prescalevalue = 0;
   trigger_name = 0;
   trigger_decision = 0;
   trigger_lastfiltername = 0;
   triggerobject_pt = 0;
   triggerobject_px = 0;
   triggerobject_py = 0;
   triggerobject_pz = 0;
   triggerobject_et = 0;
   triggerobject_energy = 0;
   triggerobject_phi = 0;
   triggerobject_eta = 0;
   triggerobject_collectionname = 0;
   standalone_triggerobject_pt = 0;
   standalone_triggerobject_px = 0;
   standalone_triggerobject_py = 0;
   standalone_triggerobject_pz = 0;
   standalone_triggerobject_et = 0;
   standalone_triggerobject_energy = 0;
   standalone_triggerobject_phi = 0;
   standalone_triggerobject_eta = 0;
   standalone_triggerobject_collectionname = 0;
   L1trigger_bit = 0;
   L1trigger_techTrigger = 0;
   L1trigger_prescalevalue = 0;
   L1trigger_name = 0;
   L1trigger_alias = 0;
   L1trigger_decision = 0;
   L1trigger_decision_nomask = 0;
   els_conversion_dist = 0;
   els_conversion_dcot = 0;
   els_PFchargedHadronIsoR03 = 0;
   els_PFphotonIsoR03 = 0;
   els_PFneutralHadronIsoR03 = 0;
   els_hasMatchedConversion = 0;
   pf_els_PFchargedHadronIsoR03 = 0;
   pf_els_PFphotonIsoR03 = 0;
   pf_els_PFneutralHadronIsoR03 = 0;
   pf_els_hasMatchedConversion = 0;
   jets_AK5PFclean_corrL2L3 = 0;
   jets_AK5PFclean_corrL2L3Residual = 0;
   jets_AK5PFclean_corrL1FastL2L3 = 0;
   jets_AK5PFclean_corrL1L2L3 = 0;
   jets_AK5PFclean_corrL1FastL2L3Residual = 0;
   jets_AK5PFclean_corrL1L2L3Residual = 0;
   jets_AK5PFclean_Uncert = 0;
   PU_zpositions = 0;
   PU_sumpT_lowpT = 0;
   PU_sumpT_highpT = 0;
   PU_ntrks_lowpT = 0;
   PU_ntrks_highpT = 0;
   PU_NumInteractions = 0;
   PU_bunchCrossing = 0;
   PU_TrueNumInteractions = 0;

   fChain->SetBranchAddress("trigger_prescalevalue", &trigger_prescalevalue, &b_trigger_prescalevalue);
   fChain->SetBranchAddress("trigger_name", &trigger_name, &b_trigger_name);
   fChain->SetBranchAddress("trigger_decision", &trigger_decision, &b_trigger_decision);
   fChain->SetBranchAddress("trigger_lastfiltername", &trigger_lastfiltername, &b_trigger_lastfiltername);
   fChain->SetBranchAddress("triggerobject_pt", &triggerobject_pt, &b_triggerobject_pt);
   fChain->SetBranchAddress("triggerobject_px", &triggerobject_px, &b_triggerobject_px);
   fChain->SetBranchAddress("triggerobject_py", &triggerobject_py, &b_triggerobject_py);
   fChain->SetBranchAddress("triggerobject_pz", &triggerobject_pz, &b_triggerobject_pz);
   fChain->SetBranchAddress("triggerobject_et", &triggerobject_et, &b_triggerobject_et);
   fChain->SetBranchAddress("triggerobject_energy", &triggerobject_energy, &b_triggerobject_energy);
   fChain->SetBranchAddress("triggerobject_phi", &triggerobject_phi, &b_triggerobject_phi);
   fChain->SetBranchAddress("triggerobject_eta", &triggerobject_eta, &b_triggerobject_eta);
   fChain->SetBranchAddress("triggerobject_collectionname", &triggerobject_collectionname, &b_triggerobject_collectionname);
   fChain->SetBranchAddress("standalone_triggerobject_pt", &standalone_triggerobject_pt, &b_standalone_triggerobject_pt);
   fChain->SetBranchAddress("standalone_triggerobject_px", &standalone_triggerobject_px, &b_standalone_triggerobject_px);
   fChain->SetBranchAddress("standalone_triggerobject_py", &standalone_triggerobject_py, &b_standalone_triggerobject_py);
   fChain->SetBranchAddress("standalone_triggerobject_pz", &standalone_triggerobject_pz, &b_standalone_triggerobject_pz);
   fChain->SetBranchAddress("standalone_triggerobject_et", &standalone_triggerobject_et, &b_standalone_triggerobject_et);
   fChain->SetBranchAddress("standalone_triggerobject_energy", &standalone_triggerobject_energy, &b_standalone_triggerobject_energy);
   fChain->SetBranchAddress("standalone_triggerobject_phi", &standalone_triggerobject_phi, &b_standalone_triggerobject_phi);
   fChain->SetBranchAddress("standalone_triggerobject_eta", &standalone_triggerobject_eta, &b_standalone_triggerobject_eta);
   fChain->SetBranchAddress("standalone_triggerobject_collectionname", &standalone_triggerobject_collectionname, &b_standalone_triggerobject_collectionname);
   fChain->SetBranchAddress("L1trigger_bit", &L1trigger_bit, &b_L1trigger_bit);
   fChain->SetBranchAddress("L1trigger_techTrigger", &L1trigger_techTrigger, &b_L1trigger_techTrigger);
   fChain->SetBranchAddress("L1trigger_prescalevalue", &L1trigger_prescalevalue, &b_L1trigger_prescalevalue);
   fChain->SetBranchAddress("L1trigger_name", &L1trigger_name, &b_L1trigger_name);
   fChain->SetBranchAddress("L1trigger_alias", &L1trigger_alias, &b_L1trigger_alias);
   fChain->SetBranchAddress("L1trigger_decision", &L1trigger_decision, &b_L1trigger_decision);
   fChain->SetBranchAddress("L1trigger_decision_nomask", &L1trigger_decision_nomask, &b_L1trigger_decision_nomask);
   fChain->SetBranchAddress("els_conversion_dist", &els_conversion_dist, &b_els_conversion_dist);
   fChain->SetBranchAddress("els_conversion_dcot", &els_conversion_dcot, &b_els_conversion_dcot);
   fChain->SetBranchAddress("els_PFchargedHadronIsoR03", &els_PFchargedHadronIsoR03, &b_els_PFchargedHadronIsoR03);
   fChain->SetBranchAddress("els_PFphotonIsoR03", &els_PFphotonIsoR03, &b_els_PFphotonIsoR03);
   fChain->SetBranchAddress("els_PFneutralHadronIsoR03", &els_PFneutralHadronIsoR03, &b_els_PFneutralHadronIsoR03);
   fChain->SetBranchAddress("els_hasMatchedConversion", &els_hasMatchedConversion, &b_els_hasMatchedConversion);
   fChain->SetBranchAddress("pf_els_PFchargedHadronIsoR03", &pf_els_PFchargedHadronIsoR03, &b_pf_els_PFchargedHadronIsoR03);
   fChain->SetBranchAddress("pf_els_PFphotonIsoR03", &pf_els_PFphotonIsoR03, &b_pf_els_PFphotonIsoR03);
   fChain->SetBranchAddress("pf_els_PFneutralHadronIsoR03", &pf_els_PFneutralHadronIsoR03, &b_pf_els_PFneutralHadronIsoR03);
   fChain->SetBranchAddress("pf_els_hasMatchedConversion", &pf_els_hasMatchedConversion, &b_pf_els_hasMatchedConversion);
   fChain->SetBranchAddress("hbhefilter_decision", &hbhefilter_decision, &b_hbhefilter_decision);
   fChain->SetBranchAddress("trackingfailurefilter_decision", &trackingfailurefilter_decision, &b_trackingfailurefilter_decision);
   fChain->SetBranchAddress("cschalofilter_decision", &cschalofilter_decision, &b_cschalofilter_decision);
   fChain->SetBranchAddress("ecalTPfilter_decision", &ecalTPfilter_decision, &b_ecalTPfilter_decision);
   fChain->SetBranchAddress("ecalBEfilter_decision", &ecalBEfilter_decision, &b_ecalBEfilter_decision);
   fChain->SetBranchAddress("scrapingVeto_decision", &scrapingVeto_decision, &b_scrapingVeto_decision);
   fChain->SetBranchAddress("greedymuonfilter_decision", &greedymuonfilter_decision, &b_greedymuonfilter_decision);
   fChain->SetBranchAddress("inconsistentPFmuonfilter_decision", &inconsistentPFmuonfilter_decision, &b_inconsistentPFmuonfilter_decision);
   fChain->SetBranchAddress("hcallaserfilter_decision", &hcallaserfilter_decision, &b_hcallaserfilter_decision);
   fChain->SetBranchAddress("eenoisefilter_decision", &eenoisefilter_decision, &b_eenoisefilter_decision);
   fChain->SetBranchAddress("eebadscfilter_decision", &eebadscfilter_decision, &b_eebadscfilter_decision);
   fChain->SetBranchAddress("passprescalePFHT350filter_decision", &passprescalePFHT350filter_decision, &b_passprescalePFHT350filter_decision);
   fChain->SetBranchAddress("passprescaleHT250filter_decision", &passprescaleHT250filter_decision, &b_passprescaleHT250filter_decision);
   fChain->SetBranchAddress("passprescaleHT300filter_decision", &passprescaleHT300filter_decision, &b_passprescaleHT300filter_decision);
   fChain->SetBranchAddress("passprescaleHT350filter_decision", &passprescaleHT350filter_decision, &b_passprescaleHT350filter_decision);
   fChain->SetBranchAddress("passprescaleHT400filter_decision", &passprescaleHT400filter_decision, &b_passprescaleHT400filter_decision);
   fChain->SetBranchAddress("passprescaleHT450filter_decision", &passprescaleHT450filter_decision, &b_passprescaleHT450filter_decision);
   fChain->SetBranchAddress("MPT", &MPT, &b_MPT);
   fChain->SetBranchAddress("genHT", &genHT, &b_genHT);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL2L3", &jets_AK5PFclean_corrL2L3, &b_jets_AK5PFclean_corrL2L3);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL2L3Residual", &jets_AK5PFclean_corrL2L3Residual, &b_jets_AK5PFclean_corrL2L3Residual);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL1FastL2L3", &jets_AK5PFclean_corrL1FastL2L3, &b_jets_AK5PFclean_corrL1FastL2L3);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL1L2L3", &jets_AK5PFclean_corrL1L2L3, &b_jets_AK5PFclean_corrL1L2L3);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL1FastL2L3Residual", &jets_AK5PFclean_corrL1FastL2L3Residual, &b_jets_AK5PFclean_corrL1FastL2L3Residual);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL1L2L3Residual", &jets_AK5PFclean_corrL1L2L3Residual, &b_jets_AK5PFclean_corrL1L2L3Residual);
   fChain->SetBranchAddress("jets_AK5PFclean_Uncert", &jets_AK5PFclean_Uncert, &b_jets_AK5PFclean_Uncert);
   fChain->SetBranchAddress("PU_zpositions", &PU_zpositions, &b_PU_zpositions);
   fChain->SetBranchAddress("PU_sumpT_lowpT", &PU_sumpT_lowpT, &b_PU_sumpT_lowpT);
   fChain->SetBranchAddress("PU_sumpT_highpT", &PU_sumpT_highpT, &b_PU_sumpT_highpT);
   fChain->SetBranchAddress("PU_ntrks_lowpT", &PU_ntrks_lowpT, &b_PU_ntrks_lowpT);
   fChain->SetBranchAddress("PU_ntrks_highpT", &PU_ntrks_highpT, &b_PU_ntrks_highpT);
   fChain->SetBranchAddress("PU_NumInteractions", &PU_NumInteractions, &b_PU_NumInteractions);
   fChain->SetBranchAddress("PU_bunchCrossing", &PU_bunchCrossing, &b_PU_bunchCrossing);
   fChain->SetBranchAddress("PU_TrueNumInteractions", &PU_TrueNumInteractions, &b_PU_TrueNumInteractions);
   fChain->SetBranchAddress("rho_kt6PFJetsForIsolation2011", &rho_kt6PFJetsForIsolation2011, &b_rho_kt6PFJetsForIsolation2011);
   fChain->SetBranchAddress("rho_kt6PFJetsForIsolation2012", &rho_kt6PFJetsForIsolation2012, &b_rho_kt6PFJetsForIsolation2012);
   fChain->SetBranchAddress("pfmets_fullSignif", &pfmets_fullSignif, &b_pfmets_fullSignif);
   fChain->SetBranchAddress("pfmets_fullSignifCov00", &pfmets_fullSignifCov00, &b_pfmets_fullSignifCov00);
   fChain->SetBranchAddress("pfmets_fullSignifCov10", &pfmets_fullSignifCov10, &b_pfmets_fullSignifCov10);
   fChain->SetBranchAddress("pfmets_fullSignifCov11", &pfmets_fullSignifCov11, &b_pfmets_fullSignifCov11);
   fChain->SetBranchAddress("softjetUp_dMEx", &softjetUp_dMEx, &b_softjetUp_dMEx);
   fChain->SetBranchAddress("softjetUp_dMEy", &softjetUp_dMEy, &b_softjetUp_dMEy);


   //special! added by hand for signal
   pdfweights_cteq=0;
   pdfweights_mstw=0;
   pdfweights_nnpdf=0;
   if ( sampleIsSignal_) {
     cout<<"Signal sample -- setting branch addresses for pdf weights ..."<<endl;
     fChain->SetBranchAddress("pdfweights_cteq", &pdfweights_cteq, &b_pdfweights_cteq);
     fChain->SetBranchAddress("pdfweights_mstw", &pdfweights_mstw, &b_pdfweights_mstw);
     fChain->SetBranchAddress("pdfweights_nnpdf", &pdfweights_nnpdf, &b_pdfweights_nnpdf);
   }

}

void EventCalculator::InitializeB(TChain *fChain)
{

   // Set object pointer
   beamSpot_x = 0;
   beamSpot_y = 0;
   beamSpot_z = 0;
   beamSpot_x0Error = 0;
   beamSpot_y0Error = 0;
   beamSpot_z0Error = 0;
   beamSpot_sigmaZ = 0;
   beamSpot_sigmaZ0Error = 0;
   beamSpot_dxdz = 0;
   beamSpot_dxdzError = 0;
   beamSpot_dydz = 0;
   beamSpot_dydzError = 0;
   beamSpot_beamWidthX = 0;
   beamSpot_beamWidthY = 0;
   beamSpot_beamWidthXError = 0;
   beamSpot_beamWidthYError = 0;
   els_energy = 0;
   els_et = 0;
   els_eta = 0;
   els_phi = 0;
   els_pt = 0;
   els_px = 0;
   els_py = 0;
   els_pz = 0;
   els_status = 0;
   els_theta = 0;
   els_gen_id = 0;
   els_gen_phi = 0;
   els_gen_pt = 0;
   els_gen_pz = 0;
   els_gen_px = 0;
   els_gen_py = 0;
   els_gen_eta = 0;
   els_gen_theta = 0;
   els_gen_et = 0;
   els_gen_mother_id = 0;
   els_gen_mother_phi = 0;
   els_gen_mother_pt = 0;
   els_gen_mother_pz = 0;
   els_gen_mother_px = 0;
   els_gen_mother_py = 0;
   els_gen_mother_eta = 0;
   els_gen_mother_theta = 0;
   els_gen_mother_et = 0;
   els_tightId = 0;
   els_looseId = 0;
   els_robustTightId = 0;
   els_robustLooseId = 0;
   els_robustHighEnergyId = 0;
   els_simpleEleId95relIso = 0;
   els_simpleEleId90relIso = 0;
   els_simpleEleId85relIso = 0;
   els_simpleEleId80relIso = 0;
   els_simpleEleId70relIso = 0;
   els_simpleEleId60relIso = 0;
   els_simpleEleId95cIso = 0;
   els_simpleEleId90cIso = 0;
   els_simpleEleId85cIso = 0;
   els_simpleEleId80cIso = 0;
   els_simpleEleId70cIso = 0;
   els_simpleEleId60cIso = 0;
   els_cIso = 0;
   els_tIso = 0;
   els_ecalIso = 0;
   els_hcalIso = 0;
   els_chi2 = 0;
   els_charge = 0;
   els_caloEnergy = 0;
   els_hadOverEm = 0;
   els_hcalOverEcalBc = 0;
   els_eOverPIn = 0;
   els_eSeedOverPOut = 0;
   els_sigmaEtaEta = 0;
   els_sigmaIEtaIEta = 0;
   els_scEnergy = 0;
   els_scRawEnergy = 0;
   els_scSeedEnergy = 0;
   els_scEta = 0;
   els_scPhi = 0;
   els_scEtaWidth = 0;
   els_scPhiWidth = 0;
   els_scE1x5 = 0;
   els_scE2x5Max = 0;
   els_scE5x5 = 0;
   els_isEB = 0;
   els_isEE = 0;
   els_dEtaIn = 0;
   els_dPhiIn = 0;
   els_dEtaOut = 0;
   els_dPhiOut = 0;
   els_numvalhits = 0;
   els_numlosthits = 0;
   els_basicClustersSize = 0;
   els_tk_pz = 0;
   els_tk_pt = 0;
   els_tk_phi = 0;
   els_tk_eta = 0;
   els_d0dum = 0;
   els_dz = 0;
   els_vx = 0;
   els_vy = 0;
   els_vz = 0;
   els_ndof = 0;
   els_ptError = 0;
   els_d0dumError = 0;
   els_dzError = 0;
   els_etaError = 0;
   els_phiError = 0;
   els_tk_charge = 0;
   els_core_ecalDrivenSeed = 0;
   els_n_inner_layer = 0;
   els_n_outer_layer = 0;
   els_ctf_tk_id = 0;
   els_ctf_tk_charge = 0;
   els_ctf_tk_eta = 0;
   els_ctf_tk_phi = 0;
   els_fbrem = 0;
   els_shFracInnerHits = 0;
   els_dr03EcalRecHitSumEt = 0;
   els_dr03HcalTowerSumEt = 0;
   els_dr03HcalDepth1TowerSumEt = 0;
   els_dr03HcalDepth2TowerSumEt = 0;
   els_dr03TkSumPt = 0;
   els_dr04EcalRecHitSumEt = 0;
   els_dr04HcalTowerSumEt = 0;
   els_dr04HcalDepth1TowerSumEt = 0;
   els_dr04HcalDepth2TowerSumEt = 0;
   els_dr04TkSumPt = 0;
   els_cpx = 0;
   els_cpy = 0;
   els_cpz = 0;
   els_vpx = 0;
   els_vpy = 0;
   els_vpz = 0;
   els_cx = 0;
   els_cy = 0;
   els_cz = 0;
   els_PATpassConversionVeto = 0;
   jets_AK5PF_status = 0;
   jets_AK5PF_phi = 0;
   jets_AK5PF_pt = 0;
   jets_AK5PF_pz = 0;
   jets_AK5PF_px = 0;
   jets_AK5PF_py = 0;
   jets_AK5PF_eta = 0;
   jets_AK5PF_theta = 0;
   jets_AK5PF_et = 0;
   jets_AK5PF_energy = 0;
   jets_AK5PF_parton_Id = 0;
   jets_AK5PF_parton_motherId = 0;
   jets_AK5PF_parton_pt = 0;
   jets_AK5PF_parton_phi = 0;
   jets_AK5PF_parton_eta = 0;
   jets_AK5PF_parton_Energy = 0;
   jets_AK5PF_parton_mass = 0;
   jets_AK5PF_gen_et = 0;
   jets_AK5PF_gen_pt = 0;
   jets_AK5PF_gen_eta = 0;
   jets_AK5PF_gen_phi = 0;
   jets_AK5PF_gen_mass = 0;
   jets_AK5PF_gen_Energy = 0;
   jets_AK5PF_gen_Id = 0;
   jets_AK5PF_gen_motherID = 0;
   jets_AK5PF_gen_threeCharge = 0;
   jets_AK5PF_partonFlavour = 0;
   jets_AK5PF_btag_TC_highPur = 0;
   jets_AK5PF_btag_TC_highEff = 0;
   jets_AK5PF_btag_jetProb = 0;
   jets_AK5PF_btag_jetBProb = 0;
   jets_AK5PF_btag_softEle = 0;
   jets_AK5PF_btag_softMuon = 0;
   jets_AK5PF_btag_secVertexHighPur = 0;
   jets_AK5PF_btag_secVertexHighEff = 0;
   jets_AK5PF_btag_secVertexCombined = 0;
   jets_AK5PF_jetCharge = 0;
   jets_AK5PF_chgEmE = 0;
   jets_AK5PF_chgHadE = 0;
   jets_AK5PF_photonEnergy = 0;
   jets_AK5PF_chgMuE = 0;
   jets_AK5PF_chg_Mult = 0;
   jets_AK5PF_neutralEmE = 0;
   jets_AK5PF_neutralHadE = 0;
   jets_AK5PF_neutral_Mult = 0;
   jets_AK5PF_mu_Mult = 0;
   jets_AK5PF_emf = 0;
   jets_AK5PF_ehf = 0;
   jets_AK5PF_n60 = 0;
   jets_AK5PF_n90 = 0;
   jets_AK5PF_etaetaMoment = 0;
   jets_AK5PF_etaphiMoment = 0;
   jets_AK5PF_phiphiMoment = 0;
   jets_AK5PF_n90Hits = 0;
   jets_AK5PF_fHPD = 0;
   jets_AK5PF_fRBX = 0;
   jets_AK5PF_hitsInN90 = 0;
   jets_AK5PF_nECALTowers = 0;
   jets_AK5PF_nHCALTowers = 0;
   jets_AK5PF_fSubDetector1 = 0;
   jets_AK5PF_fSubDetector2 = 0;
   jets_AK5PF_fSubDetector3 = 0;
   jets_AK5PF_fSubDetector4 = 0;
   jets_AK5PF_area = 0;
   jets_AK5PF_corrFactorRaw = 0;
   jets_AK5PF_rawPt = 0;
   jets_AK5PF_mass = 0;
   jets_AK5PFclean_status = 0;
   jets_AK5PFclean_phi = 0;
   jets_AK5PFclean_pt = 0;
   jets_AK5PFclean_pz = 0;
   jets_AK5PFclean_px = 0;
   jets_AK5PFclean_py = 0;
   jets_AK5PFclean_eta = 0;
   jets_AK5PFclean_theta = 0;
   jets_AK5PFclean_et = 0;
   jets_AK5PFclean_energy = 0;
   jets_AK5PFclean_parton_Id = 0;
   jets_AK5PFclean_parton_motherId = 0;
   jets_AK5PFclean_parton_pt = 0;
   jets_AK5PFclean_parton_phi = 0;
   jets_AK5PFclean_parton_eta = 0;
   jets_AK5PFclean_parton_Energy = 0;
   jets_AK5PFclean_parton_mass = 0;
   jets_AK5PFclean_gen_et = 0;
   jets_AK5PFclean_gen_pt = 0;
   jets_AK5PFclean_gen_eta = 0;
   jets_AK5PFclean_gen_phi = 0;
   jets_AK5PFclean_gen_mass = 0;
   jets_AK5PFclean_gen_Energy = 0;
   jets_AK5PFclean_gen_Id = 0;
   jets_AK5PFclean_partonFlavour = 0;
   jets_AK5PFclean_btag_TC_highPur = 0;
   jets_AK5PFclean_btag_TC_highEff = 0;
   jets_AK5PFclean_btag_jetProb = 0;
   jets_AK5PFclean_btag_jetBProb = 0;
   jets_AK5PFclean_btag_softEle = 0;
   jets_AK5PFclean_btag_softMuon = 0;
   jets_AK5PFclean_btag_secVertexHighPur = 0;
   jets_AK5PFclean_btag_secVertexHighEff = 0;
   jets_AK5PFclean_btag_secVertexCombined = 0;
   jets_AK5PFclean_jetCharge = 0;
   jets_AK5PFclean_chgEmE = 0;
   jets_AK5PFclean_chgHadE = 0;
   jets_AK5PFclean_photonEnergy = 0;
   jets_AK5PFclean_chgMuE = 0;
   jets_AK5PFclean_chg_Mult = 0;
   jets_AK5PFclean_neutralEmE = 0;
   jets_AK5PFclean_neutralHadE = 0;
   jets_AK5PFclean_neutral_Mult = 0;
   jets_AK5PFclean_mu_Mult = 0;
   jets_AK5PFclean_emf = 0;
   jets_AK5PFclean_ehf = 0;
   jets_AK5PFclean_n60 = 0;
   jets_AK5PFclean_n90 = 0;
   jets_AK5PFclean_etaetaMoment = 0;
   jets_AK5PFclean_etaphiMoment = 0;
   jets_AK5PFclean_phiphiMoment = 0;
   jets_AK5PFclean_n90Hits = 0;
   jets_AK5PFclean_fHPD = 0;
   jets_AK5PFclean_fRBX = 0;
   jets_AK5PFclean_hitsInN90 = 0;
   jets_AK5PFclean_nECALTowers = 0;
   jets_AK5PFclean_nHCALTowers = 0;
   jets_AK5PFclean_fSubDetector1 = 0;
   jets_AK5PFclean_fSubDetector2 = 0;
   jets_AK5PFclean_fSubDetector3 = 0;
   jets_AK5PFclean_fSubDetector4 = 0;
   jets_AK5PFclean_area = 0;
   jets_AK5PFclean_corrFactorRaw = 0;
   jets_AK5PFclean_rawPt = 0;
   jets_AK5PFclean_mass = 0;
   mc_doc_id = 0;
   mc_doc_pt = 0;
   mc_doc_px = 0;
   mc_doc_py = 0;
   mc_doc_pz = 0;
   mc_doc_eta = 0;
   mc_doc_phi = 0;
   mc_doc_theta = 0;
   mc_doc_energy = 0;
   mc_doc_status = 0;
   mc_doc_charge = 0;
   mc_doc_mother_id = 0;
   mc_doc_grandmother_id = 0;
   mc_doc_ggrandmother_id = 0;
   mc_doc_mother_pt = 0;
   mc_doc_vertex_x = 0;
   mc_doc_vertex_y = 0;
   mc_doc_vertex_z = 0;
   mc_doc_mass = 0;
   mc_doc_numOfDaughters = 0;
   mc_doc_numOfMothers = 0;
   mc_electrons_id = 0;
   mc_electrons_pt = 0;
   mc_electrons_px = 0;
   mc_electrons_py = 0;
   mc_electrons_pz = 0;
   mc_electrons_eta = 0;
   mc_electrons_phi = 0;
   mc_electrons_theta = 0;
   mc_electrons_status = 0;
   mc_electrons_energy = 0;
   mc_electrons_charge = 0;
   mc_electrons_mother_id = 0;
   mc_electrons_mother_pt = 0;
   mc_electrons_grandmother_id = 0;
   mc_electrons_ggrandmother_id = 0;
   mc_electrons_vertex_x = 0;
   mc_electrons_vertex_y = 0;
   mc_electrons_vertex_z = 0;
   mc_electrons_mass = 0;
   mc_electrons_numOfDaughters = 0;
   mc_mus_id = 0;
   mc_mus_pt = 0;
   mc_mus_px = 0;
   mc_mus_py = 0;
   mc_mus_pz = 0;
   mc_mus_eta = 0;
   mc_mus_phi = 0;
   mc_mus_theta = 0;
   mc_mus_status = 0;
   mc_mus_energy = 0;
   mc_mus_charge = 0;
   mc_mus_mother_id = 0;
   mc_mus_mother_pt = 0;
   mc_mus_grandmother_id = 0;
   mc_mus_ggrandmother_id = 0;
   mc_mus_vertex_x = 0;
   mc_mus_vertex_y = 0;
   mc_mus_vertex_z = 0;
   mc_mus_mass = 0;
   mc_mus_numOfDaughters = 0;
   mc_nues_id = 0;
   mc_nues_pt = 0;
   mc_nues_px = 0;
   mc_nues_py = 0;
   mc_nues_pz = 0;
   mc_nues_eta = 0;
   mc_nues_phi = 0;
   mc_nues_theta = 0;
   mc_nues_status = 0;
   mc_nues_energy = 0;
   mc_nues_charge = 0;
   mc_nues_mother_id = 0;
   mc_nues_mother_pt = 0;
   mc_nues_grandmother_id = 0;
   mc_nues_ggrandmother_id = 0;
   mc_nues_vertex_x = 0;
   mc_nues_vertex_y = 0;
   mc_nues_vertex_z = 0;
   mc_nues_mass = 0;
   mc_nues_numOfDaughters = 0;
   mc_numus_id = 0;
   mc_numus_pt = 0;
   mc_numus_px = 0;
   mc_numus_py = 0;
   mc_numus_pz = 0;
   mc_numus_eta = 0;
   mc_numus_phi = 0;
   mc_numus_theta = 0;
   mc_numus_status = 0;
   mc_numus_energy = 0;
   mc_numus_charge = 0;
   mc_numus_mother_id = 0;
   mc_numus_mother_pt = 0;
   mc_numus_grandmother_id = 0;
   mc_numus_ggrandmother_id = 0;
   mc_numus_vertex_x = 0;
   mc_numus_vertex_y = 0;
   mc_numus_vertex_z = 0;
   mc_numus_mass = 0;
   mc_numus_numOfDaughters = 0;
   mc_nutaus_id = 0;
   mc_nutaus_pt = 0;
   mc_nutaus_px = 0;
   mc_nutaus_py = 0;
   mc_nutaus_pz = 0;
   mc_nutaus_eta = 0;
   mc_nutaus_phi = 0;
   mc_nutaus_theta = 0;
   mc_nutaus_status = 0;
   mc_nutaus_energy = 0;
   mc_nutaus_charge = 0;
   mc_nutaus_mother_id = 0;
   mc_nutaus_mother_pt = 0;
   mc_nutaus_grandmother_id = 0;
   mc_nutaus_ggrandmother_id = 0;
   mc_nutaus_vertex_x = 0;
   mc_nutaus_vertex_y = 0;
   mc_nutaus_vertex_z = 0;
   mc_nutaus_mass = 0;
   mc_nutaus_numOfDaughters = 0;
   mc_photons_id = 0;
   mc_photons_pt = 0;
   mc_photons_px = 0;
   mc_photons_py = 0;
   mc_photons_pz = 0;
   mc_photons_eta = 0;
   mc_photons_phi = 0;
   mc_photons_theta = 0;
   mc_photons_status = 0;
   mc_photons_energy = 0;
   mc_photons_charge = 0;
   mc_photons_mother_id = 0;
   mc_photons_mother_pt = 0;
   mc_photons_grandmother_id = 0;
   mc_photons_ggrandmother_id = 0;
   mc_photons_vertex_x = 0;
   mc_photons_vertex_y = 0;
   mc_photons_vertex_z = 0;
   mc_photons_mass = 0;
   mc_photons_numOfDaughters = 0;
   mc_taus_id = 0;
   mc_taus_pt = 0;
   mc_taus_px = 0;
   mc_taus_py = 0;
   mc_taus_pz = 0;
   mc_taus_eta = 0;
   mc_taus_phi = 0;
   mc_taus_theta = 0;
   mc_taus_status = 0;
   mc_taus_energy = 0;
   mc_taus_charge = 0;
   mc_taus_mother_id = 0;
   mc_taus_mother_pt = 0;
   mc_taus_grandmother_id = 0;
   mc_taus_ggrandmother_id = 0;
   mc_taus_vertex_x = 0;
   mc_taus_vertex_y = 0;
   mc_taus_vertex_z = 0;
   mc_taus_mass = 0;
   mc_taus_numOfDaughters = 0;
   mets_AK5_et = 0;
   mets_AK5_phi = 0;
   mets_AK5_ex = 0;
   mets_AK5_ey = 0;
   mets_AK5_gen_et = 0;
   mets_AK5_gen_phi = 0;
   mets_AK5_sign = 0;
   mets_AK5_sumEt = 0;
   mets_AK5_unCPhi = 0;
   mets_AK5_unCPt = 0;
   mus_energy = 0;
   mus_et = 0;
   mus_eta = 0;
   mus_phi = 0;
   mus_pt = 0;
   mus_px = 0;
   mus_py = 0;
   mus_pz = 0;
   mus_status = 0;
   mus_theta = 0;
   mus_gen_id = 0;
   mus_gen_phi = 0;
   mus_gen_pt = 0;
   mus_gen_pz = 0;
   mus_gen_px = 0;
   mus_gen_py = 0;
   mus_gen_eta = 0;
   mus_gen_theta = 0;
   mus_gen_et = 0;
   mus_gen_mother_id = 0;
   mus_gen_mother_phi = 0;
   mus_gen_mother_pt = 0;
   mus_gen_mother_pz = 0;
   mus_gen_mother_px = 0;
   mus_gen_mother_py = 0;
   mus_gen_mother_eta = 0;
   mus_gen_mother_theta = 0;
   mus_gen_mother_et = 0;
   mus_tkHits = 0;
   mus_cIso = 0;
   mus_tIso = 0;
   mus_ecalIso = 0;
   mus_hcalIso = 0;
   mus_ecalvetoDep = 0;
   mus_hcalvetoDep = 0;
   mus_calEnergyEm = 0;
   mus_calEnergyHad = 0;
   mus_calEnergyHo = 0;
   mus_calEnergyEmS9 = 0;
   mus_calEnergyHadS9 = 0;
   mus_calEnergyHoS9 = 0;
   mus_iso03_emVetoEt = 0;
   mus_iso03_hadVetoEt = 0;
   mus_iso03_sumPt = 0;
   mus_iso03_emEt = 0;
   mus_iso03_hadEt = 0;
   mus_iso03_hoEt = 0;
   mus_iso03_nTracks = 0;
   mus_iso05_sumPt = 0;
   mus_iso05_emEt = 0;
   mus_iso05_hadEt = 0;
   mus_iso05_hoEt = 0;
   mus_iso05_nTracks = 0;
   mus_pfIsolationR03_sumChargedHadronPt = 0;
   mus_pfIsolationR03_sumChargedParticlePt = 0;
   mus_pfIsolationR03_sumNeutralHadronEt = 0;
   mus_pfIsolationR03_sumNeutralHadronEtHighThreshold = 0;
   mus_pfIsolationR03_sumPhotonEt = 0;
   mus_pfIsolationR03_sumPhotonEtHighThreshold = 0;
   mus_pfIsolationR03_sumPUPt = 0;
   mus_pfIsolationR04_sumChargedHadronPt = 0;
   mus_pfIsolationR04_sumChargedParticlePt = 0;
   mus_pfIsolationR04_sumNeutralHadronEt = 0;
   mus_pfIsolationR04_sumNeutralHadronEtHighThreshold = 0;
   mus_pfIsolationR04_sumPhotonEt = 0;
   mus_pfIsolationR04_sumPhotonEtHighThreshold = 0;
   mus_pfIsolationR04_sumPUPt = 0;
   mus_charge = 0;
   mus_cm_chi2 = 0;
   mus_cm_ndof = 0;
   mus_cm_chg = 0;
   mus_cm_pt = 0;
   mus_cm_px = 0;
   mus_cm_py = 0;
   mus_cm_pz = 0;
   mus_cm_eta = 0;
   mus_cm_phi = 0;
   mus_cm_theta = 0;
   mus_cm_d0dum = 0;
   mus_cm_dz = 0;
   mus_cm_vx = 0;
   mus_cm_vy = 0;
   mus_cm_vz = 0;
   mus_cm_numvalhits = 0;
   mus_cm_numlosthits = 0;
   mus_cm_numvalMuonhits = 0;
   mus_cm_d0dumErr = 0;
   mus_cm_dzErr = 0;
   mus_cm_ptErr = 0;
   mus_cm_etaErr = 0;
   mus_cm_phiErr = 0;
   mus_tk_id = 0;
   mus_tk_chi2 = 0;
   mus_tk_ndof = 0;
   mus_tk_chg = 0;
   mus_tk_pt = 0;
   mus_tk_px = 0;
   mus_tk_py = 0;
   mus_tk_pz = 0;
   mus_tk_eta = 0;
   mus_tk_phi = 0;
   mus_tk_theta = 0;
   mus_tk_d0dum = 0;
   mus_tk_dz = 0;
   mus_tk_vx = 0;
   mus_tk_vy = 0;
   mus_tk_vz = 0;
   mus_tk_numvalhits = 0;
   mus_tk_numlosthits = 0;
   mus_tk_d0dumErr = 0;
   mus_tk_dzErr = 0;
   mus_tk_ptErr = 0;
   mus_tk_etaErr = 0;
   mus_tk_phiErr = 0;
   mus_tk_numvalPixelhits = 0;
   mus_tk_numpixelWthMeasr = 0;
   mus_stamu_chi2 = 0;
   mus_stamu_ndof = 0;
   mus_stamu_chg = 0;
   mus_stamu_pt = 0;
   mus_stamu_px = 0;
   mus_stamu_py = 0;
   mus_stamu_pz = 0;
   mus_stamu_eta = 0;
   mus_stamu_phi = 0;
   mus_stamu_theta = 0;
   mus_stamu_d0dum = 0;
   mus_stamu_dz = 0;
   mus_stamu_vx = 0;
   mus_stamu_vy = 0;
   mus_stamu_vz = 0;
   mus_stamu_numvalhits = 0;
   mus_stamu_numlosthits = 0;
   mus_stamu_d0dumErr = 0;
   mus_stamu_dzErr = 0;
   mus_stamu_ptErr = 0;
   mus_stamu_etaErr = 0;
   mus_stamu_phiErr = 0;
   mus_num_matches = 0;
   mus_isPFMuon = 0;
   mus_isTrackerMuon = 0;
   mus_isStandAloneMuon = 0;
   mus_isCaloMuon = 0;
   mus_isGlobalMuon = 0;
   mus_isElectron = 0;
   mus_isConvertedPhoton = 0;
   mus_isPhoton = 0;
   mus_id_All = 0;
   mus_id_AllGlobalMuons = 0;
   mus_id_AllStandAloneMuons = 0;
   mus_id_AllTrackerMuons = 0;
   mus_id_TrackerMuonArbitrated = 0;
   mus_id_AllArbitrated = 0;
   mus_id_GlobalMuonPromptTight = 0;
   mus_id_TMLastStationLoose = 0;
   mus_id_TMLastStationTight = 0;
   mus_id_TM2DCompatibilityLoose = 0;
   mus_id_TM2DCompatibilityTight = 0;
   mus_id_TMOneStationLoose = 0;
   mus_id_TMOneStationTight = 0;
   mus_id_TMLastStationOptimizedLowPtLoose = 0;
   mus_id_TMLastStationOptimizedLowPtTight = 0;
   mus_tk_LayersWithMeasurement = 0;
   mus_tk_PixelLayersWithMeasurement = 0;
   mus_tk_ValidStripLayersWithMonoAndStereoHit = 0;
   mus_tk_LayersWithoutMeasurement = 0;
   mus_tk_ExpectedHitsInner = 0;
   mus_tk_ExpectedHitsOuter = 0;
   mus_cm_LayersWithMeasurement = 0;
   mus_cm_PixelLayersWithMeasurement = 0;
   mus_cm_ValidStripLayersWithMonoAndStereoHit = 0;
   mus_cm_LayersWithoutMeasurement = 0;
   mus_cm_ExpectedHitsInner = 0;
   mus_cm_ExpectedHitsOuter = 0;
   mus_picky_LayersWithMeasurement = 0;
   mus_picky_PixelLayersWithMeasurement = 0;
   mus_picky_ValidStripLayersWithMonoAndStereoHit = 0;
   mus_picky_LayersWithoutMeasurement = 0;
   mus_picky_ExpectedHitsInner = 0;
   mus_picky_ExpectedHitsOuter = 0;
   mus_tpfms_LayersWithMeasurement = 0;
   mus_tpfms_PixelLayersWithMeasurement = 0;
   mus_tpfms_ValidStripLayersWithMonoAndStereoHit = 0;
   mus_tpfms_LayersWithoutMeasurement = 0;
   mus_tpfms_ExpectedHitsInner = 0;
   mus_tpfms_ExpectedHitsOuter = 0;
   mus_picky_id = 0;
   mus_picky_chi2 = 0;
   mus_picky_ndof = 0;
   mus_picky_chg = 0;
   mus_picky_pt = 0;
   mus_picky_px = 0;
   mus_picky_py = 0;
   mus_picky_pz = 0;
   mus_picky_eta = 0;
   mus_picky_phi = 0;
   mus_picky_theta = 0;
   mus_picky_d0dum = 0;
   mus_picky_dz = 0;
   mus_picky_vx = 0;
   mus_picky_vy = 0;
   mus_picky_vz = 0;
   mus_picky_numvalhits = 0;
   mus_picky_numlosthits = 0;
   mus_picky_d0dumErr = 0;
   mus_picky_dzErr = 0;
   mus_picky_ptErr = 0;
   mus_picky_etaErr = 0;
   mus_picky_phiErr = 0;
   mus_picky_numvalPixelhits = 0;
   mus_tpfms_id = 0;
   mus_tpfms_chi2 = 0;
   mus_tpfms_ndof = 0;
   mus_tpfms_chg = 0;
   mus_tpfms_pt = 0;
   mus_tpfms_px = 0;
   mus_tpfms_py = 0;
   mus_tpfms_pz = 0;
   mus_tpfms_eta = 0;
   mus_tpfms_phi = 0;
   mus_tpfms_theta = 0;
   mus_tpfms_d0dum = 0;
   mus_tpfms_dz = 0;
   mus_tpfms_vx = 0;
   mus_tpfms_vy = 0;
   mus_tpfms_vz = 0;
   mus_tpfms_numvalhits = 0;
   mus_tpfms_numlosthits = 0;
   mus_tpfms_d0dumErr = 0;
   mus_tpfms_dzErr = 0;
   mus_tpfms_ptErr = 0;
   mus_tpfms_etaErr = 0;
   mus_tpfms_phiErr = 0;
   mus_tpfms_numvalPixelhits = 0;
   mus_dB = 0;
   mus_numberOfMatchedStations = 0;
   pfTypeINoXYCorrmets_et = 0;
   pfTypeINoXYCorrmets_phi = 0;
   pfTypeINoXYCorrmets_ex = 0;
   pfTypeINoXYCorrmets_ey = 0;
   pfTypeINoXYCorrmets_gen_et = 0;
   pfTypeINoXYCorrmets_gen_phi = 0;
   pfTypeINoXYCorrmets_sign = 0;
   pfTypeINoXYCorrmets_sumEt = 0;
   pfTypeINoXYCorrmets_unCPhi = 0;
   pfTypeINoXYCorrmets_unCPt = 0;
   pfTypeIType0mets_et = 0;
   pfTypeIType0mets_phi = 0;
   pfTypeIType0mets_ex = 0;
   pfTypeIType0mets_ey = 0;
   pfTypeIType0mets_gen_et = 0;
   pfTypeIType0mets_gen_phi = 0;
   pfTypeIType0mets_sign = 0;
   pfTypeIType0mets_sumEt = 0;
   pfTypeIType0mets_unCPhi = 0;
   pfTypeIType0mets_unCPt = 0;
   pfTypeImets_et = 0;
   pfTypeImets_phi = 0;
   pfTypeImets_ex = 0;
   pfTypeImets_ey = 0;
   pfTypeImets_gen_et = 0;
   pfTypeImets_gen_phi = 0;
   pfTypeImets_sign = 0;
   pfTypeImets_sumEt = 0;
   pfTypeImets_unCPhi = 0;
   pfTypeImets_unCPt = 0;
   pf_els_energy = 0;
   pf_els_et = 0;
   pf_els_eta = 0;
   pf_els_phi = 0;
   pf_els_pt = 0;
   pf_els_px = 0;
   pf_els_py = 0;
   pf_els_pz = 0;
   pf_els_status = 0;
   pf_els_theta = 0;
   pf_els_gen_id = 0;
   pf_els_gen_phi = 0;
   pf_els_gen_pt = 0;
   pf_els_gen_pz = 0;
   pf_els_gen_px = 0;
   pf_els_gen_py = 0;
   pf_els_gen_eta = 0;
   pf_els_gen_theta = 0;
   pf_els_gen_et = 0;
   pf_els_gen_mother_id = 0;
   pf_els_gen_mother_phi = 0;
   pf_els_gen_mother_pt = 0;
   pf_els_gen_mother_pz = 0;
   pf_els_gen_mother_px = 0;
   pf_els_gen_mother_py = 0;
   pf_els_gen_mother_eta = 0;
   pf_els_gen_mother_theta = 0;
   pf_els_gen_mother_et = 0;
   pf_els_tightId = 0;
   pf_els_looseId = 0;
   pf_els_robustTightId = 0;
   pf_els_robustLooseId = 0;
   pf_els_robustHighEnergyId = 0;
   pf_els_simpleEleId95relIso = 0;
   pf_els_simpleEleId90relIso = 0;
   pf_els_simpleEleId85relIso = 0;
   pf_els_simpleEleId80relIso = 0;
   pf_els_simpleEleId70relIso = 0;
   pf_els_simpleEleId60relIso = 0;
   pf_els_simpleEleId95cIso = 0;
   pf_els_simpleEleId90cIso = 0;
   pf_els_simpleEleId85cIso = 0;
   pf_els_simpleEleId80cIso = 0;
   pf_els_simpleEleId70cIso = 0;
   pf_els_simpleEleId60cIso = 0;
   pf_els_cIso = 0;
   pf_els_tIso = 0;
   pf_els_ecalIso = 0;
   pf_els_hcalIso = 0;
   pf_els_chargedHadronIso = 0;
   pf_els_photonIso = 0;
   pf_els_neutralHadronIso = 0;
   pf_els_chi2 = 0;
   pf_els_charge = 0;
   pf_els_caloEnergy = 0;
   pf_els_hadOverEm = 0;
   pf_els_hcalOverEcalBc = 0;
   pf_els_eOverPIn = 0;
   pf_els_eSeedOverPOut = 0;
   pf_els_sigmaEtaEta = 0;
   pf_els_sigmaIEtaIEta = 0;
   pf_els_scEnergy = 0;
   pf_els_scRawEnergy = 0;
   pf_els_scSeedEnergy = 0;
   pf_els_scEta = 0;
   pf_els_scPhi = 0;
   pf_els_scEtaWidth = 0;
   pf_els_scPhiWidth = 0;
   pf_els_scE1x5 = 0;
   pf_els_scE2x5Max = 0;
   pf_els_scE5x5 = 0;
   pf_els_isEB = 0;
   pf_els_isEE = 0;
   pf_els_dEtaIn = 0;
   pf_els_dPhiIn = 0;
   pf_els_dEtaOut = 0;
   pf_els_dPhiOut = 0;
   pf_els_numvalhits = 0;
   pf_els_numlosthits = 0;
   pf_els_basicClustersSize = 0;
   pf_els_tk_pz = 0;
   pf_els_tk_pt = 0;
   pf_els_tk_phi = 0;
   pf_els_tk_eta = 0;
   pf_els_d0dum = 0;
   pf_els_dz = 0;
   pf_els_vx = 0;
   pf_els_vy = 0;
   pf_els_vz = 0;
   pf_els_ndof = 0;
   pf_els_ptError = 0;
   pf_els_d0dumError = 0;
   pf_els_dzError = 0;
   pf_els_etaError = 0;
   pf_els_phiError = 0;
   pf_els_tk_charge = 0;
   pf_els_core_ecalDrivenSeed = 0;
   pf_els_n_inner_layer = 0;
   pf_els_n_outer_layer = 0;
   pf_els_ctf_tk_id = 0;
   pf_els_ctf_tk_charge = 0;
   pf_els_ctf_tk_eta = 0;
   pf_els_ctf_tk_phi = 0;
   pf_els_fbrem = 0;
   pf_els_shFracInnerHits = 0;
   pf_els_dr03EcalRecHitSumEt = 0;
   pf_els_dr03HcalTowerSumEt = 0;
   pf_els_dr03HcalDepth1TowerSumEt = 0;
   pf_els_dr03HcalDepth2TowerSumEt = 0;
   pf_els_dr03TkSumPt = 0;
   pf_els_dr04EcalRecHitSumEt = 0;
   pf_els_dr04HcalTowerSumEt = 0;
   pf_els_dr04HcalDepth1TowerSumEt = 0;
   pf_els_dr04HcalDepth2TowerSumEt = 0;
   pf_els_dr04TkSumPt = 0;
   pf_els_cpx = 0;
   pf_els_cpy = 0;
   pf_els_cpz = 0;
   pf_els_vpx = 0;
   pf_els_vpy = 0;
   pf_els_vpz = 0;
   pf_els_cx = 0;
   pf_els_cy = 0;
   pf_els_cz = 0;
   pf_els_PATpassConversionVeto = 0;
   pf_mus_energy = 0;
   pf_mus_et = 0;
   pf_mus_eta = 0;
   pf_mus_phi = 0;
   pf_mus_pt = 0;
   pf_mus_px = 0;
   pf_mus_py = 0;
   pf_mus_pz = 0;
   pf_mus_status = 0;
   pf_mus_theta = 0;
   pf_mus_gen_id = 0;
   pf_mus_gen_phi = 0;
   pf_mus_gen_pt = 0;
   pf_mus_gen_pz = 0;
   pf_mus_gen_px = 0;
   pf_mus_gen_py = 0;
   pf_mus_gen_eta = 0;
   pf_mus_gen_theta = 0;
   pf_mus_gen_et = 0;
   pf_mus_gen_mother_id = 0;
   pf_mus_gen_mother_phi = 0;
   pf_mus_gen_mother_pt = 0;
   pf_mus_gen_mother_pz = 0;
   pf_mus_gen_mother_px = 0;
   pf_mus_gen_mother_py = 0;
   pf_mus_gen_mother_eta = 0;
   pf_mus_gen_mother_theta = 0;
   pf_mus_gen_mother_et = 0;
   pf_mus_tkHits = 0;
   pf_mus_cIso = 0;
   pf_mus_tIso = 0;
   pf_mus_ecalIso = 0;
   pf_mus_hcalIso = 0;
   pf_mus_iso03_emVetoEt = 0;
   pf_mus_iso03_hadVetoEt = 0;
   pf_mus_calEnergyEm = 0;
   pf_mus_calEnergyHad = 0;
   pf_mus_calEnergyHo = 0;
   pf_mus_calEnergyEmS9 = 0;
   pf_mus_calEnergyHadS9 = 0;
   pf_mus_calEnergyHoS9 = 0;
   pf_mus_iso03_sumPt = 0;
   pf_mus_iso03_emEt = 0;
   pf_mus_iso03_hadEt = 0;
   pf_mus_iso03_hoEt = 0;
   pf_mus_iso03_nTracks = 0;
   pf_mus_iso05_sumPt = 0;
   pf_mus_iso05_emEt = 0;
   pf_mus_iso05_hadEt = 0;
   pf_mus_iso05_hoEt = 0;
   pf_mus_iso05_nTracks = 0;
   pf_mus_neutralHadronIso = 0;
   pf_mus_chargedHadronIso = 0;
   pf_mus_photonIso = 0;
   pf_mus_charge = 0;
   pf_mus_cm_chi2 = 0;
   pf_mus_cm_ndof = 0;
   pf_mus_cm_chg = 0;
   pf_mus_cm_pt = 0;
   pf_mus_cm_px = 0;
   pf_mus_cm_py = 0;
   pf_mus_cm_pz = 0;
   pf_mus_cm_eta = 0;
   pf_mus_cm_phi = 0;
   pf_mus_cm_theta = 0;
   pf_mus_cm_d0dum = 0;
   pf_mus_cm_dz = 0;
   pf_mus_cm_vx = 0;
   pf_mus_cm_vy = 0;
   pf_mus_cm_vz = 0;
   pf_mus_cm_numvalhits = 0;
   pf_mus_cm_numlosthits = 0;
   pf_mus_cm_numvalMuonhits = 0;
   pf_mus_cm_d0dumErr = 0;
   pf_mus_cm_dzErr = 0;
   pf_mus_cm_ptErr = 0;
   pf_mus_cm_etaErr = 0;
   pf_mus_cm_phiErr = 0;
   pf_mus_tk_id = 0;
   pf_mus_tk_chi2 = 0;
   pf_mus_tk_ndof = 0;
   pf_mus_tk_chg = 0;
   pf_mus_tk_pt = 0;
   pf_mus_tk_px = 0;
   pf_mus_tk_py = 0;
   pf_mus_tk_pz = 0;
   pf_mus_tk_eta = 0;
   pf_mus_tk_phi = 0;
   pf_mus_tk_theta = 0;
   pf_mus_tk_d0dum = 0;
   pf_mus_tk_dz = 0;
   pf_mus_tk_vx = 0;
   pf_mus_tk_vy = 0;
   pf_mus_tk_vz = 0;
   pf_mus_tk_numvalhits = 0;
   pf_mus_tk_numlosthits = 0;
   pf_mus_tk_d0dumErr = 0;
   pf_mus_tk_dzErr = 0;
   pf_mus_tk_ptErr = 0;
   pf_mus_tk_etaErr = 0;
   pf_mus_tk_phiErr = 0;
   pf_mus_tk_numvalPixelhits = 0;
   pf_mus_tk_numpixelWthMeasr = 0;
   pf_mus_stamu_chi2 = 0;
   pf_mus_stamu_ndof = 0;
   pf_mus_stamu_chg = 0;
   pf_mus_stamu_pt = 0;
   pf_mus_stamu_px = 0;
   pf_mus_stamu_py = 0;
   pf_mus_stamu_pz = 0;
   pf_mus_stamu_eta = 0;
   pf_mus_stamu_phi = 0;
   pf_mus_stamu_theta = 0;
   pf_mus_stamu_d0dum = 0;
   pf_mus_stamu_dz = 0;
   pf_mus_stamu_vx = 0;
   pf_mus_stamu_vy = 0;
   pf_mus_stamu_vz = 0;
   pf_mus_stamu_numvalhits = 0;
   pf_mus_stamu_numlosthits = 0;
   pf_mus_stamu_d0dumErr = 0;
   pf_mus_stamu_dzErr = 0;
   pf_mus_stamu_ptErr = 0;
   pf_mus_stamu_etaErr = 0;
   pf_mus_stamu_phiErr = 0;
   pf_mus_num_matches = 0;
   pf_mus_isTrackerMuon = 0;
   pf_mus_isStandAloneMuon = 0;
   pf_mus_isCaloMuon = 0;
   pf_mus_isGlobalMuon = 0;
   pf_mus_isElectron = 0;
   pf_mus_isConvertedPhoton = 0;
   pf_mus_isPhoton = 0;
   pf_mus_id_All = 0;
   pf_mus_id_AllGlobalMuons = 0;
   pf_mus_id_AllStandAloneMuons = 0;
   pf_mus_id_AllTrackerMuons = 0;
   pf_mus_id_TrackerMuonArbitrated = 0;
   pf_mus_id_AllArbitrated = 0;
   pf_mus_id_GlobalMuonPromptTight = 0;
   pf_mus_id_TMLastStationLoose = 0;
   pf_mus_id_TMLastStationTight = 0;
   pf_mus_id_TM2DCompatibilityLoose = 0;
   pf_mus_id_TM2DCompatibilityTight = 0;
   pf_mus_id_TMOneStationLoose = 0;
   pf_mus_id_TMOneStationTight = 0;
   pf_mus_id_TMLastStationOptimizedLowPtLoose = 0;
   pf_mus_id_TMLastStationOptimizedLowPtTight = 0;
   pf_mus_tk_LayersWithMeasurement = 0;
   pf_mus_tk_PixelLayersWithMeasurement = 0;
   pf_mus_tk_ValidStripLayersWithMonoAndStereoHit = 0;
   pf_mus_tk_LayersWithoutMeasurement = 0;
   pf_mus_tk_ExpectedHitsInner = 0;
   pf_mus_tk_ExpectedHitsOuter = 0;
   pf_mus_cm_LayersWithMeasurement = 0;
   pf_mus_cm_PixelLayersWithMeasurement = 0;
   pf_mus_cm_ValidStripLayersWithMonoAndStereoHit = 0;
   pf_mus_cm_LayersWithoutMeasurement = 0;
   pf_mus_cm_ExpectedHitsInner = 0;
   pf_mus_cm_ExpectedHitsOuter = 0;
   pf_mus_picky_LayersWithMeasurement = 0;
   pf_mus_picky_PixelLayersWithMeasurement = 0;
   pf_mus_picky_ValidStripLayersWithMonoAndStereoHit = 0;
   pf_mus_picky_LayersWithoutMeasurement = 0;
   pf_mus_picky_ExpectedHitsInner = 0;
   pf_mus_picky_ExpectedHitsOuter = 0;
   pf_mus_tpfms_LayersWithMeasurement = 0;
   pf_mus_tpfms_PixelLayersWithMeasurement = 0;
   pf_mus_tpfms_ValidStripLayersWithMonoAndStereoHit = 0;
   pf_mus_tpfms_LayersWithoutMeasurement = 0;
   pf_mus_tpfms_ExpectedHitsInner = 0;
   pf_mus_tpfms_ExpectedHitsOuter = 0;
   pf_mus_picky_id = 0;
   pf_mus_picky_chi2 = 0;
   pf_mus_picky_ndof = 0;
   pf_mus_picky_chg = 0;
   pf_mus_picky_pt = 0;
   pf_mus_picky_px = 0;
   pf_mus_picky_py = 0;
   pf_mus_picky_pz = 0;
   pf_mus_picky_eta = 0;
   pf_mus_picky_phi = 0;
   pf_mus_picky_theta = 0;
   pf_mus_picky_d0dum = 0;
   pf_mus_picky_dz = 0;
   pf_mus_picky_vx = 0;
   pf_mus_picky_vy = 0;
   pf_mus_picky_vz = 0;
   pf_mus_picky_numvalhits = 0;
   pf_mus_picky_numlosthits = 0;
   pf_mus_picky_d0dumErr = 0;
   pf_mus_picky_dzErr = 0;
   pf_mus_picky_ptErr = 0;
   pf_mus_picky_etaErr = 0;
   pf_mus_picky_phiErr = 0;
   pf_mus_picky_numvalPixelhits = 0;
   pf_mus_tpfms_id = 0;
   pf_mus_tpfms_chi2 = 0;
   pf_mus_tpfms_ndof = 0;
   pf_mus_tpfms_chg = 0;
   pf_mus_tpfms_pt = 0;
   pf_mus_tpfms_px = 0;
   pf_mus_tpfms_py = 0;
   pf_mus_tpfms_pz = 0;
   pf_mus_tpfms_eta = 0;
   pf_mus_tpfms_phi = 0;
   pf_mus_tpfms_theta = 0;
   pf_mus_tpfms_d0dum = 0;
   pf_mus_tpfms_dz = 0;
   pf_mus_tpfms_vx = 0;
   pf_mus_tpfms_vy = 0;
   pf_mus_tpfms_vz = 0;
   pf_mus_tpfms_numvalhits = 0;
   pf_mus_tpfms_numlosthits = 0;
   pf_mus_tpfms_d0dumErr = 0;
   pf_mus_tpfms_dzErr = 0;
   pf_mus_tpfms_ptErr = 0;
   pf_mus_tpfms_etaErr = 0;
   pf_mus_tpfms_phiErr = 0;
   pf_mus_tpfms_numvalPixelhits = 0;
   pf_mus_pfIsolationR03_sumChargedHadronPt = 0;
   pf_mus_pfIsolationR03_sumChargedParticlePt = 0;
   pf_mus_pfIsolationR03_sumNeutralHadronEt = 0;
   pf_mus_pfIsolationR03_sumNeutralHadronEtHighThreshold = 0;
   pf_mus_pfIsolationR03_sumPhotonEt = 0;
   pf_mus_pfIsolationR03_sumPhotonEtHighThreshold = 0;
   pf_mus_pfIsolationR03_sumPUPt = 0;
   pf_mus_pfIsolationR04_sumChargedHadronPt = 0;
   pf_mus_pfIsolationR04_sumChargedParticlePt = 0;
   pf_mus_pfIsolationR04_sumNeutralHadronEt = 0;
   pf_mus_pfIsolationR04_sumNeutralHadronEtHighThreshold = 0;
   pf_mus_pfIsolationR04_sumPhotonEt = 0;
   pf_mus_pfIsolationR04_sumPhotonEtHighThreshold = 0;
   pf_mus_pfIsolationR04_sumPUPt = 0;
   pf_mus_dB = 0;
   pf_mus_numberOfMatchedStations = 0;
   pf_mus_isPFMuon = 0;
   pfcand_pdgId = 0;
   pfcand_particleId = 0;
   pfcand_pt = 0;
   pfcand_pz = 0;
   pfcand_px = 0;
   pfcand_py = 0;
   pfcand_eta = 0;
   pfcand_phi = 0;
   pfcand_theta = 0;
   pfcand_energy = 0;
   pfcand_charge = 0;
   pfmets_et = 0;
   pfmets_phi = 0;
   pfmets_ex = 0;
   pfmets_ey = 0;
   pfmets_gen_et = 0;
   pfmets_gen_phi = 0;
   pfmets_sign = 0;
   pfmets_sumEt = 0;
   pfmets_unCPhi = 0;
   pfmets_unCPt = 0;
   photons_energy = 0;
   photons_et = 0;
   photons_eta = 0;
   photons_phi = 0;
   photons_pt = 0;
   photons_px = 0;
   photons_py = 0;
   photons_pz = 0;
   photons_status = 0;
   photons_theta = 0;
   photons_hadOverEM = 0;
   photons_scEnergy = 0;
   photons_scRawEnergy = 0;
   photons_scEta = 0;
   photons_scPhi = 0;
   photons_scEtaWidth = 0;
   photons_scPhiWidth = 0;
   photons_tIso = 0;
   photons_ecalIso = 0;
   photons_hcalIso = 0;
   photons_isoEcalRecHitDR04 = 0;
   photons_isoHcalRecHitDR04 = 0;
   photons_isoSolidTrkConeDR04 = 0;
   photons_isoHollowTrkConeDR04 = 0;
   photons_nTrkSolidConeDR04 = 0;
   photons_nTrkHollowConeDR04 = 0;
   photons_isoEcalRecHitDR03 = 0;
   photons_isoHcalRecHitDR03 = 0;
   photons_isoSolidTrkConeDR03 = 0;
   photons_isoHollowTrkConeDR03 = 0;
   photons_nTrkSolidConeDR03 = 0;
   photons_nTrkHollowConeDR03 = 0;
   photons_isAlsoElectron = 0;
   photons_hasPixelSeed = 0;
   photons_isConverted = 0;
   photons_isEBGap = 0;
   photons_isEEGap = 0;
   photons_isEBEEGap = 0;
   photons_isEBPho = 0;
   photons_isEEPho = 0;
   photons_isLoosePhoton = 0;
   photons_isTightPhoton = 0;
   photons_maxEnergyXtal = 0;
   photons_e1x5 = 0;
   photons_e2x5 = 0;
   photons_e3x3 = 0;
   photons_e5x5 = 0;
   photons_sigmaEtaEta = 0;
   photons_sigmaIetaIeta = 0;
   photons_r9 = 0;
   photons_gen_et = 0;
   photons_gen_eta = 0;
   photons_gen_phi = 0;
   photons_gen_id = 0;
   pv_x = 0;
   pv_y = 0;
   pv_z = 0;
   pv_xErr = 0;
   pv_yErr = 0;
   pv_zErr = 0;
   pv_chi2 = 0;
   pv_ndof = 0;
   pv_isFake = 0;
   pv_isValid = 0;
   pv_tracksSize = 0;
   taus_status = 0;
   taus_phi = 0;
   taus_pt = 0;
   taus_pz = 0;
   taus_px = 0;
   taus_py = 0;
   taus_eta = 0;
   taus_theta = 0;
   taus_et = 0;
   taus_energy = 0;
   taus_charge = 0;
   taus_emf = 0;
   taus_hcalTotOverPLead = 0;
   taus_hcalMaxOverPLead = 0;
   taus_hcal3x3OverPLead = 0;
   taus_ecalStripSumEOverPLead = 0;
   taus_elecPreIdOutput = 0;
   taus_elecPreIdDecision = 0;
   taus_leadPFChargedHadrCand_pt = 0;
   taus_leadPFChargedHadrCand_charge = 0;
   taus_leadPFChargedHadrCand_eta = 0;
   taus_leadPFChargedHadrCand_ECAL_eta = 0;
   taus_leadPFChargedHadrCand_phi = 0;
   taus_isoPFGammaCandsEtSum = 0;
   taus_isoPFChargedHadrCandsPtSum = 0;
   taus_leadingTrackFinding = 0;
   taus_leadingTrackPtCut = 0;
   taus_trackIsolation = 0;
   taus_ecalIsolation = 0;
   taus_byIsolation = 0;
   taus_againstElectron = 0;
   taus_againstMuon = 0;
   taus_taNC_quarter = 0;
   taus_taNC_one = 0;
   taus_taNC_half = 0;
   taus_taNC_tenth = 0;
   taus_taNC = 0;
   taus_byIsoUsingLeadingPi = 0;
   taus_tkIsoUsingLeadingPi = 0;
   taus_ecalIsoUsingLeadingPi = 0;
   taus_againstElectronLoose = 0;
   taus_againstElectronMedium = 0;
   taus_againstElectronTight = 0;
   taus_againstElectronMVA = 0;
   taus_againstMuonLoose = 0;
   taus_againstMuonMedium = 0;
   taus_againstMuonTight = 0;
   taus_decayModeFinding = 0;
   taus_byVLooseIsolation = 0;
   taus_byLooseIsolation = 0;
   taus_byMediumIsolation = 0;
   taus_byTightIsolation = 0;
   taus_byVLooseIsolationDeltaBetaCorr = 0;
   taus_byLooseIsolationDeltaBetaCorr = 0;
   taus_byMediumIsolationDeltaBetaCorr = 0;
   taus_byTightIsolationDeltaBetaCorr = 0;
   taus_signalPFChargedHadrCandsSize = 0;
   taus_muDecision = 0;
   taus_Nprongs = 0;
   tcmets_et = 0;
   tcmets_phi = 0;
   tcmets_ex = 0;
   tcmets_ey = 0;
   tcmets_sumEt = 0;
   tracks_chi2 = 0;
   tracks_ndof = 0;
   tracks_chg = 0;
   tracks_pt = 0;
   tracks_px = 0;
   tracks_py = 0;
   tracks_pz = 0;
   tracks_eta = 0;
   tracks_phi = 0;
   tracks_d0dum = 0;
   tracks_dz = 0;
   tracks_vx = 0;
   tracks_vy = 0;
   tracks_vz = 0;
   tracks_numvalhits = 0;
   tracks_numlosthits = 0;
   tracks_d0dumErr = 0;
   tracks_dzErr = 0;
   tracks_ptErr = 0;
   tracks_etaErr = 0;
   tracks_phiErr = 0;
   tracks_highPurity = 0;
   model_params = 0;

   fChain->SetBranchAddress("NbeamSpot", &NbeamSpot, &b_NbeamSpot);
   fChain->SetBranchAddress("beamSpot_x", &beamSpot_x, &b_beamSpot_x);
   fChain->SetBranchAddress("beamSpot_y", &beamSpot_y, &b_beamSpot_y);
   fChain->SetBranchAddress("beamSpot_z", &beamSpot_z, &b_beamSpot_z);
   fChain->SetBranchAddress("beamSpot_x0Error", &beamSpot_x0Error, &b_beamSpot_x0Error);
   fChain->SetBranchAddress("beamSpot_y0Error", &beamSpot_y0Error, &b_beamSpot_y0Error);
   fChain->SetBranchAddress("beamSpot_z0Error", &beamSpot_z0Error, &b_beamSpot_z0Error);
   fChain->SetBranchAddress("beamSpot_sigmaZ", &beamSpot_sigmaZ, &b_beamSpot_sigmaZ);
   fChain->SetBranchAddress("beamSpot_sigmaZ0Error", &beamSpot_sigmaZ0Error, &b_beamSpot_sigmaZ0Error);
   fChain->SetBranchAddress("beamSpot_dxdz", &beamSpot_dxdz, &b_beamSpot_dxdz);
   fChain->SetBranchAddress("beamSpot_dxdzError", &beamSpot_dxdzError, &b_beamSpot_dxdzError);
   fChain->SetBranchAddress("beamSpot_dydz", &beamSpot_dydz, &b_beamSpot_dydz);
   fChain->SetBranchAddress("beamSpot_dydzError", &beamSpot_dydzError, &b_beamSpot_dydzError);
   fChain->SetBranchAddress("beamSpot_beamWidthX", &beamSpot_beamWidthX, &b_beamSpot_beamWidthX);
   fChain->SetBranchAddress("beamSpot_beamWidthY", &beamSpot_beamWidthY, &b_beamSpot_beamWidthY);
   fChain->SetBranchAddress("beamSpot_beamWidthXError", &beamSpot_beamWidthXError, &b_beamSpot_beamWidthXError);
   fChain->SetBranchAddress("beamSpot_beamWidthYError", &beamSpot_beamWidthYError, &b_beamSpot_beamWidthYError);
   fChain->SetBranchAddress("Nels", &Nels, &b_Nels);
   fChain->SetBranchAddress("els_energy", &els_energy, &b_els_energy);
   fChain->SetBranchAddress("els_et", &els_et, &b_els_et);
   fChain->SetBranchAddress("els_eta", &els_eta, &b_els_eta);
   fChain->SetBranchAddress("els_phi", &els_phi, &b_els_phi);
   fChain->SetBranchAddress("els_pt", &els_pt, &b_els_pt);
   fChain->SetBranchAddress("els_px", &els_px, &b_els_px);
   fChain->SetBranchAddress("els_py", &els_py, &b_els_py);
   fChain->SetBranchAddress("els_pz", &els_pz, &b_els_pz);
   fChain->SetBranchAddress("els_status", &els_status, &b_els_status);
   fChain->SetBranchAddress("els_theta", &els_theta, &b_els_theta);
   fChain->SetBranchAddress("els_gen_id", &els_gen_id, &b_els_gen_id);
   fChain->SetBranchAddress("els_gen_phi", &els_gen_phi, &b_els_gen_phi);
   fChain->SetBranchAddress("els_gen_pt", &els_gen_pt, &b_els_gen_pt);
   fChain->SetBranchAddress("els_gen_pz", &els_gen_pz, &b_els_gen_pz);
   fChain->SetBranchAddress("els_gen_px", &els_gen_px, &b_els_gen_px);
   fChain->SetBranchAddress("els_gen_py", &els_gen_py, &b_els_gen_py);
   fChain->SetBranchAddress("els_gen_eta", &els_gen_eta, &b_els_gen_eta);
   fChain->SetBranchAddress("els_gen_theta", &els_gen_theta, &b_els_gen_theta);
   fChain->SetBranchAddress("els_gen_et", &els_gen_et, &b_els_gen_et);
   fChain->SetBranchAddress("els_gen_mother_id", &els_gen_mother_id, &b_els_gen_mother_id);
   fChain->SetBranchAddress("els_gen_mother_phi", &els_gen_mother_phi, &b_els_gen_mother_phi);
   fChain->SetBranchAddress("els_gen_mother_pt", &els_gen_mother_pt, &b_els_gen_mother_pt);
   fChain->SetBranchAddress("els_gen_mother_pz", &els_gen_mother_pz, &b_els_gen_mother_pz);
   fChain->SetBranchAddress("els_gen_mother_px", &els_gen_mother_px, &b_els_gen_mother_px);
   fChain->SetBranchAddress("els_gen_mother_py", &els_gen_mother_py, &b_els_gen_mother_py);
   fChain->SetBranchAddress("els_gen_mother_eta", &els_gen_mother_eta, &b_els_gen_mother_eta);
   fChain->SetBranchAddress("els_gen_mother_theta", &els_gen_mother_theta, &b_els_gen_mother_theta);
   fChain->SetBranchAddress("els_gen_mother_et", &els_gen_mother_et, &b_els_gen_mother_et);
   fChain->SetBranchAddress("els_tightId", &els_tightId, &b_els_tightId);
   fChain->SetBranchAddress("els_looseId", &els_looseId, &b_els_looseId);
   fChain->SetBranchAddress("els_robustTightId", &els_robustTightId, &b_els_robustTightId);
   fChain->SetBranchAddress("els_robustLooseId", &els_robustLooseId, &b_els_robustLooseId);
   fChain->SetBranchAddress("els_robustHighEnergyId", &els_robustHighEnergyId, &b_els_robustHighEnergyId);
   fChain->SetBranchAddress("els_simpleEleId95relIso", &els_simpleEleId95relIso, &b_els_simpleEleId95relIso);
   fChain->SetBranchAddress("els_simpleEleId90relIso", &els_simpleEleId90relIso, &b_els_simpleEleId90relIso);
   fChain->SetBranchAddress("els_simpleEleId85relIso", &els_simpleEleId85relIso, &b_els_simpleEleId85relIso);
   fChain->SetBranchAddress("els_simpleEleId80relIso", &els_simpleEleId80relIso, &b_els_simpleEleId80relIso);
   fChain->SetBranchAddress("els_simpleEleId70relIso", &els_simpleEleId70relIso, &b_els_simpleEleId70relIso);
   fChain->SetBranchAddress("els_simpleEleId60relIso", &els_simpleEleId60relIso, &b_els_simpleEleId60relIso);
   fChain->SetBranchAddress("els_simpleEleId95cIso", &els_simpleEleId95cIso, &b_els_simpleEleId95cIso);
   fChain->SetBranchAddress("els_simpleEleId90cIso", &els_simpleEleId90cIso, &b_els_simpleEleId90cIso);
   fChain->SetBranchAddress("els_simpleEleId85cIso", &els_simpleEleId85cIso, &b_els_simpleEleId85cIso);
   fChain->SetBranchAddress("els_simpleEleId80cIso", &els_simpleEleId80cIso, &b_els_simpleEleId80cIso);
   fChain->SetBranchAddress("els_simpleEleId70cIso", &els_simpleEleId70cIso, &b_els_simpleEleId70cIso);
   fChain->SetBranchAddress("els_simpleEleId60cIso", &els_simpleEleId60cIso, &b_els_simpleEleId60cIso);
   fChain->SetBranchAddress("els_cIso", &els_cIso, &b_els_cIso);
   fChain->SetBranchAddress("els_tIso", &els_tIso, &b_els_tIso);
   fChain->SetBranchAddress("els_ecalIso", &els_ecalIso, &b_els_ecalIso);
   fChain->SetBranchAddress("els_hcalIso", &els_hcalIso, &b_els_hcalIso);
   fChain->SetBranchAddress("els_chi2", &els_chi2, &b_els_chi2);
   fChain->SetBranchAddress("els_charge", &els_charge, &b_els_charge);
   fChain->SetBranchAddress("els_caloEnergy", &els_caloEnergy, &b_els_caloEnergy);
   fChain->SetBranchAddress("els_hadOverEm", &els_hadOverEm, &b_els_hadOverEm);
   fChain->SetBranchAddress("els_hcalOverEcalBc", &els_hcalOverEcalBc, &b_els_hcalOverEcalBc);
   fChain->SetBranchAddress("els_eOverPIn", &els_eOverPIn, &b_els_eOverPIn);
   fChain->SetBranchAddress("els_eSeedOverPOut", &els_eSeedOverPOut, &b_els_eSeedOverPOut);
   fChain->SetBranchAddress("els_sigmaEtaEta", &els_sigmaEtaEta, &b_els_sigmaEtaEta);
   fChain->SetBranchAddress("els_sigmaIEtaIEta", &els_sigmaIEtaIEta, &b_els_sigmaIEtaIEta);
   fChain->SetBranchAddress("els_scEnergy", &els_scEnergy, &b_els_scEnergy);
   fChain->SetBranchAddress("els_scRawEnergy", &els_scRawEnergy, &b_els_scRawEnergy);
   fChain->SetBranchAddress("els_scSeedEnergy", &els_scSeedEnergy, &b_els_scSeedEnergy);
   fChain->SetBranchAddress("els_scEta", &els_scEta, &b_els_scEta);
   fChain->SetBranchAddress("els_scPhi", &els_scPhi, &b_els_scPhi);
   fChain->SetBranchAddress("els_scEtaWidth", &els_scEtaWidth, &b_els_scEtaWidth);
   fChain->SetBranchAddress("els_scPhiWidth", &els_scPhiWidth, &b_els_scPhiWidth);
   fChain->SetBranchAddress("els_scE1x5", &els_scE1x5, &b_els_scE1x5);
   fChain->SetBranchAddress("els_scE2x5Max", &els_scE2x5Max, &b_els_scE2x5Max);
   fChain->SetBranchAddress("els_scE5x5", &els_scE5x5, &b_els_scE5x5);
   fChain->SetBranchAddress("els_isEB", &els_isEB, &b_els_isEB);
   fChain->SetBranchAddress("els_isEE", &els_isEE, &b_els_isEE);
   fChain->SetBranchAddress("els_dEtaIn", &els_dEtaIn, &b_els_dEtaIn);
   fChain->SetBranchAddress("els_dPhiIn", &els_dPhiIn, &b_els_dPhiIn);
   fChain->SetBranchAddress("els_dEtaOut", &els_dEtaOut, &b_els_dEtaOut);
   fChain->SetBranchAddress("els_dPhiOut", &els_dPhiOut, &b_els_dPhiOut);
   fChain->SetBranchAddress("els_numvalhits", &els_numvalhits, &b_els_numvalhits);
   fChain->SetBranchAddress("els_numlosthits", &els_numlosthits, &b_els_numlosthits);
   fChain->SetBranchAddress("els_basicClustersSize", &els_basicClustersSize, &b_els_basicClustersSize);
   fChain->SetBranchAddress("els_tk_pz", &els_tk_pz, &b_els_tk_pz);
   fChain->SetBranchAddress("els_tk_pt", &els_tk_pt, &b_els_tk_pt);
   fChain->SetBranchAddress("els_tk_phi", &els_tk_phi, &b_els_tk_phi);
   fChain->SetBranchAddress("els_tk_eta", &els_tk_eta, &b_els_tk_eta);
   fChain->SetBranchAddress("els_d0dum", &els_d0dum, &b_els_d0dum);
   fChain->SetBranchAddress("els_dz", &els_dz, &b_els_dz);
   fChain->SetBranchAddress("els_vx", &els_vx, &b_els_vx);
   fChain->SetBranchAddress("els_vy", &els_vy, &b_els_vy);
   fChain->SetBranchAddress("els_vz", &els_vz, &b_els_vz);
   fChain->SetBranchAddress("els_ndof", &els_ndof, &b_els_ndof);
   fChain->SetBranchAddress("els_ptError", &els_ptError, &b_els_ptError);
   fChain->SetBranchAddress("els_d0dumError", &els_d0dumError, &b_els_d0dumError);
   fChain->SetBranchAddress("els_dzError", &els_dzError, &b_els_dzError);
   fChain->SetBranchAddress("els_etaError", &els_etaError, &b_els_etaError);
   fChain->SetBranchAddress("els_phiError", &els_phiError, &b_els_phiError);
   fChain->SetBranchAddress("els_tk_charge", &els_tk_charge, &b_els_tk_charge);
   fChain->SetBranchAddress("els_core_ecalDrivenSeed", &els_core_ecalDrivenSeed, &b_els_core_ecalDrivenSeed);
   fChain->SetBranchAddress("els_n_inner_layer", &els_n_inner_layer, &b_els_n_inner_layer);
   fChain->SetBranchAddress("els_n_outer_layer", &els_n_outer_layer, &b_els_n_outer_layer);
   fChain->SetBranchAddress("els_ctf_tk_id", &els_ctf_tk_id, &b_els_ctf_tk_id);
   fChain->SetBranchAddress("els_ctf_tk_charge", &els_ctf_tk_charge, &b_els_ctf_tk_charge);
   fChain->SetBranchAddress("els_ctf_tk_eta", &els_ctf_tk_eta, &b_els_ctf_tk_eta);
   fChain->SetBranchAddress("els_ctf_tk_phi", &els_ctf_tk_phi, &b_els_ctf_tk_phi);
   fChain->SetBranchAddress("els_fbrem", &els_fbrem, &b_els_fbrem);
   fChain->SetBranchAddress("els_shFracInnerHits", &els_shFracInnerHits, &b_els_shFracInnerHits);
   fChain->SetBranchAddress("els_dr03EcalRecHitSumEt", &els_dr03EcalRecHitSumEt, &b_els_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("els_dr03HcalTowerSumEt", &els_dr03HcalTowerSumEt, &b_els_dr03HcalTowerSumEt);
   fChain->SetBranchAddress("els_dr03HcalDepth1TowerSumEt", &els_dr03HcalDepth1TowerSumEt, &b_els_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("els_dr03HcalDepth2TowerSumEt", &els_dr03HcalDepth2TowerSumEt, &b_els_dr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("els_dr03TkSumPt", &els_dr03TkSumPt, &b_els_dr03TkSumPt);
   fChain->SetBranchAddress("els_dr04EcalRecHitSumEt", &els_dr04EcalRecHitSumEt, &b_els_dr04EcalRecHitSumEt);
   fChain->SetBranchAddress("els_dr04HcalTowerSumEt", &els_dr04HcalTowerSumEt, &b_els_dr04HcalTowerSumEt);
   fChain->SetBranchAddress("els_dr04HcalDepth1TowerSumEt", &els_dr04HcalDepth1TowerSumEt, &b_els_dr04HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("els_dr04HcalDepth2TowerSumEt", &els_dr04HcalDepth2TowerSumEt, &b_els_dr04HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("els_dr04TkSumPt", &els_dr04TkSumPt, &b_els_dr04TkSumPt);
   fChain->SetBranchAddress("els_cpx", &els_cpx, &b_els_cpx);
   fChain->SetBranchAddress("els_cpy", &els_cpy, &b_els_cpy);
   fChain->SetBranchAddress("els_cpz", &els_cpz, &b_els_cpz);
   fChain->SetBranchAddress("els_vpx", &els_vpx, &b_els_vpx);
   fChain->SetBranchAddress("els_vpy", &els_vpy, &b_els_vpy);
   fChain->SetBranchAddress("els_vpz", &els_vpz, &b_els_vpz);
   fChain->SetBranchAddress("els_cx", &els_cx, &b_els_cx);
   fChain->SetBranchAddress("els_cy", &els_cy, &b_els_cy);
   fChain->SetBranchAddress("els_cz", &els_cz, &b_els_cz);
   fChain->SetBranchAddress("els_PATpassConversionVeto", &els_PATpassConversionVeto, &b_els_PATpassConversionVeto);
   fChain->SetBranchAddress("Njets_AK5PF", &Njets_AK5PF, &b_Njets_AK5PF);
   fChain->SetBranchAddress("jets_AK5PF_status", &jets_AK5PF_status, &b_jets_AK5PF_status);
   fChain->SetBranchAddress("jets_AK5PF_phi", &jets_AK5PF_phi, &b_jets_AK5PF_phi);
   fChain->SetBranchAddress("jets_AK5PF_pt", &jets_AK5PF_pt, &b_jets_AK5PF_pt);
   fChain->SetBranchAddress("jets_AK5PF_pz", &jets_AK5PF_pz, &b_jets_AK5PF_pz);
   fChain->SetBranchAddress("jets_AK5PF_px", &jets_AK5PF_px, &b_jets_AK5PF_px);
   fChain->SetBranchAddress("jets_AK5PF_py", &jets_AK5PF_py, &b_jets_AK5PF_py);
   fChain->SetBranchAddress("jets_AK5PF_eta", &jets_AK5PF_eta, &b_jets_AK5PF_eta);
   fChain->SetBranchAddress("jets_AK5PF_theta", &jets_AK5PF_theta, &b_jets_AK5PF_theta);
   fChain->SetBranchAddress("jets_AK5PF_et", &jets_AK5PF_et, &b_jets_AK5PF_et);
   fChain->SetBranchAddress("jets_AK5PF_energy", &jets_AK5PF_energy, &b_jets_AK5PF_energy);
   fChain->SetBranchAddress("jets_AK5PF_parton_Id", &jets_AK5PF_parton_Id, &b_jets_AK5PF_parton_Id);
   fChain->SetBranchAddress("jets_AK5PF_parton_motherId", &jets_AK5PF_parton_motherId, &b_jets_AK5PF_parton_motherId);
   fChain->SetBranchAddress("jets_AK5PF_parton_pt", &jets_AK5PF_parton_pt, &b_jets_AK5PF_parton_pt);
   fChain->SetBranchAddress("jets_AK5PF_parton_phi", &jets_AK5PF_parton_phi, &b_jets_AK5PF_parton_phi);
   fChain->SetBranchAddress("jets_AK5PF_parton_eta", &jets_AK5PF_parton_eta, &b_jets_AK5PF_parton_eta);
   fChain->SetBranchAddress("jets_AK5PF_parton_Energy", &jets_AK5PF_parton_Energy, &b_jets_AK5PF_parton_Energy);
   fChain->SetBranchAddress("jets_AK5PF_parton_mass", &jets_AK5PF_parton_mass, &b_jets_AK5PF_parton_mass);
   fChain->SetBranchAddress("jets_AK5PF_gen_et", &jets_AK5PF_gen_et, &b_jets_AK5PF_gen_et);
   fChain->SetBranchAddress("jets_AK5PF_gen_pt", &jets_AK5PF_gen_pt, &b_jets_AK5PF_gen_pt);
   fChain->SetBranchAddress("jets_AK5PF_gen_eta", &jets_AK5PF_gen_eta, &b_jets_AK5PF_gen_eta);
   fChain->SetBranchAddress("jets_AK5PF_gen_phi", &jets_AK5PF_gen_phi, &b_jets_AK5PF_gen_phi);
   fChain->SetBranchAddress("jets_AK5PF_gen_mass", &jets_AK5PF_gen_mass, &b_jets_AK5PF_gen_mass);
   fChain->SetBranchAddress("jets_AK5PF_gen_Energy", &jets_AK5PF_gen_Energy, &b_jets_AK5PF_gen_Energy);
   fChain->SetBranchAddress("jets_AK5PF_gen_Id", &jets_AK5PF_gen_Id, &b_jets_AK5PF_gen_Id);
   fChain->SetBranchAddress("jets_AK5PF_gen_motherID", &jets_AK5PF_gen_motherID, &b_jets_AK5PF_gen_motherID);
   fChain->SetBranchAddress("jets_AK5PF_gen_threeCharge", &jets_AK5PF_gen_threeCharge, &b_jets_AK5PF_gen_threeCharge);
   fChain->SetBranchAddress("jets_AK5PF_partonFlavour", &jets_AK5PF_partonFlavour, &b_jets_AK5PF_partonFlavour);
   fChain->SetBranchAddress("jets_AK5PF_btag_TC_highPur", &jets_AK5PF_btag_TC_highPur, &b_jets_AK5PF_btag_TC_highPur);
   fChain->SetBranchAddress("jets_AK5PF_btag_TC_highEff", &jets_AK5PF_btag_TC_highEff, &b_jets_AK5PF_btag_TC_highEff);
   fChain->SetBranchAddress("jets_AK5PF_btag_jetProb", &jets_AK5PF_btag_jetProb, &b_jets_AK5PF_btag_jetProb);
   fChain->SetBranchAddress("jets_AK5PF_btag_jetBProb", &jets_AK5PF_btag_jetBProb, &b_jets_AK5PF_btag_jetBProb);
   fChain->SetBranchAddress("jets_AK5PF_btag_softEle", &jets_AK5PF_btag_softEle, &b_jets_AK5PF_btag_softEle);
   fChain->SetBranchAddress("jets_AK5PF_btag_softMuon", &jets_AK5PF_btag_softMuon, &b_jets_AK5PF_btag_softMuon);
   fChain->SetBranchAddress("jets_AK5PF_btag_secVertexHighPur", &jets_AK5PF_btag_secVertexHighPur, &b_jets_AK5PF_btag_secVertexHighPur);
   fChain->SetBranchAddress("jets_AK5PF_btag_secVertexHighEff", &jets_AK5PF_btag_secVertexHighEff, &b_jets_AK5PF_btag_secVertexHighEff);
   fChain->SetBranchAddress("jets_AK5PF_btag_secVertexCombined", &jets_AK5PF_btag_secVertexCombined, &b_jets_AK5PF_btag_secVertexCombined);
   fChain->SetBranchAddress("jets_AK5PF_jetCharge", &jets_AK5PF_jetCharge, &b_jets_AK5PF_jetCharge);
   fChain->SetBranchAddress("jets_AK5PF_chgEmE", &jets_AK5PF_chgEmE, &b_jets_AK5PF_chgEmE);
   fChain->SetBranchAddress("jets_AK5PF_chgHadE", &jets_AK5PF_chgHadE, &b_jets_AK5PF_chgHadE);
   fChain->SetBranchAddress("jets_AK5PF_photonEnergy", &jets_AK5PF_photonEnergy, &b_jets_AK5PF_photonEnergy);
   fChain->SetBranchAddress("jets_AK5PF_chgMuE", &jets_AK5PF_chgMuE, &b_jets_AK5PF_chgMuE);
   fChain->SetBranchAddress("jets_AK5PF_chg_Mult", &jets_AK5PF_chg_Mult, &b_jets_AK5PF_chg_Mult);
   fChain->SetBranchAddress("jets_AK5PF_neutralEmE", &jets_AK5PF_neutralEmE, &b_jets_AK5PF_neutralEmE);
   fChain->SetBranchAddress("jets_AK5PF_neutralHadE", &jets_AK5PF_neutralHadE, &b_jets_AK5PF_neutralHadE);
   fChain->SetBranchAddress("jets_AK5PF_neutral_Mult", &jets_AK5PF_neutral_Mult, &b_jets_AK5PF_neutral_Mult);
   fChain->SetBranchAddress("jets_AK5PF_mu_Mult", &jets_AK5PF_mu_Mult, &b_jets_AK5PF_mu_Mult);
   fChain->SetBranchAddress("jets_AK5PF_emf", &jets_AK5PF_emf, &b_jets_AK5PF_emf);
   fChain->SetBranchAddress("jets_AK5PF_ehf", &jets_AK5PF_ehf, &b_jets_AK5PF_ehf);
   fChain->SetBranchAddress("jets_AK5PF_n60", &jets_AK5PF_n60, &b_jets_AK5PF_n60);
   fChain->SetBranchAddress("jets_AK5PF_n90", &jets_AK5PF_n90, &b_jets_AK5PF_n90);
   fChain->SetBranchAddress("jets_AK5PF_etaetaMoment", &jets_AK5PF_etaetaMoment, &b_jets_AK5PF_etaetaMoment);
   fChain->SetBranchAddress("jets_AK5PF_etaphiMoment", &jets_AK5PF_etaphiMoment, &b_jets_AK5PF_etaphiMoment);
   fChain->SetBranchAddress("jets_AK5PF_phiphiMoment", &jets_AK5PF_phiphiMoment, &b_jets_AK5PF_phiphiMoment);
   fChain->SetBranchAddress("jets_AK5PF_n90Hits", &jets_AK5PF_n90Hits, &b_jets_AK5PF_n90Hits);
   fChain->SetBranchAddress("jets_AK5PF_fHPD", &jets_AK5PF_fHPD, &b_jets_AK5PF_fHPD);
   fChain->SetBranchAddress("jets_AK5PF_fRBX", &jets_AK5PF_fRBX, &b_jets_AK5PF_fRBX);
   fChain->SetBranchAddress("jets_AK5PF_hitsInN90", &jets_AK5PF_hitsInN90, &b_jets_AK5PF_hitsInN90);
   fChain->SetBranchAddress("jets_AK5PF_nECALTowers", &jets_AK5PF_nECALTowers, &b_jets_AK5PF_nECALTowers);
   fChain->SetBranchAddress("jets_AK5PF_nHCALTowers", &jets_AK5PF_nHCALTowers, &b_jets_AK5PF_nHCALTowers);
   fChain->SetBranchAddress("jets_AK5PF_fSubDetector1", &jets_AK5PF_fSubDetector1, &b_jets_AK5PF_fSubDetector1);
   fChain->SetBranchAddress("jets_AK5PF_fSubDetector2", &jets_AK5PF_fSubDetector2, &b_jets_AK5PF_fSubDetector2);
   fChain->SetBranchAddress("jets_AK5PF_fSubDetector3", &jets_AK5PF_fSubDetector3, &b_jets_AK5PF_fSubDetector3);
   fChain->SetBranchAddress("jets_AK5PF_fSubDetector4", &jets_AK5PF_fSubDetector4, &b_jets_AK5PF_fSubDetector4);
   fChain->SetBranchAddress("jets_AK5PF_area", &jets_AK5PF_area, &b_jets_AK5PF_area);
   fChain->SetBranchAddress("jets_AK5PF_corrFactorRaw", &jets_AK5PF_corrFactorRaw, &b_jets_AK5PF_corrFactorRaw);
   fChain->SetBranchAddress("jets_AK5PF_rawPt", &jets_AK5PF_rawPt, &b_jets_AK5PF_rawPt);
   fChain->SetBranchAddress("jets_AK5PF_mass", &jets_AK5PF_mass, &b_jets_AK5PF_mass);
   fChain->SetBranchAddress("Njets_AK5PFclean", &Njets_AK5PFclean, &b_Njets_AK5PFclean);
   fChain->SetBranchAddress("jets_AK5PFclean_status", &jets_AK5PFclean_status, &b_jets_AK5PFclean_status);
   fChain->SetBranchAddress("jets_AK5PFclean_phi", &jets_AK5PFclean_phi, &b_jets_AK5PFclean_phi);
   fChain->SetBranchAddress("jets_AK5PFclean_pt", &jets_AK5PFclean_pt, &b_jets_AK5PFclean_pt);
   fChain->SetBranchAddress("jets_AK5PFclean_pz", &jets_AK5PFclean_pz, &b_jets_AK5PFclean_pz);
   fChain->SetBranchAddress("jets_AK5PFclean_px", &jets_AK5PFclean_px, &b_jets_AK5PFclean_px);
   fChain->SetBranchAddress("jets_AK5PFclean_py", &jets_AK5PFclean_py, &b_jets_AK5PFclean_py);
   fChain->SetBranchAddress("jets_AK5PFclean_eta", &jets_AK5PFclean_eta, &b_jets_AK5PFclean_eta);
   fChain->SetBranchAddress("jets_AK5PFclean_theta", &jets_AK5PFclean_theta, &b_jets_AK5PFclean_theta);
   fChain->SetBranchAddress("jets_AK5PFclean_et", &jets_AK5PFclean_et, &b_jets_AK5PFclean_et);
   fChain->SetBranchAddress("jets_AK5PFclean_energy", &jets_AK5PFclean_energy, &b_jets_AK5PFclean_energy);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_Id", &jets_AK5PFclean_parton_Id, &b_jets_AK5PFclean_parton_Id);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_motherId", &jets_AK5PFclean_parton_motherId, &b_jets_AK5PFclean_parton_motherId);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_pt", &jets_AK5PFclean_parton_pt, &b_jets_AK5PFclean_parton_pt);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_phi", &jets_AK5PFclean_parton_phi, &b_jets_AK5PFclean_parton_phi);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_eta", &jets_AK5PFclean_parton_eta, &b_jets_AK5PFclean_parton_eta);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_Energy", &jets_AK5PFclean_parton_Energy, &b_jets_AK5PFclean_parton_Energy);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_mass", &jets_AK5PFclean_parton_mass, &b_jets_AK5PFclean_parton_mass);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_et", &jets_AK5PFclean_gen_et, &b_jets_AK5PFclean_gen_et);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_pt", &jets_AK5PFclean_gen_pt, &b_jets_AK5PFclean_gen_pt);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_eta", &jets_AK5PFclean_gen_eta, &b_jets_AK5PFclean_gen_eta);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_phi", &jets_AK5PFclean_gen_phi, &b_jets_AK5PFclean_gen_phi);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_mass", &jets_AK5PFclean_gen_mass, &b_jets_AK5PFclean_gen_mass);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_Energy", &jets_AK5PFclean_gen_Energy, &b_jets_AK5PFclean_gen_Energy);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_Id", &jets_AK5PFclean_gen_Id, &b_jets_AK5PFclean_gen_Id);
   fChain->SetBranchAddress("jets_AK5PFclean_partonFlavour", &jets_AK5PFclean_partonFlavour, &b_jets_AK5PFclean_partonFlavour);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_TC_highPur", &jets_AK5PFclean_btag_TC_highPur, &b_jets_AK5PFclean_btag_TC_highPur);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_TC_highEff", &jets_AK5PFclean_btag_TC_highEff, &b_jets_AK5PFclean_btag_TC_highEff);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_jetProb", &jets_AK5PFclean_btag_jetProb, &b_jets_AK5PFclean_btag_jetProb);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_jetBProb", &jets_AK5PFclean_btag_jetBProb, &b_jets_AK5PFclean_btag_jetBProb);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_softEle", &jets_AK5PFclean_btag_softEle, &b_jets_AK5PFclean_btag_softEle);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_softMuon", &jets_AK5PFclean_btag_softMuon, &b_jets_AK5PFclean_btag_softMuon);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_secVertexHighPur", &jets_AK5PFclean_btag_secVertexHighPur, &b_jets_AK5PFclean_btag_secVertexHighPur);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_secVertexHighEff", &jets_AK5PFclean_btag_secVertexHighEff, &b_jets_AK5PFclean_btag_secVertexHighEff);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_secVertexCombined", &jets_AK5PFclean_btag_secVertexCombined, &b_jets_AK5PFclean_btag_secVertexCombined);
   fChain->SetBranchAddress("jets_AK5PFclean_jetCharge", &jets_AK5PFclean_jetCharge, &b_jets_AK5PFclean_jetCharge);
   fChain->SetBranchAddress("jets_AK5PFclean_chgEmE", &jets_AK5PFclean_chgEmE, &b_jets_AK5PFclean_chgEmE);
   fChain->SetBranchAddress("jets_AK5PFclean_chgHadE", &jets_AK5PFclean_chgHadE, &b_jets_AK5PFclean_chgHadE);
   fChain->SetBranchAddress("jets_AK5PFclean_photonEnergy", &jets_AK5PFclean_photonEnergy, &b_jets_AK5PFclean_photonEnergy);
   fChain->SetBranchAddress("jets_AK5PFclean_chgMuE", &jets_AK5PFclean_chgMuE, &b_jets_AK5PFclean_chgMuE);
   fChain->SetBranchAddress("jets_AK5PFclean_chg_Mult", &jets_AK5PFclean_chg_Mult, &b_jets_AK5PFclean_chg_Mult);
   fChain->SetBranchAddress("jets_AK5PFclean_neutralEmE", &jets_AK5PFclean_neutralEmE, &b_jets_AK5PFclean_neutralEmE);
   fChain->SetBranchAddress("jets_AK5PFclean_neutralHadE", &jets_AK5PFclean_neutralHadE, &b_jets_AK5PFclean_neutralHadE);
   fChain->SetBranchAddress("jets_AK5PFclean_neutral_Mult", &jets_AK5PFclean_neutral_Mult, &b_jets_AK5PFclean_neutral_Mult);
   fChain->SetBranchAddress("jets_AK5PFclean_mu_Mult", &jets_AK5PFclean_mu_Mult, &b_jets_AK5PFclean_mu_Mult);
   fChain->SetBranchAddress("jets_AK5PFclean_emf", &jets_AK5PFclean_emf, &b_jets_AK5PFclean_emf);
   fChain->SetBranchAddress("jets_AK5PFclean_ehf", &jets_AK5PFclean_ehf, &b_jets_AK5PFclean_ehf);
   fChain->SetBranchAddress("jets_AK5PFclean_n60", &jets_AK5PFclean_n60, &b_jets_AK5PFclean_n60);
   fChain->SetBranchAddress("jets_AK5PFclean_n90", &jets_AK5PFclean_n90, &b_jets_AK5PFclean_n90);
   fChain->SetBranchAddress("jets_AK5PFclean_etaetaMoment", &jets_AK5PFclean_etaetaMoment, &b_jets_AK5PFclean_etaetaMoment);
   fChain->SetBranchAddress("jets_AK5PFclean_etaphiMoment", &jets_AK5PFclean_etaphiMoment, &b_jets_AK5PFclean_etaphiMoment);
   fChain->SetBranchAddress("jets_AK5PFclean_phiphiMoment", &jets_AK5PFclean_phiphiMoment, &b_jets_AK5PFclean_phiphiMoment);
   fChain->SetBranchAddress("jets_AK5PFclean_n90Hits", &jets_AK5PFclean_n90Hits, &b_jets_AK5PFclean_n90Hits);
   fChain->SetBranchAddress("jets_AK5PFclean_fHPD", &jets_AK5PFclean_fHPD, &b_jets_AK5PFclean_fHPD);
   fChain->SetBranchAddress("jets_AK5PFclean_fRBX", &jets_AK5PFclean_fRBX, &b_jets_AK5PFclean_fRBX);
   fChain->SetBranchAddress("jets_AK5PFclean_hitsInN90", &jets_AK5PFclean_hitsInN90, &b_jets_AK5PFclean_hitsInN90);
   fChain->SetBranchAddress("jets_AK5PFclean_nECALTowers", &jets_AK5PFclean_nECALTowers, &b_jets_AK5PFclean_nECALTowers);
   fChain->SetBranchAddress("jets_AK5PFclean_nHCALTowers", &jets_AK5PFclean_nHCALTowers, &b_jets_AK5PFclean_nHCALTowers);
   fChain->SetBranchAddress("jets_AK5PFclean_fSubDetector1", &jets_AK5PFclean_fSubDetector1, &b_jets_AK5PFclean_fSubDetector1);
   fChain->SetBranchAddress("jets_AK5PFclean_fSubDetector2", &jets_AK5PFclean_fSubDetector2, &b_jets_AK5PFclean_fSubDetector2);
   fChain->SetBranchAddress("jets_AK5PFclean_fSubDetector3", &jets_AK5PFclean_fSubDetector3, &b_jets_AK5PFclean_fSubDetector3);
   fChain->SetBranchAddress("jets_AK5PFclean_fSubDetector4", &jets_AK5PFclean_fSubDetector4, &b_jets_AK5PFclean_fSubDetector4);
   fChain->SetBranchAddress("jets_AK5PFclean_area", &jets_AK5PFclean_area, &b_jets_AK5PFclean_area);
   fChain->SetBranchAddress("jets_AK5PFclean_corrFactorRaw", &jets_AK5PFclean_corrFactorRaw, &b_jets_AK5PFclean_corrFactorRaw);
   fChain->SetBranchAddress("jets_AK5PFclean_rawPt", &jets_AK5PFclean_rawPt, &b_jets_AK5PFclean_rawPt);
   fChain->SetBranchAddress("jets_AK5PFclean_mass", &jets_AK5PFclean_mass, &b_jets_AK5PFclean_mass);
   fChain->SetBranchAddress("Nmc_doc", &Nmc_doc, &b_Nmc_doc);
   fChain->SetBranchAddress("mc_doc_id", &mc_doc_id, &b_mc_doc_id);
   fChain->SetBranchAddress("mc_doc_pt", &mc_doc_pt, &b_mc_doc_pt);
   fChain->SetBranchAddress("mc_doc_px", &mc_doc_px, &b_mc_doc_px);
   fChain->SetBranchAddress("mc_doc_py", &mc_doc_py, &b_mc_doc_py);
   fChain->SetBranchAddress("mc_doc_pz", &mc_doc_pz, &b_mc_doc_pz);
   fChain->SetBranchAddress("mc_doc_eta", &mc_doc_eta, &b_mc_doc_eta);
   fChain->SetBranchAddress("mc_doc_phi", &mc_doc_phi, &b_mc_doc_phi);
   fChain->SetBranchAddress("mc_doc_theta", &mc_doc_theta, &b_mc_doc_theta);
   fChain->SetBranchAddress("mc_doc_energy", &mc_doc_energy, &b_mc_doc_energy);
   fChain->SetBranchAddress("mc_doc_status", &mc_doc_status, &b_mc_doc_status);
   fChain->SetBranchAddress("mc_doc_charge", &mc_doc_charge, &b_mc_doc_charge);
   fChain->SetBranchAddress("mc_doc_mother_id", &mc_doc_mother_id, &b_mc_doc_mother_id);
   fChain->SetBranchAddress("mc_doc_grandmother_id", &mc_doc_grandmother_id, &b_mc_doc_grandmother_id);
   fChain->SetBranchAddress("mc_doc_ggrandmother_id", &mc_doc_ggrandmother_id, &b_mc_doc_ggrandmother_id);
   fChain->SetBranchAddress("mc_doc_mother_pt", &mc_doc_mother_pt, &b_mc_doc_mother_pt);
   fChain->SetBranchAddress("mc_doc_vertex_x", &mc_doc_vertex_x, &b_mc_doc_vertex_x);
   fChain->SetBranchAddress("mc_doc_vertex_y", &mc_doc_vertex_y, &b_mc_doc_vertex_y);
   fChain->SetBranchAddress("mc_doc_vertex_z", &mc_doc_vertex_z, &b_mc_doc_vertex_z);
   fChain->SetBranchAddress("mc_doc_mass", &mc_doc_mass, &b_mc_doc_mass);
   fChain->SetBranchAddress("mc_doc_numOfDaughters", &mc_doc_numOfDaughters, &b_mc_doc_numOfDaughters);
   fChain->SetBranchAddress("mc_doc_numOfMothers", &mc_doc_numOfMothers, &b_mc_doc_numOfMothers);
   fChain->SetBranchAddress("Nmc_electrons", &Nmc_electrons, &b_Nmc_electrons);
   fChain->SetBranchAddress("mc_electrons_id", &mc_electrons_id, &b_mc_electrons_id);
   fChain->SetBranchAddress("mc_electrons_pt", &mc_electrons_pt, &b_mc_electrons_pt);
   fChain->SetBranchAddress("mc_electrons_px", &mc_electrons_px, &b_mc_electrons_px);
   fChain->SetBranchAddress("mc_electrons_py", &mc_electrons_py, &b_mc_electrons_py);
   fChain->SetBranchAddress("mc_electrons_pz", &mc_electrons_pz, &b_mc_electrons_pz);
   fChain->SetBranchAddress("mc_electrons_eta", &mc_electrons_eta, &b_mc_electrons_eta);
   fChain->SetBranchAddress("mc_electrons_phi", &mc_electrons_phi, &b_mc_electrons_phi);
   fChain->SetBranchAddress("mc_electrons_theta", &mc_electrons_theta, &b_mc_electrons_theta);
   fChain->SetBranchAddress("mc_electrons_status", &mc_electrons_status, &b_mc_electrons_status);
   fChain->SetBranchAddress("mc_electrons_energy", &mc_electrons_energy, &b_mc_electrons_energy);
   fChain->SetBranchAddress("mc_electrons_charge", &mc_electrons_charge, &b_mc_electrons_charge);
   fChain->SetBranchAddress("mc_electrons_mother_id", &mc_electrons_mother_id, &b_mc_electrons_mother_id);
   fChain->SetBranchAddress("mc_electrons_mother_pt", &mc_electrons_mother_pt, &b_mc_electrons_mother_pt);
   fChain->SetBranchAddress("mc_electrons_grandmother_id", &mc_electrons_grandmother_id, &b_mc_electrons_grandmother_id);
   fChain->SetBranchAddress("mc_electrons_ggrandmother_id", &mc_electrons_ggrandmother_id, &b_mc_electrons_ggrandmother_id);
   fChain->SetBranchAddress("mc_electrons_vertex_x", &mc_electrons_vertex_x, &b_mc_electrons_vertex_x);
   fChain->SetBranchAddress("mc_electrons_vertex_y", &mc_electrons_vertex_y, &b_mc_electrons_vertex_y);
   fChain->SetBranchAddress("mc_electrons_vertex_z", &mc_electrons_vertex_z, &b_mc_electrons_vertex_z);
   fChain->SetBranchAddress("mc_electrons_mass", &mc_electrons_mass, &b_mc_electrons_mass);
   fChain->SetBranchAddress("mc_electrons_numOfDaughters", &mc_electrons_numOfDaughters, &b_mc_electrons_numOfDaughters);
   fChain->SetBranchAddress("Nmc_mus", &Nmc_mus, &b_Nmc_mus);
   fChain->SetBranchAddress("mc_mus_id", &mc_mus_id, &b_mc_mus_id);
   fChain->SetBranchAddress("mc_mus_pt", &mc_mus_pt, &b_mc_mus_pt);
   fChain->SetBranchAddress("mc_mus_px", &mc_mus_px, &b_mc_mus_px);
   fChain->SetBranchAddress("mc_mus_py", &mc_mus_py, &b_mc_mus_py);
   fChain->SetBranchAddress("mc_mus_pz", &mc_mus_pz, &b_mc_mus_pz);
   fChain->SetBranchAddress("mc_mus_eta", &mc_mus_eta, &b_mc_mus_eta);
   fChain->SetBranchAddress("mc_mus_phi", &mc_mus_phi, &b_mc_mus_phi);
   fChain->SetBranchAddress("mc_mus_theta", &mc_mus_theta, &b_mc_mus_theta);
   fChain->SetBranchAddress("mc_mus_status", &mc_mus_status, &b_mc_mus_status);
   fChain->SetBranchAddress("mc_mus_energy", &mc_mus_energy, &b_mc_mus_energy);
   fChain->SetBranchAddress("mc_mus_charge", &mc_mus_charge, &b_mc_mus_charge);
   fChain->SetBranchAddress("mc_mus_mother_id", &mc_mus_mother_id, &b_mc_mus_mother_id);
   fChain->SetBranchAddress("mc_mus_mother_pt", &mc_mus_mother_pt, &b_mc_mus_mother_pt);
   fChain->SetBranchAddress("mc_mus_grandmother_id", &mc_mus_grandmother_id, &b_mc_mus_grandmother_id);
   fChain->SetBranchAddress("mc_mus_ggrandmother_id", &mc_mus_ggrandmother_id, &b_mc_mus_ggrandmother_id);
   fChain->SetBranchAddress("mc_mus_vertex_x", &mc_mus_vertex_x, &b_mc_mus_vertex_x);
   fChain->SetBranchAddress("mc_mus_vertex_y", &mc_mus_vertex_y, &b_mc_mus_vertex_y);
   fChain->SetBranchAddress("mc_mus_vertex_z", &mc_mus_vertex_z, &b_mc_mus_vertex_z);
   fChain->SetBranchAddress("mc_mus_mass", &mc_mus_mass, &b_mc_mus_mass);
   fChain->SetBranchAddress("mc_mus_numOfDaughters", &mc_mus_numOfDaughters, &b_mc_mus_numOfDaughters);
   fChain->SetBranchAddress("Nmc_nues", &Nmc_nues, &b_Nmc_nues);
   fChain->SetBranchAddress("mc_nues_id", &mc_nues_id, &b_mc_nues_id);
   fChain->SetBranchAddress("mc_nues_pt", &mc_nues_pt, &b_mc_nues_pt);
   fChain->SetBranchAddress("mc_nues_px", &mc_nues_px, &b_mc_nues_px);
   fChain->SetBranchAddress("mc_nues_py", &mc_nues_py, &b_mc_nues_py);
   fChain->SetBranchAddress("mc_nues_pz", &mc_nues_pz, &b_mc_nues_pz);
   fChain->SetBranchAddress("mc_nues_eta", &mc_nues_eta, &b_mc_nues_eta);
   fChain->SetBranchAddress("mc_nues_phi", &mc_nues_phi, &b_mc_nues_phi);
   fChain->SetBranchAddress("mc_nues_theta", &mc_nues_theta, &b_mc_nues_theta);
   fChain->SetBranchAddress("mc_nues_status", &mc_nues_status, &b_mc_nues_status);
   fChain->SetBranchAddress("mc_nues_energy", &mc_nues_energy, &b_mc_nues_energy);
   fChain->SetBranchAddress("mc_nues_charge", &mc_nues_charge, &b_mc_nues_charge);
   fChain->SetBranchAddress("mc_nues_mother_id", &mc_nues_mother_id, &b_mc_nues_mother_id);
   fChain->SetBranchAddress("mc_nues_mother_pt", &mc_nues_mother_pt, &b_mc_nues_mother_pt);
   fChain->SetBranchAddress("mc_nues_grandmother_id", &mc_nues_grandmother_id, &b_mc_nues_grandmother_id);
   fChain->SetBranchAddress("mc_nues_ggrandmother_id", &mc_nues_ggrandmother_id, &b_mc_nues_ggrandmother_id);
   fChain->SetBranchAddress("mc_nues_vertex_x", &mc_nues_vertex_x, &b_mc_nues_vertex_x);
   fChain->SetBranchAddress("mc_nues_vertex_y", &mc_nues_vertex_y, &b_mc_nues_vertex_y);
   fChain->SetBranchAddress("mc_nues_vertex_z", &mc_nues_vertex_z, &b_mc_nues_vertex_z);
   fChain->SetBranchAddress("mc_nues_mass", &mc_nues_mass, &b_mc_nues_mass);
   fChain->SetBranchAddress("mc_nues_numOfDaughters", &mc_nues_numOfDaughters, &b_mc_nues_numOfDaughters);
   fChain->SetBranchAddress("Nmc_numus", &Nmc_numus, &b_Nmc_numus);
   fChain->SetBranchAddress("mc_numus_id", &mc_numus_id, &b_mc_numus_id);
   fChain->SetBranchAddress("mc_numus_pt", &mc_numus_pt, &b_mc_numus_pt);
   fChain->SetBranchAddress("mc_numus_px", &mc_numus_px, &b_mc_numus_px);
   fChain->SetBranchAddress("mc_numus_py", &mc_numus_py, &b_mc_numus_py);
   fChain->SetBranchAddress("mc_numus_pz", &mc_numus_pz, &b_mc_numus_pz);
   fChain->SetBranchAddress("mc_numus_eta", &mc_numus_eta, &b_mc_numus_eta);
   fChain->SetBranchAddress("mc_numus_phi", &mc_numus_phi, &b_mc_numus_phi);
   fChain->SetBranchAddress("mc_numus_theta", &mc_numus_theta, &b_mc_numus_theta);
   fChain->SetBranchAddress("mc_numus_status", &mc_numus_status, &b_mc_numus_status);
   fChain->SetBranchAddress("mc_numus_energy", &mc_numus_energy, &b_mc_numus_energy);
   fChain->SetBranchAddress("mc_numus_charge", &mc_numus_charge, &b_mc_numus_charge);
   fChain->SetBranchAddress("mc_numus_mother_id", &mc_numus_mother_id, &b_mc_numus_mother_id);
   fChain->SetBranchAddress("mc_numus_mother_pt", &mc_numus_mother_pt, &b_mc_numus_mother_pt);
   fChain->SetBranchAddress("mc_numus_grandmother_id", &mc_numus_grandmother_id, &b_mc_numus_grandmother_id);
   fChain->SetBranchAddress("mc_numus_ggrandmother_id", &mc_numus_ggrandmother_id, &b_mc_numus_ggrandmother_id);
   fChain->SetBranchAddress("mc_numus_vertex_x", &mc_numus_vertex_x, &b_mc_numus_vertex_x);
   fChain->SetBranchAddress("mc_numus_vertex_y", &mc_numus_vertex_y, &b_mc_numus_vertex_y);
   fChain->SetBranchAddress("mc_numus_vertex_z", &mc_numus_vertex_z, &b_mc_numus_vertex_z);
   fChain->SetBranchAddress("mc_numus_mass", &mc_numus_mass, &b_mc_numus_mass);
   fChain->SetBranchAddress("mc_numus_numOfDaughters", &mc_numus_numOfDaughters, &b_mc_numus_numOfDaughters);
   fChain->SetBranchAddress("Nmc_nutaus", &Nmc_nutaus, &b_Nmc_nutaus);
   fChain->SetBranchAddress("mc_nutaus_id", &mc_nutaus_id, &b_mc_nutaus_id);
   fChain->SetBranchAddress("mc_nutaus_pt", &mc_nutaus_pt, &b_mc_nutaus_pt);
   fChain->SetBranchAddress("mc_nutaus_px", &mc_nutaus_px, &b_mc_nutaus_px);
   fChain->SetBranchAddress("mc_nutaus_py", &mc_nutaus_py, &b_mc_nutaus_py);
   fChain->SetBranchAddress("mc_nutaus_pz", &mc_nutaus_pz, &b_mc_nutaus_pz);
   fChain->SetBranchAddress("mc_nutaus_eta", &mc_nutaus_eta, &b_mc_nutaus_eta);
   fChain->SetBranchAddress("mc_nutaus_phi", &mc_nutaus_phi, &b_mc_nutaus_phi);
   fChain->SetBranchAddress("mc_nutaus_theta", &mc_nutaus_theta, &b_mc_nutaus_theta);
   fChain->SetBranchAddress("mc_nutaus_status", &mc_nutaus_status, &b_mc_nutaus_status);
   fChain->SetBranchAddress("mc_nutaus_energy", &mc_nutaus_energy, &b_mc_nutaus_energy);
   fChain->SetBranchAddress("mc_nutaus_charge", &mc_nutaus_charge, &b_mc_nutaus_charge);
   fChain->SetBranchAddress("mc_nutaus_mother_id", &mc_nutaus_mother_id, &b_mc_nutaus_mother_id);
   fChain->SetBranchAddress("mc_nutaus_mother_pt", &mc_nutaus_mother_pt, &b_mc_nutaus_mother_pt);
   fChain->SetBranchAddress("mc_nutaus_grandmother_id", &mc_nutaus_grandmother_id, &b_mc_nutaus_grandmother_id);
   fChain->SetBranchAddress("mc_nutaus_ggrandmother_id", &mc_nutaus_ggrandmother_id, &b_mc_nutaus_ggrandmother_id);
   fChain->SetBranchAddress("mc_nutaus_vertex_x", &mc_nutaus_vertex_x, &b_mc_nutaus_vertex_x);
   fChain->SetBranchAddress("mc_nutaus_vertex_y", &mc_nutaus_vertex_y, &b_mc_nutaus_vertex_y);
   fChain->SetBranchAddress("mc_nutaus_vertex_z", &mc_nutaus_vertex_z, &b_mc_nutaus_vertex_z);
   fChain->SetBranchAddress("mc_nutaus_mass", &mc_nutaus_mass, &b_mc_nutaus_mass);
   fChain->SetBranchAddress("mc_nutaus_numOfDaughters", &mc_nutaus_numOfDaughters, &b_mc_nutaus_numOfDaughters);
   fChain->SetBranchAddress("Nmc_photons", &Nmc_photons, &b_Nmc_photons);
   fChain->SetBranchAddress("mc_photons_id", &mc_photons_id, &b_mc_photons_id);
   fChain->SetBranchAddress("mc_photons_pt", &mc_photons_pt, &b_mc_photons_pt);
   fChain->SetBranchAddress("mc_photons_px", &mc_photons_px, &b_mc_photons_px);
   fChain->SetBranchAddress("mc_photons_py", &mc_photons_py, &b_mc_photons_py);
   fChain->SetBranchAddress("mc_photons_pz", &mc_photons_pz, &b_mc_photons_pz);
   fChain->SetBranchAddress("mc_photons_eta", &mc_photons_eta, &b_mc_photons_eta);
   fChain->SetBranchAddress("mc_photons_phi", &mc_photons_phi, &b_mc_photons_phi);
   fChain->SetBranchAddress("mc_photons_theta", &mc_photons_theta, &b_mc_photons_theta);
   fChain->SetBranchAddress("mc_photons_status", &mc_photons_status, &b_mc_photons_status);
   fChain->SetBranchAddress("mc_photons_energy", &mc_photons_energy, &b_mc_photons_energy);
   fChain->SetBranchAddress("mc_photons_charge", &mc_photons_charge, &b_mc_photons_charge);
   fChain->SetBranchAddress("mc_photons_mother_id", &mc_photons_mother_id, &b_mc_photons_mother_id);
   fChain->SetBranchAddress("mc_photons_mother_pt", &mc_photons_mother_pt, &b_mc_photons_mother_pt);
   fChain->SetBranchAddress("mc_photons_grandmother_id", &mc_photons_grandmother_id, &b_mc_photons_grandmother_id);
   fChain->SetBranchAddress("mc_photons_ggrandmother_id", &mc_photons_ggrandmother_id, &b_mc_photons_ggrandmother_id);
   fChain->SetBranchAddress("mc_photons_vertex_x", &mc_photons_vertex_x, &b_mc_photons_vertex_x);
   fChain->SetBranchAddress("mc_photons_vertex_y", &mc_photons_vertex_y, &b_mc_photons_vertex_y);
   fChain->SetBranchAddress("mc_photons_vertex_z", &mc_photons_vertex_z, &b_mc_photons_vertex_z);
   fChain->SetBranchAddress("mc_photons_mass", &mc_photons_mass, &b_mc_photons_mass);
   fChain->SetBranchAddress("mc_photons_numOfDaughters", &mc_photons_numOfDaughters, &b_mc_photons_numOfDaughters);
   fChain->SetBranchAddress("Nmc_taus", &Nmc_taus, &b_Nmc_taus);
   fChain->SetBranchAddress("mc_taus_id", &mc_taus_id, &b_mc_taus_id);
   fChain->SetBranchAddress("mc_taus_pt", &mc_taus_pt, &b_mc_taus_pt);
   fChain->SetBranchAddress("mc_taus_px", &mc_taus_px, &b_mc_taus_px);
   fChain->SetBranchAddress("mc_taus_py", &mc_taus_py, &b_mc_taus_py);
   fChain->SetBranchAddress("mc_taus_pz", &mc_taus_pz, &b_mc_taus_pz);
   fChain->SetBranchAddress("mc_taus_eta", &mc_taus_eta, &b_mc_taus_eta);
   fChain->SetBranchAddress("mc_taus_phi", &mc_taus_phi, &b_mc_taus_phi);
   fChain->SetBranchAddress("mc_taus_theta", &mc_taus_theta, &b_mc_taus_theta);
   fChain->SetBranchAddress("mc_taus_status", &mc_taus_status, &b_mc_taus_status);
   fChain->SetBranchAddress("mc_taus_energy", &mc_taus_energy, &b_mc_taus_energy);
   fChain->SetBranchAddress("mc_taus_charge", &mc_taus_charge, &b_mc_taus_charge);
   fChain->SetBranchAddress("mc_taus_mother_id", &mc_taus_mother_id, &b_mc_taus_mother_id);
   fChain->SetBranchAddress("mc_taus_mother_pt", &mc_taus_mother_pt, &b_mc_taus_mother_pt);
   fChain->SetBranchAddress("mc_taus_grandmother_id", &mc_taus_grandmother_id, &b_mc_taus_grandmother_id);
   fChain->SetBranchAddress("mc_taus_ggrandmother_id", &mc_taus_ggrandmother_id, &b_mc_taus_ggrandmother_id);
   fChain->SetBranchAddress("mc_taus_vertex_x", &mc_taus_vertex_x, &b_mc_taus_vertex_x);
   fChain->SetBranchAddress("mc_taus_vertex_y", &mc_taus_vertex_y, &b_mc_taus_vertex_y);
   fChain->SetBranchAddress("mc_taus_vertex_z", &mc_taus_vertex_z, &b_mc_taus_vertex_z);
   fChain->SetBranchAddress("mc_taus_mass", &mc_taus_mass, &b_mc_taus_mass);
   fChain->SetBranchAddress("mc_taus_numOfDaughters", &mc_taus_numOfDaughters, &b_mc_taus_numOfDaughters);
   fChain->SetBranchAddress("Nmets_AK5", &Nmets_AK5, &b_Nmets_AK5);
   fChain->SetBranchAddress("mets_AK5_et", &mets_AK5_et, &b_mets_AK5_et);
   fChain->SetBranchAddress("mets_AK5_phi", &mets_AK5_phi, &b_mets_AK5_phi);
   fChain->SetBranchAddress("mets_AK5_ex", &mets_AK5_ex, &b_mets_AK5_ex);
   fChain->SetBranchAddress("mets_AK5_ey", &mets_AK5_ey, &b_mets_AK5_ey);
   fChain->SetBranchAddress("mets_AK5_gen_et", &mets_AK5_gen_et, &b_mets_AK5_gen_et);
   fChain->SetBranchAddress("mets_AK5_gen_phi", &mets_AK5_gen_phi, &b_mets_AK5_gen_phi);
   fChain->SetBranchAddress("mets_AK5_sign", &mets_AK5_sign, &b_mets_AK5_sign);
   fChain->SetBranchAddress("mets_AK5_sumEt", &mets_AK5_sumEt, &b_mets_AK5_sumEt);
   fChain->SetBranchAddress("mets_AK5_unCPhi", &mets_AK5_unCPhi, &b_mets_AK5_unCPhi);
   fChain->SetBranchAddress("mets_AK5_unCPt", &mets_AK5_unCPt, &b_mets_AK5_unCPt);
   fChain->SetBranchAddress("Nmus", &Nmus, &b_Nmus);
   fChain->SetBranchAddress("mus_energy", &mus_energy, &b_mus_energy);
   fChain->SetBranchAddress("mus_et", &mus_et, &b_mus_et);
   fChain->SetBranchAddress("mus_eta", &mus_eta, &b_mus_eta);
   fChain->SetBranchAddress("mus_phi", &mus_phi, &b_mus_phi);
   fChain->SetBranchAddress("mus_pt", &mus_pt, &b_mus_pt);
   fChain->SetBranchAddress("mus_px", &mus_px, &b_mus_px);
   fChain->SetBranchAddress("mus_py", &mus_py, &b_mus_py);
   fChain->SetBranchAddress("mus_pz", &mus_pz, &b_mus_pz);
   fChain->SetBranchAddress("mus_status", &mus_status, &b_mus_status);
   fChain->SetBranchAddress("mus_theta", &mus_theta, &b_mus_theta);
   fChain->SetBranchAddress("mus_gen_id", &mus_gen_id, &b_mus_gen_id);
   fChain->SetBranchAddress("mus_gen_phi", &mus_gen_phi, &b_mus_gen_phi);
   fChain->SetBranchAddress("mus_gen_pt", &mus_gen_pt, &b_mus_gen_pt);
   fChain->SetBranchAddress("mus_gen_pz", &mus_gen_pz, &b_mus_gen_pz);
   fChain->SetBranchAddress("mus_gen_px", &mus_gen_px, &b_mus_gen_px);
   fChain->SetBranchAddress("mus_gen_py", &mus_gen_py, &b_mus_gen_py);
   fChain->SetBranchAddress("mus_gen_eta", &mus_gen_eta, &b_mus_gen_eta);
   fChain->SetBranchAddress("mus_gen_theta", &mus_gen_theta, &b_mus_gen_theta);
   fChain->SetBranchAddress("mus_gen_et", &mus_gen_et, &b_mus_gen_et);
   fChain->SetBranchAddress("mus_gen_mother_id", &mus_gen_mother_id, &b_mus_gen_mother_id);
   fChain->SetBranchAddress("mus_gen_mother_phi", &mus_gen_mother_phi, &b_mus_gen_mother_phi);
   fChain->SetBranchAddress("mus_gen_mother_pt", &mus_gen_mother_pt, &b_mus_gen_mother_pt);
   fChain->SetBranchAddress("mus_gen_mother_pz", &mus_gen_mother_pz, &b_mus_gen_mother_pz);
   fChain->SetBranchAddress("mus_gen_mother_px", &mus_gen_mother_px, &b_mus_gen_mother_px);
   fChain->SetBranchAddress("mus_gen_mother_py", &mus_gen_mother_py, &b_mus_gen_mother_py);
   fChain->SetBranchAddress("mus_gen_mother_eta", &mus_gen_mother_eta, &b_mus_gen_mother_eta);
   fChain->SetBranchAddress("mus_gen_mother_theta", &mus_gen_mother_theta, &b_mus_gen_mother_theta);
   fChain->SetBranchAddress("mus_gen_mother_et", &mus_gen_mother_et, &b_mus_gen_mother_et);
   fChain->SetBranchAddress("mus_tkHits", &mus_tkHits, &b_mus_tkHits);
   fChain->SetBranchAddress("mus_cIso", &mus_cIso, &b_mus_cIso);
   fChain->SetBranchAddress("mus_tIso", &mus_tIso, &b_mus_tIso);
   fChain->SetBranchAddress("mus_ecalIso", &mus_ecalIso, &b_mus_ecalIso);
   fChain->SetBranchAddress("mus_hcalIso", &mus_hcalIso, &b_mus_hcalIso);
   fChain->SetBranchAddress("mus_ecalvetoDep", &mus_ecalvetoDep, &b_mus_ecalvetoDep);
   fChain->SetBranchAddress("mus_hcalvetoDep", &mus_hcalvetoDep, &b_mus_hcalvetoDep);
   fChain->SetBranchAddress("mus_calEnergyEm", &mus_calEnergyEm, &b_mus_calEnergyEm);
   fChain->SetBranchAddress("mus_calEnergyHad", &mus_calEnergyHad, &b_mus_calEnergyHad);
   fChain->SetBranchAddress("mus_calEnergyHo", &mus_calEnergyHo, &b_mus_calEnergyHo);
   fChain->SetBranchAddress("mus_calEnergyEmS9", &mus_calEnergyEmS9, &b_mus_calEnergyEmS9);
   fChain->SetBranchAddress("mus_calEnergyHadS9", &mus_calEnergyHadS9, &b_mus_calEnergyHadS9);
   fChain->SetBranchAddress("mus_calEnergyHoS9", &mus_calEnergyHoS9, &b_mus_calEnergyHoS9);
   fChain->SetBranchAddress("mus_iso03_emVetoEt", &mus_iso03_emVetoEt, &b_mus_iso03_emVetoEt);
   fChain->SetBranchAddress("mus_iso03_hadVetoEt", &mus_iso03_hadVetoEt, &b_mus_iso03_hadVetoEt);
   fChain->SetBranchAddress("mus_iso03_sumPt", &mus_iso03_sumPt, &b_mus_iso03_sumPt);
   fChain->SetBranchAddress("mus_iso03_emEt", &mus_iso03_emEt, &b_mus_iso03_emEt);
   fChain->SetBranchAddress("mus_iso03_hadEt", &mus_iso03_hadEt, &b_mus_iso03_hadEt);
   fChain->SetBranchAddress("mus_iso03_hoEt", &mus_iso03_hoEt, &b_mus_iso03_hoEt);
   fChain->SetBranchAddress("mus_iso03_nTracks", &mus_iso03_nTracks, &b_mus_iso03_nTracks);
   fChain->SetBranchAddress("mus_iso05_sumPt", &mus_iso05_sumPt, &b_mus_iso05_sumPt);
   fChain->SetBranchAddress("mus_iso05_emEt", &mus_iso05_emEt, &b_mus_iso05_emEt);
   fChain->SetBranchAddress("mus_iso05_hadEt", &mus_iso05_hadEt, &b_mus_iso05_hadEt);
   fChain->SetBranchAddress("mus_iso05_hoEt", &mus_iso05_hoEt, &b_mus_iso05_hoEt);
   fChain->SetBranchAddress("mus_iso05_nTracks", &mus_iso05_nTracks, &b_mus_iso05_nTracks);
   fChain->SetBranchAddress("mus_pfIsolationR03_sumChargedHadronPt", &mus_pfIsolationR03_sumChargedHadronPt, &b_mus_pfIsolationR03_sumChargedHadronPt);
   fChain->SetBranchAddress("mus_pfIsolationR03_sumChargedParticlePt", &mus_pfIsolationR03_sumChargedParticlePt, &b_mus_pfIsolationR03_sumChargedParticlePt);
   fChain->SetBranchAddress("mus_pfIsolationR03_sumNeutralHadronEt", &mus_pfIsolationR03_sumNeutralHadronEt, &b_mus_pfIsolationR03_sumNeutralHadronEt);
   fChain->SetBranchAddress("mus_pfIsolationR03_sumNeutralHadronEtHighThreshold", &mus_pfIsolationR03_sumNeutralHadronEtHighThreshold, &b_mus_pfIsolationR03_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("mus_pfIsolationR03_sumPhotonEt", &mus_pfIsolationR03_sumPhotonEt, &b_mus_pfIsolationR03_sumPhotonEt);
   fChain->SetBranchAddress("mus_pfIsolationR03_sumPhotonEtHighThreshold", &mus_pfIsolationR03_sumPhotonEtHighThreshold, &b_mus_pfIsolationR03_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("mus_pfIsolationR03_sumPUPt", &mus_pfIsolationR03_sumPUPt, &b_mus_pfIsolationR03_sumPUPt);
   fChain->SetBranchAddress("mus_pfIsolationR04_sumChargedHadronPt", &mus_pfIsolationR04_sumChargedHadronPt, &b_mus_pfIsolationR04_sumChargedHadronPt);
   fChain->SetBranchAddress("mus_pfIsolationR04_sumChargedParticlePt", &mus_pfIsolationR04_sumChargedParticlePt, &b_mus_pfIsolationR04_sumChargedParticlePt);
   fChain->SetBranchAddress("mus_pfIsolationR04_sumNeutralHadronEt", &mus_pfIsolationR04_sumNeutralHadronEt, &b_mus_pfIsolationR04_sumNeutralHadronEt);
   fChain->SetBranchAddress("mus_pfIsolationR04_sumNeutralHadronEtHighThreshold", &mus_pfIsolationR04_sumNeutralHadronEtHighThreshold, &b_mus_pfIsolationR04_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("mus_pfIsolationR04_sumPhotonEt", &mus_pfIsolationR04_sumPhotonEt, &b_mus_pfIsolationR04_sumPhotonEt);
   fChain->SetBranchAddress("mus_pfIsolationR04_sumPhotonEtHighThreshold", &mus_pfIsolationR04_sumPhotonEtHighThreshold, &b_mus_pfIsolationR04_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("mus_pfIsolationR04_sumPUPt", &mus_pfIsolationR04_sumPUPt, &b_mus_pfIsolationR04_sumPUPt);
   fChain->SetBranchAddress("mus_charge", &mus_charge, &b_mus_charge);
   fChain->SetBranchAddress("mus_cm_chi2", &mus_cm_chi2, &b_mus_cm_chi2);
   fChain->SetBranchAddress("mus_cm_ndof", &mus_cm_ndof, &b_mus_cm_ndof);
   fChain->SetBranchAddress("mus_cm_chg", &mus_cm_chg, &b_mus_cm_chg);
   fChain->SetBranchAddress("mus_cm_pt", &mus_cm_pt, &b_mus_cm_pt);
   fChain->SetBranchAddress("mus_cm_px", &mus_cm_px, &b_mus_cm_px);
   fChain->SetBranchAddress("mus_cm_py", &mus_cm_py, &b_mus_cm_py);
   fChain->SetBranchAddress("mus_cm_pz", &mus_cm_pz, &b_mus_cm_pz);
   fChain->SetBranchAddress("mus_cm_eta", &mus_cm_eta, &b_mus_cm_eta);
   fChain->SetBranchAddress("mus_cm_phi", &mus_cm_phi, &b_mus_cm_phi);
   fChain->SetBranchAddress("mus_cm_theta", &mus_cm_theta, &b_mus_cm_theta);
   fChain->SetBranchAddress("mus_cm_d0dum", &mus_cm_d0dum, &b_mus_cm_d0dum);
   fChain->SetBranchAddress("mus_cm_dz", &mus_cm_dz, &b_mus_cm_dz);
   fChain->SetBranchAddress("mus_cm_vx", &mus_cm_vx, &b_mus_cm_vx);
   fChain->SetBranchAddress("mus_cm_vy", &mus_cm_vy, &b_mus_cm_vy);
   fChain->SetBranchAddress("mus_cm_vz", &mus_cm_vz, &b_mus_cm_vz);
   fChain->SetBranchAddress("mus_cm_numvalhits", &mus_cm_numvalhits, &b_mus_cm_numvalhits);
   fChain->SetBranchAddress("mus_cm_numlosthits", &mus_cm_numlosthits, &b_mus_cm_numlosthits);
   fChain->SetBranchAddress("mus_cm_numvalMuonhits", &mus_cm_numvalMuonhits, &b_mus_cm_numvalMuonhits);
   fChain->SetBranchAddress("mus_cm_d0dumErr", &mus_cm_d0dumErr, &b_mus_cm_d0dumErr);
   fChain->SetBranchAddress("mus_cm_dzErr", &mus_cm_dzErr, &b_mus_cm_dzErr);
   fChain->SetBranchAddress("mus_cm_ptErr", &mus_cm_ptErr, &b_mus_cm_ptErr);
   fChain->SetBranchAddress("mus_cm_etaErr", &mus_cm_etaErr, &b_mus_cm_etaErr);
   fChain->SetBranchAddress("mus_cm_phiErr", &mus_cm_phiErr, &b_mus_cm_phiErr);
   fChain->SetBranchAddress("mus_tk_id", &mus_tk_id, &b_mus_tk_id);
   fChain->SetBranchAddress("mus_tk_chi2", &mus_tk_chi2, &b_mus_tk_chi2);
   fChain->SetBranchAddress("mus_tk_ndof", &mus_tk_ndof, &b_mus_tk_ndof);
   fChain->SetBranchAddress("mus_tk_chg", &mus_tk_chg, &b_mus_tk_chg);
   fChain->SetBranchAddress("mus_tk_pt", &mus_tk_pt, &b_mus_tk_pt);
   fChain->SetBranchAddress("mus_tk_px", &mus_tk_px, &b_mus_tk_px);
   fChain->SetBranchAddress("mus_tk_py", &mus_tk_py, &b_mus_tk_py);
   fChain->SetBranchAddress("mus_tk_pz", &mus_tk_pz, &b_mus_tk_pz);
   fChain->SetBranchAddress("mus_tk_eta", &mus_tk_eta, &b_mus_tk_eta);
   fChain->SetBranchAddress("mus_tk_phi", &mus_tk_phi, &b_mus_tk_phi);
   fChain->SetBranchAddress("mus_tk_theta", &mus_tk_theta, &b_mus_tk_theta);
   fChain->SetBranchAddress("mus_tk_d0dum", &mus_tk_d0dum, &b_mus_tk_d0dum);
   fChain->SetBranchAddress("mus_tk_dz", &mus_tk_dz, &b_mus_tk_dz);
   fChain->SetBranchAddress("mus_tk_vx", &mus_tk_vx, &b_mus_tk_vx);
   fChain->SetBranchAddress("mus_tk_vy", &mus_tk_vy, &b_mus_tk_vy);
   fChain->SetBranchAddress("mus_tk_vz", &mus_tk_vz, &b_mus_tk_vz);
   fChain->SetBranchAddress("mus_tk_numvalhits", &mus_tk_numvalhits, &b_mus_tk_numvalhits);
   fChain->SetBranchAddress("mus_tk_numlosthits", &mus_tk_numlosthits, &b_mus_tk_numlosthits);
   fChain->SetBranchAddress("mus_tk_d0dumErr", &mus_tk_d0dumErr, &b_mus_tk_d0dumErr);
   fChain->SetBranchAddress("mus_tk_dzErr", &mus_tk_dzErr, &b_mus_tk_dzErr);
   fChain->SetBranchAddress("mus_tk_ptErr", &mus_tk_ptErr, &b_mus_tk_ptErr);
   fChain->SetBranchAddress("mus_tk_etaErr", &mus_tk_etaErr, &b_mus_tk_etaErr);
   fChain->SetBranchAddress("mus_tk_phiErr", &mus_tk_phiErr, &b_mus_tk_phiErr);
   fChain->SetBranchAddress("mus_tk_numvalPixelhits", &mus_tk_numvalPixelhits, &b_mus_tk_numvalPixelhits);
   fChain->SetBranchAddress("mus_tk_numpixelWthMeasr", &mus_tk_numpixelWthMeasr, &b_mus_tk_numpixelWthMeasr);
   fChain->SetBranchAddress("mus_stamu_chi2", &mus_stamu_chi2, &b_mus_stamu_chi2);
   fChain->SetBranchAddress("mus_stamu_ndof", &mus_stamu_ndof, &b_mus_stamu_ndof);
   fChain->SetBranchAddress("mus_stamu_chg", &mus_stamu_chg, &b_mus_stamu_chg);
   fChain->SetBranchAddress("mus_stamu_pt", &mus_stamu_pt, &b_mus_stamu_pt);
   fChain->SetBranchAddress("mus_stamu_px", &mus_stamu_px, &b_mus_stamu_px);
   fChain->SetBranchAddress("mus_stamu_py", &mus_stamu_py, &b_mus_stamu_py);
   fChain->SetBranchAddress("mus_stamu_pz", &mus_stamu_pz, &b_mus_stamu_pz);
   fChain->SetBranchAddress("mus_stamu_eta", &mus_stamu_eta, &b_mus_stamu_eta);
   fChain->SetBranchAddress("mus_stamu_phi", &mus_stamu_phi, &b_mus_stamu_phi);
   fChain->SetBranchAddress("mus_stamu_theta", &mus_stamu_theta, &b_mus_stamu_theta);
   fChain->SetBranchAddress("mus_stamu_d0dum", &mus_stamu_d0dum, &b_mus_stamu_d0dum);
   fChain->SetBranchAddress("mus_stamu_dz", &mus_stamu_dz, &b_mus_stamu_dz);
   fChain->SetBranchAddress("mus_stamu_vx", &mus_stamu_vx, &b_mus_stamu_vx);
   fChain->SetBranchAddress("mus_stamu_vy", &mus_stamu_vy, &b_mus_stamu_vy);
   fChain->SetBranchAddress("mus_stamu_vz", &mus_stamu_vz, &b_mus_stamu_vz);
   fChain->SetBranchAddress("mus_stamu_numvalhits", &mus_stamu_numvalhits, &b_mus_stamu_numvalhits);
   fChain->SetBranchAddress("mus_stamu_numlosthits", &mus_stamu_numlosthits, &b_mus_stamu_numlosthits);
   fChain->SetBranchAddress("mus_stamu_d0dumErr", &mus_stamu_d0dumErr, &b_mus_stamu_d0dumErr);
   fChain->SetBranchAddress("mus_stamu_dzErr", &mus_stamu_dzErr, &b_mus_stamu_dzErr);
   fChain->SetBranchAddress("mus_stamu_ptErr", &mus_stamu_ptErr, &b_mus_stamu_ptErr);
   fChain->SetBranchAddress("mus_stamu_etaErr", &mus_stamu_etaErr, &b_mus_stamu_etaErr);
   fChain->SetBranchAddress("mus_stamu_phiErr", &mus_stamu_phiErr, &b_mus_stamu_phiErr);
   fChain->SetBranchAddress("mus_num_matches", &mus_num_matches, &b_mus_num_matches);
   fChain->SetBranchAddress("mus_isPFMuon", &mus_isPFMuon, &b_mus_isPFMuon);
   fChain->SetBranchAddress("mus_isTrackerMuon", &mus_isTrackerMuon, &b_mus_isTrackerMuon);
   fChain->SetBranchAddress("mus_isStandAloneMuon", &mus_isStandAloneMuon, &b_mus_isStandAloneMuon);
   fChain->SetBranchAddress("mus_isCaloMuon", &mus_isCaloMuon, &b_mus_isCaloMuon);
   fChain->SetBranchAddress("mus_isGlobalMuon", &mus_isGlobalMuon, &b_mus_isGlobalMuon);
   fChain->SetBranchAddress("mus_isElectron", &mus_isElectron, &b_mus_isElectron);
   fChain->SetBranchAddress("mus_isConvertedPhoton", &mus_isConvertedPhoton, &b_mus_isConvertedPhoton);
   fChain->SetBranchAddress("mus_isPhoton", &mus_isPhoton, &b_mus_isPhoton);
   fChain->SetBranchAddress("mus_id_All", &mus_id_All, &b_mus_id_All);
   fChain->SetBranchAddress("mus_id_AllGlobalMuons", &mus_id_AllGlobalMuons, &b_mus_id_AllGlobalMuons);
   fChain->SetBranchAddress("mus_id_AllStandAloneMuons", &mus_id_AllStandAloneMuons, &b_mus_id_AllStandAloneMuons);
   fChain->SetBranchAddress("mus_id_AllTrackerMuons", &mus_id_AllTrackerMuons, &b_mus_id_AllTrackerMuons);
   fChain->SetBranchAddress("mus_id_TrackerMuonArbitrated", &mus_id_TrackerMuonArbitrated, &b_mus_id_TrackerMuonArbitrated);
   fChain->SetBranchAddress("mus_id_AllArbitrated", &mus_id_AllArbitrated, &b_mus_id_AllArbitrated);
   fChain->SetBranchAddress("mus_id_GlobalMuonPromptTight", &mus_id_GlobalMuonPromptTight, &b_mus_id_GlobalMuonPromptTight);
   fChain->SetBranchAddress("mus_id_TMLastStationLoose", &mus_id_TMLastStationLoose, &b_mus_id_TMLastStationLoose);
   fChain->SetBranchAddress("mus_id_TMLastStationTight", &mus_id_TMLastStationTight, &b_mus_id_TMLastStationTight);
   fChain->SetBranchAddress("mus_id_TM2DCompatibilityLoose", &mus_id_TM2DCompatibilityLoose, &b_mus_id_TM2DCompatibilityLoose);
   fChain->SetBranchAddress("mus_id_TM2DCompatibilityTight", &mus_id_TM2DCompatibilityTight, &b_mus_id_TM2DCompatibilityTight);
   fChain->SetBranchAddress("mus_id_TMOneStationLoose", &mus_id_TMOneStationLoose, &b_mus_id_TMOneStationLoose);
   fChain->SetBranchAddress("mus_id_TMOneStationTight", &mus_id_TMOneStationTight, &b_mus_id_TMOneStationTight);
   fChain->SetBranchAddress("mus_id_TMLastStationOptimizedLowPtLoose", &mus_id_TMLastStationOptimizedLowPtLoose, &b_mus_id_TMLastStationOptimizedLowPtLoose);
   fChain->SetBranchAddress("mus_id_TMLastStationOptimizedLowPtTight", &mus_id_TMLastStationOptimizedLowPtTight, &b_mus_id_TMLastStationOptimizedLowPtTight);
   fChain->SetBranchAddress("mus_tk_LayersWithMeasurement", &mus_tk_LayersWithMeasurement, &b_mus_tk_LayersWithMeasurement);
   fChain->SetBranchAddress("mus_tk_PixelLayersWithMeasurement", &mus_tk_PixelLayersWithMeasurement, &b_mus_tk_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("mus_tk_ValidStripLayersWithMonoAndStereoHit", &mus_tk_ValidStripLayersWithMonoAndStereoHit, &b_mus_tk_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("mus_tk_LayersWithoutMeasurement", &mus_tk_LayersWithoutMeasurement, &b_mus_tk_LayersWithoutMeasurement);
   fChain->SetBranchAddress("mus_tk_ExpectedHitsInner", &mus_tk_ExpectedHitsInner, &b_mus_tk_ExpectedHitsInner);
   fChain->SetBranchAddress("mus_tk_ExpectedHitsOuter", &mus_tk_ExpectedHitsOuter, &b_mus_tk_ExpectedHitsOuter);
   fChain->SetBranchAddress("mus_cm_LayersWithMeasurement", &mus_cm_LayersWithMeasurement, &b_mus_cm_LayersWithMeasurement);
   fChain->SetBranchAddress("mus_cm_PixelLayersWithMeasurement", &mus_cm_PixelLayersWithMeasurement, &b_mus_cm_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("mus_cm_ValidStripLayersWithMonoAndStereoHit", &mus_cm_ValidStripLayersWithMonoAndStereoHit, &b_mus_cm_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("mus_cm_LayersWithoutMeasurement", &mus_cm_LayersWithoutMeasurement, &b_mus_cm_LayersWithoutMeasurement);
   fChain->SetBranchAddress("mus_cm_ExpectedHitsInner", &mus_cm_ExpectedHitsInner, &b_mus_cm_ExpectedHitsInner);
   fChain->SetBranchAddress("mus_cm_ExpectedHitsOuter", &mus_cm_ExpectedHitsOuter, &b_mus_cm_ExpectedHitsOuter);
   fChain->SetBranchAddress("mus_picky_LayersWithMeasurement", &mus_picky_LayersWithMeasurement, &b_mus_picky_LayersWithMeasurement);
   fChain->SetBranchAddress("mus_picky_PixelLayersWithMeasurement", &mus_picky_PixelLayersWithMeasurement, &b_mus_picky_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("mus_picky_ValidStripLayersWithMonoAndStereoHit", &mus_picky_ValidStripLayersWithMonoAndStereoHit, &b_mus_picky_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("mus_picky_LayersWithoutMeasurement", &mus_picky_LayersWithoutMeasurement, &b_mus_picky_LayersWithoutMeasurement);
   fChain->SetBranchAddress("mus_picky_ExpectedHitsInner", &mus_picky_ExpectedHitsInner, &b_mus_picky_ExpectedHitsInner);
   fChain->SetBranchAddress("mus_picky_ExpectedHitsOuter", &mus_picky_ExpectedHitsOuter, &b_mus_picky_ExpectedHitsOuter);
   fChain->SetBranchAddress("mus_tpfms_LayersWithMeasurement", &mus_tpfms_LayersWithMeasurement, &b_mus_tpfms_LayersWithMeasurement);
   fChain->SetBranchAddress("mus_tpfms_PixelLayersWithMeasurement", &mus_tpfms_PixelLayersWithMeasurement, &b_mus_tpfms_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("mus_tpfms_ValidStripLayersWithMonoAndStereoHit", &mus_tpfms_ValidStripLayersWithMonoAndStereoHit, &b_mus_tpfms_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("mus_tpfms_LayersWithoutMeasurement", &mus_tpfms_LayersWithoutMeasurement, &b_mus_tpfms_LayersWithoutMeasurement);
   fChain->SetBranchAddress("mus_tpfms_ExpectedHitsInner", &mus_tpfms_ExpectedHitsInner, &b_mus_tpfms_ExpectedHitsInner);
   fChain->SetBranchAddress("mus_tpfms_ExpectedHitsOuter", &mus_tpfms_ExpectedHitsOuter, &b_mus_tpfms_ExpectedHitsOuter);
   fChain->SetBranchAddress("mus_picky_id", &mus_picky_id, &b_mus_picky_id);
   fChain->SetBranchAddress("mus_picky_chi2", &mus_picky_chi2, &b_mus_picky_chi2);
   fChain->SetBranchAddress("mus_picky_ndof", &mus_picky_ndof, &b_mus_picky_ndof);
   fChain->SetBranchAddress("mus_picky_chg", &mus_picky_chg, &b_mus_picky_chg);
   fChain->SetBranchAddress("mus_picky_pt", &mus_picky_pt, &b_mus_picky_pt);
   fChain->SetBranchAddress("mus_picky_px", &mus_picky_px, &b_mus_picky_px);
   fChain->SetBranchAddress("mus_picky_py", &mus_picky_py, &b_mus_picky_py);
   fChain->SetBranchAddress("mus_picky_pz", &mus_picky_pz, &b_mus_picky_pz);
   fChain->SetBranchAddress("mus_picky_eta", &mus_picky_eta, &b_mus_picky_eta);
   fChain->SetBranchAddress("mus_picky_phi", &mus_picky_phi, &b_mus_picky_phi);
   fChain->SetBranchAddress("mus_picky_theta", &mus_picky_theta, &b_mus_picky_theta);
   fChain->SetBranchAddress("mus_picky_d0dum", &mus_picky_d0dum, &b_mus_picky_d0dum);
   fChain->SetBranchAddress("mus_picky_dz", &mus_picky_dz, &b_mus_picky_dz);
   fChain->SetBranchAddress("mus_picky_vx", &mus_picky_vx, &b_mus_picky_vx);
   fChain->SetBranchAddress("mus_picky_vy", &mus_picky_vy, &b_mus_picky_vy);
   fChain->SetBranchAddress("mus_picky_vz", &mus_picky_vz, &b_mus_picky_vz);
   fChain->SetBranchAddress("mus_picky_numvalhits", &mus_picky_numvalhits, &b_mus_picky_numvalhits);
   fChain->SetBranchAddress("mus_picky_numlosthits", &mus_picky_numlosthits, &b_mus_picky_numlosthits);
   fChain->SetBranchAddress("mus_picky_d0dumErr", &mus_picky_d0dumErr, &b_mus_picky_d0dumErr);
   fChain->SetBranchAddress("mus_picky_dzErr", &mus_picky_dzErr, &b_mus_picky_dzErr);
   fChain->SetBranchAddress("mus_picky_ptErr", &mus_picky_ptErr, &b_mus_picky_ptErr);
   fChain->SetBranchAddress("mus_picky_etaErr", &mus_picky_etaErr, &b_mus_picky_etaErr);
   fChain->SetBranchAddress("mus_picky_phiErr", &mus_picky_phiErr, &b_mus_picky_phiErr);
   fChain->SetBranchAddress("mus_picky_numvalPixelhits", &mus_picky_numvalPixelhits, &b_mus_picky_numvalPixelhits);
   fChain->SetBranchAddress("mus_tpfms_id", &mus_tpfms_id, &b_mus_tpfms_id);
   fChain->SetBranchAddress("mus_tpfms_chi2", &mus_tpfms_chi2, &b_mus_tpfms_chi2);
   fChain->SetBranchAddress("mus_tpfms_ndof", &mus_tpfms_ndof, &b_mus_tpfms_ndof);
   fChain->SetBranchAddress("mus_tpfms_chg", &mus_tpfms_chg, &b_mus_tpfms_chg);
   fChain->SetBranchAddress("mus_tpfms_pt", &mus_tpfms_pt, &b_mus_tpfms_pt);
   fChain->SetBranchAddress("mus_tpfms_px", &mus_tpfms_px, &b_mus_tpfms_px);
   fChain->SetBranchAddress("mus_tpfms_py", &mus_tpfms_py, &b_mus_tpfms_py);
   fChain->SetBranchAddress("mus_tpfms_pz", &mus_tpfms_pz, &b_mus_tpfms_pz);
   fChain->SetBranchAddress("mus_tpfms_eta", &mus_tpfms_eta, &b_mus_tpfms_eta);
   fChain->SetBranchAddress("mus_tpfms_phi", &mus_tpfms_phi, &b_mus_tpfms_phi);
   fChain->SetBranchAddress("mus_tpfms_theta", &mus_tpfms_theta, &b_mus_tpfms_theta);
   fChain->SetBranchAddress("mus_tpfms_d0dum", &mus_tpfms_d0dum, &b_mus_tpfms_d0dum);
   fChain->SetBranchAddress("mus_tpfms_dz", &mus_tpfms_dz, &b_mus_tpfms_dz);
   fChain->SetBranchAddress("mus_tpfms_vx", &mus_tpfms_vx, &b_mus_tpfms_vx);
   fChain->SetBranchAddress("mus_tpfms_vy", &mus_tpfms_vy, &b_mus_tpfms_vy);
   fChain->SetBranchAddress("mus_tpfms_vz", &mus_tpfms_vz, &b_mus_tpfms_vz);
   fChain->SetBranchAddress("mus_tpfms_numvalhits", &mus_tpfms_numvalhits, &b_mus_tpfms_numvalhits);
   fChain->SetBranchAddress("mus_tpfms_numlosthits", &mus_tpfms_numlosthits, &b_mus_tpfms_numlosthits);
   fChain->SetBranchAddress("mus_tpfms_d0dumErr", &mus_tpfms_d0dumErr, &b_mus_tpfms_d0dumErr);
   fChain->SetBranchAddress("mus_tpfms_dzErr", &mus_tpfms_dzErr, &b_mus_tpfms_dzErr);
   fChain->SetBranchAddress("mus_tpfms_ptErr", &mus_tpfms_ptErr, &b_mus_tpfms_ptErr);
   fChain->SetBranchAddress("mus_tpfms_etaErr", &mus_tpfms_etaErr, &b_mus_tpfms_etaErr);
   fChain->SetBranchAddress("mus_tpfms_phiErr", &mus_tpfms_phiErr, &b_mus_tpfms_phiErr);
   fChain->SetBranchAddress("mus_tpfms_numvalPixelhits", &mus_tpfms_numvalPixelhits, &b_mus_tpfms_numvalPixelhits);
   fChain->SetBranchAddress("mus_dB", &mus_dB, &b_mus_dB);
   fChain->SetBranchAddress("mus_numberOfMatchedStations", &mus_numberOfMatchedStations, &b_mus_numberOfMatchedStations);
   fChain->SetBranchAddress("NpfTypeINoXYCorrmets", &NpfTypeINoXYCorrmets, &b_NpfTypeINoXYCorrmets);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_et", &pfTypeINoXYCorrmets_et, &b_pfTypeINoXYCorrmets_et);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_phi", &pfTypeINoXYCorrmets_phi, &b_pfTypeINoXYCorrmets_phi);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_ex", &pfTypeINoXYCorrmets_ex, &b_pfTypeINoXYCorrmets_ex);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_ey", &pfTypeINoXYCorrmets_ey, &b_pfTypeINoXYCorrmets_ey);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_gen_et", &pfTypeINoXYCorrmets_gen_et, &b_pfTypeINoXYCorrmets_gen_et);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_gen_phi", &pfTypeINoXYCorrmets_gen_phi, &b_pfTypeINoXYCorrmets_gen_phi);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_sign", &pfTypeINoXYCorrmets_sign, &b_pfTypeINoXYCorrmets_sign);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_sumEt", &pfTypeINoXYCorrmets_sumEt, &b_pfTypeINoXYCorrmets_sumEt);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_unCPhi", &pfTypeINoXYCorrmets_unCPhi, &b_pfTypeINoXYCorrmets_unCPhi);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_unCPt", &pfTypeINoXYCorrmets_unCPt, &b_pfTypeINoXYCorrmets_unCPt);
   fChain->SetBranchAddress("NpfTypeIType0mets", &NpfTypeIType0mets, &b_NpfTypeIType0mets);
   fChain->SetBranchAddress("pfTypeIType0mets_et", &pfTypeIType0mets_et, &b_pfTypeIType0mets_et);
   fChain->SetBranchAddress("pfTypeIType0mets_phi", &pfTypeIType0mets_phi, &b_pfTypeIType0mets_phi);
   fChain->SetBranchAddress("pfTypeIType0mets_ex", &pfTypeIType0mets_ex, &b_pfTypeIType0mets_ex);
   fChain->SetBranchAddress("pfTypeIType0mets_ey", &pfTypeIType0mets_ey, &b_pfTypeIType0mets_ey);
   fChain->SetBranchAddress("pfTypeIType0mets_gen_et", &pfTypeIType0mets_gen_et, &b_pfTypeIType0mets_gen_et);
   fChain->SetBranchAddress("pfTypeIType0mets_gen_phi", &pfTypeIType0mets_gen_phi, &b_pfTypeIType0mets_gen_phi);
   fChain->SetBranchAddress("pfTypeIType0mets_sign", &pfTypeIType0mets_sign, &b_pfTypeIType0mets_sign);
   fChain->SetBranchAddress("pfTypeIType0mets_sumEt", &pfTypeIType0mets_sumEt, &b_pfTypeIType0mets_sumEt);
   fChain->SetBranchAddress("pfTypeIType0mets_unCPhi", &pfTypeIType0mets_unCPhi, &b_pfTypeIType0mets_unCPhi);
   fChain->SetBranchAddress("pfTypeIType0mets_unCPt", &pfTypeIType0mets_unCPt, &b_pfTypeIType0mets_unCPt);
   fChain->SetBranchAddress("NpfTypeImets", &NpfTypeImets, &b_NpfTypeImets);
   fChain->SetBranchAddress("pfTypeImets_et", &pfTypeImets_et, &b_pfTypeImets_et);
   fChain->SetBranchAddress("pfTypeImets_phi", &pfTypeImets_phi, &b_pfTypeImets_phi);
   fChain->SetBranchAddress("pfTypeImets_ex", &pfTypeImets_ex, &b_pfTypeImets_ex);
   fChain->SetBranchAddress("pfTypeImets_ey", &pfTypeImets_ey, &b_pfTypeImets_ey);
   fChain->SetBranchAddress("pfTypeImets_gen_et", &pfTypeImets_gen_et, &b_pfTypeImets_gen_et);
   fChain->SetBranchAddress("pfTypeImets_gen_phi", &pfTypeImets_gen_phi, &b_pfTypeImets_gen_phi);
   fChain->SetBranchAddress("pfTypeImets_sign", &pfTypeImets_sign, &b_pfTypeImets_sign);
   fChain->SetBranchAddress("pfTypeImets_sumEt", &pfTypeImets_sumEt, &b_pfTypeImets_sumEt);
   fChain->SetBranchAddress("pfTypeImets_unCPhi", &pfTypeImets_unCPhi, &b_pfTypeImets_unCPhi);
   fChain->SetBranchAddress("pfTypeImets_unCPt", &pfTypeImets_unCPt, &b_pfTypeImets_unCPt);
   fChain->SetBranchAddress("Npf_els", &Npf_els, &b_Npf_els);
   fChain->SetBranchAddress("pf_els_energy", &pf_els_energy, &b_pf_els_energy);
   fChain->SetBranchAddress("pf_els_et", &pf_els_et, &b_pf_els_et);
   fChain->SetBranchAddress("pf_els_eta", &pf_els_eta, &b_pf_els_eta);
   fChain->SetBranchAddress("pf_els_phi", &pf_els_phi, &b_pf_els_phi);
   fChain->SetBranchAddress("pf_els_pt", &pf_els_pt, &b_pf_els_pt);
   fChain->SetBranchAddress("pf_els_px", &pf_els_px, &b_pf_els_px);
   fChain->SetBranchAddress("pf_els_py", &pf_els_py, &b_pf_els_py);
   fChain->SetBranchAddress("pf_els_pz", &pf_els_pz, &b_pf_els_pz);
   fChain->SetBranchAddress("pf_els_status", &pf_els_status, &b_pf_els_status);
   fChain->SetBranchAddress("pf_els_theta", &pf_els_theta, &b_pf_els_theta);
   fChain->SetBranchAddress("pf_els_gen_id", &pf_els_gen_id, &b_pf_els_gen_id);
   fChain->SetBranchAddress("pf_els_gen_phi", &pf_els_gen_phi, &b_pf_els_gen_phi);
   fChain->SetBranchAddress("pf_els_gen_pt", &pf_els_gen_pt, &b_pf_els_gen_pt);
   fChain->SetBranchAddress("pf_els_gen_pz", &pf_els_gen_pz, &b_pf_els_gen_pz);
   fChain->SetBranchAddress("pf_els_gen_px", &pf_els_gen_px, &b_pf_els_gen_px);
   fChain->SetBranchAddress("pf_els_gen_py", &pf_els_gen_py, &b_pf_els_gen_py);
   fChain->SetBranchAddress("pf_els_gen_eta", &pf_els_gen_eta, &b_pf_els_gen_eta);
   fChain->SetBranchAddress("pf_els_gen_theta", &pf_els_gen_theta, &b_pf_els_gen_theta);
   fChain->SetBranchAddress("pf_els_gen_et", &pf_els_gen_et, &b_pf_els_gen_et);
   fChain->SetBranchAddress("pf_els_gen_mother_id", &pf_els_gen_mother_id, &b_pf_els_gen_mother_id);
   fChain->SetBranchAddress("pf_els_gen_mother_phi", &pf_els_gen_mother_phi, &b_pf_els_gen_mother_phi);
   fChain->SetBranchAddress("pf_els_gen_mother_pt", &pf_els_gen_mother_pt, &b_pf_els_gen_mother_pt);
   fChain->SetBranchAddress("pf_els_gen_mother_pz", &pf_els_gen_mother_pz, &b_pf_els_gen_mother_pz);
   fChain->SetBranchAddress("pf_els_gen_mother_px", &pf_els_gen_mother_px, &b_pf_els_gen_mother_px);
   fChain->SetBranchAddress("pf_els_gen_mother_py", &pf_els_gen_mother_py, &b_pf_els_gen_mother_py);
   fChain->SetBranchAddress("pf_els_gen_mother_eta", &pf_els_gen_mother_eta, &b_pf_els_gen_mother_eta);
   fChain->SetBranchAddress("pf_els_gen_mother_theta", &pf_els_gen_mother_theta, &b_pf_els_gen_mother_theta);
   fChain->SetBranchAddress("pf_els_gen_mother_et", &pf_els_gen_mother_et, &b_pf_els_gen_mother_et);
   fChain->SetBranchAddress("pf_els_tightId", &pf_els_tightId, &b_pf_els_tightId);
   fChain->SetBranchAddress("pf_els_looseId", &pf_els_looseId, &b_pf_els_looseId);
   fChain->SetBranchAddress("pf_els_robustTightId", &pf_els_robustTightId, &b_pf_els_robustTightId);
   fChain->SetBranchAddress("pf_els_robustLooseId", &pf_els_robustLooseId, &b_pf_els_robustLooseId);
   fChain->SetBranchAddress("pf_els_robustHighEnergyId", &pf_els_robustHighEnergyId, &b_pf_els_robustHighEnergyId);
   fChain->SetBranchAddress("pf_els_simpleEleId95relIso", &pf_els_simpleEleId95relIso, &b_pf_els_simpleEleId95relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId90relIso", &pf_els_simpleEleId90relIso, &b_pf_els_simpleEleId90relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId85relIso", &pf_els_simpleEleId85relIso, &b_pf_els_simpleEleId85relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId80relIso", &pf_els_simpleEleId80relIso, &b_pf_els_simpleEleId80relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId70relIso", &pf_els_simpleEleId70relIso, &b_pf_els_simpleEleId70relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId60relIso", &pf_els_simpleEleId60relIso, &b_pf_els_simpleEleId60relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId95cIso", &pf_els_simpleEleId95cIso, &b_pf_els_simpleEleId95cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId90cIso", &pf_els_simpleEleId90cIso, &b_pf_els_simpleEleId90cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId85cIso", &pf_els_simpleEleId85cIso, &b_pf_els_simpleEleId85cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId80cIso", &pf_els_simpleEleId80cIso, &b_pf_els_simpleEleId80cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId70cIso", &pf_els_simpleEleId70cIso, &b_pf_els_simpleEleId70cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId60cIso", &pf_els_simpleEleId60cIso, &b_pf_els_simpleEleId60cIso);
   fChain->SetBranchAddress("pf_els_cIso", &pf_els_cIso, &b_pf_els_cIso);
   fChain->SetBranchAddress("pf_els_tIso", &pf_els_tIso, &b_pf_els_tIso);
   fChain->SetBranchAddress("pf_els_ecalIso", &pf_els_ecalIso, &b_pf_els_ecalIso);
   fChain->SetBranchAddress("pf_els_hcalIso", &pf_els_hcalIso, &b_pf_els_hcalIso);
   fChain->SetBranchAddress("pf_els_chargedHadronIso", &pf_els_chargedHadronIso, &b_pf_els_chargedHadronIso);
   fChain->SetBranchAddress("pf_els_photonIso", &pf_els_photonIso, &b_pf_els_photonIso);
   fChain->SetBranchAddress("pf_els_neutralHadronIso", &pf_els_neutralHadronIso, &b_pf_els_neutralHadronIso);
   fChain->SetBranchAddress("pf_els_chi2", &pf_els_chi2, &b_pf_els_chi2);
   fChain->SetBranchAddress("pf_els_charge", &pf_els_charge, &b_pf_els_charge);
   fChain->SetBranchAddress("pf_els_caloEnergy", &pf_els_caloEnergy, &b_pf_els_caloEnergy);
   fChain->SetBranchAddress("pf_els_hadOverEm", &pf_els_hadOverEm, &b_pf_els_hadOverEm);
   fChain->SetBranchAddress("pf_els_hcalOverEcalBc", &pf_els_hcalOverEcalBc, &b_pf_els_hcalOverEcalBc);
   fChain->SetBranchAddress("pf_els_eOverPIn", &pf_els_eOverPIn, &b_pf_els_eOverPIn);
   fChain->SetBranchAddress("pf_els_eSeedOverPOut", &pf_els_eSeedOverPOut, &b_pf_els_eSeedOverPOut);
   fChain->SetBranchAddress("pf_els_sigmaEtaEta", &pf_els_sigmaEtaEta, &b_pf_els_sigmaEtaEta);
   fChain->SetBranchAddress("pf_els_sigmaIEtaIEta", &pf_els_sigmaIEtaIEta, &b_pf_els_sigmaIEtaIEta);
   fChain->SetBranchAddress("pf_els_scEnergy", &pf_els_scEnergy, &b_pf_els_scEnergy);
   fChain->SetBranchAddress("pf_els_scRawEnergy", &pf_els_scRawEnergy, &b_pf_els_scRawEnergy);
   fChain->SetBranchAddress("pf_els_scSeedEnergy", &pf_els_scSeedEnergy, &b_pf_els_scSeedEnergy);
   fChain->SetBranchAddress("pf_els_scEta", &pf_els_scEta, &b_pf_els_scEta);
   fChain->SetBranchAddress("pf_els_scPhi", &pf_els_scPhi, &b_pf_els_scPhi);
   fChain->SetBranchAddress("pf_els_scEtaWidth", &pf_els_scEtaWidth, &b_pf_els_scEtaWidth);
   fChain->SetBranchAddress("pf_els_scPhiWidth", &pf_els_scPhiWidth, &b_pf_els_scPhiWidth);
   fChain->SetBranchAddress("pf_els_scE1x5", &pf_els_scE1x5, &b_pf_els_scE1x5);
   fChain->SetBranchAddress("pf_els_scE2x5Max", &pf_els_scE2x5Max, &b_pf_els_scE2x5Max);
   fChain->SetBranchAddress("pf_els_scE5x5", &pf_els_scE5x5, &b_pf_els_scE5x5);
   fChain->SetBranchAddress("pf_els_isEB", &pf_els_isEB, &b_pf_els_isEB);
   fChain->SetBranchAddress("pf_els_isEE", &pf_els_isEE, &b_pf_els_isEE);
   fChain->SetBranchAddress("pf_els_dEtaIn", &pf_els_dEtaIn, &b_pf_els_dEtaIn);
   fChain->SetBranchAddress("pf_els_dPhiIn", &pf_els_dPhiIn, &b_pf_els_dPhiIn);
   fChain->SetBranchAddress("pf_els_dEtaOut", &pf_els_dEtaOut, &b_pf_els_dEtaOut);
   fChain->SetBranchAddress("pf_els_dPhiOut", &pf_els_dPhiOut, &b_pf_els_dPhiOut);
   fChain->SetBranchAddress("pf_els_numvalhits", &pf_els_numvalhits, &b_pf_els_numvalhits);
   fChain->SetBranchAddress("pf_els_numlosthits", &pf_els_numlosthits, &b_pf_els_numlosthits);
   fChain->SetBranchAddress("pf_els_basicClustersSize", &pf_els_basicClustersSize, &b_pf_els_basicClustersSize);
   fChain->SetBranchAddress("pf_els_tk_pz", &pf_els_tk_pz, &b_pf_els_tk_pz);
   fChain->SetBranchAddress("pf_els_tk_pt", &pf_els_tk_pt, &b_pf_els_tk_pt);
   fChain->SetBranchAddress("pf_els_tk_phi", &pf_els_tk_phi, &b_pf_els_tk_phi);
   fChain->SetBranchAddress("pf_els_tk_eta", &pf_els_tk_eta, &b_pf_els_tk_eta);
   fChain->SetBranchAddress("pf_els_d0dum", &pf_els_d0dum, &b_pf_els_d0dum);
   fChain->SetBranchAddress("pf_els_dz", &pf_els_dz, &b_pf_els_dz);
   fChain->SetBranchAddress("pf_els_vx", &pf_els_vx, &b_pf_els_vx);
   fChain->SetBranchAddress("pf_els_vy", &pf_els_vy, &b_pf_els_vy);
   fChain->SetBranchAddress("pf_els_vz", &pf_els_vz, &b_pf_els_vz);
   fChain->SetBranchAddress("pf_els_ndof", &pf_els_ndof, &b_pf_els_ndof);
   fChain->SetBranchAddress("pf_els_ptError", &pf_els_ptError, &b_pf_els_ptError);
   fChain->SetBranchAddress("pf_els_d0dumError", &pf_els_d0dumError, &b_pf_els_d0dumError);
   fChain->SetBranchAddress("pf_els_dzError", &pf_els_dzError, &b_pf_els_dzError);
   fChain->SetBranchAddress("pf_els_etaError", &pf_els_etaError, &b_pf_els_etaError);
   fChain->SetBranchAddress("pf_els_phiError", &pf_els_phiError, &b_pf_els_phiError);
   fChain->SetBranchAddress("pf_els_tk_charge", &pf_els_tk_charge, &b_pf_els_tk_charge);
   fChain->SetBranchAddress("pf_els_core_ecalDrivenSeed", &pf_els_core_ecalDrivenSeed, &b_pf_els_core_ecalDrivenSeed);
   fChain->SetBranchAddress("pf_els_n_inner_layer", &pf_els_n_inner_layer, &b_pf_els_n_inner_layer);
   fChain->SetBranchAddress("pf_els_n_outer_layer", &pf_els_n_outer_layer, &b_pf_els_n_outer_layer);
   fChain->SetBranchAddress("pf_els_ctf_tk_id", &pf_els_ctf_tk_id, &b_pf_els_ctf_tk_id);
   fChain->SetBranchAddress("pf_els_ctf_tk_charge", &pf_els_ctf_tk_charge, &b_pf_els_ctf_tk_charge);
   fChain->SetBranchAddress("pf_els_ctf_tk_eta", &pf_els_ctf_tk_eta, &b_pf_els_ctf_tk_eta);
   fChain->SetBranchAddress("pf_els_ctf_tk_phi", &pf_els_ctf_tk_phi, &b_pf_els_ctf_tk_phi);
   fChain->SetBranchAddress("pf_els_fbrem", &pf_els_fbrem, &b_pf_els_fbrem);
   fChain->SetBranchAddress("pf_els_shFracInnerHits", &pf_els_shFracInnerHits, &b_pf_els_shFracInnerHits);
   fChain->SetBranchAddress("pf_els_dr03EcalRecHitSumEt", &pf_els_dr03EcalRecHitSumEt, &b_pf_els_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("pf_els_dr03HcalTowerSumEt", &pf_els_dr03HcalTowerSumEt, &b_pf_els_dr03HcalTowerSumEt);
   fChain->SetBranchAddress("pf_els_dr03HcalDepth1TowerSumEt", &pf_els_dr03HcalDepth1TowerSumEt, &b_pf_els_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("pf_els_dr03HcalDepth2TowerSumEt", &pf_els_dr03HcalDepth2TowerSumEt, &b_pf_els_dr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("pf_els_dr03TkSumPt", &pf_els_dr03TkSumPt, &b_pf_els_dr03TkSumPt);
   fChain->SetBranchAddress("pf_els_dr04EcalRecHitSumEt", &pf_els_dr04EcalRecHitSumEt, &b_pf_els_dr04EcalRecHitSumEt);
   fChain->SetBranchAddress("pf_els_dr04HcalTowerSumEt", &pf_els_dr04HcalTowerSumEt, &b_pf_els_dr04HcalTowerSumEt);
   fChain->SetBranchAddress("pf_els_dr04HcalDepth1TowerSumEt", &pf_els_dr04HcalDepth1TowerSumEt, &b_pf_els_dr04HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("pf_els_dr04HcalDepth2TowerSumEt", &pf_els_dr04HcalDepth2TowerSumEt, &b_pf_els_dr04HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("pf_els_dr04TkSumPt", &pf_els_dr04TkSumPt, &b_pf_els_dr04TkSumPt);
   fChain->SetBranchAddress("pf_els_cpx", &pf_els_cpx, &b_pf_els_cpx);
   fChain->SetBranchAddress("pf_els_cpy", &pf_els_cpy, &b_pf_els_cpy);
   fChain->SetBranchAddress("pf_els_cpz", &pf_els_cpz, &b_pf_els_cpz);
   fChain->SetBranchAddress("pf_els_vpx", &pf_els_vpx, &b_pf_els_vpx);
   fChain->SetBranchAddress("pf_els_vpy", &pf_els_vpy, &b_pf_els_vpy);
   fChain->SetBranchAddress("pf_els_vpz", &pf_els_vpz, &b_pf_els_vpz);
   fChain->SetBranchAddress("pf_els_cx", &pf_els_cx, &b_pf_els_cx);
   fChain->SetBranchAddress("pf_els_cy", &pf_els_cy, &b_pf_els_cy);
   fChain->SetBranchAddress("pf_els_cz", &pf_els_cz, &b_pf_els_cz);
   fChain->SetBranchAddress("pf_els_PATpassConversionVeto", &pf_els_PATpassConversionVeto, &b_pf_els_PATpassConversionVeto);
   fChain->SetBranchAddress("Npf_mus", &Npf_mus, &b_Npf_mus);
   fChain->SetBranchAddress("pf_mus_energy", &pf_mus_energy, &b_pf_mus_energy);
   fChain->SetBranchAddress("pf_mus_et", &pf_mus_et, &b_pf_mus_et);
   fChain->SetBranchAddress("pf_mus_eta", &pf_mus_eta, &b_pf_mus_eta);
   fChain->SetBranchAddress("pf_mus_phi", &pf_mus_phi, &b_pf_mus_phi);
   fChain->SetBranchAddress("pf_mus_pt", &pf_mus_pt, &b_pf_mus_pt);
   fChain->SetBranchAddress("pf_mus_px", &pf_mus_px, &b_pf_mus_px);
   fChain->SetBranchAddress("pf_mus_py", &pf_mus_py, &b_pf_mus_py);
   fChain->SetBranchAddress("pf_mus_pz", &pf_mus_pz, &b_pf_mus_pz);
   fChain->SetBranchAddress("pf_mus_status", &pf_mus_status, &b_pf_mus_status);
   fChain->SetBranchAddress("pf_mus_theta", &pf_mus_theta, &b_pf_mus_theta);
   fChain->SetBranchAddress("pf_mus_gen_id", &pf_mus_gen_id, &b_pf_mus_gen_id);
   fChain->SetBranchAddress("pf_mus_gen_phi", &pf_mus_gen_phi, &b_pf_mus_gen_phi);
   fChain->SetBranchAddress("pf_mus_gen_pt", &pf_mus_gen_pt, &b_pf_mus_gen_pt);
   fChain->SetBranchAddress("pf_mus_gen_pz", &pf_mus_gen_pz, &b_pf_mus_gen_pz);
   fChain->SetBranchAddress("pf_mus_gen_px", &pf_mus_gen_px, &b_pf_mus_gen_px);
   fChain->SetBranchAddress("pf_mus_gen_py", &pf_mus_gen_py, &b_pf_mus_gen_py);
   fChain->SetBranchAddress("pf_mus_gen_eta", &pf_mus_gen_eta, &b_pf_mus_gen_eta);
   fChain->SetBranchAddress("pf_mus_gen_theta", &pf_mus_gen_theta, &b_pf_mus_gen_theta);
   fChain->SetBranchAddress("pf_mus_gen_et", &pf_mus_gen_et, &b_pf_mus_gen_et);
   fChain->SetBranchAddress("pf_mus_gen_mother_id", &pf_mus_gen_mother_id, &b_pf_mus_gen_mother_id);
   fChain->SetBranchAddress("pf_mus_gen_mother_phi", &pf_mus_gen_mother_phi, &b_pf_mus_gen_mother_phi);
   fChain->SetBranchAddress("pf_mus_gen_mother_pt", &pf_mus_gen_mother_pt, &b_pf_mus_gen_mother_pt);
   fChain->SetBranchAddress("pf_mus_gen_mother_pz", &pf_mus_gen_mother_pz, &b_pf_mus_gen_mother_pz);
   fChain->SetBranchAddress("pf_mus_gen_mother_px", &pf_mus_gen_mother_px, &b_pf_mus_gen_mother_px);
   fChain->SetBranchAddress("pf_mus_gen_mother_py", &pf_mus_gen_mother_py, &b_pf_mus_gen_mother_py);
   fChain->SetBranchAddress("pf_mus_gen_mother_eta", &pf_mus_gen_mother_eta, &b_pf_mus_gen_mother_eta);
   fChain->SetBranchAddress("pf_mus_gen_mother_theta", &pf_mus_gen_mother_theta, &b_pf_mus_gen_mother_theta);
   fChain->SetBranchAddress("pf_mus_gen_mother_et", &pf_mus_gen_mother_et, &b_pf_mus_gen_mother_et);
   fChain->SetBranchAddress("pf_mus_tkHits", &pf_mus_tkHits, &b_pf_mus_tkHits);
   fChain->SetBranchAddress("pf_mus_cIso", &pf_mus_cIso, &b_pf_mus_cIso);
   fChain->SetBranchAddress("pf_mus_tIso", &pf_mus_tIso, &b_pf_mus_tIso);
   fChain->SetBranchAddress("pf_mus_ecalIso", &pf_mus_ecalIso, &b_pf_mus_ecalIso);
   fChain->SetBranchAddress("pf_mus_hcalIso", &pf_mus_hcalIso, &b_pf_mus_hcalIso);
   fChain->SetBranchAddress("pf_mus_iso03_emVetoEt", &pf_mus_iso03_emVetoEt, &b_pf_mus_iso03_emVetoEt);
   fChain->SetBranchAddress("pf_mus_iso03_hadVetoEt", &pf_mus_iso03_hadVetoEt, &b_pf_mus_iso03_hadVetoEt);
   fChain->SetBranchAddress("pf_mus_calEnergyEm", &pf_mus_calEnergyEm, &b_pf_mus_calEnergyEm);
   fChain->SetBranchAddress("pf_mus_calEnergyHad", &pf_mus_calEnergyHad, &b_pf_mus_calEnergyHad);
   fChain->SetBranchAddress("pf_mus_calEnergyHo", &pf_mus_calEnergyHo, &b_pf_mus_calEnergyHo);
   fChain->SetBranchAddress("pf_mus_calEnergyEmS9", &pf_mus_calEnergyEmS9, &b_pf_mus_calEnergyEmS9);
   fChain->SetBranchAddress("pf_mus_calEnergyHadS9", &pf_mus_calEnergyHadS9, &b_pf_mus_calEnergyHadS9);
   fChain->SetBranchAddress("pf_mus_calEnergyHoS9", &pf_mus_calEnergyHoS9, &b_pf_mus_calEnergyHoS9);
   fChain->SetBranchAddress("pf_mus_iso03_sumPt", &pf_mus_iso03_sumPt, &b_pf_mus_iso03_sumPt);
   fChain->SetBranchAddress("pf_mus_iso03_emEt", &pf_mus_iso03_emEt, &b_pf_mus_iso03_emEt);
   fChain->SetBranchAddress("pf_mus_iso03_hadEt", &pf_mus_iso03_hadEt, &b_pf_mus_iso03_hadEt);
   fChain->SetBranchAddress("pf_mus_iso03_hoEt", &pf_mus_iso03_hoEt, &b_pf_mus_iso03_hoEt);
   fChain->SetBranchAddress("pf_mus_iso03_nTracks", &pf_mus_iso03_nTracks, &b_pf_mus_iso03_nTracks);
   fChain->SetBranchAddress("pf_mus_iso05_sumPt", &pf_mus_iso05_sumPt, &b_pf_mus_iso05_sumPt);
   fChain->SetBranchAddress("pf_mus_iso05_emEt", &pf_mus_iso05_emEt, &b_pf_mus_iso05_emEt);
   fChain->SetBranchAddress("pf_mus_iso05_hadEt", &pf_mus_iso05_hadEt, &b_pf_mus_iso05_hadEt);
   fChain->SetBranchAddress("pf_mus_iso05_hoEt", &pf_mus_iso05_hoEt, &b_pf_mus_iso05_hoEt);
   fChain->SetBranchAddress("pf_mus_iso05_nTracks", &pf_mus_iso05_nTracks, &b_pf_mus_iso05_nTracks);
   fChain->SetBranchAddress("pf_mus_neutralHadronIso", &pf_mus_neutralHadronIso, &b_pf_mus_neutralHadronIso);
   fChain->SetBranchAddress("pf_mus_chargedHadronIso", &pf_mus_chargedHadronIso, &b_pf_mus_chargedHadronIso);
   fChain->SetBranchAddress("pf_mus_photonIso", &pf_mus_photonIso, &b_pf_mus_photonIso);
   fChain->SetBranchAddress("pf_mus_charge", &pf_mus_charge, &b_pf_mus_charge);
   fChain->SetBranchAddress("pf_mus_cm_chi2", &pf_mus_cm_chi2, &b_pf_mus_cm_chi2);
   fChain->SetBranchAddress("pf_mus_cm_ndof", &pf_mus_cm_ndof, &b_pf_mus_cm_ndof);
   fChain->SetBranchAddress("pf_mus_cm_chg", &pf_mus_cm_chg, &b_pf_mus_cm_chg);
   fChain->SetBranchAddress("pf_mus_cm_pt", &pf_mus_cm_pt, &b_pf_mus_cm_pt);
   fChain->SetBranchAddress("pf_mus_cm_px", &pf_mus_cm_px, &b_pf_mus_cm_px);
   fChain->SetBranchAddress("pf_mus_cm_py", &pf_mus_cm_py, &b_pf_mus_cm_py);
   fChain->SetBranchAddress("pf_mus_cm_pz", &pf_mus_cm_pz, &b_pf_mus_cm_pz);
   fChain->SetBranchAddress("pf_mus_cm_eta", &pf_mus_cm_eta, &b_pf_mus_cm_eta);
   fChain->SetBranchAddress("pf_mus_cm_phi", &pf_mus_cm_phi, &b_pf_mus_cm_phi);
   fChain->SetBranchAddress("pf_mus_cm_theta", &pf_mus_cm_theta, &b_pf_mus_cm_theta);
   fChain->SetBranchAddress("pf_mus_cm_d0dum", &pf_mus_cm_d0dum, &b_pf_mus_cm_d0dum);
   fChain->SetBranchAddress("pf_mus_cm_dz", &pf_mus_cm_dz, &b_pf_mus_cm_dz);
   fChain->SetBranchAddress("pf_mus_cm_vx", &pf_mus_cm_vx, &b_pf_mus_cm_vx);
   fChain->SetBranchAddress("pf_mus_cm_vy", &pf_mus_cm_vy, &b_pf_mus_cm_vy);
   fChain->SetBranchAddress("pf_mus_cm_vz", &pf_mus_cm_vz, &b_pf_mus_cm_vz);
   fChain->SetBranchAddress("pf_mus_cm_numvalhits", &pf_mus_cm_numvalhits, &b_pf_mus_cm_numvalhits);
   fChain->SetBranchAddress("pf_mus_cm_numlosthits", &pf_mus_cm_numlosthits, &b_pf_mus_cm_numlosthits);
   fChain->SetBranchAddress("pf_mus_cm_numvalMuonhits", &pf_mus_cm_numvalMuonhits, &b_pf_mus_cm_numvalMuonhits);
   fChain->SetBranchAddress("pf_mus_cm_d0dumErr", &pf_mus_cm_d0dumErr, &b_pf_mus_cm_d0dumErr);
   fChain->SetBranchAddress("pf_mus_cm_dzErr", &pf_mus_cm_dzErr, &b_pf_mus_cm_dzErr);
   fChain->SetBranchAddress("pf_mus_cm_ptErr", &pf_mus_cm_ptErr, &b_pf_mus_cm_ptErr);
   fChain->SetBranchAddress("pf_mus_cm_etaErr", &pf_mus_cm_etaErr, &b_pf_mus_cm_etaErr);
   fChain->SetBranchAddress("pf_mus_cm_phiErr", &pf_mus_cm_phiErr, &b_pf_mus_cm_phiErr);
   fChain->SetBranchAddress("pf_mus_tk_id", &pf_mus_tk_id, &b_pf_mus_tk_id);
   fChain->SetBranchAddress("pf_mus_tk_chi2", &pf_mus_tk_chi2, &b_pf_mus_tk_chi2);
   fChain->SetBranchAddress("pf_mus_tk_ndof", &pf_mus_tk_ndof, &b_pf_mus_tk_ndof);
   fChain->SetBranchAddress("pf_mus_tk_chg", &pf_mus_tk_chg, &b_pf_mus_tk_chg);
   fChain->SetBranchAddress("pf_mus_tk_pt", &pf_mus_tk_pt, &b_pf_mus_tk_pt);
   fChain->SetBranchAddress("pf_mus_tk_px", &pf_mus_tk_px, &b_pf_mus_tk_px);
   fChain->SetBranchAddress("pf_mus_tk_py", &pf_mus_tk_py, &b_pf_mus_tk_py);
   fChain->SetBranchAddress("pf_mus_tk_pz", &pf_mus_tk_pz, &b_pf_mus_tk_pz);
   fChain->SetBranchAddress("pf_mus_tk_eta", &pf_mus_tk_eta, &b_pf_mus_tk_eta);
   fChain->SetBranchAddress("pf_mus_tk_phi", &pf_mus_tk_phi, &b_pf_mus_tk_phi);
   fChain->SetBranchAddress("pf_mus_tk_theta", &pf_mus_tk_theta, &b_pf_mus_tk_theta);
   fChain->SetBranchAddress("pf_mus_tk_d0dum", &pf_mus_tk_d0dum, &b_pf_mus_tk_d0dum);
   fChain->SetBranchAddress("pf_mus_tk_dz", &pf_mus_tk_dz, &b_pf_mus_tk_dz);
   fChain->SetBranchAddress("pf_mus_tk_vx", &pf_mus_tk_vx, &b_pf_mus_tk_vx);
   fChain->SetBranchAddress("pf_mus_tk_vy", &pf_mus_tk_vy, &b_pf_mus_tk_vy);
   fChain->SetBranchAddress("pf_mus_tk_vz", &pf_mus_tk_vz, &b_pf_mus_tk_vz);
   fChain->SetBranchAddress("pf_mus_tk_numvalhits", &pf_mus_tk_numvalhits, &b_pf_mus_tk_numvalhits);
   fChain->SetBranchAddress("pf_mus_tk_numlosthits", &pf_mus_tk_numlosthits, &b_pf_mus_tk_numlosthits);
   fChain->SetBranchAddress("pf_mus_tk_d0dumErr", &pf_mus_tk_d0dumErr, &b_pf_mus_tk_d0dumErr);
   fChain->SetBranchAddress("pf_mus_tk_dzErr", &pf_mus_tk_dzErr, &b_pf_mus_tk_dzErr);
   fChain->SetBranchAddress("pf_mus_tk_ptErr", &pf_mus_tk_ptErr, &b_pf_mus_tk_ptErr);
   fChain->SetBranchAddress("pf_mus_tk_etaErr", &pf_mus_tk_etaErr, &b_pf_mus_tk_etaErr);
   fChain->SetBranchAddress("pf_mus_tk_phiErr", &pf_mus_tk_phiErr, &b_pf_mus_tk_phiErr);
   fChain->SetBranchAddress("pf_mus_tk_numvalPixelhits", &pf_mus_tk_numvalPixelhits, &b_pf_mus_tk_numvalPixelhits);
   fChain->SetBranchAddress("pf_mus_tk_numpixelWthMeasr", &pf_mus_tk_numpixelWthMeasr, &b_pf_mus_tk_numpixelWthMeasr);
   fChain->SetBranchAddress("pf_mus_stamu_chi2", &pf_mus_stamu_chi2, &b_pf_mus_stamu_chi2);
   fChain->SetBranchAddress("pf_mus_stamu_ndof", &pf_mus_stamu_ndof, &b_pf_mus_stamu_ndof);
   fChain->SetBranchAddress("pf_mus_stamu_chg", &pf_mus_stamu_chg, &b_pf_mus_stamu_chg);
   fChain->SetBranchAddress("pf_mus_stamu_pt", &pf_mus_stamu_pt, &b_pf_mus_stamu_pt);
   fChain->SetBranchAddress("pf_mus_stamu_px", &pf_mus_stamu_px, &b_pf_mus_stamu_px);
   fChain->SetBranchAddress("pf_mus_stamu_py", &pf_mus_stamu_py, &b_pf_mus_stamu_py);
   fChain->SetBranchAddress("pf_mus_stamu_pz", &pf_mus_stamu_pz, &b_pf_mus_stamu_pz);
   fChain->SetBranchAddress("pf_mus_stamu_eta", &pf_mus_stamu_eta, &b_pf_mus_stamu_eta);
   fChain->SetBranchAddress("pf_mus_stamu_phi", &pf_mus_stamu_phi, &b_pf_mus_stamu_phi);
   fChain->SetBranchAddress("pf_mus_stamu_theta", &pf_mus_stamu_theta, &b_pf_mus_stamu_theta);
   fChain->SetBranchAddress("pf_mus_stamu_d0dum", &pf_mus_stamu_d0dum, &b_pf_mus_stamu_d0dum);
   fChain->SetBranchAddress("pf_mus_stamu_dz", &pf_mus_stamu_dz, &b_pf_mus_stamu_dz);
   fChain->SetBranchAddress("pf_mus_stamu_vx", &pf_mus_stamu_vx, &b_pf_mus_stamu_vx);
   fChain->SetBranchAddress("pf_mus_stamu_vy", &pf_mus_stamu_vy, &b_pf_mus_stamu_vy);
   fChain->SetBranchAddress("pf_mus_stamu_vz", &pf_mus_stamu_vz, &b_pf_mus_stamu_vz);
   fChain->SetBranchAddress("pf_mus_stamu_numvalhits", &pf_mus_stamu_numvalhits, &b_pf_mus_stamu_numvalhits);
   fChain->SetBranchAddress("pf_mus_stamu_numlosthits", &pf_mus_stamu_numlosthits, &b_pf_mus_stamu_numlosthits);
   fChain->SetBranchAddress("pf_mus_stamu_d0dumErr", &pf_mus_stamu_d0dumErr, &b_pf_mus_stamu_d0dumErr);
   fChain->SetBranchAddress("pf_mus_stamu_dzErr", &pf_mus_stamu_dzErr, &b_pf_mus_stamu_dzErr);
   fChain->SetBranchAddress("pf_mus_stamu_ptErr", &pf_mus_stamu_ptErr, &b_pf_mus_stamu_ptErr);
   fChain->SetBranchAddress("pf_mus_stamu_etaErr", &pf_mus_stamu_etaErr, &b_pf_mus_stamu_etaErr);
   fChain->SetBranchAddress("pf_mus_stamu_phiErr", &pf_mus_stamu_phiErr, &b_pf_mus_stamu_phiErr);
   fChain->SetBranchAddress("pf_mus_num_matches", &pf_mus_num_matches, &b_pf_mus_num_matches);
   fChain->SetBranchAddress("pf_mus_isTrackerMuon", &pf_mus_isTrackerMuon, &b_pf_mus_isTrackerMuon);
   fChain->SetBranchAddress("pf_mus_isStandAloneMuon", &pf_mus_isStandAloneMuon, &b_pf_mus_isStandAloneMuon);
   fChain->SetBranchAddress("pf_mus_isCaloMuon", &pf_mus_isCaloMuon, &b_pf_mus_isCaloMuon);
   fChain->SetBranchAddress("pf_mus_isGlobalMuon", &pf_mus_isGlobalMuon, &b_pf_mus_isGlobalMuon);
   fChain->SetBranchAddress("pf_mus_isElectron", &pf_mus_isElectron, &b_pf_mus_isElectron);
   fChain->SetBranchAddress("pf_mus_isConvertedPhoton", &pf_mus_isConvertedPhoton, &b_pf_mus_isConvertedPhoton);
   fChain->SetBranchAddress("pf_mus_isPhoton", &pf_mus_isPhoton, &b_pf_mus_isPhoton);
   fChain->SetBranchAddress("pf_mus_id_All", &pf_mus_id_All, &b_pf_mus_id_All);
   fChain->SetBranchAddress("pf_mus_id_AllGlobalMuons", &pf_mus_id_AllGlobalMuons, &b_pf_mus_id_AllGlobalMuons);
   fChain->SetBranchAddress("pf_mus_id_AllStandAloneMuons", &pf_mus_id_AllStandAloneMuons, &b_pf_mus_id_AllStandAloneMuons);
   fChain->SetBranchAddress("pf_mus_id_AllTrackerMuons", &pf_mus_id_AllTrackerMuons, &b_pf_mus_id_AllTrackerMuons);
   fChain->SetBranchAddress("pf_mus_id_TrackerMuonArbitrated", &pf_mus_id_TrackerMuonArbitrated, &b_pf_mus_id_TrackerMuonArbitrated);
   fChain->SetBranchAddress("pf_mus_id_AllArbitrated", &pf_mus_id_AllArbitrated, &b_pf_mus_id_AllArbitrated);
   fChain->SetBranchAddress("pf_mus_id_GlobalMuonPromptTight", &pf_mus_id_GlobalMuonPromptTight, &b_pf_mus_id_GlobalMuonPromptTight);
   fChain->SetBranchAddress("pf_mus_id_TMLastStationLoose", &pf_mus_id_TMLastStationLoose, &b_pf_mus_id_TMLastStationLoose);
   fChain->SetBranchAddress("pf_mus_id_TMLastStationTight", &pf_mus_id_TMLastStationTight, &b_pf_mus_id_TMLastStationTight);
   fChain->SetBranchAddress("pf_mus_id_TM2DCompatibilityLoose", &pf_mus_id_TM2DCompatibilityLoose, &b_pf_mus_id_TM2DCompatibilityLoose);
   fChain->SetBranchAddress("pf_mus_id_TM2DCompatibilityTight", &pf_mus_id_TM2DCompatibilityTight, &b_pf_mus_id_TM2DCompatibilityTight);
   fChain->SetBranchAddress("pf_mus_id_TMOneStationLoose", &pf_mus_id_TMOneStationLoose, &b_pf_mus_id_TMOneStationLoose);
   fChain->SetBranchAddress("pf_mus_id_TMOneStationTight", &pf_mus_id_TMOneStationTight, &b_pf_mus_id_TMOneStationTight);
   fChain->SetBranchAddress("pf_mus_id_TMLastStationOptimizedLowPtLoose", &pf_mus_id_TMLastStationOptimizedLowPtLoose, &b_pf_mus_id_TMLastStationOptimizedLowPtLoose);
   fChain->SetBranchAddress("pf_mus_id_TMLastStationOptimizedLowPtTight", &pf_mus_id_TMLastStationOptimizedLowPtTight, &b_pf_mus_id_TMLastStationOptimizedLowPtTight);
   fChain->SetBranchAddress("pf_mus_tk_LayersWithMeasurement", &pf_mus_tk_LayersWithMeasurement, &b_pf_mus_tk_LayersWithMeasurement);
   fChain->SetBranchAddress("pf_mus_tk_PixelLayersWithMeasurement", &pf_mus_tk_PixelLayersWithMeasurement, &b_pf_mus_tk_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("pf_mus_tk_ValidStripLayersWithMonoAndStereoHit", &pf_mus_tk_ValidStripLayersWithMonoAndStereoHit, &b_pf_mus_tk_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("pf_mus_tk_LayersWithoutMeasurement", &pf_mus_tk_LayersWithoutMeasurement, &b_pf_mus_tk_LayersWithoutMeasurement);
   fChain->SetBranchAddress("pf_mus_tk_ExpectedHitsInner", &pf_mus_tk_ExpectedHitsInner, &b_pf_mus_tk_ExpectedHitsInner);
   fChain->SetBranchAddress("pf_mus_tk_ExpectedHitsOuter", &pf_mus_tk_ExpectedHitsOuter, &b_pf_mus_tk_ExpectedHitsOuter);
   fChain->SetBranchAddress("pf_mus_cm_LayersWithMeasurement", &pf_mus_cm_LayersWithMeasurement, &b_pf_mus_cm_LayersWithMeasurement);
   fChain->SetBranchAddress("pf_mus_cm_PixelLayersWithMeasurement", &pf_mus_cm_PixelLayersWithMeasurement, &b_pf_mus_cm_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("pf_mus_cm_ValidStripLayersWithMonoAndStereoHit", &pf_mus_cm_ValidStripLayersWithMonoAndStereoHit, &b_pf_mus_cm_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("pf_mus_cm_LayersWithoutMeasurement", &pf_mus_cm_LayersWithoutMeasurement, &b_pf_mus_cm_LayersWithoutMeasurement);
   fChain->SetBranchAddress("pf_mus_cm_ExpectedHitsInner", &pf_mus_cm_ExpectedHitsInner, &b_pf_mus_cm_ExpectedHitsInner);
   fChain->SetBranchAddress("pf_mus_cm_ExpectedHitsOuter", &pf_mus_cm_ExpectedHitsOuter, &b_pf_mus_cm_ExpectedHitsOuter);
   fChain->SetBranchAddress("pf_mus_picky_LayersWithMeasurement", &pf_mus_picky_LayersWithMeasurement, &b_pf_mus_picky_LayersWithMeasurement);
   fChain->SetBranchAddress("pf_mus_picky_PixelLayersWithMeasurement", &pf_mus_picky_PixelLayersWithMeasurement, &b_pf_mus_picky_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("pf_mus_picky_ValidStripLayersWithMonoAndStereoHit", &pf_mus_picky_ValidStripLayersWithMonoAndStereoHit, &b_pf_mus_picky_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("pf_mus_picky_LayersWithoutMeasurement", &pf_mus_picky_LayersWithoutMeasurement, &b_pf_mus_picky_LayersWithoutMeasurement);
   fChain->SetBranchAddress("pf_mus_picky_ExpectedHitsInner", &pf_mus_picky_ExpectedHitsInner, &b_pf_mus_picky_ExpectedHitsInner);
   fChain->SetBranchAddress("pf_mus_picky_ExpectedHitsOuter", &pf_mus_picky_ExpectedHitsOuter, &b_pf_mus_picky_ExpectedHitsOuter);
   fChain->SetBranchAddress("pf_mus_tpfms_LayersWithMeasurement", &pf_mus_tpfms_LayersWithMeasurement, &b_pf_mus_tpfms_LayersWithMeasurement);
   fChain->SetBranchAddress("pf_mus_tpfms_PixelLayersWithMeasurement", &pf_mus_tpfms_PixelLayersWithMeasurement, &b_pf_mus_tpfms_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("pf_mus_tpfms_ValidStripLayersWithMonoAndStereoHit", &pf_mus_tpfms_ValidStripLayersWithMonoAndStereoHit, &b_pf_mus_tpfms_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("pf_mus_tpfms_LayersWithoutMeasurement", &pf_mus_tpfms_LayersWithoutMeasurement, &b_pf_mus_tpfms_LayersWithoutMeasurement);
   fChain->SetBranchAddress("pf_mus_tpfms_ExpectedHitsInner", &pf_mus_tpfms_ExpectedHitsInner, &b_pf_mus_tpfms_ExpectedHitsInner);
   fChain->SetBranchAddress("pf_mus_tpfms_ExpectedHitsOuter", &pf_mus_tpfms_ExpectedHitsOuter, &b_pf_mus_tpfms_ExpectedHitsOuter);
   fChain->SetBranchAddress("pf_mus_picky_id", &pf_mus_picky_id, &b_pf_mus_picky_id);
   fChain->SetBranchAddress("pf_mus_picky_chi2", &pf_mus_picky_chi2, &b_pf_mus_picky_chi2);
   fChain->SetBranchAddress("pf_mus_picky_ndof", &pf_mus_picky_ndof, &b_pf_mus_picky_ndof);
   fChain->SetBranchAddress("pf_mus_picky_chg", &pf_mus_picky_chg, &b_pf_mus_picky_chg);
   fChain->SetBranchAddress("pf_mus_picky_pt", &pf_mus_picky_pt, &b_pf_mus_picky_pt);
   fChain->SetBranchAddress("pf_mus_picky_px", &pf_mus_picky_px, &b_pf_mus_picky_px);
   fChain->SetBranchAddress("pf_mus_picky_py", &pf_mus_picky_py, &b_pf_mus_picky_py);
   fChain->SetBranchAddress("pf_mus_picky_pz", &pf_mus_picky_pz, &b_pf_mus_picky_pz);
   fChain->SetBranchAddress("pf_mus_picky_eta", &pf_mus_picky_eta, &b_pf_mus_picky_eta);
   fChain->SetBranchAddress("pf_mus_picky_phi", &pf_mus_picky_phi, &b_pf_mus_picky_phi);
   fChain->SetBranchAddress("pf_mus_picky_theta", &pf_mus_picky_theta, &b_pf_mus_picky_theta);
   fChain->SetBranchAddress("pf_mus_picky_d0dum", &pf_mus_picky_d0dum, &b_pf_mus_picky_d0dum);
   fChain->SetBranchAddress("pf_mus_picky_dz", &pf_mus_picky_dz, &b_pf_mus_picky_dz);
   fChain->SetBranchAddress("pf_mus_picky_vx", &pf_mus_picky_vx, &b_pf_mus_picky_vx);
   fChain->SetBranchAddress("pf_mus_picky_vy", &pf_mus_picky_vy, &b_pf_mus_picky_vy);
   fChain->SetBranchAddress("pf_mus_picky_vz", &pf_mus_picky_vz, &b_pf_mus_picky_vz);
   fChain->SetBranchAddress("pf_mus_picky_numvalhits", &pf_mus_picky_numvalhits, &b_pf_mus_picky_numvalhits);
   fChain->SetBranchAddress("pf_mus_picky_numlosthits", &pf_mus_picky_numlosthits, &b_pf_mus_picky_numlosthits);
   fChain->SetBranchAddress("pf_mus_picky_d0dumErr", &pf_mus_picky_d0dumErr, &b_pf_mus_picky_d0dumErr);
   fChain->SetBranchAddress("pf_mus_picky_dzErr", &pf_mus_picky_dzErr, &b_pf_mus_picky_dzErr);
   fChain->SetBranchAddress("pf_mus_picky_ptErr", &pf_mus_picky_ptErr, &b_pf_mus_picky_ptErr);
   fChain->SetBranchAddress("pf_mus_picky_etaErr", &pf_mus_picky_etaErr, &b_pf_mus_picky_etaErr);
   fChain->SetBranchAddress("pf_mus_picky_phiErr", &pf_mus_picky_phiErr, &b_pf_mus_picky_phiErr);
   fChain->SetBranchAddress("pf_mus_picky_numvalPixelhits", &pf_mus_picky_numvalPixelhits, &b_pf_mus_picky_numvalPixelhits);
   fChain->SetBranchAddress("pf_mus_tpfms_id", &pf_mus_tpfms_id, &b_pf_mus_tpfms_id);
   fChain->SetBranchAddress("pf_mus_tpfms_chi2", &pf_mus_tpfms_chi2, &b_pf_mus_tpfms_chi2);
   fChain->SetBranchAddress("pf_mus_tpfms_ndof", &pf_mus_tpfms_ndof, &b_pf_mus_tpfms_ndof);
   fChain->SetBranchAddress("pf_mus_tpfms_chg", &pf_mus_tpfms_chg, &b_pf_mus_tpfms_chg);
   fChain->SetBranchAddress("pf_mus_tpfms_pt", &pf_mus_tpfms_pt, &b_pf_mus_tpfms_pt);
   fChain->SetBranchAddress("pf_mus_tpfms_px", &pf_mus_tpfms_px, &b_pf_mus_tpfms_px);
   fChain->SetBranchAddress("pf_mus_tpfms_py", &pf_mus_tpfms_py, &b_pf_mus_tpfms_py);
   fChain->SetBranchAddress("pf_mus_tpfms_pz", &pf_mus_tpfms_pz, &b_pf_mus_tpfms_pz);
   fChain->SetBranchAddress("pf_mus_tpfms_eta", &pf_mus_tpfms_eta, &b_pf_mus_tpfms_eta);
   fChain->SetBranchAddress("pf_mus_tpfms_phi", &pf_mus_tpfms_phi, &b_pf_mus_tpfms_phi);
   fChain->SetBranchAddress("pf_mus_tpfms_theta", &pf_mus_tpfms_theta, &b_pf_mus_tpfms_theta);
   fChain->SetBranchAddress("pf_mus_tpfms_d0dum", &pf_mus_tpfms_d0dum, &b_pf_mus_tpfms_d0dum);
   fChain->SetBranchAddress("pf_mus_tpfms_dz", &pf_mus_tpfms_dz, &b_pf_mus_tpfms_dz);
   fChain->SetBranchAddress("pf_mus_tpfms_vx", &pf_mus_tpfms_vx, &b_pf_mus_tpfms_vx);
   fChain->SetBranchAddress("pf_mus_tpfms_vy", &pf_mus_tpfms_vy, &b_pf_mus_tpfms_vy);
   fChain->SetBranchAddress("pf_mus_tpfms_vz", &pf_mus_tpfms_vz, &b_pf_mus_tpfms_vz);
   fChain->SetBranchAddress("pf_mus_tpfms_numvalhits", &pf_mus_tpfms_numvalhits, &b_pf_mus_tpfms_numvalhits);
   fChain->SetBranchAddress("pf_mus_tpfms_numlosthits", &pf_mus_tpfms_numlosthits, &b_pf_mus_tpfms_numlosthits);
   fChain->SetBranchAddress("pf_mus_tpfms_d0dumErr", &pf_mus_tpfms_d0dumErr, &b_pf_mus_tpfms_d0dumErr);
   fChain->SetBranchAddress("pf_mus_tpfms_dzErr", &pf_mus_tpfms_dzErr, &b_pf_mus_tpfms_dzErr);
   fChain->SetBranchAddress("pf_mus_tpfms_ptErr", &pf_mus_tpfms_ptErr, &b_pf_mus_tpfms_ptErr);
   fChain->SetBranchAddress("pf_mus_tpfms_etaErr", &pf_mus_tpfms_etaErr, &b_pf_mus_tpfms_etaErr);
   fChain->SetBranchAddress("pf_mus_tpfms_phiErr", &pf_mus_tpfms_phiErr, &b_pf_mus_tpfms_phiErr);
   fChain->SetBranchAddress("pf_mus_tpfms_numvalPixelhits", &pf_mus_tpfms_numvalPixelhits, &b_pf_mus_tpfms_numvalPixelhits);
   fChain->SetBranchAddress("pf_mus_pfIsolationR03_sumChargedHadronPt", &pf_mus_pfIsolationR03_sumChargedHadronPt, &b_pf_mus_pfIsolationR03_sumChargedHadronPt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR03_sumChargedParticlePt", &pf_mus_pfIsolationR03_sumChargedParticlePt, &b_pf_mus_pfIsolationR03_sumChargedParticlePt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR03_sumNeutralHadronEt", &pf_mus_pfIsolationR03_sumNeutralHadronEt, &b_pf_mus_pfIsolationR03_sumNeutralHadronEt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR03_sumNeutralHadronEtHighThreshold", &pf_mus_pfIsolationR03_sumNeutralHadronEtHighThreshold, &b_pf_mus_pfIsolationR03_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("pf_mus_pfIsolationR03_sumPhotonEt", &pf_mus_pfIsolationR03_sumPhotonEt, &b_pf_mus_pfIsolationR03_sumPhotonEt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR03_sumPhotonEtHighThreshold", &pf_mus_pfIsolationR03_sumPhotonEtHighThreshold, &b_pf_mus_pfIsolationR03_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("pf_mus_pfIsolationR03_sumPUPt", &pf_mus_pfIsolationR03_sumPUPt, &b_pf_mus_pfIsolationR03_sumPUPt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR04_sumChargedHadronPt", &pf_mus_pfIsolationR04_sumChargedHadronPt, &b_pf_mus_pfIsolationR04_sumChargedHadronPt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR04_sumChargedParticlePt", &pf_mus_pfIsolationR04_sumChargedParticlePt, &b_pf_mus_pfIsolationR04_sumChargedParticlePt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR04_sumNeutralHadronEt", &pf_mus_pfIsolationR04_sumNeutralHadronEt, &b_pf_mus_pfIsolationR04_sumNeutralHadronEt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR04_sumNeutralHadronEtHighThreshold", &pf_mus_pfIsolationR04_sumNeutralHadronEtHighThreshold, &b_pf_mus_pfIsolationR04_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("pf_mus_pfIsolationR04_sumPhotonEt", &pf_mus_pfIsolationR04_sumPhotonEt, &b_pf_mus_pfIsolationR04_sumPhotonEt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR04_sumPhotonEtHighThreshold", &pf_mus_pfIsolationR04_sumPhotonEtHighThreshold, &b_pf_mus_pfIsolationR04_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("pf_mus_pfIsolationR04_sumPUPt", &pf_mus_pfIsolationR04_sumPUPt, &b_pf_mus_pfIsolationR04_sumPUPt);
   fChain->SetBranchAddress("pf_mus_dB", &pf_mus_dB, &b_pf_mus_dB);
   fChain->SetBranchAddress("pf_mus_numberOfMatchedStations", &pf_mus_numberOfMatchedStations, &b_pf_mus_numberOfMatchedStations);
   fChain->SetBranchAddress("pf_mus_isPFMuon", &pf_mus_isPFMuon, &b_pf_mus_isPFMuon);
   fChain->SetBranchAddress("Npfcand", &Npfcand, &b_Npfcand);
   fChain->SetBranchAddress("pfcand_pdgId", &pfcand_pdgId, &b_pfcand_pdgId);
   fChain->SetBranchAddress("pfcand_particleId", &pfcand_particleId, &b_pfcand_particleId);
   fChain->SetBranchAddress("pfcand_pt", &pfcand_pt, &b_pfcand_pt);
   fChain->SetBranchAddress("pfcand_pz", &pfcand_pz, &b_pfcand_pz);
   fChain->SetBranchAddress("pfcand_px", &pfcand_px, &b_pfcand_px);
   fChain->SetBranchAddress("pfcand_py", &pfcand_py, &b_pfcand_py);
   fChain->SetBranchAddress("pfcand_eta", &pfcand_eta, &b_pfcand_eta);
   fChain->SetBranchAddress("pfcand_phi", &pfcand_phi, &b_pfcand_phi);
   fChain->SetBranchAddress("pfcand_theta", &pfcand_theta, &b_pfcand_theta);
   fChain->SetBranchAddress("pfcand_energy", &pfcand_energy, &b_pfcand_energy);
   fChain->SetBranchAddress("pfcand_charge", &pfcand_charge, &b_pfcand_charge);
   fChain->SetBranchAddress("Npfmets", &Npfmets, &b_Npfmets);
   fChain->SetBranchAddress("pfmets_et", &pfmets_et, &b_pfmets_et);
   fChain->SetBranchAddress("pfmets_phi", &pfmets_phi, &b_pfmets_phi);
   fChain->SetBranchAddress("pfmets_ex", &pfmets_ex, &b_pfmets_ex);
   fChain->SetBranchAddress("pfmets_ey", &pfmets_ey, &b_pfmets_ey);
   fChain->SetBranchAddress("pfmets_gen_et", &pfmets_gen_et, &b_pfmets_gen_et);
   fChain->SetBranchAddress("pfmets_gen_phi", &pfmets_gen_phi, &b_pfmets_gen_phi);
   fChain->SetBranchAddress("pfmets_sign", &pfmets_sign, &b_pfmets_sign);
   fChain->SetBranchAddress("pfmets_sumEt", &pfmets_sumEt, &b_pfmets_sumEt);
   fChain->SetBranchAddress("pfmets_unCPhi", &pfmets_unCPhi, &b_pfmets_unCPhi);
   fChain->SetBranchAddress("pfmets_unCPt", &pfmets_unCPt, &b_pfmets_unCPt);
   fChain->SetBranchAddress("Nphotons", &Nphotons, &b_Nphotons);
   fChain->SetBranchAddress("photons_energy", &photons_energy, &b_photons_energy);
   fChain->SetBranchAddress("photons_et", &photons_et, &b_photons_et);
   fChain->SetBranchAddress("photons_eta", &photons_eta, &b_photons_eta);
   fChain->SetBranchAddress("photons_phi", &photons_phi, &b_photons_phi);
   fChain->SetBranchAddress("photons_pt", &photons_pt, &b_photons_pt);
   fChain->SetBranchAddress("photons_px", &photons_px, &b_photons_px);
   fChain->SetBranchAddress("photons_py", &photons_py, &b_photons_py);
   fChain->SetBranchAddress("photons_pz", &photons_pz, &b_photons_pz);
   fChain->SetBranchAddress("photons_status", &photons_status, &b_photons_status);
   fChain->SetBranchAddress("photons_theta", &photons_theta, &b_photons_theta);
   fChain->SetBranchAddress("photons_hadOverEM", &photons_hadOverEM, &b_photons_hadOverEM);
   fChain->SetBranchAddress("photons_scEnergy", &photons_scEnergy, &b_photons_scEnergy);
   fChain->SetBranchAddress("photons_scRawEnergy", &photons_scRawEnergy, &b_photons_scRawEnergy);
   fChain->SetBranchAddress("photons_scEta", &photons_scEta, &b_photons_scEta);
   fChain->SetBranchAddress("photons_scPhi", &photons_scPhi, &b_photons_scPhi);
   fChain->SetBranchAddress("photons_scEtaWidth", &photons_scEtaWidth, &b_photons_scEtaWidth);
   fChain->SetBranchAddress("photons_scPhiWidth", &photons_scPhiWidth, &b_photons_scPhiWidth);
   fChain->SetBranchAddress("photons_tIso", &photons_tIso, &b_photons_tIso);
   fChain->SetBranchAddress("photons_ecalIso", &photons_ecalIso, &b_photons_ecalIso);
   fChain->SetBranchAddress("photons_hcalIso", &photons_hcalIso, &b_photons_hcalIso);
   fChain->SetBranchAddress("photons_isoEcalRecHitDR04", &photons_isoEcalRecHitDR04, &b_photons_isoEcalRecHitDR04);
   fChain->SetBranchAddress("photons_isoHcalRecHitDR04", &photons_isoHcalRecHitDR04, &b_photons_isoHcalRecHitDR04);
   fChain->SetBranchAddress("photons_isoSolidTrkConeDR04", &photons_isoSolidTrkConeDR04, &b_photons_isoSolidTrkConeDR04);
   fChain->SetBranchAddress("photons_isoHollowTrkConeDR04", &photons_isoHollowTrkConeDR04, &b_photons_isoHollowTrkConeDR04);
   fChain->SetBranchAddress("photons_nTrkSolidConeDR04", &photons_nTrkSolidConeDR04, &b_photons_nTrkSolidConeDR04);
   fChain->SetBranchAddress("photons_nTrkHollowConeDR04", &photons_nTrkHollowConeDR04, &b_photons_nTrkHollowConeDR04);
   fChain->SetBranchAddress("photons_isoEcalRecHitDR03", &photons_isoEcalRecHitDR03, &b_photons_isoEcalRecHitDR03);
   fChain->SetBranchAddress("photons_isoHcalRecHitDR03", &photons_isoHcalRecHitDR03, &b_photons_isoHcalRecHitDR03);
   fChain->SetBranchAddress("photons_isoSolidTrkConeDR03", &photons_isoSolidTrkConeDR03, &b_photons_isoSolidTrkConeDR03);
   fChain->SetBranchAddress("photons_isoHollowTrkConeDR03", &photons_isoHollowTrkConeDR03, &b_photons_isoHollowTrkConeDR03);
   fChain->SetBranchAddress("photons_nTrkSolidConeDR03", &photons_nTrkSolidConeDR03, &b_photons_nTrkSolidConeDR03);
   fChain->SetBranchAddress("photons_nTrkHollowConeDR03", &photons_nTrkHollowConeDR03, &b_photons_nTrkHollowConeDR03);
   fChain->SetBranchAddress("photons_isAlsoElectron", &photons_isAlsoElectron, &b_photons_isAlsoElectron);
   fChain->SetBranchAddress("photons_hasPixelSeed", &photons_hasPixelSeed, &b_photons_hasPixelSeed);
   fChain->SetBranchAddress("photons_isConverted", &photons_isConverted, &b_photons_isConverted);
   fChain->SetBranchAddress("photons_isEBGap", &photons_isEBGap, &b_photons_isEBGap);
   fChain->SetBranchAddress("photons_isEEGap", &photons_isEEGap, &b_photons_isEEGap);
   fChain->SetBranchAddress("photons_isEBEEGap", &photons_isEBEEGap, &b_photons_isEBEEGap);
   fChain->SetBranchAddress("photons_isEBPho", &photons_isEBPho, &b_photons_isEBPho);
   fChain->SetBranchAddress("photons_isEEPho", &photons_isEEPho, &b_photons_isEEPho);
   fChain->SetBranchAddress("photons_isLoosePhoton", &photons_isLoosePhoton, &b_photons_isLoosePhoton);
   fChain->SetBranchAddress("photons_isTightPhoton", &photons_isTightPhoton, &b_photons_isTightPhoton);
   fChain->SetBranchAddress("photons_maxEnergyXtal", &photons_maxEnergyXtal, &b_photons_maxEnergyXtal);
   fChain->SetBranchAddress("photons_e1x5", &photons_e1x5, &b_photons_e1x5);
   fChain->SetBranchAddress("photons_e2x5", &photons_e2x5, &b_photons_e2x5);
   fChain->SetBranchAddress("photons_e3x3", &photons_e3x3, &b_photons_e3x3);
   fChain->SetBranchAddress("photons_e5x5", &photons_e5x5, &b_photons_e5x5);
   fChain->SetBranchAddress("photons_sigmaEtaEta", &photons_sigmaEtaEta, &b_photons_sigmaEtaEta);
   fChain->SetBranchAddress("photons_sigmaIetaIeta", &photons_sigmaIetaIeta, &b_photons_sigmaIetaIeta);
   fChain->SetBranchAddress("photons_r9", &photons_r9, &b_photons_r9);
   fChain->SetBranchAddress("photons_gen_et", &photons_gen_et, &b_photons_gen_et);
   fChain->SetBranchAddress("photons_gen_eta", &photons_gen_eta, &b_photons_gen_eta);
   fChain->SetBranchAddress("photons_gen_phi", &photons_gen_phi, &b_photons_gen_phi);
   fChain->SetBranchAddress("photons_gen_id", &photons_gen_id, &b_photons_gen_id);
   fChain->SetBranchAddress("Npv", &Npv, &b_Npv);
   fChain->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
   fChain->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
   fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   fChain->SetBranchAddress("pv_xErr", &pv_xErr, &b_pv_xErr);
   fChain->SetBranchAddress("pv_yErr", &pv_yErr, &b_pv_yErr);
   fChain->SetBranchAddress("pv_zErr", &pv_zErr, &b_pv_zErr);
   fChain->SetBranchAddress("pv_chi2", &pv_chi2, &b_pv_chi2);
   fChain->SetBranchAddress("pv_ndof", &pv_ndof, &b_pv_ndof);
   fChain->SetBranchAddress("pv_isFake", &pv_isFake, &b_pv_isFake);
   fChain->SetBranchAddress("pv_isValid", &pv_isValid, &b_pv_isValid);
   fChain->SetBranchAddress("pv_tracksSize", &pv_tracksSize, &b_pv_tracksSize);
   fChain->SetBranchAddress("Ntaus", &Ntaus, &b_Ntaus);
   fChain->SetBranchAddress("taus_status", &taus_status, &b_taus_status);
   fChain->SetBranchAddress("taus_phi", &taus_phi, &b_taus_phi);
   fChain->SetBranchAddress("taus_pt", &taus_pt, &b_taus_pt);
   fChain->SetBranchAddress("taus_pz", &taus_pz, &b_taus_pz);
   fChain->SetBranchAddress("taus_px", &taus_px, &b_taus_px);
   fChain->SetBranchAddress("taus_py", &taus_py, &b_taus_py);
   fChain->SetBranchAddress("taus_eta", &taus_eta, &b_taus_eta);
   fChain->SetBranchAddress("taus_theta", &taus_theta, &b_taus_theta);
   fChain->SetBranchAddress("taus_et", &taus_et, &b_taus_et);
   fChain->SetBranchAddress("taus_energy", &taus_energy, &b_taus_energy);
   fChain->SetBranchAddress("taus_charge", &taus_charge, &b_taus_charge);
   fChain->SetBranchAddress("taus_emf", &taus_emf, &b_taus_emf);
   fChain->SetBranchAddress("taus_hcalTotOverPLead", &taus_hcalTotOverPLead, &b_taus_hcalTotOverPLead);
   fChain->SetBranchAddress("taus_hcalMaxOverPLead", &taus_hcalMaxOverPLead, &b_taus_hcalMaxOverPLead);
   fChain->SetBranchAddress("taus_hcal3x3OverPLead", &taus_hcal3x3OverPLead, &b_taus_hcal3x3OverPLead);
   fChain->SetBranchAddress("taus_ecalStripSumEOverPLead", &taus_ecalStripSumEOverPLead, &b_taus_ecalStripSumEOverPLead);
   fChain->SetBranchAddress("taus_elecPreIdOutput", &taus_elecPreIdOutput, &b_taus_elecPreIdOutput);
   fChain->SetBranchAddress("taus_elecPreIdDecision", &taus_elecPreIdDecision, &b_taus_elecPreIdDecision);
   fChain->SetBranchAddress("taus_leadPFChargedHadrCand_pt", &taus_leadPFChargedHadrCand_pt, &b_taus_leadPFChargedHadrCand_pt);
   fChain->SetBranchAddress("taus_leadPFChargedHadrCand_charge", &taus_leadPFChargedHadrCand_charge, &b_taus_leadPFChargedHadrCand_charge);
   fChain->SetBranchAddress("taus_leadPFChargedHadrCand_eta", &taus_leadPFChargedHadrCand_eta, &b_taus_leadPFChargedHadrCand_eta);
   fChain->SetBranchAddress("taus_leadPFChargedHadrCand_ECAL_eta", &taus_leadPFChargedHadrCand_ECAL_eta, &b_taus_leadPFChargedHadrCand_ECAL_eta);
   fChain->SetBranchAddress("taus_leadPFChargedHadrCand_phi", &taus_leadPFChargedHadrCand_phi, &b_taus_leadPFChargedHadrCand_phi);
   fChain->SetBranchAddress("taus_isoPFGammaCandsEtSum", &taus_isoPFGammaCandsEtSum, &b_taus_isoPFGammaCandsEtSum);
   fChain->SetBranchAddress("taus_isoPFChargedHadrCandsPtSum", &taus_isoPFChargedHadrCandsPtSum, &b_taus_isoPFChargedHadrCandsPtSum);
   fChain->SetBranchAddress("taus_leadingTrackFinding", &taus_leadingTrackFinding, &b_taus_leadingTrackFinding);
   fChain->SetBranchAddress("taus_leadingTrackPtCut", &taus_leadingTrackPtCut, &b_taus_leadingTrackPtCut);
   fChain->SetBranchAddress("taus_trackIsolation", &taus_trackIsolation, &b_taus_trackIsolation);
   fChain->SetBranchAddress("taus_ecalIsolation", &taus_ecalIsolation, &b_taus_ecalIsolation);
   fChain->SetBranchAddress("taus_byIsolation", &taus_byIsolation, &b_taus_byIsolation);
   fChain->SetBranchAddress("taus_againstElectron", &taus_againstElectron, &b_taus_againstElectron);
   fChain->SetBranchAddress("taus_againstMuon", &taus_againstMuon, &b_taus_againstMuon);
   fChain->SetBranchAddress("taus_taNC_quarter", &taus_taNC_quarter, &b_taus_taNC_quarter);
   fChain->SetBranchAddress("taus_taNC_one", &taus_taNC_one, &b_taus_taNC_one);
   fChain->SetBranchAddress("taus_taNC_half", &taus_taNC_half, &b_taus_taNC_half);
   fChain->SetBranchAddress("taus_taNC_tenth", &taus_taNC_tenth, &b_taus_taNC_tenth);
   fChain->SetBranchAddress("taus_taNC", &taus_taNC, &b_taus_taNC);
   fChain->SetBranchAddress("taus_byIsoUsingLeadingPi", &taus_byIsoUsingLeadingPi, &b_taus_byIsoUsingLeadingPi);
   fChain->SetBranchAddress("taus_tkIsoUsingLeadingPi", &taus_tkIsoUsingLeadingPi, &b_taus_tkIsoUsingLeadingPi);
   fChain->SetBranchAddress("taus_ecalIsoUsingLeadingPi", &taus_ecalIsoUsingLeadingPi, &b_taus_ecalIsoUsingLeadingPi);
   fChain->SetBranchAddress("taus_againstElectronLoose", &taus_againstElectronLoose, &b_taus_againstElectronLoose);
   fChain->SetBranchAddress("taus_againstElectronMedium", &taus_againstElectronMedium, &b_taus_againstElectronMedium);
   fChain->SetBranchAddress("taus_againstElectronTight", &taus_againstElectronTight, &b_taus_againstElectronTight);
   fChain->SetBranchAddress("taus_againstElectronMVA", &taus_againstElectronMVA, &b_taus_againstElectronMVA);
   fChain->SetBranchAddress("taus_againstMuonLoose", &taus_againstMuonLoose, &b_taus_againstMuonLoose);
   fChain->SetBranchAddress("taus_againstMuonMedium", &taus_againstMuonMedium, &b_taus_againstMuonMedium);
   fChain->SetBranchAddress("taus_againstMuonTight", &taus_againstMuonTight, &b_taus_againstMuonTight);
   fChain->SetBranchAddress("taus_decayModeFinding", &taus_decayModeFinding, &b_taus_decayModeFinding);
   fChain->SetBranchAddress("taus_byVLooseIsolation", &taus_byVLooseIsolation, &b_taus_byVLooseIsolation);
   fChain->SetBranchAddress("taus_byLooseIsolation", &taus_byLooseIsolation, &b_taus_byLooseIsolation);
   fChain->SetBranchAddress("taus_byMediumIsolation", &taus_byMediumIsolation, &b_taus_byMediumIsolation);
   fChain->SetBranchAddress("taus_byTightIsolation", &taus_byTightIsolation, &b_taus_byTightIsolation);
   fChain->SetBranchAddress("taus_byVLooseIsolationDeltaBetaCorr", &taus_byVLooseIsolationDeltaBetaCorr, &b_taus_byVLooseIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("taus_byLooseIsolationDeltaBetaCorr", &taus_byLooseIsolationDeltaBetaCorr, &b_taus_byLooseIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("taus_byMediumIsolationDeltaBetaCorr", &taus_byMediumIsolationDeltaBetaCorr, &b_taus_byMediumIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("taus_byTightIsolationDeltaBetaCorr", &taus_byTightIsolationDeltaBetaCorr, &b_taus_byTightIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("taus_signalPFChargedHadrCandsSize", &taus_signalPFChargedHadrCandsSize, &b_taus_signalPFChargedHadrCandsSize);
   fChain->SetBranchAddress("taus_muDecision", &taus_muDecision, &b_taus_muDecision);
   fChain->SetBranchAddress("taus_Nprongs", &taus_Nprongs, &b_taus_Nprongs);
   fChain->SetBranchAddress("Ntcmets", &Ntcmets, &b_Ntcmets);
   fChain->SetBranchAddress("tcmets_et", &tcmets_et, &b_tcmets_et);
   fChain->SetBranchAddress("tcmets_phi", &tcmets_phi, &b_tcmets_phi);
   fChain->SetBranchAddress("tcmets_ex", &tcmets_ex, &b_tcmets_ex);
   fChain->SetBranchAddress("tcmets_ey", &tcmets_ey, &b_tcmets_ey);
   fChain->SetBranchAddress("tcmets_sumEt", &tcmets_sumEt, &b_tcmets_sumEt);
   fChain->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
   fChain->SetBranchAddress("tracks_chi2", &tracks_chi2, &b_tracks_chi2);
   fChain->SetBranchAddress("tracks_ndof", &tracks_ndof, &b_tracks_ndof);
   fChain->SetBranchAddress("tracks_chg", &tracks_chg, &b_tracks_chg);
   fChain->SetBranchAddress("tracks_pt", &tracks_pt, &b_tracks_pt);
   fChain->SetBranchAddress("tracks_px", &tracks_px, &b_tracks_px);
   fChain->SetBranchAddress("tracks_py", &tracks_py, &b_tracks_py);
   fChain->SetBranchAddress("tracks_pz", &tracks_pz, &b_tracks_pz);
   fChain->SetBranchAddress("tracks_eta", &tracks_eta, &b_tracks_eta);
   fChain->SetBranchAddress("tracks_phi", &tracks_phi, &b_tracks_phi);
   fChain->SetBranchAddress("tracks_d0dum", &tracks_d0dum, &b_tracks_d0dum);
   fChain->SetBranchAddress("tracks_dz", &tracks_dz, &b_tracks_dz);
   fChain->SetBranchAddress("tracks_vx", &tracks_vx, &b_tracks_vx);
   fChain->SetBranchAddress("tracks_vy", &tracks_vy, &b_tracks_vy);
   fChain->SetBranchAddress("tracks_vz", &tracks_vz, &b_tracks_vz);
   fChain->SetBranchAddress("tracks_numvalhits", &tracks_numvalhits, &b_tracks_numvalhits);
   fChain->SetBranchAddress("tracks_numlosthits", &tracks_numlosthits, &b_tracks_numlosthits);
   fChain->SetBranchAddress("tracks_d0dumErr", &tracks_d0dumErr, &b_tracks_d0dumErr);
   fChain->SetBranchAddress("tracks_dzErr", &tracks_dzErr, &b_tracks_dzErr);
   fChain->SetBranchAddress("tracks_ptErr", &tracks_ptErr, &b_tracks_ptErr);
   fChain->SetBranchAddress("tracks_etaErr", &tracks_etaErr, &b_tracks_etaErr);
   fChain->SetBranchAddress("tracks_phiErr", &tracks_phiErr, &b_tracks_phiErr);
   fChain->SetBranchAddress("tracks_highPurity", &tracks_highPurity, &b_tracks_highPurity);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("experimentType", &experimentType, &b_experimentType);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("orbitNumber", &orbitNumber, &b_orbitNumber);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("model_params", &model_params, &b_model_params);

}


double EventCalculator::getDeltaPhi(double a, double b) {
  return jmt::deltaPhi(a,b);
}

cellECAL::cellECAL(double a, double b, int c) {
  eta=a;
  phi=b;
  status=c;
}

triggerData::triggerData() : pass(false), prescale(0), version(0) 
{
}
