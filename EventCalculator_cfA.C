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

#include "inJSON2012.h"

#include "MiscUtil.cxx"

//#include <RooRealVar.h>
//#include <RooMinuit.h>
//#include <RooTransverseThrustVar.h>

#include <cassert>
#include <fstream>

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
  //default settings for the collection pointers
  //this is what used to be in 'InitializeStuff()'
/* cfA changes
  myJetsPF( &jet2),    //selectedPatJetsPF
  myJetsPFhelper( &jethelper2),
  myElectronsPF( &electron1),
  myElectronsPFhelper( &electronhelper1),
  myMuonsPF(&muon1),
  myMuonsRECO(&muon),
  myMuonsPFhelper(&muonhelper1),
  myMuonsRECOhelper(&muonhelper),
  myTausPF(&tau),
  myMETPF(&met1),
  myMETPFType1(&met4),
  myMETcalo(&met),
  myVertex(&vertex),
  myGenParticles(&genparticlehelperra2),
  myGenWeight(&geneventinfoproduct_weight),
  myJetsPF_temp(0),
  myMETPF_temp(0),
  myEDM_bunchCrossing ( &eventhelper_bunchCrossing),
  myEDM_event ( &eventhelper_event),
  myEDM_isRealData ( &eventhelper_isRealData),
  myEDM_luminosityBlock ( &eventhelper_luminosityBlock),
  myEDM_run ( &eventhelper_run),
*/
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
  fDataJetRes_(0),
  //FIXME CFA
  ResJetPar_(0), //new JetCorrectorParameters("START42_V13_AK5PFchs_L2L3Residual.txt") ),
  L3JetPar_(0),// new JetCorrectorParameters("START42_V13_AK5PFchs_L3Absolute.txt") ),
  L2JetPar_(0),//new JetCorrectorParameters("START42_V13_AK5PFchs_L2Relative.txt") ),
  L1JetPar_(0),//new JetCorrectorParameters("START42_V13_AK5PFchs_L1FastJet.txt") ),
  JetCorrector_(0 ),
  jecUnc_(0),//new JetCorrectionUncertainty("START42_V13_AK5PFchs_Uncertainty.txt")),
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

  //FIXME CFA
//   vPar_.clear(); //should be overkill
//   vPar_.push_back(*L1JetPar_);
//   vPar_.push_back(*L2JetPar_);
//   vPar_.push_back(*L3JetPar_);
//   vPar_.push_back(*ResJetPar_);
//   JetCorrector_ =  new FactorizedJetCorrector(vPar_);

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
  theJESNames_[kJESFLY]="JESFLY";

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


// Don wrote this...i hope it is not needed anymore [cfA conversion]
/*
void EventCalculator::checkConsistency() {

 //this is probably bad coding, but the only thing I can think of for this already horribly-constructed code
  std::string jettype = typeid( *myJetsPF ).name();
  std::string recopfjet = "jet_s";
  //std::string pf2patjet = "jet1_s";
  std::string pf2patjet = "jet2_s";
  if( theJetType_==kPF2PAT && std::string::npos != jettype.find(recopfjet)) { 
    //std::cout << "myJetsPF is pointing to recopfjet" << std::endl;
    std::cout << "ERROR: theJetType_ is set to kPF2PAT jets, while myJetsPF is pointing to something else. "
	      << "This will screw up (at least) the cleaning selection.  Aborting." << std::endl;
    assert(0);
  }
  if( theJetType_==kRECOPF &&  std::string::npos != jettype.find(pf2patjet)){
    //std::cout << "myJetsPF is pointing to pf2patjet" << std::endl;
    std::cout << "ERROR: theJetType_ is set to kRECOPF jets, while myJetsPF is pointing to something else. "
	      << "This will screw up (at least) the cleaning selection.  Aborting." << std::endl;
    assert(0);
  }

  std::string mettype = typeid( *myMETPF ).name();
  std::string rawpfmet = "met1_s";
  //std::string type1pfmet = "met2_s";
  std::string type1pfmet = "met4_s";
  if( theMETType_ == kPFMETTypeI &&  std::string::npos != jettype.find(type1pfmet)){ //FIXME why a reference to jettype here?
    std::cout << "ERROR: theMETType_ is set to PFMETTypeI, while myMETsPF is also pointing to type1pfmet.  "
	      << "This will double-correct the MET!  Aborting." << std::endl;
    assert(0);
  }

  //theScanType_ is set automatically now, so this is basically redundant
  if (theScanType_!=kNotScan && sampleName_.Contains("LM")) {cout<<"LM point is not a scan! Check theScanType_!"<<endl; assert(0);}
  if (theScanType_!=kmSugra && sampleName_.Contains("SUGRA")) {cout<<"mSugra is a scan! Check theScanType_!"<<endl; assert(0);}
  if (theScanType_!=kSMS && sampleName_.Contains("T1bbbb")) {cout<<"T1bbbb is a scan! Check theScanType_!"<<endl; assert(0);}
  if (theScanType_!=kSMS && sampleName_.Contains("T2bb")) {cout<<"T2bb is a scan! Check theScanType_!"<<endl; assert(0);}
  if (theScanType_!=kSMS && sampleName_.Contains("T2tt")) {cout<<"T2tt is a scan! Check theScanType_!"<<endl; assert(0);}


}
*/

void EventCalculator::stopTimer(const Long64_t ntotal) {
  TDatime stoptime; //default ctor is for current time
  double elapsed= stoptime.Convert() - starttime_->Convert();
  std::cout<<"events / time = "<<ntotal<<" / "<<elapsed<<" = "<<double(ntotal)/double(elapsed)<<" Hz"<<std::endl;
  
  delete starttime_;
  starttime_=0;
}

bool EventCalculator::isSampleRealData() {

  if (sampleName_.BeginsWith("HT_Run2012A-PromptReco")) return true;
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
  else  if ( getOptPiece("JES",opt)== theJESNames_[kJESFLY]) theJESType_ = kJESFLY;
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

  if (sampleName_.Contains("QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6")) assert(0); //FIXME CFA
//w *= (*myGenWeight);

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

bool EventCalculator::isGoodMuon(const unsigned int imuon, const bool disableRelIso, const float ptthreshold) {

  //stealing code from Keith
  const float beamx = beamSpot_x->at(0);   
  const float beamy = beamSpot_y->at(0);   
  const float d0 = pf_mus_tk_d0dum->at(imuon) - beamx*sin(pf_mus_tk_phi->at(imuon)) + beamy*cos(pf_mus_tk_phi->at(imuon));

  const float muRelIso = (pf_mus_chargedHadronIso->at(imuon) + pf_mus_photonIso->at(imuon) + pf_mus_neutralHadronIso->at(imuon) ) / pf_mus_pt->at(imuon) ;

  if (pf_mus_pt->at(imuon) >= ptthreshold
      && fabs(pf_mus_eta->at(imuon))<2.4
      && TMath::Nint(pf_mus_id_GlobalMuonPromptTight->at(imuon)) == 1
      && TMath::Nint(pf_mus_isTrackerMuon->at(imuon)) == 1 
      && TMath::Nint(pf_mus_tk_numvalhits->at(imuon)) >=11
      && TMath::Nint(pf_mus_tk_numvalPixelhits->at(imuon)) >= 1
      //&& fabs(myMuonsPF->at(imuon).dB) < 0.02
      && fabs(d0) < 0.02
      && fabs(pf_mus_tk_vz->at(imuon) - pv_z->at(0) ) <= 1.0
      && ((muRelIso <0.2 ||disableRelIso)
	  )) {
    return true;
  }
  
  return false;
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
bool EventCalculator::isGoodElectron(const unsigned int iele, const bool disableRelIso, const float ptthreshold) {

  const float  beamx = beamSpot_x->at(0);
  const float  beamy = beamSpot_y->at(0); 
  const float   d0 = pf_els_d0dum->at(iele) - beamx*sin(pf_els_phi->at(iele)) + beamy*cos(pf_els_phi->at(iele));
  const float  elRelIso = (pf_els_chargedHadronIso->at(iele) + pf_els_photonIso->at(iele) + pf_els_neutralHadronIso->at(iele) ) / pf_els_pt->at(iele) ;

  const float sceta=pf_els_scEta->at(iele);
  
  if ( pf_els_pt->at(iele) >= ptthreshold
       && fabs( sceta) < 2.5 
       && !(fabs( sceta ) > 1.4442 
	    && fabs( sceta ) < 1.566)
       && TMath::Nint(pf_els_numlosthits->at(iele)) <= 1
       //&& fabs(myElectronsPF->at(iele).dB) < 0.02
       && fabs(d0) < 0.02
       && fabs(pf_els_vz->at(iele) - pv_z->at(0) ) < 1.0
       && (elRelIso <0.2 || disableRelIso)
       ) {
    return true;
  }
  
  return false;
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


bool EventCalculator::passHLT() {

  //for MC just return true
  if ( !isSampleRealData() ) return true;

  bool passtrigger=false;
  for (unsigned int itrig=0; itrig<trigger_name->size(); itrig++) {

    if ( TString(trigger_name->at(itrig).c_str()).BeginsWith("HLT_PFHT350_PFMET100_v")) {
      if ( TMath::Nint(trigger_prescalevalue->at(itrig)) == 1) {
	if ( TMath::Nint(trigger_decision->at(itrig))==1 ) passtrigger=true;
      }
    }
  }
  return passtrigger;
  
}

bool EventCalculator::passUtilityHLT(int &version, int &prescale) {

  //for MC just return true
  if ( !isSampleRealData() ) return true;

  bool passtrigger=false;
  for (unsigned int itrig=0; itrig<trigger_name->size(); itrig++) {
    //the flaw in this code is that it will return the results of only the *last* matching trigger
    if ( TString(trigger_name->at(itrig).c_str()).BeginsWith("HLT_PFHT350_v")) {
      TString vnumber=   TString(trigger_name->at(itrig).c_str()).Tokenize("_v")->At(2)->GetName();
      version=vnumber.Atoi();
      prescale =  TMath::Nint(trigger_prescalevalue->at(itrig));
      if ( TMath::Nint(trigger_decision->at(itrig))==1 ) passtrigger=true;
      //      cout<<version<<" "<<prescale<<" "<<passtrigger<<endl;
    }
  }

  return passtrigger;
 
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

int EventCalculator::countGoodPV() {
  
  int ngood=0;
  for (unsigned int ipv = 0; ipv<pv_x->size(); ipv++) {
    if ( ! pv_isFake->at(ipv) ) {
      if ( fabs(pv_z->at(ipv)) < 24 ) {
	if (sqrt(pow(pv_x->at(ipv),2.0)+pow(pv_y->at(ipv),2.0))  < 2 ) {
	  if ( pv_ndof->at(ipv) > 4 ) {
	    ++ngood;
	  }
	}
      }
    }
  }

  return ngood; 
}

bool EventCalculator::passPV() {
  //consider only the first PV

  bool npass=false;
  const unsigned int ipv=0;

  if ( ! pv_isFake->at(ipv) ) {
    if ( fabs(pv_z->at(ipv)) < 24 ) {
      if (sqrt(pow(pv_x->at(ipv),2.0)+pow(pv_y->at(ipv),2.0))  < 2 ) {
	if ( pv_ndof->at(ipv) > 4 ) {
	  npass=true;
	}
      }
    }
  }
    
  return npass; 

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
    assert( theMETType_ != kPFMETTypeI); //not implemented (or at least not checked)
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
    assert( theMETType_ != kPFMETTypeI); //not implemented (or at least not checked)
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


// JMT -- do we need this anymore????
bool EventCalculator::passCut(const TString cutTag) {

  if (cutTag=="cutInclusive") return true;
  if (cutTag=="cutTrigger" ) return passHLT();
  int dummyver = 0, dummyprs;
  if (cutTag=="cutUtilityTrigger" ) return ( passUtilityHLT(dummyver,dummyprs)>0  );

  if (cutTag=="cutPV") return passPV();
  
  if (cutTag=="cutHT") return getHT()>=400;
  
  if (cutTag=="cut3Jets") return (nGoodJets() >= 3);
  
  if (cutTag=="cutMuVeto") return (countMu()==0);
  if (cutTag=="cutEleVeto") return (countEle()==0);

  if (cutTag == "cutMET") {
    float mymet = getMET();
    return mymet >= 250;
  }

  if (cutTag == "cutDeltaPhi") {
    return ( getMinDeltaPhiMET(3) >= 0.3 );
  }
  if (cutTag== "cutDeltaPhiN"){
    return ( getMinDeltaPhiMETN(3) >= 4.0 );
  }
  if (cutTag=="cutDeltaPhiTaus"){
    return ( getMinDeltaPhiMETTaus() > 0.3 );
  }

  int nbtags = nGoodBJets();
  if (cutTag == "cut1b") return nbtags >=1;
  if (cutTag == "cut2b") return nbtags >=2;
  if (cutTag == "cut3b") return nbtags >=3;
  if (cutTag == "cutEq1b") return nbtags == 1;
  
  //jmt -- if we get here then i think something is wrong!
  assert(0);

  return true;
}

// JMT -- do we need this anymore?
bool EventCalculator::setCutScheme() {

  //if (cutscheme == nCutSchemes) return false;
  //theCutScheme_ = cutscheme;

  cutNames_.clear();
  cutTags_.clear();

  cutTags_.push_back("cutInclusive");cutNames_[ cutTags_.back()] = "Inclusive";
  
  //cutTags_.push_back("cut2SUSYb"); cutNames_[ cutTags_.back()] = "==2SUSYb";
  //cutTags_.push_back("cut4SUSYb"); cutNames_[ cutTags_.back()] = "==4SUSYb";
  
  cutTags_.push_back("cutTrigger"); cutNames_[cutTags_.back()]="Trigger";
  cutTags_.push_back("cutPV"); cutNames_[cutTags_.back()]="PV";
  cutTags_.push_back("cutHT"); cutNames_[cutTags_.back()]="HT";
  cutTags_.push_back("cut3Jets"); cutNames_[cutTags_.back()]=">=3Jets";
  //cutTags_.push_back("cutJetPt1");  cutNames_[cutTags_.back()]="JetPt1";
  //cutTags_.push_back("cutJetPt2");  cutNames_[cutTags_.back()]="JetPt2";
  //cutTags_.push_back("cutJetPt3");  cutNames_[cutTags_.back()]="JetPt3";
  
  cutTags_.push_back("cutEleVeto");  cutNames_[cutTags_.back()]="EleVeto";
  cutTags_.push_back("cutMuVeto");  cutNames_[cutTags_.back()]="MuVeto";
  //cutTags_.push_back("cutTauVeto");  cutNames_[cutTags_.back()]="TauVeto";
  
  cutTags_.push_back("cutMET");  cutNames_[cutTags_.back()]="MET";
  
  //cutTags_.push_back("cutDeltaPhi"); cutNames_[cutTags_.back()]="DeltaPhi";
  cutTags_.push_back("cutDeltaPhiN"); cutNames_[cutTags_.back()]="DeltaPhiN";
  //cutTags_.push_back("cutCleaning"); cutNames_[cutTags_.back()]="TailCleaning";
  
  cutTags_.push_back("cut1b"); cutNames_[cutTags_.back()]=">=1b";
  cutTags_.push_back("cut2b"); cutNames_[cutTags_.back()]=">=2b";
  cutTags_.push_back("cut3b"); cutNames_[cutTags_.back()]=">=3b";

  //cutTags_.push_back("cutCleaning"); cutNames_[cutTags_.back()]="TailCleaning";

  return true;
}

//JMT -- do we need this anymore?
void EventCalculator::setIgnoredCut(const TString cutTag) {

  bool ok=false;
  for (unsigned int i=0; i<cutTags_.size() ; i++) {
    if (cutTags_[i] == cutTag) {ok=true; break;}
  }
  if (!ok) {cout<<"Invalid cutIndex"<<endl; return;}

  ignoredCut_.push_back(cutTag);

}

//JMT -- do we need this anymore?
void EventCalculator::setRequiredCut(const TString cutTag) {

  bool ok=false;
  for (unsigned int i=0; i<cutTags_.size() ; i++) {
    if (cutTags_[i] == cutTag) {ok=true; break;}
  }
  if (!ok) {cout<<"Invalid cutIndex"<<endl; return;}

  requiredCut_.push_back(cutTag);

}

//JMT -- do we need this anymore?
bool EventCalculator::cutRequired(const TString cutTag) { //should put an & in here to improve performance

  //check if we are *ignoring* this cut (for N-1 studies, etc)
  for (unsigned int i = 0; i< ignoredCut_.size() ; i++) {
    if ( cutTag == ignoredCut_.at(i) ) return false;
  }

  //check if we are *requiring* this cut (special from the normal scheme)
  for (unsigned int i = 0; i< requiredCut_.size() ; i++) {
    if ( cutTag == requiredCut_.at(i) ) return true;
  }

  bool cutIsRequired=false;


  if      (cutTag == "cutInclusive")  cutIsRequired =  true;

  //else if (cutTag == "cut2SUSYb")  cutIsRequired = false;
  //else if (cutTag == "cut4SUSYb")  cutIsRequired = false;

  else if (cutTag == "cutTrigger")  cutIsRequired = true;
  else if (cutTag == "cutPV")  cutIsRequired =  true;
  else if (cutTag == "cutHT")  cutIsRequired =  true;
  else if (cutTag == "cut3Jets")  cutIsRequired =  true;
  //else if (cutTag == "cutJetPt1")  cutIsRequired =  false;
  //else if (cutTag == "cutJetPt2")  cutIsRequired =  false;
  //else if (cutTag == "cutJetPt3")  cutIsRequired =  false;
  else if (cutTag == "cutMET")  cutIsRequired =  true; 
  else if (cutTag == "cutMuVeto") cutIsRequired =  true;
  else if (cutTag == "cutEleVeto") cutIsRequired =  true;
  //else if (cutTag == "cutTauVeto") cutIsRequired =  false; //not required for now
  //else if (cutTag == "cutDeltaPhi") cutIsRequired =  true;
  else if (cutTag == "cutDeltaPhiN") cutIsRequired =  true;
  //else if (cutTag == "cutDeltaPhiTaus") cutIsRequired =  false;
  else if (cutTag == "cut1b") cutIsRequired =  true;
  else if (cutTag == "cut2b") cutIsRequired =  true;
  else if (cutTag == "cut3b") cutIsRequired =  true;
  //else if (cutTag == "cutCleaning") cutIsRequired = true;
  //else assert(0);

  return cutIsRequired;
}

//JMT -- do we need this anymore?
int EventCalculator::Cut()
{
  //be careful...not every function that makes cuts uses this! e.g. ::cutflow()

  for (unsigned int i=0; i< cutTags_.size(); i++) {
    if (cutRequired( cutTags_[i] ) && !passCut( cutTags_[i]) ) return -1;
  }
  
  return 1;
}

//JMT -- do we need this anymore?
void EventCalculator::resetIgnoredCut() {
  ignoredCut_.clear();
}

//JMT -- do we need this anymore?
void EventCalculator::resetRequiredCut() {
  requiredCut_.clear();
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
						  float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith) {
  double deltaT = getDeltaPhiMETN_electron_deltaT(ielectron, otherpt, othereta, otherid, dataJetRes, keith);

  //calculate deltaPhiMETN
  double dp =  getDeltaPhi(pf_els_phi->at(ielectron), getMETphi());
  double dpN = dp / atan2(deltaT, getMET());
  
  return dpN;
}
//the imuon is a different logic than is used for the jet index in getDeltaPhiMETN()
//it usable directly on the muon list. we have already checked that it is a good muon
double EventCalculator::getDeltaPhiMETN_muon( const unsigned int imuon, 
						  float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith) {
  double deltaT = getDeltaPhiMETN_muon_deltaT(imuon, otherpt, othereta, otherid, dataJetRes, keith);

  //calculate deltaPhiMETN
  double dp =  getDeltaPhi(pf_mus_phi->at(imuon), getMETphi());
  double dpN = dp / atan2(deltaT, getMET());
  
  return dpN;
}

double EventCalculator::getDeltaPhiMETN( unsigned int goodJetN, float mainpt, float maineta, bool mainid,
					 float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith) {//Ben
  
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
  double dpN = dp / atan2(deltaT, getMET());
  
  return dpN;
}


double EventCalculator::getMinDeltaPhiMETN(unsigned int maxjets, float mainpt, float maineta, bool mainid, 
					   float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith, bool includeLeptons) {//Ben
  
  double mdpN=1E12;
  
  for (unsigned int i=0; i<maxjets; i++) {
    if(i>=nGoodJets()) break;
    double dpN =  getDeltaPhiMETN(i, mainpt, maineta, mainid, otherpt, othereta, otherid, dataJetRes, keith);//i is for i'th *good* jet, starting at i=0. returns -99 if bad jet.
    if (dpN>=0 && dpN<mdpN) mdpN=dpN;//checking that dpN>=0 shouldn't be necessary after break statement above, but add it anyway 
  }

  //this option causes leptons to be treated like jets in calculating DeltaPhi(lepton, MET). the resolution term remains jets only
  if (includeLeptons) { 

    double  mdpN_e=2e12;
    double  mdpN_m=2e12;
    //need to get DeltaPhiMETN for any good leptons in the event
    for (unsigned int i=0; i <pf_els_pt->size() ; i++) {
      if (isGoodElectron(i,false,10)) { //pt threshold hard-coded to 10
	double  dPN_e=  getDeltaPhiMETN_electron(i, otherpt, othereta, otherid, dataJetRes, keith);
	if (dPN_e>=0 && dPN_e<mdpN_e) mdpN_e=dPN_e;
      }
    }
    for (unsigned int i=0; i <pf_mus_pt->size() ; i++) {
      if (isCleanMuon(i,10)) { //pt threshold hard-coded to 10
	double  dPN_m=  getDeltaPhiMETN_muon(i, otherpt, othereta, otherid, dataJetRes, keith);
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

//FIXME CFA
int EventCalculator::doPBNR() {
  return 1;

  /*
  bool nhBad=false;
  bool phBad=false;
  
  for (unsigned int i=0; i< jets_AK5PF_pt->size(); i++) {
    if (getJetPt(i) >30 ) {
      if (myJetsPF->at(i).neutralHadronEnergyFraction >0.90) nhBad=true;
      if (myJetsPF->at(i).photonEnergy / myJetsPF->at(i).energy > 0.95) phBad=true;
    }
  }
  
  if (nhBad && phBad) return -3;
  else if (phBad) return -2;
  else if (nhBad) return -1;
  return 1;
  */
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

unsigned int EventCalculator::nGoodBJets( BTaggerType btagger) {
  unsigned int nb=0;
  for (unsigned int i = 0; i < jets_AK5PF_pt->size(); ++i) {
    if (isGoodJet(i,30) ) {
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
  assert(theJetType_ == kPF2PAT);

  //start with uncorrected MET
  if (theMETType_== kPFMET ) {
    myMET = pfmets_et->at(0);
    myMETphi = pfmets_phi->at(0);
  }
  else {assert(0);} //this code isn't set up to handle corrected MET

  float myMETx = myMET * cos(myMETphi);
  float myMETy = myMET * sin(myMETphi);

  //  float a,b; getCorrectedMET(a,b);
  //  cout<<a<<"\t";
 
  //first remove jets from MET
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
    //    if (isCleanJet(i) ){//this only has an effect for recopfjets    
      //remove the uncorrected jet pT from the raw MET
    if ( jets_AK5PF_pt->at(i) >10 ) {
      myMETx += getUncorrectedJetPt(i) * cos(jets_AK5PF_phi->at(i));
      myMETy += getUncorrectedJetPt(i) * sin(jets_AK5PF_phi->at(i));
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

  //  if (sqrt(myMETx*myMETx + myMETy*myMETy) >50) {
  //    cout<<sqrt(myMETx*myMETx + myMETy*myMETy)<<"\t"<<atan2(myMETy,myMETx)<<endl; dumpEvent();
  //  }

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
    //    if (isCleanJet(i) ){//this only has an effect for recopfjets    
      //remove the corrected jet pT from MET
    if (jets_AK5PF_pt->at(i) >10 ) {
      myMETx -= getUncorrectedJetPt(i) * cos(jets_AK5PF_phi->at(i));
      myMETy -= getUncorrectedJetPt(i) * sin(jets_AK5PF_phi->at(i));
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
    //numbers from JME-10-014-pas
    //i have done my own averaging of Table 1
    //these numbers are only valid for PF Jets
    float m = 0;
    if (theJERType_ == kJERup) m=1;
    else if (theJERType_==kJERdown) m=-1;

    if (abseta < 1.1) {
      return 0.06303 + m*0.02478;
    }
    else if (abseta <1.7 && abseta>=1.1) {
      return 0.08467 + m*0.033169;
    }
    else if (abseta <2.3 && abseta>=1.7) {
      return  0.02557 + m*0.049505;
    }
    else if (abseta<=5 && abseta>= 2.3) {
      return 0.15819 + m*0.0663391;
    }
    else { //eta>5 not given
      return 0;
    }
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

  /*  float uncertainty = -1000; */
  return 0; //FIXME CFA


/* for timing tests only
  if ( watch_==0 ) watch_ = new TStopwatch(); //ctor starts the timer
  else watch_->Continue(); //kFALSE == don't reset the timer
*/

/* FIXME CFA
  if ( theJESType_ == kJESup ) {
    uncertainty = myJetsPFhelper->at(ijet).jetUncPlus; //ntuple values are much faster than running this on the fly
//     if( fabs(myJetsPF->at(ijet).eta) <5 ){
//       jecUnc_->setJetEta(myJetsPF->at(ijet).eta);
//       jecUnc_->setJetPt(myJetsPF->at(ijet).pt); // here you must use the CORRECTED jet pt
//       uncertainty = jecUnc_->getUncertainty(true);

//       //      cout<<"[EventCalculator::getJESUncertainty] "<<myJetsPFhelper->at(ijet).jetUncPlus<<"\t"<<uncertainty<<endl;
//       //      if ( fabs(  myJetsPFhelper->at(ijet).jetUncPlus - uncertainty) >0.01) cout<<"[EventCalculator::getJESUncertainty] something Bad!"<<endl;
      
//     }
  }
  else if (theJESType_ == kJESdown) {
    uncertainty = myJetsPFhelper->at(ijet).jetUncMinus; //ntuple values are much faster than running this on the fly
//     if( fabs(myJetsPF->at(ijet).eta) <5 ){
//       jecUnc_->setJetEta(myJetsPF->at(ijet).eta);
//       jecUnc_->setJetPt(myJetsPF->at(ijet).pt); // here you must use the CORRECTED jet pt
//       uncertainty = jecUnc_->getUncertainty(false);
//     }
  }

  //  watch_->Stop(); //for timing tests only

  if (fabs(uncertainty) > 1) uncertainty=0; //important sanity check

  //  cout<<"\t JES unc ("<< myJetsPF->at(ijet).pt <<") = "<<uncertainty;//JMT DEBUG

  if (addL2L3toJES ) {
    JetCorrector_->setJetEta( myJetsPF->at(ijet).eta);
    JetCorrector_->setJetPt( myJetsPF->at(ijet).uncor_pt );
    JetCorrector_->setJetA( myJetsPF->at(ijet).jetArea );
    JetCorrector_->setRho( sdouble_kt6pfjets_rho_value ); 
    std::vector<float> factors = JetCorrector_->getSubCorrections();
    float L2L3Residual = factors[3]/factors[2] - 1;

    uncertainty = sqrt(uncertainty*uncertainty + pow(L2L3Residual,2) );
    //    cout<<" ++ "<<L2L3Residual<<" = "<<uncertainty<<endl; //JMT DEBUG
  }

  return uncertainty;
*/
}

float EventCalculator::getJetPt( unsigned int ijet, bool addL2L3toJES ) {

  //i am going to write this to allow simultaneous use of JER and JES
  //(hence use if repeated 'if' instead of 'else if')

  float pt = jets_AK5PF_pt->at(ijet);
  //if ( theJESType_ == kJES0 && theJERType_ == kJER0) return pt;

  if (theJESType_ == kJESFLY) {
    assert(0); //FIXME CFA
    /* FIXME CFA
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
    */
    //for study of L2L3Res
	    //if (pt>30)  cout<<ijet<<" pT L2L3Residual Unc+ Unc- : "<<pt<<"\t"<<1-L2L3Residualonly<<" "<<myJetsPFhelper->at(ijet).jetUncPlus<<" "<<myJetsPFhelper->at(ijet).jetUncMinus <<" "<< sqrt( pow(1-L2L3Residualonly,2) + pow(myJetsPFhelper->at(ijet).jetUncPlus,2) ) <<endl;

  }

  
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
  //FIXME CFA
  unsigned int ngood=0;
/*
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {

    if( isSampleRealData() ) return 0;

    bool pass=false;
    pass = isGoodJet(i);

    if (pass ) {
      ngood++;
      if (ngood==n) return myJetsPF->at(i).genJet_phi;
    }
  }
*/
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

 return 99;//FIXME CFA
/*
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
*/
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

  if (jets_AK5PF_partonFlavour==0) return 0; //FIXME CFA

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

//FIXME CFA
int  EventCalculator::jetChargedHadronMultOfN(unsigned int n) {

  unsigned int ngood=0;
/*
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
    bool pass=false;
    pass = isGoodJet(i);

    if (pass ) {
      ++ngood;
      if (ngood==n) return  TMath::Nint(myJetsPF->at(i).chargedHadronMultiplicity);
    }
  }
*/
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

//FIXME CFA
float  EventCalculator::bjetChargedHadronFracOfN(unsigned int n) {

  unsigned int ngood=0;
/*
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {

    bool pass=false;
    pass = (isGoodJet30(i) && passBTagger(i));

    if (pass ) {
      ngood++;
      if (ngood==n) return   myJetsPF->at(i).chargedHadronEnergyFraction;
    }
  }
*/
  return 0;
}

int  EventCalculator::bjetChargedHadronMultOfN(unsigned int n) {
  //FIXME CFA
/*
  unsigned int ngood=0;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {

    bool pass=false;
    pass = (isGoodJet30(i) && passBTagger(i));

    if (pass ) {
      ngood++;
      if (ngood==n) return   TMath::Nint(myJetsPF->at(i).chargedHadronMultiplicity);
    }
  }
*/
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

//FIXME CFA
unsigned int EventCalculator::findSUSYMaternity( unsigned int k ) {

  int numSUSY=0;

  return numSUSY; //FIXME CFA

//   if (TMath::Nint(myGenParticles->at(k).firstMother) <0) {/*cout<<"got no mother!"<<endl;*/ return numSUSY;}

//   int motherId = abs(TMath::Nint( myGenParticles->at(  TMath::Nint(myGenParticles->at(k).firstMother)).pdgId )); //get rid of minus signs!
//   if ( (motherId>= 1000001 && motherId <=1000039) || (motherId>= 2000001 && motherId <=2000015)) {numSUSY++;}
//   else if (motherId == 21) {/*cout<<"Found gluon splitting!"<<endl;*/}
//   else if ( int(k) == TMath::Nint(myGenParticles->at(k).firstMother)) {} //particle is its own mother. do nothing. (avoids infinite loop)
//   else { numSUSY+= findSUSYMaternity(TMath::Nint(myGenParticles->at(k).firstMother)  ); }

//   return numSUSY;

}

unsigned int EventCalculator::getSUSYnb(std::vector<unsigned int> &susyb_index) {

  //  for (unsigned int k = 0; k<myGenParticles->size(); k++) {
    //for debugging
  //  cout<<k<<"\t"<<TMath::Nint(myGenParticles->at(k).pdgId )<<"\t"<< TMath::Nint(myGenParticles->at(k).firstMother)<<"\t"<<TMath::Nint(myGenParticles->at(k).lastMother)<<"\t"<<TMath::Nint(myGenParticles->at(k).status ) <<endl;
  //  }
  
  unsigned int SUSY_nb=0;
  /* FIXME CFA
  for (unsigned int k = 0; k<myGenParticles->size(); k++) {
    //important to require status 3!
    if ( abs(TMath::Nint(myGenParticles->at(k).pdgId )) == 5 && TMath::Nint(myGenParticles->at(k).status)==3 ) { //find b quark
      unsigned int nSUSY = findSUSYMaternity( k );      
      if (nSUSY>0) {SUSY_nb++; susyb_index.push_back(k);}
    }
  }
  */
  //cout<<"SUSY_nb = "<<SUSY_nb<<endl;
  return SUSY_nb;

}

//in case you don't care about the indices of the b quarks
unsigned int EventCalculator::getSUSYnb() {
  std::vector<unsigned int> dummy;
  return  getSUSYnb(dummy);

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


//FIXME CFA
std::pair<int,int> EventCalculator::getSMSmasses() {
  assert(theScanType_==kSMS);

  //our first pass at the official T1bbbb sample had bogus mGL, mLSP. we're not going to deal with that here.

  //Don says that the squark mass for T2qq is stored in the mGL field

  return  make_pair(0,0);  //FIXME CFA
  //  return make_pair( TMath::Nint(eventlhehelperextra_mGL), TMath::Nint(eventlhehelperextra_mLSP));
}

double EventCalculator::checkPdfWeightSanity( double a) {

  if ( std::isnan(a) ) {cout<<"Found PDF NaN!"<<endl; return 1;}
  if ( a<0 ) return 1;
  if ( a>10) return 1; //i'm not sure if this is a good idea at all, nor do i know what threshold to set

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

double EventCalculator::getCrossSection(){
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV

  //const double bf = 0.32442;

  //7 TeV to be updated
  if (sampleName_.Contains("QCD_Pt-0to5_TuneZ2_7TeV_pythia6") )                    return 4.844e10;
  if (sampleName_.Contains("QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6") )                return .3321;
  if (sampleName_.Contains("QCD_Pt-120to170_TuneZ2_7TeV_pythia6") )                  return 1.151e5;
  if (sampleName_.Contains("QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6") )              return .01087;
  if (sampleName_.Contains("QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6") )             return 2.213e10;
  if (sampleName_.Contains("QCD_Pt-15to30_TuneZ2_7TeV_pythia6") )                    return 8.159e8;
  if (sampleName_.Contains("QCD_Pt-170to300_TuneZ2_7TeV_pythia6") )                  return 2.426e4;
  if (sampleName_.Contains("QCD_Pt-1800_TuneZ2_7TeV_pythia6") )                      return .0003575;
  if (sampleName_.Contains("QCD_Pt-300to470_TuneZ2_7TeV_pythia6") )                  return 1.168e3;
  if (sampleName_.Contains("QCD_Pt-30to50_TuneZ2_7TeV_pythia6") )                    return 5.312e7;
  if (sampleName_.Contains("QCD_Pt-470to600_TuneZ2_7TeV_pythia6") )                  return 7.022e1;
  if (sampleName_.Contains("QCD_Pt-50to80_TuneZ2_7TeV_pythia6") )                    return 6.359e6;
  if (sampleName_.Contains("QCD_Pt-5to15_TuneZ2_7TeV_pythia6") )                     return 3.675e10;
  if (sampleName_.Contains("QCD_Pt-600to800_TuneZ2_7TeV_pythia6") )                  return 1.555e1;
  if (sampleName_.Contains("QCD_Pt-800to1000_TuneZ2_7TeV_pythia6") )                 return 1.844;
  if (sampleName_.Contains("QCD_Pt-80to120_TuneZ2_7TeV_pythia6") )                   return 7.843e5;
  if (sampleName_.Contains("ww") )                                                   return 27.83;//from PREP
  if (sampleName_.Contains("wz") )                                                   return 10.47;//from PREP
  if (sampleName_.Contains("zz") )                                                   return 4.287;//from PREP
  //7TeV
  if (sampleName_.Contains("T_TuneZ2_s-channel_7TeV-powheg-tauola") )                return 3.19; 
  if (sampleName_.Contains("Tbar_TuneZ2_s-channel_7TeV-powheg-tauola") )             return 1.44; 
  if (sampleName_.Contains("T_TuneZ2_t-channel_7TeV-powheg-tauola") )                return 41.92;
  if (sampleName_.Contains("Tbar_TuneZ2_t-channel_7TeV-powheg-tauola") )             return 22.65;
  if (sampleName_.Contains("T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola") )            return 7.87; 
  if (sampleName_.Contains("Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola") )         return 7.87; 

  if (sampleName_.Contains("WJetsToLNu_TuneZ2Star_8TeV-madgraph") )               return 36257.2; //NNLO 8 TeV

  //7TeV
  if (sampleName_.Contains("WJetsToLNu_250_HT_300_TuneZ2_7TeV-madgraph-tauola") )    return 34.8;
  if (sampleName_.Contains("WJetsToLNu_300_HT_inf_TuneZ2_7TeV-madgraph-tauola") )    return 48.49;
  if (sampleName_.Contains("DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola") )          return 3048;
  if (sampleName_.Contains("LM9_SUSY_sftsht_7TeV-pythia6") )                         return 7.134 * 1.48; //from ProductionSpring2011 twiki

  //8 TeV  
  if (sampleName_.Contains("SUSY_LM9_sftsht_8TeV"))          return 9.287; //LO 8TeV from PREP
  if (sampleName_.Contains("TTJets_TuneZ2star_8TeV-madgraph-tauola") )         return 234;  //approx NNLO
 

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

  //8TeV
  if (sampleName_.Contains("TTJets_TuneZ2star_8TeV-madgraph-tauola") )              return "TTbarJets";

  if (sampleName_.Contains("zjets") )                                                return "Zinvisible";
  //if (sampleName_.Contains("wjets") )                                                return "WJets";
  if (sampleName_.Contains("WJetsToLNu_TuneZ2_7TeV-madgraph-tauola") )               return "WJets";
  if (sampleName_.Contains("WJetsToLNu_250_HT_300_TuneZ2_7TeV-madgraph-tauola") )    return "WJetsHT250";
  if (sampleName_.Contains("WJetsToLNu_300_HT_inf_TuneZ2_7TeV-madgraph-tauola") )    return "WJetsHT300";
  if (sampleName_.Contains("DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola") )          return "ZJets";
  if (sampleName_.Contains("ww") )                                                   return "WW";
  if (sampleName_.Contains("wz") )                                                   return "WZ";
  if (sampleName_.Contains("zz") )                                                   return "ZZ";
  //if (sampleName_.Contains("t_s-channel") )                                          return "SingleTop-sChannel";
  //if (sampleName_.Contains("tbar_s-channel") )                                       return "SingleTopBar-sChannel";
  //if (sampleName_.Contains("t_t-channel") )                                          return "SingleTop-tChannel";
  //if (sampleName_.Contains("tbar_t-channel") )                                       return "SingleTopBar-tChannel";
  //if (sampleName_.Contains("t_tW-channel") )                                         return "SingleTop-tWChannel";
  //if (sampleName_.Contains("tbar_tW-channel") )                                      return "SingleTopBar-tWChannel";
  if (sampleName_.Contains("T_TuneZ2_s-channel_7TeV-powheg-tauola") )                return "SingleTop-sChannel";
  if (sampleName_.Contains("Tbar_TuneZ2_s-channel_7TeV-powheg-tauola") )             return "SingleTopBar-sChannel";
  if (sampleName_.Contains("T_TuneZ2_t-channel_7TeV-powheg-tauola") )                return "SingleTop-tChannel";
  if (sampleName_.Contains("Tbar_TuneZ2_t-channel_7TeV-powheg-tauola") )             return "SingleTopBar-tChannel";
  if (sampleName_.Contains("T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola") )            return "SingleTop-tWChannel";
  if (sampleName_.Contains("Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola") )         return "SingleTopBar-tWChannel";



  if (sampleName_.Contains("LM9_SUSY_sftsht_7TeV-pythia6") )                         return "LM9";

  //if it isn't found, just use the full name 
  return sampleName_;

}

TString EventCalculator::getCutDescriptionString(){
  
  TString cut = "";
  cut += theBTaggerNames_[theBTaggerType_];
  cut += "_";

  if (theMETType_ ==kPFMET) {
  }
  else if (theMETType_ == kPFMETTypeI) {
    cut += "PFMETTypeI_";
  }
  else assert(0);

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

void EventCalculator::jjResonanceFinder(float & mjj1, float & mjj2) {//simple first try
  mjj1=0;
  mjj2=0;

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
    }

  }
  
  
}

double EventCalculator::getSMSScanCrossSection( const double mgluino) {
  //the SMS scan cross sections are not very important, so I have not bothered to fully implement this

  return 1; //FIXME need to implement T2bb etc
  //I don't plan to ever fix this...load these in drawReducedTrees now
}


void EventCalculator::loadJetTagEffMaps() {
  //load the sample's efficiency maps
  if (sampleName_.Contains("QCD"))
    f_tageff_ = new TFile("histos_btageff_csvm_qcd.root","READ");
  //else if(sampleName_.Contains("t_s-channel"))
  else if (sampleName_.Contains("T_TuneZ2_s-channel_7TeV-powheg-tauola") )
    f_tageff_ = new TFile("histos_btageff_csvm_singletop_s.root","READ");
  //else if(sampleName_.Contains("tbar_s-channel"))
  else if (sampleName_.Contains("Tbar_TuneZ2_s-channel_7TeV-powheg-tauola") )
    f_tageff_ = new TFile("histos_btageff_csvm_singletopbar_s.root","READ");
  //else if(sampleName_.Contains("t_t-channel"))
  else if (sampleName_.Contains("T_TuneZ2_t-channel_7TeV-powheg-tauola") )   
    f_tageff_ = new TFile("histos_btageff_csvm_singletop_t.root","READ");
  //else if(sampleName_.Contains("tbar_t-channel"))
  else if (sampleName_.Contains("Tbar_TuneZ2_t-channel_7TeV-powheg-tauola") )
    f_tageff_ = new TFile("histos_btageff_csvm_singletopbar_t.root","READ");
  //else if(sampleName_.Contains("t_tW-channel"))
  else if (sampleName_.Contains("T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola") )
    f_tageff_ = new TFile("histos_btageff_csvm_singletop_tW.root","READ");
  //else if(sampleName_.Contains("tbar_tW-channel"))
  else if (sampleName_.Contains("Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola") )
    f_tageff_ = new TFile("histos_btageff_csvm_singletopbar_tW.root","READ");
  else if (sampleName_.Contains("zjets"))   
    f_tageff_ = new TFile("histos_btageff_csvm_zinvisible.root","READ");
  else if (sampleName_.Contains("WJetsToLNu_250_HT_300_TuneZ2_7TeV-madgraph-tauola") )
    f_tageff_ = new TFile("histos_btageff_csvm_wjetsHT250.root","READ");
  else if (sampleName_.Contains("WJetsToLNu_300_HT_inf_TuneZ2_7TeV-madgraph-tauola") )    
    f_tageff_ = new TFile("histos_btageff_csvm_wjetsHT300.root","READ");
  else if (sampleName_.Contains("ww") )  
    f_tageff_ = new TFile("histos_btageff_csvm_ww.root","READ");
  else if (sampleName_.Contains("wz") )  
    f_tageff_ = new TFile("histos_btageff_csvm_wz.root","READ");
  else if (sampleName_.Contains("zz") )  
    f_tageff_ = new TFile("histos_btageff_csvm_zz.root","READ");
  else if (sampleName_.Contains("DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola") )
    f_tageff_ = new TFile("histos_btageff_csvm_dyjets.root","READ");
  else if(sampleName_.Contains("LM9")	
	  ||sampleName_.Contains("SUGRA")
	  ||sampleName_.Contains("T1bbbb") 
	  ||sampleName_.Contains("T2bb") 
	  ||sampleName_.Contains("T2tt") 
	  ||sampleName_.Contains("T1tttt") 
	  )   
    f_tageff_ = new TFile("histos_btageff_csvm_lm9.root","READ");
  else if (sampleName_.Contains("TTJets_TuneZ2_7TeV-madgraph-tauola_Fall11_v2") )
    f_tageff_ = new TFile("histos_btageff_csvm_ttbarjets_fall11.root","READ");
  else{ //if all else fails, use fall11 ttbar
    std::cout << "loadJetTagEffMaps: could not find sample - using Fall11 ttbar efficiencies." << std::endl;
    //f_tageff_ = new TFile("histos_btageff_csvm.root","READ");
    f_tageff_ = new TFile("histos_btageff_csvm_ttbarjets_fall11.root","READ");
  }
}

void EventCalculator::calculateTagProb(float &Prob0, float &ProbGEQ1, float &Prob1, float &ProbGEQ2, float &Prob2, float &ProbGEQ3,
				       const float extraSFb, const float extraSFc, const float extraSFl, BTagEffModifier modifier) {

  char btageffname[200], ctageffname[200], ltageffname[200];
  std::string sbtageff = "h_btageff";  std::string sctageff = "h_ctageff";  std::string sltageff = "h_ltageff";
  sprintf(btageffname,"%s",sbtageff.c_str());   
  sprintf(ctageffname,"%s",sctageff.c_str());   
  sprintf(ltageffname,"%s",sltageff.c_str());   
  TH1F * h_btageff  = (TH1F *)f_tageff_->Get(btageffname);
  TH1F * h_ctageff  = (TH1F *)f_tageff_->Get(ctageffname);
  TH1F * h_ltageff  = (TH1F *)f_tageff_->Get(ltageffname);

  //must initialize correctly
  Prob2 = 0;
  Prob1 = 0; ProbGEQ1 = 1; Prob0 = 1; ProbGEQ2 = 0;

  for (unsigned int ijet=0; ijet<jets_AK5PF_pt->size(); ++ijet) {
    float subprob1=0;
    if(isGoodJet30(ijet)){

      float effi = jetTagEff(ijet, h_btageff, h_ctageff, h_ltageff, extraSFb,extraSFc,extraSFl,modifier);
      //      cout<<"effi = "<<effi<<endl;
      Prob0 = Prob0* ( 1 - effi);
      
      double product = 1;
      for (unsigned int kjet=0; kjet<jets_AK5PF_pt->size(); ++kjet) {
	if(isGoodJet30(kjet)){
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
		if(isGoodJet30(jjet)){
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


//removed unused code getBTagIPWeight()

//get MC btag efficiency
float EventCalculator::jetTagEff(unsigned int ijet, TH1F* h_btageff, TH1F* h_ctageff, TH1F* h_ltageff,
				 const float extraSFb, const float extraSFc, const float extraSFl,//options that were added but aren't used much
				 BTagEffModifier modifier ) { 

  float tageff=0;
  const float pt = getJetPt(ijet);
  const float eta = fabs(jets_AK5PF_eta->at(ijet));
  int flavor = jets_AK5PF_partonFlavour->at(ijet);

  float noEfficiencyThreshold=350; //threshold for the highest SF bin, which we sometimes set to 0 efficiency

  if(isGoodJet30(ijet)) {

    if (theBTagEffType_ == kBTagEff04 || theBTagEffType_ == kBTagEffup4 || theBTagEffType_ == kBTagEffdown4) { //new 2012 BTV prescription 

      if (abs(flavor) ==4 || abs(flavor)==5) { //heavy flavor
	float errFactor = 1;
	float x=pt;
	// Tagger: CSVM within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
	if (x >670) { //use SF for 670 with twice the errors
	  x=670;
	  errFactor=2;
	}
	if (abs(flavor) == 4)  errFactor=2; //charm has double the errors   "SFc = SFb with twice the quoted uncertainty"

	float  SFb = 0.6981*((1.+(0.414063*x))/(1.+(0.300155*x)));

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
	  myErr *= errFactor; //high pT and charm get scaled up

	  if ((theBTagEffType_ == kBTagEffup4) || (modifier==kHFup)) SFb += myErr;
	  else if ((theBTagEffType_ == kBTagEffdown4) || (modifier==kHFdown)) SFb -= myErr;
	  else assert(0);
	}

	//	cout<<"jet flavor, pt, SF = "<<abs(flavor)<<" "<<pt<<" "<<SFb<<endl;
	if      (abs(flavor) == 5) tageff = extraSFb * SFb * h_btageff->GetBinContent( h_btageff->FindBin( pt ) );
	else if (abs(flavor) == 4) tageff = extraSFc * SFb * h_ctageff->GetBinContent( h_ctageff->FindBin( pt ) );
	else assert(0);
      }
      else { //light flavor [ see https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs.C ]
	float SF=0;
	if ( eta < 0.8 ) {
	  if       (theBTagEffType_ == kBTagEff04 && (modifier==kBTagModifier0 || modifier==kHFdown||modifier==kHFup))    SF = ((1.06182+(0.000617034*pt))+(-1.5732e-06*(pt*pt)))+(3.02909e-10*(pt*(pt*pt)));
	  else  if (theBTagEffType_ == kBTagEffdown4 || (modifier==kLFdown)) SF = ((0.972455+(7.51396e-06*pt))+(4.91857e-07*(pt*pt)))+(-1.47661e-09*(pt*(pt*pt)));
	  else  if (theBTagEffType_ == kBTagEffup4 || (modifier==kLFup))   SF = ((1.15116+(0.00122657*pt))+(-3.63826e-06*(pt*pt)))+(2.08242e-09*(pt*(pt*pt)));
	}
	else if (eta>=0.8 && eta<1.6) {
	  if       (theBTagEffType_ == kBTagEff04 && (modifier==kBTagModifier0 || modifier==kHFdown||modifier==kHFup))    SF = ((1.111+(-9.64191e-06*pt))+(1.80811e-07*(pt*pt)))+(-5.44868e-10*(pt*(pt*pt)));
	  else  if (theBTagEffType_ == kBTagEffdown4 || (modifier==kLFdown)) SF = ((1.02055+(-0.000378856*pt))+(1.49029e-06*(pt*pt)))+(-1.74966e-09*(pt*(pt*pt)));
	  else  if (theBTagEffType_ == kBTagEffup4 || (modifier==kLFup))   SF = ((1.20146+(0.000359543*pt))+(-1.12866e-06*(pt*pt)))+(6.59918e-10*(pt*(pt*pt)));
	}
	else if (eta>=1.6 && eta<=2.4) {
	  if       (theBTagEffType_ == kBTagEff04 && (modifier==kBTagModifier0 || modifier==kHFdown||modifier==kHFup))    SF = ((1.08498+(-0.000701422*pt))+(3.43612e-06*(pt*pt)))+(-4.11794e-09*(pt*(pt*pt)));
	  else  if (theBTagEffType_ == kBTagEffdown4 || (modifier==kLFdown)) SF = ((0.983476+(-0.000607242*pt))+(3.17997e-06*(pt*pt)))+(-4.01242e-09*(pt*(pt*pt)));
	  else  if (theBTagEffType_ == kBTagEffup4 || (modifier==kLFup))   SF = ((1.18654+(-0.000795808*pt))+(3.69226e-06*(pt*pt)))+(-4.22347e-09*(pt*(pt*pt)));
	}
	//design question -- what do to for jets at eta>2.4? assert? or return tageff=0?
	//i guess tageff 0 makes the most sense

	tageff = SF * extraSFl * h_ltageff->GetBinContent( h_ltageff->FindBin( pt )); 
	//	cout<<"jet flavor, pt, SF = "<<abs(flavor)<<" "<<pt<<" "<<SF<<endl;

      }

    }
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

/* FIXME CFA
void EventCalculator::dumpEvent () {

  cout<<" == jets"<<endl;
  for (unsigned int i=0; i<jets_AK5PF_pt->size(); i++) {
    //    if (isCleanJet(i) ){//this only has an effect for recopfjets    
    cout<<"\t"<<myJetsPF->at(i).pt<<" "<<myJetsPF->at(i).eta<<" "<<myJetsPF->at(i).phi<<"\t"<<isCleanJet(i)<<endl;
  }
  cout<<" == muons"<<endl;
  for ( unsigned int i = 0; i<myMuonsPF->size() ; i++) {
    cout<<"\t"<<myMuonsPF->at(i).pt<<" "<<myMuonsPF->at(i).eta<<" "<<myMuonsPF->at(i).phi<<"\t"<<isCleanMuon(i)<<endl;
  }
  cout<<" == electrons"<<endl;
  for ( unsigned int i = 0; i<myElectronsPF->size() ; i++) {
    cout<<"\t"<<myElectronsPF->at(i).pt<<" "<<myElectronsPF->at(i).eta<<" "<<myElectronsPF->at(i).phi<<"\t"<<isGoodElectron(i)<<endl;
  }

}
*/


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

void EventCalculator::reducedTree(TString outputpath) {

  //JSON file reading. For now the JSON file name must be hard-coded here
  vector< vector<int> > VRunLumi = MakeVRunLumi("Cert_190456-194479_8TeV_PromptReco_Collisions12_JSON.txt");

  //open output file
  TString outfilename="reducedTree";
  outfilename += ".";
  outfilename += getCutDescriptionString();
  outfilename += ".";

  outfilename += getSampleNameOutputString();
  outfilename+=".root";
  if (outputpath[outputpath.Length()-1] != '/') outputpath += "/";
  outfilename.Prepend(outputpath);
  TFile fout(outfilename,"RECREATE");

  // define the TTree
  TTree reducedTree("reducedTree","tree with minimal cuts");
  reducedTree.SetMaxTreeSize(9900000000LL); //increase maximum tree size to 9.9GB

  TRandom3 random(getSeed());

  // ~~~~~~~~~ declare reducedTree variables ~~~~~~~~~~
  //we're making an ntuple, so size matters -- use float not double		     
  double weight; //one exception to the float rule				     
  //  float btagIPweight;//, pfmhtweight; //
  float PUweight;
  float hltHTeff;
  float hltMHTeff;
  float hltMHTeffBNN,hltMHTeffBNNUp,hltMHTeffBNNDown;

  double scanCrossSection,scanCrossSectionPlus,scanCrossSectionMinus;
  int m0,m12;

  int W1decayType = -1, W2decayType = -1;
  int decayType = -1;;

  ULong64_t lumiSection, eventNumber, runNumber;
  //  float METsig;
  float ST, STeff, HT, MHT, MET, METphi, minDeltaPhi, minDeltaPhiAll, minDeltaPhiAll30,minDeltaPhi30_eta5_noIdAll;
  //  float correctedMET, correctedMETphi
  float caloMET;
  float METPFType1,METPFType1phi;
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

  float maxJetMis, max2JetMis, maxJetMisAll30, max2JetMisAll30;
  float maxJetFracMis, max2JetFracMis, maxJetFracMisAll30, max2JetFracMisAll30;
  float deltaPhiMETJetMaxMis, deltaPhiMETJetMaxMis30;

  float deltaPhiStar, deltaPhiStar_badjet_pt, deltaPhiStar_badjet_eta, deltaPhiStar_badjet_phi;

  float minDeltaPhiMetTau;

  //items to investigate possible heavy flavor understimate in SIG -- looking at semi leptonic decays, bs, etc
  float CSVout1, CSVout2, CSVout3; //CSV tagger output for the lead three b-tagged jets
  float minDeltaPhiAllb30, deltaPhib1, deltaPhib2, deltaPhib3;
  float minDeltaPhiMETMuonsAll;

  bool cutHT,cutPV,cutTrigger;
  bool cut3Jets,cutEleVeto,cutMuVeto,cutMET,cutDeltaPhi;

  bool	csctighthaloFilter;
  bool	eenoiseFilter;
  bool	greedymuonFilter;
  bool	hbhenoiseFilter;
  bool	inconsistentmuonFilter;
  bool	ra2ecaltpFilter;
  bool	scrapingvetoFilter;
  bool	trackingfailureFilter;
  bool	trackingfailureFilterPFLOW;
  bool  recovrechitFilter;
  //  bool	ra2ecalbeFilter;

  bool passCleaning;
  int PBNRcode;

  bool buggyEvent;

  bool isRealData;
  bool pass_utilityHLT;
  UInt_t prescaleUtilityHLT;
  UInt_t versionUtilityHLT;
  bool pass_utilityPrescaleModuleHLT;


  int njets, njets30, nbjets, ntruebjets, nElectrons, nMuons, nTaus;
  int nbjetsSSVM,nbjetsTCHET,nbjetsSSVHPT,nbjetsTCHPT,nbjetsTCHPM,nbjetsCSVM,nbjetsCSVL;
  int nElectrons5, nElectrons15,nElectrons20;
  int nMuons5, nMuons15,nMuons20;

  float bestZmass;
  float mjj1,mjj2,mjjdiff;

  float jetpt1, jetgenpt1, jetphi1, jetgenphi1, jeteta1, jetgeneta1, jetenergy1, bjetpt1, bjetphi1, bjeteta1, bjetenergy1;
  float jetpt2, jetgenpt2, jetphi2, jetgenphi2, jeteta2, jetgeneta2, jetenergy2, bjetpt2, bjetphi2, bjeteta2, bjetenergy2;
  float jetpt3, jetgenpt3, jetphi3, jetgenphi3, jeteta3, jetgeneta3, jetenergy3, bjetpt3, bjetphi3, bjeteta3, bjetenergy3;
  int jetflavor1, jetflavor2, jetflavor3, bjetflavor1, bjetflavor2, bjetflavor3;

  float jetchargedhadronfrac1, jetchargedhadronfrac2, jetchargedhadronfrac3, bjetchargedhadronfrac1, bjetchargedhadronfrac2, bjetchargedhadronfrac3;
  int jetchargedhadronmult1, jetchargedhadronmult2, jetchargedhadronmult3, bjetchargedhadronmult1, bjetchargedhadronmult2, bjetchargedhadronmult3;

  float eleet1, elephi1, eleeta1, muonpt1, muonphi1, muoneta1;
  float eleet2, elephi2, eleeta2, muonpt2, muonphi2, muoneta2;
  float muoniso1,muonchhadiso1,muonphotoniso1,muonneutralhadiso1;
  float taupt1, taueta1;
  float eleRelIso,muonRelIso;

  float recomuonpt1, recomuonphi1, recomuoneta1;
  float recomuonpt2, recomuonphi2, recomuoneta2;
  float recomuoniso1, recomuoniso2;
  float recomuonmindphijet1, recomuonmindphijet2;

  float rl,rMET;

  float MT_Wlep;
  float MT_Wlep5,MT_Wlep15;
  float wMass, topMass, wCosHel, topCosHel;

  int nGoodPV;

  int SUSY_nb;
  int SUSY_process;
  float SUSY_recoilPt;

  //new variables from Luke
/* jmt -- remove for now
  float lambda1_allJets;
  float lambda2_allJets;
  float determinant_allJets;
  float lambda1_allJetsPlusMET;
  float lambda2_allJetsPlusMET;
  float determinant_allJetsPlusMET;
  float lambda1_topThreeJets;
  float lambda2_topThreeJets;
  float determinant_topThreeJets;
  float lambda1_topThreeJetsPlusMET;
  float lambda2_topThreeJetsPlusMET;
  float determinant_topThreeJetsPlusMET;
*/
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
    for (int im0= 0; im0<=1500; im0+=5) {
      for (int im12= 0; im12<=1500; im12+=5) {
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
  if (theScanType_==kSMS) scanSMSngen = new TH2D("scanSMSngen","number of generated events",150,0,1500,150,0,1500); //mgluino,mLSP

  const  bool puReweightIs1D = true;//((theScanType_!=kNotScan) || sampleName_.Contains("QCD")); //CFA -- no 3D for now

  //initialize PU things
  std::vector< float > DataDist;
  std::vector< float > MCDist;
  for( int i=0; i<60; ++i) {
    //for PU_S3 (only in-time)
    //1D reweighting -- for QCD and Fastsim scans
    if (puReweightIs1D)    DataDist.push_back(pu::FirstData2012[i]);
    //for 3dPU reweighting 
    //    else DataDist2011.push_back(pu::TrueDist2011_f[i]);
    MCDist.push_back(pu::Summer2012[i]);
  } 

  reweight::LumiReWeighting * LumiWeights=0;
  // FIXME CFA
  //Lumi3DReWeighting * LumiWeights3D=0;
  if (puReweightIs1D)  LumiWeights = new reweight::LumiReWeighting( MCDist, DataDist );
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
  reducedTree.Branch("weight",&weight,"weight/D");
  reducedTree.Branch("scanCrossSection",&scanCrossSection,"scanCrossSection/D");
  reducedTree.Branch("scanCrossSectionPlus",&scanCrossSectionPlus,"scanCrossSectionPlus/D");
  reducedTree.Branch("scanCrossSectionMinus",&scanCrossSectionMinus,"scanCrossSectionMinus/D");
  reducedTree.Branch("runNumber",&runNumber,"runNumber/l");
  reducedTree.Branch("lumiSection",&lumiSection,"lumiSection/l");
  reducedTree.Branch("eventNumber",&eventNumber,"eventNumber/l");

  reducedTree.Branch("m0",&m0,"m0/I");
  reducedTree.Branch("m12",&m12,"m12/I");

  reducedTree.Branch("W1decayType",&W1decayType,"W1decayType/I");
  reducedTree.Branch("W2decayType",&W2decayType,"W2decayType/I");
  reducedTree.Branch("decayType",&decayType,"decayType/I");

  //  reducedTree.Branch("btagIPweight",&btagIPweight,"btagIPweight/F");
  //  reducedTree.Branch("pfmhtweight",&pfmhtweight,"pfmhtweight/F");
  reducedTree.Branch("PUweight",&PUweight,"PUweight/F");
  reducedTree.Branch("hltHTeff",&hltHTeff,"hltHTeff/F");
  reducedTree.Branch("hltMHTeff",&hltMHTeff,"hltMHTeff/F");
  reducedTree.Branch("hltMHTeffBNN",&hltMHTeffBNN,"hltMHTeffBNN/F");
  reducedTree.Branch("hltMHTeffBNNUp",&hltMHTeffBNNUp,"hltMHTeffBNNUp/F");
  reducedTree.Branch("hltMHTeffBNNDown",&hltMHTeffBNNDown,"hltMHTeffBNNDown/F");

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


  //should consider whether some of these should be killed off
  reducedTree.Branch("cutHT",&cutHT,"cutHT/O");
  reducedTree.Branch("cutPV",&cutPV,"cutPV/O");
  reducedTree.Branch("cutTrigger",&cutTrigger,"cutTrigger/O");
  reducedTree.Branch("cut3Jets",&cut3Jets,"cut3Jets/O");
  reducedTree.Branch("cutEleVeto",&cutEleVeto,"cutEleVeto/O");
  reducedTree.Branch("cutMuVeto",&cutMuVeto,"cutMuVeto/O");
  reducedTree.Branch("cutMET",&cutMET,"cutMET/O");
  reducedTree.Branch("cutDeltaPhi",&cutDeltaPhi,"cutDeltaPhi/O");
  
  reducedTree.Branch("csctighthaloFilter",&csctighthaloFilter,"csctighthaloFilter/O");
  reducedTree.Branch("eenoiseFilter",&eenoiseFilter,"eenoiseFilter/O");
  reducedTree.Branch("greedymuonFilter",&greedymuonFilter,"greedymuonFilter/O");
  reducedTree.Branch("hbhenoiseFilter",&hbhenoiseFilter,"hbhenoiseFilter/O");
  reducedTree.Branch("inconsistentmuonFilter",&inconsistentmuonFilter,"inconsistentmuonFilter/O");
  reducedTree.Branch("ra2ecaltpFilter",&ra2ecaltpFilter,"ra2ecaltpFilter/O");
  reducedTree.Branch("scrapingvetoFilter",&scrapingvetoFilter,"scrapingvetoFilter/O");
  reducedTree.Branch("trackingfailureFilter",&trackingfailureFilter,"trackingfailureFilter/O");
  reducedTree.Branch("trackingfailureFilterPFLOW",&trackingfailureFilterPFLOW,"trackingfailureFilterPFLOW/O");
  reducedTree.Branch("recovrechitFilter",&recovrechitFilter,"recovrechitFilter/O");
  //  reducedTree.Branch("ra2ecalbeFilter",&ra2ecalbeFilter,"ra2ecalbeFilter/O");
  reducedTree.Branch("passCleaning",&passCleaning,"passCleaning/O");
  reducedTree.Branch("PBNRcode",&PBNRcode,"PBNRcode/I");

  reducedTree.Branch("buggyEvent",&buggyEvent,"buggyEvent/O");

  reducedTree.Branch("nGoodPV",&nGoodPV,"nGoodPV/I");

  reducedTree.Branch("SUSY_nb",&SUSY_nb,"SUSY_nb/I");
  reducedTree.Branch("SUSY_process",&SUSY_process,"SUSY_process/I");
  reducedTree.Branch("SUSY_recoilPt",&SUSY_recoilPt,"SUSY_recoilPt/F");

  reducedTree.Branch("njets",&njets,"njets/I");
  reducedTree.Branch("njets30",&njets30,"njets30/I");
  reducedTree.Branch("nbjets",&nbjets,"nbjets/I");
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

  reducedTree.Branch("HT",&HT,"HT/F");
  reducedTree.Branch("ST",&ST,"ST/F"); //includes HT + leptons
  reducedTree.Branch("STeff",&STeff,"STeff/F"); //includes HT + leptons + MET
  reducedTree.Branch("MET",&MET,"MET/F");
  //  reducedTree.Branch("METsig",&METsig,"METsig/F");
  reducedTree.Branch("METphi",&METphi,"METphi/F");
  reducedTree.Branch("MHT",&MHT,"MHT/F");

  //  reducedTree.Branch("correctedMET",&correctedMET,"correctedMET/F");
  //  reducedTree.Branch("correctedMETphi",&correctedMETphi,"correctedMETphi/F");

  reducedTree.Branch("caloMET",&caloMET,"caloMET/F");
  reducedTree.Branch("METPFType1",&METPFType1, "METPFType1/F");
  reducedTree.Branch("METPFType1phi",&METPFType1phi, "METPFType1phi/F");
  
  reducedTree.Branch("bestWMass",&wMass,"bestWMass/F");
  reducedTree.Branch("bestTopMass",&topMass,"bestTopMass/F");
  reducedTree.Branch("topCosHel",&topCosHel,"topCosHel/F");
  reducedTree.Branch("WCosHel",&wCosHel,"WCosHel/F");
  reducedTree.Branch("MT_Wlep",&MT_Wlep, "MT_Wlep/F");
  reducedTree.Branch("MT_Wlep5",&MT_Wlep5, "MT_Wlep5/F");
  reducedTree.Branch("MT_Wlep15",&MT_Wlep15, "MT_Wlep15/F");

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

  reducedTree.Branch("eleet1",&eleet1,"eleet1/F");
  reducedTree.Branch("elephi1",&elephi1,"elephi1/F");
  reducedTree.Branch("eleeta1",&eleeta1,"eleeta1/F");
  reducedTree.Branch("muonpt1",&muonpt1,"muonpt1/F");
  reducedTree.Branch("muonphi1",&muonphi1,"muonphi1/F");
  reducedTree.Branch("muoneta1",&muoneta1,"muoneta1/F");
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

/*
  reducedTree.Branch("lambda1_allJets",&lambda1_allJets,"lambda1_allJets/F");
  reducedTree.Branch("lambda2_allJets",&lambda2_allJets,"lambda2_allJets/F");
  reducedTree.Branch("determinant_allJets",&determinant_allJets,"determinant_allJets/F");
  reducedTree.Branch("lambda1_allJetsPlusMET",&lambda1_allJetsPlusMET,"lambda1_allJetsPlusMET/F");
  reducedTree.Branch("lambda2_allJetsPlusMET",&lambda2_allJetsPlusMET,"lambda2_allJetsPlusMET/F");
  reducedTree.Branch("determinant_allJetsPlusMET",&determinant_allJetsPlusMET,"determinant_allJetsPlusMET/F");
  reducedTree.Branch("lambda1_topThreeJets",&lambda1_topThreeJets,"lambda1_topThreeJets/F");
  reducedTree.Branch("lambda2_topThreeJets",&lambda2_topThreeJets,"lambda2_topThreeJets/F");
  reducedTree.Branch("determinant_topThreeJets",&determinant_topThreeJets,"determinant_topThreeJets/F");
  reducedTree.Branch("lambda1_topThreeJetsPlusMET",&lambda1_topThreeJetsPlusMET,"lambda1_topThreeJetsPlusMET/F");
  reducedTree.Branch("lambda2_topThreeJetsPlusMET",&lambda2_topThreeJetsPlusMET,"lambda2_topThreeJetsPlusMET/F");
  reducedTree.Branch("determinant_topThreeJetsPlusMET",&determinant_topThreeJetsPlusMET,"determinant_topThreeJetsPlusMET/F");
*/
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

  const Long64_t nevents = chainA->GetEntries();
  const Long64_t neventsB = chainB->GetEntries();
  assert(nevents==neventsB);

  startTimer();
  // ~~~~ now start the real event loop
  for(Long64_t entry=0; entry < nevents; ++entry) {
    chainB->GetEntry(entry);
    chainA->GetEntry(entry);

    //some output to watch while it's running
    if (entry==0) {
      if ( !isSampleRealData() ) {
	cout << "MC xsec: " << getCrossSection() << endl;
      }
      else{
	cout << "This is data!"<< endl;
      }
      cout << "weight: "<< getWeight(nevents) << endl;
    }
    if (entry%100000==0 ) cout << "  entry: " << entry << ", percent done=" << (int)(entry/(double)nevents*100.)<<  endl;


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
      assert(0); // FIXME CFA

      //increment a 2d histogram of mGL, mLSP
      //we know 10k were generated everywhere, but what if we have failed jobs?
      scanSMSngen->Fill(m0,m12); //we can probably kill this 

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


    //very loose skim for reducedTrees (HT, trigger)
    ST = getST(30,10); //ST is always bigger than HT; reducing the jet cut to 30 will only make it larger
    if ( (passCut("cutTrigger") || passCut("cutUtilityTrigger")) && (ST>=400) ) {
      cutHT = passCut("cutHT"); //should always be true

      weight = getWeight(nevents);

      
      if (theScanType_!=kSMS) {
	scanCrossSection = getScanCrossSection(prodprocess,"");
	scanCrossSectionPlus = getScanCrossSection(prodprocess,"Plus");
	scanCrossSectionMinus = getScanCrossSection(prodprocess,"Minus");
      }
      else {
	scanCrossSection = getSMSScanCrossSection(m0); //gluino mass
	scanCrossSectionPlus = scanCrossSection;
	scanCrossSectionMinus = scanCrossSection; //there are no cross section errors for SMS
      }

      runNumber = getRunNumber();
      lumiSection = getLumiSection();
      eventNumber = getEventNumber();

      //stuff that we only do on data....
      bool      passJSON=true;
      if (isSampleRealData()) {
	//check the JSON file
	passJSON=  inJSON(VRunLumi,runNumber,lumiSection);
      }
      if (!passJSON) continue; //don't put in event that aren't in the json

      //if we are running over ttbar, fill info on decay mode
/* FIXME CFA
      if (sampleName_.Contains("ttjets_madgraph") || sampleName_.Contains("TTJets_TuneZ2_7TeV-madgraph-tauola_Fall11_v2")){
	int W1, W1daughter, W2, W2daughter;
	decayType = getTTbarDecayType(W1decayType, W2decayType, W1, W1daughter, W2, W2daughter,false); //FIXME CFA
      }
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

      HT=getHT();
      hltHTeff = getHLTHTeff(HT);
      PUweight =  puReweightIs1D ? getPUWeight(*LumiWeights) : 1;// FIXME CFA getPUWeight(*LumiWeights3D);
      //      btagIPweight = getBTagIPWeight();
      //      pfmhtweight = getPFMHTWeight();

      cutTrigger = passCut("cutTrigger");
      cutPV = passCut("cutPV");
      cut3Jets = passCut("cut3Jets");
      cutEleVeto = passCut("cutEleVeto");
      cutMuVeto = passCut("cutMuVeto");
      cutMET = passCut("cutMET");
      cutDeltaPhi = passCut("cutDeltaPhi");

      //FIXME CFA
/*
      calculateTagProb(prob0,probge1,prob1,probge2,prob2,probge3);

      calculateTagProb(prob0_HFplus,probge1_HFplus,prob1_HFplus,probge2_HFplus,prob2_HFplus,probge3_HFplus,1,1,1,kHFup);
      calculateTagProb(prob0_HFminus,probge1_HFminus,prob1_HFminus,probge2_HFminus,prob2_HFminus,probge3_HFminus,1,1,1,kHFdown);

      calculateTagProb(prob0_LFplus,probge1_LFplus,prob1_LFplus,probge2_LFplus,prob2_LFplus,probge3_LFplus,1,1,1,kLFup);
      calculateTagProb(prob0_LFminus,probge1_LFminus,prob1_LFminus,probge2_LFminus,prob2_LFminus,probge3_LFminus,1,1,1,kLFdown);
*/

      isRealData = isSampleRealData();
      int version = 0, prescale = 0;
      pass_utilityHLT = passUtilityHLT(version, prescale);
      prescaleUtilityHLT = prescale;
      versionUtilityHLT = version;

      pass_utilityPrescaleModuleHLT = passUtilityPrescaleModuleHLT();

      nGoodPV = countGoodPV();

      SUSY_nb = sampleIsSignal_ ? getSUSYnb() : 0;
      //      bjetSumSUSY[thispoint] += SUSY_nb;
      //if(SUSY_process==NotFound) cout<<"SUSY_nb = "<<SUSY_nb<<endl;

      njets = nGoodJets();
      njets30 = nGoodJets30();

      jjResonanceFinder(mjj1,mjj2);
      mjjdiff=fabs(mjj1-mjj2);

      ntruebjets = nTrueBJets();
      nbjets = nGoodBJets();
      nbjetsSSVM = nGoodBJets( kSSVM);
      nbjetsSSVHPT = nGoodBJets( kSSVHPT);
      nbjetsTCHET = nGoodBJets( kTCHET);
      nbjetsTCHPT = nGoodBJets( kTCHPT);
      nbjetsTCHPM = nGoodBJets( kTCHPM);
      nbjetsCSVM = nGoodBJets( kCSVM);
      nbjetsCSVL = nGoodBJets( kCSVL);

      CSVout1=bjetCSVOfN(1);
      CSVout2=bjetCSVOfN(2);
      CSVout3=bjetCSVOfN(3);
      minDeltaPhiAllb30=getMinDeltaPhiMET30(99,true);
      deltaPhib1=getDeltaPhiMET(1,30,true);
      deltaPhib2=getDeltaPhiMET(2,30,true); 
      deltaPhib3=getDeltaPhiMET(3,30,true);
      minDeltaPhiMETMuonsAll=getMinDeltaPhiMETMuons(99);

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
      MET=getMET();
      STeff = ST + MET;
      //      METsig = myMETPF->at(0).mEtSig; //FIXME hard coded for PF
      MHT=getMHT();
      METphi = getMETphi();

      //      getCorrectedMET(correctedMET, correctedMETphi);

      caloMET = mets_AK5_et->at(0); //is this right?
      METPFType1 = pfTypeImets_et->at(0);
      METPFType1phi = pfTypeImets_phi->at(0);

      minDeltaPhi = getMinDeltaPhiMET(3);
      minDeltaPhiAll = getMinDeltaPhiMET(99);
      minDeltaPhiAll30 = getMinDeltaPhiMET30(99);
      minDeltaPhi30_eta5_noIdAll = getMinDeltaPhiMET30_eta5_noId(99);
      maxDeltaPhi = getMaxDeltaPhiMET(3);
      maxDeltaPhiAll = getMaxDeltaPhiMET(99);
      maxDeltaPhiAll30 = getMaxDeltaPhiMET30(99);
      maxDeltaPhi30_eta5_noIdAll = getMaxDeltaPhiMET30_eta5_noId(99);
      
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

      int badjet;
      deltaPhiStar = getDeltaPhiStar(badjet);
      deltaPhiStar_badjet_pt = (badjet>=0) ? jets_AK5PF_pt->at(badjet) : -1;
      deltaPhiStar_badjet_phi = (badjet>=0) ? jets_AK5PF_phi->at(badjet) : -1;
      deltaPhiStar_badjet_eta = (badjet>=0) ? jets_AK5PF_eta->at(badjet) : -1;

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
      deltaPhiN1 = getDeltaPhiMETN(0);
      deltaPhiN2 = getDeltaPhiMETN(1);
      deltaPhiN3 = getDeltaPhiMETN(2);

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
      double effUp, effDown;
      hltMHTeffBNN = getHLTMHTeffBNN(MET, HT, nElectrons, nMuons, minDeltaPhiN, effUp, effDown);
      hltMHTeffBNNUp = effUp;
      hltMHTeffBNNDown = effDown;

      //alternatives for testing
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
      
      minDeltaPhiN_Luke = getMinDeltaPhiNMET(3);
      maxDeltaPhiN_Luke = getMaxDeltaPhiNMET(3);
      deltaPhiN1_Luke = getDeltaPhiNMET(0);
      deltaPhiN2_Luke = getDeltaPhiNMET(1);
      deltaPhiN3_Luke = getDeltaPhiNMET(2);
      
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
      
      minTransverseMETSignificance = getMinTransverseMETSignificance(3);
      maxTransverseMETSignificance = getMaxTransverseMETSignificance(3);
      transverseMETSignificance1 = getTransverseMETSignificance(0);
      transverseMETSignificance2 = getTransverseMETSignificance(1);
      transverseMETSignificance3 = getTransverseMETSignificance(2);
      
      MT_Wlep = getMT_Wlep();
      MT_Wlep5 = getMT_Wlep(5);
      MT_Wlep15 = getMT_Wlep(15);
      
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
      muonpt1 = muonPtOfN(1,5);
      muonphi1 = muonPhiOfN(1,5);
      muoneta1 = muonEtaOfN(1,5);
      muoniso1 = muonIsoOfN(1,5);
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


      csctighthaloFilter = cschalofilter_decision==1 ? true : false;
      eenoiseFilter = eenoisefilter_decision==1 ? true : false;
      greedymuonFilter = greedymuonfilter_decision==1 ? true : false;
      hbhenoiseFilter = (theScanType_!=kNotScan || (hbhefilter_decision==1)) ? true : false;
      inconsistentmuonFilter = inconsistentPFmuonfilter_decision==1 ? true : false;
      ra2ecaltpFilter =ecalTPfilter_decision==1 ? true:false; 
      scrapingvetoFilter = scrapingVeto_decision==1?true:false;
      trackingfailureFilter = true; //FIXME CFA (doesn't exist in cfA?)
      trackingfailureFilterPFLOW = trackingfailurefilter_decision ==1 ? true:false;
      recovrechitFilter = true; //FIXME CFA

      passCleaning = csctighthaloFilter && eenoiseFilter && greedymuonFilter && hbhenoiseFilter && inconsistentmuonFilter && ra2ecaltpFilter && scrapingvetoFilter && trackingfailureFilterPFLOW;

      PBNRcode = doPBNR();

      buggyEvent = inEventList(vrun, vlumi, vevent);
      
      //fill new variables from Luke
/*
      getSphericityJetMET(lambda1_allJets,lambda2_allJets,determinant_allJets,99,false);
      getSphericityJetMET(lambda1_allJetsPlusMET,lambda2_allJetsPlusMET,determinant_allJetsPlusMET,99,true);
      getSphericityJetMET(lambda1_topThreeJets,lambda2_topThreeJets,determinant_topThreeJets,3,false);
      getSphericityJetMET(lambda1_topThreeJetsPlusMET,lambda2_topThreeJetsPlusMET,determinant_topThreeJetsPlusMET,3,true);
*/
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


void EventCalculator::changeVariablesGenSmearing(TRandom3* random) {
  //this code is very much a work in progress
  //the goal is still far away...

  //goal: starting with gen jets, fill jet list from scratch with smeared gen jets
  //two ways to proceed:
  // (1) start with gen particle list, looking for (status 3?) quarks and gluons
  // (2) start with genPartons associated with reco'd jets
  //version (2) is where i will start
  /*
Let me explain it with a pseudo-algorithm.
1) take truly balanced QCD events (gen-level partons, or even better
rebalanced events like in R+S)

2) smear exactly 1 jet in a way that can give a large mismeasurement
[maybe a double or triple Gaussian smearing, with the tail functions having
a very large width]

3) smear the other jets with a 0.10*pT Gaussian
[simple Gaussian smearing, so that virtually no jets get a large
mismeasurement]

4) compute MET, DeltaPhiN and study the r(MET) plot

  */

 //  if(recalculatedVariables_) return;
//   recalculatedVariables_=true;
//   myJetsPF_temp                   = myJetsPF;		  
//   myMETPF_temp		          = myMETPF;		  

//   myJetsPF                = new std::vector<jet2_s>;	   //jmt -- switch to PF2PAT
//   myMETPF		  = new std::vector<met1_s>(*myMETPF_temp);	  


/*
  //as a cross-check, let's calculate a jet-only gen-level MET
//this didn't work very well....
  double genMETx=0;
  double genMETy=0;
  for (unsigned int i=0; i<myJetsPF_temp->size(); i++) {
    double genpt =  jet2_genParton_pt .at( i);
    double geneta=  jet2_genParton_eta.at( i);
    double genphi=  jet2_genParton_phi.at( i);

    if (fabs(geneta)<5) {
      genMETx -= genpt * cos(genphi);
      genMETy -= genpt * sin(genphi);

      myJetsPF->push_back( myJetsPF_temp->at(i)); //temp
    }
  }
  cout<<"genMET = "<<sqrt(genMETx*genMETx + genMETy*genMETy)<<" "<<met1_genMET_et.at(0)<<endl;
*/

  cout<<"=="<<endl;
/*
//start with MET
  double METx = myMETPF_temp->at(0).pt * cos(myMETPF_temp->at(0).phi);
  double METy = myMETPF_temp->at(0).pt * sin(myMETPF_temp->at(0).phi);
//remove the jets from the MET  
  for (unsigned int i=0; i<myJetsPF_temp->size(); i++) {
    METx += myJetsPF_temp->at(i).uncor_pt * cos(myJetsPF_temp->at(i).uncor_phi);
    METy += myJetsPF_temp->at(i).uncor_pt * sin(myJetsPF_temp->at(i).uncor_phi);
  }
  cout<<myMETPF_temp->at(0).pt<<" "<<sqrt(METx*METx + METy*METy)<<flush;
  for (unsigned int i=0; i<myElectronsPF->size(); i++) {
    METx += myElectronsPF->at(i).pt * cos(myElectronsPF->at(i).phi);
    METy += myElectronsPF->at(i).pt * sin(myElectronsPF->at(i).phi);
  }
  for (unsigned int i=0; i<myMuonsPF->size(); i++) {
    METx += myMuonsPF->at(i).pt * cos(myMuonsPF->at(i).phi);
    METy += myMuonsPF->at(i).pt * sin(myMuonsPF->at(i).phi);
  }

  cout<<" "<<sqrt(METx*METx + METy*METy)<<endl;
  
  //now try the hail mary approach

  for (unsigned int i=0 ; i<myGenParticles->size(); i++) {
    if (TMath::Nint(   myGenParticles->at(i).status)==3) cout<<i<<"  "<<myGenParticles->at(i).firstMother<<"\t"<<myGenParticles->at(i).pdgId
							     <<"\t"<<myGenParticles->at(i).firstDaughter<<endl;
  }
*/

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


void EventCalculator::getSphericityJetMET(float & lambda1, float & lambda2, float & det,
			 const int jetmax, bool addMET) {

  TMatrixD top3Jets(2,2);
  double top3JetsScale = 0;
  
  unsigned int njets=  jets_AK5PF_pt->size();
  int ngoodj=0;
  for (unsigned int i=0; i<njets ; i++) {
    if ( isGoodJet(i) ) {
      ++ngoodj;
      double phi =jets_AK5PF_phi->at(i);
      double eT = getJetPt(i);
      double eX = eT*cos(phi);
      double eY = eT*sin(phi);
      if(ngoodj <= jetmax) {
	top3Jets[0][0] += eX*eX/eT;
	top3Jets[0][1] += eX*eY/eT;
	top3Jets[1][0] += eX*eY/eT;
	top3Jets[1][1] += eY*eY/eT;
	top3JetsScale += eT;
      }
      else break;
    }
  }

  if (addMET) {
    double phi = getMETphi();
    double eT = getMET();
    double eX = eT*cos(phi);
    double eY = eT*sin(phi);
  
    top3Jets[0][0] += eX*eX/eT;
    top3Jets[0][1] += eX*eY/eT;
    top3Jets[1][0] += eX*eY/eT;
    top3Jets[1][1] += eY*eY/eT;
    top3JetsScale += eT;
  }
  top3Jets*=1/top3JetsScale;
  TMatrixDEigen top3JetsEigen(top3Jets);
  lambda1 = top3JetsEigen.GetEigenValuesRe()[0];
  lambda2 = top3JetsEigen.GetEigenValuesRe()[1];
  det = top3Jets.Determinant();
  
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
  if (sampleName_.Contains("ttjets_madgraph") )
    inFile.open("badeventlist_pythia6bug_ttjets_madgraph.txt");
  //else if (sampleName_.Contains("t_s-channel") )
  else if (sampleName_.Contains("T_TuneZ2_s-channel_7TeV-powheg-tauola") )            
    inFile.open("badeventlist_pythia6bug_t_schannel_madgraph.txt");
  //else if (sampleName_.Contains("tbar_s-channel") )
  else if (sampleName_.Contains("Tbar_TuneZ2_s-channel_7TeV-powheg-tauola") )         
    inFile.open("badeventlist_pythia6bug_tbar_schannel_madgraph.txt");
  //else if (sampleName_.Contains("t_t-channel") )
  else if (sampleName_.Contains("T_TuneZ2_t-channel_7TeV-powheg-tauola") )            
    inFile.open("badeventlist_pythia6bug_t_tchannel_madgraph.txt");
  //else if (sampleName_.Contains("tbar_t-channel") )
  else if (sampleName_.Contains("Tbar_TuneZ2_t-channel_7TeV-powheg-tauola") )         
    inFile.open("badeventlist_pythia6bug_tbar_tchannel_madgraph.txt");
  //else if (sampleName_.Contains("t_tW-channel") )
  else if (sampleName_.Contains("T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola") )        
    inFile.open("badeventlist_pythia6bug_t_tWchannel_madgraph.txt");
  //else if (sampleName_.Contains("tbar_tW-channel") )
  else if (sampleName_.Contains("Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola") )     
    inFile.open("badeventlist_pythia6bug_tbar_tWchannel_madgraph.txt");
  else if (sampleName_.Contains("WJetsToLNu_300_HT_inf_TuneZ2_7TeV-madgraph-tauola") )
    inFile.open("badeventlist_pythia6bug_wjets_ht300_madgraph.txt");
  else if (sampleName_.Contains("DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola") )
    inFile.open("badeventlist_pythia6bug_dyjets_madgraph.txt");
  else if (sampleName_.Contains("zjets") )
    inFile.open("badeventlist_pythia6bug_zinvjets_madgraph.txt");
  else return;

  if(!inFile) {std::cout << "ERROR: can't open event list" << std::endl;  assert(0);}
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
  setCutScheme();
  setIgnoredCut("cut1b");
  setIgnoredCut("cut2b");
  setIgnoredCut("cut3b");

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

//separate out the jettag efficiency code from sampleAnalyzer for ease of use
void EventCalculator::plotBTagEffMC( ) {


  TFile fout("histos.root","RECREATE");
  
  //check the MC efficiencies from Summer11 TTbar
  //this histogram has bin edges {30,50,75,100,150,200,240,500,1000}
  //double bins[10] = {30,50,75,100,150,200,240,350,500,1000};
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

  setCutScheme();
  setIgnoredCut("cut1b");
  setIgnoredCut("cut2b");
  setIgnoredCut("cut3b");

  cout<<"Running..."<<endl;  
  int npass = 0;

  int ntaggedjets = 0;
  int ntaggedjets_b = 0;
  int ntaggedjets_c = 0;
  int ntaggedjets_l = 0;


  startTimer();
  for(Long64_t entry=0; entry < nevents; ++entry){
    chainB->GetEntry(entry);
    chainA->GetEntry(entry);

    if(entry%10000==0) cout << "entry: " << entry << ", percent done=" << (int)(entry/(double)nevents*100.)<<  endl;
    
    npass++;
    
    for (unsigned int i = 0; i < jets_AK5PF_pt->size(); ++i) {
      int flavor = jets_AK5PF_partonFlavour->at(i);
      
      if(isGoodJet(i,30)){
      
	  
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
	
      if (isGoodJet(i,30) && passBTagger(i) ){
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

cout << "initializing A" << endl;
   // Set object pointer
   trigger_prescalevalue = 0;
   trigger_name = 0;
   trigger_decision = 0;
   triggerobject_pt = 0;
   triggerobject_phi = 0;
   triggerobject_eta = 0;
   L1trigger_prescalevalue = 0;
   L1trigger_name = 0;
   L1trigger_decision = 0;
   els_conversion_dist = 0;
   els_conversion_dcot = 0;
   jets_AK5PFclean_corrL2L3 = 0;
   jets_AK5PFclean_corrL2L3Residual = 0;
   jets_AK5PFclean_corrL1FastL2L3 = 0;
   jets_AK5PFclean_corrL1L2L3 = 0;
   jets_AK5PFclean_corrL1FastL2L3Residual = 0;
   jets_AK5PFclean_corrL1L2L3Residual = 0;
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
   fChain->SetBranchAddress("triggerobject_pt", &triggerobject_pt, &b_triggerobject_pt);
   fChain->SetBranchAddress("triggerobject_phi", &triggerobject_phi, &b_triggerobject_phi);
   fChain->SetBranchAddress("triggerobject_eta", &triggerobject_eta, &b_triggerobject_eta);
   fChain->SetBranchAddress("L1trigger_prescalevalue", &L1trigger_prescalevalue, &b_L1trigger_prescalevalue);
   fChain->SetBranchAddress("L1trigger_name", &L1trigger_name, &b_L1trigger_name);
   fChain->SetBranchAddress("L1trigger_decision", &L1trigger_decision, &b_L1trigger_decision);
   fChain->SetBranchAddress("els_conversion_dist", &els_conversion_dist, &b_els_conversion_dist);
   fChain->SetBranchAddress("els_conversion_dcot", &els_conversion_dcot, &b_els_conversion_dcot);
   fChain->SetBranchAddress("hbhefilter_decision", &hbhefilter_decision, &b_hbhefilter_decision);
   fChain->SetBranchAddress("trackingfailurefilter_decision", &trackingfailurefilter_decision, &b_trackingfailurefilter_decision);
   fChain->SetBranchAddress("cschalofilter_decision", &cschalofilter_decision, &b_cschalofilter_decision);
   fChain->SetBranchAddress("ecalTPfilter_decision", &ecalTPfilter_decision, &b_ecalTPfilter_decision);
   fChain->SetBranchAddress("ecalBEfilter_decision", &ecalBEfilter_decision, &b_ecalBEfilter_decision);
   fChain->SetBranchAddress("scrapingVeto_decision", &scrapingVeto_decision, &b_scrapingVeto_decision);
   fChain->SetBranchAddress("greedymuonfilter_decision", &greedymuonfilter_decision, &b_greedymuonfilter_decision);
   fChain->SetBranchAddress("inconsistentPFmuonfilter_decision", &inconsistentPFmuonfilter_decision, &b_inconsistentPFmuonfilter_decision);
   fChain->SetBranchAddress("eenoisefilter_decision", &eenoisefilter_decision, &b_eenoisefilter_decision);
   fChain->SetBranchAddress("passprescalePFHT350filter_decision", &passprescalePFHT350filter_decision, &b_passprescalePFHT350filter_decision);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL2L3", &jets_AK5PFclean_corrL2L3, &b_jets_AK5PFclean_corrL2L3);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL2L3Residual", &jets_AK5PFclean_corrL2L3Residual, &b_jets_AK5PFclean_corrL2L3Residual);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL1FastL2L3", &jets_AK5PFclean_corrL1FastL2L3, &b_jets_AK5PFclean_corrL1FastL2L3);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL1L2L3", &jets_AK5PFclean_corrL1L2L3, &b_jets_AK5PFclean_corrL1L2L3);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL1FastL2L3Residual", &jets_AK5PFclean_corrL1FastL2L3Residual, &b_jets_AK5PFclean_corrL1FastL2L3Residual);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL1L2L3Residual", &jets_AK5PFclean_corrL1L2L3Residual, &b_jets_AK5PFclean_corrL1L2L3Residual);
   fChain->SetBranchAddress("PU_zpositions", &PU_zpositions, &b_PU_zpositions);
   fChain->SetBranchAddress("PU_sumpT_lowpT", &PU_sumpT_lowpT, &b_PU_sumpT_lowpT);
   fChain->SetBranchAddress("PU_sumpT_highpT", &PU_sumpT_highpT, &b_PU_sumpT_highpT);
   fChain->SetBranchAddress("PU_ntrks_lowpT", &PU_ntrks_lowpT, &b_PU_ntrks_lowpT);
   fChain->SetBranchAddress("PU_ntrks_highpT", &PU_ntrks_highpT, &b_PU_ntrks_highpT);
   fChain->SetBranchAddress("PU_NumInteractions", &PU_NumInteractions, &b_PU_NumInteractions);
   fChain->SetBranchAddress("PU_bunchCrossing", &PU_bunchCrossing, &b_PU_bunchCrossing);
   fChain->SetBranchAddress("PU_TrueNumInteractions", &PU_TrueNumInteractions, &b_PU_TrueNumInteractions);

   cout << "end initialize a" << endl;
}

void EventCalculator::InitializeB(TChain *fChain)
{
   // Set object pointer
   beamSpot_x = 0;
   beamSpot_y = 0;
   els_energy = 0;
   els_et = 0;
   els_eta = 0;
   els_phi = 0;
   els_pt = 0;
   els_px = 0;
   els_py = 0;
   els_pz = 0;
   els_robustTightId = 0;
   els_simpleEleId95relIso = 0;
   els_simpleEleId90relIso = 0;
   els_simpleEleId85relIso = 0;
   els_simpleEleId80relIso = 0;
   els_simpleEleId70relIso = 0;
   els_simpleEleId95cIso = 0;
   els_simpleEleId90cIso = 0;
   els_simpleEleId85cIso = 0;
   els_simpleEleId80cIso = 0;
   els_simpleEleId70cIso = 0;
   els_cIso = 0;
   els_tIso = 0;
   els_ecalIso = 0;
   els_hcalIso = 0;
   els_charge = 0;
   els_hadOverEm = 0;
   els_eOverPIn = 0;
   els_sigmaIEtaIEta = 0;
   els_scEnergy = 0;
   els_scEta = 0;
   els_scE1x5 = 0;
   els_scE2x5Max = 0;
   els_scE5x5 = 0;
   els_isEB = 0;
   els_isEE = 0;
   els_dEtaIn = 0;
   els_dPhiIn = 0;
   els_dEtaOut = 0;
   els_dPhiOut = 0;
   els_numlosthits = 0;
   els_tk_phi = 0;
   els_d0dum = 0;
   els_vx = 0;
   els_vy = 0;
   els_vz = 0;
   els_ptError = 0;
   els_n_inner_layer = 0;
   els_dr03EcalRecHitSumEt = 0;
   els_dr03HcalTowerSumEt = 0;
   els_dr03HcalDepth1TowerSumEt = 0;
   els_dr03HcalDepth2TowerSumEt = 0;
   els_dr03TkSumPt = 0;
   jets_AK5PF_phi = 0;
   jets_AK5PF_pt = 0;
   jets_AK5PF_pz = 0;
   jets_AK5PF_px = 0;
   jets_AK5PF_py = 0;
   jets_AK5PF_eta = 0;
   jets_AK5PF_et = 0;
   jets_AK5PF_energy = 0;
   jets_AK5PF_parton_Id = 0;
   jets_AK5PF_parton_motherId = 0;
   jets_AK5PF_gen_pt = 0;
   jets_AK5PF_gen_phi = 0;
   jets_AK5PF_partonFlavour = 0;
   jets_AK5PF_btag_TC_highPur = 0;
   jets_AK5PF_btag_TC_highEff = 0;
   jets_AK5PF_btag_jetProb = 0;
   jets_AK5PF_btag_jetBProb = 0;
   jets_AK5PF_btag_secVertexHighPur = 0;
   jets_AK5PF_btag_secVertexHighEff = 0;
   jets_AK5PF_btag_secVertexCombined = 0;
   jets_AK5PF_jetCharge = 0;
   jets_AK5PF_chgEmE = 0;
   jets_AK5PF_chgHadE = 0;
   jets_AK5PF_photonEnergy = 0;
   jets_AK5PF_chg_Mult = 0;
   jets_AK5PF_neutralEmE = 0;
   jets_AK5PF_neutralHadE = 0;
   jets_AK5PF_neutral_Mult = 0;
   jets_AK5PF_mu_Mult = 0;
   jets_AK5PF_ehf = 0;
   jets_AK5PF_corrFactorRaw = 0;
   jets_AK5PFclean_phi = 0;
   jets_AK5PFclean_pt = 0;
   jets_AK5PFclean_pz = 0;
   jets_AK5PFclean_px = 0;
   jets_AK5PFclean_py = 0;
   jets_AK5PFclean_eta = 0;
   jets_AK5PFclean_et = 0;
   jets_AK5PFclean_energy = 0;
   jets_AK5PFclean_parton_Id = 0;
   jets_AK5PFclean_parton_motherId = 0;
   jets_AK5PFclean_gen_pt = 0;
   jets_AK5PFclean_gen_phi = 0;
   jets_AK5PFclean_partonFlavour = 0;
   jets_AK5PFclean_btag_TC_highPur = 0;
   jets_AK5PFclean_btag_TC_highEff = 0;
   jets_AK5PFclean_btag_jetProb = 0;
   jets_AK5PFclean_btag_jetBProb = 0;
   jets_AK5PFclean_btag_secVertexHighPur = 0;
   jets_AK5PFclean_btag_secVertexHighEff = 0;
   jets_AK5PFclean_btag_secVertexCombined = 0;
   jets_AK5PFclean_jetCharge = 0;
   jets_AK5PFclean_chgEmE = 0;
   jets_AK5PFclean_chgHadE = 0;
   jets_AK5PFclean_photonEnergy = 0;
   jets_AK5PFclean_chg_Mult = 0;
   jets_AK5PFclean_neutralEmE = 0;
   jets_AK5PFclean_neutralHadE = 0;
   jets_AK5PFclean_neutral_Mult = 0;
   jets_AK5PFclean_mu_Mult = 0;
   jets_AK5PFclean_ehf = 0;
   jets_AK5PFclean_corrFactorRaw = 0;
   mc_doc_id = 0;
   mc_doc_pt = 0;
   mc_doc_px = 0;
   mc_doc_py = 0;
   mc_doc_pz = 0;
   mc_doc_eta = 0;
   mc_doc_phi = 0;
   mc_doc_mother_id = 0;
   mc_doc_grandmother_id = 0;
   mc_doc_ggrandmother_id = 0;
   mc_doc_mass = 0;
   mc_electrons_pt = 0;
   mc_electrons_px = 0;
   mc_electrons_py = 0;
   mc_electrons_pz = 0;
   mc_electrons_eta = 0;
   mc_electrons_phi = 0;
   mc_electrons_energy = 0;
   mc_electrons_mother_id = 0;
   mc_electrons_grandmother_id = 0;
   mc_mus_pt = 0;
   mc_mus_px = 0;
   mc_mus_py = 0;
   mc_mus_pz = 0;
   mc_mus_eta = 0;
   mc_mus_phi = 0;
   mc_mus_energy = 0;
   mc_mus_mother_id = 0;
   mc_mus_grandmother_id = 0;
   mc_nues_pt = 0;
   mc_nues_px = 0;
   mc_nues_py = 0;
   mc_nues_pz = 0;
   mc_nues_eta = 0;
   mc_nues_phi = 0;
   mc_nues_energy = 0;
   mc_nues_mother_id = 0;
   mc_nues_grandmother_id = 0;
   mc_numus_pt = 0;
   mc_numus_px = 0;
   mc_numus_py = 0;
   mc_numus_pz = 0;
   mc_numus_eta = 0;
   mc_numus_phi = 0;
   mc_numus_energy = 0;
   mc_numus_mother_id = 0;
   mc_numus_grandmother_id = 0;
   mc_nutaus_pt = 0;
   mc_nutaus_px = 0;
   mc_nutaus_py = 0;
   mc_nutaus_pz = 0;
   mc_nutaus_eta = 0;
   mc_nutaus_phi = 0;
   mc_nutaus_energy = 0;
   mc_nutaus_mother_id = 0;
   mc_nutaus_grandmother_id = 0;
   mets_AK5_et = 0;
   mets_AK5_phi = 0;
   mus_energy = 0;
   mus_et = 0;
   mus_eta = 0;
   mus_phi = 0;
   mus_pt = 0;
   mus_px = 0;
   mus_py = 0;
   mus_pz = 0;
   mus_numberOfMatchedStations = 0;
   mus_cIso = 0;
   mus_tIso = 0;
   mus_ecalIso = 0;
   mus_hcalIso = 0;
   mus_ecalvetoDep = 0;
   mus_hcalvetoDep = 0;
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
   mus_cm_pt = 0;
   mus_cm_ptErr = 0;
   mus_tk_chi2 = 0;
   mus_tk_ndof = 0;
   mus_tk_pt = 0;
   mus_tk_phi = 0;
   mus_tk_d0dum = 0;
   mus_tk_vz = 0;
   mus_tk_numvalhits = 0;
   mus_tk_ptErr = 0;
   mus_tk_numvalPixelhits = 0;
   mus_tk_numpixelWthMeasr = 0;
   mus_stamu_pt = 0;
   mus_stamu_ptErr = 0;
   mus_num_matches = 0;
   mus_isPFMuon = 0;
   mus_isTrackerMuon = 0;
   mus_isGlobalMuon = 0;
   mus_id_AllGlobalMuons = 0;
   mus_id_AllTrackerMuons = 0;
   mus_id_GlobalMuonPromptTight = 0;
   pfTypeImets_et = 0;
   pfTypeImets_phi = 0;
   pfTypeImets_ex = 0;
   pfTypeImets_ey = 0;
   pfTypeImets_gen_et = 0;
   pfTypeImets_gen_phi = 0;
   pfTypeImets_sumEt = 0;
   pf_els_energy = 0;
   pf_els_et = 0;
   pf_els_eta = 0;
   pf_els_phi = 0;
   pf_els_pt = 0;
   pf_els_px = 0;
   pf_els_py = 0;
   pf_els_pz = 0;
   pf_els_robustTightId = 0;
   pf_els_simpleEleId95relIso = 0;
   pf_els_simpleEleId90relIso = 0;
   pf_els_simpleEleId85relIso = 0;
   pf_els_simpleEleId80relIso = 0;
   pf_els_simpleEleId70relIso = 0;
   pf_els_simpleEleId95cIso = 0;
   pf_els_simpleEleId90cIso = 0;
   pf_els_simpleEleId85cIso = 0;
   pf_els_simpleEleId80cIso = 0;
   pf_els_simpleEleId70cIso = 0;
   pf_els_cIso = 0;
   pf_els_tIso = 0;
   pf_els_ecalIso = 0;
   pf_els_hcalIso = 0;
   pf_els_chargedHadronIso = 0;
   pf_els_photonIso = 0;
   pf_els_neutralHadronIso = 0;
   pf_els_charge = 0;
   pf_els_hadOverEm = 0;
   pf_els_eOverPIn = 0;
   pf_els_sigmaIEtaIEta = 0;
   pf_els_scEnergy = 0;
   pf_els_scEta = 0;
   pf_els_scE1x5 = 0;
   pf_els_scE2x5Max = 0;
   pf_els_scE5x5 = 0;
   pf_els_isEB = 0;
   pf_els_isEE = 0;
   pf_els_dEtaIn = 0;
   pf_els_dPhiIn = 0;
   pf_els_dEtaOut = 0;
   pf_els_dPhiOut = 0;
   pf_els_numlosthits = 0;
   pf_els_tk_phi = 0;
   pf_els_d0dum = 0;
   pf_els_vx = 0;
   pf_els_vy = 0;
   pf_els_vz = 0;
   pf_els_ptError = 0;
   pf_els_n_inner_layer = 0;
   pf_els_dr03EcalRecHitSumEt = 0;
   pf_els_dr03HcalTowerSumEt = 0;
   pf_els_dr03HcalDepth1TowerSumEt = 0;
   pf_els_dr03HcalDepth2TowerSumEt = 0;
   pf_els_dr03TkSumPt = 0;
   pf_mus_energy = 0;
   pf_mus_et = 0;
   pf_mus_eta = 0;
   pf_mus_phi = 0;
   pf_mus_pt = 0;
   pf_mus_px = 0;
   pf_mus_py = 0;
   pf_mus_pz = 0;
   pf_mus_cIso = 0;
   pf_mus_tIso = 0;
   pf_mus_ecalIso = 0;
   pf_mus_hcalIso = 0;
   pf_mus_neutralHadronIso = 0;
   pf_mus_chargedHadronIso = 0;
   pf_mus_photonIso = 0;
   pf_mus_charge = 0;
   pf_mus_cm_chi2 = 0;
   pf_mus_cm_ndof = 0;
   pf_mus_cm_pt = 0;
   pf_mus_cm_ptErr = 0;
   pf_mus_tk_chi2 = 0;
   pf_mus_tk_ndof = 0;
   pf_mus_tk_pt = 0;
   pf_mus_tk_phi = 0;
   pf_mus_tk_d0dum = 0;
   pf_mus_tk_vz = 0;
   pf_mus_tk_numvalhits = 0;
   pf_mus_tk_ptErr = 0;
   pf_mus_tk_numvalPixelhits = 0;
   pf_mus_tk_numpixelWthMeasr = 0;
   pf_mus_stamu_pt = 0;
   pf_mus_stamu_ptErr = 0;
   pf_mus_num_matches = 0;
   pf_mus_isTrackerMuon = 0;
   pf_mus_isGlobalMuon = 0;
   pf_mus_id_GlobalMuonPromptTight = 0;
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
   pfmets_sumEt = 0;
   pv_x = 0;
   pv_y = 0;
   pv_z = 0;
   pv_ndof = 0;
   pv_isFake = 0;
   pv_tracksSize = 0;
   model_params = 0;
   fChain->SetBranchAddress("beamSpot_x", &beamSpot_x, &b_beamSpot_x);
   fChain->SetBranchAddress("beamSpot_y", &beamSpot_y, &b_beamSpot_y);
   fChain->SetBranchAddress("els_energy", &els_energy, &b_els_energy);
   fChain->SetBranchAddress("els_et", &els_et, &b_els_et);
   fChain->SetBranchAddress("els_eta", &els_eta, &b_els_eta);
   fChain->SetBranchAddress("els_phi", &els_phi, &b_els_phi);
   fChain->SetBranchAddress("els_pt", &els_pt, &b_els_pt);
   fChain->SetBranchAddress("els_px", &els_px, &b_els_px);
   fChain->SetBranchAddress("els_py", &els_py, &b_els_py);
   fChain->SetBranchAddress("els_pz", &els_pz, &b_els_pz);
   fChain->SetBranchAddress("els_robustTightId", &els_robustTightId, &b_els_robustTightId);
   fChain->SetBranchAddress("els_simpleEleId95relIso", &els_simpleEleId95relIso, &b_els_simpleEleId95relIso);
   fChain->SetBranchAddress("els_simpleEleId90relIso", &els_simpleEleId90relIso, &b_els_simpleEleId90relIso);
   fChain->SetBranchAddress("els_simpleEleId85relIso", &els_simpleEleId85relIso, &b_els_simpleEleId85relIso);
   fChain->SetBranchAddress("els_simpleEleId80relIso", &els_simpleEleId80relIso, &b_els_simpleEleId80relIso);
   fChain->SetBranchAddress("els_simpleEleId70relIso", &els_simpleEleId70relIso, &b_els_simpleEleId70relIso);
   fChain->SetBranchAddress("els_simpleEleId95cIso", &els_simpleEleId95cIso, &b_els_simpleEleId95cIso);
   fChain->SetBranchAddress("els_simpleEleId90cIso", &els_simpleEleId90cIso, &b_els_simpleEleId90cIso);
   fChain->SetBranchAddress("els_simpleEleId85cIso", &els_simpleEleId85cIso, &b_els_simpleEleId85cIso);
   fChain->SetBranchAddress("els_simpleEleId80cIso", &els_simpleEleId80cIso, &b_els_simpleEleId80cIso);
   fChain->SetBranchAddress("els_simpleEleId70cIso", &els_simpleEleId70cIso, &b_els_simpleEleId70cIso);
   fChain->SetBranchAddress("els_cIso", &els_cIso, &b_els_cIso);
   fChain->SetBranchAddress("els_tIso", &els_tIso, &b_els_tIso);
   fChain->SetBranchAddress("els_ecalIso", &els_ecalIso, &b_els_ecalIso);
   fChain->SetBranchAddress("els_hcalIso", &els_hcalIso, &b_els_hcalIso);
   fChain->SetBranchAddress("els_charge", &els_charge, &b_els_charge);
   fChain->SetBranchAddress("els_hadOverEm", &els_hadOverEm, &b_els_hadOverEm);
   fChain->SetBranchAddress("els_eOverPIn", &els_eOverPIn, &b_els_eOverPIn);
   fChain->SetBranchAddress("els_sigmaIEtaIEta", &els_sigmaIEtaIEta, &b_els_sigmaIEtaIEta);
   fChain->SetBranchAddress("els_scEnergy", &els_scEnergy, &b_els_scEnergy);
   fChain->SetBranchAddress("els_scEta", &els_scEta, &b_els_scEta);
   fChain->SetBranchAddress("els_scE1x5", &els_scE1x5, &b_els_scE1x5);
   fChain->SetBranchAddress("els_scE2x5Max", &els_scE2x5Max, &b_els_scE2x5Max);
   fChain->SetBranchAddress("els_scE5x5", &els_scE5x5, &b_els_scE5x5);
   fChain->SetBranchAddress("els_isEB", &els_isEB, &b_els_isEB);
   fChain->SetBranchAddress("els_isEE", &els_isEE, &b_els_isEE);
   fChain->SetBranchAddress("els_dEtaIn", &els_dEtaIn, &b_els_dEtaIn);
   fChain->SetBranchAddress("els_dPhiIn", &els_dPhiIn, &b_els_dPhiIn);
   fChain->SetBranchAddress("els_dEtaOut", &els_dEtaOut, &b_els_dEtaOut);
   fChain->SetBranchAddress("els_dPhiOut", &els_dPhiOut, &b_els_dPhiOut);
   fChain->SetBranchAddress("els_numlosthits", &els_numlosthits, &b_els_numlosthits);
   fChain->SetBranchAddress("els_tk_phi", &els_tk_phi, &b_els_tk_phi);
   fChain->SetBranchAddress("els_d0dum", &els_d0dum, &b_els_d0dum);
   fChain->SetBranchAddress("els_vx", &els_vx, &b_els_vx);
   fChain->SetBranchAddress("els_vy", &els_vy, &b_els_vy);
   fChain->SetBranchAddress("els_vz", &els_vz, &b_els_vz);
   fChain->SetBranchAddress("els_ptError", &els_ptError, &b_els_ptError);
   fChain->SetBranchAddress("els_n_inner_layer", &els_n_inner_layer, &b_els_n_inner_layer);
   fChain->SetBranchAddress("els_dr03EcalRecHitSumEt", &els_dr03EcalRecHitSumEt, &b_els_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("els_dr03HcalTowerSumEt", &els_dr03HcalTowerSumEt, &b_els_dr03HcalTowerSumEt);
   fChain->SetBranchAddress("els_dr03HcalDepth1TowerSumEt", &els_dr03HcalDepth1TowerSumEt, &b_els_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("els_dr03HcalDepth2TowerSumEt", &els_dr03HcalDepth2TowerSumEt, &b_els_dr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("els_dr03TkSumPt", &els_dr03TkSumPt, &b_els_dr03TkSumPt);
   fChain->SetBranchAddress("jets_AK5PF_phi", &jets_AK5PF_phi, &b_jets_AK5PF_phi);
   fChain->SetBranchAddress("jets_AK5PF_pt", &jets_AK5PF_pt, &b_jets_AK5PF_pt);
   fChain->SetBranchAddress("jets_AK5PF_pz", &jets_AK5PF_pz, &b_jets_AK5PF_pz);
   fChain->SetBranchAddress("jets_AK5PF_px", &jets_AK5PF_px, &b_jets_AK5PF_px);
   fChain->SetBranchAddress("jets_AK5PF_py", &jets_AK5PF_py, &b_jets_AK5PF_py);
   fChain->SetBranchAddress("jets_AK5PF_eta", &jets_AK5PF_eta, &b_jets_AK5PF_eta);
   fChain->SetBranchAddress("jets_AK5PF_et", &jets_AK5PF_et, &b_jets_AK5PF_et);
   fChain->SetBranchAddress("jets_AK5PF_energy", &jets_AK5PF_energy, &b_jets_AK5PF_energy);
   fChain->SetBranchAddress("jets_AK5PF_parton_Id", &jets_AK5PF_parton_Id, &b_jets_AK5PF_parton_Id);
   fChain->SetBranchAddress("jets_AK5PF_parton_motherId", &jets_AK5PF_parton_motherId, &b_jets_AK5PF_parton_motherId);
   fChain->SetBranchAddress("jets_AK5PF_gen_pt", &jets_AK5PF_gen_pt, &b_jets_AK5PF_gen_pt);
   fChain->SetBranchAddress("jets_AK5PF_gen_phi", &jets_AK5PF_gen_phi, &b_jets_AK5PF_gen_phi);
   fChain->SetBranchAddress("jets_AK5PF_partonFlavour", &jets_AK5PF_partonFlavour, &b_jets_AK5PF_partonFlavour);
   fChain->SetBranchAddress("jets_AK5PF_btag_TC_highPur", &jets_AK5PF_btag_TC_highPur, &b_jets_AK5PF_btag_TC_highPur);
   fChain->SetBranchAddress("jets_AK5PF_btag_TC_highEff", &jets_AK5PF_btag_TC_highEff, &b_jets_AK5PF_btag_TC_highEff);
   fChain->SetBranchAddress("jets_AK5PF_btag_jetProb", &jets_AK5PF_btag_jetProb, &b_jets_AK5PF_btag_jetProb);
   fChain->SetBranchAddress("jets_AK5PF_btag_jetBProb", &jets_AK5PF_btag_jetBProb, &b_jets_AK5PF_btag_jetBProb);
   fChain->SetBranchAddress("jets_AK5PF_btag_secVertexHighPur", &jets_AK5PF_btag_secVertexHighPur, &b_jets_AK5PF_btag_secVertexHighPur);
   fChain->SetBranchAddress("jets_AK5PF_btag_secVertexHighEff", &jets_AK5PF_btag_secVertexHighEff, &b_jets_AK5PF_btag_secVertexHighEff);
   fChain->SetBranchAddress("jets_AK5PF_btag_secVertexCombined", &jets_AK5PF_btag_secVertexCombined, &b_jets_AK5PF_btag_secVertexCombined);
   fChain->SetBranchAddress("jets_AK5PF_jetCharge", &jets_AK5PF_jetCharge, &b_jets_AK5PF_jetCharge);
   fChain->SetBranchAddress("jets_AK5PF_chgEmE", &jets_AK5PF_chgEmE, &b_jets_AK5PF_chgEmE);
   fChain->SetBranchAddress("jets_AK5PF_chgHadE", &jets_AK5PF_chgHadE, &b_jets_AK5PF_chgHadE);
   fChain->SetBranchAddress("jets_AK5PF_photonEnergy", &jets_AK5PF_photonEnergy, &b_jets_AK5PF_photonEnergy);
   fChain->SetBranchAddress("jets_AK5PF_chg_Mult", &jets_AK5PF_chg_Mult, &b_jets_AK5PF_chg_Mult);
   fChain->SetBranchAddress("jets_AK5PF_neutralEmE", &jets_AK5PF_neutralEmE, &b_jets_AK5PF_neutralEmE);
   fChain->SetBranchAddress("jets_AK5PF_neutralHadE", &jets_AK5PF_neutralHadE, &b_jets_AK5PF_neutralHadE);
   fChain->SetBranchAddress("jets_AK5PF_neutral_Mult", &jets_AK5PF_neutral_Mult, &b_jets_AK5PF_neutral_Mult);
   fChain->SetBranchAddress("jets_AK5PF_mu_Mult", &jets_AK5PF_mu_Mult, &b_jets_AK5PF_mu_Mult);
   fChain->SetBranchAddress("jets_AK5PF_ehf", &jets_AK5PF_ehf, &b_jets_AK5PF_ehf);
   fChain->SetBranchAddress("jets_AK5PF_corrFactorRaw", &jets_AK5PF_corrFactorRaw, &b_jets_AK5PF_corrFactorRaw);
   fChain->SetBranchAddress("jets_AK5PFclean_phi", &jets_AK5PFclean_phi, &b_jets_AK5PFclean_phi);
   fChain->SetBranchAddress("jets_AK5PFclean_pt", &jets_AK5PFclean_pt, &b_jets_AK5PFclean_pt);
   fChain->SetBranchAddress("jets_AK5PFclean_pz", &jets_AK5PFclean_pz, &b_jets_AK5PFclean_pz);
   fChain->SetBranchAddress("jets_AK5PFclean_px", &jets_AK5PFclean_px, &b_jets_AK5PFclean_px);
   fChain->SetBranchAddress("jets_AK5PFclean_py", &jets_AK5PFclean_py, &b_jets_AK5PFclean_py);
   fChain->SetBranchAddress("jets_AK5PFclean_eta", &jets_AK5PFclean_eta, &b_jets_AK5PFclean_eta);
   fChain->SetBranchAddress("jets_AK5PFclean_et", &jets_AK5PFclean_et, &b_jets_AK5PFclean_et);
   fChain->SetBranchAddress("jets_AK5PFclean_energy", &jets_AK5PFclean_energy, &b_jets_AK5PFclean_energy);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_Id", &jets_AK5PFclean_parton_Id, &b_jets_AK5PFclean_parton_Id);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_motherId", &jets_AK5PFclean_parton_motherId, &b_jets_AK5PFclean_parton_motherId);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_pt", &jets_AK5PFclean_gen_pt, &b_jets_AK5PFclean_gen_pt);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_phi", &jets_AK5PFclean_gen_phi, &b_jets_AK5PFclean_gen_phi);
   fChain->SetBranchAddress("jets_AK5PFclean_partonFlavour", &jets_AK5PFclean_partonFlavour, &b_jets_AK5PFclean_partonFlavour);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_TC_highPur", &jets_AK5PFclean_btag_TC_highPur, &b_jets_AK5PFclean_btag_TC_highPur);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_TC_highEff", &jets_AK5PFclean_btag_TC_highEff, &b_jets_AK5PFclean_btag_TC_highEff);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_jetProb", &jets_AK5PFclean_btag_jetProb, &b_jets_AK5PFclean_btag_jetProb);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_jetBProb", &jets_AK5PFclean_btag_jetBProb, &b_jets_AK5PFclean_btag_jetBProb);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_secVertexHighPur", &jets_AK5PFclean_btag_secVertexHighPur, &b_jets_AK5PFclean_btag_secVertexHighPur);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_secVertexHighEff", &jets_AK5PFclean_btag_secVertexHighEff, &b_jets_AK5PFclean_btag_secVertexHighEff);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_secVertexCombined", &jets_AK5PFclean_btag_secVertexCombined, &b_jets_AK5PFclean_btag_secVertexCombined);
   fChain->SetBranchAddress("jets_AK5PFclean_jetCharge", &jets_AK5PFclean_jetCharge, &b_jets_AK5PFclean_jetCharge);
   fChain->SetBranchAddress("jets_AK5PFclean_chgEmE", &jets_AK5PFclean_chgEmE, &b_jets_AK5PFclean_chgEmE);
   fChain->SetBranchAddress("jets_AK5PFclean_chgHadE", &jets_AK5PFclean_chgHadE, &b_jets_AK5PFclean_chgHadE);
   fChain->SetBranchAddress("jets_AK5PFclean_photonEnergy", &jets_AK5PFclean_photonEnergy, &b_jets_AK5PFclean_photonEnergy);
   fChain->SetBranchAddress("jets_AK5PFclean_chg_Mult", &jets_AK5PFclean_chg_Mult, &b_jets_AK5PFclean_chg_Mult);
   fChain->SetBranchAddress("jets_AK5PFclean_neutralEmE", &jets_AK5PFclean_neutralEmE, &b_jets_AK5PFclean_neutralEmE);
   fChain->SetBranchAddress("jets_AK5PFclean_neutralHadE", &jets_AK5PFclean_neutralHadE, &b_jets_AK5PFclean_neutralHadE);
   fChain->SetBranchAddress("jets_AK5PFclean_neutral_Mult", &jets_AK5PFclean_neutral_Mult, &b_jets_AK5PFclean_neutral_Mult);
   fChain->SetBranchAddress("jets_AK5PFclean_mu_Mult", &jets_AK5PFclean_mu_Mult, &b_jets_AK5PFclean_mu_Mult);
   fChain->SetBranchAddress("jets_AK5PFclean_ehf", &jets_AK5PFclean_ehf, &b_jets_AK5PFclean_ehf);
   fChain->SetBranchAddress("jets_AK5PFclean_corrFactorRaw", &jets_AK5PFclean_corrFactorRaw, &b_jets_AK5PFclean_corrFactorRaw);
   fChain->SetBranchAddress("mc_doc_id", &mc_doc_id, &b_mc_doc_id);
   fChain->SetBranchAddress("mc_doc_pt", &mc_doc_pt, &b_mc_doc_pt);
   fChain->SetBranchAddress("mc_doc_px", &mc_doc_px, &b_mc_doc_px);
   fChain->SetBranchAddress("mc_doc_py", &mc_doc_py, &b_mc_doc_py);
   fChain->SetBranchAddress("mc_doc_pz", &mc_doc_pz, &b_mc_doc_pz);
   fChain->SetBranchAddress("mc_doc_eta", &mc_doc_eta, &b_mc_doc_eta);
   fChain->SetBranchAddress("mc_doc_phi", &mc_doc_phi, &b_mc_doc_phi);
   fChain->SetBranchAddress("mc_doc_mother_id", &mc_doc_mother_id, &b_mc_doc_mother_id);
   fChain->SetBranchAddress("mc_doc_grandmother_id", &mc_doc_grandmother_id, &b_mc_doc_grandmother_id);
   fChain->SetBranchAddress("mc_doc_ggrandmother_id", &mc_doc_ggrandmother_id, &b_mc_doc_ggrandmother_id);
   fChain->SetBranchAddress("mc_doc_mass", &mc_doc_mass, &b_mc_doc_mass);
   fChain->SetBranchAddress("mc_electrons_pt", &mc_electrons_pt, &b_mc_electrons_pt);
   fChain->SetBranchAddress("mc_electrons_px", &mc_electrons_px, &b_mc_electrons_px);
   fChain->SetBranchAddress("mc_electrons_py", &mc_electrons_py, &b_mc_electrons_py);
   fChain->SetBranchAddress("mc_electrons_pz", &mc_electrons_pz, &b_mc_electrons_pz);
   fChain->SetBranchAddress("mc_electrons_eta", &mc_electrons_eta, &b_mc_electrons_eta);
   fChain->SetBranchAddress("mc_electrons_phi", &mc_electrons_phi, &b_mc_electrons_phi);
   fChain->SetBranchAddress("mc_electrons_energy", &mc_electrons_energy, &b_mc_electrons_energy);
   fChain->SetBranchAddress("mc_electrons_mother_id", &mc_electrons_mother_id, &b_mc_electrons_mother_id);
   fChain->SetBranchAddress("mc_electrons_grandmother_id", &mc_electrons_grandmother_id, &b_mc_electrons_grandmother_id);
   fChain->SetBranchAddress("mc_mus_pt", &mc_mus_pt, &b_mc_mus_pt);
   fChain->SetBranchAddress("mc_mus_px", &mc_mus_px, &b_mc_mus_px);
   fChain->SetBranchAddress("mc_mus_py", &mc_mus_py, &b_mc_mus_py);
   fChain->SetBranchAddress("mc_mus_pz", &mc_mus_pz, &b_mc_mus_pz);
   fChain->SetBranchAddress("mc_mus_eta", &mc_mus_eta, &b_mc_mus_eta);
   fChain->SetBranchAddress("mc_mus_phi", &mc_mus_phi, &b_mc_mus_phi);
   fChain->SetBranchAddress("mc_mus_energy", &mc_mus_energy, &b_mc_mus_energy);
   fChain->SetBranchAddress("mc_mus_mother_id", &mc_mus_mother_id, &b_mc_mus_mother_id);
   fChain->SetBranchAddress("mc_mus_grandmother_id", &mc_mus_grandmother_id, &b_mc_mus_grandmother_id);
   fChain->SetBranchAddress("mc_nues_pt", &mc_nues_pt, &b_mc_nues_pt);
   fChain->SetBranchAddress("mc_nues_px", &mc_nues_px, &b_mc_nues_px);
   fChain->SetBranchAddress("mc_nues_py", &mc_nues_py, &b_mc_nues_py);
   fChain->SetBranchAddress("mc_nues_pz", &mc_nues_pz, &b_mc_nues_pz);
   fChain->SetBranchAddress("mc_nues_eta", &mc_nues_eta, &b_mc_nues_eta);
   fChain->SetBranchAddress("mc_nues_phi", &mc_nues_phi, &b_mc_nues_phi);
   fChain->SetBranchAddress("mc_nues_energy", &mc_nues_energy, &b_mc_nues_energy);
   fChain->SetBranchAddress("mc_nues_mother_id", &mc_nues_mother_id, &b_mc_nues_mother_id);
   fChain->SetBranchAddress("mc_nues_grandmother_id", &mc_nues_grandmother_id, &b_mc_nues_grandmother_id);
   fChain->SetBranchAddress("mc_numus_pt", &mc_numus_pt, &b_mc_numus_pt);
   fChain->SetBranchAddress("mc_numus_px", &mc_numus_px, &b_mc_numus_px);
   fChain->SetBranchAddress("mc_numus_py", &mc_numus_py, &b_mc_numus_py);
   fChain->SetBranchAddress("mc_numus_pz", &mc_numus_pz, &b_mc_numus_pz);
   fChain->SetBranchAddress("mc_numus_eta", &mc_numus_eta, &b_mc_numus_eta);
   fChain->SetBranchAddress("mc_numus_phi", &mc_numus_phi, &b_mc_numus_phi);
   fChain->SetBranchAddress("mc_numus_energy", &mc_numus_energy, &b_mc_numus_energy);
   fChain->SetBranchAddress("mc_numus_mother_id", &mc_numus_mother_id, &b_mc_numus_mother_id);
   fChain->SetBranchAddress("mc_numus_grandmother_id", &mc_numus_grandmother_id, &b_mc_numus_grandmother_id);
   fChain->SetBranchAddress("mc_nutaus_pt", &mc_nutaus_pt, &b_mc_nutaus_pt);
   fChain->SetBranchAddress("mc_nutaus_px", &mc_nutaus_px, &b_mc_nutaus_px);
   fChain->SetBranchAddress("mc_nutaus_py", &mc_nutaus_py, &b_mc_nutaus_py);
   fChain->SetBranchAddress("mc_nutaus_pz", &mc_nutaus_pz, &b_mc_nutaus_pz);
   fChain->SetBranchAddress("mc_nutaus_eta", &mc_nutaus_eta, &b_mc_nutaus_eta);
   fChain->SetBranchAddress("mc_nutaus_phi", &mc_nutaus_phi, &b_mc_nutaus_phi);
   fChain->SetBranchAddress("mc_nutaus_energy", &mc_nutaus_energy, &b_mc_nutaus_energy);
   fChain->SetBranchAddress("mc_nutaus_mother_id", &mc_nutaus_mother_id, &b_mc_nutaus_mother_id);
   fChain->SetBranchAddress("mc_nutaus_grandmother_id", &mc_nutaus_grandmother_id, &b_mc_nutaus_grandmother_id);
   fChain->SetBranchAddress("mets_AK5_et", &mets_AK5_et, &b_mets_AK5_et);
   fChain->SetBranchAddress("mets_AK5_phi", &mets_AK5_phi, &b_mets_AK5_phi);
   fChain->SetBranchAddress("mus_energy", &mus_energy, &b_mus_energy);
   fChain->SetBranchAddress("mus_et", &mus_et, &b_mus_et);
   fChain->SetBranchAddress("mus_eta", &mus_eta, &b_mus_eta);
   fChain->SetBranchAddress("mus_phi", &mus_phi, &b_mus_phi);
   fChain->SetBranchAddress("mus_pt", &mus_pt, &b_mus_pt);
   fChain->SetBranchAddress("mus_px", &mus_px, &b_mus_px);
   fChain->SetBranchAddress("mus_py", &mus_py, &b_mus_py);
   fChain->SetBranchAddress("mus_pz", &mus_pz, &b_mus_pz);
   fChain->SetBranchAddress("mus_numberOfMatchedStations", &mus_numberOfMatchedStations, &b_mus_numberOfMatchedStations);
   fChain->SetBranchAddress("mus_cIso", &mus_cIso, &b_mus_cIso);
   fChain->SetBranchAddress("mus_tIso", &mus_tIso, &b_mus_tIso);
   fChain->SetBranchAddress("mus_ecalIso", &mus_ecalIso, &b_mus_ecalIso);
   fChain->SetBranchAddress("mus_hcalIso", &mus_hcalIso, &b_mus_hcalIso);
   fChain->SetBranchAddress("mus_ecalvetoDep", &mus_ecalvetoDep, &b_mus_ecalvetoDep);
   fChain->SetBranchAddress("mus_hcalvetoDep", &mus_hcalvetoDep, &b_mus_hcalvetoDep);
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
   fChain->SetBranchAddress("mus_cm_pt", &mus_cm_pt, &b_mus_cm_pt);
   fChain->SetBranchAddress("mus_cm_ptErr", &mus_cm_ptErr, &b_mus_cm_ptErr);
   fChain->SetBranchAddress("mus_tk_chi2", &mus_tk_chi2, &b_mus_tk_chi2);
   fChain->SetBranchAddress("mus_tk_ndof", &mus_tk_ndof, &b_mus_tk_ndof);
   fChain->SetBranchAddress("mus_tk_pt", &mus_tk_pt, &b_mus_tk_pt);
   fChain->SetBranchAddress("mus_tk_phi", &mus_tk_phi, &b_mus_tk_phi);
   fChain->SetBranchAddress("mus_tk_d0dum", &mus_tk_d0dum, &b_mus_tk_d0dum);
   fChain->SetBranchAddress("mus_tk_vz", &mus_tk_vz, &b_mus_tk_vz);
   fChain->SetBranchAddress("mus_tk_numvalhits", &mus_tk_numvalhits, &b_mus_tk_numvalhits);
   fChain->SetBranchAddress("mus_tk_ptErr", &mus_tk_ptErr, &b_mus_tk_ptErr);
   fChain->SetBranchAddress("mus_tk_numvalPixelhits", &mus_tk_numvalPixelhits, &b_mus_tk_numvalPixelhits);
   fChain->SetBranchAddress("mus_tk_numpixelWthMeasr", &mus_tk_numpixelWthMeasr, &b_mus_tk_numpixelWthMeasr);
   fChain->SetBranchAddress("mus_stamu_pt", &mus_stamu_pt, &b_mus_stamu_pt);
   fChain->SetBranchAddress("mus_stamu_ptErr", &mus_stamu_ptErr, &b_mus_stamu_ptErr);
   fChain->SetBranchAddress("mus_num_matches", &mus_num_matches, &b_mus_num_matches);
   fChain->SetBranchAddress("mus_isPFMuon", &mus_isPFMuon, &b_mus_isPFMuon);
   fChain->SetBranchAddress("mus_isTrackerMuon", &mus_isTrackerMuon, &b_mus_isTrackerMuon);
   fChain->SetBranchAddress("mus_isGlobalMuon", &mus_isGlobalMuon, &b_mus_isGlobalMuon);
   fChain->SetBranchAddress("mus_id_AllGlobalMuons", &mus_id_AllGlobalMuons, &b_mus_id_AllGlobalMuons);
   fChain->SetBranchAddress("mus_id_AllTrackerMuons", &mus_id_AllTrackerMuons, &b_mus_id_AllTrackerMuons);
   fChain->SetBranchAddress("mus_id_GlobalMuonPromptTight", &mus_id_GlobalMuonPromptTight, &b_mus_id_GlobalMuonPromptTight);
   fChain->SetBranchAddress("pfTypeImets_et", &pfTypeImets_et, &b_pfTypeImets_et);
   fChain->SetBranchAddress("pfTypeImets_phi", &pfTypeImets_phi, &b_pfTypeImets_phi);
   fChain->SetBranchAddress("pfTypeImets_ex", &pfTypeImets_ex, &b_pfTypeImets_ex);
   fChain->SetBranchAddress("pfTypeImets_ey", &pfTypeImets_ey, &b_pfTypeImets_ey);
   fChain->SetBranchAddress("pfTypeImets_gen_et", &pfTypeImets_gen_et, &b_pfTypeImets_gen_et);
   fChain->SetBranchAddress("pfTypeImets_gen_phi", &pfTypeImets_gen_phi, &b_pfTypeImets_gen_phi);
   fChain->SetBranchAddress("pfTypeImets_sumEt", &pfTypeImets_sumEt, &b_pfTypeImets_sumEt);
   fChain->SetBranchAddress("pf_els_energy", &pf_els_energy, &b_pf_els_energy);
   fChain->SetBranchAddress("pf_els_et", &pf_els_et, &b_pf_els_et);
   fChain->SetBranchAddress("pf_els_eta", &pf_els_eta, &b_pf_els_eta);
   fChain->SetBranchAddress("pf_els_phi", &pf_els_phi, &b_pf_els_phi);
   fChain->SetBranchAddress("pf_els_pt", &pf_els_pt, &b_pf_els_pt);
   fChain->SetBranchAddress("pf_els_px", &pf_els_px, &b_pf_els_px);
   fChain->SetBranchAddress("pf_els_py", &pf_els_py, &b_pf_els_py);
   fChain->SetBranchAddress("pf_els_pz", &pf_els_pz, &b_pf_els_pz);
   fChain->SetBranchAddress("pf_els_robustTightId", &pf_els_robustTightId, &b_pf_els_robustTightId);
   fChain->SetBranchAddress("pf_els_simpleEleId95relIso", &pf_els_simpleEleId95relIso, &b_pf_els_simpleEleId95relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId90relIso", &pf_els_simpleEleId90relIso, &b_pf_els_simpleEleId90relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId85relIso", &pf_els_simpleEleId85relIso, &b_pf_els_simpleEleId85relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId80relIso", &pf_els_simpleEleId80relIso, &b_pf_els_simpleEleId80relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId70relIso", &pf_els_simpleEleId70relIso, &b_pf_els_simpleEleId70relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId95cIso", &pf_els_simpleEleId95cIso, &b_pf_els_simpleEleId95cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId90cIso", &pf_els_simpleEleId90cIso, &b_pf_els_simpleEleId90cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId85cIso", &pf_els_simpleEleId85cIso, &b_pf_els_simpleEleId85cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId80cIso", &pf_els_simpleEleId80cIso, &b_pf_els_simpleEleId80cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId70cIso", &pf_els_simpleEleId70cIso, &b_pf_els_simpleEleId70cIso);
   fChain->SetBranchAddress("pf_els_cIso", &pf_els_cIso, &b_pf_els_cIso);
   fChain->SetBranchAddress("pf_els_tIso", &pf_els_tIso, &b_pf_els_tIso);
   fChain->SetBranchAddress("pf_els_ecalIso", &pf_els_ecalIso, &b_pf_els_ecalIso);
   fChain->SetBranchAddress("pf_els_hcalIso", &pf_els_hcalIso, &b_pf_els_hcalIso);
   fChain->SetBranchAddress("pf_els_chargedHadronIso", &pf_els_chargedHadronIso, &b_pf_els_chargedHadronIso);
   fChain->SetBranchAddress("pf_els_photonIso", &pf_els_photonIso, &b_pf_els_photonIso);
   fChain->SetBranchAddress("pf_els_neutralHadronIso", &pf_els_neutralHadronIso, &b_pf_els_neutralHadronIso);
   fChain->SetBranchAddress("pf_els_charge", &pf_els_charge, &b_pf_els_charge);
   fChain->SetBranchAddress("pf_els_hadOverEm", &pf_els_hadOverEm, &b_pf_els_hadOverEm);
   fChain->SetBranchAddress("pf_els_eOverPIn", &pf_els_eOverPIn, &b_pf_els_eOverPIn);
   fChain->SetBranchAddress("pf_els_sigmaIEtaIEta", &pf_els_sigmaIEtaIEta, &b_pf_els_sigmaIEtaIEta);
   fChain->SetBranchAddress("pf_els_scEnergy", &pf_els_scEnergy, &b_pf_els_scEnergy);
   fChain->SetBranchAddress("pf_els_scEta", &pf_els_scEta, &b_pf_els_scEta);
   fChain->SetBranchAddress("pf_els_scE1x5", &pf_els_scE1x5, &b_pf_els_scE1x5);
   fChain->SetBranchAddress("pf_els_scE2x5Max", &pf_els_scE2x5Max, &b_pf_els_scE2x5Max);
   fChain->SetBranchAddress("pf_els_scE5x5", &pf_els_scE5x5, &b_pf_els_scE5x5);
   fChain->SetBranchAddress("pf_els_isEB", &pf_els_isEB, &b_pf_els_isEB);
   fChain->SetBranchAddress("pf_els_isEE", &pf_els_isEE, &b_pf_els_isEE);
   fChain->SetBranchAddress("pf_els_dEtaIn", &pf_els_dEtaIn, &b_pf_els_dEtaIn);
   fChain->SetBranchAddress("pf_els_dPhiIn", &pf_els_dPhiIn, &b_pf_els_dPhiIn);
   fChain->SetBranchAddress("pf_els_dEtaOut", &pf_els_dEtaOut, &b_pf_els_dEtaOut);
   fChain->SetBranchAddress("pf_els_dPhiOut", &pf_els_dPhiOut, &b_pf_els_dPhiOut);
   fChain->SetBranchAddress("pf_els_numlosthits", &pf_els_numlosthits, &b_pf_els_numlosthits);
   fChain->SetBranchAddress("pf_els_tk_phi", &pf_els_tk_phi, &b_pf_els_tk_phi);
   fChain->SetBranchAddress("pf_els_d0dum", &pf_els_d0dum, &b_pf_els_d0dum);
   fChain->SetBranchAddress("pf_els_vx", &pf_els_vx, &b_pf_els_vx);
   fChain->SetBranchAddress("pf_els_vy", &pf_els_vy, &b_pf_els_vy);
   fChain->SetBranchAddress("pf_els_vz", &pf_els_vz, &b_pf_els_vz);
   fChain->SetBranchAddress("pf_els_ptError", &pf_els_ptError, &b_pf_els_ptError);
   fChain->SetBranchAddress("pf_els_n_inner_layer", &pf_els_n_inner_layer, &b_pf_els_n_inner_layer);
   fChain->SetBranchAddress("pf_els_dr03EcalRecHitSumEt", &pf_els_dr03EcalRecHitSumEt, &b_pf_els_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("pf_els_dr03HcalTowerSumEt", &pf_els_dr03HcalTowerSumEt, &b_pf_els_dr03HcalTowerSumEt);
   fChain->SetBranchAddress("pf_els_dr03HcalDepth1TowerSumEt", &pf_els_dr03HcalDepth1TowerSumEt, &b_pf_els_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("pf_els_dr03HcalDepth2TowerSumEt", &pf_els_dr03HcalDepth2TowerSumEt, &b_pf_els_dr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("pf_els_dr03TkSumPt", &pf_els_dr03TkSumPt, &b_pf_els_dr03TkSumPt);
   fChain->SetBranchAddress("pf_mus_energy", &pf_mus_energy, &b_pf_mus_energy);
   fChain->SetBranchAddress("pf_mus_et", &pf_mus_et, &b_pf_mus_et);
   fChain->SetBranchAddress("pf_mus_eta", &pf_mus_eta, &b_pf_mus_eta);
   fChain->SetBranchAddress("pf_mus_phi", &pf_mus_phi, &b_pf_mus_phi);
   fChain->SetBranchAddress("pf_mus_pt", &pf_mus_pt, &b_pf_mus_pt);
   fChain->SetBranchAddress("pf_mus_px", &pf_mus_px, &b_pf_mus_px);
   fChain->SetBranchAddress("pf_mus_py", &pf_mus_py, &b_pf_mus_py);
   fChain->SetBranchAddress("pf_mus_pz", &pf_mus_pz, &b_pf_mus_pz);
   fChain->SetBranchAddress("pf_mus_cIso", &pf_mus_cIso, &b_pf_mus_cIso);
   fChain->SetBranchAddress("pf_mus_tIso", &pf_mus_tIso, &b_pf_mus_tIso);
   fChain->SetBranchAddress("pf_mus_ecalIso", &pf_mus_ecalIso, &b_pf_mus_ecalIso);
   fChain->SetBranchAddress("pf_mus_hcalIso", &pf_mus_hcalIso, &b_pf_mus_hcalIso);
   fChain->SetBranchAddress("pf_mus_neutralHadronIso", &pf_mus_neutralHadronIso, &b_pf_mus_neutralHadronIso);
   fChain->SetBranchAddress("pf_mus_chargedHadronIso", &pf_mus_chargedHadronIso, &b_pf_mus_chargedHadronIso);
   fChain->SetBranchAddress("pf_mus_photonIso", &pf_mus_photonIso, &b_pf_mus_photonIso);
   fChain->SetBranchAddress("pf_mus_charge", &pf_mus_charge, &b_pf_mus_charge);
   fChain->SetBranchAddress("pf_mus_cm_chi2", &pf_mus_cm_chi2, &b_pf_mus_cm_chi2);
   fChain->SetBranchAddress("pf_mus_cm_ndof", &pf_mus_cm_ndof, &b_pf_mus_cm_ndof);
   fChain->SetBranchAddress("pf_mus_cm_pt", &pf_mus_cm_pt, &b_pf_mus_cm_pt);
   fChain->SetBranchAddress("pf_mus_cm_ptErr", &pf_mus_cm_ptErr, &b_pf_mus_cm_ptErr);
   fChain->SetBranchAddress("pf_mus_tk_chi2", &pf_mus_tk_chi2, &b_pf_mus_tk_chi2);
   fChain->SetBranchAddress("pf_mus_tk_ndof", &pf_mus_tk_ndof, &b_pf_mus_tk_ndof);
   fChain->SetBranchAddress("pf_mus_tk_pt", &pf_mus_tk_pt, &b_pf_mus_tk_pt);
   fChain->SetBranchAddress("pf_mus_tk_phi", &pf_mus_tk_phi, &b_pf_mus_tk_phi);
   fChain->SetBranchAddress("pf_mus_tk_d0dum", &pf_mus_tk_d0dum, &b_pf_mus_tk_d0dum);
   fChain->SetBranchAddress("pf_mus_tk_vz", &pf_mus_tk_vz, &b_pf_mus_tk_vz);
   fChain->SetBranchAddress("pf_mus_tk_numvalhits", &pf_mus_tk_numvalhits, &b_pf_mus_tk_numvalhits);
   fChain->SetBranchAddress("pf_mus_tk_ptErr", &pf_mus_tk_ptErr, &b_pf_mus_tk_ptErr);
   fChain->SetBranchAddress("pf_mus_tk_numvalPixelhits", &pf_mus_tk_numvalPixelhits, &b_pf_mus_tk_numvalPixelhits);
   fChain->SetBranchAddress("pf_mus_tk_numpixelWthMeasr", &pf_mus_tk_numpixelWthMeasr, &b_pf_mus_tk_numpixelWthMeasr);
   fChain->SetBranchAddress("pf_mus_stamu_pt", &pf_mus_stamu_pt, &b_pf_mus_stamu_pt);
   fChain->SetBranchAddress("pf_mus_stamu_ptErr", &pf_mus_stamu_ptErr, &b_pf_mus_stamu_ptErr);
   fChain->SetBranchAddress("pf_mus_num_matches", &pf_mus_num_matches, &b_pf_mus_num_matches);
   fChain->SetBranchAddress("pf_mus_isTrackerMuon", &pf_mus_isTrackerMuon, &b_pf_mus_isTrackerMuon);
   fChain->SetBranchAddress("pf_mus_isGlobalMuon", &pf_mus_isGlobalMuon, &b_pf_mus_isGlobalMuon);
   fChain->SetBranchAddress("pf_mus_id_GlobalMuonPromptTight", &pf_mus_id_GlobalMuonPromptTight, &b_pf_mus_id_GlobalMuonPromptTight);
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
   fChain->SetBranchAddress("pfmets_et", &pfmets_et, &b_pfmets_et);
   fChain->SetBranchAddress("pfmets_phi", &pfmets_phi, &b_pfmets_phi);
   fChain->SetBranchAddress("pfmets_ex", &pfmets_ex, &b_pfmets_ex);
   fChain->SetBranchAddress("pfmets_ey", &pfmets_ey, &b_pfmets_ey);
   fChain->SetBranchAddress("pfmets_gen_et", &pfmets_gen_et, &b_pfmets_gen_et);
   fChain->SetBranchAddress("pfmets_gen_phi", &pfmets_gen_phi, &b_pfmets_gen_phi);
   fChain->SetBranchAddress("pfmets_sumEt", &pfmets_sumEt, &b_pfmets_sumEt);
   fChain->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
   fChain->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
   fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   fChain->SetBranchAddress("pv_ndof", &pv_ndof, &b_pv_ndof);
   fChain->SetBranchAddress("pv_isFake", &pv_isFake, &b_pv_isFake);
   fChain->SetBranchAddress("pv_tracksSize", &pv_tracksSize, &b_pv_tracksSize);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
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
