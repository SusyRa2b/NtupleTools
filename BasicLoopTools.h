//
// Cornell University
//

//#define isMC //uncomment to run on V00-01-03 (spring11) MC
#ifdef isMC
#include "BasicLoopCU_MC.h"
#else
#include "BasicLoopCU.h"
#endif

#include <TLorentzVector.h>
#include <TVector3.h>
#include "TMatrixT.h"
#include "TRandom3.h"
#include "TMatrixDEigen.h"
#include <RooRealVar.h>
#include <RooMinuit.h>
#include <RooTransverseThrustVar.h>
#include <cassert>

#include <typeinfo>

#include "TGraphAsymmErrors.h"

//Pile-up reweighting stuff
//NOTE: Grab this header from PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h
#include "LumiReweightingStandAlone.h"

//the data histogram obtained from: 
// /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Pileup_2011_EPS_8_jul.root
float TrueDist2011_f[25] = {
  1.45417e+07,
  3.47743e+07,
  7.89247e+07,
  1.26467e+08,
  1.59329e+08,
  1.67603e+08,
  1.52684e+08,
  1.23794e+08,
  9.09462e+07,
  6.13973e+07,
  3.8505e+07,
  2.2628e+07,
  1.25503e+07,
  6.61051e+06,
  3.32403e+06,
  1.60286e+06,
  743920,
  333477,
  144861,
  61112.7,
  25110.2,
  10065.1,
  3943.98,
  1513.54,
  896.161
};
//Summer11 PU_S4 distribution
//obtained from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupMCReweightingUtilities
Double_t PoissonIntDist_f[25] = {
  0.104109,
  0.0703573,
  0.0698445,
  0.0698254,
  0.0697054,
  0.0697907,
  0.0696751,
  0.0694486,
  0.0680332,
  0.0651044,
  0.0598036,
  0.0527395,
  0.0439513,
  0.0352202,
  0.0266714,
  0.019411,
  0.0133974,
  0.00898536,
  0.0057516,
  0.00351493,
  0.00212087,
  0.00122891,
  0.00070592,
  0.000384744,
  0.000219377
};


////the LM9 PU histogram
//float  LM9DistSummer2011_f[25] = {
//  8374,
//  5764,
//  5861,
//  4747,
//  5519,
//  5914,
//  5848,
//  4697,
//  5198,
//  5759,
//  4821,
//  4314,
//  3746,
//  2805,
//  1953,
//  1457,
//  1124,
//  721,
//  443,
//  242,
//  180,
//  96,
//  48,
//  9,
//  9
//}


//#include <string>
//#include <sys/stat.h>
//#include "CondFormats/JetMETObjects/interface/JetResolution.h"

// Declare the BNN functions:
double pfmht_efficiency(std::vector<double>& inputvars,
			int first=0,
			int last=99);


using namespace std;

//should really make this a class so we don't have to do this...
const double mW_ = 80.399;
const double mtop_ = 172.0;
const double lumi_ = 1.; //fix to 1/pb and scale MC later (e.g. in drawReducedTrees)
TString sampleName_ = ""; 

//jmt -- commenting out the ones that are not currently used
// Note: when adding a new type type, make sure to give it a name is the PseudoConstructor
//enum CutScheme {kBaseline2010};
//CutScheme theCutScheme_;
enum METType {kPFMET=0,kPFMETTypeI};
METType theMETType_; std::map<METType, TString> theMETNames_;
//enum METRange {kMedium=0, kHigh, kWide, kMedhigh, kLSB};
//METRange theMETRange_;
enum jetType {kPF2PAT=0, kRECOPF};
jetType theJetType_; std::map<jetType, TString> theJetNames_;
//enum leptonType {kNormal=0, kPFLeptons, kPFLeptonsRA2};
//leptonType theLeptonType_;
//enum dpType {kDeltaPhi=0, kminDP, kMPT, kDPSync1, kminDPinv, kminDPAll30};
//dpType theDPType_;
enum JESType {kJES0=0,kJESup,kJESdown};
JESType theJESType_; std::map<JESType, TString> theJESNames_;
enum JERType {kJER0=0,kJERbias,kJERup,kJERdown,kJERra2};
JERType theJERType_; std::map<JERType, TString> theJERNames_;
enum METuncType {kMETunc0=0,kMETuncDown,kMETuncUp};
METuncType theMETuncType_; std::map<METuncType, TString> theMETuncNames_;
enum PUuncType {kPUunc0=0,kPUuncDown,kPUuncUp};
PUuncType thePUuncType_; std::map<PUuncType, TString> thePUuncNames_;
//enum BTagEffType {kBTagEff0=0,kBTagEffup,kBTagEffdown};
//BTagEffType theBTagEffType_;
//enum tailCleaningType {kNoCleaning=0, kMuonCleaning, kMuonEcalCleaning};
//tailCleaningType theCleaningType_;
//enum flavorHistoryType {kNoFlvHist=0, kFlvHist};
//flavorHistoryType theFlavorHistoryType_;
enum BTaggerType {kSSVM=0, kTCHET, kSSVHPT, kTCHPT, kTCHPM, Nbtaggers};
BTaggerType theBTaggerType_; std::map<BTaggerType, TString> theBTaggerNames_;

bool recalculatedVariables = false;

//-Define some pointers so we can use more intuitive names.
//-This also isolates some of the dependency on the configuration of the ntuple.
//--These should be checked each time we make new ntuples.
#ifdef isMC
std::vector<jet3_s> * myJetsPF;
std::vector<jet1_s> * myJetsPFhelper;
std::vector<electron3_s> * myElectronsPF;
std::vector<muon3_s> * myMuonsPF;
std::vector<tau1_s> * myTausPF;
std::vector<met1_s> * myMETPF;
std::vector<vertex_s> * myVertex;
std::vector<genparticlera2_s> * myGenParticles; //jmt
double* myGenWeight;
double* myEDM_bunchCrossing;
double* myEDM_event;
double* myEDM_isRealData;
double* myEDM_luminosityBlock;
double* myEDM_run;

void InitializeStuff(){
  myJetsPF = &jet3;
  myJetsPFhelper = &jet1;
  myElectronsPF = &electron3;
  myMuonsPF = &muon3;
  myTausPF = &tau1;
  myMETPF = &met1;
  myVertex = &vertex;
  myGenParticles = &genparticlera2; //jmt
  myGenWeight = &geneventinfoproduct1_weight;
  myEDM_bunchCrossing = &edmevent_bunchCrossing;
  myEDM_event = &edmevent_event;
  myEDM_isRealData = &edmevent_isRealData;
  myEDM_luminosityBlock = &edmevent_luminosityBlock;
  myEDM_run = &edmevent_run;
}
#else
//data
//std::vector<jet1_s> * myJetsPF;
//std::vector<jethelper_s> * myJetsPFhelper;
std::vector<jet_s> * myJetsPF;
std::vector<jetcleanhelper_s> * myJetsPFhelper;

std::vector<electron1_s> * myElectronsPF;
std::vector<muon1_s> * myMuonsPF;
std::vector<tau_s> * myTausPF;
std::vector<met1_s> * myMETPF;
//std::vector<met2_s> * myMETPF;
std::vector<vertex_s> * myVertex;
std::vector<genparticlehelperra2_s> * myGenParticles; //jmt
double * myGenWeight;
double* myEDM_bunchCrossing;
double* myEDM_event;
double* myEDM_isRealData;
double* myEDM_luminosityBlock;
double* myEDM_run;

void InitializeStuff(){
  //myJetsPF = &jet1;            //selectedPatJetsPF
  //myJetsPFhelper = &jethelper; //selectedPatJetsPF helper

  myJetsPF = &jet;            //cleanPatJetsAK5PF
  myJetsPFhelper = &jetcleanhelper; //cleanPatJetsAK5PF helper

  myElectronsPF = &electron1;  //selectedPatElectronsPF
  myMuonsPF = &muon1;          //selectedPatMuonsPF
  myTausPF = &tau;            //selectedPatTausPF
  myMETPF = &met1;             //patMETsPF
  //myMETPF = &met2;             //patMETsTypeIPF
  myVertex = &vertex;         //offlinePrimaryVertices
  myGenParticles = &genparticlehelperra2; //jmt
  myGenWeight = &geneventinfoproduct_weight;
  myEDM_bunchCrossing = &eventhelper_bunchCrossing;
  myEDM_event = &eventhelper_event;
  myEDM_isRealData = &eventhelper_isRealData;
  myEDM_luminosityBlock = &eventhelper_luminosityBlock;
  myEDM_run = &eventhelper_run;


}
#endif


float eleet1_,muonpt1_,elephi1_,muonphi1_; //FIXME this is a hack that I don't like!
bool isGoodMuon(const unsigned int imuon) {

  if(myMuonsPF->at(imuon).pt > 10
     && fabs(myMuonsPF->at(imuon).eta)<2.4
     && myMuonsPF->at(imuon).GlobalMuonPromptTight == 1
     && myMuonsPF->at(imuon).innerTrack_numberOfValidHits >=11
     && myMuonsPF->at(imuon).track_hitPattern_numberOfValidPixelHits >= 1
     && fabs(myMuonsPF->at(imuon).dB) < 0.02
     && fabs(myMuonsPF->at(imuon).vz - myVertex->at(0).z ) <1
     && (myMuonsPF->at(imuon).chargedHadronIso 
	 + myMuonsPF->at(imuon).photonIso 
	 + myMuonsPF->at(imuon).neutralHadronIso)/myMuonsPF->at(imuon).pt <0.2 
     ){
    return true;
  }
  
  return false;
}

bool isGoodElectron(const unsigned int iele) {
  
  if(myElectronsPF->at(iele).pt > 10
     && fabs(myElectronsPF->at(iele).superCluster_eta) < 2.5 
     && !(fabs(myElectronsPF->at(iele).superCluster_eta) > 1.4442 
	  && fabs(myElectronsPF->at(iele).superCluster_eta) < 1.566)
     && myElectronsPF->at(iele).gsfTrack_trackerExpectedHitsInner_numberOfLostHits <= 1
     && fabs(myElectronsPF->at(iele).dB) < 0.02
     && fabs(myElectronsPF->at(iele).vz - myVertex->at(0).z ) <1
     && (myElectronsPF->at(iele).chargedHadronIso 
	 + myElectronsPF->at(iele).photonIso 
	 + myElectronsPF->at(iele).neutralHadronIso)/myElectronsPF->at(iele).pt <0.2 
       ){
    return true;
  }
  
  return false;
}



uint countEle() {

  int ngoodele=0;
  unsigned int nele = 0;
  nele = myElectronsPF->size();
  for (unsigned int i=0; i < nele; i++) {
    if(isGoodElectron(i)){
      //FIXME adding this as a hack...i would like to do this more elegantly, but for now this will work
      if (ngoodele==0){
	eleet1_ = myElectronsPF->at(i).et;
	elephi1_ = myElectronsPF->at(i).phi;
      }
      ++ngoodele;
    }
  }
  return ngoodele;
}



uint nElecut_;
void setEleReq(unsigned int ne) {
  nElecut_=ne;
}

bool passEleVeto() {
  return (countEle() == nElecut_);
}


double deltaPhi(double phi1, double phi2) { 
  double result = phi1 - phi2;
  while (result > M_PI) result -= 2*M_PI;
  while (result <= -M_PI) result += 2*M_PI;
  return result;
}

double deltaR2(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return deta*deta + dphi*dphi;
}
double deltaR(double eta1, double phi1, double eta2, double phi2) {
  return std::sqrt(deltaR2 (eta1, phi1, eta2, phi2));
}

bool isCleanJet(const unsigned int ijet){
    
  if(theJetType_ == kPF2PAT) return true; //pf2pat jets are clean by construction

  //if it's near a good muon, it's not clean
  bool isNearMuon = false;
  for ( unsigned int j = 0; j< myMuonsPF->size(); j++) {
    if( isGoodMuon(j) && deltaR( myJetsPF->at(ijet).eta, myJetsPF->at(ijet).phi,myMuonsPF->at(j).eta, myMuonsPF->at(j).phi)<0.1 ){
      isNearMuon = true; break;
    }
  }
  
  if(isNearMuon) return false;
  
  //if it's near a good electron, it's not clean
  bool isNearElectron = false;
  for ( unsigned int j = 0; j< myElectronsPF->size(); j++) {
    if( isGoodElectron(j) && deltaR( myJetsPF->at(ijet).eta, myJetsPF->at(ijet).phi,myElectronsPF->at(j).eta, myElectronsPF->at(j).phi)<0.3 ){
      isNearElectron = true; break;
    }
  }
    
  if(isNearElectron) return false;

  return true;
}



#ifdef isMC
std::vector<jet3_s> * myJetsPF_temp;
//std::vector<electron3_s> * myElectronsPF_temp;
//std::vector<muon3_s> * myMuonsPF_temp;
//std::vector<tau1_s> * myTausPF_temp;
std::vector<met1_s> * myMETPF_temp;
//std::vector<vertex_s> * myVertex_temp;
//double* myGenWeight_temp;
//double* myEDM_bunchCrossing_temp;
//double* myEDM_event_temp;
//double* myEDM_isRealData_temp;
//double* myEDM_luminosityBlock_temp;
//double* myEDM_run_temp;


void changeVariables(TRandom* random, double jetLossProbability, int& nLostJets)
{
  if(recalculatedVariables) return;
  recalculatedVariables=true;
  myJetsPF_temp                   = myJetsPF;		  
  //myElectronsPF_temp	          = myElectronsPF;	  
  //myMuonsPF_temp		  = myMuonsPF;		  
  //myTausPF_temp		          = myTausPF;		  
  myMETPF_temp		          = myMETPF;		  
  //myVertex_temp		          = myVertex;		  
  //myGenWeight_temp	          = myGenWeight;	  
  //myEDM_bunchCrossing_temp        = myEDM_bunchCrossing;  
  //myEDM_event_temp	          = myEDM_event;	  
  //myEDM_isRealData_temp	          = myEDM_isRealData;	  
  //myEDM_luminosityBlock_temp      = myEDM_luminosityBlock;
  //myEDM_run_temp                  = myEDM_run;      

  myJetsPF                = new std::vector<jet3_s>;	  
  //myElectronsPF	          = new std::vector<electron3_s>;	  
  //myMuonsPF		  = new std::vector<muon3_s>;	  
  //myTausPF		  = new std::vector<tau1_s>;	  
  myMETPF		  = new std::vector<met1_s>(*myMETPF_temp);;	  
  //myVertex		  = new std::vector<vertex_s>;	  
  //myGenWeight	          = new double;	  
  //myEDM_bunchCrossing     = new double;
  //myEDM_event	          = new double;	  
  //myEDM_isRealData	  = new double;
  //myEDM_luminosityBlock   = new double;
  //myEDM_run               = new double;

  for(vector<jet3_s>::iterator thisJet = myJetsPF_temp->begin(); thisJet != myJetsPF_temp->end(); thisJet++)
    {
      if(random->Rndm() > jetLossProbability)
	{
	  myJetsPF->push_back(*thisJet);
	}
      else
	{
	  double jetPx = thisJet->uncor_pt * cos(thisJet->uncor_phi);
	  double jetPy = thisJet->uncor_pt * sin(thisJet->uncor_phi);
	  double METx = myMETPF->at(0).pt * cos(myMETPF->at(0).phi) - jetPx;
	  double METy = myMETPF->at(0).pt * sin(myMETPF->at(0).phi) - jetPy;
	  myMETPF->at(0).pt = sqrt(METx*METx + METy*METy);
	  myMETPF->at(0).phi = atan2(METy,METx);
	}
    }
  nLostJets = myJetsPF_temp->size() - myJetsPF->size();
}

#else
std::vector<jet_s> * myJetsPF_temp;
//std::vector<electron1_s> * myElectronsPF_temp;
//std::vector<muon1_s> * myMuonsPF_temp;
//std::vector<tau1_s> * myTausPF_temp;
std::vector<met1_s> * myMETPF_temp;
//std::vector<vertex_s> * myVertex_temp;
//double* myGenWeight_temp;
//double* myEDM_bunchCrossing_temp;
//double* myEDM_event_temp;
//double* myEDM_isRealData_temp;
//double* myEDM_luminosityBlock_temp;
//double* myEDM_run_temp;

void changeVariables(TRandom* random, double jetLossProbability, int& nLostJets)
{
  if(recalculatedVariables) return;
  recalculatedVariables=true;
  myJetsPF_temp                   = myJetsPF;		  
  //myElectronsPF_temp	          = myElectronsPF;	  
  //myMuonsPF_temp		  = myMuonsPF;		  
  //myTausPF_temp		          = myTausPF;		  
  myMETPF_temp		          = myMETPF;		  
  //myVertex_temp		          = myVertex;		  
  //myGenWeight_temp	          = myGenWeight;	  
  //myEDM_bunchCrossing_temp        = myEDM_bunchCrossing;  
  //myEDM_event_temp	          = myEDM_event;	  
  //myEDM_isRealData_temp	          = myEDM_isRealData;	  
  //myEDM_luminosityBlock_temp      = myEDM_luminosityBlock;
  //myEDM_run_temp                  = myEDM_run;         

  myJetsPF                = new std::vector<jet_s>;	  
  //myElectronsPF	          = new std::vector<electron1_s>;	  
  //myMuonsPF		  = new std::vector<muon1_s>;	  
  //myTausPF		  = new std::vector<tau1_s>;	  
  myMETPF		  = new std::vector<met1_s>(*myMETPF_temp);	  
  //myVertex		  = new std::vector<vertex_s>;	  
  //myGenWeight	          = new double;	  
  //myEDM_bunchCrossing     = new double;
  //myEDM_event	          = new double;	  
  //myEDM_isRealData	  = new double;
  //myEDM_luminosityBlock   = new double;
  //myEDM_run               = new double;
  for(vector<jet_s>::iterator thisJet = myJetsPF_temp->begin(); thisJet != myJetsPF_temp->end(); thisJet++)
    {
      if(random->Rndm() > jetLossProbability)
	{
	  myJetsPF->push_back(*thisJet);
	}
      else
	{
	  double jetPx = thisJet->uncor_pt * cos(thisJet->uncor_phi);
	  double jetPy = thisJet->uncor_pt * sin(thisJet->uncor_phi);
	  double METx = myMETPF->at(0).pt * cos(myMETPF->at(0).phi) - jetPx;
	  double METy = myMETPF->at(0).pt * sin(myMETPF->at(0).phi) - jetPy;
	  myMETPF->at(0).pt = sqrt(METx*METx + METy*METy);
	  myMETPF->at(0).phi = atan2(METy,METx);
	}
    }
  nLostJets = myJetsPF_temp->size() - myJetsPF->size();
}

#endif

void resetVariables()
{
  if(!recalculatedVariables) return;
  recalculatedVariables = false;
  if(myJetsPF != 0) delete myJetsPF;
  else cout << "you've done something wrong!" << endl;
  //delete myElectronsPF;	
  //delete myMuonsPF;		
  //delete myTausPF;		
  if(myMETPF != 0)delete myMETPF;		
  else cout << "you've done something wrong!" << endl;
  //delete myVertex;		
  //delete myGenWeight;	  	
  //delete myEDM_bunchCrossing;  
  //delete myEDM_event;	  	
  //delete myEDM_isRealData;	
  //delete myEDM_luminosityBlock;
  //delete myEDM_run;        

  myJetsPF                = myJetsPF_temp;		  
  //myElectronsPF	          = myElectronsPF_temp;	  
  //myMuonsPF		  = myMuonsPF_temp;		  
  //myTausPF		  = myTausPF_temp;		  
  myMETPF		  = myMETPF_temp;		  
  //myVertex		  = myVertex_temp;		  
  //myGenWeight	          = myGenWeight_temp;	  
  //myEDM_bunchCrossing     = myEDM_bunchCrossing_temp;  
  //myEDM_event	          = myEDM_event_temp;	  
  //myEDM_isRealData	  = myEDM_isRealData_temp;	  
  //myEDM_luminosityBlock   = myEDM_luminosityBlock_temp;
  //myEDM_run               = myEDM_run_temp;      
    
}

void checkConsistency() {
 //this is probably bad coding, but the only thing I can think of for this already horribly-constructed code
  std::string jettype = typeid( *myJetsPF ).name();
  std::string recopfjet = "jet_s";
  std::string pf2patjet = "jet1_s";
  if( theJetType_==kPF2PAT && std::string::npos != jettype.find(recopfjet)){ 
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
  std::string type1pfmet = "met2_s";
  if( theMETType_ == kPFMETTypeI &&  std::string::npos != jettype.find(type1pfmet)){
    std::cout << "ERROR: theMETType_ is set to PFMETTypeI, while myMETsPF is also pointing to type1pfmet.  "
	      << "This will double-correct the MET!  Aborting." << std::endl;
    assert(0);
  }

}


void PseudoConstructor() {
  //  theCutScheme_ = kBaseline2010;
  theMETType_=kPFMETTypeI;
  //  theMETType_=kPFMET;
  //  theMETRange_=kHigh;
  theJetType_=kRECOPF;
  //theJetType_=kPF2PAT; 

  checkConsistency();

  //  theLeptonType_=kPFLeptonsRA2;
  //  theDPType_=kminDP;
  theJESType_=kJES0;
  theJERType_=kJER0;
  theMETuncType_=kMETunc0;
  thePUuncType_=kPUunc0;
  //  theBTagEffType_=kBTagErr0;
  //  theCleaningType_=kNoCleaning;
  //  theFlavorHistoryType_=kNoFlvHist;
  theBTaggerType_=kSSVHPT;

  theBTaggerNames_[kSSVM]="SSVHEM";
  theBTaggerNames_[kTCHET]="TCHET";
  theBTaggerNames_[kSSVHPT]="SSVHPT";
  theBTaggerNames_[kTCHPT]="TCHPT";
  theBTaggerNames_[kTCHPM]="TCHPM";

  theMETNames_[kPFMET] = "PFMET";
  theMETNames_[kPFMETTypeI] = "PFMETTypeI";

  theJetNames_[kPF2PAT]="PF2PATjets";
  theJetNames_[kRECOPF]="RecoPFjets";

  theJESNames_[kJES0]="JES0";
  theJESNames_[kJESup]="JESup";
  theJESNames_[kJESdown]="JESdown";

  theMETuncNames_[kMETunc0]="METunc0";
  theMETuncNames_[kMETuncUp]="METuncUp";
  theMETuncNames_[kMETuncDown]="METuncDown";

  thePUuncNames_[kPUunc0]="PUunc0";
  thePUuncNames_[kPUuncUp]="PUuncUp";
  thePUuncNames_[kPUuncDown]="PUuncDown";

  theJERNames_[kJER0]="JER0";
  theJERNames_[kJERup]="JERup";
  theJERNames_[kJERbias]="JERbias";
  theJERNames_[kJERra2]="JERra2";
  theJERNames_[kJERdown]="JERdown";

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//#ifdef isMC
unsigned int findSUSYMaternity( unsigned int k ) {

  int numSUSY=0;
  //cout<<"[findSUSYMaternity]"<<endl;
  //  if  (TMath::Nint(myGenParticles->at(k).firstMother) != TMath::Nint(myGenParticles->at(k).lastMother )) cout<<"Multiple mothers!"<<endl;

  if (TMath::Nint(myGenParticles->at(k).firstMother) <0) {/*cout<<"got no mother!"<<endl;*/ return numSUSY;}

  int motherId = abs(TMath::Nint( myGenParticles->at(  TMath::Nint(myGenParticles->at(k).firstMother)).pdgId )); //get rid of minus signs!
  if ( (motherId>= 1000001 && motherId <=1000037) || (motherId>= 2000001 && motherId <=2000015)) {numSUSY++;}
  else if (motherId == 21) {/*cout<<"Found gluon splitting!"<<endl;*/}
  else { numSUSY+= findSUSYMaternity(TMath::Nint(myGenParticles->at(k).firstMother)  ); }

  return numSUSY;

}

unsigned int getSUSYnb(std::vector<uint> &susyb_index) {

  unsigned int SUSY_nb=0;
  for (unsigned int k = 0; k<myGenParticles->size(); k++) {
    //         cout<<k<<"\t"<<TMath::Nint(myGenParticles->at(k).pdgId )<<"\t"<< TMath::Nint(myGenParticles->at(k).firstMother)<<"\t"<<TMath::Nint(myGenParticles->at(k).lastMother)<<"\t"<<TMath::Nint(myGenParticles->at(k).status ) <<endl;

    if ( abs(TMath::Nint(myGenParticles->at(k).pdgId )) == 5 ) { //find b quark
      unsigned int nSUSY = findSUSYMaternity( k );      
      if (nSUSY>0) {SUSY_nb++; susyb_index.push_back(k);}
    }
  }
  
  return SUSY_nb;

}

////#else
//
//unsigned int findSUSYMaternity( unsigned int k ) {
//  return 0;
//}
//
//unsigned int getSUSYnb() {
//  return 0;
//}
//
////#endif
/*
void testLoop(itreestream& stream) {
  InitializeStuff();//required!


  int nevents = stream.size();

  for(int entry=0; entry < nevents; ++entry){
    // Read event into memory
    stream.read(entry);
    fillObjects();

    cout<<" ==== "<<endl;
//    analyzeDecayTree();
  }
}
*/

TString getCutDescriptionString(){
  
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

  return cut;
}

void setSampleName_(TString name){
  sampleName_ = name;
}

void setBTaggerType(BTaggerType btaggertype){
  theBTaggerType_ = btaggertype;
}

int isRealDataInt(double * isRealDataDouble){
  return ((int)((*isRealDataDouble)+0.5));
}


//This is a function which is temporarily hard-coded to remove data that has become uncertified
bool passLumiMask(){

  int lumi =  TMath::Nint( *myEDM_luminosityBlock );
  int run = TMath::Nint( *myEDM_run );

  if(run == 165525 && lumi>=1 && lumi<=31) return false;
  if(run == 165537 && lumi>=1 && lumi<=20) return false;
  if(run == 165537 && lumi>=22 && lumi<295) return false;

  return true;

}


bool passHLT() { 

  long runnumber = (long)((*myEDM_run) +0.5); //just in case some crazy thing happens to myEDM_run being saved as a double

  //RA2b - 2011 Triggers
  bool passTrig = false;

  if(isRealDataInt(myEDM_isRealData)){
    #ifndef isMC
    //edmtriggerresults is a double - possible values are 0 (failed trigger), 1 (passed trigger), -9999 (trig result not available)
    /* using BTagIP triggers
    if(runnumber >= 160431 && runnumber < 161205) passTrig = (edmtriggerresults_HLT_HT260_MHT60_v2 > 0);
    else if (runnumber >= 161205 && runnumber < 163269) passTrig = (edmtriggerresults_HLT_HT250_MHT60_v2 > 0);
    else if (runnumber >= 163269 && runnumber < 164924) passTrig = (edmtriggerresults_HLT_HT250_MHT60_v3 > 0);
    else if (runnumber >= 164924 && runnumber < 165922) passTrig = edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v2;
    else if (runnumber >= 165922 && runnumber < 166301) passTrig = edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3;
    //else if (runnumber >= 166301 && runnumber < 166374 ) passTrig = edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v4; 
    else if (runnumber >= 166374 && runnumber < 167078) passTrig = edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3;
    //else if (runnumber >= 167078) passTrig = edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v5; 
    */
    //jmt -- 6 July 2011 -- use agreed RA2b triggers for Summer 2011
    if      (runnumber >= 160431 && runnumber <= 161204)  passTrig = (triggerresultshelper_HLT_HT260_MHT60_v2 > 0);
    else if (runnumber >= 161205 && runnumber <= 163268)  passTrig = (triggerresultshelper_HLT_HT250_MHT60_v2 > 0);
    else if (runnumber >= 163269 && runnumber <= 164923)  passTrig = (triggerresultshelper_HLT_HT250_MHT60_v3 > 0);   
    else if (runnumber >= 164924 && runnumber <= 165921)  passTrig = (triggerresultshelper_HLT_HT250_MHT70_v1 > 0);
    else if (runnumber >= 165922 && runnumber <= 166300)  passTrig = (triggerresultshelper_HLT_HT300_MHT75_v7 > 0);
    else if (runnumber >= 166301 && runnumber <= 166373)  passTrig = (triggerresultshelper_HLT_HT300_MHT75_v8 > 0);
    else if (runnumber >= 166374 && runnumber <= 166978)  passTrig = (triggerresultshelper_HLT_HT300_MHT75_v7 > 0);
    //168941 is an arbitrary run during July Tech Stop		      
    else if (runnumber >= 166979 && runnumber <= 168941)  passTrig = (triggerresultshelper_HLT_HT300_MHT80_v1 > 0);
    
    else {cout<<"No trigger assigned for run = "<<runnumber<<endl; assert(0);}
    #endif
  }
  else passTrig = true;   //use no trigger for MC   

  return passTrig;
}

unsigned int utilityHLT_HT300(){
  if(!isRealDataInt(myEDM_isRealData)) return 999;
  
  unsigned int passTrig = 0;
  #ifndef isMC
  if(triggerresultshelper_HLT_HT300_v1>0) passTrig = 1;
  if(triggerresultshelper_HLT_HT300_v2>0) passTrig = 2;
  if(triggerresultshelper_HLT_HT300_v3>0) passTrig = 3;
  if(triggerresultshelper_HLT_HT300_v4>0) passTrig = 4;
  if(triggerresultshelper_HLT_HT300_v5>0) passTrig = 5;
  if(triggerresultshelper_HLT_HT300_v6>0) passTrig = 6;
  if(triggerresultshelper_HLT_HT300_v7>0) passTrig = 7;
  if(triggerresultshelper_HLT_HT300_v8>0) passTrig = 8;
  #endif
  return passTrig;
}

unsigned int utilityHLT_HT300_CentralJet30_BTagIP(){
  if(!isRealDataInt(myEDM_isRealData)) return 999;

  unsigned int passTrig = 0;
  #ifndef isMC
  if(triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v2>0) passTrig = 2;
  if(triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v3>0) passTrig = 3; 
  if(triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v4>0) passTrig = 4; 
  if(triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v5>0) passTrig = 5; 
  #endif
  return passTrig;
}


float getHLTHTeff(float offHT) {

  float eff=100;
  //get this file from ~/public/wteo/
  TFile f_eff("ht300_out_166864_eff_histogram.root","READ");

  char graphname[200];
  std::string grname = "myeff";
  sprintf(graphname,"%s",grname.c_str());   
  TGraphAsymmErrors * htgraph  = (TGraphAsymmErrors *)f_eff.Get(graphname);

  eff = htgraph->Eval(offHT);

  return eff;
}

int countGoodPV() {
  
  int ngood=0;
  for (unsigned int ipv = 0; ipv<myVertex->size(); ipv++) {
    if ( !myVertex->at(ipv).isFake ) {
      if ( fabs(myVertex->at(ipv).z) < 24 ){
	if ( fabs(myVertex->at(ipv).position_Rho) < 2 ){
	  if ( myVertex->at(ipv).ndof > 4 ){
	    ++ngood;
	  }
	}
      }
    }
  }

  return ngood; 
}

bool passPV() {
  
  bool npass=false;
  //for (unsigned int ipv = 0; ipv<myVertex->size(); ipv++) {
  if ( !myVertex->at(0).isFake ){
    if ( fabs(myVertex->at(0).z) < 24 ){
      if ( fabs(myVertex->at(0).position_Rho) < 2 ){
	if ( myVertex->at(0).ndof > 4 ){
	  npass = true;
	}
      }
    }
  }
  
  return npass; 

}

float getJERbiasFactor(unsigned int ijet) {
  float abseta = fabs(myJetsPF->at(ijet).eta);

 
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

float getJetPt( unsigned int ijet ) {

  //i am going to write this to allow simultaneous use of JER and JES
  //(hence use if repeated 'if' instead of 'else if')
  //I hope this is sensible!

  float pt = myJetsPF->at(ijet).pt;
  //  if ( theJESType_ == kJES0 && theJERType_ == kJER0) return pt;

  //first JER
  if ( theJERType_ != kJER0 ) {
    float genpt = myJetsPF->at(ijet).genJet_pt; //not sure what it will be called
    if (genpt > 15) {
      float factor = getJERbiasFactor(ijet);
      float deltapt = (pt - genpt) * factor;
      float frac = (pt+deltapt)/pt;
      float ptscale = frac>0 ? frac : 0;
      pt *= ptscale;
    }

  } 
  //then JES
  if ( theJESType_ == kJESup ) {
    //in 2010 there was an extra term added in quadrature. 
    //i will not implement that because i don't know if that term should exist in 2011
    //note that this if statement is important because there are dummy values like -1000 sometimes
    if (  fabs(myJetsPFhelper->at(ijet).jetUncPlus)<1 )   pt *= (1+myJetsPFhelper->at(ijet).jetUncPlus);
  }
  else if (theJESType_ == kJESdown) {
    if (  fabs(myJetsPFhelper->at(ijet).jetUncMinus)<1 ) pt *= (1-myJetsPFhelper->at(ijet).jetUncMinus);
  }

  return pt;

}

float getJetPx( unsigned int ijet ) {
  return getJetPt(ijet) * cos(myJetsPF->at(ijet).phi);
}

float getJetPy( unsigned int ijet ) {
  return getJetPt(ijet) * sin(myJetsPF->at(ijet).phi);
}

float getJetPz( unsigned int ijet ) {
  return getJetPt(ijet) * sinh(myJetsPF->at(ijet).eta);
}

float getJetEnergy( unsigned int ijet ) {
  //need to move to using jecFactor 
  return myJetsPF->at(ijet).energy;
}

bool jetPassLooseID( unsigned int ijet ) {

  if(myJetsPF->at(ijet).neutralHadronEnergyFraction < 0.99
     && myJetsPF->at(ijet).neutralEmEnergyFraction < 0.99
     && myJetsPF->at(ijet).numberOfDaughters > 1
     && ( fabs(myJetsPF->at(ijet).uncor_eta)>=2.4 
	  || (fabs(myJetsPF->at(ijet).uncor_eta)<2.4 && myJetsPF->at(ijet).chargedHadronEnergyFraction>0))
     && ( fabs(myJetsPF->at(ijet).uncor_eta)>=2.4 
	  || (fabs(myJetsPF->at(ijet).uncor_eta)<2.4 && myJetsPF->at(ijet).chargedEmEnergyFraction<0.99))
     && ( fabs(myJetsPF->at(ijet).uncor_eta)>=2.4 
	  || (fabs(myJetsPF->at(ijet).uncor_eta)<2.4 && myJetsPF->at(ijet).chargedMultiplicity>0))    
     ){
    return true;
  }

  return false;
}


bool isGoodJet(const unsigned int ijet, const float pTthreshold=50, const float etaMax=2.4) {

  if ( getJetPt(ijet) <pTthreshold) return false;
  if ( fabs(myJetsPF->at(ijet).eta) > etaMax) return false;
  if ( !jetPassLooseID(ijet) ) return false;

  //do manual cleaning for reco pfjets
  if(theJetType_ == kRECOPF){
    if ( !isCleanJet(ijet) ) return false;
  }

  return true;
}

bool isGoodJet10(unsigned int ijet) {

  return isGoodJet(ijet,10,2.4);
}

bool isGoodJet30(unsigned int ijet) {
  return isGoodJet(ijet,30,2.4);
}

bool isGoodJetMHT(unsigned int ijet) {

  if ( getJetPt(ijet) <30) return false;
  if ( fabs(myJetsPF->at(ijet).eta) > 5) return false;
  //no jet id for MHT

  return true;
}


bool passBTagger(int ijet, BTaggerType btagger=Nbtaggers ) {

  //if we are passed the dummy value, then set to the default type
  if (btagger==Nbtaggers) btagger = theBTaggerType_;

  if (btagger==kSSVM){
    return ( myJetsPF->at(ijet).simpleSecondaryVertexBJetTags >= 1.74 
	     || myJetsPF->at(ijet).simpleSecondaryVertexHighEffBJetTags >= 1.74);
  }
  else if (btagger==kTCHET){
    return (myJetsPF->at(ijet).trackCountingHighEffBJetTags >= 10.2);
  }
  else if (btagger==kSSVHPT) return myJetsPF->at(ijet).simpleSecondaryVertexHighPurBJetTags >= 2;
  else if (btagger==kTCHPT ) return myJetsPF->at(ijet).trackCountingHighPurBJetTags >= 3.41;
  else if (btagger==kTCHPM ) return myJetsPF->at(ijet).trackCountingHighPurBJetTags >= 1.93;
  else{
    cout << "Invalid b tagger!" << endl;
    assert(0);
    return false;
  }
}

uint nGoodJets() {
  
  uint njets=0;
  for (unsigned int i=0; i < myJetsPF->size(); ++i) {
    if (isGoodJet(i) )   njets++;
  }
  return njets;
}


uint nGoodBJets( BTaggerType btagger=Nbtaggers) {
  uint nb=0;
  for (uint i = 0; i < myJetsPF->size(); ++i) {
    if (isGoodJet(i,30) ) {
      if ( passBTagger(i,btagger) ) nb++;
    }
  }
  return nb;
}



bool isCleanMuon(const unsigned int imuon){

  if(!isGoodMuon(imuon)) return false;

  //clean muons if using reco-pfjets
  if(theJetType_ == kRECOPF){    
    bool isNearJet = false;
    for ( unsigned int j = 0; j< myJetsPF->size(); j++) {
      if( isGoodJet(j) && deltaR( myJetsPF->at(j).eta, myJetsPF->at(j).phi,myMuonsPF->at(imuon).eta, myMuonsPF->at(imuon).phi)<0.3 ){
	isNearJet = true ; break;
      }
    }
    if(isNearJet) return false;
  }

  return true;

}


uint countMu() {

  int ngoodmu=0;
  unsigned int nmu = 0;
  nmu = myMuonsPF->size();
  for ( unsigned int i = 0; i< nmu; i++) {
    if(isCleanMuon(i)){
      //once we reach here we've got a good muon in hand
      //FIXME adding this as a hack...i would like to do this more elegantly, but for now this will work
      if (ngoodmu==0){
	muonpt1_ = myMuonsPF->at(i).pt;
	muonphi1_ = myMuonsPF->at(i).phi;
      }
      ++ngoodmu;   
    }
  }
  
  return ngoodmu;
}

uint nMucut_;
void setMuonReq(unsigned int nmu) {
  nMucut_=nmu;
}

bool passMuVeto() {
  return (countMu() == nMucut_);
}

bool passBadPFMuonFilter(){
  bool hasBadMuon = false;

  //unsigned int nmu = 0;
  //nmu = myMuonsPF->size();
  //
  //for (unsigned int i=0; i < nmu; i++) {
  //  //if(muonhelper.at(i).badPFmuon) hasBadMuon = true; NEED TO GET THIS BACK IN BEN FIXME
  //}

  if(hasBadMuon) return false;
  return true;
}

bool passInconsistentMuonFilter(){

  bool hasInconsistentMuon = false;

  //this should not be looping over selectedPatMuonsPF (which has some cuts applied), it should be looping over all PFCandidate muons
  //use the ra2 code
  /*
  unsigned int nmu = 0;
  nmu = myMuonsPF->size();

  for (unsigned int i=0; i < nmu; i++) {
    if(myMuonsPF->at(i).pt > 100){
      if(myMuonsPF->at(i).isTrackerMuon &&  myMuonsPF->at(i).isGlobalMuon){
	if( fabs( myMuonsPF->at(i).innerTrack_pt/myMuonsPF->at(i).globalTrack_pt - 1 ) > 0.1 )
	  hasInconsistentMuon = true;
      }
    }
  }
  */

  if(hasInconsistentMuon) return false;
  return true;

}


float getBTagIPWeight() {//this function should be called *after* offline tagging 
  float w_event = 1;

  float ProbGEQ1 = 1, Prob0 = 1;
  //float Prob1 = 0;

  //From TTbarMC, the BTagIP efficiency on offline tagged b-jets (non b-jets) 
  //is roughly constant and equal to:
  float w_bjet = 0.88;
  float w_nonbjet = 0.61;

  //Count the number of offline-tagged b-jets and nonb-jets
  uint nb=0, nnonb=0;
  for (uint i = 0; i < myJetsPF->size(); ++i) {
    if (isGoodJet(i) ) {

      float effi  = 0;

      int flavor = myJetsPF->at(i).partonFlavour;
      if ( passBTagger(i) && abs(flavor)==5){
	nb++;
	effi = w_bjet;
      }
      else if ( passBTagger(i) && 
		(abs(flavor)==4 || abs(flavor)==3 || abs(flavor)==2 
		 || abs(flavor)==1 || abs(flavor)==21)){
	nnonb++;
	effi = w_nonbjet;
      }

      Prob0 = Prob0* ( 1 - effi);
      
      //double product = 1;
      //for (uint j=0; j<myJetsPF->size(); ++j) {
      //	if(isGoodJet(j)){
      //	  float effj = 0;
      //	  int flavorj = myJetsPF->at(j).partonFlavour;
      //	  if ( passBTagger(j) && abs(flavor)==5){
      //	    effj = w_bjet;
      //	  }
      //	  else if ( passBTagger(j) && 
      //		    (abs(flavor)==4 || abs(flavor)==3 || abs(flavor)==2 
      //		     || abs(flavor)==1 || abs(flavor)==21)){
      //	    effj = w_nonbjet;
      //	  }
      //	  if(j != i) product = product*(1-effj);	  
      //	}
      //}
      //Prob1 += effi*product;
    }
  }

  ProbGEQ1 = 1 - Prob0;
  w_event = ProbGEQ1;

  if(nb==0 && nnonb==0){
    return 0; //this event is not even offline tagged
  }

  return w_event;
}


float jetPtOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<myJetsPF->size(); i++) {

    bool pass=false;
    pass = isGoodJet(i);

    if (pass ) {
      ngood++;
      if (ngood==n) return  getJetPt(i);
    }
  }
  return 0;
}

float jetPhiOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<myJetsPF->size(); i++) {

    bool pass=false;
    pass = isGoodJet(i);

    if (pass ) {
      ngood++;
      if (ngood==n) return  myJetsPF->at(i).phi;
    }
  }
  return 0;
}

float jetEtaOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<myJetsPF->size(); i++) {

    bool pass=false;
    pass = isGoodJet(i);

    if (pass ) {
      ngood++;
      if (ngood==n) return  myJetsPF->at(i).eta;
    }
  }
  return 0;
}

float jetEnergyOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<myJetsPF->size(); i++) {

    bool pass=false;
    pass = isGoodJet(i);

    if (pass ) {
      ngood++;
      if (ngood==n) return  getJetEnergy(i);
    }
  }
  return 0;
}


float bjetPtOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<myJetsPF->size(); i++) {

    bool pass=false;
    pass = (isGoodJet(i) && passBTagger(i));

    if (pass ) {
      ngood++;
      if (ngood==n) return  getJetPt(i);
    }
  }
  return 0;
}

float bjetPhiOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<myJetsPF->size(); i++) {

    bool pass=false;
    pass = (isGoodJet(i) && passBTagger(i));

    if (pass ) {
      ngood++;
      if (ngood==n) return  myJetsPF->at(i).phi;
    }
  }
  return 0;
}

float bjetEtaOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<myJetsPF->size(); i++) {

    bool pass=false;
    pass = (isGoodJet(i) && passBTagger(i));

    if (pass ) {
      ngood++;
      if (ngood==n) return  myJetsPF->at(i).eta;
    }
  }
  return 0;
}

float bjetEnergyOfN(unsigned int n) {

  unsigned int ngood=0;
  for (unsigned int i=0; i<myJetsPF->size(); i++) {

    bool pass=false;
    pass = (isGoodJet(i) && passBTagger(i));

    if (pass ) {
      ngood++;
      if (ngood==n) return  getJetEnergy(i);
    }
  }
  return 0;
}

float getHT() {
  float ht=0;
  for (unsigned int i=0; i<myJetsPF->size(); i++) {
    if (isGoodJet( i ) ) ht+= getJetPt(i);
  }
  return ht;
}

void getCorrectedMET(float& correctedMET, float& correctedMETPhi) {

  double METx = myMETPF->at(0).pt * cos(myMETPF->at(0).phi);
  double METy = myMETPF->at(0).pt * sin(myMETPF->at(0).phi);

  for(unsigned int thisJet = 0; thisJet < myJetsPF->size(); thisJet++)
    {
      if (isCleanJet(thisJet) ){//this only has an effect for recopfjets
	METx += myJetsPF->at(thisJet).uncor_pt * cos(myJetsPF->at(thisJet).uncor_phi);
	METx -= getJetPt(thisJet) * cos(myJetsPF->at(thisJet).phi);
	METy += myJetsPF->at(thisJet).uncor_pt * sin(myJetsPF->at(thisJet).uncor_phi);
	METy -= getJetPt(thisJet) * sin(myJetsPF->at(thisJet).phi);
      }
    }
  correctedMET = sqrt(METx*METx + METy*METy);
  correctedMETPhi = atan2(METy,METx);
}


void dumpEvent () {

  cout<<" == jets"<<endl;
  for (unsigned int i=0; i<myJetsPF->size(); i++) {
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


std::pair<float,float> getSmearedUnclusteredMETxy() {

  float myMET=0;
  float myMETphi=0;

  //in both cases start with uncorrected MET
  if (theMETType_== kPFMET || theMETType_==kPFMETTypeI) {
    myMET = myMETPF->at(0).pt;
    myMETphi = myMETPF->at(0).phi;
  }
  else {assert(0);}

  float myMETx = myMET * cos(myMETphi);
  float myMETy = myMET * sin(myMETphi);

  //  float a,b; getCorrectedMET(a,b);
  //  cout<<a<<"\t";
 
  //first remove jets from MET
  for (unsigned int i=0; i<myJetsPF->size(); i++) {
    if (isCleanJet(i) ){//this only has an effect for recopfjets    
      //remove the uncorrected jet pT from the raw MET
      myMETx += myJetsPF->at(i).uncor_pt * cos(myJetsPF->at(i).uncor_phi);
      myMETy += myJetsPF->at(i).uncor_pt * sin(myJetsPF->at(i).uncor_phi);
    }
  }
  //then muons
  for ( unsigned int i = 0; i<myMuonsPF->size() ; i++) {
    if (isCleanMuon(i)) {
      myMETx += myMuonsPF->at(i).pt * cos(myMuonsPF->at(i).phi);
      myMETy += myMuonsPF->at(i).pt * sin(myMuonsPF->at(i).phi);
    }
  }
  //electrons
  for ( unsigned int i = 0; i<myElectronsPF->size() ; i++) {
    if (isGoodElectron(i)) {
      myMETx += myElectronsPF->at(i).pt * cos(myElectronsPF->at(i).phi);
      myMETy += myElectronsPF->at(i).pt * sin(myElectronsPF->at(i).phi);
    }
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

  //note that this time we use corrected jets in order to get corrected MET

  //jets
  for (unsigned int i=0; i<myJetsPF->size(); i++) {
    if (isCleanJet(i) ){//this only has an effect for recopfjets    
      //remove the corrected jet pT from MET
      myMETx -= myJetsPF->at(i).pt * cos(myJetsPF->at(i).phi);
      myMETy -= myJetsPF->at(i).pt * sin(myJetsPF->at(i).phi);
    }
  }
  //muons
  for ( unsigned int i = 0; i<myMuonsPF->size() ; i++) {
    if (isCleanMuon(i)) {
      myMETx -= myMuonsPF->at(i).pt * cos(myMuonsPF->at(i).phi);
      myMETy -= myMuonsPF->at(i).pt * sin(myMuonsPF->at(i).phi);
    }
  }
  //electrons
  for ( unsigned int i = 0; i<myElectronsPF->size() ; i++) {
    if (isGoodElectron(i)) {
      myMETx -= myElectronsPF->at(i).pt * cos(myElectronsPF->at(i).phi);
      myMETy -= myElectronsPF->at(i).pt * sin(myElectronsPF->at(i).phi);
    }
  }

  //  cout<<sqrt(myMETx*myMETx + myMETy*myMETy)<<endl;
  return make_pair(myMETx,myMETy);
}


float getMET() {
  float myMET=-1;
  float myMETphi =-1000;

  if(theMETType_== kPFMET){ 
    assert(theJESType_==kJES0 && theJERType_==kJER0);
    myMET = myMETPF->at(0).pt;
  }
  else if(theMETType_ == kPFMETTypeI){
    getCorrectedMET(myMET,myMETphi);
  }

  //JER and JES are automatically handled inside getCorrectedMET()
  if (theMETuncType_!=kMETunc0) {
    std::pair<float, float> metxy = getSmearedUnclusteredMETxy();
    myMET = sqrt( metxy.first*metxy.first + metxy.second*metxy.second);
  }
  return myMET;
}

float getMETphi() {
  float myMET=-1;
  float myMETphi=-99;

  if(theMETType_== kPFMET){
    myMETphi = myMETPF->at(0).phi;
  }
  else if(theMETType_ == kPFMETTypeI){
    getCorrectedMET(myMET,myMETphi);
  }

  if (theMETuncType_!=kMETunc0) {
    
    std::pair<float, float> metxy = getSmearedUnclusteredMETxy();
    myMETphi = atan2( metxy.second, metxy.first);
  }

  return myMETphi;
}

double getDeltaPhi(double phi1, double phi2) {
  return acos(cos(phi1-phi2));
}

double getMinDeltaPhiMET(unsigned int maxjets) {

  double mindp=99;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< myJetsPF->size(); i++) {
    
    if (isGoodJet(i)) {
      ++ngood;
      double dp =  getDeltaPhi( myJetsPF->at(i).phi , getMETphi());
      if (dp<mindp) mindp=dp;
      if (ngood >= maxjets) break;
    }
  }
  
  return mindp;
}

double getTransverseMETError(unsigned int thisJet){

  if(!(thisJet<myJetsPF->size())) return -99;

  double deltaT=0;

  for (unsigned int i=0; i< myJetsPF->size(); i++) {
    
    if (isGoodJet30(i) && i!=thisJet) {
      double dp = getDeltaPhi( myJetsPF->at(i).phi , myJetsPF->at(thisJet).phi );
      deltaT+=pow(getJetPt(i)*sin(dp),2);
    }
  }
  return 0.1*sqrt(deltaT);
}

//double getTransverseMETErrorWithCorrections(unsigned int thisJet){
//
//  double deltaT=0;
//
//  string   era("Spring10");
//  string   alg("AK5PF");
//  bool     doGaussian(true);
//
//  string cmssw_base(getenv("CMSSW_BASE"));
//  string cmssw_release_base(getenv("CMSSW_RELEASE_BASE"));
//  string path = cmssw_base + "/src/CondFormats/JetMETObjects/data";
//  struct stat st;
//  if (stat(path.c_str(),&st)!=0)
//    path = cmssw_release_base + "/src/CondFormats/JetMETObjects/data";
//  if (stat(path.c_str(),&st)!=0) {
//    cerr<<"ERROR: tried to set path but failed, abort."<<endl;
//    return 0;
//  }
//  
//  string ptFileName  = path + "/" + era + "_PtResolution_" +alg+".txt";
//
//  JetResolution ptResol(ptFileName,doGaussian);
//
//  for (unsigned int i=0; i< myJetsPF->size(); i++) {
//    
//    if (isGoodJet(i) && i!=thisJet) {
//      double dp = getDeltaPhi( myJetsPF->at(i).phi , myJetsPF->at(thisJet).phi );
//      //Put in the correct energy corrections!
//      double ptResolution = ptResol.parameterEta( "sigma", myJetsPF->at(i).eta )->Eval( getJetPt(i) );
//      deltaT+=pow(ptResolution*getJetPt(i)*sin(dp),2);
//    }
//  }
//  return sqrt(deltaT);
//}

double getDeltaPhiNMET(unsigned int thisJet){

  if(!(thisJet<myJetsPF->size())) return 1E12;

  double dp =  getDeltaPhi( myJetsPF->at(thisJet).phi , getMETphi() );
  double deltaT = getTransverseMETError(thisJet);
  return dp/atan2(deltaT,getMET());

}

double getMinDeltaPhiNMET(unsigned int maxjets){
  
  double mindpn=1E12;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< myJetsPF->size(); i++) {
    
    if (isGoodJet(i)) {
      ++ngood;
      double dpn =  getDeltaPhiNMET( i );
      if (dpn<mindpn) mindpn=dpn;
      if (ngood >= maxjets) break;
    }
  }
  return mindpn;
}

double getMaxDeltaPhiNMET(unsigned int maxjets){
  
  double maxdpn=-99.;

  unsigned int ngood=0;
  //get the maximum angle between the first n jets and MET
  for (unsigned int i=0; i< myJetsPF->size(); i++) {
    
    if (isGoodJet(i)) {
      ++ngood;
      double dpn =  getDeltaPhiNMET( i );
      if (dpn>maxdpn) maxdpn=dpn;
      if (ngood >= maxjets) break;
    }
  }
  return maxdpn;
}

double getTransverseMETSignificance(unsigned int thisJet){

  if(!(thisJet<myJetsPF->size())) return -99;

  double deltaT=getTransverseMETError(thisJet);
  double dp = getDeltaPhi( myJetsPF->at(thisJet).phi , getMETphi() );
  return getMET()*sin(dp)/deltaT;
}

double getMaxTransverseMETSignificance(unsigned int maxjets){

  double maxTransverseMETSignificance=-100;

  unsigned int ngood=0;
  //get the maximum Transverse MET Significance
  for (unsigned int i=0; i< myJetsPF->size(); i++) {
    
    if (isGoodJet(i)) {
      ++ngood;
      double transverseMETSignificance = getTransverseMETSignificance(i);
      if (transverseMETSignificance>maxTransverseMETSignificance) maxTransverseMETSignificance=transverseMETSignificance;
      if (ngood >= maxjets) break;
    }
  }
  
  return maxTransverseMETSignificance;
}

double getMinTransverseMETSignificance(unsigned int maxjets){

  double minTransverseMETSignificance=1E12;

  unsigned int ngood=0;
  //get the minimum Transverse MET Significance
  for (unsigned int i=0; i< myJetsPF->size(); i++) {
    
    if (isGoodJet(i)) {
      ++ngood;
      double transverseMETSignificance = getTransverseMETSignificance(i);
      if (transverseMETSignificance<minTransverseMETSignificance) minTransverseMETSignificance=transverseMETSignificance;
      if (ngood >= maxjets) break;
    }
  }
  
  return minTransverseMETSignificance;
}


double getMinDeltaPhiMET30(unsigned int maxjets) {

  double mindp=99;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< myJetsPF->size(); i++) {
    
    if (isGoodJet30(i)) {
      ++ngood;
      double dp =  getDeltaPhi( myJetsPF->at(i).phi , getMETphi());
      if (dp<mindp) mindp=dp;
      if (ngood >= maxjets) break;
    }
  }
  
  return mindp;
}

double getMinDeltaPhiMET30_eta5(unsigned int maxjets) {

  double mindp=99;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< myJetsPF->size(); i++) {
    
    if (getJetPt(i) > 30 && fabs(myJetsPF->at(i).eta) < 5 && jetPassLooseID(i)) {
      ++ngood;
      double dp =  getDeltaPhi( myJetsPF->at(i).phi , getMETphi());
      if (dp<mindp) mindp=dp;
      if (ngood >= maxjets) break;
    }
  }
  
  return mindp;
}

double getMinDeltaPhiMET30_eta5_noId(unsigned int maxjets) {

  double mindp=99;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< myJetsPF->size(); i++) {
    
    if (getJetPt(i) > 30 && fabs(myJetsPF->at(i).eta) < 5) {
      ++ngood;
      double dp =  getDeltaPhi( myJetsPF->at(i).phi , getMETphi());
      if (dp<mindp) mindp=dp;
      if (ngood >= maxjets) break;
    }
  }
  
  return mindp;
}

double getMaxDeltaPhiMET(unsigned int maxjets) {

  double maxdp=0;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< myJetsPF->size(); i++) {
    
    if (isGoodJet(i)) {
      ++ngood;
      double dp =  getDeltaPhi( myJetsPF->at(i).phi , getMETphi());
      if (dp>maxdp) maxdp=dp;
      if (ngood >= maxjets) break;
    }
  }
  
  return maxdp;
}
double getMaxDeltaPhiMET30(unsigned int maxjets) {

  double maxdp=0;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< myJetsPF->size(); i++) {
    
    if (isGoodJet30(i)) {
      ++ngood;
      double dp =  getDeltaPhi( myJetsPF->at(i).phi , getMETphi());
      if (dp>maxdp) maxdp=dp;
      if (ngood >= maxjets) break;
    }
  }
  
  return maxdp;
}

double getMaxDeltaPhiMET30_eta5(unsigned int maxjets) {

  double maxdp=0;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< myJetsPF->size(); i++) {
    
    if (getJetPt(i) > 30 && fabs(myJetsPF->at(i).eta) < 5 && jetPassLooseID(i)) {
      ++ngood;
      double dp =  getDeltaPhi( myJetsPF->at(i).phi , getMETphi());
      if (dp>maxdp) maxdp=dp;
      if (ngood >= maxjets) break;
    }
  }
  
  return maxdp;
}
double getMaxDeltaPhiMET30_eta5_noId(unsigned int maxjets) {

  double maxdp=0;

  unsigned int ngood=0;
  //get the minimum angle between the first n jets and MET
  for (unsigned int i=0; i< myJetsPF->size(); i++) {
    
    if (getJetPt(i) > 30 && fabs(myJetsPF->at(i).eta) < 5) {
      ++ngood;
      double dp =  getDeltaPhi( myJetsPF->at(i).phi , getMETphi());
      if (dp>maxdp) maxdp=dp;
      if (ngood >= maxjets) break;
    }
  }
  
  return maxdp;
}


double getDeltaPhiMETN( unsigned int goodJetN ){
  
  //find the goodJetN-th good jet -- this is the jet deltaPhiN will be calculated for
  unsigned int ijet = 999999;
  unsigned int goodJetI=0;
  for (unsigned int i=0; i< myJetsPF->size(); i++) {
    if (isGoodJet(i)) {
      if(goodJetI == goodJetN){
	ijet = i;
	break;
      }
      goodJetI++;
    }
  }
  if(ijet == 999999) return -99;
  
  //get sum for deltaT
  double sum = 0;
  for (unsigned int i=0; i< myJetsPF->size(); i++) {
    if(i==ijet) continue; //not really needed since it adds zero to the sum
    if(isGoodJet30(i)){
      sum += pow( (getJetPx(ijet)*getJetPy(i) - getJetPy(ijet)*getJetPx(i)), 2);      
    }//is good jet
  }//i
  
  //get deltaT
  double deltaT = (0.1 / getJetPt(ijet))*sqrt(sum);
  
  //calculate deltaPhiMETN
  double dp =  getDeltaPhi( myJetsPF->at(ijet).phi , getMETphi());
  double dpN = dp / atan2(deltaT, getMET());
  
  return dpN;
}



double getMinDeltaPhiMETN(unsigned int maxjets){
  
  double mdpN=1E12;
  
  for (unsigned int i=0; i<maxjets; i++) {//could really start at i=1
    if(i>=nGoodJets()) break;
    double dpN =  getDeltaPhiMETN(i);//i is for i'th *good* jet, starting at i=0. returns -99 if bad jet.
    if (dpN>=0 && dpN<mdpN) mdpN=dpN;//checking that dpN>=0 shouldn't be necessary after break statement above, but add it anyway 
  }
  return mdpN;
}

double getMinDeltaPhiMETTaus() {
  double mindp=99;
  for (unsigned int i=0; i< myTausPF->size(); i++) {
    if (myTausPF->at(i).pt > 30) {
      double dp =  getDeltaPhi( myTausPF->at(i).phi , getMETphi());
      if (dp<mindp) mindp=dp;
    }
  }
  return mindp;
}


std::pair<float,float> getJERAdjustedMHTxy() {
  //no difference from MHT calculation -- just return both pieces
  float mhtx=0;
  float mhty=0;

  for (unsigned int i=0; i<myJetsPF->size(); i++) {
    if (isGoodJetMHT( i ) ) {

      mhtx -= getJetPx(i);
      mhty -= getJetPy(i);
    }
  }

  //this is an experiment! add taus in!
  //-->experiment successful! leave taus in.
  //FIXME hard-coded for PF
  for (unsigned int i=0; i<myTausPF->size(); i++) {
    //start by using same pT threshold as jets
    if ( myTausPF->at(i).pt > 30 && fabs(myTausPF->at(i).eta)<5 ) {

      mhtx -= myTausPF->at(i).pt * cos(myTausPF->at(i).phi);
      mhty -= myTausPF->at(i).pt * sin(myTausPF->at(i).phi);
    }
  }
  
  return make_pair(mhtx,mhty);
}

float getMHT() {
  std::pair<float,float> mht=getJERAdjustedMHTxy();
  
  return sqrt(mht.first*mht.first + mht.second*mht.second);
}

float getMHTphi() {
  std::pair<float,float> mht=getJERAdjustedMHTxy();

  return atan2(mht.second,mht.first);
}


double getPFMHTWeight() {
  double pfmhtweight = 1;

  // BNNs:
  std::vector<double> pfmet;
  pfmet.push_back( (double) getMET() );

  // Get the efficiency value from the BNNs:
  pfmhtweight = pfmht_efficiency( pfmet ); 


  return pfmhtweight;
}


//std::pair<float,float> getJESAdjustedMETxy() {
//
//  float myMET=-1;
//  float myMETphi=-99;
//
//  myMET    = myMETPF->at(0).pt;
//  myMETphi = myMETPF->at(0).phi;
//
//  float myMETx = myMET * cos(myMETphi);
//  float myMETy = myMET * sin(myMETphi);
//
//  if (theJESType_ != kJES0){
//  
//    //loop over all jets
//    for (unsigned int ijet=0; ijet<myJetsPF->size(); ++ijet) {
//      float jetUx = myJetsPF->at(ijet).uncor_pt * cos(myJetsPF->at(ijet).uncor_phi);
//      float jetUy = myJetsPF->at(ijet).uncor_pt * sin(myJetsPF->at(ijet).uncor_phi);
//
//      myMETx += jetUx;
//      myMETy += jetUy;
//      float jes_factor=1;
//      if (theJESType_ == kJESup) {
//	float unc =   loosejetJECUncPlus->at(ijet);
//	float cor = getJESExtraUnc(ijet);
//	unc = sqrt(unc*unc + cor*cor);
//	jes_factor += unc;
//      }
//      else if (theJESType_ == kJESdown) {
//	float unc =  loosejetJECUncMinus->at(ijet);
//	float cor = getJESExtraUnc(ijet);
//	jes_factor -= sqrt(unc*unc + cor*cor);
//      }
//
//      jetUx *= jes_factor;
//      jetUy *= jes_factor;
//
//      myMETx -= jetUx;
//      myMETy -= jetUy;
//    }
//
//  }
//
//  return make_pair(myMETx, myMETy);
//
//}










//double calc_theta( double eta){
//  return 2*atan(exp(-eta));
//}
float WCosHel_,topCosHel_,bestWMass_,bestTopMass_;


double calc_mNj( std::vector<unsigned int> jNi ) {

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
    sumE += myJetsPF->at( j1i ).energy; 

    sumPx += getJetPx(j1i);
    sumPy += getJetPy(j1i);
    sumPz += getJetPz(j1i);
  }

  double sumP2 = sumPx*sumPx + sumPy*sumPy + sumPz*sumPz ;

  if ( (sumE*sumE) < sumP2 ) return -2. ;

  return sqrt( sumE*sumE - sumP2 );
}

double calc_mNj( unsigned int j1i, unsigned int j2i) {
  std::vector<unsigned int> v;

  v.push_back(j1i);
  v.push_back(j2i);
  return calc_mNj(v);
}

double calc_mNj( unsigned int j1i, unsigned int j2i, unsigned int j3i) {
  std::vector<unsigned int> v;

  v.push_back(j1i);
  v.push_back(j2i);
  v.push_back(j3i);
  return calc_mNj(v);
}


//-- first two jets are from W.  third is b jet
void calcCosHel( unsigned int j1i, unsigned int j2i, unsigned int j3i) {

  //bool verb(true) ;
  bool verb(false) ;
  if ( verb ) { printf( "\n" ) ; }

  if ( j1i == j2i || j1i == j3i || j2i == j3i) {
    WCosHel_ = -99;
    topCosHel_ = -99;
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


  topCosHel_ = j3p4trf.Vect().Dot( betatopvec ) / ( betatopvec.Mag() * j3p4trf.Vect().Mag() ) ;
  if ( verb ) { printf(" calc_wcoshel: top cos hel = %6.3f\n", topCosHel_ ) ; }

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
  WCosHel_ = j1p4wrf.Vect().Dot( betawvec ) / ( betawvec.Mag() * j1p4wrf.Vect().Mag() ) ;
  if (verb) {
    printf(" calc_wcoshel: j1p4 in W rf: %6.1f,  %6.1f,  %6.1f,   %6.1f\n",
	   j1p4wrf.Px(), j1p4wrf.Py(), j1p4wrf.Pz(), j1p4wrf.E() ) ;
    printf(" calc_wcoshel: j2p4 in W rf: %6.1f,  %6.1f,  %6.1f,   %6.1f\n",
	   j2p4wrf.Px(), j2p4wrf.Py(), j2p4wrf.Pz(), j2p4wrf.E() ) ;
    TLorentzVector j1j2wrf = j1p4wrf + j2p4wrf ;
    printf(" calc_wcoshel: j1j2 mass in W rf: %7.2f\n", j1j2wrf.M() ) ;
    printf(" calc_wcoshel: w cos hel = = %6.3f\n", WCosHel_ ) ;
  }

  //   cout<<"cosHel calculator (t,W): "<<tcoshel<<" "<<wcoshel<<endl;
}


void fillWTop() {

  // cout<<" == event =="<<endl;

  double bestM2j=-9999;//, bestM2j_j1pt=0, bestM2j_j2pt=0;
  double bestM3j=-9999;//, bestM3j_j3pt=0;

  //adopting this code from Owen -- note the loop goes to the second to last jet only
  for (unsigned int j1i = 0; j1i < myJetsPF->size() -1; j1i++) {

    if ( isGoodJet30(j1i)) { //owen is using pT>30 cut

      //use exactly the same logic as Owen, to avoid bugs
      if (passBTagger(j1i) ) continue; //veto b jets

      //note how owen does the loop indexing here
      for (unsigned int j2i =j1i+1; j2i<myJetsPF->size(); j2i++) {
	if ( isGoodJet10(j2i)) { //owen is using a pT>10 cut here!

	  if (isGoodJet30(j2i) && passBTagger(j2i)) continue; //veto b jets with >30 gev

	  double m2j = calc_mNj(j1i,j2i);
	  if ( fabs(m2j- mW_) < fabs(bestM2j - mW_) ) {

	    bestM2j = m2j;
	    // bestM2j_j1pt = getLooseJetPt(j1i);
	    //bestM2j_j2pt = getLooseJetPt(j2i);

	    for ( unsigned int j3i=0; j3i<myJetsPF->size(); j3i++) {

	      if (j3i==j1i || j3i==j2i) continue;

	      if ( isGoodJet30(j3i) && passBTagger(j3i)) { //owen uses 30 GeV pT cut

		double m3j = calc_mNj(j1i,j2i,j3i);

		if ( fabs(m3j-mtop_) < fabs(bestM3j-mtop_) ) {
		  bestM3j=m3j;
		  calcCosHel(j1i,j2i,j3i); //update helicity angles
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

  bestWMass_ = bestM2j;
  bestTopMass_ = bestM3j;

}






bool passCleaning() {

  //if (theCleaningType_ == kNoCleaning) return true;
  //else if (theCleaningType_ == kMuonCleaning) {
    //these bools have now been simplified in the ntuple.
    //there are no longer extraneous _PF/not _PF copies
  if(passBadPFMuonFilter() && passInconsistentMuonFilter() ) return true;
    //}
  //else if (theCleaningType_ ==kMuonEcalCleaning) {
  //  bool passMuon = passesBadPFMuonFilter && passesInconsistentMuonPFCandidateFilter;
  //  bool passEcal = passEcalDeadCellCleaning();
  //  return passMuon && passEcal;
  //}
  //else {assert(0);}

  return false;
}




float getMT_Wlep() {

  float MT = -1.;
  float MT2;
  int nE = countEle();
  int nM = countMu();
  
  float myMET = getMET();
  float myMETphi = getMETphi();
  float myMETx = myMET * cos(myMETphi);
  float myMETy = myMET * sin(myMETphi);

  float myP=0., myPphi=0.;
  bool newMT = false;
  if(nE==1 && nM==0){
    newMT = true;
    myP = eleet1_;
    myPphi = elephi1_; 
  }
  else if(nE==0 && nM==1){
    newMT=true;
    myP = muonpt1_;
    myPphi = muonphi1_;
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



void getTransverseThrustVariables(float & thrust, float & thrustPhi, bool addMET)
{
    //Begin Thrust Calculations
    //Copy Formula from http://home.thep.lu.se/~torbjorn/pythia/lutp0613man2.pdf 

    unsigned int njets = myJetsPF->size();
     
    vector<double> goodJetPx;
    vector<double> goodJetPy;


    for (unsigned int i = 0; i<njets; i++) {
      if ( !isGoodJet(i) ) continue;
      //ngoodjets++;
      double phi = myJetsPF->at(i).phi;
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

    double phiGuess = myJetsPF->at(0).phi;
    RooRealVar phi("phi","phi",phiGuess,-1.0*TMath::Pi(),TMath::Pi());
    RooTransverseThrustVar negThrustValue("negThrustValue","negThrustValue",phi,goodJetPx,goodJetPy);
    RooMinuit rm(negThrustValue);
    rm.migrad();
    thrust = -1*negThrustValue.getVal();
    thrustPhi = (phi.getVal()>0) ? phi.getVal()-TMath::Pi() : phi.getVal()+TMath::Pi();
}

//event shape variables from Luke
void getSphericityJetMET(float & lambda1, float & lambda2, float & det,
				    const int jetmax, bool addMET) {


  TMatrixD top3Jets(2,2);
  double top3JetsScale = 0;
  
  unsigned int njets=  myJetsPF->size();
  int ngoodj=0;
  for (unsigned int i=0; i<njets ; i++) {
    if ( isGoodJet(i) ) {
      ++ngoodj;
      double phi = myJetsPF->at(i).phi;
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



bool passCut(const TString cutTag) {

  if (cutTag=="cutLumiMask" ) return passLumiMask();
  if (cutTag=="cutTrigger" ) return passHLT();
  if (cutTag=="cutUtilityTrigger" ) return ( utilityHLT_HT300()>0 || utilityHLT_HT300_CentralJet30_BTagIP()>0 );

  if (cutTag=="cutPV") return passPV();
  
  if (cutTag=="cutHT") return getHT()>350;
  
  if (cutTag=="cut3Jets") return (nGoodJets() >= 3);
  
  if (cutTag=="cutMuVeto") return passMuVeto();
  if (cutTag=="cutEleVeto") return passEleVeto();

  if (cutTag == "cutMET") {
    float mymet = getMET();
    return mymet >= 200;
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
  
  if (cutTag == "cutCleaning") return passCleaning();

  return true;
}


std::vector<TString> cutTags_;
std::map<TString, TString> cutNames_; //key is a cutTag. these should be "human readable" but should not have any spaces
std::vector<TString> ignoredCut_; //allow more than 1 ignored cut!
std::vector<TString> requiredCut_; //new feature to *turn on* a cut that is usually not required by a given cut scheme

bool setCutScheme() {

  //if (cutscheme == nCutSchemes) return false;
  //theCutScheme_ = cutscheme;

  cutNames_.clear();
  cutTags_.clear();

  cutTags_.push_back("cutInclusive");cutNames_[ cutTags_.back()] = "Inclusive";
  
  //cutTags_.push_back("cut2SUSYb"); cutNames_[ cutTags_.back()] = "==2SUSYb";
  //cutTags_.push_back("cut4SUSYb"); cutNames_[ cutTags_.back()] = "==4SUSYb";
  
  cutTags_.push_back("cutLumiMask"); cutNames_[cutTags_.back()]="LumiMask";
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
  cutTags_.push_back("cutDeltaPhiTaus"); cutNames_[cutTags_.back()]="DeltaPhiMETTaus";
  //cutTags_.push_back("cutCleaning"); cutNames_[cutTags_.back()]="TailCleaning";
  
  cutTags_.push_back("cut1b"); cutNames_[cutTags_.back()]=">=1b";
  cutTags_.push_back("cut2b"); cutNames_[cutTags_.back()]=">=2b";
  cutTags_.push_back("cut3b"); cutNames_[cutTags_.back()]=">=3b";

  cutTags_.push_back("cutCleaning"); cutNames_[cutTags_.back()]="TailCleaning";


  return true;
}

uint nBcut_;
void setBCut(unsigned int nb) {
  nBcut_=nb;
}

void setIgnoredCut(const TString cutTag) {

  bool ok=false;
  for (unsigned int i=0; i<cutTags_.size() ; i++) {
    if (cutTags_[i] == cutTag) {ok=true; break;}
  }
  if (!ok) {cout<<"Invalid cutIndex"<<endl; return;}

  ignoredCut_.push_back(cutTag);

}

void setRequiredCut(const TString cutTag) {

  bool ok=false;
  for (unsigned int i=0; i<cutTags_.size() ; i++) {
    if (cutTags_[i] == cutTag) {ok=true; break;}
  }
  if (!ok) {cout<<"Invalid cutIndex"<<endl; return;}

  requiredCut_.push_back(cutTag);

}

void resetIgnoredCut() {
  ignoredCut_.clear();
}

void resetRequiredCut() {
  requiredCut_.clear();
}




bool cutRequired(const TString cutTag) { //should put an & in here to improve performance

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

  else if (cutTag == "cutLumiMask")  cutIsRequired = true;
  else if (cutTag == "cutTrigger")  cutIsRequired = true;
  else if (cutTag == "cutPV")  cutIsRequired =  true;
  else if (cutTag == "cutHT")  cutIsRequired =  true;
  else if (cutTag == "cut3Jets")  cutIsRequired =  true;
  //else if (cutTag == "cutJetPt1")  cutIsRequired =  false;
  //else if (cutTag == "cutJetPt2")  cutIsRequired =  false;
  //else if (cutTag == "cutJetPt3")  cutIsRequired =  false;
  //roll MHT into MET
  //else if (cutTag == "cutMET")  cutIsRequired =  (theMETType_== kMET || 
  //						  theMETType_==ktcMET ||
  //					  theMETType_==kpfMET ||
  //					  theMETType_==kMHT);
  else if (cutTag == "cutMET")  cutIsRequired =  true; 
  else if (cutTag == "cutMuVeto") cutIsRequired =  true;
  else if (cutTag == "cutEleVeto") cutIsRequired =  true;
  //else if (cutTag == "cutTauVeto") cutIsRequired =  false; //not required for now
  //else if (cutTag == "cutDeltaPhi") cutIsRequired =  true;
  else if (cutTag == "cutDeltaPhiN") cutIsRequired =  true;
  else if (cutTag == "cutDeltaPhiTaus") cutIsRequired =  true;
  else if (cutTag == "cut1b") cutIsRequired =  nBcut_ >=1;
  else if (cutTag == "cut2b") cutIsRequired =  nBcut_ >=2;
  else if (cutTag == "cut3b") cutIsRequired =  nBcut_ >=3;
  else if (cutTag == "cutCleaning") cutIsRequired = true;
  //else assert(0);

  return cutIsRequired;
}


int Cut(unsigned int entry)
{
  //be careful...not every function that makes cuts uses this! e.g. ::cutflow()

  for (unsigned int i=0; i< cutTags_.size(); i++) {
    if (cutRequired( cutTags_[i] ) && !passCut( cutTags_[i]) ) return -1;
  }
  
  return 1;
  
}



double getCrossSection(TString inname){
  //from https://twiki.cern.ch/twiki/bin/view/CMS/SusyRA2BJets2011?rev=36

  const double bf = 0.32442;

  if (inname.Contains("DYJetsToLL_TuneD6T_M-10To50_7TeV-madgraph-tauola") )     return 310;
  if (inname.Contains("DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola") )         return 3048;
  if (inname.Contains("LM13_SUSY_sftsht_7TeV-pythia6") )                        return 6.899;
  if (inname.Contains("LM9_SUSY_sftsht_7TeV-pythia6_2") )                       return 7.134; //from ProductionSpring2011 twiki
  if (inname.Contains("QCD_Pt_0to5_TuneZ2_7TeV_pythia6_3") )                    return 4.844e10;
  if (inname.Contains("QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6") )               return .3321;
  if (inname.Contains("QCD_Pt_120to170_TuneZ2_7TeV_pythia6") )                  return 1.151e5;
  if (inname.Contains("QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6_3") )              return .01087;
  if (inname.Contains("QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6") )             return 2.213e10;
  if (inname.Contains("QCD_Pt_15to30_TuneZ2_7TeV_pythia6") )                    return 8.159e8;
  if (inname.Contains("QCD_Pt_170to300_TuneZ2_7TeV_pythia6") )                  return 2.426e4;
  if (inname.Contains("QCD_Pt_1800_TuneZ2_7TeV_pythia6") )                      return .0003575;
  if (inname.Contains("QCD_Pt_300to470_TuneZ2_7TeV_pythia6") )                  return 1.168e3;
  if (inname.Contains("QCD_Pt_30to50_TuneZ2_7TeV_pythia6") )                    return 5.312e7;
  if (inname.Contains("QCD_Pt_470to600_TuneZ2_7TeV_pythia6") )                  return 7.022e1;
  if (inname.Contains("QCD_Pt_50to80_TuneZ2_7TeV_pythia6") )                    return 6.359e6;
  if (inname.Contains("QCD_Pt_5to15_TuneZ2_7TeV_pythia6") )                     return 3.675e10;
  if (inname.Contains("QCD_Pt_600to800_TuneZ2_7TeV_pythia6") )                  return 1.555e1;
  if (inname.Contains("QCD_Pt_800to1000_TuneZ2_7TeV_pythia6") )                 return 1.844;
  if (inname.Contains("QCD_Pt_80to120_TuneZ2_7TeV_pythia6_3") )                 return 7.843e5;
  if (inname.Contains("TTJets_TuneD6T_7TeV-madgraph-tauola") )                  return 158; // +/- 10 +/- 15 //CMS PAS TOP-11-001
  if (inname.Contains("TToBLNu_TuneZ2_s-channel_7TeV-madgraph") )               return bf*4.6; //note BF factor
  if (inname.Contains("TToBLNu_TuneZ2_t-channel_7TeV-madgraph") )               return bf*64.6; //note BF factor
  if (inname.Contains("TToBLNu_TuneZ2_tW-channel_7TeV-madgraph") )              return 10.6;
  if (inname.Contains("WJetsToLNu_TuneZ2_7TeV-madgraph-tauola") )               return 31314;
  if (inname.Contains("WWtoAnything_TuneZ2_7TeV-pythia6-tauola") )              return 43;
  if (inname.Contains("WZtoAnything_TuneZ2_7TeV-pythia6-tauola") )              return 10.4;
  if (inname.Contains("ZinvisibleJets_7TeV-madgraph") )                         return 5760; //NNLO, taken from RA2 note //RA1 uses 5715
  if (inname.Contains("ZZtoAnything_TuneZ2_7TeV-pythia6-tauola") )              return 4.297;

  //Summer11 samples
  if (inname.Contains("qcd_tunez2_pt0to5_summer11") )                           return 4.844e10;
  if (inname.Contains("qcd_tunez2_pt1000to1400_summer11") )                     return .3321;
  if (inname.Contains("qcd_tunez2_pt120to170_summer11") )                       return 1.151e5;
  if (inname.Contains("qcd_tunez2_pt1400to1800_summer11") )                     return .01087;
  if (inname.Contains("qcd_tunez2_pt15to3000_summer11") )                       return 2.213e10;
  if (inname.Contains("qcd_tunez2_pt15to30_summer11") )                         return 8.159e8;
  if (inname.Contains("qcd_tunez2_pt170to300_summer11") )                       return 2.426e4;
  if (inname.Contains("qcd_tunez2_pt1800_summer11") )                           return .0003575;
  if (inname.Contains("qcd_tunez2_pt300to470_summer11") )                       return 1.168e3;
  if (inname.Contains("qcd_tunez2_pt30to50_summer11") )                         return 5.312e7;
  if (inname.Contains("qcd_tunez2_pt470to600_summer11") )                       return 7.022e1;
  if (inname.Contains("qcd_tunez2_pt50to80_summer11") )                         return 6.359e6;
  if (inname.Contains("qcd_tunez2_pt5to15_summer11") )                          return 3.675e10;
  if (inname.Contains("qcd_tunez2_pt600to800_summer11") )                       return 1.555e1;
  if (inname.Contains("qcd_tunez2_pt800to1000_summer11") )                      return 1.844;
  if (inname.Contains("qcd_tunez2_pt80to120_summer11") )                        return 7.843e5;
  if (inname.Contains("ttjets_tunez2_madgraph_tauola_summer11") )               return 158; // +/- 10 +/- 15 //CMS PAS TOP-11-001

  if (inname.Contains("LM9_SUSY_sftsht_7TeV-pythia6") )                         return 7.134; //from ProductionSpring2011 twiki
  if (inname.Contains("TTJets_TuneZ2_7TeV-madgraph-tauola") )                   return 158; // +/- 10 +/- 15 //CMS PAS TOP-11-001
    
  std::cout<<"Cannot find cross section for this sample!"<<std::endl;
  assert(0); 
  return -1;
}


TString getSampleNameOutputString(TString inname){

  //strategy: as much as possible, give the name that drawReducedTrees expects,
  //and return sampleName for samples that have to be 'hadd'ed afterwards anyway

  if (inname.Contains("DYJetsToLL_TuneD6T_M-10To50_7TeV-madgraph-tauola") )     return inname;
  if (inname.Contains("DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola") )         return inname;
  if (inname.Contains("LM13_SUSY_sftsht_7TeV-pythia6") )                        return "LM13";
  if (inname.Contains("LM9_SUSY_sftsht_7TeV-pythia6") )                       return "LM9";
  if (inname.Contains("QCD_Pt_0to5_TuneZ2_7TeV_pythia6") )                    return inname;
  if (inname.Contains("QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6") )                return inname;
  if (inname.Contains("QCD_Pt_120to170_TuneZ2_7TeV_pythia6") )                  return inname;
  if (inname.Contains("QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6") )              return inname;
  if (inname.Contains("QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6") )             return "PythiaPUQCDFlat";
  if (inname.Contains("QCD_Pt_15to30_TuneZ2_7TeV_pythia6") )                    return inname;
  if (inname.Contains("QCD_Pt_170to300_TuneZ2_7TeV_pythia6") )                  return inname;
  if (inname.Contains("QCD_Pt_1800_TuneZ2_7TeV_pythia6") )                      return inname;
  if (inname.Contains("QCD_Pt_300to470_TuneZ2_7TeV_pythia6") )                  return inname;
  if (inname.Contains("QCD_Pt_30to50_TuneZ2_7TeV_pythia6") )                    return inname;
  if (inname.Contains("QCD_Pt_470to600_TuneZ2_7TeV_pythia6") )                  return inname;
  if (inname.Contains("QCD_Pt_50to80_TuneZ2_7TeV_pythia6") )                    return inname;
  if (inname.Contains("QCD_Pt_5to15_TuneZ2_7TeV_pythia6") )                     return inname;
  if (inname.Contains("QCD_Pt_600to800_TuneZ2_7TeV_pythia6") )                  return inname;
  if (inname.Contains("QCD_Pt_800to1000_TuneZ2_7TeV_pythia6") )                 return inname;
  if (inname.Contains("QCD_Pt_80to120_TuneZ2_7TeV_pythia6") )                 return inname;
  if (inname.Contains("TTJets_TuneD6T_7TeV-madgraph-tauola") )                  return "TTbarJets";
  if (inname.Contains("TToBLNu_TuneZ2_s-channel_7TeV-madgraph") )               return "SingleTop-sChannel";
  if (inname.Contains("TToBLNu_TuneZ2_t-channel_7TeV-madgraph") )               return "SingleTop-tChannel";
  if (inname.Contains("TToBLNu_TuneZ2_tW-channel_7TeV-madgraph") )              return "SingleTop-tWChannel";
  if (inname.Contains("WJetsToLNu_TuneZ2_7TeV-madgraph-tauola") )               return "WJets";
  if (inname.Contains("WWtoAnything_TuneZ2_7TeV-pythia6-tauola") )              return "WW";
  if (inname.Contains("WZtoAnything_TuneZ2_7TeV-pythia6-tauola") )              return "WZ";
  if (inname.Contains("ZinvisibleJets_7TeV-madgraph") )                         return "Zinvisible";
  if (inname.Contains("ZZtoAnything_TuneZ2_7TeV-pythia6-tauola") )              return "ZZ";


  //Summer11 samples
  if (inname.Contains("qcd_tunez2_pt0to5_summer11") )                           return inname;
  if (inname.Contains("qcd_tunez2_pt1000to1400_summer11") )                     return inname;
  if (inname.Contains("qcd_tunez2_pt120to170_summer11") )                       return inname;
  if (inname.Contains("qcd_tunez2_pt1400to1800_summer11") )                     return inname;
  if (inname.Contains("qcd_tunez2_pt15to3000_summer11") )                       return inname;
  if (inname.Contains("qcd_tunez2_pt15to30_summer11") )                         return inname;
  if (inname.Contains("qcd_tunez2_pt170to300_summer11") )                       return inname;
  if (inname.Contains("qcd_tunez2_pt1800_summer11") )                           return inname;
  if (inname.Contains("qcd_tunez2_pt300to470_summer11") )                       return inname;
  if (inname.Contains("qcd_tunez2_pt30to50_summer11") )                         return inname;
  if (inname.Contains("qcd_tunez2_pt470to600_summer11") )                       return inname;
  if (inname.Contains("qcd_tunez2_pt50to80_summer11") )                         return inname;
  if (inname.Contains("qcd_tunez2_pt5to15_summer11") )                          return inname;
  if (inname.Contains("qcd_tunez2_pt600to800_summer11") )                       return inname;
  if (inname.Contains("qcd_tunez2_pt800to1000_summer11") )                      return inname;
  if (inname.Contains("qcd_tunez2_pt80to120_summer11") )                        return inname;
  if (inname.Contains("ttjets_tunez2_madgraph_tauola_summer11") )               return "TTbarJets";


  //if it isn't found, just use the full name 
  //  -- data will fall into this category, which is fine because we have to hadd the files after anyway
  return inname;

}


double getWeight(Long64_t nentries) {
  if(isRealDataInt(myEDM_isRealData)) return 1;

  double sigma = getCrossSection(sampleName_);
  double w = lumi_ * sigma / double(nentries);

  if (sampleName_.Contains("QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6")) w *= (*myGenWeight);

  return  w;
}

float getPUWeight(reweight::LumiReWeighting lumiWeights_){

  float weight;
  float sum_nvtx = 0;
  int npv = 0;
  for ( unsigned int i = 0; i<pileupsummaryinfo.size() ; i++) {
    npv = pileupsummaryinfo.at(i).addpileupinfo_getPU_NumInteractions;
    sum_nvtx += float(npv);
  }
      
  float ave_nvtx = sum_nvtx/3.;
  weight = lumiWeights_.ITweight3BX( ave_nvtx );


  //following: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupSystematicErrors
  reweight::PoissonMeanShifter PShift_;
  if(thePUuncType_ == kPUuncDown){
    PShift_ = reweight::PoissonMeanShifter(-0.6);
    weight = weight * PShift_.ShiftWeight( ave_nvtx );
  }
  else if(thePUuncType_ == kPUuncUp){
    PShift_ = reweight::PoissonMeanShifter(0.6);
    weight = weight * PShift_.ShiftWeight( ave_nvtx );
  }

  return weight;
}



TDatime* starttime_;
void startTimer() {
  starttime_ = new TDatime();
}

//void checkTimer(const Long64_t ndone, const Long64_t ntotal) {
//  if (ndone==0) return;
//  double fracdone = double(ndone)/double(ntotal);
//  
//  TDatime timenow;
//  UInt_t elapsed= timenow.Convert() - starttime_->Convert();
//  
//  std::cout<<int(100*fracdone) << " ["<<(1-fracdone)*double(elapsed)/fracdone<<" seconds remaining]"<<std::endl;
//  
//}

void stopTimer(const Long64_t ntotal) {
  TDatime stoptime; //default ctor is for current time
  double elapsed= stoptime.Convert() - starttime_->Convert();
  std::cout<<"events / time = "<<ntotal<<" / "<<elapsed<<" = "<<double(ntotal)/double(elapsed)<<" Hz"<<std::endl;
  
  delete starttime_;
  starttime_=0;
}


void cutflow(itreestream& stream, int maxevents=-1){

  std::vector<int> npass;
  std::vector<double> sumw; //sum of weights
  std::vector<double> sumw2; //sum of weights squared
  int npass_eq1b=0;
  double sumw_eq1b=0,sumw2_eq1b=0; //special kludge for eq1b case

  int nevents = stream.size();

  setMuonReq(0);
  setEleReq(0);
  setBCut(3);
  setCutScheme();
  setIgnoredCut("cutDeltaPhiTaus");

  for (unsigned int i=0 ; i<cutTags_.size(); i++) {
    //std::cout << "cutTags_["<< i<< "]=" << cutTags_.at(i) << std::endl;
    npass.push_back(0);
    sumw.push_back(0);
    sumw2.push_back(0);
  }

  cout<<"Running..."<<endl;  

  InitializeStuff();//BEN


  startTimer();
  for(int entry=0; entry < nevents; ++entry){

    if(maxevents>0 && entry>=maxevents) break;

    // Read event into memory
    stream.read(entry);
    fillObjects();

    //some output to watch while it's running
    if(entry==0){
      if(!isRealDataInt(myEDM_isRealData)){
        cout << "MC xsec: " << getCrossSection(sampleName_) << endl;
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
	if (passCut("cutEq1b"))	{
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
      	//	error = sqrt(npass_eq1b);
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



void reducedTree(TString outputpath, itreestream& stream)
{

  //open output file
  TString outfilename="reducedTree";
  outfilename += ".";
  outfilename += getCutDescriptionString();
  outfilename += ".";

  outfilename += getSampleNameOutputString(sampleName_);
  outfilename+=".root";
  if (outputpath[outputpath.Length()-1] != '/') outputpath += "/";
  outfilename.Prepend(outputpath);
  TFile fout(outfilename,"RECREATE");

  //TFile fout("testing.root","RECREATE");

  // define the TTree
  TTree reducedTree("reducedTree","tree with minimal cuts");
  reducedTree.SetMaxTreeSize(4900000000LL); //hopefully this cures crash at 1.9GB
  //if more variables keep getting added, I will have to split the reduced trees into pieces
  //hopefully that can be avoided

  unsigned int seed = 4357;

  if (sampleName_.Contains("DYJetsToLL_TuneD6T_M-10To50_7TeV-madgraph-tauola") )     seed = 4357;
  if (sampleName_.Contains("DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola") )         seed = 4358;
  if (sampleName_.Contains("LM13_SUSY_sftsht_7TeV-pythia6") )                        seed = 4359;
  if (sampleName_.Contains("LM9_SUSY_sftsht_7TeV-pythia6_2") )                       seed = 4360;
  if (sampleName_.Contains("QCD_Pt_0to5_TuneZ2_7TeV_pythia6_3") )                    seed = 4361;
  if (sampleName_.Contains("QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6") )                seed = 4362;
  if (sampleName_.Contains("QCD_Pt_120to170_TuneZ2_7TeV_pythia6") )                  seed = 4363;
  if (sampleName_.Contains("QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6_3") )              seed = 4364;
  if (sampleName_.Contains("QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6") )             seed = 4365;
  if (sampleName_.Contains("QCD_Pt_15to30_TuneZ2_7TeV_pythia6") )                    seed = 4366;
  if (sampleName_.Contains("QCD_Pt_170to300_TuneZ2_7TeV_pythia6") )                  seed = 4367;
  if (sampleName_.Contains("QCD_Pt_1800_TuneZ2_7TeV_pythia6") )                      seed = 4368;
  if (sampleName_.Contains("QCD_Pt_300to470_TuneZ2_7TeV_pythia6") )                  seed = 4369;
  if (sampleName_.Contains("QCD_Pt_30to50_TuneZ2_7TeV_pythia6") )                    seed = 4370;
  if (sampleName_.Contains("QCD_Pt_470to600_TuneZ2_7TeV_pythia6") )                  seed = 4371;
  if (sampleName_.Contains("QCD_Pt_50to80_TuneZ2_7TeV_pythia6") )                    seed = 4372;
  if (sampleName_.Contains("QCD_Pt_5to15_TuneZ2_7TeV_pythia6") )                     seed = 4373;
  if (sampleName_.Contains("QCD_Pt_600to800_TuneZ2_7TeV_pythia6") )                  seed = 4374;
  if (sampleName_.Contains("QCD_Pt_800to1000_TuneZ2_7TeV_pythia6") )                 seed = 4375;
  if (sampleName_.Contains("QCD_Pt_80to120_TuneZ2_7TeV_pythia6_3") )                 seed = 4376;
  if (sampleName_.Contains("TTJets_TuneD6T_7TeV-madgraph-tauola") )                  seed = 4377;
  if (sampleName_.Contains("TToBLNu_TuneZ2_s-channel_7TeV-madgraph") )               seed = 4378;
  if (sampleName_.Contains("TToBLNu_TuneZ2_t-channel_7TeV-madgraph") )               seed = 4379;
  if (sampleName_.Contains("TToBLNu_TuneZ2_tW-channel_7TeV-madgraph") )              seed = 4380;
  if (sampleName_.Contains("WJetsToLNu_TuneZ2_7TeV-madgraph-tauola") )               seed = 4381;
  if (sampleName_.Contains("WWtoAnything_TuneZ2_7TeV-pythia6-tauola") )              seed = 4382;
  if (sampleName_.Contains("WZtoAnything_TuneZ2_7TeV-pythia6-tauola") )              seed = 4383;
  if (sampleName_.Contains("ZinvisibleJets_7TeV-madgraph") )                         seed = 4384;
  if (sampleName_.Contains("ZZtoAnything_TuneZ2_7TeV-pythia6-tauola") )              seed = 4385;  

  if (sampleName_.Contains("qcd_tunez2_pt0to5_summer11") )                           seed = 4384;
  if (sampleName_.Contains("qcd_tunez2_pt1000to1400_summer11") )                     seed = 4385;
  if (sampleName_.Contains("qcd_tunez2_pt120to170_summer11") )                       seed = 4386;
  if (sampleName_.Contains("qcd_tunez2_pt1400to1800_summer11") )                     seed = 4387;
  if (sampleName_.Contains("qcd_tunez2_pt15to3000_summer11") )                       seed = 4388;
  if (sampleName_.Contains("qcd_tunez2_pt15to30_summer11") )                         seed = 4389;
  if (sampleName_.Contains("qcd_tunez2_pt170to300_summer11") )                       seed = 4390;
  if (sampleName_.Contains("qcd_tunez2_pt1800_summer11") )                           seed = 4391;
  if (sampleName_.Contains("qcd_tunez2_pt300to470_summer11") )                       seed = 4392;
  if (sampleName_.Contains("qcd_tunez2_pt30to50_summer11") )                         seed = 4393;
  if (sampleName_.Contains("qcd_tunez2_pt470to600_summer11") )                       seed = 4394;
  if (sampleName_.Contains("qcd_tunez2_pt50to80_summer11") )                         seed = 4395;
  if (sampleName_.Contains("qcd_tunez2_pt5to15_summer11") )                          seed = 4396;
  if (sampleName_.Contains("qcd_tunez2_pt600to800_summer11") )                       seed = 4397;
  if (sampleName_.Contains("qcd_tunez2_pt800to1000_summer11") )                      seed = 4398;
  if (sampleName_.Contains("qcd_tunez2_pt80to120_summer11") )                        seed = 4399;
  if (sampleName_.Contains("ttjets_tunez2_madgraph_tauola_summer11") )               seed = 4400;

  if (sampleName_.Contains("ht_run2011a_423_jul1.txt") )                             seed = 4401;
  if (sampleName_.Contains("ht_run2011a_423_jul6.txt") )                             seed = 4402;
  if (sampleName_.Contains("ht_run2011a_423_jun10.txt") )                            seed = 4403;
  if (sampleName_.Contains("ht_run2011a_423_jun17.txt") )                            seed = 4404;
  if (sampleName_.Contains("ht_run2011a_423_jun24.txt") )                            seed = 4405;
  if (sampleName_.Contains("ht_run2011a_423_jun3.txt") )                             seed = 4406;
  if (sampleName_.Contains("ht_run2011a_423_may10rereco.txt") )                      seed = 4407;
  if (sampleName_.Contains("ht_run2011a_423_may10rereco_added.txt") )                seed = 4408;
  if (sampleName_.Contains("ht_run2011a_423_may10rereco_added_b.txt") )              seed = 4409;
  if (sampleName_.Contains("ht_run2011a_423_may24.txt") )                            seed = 4410;
  if (sampleName_.Contains("ht_run2011a_423_may27.txt") )                            seed = 4411;
										     
  TRandom3 random(seed);							     
										     
  //we're making an ntuple, so size matters -- use float not double		     
  double weight; //one exception to the float rule				     
  float btagIPweight, pfmhtweight;
  float PUweight;
  float hltHTeff;
  ULong64_t lumiSection, eventNumber, runNumber;
  float METsig;
  float HT, MHT, MET, METphi, minDeltaPhi, minDeltaPhiAll, minDeltaPhiAll30,minDeltaPhi30_eta5_noIdAll;
  float correctedMET, correctedMETphi;
  float maxDeltaPhi, maxDeltaPhiAll, maxDeltaPhiAll30, maxDeltaPhi30_eta5_noIdAll;
  float sumDeltaPhi, diffDeltaPhi;
  float minDeltaPhiN, deltaPhiN1, deltaPhiN2, deltaPhiN3;
  float minDeltaPhiMetTau;

  bool cutHT,cutPV,cutTrigger;
  bool cut3Jets,cutEleVeto,cutMuVeto,cutMET,cutDeltaPhi, cutCleaning;

  bool passBadPFMuon, passInconsistentMuon;

  bool isRealData;
  UInt_t pass_utilityHLT_HT300;
  UInt_t prescale_utilityHLT_HT300;
  UInt_t pass_utilityHLT_HT300_CentralJet30_BTagIP;
  UInt_t prescale_utilityHLT_HT300_CentralJet30_BTagIP;

  int njets, nbjets, nElectrons, nMuons;
  int nbjetsSSVM,nbjetsTCHET,nbjetsSSVHPT,nbjetsTCHPT,nbjetsTCHPM;

  float jetpt1,jetphi1, jeteta1, jetenergy1, bjetpt1, bjetphi1, bjeteta1, bjetenergy1;
  float jetpt2,jetphi2, jeteta2, jetenergy2, bjetpt2, bjetphi2, bjeteta2, bjetenergy2;
  float jetpt3,jetphi3, jeteta3, jetenergy3, bjetpt3, bjetphi3, bjeteta3, bjetenergy3;
  float eleet1;
  float muonpt1;

  float MT_Wlep;

  int nGoodPV;

  int SUSY_nb;

  //new variables from Luke
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

  float transverseThrust,transverseThrustPhi;
  float transverseThrustWithMET,transverseThrustWithMETPhi;

  float minDeltaPhiN_Luke, maxDeltaPhiN_Luke, deltaPhiN1_Luke, deltaPhiN2_Luke, deltaPhiN3_Luke;
  float minTransverseMETSignificance, maxTransverseMETSignificance, transverseMETSignificance1, transverseMETSignificance2, transverseMETSignificance3;

  int nLostJet;
  int njets_lostJet, nbjets_lostJet;
  float minDeltaPhiN_lostJet, deltaPhiN1_lostJet, deltaPhiN2_lostJet, deltaPhiN3_lostJet;
  float minDeltaPhiN_Luke_lostJet, maxDeltaPhiN_Luke_lostJet, deltaPhiN1_Luke_lostJet, deltaPhiN2_Luke_lostJet, deltaPhiN3_Luke_lostJet;
  float minTransverseMETSignificance_lostJet, maxTransverseMETSignificance_lostJet, transverseMETSignificance1_lostJet, transverseMETSignificance2_lostJet, transverseMETSignificance3_lostJet;


  //initialize PU things
  std::vector< float > TrueDist2011;
  std::vector< float > MCDist2011;
  for( int i=0; i<25; ++i) {
    TrueDist2011.push_back(TrueDist2011_f[i]);
    MCDist2011.push_back(PoissonIntDist_f[i]);
  }  
  reweight::LumiReWeighting LumiWeights_ = reweight::LumiReWeighting( MCDist2011, TrueDist2011 );


  reducedTree.Branch("weight",&weight,"weight/D");
  reducedTree.Branch("runNumber",&runNumber,"runNumber/l");
  reducedTree.Branch("lumiSection",&lumiSection,"lumiSection/l");
  reducedTree.Branch("eventNumber",&eventNumber,"eventNumber/l");

  reducedTree.Branch("btagIPweight",&btagIPweight,"btagIPweight/F");
  reducedTree.Branch("pfmhtweight",&pfmhtweight,"pfmhtweight/F");
  reducedTree.Branch("PUweight",&PUweight,"PUweight/F");
  reducedTree.Branch("hltHTeff",&hltHTeff,"hltHTeff/F");

  reducedTree.Branch("cutHT",&cutHT,"cutHT/O");
  reducedTree.Branch("cutPV",&cutPV,"cutPV/O");
  reducedTree.Branch("cutTrigger",&cutTrigger,"cutTrigger/O");
  reducedTree.Branch("cut3Jets",&cut3Jets,"cut3Jets/O");
  reducedTree.Branch("cutEleVeto",&cutEleVeto,"cutEleVeto/O");
  reducedTree.Branch("cutMuVeto",&cutMuVeto,"cutMuVeto/O");
  reducedTree.Branch("cutMET",&cutMET,"cutMET/O");
  reducedTree.Branch("cutDeltaPhi",&cutDeltaPhi,"cutDeltaPhi/O");
  reducedTree.Branch("cutCleaning",&cutCleaning,"cutCleaning/O");

  reducedTree.Branch("nGoodPV",&nGoodPV,"nGoodPV/I");

  reducedTree.Branch("SUSY_nb",&SUSY_nb,"SUSY_nb/I");

  reducedTree.Branch("passBadPFMuon",&passBadPFMuon,"passBadPFMuon/O");
  reducedTree.Branch("passInconsistentMuon",&passInconsistentMuon,"passInconsistentMuon/O");

  reducedTree.Branch("njets",&njets,"njets/I");
  reducedTree.Branch("nbjets",&nbjets,"nbjets/I");
  reducedTree.Branch("nElectrons",&nElectrons,"nElectrons/I");
  reducedTree.Branch("nMuons",&nMuons,"nMuons/I");

  reducedTree.Branch("nbjetsSSVM",&nbjetsSSVM,"nbjetsSSVM/I");
  reducedTree.Branch("nbjetsSSVHPT",&nbjetsSSVHPT,"nbjetsSSVHPT/I");
  reducedTree.Branch("nbjetsTCHPT",&nbjetsTCHPT,"nbjetsTCHPT/I");
  reducedTree.Branch("nbjetsTCHET",&nbjetsTCHET,"nbjetsTCHET/I");
  reducedTree.Branch("nbjetsTCHPM",&nbjetsTCHPM,"nbjetsTCHPM/I");

  reducedTree.Branch("isRealData",&isRealData,"isRealData/O");
  reducedTree.Branch("pass_utilityHLT_HT300",&pass_utilityHLT_HT300,"pass_utilityHLT_HT300/i");
  reducedTree.Branch("prescale_utilityHLT_HT300", &prescale_utilityHLT_HT300, "prescale_utilityHLT_HT300/i");
  reducedTree.Branch("pass_utilityHLT_HT300_CentralJet30_BTagIP",&pass_utilityHLT_HT300_CentralJet30_BTagIP,"pass_utilityHLT_HT300_CentralJet30_BTagIP/i");
  reducedTree.Branch("prescale_utilityHLT_HT300_CentralJet30_BTagIP", &prescale_utilityHLT_HT300_CentralJet30_BTagIP, "prescale_utilityHLT_HT300_CentralJet30_BTagIP/i");  

  reducedTree.Branch("HT",&HT,"HT/F");
  reducedTree.Branch("MET",&MET,"MET/F");
  reducedTree.Branch("METsig",&METsig,"METsig/F");
  reducedTree.Branch("METphi",&METphi,"METphi/F");
  reducedTree.Branch("MHT",&MHT,"MHT/F");

  reducedTree.Branch("correctedMET",&correctedMET,"correctedMET/F");
  reducedTree.Branch("correctedMETphi",&correctedMETphi,"correctedMETphi/F");

  reducedTree.Branch("bestWMass",&bestWMass_,"bestWMass/F");
  reducedTree.Branch("bestTopMass",&bestTopMass_,"bestTopMass/F");
  reducedTree.Branch("topCosHel",&topCosHel_,"topCosHel/F");
  reducedTree.Branch("WCosHel",&WCosHel_,"WCosHel/F");
  reducedTree.Branch("MT_Wlep",&MT_Wlep, "MT_Wlep/F");

  reducedTree.Branch("minDeltaPhi",&minDeltaPhi,"minDeltaPhi/F");
  reducedTree.Branch("minDeltaPhiAll",&minDeltaPhiAll,"minDeltaPhiAll/F");
  reducedTree.Branch("minDeltaPhiAll30",&minDeltaPhiAll30,"minDeltaPhiAll30/F");
  reducedTree.Branch("minDeltaPhi30_eta5_noIdAll",&minDeltaPhi30_eta5_noIdAll,"minDeltaPhi30_eta5_noIdAll/F");

  reducedTree.Branch("minDeltaPhiMetTau",&minDeltaPhiMetTau,"minDeltaPhiMetTau/F");

  reducedTree.Branch("maxDeltaPhi",&maxDeltaPhi,"maxDeltaPhi/F");
  reducedTree.Branch("maxDeltaPhiAll",&maxDeltaPhiAll,"maxDeltaPhiAll/F");
  reducedTree.Branch("maxDeltaPhiAll30",&maxDeltaPhiAll30,"maxDeltaPhiAll30/F");
  reducedTree.Branch("maxDeltaPhi30_eta5_noIdAll",&maxDeltaPhi30_eta5_noIdAll,"maxDeltaPhi30_eta5_noIdAll/F");

  reducedTree.Branch("sumDeltaPhi",&sumDeltaPhi,"sumDeltaPhi/F");
  reducedTree.Branch("diffDeltaPhi",&diffDeltaPhi,"diffDeltaPhi/F");

  reducedTree.Branch("minDeltaPhiN", &minDeltaPhiN, "minDeltaPhiN/F");
  reducedTree.Branch("deltaPhiN1", &deltaPhiN1, "deltaPhiN1/F");
  reducedTree.Branch("deltaPhiN2", &deltaPhiN2, "deltaPhiN2/F");
  reducedTree.Branch("deltaPhiN3", &deltaPhiN3, "deltaPhiN3/F");

  reducedTree.Branch("minDeltaPhiN_lostJet", &minDeltaPhiN_lostJet, "minDeltaPhiN_lostJet/F");
  reducedTree.Branch("deltaPhiN1_lostJet", &deltaPhiN1_lostJet, "deltaPhiN1_lostJet/F");
  reducedTree.Branch("deltaPhiN2_lostJet", &deltaPhiN2_lostJet, "deltaPhiN2_lostJet/F");
  reducedTree.Branch("deltaPhiN3_lostJet", &deltaPhiN3_lostJet, "deltaPhiN3_lostJet/F");

  reducedTree.Branch("jetpt1",&jetpt1,"jetpt1/F");
  reducedTree.Branch("jeteta1",&jeteta1,"jeteta1/F");
  reducedTree.Branch("jetphi1",&jetphi1,"jetphi1/F");
  reducedTree.Branch("jetenergy1",&jetenergy1,"jetenergy1/F");

  reducedTree.Branch("jetpt2",&jetpt2,"jetpt2/F");
  reducedTree.Branch("jeteta2",&jeteta2,"jeteta2/F");
  reducedTree.Branch("jetphi2",&jetphi2,"jetphi2/F");
  reducedTree.Branch("jetenergy2",&jetenergy2,"jetenergy2/F");

  reducedTree.Branch("jetpt3",&jetpt3,"jetpt3/F");
  reducedTree.Branch("jeteta3",&jeteta3,"jeteta3/F");
  reducedTree.Branch("jetphi3",&jetphi3,"jetphi3/F");
  reducedTree.Branch("jetenergy3",&jetenergy3,"jetenergy3/F");

  reducedTree.Branch("bjetpt1",&bjetpt1,"bjetpt1/F");
  reducedTree.Branch("bjeteta1",&bjeteta1,"bjeteta1/F");
  reducedTree.Branch("bjetphi1",&bjetphi1,"bjetphi1/F");
  reducedTree.Branch("bjetenergy1",&bjetenergy1,"bjetenergy1/F");  

  reducedTree.Branch("bjetpt2",&bjetpt2,"bjetpt2/F");
  reducedTree.Branch("bjeteta2",&bjeteta2,"bjeteta2/F");
  reducedTree.Branch("bjetphi2",&bjetphi2,"bjetphi2/F");
  reducedTree.Branch("bjetenergy2",&bjetenergy2,"bjetenergy2/F");  
  
  reducedTree.Branch("bjetpt3",&bjetpt3,"bjetpt3/F");
  reducedTree.Branch("bjeteta3",&bjeteta3,"bjeteta3/F");
  reducedTree.Branch("bjetphi3",&bjetphi3,"bjetphi3/F");
  reducedTree.Branch("bjetenergy3",&bjetenergy3,"bjetenergy3/F");  

  reducedTree.Branch("eleet1",&eleet1,"eleet1/F");
  reducedTree.Branch("muonpt1",&muonpt1,"muonpt1/F");

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

  reducedTree.Branch("transverseThrust",&transverseThrust,"transverseThrust/F");
  reducedTree.Branch("transverseThrustPhi",&transverseThrustPhi,"transverseThrustPhi/F");

  reducedTree.Branch("transverseThrustWithMET",&transverseThrustWithMET,"transverseThrustWithMET/F");
  reducedTree.Branch("transverseThrustWithMETPhi",&transverseThrustWithMETPhi,"transverseThrustWithMETPhi/F");

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

  int nevents = stream.size();

  InitializeStuff();//BEN

  startTimer();
  for(int entry=0; entry < nevents; ++entry){
    // Read event into memory
    stream.read(entry);
    fillObjects();

    //some output to watch while it's running
    if(entry==0){
      if(!isRealDataInt(myEDM_isRealData)){
	cout << "MC xsec: " << getCrossSection(sampleName_) << endl;
      }
      else{
	cout << "This is data!"<< endl;
      }
      cout << "weight: "<< getWeight(nevents) << endl;
    }
    if(entry%100000==0) cout << "  entry: " << entry << ", percent done=" << (int)(entry/(double)nevents*100.)<<  endl;
      

    if ( passCut("cutLumiMask") && (passCut("cutTrigger") || passCut("cutUtilityTrigger")) && passCut("cutHT") ) {
    
           
      //if (entry%1000000==0) checkTimer(entry,nevents);
      weight = getWeight(nevents);
      
      // cast these as long ints, with rounding, assumes they are positive to begin with 
      runNumber = (ULong64_t)((*myEDM_run)+0.5);
      lumiSection = (ULong64_t)((*myEDM_luminosityBlock)+0.5);
      eventNumber = (ULong64_t)((*myEDM_event)+0.5);
            
      HT=getHT();
      hltHTeff = getHLTHTeff(HT);
      PUweight = getPUWeight(LumiWeights_);
      btagIPweight = getBTagIPWeight();
      pfmhtweight = getPFMHTWeight();

      cutHT = true; 
      //cutHT = passCut("cutHT");

      cutTrigger = passCut("cutTrigger");
      cutPV = passCut("cutPV");
      cut3Jets = passCut("cut3Jets");
      cutEleVeto = passCut("cutEleVeto");
      cutMuVeto = passCut("cutMuVeto");
      cutMET = passCut("cutMET");
      cutDeltaPhi = passCut("cutDeltaPhi");
      cutCleaning = passCut("cutCleaning");
      
      passBadPFMuon = passBadPFMuonFilter();
      passInconsistentMuon = passInconsistentMuonFilter();
      
      isRealData = isRealDataInt(myEDM_isRealData);
      pass_utilityHLT_HT300 = utilityHLT_HT300();
      prescale_utilityHLT_HT300 = 0;
      pass_utilityHLT_HT300_CentralJet30_BTagIP = utilityHLT_HT300_CentralJet30_BTagIP();
      prescale_utilityHLT_HT300_CentralJet30_BTagIP = 0;

      nGoodPV = countGoodPV();

      //SUSY_nb = getSUSYnb();

      njets = nGoodJets();
      nbjets = nGoodBJets();

      //  int nbjetsSSVM,nbjetsTCHET,nbjetsSSVHPT,nbjetsTCHPT,nbjetsTCHPM;
      nbjetsSSVM = nGoodBJets( kSSVM);
      nbjetsSSVHPT = nGoodBJets( kSSVHPT);
      nbjetsTCHET = nGoodBJets( kTCHET);
      nbjetsTCHPT = nGoodBJets( kTCHPT);
      nbjetsTCHPM = nGoodBJets( kTCHPM);

      nElectrons = countEle();
      nMuons = countMu();
      MET=getMET();
      METsig = myMETPF->at(0).mEtSig; //FIXME hard coded for PF
      MHT=getMHT();
      METphi = getMETphi();

      getCorrectedMET(correctedMET, correctedMETphi);

      minDeltaPhi = getMinDeltaPhiMET(3);
      minDeltaPhiAll = getMinDeltaPhiMET(99);
      minDeltaPhiAll30 = getMinDeltaPhiMET30(99);
      minDeltaPhi30_eta5_noIdAll = getMinDeltaPhiMET30_eta5_noId(99);
      maxDeltaPhi = getMaxDeltaPhiMET(3);
      maxDeltaPhiAll = getMaxDeltaPhiMET(99);
      maxDeltaPhiAll30 = getMaxDeltaPhiMET30(99);
      maxDeltaPhi30_eta5_noIdAll = getMaxDeltaPhiMET30_eta5_noId(99);
      
      minDeltaPhiMetTau = getMinDeltaPhiMETTaus();

      sumDeltaPhi = maxDeltaPhi + minDeltaPhi;
      diffDeltaPhi = maxDeltaPhi - minDeltaPhi;
      
      minDeltaPhiN = getMinDeltaPhiMETN(3);
      deltaPhiN1 = getDeltaPhiMETN(0);
      deltaPhiN2 = getDeltaPhiMETN(1);
      deltaPhiN3 = getDeltaPhiMETN(2);

      minDeltaPhiN_Luke = getMinDeltaPhiNMET(3);
      maxDeltaPhiN_Luke = getMaxDeltaPhiNMET(3);
      deltaPhiN1_Luke = getDeltaPhiNMET(0);
      deltaPhiN2_Luke = getDeltaPhiNMET(1);
      deltaPhiN3_Luke = getDeltaPhiNMET(2);

      minTransverseMETSignificance = getMinTransverseMETSignificance(3);
      maxTransverseMETSignificance = getMaxTransverseMETSignificance(3);
      transverseMETSignificance1 = getTransverseMETSignificance(0);
      transverseMETSignificance2 = getTransverseMETSignificance(1);
      transverseMETSignificance3 = getTransverseMETSignificance(2);

      MT_Wlep = getMT_Wlep();
      
      jetpt1 = jetPtOfN(1);
      jetphi1 = jetPhiOfN(1);
      jeteta1 = jetEtaOfN(1);
      jetenergy1 = jetEnergyOfN(1);

      jetpt2 = jetPtOfN(2);
      jetphi2 = jetPhiOfN(2);
      jeteta2 = jetEtaOfN(2);
      jetenergy2 = jetEnergyOfN(2);
      
      jetpt3 = jetPtOfN(3);
      jetphi3 = jetPhiOfN(3);
      jeteta3 = jetEtaOfN(3);
      jetenergy3 = jetEnergyOfN(3);      

      bjetpt1 = bjetPtOfN(1);
      bjetphi1 = bjetPhiOfN(1);
      bjeteta1 = bjetEtaOfN(1);
      bjetenergy1 = bjetEnergyOfN(1); 
      
      bjetpt2 = bjetPtOfN(2);
      bjetphi2 = bjetPhiOfN(2);
      bjeteta2 = bjetEtaOfN(2);
      bjetenergy2 = bjetEnergyOfN(2); 
      
      bjetpt3 = bjetPtOfN(3);
      bjetphi3 = bjetPhiOfN(3);
      bjeteta3 = bjetEtaOfN(3);
      bjetenergy3 = bjetEnergyOfN(3);       

      if (nElectrons>=1) {
	eleet1 = eleet1_; //this is filled when i call countEle(). a hack to be cleaned up later
      }
      else {eleet1=-1;}
      if (nMuons>=1) {
	muonpt1 = muonpt1_; //this is filled when i call countMu(). a hack to be cleaned up later
      }
      else {muonpt1=-1;}
      
      
      fillWTop(); //fill W,top masses and helicity angles
      
      
      //fill new variables from Luke
      getSphericityJetMET(lambda1_allJets,lambda2_allJets,determinant_allJets,99,false);
      getSphericityJetMET(lambda1_allJetsPlusMET,lambda2_allJetsPlusMET,determinant_allJetsPlusMET,99,true);
      getSphericityJetMET(lambda1_topThreeJets,lambda2_topThreeJets,determinant_topThreeJets,3,false);
      getSphericityJetMET(lambda1_topThreeJetsPlusMET,lambda2_topThreeJetsPlusMET,determinant_topThreeJetsPlusMET,3,true);
      
      //Uncomment next two lines to do thrust calculations
      //getTransverseThrustVariables(transverseThrust, transverseThrustPhi, false);
      //getTransverseThrustVariables(transverseThrustWithMET, transverseThrustWithMETPhi, true);

      changeVariables(&random,0.05,nLostJet);

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

      reducedTree.Fill();
      
    }
  }

  stopTimer(nevents);
  fout.Write();
  fout.Close();
  
  return;
}


void sampleAnalyzer(itreestream& stream){

  int nevents = stream.size();

  TFile fout("histos.root","RECREATE");
  //TH1F * h_taupt = new TH1F("h_taupt","tau pt",200,0,200);

  InitializeStuff();//BEN

  //setMuonReq(0);
  //setEleReq(0);
  setCutScheme();
  //setBCut(2);

  setIgnoredCut("cutTrigger");
  setIgnoredCut("cutPV");
  setIgnoredCut("cutHT");
  setIgnoredCut("cut3Jets");
  setIgnoredCut("cutEleVeto");
  setIgnoredCut("cutMuVeto");
  setIgnoredCut("cutMET");
  setIgnoredCut("cutDeltaPhiN");
  setIgnoredCut("cutDeltaPhiTaus");
  setBCut(0);

  int count = 0;
  startTimer();
  for(int entry=0; entry < nevents; ++entry){
    // Read event into memory
    stream.read(entry);
    fillObjects();
    if(entry%100000==0) cout << "  entry: " << entry << ", percent done=" << (int)(entry/(double)nevents*100.)<<  endl;

    if (Cut(entry) < 0) continue;
    count++;

    //std::cout << "genweight = " << (*myGenWeight) << std::endl;
    //std::cout << "btagIP weight = " << getBTagIPWeight() << std::endl;
    //std::cout << "PFMHT weight = " << getPFMHTWeight() << std::endl;

  }
  stopTimer(nevents);

  fout.Write();
  fout.Close();

  return;
}
