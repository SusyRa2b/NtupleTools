//
// Cornell University
//


#include "BasicLoopCU.h"
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TMatrixT.h"
#include "TMatrixDEigen.h"
#include <cassert>

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

enum CutScheme {kBaseline2010};
CutScheme theCutScheme_;
enum METType {kMHT=0, kMET, ktcMET, kpfMET};
METType theMETType_;
enum METRange {kMedium=0, kHigh, kWide, kMedhigh, kLSB};
METRange theMETRange_;
enum jetType {kCalo=0, kPF, kJPT};
jetType theJetType_;
enum leptonType {kNormal=0, kPFLeptons, kPFLeptonsRA2};
leptonType theLeptonType_;
enum dpType {kDeltaPhi=0, kminDP, kMPT, kDPSync1, kminDPinv, kminDPAll30};
dpType theDPType_;
enum JESType {kJES0=0,kJESup,kJESdown};
JESType theJESType_;
enum JERType {kJER0=0,kJERbias,kJERup,kJERbias6};
JERType theJERType_;
enum METuncType {kMETunc0=0,kMETuncDown,kMETuncUp};
METuncType theMETuncType_;
enum BTagEffType {kBTagEff0=0,kBTagEffup,kBTagEffdown};
BTagEffType theBTagEffType_;
enum tailCleaningType {kNoCleaning=0, kMuonCleaning, kMuonEcalCleaning};
tailCleaningType theCleaningType_;
enum flavorHistoryType {kNoFlvHist=0, kFlvHist};
flavorHistoryType theFlavorHistoryType_;
enum BTaggerType {kSSVM=0, kTCHET};
BTaggerType theBTaggerType_;

//-Define some pointers so we can use more intuitive names.
//-This also isolates some of the dependency on the configuration of the ntuple.
//--These should be checked each time we make new ntuples.
std::vector<jet2_s> * myJetsPF; 
std::vector<electron1_s> * myElectronsPF;
std::vector<muon1_s> * myMuonsPF;
std::vector<tau1_s> * myTausPF;
std::vector<met1_s> * myMETPF;
std::vector<vertex_s> * myVertex;

void InitializeStuff(){       //If we resurrect the class structure, this will go into the constructor.
  myJetsPF = &jet2;            //selectedPatJetsPF
  myElectronsPF = &electron1;  //selectedPatElectronsPF
  myMuonsPF = &muon1;          //selectedPatMuonsPF
  myTausPF = &tau1;            //selectedPatTausPF
  myMETPF = &met1;             //patMETsPF
  myVertex = &vertex;         //offlinePrimaryVertices
}

TString getCutDescriptionString(){
  
  TString cut = "";
  if(theBTaggerType_==kTCHET){
    cut+="TCHET";
  }  
  return cut;
}

void setSampleName_(TString name){
  sampleName_ = name;
}

void setBTaggerType(BTaggerType btaggertype){
  theBTaggerType_ = btaggertype;
}

bool passHLT() { 

  long runnumber = (long)(edmevent_run +0.5); //just in case some crazy thing happens to edmevent_run being saved as a double

  //RA2b - 2011 Triggers
  bool passTrig = false;

  if(edmevent_isRealData){
    //edmtriggerresults is a double - possible values are 0 (failed trigger), 1 (passed trigger), -9999 (trig result not available)
    if(runnumber >= 160431 && runnumber < 161205) passTrig = (edmtriggerresults_HLT_HT260_MHT60_v2 > 0);
    else if (runnumber >= 161205 && runnumber < 163269) passTrig = (edmtriggerresults_HLT_HT250_MHT60_v2 > 0);
    else if (runnumber >= 163269 && runnumber < 164924) passTrig = (edmtriggerresults_HLT_HT250_MHT60_v3 > 0);
    else if (runnumber >= 164924 && runnumber < 165922) passTrig = edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v2;
    else if (runnumber >= 165922 && runnumber < 166301) passTrig = edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3;
    //else if (runnumber >= 166301 && runnumber < 166374 ) passTrig = edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v4; 
    else if (runnumber >= 166374 && runnumber < 167078) passTrig = edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3;
    //else if (runnumber >= 167078) passTrig = edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v5; 

  }
  else passTrig = true;   //use no trigger for MC   

  return passTrig;
}

unsigned int utilityHLT_HT300(){
  //for data, this function will return 0 if the utility trigger was not passed
  //and the prescale value if it was passed. For MC, it returns 1.
  
  if(edmevent_isRealData) return 1;
 
  unsigned int myPrescale = 0;
  if(edmtriggerresults_HLT_HT300_v1>0) myPrescale = edmtriggerresults_HLT_HT300_v1_prs;
  if(edmtriggerresults_HLT_HT300_v2>0) myPrescale = edmtriggerresults_HLT_HT300_v2_prs;
  if(edmtriggerresults_HLT_HT300_v3>0) myPrescale = edmtriggerresults_HLT_HT300_v3_prs;
  if(edmtriggerresults_HLT_HT300_v4>0) myPrescale = edmtriggerresults_HLT_HT300_v4_prs;
  if(edmtriggerresults_HLT_HT300_v5>0) myPrescale = edmtriggerresults_HLT_HT300_v5_prs;
  return myPrescale;
}

unsigned int utilityHLT_HT300_CentralJet30_BTagIP(){
  //this function will return 0 if the utility trigger was not passed
  //and the prescale value if it was passed
  
  if(edmevent_isRealData) return 1;

  unsigned int myPrescale = 0;
  if(edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v2>0) myPrescale = edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v2_prs;
  if(edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v3>0) myPrescale = edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v3_prs;
  return myPrescale;
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

float getJetPt( unsigned int ijet ) {
  //need to move to using jecFactor 
  return myJetsPF->at(ijet).pt;
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


bool isGoodJet(unsigned int ijet) {

  if ( getJetPt(ijet) <50) return false;
  if ( fabs(myJetsPF->at(ijet).eta) > 2.4) return false;
  if ( !jetPassLooseID(ijet) ) return false;

  return true;
}

bool isGoodJet10(unsigned int ijet) {

  if ( getJetPt(ijet) <10) return false;
  if ( fabs(myJetsPF->at(ijet).eta) > 2.4) return false;
  if ( !jetPassLooseID(ijet) ) return false;

  return true;
}

bool isGoodJet30(unsigned int ijet) {

  if ( getJetPt(ijet) <30) return false;
  if ( fabs(myJetsPF->at(ijet).eta) > 2.4) return false;
  if ( !jetPassLooseID(ijet) ) return false;


  return true;
}

bool isGoodJetMHT(unsigned int ijet) {

  if ( getJetPt(ijet) <30) return false;
  if ( fabs(myJetsPF->at(ijet).eta) > 5) return false;
  //no jet id for MHT

  return true;
}


bool passBTagger(int ijet) {

  if(theBTaggerType_==kSSVM){
    return ( myJetsPF->at(ijet).simpleSecondaryVertexBJetTags >= 1.74 
	     || myJetsPF->at(ijet).simpleSecondaryVertexHighEffBJetTags >= 1.74);
  }
  else if(theBTaggerType_==kTCHET){
    return (myJetsPF->at(ijet).trackCountingHighEffBJetTags >= 10.2);
  }
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


uint nGoodBJets() {
  uint nb=0;
  for (uint i = 0; i < myJetsPF->size(); ++i) {
    if (isGoodJet(i) ) {
      if ( passBTagger(i) ) nb++;
    }
  }
  return nb;
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

  /*
  //weight = Prob( >=1 HLT Tag | m offline true tags and n offline fake tags)

  //should find a more elegant way to do this... (implemented above)

  //weight for offline-tagged b-jets
  if(nb==1 && nnonb==0) 
    w_event = w_bjet;
  else if (nb==2 && nnonb==0) 
    w_event = w_bjet*(1 - w_bjet)*2 + w_bjet*w_bjet;
  else if (nb==3 && nnonb==0)
    w_event = w_bjet*(1 - w_bjet)*(1 - w_bjet)*3
      + w_bjet*w_bjet*(1 - w_bjet)*3 
      + pow(w_bjet,3);
  else if (nb==4 && nnonb==0)
    w_event = w_bjet*pow( (1-w_bjet),3)*4 
      + pow( w_bjet,2)*pow( (1-w_bjet),2)*6 
      + pow(w_bjet,3)*(1-w_bjet)*4 
      + pow(w_bjet,4);

  //weight for offline-tagged nonb-jets
  else if (nnonb==1 && nb==0) 
    w_event = w_nonbjet;
  else if (nnonb==2 && nb==0) 
    w_event = w_nonbjet*(1 - w_nonbjet)*2 + w_nonbjet*w_nonbjet;
  else if (nnonb==3 && nb==0)
    w_event = w_nonbjet*(1 - w_nonbjet)*(1 - w_nonbjet)*3
      + w_nonbjet*w_nonbjet*(1 - w_nonbjet)*3 
      + pow(w_nonbjet,3);
  else if (nnonb==4 && nb==0)
    w_event = w_nonbjet*pow( (1-w_nonbjet),3)*4 
      + pow( w_nonbjet,2)*pow( (1-w_nonbjet),2)*6 
      + pow(w_nonbjet,3)*(1-w_nonbjet)*4 
      + pow(w_nonbjet,4);

  //what if there are mistags and real tags?
  else if(nb==1 && nnonb==1)
    w_event = w_bjet*(1-w_nonbjet) + w_nonbjet*(1-w_bjet) + w_bjet*w_nonbjet;
  else if(nb==2 && nnonb==1)
    w_event = w_bjet*(1-w_bjet)*(1-w_nonbjet)*2 
      + w_bjet*w_bjet*(1-w_nonbjet) 
      + w_nonbjet*(1-w_bjet)*(1-w_bjet) 
      * w_nonbjet*w_bjet*(1-w_bjet)*2 
      + w_bjet*w_bjet*w_nonbjet; 
  else if (nb==1 && nnonb==2)
    w_event = w_nonbjet*(1-w_nonbjet)*(1-w_bjet)*2 
      + w_nonbjet*w_nonbjet*(1-w_bjet) 
      + w_bjet*(1-w_nonbjet)*(1-w_nonbjet) 
      * w_bjet*w_nonbjet*(1-w_nonbjet)*2 
      + w_nonbjet*w_nonbjet*w_bjet; 
  else if (nb==3 && nnonb==1)
    w_event = w_bjet*pow((1-w_bjet),2)*(1-w_nonbjet)*3 + w_bjet*pow((1-w_bjet),2)*w_nonbjet*3
      + w_bjet*w_bjet*(1-w_bjet)*(1-w_nonbjet)*3 + w_bjet*w_bjet*(1-w_bjet)*w_nonbjet*3
      + pow(w_bjet,3)*(1-w_nonbjet) + pow(w_bjet,3)*w_nonbjet
      + pow((1-w_bjet),3)*w_nonbjet;
  else if (nb==1 && nnonb==3)
    w_event = w_nonbjet*pow((1-w_nonbjet),2)*(1-w_bjet)*3 + w_nonbjet*pow((1-w_nonbjet),2)*w_bjet*3
      + w_nonbjet*w_nonbjet*(1-w_nonbjet)*(1-w_bjet)*3 + w_nonbjet*w_nonbjet*(1-w_nonbjet)*w_bjet*3
      + pow(w_nonbjet,3)*(1-w_bjet) + pow(w_nonbjet,3)*w_bjet
      + pow((1-w_nonbjet),3)*w_bjet;
  else if (nb==2 && nnonb==2)
    w_event = w_bjet*(1-w_bjet)*pow( (1-w_nonbjet),2)*2 
      + w_nonbjet*(1-w_nonbjet)*pow( (1- w_bjet),2)*2
      + w_bjet*w_bjet*pow( (1-w_nonbjet),2)
      + w_nonbjet*w_nonbjet*pow( (1-w_bjet),2)
      + w_bjet*(1-w_bjet)*w_nonbjet*(1-w_nonbjet)*4
      + w_bjet*w_bjet*w_nonbjet*(1-w_nonbjet)*2
      + w_nonbjet*w_nonbjet*w_bjet*(1-w_bjet)*2
      + w_bjet*w_bjet*w_nonbjet*w_nonbjet;

  //more than 4 tags should be unlikely....
  else{
    std::cout << "WARNING, BTagIP weight was not computed!" << std::endl;
    std::cout << "\t This event has " << nb << " true tags and " << nnonb << " mistags." << std::endl;
  }
  */

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

float getMET() {
  float myMET=-1;

  myMET = myMETPF->at(0).pt;

  return myMET;
}

float getMETphi() {
  float myMETphi=-99;

  myMETphi = myMETPF->at(0).phi;

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


float getJetPx( unsigned int ijet ) {
  return getJetPt(ijet) * cos(myJetsPF->at(ijet).phi);
}
float getJetPy( unsigned int ijet ) {
  return getJetPt(ijet) * sin(myJetsPF->at(ijet).phi);
}
float getJetPz( unsigned int ijet ) {
  return getJetPt(ijet) * sinh(myJetsPF->at(ijet).eta);

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
float eleet1_,muonpt1_,elephi1_,muonphi1_; //FIXME this is a hack that I don't like!


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





uint countMu() {

  int ngoodmu=0;

  unsigned int nmu = 0;
  nmu = myMuonsPF->size();

  for ( unsigned int i = 0; i< nmu; i++) {

    if(myMuonsPF->at(i).pt > 10
       && fabs(myMuonsPF->at(i).eta)<2.4
       && myMuonsPF->at(i).GlobalMuonPromptTight == 1
       && myMuonsPF->at(i).innerTrack_numberOfValidHits >=11
       && myMuonsPF->at(i).track_hitPattern_numberOfValidPixelHits >= 1
       && fabs(myMuonsPF->at(i).dB) < 0.02
       && fabs(myMuonsPF->at(i).vz - myVertex->at(0).z ) <1
       && (myMuonsPF->at(i).chargedHadronIso 
	   + myMuonsPF->at(i).photonIso 
	   + myMuonsPF->at(i).neutralHadronIso)/myMuonsPF->at(i).pt <0.2 
       ){


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

  unsigned int nmu = 0;
  nmu = myMuonsPF->size();

  for (unsigned int i=0; i < nmu; i++) {
    //if(muonhelper.at(i).badPFmuon) hasBadMuon = true; NEED TO GET THIS BACK IN BEN FIXME
  }

  if(hasBadMuon) return false;
  return true;
}

bool passInconsistentMuonFilter(){

  bool hasInconsistentMuon = false;

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

  if(hasInconsistentMuon) return false;
  return true;

}



uint countEle() {

  int ngoodele=0;

  unsigned int nele = 0;
  nele = myElectronsPF->size();

  for (unsigned int i=0; i < nele; i++) {

    if(myElectronsPF->at(i).pt > 10
       && fabs(myElectronsPF->at(i).superCluster_eta) < 2.5 
       && !(fabs(myElectronsPF->at(i).superCluster_eta) > 1.4442 
	    && fabs(myElectronsPF->at(i).superCluster_eta) < 1.566)
       && myElectronsPF->at(i).gsfTrack_trackerExpectedHitsInner_numberOfLostHits <= 1
       && fabs(myElectronsPF->at(i).dB) < 0.02
       && fabs(myElectronsPF->at(i).vz - myVertex->at(0).z ) <1
       && (myElectronsPF->at(i).chargedHadronIso 
	   + myElectronsPF->at(i).photonIso 
	   + myElectronsPF->at(i).neutralHadronIso)/myElectronsPF->at(i).pt <0.2 
       ){


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

  if (cutTag=="cutTrigger" ) return passHLT();
  if (cutTag=="cutUtilityTrigger" ) return ( utilityHLT_HT300()>0 || utilityHLT_HT300_CentralJet30_BTagIP()>0 );

  if (cutTag=="cutPV") return passPV();
  
  if (cutTag=="cutHT") return getHT()>350;
  
  if (cutTag=="cut3Jets") return (nGoodJets() >= 3);
  
  if (cutTag=="cutMuVeto") return passMuVeto();
  if (cutTag=="cutEleVeto") return passEleVeto();

  if (cutTag == "cutMET") {
    float mymet = getMET();
    return mymet >= 150;
  }

  if (cutTag == "cutDeltaPhi") {
    return ( getMinDeltaPhiMET(3) >= 0.3 );
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
  
  cutTags_.push_back("cutDeltaPhi"); cutNames_[cutTags_.back()]="DeltaPhi";
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
  else if (cutTag == "cutDeltaPhi") cutIsRequired =  true;
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
  if (inname.Contains("TTJets_TuneD6T_7TeV-madgraph-tauola") )                  return 165;
  if (inname.Contains("TToBLNu_TuneZ2_s-channel_7TeV-madgraph") )               return 4.6;
  if (inname.Contains("TToBLNu_TuneZ2_t-channel_7TeV-madgraph") )               return 64.6;
  if (inname.Contains("TToBLNu_TuneZ2_tW-channel_7TeV-madgraph") )              return 10.6;
  if (inname.Contains("WJetsToLNu_TuneZ2_7TeV-madgraph-tauola") )               return 31314;
  if (inname.Contains("WWtoAnything_TuneZ2_7TeV-pythia6-tauola") )              return 43;
  if (inname.Contains("WZtoAnything_TuneZ2_7TeV-pythia6-tauola") )              return 10.4;
  if (inname.Contains("ZinvisibleJets_7TeV-madgraph") )                         return 4500;
  if (inname.Contains("ZZtoAnything_TuneZ2_7TeV-pythia6-tauola") )              return 4.297;
    
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
  if (inname.Contains("LM9_SUSY_sftsht_7TeV-pythia6_2") )                       return "LM9";
  if (inname.Contains("QCD_Pt_0to5_TuneZ2_7TeV_pythia6_3") )                    return inname;
  if (inname.Contains("QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6") )                return inname;
  if (inname.Contains("QCD_Pt_120to170_TuneZ2_7TeV_pythia6") )                  return inname;
  if (inname.Contains("QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6_3") )              return inname;
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
  if (inname.Contains("QCD_Pt_80to120_TuneZ2_7TeV_pythia6_3") )                 return inname;
  if (inname.Contains("TTJets_TuneD6T_7TeV-madgraph-tauola") )                  return "TTbarJets";
  if (inname.Contains("TToBLNu_TuneZ2_s-channel_7TeV-madgraph") )               return "SingleTop-sChannel";
  if (inname.Contains("TToBLNu_TuneZ2_t-channel_7TeV-madgraph") )               return "SingleTop-tChannel";
  if (inname.Contains("TToBLNu_TuneZ2_tW-channel_7TeV-madgraph") )              return "SingleTop-tWChannel";
  if (inname.Contains("WJetsToLNu_TuneZ2_7TeV-madgraph-tauola") )               return "WJets";
  if (inname.Contains("WWtoAnything_TuneZ2_7TeV-pythia6-tauola") )              return "WW";
  if (inname.Contains("WZtoAnything_TuneZ2_7TeV-pythia6-tauola") )              return "WZ";
  if (inname.Contains("ZinvisibleJets_7TeV-madgraph") )                         return "Zinvisible";
  if (inname.Contains("ZZtoAnything_TuneZ2_7TeV-pythia6-tauola") )              return "ZZ";

  //if it isn't found, just use the full name 
  //  -- data will fall into this category, which is fine because we have to hadd the files after anyway
  return inname;

}


double getWeight(Long64_t nentries) {
  if(edmevent_isRealData) return 1;

  double sigma = getCrossSection(sampleName_);
  double w = lumi_ * sigma / double(nentries);

  if (sampleName_.Contains("QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6")) w *=  geneventinfoproduct_weight;

  return  w;
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


void cutflow(itreestream& stream){

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
    // Read event into memory
    stream.read(entry);
    fillObjects();

    //some output to watch while it's running
    if(entry==0){
      if(!edmevent_isRealData){
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

  //we're making an ntuple, so size matters -- use float not double
  double weight; //one exception to the float rule
  float btagIPweight, pfmhtweight;
  ULong64_t lumiSection, eventNumber, runNumber;
  float HT, MHT, MET, METphi, minDeltaPhi, minDeltaPhiAll, minDeltaPhiAll30,minDeltaPhi30_eta5_noIdAll;
  float maxDeltaPhi, maxDeltaPhiAll, maxDeltaPhiAll30, maxDeltaPhi30_eta5_noIdAll;
  float sumDeltaPhi, diffDeltaPhi;

  bool cutHT,cutPV,cutTrigger;
  bool cut3Jets,cutEleVeto,cutMuVeto,cutMET,cutDeltaPhi, cutCleaning;

  bool passBadPFMuon, passInconsistentMuon;

  bool pass_utilityHLT_HT300;
  UInt_t prescale_utilityHLT_HT300;
  bool pass_utilityHLT_HT300_CentralJet30_BTagIP;
  UInt_t prescale_utilityHLT_HT300_CentralJet30_BTagIP;

  int njets, nbjets, nElectrons, nMuons;

  float jetpt1,jetphi1, jeteta1, jetenergy1, bjetpt1, bjetphi1, bjeteta1, bjetenergy1;
  float jetpt2,jetphi2, jeteta2, jetenergy2, bjetpt2, bjetphi2, bjeteta2, bjetenergy2;
  float jetpt3,jetphi3, jeteta3, jetenergy3, bjetpt3, bjetphi3, bjeteta3, bjetenergy3;
  float eleet1;
  float muonpt1;

  float MT_Wlep;

  int nGoodPV;

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

  reducedTree.Branch("weight",&weight,"weight/D");
  reducedTree.Branch("runNumber",&runNumber,"runNumber/l");
  reducedTree.Branch("lumiSection",&lumiSection,"lumiSection/l");
  reducedTree.Branch("eventNumber",&eventNumber,"eventNumber/l");

  reducedTree.Branch("btagIPweight",&btagIPweight,"btagIPweight/F");
  reducedTree.Branch("pfmhtweight",&pfmhtweight,"pfmhtweight/F");

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

  reducedTree.Branch("passBadPFMuon",&passBadPFMuon,"passBadPFMuon/O");
  reducedTree.Branch("passInconsistentMuon",&passInconsistentMuon,"passInconsistentMuon/O");

  reducedTree.Branch("njets",&njets,"njets/I");
  reducedTree.Branch("nbjets",&nbjets,"nbjets/I");
  reducedTree.Branch("nElectrons",&nElectrons,"nElectrons/I");
  reducedTree.Branch("nMuons",&nMuons,"nMuons/I");

  reducedTree.Branch("pass_utilityHLT_HT300",&pass_utilityHLT_HT300,"pass_utilityHLT_HT300/O");
  reducedTree.Branch("prescale_utilityHLT_HT300", &prescale_utilityHLT_HT300, "prescale_utilityHLT_HT300/i");
  reducedTree.Branch("pass_utilityHLT_HT300_CentralJet30_BTagIP",&pass_utilityHLT_HT300_CentralJet30_BTagIP,"pass_utilityHLT_HT300_CentralJet30_BTagIP/O");
  reducedTree.Branch("prescale_utilityHLT_HT300_CentralJet30_BTagIP", &prescale_utilityHLT_HT300_CentralJet30_BTagIP, "prescale_utilityHLT_HT300_CentralJet30_BTagIP/i");  

  reducedTree.Branch("HT",&HT,"HT/F");
  reducedTree.Branch("MET",&MET,"MET/F");
  reducedTree.Branch("METphi",&METphi,"METphi/F");
  reducedTree.Branch("MHT",&MHT,"MHT/F");

  reducedTree.Branch("bestWMass",&bestWMass_,"bestWMass/F");
  reducedTree.Branch("bestTopMass",&bestTopMass_,"bestTopMass/F");
  reducedTree.Branch("topCosHel",&topCosHel_,"topCosHel/F");
  reducedTree.Branch("WCosHel",&WCosHel_,"WCosHel/F");
  reducedTree.Branch("MT_Wlep",&MT_Wlep, "MT_Wlep/F");

  reducedTree.Branch("minDeltaPhi",&minDeltaPhi,"minDeltaPhi/F");
  reducedTree.Branch("minDeltaPhiAll",&minDeltaPhiAll,"minDeltaPhiAll/F");
  reducedTree.Branch("minDeltaPhiAll30",&minDeltaPhiAll30,"minDeltaPhiAll30/F");
  reducedTree.Branch("minDeltaPhi30_eta5_noIdAll",&minDeltaPhi30_eta5_noIdAll,"minDeltaPhi30_eta5_noIdAll/F");

  reducedTree.Branch("maxDeltaPhi",&maxDeltaPhi,"maxDeltaPhi/F");
  reducedTree.Branch("maxDeltaPhiAll",&maxDeltaPhiAll,"maxDeltaPhiAll/F");
  reducedTree.Branch("maxDeltaPhiAll30",&maxDeltaPhiAll30,"maxDeltaPhiAll30/F");
  reducedTree.Branch("maxDeltaPhi30_eta5_noIdAll",&maxDeltaPhi30_eta5_noIdAll,"maxDeltaPhi30_eta5_noIdAll/F");

  reducedTree.Branch("sumDeltaPhi",&sumDeltaPhi,"sumDeltaPhi/F");
  reducedTree.Branch("diffDeltaPhi",&diffDeltaPhi,"diffDeltaPhi/F");

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


  int nevents = stream.size();

  InitializeStuff();//BEN
   
  startTimer();
  for(int entry=0; entry < nevents; ++entry){
    // Read event into memory
    stream.read(entry);
    fillObjects();

    //some output to watch while it's running
    if(entry==0){
      if(!edmevent_isRealData){
	cout << "MC xsec: " << getCrossSection(sampleName_) << endl;
      }
      else{
	cout << "This is data!"<< endl;
      }
      cout << "weight: "<< getWeight(nevents) << endl;
    }
    if(entry%100000==0) cout << "  entry: " << entry << ", percent done=" << (int)(entry/(double)nevents*100.)<<  endl;
      

    if ((passCut("cutTrigger") || passCut("cutUtilityTrigger")) && passCut("cutHT") ) {
           
      //if (entry%1000000==0) checkTimer(entry,nevents);
      weight = getWeight(nevents);
      
      // cast these as long ints, with rounding, assumes they are positive to begin with 
      runNumber = (ULong64_t)(edmevent_run+0.5);
      lumiSection = (ULong64_t)(edmevent_luminosityBlock+0.5);
      eventNumber = (ULong64_t)(edmevent_event+0.5);
      
      btagIPweight = getBTagIPWeight();
      pfmhtweight = getPFMHTWeight();

      cutHT = true; 
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
      
      pass_utilityHLT_HT300 = (utilityHLT_HT300()>0);
      prescale_utilityHLT_HT300 = (UInt_t)utilityHLT_HT300();
      pass_utilityHLT_HT300_CentralJet30_BTagIP = (utilityHLT_HT300_CentralJet30_BTagIP()>0);
      prescale_utilityHLT_HT300_CentralJet30_BTagIP = (UInt_t)utilityHLT_HT300_CentralJet30_BTagIP();

      nGoodPV = countGoodPV();

      njets = nGoodJets();
      nbjets = nGoodBJets();
      nElectrons = countEle();
      nMuons = countMu();
      HT=getHT();
      MET=getMET();
      MHT=getMHT();
      METphi = getMETphi();
      minDeltaPhi = getMinDeltaPhiMET(3);
      minDeltaPhiAll = getMinDeltaPhiMET(99);
      minDeltaPhiAll30 = getMinDeltaPhiMET30(99);
      minDeltaPhi30_eta5_noIdAll = getMinDeltaPhiMET30_eta5_noId(99);
      maxDeltaPhi = getMaxDeltaPhiMET(3);
      maxDeltaPhiAll = getMaxDeltaPhiMET(99);
      maxDeltaPhiAll30 = getMaxDeltaPhiMET30(99);
      maxDeltaPhi30_eta5_noIdAll = getMaxDeltaPhiMET30_eta5_noId(99);
      
      sumDeltaPhi = maxDeltaPhi + minDeltaPhi;
      diffDeltaPhi = maxDeltaPhi - minDeltaPhi;
      
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

  setMuonReq(0);
  setEleReq(0);
  setCutScheme();
  InitializeStuff();//BEN
  setBCut(2);

  int count = 0;
  startTimer();
  for(int entry=0; entry < nevents; ++entry){
    // Read event into memory
    stream.read(entry);
    fillObjects();

    
    if (Cut(entry) < 0) continue;
    count++;
    std::cout << "genweight = " << geneventinfoproduct_weight << std::endl;
    std::cout << "btagIP weight = " << getBTagIPWeight() << std::endl;
    //std::cout << "PFMHT weight = " << getPFMHTWeight() << std::endl;

  }
  stopTimer(nevents);
  std::cout << "count = " << count << std::endl;

  return;
}
