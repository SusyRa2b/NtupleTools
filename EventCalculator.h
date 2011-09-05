// -*- C++ -*-

#ifndef EVENTCALCULATOR_H
#define EVENTCALCULATOR_H

#include "CrossSectionTable.h"

#include "MiscUtil.cxx"

//Pile-up reweighting stuff
//NOTE: Grab this header from PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h
#include "LumiReweightingStandAlone.h"
//JES on the fly
#include "JetCorrectorParameters.h"
#include "FactorizedJetCorrector.h"
#include "JetCorrectionUncertainty.h"

#include "TMath.h"
#include "TGraphAsymmErrors.h"

/*
left undone -- have not tried to port luke's jet killing stuff (should be done)
*/
#include "BasicLoopCU.h" //get all of the tree-related variables



//global constants
const double mW_ = 80.399;
const double mtop_ = 172.0;
const double lumi_ = 1.; //fix to 1/pb and scale MC later (e.g. in drawReducedTrees)

class EventCalculator {
public:
  //enums for configuration
  enum scanType {kNotScan=0, kmSugra, kSMS};
  enum METType {kPFMET=0,kPFMETTypeI};
  enum jetType {kPF2PAT=0, kRECOPF};
  enum JESType {kJES0=0,kJESup,kJESdown,kJESFLY};
  enum JERType {kJER0=0,kJERbias,kJERup,kJERdown,kJERra2};
  enum METuncType {kMETunc0=0,kMETuncDown,kMETuncUp};
  enum PUuncType {kPUunc0=0,kPUuncDown,kPUuncUp};
  enum BTagEffType {kBTagEff0=0,kBTagEffup,kBTagEffdown,kBTagEff02,kBTagEffup2,kBTagEffdown2,kBTagEff03,kBTagEffup3,kBTagEffdown3};
  enum HLTEffType {kHLTEff0=0,kHLTEffup,kHLTEffdown};
  enum BTaggerType {kSSVM=0, kTCHET, kSSVHPT, kTCHPT, kTCHPM, Nbtaggers};

  EventCalculator(const TString & sampleName, jetType theJetType, METType theMETType);
  ~EventCalculator();

  //setters
  void setBTaggerType(BTaggerType btaggertype) {theBTaggerType_ = btaggertype;}
  void setOptions( const TString & opt);

  //loop over events
  void reducedTree(TString outputpath, itreestream& stream);

  // functions that calculate stuff

  double getWeight(Long64_t nentries);
  bool noPUWeight();
  float getPUWeight(reweight::LumiReWeighting lumiWeights);

  bool isGoodMuon(const unsigned int imuon);
  bool isGoodElectron(const unsigned int iele);
  unsigned int countEle() ;
  bool isCleanMuon(const unsigned int imuon);
  unsigned int countMu();

  //  bool passBEFilter(); //TODO migrate this code
  bool passHLT();
  unsigned int utilityHLT_HT300();
  unsigned int utilityHLT_HT300_CentralJet30_BTagIP();
  float getHLTHTeff(float offHT);
  float getHLTMHTeff(float offMET);

  int countGoodPV();
  bool passPV() ;

  float getHT();

  float getMET();
  float getMETphi();

  float getMHT();
  float getMHTphi();

  //  float getMT_Wlep(); //SKIPping this for now

  /* stuff that i'm skipping for now
    void getTransverseThrustVariables(float & thrust, float & thrustPhi, bool addMET);
void getSphericityJetMET(float & lambda1, float & lambda2, float & det,const int jetmax, bool addMET);

*/

  bool passCut(const TString cutTag);

  double getMinDeltaPhiMET(unsigned int maxjets);
  double getTransverseMETError(unsigned int thisJet);
  double getDeltaPhiNMET(unsigned int thisJet);
  double getDeltaPhiMETN( unsigned int goodJetN ); //yes, we have two !
  double getMinDeltaPhiNMET(unsigned int maxjets);
  double getMinDeltaPhiMETN(unsigned int maxjets); //again a duplicate!

  double getMaxDeltaPhiNMET(unsigned int maxjets);
  double getTransverseMETSignificance(unsigned int thisJet);
  double getMaxTransverseMETSignificance(unsigned int maxjets);
  double getMinTransverseMETSignificance(unsigned int maxjets);
  double getMinDeltaPhiMET30(unsigned int maxjets);
  double getMinDeltaPhiMET30_eta5(unsigned int maxjets);
  double getMinDeltaPhiMET30_eta5_noId(unsigned int maxjets);
  double getMaxDeltaPhiMET(unsigned int maxjets);
  double getMaxDeltaPhiMET30(unsigned int maxjets);
  double getMaxDeltaPhiMET30_eta5(unsigned int maxjets);
  double getMaxDeltaPhiMET30_eta5_noId(unsigned int maxjets);

  double getMinDeltaPhiMETTaus();

  void getCorrectedMET(float& correctedMET, float& correctedMETPhi);
  void getUncorrectedMET(float& uncorrectedMET, float& uncorrectedMETPhi);

  int doPBNR(); //particle based noise rejection

  //functions that are mostly used internally
  ULong64_t getRunNumber() {return TMath::Nint( *myEDM_run );}
  ULong64_t getLumiSection() {return TMath::Nint( *myEDM_luminosityBlock );}
  ULong64_t getEventNumber() {return TMath::Nint( *myEDM_event );}

  bool isCleanJet(const unsigned int ijet);
  bool isGoodJet(const unsigned int ijet, const float pTthreshold=50, const float etaMax=2.4); //subset of isCleanJet
  bool isGoodJet10(unsigned int ijet) {return isGoodJet(ijet,10,2.4);}
  bool isGoodJet30(unsigned int ijet) {return isGoodJet(ijet,30,2.4);}
  bool isGoodJetMHT(unsigned int ijet);
  bool passBTagger(int ijet, BTaggerType btagger=Nbtaggers );

  unsigned int nGoodJets();
  unsigned int nGoodJets30();
  unsigned int nGoodBJets( BTaggerType btagger=Nbtaggers);

  void getSmearedUnclusteredMET(float & myMET, float & myMETphi);


  float getJERbiasFactor(unsigned int ijet);
  float getJESUncertainty( unsigned int ijet );

  float getJetPt( unsigned int ijet );
  float getUncorrectedJetPt( unsigned int ijet );
  float getJetPx( unsigned int ijet ) ;
  float getJetPy( unsigned int ijet ) ;
  float getJetPz( unsigned int ijet ) ;
  float getJetEnergy( unsigned int ijet ) ;
  bool jetPassLooseID( unsigned int ijet );

  float jetPtOfN(unsigned int n);
  float jetPhiOfN(unsigned int n);
  float jetEtaOfN(unsigned int n);
  float jetEnergyOfN(unsigned int n);
  float bjetPtOfN(unsigned int n);
  float bjetPhiOfN(unsigned int n);
  float bjetEtaOfN(unsigned int n);
  float bjetEnergyOfN(unsigned int n);

  float elePtOfN(unsigned int n);
  float elePhiOfN(unsigned int n);

  float muonPtOfN(unsigned int n);
  float muonPhiOfN(unsigned int n);

  float getMT_Wlep();
  void calcTopDecayVariables(float & wmass, float & tmass, float & wcoshel, float & tcoshel);
  void calcCosHel(unsigned int j1i, unsigned int j2i, unsigned int j3i, float & wcoshel,float &tcoshel);

  std::pair<float,float> getJERAdjustedMHTxy();

  //MC tools
  unsigned int findSUSYMaternity( unsigned int k );
  unsigned int getSUSYnb(std::vector<unsigned int> &susyb_index);
  unsigned int getSUSYnb();
  SUSYProcess getSUSYProcess();
  std::pair<int,int> getSMSmasses();
  double checkPdfWeightSanity( double a) ;
  void getPdfWeights(const TString & pdfset, Float_t * pdfWeights, TH1D * sumofweights) ;
  double getCrossSection();

  //  double getPFMHTWeight(); //TODO

  double getScanCrossSection( SUSYProcess p, const TString & variation );
  double getSMSScanCrossSection( const double mgluino);
  void calculateTagProb(float &Prob0, float &ProbGEQ1, float &Prob1, float &ProbGEQ2);

  //btag stuff
  float getBTagIPWeight(); //this function should be called *after* offline tagging 
  float jetTagEff(unsigned int ijet, TH1F* h_btageff, TH1F* h_ctageff, TH1F* h_ltageff);

  //other stuff
  double getDeltaPhi(double a, double b) { return jmt::deltaPhi(a,b);}
  bool isSampleRealData() {return jmt::doubleToBool( *myEDM_isRealData);}
  bool passLumiMask();
  //  void loadBEFailEvents(std::vector<int> &vrun, std::vector<int> &vlumi, std::vector<int> &vevent);
  void dumpEvent ();
  TString getSampleNameOutputString();


  //bookkeeping
  TString getCutDescriptionString();

protected:
  double calc_mNj( std::vector<unsigned int> jNi );
  double calc_mNj( unsigned int j1i, unsigned int j2i);
  double calc_mNj( unsigned int j1i, unsigned int j2i, unsigned int j3i);

private:
  TString sampleName_;
  bool sampleIsSignal_;

  //configuration options
  scanType theScanType_;
  METType theMETType_; 
  jetType theJetType_;
  JESType theJESType_;
  JERType theJERType_;
  METuncType theMETuncType_;
  PUuncType thePUuncType_;
  BTagEffType theBTagEffType_;
  HLTEffType theHLTEffType_;
  BTaggerType theBTaggerType_;

  //pointers to the critical jet, etc collections
  std::vector<jet1_s> * myJetsPF;
  std::vector<jethelper_s> * myJetsPFhelper;
  
  std::vector<electron1_s> * myElectronsPF;
  std::vector<muon1_s> * myMuonsPF;
  std::vector<tau_s> * myTausPF;
  std::vector<met1_s> * myMETPF;
  std::vector<vertex_s> * myVertex;
  std::vector<genparticlehelperra2_s> * myGenParticles;
  double* myGenWeight;
  double* myEDM_bunchCrossing;
  double* myEDM_event;
  double* myEDM_isRealData;
  double* myEDM_luminosityBlock;
  double* myEDM_run;

  CrossSectionTable * crossSectionTanb40_10_;
  CrossSectionTable * crossSectionTanb40_05_;
  CrossSectionTable * crossSectionTanb40_20_;

  TFile *smsCrossSectionFile_;
  //stuff for the HT turn-on
  TFile *f_eff_;
  TGraphAsymmErrors *htgraph_;
  TGraphAsymmErrors *htgraphPlus_;
  TGraphAsymmErrors *htgraphMinus_;
  TFile *f_eff_mht_;
  TGraphAsymmErrors *mhtgraph_;
  TGraphAsymmErrors *mhtgraphPlus_;
  TGraphAsymmErrors *mhtgraphMinus_;

  void loadHLTHTeff();
  void loadHLTMHTeff();

  TString getOptPiece(const TString &key, const TString & opt);

  //JES on-the-fly stuff
  JetCorrectorParameters *ResJetPar_;
  JetCorrectorParameters *L3JetPar_;
  JetCorrectorParameters *L2JetPar_; 
  JetCorrectorParameters *L1JetPar_; 
  std::vector<JetCorrectorParameters> vPar_;
  FactorizedJetCorrector *JetCorrector_;
  JetCorrectionUncertainty *jecUnc_;

  TDatime* starttime_;
  void startTimer() {  starttime_ = new TDatime();}
  void stopTimer(const Long64_t ntotal);

  //hold human-readable names for configuration options
  void initEnumNames();
  std::map<METType, TString> theMETNames_;
  std::map<jetType, TString> theJetNames_;
  std::map<JESType, TString> theJESNames_;
  std::map<JERType, TString> theJERNames_;
  std::map<METuncType, TString> theMETuncNames_;
  std::map<PUuncType, TString> thePUuncNames_;
  std::map<BTagEffType, TString> theBTagEffNames_;
  std::map<HLTEffType, TString> theHLTEffNames_;
  std::map<BTaggerType, TString> theBTaggerNames_;



  void checkConsistency();

  void  loadSusyScanCrossSections();

  unsigned int getSeed();

};

#endif
