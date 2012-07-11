// -*- C++ -*-

#ifndef EVENTCALCULATOR_H
#define EVENTCALCULATOR_H

#include "CrossSectionTable.h"


//Pile-up reweighting stuff
//NOTE: Grab this header from PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h
#include "LumiReweightingStandAlone.h"
#include "Lumi3DReWeighting.h"
//JES on the fly
#include "JetCorrectorParameters.h"
#include "FactorizedJetCorrector.h"
#include "JetCorrectionUncertainty.h"

#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TBranch.h"
#include "TChain.h"

#include <set>
//#include "BasicLoopCU.h" //get all of the tree-related variables



class TRandom3;

//global constants
const double mW_ = 80.399;
const double mZ_ = 91.2;
const double mtop_ = 172.0;
const double lumi_ = 1.; //fix to 1/pb and scale MC later (e.g. in drawReducedTrees)


class cellECAL {
public:
  cellECAL(double,double,int);
  double eta;
  double phi;
  int status;
};

class triggerData {
public:
  triggerData() ;
  bool pass;
  int prescale;
  int version;
};

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
  enum BTagEffType {kBTagEff0=0,kBTagEffup,kBTagEffdown,kBTagEff02,kBTagEffup2,kBTagEffdown2,kBTagEff03,kBTagEffup3,kBTagEffdown3,kBTagEff04,kBTagEffup4,kBTagEffdown4};
  enum HLTEffType {kHLTEff0=0,kHLTEffup,kHLTEffdown};
  enum BTaggerType {kSSVM=0, kTCHET, kSSVHPT, kTCHPT, kTCHPM, kCSVM, kCSVL,Nbtaggers};

  enum BTagEffModifier {kBTagModifier0=0,kLFdown,kLFup,kHFdown,kHFup}; //ugg...not in love with this design

  EventCalculator(const TString & sampleName, const std::vector<std::string>  inputFiles, jetType theJetType, METType theMETType);
  ~EventCalculator();

  //setters
  void setBTaggerType(BTaggerType btaggertype) {theBTaggerType_ = btaggertype;}
  void setOptions( const TString & opt);

  //loop over events
  void reducedTree(TString outputpath);

  //  void cutflow(itreestream& stream, int maxevents);
  void sampleAnalyzer();
  void plotBTagEffMC();

  //load external list of event ID's
  void loadEventList(std::vector<int> &vrun, std::vector<int> &vlumi, std::vector<int> &vevent);
  bool inEventList(std::vector<int> &vrun, std::vector<int> &vlumi, std::vector<int> &vevent);

  // functions that calculate stuff
  double getWeight(Long64_t nentries);
  Long64_t getNEventsGenerated();
  float getPUWeight(reweight::LumiReWeighting lumiWeights);
  //  float getPUWeight(/* Lumi3DReWeighting lumiWeights*/); //FIXME CFA

  bool isGoodMuon(const unsigned int imuon, const bool disableRelIso=false, const float ptthreshold=10);
  bool isGoodRecoMuon(const unsigned int imuon, const bool disableRelIso=false, const float ptthreshold=10);
  bool isGoodElectron(const unsigned int iele, const bool disableRelIso=false, const float ptthreshold=10);
  unsigned int countEle(const float ptthreshold=10) ;
  bool isCleanMuon(const unsigned int imuon, const float ptthreshold=10);
  unsigned int countMu(const float ptthreshold=10);

  bool isGoodTau(const unsigned int itau, const float pTthreshold=20, const float etaMax=2.4);
  float getTauPt( unsigned int itau );
  unsigned int countTau();


  TString stripTriggerVersion(const TString & fullname, int & version);
  bool passHLT(std::map<TString,triggerData> & triggerlist);
  //  bool passUtilityHLT(int &version, int &prescale); //deprecated
  bool passUtilityPrescaleModuleHLT();
  float getHLTHTeff(float offHT);
  float getHLTMHTeff(float offMET, float offHT, uint nElectrons, uint nMuons, double mindphin);
  //  double getPFMHTWeight(); //TODO
  double getHLTMHTeffBNN(float offMET, float offHT, uint nElectrons, uint nMuons, double mindphin, 
			 double& effUp, double& effDown);

  int countGoodPV();
  bool passPV() ;

  float getHT(float ptthreshold=50);
  float getST(float jetthreshold=50,float leptonthreshold=10);

  float getMET();
  float getMETphi();

  float getMHT();
  float getMHTphi();
  float getMHTphi(int ignoredJet);

  //  void getTransverseThrustVariables(float & thrust, float & thrustPhi, bool addMET);
  void getSphericityJetMET(float & lambda1, float & lambda2, float & det,const int jetmax, bool addMET);

  std::vector<unsigned int> jetsetToVector(const std::vector<unsigned int> & goodjets, const std::set<unsigned int> & myset) ;
  void jjResonanceFinder(float & mjj1, float & mjj2);//simple first try
  void jjResonanceFinder5(float & mjj1, float & mjj2);

  unsigned int getNthGoodJet(unsigned int goodJetN, float mainpt, float maineta, bool mainid);
  double getMinDeltaPhiMET(unsigned int maxjets);
  double getTransverseMETError(unsigned int thisJet);
  double getDeltaPhiMET(unsigned int n, float ptThreshold = 50, bool bjetsonly = false);
  double getDeltaPhiNMET(unsigned int thisJet); //Luke
  double getDeltaPhiMETN_deltaT(unsigned int ijet, float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith);
  double getDeltaPhiMETN_deltaT(unsigned int ijet) { return getDeltaPhiMETN_deltaT(ijet,30,2.4,true,false,false);  } //overloaded
  double getDeltaPhiMETN( unsigned int goodJetN, float mainpt, float maineta, bool mainid, float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith ); //Ben
  double getDeltaPhiMETN( unsigned int goodJetN ) {return getDeltaPhiMETN(goodJetN,50,2.4,true,30,2.4,true,false,false); }; //Ben, overloaded
  double getMinDeltaPhiNMET(unsigned int maxjets); //Luke
  double getMinDeltaPhiMETN(unsigned int maxjets, float mainmt, float maineta, bool mainid, float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith, bool includeLeptons=false ); //Ben
  double getMinDeltaPhiMETN(unsigned int maxjets) {return getMinDeltaPhiMETN(maxjets,50,2.4,true,30,2.4,true,false,false); }; //Ben, overloaded

  double getMaxDeltaPhiNMET(unsigned int maxjets);
  double getTransverseMETSignificance(unsigned int thisJet);
  double getMaxTransverseMETSignificance(unsigned int maxjets);
  double getMinTransverseMETSignificance(unsigned int maxjets);
  double getMinDeltaPhiMET30(unsigned int maxjets, bool bjetsonly = false);
  double getMinDeltaPhiMET30_eta5(unsigned int maxjets);
  double getMinDeltaPhiMET30_eta5_noId(unsigned int maxjets);
  double getMaxDeltaPhiMET(unsigned int maxjets);
  double getMaxDeltaPhiMET30(unsigned int maxjets);
  double getMaxDeltaPhiMET30_eta5(unsigned int maxjets);
  double getMaxDeltaPhiMET30_eta5_noId(unsigned int maxjets);

  double getDeltaPhiMETN_electron_deltaT(unsigned int ielectron, float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith);
  double getDeltaPhiMETN_muon_deltaT(unsigned int imuon, float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith);
  double getDeltaPhiMETN_electron( const unsigned int ielectron, float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith) ;
  double getDeltaPhiMETN_muon( const unsigned int imuon, float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith) ;

  float getMaxJetMis(unsigned int rank, unsigned int maxjets, float jetpt);
  float getMaxJetFracMis(unsigned int rank, unsigned int maxjets, float jetpt);
  float getDeltaPhiMETJetMaxMis(float jetpt);

  float getDeltaPhiStar(int & badjet);

  double getMinDeltaPhiMETTaus();
  double getMinDeltaPhiMETMuons(unsigned int maxmuons);

  void getCorrectedMET(float& correctedMET, float& correctedMETPhi);
  void getUncorrectedMET(float& uncorrectedMET, float& uncorrectedMETPhi);

  int doPBNR(); //particle based noise rejection

  //functions that are mostly used internally
  ULong64_t getRunNumber() {return run;}
  ULong64_t getLumiSection() {return lumiblock;}
  ULong64_t getEventNumber() {return event;}

  bool isCleanJet(const unsigned int ijet);
  bool isGoodJet(const unsigned int ijet, const float pTthreshold=50, const float etaMax=2.4, const bool jetid=true); //subset of isCleanJet
  bool isGoodJet10(unsigned int ijet) {return isGoodJet(ijet,10,2.4,true);}
  bool isGoodJet30(unsigned int ijet) {return isGoodJet(ijet,30,2.4,true);}
  bool isGoodJetMHT(unsigned int ijet);
  bool passBTagger(int ijet, BTaggerType btagger=Nbtaggers );

  unsigned int nGoodJets();
  //  unsigned int nGoodJets(TH2D* count,TH2D* unc,TH2D* l2l3); //for a test
  unsigned int nGoodJets30();
  unsigned int nGoodBJets( BTaggerType btagger=Nbtaggers);
  unsigned int nTrueBJets();

  void getSmearedUnclusteredMET(float & myMET, float & myMETphi);


  float getJERbiasFactor(unsigned int ijet);
  float getJESUncertainty( unsigned int ijet, bool addL2L3toJES );

  float getJetPt( unsigned int ijet, bool addL2L3toJES=false );
  float getUncorrectedJetPt( unsigned int ijet, bool addL2L3toJES=false );
  float getJetPx( unsigned int ijet ) ;
  float getJetPy( unsigned int ijet ) ;
  float getJetPz( unsigned int ijet ) ;
  float getJetEnergy( unsigned int ijet ) ;
  bool jetPassLooseID( unsigned int ijet );
  float getJetCSV( unsigned int ijet );

  float jetPtOfN(unsigned int n);
  float jetGenPtOfN(unsigned int n);
  float jetPhiOfN(unsigned int n);
  float jetGenPhiOfN(unsigned int n);
  float jetEtaOfN(unsigned int n);
  float jetGenEtaOfN(unsigned int n);
  float jetEnergyOfN(unsigned int n);
  int jetFlavorOfN(unsigned int n);
  float jetChargedHadronFracOfN(unsigned int n);
  int jetChargedHadronMultOfN(unsigned int n);

  float bjetPtOfN(unsigned int n);
  float bjetPhiOfN(unsigned int n);
  float bjetEtaOfN(unsigned int n);
  float bjetEnergyOfN(unsigned int n);
  float bjetCSVOfN(unsigned int n);
  int bjetFlavorOfN(unsigned int n);
  float bjetChargedHadronFracOfN(unsigned int n);
  int bjetChargedHadronMultOfN(unsigned int n);

  float elePtOfN(unsigned int n, const float ptthreshold=10);
  float eleEtaOfN(unsigned int n, const float ptthreshold=10);
  float elePhiOfN(unsigned int n, const float ptthreshold=10);

  float muonPtOfN(unsigned int n, const float ptthreshold=10);
  float muonEtaOfN(unsigned int n, const float ptthreshold=10);
  float muonPhiOfN(unsigned int n, const float ptthreshold=10);
  float muonIsoOfN(unsigned int n, const float ptthreshold=10);
  float muonChHadIsoOfN(unsigned int n, const float ptthreshold=10);
  float muonPhotonIsoOfN(unsigned int n, const float ptthreshold=10);
  float muonNeutralHadIsoOfN(unsigned int n, const float ptthreshold=10);

  float getBestZCandidate(const float pt_threshold1,const float pt_threshold2);
  float getBestZeeCandidate(const float pt_threshold1,const float pt_threshold2);
  float getBestZmmCandidate(const float pt_threshold1,const float pt_threshold2);

  float tauPtOfN(unsigned int n);
  float tauEtaOfN(unsigned int n);

  float recoMuonPtOfN(unsigned int n, const float ptthreshold=10);
  float recoMuonEtaOfN(unsigned int n, const float ptthreshold=10);
  float recoMuonPhiOfN(unsigned int n, const float ptthreshold=10);
  float recoMuonIsoOfN(unsigned int n, const float ptthreshold=10);
  float recoMuonMinDeltaPhiJetOfN(unsigned int n, const float ptthreshold=10);

  float getRelIsoForIsolationStudyEle();
  float getRelIsoForIsolationStudyMuon();

  float getMT_Wlep(const float pttreshold=10);
  void calcTopDecayVariables(float & wmass, float & tmass, float & wcoshel, float & tcoshel);
  void calcCosHel(unsigned int j1i, unsigned int j2i, unsigned int j3i, float & wcoshel,float &tcoshel);

  std::pair<float,float> getJERAdjustedMHTxy(int ignoredJet=-1);

  //MC tools
  unsigned int findSUSYMaternity( unsigned int k );
  unsigned int getSUSYnb(std::vector<unsigned int> &susyb_index);
  unsigned int getSUSYnb();
  SUSYProcess getSUSYProcess(float & pt1, float & phi1, float & pt2, float & phi2); //momenta of the SUSY mothers
  std::pair<int,int> getSMSmasses();
  double checkPdfWeightSanity( double a) ;
  //  void getPdfWeights(const TString & pdfset, Float_t * pdfWeights, TH1D * sumofweights) ;
  int findTop(int& top1, int& top2);
  
  int WDecayType(const int Wparent,int& Wdaughter);
  int findW(int& W, int& Wdaughter, int parent, bool fromtop);
  int muonMatch(const int trueMuon);
  int electronMatch(const int trueElectron);
  int tauMatch(const int trueTau);
  int daughterMatch(const int Wdaughter, const int WdecayType);
  //  int getTTbarDecayType(int& W1decayType, int& W2decayType, int& W1, int& W1daughter, int& W2, int& W2daughter, bool passW2info);
  //  int getWDecayType(int& WdecayType, int& W, int& Wdaughter, bool fromtop);

  double getCrossSection();

  double getScanCrossSection( SUSYProcess p, const TString & variation );
  double getSMSScanCrossSection( const double mgluino);
  void calculateTagProb(float &Prob0, float &ProbGEQ1, float &Prob1, float &ProbGEQ2, float & Prob2, float &ProbGEQ3, 
			float extraSFb=1, float extraSFc=1, float extraSFl=1, BTagEffModifier modifier=kBTagModifier0);
  //  void averageBeff(double & bjetEffSum);// , Long64_t & bjetSum);

  //btag stuff
  //  float getBTagIPWeight(); //this function should be called *after* offline tagging 
  float jetTagEff(unsigned int ijet, TH1F* h_btageff, TH1F* h_ctageff, TH1F* h_ltageff,
		  const float extraSFb, const float extraSFc, const float extraSFl,const BTagEffModifier modifier=kBTagModifier0);

  //other stuff
  double getDeltaPhi(double a, double b);
  bool isSampleRealData();
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


  CrossSectionTable * crossSectionTanb40_10_;
  CrossSectionTable * crossSectionTanb40_05_;
  CrossSectionTable * crossSectionTanb40_20_;

  TFile *smsCrossSectionFile_;
  //stuff for the HT turn-on
  TFile *f_eff_ht300_;
  TFile *f_eff_ht350_;
  TGraphAsymmErrors *htgraph_ht300_;
  TGraphAsymmErrors *htgraphPlus_ht300_;
  TGraphAsymmErrors *htgraphMinus_ht300_;
  TGraphAsymmErrors *htgraph_ht350_;
  TGraphAsymmErrors *htgraphPlus_ht350_;
  TGraphAsymmErrors *htgraphMinus_ht350_;
  TFile *f_eff_mht_;
  TGraphAsymmErrors *mhtgraph_;
  TGraphAsymmErrors *mhtgraphPlus_;
  TGraphAsymmErrors *mhtgraphMinus_;

  //stuff for the btag probability 
  TFile *f_tageff_;
  void loadJetTagEffMaps();

  void loadHLTHTeff();
  void loadHLTMHTeff();

  void loadDataJetRes();
  float getDataJetRes(float pt, float eta);
  TFile *fDataJetRes_;
  TH1D *hEta0_0p5_;
  TH1D *hEta0p5_1_;
  TH1D *hEta1_1p5_;
  TH1D *hEta1p5_2_;
  TH1D *hEta2_2p5_;
  TH1D *hEta2p5_3_;
  TH1D *hEta3_5_;

  //ecal dead cell
  std::vector<cellECAL> badECAL_;
  void loadECALStatus();
  bool passBadECALFilter(TString nearmettype = "METphi", double nearmetphicut=0.5,TString nearecaltype ="deadCell", double nearecalRcut=0.3,int nbadcellsthr=10, int badcellstatusth=10);
  float passBadECALFilter_worstMis(TString nearmettype = "METphi", double nearmetphicut=0.5,TString nearecaltype ="deadCell", double nearecalRcut=0.3, int nbadcellsthr=10, int badcellstatusth=10, TString mistype ="abs");
  bool jetNearMET(unsigned int i, TString type, double nearmetphicut);
  bool jetNearBadECALCell(unsigned int i, TString type, double nearecalRcut, int nbadcellsthr, int badcellstatusthr);

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



  //  void checkConsistency();

  void  loadSusyScanCrossSections();

  unsigned int getSeed();
  //void changeVariables(TRandom3* random, double jetLossProbability, int& nLostJets); //FIXME CFA
  //void resetVariables(); //FIXME CFA
  bool  recalculatedVariables_;

  TStopwatch* watch_; //not used for everyday running; just for testing


  // ==== BEGIN giant copy/paste of cfA variables =============================================================
  // generated from slimmed v60 --
  //source: http://cms2.physics.ucsb.edu/cfA/WJetsToLNu_HT-400ToInf_8TeV-madgraph_Summer12-PU_S7_START52_V9-v1_AODSIM_UCSB1229_v60s/WJetsToLNu_HT-400ToInf_8TeV-madgraph_Summer12-PU_S7_START52_V9-v1_AODSIM_UCSB1229_v60s.h

   // Declaration of leaf types
   std::vector<float>   *trigger_prescalevalue;
   std::vector<std::string>  *trigger_name;
   std::vector<float>   *trigger_decision;
   std::vector<std::vector<float> > *triggerobject_pt;
   std::vector<std::vector<float> > *triggerobject_phi;
   std::vector<std::vector<float> > *triggerobject_eta;
   std::vector<float>   *L1trigger_prescalevalue;
   std::vector<std::string>  *L1trigger_name;
   std::vector<float>   *L1trigger_decision;
   std::vector<float>   *els_conversion_dist;
   std::vector<float>   *els_conversion_dcot;
   Int_t           hbhefilter_decision;
   Int_t           trackingfailurefilter_decision;
   Int_t           cschalofilter_decision;
   Int_t           ecalTPfilter_decision;
   Int_t           ecalBEfilter_decision;
   Int_t           scrapingVeto_decision;
   Int_t           greedymuonfilter_decision;
   Int_t           inconsistentPFmuonfilter_decision;
   Int_t           eenoisefilter_decision;
   Int_t           passprescalePFHT350filter_decision;
   std::vector<float>   *jets_AK5PFclean_corrL2L3;
   std::vector<float>   *jets_AK5PFclean_corrL2L3Residual;
   std::vector<float>   *jets_AK5PFclean_corrL1FastL2L3;
   std::vector<float>   *jets_AK5PFclean_corrL1L2L3;
   std::vector<float>   *jets_AK5PFclean_corrL1FastL2L3Residual;
   std::vector<float>   *jets_AK5PFclean_corrL1L2L3Residual;
   std::vector<std::vector<float> > *PU_zpositions;
   std::vector<std::vector<float> > *PU_sumpT_lowpT;
   std::vector<std::vector<float> > *PU_sumpT_highpT;
   std::vector<std::vector<int> > *PU_ntrks_lowpT;
   std::vector<std::vector<int> > *PU_ntrks_highpT;
   std::vector<int>     *PU_NumInteractions;
   std::vector<int>     *PU_bunchCrossing;
   std::vector<int>     *PU_TrueNumInteractions;

   // List of branches
   TBranch        *b_trigger_prescalevalue;   //!
   TBranch        *b_trigger_name;   //!
   TBranch        *b_trigger_decision;   //!
   TBranch        *b_triggerobject_pt;   //!
   TBranch        *b_triggerobject_phi;   //!
   TBranch        *b_triggerobject_eta;   //!
   TBranch        *b_L1trigger_prescalevalue;   //!
   TBranch        *b_L1trigger_name;   //!
   TBranch        *b_L1trigger_decision;   //!
   TBranch        *b_els_conversion_dist;   //!
   TBranch        *b_els_conversion_dcot;   //!
   TBranch        *b_hbhefilter_decision;   //!
   TBranch        *b_trackingfailurefilter_decision;   //!
   TBranch        *b_cschalofilter_decision;   //!
   TBranch        *b_ecalTPfilter_decision;   //!
   TBranch        *b_ecalBEfilter_decision;   //!
   TBranch        *b_scrapingVeto_decision;   //!
   TBranch        *b_greedymuonfilter_decision;   //!
   TBranch        *b_inconsistentPFmuonfilter_decision;   //!
   TBranch        *b_eenoisefilter_decision;   //!
   TBranch        *b_passprescalePFHT350filter_decision;   //!
   TBranch        *b_jets_AK5PFclean_corrL2L3;   //!
   TBranch        *b_jets_AK5PFclean_corrL2L3Residual;   //!
   TBranch        *b_jets_AK5PFclean_corrL1FastL2L3;   //!
   TBranch        *b_jets_AK5PFclean_corrL1L2L3;   //!
   TBranch        *b_jets_AK5PFclean_corrL1FastL2L3Residual;   //!
   TBranch        *b_jets_AK5PFclean_corrL1L2L3Residual;   //!
   TBranch        *b_PU_zpositions;   //!
   TBranch        *b_PU_sumpT_lowpT;   //!
   TBranch        *b_PU_sumpT_highpT;   //!
   TBranch        *b_PU_ntrks_lowpT;   //!
   TBranch        *b_PU_ntrks_highpT;   //!
   TBranch        *b_PU_NumInteractions;   //!
   TBranch        *b_PU_bunchCrossing;   //!
   TBranch        *b_PU_TrueNumInteractions;   //!

   // Declaration of leaf types
   std::vector<float>   *beamSpot_x;
   std::vector<float>   *beamSpot_y;
   std::vector<float>   *els_energy;
   std::vector<float>   *els_et;
   std::vector<float>   *els_eta;
   std::vector<float>   *els_phi;
   std::vector<float>   *els_pt;
   std::vector<float>   *els_px;
   std::vector<float>   *els_py;
   std::vector<float>   *els_pz;
   std::vector<float>   *els_robustTightId;
   std::vector<float>   *els_simpleEleId95relIso;
   std::vector<float>   *els_simpleEleId90relIso;
   std::vector<float>   *els_simpleEleId85relIso;
   std::vector<float>   *els_simpleEleId80relIso;
   std::vector<float>   *els_simpleEleId70relIso;
   std::vector<float>   *els_simpleEleId95cIso;
   std::vector<float>   *els_simpleEleId90cIso;
   std::vector<float>   *els_simpleEleId85cIso;
   std::vector<float>   *els_simpleEleId80cIso;
   std::vector<float>   *els_simpleEleId70cIso;
   std::vector<float>   *els_cIso;
   std::vector<float>   *els_tIso;
   std::vector<float>   *els_ecalIso;
   std::vector<float>   *els_hcalIso;
   std::vector<float>   *els_charge;
   std::vector<float>   *els_hadOverEm;
   std::vector<float>   *els_eOverPIn;
   std::vector<float>   *els_sigmaIEtaIEta;
   std::vector<float>   *els_scEnergy;
   std::vector<float>   *els_scEta;
   std::vector<float>   *els_scE1x5;
   std::vector<float>   *els_scE2x5Max;
   std::vector<float>   *els_scE5x5;
   std::vector<float>   *els_isEB;
   std::vector<float>   *els_isEE;
   std::vector<float>   *els_dEtaIn;
   std::vector<float>   *els_dPhiIn;
   std::vector<float>   *els_dEtaOut;
   std::vector<float>   *els_dPhiOut;
   std::vector<float>   *els_numlosthits;
   std::vector<float>   *els_tk_phi;
   std::vector<float>   *els_d0dum;
   std::vector<float>   *els_vx;
   std::vector<float>   *els_vy;
   std::vector<float>   *els_vz;
   std::vector<float>   *els_ptError;
   std::vector<float>   *els_n_inner_layer;
   std::vector<float>   *els_dr03EcalRecHitSumEt;
   std::vector<float>   *els_dr03HcalTowerSumEt;
   std::vector<float>   *els_dr03HcalDepth1TowerSumEt;
   std::vector<float>   *els_dr03HcalDepth2TowerSumEt;
   std::vector<float>   *els_dr03TkSumPt;
   std::vector<float>   *jets_AK5PF_phi;
   std::vector<float>   *jets_AK5PF_pt;
   std::vector<float>   *jets_AK5PF_pz;
   std::vector<float>   *jets_AK5PF_px;
   std::vector<float>   *jets_AK5PF_py;
   std::vector<float>   *jets_AK5PF_eta;
   std::vector<float>   *jets_AK5PF_et;
   std::vector<float>   *jets_AK5PF_energy;
   std::vector<float>   *jets_AK5PF_parton_Id;
   std::vector<float>   *jets_AK5PF_parton_motherId;
   std::vector<float>   *jets_AK5PF_gen_pt;
   std::vector<float>   *jets_AK5PF_gen_phi;
   std::vector<float>   *jets_AK5PF_partonFlavour;
   std::vector<float>   *jets_AK5PF_btag_TC_highPur;
   std::vector<float>   *jets_AK5PF_btag_TC_highEff;
   std::vector<float>   *jets_AK5PF_btag_jetProb;
   std::vector<float>   *jets_AK5PF_btag_jetBProb;
   std::vector<float>   *jets_AK5PF_btag_secVertexHighPur;
   std::vector<float>   *jets_AK5PF_btag_secVertexHighEff;
   std::vector<float>   *jets_AK5PF_btag_secVertexCombined;
   std::vector<float>   *jets_AK5PF_jetCharge;
   std::vector<float>   *jets_AK5PF_chgEmE;
   std::vector<float>   *jets_AK5PF_chgHadE;
   std::vector<float>   *jets_AK5PF_photonEnergy;
   std::vector<float>   *jets_AK5PF_chg_Mult;
   std::vector<float>   *jets_AK5PF_neutralEmE;
   std::vector<float>   *jets_AK5PF_neutralHadE;
   std::vector<float>   *jets_AK5PF_neutral_Mult;
   std::vector<float>   *jets_AK5PF_mu_Mult;
   std::vector<float>   *jets_AK5PF_ehf;
   std::vector<float>   *jets_AK5PF_corrFactorRaw;
   std::vector<float>   *jets_AK5PFclean_phi;
   std::vector<float>   *jets_AK5PFclean_pt;
   std::vector<float>   *jets_AK5PFclean_pz;
   std::vector<float>   *jets_AK5PFclean_px;
   std::vector<float>   *jets_AK5PFclean_py;
   std::vector<float>   *jets_AK5PFclean_eta;
   std::vector<float>   *jets_AK5PFclean_et;
   std::vector<float>   *jets_AK5PFclean_energy;
   std::vector<float>   *jets_AK5PFclean_parton_Id;
   std::vector<float>   *jets_AK5PFclean_parton_motherId;
   std::vector<float>   *jets_AK5PFclean_gen_pt;
   std::vector<float>   *jets_AK5PFclean_gen_phi;
   std::vector<float>   *jets_AK5PFclean_partonFlavour;
   std::vector<float>   *jets_AK5PFclean_btag_TC_highPur;
   std::vector<float>   *jets_AK5PFclean_btag_TC_highEff;
   std::vector<float>   *jets_AK5PFclean_btag_jetProb;
   std::vector<float>   *jets_AK5PFclean_btag_jetBProb;
   std::vector<float>   *jets_AK5PFclean_btag_secVertexHighPur;
   std::vector<float>   *jets_AK5PFclean_btag_secVertexHighEff;
   std::vector<float>   *jets_AK5PFclean_btag_secVertexCombined;
   std::vector<float>   *jets_AK5PFclean_jetCharge;
   std::vector<float>   *jets_AK5PFclean_chgEmE;
   std::vector<float>   *jets_AK5PFclean_chgHadE;
   std::vector<float>   *jets_AK5PFclean_photonEnergy;
   std::vector<float>   *jets_AK5PFclean_chg_Mult;
   std::vector<float>   *jets_AK5PFclean_neutralEmE;
   std::vector<float>   *jets_AK5PFclean_neutralHadE;
   std::vector<float>   *jets_AK5PFclean_neutral_Mult;
   std::vector<float>   *jets_AK5PFclean_mu_Mult;
   std::vector<float>   *jets_AK5PFclean_ehf;
   std::vector<float>   *jets_AK5PFclean_corrFactorRaw;
   std::vector<float>   *mc_doc_id;
   std::vector<float>   *mc_doc_pt;
   std::vector<float>   *mc_doc_px;
   std::vector<float>   *mc_doc_py;
   std::vector<float>   *mc_doc_pz;
   std::vector<float>   *mc_doc_eta;
   std::vector<float>   *mc_doc_phi;
   std::vector<float>   *mc_doc_mother_id;
   std::vector<float>   *mc_doc_grandmother_id;
   std::vector<float>   *mc_doc_ggrandmother_id;
   std::vector<float>   *mc_doc_mass;
   std::vector<float>   *mc_electrons_pt;
   std::vector<float>   *mc_electrons_px;
   std::vector<float>   *mc_electrons_py;
   std::vector<float>   *mc_electrons_pz;
   std::vector<float>   *mc_electrons_eta;
   std::vector<float>   *mc_electrons_phi;
   std::vector<float>   *mc_electrons_energy;
   std::vector<float>   *mc_electrons_mother_id;
   std::vector<float>   *mc_electrons_grandmother_id;
   std::vector<float>   *mc_mus_pt;
   std::vector<float>   *mc_mus_px;
   std::vector<float>   *mc_mus_py;
   std::vector<float>   *mc_mus_pz;
   std::vector<float>   *mc_mus_eta;
   std::vector<float>   *mc_mus_phi;
   std::vector<float>   *mc_mus_energy;
   std::vector<float>   *mc_mus_mother_id;
   std::vector<float>   *mc_mus_grandmother_id;
   std::vector<float>   *mc_nues_pt;
   std::vector<float>   *mc_nues_px;
   std::vector<float>   *mc_nues_py;
   std::vector<float>   *mc_nues_pz;
   std::vector<float>   *mc_nues_eta;
   std::vector<float>   *mc_nues_phi;
   std::vector<float>   *mc_nues_energy;
   std::vector<float>   *mc_nues_mother_id;
   std::vector<float>   *mc_nues_grandmother_id;
   std::vector<float>   *mc_numus_pt;
   std::vector<float>   *mc_numus_px;
   std::vector<float>   *mc_numus_py;
   std::vector<float>   *mc_numus_pz;
   std::vector<float>   *mc_numus_eta;
   std::vector<float>   *mc_numus_phi;
   std::vector<float>   *mc_numus_energy;
   std::vector<float>   *mc_numus_mother_id;
   std::vector<float>   *mc_numus_grandmother_id;
   std::vector<float>   *mc_nutaus_pt;
   std::vector<float>   *mc_nutaus_px;
   std::vector<float>   *mc_nutaus_py;
   std::vector<float>   *mc_nutaus_pz;
   std::vector<float>   *mc_nutaus_eta;
   std::vector<float>   *mc_nutaus_phi;
   std::vector<float>   *mc_nutaus_energy;
   std::vector<float>   *mc_nutaus_mother_id;
   std::vector<float>   *mc_nutaus_grandmother_id;
   std::vector<float>   *mets_AK5_et;
   std::vector<float>   *mets_AK5_phi;
   std::vector<float>   *mus_energy;
   std::vector<float>   *mus_et;
   std::vector<float>   *mus_eta;
   std::vector<float>   *mus_phi;
   std::vector<float>   *mus_pt;
   std::vector<float>   *mus_px;
   std::vector<float>   *mus_py;
   std::vector<float>   *mus_pz;
   std::vector<float>   *mus_numberOfMatchedStations;
   std::vector<float>   *mus_cIso;
   std::vector<float>   *mus_tIso;
   std::vector<float>   *mus_ecalIso;
   std::vector<float>   *mus_hcalIso;
   std::vector<float>   *mus_ecalvetoDep;
   std::vector<float>   *mus_hcalvetoDep;
   std::vector<float>   *mus_pfIsolationR03_sumChargedHadronPt;
   std::vector<float>   *mus_pfIsolationR03_sumChargedParticlePt;
   std::vector<float>   *mus_pfIsolationR03_sumNeutralHadronEt;
   std::vector<float>   *mus_pfIsolationR03_sumNeutralHadronEtHighThreshold;
   std::vector<float>   *mus_pfIsolationR03_sumPhotonEt;
   std::vector<float>   *mus_pfIsolationR03_sumPhotonEtHighThreshold;
   std::vector<float>   *mus_pfIsolationR03_sumPUPt;
   std::vector<float>   *mus_pfIsolationR04_sumChargedHadronPt;
   std::vector<float>   *mus_pfIsolationR04_sumChargedParticlePt;
   std::vector<float>   *mus_pfIsolationR04_sumNeutralHadronEt;
   std::vector<float>   *mus_pfIsolationR04_sumNeutralHadronEtHighThreshold;
   std::vector<float>   *mus_pfIsolationR04_sumPhotonEt;
   std::vector<float>   *mus_pfIsolationR04_sumPhotonEtHighThreshold;
   std::vector<float>   *mus_pfIsolationR04_sumPUPt;
   std::vector<float>   *mus_charge;
   std::vector<float>   *mus_cm_chi2;
   std::vector<float>   *mus_cm_ndof;
   std::vector<float>   *mus_cm_pt;
   std::vector<float>   *mus_cm_ptErr;
   std::vector<float>   *mus_tk_chi2;
   std::vector<float>   *mus_tk_ndof;
   std::vector<float>   *mus_tk_pt;
   std::vector<float>   *mus_tk_phi;
   std::vector<float>   *mus_tk_d0dum;
   std::vector<float>   *mus_tk_vz;
   std::vector<float>   *mus_tk_numvalhits;
   std::vector<float>   *mus_tk_ptErr;
   std::vector<float>   *mus_tk_numvalPixelhits;
   std::vector<float>   *mus_tk_numpixelWthMeasr;
   std::vector<float>   *mus_stamu_pt;
   std::vector<float>   *mus_stamu_ptErr;
   std::vector<float>   *mus_num_matches;
   std::vector<float>   *mus_isPFMuon;
   std::vector<float>   *mus_isTrackerMuon;
   std::vector<float>   *mus_isGlobalMuon;
   std::vector<float>   *mus_id_AllGlobalMuons;
   std::vector<float>   *mus_id_AllTrackerMuons;
   std::vector<float>   *mus_id_GlobalMuonPromptTight;
   std::vector<float>   *pfTypeImets_et;
   std::vector<float>   *pfTypeImets_phi;
   std::vector<float>   *pfTypeImets_ex;
   std::vector<float>   *pfTypeImets_ey;
   std::vector<float>   *pfTypeImets_gen_et;
   std::vector<float>   *pfTypeImets_gen_phi;
   std::vector<float>   *pfTypeImets_sumEt;
   std::vector<float>   *pf_els_energy;
   std::vector<float>   *pf_els_et;
   std::vector<float>   *pf_els_eta;
   std::vector<float>   *pf_els_phi;
   std::vector<float>   *pf_els_pt;
   std::vector<float>   *pf_els_px;
   std::vector<float>   *pf_els_py;
   std::vector<float>   *pf_els_pz;
   std::vector<float>   *pf_els_robustTightId;
   std::vector<float>   *pf_els_simpleEleId95relIso;
   std::vector<float>   *pf_els_simpleEleId90relIso;
   std::vector<float>   *pf_els_simpleEleId85relIso;
   std::vector<float>   *pf_els_simpleEleId80relIso;
   std::vector<float>   *pf_els_simpleEleId70relIso;
   std::vector<float>   *pf_els_simpleEleId95cIso;
   std::vector<float>   *pf_els_simpleEleId90cIso;
   std::vector<float>   *pf_els_simpleEleId85cIso;
   std::vector<float>   *pf_els_simpleEleId80cIso;
   std::vector<float>   *pf_els_simpleEleId70cIso;
   std::vector<float>   *pf_els_cIso;
   std::vector<float>   *pf_els_tIso;
   std::vector<float>   *pf_els_ecalIso;
   std::vector<float>   *pf_els_hcalIso;
   std::vector<float>   *pf_els_chargedHadronIso;
   std::vector<float>   *pf_els_photonIso;
   std::vector<float>   *pf_els_neutralHadronIso;
   std::vector<float>   *pf_els_charge;
   std::vector<float>   *pf_els_hadOverEm;
   std::vector<float>   *pf_els_eOverPIn;
   std::vector<float>   *pf_els_sigmaIEtaIEta;
   std::vector<float>   *pf_els_scEnergy;
   std::vector<float>   *pf_els_scEta;
   std::vector<float>   *pf_els_scE1x5;
   std::vector<float>   *pf_els_scE2x5Max;
   std::vector<float>   *pf_els_scE5x5;
   std::vector<float>   *pf_els_isEB;
   std::vector<float>   *pf_els_isEE;
   std::vector<float>   *pf_els_dEtaIn;
   std::vector<float>   *pf_els_dPhiIn;
   std::vector<float>   *pf_els_dEtaOut;
   std::vector<float>   *pf_els_dPhiOut;
   std::vector<float>   *pf_els_numlosthits;
   std::vector<float>   *pf_els_tk_phi;
   std::vector<float>   *pf_els_d0dum;
   std::vector<float>   *pf_els_vx;
   std::vector<float>   *pf_els_vy;
   std::vector<float>   *pf_els_vz;
   std::vector<float>   *pf_els_ptError;
   std::vector<float>   *pf_els_n_inner_layer;
   std::vector<float>   *pf_els_dr03EcalRecHitSumEt;
   std::vector<float>   *pf_els_dr03HcalTowerSumEt;
   std::vector<float>   *pf_els_dr03HcalDepth1TowerSumEt;
   std::vector<float>   *pf_els_dr03HcalDepth2TowerSumEt;
   std::vector<float>   *pf_els_dr03TkSumPt;
   std::vector<float>   *pf_mus_energy;
   std::vector<float>   *pf_mus_et;
   std::vector<float>   *pf_mus_eta;
   std::vector<float>   *pf_mus_phi;
   std::vector<float>   *pf_mus_pt;
   std::vector<float>   *pf_mus_px;
   std::vector<float>   *pf_mus_py;
   std::vector<float>   *pf_mus_pz;
   std::vector<float>   *pf_mus_cIso;
   std::vector<float>   *pf_mus_tIso;
   std::vector<float>   *pf_mus_ecalIso;
   std::vector<float>   *pf_mus_hcalIso;
   std::vector<float>   *pf_mus_neutralHadronIso;
   std::vector<float>   *pf_mus_chargedHadronIso;
   std::vector<float>   *pf_mus_photonIso;
   std::vector<float>   *pf_mus_charge;
   std::vector<float>   *pf_mus_cm_chi2;
   std::vector<float>   *pf_mus_cm_ndof;
   std::vector<float>   *pf_mus_cm_pt;
   std::vector<float>   *pf_mus_cm_ptErr;
   std::vector<float>   *pf_mus_tk_chi2;
   std::vector<float>   *pf_mus_tk_ndof;
   std::vector<float>   *pf_mus_tk_pt;
   std::vector<float>   *pf_mus_tk_phi;
   std::vector<float>   *pf_mus_tk_d0dum;
   std::vector<float>   *pf_mus_tk_vz;
   std::vector<float>   *pf_mus_tk_numvalhits;
   std::vector<float>   *pf_mus_tk_ptErr;
   std::vector<float>   *pf_mus_tk_numvalPixelhits;
   std::vector<float>   *pf_mus_tk_numpixelWthMeasr;
   std::vector<float>   *pf_mus_stamu_pt;
   std::vector<float>   *pf_mus_stamu_ptErr;
   std::vector<float>   *pf_mus_num_matches;
   std::vector<float>   *pf_mus_isTrackerMuon;
   std::vector<float>   *pf_mus_isGlobalMuon;
   std::vector<float>   *pf_mus_id_GlobalMuonPromptTight;
   std::vector<float>   *pfcand_pdgId;
   std::vector<float>   *pfcand_particleId;
   std::vector<float>   *pfcand_pt;
   std::vector<float>   *pfcand_pz;
   std::vector<float>   *pfcand_px;
   std::vector<float>   *pfcand_py;
   std::vector<float>   *pfcand_eta;
   std::vector<float>   *pfcand_phi;
   std::vector<float>   *pfcand_theta;
   std::vector<float>   *pfcand_energy;
   std::vector<float>   *pfcand_charge;
   std::vector<float>   *pfmets_et;
   std::vector<float>   *pfmets_phi;
   std::vector<float>   *pfmets_ex;
   std::vector<float>   *pfmets_ey;
   std::vector<float>   *pfmets_gen_et;
   std::vector<float>   *pfmets_gen_phi;
   std::vector<float>   *pfmets_sumEt;
   std::vector<float>   *pv_x;
   std::vector<float>   *pv_y;
   std::vector<float>   *pv_z;
   std::vector<float>   *pv_ndof;
   std::vector<float>   *pv_isFake;
   std::vector<float>   *pv_tracksSize;
   UInt_t          run;
   UInt_t          event;
   UInt_t          lumiblock;
   UInt_t          bunchCrossing;
   Float_t         weight;
   std::string          *model_params;

   // List of branches
   TBranch        *b_beamSpot_x;   //!
   TBranch        *b_beamSpot_y;   //!
   TBranch        *b_els_energy;   //!
   TBranch        *b_els_et;   //!
   TBranch        *b_els_eta;   //!
   TBranch        *b_els_phi;   //!
   TBranch        *b_els_pt;   //!
   TBranch        *b_els_px;   //!
   TBranch        *b_els_py;   //!
   TBranch        *b_els_pz;   //!
   TBranch        *b_els_robustTightId;   //!
   TBranch        *b_els_simpleEleId95relIso;   //!
   TBranch        *b_els_simpleEleId90relIso;   //!
   TBranch        *b_els_simpleEleId85relIso;   //!
   TBranch        *b_els_simpleEleId80relIso;   //!
   TBranch        *b_els_simpleEleId70relIso;   //!
   TBranch        *b_els_simpleEleId95cIso;   //!
   TBranch        *b_els_simpleEleId90cIso;   //!
   TBranch        *b_els_simpleEleId85cIso;   //!
   TBranch        *b_els_simpleEleId80cIso;   //!
   TBranch        *b_els_simpleEleId70cIso;   //!
   TBranch        *b_els_cIso;   //!
   TBranch        *b_els_tIso;   //!
   TBranch        *b_els_ecalIso;   //!
   TBranch        *b_els_hcalIso;   //!
   TBranch        *b_els_charge;   //!
   TBranch        *b_els_hadOverEm;   //!
   TBranch        *b_els_eOverPIn;   //!
   TBranch        *b_els_sigmaIEtaIEta;   //!
   TBranch        *b_els_scEnergy;   //!
   TBranch        *b_els_scEta;   //!
   TBranch        *b_els_scE1x5;   //!
   TBranch        *b_els_scE2x5Max;   //!
   TBranch        *b_els_scE5x5;   //!
   TBranch        *b_els_isEB;   //!
   TBranch        *b_els_isEE;   //!
   TBranch        *b_els_dEtaIn;   //!
   TBranch        *b_els_dPhiIn;   //!
   TBranch        *b_els_dEtaOut;   //!
   TBranch        *b_els_dPhiOut;   //!
   TBranch        *b_els_numlosthits;   //!
   TBranch        *b_els_tk_phi;   //!
   TBranch        *b_els_d0dum;   //!
   TBranch        *b_els_vx;   //!
   TBranch        *b_els_vy;   //!
   TBranch        *b_els_vz;   //!
   TBranch        *b_els_ptError;   //!
   TBranch        *b_els_n_inner_layer;   //!
   TBranch        *b_els_dr03EcalRecHitSumEt;   //!
   TBranch        *b_els_dr03HcalTowerSumEt;   //!
   TBranch        *b_els_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_els_dr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_els_dr03TkSumPt;   //!
   TBranch        *b_jets_AK5PF_phi;   //!
   TBranch        *b_jets_AK5PF_pt;   //!
   TBranch        *b_jets_AK5PF_pz;   //!
   TBranch        *b_jets_AK5PF_px;   //!
   TBranch        *b_jets_AK5PF_py;   //!
   TBranch        *b_jets_AK5PF_eta;   //!
   TBranch        *b_jets_AK5PF_et;   //!
   TBranch        *b_jets_AK5PF_energy;   //!
   TBranch        *b_jets_AK5PF_parton_Id;   //!
   TBranch        *b_jets_AK5PF_parton_motherId;   //!
   TBranch        *b_jets_AK5PF_gen_pt;   //!
   TBranch        *b_jets_AK5PF_gen_phi;   //!
   TBranch        *b_jets_AK5PF_partonFlavour;   //!
   TBranch        *b_jets_AK5PF_btag_TC_highPur;   //!
   TBranch        *b_jets_AK5PF_btag_TC_highEff;   //!
   TBranch        *b_jets_AK5PF_btag_jetProb;   //!
   TBranch        *b_jets_AK5PF_btag_jetBProb;   //!
   TBranch        *b_jets_AK5PF_btag_secVertexHighPur;   //!
   TBranch        *b_jets_AK5PF_btag_secVertexHighEff;   //!
   TBranch        *b_jets_AK5PF_btag_secVertexCombined;   //!
   TBranch        *b_jets_AK5PF_jetCharge;   //!
   TBranch        *b_jets_AK5PF_chgEmE;   //!
   TBranch        *b_jets_AK5PF_chgHadE;   //!
   TBranch        *b_jets_AK5PF_photonEnergy;   //!
   TBranch        *b_jets_AK5PF_chg_Mult;   //!
   TBranch        *b_jets_AK5PF_neutralEmE;   //!
   TBranch        *b_jets_AK5PF_neutralHadE;   //!
   TBranch        *b_jets_AK5PF_neutral_Mult;   //!
   TBranch        *b_jets_AK5PF_mu_Mult;   //!
   TBranch        *b_jets_AK5PF_ehf;   //!
   TBranch        *b_jets_AK5PF_corrFactorRaw;   //!
   TBranch        *b_jets_AK5PFclean_phi;   //!
   TBranch        *b_jets_AK5PFclean_pt;   //!
   TBranch        *b_jets_AK5PFclean_pz;   //!
   TBranch        *b_jets_AK5PFclean_px;   //!
   TBranch        *b_jets_AK5PFclean_py;   //!
   TBranch        *b_jets_AK5PFclean_eta;   //!
   TBranch        *b_jets_AK5PFclean_et;   //!
   TBranch        *b_jets_AK5PFclean_energy;   //!
   TBranch        *b_jets_AK5PFclean_parton_Id;   //!
   TBranch        *b_jets_AK5PFclean_parton_motherId;   //!
   TBranch        *b_jets_AK5PFclean_gen_pt;   //!
   TBranch        *b_jets_AK5PFclean_gen_phi;   //!
   TBranch        *b_jets_AK5PFclean_partonFlavour;   //!
   TBranch        *b_jets_AK5PFclean_btag_TC_highPur;   //!
   TBranch        *b_jets_AK5PFclean_btag_TC_highEff;   //!
   TBranch        *b_jets_AK5PFclean_btag_jetProb;   //!
   TBranch        *b_jets_AK5PFclean_btag_jetBProb;   //!
   TBranch        *b_jets_AK5PFclean_btag_secVertexHighPur;   //!
   TBranch        *b_jets_AK5PFclean_btag_secVertexHighEff;   //!
   TBranch        *b_jets_AK5PFclean_btag_secVertexCombined;   //!
   TBranch        *b_jets_AK5PFclean_jetCharge;   //!
   TBranch        *b_jets_AK5PFclean_chgEmE;   //!
   TBranch        *b_jets_AK5PFclean_chgHadE;   //!
   TBranch        *b_jets_AK5PFclean_photonEnergy;   //!
   TBranch        *b_jets_AK5PFclean_chg_Mult;   //!
   TBranch        *b_jets_AK5PFclean_neutralEmE;   //!
   TBranch        *b_jets_AK5PFclean_neutralHadE;   //!
   TBranch        *b_jets_AK5PFclean_neutral_Mult;   //!
   TBranch        *b_jets_AK5PFclean_mu_Mult;   //!
   TBranch        *b_jets_AK5PFclean_ehf;   //!
   TBranch        *b_jets_AK5PFclean_corrFactorRaw;   //!
   TBranch        *b_mc_doc_id;   //!
   TBranch        *b_mc_doc_pt;   //!
   TBranch        *b_mc_doc_px;   //!
   TBranch        *b_mc_doc_py;   //!
   TBranch        *b_mc_doc_pz;   //!
   TBranch        *b_mc_doc_eta;   //!
   TBranch        *b_mc_doc_phi;   //!
   TBranch        *b_mc_doc_mother_id;   //!
   TBranch        *b_mc_doc_grandmother_id;   //!
   TBranch        *b_mc_doc_ggrandmother_id;   //!
   TBranch        *b_mc_doc_mass;   //!
   TBranch        *b_mc_electrons_pt;   //!
   TBranch        *b_mc_electrons_px;   //!
   TBranch        *b_mc_electrons_py;   //!
   TBranch        *b_mc_electrons_pz;   //!
   TBranch        *b_mc_electrons_eta;   //!
   TBranch        *b_mc_electrons_phi;   //!
   TBranch        *b_mc_electrons_energy;   //!
   TBranch        *b_mc_electrons_mother_id;   //!
   TBranch        *b_mc_electrons_grandmother_id;   //!
   TBranch        *b_mc_mus_pt;   //!
   TBranch        *b_mc_mus_px;   //!
   TBranch        *b_mc_mus_py;   //!
   TBranch        *b_mc_mus_pz;   //!
   TBranch        *b_mc_mus_eta;   //!
   TBranch        *b_mc_mus_phi;   //!
   TBranch        *b_mc_mus_energy;   //!
   TBranch        *b_mc_mus_mother_id;   //!
   TBranch        *b_mc_mus_grandmother_id;   //!
   TBranch        *b_mc_nues_pt;   //!
   TBranch        *b_mc_nues_px;   //!
   TBranch        *b_mc_nues_py;   //!
   TBranch        *b_mc_nues_pz;   //!
   TBranch        *b_mc_nues_eta;   //!
   TBranch        *b_mc_nues_phi;   //!
   TBranch        *b_mc_nues_energy;   //!
   TBranch        *b_mc_nues_mother_id;   //!
   TBranch        *b_mc_nues_grandmother_id;   //!
   TBranch        *b_mc_numus_pt;   //!
   TBranch        *b_mc_numus_px;   //!
   TBranch        *b_mc_numus_py;   //!
   TBranch        *b_mc_numus_pz;   //!
   TBranch        *b_mc_numus_eta;   //!
   TBranch        *b_mc_numus_phi;   //!
   TBranch        *b_mc_numus_energy;   //!
   TBranch        *b_mc_numus_mother_id;   //!
   TBranch        *b_mc_numus_grandmother_id;   //!
   TBranch        *b_mc_nutaus_pt;   //!
   TBranch        *b_mc_nutaus_px;   //!
   TBranch        *b_mc_nutaus_py;   //!
   TBranch        *b_mc_nutaus_pz;   //!
   TBranch        *b_mc_nutaus_eta;   //!
   TBranch        *b_mc_nutaus_phi;   //!
   TBranch        *b_mc_nutaus_energy;   //!
   TBranch        *b_mc_nutaus_mother_id;   //!
   TBranch        *b_mc_nutaus_grandmother_id;   //!
   TBranch        *b_mets_AK5_et;   //!
   TBranch        *b_mets_AK5_phi;   //!
   TBranch        *b_mus_energy;   //!
   TBranch        *b_mus_et;   //!
   TBranch        *b_mus_eta;   //!
   TBranch        *b_mus_phi;   //!
   TBranch        *b_mus_pt;   //!
   TBranch        *b_mus_px;   //!
   TBranch        *b_mus_py;   //!
   TBranch        *b_mus_pz;   //!
   TBranch        *b_mus_numberOfMatchedStations;   //!
   TBranch        *b_mus_cIso;   //!
   TBranch        *b_mus_tIso;   //!
   TBranch        *b_mus_ecalIso;   //!
   TBranch        *b_mus_hcalIso;   //!
   TBranch        *b_mus_ecalvetoDep;   //!
   TBranch        *b_mus_hcalvetoDep;   //!
   TBranch        *b_mus_pfIsolationR03_sumChargedHadronPt;   //!
   TBranch        *b_mus_pfIsolationR03_sumChargedParticlePt;   //!
   TBranch        *b_mus_pfIsolationR03_sumNeutralHadronEt;   //!
   TBranch        *b_mus_pfIsolationR03_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_mus_pfIsolationR03_sumPhotonEt;   //!
   TBranch        *b_mus_pfIsolationR03_sumPhotonEtHighThreshold;   //!
   TBranch        *b_mus_pfIsolationR03_sumPUPt;   //!
   TBranch        *b_mus_pfIsolationR04_sumChargedHadronPt;   //!
   TBranch        *b_mus_pfIsolationR04_sumChargedParticlePt;   //!
   TBranch        *b_mus_pfIsolationR04_sumNeutralHadronEt;   //!
   TBranch        *b_mus_pfIsolationR04_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_mus_pfIsolationR04_sumPhotonEt;   //!
   TBranch        *b_mus_pfIsolationR04_sumPhotonEtHighThreshold;   //!
   TBranch        *b_mus_pfIsolationR04_sumPUPt;   //!
   TBranch        *b_mus_charge;   //!
   TBranch        *b_mus_cm_chi2;   //!
   TBranch        *b_mus_cm_ndof;   //!
   TBranch        *b_mus_cm_pt;   //!
   TBranch        *b_mus_cm_ptErr;   //!
   TBranch        *b_mus_tk_chi2;   //!
   TBranch        *b_mus_tk_ndof;   //!
   TBranch        *b_mus_tk_pt;   //!
   TBranch        *b_mus_tk_phi;   //!
   TBranch        *b_mus_tk_d0dum;   //!
   TBranch        *b_mus_tk_vz;   //!
   TBranch        *b_mus_tk_numvalhits;   //!
   TBranch        *b_mus_tk_ptErr;   //!
   TBranch        *b_mus_tk_numvalPixelhits;   //!
   TBranch        *b_mus_tk_numpixelWthMeasr;   //!
   TBranch        *b_mus_stamu_pt;   //!
   TBranch        *b_mus_stamu_ptErr;   //!
   TBranch        *b_mus_num_matches;   //!
   TBranch        *b_mus_isPFMuon;   //!
   TBranch        *b_mus_isTrackerMuon;   //!
   TBranch        *b_mus_isGlobalMuon;   //!
   TBranch        *b_mus_id_AllGlobalMuons;   //!
   TBranch        *b_mus_id_AllTrackerMuons;   //!
   TBranch        *b_mus_id_GlobalMuonPromptTight;   //!
   TBranch        *b_pfTypeImets_et;   //!
   TBranch        *b_pfTypeImets_phi;   //!
   TBranch        *b_pfTypeImets_ex;   //!
   TBranch        *b_pfTypeImets_ey;   //!
   TBranch        *b_pfTypeImets_gen_et;   //!
   TBranch        *b_pfTypeImets_gen_phi;   //!
   TBranch        *b_pfTypeImets_sumEt;   //!
   TBranch        *b_pf_els_energy;   //!
   TBranch        *b_pf_els_et;   //!
   TBranch        *b_pf_els_eta;   //!
   TBranch        *b_pf_els_phi;   //!
   TBranch        *b_pf_els_pt;   //!
   TBranch        *b_pf_els_px;   //!
   TBranch        *b_pf_els_py;   //!
   TBranch        *b_pf_els_pz;   //!
   TBranch        *b_pf_els_robustTightId;   //!
   TBranch        *b_pf_els_simpleEleId95relIso;   //!
   TBranch        *b_pf_els_simpleEleId90relIso;   //!
   TBranch        *b_pf_els_simpleEleId85relIso;   //!
   TBranch        *b_pf_els_simpleEleId80relIso;   //!
   TBranch        *b_pf_els_simpleEleId70relIso;   //!
   TBranch        *b_pf_els_simpleEleId95cIso;   //!
   TBranch        *b_pf_els_simpleEleId90cIso;   //!
   TBranch        *b_pf_els_simpleEleId85cIso;   //!
   TBranch        *b_pf_els_simpleEleId80cIso;   //!
   TBranch        *b_pf_els_simpleEleId70cIso;   //!
   TBranch        *b_pf_els_cIso;   //!
   TBranch        *b_pf_els_tIso;   //!
   TBranch        *b_pf_els_ecalIso;   //!
   TBranch        *b_pf_els_hcalIso;   //!
   TBranch        *b_pf_els_chargedHadronIso;   //!
   TBranch        *b_pf_els_photonIso;   //!
   TBranch        *b_pf_els_neutralHadronIso;   //!
   TBranch        *b_pf_els_charge;   //!
   TBranch        *b_pf_els_hadOverEm;   //!
   TBranch        *b_pf_els_eOverPIn;   //!
   TBranch        *b_pf_els_sigmaIEtaIEta;   //!
   TBranch        *b_pf_els_scEnergy;   //!
   TBranch        *b_pf_els_scEta;   //!
   TBranch        *b_pf_els_scE1x5;   //!
   TBranch        *b_pf_els_scE2x5Max;   //!
   TBranch        *b_pf_els_scE5x5;   //!
   TBranch        *b_pf_els_isEB;   //!
   TBranch        *b_pf_els_isEE;   //!
   TBranch        *b_pf_els_dEtaIn;   //!
   TBranch        *b_pf_els_dPhiIn;   //!
   TBranch        *b_pf_els_dEtaOut;   //!
   TBranch        *b_pf_els_dPhiOut;   //!
   TBranch        *b_pf_els_numlosthits;   //!
   TBranch        *b_pf_els_tk_phi;   //!
   TBranch        *b_pf_els_d0dum;   //!
   TBranch        *b_pf_els_vx;   //!
   TBranch        *b_pf_els_vy;   //!
   TBranch        *b_pf_els_vz;   //!
   TBranch        *b_pf_els_ptError;   //!
   TBranch        *b_pf_els_n_inner_layer;   //!
   TBranch        *b_pf_els_dr03EcalRecHitSumEt;   //!
   TBranch        *b_pf_els_dr03HcalTowerSumEt;   //!
   TBranch        *b_pf_els_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_pf_els_dr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_pf_els_dr03TkSumPt;   //!
   TBranch        *b_pf_mus_energy;   //!
   TBranch        *b_pf_mus_et;   //!
   TBranch        *b_pf_mus_eta;   //!
   TBranch        *b_pf_mus_phi;   //!
   TBranch        *b_pf_mus_pt;   //!
   TBranch        *b_pf_mus_px;   //!
   TBranch        *b_pf_mus_py;   //!
   TBranch        *b_pf_mus_pz;   //!
   TBranch        *b_pf_mus_cIso;   //!
   TBranch        *b_pf_mus_tIso;   //!
   TBranch        *b_pf_mus_ecalIso;   //!
   TBranch        *b_pf_mus_hcalIso;   //!
   TBranch        *b_pf_mus_neutralHadronIso;   //!
   TBranch        *b_pf_mus_chargedHadronIso;   //!
   TBranch        *b_pf_mus_photonIso;   //!
   TBranch        *b_pf_mus_charge;   //!
   TBranch        *b_pf_mus_cm_chi2;   //!
   TBranch        *b_pf_mus_cm_ndof;   //!
   TBranch        *b_pf_mus_cm_pt;   //!
   TBranch        *b_pf_mus_cm_ptErr;   //!
   TBranch        *b_pf_mus_tk_chi2;   //!
   TBranch        *b_pf_mus_tk_ndof;   //!
   TBranch        *b_pf_mus_tk_pt;   //!
   TBranch        *b_pf_mus_tk_phi;   //!
   TBranch        *b_pf_mus_tk_d0dum;   //!
   TBranch        *b_pf_mus_tk_vz;   //!
   TBranch        *b_pf_mus_tk_numvalhits;   //!
   TBranch        *b_pf_mus_tk_ptErr;   //!
   TBranch        *b_pf_mus_tk_numvalPixelhits;   //!
   TBranch        *b_pf_mus_tk_numpixelWthMeasr;   //!
   TBranch        *b_pf_mus_stamu_pt;   //!
   TBranch        *b_pf_mus_stamu_ptErr;   //!
   TBranch        *b_pf_mus_num_matches;   //!
   TBranch        *b_pf_mus_isTrackerMuon;   //!
   TBranch        *b_pf_mus_isGlobalMuon;   //!
   TBranch        *b_pf_mus_id_GlobalMuonPromptTight;   //!
   TBranch        *b_pfcand_pdgId;   //!
   TBranch        *b_pfcand_particleId;   //!
   TBranch        *b_pfcand_pt;   //!
   TBranch        *b_pfcand_pz;   //!
   TBranch        *b_pfcand_px;   //!
   TBranch        *b_pfcand_py;   //!
   TBranch        *b_pfcand_eta;   //!
   TBranch        *b_pfcand_phi;   //!
   TBranch        *b_pfcand_theta;   //!
   TBranch        *b_pfcand_energy;   //!
   TBranch        *b_pfcand_charge;   //!
   TBranch        *b_pfmets_et;   //!
   TBranch        *b_pfmets_phi;   //!
   TBranch        *b_pfmets_ex;   //!
   TBranch        *b_pfmets_ey;   //!
   TBranch        *b_pfmets_gen_et;   //!
   TBranch        *b_pfmets_gen_phi;   //!
   TBranch        *b_pfmets_sumEt;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_ndof;   //!
   TBranch        *b_pv_isFake;   //!
   TBranch        *b_pv_tracksSize;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_model_params;   //!


   TChain* chainB;
   TChain* chainV;
   TChain* chainA;
   void InitializeA(TChain *fChain);
   void InitializeB(TChain *fChain);
   void InitializeV(TChain *fChain);

  // ==== END giant copy/paste of cfA variable ================================================================



};

#endif
