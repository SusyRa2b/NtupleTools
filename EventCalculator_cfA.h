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

//#include "BasicLoopCU.h" //get all of the tree-related variables



class TRandom3;

//global constants
const double mW_ = 80.399;
const double mtop_ = 172.0;
const double lumi_ = 1.; //fix to 1/pb and scale MC later (e.g. in drawReducedTrees)

//define the BNN functions
/*
double HT260_MHT60_v2_HT260_v2_160431_161204_sel1(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT250_MHT60_v2_HT250_v2_161205_163268_sel1(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT250_MHT60_v3_HT250_v3_163269_164923_sel1(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT250_MHT70_v1_HT250_v4_164924_165921_sel1(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT300_MHT75_v7_HT300_v6_165922_166300_166374_166978_sel1(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT300_MHT75_v8_HT300_v7_166301_166373_sel1(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT300_MHT80_v1_HT300_v6_166979_170064_sel1(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT300_MHT80_v2_HT300_v9_170065_173211_sel1(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT300_MHT90_v2_HT300_v9_173212_176544_sel1(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT350_MHT90_v1_HT350_v8_176545_178410_sel1(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT350_MHT110_v3_HT350_v11_178411_180252_sel1(std::vector<double>& inputvars, int first=0, int last=100-1);

double HT260_MHT60_v2_HT260_v2_160431_161204_sel2(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT250_MHT60_v2_HT250_v2_161205_163268_sel2(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT250_MHT60_v3_HT250_v3_163269_164923_sel2(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT250_MHT70_v1_HT250_v4_164924_165921_sel2(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT300_MHT75_v7_HT300_v6_165922_166300_166374_166978_sel2(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT300_MHT75_v8_HT300_v7_166301_166373_sel2(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT300_MHT80_v1_HT300_v6_166979_170064_sel2(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT300_MHT80_v2_HT300_v9_170065_173211_sel2(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT300_MHT90_v2_HT300_v9_173212_176544_sel2(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT350_MHT90_v1_HT350_v8_176545_178410_sel2(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT350_MHT110_v3_HT350_v11_178411_180252_sel2(std::vector<double>& inputvars, int first=0, int last=100-1);

double ElHad_300_75(std::vector<double>& inputvars, int first=0, int last=100-1);
double MuHad_300_75(std::vector<double>& inputvars, int first=0, int last=100-1);
double ElHad_300_80(std::vector<double>& inputvars, int first=0, int last=100-1);
double MuHad_300_80(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT300MHT90elNew(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT300MHT90muNew(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT350MHT90elNew(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT350MHT90muNew(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT350MHT110elNew(std::vector<double>& inputvars, int first=0, int last=100-1);
double HT350MHT110muNew(std::vector<double>& inputvars, int first=0, int last=100-1);
*/

class cellECAL {
public:
  cellECAL(double,double,int);
  double eta;
  double phi;
  int status;
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
  float getPUWeight(/*reweight::LumiReWeighting lumiWeights*/);
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

  //  bool passBEFilter(); //TODO migrate this code
  bool passHLT();
  bool passUtilityHLT(int &version, int &prescale);
  bool passUtilityPrescaleModuleHLT();
  unsigned int utilityHLT_HT300_CentralJet30_BTagIP();
  float getHLTHTeff(float offHT);
  float getHLTMHTeff(float offMET, float offHT, uint nElectrons, uint nMuons, double mindphin);
  //  double getPFMHTWeight(); //TODO
  double getHLTMHTeffBNN(float offMET, float offHT, uint nElectrons, uint nMuons, double mindphin, 
			 double& effUp, double& effDown);

  int countGoodPV();
  bool passPV() ;

  float getHT();
  float getST();

  float getMET();
  float getMETphi();

  float getMHT();
  float getMHTphi();
  float getMHTphi(int ignoredJet);

  //  void getTransverseThrustVariables(float & thrust, float & thrustPhi, bool addMET);
  void getSphericityJetMET(float & lambda1, float & lambda2, float & det,const int jetmax, bool addMET);


  bool passCut(const TString cutTag);
  bool setCutScheme();
  void setIgnoredCut(const TString cutTag);
  void setRequiredCut(const TString cutTag);
  void resetIgnoredCut();
  void resetRequiredCut();
  bool cutRequired(const TString cutTag);
  int Cut();

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
  float getBTagIPWeight(); //this function should be called *after* offline tagging 
  float jetTagEff(unsigned int ijet, TH1F* h_btageff, TH1F* h_ctageff, TH1F* h_ltageff,
		  const float extraSFb, const float extraSFc, const float extraSFl,const BTagEffModifier modifier=kBTagModifier0);

  //other stuff
  double getDeltaPhi(double a, double b);
  bool isSampleRealData();
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
//   std::vector<jet2_s> * myJetsPF;
//   std::vector<jethelper2_s> * myJetsPFhelper;

//   std::vector<electron1_s> * myElectronsPF;
//   std::vector<electronhelper1_s> * myElectronsPFhelper;
//   std::vector<muon1_s> * myMuonsPF;
//   std::vector<muon_s> * myMuonsRECO;
//   std::vector<muonhelper1_s> * myMuonsPFhelper;
//   std::vector<muonhelper_s> * myMuonsRECOhelper;
//   std::vector<tau_s> * myTausPF;
//   std::vector<met1_s> * myMETPF;
//   std::vector<met4_s> * myMETPFType1;
//   std::vector<met_s> * myMETcalo;
//   std::vector<vertex_s> * myVertex;
//   std::vector<genparticlehelperra2_s> * myGenParticles;
//   double* myGenWeight;

  //versions for jet-loss study
//   std::vector<jet2_s> * myJetsPF_temp;
//   std::vector<met1_s> * myMETPF_temp;

//   double* myEDM_bunchCrossing;
//   double* myEDM_event;
//   double* myEDM_isRealData;
//   double* myEDM_luminosityBlock;
//   double* myEDM_run;

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
  void changeVariablesGenSmearing(TRandom3* random);
  //void resetVariables(); //FIXME CFA
  bool  recalculatedVariables_;

  TStopwatch* watch_; //not used for everyday running; just for testing

  std::vector<TString> cutTags_;
  std::map<TString, TString> cutNames_; //key is a cutTag. these should be "human readable" but should not have any spaces
  std::vector<TString> ignoredCut_; //allow more than 1 ignored cut!
  std::vector<TString> requiredCut_; //new feature to *turn on* a cut that is usually not required by a given cut scheme


  // ==== BEGIN giant copy/paste of cfA variables =============================================================
   std::vector<float>   *trigger_prescalevalue;
   std::vector<std::string>  *trigger_name;
   std::vector<float>   *trigger_decision;
   std::vector<std::string>  *trigger_lastfiltername;
   std::vector<std::vector<float> > *triggerobject_pt;
   std::vector<std::vector<float> > *triggerobject_px;
   std::vector<std::vector<float> > *triggerobject_py;
   std::vector<std::vector<float> > *triggerobject_pz;
   std::vector<std::vector<float> > *triggerobject_et;
   std::vector<std::vector<float> > *triggerobject_energy;
   std::vector<std::vector<float> > *triggerobject_phi;
   std::vector<std::vector<float> > *triggerobject_eta;
   std::vector<std::vector<TString> > *triggerobject_collectionname;
   std::vector<float>   *L1trigger_bit;
   std::vector<float>   *L1trigger_techTrigger;
   std::vector<float>   *L1trigger_prescalevalue;
   std::vector<std::string>  *L1trigger_name;
   std::vector<std::string>  *L1trigger_alias;
   std::vector<float>   *L1trigger_decision;
   std::vector<float>   *L1trigger_decision_nomask;
   std::vector<float>   *els_conversion_dist;
   std::vector<float>   *els_conversion_dcot;
   Int_t           hbhefilter_decision;
   Int_t           trackingfailurefilter_decision;
   Int_t           cschalofilter_decision;
   Int_t           ecalTPfilter_decision;
   Int_t           scrapingVeto_decision;
   Int_t           inconsistentPFmuonfilter_decision;
   Int_t           greedymuonfilter_decision;
   Int_t           eenoisefilter_decision;
 
   // Declaration of leave types
   UInt_t          NbeamSpot;
   std::vector<float>   *beamSpot_x;
   std::vector<float>   *beamSpot_y;
   std::vector<float>   *beamSpot_z;
   std::vector<float>   *beamSpot_x0Error;
   std::vector<float>   *beamSpot_y0Error;
   std::vector<float>   *beamSpot_z0Error;
   std::vector<float>   *beamSpot_sigmaZ;
   std::vector<float>   *beamSpot_sigmaZ0Error;
   std::vector<float>   *beamSpot_dxdz;
   std::vector<float>   *beamSpot_dxdzError;
   std::vector<float>   *beamSpot_dydz;
   std::vector<float>   *beamSpot_dydzError;
   std::vector<float>   *beamSpot_beamWidthX;
   std::vector<float>   *beamSpot_beamWidthY;
   std::vector<float>   *beamSpot_beamWidthXError;
   std::vector<float>   *beamSpot_beamWidthYError;
   UInt_t          Nels;
   std::vector<float>   *pf_els_energy;
   std::vector<float>   *pf_els_et;
   std::vector<float>   *pf_els_eta;
   std::vector<float>   *pf_els_phi;
   std::vector<float>   *pf_els_pt;
   std::vector<float>   *pf_els_px;
   std::vector<float>   *pf_els_py;
   std::vector<float>   *pf_els_pz;
   std::vector<float>   *pf_els_status;
   std::vector<float>   *pf_els_theta;
   std::vector<float>   *els_pt;
   std::vector<float>   *els_eta;
   std::vector<float>   *els_phi;
   std::vector<float>   *els_charge;
   std::vector<float>   *els_gen_id;
   std::vector<float>   *els_gen_phi;
   std::vector<float>   *els_gen_pt;
   std::vector<float>   *els_gen_pz;
   std::vector<float>   *els_gen_px;
   std::vector<float>   *els_gen_py;
   std::vector<float>   *els_gen_eta;
   std::vector<float>   *els_gen_theta;
   std::vector<float>   *els_gen_et;
   std::vector<float>   *els_gen_mother_id;
   std::vector<float>   *els_gen_mother_phi;
   std::vector<float>   *els_gen_mother_pt;
   std::vector<float>   *els_gen_mother_pz;
   std::vector<float>   *els_gen_mother_px;
   std::vector<float>   *els_gen_mother_py;
   std::vector<float>   *els_gen_mother_eta;
   std::vector<float>   *els_gen_mother_theta;
   std::vector<float>   *els_gen_mother_et;
   std::vector<float>   *pf_els_tightId;
   std::vector<float>   *pf_els_looseId;
   std::vector<float>   *pf_els_robustTightId;
   std::vector<float>   *pf_els_robustLooseId;
   std::vector<float>   *pf_els_robustHighEnergyId;
   std::vector<float>   *pf_els_simpleEleId95relIso;
   std::vector<float>   *pf_els_simpleEleId90relIso;
   std::vector<float>   *pf_els_simpleEleId85relIso;
   std::vector<float>   *pf_els_simpleEleId80relIso;
   std::vector<float>   *pf_els_simpleEleId70relIso;
   std::vector<float>   *pf_els_simpleEleId60relIso;
   std::vector<float>   *pf_els_simpleEleId95cIso;
   std::vector<float>   *pf_els_simpleEleId90cIso;
   std::vector<float>   *pf_els_simpleEleId85cIso;
   std::vector<float>   *pf_els_simpleEleId80cIso;
   std::vector<float>   *pf_els_simpleEleId70cIso;
   std::vector<float>   *pf_els_simpleEleId60cIso;
   std::vector<float>   *pf_els_cIso;
   std::vector<float>   *pf_els_tIso;
   std::vector<float>   *pf_els_ecalIso;
   std::vector<float>   *pf_els_hcalIso;
   std::vector<float>   *pf_els_chi2;
   std::vector<float>   *pf_els_charge;
   std::vector<float>   *pf_els_caloEnergy;
   std::vector<float>   *pf_els_hadOverEm;
   std::vector<float>   *pf_els_eOverPIn;
   std::vector<float>   *pf_els_eSeedOverPOut;
   std::vector<float>   *pf_els_sigmaEtaEta;
   std::vector<float>   *pf_els_sigmaIEtaIEta;
   std::vector<float>   *pf_els_scEnergy;
   std::vector<float>   *pf_els_scRawEnergy;
   std::vector<float>   *pf_els_scSeedEnergy;
   std::vector<float>   *pf_els_scEta;
   std::vector<float>   *pf_els_scPhi;
   std::vector<float>   *pf_els_scEtaWidth;
   std::vector<float>   *pf_els_scPhiWidth;
   std::vector<float>   *pf_els_scE1x5;
   std::vector<float>   *pf_els_scE2x5Max;
   std::vector<float>   *pf_els_scE5x5;
   std::vector<float>   *pf_els_isEB;
   std::vector<float>   *pf_els_isEE;
   std::vector<float>   *pf_els_dEtaIn;
   std::vector<float>   *pf_els_dPhiIn;
   std::vector<float>   *pf_els_dEtaOut;
   std::vector<float>   *pf_els_dPhiOut;
   std::vector<float>   *pf_els_numvalhits;
   std::vector<float>   *pf_els_numlosthits;
   std::vector<float>   *pf_els_basicClustersSize;
   std::vector<float>   *pf_els_tk_pz;
   std::vector<float>   *pf_els_tk_pt;
   std::vector<float>   *pf_els_tk_phi;
   std::vector<float>   *pf_els_tk_eta;
   std::vector<float>   *pf_els_d0dum;
   std::vector<float>   *pf_els_dz;
   std::vector<float>   *pf_els_vx;
   std::vector<float>   *pf_els_vy;
   std::vector<float>   *pf_els_vz;
   std::vector<float>   *pf_els_ndof;
   std::vector<float>   *pf_els_ptError;
   std::vector<float>   *pf_els_d0dumError;
   std::vector<float>   *pf_els_dzError;
   std::vector<float>   *pf_els_etaError;
   std::vector<float>   *pf_els_phiError;
   std::vector<float>   *pf_els_tk_charge;
   std::vector<float>   *pf_els_core_ecalDrivenSeed;
   std::vector<float>   *pf_els_n_inner_layer;
   std::vector<float>   *pf_els_n_outer_layer;
   std::vector<float>   *pf_els_ctf_tk_id;
   std::vector<float>   *pf_els_ctf_tk_charge;
   std::vector<float>   *pf_els_ctf_tk_eta;
   std::vector<float>   *pf_els_ctf_tk_phi;
   std::vector<float>   *pf_els_fbrem;
   std::vector<float>   *pf_els_shFracInnerHits;
   std::vector<float>   *pf_els_dr03EcalRecHitSumEt;
   std::vector<float>   *pf_els_dr03HcalTowerSumEt;
   std::vector<float>   *pf_els_dr03HcalDepth1TowerSumEt;
   std::vector<float>   *pf_els_dr03HcalDepth2TowerSumEt;
   std::vector<float>   *pf_els_dr03TkSumPt;
   std::vector<float>   *pf_els_dr04EcalRecHitSumEt;
   std::vector<float>   *pf_els_dr04HcalTowerSumEt;
   std::vector<float>   *pf_els_dr04HcalDepth1TowerSumEt;
   std::vector<float>   *pf_els_dr04HcalDepth2TowerSumEt;
   std::vector<float>   *pf_els_dr04TkSumPt;
   std::vector<float>   *pf_els_chargedHadronIso;
   std::vector<float>   *pf_els_photonIso;
   std::vector<float>   *pf_els_neutralHadronIso;
   std::vector<float>   *pf_els_cpx;
   std::vector<float>   *pf_els_cpy;
   std::vector<float>   *pf_els_cpz;
   std::vector<float>   *pf_els_vpx;
   std::vector<float>   *pf_els_vpy;
   std::vector<float>   *pf_els_vpz;
   std::vector<float>   *pf_els_cx;
   std::vector<float>   *pf_els_cy;
   std::vector<float>   *pf_els_cz;
   UInt_t          NhcalNoiseSummary;
   std::vector<float>   *hcalNoiseSummary_passLooseNoiseFilter;
   std::vector<float>   *hcalNoiseSummary_passTightNoiseFilter;
   std::vector<float>   *hcalNoiseSummary_passHighLevelNoiseFilter;
   std::vector<float>   *hcalNoiseSummary_noiseFilterStatus;
   std::vector<float>   *hcalNoiseSummary_noiseType;
   std::vector<float>   *hcalNoiseSummary_eventEMEnergy;
   std::vector<float>   *hcalNoiseSummary_eventHadEnergy;
   std::vector<float>   *hcalNoiseSummary_eventTrackEnergy;
   std::vector<float>   *hcalNoiseSummary_eventEMFraction;
   std::vector<float>   *hcalNoiseSummary_eventChargeFraction;
   std::vector<float>   *hcalNoiseSummary_min10GeVHitTime;
   std::vector<float>   *hcalNoiseSummary_max10GeVHitTime;
   std::vector<float>   *hcalNoiseSummary_rms10GeVHitTime;
   std::vector<float>   *hcalNoiseSummary_min25GeVHitTime;
   std::vector<float>   *hcalNoiseSummary_max25GeVHitTime;
   std::vector<float>   *hcalNoiseSummary_rms25GeVHitTime;
   std::vector<float>   *hcalNoiseSummary_num10GeVHits;
   std::vector<float>   *hcalNoiseSummary_num25GeVHits;
   std::vector<float>   *hcalNoiseSummary_minE2TS;
   std::vector<float>   *hcalNoiseSummary_minE10TS;
   std::vector<float>   *hcalNoiseSummary_minE2Over10TS;
   std::vector<float>   *hcalNoiseSummary_maxZeros;
   std::vector<float>   *hcalNoiseSummary_maxHPDHits;
   std::vector<float>   *hcalNoiseSummary_maxRBXHits;
   std::vector<float>   *hcalNoiseSummary_minHPDEMF;
   std::vector<float>   *hcalNoiseSummary_minRBXEMF;
   std::vector<float>   *hcalNoiseSummary_numProblematicRBXs;
   std::vector<float>   *hcalNoiseSummary_maxE2Over10TS;
   std::vector<float>   *hcalNoiseSummary_maxHPDNoOtherHits;
   UInt_t          NhybridBBC;
   std::vector<float>   *hybridBBC_energy;
   std::vector<float>   *hybridBBC_x;
   std::vector<float>   *hybridBBC_y;
   std::vector<float>   *hybridBBC_z;
   std::vector<float>   *hybridBBC_rho;
   std::vector<float>   *hybridBBC_phi;
   std::vector<float>   *hybridBBC_eta;
   std::vector<float>   *hybridBBC_theta;
   UInt_t          Njets_AK5;
   std::vector<float>   *jets_AK5_energy;
   std::vector<float>   *jets_AK5_et;
   std::vector<float>   *jets_AK5_eta;
   std::vector<float>   *jets_AK5_phi;
   std::vector<float>   *jets_AK5_pt;
   std::vector<float>   *jets_AK5_px;
   std::vector<float>   *jets_AK5_py;
   std::vector<float>   *jets_AK5_pz;
   std::vector<float>   *jets_AK5_status;
   std::vector<float>   *jets_AK5_theta;
   std::vector<float>   *jets_AK5_parton_Id;
   std::vector<float>   *jets_AK5_parton_motherId;
   std::vector<float>   *jets_AK5_parton_pt;
   std::vector<float>   *jets_AK5_parton_phi;
   std::vector<float>   *jets_AK5_parton_eta;
   std::vector<float>   *jets_AK5_parton_Energy;
   std::vector<float>   *jets_AK5_parton_mass;
   std::vector<float>   *jets_AK5_gen_et;
   std::vector<float>   *jets_AK5_gen_pt;
   std::vector<float>   *jets_AK5_gen_eta;
   std::vector<float>   *jets_AK5_gen_phi;
   std::vector<float>   *jets_AK5_gen_mass;
   std::vector<float>   *jets_AK5_gen_Energy;
   std::vector<float>   *jets_AK5_gen_Id;
   std::vector<float>   *jets_AK5_gen_motherID;
   std::vector<float>   *jets_AK5_gen_threeCharge;
   std::vector<float>   *jets_AK5_partonFlavour;
   std::vector<float>   *jets_AK5_btag_TC_highPur;
   std::vector<float>   *jets_AK5_btag_TC_highEff;
   std::vector<float>   *jets_AK5_btag_jetProb;
   std::vector<float>   *jets_AK5_btag_jetBProb;
   std::vector<float>   *jets_AK5_btag_softEle;
   std::vector<float>   *jets_AK5_btag_softMuon;
   std::vector<float>   *jets_AK5_btag_secVertexHighPur;
   std::vector<float>   *jets_AK5_btag_secVertexHighEff;
   std::vector<float>   *jets_AK5_btag_secVertexCombined;
   std::vector<float>   *jets_AK5_jetCharge;
   std::vector<float>   *jets_AK5_chgEmE;
   std::vector<float>   *jets_AK5_chgHadE;
   std::vector<float>   *jets_AK5_chgMuE;
   std::vector<float>   *jets_AK5_chg_Mult;
   std::vector<float>   *jets_AK5_neutralEmE;
   std::vector<float>   *jets_AK5_neutralHadE;
   std::vector<float>   *jets_AK5_neutral_Mult;
   std::vector<float>   *jets_AK5_mu_Mult;
   std::vector<float>   *jets_AK5_emf;
   std::vector<float>   *jets_AK5_ehf;
   std::vector<float>   *jets_AK5_n60;
   std::vector<float>   *jets_AK5_n90;
   std::vector<float>   *jets_AK5_etaetaMoment;
   std::vector<float>   *jets_AK5_etaphiMoment;
   std::vector<float>   *jets_AK5_phiphiMoment;
   std::vector<float>   *jets_AK5_n90Hits;
   std::vector<float>   *jets_AK5_fHPD;
   std::vector<float>   *jets_AK5_fRBX;
   std::vector<float>   *jets_AK5_hitsInN90;
   std::vector<float>   *jets_AK5_nECALTowers;
   std::vector<float>   *jets_AK5_nHCALTowers;
   std::vector<float>   *jets_AK5_fSubDetector1;
   std::vector<float>   *jets_AK5_fSubDetector2;
   std::vector<float>   *jets_AK5_fSubDetector3;
   std::vector<float>   *jets_AK5_fSubDetector4;
   std::vector<float>   *jets_AK5_area;
   std::vector<float>   *jets_AK5_corrFactorRaw;
   std::vector<float>   *jets_AK5_mass;
   UInt_t          Njets_AK5JPT;
   std::vector<float>   *jets_AK5JPT_energy;
   std::vector<float>   *jets_AK5JPT_et;
   std::vector<float>   *jets_AK5JPT_eta;
   std::vector<float>   *jets_AK5JPT_phi;
   std::vector<float>   *jets_AK5JPT_pt;
   std::vector<float>   *jets_AK5JPT_px;
   std::vector<float>   *jets_AK5JPT_py;
   std::vector<float>   *jets_AK5JPT_pz;
   std::vector<float>   *jets_AK5JPT_status;
   std::vector<float>   *jets_AK5JPT_theta;
   std::vector<float>   *jets_AK5JPT_parton_Id;
   std::vector<float>   *jets_AK5JPT_parton_motherId;
   std::vector<float>   *jets_AK5JPT_parton_pt;
   std::vector<float>   *jets_AK5JPT_parton_phi;
   std::vector<float>   *jets_AK5JPT_parton_eta;
   std::vector<float>   *jets_AK5JPT_parton_Energy;
   std::vector<float>   *jets_AK5JPT_parton_mass;
   std::vector<float>   *jets_AK5JPT_gen_et;
   std::vector<float>   *jets_AK5JPT_gen_pt;
   std::vector<float>   *jets_AK5JPT_gen_eta;
   std::vector<float>   *jets_AK5JPT_gen_phi;
   std::vector<float>   *jets_AK5JPT_gen_mass;
   std::vector<float>   *jets_AK5JPT_gen_Energy;
   std::vector<float>   *jets_AK5JPT_gen_Id;
   std::vector<float>   *jets_AK5JPT_gen_motherID;
   std::vector<float>   *jets_AK5JPT_gen_threeCharge;
   std::vector<float>   *jets_AK5JPT_partonFlavour;
   std::vector<float>   *jets_AK5JPT_btag_TC_highPur;
   std::vector<float>   *jets_AK5JPT_btag_TC_highEff;
   std::vector<float>   *jets_AK5JPT_btag_jetProb;
   std::vector<float>   *jets_AK5JPT_btag_jetBProb;
   std::vector<float>   *jets_AK5JPT_btag_softEle;
   std::vector<float>   *jets_AK5JPT_btag_softMuon;
   std::vector<float>   *jets_AK5JPT_btag_secVertexHighPur;
   std::vector<float>   *jets_AK5JPT_btag_secVertexHighEff;
   std::vector<float>   *jets_AK5JPT_btag_secVertexCombined;
   std::vector<float>   *jets_AK5JPT_jetCharge;
   std::vector<float>   *jets_AK5JPT_chgEmE;
   std::vector<float>   *jets_AK5JPT_chgHadE;
   std::vector<float>   *jets_AK5JPT_chgMuE;
   std::vector<float>   *jets_AK5JPT_chg_Mult;
   std::vector<float>   *jets_AK5JPT_neutralEmE;
   std::vector<float>   *jets_AK5JPT_neutralHadE;
   std::vector<float>   *jets_AK5JPT_neutral_Mult;
   std::vector<float>   *jets_AK5JPT_mu_Mult;
   std::vector<float>   *jets_AK5JPT_emf;
   std::vector<float>   *jets_AK5JPT_ehf;
   std::vector<float>   *jets_AK5JPT_n60;
   std::vector<float>   *jets_AK5JPT_n90;
   std::vector<float>   *jets_AK5JPT_etaetaMoment;
   std::vector<float>   *jets_AK5JPT_etaphiMoment;
   std::vector<float>   *jets_AK5JPT_phiphiMoment;
   std::vector<float>   *jets_AK5JPT_n90Hits;
   std::vector<float>   *jets_AK5JPT_fHPD;
   std::vector<float>   *jets_AK5JPT_fRBX;
   std::vector<float>   *jets_AK5JPT_hitsInN90;
   std::vector<float>   *jets_AK5JPT_nECALTowers;
   std::vector<float>   *jets_AK5JPT_nHCALTowers;
   std::vector<float>   *jets_AK5JPT_fSubDetector1;
   std::vector<float>   *jets_AK5JPT_fSubDetector2;
   std::vector<float>   *jets_AK5JPT_fSubDetector3;
   std::vector<float>   *jets_AK5JPT_fSubDetector4;
   std::vector<float>   *jets_AK5JPT_area;
   std::vector<float>   *jets_AK5JPT_corrFactorRaw;
   std::vector<float>   *jets_AK5JPT_mass;
   UInt_t          Njets_AK5PF;
   std::vector<float>   *jets_AK5PF_energy;
   std::vector<float>   *jets_AK5PF_et;
   std::vector<float>   *jets_AK5PF_eta;
   std::vector<float>   *jets_AK5PF_phi;
   std::vector<float>   *jets_AK5PF_pt;
   std::vector<float>   *jets_AK5PF_px;
   std::vector<float>   *jets_AK5PF_py;
   std::vector<float>   *jets_AK5PF_pz;
   std::vector<float>   *jets_AK5PF_status;
   std::vector<float>   *jets_AK5PF_theta;
   std::vector<float>   *jets_AK5PF_parton_Id;
   std::vector<float>   *jets_AK5PF_parton_motherId;
   std::vector<float>   *jets_AK5PF_parton_pt;
   std::vector<float>   *jets_AK5PF_parton_phi;
   std::vector<float>   *jets_AK5PF_parton_eta;
   std::vector<float>   *jets_AK5PF_parton_Energy;
   std::vector<float>   *jets_AK5PF_parton_mass;
   std::vector<float>   *jets_AK5PF_gen_et;
   std::vector<float>   *jets_AK5PF_gen_pt;
   std::vector<float>   *jets_AK5PF_gen_eta;
   std::vector<float>   *jets_AK5PF_gen_phi;
   std::vector<float>   *jets_AK5PF_gen_mass;
   std::vector<float>   *jets_AK5PF_gen_Energy;
   std::vector<float>   *jets_AK5PF_gen_Id;
   std::vector<float>   *jets_AK5PF_gen_motherID;
   std::vector<float>   *jets_AK5PF_gen_threeCharge;
   std::vector<float>   *jets_AK5PF_partonFlavour;
   std::vector<float>   *jets_AK5PF_btag_TC_highPur;
   std::vector<float>   *jets_AK5PF_btag_TC_highEff;
   std::vector<float>   *jets_AK5PF_btag_jetProb;
   std::vector<float>   *jets_AK5PF_btag_jetBProb;
   std::vector<float>   *jets_AK5PF_btag_softEle;
   std::vector<float>   *jets_AK5PF_btag_softMuon;
   std::vector<float>   *jets_AK5PF_btag_secVertexHighPur;
   std::vector<float>   *jets_AK5PF_btag_secVertexHighEff;
   std::vector<float>   *jets_AK5PF_btag_secVertexCombined;
   std::vector<float>   *jets_AK5PF_jetCharge;
   std::vector<float>   *jets_AK5PF_chgEmE;
   std::vector<float>   *jets_AK5PF_chgHadE;
   std::vector<float>   *jets_AK5PF_chgMuE;
   std::vector<float>   *jets_AK5PF_chg_Mult;
   std::vector<float>   *jets_AK5PF_neutralEmE;
   std::vector<float>   *jets_AK5PF_neutralHadE;
   std::vector<float>   *jets_AK5PF_neutral_Mult;
   std::vector<float>   *jets_AK5PF_mu_Mult;
   std::vector<float>   *jets_AK5PF_emf;
   std::vector<float>   *jets_AK5PF_ehf;
   std::vector<float>   *jets_AK5PF_n60;
   std::vector<float>   *jets_AK5PF_n90;
   std::vector<float>   *jets_AK5PF_etaetaMoment;
   std::vector<float>   *jets_AK5PF_etaphiMoment;
   std::vector<float>   *jets_AK5PF_phiphiMoment;
   std::vector<float>   *jets_AK5PF_n90Hits;
   std::vector<float>   *jets_AK5PF_fHPD;
   std::vector<float>   *jets_AK5PF_fRBX;
   std::vector<float>   *jets_AK5PF_hitsInN90;
   std::vector<float>   *jets_AK5PF_nECALTowers;
   std::vector<float>   *jets_AK5PF_nHCALTowers;
   std::vector<float>   *jets_AK5PF_fSubDetector1;
   std::vector<float>   *jets_AK5PF_fSubDetector2;
   std::vector<float>   *jets_AK5PF_fSubDetector3;
   std::vector<float>   *jets_AK5PF_fSubDetector4;
   std::vector<float>   *jets_AK5PF_area;
   std::vector<float>   *jets_AK5PF_corrFactorRaw;
   std::vector<float>   *jets_AK5PF_mass;
   UInt_t          Njets_AK5PFclean;
   std::vector<float>   *jets_AK5PFclean_energy;
   std::vector<float>   *jets_AK5PFclean_et;
   std::vector<float>   *jets_AK5PFclean_eta;
   std::vector<float>   *jets_AK5PFclean_phi;
   std::vector<float>   *jets_AK5PFclean_pt;
   std::vector<float>   *jets_AK5PFclean_px;
   std::vector<float>   *jets_AK5PFclean_py;
   std::vector<float>   *jets_AK5PFclean_pz;
   std::vector<float>   *jets_AK5PFclean_status;
   std::vector<float>   *jets_AK5PFclean_theta;
   std::vector<float>   *jets_AK5PFclean_parton_Id;
   std::vector<float>   *jets_AK5PFclean_parton_motherId;
   std::vector<float>   *jets_AK5PFclean_parton_pt;
   std::vector<float>   *jets_AK5PFclean_parton_phi;
   std::vector<float>   *jets_AK5PFclean_parton_eta;
   std::vector<float>   *jets_AK5PFclean_parton_Energy;
   std::vector<float>   *jets_AK5PFclean_parton_mass;
   std::vector<float>   *jets_AK5PFclean_gen_et;
   std::vector<float>   *jets_AK5PFclean_gen_pt;
   std::vector<float>   *jets_AK5PFclean_gen_eta;
   std::vector<float>   *jets_AK5PFclean_gen_phi;
   std::vector<float>   *jets_AK5PFclean_gen_mass;
   std::vector<float>   *jets_AK5PFclean_gen_Energy;
   std::vector<float>   *jets_AK5PFclean_gen_Id;
   std::vector<float>   *jets_AK5PFclean_partonFlavour;
   std::vector<float>   *jets_AK5PFclean_btag_TC_highPur;
   std::vector<float>   *jets_AK5PFclean_btag_TC_highEff;
   std::vector<float>   *jets_AK5PFclean_btag_jetProb;
   std::vector<float>   *jets_AK5PFclean_btag_jetBProb;
   std::vector<float>   *jets_AK5PFclean_btag_softEle;
   std::vector<float>   *jets_AK5PFclean_btag_softMuon;
   std::vector<float>   *jets_AK5PFclean_btag_secVertexHighPur;
   std::vector<float>   *jets_AK5PFclean_btag_secVertexHighEff;
   std::vector<float>   *jets_AK5PFclean_btag_secVertexCombined;
   std::vector<float>   *jets_AK5PFclean_jetCharge;
   std::vector<float>   *jets_AK5PFclean_chgEmE;
   std::vector<float>   *jets_AK5PFclean_chgHadE;
   std::vector<float>   *jets_AK5PFclean_chgMuE;
   std::vector<float>   *jets_AK5PFclean_chg_Mult;
   std::vector<float>   *jets_AK5PFclean_neutralEmE;
   std::vector<float>   *jets_AK5PFclean_neutralHadE;
   std::vector<float>   *jets_AK5PFclean_neutral_Mult;
   std::vector<float>   *jets_AK5PFclean_mu_Mult;
   std::vector<float>   *jets_AK5PFclean_emf;
   std::vector<float>   *jets_AK5PFclean_ehf;
   std::vector<float>   *jets_AK5PFclean_n60;
   std::vector<float>   *jets_AK5PFclean_n90;
   std::vector<float>   *jets_AK5PFclean_etaetaMoment;
   std::vector<float>   *jets_AK5PFclean_etaphiMoment;
   std::vector<float>   *jets_AK5PFclean_phiphiMoment;
   std::vector<float>   *jets_AK5PFclean_n90Hits;
   std::vector<float>   *jets_AK5PFclean_fHPD;
   std::vector<float>   *jets_AK5PFclean_fRBX;
   std::vector<float>   *jets_AK5PFclean_hitsInN90;
   std::vector<float>   *jets_AK5PFclean_nECALTowers;
   std::vector<float>   *jets_AK5PFclean_nHCALTowers;
   std::vector<float>   *jets_AK5PFclean_fSubDetector1;
   std::vector<float>   *jets_AK5PFclean_fSubDetector2;
   std::vector<float>   *jets_AK5PFclean_fSubDetector3;
   std::vector<float>   *jets_AK5PFclean_fSubDetector4;
   std::vector<float>   *jets_AK5PFclean_area;
   std::vector<float>   *jets_AK5PFclean_corrFactorRaw;
   std::vector<float>   *jets_AK5PFclean_mass;
   UInt_t          Nmc_doc;
   std::vector<float>   *mc_doc_id;
   std::vector<float>   *mc_doc_pt;
   std::vector<float>   *mc_doc_px;
   std::vector<float>   *mc_doc_py;
   std::vector<float>   *mc_doc_pz;
   std::vector<float>   *mc_doc_eta;
   std::vector<float>   *mc_doc_phi;
   std::vector<float>   *mc_doc_theta;
   std::vector<float>   *mc_doc_energy;
   std::vector<float>   *mc_doc_status;
   std::vector<float>   *mc_doc_charge;
   std::vector<float>   *mc_doc_mother_id;
   std::vector<float>   *mc_doc_grandmother_id;
   std::vector<float>   *mc_doc_ggrandmother_id;
   std::vector<float>   *mc_doc_mother_pt;
   std::vector<float>   *mc_doc_vertex_x;
   std::vector<float>   *mc_doc_vertex_y;
   std::vector<float>   *mc_doc_vertex_z;
   std::vector<float>   *mc_doc_mass;
   std::vector<float>   *mc_doc_numOfDaughters;
   std::vector<float>   *mc_doc_numOfMothers;
   UInt_t          Nmc_electrons;
   std::vector<float>   *mc_electrons_id;
   std::vector<float>   *mc_electrons_pt;
   std::vector<float>   *mc_electrons_px;
   std::vector<float>   *mc_electrons_py;
   std::vector<float>   *mc_electrons_pz;
   std::vector<float>   *mc_electrons_eta;
   std::vector<float>   *mc_electrons_phi;
   std::vector<float>   *mc_electrons_theta;
   std::vector<float>   *mc_electrons_status;
   std::vector<float>   *mc_electrons_energy;
   std::vector<float>   *mc_electrons_charge;
   std::vector<float>   *mc_electrons_mother_id;
   std::vector<float>   *mc_electrons_mother_pt;
   std::vector<float>   *mc_electrons_grandmother_id;
   std::vector<float>   *mc_electrons_ggrandmother_id;
   std::vector<float>   *mc_electrons_vertex_x;
   std::vector<float>   *mc_electrons_vertex_y;
   std::vector<float>   *mc_electrons_vertex_z;
   std::vector<float>   *mc_electrons_mass;
   std::vector<float>   *mc_electrons_numOfDaughters;
   UInt_t          Nmc_mus;
   std::vector<float>   *mc_mus_id;
   std::vector<float>   *mc_mus_pt;
   std::vector<float>   *mc_mus_px;
   std::vector<float>   *mc_mus_py;
   std::vector<float>   *mc_mus_pz;
   std::vector<float>   *mc_mus_eta;
   std::vector<float>   *mc_mus_phi;
   std::vector<float>   *mc_mus_theta;
   std::vector<float>   *mc_mus_status;
   std::vector<float>   *mc_mus_energy;
   std::vector<float>   *mc_mus_charge;
   std::vector<float>   *mc_mus_mother_id;
   std::vector<float>   *mc_mus_mother_pt;
   std::vector<float>   *mc_mus_grandmother_id;
   std::vector<float>   *mc_mus_ggrandmother_id;
   std::vector<float>   *mc_mus_vertex_x;
   std::vector<float>   *mc_mus_vertex_y;
   std::vector<float>   *mc_mus_vertex_z;
   std::vector<float>   *mc_mus_mass;
   std::vector<float>   *mc_mus_numOfDaughters;
   UInt_t          Nmc_photons;
   std::vector<float>   *mc_photons_id;
   std::vector<float>   *mc_photons_pt;
   std::vector<float>   *mc_photons_px;
   std::vector<float>   *mc_photons_py;
   std::vector<float>   *mc_photons_pz;
   std::vector<float>   *mc_photons_eta;
   std::vector<float>   *mc_photons_phi;
   std::vector<float>   *mc_photons_theta;
   std::vector<float>   *mc_photons_status;
   std::vector<float>   *mc_photons_energy;
   std::vector<float>   *mc_photons_charge;
   std::vector<float>   *mc_photons_mother_id;
   std::vector<float>   *mc_photons_mother_pt;
   std::vector<float>   *mc_photons_grandmother_id;
   std::vector<float>   *mc_photons_ggrandmother_id;
   std::vector<float>   *mc_photons_vertex_x;
   std::vector<float>   *mc_photons_vertex_y;
   std::vector<float>   *mc_photons_vertex_z;
   std::vector<float>   *mc_photons_mass;
   std::vector<float>   *mc_photons_numOfDaughters;
   UInt_t          Nmc_taus;
   std::vector<float>   *mc_taus_id;
   std::vector<float>   *mc_taus_pt;
   std::vector<float>   *mc_taus_px;
   std::vector<float>   *mc_taus_py;
   std::vector<float>   *mc_taus_pz;
   std::vector<float>   *mc_taus_eta;
   std::vector<float>   *mc_taus_phi;
   std::vector<float>   *mc_taus_theta;
   std::vector<float>   *mc_taus_status;
   std::vector<float>   *mc_taus_energy;
   std::vector<float>   *mc_taus_charge;
   std::vector<float>   *mc_taus_mother_id;
   std::vector<float>   *mc_taus_mother_pt;
   std::vector<float>   *mc_taus_grandmother_id;
   std::vector<float>   *mc_taus_ggrandmother_id;
   std::vector<float>   *mc_taus_vertex_x;
   std::vector<float>   *mc_taus_vertex_y;
   std::vector<float>   *mc_taus_vertex_z;
   std::vector<float>   *mc_taus_mass;
   std::vector<float>   *mc_taus_numOfDaughters;
   UInt_t          Nmets_AK5;
   std::vector<float>   *mets_AK5_et;
   std::vector<float>   *mets_AK5_phi;
   std::vector<float>   *mets_AK5_ex;
   std::vector<float>   *mets_AK5_ey;
   std::vector<float>   *mets_AK5_gen_et;
   std::vector<float>   *mets_AK5_gen_phi;
   std::vector<float>   *mets_AK5_sign;
   std::vector<float>   *mets_AK5_sumEt;
   std::vector<float>   *mets_AK5_unCPhi;
   std::vector<float>   *mets_AK5_unCPt;
   UInt_t          Nmulti5x5EBC;
   std::vector<float>   *multi5x5EBC_energy;
   std::vector<float>   *multi5x5EBC_x;
   std::vector<float>   *multi5x5EBC_y;
   std::vector<float>   *multi5x5EBC_z;
   std::vector<float>   *multi5x5EBC_rho;
   std::vector<float>   *multi5x5EBC_phi;
   std::vector<float>   *multi5x5EBC_eta;
   std::vector<float>   *multi5x5EBC_theta;
   UInt_t          Nmus;
   std::vector<float>   *pf_mus_energy;
   std::vector<float>   *pf_mus_et;
   std::vector<float>   *pf_mus_eta;
   std::vector<float>   *pf_mus_phi;
   std::vector<float>   *pf_mus_pt;
   std::vector<float>   *pf_mus_px;
   std::vector<float>   *pf_mus_py;
   std::vector<float>   *pf_mus_pz;
   std::vector<float>   *pf_mus_status;
   std::vector<float>   *pf_mus_theta;
   std::vector<float>   *mus_eta;
   std::vector<float>   *mus_phi;
   std::vector<float>   *mus_pt;   
   std::vector<float>   *mus_px;   
   std::vector<float>   *mus_py;   
   std::vector<float>   *mus_pz; 
   std::vector<float>   *mus_charge;   
   std::vector<float>   *mus_energy;
   std::vector<float>   *mus_gen_id;
   std::vector<float>   *mus_gen_phi;
   std::vector<float>   *mus_gen_pt;
   std::vector<float>   *mus_gen_pz;
   std::vector<float>   *mus_gen_px;
   std::vector<float>   *mus_gen_py;
   std::vector<float>   *mus_gen_eta;
   std::vector<float>   *mus_gen_theta;
   std::vector<float>   *mus_gen_et;
   std::vector<float>   *mus_gen_mother_id;
   std::vector<float>   *mus_gen_mother_phi;
   std::vector<float>   *mus_gen_mother_pt;
   std::vector<float>   *mus_gen_mother_pz;
   std::vector<float>   *mus_gen_mother_px;
   std::vector<float>   *mus_gen_mother_py;
   std::vector<float>   *mus_gen_mother_eta;
   std::vector<float>   *mus_gen_mother_theta;
   std::vector<float>   *mus_gen_mother_et;
   std::vector<float>   *pf_mus_tkHits;
   std::vector<float>   *pf_mus_cIso;
   std::vector<float>   *pf_mus_tIso;
   std::vector<float>   *pf_mus_ecalIso;
   std::vector<float>   *pf_mus_hcalIso;
   std::vector<float>   *pf_mus_ecalvetoDep;
   std::vector<float>   *pf_mus_hcalvetoDep;
   std::vector<float>   *pf_mus_calEnergyEm;
   std::vector<float>   *pf_mus_calEnergyHad;
   std::vector<float>   *pf_mus_calEnergyHo;
   std::vector<float>   *pf_mus_calEnergyEmS9;
   std::vector<float>   *pf_mus_calEnergyHadS9;
   std::vector<float>   *pf_mus_calEnergyHoS9;
   std::vector<float>   *pf_mus_iso03_sumPt;
   std::vector<float>   *pf_mus_iso03_emEt;
   std::vector<float>   *pf_mus_iso03_hadEt;
   std::vector<float>   *pf_mus_iso03_hoEt;
   std::vector<float>   *pf_mus_iso03_nTracks;
   std::vector<float>   *pf_mus_iso05_sumPt;
   std::vector<float>   *pf_mus_iso05_emEt;
   std::vector<float>   *pf_mus_iso05_hadEt;
   std::vector<float>   *pf_mus_iso05_hoEt;
   std::vector<float>   *pf_mus_iso05_nTracks;
   std::vector<float>   *pf_mus_chargedHadronIso;
   std::vector<float>   *pf_mus_photonIso;
   std::vector<float>   *pf_mus_neutralHadronIso;
   std::vector<float>   *pf_mus_charge;
   std::vector<float>   *pf_mus_cm_chi2;
   std::vector<float>   *pf_mus_cm_ndof;
   std::vector<float>   *pf_mus_cm_chg;
   std::vector<float>   *pf_mus_cm_pt;
   std::vector<float>   *pf_mus_cm_px;
   std::vector<float>   *pf_mus_cm_py;
   std::vector<float>   *pf_mus_cm_pz;
   std::vector<float>   *pf_mus_cm_eta;
   std::vector<float>   *pf_mus_cm_phi;
   std::vector<float>   *pf_mus_cm_theta;
   std::vector<float>   *pf_mus_cm_d0dum;
   std::vector<float>   *pf_mus_cm_dz;
   std::vector<float>   *pf_mus_cm_vx;
   std::vector<float>   *pf_mus_cm_vy;
   std::vector<float>   *pf_mus_cm_vz;
   std::vector<float>   *pf_mus_cm_numvalhits;
   std::vector<float>   *pf_mus_cm_numlosthits;
   std::vector<float>   *pf_mus_cm_numvalMuonhits;
   std::vector<float>   *pf_mus_cm_d0dumErr;
   std::vector<float>   *pf_mus_cm_dzErr;
   std::vector<float>   *pf_mus_cm_ptErr;
   std::vector<float>   *pf_mus_cm_etaErr;
   std::vector<float>   *pf_mus_cm_phiErr;
   std::vector<float>   *pf_mus_tk_id;
   std::vector<float>   *pf_mus_tk_chi2;
   std::vector<float>   *pf_mus_tk_ndof;
   std::vector<float>   *pf_mus_tk_chg;
   std::vector<float>   *pf_mus_tk_pt;
   std::vector<float>   *pf_mus_tk_px;
   std::vector<float>   *pf_mus_tk_py;
   std::vector<float>   *pf_mus_tk_pz;
   std::vector<float>   *pf_mus_tk_eta;
   std::vector<float>   *pf_mus_tk_phi;
   std::vector<float>   *pf_mus_tk_theta;
   std::vector<float>   *pf_mus_tk_d0dum;
   std::vector<float>   *pf_mus_tk_dz;
   std::vector<float>   *pf_mus_tk_vx;
   std::vector<float>   *pf_mus_tk_vy;
   std::vector<float>   *pf_mus_tk_vz;
   std::vector<float>   *pf_mus_tk_numvalhits;
   std::vector<float>   *pf_mus_tk_numlosthits;
   std::vector<float>   *pf_mus_tk_d0dumErr;
   std::vector<float>   *pf_mus_tk_dzErr;
   std::vector<float>   *pf_mus_tk_ptErr;
   std::vector<float>   *pf_mus_tk_etaErr;
   std::vector<float>   *pf_mus_tk_phiErr;
   std::vector<float>   *pf_mus_tk_numvalPixelhits;
   std::vector<float>   *pf_mus_tk_numpixelWthMeasr;
   std::vector<float>   *pf_mus_stamu_chi2;
   std::vector<float>   *pf_mus_stamu_ndof;
   std::vector<float>   *pf_mus_stamu_chg;
   std::vector<float>   *pf_mus_stamu_pt;
   std::vector<float>   *pf_mus_stamu_px;
   std::vector<float>   *pf_mus_stamu_py;
   std::vector<float>   *pf_mus_stamu_pz;
   std::vector<float>   *pf_mus_stamu_eta;
   std::vector<float>   *pf_mus_stamu_phi;
   std::vector<float>   *pf_mus_stamu_theta;
   std::vector<float>   *pf_mus_stamu_d0dum;
   std::vector<float>   *pf_mus_stamu_dz;
   std::vector<float>   *pf_mus_stamu_vx;
   std::vector<float>   *pf_mus_stamu_vy;
   std::vector<float>   *pf_mus_stamu_vz;
   std::vector<float>   *pf_mus_stamu_numvalhits;
   std::vector<float>   *pf_mus_stamu_numlosthits;
   std::vector<float>   *pf_mus_stamu_d0dumErr;
   std::vector<float>   *pf_mus_stamu_dzErr;
   std::vector<float>   *pf_mus_stamu_ptErr;
   std::vector<float>   *pf_mus_stamu_etaErr;
   std::vector<float>   *pf_mus_stamu_phiErr;
   std::vector<float>   *pf_mus_num_matches;
   std::vector<float>   *pf_mus_isTrackerMuon;
   std::vector<float>   *pf_mus_isStandAloneMuon;
   std::vector<float>   *pf_mus_isCaloMuon;
   std::vector<float>   *pf_mus_isGlobalMuon;
   std::vector<float>   *pf_mus_isElectron;
   std::vector<float>   *pf_mus_isConvertedPhoton;
   std::vector<float>   *pf_mus_isPhoton;
   std::vector<float>   *pf_mus_id_All;
   std::vector<float>   *pf_mus_id_AllGlobalMuons;
   std::vector<float>   *pf_mus_id_AllStandAloneMuons;
   std::vector<float>   *pf_mus_id_AllTrackerMuons;
   std::vector<float>   *pf_mus_id_TrackerMuonArbitrated;
   std::vector<float>   *pf_mus_id_AllArbitrated;
   std::vector<float>   *pf_mus_id_GlobalMuonPromptTight;
   std::vector<float>   *pf_mus_id_TMLastStationLoose;
   std::vector<float>   *pf_mus_id_TMLastStationTight;
   std::vector<float>   *pf_mus_id_TM2DCompatibilityLoose;
   std::vector<float>   *pf_mus_id_TM2DCompatibilityTight;
   std::vector<float>   *pf_mus_id_TMOneStationLoose;
   std::vector<float>   *pf_mus_id_TMOneStationTight;
   std::vector<float>   *pf_mus_id_TMLastStationOptimizedLowPtLoose;
   std::vector<float>   *pf_mus_id_TMLastStationOptimizedLowPtTight;
   UInt_t          NpfTypeImets;
   std::vector<float>   *pfTypeImets_et;
   std::vector<float>   *pfTypeImets_phi;
   std::vector<float>   *pfTypeImets_ex;
   std::vector<float>   *pfTypeImets_ey;
   std::vector<float>   *pfTypeImets_gen_et;
   std::vector<float>   *pfTypeImets_gen_phi;
   std::vector<float>   *pfTypeImets_sign;
   std::vector<float>   *pfTypeImets_sumEt;
   std::vector<float>   *pfTypeImets_unCPhi;
   std::vector<float>   *pfTypeImets_unCPt;
   UInt_t          Npf_els;

   UInt_t          Npf_mus;
   std::vector<float>   *pf_mus_gen_id;
   std::vector<float>   *pf_mus_gen_phi;
   std::vector<float>   *pf_mus_gen_pt;
   std::vector<float>   *pf_mus_gen_pz;
   std::vector<float>   *pf_mus_gen_px;
   std::vector<float>   *pf_mus_gen_py;
   std::vector<float>   *pf_mus_gen_eta;
   std::vector<float>   *pf_mus_gen_theta;
   std::vector<float>   *pf_mus_gen_et;
   std::vector<float>   *pf_mus_gen_mother_id;
   std::vector<float>   *pf_mus_gen_mother_phi;
   std::vector<float>   *pf_mus_gen_mother_pt;
   std::vector<float>   *pf_mus_gen_mother_pz;
   std::vector<float>   *pf_mus_gen_mother_px;
   std::vector<float>   *pf_mus_gen_mother_py;
   std::vector<float>   *pf_mus_gen_mother_eta;
   std::vector<float>   *pf_mus_gen_mother_theta;
   std::vector<float>   *pf_mus_gen_mother_et;
   UInt_t          Npfcand_els;
   UInt_t          Npfcand_mus;
   std::vector<float>   *pf_mus_particleId;
   UInt_t          Npfmets;
   std::vector<float>   *pfmets_et;
   std::vector<float>   *pfmets_phi;
   std::vector<float>   *pfmets_ex;
   std::vector<float>   *pfmets_ey;
   std::vector<float>   *pfmets_gen_et;
   std::vector<float>   *pfmets_gen_phi;
   std::vector<float>   *pfmets_sign;
   std::vector<float>   *pfmets_sumEt;
   std::vector<float>   *pfmets_unCPhi;
   std::vector<float>   *pfmets_unCPt;
   UInt_t          Nphotons;
   std::vector<float>   *photons_energy;
   std::vector<float>   *photons_et;
   std::vector<float>   *photons_eta;
   std::vector<float>   *photons_phi;
   std::vector<float>   *photons_pt;
   std::vector<float>   *photons_px;
   std::vector<float>   *photons_py;
   std::vector<float>   *photons_pz;
   std::vector<float>   *photons_status;
   std::vector<float>   *photons_theta;
   std::vector<float>   *photons_hadOverEM;
   std::vector<float>   *photons_scEnergy;
   std::vector<float>   *photons_scRawEnergy;
   std::vector<float>   *photons_scEta;
   std::vector<float>   *photons_scPhi;
   std::vector<float>   *photons_scEtaWidth;
   std::vector<float>   *photons_scPhiWidth;
   std::vector<float>   *photons_tIso;
   std::vector<float>   *photons_ecalIso;
   std::vector<float>   *photons_hcalIso;
   std::vector<float>   *photons_isoEcalRecHitDR04;
   std::vector<float>   *photons_isoHcalRecHitDR04;
   std::vector<float>   *photons_isoSolidTrkConeDR04;
   std::vector<float>   *photons_isoHollowTrkConeDR04;
   std::vector<float>   *photons_nTrkSolidConeDR04;
   std::vector<float>   *photons_nTrkHollowConeDR04;
   std::vector<float>   *photons_isoEcalRecHitDR03;
   std::vector<float>   *photons_isoHcalRecHitDR03;
   std::vector<float>   *photons_isoSolidTrkConeDR03;
   std::vector<float>   *photons_isoHollowTrkConeDR03;
   std::vector<float>   *photons_nTrkSolidConeDR03;
   std::vector<float>   *photons_nTrkHollowConeDR03;
   std::vector<float>   *photons_isAlsoElectron;
   std::vector<float>   *photons_hasPixelSeed;
   std::vector<float>   *photons_isConverted;
   std::vector<float>   *photons_isEBGap;
   std::vector<float>   *photons_isEEGap;
   std::vector<float>   *photons_isEBEEGap;
   std::vector<float>   *photons_isEBPho;
   std::vector<float>   *photons_isEEPho;
   std::vector<float>   *photons_isLoosePhoton;
   std::vector<float>   *photons_isTightPhoton;
   std::vector<float>   *photons_maxEnergyXtal;
   std::vector<float>   *photons_e1x5;
   std::vector<float>   *photons_e2x5;
   std::vector<float>   *photons_e3x3;
   std::vector<float>   *photons_e5x5;
   std::vector<float>   *photons_sigmaEtaEta;
   std::vector<float>   *photons_sigmaIetaIeta;
   std::vector<float>   *photons_r9;
   std::vector<float>   *photons_gen_et;
   std::vector<float>   *photons_gen_eta;
   std::vector<float>   *photons_gen_phi;
   std::vector<float>   *photons_gen_id;
   UInt_t          Npv;
   std::vector<float>   *pv_x;
   std::vector<float>   *pv_y;
   std::vector<float>   *pv_z;
   std::vector<float>   *pv_xErr;
   std::vector<float>   *pv_yErr;
   std::vector<float>   *pv_zErr;
   std::vector<float>   *pv_chi2;
   std::vector<float>   *pv_ndof;
   std::vector<float>   *pv_isFake;
   std::vector<float>   *pv_tracksSize;
   UInt_t          Ntaus;
   std::vector<float>   *taus_energy;
   std::vector<float>   *taus_et;
   std::vector<float>   *taus_eta;
   std::vector<float>   *taus_phi;
   std::vector<float>   *taus_pt;
   std::vector<float>   *taus_px;
   std::vector<float>   *taus_py;
   std::vector<float>   *taus_pz;
   std::vector<float>   *taus_status;
   std::vector<float>   *taus_theta;
   std::vector<float>   *taus_charge;
   std::vector<float>   *taus_emf;
   std::vector<float>   *taus_hcalTotOverPLead;
   std::vector<float>   *taus_hcalMaxOverPLead;
   std::vector<float>   *taus_hcal3x3OverPLead;
   std::vector<float>   *taus_ecalStripSumEOverPLead;
   std::vector<float>   *taus_elecPreIdOutput;
   std::vector<float>   *taus_elecPreIdDecision;
   std::vector<float>   *taus_leadPFChargedHadrCand_pt;
   std::vector<float>   *taus_leadPFChargedHadrCand_charge;
   std::vector<float>   *taus_leadPFChargedHadrCand_eta;
   std::vector<float>   *taus_leadPFChargedHadrCand_ECAL_eta;
   std::vector<float>   *taus_leadPFChargedHadrCand_phi;
   std::vector<float>   *taus_isoPFGammaCandsEtSum;
   std::vector<float>   *taus_isoPFChargedHadrCandsPtSum;
   std::vector<float>   *taus_leadingTrackFinding;
   std::vector<float>   *taus_leadingTrackPtCut;
   std::vector<float>   *taus_trackIsolation;
   std::vector<float>   *taus_ecalIsolation;
   std::vector<float>   *taus_byIsolation;
   std::vector<float>   *taus_againstElectron;
   std::vector<float>   *taus_againstMuon;
   std::vector<float>   *taus_taNC_quarter;
   std::vector<float>   *taus_taNC_one;
   std::vector<float>   *taus_taNC_half;
   std::vector<float>   *taus_taNC_tenth;
   std::vector<float>   *taus_taNC;
   std::vector<float>   *taus_byIsoUsingLeadingPi;
   std::vector<float>   *taus_tkIsoUsingLeadingPi;
   std::vector<float>   *taus_ecalIsoUsingLeadingPi;
   std::vector<float>   *taus_signalPFChargedHadrCandsSize;
   std::vector<float>   *taus_muDecision;
   std::vector<float>   *taus_Nprongs;
   UInt_t          Ntcmets;
   std::vector<float>   *tcmets_et;
   std::vector<float>   *tcmets_phi;
   std::vector<float>   *tcmets_ex;
   std::vector<float>   *tcmets_ey;
   std::vector<float>   *tcmets_sumEt;
   UInt_t          Ntracks;
   std::vector<float>   *tracks_chi2;
   std::vector<float>   *tracks_trkExptHitsInner;
   std::vector<float>   *tracks_trkExptHitsOuter;
   std::vector<float>   *tracks_trks_nlayerslost;
   std::vector<float>   *tracks_trks_nlayers;
   std::vector<float>   *tracks_trksvalidpixelhits;
   std::vector<float>   *tracks_trkslostpixelhits;
   std::vector<float>   *tracks_ndof;
   std::vector<float>   *tracks_chg;
   std::vector<float>   *tracks_pt;
   std::vector<float>   *tracks_px;
   std::vector<float>   *tracks_py;
   std::vector<float>   *tracks_pz;
   std::vector<float>   *tracks_eta;
   std::vector<float>   *tracks_phi;
   std::vector<float>   *tracks_theta;
   std::vector<float>   *tracks_d0dum;
   std::vector<float>   *tracks_dz;
   std::vector<float>   *tracks_vx;
   std::vector<float>   *tracks_vy;
   std::vector<float>   *tracks_vz;
   std::vector<float>   *tracks_numvalhits;
   std::vector<float>   *tracks_numlosthits;
   std::vector<float>   *tracks_d0dumErr;
   std::vector<float>   *tracks_dzErr;
   std::vector<float>   *tracks_ptErr;
   std::vector<float>   *tracks_etaErr;
   std::vector<float>   *tracks_phiErr;
   std::vector<float>   *tracks_Nrechits;
   std::vector<float>   *tracks_innerHitX;
   std::vector<float>   *tracks_innerHitY;
   std::vector<float>   *tracks_innerHitZ;
   std::vector<float>   *tracks_outerHitX;
   std::vector<float>   *tracks_outerHitY;
   std::vector<float>   *tracks_outerHitZ;
   std::vector<float>   *tracks_outerPx;
   std::vector<float>   *tracks_outerPy;
   std::vector<float>   *tracks_outerPz;
   std::vector<float>   *tracks_algo;
   std::vector<float>   *tracks_highPurity;
   UInt_t          run;
   UInt_t          event;
   UInt_t          lumiblock;
   UInt_t          experimentType;
   UInt_t          bunchCrossing;
   UInt_t          orbitNumber;
   
   // List of branches

   TBranch        *b_trigger_prescalevalue;   //!
   TBranch        *b_trigger_name;   //!
   TBranch        *b_trigger_decision;   //!
   TBranch        *b_trigger_lastfiltername;   //!
   TBranch        *b_triggerobject_pt;   //!
   TBranch        *b_triggerobject_px;   //!
   TBranch        *b_triggerobject_py;   //!
   TBranch        *b_triggerobject_pz;   //!
   TBranch        *b_triggerobject_et;   //!
   TBranch        *b_triggerobject_energy;   //!
   TBranch        *b_triggerobject_phi;   //!
   TBranch        *b_triggerobject_eta;   //!
   TBranch        *b_triggerobject_collectionname;   //!
   TBranch        *b_L1trigger_bit;   //!
   TBranch        *b_L1trigger_techTrigger;   //!
   TBranch        *b_L1trigger_prescalevalue;   //!
   TBranch        *b_L1trigger_name;   //!
   TBranch        *b_L1trigger_alias;   //!
   TBranch        *b_L1trigger_decision;   //!
   TBranch        *b_L1trigger_decision_nomask;   //!
   TBranch        *b_els_conversion_dist;   //!
   TBranch        *b_els_conversion_dcot;   //!
   TBranch        *b_hbhefilter_decision;   //!
   TBranch        *b_trackingfailurefilter_decision;   //!
   TBranch	  *b_cschalofilter_decision;   //!
   TBranch	  *b_ecalTPfilter_decision;   //!
   TBranch	  *b_scrapingVeto_decision;   //!
   TBranch	  *b_inconsistentPFmuonfilter_decision;   //!
   TBranch	  *b_greedymuonfilter_decision;   //!
   TBranch	  *b_eenoisefilter_decision;   //!
   TBranch        *b_NbeamSpot;   //!
   TBranch        *b_beamSpot_x;   //!
   TBranch        *b_beamSpot_y;   //!
   TBranch        *b_beamSpot_z;   //!
   TBranch        *b_beamSpot_x0Error;   //!
   TBranch        *b_beamSpot_y0Error;   //!
   TBranch        *b_beamSpot_z0Error;   //!
   TBranch        *b_beamSpot_sigmaZ;   //!
   TBranch        *b_beamSpot_sigmaZ0Error;   //!
   TBranch        *b_beamSpot_dxdz;   //!
   TBranch        *b_beamSpot_dxdzError;   //!
   TBranch        *b_beamSpot_dydz;   //!
   TBranch        *b_beamSpot_dydzError;   //!
   TBranch        *b_beamSpot_beamWidthX;   //!
   TBranch        *b_beamSpot_beamWidthY;   //!
   TBranch        *b_beamSpot_beamWidthXError;   //!
   TBranch        *b_beamSpot_beamWidthYError;   //!
   TBranch        *b_Nels;   //!
   TBranch        *b_pf_els_energy;   //!
   TBranch        *b_pf_els_et;   //!
   TBranch        *b_pf_els_eta;   //!
   TBranch        *b_pf_els_phi;   //!
   TBranch        *b_pf_els_pt;   //!
   TBranch        *b_pf_els_px;   //!
   TBranch        *b_pf_els_py;   //!
   TBranch        *b_pf_els_pz;   //!
   TBranch        *b_pf_els_status;   //!
   TBranch        *b_pf_els_theta;   //!
   TBranch        *b_els_pt;   //!
   TBranch        *b_els_eta;   //!
   TBranch        *b_els_phi;   //!
   TBranch        *b_els_charge;   //!
   TBranch        *b_els_gen_id;   //!
   TBranch        *b_els_gen_phi;   //!
   TBranch        *b_els_gen_pt;   //!
   TBranch        *b_els_gen_pz;   //!
   TBranch        *b_els_gen_px;   //!
   TBranch        *b_els_gen_py;   //!
   TBranch        *b_els_gen_eta;   //!
   TBranch        *b_els_gen_theta;   //!
   TBranch        *b_els_gen_et;   //!
   TBranch        *b_els_gen_mother_id;   //!
   TBranch        *b_els_gen_mother_phi;   //!
   TBranch        *b_els_gen_mother_pt;   //!
   TBranch        *b_els_gen_mother_pz;   //!
   TBranch        *b_els_gen_mother_px;   //!
   TBranch        *b_els_gen_mother_py;   //!
   TBranch        *b_els_gen_mother_eta;   //!
   TBranch        *b_els_gen_mother_theta;   //!
   TBranch        *b_els_gen_mother_et;   //!
   TBranch        *b_pf_els_tightId;   //!
   TBranch        *b_pf_els_looseId;   //!
   TBranch        *b_pf_els_robustTightId;   //!
   TBranch        *b_pf_els_robustLooseId;   //!
   TBranch        *b_pf_els_robustHighEnergyId;   //!
   TBranch        *b_pf_els_simpleEleId95relIso;   //!
   TBranch        *b_pf_els_simpleEleId90relIso;   //!
   TBranch        *b_pf_els_simpleEleId85relIso;   //!
   TBranch        *b_pf_els_simpleEleId80relIso;   //!
   TBranch        *b_pf_els_simpleEleId70relIso;   //!
   TBranch        *b_pf_els_simpleEleId60relIso;   //!
   TBranch        *b_pf_els_simpleEleId95cIso;   //!
   TBranch        *b_pf_els_simpleEleId90cIso;   //!
   TBranch        *b_pf_els_simpleEleId85cIso;   //!
   TBranch        *b_pf_els_simpleEleId80cIso;   //!
   TBranch        *b_pf_els_simpleEleId70cIso;   //!
   TBranch        *b_pf_els_simpleEleId60cIso;   //!
   TBranch        *b_pf_els_cIso;   //!
   TBranch        *b_pf_els_tIso;   //!
   TBranch        *b_pf_els_ecalIso;   //!
   TBranch        *b_pf_els_hcalIso;   //!
   TBranch        *b_pf_els_chi2;   //!
   TBranch        *b_pf_els_charge;   //!
   TBranch        *b_pf_els_caloEnergy;   //!
   TBranch        *b_pf_els_hadOverEm;   //!
   TBranch        *b_pf_els_eOverPIn;   //!
   TBranch        *b_pf_els_eSeedOverPOut;   //!
   TBranch        *b_pf_els_sigmaEtaEta;   //!
   TBranch        *b_pf_els_sigmaIEtaIEta;   //!
   TBranch        *b_pf_els_scEnergy;   //!
   TBranch        *b_pf_els_scRawEnergy;   //!
   TBranch        *b_pf_els_scSeedEnergy;   //!
   TBranch        *b_pf_els_scEta;   //!
   TBranch        *b_pf_els_scPhi;   //!
   TBranch        *b_pf_els_scEtaWidth;   //!
   TBranch        *b_pf_els_scPhiWidth;   //!
   TBranch        *b_pf_els_scE1x5;   //!
   TBranch        *b_pf_els_scE2x5Max;   //!
   TBranch        *b_pf_els_scE5x5;   //!
   TBranch        *b_pf_els_isEB;   //!
   TBranch        *b_pf_els_isEE;   //!
   TBranch        *b_pf_els_dEtaIn;   //!
   TBranch        *b_pf_els_dPhiIn;   //!
   TBranch        *b_pf_els_dEtaOut;   //!
   TBranch        *b_pf_els_dPhiOut;   //!
   TBranch        *b_pf_els_numvalhits;   //!
   TBranch        *b_pf_els_numlosthits;   //!
   TBranch        *b_pf_els_basicClustersSize;   //!
   TBranch        *b_pf_els_tk_pz;   //!
   TBranch        *b_pf_els_tk_pt;   //!
   TBranch        *b_pf_els_tk_phi;   //!
   TBranch        *b_pf_els_tk_eta;   //!
   TBranch        *b_pf_els_d0dum;   //!
   TBranch        *b_pf_els_dz;   //!
   TBranch        *b_pf_els_vx;   //!
   TBranch        *b_pf_els_vy;   //!
   TBranch        *b_pf_els_vz;   //!
   TBranch        *b_pf_els_ndof;   //!
   TBranch        *b_pf_els_ptError;   //!
   TBranch        *b_pf_els_d0dumError;   //!
   TBranch        *b_pf_els_dzError;   //!
   TBranch        *b_pf_els_etaError;   //!
   TBranch        *b_pf_els_phiError;   //!
   TBranch        *b_pf_els_tk_charge;   //!
   TBranch        *b_pf_els_core_ecalDrivenSeed;   //!
   TBranch        *b_pf_els_n_inner_layer;   //!
   TBranch        *b_pf_els_n_outer_layer;   //!
   TBranch        *b_pf_els_ctf_tk_id;   //!
   TBranch        *b_pf_els_ctf_tk_charge;   //!
   TBranch        *b_pf_els_ctf_tk_eta;   //!
   TBranch        *b_pf_els_ctf_tk_phi;   //!
   TBranch        *b_pf_els_fbrem;   //!
   TBranch        *b_pf_els_shFracInnerHits;   //!
   TBranch        *b_pf_els_dr03EcalRecHitSumEt;   //!
   TBranch        *b_pf_els_dr03HcalTowerSumEt;   //!
   TBranch        *b_pf_els_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_pf_els_dr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_pf_els_dr03TkSumPt;   //!
   TBranch        *b_pf_els_dr04EcalRecHitSumEt;   //!
   TBranch        *b_pf_els_dr04HcalTowerSumEt;   //!
   TBranch        *b_pf_els_dr04HcalDepth1TowerSumEt;   //!
   TBranch        *b_pf_els_dr04HcalDepth2TowerSumEt;   //!
   TBranch        *b_pf_els_dr04TkSumPt;   //!
   TBranch        *b_pf_els_chargedHadronIso;   //!
   TBranch        *b_pf_els_photonIso;   //!
   TBranch        *b_pf_els_neutralHadronIso;   //!
   TBranch        *b_pf_els_cpx;   //!
   TBranch        *b_pf_els_cpy;   //!
   TBranch        *b_pf_els_cpz;   //!
   TBranch        *b_pf_els_vpx;   //!
   TBranch        *b_pf_els_vpy;   //!
   TBranch        *b_pf_els_vpz;   //!
   TBranch        *b_pf_els_cx;   //!
   TBranch        *b_pf_els_cy;   //!
   TBranch        *b_pf_els_cz;   //!
   TBranch        *b_NhcalNoiseSummary;   //!
   TBranch        *b_hcalNoiseSummary_passLooseNoiseFilter;   //!
   TBranch        *b_hcalNoiseSummary_passTightNoiseFilter;   //!
   TBranch        *b_hcalNoiseSummary_passHighLevelNoiseFilter;   //!
   TBranch        *b_hcalNoiseSummary_noiseFilterStatus;   //!
   TBranch        *b_hcalNoiseSummary_noiseType;   //!
   TBranch        *b_hcalNoiseSummary_eventEMEnergy;   //!
   TBranch        *b_hcalNoiseSummary_eventHadEnergy;   //!
   TBranch        *b_hcalNoiseSummary_eventTrackEnergy;   //!
   TBranch        *b_hcalNoiseSummary_eventEMFraction;   //!
   TBranch        *b_hcalNoiseSummary_eventChargeFraction;   //!
   TBranch        *b_hcalNoiseSummary_min10GeVHitTime;   //!
   TBranch        *b_hcalNoiseSummary_max10GeVHitTime;   //!
   TBranch        *b_hcalNoiseSummary_rms10GeVHitTime;   //!
   TBranch        *b_hcalNoiseSummary_min25GeVHitTime;   //!
   TBranch        *b_hcalNoiseSummary_max25GeVHitTime;   //!
   TBranch        *b_hcalNoiseSummary_rms25GeVHitTime;   //!
   TBranch        *b_hcalNoiseSummary_num10GeVHits;   //!
   TBranch        *b_hcalNoiseSummary_num25GeVHits;   //!
   TBranch        *b_hcalNoiseSummary_minE2TS;   //!
   TBranch        *b_hcalNoiseSummary_minE10TS;   //!
   TBranch        *b_hcalNoiseSummary_minE2Over10TS;   //!
   TBranch        *b_hcalNoiseSummary_maxZeros;   //!
   TBranch        *b_hcalNoiseSummary_maxHPDHits;   //!
   TBranch        *b_hcalNoiseSummary_maxRBXHits;   //!
   TBranch        *b_hcalNoiseSummary_minHPDEMF;   //!
   TBranch        *b_hcalNoiseSummary_minRBXEMF;   //!
   TBranch        *b_hcalNoiseSummary_numProblematicRBXs;   //!
   TBranch        *b_hcalNoiseSummary_maxE2Over10TS;   //!
   TBranch        *b_hcalNoiseSummary_maxHPDNoOtherHits;   //!
   TBranch        *b_NhybridBBC;   //!
   TBranch        *b_hybridBBC_energy;   //!
   TBranch        *b_hybridBBC_x;   //!
   TBranch        *b_hybridBBC_y;   //!
   TBranch        *b_hybridBBC_z;   //!
   TBranch        *b_hybridBBC_rho;   //!
   TBranch        *b_hybridBBC_phi;   //!
   TBranch        *b_hybridBBC_eta;   //!
   TBranch        *b_hybridBBC_theta;   //!
   TBranch        *b_Njets_AK5;   //!
   TBranch        *b_jets_AK5_energy;   //!
   TBranch        *b_jets_AK5_et;   //!
   TBranch        *b_jets_AK5_eta;   //!
   TBranch        *b_jets_AK5_phi;   //!
   TBranch        *b_jets_AK5_pt;   //!
   TBranch        *b_jets_AK5_px;   //!
   TBranch        *b_jets_AK5_py;   //!
   TBranch        *b_jets_AK5_pz;   //!
   TBranch        *b_jets_AK5_status;   //!
   TBranch        *b_jets_AK5_theta;   //!
   TBranch        *b_jets_AK5_parton_Id;   //!
   TBranch        *b_jets_AK5_parton_motherId;   //!
   TBranch        *b_jets_AK5_parton_pt;   //!
   TBranch        *b_jets_AK5_parton_phi;   //!
   TBranch        *b_jets_AK5_parton_eta;   //!
   TBranch        *b_jets_AK5_parton_Energy;   //!
   TBranch        *b_jets_AK5_parton_mass;   //!
   TBranch        *b_jets_AK5_gen_et;   //!
   TBranch        *b_jets_AK5_gen_pt;   //!
   TBranch        *b_jets_AK5_gen_eta;   //!
   TBranch        *b_jets_AK5_gen_phi;   //!
   TBranch        *b_jets_AK5_gen_mass;   //!
   TBranch        *b_jets_AK5_gen_Energy;   //!
   TBranch        *b_jets_AK5_gen_Id;   //!
   TBranch        *b_jets_AK5_gen_motherID;   //!
   TBranch        *b_jets_AK5_gen_threeCharge;   //!
   TBranch        *b_jets_AK5_partonFlavour;   //!
   TBranch        *b_jets_AK5_btag_TC_highPur;   //!
   TBranch        *b_jets_AK5_btag_TC_highEff;   //!
   TBranch        *b_jets_AK5_btag_jetProb;   //!
   TBranch        *b_jets_AK5_btag_jetBProb;   //!
   TBranch        *b_jets_AK5_btag_softEle;   //!
   TBranch        *b_jets_AK5_btag_softMuon;   //!
   TBranch        *b_jets_AK5_btag_secVertexHighPur;   //!
   TBranch        *b_jets_AK5_btag_secVertexHighEff;   //!
   TBranch        *b_jets_AK5_btag_secVertexCombined;   //!
   TBranch        *b_jets_AK5_jetCharge;   //!
   TBranch        *b_jets_AK5_chgEmE;   //!
   TBranch        *b_jets_AK5_chgHadE;   //!
   TBranch        *b_jets_AK5_chgMuE;   //!
   TBranch        *b_jets_AK5_chg_Mult;   //!
   TBranch        *b_jets_AK5_neutralEmE;   //!
   TBranch        *b_jets_AK5_neutralHadE;   //!
   TBranch        *b_jets_AK5_neutral_Mult;   //!
   TBranch        *b_jets_AK5_mu_Mult;   //!
   TBranch        *b_jets_AK5_emf;   //!
   TBranch        *b_jets_AK5_ehf;   //!
   TBranch        *b_jets_AK5_n60;   //!
   TBranch        *b_jets_AK5_n90;   //!
   TBranch        *b_jets_AK5_etaetaMoment;   //!
   TBranch        *b_jets_AK5_etaphiMoment;   //!
   TBranch        *b_jets_AK5_phiphiMoment;   //!
   TBranch        *b_jets_AK5_n90Hits;   //!
   TBranch        *b_jets_AK5_fHPD;   //!
   TBranch        *b_jets_AK5_fRBX;   //!
   TBranch        *b_jets_AK5_hitsInN90;   //!
   TBranch        *b_jets_AK5_nECALTowers;   //!
   TBranch        *b_jets_AK5_nHCALTowers;   //!
   TBranch        *b_jets_AK5_fSubDetector1;   //!
   TBranch        *b_jets_AK5_fSubDetector2;   //!
   TBranch        *b_jets_AK5_fSubDetector3;   //!
   TBranch        *b_jets_AK5_fSubDetector4;   //!
   TBranch        *b_jets_AK5_area;   //!
   TBranch        *b_jets_AK5_corrFactorRaw;   //!
   TBranch        *b_jets_AK5_mass;   //!
   TBranch        *b_Njets_AK5JPT;   //!
   TBranch        *b_jets_AK5JPT_energy;   //!
   TBranch        *b_jets_AK5JPT_et;   //!
   TBranch        *b_jets_AK5JPT_eta;   //!
   TBranch        *b_jets_AK5JPT_phi;   //!
   TBranch        *b_jets_AK5JPT_pt;   //!
   TBranch        *b_jets_AK5JPT_px;   //!
   TBranch        *b_jets_AK5JPT_py;   //!
   TBranch        *b_jets_AK5JPT_pz;   //!
   TBranch        *b_jets_AK5JPT_status;   //!
   TBranch        *b_jets_AK5JPT_theta;   //!
   TBranch        *b_jets_AK5JPT_parton_Id;   //!
   TBranch        *b_jets_AK5JPT_parton_motherId;   //!
   TBranch        *b_jets_AK5JPT_parton_pt;   //!
   TBranch        *b_jets_AK5JPT_parton_phi;   //!
   TBranch        *b_jets_AK5JPT_parton_eta;   //!
   TBranch        *b_jets_AK5JPT_parton_Energy;   //!
   TBranch        *b_jets_AK5JPT_parton_mass;   //!
   TBranch        *b_jets_AK5JPT_gen_et;   //!
   TBranch        *b_jets_AK5JPT_gen_pt;   //!
   TBranch        *b_jets_AK5JPT_gen_eta;   //!
   TBranch        *b_jets_AK5JPT_gen_phi;   //!
   TBranch        *b_jets_AK5JPT_gen_mass;   //!
   TBranch        *b_jets_AK5JPT_gen_Energy;   //!
   TBranch        *b_jets_AK5JPT_gen_Id;   //!
   TBranch        *b_jets_AK5JPT_gen_motherID;   //!
   TBranch        *b_jets_AK5JPT_gen_threeCharge;   //!
   TBranch        *b_jets_AK5JPT_partonFlavour;   //!
   TBranch        *b_jets_AK5JPT_btag_TC_highPur;   //!
   TBranch        *b_jets_AK5JPT_btag_TC_highEff;   //!
   TBranch        *b_jets_AK5JPT_btag_jetProb;   //!
   TBranch        *b_jets_AK5JPT_btag_jetBProb;   //!
   TBranch        *b_jets_AK5JPT_btag_softEle;   //!
   TBranch        *b_jets_AK5JPT_btag_softMuon;   //!
   TBranch        *b_jets_AK5JPT_btag_secVertexHighPur;   //!
   TBranch        *b_jets_AK5JPT_btag_secVertexHighEff;   //!
   TBranch        *b_jets_AK5JPT_btag_secVertexCombined;   //!
   TBranch        *b_jets_AK5JPT_jetCharge;   //!
   TBranch        *b_jets_AK5JPT_chgEmE;   //!
   TBranch        *b_jets_AK5JPT_chgHadE;   //!
   TBranch        *b_jets_AK5JPT_chgMuE;   //!
   TBranch        *b_jets_AK5JPT_chg_Mult;   //!
   TBranch        *b_jets_AK5JPT_neutralEmE;   //!
   TBranch        *b_jets_AK5JPT_neutralHadE;   //!
   TBranch        *b_jets_AK5JPT_neutral_Mult;   //!
   TBranch        *b_jets_AK5JPT_mu_Mult;   //!
   TBranch        *b_jets_AK5JPT_emf;   //!
   TBranch        *b_jets_AK5JPT_ehf;   //!
   TBranch        *b_jets_AK5JPT_n60;   //!
   TBranch        *b_jets_AK5JPT_n90;   //!
   TBranch        *b_jets_AK5JPT_etaetaMoment;   //!
   TBranch        *b_jets_AK5JPT_etaphiMoment;   //!
   TBranch        *b_jets_AK5JPT_phiphiMoment;   //!
   TBranch        *b_jets_AK5JPT_n90Hits;   //!
   TBranch        *b_jets_AK5JPT_fHPD;   //!
   TBranch        *b_jets_AK5JPT_fRBX;   //!
   TBranch        *b_jets_AK5JPT_hitsInN90;   //!
   TBranch        *b_jets_AK5JPT_nECALTowers;   //!
   TBranch        *b_jets_AK5JPT_nHCALTowers;   //!
   TBranch        *b_jets_AK5JPT_fSubDetector1;   //!
   TBranch        *b_jets_AK5JPT_fSubDetector2;   //!
   TBranch        *b_jets_AK5JPT_fSubDetector3;   //!
   TBranch        *b_jets_AK5JPT_fSubDetector4;   //!
   TBranch        *b_jets_AK5JPT_area;   //!
   TBranch        *b_jets_AK5JPT_corrFactorRaw;   //!
   TBranch        *b_jets_AK5JPT_mass;   //!
   TBranch        *b_Njets_AK5PF;   //!
   TBranch        *b_jets_AK5PF_energy;   //!
   TBranch        *b_jets_AK5PF_et;   //!
   TBranch        *b_jets_AK5PF_eta;   //!
   TBranch        *b_jets_AK5PF_phi;   //!
   TBranch        *b_jets_AK5PF_pt;   //!
   TBranch        *b_jets_AK5PF_px;   //!
   TBranch        *b_jets_AK5PF_py;   //!
   TBranch        *b_jets_AK5PF_pz;   //!
   TBranch        *b_jets_AK5PF_status;   //!
   TBranch        *b_jets_AK5PF_theta;   //!
   TBranch        *b_jets_AK5PF_parton_Id;   //!
   TBranch        *b_jets_AK5PF_parton_motherId;   //!
   TBranch        *b_jets_AK5PF_parton_pt;   //!
   TBranch        *b_jets_AK5PF_parton_phi;   //!
   TBranch        *b_jets_AK5PF_parton_eta;   //!
   TBranch        *b_jets_AK5PF_parton_Energy;   //!
   TBranch        *b_jets_AK5PF_parton_mass;   //!
   TBranch        *b_jets_AK5PF_gen_et;   //!
   TBranch        *b_jets_AK5PF_gen_pt;   //!
   TBranch        *b_jets_AK5PF_gen_eta;   //!
   TBranch        *b_jets_AK5PF_gen_phi;   //!
   TBranch        *b_jets_AK5PF_gen_mass;   //!
   TBranch        *b_jets_AK5PF_gen_Energy;   //!
   TBranch        *b_jets_AK5PF_gen_Id;   //!
   TBranch        *b_jets_AK5PF_gen_motherID;   //!
   TBranch        *b_jets_AK5PF_gen_threeCharge;   //!
   TBranch        *b_jets_AK5PF_partonFlavour;   //!
   TBranch        *b_jets_AK5PF_btag_TC_highPur;   //!
   TBranch        *b_jets_AK5PF_btag_TC_highEff;   //!
   TBranch        *b_jets_AK5PF_btag_jetProb;   //!
   TBranch        *b_jets_AK5PF_btag_jetBProb;   //!
   TBranch        *b_jets_AK5PF_btag_softEle;   //!
   TBranch        *b_jets_AK5PF_btag_softMuon;   //!
   TBranch        *b_jets_AK5PF_btag_secVertexHighPur;   //!
   TBranch        *b_jets_AK5PF_btag_secVertexHighEff;   //!
   TBranch        *b_jets_AK5PF_btag_secVertexCombined;   //!
   TBranch        *b_jets_AK5PF_jetCharge;   //!
   TBranch        *b_jets_AK5PF_chgEmE;   //!
   TBranch        *b_jets_AK5PF_chgHadE;   //!
   TBranch        *b_jets_AK5PF_chgMuE;   //!
   TBranch        *b_jets_AK5PF_chg_Mult;   //!
   TBranch        *b_jets_AK5PF_neutralEmE;   //!
   TBranch        *b_jets_AK5PF_neutralHadE;   //!
   TBranch        *b_jets_AK5PF_neutral_Mult;   //!
   TBranch        *b_jets_AK5PF_mu_Mult;   //!
   TBranch        *b_jets_AK5PF_emf;   //!
   TBranch        *b_jets_AK5PF_ehf;   //!
   TBranch        *b_jets_AK5PF_n60;   //!
   TBranch        *b_jets_AK5PF_n90;   //!
   TBranch        *b_jets_AK5PF_etaetaMoment;   //!
   TBranch        *b_jets_AK5PF_etaphiMoment;   //!
   TBranch        *b_jets_AK5PF_phiphiMoment;   //!
   TBranch        *b_jets_AK5PF_n90Hits;   //!
   TBranch        *b_jets_AK5PF_fHPD;   //!
   TBranch        *b_jets_AK5PF_fRBX;   //!
   TBranch        *b_jets_AK5PF_hitsInN90;   //!
   TBranch        *b_jets_AK5PF_nECALTowers;   //!
   TBranch        *b_jets_AK5PF_nHCALTowers;   //!
   TBranch        *b_jets_AK5PF_fSubDetector1;   //!
   TBranch        *b_jets_AK5PF_fSubDetector2;   //!
   TBranch        *b_jets_AK5PF_fSubDetector3;   //!
   TBranch        *b_jets_AK5PF_fSubDetector4;   //!
   TBranch        *b_jets_AK5PF_area;   //!
   TBranch        *b_jets_AK5PF_corrFactorRaw;   //!
   TBranch        *b_jets_AK5PF_mass;   //!
   TBranch        *b_Njets_AK5PFclean;   //!
   TBranch        *b_jets_AK5PFclean_energy;   //!
   TBranch        *b_jets_AK5PFclean_et;   //!
   TBranch        *b_jets_AK5PFclean_eta;   //!
   TBranch        *b_jets_AK5PFclean_phi;   //!
   TBranch        *b_jets_AK5PFclean_pt;   //!
   TBranch        *b_jets_AK5PFclean_px;   //!
   TBranch        *b_jets_AK5PFclean_py;   //!
   TBranch        *b_jets_AK5PFclean_pz;   //!
   TBranch        *b_jets_AK5PFclean_status;   //!
   TBranch        *b_jets_AK5PFclean_theta;   //!
   TBranch        *b_jets_AK5PFclean_parton_Id;   //!
   TBranch        *b_jets_AK5PFclean_parton_motherId;   //!
   TBranch        *b_jets_AK5PFclean_parton_pt;   //!
   TBranch        *b_jets_AK5PFclean_parton_phi;   //!
   TBranch        *b_jets_AK5PFclean_parton_eta;   //!
   TBranch        *b_jets_AK5PFclean_parton_Energy;   //!
   TBranch        *b_jets_AK5PFclean_parton_mass;   //!
   TBranch        *b_jets_AK5PFclean_gen_et;   //!
   TBranch        *b_jets_AK5PFclean_gen_pt;   //!
   TBranch        *b_jets_AK5PFclean_gen_eta;   //!
   TBranch        *b_jets_AK5PFclean_gen_phi;   //!
   TBranch        *b_jets_AK5PFclean_gen_mass;   //!
   TBranch        *b_jets_AK5PFclean_gen_Energy;   //!
   TBranch        *b_jets_AK5PFclean_gen_Id;   //!
   TBranch        *b_jets_AK5PFclean_partonFlavour;   //!
   TBranch        *b_jets_AK5PFclean_btag_TC_highPur;   //!
   TBranch        *b_jets_AK5PFclean_btag_TC_highEff;   //!
   TBranch        *b_jets_AK5PFclean_btag_jetProb;   //!
   TBranch        *b_jets_AK5PFclean_btag_jetBProb;   //!
   TBranch        *b_jets_AK5PFclean_btag_softEle;   //!
   TBranch        *b_jets_AK5PFclean_btag_softMuon;   //!
   TBranch        *b_jets_AK5PFclean_btag_secVertexHighPur;   //!
   TBranch        *b_jets_AK5PFclean_btag_secVertexHighEff;   //!
   TBranch        *b_jets_AK5PFclean_btag_secVertexCombined;   //!
   TBranch        *b_jets_AK5PFclean_jetCharge;   //!
   TBranch        *b_jets_AK5PFclean_chgEmE;   //!
   TBranch        *b_jets_AK5PFclean_chgHadE;   //!
   TBranch        *b_jets_AK5PFclean_chgMuE;   //!
   TBranch        *b_jets_AK5PFclean_chg_Mult;   //!
   TBranch        *b_jets_AK5PFclean_neutralEmE;   //!
   TBranch        *b_jets_AK5PFclean_neutralHadE;   //!
   TBranch        *b_jets_AK5PFclean_neutral_Mult;   //!
   TBranch        *b_jets_AK5PFclean_mu_Mult;   //!
   TBranch        *b_jets_AK5PFclean_emf;   //!
   TBranch        *b_jets_AK5PFclean_ehf;   //!
   TBranch        *b_jets_AK5PFclean_n60;   //!
   TBranch        *b_jets_AK5PFclean_n90;   //!
   TBranch        *b_jets_AK5PFclean_etaetaMoment;   //!
   TBranch        *b_jets_AK5PFclean_etaphiMoment;   //!
   TBranch        *b_jets_AK5PFclean_phiphiMoment;   //!
   TBranch        *b_jets_AK5PFclean_n90Hits;   //!
   TBranch        *b_jets_AK5PFclean_fHPD;   //!
   TBranch        *b_jets_AK5PFclean_fRBX;   //!
   TBranch        *b_jets_AK5PFclean_hitsInN90;   //!
   TBranch        *b_jets_AK5PFclean_nECALTowers;   //!
   TBranch        *b_jets_AK5PFclean_nHCALTowers;   //!
   TBranch        *b_jets_AK5PFclean_fSubDetector1;   //!
   TBranch        *b_jets_AK5PFclean_fSubDetector2;   //!
   TBranch        *b_jets_AK5PFclean_fSubDetector3;   //!
   TBranch        *b_jets_AK5PFclean_fSubDetector4;   //!
   TBranch        *b_jets_AK5PFclean_area;   //!
   TBranch        *b_jets_AK5PFclean_corrFactorRaw;   //!
   TBranch        *b_jets_AK5PFclean_mass;   //!
   TBranch        *b_Nmc_doc;   //!
   TBranch        *b_mc_doc_id;   //!
   TBranch        *b_mc_doc_pt;   //!
   TBranch        *b_mc_doc_px;   //!
   TBranch        *b_mc_doc_py;   //!
   TBranch        *b_mc_doc_pz;   //!
   TBranch        *b_mc_doc_eta;   //!
   TBranch        *b_mc_doc_phi;   //!
   TBranch        *b_mc_doc_theta;   //!
   TBranch        *b_mc_doc_energy;   //!
   TBranch        *b_mc_doc_status;   //!
   TBranch        *b_mc_doc_charge;   //!
   TBranch        *b_mc_doc_mother_id;   //!
   TBranch        *b_mc_doc_grandmother_id;   //!
   TBranch        *b_mc_doc_ggrandmother_id;   //!
   TBranch        *b_mc_doc_mother_pt;   //!
   TBranch        *b_mc_doc_vertex_x;   //!
   TBranch        *b_mc_doc_vertex_y;   //!
   TBranch        *b_mc_doc_vertex_z;   //!
   TBranch        *b_mc_doc_mass;   //!
   TBranch        *b_mc_doc_numOfDaughters;   //!
   TBranch        *b_mc_doc_numOfMothers;   //!
   TBranch        *b_Nmc_electrons;   //!
   TBranch        *b_mc_electrons_id;   //!
   TBranch        *b_mc_electrons_pt;   //!
   TBranch        *b_mc_electrons_px;   //!
   TBranch        *b_mc_electrons_py;   //!
   TBranch        *b_mc_electrons_pz;   //!
   TBranch        *b_mc_electrons_eta;   //!
   TBranch        *b_mc_electrons_phi;   //!
   TBranch        *b_mc_electrons_theta;   //!
   TBranch        *b_mc_electrons_status;   //!
   TBranch        *b_mc_electrons_energy;   //!
   TBranch        *b_mc_electrons_charge;   //!
   TBranch        *b_mc_electrons_mother_id;   //!
   TBranch        *b_mc_electrons_mother_pt;   //!
   TBranch        *b_mc_electrons_grandmother_id;   //!
   TBranch        *b_mc_electrons_ggrandmother_id;   //!
   TBranch        *b_mc_electrons_vertex_x;   //!
   TBranch        *b_mc_electrons_vertex_y;   //!
   TBranch        *b_mc_electrons_vertex_z;   //!
   TBranch        *b_mc_electrons_mass;   //!
   TBranch        *b_mc_electrons_numOfDaughters;   //!
   TBranch        *b_Nmc_mus;   //!
   TBranch        *b_mc_mus_id;   //!
   TBranch        *b_mc_mus_pt;   //!
   TBranch        *b_mc_mus_px;   //!
   TBranch        *b_mc_mus_py;   //!
   TBranch        *b_mc_mus_pz;   //!
   TBranch        *b_mc_mus_eta;   //!
   TBranch        *b_mc_mus_phi;   //!
   TBranch        *b_mc_mus_theta;   //!
   TBranch        *b_mc_mus_status;   //!
   TBranch        *b_mc_mus_energy;   //!
   TBranch        *b_mc_mus_charge;   //!
   TBranch        *b_mc_mus_mother_id;   //!
   TBranch        *b_mc_mus_mother_pt;   //!
   TBranch        *b_mc_mus_grandmother_id;   //!
   TBranch        *b_mc_mus_ggrandmother_id;   //!
   TBranch        *b_mc_mus_vertex_x;   //!
   TBranch        *b_mc_mus_vertex_y;   //!
   TBranch        *b_mc_mus_vertex_z;   //!
   TBranch        *b_mc_mus_mass;   //!
   TBranch        *b_mc_mus_numOfDaughters;   //!
   TBranch        *b_Nmc_photons;   //!
   TBranch        *b_mc_photons_id;   //!
   TBranch        *b_mc_photons_pt;   //!
   TBranch        *b_mc_photons_px;   //!
   TBranch        *b_mc_photons_py;   //!
   TBranch        *b_mc_photons_pz;   //!
   TBranch        *b_mc_photons_eta;   //!
   TBranch        *b_mc_photons_phi;   //!
   TBranch        *b_mc_photons_theta;   //!
   TBranch        *b_mc_photons_status;   //!
   TBranch        *b_mc_photons_energy;   //!
   TBranch        *b_mc_photons_charge;   //!
   TBranch        *b_mc_photons_mother_id;   //!
   TBranch        *b_mc_photons_mother_pt;   //!
   TBranch        *b_mc_photons_grandmother_id;   //!
   TBranch        *b_mc_photons_ggrandmother_id;   //!
   TBranch        *b_mc_photons_vertex_x;   //!
   TBranch        *b_mc_photons_vertex_y;   //!
   TBranch        *b_mc_photons_vertex_z;   //!
   TBranch        *b_mc_photons_mass;   //!
   TBranch        *b_mc_photons_numOfDaughters;   //!
   TBranch        *b_Nmc_taus;   //!
   TBranch        *b_mc_taus_id;   //!
   TBranch        *b_mc_taus_pt;   //!
   TBranch        *b_mc_taus_px;   //!
   TBranch        *b_mc_taus_py;   //!
   TBranch        *b_mc_taus_pz;   //!
   TBranch        *b_mc_taus_eta;   //!
   TBranch        *b_mc_taus_phi;   //!
   TBranch        *b_mc_taus_theta;   //!
   TBranch        *b_mc_taus_status;   //!
   TBranch        *b_mc_taus_energy;   //!
   TBranch        *b_mc_taus_charge;   //!
   TBranch        *b_mc_taus_mother_id;   //!
   TBranch        *b_mc_taus_mother_pt;   //!
   TBranch        *b_mc_taus_grandmother_id;   //!
   TBranch        *b_mc_taus_ggrandmother_id;   //!
   TBranch        *b_mc_taus_vertex_x;   //!
   TBranch        *b_mc_taus_vertex_y;   //!
   TBranch        *b_mc_taus_vertex_z;   //!
   TBranch        *b_mc_taus_mass;   //!
   TBranch        *b_mc_taus_numOfDaughters;   //!
   TBranch        *b_Nmets_AK5;   //!
   TBranch        *b_mets_AK5_et;   //!
   TBranch        *b_mets_AK5_phi;   //!
   TBranch        *b_mets_AK5_ex;   //!
   TBranch        *b_mets_AK5_ey;   //!
   TBranch        *b_mets_AK5_gen_et;   //!
   TBranch        *b_mets_AK5_gen_phi;   //!
   TBranch        *b_mets_AK5_sign;   //!
   TBranch        *b_mets_AK5_sumEt;   //!
   TBranch        *b_mets_AK5_unCPhi;   //!
   TBranch        *b_mets_AK5_unCPt;   //!
   TBranch        *b_Nmulti5x5EBC;   //!
   TBranch        *b_multi5x5EBC_energy;   //!
   TBranch        *b_multi5x5EBC_x;   //!
   TBranch        *b_multi5x5EBC_y;   //!
   TBranch        *b_multi5x5EBC_z;   //!
   TBranch        *b_multi5x5EBC_rho;   //!
   TBranch        *b_multi5x5EBC_phi;   //!
   TBranch        *b_multi5x5EBC_eta;   //!
   TBranch        *b_multi5x5EBC_theta;   //!
   TBranch        *b_Nmus;   //!
   TBranch        *b_pf_mus_energy;   //!
   TBranch        *b_pf_mus_et;   //!
   TBranch        *b_pf_mus_eta;   //!
   TBranch        *b_pf_mus_phi;   //!
   TBranch        *b_pf_mus_pt;   //!
   TBranch        *b_pf_mus_px;   //!
   TBranch        *b_pf_mus_py;   //!
   TBranch        *b_pf_mus_pz;   //!
   TBranch        *b_pf_mus_status;   //!
   TBranch        *b_pf_mus_theta;   //!
   TBranch        *b_mus_eta;	//!
   TBranch        *b_mus_phi;	//!
   TBranch        *b_mus_pt;   //!
   TBranch        *b_mus_px;   //!   
   TBranch        *b_mus_py;   //!   
   TBranch        *b_mus_pz;   //!
   TBranch        *b_mus_charge;   //!   
   TBranch        *b_mus_energy;   //!
   TBranch        *b_mus_gen_id;   //!
   TBranch        *b_mus_gen_phi;   //!
   TBranch        *b_mus_gen_pt;   //!
   TBranch        *b_mus_gen_pz;   //!
   TBranch        *b_mus_gen_px;   //!
   TBranch        *b_mus_gen_py;   //!
   TBranch        *b_mus_gen_eta;   //!
   TBranch        *b_mus_gen_theta;   //!
   TBranch        *b_mus_gen_et;   //!
   TBranch        *b_mus_gen_mother_id;   //!
   TBranch        *b_mus_gen_mother_phi;   //!
   TBranch        *b_mus_gen_mother_pt;   //!
   TBranch        *b_mus_gen_mother_pz;   //!
   TBranch        *b_mus_gen_mother_px;   //!
   TBranch        *b_mus_gen_mother_py;   //!
   TBranch        *b_mus_gen_mother_eta;   //!
   TBranch        *b_mus_gen_mother_theta;   //!
   TBranch        *b_mus_gen_mother_et;   //!
   TBranch        *b_pf_mus_tkHits;   //!
   TBranch        *b_pf_mus_cIso;   //!
   TBranch        *b_pf_mus_tIso;   //!
   TBranch        *b_pf_mus_ecalIso;   //!
   TBranch        *b_pf_mus_hcalIso;   //!
   TBranch        *b_pf_mus_ecalvetoDep;   //!
   TBranch        *b_pf_mus_hcalvetoDep;   //!
   TBranch        *b_pf_mus_calEnergyEm;   //!
   TBranch        *b_pf_mus_calEnergyHad;   //!
   TBranch        *b_pf_mus_calEnergyHo;   //!
   TBranch        *b_pf_mus_calEnergyEmS9;   //!
   TBranch        *b_pf_mus_calEnergyHadS9;   //!
   TBranch        *b_pf_mus_calEnergyHoS9;   //!
   TBranch        *b_pf_mus_iso03_sumPt;   //!
   TBranch        *b_pf_mus_iso03_emEt;   //!
   TBranch        *b_pf_mus_iso03_hadEt;   //!
   TBranch        *b_pf_mus_iso03_hoEt;   //!
   TBranch        *b_pf_mus_iso03_nTracks;   //!
   TBranch        *b_pf_mus_iso05_sumPt;   //!
   TBranch        *b_pf_mus_iso05_emEt;   //!
   TBranch        *b_pf_mus_iso05_hadEt;   //!
   TBranch        *b_pf_mus_iso05_hoEt;   //!
   TBranch        *b_pf_mus_iso05_nTracks;   //!
   TBranch        *b_pf_mus_chargedHadronIso;   //!
   TBranch        *b_pf_mus_photonIso;   //!
   TBranch        *b_pf_mus_neutralHadronIso;   //!
   TBranch        *b_pf_mus_charge;   //!
   TBranch        *b_pf_mus_cm_chi2;   //!
   TBranch        *b_pf_mus_cm_ndof;   //!
   TBranch        *b_pf_mus_cm_chg;   //!
   TBranch        *b_pf_mus_cm_pt;   //!
   TBranch        *b_pf_mus_cm_px;   //!
   TBranch        *b_pf_mus_cm_py;   //!
   TBranch        *b_pf_mus_cm_pz;   //!
   TBranch        *b_pf_mus_cm_eta;   //!
   TBranch        *b_pf_mus_cm_phi;   //!
   TBranch        *b_pf_mus_cm_theta;   //!
   TBranch        *b_pf_mus_cm_d0dum;   //!
   TBranch        *b_pf_mus_cm_dz;   //!
   TBranch        *b_pf_mus_cm_vx;   //!
   TBranch        *b_pf_mus_cm_vy;   //!
   TBranch        *b_pf_mus_cm_vz;   //!
   TBranch        *b_pf_mus_cm_numvalhits;   //!
   TBranch        *b_pf_mus_cm_numlosthits;   //!
   TBranch        *b_pf_mus_cm_numvalMuonhits;   //!
   TBranch        *b_pf_mus_cm_d0dumErr;   //!
   TBranch        *b_pf_mus_cm_dzErr;   //!
   TBranch        *b_pf_mus_cm_ptErr;   //!
   TBranch        *b_pf_mus_cm_etaErr;   //!
   TBranch        *b_pf_mus_cm_phiErr;   //!
   TBranch        *b_pf_mus_tk_id;   //!
   TBranch        *b_pf_mus_tk_chi2;   //!
   TBranch        *b_pf_mus_tk_ndof;   //!
   TBranch        *b_pf_mus_tk_chg;   //!
   TBranch        *b_pf_mus_tk_pt;   //!
   TBranch        *b_pf_mus_tk_px;   //!
   TBranch        *b_pf_mus_tk_py;   //!
   TBranch        *b_pf_mus_tk_pz;   //!
   TBranch        *b_pf_mus_tk_eta;   //!
   TBranch        *b_pf_mus_tk_phi;   //!
   TBranch        *b_pf_mus_tk_theta;   //!
   TBranch        *b_pf_mus_tk_d0dum;   //!
   TBranch        *b_pf_mus_tk_dz;   //!
   TBranch        *b_pf_mus_tk_vx;   //!
   TBranch        *b_pf_mus_tk_vy;   //!
   TBranch        *b_pf_mus_tk_vz;   //!
   TBranch        *b_pf_mus_tk_numvalhits;   //!
   TBranch        *b_pf_mus_tk_numlosthits;   //!
   TBranch        *b_pf_mus_tk_d0dumErr;   //!
   TBranch        *b_pf_mus_tk_dzErr;   //!
   TBranch        *b_pf_mus_tk_ptErr;   //!
   TBranch        *b_pf_mus_tk_etaErr;   //!
   TBranch        *b_pf_mus_tk_phiErr;   //!
   TBranch        *b_pf_mus_tk_numvalPixelhits;   //!
   TBranch        *b_pf_mus_stamu_chi2;   //!
   TBranch        *b_pf_mus_stamu_ndof;   //!
   TBranch        *b_pf_mus_stamu_chg;   //!
   TBranch        *b_pf_mus_stamu_pt;   //!
   TBranch        *b_pf_mus_stamu_px;   //!
   TBranch        *b_pf_mus_stamu_py;   //!
   TBranch        *b_pf_mus_stamu_pz;   //!
   TBranch        *b_pf_mus_stamu_eta;   //!
   TBranch        *b_pf_mus_stamu_phi;   //!
   TBranch        *b_pf_mus_stamu_theta;   //!
   TBranch        *b_pf_mus_stamu_d0dum;   //!
   TBranch        *b_pf_mus_stamu_dz;   //!
   TBranch        *b_pf_mus_stamu_vx;   //!
   TBranch        *b_pf_mus_stamu_vy;   //!
   TBranch        *b_pf_mus_stamu_vz;   //!
   TBranch        *b_pf_mus_stamu_numvalhits;   //!
   TBranch        *b_pf_mus_stamu_numlosthits;   //!
   TBranch        *b_pf_mus_stamu_d0dumErr;   //!
   TBranch        *b_pf_mus_stamu_dzErr;   //!
   TBranch        *b_pf_mus_stamu_ptErr;   //!
   TBranch        *b_pf_mus_stamu_etaErr;   //!
   TBranch        *b_pf_mus_stamu_phiErr;   //!
   TBranch        *b_pf_mus_num_matches;   //!
   TBranch        *b_pf_mus_isTrackerMuon;   //!
   TBranch        *b_pf_mus_isStandAloneMuon;   //!
   TBranch        *b_pf_mus_isCaloMuon;   //!
   TBranch        *b_pf_mus_isGlobalMuon;   //!
   TBranch        *b_pf_mus_isElectron;   //!
   TBranch        *b_pf_mus_isConvertedPhoton;   //!
   TBranch        *b_pf_mus_isPhoton;   //!
   TBranch        *b_pf_mus_id_All;   //!
   TBranch        *b_pf_mus_id_AllGlobalMuons;   //!
   TBranch        *b_pf_mus_id_AllStandAloneMuons;   //!
   TBranch        *b_pf_mus_id_AllTrackerMuons;   //!
   TBranch        *b_pf_mus_id_TrackerMuonArbitrated;   //!
   TBranch        *b_pf_mus_id_AllArbitrated;   //!
   TBranch        *b_pf_mus_id_GlobalMuonPromptTight;   //!
   TBranch        *b_pf_mus_id_TMLastStationLoose;   //!
   TBranch        *b_pf_mus_id_TMLastStationTight;   //!
   TBranch        *b_pf_mus_id_TM2DCompatibilityLoose;   //!
   TBranch        *b_pf_mus_id_TM2DCompatibilityTight;   //!
   TBranch        *b_pf_mus_id_TMOneStationLoose;   //!
   TBranch        *b_pf_mus_id_TMOneStationTight;   //!
   TBranch        *b_pf_mus_id_TMLastStationOptimizedLowPtLoose;   //!
   TBranch        *b_pf_mus_id_TMLastStationOptimizedLowPtTight;   //!
   TBranch        *b_NpfTypeImets;   //!
   TBranch        *b_pfTypeImets_et;   //!
   TBranch        *b_pfTypeImets_phi;   //!
   TBranch        *b_pfTypeImets_ex;   //!
   TBranch        *b_pfTypeImets_ey;   //!
   TBranch        *b_pfTypeImets_gen_et;   //!
   TBranch        *b_pfTypeImets_gen_phi;   //!
   TBranch        *b_pfTypeImets_sign;   //!
   TBranch        *b_pfTypeImets_sumEt;   //!
   TBranch        *b_pfTypeImets_unCPhi;   //!
   TBranch        *b_pfTypeImets_unCPt;   //!
   TBranch        *b_Npf_els;   //!
   TBranch        *b_Npf_mus;   //!
   TBranch        *b_pf_mus_gen_id;   //!
   TBranch        *b_pf_mus_gen_phi;   //!
   TBranch        *b_pf_mus_gen_pt;   //!
   TBranch        *b_pf_mus_gen_pz;   //!
   TBranch        *b_pf_mus_gen_px;   //!
   TBranch        *b_pf_mus_gen_py;   //!
   TBranch        *b_pf_mus_gen_eta;   //!
   TBranch        *b_pf_mus_gen_theta;   //!
   TBranch        *b_pf_mus_gen_et;   //!
   TBranch        *b_pf_mus_gen_mother_id;   //!
   TBranch        *b_pf_mus_gen_mother_phi;   //!
   TBranch        *b_pf_mus_gen_mother_pt;   //!
   TBranch        *b_pf_mus_gen_mother_pz;   //!
   TBranch        *b_pf_mus_gen_mother_px;   //!
   TBranch        *b_pf_mus_gen_mother_py;   //!
   TBranch        *b_pf_mus_gen_mother_eta;   //!
   TBranch        *b_pf_mus_gen_mother_theta;   //!
   TBranch        *b_pf_mus_gen_mother_et;   //!
   TBranch        *b_Npfcand_els;   //!
   TBranch        *b_pf_els_particleId;   //!
   TBranch        *b_Npfcand_mus;   //!
   TBranch        *b_pf_mus_particleId;   //!
   TBranch        *b_Npfmets;   //!
   TBranch        *b_pfmets_et;   //!
   TBranch        *b_pfmets_phi;   //!
   TBranch        *b_pfmets_ex;   //!
   TBranch        *b_pfmets_ey;   //!
   TBranch        *b_pfmets_gen_et;   //!
   TBranch        *b_pfmets_gen_phi;   //!
   TBranch        *b_pfmets_sign;   //!
   TBranch        *b_pfmets_sumEt;   //!
   TBranch        *b_pfmets_unCPhi;   //!
   TBranch        *b_pfmets_unCPt;   //!
   TBranch        *b_Nphotons;   //!
   TBranch        *b_photons_energy;   //!
   TBranch        *b_photons_et;   //!
   TBranch        *b_photons_eta;   //!
   TBranch        *b_photons_phi;   //!
   TBranch        *b_photons_pt;   //!
   TBranch        *b_photons_px;   //!
   TBranch        *b_photons_py;   //!
   TBranch        *b_photons_pz;   //!
   TBranch        *b_photons_status;   //!
   TBranch        *b_photons_theta;   //!
   TBranch        *b_photons_hadOverEM;   //!
   TBranch        *b_photons_scEnergy;   //!
   TBranch        *b_photons_scRawEnergy;   //!
   TBranch        *b_photons_scEta;   //!
   TBranch        *b_photons_scPhi;   //!
   TBranch        *b_photons_scEtaWidth;   //!
   TBranch        *b_photons_scPhiWidth;   //!
   TBranch        *b_photons_tIso;   //!
   TBranch        *b_photons_ecalIso;   //!
   TBranch        *b_photons_hcalIso;   //!
   TBranch        *b_photons_isoEcalRecHitDR04;   //!
   TBranch        *b_photons_isoHcalRecHitDR04;   //!
   TBranch        *b_photons_isoSolidTrkConeDR04;   //!
   TBranch        *b_photons_isoHollowTrkConeDR04;   //!
   TBranch        *b_photons_nTrkSolidConeDR04;   //!
   TBranch        *b_photons_nTrkHollowConeDR04;   //!
   TBranch        *b_photons_isoEcalRecHitDR03;   //!
   TBranch        *b_photons_isoHcalRecHitDR03;   //!
   TBranch        *b_photons_isoSolidTrkConeDR03;   //!
   TBranch        *b_photons_isoHollowTrkConeDR03;   //!
   TBranch        *b_photons_nTrkSolidConeDR03;   //!
   TBranch        *b_photons_nTrkHollowConeDR03;   //!
   TBranch        *b_photons_isAlsoElectron;   //!
   TBranch        *b_photons_hasPixelSeed;   //!
   TBranch        *b_photons_isConverted;   //!
   TBranch        *b_photons_isEBGap;   //!
   TBranch        *b_photons_isEEGap;   //!
   TBranch        *b_photons_isEBEEGap;   //!
   TBranch        *b_photons_isEBPho;   //!
   TBranch        *b_photons_isEEPho;   //!
   TBranch        *b_photons_isLoosePhoton;   //!
   TBranch        *b_photons_isTightPhoton;   //!
   TBranch        *b_photons_maxEnergyXtal;   //!
   TBranch        *b_photons_e1x5;   //!
   TBranch        *b_photons_e2x5;   //!
   TBranch        *b_photons_e3x3;   //!
   TBranch        *b_photons_e5x5;   //!
   TBranch        *b_photons_sigmaEtaEta;   //!
   TBranch        *b_photons_sigmaIetaIeta;   //!
   TBranch        *b_photons_r9;   //!
   TBranch        *b_photons_gen_et;   //!
   TBranch        *b_photons_gen_eta;   //!
   TBranch        *b_photons_gen_phi;   //!
   TBranch        *b_photons_gen_id;   //!
   TBranch        *b_Npv;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_xErr;   //!
   TBranch        *b_pv_yErr;   //!
   TBranch        *b_pv_zErr;   //!
   TBranch        *b_pv_chi2;   //!
   TBranch        *b_pv_ndof;   //!
   TBranch        *b_pv_isFake;   //!
   TBranch        *b_pv_tracksSize;   //!
   TBranch        *b_Ntaus;   //!
   TBranch        *b_taus_energy;   //!
   TBranch        *b_taus_et;   //!
   TBranch        *b_taus_eta;   //!
   TBranch        *b_taus_phi;   //!
   TBranch        *b_taus_pt;   //!
   TBranch        *b_taus_px;   //!
   TBranch        *b_taus_py;   //!
   TBranch        *b_taus_pz;   //!
   TBranch        *b_taus_status;   //!
   TBranch        *b_taus_theta;   //!
   TBranch        *b_taus_charge;   //!
   TBranch        *b_taus_emf;   //!
   TBranch        *b_taus_hcalTotOverPLead;   //!
   TBranch        *b_taus_hcalMaxOverPLead;   //!
   TBranch        *b_taus_hcal3x3OverPLead;   //!
   TBranch        *b_taus_ecalStripSumEOverPLead;   //!
   TBranch        *b_taus_elecPreIdOutput;   //!
   TBranch        *b_taus_elecPreIdDecision;   //!
   TBranch        *b_taus_leadPFChargedHadrCand_pt;   //!
   TBranch        *b_taus_leadPFChargedHadrCand_charge;   //!
   TBranch        *b_taus_leadPFChargedHadrCand_eta;   //!
   TBranch        *b_taus_leadPFChargedHadrCand_ECAL_eta;   //!
   TBranch        *b_taus_leadPFChargedHadrCand_phi;   //!
   TBranch        *b_taus_isoPFGammaCandsEtSum;   //!
   TBranch        *b_taus_isoPFChargedHadrCandsPtSum;   //!
   TBranch        *b_taus_leadingTrackFinding;   //!
   TBranch        *b_taus_leadingTrackPtCut;   //!
   TBranch        *b_taus_trackIsolation;   //!
   TBranch        *b_taus_ecalIsolation;   //!
   TBranch        *b_taus_byIsolation;   //!
   TBranch        *b_taus_againstElectron;   //!
   TBranch        *b_taus_againstMuon;   //!
   TBranch        *b_taus_taNC_quarter;   //!
   TBranch        *b_taus_taNC_one;   //!
   TBranch        *b_taus_taNC_half;   //!
   TBranch        *b_taus_taNC_tenth;   //!
   TBranch        *b_taus_taNC;   //!
   TBranch        *b_taus_byIsoUsingLeadingPi;   //!
   TBranch        *b_taus_tkIsoUsingLeadingPi;   //!
   TBranch        *b_taus_ecalIsoUsingLeadingPi;   //!
   TBranch        *b_taus_signalPFChargedHadrCandsSize;   //!
   TBranch        *b_taus_muDecision;   //!
   TBranch        *b_taus_Nprongs;   //!
   TBranch        *b_Ntcmets;   //!
   TBranch        *b_tcmets_et;   //!
   TBranch        *b_tcmets_phi;   //!
   TBranch        *b_tcmets_ex;   //!
   TBranch        *b_tcmets_ey;   //!
   TBranch        *b_tcmets_sumEt;   //!
   TBranch        *b_Ntracks;   //!
   TBranch        *b_tracks_chi2;   //!
   TBranch        *b_tracks_trkExptHitsInner;   //!
   TBranch        *b_tracks_trkExptHitsOuter;   //!
   TBranch        *b_tracks_trks_nlayerslost;   //!
   TBranch        *b_tracks_trks_nlayers;   //!
   TBranch        *b_tracks_trksvalidpixelhits;   //!
   TBranch        *b_tracks_trkslostpixelhits;   //!
   TBranch        *b_tracks_ndof;   //!
   TBranch        *b_tracks_chg;   //!
   TBranch        *b_tracks_pt;   //!
   TBranch        *b_tracks_px;   //!
   TBranch        *b_tracks_py;   //!
   TBranch        *b_tracks_pz;   //!
   TBranch        *b_tracks_eta;   //!
   TBranch        *b_tracks_phi;   //!
   TBranch        *b_tracks_theta;   //!
   TBranch        *b_tracks_d0dum;   //!
   TBranch        *b_tracks_dz;   //!
   TBranch        *b_tracks_vx;   //!
   TBranch        *b_tracks_vy;   //!
   TBranch        *b_tracks_vz;   //!
   TBranch        *b_tracks_numvalhits;   //!
   TBranch        *b_tracks_numlosthits;   //!
   TBranch        *b_tracks_d0dumErr;   //!
   TBranch        *b_tracks_dzErr;   //!
   TBranch        *b_tracks_ptErr;   //!
   TBranch        *b_tracks_etaErr;   //!
   TBranch        *b_tracks_phiErr;   //!
   TBranch        *b_tracks_Nrechits;   //!
   TBranch        *b_tracks_innerHitX;   //!
   TBranch        *b_tracks_innerHitY;   //!
   TBranch        *b_tracks_innerHitZ;   //!
   TBranch        *b_tracks_outerHitX;   //!
   TBranch        *b_tracks_outerHitY;   //!
   TBranch        *b_tracks_outerHitZ;   //!
   TBranch        *b_tracks_outerPx;   //!
   TBranch        *b_tracks_outerPy;   //!
   TBranch        *b_tracks_outerPz;   //!
   TBranch        *b_tracks_algo;   //!
   TBranch        *b_tracks_highPurity;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_experimentType;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_orbitNumber;   //!


   // Declaration of leaf types
   Double_t        AlCa_EcalEta_8E29;
   Double_t        AlCa_EcalPhiSym;
   Double_t        AlCa_EcalPi0_8E29;
   Double_t        AlCa_HcalPhiSym;
   Double_t        AlCa_RPCMuonNoHits;
   Double_t        AlCa_RPCMuonNormalisation;
   Double_t        HLTAnalyzerEndpath;
   Double_t        HLT_BTagIP_Jet50U;
   Double_t        HLT_BTagMu_Jet10U;
   Double_t        HLT_BackwardBSC;
   Double_t        HLT_CSCBeamHalo;
   Double_t        HLT_CSCBeamHaloOverlapRing1;
   Double_t        HLT_CSCBeamHaloOverlapRing2;
   Double_t        HLT_CSCBeamHaloRing2or3;
   Double_t        HLT_DiJetAve15U_8E29;
   Double_t        HLT_DiJetAve30U_8E29;
   Double_t        HLT_DiJetAve50U;
   Double_t        HLT_DoubleEle5_SW_L1R;
   Double_t        HLT_DoubleLooseIsoTau15;
   Double_t        HLT_DoubleMu0;
   Double_t        HLT_DoubleMu3;
   Double_t        HLT_DoublePhoton10_L1R;
   Double_t        HLT_DoublePhoton5_Jpsi_L1R;
   Double_t        HLT_DoublePhoton5_Upsilon_L1R;
   Double_t        HLT_DoublePhoton5_eeRes_L1R;
   Double_t        HLT_Ele10_LW_EleId_L1R;
   Double_t        HLT_Ele10_LW_L1R;
   Double_t        HLT_Ele15_LW_L1R;
   Double_t        HLT_Ele15_SC10_LW_L1R;
   Double_t        HLT_Ele15_SiStrip_L1R;
   Double_t        HLT_Ele20_LW_L1R;
   Double_t        HLT_ForwardBSC;
   Double_t        HLT_FwdJet20U;
   Double_t        HLT_HT100U;
   Double_t        HLT_IsoMu3;
   Double_t        HLT_IsoTrack_8E29;
   Double_t        HLT_Jet100U;
   Double_t        HLT_Jet15U;
   Double_t        HLT_Jet30U;
   Double_t        HLT_Jet50U;
   Double_t        HLT_Jet70U;
   Double_t        HLT_L1DoubleEG5;
   Double_t        HLT_L1DoubleMuOpen;
   Double_t        HLT_L1Jet6U;
   Double_t        HLT_L1MET20;
   Double_t        HLT_L1Mu;
   Double_t        HLT_L1Mu14_L1ETM30;
   Double_t        HLT_L1Mu14_L1SingleEG10;
   Double_t        HLT_L1Mu14_L1SingleJet6U;
   Double_t        HLT_L1Mu20;
   Double_t        HLT_L1MuOpen;
   Double_t        HLT_L1SingleEG5;
   Double_t        HLT_L1SingleEG8;
   Double_t        HLT_L2Mu11;
   Double_t        HLT_L2Mu9;
   Double_t        HLT_MET100;
   Double_t        HLT_MET45;
   Double_t        HLT_MinBiasEcal;
   Double_t        HLT_MinBiasHcal;
   Double_t        HLT_MinBiasPixel;
   Double_t        HLT_MinBiasPixel_Trk5;
   Double_t        HLT_Mu3;
   Double_t        HLT_Mu5;
   Double_t        HLT_Mu9;
   Double_t        HLT_Photon10_L1R;
   Double_t        HLT_Photon15_L1R;
   Double_t        HLT_Photon15_LooseEcalIso_L1R;
   Double_t        HLT_Photon15_TrackIso_L1R;
   Double_t        HLT_Photon20_L1R;
   Double_t        HLT_Photon30_L1R_8E29;
   Double_t        HLT_QuadJet15U;
   Double_t        HLT_SingleLooseIsoTau20;
   Double_t        HLT_StoppedHSCP_8E29;
   Double_t        HLT_TrackerCosmics;
   Double_t        HLT_ZeroBias;
   Double_t        HLTriggerFinalPath;
   Double_t        HLTriggerFirstPath;
   Double_t        L1Bit_0;
   Double_t        L1Bit_1;
   Double_t        L1Bit_10;
   Double_t        L1Bit_100;
   Double_t        L1Bit_101;
   Double_t        L1Bit_102;
   Double_t        L1Bit_103;
   Double_t        L1Bit_104;
   Double_t        L1Bit_105;
   Double_t        L1Bit_106;
   Double_t        L1Bit_107;
   Double_t        L1Bit_108;
   Double_t        L1Bit_109;
   Double_t        L1Bit_11;
   Double_t        L1Bit_110;
   Double_t        L1Bit_111;
   Double_t        L1Bit_112;
   Double_t        L1Bit_113;
   Double_t        L1Bit_114;
   Double_t        L1Bit_115;
   Double_t        L1Bit_116;
   Double_t        L1Bit_117;
   Double_t        L1Bit_118;
   Double_t        L1Bit_119;
   Double_t        L1Bit_12;
   Double_t        L1Bit_120;
   Double_t        L1Bit_121;
   Double_t        L1Bit_122;
   Double_t        L1Bit_123;
   Double_t        L1Bit_124;
   Double_t        L1Bit_125;
   Double_t        L1Bit_126;
   Double_t        L1Bit_127;
   Double_t        L1Bit_13;
   Double_t        L1Bit_14;
   Double_t        L1Bit_15;
   Double_t        L1Bit_16;
   Double_t        L1Bit_17;
   Double_t        L1Bit_18;
   Double_t        L1Bit_19;
   Double_t        L1Bit_2;
   Double_t        L1Bit_20;
   Double_t        L1Bit_21;
   Double_t        L1Bit_22;
   Double_t        L1Bit_23;
   Double_t        L1Bit_24;
   Double_t        L1Bit_25;
   Double_t        L1Bit_26;
   Double_t        L1Bit_27;
   Double_t        L1Bit_28;
   Double_t        L1Bit_29;
   Double_t        L1Bit_3;
   Double_t        L1Bit_30;
   Double_t        L1Bit_31;
   Double_t        L1Bit_32;
   Double_t        L1Bit_33;
   Double_t        L1Bit_34;
   Double_t        L1Bit_35;
   Double_t        L1Bit_36;
   Double_t        L1Bit_37;
   Double_t        L1Bit_38;
   Double_t        L1Bit_39;
   Double_t        L1Bit_4;
   Double_t        L1Bit_40;
   Double_t        L1Bit_41;
   Double_t        L1Bit_42;
   Double_t        L1Bit_43;
   Double_t        L1Bit_44;
   Double_t        L1Bit_45;
   Double_t        L1Bit_46;
   Double_t        L1Bit_47;
   Double_t        L1Bit_48;
   Double_t        L1Bit_49;
   Double_t        L1Bit_5;
   Double_t        L1Bit_50;
   Double_t        L1Bit_51;
   Double_t        L1Bit_52;
   Double_t        L1Bit_53;
   Double_t        L1Bit_54;
   Double_t        L1Bit_55;
   Double_t        L1Bit_56;
   Double_t        L1Bit_57;
   Double_t        L1Bit_58;
   Double_t        L1Bit_59;
   Double_t        L1Bit_6;
   Double_t        L1Bit_60;
   Double_t        L1Bit_61;
   Double_t        L1Bit_62;
   Double_t        L1Bit_63;
   Double_t        L1Bit_64;
   Double_t        L1Bit_65;
   Double_t        L1Bit_66;
   Double_t        L1Bit_67;
   Double_t        L1Bit_68;
   Double_t        L1Bit_69;
   Double_t        L1Bit_7;
   Double_t        L1Bit_70;
   Double_t        L1Bit_71;
   Double_t        L1Bit_72;
   Double_t        L1Bit_73;
   Double_t        L1Bit_74;
   Double_t        L1Bit_75;
   Double_t        L1Bit_76;
   Double_t        L1Bit_77;
   Double_t        L1Bit_78;
   Double_t        L1Bit_79;
   Double_t        L1Bit_8;
   Double_t        L1Bit_80;
   Double_t        L1Bit_81;
   Double_t        L1Bit_82;
   Double_t        L1Bit_83;
   Double_t        L1Bit_84;
   Double_t        L1Bit_85;
   Double_t        L1Bit_86;
   Double_t        L1Bit_87;
   Double_t        L1Bit_88;
   Double_t        L1Bit_89;
   Double_t        L1Bit_9;
   Double_t        L1Bit_90;
   Double_t        L1Bit_91;
   Double_t        L1Bit_92;
   Double_t        L1Bit_93;
   Double_t        L1Bit_94;
   Double_t        L1Bit_95;
   Double_t        L1Bit_96;
   Double_t        L1Bit_97;
   Double_t        L1Bit_98;
   Double_t        L1Bit_99;
   Double_t        L1Bit_TechTrig_0;
   Double_t        L1Bit_TechTrig_1;
   Double_t        L1Bit_TechTrig_10;
   Double_t        L1Bit_TechTrig_11;
   Double_t        L1Bit_TechTrig_12;
   Double_t        L1Bit_TechTrig_13;
   Double_t        L1Bit_TechTrig_14;
   Double_t        L1Bit_TechTrig_15;
   Double_t        L1Bit_TechTrig_16;
   Double_t        L1Bit_TechTrig_17;
   Double_t        L1Bit_TechTrig_18;
   Double_t        L1Bit_TechTrig_19;
   Double_t        L1Bit_TechTrig_2;
   Double_t        L1Bit_TechTrig_20;
   Double_t        L1Bit_TechTrig_21;
   Double_t        L1Bit_TechTrig_22;
   Double_t        L1Bit_TechTrig_23;
   Double_t        L1Bit_TechTrig_24;
   Double_t        L1Bit_TechTrig_25;
   Double_t        L1Bit_TechTrig_26;
   Double_t        L1Bit_TechTrig_27;
   Double_t        L1Bit_TechTrig_28;
   Double_t        L1Bit_TechTrig_29;
   Double_t        L1Bit_TechTrig_3;
   Double_t        L1Bit_TechTrig_30;
   Double_t        L1Bit_TechTrig_31;
   Double_t        L1Bit_TechTrig_32;
   Double_t        L1Bit_TechTrig_33;
   Double_t        L1Bit_TechTrig_34;
   Double_t        L1Bit_TechTrig_35;
   Double_t        L1Bit_TechTrig_36;
   Double_t        L1Bit_TechTrig_37;
   Double_t        L1Bit_TechTrig_38;
   Double_t        L1Bit_TechTrig_39;
   Double_t        L1Bit_TechTrig_4;
   Double_t        L1Bit_TechTrig_40;
   Double_t        L1Bit_TechTrig_41;
   Double_t        L1Bit_TechTrig_42;
   Double_t        L1Bit_TechTrig_43;
   Double_t        L1Bit_TechTrig_44;
   Double_t        L1Bit_TechTrig_45;
   Double_t        L1Bit_TechTrig_46;
   Double_t        L1Bit_TechTrig_47;
   Double_t        L1Bit_TechTrig_48;
   Double_t        L1Bit_TechTrig_49;
   Double_t        L1Bit_TechTrig_5;
   Double_t        L1Bit_TechTrig_50;
   Double_t        L1Bit_TechTrig_51;
   Double_t        L1Bit_TechTrig_52;
   Double_t        L1Bit_TechTrig_53;
   Double_t        L1Bit_TechTrig_54;
   Double_t        L1Bit_TechTrig_55;
   Double_t        L1Bit_TechTrig_56;
   Double_t        L1Bit_TechTrig_57;
   Double_t        L1Bit_TechTrig_58;
   Double_t        L1Bit_TechTrig_59;
   Double_t        L1Bit_TechTrig_6;
   Double_t        L1Bit_TechTrig_60;
   Double_t        L1Bit_TechTrig_61;
   Double_t        L1Bit_TechTrig_62;
   Double_t        L1Bit_TechTrig_63;
   Double_t        L1Bit_TechTrig_7;
   Double_t        L1Bit_TechTrig_8;
   Double_t        L1Bit_TechTrig_9;


   // List of branches
   TBranch        *b_AlCa_EcalEta_8E29;   //!
   TBranch        *b_AlCa_EcalPhiSym;   //!
   TBranch        *b_AlCa_EcalPi0_8E29;   //!
   TBranch        *b_AlCa_HcalPhiSym;   //!
   TBranch        *b_AlCa_RPCMuonNoHits;   //!
   TBranch        *b_AlCa_RPCMuonNormalisation;   //!
   TBranch        *b_HLTAnalyzerEndpath;   //!
   TBranch        *b_HLT_BTagIP_Jet50U;   //!
   TBranch        *b_HLT_BTagMu_Jet10U;   //!
   TBranch        *b_HLT_BackwardBSC;   //!
   TBranch        *b_HLT_CSCBeamHalo;   //!
   TBranch        *b_HLT_CSCBeamHaloOverlapRing1;   //!
   TBranch        *b_HLT_CSCBeamHaloOverlapRing2;   //!
   TBranch        *b_HLT_CSCBeamHaloRing2or3;   //!
   TBranch        *b_HLT_DiJetAve15U_8E29;   //!
   TBranch        *b_HLT_DiJetAve30U_8E29;   //!
   TBranch        *b_HLT_DiJetAve50U;   //!
   TBranch        *b_HLT_DoubleEle5_SW_L1R;   //!
   TBranch        *b_HLT_DoubleLooseIsoTau15;   //!
   TBranch        *b_HLT_DoubleMu0;   //!
   TBranch        *b_HLT_DoubleMu3;   //!
   TBranch        *b_HLT_DoublePhoton10_L1R;   //!
   TBranch        *b_HLT_DoublePhoton5_Jpsi_L1R;   //!
   TBranch        *b_HLT_DoublePhoton5_Upsilon_L1R;   //!
   TBranch        *b_HLT_DoublePhoton5_eeRes_L1R;   //!
   TBranch        *b_HLT_Ele10_LW_EleId_L1R;   //!
   TBranch        *b_HLT_Ele10_LW_L1R;   //!
   TBranch        *b_HLT_Ele15_LW_L1R;   //!
   TBranch        *b_HLT_Ele15_SC10_LW_L1R;   //!
   TBranch        *b_HLT_Ele15_SiStrip_L1R;   //!
   TBranch        *b_HLT_Ele20_LW_L1R;   //!
   TBranch        *b_HLT_ForwardBSC;   //!
   TBranch        *b_HLT_FwdJet20U;   //!
   TBranch        *b_HLT_HT100U;   //!
   TBranch        *b_HLT_IsoMu3;   //!
   TBranch        *b_HLT_IsoTrack_8E29;   //!
   TBranch        *b_HLT_Jet100U;   //!
   TBranch        *b_HLT_Jet15U;   //!
   TBranch        *b_HLT_Jet30U;   //!
   TBranch        *b_HLT_Jet50U;   //!
   TBranch        *b_HLT_Jet70U;   //!
   TBranch        *b_HLT_L1DoubleEG5;   //!
   TBranch        *b_HLT_L1DoubleMuOpen;   //!
   TBranch        *b_HLT_L1Jet6U;   //!
   TBranch        *b_HLT_L1MET20;   //!
   TBranch        *b_HLT_L1Mu;   //!
   TBranch        *b_HLT_L1Mu14_L1ETM30;   //!
   TBranch        *b_HLT_L1Mu14_L1SingleEG10;   //!
   TBranch        *b_HLT_L1Mu14_L1SingleJet6U;   //!
   TBranch        *b_HLT_L1Mu20;   //!
   TBranch        *b_HLT_L1MuOpen;   //!
   TBranch        *b_HLT_L1SingleEG5;   //!
   TBranch        *b_HLT_L1SingleEG8;   //!
   TBranch        *b_HLT_L2Mu11;   //!
   TBranch        *b_HLT_L2Mu9;   //!
   TBranch        *b_HLT_MET100;   //!
   TBranch        *b_HLT_MET45;   //!
   TBranch        *b_HLT_MinBiasEcal;   //!
   TBranch        *b_HLT_MinBiasHcal;   //!
   TBranch        *b_HLT_MinBiasPixel;   //!
   TBranch        *b_HLT_MinBiasPixel_Trk5;   //!
   TBranch        *b_HLT_Mu3;   //!
   TBranch        *b_HLT_Mu5;   //!
   TBranch        *b_HLT_Mu9;   //!
   TBranch        *b_HLT_Photon10_L1R;   //!
   TBranch        *b_HLT_Photon15_L1R;   //!
   TBranch        *b_HLT_Photon15_LooseEcalIso_L1R;   //!
   TBranch        *b_HLT_Photon15_TrackIso_L1R;   //!
   TBranch        *b_HLT_Photon20_L1R;   //!
   TBranch        *b_HLT_Photon30_L1R_8E29;   //!
   TBranch        *b_HLT_QuadJet15U;   //!
   TBranch        *b_HLT_SingleLooseIsoTau20;   //!
   TBranch        *b_HLT_StoppedHSCP_8E29;   //!
   TBranch        *b_HLT_TrackerCosmics;   //!
   TBranch        *b_HLT_ZeroBias;   //!
   TBranch        *b_HLTriggerFinalPath;   //!
   TBranch        *b_HLTriggerFirstPath;   //!
   TBranch        *b_L1Bit_0;   //!
   TBranch        *b_L1Bit_1;   //!
   TBranch        *b_L1Bit_10;   //!
   TBranch        *b_L1Bit_100;   //!
   TBranch        *b_L1Bit_101;   //!
   TBranch        *b_L1Bit_102;   //!
   TBranch        *b_L1Bit_103;   //!
   TBranch        *b_L1Bit_104;   //!
   TBranch        *b_L1Bit_105;   //!
   TBranch        *b_L1Bit_106;   //!
   TBranch        *b_L1Bit_107;   //!
   TBranch        *b_L1Bit_108;   //!
   TBranch        *b_L1Bit_109;   //!
   TBranch        *b_L1Bit_11;   //!
   TBranch        *b_L1Bit_110;   //!
   TBranch        *b_L1Bit_111;   //!
   TBranch        *b_L1Bit_112;   //!
   TBranch        *b_L1Bit_113;   //!
   TBranch        *b_L1Bit_114;   //!
   TBranch        *b_L1Bit_115;   //!
   TBranch        *b_L1Bit_116;   //!
   TBranch        *b_L1Bit_117;   //!
   TBranch        *b_L1Bit_118;   //!
   TBranch        *b_L1Bit_119;   //!
   TBranch        *b_L1Bit_12;   //!
   TBranch        *b_L1Bit_120;   //!
   TBranch        *b_L1Bit_121;   //!
   TBranch        *b_L1Bit_122;   //!
   TBranch        *b_L1Bit_123;   //!
   TBranch        *b_L1Bit_124;   //!
   TBranch        *b_L1Bit_125;   //!
   TBranch        *b_L1Bit_126;   //!
   TBranch        *b_L1Bit_127;   //!
   TBranch        *b_L1Bit_13;   //!
   TBranch        *b_L1Bit_14;   //!
   TBranch        *b_L1Bit_15;   //!
   TBranch        *b_L1Bit_16;   //!
   TBranch        *b_L1Bit_17;   //!
   TBranch        *b_L1Bit_18;   //!
   TBranch        *b_L1Bit_19;   //!
   TBranch        *b_L1Bit_2;   //!
   TBranch        *b_L1Bit_20;   //!
   TBranch        *b_L1Bit_21;   //!
   TBranch        *b_L1Bit_22;   //!
   TBranch        *b_L1Bit_23;   //!
   TBranch        *b_L1Bit_24;   //!
   TBranch        *b_L1Bit_25;   //!
   TBranch        *b_L1Bit_26;   //!
   TBranch        *b_L1Bit_27;   //!
   TBranch        *b_L1Bit_28;   //!
   TBranch        *b_L1Bit_29;   //!
   TBranch        *b_L1Bit_3;   //!
   TBranch        *b_L1Bit_30;   //!
   TBranch        *b_L1Bit_31;   //!
   TBranch        *b_L1Bit_32;   //!
   TBranch        *b_L1Bit_33;   //!
   TBranch        *b_L1Bit_34;   //!
   TBranch        *b_L1Bit_35;   //!
   TBranch        *b_L1Bit_36;   //!
   TBranch        *b_L1Bit_37;   //!
   TBranch        *b_L1Bit_38;   //!
   TBranch        *b_L1Bit_39;   //!
   TBranch        *b_L1Bit_4;   //!
   TBranch        *b_L1Bit_40;   //!
   TBranch        *b_L1Bit_41;   //!
   TBranch        *b_L1Bit_42;   //!
   TBranch        *b_L1Bit_43;   //!
   TBranch        *b_L1Bit_44;   //!
   TBranch        *b_L1Bit_45;   //!
   TBranch        *b_L1Bit_46;   //!
   TBranch        *b_L1Bit_47;   //!
   TBranch        *b_L1Bit_48;   //!
   TBranch        *b_L1Bit_49;   //!
   TBranch        *b_L1Bit_5;   //!
   TBranch        *b_L1Bit_50;   //!
   TBranch        *b_L1Bit_51;   //!
   TBranch        *b_L1Bit_52;   //!
   TBranch        *b_L1Bit_53;   //!
   TBranch        *b_L1Bit_54;   //!
   TBranch        *b_L1Bit_55;   //!
   TBranch        *b_L1Bit_56;   //!
   TBranch        *b_L1Bit_57;   //!
   TBranch        *b_L1Bit_58;   //!
   TBranch        *b_L1Bit_59;   //!
   TBranch        *b_L1Bit_6;   //!
   TBranch        *b_L1Bit_60;   //!
   TBranch        *b_L1Bit_61;   //!
   TBranch        *b_L1Bit_62;   //!
   TBranch        *b_L1Bit_63;   //!
   TBranch        *b_L1Bit_64;   //!
   TBranch        *b_L1Bit_65;   //!
   TBranch        *b_L1Bit_66;   //!
   TBranch        *b_L1Bit_67;   //!
   TBranch        *b_L1Bit_68;   //!
   TBranch        *b_L1Bit_69;   //!
   TBranch        *b_L1Bit_7;   //!
   TBranch        *b_L1Bit_70;   //!
   TBranch        *b_L1Bit_71;   //!
   TBranch        *b_L1Bit_72;   //!
   TBranch        *b_L1Bit_73;   //!
   TBranch        *b_L1Bit_74;   //!
   TBranch        *b_L1Bit_75;   //!
   TBranch        *b_L1Bit_76;   //!
   TBranch        *b_L1Bit_77;   //!
   TBranch        *b_L1Bit_78;   //!
   TBranch        *b_L1Bit_79;   //!
   TBranch        *b_L1Bit_8;   //!
   TBranch        *b_L1Bit_80;   //!
   TBranch        *b_L1Bit_81;   //!
   TBranch        *b_L1Bit_82;   //!
   TBranch        *b_L1Bit_83;   //!
   TBranch        *b_L1Bit_84;   //!
   TBranch        *b_L1Bit_85;   //!
   TBranch        *b_L1Bit_86;   //!
   TBranch        *b_L1Bit_87;   //!
   TBranch        *b_L1Bit_88;   //!
   TBranch        *b_L1Bit_89;   //!
   TBranch        *b_L1Bit_9;   //!
   TBranch        *b_L1Bit_90;   //!
   TBranch        *b_L1Bit_91;   //!
   TBranch        *b_L1Bit_92;   //!
   TBranch        *b_L1Bit_93;   //!
   TBranch        *b_L1Bit_94;   //!
   TBranch        *b_L1Bit_95;   //!
   TBranch        *b_L1Bit_96;   //!
   TBranch        *b_L1Bit_97;   //!
   TBranch        *b_L1Bit_98;   //!
   TBranch        *b_L1Bit_99;   //!
   TBranch        *b_L1Bit_TechTrig_0;   //!
   TBranch        *b_L1Bit_TechTrig_1;   //!
   TBranch        *b_L1Bit_TechTrig_10;   //!
   TBranch        *b_L1Bit_TechTrig_11;   //!
   TBranch        *b_L1Bit_TechTrig_12;   //!
   TBranch        *b_L1Bit_TechTrig_13;   //!
   TBranch        *b_L1Bit_TechTrig_14;   //!
   TBranch        *b_L1Bit_TechTrig_15;   //!
   TBranch        *b_L1Bit_TechTrig_16;   //!
   TBranch        *b_L1Bit_TechTrig_17;   //!
   TBranch        *b_L1Bit_TechTrig_18;   //!
   TBranch        *b_L1Bit_TechTrig_19;   //!
   TBranch        *b_L1Bit_TechTrig_2;   //!
   TBranch        *b_L1Bit_TechTrig_20;   //!
   TBranch        *b_L1Bit_TechTrig_21;   //!
   TBranch        *b_L1Bit_TechTrig_22;   //!
   TBranch        *b_L1Bit_TechTrig_23;   //!
   TBranch        *b_L1Bit_TechTrig_24;   //!
   TBranch        *b_L1Bit_TechTrig_25;   //!
   TBranch        *b_L1Bit_TechTrig_26;   //!
   TBranch        *b_L1Bit_TechTrig_27;   //!
   TBranch        *b_L1Bit_TechTrig_28;   //!
   TBranch        *b_L1Bit_TechTrig_29;   //!
   TBranch        *b_L1Bit_TechTrig_3;   //!
   TBranch        *b_L1Bit_TechTrig_30;   //!
   TBranch        *b_L1Bit_TechTrig_31;   //!
   TBranch        *b_L1Bit_TechTrig_32;   //!
   TBranch        *b_L1Bit_TechTrig_33;   //!
   TBranch        *b_L1Bit_TechTrig_34;   //!
   TBranch        *b_L1Bit_TechTrig_35;   //!
   TBranch        *b_L1Bit_TechTrig_36;   //!
   TBranch        *b_L1Bit_TechTrig_37;   //!
   TBranch        *b_L1Bit_TechTrig_38;   //!
   TBranch        *b_L1Bit_TechTrig_39;   //!
   TBranch        *b_L1Bit_TechTrig_4;   //!
   TBranch        *b_L1Bit_TechTrig_40;   //!
   TBranch        *b_L1Bit_TechTrig_41;   //!
   TBranch        *b_L1Bit_TechTrig_42;   //!
   TBranch        *b_L1Bit_TechTrig_43;   //!
   TBranch        *b_L1Bit_TechTrig_44;   //!
   TBranch        *b_L1Bit_TechTrig_45;   //!
   TBranch        *b_L1Bit_TechTrig_46;   //!
   TBranch        *b_L1Bit_TechTrig_47;   //!
   TBranch        *b_L1Bit_TechTrig_48;   //!
   TBranch        *b_L1Bit_TechTrig_49;   //!
   TBranch        *b_L1Bit_TechTrig_5;   //!
   TBranch        *b_L1Bit_TechTrig_50;   //!
   TBranch        *b_L1Bit_TechTrig_51;   //!
   TBranch        *b_L1Bit_TechTrig_52;   //!
   TBranch        *b_L1Bit_TechTrig_53;   //!
   TBranch        *b_L1Bit_TechTrig_54;   //!
   TBranch        *b_L1Bit_TechTrig_55;   //!
   TBranch        *b_L1Bit_TechTrig_56;   //!
   TBranch        *b_L1Bit_TechTrig_57;   //!
   TBranch        *b_L1Bit_TechTrig_58;   //!
   TBranch        *b_L1Bit_TechTrig_59;   //!
   TBranch        *b_L1Bit_TechTrig_6;   //!
   TBranch        *b_L1Bit_TechTrig_60;   //!
   TBranch        *b_L1Bit_TechTrig_61;   //!
   TBranch        *b_L1Bit_TechTrig_62;   //!
   TBranch        *b_L1Bit_TechTrig_63;   //!
   TBranch        *b_L1Bit_TechTrig_7;   //!
   TBranch        *b_L1Bit_TechTrig_8;   //!
   TBranch        *b_L1Bit_TechTrig_9;   //!
   
   TChain* chainB;
   TChain* chainV;
   TChain* chainA;
   void InitializeA(TChain *fChain);
   void InitializeB(TChain *fChain);
   void InitializeV(TChain *fChain);

  // ==== END giant copy/paste of cfA variable ================================================================



};

#endif
