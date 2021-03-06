// -*- C++ -*-

#ifndef EVENTCALCULATOR_H
#define EVENTCALCULATOR_H

#include "CrossSectionTable.h"

#include "MiscUtil.cxx"

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

/*
left undone -- have not tried to port luke's jet killing stuff (should be done)
*/
#include "BasicLoopCU.h" //get all of the tree-related variables

class TRandom3;

//global constants
const double mW_ = 80.399;
const double mtop_ = 172.0;
const double lumi_ = 1.; //fix to 1/pb and scale MC later (e.g. in drawReducedTrees)
const double mZ_ = 91.188;

//define the BNN functions
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


class cellECAL {
public:
  cellECAL(double,double,int);
  double eta;
  double phi;
  int status;
};
cellECAL::cellECAL(double a, double b, int c) {
  eta=a;
  phi=b;
  status=c;
}
  

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

  EventCalculator(const TString & sampleName, jetType theJetType, METType theMETType);
  ~EventCalculator();

  //setters
  void setBTaggerType(BTaggerType btaggertype) {theBTaggerType_ = btaggertype;}
  void setOptions( const TString & opt);

  //loop over events
  void reducedTree(TString outputpath, itreestream& stream);

  void cutflow(itreestream& stream, int maxevents);
  void sampleAnalyzer(itreestream& stream);
  void plotBTagEffMC(itreestream& stream);

  //load external list of event ID's
  void loadEventList(std::vector<int> &vrun, std::vector<int> &vlumi, std::vector<int> &vevent);
  bool inEventList(std::vector<int> &vrun, std::vector<int> &vlumi, std::vector<int> &vevent);

  // functions that calculate stuff
  double getWeight(Long64_t nentries);
  float getPUWeight(reweight::LumiReWeighting lumiWeights);
  float getPUWeight(Lumi3DReWeighting lumiWeights);

  bool isGoodMuon(const unsigned int imuon, const bool disableRelIso=false, const float ptthreshold=10);
  bool isGoodRecoMuon(const unsigned int imuon, const bool disableRelIso=false, const float ptthreshold=10);
  bool isInAccRecoMuon(const unsigned int imuon, const float ptthreshold);
  bool isGoodElectron(const unsigned int iele, const bool disableRelIso=false, const float ptthreshold=10);
  bool isTightElectron(const unsigned int iele, const float isothreshold, const float ptthreshold=10);
  bool isInAccRecoElectron(const unsigned int iele, const float ptthreshold);
  bool isIDRecoElectron(const unsigned int iele, const float ptthreshold);
  unsigned int countEle(const float ptthreshold=10) ;
  bool isCleanMuon(const unsigned int imuon, const float ptthreshold=10);
  unsigned int countMu(const float ptthreshold=10);

  bool isZmumuCandidateEvent(const float ptthreshold, const float mllthreshold, float& m_ll);
  bool isZeeCandidateEvent(const float ptthreshold, const float mllthreshold, float& m_ll);


  bool isGoodTau(const unsigned int itau, const float pTthreshold=20, const float etaMax=2.4);
  float getTauPt( unsigned int itau );
  unsigned int countTau();

  //  bool passBEFilter(); //TODO migrate this code
  bool passHLT();
  bool passUtilityHLT(int &version, int &prescale);
  bool passUtilityPrescaleModuleHLT();
  bool passZmumuHLT();
  bool passZeeHLT();
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

  float getZllMET(bool isMuMu, bool isHybridMC, float& ZllMETphi);

  void getTransverseThrustVariables(float & thrust, float & thrustPhi, bool addMET);
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
 double getDeltaPhiMET(unsigned int n, float ptThreshold = 50, bool bjetsonly = false);
  
  double getTransverseMETError(unsigned int thisJet);
  double getDeltaPhiNMET(unsigned int thisJet); //Luke
  double getMinDeltaPhiNMET(unsigned int maxjets); //Luke
  
  double getDeltaPhiMETN_deltaT(unsigned int ijet, float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith, bool includeMismeasuredJet=false);
  double getDeltaPhiMETN_deltaT(unsigned int ijet) { return getDeltaPhiMETN_deltaT(ijet,30,2.4,true,false,false);  } //overloaded
  double getDeltaPhiMETN( unsigned int goodJetN, float mainpt, float maineta, bool mainid, float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith, bool includeMismeasuredJet=false, bool addPi=false); //Ben
  double getDeltaPhiMETN( unsigned int goodJetN ) {return getDeltaPhiMETN(goodJetN,50,2.4,true,30,2.4,true,false,false); }; //Ben, overloaded
  double getMinDeltaPhiMETN(unsigned int maxjets, float mainmt, float maineta, bool mainid, float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith, bool includeLeptons=false, bool includeMismeasuredJet=false, bool addPi=false ); //Ben
  double getMinDeltaPhiMETN(unsigned int maxjets) {return getMinDeltaPhiMETN(maxjets,50,2.4,true,30,2.4,true,false,false); }; //Ben, overloaded

  //mc truth mdpn
  double getDeltaT_MC(unsigned int ijet, bool addquad, bool keith);
  double getDeltaPhiMETN_MC(unsigned int goodJetN, bool addquad, bool keith);
  double getMinDeltaPhiMETN_MC(bool addquad, bool keith);

  void minDeltaPhiMETN_diySmear(TString smearingType, float &MET_g, float &MET_s, float &HT_s, float &mdpN_s, int &chosenJet, int &njets);


  //for Zll-modified MET 
  double getDeltaPhiZllMETN_deltaT(unsigned int ijet, float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith, bool ismumu, bool isHybridMC);
  double getDeltaPhiZllMETN_deltaT(unsigned int ijet, bool ismumu, bool isHybridMC) { return getDeltaPhiZllMETN_deltaT(ijet,30,2.4,true,false,false,ismumu, isHybridMC);  } 
  double getDeltaPhiZllMETN( unsigned int goodJetN, float mainpt, float maineta, bool mainid, float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith, bool ismumu, bool isHybridMC ); 
  double getDeltaPhiZllMETN( unsigned int goodJetN, bool ismumu, bool isHybridMC) {return getDeltaPhiZllMETN(goodJetN,50,2.4,true,30,2.4,true,false,false,ismumu, isHybridMC); }; 
  double getMinDeltaPhiZllMETN(unsigned int maxjets, float mainmt, float maineta, bool mainid, float otherpt, float othereta, bool otherid, bool dataJetRes, bool keith, bool ismumu, bool isHybridMC);
  double getMinDeltaPhiZllMETN(unsigned int maxjets, bool ismumu, bool isHybridMC = false) {return getMinDeltaPhiZllMETN(maxjets,50,2.4,true,30,2.4,true,false,false,ismumu, isHybridMC); }; 

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

  float getMaxJetMis(unsigned int rank, unsigned int maxjets, float jetpt, bool abs);
  float getMaxJetFracMis(unsigned int rank, unsigned int maxjets, float jetpt);
  float getDeltaPhiMETJetMaxMis(float jetpt);
  unsigned int getRankJetMaxMis(unsigned int maxjets, float jetpt);
  int getJetMisCategoryType1(bool skipMETfraction);
  int getJetMisCategoryType2(unsigned int targetJet, bool skipMETfraction);

  float getDeltaPhiStar(int & badjet);

  double getMinDeltaPhiMETTaus();
  double getMinDeltaPhiMETMuons(unsigned int maxmuons);

  void getCorrectedMET(float& correctedMET, float& correctedMETPhi);
  void getUncorrectedMET(float& uncorrectedMET, float& uncorrectedMETPhi);

  int doPBNR(); //particle based noise rejection

  //functions that are mostly used internally
  ULong64_t getRunNumber() {return TMath::Nint( *myEDM_run );}
  ULong64_t getLumiSection() {return TMath::Nint( *myEDM_luminosityBlock );}
  ULong64_t getEventNumber() {return TMath::Nint( *myEDM_event );}

  void hadronicTopFinder_DeltaR(float & mjjb1, float & mjjb2 , float & topPT1, float & topPT2);

  bool isCleanJet(const unsigned int ijet);
  bool isGoodJet(const unsigned int ijet, const float pTthreshold=50, const float etaMax=2.4, const bool jetid=true); //subset of isCleanJet
  bool isGoodJet10(unsigned int ijet) {return isGoodJet(ijet,10,2.4,true);}
  bool isGoodJet30(unsigned int ijet) {return isGoodJet(ijet,30,2.4,true);}
  bool isGoodJetMHT(unsigned int ijet);
  bool passBTagger(int ijet, BTaggerType btagger=Nbtaggers );
  bool passVLBTagger(int ijet);

  unsigned int nGoodJets();
  //  unsigned int nGoodJets(TH2D* count,TH2D* unc,TH2D* l2l3); //for a test
  unsigned int nGoodJets30();
  unsigned int nGoodBJets( BTaggerType btagger=Nbtaggers);
  unsigned int nTrueBJets();
  unsigned int nGoodVLBJets();

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
  float getJetTCHE( unsigned int ijet );

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

  float bjetPtOfN(unsigned int n, BTaggerType thebtaggertype=Nbtaggers);
  float bjetPhiOfN(unsigned int n, BTaggerType thebtaggertype=Nbtaggers);
  float bjetEtaOfN(unsigned int n, BTaggerType thebtaggertype=Nbtaggers);
  float bjetEnergyOfN(unsigned int n, BTaggerType thebtaggertype=Nbtaggers);
  float bjetCSVOfN(unsigned int n, BTaggerType thebtaggertype=Nbtaggers);
  int bjetFlavorOfN(unsigned int n, BTaggerType thebtaggertype=Nbtaggers);
  float bjetChargedHadronFracOfN(unsigned int n, BTaggerType thebtaggertype=Nbtaggers);
  int bjetChargedHadronMultOfN(unsigned int n, BTaggerType thebtaggertype=Nbtaggers);

  float vlBjetCSVOfN(unsigned int n);
  float vlBjetTCHEOfN(unsigned int n);

  float elePtOfN(unsigned int n, const float ptthreshold=10);
  float eleEtaOfN(unsigned int n, const float ptthreshold=10);
  float elePhiOfN(unsigned int n, const float ptthreshold=10);
  int eleChargeOfN(unsigned int n, const float ptthreshold=10);
  float eleHOverEOfN(unsigned int n, const float ptthreshold=10);
  float eleDphiOfN(unsigned int n, const float ptthreshold=10);
  float eleDetaOfN(unsigned int n, const float ptthreshold=10);
  float eleSigmaIetaIetaOfN(unsigned int n, const float ptthreshold=10);

  float muonPtOfN(unsigned int n, const float ptthreshold=10);
  float muonEtaOfN(unsigned int n, const float ptthreshold=10);
  float muonPhiOfN(unsigned int n, const float ptthreshold=10);
  float muonIsoOfN(unsigned int n, const float ptthreshold=10);
  float muonChHadIsoOfN(unsigned int n, const float ptthreshold=10);
  float muonPhotonIsoOfN(unsigned int n, const float ptthreshold=10);
  float muonNeutralHadIsoOfN(unsigned int n, const float ptthreshold=10);
  int muonChargeOfN(unsigned int n, const float ptthreshold=10);

  float tauPtOfN(unsigned int n);
  float tauEtaOfN(unsigned int n);

  float recoMuonPtOfN(unsigned int n, const float ptthreshold=10);
  float recoMuonEtaOfN(unsigned int n, const float ptthreshold=10);
  float recoMuonPhiOfN(unsigned int n, const float ptthreshold=10);
  float recoMuonIsoOfN(unsigned int n, const float ptthreshold=10);
  float recoMuonMinDeltaPhiJetOfN(unsigned int n, const float ptthreshold=10);

  float getRelIsoForIsolationStudyEle();
  float getRelIsoForIsolationStudyMuon();

  std::vector<bool> passTagAndProbeMuon(std::vector<double>& probes_mll);
  std::vector<bool> passTagAndProbeElectron(std::vector<double>& probes_mll, TString tagMode="nominal", bool relIDmode=false);

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
  int muonMatch(const int trueMuon, const float ptthreshold=10, bool mustMatchToGood=false, bool matchtoRECO=false);
  int electronMatch(const int trueElectron, const float ptthreshold=10, bool mustMatchToGood=false, bool matchtoRECO=false);
  int tauMatch(const int trueTau);
  int daughterMatch(const int Wdaughter, const int WdecayType);
  int getTTbarDecayType(int& W1decayType, int& W2decayType, int& W1, int& W1daughter, int& W2, int& W2daughter, bool passW2info);
  int getWDecayType(int& WdecayType, int& W, int& Wdaughter, bool fromtop);

  bool findZ(int& decaymode, int& lepton1_index, int& lepton2_index);
  bool genZllInAcceptance(int lepton1_index, int lepton2_index, float& zllMass);

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

  //FASTSIM b-tagging efficiencies
  float bJetFastsimSF(const TString & what, int flavor,float pt);
  float get_AN_12_175_Table2_Value(const float pt);
  float get_AN_12_175_Table3_Value(const float pt);
  float get_AN_12_175_Table4_Value(const float pt);
  float get_AN_12_175_Table5_Value(const float pt);
  float get_AN_12_175_Table6_Value(const float pt);
  float get_AN_12_175_Table8_Value(const float pt);
  float get_AN_12_175_Table2_Error(const float pt);
  int getPtBinIndex(const float pt) ;


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
  std::vector<jet2_s> * myJetsPF;
  std::vector<jethelper2_s> * myJetsPFhelper;

  std::vector<electron1_s> * myElectronsPF;
  std::vector<electron_s> * myElectronsRECO;
  std::vector<electronhelper1_s> * myElectronsPFhelper;
  std::vector<electronhelper_s> * myElectronsRECOhelper;
  std::vector<muon1_s> * myMuonsPF;
  std::vector<muon_s> * myMuonsRECO;
  std::vector<muonhelper1_s> * myMuonsPFhelper;
  std::vector<muonhelper_s> * myMuonsRECOhelper;
  std::vector<tau_s> * myTausPF;
  std::vector<met1_s> * myMETPF;
  std::vector<met4_s> * myMETPFType1;
  std::vector<met_s> * myMETcalo;
  std::vector<vertex_s> * myVertex;
  std::vector<genparticlehelperra2_s> * myGenParticles;
  double* myGenWeight;

  //versions for jet-loss study
  std::vector<jet2_s> * myJetsPF_temp;
  std::vector<met1_s> * myMETPF_temp;

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

  void loadDiySmear();
  TF1 *hDiyGaus_;
  TF1 *hDiy3Gaus_;

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



  void checkConsistency();

  void  loadSusyScanCrossSections();

  unsigned int getSeed();
  void changeVariables(TRandom3* random, double jetLossProbability, int& nLostJets);
  void changeVariablesGenSmearing(TRandom3* random);
  void resetVariables();
  bool  recalculatedVariables_;

  TStopwatch* watch_; //not used for everyday running; just for testing

  std::vector<TString> cutTags_;
  std::map<TString, TString> cutNames_; //key is a cutTag. these should be "human readable" but should not have any spaces
  std::vector<TString> ignoredCut_; //allow more than 1 ignored cut!
  std::vector<TString> requiredCut_; //new feature to *turn on* a cut that is usually not required by a given cut scheme


  int ZmumuCand1_, ZmumuCand2_;
  int ZeeCand1_, ZeeCand2_;
  int ZllMCCand1_, ZllMCCand2_;
  
  struct sort_pairF {
    bool operator()(const std::pair<double,double> &left, const std::pair<double,double> &right) {
      return left.first < right.first;
    }
  };
  
};

#endif
