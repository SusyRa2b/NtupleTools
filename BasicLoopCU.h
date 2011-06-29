#ifndef BASICLOOPCU_H
#define BASICLOOPCU_H
//-----------------------------------------------------------------------------
// File:        ra2banalyzer.h
// Description: Analyzer header for ntuples created by TheNtupleMaker
// Created:     Wed Jun 29 04:27:16 2011 by mkntanalyzer.py
// Author:      Wee Teo
//-----------------------------------------------------------------------------

// -- System

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

#ifdef PROJECT_NAME

// --- CMSSW

#include "PhysicsTools/TheNtupleMaker/interface/treestream.h"
#include "PhysicsTools/TheNtupleMaker/interface/pdg.h"

#else

#include "treestream.h"
#include "pdg.h"

#endif

// -- Root

#include "TROOT.h"
#include "TApplication.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1F.h"
#include "TH2F.h"
//-----------------------------------------------------------------------------
// -- Utilities
//-----------------------------------------------------------------------------
void
error(std::string message)
{
  std::cout << "** error ** " << message << std::endl;
  exit(0);
}

std::string 
strip(std::string line)
{
  int l = line.size();
  if ( l == 0 ) return std::string("");
  int n = 0;
  while (((line[n] == 0)    ||
	  (line[n] == ' ' ) ||
	  (line[n] == '\n') ||
	  (line[n] == '\t')) && n < l) n++;

  int m = l-1;
  while (((line[m] == 0)    ||
	  (line[m] == ' ')  ||
	  (line[m] == '\n') ||
	  (line[m] == '\t')) && m > 0) m--;
  return line.substr(n,m-n+1);
}

std::string
nameonly(std::string filename)
{
  int i = filename.rfind("/");
  int j = filename.rfind(".");
  if ( j < 0 ) j = filename.size();
  return filename.substr(i+1,j-i-1);
}
//-----------------------------------------------------------------------------
struct outputFile
{
  outputFile(std::string filename)
   : filename_(filename),
	 file_(new TFile(filename_.c_str(), "recreate")),
	 tree_(0),
	 b_weight_(0),
	 entry_(0),
	 SAVECOUNT_(50000)
  {
	file_->cd();
	hist_ = new TH1F("counts", "", 1,0,1);
	hist_->SetBit(TH1::kCanRebin);
	hist_->SetStats(0);
  }

  outputFile(std::string filename, itreestream& stream, int savecount=50000) 
   : filename_(filename),
	 file_(new TFile(filename.c_str(), "recreate")),
	 tree_(stream.tree()->CloneTree(0)),
	 b_weight_(tree_->Branch("eventWeight", &weight_, "eventWeight/D")),
	 entry_(0),
	 SAVECOUNT_(savecount)
  {
	std::cout << "events will be skimmed to file "
			  << filename_ << std::endl;
	file_->cd();
	hist_ = new TH1F("counts", "", 1,0,1);
	hist_->SetBit(TH1::kCanRebin);
	hist_->SetStats(0);
  }

  void addEvent(double weight=1)
  {
    if ( tree_ == 0 ) return;
	
    weight_ = weight;
	file_   = tree_->GetCurrentFile();
	file_->cd();
	tree_->Fill();

	entry_++;
	if ( entry_ % SAVECOUNT_ == 0 )
	  tree_->AutoSave("SaveSelf");
  }

  void count(std::string cond, int w=1)
  {
    hist_->Fill(cond.c_str(), w);
  }
  
  void close()
  {
  	std::cout << "==> histograms saved to file " << filename_ << std::endl;
    if ( tree_ != 0 )
	  {
	    std::cout << "==> events skimmed to file " << filename_ << std::endl;
	    file_ = tree_->GetCurrentFile();
	  }
	file_->cd();
	//file_->Write("", TObject::kWriteDelete);
	file_->Write();
	file_->ls();
	file_->Close();
  }

  std::string filename_;  
  TFile* file_;
  TTree* tree_;
  TH1F*  hist_;
  TBranch* b_weight_;
  double     weight_;
  int    entry_;
  int    SAVECOUNT_;
};

struct commandLine
{
  std::string progname;
  std::string filelist;
  std::string outputfilename;
};


void
decodeCommandLine(int argc, char** argv, commandLine& cl)
{
  cl.progname = std::string(argv[0]);

  // 1st (optional) argument
  if ( argc > 1 )
	cl.filelist = std::string(argv[1]);
  else
	cl.filelist = std::string("filelist.txt");

  // 2nd (optional) command line argument
  if ( argc > 2 ) 
	cl.outputfilename = std::string(argv[2]);
  else
	cl.outputfilename = cl.progname + std::string("_histograms");

  // Make sure extension is ".root"
  std::string name = cl.outputfilename;
  if ( name.substr(name.size()-5, 5) != std::string(".root") )
    cl.outputfilename += std::string(".root");
}

// Read ntuple filenames from file list

std::vector<std::string>
getFilenames(std::string filelist)
{
  std::ifstream stream(filelist.c_str());
  if ( !stream.good() ) error("unable to open file: " + filelist);

  // Get list of ntuple files to be processed

  std::vector<std::string> v;
  std::string filename;
  while ( stream >> filename )
	if ( strip(filename) != "" ) v.push_back(filename);
  return v;
}
//-----------------------------------------------------------------------------
// -- Declare variables to be read
//-----------------------------------------------------------------------------
//double	basictype<double>2_value;
//double	basictype<double>_value;
//int	basictype<unsigned int>_value;
double	sigma;
double	rho;
int	flavorhistorypath;
double	beamspot_x0;
double	beamspot_y0;
double	beamspot_z0;
std::vector<float>	electron_caloIso(50,0);
std::vector<int>	electron_charge(50,0);
std::vector<float>	electron_chargedHadronIso(50,0);
std::vector<double>	electron_dB(50,0);
std::vector<double>	electron_deltaEtaSuperClusterTrackAtVtx(50,0);
std::vector<double>	electron_deltaPhiSuperClusterTrackAtVtx(50,0);
std::vector<float>	electron_dr03EcalRecHitSumEt(50,0);
std::vector<float>	electron_dr03HcalTowerSumEt(50,0);
std::vector<float>	electron_dr03TkSumPt(50,0);
std::vector<float>	electron_ecalIso(50,0);
std::vector<float>	electron_eidRobustTight(50,0);
std::vector<double>	electron_energy(50,0);
std::vector<double>	electron_et(50,0);
std::vector<double>	electron_eta(50,0);
std::vector<double>	electron_genLepton_eta(50,0);
std::vector<int>	electron_genLepton_pdgId(50,0);
std::vector<double>	electron_genLepton_phi(50,0);
std::vector<double>	electron_genLepton_pt(50,0);
std::vector<double>	electron_gsfTrack_d0(50,0);
std::vector<int>	electron_gsfTrack_trackerExpectedHitsInner_numberOfLostHits(50,0);
std::vector<double>	electron_hadronicOverEm(50,0);
std::vector<float>	electron_hcalIso(50,0);
std::vector<float>	electron_neutralHadronIso(50,0);
std::vector<double>	electron_phi(50,0);
std::vector<float>	electron_photonIso(50,0);
std::vector<double>	electron_pt(50,0);
std::vector<double>	electron_px(50,0);
std::vector<double>	electron_py(50,0);
std::vector<double>	electron_pz(50,0);
std::vector<double>	electron_sigmaIetaIeta(50,0);
std::vector<float>	electron_simpleEleId95relIso(50,0);
std::vector<double>	electron_superCluster_eta(50,0);
std::vector<float>	electron_trackIso(50,0);
std::vector<double>	electron_vz(50,0);
std::vector<double>	electronhelper_dxywrtBeamSpot(50,0);
int	eventhelper_bunchCrossing;
int	eventhelper_event;
int	eventhelper_isRealData;
int	eventhelper_luminosityBlock;
int	eventhelper_run;
double	geneventinfoproduct_pdf1;
double	geneventinfoproduct_pdf2;
double	geneventinfoproduct_scalePDF;
double	geneventinfoproduct_weight;
double	geneventinfoproduct_x1;
double	geneventinfoproduct_x2;
std::vector<double>	geneventinfoproducthelper_pdf1(50,0);
std::vector<double>	geneventinfoproducthelper_pdf2(50,0);
std::vector<double>	geneventinfoproducthelper_pdfweight(50,0);
std::vector<double>	geneventinfoproducthelper_pdfweightsum(50,0);
std::vector<int>	genparticlehelperra2b_charge(200,0);
std::vector<double>	genparticlehelperra2b_eta(200,0);
std::vector<int>	genparticlehelperra2b_firstDaughter(200,0);
std::vector<int>	genparticlehelperra2b_firstMother(200,0);
std::vector<int>	genparticlehelperra2b_lastDaughter(200,0);
std::vector<int>	genparticlehelperra2b_lastMother(200,0);
std::vector<double>	genparticlehelperra2b_mass(200,0);
std::vector<int>	genparticlehelperra2b_pdgId(200,0);
std::vector<double>	genparticlehelperra2b_phi(200,0);
std::vector<double>	genparticlehelperra2b_pt(200,0);
std::vector<int>	genparticlehelperra2b_status(200,0);
double	genruninfoproduct_externalXSecLO_error;
double	genruninfoproduct_externalXSecLO_value;
double	genruninfoproduct_filterEfficiency;
double	genruninfoproduct_internalXSec_error;
double	genruninfoproduct_internalXSec_value;
std::vector<float>	jet2_combinedSecondaryVertexBJetTags(500,0);
std::vector<float>	jet2_combinedSecondaryVertexMVABJetTags(500,0);
std::vector<double>	jet2_energy(500,0);
std::vector<double>	jet2_eta(500,0);
std::vector<float>	jet2_jetBProbabilityBJetTags(500,0);
std::vector<double>	jet2_jetID_fHPD(500,0);
std::vector<int>	jet2_jetID_n90Hits(500,0);
std::vector<float>	jet2_jetProbabilityBJetTags(500,0);
std::vector<int>	jet2_partonFlavour(500,0);
std::vector<double>	jet2_phi(500,0);
std::vector<double>	jet2_pt(500,0);
std::vector<float>	jet2_simpleSecondaryVertexBJetTags(500,0);
std::vector<float>	jet2_simpleSecondaryVertexHighEffBJetTags(500,0);
std::vector<float>	jet2_simpleSecondaryVertexHighPurBJetTags(500,0);
std::vector<float>	jet2_softElectronByIP3dBJetTags(500,0);
std::vector<float>	jet2_softElectronByPtBJetTags(500,0);
std::vector<float>	jet2_softMuonBJetTags(500,0);
std::vector<float>	jet2_softMuonByIP3dBJetTags(500,0);
std::vector<float>	jet2_softMuonByPtBJetTags(500,0);
std::vector<float>	jet2_trackCountingHighEffBJetTags(500,0);
std::vector<float>	jet2_trackCountingHighPurBJetTags(500,0);
std::vector<double>	jet2_uncor_energy(500,0);
std::vector<double>	jet2_uncor_eta(500,0);
std::vector<double>	jet2_uncor_phi(500,0);
std::vector<double>	jet2_uncor_pt(500,0);
std::vector<float>	jet_HFHadronEnergy(500,0);
std::vector<float>	jet_chargedEmEnergyFraction(500,0);
std::vector<float>	jet_chargedHadronEnergyFraction(500,0);
std::vector<int>	jet_chargedMultiplicity(500,0);
std::vector<float>	jet_combinedSecondaryVertexBJetTags(500,0);
std::vector<float>	jet_combinedSecondaryVertexMVABJetTags(500,0);
std::vector<float>	jet_emEnergyFraction(500,0);
std::vector<double>	jet_energy(500,0);
std::vector<double>	jet_eta(500,0);
std::vector<float>	jet_jetArea(500,0);
std::vector<float>	jet_jetBProbabilityBJetTags(500,0);
std::vector<double>	jet_jetID_fHPD(500,0);
std::vector<int>	jet_jetID_n90Hits(500,0);
std::vector<float>	jet_jetProbabilityBJetTags(500,0);
std::vector<float>	jet_neutralEmEnergyFraction(500,0);
std::vector<float>	jet_neutralHadronEnergy(500,0);
std::vector<float>	jet_neutralHadronEnergyFraction(500,0);
std::vector<int>	jet_numberOfDaughters(500,0);
std::vector<int>	jet_partonFlavour(500,0);
std::vector<double>	jet_phi(500,0);
std::vector<double>	jet_pt(500,0);
std::vector<float>	jet_simpleSecondaryVertexBJetTags(500,0);
std::vector<float>	jet_simpleSecondaryVertexHighEffBJetTags(500,0);
std::vector<float>	jet_simpleSecondaryVertexHighPurBJetTags(500,0);
std::vector<float>	jet_softElectronByIP3dBJetTags(500,0);
std::vector<float>	jet_softElectronByPtBJetTags(500,0);
std::vector<float>	jet_softMuonBJetTags(500,0);
std::vector<float>	jet_softMuonByIP3dBJetTags(500,0);
std::vector<float>	jet_softMuonByPtBJetTags(500,0);
std::vector<float>	jet_trackCountingHighEffBJetTags(500,0);
std::vector<float>	jet_trackCountingHighPurBJetTags(500,0);
std::vector<double>	jet_uncor_energy(500,0);
std::vector<double>	jet_uncor_eta(500,0);
std::vector<double>	jet_uncor_phi(500,0);
std::vector<double>	jet_uncor_pt(500,0);
std::vector<double>	jethelper_jecFactor(500,0);
std::vector<double>	jethelper_jecFactorNoL1Fast(500,0);
std::vector<double>	jethelper_jetUncMinus(500,0);
std::vector<double>	jethelper_jetUncPlus(500,0);
std::vector<double>	met2_energy(10,0);
std::vector<double>	met2_et(10,0);
std::vector<double>	met2_mEtSig(10,0);
std::vector<double>	met2_phi(10,0);
std::vector<double>	met2_pt(10,0);
std::vector<double>	met2_sumEt(10,0);
std::vector<double>	met3_energy(10,0);
std::vector<double>	met3_et(10,0);
std::vector<double>	met3_mEtSig(10,0);
std::vector<double>	met3_phi(10,0);
std::vector<double>	met3_pt(10,0);
std::vector<double>	met3_sumEt(10,0);
std::vector<double>	met_energy(10,0);
std::vector<double>	met_et(10,0);
std::vector<double>	met_mEtSig(10,0);
std::vector<double>	met_phi(10,0);
std::vector<double>	met_pt(10,0);
std::vector<double>	met_sumEt(10,0);
std::vector<int>	muon_AllGlobalMuons(50,0);
std::vector<int>	muon_GlobalMuonPromptTight(50,0);
std::vector<int>	muon_charge(50,0);
std::vector<float>	muon_chargedHadronIso(50,0);
std::vector<double>	muon_combinedMuon_chi2(50,0);
std::vector<double>	muon_combinedMuon_ndof(50,0);
std::vector<unsigned short>	muon_combinedMuon_numberOfValidHits(50,0);
std::vector<double>	muon_dB(50,0);
std::vector<float>	muon_ecalIso(50,0);
std::vector<double>	muon_energy(50,0);
std::vector<double>	muon_eta(50,0);
std::vector<double>	muon_genLepton_eta(50,0);
std::vector<int>	muon_genLepton_pdgId(50,0);
std::vector<double>	muon_genLepton_phi(50,0);
std::vector<double>	muon_genLepton_pt(50,0);
std::vector<float>	muon_globalTrack_pt(50,0);
std::vector<float>	muon_hcalIso(50,0);
std::vector<unsigned short>	muon_innerTrack_numberOfValidHits(50,0);
std::vector<float>	muon_innerTrack_pt(50,0);
std::vector<int>	muon_isGlobalMuon(50,0);
std::vector<int>	muon_isTrackerMuon(50,0);
std::vector<float>	muon_neutralHadronIso(50,0);
std::vector<double>	muon_phi(50,0);
std::vector<float>	muon_photonIso(50,0);
std::vector<double>	muon_pt(50,0);
std::vector<double>	muon_px(50,0);
std::vector<double>	muon_py(50,0);
std::vector<double>	muon_pz(50,0);
std::vector<float>	muon_trackIso(50,0);
std::vector<double>	muon_track_d0(50,0);
std::vector<unsigned short>	muon_track_hitPattern_numberOfValidMuonHits(50,0);
std::vector<unsigned short>	muon_track_hitPattern_numberOfValidPixelHits(50,0);
std::vector<double>	muon_track_normalizedChi2(50,0);
std::vector<double>	muon_vz(50,0);
std::vector<double>	muonhelper_dxywrtBeamSpot(50,0);
std::vector<double>	muonhelper_muonPFcandptdiff(50,0);
int	nelectron;
int	nelectronhelper;
int	ngeneventinfoproducthelper;
int	ngenparticlehelperra2b;
int	njet;
int	njet2;
int	njethelper;
int	nmet;
int	nmet2;
int	nmet3;
int	nmuon;
int	nmuonhelper;
int	ntau;
int	ntauhelper;
int	nvertex;
std::vector<float>	tau_caloIso(50,0);
std::vector<float>	tau_ecalIso(50,0);
std::vector<double>	tau_energy(50,0);
std::vector<double>	tau_et(50,0);
std::vector<double>	tau_eta(50,0);
std::vector<float>	tau_hcalIso(50,0);
std::vector<double>	tau_phi(50,0);
std::vector<double>	tau_pt(50,0);
std::vector<double>	tau_px(50,0);
std::vector<double>	tau_py(50,0);
std::vector<double>	tau_pz(50,0);
std::vector<float>	tau_tauID_againstElectron(50,0);
std::vector<float>	tau_tauID_againstMuon(50,0);
std::vector<float>	tau_tauID_byIsolation(50,0);
std::vector<float>	tau_tauID_byTaNC(50,0);
std::vector<float>	tau_tauID_byTaNCfrHalfPercent(50,0);
std::vector<float>	tau_tauID_byTaNCfrQuarterPercent(50,0);
std::vector<float>	tau_trackIso(50,0);
std::vector<int>	tauhelper_genTauDecayModeID(50,0);
int	triggerresultshelper_HLT_HT100U;
unsigned int	triggerresultshelper_HLT_HT100U_prs;
int	triggerresultshelper_HLT_HT120U;
unsigned int	triggerresultshelper_HLT_HT120U_prs;
int	triggerresultshelper_HLT_HT140U;
unsigned int	triggerresultshelper_HLT_HT140U_prs;
int	triggerresultshelper_HLT_HT150U_v3;
unsigned int	triggerresultshelper_HLT_HT150U_v3_prs;
int	triggerresultshelper_HLT_HT150_v1;
unsigned int	triggerresultshelper_HLT_HT150_v1_prs;
int	triggerresultshelper_HLT_HT150_v2;
unsigned int	triggerresultshelper_HLT_HT150_v2_prs;
int	triggerresultshelper_HLT_HT150_v3;
unsigned int	triggerresultshelper_HLT_HT150_v3_prs;
int	triggerresultshelper_HLT_HT150_v4;
unsigned int	triggerresultshelper_HLT_HT150_v4_prs;
int	triggerresultshelper_HLT_HT150_v5;
unsigned int	triggerresultshelper_HLT_HT150_v5_prs;
int	triggerresultshelper_HLT_HT250_MHT60_v2;
unsigned int	triggerresultshelper_HLT_HT250_MHT60_v2_prs;
int	triggerresultshelper_HLT_HT250_MHT60_v3;
unsigned int	triggerresultshelper_HLT_HT250_MHT60_v3_prs;
int	triggerresultshelper_HLT_HT250_v1;
unsigned int	triggerresultshelper_HLT_HT250_v1_prs;
int	triggerresultshelper_HLT_HT250_v2;
unsigned int	triggerresultshelper_HLT_HT250_v2_prs;
int	triggerresultshelper_HLT_HT250_v3;
unsigned int	triggerresultshelper_HLT_HT250_v3_prs;
int	triggerresultshelper_HLT_HT250_v4;
unsigned int	triggerresultshelper_HLT_HT250_v4_prs;
int	triggerresultshelper_HLT_HT250_v5;
unsigned int	triggerresultshelper_HLT_HT250_v5_prs;
int	triggerresultshelper_HLT_HT260_MHT60_v2;
unsigned int	triggerresultshelper_HLT_HT260_MHT60_v2_prs;
int	triggerresultshelper_HLT_HT260_v2;
unsigned int	triggerresultshelper_HLT_HT260_v2_prs;
int	triggerresultshelper_HLT_HT260_v3;
unsigned int	triggerresultshelper_HLT_HT260_v3_prs;
int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v2;
unsigned int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v2_prs;
int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3;
unsigned int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3_prs;
int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v4;
unsigned int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v4_prs;
int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v5;
unsigned int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v5_prs;
int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2;
unsigned int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2_prs;
int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3;
unsigned int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3_prs;
int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4;
unsigned int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4_prs;
int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5;
unsigned int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5_prs;
int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v2;
unsigned int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v2_prs;
int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v3;
unsigned int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v3_prs;
int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v4;
unsigned int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v4_prs;
int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v5;
unsigned int	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v5_prs;
int	triggerresultshelper_HLT_HT300_MHT75_v2;
unsigned int	triggerresultshelper_HLT_HT300_MHT75_v2_prs;
int	triggerresultshelper_HLT_HT300_MHT75_v3;
unsigned int	triggerresultshelper_HLT_HT300_MHT75_v3_prs;
int	triggerresultshelper_HLT_HT300_MHT75_v4;
unsigned int	triggerresultshelper_HLT_HT300_MHT75_v4_prs;
int	triggerresultshelper_HLT_HT300_PFMHT55_v2;
unsigned int	triggerresultshelper_HLT_HT300_PFMHT55_v2_prs;
int	triggerresultshelper_HLT_HT300_PFMHT55_v3;
unsigned int	triggerresultshelper_HLT_HT300_PFMHT55_v3_prs;
int	triggerresultshelper_HLT_HT300_PFMHT55_v4;
unsigned int	triggerresultshelper_HLT_HT300_PFMHT55_v4_prs;
int	triggerresultshelper_HLT_HT300_PFMHT55_v5;
unsigned int	triggerresultshelper_HLT_HT300_PFMHT55_v5_prs;
int	triggerresultshelper_HLT_HT300_v1;
unsigned int	triggerresultshelper_HLT_HT300_v1_prs;
int	triggerresultshelper_HLT_HT300_v2;
unsigned int	triggerresultshelper_HLT_HT300_v2_prs;
int	triggerresultshelper_HLT_HT300_v3;
unsigned int	triggerresultshelper_HLT_HT300_v3_prs;
int	triggerresultshelper_HLT_HT300_v4;
unsigned int	triggerresultshelper_HLT_HT300_v4_prs;
int	triggerresultshelper_HLT_HT300_v5;
unsigned int	triggerresultshelper_HLT_HT300_v5_prs;
int	triggerresultshelper_HLT_HT300_v6;
unsigned int	triggerresultshelper_HLT_HT300_v6_prs;
int	triggerresultshelper_HLT_HT300_v7;
unsigned int	triggerresultshelper_HLT_HT300_v7_prs;
std::vector<int>	vertex_isFake(50,0);
std::vector<int>	vertex_isValid(50,0);
std::vector<double>	vertex_ndof(50,0);
std::vector<double>	vertex_normalizedChi2(50,0);
std::vector<double>	vertex_position_Rho(50,0);
std::vector<int>	vertex_tracksSize(50,0);
std::vector<double>	vertex_x(50,0);
std::vector<double>	vertex_xError(50,0);
std::vector<double>	vertex_y(50,0);
std::vector<double>	vertex_yError(50,0);
std::vector<double>	vertex_z(50,0);
std::vector<double>	vertex_zError(50,0);


//-----------------------------------------------------------------------------
// -- Select variables to be read
//-----------------------------------------------------------------------------
void selectVariables(itreestream& stream)
{
  stream.select("sdouble_kt6PFJets_sigma.value", sigma);
  stream.select("sdouble_kt6PFJets_rho.value", rho);
  stream.select("suint_flavorHistoryFilter.value", flavorhistorypath);
  stream.select("recoBeamSpot_offlineBeamSpot.x0", beamspot_x0);
  stream.select("recoBeamSpot_offlineBeamSpot.y0", beamspot_y0);
  stream.select("recoBeamSpot_offlineBeamSpot.z0", beamspot_z0);
  stream.select("patElectron_selectedPatElectronsPF.caloIso", electron_caloIso);
  stream.select("patElectron_selectedPatElectronsPF.charge", electron_charge);
  stream.select("patElectron_selectedPatElectronsPF.chargedHadronIso", electron_chargedHadronIso);
  stream.select("patElectron_selectedPatElectronsPF.dB", electron_dB);
  stream.select("patElectron_selectedPatElectronsPF.deltaEtaSuperClusterTrackAtVtx", electron_deltaEtaSuperClusterTrackAtVtx);
  stream.select("patElectron_selectedPatElectronsPF.deltaPhiSuperClusterTrackAtVtx", electron_deltaPhiSuperClusterTrackAtVtx);
  stream.select("patElectron_selectedPatElectronsPF.dr03EcalRecHitSumEt", electron_dr03EcalRecHitSumEt);
  stream.select("patElectron_selectedPatElectronsPF.dr03HcalTowerSumEt", electron_dr03HcalTowerSumEt);
  stream.select("patElectron_selectedPatElectronsPF.dr03TkSumPt", electron_dr03TkSumPt);
  stream.select("patElectron_selectedPatElectronsPF.ecalIso", electron_ecalIso);
  stream.select("patElectron_selectedPatElectronsPF.eidRobustTight", electron_eidRobustTight);
  stream.select("patElectron_selectedPatElectronsPF.energy", electron_energy);
  stream.select("patElectron_selectedPatElectronsPF.et", electron_et);
  stream.select("patElectron_selectedPatElectronsPF.eta", electron_eta);
  stream.select("patElectron_selectedPatElectronsPF.genLepton_eta", electron_genLepton_eta);
  stream.select("patElectron_selectedPatElectronsPF.genLepton_pdgId", electron_genLepton_pdgId);
  stream.select("patElectron_selectedPatElectronsPF.genLepton_phi", electron_genLepton_phi);
  stream.select("patElectron_selectedPatElectronsPF.genLepton_pt", electron_genLepton_pt);
  stream.select("patElectron_selectedPatElectronsPF.gsfTrack_d0", electron_gsfTrack_d0);
  stream.select("patElectron_selectedPatElectronsPF.gsfTrack_trackerExpectedHitsInner_numberOfLostHits", electron_gsfTrack_trackerExpectedHitsInner_numberOfLostHits);
  stream.select("patElectron_selectedPatElectronsPF.hadronicOverEm", electron_hadronicOverEm);
  stream.select("patElectron_selectedPatElectronsPF.hcalIso", electron_hcalIso);
  stream.select("patElectron_selectedPatElectronsPF.neutralHadronIso", electron_neutralHadronIso);
  stream.select("patElectron_selectedPatElectronsPF.phi", electron_phi);
  stream.select("patElectron_selectedPatElectronsPF.photonIso", electron_photonIso);
  stream.select("patElectron_selectedPatElectronsPF.pt", electron_pt);
  stream.select("patElectron_selectedPatElectronsPF.px", electron_px);
  stream.select("patElectron_selectedPatElectronsPF.py", electron_py);
  stream.select("patElectron_selectedPatElectronsPF.pz", electron_pz);
  stream.select("patElectron_selectedPatElectronsPF.sigmaIetaIeta", electron_sigmaIetaIeta);
  stream.select("patElectron_selectedPatElectronsPF.simpleEleId95relIso", electron_simpleEleId95relIso);
  stream.select("patElectron_selectedPatElectronsPF.superCluster_eta", electron_superCluster_eta);
  stream.select("patElectron_selectedPatElectronsPF.trackIso", electron_trackIso);
  stream.select("patElectron_selectedPatElectronsPF.vz", electron_vz);
  stream.select("patElectronHelper_selectedPatElectronsPF.dxywrtBeamSpot", electronhelper_dxywrtBeamSpot);
  stream.select("edmEventHelper_info.bunchCrossing", eventhelper_bunchCrossing);
  stream.select("edmEventHelper_info.event", eventhelper_event);
  stream.select("edmEventHelper_info.isRealData", eventhelper_isRealData);
  stream.select("edmEventHelper_info.luminosityBlock", eventhelper_luminosityBlock);
  stream.select("edmEventHelper_info.run", eventhelper_run);
  stream.select("GenEventInfoProduct_generator.pdf1", geneventinfoproduct_pdf1);
  stream.select("GenEventInfoProduct_generator.pdf2", geneventinfoproduct_pdf2);
  stream.select("GenEventInfoProduct_generator.scalePDF", geneventinfoproduct_scalePDF);
  stream.select("GenEventInfoProduct_generator.weight", geneventinfoproduct_weight);
  stream.select("GenEventInfoProduct_generator.x1", geneventinfoproduct_x1);
  stream.select("GenEventInfoProduct_generator.x2", geneventinfoproduct_x2);
  stream.select("GenEventInfoProductHelper_generator.pdf1", geneventinfoproducthelper_pdf1);
  stream.select("GenEventInfoProductHelper_generator.pdf2", geneventinfoproducthelper_pdf2);
  stream.select("GenEventInfoProductHelper_generator.pdfweight", geneventinfoproducthelper_pdfweight);
  stream.select("GenEventInfoProductHelper_generator.pdfweightsum", geneventinfoproducthelper_pdfweightsum);
  stream.select("recoGenParticleHelperRA2b_genParticles.charge", genparticlehelperra2b_charge);
  stream.select("recoGenParticleHelperRA2b_genParticles.eta", genparticlehelperra2b_eta);
  stream.select("recoGenParticleHelperRA2b_genParticles.firstDaughter", genparticlehelperra2b_firstDaughter);
  stream.select("recoGenParticleHelperRA2b_genParticles.firstMother", genparticlehelperra2b_firstMother);
  stream.select("recoGenParticleHelperRA2b_genParticles.lastDaughter", genparticlehelperra2b_lastDaughter);
  stream.select("recoGenParticleHelperRA2b_genParticles.lastMother", genparticlehelperra2b_lastMother);
  stream.select("recoGenParticleHelperRA2b_genParticles.mass", genparticlehelperra2b_mass);
  stream.select("recoGenParticleHelperRA2b_genParticles.pdgId", genparticlehelperra2b_pdgId);
  stream.select("recoGenParticleHelperRA2b_genParticles.phi", genparticlehelperra2b_phi);
  stream.select("recoGenParticleHelperRA2b_genParticles.pt", genparticlehelperra2b_pt);
  stream.select("recoGenParticleHelperRA2b_genParticles.status", genparticlehelperra2b_status);
  stream.select("GenRunInfoProduct_generator.externalXSecLO_error", genruninfoproduct_externalXSecLO_error);
  stream.select("GenRunInfoProduct_generator.externalXSecLO_value", genruninfoproduct_externalXSecLO_value);
  stream.select("GenRunInfoProduct_generator.filterEfficiency", genruninfoproduct_filterEfficiency);
  stream.select("GenRunInfoProduct_generator.internalXSec_error", genruninfoproduct_internalXSec_error);
  stream.select("GenRunInfoProduct_generator.internalXSec_value", genruninfoproduct_internalXSec_value);
  stream.select("patJet_cleanPatJetsAK5Calo.combinedSecondaryVertexBJetTags", jet2_combinedSecondaryVertexBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.combinedSecondaryVertexMVABJetTags", jet2_combinedSecondaryVertexMVABJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.energy", jet2_energy);
  stream.select("patJet_cleanPatJetsAK5Calo.eta", jet2_eta);
  stream.select("patJet_cleanPatJetsAK5Calo.jetBProbabilityBJetTags", jet2_jetBProbabilityBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.jetID_fHPD", jet2_jetID_fHPD);
  stream.select("patJet_cleanPatJetsAK5Calo.jetID_n90Hits", jet2_jetID_n90Hits);
  stream.select("patJet_cleanPatJetsAK5Calo.jetProbabilityBJetTags", jet2_jetProbabilityBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.partonFlavour", jet2_partonFlavour);
  stream.select("patJet_cleanPatJetsAK5Calo.phi", jet2_phi);
  stream.select("patJet_cleanPatJetsAK5Calo.pt", jet2_pt);
  stream.select("patJet_cleanPatJetsAK5Calo.simpleSecondaryVertexBJetTags", jet2_simpleSecondaryVertexBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.simpleSecondaryVertexHighEffBJetTags", jet2_simpleSecondaryVertexHighEffBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.simpleSecondaryVertexHighPurBJetTags", jet2_simpleSecondaryVertexHighPurBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.softElectronByIP3dBJetTags", jet2_softElectronByIP3dBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.softElectronByPtBJetTags", jet2_softElectronByPtBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.softMuonBJetTags", jet2_softMuonBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.softMuonByIP3dBJetTags", jet2_softMuonByIP3dBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.softMuonByPtBJetTags", jet2_softMuonByPtBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.trackCountingHighEffBJetTags", jet2_trackCountingHighEffBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.trackCountingHighPurBJetTags", jet2_trackCountingHighPurBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.uncor_energy", jet2_uncor_energy);
  stream.select("patJet_cleanPatJetsAK5Calo.uncor_eta", jet2_uncor_eta);
  stream.select("patJet_cleanPatJetsAK5Calo.uncor_phi", jet2_uncor_phi);
  stream.select("patJet_cleanPatJetsAK5Calo.uncor_pt", jet2_uncor_pt);
  stream.select("patJet_selectedPatJetsPF.HFHadronEnergy", jet_HFHadronEnergy);
  stream.select("patJet_selectedPatJetsPF.chargedEmEnergyFraction", jet_chargedEmEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.chargedHadronEnergyFraction", jet_chargedHadronEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.chargedMultiplicity", jet_chargedMultiplicity);
  stream.select("patJet_selectedPatJetsPF.combinedSecondaryVertexBJetTags", jet_combinedSecondaryVertexBJetTags);
  stream.select("patJet_selectedPatJetsPF.combinedSecondaryVertexMVABJetTags", jet_combinedSecondaryVertexMVABJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.emEnergyFraction", jet_emEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.energy", jet_energy);
  stream.select("patJet_selectedPatJetsPF.eta", jet_eta);
  stream.select("patJet_selectedPatJetsPF.jetArea", jet_jetArea);
  stream.select("patJet_selectedPatJetsPF.jetBProbabilityBJetTags", jet_jetBProbabilityBJetTags);
  stream.select("patJet_selectedPatJetsPF.jetID_fHPD", jet_jetID_fHPD);
  stream.select("patJet_selectedPatJetsPF.jetID_n90Hits", jet_jetID_n90Hits);
  stream.select("patJet_selectedPatJetsPF.jetProbabilityBJetTags", jet_jetProbabilityBJetTags);
  stream.select("patJet_selectedPatJetsPF.neutralEmEnergyFraction", jet_neutralEmEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.neutralHadronEnergy", jet_neutralHadronEnergy);
  stream.select("patJet_selectedPatJetsPF.neutralHadronEnergyFraction", jet_neutralHadronEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.numberOfDaughters", jet_numberOfDaughters);
  stream.select("patJet_selectedPatJetsPF.partonFlavour", jet_partonFlavour);
  stream.select("patJet_selectedPatJetsPF.phi", jet_phi);
  stream.select("patJet_selectedPatJetsPF.pt", jet_pt);
  stream.select("patJet_selectedPatJetsPF.simpleSecondaryVertexBJetTags", jet_simpleSecondaryVertexBJetTags);
  stream.select("patJet_selectedPatJetsPF.simpleSecondaryVertexHighEffBJetTags", jet_simpleSecondaryVertexHighEffBJetTags);
  stream.select("patJet_selectedPatJetsPF.simpleSecondaryVertexHighPurBJetTags", jet_simpleSecondaryVertexHighPurBJetTags);
  stream.select("patJet_selectedPatJetsPF.softElectronByIP3dBJetTags", jet_softElectronByIP3dBJetTags);
  stream.select("patJet_selectedPatJetsPF.softElectronByPtBJetTags", jet_softElectronByPtBJetTags);
  stream.select("patJet_selectedPatJetsPF.softMuonBJetTags", jet_softMuonBJetTags);
  stream.select("patJet_selectedPatJetsPF.softMuonByIP3dBJetTags", jet_softMuonByIP3dBJetTags);
  stream.select("patJet_selectedPatJetsPF.softMuonByPtBJetTags", jet_softMuonByPtBJetTags);
  stream.select("patJet_selectedPatJetsPF.trackCountingHighEffBJetTags", jet_trackCountingHighEffBJetTags);
  stream.select("patJet_selectedPatJetsPF.trackCountingHighPurBJetTags", jet_trackCountingHighPurBJetTags);
  stream.select("patJet_selectedPatJetsPF.uncor_energy", jet_uncor_energy);
  stream.select("patJet_selectedPatJetsPF.uncor_eta", jet_uncor_eta);
  stream.select("patJet_selectedPatJetsPF.uncor_phi", jet_uncor_phi);
  stream.select("patJet_selectedPatJetsPF.uncor_pt", jet_uncor_pt);
  stream.select("patJetHelper_selectedPatJetsPF.jecFactor", jethelper_jecFactor);
  stream.select("patJetHelper_selectedPatJetsPF.jecFactorNoL1Fast", jethelper_jecFactorNoL1Fast);
  stream.select("patJetHelper_selectedPatJetsPF.jetUncMinus", jethelper_jetUncMinus);
  stream.select("patJetHelper_selectedPatJetsPF.jetUncPlus", jethelper_jetUncPlus);
  stream.select("patMET_patMETsTypeIPF.energy", met2_energy);
  stream.select("patMET_patMETsTypeIPF.et", met2_et);
  stream.select("patMET_patMETsTypeIPF.mEtSig", met2_mEtSig);
  stream.select("patMET_patMETsTypeIPF.phi", met2_phi);
  stream.select("patMET_patMETsTypeIPF.pt", met2_pt);
  stream.select("patMET_patMETsTypeIPF.sumEt", met2_sumEt);
  stream.select("patMET_patMETsAK5Calo.energy", met3_energy);
  stream.select("patMET_patMETsAK5Calo.et", met3_et);
  stream.select("patMET_patMETsAK5Calo.mEtSig", met3_mEtSig);
  stream.select("patMET_patMETsAK5Calo.phi", met3_phi);
  stream.select("patMET_patMETsAK5Calo.pt", met3_pt);
  stream.select("patMET_patMETsAK5Calo.sumEt", met3_sumEt);
  stream.select("patMET_patMETsPF.energy", met_energy);
  stream.select("patMET_patMETsPF.et", met_et);
  stream.select("patMET_patMETsPF.mEtSig", met_mEtSig);
  stream.select("patMET_patMETsPF.phi", met_phi);
  stream.select("patMET_patMETsPF.pt", met_pt);
  stream.select("patMET_patMETsPF.sumEt", met_sumEt);
  stream.select("patMuon_selectedPatMuonsPF.AllGlobalMuons", muon_AllGlobalMuons);
  stream.select("patMuon_selectedPatMuonsPF.GlobalMuonPromptTight", muon_GlobalMuonPromptTight);
  stream.select("patMuon_selectedPatMuonsPF.charge", muon_charge);
  stream.select("patMuon_selectedPatMuonsPF.chargedHadronIso", muon_chargedHadronIso);
  stream.select("patMuon_selectedPatMuonsPF.combinedMuon_chi2", muon_combinedMuon_chi2);
  stream.select("patMuon_selectedPatMuonsPF.combinedMuon_ndof", muon_combinedMuon_ndof);
  stream.select("patMuon_selectedPatMuonsPF.combinedMuon_numberOfValidHits", muon_combinedMuon_numberOfValidHits);
  stream.select("patMuon_selectedPatMuonsPF.dB", muon_dB);
  stream.select("patMuon_selectedPatMuonsPF.ecalIso", muon_ecalIso);
  stream.select("patMuon_selectedPatMuonsPF.energy", muon_energy);
  stream.select("patMuon_selectedPatMuonsPF.eta", muon_eta);
  stream.select("patMuon_selectedPatMuonsPF.genLepton_eta", muon_genLepton_eta);
  stream.select("patMuon_selectedPatMuonsPF.genLepton_pdgId", muon_genLepton_pdgId);
  stream.select("patMuon_selectedPatMuonsPF.genLepton_phi", muon_genLepton_phi);
  stream.select("patMuon_selectedPatMuonsPF.genLepton_pt", muon_genLepton_pt);
  stream.select("patMuon_selectedPatMuonsPF.globalTrack_pt", muon_globalTrack_pt);
  stream.select("patMuon_selectedPatMuonsPF.hcalIso", muon_hcalIso);
  stream.select("patMuon_selectedPatMuonsPF.innerTrack_numberOfValidHits", muon_innerTrack_numberOfValidHits);
  stream.select("patMuon_selectedPatMuonsPF.innerTrack_pt", muon_innerTrack_pt);
  stream.select("patMuon_selectedPatMuonsPF.isGlobalMuon", muon_isGlobalMuon);
  stream.select("patMuon_selectedPatMuonsPF.isTrackerMuon", muon_isTrackerMuon);
  stream.select("patMuon_selectedPatMuonsPF.neutralHadronIso", muon_neutralHadronIso);
  stream.select("patMuon_selectedPatMuonsPF.phi", muon_phi);
  stream.select("patMuon_selectedPatMuonsPF.photonIso", muon_photonIso);
  stream.select("patMuon_selectedPatMuonsPF.pt", muon_pt);
  stream.select("patMuon_selectedPatMuonsPF.px", muon_px);
  stream.select("patMuon_selectedPatMuonsPF.py", muon_py);
  stream.select("patMuon_selectedPatMuonsPF.pz", muon_pz);
  stream.select("patMuon_selectedPatMuonsPF.trackIso", muon_trackIso);
  stream.select("patMuon_selectedPatMuonsPF.track_d0", muon_track_d0);
  stream.select("patMuon_selectedPatMuonsPF.track_hitPattern_numberOfValidMuonHits", muon_track_hitPattern_numberOfValidMuonHits);
  stream.select("patMuon_selectedPatMuonsPF.track_hitPattern_numberOfValidPixelHits", muon_track_hitPattern_numberOfValidPixelHits);
  stream.select("patMuon_selectedPatMuonsPF.track_normalizedChi2", muon_track_normalizedChi2);
  stream.select("patMuon_selectedPatMuonsPF.vz", muon_vz);
  stream.select("patMuonHelper_selectedPatMuonsPF.dxywrtBeamSpot", muonhelper_dxywrtBeamSpot);
  stream.select("patMuonHelper_selectedPatMuonsPF.muonPFcandptdiff", muonhelper_muonPFcandptdiff);
  stream.select("npatElectron_selectedPatElectronsPF", nelectron);
  stream.select("npatElectronHelper_selectedPatElectronsPF", nelectronhelper);
  stream.select("nGenEventInfoProductHelper_generator", ngeneventinfoproducthelper);
  stream.select("nrecoGenParticleHelperRA2b_genParticles", ngenparticlehelperra2b);
  stream.select("npatJet_selectedPatJetsPF", njet);
  stream.select("npatJet_cleanPatJetsAK5Calo", njet2);
  stream.select("npatJetHelper_selectedPatJetsPF", njethelper);
  stream.select("npatMET_patMETsPF", nmet);
  stream.select("npatMET_patMETsTypeIPF", nmet2);
  stream.select("npatMET_patMETsAK5Calo", nmet3);
  stream.select("npatMuon_selectedPatMuonsPF", nmuon);
  stream.select("npatMuonHelper_selectedPatMuonsPF", nmuonhelper);
  stream.select("npatTau_selectedPatTausPF", ntau);
  stream.select("npatTauHelper_selectedPatTausPF", ntauhelper);
  stream.select("nrecoVertex_offlinePrimaryVertices", nvertex);
  stream.select("patTau_selectedPatTausPF.caloIso", tau_caloIso);
  stream.select("patTau_selectedPatTausPF.ecalIso", tau_ecalIso);
  stream.select("patTau_selectedPatTausPF.energy", tau_energy);
  stream.select("patTau_selectedPatTausPF.et", tau_et);
  stream.select("patTau_selectedPatTausPF.eta", tau_eta);
  stream.select("patTau_selectedPatTausPF.hcalIso", tau_hcalIso);
  stream.select("patTau_selectedPatTausPF.phi", tau_phi);
  stream.select("patTau_selectedPatTausPF.pt", tau_pt);
  stream.select("patTau_selectedPatTausPF.px", tau_px);
  stream.select("patTau_selectedPatTausPF.py", tau_py);
  stream.select("patTau_selectedPatTausPF.pz", tau_pz);
  stream.select("patTau_selectedPatTausPF.tauID_againstElectron", tau_tauID_againstElectron);
  stream.select("patTau_selectedPatTausPF.tauID_againstMuon", tau_tauID_againstMuon);
  stream.select("patTau_selectedPatTausPF.tauID_byIsolation", tau_tauID_byIsolation);
  stream.select("patTau_selectedPatTausPF.tauID_byTaNC", tau_tauID_byTaNC);
  stream.select("patTau_selectedPatTausPF.tauID_byTaNCfrHalfPercent", tau_tauID_byTaNCfrHalfPercent);
  stream.select("patTau_selectedPatTausPF.tauID_byTaNCfrQuarterPercent", tau_tauID_byTaNCfrQuarterPercent);
  stream.select("patTau_selectedPatTausPF.trackIso", tau_trackIso);
  stream.select("patTauHelper_selectedPatTausPF.genTauDecayModeID", tauhelper_genTauDecayModeID);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT100U", triggerresultshelper_HLT_HT100U);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT100U_prs", triggerresultshelper_HLT_HT100U_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT120U", triggerresultshelper_HLT_HT120U);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT120U_prs", triggerresultshelper_HLT_HT120U_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT140U", triggerresultshelper_HLT_HT140U);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT140U_prs", triggerresultshelper_HLT_HT140U_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150U_v3", triggerresultshelper_HLT_HT150U_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150U_v3_prs", triggerresultshelper_HLT_HT150U_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v1", triggerresultshelper_HLT_HT150_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v1_prs", triggerresultshelper_HLT_HT150_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v2", triggerresultshelper_HLT_HT150_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v2_prs", triggerresultshelper_HLT_HT150_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v3", triggerresultshelper_HLT_HT150_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v3_prs", triggerresultshelper_HLT_HT150_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v4", triggerresultshelper_HLT_HT150_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v4_prs", triggerresultshelper_HLT_HT150_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v5", triggerresultshelper_HLT_HT150_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v5_prs", triggerresultshelper_HLT_HT150_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT60_v2", triggerresultshelper_HLT_HT250_MHT60_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT60_v2_prs", triggerresultshelper_HLT_HT250_MHT60_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT60_v3", triggerresultshelper_HLT_HT250_MHT60_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT60_v3_prs", triggerresultshelper_HLT_HT250_MHT60_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v1", triggerresultshelper_HLT_HT250_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v1_prs", triggerresultshelper_HLT_HT250_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v2", triggerresultshelper_HLT_HT250_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v2_prs", triggerresultshelper_HLT_HT250_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v3", triggerresultshelper_HLT_HT250_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v3_prs", triggerresultshelper_HLT_HT250_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v4", triggerresultshelper_HLT_HT250_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v4_prs", triggerresultshelper_HLT_HT250_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v5", triggerresultshelper_HLT_HT250_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v5_prs", triggerresultshelper_HLT_HT250_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT260_MHT60_v2", triggerresultshelper_HLT_HT260_MHT60_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT260_MHT60_v2_prs", triggerresultshelper_HLT_HT260_MHT60_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT260_v2", triggerresultshelper_HLT_HT260_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT260_v2_prs", triggerresultshelper_HLT_HT260_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT260_v3", triggerresultshelper_HLT_HT260_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT260_v3_prs", triggerresultshelper_HLT_HT260_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v2", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v2_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v4", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v4_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v5", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v5_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v2", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v2_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v3", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v3_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v4", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v4_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v5", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v5_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v2", triggerresultshelper_HLT_HT300_MHT75_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v2_prs", triggerresultshelper_HLT_HT300_MHT75_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v3", triggerresultshelper_HLT_HT300_MHT75_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v3_prs", triggerresultshelper_HLT_HT300_MHT75_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v4", triggerresultshelper_HLT_HT300_MHT75_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v4_prs", triggerresultshelper_HLT_HT300_MHT75_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v2", triggerresultshelper_HLT_HT300_PFMHT55_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v2_prs", triggerresultshelper_HLT_HT300_PFMHT55_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v3", triggerresultshelper_HLT_HT300_PFMHT55_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v3_prs", triggerresultshelper_HLT_HT300_PFMHT55_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v4", triggerresultshelper_HLT_HT300_PFMHT55_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v4_prs", triggerresultshelper_HLT_HT300_PFMHT55_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v5", triggerresultshelper_HLT_HT300_PFMHT55_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v5_prs", triggerresultshelper_HLT_HT300_PFMHT55_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v1", triggerresultshelper_HLT_HT300_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v1_prs", triggerresultshelper_HLT_HT300_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v2", triggerresultshelper_HLT_HT300_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v2_prs", triggerresultshelper_HLT_HT300_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v3", triggerresultshelper_HLT_HT300_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v3_prs", triggerresultshelper_HLT_HT300_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v4", triggerresultshelper_HLT_HT300_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v4_prs", triggerresultshelper_HLT_HT300_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v5", triggerresultshelper_HLT_HT300_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v5_prs", triggerresultshelper_HLT_HT300_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v6", triggerresultshelper_HLT_HT300_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v6_prs", triggerresultshelper_HLT_HT300_v6_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v7", triggerresultshelper_HLT_HT300_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v7_prs", triggerresultshelper_HLT_HT300_v7_prs);
  stream.select("recoVertex_offlinePrimaryVertices.isFake", vertex_isFake);
  stream.select("recoVertex_offlinePrimaryVertices.isValid", vertex_isValid);
  stream.select("recoVertex_offlinePrimaryVertices.ndof", vertex_ndof);
  stream.select("recoVertex_offlinePrimaryVertices.normalizedChi2", vertex_normalizedChi2);
  stream.select("recoVertex_offlinePrimaryVertices.position_Rho", vertex_position_Rho);
  stream.select("recoVertex_offlinePrimaryVertices.tracksSize", vertex_tracksSize);
  stream.select("recoVertex_offlinePrimaryVertices.x", vertex_x);
  stream.select("recoVertex_offlinePrimaryVertices.xError", vertex_xError);
  stream.select("recoVertex_offlinePrimaryVertices.y", vertex_y);
  stream.select("recoVertex_offlinePrimaryVertices.yError", vertex_yError);
  stream.select("recoVertex_offlinePrimaryVertices.z", vertex_z);
  stream.select("recoVertex_offlinePrimaryVertices.zError", vertex_zError);

}

//-----------------------------------------------------------------------------
// --- These structs can be filled by calling fillObjects()
// --- after the call to stream.read(...)
//-----------------------------------------------------------------------------
struct electron_s
{
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
  double	px;
  double	py;
  double	pz;
  float	trackIso;
  float	ecalIso;
  float	hcalIso;
  float	caloIso;
  float	eidRobustTight;
  float	simpleEleId95relIso;
  double	gsfTrack_d0;
  float	dr03TkSumPt;
  float	dr03EcalRecHitSumEt;
  float	dr03HcalTowerSumEt;
  double	vz;
  double	dB;
  float	chargedHadronIso;
  float	photonIso;
  float	neutralHadronIso;
  int	gsfTrack_trackerExpectedHitsInner_numberOfLostHits;
  double	superCluster_eta;
  int	charge;
  double	genLepton_pt;
  double	genLepton_eta;
  double	genLepton_phi;
  int	genLepton_pdgId;
  double	hadronicOverEm;
  double	deltaPhiSuperClusterTrackAtVtx;
  double	deltaEtaSuperClusterTrackAtVtx;
  double	sigmaIetaIeta;
};
std::vector<electron_s> electron(50);

std::ostream& operator<<(std::ostream& os, const electron_s& o)
{
  char r[1024];
  os << "electron" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "px", (double)o.px); os << r;
  sprintf(r, "  %-32s: %f\n", "py", (double)o.py); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "trackIso", (double)o.trackIso); os << r;
  sprintf(r, "  %-32s: %f\n", "ecalIso", (double)o.ecalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "hcalIso", (double)o.hcalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "caloIso", (double)o.caloIso); os << r;
  sprintf(r, "  %-32s: %f\n", "eidRobustTight", (double)o.eidRobustTight); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleEleId95relIso", (double)o.simpleEleId95relIso); os << r;
  sprintf(r, "  %-32s: %f\n", "gsfTrack_d0", (double)o.gsfTrack_d0); os << r;
  sprintf(r, "  %-32s: %f\n", "dr03TkSumPt", (double)o.dr03TkSumPt); os << r;
  sprintf(r, "  %-32s: %f\n", "dr03EcalRecHitSumEt", (double)o.dr03EcalRecHitSumEt); os << r;
  sprintf(r, "  %-32s: %f\n", "dr03HcalTowerSumEt", (double)o.dr03HcalTowerSumEt); os << r;
  sprintf(r, "  %-32s: %f\n", "vz", (double)o.vz); os << r;
  sprintf(r, "  %-32s: %f\n", "dB", (double)o.dB); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronIso", (double)o.chargedHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "photonIso", (double)o.photonIso); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronIso", (double)o.neutralHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "gsfTrack_trackerExpectedHitsInner_numberOfLostHits", (double)o.gsfTrack_trackerExpectedHitsInner_numberOfLostHits); os << r;
  sprintf(r, "  %-32s: %f\n", "superCluster_eta", (double)o.superCluster_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_pt", (double)o.genLepton_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_eta", (double)o.genLepton_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_phi", (double)o.genLepton_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_pdgId", (double)o.genLepton_pdgId); os << r;
  sprintf(r, "  %-32s: %f\n", "hadronicOverEm", (double)o.hadronicOverEm); os << r;
  sprintf(r, "  %-32s: %f\n", "deltaPhiSuperClusterTrackAtVtx", (double)o.deltaPhiSuperClusterTrackAtVtx); os << r;
  sprintf(r, "  %-32s: %f\n", "deltaEtaSuperClusterTrackAtVtx", (double)o.deltaEtaSuperClusterTrackAtVtx); os << r;
  sprintf(r, "  %-32s: %f\n", "sigmaIetaIeta", (double)o.sigmaIetaIeta); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct electronhelper_s
{
  double	dxywrtBeamSpot;
};
std::vector<electronhelper_s> electronhelper(50);

std::ostream& operator<<(std::ostream& os, const electronhelper_s& o)
{
  char r[1024];
  os << "electronhelper" << std::endl;
  sprintf(r, "  %-32s: %f\n", "dxywrtBeamSpot", (double)o.dxywrtBeamSpot); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct geneventinfoproducthelper_s
{
  double	pdf1;
  double	pdf2;
  double	pdfweight;
  double	pdfweightsum;
};
std::vector<geneventinfoproducthelper_s> geneventinfoproducthelper(50);

std::ostream& operator<<(std::ostream& os, const geneventinfoproducthelper_s& o)
{
  char r[1024];
  os << "geneventinfoproducthelper" << std::endl;
  sprintf(r, "  %-32s: %f\n", "pdf1", (double)o.pdf1); os << r;
  sprintf(r, "  %-32s: %f\n", "pdf2", (double)o.pdf2); os << r;
  sprintf(r, "  %-32s: %f\n", "pdfweight", (double)o.pdfweight); os << r;
  sprintf(r, "  %-32s: %f\n", "pdfweightsum", (double)o.pdfweightsum); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct genparticlehelperra2b_s
{
  int	firstMother;
  int	lastMother;
  int	firstDaughter;
  int	lastDaughter;
  int	charge;
  int	pdgId;
  int	status;
  double	pt;
  double	eta;
  double	phi;
  double	mass;
};
std::vector<genparticlehelperra2b_s> genparticlehelperra2b(200);

std::ostream& operator<<(std::ostream& os, const genparticlehelperra2b_s& o)
{
  char r[1024];
  os << "genparticlehelperra2b" << std::endl;
  sprintf(r, "  %-32s: %f\n", "firstMother", (double)o.firstMother); os << r;
  sprintf(r, "  %-32s: %f\n", "lastMother", (double)o.lastMother); os << r;
  sprintf(r, "  %-32s: %f\n", "firstDaughter", (double)o.firstDaughter); os << r;
  sprintf(r, "  %-32s: %f\n", "lastDaughter", (double)o.lastDaughter); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "pdgId", (double)o.pdgId); os << r;
  sprintf(r, "  %-32s: %f\n", "status", (double)o.status); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "mass", (double)o.mass); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jet_s
{
  double	energy;
  double	uncor_energy;
  double	pt;
  double	uncor_pt;
  double	phi;
  double	uncor_phi;
  double	eta;
  double	uncor_eta;
  float	jetArea;
  double	jetID_fHPD;
  int	jetID_n90Hits;
  float	HFHadronEnergy;
  float	neutralHadronEnergy;
  float	neutralEmEnergyFraction;
  float	neutralHadronEnergyFraction;
  int	numberOfDaughters;
  float	chargedHadronEnergyFraction;
  int	chargedMultiplicity;
  float	chargedEmEnergyFraction;
  float	trackCountingHighEffBJetTags;
  float	trackCountingHighPurBJetTags;
  float	simpleSecondaryVertexBJetTags;
  float	simpleSecondaryVertexHighEffBJetTags;
  float	simpleSecondaryVertexHighPurBJetTags;
  float	combinedSecondaryVertexBJetTags;
  float	combinedSecondaryVertexMVABJetTags;
  float	jetBProbabilityBJetTags;
  float	jetProbabilityBJetTags;
  float	softElectronByIP3dBJetTags;
  float	softElectronByPtBJetTags;
  float	softMuonBJetTags;
  float	softMuonByIP3dBJetTags;
  float	softMuonByPtBJetTags;
  int	partonFlavour;
  float	emEnergyFraction;
};
std::vector<jet_s> jet(500);

std::ostream& operator<<(std::ostream& os, const jet_s& o)
{
  char r[1024];
  os << "jet" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_energy", (double)o.uncor_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_pt", (double)o.uncor_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_phi", (double)o.uncor_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_eta", (double)o.uncor_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "jetArea", (double)o.jetArea); os << r;
  sprintf(r, "  %-32s: %f\n", "jetID_fHPD", (double)o.jetID_fHPD); os << r;
  sprintf(r, "  %-32s: %f\n", "jetID_n90Hits", (double)o.jetID_n90Hits); os << r;
  sprintf(r, "  %-32s: %f\n", "HFHadronEnergy", (double)o.HFHadronEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergy", (double)o.neutralHadronEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralEmEnergyFraction", (double)o.neutralEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergyFraction", (double)o.neutralHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfDaughters", (double)o.numberOfDaughters); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronEnergyFraction", (double)o.chargedHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedMultiplicity", (double)o.chargedMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedEmEnergyFraction", (double)o.chargedEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighEffBJetTags", (double)o.trackCountingHighEffBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighPurBJetTags", (double)o.trackCountingHighPurBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleSecondaryVertexBJetTags", (double)o.simpleSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleSecondaryVertexHighEffBJetTags", (double)o.simpleSecondaryVertexHighEffBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleSecondaryVertexHighPurBJetTags", (double)o.simpleSecondaryVertexHighPurBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexBJetTags", (double)o.combinedSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexMVABJetTags", (double)o.combinedSecondaryVertexMVABJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetBProbabilityBJetTags", (double)o.jetBProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetProbabilityBJetTags", (double)o.jetProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softElectronByIP3dBJetTags", (double)o.softElectronByIP3dBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softElectronByPtBJetTags", (double)o.softElectronByPtBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softMuonBJetTags", (double)o.softMuonBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softMuonByIP3dBJetTags", (double)o.softMuonByIP3dBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softMuonByPtBJetTags", (double)o.softMuonByPtBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "partonFlavour", (double)o.partonFlavour); os << r;
  sprintf(r, "  %-32s: %f\n", "emEnergyFraction", (double)o.emEnergyFraction); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jet2_s
{
  double	energy;
  double	uncor_energy;
  double	pt;
  double	uncor_pt;
  double	phi;
  double	uncor_phi;
  double	eta;
  double	uncor_eta;
  double	jetID_fHPD;
  int	jetID_n90Hits;
  float	trackCountingHighEffBJetTags;
  float	trackCountingHighPurBJetTags;
  float	simpleSecondaryVertexBJetTags;
  float	simpleSecondaryVertexHighEffBJetTags;
  float	simpleSecondaryVertexHighPurBJetTags;
  float	combinedSecondaryVertexBJetTags;
  float	combinedSecondaryVertexMVABJetTags;
  float	jetBProbabilityBJetTags;
  float	jetProbabilityBJetTags;
  float	softElectronByIP3dBJetTags;
  float	softElectronByPtBJetTags;
  float	softMuonBJetTags;
  float	softMuonByIP3dBJetTags;
  float	softMuonByPtBJetTags;
  int	partonFlavour;
};
std::vector<jet2_s> jet2(500);

std::ostream& operator<<(std::ostream& os, const jet2_s& o)
{
  char r[1024];
  os << "jet2" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_energy", (double)o.uncor_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_pt", (double)o.uncor_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_phi", (double)o.uncor_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_eta", (double)o.uncor_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "jetID_fHPD", (double)o.jetID_fHPD); os << r;
  sprintf(r, "  %-32s: %f\n", "jetID_n90Hits", (double)o.jetID_n90Hits); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighEffBJetTags", (double)o.trackCountingHighEffBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighPurBJetTags", (double)o.trackCountingHighPurBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleSecondaryVertexBJetTags", (double)o.simpleSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleSecondaryVertexHighEffBJetTags", (double)o.simpleSecondaryVertexHighEffBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleSecondaryVertexHighPurBJetTags", (double)o.simpleSecondaryVertexHighPurBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexBJetTags", (double)o.combinedSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexMVABJetTags", (double)o.combinedSecondaryVertexMVABJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetBProbabilityBJetTags", (double)o.jetBProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetProbabilityBJetTags", (double)o.jetProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softElectronByIP3dBJetTags", (double)o.softElectronByIP3dBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softElectronByPtBJetTags", (double)o.softElectronByPtBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softMuonBJetTags", (double)o.softMuonBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softMuonByIP3dBJetTags", (double)o.softMuonByIP3dBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softMuonByPtBJetTags", (double)o.softMuonByPtBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "partonFlavour", (double)o.partonFlavour); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jethelper_s
{
  double	jetUncPlus;
  double	jetUncMinus;
  double	jecFactor;
  double	jecFactorNoL1Fast;
};
std::vector<jethelper_s> jethelper(500);

std::ostream& operator<<(std::ostream& os, const jethelper_s& o)
{
  char r[1024];
  os << "jethelper" << std::endl;
  sprintf(r, "  %-32s: %f\n", "jetUncPlus", (double)o.jetUncPlus); os << r;
  sprintf(r, "  %-32s: %f\n", "jetUncMinus", (double)o.jetUncMinus); os << r;
  sprintf(r, "  %-32s: %f\n", "jecFactor", (double)o.jecFactor); os << r;
  sprintf(r, "  %-32s: %f\n", "jecFactorNoL1Fast", (double)o.jecFactorNoL1Fast); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct met_s
{
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	sumEt;
  double	mEtSig;
};
std::vector<met_s> met(10);

std::ostream& operator<<(std::ostream& os, const met_s& o)
{
  char r[1024];
  os << "met" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct met2_s
{
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	sumEt;
  double	mEtSig;
};
std::vector<met2_s> met2(10);

std::ostream& operator<<(std::ostream& os, const met2_s& o)
{
  char r[1024];
  os << "met2" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct met3_s
{
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	sumEt;
  double	mEtSig;
};
std::vector<met3_s> met3(10);

std::ostream& operator<<(std::ostream& os, const met3_s& o)
{
  char r[1024];
  os << "met3" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct muon_s
{
  double	energy;
  double	pt;
  double	px;
  double	py;
  double	pz;
  double	phi;
  double	eta;
  float	trackIso;
  float	ecalIso;
  float	hcalIso;
  int	AllGlobalMuons;
  int	isGlobalMuon;
  int	isTrackerMuon;
  int	GlobalMuonPromptTight;
  double	track_normalizedChi2;
  unsigned short	track_hitPattern_numberOfValidMuonHits;
  unsigned short	track_hitPattern_numberOfValidPixelHits;
  unsigned short	innerTrack_numberOfValidHits;
  double	track_d0;
  float	innerTrack_pt;
  float	globalTrack_pt;
  double	dB;
  double	vz;
  float	chargedHadronIso;
  float	photonIso;
  float	neutralHadronIso;
  int	charge;
  double	combinedMuon_ndof;
  double	combinedMuon_chi2;
  unsigned short	combinedMuon_numberOfValidHits;
  double	genLepton_pt;
  double	genLepton_eta;
  double	genLepton_phi;
  int	genLepton_pdgId;
};
std::vector<muon_s> muon(50);

std::ostream& operator<<(std::ostream& os, const muon_s& o)
{
  char r[1024];
  os << "muon" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "px", (double)o.px); os << r;
  sprintf(r, "  %-32s: %f\n", "py", (double)o.py); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "trackIso", (double)o.trackIso); os << r;
  sprintf(r, "  %-32s: %f\n", "ecalIso", (double)o.ecalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "hcalIso", (double)o.hcalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "AllGlobalMuons", (double)o.AllGlobalMuons); os << r;
  sprintf(r, "  %-32s: %f\n", "isGlobalMuon", (double)o.isGlobalMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "isTrackerMuon", (double)o.isTrackerMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "GlobalMuonPromptTight", (double)o.GlobalMuonPromptTight); os << r;
  sprintf(r, "  %-32s: %f\n", "track_normalizedChi2", (double)o.track_normalizedChi2); os << r;
  sprintf(r, "  %-32s: %f\n", "track_hitPattern_numberOfValidMuonHits", (double)o.track_hitPattern_numberOfValidMuonHits); os << r;
  sprintf(r, "  %-32s: %f\n", "track_hitPattern_numberOfValidPixelHits", (double)o.track_hitPattern_numberOfValidPixelHits); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_numberOfValidHits", (double)o.innerTrack_numberOfValidHits); os << r;
  sprintf(r, "  %-32s: %f\n", "track_d0", (double)o.track_d0); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_pt", (double)o.innerTrack_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_pt", (double)o.globalTrack_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "dB", (double)o.dB); os << r;
  sprintf(r, "  %-32s: %f\n", "vz", (double)o.vz); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronIso", (double)o.chargedHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "photonIso", (double)o.photonIso); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronIso", (double)o.neutralHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedMuon_ndof", (double)o.combinedMuon_ndof); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedMuon_chi2", (double)o.combinedMuon_chi2); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedMuon_numberOfValidHits", (double)o.combinedMuon_numberOfValidHits); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_pt", (double)o.genLepton_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_eta", (double)o.genLepton_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_phi", (double)o.genLepton_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_pdgId", (double)o.genLepton_pdgId); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct muonhelper_s
{
  double	dxywrtBeamSpot;
  double	muonPFcandptdiff;
};
std::vector<muonhelper_s> muonhelper(50);

std::ostream& operator<<(std::ostream& os, const muonhelper_s& o)
{
  char r[1024];
  os << "muonhelper" << std::endl;
  sprintf(r, "  %-32s: %f\n", "dxywrtBeamSpot", (double)o.dxywrtBeamSpot); os << r;
  sprintf(r, "  %-32s: %f\n", "muonPFcandptdiff", (double)o.muonPFcandptdiff); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct tau_s
{
  double	energy;
  double	pt;
  double	px;
  double	py;
  double	pz;
  double	et;
  double	phi;
  double	eta;
  float	tauID_againstElectron;
  float	tauID_againstMuon;
  float	tauID_byIsolation;
  float	tauID_byTaNC;
  float	tauID_byTaNCfrHalfPercent;
  float	tauID_byTaNCfrQuarterPercent;
  float	trackIso;
  float	ecalIso;
  float	hcalIso;
  float	caloIso;
};
std::vector<tau_s> tau(50);

std::ostream& operator<<(std::ostream& os, const tau_s& o)
{
  char r[1024];
  os << "tau" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "px", (double)o.px); os << r;
  sprintf(r, "  %-32s: %f\n", "py", (double)o.py); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "tauID_againstElectron", (double)o.tauID_againstElectron); os << r;
  sprintf(r, "  %-32s: %f\n", "tauID_againstMuon", (double)o.tauID_againstMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "tauID_byIsolation", (double)o.tauID_byIsolation); os << r;
  sprintf(r, "  %-32s: %f\n", "tauID_byTaNC", (double)o.tauID_byTaNC); os << r;
  sprintf(r, "  %-32s: %f\n", "tauID_byTaNCfrHalfPercent", (double)o.tauID_byTaNCfrHalfPercent); os << r;
  sprintf(r, "  %-32s: %f\n", "tauID_byTaNCfrQuarterPercent", (double)o.tauID_byTaNCfrQuarterPercent); os << r;
  sprintf(r, "  %-32s: %f\n", "trackIso", (double)o.trackIso); os << r;
  sprintf(r, "  %-32s: %f\n", "ecalIso", (double)o.ecalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "hcalIso", (double)o.hcalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "caloIso", (double)o.caloIso); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct tauhelper_s
{
  int	genTauDecayModeID;
};
std::vector<tauhelper_s> tauhelper(50);

std::ostream& operator<<(std::ostream& os, const tauhelper_s& o)
{
  char r[1024];
  os << "tauhelper" << std::endl;
  sprintf(r, "  %-32s: %f\n", "genTauDecayModeID", (double)o.genTauDecayModeID); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct vertex_s
{
  int	isValid;
  int	isFake;
  double	ndof;
  double	x;
  double	y;
  double	z;
  double	position_Rho;
  double	normalizedChi2;
  int	tracksSize;
  double	xError;
  double	yError;
  double	zError;
};
std::vector<vertex_s> vertex(50);

std::ostream& operator<<(std::ostream& os, const vertex_s& o)
{
  char r[1024];
  os << "vertex" << std::endl;
  sprintf(r, "  %-32s: %f\n", "isValid", (double)o.isValid); os << r;
  sprintf(r, "  %-32s: %f\n", "isFake", (double)o.isFake); os << r;
  sprintf(r, "  %-32s: %f\n", "ndof", (double)o.ndof); os << r;
  sprintf(r, "  %-32s: %f\n", "x", (double)o.x); os << r;
  sprintf(r, "  %-32s: %f\n", "y", (double)o.y); os << r;
  sprintf(r, "  %-32s: %f\n", "z", (double)o.z); os << r;
  sprintf(r, "  %-32s: %f\n", "position_Rho", (double)o.position_Rho); os << r;
  sprintf(r, "  %-32s: %f\n", "normalizedChi2", (double)o.normalizedChi2); os << r;
  sprintf(r, "  %-32s: %f\n", "tracksSize", (double)o.tracksSize); os << r;
  sprintf(r, "  %-32s: %f\n", "xError", (double)o.xError); os << r;
  sprintf(r, "  %-32s: %f\n", "yError", (double)o.yError); os << r;
  sprintf(r, "  %-32s: %f\n", "zError", (double)o.zError); os << r;
  return os;
}
//-----------------------------------------------------------------------------

void fillObjects()
{

  electron.resize(electron_energy.size());
  for(unsigned int i=0; i < electron_energy.size(); ++i)
    {
      electron[i].energy	= electron_energy[i];
      electron[i].et	= electron_et[i];
      electron[i].pt	= electron_pt[i];
      electron[i].phi	= electron_phi[i];
      electron[i].eta	= electron_eta[i];
      electron[i].px	= electron_px[i];
      electron[i].py	= electron_py[i];
      electron[i].pz	= electron_pz[i];
      electron[i].trackIso	= electron_trackIso[i];
      electron[i].ecalIso	= electron_ecalIso[i];
      electron[i].hcalIso	= electron_hcalIso[i];
      electron[i].caloIso	= electron_caloIso[i];
      electron[i].eidRobustTight	= electron_eidRobustTight[i];
      electron[i].simpleEleId95relIso	= electron_simpleEleId95relIso[i];
      electron[i].gsfTrack_d0	= electron_gsfTrack_d0[i];
      electron[i].dr03TkSumPt	= electron_dr03TkSumPt[i];
      electron[i].dr03EcalRecHitSumEt	= electron_dr03EcalRecHitSumEt[i];
      electron[i].dr03HcalTowerSumEt	= electron_dr03HcalTowerSumEt[i];
      electron[i].vz	= electron_vz[i];
      electron[i].dB	= electron_dB[i];
      electron[i].chargedHadronIso	= electron_chargedHadronIso[i];
      electron[i].photonIso	= electron_photonIso[i];
      electron[i].neutralHadronIso	= electron_neutralHadronIso[i];
      electron[i].gsfTrack_trackerExpectedHitsInner_numberOfLostHits	= electron_gsfTrack_trackerExpectedHitsInner_numberOfLostHits[i];
      electron[i].superCluster_eta	= electron_superCluster_eta[i];
      electron[i].charge	= electron_charge[i];
      electron[i].genLepton_pt	= electron_genLepton_pt[i];
      electron[i].genLepton_eta	= electron_genLepton_eta[i];
      electron[i].genLepton_phi	= electron_genLepton_phi[i];
      electron[i].genLepton_pdgId	= electron_genLepton_pdgId[i];
      electron[i].hadronicOverEm	= electron_hadronicOverEm[i];
      electron[i].deltaPhiSuperClusterTrackAtVtx	= electron_deltaPhiSuperClusterTrackAtVtx[i];
      electron[i].deltaEtaSuperClusterTrackAtVtx	= electron_deltaEtaSuperClusterTrackAtVtx[i];
      electron[i].sigmaIetaIeta	= electron_sigmaIetaIeta[i];
    }

  electronhelper.resize(electronhelper_dxywrtBeamSpot.size());
  for(unsigned int i=0; i < electronhelper_dxywrtBeamSpot.size(); ++i)
    {
      electronhelper[i].dxywrtBeamSpot	= electronhelper_dxywrtBeamSpot[i];
    }

  geneventinfoproducthelper.resize(geneventinfoproducthelper_pdf1.size());
  for(unsigned int i=0; i < geneventinfoproducthelper_pdf1.size(); ++i)
    {
      geneventinfoproducthelper[i].pdf1	= geneventinfoproducthelper_pdf1[i];
      geneventinfoproducthelper[i].pdf2	= geneventinfoproducthelper_pdf2[i];
      geneventinfoproducthelper[i].pdfweight	= geneventinfoproducthelper_pdfweight[i];
      geneventinfoproducthelper[i].pdfweightsum	= geneventinfoproducthelper_pdfweightsum[i];
    }

  genparticlehelperra2b.resize(genparticlehelperra2b_firstMother.size());
  for(unsigned int i=0; i < genparticlehelperra2b_firstMother.size(); ++i)
    {
      genparticlehelperra2b[i].firstMother	= genparticlehelperra2b_firstMother[i];
      genparticlehelperra2b[i].lastMother	= genparticlehelperra2b_lastMother[i];
      genparticlehelperra2b[i].firstDaughter	= genparticlehelperra2b_firstDaughter[i];
      genparticlehelperra2b[i].lastDaughter	= genparticlehelperra2b_lastDaughter[i];
      genparticlehelperra2b[i].charge	= genparticlehelperra2b_charge[i];
      genparticlehelperra2b[i].pdgId	= genparticlehelperra2b_pdgId[i];
      genparticlehelperra2b[i].status	= genparticlehelperra2b_status[i];
      genparticlehelperra2b[i].pt	= genparticlehelperra2b_pt[i];
      genparticlehelperra2b[i].eta	= genparticlehelperra2b_eta[i];
      genparticlehelperra2b[i].phi	= genparticlehelperra2b_phi[i];
      genparticlehelperra2b[i].mass	= genparticlehelperra2b_mass[i];
    }

  jet.resize(jet_energy.size());
  for(unsigned int i=0; i < jet_energy.size(); ++i)
    {
      jet[i].energy	= jet_energy[i];
      jet[i].uncor_energy	= jet_uncor_energy[i];
      jet[i].pt	= jet_pt[i];
      jet[i].uncor_pt	= jet_uncor_pt[i];
      jet[i].phi	= jet_phi[i];
      jet[i].uncor_phi	= jet_uncor_phi[i];
      jet[i].eta	= jet_eta[i];
      jet[i].uncor_eta	= jet_uncor_eta[i];
      jet[i].jetArea	= jet_jetArea[i];
      jet[i].jetID_fHPD	= jet_jetID_fHPD[i];
      jet[i].jetID_n90Hits	= jet_jetID_n90Hits[i];
      jet[i].HFHadronEnergy	= jet_HFHadronEnergy[i];
      jet[i].neutralHadronEnergy	= jet_neutralHadronEnergy[i];
      jet[i].neutralEmEnergyFraction	= jet_neutralEmEnergyFraction[i];
      jet[i].neutralHadronEnergyFraction	= jet_neutralHadronEnergyFraction[i];
      jet[i].numberOfDaughters	= jet_numberOfDaughters[i];
      jet[i].chargedHadronEnergyFraction	= jet_chargedHadronEnergyFraction[i];
      jet[i].chargedMultiplicity	= jet_chargedMultiplicity[i];
      jet[i].chargedEmEnergyFraction	= jet_chargedEmEnergyFraction[i];
      jet[i].trackCountingHighEffBJetTags	= jet_trackCountingHighEffBJetTags[i];
      jet[i].trackCountingHighPurBJetTags	= jet_trackCountingHighPurBJetTags[i];
      jet[i].simpleSecondaryVertexBJetTags	= jet_simpleSecondaryVertexBJetTags[i];
      jet[i].simpleSecondaryVertexHighEffBJetTags	= jet_simpleSecondaryVertexHighEffBJetTags[i];
      jet[i].simpleSecondaryVertexHighPurBJetTags	= jet_simpleSecondaryVertexHighPurBJetTags[i];
      jet[i].combinedSecondaryVertexBJetTags	= jet_combinedSecondaryVertexBJetTags[i];
      jet[i].combinedSecondaryVertexMVABJetTags	= jet_combinedSecondaryVertexMVABJetTags[i];
      jet[i].jetBProbabilityBJetTags	= jet_jetBProbabilityBJetTags[i];
      jet[i].jetProbabilityBJetTags	= jet_jetProbabilityBJetTags[i];
      jet[i].softElectronByIP3dBJetTags	= jet_softElectronByIP3dBJetTags[i];
      jet[i].softElectronByPtBJetTags	= jet_softElectronByPtBJetTags[i];
      jet[i].softMuonBJetTags	= jet_softMuonBJetTags[i];
      jet[i].softMuonByIP3dBJetTags	= jet_softMuonByIP3dBJetTags[i];
      jet[i].softMuonByPtBJetTags	= jet_softMuonByPtBJetTags[i];
      jet[i].partonFlavour	= jet_partonFlavour[i];
      jet[i].emEnergyFraction	= jet_emEnergyFraction[i];
    }

  jet2.resize(jet2_energy.size());
  for(unsigned int i=0; i < jet2_energy.size(); ++i)
    {
      jet2[i].energy	= jet2_energy[i];
      jet2[i].uncor_energy	= jet2_uncor_energy[i];
      jet2[i].pt	= jet2_pt[i];
      jet2[i].uncor_pt	= jet2_uncor_pt[i];
      jet2[i].phi	= jet2_phi[i];
      jet2[i].uncor_phi	= jet2_uncor_phi[i];
      jet2[i].eta	= jet2_eta[i];
      jet2[i].uncor_eta	= jet2_uncor_eta[i];
      jet2[i].jetID_fHPD	= jet2_jetID_fHPD[i];
      jet2[i].jetID_n90Hits	= jet2_jetID_n90Hits[i];
      jet2[i].trackCountingHighEffBJetTags	= jet2_trackCountingHighEffBJetTags[i];
      jet2[i].trackCountingHighPurBJetTags	= jet2_trackCountingHighPurBJetTags[i];
      jet2[i].simpleSecondaryVertexBJetTags	= jet2_simpleSecondaryVertexBJetTags[i];
      jet2[i].simpleSecondaryVertexHighEffBJetTags	= jet2_simpleSecondaryVertexHighEffBJetTags[i];
      jet2[i].simpleSecondaryVertexHighPurBJetTags	= jet2_simpleSecondaryVertexHighPurBJetTags[i];
      jet2[i].combinedSecondaryVertexBJetTags	= jet2_combinedSecondaryVertexBJetTags[i];
      jet2[i].combinedSecondaryVertexMVABJetTags	= jet2_combinedSecondaryVertexMVABJetTags[i];
      jet2[i].jetBProbabilityBJetTags	= jet2_jetBProbabilityBJetTags[i];
      jet2[i].jetProbabilityBJetTags	= jet2_jetProbabilityBJetTags[i];
      jet2[i].softElectronByIP3dBJetTags	= jet2_softElectronByIP3dBJetTags[i];
      jet2[i].softElectronByPtBJetTags	= jet2_softElectronByPtBJetTags[i];
      jet2[i].softMuonBJetTags	= jet2_softMuonBJetTags[i];
      jet2[i].softMuonByIP3dBJetTags	= jet2_softMuonByIP3dBJetTags[i];
      jet2[i].softMuonByPtBJetTags	= jet2_softMuonByPtBJetTags[i];
      jet2[i].partonFlavour	= jet2_partonFlavour[i];
    }

  jethelper.resize(jethelper_jetUncPlus.size());
  for(unsigned int i=0; i < jethelper_jetUncPlus.size(); ++i)
    {
      jethelper[i].jetUncPlus	= jethelper_jetUncPlus[i];
      jethelper[i].jetUncMinus	= jethelper_jetUncMinus[i];
      jethelper[i].jecFactor	= jethelper_jecFactor[i];
      jethelper[i].jecFactorNoL1Fast	= jethelper_jecFactorNoL1Fast[i];
    }

  met.resize(met_energy.size());
  for(unsigned int i=0; i < met_energy.size(); ++i)
    {
      met[i].energy	= met_energy[i];
      met[i].et	= met_et[i];
      met[i].pt	= met_pt[i];
      met[i].phi	= met_phi[i];
      met[i].sumEt	= met_sumEt[i];
      met[i].mEtSig	= met_mEtSig[i];
    }

  met2.resize(met2_energy.size());
  for(unsigned int i=0; i < met2_energy.size(); ++i)
    {
      met2[i].energy	= met2_energy[i];
      met2[i].et	= met2_et[i];
      met2[i].pt	= met2_pt[i];
      met2[i].phi	= met2_phi[i];
      met2[i].sumEt	= met2_sumEt[i];
      met2[i].mEtSig	= met2_mEtSig[i];
    }

  met3.resize(met3_energy.size());
  for(unsigned int i=0; i < met3_energy.size(); ++i)
    {
      met3[i].energy	= met3_energy[i];
      met3[i].et	= met3_et[i];
      met3[i].pt	= met3_pt[i];
      met3[i].phi	= met3_phi[i];
      met3[i].sumEt	= met3_sumEt[i];
      met3[i].mEtSig	= met3_mEtSig[i];
    }

  muon.resize(muon_energy.size());
  for(unsigned int i=0; i < muon_energy.size(); ++i)
    {
      muon[i].energy	= muon_energy[i];
      muon[i].pt	= muon_pt[i];
      muon[i].px	= muon_px[i];
      muon[i].py	= muon_py[i];
      muon[i].pz	= muon_pz[i];
      muon[i].phi	= muon_phi[i];
      muon[i].eta	= muon_eta[i];
      muon[i].trackIso	= muon_trackIso[i];
      muon[i].ecalIso	= muon_ecalIso[i];
      muon[i].hcalIso	= muon_hcalIso[i];
      muon[i].AllGlobalMuons	= muon_AllGlobalMuons[i];
      muon[i].isGlobalMuon	= muon_isGlobalMuon[i];
      muon[i].isTrackerMuon	= muon_isTrackerMuon[i];
      muon[i].GlobalMuonPromptTight	= muon_GlobalMuonPromptTight[i];
      muon[i].track_normalizedChi2	= muon_track_normalizedChi2[i];
      muon[i].track_hitPattern_numberOfValidMuonHits	= muon_track_hitPattern_numberOfValidMuonHits[i];
      muon[i].track_hitPattern_numberOfValidPixelHits	= muon_track_hitPattern_numberOfValidPixelHits[i];
      muon[i].innerTrack_numberOfValidHits	= muon_innerTrack_numberOfValidHits[i];
      muon[i].track_d0	= muon_track_d0[i];
      muon[i].innerTrack_pt	= muon_innerTrack_pt[i];
      muon[i].globalTrack_pt	= muon_globalTrack_pt[i];
      muon[i].dB	= muon_dB[i];
      muon[i].vz	= muon_vz[i];
      muon[i].chargedHadronIso	= muon_chargedHadronIso[i];
      muon[i].photonIso	= muon_photonIso[i];
      muon[i].neutralHadronIso	= muon_neutralHadronIso[i];
      muon[i].charge	= muon_charge[i];
      muon[i].combinedMuon_ndof	= muon_combinedMuon_ndof[i];
      muon[i].combinedMuon_chi2	= muon_combinedMuon_chi2[i];
      muon[i].combinedMuon_numberOfValidHits	= muon_combinedMuon_numberOfValidHits[i];
      muon[i].genLepton_pt	= muon_genLepton_pt[i];
      muon[i].genLepton_eta	= muon_genLepton_eta[i];
      muon[i].genLepton_phi	= muon_genLepton_phi[i];
      muon[i].genLepton_pdgId	= muon_genLepton_pdgId[i];
    }

  muonhelper.resize(muonhelper_dxywrtBeamSpot.size());
  for(unsigned int i=0; i < muonhelper_dxywrtBeamSpot.size(); ++i)
    {
      muonhelper[i].dxywrtBeamSpot	= muonhelper_dxywrtBeamSpot[i];
      muonhelper[i].muonPFcandptdiff	= muonhelper_muonPFcandptdiff[i];
    }

  tau.resize(tau_energy.size());
  for(unsigned int i=0; i < tau_energy.size(); ++i)
    {
      tau[i].energy	= tau_energy[i];
      tau[i].pt	= tau_pt[i];
      tau[i].px	= tau_px[i];
      tau[i].py	= tau_py[i];
      tau[i].pz	= tau_pz[i];
      tau[i].et	= tau_et[i];
      tau[i].phi	= tau_phi[i];
      tau[i].eta	= tau_eta[i];
      tau[i].tauID_againstElectron	= tau_tauID_againstElectron[i];
      tau[i].tauID_againstMuon	= tau_tauID_againstMuon[i];
      tau[i].tauID_byIsolation	= tau_tauID_byIsolation[i];
      tau[i].tauID_byTaNC	= tau_tauID_byTaNC[i];
      tau[i].tauID_byTaNCfrHalfPercent	= tau_tauID_byTaNCfrHalfPercent[i];
      tau[i].tauID_byTaNCfrQuarterPercent	= tau_tauID_byTaNCfrQuarterPercent[i];
      tau[i].trackIso	= tau_trackIso[i];
      tau[i].ecalIso	= tau_ecalIso[i];
      tau[i].hcalIso	= tau_hcalIso[i];
      tau[i].caloIso	= tau_caloIso[i];
    }

  tauhelper.resize(tauhelper_genTauDecayModeID.size());
  for(unsigned int i=0; i < tauhelper_genTauDecayModeID.size(); ++i)
    {
      tauhelper[i].genTauDecayModeID	= tauhelper_genTauDecayModeID[i];
    }

  vertex.resize(vertex_isValid.size());
  for(unsigned int i=0; i < vertex_isValid.size(); ++i)
    {
      vertex[i].isValid	= vertex_isValid[i];
      vertex[i].isFake	= vertex_isFake[i];
      vertex[i].ndof	= vertex_ndof[i];
      vertex[i].x	= vertex_x[i];
      vertex[i].y	= vertex_y[i];
      vertex[i].z	= vertex_z[i];
      vertex[i].position_Rho	= vertex_position_Rho[i];
      vertex[i].normalizedChi2	= vertex_normalizedChi2[i];
      vertex[i].tracksSize	= vertex_tracksSize[i];
      vertex[i].xError	= vertex_xError[i];
      vertex[i].yError	= vertex_yError[i];
      vertex[i].zError	= vertex_zError[i];
    }
}


#endif
