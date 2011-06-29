#ifndef BASICLOOPCU_H
#define BASICLOOPCU_H
//-----------------------------------------------------------------------------
// File:        BasicLoopCU.h
// Description: Analyzer header for ntuples created by TheNtupleMaker
// Created:     Tue Jun 28 14:01:46 2011 by mkntanalyzer.py
// Author:      Ben Kreis
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
double	beamspot_x0;
double	beamspot_y0;
double	beamspot_z0;
double	edmevent_bunchCrossing;
double	edmevent_event;
double	edmevent_isRealData;
double	edmevent_luminosityBlock;
double	edmevent_run;
double	edmtriggerresults_HLT_HT100U;
double	edmtriggerresults_HLT_HT100U_prs;
double	edmtriggerresults_HLT_HT120U;
double	edmtriggerresults_HLT_HT120U_prs;
double	edmtriggerresults_HLT_HT140U;
double	edmtriggerresults_HLT_HT140U_prs;
double	edmtriggerresults_HLT_HT150U_v3;
double	edmtriggerresults_HLT_HT150U_v3_prs;
double	edmtriggerresults_HLT_HT150_v1;
double	edmtriggerresults_HLT_HT150_v1_prs;
double	edmtriggerresults_HLT_HT150_v2;
double	edmtriggerresults_HLT_HT150_v2_prs;
double	edmtriggerresults_HLT_HT150_v3;
double	edmtriggerresults_HLT_HT150_v3_prs;
double	edmtriggerresults_HLT_HT150_v4;
double	edmtriggerresults_HLT_HT150_v4_prs;
double	edmtriggerresults_HLT_HT150_v5;
double	edmtriggerresults_HLT_HT150_v5_prs;
double	edmtriggerresults_HLT_HT250_MHT60_v2;
double	edmtriggerresults_HLT_HT250_MHT60_v2_prs;
double	edmtriggerresults_HLT_HT250_MHT60_v3;
double	edmtriggerresults_HLT_HT250_MHT60_v3_prs;
double	edmtriggerresults_HLT_HT250_v1;
double	edmtriggerresults_HLT_HT250_v1_prs;
double	edmtriggerresults_HLT_HT250_v2;
double	edmtriggerresults_HLT_HT250_v2_prs;
double	edmtriggerresults_HLT_HT250_v3;
double	edmtriggerresults_HLT_HT250_v3_prs;
double	edmtriggerresults_HLT_HT250_v4;
double	edmtriggerresults_HLT_HT250_v4_prs;
double	edmtriggerresults_HLT_HT250_v5;
double	edmtriggerresults_HLT_HT250_v5_prs;
double	edmtriggerresults_HLT_HT260_MHT60_v2;
double	edmtriggerresults_HLT_HT260_MHT60_v2_prs;
double	edmtriggerresults_HLT_HT260_v2;
double	edmtriggerresults_HLT_HT260_v2_prs;
double	edmtriggerresults_HLT_HT260_v3;
double	edmtriggerresults_HLT_HT260_v3_prs;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v2;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v2_prs;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3_prs;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v4;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v4_prs;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v5;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v5_prs;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2_prs;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3_prs;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4_prs;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5_prs;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v2;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v2_prs;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v3;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v3_prs;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v4;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v4_prs;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v5;
double	edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v5_prs;
double	edmtriggerresults_HLT_HT300_MHT75_v2;
double	edmtriggerresults_HLT_HT300_MHT75_v2_prs;
double	edmtriggerresults_HLT_HT300_MHT75_v3;
double	edmtriggerresults_HLT_HT300_MHT75_v3_prs;
double	edmtriggerresults_HLT_HT300_MHT75_v4;
double	edmtriggerresults_HLT_HT300_MHT75_v4_prs;
double	edmtriggerresults_HLT_HT300_PFMHT55_v2;
double	edmtriggerresults_HLT_HT300_PFMHT55_v2_prs;
double	edmtriggerresults_HLT_HT300_PFMHT55_v3;
double	edmtriggerresults_HLT_HT300_PFMHT55_v3_prs;
double	edmtriggerresults_HLT_HT300_PFMHT55_v4;
double	edmtriggerresults_HLT_HT300_PFMHT55_v4_prs;
double	edmtriggerresults_HLT_HT300_PFMHT55_v5;
double	edmtriggerresults_HLT_HT300_PFMHT55_v5_prs;
double	edmtriggerresults_HLT_HT300_v1;
double	edmtriggerresults_HLT_HT300_v1_prs;
double	edmtriggerresults_HLT_HT300_v2;
double	edmtriggerresults_HLT_HT300_v2_prs;
double	edmtriggerresults_HLT_HT300_v3;
double	edmtriggerresults_HLT_HT300_v3_prs;
double	edmtriggerresults_HLT_HT300_v4;
double	edmtriggerresults_HLT_HT300_v4_prs;
double	edmtriggerresults_HLT_HT300_v5;
double	edmtriggerresults_HLT_HT300_v5_prs;
double	edmtriggerresults_HLT_HT300_v6;
double	edmtriggerresults_HLT_HT300_v6_prs;
double	edmtriggerresults_HLT_HT300_v7;
double	edmtriggerresults_HLT_HT300_v7_prs;
std::vector<double>	electron1_caloIso(7,0);
std::vector<double>	electron1_charge(7,0);
std::vector<double>	electron1_chargedHadronIso(7,0);
std::vector<double>	electron1_dB(7,0);
std::vector<double>	electron1_deltaEtaSuperClusterTrackAtVtx(7,0);
std::vector<double>	electron1_deltaPhiSuperClusterTrackAtVtx(7,0);
std::vector<double>	electron1_dr03EcalRecHitSumEt(7,0);
std::vector<double>	electron1_dr03HcalTowerSumEt(7,0);
std::vector<double>	electron1_dr03TkSumPt(7,0);
std::vector<double>	electron1_ecalIso(7,0);
std::vector<double>	electron1_eidRobustTight(7,0);
std::vector<double>	electron1_energy(7,0);
std::vector<double>	electron1_et(7,0);
std::vector<double>	electron1_eta(7,0);
std::vector<double>	electron1_genLepton_eta(7,0);
std::vector<double>	electron1_genLepton_pdgId(7,0);
std::vector<double>	electron1_genLepton_phi(7,0);
std::vector<double>	electron1_genLepton_pt(7,0);
std::vector<double>	electron1_gsfTrack_d0(7,0);
std::vector<double>	electron1_gsfTrack_trackerExpectedHitsInner_numberOfLostHits(7,0);
std::vector<double>	electron1_hadronicOverEm(7,0);
std::vector<double>	electron1_hcalIso(7,0);
std::vector<double>	electron1_neutralHadronIso(7,0);
std::vector<double>	electron1_phi(7,0);
std::vector<double>	electron1_photonIso(7,0);
std::vector<double>	electron1_pt(7,0);
std::vector<double>	electron1_px(7,0);
std::vector<double>	electron1_py(7,0);
std::vector<double>	electron1_pz(7,0);
std::vector<double>	electron1_sigmaIetaIeta(7,0);
std::vector<double>	electron1_simpleEleId95relIso(7,0);
std::vector<double>	electron1_superCluster_eta(7,0);
std::vector<double>	electron1_trackIso(7,0);
std::vector<double>	electron1_vz(7,0);
std::vector<double>	electron_dxywrtBeamSpot(7,0);
double	geneventinfoproduct1_pdf1;
double	geneventinfoproduct1_pdf2;
double	geneventinfoproduct1_scalePDF;
double	geneventinfoproduct1_weight;
double	geneventinfoproduct1_x1;
double	geneventinfoproduct1_x2;
double	geneventinfoproduct_pdf1;
double	geneventinfoproduct_pdf2;
double	geneventinfoproduct_pdfweight;
double	geneventinfoproduct_pdfweightsum;
double	genparticlera2_charge;
double	genparticlera2_eta;
double	genparticlera2_firstDaughter;
double	genparticlera2_firstMother;
double	genparticlera2_lastDaughter;
double	genparticlera2_lastMother;
double	genparticlera2_mass;
double	genparticlera2_pdgId;
double	genparticlera2_phi;
double	genparticlera2_pt;
double	genparticlera2_status;
double	genruninfoproduct_externalXSecLO_error;
double	genruninfoproduct_externalXSecLO_value;
double	genruninfoproduct_filterEfficiency;
double	genruninfoproduct_internalXSec_error;
double	genruninfoproduct_internalXSec_value;
std::vector<double>	jet1_combinedSecondaryVertexBJetTags(177,0);
std::vector<double>	jet1_combinedSecondaryVertexMVABJetTags(177,0);
std::vector<double>	jet1_emEnergyFraction(177,0);
std::vector<double>	jet1_energy(177,0);
std::vector<double>	jet1_eta(177,0);
std::vector<double>	jet1_jetBProbabilityBJetTags(177,0);
std::vector<double>	jet1_jetID_fHPD(177,0);
std::vector<double>	jet1_jetID_n90Hits(177,0);
std::vector<double>	jet1_jetProbabilityBJetTags(177,0);
std::vector<double>	jet1_partonFlavour(177,0);
std::vector<double>	jet1_phi(177,0);
std::vector<double>	jet1_pt(177,0);
std::vector<double>	jet1_simpleSecondaryVertexBJetTags(177,0);
std::vector<double>	jet1_simpleSecondaryVertexHighEffBJetTags(177,0);
std::vector<double>	jet1_simpleSecondaryVertexHighPurBJetTags(177,0);
std::vector<double>	jet1_softElectronByIP3dBJetTags(177,0);
std::vector<double>	jet1_softElectronByPtBJetTags(177,0);
std::vector<double>	jet1_softMuonBJetTags(177,0);
std::vector<double>	jet1_softMuonByIP3dBJetTags(177,0);
std::vector<double>	jet1_softMuonByPtBJetTags(177,0);
std::vector<double>	jet1_trackCountingHighEffBJetTags(177,0);
std::vector<double>	jet1_trackCountingHighPurBJetTags(177,0);
std::vector<double>	jet1_uncor_energy(177,0);
std::vector<double>	jet1_uncor_eta(177,0);
std::vector<double>	jet1_uncor_phi(177,0);
std::vector<double>	jet1_uncor_pt(177,0);
std::vector<double>	jet2_HFHadronEnergy(173,0);
std::vector<double>	jet2_chargedEmEnergyFraction(173,0);
std::vector<double>	jet2_chargedHadronEnergyFraction(173,0);
std::vector<double>	jet2_chargedMultiplicity(173,0);
std::vector<double>	jet2_combinedSecondaryVertexBJetTags(173,0);
std::vector<double>	jet2_combinedSecondaryVertexMVABJetTags(173,0);
std::vector<double>	jet2_energy(173,0);
std::vector<double>	jet2_eta(173,0);
std::vector<double>	jet2_jetArea(173,0);
std::vector<double>	jet2_jetBProbabilityBJetTags(173,0);
std::vector<double>	jet2_jetID_fHPD(173,0);
std::vector<double>	jet2_jetID_n90Hits(173,0);
std::vector<double>	jet2_jetProbabilityBJetTags(173,0);
std::vector<double>	jet2_neutralEmEnergyFraction(173,0);
std::vector<double>	jet2_neutralHadronEnergy(173,0);
std::vector<double>	jet2_neutralHadronEnergyFraction(173,0);
std::vector<double>	jet2_numberOfDaughters(173,0);
std::vector<double>	jet2_partonFlavour(173,0);
std::vector<double>	jet2_phi(173,0);
std::vector<double>	jet2_pt(173,0);
std::vector<double>	jet2_simpleSecondaryVertexBJetTags(173,0);
std::vector<double>	jet2_simpleSecondaryVertexHighEffBJetTags(173,0);
std::vector<double>	jet2_simpleSecondaryVertexHighPurBJetTags(173,0);
std::vector<double>	jet2_softElectronByIP3dBJetTags(173,0);
std::vector<double>	jet2_softElectronByPtBJetTags(173,0);
std::vector<double>	jet2_softMuonBJetTags(173,0);
std::vector<double>	jet2_softMuonByIP3dBJetTags(173,0);
std::vector<double>	jet2_softMuonByPtBJetTags(173,0);
std::vector<double>	jet2_trackCountingHighEffBJetTags(173,0);
std::vector<double>	jet2_trackCountingHighPurBJetTags(173,0);
std::vector<double>	jet2_uncor_energy(173,0);
std::vector<double>	jet2_uncor_eta(173,0);
std::vector<double>	jet2_uncor_phi(173,0);
std::vector<double>	jet2_uncor_pt(173,0);
std::vector<double>	jet_jecFactor(173,0);
std::vector<double>	jet_jecFactorNoL1Fast(173,0);
std::vector<double>	jet_jetUncMinus(173,0);
std::vector<double>	jet_jetUncPlus(173,0);
std::vector<double>	met1_energy(3,0);
std::vector<double>	met1_et(3,0);
std::vector<double>	met1_mEtSig(3,0);
std::vector<double>	met1_phi(3,0);
std::vector<double>	met1_pt(3,0);
std::vector<double>	met1_sumEt(3,0);
std::vector<double>	met2_energy(3,0);
std::vector<double>	met2_et(3,0);
std::vector<double>	met2_mEtSig(3,0);
std::vector<double>	met2_phi(3,0);
std::vector<double>	met2_pt(3,0);
std::vector<double>	met2_sumEt(3,0);
std::vector<double>	met_energy(3,0);
std::vector<double>	met_et(3,0);
std::vector<double>	met_mEtSig(3,0);
std::vector<double>	met_phi(3,0);
std::vector<double>	met_pt(3,0);
std::vector<double>	met_sumEt(3,0);
std::vector<double>	muon1_AllGlobalMuons(5,0);
std::vector<double>	muon1_GlobalMuonPromptTight(5,0);
std::vector<double>	muon1_charge(5,0);
std::vector<double>	muon1_chargedHadronIso(5,0);
std::vector<double>	muon1_combinedMuon_chi2(5,0);
std::vector<double>	muon1_combinedMuon_ndof(5,0);
std::vector<double>	muon1_combinedMuon_numberOfValidHits(5,0);
std::vector<double>	muon1_dB(5,0);
std::vector<double>	muon1_ecalIso(5,0);
std::vector<double>	muon1_energy(5,0);
std::vector<double>	muon1_eta(5,0);
std::vector<double>	muon1_genLepton_eta(5,0);
std::vector<double>	muon1_genLepton_pdgId(5,0);
std::vector<double>	muon1_genLepton_phi(5,0);
std::vector<double>	muon1_genLepton_pt(5,0);
std::vector<double>	muon1_globalTrack_pt(5,0);
std::vector<double>	muon1_hcalIso(5,0);
std::vector<double>	muon1_innerTrack_numberOfValidHits(5,0);
std::vector<double>	muon1_innerTrack_pt(5,0);
std::vector<double>	muon1_isGlobalMuon(5,0);
std::vector<double>	muon1_isTrackerMuon(5,0);
std::vector<double>	muon1_neutralHadronIso(5,0);
std::vector<double>	muon1_phi(5,0);
std::vector<double>	muon1_photonIso(5,0);
std::vector<double>	muon1_pt(5,0);
std::vector<double>	muon1_px(5,0);
std::vector<double>	muon1_py(5,0);
std::vector<double>	muon1_pz(5,0);
std::vector<double>	muon1_trackIso(5,0);
std::vector<double>	muon1_track_d0(5,0);
std::vector<double>	muon1_track_hitPattern_numberOfValidMuonHits(5,0);
std::vector<double>	muon1_track_hitPattern_numberOfValidPixelHits(5,0);
std::vector<double>	muon1_track_normalizedChi2(5,0);
std::vector<double>	muon1_vz(5,0);
std::vector<double>	muon_dxywrtBeamSpot(5,0);
std::vector<double>	muon_muonPFcandptdiff(5,0);
int	nGenEventInfoProductHelper_generator;
int	nelectron;
int	nelectron1;
int	ngenparticlera2b;
int	njet;
int	njet1;
int	njet2;
int	nmet;
int	nmet1;
int	nmet2;
int	nmuon;
int	nmuon1;
int	ntau;
int	ntau1;
int	nvertex;
double	sdouble_kt6pfjets_rho_value;
double	sdouble_kt6pfjets_sigma_value;
std::vector<double>	tau1_caloIso(9,0);
std::vector<double>	tau1_ecalIso(9,0);
std::vector<double>	tau1_energy(9,0);
std::vector<double>	tau1_et(9,0);
std::vector<double>	tau1_eta(9,0);
std::vector<double>	tau1_hcalIso(9,0);
std::vector<double>	tau1_phi(9,0);
std::vector<double>	tau1_pt(9,0);
std::vector<double>	tau1_px(9,0);
std::vector<double>	tau1_py(9,0);
std::vector<double>	tau1_pz(9,0);
std::vector<double>	tau1_tauID_againstElectron(9,0);
std::vector<double>	tau1_tauID_againstMuon(9,0);
std::vector<double>	tau1_tauID_byIsolation(9,0);
std::vector<double>	tau1_tauID_byTaNC(9,0);
std::vector<double>	tau1_tauID_byTaNCfrHalfPercent(9,0);
std::vector<double>	tau1_tauID_byTaNCfrQuarterPercent(9,0);
std::vector<double>	tau1_trackIso(9,0);
std::vector<double>	tau_genTauDecayModeID(9,0);
std::vector<double>	vertex_isFake(35,0);
std::vector<double>	vertex_isValid(35,0);
std::vector<double>	vertex_ndof(35,0);
std::vector<double>	vertex_normalizedChi2(35,0);
std::vector<double>	vertex_position_Rho(35,0);
std::vector<double>	vertex_tracksSize(35,0);
std::vector<double>	vertex_x(35,0);
std::vector<double>	vertex_xError(35,0);
std::vector<double>	vertex_y(35,0);
std::vector<double>	vertex_yError(35,0);
std::vector<double>	vertex_z(35,0);
std::vector<double>	vertex_zError(35,0);

//-----------------------------------------------------------------------------
// -- Select variables to be read
//-----------------------------------------------------------------------------
void selectVariables(itreestream& stream)
{
  stream.select("recoBeamSpot_offlineBeamSpot.x0", beamspot_x0);
  stream.select("recoBeamSpot_offlineBeamSpot.y0", beamspot_y0);
  stream.select("recoBeamSpot_offlineBeamSpot.z0", beamspot_z0);
  stream.select("edmEventHelper_info.bunchCrossing", edmevent_bunchCrossing);
  stream.select("edmEventHelper_info.event", edmevent_event);
  stream.select("edmEventHelper_info.isRealData", edmevent_isRealData);
  stream.select("edmEventHelper_info.luminosityBlock", edmevent_luminosityBlock);
  stream.select("edmEventHelper_info.run", edmevent_run);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT100U", edmtriggerresults_HLT_HT100U);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT100U_prs", edmtriggerresults_HLT_HT100U_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT120U", edmtriggerresults_HLT_HT120U);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT120U_prs", edmtriggerresults_HLT_HT120U_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT140U", edmtriggerresults_HLT_HT140U);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT140U_prs", edmtriggerresults_HLT_HT140U_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150U_v3", edmtriggerresults_HLT_HT150U_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150U_v3_prs", edmtriggerresults_HLT_HT150U_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v1", edmtriggerresults_HLT_HT150_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v1_prs", edmtriggerresults_HLT_HT150_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v2", edmtriggerresults_HLT_HT150_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v2_prs", edmtriggerresults_HLT_HT150_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v3", edmtriggerresults_HLT_HT150_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v3_prs", edmtriggerresults_HLT_HT150_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v4", edmtriggerresults_HLT_HT150_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v4_prs", edmtriggerresults_HLT_HT150_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v5", edmtriggerresults_HLT_HT150_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v5_prs", edmtriggerresults_HLT_HT150_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT60_v2", edmtriggerresults_HLT_HT250_MHT60_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT60_v2_prs", edmtriggerresults_HLT_HT250_MHT60_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT60_v3", edmtriggerresults_HLT_HT250_MHT60_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT60_v3_prs", edmtriggerresults_HLT_HT250_MHT60_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v1", edmtriggerresults_HLT_HT250_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v1_prs", edmtriggerresults_HLT_HT250_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v2", edmtriggerresults_HLT_HT250_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v2_prs", edmtriggerresults_HLT_HT250_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v3", edmtriggerresults_HLT_HT250_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v3_prs", edmtriggerresults_HLT_HT250_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v4", edmtriggerresults_HLT_HT250_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v4_prs", edmtriggerresults_HLT_HT250_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v5", edmtriggerresults_HLT_HT250_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v5_prs", edmtriggerresults_HLT_HT250_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT260_MHT60_v2", edmtriggerresults_HLT_HT260_MHT60_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT260_MHT60_v2_prs", edmtriggerresults_HLT_HT260_MHT60_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT260_v2", edmtriggerresults_HLT_HT260_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT260_v2_prs", edmtriggerresults_HLT_HT260_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT260_v3", edmtriggerresults_HLT_HT260_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT260_v3_prs", edmtriggerresults_HLT_HT260_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v2", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v2_prs", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3_prs", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v4", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v4_prs", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v5", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v5_prs", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2_prs", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3_prs", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4_prs", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5_prs", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v2", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v2_prs", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v3", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v3_prs", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v4", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v4_prs", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v5", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v5_prs", edmtriggerresults_HLT_HT300_CentralJet30_BTagIP_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v2", edmtriggerresults_HLT_HT300_MHT75_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v2_prs", edmtriggerresults_HLT_HT300_MHT75_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v3", edmtriggerresults_HLT_HT300_MHT75_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v3_prs", edmtriggerresults_HLT_HT300_MHT75_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v4", edmtriggerresults_HLT_HT300_MHT75_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v4_prs", edmtriggerresults_HLT_HT300_MHT75_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v2", edmtriggerresults_HLT_HT300_PFMHT55_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v2_prs", edmtriggerresults_HLT_HT300_PFMHT55_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v3", edmtriggerresults_HLT_HT300_PFMHT55_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v3_prs", edmtriggerresults_HLT_HT300_PFMHT55_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v4", edmtriggerresults_HLT_HT300_PFMHT55_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v4_prs", edmtriggerresults_HLT_HT300_PFMHT55_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v5", edmtriggerresults_HLT_HT300_PFMHT55_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v5_prs", edmtriggerresults_HLT_HT300_PFMHT55_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v1", edmtriggerresults_HLT_HT300_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v1_prs", edmtriggerresults_HLT_HT300_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v2", edmtriggerresults_HLT_HT300_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v2_prs", edmtriggerresults_HLT_HT300_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v3", edmtriggerresults_HLT_HT300_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v3_prs", edmtriggerresults_HLT_HT300_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v4", edmtriggerresults_HLT_HT300_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v4_prs", edmtriggerresults_HLT_HT300_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v5", edmtriggerresults_HLT_HT300_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v5_prs", edmtriggerresults_HLT_HT300_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v6", edmtriggerresults_HLT_HT300_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v6_prs", edmtriggerresults_HLT_HT300_v6_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v7", edmtriggerresults_HLT_HT300_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v7_prs", edmtriggerresults_HLT_HT300_v7_prs);
  stream.select("patElectron_selectedPatElectronsPF.caloIso", electron1_caloIso);
  stream.select("patElectron_selectedPatElectronsPF.charge", electron1_charge);
  stream.select("patElectron_selectedPatElectronsPF.chargedHadronIso", electron1_chargedHadronIso);
  stream.select("patElectron_selectedPatElectronsPF.dB", electron1_dB);
  stream.select("patElectron_selectedPatElectronsPF.deltaEtaSuperClusterTrackAtVtx", electron1_deltaEtaSuperClusterTrackAtVtx);
  stream.select("patElectron_selectedPatElectronsPF.deltaPhiSuperClusterTrackAtVtx", electron1_deltaPhiSuperClusterTrackAtVtx);
  stream.select("patElectron_selectedPatElectronsPF.dr03EcalRecHitSumEt", electron1_dr03EcalRecHitSumEt);
  stream.select("patElectron_selectedPatElectronsPF.dr03HcalTowerSumEt", electron1_dr03HcalTowerSumEt);
  stream.select("patElectron_selectedPatElectronsPF.dr03TkSumPt", electron1_dr03TkSumPt);
  stream.select("patElectron_selectedPatElectronsPF.ecalIso", electron1_ecalIso);
  stream.select("patElectron_selectedPatElectronsPF.eidRobustTight", electron1_eidRobustTight);
  stream.select("patElectron_selectedPatElectronsPF.energy", electron1_energy);
  stream.select("patElectron_selectedPatElectronsPF.et", electron1_et);
  stream.select("patElectron_selectedPatElectronsPF.eta", electron1_eta);
  stream.select("patElectron_selectedPatElectronsPF.genLepton_eta", electron1_genLepton_eta);
  stream.select("patElectron_selectedPatElectronsPF.genLepton_pdgId", electron1_genLepton_pdgId);
  stream.select("patElectron_selectedPatElectronsPF.genLepton_phi", electron1_genLepton_phi);
  stream.select("patElectron_selectedPatElectronsPF.genLepton_pt", electron1_genLepton_pt);
  stream.select("patElectron_selectedPatElectronsPF.gsfTrack_d0", electron1_gsfTrack_d0);
  stream.select("patElectron_selectedPatElectronsPF.gsfTrack_trackerExpectedHitsInner_numberOfLostHits", electron1_gsfTrack_trackerExpectedHitsInner_numberOfLostHits);
  stream.select("patElectron_selectedPatElectronsPF.hadronicOverEm", electron1_hadronicOverEm);
  stream.select("patElectron_selectedPatElectronsPF.hcalIso", electron1_hcalIso);
  stream.select("patElectron_selectedPatElectronsPF.neutralHadronIso", electron1_neutralHadronIso);
  stream.select("patElectron_selectedPatElectronsPF.phi", electron1_phi);
  stream.select("patElectron_selectedPatElectronsPF.photonIso", electron1_photonIso);
  stream.select("patElectron_selectedPatElectronsPF.pt", electron1_pt);
  stream.select("patElectron_selectedPatElectronsPF.px", electron1_px);
  stream.select("patElectron_selectedPatElectronsPF.py", electron1_py);
  stream.select("patElectron_selectedPatElectronsPF.pz", electron1_pz);
  stream.select("patElectron_selectedPatElectronsPF.sigmaIetaIeta", electron1_sigmaIetaIeta);
  stream.select("patElectron_selectedPatElectronsPF.simpleEleId95relIso", electron1_simpleEleId95relIso);
  stream.select("patElectron_selectedPatElectronsPF.superCluster_eta", electron1_superCluster_eta);
  stream.select("patElectron_selectedPatElectronsPF.trackIso", electron1_trackIso);
  stream.select("patElectron_selectedPatElectronsPF.vz", electron1_vz);
  stream.select("patElectronHelper_selectedPatElectronsPF.dxywrtBeamSpot", electron_dxywrtBeamSpot);
  stream.select("GenEventInfoProduct_generator.pdf1", geneventinfoproduct1_pdf1);
  stream.select("GenEventInfoProduct_generator.pdf2", geneventinfoproduct1_pdf2);
  stream.select("GenEventInfoProduct_generator.scalePDF", geneventinfoproduct1_scalePDF);
  stream.select("GenEventInfoProduct_generator.weight", geneventinfoproduct1_weight);
  stream.select("GenEventInfoProduct_generator.x1", geneventinfoproduct1_x1);
  stream.select("GenEventInfoProduct_generator.x2", geneventinfoproduct1_x2);
  stream.select("GenEventInfoProductHelper_generator.pdf1", geneventinfoproduct_pdf1);
  stream.select("GenEventInfoProductHelper_generator.pdf2", geneventinfoproduct_pdf2);
  stream.select("GenEventInfoProductHelper_generator.pdfweight", geneventinfoproduct_pdfweight);
  stream.select("GenEventInfoProductHelper_generator.pdfweightsum", geneventinfoproduct_pdfweightsum);
  stream.select("recoGenParticleHelperRA2b_genParticles.charge", genparticlera2_charge);
  stream.select("recoGenParticleHelperRA2b_genParticles.eta", genparticlera2_eta);
  stream.select("recoGenParticleHelperRA2b_genParticles.firstDaughter", genparticlera2_firstDaughter);
  stream.select("recoGenParticleHelperRA2b_genParticles.firstMother", genparticlera2_firstMother);
  stream.select("recoGenParticleHelperRA2b_genParticles.lastDaughter", genparticlera2_lastDaughter);
  stream.select("recoGenParticleHelperRA2b_genParticles.lastMother", genparticlera2_lastMother);
  stream.select("recoGenParticleHelperRA2b_genParticles.mass", genparticlera2_mass);
  stream.select("recoGenParticleHelperRA2b_genParticles.pdgId", genparticlera2_pdgId);
  stream.select("recoGenParticleHelperRA2b_genParticles.phi", genparticlera2_phi);
  stream.select("recoGenParticleHelperRA2b_genParticles.pt", genparticlera2_pt);
  stream.select("recoGenParticleHelperRA2b_genParticles.status", genparticlera2_status);
  stream.select("GenRunInfoProduct_generator.externalXSecLO_error", genruninfoproduct_externalXSecLO_error);
  stream.select("GenRunInfoProduct_generator.externalXSecLO_value", genruninfoproduct_externalXSecLO_value);
  stream.select("GenRunInfoProduct_generator.filterEfficiency", genruninfoproduct_filterEfficiency);
  stream.select("GenRunInfoProduct_generator.internalXSec_error", genruninfoproduct_internalXSec_error);
  stream.select("GenRunInfoProduct_generator.internalXSec_value", genruninfoproduct_internalXSec_value);
  stream.select("patJet_cleanPatJetsAK5Calo.combinedSecondaryVertexBJetTags", jet1_combinedSecondaryVertexBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.combinedSecondaryVertexMVABJetTags", jet1_combinedSecondaryVertexMVABJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.emEnergyFraction", jet1_emEnergyFraction);
  stream.select("patJet_cleanPatJetsAK5Calo.energy", jet1_energy);
  stream.select("patJet_cleanPatJetsAK5Calo.eta", jet1_eta);
  stream.select("patJet_cleanPatJetsAK5Calo.jetBProbabilityBJetTags", jet1_jetBProbabilityBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.jetID_fHPD", jet1_jetID_fHPD);
  stream.select("patJet_cleanPatJetsAK5Calo.jetID_n90Hits", jet1_jetID_n90Hits);
  stream.select("patJet_cleanPatJetsAK5Calo.jetProbabilityBJetTags", jet1_jetProbabilityBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.partonFlavour", jet1_partonFlavour);
  stream.select("patJet_cleanPatJetsAK5Calo.phi", jet1_phi);
  stream.select("patJet_cleanPatJetsAK5Calo.pt", jet1_pt);
  stream.select("patJet_cleanPatJetsAK5Calo.simpleSecondaryVertexBJetTags", jet1_simpleSecondaryVertexBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.simpleSecondaryVertexHighEffBJetTags", jet1_simpleSecondaryVertexHighEffBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.simpleSecondaryVertexHighPurBJetTags", jet1_simpleSecondaryVertexHighPurBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.softElectronByIP3dBJetTags", jet1_softElectronByIP3dBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.softElectronByPtBJetTags", jet1_softElectronByPtBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.softMuonBJetTags", jet1_softMuonBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.softMuonByIP3dBJetTags", jet1_softMuonByIP3dBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.softMuonByPtBJetTags", jet1_softMuonByPtBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.trackCountingHighEffBJetTags", jet1_trackCountingHighEffBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.trackCountingHighPurBJetTags", jet1_trackCountingHighPurBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.uncor_energy", jet1_uncor_energy);
  stream.select("patJet_cleanPatJetsAK5Calo.uncor_eta", jet1_uncor_eta);
  stream.select("patJet_cleanPatJetsAK5Calo.uncor_phi", jet1_uncor_phi);
  stream.select("patJet_cleanPatJetsAK5Calo.uncor_pt", jet1_uncor_pt);
  stream.select("patJet_selectedPatJetsPF.HFHadronEnergy", jet2_HFHadronEnergy);
  stream.select("patJet_selectedPatJetsPF.chargedEmEnergyFraction", jet2_chargedEmEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.chargedHadronEnergyFraction", jet2_chargedHadronEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.chargedMultiplicity", jet2_chargedMultiplicity);
  stream.select("patJet_selectedPatJetsPF.combinedSecondaryVertexBJetTags", jet2_combinedSecondaryVertexBJetTags);
  stream.select("patJet_selectedPatJetsPF.combinedSecondaryVertexMVABJetTags", jet2_combinedSecondaryVertexMVABJetTags);
  stream.select("patJet_selectedPatJetsPF.energy", jet2_energy);
  stream.select("patJet_selectedPatJetsPF.eta", jet2_eta);
  stream.select("patJet_selectedPatJetsPF.jetArea", jet2_jetArea);
  stream.select("patJet_selectedPatJetsPF.jetBProbabilityBJetTags", jet2_jetBProbabilityBJetTags);
  stream.select("patJet_selectedPatJetsPF.jetID_fHPD", jet2_jetID_fHPD);
  stream.select("patJet_selectedPatJetsPF.jetID_n90Hits", jet2_jetID_n90Hits);
  stream.select("patJet_selectedPatJetsPF.jetProbabilityBJetTags", jet2_jetProbabilityBJetTags);
  stream.select("patJet_selectedPatJetsPF.neutralEmEnergyFraction", jet2_neutralEmEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.neutralHadronEnergy", jet2_neutralHadronEnergy);
  stream.select("patJet_selectedPatJetsPF.neutralHadronEnergyFraction", jet2_neutralHadronEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.numberOfDaughters", jet2_numberOfDaughters);
  stream.select("patJet_selectedPatJetsPF.partonFlavour", jet2_partonFlavour);
  stream.select("patJet_selectedPatJetsPF.phi", jet2_phi);
  stream.select("patJet_selectedPatJetsPF.pt", jet2_pt);
  stream.select("patJet_selectedPatJetsPF.simpleSecondaryVertexBJetTags", jet2_simpleSecondaryVertexBJetTags);
  stream.select("patJet_selectedPatJetsPF.simpleSecondaryVertexHighEffBJetTags", jet2_simpleSecondaryVertexHighEffBJetTags);
  stream.select("patJet_selectedPatJetsPF.simpleSecondaryVertexHighPurBJetTags", jet2_simpleSecondaryVertexHighPurBJetTags);
  stream.select("patJet_selectedPatJetsPF.softElectronByIP3dBJetTags", jet2_softElectronByIP3dBJetTags);
  stream.select("patJet_selectedPatJetsPF.softElectronByPtBJetTags", jet2_softElectronByPtBJetTags);
  stream.select("patJet_selectedPatJetsPF.softMuonBJetTags", jet2_softMuonBJetTags);
  stream.select("patJet_selectedPatJetsPF.softMuonByIP3dBJetTags", jet2_softMuonByIP3dBJetTags);
  stream.select("patJet_selectedPatJetsPF.softMuonByPtBJetTags", jet2_softMuonByPtBJetTags);
  stream.select("patJet_selectedPatJetsPF.trackCountingHighEffBJetTags", jet2_trackCountingHighEffBJetTags);
  stream.select("patJet_selectedPatJetsPF.trackCountingHighPurBJetTags", jet2_trackCountingHighPurBJetTags);
  stream.select("patJet_selectedPatJetsPF.uncor_energy", jet2_uncor_energy);
  stream.select("patJet_selectedPatJetsPF.uncor_eta", jet2_uncor_eta);
  stream.select("patJet_selectedPatJetsPF.uncor_phi", jet2_uncor_phi);
  stream.select("patJet_selectedPatJetsPF.uncor_pt", jet2_uncor_pt);
  stream.select("patJetHelper_selectedPatJetsPF.jecFactor", jet_jecFactor);
  stream.select("patJetHelper_selectedPatJetsPF.jecFactorNoL1Fast", jet_jecFactorNoL1Fast);
  stream.select("patJetHelper_selectedPatJetsPF.jetUncMinus", jet_jetUncMinus);
  stream.select("patJetHelper_selectedPatJetsPF.jetUncPlus", jet_jetUncPlus);
  stream.select("patMET_patMETsPF.energy", met1_energy);
  stream.select("patMET_patMETsPF.et", met1_et);
  stream.select("patMET_patMETsPF.mEtSig", met1_mEtSig);
  stream.select("patMET_patMETsPF.phi", met1_phi);
  stream.select("patMET_patMETsPF.pt", met1_pt);
  stream.select("patMET_patMETsPF.sumEt", met1_sumEt);
  stream.select("patMET_patMETsTypeIPF.energy", met2_energy);
  stream.select("patMET_patMETsTypeIPF.et", met2_et);
  stream.select("patMET_patMETsTypeIPF.mEtSig", met2_mEtSig);
  stream.select("patMET_patMETsTypeIPF.phi", met2_phi);
  stream.select("patMET_patMETsTypeIPF.pt", met2_pt);
  stream.select("patMET_patMETsTypeIPF.sumEt", met2_sumEt);
  stream.select("patMET_patMETsAK5Calo.energy", met_energy);
  stream.select("patMET_patMETsAK5Calo.et", met_et);
  stream.select("patMET_patMETsAK5Calo.mEtSig", met_mEtSig);
  stream.select("patMET_patMETsAK5Calo.phi", met_phi);
  stream.select("patMET_patMETsAK5Calo.pt", met_pt);
  stream.select("patMET_patMETsAK5Calo.sumEt", met_sumEt);
  stream.select("patMuon_selectedPatMuonsPF.AllGlobalMuons", muon1_AllGlobalMuons);
  stream.select("patMuon_selectedPatMuonsPF.GlobalMuonPromptTight", muon1_GlobalMuonPromptTight);
  stream.select("patMuon_selectedPatMuonsPF.charge", muon1_charge);
  stream.select("patMuon_selectedPatMuonsPF.chargedHadronIso", muon1_chargedHadronIso);
  stream.select("patMuon_selectedPatMuonsPF.combinedMuon_chi2", muon1_combinedMuon_chi2);
  stream.select("patMuon_selectedPatMuonsPF.combinedMuon_ndof", muon1_combinedMuon_ndof);
  stream.select("patMuon_selectedPatMuonsPF.combinedMuon_numberOfValidHits", muon1_combinedMuon_numberOfValidHits);
  stream.select("patMuon_selectedPatMuonsPF.dB", muon1_dB);
  stream.select("patMuon_selectedPatMuonsPF.ecalIso", muon1_ecalIso);
  stream.select("patMuon_selectedPatMuonsPF.energy", muon1_energy);
  stream.select("patMuon_selectedPatMuonsPF.eta", muon1_eta);
  stream.select("patMuon_selectedPatMuonsPF.genLepton_eta", muon1_genLepton_eta);
  stream.select("patMuon_selectedPatMuonsPF.genLepton_pdgId", muon1_genLepton_pdgId);
  stream.select("patMuon_selectedPatMuonsPF.genLepton_phi", muon1_genLepton_phi);
  stream.select("patMuon_selectedPatMuonsPF.genLepton_pt", muon1_genLepton_pt);
  stream.select("patMuon_selectedPatMuonsPF.globalTrack_pt", muon1_globalTrack_pt);
  stream.select("patMuon_selectedPatMuonsPF.hcalIso", muon1_hcalIso);
  stream.select("patMuon_selectedPatMuonsPF.innerTrack_numberOfValidHits", muon1_innerTrack_numberOfValidHits);
  stream.select("patMuon_selectedPatMuonsPF.innerTrack_pt", muon1_innerTrack_pt);
  stream.select("patMuon_selectedPatMuonsPF.isGlobalMuon", muon1_isGlobalMuon);
  stream.select("patMuon_selectedPatMuonsPF.isTrackerMuon", muon1_isTrackerMuon);
  stream.select("patMuon_selectedPatMuonsPF.neutralHadronIso", muon1_neutralHadronIso);
  stream.select("patMuon_selectedPatMuonsPF.phi", muon1_phi);
  stream.select("patMuon_selectedPatMuonsPF.photonIso", muon1_photonIso);
  stream.select("patMuon_selectedPatMuonsPF.pt", muon1_pt);
  stream.select("patMuon_selectedPatMuonsPF.px", muon1_px);
  stream.select("patMuon_selectedPatMuonsPF.py", muon1_py);
  stream.select("patMuon_selectedPatMuonsPF.pz", muon1_pz);
  stream.select("patMuon_selectedPatMuonsPF.trackIso", muon1_trackIso);
  stream.select("patMuon_selectedPatMuonsPF.track_d0", muon1_track_d0);
  stream.select("patMuon_selectedPatMuonsPF.track_hitPattern_numberOfValidMuonHits", muon1_track_hitPattern_numberOfValidMuonHits);
  stream.select("patMuon_selectedPatMuonsPF.track_hitPattern_numberOfValidPixelHits", muon1_track_hitPattern_numberOfValidPixelHits);
  stream.select("patMuon_selectedPatMuonsPF.track_normalizedChi2", muon1_track_normalizedChi2);
  stream.select("patMuon_selectedPatMuonsPF.vz", muon1_vz);
  stream.select("patMuonHelper_selectedPatMuonsPF.dxywrtBeamSpot", muon_dxywrtBeamSpot);
  stream.select("patMuonHelper_selectedPatMuonsPF.muonPFcandptdiff", muon_muonPFcandptdiff);
  stream.select("nGenEventInfoProductHelper_generator", nGenEventInfoProductHelper_generator);
  stream.select("npatElectronHelper_selectedPatElectronsPF", nelectron);
  stream.select("npatElectron_selectedPatElectronsPF", nelectron1);
  stream.select("nrecoGenParticleHelperRA2b_genParticles", ngenparticlera2b);
  stream.select("npatJetHelper_selectedPatJetsPF", njet);
  stream.select("npatJet_cleanPatJetsAK5Calo", njet1);
  stream.select("npatJet_selectedPatJetsPF", njet2);
  stream.select("npatMET_patMETsAK5Calo", nmet);
  stream.select("npatMET_patMETsPF", nmet1);
  stream.select("npatMET_patMETsTypeIPF", nmet2);
  stream.select("npatMuonHelper_selectedPatMuonsPF", nmuon);
  stream.select("npatMuon_selectedPatMuonsPF", nmuon1);
  stream.select("npatTauHelper_selectedPatTausPF", ntau);
  stream.select("npatTau_selectedPatTausPF", ntau1);
  stream.select("nrecoVertex_offlinePrimaryVertices", nvertex);
  stream.select("sdouble_kt6PFJets_rho.value", sdouble_kt6pfjets_rho_value);
  stream.select("sdouble_kt6PFJets_sigma.value", sdouble_kt6pfjets_sigma_value);
  stream.select("patTau_selectedPatTausPF.caloIso", tau1_caloIso);
  stream.select("patTau_selectedPatTausPF.ecalIso", tau1_ecalIso);
  stream.select("patTau_selectedPatTausPF.energy", tau1_energy);
  stream.select("patTau_selectedPatTausPF.et", tau1_et);
  stream.select("patTau_selectedPatTausPF.eta", tau1_eta);
  stream.select("patTau_selectedPatTausPF.hcalIso", tau1_hcalIso);
  stream.select("patTau_selectedPatTausPF.phi", tau1_phi);
  stream.select("patTau_selectedPatTausPF.pt", tau1_pt);
  stream.select("patTau_selectedPatTausPF.px", tau1_px);
  stream.select("patTau_selectedPatTausPF.py", tau1_py);
  stream.select("patTau_selectedPatTausPF.pz", tau1_pz);
  stream.select("patTau_selectedPatTausPF.tauID_againstElectron", tau1_tauID_againstElectron);
  stream.select("patTau_selectedPatTausPF.tauID_againstMuon", tau1_tauID_againstMuon);
  stream.select("patTau_selectedPatTausPF.tauID_byIsolation", tau1_tauID_byIsolation);
  stream.select("patTau_selectedPatTausPF.tauID_byTaNC", tau1_tauID_byTaNC);
  stream.select("patTau_selectedPatTausPF.tauID_byTaNCfrHalfPercent", tau1_tauID_byTaNCfrHalfPercent);
  stream.select("patTau_selectedPatTausPF.tauID_byTaNCfrQuarterPercent", tau1_tauID_byTaNCfrQuarterPercent);
  stream.select("patTau_selectedPatTausPF.trackIso", tau1_trackIso);
  stream.select("patTauHelper_selectedPatTausPF.genTauDecayModeID", tau_genTauDecayModeID);
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
  double	dxywrtBeamSpot;
};
std::vector<electron_s> electron(7);

std::ostream& operator<<(std::ostream& os, const electron_s& o)
{
  char r[1024];
  os << "electron" << std::endl;
  sprintf(r, "  %-32s: %f\n", "dxywrtBeamSpot", (double)o.dxywrtBeamSpot); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct electron1_s
{
  double	caloIso;
  double	charge;
  double	chargedHadronIso;
  double	dB;
  double	deltaEtaSuperClusterTrackAtVtx;
  double	deltaPhiSuperClusterTrackAtVtx;
  double	dr03EcalRecHitSumEt;
  double	dr03HcalTowerSumEt;
  double	dr03TkSumPt;
  double	ecalIso;
  double	eidRobustTight;
  double	energy;
  double	et;
  double	eta;
  double	genLepton_eta;
  double	genLepton_pdgId;
  double	genLepton_phi;
  double	genLepton_pt;
  double	gsfTrack_d0;
  double	gsfTrack_trackerExpectedHitsInner_numberOfLostHits;
  double	hadronicOverEm;
  double	hcalIso;
  double	neutralHadronIso;
  double	phi;
  double	photonIso;
  double	pt;
  double	px;
  double	py;
  double	pz;
  double	sigmaIetaIeta;
  double	simpleEleId95relIso;
  double	superCluster_eta;
  double	trackIso;
  double	vz;
};
std::vector<electron1_s> electron1(7);

std::ostream& operator<<(std::ostream& os, const electron1_s& o)
{
  char r[1024];
  os << "electron1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "caloIso", (double)o.caloIso); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronIso", (double)o.chargedHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "dB", (double)o.dB); os << r;
  sprintf(r, "  %-32s: %f\n", "deltaEtaSuperClusterTrackAtVtx", (double)o.deltaEtaSuperClusterTrackAtVtx); os << r;
  sprintf(r, "  %-32s: %f\n", "deltaPhiSuperClusterTrackAtVtx", (double)o.deltaPhiSuperClusterTrackAtVtx); os << r;
  sprintf(r, "  %-32s: %f\n", "dr03EcalRecHitSumEt", (double)o.dr03EcalRecHitSumEt); os << r;
  sprintf(r, "  %-32s: %f\n", "dr03HcalTowerSumEt", (double)o.dr03HcalTowerSumEt); os << r;
  sprintf(r, "  %-32s: %f\n", "dr03TkSumPt", (double)o.dr03TkSumPt); os << r;
  sprintf(r, "  %-32s: %f\n", "ecalIso", (double)o.ecalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "eidRobustTight", (double)o.eidRobustTight); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_eta", (double)o.genLepton_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_pdgId", (double)o.genLepton_pdgId); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_phi", (double)o.genLepton_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_pt", (double)o.genLepton_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "gsfTrack_d0", (double)o.gsfTrack_d0); os << r;
  sprintf(r, "  %-32s: %f\n", "gsfTrack_trackerExpectedHitsInner_numberOfLostHits", (double)o.gsfTrack_trackerExpectedHitsInner_numberOfLostHits); os << r;
  sprintf(r, "  %-32s: %f\n", "hadronicOverEm", (double)o.hadronicOverEm); os << r;
  sprintf(r, "  %-32s: %f\n", "hcalIso", (double)o.hcalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronIso", (double)o.neutralHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "photonIso", (double)o.photonIso); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "px", (double)o.px); os << r;
  sprintf(r, "  %-32s: %f\n", "py", (double)o.py); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "sigmaIetaIeta", (double)o.sigmaIetaIeta); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleEleId95relIso", (double)o.simpleEleId95relIso); os << r;
  sprintf(r, "  %-32s: %f\n", "superCluster_eta", (double)o.superCluster_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "trackIso", (double)o.trackIso); os << r;
  sprintf(r, "  %-32s: %f\n", "vz", (double)o.vz); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jet_s
{
  double	jecFactor;
  double	jecFactorNoL1Fast;
  double	jetUncMinus;
  double	jetUncPlus;
};
std::vector<jet_s> jet(173);

std::ostream& operator<<(std::ostream& os, const jet_s& o)
{
  char r[1024];
  os << "jet" << std::endl;
  sprintf(r, "  %-32s: %f\n", "jecFactor", (double)o.jecFactor); os << r;
  sprintf(r, "  %-32s: %f\n", "jecFactorNoL1Fast", (double)o.jecFactorNoL1Fast); os << r;
  sprintf(r, "  %-32s: %f\n", "jetUncMinus", (double)o.jetUncMinus); os << r;
  sprintf(r, "  %-32s: %f\n", "jetUncPlus", (double)o.jetUncPlus); os << r;
  return os;
}
//-----------------------------------------------------------------------------

struct jet1_s
{
  double	combinedSecondaryVertexBJetTags;
  double	combinedSecondaryVertexMVABJetTags;
  double	emEnergyFraction;
  double	energy;
  double	eta;
  double	jetBProbabilityBJetTags;
  double	jetID_fHPD;
  double	jetID_n90Hits;
  double	jetProbabilityBJetTags;
  double	partonFlavour;
  double	phi;
  double	pt;
  double	simpleSecondaryVertexBJetTags;
  double	simpleSecondaryVertexHighEffBJetTags;
  double	simpleSecondaryVertexHighPurBJetTags;
  double	softElectronByIP3dBJetTags;
  double	softElectronByPtBJetTags;
  double	softMuonBJetTags;
  double	softMuonByIP3dBJetTags;
  double	softMuonByPtBJetTags;
  double	trackCountingHighEffBJetTags;
  double	trackCountingHighPurBJetTags;
  double	uncor_energy;
  double	uncor_eta;
  double	uncor_phi;
  double	uncor_pt;
};
std::vector<jet1_s> jet1(177);

std::ostream& operator<<(std::ostream& os, const jet1_s& o)
{
  char r[1024];
  os << "jet1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexBJetTags", (double)o.combinedSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexMVABJetTags", (double)o.combinedSecondaryVertexMVABJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "emEnergyFraction", (double)o.emEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "jetBProbabilityBJetTags", (double)o.jetBProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetID_fHPD", (double)o.jetID_fHPD); os << r;
  sprintf(r, "  %-32s: %f\n", "jetID_n90Hits", (double)o.jetID_n90Hits); os << r;
  sprintf(r, "  %-32s: %f\n", "jetProbabilityBJetTags", (double)o.jetProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "partonFlavour", (double)o.partonFlavour); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleSecondaryVertexBJetTags", (double)o.simpleSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleSecondaryVertexHighEffBJetTags", (double)o.simpleSecondaryVertexHighEffBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleSecondaryVertexHighPurBJetTags", (double)o.simpleSecondaryVertexHighPurBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softElectronByIP3dBJetTags", (double)o.softElectronByIP3dBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softElectronByPtBJetTags", (double)o.softElectronByPtBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softMuonBJetTags", (double)o.softMuonBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softMuonByIP3dBJetTags", (double)o.softMuonByIP3dBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softMuonByPtBJetTags", (double)o.softMuonByPtBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighEffBJetTags", (double)o.trackCountingHighEffBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighPurBJetTags", (double)o.trackCountingHighPurBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_energy", (double)o.uncor_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_eta", (double)o.uncor_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_phi", (double)o.uncor_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_pt", (double)o.uncor_pt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jet2_s
{
  double	HFHadronEnergy;
  double	chargedEmEnergyFraction;
  double	chargedHadronEnergyFraction;
  double	chargedMultiplicity;
  double	combinedSecondaryVertexBJetTags;
  double	combinedSecondaryVertexMVABJetTags;
  double	energy;
  double	eta;
  double	jetArea;
  double	jetBProbabilityBJetTags;
  double	jetID_fHPD;
  double	jetID_n90Hits;
  double	jetProbabilityBJetTags;
  double	neutralEmEnergyFraction;
  double	neutralHadronEnergy;
  double	neutralHadronEnergyFraction;
  double	numberOfDaughters;
  double	partonFlavour;
  double	phi;
  double	pt;
  double	simpleSecondaryVertexBJetTags;
  double	simpleSecondaryVertexHighEffBJetTags;
  double	simpleSecondaryVertexHighPurBJetTags;
  double	softElectronByIP3dBJetTags;
  double	softElectronByPtBJetTags;
  double	softMuonBJetTags;
  double	softMuonByIP3dBJetTags;
  double	softMuonByPtBJetTags;
  double	trackCountingHighEffBJetTags;
  double	trackCountingHighPurBJetTags;
  double	uncor_energy;
  double	uncor_eta;
  double	uncor_phi;
  double	uncor_pt;
};
std::vector<jet2_s> jet2(173);

std::ostream& operator<<(std::ostream& os, const jet2_s& o)
{
  char r[1024];
  os << "jet2" << std::endl;
  sprintf(r, "  %-32s: %f\n", "HFHadronEnergy", (double)o.HFHadronEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedEmEnergyFraction", (double)o.chargedEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronEnergyFraction", (double)o.chargedHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedMultiplicity", (double)o.chargedMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexBJetTags", (double)o.combinedSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexMVABJetTags", (double)o.combinedSecondaryVertexMVABJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "jetArea", (double)o.jetArea); os << r;
  sprintf(r, "  %-32s: %f\n", "jetBProbabilityBJetTags", (double)o.jetBProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetID_fHPD", (double)o.jetID_fHPD); os << r;
  sprintf(r, "  %-32s: %f\n", "jetID_n90Hits", (double)o.jetID_n90Hits); os << r;
  sprintf(r, "  %-32s: %f\n", "jetProbabilityBJetTags", (double)o.jetProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralEmEnergyFraction", (double)o.neutralEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergy", (double)o.neutralHadronEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergyFraction", (double)o.neutralHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfDaughters", (double)o.numberOfDaughters); os << r;
  sprintf(r, "  %-32s: %f\n", "partonFlavour", (double)o.partonFlavour); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleSecondaryVertexBJetTags", (double)o.simpleSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleSecondaryVertexHighEffBJetTags", (double)o.simpleSecondaryVertexHighEffBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleSecondaryVertexHighPurBJetTags", (double)o.simpleSecondaryVertexHighPurBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softElectronByIP3dBJetTags", (double)o.softElectronByIP3dBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softElectronByPtBJetTags", (double)o.softElectronByPtBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softMuonBJetTags", (double)o.softMuonBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softMuonByIP3dBJetTags", (double)o.softMuonByIP3dBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "softMuonByPtBJetTags", (double)o.softMuonByPtBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighEffBJetTags", (double)o.trackCountingHighEffBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighPurBJetTags", (double)o.trackCountingHighPurBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_energy", (double)o.uncor_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_eta", (double)o.uncor_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_phi", (double)o.uncor_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_pt", (double)o.uncor_pt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct met_s
{
  double	energy;
  double	et;
  double	mEtSig;
  double	phi;
  double	pt;
  double	sumEt;
};
std::vector<met_s> met(3);

std::ostream& operator<<(std::ostream& os, const met_s& o)
{
  char r[1024];
  os << "met" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct met1_s
{
  double	energy;
  double	et;
  double	mEtSig;
  double	phi;
  double	pt;
  double	sumEt;
};
std::vector<met1_s> met1(3);

std::ostream& operator<<(std::ostream& os, const met1_s& o)
{
  char r[1024];
  os << "met1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct met2_s
{
  double	energy;
  double	et;
  double	mEtSig;
  double	phi;
  double	pt;
  double	sumEt;
};
std::vector<met2_s> met2(3);

std::ostream& operator<<(std::ostream& os, const met2_s& o)
{
  char r[1024];
  os << "met2" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct muon_s
{
  double	dxywrtBeamSpot;
  double	muonPFcandptdiff;
};
std::vector<muon_s> muon(5);

std::ostream& operator<<(std::ostream& os, const muon_s& o)
{
  char r[1024];
  os << "muon" << std::endl;
  sprintf(r, "  %-32s: %f\n", "dxywrtBeamSpot", (double)o.dxywrtBeamSpot); os << r;
  sprintf(r, "  %-32s: %f\n", "muonPFcandptdiff", (double)o.muonPFcandptdiff); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct muon1_s
{
  double	AllGlobalMuons;
  double	GlobalMuonPromptTight;
  double	charge;
  double	chargedHadronIso;
  double	combinedMuon_chi2;
  double	combinedMuon_ndof;
  double	combinedMuon_numberOfValidHits;
  double	dB;
  double	ecalIso;
  double	energy;
  double	eta;
  double	genLepton_eta;
  double	genLepton_pdgId;
  double	genLepton_phi;
  double	genLepton_pt;
  double	globalTrack_pt;
  double	hcalIso;
  double	innerTrack_numberOfValidHits;
  double	innerTrack_pt;
  double	isGlobalMuon;
  double	isTrackerMuon;
  double	neutralHadronIso;
  double	phi;
  double	photonIso;
  double	pt;
  double	px;
  double	py;
  double	pz;
  double	trackIso;
  double	track_d0;
  double	track_hitPattern_numberOfValidMuonHits;
  double	track_hitPattern_numberOfValidPixelHits;
  double	track_normalizedChi2;
  double	vz;
};
std::vector<muon1_s> muon1(5);

std::ostream& operator<<(std::ostream& os, const muon1_s& o)
{
  char r[1024];
  os << "muon1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "AllGlobalMuons", (double)o.AllGlobalMuons); os << r;
  sprintf(r, "  %-32s: %f\n", "GlobalMuonPromptTight", (double)o.GlobalMuonPromptTight); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronIso", (double)o.chargedHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedMuon_chi2", (double)o.combinedMuon_chi2); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedMuon_ndof", (double)o.combinedMuon_ndof); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedMuon_numberOfValidHits", (double)o.combinedMuon_numberOfValidHits); os << r;
  sprintf(r, "  %-32s: %f\n", "dB", (double)o.dB); os << r;
  sprintf(r, "  %-32s: %f\n", "ecalIso", (double)o.ecalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_eta", (double)o.genLepton_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_pdgId", (double)o.genLepton_pdgId); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_phi", (double)o.genLepton_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_pt", (double)o.genLepton_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_pt", (double)o.globalTrack_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "hcalIso", (double)o.hcalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_numberOfValidHits", (double)o.innerTrack_numberOfValidHits); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_pt", (double)o.innerTrack_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "isGlobalMuon", (double)o.isGlobalMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "isTrackerMuon", (double)o.isTrackerMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronIso", (double)o.neutralHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "photonIso", (double)o.photonIso); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "px", (double)o.px); os << r;
  sprintf(r, "  %-32s: %f\n", "py", (double)o.py); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "trackIso", (double)o.trackIso); os << r;
  sprintf(r, "  %-32s: %f\n", "track_d0", (double)o.track_d0); os << r;
  sprintf(r, "  %-32s: %f\n", "track_hitPattern_numberOfValidMuonHits", (double)o.track_hitPattern_numberOfValidMuonHits); os << r;
  sprintf(r, "  %-32s: %f\n", "track_hitPattern_numberOfValidPixelHits", (double)o.track_hitPattern_numberOfValidPixelHits); os << r;
  sprintf(r, "  %-32s: %f\n", "track_normalizedChi2", (double)o.track_normalizedChi2); os << r;
  sprintf(r, "  %-32s: %f\n", "vz", (double)o.vz); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct tau_s
{
  double	genTauDecayModeID;
};
std::vector<tau_s> tau(9);

std::ostream& operator<<(std::ostream& os, const tau_s& o)
{
  char r[1024];
  os << "tau" << std::endl;
  sprintf(r, "  %-32s: %f\n", "genTauDecayModeID", (double)o.genTauDecayModeID); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct tau1_s
{
  double	caloIso;
  double	ecalIso;
  double	energy;
  double	et;
  double	eta;
  double	hcalIso;
  double	phi;
  double	pt;
  double	px;
  double	py;
  double	pz;
  double	tauID_againstElectron;
  double	tauID_againstMuon;
  double	tauID_byIsolation;
  double	tauID_byTaNC;
  double	tauID_byTaNCfrHalfPercent;
  double	tauID_byTaNCfrQuarterPercent;
  double	trackIso;
};
std::vector<tau1_s> tau1(9);

std::ostream& operator<<(std::ostream& os, const tau1_s& o)
{
  char r[1024];
  os << "tau1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "caloIso", (double)o.caloIso); os << r;
  sprintf(r, "  %-32s: %f\n", "ecalIso", (double)o.ecalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "hcalIso", (double)o.hcalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "px", (double)o.px); os << r;
  sprintf(r, "  %-32s: %f\n", "py", (double)o.py); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "tauID_againstElectron", (double)o.tauID_againstElectron); os << r;
  sprintf(r, "  %-32s: %f\n", "tauID_againstMuon", (double)o.tauID_againstMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "tauID_byIsolation", (double)o.tauID_byIsolation); os << r;
  sprintf(r, "  %-32s: %f\n", "tauID_byTaNC", (double)o.tauID_byTaNC); os << r;
  sprintf(r, "  %-32s: %f\n", "tauID_byTaNCfrHalfPercent", (double)o.tauID_byTaNCfrHalfPercent); os << r;
  sprintf(r, "  %-32s: %f\n", "tauID_byTaNCfrQuarterPercent", (double)o.tauID_byTaNCfrQuarterPercent); os << r;
  sprintf(r, "  %-32s: %f\n", "trackIso", (double)o.trackIso); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct vertex_s
{
  double	isFake;
  double	isValid;
  double	ndof;
  double	normalizedChi2;
  double	position_Rho;
  double	tracksSize;
  double	x;
  double	xError;
  double	y;
  double	yError;
  double	z;
  double	zError;
};
std::vector<vertex_s> vertex(35);

std::ostream& operator<<(std::ostream& os, const vertex_s& o)
{
  char r[1024];
  os << "vertex" << std::endl;
  sprintf(r, "  %-32s: %f\n", "isFake", (double)o.isFake); os << r;
  sprintf(r, "  %-32s: %f\n", "isValid", (double)o.isValid); os << r;
  sprintf(r, "  %-32s: %f\n", "ndof", (double)o.ndof); os << r;
  sprintf(r, "  %-32s: %f\n", "normalizedChi2", (double)o.normalizedChi2); os << r;
  sprintf(r, "  %-32s: %f\n", "position_Rho", (double)o.position_Rho); os << r;
  sprintf(r, "  %-32s: %f\n", "tracksSize", (double)o.tracksSize); os << r;
  sprintf(r, "  %-32s: %f\n", "x", (double)o.x); os << r;
  sprintf(r, "  %-32s: %f\n", "xError", (double)o.xError); os << r;
  sprintf(r, "  %-32s: %f\n", "y", (double)o.y); os << r;
  sprintf(r, "  %-32s: %f\n", "yError", (double)o.yError); os << r;
  sprintf(r, "  %-32s: %f\n", "z", (double)o.z); os << r;
  sprintf(r, "  %-32s: %f\n", "zError", (double)o.zError); os << r;
  return os;
}
//-----------------------------------------------------------------------------

void fillObjects()
{

  electron.resize(electron_dxywrtBeamSpot.size());
  for(unsigned int i=0; i < electron_dxywrtBeamSpot.size(); ++i)
    {
      electron[i].dxywrtBeamSpot	= electron_dxywrtBeamSpot[i];
    }

  electron1.resize(electron1_caloIso.size());
  for(unsigned int i=0; i < electron1_caloIso.size(); ++i)
    {
      electron1[i].caloIso	= electron1_caloIso[i];
      electron1[i].charge	= electron1_charge[i];
      electron1[i].chargedHadronIso	= electron1_chargedHadronIso[i];
      electron1[i].dB	= electron1_dB[i];
      electron1[i].deltaEtaSuperClusterTrackAtVtx	= electron1_deltaEtaSuperClusterTrackAtVtx[i];
      electron1[i].deltaPhiSuperClusterTrackAtVtx	= electron1_deltaPhiSuperClusterTrackAtVtx[i];
      electron1[i].dr03EcalRecHitSumEt	= electron1_dr03EcalRecHitSumEt[i];
      electron1[i].dr03HcalTowerSumEt	= electron1_dr03HcalTowerSumEt[i];
      electron1[i].dr03TkSumPt	= electron1_dr03TkSumPt[i];
      electron1[i].ecalIso	= electron1_ecalIso[i];
      electron1[i].eidRobustTight	= electron1_eidRobustTight[i];
      electron1[i].energy	= electron1_energy[i];
      electron1[i].et	= electron1_et[i];
      electron1[i].eta	= electron1_eta[i];
      electron1[i].genLepton_eta	= electron1_genLepton_eta[i];
      electron1[i].genLepton_pdgId	= electron1_genLepton_pdgId[i];
      electron1[i].genLepton_phi	= electron1_genLepton_phi[i];
      electron1[i].genLepton_pt	= electron1_genLepton_pt[i];
      electron1[i].gsfTrack_d0	= electron1_gsfTrack_d0[i];
      electron1[i].gsfTrack_trackerExpectedHitsInner_numberOfLostHits	= electron1_gsfTrack_trackerExpectedHitsInner_numberOfLostHits[i];
      electron1[i].hadronicOverEm	= electron1_hadronicOverEm[i];
      electron1[i].hcalIso	= electron1_hcalIso[i];
      electron1[i].neutralHadronIso	= electron1_neutralHadronIso[i];
      electron1[i].phi	= electron1_phi[i];
      electron1[i].photonIso	= electron1_photonIso[i];
      electron1[i].pt	= electron1_pt[i];
      electron1[i].px	= electron1_px[i];
      electron1[i].py	= electron1_py[i];
      electron1[i].pz	= electron1_pz[i];
      electron1[i].sigmaIetaIeta	= electron1_sigmaIetaIeta[i];
      electron1[i].simpleEleId95relIso	= electron1_simpleEleId95relIso[i];
      electron1[i].superCluster_eta	= electron1_superCluster_eta[i];
      electron1[i].trackIso	= electron1_trackIso[i];
      electron1[i].vz	= electron1_vz[i];
    }

  jet.resize(jet_jecFactor.size());
  for(unsigned int i=0; i < jet_jecFactor.size(); ++i)
    {
      jet[i].jecFactor	= jet_jecFactor[i];
      jet[i].jecFactorNoL1Fast	= jet_jecFactorNoL1Fast[i];
      jet[i].jetUncMinus	= jet_jetUncMinus[i];
      jet[i].jetUncPlus	= jet_jetUncPlus[i];
    }

  jet1.resize(jet1_combinedSecondaryVertexBJetTags.size());
  for(unsigned int i=0; i < jet1_combinedSecondaryVertexBJetTags.size(); ++i)
    {
      jet1[i].combinedSecondaryVertexBJetTags	= jet1_combinedSecondaryVertexBJetTags[i];
      jet1[i].combinedSecondaryVertexMVABJetTags	= jet1_combinedSecondaryVertexMVABJetTags[i];
      jet1[i].emEnergyFraction	= jet1_emEnergyFraction[i];
      jet1[i].energy	= jet1_energy[i];
      jet1[i].eta	= jet1_eta[i];
      jet1[i].jetBProbabilityBJetTags	= jet1_jetBProbabilityBJetTags[i];
      jet1[i].jetID_fHPD	= jet1_jetID_fHPD[i];
      jet1[i].jetID_n90Hits	= jet1_jetID_n90Hits[i];
      jet1[i].jetProbabilityBJetTags	= jet1_jetProbabilityBJetTags[i];
      jet1[i].partonFlavour	= jet1_partonFlavour[i];
      jet1[i].phi	= jet1_phi[i];
      jet1[i].pt	= jet1_pt[i];
      jet1[i].simpleSecondaryVertexBJetTags	= jet1_simpleSecondaryVertexBJetTags[i];
      jet1[i].simpleSecondaryVertexHighEffBJetTags	= jet1_simpleSecondaryVertexHighEffBJetTags[i];
      jet1[i].simpleSecondaryVertexHighPurBJetTags	= jet1_simpleSecondaryVertexHighPurBJetTags[i];
      jet1[i].softElectronByIP3dBJetTags	= jet1_softElectronByIP3dBJetTags[i];
      jet1[i].softElectronByPtBJetTags	= jet1_softElectronByPtBJetTags[i];
      jet1[i].softMuonBJetTags	= jet1_softMuonBJetTags[i];
      jet1[i].softMuonByIP3dBJetTags	= jet1_softMuonByIP3dBJetTags[i];
      jet1[i].softMuonByPtBJetTags	= jet1_softMuonByPtBJetTags[i];
      jet1[i].trackCountingHighEffBJetTags	= jet1_trackCountingHighEffBJetTags[i];
      jet1[i].trackCountingHighPurBJetTags	= jet1_trackCountingHighPurBJetTags[i];
      jet1[i].uncor_energy	= jet1_uncor_energy[i];
      jet1[i].uncor_eta	= jet1_uncor_eta[i];
      jet1[i].uncor_phi	= jet1_uncor_phi[i];
      jet1[i].uncor_pt	= jet1_uncor_pt[i];
    }

  jet2.resize(jet2_HFHadronEnergy.size());
  for(unsigned int i=0; i < jet2_HFHadronEnergy.size(); ++i)
    {
      jet2[i].HFHadronEnergy	= jet2_HFHadronEnergy[i];
      jet2[i].chargedEmEnergyFraction	= jet2_chargedEmEnergyFraction[i];
      jet2[i].chargedHadronEnergyFraction	= jet2_chargedHadronEnergyFraction[i];
      jet2[i].chargedMultiplicity	= jet2_chargedMultiplicity[i];
      jet2[i].combinedSecondaryVertexBJetTags	= jet2_combinedSecondaryVertexBJetTags[i];
      jet2[i].combinedSecondaryVertexMVABJetTags	= jet2_combinedSecondaryVertexMVABJetTags[i];
      jet2[i].energy	= jet2_energy[i];
      jet2[i].eta	= jet2_eta[i];
      jet2[i].jetArea	= jet2_jetArea[i];
      jet2[i].jetBProbabilityBJetTags	= jet2_jetBProbabilityBJetTags[i];
      jet2[i].jetID_fHPD	= jet2_jetID_fHPD[i];
      jet2[i].jetID_n90Hits	= jet2_jetID_n90Hits[i];
      jet2[i].jetProbabilityBJetTags	= jet2_jetProbabilityBJetTags[i];
      jet2[i].neutralEmEnergyFraction	= jet2_neutralEmEnergyFraction[i];
      jet2[i].neutralHadronEnergy	= jet2_neutralHadronEnergy[i];
      jet2[i].neutralHadronEnergyFraction	= jet2_neutralHadronEnergyFraction[i];
      jet2[i].numberOfDaughters	= jet2_numberOfDaughters[i];
      jet2[i].partonFlavour	= jet2_partonFlavour[i];
      jet2[i].phi	= jet2_phi[i];
      jet2[i].pt	= jet2_pt[i];
      jet2[i].simpleSecondaryVertexBJetTags	= jet2_simpleSecondaryVertexBJetTags[i];
      jet2[i].simpleSecondaryVertexHighEffBJetTags	= jet2_simpleSecondaryVertexHighEffBJetTags[i];
      jet2[i].simpleSecondaryVertexHighPurBJetTags	= jet2_simpleSecondaryVertexHighPurBJetTags[i];
      jet2[i].softElectronByIP3dBJetTags	= jet2_softElectronByIP3dBJetTags[i];
      jet2[i].softElectronByPtBJetTags	= jet2_softElectronByPtBJetTags[i];
      jet2[i].softMuonBJetTags	= jet2_softMuonBJetTags[i];
      jet2[i].softMuonByIP3dBJetTags	= jet2_softMuonByIP3dBJetTags[i];
      jet2[i].softMuonByPtBJetTags	= jet2_softMuonByPtBJetTags[i];
      jet2[i].trackCountingHighEffBJetTags	= jet2_trackCountingHighEffBJetTags[i];
      jet2[i].trackCountingHighPurBJetTags	= jet2_trackCountingHighPurBJetTags[i];
      jet2[i].uncor_energy	= jet2_uncor_energy[i];
      jet2[i].uncor_eta	= jet2_uncor_eta[i];
      jet2[i].uncor_phi	= jet2_uncor_phi[i];
      jet2[i].uncor_pt	= jet2_uncor_pt[i];
    }

  met.resize(met_energy.size());
  for(unsigned int i=0; i < met_energy.size(); ++i)
    {
      met[i].energy	= met_energy[i];
      met[i].et	= met_et[i];
      met[i].mEtSig	= met_mEtSig[i];
      met[i].phi	= met_phi[i];
      met[i].pt	= met_pt[i];
      met[i].sumEt	= met_sumEt[i];
    }

  met1.resize(met1_energy.size());
  for(unsigned int i=0; i < met1_energy.size(); ++i)
    {
      met1[i].energy	= met1_energy[i];
      met1[i].et	= met1_et[i];
      met1[i].mEtSig	= met1_mEtSig[i];
      met1[i].phi	= met1_phi[i];
      met1[i].pt	= met1_pt[i];
      met1[i].sumEt	= met1_sumEt[i];
    }

  met2.resize(met2_energy.size());
  for(unsigned int i=0; i < met2_energy.size(); ++i)
    {
      met2[i].energy	= met2_energy[i];
      met2[i].et	= met2_et[i];
      met2[i].mEtSig	= met2_mEtSig[i];
      met2[i].phi	= met2_phi[i];
      met2[i].pt	= met2_pt[i];
      met2[i].sumEt	= met2_sumEt[i];
    }

  muon.resize(muon_dxywrtBeamSpot.size());
  for(unsigned int i=0; i < muon_dxywrtBeamSpot.size(); ++i)
    {
      muon[i].dxywrtBeamSpot	= muon_dxywrtBeamSpot[i];
      muon[i].muonPFcandptdiff	= muon_muonPFcandptdiff[i];
    }

  muon1.resize(muon1_AllGlobalMuons.size());
  for(unsigned int i=0; i < muon1_AllGlobalMuons.size(); ++i)
    {
      muon1[i].AllGlobalMuons	= muon1_AllGlobalMuons[i];
      muon1[i].GlobalMuonPromptTight	= muon1_GlobalMuonPromptTight[i];
      muon1[i].charge	= muon1_charge[i];
      muon1[i].chargedHadronIso	= muon1_chargedHadronIso[i];
      muon1[i].combinedMuon_chi2	= muon1_combinedMuon_chi2[i];
      muon1[i].combinedMuon_ndof	= muon1_combinedMuon_ndof[i];
      muon1[i].combinedMuon_numberOfValidHits	= muon1_combinedMuon_numberOfValidHits[i];
      muon1[i].dB	= muon1_dB[i];
      muon1[i].ecalIso	= muon1_ecalIso[i];
      muon1[i].energy	= muon1_energy[i];
      muon1[i].eta	= muon1_eta[i];
      muon1[i].genLepton_eta	= muon1_genLepton_eta[i];
      muon1[i].genLepton_pdgId	= muon1_genLepton_pdgId[i];
      muon1[i].genLepton_phi	= muon1_genLepton_phi[i];
      muon1[i].genLepton_pt	= muon1_genLepton_pt[i];
      muon1[i].globalTrack_pt	= muon1_globalTrack_pt[i];
      muon1[i].hcalIso	= muon1_hcalIso[i];
      muon1[i].innerTrack_numberOfValidHits	= muon1_innerTrack_numberOfValidHits[i];
      muon1[i].innerTrack_pt	= muon1_innerTrack_pt[i];
      muon1[i].isGlobalMuon	= muon1_isGlobalMuon[i];
      muon1[i].isTrackerMuon	= muon1_isTrackerMuon[i];
      muon1[i].neutralHadronIso	= muon1_neutralHadronIso[i];
      muon1[i].phi	= muon1_phi[i];
      muon1[i].photonIso	= muon1_photonIso[i];
      muon1[i].pt	= muon1_pt[i];
      muon1[i].px	= muon1_px[i];
      muon1[i].py	= muon1_py[i];
      muon1[i].pz	= muon1_pz[i];
      muon1[i].trackIso	= muon1_trackIso[i];
      muon1[i].track_d0	= muon1_track_d0[i];
      muon1[i].track_hitPattern_numberOfValidMuonHits	= muon1_track_hitPattern_numberOfValidMuonHits[i];
      muon1[i].track_hitPattern_numberOfValidPixelHits	= muon1_track_hitPattern_numberOfValidPixelHits[i];
      muon1[i].track_normalizedChi2	= muon1_track_normalizedChi2[i];
      muon1[i].vz	= muon1_vz[i];
    }

  tau.resize(tau_genTauDecayModeID.size());
  for(unsigned int i=0; i < tau_genTauDecayModeID.size(); ++i)
    {
      tau[i].genTauDecayModeID	= tau_genTauDecayModeID[i];
    }

  tau1.resize(tau1_caloIso.size());
  for(unsigned int i=0; i < tau1_caloIso.size(); ++i)
    {
      tau1[i].caloIso	= tau1_caloIso[i];
      tau1[i].ecalIso	= tau1_ecalIso[i];
      tau1[i].energy	= tau1_energy[i];
      tau1[i].et	= tau1_et[i];
      tau1[i].eta	= tau1_eta[i];
      tau1[i].hcalIso	= tau1_hcalIso[i];
      tau1[i].phi	= tau1_phi[i];
      tau1[i].pt	= tau1_pt[i];
      tau1[i].px	= tau1_px[i];
      tau1[i].py	= tau1_py[i];
      tau1[i].pz	= tau1_pz[i];
      tau1[i].tauID_againstElectron	= tau1_tauID_againstElectron[i];
      tau1[i].tauID_againstMuon	= tau1_tauID_againstMuon[i];
      tau1[i].tauID_byIsolation	= tau1_tauID_byIsolation[i];
      tau1[i].tauID_byTaNC	= tau1_tauID_byTaNC[i];
      tau1[i].tauID_byTaNCfrHalfPercent	= tau1_tauID_byTaNCfrHalfPercent[i];
      tau1[i].tauID_byTaNCfrQuarterPercent	= tau1_tauID_byTaNCfrQuarterPercent[i];
      tau1[i].trackIso	= tau1_trackIso[i];
    }

  vertex.resize(vertex_isFake.size());
  for(unsigned int i=0; i < vertex_isFake.size(); ++i)
    {
      vertex[i].isFake	= vertex_isFake[i];
      vertex[i].isValid	= vertex_isValid[i];
      vertex[i].ndof	= vertex_ndof[i];
      vertex[i].normalizedChi2	= vertex_normalizedChi2[i];
      vertex[i].position_Rho	= vertex_position_Rho[i];
      vertex[i].tracksSize	= vertex_tracksSize[i];
      vertex[i].x	= vertex_x[i];
      vertex[i].xError	= vertex_xError[i];
      vertex[i].y	= vertex_y[i];
      vertex[i].yError	= vertex_yError[i];
      vertex[i].z	= vertex_z[i];
      vertex[i].zError	= vertex_zError[i];
    }
}


#endif
