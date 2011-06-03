#ifndef BASICLOOPCU_H
#define BASICLOOPCU_H
//-----------------------------------------------------------------------------
// File:        BasicLoopCU.h
// Description: Analyzer header for ntuples created by TheNtupleMaker
// Created:     Fri Jun  3 16:38:31 2011 by mkntanalyzer.py
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
double	edmtriggerresults_HLT_HT260_MHT60_v2;
double	edmtriggerresults_HLT_HT260_MHT60_v2_prs;
double	edmtriggerresults_HLT_HT300_MHT75_v2;
double	edmtriggerresults_HLT_HT300_MHT75_v2_prs;
double	edmtriggerresults_HLT_HT300_MHT75_v3;
double	edmtriggerresults_HLT_HT300_MHT75_v3_prs;
double	edmtriggerresults_HLT_HT300_MHT75_v4;
double	edmtriggerresults_HLT_HT300_MHT75_v4_prs;
double	edmtriggerresults_HLT_HT300_v1;
double	edmtriggerresults_HLT_HT300_v1_prs;
double	edmtriggerresults_HLT_HT300_v2;
double	edmtriggerresults_HLT_HT300_v2_prs;
double	edmtriggerresults_HLT_HT300_v3;
double	edmtriggerresults_HLT_HT300_v3_prs;
double	edmtriggerresults_HLT_HT300_v4;
double	edmtriggerresults_HLT_HT300_v4_prs;
std::vector<double>	electron1_dxywrtBeamSpot(3,0);
std::vector<double>	electron2_caloIso(9,0);
std::vector<double>	electron2_charge(9,0);
std::vector<double>	electron2_chargedHadronIso(9,0);
std::vector<double>	electron2_dB(9,0);
std::vector<double>	electron2_deltaEtaSuperClusterTrackAtVtx(9,0);
std::vector<double>	electron2_deltaPhiSuperClusterTrackAtVtx(9,0);
std::vector<double>	electron2_dr03EcalRecHitSumEt(9,0);
std::vector<double>	electron2_dr03HcalTowerSumEt(9,0);
std::vector<double>	electron2_dr03TkSumPt(9,0);
std::vector<double>	electron2_ecalIso(9,0);
std::vector<double>	electron2_eidRobustTight(9,0);
std::vector<double>	electron2_energy(9,0);
std::vector<double>	electron2_et(9,0);
std::vector<double>	electron2_eta(9,0);
std::vector<double>	electron2_genLepton_eta(9,0);
std::vector<double>	electron2_genLepton_pdgId(9,0);
std::vector<double>	electron2_genLepton_phi(9,0);
std::vector<double>	electron2_genLepton_pt(9,0);
std::vector<double>	electron2_gsfTrack_d0(9,0);
std::vector<double>	electron2_gsfTrack_trackerExpectedHitsInner_numberOfLostHits(9,0);
std::vector<double>	electron2_hadronicOverEm(9,0);
std::vector<double>	electron2_hcalIso(9,0);
std::vector<double>	electron2_neutralHadronIso(9,0);
std::vector<double>	electron2_phi(9,0);
std::vector<double>	electron2_photonIso(9,0);
std::vector<double>	electron2_pt(9,0);
std::vector<double>	electron2_px(9,0);
std::vector<double>	electron2_py(9,0);
std::vector<double>	electron2_pz(9,0);
std::vector<double>	electron2_sigmaIetaIeta(9,0);
std::vector<double>	electron2_simpleEleId95relIso(9,0);
std::vector<double>	electron2_superCluster_eta(9,0);
std::vector<double>	electron2_trackIso(9,0);
std::vector<double>	electron2_vz(9,0);
std::vector<double>	electron3_caloIso(3,0);
std::vector<double>	electron3_charge(3,0);
std::vector<double>	electron3_chargedHadronIso(3,0);
std::vector<double>	electron3_dB(3,0);
std::vector<double>	electron3_deltaEtaSuperClusterTrackAtVtx(3,0);
std::vector<double>	electron3_deltaPhiSuperClusterTrackAtVtx(3,0);
std::vector<double>	electron3_dr03EcalRecHitSumEt(3,0);
std::vector<double>	electron3_dr03HcalTowerSumEt(3,0);
std::vector<double>	electron3_dr03TkSumPt(3,0);
std::vector<double>	electron3_ecalIso(3,0);
std::vector<double>	electron3_eidRobustTight(3,0);
std::vector<double>	electron3_energy(3,0);
std::vector<double>	electron3_et(3,0);
std::vector<double>	electron3_eta(3,0);
std::vector<double>	electron3_genLepton_eta(3,0);
std::vector<double>	electron3_genLepton_pdgId(3,0);
std::vector<double>	electron3_genLepton_phi(3,0);
std::vector<double>	electron3_genLepton_pt(3,0);
std::vector<double>	electron3_gsfTrack_d0(3,0);
std::vector<double>	electron3_gsfTrack_trackerExpectedHitsInner_numberOfLostHits(3,0);
std::vector<double>	electron3_hadronicOverEm(3,0);
std::vector<double>	electron3_hcalIso(3,0);
std::vector<double>	electron3_neutralHadronIso(3,0);
std::vector<double>	electron3_phi(3,0);
std::vector<double>	electron3_photonIso(3,0);
std::vector<double>	electron3_pt(3,0);
std::vector<double>	electron3_px(3,0);
std::vector<double>	electron3_py(3,0);
std::vector<double>	electron3_pz(3,0);
std::vector<double>	electron3_sigmaIetaIeta(3,0);
std::vector<double>	electron3_simpleEleId95relIso(3,0);
std::vector<double>	electron3_superCluster_eta(3,0);
std::vector<double>	electron3_trackIso(3,0);
std::vector<double>	electron3_vz(3,0);
std::vector<double>	electron_dxywrtBeamSpot(9,0);
double	geneventinfoproduct1_pdf1;
double	geneventinfoproduct1_pdf2;
double	geneventinfoproduct1_scalePDF;
double	geneventinfoproduct1_weight;
double	geneventinfoproduct1_x1;
double	geneventinfoproduct1_x2;
std::vector<double>	geneventinfoproduct_pdf1(91,0);
std::vector<double>	geneventinfoproduct_pdf2(91,0);
std::vector<double>	geneventinfoproduct_pdfweight(91,0);
std::vector<double>	geneventinfoproduct_pdfweightsum(91,0);
std::vector<double>	genparticlera2_charge(31,0);
std::vector<double>	genparticlera2_eta(31,0);
std::vector<double>	genparticlera2_firstDaughter(31,0);
std::vector<double>	genparticlera2_firstMother(31,0);
std::vector<double>	genparticlera2_lastDaughter(31,0);
std::vector<double>	genparticlera2_lastMother(31,0);
std::vector<double>	genparticlera2_mass(31,0);
std::vector<double>	genparticlera2_pdgId(31,0);
std::vector<double>	genparticlera2_phi(31,0);
std::vector<double>	genparticlera2_pt(31,0);
std::vector<double>	genparticlera2_status(31,0);
double	genruninfoproduct_externalXSecLO_error;
double	genruninfoproduct_externalXSecLO_value;
double	genruninfoproduct_filterEfficiency;
double	genruninfoproduct_internalXSec_error;
double	genruninfoproduct_internalXSec_value;
std::vector<double>	jet1_jetUncMinus(175,0);
std::vector<double>	jet1_jetUncPlus(175,0);
std::vector<double>	jet2_combinedSecondaryVertexBJetTags(167,0);
std::vector<double>	jet2_combinedSecondaryVertexMVABJetTags(167,0);
std::vector<double>	jet2_emEnergyFraction(167,0);
std::vector<double>	jet2_energy(167,0);
std::vector<double>	jet2_eta(167,0);
std::vector<double>	jet2_jetBProbabilityBJetTags(167,0);
std::vector<double>	jet2_jetID_fHPD(167,0);
std::vector<double>	jet2_jetID_n90Hits(167,0);
std::vector<double>	jet2_jetProbabilityBJetTags(167,0);
std::vector<double>	jet2_partonFlavour(167,0);
std::vector<double>	jet2_phi(167,0);
std::vector<double>	jet2_pt(167,0);
std::vector<double>	jet2_simpleSecondaryVertexBJetTags(167,0);
std::vector<double>	jet2_simpleSecondaryVertexHighEffBJetTags(167,0);
std::vector<double>	jet2_simpleSecondaryVertexHighPurBJetTags(167,0);
std::vector<double>	jet2_softElectronByIP3dBJetTags(167,0);
std::vector<double>	jet2_softElectronByPtBJetTags(167,0);
std::vector<double>	jet2_softMuonBJetTags(167,0);
std::vector<double>	jet2_softMuonByIP3dBJetTags(167,0);
std::vector<double>	jet2_softMuonByPtBJetTags(167,0);
std::vector<double>	jet2_trackCountingHighEffBJetTags(167,0);
std::vector<double>	jet2_trackCountingHighPurBJetTags(167,0);
std::vector<double>	jet2_uncor_energy(167,0);
std::vector<double>	jet2_uncor_eta(167,0);
std::vector<double>	jet2_uncor_phi(167,0);
std::vector<double>	jet2_uncor_pt(167,0);
std::vector<double>	jet3_HFHadronEnergy(175,0);
std::vector<double>	jet3_chargedEmEnergyFraction(175,0);
std::vector<double>	jet3_chargedHadronEnergyFraction(175,0);
std::vector<double>	jet3_chargedMultiplicity(175,0);
std::vector<double>	jet3_combinedSecondaryVertexBJetTags(175,0);
std::vector<double>	jet3_combinedSecondaryVertexMVABJetTags(175,0);
std::vector<double>	jet3_energy(175,0);
std::vector<double>	jet3_eta(175,0);
std::vector<double>	jet3_jetBProbabilityBJetTags(175,0);
std::vector<double>	jet3_jetID_fHPD(175,0);
std::vector<double>	jet3_jetID_n90Hits(175,0);
std::vector<double>	jet3_jetProbabilityBJetTags(175,0);
std::vector<double>	jet3_neutralEmEnergyFraction(175,0);
std::vector<double>	jet3_neutralHadronEnergy(175,0);
std::vector<double>	jet3_neutralHadronEnergyFraction(175,0);
std::vector<double>	jet3_numberOfDaughters(175,0);
std::vector<double>	jet3_partonFlavour(175,0);
std::vector<double>	jet3_phi(175,0);
std::vector<double>	jet3_pt(175,0);
std::vector<double>	jet3_simpleSecondaryVertexBJetTags(175,0);
std::vector<double>	jet3_simpleSecondaryVertexHighEffBJetTags(175,0);
std::vector<double>	jet3_simpleSecondaryVertexHighPurBJetTags(175,0);
std::vector<double>	jet3_softElectronByIP3dBJetTags(175,0);
std::vector<double>	jet3_softElectronByPtBJetTags(175,0);
std::vector<double>	jet3_softMuonBJetTags(175,0);
std::vector<double>	jet3_softMuonByIP3dBJetTags(175,0);
std::vector<double>	jet3_softMuonByPtBJetTags(175,0);
std::vector<double>	jet3_trackCountingHighEffBJetTags(175,0);
std::vector<double>	jet3_trackCountingHighPurBJetTags(175,0);
std::vector<double>	jet3_uncor_energy(175,0);
std::vector<double>	jet3_uncor_eta(175,0);
std::vector<double>	jet3_uncor_phi(175,0);
std::vector<double>	jet3_uncor_pt(175,0);
std::vector<double>	jet_jetUncMinus(167,0);
std::vector<double>	jet_jetUncPlus(167,0);
std::vector<double>	met1_energy(3,0);
std::vector<double>	met1_et(3,0);
std::vector<double>	met1_mEtSig(3,0);
std::vector<double>	met1_phi(3,0);
std::vector<double>	met1_pt(3,0);
std::vector<double>	met1_sumEt(3,0);
std::vector<double>	met_energy(3,0);
std::vector<double>	met_et(3,0);
std::vector<double>	met_mEtSig(3,0);
std::vector<double>	met_phi(3,0);
std::vector<double>	met_pt(3,0);
std::vector<double>	met_sumEt(3,0);
std::vector<double>	muon1_dxywrtBeamSpot(3,0);
std::vector<double>	muon2_AllGlobalMuons(9,0);
std::vector<double>	muon2_GlobalMuonPromptTight(9,0);
std::vector<double>	muon2_charge(9,0);
std::vector<double>	muon2_chargedHadronIso(9,0);
std::vector<double>	muon2_combinedMuon_chi2(9,0);
std::vector<double>	muon2_combinedMuon_ndof(9,0);
std::vector<double>	muon2_combinedMuon_numberOfValidHits(9,0);
std::vector<double>	muon2_dB(9,0);
std::vector<double>	muon2_ecalIso(9,0);
std::vector<double>	muon2_energy(9,0);
std::vector<double>	muon2_eta(9,0);
std::vector<double>	muon2_genLepton_eta(9,0);
std::vector<double>	muon2_genLepton_pdgId(9,0);
std::vector<double>	muon2_genLepton_phi(9,0);
std::vector<double>	muon2_genLepton_pt(9,0);
std::vector<double>	muon2_globalTrack_pt(9,0);
std::vector<double>	muon2_hcalIso(9,0);
std::vector<double>	muon2_innerTrack_numberOfValidHits(9,0);
std::vector<double>	muon2_innerTrack_pt(9,0);
std::vector<double>	muon2_isGlobalMuon(9,0);
std::vector<double>	muon2_isTrackerMuon(9,0);
std::vector<double>	muon2_neutralHadronIso(9,0);
std::vector<double>	muon2_phi(9,0);
std::vector<double>	muon2_photonIso(9,0);
std::vector<double>	muon2_pt(9,0);
std::vector<double>	muon2_px(9,0);
std::vector<double>	muon2_py(9,0);
std::vector<double>	muon2_pz(9,0);
std::vector<double>	muon2_trackIso(9,0);
std::vector<double>	muon2_track_d0(9,0);
std::vector<double>	muon2_track_hitPattern_numberOfValidMuonHits(9,0);
std::vector<double>	muon2_track_hitPattern_numberOfValidPixelHits(9,0);
std::vector<double>	muon2_track_normalizedChi2(9,0);
std::vector<double>	muon2_vz(9,0);
std::vector<double>	muon3_AllGlobalMuons(3,0);
std::vector<double>	muon3_GlobalMuonPromptTight(3,0);
std::vector<double>	muon3_charge(3,0);
std::vector<double>	muon3_chargedHadronIso(3,0);
std::vector<double>	muon3_combinedMuon_chi2(3,0);
std::vector<double>	muon3_combinedMuon_ndof(3,0);
std::vector<double>	muon3_combinedMuon_numberOfValidHits(3,0);
std::vector<double>	muon3_dB(3,0);
std::vector<double>	muon3_ecalIso(3,0);
std::vector<double>	muon3_energy(3,0);
std::vector<double>	muon3_eta(3,0);
std::vector<double>	muon3_genLepton_eta(3,0);
std::vector<double>	muon3_genLepton_pdgId(3,0);
std::vector<double>	muon3_genLepton_phi(3,0);
std::vector<double>	muon3_genLepton_pt(3,0);
std::vector<double>	muon3_globalTrack_pt(3,0);
std::vector<double>	muon3_hcalIso(3,0);
std::vector<double>	muon3_innerTrack_numberOfValidHits(3,0);
std::vector<double>	muon3_innerTrack_pt(3,0);
std::vector<double>	muon3_isGlobalMuon(3,0);
std::vector<double>	muon3_isTrackerMuon(3,0);
std::vector<double>	muon3_neutralHadronIso(3,0);
std::vector<double>	muon3_phi(3,0);
std::vector<double>	muon3_photonIso(3,0);
std::vector<double>	muon3_pt(3,0);
std::vector<double>	muon3_px(3,0);
std::vector<double>	muon3_py(3,0);
std::vector<double>	muon3_pz(3,0);
std::vector<double>	muon3_trackIso(3,0);
std::vector<double>	muon3_track_d0(3,0);
std::vector<double>	muon3_track_hitPattern_numberOfValidMuonHits(3,0);
std::vector<double>	muon3_track_hitPattern_numberOfValidPixelHits(3,0);
std::vector<double>	muon3_track_normalizedChi2(3,0);
std::vector<double>	muon3_vz(3,0);
std::vector<double>	muon_dxywrtBeamSpot(9,0);
int	nGenEventInfoProductHelper_generator;
int	nelectron;
int	nelectron1;
int	nelectron2;
int	nelectron3;
int	ngenparticlera2b;
int	njet;
int	njet1;
int	njet2;
int	njet3;
int	nmet;
int	nmet1;
int	nmuon;
int	nmuon1;
int	nmuon2;
int	nmuon3;
int	ntau;
int	ntau1;
int	nvertex;
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
std::vector<double>	vertex_isFake(39,0);
std::vector<double>	vertex_isValid(39,0);
std::vector<double>	vertex_ndof(39,0);
std::vector<double>	vertex_normalizedChi2(39,0);
std::vector<double>	vertex_position_Rho(39,0);
std::vector<double>	vertex_tracksSize(39,0);
std::vector<double>	vertex_x(39,0);
std::vector<double>	vertex_xError(39,0);
std::vector<double>	vertex_y(39,0);
std::vector<double>	vertex_yError(39,0);
std::vector<double>	vertex_z(39,0);
std::vector<double>	vertex_zError(39,0);


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
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT100U", edmtriggerresults_HLT_HT100U);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT100U_prs", edmtriggerresults_HLT_HT100U_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT120U", edmtriggerresults_HLT_HT120U);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT120U_prs", edmtriggerresults_HLT_HT120U_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT140U", edmtriggerresults_HLT_HT140U);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT140U_prs", edmtriggerresults_HLT_HT140U_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT150U_v3", edmtriggerresults_HLT_HT150U_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT150U_v3_prs", edmtriggerresults_HLT_HT150U_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT150_v1", edmtriggerresults_HLT_HT150_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT150_v1_prs", edmtriggerresults_HLT_HT150_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT150_v2", edmtriggerresults_HLT_HT150_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT150_v2_prs", edmtriggerresults_HLT_HT150_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT150_v3", edmtriggerresults_HLT_HT150_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT150_v3_prs", edmtriggerresults_HLT_HT150_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT250_MHT60_v2", edmtriggerresults_HLT_HT250_MHT60_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT250_MHT60_v2_prs", edmtriggerresults_HLT_HT250_MHT60_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT250_MHT60_v3", edmtriggerresults_HLT_HT250_MHT60_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT250_MHT60_v3_prs", edmtriggerresults_HLT_HT250_MHT60_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT250_v1", edmtriggerresults_HLT_HT250_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT250_v1_prs", edmtriggerresults_HLT_HT250_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT250_v2", edmtriggerresults_HLT_HT250_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT250_v2_prs", edmtriggerresults_HLT_HT250_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT250_v3", edmtriggerresults_HLT_HT250_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT250_v3_prs", edmtriggerresults_HLT_HT250_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT260_MHT60_v2", edmtriggerresults_HLT_HT260_MHT60_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT260_MHT60_v2_prs", edmtriggerresults_HLT_HT260_MHT60_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT300_MHT75_v2", edmtriggerresults_HLT_HT300_MHT75_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT300_MHT75_v2_prs", edmtriggerresults_HLT_HT300_MHT75_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT300_MHT75_v3", edmtriggerresults_HLT_HT300_MHT75_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT300_MHT75_v3_prs", edmtriggerresults_HLT_HT300_MHT75_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT300_MHT75_v4", edmtriggerresults_HLT_HT300_MHT75_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT300_MHT75_v4_prs", edmtriggerresults_HLT_HT300_MHT75_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT300_v1", edmtriggerresults_HLT_HT300_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT300_v1_prs", edmtriggerresults_HLT_HT300_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT300_v2", edmtriggerresults_HLT_HT300_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT300_v2_prs", edmtriggerresults_HLT_HT300_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT300_v3", edmtriggerresults_HLT_HT300_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT300_v3_prs", edmtriggerresults_HLT_HT300_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT300_v4", edmtriggerresults_HLT_HT300_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults.HLT_HT300_v4_prs", edmtriggerresults_HLT_HT300_v4_prs);
  stream.select("patElectronHelper_selectedPatElectronsPF.dxywrtBeamSpot", electron1_dxywrtBeamSpot);
  stream.select("patElectron_patElectronsNoIsoPF.caloIso", electron2_caloIso);
  stream.select("patElectron_patElectronsNoIsoPF.charge", electron2_charge);
  stream.select("patElectron_patElectronsNoIsoPF.chargedHadronIso", electron2_chargedHadronIso);
  stream.select("patElectron_patElectronsNoIsoPF.dB", electron2_dB);
  stream.select("patElectron_patElectronsNoIsoPF.deltaEtaSuperClusterTrackAtVtx", electron2_deltaEtaSuperClusterTrackAtVtx);
  stream.select("patElectron_patElectronsNoIsoPF.deltaPhiSuperClusterTrackAtVtx", electron2_deltaPhiSuperClusterTrackAtVtx);
  stream.select("patElectron_patElectronsNoIsoPF.dr03EcalRecHitSumEt", electron2_dr03EcalRecHitSumEt);
  stream.select("patElectron_patElectronsNoIsoPF.dr03HcalTowerSumEt", electron2_dr03HcalTowerSumEt);
  stream.select("patElectron_patElectronsNoIsoPF.dr03TkSumPt", electron2_dr03TkSumPt);
  stream.select("patElectron_patElectronsNoIsoPF.ecalIso", electron2_ecalIso);
  stream.select("patElectron_patElectronsNoIsoPF.eidRobustTight", electron2_eidRobustTight);
  stream.select("patElectron_patElectronsNoIsoPF.energy", electron2_energy);
  stream.select("patElectron_patElectronsNoIsoPF.et", electron2_et);
  stream.select("patElectron_patElectronsNoIsoPF.eta", electron2_eta);
  stream.select("patElectron_patElectronsNoIsoPF.genLepton_eta", electron2_genLepton_eta);
  stream.select("patElectron_patElectronsNoIsoPF.genLepton_pdgId", electron2_genLepton_pdgId);
  stream.select("patElectron_patElectronsNoIsoPF.genLepton_phi", electron2_genLepton_phi);
  stream.select("patElectron_patElectronsNoIsoPF.genLepton_pt", electron2_genLepton_pt);
  stream.select("patElectron_patElectronsNoIsoPF.gsfTrack_d0", electron2_gsfTrack_d0);
  stream.select("patElectron_patElectronsNoIsoPF.gsfTrack_trackerExpectedHitsInner_numberOfLostHits", electron2_gsfTrack_trackerExpectedHitsInner_numberOfLostHits);
  stream.select("patElectron_patElectronsNoIsoPF.hadronicOverEm", electron2_hadronicOverEm);
  stream.select("patElectron_patElectronsNoIsoPF.hcalIso", electron2_hcalIso);
  stream.select("patElectron_patElectronsNoIsoPF.neutralHadronIso", electron2_neutralHadronIso);
  stream.select("patElectron_patElectronsNoIsoPF.phi", electron2_phi);
  stream.select("patElectron_patElectronsNoIsoPF.photonIso", electron2_photonIso);
  stream.select("patElectron_patElectronsNoIsoPF.pt", electron2_pt);
  stream.select("patElectron_patElectronsNoIsoPF.px", electron2_px);
  stream.select("patElectron_patElectronsNoIsoPF.py", electron2_py);
  stream.select("patElectron_patElectronsNoIsoPF.pz", electron2_pz);
  stream.select("patElectron_patElectronsNoIsoPF.sigmaIetaIeta", electron2_sigmaIetaIeta);
  stream.select("patElectron_patElectronsNoIsoPF.simpleEleId95relIso", electron2_simpleEleId95relIso);
  stream.select("patElectron_patElectronsNoIsoPF.superCluster_eta", electron2_superCluster_eta);
  stream.select("patElectron_patElectronsNoIsoPF.trackIso", electron2_trackIso);
  stream.select("patElectron_patElectronsNoIsoPF.vz", electron2_vz);
  stream.select("patElectron_selectedPatElectronsPF.caloIso", electron3_caloIso);
  stream.select("patElectron_selectedPatElectronsPF.charge", electron3_charge);
  stream.select("patElectron_selectedPatElectronsPF.chargedHadronIso", electron3_chargedHadronIso);
  stream.select("patElectron_selectedPatElectronsPF.dB", electron3_dB);
  stream.select("patElectron_selectedPatElectronsPF.deltaEtaSuperClusterTrackAtVtx", electron3_deltaEtaSuperClusterTrackAtVtx);
  stream.select("patElectron_selectedPatElectronsPF.deltaPhiSuperClusterTrackAtVtx", electron3_deltaPhiSuperClusterTrackAtVtx);
  stream.select("patElectron_selectedPatElectronsPF.dr03EcalRecHitSumEt", electron3_dr03EcalRecHitSumEt);
  stream.select("patElectron_selectedPatElectronsPF.dr03HcalTowerSumEt", electron3_dr03HcalTowerSumEt);
  stream.select("patElectron_selectedPatElectronsPF.dr03TkSumPt", electron3_dr03TkSumPt);
  stream.select("patElectron_selectedPatElectronsPF.ecalIso", electron3_ecalIso);
  stream.select("patElectron_selectedPatElectronsPF.eidRobustTight", electron3_eidRobustTight);
  stream.select("patElectron_selectedPatElectronsPF.energy", electron3_energy);
  stream.select("patElectron_selectedPatElectronsPF.et", electron3_et);
  stream.select("patElectron_selectedPatElectronsPF.eta", electron3_eta);
  stream.select("patElectron_selectedPatElectronsPF.genLepton_eta", electron3_genLepton_eta);
  stream.select("patElectron_selectedPatElectronsPF.genLepton_pdgId", electron3_genLepton_pdgId);
  stream.select("patElectron_selectedPatElectronsPF.genLepton_phi", electron3_genLepton_phi);
  stream.select("patElectron_selectedPatElectronsPF.genLepton_pt", electron3_genLepton_pt);
  stream.select("patElectron_selectedPatElectronsPF.gsfTrack_d0", electron3_gsfTrack_d0);
  stream.select("patElectron_selectedPatElectronsPF.gsfTrack_trackerExpectedHitsInner_numberOfLostHits", electron3_gsfTrack_trackerExpectedHitsInner_numberOfLostHits);
  stream.select("patElectron_selectedPatElectronsPF.hadronicOverEm", electron3_hadronicOverEm);
  stream.select("patElectron_selectedPatElectronsPF.hcalIso", electron3_hcalIso);
  stream.select("patElectron_selectedPatElectronsPF.neutralHadronIso", electron3_neutralHadronIso);
  stream.select("patElectron_selectedPatElectronsPF.phi", electron3_phi);
  stream.select("patElectron_selectedPatElectronsPF.photonIso", electron3_photonIso);
  stream.select("patElectron_selectedPatElectronsPF.pt", electron3_pt);
  stream.select("patElectron_selectedPatElectronsPF.px", electron3_px);
  stream.select("patElectron_selectedPatElectronsPF.py", electron3_py);
  stream.select("patElectron_selectedPatElectronsPF.pz", electron3_pz);
  stream.select("patElectron_selectedPatElectronsPF.sigmaIetaIeta", electron3_sigmaIetaIeta);
  stream.select("patElectron_selectedPatElectronsPF.simpleEleId95relIso", electron3_simpleEleId95relIso);
  stream.select("patElectron_selectedPatElectronsPF.superCluster_eta", electron3_superCluster_eta);
  stream.select("patElectron_selectedPatElectronsPF.trackIso", electron3_trackIso);
  stream.select("patElectron_selectedPatElectronsPF.vz", electron3_vz);
  stream.select("patElectronHelper_patElectronsNoIsoPF.dxywrtBeamSpot", electron_dxywrtBeamSpot);
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
  stream.select("patJetHelper_selectedPatJetsPF.jetUncMinus", jet1_jetUncMinus);
  stream.select("patJetHelper_selectedPatJetsPF.jetUncPlus", jet1_jetUncPlus);
  stream.select("patJet_cleanPatJetsAK5Calo.combinedSecondaryVertexBJetTags", jet2_combinedSecondaryVertexBJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.combinedSecondaryVertexMVABJetTags", jet2_combinedSecondaryVertexMVABJetTags);
  stream.select("patJet_cleanPatJetsAK5Calo.emEnergyFraction", jet2_emEnergyFraction);
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
  stream.select("patJet_selectedPatJetsPF.HFHadronEnergy", jet3_HFHadronEnergy);
  stream.select("patJet_selectedPatJetsPF.chargedEmEnergyFraction", jet3_chargedEmEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.chargedHadronEnergyFraction", jet3_chargedHadronEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.chargedMultiplicity", jet3_chargedMultiplicity);
  stream.select("patJet_selectedPatJetsPF.combinedSecondaryVertexBJetTags", jet3_combinedSecondaryVertexBJetTags);
  stream.select("patJet_selectedPatJetsPF.combinedSecondaryVertexMVABJetTags", jet3_combinedSecondaryVertexMVABJetTags);
  stream.select("patJet_selectedPatJetsPF.energy", jet3_energy);
  stream.select("patJet_selectedPatJetsPF.eta", jet3_eta);
  stream.select("patJet_selectedPatJetsPF.jetBProbabilityBJetTags", jet3_jetBProbabilityBJetTags);
  stream.select("patJet_selectedPatJetsPF.jetID_fHPD", jet3_jetID_fHPD);
  stream.select("patJet_selectedPatJetsPF.jetID_n90Hits", jet3_jetID_n90Hits);
  stream.select("patJet_selectedPatJetsPF.jetProbabilityBJetTags", jet3_jetProbabilityBJetTags);
  stream.select("patJet_selectedPatJetsPF.neutralEmEnergyFraction", jet3_neutralEmEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.neutralHadronEnergy", jet3_neutralHadronEnergy);
  stream.select("patJet_selectedPatJetsPF.neutralHadronEnergyFraction", jet3_neutralHadronEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.numberOfDaughters", jet3_numberOfDaughters);
  stream.select("patJet_selectedPatJetsPF.partonFlavour", jet3_partonFlavour);
  stream.select("patJet_selectedPatJetsPF.phi", jet3_phi);
  stream.select("patJet_selectedPatJetsPF.pt", jet3_pt);
  stream.select("patJet_selectedPatJetsPF.simpleSecondaryVertexBJetTags", jet3_simpleSecondaryVertexBJetTags);
  stream.select("patJet_selectedPatJetsPF.simpleSecondaryVertexHighEffBJetTags", jet3_simpleSecondaryVertexHighEffBJetTags);
  stream.select("patJet_selectedPatJetsPF.simpleSecondaryVertexHighPurBJetTags", jet3_simpleSecondaryVertexHighPurBJetTags);
  stream.select("patJet_selectedPatJetsPF.softElectronByIP3dBJetTags", jet3_softElectronByIP3dBJetTags);
  stream.select("patJet_selectedPatJetsPF.softElectronByPtBJetTags", jet3_softElectronByPtBJetTags);
  stream.select("patJet_selectedPatJetsPF.softMuonBJetTags", jet3_softMuonBJetTags);
  stream.select("patJet_selectedPatJetsPF.softMuonByIP3dBJetTags", jet3_softMuonByIP3dBJetTags);
  stream.select("patJet_selectedPatJetsPF.softMuonByPtBJetTags", jet3_softMuonByPtBJetTags);
  stream.select("patJet_selectedPatJetsPF.trackCountingHighEffBJetTags", jet3_trackCountingHighEffBJetTags);
  stream.select("patJet_selectedPatJetsPF.trackCountingHighPurBJetTags", jet3_trackCountingHighPurBJetTags);
  stream.select("patJet_selectedPatJetsPF.uncor_energy", jet3_uncor_energy);
  stream.select("patJet_selectedPatJetsPF.uncor_eta", jet3_uncor_eta);
  stream.select("patJet_selectedPatJetsPF.uncor_phi", jet3_uncor_phi);
  stream.select("patJet_selectedPatJetsPF.uncor_pt", jet3_uncor_pt);
  stream.select("patJetHelper_cleanPatJetsAK5Calo.jetUncMinus", jet_jetUncMinus);
  stream.select("patJetHelper_cleanPatJetsAK5Calo.jetUncPlus", jet_jetUncPlus);
  stream.select("patMET_patMETsPF.energy", met1_energy);
  stream.select("patMET_patMETsPF.et", met1_et);
  stream.select("patMET_patMETsPF.mEtSig", met1_mEtSig);
  stream.select("patMET_patMETsPF.phi", met1_phi);
  stream.select("patMET_patMETsPF.pt", met1_pt);
  stream.select("patMET_patMETsPF.sumEt", met1_sumEt);
  stream.select("patMET_patMETsAK5Calo.energy", met_energy);
  stream.select("patMET_patMETsAK5Calo.et", met_et);
  stream.select("patMET_patMETsAK5Calo.mEtSig", met_mEtSig);
  stream.select("patMET_patMETsAK5Calo.phi", met_phi);
  stream.select("patMET_patMETsAK5Calo.pt", met_pt);
  stream.select("patMET_patMETsAK5Calo.sumEt", met_sumEt);
  stream.select("patMuonHelper_selectedPatMuonsPF.dxywrtBeamSpot", muon1_dxywrtBeamSpot);
  stream.select("patMuon_patMuonsNoIsoPF.AllGlobalMuons", muon2_AllGlobalMuons);
  stream.select("patMuon_patMuonsNoIsoPF.GlobalMuonPromptTight", muon2_GlobalMuonPromptTight);
  stream.select("patMuon_patMuonsNoIsoPF.charge", muon2_charge);
  stream.select("patMuon_patMuonsNoIsoPF.chargedHadronIso", muon2_chargedHadronIso);
  stream.select("patMuon_patMuonsNoIsoPF.combinedMuon_chi2", muon2_combinedMuon_chi2);
  stream.select("patMuon_patMuonsNoIsoPF.combinedMuon_ndof", muon2_combinedMuon_ndof);
  stream.select("patMuon_patMuonsNoIsoPF.combinedMuon_numberOfValidHits", muon2_combinedMuon_numberOfValidHits);
  stream.select("patMuon_patMuonsNoIsoPF.dB", muon2_dB);
  stream.select("patMuon_patMuonsNoIsoPF.ecalIso", muon2_ecalIso);
  stream.select("patMuon_patMuonsNoIsoPF.energy", muon2_energy);
  stream.select("patMuon_patMuonsNoIsoPF.eta", muon2_eta);
  stream.select("patMuon_patMuonsNoIsoPF.genLepton_eta", muon2_genLepton_eta);
  stream.select("patMuon_patMuonsNoIsoPF.genLepton_pdgId", muon2_genLepton_pdgId);
  stream.select("patMuon_patMuonsNoIsoPF.genLepton_phi", muon2_genLepton_phi);
  stream.select("patMuon_patMuonsNoIsoPF.genLepton_pt", muon2_genLepton_pt);
  stream.select("patMuon_patMuonsNoIsoPF.globalTrack_pt", muon2_globalTrack_pt);
  stream.select("patMuon_patMuonsNoIsoPF.hcalIso", muon2_hcalIso);
  stream.select("patMuon_patMuonsNoIsoPF.innerTrack_numberOfValidHits", muon2_innerTrack_numberOfValidHits);
  stream.select("patMuon_patMuonsNoIsoPF.innerTrack_pt", muon2_innerTrack_pt);
  stream.select("patMuon_patMuonsNoIsoPF.isGlobalMuon", muon2_isGlobalMuon);
  stream.select("patMuon_patMuonsNoIsoPF.isTrackerMuon", muon2_isTrackerMuon);
  stream.select("patMuon_patMuonsNoIsoPF.neutralHadronIso", muon2_neutralHadronIso);
  stream.select("patMuon_patMuonsNoIsoPF.phi", muon2_phi);
  stream.select("patMuon_patMuonsNoIsoPF.photonIso", muon2_photonIso);
  stream.select("patMuon_patMuonsNoIsoPF.pt", muon2_pt);
  stream.select("patMuon_patMuonsNoIsoPF.px", muon2_px);
  stream.select("patMuon_patMuonsNoIsoPF.py", muon2_py);
  stream.select("patMuon_patMuonsNoIsoPF.pz", muon2_pz);
  stream.select("patMuon_patMuonsNoIsoPF.trackIso", muon2_trackIso);
  stream.select("patMuon_patMuonsNoIsoPF.track_d0", muon2_track_d0);
  stream.select("patMuon_patMuonsNoIsoPF.track_hitPattern_numberOfValidMuonHits", muon2_track_hitPattern_numberOfValidMuonHits);
  stream.select("patMuon_patMuonsNoIsoPF.track_hitPattern_numberOfValidPixelHits", muon2_track_hitPattern_numberOfValidPixelHits);
  stream.select("patMuon_patMuonsNoIsoPF.track_normalizedChi2", muon2_track_normalizedChi2);
  stream.select("patMuon_patMuonsNoIsoPF.vz", muon2_vz);
  stream.select("patMuon_selectedPatMuonsPF.AllGlobalMuons", muon3_AllGlobalMuons);
  stream.select("patMuon_selectedPatMuonsPF.GlobalMuonPromptTight", muon3_GlobalMuonPromptTight);
  stream.select("patMuon_selectedPatMuonsPF.charge", muon3_charge);
  stream.select("patMuon_selectedPatMuonsPF.chargedHadronIso", muon3_chargedHadronIso);
  stream.select("patMuon_selectedPatMuonsPF.combinedMuon_chi2", muon3_combinedMuon_chi2);
  stream.select("patMuon_selectedPatMuonsPF.combinedMuon_ndof", muon3_combinedMuon_ndof);
  stream.select("patMuon_selectedPatMuonsPF.combinedMuon_numberOfValidHits", muon3_combinedMuon_numberOfValidHits);
  stream.select("patMuon_selectedPatMuonsPF.dB", muon3_dB);
  stream.select("patMuon_selectedPatMuonsPF.ecalIso", muon3_ecalIso);
  stream.select("patMuon_selectedPatMuonsPF.energy", muon3_energy);
  stream.select("patMuon_selectedPatMuonsPF.eta", muon3_eta);
  stream.select("patMuon_selectedPatMuonsPF.genLepton_eta", muon3_genLepton_eta);
  stream.select("patMuon_selectedPatMuonsPF.genLepton_pdgId", muon3_genLepton_pdgId);
  stream.select("patMuon_selectedPatMuonsPF.genLepton_phi", muon3_genLepton_phi);
  stream.select("patMuon_selectedPatMuonsPF.genLepton_pt", muon3_genLepton_pt);
  stream.select("patMuon_selectedPatMuonsPF.globalTrack_pt", muon3_globalTrack_pt);
  stream.select("patMuon_selectedPatMuonsPF.hcalIso", muon3_hcalIso);
  stream.select("patMuon_selectedPatMuonsPF.innerTrack_numberOfValidHits", muon3_innerTrack_numberOfValidHits);
  stream.select("patMuon_selectedPatMuonsPF.innerTrack_pt", muon3_innerTrack_pt);
  stream.select("patMuon_selectedPatMuonsPF.isGlobalMuon", muon3_isGlobalMuon);
  stream.select("patMuon_selectedPatMuonsPF.isTrackerMuon", muon3_isTrackerMuon);
  stream.select("patMuon_selectedPatMuonsPF.neutralHadronIso", muon3_neutralHadronIso);
  stream.select("patMuon_selectedPatMuonsPF.phi", muon3_phi);
  stream.select("patMuon_selectedPatMuonsPF.photonIso", muon3_photonIso);
  stream.select("patMuon_selectedPatMuonsPF.pt", muon3_pt);
  stream.select("patMuon_selectedPatMuonsPF.px", muon3_px);
  stream.select("patMuon_selectedPatMuonsPF.py", muon3_py);
  stream.select("patMuon_selectedPatMuonsPF.pz", muon3_pz);
  stream.select("patMuon_selectedPatMuonsPF.trackIso", muon3_trackIso);
  stream.select("patMuon_selectedPatMuonsPF.track_d0", muon3_track_d0);
  stream.select("patMuon_selectedPatMuonsPF.track_hitPattern_numberOfValidMuonHits", muon3_track_hitPattern_numberOfValidMuonHits);
  stream.select("patMuon_selectedPatMuonsPF.track_hitPattern_numberOfValidPixelHits", muon3_track_hitPattern_numberOfValidPixelHits);
  stream.select("patMuon_selectedPatMuonsPF.track_normalizedChi2", muon3_track_normalizedChi2);
  stream.select("patMuon_selectedPatMuonsPF.vz", muon3_vz);
  stream.select("patMuonHelper_patMuonsNoIsoPF.dxywrtBeamSpot", muon_dxywrtBeamSpot);
  stream.select("nGenEventInfoProductHelper_generator", nGenEventInfoProductHelper_generator);
  stream.select("npatElectronHelper_patElectronsNoIsoPF", nelectron);
  stream.select("npatElectronHelper_selectedPatElectronsPF", nelectron1);
  stream.select("npatElectron_patElectronsNoIsoPF", nelectron2);
  stream.select("npatElectron_selectedPatElectronsPF", nelectron3);
  stream.select("nrecoGenParticleHelperRA2b_genParticles", ngenparticlera2b);
  stream.select("npatJetHelper_cleanPatJetsAK5Calo", njet);
  stream.select("npatJetHelper_selectedPatJetsPF", njet1);
  stream.select("npatJet_cleanPatJetsAK5Calo", njet2);
  stream.select("npatJet_selectedPatJetsPF", njet3);
  stream.select("npatMET_patMETsAK5Calo", nmet);
  stream.select("npatMET_patMETsPF", nmet1);
  stream.select("npatMuonHelper_patMuonsNoIsoPF", nmuon);
  stream.select("npatMuonHelper_selectedPatMuonsPF", nmuon1);
  stream.select("npatMuon_patMuonsNoIsoPF", nmuon2);
  stream.select("npatMuon_selectedPatMuonsPF", nmuon3);
  stream.select("npatTauHelper_selectedPatTausPF", ntau);
  stream.select("npatTau_selectedPatTausPF", ntau1);
  stream.select("nrecoVertex_offlinePrimaryVertices", nvertex);
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
std::vector<electron_s> electron(9);

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
  double	dxywrtBeamSpot;
};
std::vector<electron1_s> electron1(3);

std::ostream& operator<<(std::ostream& os, const electron1_s& o)
{
  char r[1024];
  os << "electron1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "dxywrtBeamSpot", (double)o.dxywrtBeamSpot); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct electron2_s
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
std::vector<electron2_s> electron2(9);

std::ostream& operator<<(std::ostream& os, const electron2_s& o)
{
  char r[1024];
  os << "electron2" << std::endl;
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
struct electron3_s
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
std::vector<electron3_s> electron3(3);

std::ostream& operator<<(std::ostream& os, const electron3_s& o)
{
  char r[1024];
  os << "electron3" << std::endl;
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
struct geneventinfoproduct_s
{
  double	pdf1;
  double	pdf2;
  double	pdfweight;
  double	pdfweightsum;
};
std::vector<geneventinfoproduct_s> geneventinfoproduct(91);

std::ostream& operator<<(std::ostream& os, const geneventinfoproduct_s& o)
{
  char r[1024];
  os << "geneventinfoproduct" << std::endl;
  sprintf(r, "  %-32s: %f\n", "pdf1", (double)o.pdf1); os << r;
  sprintf(r, "  %-32s: %f\n", "pdf2", (double)o.pdf2); os << r;
  sprintf(r, "  %-32s: %f\n", "pdfweight", (double)o.pdfweight); os << r;
  sprintf(r, "  %-32s: %f\n", "pdfweightsum", (double)o.pdfweightsum); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct genparticlera2_s
{
  double	charge;
  double	eta;
  double	firstDaughter;
  double	firstMother;
  double	lastDaughter;
  double	lastMother;
  double	mass;
  double	pdgId;
  double	phi;
  double	pt;
  double	status;
};
std::vector<genparticlera2_s> genparticlera2(31);

std::ostream& operator<<(std::ostream& os, const genparticlera2_s& o)
{
  char r[1024];
  os << "genparticlera2" << std::endl;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "firstDaughter", (double)o.firstDaughter); os << r;
  sprintf(r, "  %-32s: %f\n", "firstMother", (double)o.firstMother); os << r;
  sprintf(r, "  %-32s: %f\n", "lastDaughter", (double)o.lastDaughter); os << r;
  sprintf(r, "  %-32s: %f\n", "lastMother", (double)o.lastMother); os << r;
  sprintf(r, "  %-32s: %f\n", "mass", (double)o.mass); os << r;
  sprintf(r, "  %-32s: %f\n", "pdgId", (double)o.pdgId); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "status", (double)o.status); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jet_s
{
  double	jetUncMinus;
  double	jetUncPlus;
};
std::vector<jet_s> jet(167);

std::ostream& operator<<(std::ostream& os, const jet_s& o)
{
  char r[1024];
  os << "jet" << std::endl;
  sprintf(r, "  %-32s: %f\n", "jetUncMinus", (double)o.jetUncMinus); os << r;
  sprintf(r, "  %-32s: %f\n", "jetUncPlus", (double)o.jetUncPlus); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jet1_s
{
  double	jetUncMinus;
  double	jetUncPlus;
};
std::vector<jet1_s> jet1(175);

std::ostream& operator<<(std::ostream& os, const jet1_s& o)
{
  char r[1024];
  os << "jet1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "jetUncMinus", (double)o.jetUncMinus); os << r;
  sprintf(r, "  %-32s: %f\n", "jetUncPlus", (double)o.jetUncPlus); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jet2_s
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
std::vector<jet2_s> jet2(167);

std::ostream& operator<<(std::ostream& os, const jet2_s& o)
{
  char r[1024];
  os << "jet2" << std::endl;
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
struct jet3_s
{
  double	HFHadronEnergy;
  double	chargedEmEnergyFraction;
  double	chargedHadronEnergyFraction;
  double	chargedMultiplicity;
  double	combinedSecondaryVertexBJetTags;
  double	combinedSecondaryVertexMVABJetTags;
  double	energy;
  double	eta;
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
std::vector<jet3_s> jet3(175);

std::ostream& operator<<(std::ostream& os, const jet3_s& o)
{
  char r[1024];
  os << "jet3" << std::endl;
  sprintf(r, "  %-32s: %f\n", "HFHadronEnergy", (double)o.HFHadronEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedEmEnergyFraction", (double)o.chargedEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronEnergyFraction", (double)o.chargedHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedMultiplicity", (double)o.chargedMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexBJetTags", (double)o.combinedSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexMVABJetTags", (double)o.combinedSecondaryVertexMVABJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
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
struct muon_s
{
  double	dxywrtBeamSpot;
};
std::vector<muon_s> muon(9);

std::ostream& operator<<(std::ostream& os, const muon_s& o)
{
  char r[1024];
  os << "muon" << std::endl;
  sprintf(r, "  %-32s: %f\n", "dxywrtBeamSpot", (double)o.dxywrtBeamSpot); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct muon1_s
{
  double	dxywrtBeamSpot;
};
std::vector<muon1_s> muon1(3);

std::ostream& operator<<(std::ostream& os, const muon1_s& o)
{
  char r[1024];
  os << "muon1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "dxywrtBeamSpot", (double)o.dxywrtBeamSpot); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct muon2_s
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
std::vector<muon2_s> muon2(9);

std::ostream& operator<<(std::ostream& os, const muon2_s& o)
{
  char r[1024];
  os << "muon2" << std::endl;
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
struct muon3_s
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
std::vector<muon3_s> muon3(3);

std::ostream& operator<<(std::ostream& os, const muon3_s& o)
{
  char r[1024];
  os << "muon3" << std::endl;
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
std::vector<vertex_s> vertex(39);

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

  electron1.resize(electron1_dxywrtBeamSpot.size());
  for(unsigned int i=0; i < electron1_dxywrtBeamSpot.size(); ++i)
    {
      electron1[i].dxywrtBeamSpot	= electron1_dxywrtBeamSpot[i];
    }

  electron2.resize(electron2_caloIso.size());
  for(unsigned int i=0; i < electron2_caloIso.size(); ++i)
    {
      electron2[i].caloIso	= electron2_caloIso[i];
      electron2[i].charge	= electron2_charge[i];
      electron2[i].chargedHadronIso	= electron2_chargedHadronIso[i];
      electron2[i].dB	= electron2_dB[i];
      electron2[i].deltaEtaSuperClusterTrackAtVtx	= electron2_deltaEtaSuperClusterTrackAtVtx[i];
      electron2[i].deltaPhiSuperClusterTrackAtVtx	= electron2_deltaPhiSuperClusterTrackAtVtx[i];
      electron2[i].dr03EcalRecHitSumEt	= electron2_dr03EcalRecHitSumEt[i];
      electron2[i].dr03HcalTowerSumEt	= electron2_dr03HcalTowerSumEt[i];
      electron2[i].dr03TkSumPt	= electron2_dr03TkSumPt[i];
      electron2[i].ecalIso	= electron2_ecalIso[i];
      electron2[i].eidRobustTight	= electron2_eidRobustTight[i];
      electron2[i].energy	= electron2_energy[i];
      electron2[i].et	= electron2_et[i];
      electron2[i].eta	= electron2_eta[i];
      electron2[i].genLepton_eta	= electron2_genLepton_eta[i];
      electron2[i].genLepton_pdgId	= electron2_genLepton_pdgId[i];
      electron2[i].genLepton_phi	= electron2_genLepton_phi[i];
      electron2[i].genLepton_pt	= electron2_genLepton_pt[i];
      electron2[i].gsfTrack_d0	= electron2_gsfTrack_d0[i];
      electron2[i].gsfTrack_trackerExpectedHitsInner_numberOfLostHits	= electron2_gsfTrack_trackerExpectedHitsInner_numberOfLostHits[i];
      electron2[i].hadronicOverEm	= electron2_hadronicOverEm[i];
      electron2[i].hcalIso	= electron2_hcalIso[i];
      electron2[i].neutralHadronIso	= electron2_neutralHadronIso[i];
      electron2[i].phi	= electron2_phi[i];
      electron2[i].photonIso	= electron2_photonIso[i];
      electron2[i].pt	= electron2_pt[i];
      electron2[i].px	= electron2_px[i];
      electron2[i].py	= electron2_py[i];
      electron2[i].pz	= electron2_pz[i];
      electron2[i].sigmaIetaIeta	= electron2_sigmaIetaIeta[i];
      electron2[i].simpleEleId95relIso	= electron2_simpleEleId95relIso[i];
      electron2[i].superCluster_eta	= electron2_superCluster_eta[i];
      electron2[i].trackIso	= electron2_trackIso[i];
      electron2[i].vz	= electron2_vz[i];
    }

  electron3.resize(electron3_caloIso.size());
  for(unsigned int i=0; i < electron3_caloIso.size(); ++i)
    {
      electron3[i].caloIso	= electron3_caloIso[i];
      electron3[i].charge	= electron3_charge[i];
      electron3[i].chargedHadronIso	= electron3_chargedHadronIso[i];
      electron3[i].dB	= electron3_dB[i];
      electron3[i].deltaEtaSuperClusterTrackAtVtx	= electron3_deltaEtaSuperClusterTrackAtVtx[i];
      electron3[i].deltaPhiSuperClusterTrackAtVtx	= electron3_deltaPhiSuperClusterTrackAtVtx[i];
      electron3[i].dr03EcalRecHitSumEt	= electron3_dr03EcalRecHitSumEt[i];
      electron3[i].dr03HcalTowerSumEt	= electron3_dr03HcalTowerSumEt[i];
      electron3[i].dr03TkSumPt	= electron3_dr03TkSumPt[i];
      electron3[i].ecalIso	= electron3_ecalIso[i];
      electron3[i].eidRobustTight	= electron3_eidRobustTight[i];
      electron3[i].energy	= electron3_energy[i];
      electron3[i].et	= electron3_et[i];
      electron3[i].eta	= electron3_eta[i];
      electron3[i].genLepton_eta	= electron3_genLepton_eta[i];
      electron3[i].genLepton_pdgId	= electron3_genLepton_pdgId[i];
      electron3[i].genLepton_phi	= electron3_genLepton_phi[i];
      electron3[i].genLepton_pt	= electron3_genLepton_pt[i];
      electron3[i].gsfTrack_d0	= electron3_gsfTrack_d0[i];
      electron3[i].gsfTrack_trackerExpectedHitsInner_numberOfLostHits	= electron3_gsfTrack_trackerExpectedHitsInner_numberOfLostHits[i];
      electron3[i].hadronicOverEm	= electron3_hadronicOverEm[i];
      electron3[i].hcalIso	= electron3_hcalIso[i];
      electron3[i].neutralHadronIso	= electron3_neutralHadronIso[i];
      electron3[i].phi	= electron3_phi[i];
      electron3[i].photonIso	= electron3_photonIso[i];
      electron3[i].pt	= electron3_pt[i];
      electron3[i].px	= electron3_px[i];
      electron3[i].py	= electron3_py[i];
      electron3[i].pz	= electron3_pz[i];
      electron3[i].sigmaIetaIeta	= electron3_sigmaIetaIeta[i];
      electron3[i].simpleEleId95relIso	= electron3_simpleEleId95relIso[i];
      electron3[i].superCluster_eta	= electron3_superCluster_eta[i];
      electron3[i].trackIso	= electron3_trackIso[i];
      electron3[i].vz	= electron3_vz[i];
    }

  geneventinfoproduct.resize(geneventinfoproduct_pdf1.size());
  for(unsigned int i=0; i < geneventinfoproduct_pdf1.size(); ++i)
    {
      geneventinfoproduct[i].pdf1	= geneventinfoproduct_pdf1[i];
      geneventinfoproduct[i].pdf2	= geneventinfoproduct_pdf2[i];
      geneventinfoproduct[i].pdfweight	= geneventinfoproduct_pdfweight[i];
      geneventinfoproduct[i].pdfweightsum	= geneventinfoproduct_pdfweightsum[i];
    }

  genparticlera2.resize(genparticlera2_charge.size());
  for(unsigned int i=0; i < genparticlera2_charge.size(); ++i)
    {
      genparticlera2[i].charge	= genparticlera2_charge[i];
      genparticlera2[i].eta	= genparticlera2_eta[i];
      genparticlera2[i].firstDaughter	= genparticlera2_firstDaughter[i];
      genparticlera2[i].firstMother	= genparticlera2_firstMother[i];
      genparticlera2[i].lastDaughter	= genparticlera2_lastDaughter[i];
      genparticlera2[i].lastMother	= genparticlera2_lastMother[i];
      genparticlera2[i].mass	= genparticlera2_mass[i];
      genparticlera2[i].pdgId	= genparticlera2_pdgId[i];
      genparticlera2[i].phi	= genparticlera2_phi[i];
      genparticlera2[i].pt	= genparticlera2_pt[i];
      genparticlera2[i].status	= genparticlera2_status[i];
    }

  jet.resize(jet_jetUncMinus.size());
  for(unsigned int i=0; i < jet_jetUncMinus.size(); ++i)
    {
      jet[i].jetUncMinus	= jet_jetUncMinus[i];
      jet[i].jetUncPlus	= jet_jetUncPlus[i];
    }

  jet1.resize(jet1_jetUncMinus.size());
  for(unsigned int i=0; i < jet1_jetUncMinus.size(); ++i)
    {
      jet1[i].jetUncMinus	= jet1_jetUncMinus[i];
      jet1[i].jetUncPlus	= jet1_jetUncPlus[i];
    }

  jet2.resize(jet2_combinedSecondaryVertexBJetTags.size());
  for(unsigned int i=0; i < jet2_combinedSecondaryVertexBJetTags.size(); ++i)
    {
      jet2[i].combinedSecondaryVertexBJetTags	= jet2_combinedSecondaryVertexBJetTags[i];
      jet2[i].combinedSecondaryVertexMVABJetTags	= jet2_combinedSecondaryVertexMVABJetTags[i];
      jet2[i].emEnergyFraction	= jet2_emEnergyFraction[i];
      jet2[i].energy	= jet2_energy[i];
      jet2[i].eta	= jet2_eta[i];
      jet2[i].jetBProbabilityBJetTags	= jet2_jetBProbabilityBJetTags[i];
      jet2[i].jetID_fHPD	= jet2_jetID_fHPD[i];
      jet2[i].jetID_n90Hits	= jet2_jetID_n90Hits[i];
      jet2[i].jetProbabilityBJetTags	= jet2_jetProbabilityBJetTags[i];
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

  jet3.resize(jet3_HFHadronEnergy.size());
  for(unsigned int i=0; i < jet3_HFHadronEnergy.size(); ++i)
    {
      jet3[i].HFHadronEnergy	= jet3_HFHadronEnergy[i];
      jet3[i].chargedEmEnergyFraction	= jet3_chargedEmEnergyFraction[i];
      jet3[i].chargedHadronEnergyFraction	= jet3_chargedHadronEnergyFraction[i];
      jet3[i].chargedMultiplicity	= jet3_chargedMultiplicity[i];
      jet3[i].combinedSecondaryVertexBJetTags	= jet3_combinedSecondaryVertexBJetTags[i];
      jet3[i].combinedSecondaryVertexMVABJetTags	= jet3_combinedSecondaryVertexMVABJetTags[i];
      jet3[i].energy	= jet3_energy[i];
      jet3[i].eta	= jet3_eta[i];
      jet3[i].jetBProbabilityBJetTags	= jet3_jetBProbabilityBJetTags[i];
      jet3[i].jetID_fHPD	= jet3_jetID_fHPD[i];
      jet3[i].jetID_n90Hits	= jet3_jetID_n90Hits[i];
      jet3[i].jetProbabilityBJetTags	= jet3_jetProbabilityBJetTags[i];
      jet3[i].neutralEmEnergyFraction	= jet3_neutralEmEnergyFraction[i];
      jet3[i].neutralHadronEnergy	= jet3_neutralHadronEnergy[i];
      jet3[i].neutralHadronEnergyFraction	= jet3_neutralHadronEnergyFraction[i];
      jet3[i].numberOfDaughters	= jet3_numberOfDaughters[i];
      jet3[i].partonFlavour	= jet3_partonFlavour[i];
      jet3[i].phi	= jet3_phi[i];
      jet3[i].pt	= jet3_pt[i];
      jet3[i].simpleSecondaryVertexBJetTags	= jet3_simpleSecondaryVertexBJetTags[i];
      jet3[i].simpleSecondaryVertexHighEffBJetTags	= jet3_simpleSecondaryVertexHighEffBJetTags[i];
      jet3[i].simpleSecondaryVertexHighPurBJetTags	= jet3_simpleSecondaryVertexHighPurBJetTags[i];
      jet3[i].softElectronByIP3dBJetTags	= jet3_softElectronByIP3dBJetTags[i];
      jet3[i].softElectronByPtBJetTags	= jet3_softElectronByPtBJetTags[i];
      jet3[i].softMuonBJetTags	= jet3_softMuonBJetTags[i];
      jet3[i].softMuonByIP3dBJetTags	= jet3_softMuonByIP3dBJetTags[i];
      jet3[i].softMuonByPtBJetTags	= jet3_softMuonByPtBJetTags[i];
      jet3[i].trackCountingHighEffBJetTags	= jet3_trackCountingHighEffBJetTags[i];
      jet3[i].trackCountingHighPurBJetTags	= jet3_trackCountingHighPurBJetTags[i];
      jet3[i].uncor_energy	= jet3_uncor_energy[i];
      jet3[i].uncor_eta	= jet3_uncor_eta[i];
      jet3[i].uncor_phi	= jet3_uncor_phi[i];
      jet3[i].uncor_pt	= jet3_uncor_pt[i];
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

  muon.resize(muon_dxywrtBeamSpot.size());
  for(unsigned int i=0; i < muon_dxywrtBeamSpot.size(); ++i)
    {
      muon[i].dxywrtBeamSpot	= muon_dxywrtBeamSpot[i];
    }

  muon1.resize(muon1_dxywrtBeamSpot.size());
  for(unsigned int i=0; i < muon1_dxywrtBeamSpot.size(); ++i)
    {
      muon1[i].dxywrtBeamSpot	= muon1_dxywrtBeamSpot[i];
    }

  muon2.resize(muon2_AllGlobalMuons.size());
  for(unsigned int i=0; i < muon2_AllGlobalMuons.size(); ++i)
    {
      muon2[i].AllGlobalMuons	= muon2_AllGlobalMuons[i];
      muon2[i].GlobalMuonPromptTight	= muon2_GlobalMuonPromptTight[i];
      muon2[i].charge	= muon2_charge[i];
      muon2[i].chargedHadronIso	= muon2_chargedHadronIso[i];
      muon2[i].combinedMuon_chi2	= muon2_combinedMuon_chi2[i];
      muon2[i].combinedMuon_ndof	= muon2_combinedMuon_ndof[i];
      muon2[i].combinedMuon_numberOfValidHits	= muon2_combinedMuon_numberOfValidHits[i];
      muon2[i].dB	= muon2_dB[i];
      muon2[i].ecalIso	= muon2_ecalIso[i];
      muon2[i].energy	= muon2_energy[i];
      muon2[i].eta	= muon2_eta[i];
      muon2[i].genLepton_eta	= muon2_genLepton_eta[i];
      muon2[i].genLepton_pdgId	= muon2_genLepton_pdgId[i];
      muon2[i].genLepton_phi	= muon2_genLepton_phi[i];
      muon2[i].genLepton_pt	= muon2_genLepton_pt[i];
      muon2[i].globalTrack_pt	= muon2_globalTrack_pt[i];
      muon2[i].hcalIso	= muon2_hcalIso[i];
      muon2[i].innerTrack_numberOfValidHits	= muon2_innerTrack_numberOfValidHits[i];
      muon2[i].innerTrack_pt	= muon2_innerTrack_pt[i];
      muon2[i].isGlobalMuon	= muon2_isGlobalMuon[i];
      muon2[i].isTrackerMuon	= muon2_isTrackerMuon[i];
      muon2[i].neutralHadronIso	= muon2_neutralHadronIso[i];
      muon2[i].phi	= muon2_phi[i];
      muon2[i].photonIso	= muon2_photonIso[i];
      muon2[i].pt	= muon2_pt[i];
      muon2[i].px	= muon2_px[i];
      muon2[i].py	= muon2_py[i];
      muon2[i].pz	= muon2_pz[i];
      muon2[i].trackIso	= muon2_trackIso[i];
      muon2[i].track_d0	= muon2_track_d0[i];
      muon2[i].track_hitPattern_numberOfValidMuonHits	= muon2_track_hitPattern_numberOfValidMuonHits[i];
      muon2[i].track_hitPattern_numberOfValidPixelHits	= muon2_track_hitPattern_numberOfValidPixelHits[i];
      muon2[i].track_normalizedChi2	= muon2_track_normalizedChi2[i];
      muon2[i].vz	= muon2_vz[i];
    }

  muon3.resize(muon3_AllGlobalMuons.size());
  for(unsigned int i=0; i < muon3_AllGlobalMuons.size(); ++i)
    {
      muon3[i].AllGlobalMuons	= muon3_AllGlobalMuons[i];
      muon3[i].GlobalMuonPromptTight	= muon3_GlobalMuonPromptTight[i];
      muon3[i].charge	= muon3_charge[i];
      muon3[i].chargedHadronIso	= muon3_chargedHadronIso[i];
      muon3[i].combinedMuon_chi2	= muon3_combinedMuon_chi2[i];
      muon3[i].combinedMuon_ndof	= muon3_combinedMuon_ndof[i];
      muon3[i].combinedMuon_numberOfValidHits	= muon3_combinedMuon_numberOfValidHits[i];
      muon3[i].dB	= muon3_dB[i];
      muon3[i].ecalIso	= muon3_ecalIso[i];
      muon3[i].energy	= muon3_energy[i];
      muon3[i].eta	= muon3_eta[i];
      muon3[i].genLepton_eta	= muon3_genLepton_eta[i];
      muon3[i].genLepton_pdgId	= muon3_genLepton_pdgId[i];
      muon3[i].genLepton_phi	= muon3_genLepton_phi[i];
      muon3[i].genLepton_pt	= muon3_genLepton_pt[i];
      muon3[i].globalTrack_pt	= muon3_globalTrack_pt[i];
      muon3[i].hcalIso	= muon3_hcalIso[i];
      muon3[i].innerTrack_numberOfValidHits	= muon3_innerTrack_numberOfValidHits[i];
      muon3[i].innerTrack_pt	= muon3_innerTrack_pt[i];
      muon3[i].isGlobalMuon	= muon3_isGlobalMuon[i];
      muon3[i].isTrackerMuon	= muon3_isTrackerMuon[i];
      muon3[i].neutralHadronIso	= muon3_neutralHadronIso[i];
      muon3[i].phi	= muon3_phi[i];
      muon3[i].photonIso	= muon3_photonIso[i];
      muon3[i].pt	= muon3_pt[i];
      muon3[i].px	= muon3_px[i];
      muon3[i].py	= muon3_py[i];
      muon3[i].pz	= muon3_pz[i];
      muon3[i].trackIso	= muon3_trackIso[i];
      muon3[i].track_d0	= muon3_track_d0[i];
      muon3[i].track_hitPattern_numberOfValidMuonHits	= muon3_track_hitPattern_numberOfValidMuonHits[i];
      muon3[i].track_hitPattern_numberOfValidPixelHits	= muon3_track_hitPattern_numberOfValidPixelHits[i];
      muon3[i].track_normalizedChi2	= muon3_track_normalizedChi2[i];
      muon3[i].vz	= muon3_vz[i];
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
