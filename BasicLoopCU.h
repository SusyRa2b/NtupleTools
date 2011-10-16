#ifndef BASICLOOPCU_H
#define BASICLOOPCU_H
//-----------------------------------------------------------------------------
// File:        BasicLoopCU.h
// Description: Analyzer header for ntuples created by TheNtupleMaker
// Created:     Wed Oct 12 11:42:35 2011 by mkntanalyzer.py
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
double	beamspot_x0;
double	beamspot_y0;
double	beamspot_z0;
std::vector<double>	electron1_caloIso(7,0);
std::vector<double>	electron1_charge(7,0);
std::vector<double>	electron1_chargedHadronIso(7,0);
std::vector<double>	electron1_convDcot(7,0);
std::vector<double>	electron1_convDist(7,0);
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
std::vector<double>	electron1_sigmaIetaIeta(7,0);
std::vector<double>	electron1_superCluster_eta(7,0);
std::vector<double>	electron1_trackIso(7,0);
std::vector<double>	electron1_vertex_z(7,0);
std::vector<double>	electron1_vz(7,0);
std::vector<double>	electron_caloIso(15,0);
std::vector<double>	electron_charge(15,0);
std::vector<double>	electron_chargedHadronIso(15,0);
std::vector<double>	electron_convDcot(15,0);
std::vector<double>	electron_convDist(15,0);
std::vector<double>	electron_dB(15,0);
std::vector<double>	electron_deltaEtaSuperClusterTrackAtVtx(15,0);
std::vector<double>	electron_deltaPhiSuperClusterTrackAtVtx(15,0);
std::vector<double>	electron_dr03EcalRecHitSumEt(15,0);
std::vector<double>	electron_dr03HcalTowerSumEt(15,0);
std::vector<double>	electron_dr03TkSumPt(15,0);
std::vector<double>	electron_ecalIso(15,0);
std::vector<double>	electron_eidRobustTight(15,0);
std::vector<double>	electron_energy(15,0);
std::vector<double>	electron_et(15,0);
std::vector<double>	electron_eta(15,0);
std::vector<double>	electron_genLepton_eta(15,0);
std::vector<double>	electron_genLepton_pdgId(15,0);
std::vector<double>	electron_genLepton_phi(15,0);
std::vector<double>	electron_genLepton_pt(15,0);
std::vector<double>	electron_gsfTrack_d0(15,0);
std::vector<double>	electron_gsfTrack_trackerExpectedHitsInner_numberOfLostHits(15,0);
std::vector<double>	electron_hadronicOverEm(15,0);
std::vector<double>	electron_hcalIso(15,0);
std::vector<double>	electron_neutralHadronIso(15,0);
std::vector<double>	electron_phi(15,0);
std::vector<double>	electron_photonIso(15,0);
std::vector<double>	electron_pt(15,0);
std::vector<double>	electron_sigmaIetaIeta(15,0);
std::vector<double>	electron_simpleEleId80cIso(15,0);
std::vector<double>	electron_simpleEleId80relIso(15,0);
std::vector<double>	electron_simpleEleId95cIso(15,0);
std::vector<double>	electron_simpleEleId95relIso(15,0);
std::vector<double>	electron_superCluster_eta(15,0);
std::vector<double>	electron_trackIso(15,0);
std::vector<double>	electron_vertex_z(15,0);
std::vector<double>	electron_vz(15,0);
std::vector<double>	electronhelper1_dxywrtBeamSpot(7,0);
std::vector<double>	electronhelper_dxywrtBeamSpot(15,0);
double	eventhelper_bunchCrossing;
double	eventhelper_event;
double	eventhelper_isRealData;
double	eventhelper_luminosityBlock;
double	eventhelper_run;
double	eventhelperextra_trackMPT1phi;
double	eventhelperextra_trackMPT1pt;
double	eventhelperextra_trackMPT5phi;
double	eventhelperextra_trackMPT5pt;
double	eventlhehelperextra_m0;
double	eventlhehelperextra_m12;
double	eventlhehelperextra_mGL;
double	eventlhehelperextra_mLSP;
double	eventlhehelperextra_xCHI;
double	geneventinfoproduct_pdf1;
double	geneventinfoproduct_pdf2;
double	geneventinfoproduct_scalePDF;
double	geneventinfoproduct_weight;
double	geneventinfoproduct_x1;
double	geneventinfoproduct_x2;
std::vector<double>	geneventinfoproducthelper1_pdf1(91,0);
std::vector<double>	geneventinfoproducthelper1_pdf2(91,0);
std::vector<double>	geneventinfoproducthelper1_pdfweight(91,0);
std::vector<double>	geneventinfoproducthelper1_pdfweightsum(91,0);
std::vector<double>	geneventinfoproducthelper2_pdf1(83,0);
std::vector<double>	geneventinfoproducthelper2_pdf2(83,0);
std::vector<double>	geneventinfoproducthelper2_pdfweight(83,0);
std::vector<double>	geneventinfoproducthelper2_pdfweightsum(83,0);
std::vector<double>	geneventinfoproducthelper_pdf1(201,0);
std::vector<double>	geneventinfoproducthelper_pdf2(201,0);
std::vector<double>	geneventinfoproducthelper_pdfweight(201,0);
std::vector<double>	geneventinfoproducthelper_pdfweightsum(201,0);
std::vector<double>	genparticlehelperra2_charge(73,0);
std::vector<double>	genparticlehelperra2_eta(73,0);
std::vector<double>	genparticlehelperra2_firstDaughter(73,0);
std::vector<double>	genparticlehelperra2_firstMother(73,0);
std::vector<double>	genparticlehelperra2_lastDaughter(73,0);
std::vector<double>	genparticlehelperra2_lastMother(73,0);
std::vector<double>	genparticlehelperra2_mass(73,0);
std::vector<double>	genparticlehelperra2_pdgId(73,0);
std::vector<double>	genparticlehelperra2_phi(73,0);
std::vector<double>	genparticlehelperra2_pt(73,0);
std::vector<double>	genparticlehelperra2_status(73,0);
double	genruninfoproduct_externalXSecLO_error;
double	genruninfoproduct_externalXSecLO_value;
double	genruninfoproduct_filterEfficiency;
double	genruninfoproduct_internalXSec_error;
double	genruninfoproduct_internalXSec_value;
std::vector<double>	jet1_HFEMEnergy(177,0);
std::vector<double>	jet1_HFEMMultiplicity(177,0);
std::vector<double>	jet1_HFHadronEnergy(177,0);
std::vector<double>	jet1_HFHadronMultiplicity(177,0);
std::vector<double>	jet1_chargedEmEnergyFraction(177,0);
std::vector<double>	jet1_chargedHadronEnergyFraction(177,0);
std::vector<double>	jet1_chargedHadronMultiplicity(177,0);
std::vector<double>	jet1_chargedMuEnergy(177,0);
std::vector<double>	jet1_chargedMultiplicity(177,0);
std::vector<double>	jet1_combinedSecondaryVertexBJetTags(177,0);
std::vector<double>	jet1_combinedSecondaryVertexMVABJetTags(177,0);
std::vector<double>	jet1_electronEnergy(177,0);
std::vector<double>	jet1_electronMultiplicity(177,0);
std::vector<double>	jet1_energy(177,0);
std::vector<double>	jet1_et(177,0);
std::vector<double>	jet1_eta(177,0);
std::vector<double>	jet1_genJet_energy(177,0);
std::vector<double>	jet1_genJet_eta(177,0);
std::vector<double>	jet1_genJet_invisibleEnergy(177,0);
std::vector<double>	jet1_genJet_phi(177,0);
std::vector<double>	jet1_genJet_pt(177,0);
std::vector<double>	jet1_genParton_energy(177,0);
std::vector<double>	jet1_genParton_eta(177,0);
std::vector<double>	jet1_genParton_pdgId(177,0);
std::vector<double>	jet1_genParton_phi(177,0);
std::vector<double>	jet1_genParton_pt(177,0);
std::vector<double>	jet1_jetArea(177,0);
std::vector<double>	jet1_jetBProbabilityBJetTags(177,0);
std::vector<double>	jet1_jetID_fHPD(177,0);
std::vector<double>	jet1_jetID_n90Hits(177,0);
std::vector<double>	jet1_jetProbabilityBJetTags(177,0);
std::vector<double>	jet1_muonEnergy(177,0);
std::vector<double>	jet1_neutralEmEnergyFraction(177,0);
std::vector<double>	jet1_neutralHadronEnergy(177,0);
std::vector<double>	jet1_neutralHadronEnergyFraction(177,0);
std::vector<double>	jet1_neutralHadronMultiplicity(177,0);
std::vector<double>	jet1_neutralMultiplicity(177,0);
std::vector<double>	jet1_numberOfDaughters(177,0);
std::vector<double>	jet1_partonFlavour(177,0);
std::vector<double>	jet1_phi(177,0);
std::vector<double>	jet1_photonEnergy(177,0);
std::vector<double>	jet1_photonMultiplicity(177,0);
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
std::vector<double>	jet1_uncor_et(177,0);
std::vector<double>	jet1_uncor_eta(177,0);
std::vector<double>	jet1_uncor_phi(177,0);
std::vector<double>	jet1_uncor_pt(177,0);
std::vector<double>	jet2_HFEMEnergy(177,0);
std::vector<double>	jet2_HFEMMultiplicity(177,0);
std::vector<double>	jet2_HFHadronEnergy(177,0);
std::vector<double>	jet2_HFHadronMultiplicity(177,0);
std::vector<double>	jet2_chargedEmEnergyFraction(177,0);
std::vector<double>	jet2_chargedHadronEnergyFraction(177,0);
std::vector<double>	jet2_chargedHadronMultiplicity(177,0);
std::vector<double>	jet2_chargedMuEnergy(177,0);
std::vector<double>	jet2_chargedMultiplicity(177,0);
std::vector<double>	jet2_combinedSecondaryVertexBJetTags(177,0);
std::vector<double>	jet2_combinedSecondaryVertexMVABJetTags(177,0);
std::vector<double>	jet2_electronEnergy(177,0);
std::vector<double>	jet2_electronMultiplicity(177,0);
std::vector<double>	jet2_energy(177,0);
std::vector<double>	jet2_et(177,0);
std::vector<double>	jet2_eta(177,0);
std::vector<double>	jet2_genJet_energy(177,0);
std::vector<double>	jet2_genJet_eta(177,0);
std::vector<double>	jet2_genJet_invisibleEnergy(177,0);
std::vector<double>	jet2_genJet_phi(177,0);
std::vector<double>	jet2_genJet_pt(177,0);
std::vector<double>	jet2_genParton_energy(177,0);
std::vector<double>	jet2_genParton_eta(177,0);
std::vector<double>	jet2_genParton_pdgId(177,0);
std::vector<double>	jet2_genParton_phi(177,0);
std::vector<double>	jet2_genParton_pt(177,0);
std::vector<double>	jet2_jetArea(177,0);
std::vector<double>	jet2_jetBProbabilityBJetTags(177,0);
std::vector<double>	jet2_jetID_fHPD(177,0);
std::vector<double>	jet2_jetID_n90Hits(177,0);
std::vector<double>	jet2_jetProbabilityBJetTags(177,0);
std::vector<double>	jet2_muonEnergy(177,0);
std::vector<double>	jet2_neutralEmEnergyFraction(177,0);
std::vector<double>	jet2_neutralHadronEnergy(177,0);
std::vector<double>	jet2_neutralHadronEnergyFraction(177,0);
std::vector<double>	jet2_neutralHadronMultiplicity(177,0);
std::vector<double>	jet2_neutralMultiplicity(177,0);
std::vector<double>	jet2_numberOfDaughters(177,0);
std::vector<double>	jet2_partonFlavour(177,0);
std::vector<double>	jet2_phi(177,0);
std::vector<double>	jet2_photonEnergy(177,0);
std::vector<double>	jet2_photonMultiplicity(177,0);
std::vector<double>	jet2_pt(177,0);
std::vector<double>	jet2_simpleSecondaryVertexBJetTags(177,0);
std::vector<double>	jet2_simpleSecondaryVertexHighEffBJetTags(177,0);
std::vector<double>	jet2_simpleSecondaryVertexHighPurBJetTags(177,0);
std::vector<double>	jet2_softElectronByIP3dBJetTags(177,0);
std::vector<double>	jet2_softElectronByPtBJetTags(177,0);
std::vector<double>	jet2_softMuonBJetTags(177,0);
std::vector<double>	jet2_softMuonByIP3dBJetTags(177,0);
std::vector<double>	jet2_softMuonByPtBJetTags(177,0);
std::vector<double>	jet2_trackCountingHighEffBJetTags(177,0);
std::vector<double>	jet2_trackCountingHighPurBJetTags(177,0);
std::vector<double>	jet2_uncor_energy(177,0);
std::vector<double>	jet2_uncor_et(177,0);
std::vector<double>	jet2_uncor_eta(177,0);
std::vector<double>	jet2_uncor_phi(177,0);
std::vector<double>	jet2_uncor_pt(177,0);
std::vector<double>	jet_HFEMEnergy(205,0);
std::vector<double>	jet_HFEMMultiplicity(205,0);
std::vector<double>	jet_HFHadronEnergy(205,0);
std::vector<double>	jet_HFHadronMultiplicity(205,0);
std::vector<double>	jet_chargedEmEnergyFraction(205,0);
std::vector<double>	jet_chargedHadronEnergyFraction(205,0);
std::vector<double>	jet_chargedHadronMultiplicity(205,0);
std::vector<double>	jet_chargedMuEnergy(205,0);
std::vector<double>	jet_chargedMultiplicity(205,0);
std::vector<double>	jet_combinedSecondaryVertexBJetTags(205,0);
std::vector<double>	jet_combinedSecondaryVertexMVABJetTags(205,0);
std::vector<double>	jet_electronEnergy(205,0);
std::vector<double>	jet_electronMultiplicity(205,0);
std::vector<double>	jet_energy(205,0);
std::vector<double>	jet_et(205,0);
std::vector<double>	jet_eta(205,0);
std::vector<double>	jet_genJet_energy(205,0);
std::vector<double>	jet_genJet_eta(205,0);
std::vector<double>	jet_genJet_invisibleEnergy(205,0);
std::vector<double>	jet_genJet_phi(205,0);
std::vector<double>	jet_genJet_pt(205,0);
std::vector<double>	jet_genParton_energy(205,0);
std::vector<double>	jet_genParton_eta(205,0);
std::vector<double>	jet_genParton_pdgId(205,0);
std::vector<double>	jet_genParton_phi(205,0);
std::vector<double>	jet_genParton_pt(205,0);
std::vector<double>	jet_jetArea(205,0);
std::vector<double>	jet_jetBProbabilityBJetTags(205,0);
std::vector<double>	jet_jetID_fHPD(205,0);
std::vector<double>	jet_jetID_n90Hits(205,0);
std::vector<double>	jet_jetProbabilityBJetTags(205,0);
std::vector<double>	jet_muonEnergy(205,0);
std::vector<double>	jet_neutralEmEnergyFraction(205,0);
std::vector<double>	jet_neutralHadronEnergy(205,0);
std::vector<double>	jet_neutralHadronEnergyFraction(205,0);
std::vector<double>	jet_neutralHadronMultiplicity(205,0);
std::vector<double>	jet_neutralMultiplicity(205,0);
std::vector<double>	jet_numberOfDaughters(205,0);
std::vector<double>	jet_partonFlavour(205,0);
std::vector<double>	jet_phi(205,0);
std::vector<double>	jet_photonEnergy(205,0);
std::vector<double>	jet_photonMultiplicity(205,0);
std::vector<double>	jet_pt(205,0);
std::vector<double>	jet_simpleSecondaryVertexBJetTags(205,0);
std::vector<double>	jet_simpleSecondaryVertexHighEffBJetTags(205,0);
std::vector<double>	jet_simpleSecondaryVertexHighPurBJetTags(205,0);
std::vector<double>	jet_softElectronByIP3dBJetTags(205,0);
std::vector<double>	jet_softElectronByPtBJetTags(205,0);
std::vector<double>	jet_softMuonBJetTags(205,0);
std::vector<double>	jet_softMuonByIP3dBJetTags(205,0);
std::vector<double>	jet_softMuonByPtBJetTags(205,0);
std::vector<double>	jet_trackCountingHighEffBJetTags(205,0);
std::vector<double>	jet_trackCountingHighPurBJetTags(205,0);
std::vector<double>	jet_uncor_energy(205,0);
std::vector<double>	jet_uncor_et(205,0);
std::vector<double>	jet_uncor_eta(205,0);
std::vector<double>	jet_uncor_phi(205,0);
std::vector<double>	jet_uncor_pt(205,0);
std::vector<double>	jethelper1_jecFactor(177,0);
std::vector<double>	jethelper1_jecFactorNoL1Fast(177,0);
std::vector<double>	jethelper1_jetUncMinus(177,0);
std::vector<double>	jethelper1_jetUncPlus(177,0);
std::vector<double>	jethelper1_nSecondaryVertices(177,0);
std::vector<double>	jethelper1_secondaryVertex0DirectionX(177,0);
std::vector<double>	jethelper1_secondaryVertex0DirectionY(177,0);
std::vector<double>	jethelper1_secondaryVertex0DirectionZ(177,0);
std::vector<double>	jethelper1_secondaryVertex0FlightDistance(177,0);
std::vector<double>	jethelper1_secondaryVertex0FlightDistanceSig(177,0);
std::vector<double>	jethelper1_secondaryVertex0Mass(177,0);
std::vector<double>	jethelper1_secondaryVertex1DirectionX(177,0);
std::vector<double>	jethelper1_secondaryVertex1DirectionY(177,0);
std::vector<double>	jethelper1_secondaryVertex1DirectionZ(177,0);
std::vector<double>	jethelper1_secondaryVertex1FlightDistance(177,0);
std::vector<double>	jethelper1_secondaryVertex1FlightDistanceSig(177,0);
std::vector<double>	jethelper1_secondaryVertex1Mass(177,0);
std::vector<double>	jethelper2_jecFactor(177,0);
std::vector<double>	jethelper2_jecFactorNoL1Fast(177,0);
std::vector<double>	jethelper2_jetUncMinus(177,0);
std::vector<double>	jethelper2_jetUncPlus(177,0);
std::vector<double>	jethelper2_nSecondaryVertices(177,0);
std::vector<double>	jethelper2_secondaryVertex0DirectionX(177,0);
std::vector<double>	jethelper2_secondaryVertex0DirectionY(177,0);
std::vector<double>	jethelper2_secondaryVertex0DirectionZ(177,0);
std::vector<double>	jethelper2_secondaryVertex0FlightDistance(177,0);
std::vector<double>	jethelper2_secondaryVertex0FlightDistanceSig(177,0);
std::vector<double>	jethelper2_secondaryVertex0Mass(177,0);
std::vector<double>	jethelper2_secondaryVertex1DirectionX(177,0);
std::vector<double>	jethelper2_secondaryVertex1DirectionY(177,0);
std::vector<double>	jethelper2_secondaryVertex1DirectionZ(177,0);
std::vector<double>	jethelper2_secondaryVertex1FlightDistance(177,0);
std::vector<double>	jethelper2_secondaryVertex1FlightDistanceSig(177,0);
std::vector<double>	jethelper2_secondaryVertex1Mass(177,0);
std::vector<double>	jethelper_jecFactor(205,0);
std::vector<double>	jethelper_jecFactorNoL1Fast(205,0);
std::vector<double>	jethelper_jetUncMinus(205,0);
std::vector<double>	jethelper_jetUncPlus(205,0);
std::vector<double>	jethelper_nSecondaryVertices(205,0);
std::vector<double>	jethelper_secondaryVertex0DirectionX(205,0);
std::vector<double>	jethelper_secondaryVertex0DirectionY(205,0);
std::vector<double>	jethelper_secondaryVertex0DirectionZ(205,0);
std::vector<double>	jethelper_secondaryVertex0FlightDistance(205,0);
std::vector<double>	jethelper_secondaryVertex0FlightDistanceSig(205,0);
std::vector<double>	jethelper_secondaryVertex0Mass(205,0);
std::vector<double>	jethelper_secondaryVertex1DirectionX(205,0);
std::vector<double>	jethelper_secondaryVertex1DirectionY(205,0);
std::vector<double>	jethelper_secondaryVertex1DirectionZ(205,0);
std::vector<double>	jethelper_secondaryVertex1FlightDistance(205,0);
std::vector<double>	jethelper_secondaryVertex1FlightDistanceSig(205,0);
std::vector<double>	jethelper_secondaryVertex1Mass(205,0);
std::vector<double>	met1_energy(3,0);
std::vector<double>	met1_et(3,0);
std::vector<double>	met1_genMET_energy(3,0);
std::vector<double>	met1_genMET_et(3,0);
std::vector<double>	met1_genMET_phi(3,0);
std::vector<double>	met1_genMET_pt(3,0);
std::vector<double>	met1_mEtSig(3,0);
std::vector<double>	met1_phi(3,0);
std::vector<double>	met1_pt(3,0);
std::vector<double>	met1_significance(3,0);
std::vector<double>	met1_sumEt(3,0);
std::vector<double>	met2_energy(3,0);
std::vector<double>	met2_et(3,0);
std::vector<double>	met2_genMET_energy(3,0);
std::vector<double>	met2_genMET_et(3,0);
std::vector<double>	met2_genMET_phi(3,0);
std::vector<double>	met2_genMET_pt(3,0);
std::vector<double>	met2_mEtSig(3,0);
std::vector<double>	met2_phi(3,0);
std::vector<double>	met2_pt(3,0);
std::vector<double>	met2_significance(3,0);
std::vector<double>	met2_sumEt(3,0);
std::vector<double>	met3_energy(3,0);
std::vector<double>	met3_et(3,0);
std::vector<double>	met3_genMET_energy(3,0);
std::vector<double>	met3_genMET_et(3,0);
std::vector<double>	met3_genMET_phi(3,0);
std::vector<double>	met3_genMET_pt(3,0);
std::vector<double>	met3_mEtSig(3,0);
std::vector<double>	met3_phi(3,0);
std::vector<double>	met3_pt(3,0);
std::vector<double>	met3_significance(3,0);
std::vector<double>	met3_sumEt(3,0);
std::vector<double>	met4_energy(3,0);
std::vector<double>	met4_et(3,0);
std::vector<double>	met4_genMET_energy(3,0);
std::vector<double>	met4_genMET_et(3,0);
std::vector<double>	met4_genMET_phi(3,0);
std::vector<double>	met4_genMET_pt(3,0);
std::vector<double>	met4_mEtSig(3,0);
std::vector<double>	met4_phi(3,0);
std::vector<double>	met4_pt(3,0);
std::vector<double>	met4_significance(3,0);
std::vector<double>	met4_sumEt(3,0);
std::vector<double>	met5_energy(3,0);
std::vector<double>	met5_et(3,0);
std::vector<double>	met5_genMET_energy(3,0);
std::vector<double>	met5_genMET_et(3,0);
std::vector<double>	met5_genMET_phi(3,0);
std::vector<double>	met5_genMET_pt(3,0);
std::vector<double>	met5_mEtSig(3,0);
std::vector<double>	met5_phi(3,0);
std::vector<double>	met5_pt(3,0);
std::vector<double>	met5_significance(3,0);
std::vector<double>	met5_sumEt(3,0);
std::vector<double>	met6_energy(3,0);
std::vector<double>	met6_et(3,0);
std::vector<double>	met6_genMET_energy(3,0);
std::vector<double>	met6_genMET_et(3,0);
std::vector<double>	met6_genMET_phi(3,0);
std::vector<double>	met6_genMET_pt(3,0);
std::vector<double>	met6_mEtSig(3,0);
std::vector<double>	met6_phi(3,0);
std::vector<double>	met6_pt(3,0);
std::vector<double>	met6_significance(3,0);
std::vector<double>	met6_sumEt(3,0);
std::vector<double>	met7_energy(3,0);
std::vector<double>	met7_et(3,0);
std::vector<double>	met7_genMET_energy(3,0);
std::vector<double>	met7_genMET_et(3,0);
std::vector<double>	met7_genMET_phi(3,0);
std::vector<double>	met7_genMET_pt(3,0);
std::vector<double>	met7_mEtSig(3,0);
std::vector<double>	met7_phi(3,0);
std::vector<double>	met7_pt(3,0);
std::vector<double>	met7_significance(3,0);
std::vector<double>	met7_sumEt(3,0);
std::vector<double>	met_energy(3,0);
std::vector<double>	met_et(3,0);
std::vector<double>	met_genMET_energy(3,0);
std::vector<double>	met_genMET_et(3,0);
std::vector<double>	met_genMET_phi(3,0);
std::vector<double>	met_genMET_pt(3,0);
std::vector<double>	met_mEtSig(3,0);
std::vector<double>	met_metSignificance(3,0);
std::vector<double>	met_phi(3,0);
std::vector<double>	met_pt(3,0);
std::vector<double>	met_significance(3,0);
std::vector<double>	met_sumEt(3,0);
std::vector<double>	methelper1_significance_dxx(3,0);
std::vector<double>	methelper1_significance_dxy(3,0);
std::vector<double>	methelper1_significance_dyx(3,0);
std::vector<double>	methelper1_significance_dyy(3,0);
std::vector<double>	methelper_significance_dxx(3,0);
std::vector<double>	methelper_significance_dxy(3,0);
std::vector<double>	methelper_significance_dyx(3,0);
std::vector<double>	methelper_significance_dyy(3,0);
std::vector<double>	muon1_AllGlobalMuons(7,0);
std::vector<double>	muon1_AllTrackerMuons(7,0);
std::vector<double>	muon1_GlobalMuonPromptTight(7,0);
std::vector<double>	muon1_charge(7,0);
std::vector<double>	muon1_chargedHadronIso(7,0);
std::vector<double>	muon1_combinedMuon_chi2(7,0);
std::vector<double>	muon1_combinedMuon_ndof(7,0);
std::vector<double>	muon1_combinedMuon_numberOfValidHits(7,0);
std::vector<double>	muon1_dB(7,0);
std::vector<double>	muon1_ecalIso(7,0);
std::vector<double>	muon1_energy(7,0);
std::vector<double>	muon1_et(7,0);
std::vector<double>	muon1_eta(7,0);
std::vector<double>	muon1_genLepton_eta(7,0);
std::vector<double>	muon1_genLepton_pdgId(7,0);
std::vector<double>	muon1_genLepton_phi(7,0);
std::vector<double>	muon1_genLepton_pt(7,0);
std::vector<double>	muon1_globalTrack_chi2(7,0);
std::vector<double>	muon1_globalTrack_hitPattern_numberOfValidTrackerHits(7,0);
std::vector<double>	muon1_globalTrack_ndof(7,0);
std::vector<double>	muon1_globalTrack_pt(7,0);
std::vector<double>	muon1_globalTrack_ptError(7,0);
std::vector<double>	muon1_hcalIso(7,0);
std::vector<double>	muon1_innerTrack_d0(7,0);
std::vector<double>	muon1_innerTrack_hitPattern_pixelLayersWithMeasurement(7,0);
std::vector<double>	muon1_innerTrack_numberOfValidHits(7,0);
std::vector<double>	muon1_innerTrack_phi(7,0);
std::vector<double>	muon1_innerTrack_pt(7,0);
std::vector<double>	muon1_innerTrack_vertex_z(7,0);
std::vector<double>	muon1_isGlobalMuon(7,0);
std::vector<double>	muon1_isTrackerMuon(7,0);
std::vector<double>	muon1_neutralHadronIso(7,0);
std::vector<double>	muon1_numberOfMatches(7,0);
std::vector<double>	muon1_phi(7,0);
std::vector<double>	muon1_photonIso(7,0);
std::vector<double>	muon1_pt(7,0);
std::vector<double>	muon1_trackIso(7,0);
std::vector<double>	muon1_track_d0(7,0);
std::vector<double>	muon1_track_hitPattern_numberOfValidMuonHits(7,0);
std::vector<double>	muon1_track_hitPattern_numberOfValidPixelHits(7,0);
std::vector<double>	muon1_track_normalizedChi2(7,0);
std::vector<double>	muon1_vz(7,0);
std::vector<double>	muon_AllGlobalMuons(29,0);
std::vector<double>	muon_AllTrackerMuons(29,0);
std::vector<double>	muon_GlobalMuonPromptTight(29,0);
std::vector<double>	muon_charge(29,0);
std::vector<double>	muon_chargedHadronIso(29,0);
std::vector<double>	muon_combinedMuon_chi2(29,0);
std::vector<double>	muon_combinedMuon_ndof(29,0);
std::vector<double>	muon_combinedMuon_numberOfValidHits(29,0);
std::vector<double>	muon_dB(29,0);
std::vector<double>	muon_ecalIso(29,0);
std::vector<double>	muon_energy(29,0);
std::vector<double>	muon_et(29,0);
std::vector<double>	muon_eta(29,0);
std::vector<double>	muon_genLepton_eta(29,0);
std::vector<double>	muon_genLepton_pdgId(29,0);
std::vector<double>	muon_genLepton_phi(29,0);
std::vector<double>	muon_genLepton_pt(29,0);
std::vector<double>	muon_globalTrack_chi2(29,0);
std::vector<double>	muon_globalTrack_hitPattern_numberOfValidTrackerHits(29,0);
std::vector<double>	muon_globalTrack_ndof(29,0);
std::vector<double>	muon_globalTrack_pt(29,0);
std::vector<double>	muon_globalTrack_ptError(29,0);
std::vector<double>	muon_hcalIso(29,0);
std::vector<double>	muon_innerTrack_d0(29,0);
std::vector<double>	muon_innerTrack_hitPattern_pixelLayersWithMeasurement(29,0);
std::vector<double>	muon_innerTrack_numberOfValidHits(29,0);
std::vector<double>	muon_innerTrack_phi(29,0);
std::vector<double>	muon_innerTrack_pt(29,0);
std::vector<double>	muon_innerTrack_vertex_z(29,0);
std::vector<double>	muon_isGlobalMuon(29,0);
std::vector<double>	muon_isTrackerMuon(29,0);
std::vector<double>	muon_neutralHadronIso(29,0);
std::vector<double>	muon_numberOfMatches(29,0);
std::vector<double>	muon_phi(29,0);
std::vector<double>	muon_photonIso(29,0);
std::vector<double>	muon_pt(29,0);
std::vector<double>	muon_trackIso(29,0);
std::vector<double>	muon_track_d0(29,0);
std::vector<double>	muon_track_hitPattern_numberOfValidMuonHits(29,0);
std::vector<double>	muon_track_hitPattern_numberOfValidPixelHits(29,0);
std::vector<double>	muon_track_normalizedChi2(29,0);
std::vector<double>	muon_vz(29,0);
std::vector<double>	muonhelper1_dxywrtBeamSpot(7,0);
std::vector<double>	muonhelper1_muonPFcandptdiff(7,0);
std::vector<double>	muonhelper_dxywrtBeamSpot(29,0);
std::vector<double>	muonhelper_muonPFcandptdiff(29,0);
int	nGenEventInfoProductHelper_generator_NNPDF;
int	nGenEventInfoProductHelper_generator_cteq;
int	nGenEventInfoProductHelper_generator_mstw;
int	nPileupSummaryInfo_addPileupInfo;
int	nelectron;
int	nelectron1;
int	nelectronhelper;
int	nelectronhelper1;
int	ngenparticlehelperra2b;
int	njet;
int	njet1;
int	njet2;
int	njethelper;
int	njethelper1;
int	njethelper2;
int	nmet;
int	nmet1;
int	nmet2;
int	nmet3;
int	nmet4;
int	nmet5;
int	nmet6;
int	nmet7;
int	nmethelper;
int	nmethelper1;
int	nmuon;
int	nmuon1;
int	nmuonhelper;
int	nmuonhelper1;
int	ntau;
int	ntauhelper;
int	nvertex;
std::vector<double>	pileupsummaryinfo_addpileupinfo_getBunchCrossing(7,0);
std::vector<double>	pileupsummaryinfo_addpileupinfo_getPU_NumInteractions(7,0);
double	sdouble_kt6pfjets_rho_value;
double	sdouble_kt6pfjets_sigma_value;
double	suint_flavorhistoryfilter_value;
std::vector<double>	tau_caloIso(7,0);
std::vector<double>	tau_ecalIso(7,0);
std::vector<double>	tau_energy(7,0);
std::vector<double>	tau_et(7,0);
std::vector<double>	tau_eta(7,0);
std::vector<double>	tau_hcalIso(7,0);
std::vector<double>	tau_phi(7,0);
std::vector<double>	tau_pt(7,0);
std::vector<double>	tau_tauID_againstElectron(7,0);
std::vector<double>	tau_tauID_againstMuon(7,0);
std::vector<double>	tau_tauID_byIsolation(7,0);
std::vector<double>	tau_tauID_byTaNC(7,0);
std::vector<double>	tau_tauID_byTaNCfrHalfPercent(7,0);
std::vector<double>	tau_tauID_byTaNCfrQuarterPercent(7,0);
std::vector<double>	tau_trackIso(7,0);
std::vector<double>	tauhelper_genTauDecayModeID(7,0);
double	triggerresultshelper1_csctighthaloFilter;
double	triggerresultshelper1_eenoiseFilter;
double	triggerresultshelper1_greedymuonFilter;
double	triggerresultshelper1_hbhenoiseFilter;
double	triggerresultshelper1_inconsistentmuonFilter;
double	triggerresultshelper1_ra2ecaltpFilter;
double	triggerresultshelper1_recovrechitFilter;
double	triggerresultshelper1_scrapingvetoFilter;
double	triggerresultshelper1_trackingfailureFilter;
double	triggerresultshelper1_trackingfailureFilterPFLOW;
double	triggerresultshelper_HLT_HT150_v1;
double	triggerresultshelper_HLT_HT150_v10;
double	triggerresultshelper_HLT_HT150_v10_prs;
double	triggerresultshelper_HLT_HT150_v1_prs;
double	triggerresultshelper_HLT_HT150_v2;
double	triggerresultshelper_HLT_HT150_v2_prs;
double	triggerresultshelper_HLT_HT150_v3;
double	triggerresultshelper_HLT_HT150_v3_prs;
double	triggerresultshelper_HLT_HT150_v4;
double	triggerresultshelper_HLT_HT150_v4_prs;
double	triggerresultshelper_HLT_HT150_v5;
double	triggerresultshelper_HLT_HT150_v5_prs;
double	triggerresultshelper_HLT_HT150_v6;
double	triggerresultshelper_HLT_HT150_v6_prs;
double	triggerresultshelper_HLT_HT150_v7;
double	triggerresultshelper_HLT_HT150_v7_prs;
double	triggerresultshelper_HLT_HT150_v8;
double	triggerresultshelper_HLT_HT150_v8_prs;
double	triggerresultshelper_HLT_HT150_v9;
double	triggerresultshelper_HLT_HT150_v9_prs;
double	triggerresultshelper_HLT_HT2000_v1;
double	triggerresultshelper_HLT_HT2000_v1_prs;
double	triggerresultshelper_HLT_HT2000_v2;
double	triggerresultshelper_HLT_HT2000_v2_prs;
double	triggerresultshelper_HLT_HT2000_v3;
double	triggerresultshelper_HLT_HT2000_v3_prs;
double	triggerresultshelper_HLT_HT2000_v4;
double	triggerresultshelper_HLT_HT2000_v4_prs;
double	triggerresultshelper_HLT_HT200_v1;
double	triggerresultshelper_HLT_HT200_v10;
double	triggerresultshelper_HLT_HT200_v10_prs;
double	triggerresultshelper_HLT_HT200_v1_prs;
double	triggerresultshelper_HLT_HT200_v2;
double	triggerresultshelper_HLT_HT200_v2_prs;
double	triggerresultshelper_HLT_HT200_v3;
double	triggerresultshelper_HLT_HT200_v3_prs;
double	triggerresultshelper_HLT_HT200_v4;
double	triggerresultshelper_HLT_HT200_v4_prs;
double	triggerresultshelper_HLT_HT200_v5;
double	triggerresultshelper_HLT_HT200_v5_prs;
double	triggerresultshelper_HLT_HT200_v6;
double	triggerresultshelper_HLT_HT200_v6_prs;
double	triggerresultshelper_HLT_HT200_v7;
double	triggerresultshelper_HLT_HT200_v7_prs;
double	triggerresultshelper_HLT_HT200_v8;
double	triggerresultshelper_HLT_HT200_v8_prs;
double	triggerresultshelper_HLT_HT200_v9;
double	triggerresultshelper_HLT_HT200_v9_prs;
double	triggerresultshelper_HLT_HT250_MHT100_v2;
double	triggerresultshelper_HLT_HT250_MHT100_v2_prs;
double	triggerresultshelper_HLT_HT250_MHT60_v2;
double	triggerresultshelper_HLT_HT250_MHT60_v2_prs;
double	triggerresultshelper_HLT_HT250_MHT60_v3;
double	triggerresultshelper_HLT_HT250_MHT60_v3_prs;
double	triggerresultshelper_HLT_HT250_MHT70_v1;
double	triggerresultshelper_HLT_HT250_MHT70_v1_prs;
double	triggerresultshelper_HLT_HT250_MHT80_v3;
double	triggerresultshelper_HLT_HT250_MHT80_v3_prs;
double	triggerresultshelper_HLT_HT250_MHT80_v4;
double	triggerresultshelper_HLT_HT250_MHT80_v4_prs;
double	triggerresultshelper_HLT_HT250_MHT90_v1;
double	triggerresultshelper_HLT_HT250_MHT90_v1_prs;
double	triggerresultshelper_HLT_HT250_MHT90_v2;
double	triggerresultshelper_HLT_HT250_MHT90_v2_prs;
double	triggerresultshelper_HLT_HT250_v1;
double	triggerresultshelper_HLT_HT250_v10;
double	triggerresultshelper_HLT_HT250_v10_prs;
double	triggerresultshelper_HLT_HT250_v1_prs;
double	triggerresultshelper_HLT_HT250_v2;
double	triggerresultshelper_HLT_HT250_v2_prs;
double	triggerresultshelper_HLT_HT250_v3;
double	triggerresultshelper_HLT_HT250_v3_prs;
double	triggerresultshelper_HLT_HT250_v4;
double	triggerresultshelper_HLT_HT250_v4_prs;
double	triggerresultshelper_HLT_HT250_v5;
double	triggerresultshelper_HLT_HT250_v5_prs;
double	triggerresultshelper_HLT_HT250_v6;
double	triggerresultshelper_HLT_HT250_v6_prs;
double	triggerresultshelper_HLT_HT250_v7;
double	triggerresultshelper_HLT_HT250_v7_prs;
double	triggerresultshelper_HLT_HT250_v8;
double	triggerresultshelper_HLT_HT250_v8_prs;
double	triggerresultshelper_HLT_HT250_v9;
double	triggerresultshelper_HLT_HT250_v9_prs;
double	triggerresultshelper_HLT_HT260_MHT60_v2;
double	triggerresultshelper_HLT_HT260_MHT60_v2_prs;
double	triggerresultshelper_HLT_HT260_v2;
double	triggerresultshelper_HLT_HT260_v2_prs;
double	triggerresultshelper_HLT_HT260_v3;
double	triggerresultshelper_HLT_HT260_v3_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v2;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v2_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v4;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v4_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v5;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v5_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v7;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v7_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v8;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v8_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT65_v1;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT65_v1_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v7;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v7_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v2;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v2_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v3;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v3_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v4;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v4_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v5;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v5_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v6;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v6_prs;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v7;
double	triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v7_prs;
double	triggerresultshelper_HLT_HT300_MHT75_v2;
double	triggerresultshelper_HLT_HT300_MHT75_v2_prs;
double	triggerresultshelper_HLT_HT300_MHT75_v3;
double	triggerresultshelper_HLT_HT300_MHT75_v3_prs;
double	triggerresultshelper_HLT_HT300_MHT75_v4;
double	triggerresultshelper_HLT_HT300_MHT75_v4_prs;
double	triggerresultshelper_HLT_HT300_MHT75_v7;
double	triggerresultshelper_HLT_HT300_MHT75_v7_prs;
double	triggerresultshelper_HLT_HT300_MHT75_v8;
double	triggerresultshelper_HLT_HT300_MHT75_v8_prs;
double	triggerresultshelper_HLT_HT300_MHT80_v1;
double	triggerresultshelper_HLT_HT300_MHT80_v1_prs;
double	triggerresultshelper_HLT_HT300_MHT80_v2;
double	triggerresultshelper_HLT_HT300_MHT80_v2_prs;
double	triggerresultshelper_HLT_HT300_MHT90_v2;
double	triggerresultshelper_HLT_HT300_MHT90_v2_prs;
double	triggerresultshelper_HLT_HT300_PFMHT55_v2;
double	triggerresultshelper_HLT_HT300_PFMHT55_v2_prs;
double	triggerresultshelper_HLT_HT300_PFMHT55_v3;
double	triggerresultshelper_HLT_HT300_PFMHT55_v3_prs;
double	triggerresultshelper_HLT_HT300_PFMHT55_v4;
double	triggerresultshelper_HLT_HT300_PFMHT55_v4_prs;
double	triggerresultshelper_HLT_HT300_PFMHT55_v5;
double	triggerresultshelper_HLT_HT300_PFMHT55_v5_prs;
double	triggerresultshelper_HLT_HT300_PFMHT55_v7;
double	triggerresultshelper_HLT_HT300_PFMHT55_v7_prs;
double	triggerresultshelper_HLT_HT300_PFMHT55_v8;
double	triggerresultshelper_HLT_HT300_PFMHT55_v8_prs;
double	triggerresultshelper_HLT_HT300_v1;
double	triggerresultshelper_HLT_HT300_v10;
double	triggerresultshelper_HLT_HT300_v10_prs;
double	triggerresultshelper_HLT_HT300_v11;
double	triggerresultshelper_HLT_HT300_v11_prs;
double	triggerresultshelper_HLT_HT300_v1_prs;
double	triggerresultshelper_HLT_HT300_v2;
double	triggerresultshelper_HLT_HT300_v2_prs;
double	triggerresultshelper_HLT_HT300_v3;
double	triggerresultshelper_HLT_HT300_v3_prs;
double	triggerresultshelper_HLT_HT300_v4;
double	triggerresultshelper_HLT_HT300_v4_prs;
double	triggerresultshelper_HLT_HT300_v5;
double	triggerresultshelper_HLT_HT300_v5_prs;
double	triggerresultshelper_HLT_HT300_v6;
double	triggerresultshelper_HLT_HT300_v6_prs;
double	triggerresultshelper_HLT_HT300_v7;
double	triggerresultshelper_HLT_HT300_v7_prs;
double	triggerresultshelper_HLT_HT300_v8;
double	triggerresultshelper_HLT_HT300_v8_prs;
double	triggerresultshelper_HLT_HT300_v9;
double	triggerresultshelper_HLT_HT300_v9_prs;
double	triggerresultshelper_HLT_HT350_MHT70_v1;
double	triggerresultshelper_HLT_HT350_MHT70_v1_prs;
double	triggerresultshelper_HLT_HT350_MHT70_v2;
double	triggerresultshelper_HLT_HT350_MHT70_v2_prs;
double	triggerresultshelper_HLT_HT350_MHT80_v2;
double	triggerresultshelper_HLT_HT350_MHT80_v2_prs;
double	triggerresultshelper_HLT_HT350_MHT90_v1;
double	triggerresultshelper_HLT_HT350_MHT90_v1_prs;
double	triggerresultshelper_HLT_HT350_v1;
double	triggerresultshelper_HLT_HT350_v10;
double	triggerresultshelper_HLT_HT350_v10_prs;
double	triggerresultshelper_HLT_HT350_v1_prs;
double	triggerresultshelper_HLT_HT350_v2;
double	triggerresultshelper_HLT_HT350_v2_prs;
double	triggerresultshelper_HLT_HT350_v3;
double	triggerresultshelper_HLT_HT350_v3_prs;
double	triggerresultshelper_HLT_HT350_v4;
double	triggerresultshelper_HLT_HT350_v4_prs;
double	triggerresultshelper_HLT_HT350_v5;
double	triggerresultshelper_HLT_HT350_v5_prs;
double	triggerresultshelper_HLT_HT350_v6;
double	triggerresultshelper_HLT_HT350_v6_prs;
double	triggerresultshelper_HLT_HT350_v7;
double	triggerresultshelper_HLT_HT350_v7_prs;
double	triggerresultshelper_HLT_HT350_v8;
double	triggerresultshelper_HLT_HT350_v8_prs;
double	triggerresultshelper_HLT_HT350_v9;
double	triggerresultshelper_HLT_HT350_v9_prs;
double	triggerresultshelper_HLT_HT360_v2;
double	triggerresultshelper_HLT_HT360_v2_prs;
double	triggerresultshelper_HLT_HT400_MHT80_v1;
double	triggerresultshelper_HLT_HT400_MHT80_v1_prs;
double	triggerresultshelper_HLT_HT400_v1;
double	triggerresultshelper_HLT_HT400_v10;
double	triggerresultshelper_HLT_HT400_v10_prs;
double	triggerresultshelper_HLT_HT400_v1_prs;
double	triggerresultshelper_HLT_HT400_v2;
double	triggerresultshelper_HLT_HT400_v2_prs;
double	triggerresultshelper_HLT_HT400_v3;
double	triggerresultshelper_HLT_HT400_v3_prs;
double	triggerresultshelper_HLT_HT400_v4;
double	triggerresultshelper_HLT_HT400_v4_prs;
double	triggerresultshelper_HLT_HT400_v5;
double	triggerresultshelper_HLT_HT400_v5_prs;
double	triggerresultshelper_HLT_HT400_v6;
double	triggerresultshelper_HLT_HT400_v6_prs;
double	triggerresultshelper_HLT_HT400_v7;
double	triggerresultshelper_HLT_HT400_v7_prs;
double	triggerresultshelper_HLT_HT400_v8;
double	triggerresultshelper_HLT_HT400_v8_prs;
double	triggerresultshelper_HLT_HT400_v9;
double	triggerresultshelper_HLT_HT400_v9_prs;
double	triggerresultshelper_HLT_HT450_v1;
double	triggerresultshelper_HLT_HT450_v10;
double	triggerresultshelper_HLT_HT450_v10_prs;
double	triggerresultshelper_HLT_HT450_v1_prs;
double	triggerresultshelper_HLT_HT450_v2;
double	triggerresultshelper_HLT_HT450_v2_prs;
double	triggerresultshelper_HLT_HT450_v3;
double	triggerresultshelper_HLT_HT450_v3_prs;
double	triggerresultshelper_HLT_HT450_v4;
double	triggerresultshelper_HLT_HT450_v4_prs;
double	triggerresultshelper_HLT_HT450_v5;
double	triggerresultshelper_HLT_HT450_v5_prs;
double	triggerresultshelper_HLT_HT450_v6;
double	triggerresultshelper_HLT_HT450_v6_prs;
double	triggerresultshelper_HLT_HT450_v7;
double	triggerresultshelper_HLT_HT450_v7_prs;
double	triggerresultshelper_HLT_HT450_v8;
double	triggerresultshelper_HLT_HT450_v8_prs;
double	triggerresultshelper_HLT_HT450_v9;
double	triggerresultshelper_HLT_HT450_v9_prs;
double	triggerresultshelper_HLT_HT500_v10;
double	triggerresultshelper_HLT_HT500_v10_prs;
double	triggerresultshelper_HLT_HT500_v3;
double	triggerresultshelper_HLT_HT500_v3_prs;
double	triggerresultshelper_HLT_HT500_v4;
double	triggerresultshelper_HLT_HT500_v4_prs;
double	triggerresultshelper_HLT_HT500_v5;
double	triggerresultshelper_HLT_HT500_v5_prs;
double	triggerresultshelper_HLT_HT500_v6;
double	triggerresultshelper_HLT_HT500_v6_prs;
double	triggerresultshelper_HLT_HT500_v7;
double	triggerresultshelper_HLT_HT500_v7_prs;
double	triggerresultshelper_HLT_HT500_v8;
double	triggerresultshelper_HLT_HT500_v8_prs;
double	triggerresultshelper_HLT_HT500_v9;
double	triggerresultshelper_HLT_HT500_v9_prs;
double	triggerresultshelper_HLT_HT550_v5;
double	triggerresultshelper_HLT_HT550_v5_prs;
double	triggerresultshelper_HLT_HT550_v6;
double	triggerresultshelper_HLT_HT550_v6_prs;
double	triggerresultshelper_HLT_HT550_v7;
double	triggerresultshelper_HLT_HT550_v7_prs;
double	triggerresultshelper_HLT_HT550_v8;
double	triggerresultshelper_HLT_HT550_v8_prs;
double	triggerresultshelper_HLT_HT600_v1;
double	triggerresultshelper_HLT_HT600_v1_prs;
double	triggerresultshelper_HLT_HT600_v2;
double	triggerresultshelper_HLT_HT600_v2_prs;
double	triggerresultshelper_HLT_HT600_v3;
double	triggerresultshelper_HLT_HT600_v3_prs;
double	triggerresultshelper_HLT_HT650_v1;
double	triggerresultshelper_HLT_HT650_v1_prs;
double	triggerresultshelper_HLT_HT650_v2;
double	triggerresultshelper_HLT_HT650_v2_prs;
double	triggerresultshelper_HLT_HT650_v3;
double	triggerresultshelper_HLT_HT650_v3_prs;
double	triggerresultshelper_HLT_HT750_v1;
double	triggerresultshelper_HLT_HT750_v1_prs;
double	triggerresultshelper_HLT_HT750_v2;
double	triggerresultshelper_HLT_HT750_v2_prs;
double	triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetEta1;
double	triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetEta2;
double	triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetEta3;
double	triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetEta4;
double	triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetPhi1;
double	triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetPhi2;
double	triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetPhi3;
double	triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetPhi4;
double	triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetPt1;
double	triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetPt2;
double	triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetPt3;
double	triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetPt4;
std::vector<double>	vertex_isFake(41,0);
std::vector<double>	vertex_isValid(41,0);
std::vector<double>	vertex_ndof(41,0);
std::vector<double>	vertex_normalizedChi2(41,0);
std::vector<double>	vertex_position_Rho(41,0);
std::vector<double>	vertex_tracksSize(41,0);
std::vector<double>	vertex_x(41,0);
std::vector<double>	vertex_xError(41,0);
std::vector<double>	vertex_y(41,0);
std::vector<double>	vertex_yError(41,0);
std::vector<double>	vertex_z(41,0);
std::vector<double>	vertex_zError(41,0);


//-----------------------------------------------------------------------------
// -- Select variables to be read
//-----------------------------------------------------------------------------
void selectVariables(itreestream& stream)
{
  stream.select("recoBeamSpot_offlineBeamSpot.x0", beamspot_x0);
  stream.select("recoBeamSpot_offlineBeamSpot.y0", beamspot_y0);
  stream.select("recoBeamSpot_offlineBeamSpot.z0", beamspot_z0);
  stream.select("patElectron_selectedPatElectronsPF.caloIso", electron1_caloIso);
  stream.select("patElectron_selectedPatElectronsPF.charge", electron1_charge);
  stream.select("patElectron_selectedPatElectronsPF.chargedHadronIso", electron1_chargedHadronIso);
  stream.select("patElectron_selectedPatElectronsPF.convDcot", electron1_convDcot);
  stream.select("patElectron_selectedPatElectronsPF.convDist", electron1_convDist);
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
  stream.select("patElectron_selectedPatElectronsPF.sigmaIetaIeta", electron1_sigmaIetaIeta);
  stream.select("patElectron_selectedPatElectronsPF.superCluster_eta", electron1_superCluster_eta);
  stream.select("patElectron_selectedPatElectronsPF.trackIso", electron1_trackIso);
  stream.select("patElectron_selectedPatElectronsPF.vertex_z", electron1_vertex_z);
  stream.select("patElectron_selectedPatElectronsPF.vz", electron1_vz);
  stream.select("patElectron_cleanPatElectrons.caloIso", electron_caloIso);
  stream.select("patElectron_cleanPatElectrons.charge", electron_charge);
  stream.select("patElectron_cleanPatElectrons.chargedHadronIso", electron_chargedHadronIso);
  stream.select("patElectron_cleanPatElectrons.convDcot", electron_convDcot);
  stream.select("patElectron_cleanPatElectrons.convDist", electron_convDist);
  stream.select("patElectron_cleanPatElectrons.dB", electron_dB);
  stream.select("patElectron_cleanPatElectrons.deltaEtaSuperClusterTrackAtVtx", electron_deltaEtaSuperClusterTrackAtVtx);
  stream.select("patElectron_cleanPatElectrons.deltaPhiSuperClusterTrackAtVtx", electron_deltaPhiSuperClusterTrackAtVtx);
  stream.select("patElectron_cleanPatElectrons.dr03EcalRecHitSumEt", electron_dr03EcalRecHitSumEt);
  stream.select("patElectron_cleanPatElectrons.dr03HcalTowerSumEt", electron_dr03HcalTowerSumEt);
  stream.select("patElectron_cleanPatElectrons.dr03TkSumPt", electron_dr03TkSumPt);
  stream.select("patElectron_cleanPatElectrons.ecalIso", electron_ecalIso);
  stream.select("patElectron_cleanPatElectrons.eidRobustTight", electron_eidRobustTight);
  stream.select("patElectron_cleanPatElectrons.energy", electron_energy);
  stream.select("patElectron_cleanPatElectrons.et", electron_et);
  stream.select("patElectron_cleanPatElectrons.eta", electron_eta);
  stream.select("patElectron_cleanPatElectrons.genLepton_eta", electron_genLepton_eta);
  stream.select("patElectron_cleanPatElectrons.genLepton_pdgId", electron_genLepton_pdgId);
  stream.select("patElectron_cleanPatElectrons.genLepton_phi", electron_genLepton_phi);
  stream.select("patElectron_cleanPatElectrons.genLepton_pt", electron_genLepton_pt);
  stream.select("patElectron_cleanPatElectrons.gsfTrack_d0", electron_gsfTrack_d0);
  stream.select("patElectron_cleanPatElectrons.gsfTrack_trackerExpectedHitsInner_numberOfLostHits", electron_gsfTrack_trackerExpectedHitsInner_numberOfLostHits);
  stream.select("patElectron_cleanPatElectrons.hadronicOverEm", electron_hadronicOverEm);
  stream.select("patElectron_cleanPatElectrons.hcalIso", electron_hcalIso);
  stream.select("patElectron_cleanPatElectrons.neutralHadronIso", electron_neutralHadronIso);
  stream.select("patElectron_cleanPatElectrons.phi", electron_phi);
  stream.select("patElectron_cleanPatElectrons.photonIso", electron_photonIso);
  stream.select("patElectron_cleanPatElectrons.pt", electron_pt);
  stream.select("patElectron_cleanPatElectrons.sigmaIetaIeta", electron_sigmaIetaIeta);
  stream.select("patElectron_cleanPatElectrons.simpleEleId80cIso", electron_simpleEleId80cIso);
  stream.select("patElectron_cleanPatElectrons.simpleEleId80relIso", electron_simpleEleId80relIso);
  stream.select("patElectron_cleanPatElectrons.simpleEleId95cIso", electron_simpleEleId95cIso);
  stream.select("patElectron_cleanPatElectrons.simpleEleId95relIso", electron_simpleEleId95relIso);
  stream.select("patElectron_cleanPatElectrons.superCluster_eta", electron_superCluster_eta);
  stream.select("patElectron_cleanPatElectrons.trackIso", electron_trackIso);
  stream.select("patElectron_cleanPatElectrons.vertex_z", electron_vertex_z);
  stream.select("patElectron_cleanPatElectrons.vz", electron_vz);
  stream.select("patElectronHelper_selectedPatElectronsPF.dxywrtBeamSpot", electronhelper1_dxywrtBeamSpot);
  stream.select("patElectronHelper_cleanPatElectrons.dxywrtBeamSpot", electronhelper_dxywrtBeamSpot);
  stream.select("edmEventHelper_info.bunchCrossing", eventhelper_bunchCrossing);
  stream.select("edmEventHelper_info.event", eventhelper_event);
  stream.select("edmEventHelper_info.isRealData", eventhelper_isRealData);
  stream.select("edmEventHelper_info.luminosityBlock", eventhelper_luminosityBlock);
  stream.select("edmEventHelper_info.run", eventhelper_run);
  stream.select("edmEventHelperExtra_info.trackMPT1phi", eventhelperextra_trackMPT1phi);
  stream.select("edmEventHelperExtra_info.trackMPT1pt", eventhelperextra_trackMPT1pt);
  stream.select("edmEventHelperExtra_info.trackMPT5phi", eventhelperextra_trackMPT5phi);
  stream.select("edmEventHelperExtra_info.trackMPT5pt", eventhelperextra_trackMPT5pt);
  stream.select("edmEventLHEHelperExtra_info.m0", eventlhehelperextra_m0);
  stream.select("edmEventLHEHelperExtra_info.m12", eventlhehelperextra_m12);
  stream.select("edmEventLHEHelperExtra_info.mGL", eventlhehelperextra_mGL);
  stream.select("edmEventLHEHelperExtra_info.mLSP", eventlhehelperextra_mLSP);
  stream.select("edmEventLHEHelperExtra_info.xCHI", eventlhehelperextra_xCHI);
  stream.select("GenEventInfoProduct_generator.pdf1", geneventinfoproduct_pdf1);
  stream.select("GenEventInfoProduct_generator.pdf2", geneventinfoproduct_pdf2);
  stream.select("GenEventInfoProduct_generator.scalePDF", geneventinfoproduct_scalePDF);
  stream.select("GenEventInfoProduct_generator.weight", geneventinfoproduct_weight);
  stream.select("GenEventInfoProduct_generator.x1", geneventinfoproduct_x1);
  stream.select("GenEventInfoProduct_generator.x2", geneventinfoproduct_x2);
  stream.select("GenEventInfoProductHelper_generator_cteq.pdf1", geneventinfoproducthelper1_pdf1);
  stream.select("GenEventInfoProductHelper_generator_cteq.pdf2", geneventinfoproducthelper1_pdf2);
  stream.select("GenEventInfoProductHelper_generator_cteq.pdfweight", geneventinfoproducthelper1_pdfweight);
  stream.select("GenEventInfoProductHelper_generator_cteq.pdfweightsum", geneventinfoproducthelper1_pdfweightsum);
  stream.select("GenEventInfoProductHelper_generator_mstw.pdf1", geneventinfoproducthelper2_pdf1);
  stream.select("GenEventInfoProductHelper_generator_mstw.pdf2", geneventinfoproducthelper2_pdf2);
  stream.select("GenEventInfoProductHelper_generator_mstw.pdfweight", geneventinfoproducthelper2_pdfweight);
  stream.select("GenEventInfoProductHelper_generator_mstw.pdfweightsum", geneventinfoproducthelper2_pdfweightsum);
  stream.select("GenEventInfoProductHelper_generator_NNPDF.pdf1", geneventinfoproducthelper_pdf1);
  stream.select("GenEventInfoProductHelper_generator_NNPDF.pdf2", geneventinfoproducthelper_pdf2);
  stream.select("GenEventInfoProductHelper_generator_NNPDF.pdfweight", geneventinfoproducthelper_pdfweight);
  stream.select("GenEventInfoProductHelper_generator_NNPDF.pdfweightsum", geneventinfoproducthelper_pdfweightsum);
  stream.select("recoGenParticleHelperRA2b_genParticles.charge", genparticlehelperra2_charge);
  stream.select("recoGenParticleHelperRA2b_genParticles.eta", genparticlehelperra2_eta);
  stream.select("recoGenParticleHelperRA2b_genParticles.firstDaughter", genparticlehelperra2_firstDaughter);
  stream.select("recoGenParticleHelperRA2b_genParticles.firstMother", genparticlehelperra2_firstMother);
  stream.select("recoGenParticleHelperRA2b_genParticles.lastDaughter", genparticlehelperra2_lastDaughter);
  stream.select("recoGenParticleHelperRA2b_genParticles.lastMother", genparticlehelperra2_lastMother);
  stream.select("recoGenParticleHelperRA2b_genParticles.mass", genparticlehelperra2_mass);
  stream.select("recoGenParticleHelperRA2b_genParticles.pdgId", genparticlehelperra2_pdgId);
  stream.select("recoGenParticleHelperRA2b_genParticles.phi", genparticlehelperra2_phi);
  stream.select("recoGenParticleHelperRA2b_genParticles.pt", genparticlehelperra2_pt);
  stream.select("recoGenParticleHelperRA2b_genParticles.status", genparticlehelperra2_status);
  stream.select("GenRunInfoProduct_generator.externalXSecLO_error", genruninfoproduct_externalXSecLO_error);
  stream.select("GenRunInfoProduct_generator.externalXSecLO_value", genruninfoproduct_externalXSecLO_value);
  stream.select("GenRunInfoProduct_generator.filterEfficiency", genruninfoproduct_filterEfficiency);
  stream.select("GenRunInfoProduct_generator.internalXSec_error", genruninfoproduct_internalXSec_error);
  stream.select("GenRunInfoProduct_generator.internalXSec_value", genruninfoproduct_internalXSec_value);
  stream.select("patJet_selectedPatJetsPF.HFEMEnergy", jet1_HFEMEnergy);
  stream.select("patJet_selectedPatJetsPF.HFEMMultiplicity", jet1_HFEMMultiplicity);
  stream.select("patJet_selectedPatJetsPF.HFHadronEnergy", jet1_HFHadronEnergy);
  stream.select("patJet_selectedPatJetsPF.HFHadronMultiplicity", jet1_HFHadronMultiplicity);
  stream.select("patJet_selectedPatJetsPF.chargedEmEnergyFraction", jet1_chargedEmEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.chargedHadronEnergyFraction", jet1_chargedHadronEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.chargedHadronMultiplicity", jet1_chargedHadronMultiplicity);
  stream.select("patJet_selectedPatJetsPF.chargedMuEnergy", jet1_chargedMuEnergy);
  stream.select("patJet_selectedPatJetsPF.chargedMultiplicity", jet1_chargedMultiplicity);
  stream.select("patJet_selectedPatJetsPF.combinedSecondaryVertexBJetTags", jet1_combinedSecondaryVertexBJetTags);
  stream.select("patJet_selectedPatJetsPF.combinedSecondaryVertexMVABJetTags", jet1_combinedSecondaryVertexMVABJetTags);
  stream.select("patJet_selectedPatJetsPF.electronEnergy", jet1_electronEnergy);
  stream.select("patJet_selectedPatJetsPF.electronMultiplicity", jet1_electronMultiplicity);
  stream.select("patJet_selectedPatJetsPF.energy", jet1_energy);
  stream.select("patJet_selectedPatJetsPF.et", jet1_et);
  stream.select("patJet_selectedPatJetsPF.eta", jet1_eta);
  stream.select("patJet_selectedPatJetsPF.genJet_energy", jet1_genJet_energy);
  stream.select("patJet_selectedPatJetsPF.genJet_eta", jet1_genJet_eta);
  stream.select("patJet_selectedPatJetsPF.genJet_invisibleEnergy", jet1_genJet_invisibleEnergy);
  stream.select("patJet_selectedPatJetsPF.genJet_phi", jet1_genJet_phi);
  stream.select("patJet_selectedPatJetsPF.genJet_pt", jet1_genJet_pt);
  stream.select("patJet_selectedPatJetsPF.genParton_energy", jet1_genParton_energy);
  stream.select("patJet_selectedPatJetsPF.genParton_eta", jet1_genParton_eta);
  stream.select("patJet_selectedPatJetsPF.genParton_pdgId", jet1_genParton_pdgId);
  stream.select("patJet_selectedPatJetsPF.genParton_phi", jet1_genParton_phi);
  stream.select("patJet_selectedPatJetsPF.genParton_pt", jet1_genParton_pt);
  stream.select("patJet_selectedPatJetsPF.jetArea", jet1_jetArea);
  stream.select("patJet_selectedPatJetsPF.jetBProbabilityBJetTags", jet1_jetBProbabilityBJetTags);
  stream.select("patJet_selectedPatJetsPF.jetID_fHPD", jet1_jetID_fHPD);
  stream.select("patJet_selectedPatJetsPF.jetID_n90Hits", jet1_jetID_n90Hits);
  stream.select("patJet_selectedPatJetsPF.jetProbabilityBJetTags", jet1_jetProbabilityBJetTags);
  stream.select("patJet_selectedPatJetsPF.muonEnergy", jet1_muonEnergy);
  stream.select("patJet_selectedPatJetsPF.neutralEmEnergyFraction", jet1_neutralEmEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.neutralHadronEnergy", jet1_neutralHadronEnergy);
  stream.select("patJet_selectedPatJetsPF.neutralHadronEnergyFraction", jet1_neutralHadronEnergyFraction);
  stream.select("patJet_selectedPatJetsPF.neutralHadronMultiplicity", jet1_neutralHadronMultiplicity);
  stream.select("patJet_selectedPatJetsPF.neutralMultiplicity", jet1_neutralMultiplicity);
  stream.select("patJet_selectedPatJetsPF.numberOfDaughters", jet1_numberOfDaughters);
  stream.select("patJet_selectedPatJetsPF.partonFlavour", jet1_partonFlavour);
  stream.select("patJet_selectedPatJetsPF.phi", jet1_phi);
  stream.select("patJet_selectedPatJetsPF.photonEnergy", jet1_photonEnergy);
  stream.select("patJet_selectedPatJetsPF.photonMultiplicity", jet1_photonMultiplicity);
  stream.select("patJet_selectedPatJetsPF.pt", jet1_pt);
  stream.select("patJet_selectedPatJetsPF.simpleSecondaryVertexBJetTags", jet1_simpleSecondaryVertexBJetTags);
  stream.select("patJet_selectedPatJetsPF.simpleSecondaryVertexHighEffBJetTags", jet1_simpleSecondaryVertexHighEffBJetTags);
  stream.select("patJet_selectedPatJetsPF.simpleSecondaryVertexHighPurBJetTags", jet1_simpleSecondaryVertexHighPurBJetTags);
  stream.select("patJet_selectedPatJetsPF.softElectronByIP3dBJetTags", jet1_softElectronByIP3dBJetTags);
  stream.select("patJet_selectedPatJetsPF.softElectronByPtBJetTags", jet1_softElectronByPtBJetTags);
  stream.select("patJet_selectedPatJetsPF.softMuonBJetTags", jet1_softMuonBJetTags);
  stream.select("patJet_selectedPatJetsPF.softMuonByIP3dBJetTags", jet1_softMuonByIP3dBJetTags);
  stream.select("patJet_selectedPatJetsPF.softMuonByPtBJetTags", jet1_softMuonByPtBJetTags);
  stream.select("patJet_selectedPatJetsPF.trackCountingHighEffBJetTags", jet1_trackCountingHighEffBJetTags);
  stream.select("patJet_selectedPatJetsPF.trackCountingHighPurBJetTags", jet1_trackCountingHighPurBJetTags);
  stream.select("patJet_selectedPatJetsPF.uncor_energy", jet1_uncor_energy);
  stream.select("patJet_selectedPatJetsPF.uncor_et", jet1_uncor_et);
  stream.select("patJet_selectedPatJetsPF.uncor_eta", jet1_uncor_eta);
  stream.select("patJet_selectedPatJetsPF.uncor_phi", jet1_uncor_phi);
  stream.select("patJet_selectedPatJetsPF.uncor_pt", jet1_uncor_pt);
  stream.select("patJet_selectedPatJetsPFLOW.HFEMEnergy", jet2_HFEMEnergy);
  stream.select("patJet_selectedPatJetsPFLOW.HFEMMultiplicity", jet2_HFEMMultiplicity);
  stream.select("patJet_selectedPatJetsPFLOW.HFHadronEnergy", jet2_HFHadronEnergy);
  stream.select("patJet_selectedPatJetsPFLOW.HFHadronMultiplicity", jet2_HFHadronMultiplicity);
  stream.select("patJet_selectedPatJetsPFLOW.chargedEmEnergyFraction", jet2_chargedEmEnergyFraction);
  stream.select("patJet_selectedPatJetsPFLOW.chargedHadronEnergyFraction", jet2_chargedHadronEnergyFraction);
  stream.select("patJet_selectedPatJetsPFLOW.chargedHadronMultiplicity", jet2_chargedHadronMultiplicity);
  stream.select("patJet_selectedPatJetsPFLOW.chargedMuEnergy", jet2_chargedMuEnergy);
  stream.select("patJet_selectedPatJetsPFLOW.chargedMultiplicity", jet2_chargedMultiplicity);
  stream.select("patJet_selectedPatJetsPFLOW.combinedSecondaryVertexBJetTags", jet2_combinedSecondaryVertexBJetTags);
  stream.select("patJet_selectedPatJetsPFLOW.combinedSecondaryVertexMVABJetTags", jet2_combinedSecondaryVertexMVABJetTags);
  stream.select("patJet_selectedPatJetsPFLOW.electronEnergy", jet2_electronEnergy);
  stream.select("patJet_selectedPatJetsPFLOW.electronMultiplicity", jet2_electronMultiplicity);
  stream.select("patJet_selectedPatJetsPFLOW.energy", jet2_energy);
  stream.select("patJet_selectedPatJetsPFLOW.et", jet2_et);
  stream.select("patJet_selectedPatJetsPFLOW.eta", jet2_eta);
  stream.select("patJet_selectedPatJetsPFLOW.genJet_energy", jet2_genJet_energy);
  stream.select("patJet_selectedPatJetsPFLOW.genJet_eta", jet2_genJet_eta);
  stream.select("patJet_selectedPatJetsPFLOW.genJet_invisibleEnergy", jet2_genJet_invisibleEnergy);
  stream.select("patJet_selectedPatJetsPFLOW.genJet_phi", jet2_genJet_phi);
  stream.select("patJet_selectedPatJetsPFLOW.genJet_pt", jet2_genJet_pt);
  stream.select("patJet_selectedPatJetsPFLOW.genParton_energy", jet2_genParton_energy);
  stream.select("patJet_selectedPatJetsPFLOW.genParton_eta", jet2_genParton_eta);
  stream.select("patJet_selectedPatJetsPFLOW.genParton_pdgId", jet2_genParton_pdgId);
  stream.select("patJet_selectedPatJetsPFLOW.genParton_phi", jet2_genParton_phi);
  stream.select("patJet_selectedPatJetsPFLOW.genParton_pt", jet2_genParton_pt);
  stream.select("patJet_selectedPatJetsPFLOW.jetArea", jet2_jetArea);
  stream.select("patJet_selectedPatJetsPFLOW.jetBProbabilityBJetTags", jet2_jetBProbabilityBJetTags);
  stream.select("patJet_selectedPatJetsPFLOW.jetID_fHPD", jet2_jetID_fHPD);
  stream.select("patJet_selectedPatJetsPFLOW.jetID_n90Hits", jet2_jetID_n90Hits);
  stream.select("patJet_selectedPatJetsPFLOW.jetProbabilityBJetTags", jet2_jetProbabilityBJetTags);
  stream.select("patJet_selectedPatJetsPFLOW.muonEnergy", jet2_muonEnergy);
  stream.select("patJet_selectedPatJetsPFLOW.neutralEmEnergyFraction", jet2_neutralEmEnergyFraction);
  stream.select("patJet_selectedPatJetsPFLOW.neutralHadronEnergy", jet2_neutralHadronEnergy);
  stream.select("patJet_selectedPatJetsPFLOW.neutralHadronEnergyFraction", jet2_neutralHadronEnergyFraction);
  stream.select("patJet_selectedPatJetsPFLOW.neutralHadronMultiplicity", jet2_neutralHadronMultiplicity);
  stream.select("patJet_selectedPatJetsPFLOW.neutralMultiplicity", jet2_neutralMultiplicity);
  stream.select("patJet_selectedPatJetsPFLOW.numberOfDaughters", jet2_numberOfDaughters);
  stream.select("patJet_selectedPatJetsPFLOW.partonFlavour", jet2_partonFlavour);
  stream.select("patJet_selectedPatJetsPFLOW.phi", jet2_phi);
  stream.select("patJet_selectedPatJetsPFLOW.photonEnergy", jet2_photonEnergy);
  stream.select("patJet_selectedPatJetsPFLOW.photonMultiplicity", jet2_photonMultiplicity);
  stream.select("patJet_selectedPatJetsPFLOW.pt", jet2_pt);
  stream.select("patJet_selectedPatJetsPFLOW.simpleSecondaryVertexBJetTags", jet2_simpleSecondaryVertexBJetTags);
  stream.select("patJet_selectedPatJetsPFLOW.simpleSecondaryVertexHighEffBJetTags", jet2_simpleSecondaryVertexHighEffBJetTags);
  stream.select("patJet_selectedPatJetsPFLOW.simpleSecondaryVertexHighPurBJetTags", jet2_simpleSecondaryVertexHighPurBJetTags);
  stream.select("patJet_selectedPatJetsPFLOW.softElectronByIP3dBJetTags", jet2_softElectronByIP3dBJetTags);
  stream.select("patJet_selectedPatJetsPFLOW.softElectronByPtBJetTags", jet2_softElectronByPtBJetTags);
  stream.select("patJet_selectedPatJetsPFLOW.softMuonBJetTags", jet2_softMuonBJetTags);
  stream.select("patJet_selectedPatJetsPFLOW.softMuonByIP3dBJetTags", jet2_softMuonByIP3dBJetTags);
  stream.select("patJet_selectedPatJetsPFLOW.softMuonByPtBJetTags", jet2_softMuonByPtBJetTags);
  stream.select("patJet_selectedPatJetsPFLOW.trackCountingHighEffBJetTags", jet2_trackCountingHighEffBJetTags);
  stream.select("patJet_selectedPatJetsPFLOW.trackCountingHighPurBJetTags", jet2_trackCountingHighPurBJetTags);
  stream.select("patJet_selectedPatJetsPFLOW.uncor_energy", jet2_uncor_energy);
  stream.select("patJet_selectedPatJetsPFLOW.uncor_et", jet2_uncor_et);
  stream.select("patJet_selectedPatJetsPFLOW.uncor_eta", jet2_uncor_eta);
  stream.select("patJet_selectedPatJetsPFLOW.uncor_phi", jet2_uncor_phi);
  stream.select("patJet_selectedPatJetsPFLOW.uncor_pt", jet2_uncor_pt);
  stream.select("patJet_cleanPatJetsAK5PF.HFEMEnergy", jet_HFEMEnergy);
  stream.select("patJet_cleanPatJetsAK5PF.HFEMMultiplicity", jet_HFEMMultiplicity);
  stream.select("patJet_cleanPatJetsAK5PF.HFHadronEnergy", jet_HFHadronEnergy);
  stream.select("patJet_cleanPatJetsAK5PF.HFHadronMultiplicity", jet_HFHadronMultiplicity);
  stream.select("patJet_cleanPatJetsAK5PF.chargedEmEnergyFraction", jet_chargedEmEnergyFraction);
  stream.select("patJet_cleanPatJetsAK5PF.chargedHadronEnergyFraction", jet_chargedHadronEnergyFraction);
  stream.select("patJet_cleanPatJetsAK5PF.chargedHadronMultiplicity", jet_chargedHadronMultiplicity);
  stream.select("patJet_cleanPatJetsAK5PF.chargedMuEnergy", jet_chargedMuEnergy);
  stream.select("patJet_cleanPatJetsAK5PF.chargedMultiplicity", jet_chargedMultiplicity);
  stream.select("patJet_cleanPatJetsAK5PF.combinedSecondaryVertexBJetTags", jet_combinedSecondaryVertexBJetTags);
  stream.select("patJet_cleanPatJetsAK5PF.combinedSecondaryVertexMVABJetTags", jet_combinedSecondaryVertexMVABJetTags);
  stream.select("patJet_cleanPatJetsAK5PF.electronEnergy", jet_electronEnergy);
  stream.select("patJet_cleanPatJetsAK5PF.electronMultiplicity", jet_electronMultiplicity);
  stream.select("patJet_cleanPatJetsAK5PF.energy", jet_energy);
  stream.select("patJet_cleanPatJetsAK5PF.et", jet_et);
  stream.select("patJet_cleanPatJetsAK5PF.eta", jet_eta);
  stream.select("patJet_cleanPatJetsAK5PF.genJet_energy", jet_genJet_energy);
  stream.select("patJet_cleanPatJetsAK5PF.genJet_eta", jet_genJet_eta);
  stream.select("patJet_cleanPatJetsAK5PF.genJet_invisibleEnergy", jet_genJet_invisibleEnergy);
  stream.select("patJet_cleanPatJetsAK5PF.genJet_phi", jet_genJet_phi);
  stream.select("patJet_cleanPatJetsAK5PF.genJet_pt", jet_genJet_pt);
  stream.select("patJet_cleanPatJetsAK5PF.genParton_energy", jet_genParton_energy);
  stream.select("patJet_cleanPatJetsAK5PF.genParton_eta", jet_genParton_eta);
  stream.select("patJet_cleanPatJetsAK5PF.genParton_pdgId", jet_genParton_pdgId);
  stream.select("patJet_cleanPatJetsAK5PF.genParton_phi", jet_genParton_phi);
  stream.select("patJet_cleanPatJetsAK5PF.genParton_pt", jet_genParton_pt);
  stream.select("patJet_cleanPatJetsAK5PF.jetArea", jet_jetArea);
  stream.select("patJet_cleanPatJetsAK5PF.jetBProbabilityBJetTags", jet_jetBProbabilityBJetTags);
  stream.select("patJet_cleanPatJetsAK5PF.jetID_fHPD", jet_jetID_fHPD);
  stream.select("patJet_cleanPatJetsAK5PF.jetID_n90Hits", jet_jetID_n90Hits);
  stream.select("patJet_cleanPatJetsAK5PF.jetProbabilityBJetTags", jet_jetProbabilityBJetTags);
  stream.select("patJet_cleanPatJetsAK5PF.muonEnergy", jet_muonEnergy);
  stream.select("patJet_cleanPatJetsAK5PF.neutralEmEnergyFraction", jet_neutralEmEnergyFraction);
  stream.select("patJet_cleanPatJetsAK5PF.neutralHadronEnergy", jet_neutralHadronEnergy);
  stream.select("patJet_cleanPatJetsAK5PF.neutralHadronEnergyFraction", jet_neutralHadronEnergyFraction);
  stream.select("patJet_cleanPatJetsAK5PF.neutralHadronMultiplicity", jet_neutralHadronMultiplicity);
  stream.select("patJet_cleanPatJetsAK5PF.neutralMultiplicity", jet_neutralMultiplicity);
  stream.select("patJet_cleanPatJetsAK5PF.numberOfDaughters", jet_numberOfDaughters);
  stream.select("patJet_cleanPatJetsAK5PF.partonFlavour", jet_partonFlavour);
  stream.select("patJet_cleanPatJetsAK5PF.phi", jet_phi);
  stream.select("patJet_cleanPatJetsAK5PF.photonEnergy", jet_photonEnergy);
  stream.select("patJet_cleanPatJetsAK5PF.photonMultiplicity", jet_photonMultiplicity);
  stream.select("patJet_cleanPatJetsAK5PF.pt", jet_pt);
  stream.select("patJet_cleanPatJetsAK5PF.simpleSecondaryVertexBJetTags", jet_simpleSecondaryVertexBJetTags);
  stream.select("patJet_cleanPatJetsAK5PF.simpleSecondaryVertexHighEffBJetTags", jet_simpleSecondaryVertexHighEffBJetTags);
  stream.select("patJet_cleanPatJetsAK5PF.simpleSecondaryVertexHighPurBJetTags", jet_simpleSecondaryVertexHighPurBJetTags);
  stream.select("patJet_cleanPatJetsAK5PF.softElectronByIP3dBJetTags", jet_softElectronByIP3dBJetTags);
  stream.select("patJet_cleanPatJetsAK5PF.softElectronByPtBJetTags", jet_softElectronByPtBJetTags);
  stream.select("patJet_cleanPatJetsAK5PF.softMuonBJetTags", jet_softMuonBJetTags);
  stream.select("patJet_cleanPatJetsAK5PF.softMuonByIP3dBJetTags", jet_softMuonByIP3dBJetTags);
  stream.select("patJet_cleanPatJetsAK5PF.softMuonByPtBJetTags", jet_softMuonByPtBJetTags);
  stream.select("patJet_cleanPatJetsAK5PF.trackCountingHighEffBJetTags", jet_trackCountingHighEffBJetTags);
  stream.select("patJet_cleanPatJetsAK5PF.trackCountingHighPurBJetTags", jet_trackCountingHighPurBJetTags);
  stream.select("patJet_cleanPatJetsAK5PF.uncor_energy", jet_uncor_energy);
  stream.select("patJet_cleanPatJetsAK5PF.uncor_et", jet_uncor_et);
  stream.select("patJet_cleanPatJetsAK5PF.uncor_eta", jet_uncor_eta);
  stream.select("patJet_cleanPatJetsAK5PF.uncor_phi", jet_uncor_phi);
  stream.select("patJet_cleanPatJetsAK5PF.uncor_pt", jet_uncor_pt);
  stream.select("patJetHelper_selectedPatJetsPF.jecFactor", jethelper1_jecFactor);
  stream.select("patJetHelper_selectedPatJetsPF.jecFactorNoL1Fast", jethelper1_jecFactorNoL1Fast);
  stream.select("patJetHelper_selectedPatJetsPF.jetUncMinus", jethelper1_jetUncMinus);
  stream.select("patJetHelper_selectedPatJetsPF.jetUncPlus", jethelper1_jetUncPlus);
  stream.select("patJetHelper_selectedPatJetsPF.nSecondaryVertices", jethelper1_nSecondaryVertices);
  stream.select("patJetHelper_selectedPatJetsPF.secondaryVertex0DirectionX", jethelper1_secondaryVertex0DirectionX);
  stream.select("patJetHelper_selectedPatJetsPF.secondaryVertex0DirectionY", jethelper1_secondaryVertex0DirectionY);
  stream.select("patJetHelper_selectedPatJetsPF.secondaryVertex0DirectionZ", jethelper1_secondaryVertex0DirectionZ);
  stream.select("patJetHelper_selectedPatJetsPF.secondaryVertex0FlightDistance", jethelper1_secondaryVertex0FlightDistance);
  stream.select("patJetHelper_selectedPatJetsPF.secondaryVertex0FlightDistanceSig", jethelper1_secondaryVertex0FlightDistanceSig);
  stream.select("patJetHelper_selectedPatJetsPF.secondaryVertex0Mass", jethelper1_secondaryVertex0Mass);
  stream.select("patJetHelper_selectedPatJetsPF.secondaryVertex1DirectionX", jethelper1_secondaryVertex1DirectionX);
  stream.select("patJetHelper_selectedPatJetsPF.secondaryVertex1DirectionY", jethelper1_secondaryVertex1DirectionY);
  stream.select("patJetHelper_selectedPatJetsPF.secondaryVertex1DirectionZ", jethelper1_secondaryVertex1DirectionZ);
  stream.select("patJetHelper_selectedPatJetsPF.secondaryVertex1FlightDistance", jethelper1_secondaryVertex1FlightDistance);
  stream.select("patJetHelper_selectedPatJetsPF.secondaryVertex1FlightDistanceSig", jethelper1_secondaryVertex1FlightDistanceSig);
  stream.select("patJetHelper_selectedPatJetsPF.secondaryVertex1Mass", jethelper1_secondaryVertex1Mass);
  stream.select("patJetHelper_selectedPatJetsPFLOW.jecFactor", jethelper2_jecFactor);
  stream.select("patJetHelper_selectedPatJetsPFLOW.jecFactorNoL1Fast", jethelper2_jecFactorNoL1Fast);
  stream.select("patJetHelper_selectedPatJetsPFLOW.jetUncMinus", jethelper2_jetUncMinus);
  stream.select("patJetHelper_selectedPatJetsPFLOW.jetUncPlus", jethelper2_jetUncPlus);
  stream.select("patJetHelper_selectedPatJetsPFLOW.nSecondaryVertices", jethelper2_nSecondaryVertices);
  stream.select("patJetHelper_selectedPatJetsPFLOW.secondaryVertex0DirectionX", jethelper2_secondaryVertex0DirectionX);
  stream.select("patJetHelper_selectedPatJetsPFLOW.secondaryVertex0DirectionY", jethelper2_secondaryVertex0DirectionY);
  stream.select("patJetHelper_selectedPatJetsPFLOW.secondaryVertex0DirectionZ", jethelper2_secondaryVertex0DirectionZ);
  stream.select("patJetHelper_selectedPatJetsPFLOW.secondaryVertex0FlightDistance", jethelper2_secondaryVertex0FlightDistance);
  stream.select("patJetHelper_selectedPatJetsPFLOW.secondaryVertex0FlightDistanceSig", jethelper2_secondaryVertex0FlightDistanceSig);
  stream.select("patJetHelper_selectedPatJetsPFLOW.secondaryVertex0Mass", jethelper2_secondaryVertex0Mass);
  stream.select("patJetHelper_selectedPatJetsPFLOW.secondaryVertex1DirectionX", jethelper2_secondaryVertex1DirectionX);
  stream.select("patJetHelper_selectedPatJetsPFLOW.secondaryVertex1DirectionY", jethelper2_secondaryVertex1DirectionY);
  stream.select("patJetHelper_selectedPatJetsPFLOW.secondaryVertex1DirectionZ", jethelper2_secondaryVertex1DirectionZ);
  stream.select("patJetHelper_selectedPatJetsPFLOW.secondaryVertex1FlightDistance", jethelper2_secondaryVertex1FlightDistance);
  stream.select("patJetHelper_selectedPatJetsPFLOW.secondaryVertex1FlightDistanceSig", jethelper2_secondaryVertex1FlightDistanceSig);
  stream.select("patJetHelper_selectedPatJetsPFLOW.secondaryVertex1Mass", jethelper2_secondaryVertex1Mass);
  stream.select("patJetHelper_cleanPatJetsAK5PF.jecFactor", jethelper_jecFactor);
  stream.select("patJetHelper_cleanPatJetsAK5PF.jecFactorNoL1Fast", jethelper_jecFactorNoL1Fast);
  stream.select("patJetHelper_cleanPatJetsAK5PF.jetUncMinus", jethelper_jetUncMinus);
  stream.select("patJetHelper_cleanPatJetsAK5PF.jetUncPlus", jethelper_jetUncPlus);
  stream.select("patJetHelper_cleanPatJetsAK5PF.nSecondaryVertices", jethelper_nSecondaryVertices);
  stream.select("patJetHelper_cleanPatJetsAK5PF.secondaryVertex0DirectionX", jethelper_secondaryVertex0DirectionX);
  stream.select("patJetHelper_cleanPatJetsAK5PF.secondaryVertex0DirectionY", jethelper_secondaryVertex0DirectionY);
  stream.select("patJetHelper_cleanPatJetsAK5PF.secondaryVertex0DirectionZ", jethelper_secondaryVertex0DirectionZ);
  stream.select("patJetHelper_cleanPatJetsAK5PF.secondaryVertex0FlightDistance", jethelper_secondaryVertex0FlightDistance);
  stream.select("patJetHelper_cleanPatJetsAK5PF.secondaryVertex0FlightDistanceSig", jethelper_secondaryVertex0FlightDistanceSig);
  stream.select("patJetHelper_cleanPatJetsAK5PF.secondaryVertex0Mass", jethelper_secondaryVertex0Mass);
  stream.select("patJetHelper_cleanPatJetsAK5PF.secondaryVertex1DirectionX", jethelper_secondaryVertex1DirectionX);
  stream.select("patJetHelper_cleanPatJetsAK5PF.secondaryVertex1DirectionY", jethelper_secondaryVertex1DirectionY);
  stream.select("patJetHelper_cleanPatJetsAK5PF.secondaryVertex1DirectionZ", jethelper_secondaryVertex1DirectionZ);
  stream.select("patJetHelper_cleanPatJetsAK5PF.secondaryVertex1FlightDistance", jethelper_secondaryVertex1FlightDistance);
  stream.select("patJetHelper_cleanPatJetsAK5PF.secondaryVertex1FlightDistanceSig", jethelper_secondaryVertex1FlightDistanceSig);
  stream.select("patJetHelper_cleanPatJetsAK5PF.secondaryVertex1Mass", jethelper_secondaryVertex1Mass);
  stream.select("patMET_patMETsPF.energy", met1_energy);
  stream.select("patMET_patMETsPF.et", met1_et);
  stream.select("patMET_patMETsPF.genMET_energy", met1_genMET_energy);
  stream.select("patMET_patMETsPF.genMET_et", met1_genMET_et);
  stream.select("patMET_patMETsPF.genMET_phi", met1_genMET_phi);
  stream.select("patMET_patMETsPF.genMET_pt", met1_genMET_pt);
  stream.select("patMET_patMETsPF.mEtSig", met1_mEtSig);
  stream.select("patMET_patMETsPF.phi", met1_phi);
  stream.select("patMET_patMETsPF.pt", met1_pt);
  stream.select("patMET_patMETsPF.significance", met1_significance);
  stream.select("patMET_patMETsPF.sumEt", met1_sumEt);
  stream.select("patMET_patPFMETsTypeIcorrected.energy", met2_energy);
  stream.select("patMET_patPFMETsTypeIcorrected.et", met2_et);
  stream.select("patMET_patPFMETsTypeIcorrected.genMET_energy", met2_genMET_energy);
  stream.select("patMET_patPFMETsTypeIcorrected.genMET_et", met2_genMET_et);
  stream.select("patMET_patPFMETsTypeIcorrected.genMET_phi", met2_genMET_phi);
  stream.select("patMET_patPFMETsTypeIcorrected.genMET_pt", met2_genMET_pt);
  stream.select("patMET_patPFMETsTypeIcorrected.mEtSig", met2_mEtSig);
  stream.select("patMET_patPFMETsTypeIcorrected.phi", met2_phi);
  stream.select("patMET_patPFMETsTypeIcorrected.pt", met2_pt);
  stream.select("patMET_patPFMETsTypeIcorrected.significance", met2_significance);
  stream.select("patMET_patPFMETsTypeIcorrected.sumEt", met2_sumEt);
  stream.select("patMET_patPFMETsTypeIcorrectedPF.energy", met3_energy);
  stream.select("patMET_patPFMETsTypeIcorrectedPF.et", met3_et);
  stream.select("patMET_patPFMETsTypeIcorrectedPF.genMET_energy", met3_genMET_energy);
  stream.select("patMET_patPFMETsTypeIcorrectedPF.genMET_et", met3_genMET_et);
  stream.select("patMET_patPFMETsTypeIcorrectedPF.genMET_phi", met3_genMET_phi);
  stream.select("patMET_patPFMETsTypeIcorrectedPF.genMET_pt", met3_genMET_pt);
  stream.select("patMET_patPFMETsTypeIcorrectedPF.mEtSig", met3_mEtSig);
  stream.select("patMET_patPFMETsTypeIcorrectedPF.phi", met3_phi);
  stream.select("patMET_patPFMETsTypeIcorrectedPF.pt", met3_pt);
  stream.select("patMET_patPFMETsTypeIcorrectedPF.significance", met3_significance);
  stream.select("patMET_patPFMETsTypeIcorrectedPF.sumEt", met3_sumEt);
  stream.select("patMET_patPFMETsTypeIcorrectedPFLOW.energy", met4_energy);
  stream.select("patMET_patPFMETsTypeIcorrectedPFLOW.et", met4_et);
  stream.select("patMET_patPFMETsTypeIcorrectedPFLOW.genMET_energy", met4_genMET_energy);
  stream.select("patMET_patPFMETsTypeIcorrectedPFLOW.genMET_et", met4_genMET_et);
  stream.select("patMET_patPFMETsTypeIcorrectedPFLOW.genMET_phi", met4_genMET_phi);
  stream.select("patMET_patPFMETsTypeIcorrectedPFLOW.genMET_pt", met4_genMET_pt);
  stream.select("patMET_patPFMETsTypeIcorrectedPFLOW.mEtSig", met4_mEtSig);
  stream.select("patMET_patPFMETsTypeIcorrectedPFLOW.phi", met4_phi);
  stream.select("patMET_patPFMETsTypeIcorrectedPFLOW.pt", met4_pt);
  stream.select("patMET_patPFMETsTypeIcorrectedPFLOW.significance", met4_significance);
  stream.select("patMET_patPFMETsTypeIcorrectedPFLOW.sumEt", met4_sumEt);
  stream.select("patMET_patPFMETsTypeIpIIcorrected.energy", met5_energy);
  stream.select("patMET_patPFMETsTypeIpIIcorrected.et", met5_et);
  stream.select("patMET_patPFMETsTypeIpIIcorrected.genMET_energy", met5_genMET_energy);
  stream.select("patMET_patPFMETsTypeIpIIcorrected.genMET_et", met5_genMET_et);
  stream.select("patMET_patPFMETsTypeIpIIcorrected.genMET_phi", met5_genMET_phi);
  stream.select("patMET_patPFMETsTypeIpIIcorrected.genMET_pt", met5_genMET_pt);
  stream.select("patMET_patPFMETsTypeIpIIcorrected.mEtSig", met5_mEtSig);
  stream.select("patMET_patPFMETsTypeIpIIcorrected.phi", met5_phi);
  stream.select("patMET_patPFMETsTypeIpIIcorrected.pt", met5_pt);
  stream.select("patMET_patPFMETsTypeIpIIcorrected.significance", met5_significance);
  stream.select("patMET_patPFMETsTypeIpIIcorrected.sumEt", met5_sumEt);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPF.energy", met6_energy);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPF.et", met6_et);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPF.genMET_energy", met6_genMET_energy);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPF.genMET_et", met6_genMET_et);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPF.genMET_phi", met6_genMET_phi);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPF.genMET_pt", met6_genMET_pt);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPF.mEtSig", met6_mEtSig);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPF.phi", met6_phi);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPF.pt", met6_pt);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPF.significance", met6_significance);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPF.sumEt", met6_sumEt);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPFLOW.energy", met7_energy);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPFLOW.et", met7_et);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPFLOW.genMET_energy", met7_genMET_energy);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPFLOW.genMET_et", met7_genMET_et);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPFLOW.genMET_phi", met7_genMET_phi);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPFLOW.genMET_pt", met7_genMET_pt);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPFLOW.mEtSig", met7_mEtSig);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPFLOW.phi", met7_phi);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPFLOW.pt", met7_pt);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPFLOW.significance", met7_significance);
  stream.select("patMET_patPFMETsTypeIpIIcorrectedPFLOW.sumEt", met7_sumEt);
  stream.select("patMET_patMETsAK5Calo.energy", met_energy);
  stream.select("patMET_patMETsAK5Calo.et", met_et);
  stream.select("patMET_patMETsAK5Calo.genMET_energy", met_genMET_energy);
  stream.select("patMET_patMETsAK5Calo.genMET_et", met_genMET_et);
  stream.select("patMET_patMETsAK5Calo.genMET_phi", met_genMET_phi);
  stream.select("patMET_patMETsAK5Calo.genMET_pt", met_genMET_pt);
  stream.select("patMET_patMETsAK5Calo.mEtSig", met_mEtSig);
  stream.select("patMET_patMETsAK5Calo.metSignificance", met_metSignificance);
  stream.select("patMET_patMETsAK5Calo.phi", met_phi);
  stream.select("patMET_patMETsAK5Calo.pt", met_pt);
  stream.select("patMET_patMETsAK5Calo.significance", met_significance);
  stream.select("patMET_patMETsAK5Calo.sumEt", met_sumEt);
  stream.select("patMETHelper_patMETsPFLOW.significance_dxx", methelper1_significance_dxx);
  stream.select("patMETHelper_patMETsPFLOW.significance_dxy", methelper1_significance_dxy);
  stream.select("patMETHelper_patMETsPFLOW.significance_dyx", methelper1_significance_dyx);
  stream.select("patMETHelper_patMETsPFLOW.significance_dyy", methelper1_significance_dyy);
  stream.select("patMETHelper_patMETsPF.significance_dxx", methelper_significance_dxx);
  stream.select("patMETHelper_patMETsPF.significance_dxy", methelper_significance_dxy);
  stream.select("patMETHelper_patMETsPF.significance_dyx", methelper_significance_dyx);
  stream.select("patMETHelper_patMETsPF.significance_dyy", methelper_significance_dyy);
  stream.select("patMuon_selectedPatMuonsPF.AllGlobalMuons", muon1_AllGlobalMuons);
  stream.select("patMuon_selectedPatMuonsPF.AllTrackerMuons", muon1_AllTrackerMuons);
  stream.select("patMuon_selectedPatMuonsPF.GlobalMuonPromptTight", muon1_GlobalMuonPromptTight);
  stream.select("patMuon_selectedPatMuonsPF.charge", muon1_charge);
  stream.select("patMuon_selectedPatMuonsPF.chargedHadronIso", muon1_chargedHadronIso);
  stream.select("patMuon_selectedPatMuonsPF.combinedMuon_chi2", muon1_combinedMuon_chi2);
  stream.select("patMuon_selectedPatMuonsPF.combinedMuon_ndof", muon1_combinedMuon_ndof);
  stream.select("patMuon_selectedPatMuonsPF.combinedMuon_numberOfValidHits", muon1_combinedMuon_numberOfValidHits);
  stream.select("patMuon_selectedPatMuonsPF.dB", muon1_dB);
  stream.select("patMuon_selectedPatMuonsPF.ecalIso", muon1_ecalIso);
  stream.select("patMuon_selectedPatMuonsPF.energy", muon1_energy);
  stream.select("patMuon_selectedPatMuonsPF.et", muon1_et);
  stream.select("patMuon_selectedPatMuonsPF.eta", muon1_eta);
  stream.select("patMuon_selectedPatMuonsPF.genLepton_eta", muon1_genLepton_eta);
  stream.select("patMuon_selectedPatMuonsPF.genLepton_pdgId", muon1_genLepton_pdgId);
  stream.select("patMuon_selectedPatMuonsPF.genLepton_phi", muon1_genLepton_phi);
  stream.select("patMuon_selectedPatMuonsPF.genLepton_pt", muon1_genLepton_pt);
  stream.select("patMuon_selectedPatMuonsPF.globalTrack_chi2", muon1_globalTrack_chi2);
  stream.select("patMuon_selectedPatMuonsPF.globalTrack_hitPattern_numberOfValidTrackerHits", muon1_globalTrack_hitPattern_numberOfValidTrackerHits);
  stream.select("patMuon_selectedPatMuonsPF.globalTrack_ndof", muon1_globalTrack_ndof);
  stream.select("patMuon_selectedPatMuonsPF.globalTrack_pt", muon1_globalTrack_pt);
  stream.select("patMuon_selectedPatMuonsPF.globalTrack_ptError", muon1_globalTrack_ptError);
  stream.select("patMuon_selectedPatMuonsPF.hcalIso", muon1_hcalIso);
  stream.select("patMuon_selectedPatMuonsPF.innerTrack_d0", muon1_innerTrack_d0);
  stream.select("patMuon_selectedPatMuonsPF.innerTrack_hitPattern_pixelLayersWithMeasurement", muon1_innerTrack_hitPattern_pixelLayersWithMeasurement);
  stream.select("patMuon_selectedPatMuonsPF.innerTrack_numberOfValidHits", muon1_innerTrack_numberOfValidHits);
  stream.select("patMuon_selectedPatMuonsPF.innerTrack_phi", muon1_innerTrack_phi);
  stream.select("patMuon_selectedPatMuonsPF.innerTrack_pt", muon1_innerTrack_pt);
  stream.select("patMuon_selectedPatMuonsPF.innerTrack_vertex_z", muon1_innerTrack_vertex_z);
  stream.select("patMuon_selectedPatMuonsPF.isGlobalMuon", muon1_isGlobalMuon);
  stream.select("patMuon_selectedPatMuonsPF.isTrackerMuon", muon1_isTrackerMuon);
  stream.select("patMuon_selectedPatMuonsPF.neutralHadronIso", muon1_neutralHadronIso);
  stream.select("patMuon_selectedPatMuonsPF.numberOfMatches", muon1_numberOfMatches);
  stream.select("patMuon_selectedPatMuonsPF.phi", muon1_phi);
  stream.select("patMuon_selectedPatMuonsPF.photonIso", muon1_photonIso);
  stream.select("patMuon_selectedPatMuonsPF.pt", muon1_pt);
  stream.select("patMuon_selectedPatMuonsPF.trackIso", muon1_trackIso);
  stream.select("patMuon_selectedPatMuonsPF.track_d0", muon1_track_d0);
  stream.select("patMuon_selectedPatMuonsPF.track_hitPattern_numberOfValidMuonHits", muon1_track_hitPattern_numberOfValidMuonHits);
  stream.select("patMuon_selectedPatMuonsPF.track_hitPattern_numberOfValidPixelHits", muon1_track_hitPattern_numberOfValidPixelHits);
  stream.select("patMuon_selectedPatMuonsPF.track_normalizedChi2", muon1_track_normalizedChi2);
  stream.select("patMuon_selectedPatMuonsPF.vz", muon1_vz);
  stream.select("patMuon_cleanPatMuons.AllGlobalMuons", muon_AllGlobalMuons);
  stream.select("patMuon_cleanPatMuons.AllTrackerMuons", muon_AllTrackerMuons);
  stream.select("patMuon_cleanPatMuons.GlobalMuonPromptTight", muon_GlobalMuonPromptTight);
  stream.select("patMuon_cleanPatMuons.charge", muon_charge);
  stream.select("patMuon_cleanPatMuons.chargedHadronIso", muon_chargedHadronIso);
  stream.select("patMuon_cleanPatMuons.combinedMuon_chi2", muon_combinedMuon_chi2);
  stream.select("patMuon_cleanPatMuons.combinedMuon_ndof", muon_combinedMuon_ndof);
  stream.select("patMuon_cleanPatMuons.combinedMuon_numberOfValidHits", muon_combinedMuon_numberOfValidHits);
  stream.select("patMuon_cleanPatMuons.dB", muon_dB);
  stream.select("patMuon_cleanPatMuons.ecalIso", muon_ecalIso);
  stream.select("patMuon_cleanPatMuons.energy", muon_energy);
  stream.select("patMuon_cleanPatMuons.et", muon_et);
  stream.select("patMuon_cleanPatMuons.eta", muon_eta);
  stream.select("patMuon_cleanPatMuons.genLepton_eta", muon_genLepton_eta);
  stream.select("patMuon_cleanPatMuons.genLepton_pdgId", muon_genLepton_pdgId);
  stream.select("patMuon_cleanPatMuons.genLepton_phi", muon_genLepton_phi);
  stream.select("patMuon_cleanPatMuons.genLepton_pt", muon_genLepton_pt);
  stream.select("patMuon_cleanPatMuons.globalTrack_chi2", muon_globalTrack_chi2);
  stream.select("patMuon_cleanPatMuons.globalTrack_hitPattern_numberOfValidTrackerHits", muon_globalTrack_hitPattern_numberOfValidTrackerHits);
  stream.select("patMuon_cleanPatMuons.globalTrack_ndof", muon_globalTrack_ndof);
  stream.select("patMuon_cleanPatMuons.globalTrack_pt", muon_globalTrack_pt);
  stream.select("patMuon_cleanPatMuons.globalTrack_ptError", muon_globalTrack_ptError);
  stream.select("patMuon_cleanPatMuons.hcalIso", muon_hcalIso);
  stream.select("patMuon_cleanPatMuons.innerTrack_d0", muon_innerTrack_d0);
  stream.select("patMuon_cleanPatMuons.innerTrack_hitPattern_pixelLayersWithMeasurement", muon_innerTrack_hitPattern_pixelLayersWithMeasurement);
  stream.select("patMuon_cleanPatMuons.innerTrack_numberOfValidHits", muon_innerTrack_numberOfValidHits);
  stream.select("patMuon_cleanPatMuons.innerTrack_phi", muon_innerTrack_phi);
  stream.select("patMuon_cleanPatMuons.innerTrack_pt", muon_innerTrack_pt);
  stream.select("patMuon_cleanPatMuons.innerTrack_vertex_z", muon_innerTrack_vertex_z);
  stream.select("patMuon_cleanPatMuons.isGlobalMuon", muon_isGlobalMuon);
  stream.select("patMuon_cleanPatMuons.isTrackerMuon", muon_isTrackerMuon);
  stream.select("patMuon_cleanPatMuons.neutralHadronIso", muon_neutralHadronIso);
  stream.select("patMuon_cleanPatMuons.numberOfMatches", muon_numberOfMatches);
  stream.select("patMuon_cleanPatMuons.phi", muon_phi);
  stream.select("patMuon_cleanPatMuons.photonIso", muon_photonIso);
  stream.select("patMuon_cleanPatMuons.pt", muon_pt);
  stream.select("patMuon_cleanPatMuons.trackIso", muon_trackIso);
  stream.select("patMuon_cleanPatMuons.track_d0", muon_track_d0);
  stream.select("patMuon_cleanPatMuons.track_hitPattern_numberOfValidMuonHits", muon_track_hitPattern_numberOfValidMuonHits);
  stream.select("patMuon_cleanPatMuons.track_hitPattern_numberOfValidPixelHits", muon_track_hitPattern_numberOfValidPixelHits);
  stream.select("patMuon_cleanPatMuons.track_normalizedChi2", muon_track_normalizedChi2);
  stream.select("patMuon_cleanPatMuons.vz", muon_vz);
  stream.select("patMuonHelper_selectedPatMuonsPF.dxywrtBeamSpot", muonhelper1_dxywrtBeamSpot);
  stream.select("patMuonHelper_selectedPatMuonsPF.muonPFcandptdiff", muonhelper1_muonPFcandptdiff);
  stream.select("patMuonHelper_cleanPatMuons.dxywrtBeamSpot", muonhelper_dxywrtBeamSpot);
  stream.select("patMuonHelper_cleanPatMuons.muonPFcandptdiff", muonhelper_muonPFcandptdiff);
  stream.select("nGenEventInfoProductHelper_generator_NNPDF", nGenEventInfoProductHelper_generator_NNPDF);
  stream.select("nGenEventInfoProductHelper_generator_cteq", nGenEventInfoProductHelper_generator_cteq);
  stream.select("nGenEventInfoProductHelper_generator_mstw", nGenEventInfoProductHelper_generator_mstw);
  stream.select("nPileupSummaryInfo_addPileupInfo", nPileupSummaryInfo_addPileupInfo);
  stream.select("npatElectron_cleanPatElectrons", nelectron);
  stream.select("npatElectron_selectedPatElectronsPF", nelectron1);
  stream.select("npatElectronHelper_cleanPatElectrons", nelectronhelper);
  stream.select("npatElectronHelper_selectedPatElectronsPF", nelectronhelper1);
  stream.select("nrecoGenParticleHelperRA2b_genParticles", ngenparticlehelperra2b);
  stream.select("npatJet_cleanPatJetsAK5PF", njet);
  stream.select("npatJet_selectedPatJetsPF", njet1);
  stream.select("npatJet_selectedPatJetsPFLOW", njet2);
  stream.select("npatJetHelper_cleanPatJetsAK5PF", njethelper);
  stream.select("npatJetHelper_selectedPatJetsPF", njethelper1);
  stream.select("npatJetHelper_selectedPatJetsPFLOW", njethelper2);
  stream.select("npatMET_patMETsAK5Calo", nmet);
  stream.select("npatMET_patMETsPF", nmet1);
  stream.select("npatMET_patPFMETsTypeIcorrected", nmet2);
  stream.select("npatMET_patPFMETsTypeIcorrectedPF", nmet3);
  stream.select("npatMET_patPFMETsTypeIcorrectedPFLOW", nmet4);
  stream.select("npatMET_patPFMETsTypeIpIIcorrected", nmet5);
  stream.select("npatMET_patPFMETsTypeIpIIcorrectedPF", nmet6);
  stream.select("npatMET_patPFMETsTypeIpIIcorrectedPFLOW", nmet7);
  stream.select("npatMETHelper_patMETsPF", nmethelper);
  stream.select("npatMETHelper_patMETsPFLOW", nmethelper1);
  stream.select("npatMuon_cleanPatMuons", nmuon);
  stream.select("npatMuon_selectedPatMuonsPF", nmuon1);
  stream.select("npatMuonHelper_cleanPatMuons", nmuonhelper);
  stream.select("npatMuonHelper_selectedPatMuonsPF", nmuonhelper1);
  stream.select("npatTau_selectedPatTausPF", ntau);
  stream.select("npatTauHelper_selectedPatTausPF", ntauhelper);
  stream.select("nrecoVertex_offlinePrimaryVertices", nvertex);
  stream.select("PileupSummaryInfo_addPileupInfo.getBunchCrossing", pileupsummaryinfo_addpileupinfo_getBunchCrossing);
  stream.select("PileupSummaryInfo_addPileupInfo.getPU_NumInteractions", pileupsummaryinfo_addpileupinfo_getPU_NumInteractions);
  stream.select("sdouble_kt6PFJets_rho.value", sdouble_kt6pfjets_rho_value);
  stream.select("sdouble_kt6PFJets_sigma.value", sdouble_kt6pfjets_sigma_value);
  stream.select("suint_flavorHistoryFilter.value", suint_flavorhistoryfilter_value);
  stream.select("patTau_selectedPatTausPF.caloIso", tau_caloIso);
  stream.select("patTau_selectedPatTausPF.ecalIso", tau_ecalIso);
  stream.select("patTau_selectedPatTausPF.energy", tau_energy);
  stream.select("patTau_selectedPatTausPF.et", tau_et);
  stream.select("patTau_selectedPatTausPF.eta", tau_eta);
  stream.select("patTau_selectedPatTausPF.hcalIso", tau_hcalIso);
  stream.select("patTau_selectedPatTausPF.phi", tau_phi);
  stream.select("patTau_selectedPatTausPF.pt", tau_pt);
  stream.select("patTau_selectedPatTausPF.tauID_againstElectron", tau_tauID_againstElectron);
  stream.select("patTau_selectedPatTausPF.tauID_againstMuon", tau_tauID_againstMuon);
  stream.select("patTau_selectedPatTausPF.tauID_byIsolation", tau_tauID_byIsolation);
  stream.select("patTau_selectedPatTausPF.tauID_byTaNC", tau_tauID_byTaNC);
  stream.select("patTau_selectedPatTausPF.tauID_byTaNCfrHalfPercent", tau_tauID_byTaNCfrHalfPercent);
  stream.select("patTau_selectedPatTausPF.tauID_byTaNCfrQuarterPercent", tau_tauID_byTaNCfrQuarterPercent);
  stream.select("patTau_selectedPatTausPF.trackIso", tau_trackIso);
  stream.select("patTauHelper_selectedPatTausPF.genTauDecayModeID", tauhelper_genTauDecayModeID);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.csctighthaloFilter", triggerresultshelper1_csctighthaloFilter);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.eenoiseFilter", triggerresultshelper1_eenoiseFilter);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.greedymuonFilter", triggerresultshelper1_greedymuonFilter);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.hbhenoiseFilter", triggerresultshelper1_hbhenoiseFilter);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.inconsistentmuonFilter", triggerresultshelper1_inconsistentmuonFilter);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.ra2ecaltpFilter", triggerresultshelper1_ra2ecaltpFilter);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.recovrechitFilter", triggerresultshelper1_recovrechitFilter);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.scrapingvetoFilter", triggerresultshelper1_scrapingvetoFilter);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.trackingfailureFilter", triggerresultshelper1_trackingfailureFilter);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.trackingfailureFilterPFLOW", triggerresultshelper1_trackingfailureFilterPFLOW);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v1", triggerresultshelper_HLT_HT150_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v10", triggerresultshelper_HLT_HT150_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v10_prs", triggerresultshelper_HLT_HT150_v10_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v1_prs", triggerresultshelper_HLT_HT150_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v2", triggerresultshelper_HLT_HT150_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v2_prs", triggerresultshelper_HLT_HT150_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v3", triggerresultshelper_HLT_HT150_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v3_prs", triggerresultshelper_HLT_HT150_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v4", triggerresultshelper_HLT_HT150_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v4_prs", triggerresultshelper_HLT_HT150_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v5", triggerresultshelper_HLT_HT150_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v5_prs", triggerresultshelper_HLT_HT150_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v6", triggerresultshelper_HLT_HT150_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v6_prs", triggerresultshelper_HLT_HT150_v6_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v7", triggerresultshelper_HLT_HT150_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v7_prs", triggerresultshelper_HLT_HT150_v7_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v8", triggerresultshelper_HLT_HT150_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v8_prs", triggerresultshelper_HLT_HT150_v8_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v9", triggerresultshelper_HLT_HT150_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT150_v9_prs", triggerresultshelper_HLT_HT150_v9_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT2000_v1", triggerresultshelper_HLT_HT2000_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT2000_v1_prs", triggerresultshelper_HLT_HT2000_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT2000_v2", triggerresultshelper_HLT_HT2000_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT2000_v2_prs", triggerresultshelper_HLT_HT2000_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT2000_v3", triggerresultshelper_HLT_HT2000_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT2000_v3_prs", triggerresultshelper_HLT_HT2000_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT2000_v4", triggerresultshelper_HLT_HT2000_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT2000_v4_prs", triggerresultshelper_HLT_HT2000_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v1", triggerresultshelper_HLT_HT200_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v10", triggerresultshelper_HLT_HT200_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v10_prs", triggerresultshelper_HLT_HT200_v10_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v1_prs", triggerresultshelper_HLT_HT200_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v2", triggerresultshelper_HLT_HT200_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v2_prs", triggerresultshelper_HLT_HT200_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v3", triggerresultshelper_HLT_HT200_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v3_prs", triggerresultshelper_HLT_HT200_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v4", triggerresultshelper_HLT_HT200_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v4_prs", triggerresultshelper_HLT_HT200_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v5", triggerresultshelper_HLT_HT200_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v5_prs", triggerresultshelper_HLT_HT200_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v6", triggerresultshelper_HLT_HT200_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v6_prs", triggerresultshelper_HLT_HT200_v6_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v7", triggerresultshelper_HLT_HT200_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v7_prs", triggerresultshelper_HLT_HT200_v7_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v8", triggerresultshelper_HLT_HT200_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v8_prs", triggerresultshelper_HLT_HT200_v8_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v9", triggerresultshelper_HLT_HT200_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT200_v9_prs", triggerresultshelper_HLT_HT200_v9_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT100_v2", triggerresultshelper_HLT_HT250_MHT100_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT100_v2_prs", triggerresultshelper_HLT_HT250_MHT100_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT60_v2", triggerresultshelper_HLT_HT250_MHT60_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT60_v2_prs", triggerresultshelper_HLT_HT250_MHT60_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT60_v3", triggerresultshelper_HLT_HT250_MHT60_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT60_v3_prs", triggerresultshelper_HLT_HT250_MHT60_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT70_v1", triggerresultshelper_HLT_HT250_MHT70_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT70_v1_prs", triggerresultshelper_HLT_HT250_MHT70_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT80_v3", triggerresultshelper_HLT_HT250_MHT80_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT80_v3_prs", triggerresultshelper_HLT_HT250_MHT80_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT80_v4", triggerresultshelper_HLT_HT250_MHT80_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT80_v4_prs", triggerresultshelper_HLT_HT250_MHT80_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT90_v1", triggerresultshelper_HLT_HT250_MHT90_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT90_v1_prs", triggerresultshelper_HLT_HT250_MHT90_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT90_v2", triggerresultshelper_HLT_HT250_MHT90_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_MHT90_v2_prs", triggerresultshelper_HLT_HT250_MHT90_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v1", triggerresultshelper_HLT_HT250_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v10", triggerresultshelper_HLT_HT250_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v10_prs", triggerresultshelper_HLT_HT250_v10_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v1_prs", triggerresultshelper_HLT_HT250_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v2", triggerresultshelper_HLT_HT250_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v2_prs", triggerresultshelper_HLT_HT250_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v3", triggerresultshelper_HLT_HT250_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v3_prs", triggerresultshelper_HLT_HT250_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v4", triggerresultshelper_HLT_HT250_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v4_prs", triggerresultshelper_HLT_HT250_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v5", triggerresultshelper_HLT_HT250_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v5_prs", triggerresultshelper_HLT_HT250_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v6", triggerresultshelper_HLT_HT250_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v6_prs", triggerresultshelper_HLT_HT250_v6_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v7", triggerresultshelper_HLT_HT250_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v7_prs", triggerresultshelper_HLT_HT250_v7_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v8", triggerresultshelper_HLT_HT250_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v8_prs", triggerresultshelper_HLT_HT250_v8_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v9", triggerresultshelper_HLT_HT250_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT250_v9_prs", triggerresultshelper_HLT_HT250_v9_prs);
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
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v7", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v7_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v7_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v8", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT55_v8_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT55_v8_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT65_v1", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT65_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT65_v1_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT65_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v7", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_PFMHT75_v7_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_PFMHT75_v7_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v2", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v2_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v3", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v3_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v4", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v4_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v5", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v5_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v6", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v6_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v6_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v7", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_CentralJet30_BTagIP_v7_prs", triggerresultshelper_HLT_HT300_CentralJet30_BTagIP_v7_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v2", triggerresultshelper_HLT_HT300_MHT75_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v2_prs", triggerresultshelper_HLT_HT300_MHT75_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v3", triggerresultshelper_HLT_HT300_MHT75_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v3_prs", triggerresultshelper_HLT_HT300_MHT75_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v4", triggerresultshelper_HLT_HT300_MHT75_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v4_prs", triggerresultshelper_HLT_HT300_MHT75_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v7", triggerresultshelper_HLT_HT300_MHT75_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v7_prs", triggerresultshelper_HLT_HT300_MHT75_v7_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v8", triggerresultshelper_HLT_HT300_MHT75_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT75_v8_prs", triggerresultshelper_HLT_HT300_MHT75_v8_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT80_v1", triggerresultshelper_HLT_HT300_MHT80_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT80_v1_prs", triggerresultshelper_HLT_HT300_MHT80_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT80_v2", triggerresultshelper_HLT_HT300_MHT80_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT80_v2_prs", triggerresultshelper_HLT_HT300_MHT80_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT90_v2", triggerresultshelper_HLT_HT300_MHT90_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_MHT90_v2_prs", triggerresultshelper_HLT_HT300_MHT90_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v2", triggerresultshelper_HLT_HT300_PFMHT55_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v2_prs", triggerresultshelper_HLT_HT300_PFMHT55_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v3", triggerresultshelper_HLT_HT300_PFMHT55_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v3_prs", triggerresultshelper_HLT_HT300_PFMHT55_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v4", triggerresultshelper_HLT_HT300_PFMHT55_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v4_prs", triggerresultshelper_HLT_HT300_PFMHT55_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v5", triggerresultshelper_HLT_HT300_PFMHT55_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v5_prs", triggerresultshelper_HLT_HT300_PFMHT55_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v7", triggerresultshelper_HLT_HT300_PFMHT55_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v7_prs", triggerresultshelper_HLT_HT300_PFMHT55_v7_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v8", triggerresultshelper_HLT_HT300_PFMHT55_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_PFMHT55_v8_prs", triggerresultshelper_HLT_HT300_PFMHT55_v8_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v1", triggerresultshelper_HLT_HT300_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v10", triggerresultshelper_HLT_HT300_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v10_prs", triggerresultshelper_HLT_HT300_v10_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v11", triggerresultshelper_HLT_HT300_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v11_prs", triggerresultshelper_HLT_HT300_v11_prs);
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
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v8", triggerresultshelper_HLT_HT300_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v8_prs", triggerresultshelper_HLT_HT300_v8_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v9", triggerresultshelper_HLT_HT300_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT300_v9_prs", triggerresultshelper_HLT_HT300_v9_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_MHT70_v1", triggerresultshelper_HLT_HT350_MHT70_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_MHT70_v1_prs", triggerresultshelper_HLT_HT350_MHT70_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_MHT70_v2", triggerresultshelper_HLT_HT350_MHT70_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_MHT70_v2_prs", triggerresultshelper_HLT_HT350_MHT70_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_MHT80_v2", triggerresultshelper_HLT_HT350_MHT80_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_MHT80_v2_prs", triggerresultshelper_HLT_HT350_MHT80_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_MHT90_v1", triggerresultshelper_HLT_HT350_MHT90_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_MHT90_v1_prs", triggerresultshelper_HLT_HT350_MHT90_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v1", triggerresultshelper_HLT_HT350_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v10", triggerresultshelper_HLT_HT350_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v10_prs", triggerresultshelper_HLT_HT350_v10_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v1_prs", triggerresultshelper_HLT_HT350_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v2", triggerresultshelper_HLT_HT350_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v2_prs", triggerresultshelper_HLT_HT350_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v3", triggerresultshelper_HLT_HT350_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v3_prs", triggerresultshelper_HLT_HT350_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v4", triggerresultshelper_HLT_HT350_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v4_prs", triggerresultshelper_HLT_HT350_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v5", triggerresultshelper_HLT_HT350_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v5_prs", triggerresultshelper_HLT_HT350_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v6", triggerresultshelper_HLT_HT350_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v6_prs", triggerresultshelper_HLT_HT350_v6_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v7", triggerresultshelper_HLT_HT350_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v7_prs", triggerresultshelper_HLT_HT350_v7_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v8", triggerresultshelper_HLT_HT350_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v8_prs", triggerresultshelper_HLT_HT350_v8_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v9", triggerresultshelper_HLT_HT350_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT350_v9_prs", triggerresultshelper_HLT_HT350_v9_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT360_v2", triggerresultshelper_HLT_HT360_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT360_v2_prs", triggerresultshelper_HLT_HT360_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_MHT80_v1", triggerresultshelper_HLT_HT400_MHT80_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_MHT80_v1_prs", triggerresultshelper_HLT_HT400_MHT80_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v1", triggerresultshelper_HLT_HT400_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v10", triggerresultshelper_HLT_HT400_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v10_prs", triggerresultshelper_HLT_HT400_v10_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v1_prs", triggerresultshelper_HLT_HT400_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v2", triggerresultshelper_HLT_HT400_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v2_prs", triggerresultshelper_HLT_HT400_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v3", triggerresultshelper_HLT_HT400_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v3_prs", triggerresultshelper_HLT_HT400_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v4", triggerresultshelper_HLT_HT400_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v4_prs", triggerresultshelper_HLT_HT400_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v5", triggerresultshelper_HLT_HT400_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v5_prs", triggerresultshelper_HLT_HT400_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v6", triggerresultshelper_HLT_HT400_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v6_prs", triggerresultshelper_HLT_HT400_v6_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v7", triggerresultshelper_HLT_HT400_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v7_prs", triggerresultshelper_HLT_HT400_v7_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v8", triggerresultshelper_HLT_HT400_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v8_prs", triggerresultshelper_HLT_HT400_v8_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v9", triggerresultshelper_HLT_HT400_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT400_v9_prs", triggerresultshelper_HLT_HT400_v9_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v1", triggerresultshelper_HLT_HT450_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v10", triggerresultshelper_HLT_HT450_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v10_prs", triggerresultshelper_HLT_HT450_v10_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v1_prs", triggerresultshelper_HLT_HT450_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v2", triggerresultshelper_HLT_HT450_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v2_prs", triggerresultshelper_HLT_HT450_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v3", triggerresultshelper_HLT_HT450_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v3_prs", triggerresultshelper_HLT_HT450_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v4", triggerresultshelper_HLT_HT450_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v4_prs", triggerresultshelper_HLT_HT450_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v5", triggerresultshelper_HLT_HT450_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v5_prs", triggerresultshelper_HLT_HT450_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v6", triggerresultshelper_HLT_HT450_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v6_prs", triggerresultshelper_HLT_HT450_v6_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v7", triggerresultshelper_HLT_HT450_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v7_prs", triggerresultshelper_HLT_HT450_v7_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v8", triggerresultshelper_HLT_HT450_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v8_prs", triggerresultshelper_HLT_HT450_v8_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v9", triggerresultshelper_HLT_HT450_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v9_prs", triggerresultshelper_HLT_HT450_v9_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v10", triggerresultshelper_HLT_HT500_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v10_prs", triggerresultshelper_HLT_HT500_v10_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v3", triggerresultshelper_HLT_HT500_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v3_prs", triggerresultshelper_HLT_HT500_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v4", triggerresultshelper_HLT_HT500_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v4_prs", triggerresultshelper_HLT_HT500_v4_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v5", triggerresultshelper_HLT_HT500_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v5_prs", triggerresultshelper_HLT_HT500_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v6", triggerresultshelper_HLT_HT500_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v6_prs", triggerresultshelper_HLT_HT500_v6_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v7", triggerresultshelper_HLT_HT500_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v7_prs", triggerresultshelper_HLT_HT500_v7_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v8", triggerresultshelper_HLT_HT500_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v8_prs", triggerresultshelper_HLT_HT500_v8_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v9", triggerresultshelper_HLT_HT500_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v9_prs", triggerresultshelper_HLT_HT500_v9_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v5", triggerresultshelper_HLT_HT550_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v5_prs", triggerresultshelper_HLT_HT550_v5_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v6", triggerresultshelper_HLT_HT550_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v6_prs", triggerresultshelper_HLT_HT550_v6_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v7", triggerresultshelper_HLT_HT550_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v7_prs", triggerresultshelper_HLT_HT550_v7_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v8", triggerresultshelper_HLT_HT550_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v8_prs", triggerresultshelper_HLT_HT550_v8_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT600_v1", triggerresultshelper_HLT_HT600_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT600_v1_prs", triggerresultshelper_HLT_HT600_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT600_v2", triggerresultshelper_HLT_HT600_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT600_v2_prs", triggerresultshelper_HLT_HT600_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT600_v3", triggerresultshelper_HLT_HT600_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT600_v3_prs", triggerresultshelper_HLT_HT600_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT650_v1", triggerresultshelper_HLT_HT650_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT650_v1_prs", triggerresultshelper_HLT_HT650_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT650_v2", triggerresultshelper_HLT_HT650_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT650_v2_prs", triggerresultshelper_HLT_HT650_v2_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT650_v3", triggerresultshelper_HLT_HT650_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT650_v3_prs", triggerresultshelper_HLT_HT650_v3_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT750_v1", triggerresultshelper_HLT_HT750_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT750_v1_prs", triggerresultshelper_HLT_HT750_v1_prs);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT750_v2", triggerresultshelper_HLT_HT750_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT750_v2_prs", triggerresultshelper_HLT_HT750_v2_prs);
  stream.select("triggerTriggerEventHelper_hltTriggerSummaryAOD_HLT.btagIPjetEta1", triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetEta1);
  stream.select("triggerTriggerEventHelper_hltTriggerSummaryAOD_HLT.btagIPjetEta2", triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetEta2);
  stream.select("triggerTriggerEventHelper_hltTriggerSummaryAOD_HLT.btagIPjetEta3", triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetEta3);
  stream.select("triggerTriggerEventHelper_hltTriggerSummaryAOD_HLT.btagIPjetEta4", triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetEta4);
  stream.select("triggerTriggerEventHelper_hltTriggerSummaryAOD_HLT.btagIPjetPhi1", triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetPhi1);
  stream.select("triggerTriggerEventHelper_hltTriggerSummaryAOD_HLT.btagIPjetPhi2", triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetPhi2);
  stream.select("triggerTriggerEventHelper_hltTriggerSummaryAOD_HLT.btagIPjetPhi3", triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetPhi3);
  stream.select("triggerTriggerEventHelper_hltTriggerSummaryAOD_HLT.btagIPjetPhi4", triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetPhi4);
  stream.select("triggerTriggerEventHelper_hltTriggerSummaryAOD_HLT.btagIPjetPt1", triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetPt1);
  stream.select("triggerTriggerEventHelper_hltTriggerSummaryAOD_HLT.btagIPjetPt2", triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetPt2);
  stream.select("triggerTriggerEventHelper_hltTriggerSummaryAOD_HLT.btagIPjetPt3", triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetPt3);
  stream.select("triggerTriggerEventHelper_hltTriggerSummaryAOD_HLT.btagIPjetPt4", triggertriggereventhelper_hlttriggersummaryaod_hlt_btagIPjetPt4);
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
  double	caloIso;
  double	charge;
  double	chargedHadronIso;
  double	convDcot;
  double	convDist;
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
  double	sigmaIetaIeta;
  double	simpleEleId80cIso;
  double	simpleEleId80relIso;
  double	simpleEleId95cIso;
  double	simpleEleId95relIso;
  double	superCluster_eta;
  double	trackIso;
  double	vertex_z;
  double	vz;
};
std::vector<electron_s> electron(15);

std::ostream& operator<<(std::ostream& os, const electron_s& o)
{
  char r[1024];
  os << "electron" << std::endl;
  sprintf(r, "  %-32s: %f\n", "caloIso", (double)o.caloIso); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronIso", (double)o.chargedHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "convDcot", (double)o.convDcot); os << r;
  sprintf(r, "  %-32s: %f\n", "convDist", (double)o.convDist); os << r;
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
  sprintf(r, "  %-32s: %f\n", "sigmaIetaIeta", (double)o.sigmaIetaIeta); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleEleId80cIso", (double)o.simpleEleId80cIso); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleEleId80relIso", (double)o.simpleEleId80relIso); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleEleId95cIso", (double)o.simpleEleId95cIso); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleEleId95relIso", (double)o.simpleEleId95relIso); os << r;
  sprintf(r, "  %-32s: %f\n", "superCluster_eta", (double)o.superCluster_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "trackIso", (double)o.trackIso); os << r;
  sprintf(r, "  %-32s: %f\n", "vertex_z", (double)o.vertex_z); os << r;
  sprintf(r, "  %-32s: %f\n", "vz", (double)o.vz); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct electron1_s
{
  double	caloIso;
  double	charge;
  double	chargedHadronIso;
  double	convDcot;
  double	convDist;
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
  double	sigmaIetaIeta;
  double	superCluster_eta;
  double	trackIso;
  double	vertex_z;
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
  sprintf(r, "  %-32s: %f\n", "convDcot", (double)o.convDcot); os << r;
  sprintf(r, "  %-32s: %f\n", "convDist", (double)o.convDist); os << r;
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
  sprintf(r, "  %-32s: %f\n", "sigmaIetaIeta", (double)o.sigmaIetaIeta); os << r;
  sprintf(r, "  %-32s: %f\n", "superCluster_eta", (double)o.superCluster_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "trackIso", (double)o.trackIso); os << r;
  sprintf(r, "  %-32s: %f\n", "vertex_z", (double)o.vertex_z); os << r;
  sprintf(r, "  %-32s: %f\n", "vz", (double)o.vz); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct electronhelper_s
{
  double	dxywrtBeamSpot;
};
std::vector<electronhelper_s> electronhelper(15);

std::ostream& operator<<(std::ostream& os, const electronhelper_s& o)
{
  char r[1024];
  os << "electronhelper" << std::endl;
  sprintf(r, "  %-32s: %f\n", "dxywrtBeamSpot", (double)o.dxywrtBeamSpot); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct electronhelper1_s
{
  double	dxywrtBeamSpot;
};
std::vector<electronhelper1_s> electronhelper1(7);

std::ostream& operator<<(std::ostream& os, const electronhelper1_s& o)
{
  char r[1024];
  os << "electronhelper1" << std::endl;
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
std::vector<geneventinfoproducthelper_s> geneventinfoproducthelper(201);

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
struct geneventinfoproducthelper1_s
{
  double	pdf1;
  double	pdf2;
  double	pdfweight;
  double	pdfweightsum;
};
std::vector<geneventinfoproducthelper1_s> geneventinfoproducthelper1(91);

std::ostream& operator<<(std::ostream& os, const geneventinfoproducthelper1_s& o)
{
  char r[1024];
  os << "geneventinfoproducthelper1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "pdf1", (double)o.pdf1); os << r;
  sprintf(r, "  %-32s: %f\n", "pdf2", (double)o.pdf2); os << r;
  sprintf(r, "  %-32s: %f\n", "pdfweight", (double)o.pdfweight); os << r;
  sprintf(r, "  %-32s: %f\n", "pdfweightsum", (double)o.pdfweightsum); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct geneventinfoproducthelper2_s
{
  double	pdf1;
  double	pdf2;
  double	pdfweight;
  double	pdfweightsum;
};
std::vector<geneventinfoproducthelper2_s> geneventinfoproducthelper2(83);

std::ostream& operator<<(std::ostream& os, const geneventinfoproducthelper2_s& o)
{
  char r[1024];
  os << "geneventinfoproducthelper2" << std::endl;
  sprintf(r, "  %-32s: %f\n", "pdf1", (double)o.pdf1); os << r;
  sprintf(r, "  %-32s: %f\n", "pdf2", (double)o.pdf2); os << r;
  sprintf(r, "  %-32s: %f\n", "pdfweight", (double)o.pdfweight); os << r;
  sprintf(r, "  %-32s: %f\n", "pdfweightsum", (double)o.pdfweightsum); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct genparticlehelperra2_s
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
std::vector<genparticlehelperra2_s> genparticlehelperra2(73);

std::ostream& operator<<(std::ostream& os, const genparticlehelperra2_s& o)
{
  char r[1024];
  os << "genparticlehelperra2" << std::endl;
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
  double	HFEMEnergy;
  double	HFEMMultiplicity;
  double	HFHadronEnergy;
  double	HFHadronMultiplicity;
  double	chargedEmEnergyFraction;
  double	chargedHadronEnergyFraction;
  double	chargedHadronMultiplicity;
  double	chargedMuEnergy;
  double	chargedMultiplicity;
  double	combinedSecondaryVertexBJetTags;
  double	combinedSecondaryVertexMVABJetTags;
  double	electronEnergy;
  double	electronMultiplicity;
  double	energy;
  double	et;
  double	eta;
  double	genJet_energy;
  double	genJet_eta;
  double	genJet_invisibleEnergy;
  double	genJet_phi;
  double	genJet_pt;
  double	genParton_energy;
  double	genParton_eta;
  double	genParton_pdgId;
  double	genParton_phi;
  double	genParton_pt;
  double	jetArea;
  double	jetBProbabilityBJetTags;
  double	jetID_fHPD;
  double	jetID_n90Hits;
  double	jetProbabilityBJetTags;
  double	muonEnergy;
  double	neutralEmEnergyFraction;
  double	neutralHadronEnergy;
  double	neutralHadronEnergyFraction;
  double	neutralHadronMultiplicity;
  double	neutralMultiplicity;
  double	numberOfDaughters;
  double	partonFlavour;
  double	phi;
  double	photonEnergy;
  double	photonMultiplicity;
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
  double	uncor_et;
  double	uncor_eta;
  double	uncor_phi;
  double	uncor_pt;
};
std::vector<jet_s> jet(205);

std::ostream& operator<<(std::ostream& os, const jet_s& o)
{
  char r[1024];
  os << "jet" << std::endl;
  sprintf(r, "  %-32s: %f\n", "HFEMEnergy", (double)o.HFEMEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "HFEMMultiplicity", (double)o.HFEMMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "HFHadronEnergy", (double)o.HFHadronEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "HFHadronMultiplicity", (double)o.HFHadronMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedEmEnergyFraction", (double)o.chargedEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronEnergyFraction", (double)o.chargedHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronMultiplicity", (double)o.chargedHadronMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedMuEnergy", (double)o.chargedMuEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedMultiplicity", (double)o.chargedMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexBJetTags", (double)o.combinedSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexMVABJetTags", (double)o.combinedSecondaryVertexMVABJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "electronEnergy", (double)o.electronEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "electronMultiplicity", (double)o.electronMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genJet_energy", (double)o.genJet_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "genJet_eta", (double)o.genJet_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genJet_invisibleEnergy", (double)o.genJet_invisibleEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "genJet_phi", (double)o.genJet_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genJet_pt", (double)o.genJet_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "genParton_energy", (double)o.genParton_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "genParton_eta", (double)o.genParton_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genParton_pdgId", (double)o.genParton_pdgId); os << r;
  sprintf(r, "  %-32s: %f\n", "genParton_phi", (double)o.genParton_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genParton_pt", (double)o.genParton_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "jetArea", (double)o.jetArea); os << r;
  sprintf(r, "  %-32s: %f\n", "jetBProbabilityBJetTags", (double)o.jetBProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetID_fHPD", (double)o.jetID_fHPD); os << r;
  sprintf(r, "  %-32s: %f\n", "jetID_n90Hits", (double)o.jetID_n90Hits); os << r;
  sprintf(r, "  %-32s: %f\n", "jetProbabilityBJetTags", (double)o.jetProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "muonEnergy", (double)o.muonEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralEmEnergyFraction", (double)o.neutralEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergy", (double)o.neutralHadronEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergyFraction", (double)o.neutralHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronMultiplicity", (double)o.neutralHadronMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralMultiplicity", (double)o.neutralMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfDaughters", (double)o.numberOfDaughters); os << r;
  sprintf(r, "  %-32s: %f\n", "partonFlavour", (double)o.partonFlavour); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "photonEnergy", (double)o.photonEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "photonMultiplicity", (double)o.photonMultiplicity); os << r;
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
  sprintf(r, "  %-32s: %f\n", "uncor_et", (double)o.uncor_et); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_eta", (double)o.uncor_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_phi", (double)o.uncor_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_pt", (double)o.uncor_pt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jet1_s
{
  double	HFEMEnergy;
  double	HFEMMultiplicity;
  double	HFHadronEnergy;
  double	HFHadronMultiplicity;
  double	chargedEmEnergyFraction;
  double	chargedHadronEnergyFraction;
  double	chargedHadronMultiplicity;
  double	chargedMuEnergy;
  double	chargedMultiplicity;
  double	combinedSecondaryVertexBJetTags;
  double	combinedSecondaryVertexMVABJetTags;
  double	electronEnergy;
  double	electronMultiplicity;
  double	energy;
  double	et;
  double	eta;
  double	genJet_energy;
  double	genJet_eta;
  double	genJet_invisibleEnergy;
  double	genJet_phi;
  double	genJet_pt;
  double	genParton_energy;
  double	genParton_eta;
  double	genParton_pdgId;
  double	genParton_phi;
  double	genParton_pt;
  double	jetArea;
  double	jetBProbabilityBJetTags;
  double	jetID_fHPD;
  double	jetID_n90Hits;
  double	jetProbabilityBJetTags;
  double	muonEnergy;
  double	neutralEmEnergyFraction;
  double	neutralHadronEnergy;
  double	neutralHadronEnergyFraction;
  double	neutralHadronMultiplicity;
  double	neutralMultiplicity;
  double	numberOfDaughters;
  double	partonFlavour;
  double	phi;
  double	photonEnergy;
  double	photonMultiplicity;
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
  double	uncor_et;
  double	uncor_eta;
  double	uncor_phi;
  double	uncor_pt;
};
std::vector<jet1_s> jet1(177);

std::ostream& operator<<(std::ostream& os, const jet1_s& o)
{
  char r[1024];
  os << "jet1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "HFEMEnergy", (double)o.HFEMEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "HFEMMultiplicity", (double)o.HFEMMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "HFHadronEnergy", (double)o.HFHadronEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "HFHadronMultiplicity", (double)o.HFHadronMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedEmEnergyFraction", (double)o.chargedEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronEnergyFraction", (double)o.chargedHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronMultiplicity", (double)o.chargedHadronMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedMuEnergy", (double)o.chargedMuEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedMultiplicity", (double)o.chargedMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexBJetTags", (double)o.combinedSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexMVABJetTags", (double)o.combinedSecondaryVertexMVABJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "electronEnergy", (double)o.electronEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "electronMultiplicity", (double)o.electronMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genJet_energy", (double)o.genJet_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "genJet_eta", (double)o.genJet_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genJet_invisibleEnergy", (double)o.genJet_invisibleEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "genJet_phi", (double)o.genJet_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genJet_pt", (double)o.genJet_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "genParton_energy", (double)o.genParton_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "genParton_eta", (double)o.genParton_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genParton_pdgId", (double)o.genParton_pdgId); os << r;
  sprintf(r, "  %-32s: %f\n", "genParton_phi", (double)o.genParton_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genParton_pt", (double)o.genParton_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "jetArea", (double)o.jetArea); os << r;
  sprintf(r, "  %-32s: %f\n", "jetBProbabilityBJetTags", (double)o.jetBProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetID_fHPD", (double)o.jetID_fHPD); os << r;
  sprintf(r, "  %-32s: %f\n", "jetID_n90Hits", (double)o.jetID_n90Hits); os << r;
  sprintf(r, "  %-32s: %f\n", "jetProbabilityBJetTags", (double)o.jetProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "muonEnergy", (double)o.muonEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralEmEnergyFraction", (double)o.neutralEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergy", (double)o.neutralHadronEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergyFraction", (double)o.neutralHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronMultiplicity", (double)o.neutralHadronMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralMultiplicity", (double)o.neutralMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfDaughters", (double)o.numberOfDaughters); os << r;
  sprintf(r, "  %-32s: %f\n", "partonFlavour", (double)o.partonFlavour); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "photonEnergy", (double)o.photonEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "photonMultiplicity", (double)o.photonMultiplicity); os << r;
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
  sprintf(r, "  %-32s: %f\n", "uncor_et", (double)o.uncor_et); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_eta", (double)o.uncor_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_phi", (double)o.uncor_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_pt", (double)o.uncor_pt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jet2_s
{
  double	HFEMEnergy;
  double	HFEMMultiplicity;
  double	HFHadronEnergy;
  double	HFHadronMultiplicity;
  double	chargedEmEnergyFraction;
  double	chargedHadronEnergyFraction;
  double	chargedHadronMultiplicity;
  double	chargedMuEnergy;
  double	chargedMultiplicity;
  double	combinedSecondaryVertexBJetTags;
  double	combinedSecondaryVertexMVABJetTags;
  double	electronEnergy;
  double	electronMultiplicity;
  double	energy;
  double	et;
  double	eta;
  double	genJet_energy;
  double	genJet_eta;
  double	genJet_invisibleEnergy;
  double	genJet_phi;
  double	genJet_pt;
  double	genParton_energy;
  double	genParton_eta;
  double	genParton_pdgId;
  double	genParton_phi;
  double	genParton_pt;
  double	jetArea;
  double	jetBProbabilityBJetTags;
  double	jetID_fHPD;
  double	jetID_n90Hits;
  double	jetProbabilityBJetTags;
  double	muonEnergy;
  double	neutralEmEnergyFraction;
  double	neutralHadronEnergy;
  double	neutralHadronEnergyFraction;
  double	neutralHadronMultiplicity;
  double	neutralMultiplicity;
  double	numberOfDaughters;
  double	partonFlavour;
  double	phi;
  double	photonEnergy;
  double	photonMultiplicity;
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
  double	uncor_et;
  double	uncor_eta;
  double	uncor_phi;
  double	uncor_pt;
};
std::vector<jet2_s> jet2(177);

std::ostream& operator<<(std::ostream& os, const jet2_s& o)
{
  char r[1024];
  os << "jet2" << std::endl;
  sprintf(r, "  %-32s: %f\n", "HFEMEnergy", (double)o.HFEMEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "HFEMMultiplicity", (double)o.HFEMMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "HFHadronEnergy", (double)o.HFHadronEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "HFHadronMultiplicity", (double)o.HFHadronMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedEmEnergyFraction", (double)o.chargedEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronEnergyFraction", (double)o.chargedHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronMultiplicity", (double)o.chargedHadronMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedMuEnergy", (double)o.chargedMuEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedMultiplicity", (double)o.chargedMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexBJetTags", (double)o.combinedSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexMVABJetTags", (double)o.combinedSecondaryVertexMVABJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "electronEnergy", (double)o.electronEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "electronMultiplicity", (double)o.electronMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genJet_energy", (double)o.genJet_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "genJet_eta", (double)o.genJet_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genJet_invisibleEnergy", (double)o.genJet_invisibleEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "genJet_phi", (double)o.genJet_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genJet_pt", (double)o.genJet_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "genParton_energy", (double)o.genParton_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "genParton_eta", (double)o.genParton_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genParton_pdgId", (double)o.genParton_pdgId); os << r;
  sprintf(r, "  %-32s: %f\n", "genParton_phi", (double)o.genParton_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genParton_pt", (double)o.genParton_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "jetArea", (double)o.jetArea); os << r;
  sprintf(r, "  %-32s: %f\n", "jetBProbabilityBJetTags", (double)o.jetBProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetID_fHPD", (double)o.jetID_fHPD); os << r;
  sprintf(r, "  %-32s: %f\n", "jetID_n90Hits", (double)o.jetID_n90Hits); os << r;
  sprintf(r, "  %-32s: %f\n", "jetProbabilityBJetTags", (double)o.jetProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "muonEnergy", (double)o.muonEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralEmEnergyFraction", (double)o.neutralEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergy", (double)o.neutralHadronEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergyFraction", (double)o.neutralHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronMultiplicity", (double)o.neutralHadronMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralMultiplicity", (double)o.neutralMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfDaughters", (double)o.numberOfDaughters); os << r;
  sprintf(r, "  %-32s: %f\n", "partonFlavour", (double)o.partonFlavour); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "photonEnergy", (double)o.photonEnergy); os << r;
  sprintf(r, "  %-32s: %f\n", "photonMultiplicity", (double)o.photonMultiplicity); os << r;
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
  sprintf(r, "  %-32s: %f\n", "uncor_et", (double)o.uncor_et); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_eta", (double)o.uncor_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_phi", (double)o.uncor_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_pt", (double)o.uncor_pt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jethelper_s
{
  double	jecFactor;
  double	jecFactorNoL1Fast;
  double	jetUncMinus;
  double	jetUncPlus;
  double	nSecondaryVertices;
  double	secondaryVertex0DirectionX;
  double	secondaryVertex0DirectionY;
  double	secondaryVertex0DirectionZ;
  double	secondaryVertex0FlightDistance;
  double	secondaryVertex0FlightDistanceSig;
  double	secondaryVertex0Mass;
  double	secondaryVertex1DirectionX;
  double	secondaryVertex1DirectionY;
  double	secondaryVertex1DirectionZ;
  double	secondaryVertex1FlightDistance;
  double	secondaryVertex1FlightDistanceSig;
  double	secondaryVertex1Mass;
};
std::vector<jethelper_s> jethelper(205);

std::ostream& operator<<(std::ostream& os, const jethelper_s& o)
{
  char r[1024];
  os << "jethelper" << std::endl;
  sprintf(r, "  %-32s: %f\n", "jecFactor", (double)o.jecFactor); os << r;
  sprintf(r, "  %-32s: %f\n", "jecFactorNoL1Fast", (double)o.jecFactorNoL1Fast); os << r;
  sprintf(r, "  %-32s: %f\n", "jetUncMinus", (double)o.jetUncMinus); os << r;
  sprintf(r, "  %-32s: %f\n", "jetUncPlus", (double)o.jetUncPlus); os << r;
  sprintf(r, "  %-32s: %f\n", "nSecondaryVertices", (double)o.nSecondaryVertices); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0DirectionX", (double)o.secondaryVertex0DirectionX); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0DirectionY", (double)o.secondaryVertex0DirectionY); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0DirectionZ", (double)o.secondaryVertex0DirectionZ); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0FlightDistance", (double)o.secondaryVertex0FlightDistance); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0FlightDistanceSig", (double)o.secondaryVertex0FlightDistanceSig); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0Mass", (double)o.secondaryVertex0Mass); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1DirectionX", (double)o.secondaryVertex1DirectionX); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1DirectionY", (double)o.secondaryVertex1DirectionY); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1DirectionZ", (double)o.secondaryVertex1DirectionZ); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1FlightDistance", (double)o.secondaryVertex1FlightDistance); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1FlightDistanceSig", (double)o.secondaryVertex1FlightDistanceSig); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1Mass", (double)o.secondaryVertex1Mass); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jethelper1_s
{
  double	jecFactor;
  double	jecFactorNoL1Fast;
  double	jetUncMinus;
  double	jetUncPlus;
  double	nSecondaryVertices;
  double	secondaryVertex0DirectionX;
  double	secondaryVertex0DirectionY;
  double	secondaryVertex0DirectionZ;
  double	secondaryVertex0FlightDistance;
  double	secondaryVertex0FlightDistanceSig;
  double	secondaryVertex0Mass;
  double	secondaryVertex1DirectionX;
  double	secondaryVertex1DirectionY;
  double	secondaryVertex1DirectionZ;
  double	secondaryVertex1FlightDistance;
  double	secondaryVertex1FlightDistanceSig;
  double	secondaryVertex1Mass;
};
std::vector<jethelper1_s> jethelper1(177);

std::ostream& operator<<(std::ostream& os, const jethelper1_s& o)
{
  char r[1024];
  os << "jethelper1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "jecFactor", (double)o.jecFactor); os << r;
  sprintf(r, "  %-32s: %f\n", "jecFactorNoL1Fast", (double)o.jecFactorNoL1Fast); os << r;
  sprintf(r, "  %-32s: %f\n", "jetUncMinus", (double)o.jetUncMinus); os << r;
  sprintf(r, "  %-32s: %f\n", "jetUncPlus", (double)o.jetUncPlus); os << r;
  sprintf(r, "  %-32s: %f\n", "nSecondaryVertices", (double)o.nSecondaryVertices); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0DirectionX", (double)o.secondaryVertex0DirectionX); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0DirectionY", (double)o.secondaryVertex0DirectionY); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0DirectionZ", (double)o.secondaryVertex0DirectionZ); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0FlightDistance", (double)o.secondaryVertex0FlightDistance); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0FlightDistanceSig", (double)o.secondaryVertex0FlightDistanceSig); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0Mass", (double)o.secondaryVertex0Mass); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1DirectionX", (double)o.secondaryVertex1DirectionX); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1DirectionY", (double)o.secondaryVertex1DirectionY); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1DirectionZ", (double)o.secondaryVertex1DirectionZ); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1FlightDistance", (double)o.secondaryVertex1FlightDistance); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1FlightDistanceSig", (double)o.secondaryVertex1FlightDistanceSig); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1Mass", (double)o.secondaryVertex1Mass); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jethelper2_s
{
  double	jecFactor;
  double	jecFactorNoL1Fast;
  double	jetUncMinus;
  double	jetUncPlus;
  double	nSecondaryVertices;
  double	secondaryVertex0DirectionX;
  double	secondaryVertex0DirectionY;
  double	secondaryVertex0DirectionZ;
  double	secondaryVertex0FlightDistance;
  double	secondaryVertex0FlightDistanceSig;
  double	secondaryVertex0Mass;
  double	secondaryVertex1DirectionX;
  double	secondaryVertex1DirectionY;
  double	secondaryVertex1DirectionZ;
  double	secondaryVertex1FlightDistance;
  double	secondaryVertex1FlightDistanceSig;
  double	secondaryVertex1Mass;
};
std::vector<jethelper2_s> jethelper2(177);

std::ostream& operator<<(std::ostream& os, const jethelper2_s& o)
{
  char r[1024];
  os << "jethelper2" << std::endl;
  sprintf(r, "  %-32s: %f\n", "jecFactor", (double)o.jecFactor); os << r;
  sprintf(r, "  %-32s: %f\n", "jecFactorNoL1Fast", (double)o.jecFactorNoL1Fast); os << r;
  sprintf(r, "  %-32s: %f\n", "jetUncMinus", (double)o.jetUncMinus); os << r;
  sprintf(r, "  %-32s: %f\n", "jetUncPlus", (double)o.jetUncPlus); os << r;
  sprintf(r, "  %-32s: %f\n", "nSecondaryVertices", (double)o.nSecondaryVertices); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0DirectionX", (double)o.secondaryVertex0DirectionX); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0DirectionY", (double)o.secondaryVertex0DirectionY); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0DirectionZ", (double)o.secondaryVertex0DirectionZ); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0FlightDistance", (double)o.secondaryVertex0FlightDistance); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0FlightDistanceSig", (double)o.secondaryVertex0FlightDistanceSig); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex0Mass", (double)o.secondaryVertex0Mass); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1DirectionX", (double)o.secondaryVertex1DirectionX); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1DirectionY", (double)o.secondaryVertex1DirectionY); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1DirectionZ", (double)o.secondaryVertex1DirectionZ); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1FlightDistance", (double)o.secondaryVertex1FlightDistance); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1FlightDistanceSig", (double)o.secondaryVertex1FlightDistanceSig); os << r;
  sprintf(r, "  %-32s: %f\n", "secondaryVertex1Mass", (double)o.secondaryVertex1Mass); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct met_s
{
  double	energy;
  double	et;
  double	genMET_energy;
  double	genMET_et;
  double	genMET_phi;
  double	genMET_pt;
  double	mEtSig;
  double	metSignificance;
  double	phi;
  double	pt;
  double	significance;
  double	sumEt;
};
std::vector<met_s> met(3);

std::ostream& operator<<(std::ostream& os, const met_s& o)
{
  char r[1024];
  os << "met" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_energy", (double)o.genMET_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_et", (double)o.genMET_et); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_phi", (double)o.genMET_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_pt", (double)o.genMET_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  sprintf(r, "  %-32s: %f\n", "metSignificance", (double)o.metSignificance); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "significance", (double)o.significance); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct met1_s
{
  double	energy;
  double	et;
  double	genMET_energy;
  double	genMET_et;
  double	genMET_phi;
  double	genMET_pt;
  double	mEtSig;
  double	phi;
  double	pt;
  double	significance;
  double	sumEt;
};
std::vector<met1_s> met1(3);

std::ostream& operator<<(std::ostream& os, const met1_s& o)
{
  char r[1024];
  os << "met1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_energy", (double)o.genMET_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_et", (double)o.genMET_et); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_phi", (double)o.genMET_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_pt", (double)o.genMET_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "significance", (double)o.significance); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct met2_s
{
  double	energy;
  double	et;
  double	genMET_energy;
  double	genMET_et;
  double	genMET_phi;
  double	genMET_pt;
  double	mEtSig;
  double	phi;
  double	pt;
  double	significance;
  double	sumEt;
};
std::vector<met2_s> met2(3);

std::ostream& operator<<(std::ostream& os, const met2_s& o)
{
  char r[1024];
  os << "met2" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_energy", (double)o.genMET_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_et", (double)o.genMET_et); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_phi", (double)o.genMET_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_pt", (double)o.genMET_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "significance", (double)o.significance); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct met3_s
{
  double	energy;
  double	et;
  double	genMET_energy;
  double	genMET_et;
  double	genMET_phi;
  double	genMET_pt;
  double	mEtSig;
  double	phi;
  double	pt;
  double	significance;
  double	sumEt;
};
std::vector<met3_s> met3(3);

std::ostream& operator<<(std::ostream& os, const met3_s& o)
{
  char r[1024];
  os << "met3" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_energy", (double)o.genMET_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_et", (double)o.genMET_et); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_phi", (double)o.genMET_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_pt", (double)o.genMET_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "significance", (double)o.significance); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct met4_s
{
  double	energy;
  double	et;
  double	genMET_energy;
  double	genMET_et;
  double	genMET_phi;
  double	genMET_pt;
  double	mEtSig;
  double	phi;
  double	pt;
  double	significance;
  double	sumEt;
};
std::vector<met4_s> met4(3);

std::ostream& operator<<(std::ostream& os, const met4_s& o)
{
  char r[1024];
  os << "met4" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_energy", (double)o.genMET_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_et", (double)o.genMET_et); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_phi", (double)o.genMET_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_pt", (double)o.genMET_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "significance", (double)o.significance); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct met5_s
{
  double	energy;
  double	et;
  double	genMET_energy;
  double	genMET_et;
  double	genMET_phi;
  double	genMET_pt;
  double	mEtSig;
  double	phi;
  double	pt;
  double	significance;
  double	sumEt;
};
std::vector<met5_s> met5(3);

std::ostream& operator<<(std::ostream& os, const met5_s& o)
{
  char r[1024];
  os << "met5" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_energy", (double)o.genMET_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_et", (double)o.genMET_et); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_phi", (double)o.genMET_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_pt", (double)o.genMET_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "significance", (double)o.significance); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct met6_s
{
  double	energy;
  double	et;
  double	genMET_energy;
  double	genMET_et;
  double	genMET_phi;
  double	genMET_pt;
  double	mEtSig;
  double	phi;
  double	pt;
  double	significance;
  double	sumEt;
};
std::vector<met6_s> met6(3);

std::ostream& operator<<(std::ostream& os, const met6_s& o)
{
  char r[1024];
  os << "met6" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_energy", (double)o.genMET_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_et", (double)o.genMET_et); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_phi", (double)o.genMET_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_pt", (double)o.genMET_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "significance", (double)o.significance); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct met7_s
{
  double	energy;
  double	et;
  double	genMET_energy;
  double	genMET_et;
  double	genMET_phi;
  double	genMET_pt;
  double	mEtSig;
  double	phi;
  double	pt;
  double	significance;
  double	sumEt;
};
std::vector<met7_s> met7(3);

std::ostream& operator<<(std::ostream& os, const met7_s& o)
{
  char r[1024];
  os << "met7" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_energy", (double)o.genMET_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_et", (double)o.genMET_et); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_phi", (double)o.genMET_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genMET_pt", (double)o.genMET_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "significance", (double)o.significance); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct methelper_s
{
  double	significance_dxx;
  double	significance_dxy;
  double	significance_dyx;
  double	significance_dyy;
};
std::vector<methelper_s> methelper(3);

std::ostream& operator<<(std::ostream& os, const methelper_s& o)
{
  char r[1024];
  os << "methelper" << std::endl;
  sprintf(r, "  %-32s: %f\n", "significance_dxx", (double)o.significance_dxx); os << r;
  sprintf(r, "  %-32s: %f\n", "significance_dxy", (double)o.significance_dxy); os << r;
  sprintf(r, "  %-32s: %f\n", "significance_dyx", (double)o.significance_dyx); os << r;
  sprintf(r, "  %-32s: %f\n", "significance_dyy", (double)o.significance_dyy); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct methelper1_s
{
  double	significance_dxx;
  double	significance_dxy;
  double	significance_dyx;
  double	significance_dyy;
};
std::vector<methelper1_s> methelper1(3);

std::ostream& operator<<(std::ostream& os, const methelper1_s& o)
{
  char r[1024];
  os << "methelper1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "significance_dxx", (double)o.significance_dxx); os << r;
  sprintf(r, "  %-32s: %f\n", "significance_dxy", (double)o.significance_dxy); os << r;
  sprintf(r, "  %-32s: %f\n", "significance_dyx", (double)o.significance_dyx); os << r;
  sprintf(r, "  %-32s: %f\n", "significance_dyy", (double)o.significance_dyy); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct muon_s
{
  double	AllGlobalMuons;
  double	AllTrackerMuons;
  double	GlobalMuonPromptTight;
  double	charge;
  double	chargedHadronIso;
  double	combinedMuon_chi2;
  double	combinedMuon_ndof;
  double	combinedMuon_numberOfValidHits;
  double	dB;
  double	ecalIso;
  double	energy;
  double	et;
  double	eta;
  double	genLepton_eta;
  double	genLepton_pdgId;
  double	genLepton_phi;
  double	genLepton_pt;
  double	globalTrack_chi2;
  double	globalTrack_hitPattern_numberOfValidTrackerHits;
  double	globalTrack_ndof;
  double	globalTrack_pt;
  double	globalTrack_ptError;
  double	hcalIso;
  double	innerTrack_d0;
  double	innerTrack_hitPattern_pixelLayersWithMeasurement;
  double	innerTrack_numberOfValidHits;
  double	innerTrack_phi;
  double	innerTrack_pt;
  double	innerTrack_vertex_z;
  double	isGlobalMuon;
  double	isTrackerMuon;
  double	neutralHadronIso;
  double	numberOfMatches;
  double	phi;
  double	photonIso;
  double	pt;
  double	trackIso;
  double	track_d0;
  double	track_hitPattern_numberOfValidMuonHits;
  double	track_hitPattern_numberOfValidPixelHits;
  double	track_normalizedChi2;
  double	vz;
};
std::vector<muon_s> muon(29);

std::ostream& operator<<(std::ostream& os, const muon_s& o)
{
  char r[1024];
  os << "muon" << std::endl;
  sprintf(r, "  %-32s: %f\n", "AllGlobalMuons", (double)o.AllGlobalMuons); os << r;
  sprintf(r, "  %-32s: %f\n", "AllTrackerMuons", (double)o.AllTrackerMuons); os << r;
  sprintf(r, "  %-32s: %f\n", "GlobalMuonPromptTight", (double)o.GlobalMuonPromptTight); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronIso", (double)o.chargedHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedMuon_chi2", (double)o.combinedMuon_chi2); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedMuon_ndof", (double)o.combinedMuon_ndof); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedMuon_numberOfValidHits", (double)o.combinedMuon_numberOfValidHits); os << r;
  sprintf(r, "  %-32s: %f\n", "dB", (double)o.dB); os << r;
  sprintf(r, "  %-32s: %f\n", "ecalIso", (double)o.ecalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_eta", (double)o.genLepton_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_pdgId", (double)o.genLepton_pdgId); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_phi", (double)o.genLepton_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_pt", (double)o.genLepton_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_chi2", (double)o.globalTrack_chi2); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_hitPattern_numberOfValidTrackerHits", (double)o.globalTrack_hitPattern_numberOfValidTrackerHits); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_ndof", (double)o.globalTrack_ndof); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_pt", (double)o.globalTrack_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_ptError", (double)o.globalTrack_ptError); os << r;
  sprintf(r, "  %-32s: %f\n", "hcalIso", (double)o.hcalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_d0", (double)o.innerTrack_d0); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_hitPattern_pixelLayersWithMeasurement", (double)o.innerTrack_hitPattern_pixelLayersWithMeasurement); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_numberOfValidHits", (double)o.innerTrack_numberOfValidHits); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_phi", (double)o.innerTrack_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_pt", (double)o.innerTrack_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_vertex_z", (double)o.innerTrack_vertex_z); os << r;
  sprintf(r, "  %-32s: %f\n", "isGlobalMuon", (double)o.isGlobalMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "isTrackerMuon", (double)o.isTrackerMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronIso", (double)o.neutralHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfMatches", (double)o.numberOfMatches); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "photonIso", (double)o.photonIso); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "trackIso", (double)o.trackIso); os << r;
  sprintf(r, "  %-32s: %f\n", "track_d0", (double)o.track_d0); os << r;
  sprintf(r, "  %-32s: %f\n", "track_hitPattern_numberOfValidMuonHits", (double)o.track_hitPattern_numberOfValidMuonHits); os << r;
  sprintf(r, "  %-32s: %f\n", "track_hitPattern_numberOfValidPixelHits", (double)o.track_hitPattern_numberOfValidPixelHits); os << r;
  sprintf(r, "  %-32s: %f\n", "track_normalizedChi2", (double)o.track_normalizedChi2); os << r;
  sprintf(r, "  %-32s: %f\n", "vz", (double)o.vz); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct muon1_s
{
  double	AllGlobalMuons;
  double	AllTrackerMuons;
  double	GlobalMuonPromptTight;
  double	charge;
  double	chargedHadronIso;
  double	combinedMuon_chi2;
  double	combinedMuon_ndof;
  double	combinedMuon_numberOfValidHits;
  double	dB;
  double	ecalIso;
  double	energy;
  double	et;
  double	eta;
  double	genLepton_eta;
  double	genLepton_pdgId;
  double	genLepton_phi;
  double	genLepton_pt;
  double	globalTrack_chi2;
  double	globalTrack_hitPattern_numberOfValidTrackerHits;
  double	globalTrack_ndof;
  double	globalTrack_pt;
  double	globalTrack_ptError;
  double	hcalIso;
  double	innerTrack_d0;
  double	innerTrack_hitPattern_pixelLayersWithMeasurement;
  double	innerTrack_numberOfValidHits;
  double	innerTrack_phi;
  double	innerTrack_pt;
  double	innerTrack_vertex_z;
  double	isGlobalMuon;
  double	isTrackerMuon;
  double	neutralHadronIso;
  double	numberOfMatches;
  double	phi;
  double	photonIso;
  double	pt;
  double	trackIso;
  double	track_d0;
  double	track_hitPattern_numberOfValidMuonHits;
  double	track_hitPattern_numberOfValidPixelHits;
  double	track_normalizedChi2;
  double	vz;
};
std::vector<muon1_s> muon1(7);

std::ostream& operator<<(std::ostream& os, const muon1_s& o)
{
  char r[1024];
  os << "muon1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "AllGlobalMuons", (double)o.AllGlobalMuons); os << r;
  sprintf(r, "  %-32s: %f\n", "AllTrackerMuons", (double)o.AllTrackerMuons); os << r;
  sprintf(r, "  %-32s: %f\n", "GlobalMuonPromptTight", (double)o.GlobalMuonPromptTight); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronIso", (double)o.chargedHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedMuon_chi2", (double)o.combinedMuon_chi2); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedMuon_ndof", (double)o.combinedMuon_ndof); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedMuon_numberOfValidHits", (double)o.combinedMuon_numberOfValidHits); os << r;
  sprintf(r, "  %-32s: %f\n", "dB", (double)o.dB); os << r;
  sprintf(r, "  %-32s: %f\n", "ecalIso", (double)o.ecalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_eta", (double)o.genLepton_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_pdgId", (double)o.genLepton_pdgId); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_phi", (double)o.genLepton_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "genLepton_pt", (double)o.genLepton_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_chi2", (double)o.globalTrack_chi2); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_hitPattern_numberOfValidTrackerHits", (double)o.globalTrack_hitPattern_numberOfValidTrackerHits); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_ndof", (double)o.globalTrack_ndof); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_pt", (double)o.globalTrack_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_ptError", (double)o.globalTrack_ptError); os << r;
  sprintf(r, "  %-32s: %f\n", "hcalIso", (double)o.hcalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_d0", (double)o.innerTrack_d0); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_hitPattern_pixelLayersWithMeasurement", (double)o.innerTrack_hitPattern_pixelLayersWithMeasurement); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_numberOfValidHits", (double)o.innerTrack_numberOfValidHits); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_phi", (double)o.innerTrack_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_pt", (double)o.innerTrack_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_vertex_z", (double)o.innerTrack_vertex_z); os << r;
  sprintf(r, "  %-32s: %f\n", "isGlobalMuon", (double)o.isGlobalMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "isTrackerMuon", (double)o.isTrackerMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronIso", (double)o.neutralHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfMatches", (double)o.numberOfMatches); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "photonIso", (double)o.photonIso); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "trackIso", (double)o.trackIso); os << r;
  sprintf(r, "  %-32s: %f\n", "track_d0", (double)o.track_d0); os << r;
  sprintf(r, "  %-32s: %f\n", "track_hitPattern_numberOfValidMuonHits", (double)o.track_hitPattern_numberOfValidMuonHits); os << r;
  sprintf(r, "  %-32s: %f\n", "track_hitPattern_numberOfValidPixelHits", (double)o.track_hitPattern_numberOfValidPixelHits); os << r;
  sprintf(r, "  %-32s: %f\n", "track_normalizedChi2", (double)o.track_normalizedChi2); os << r;
  sprintf(r, "  %-32s: %f\n", "vz", (double)o.vz); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct muonhelper_s
{
  double	dxywrtBeamSpot;
  double	muonPFcandptdiff;
};
std::vector<muonhelper_s> muonhelper(29);

std::ostream& operator<<(std::ostream& os, const muonhelper_s& o)
{
  char r[1024];
  os << "muonhelper" << std::endl;
  sprintf(r, "  %-32s: %f\n", "dxywrtBeamSpot", (double)o.dxywrtBeamSpot); os << r;
  sprintf(r, "  %-32s: %f\n", "muonPFcandptdiff", (double)o.muonPFcandptdiff); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct muonhelper1_s
{
  double	dxywrtBeamSpot;
  double	muonPFcandptdiff;
};
std::vector<muonhelper1_s> muonhelper1(7);

std::ostream& operator<<(std::ostream& os, const muonhelper1_s& o)
{
  char r[1024];
  os << "muonhelper1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "dxywrtBeamSpot", (double)o.dxywrtBeamSpot); os << r;
  sprintf(r, "  %-32s: %f\n", "muonPFcandptdiff", (double)o.muonPFcandptdiff); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct pileupsummaryinfo_s
{
  double	addpileupinfo_getBunchCrossing;
  double	addpileupinfo_getPU_NumInteractions;
};
std::vector<pileupsummaryinfo_s> pileupsummaryinfo(7);

std::ostream& operator<<(std::ostream& os, const pileupsummaryinfo_s& o)
{
  char r[1024];
  os << "pileupsummaryinfo" << std::endl;
  sprintf(r, "  %-32s: %f\n", "addpileupinfo_getBunchCrossing", (double)o.addpileupinfo_getBunchCrossing); os << r;
  sprintf(r, "  %-32s: %f\n", "addpileupinfo_getPU_NumInteractions", (double)o.addpileupinfo_getPU_NumInteractions); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct tau_s
{
  double	caloIso;
  double	ecalIso;
  double	energy;
  double	et;
  double	eta;
  double	hcalIso;
  double	phi;
  double	pt;
  double	tauID_againstElectron;
  double	tauID_againstMuon;
  double	tauID_byIsolation;
  double	tauID_byTaNC;
  double	tauID_byTaNCfrHalfPercent;
  double	tauID_byTaNCfrQuarterPercent;
  double	trackIso;
};
std::vector<tau_s> tau(7);

std::ostream& operator<<(std::ostream& os, const tau_s& o)
{
  char r[1024];
  os << "tau" << std::endl;
  sprintf(r, "  %-32s: %f\n", "caloIso", (double)o.caloIso); os << r;
  sprintf(r, "  %-32s: %f\n", "ecalIso", (double)o.ecalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "hcalIso", (double)o.hcalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
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
struct tauhelper_s
{
  double	genTauDecayModeID;
};
std::vector<tauhelper_s> tauhelper(7);

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
std::vector<vertex_s> vertex(41);

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

  electron.resize(electron_caloIso.size());
  for(unsigned int i=0; i < electron_caloIso.size(); ++i)
    {
      electron[i].caloIso	= electron_caloIso[i];
      electron[i].charge	= electron_charge[i];
      electron[i].chargedHadronIso	= electron_chargedHadronIso[i];
      electron[i].convDcot	= electron_convDcot[i];
      electron[i].convDist	= electron_convDist[i];
      electron[i].dB	= electron_dB[i];
      electron[i].deltaEtaSuperClusterTrackAtVtx	= electron_deltaEtaSuperClusterTrackAtVtx[i];
      electron[i].deltaPhiSuperClusterTrackAtVtx	= electron_deltaPhiSuperClusterTrackAtVtx[i];
      electron[i].dr03EcalRecHitSumEt	= electron_dr03EcalRecHitSumEt[i];
      electron[i].dr03HcalTowerSumEt	= electron_dr03HcalTowerSumEt[i];
      electron[i].dr03TkSumPt	= electron_dr03TkSumPt[i];
      electron[i].ecalIso	= electron_ecalIso[i];
      electron[i].eidRobustTight	= electron_eidRobustTight[i];
      electron[i].energy	= electron_energy[i];
      electron[i].et	= electron_et[i];
      electron[i].eta	= electron_eta[i];
      electron[i].genLepton_eta	= electron_genLepton_eta[i];
      electron[i].genLepton_pdgId	= electron_genLepton_pdgId[i];
      electron[i].genLepton_phi	= electron_genLepton_phi[i];
      electron[i].genLepton_pt	= electron_genLepton_pt[i];
      electron[i].gsfTrack_d0	= electron_gsfTrack_d0[i];
      electron[i].gsfTrack_trackerExpectedHitsInner_numberOfLostHits	= electron_gsfTrack_trackerExpectedHitsInner_numberOfLostHits[i];
      electron[i].hadronicOverEm	= electron_hadronicOverEm[i];
      electron[i].hcalIso	= electron_hcalIso[i];
      electron[i].neutralHadronIso	= electron_neutralHadronIso[i];
      electron[i].phi	= electron_phi[i];
      electron[i].photonIso	= electron_photonIso[i];
      electron[i].pt	= electron_pt[i];
      electron[i].sigmaIetaIeta	= electron_sigmaIetaIeta[i];
      electron[i].simpleEleId80cIso	= electron_simpleEleId80cIso[i];
      electron[i].simpleEleId80relIso	= electron_simpleEleId80relIso[i];
      electron[i].simpleEleId95cIso	= electron_simpleEleId95cIso[i];
      electron[i].simpleEleId95relIso	= electron_simpleEleId95relIso[i];
      electron[i].superCluster_eta	= electron_superCluster_eta[i];
      electron[i].trackIso	= electron_trackIso[i];
      electron[i].vertex_z	= electron_vertex_z[i];
      electron[i].vz	= electron_vz[i];
    }

  electron1.resize(electron1_caloIso.size());
  for(unsigned int i=0; i < electron1_caloIso.size(); ++i)
    {
      electron1[i].caloIso	= electron1_caloIso[i];
      electron1[i].charge	= electron1_charge[i];
      electron1[i].chargedHadronIso	= electron1_chargedHadronIso[i];
      electron1[i].convDcot	= electron1_convDcot[i];
      electron1[i].convDist	= electron1_convDist[i];
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
      electron1[i].sigmaIetaIeta	= electron1_sigmaIetaIeta[i];
      electron1[i].superCluster_eta	= electron1_superCluster_eta[i];
      electron1[i].trackIso	= electron1_trackIso[i];
      electron1[i].vertex_z	= electron1_vertex_z[i];
      electron1[i].vz	= electron1_vz[i];
    }

  electronhelper.resize(electronhelper_dxywrtBeamSpot.size());
  for(unsigned int i=0; i < electronhelper_dxywrtBeamSpot.size(); ++i)
    {
      electronhelper[i].dxywrtBeamSpot	= electronhelper_dxywrtBeamSpot[i];
    }

  electronhelper1.resize(electronhelper1_dxywrtBeamSpot.size());
  for(unsigned int i=0; i < electronhelper1_dxywrtBeamSpot.size(); ++i)
    {
      electronhelper1[i].dxywrtBeamSpot	= electronhelper1_dxywrtBeamSpot[i];
    }

  geneventinfoproducthelper.resize(geneventinfoproducthelper_pdf1.size());
  for(unsigned int i=0; i < geneventinfoproducthelper_pdf1.size(); ++i)
    {
      geneventinfoproducthelper[i].pdf1	= geneventinfoproducthelper_pdf1[i];
      geneventinfoproducthelper[i].pdf2	= geneventinfoproducthelper_pdf2[i];
      geneventinfoproducthelper[i].pdfweight	= geneventinfoproducthelper_pdfweight[i];
      geneventinfoproducthelper[i].pdfweightsum	= geneventinfoproducthelper_pdfweightsum[i];
    }

  geneventinfoproducthelper1.resize(geneventinfoproducthelper1_pdf1.size());
  for(unsigned int i=0; i < geneventinfoproducthelper1_pdf1.size(); ++i)
    {
      geneventinfoproducthelper1[i].pdf1	= geneventinfoproducthelper1_pdf1[i];
      geneventinfoproducthelper1[i].pdf2	= geneventinfoproducthelper1_pdf2[i];
      geneventinfoproducthelper1[i].pdfweight	= geneventinfoproducthelper1_pdfweight[i];
      geneventinfoproducthelper1[i].pdfweightsum	= geneventinfoproducthelper1_pdfweightsum[i];
    }

  geneventinfoproducthelper2.resize(geneventinfoproducthelper2_pdf1.size());
  for(unsigned int i=0; i < geneventinfoproducthelper2_pdf1.size(); ++i)
    {
      geneventinfoproducthelper2[i].pdf1	= geneventinfoproducthelper2_pdf1[i];
      geneventinfoproducthelper2[i].pdf2	= geneventinfoproducthelper2_pdf2[i];
      geneventinfoproducthelper2[i].pdfweight	= geneventinfoproducthelper2_pdfweight[i];
      geneventinfoproducthelper2[i].pdfweightsum	= geneventinfoproducthelper2_pdfweightsum[i];
    }

  genparticlehelperra2.resize(genparticlehelperra2_charge.size());
  for(unsigned int i=0; i < genparticlehelperra2_charge.size(); ++i)
    {
      genparticlehelperra2[i].charge	= genparticlehelperra2_charge[i];
      genparticlehelperra2[i].eta	= genparticlehelperra2_eta[i];
      genparticlehelperra2[i].firstDaughter	= genparticlehelperra2_firstDaughter[i];
      genparticlehelperra2[i].firstMother	= genparticlehelperra2_firstMother[i];
      genparticlehelperra2[i].lastDaughter	= genparticlehelperra2_lastDaughter[i];
      genparticlehelperra2[i].lastMother	= genparticlehelperra2_lastMother[i];
      genparticlehelperra2[i].mass	= genparticlehelperra2_mass[i];
      genparticlehelperra2[i].pdgId	= genparticlehelperra2_pdgId[i];
      genparticlehelperra2[i].phi	= genparticlehelperra2_phi[i];
      genparticlehelperra2[i].pt	= genparticlehelperra2_pt[i];
      genparticlehelperra2[i].status	= genparticlehelperra2_status[i];
    }

  jet.resize(jet_HFEMEnergy.size());
  for(unsigned int i=0; i < jet_HFEMEnergy.size(); ++i)
    {
      jet[i].HFEMEnergy	= jet_HFEMEnergy[i];
      jet[i].HFEMMultiplicity	= jet_HFEMMultiplicity[i];
      jet[i].HFHadronEnergy	= jet_HFHadronEnergy[i];
      jet[i].HFHadronMultiplicity	= jet_HFHadronMultiplicity[i];
      jet[i].chargedEmEnergyFraction	= jet_chargedEmEnergyFraction[i];
      jet[i].chargedHadronEnergyFraction	= jet_chargedHadronEnergyFraction[i];
      jet[i].chargedHadronMultiplicity	= jet_chargedHadronMultiplicity[i];
      jet[i].chargedMuEnergy	= jet_chargedMuEnergy[i];
      jet[i].chargedMultiplicity	= jet_chargedMultiplicity[i];
      jet[i].combinedSecondaryVertexBJetTags	= jet_combinedSecondaryVertexBJetTags[i];
      jet[i].combinedSecondaryVertexMVABJetTags	= jet_combinedSecondaryVertexMVABJetTags[i];
      jet[i].electronEnergy	= jet_electronEnergy[i];
      jet[i].electronMultiplicity	= jet_electronMultiplicity[i];
      jet[i].energy	= jet_energy[i];
      jet[i].et	= jet_et[i];
      jet[i].eta	= jet_eta[i];
      jet[i].genJet_energy	= jet_genJet_energy[i];
      jet[i].genJet_eta	= jet_genJet_eta[i];
      jet[i].genJet_invisibleEnergy	= jet_genJet_invisibleEnergy[i];
      jet[i].genJet_phi	= jet_genJet_phi[i];
      jet[i].genJet_pt	= jet_genJet_pt[i];
      jet[i].genParton_energy	= jet_genParton_energy[i];
      jet[i].genParton_eta	= jet_genParton_eta[i];
      jet[i].genParton_pdgId	= jet_genParton_pdgId[i];
      jet[i].genParton_phi	= jet_genParton_phi[i];
      jet[i].genParton_pt	= jet_genParton_pt[i];
      jet[i].jetArea	= jet_jetArea[i];
      jet[i].jetBProbabilityBJetTags	= jet_jetBProbabilityBJetTags[i];
      jet[i].jetID_fHPD	= jet_jetID_fHPD[i];
      jet[i].jetID_n90Hits	= jet_jetID_n90Hits[i];
      jet[i].jetProbabilityBJetTags	= jet_jetProbabilityBJetTags[i];
      jet[i].muonEnergy	= jet_muonEnergy[i];
      jet[i].neutralEmEnergyFraction	= jet_neutralEmEnergyFraction[i];
      jet[i].neutralHadronEnergy	= jet_neutralHadronEnergy[i];
      jet[i].neutralHadronEnergyFraction	= jet_neutralHadronEnergyFraction[i];
      jet[i].neutralHadronMultiplicity	= jet_neutralHadronMultiplicity[i];
      jet[i].neutralMultiplicity	= jet_neutralMultiplicity[i];
      jet[i].numberOfDaughters	= jet_numberOfDaughters[i];
      jet[i].partonFlavour	= jet_partonFlavour[i];
      jet[i].phi	= jet_phi[i];
      jet[i].photonEnergy	= jet_photonEnergy[i];
      jet[i].photonMultiplicity	= jet_photonMultiplicity[i];
      jet[i].pt	= jet_pt[i];
      jet[i].simpleSecondaryVertexBJetTags	= jet_simpleSecondaryVertexBJetTags[i];
      jet[i].simpleSecondaryVertexHighEffBJetTags	= jet_simpleSecondaryVertexHighEffBJetTags[i];
      jet[i].simpleSecondaryVertexHighPurBJetTags	= jet_simpleSecondaryVertexHighPurBJetTags[i];
      jet[i].softElectronByIP3dBJetTags	= jet_softElectronByIP3dBJetTags[i];
      jet[i].softElectronByPtBJetTags	= jet_softElectronByPtBJetTags[i];
      jet[i].softMuonBJetTags	= jet_softMuonBJetTags[i];
      jet[i].softMuonByIP3dBJetTags	= jet_softMuonByIP3dBJetTags[i];
      jet[i].softMuonByPtBJetTags	= jet_softMuonByPtBJetTags[i];
      jet[i].trackCountingHighEffBJetTags	= jet_trackCountingHighEffBJetTags[i];
      jet[i].trackCountingHighPurBJetTags	= jet_trackCountingHighPurBJetTags[i];
      jet[i].uncor_energy	= jet_uncor_energy[i];
      jet[i].uncor_et	= jet_uncor_et[i];
      jet[i].uncor_eta	= jet_uncor_eta[i];
      jet[i].uncor_phi	= jet_uncor_phi[i];
      jet[i].uncor_pt	= jet_uncor_pt[i];
    }

  jet1.resize(jet1_HFEMEnergy.size());
  for(unsigned int i=0; i < jet1_HFEMEnergy.size(); ++i)
    {
      jet1[i].HFEMEnergy	= jet1_HFEMEnergy[i];
      jet1[i].HFEMMultiplicity	= jet1_HFEMMultiplicity[i];
      jet1[i].HFHadronEnergy	= jet1_HFHadronEnergy[i];
      jet1[i].HFHadronMultiplicity	= jet1_HFHadronMultiplicity[i];
      jet1[i].chargedEmEnergyFraction	= jet1_chargedEmEnergyFraction[i];
      jet1[i].chargedHadronEnergyFraction	= jet1_chargedHadronEnergyFraction[i];
      jet1[i].chargedHadronMultiplicity	= jet1_chargedHadronMultiplicity[i];
      jet1[i].chargedMuEnergy	= jet1_chargedMuEnergy[i];
      jet1[i].chargedMultiplicity	= jet1_chargedMultiplicity[i];
      jet1[i].combinedSecondaryVertexBJetTags	= jet1_combinedSecondaryVertexBJetTags[i];
      jet1[i].combinedSecondaryVertexMVABJetTags	= jet1_combinedSecondaryVertexMVABJetTags[i];
      jet1[i].electronEnergy	= jet1_electronEnergy[i];
      jet1[i].electronMultiplicity	= jet1_electronMultiplicity[i];
      jet1[i].energy	= jet1_energy[i];
      jet1[i].et	= jet1_et[i];
      jet1[i].eta	= jet1_eta[i];
      jet1[i].genJet_energy	= jet1_genJet_energy[i];
      jet1[i].genJet_eta	= jet1_genJet_eta[i];
      jet1[i].genJet_invisibleEnergy	= jet1_genJet_invisibleEnergy[i];
      jet1[i].genJet_phi	= jet1_genJet_phi[i];
      jet1[i].genJet_pt	= jet1_genJet_pt[i];
      jet1[i].genParton_energy	= jet1_genParton_energy[i];
      jet1[i].genParton_eta	= jet1_genParton_eta[i];
      jet1[i].genParton_pdgId	= jet1_genParton_pdgId[i];
      jet1[i].genParton_phi	= jet1_genParton_phi[i];
      jet1[i].genParton_pt	= jet1_genParton_pt[i];
      jet1[i].jetArea	= jet1_jetArea[i];
      jet1[i].jetBProbabilityBJetTags	= jet1_jetBProbabilityBJetTags[i];
      jet1[i].jetID_fHPD	= jet1_jetID_fHPD[i];
      jet1[i].jetID_n90Hits	= jet1_jetID_n90Hits[i];
      jet1[i].jetProbabilityBJetTags	= jet1_jetProbabilityBJetTags[i];
      jet1[i].muonEnergy	= jet1_muonEnergy[i];
      jet1[i].neutralEmEnergyFraction	= jet1_neutralEmEnergyFraction[i];
      jet1[i].neutralHadronEnergy	= jet1_neutralHadronEnergy[i];
      jet1[i].neutralHadronEnergyFraction	= jet1_neutralHadronEnergyFraction[i];
      jet1[i].neutralHadronMultiplicity	= jet1_neutralHadronMultiplicity[i];
      jet1[i].neutralMultiplicity	= jet1_neutralMultiplicity[i];
      jet1[i].numberOfDaughters	= jet1_numberOfDaughters[i];
      jet1[i].partonFlavour	= jet1_partonFlavour[i];
      jet1[i].phi	= jet1_phi[i];
      jet1[i].photonEnergy	= jet1_photonEnergy[i];
      jet1[i].photonMultiplicity	= jet1_photonMultiplicity[i];
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
      jet1[i].uncor_et	= jet1_uncor_et[i];
      jet1[i].uncor_eta	= jet1_uncor_eta[i];
      jet1[i].uncor_phi	= jet1_uncor_phi[i];
      jet1[i].uncor_pt	= jet1_uncor_pt[i];
    }

  jet2.resize(jet2_HFEMEnergy.size());
  for(unsigned int i=0; i < jet2_HFEMEnergy.size(); ++i)
    {
      jet2[i].HFEMEnergy	= jet2_HFEMEnergy[i];
      jet2[i].HFEMMultiplicity	= jet2_HFEMMultiplicity[i];
      jet2[i].HFHadronEnergy	= jet2_HFHadronEnergy[i];
      jet2[i].HFHadronMultiplicity	= jet2_HFHadronMultiplicity[i];
      jet2[i].chargedEmEnergyFraction	= jet2_chargedEmEnergyFraction[i];
      jet2[i].chargedHadronEnergyFraction	= jet2_chargedHadronEnergyFraction[i];
      jet2[i].chargedHadronMultiplicity	= jet2_chargedHadronMultiplicity[i];
      jet2[i].chargedMuEnergy	= jet2_chargedMuEnergy[i];
      jet2[i].chargedMultiplicity	= jet2_chargedMultiplicity[i];
      jet2[i].combinedSecondaryVertexBJetTags	= jet2_combinedSecondaryVertexBJetTags[i];
      jet2[i].combinedSecondaryVertexMVABJetTags	= jet2_combinedSecondaryVertexMVABJetTags[i];
      jet2[i].electronEnergy	= jet2_electronEnergy[i];
      jet2[i].electronMultiplicity	= jet2_electronMultiplicity[i];
      jet2[i].energy	= jet2_energy[i];
      jet2[i].et	= jet2_et[i];
      jet2[i].eta	= jet2_eta[i];
      jet2[i].genJet_energy	= jet2_genJet_energy[i];
      jet2[i].genJet_eta	= jet2_genJet_eta[i];
      jet2[i].genJet_invisibleEnergy	= jet2_genJet_invisibleEnergy[i];
      jet2[i].genJet_phi	= jet2_genJet_phi[i];
      jet2[i].genJet_pt	= jet2_genJet_pt[i];
      jet2[i].genParton_energy	= jet2_genParton_energy[i];
      jet2[i].genParton_eta	= jet2_genParton_eta[i];
      jet2[i].genParton_pdgId	= jet2_genParton_pdgId[i];
      jet2[i].genParton_phi	= jet2_genParton_phi[i];
      jet2[i].genParton_pt	= jet2_genParton_pt[i];
      jet2[i].jetArea	= jet2_jetArea[i];
      jet2[i].jetBProbabilityBJetTags	= jet2_jetBProbabilityBJetTags[i];
      jet2[i].jetID_fHPD	= jet2_jetID_fHPD[i];
      jet2[i].jetID_n90Hits	= jet2_jetID_n90Hits[i];
      jet2[i].jetProbabilityBJetTags	= jet2_jetProbabilityBJetTags[i];
      jet2[i].muonEnergy	= jet2_muonEnergy[i];
      jet2[i].neutralEmEnergyFraction	= jet2_neutralEmEnergyFraction[i];
      jet2[i].neutralHadronEnergy	= jet2_neutralHadronEnergy[i];
      jet2[i].neutralHadronEnergyFraction	= jet2_neutralHadronEnergyFraction[i];
      jet2[i].neutralHadronMultiplicity	= jet2_neutralHadronMultiplicity[i];
      jet2[i].neutralMultiplicity	= jet2_neutralMultiplicity[i];
      jet2[i].numberOfDaughters	= jet2_numberOfDaughters[i];
      jet2[i].partonFlavour	= jet2_partonFlavour[i];
      jet2[i].phi	= jet2_phi[i];
      jet2[i].photonEnergy	= jet2_photonEnergy[i];
      jet2[i].photonMultiplicity	= jet2_photonMultiplicity[i];
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
      jet2[i].uncor_et	= jet2_uncor_et[i];
      jet2[i].uncor_eta	= jet2_uncor_eta[i];
      jet2[i].uncor_phi	= jet2_uncor_phi[i];
      jet2[i].uncor_pt	= jet2_uncor_pt[i];
    }

  jethelper.resize(jethelper_jecFactor.size());
  for(unsigned int i=0; i < jethelper_jecFactor.size(); ++i)
    {
      jethelper[i].jecFactor	= jethelper_jecFactor[i];
      jethelper[i].jecFactorNoL1Fast	= jethelper_jecFactorNoL1Fast[i];
      jethelper[i].jetUncMinus	= jethelper_jetUncMinus[i];
      jethelper[i].jetUncPlus	= jethelper_jetUncPlus[i];
      jethelper[i].nSecondaryVertices	= jethelper_nSecondaryVertices[i];
      jethelper[i].secondaryVertex0DirectionX	= jethelper_secondaryVertex0DirectionX[i];
      jethelper[i].secondaryVertex0DirectionY	= jethelper_secondaryVertex0DirectionY[i];
      jethelper[i].secondaryVertex0DirectionZ	= jethelper_secondaryVertex0DirectionZ[i];
      jethelper[i].secondaryVertex0FlightDistance	= jethelper_secondaryVertex0FlightDistance[i];
      jethelper[i].secondaryVertex0FlightDistanceSig	= jethelper_secondaryVertex0FlightDistanceSig[i];
      jethelper[i].secondaryVertex0Mass	= jethelper_secondaryVertex0Mass[i];
      jethelper[i].secondaryVertex1DirectionX	= jethelper_secondaryVertex1DirectionX[i];
      jethelper[i].secondaryVertex1DirectionY	= jethelper_secondaryVertex1DirectionY[i];
      jethelper[i].secondaryVertex1DirectionZ	= jethelper_secondaryVertex1DirectionZ[i];
      jethelper[i].secondaryVertex1FlightDistance	= jethelper_secondaryVertex1FlightDistance[i];
      jethelper[i].secondaryVertex1FlightDistanceSig	= jethelper_secondaryVertex1FlightDistanceSig[i];
      jethelper[i].secondaryVertex1Mass	= jethelper_secondaryVertex1Mass[i];
    }

  jethelper1.resize(jethelper1_jecFactor.size());
  for(unsigned int i=0; i < jethelper1_jecFactor.size(); ++i)
    {
      jethelper1[i].jecFactor	= jethelper1_jecFactor[i];
      jethelper1[i].jecFactorNoL1Fast	= jethelper1_jecFactorNoL1Fast[i];
      jethelper1[i].jetUncMinus	= jethelper1_jetUncMinus[i];
      jethelper1[i].jetUncPlus	= jethelper1_jetUncPlus[i];
      jethelper1[i].nSecondaryVertices	= jethelper1_nSecondaryVertices[i];
      jethelper1[i].secondaryVertex0DirectionX	= jethelper1_secondaryVertex0DirectionX[i];
      jethelper1[i].secondaryVertex0DirectionY	= jethelper1_secondaryVertex0DirectionY[i];
      jethelper1[i].secondaryVertex0DirectionZ	= jethelper1_secondaryVertex0DirectionZ[i];
      jethelper1[i].secondaryVertex0FlightDistance	= jethelper1_secondaryVertex0FlightDistance[i];
      jethelper1[i].secondaryVertex0FlightDistanceSig	= jethelper1_secondaryVertex0FlightDistanceSig[i];
      jethelper1[i].secondaryVertex0Mass	= jethelper1_secondaryVertex0Mass[i];
      jethelper1[i].secondaryVertex1DirectionX	= jethelper1_secondaryVertex1DirectionX[i];
      jethelper1[i].secondaryVertex1DirectionY	= jethelper1_secondaryVertex1DirectionY[i];
      jethelper1[i].secondaryVertex1DirectionZ	= jethelper1_secondaryVertex1DirectionZ[i];
      jethelper1[i].secondaryVertex1FlightDistance	= jethelper1_secondaryVertex1FlightDistance[i];
      jethelper1[i].secondaryVertex1FlightDistanceSig	= jethelper1_secondaryVertex1FlightDistanceSig[i];
      jethelper1[i].secondaryVertex1Mass	= jethelper1_secondaryVertex1Mass[i];
    }

  jethelper2.resize(jethelper2_jecFactor.size());
  for(unsigned int i=0; i < jethelper2_jecFactor.size(); ++i)
    {
      jethelper2[i].jecFactor	= jethelper2_jecFactor[i];
      jethelper2[i].jecFactorNoL1Fast	= jethelper2_jecFactorNoL1Fast[i];
      jethelper2[i].jetUncMinus	= jethelper2_jetUncMinus[i];
      jethelper2[i].jetUncPlus	= jethelper2_jetUncPlus[i];
      jethelper2[i].nSecondaryVertices	= jethelper2_nSecondaryVertices[i];
      jethelper2[i].secondaryVertex0DirectionX	= jethelper2_secondaryVertex0DirectionX[i];
      jethelper2[i].secondaryVertex0DirectionY	= jethelper2_secondaryVertex0DirectionY[i];
      jethelper2[i].secondaryVertex0DirectionZ	= jethelper2_secondaryVertex0DirectionZ[i];
      jethelper2[i].secondaryVertex0FlightDistance	= jethelper2_secondaryVertex0FlightDistance[i];
      jethelper2[i].secondaryVertex0FlightDistanceSig	= jethelper2_secondaryVertex0FlightDistanceSig[i];
      jethelper2[i].secondaryVertex0Mass	= jethelper2_secondaryVertex0Mass[i];
      jethelper2[i].secondaryVertex1DirectionX	= jethelper2_secondaryVertex1DirectionX[i];
      jethelper2[i].secondaryVertex1DirectionY	= jethelper2_secondaryVertex1DirectionY[i];
      jethelper2[i].secondaryVertex1DirectionZ	= jethelper2_secondaryVertex1DirectionZ[i];
      jethelper2[i].secondaryVertex1FlightDistance	= jethelper2_secondaryVertex1FlightDistance[i];
      jethelper2[i].secondaryVertex1FlightDistanceSig	= jethelper2_secondaryVertex1FlightDistanceSig[i];
      jethelper2[i].secondaryVertex1Mass	= jethelper2_secondaryVertex1Mass[i];
    }

  met.resize(met_energy.size());
  for(unsigned int i=0; i < met_energy.size(); ++i)
    {
      met[i].energy	= met_energy[i];
      met[i].et	= met_et[i];
      met[i].genMET_energy	= met_genMET_energy[i];
      met[i].genMET_et	= met_genMET_et[i];
      met[i].genMET_phi	= met_genMET_phi[i];
      met[i].genMET_pt	= met_genMET_pt[i];
      met[i].mEtSig	= met_mEtSig[i];
      met[i].metSignificance	= met_metSignificance[i];
      met[i].phi	= met_phi[i];
      met[i].pt	= met_pt[i];
      met[i].significance	= met_significance[i];
      met[i].sumEt	= met_sumEt[i];
    }

  met1.resize(met1_energy.size());
  for(unsigned int i=0; i < met1_energy.size(); ++i)
    {
      met1[i].energy	= met1_energy[i];
      met1[i].et	= met1_et[i];
      met1[i].genMET_energy	= met1_genMET_energy[i];
      met1[i].genMET_et	= met1_genMET_et[i];
      met1[i].genMET_phi	= met1_genMET_phi[i];
      met1[i].genMET_pt	= met1_genMET_pt[i];
      met1[i].mEtSig	= met1_mEtSig[i];
      met1[i].phi	= met1_phi[i];
      met1[i].pt	= met1_pt[i];
      met1[i].significance	= met1_significance[i];
      met1[i].sumEt	= met1_sumEt[i];
    }

  met2.resize(met2_energy.size());
  for(unsigned int i=0; i < met2_energy.size(); ++i)
    {
      met2[i].energy	= met2_energy[i];
      met2[i].et	= met2_et[i];
      met2[i].genMET_energy	= met2_genMET_energy[i];
      met2[i].genMET_et	= met2_genMET_et[i];
      met2[i].genMET_phi	= met2_genMET_phi[i];
      met2[i].genMET_pt	= met2_genMET_pt[i];
      met2[i].mEtSig	= met2_mEtSig[i];
      met2[i].phi	= met2_phi[i];
      met2[i].pt	= met2_pt[i];
      met2[i].significance	= met2_significance[i];
      met2[i].sumEt	= met2_sumEt[i];
    }

  met3.resize(met3_energy.size());
  for(unsigned int i=0; i < met3_energy.size(); ++i)
    {
      met3[i].energy	= met3_energy[i];
      met3[i].et	= met3_et[i];
      met3[i].genMET_energy	= met3_genMET_energy[i];
      met3[i].genMET_et	= met3_genMET_et[i];
      met3[i].genMET_phi	= met3_genMET_phi[i];
      met3[i].genMET_pt	= met3_genMET_pt[i];
      met3[i].mEtSig	= met3_mEtSig[i];
      met3[i].phi	= met3_phi[i];
      met3[i].pt	= met3_pt[i];
      met3[i].significance	= met3_significance[i];
      met3[i].sumEt	= met3_sumEt[i];
    }

  met4.resize(met4_energy.size());
  for(unsigned int i=0; i < met4_energy.size(); ++i)
    {
      met4[i].energy	= met4_energy[i];
      met4[i].et	= met4_et[i];
      met4[i].genMET_energy	= met4_genMET_energy[i];
      met4[i].genMET_et	= met4_genMET_et[i];
      met4[i].genMET_phi	= met4_genMET_phi[i];
      met4[i].genMET_pt	= met4_genMET_pt[i];
      met4[i].mEtSig	= met4_mEtSig[i];
      met4[i].phi	= met4_phi[i];
      met4[i].pt	= met4_pt[i];
      met4[i].significance	= met4_significance[i];
      met4[i].sumEt	= met4_sumEt[i];
    }

  met5.resize(met5_energy.size());
  for(unsigned int i=0; i < met5_energy.size(); ++i)
    {
      met5[i].energy	= met5_energy[i];
      met5[i].et	= met5_et[i];
      met5[i].genMET_energy	= met5_genMET_energy[i];
      met5[i].genMET_et	= met5_genMET_et[i];
      met5[i].genMET_phi	= met5_genMET_phi[i];
      met5[i].genMET_pt	= met5_genMET_pt[i];
      met5[i].mEtSig	= met5_mEtSig[i];
      met5[i].phi	= met5_phi[i];
      met5[i].pt	= met5_pt[i];
      met5[i].significance	= met5_significance[i];
      met5[i].sumEt	= met5_sumEt[i];
    }

  met6.resize(met6_energy.size());
  for(unsigned int i=0; i < met6_energy.size(); ++i)
    {
      met6[i].energy	= met6_energy[i];
      met6[i].et	= met6_et[i];
      met6[i].genMET_energy	= met6_genMET_energy[i];
      met6[i].genMET_et	= met6_genMET_et[i];
      met6[i].genMET_phi	= met6_genMET_phi[i];
      met6[i].genMET_pt	= met6_genMET_pt[i];
      met6[i].mEtSig	= met6_mEtSig[i];
      met6[i].phi	= met6_phi[i];
      met6[i].pt	= met6_pt[i];
      met6[i].significance	= met6_significance[i];
      met6[i].sumEt	= met6_sumEt[i];
    }

  met7.resize(met7_energy.size());
  for(unsigned int i=0; i < met7_energy.size(); ++i)
    {
      met7[i].energy	= met7_energy[i];
      met7[i].et	= met7_et[i];
      met7[i].genMET_energy	= met7_genMET_energy[i];
      met7[i].genMET_et	= met7_genMET_et[i];
      met7[i].genMET_phi	= met7_genMET_phi[i];
      met7[i].genMET_pt	= met7_genMET_pt[i];
      met7[i].mEtSig	= met7_mEtSig[i];
      met7[i].phi	= met7_phi[i];
      met7[i].pt	= met7_pt[i];
      met7[i].significance	= met7_significance[i];
      met7[i].sumEt	= met7_sumEt[i];
    }

  methelper.resize(methelper_significance_dxx.size());
  for(unsigned int i=0; i < methelper_significance_dxx.size(); ++i)
    {
      methelper[i].significance_dxx	= methelper_significance_dxx[i];
      methelper[i].significance_dxy	= methelper_significance_dxy[i];
      methelper[i].significance_dyx	= methelper_significance_dyx[i];
      methelper[i].significance_dyy	= methelper_significance_dyy[i];
    }

  methelper1.resize(methelper1_significance_dxx.size());
  for(unsigned int i=0; i < methelper1_significance_dxx.size(); ++i)
    {
      methelper1[i].significance_dxx	= methelper1_significance_dxx[i];
      methelper1[i].significance_dxy	= methelper1_significance_dxy[i];
      methelper1[i].significance_dyx	= methelper1_significance_dyx[i];
      methelper1[i].significance_dyy	= methelper1_significance_dyy[i];
    }

  muon.resize(muon_AllGlobalMuons.size());
  for(unsigned int i=0; i < muon_AllGlobalMuons.size(); ++i)
    {
      muon[i].AllGlobalMuons	= muon_AllGlobalMuons[i];
      muon[i].AllTrackerMuons	= muon_AllTrackerMuons[i];
      muon[i].GlobalMuonPromptTight	= muon_GlobalMuonPromptTight[i];
      muon[i].charge	= muon_charge[i];
      muon[i].chargedHadronIso	= muon_chargedHadronIso[i];
      muon[i].combinedMuon_chi2	= muon_combinedMuon_chi2[i];
      muon[i].combinedMuon_ndof	= muon_combinedMuon_ndof[i];
      muon[i].combinedMuon_numberOfValidHits	= muon_combinedMuon_numberOfValidHits[i];
      muon[i].dB	= muon_dB[i];
      muon[i].ecalIso	= muon_ecalIso[i];
      muon[i].energy	= muon_energy[i];
      muon[i].et	= muon_et[i];
      muon[i].eta	= muon_eta[i];
      muon[i].genLepton_eta	= muon_genLepton_eta[i];
      muon[i].genLepton_pdgId	= muon_genLepton_pdgId[i];
      muon[i].genLepton_phi	= muon_genLepton_phi[i];
      muon[i].genLepton_pt	= muon_genLepton_pt[i];
      muon[i].globalTrack_chi2	= muon_globalTrack_chi2[i];
      muon[i].globalTrack_hitPattern_numberOfValidTrackerHits	= muon_globalTrack_hitPattern_numberOfValidTrackerHits[i];
      muon[i].globalTrack_ndof	= muon_globalTrack_ndof[i];
      muon[i].globalTrack_pt	= muon_globalTrack_pt[i];
      muon[i].globalTrack_ptError	= muon_globalTrack_ptError[i];
      muon[i].hcalIso	= muon_hcalIso[i];
      muon[i].innerTrack_d0	= muon_innerTrack_d0[i];
      muon[i].innerTrack_hitPattern_pixelLayersWithMeasurement	= muon_innerTrack_hitPattern_pixelLayersWithMeasurement[i];
      muon[i].innerTrack_numberOfValidHits	= muon_innerTrack_numberOfValidHits[i];
      muon[i].innerTrack_phi	= muon_innerTrack_phi[i];
      muon[i].innerTrack_pt	= muon_innerTrack_pt[i];
      muon[i].innerTrack_vertex_z	= muon_innerTrack_vertex_z[i];
      muon[i].isGlobalMuon	= muon_isGlobalMuon[i];
      muon[i].isTrackerMuon	= muon_isTrackerMuon[i];
      muon[i].neutralHadronIso	= muon_neutralHadronIso[i];
      muon[i].numberOfMatches	= muon_numberOfMatches[i];
      muon[i].phi	= muon_phi[i];
      muon[i].photonIso	= muon_photonIso[i];
      muon[i].pt	= muon_pt[i];
      muon[i].trackIso	= muon_trackIso[i];
      muon[i].track_d0	= muon_track_d0[i];
      muon[i].track_hitPattern_numberOfValidMuonHits	= muon_track_hitPattern_numberOfValidMuonHits[i];
      muon[i].track_hitPattern_numberOfValidPixelHits	= muon_track_hitPattern_numberOfValidPixelHits[i];
      muon[i].track_normalizedChi2	= muon_track_normalizedChi2[i];
      muon[i].vz	= muon_vz[i];
    }

  muon1.resize(muon1_AllGlobalMuons.size());
  for(unsigned int i=0; i < muon1_AllGlobalMuons.size(); ++i)
    {
      muon1[i].AllGlobalMuons	= muon1_AllGlobalMuons[i];
      muon1[i].AllTrackerMuons	= muon1_AllTrackerMuons[i];
      muon1[i].GlobalMuonPromptTight	= muon1_GlobalMuonPromptTight[i];
      muon1[i].charge	= muon1_charge[i];
      muon1[i].chargedHadronIso	= muon1_chargedHadronIso[i];
      muon1[i].combinedMuon_chi2	= muon1_combinedMuon_chi2[i];
      muon1[i].combinedMuon_ndof	= muon1_combinedMuon_ndof[i];
      muon1[i].combinedMuon_numberOfValidHits	= muon1_combinedMuon_numberOfValidHits[i];
      muon1[i].dB	= muon1_dB[i];
      muon1[i].ecalIso	= muon1_ecalIso[i];
      muon1[i].energy	= muon1_energy[i];
      muon1[i].et	= muon1_et[i];
      muon1[i].eta	= muon1_eta[i];
      muon1[i].genLepton_eta	= muon1_genLepton_eta[i];
      muon1[i].genLepton_pdgId	= muon1_genLepton_pdgId[i];
      muon1[i].genLepton_phi	= muon1_genLepton_phi[i];
      muon1[i].genLepton_pt	= muon1_genLepton_pt[i];
      muon1[i].globalTrack_chi2	= muon1_globalTrack_chi2[i];
      muon1[i].globalTrack_hitPattern_numberOfValidTrackerHits	= muon1_globalTrack_hitPattern_numberOfValidTrackerHits[i];
      muon1[i].globalTrack_ndof	= muon1_globalTrack_ndof[i];
      muon1[i].globalTrack_pt	= muon1_globalTrack_pt[i];
      muon1[i].globalTrack_ptError	= muon1_globalTrack_ptError[i];
      muon1[i].hcalIso	= muon1_hcalIso[i];
      muon1[i].innerTrack_d0	= muon1_innerTrack_d0[i];
      muon1[i].innerTrack_hitPattern_pixelLayersWithMeasurement	= muon1_innerTrack_hitPattern_pixelLayersWithMeasurement[i];
      muon1[i].innerTrack_numberOfValidHits	= muon1_innerTrack_numberOfValidHits[i];
      muon1[i].innerTrack_phi	= muon1_innerTrack_phi[i];
      muon1[i].innerTrack_pt	= muon1_innerTrack_pt[i];
      muon1[i].innerTrack_vertex_z	= muon1_innerTrack_vertex_z[i];
      muon1[i].isGlobalMuon	= muon1_isGlobalMuon[i];
      muon1[i].isTrackerMuon	= muon1_isTrackerMuon[i];
      muon1[i].neutralHadronIso	= muon1_neutralHadronIso[i];
      muon1[i].numberOfMatches	= muon1_numberOfMatches[i];
      muon1[i].phi	= muon1_phi[i];
      muon1[i].photonIso	= muon1_photonIso[i];
      muon1[i].pt	= muon1_pt[i];
      muon1[i].trackIso	= muon1_trackIso[i];
      muon1[i].track_d0	= muon1_track_d0[i];
      muon1[i].track_hitPattern_numberOfValidMuonHits	= muon1_track_hitPattern_numberOfValidMuonHits[i];
      muon1[i].track_hitPattern_numberOfValidPixelHits	= muon1_track_hitPattern_numberOfValidPixelHits[i];
      muon1[i].track_normalizedChi2	= muon1_track_normalizedChi2[i];
      muon1[i].vz	= muon1_vz[i];
    }

  muonhelper.resize(muonhelper_dxywrtBeamSpot.size());
  for(unsigned int i=0; i < muonhelper_dxywrtBeamSpot.size(); ++i)
    {
      muonhelper[i].dxywrtBeamSpot	= muonhelper_dxywrtBeamSpot[i];
      muonhelper[i].muonPFcandptdiff	= muonhelper_muonPFcandptdiff[i];
    }

  muonhelper1.resize(muonhelper1_dxywrtBeamSpot.size());
  for(unsigned int i=0; i < muonhelper1_dxywrtBeamSpot.size(); ++i)
    {
      muonhelper1[i].dxywrtBeamSpot	= muonhelper1_dxywrtBeamSpot[i];
      muonhelper1[i].muonPFcandptdiff	= muonhelper1_muonPFcandptdiff[i];
    }

  pileupsummaryinfo.resize(pileupsummaryinfo_addpileupinfo_getBunchCrossing.size());
  for(unsigned int i=0; i < pileupsummaryinfo_addpileupinfo_getBunchCrossing.size(); ++i)
    {
      pileupsummaryinfo[i].addpileupinfo_getBunchCrossing	= pileupsummaryinfo_addpileupinfo_getBunchCrossing[i];
      pileupsummaryinfo[i].addpileupinfo_getPU_NumInteractions	= pileupsummaryinfo_addpileupinfo_getPU_NumInteractions[i];
    }

  tau.resize(tau_caloIso.size());
  for(unsigned int i=0; i < tau_caloIso.size(); ++i)
    {
      tau[i].caloIso	= tau_caloIso[i];
      tau[i].ecalIso	= tau_ecalIso[i];
      tau[i].energy	= tau_energy[i];
      tau[i].et	= tau_et[i];
      tau[i].eta	= tau_eta[i];
      tau[i].hcalIso	= tau_hcalIso[i];
      tau[i].phi	= tau_phi[i];
      tau[i].pt	= tau_pt[i];
      tau[i].tauID_againstElectron	= tau_tauID_againstElectron[i];
      tau[i].tauID_againstMuon	= tau_tauID_againstMuon[i];
      tau[i].tauID_byIsolation	= tau_tauID_byIsolation[i];
      tau[i].tauID_byTaNC	= tau_tauID_byTaNC[i];
      tau[i].tauID_byTaNCfrHalfPercent	= tau_tauID_byTaNCfrHalfPercent[i];
      tau[i].tauID_byTaNCfrQuarterPercent	= tau_tauID_byTaNCfrQuarterPercent[i];
      tau[i].trackIso	= tau_trackIso[i];
    }

  tauhelper.resize(tauhelper_genTauDecayModeID.size());
  for(unsigned int i=0; i < tauhelper_genTauDecayModeID.size(); ++i)
    {
      tauhelper[i].genTauDecayModeID	= tauhelper_genTauDecayModeID[i];
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
