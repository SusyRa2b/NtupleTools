
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TMath.h"

#include <iostream>


  void doSkimSlim( const char* infile_name, bool doSlim = false  ) {

      TFile* infile = new TFile( infile_name ) ;
      if ( ! (infile->IsOpen()) ) return ;

      TTree* inReducedTree = (TTree*) infile->Get("reducedTree") ;

      Long64_t nentries = inReducedTree -> GetEntries() ;

      printf("\n\n Number of entries: %llu\n\n", nentries ) ;

      if ( doSlim ) {

         inReducedTree -> SetBranchStatus("*",0) ; // disable all branches.

        //--- Now enable only the ones we want to save in the output file.
        //    Add anything you care about here.

        //--- basic cut vars
         inReducedTree -> SetBranchStatus("cutPV",1) ;
         inReducedTree -> SetBranchStatus("passCleaning",1) ;
         inReducedTree -> SetBranchStatus("buggyEvent",1) ;

        //--- trigger
         inReducedTree -> SetBranchStatus("passMC_DiCentralPFJet30_PFMET80_BTagCSV07",1) ;
         inReducedTree -> SetBranchStatus("passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05",1) ;
         inReducedTree -> SetBranchStatus("passMC_PFMET150",1) ;

        //--- lepton veto vars
         inReducedTree -> SetBranchStatus("nMuons",1) ;
         inReducedTree -> SetBranchStatus("nElectrons",1) ;
         inReducedTree -> SetBranchStatus("nIsoTracks15_005_03",1) ;
         inReducedTree -> SetBranchStatus("nIsoTracks5_005_03",1) ;
         inReducedTree -> SetBranchStatus("nTausLoose",1) ;

        //--- njets
         inReducedTree -> SetBranchStatus("njets20",1) ;

        //--- btagging
         inReducedTree -> SetBranchStatus("CSVbest1",1) ;
         inReducedTree -> SetBranchStatus("CSVbest2",1) ;
         inReducedTree -> SetBranchStatus("CSVbest3",1) ;
         inReducedTree -> SetBranchStatus("CSVbest4",1) ;

        //--- Essential Higgs vars
         inReducedTree -> SetBranchStatus("higgsMbb1MassDiff",1) ;
         inReducedTree -> SetBranchStatus("higgsMbb2MassDiff",1) ;
         inReducedTree -> SetBranchStatus("higgsMbb1MassDiff_correct",1) ;
         inReducedTree -> SetBranchStatus("maxDeltaPhi_bb_bb",1) ;
         inReducedTree -> SetBranchStatus("deltaRmax_hh",1) ;
         inReducedTree -> SetBranchStatus("deltaRmin_hh",1) ;

        //--- Kinematic vars
         inReducedTree -> SetBranchStatus("MET",1) ;
         inReducedTree -> SetBranchStatus("METsig",1) ;

        //--- Event weights
         inReducedTree -> SetBranchStatus("weight",1) ;
         inReducedTree -> SetBranchStatus("weight2",1) ;
         inReducedTree -> SetBranchStatus("weight3",1) ;

      } else {

         inReducedTree -> SetBranchStatus("*",1) ; // enable all branches.

      }






     //--- Vars needed to decide whether or not to save the event.
      bool trig1, trig2, trig3 ;
      inReducedTree -> SetBranchAddress("passMC_DiCentralPFJet30_PFMET80_BTagCSV07", &trig1 ) ;
      inReducedTree -> SetBranchAddress("passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05", &trig2 ) ;
      inReducedTree -> SetBranchAddress("passMC_PFMET150", &trig3 ) ;

      int njets20 ;
      inReducedTree -> SetBranchAddress("njets20", &njets20 ) ;

      float CSVbest2 ;
      inReducedTree -> SetBranchAddress("CSVbest2", &CSVbest2) ;





     //--- Open output file
      TString outfile_name( infile_name ) ;
      if ( doSlim ) {
         if ( outfile_name.Contains("-skim.root") ) {
            outfile_name.ReplaceAll( "-skim.root", "-slimskim.root" ) ;
         } else {
            outfile_name.ReplaceAll( ".root", "-slimskim.root" ) ;
         }
      } else {
         outfile_name.ReplaceAll( ".root", "-skim.root" ) ;
      }
      if ( outfile_name.CompareTo( infile_name ) == 0 ) {
         printf("\n\n *** Input and output file names same.  Input doesn't contain .root in name?\n") ;
         printf("    input: %s \n", infile_name ) ;
         printf("    output: %s \n", outfile_name.Data() ) ;
         return ;
      }
      printf("\n\n Output file: %s\n\n", outfile_name.Data() ) ;
      char command[10000] ;
      sprintf( command, "ls %s >& /dev/null", outfile_name.Data() ) ;
      int returnstat = gSystem->Exec( command ) ;
      if ( returnstat == 0 ) {
         char mvfile[10000] ;
         sprintf( mvfile, "%s-old", outfile_name.Data() ) ;
         printf("\n\n *** Output file already exists.  Moving it to %s\n\n", mvfile ) ;
         sprintf( command, "mv %s %s", outfile_name.Data(), mvfile ) ;
         gSystem->Exec( command ) ;
      }
      TFile* outfile = new TFile( outfile_name, "recreate" ) ;
      TTree* outReducedTree = inReducedTree->CloneTree(0) ;







     //--- Loop over the events.
      TStopwatch sw ;
      sw.Start() ;
      int time(0) ;
      for ( Long64_t ievt=0; ievt<nentries; ievt++ ) {

         if ( ievt%1000 == 0 ) {
            int thistime = sw.RealTime() ;
            sw.Continue() ;
            float projected_remaining ;
            if ( thistime < 2 ) {
               printf("   %10llu out of %10llu  (%6.1f%%) \r", ievt, nentries, 100.*ievt/(1.*nentries) ) ;
            } else {
               if ( thistime > time ) projected_remaining = (1.*thistime)/(1.*ievt)*(nentries-ievt) ;
               if ( projected_remaining < 100 ) {
                  printf("   %10llu out of %10llu  (%6.1f%%)    seconds remaining %4.0f                       \r", ievt, nentries, 100.*ievt/(1.*nentries), projected_remaining ) ;
               } else if ( projected_remaining < 3600 ) {
                  printf("   %10llu out of %10llu  (%6.1f%%)    time remaining  %2d:%02d   \r", ievt, nentries, 100.*ievt/(1.*nentries), TMath::Nint(projected_remaining)/60, TMath::Nint(projected_remaining)%60 ) ;
               }
            }
            cout << flush ;
            time = thistime ;
         }

         inReducedTree -> GetEntry(ievt) ;

         if ( !(trig1 || trig2 || trig3) ) continue ;
         if ( njets20<4 || njets20>5 ) continue ;
         if ( CSVbest2 < 0.898 ) continue ;

         outReducedTree->Fill() ;

      } // ievt.

      printf("\n\n\n Done.\n\n\n") ;

      printf("\n\n Output file:  %s\n\n\n", outfile_name.Data() ) ;

      outReducedTree->AutoSave() ;

      delete outfile ;

  }

