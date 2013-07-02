
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TMath.h"

#include <iostream>

//
// Note: If you are using this on a SLC6 machine, you need to do
//       this in your interactive root session
//
//          gSystem->AddIncludePath(" -D__USE_XOPEN2K8 ")
//
//       before doing this
//
//          .L doSkimSlim.C+
//

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
         inReducedTree -> SetBranchStatus("njets30",1) ;

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
         inReducedTree -> SetBranchStatus("PUweight",1) ;
         inReducedTree -> SetBranchStatus("nGoodPV",1) ;

        //--- jet pt
         inReducedTree -> SetBranchStatus("jetpt1",1) ;
         inReducedTree -> SetBranchStatus("jetpt2",1) ;
         inReducedTree -> SetBranchStatus("jetpt3",1) ;
         inReducedTree -> SetBranchStatus("jetpt4",1) ;

        //--- new vars from Josh for trying to understand ABCD non-closure.
         inReducedTree -> SetBranchStatus("higgs1jetpt1",1) ;
         inReducedTree -> SetBranchStatus("higgs1jetpt2",1) ;
         inReducedTree -> SetBranchStatus("higgs2jetpt1",1) ;
         inReducedTree -> SetBranchStatus("higgs2jetpt2",1) ;
         inReducedTree -> SetBranchStatus("higgs1CSV1",1) ;
         inReducedTree -> SetBranchStatus("higgs1CSV2",1) ;
         inReducedTree -> SetBranchStatus("higgs2CSV1",1) ;
         inReducedTree -> SetBranchStatus("higgs2CSV2",1) ;
         inReducedTree -> SetBranchStatus("higgs1partonId1",1) ;
         inReducedTree -> SetBranchStatus("higgs1partonId2",1) ;
         inReducedTree -> SetBranchStatus("higgs2partonId1",1) ;
         inReducedTree -> SetBranchStatus("higgs2partonId2",1) ;
         inReducedTree -> SetBranchStatus("higgs1partonMomId1",1) ;
         inReducedTree -> SetBranchStatus("higgs1partonMomId2",1) ;
         inReducedTree -> SetBranchStatus("higgs2partonMomId1",1) ;
         inReducedTree -> SetBranchStatus("higgs2partonMomId2",1) ;
         inReducedTree -> SetBranchStatus("nGluonsSplitToLF",1) ;
         inReducedTree -> SetBranchStatus("nGluonsSplitToC",1) ;
         inReducedTree -> SetBranchStatus("nGluonsSplitToB",1) ;



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
      int njets30 ;
      inReducedTree -> SetBranchAddress("njets30", &njets30 ) ;

      float CSVbest2 ;
      inReducedTree -> SetBranchAddress("CSVbest2", &CSVbest2) ;

      bool cutPV, passCleaning, buggyEvent ;
      inReducedTree -> SetBranchAddress("cutPV", &cutPV ) ;
      inReducedTree -> SetBranchAddress("passCleaning", &passCleaning ) ;
      inReducedTree -> SetBranchAddress("buggyEvent", &buggyEvent ) ;




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
      float projected_remaining(999999.) ;
      for ( Long64_t ievt=0; ievt<nentries; ievt++ ) {

         if ( ievt%1000 == 0 ) {
            int thistime = sw.RealTime() ;
            sw.Continue() ;
            if ( thistime < 2 ) {
               printf("   %10llu out of %10llu  (%6.1f%%) \r", ievt, nentries, 100.*ievt/(1.*nentries) ) ;
            } else {
               if ( thistime > time ) projected_remaining = (1.*thistime)/(1.*ievt)*(nentries-ievt) ;
               if ( projected_remaining < 100 ) {
                  printf("   %10llu out of %10llu  (%6.1f%%)    seconds remaining %4.0f                       \r", ievt, nentries, 100.*ievt/(1.*nentries), projected_remaining ) ;
               } else if ( projected_remaining < 3600 ) {
                  printf("   %10llu out of %10llu  (%6.1f%%)    time remaining     %2d:%02d   \r", ievt, nentries, 100.*ievt/(1.*nentries),
                       TMath::Nint(projected_remaining)/60, TMath::Nint(projected_remaining)%60 ) ;
               } else {
                  printf("   %10llu out of %10llu  (%6.1f%%)    time remaining  %2d:%02d:%02d   \r", ievt, nentries, 100.*ievt/(1.*nentries),
                       TMath::Nint(projected_remaining)/3600, (TMath::Nint(projected_remaining)%3600)/60, TMath::Nint(projected_remaining)%60 ) ;
               }
            }
            cout << flush ;
            time = thistime ;
         }

         inReducedTree -> GetEntry(ievt) ;

         if ( !cutPV ) continue ;
         if ( !passCleaning ) continue ;
         if ( buggyEvent ) continue ;
         if ( !(trig1 || trig2 || trig3) ) continue ;
         //-----------
         // Owen: make njets cut safe for pt>20 or pt>30.
         //       this means you must cut on the appropriate njets variable when using the skim output.
         if ( njets20<4 || njets30>5 ) continue ;
         //-----------
         if ( CSVbest2 < 0.898 ) continue ;

         outReducedTree->Fill() ;

      } // ievt.

      printf("\n\n\n Done.\n\n\n") ;

      printf("\n\n Output file:  %s\n\n\n", outfile_name.Data() ) ;

      if ( time > 3600 ) {
         printf( "   Total time:  %2d:%02d:%02d \n\n\n", time/3600, (time%3600)/60, time%60 ) ;
      } else if ( time > 100 ) {
         printf( "   Total time:     %02d:%02d \n\n\n", time/60, time%60 ) ;
      } else {
         printf( "   Total time:     %d seconds \n\n\n", time ) ;
      }

      outReducedTree->AutoSave() ;

      delete outfile ;

  }

