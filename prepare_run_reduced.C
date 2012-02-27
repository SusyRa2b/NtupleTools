#include "TString.h"
#include "TChain.h"
#include "TObjArray.h"
#include <fstream> 
#include <iostream>

using namespace std;

void prepare_run_reduced(){
  //specify where the ntuples are
  
  //NEWER 42X data and MC
  TString dir1 = "/cu2/ra2b/Ftuples/";
  TString inputVersion = "V00-02-35_example"; 

  //NEW 42X data and MC - don't use for SUS-11-006!!
  //TString dir1 = "/cu2/ra2b/Ftuples/";
  //TString inputVersion = "V00-02-24";

  ////42X data and MC
  //TString dir1 = "/cu1/kreis/Ftuples/";
  //TString inputVersion = "V00-02-05";

  ////41X MC
  //TString dir1 = "/cu3/kreis/Ftuples/";
  //TString inputVersion = "V00-01-03";

  ofstream runstream;
  runstream.open("run_reducedTrees.sh", ios::trunc);
  
  dir1+=inputVersion;
  dir1+="/*";
  TChain dummy1("dummy1");
  dummy1.Add(dir1);
  TObjArray* dir1list = dummy1.GetListOfFiles();
  int nfiles1=dir1list->GetEntries();
  
  for (int ifile1=0; ifile1<nfiles1; ifile1++) {
    TString samplefiles1 = dir1list->At(ifile1)->GetTitle();
    TString sampleName =  samplefiles1( samplefiles1.Last('/')+1, samplefiles1.Length() );
    //if(!sampleName.Contains("QCD")) continue;

    cout << sampleName << endl;
    ofstream fileliststream;
    fileliststream.open(sampleName+".txt", ios::trunc);
    
    TString dir2 = samplefiles1;
    dir2+="/*.root";
    TChain dummy2("dummy3");
    dummy2.Add(dir2);
    TObjArray* dir2list = dummy2.GetListOfFiles();
    int nfiles2=dir2list->GetEntries();
    for (int ifile2=0; ifile2<nfiles2; ifile2++) {
      TString samplefiles2 = dir2list->At(ifile2)->GetTitle();
      fileliststream << samplefiles2 << endl;
    }

    fileliststream.close();
    runstream << "./BasicLoopCU " << sampleName+".txt"<< endl;

  }
  runstream.close();
}
