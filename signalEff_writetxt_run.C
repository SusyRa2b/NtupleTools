void signalEff_writetxt_run(TString what, TString sample) {

  gSystem->Load("CrossSectionTable_cxx.so");

  gSystem->Load("signalEff2012_writetxt_C.so");

  TString prefix = "eventcounts.";
  TString prefix2="eventcounts.mergebbins.";

  if (what=="counts" || what=="PU") writetxt(what,sample,prefix,true,false);
  else {
    writetxt(what,sample,prefix2,true,false);
  }

}
