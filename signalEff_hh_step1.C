void signalEff_hh_step1(TString file,TString options) {

  gSystem->Load("signalEff_hbbhbb_C.so");

  //split off the path
  TString path = file(0,file.Index("reducedTree."));

  //chop off .root
  TString stub= file(file.Index("reducedTree.")+12,file.Length()-file.Index("reducedTree.")-5-12);

  bool joinb = options.Contains("joinbtag");
  bool usebtagsf = !joinb; //can't use btagsf if we're joining b-tag bins
  bool dopdfs=options.Contains("pdfs");
  bool dopu=options.Contains("pushift");

  IsrMode isrmode = kIsr0;
  if      ( options.Contains("isrup"))   isrmode = kIsrUp;
  else if ( options.Contains("isrdown")) isrmode = kIsrDown;

  BTagMode btagmode = kBTag0;
  if (options.Contains("btagup")) btagmode = kBTagHfUp;
  else if (options.Contains("btagdown"))  btagmode = kBTagHfDown;

  signalEff_hbbhbb looper(path,stub,joinb,usebtagsf,dopdfs,dopu,isrmode,btagmode);
  looper.Loop();


}
