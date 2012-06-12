
namespace pu {

  //See this twiki for more info:
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupMCReweightingUtilities

  /*
  this data distribution is an hadd of the following histograms:  
  (all in /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/)

  Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.pileupTruth_v2.root
  Cert_165088-167913_7TeV_PromptReco_JSON.pileupTruth_v2.root
  Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v2.pileupTruth_v2.root
  Cert_172620-173692_PromptReco_JSON.pileupTruth_v2.root
  Cert_175832-177515_PromptReco_JSON.pileupTruth_v2.root
  Cert_177718_178078_7TeV_PromptReco_Collisons11_JSON.pileupTruth_v2.root
  Cert_178098-180252_7TeV_PromptReco_Collisions11_JSON.pileupTruth.root
  */

  float TrueDist2011_f[35] = {
    0,
    304761,
    5.78132e+06,
    5.06699e+07,
    2.69238e+08,
    5.18854e+08,
    5.76628e+08,
    5.33043e+08,
    4.66102e+08,
    4.26606e+08,
    3.82425e+08,
    3.51807e+08,
    3.30761e+08,
    2.94186e+08,
    2.29225e+08,
    1.45499e+08,
    7.43711e+07,
    3.08969e+07,
    1.0913e+07,
    3.68176e+06,
    1.13003e+06,
    288668,
    64020.4,
    2350.83,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  };

  /*
  this data distribution is an hadd of the following histograms:  
  (all in /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/)

  Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.pileup_v2.root
  Cert_165088-167913_7TeV_PromptReco_JSON.pileup_v2.root
  Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v2.pileup_v2.root
  Cert_172620-173692_PromptReco_JSON.pileup_v2.root
  Cert_175832-177515_PromptReco_JSON.pileup_v2.root
  Cert_177718_178078_7TeV_PromptReco_Collisons11_JSON.pileup_v2.root
  Cert_178098-180252_7TeV_PromptReco_Collisions11_JSON.pileup_v2.root
  */

  float ObsDist2011_f[35] = {
    1.34465e+07,
    5.90653e+07,
    1.40903e+08,
    2.41301e+08,
    3.33745e+08,
    3.98711e+08,
    4.30106e+08,
    4.32283e+08,
    4.1382e+08,
    3.82846e+08,
    3.45164e+08,
    3.04344e+08,
    2.62555e+08,
    2.21331e+08,
    1.81983e+08,
    1.4569e+08,
    1.13413e+08,
    8.57789e+07,
    6.30124e+07,
    4.49596e+07,
    3.1169e+07,
    2.10079e+07,
    1.37759e+07,
    8.79641e+06,
    5.47442e+06,
    3.32378e+06,
    1.97064e+06,
    1.14204e+06,
    647539,
    359547,
    195673,
    104460,
    54745.2,
    28185.6,
    28005.5
  };


  // Flat10+Tail distribution taken directly from MixingModule input:  
  //(Can be used for Spring11 and Summer11 if you don't worry about small shifts in the mean) 
  //SHOULD be used for 3-D Reweighting, as this is the "true" input for all Summer11 samples.

  Double_t probdistFlat10_f[35] = {
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0630151648,
    0.0526654164,
    0.0402754482,
    0.0292988928,
    0.0194384503,
    0.0122016783,
    0.007207042,
    0.004003637,
    0.0020278322,
    0.0010739954,
    0.0004595759,
    0.0002229748,
    0.0001028162,
    4.58337152809607E-05,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  };


  // Summer11 PU_S4, distribution obtained by only looking at the in-time crossing.  This is the "spike+smear" distribution
  // RECOMMENDED FOR REWEIGHTING (if ignoring out-of-time PU)
  float PoissonOneXDist_f[35] = {
    1.45346E-01,
    6.42802E-02,
    6.95255E-02,
    6.96747E-02,
    6.92955E-02,
    6.84997E-02,
    6.69528E-02,
    6.45515E-02,
    6.09865E-02,
    5.63323E-02,
    5.07322E-02,
    4.44681E-02,
    3.79205E-02,
    3.15131E-02,
    2.54220E-02,
    2.00184E-02,
    1.53776E-02,
    1.15387E-02,
    8.47608E-03,
    6.08715E-03,
    4.28255E-03,
    2.97185E-03,
    2.01918E-03,
    1.34490E-03,
    8.81587E-04,
    5.69954E-04,
    3.61493E-04,
    2.28692E-04,
    1.40791E-04,
    8.44606E-05,
    5.10204E-05,
    3.07802E-05,
    1.81401E-05,
    1.00201E-05,
    5.80004E-06
  };


// Distribution used for Summer2012 MC.
//from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupMCReweightingUtilities
  float Summer2012[60] = {
    2.344E-05,
    2.344E-05,
    2.344E-05,
    2.344E-05,
    4.687E-04,
    4.687E-04,
    7.032E-04,
    9.414E-04,
    1.234E-03,
    1.603E-03,
    2.464E-03,
    3.250E-03,
    5.021E-03,
    6.644E-03,
    8.502E-03,
    1.121E-02,
    1.518E-02,
    2.033E-02,
    2.608E-02,
    3.171E-02,
    3.667E-02,
    4.060E-02,
    4.338E-02,
    4.520E-02,
    4.641E-02,
    4.735E-02,
    4.816E-02,
    4.881E-02,
    4.917E-02,
    4.909E-02,
    4.842E-02,
    4.707E-02,
    4.501E-02,
    4.228E-02,
    3.896E-02,
    3.521E-02,
    3.118E-02,
    2.702E-02,
    2.287E-02,
    1.885E-02,
    1.508E-02,
    1.166E-02,
    8.673E-03,
    6.190E-03,
    4.222E-03,
    2.746E-03,
    1.698E-03,
    9.971E-04,
    5.549E-04,
    2.924E-04,
    1.457E-04,
    6.864E-05,
    3.054E-05,
    1.282E-05,
    5.081E-06,
    1.898E-06,
    6.688E-07,
    2.221E-07,
    6.947E-08,
    2.047E-08
   };  

  //https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
  //very crude first attempt...just trying to get something
  //at command line:
/*
pileupCalc.py -i Cert_190456-194479_8TeV_PromptReco_Collisions12_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-194533.txt --calcMode true --minBiasXsec 69400 --maxPileupBin 60 --numPileupBins 60 MyDataPileupHistogram.root

note that the json files don't exactly correspond
Then in ROOT:
 for (int i=1; i<=60; i++) cout<<pileup->GetBinContent(i)<<","<<endl;
Then copy/paste here

Why the bump at low PU?
Why the large number of overflows in bin 60?
*/

  float FirstData2012[60] = {
    568.943,
    1358.19,
    9475.99,
    2.79881e+06,
    7.33073e+06,
    133124,
    318824,
    3.61796e+06,
    1.48302e+07,
    3.16554e+07,
    5.15906e+07,
    7.65019e+07,
    1.05187e+08,
    1.30138e+08,
    1.51041e+08,
    1.68285e+08,
    1.70619e+08,
    1.52533e+08,
    1.23968e+08,
    9.62378e+07,
    7.39321e+07,
    5.72667e+07,
    4.48809e+07,
    3.5617e+07,
    2.84809e+07,
    2.27126e+07,
    1.79169e+07,
    1.39221e+07,
    1.0631e+07,
    7.96323e+06,
    5.84208e+06,
    4.19264e+06,
    2.94107e+06,
    2.01569e+06,
    1.34937e+06,
    882181,
    563154,
    350949,
    213448,
    126660,
    73306.8,
    41369.1,
    22757.1,
    12200.3,
    6373.27,
    3243.71,
    1608.35,
    776.896,
    365.588,
    167.601,
    74.8564,
    32.5737,
    13.8103,
    5.70494,
    2.29621,
    0.900496,
    0.344075,
    0.128088,
    0.0464534,
    414736};

  ///////////////////////////////////////////////////
  //The below PU distributions were used for the Summer 2011 result
  ///////////////////////////////////////////////////

  ////the data histogram obtained from: 
  //// /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Pileup_2011_EPS_8_jul.root
  //float TrueDist2011_f[25] = {
  //  1.45417e+07,
  //  3.47743e+07,
  //  7.89247e+07,
  //  1.26467e+08,
  //  1.59329e+08,
  //  1.67603e+08,
  //  1.52684e+08,
  //  1.23794e+08,
  //  9.09462e+07,
  //  6.13973e+07,
  //  3.8505e+07,
  //  2.2628e+07,
  //  1.25503e+07,
  //  6.61051e+06,
  //  3.32403e+06,
  //  1.60286e+06,
  //  743920,
  //  333477,
  //  144861,
  //  61112.7,
  //  25110.2,
  //  10065.1,
  //  3943.98,
  //  1513.54,
  //  896.161
  //};
  ////Summer11 PU_S4 distribution
  ////obtained from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupMCReweightingUtilities
  //Double_t PoissonIntDist_f[25] = {
  //  0.104109,
  //  0.0703573,
  //  0.0698445,
  //  0.0698254,
  //  0.0697054,
  //  0.0697907,
  //  0.0696751,
  //  0.0694486,
  //  0.0680332,
  //  0.0651044,
  //  0.0598036,
  //  0.0527395,
  //  0.0439513,
  //  0.0352202,
  //  0.0266714,
  //  0.019411,
  //  0.0133974,
  //  0.00898536,
  //  0.0057516,
  //  0.00351493,
  //  0.00212087,
  //  0.00122891,
  //  0.00070592,
  //  0.000384744,
  //  0.000219377
  //};

}
