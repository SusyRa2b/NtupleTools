
namespace pu {

  //See this twiki for more info:
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupMCReweightingUtilities

  /*
  this data distribution is an hadd of the following histograms:  
  (all in /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/)

  Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.pileup_v2.root
  Cert_165088-167913_7TeV_PromptReco_JSON.pileup_v2.root
  Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v2.pileup_v2.root
  Cert_172620-173692_PromptReco_JSON.pileup_v2.root
  Cert_175832-177515_PromptReco_JSON.pileup_v2.root
  Cert_177718_178078_7TeV_PromptReco_Collisons11_JSON.pileup_v2.root

  */
  float ObsDist2011_f[35] = {
    1.32542e+07,
    5.77956e+07,
    1.36378e+08,
    2.3e+08,
    3.11564e+08,
    3.62173e+08,
    3.77367e+08,
    3.63546e+08,
    3.31177e+08,
    2.89796e+08,
    2.46047e+08,
    2.03797e+08,
    1.65016e+08,
    1.30608e+08,
    1.00934e+08,
    7.60516e+07,
    5.57973e+07,
    3.98224e+07,
    2.76309e+07,
    1.86346e+07,
    1.22163e+07,
    7.78765e+06,
    4.8301e+06,
    2.91667e+06,
    1.71611e+06,
    984715,
    551540,
    301823,
    161527,
    84616.8,
    43429,
    21857,
    10795.3,
    5236.56,
    4645.72
  };


  // Summer11 PU_S4, distribution obtained by only looking at the in-time crossing.  This is the "spike+smear" distribution, RECOMMENDED FOR REWEIGHTING.
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
