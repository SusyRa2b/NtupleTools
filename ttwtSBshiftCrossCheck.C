#include "MiscUtil.cxx"

void docalc(double nsb,double nsb_err,double sbsub,double sbsub_err,double nsbsl,double nsbsl_err,double &r, double &r_err) {

  double sb = nsb - sbsub;
  double sb_err = sqrt(nsb_err*nsb_err + sbsub_err*sbsub_err);
  r = sb/nsbsl;
  r_err = jmt::errAoverB( sb,sb_err,nsbsl,nsbsl_err);

}

void sb()
{

  double r,r_err;
  const  double trigeff=0.841;
  const  double trigeff_perr =5.9e-2 ;
  const  double trigeff_merr =9.0e-2 ;

  //1BT
  double nsb = 173;
  double sbsub = 41;
  double sbsub_err = 9;
  double nsbsl = 81;

  //calc sb / sbsl
  docalc(nsb,sqrt(nsb),sbsub,sbsub_err,nsbsl,sqrt(nsbsl),r,r_err);
  cout<<"shifted SB r = "<<r<<" +/- "<<r_err<<endl;

  //1BT nominal with eff correction
  nsb = 328;
  sbsub = 108;
  sbsub_err = 13;
  nsbsl = 173;

  docalc(nsb / trigeff, sqrt(nsb)/trigeff, sbsub,sbsub_err,nsbsl,sqrt(nsbsl),r,r_err);
  double rnom=r; double rnom_err=r_err;
  docalc(nsb / (trigeff+trigeff_perr), sqrt(nsb)/(trigeff+trigeff_perr), sbsub,sbsub_err,nsbsl,sqrt(nsbsl),r,r_err);
  double rp=r;
  docalc(nsb / (trigeff-trigeff_merr), sqrt(nsb)/(trigeff-trigeff_merr), sbsub,sbsub_err,nsbsl,sqrt(nsbsl),r,r_err);
  double rm=r;
  cout<<"nominal    r = "<<rnom<<" +/- "<<rnom_err<< " - "<<rnom-rp<<" + "<<rm-rnom<<endl;

  cout<<" == 2BT =="<<endl;
  //2BT shifted
  nsb = 34;
  sbsub = 4;
  sbsub_err = 2;
  nsbsl = 18;

  //calc sb / sbsl
  docalc(nsb,sqrt(nsb),sbsub,sbsub_err,nsbsl,sqrt(nsbsl),r,r_err);
  cout<<"shifted SB r = "<<r<<" +/- "<<r_err<<endl;
  

  //2BT nominal with eff correction
  nsb = 52 ;
  sbsub = 10;
  sbsub_err = 3;
  nsbsl = 41;

  docalc(nsb / trigeff, sqrt(nsb)/trigeff, sbsub,sbsub_err,nsbsl,sqrt(nsbsl),r,r_err);
  rnom=r; rnom_err=r_err;
  docalc(nsb / (trigeff+trigeff_perr), sqrt(nsb)/(trigeff+trigeff_perr), sbsub,sbsub_err,nsbsl,sqrt(nsbsl),r,r_err);
  rp=r;
  docalc(nsb / (trigeff-trigeff_merr), sqrt(nsb)/(trigeff-trigeff_merr), sbsub,sbsub_err,nsbsl,sqrt(nsbsl),r,r_err);
  rm=r;
  cout<<"nominal    r = "<<rnom<<" +/- "<<rnom_err<< " - "<<rnom-rp<<" + "<<rm-rnom<<endl;

}
