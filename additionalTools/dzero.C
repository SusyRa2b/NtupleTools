#include <iostream>
#include "MiscUtil.cxx"

using namespace std;

void dzero() {
  
  //inputs for AN, April 13 
  double myR = 0.115263;
  double myRerr = 0.00447743;
  double D = 2.03874;
  double Derr = 1.41421;
  double Dsub = 2.05436;
  double Dsuberr = 0.223389;
  double eff_MHT = 0.981;
  double MHTLDPDeltaDown = -0.00278582;
  double MHTLDPDeltaUp = 0.00878194;
  double MHTDeltaDown = 0.945;
  double MHTDeltaUp = 0.993;
  double closure = 318.1*.01;
  double estm40 = 0.0911502;
  //

  //crack noF

  //crack F

  

  //////////////////////////////////////////
  
  double estimate = 0;
  double estimateerr = jmt::errAtimesB(myR, myRerr, D-Dsub, sqrt(Derr*Derr+Dsuberr*Dsuberr));
  

  double delta_down = MHTLDPDeltaDown;
  double delta_up = MHTLDPDeltaUp;
  double estFluctuate = estimateerr;
  double delta_down2 = (MHTDeltaDown-eff_MHT)*estFluctuate;
  double delta_up2 = (MHTDeltaUp-eff_MHT)*estFluctuate;
  
  double estimateerrTrigPlus  = sqrt(delta_up*delta_up + delta_up2*delta_up2);
  double estimateerrTrigMinus = sqrt(delta_down*delta_down + delta_down2*delta_down2);
  
  double systm = sqrt(estm40*estm40*closure*closure + estimateerrTrigMinus*estimateerrTrigMinus);
  double systp = sqrt(estm40*estm40*closure*closure + estimateerrTrigPlus*estimateerrTrigPlus);
  
  cout << "Estimate for QCD results table: " << endl;
  cout << estimate << " +- " << estimateerr << " - " << estimateerrTrigMinus << " + " << estimateerrTrigPlus << endl;
  
  cout << "Estimate for results summary: " << endl;
  cout << estimate << " +- " << estimateerr << " - " << systm << " + " << systp << endl;
    

}
