#include <iostream>
#include "MiscUtil.cxx"

using namespace std;


void contlsb(double A, double Acont, double Aconterr, double B, double Bcont, double Bconterr) {

  double scale = 30.15471/4683.719;

  Acont = scale*Acont;
  Aconterr = scale*Aconterr;
  Bcont = scale*Bcont;
  Bconterr = scale*Bconterr;

  double Aerr = sqrt(A);
  double Berr = sqrt(B);

  double r = A/B;
  double rerr = jmt::errAoverB(A,Aerr,B,Berr);
  cout << "without contamination sub: " << r << " +- " << rerr << endl;

  cout << "A: " << A << "->" << A-Acont << endl;
  cout << "B: " << B << "->" << B-Bcont << endl;

  double rc = (A-Acont)/(B-Bcont);
  double rcerr = jmt::errAoverB( A-Acont, sqrt(Aerr*Aerr+Aconterr*Aconterr) , B-Bcont, sqrt(Berr*Berr+Bconterr*Bconterr));
  cout << "with contamination sub: " << rc << " +- " << rcerr << endl;

  double a=r;
  double aerr = rerr;
  double b=rc;
  double berr = rcerr;

  double perc = 100.0*(a-b)/( (a+b)/2.0 );
  double nsig = (a-b)/( (aerr+berr)/2.0 );


  cout << "perc: " << perc << ", nsig: " << nsig << endl;

    

}
