#include <iostream>

using namespace std;


void calclsb(double a, double aerr, double b, double berr) {

  double perc = 100.0*(a-b)/( (a+b)/2.0 );
  double nsig = (a-b)/( sqrt(aerr*aerr+berr*berr) );


  cout << "perc: " << perc << ", nsig: " << nsig << endl; 

}
