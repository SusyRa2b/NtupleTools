#include "CrossSectionTable.h"

/*
Class that loads a text file containing NLO cross sections across the mSugra plane
format is:

m0 m12 s0 s1 s2 s3 s4 s5 s6 s7 s8 s9

where the sN are the cross-sections for the production sub-processes given in the 
enum defined in the header file

*/

CrossSectionTable::CrossSectionTable(const TString & inputFile) //:
//  filename_(inputFile) {
{
  //load the contents of file into the database
  loadFileToDatabase(inputFile);
}

CrossSectionTable::~CrossSectionTable() {}

std::map<SUSYProcess,double>  CrossSectionTable::operator[] (const std::pair<int,int> & scanpoint) const {

  std::map<std::pair<int, int>, std::map<SUSYProcess, double> >::const_iterator mapval = database_.find(scanpoint);

  if (mapval == database_.end() ) {
    std::cout<<"[CrossSectionTable] did not find cross section for scan point "<<scanpoint.first<<" "<<scanpoint.second<<std::endl;
  }
  return mapval->second;

}

double CrossSectionTable::getCrossSection(const int m0, const int m12, const SUSYProcess process) const {

  std::map<SUSYProcess,double> forThisPoint = (*this)[std::make_pair(m0,m12)];
  return forThisPoint[process];

}

void CrossSectionTable::loadFileToDatabase(const TString & filename) {
  using namespace std;

  cout<<"Loading file "<<filename<<" into CrossSectionTable"<<endl;


  int m0,m12;
  double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10;
  ifstream file10(filename.Data());
  if (!file10.good() ) {
    cout<<"Problem with file! Terminating..."<<endl;
    assert(0);
  }
  database_.clear(); //just to be sure
  while ( file10>>m0>>m12>>c1>>c2>>c3>>c4>>c5>>c6>>c7>>c8>>c9>>c10 ) {
    pair<int, int> scanpoint = make_pair(m0,m12);

    map<SUSYProcess, double> theseCrossSections;
    theseCrossSections[ng] = c1;
    theseCrossSections[ns] = c2;
    theseCrossSections[nn] = c3;
    theseCrossSections[ll] = c4;
    theseCrossSections[sb] = c5;
    theseCrossSections[ss] = c6;
    theseCrossSections[tb] = c7;
    theseCrossSections[bb] = c8;
    theseCrossSections[gg] = c9;
    theseCrossSections[sg] = c10;
    database_[scanpoint] = theseCrossSections;
  }
  file10.close();
  cout<<"Loaded cross sections for "<<database_.size()<<" scan points"<<endl;

}
