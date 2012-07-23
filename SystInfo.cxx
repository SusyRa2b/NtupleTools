#include "SystInfo.h"


void SystInfo::write(std::ofstream* outfile) const {

  (*outfile)<<status<<" "<<minus<<" "<<plus<<std::endl;

}

SystInfo::SystInfo(std::ifstream* input) :
  plus(0),minus(0),status(99)
{

  (*input)>>status>>minus>>plus;
  std::cout<< "  Loaded status,minus,plus = "<<status<<" "<<minus<<" "<<plus<<std::endl; //DEBUG
}

SystInfo::SystInfo(float p, float m, int s) :
  plus(p),  minus(m),  status(s) {}
SystInfo::~SystInfo() {}
