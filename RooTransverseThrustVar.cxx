/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id: RooEffProd.cxx 25184 2008-08-20 13:59:55Z wouter $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, NIKHEF
 *   GR, Gerhard Raven, NIKHEF/VU                                            *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/


/////////////////////////////////////////////////////////////////////////////////////
// BEGIN_HTML
// The class RooTransverseThrustVar implements the product of a PDF with an efficiency function.
// The normalization integral of the product is calculated numerically, but the
// event generation is handled by a specialized generator context that implements
// the event generation in a more efficient for cases where the PDF has an internal
// generator that is smarter than accept reject. 
// END_HTML
//

#include "RooFit.h"
#include "RooTransverseThrustVar.h"
#include <cmath>

ClassImp(RooTransverseThrustVar)
  ;



//_____________________________________________________________________________
RooTransverseThrustVar::RooTransverseThrustVar(const char *name, const char *title, 
			   RooAbsReal& phi,
			   vector<double> jetPx, vector<double> jetPy):
  RooAbsReal(name,title),
  _phi("phi","phi",this,phi),
  _jetPx(jetPx),
  _jetPy(jetPy)
{  
  unsigned int njets=  _jetPy.size();
  for (unsigned int i = 0; i<njets; i++) {
    _jetPt.push_back(sqrt(_jetPx[i]*_jetPx[i] + _jetPy[i]*_jetPy[i]));
  }
}




//_____________________________________________________________________________
RooTransverseThrustVar::RooTransverseThrustVar(const RooTransverseThrustVar& other, const char* name) : 
  RooAbsReal(other, name),
  _phi("phi",this,other._phi),
  _jetPx(other._jetPx),
  _jetPy(other._jetPy),
  _jetPt(other._jetPt)
{
  // Copy constructor
}




//_____________________________________________________________________________
RooTransverseThrustVar::~RooTransverseThrustVar() 
{
  // Destructor
}



//_____________________________________________________________________________
Double_t RooTransverseThrustVar::evaluate() const
{
  double thrust = 0;
  double sumMomentum = 0;
  double nx = cos(_phi);
  double ny = sin(_phi);
  unsigned int njets=  _jetPt.size();
  if(_jetPt.size() != _jetPx.size()) cout << "pT and pX do not match" << endl;
  if(_jetPt.size() != _jetPy.size()) cout << "pT and pY do not match" << endl;
  if(_jetPy.size() != _jetPx.size()) cout << "pX and pY do not match" << endl;
  for (unsigned int i = 0; i<njets; i++) {
    sumMomentum += _jetPt[i];
    thrust -= abs(nx*_jetPx[i] + ny*_jetPy[i]);
  }
  return thrust/sumMomentum;
}
