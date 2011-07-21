/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooTransverseThrustVar.h,v 1.2 2007/05/11 10:14:56 verkerke Exp $
 * Authors:                                                                  *
 *   GR, Gerhard Raven, NIKHEF/VU                                            *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_TRANSVERSE_THRUST_VAR
#define ROO_TRANSVERSE_THRUST_VAR

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include <vector>

class RooTransverseThrustVar: public RooAbsReal {
public:
  // Constructors, assignment etc
  //  inline RooThrustVar() : _nset(0), _fixedNset(0) { };
  virtual ~RooTransverseThrustVar();
  RooTransverseThrustVar(const char *name, const char *title, RooAbsReal& phi,
			 vector<double> jetPx, vector<double> jetPy);
  RooTransverseThrustVar(const RooTransverseThrustVar& other, const char* name=0);

  virtual TObject* clone(const char* newname) const { return new RooTransverseThrustVar(*this,newname); }

  //virtual Double_t getVal(const RooArgSet* set=0) const ;
  
protected:
  
  // Function evaluation
  virtual Double_t evaluate() const ;

  // the real stuff...
  RooRealProxy _phi ;               // Probability Density function positive
  std::vector<double> _jetPx ;
  std::vector<double> _jetPy ;
  std::vector<double> _jetPt ;

  ClassDef(RooTransverseThrustVar,0)             // 
};

#endif
