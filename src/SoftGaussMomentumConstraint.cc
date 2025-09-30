/*! \file 
 *  \brief Implements class SoftGaussMomentumConstraint
 *
 * \b Changelog:
 *
 * \b CVS Log messages:
 * - $Log: SoftGaussMomentumConstraint.cc,v $
 * - Revision 1.1  2008/02/18 09:59:35  blist
 * - MomentumConstraint and SoftGaussMomentumCOnstraint added; PConstraint is obsolete
 * -
 * -
 *
 */ 

#include "SoftGaussMomentumConstraint.h"
#include "ParticleFitObject.h"

#include<iostream>
#include<cmath>
#undef NDEBUG
#include<cassert>

using std::cerr;
using std::cout;
using std::endl;

// constructor
SoftGaussMomentumConstraint::SoftGaussMomentumConstraint (double sigma_, double efact_, double pxfact_, 
                                                          double pyfact_, double pzfact_, double value_) 
: SoftGaussParticleConstraint (sigma_),
  efact (efact_),
  pxfact (pxfact_),
  pyfact (pyfact_),
  pzfact (pzfact_),
  value (value_)
{}

// destructor
SoftGaussMomentumConstraint::~SoftGaussMomentumConstraint () {
  // std::cout << "destroying SoftGaussMomentumConstraint" << std::endl;
}

// calulate current value of constraint function
double SoftGaussMomentumConstraint::getValue() const {
  double totpx = 0;
  double totpy = 0;
  double totpz = 0;
  double totE = 0;
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    if (pxfact != 0) totpx += fitobjects[i]->getPx(); 
    if (pyfact != 0) totpy += fitobjects[i]->getPy(); 
    if (pzfact != 0) totpz += fitobjects[i]->getPz(); 
    if (efact  != 0) totE  += fitobjects[i]->getE(); 
  }
  return pxfact*totpx + pyfact*totpy + pzfact*totpz + efact*totE - value;
}

// calculate vector/array of derivatives of this contraint 
// w.r.t. to ALL parameters of all fitobjects
// here: d M /d par(j) 
//          = d M /d p(i) * d p(i) /d par(j)
//          =  +-1/M * p(i) * d p(i) /d par(j)
void SoftGaussMomentumConstraint::getDerivatives(int idim, double der[]) const {
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    for (int ilocal = 0; ilocal < fitobjects[i]->getNPar(); ilocal++) {
      if (!fitobjects[i]->isParamFixed(ilocal)) {
        int iglobal = fitobjects[i]->getGlobalParNum (ilocal);
        assert (iglobal >= 0 && iglobal < idim);
        double d = 0;
        if (pxfact != 0) d += pxfact*fitobjects[i]->getDPx (ilocal);
        if (pyfact != 0) d += pyfact*fitobjects[i]->getDPy (ilocal);
        if (pzfact != 0) d += pzfact*fitobjects[i]->getDPz (ilocal);
        if (efact  != 0) d +=  efact*fitobjects[i]->getDE  (ilocal);
        der[iglobal] = d;
      }
    }
  }
}
  

bool SoftGaussMomentumConstraint::firstDerivatives (int, double *derivs) const {
  derivs[0] = efact;
  derivs[1] = pxfact;
  derivs[2] = pyfact;
  derivs[3] = pzfact;
  return true;
}

bool SoftGaussMomentumConstraint::secondDerivatives (int, int, double*) const {
  return false;
}  
  
