/*! \file 
 *  \brief Implements class MomentumConstraint
 *
 * \b Changelog:
 *
 * \b CVS Log messages:
 * - $Log: MomentumConstraint.cc,v $
 * - Revision 1.1  2008/02/18 09:59:35  blist
 * - MomentumConstraint and SoftGaussMomentumCOnstraint added; PConstraint is obsolete
 * -
 * -
 *
 */ 

#include "MomentumConstraint.h"
#include "ParticleFitObject.h"

#include<iostream>
#include<cassert>

using std::cout;
using std::endl;

MomentumConstraint::MomentumConstraint (double efact_, double pxfact_, double pyfact_, 
                                        double pzfact_, double value_) 
: efact (efact_),
  pxfact (pxfact_),
  pyfact (pyfact_),
  pzfact (pzfact_),
  value (value_),
  cachevalid(false),
  nparams(0)
{}

// destructor
MomentumConstraint::~MomentumConstraint () {
  //std::cout << "destroying MomentumConstraint" << std::endl;
}

// calculate current value of constraint function
double MomentumConstraint::getValue() const {
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
// here: d sum(px) /d par(i,j) 
//                      = d sum(px) /d px(i) * d px(i) /d par(i, j)
//                                      =  1 * d px(i) /d par(i, j)
void MomentumConstraint::getDerivatives(int idim, double der[]) const {
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

void MomentumConstraint::invalidateCache() const {
  cachevalid = false;
}

void MomentumConstraint::updateCache() const {
  nparams = 0;
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    for (int ilocal = 0; ilocal < fitobjects[i]->getNPar(); ilocal++) {
      int iglobal = fitobjects[i]->getGlobalParNum (ilocal);
      if (!fitobjects[i]->isParamFixed(ilocal)) {
        assert (iglobal >= 0);
        nparams++;
      }
    }
  }
  cachevalid = true;
}
  
bool MomentumConstraint::secondDerivatives (int i, int j, double *dderivatives) const {
  return false;
}  
  
bool MomentumConstraint::firstDerivatives (int i, double *dderivatives) const {
  dderivatives[0] = efact;
  dderivatives[1] = pxfact;
  dderivatives[2] = pyfact;
  dderivatives[3] = pzfact;
  return true;
}

int MomentumConstraint::getVarBasis() const {
  return VAR_BASIS;
}
