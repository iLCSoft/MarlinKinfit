/*! \file 
 *  \brief Implements class TrackMomentumConstraint
 *
 * \b Changelog:
 * - 23.11.04 BL: First version (taken from TrackMomentumConstraint)
 *
 */ 

#include "TrackMomentumConstraint.h"
#include "TrackFitObject.h"

#include<iostream>
#include<cassert>

using std::cout;
using std::endl;

TrackMomentumConstraint::TrackMomentumConstraint (double pxfact_, double pyfact_, double pzfact_,
                          double efact_, double value_) 
    : factor (efact_, -pxfact_, -pyfact_, -pzfact_),  
      value (value_),
      cachevalid(false)
{}

TrackMomentumConstraint::TrackMomentumConstraint (int axis, double value_) 
    : value (value_),
      cachevalid(false)
{
  switch (axis) {
    case 1: factor.setValues (0, -1, 0, 0);
      break;
    case 2: factor.setValues (0, 0, -1, 0);
      break;
    case 3: factor.setValues (0, 0, 0, -1);
      break;
    default: factor.setValues (1, 0, 0, 0);
      break;
  }
}

// destructor
TrackMomentumConstraint::~TrackMomentumConstraint () {}

// calculate current value of constraint function
double TrackMomentumConstraint::getValue() const {
  FourVector sum;
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    sum += sign[i]*(fitobjects[i]->getMomentum (flags[i]));
//     cout <<  "TrackMomentumConstraint::getValue: i=" << i << ", sign=" << sign[i]
//          << ", mom=" << fitobjects[i]->getMomentum (flags[i])
//          << ", sign*mom=" << sign[i]*fitobjects[i]->getMomentum (flags[i])
//          << ", new sum=" << sum << endl;
  }
//   cout <<  "factor=" << factor << ", result=" << factor*sum - value << endl;
  return factor*sum - value;
}

// calculate vector/array of derivatives of this contraint 
// w.r.t. to ALL parameters of all fitobjects
// here: d sum(px) /d par(i,j) 
//                      = d sum(px) /d px(i) * d px(i) /d par(i, j)
//                                      =  1 * d px(i) /d par(i, j)
void TrackMomentumConstraint::getDerivatives(int idim, double der[]) const {
  for (int iglobal = 0; iglobal < idim; iglobal++) der[iglobal] = 0;
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    for (int ilocal = 0; ilocal < fitobjects[i]->getNPar(); ilocal++) {
      if (!fitobjects[i]->isParamFixed(ilocal)) {
        int iglobal = fitobjects[i]->getGlobalParNum (ilocal);
        assert (iglobal >= 0 && iglobal < idim);
        der[iglobal] += factor*(sign[i]*fitobjects[i]->getMomentumDerivative (flags[i], ilocal));
      }
    }
  }
}
  
void TrackMomentumConstraint::add1stDerivativesToMatrix(int idim, double *M) const {

  assert (0);

  assert (M);
  int kglobal = getGlobalNum();
  assert (kglobal >= 0 && kglobal < idim);
  
  for (ConstFitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
    const TrackFitObject *fo = *i;
    assert (fo);
    for (int ilocal = 0; ilocal < fo->getNPar(); ++ilocal) {
      if (!fo->isParamFixed(ilocal)) {
        int iglobal = fo->getGlobalParNum (ilocal);
        assert (iglobal >= 0 && iglobal < idim);
        double d = 0;
        M[idim*iglobal+kglobal] += d;
        M[idim*kglobal+iglobal] += d;
      }
    }
    
  }
}
  
void TrackMomentumConstraint::add2ndDerivativesToMatrix(int idim, double *M, double lambda) const {

  assert (0);
}
    
void TrackMomentumConstraint::addToGlobalDerMatrix (double lambda, int idim, double *M) const {

  assert (0);
  // Add lambda*d^2 g / d x_i dx_j to global matrix
  
  if (lambda == 0) return;
  
  // d^2 g / (dx_i dx_j) = 
  //   = sum_k,l d^2 g/(dpx_k dpx_l) * dpx_k/dx_i dpx_l/dx_j
  //     + sum_k dg/dpx_k * d^2 px_k/(dx_i dx_j)
  //   = sum_k,l      1              * dpx_k/dx_i dpx_l/dx_j
  
  // assume here that different 4-vectors always depend on 
  // different parameters!
  
  if (!cachevalid) updateCache();
  
  int *globalParNum = new int[nparams];
  double *der = new double[nparams];
  
//   ipar = 0;
//   for (int i = 0; i < fitobjects.size(); i++) {
//     for (int ilocal = 0; ilocal < fitobjects[i]->getNPar(); ilocal++) {
//       int iglobal = fitobjects[i]->getGlobalParNum (ilocal);
//       if (iglobal >= 0) {
//         assert (ipar < nparams);
//         globalParNum[ipar] = iglobal;
//         der[ipar] = fitobjects[i]->getDPx (ilocal);
//         ipar++;
//       }
//   }
  
  for (int ipar = 0; ipar < nparams; ipar++) {
    int iglobal = globalParNum[ipar];
    double der_i = der[ipar];
    for (int jpar = ipar; jpar < nparams; jpar++) {
      int jglobal = globalParNum[ipar];
      double der_j = der[jpar];
      double l_der_ij = lambda*der_i*der_j;
      M[idim*iglobal+jglobal] += l_der_ij;
      if (ipar != jpar) M[idim*jglobal+iglobal] += l_der_ij;
    }
  }       
}

void TrackMomentumConstraint::invalidateCache() const {
  cachevalid = false;
}

void TrackMomentumConstraint::updateCache() const {
  nparams = 0;
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    for (int ilocal = 0; ilocal < fitobjects[i]->getNPar(); ilocal++) {
      if (!fitobjects[i]->isParamFixed(ilocal)) {
        int iglobal = fitobjects[i]->getGlobalParNum (ilocal);
        if (iglobal >= 0) nparams++;
      }
    }
  }
  cachevalid = true;
}
  
