/*! \file 
 *  \brief Implements class TrackFitObject
 *
 * \b Changelog:
 *
 * \b CVS Log messages:
 * - $Log: TrackFitObject.C,v $
 * - Revision 1.2  2007/09/13 13:33:06  blist
 * - Print methods return os
 * -
 *
 */ 

#include "TrackFitObject.h"
#include "ThreeVector.h"
#include "FourVector.h"
//#include "cernlib.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdio>

using std::isfinite;

#include <cstring>

TrackFitObject::TrackFitObject(const char *name_ )
: cachevalid (false), covinvvalid (false), name (0)
{
  for (int i = 0; i < NPARMAX; ++i) {
    fixed[i] = true;
    measured [i] = false;
    par[i] = mpar[ i] = err [i] = 0;
    globalParNum[i] = -1;
    for (int j = 0; j < NPARMAX; ++j) cov[i][j] = covinv[i][j] = 0; 
  }
  setName (name_);
}

TrackFitObject::TrackFitObject (const TrackFitObject& rhs) 
: cachevalid (false), covinvvalid (false), name (0)
{
  copy (rhs);
}
TrackFitObject& TrackFitObject::operator= (const TrackFitObject& rhs) {
  if (&rhs != this) copy (rhs);
  return *this;
}

void TrackFitObject::copy (const TrackFitObject& source) 
{
  cachevalid = false;
  covinvvalid = false;
  name = 0;
  for (int i = 0; i < NPARMAX; i++) {
    par[i] = source.par[i];
    mpar[i] = source.mpar[i];
    err[i] = source.err[i];
    measured[i] = source.measured[i];
    fixed[i] = source.fixed[i];
    globalParNum[i] = source.globalParNum[i];
    for (int j = 0; j < NPARMAX; j++) {
      cov[i][j] = source.cov[i][j];
    }
  } 
  setName (source.name);
}


TrackFitObject::~TrackFitObject()
{
  delete[] name;
}

bool TrackFitObject::setParam (int i, double par_, 
                                 bool measured_, bool fixed_) {
  assert (i >= 0 && i < NPARMAX);
  if (measured[i] != measured_ || fixed[i] != fixed_) invalidateCache();
  measured[i] = measured_;
  fixed[i] = fixed_;
  return setParam (i, par_);
}  

bool TrackFitObject::setParam (int i, double par_ ) {
  if (!isfinite(par_)) return false;
  assert (i >= 0 && i < NPARMAX);
  if (par[i] == par_) return true;
  invalidateCache();
  par[i] = par_;
  return true;
}  
bool TrackFitObject::setMParam (int i, double mpar_ ) {
  if (!isfinite(mpar_)) return false;
  assert (i >= 0 && i < NPARMAX);
  if (mpar[i] == mpar_) return true;
  invalidateCache();
  mpar[i] = mpar_;
  return true;
}  

bool TrackFitObject::setError (int ilocal, double err_) {
  if (!isfinite(err_)) return false;
  assert (ilocal >= 0 && ilocal < NPARMAX);
  invalidateCache();
  cov[ilocal][ilocal] = err_*err_;
  covinvvalid = false;
  return true;
}

bool TrackFitObject::setCov (int ilocal, int jlocal, double cov_) {
  if (!isfinite(cov_)) return false;
  assert (ilocal >= 0 && ilocal < NPARMAX);
  assert (jlocal >= 0 && jlocal < NPARMAX);
  invalidateCache();
  covinvvalid = false;
  cov[ilocal][jlocal] = cov[jlocal][ilocal] = cov_;
  return true;
}

const char *TrackFitObject::getName() const {
  return name ? name : "???";
}

void TrackFitObject::setName(const char *name_) {
  if (!name_) return;
  delete[] name;
  name = 0;
  
  size_t len = strlen(name_)+1;
  name = new char[len];
  strncpy (name, name_, len);
}

// B field in Tesla
double TrackFitObject::bfield = 1.14;
//double TrackFitObject::bfield = 4.;

double TrackFitObject::setBfield (double bfield_) {
  return bfield = bfield_;
}
  
bool TrackFitObject::setGlobalParNum (int ilocal, int iglobal) {
  assert (ilocal >= 0 && ilocal < NPARMAX);
  assert (!isParamFixed(ilocal));
  globalParNum[ilocal] = iglobal;
  return true;
}
bool TrackFitObject::fixParam (int ilocal, bool fix) {
  assert (ilocal >= 0 && ilocal < NPARMAX);
  return fixed [ilocal] = fix;
}
int  TrackFitObject::getGlobalParNum(int ilocal) const {
  if (ilocal < 0 || ilocal >= getNPar()) return -1;
  return globalParNum[ilocal];
}

double TrackFitObject::getParam (int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPARMAX);
  return par[ilocal];
}
double TrackFitObject::getMParam (int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPARMAX);
  assert (isParamMeasured(ilocal));
  return mpar[ilocal];
}

double TrackFitObject::getError (int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPARMAX);
  return std::sqrt(cov[ilocal][ilocal]);
}

double TrackFitObject::getCov (int ilocal, int jlocal) const {
  assert (ilocal >= 0 && ilocal < NPARMAX);
  assert (jlocal >= 0 && jlocal < NPARMAX);
  return cov[ilocal][jlocal];
}

bool TrackFitObject::isParamMeasured (int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPARMAX);
  return measured[ilocal];
}

bool TrackFitObject::isParamFixed (int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPARMAX);
  return fixed[ilocal];
}
    
std::ostream&  TrackFitObject::print (std::ostream& os) const {
  os << getName() << ":\n";
  printParams(os);
  os << "\n     4-vect(0): " << getMomentum (0);
  os << "\n     vertex(0): " << getVertex (0);
  return os;
}

void TrackFitObject::invalidateCache() {
  cachevalid = false;
}

ThreeVector TrackFitObject::getTrajectoryPoint (double s) const {
  ThreeVector result;
  getTrajectoryPointEx (s, result);
  return result;
}

ThreeVector TrackFitObject::getVertex (int i) const {
  ThreeVector result;
  getVertexEx (i, result);
  return result;
}

ThreeVector TrackFitObject::getTrajectoryDerivative (double s, 
                                                     int ilocal 
                                                    ) const {
  ThreeVector result;
  getTrajectoryDerivativeEx (s, ilocal, result);
  return result;
}
                                                    
ThreeVector TrackFitObject::getVertexDerivative (int i, 
                                                 int ilocal 
                                                ) const  {
  ThreeVector result;
  getVertexDerivativeEx (i, ilocal, result);
  return result;
}

FourVector TrackFitObject::getMomentumAtTrajectory (double s) const {
  FourVector result;
  getMomentumAtTrajectoryEx (s, result);
  return result;
}

 FourVector TrackFitObject::getMomentum (int i) const {
  FourVector result;
  getMomentumEx (i, result);
  return result;
}

FourVector TrackFitObject::getMomentumDerivativeAtTrajectory (double s, 
                                                              int ilocal 
                                                             ) const {
  FourVector result;
  getMomentumDerivativeAtTrajectoryEx (s, ilocal, result);
  return result;
}
                                                    
FourVector TrackFitObject::getMomentumDerivative (int i, 
                                                  int ilocal 
                                                 ) const  {
  FourVector result;
  getMomentumDerivativeEx (i, ilocal, result);
  return result;
}

void TrackFitObject::addToGlobCov(double *glcov, int idim) const {
  int globalnum[NPARMAX];
  bool ok [NPARMAX];
  for (int ilocal = 0; ilocal < getNPar(); ilocal++) {
    int iglobal = globalnum[ilocal] = getGlobalParNum(ilocal);
    if (ok [ilocal] = (iglobal >= 0 && !isParamFixed(ilocal) && isParamMeasured (ilocal))) {
      for (int jlocal = 0; jlocal < ilocal; jlocal++) {
        if (ok [jlocal]) {
          int jglobal = globalnum[jlocal];
          glcov[jglobal+iglobal*idim] += cov[ilocal][jlocal];
          glcov[iglobal+jglobal*idim] += cov[ilocal][jlocal];
        }
      }
      glcov[iglobal*(idim+1)] += cov[ilocal][ilocal];
    }
  }
  // printCov (std::cout);
} 


void TrackFitObject::initCov() {
  int n = getNPar();
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      cov[i][j] = static_cast<double>(i == j);
    }
  }
  covinvvalid = false;
}

void TrackFitObject::checkCov() {
  // Check that correlations are not greater than 1
  // (can happen due to resolution effects)
  int n = getNPar();
  for (int i = 0; i < n; ++i) {
    if (isParamMeasured (i)) {
      for (int j = 0; j < n; ++j) {
        if (isParamMeasured (j)) {
          double cc = cov[i][i]*cov[j][j];
          if (cov[i][j]*cov[i][j] > cc) {
            cov[i][j] = cov[j][i] = (cov[i][j] > 0 ? 0.999 : -0.999) *
                                    std::sqrt(cc);
            covinvvalid = false;
          } 
        } 
      } 
    }
  }
  int it = 0;
  while (it++ < 10) {
    if (calculateCovInv()) break;
    std::cout << "TrackFitObject::checkCov: rescaling" << std::endl;
    int imax = 0;
    int jmax = 1;
    double rmax = 0;
    for (int i = 0; i < n; ++i) {
      if (isParamMeasured (i)) {
        for (int j = i+1; j < n; ++j) {
          if (isParamMeasured (j)) {
            double r = std::abs(cov[i][j]/(std::sqrt(cov[i][i]*cov[j][j])));
            if (r > rmax) {
              rmax = r;
              imax = i;
              jmax = j;
            } 
          } 
        } 
      }
    }
    cov[imax][jmax] *= 0.99;
    cov[jmax][imax] *= 0.99;
  }
}

double TrackFitObject::getChi2() const {
  calculateChi2();
  return chi2;
}

bool TrackFitObject::calculateCovInv() const {
  int n = getNPar();
  
/*
 *   int idim = 0;
 *   for (int i = 0; i < n; ++i) {
 *     if (isParamMeasured (i)) {
 *       idim = i;
 *       for (int j = 0; j < n; ++j) {
 *         covinv[i][j] = isParamMeasured (j) ? cov[i][j] : 0;
 *       }
 *     }
 *     else {
 *       for (int j = 0; j < n; ++j) {
 *         covinv[i][j] = static_cast<double>(i == j);
 *       }
 *     }
 *   }
 *   int ierr = (idim == 0) ? 0 : dsinv(idim, &covinv[0][0], NPARMAX);
 * //   if (ierr != 0) {
 * //     std::cout << "TrackFitObject::calculateCovInv: Error "
 * //               << ierr << " from dsinv! Object " << getName() << std::endl;
 * //     printCov (std::cout);         
 * //   }
 *   return covinvvalid = (ierr == 0);
 */ 
 
  gsl_matrix *covm = gsl_matrix_alloc (n, n);
  gsl_matrix_set_identity (covm);
  
  for (int i = 0; i < n; ++i) {
    if (isParamMeasured (i)) {
      for (int j = 0; j < n; ++j) {
        if (isParamMeasured (j)) gsl_matrix_set (covm, i, j, cov[i][j]);
      }
    }
  }
  gsl_error_handler_t *e = gsl_set_error_handler_off ();
  int result = gsl_linalg_cholesky_decomp (covm);
  if (result == 0) result = gsl_linalg_cholesky_invert (covm);
  gsl_set_error_handler (e);
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      covinv[i][j] = covinv[j][i] = gsl_matrix_get (covm, i, j);
    }
    covinv[i][i] = gsl_matrix_get (covm, i, i);
  }

  gsl_matrix_free(covm);
  covinvvalid = (result == 0);
  return covinvvalid;
 
}


void TrackFitObject::calculateChi2() const {
  if (!covinvvalid) calculateCovInv();
  if (!covinvvalid) {
    chi2 = -1;
    return;
  }
  chi2 = 0;
  for (int i = 0; i < getNPar(); ++i) {
    resid[i] = par[i]-mpar[i];
    if (chi2contr[i] = isParamMeasured(i) && !isParamFixed(i)) {
      chi2 += resid[i]*covinv[i][i]*resid[i];
      for (int j = 0; j < i; ++j) {
        if (chi2contr[j]) chi2 += 2*(resid[i])*covinv[i][j]*(resid[j]);
      }
    }
  }
}

bool TrackFitObject::releaseVertexParam (int ivertex) {
  return fixVertexParam (ivertex, false);
}

void TrackFitObject::printCov (std::ostream& os) const {
  os << "Covariance matrix of " << name << "\n";
  for(int k = 0; k < getNPar(); ++k) {
    for (int l = 0; l < getNPar(); ++l) {
      static char str[20];
      if (l <= k) {
        snprintf (str, 20, " %10.6f", cov[k][l]);
      }
      else if (l == k+1) {
        snprintf (str, 20, " |%9.6f", cov[k][l]/std::sqrt(cov[k][k]*cov[l][l]));
      }
      else {
        snprintf (str, 20, " %10.6f", cov[k][l]/std::sqrt(cov[k][k]*cov[l][l]));
      }
      os << str;
    }
    os << "   " << getParamName (k) << "\n";
  }        
}
