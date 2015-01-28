/*! \file 
 *  \brief Implements class VertexFitObject
 *
 * \b Changelog:
 *
 * \b CVS Log messages:
 * - $Log: VertexFitObject.C,v $
 * - Revision 1.2  2007/09/13 13:33:06  blist
 * - Print methods return os
 * -
 *
 */ 

#include "VertexFitObject.h"
#include "TwoVector.h"
#include "ThreeVector.h"
#include "FourVector.h"
#include "VertexConstraint.h"
#include "TrackMomentumConstraint.h"
#include "TrackFitObject.h"
#include "BaseFitter.h"
#include "JBLHelix.h"
//#include "cernlib.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include <iostream>
using std::cout;
using std::endl;

#include <cassert>
#include <cmath>
using std::isfinite;

#include <cstring>

static int debug = 0;

VertexFitObject::VertexFitObject(const char *name_,
                                 double x,
                                 double y,
                                 double z
                                )
: covinvvalid (false), name (0),
  tracks (0), constraints (0)
{
  setParam (0, x, false);
  setParam (1, y, false);
  setParam (2, z, false);
  setMParam (0, x);
  setMParam (1, y);
  setMParam (2, z);
  setName (name_);
  initCov();
}

VertexFitObject::VertexFitObject (const VertexFitObject& rhs) 
: covinvvalid (false), name (0)
{
  copy (rhs);
}

VertexFitObject& VertexFitObject::operator= (const VertexFitObject& rhs) {
  if (&rhs != this) copy (rhs);
  return *this;
}

int VertexFitObject::getNPar() const {
  return NPAR;
}

void VertexFitObject::copy (const VertexFitObject& source) 
{
  covinvvalid =false;
  name = 0;
  for (int i = 0; i < NPAR; i++) {
    par[i] = source.par[i];
    mpar[i] = source.mpar[i];
    err[i] = source.err[i];
    measured[i] = source.measured[i];
    fixed[i] = source.fixed[i];
    globalParNum[i] = source.globalParNum[i];
    for (int j = 0; j < NPAR; j++) {
      cov[i][j] = source.cov[i][j];
    }
  } 
  tracks = source.tracks;
  setName (source.name);
}

VertexFitObject *VertexFitObject::copy() const {
  return new VertexFitObject (*this);
}
    
VertexFitObject& VertexFitObject::assign (const BaseFitObject& source) {
  if (const VertexFitObject *psource = dynamic_cast<const VertexFitObject *>(&source)) {
    if (psource != this) *this = *psource;
  }
  else {
    assert (0);
  }
  return *this;
}

VertexFitObject::~VertexFitObject()
{
  delete[] name;
  name = 0;
  for (CIterator it = constraints.begin(); it != constraints.end(); ++it) {
    delete (*it);
    (*it) = 0;
  }
}

bool VertexFitObject::setParam (int ilocal, double par_, 
                                 bool measured_, bool fixed_) {
  assert (ilocal >= 0 && ilocal < NPAR);
  measured[ilocal] = measured_;
  fixed[ilocal] = fixed_;
  return setParam (ilocal, par_);
} 

bool VertexFitObject::setParam (int ilocal, double par_ ) {
  if (!isfinite(par_)) return false;
  assert (ilocal >= 0 && ilocal < NPAR);
  par[ilocal] = par_;
  return true;
}  

bool VertexFitObject::setMParam (int ilocal, double mpar_ ) {
  if (!isfinite(mpar_)) return false;
  assert (ilocal >= 0 && ilocal < NPAR);
  mpar[ilocal] = mpar_;
  return true;
}  

bool VertexFitObject::setError (int ilocal, double err_) {
  if (!isfinite(err_)) return false;
  assert (ilocal >= 0 && ilocal < NPAR);
  cov[ilocal][ilocal] = err_*err_;
  covinvvalid = false;
  return true;
}

bool VertexFitObject::setCov (int ilocal, int jlocal, double cov_) {
  if (!isfinite(cov_)) return false;
  assert (ilocal >= 0 && ilocal < NPAR);
  assert (jlocal >= 0 && jlocal < NPAR);
  cov[ilocal][jlocal] = cov[jlocal][ilocal] = cov_;
  covinvvalid = false;
  return true;
}

const char *VertexFitObject::getName() const {
  return name ? name : "???";
}

const char *VertexFitObject::getParamName (int ilocal) const {
  switch (ilocal) {
    case 0: return "x";
    case 1: return "y";
    case 2: return "z";
  }
  return "undefined";
}

void VertexFitObject::setName(const char *name_) {
  if (!name_) return;
  delete[] name;
  name = 0;
  
  size_t len = strlen(name_)+1;
  name = new char[len];
  strncpy (name, name_, len);
}

  
bool VertexFitObject::setGlobalParNum (int ilocal, int iglobal) {
  assert (ilocal >= 0 && ilocal < NPAR);
  assert (!isParamFixed(ilocal));
  globalParNum[ilocal] = iglobal;
  return true;
}
bool VertexFitObject::fixParam (int ilocal, bool fix) {
  assert (ilocal >= 0 && ilocal < NPAR);
  return fixed [ilocal] = fix;
}
int  VertexFitObject::getGlobalParNum(int ilocal) const {
  if (ilocal < 0 || ilocal >= getNPar()) return -1;
  return globalParNum[ilocal];
}

double VertexFitObject::getParam (int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  return par[ilocal];
}
double VertexFitObject::getMParam (int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  assert (isParamMeasured(ilocal));
  return mpar[ilocal];
}

double VertexFitObject::getError (int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  return std::sqrt(cov[ilocal][ilocal]);
}

double VertexFitObject::getCov (int ilocal, int jlocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  assert (jlocal >= 0 && jlocal < NPAR);
  return cov[ilocal][jlocal];
}

bool VertexFitObject::isParamMeasured (int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  return measured[ilocal];
}

bool VertexFitObject::isParamFixed (int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  return fixed[ilocal];
}
    
std::ostream& VertexFitObject::print (std::ostream& os) const {
  os << getName() << ":\n";
  printParams(os);
  return os;
}


ThreeVector VertexFitObject::getVertex () const {
  ThreeVector result;
  getVertexEx (result);
  return result;
}

void VertexFitObject::getVertexEx (ThreeVector& p) const {
  p.setValues (par[0], par[1], par[2]);
} 
                                                    
ThreeVector VertexFitObject::getVertexDerivative (int ilocal) const  {
  ThreeVector result;
  getVertexDerivativeEx (ilocal, result);
  return result;
}

void VertexFitObject::getVertexDerivativeEx (int ilocal, ThreeVector& p) const {
  switch (ilocal) {
    case 0: // d / d par[0] = d / d x
      p.setValues (1, 0, 0);
      break;
    case 1: // d / d par[1] = d / d y
      p.setValues (0, 1, 0);
      break;
    case 2: // d / d par[2] = d / d z
      p.setValues (0, 0, 1);
      break;
    default: // should never happen!
      assert (0);
  }
}

void VertexFitObject::addToGlobCov(double *glcov, int idim) const {
  int globalnum[NPAR];
  bool ok [NPAR];
  for (int ilocal = 0; ilocal < getNPar(); ilocal++) {
    int iglobal = globalnum[ilocal] = getGlobalParNum(ilocal);
    assert (iglobal < idim);
    if (ok [ilocal] = (iglobal >= 0 && !isParamFixed(ilocal) && isParamMeasured (ilocal))) {
      for (int jlocal = 0; jlocal <= ilocal; jlocal++) {
        int jglobal = globalnum[jlocal];
        assert (jglobal < idim);
        if (ok [jlocal]) {
          double c = cov[ilocal][jlocal];
          glcov[jglobal+iglobal*idim] += c;
          if (jglobal != iglobal) glcov[iglobal+jglobal*idim] += c;
        }
      }
    }
  }
} 


void VertexFitObject::initCov() {
  int n = getNPar();
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      cov[i][j] = static_cast<double> (i == j);
    }
  }
  covinvvalid = false;
}

double VertexFitObject::getChi2() const {
  calculateChi2();
  return chi2;
}

bool VertexFitObject::calculateCovInv() const {
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
 *   int ierr = (idim == 0) ? 0 : dsinv(idim, &covinv[0][0], NPAR);
 *   if (ierr != 0) {
 *     std::cerr << "VertexFitObject::calculateCovInv: Error "
 *               << ierr << " from dsinv!" << std::endl;
 *   }
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


void VertexFitObject::calculateChi2() const {
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

void VertexFitObject::addVertexConstraints (BaseFitter& fitter, int axis) {
  for (TIterator it = tracks.begin(); it != tracks.end(); ++it) {
    assert(it->track); 
    VertexConstraint *con = new VertexConstraint (*this, *(it->track), (it->inbound ? 1 : 0), axis);
    fitter.addConstraint (con);
    constraints.push_back (con);
    it->track->releaseVertexParam (it->inbound ? 1 : 0);
  }
}

void VertexFitObject::addMomentumConstraint (BaseFitter& fitter, int axis) {
  TrackMomentumConstraint *con = new TrackMomentumConstraint (axis);
  fitter.addConstraint (con);
  constraints.push_back (con);
  
  // Add first inbound track
  for (TIterator it = tracks.begin(); it != tracks.end(); ++it) {
    if (it->inbound) {
      assert(it->track); 
      con->addToFOList(*(it->track), 1);
      break;
    }
  }
  // Now add outgoing tracks
  for (TIterator it = tracks.begin(); it != tracks.end(); ++it) {
    assert(it->track); 
    if (!it->inbound) con->addToFOList(*(it->track), 0);
  }
}

void VertexFitObject::addConstraints (BaseFitter& fitter, int mask) {
  if (mask & VX) addVertexConstraints (fitter, 0); else fixParam (0);
  if (mask & VY) addVertexConstraints (fitter, 1); else fixParam (1);
  if (mask & VZ) addVertexConstraints (fitter, 2); else fixParam (2);
  if (mask & PX) addMomentumConstraint (fitter, 1);
  if (mask & PY) addMomentumConstraint (fitter, 2);
  if (mask & PZ) addMomentumConstraint (fitter, 3);
  if (mask & E)  addMomentumConstraint (fitter, 0);
}

void VertexFitObject::addTrack (TrackFitObject *track, bool inbound, bool measured) {
  assert (track);
  tracks.push_back (TrackDescriptor (track, inbound, measured));
}

ThreeVector VertexFitObject::estimatePosition () {
  if (debug) cout << "VertexFitObject::estimatePosition(): starting" << endl;
  ThreeVector position (0, 0, 0);
  int n = 0;
  for (TIterator it0 = tracks.begin(); it0 != tracks.end(); ++it0) {
    if (it0->measured) {
      JBLHelix h0 = it0->track->getTangentialHelix(0);
      TIterator it1 = it0;
      for (++it1; it1 != tracks.end(); ++it1) {
        if (it1->measured) {
          if (debug)
            cout << "Intersecting " << it0->track->getName()
                 << " and " << it1->track->getName() << endl;
          JBLHelix h1 = it1->track->getTangentialHelix(0);
          double s0, s1, s02nd, s12nd;
          h0.getClosestApproach (h1, s0, s1, s02nd, s12nd);
          position += h0.getTrajectoryPoint (s0);
          position += h1.getTrajectoryPoint (s1);
          if (debug)
            cout << "  point 0: " << h0.getTrajectoryPoint (s0)
                 << ", point 1: " << h1.getTrajectoryPoint (s1) << endl;
          n += 2;
        }
      }
    }
  }
  position *= (1./n);
  if (debug) cout << "Final position estimate: " << position << endl;
  return position;
}

void VertexFitObject::initForFit() {
  if (debug) 
    cout << "VertexFitObject::initForFit(): starting for " << getName() << endl;

  // Estimate and set vertex position
  ThreeVector position = estimatePosition ();
  for (int i = 0; i < 3; i++) par[i] = position.getComponent (i);
  if (debug) 
    cout << "VertexFitObject now: " << *this << endl;
  if (debug) {
    cout << "TrackFit objects before adjustment: \n";
    for (TIterator it = tracks.begin(); it != tracks.end(); ++it) {
      cout << "   " << it->track->getName() << " = " << *(it->track)
           << (it->measured ? "  measured " : "  not meas.")
           << (it->inbound ? " decays at " : " starts at ")
           << it->track->getVertex(it->inbound ? 1 : 0)
           << endl;
    }
  }  
    
  // For all measured tracks, set s to be close to vertex position,
  // at the same time estimate total 4-momentum and total charge
  FourVector ptot;
  double chargetot = 0;
  int nunm = 0;
  for (TIterator it = tracks.begin(); it != tracks.end(); ++it) {
    if (it->measured) {
      if (it->inbound) {
        it->track->setVertex (1, TwoVector (par[0], par[1]));
        ptot -= it->track->getMomentum (1);
        chargetot -= it->track->getCharge();
      }
      else {
        it->track->setVertex (0, TwoVector (par[0], par[1]));
        ptot += it->track->getMomentum (0);
        chargetot += it->track->getCharge();
      }
    }
    else {
      nunm++;
    }
  }
  assert (nunm <= 1);
  // Now, init remaining unmeasured track
  for (TIterator it = tracks.begin(); it != tracks.end(); ++it) {
    if (!it->measured) {
      if (it->inbound) {
        assert (ptot.getE() > 0);
        it->track->setParameters (1, position, ptot, chargetot);
      }
      else {
        assert (ptot.getE() < 0);
        it->track->setParameters (0, position, -ptot, -chargetot);
      }
    }
  }
  if (debug) {
    cout << "TrackFit objects after adjustment: \n";
    for (TIterator it = tracks.begin(); it != tracks.end(); ++it) {
      cout << "   " << it->track->getName() << " = " << *(it->track)
           << (it->measured ? "  measured " : "  not meas.")
           << (it->inbound ? " decays at " : " starts at ")
           << it->track->getVertex(it->inbound ? 1 : 0)
           << endl;
    }
  }  
}
