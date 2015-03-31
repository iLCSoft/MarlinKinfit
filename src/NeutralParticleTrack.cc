/*! \file 
 *  \brief Implements class NeutralParticleTrack
 *
 * \b Changelog:
 * - 23.11.04 BL First version
 *
 */ 
#include "NeutralParticleTrack.h"
#include "TwoVector.h"
#include "ThreeVector.h"
#include "FourVector.h"
#include "JBLHelix.h"

#include <cassert>
#include <cmath>
using std::sin;
using std::cos;
using std::sqrt;
using std::tan;
using std::pow;
using std::abs;
using std::isfinite;
#include <iostream>
using std::cout;
using std::endl;


const double NeutralParticleTrack::parfact[NPARMAX] = 
  {1., 1., 1., 1., 1., 1., 1., 1.};

NeutralParticleTrack::NeutralParticleTrack (const char *name_,
                                            double pt,
                                            double phi0,
                                            double theta,
                                            double dca,
                                            double z0,
                                            double mass,
                                            double sstart,
                                            double sstop
                                            )
: TrackFitObject (name_) 
{
  setParam (0, pt, false);
  setParam (1, phi0, false);
  setParam (2, theta, false);
  setParam (3, dca, false);
  setParam (4, z0, false);
  setParam (5, mass, false, true);
  setParam (6, sstart, false);
  setParam (7, sstop, false);
  // setParam (6, sstop, false);
  setMParam (0, pt);
  setMParam (1, phi0);
  setMParam (2, theta);
  setMParam (3, dca);
  setMParam (4, z0);
  TrackFitObject::initCov();
}
NeutralParticleTrack::NeutralParticleTrack (const char *name_,
                                            double pt,
                                            double phi0,
                                            double theta,
                                            double dca,
                                            double z0,
                                            double mass,
                                            double sstart
                                            )
: TrackFitObject (name_) 
{
  setParam (0, pt, false);
  setParam (1, phi0, false);
  setParam (2, theta, false);
  setParam (3, dca, false);
  setParam (4, z0, false);
  setParam (5, mass, false, true);
  setParam (6, sstart, false);
  setParam (7, 0, false, true);
  // setParam (6, sstop, false);
  setMParam (0, pt);
  setMParam (1, phi0);
  setMParam (2, theta);
  setMParam (3, dca);
  setMParam (4, z0);
  TrackFitObject::initCov();
}
NeutralParticleTrack::NeutralParticleTrack (const char *name_,
                                            const ThreeVector& vertex,
                                            const ThreeVector& momentum,
                                            double mass
                                           )
: TrackFitObject (name_) 
{
  pt = momentum.getPt();
  theta = momentum.getTheta();
  phi0 = momentum.getPhi();
  
//   cout << "NeutralParticleTrack::NeutralParticleTrack"
//        << "\nvertex = " << vertex << ", momentum = " << momentum
//        << "\npt=" << pt << ", theta=" << theta << ", phi0=" << phi0 << endl;
  
  sphi0 = sin(phi0);
  cphi0 = cos(phi0);
  double sstop = vertex.getX()*cphi0 + vertex.getY()*sphi0;
  dca = vertex.getX()*sphi0 - vertex.getY()*cphi0; 
  z0 = vertex.getZ() - sstop*cos(theta)/sin(theta);

//   cout << "sstart=" << sstart << ", dca=" << dca << endl;
//   cout << "x(sstart)=" << sstart*cphi0 + dca*sphi0
//        << ", y(sstart)=" << +sstart*sphi0 - dca*cphi0 << endl;
  

  setParam (0, pt, false);
  setParam (1, phi0, false);
  setParam (2, theta, false);
  setParam (3, dca, false);
  setParam (4, z0, false);
  setParam (5, mass, false, true);
  setParam (6, 0, false, true);
  setParam (7, sstop, false);
  
//   printParams (cout);
//   cout << endl << endl;
  
  setMParam (0, pt);
  setMParam (1, phi0);
  setMParam (2, theta);
  setMParam (3, dca);
  setMParam (4, z0);
  TrackFitObject::initCov();
}

NeutralParticleTrack *NeutralParticleTrack::copy() const {
  return new NeutralParticleTrack (*this);
}
    
NeutralParticleTrack& NeutralParticleTrack::assign (const BaseFitObject& source) {
  if (const NeutralParticleTrack *psource = dynamic_cast<const NeutralParticleTrack *>(&source)) {
    if (psource != this) *this = *psource;
  }
  else {
    assert (0);
  }
  return *this;
}

NeutralParticleTrack::~NeutralParticleTrack ()
{}

int NeutralParticleTrack::getNPar() const {
  return NPAR;
}

const char *NeutralParticleTrack::getParamName (int ilocal) const {
  switch (ilocal) {
    case 0: return "pt";
    case 1: return "phi0";
    case 2: return "theta";
    case 3: return "dca";
    case 4: return "z0";
    case 5: return "mass";
    case 6: return "sstart";
    case 7: return "sstop";
  }
  return "undefined";
}

bool NeutralParticleTrack::setParam (int ilocal, double par_, 
                                 bool measured_, bool fixed_) {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (measured[ilocal] != measured_ || fixed[ilocal] != fixed_) invalidateCache();
  measured[ilocal] = measured_;
  fixed[ilocal] = fixed_;
  return setParam (ilocal, par_);
}  

bool NeutralParticleTrack::setParam (int ilocal, double par_ ) {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!isfinite(par_)) return false;
  // enforce parameter range restrictions
  // dont restrict phi0 to -pi ... pi, because otherwise residual calculation
  // becomes too complicated!
  switch (ilocal) {
    case 0:  if (par_ <= 0) return false;
      break;
    // case 1:  if (abs(par_) > M_PI) par_ = atan2 (sin (par_), cos (par_); // -pi <= phi0 <= pi
    //   break;
    case 2:  if (par_ < 0 || par_*parfact[2] >= M_PI) par_ = acos (cos (par_*parfact[2]))/parfact[2];
      break;
    case 5:  if (par_ < 0) par_ = abs (par_);
      break;
  }
  if (par[ilocal] == par_) return true;
  invalidateCache();
  par[ilocal] = par_;
  return true;
}  

bool NeutralParticleTrack::setParameters (int ivertex, 
                                          const ThreeVector& vertex,  
                                          const FourVector& momentum,
                                          double charge_    
                                         ) { 
  bool result = true;
  pt = momentum.getPt();
  theta = momentum.getTheta();
  phi0 = momentum.getPhi();
  
//   cout << "NeutralParticleTrack::NeutralParticleTrack"
//        << "\nvertex = " << vertex << ", momentum = " << momentum
//        << "\npt=" << pt << ", theta=" << theta << ", phi0=" << phi0 << endl;
  
  sphi0 = sin(phi0);
  cphi0 = cos(phi0);
  double s = vertex.getX()*cphi0 + vertex.getY()*sphi0;
  dca = vertex.getX()*sphi0 - vertex.getY()*cphi0; 
  z0 = vertex.getZ() - s*cos(theta)/sin(theta);

//   cout << "sstart=" << sstart << ", dca=" << dca << endl;
//   cout << "x(sstart)=" << sstart*cphi0 + dca*sphi0
//        << ", y(sstart)=" << +sstart*sphi0 - dca*cphi0 << endl;
  

  if (!isParamFixed (0)) result &= setParam (0, pt/parfact[0]);
  if (!isParamFixed (1)) result &= setParam (1, phi0/parfact[1]);
  if (!isParamFixed (2)) result &= setParam (2, theta/parfact[2]);
  if (!isParamFixed (3)) result &= setParam (3, dca/parfact[3]);
  if (!isParamFixed (4)) result &= setParam (4, z0/parfact[4]);
  if (!isParamFixed (5)) result &= setParam (5, mass/parfact[5]);
  assert (ivertex >= 0 && 6+ivertex < NPAR);
  result &= setParam (6+ivertex, s/parfact[6+ivertex]);
  return result;
}


void NeutralParticleTrack::getTrajectoryPointEx (double s, ThreeVector& p) const {
  if (!cachevalid) updateCache();
    p.setValues (x0 + cphi0*s,
                 y0 + sphi0*s,
                 z0 + s*cottheta);
}

void NeutralParticleTrack::getTrajectoryDerivativeEx (double s, int ilocal, ThreeVector& p) const {
  if (!cachevalid) updateCache();
  assert (ilocal >= 0 && ilocal < NPAR);
  switch (ilocal) {
    case 0: // d / d pt
      p.setValues (0, 0, 0);          
      break;
    case 1: // d / d phi0
      p.setValues (-y0 - sphi0*s, x0 + cphi0*s, 0);
      break;
    case 2: // d / d theta
      p.setValues (0, 0, - s/sin2theta);          
      break;
    case 3: // d / d dca
      p.setValues (sphi0, -cphi0, 0);          
      break;
    case 4: // d / d z0
      p.setValues (0, 0, 1);  
      break;
    case 5: // d / d mass
      p.setValues (0, 0, 0);  
      break;
    case 6: // d / d s 
            // fall through 
    case 7: // d / d s
        p.setValues (cphi0, sphi0, cottheta);   
      break;
    default: // should never happen!
      assert (0);
  }
  p *= parfact[ilocal];
}

void NeutralParticleTrack::getVertexEx (int ivertex, ThreeVector& p) const {
  assert (ivertex >= 0 && 6+ivertex < NPAR);
  getTrajectoryPointEx (par[6+ivertex]*parfact[6+ivertex], p);
} 

void NeutralParticleTrack::setVertex (int ivertex, const TwoVector& p) {
  assert (ivertex >= 0 && 6+ivertex < NPAR);
  if (!cachevalid) updateCache();
  // Get point of closest approach to origin
  double xc = sphi0*dca;
  double yc = -cphi0*dca;
  par[6+ivertex] = ((p.getX()-xc)*cphi0 + (p.getY()-yc)*sphi0)/parfact[6+ivertex];
}

void NeutralParticleTrack::getVertexDerivativeEx (int ivertex, 
                                                  int ilocal, 
                                                  ThreeVector& p) const {
  assert (ivertex >= 0 && 6+ivertex < NPAR);
  assert (ilocal >= 0 && ilocal < NPAR);
  // Check: Parameters 5 and 6 are sstart ans sstop.
  // Start vertex does not depend on par[7], 
  // stop vertex does not depend on par[6]
  if (ilocal >= 6 && ilocal != 6+ivertex) {
    p.setValues (0, 0, 0);
  }
  else {
    getTrajectoryDerivativeEx (par[6+ivertex]*parfact[6+ivertex], ilocal, p);
  }
} 

void NeutralParticleTrack::getMomentumAtTrajectoryEx (double s, FourVector& p) const {
  if (!cachevalid) updateCache();
  p.setValues (energy,
               px,
               py,
               pt*cottheta);
}

void NeutralParticleTrack::getMomentumEx (int ivertex, FourVector& p) const {
  assert (ivertex >= 0 && 6+ivertex < NPAR);
  getMomentumAtTrajectoryEx (par[6+ivertex]*parfact[6+ivertex], p);
}


double NeutralParticleTrack::getCharge () const {
  return 0;
} 

double NeutralParticleTrack::getMass () const {
  return  par[5]*parfact[5];
} 

void NeutralParticleTrack::getMomentumDerivativeAtTrajectoryEx (double s, int ilocal, FourVector& p) const {
  if (!cachevalid) updateCache();
  assert (ilocal >= 0 && ilocal < getNPar() && ilocal < NPAR);
  switch (ilocal) {
    case 0: // d / d pt
      {
        p.setValues (beta/sintheta, cphi0, sphi0, cottheta);
      }
      break;
    case 1: // d / d phi0
      p.setValues (0, -py, px, 0);
      break;
    case 2: // d / d theta
      p.setValues (beta*momentum*cottheta,0, 0, -pt/sin2theta);          
      break;
    case 3: // d / d dca
      p.setValues (0, 0, 0, 0);          
      break;
    case 4: // d / d z0
      p.setValues (0, 0, 0, 0);  
      break;
    case 5: // d / d mass
      p.setValues (par[5]/energy, 0, 0, 0);  
      break;
    case 6: // d / d s 
            // fall through 
    case 7: // d / d s
      p.setValues (0, 0, 0, 0);  
      break;
    default: // should never happen!
      assert (0);
  }
  p *= parfact[ilocal];
}

void NeutralParticleTrack::getMomentumDerivativeEx (int ivertex, 
                                                int ilocal, 
                                                FourVector& p) const {
  assert (ivertex >= 0 && 6+ivertex < NPAR);
  assert (ilocal >= 0 && ilocal < NPAR);
  // Check: Parameters 5 and 6 are sstart and sstop.
  // Momentum does not depend on them!
  if (ilocal >= 6) {
    p.setValues (0, 0, 0, 0);
  }
  else {
    getMomentumDerivativeAtTrajectoryEx (par[6+ivertex]*parfact[6+ivertex], ilocal, p);
  }
} 


double NeutralParticleTrack::getArcLength (int ivertex) const {
  assert (ivertex >= 0 && 6+ivertex < NPAR);
  return par[6+ivertex]*parfact[6+ivertex];
}

void NeutralParticleTrack::addToGlobCov(double *glcov, int idim) const {
  int globalnum[NPARMAX];
  bool ok [NPARMAX];
  for (int ilocal = 0; ilocal < getNPar(); ilocal++) {
    int iglobal = globalnum[ilocal] = getGlobalParNum(ilocal);
    if (ok [ilocal] = (iglobal >= 0 && !isParamFixed(ilocal) && isParamMeasured (ilocal))) {
      for (int jlocal = 0; jlocal <= ilocal; jlocal++) {
        int jglobal = globalnum[jlocal];
        if (ok [jlocal]) {
          double c = cov[ilocal][jlocal];
          glcov[jglobal+iglobal*idim] += c;
          if (jlocal != ilocal) glcov[iglobal+jglobal*idim] += c;
        }
      }
    }
  }
} 



void NeutralParticleTrack::updateCache() const {
  pt    = par[0]*parfact[0];
  phi0  = par[1]*parfact[1];
  theta = par[2]*parfact[2];
  dca   = par[3]*parfact[3];
  z0    = par[4]*parfact[4];
  mass  = par[5]*parfact[5];
  
  sphi0 = sin(phi0);
  cphi0 = cos(phi0);
  px =  pt*cphi0;
  py =  pt*sphi0;
  x0 =  dca*sphi0;
  y0 = -dca*cphi0;
  
  sintheta = sin(theta);
  cottheta = cos(theta)/sintheta;
  sin2theta = sintheta*sintheta;
  
  momentum = pt/sintheta;
     
  energy = sqrt(momentum*momentum + mass*mass);
  beta =   (energy != 0) ? momentum/energy : 0;      
    
  calculateChi2();
  cachevalid = true;
 
}

void NeutralParticleTrack::initCov(const float cov_[15]) {
  TrackFitObject::initCov();

  int i = 0;
  for (int k = 0; k < 5; ++k) {
    for (int l = k; l < 5; ++l) {
      cov[k][l] = cov[l][k] = cov_[i++]/(parfact[k]*parfact[l]);
    }
  }                         
  assert (i==15);
  TrackFitObject::checkCov();
 
  covinvvalid = false;
}

JBLHelix NeutralParticleTrack::getTangentialHelix (double) {
  return JBLHelix (0, 
                   par[1]*parfact[1], 
                   par[2]*parfact[2], 
                   par[3]*parfact[3], 
                   par[4]*parfact[4]);
}
