/*! \file 
 *  \brief Implements class ChargedParticleTrack
 *
 * \b Changelog:
 * - 17.11.04 BL First version
 *
 */ 
#include "ChargedParticleTrack.h"
#include "TwoVector.h"
#include "ThreeVector.h"
#include "FourVector.h"
#include "JBLHelix.h"

#include <cassert>
#include <cstring>
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

const double ChargedParticleTrack::parfact[NPARMAX] = 
  {0.0001, 0.001, 0.01, 0.01, 1., 1., 1., 1.};

ChargedParticleTrack::ChargedParticleTrack (const char *name_,
                                            double kappa,
                                            double phi0,
                                            double theta,
                                            double dca,
                                            double z0,
                                            double mass,
                                            double charge_,
                                            double sstart,
                                            double sstop
                                            )
: TrackFitObject (name_), charge (abs(charge_))                                          
{
  setParam (0, kappa/parfact[0], true);
  setParam (1, phi0/parfact[1], true);
  setParam (2, theta/parfact[2], true);
  setParam (3, dca/parfact[3], true);
  setParam (4, z0/parfact[4], true);
  setParam (5, mass/parfact[5], false, true);
  setParam (6, sstart/parfact[6], false);
  setParam (7, sstop/parfact[7], false);
  setMParam (0, kappa/parfact[0]);
  setMParam (1, phi0/parfact[1]);
  setMParam (2, theta/parfact[2]);
  setMParam (3, dca/parfact[3]);
  setMParam (4, z0/parfact[4]);
  TrackFitObject::initCov();
}
ChargedParticleTrack::ChargedParticleTrack (const char *name_,
                                            double kappa,
                                            double phi0,
                                            double theta,
                                            double dca,
                                            double z0,
                                            double mass,
                                            double charge_,
                                            double sstart
                                            )
: TrackFitObject (name_), charge (abs(charge_))                                          
{
  setParam (0, kappa/parfact[0], true);
  setParam (1, phi0/parfact[1], true);
  setParam (2, theta/parfact[2], true);
  setParam (3, dca/parfact[3], true);
  setParam (4, z0/parfact[4], true);
  setParam (5, mass/parfact[5], false, true);
  setParam (6, sstart/parfact[6], false);
  setParam (7, 0, false, true);
  setMParam (0, kappa/parfact[0]);
  setMParam (1, phi0/parfact[1]);
  setMParam (2, theta/parfact[2]);
  setMParam (3, dca/parfact[3]);
  setMParam (4, z0/parfact[4]);
  TrackFitObject::initCov();
}
ChargedParticleTrack::ChargedParticleTrack (const char *name_,
                                            const ThreeVector& vertex,
                                            const ThreeVector& momentum,
                                            double mass,
                                            double charge_
                                           )
: TrackFitObject (name_), charge (abs (charge_))                                          
{
  setParameters (0, vertex, FourVector (momentum, mass), charge_);
  assert (charge != 0);
  assert (bfield != 0);
  cBq = charge_*0.002998*bfield;    // charge_ is signed, charge is not!
  kappa = -cBq/momentum.getPt();
  r = 1/kappa;
  theta = momentum.getTheta();
  double phi0plkappas = momentum.getPhi();
//  double phi0plkappas = atan2 (momentum.getPy(), momentum.getPx());
  double si = sin(phi0plkappas);
  double co = cos(phi0plkappas);
  
  double xc = vertex.getX()-r*si;
  double yc = vertex.getY()+r*co;
  
  // rc is the distance of the circle's center to the origin
  double rc = sqrt(xc*xc + yc*yc);
  dca = r*(1 - abs(kappa)*rc); 
  
  if (dca-r > 0) 
    phi0 = atan2 (xc, -yc);
  else
    phi0 = atan2 (-xc, yc);
  
//   cout << "ChargedParticleTrack::ChargedParticleTrack"
//        << "\nvertex = " << vertex << ", momentum = " << momentum
//        << "\ncBq=" << cBq << ", kappa=" << kappa << ", r=" << r << ", theta=" << theta 
//        << "\nphi0plkappas=" << phi0plkappas << ", si=" << si << ", co=" << co << ", phi0=" << phi0 
//        << "\nxc=" << xc << ", yc=" << yc 
//        << endl;
       
  // if (phi0plkappas - phi0 < -1.57) phi0plkappas += 6.2831853;
  double sstart = (phi0plkappas-phi0)/kappa;
//   cout << "sstart = " << sstart;
  if (sstart < -1.57*abs(r)) sstart += 6.2831853*abs(r);
//   cout << " => sstart = " << sstart << endl;
    
  double z0 = vertex.getZ()  - sstart*(cos(theta)/sin(theta));

//   cout << "rc = " << rc << ", dca = " << dca << ", z0 = " << z0 << endl << endl;

  setParam (0, kappa/parfact[0], true);
  setParam (1, phi0/parfact[1], true);
  setParam (2, theta/parfact[2], true);
  setParam (3, dca/parfact[3], true);
  setParam (4, z0/parfact[4], true);
  setParam (5, mass/parfact[5], false, true);
  setParam (6, sstart/parfact[6], false);
  setParam (7, 0, false, true);
  setMParam (0, kappa/parfact[0]);
  setMParam (1, phi0/parfact[1]);
  setMParam (2, theta/parfact[2]);
  setMParam (3, dca/parfact[3]);
  setMParam (4, z0/parfact[4]);
  
  TrackFitObject::initCov();
}
ChargedParticleTrack::ChargedParticleTrack (const char *name_,
                                            const float par_[5],
                                            const float cov_[15],
                                            double mass,
                                            double charge_,
                                            double sstart,
                                            double sstop
                                           ) 
: TrackFitObject (name_), charge (abs(charge_))  
{                                        
  setParam  (0, par_[0]/parfact[0], true);
  setParam  (1, par_[1]/parfact[1], true);
  setParam  (2, par_[2]/parfact[2], true);
  setParam  (3, par_[3]/parfact[3], true);
  setParam  (4, par_[4]/parfact[4], true);
  setParam  (5, mass   /parfact[5], false, true);
  setParam  (6, sstart /parfact[6], false);
  setParam  (7, sstop  /parfact[7], false, true);
  setMParam (0, par_[0]/parfact[0]);
  setMParam (1, par_[1]/parfact[1]);
  setMParam (2, par_[2]/parfact[2]);
  setMParam (3, par_[3]/parfact[3]);
  setMParam (4, par_[4]/parfact[4]);
  initCov(cov_);
  calculateCovInv();
}

ChargedParticleTrack *ChargedParticleTrack::copy() const {
  return new ChargedParticleTrack (*this);
}
    
ChargedParticleTrack& ChargedParticleTrack::assign (const BaseFitObject& source) {
  if (const ChargedParticleTrack *psource = dynamic_cast<const ChargedParticleTrack *>(&source)) {
    if (psource != this) *this = *psource;
  }
  else {
    assert (0);
  }
  return *this;
}

ChargedParticleTrack::~ChargedParticleTrack ()
{}

int ChargedParticleTrack::getNPar() const {
  return NPAR;
}

const char *ChargedParticleTrack::getParamName (int ilocal) const {
  switch (ilocal) {
    case 0: return "kappa";
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

bool ChargedParticleTrack::setParam (int ilocal, double par_, 
                                 bool measured_, bool fixed_) {
  assert (ilocal >= 0 && ilocal < NPARMAX);
  if (measured[ilocal] != measured_ || fixed[ilocal] != fixed_) invalidateCache();
  measured[ilocal] = measured_;
  fixed[ilocal] = fixed_;
  return setParam (ilocal, par_);
}  

bool ChargedParticleTrack::setParameters (int ivertex, 
                                          const ThreeVector& vertex,  
                                          const FourVector& momentum,
                                          double charge_    
                                         ) { 
  bool result = true;                                       
  charge = std::abs(charge_);                                       
  assert (charge > 0);
  assert (bfield > 0);
  cBq = 0.002998*bfield*charge_;
  if (momentum.getPt() == 0) return false;
  kappa = -cBq/momentum.getPt();
  r = 1/kappa;
  theta = momentum.getTheta();
  double phi0plkappas = momentum.getPhi();
  
  mass = (isParamFixed(5)) ? par[5]*parfact[5] : momentum.getMass();
  
//  double phi0plkappas = atan2 (momentum.getPy(), momentum.getPx());
  double si = sin(phi0plkappas);
  double co = cos(phi0plkappas);
  
  double xc = vertex.getX()-r*si;
  double yc = vertex.getY()+r*co;
  
  // rc is the distance of the circle's center to the origin
  double rc = sqrt(xc*xc + yc*yc);
  dca = r*(1 - abs(kappa)*rc); 
  
  if (dca-r > 0) 
    phi0 = atan2 (xc, -yc);
  else
    phi0 = atan2 (-xc, yc);
  
//   cout << "ChargedParticleTrack::setParameters"
//        << "\nvertex = " << vertex << ", momentum = " << momentum
//        << "\ncBq=" << cBq << ", kappa=" << kappa << ", r=" << r << ", theta=" << theta 
//        << "\nphi0plkappas=" << phi0plkappas << ", si=" << si << ", co=" << co << ", phi0=" << phi0 
//        << "\nxc=" << xc << ", yc=" << yc 
//        << endl;
       
  // if (phi0plkappas - phi0 < -1.57) phi0plkappas += 6.2831853;
  double s = getNormalS((phi0plkappas-phi0)/kappa);
  z0 = vertex.getZ()  - s*(cos(theta)/sin(theta));

//   cout << "rc = " << rc << ", dca = " << dca << ", z0 = " << z0 << endl << endl;

  
  if (!isParamFixed (0)) result &= setParam (0, kappa/parfact[0]);
  if (!isParamFixed (1)) result &= setParam (1, phi0/parfact[1]);
  if (!isParamFixed (2)) result &= setParam (2, theta/parfact[2]);
  if (!isParamFixed (3)) result &= setParam (3, dca/parfact[3]);
  if (!isParamFixed (4)) result &= setParam (4, z0/parfact[4]);
  if (!isParamFixed (5)) result &= setParam (5, mass/parfact[5]);
  assert (ivertex >= 0 && 6+ivertex < NPAR);
  result &= setParam (6+ivertex, s/parfact[6+ivertex]);
  return result;
}


bool ChargedParticleTrack::setParam (int ilocal, double par_ ) {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!isfinite(par_)) return false;
  // enforce parameter range restrictions
  // dont restrict phi0 to -pi ... pi, because otherwise residual calculation
  // becomes too complicated!
  switch (ilocal) {
   // case 1:  if (abs(par_) > M_PI) par_ = atan2 (sin (par_), cos (par_)); // -pi <= phi0 <= pi
   //   break;
    case 2:  if (par_ < 0 || par_ >= M_PI/parfact[2]) par_ = acos (cos (par_*parfact[2]))/parfact[2];
      break;
    case 5:  if (par_ < 0) par_ = abs (par_);
      break;
  } 
  if (par[ilocal] == par_) return true;
  invalidateCache();
  par[ilocal] = par_;
  return true;
}  

void ChargedParticleTrack::getTrajectoryPointEx (double s, ThreeVector& p) const {
  if (!cachevalid) updateCache();
  
  if (abs (kappa*s) < 1.e-6) {
    double ssq = s*s;
    double kssq = 0.5*kappa*ssq;
    double dcamikssq = dca - kssq;
    p.setValues ( sphi0*dcamikssq + cphi0*s, // bug fixed! BL 4.1.05
                 -cphi0*dcamikssq + sphi0*s,
                  z0 + s*cottheta);
  }
  else {
    double kappas = kappa*s;
    double phi0plkappas = phi0 + kappas;
    double si = sin(phi0plkappas);
    double co = cos(phi0plkappas);
    p.setValues ( dcamir*sphi0 + r*si,
                 -dcamir*cphi0 - r*co,
                 z0 + s*cottheta);
  }
}

void ChargedParticleTrack::getTrajectoryDerivativeEx (double s, int ilocal, ThreeVector& p) const {
  if (!cachevalid) updateCache();
  assert (ilocal >= 0 && ilocal < NPAR);
  switch (ilocal) {
    case 0: // d / d par[0] = d / d kappa * kappafact
      if (abs (kappa*s) < 1.e-6) {
        double ssq = s*s;
        p.setValues (-0.5*sphi0*ssq,
                      0.5*cphi0*ssq, 
                      0);          
      }
      else {
        double kappas = kappa*s;
        double phi0plkappas = phi0 + kappas;
        double si = sin(phi0plkappas);
        double co = cos(phi0plkappas);
        p.setValues ( co*kappas + sphi0 -si,
                      si*kappas - cphi0 +co, 
                      0);          
      }
      break;
    case 1: // d / d phi0
      if (abs (kappa*s) < 1.e-6) {
        double ssq = s*s;
        double kssq = 0.5*kappa*ssq;
        double dcamikssq = dca - kssq;
      p.setValues ( cphi0*dcamikssq - sphi0*s,
                    sphi0*dcamikssq + cphi0*s,
                    0);
      }
      else {
        double kappas = kappa*s;
        double phi0plkappas = phi0 + kappas;
        p.setValues ( dcamir*cphi0 + r*cos(phi0plkappas),
                      dcamir*sphi0 + r*sin(phi0plkappas),
                      0);
      }
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
      {
        double kappas = kappa*s;
        double phi0plkappas = phi0 + kappas;
        p.setValues (cos(phi0plkappas), sin(phi0plkappas), cottheta);   
      }
      break;
    default: // should never happen!
      assert (0);
  }
  p *= parfact[ilocal];
}

void ChargedParticleTrack::getVertexEx (int ivertex, ThreeVector& p) const {
  assert (ivertex >= 0 && 6+ivertex < NPAR);
  getTrajectoryPointEx (par[6+ivertex]*parfact[6+ivertex], p);
} 

void ChargedParticleTrack::setVertex (int ivertex, const TwoVector& p) {
  assert (ivertex >= 0 && 6+ivertex < NPAR);
  if (!cachevalid) updateCache();
  if (std::abs (kappa) < 1.e-8) {
    // If helix is almost a straight line,
    // treat it as straight line here.
    // Get point of closest approach to origin
    double xc = sphi0*dca;
    double yc = -cphi0*dca;
    par[6+ivertex] = ((p.getX()-xc)*cphi0 + (p.getY()-yc)*sphi0)/parfact[6+ivertex];
  }
  else {
    // Get center of circle
    double xc =  dcamir*sphi0;
    double yc = -dcamir*cphi0;
    double psi = (r > 0) ? 
                 std::atan2(p.getX()-xc, -p.getY()+yc) - phi0 :
                 std::atan2(p.getX()-xc,  p.getY()-yc) + phi0;
    par[6+ivertex] =  getNormalS (std::abs(r)*psi)/parfact[6+ivertex];
  }
}

void ChargedParticleTrack::getVertexDerivativeEx (int ivertex, 
                                                int ilocal, 
                                                ThreeVector& p) const {
  assert (ivertex >= 0 && 6+ivertex < NPAR);
  assert (ilocal >= 0 && ilocal < NPAR);
  // Check: Parameters 6 and 7 are sstart and sstop.
  // Start vertex does not depend on par[6], 
  // stop vertex does not depend on par[5]
  if (ilocal >= 6 && ilocal != 6+ivertex) {
    p.setValues (0, 0, 0);
  }
  else {
    getTrajectoryDerivativeEx (par[6+ivertex]*parfact[6+ivertex], ilocal, p);
  }
} 

void ChargedParticleTrack::getMomentumAtTrajectoryEx (double s, FourVector& p) const {
  if (!cachevalid) updateCache();
  double phi0plkappas = phi0 + kappa*s;
  p.setValues (energy,
               pt*cos(phi0plkappas),
               pt*sin(phi0plkappas),
               pt*cottheta);
}

void ChargedParticleTrack::getMomentumEx (int ivertex, FourVector& p) const {
  assert (ivertex >= 0 && 6+ivertex < NPAR);
  getMomentumAtTrajectoryEx (par[6+ivertex]*parfact[6+ivertex], p);
}

double ChargedParticleTrack::getCharge () const {
  return (kappa<0) ? charge : -charge;
} 

double ChargedParticleTrack::getMass () const {
  return  par[5]*parfact[5];
} 

void ChargedParticleTrack::getMomentumDerivativeAtTrajectoryEx (double s, int ilocal, FourVector& p) const {
  if (!cachevalid) updateCache();
  assert (ilocal >= 0 && ilocal < NPAR);
  switch (ilocal) {
    case 0: // d / d par[0] = d / d kappa * kappafact
      {
        double kappas = kappa*s;
        double phi0plkappas = phi0 + kappas;
        double si = sin(phi0plkappas);
        double co = cos(phi0plkappas);
        p.setValues (-beta*momentum/(energy*kappa),
                      momderfact*(co + kappas*si) ,
                      momderfact*(si - kappas*co) ,
                      momderfact*cottheta         );
      }
      break;
    case 1: // d / d phi0
      {
        double kappas = kappa*s;
        double phi0plkappas = phi0 + kappas;
        p.setValues (0, -pt*sin(phi0plkappas), pt*cos(phi0plkappas), 0);
      }
      break;
    case 2: // d / d theta
      p.setValues (-beta*momentum*cottheta, 0, 0, -pt/sin2theta);          
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
      {
        double kappas = kappa*s;
        double phi0plkappas = phi0 + kappas;
        p.setValues ( 0,
                      cBq*sin(phi0plkappas), 
                     -cBq*cos(phi0plkappas), 
                      0);   
      }
      break;
    default: // should never happen!
      assert (0);
  }
  p *= parfact[ilocal];
}

void ChargedParticleTrack::getMomentumDerivativeEx (int ivertex, 
                                                int ilocal, 
                                                FourVector& p) const {
  assert (ivertex >= 0 && 6+ivertex < NPAR);
  assert (ilocal >= 0 && ilocal < NPAR);
  // Check: Parameters 6 and 7 are sstart ans sstop.
  // Start vertex does not depend on par[6], 
  // stop vertex does not depend on par[5]
  if (ilocal >= 6 && ilocal != 6+ivertex) {
    p.setValues (0, 0, 0, 0);
  }
  else {
    getMomentumDerivativeAtTrajectoryEx (par[6+ivertex]*parfact[6+ivertex], ilocal, p);
  }
} 

double ChargedParticleTrack::getArcLength (int ivertex) const {
  assert (ivertex >= 0 && 6+ivertex < NPAR);
  return par[6+ivertex]*parfact[6+ivertex];
}

void ChargedParticleTrack::updateCache() const {
  kappa = par[0]*parfact[0];
  phi0  = par[1]*parfact[1];
  theta = par[2]*parfact[2];
  dca   = par[3]*parfact[3];
  z0    = par[4]*parfact[4];
  mass  = par[5]*parfact[5];
  sphi0 = sin(phi0);
  cphi0 = cos(phi0);
  r = (kappa != 0) ? 1/kappa : 0;
  dcamir = dca - r;
  sintheta = sin(theta);
  cottheta = cos(theta)/sintheta;
  sin2theta = sintheta*sintheta;
  assert (charge > 0);
  assert (bfield > 0);
  cBq = 0.002998*bfield*charge * ((kappa < 0) ? +1 : -1);
  pt = (kappa != 0) ? -cBq/kappa : 0;
  assert (kappa == 0 || pt > 0);
  momderfact = (kappa != 0) ? -pt/kappa : 0;
  
  momentum = (sintheta != 0) ? pt/sintheta : 0;
  // assert (momentum > 0);
  energy = (kappa != 0) ? 
           sqrt(momentum*momentum + mass*mass) : 
           mass;
  beta =   (energy != 0) ? momentum/energy : 0;     
   
  calculateChi2();
  cachevalid = true;
 
}

void ChargedParticleTrack::initCov(const float cov_[15]) {
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


JBLHelix ChargedParticleTrack::getTangentialHelix (double) {
  return JBLHelix (par[0]*parfact[0], 
                   par[1]*parfact[1], 
                   par[2]*parfact[2], 
                   par[3]*parfact[3], 
                   par[4]*parfact[4]);


}

double ChargedParticleTrack::getNormalS (double s) const {
  kappa = par[0]*parfact[0];
  double kappas = par[0]*parfact[0]*s;
  if (kappas >= -M_PI &&kappas < M_PI) return s;
  return std::atan2 (std::sin(kappas), std::cos (kappas))/kappa;
}
