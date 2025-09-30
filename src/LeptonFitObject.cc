/*! \file 
 *  \brief Implements class LeptonFitObject
 *  LeptonFitObject works similiar to JetFitObject, but it uses a q/pt, theta, phi parametrization for the 
 *  leptons, which is e.g. more appropriate for muons and other particles reconstructed from track helices.
 *  The covariance matrix differs from the common E, theta, phi parametrization.
 *  Updated so that q/pt (ptinv) includes the sign of the geometric curvature (scaled version of Omega).
 */ 

#include "LeptonFitObject.h"
#include "EVENT/Track.h"
#include "lcio.h"
#include <cmath>

#undef NDEBUG
#include <cassert>

#include <iostream>

using std::sqrt;
using std::sin;
using std::cos;
using std::cout; 
using std::endl;

using namespace lcio;

// constructor
LeptonFitObject::LeptonFitObject(double ptinv, double theta, double phi,  
                           double Dptinv, double Dtheta, double Dphi, 
                           double m) 
  : ctheta(0), stheta(0), stheta2(0), cphi(0), sphi(0), cottheta(0),
    p2(0), p(0), e(0), e2(0), pt(0), pt2(0), pt3(0), px(0), py(0), pz(0), dpdptinv(0), dpdtheta(0), dptdptinv(0),
    dpxdptinv(0), dpydptinv(0), dpzdptinv(0), dpxdtheta(0), dpydtheta(0), dpzdtheta(0), dpxdphi(0), dpydphi(0), dpzdphi(0),
    chi2(0), dEdptinv(0), dEdtheta(0), dEdp(0), qsign(0), ptinv2(0)
{

  assert( int(NPAR) <= int(BaseDefs::MAXPAR) );

  initCov();                         
  setMass (m);
  adjustPtinvThetaPhi (m, ptinv, theta, phi);
  setParam (0, ptinv, true);
  setParam (1, theta, true);
  setParam (2, phi, true);
  setMParam (0, ptinv);
  setMParam (1, theta);
  setMParam (2, phi);
  setError (0, Dptinv);
  setError (1, Dtheta);
  setError (2, Dphi);

  // parameter 2 repeats every 2*pi
  paramCycl[2]=2.*M_PI;

  invalidateCache();
}

// extended constructor
LeptonFitObject::LeptonFitObject(double ptinv, double theta, double phi,  
				 double Dptinv, double Dtheta, double Dphi,
				 double Rhoptinvtheta, double Rhoptinvphi, double Rhothetaphi, 
				 double m)  
  : ctheta(0), stheta(0), stheta2(0), cphi(0), sphi(0), cottheta(0),
    p2(0), p(0), e(0), e2(0), pt(0), pt2(0), pt3(0), px(0), py(0), pz(0), dpdptinv(0), dpdtheta(0), dptdptinv(0),
    dpxdptinv(0), dpydptinv(0), dpzdptinv(0), dpxdtheta(0), dpydtheta(0), dpzdtheta(0), dpxdphi(0), dpydphi(0), dpzdphi(0),
    chi2(0), dEdptinv(0), dEdtheta(0), dEdp(0), qsign(0), ptinv2(0)
{

  assert( int(NPAR) <= int(BaseDefs::MAXPAR) );

  initCov();                         
  setMass (m);
  adjustPtinvThetaPhi (m, ptinv, theta, phi);
  setParam (0, ptinv, true);
  setParam (1, theta, true);
  setParam (2, phi, true);
  setMParam (0, ptinv);
  setMParam (1, theta);
  setMParam (2, phi);
  setError (0, Dptinv);
  setError (1, Dtheta);
  setError (2, Dphi);
  setCov (0, 1, Rhoptinvtheta*Dptinv*Dtheta);
  setCov (0, 2, Rhoptinvphi*Dptinv*Dphi);
  setCov (1, 2, Rhothetaphi*Dtheta*Dphi);

  // parameter 2 repeats every 2*pi
  paramCycl[2]=2.*M_PI;

  invalidateCache();
}

// constructor based on Track
LeptonFitObject::LeptonFitObject(Track* track, double Bfield, double m) 
  : ctheta(0), stheta(0), stheta2(0), cphi(0), sphi(0), cottheta(0),
    p2(0), p(0), e(0), e2(0), pt(0), pt2(0), pt3(0), px(0), py(0), pz(0), dpdptinv(0), dpdtheta(0), dptdptinv(0),
    dpxdptinv(0), dpydptinv(0), dpzdptinv(0), dpxdtheta(0), dpydtheta(0), dpzdtheta(0), dpxdphi(0), dpydphi(0), dpzdphi(0),
    chi2(0), dEdptinv(0), dEdtheta(0), dEdp(0), qsign(0), ptinv2(0)
{

  assert( int(NPAR) <= int(BaseDefs::MAXPAR) );

  const double c = 2.99792458e8; // m*s^-1
//  const double Bfield = 3.5;          // Tesla       should not be hard-coded here
  const double mm2m = 1e-3;
  const double eV2GeV = 1e-9;
  const double eB = Bfield*c*mm2m*eV2GeV;

  double omega = track->getOmega();
  double ptinv = omega/eB;                   // signed q/pT in GeV^-1
  double tanl = track->getTanLambda();
  double theta = std::atan(1.0/tanl);  
  if (theta<0.0) theta += M_PI;
  double phi = track->getPhi();

  double d3 = 1.0/eB;                        // d(ptinv)/dOmega
  double d5 = -(1.0/(1.0+tanl*tanl));        // d(theta)/d(tanl)  

  FloatVec covT(15, 0.0);
  covT = track->getCovMatrix();

  initCov();                         
  setMass (m);
  adjustPtinvThetaPhi (m, ptinv, theta, phi);
  setParam (0, ptinv, true);
  setParam (1, theta, true);
  setParam (2, phi, true);
  setMParam (0, ptinv);
  setMParam (1, theta);
  setMParam (2, phi);

  setError (0, d3*std::sqrt(covT[5]) );
  setError (1, d5*std::sqrt(covT[14]) );
  setError (2, std::sqrt(covT[2]) );
  setCov (0, 1, d3*d5*covT[12] );            
  setCov (0, 2, d3*covT[4] );
  setCov (1, 2, d5*covT[11] );

  // parameter 2 repeats every 2*pi
  paramCycl[2]=2.*M_PI;

  invalidateCache();
}

// constructor based on TrackState
LeptonFitObject::LeptonFitObject(const TrackState* trackstate, double Bfield, double m) 
  : ctheta(0), stheta(0), stheta2(0), cphi(0), sphi(0), cottheta(0),
    p2(0), p(0), e(0), e2(0), pt(0), pt2(0), pt3(0), px(0), py(0), pz(0), dpdptinv(0), dpdtheta(0), dptdptinv(0),
    dpxdptinv(0), dpydptinv(0), dpzdptinv(0), dpxdtheta(0), dpydtheta(0), dpzdtheta(0), dpxdphi(0), dpydphi(0), dpzdphi(0),
    chi2(0), dEdptinv(0), dEdtheta(0), dEdp(0), qsign(0), ptinv2(0)
{

  assert( int(NPAR) <= int(BaseDefs::MAXPAR) );

  const double c = 2.99792458e8; // m*s^-1
//  const double Bfield = 3.5;          // Tesla       should not be hard-coded here
  const double mm2m = 1e-3;
  const double eV2GeV = 1e-9;
  const double eB = Bfield*c*mm2m*eV2GeV;

  double omega = trackstate->getOmega();
  double ptinv = omega/eB;                   // signed q/pT in GeV^-1
  double tanl = trackstate->getTanLambda();
  double theta = std::atan(1.0/tanl);  
  if (theta<0.0) theta += M_PI;
  double phi = trackstate->getPhi();

  double d3 = 1.0/eB;                        // d(ptinv)/dOmega
  double d5 = -(1.0/(1.0+tanl*tanl));        // d(theta)/d(tanl)  

  FloatVec covT(15, 0.0);
  covT = trackstate->getCovMatrix();

  initCov();                         
  setMass (m);
  adjustPtinvThetaPhi (m, ptinv, theta, phi);
  setParam (0, ptinv, true);
  setParam (1, theta, true);
  setParam (2, phi, true);
  setMParam (0, ptinv);
  setMParam (1, theta);
  setMParam (2, phi);

  setError (0, d3*std::sqrt(covT[5]) );
  setError (1, d5*std::sqrt(covT[14]) );
  setError (2, std::sqrt(covT[2]) );
  setCov (0, 1, d3*d5*covT[12] );            
  setCov (0, 2, d3*covT[4] );
  setCov (1, 2, d5*covT[11] );

  // parameter 2 repeats every 2*pi
  paramCycl[2]=2.*M_PI;

  invalidateCache();
}


// destructor
LeptonFitObject::~LeptonFitObject() {}


// We get a warning that ParticleFitObject should be explicitly initialized
// here, but I don't want to change this part because, I think everything is
// done properly already and not changing behavior is more important.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra"
LeptonFitObject::LeptonFitObject (const LeptonFitObject& rhs)
  : ctheta(0), stheta(0), stheta2(0), cphi(0), sphi(0), cottheta(0),
    p2(0), p(0), e(0), e2(0), pt(0), pt2(0), pt3(0), px(0), py(0), pz(0), dpdptinv(0), dpdtheta(0), dptdptinv(0),
    dpxdptinv(0), dpydptinv(0), dpzdptinv(0), dpxdtheta(0), dpydtheta(0), dpzdtheta(0), dpxdphi(0), dpydphi(0), dpzdphi(0),
    chi2(0), dEdptinv(0), dEdtheta(0), dEdp(0), qsign(0), ptinv2(0)
{
  //std::cout << "copying LeptonFitObject with name" << rhs.name << std::endl;
  LeptonFitObject::assign (rhs);
}
#pragma GCC diagnostic pop

LeptonFitObject& LeptonFitObject::operator= (const LeptonFitObject& rhs) {
  if (this != &rhs) {
    assign (rhs); // calls virtual function assign of derived class
  }
  return *this;
}

LeptonFitObject *LeptonFitObject::copy() const {
  return new LeptonFitObject (*this);
}
    
LeptonFitObject& LeptonFitObject::assign (const BaseFitObject& source) {
  if (const LeptonFitObject *psource = dynamic_cast<const LeptonFitObject *>(&source)) {
    if (psource != this) {
      ParticleFitObject::assign (source);
      // only mutable data members, need not to be copied, if cache is invalid
    }
  }
  else {
    assert (0);
  }
  return *this;
}

const char *LeptonFitObject::getParamName (int ilocal) const {
  switch (ilocal) {
    case 0: return "1/Pt";
    case 1: return "theta";
    case 2: return "phi";
  }
  return "undefined";
}

 
bool LeptonFitObject::updateParams (double pp[], int idim) {
  invalidateCache();
  
  int iptinv = getGlobalParNum(0);
  int ith    = getGlobalParNum(1);
  int iph    = getGlobalParNum(2);
  assert (iptinv >= 0 && iptinv  < idim);
  assert (ith    >= 0 && ith     < idim);
  assert (iph    >= 0 && iph     < idim);
  
//  double ptinv  = std::abs(p[iptinv]);
  double ptinv  = pp[iptinv];                    // can be +ve or -ve.
  double th = pp[ith];
  double ph = pp[iph];
  
  bool result = ((ptinv -par[0])*(ptinv -par[0]) > eps2*cov[0][0]) ||
                ((th-par[1])*(th-par[1]) > eps2*cov[1][1]) ||
                ((ph-par[2])*(ph-par[2]) > eps2*cov[2][2]);
                
  par[0] = ptinv;
  par[1] = th;
  par[2] = ph;
  pp[iptinv] = par[0];         
  pp[ith]    = par[1];         
  pp[iph]    = par[2]; 

  // std::cout << "GWW hello from LeptonFitObject::updateParams()" << par[0] << " " << par[1] << " " << par[2] << std::endl;
        
  return result;
}  

// these depend on actual parametrisation!

double LeptonFitObject::getDPx(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpxdptinv;
    case 1: return 0;            // dpxdtheta = 0
    case 2: return dpxdphi;
  }
  return 0; 
}

double LeptonFitObject::getDPy(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpydptinv;
    case 1: return 0;            // dpydtheta = 0
    case 2: return dpydphi;
  }
  return 0; 
}

double LeptonFitObject::getDPz(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpzdptinv; 
    case 1: return dpzdtheta;
    case 2: return 0;            // dpzdphi = 0
  }
  return 0; 
}

double LeptonFitObject::getDE(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR); 
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dEdptinv;
    case 1: return dEdtheta;
    case 2: return 0;            // dEdphi = 0
  }
  return 0; 
}

double LeptonFitObject::getFirstDerivative_Meta_Local( int iMeta, int ilocal , int metaSet ) const {
  // iMeta = intermediate variable (i.e. E,px,py,pz)
  // ilocal = local variable (ptinv, theta, phi)
  // metaSet = which set of intermediate varlables

  assert (metaSet==0); // only defined for E,px,py,pz  

  switch ( iMeta ) {
  case 0: // E
    return getDE(ilocal);
    break;
  case 1: // Px
    return getDPx(ilocal);
    break;
  case 2: // Py
    return getDPy(ilocal);
    break;
  case 3: // Pz
    return getDPz(ilocal);
    break;
  default:
    assert(0);
  }
  return -999;
}

double LeptonFitObject::getSecondDerivative_Meta_Local( int iMeta, int ilocal, int jlocal, int metaSet ) const {
  assert ( metaSet==0 );
  if (!cachevalid) updateCache();
  if ( jlocal<ilocal ) {
    int temp=jlocal; 
    jlocal=ilocal;   
    ilocal=temp;
  }

  // Evaluated and checked with Mathematica using (px, py, pz) = q/k [cos(phi), sin(phi), cot(theta)]
  // with local parameters of (k, theta, phi)
  // q is the sign of the track curvature (+-1) and k the signed 1/pT. So pT = q/k for either charge.
  // The four non-zero off-diagonal elements are sensitive to the charge sign (qsign).
  // Graham Wilson

  switch (iMeta) {
  case 0: // E
    if      ( ilocal==0 && jlocal==0 ) return (2*e2+mass*mass)*pt2*p*p/(e2*e);
    else if ( ilocal==0 && jlocal==1 ) return qsign*pt*p2*cottheta*(e2+mass*mass)/(e2*e);
    else if ( ilocal==1 && jlocal==1 ) return (p2/(e*stheta2))*(1.0+ctheta*ctheta*(1.0+(mass*mass/e2)));
    else return 0;
    break;
  case 1: // px
    if      ( ilocal==0 && jlocal==0 ) return 2*pt3*cphi;
    else if ( ilocal==0 && jlocal==2 ) return qsign*pt2*sphi;
    else if ( ilocal==2 && jlocal==2 ) return -px;
    else return 0;
    break;
  case 2: // py
    if      ( ilocal==0 && jlocal==0 ) return 2*pt3*sphi;
    else if ( ilocal==0 && jlocal==2 ) return -qsign*pt2*cphi;
    else if ( ilocal==2 && jlocal==2 ) return -py;
    else return 0;
    break;
  case 3: // pz
    if      ( ilocal==0 && jlocal==0 ) return 2*pt3*cottheta;
    else if ( ilocal==0 && jlocal==1 ) return qsign*pt2/stheta2;
    else if ( ilocal==1 && jlocal==1 ) return 2*pt*cottheta/stheta2;
    else return 0;
    break;
  default:
    assert(0);
  }
  return -999;
}         

void LeptonFitObject::updateCache() const {
  // std::cout << "LeptonFitObject::updateCache" << std::endl;
  //  chi2 = getChi2 ();
  
  double ptinv = par[0];
  double theta = par[1];
  double phi   = par[2];
  
  qsign = 1.0;
  if(ptinv<0.0)qsign = -1.0;

  ptinv2 = ptinv*ptinv;

//  pt = 1/std::abs(ptinv);
  pt = qsign/ptinv;
  pt2 = pt*pt;
  pt3 = pt2*pt;

  ctheta = cos(theta);
  stheta = sin(theta);
  stheta2 = stheta*stheta;
  cphi   = cos(phi);
  sphi   = sin(phi);
  cottheta = ctheta/stheta;
 
  p = pt/stheta;
  assert(p!=0);
  
  p2 = p*p;
  e2 = p2+mass*mass;
  e = std::sqrt(e2);
  px = pt*cphi;
  py = pt*sphi;
  pz = pt*cottheta;

  fourMomentum.setValues( e, px, py, pz );
  
// Need to take care of charge sign for the d/dptinv terms.    Graham
  dpdptinv = -qsign/(ptinv2*stheta);
  dpxdptinv = -qsign*cphi/ptinv2;
  dpydptinv = -qsign*sphi/ptinv2;
  dpzdptinv = -qsign*cottheta/ptinv2;

  dpdtheta = (qsign/ptinv)*(-cottheta/stheta);
  dpzdtheta = (qsign/ptinv)*(-1.0/stheta2);

  dpxdphi   = (qsign/ptinv)*(-sphi);
  dpydphi   = (qsign/ptinv)*cphi;

  dEdp      = p/e;
  dEdptinv  = dEdp*dpdptinv;
  dEdtheta  = dEdp*dpdtheta;
 
  cachevalid = true;
}

//double LeptonFitObject::getChi2 () const {
//  if (!cachevalid) updateCache();
//  return chi2;
//}

bool LeptonFitObject::adjustPtinvThetaPhi (double&, double &, double& theta, double& phi) {
  bool result = false;
  
/*  Keep the sign information - off-diagonal terms of the error matrix care about this ...  Graham
  if (ptinv<0) {
    // cout << "LeptonFitObject::adjustPtinvThetaPhi: mirrored E!\n";
    ptinv  = -ptinv;
    result = true;
  }
*/
  if (theta < -M_PI || theta > M_PI) {
    while (theta < -M_PI) theta += 2*M_PI;
    while (theta >  M_PI) theta -= 2*M_PI;
    result = true;
  }
  
  if (theta<0) {
    // cout << "LeptonFitObject::adjustPtinvThetaPhi: mirrored theta!\n";
    theta = -theta;
    phi = phi > 0 ? phi-M_PI : phi+M_PI;
    result = true;
  }
  else if (theta>M_PI) {
    // cout << "LeptonFitObject::adjustPtinvThetaPhi: mirrored theta!\n";
    theta = 2*M_PI-theta;
    phi = phi > 0 ? phi-M_PI : phi+M_PI;
    result = true;
  }
  if (phi < -M_PI || phi > M_PI) {
    while (phi < -M_PI) phi += 2*M_PI;
    while (phi >  M_PI) phi -= 2*M_PI;
    result = true;
  }

  return result;
}
