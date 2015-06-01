/*! \file 
 *  \brief Implements class LeptonFitObject
 *  LeptonFitObject works similiar to JetFitObject, but it uses a 1/pt, eta, phi parametrization for the 
 *  leptons, which is e.g. more appropriate for electrons.
 *  Especially the covarianz matrix differs from a common E, theta, phi parametrization.
 */ 

#include "LeptonFitObject.h"
#include <cmath>
#include <cassert>
#include <iostream>

using std::sqrt;
using std::sin;
using std::cos;
using std::cout; 
using std::endl;

// constructor
LeptonFitObject::LeptonFitObject(double ptinv, double theta, double phi,  
                           double Dptinv, double Dtheta, double Dphi, 
                           double m) {

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

// destructor
LeptonFitObject::~LeptonFitObject() {}

LeptonFitObject::LeptonFitObject (const LeptonFitObject& rhs)
{
  //std::cout << "copying LeptonFitObject with name" << rhs.name << std::endl;
  LeptonFitObject::assign (rhs);
}

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

 
bool LeptonFitObject::updateParams (double p[], int idim) {
  invalidateCache();
  
  int iptinv = getGlobalParNum(0);
  int ith    = getGlobalParNum(1);
  int iph    = getGlobalParNum(2);
  assert (iptinv >= 0 && iptinv  < idim);
  assert (ith    >= 0 && ith     < idim);
  assert (iph    >= 0 && iph     < idim);
  
  double ptinv  = std::abs(p[iptinv]);
  double th = p[ith];
  double ph = p[iph];
  
  bool result = ((ptinv -par[0])*(ptinv -par[0]) > eps2*cov[0][0]) ||
                ((th-par[1])*(th-par[1]) > eps2*cov[1][1]) ||
                ((ph-par[2])*(ph-par[2]) > eps2*cov[2][2]);
                
  par[0] = ptinv;
  par[1] = th;
  par[2] = ph;
  p[iptinv] = par[0];         
  p[ith]    = par[1];         
  p[iph]    = par[2];         
  return result;
}  

// these depend on actual parametrisation!

double LeptonFitObject::getDPx(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpxdptinv;
      //case 1: return dpxdtheta; 
    case 1: return 0;
    case 2: return dpxdphi;
  }
  return 0; 
}

double LeptonFitObject::getDPy(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpydptinv;
      // case 1: return dpydtheta;
    case 1: return 0;
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
    case 2: return 0;
  }
  return 0; 
}

double LeptonFitObject::getDE(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR); 
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dEdptinv;
    case 1: return dEdtheta;
    case 2: return 0;
  }
  return 0; 
}

double LeptonFitObject::getFirstDerivative( int iMeta, int ilocal , int metaSet ) const {
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
}

double LeptonFitObject::getSecondDerivative( int iMeta, int ilocal , int jlocal , int metaSet ) const {
  assert ( metaSet==0 );
  if (!cachevalid) updateCache();
  if ( jlocal<ilocal ) {
    int temp=jlocal;
    ilocal=jlocal;
    ilocal=temp;
  }

  // daniel hasn't checked these, copied from orig code
  double d2Edp2 = mass*mass/(e2*e);
  switch (iMeta) {
  case 0: // E
    if      ( ilocal==0 && jlocal==0 ) return (2.*pt3/stheta)*dEdp + dpdptinv*dpdptinv*d2Edp2;
    else if ( ilocal==0 && jlocal==1 ) return pt2*cottheta/stheta*dEdp + dpdptinv*dpdtheta*d2Edp2;
    else if ( ilocal==1 && jlocal==1 ) return pt/stheta*(2*cottheta*cottheta + 1);
    else return 0;
    break;
  case 1: // px
    if      ( ilocal==0 && jlocal==0 ) return 2*pt3*cphi;
    else if ( ilocal==0 && jlocal==2 ) return pt2*sphi;
    //    else if ( ilocal==2 && jlocal==2 ) return px; // "-" in orig code, DJ fixed 2015may27
    else if ( ilocal==2 && jlocal==2 ) return -px;
    else return 0;
    break;
  case 2: // py
    if      ( ilocal==0 && jlocal==0 ) return 2*pt3*sphi;
    else if ( ilocal==0 && jlocal==2 ) return pt2*cphi;
    else if ( ilocal==2 && jlocal==2 ) return -py;
    else return 0;
    break;
  case 3: // pz
    if      ( ilocal==0 && jlocal==0 ) return 2*pt3*cottheta;
    else if ( ilocal==0 && jlocal==1 ) return pt2/stheta2;
    else if ( ilocal==1 && jlocal==1 ) return 2*pt*cottheta/stheta2;
    else return 0;
    break;
  default:
    assert(0);
  }
}
         

void LeptonFitObject::updateCache() const {
  // std::cout << "LeptonFitObject::updateCache" << std::endl;
  //  chi2 = getChi2 ();
  
  double ptinv = par[0];
  double theta = par[1];
  double phi   = par[2];
  
  pt = 1/ptinv;
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
  px = cphi*pt;
  py = sphi*pt;
  pz = cottheta*pt;

  fourMomentum.setValues( e, px, py, pz );
  
  dpdptinv = -pt2/stheta;
  dpxdptinv = -cphi*pt2;
  dpydptinv = -sphi*pt2;
  dpzdptinv = -cottheta*pt2;

  dpdtheta = -(cottheta*pt/stheta);
  dpzdtheta = -pt/stheta2;

  dpxdphi   = -sphi*pt;
  dpydphi   =  cphi*pt;

  dEdp      = p/e;
  dEdptinv  = dpdptinv*dEdp;
  dEdtheta  = dpdtheta*dEdp;
 
  cachevalid = true;
}

//double LeptonFitObject::getChi2 () const {
//  if (!cachevalid) updateCache();
//  return chi2;
//}

bool LeptonFitObject::adjustPtinvThetaPhi (double& m, double &ptinv, double& theta, double& phi) {
  bool result = false;
  
  if (ptinv<0) {
    // cout << "LeptonFitObject::adjustEThetaPhi: mirrored E!\n";
    ptinv  = -ptinv;
    result = true;
  }
  if (theta < -M_PI || theta > M_PI) {
    while (theta < -M_PI) theta += 2*M_PI;
    while (theta >  M_PI) theta -= 2*M_PI;
    result = true;
  }
  
  if (theta<0) {
    // cout << "LeptonFitObject::adjustEThetaPhi: mirrored theta!\n";
    theta = -theta;
    phi = phi > 0 ? phi-M_PI : phi+M_PI;
    result = true;
  }
  else if (theta>M_PI) {
    // cout << "LeptonFitObject::adjustEThetaPhi: mirrored theta!\n";
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


