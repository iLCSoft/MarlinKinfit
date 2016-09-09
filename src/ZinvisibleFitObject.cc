/*! \file
 *  \brief Implements clsss ZinvisibleFitObject
 *  class for Z->neutrinos with (E, theta, phi) in kinematic fit
 *  Zinvisible works similar to NeutrinoFitObject, but its mass is set to 91.1876
 *  developed for ZHH->vvHH, ZinvisibleFO represents Z->vv in fit
 *
 */

#include "ZinvisibleFitObject.h"
#include <cmath>

#undef NDEBUG
#include <cassert>

#include <algorithm>

using std::sqrt;
using std::sin;
using std::cos;
using std::cout; 
using std::endl;

// constructor
ZinvisibleFitObject::ZinvisibleFitObject(double E, double theta, double phi, 
					 double DE, double Dtheta, double Dphi, double m) 
  : cachevalid(false), ctheta(0), stheta(0), cphi(0), sphi(0),p2(0), p(0), dpdE(0), pt(0), px(0), py(0), pz(0), dptdE(0),
    dpxdE(0), dpydE(0), dpxdtheta(0), dpydtheta(0), chi2(0)

{  //hier double m

  assert( int(NPAR) <= int(BaseDefs::MAXPAR) );
  setMass (m);  
  setParam (0, E, false);
  setParam (1, theta, false);
  setParam (2, phi, false);
  setError (0, DE);
  setError (1, Dtheta);
  setError (2, Dphi);
  invalidateCache();
}

// destructor
ZinvisibleFitObject::~ZinvisibleFitObject() {}

ZinvisibleFitObject::ZinvisibleFitObject (const ZinvisibleFitObject& rhs)
  : cachevalid(false), ctheta(0), stheta(0), cphi(0), sphi(0),p2(0), p(0), dpdE(0), pt(0), px(0), py(0), pz(0), dptdE(0),
    dpxdE(0), dpydE(0), dpxdtheta(0), dpydtheta(0), chi2(0)
{
  //std::cout << "copying ZinvisibleFitObject with name" << rhs.name << std::endl;
  ZinvisibleFitObject::assign (rhs);
}

ZinvisibleFitObject& ZinvisibleFitObject::operator= (const ZinvisibleFitObject& rhs) {
  if (this != &rhs) {
    assign (rhs); // calls virtual function assign of derived class
  }
  return *this;
}

ZinvisibleFitObject *ZinvisibleFitObject::copy() const {
  return new ZinvisibleFitObject (*this);
}
    
ZinvisibleFitObject& ZinvisibleFitObject::assign (const BaseFitObject& source) {
  if (const ZinvisibleFitObject *psource = dynamic_cast<const ZinvisibleFitObject *>(&source)) {
    if (psource != this){
      ParticleFitObject::assign (source);
      // only mutable data members, need not to be copied, if cache is invalid
    }
  }
  else {
    assert (0);
  }
  return *this;
}

const char *ZinvisibleFitObject::getParamName (int ilocal) const {
  switch (ilocal) {
    case 0: return "E";
    case 1: return "theta";
    case 2: return "phi";
  }
  return "undefined";
}

bool ZinvisibleFitObject::updateParams (double pp[], int idim) {

  invalidateCache();
  int iE  = getGlobalParNum(0);
  int ith = getGlobalParNum(1);
  int iph = getGlobalParNum(2);
  assert (iE  >= 0 && iE  < idim);
  assert (ith >= 0 && ith < idim);
  assert (iph >= 0 && iph < idim);
  
  double e  = pp[iE];
  double th = pp[ith];
  double ph = pp[iph];
  if (e<0) {
    // cout << "ZinvisibleFitObject::updateParams: mirrored E!\n";
    e  = -e;
    th = M_PI-th;
    ph = M_PI+ph;
  }

  double massPlusEpsilon = mass*(1.0000001);
  if (e < massPlusEpsilon) e = massPlusEpsilon;

  if (th<0 || th>M_PI) {
    // cout << "ZinvisibleFitObject::updateParams: mirrored theta!\n";
    th = M_PI-th;
    ph = M_PI+ph;
  }
  
  bool result = (e -par[0])*(e -par[0]) > eps2*cov[0][0] ||
                (th-par[1])*(th-par[1]) > eps2*cov[1][1] ||
                (ph-par[2])*(ph-par[2]) > eps2*cov[2][2];

  par[0] = (e >= mass) ? e:mass ;
  par[1] = (th >= 0 && th < M_PI) ? 
            th : std::acos (std::cos (th));
  if (std::abs(ph) > M_PI) ph = atan2 (sin(ph), cos (ph));          
  par[2] = ph; 
  pp[iE]  = par[0];         
  pp[ith] = par[1];         
  pp[iph] = par[2];         
  return result;
}  

// these depend on actual parametrisation!
double ZinvisibleFitObject::getPx() const {
  if (!cachevalid) updateCache();
  return px;
}
double ZinvisibleFitObject::getPy() const {
  if (!cachevalid) updateCache();
  return py;
}
double ZinvisibleFitObject::getPz() const {
  if (!cachevalid) updateCache();
  return pz;
}
double ZinvisibleFitObject::getE() const {
  return par[0];
}

double ZinvisibleFitObject::getP() const {
    if (!cachevalid) updateCache();
    return p; 
}

double ZinvisibleFitObject::getP2() const {
  if (!cachevalid) updateCache();
   return p2; 
}
double ZinvisibleFitObject::getPt() const {
  if (!cachevalid) updateCache();
  return pt;
}
double ZinvisibleFitObject::getPt2() const {
  if (!cachevalid) updateCache();
  return pt*pt;
}

double ZinvisibleFitObject::getDPx(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpxdE;
    case 1: return dpxdtheta;
    case 2: return -py;
  }
  return 0; 
}

double ZinvisibleFitObject::getDPy(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpydE;
    case 1: return dpydtheta;
    case 2: return px;
  }
  return 0; 
}

double ZinvisibleFitObject::getDPz(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return ctheta;
    case 1: return -pt;
    case 2: return 0;
  }
  return 0; 
}

double ZinvisibleFitObject::getDE(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  switch (ilocal) {
    case 0: return 1;
    case 1: return 0;
    case 2: return 0;
  }
  return 0; 
}

double ZinvisibleFitObject::getFirstDerivative_Meta_Local( int iMeta, int ilocal , int metaSet ) const {
  // iMeta = intermediate variable (i.e. E,px,py,pz)
  // ilocal = local variable (E, theta, phi)
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

double ZinvisibleFitObject::getSecondDerivative_Meta_Local( int iMeta, int ilocal , int jlocal , int metaSet ) const {
  assert ( metaSet==0 );
  if (!cachevalid) updateCache();

  if ( jlocal<ilocal ) {
    int temp=jlocal;
    jlocal=ilocal;
    ilocal=temp;
  }

  // daniel hasn't checked these, copied from orig code
  switch ( iMeta ) {

  case 0:
    return 0;
    break;
  case 1:
    if      ( ilocal==0 && jlocal==1 ) return ctheta*cphi;
    else if ( ilocal==0 && jlocal==2 ) return -dpydE;
    else if ( ilocal==1 && jlocal==1 ) return -px;
    else if ( ilocal==1 && jlocal==2 ) return -dpydtheta;
    else if ( ilocal==2 && jlocal==2 ) return -px;
    else return 0;
    break;
  case 2:
    if      ( ilocal==0 && jlocal==1 ) return ctheta*sphi;
    else if ( ilocal==0 && jlocal==2 ) return dpxdE;
    else if ( ilocal==1 && jlocal==1 ) return -py;
    else if ( ilocal==1 && jlocal==2 ) return dpxdtheta;
    else if ( ilocal==2 && jlocal==2 ) return -py;
    else return 0;
    break;
  case 3:
    if      ( ilocal==0 && jlocal==1 ) return -stheta;
    else if ( ilocal==1 && jlocal==1 ) return -pz;
    else return 0;
    break;
  default:
    assert(0);
  }
  return -999;
}

void ZinvisibleFitObject::invalidateCache() const {
  cachevalid = false;
}
    
void ZinvisibleFitObject::updateCache() const {
  double e     = par[0];
  double theta = par[1];
  double phi   = par[2];

  ctheta = cos(theta);
  stheta = sin(theta);
  cphi   = cos(phi);
  sphi   = sin(phi);

  p2 = std::abs(e*e-mass*mass);
  p= std::sqrt(p2);
  dpdE = e/p;
 
  pt = p*stheta;

  px = pt*cphi;
  py = pt*sphi;
  pz = p*ctheta; 

  dpxdE = stheta*cphi;
  dpydE = stheta*sphi;
  dpxdtheta = pz*cphi;
  dpydtheta = pz*sphi;
 
  cachevalid = true;
}
