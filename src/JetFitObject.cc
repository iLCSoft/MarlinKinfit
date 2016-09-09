/*! \file 
 *  \brief Implements class JetFitObject
 *
 * \b Changelog:
 * - 26.09.2008 mbeckman: Minor bug fixes, dodged possible division by zero
 *
 * \b CVS Log messages:
 * - $Log: JetFitObject.cc,v $
 * - Revision 1.9  2011/03/03 15:03:03  blist
 * - Latest version, with NewFitterGSL
 * -
 * - Revision 1.8  2010/07/05 19:52:33  mbeckman
 * - JetFitObject.cc: commented out assert(!isnan(...))
 * -
 * - Revision 1.7  2009/02/23 12:04:05  mbeckman
 * - - PhotonFitObject:     bug fix (1/0), removed dispensable variables
 * - - PhotonFitObjectPxyg: bug fixes (1/0, order of computing variables), modified parametrization
 * - - JetFitObject:        added start parameter check (inf, nan)
 * -
 * - Revision 1.6  2009/02/18 11:56:21  mbeckman
 * - PhotonFitObject*.cc: documentation, debug output
 * - NewtonFitterGSL.cc:  bug fix (Lagrange multipliers not initialized), debug output
 * - JetFitObject.cc:     bug fix: division by 0, if energy <= mass
 * -
 * - Revision 1.5  2009/02/17 12:46:34  blist
 * - Improved version of NewtonFitterGSL, JetFitObject changed
 * -
 * - Revision 1.4  2008/10/16 08:13:44  blist
 * - New versions of OPALfitter and Newtonfitter using GSL
 * -
 * - Revision 1.3  2008/09/26 15:26:39  mbeckman
 * - Added fit objects for pseudo-measured photons, minor changes/bug fixes in fitters and JetFitObject
 * -
 * - Revision 1.2  2008/09/26 09:58:10  boehmej
 * - removed ~100 semicolons after } at end of function implementation :)
 * -
 * - Revision 1.1  2008/02/12 10:19:08  blist
 * - First version of MarlinKinfit
 * -
 * - Revision 1.10  2008/02/07 08:15:25  blist
 * - error calculation of constraints fixed
 * -
 * - Revision 1.9  2008/02/04 17:30:53  blist
 * - NewtonFitter works now!
 * -
 * - Revision 1.8  2008/01/30 21:48:02  blist
 * - Newton Fitter still doesnt work :-(
 * -
 * - Revision 1.7  2008/01/30 09:14:53  blist
 * - Preparations for NewtonFitter
 * -
 * - Revision 1.6  2008/01/29 17:22:59  blist
 * - new addTo2ndDerivatives and setParam
 * -
 * - Revision 1.5  2007/09/17 12:50:15  blist
 * - Some parameters reordered
 * -
 * - Revision 1.4  2007/09/13 08:09:50  blist
 * - Updated 2nd derivatives for px,py,pz,E constraints, improved header documentation
 * -
 *
 */ 

#include "JetFitObject.h"
#include <cmath>

#undef NDEBUG
#include <cassert>

#include <iostream>

using std::sqrt;
using std::sin;
using std::cos;
using std::cout; 
using std::endl;

// constructor
JetFitObject::JetFitObject(double E, double theta, double phi,  
                           double DE, double Dtheta, double Dphi, 
                           double m)
  : ctheta(0), stheta(0), cphi(0), sphi(0),
    p2(0), p(0), pt(0), px(0), py(0), pz(0), dpdE(0), dptdE(0), 
    dpxdE(0), dpydE(0), dpzdE(0), dpxdtheta(0), dpydtheta(0), chi2(0)
{

  assert( int(NPAR) <= int(BaseDefs::MAXPAR) );

  initCov();                         
//  assert( !isinf(E) );        assert( !isnan(E) );
//  assert( !isinf(theta) );    assert( !isnan(theta) );
//  assert( !isinf(phi) );      assert( !isnan(phi) );
//  assert( !isinf(DE) );       assert( !isnan(DE) );
//  assert( !isinf(Dtheta) );   assert( !isnan(Dtheta) );
//  assert( !isinf(Dphi) );     assert( !isnan(Dphi) );
//  assert( !isinf(m) );        assert( !isnan(m) );
  setMass (m);
  adjustEThetaPhi (m, E, theta, phi);
  setParam (0, E, true);
  setParam (1, theta, true);
  setParam (2, phi, true);
  setMParam (0, E);
  setMParam (1, theta);
  setMParam (2, phi);
  setError (0, DE);
  setError (1, Dtheta);
  setError (2, Dphi);

  // parameter 2 repeats every 2*pi
  paramCycl[2]=2.*M_PI;

  invalidateCache();
//   std::cout << "JetFitObject::JetFitObject: E = " << E << std::endl;
//   std::cout << "JetFitObject::JetFitObject: getParam(0) = " << getParam(0) << std::endl;
//   std::cout << "JetFitObject::JetFitObject: " << *this << std::endl;
//   std::cout << "mpar= " << mpar[0] << ", " << mpar[1] << ", " << mpar[2] << std::endl;
}

// destructor
JetFitObject::~JetFitObject() {}

JetFitObject::JetFitObject (const JetFitObject& rhs)
  : ctheta(0), stheta(0), cphi(0), sphi(0),
    p2(0), p(0), pt(0), px(0), py(0), pz(0), dpdE(0), dptdE(0), 
    dpxdE(0), dpydE(0), dpzdE(0), dpxdtheta(0), dpydtheta(0), chi2(0)
{
  //std::cout << "copying JetFitObject with name " << rhs.name << std::endl;
  JetFitObject::assign (rhs);
}

JetFitObject& JetFitObject::operator= (const JetFitObject& rhs) {
  if (this != &rhs) {
    assign (rhs); // calls virtual function assign of derived class
  }
  return *this;
}

JetFitObject *JetFitObject::copy() const {
  return new JetFitObject (*this);
}
    
JetFitObject& JetFitObject::assign (const BaseFitObject& source) {
  if (const JetFitObject *psource = dynamic_cast<const JetFitObject *>(&source)) {
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

const char *JetFitObject::getParamName (int ilocal) const {
  switch (ilocal) {
    case 0: return "E";
    case 1: return "theta";
    case 2: return "phi";
  }
  return "undefined";
}

 
bool JetFitObject::updateParams (double pp[], int idim) {
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
    // cout << "JetFitObject::updateParams: mirrored E!\n";
    e  = -e;
    th = M_PI-th;
    ph = M_PI+ph;
  }
  
  double massPlusEpsilon = mass*(1.0000001);
  if (e < massPlusEpsilon) e = massPlusEpsilon;
  
  bool result = ((e -par[0])*(e -par[0]) > eps2*cov[0][0]) ||
                ((th-par[1])*(th-par[1]) > eps2*cov[1][1]) ||
                ((ph-par[2])*(ph-par[2]) > eps2*cov[2][2]);
                
  par[0] = e;
  par[1] = th;
  par[2] = ph;
  pp[iE]  = par[0];         
  pp[ith] = par[1];         
  pp[iph] = par[2];         
  return result;
}  

double JetFitObject::getDPx(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpxdE;
    case 1: return dpxdtheta;
    case 2: return -py;
  }
  return 0; 
}

double JetFitObject::getDPy(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpydE;
    case 1: return dpydtheta;
    case 2: return px;
  }
  return 0; 
}

double JetFitObject::getDPz(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpzdE;
    case 1: return -pt;
    case 2: return 0;
  }
  return 0; 
}

double JetFitObject::getDE(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  switch (ilocal) {
    case 0: return 1;
    case 1: return 0;
    case 2: return 0;
  }
  return 0; 
}

double JetFitObject::getError (int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  //cout << "JetFitObject::getError (ilocal = " << ilocal << ") = " <<  std::sqrt(cov[ilocal][ilocal]) << endl;
  return std::sqrt(cov[ilocal][ilocal]);
}

double JetFitObject::getCov (int ilocal, int jlocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  assert (jlocal >= 0 && jlocal < NPAR);
  return cov[ilocal][jlocal];
} 

//void   JetFitObject::addToDerivatives (double der[], int idim, 
//                                       double efact, double pxfact, 
//                                       double pyfact, double pzfact) const {
//  int i_E     = globalParNum[0];
//  int i_theta = globalParNum[1];
//  int i_phi   = globalParNum[2];
//  assert (i_E     >= 0 && i_E     < idim);
//  assert (i_theta >= 0 && i_theta < idim);
//  assert (i_phi   >= 0 && i_phi   < idim);
//  
//  if (!cachevalid) updateCache();
//  // for numerical accuracy, add up derivatives first,
//  // then add them to global vector
//  double der_E = efact;
//  double der_theta = 0;
//  double der_phi = 0;
//  
//  if (pxfact != 0) {
//    der_E     += pxfact*dpxdE;
//    der_theta += pxfact*dpxdtheta;
//    der_phi   -= pxfact*py;



double JetFitObject::getFirstDerivative_Meta_Local( int iMeta, int ilocal , int metaSet ) const {

  assert ( metaSet==0 );
  switch ( iMeta ) {
  case 0:
    return getDE(ilocal);
    break;
  case 1:
    return getDPx(ilocal);
    break;
  case 2:
    return getDPy(ilocal);
    break;
  case 3:
    return getDPz(ilocal);
    break;
  default:
    assert(0);

  }
  return -999;
}


double JetFitObject::getSecondDerivative_Meta_Local( int iMeta, int ilocal , int jlocal , int metaSet ) const {
  assert ( metaSet==0 );
  if (!cachevalid) updateCache();

  if ( jlocal<ilocal ) {
    int temp=jlocal;
    jlocal=ilocal;
    ilocal=temp;
  }

  double d2pdE2 = (mass != 0) ? -mass*mass/(p*p*p) : 0;
  double d2ptdE2 = d2pdE2*stheta;
  
  // daniel hasn't checked these, copied from orig code
  switch ( iMeta ) {

  case 0:
    return 0;
    break;
  case 1:
    if      ( ilocal==0 && jlocal==0 ) return  d2ptdE2*cphi;
    else if ( ilocal==0 && jlocal==1 ) return  dpzdE*cphi;
    else if ( ilocal==0 && jlocal==2 ) return -dpydE;
    else if ( ilocal==1 && jlocal==1 ) return -px;
    else if ( ilocal==1 && jlocal==2 ) return -dpydtheta;
    else if ( ilocal==2 && jlocal==2 ) return -px;
    break;
  case 2:
    if      ( ilocal==0 && jlocal==0 ) return  d2ptdE2*sphi;
    else if ( ilocal==0 && jlocal==1 ) return  dpzdE*sphi;
    else if ( ilocal==0 && jlocal==2 ) return  dpxdE;
    else if ( ilocal==1 && jlocal==1 ) return -py;
    else if ( ilocal==1 && jlocal==2 ) return  dpxdtheta;
    else if ( ilocal==2 && jlocal==2 ) return -py;
    break;
  case 3:
    if      ( ilocal==0 && jlocal==0 ) return d2pdE2*ctheta;
    //    else if ( ilocal==0 && jlocal==1 ) return dptdE;
    else if ( ilocal==0 && jlocal==1 ) return -dptdE; // this is "-" in the orig JetFitObject, DJ fixed 2015may27
    else if ( ilocal==0 && jlocal==2 ) return 0;
    else if ( ilocal==1 && jlocal==1 ) return -pz;
    else if ( ilocal==1 && jlocal==2 ) return 0;
    else if ( ilocal==2 && jlocal==2 ) return 0;
    break;
  default:
    assert(0);
  }
  return -999;

}


void JetFitObject::updateCache() const {

  double e     = par[0];
  double theta = par[1];
  double phi   = par[2];

  ctheta = cos(theta);
  stheta = sin(theta);
  cphi   = cos(phi);
  sphi   = sin(phi);

  p2 = std::abs(e*e-mass*mass);
  p = std::sqrt(p2);
  assert (p != 0);

  pt = p*stheta;

  px = pt*cphi;
  py = pt*sphi;
  pz = p*ctheta;

  fourMomentum.setValues(e, px, py, pz);
  
  dpdE = e/p;
  dptdE = dpdE*stheta;
  dpxdE = dptdE*cphi;
  dpydE = dptdE*sphi;
  dpzdE = dpdE*ctheta;
  
  dpxdtheta = pz*cphi;
  dpydtheta = pz*sphi;
  
  cachevalid = true;

}

//double JetFitObject::getChi2 () const {
//  if (!cachevalid) updateCache();
//  return chi2;
//}
//

bool JetFitObject::adjustEThetaPhi (double& m, double &E, double& theta, double& phi) {
  bool result = false;
  
  if (E<0) {
    // cout << "JetFitObject::adjustEThetaPhi: mirrored E!\n";
    E  = -E;
    theta = M_PI-theta;
    phi = M_PI+phi;
    result = true;
  }
  if (E < m) {
    E = m;
    result = true;
  }
  if (theta < -M_PI || theta > M_PI) {
    while (theta < -M_PI) theta += 2*M_PI;
    while (theta >  M_PI) theta -= 2*M_PI;
    result = true;
  }
  
  if (theta<0) {
    // cout << "JetFitObject::adjustEThetaPhi: mirrored theta!\n";
    theta = -theta;
    phi = phi > 0 ? phi-M_PI : phi+M_PI;
    result = true;
  }
  else if (theta>M_PI) {
    // cout << "JetFitObject::adjustEThetaPhi: mirrored theta!\n";
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

