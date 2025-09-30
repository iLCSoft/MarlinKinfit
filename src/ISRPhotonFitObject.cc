/*! \file 
 *  \brief Implements class ISRPhotonFitObject
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: ISRPhotonFitObject.cc,v $
 * - Revision 1.4  2011/03/16 16:33:24  mbeckman
 * - Compatibility fixes with ILCSoft svn
 * -
 * - Revision 1.3  2011/03/03 15:32:32  boehmej
 * - removed obsolete PhotonFitObject.o and PhotonFitObjectPgxy.o from src/Makefile, activated NO_MARLIN in src/ISRPhotonFitObject.cc
 * -
 * - Revision 1.2  2010/07/05 20:08:43  mbeckman
 * - ISRPhotonFitObject.cc: Included flag for output via cout/marlin
 * -
 * - Revision 1.1  2010/06/11 20:32:51  mbeckman
 * - Renamed PhotonFitObjects, cleaned them up for usage
 * -
 * - Revision 1.10  2009/04/02 12:47:35  mbeckman
 * - PhotonFitObject.cc, PseudoMeasuredPhotonFitObjectPxyz.cc: bug fix (measured p = 0 instead of start value)
 * - PhotonFitObjectPxyg.cc: added assertion to catch up division by zero
 * -
 * - Revision 1.9  2009/04/01 09:00:15  mbeckman
 * - Corrected derivatives by factor sqrt(2)
 * -
 * - Revision 1.8  2009/03/26 08:47:30  mbeckman
 * - Bug fix (measured p = 0 instead of start value)
 * -
 * - Revision 1.7  2009/03/17 13:27:03  mbeckman
 * - fixed for compiling with gcc 3.2.3 (for compability with ILCSoft)
 * -
 * - Revision 1.6  2009/02/23 12:04:05  mbeckman
 * - - PhotonFitObject:     bug fix (1/0), removed dispensable variables
 * - - PhotonFitObjectPxyg: bug fixes (1/0, order of computing variables), modified parametrization
 * - - JetFitObject:        added start parameter check (inf, nan)
 * -
 * - Revision 1.5  2009/02/18 11:56:22  mbeckman
 * - PhotonFitObject*.cc: documentation, debug output
 * - NewtonFitterGSL.cc:  bug fix (Lagrange multipliers not initialized), debug output
 * - JetFitObject.cc:     bug fix: division by 0, if energy <= mass
 *
 */ 

#define NO_MARLIN		// if defined: all output via cout, Marlin inclusion not required
#include "ISRPhotonFitObject.h"
#include <cmath>

#undef NDEBUG
#include <cassert>

#include <iostream>
#ifndef NO_MARLIN
#include "marlin/Processor.h"
#endif

using std::sqrt;
using std::exp;
using std::pow;
using std::cout; 
using std::endl;
#ifndef NO_MARLIN
using namespace marlin;
#endif

static const double pi_ = M_PI, // 3.14159265358979323846264338328,
                    a   = 8./3./pi_*(pi_-3.)/(4.-pi_);    // = ca. 0.140012289

// constructor
ISRPhotonFitObject::ISRPhotonFitObject(double px, double py, double ppz,
                                         double b_, double PzMaxB_, double PzMinB_) 
  : cachevalid(false),    
    pt2(0), p2(0), p(0), pz(0),
    dpx0(0), dpy0(0), dpz0(0), dE0(0), dpx1(0), dpy1(0), dpz1(0), dE1(0),
    dpx2(0), dpy2(0), dpz2(0), dE2(0), d2pz22(0), d2E22(0),
    chi2(0), b(0), PzMinB(0), PzMaxB(0), dp2zFact(0)
{

  assert( int(NPAR) <= int(BaseDefs::MAXPAR) );

  initCov();
  b = b_;
  PzMinB = PzMinB_;
  PzMaxB = PzMaxB_;
  #ifdef DEBUG
    cout << "ISRPhotonFitObject:   b: " << b << "   PzMinB: " << PzMinB << "   PzMaxB: " << PzMaxB << endl;
  #endif

  if(b <= 0. || b >= 1.){
    cout << "ISRPhotonFitObject:   b must be from ]0,1[ "  << endl;
  }
  assert(b > 0. && b < 1.);
  if(PzMinB < 0. || PzMaxB <= PzMinB){
    cout << "ISRPhotonFitObject:   PzMinB and PzMaxB must be chosen such that 0 <= PzMinB < PzMaxB"  << endl;
  }
  assert(PzMinB >= 0.);
  assert(PzMaxB > PzMinB);
  dp2zFact = (PzMaxB-PzMinB)/b*sqrt(2./pi_);
  double pg = PgFromPz(ppz);         // using internally Gauss-distributed parameter p_g instead of p_z
  setParam (0, px, true, true);
  setParam (1, py, true, true);
  setParam (2, pg, true);
  setMParam (0, 0.);                // all measured parameters
  setMParam (1, 0.);                // are assumed to be zero
  setMParam (2, 0.);                // in this photon parametrization
  #ifdef DEBUG
    cout << "ISRPhotonFitObject:   Initial pg: " << pg << endl;
  #endif
  setError (2, 1.);
  setMass (0.);
  invalidateCache();
}

// destructor
ISRPhotonFitObject::~ISRPhotonFitObject() {}

// We get a warning that ParticleFitObject should be explicitly initialized
// here, but I don't want to change this part because, I think everything is
// done properly already and not changing behavior is more important.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra"
ISRPhotonFitObject::ISRPhotonFitObject (const ISRPhotonFitObject& rhs)
  : cachevalid(false),    
    pt2(0), p2(0), p(0), pz(0),
    dpx0(0), dpy0(0), dpz0(0), dE0(0), dpx1(0), dpy1(0), dpz1(0), dE1(0),
    dpx2(0), dpy2(0), dpz2(0), dE2(0), d2pz22(0), d2E22(0),
    chi2(0), b(0), PzMinB(0), PzMaxB(0), dp2zFact(0)
{
  //std::cout << "copying ISRPhotonFitObject with name" << rhs.name << std::endl;
  ISRPhotonFitObject::assign (rhs);
}
#pragma GCC diagnostic pop

ISRPhotonFitObject& ISRPhotonFitObject::operator= (const ISRPhotonFitObject& rhs) {
  if (this != &rhs) {
    assign (rhs); // calls virtual function assign of derived class
  }
  return *this;
}

ISRPhotonFitObject *ISRPhotonFitObject::copy() const {
  return new ISRPhotonFitObject (*this);
}
    
ISRPhotonFitObject& ISRPhotonFitObject::assign (const BaseFitObject& source) {
  if (const ISRPhotonFitObject *psource = dynamic_cast<const ISRPhotonFitObject *>(&source)) {
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

const char *ISRPhotonFitObject::getParamName (int ilocal) const {
  switch (ilocal) {
    case 0: return "P_x";
    case 1: return "P_y";
    case 2: return "P_g";
  }
  return "undefined";
}
 
bool ISRPhotonFitObject::updateParams (double pp[], int idim) {
  invalidateCache();
  int i2 = getGlobalParNum(2);
  assert (i2 >= 0 && i2 < idim);
  double pp2 = pp[i2];
  #ifdef DEBUG
    std::cout << "ISRPhotonFitObject::updateParams:   p2(new) = " << pp[i2] << "   par[2](old) = " << par[2] << endl;
  #endif
  bool result = ((pp2-par[2])*(pp2-par[2]) > eps2*cov[2][2]);
  par[2] = pp2;
  pp[i2] = par[2];
  return result;
}  

double ISRPhotonFitObject::getDPx(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpx0;
    case 1: return dpx1;
    case 2: return dpx2;
  }
  return 0; 
}

double ISRPhotonFitObject::getDPy(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpy0;
    case 1: return dpy1;
    case 2: return dpy2;
  }
  return 0; 
}

double ISRPhotonFitObject::getDPz(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpz0;
    case 1: return dpz1;
    case 2: return dpz2;
  }
  return 0; 
}

double ISRPhotonFitObject::getDE(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dE0;
    case 1: return dE1;
    case 2: return dE2;
  }
  return 0; 
}
 

double ISRPhotonFitObject::getFirstDerivative_Meta_Local( int iMeta, int ilocal , int metaSet ) const {
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

double ISRPhotonFitObject::getSecondDerivative_Meta_Local( int iMeta, int ilocal , int jlocal, int metaSet ) const {
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
    if      ( ilocal==2 && jlocal==2 ) return  d2E22;
    else return 0;
    break;
  case 1:
  case 2:
    return 0;
    break;
  case 3:
    if      ( ilocal==2 && jlocal==2 ) return  d2pz22;
    else return 0;
    break;
  default:
    return 0;
  }

}


         
double ISRPhotonFitObject::PgFromPz(double ppz){

  int sign = (ppz>0.) - (ppz<0.);
  double u = ( pow(fabs(ppz),b) - PzMinB ) / (PzMaxB-PzMinB);

   if(u<0.){
   #ifdef NO_MARLIN
     cout << 
   #else
     m_out(WARNING) << 
   #endif
     "ISRPhotonFitObject: Initial pz with abs(pz) < pzMin adjusted to zero." << std::endl;
     u = 0.;
   }
 
   if(u>=1.){
//    #ifdef NO_MARLIN
//      cout << 
//    #else
//      m_out(WARNING) << 
//    #endif
//      "ISRPhotonFitObject: Initial pz with abs(pz) >= pzMax adjusted." << std::endl;
     u = 0.99999999;
   }

  double g = std::log(1.-u*u);
  double g4pa = g + 4./pi_/a;
  return sign*sqrt( -g4pa+sqrt( g4pa*g4pa-4./a*g ) ) ;
}


void ISRPhotonFitObject::updateCache() const {
  double px = par[0];
  double py = par[1];
  double pg = par[2];

  int sign = (pg>0.) - (pg<0.);
  double pg2h = pg*pg/2.;
  double exponent = -pg2h*(4./pi_+a*pg2h)/(1.+a*pg2h);
  double u = sqrt( (exponent<-1.e-14) ? 1.-exp( exponent ) : -exponent );  // approximation to avoid numerical problem
  pz = sign*pow( ( PzMinB + (PzMaxB-PzMinB)*u ) , (1./b) );

  pt2 = px*px+py*py;
  p2  = pt2+pz*pz;
  p   = std::sqrt(p2);

  fourMomentum.setValues( p, px, py, pz );
  
  dpx0 = 1.;
  dpx1 = 0.;
  dpx2 = 0.;
  dpy0 = 0.;
  dpy1 = 1.;
  dpy2 = 0.;
  dpz0 = 0.;
  dpz1 = 0.;
  dpz2 = dp2zFact*pow(fabs(pz),(1.-b))*exp(-pg*pg/2.);
  
  // if p,pz==0, derivatives are zero (catch up 1/0)
  if(pz){
    d2pz22 = dpz2*( (1.-b)*dpz2/pz - par[2] );
  }
  else{
    d2pz22 = 0.;
  }
  if(p){
    dE0   = px/p;
    dE1   = py/p;
    dE2   = pz/p*dpz2;
    d2E22 = pz/p*d2pz22;
    if(pt2){
      d2E22 += pt2/p/p/p*dpz2*dpz2;   // NOT using /p2/p to avoid numerical problems
    }
  }
  else{
    dE0   = 0.;
    dE1   = 0.;
    dE2   = 0.;
    d2E22 = 0.;
  }

  #ifdef DEBUG
    cout << "ISRPhotonFitObject::updateCache:   pg: " << pg << "   pz: " << pz << "   p: " << p << "   p^2: " << p2 << "\n"
         << "                                    Dpz/Dpg: " << dpz2 << "   DE/Dpg: " << dE2 << "   D^2pz/Dpg^2: " << d2pz22
                                          << "   D^2E/Dpg^2: " << d2E22 << endl;
  #endif
  
  cachevalid = true;
}
