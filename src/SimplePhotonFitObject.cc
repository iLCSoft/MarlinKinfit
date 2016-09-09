/*! \file 
 *  \brief Implements class SimplePhotonFitObject
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: SimplePhotonFitObject.cc,v $
 * - Revision 1.1  2010/06/11 20:32:51  mbeckman
 * - Renamed PhotonFitObjects, cleaned them up for usage
 * -
 * - Revision 1.8  2010/05/25 13:23:56  boehmej
 * - fixed RootTracer
 * -
 * - Revision 1.7  2009/04/02 12:47:35  mbeckman
 * - PhotonFitObject.cc, PseudoMeasuredPhotonFitObjectPxyz.cc: bug fix (measured p = 0 instead of start value)
 * - PhotonFitObjectPxyg.cc: added assertion to catch up division by zero
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
 * -
 *
 */ 

#include "SimplePhotonFitObject.h"
#include <cmath>
#undef NDEBUG
#include <cassert>
#include <iostream>
//#include "marlin/Processor.h"

using std::sqrt;
using std::cout; 
using std::endl;

// constructor
SimplePhotonFitObject::SimplePhotonFitObject(double px, double py, double pz, double Dpz) : pt2(0), p2(0), p(0),dE0(0), dE1(0), dE2(0),chi2(0)
{

  assert( int(NPAR) <= int(BaseDefs::MAXPAR) );

  initCov();                         
  setParam (0, px, true, true);
  setParam (1, py, true, true);
  setParam (2, pz, true);
  setMParam (0, 0.);
  setMParam (1, 0.);
  setMParam (2, 0.);
  setError (2, Dpz);
  setMass (0.);
  invalidateCache();
}

// destructor
SimplePhotonFitObject::~SimplePhotonFitObject() {}

SimplePhotonFitObject::SimplePhotonFitObject (const SimplePhotonFitObject& rhs) : pt2(0), p2(0), p(0),dE0(0), dE1(0), dE2(0),chi2(0)
{
  //std::cout << "copying SimplePhotonFitObject with name" << rhs.name << std::endl;
  SimplePhotonFitObject::assign (rhs);
}

SimplePhotonFitObject& SimplePhotonFitObject::operator= (const SimplePhotonFitObject& rhs) {
  if (this != &rhs) {
    assign (rhs); // calls virtual function assign of derived class
  }
  return *this;
}

SimplePhotonFitObject *SimplePhotonFitObject::copy() const {
  return new SimplePhotonFitObject (*this);
}
    
SimplePhotonFitObject& SimplePhotonFitObject::assign (const BaseFitObject& source) {
  if (const SimplePhotonFitObject *psource = dynamic_cast<const SimplePhotonFitObject *>(&source)) {
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

const char *SimplePhotonFitObject::getParamName (int ilocal) const {
  switch (ilocal) {
    case 0: return "P_x";
    case 1: return "P_y";
    case 2: return "P_z";
  }
  return "undefined";
}
 
bool SimplePhotonFitObject::updateParams (double pp[], int idim) {
  invalidateCache();
  
  int i2 = getGlobalParNum(2);
// std::cout << "updateParams: i2 = " << i2 << "\n";
// std::cout << "updateParams: idim = " << idim << "\n";
  assert (i2 >= 0 && i2 < idim);
  
  double pp2 = pp[i2];
// std::cout << "updateParams: p2 = " << p[i2] << "   par[2] = " << par[2] << "\n";
  
  bool result = ((pp2-par[2])*(pp2-par[2]) > eps2*cov[2][2]);

  par[2] = pp2;
  pp[i2] = par[2];         
  return result;
}  

// these depend on actual parametrisation!
double SimplePhotonFitObject::getDPx(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
//   if (!cachevalid) updateCache();
  switch (ilocal) {
//     case 0: return dpx0;
//     case 1: return dpx1;
//     case 2: return dpx2;
    case 0: return 1.;
    case 1: return 0.;
    case 2: return 0.;
  }
  return 0; 
}

double SimplePhotonFitObject::getDPy(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
//   if (!cachevalid) updateCache();
  switch (ilocal) {
//     case 0: return dpy0;
//     case 1: return dpy1;
//     case 2: return dpy2;
    case 0: return 0.;
    case 1: return 1.;
    case 2: return 0.;
  }
  return 0; 
}

double SimplePhotonFitObject::getDPz(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
//   if (!cachevalid) updateCache();
  switch (ilocal) {
//     case 0: return dpz0;
//     case 1: return dpz1;
//     case 2: return dpz2;
    case 0: return 0.;
    case 1: return 0.;
    case 2: return 1.;
  }
  return 0; 
}

double SimplePhotonFitObject::getDE(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dE0;
    case 1: return dE1;
    case 2: return dE2;
  }
  return 0; 
}

double SimplePhotonFitObject::getFirstDerivative_Meta_Local( int iMeta, int ilocal , int metaSet ) const {

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
    assert(0); // should never get here
  }
  // should really never get here! just to get rid of compiler warning
  return -999;
}


double SimplePhotonFitObject::getSecondDerivative_Meta_Local( int iMeta, int ilocal , int jlocal , int metaSet ) const {
  assert ( metaSet==0 );
  if (!cachevalid) updateCache();

  if ( jlocal<ilocal ) {
    int temp=jlocal;
    jlocal=ilocal;
    ilocal=temp;
  }
  
  // daniel hasn't checked these, copied from orig code
  // DANIEL needs to check this. i think it's probably nonsense...
  switch ( iMeta ) {

  case 0:
    if      ( ilocal==0 && jlocal==0 ) return pt2/p/p/p;
    else  return 0;
    break;
  case 1: // fall through
  case 2: // fall through
  case 3:
  default:
    return 0;
  }

}

         

void SimplePhotonFitObject::updateCache() const {
  // std::cout << "SimplePhotonFitObject::updateCache" << std::endl;
  double px = par[0];
  double py = par[1];
  double pz = par[2];

  pt2 = px*px+py*py;
  p2  = pt2+pz*pz;
  p   = std::sqrt(p2);

  fourMomentum.setValues(p, px, py, pz);
  
//   dpx0 = 1.;
//   dpx1 = 0.;
//   dpx2 = 0.;
//   dpy0 = 0.;
//   dpy1 = 1.;
//   dpy2 = 0.;
//   dpz0 = 0.;
//   dpz1 = 0.;
//   dpz2 = 1.;
  // if p==0, derivates are zero (catch up division by zero)
  if(p){
    dE0  = px/p;
    dE1  = py/p;
    dE2  = pz/p;
  }
  else{
    dE0  = 0.;
    dE1  = 0.;
    dE2  = 0.;
  }

  cachevalid = true;
}
