/*! \file
 *  \brief Implements class TrackParticleFitObject
 *  TrackParticleFitObject takes an LCIO track or trackstate (parameters & covariance) as input.
 *  it knows about the 4-momentum and "decay plane" of the TrackParticleObject
 *   "decay plane" is the plane defined by the IP and the tangent to the track at PCA
 *  At present it considers the track only at the PCA:
 *    i.e. the momentum/decay plane are considered only at the PCA, not at a general position along the helix
 *    this means it's no good for fitting (significantly) displaced vertices 
 *       in which the tracks have bent significantly between the PCA and vertex positions 
 *    this should be improved in a future update by including extra parameters (as was done for the original Track-based classes)
 */


#include "TrackParticleFitObject.h"
#include <cmath>

#undef NDEBUG
#include <cassert>

#include <iostream>
#include <iomanip>

using std::sqrt;
using std::sin;
using std::cos;
using std::cout;
using std::endl;

const double TrackParticleFitObject::omega_pt_conv = 2.99792458e-4; // for mm, Tesla, GeV
const double TrackParticleFitObject::maxpt = 500; // GeV

// parameter scalings, to get in approx range 1-100 { iD0,iPhi0,iOmega,iZ0,iTanL }
// DJeans is not sure if we really need these: set them to 1 for now (in which case they have no effect)
//const double TrackParticleFitObject::parfact[NPAR] = {1., 1., 1., 1., 1.};
const double TrackParticleFitObject::parfact[NPAR] = {1.e-2, 1., 1.e-3, 1.e-2, 1., 1., 1.};

TrackParticleFitObject::TrackParticleFitObject( const EVENT::Track* trk, double m) 
  : trackReferencePoint( ThreeVector(0,0,0) ),
    trackPlaneNormal( ThreeVector(0,0,0) ),
    trackPcaVector( ThreeVector(0,0,0) ),
    trajectoryPointAtPCA( ThreeVector(0,0,0) ),
    trajectoryPointAtStart( ThreeVector(0,0,0) ),
    trajectoryPointAtEnd( ThreeVector(0,0,0) ),
    momentumAtPCA( ThreeVector(0,0,0) ),
    momentumAtStart( ThreeVector(0,0,0) ),
    momentumAtEnd( ThreeVector(0,0,0) ),
    phi0(0), omega(0), tanl(0), d0(0), z0(0), s_start(0), s_end(0), chi2(0)
{
  invalidateCache();

  double ppar[NPAR];
  ppar[ iD0    ] =  trk->getD0()       ;
  ppar[ iPhi0  ] =  trk->getPhi()      ;
  ppar[ iOmega ] =  trk->getOmega()    ;
  ppar[ iZ0    ] =  trk->getZ0()       ;
  ppar[ iTanL  ] =  trk->getTanLambda();
  ppar[ iStart ] =  0; // this one is not measured
  ppar[ iEnd   ] =  0; // this one is not measured

  double ccov[15];
  for (int i=0; i<15; i++) ccov[i]=trk->getCovMatrix()[i];

  trackReferencePoint.setValues( trk->getReferencePoint()[0],
                                 trk->getReferencePoint()[1],
                                 trk->getReferencePoint()[2] );

  initialise( ppar , ccov, m );
}

TrackParticleFitObject::TrackParticleFitObject( const EVENT::TrackState* trk, double m) 
  : trackReferencePoint( ThreeVector(0,0,0) ),
    trackPlaneNormal( ThreeVector(0,0,0) ),
    trackPcaVector( ThreeVector(0,0,0) ),
    trajectoryPointAtPCA( ThreeVector(0,0,0) ),
    trajectoryPointAtStart( ThreeVector(0,0,0) ),
    trajectoryPointAtEnd( ThreeVector(0,0,0) ),
    momentumAtPCA( ThreeVector(0,0,0) ),
    momentumAtStart( ThreeVector(0,0,0) ),
    momentumAtEnd( ThreeVector(0,0,0) ),
    phi0(0), omega(0), tanl(0), d0(0), z0(0), s_start(0), s_end(0), chi2(0)
{
  invalidateCache();

  double ppar[NPAR];
  ppar[ iD0    ] =  trk->getD0()       ;
  ppar[ iPhi0  ] =  trk->getPhi()      ;
  ppar[ iOmega ] =  trk->getOmega()    ;
  ppar[ iZ0    ] =  trk->getZ0()       ;
  ppar[ iTanL  ] =  trk->getTanLambda();
  ppar[ iStart ] =  0; // this one is not measured
  ppar[ iEnd   ] =  0; // this one is not measured

  double ccov[15];
  for (int i=0; i<15; i++) ccov[i]=trk->getCovMatrix()[i];

  trackReferencePoint.setValues( trk->getReferencePoint()[0],
                                 trk->getReferencePoint()[1],
                                 trk->getReferencePoint()[2] );

  initialise( ppar , ccov, m );
}

TrackParticleFitObject::TrackParticleFitObject( const double* _ppars, const double* _cov, double m, const double* refPt_) 
  : trackReferencePoint( ThreeVector(0,0,0) ),
    trackPlaneNormal( ThreeVector(0,0,0) ),
    trackPcaVector( ThreeVector(0,0,0) ),
    trajectoryPointAtPCA( ThreeVector(0,0,0) ),
    trajectoryPointAtStart( ThreeVector(0,0,0) ),
    trajectoryPointAtEnd( ThreeVector(0,0,0) ),
    momentumAtPCA( ThreeVector(0,0,0) ),
    momentumAtStart( ThreeVector(0,0,0) ),
    momentumAtEnd( ThreeVector(0,0,0) ),
    phi0(0), omega(0), tanl(0), d0(0), z0(0), s_start(0), s_end(0), chi2(0)
{
  invalidateCache();

  assert( int(NPAR) <= int(BaseDefs::MAXPAR) );
  if ( refPt_ ) trackReferencePoint.setValues(refPt_[0],refPt_[1],refPt_[2]);
  else          trackReferencePoint.setValues(0,0,0);

  initialise(_ppars, _cov, m);
}


void TrackParticleFitObject::initialise( const double* _ppars, const double* _cov, double m) {

  //  cout << "hello from  TrackParticleFitObject::initialise" << endl;

  assert( int(NPAR) <= int(BaseDefs::MAXPAR) );

  initCov();
  setMass (m);

  for (int i=0; i<NPAR; i++) {
    bool ffmeasured = i < iStart;
    setParam ( i, _ppars[i]/parfact[i], ffmeasured, false );
    if ( ffmeasured ) 
      setMParam( i, _ppars[i]/parfact[i] );
  }

  setCov( iD0   , iD0    , _cov[ 0] / (parfact[iD0   ]*parfact[iD0   ]) ); // d0  d0
  setCov( iPhi0 , iD0    , _cov[ 1] / (parfact[iPhi0 ]*parfact[iD0   ]) ); // phi d0
  setCov( iPhi0 , iPhi0  , _cov[ 2] / (parfact[iPhi0 ]*parfact[iPhi0 ]) ); // phi phi
  setCov( iOmega, iD0    , _cov[ 3] / (parfact[iOmega]*parfact[iD0   ]) ); // ome d0
  setCov( iOmega, iPhi0  , _cov[ 4] / (parfact[iOmega]*parfact[iPhi0 ]) ); // ome phi
  setCov( iOmega, iOmega , _cov[ 5] / (parfact[iOmega]*parfact[iOmega]) ); // ome ome
  setCov( iZ0   , iD0    , _cov[ 6] / (parfact[iZ0   ]*parfact[iD0   ]) ); // z0  d0
  setCov( iZ0   , iPhi0  , _cov[ 7] / (parfact[iZ0   ]*parfact[iPhi0 ]) ); // z0  phi
  setCov( iZ0   , iOmega , _cov[ 8] / (parfact[iZ0   ]*parfact[iOmega]) ); // z0  ome
  setCov( iZ0   , iZ0    , _cov[ 9] / (parfact[iZ0   ]*parfact[iZ0   ]) ); // z0  z0
  setCov( iTanL , iD0    , _cov[10] / (parfact[iTanL ]*parfact[iD0   ]) ); // tan d0
  setCov( iTanL , iPhi0  , _cov[11] / (parfact[iTanL ]*parfact[iPhi0 ]) ); // tan phi
  setCov( iTanL , iOmega , _cov[12] / (parfact[iTanL ]*parfact[iOmega]) ); // tan ome
  setCov( iTanL , iZ0    , _cov[13] / (parfact[iTanL ]*parfact[iZ0   ]) ); // tan z0
  setCov( iTanL , iTanL  , _cov[14] / (parfact[iTanL ]*parfact[iTanL ]) ); // tan tan

  // parameter iPhi0 repeats every 2*pi
  paramCycl[iPhi0]=2.*M_PI/parfact[iPhi0 ];

  // by default, fix the vertex parameters
  fixVertexParam(0, true);
  fixVertexParam(1, true);

  invalidateCache();

  return;
}


// destructor
TrackParticleFitObject::~TrackParticleFitObject() {}



TrackParticleFitObject::TrackParticleFitObject (const TrackParticleFitObject& rhs)
  : trackReferencePoint( ThreeVector(0,0,0) ),
    trackPlaneNormal( ThreeVector(0,0,0) ),
    trackPcaVector( ThreeVector(0,0,0) ),
    trajectoryPointAtPCA( ThreeVector(0,0,0) ),
    trajectoryPointAtStart( ThreeVector(0,0,0) ),
    trajectoryPointAtEnd( ThreeVector(0,0,0) ),
    momentumAtPCA( ThreeVector(0,0,0) ),
    momentumAtStart( ThreeVector(0,0,0) ),
    momentumAtEnd( ThreeVector(0,0,0) ),
    phi0(0), omega(0), tanl(0), d0(0), z0(0), s_start(0), s_end(0), chi2(0)
{
  //std::cout << "copying TrackParticleFitObject with name " << rhs.name << std::endl;
  TrackParticleFitObject::assign (rhs);
}

TrackParticleFitObject& TrackParticleFitObject::operator= (const TrackParticleFitObject& rhs) {
  if (this != &rhs) {
    assign (rhs); // calls virtual function assign of derived class
  }
  return *this;
}



TrackParticleFitObject *TrackParticleFitObject::copy() const {
  return new TrackParticleFitObject (*this);
}

TrackParticleFitObject& TrackParticleFitObject::assign (const BaseFitObject& source) {
  if (const TrackParticleFitObject *psource = dynamic_cast<const TrackParticleFitObject *>(&source)) {
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

const char *TrackParticleFitObject::getParamName (int ilocal) const {
  switch (ilocal) {
  case iD0   : return "d0" ;
  case iPhi0 : return "phi0" ;
  case iOmega: return "omega" ;
  case iZ0   : return "z0" ;
  case iTanL : return "tanLambda" ;
  case iStart: return "sStart";
  case iEnd  : return "sEnd";
  }
  return "undefined";
}

bool TrackParticleFitObject::updateParams (double p[], int idim) {
  invalidateCache();

  double tempPar[NPAR]={0};

  // check that omega is not too small (pt too large)
  double omegaMin = fabs(omega_pt_conv*getBfield()/maxpt);

  bool result=false;

  for (int i=0; i<getNPar(); i++) {
    int iglobal = getGlobalParNum(i);
    if (iglobal>=0) {
      tempPar[i] = p[ iglobal ];

      // check that pt is not unphysically large
      if ( i==iOmega ) {
        double ffomega = tempPar[iOmega]*parfact[iOmega];
        if ( fabs( ffomega ) < omegaMin ) {
          int signO = par[i]>0 ? 1 : -1;
          tempPar[i] = signO*omegaMin/parfact[iOmega];
          cout << "TrackParticleFitObject::updateParams INFO: regularising Omega to: " << tempPar[i] << endl;
        }
      }

      // check is there has been a significant parameter update
      if ( pow( tempPar[i] - par[i], 2) >  eps2 * cov[i][i] ) result=true; // check if any have been updated
      p[iglobal]=tempPar[i];  // update the global vars
      par[i]    =tempPar[i];  // update local variables
    }
  }

  return result;
}

ThreeVector TrackParticleFitObject::getTrackPlaneNormal() const {
  updateCache(); 
  return trackPlaneNormal;
}

ThreeVector TrackParticleFitObject::getTrackPcaVector() const {
  updateCache(); 
  return trackPcaVector;
}

double TrackParticleFitObject::getDPx(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  //  if (!cachevalid) 
  updateCache();
  return getMomentumFirstDerivatives(1, ilocal);
}

double TrackParticleFitObject::getDPy(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  //  if (!cachevalid) updateCache();
  updateCache();
  return getMomentumFirstDerivatives(2, ilocal);
}

double TrackParticleFitObject::getDPz(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  //if (!cachevalid) 
  updateCache();
  return getMomentumFirstDerivatives(3, ilocal);
}

double TrackParticleFitObject::getDE(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  //  if (!cachevalid) 
  updateCache();
  return getMomentumFirstDerivatives(0, ilocal);
}

double TrackParticleFitObject::getFirstDerivative_Meta_Local( int iMeta, int ilocal , int metaSet ) const {
  // iMeta = intermediate variable (i.e. E,px,py,pz)
  // ilocal = local variable (ptinv, theta, phi)
  // metaSet = which set of intermediate varlables
  //if (!cachevalid) 
  updateCache();
  switch ( metaSet ) {
  case BaseDefs::VARBASIS_EPXYZ:
    return getMomentumFirstDerivatives(iMeta, ilocal);
    break;
  case BaseDefs::VARBASIS_TRKNORMAL:
    return getNormalFirstDerivatives(iMeta, ilocal);
    break;

    // do we need the vertex stuff here?? if so, need to know which vertex we're talking about

  default:
    assert(0);
  }
  return 0; // should never get here
}

double TrackParticleFitObject::getSecondDerivative_Meta_Local( int iMeta, int ilocal , int jlocal , int metaSet ) const {
  //  if (!cachevalid) 
  updateCache();
  switch ( metaSet ) {
  case BaseDefs::VARBASIS_EPXYZ:
    return getMomentumSecondDerivatives(iMeta, ilocal, jlocal);
    break;
  case BaseDefs::VARBASIS_TRKNORMAL:
    return getNormalSecondDerivatives(iMeta, ilocal, jlocal);
    break;

    // do we need the vertex stuff here??

  default:
    assert(0);
  }
  return 0; // should never get here
}


double TrackParticleFitObject::getChi2 () const {
  //if (!cachevalid)
  updateCache();
  return chi2;
}

void TrackParticleFitObject::updateCache() const {

  if ( cachevalid ) return;

  //  std::cout << "TrackParticleFitObject::updateCache" << std::endl;

  chi2 = ParticleFitObject::getChi2 ();

  //  cout << "chisq = " << chi2 << endl;

  phi0  = getParam(iPhi0 )*parfact[iPhi0 ] ; // rescale to physical units
  omega = getParam(iOmega)*parfact[iOmega] ;
  tanl  = getParam(iTanL )*parfact[iTanL ] ;
  d0    = getParam(iD0   )*parfact[iD0   ] ;
  z0    = getParam(iZ0   )*parfact[iZ0   ] ;

  s_start= getParam(iStart)*parfact[iStart] ;

  s_end  = getParam(iEnd)*parfact[iEnd] ;

  //  cout << "TrackParticleFitObject::updateCache() " << chi2 << " " << phi0 << " " << omega << " " << tanl << " " << d0 << " " << z0 << " " << s_start << " " << s_end << endl;

  //  cout <<  getParam(iPhi0 ) << " " <<  getParam(iOmega) << " " <<  getParam(iTanL ) << " " <<  getParam(iD0   ) << " " <<  getParam(iZ0   ) << " " << getParam(iStart) << endl;

  double aB = omega_pt_conv*getBfield();
  double pt = aB/fabs( omega );
  double p  = pt * sqrt ( 1 + pow( tanl, 2 ) );

  // this is the 4mom at PCA (s=0)
  fourMomentum.setValues( sqrt ( pow( p, 2 ) + pow ( mass, 2 ) ) ,
                          pt*cos( phi0 ),
                          pt*sin( phi0 ),
                          pt*tanl );


  // momentum and point on trajectory at PCA, s_start, s_end
  //  ThreeVector trajectoryPointAtPCA;
  trajectoryPointAtPCA.setValues( trackReferencePoint.getX() + d0*sin(phi0),
				  trackReferencePoint.getY() + d0*cos(phi0),
				  trackReferencePoint.getZ() + z0 );
  momentumAtPCA.setValues( pt*cos( phi0 ), 
			   pt*sin( phi0 ), 
			   pt*tanl );

  double xx = omega*s_start/2.;
  double sincxx = fabs(xx)>1e-9 ? sin(xx)/xx : 1;
  trajectoryPointAtStart.setValues( trajectoryPointAtPCA.getX() + s_start*sincxx*cos(phi0 - xx),
				    trajectoryPointAtPCA.getY() + s_start*sincxx*sin(phi0 - xx),
				    trajectoryPointAtPCA.getZ() + s_start*tanl );
  
  double phiStart = phi0 - s_start*omega;
  //  ThreeVector momentumAtStart;
  momentumAtStart.setValues( pt*cos( phiStart ), 
			     pt*sin( phiStart ), 
			     pt*tanl );
  
  xx = omega*s_end/2.;
  sincxx = fabs(xx)>1e-9 ? sin(xx)/xx : 1;
  trajectoryPointAtEnd.setValues( trajectoryPointAtPCA.getX() + s_end*sincxx*cos(phi0 - xx),
				  trajectoryPointAtPCA.getY() + s_end*sincxx*sin(phi0 - xx),
				  trajectoryPointAtPCA.getZ() + s_end*tanl );
  
  double phiEnd = phi0 - s_end*omega;
  momentumAtEnd.setValues( pt*cos( phiEnd ), 
			   pt*sin( phiEnd ), 
			   pt*tanl );
  
  //  cout << "TrackParticleFitObject::updateCache : FourMomentum = " << fourMomentum << endl;
  //  cout << "TrackParticleFitObject::updateCache() inter1: " << chi2 << " " << phi0 << " " << omega << " " << tanl << " " << d0 << " " << z0 << endl;

  updateNormalDerivatives();

  //  cout << "TrackParticleFitObject::updateCache() inter2: " << chi2 << " " << phi0 << " " << omega << " " << tanl << " " << d0 << " " << z0 << endl;

  updateMomentumDerivatives();

  //  cout << "TrackParticleFitObject::updateCache() inter3: " << chi2 << " " << phi0 << " " << omega << " " << tanl << " " << d0 << " " << z0 << endl;

  updateTrajectoryDerivatives();

  //  cout << "TrackParticleFitObject::updateCache() done: " << chi2 << " " << phi0 << " " << omega << " " << tanl << " " << d0 << " " << z0 << endl;
  //  cout <<  getParam(iPhi0 ) << " " <<  getParam(iOmega) << " " <<  getParam(iTanL ) << " " <<  getParam(iD0   ) << " " <<  getParam(iZ0   ) << " " << getParam(iStart) << endl;
  //  cout << "Normal vector = " << trackPlaneNormal << " " << trackPlaneNormal.getMag() << endl;

  cachevalid = true;

  //  cout << "... updated cache " << cachevalid << endl;

  return;
}

void TrackParticleFitObject::updateTrajectoryDerivatives() const {
  resetTrajectoryFirstDerivatives();
  resetTrajectorySecondDerivatives();

//  cout << "hello from updateTrajectoryDerivatives" << endl;
//  cout <<  "parameters: " << getParam(iPhi0 ) << " " <<  getParam(iOmega) << " " <<  getParam(iTanL ) << " " <<  
//    getParam(iD0   ) << " " <<  getParam(iZ0   ) << " " << getParam(iStart) << " " << getParam(iEnd) << endl;

  const int nVars = 3;

  assert( nVars == BaseDefs::nMetaVars[BaseDefs::VARBASIS_VXYZ] );

  // local parameters
  //  enum {iD0=0, iPhi0, iOmega, iZ0, iTanL, iStart, NPAR};
  //  enum {iD0=0, iPhi0, iOmega, iZ0, iTanL, iStart, iEnd, NPAR}; <--- now include ending point

  // point on trajectory in terms of intermediate variables
  // nb sinc() makes "badly behaved" derivatives (with s in denominator...) when I use mathematica
  // so use param. with omega (W) - this will run into problems when W->0...but I (hope) this will not happen
  //                                if it does, may have to think of alternative approach

  // x = xr + d0*sin(phi0) + s*sinc(omega*s/2)*cos(phi0 - omega*s/2) = X + D*sin(P) + (2/W)*sin(W*S/2)*cos(P-W*S/2)
  // y = yr + d0*cos(phi0) + s*sinc(omega*s/2)*sin(phi0 - omega*s/2) = Y + D*cos(P) + (2/W)*sin(W*S/2)*sin(P-W*S/2)
  // z = zr + z0 + s*tanl = R + Z + S*T

  // intermediate/meta variables: d0, phi0, omega, z0, tanl, s
  const int int_D = 0;
  const int int_P = 1;
  const int int_W = 2;
  const int int_Z = 3;
  const int int_T = 4;
  const int int_S = 5;
  const int nInt  = 6;

  double D = d0;
  double P = phi0;
  double W = omega;
  double Z = z0;
  double T = tanl;
  double S(0);
  
  for (int iss=0; iss<2; iss++) { // start/end points of track

    //if ( iss==0 ) cout << "--start vertex" << endl;
    //else          cout << "--end vertex" << endl;

    S           = iss==0 ? s_start: s_end ;
    int iPoint  = iss==0 ? iStart : iEnd ;

    double interFirstDerivs [nInt][NPAR];
    double interSecondDerivs[nInt][NPAR][NPAR];

    for (int i=0; i<nInt; i++) {
      for (int j=0; j<NPAR; j++) {
	interFirstDerivs[i][j]=0;
	for (int k=0; k<NPAR; k++) {
	  interSecondDerivs[i][j][k]=0;
	}
      }
    }

    interFirstDerivs[int_D][iD0   ] = 1;
    interFirstDerivs[int_P][iPhi0 ] = 1;
    interFirstDerivs[int_W][iOmega] = 1;
    interFirstDerivs[int_Z][iZ0   ] = 1;
    interFirstDerivs[int_T][iTanL ] = 1;
    interFirstDerivs[int_S][iPoint] = 1;

    // cout << "interFirstDerivs" << endl;
    // for (int i=0; i<nInt; i++) {
    //   for (int j=0; j<NPAR; j++) {
    //     cout << interFirstDerivs[i][j] << " ";
    //   }
    //   cout << endl;
    // }

    double trajectoryInterFirstDerivs[nVars][nInt];
    double trajectoryInterSecondDerivs[nVars][nInt][nInt];
    for (int i=0; i<nVars; i++) {
      for (int j=0; j<nInt; j++) {
	trajectoryInterFirstDerivs[i][j]=0;
	for (int k=0; k<nInt; k++) {
	  trajectoryInterSecondDerivs[i][j][k]=0;
	}
      }
    }

    /*
x        = X + D Sin[P] + (2 Cos[P - (S W)/2] Sin[(S W)/2])/W

dx/dD    = Sin[P]
dx/dP    = (Cos[P] + D W Cos[P] - Cos[P - S W])/W
dx/dS    = Cos[P - S W]
dx/dW    = (S W Cos[P - S W] - Sin[P] + Sin[P - S W])/W^2

d2x/dDdD = 0
d2x/dDdP = Cos[P]
d2x/dDdS = 0
d2x/dDdW = 0

d2x/dPdD = Cos[P]
d2x/dPdP = (-(1 + D W) Sin[P] + Sin[P - S W])/W
d2x/dPdS = -Sin[P - S W]
d2x/dPdW = -((Cos[P] - Cos[P - S W] + S W Sin[P - S W])/W^2)

d2x/dSdD = 0
d2x/dSdP = -Sin[P - S W]
d2x/dSdS = W Sin[P - S W]
d2x/dSdW = S Sin[P - S W]

d2x/dWdD = 0
d2x/dWdP = -((Cos[P] - Cos[P - S W] + S W Sin[P - S W])/W^2)
d2x/dWdS = S Sin[P - S W]
d2x/dWdW = (-2 S W Cos[P - S W] + 2 Sin[P] + (-2 + S^2 W^2) Sin[P - S W])/W^3


    */

    double PmSW = P - S*W;

    trajectoryInterFirstDerivs[0][int_D]         = sin(P);
    trajectoryInterFirstDerivs[0][int_P]         = (cos(P) + D*W*cos(P) - cos(PmSW))/W;
    trajectoryInterFirstDerivs[0][int_S]         = cos(PmSW);
    trajectoryInterFirstDerivs[0][int_W]         = (S*W*cos(PmSW) - sin(P) + sin(PmSW))/pow(W,2);

    trajectoryInterSecondDerivs[0][int_D][int_D] = 0;
    trajectoryInterSecondDerivs[0][int_D][int_P] = cos(P);
    trajectoryInterSecondDerivs[0][int_D][int_S] = 0;
    trajectoryInterSecondDerivs[0][int_D][int_W] = 0;

    trajectoryInterSecondDerivs[0][int_P][int_D] = cos(P);
    trajectoryInterSecondDerivs[0][int_P][int_P] = (-(1 + D*W)*sin(P) + sin(PmSW))/W;
    trajectoryInterSecondDerivs[0][int_P][int_S] = -sin(PmSW);
    trajectoryInterSecondDerivs[0][int_P][int_W] = -((cos(P) - cos(PmSW) + S*W*sin(PmSW))/pow(W,2));

    trajectoryInterSecondDerivs[0][int_S][int_D] = 0;
    trajectoryInterSecondDerivs[0][int_S][int_P] =  -sin(PmSW);
    trajectoryInterSecondDerivs[0][int_S][int_S] = W*sin(PmSW);
    trajectoryInterSecondDerivs[0][int_S][int_W] = S*sin(PmSW);

    trajectoryInterSecondDerivs[0][int_W][int_D] = 0;
    trajectoryInterSecondDerivs[0][int_W][int_P] = -((cos(P) - cos(PmSW) + S*W*sin(PmSW))/pow(W,2));
    trajectoryInterSecondDerivs[0][int_W][int_S] = S*sin(PmSW);
    trajectoryInterSecondDerivs[0][int_W][int_W] = ( -2*S*W*cos(PmSW) + 2*sin(P) + (-2 + pow(S,2)*pow(W,2))*sin(PmSW) ) / pow(W,3);


    /*

y        = Y + D Cos[P] + (2 Sin[(S W)/2] Sin[P - (S W)/2])/W
           
dy/dD    = Cos[P]
dy/dP    = -(((-1 + D W) Sin[P] + Sin[P - S W])/W)
dy/dS    = Sin[P - S W]
dy/dW    = (Cos[P] - Cos[P - S W] + S W Sin[P - S W])/W^2
           
d2y/dDdD = 0
d2y/dDdP = -Sin[P]
d2y/dDdS = 0
d2y/dDdW = 0
           
d2y/dPdD = -Sin[P]
d2y/dPdP = -(((-1 + D W) Cos[P] + Cos[P - S W])/W)
d2y/dPdS = Cos[P - S W]
d2y/dPdW = (S W Cos[P - S W] - Sin[P] + Sin[P - S W])/W^2
           
d2y/dSdD = 0
d2y/dSdP = Cos[P - S W]
d2y/dSdS = -W Cos[P - S W]
d2y/dSdW = -S Cos[P - S W]
           
d2y/dWdD = 0
d2y/dWdP = (S W Cos[P - S W] - Sin[P] + Sin[P - S W])/W^2
d2y/dWdS = -S Cos[P - S W]
d2y/dWdW = (-2 Cos[P] + (2 - S^2 W^2) Cos[P - S W] - 2 S W Sin[P - S W])/W^3

     */

    trajectoryInterFirstDerivs[1][int_D]         =  cos(P);
    trajectoryInterFirstDerivs[1][int_P]         =  -(((-1 + D*W)*sin(P) + sin(PmSW))/W);
    trajectoryInterFirstDerivs[1][int_S]         =  sin(PmSW);
    trajectoryInterFirstDerivs[1][int_W]         =  (cos(P) - cos(PmSW) + S*W*sin(PmSW))/pow(W,2);
                                                                                                                 
    trajectoryInterSecondDerivs[1][int_D][int_D] =  0;
    trajectoryInterSecondDerivs[1][int_D][int_P] =  -sin(P);
    trajectoryInterSecondDerivs[1][int_D][int_S] =  0;
    trajectoryInterSecondDerivs[1][int_D][int_W] =  0;
                                                                                                                 
    trajectoryInterSecondDerivs[1][int_P][int_D] =  -sin(P);
    trajectoryInterSecondDerivs[1][int_P][int_P] =  -(((-1 + D*W)*cos(P) + cos(PmSW))/W);
    trajectoryInterSecondDerivs[1][int_P][int_S] =  cos(PmSW);
    trajectoryInterSecondDerivs[1][int_P][int_W] =  (S*W*cos(PmSW) - sin(P) + sin(PmSW))/pow(W,2);
     
    trajectoryInterSecondDerivs[1][int_S][int_D] =  0;
    trajectoryInterSecondDerivs[1][int_S][int_P] =  cos(PmSW); 
    trajectoryInterSecondDerivs[1][int_S][int_S] =  -W*cos(PmSW);
    trajectoryInterSecondDerivs[1][int_S][int_W] =  -S*cos(PmSW);

    trajectoryInterSecondDerivs[1][int_W][int_D] =  0;
    trajectoryInterSecondDerivs[1][int_W][int_P] =  (S*W*cos(PmSW) - sin(P) + sin(PmSW))/pow(W,2);
    trajectoryInterSecondDerivs[1][int_W][int_S] =  -S*cos(PmSW);
    trajectoryInterSecondDerivs[1][int_W][int_W] =  (-2*cos(P) + (2 - pow(S,2)*pow(W,2))*cos(PmSW) - 2*S*W*sin(PmSW))/pow(W,3);

    /*

z = R + S T + Z

dz/dZ = 1
dz/dS = T
dz/dT = S

d2z/dZdZ = 0
d2z/dZdS = 0
d2z/dZdT = 0

d2z/dSdZ = 0
d2z/dSdS = 0
d2z/dSdT = 1

d2z/dTdZ = 0
d2z/dTdS = 1
d2z/dTdT = 0

    */

    trajectoryInterFirstDerivs[2][int_Z] = 1;
    trajectoryInterFirstDerivs[2][int_S] = T;
    trajectoryInterFirstDerivs[2][int_T] = S;

    trajectoryInterSecondDerivs[2][int_Z][int_Z] = 0;
    trajectoryInterSecondDerivs[2][int_Z][int_S] = 0;
    trajectoryInterSecondDerivs[2][int_Z][int_T] = 0;

    trajectoryInterSecondDerivs[2][int_S][int_Z] = 0;
    trajectoryInterSecondDerivs[2][int_S][int_S] = 0;
    trajectoryInterSecondDerivs[2][int_S][int_T] = 1;

    trajectoryInterSecondDerivs[2][int_T][int_Z] = 0;
    trajectoryInterSecondDerivs[2][int_T][int_S] = 1;
    trajectoryInterSecondDerivs[2][int_T][int_T] = 0;

    // cout << "trajectoryInterFirstDerivs" << endl;
    // for (int i=0; i<nVars; i++) {
    //   for (int j=0; j<nInt; j++) {
    // 	cout << std::setw(10) << trajectoryInterFirstDerivs[i][j] << " ";
    //   }
    //   cout << endl;
    // }



    // now calculate the total derivatives of xyz wrt trk params using chain rule
    //
    for (int ipe=0; ipe<nVars; ipe++) { // the 3 xyz
      for (int ipar=0; ipar<NPAR; ipar++) { // 7 track parameters
	// the first derivs (chain rule)
	double dd(0);
	for (int j=0; j<nInt; j++) { // 6 intermediate variables
	  dd+=trajectoryInterFirstDerivs[ipe][j]*interFirstDerivs[j][ipar];
	}
	if ( iss==0 ) {
	  setTrajectoryStartFirstDerivatives(ipe, ipar, dd);
	} else if ( iss==1 ) { 
	  setTrajectoryEndFirstDerivatives(ipe, ipar, dd);
	} else {
	  assert(0);
	}

	// the second derivs
	for (int jpar=0; jpar<NPAR; jpar++) {
	  double dd2(0);
	  for (int j=0; j<nInt; j++) {
	    dd2+=trajectoryInterFirstDerivs[ipe][j]*interSecondDerivs[j][ipar][jpar];
	    for (int k=0; k<nInt; k++) {
	      dd2+=trajectoryInterSecondDerivs[ipe][j][k]*interFirstDerivs[j][ipar]*interFirstDerivs[k][jpar];
	    }
	  }
	  if ( iss==0 ) {
	    setTrajectoryStartSecondDerivatives(ipe, ipar, jpar, dd2);	
	  } else {
	    setTrajectoryEndSecondDerivatives(ipe, ipar, jpar, dd2);	
	  }
	}
      }
    }

  } // iss ; start/end vertex

  //cout << "bye from updateTrajectoryDerivatives" << endl;
  //cout << "first derivatives (START): " << endl;
  //for (int ipe=0; ipe<nVars; ipe++) {
  //  for (int ipar=0; ipar<NPAR; ipar++) {
  //    cout << std::setw(10) << getTrajectoryStartFirstDerivatives(ipe, ipar) << " ";
  //  }
  //  cout << endl;
  //}
  //cout << "first derivatives (END): " << endl;
  //for (int ipe=0; ipe<nVars; ipe++) {
  //  for (int ipar=0; ipar<NPAR; ipar++) {
  //    cout << std::setw(10) << getTrajectoryEndFirstDerivatives(ipe, ipar) << " ";
  //  }
  //  cout << endl;
  //}

  return;
}

void TrackParticleFitObject::updateMomentumDerivatives() const {

  // daniel adding start variable

  // -------------------------------
  // the momentum derivatives
  // -------------------------------
  resetMomentumFirstDerivatives();
  resetMomentumSecondDerivatives();

  /*
    pt = aB/|omega|

    p  = pt ( 1 + tanl^2 )

    e = ( p^2 + m^2 )^(1/2)
    px = p cos (phi) ( 1 + tanl^2 )^-1
    py = p sin (phi) ( 1 + tanl^2 )^-1
    pz = p tanl ( 1 + tanl^2 )^-1

    (tanl = t)

    intermediate vars: p, phi, tanl

    orig vars: omega, phi0, tanl, z0, d0, start

    p  = aB ( 1 + t^2 ) / |omega|

    dp / dt = aB * 2 * t / |omega|
    dp / domega = - sign(omega) aB ( 1 + t^2 ) / |omega|^2

    d2p / dt2 = 2 aB / |omega|
    d2p / domega dt = - sign (omega) * 2 ab t / |omega|^2
    d2p / domega2 = 2 ab ( 1 + t^2 ) / |omega|^3

  */

  // intermediate vars: p=0, phi=1, tanl=2
  const int iP = 0;
  const int iPh= 1;
  const int iT = 2;
  const int nInt = 3;

  // track vars: iD0=0, iPhi0, iOmega, iZ0, iTanL, NPAR

  double aB = omega_pt_conv*getBfield();
  double pt = aB/fabs( omega );
  double p  = pt * sqrt ( 1 + pow( tanl, 2 ) );
  double e = sqrt( p*p + mass*mass );
  double one_tan2 = 1 + pow(tanl,2);
  
  double interFirstDerivs[nInt][NPAR];
  double interSecondDerivs[nInt][NPAR][NPAR];
  for (int i=0; i<nInt; i++) {
    for (int j=0; j<NPAR; j++) {
      interFirstDerivs[i][j]=0;
      for (int k=0; k<NPAR; k++) {
	interSecondDerivs[i][j][k]=0;
      }
    }
  }

  int osign = omega>0 ? +1 : -1 ;

  interFirstDerivs[iP][iOmega] = - osign * aB * one_tan2 / pow(omega,2) ;
  interFirstDerivs[iP][iTanL]  = 2*aB*tanl/fabs(omega);

  interFirstDerivs[iPh][iOmega] = -s_start; // daniel added
  interFirstDerivs[iPh][iPhi0]  = 1;
  interFirstDerivs[iPh][iStart] = -omega; // daniel added

  interFirstDerivs[iT][iTanL] = 1;

  //  interSecondDerivs[iP][iTanL ][iTanL ] = 2*aB*(1+pow(tanl,2))/pow(fabs(omega),2);
  interSecondDerivs[iP][iTanL ][iTanL ] = 2*aB/fabs(omega); // DJeans fixed 28May2015
  interSecondDerivs[iP][iOmega][iTanL ] = interSecondDerivs[iP][iTanL][iOmega] = -osign*2*aB*tanl/pow(omega,2);
  interSecondDerivs[iP][iOmega][iOmega] = 2*aB*one_tan2/pow(fabs(omega),3);

  interSecondDerivs[iPh][iOmega][iStart] = interSecondDerivs[iPh][iStart][iOmega] = -1; // daniel added


  // intermediate vars: p=0, phi=1, tanl=2
  double momentumInterFirstDerivs[4][nInt];
  double momentumInterSecondDerivs[4][nInt][nInt];

  for (int i=0; i<4; i++) {
    for (int j=0; j<nInt; j++) {
      momentumInterFirstDerivs[i][j]=0;
      for (int k=0; k<nInt; k++) {
	momentumInterSecondDerivs[i][j][k]=0;
      }
    }
  }

  /*
    e = ( p^2 + m^2 )^(1/2)
    de / dp = (1/2) 2p ( p^2 + m^2 )^(-1/2) = p ( p^2 + m^2 )^(-1/2)
    de / dphi = 0
    de / dt = 0

    d2e / dp2 = ( p^2 + m^2 )^(-1/2) - p*(1/2)*2p*( p^2 + m^2 )^(-3/2) = ( p^2 + m^2 )^(-1/2) - p^2 ( p^2 + m^2 )^(-3/2)
  */

  momentumInterFirstDerivs[0][iP] = p/e;
  momentumInterSecondDerivs[0][iP][iP] = 1./e - pow(p,2)/pow(e,3);

  /*
    px = p cos (phi) ( 1 + t^2 )^-1
    dpx / dp = cos (phi) ( 1 + t^2 )^-1
    dpx / dphi = -p sin(phi) ( 1 + t^2 )^-1
    dpx / dt = p cos (phi) * -2 t (1+t^2)^-2

    d2px/ dp2 = 0
    d2px/ dphi dp  = -sin (phi) ( 1 + t^2 )^-1
    d2px/ dt   dp  = -2t*cos (phi) ( 1 + t^2 )^-2
    d2px/ dp dphi  = - sin(phi) ( 1 + t^2 )^-1
    d2px/ dphi2 = -p cos(phi) ( 1 + t^2 )^-1
    d2px/ dt dphi  = p sin(phi) *2t * ( 1 + t^2 )^-2
    d2px / dp dt   = - cos (phi) * 2 t (1+t^2)^-2
    d2px / dphi dt = p sin(phi) * 2 t (1+t^2)^-2
    d2px / dt2     = -2 p cos (phi) (   (1+t^2)^-2  - 4 t^2 (1+t^2)^-3 )
  */

  momentumInterFirstDerivs[1][iP] = cos(phi0) / one_tan2;
  momentumInterFirstDerivs[1][iPh] = -p*sin(phi0)/one_tan2;
  momentumInterFirstDerivs[1][iT] = -2*p*tanl*cos(phi0)/pow(one_tan2,2);

  momentumInterSecondDerivs[1][iP ][iPh] = momentumInterSecondDerivs[1][iPh][iP ] = -sin(phi0)/one_tan2;
  momentumInterSecondDerivs[1][iP ][iT ] = momentumInterSecondDerivs[1][iT ][iP ] = -2*tanl*cos(phi0)/pow(one_tan2,2);
  momentumInterSecondDerivs[1][iPh][iPh] = -p*cos(phi0)/one_tan2;
  momentumInterSecondDerivs[1][iP ][iT ] = momentumInterSecondDerivs[1][iT ][iP ] = -2*tanl*cos(phi0)/pow(one_tan2,2);
  momentumInterSecondDerivs[1][iPh][iT ] = momentumInterSecondDerivs[1][iT ][iPh] = 2*tanl*p*sin(phi0)/pow(one_tan2,2);
  momentumInterSecondDerivs[1][iT ][iT ] = -2*p*cos(phi0)*( 1./pow(one_tan2,2) - 4*pow(tanl,2)/pow(one_tan2,3) );

  /*
    py = p sin (phi) ( 1 + t^2 )^-1
    dpy / dp = sin (phi) ( 1 + t^2 )^-1
    dpy / dphi = p cos (phi) ( 1 + t^2 )^-1
    dpy / dt = p sin (phi) * -2 t (1+t^2)^-2

    d2py / dp2 = 0
    d2py / dphi dp = cos(phi) ( 1 + t^2 )^-1
    d2py / dt dp   = -2t sin (phi) ( 1 + t^2 )^-2

    d2py / dp dphi = cos (phi) ( 1 + t^2 )^-1
    d2py / dphi2   = -p sin (phi) ( 1 + t^2 )^-1
    d2py / dt dphi = -2t p cos (phi) ( 1 + t^2 )^-2

    d2py / dp dt = sin (phi) * -2 t (1+t^2)^-2
    d2py / dphi dt = p cos (phi) * -2 t (1+t^2)^-2
    d2py / dt2 = -2 p sin (phi) ( (1+t^2)^-2 - 4 t^2 (1+t^2)^-3 )
  */

  momentumInterFirstDerivs[2][iP ] = sin(phi0) / one_tan2;
  momentumInterFirstDerivs[2][iPh] = p*cos(phi0)/one_tan2;
  momentumInterFirstDerivs[2][iT ] = -2*p*tanl*sin(phi0)/pow(one_tan2,2);

  momentumInterSecondDerivs[2][iP ][iPh] = momentumInterSecondDerivs[1][iPh][iP ] = cos(phi0)/one_tan2;
  momentumInterSecondDerivs[2][iP ][iT ] = momentumInterSecondDerivs[1][iT ][iP ] = -2*tanl*sin(phi0)/pow(one_tan2,2);
  momentumInterSecondDerivs[2][iPh][iPh] = -p*sin(phi0)/one_tan2;
  momentumInterSecondDerivs[2][iP ][iT ] = momentumInterSecondDerivs[1][iT ][iP ] = -2*tanl*sin(phi0)/pow(one_tan2,2);
  momentumInterSecondDerivs[2][iPh][iT ] = momentumInterSecondDerivs[1][iT ][iPh] = -2*tanl*p*cos(phi0)/pow(one_tan2,2);
  momentumInterSecondDerivs[2][iT ][iT ] = -2*p*sin(phi0)*( 1./pow(one_tan2,2) - 4*pow(tanl,2)/pow(one_tan2,3) );

  /*

    pz = p t ( 1 + t^2 )^-1
    dpz / dp = t ( 1 + t^2 )^-1
    dpz / dphi = 0
    dpz / dt = p ( 1 + t^2 )^-1 + p*t*-2t*(1+t^2)-2 = p ( ( 1 + t^2 )^-1 - 2t^2 (1+t^2)^-2 );

    d2pz / dp2 = 0
    d2pz / dphi dp = 0
    d2pz / dt dp = ( 1 + t^2 )^-1 - 2 t^2 ( 1 + t^2 )^-2

    d2pz / dp dt = ( ( 1 + t^2 )^-1 - 2t^2 (1+t^2)^-2 )
    d2pz / dt2 = p ( -2t( 1 + t^2 )^-2 - 4t(1+t^2)^-2 + 8 t^3 (1+t^2)^-3 ) = p ( -6t (1+t^2)^-2 + 8 t^3 (1+t^2)^-3 )


  */

  momentumInterFirstDerivs[3][iP ] = tanl/one_tan2;
  momentumInterFirstDerivs[3][iT ] = p*( 1./one_tan2 - 2*pow(tanl/one_tan2,2) );

  momentumInterSecondDerivs[3][iP ][iT ] = momentumInterSecondDerivs[3][iT ][iP ] = 
    1./one_tan2 - 2.*pow(tanl/one_tan2, 2);
  momentumInterSecondDerivs[3][iT ][iT ] = p*( -6.*tanl/pow(one_tan2,2) + 8.*pow(tanl/one_tan2, 3) );

  // now calculate the total derivatives of Epxpypz wrt trk params using chain rule
  for (int ipe=0; ipe<4; ipe++) {
    for (int ipar=0; ipar<NPAR; ipar++) {
      // the first derivs (chain rule)
      double dd(0);
      for (int j=0; j<nInt; j++) {
	dd+=momentumInterFirstDerivs[ipe][j]*interFirstDerivs[j][ipar];
      }
      setMomentumFirstDerivatives(ipe, ipar, dd);
      // the second derivs
      for (int jpar=0; jpar<NPAR; jpar++) {
	double dd2(0);
	for (int j=0; j<nInt; j++) {
	  dd2+=momentumInterFirstDerivs[ipe][j]*interSecondDerivs[j][ipar][jpar];
	  for (int k=0; k<nInt; k++) {
	    dd2+=momentumInterSecondDerivs[ipe][j][k]*interFirstDerivs[j][ipar]*interFirstDerivs[k][jpar];
	  }
	}
	setMomentumSecondDerivatives(ipe, ipar, jpar, dd2);	
      }
    }
  }

  return;
}


void TrackParticleFitObject::updateNormalDerivatives() const {
  // ------------------------------
  // the derivatives of normal to track-IP plane wrt track parameters
  // ------------------------------
  // this is rather messy
  // we use the chain rule to simplify somewhat
  // with intermediate parameters (a,b,c), proportional to the
  // components of the cross product of the line from IP->PCA and the momentum vector @ PCA

  // as of June 2015, this part is not thoroughly tested. DJeans.

  resetNormalFirstDerivatives();
  resetNormalSecondDerivatives();

  // (x,y,z) is the reference point of the track parameters
  double x = trackReferencePoint.getX();
  double y = trackReferencePoint.getY();
  double z = trackReferencePoint.getZ();

  // vector from ref point -> PCA is ( -d0 sin(phi), d0 cos(phi), z0 )
  // PCA vector: PCA = (x,y,z) + ( -d0 sin(phi), d0 cos(phi), z0 )
  // momentum 3-vector at PCA: MOM = pt*( cos(phi0), sin(phi0), tanl )

  trackPcaVector.setValues( x - d0*sin(phi0) , y + d0*cos(phi0) , z + z0 );

  //  cout << "fikka: refpt " << x << " " << y << " " << z << endl;
  //  cout << "4-mom " << fourMomentum << endl;
  //  cout << "PCA   " << trackPcaVector << endl;
  //  cout << "dot prod " << fourMomentum.getThreeVector() * trackPcaVector << endl;


  //  ThreeVector pca
  //  cout << "check PCA: d0 " << d0 << " phi0 " << phi0 << " pca:" << pca << endl;

  //  cout << "TrackParticleFitObject::updateNormalDerivatives : pca " << pca << endl;

  //PCA cross MOM
  //  = pt ( (x-d0 sin(phi), y+d0 cos(phi), z+z0) cross ( cos(phi), sin(phi), tanl ) )
  //  = pt ( a , b , c ) <--- a,b,c are "intermediate" parameters, used in applying chain rule
  //
  //a = (y+d0 cos(phi))*tanl - (z+z0)*sin(phi)
  //b = (z+z0)*cos(phi) - (x-d0 sin(phi))*tanl
  //c = (x-d0 sin(phi))*sin(phi) - (y+d0 cos(phi))*cos(phi)
  //  = x sin(phi) - y cos(phi) - d0 (sin2(phi) + cos2(phi) )
  //  = x sin(phi) - y cos(phi) - d0

  double ABC[3] = { (y+d0*cos(phi0))*tanl - (z+z0)*sin(phi0),
                    (z+z0)*cos(phi0) - (x-d0*sin(phi0))*tanl,
                    x*sin(phi0) - y*cos(phi0) - d0};

  // cout << "checking normal... " << x << " " << y << " " << z << " , " << 
  //   d0 << " " << z0 << " " << phi0 << " " << tanl << " , " << 
  //   ABC[0] << " " << ABC[1] <<  " " << ABC[2] << endl;

  // this is the ThreeVector perpendicular to the plane defined by IP, PCA, and track momentum @ PCA
  trackPlaneNormal.setValues( ABC[0], ABC[1], ABC[2] );
  trackPlaneNormal*=1./trackPlaneNormal.getMag();

  //  cout << "TrackParticleFitObject::updateNormalDerivatives : trackPlaneNormal " << trackPlaneNormal << endl;
  

  //
  // first and second derivatives of intermediate parameters abc wrt track parameters
  //
  //da/d(d0)   = cos(phi)*tanl
  //da/d(z0)   = -sin(phi)
  //da/d(phi)  = -d0 sin(phi)*tanl - (z+z0)*cos(phi)
  //da/d(tanl) = (y+d0 cos(phi))
  //
  //db/d(d0)   = sin(phi)*tanl
  //bd/d(z0)   = cos(phi)
  //db/d(phi)  = -(z+z0)*sin(phi) + d0 cos(phi)*tanl
  //db/d(tanl) = - (x-d0 sin(phi))
  //
  //dc/d(d0)   = -1
  //bc/d(z0)   = 0
  //dc/d(phi)  = x cos(phi) + y sin(phi)
  //dc/d(tanl) = 0

  double ABCderivs[3][NPAR];
  double ABCsecondderivs[3][NPAR][NPAR];

  ABCderivs[0][iPhi0 ]= -d0*sin(phi0)*tanl - (z+z0)*cos(phi0);
  ABCderivs[0][iOmega]= 0;
  ABCderivs[0][iTanL ]= y + d0*cos(phi0);
  ABCderivs[0][iD0   ]= cos(phi0)*tanl;
  ABCderivs[0][iZ0   ]= -sin(phi0);

  for (int i=0; i<3; i++)
    for (int j=0; j<NPAR; j++)
      for (int k=0; k<NPAR; k++)
        ABCsecondderivs[i][j][k]=0;

  ABCsecondderivs[0][iPhi0 ][iPhi0 ] = -d0*cos(phi0)*tanl + (z+z0)*sin(phi0);
  ABCsecondderivs[0][iPhi0 ][iTanL ] = -d0*sin(phi0);
  ABCsecondderivs[0][iPhi0 ][iD0   ] = -sin(phi0)*tanl;
  ABCsecondderivs[0][iPhi0 ][iZ0   ] = -cos(phi0);

  ABCsecondderivs[0][iTanL ][iPhi0 ] = -d0*sin(phi0);
  ABCsecondderivs[0][iTanL ][iD0   ] = cos(phi0);

  ABCsecondderivs[0][iD0   ][iPhi0 ] = -sin(phi0)*tanl;
  ABCsecondderivs[0][iD0   ][iTanL ] = cos(phi0);

  ABCsecondderivs[0][iZ0   ][iPhi0 ] = -cos(phi0);

  ABCderivs[1][iPhi0 ]= -(z+z0)*sin(phi0) + d0*cos(phi0)*tanl;
  ABCderivs[1][iOmega]= 0;
  ABCderivs[1][iTanL ]= -(x-d0*sin(phi0));
  ABCderivs[1][iD0   ]= sin(phi0)*tanl;
  ABCderivs[1][iZ0   ]= cos(phi0);

  ABCsecondderivs[1][iPhi0 ][iPhi0 ] = -(z+z0)*cos(phi0) - d0*sin(phi0)*tanl;
  ABCsecondderivs[1][iPhi0 ][iTanL ] = d0*cos(phi0);
  ABCsecondderivs[1][iPhi0 ][iD0   ] = cos(phi0)*tanl;
  ABCsecondderivs[1][iPhi0 ][iZ0   ] = -sin(phi0);

  ABCsecondderivs[1][iTanL ][iPhi0 ] = d0*cos(phi0);
  ABCsecondderivs[1][iTanL ][iD0   ] = sin(phi0);

  ABCsecondderivs[1][iD0   ][iPhi0 ] = cos(phi0)*tanl;
  ABCsecondderivs[1][iD0   ][iTanL ] = sin(phi0);

  ABCsecondderivs[1][iZ0   ][iPhi0 ] = -sin(phi0);


  ABCderivs[2][iPhi0 ]= x*cos(phi0) + y*sin(phi0);
  ABCderivs[2][iOmega]= 0;
  ABCderivs[2][iTanL ]= 0;
  ABCderivs[2][iD0   ]= -1;
  ABCderivs[2][iZ0   ]= 0;

  ABCsecondderivs[2][iPhi0 ][iPhi0 ] = -x*sin(phi0) + y*cos(phi0);


  //
  // derivatives of the normal vector N wrt the intermediate parameters abc
  //

  double NderivsABC[3][3];

  double sqabc = sqrt( pow(ABC[0],2) + pow(ABC[1],2) + pow(ABC[2],2) );

  for (int j=0; j<3; j++) { // <---- a,b,c
    for (int i=0; i<3; i++) { // <--- 3-vector
      NderivsABC[j][i] = 0;
      if (i==j) NderivsABC[j][i]+=1./sqabc;
      NderivsABC[j][i]-=ABC[j]*ABC[i]/pow(sqabc,3);
    }
  }

  // the normal's second derivatives wrt ABC
  // a little messy. i think this is ok...
  double NsecondderivsABC[3][3][3];
  for (int i=0; i<3; i++) { // <---- a,b,c
    for (int j=0; j<3; j++) { // <---- a,b,c
      for (int k=0; k<3; k++) { // <--- 3-vector
        NsecondderivsABC[i][j][k]=0;
        NsecondderivsABC[i][j][k]+=3*ABC[i]*ABC[j]*ABC[k]/pow(sqabc,5);
        if ( k==i ) {
          NsecondderivsABC[i][j][k]-=ABC[j]/pow(sqabc,3);
        }
        if ( k==j ) {
          NsecondderivsABC[i][j][k]-=ABC[i]/pow(sqabc,3);
        }
        if ( i==j ) {
          NsecondderivsABC[i][j][k]-=ABC[k]/pow(sqabc,3);
        }
      }
    }
  }


  // now sum over abc to get the derivatives of the normal vector N wrt the track parameters.

  // derivatives of vector wrt to track parameters (chain rule)
  //dN/d(d0) = dN/da da/dd0 + dN/db db/dd0 + dN/dc dc/dd0
  for (int i=0; i<NPAR; i++) { // <-- the object's parameters
    for (int j=0; j<3; j++) {  // <-- the normal's parameters
      double totsum(0);
      for (int k=0; k<3; k++) { // <-- sum over intermediate ABC params
        totsum+=NderivsABC[k][j]*ABCderivs[k][i];
      }
      setNormalFirstDerivatives(j,i,totsum);
    }
  }


  for (int i=0; i<NPAR; i++) { // <-- the object's parameters1
    for (int j=0; j<NPAR; j++) { // <-- the object's parameters2
      for (int m=0; m<3; m++) { // <-- three vector
        double sumtot(0);
        for (int k=0; k<3; k++) { // <-- sum over intermediate ABC params
          sumtot+=NderivsABC[k][m]*ABCsecondderivs[k][i][j];
          for (int l=0; l<3; l++) { // <-- sum over intermediate ABC params
            sumtot += NsecondderivsABC[k][l][m]*ABCderivs[k][i]*ABCderivs[l][j];
          }
        }
        setNormalSecondDerivatives(m, i, j, sumtot);
      }
    }
  }

  return;
}


void TrackParticleFitObject::setNormalFirstDerivatives(int i, int j, double x) const {
  assert ( i>=0 && i<3 && j>=0 && j<NPAR );
  normalFirstDerivatives[i][j]=x * parfact[j];
}

void TrackParticleFitObject::setNormalSecondDerivatives(int i, int j, int k, double x) const {
  assert ( i>=0 && i<3 && j>=0 && j<NPAR && k>=0 && k<NPAR);
  normalSecondDerivatives[i][j][k]=x * parfact[j] * parfact[k];
  normalSecondDerivatives[i][k][j]=normalSecondDerivatives[i][j][k];
}

void TrackParticleFitObject::setMomentumFirstDerivatives(int i, int j, double x) const {
  assert ( i>=0 && i<4 && j>=0 && j<NPAR );
  momentumFirstDerivatives[i][j]=x * parfact[j];
}

void TrackParticleFitObject::setMomentumSecondDerivatives(int i, int j, int k, double x) const {
  assert ( i>=0 && i<4 && j>=0 && j<NPAR && k>=0 && k<NPAR);
  momentumSecondDerivatives[i][j][k]=x * parfact[j] * parfact[k];
  momentumSecondDerivatives[i][k][j]=momentumSecondDerivatives[i][j][k];
}

void TrackParticleFitObject::setTrajectoryStartFirstDerivatives(int i, int j, double x) const {
  assert ( i>=0 && i<3 && j>=0 && j<NPAR );
  trajectoryStartFirstDerivatives[i][j]=x * parfact[j];
}

void TrackParticleFitObject::setTrajectoryStartSecondDerivatives(int i, int j, int k, double x) const {
  assert ( i>=0 && i<3 && j>=0 && j<NPAR && k>=0 && k<NPAR);
  trajectoryStartSecondDerivatives[i][j][k]=x * parfact[j] * parfact[k];
  trajectoryStartSecondDerivatives[i][k][j]=trajectoryStartSecondDerivatives[i][j][k];
}

void TrackParticleFitObject::setTrajectoryEndFirstDerivatives(int i, int j, double x) const {
  assert ( i>=0 && i<3 && j>=0 && j<NPAR );
  trajectoryEndFirstDerivatives[i][j]=x * parfact[j];
}

void TrackParticleFitObject::setTrajectoryEndSecondDerivatives(int i, int j, int k, double x) const {
  assert ( i>=0 && i<3 && j>=0 && j<NPAR && k>=0 && k<NPAR);
  trajectoryEndSecondDerivatives[i][j][k]=x * parfact[j] * parfact[k];
  trajectoryEndSecondDerivatives[i][k][j]=trajectoryEndSecondDerivatives[i][j][k];
}

double TrackParticleFitObject::getNormalFirstDerivatives(int i, int j) const {
  assert ( i>=0 && i<3 && j>=0 && j<NPAR );
  return normalFirstDerivatives[i][j];
}

double TrackParticleFitObject::getNormalSecondDerivatives(int i, int j, int k) const {
  assert ( i>=0 && i<3 && j>=0 && j<NPAR && k>=0 && k<NPAR);
  return normalSecondDerivatives[i][j][k];
}

double TrackParticleFitObject::getMomentumFirstDerivatives(int i, int j) const {
  assert ( i>=0 && i<4 && j>=0 && j<NPAR );
  return momentumFirstDerivatives[i][j];
}

double TrackParticleFitObject::getMomentumSecondDerivatives(int i, int j, int k) const {
  assert ( i>=0 && i<4 && j>=0 && j<NPAR && k>=0 && k<NPAR);
  return momentumSecondDerivatives[i][j][k];
}

double TrackParticleFitObject::getTrajectoryStartFirstDerivatives(int i, int j) const {
  assert ( i>=0 && i<3 && j>=0 && j<NPAR );
  return trajectoryStartFirstDerivatives[i][j];
}

double TrackParticleFitObject::getTrajectoryStartSecondDerivatives(int i, int j, int k) const {
  assert ( i>=0 && i<3 && j>=0 && j<NPAR && k>=0 && k<NPAR);
  return trajectoryStartSecondDerivatives[i][j][k];
}

double TrackParticleFitObject::getTrajectoryEndFirstDerivatives(int i, int j) const {
  assert ( i>=0 && i<3 && j>=0 && j<NPAR );
  return trajectoryEndFirstDerivatives[i][j];
}

double TrackParticleFitObject::getTrajectoryEndSecondDerivatives(int i, int j, int k) const {
  assert ( i>=0 && i<3 && j>=0 && j<NPAR && k>=0 && k<NPAR);
  return trajectoryEndSecondDerivatives[i][j][k];
}

void TrackParticleFitObject::resetMomentumFirstDerivatives() const {
  for (int i=0; i<4; i++)
    for (int j=0; j<NPAR; j++)
      momentumFirstDerivatives[i][j]=0;
  return;
}
void TrackParticleFitObject::resetMomentumSecondDerivatives() const {
  for (int i=0; i<4; i++)
    for (int j=0; j<NPAR; j++)
      for (int k=0; k<NPAR; k++)
        momentumSecondDerivatives[i][j][k]=0;
  return;
}

void TrackParticleFitObject::resetTrajectoryFirstDerivatives() const {
  for (int i=0; i<3; i++)
    for (int j=0; j<NPAR; j++) {
      trajectoryStartFirstDerivatives[i][j]=0;
      trajectoryEndFirstDerivatives[i][j]=0;
    }
  return;
}
void TrackParticleFitObject::resetTrajectorySecondDerivatives() const {
  for (int i=0; i<3; i++)
    for (int j=0; j<NPAR; j++)
      for (int k=0; k<NPAR; k++) {
        trajectoryStartSecondDerivatives[i][j][k]=0;
        trajectoryEndSecondDerivatives[i][j][k]=0;
      }
  return;
}

void TrackParticleFitObject::resetNormalFirstDerivatives() const {
  for (int i=0; i<3; i++)
    for (int j=0; j<NPAR; j++)
      normalFirstDerivatives[i][j]=0;
  return;
}

void TrackParticleFitObject::resetNormalSecondDerivatives() const {
  for (int i=0; i<3; i++)
    for (int j=0; j<NPAR; j++)
      for (int k=0; k<NPAR; k++)
        normalSecondDerivatives[i][j][k]=0;
  return;
}

// B field in Tesla
double TrackParticleFitObject::bfield = 3.5;

double TrackParticleFitObject::setBfield (double bfield_) {
  //  invalidateCache(); <-- not possible since this is static function
  return bfield = bfield_;
}

double TrackParticleFitObject::getBfield () {
  return bfield;
}

int TrackParticleFitObject::getCharge() const {
  updateCache(); 
  return omega>0 ? 1 : -1;
}

ThreeVector TrackParticleFitObject::getVertex (int ivertex) const {
  updateCache();
  return ivertex==0 ? trajectoryPointAtStart : trajectoryPointAtEnd ;
}


JBLHelix TrackParticleFitObject::getJBLHelix () {
  updateCache();
  // this is the helix wrt the track reference point
  //   mutable double phi0  ;
  //   mutable double omega ;
  //   mutable double tanl  ;
  //   mutable double d0    ;
  //   mutable double z0    ;
  // convert ILD parameters to those of JBL
  // daniel took this from http://ilcsoft.desy.de/MarlinTPC/current/doc/html/classmarlintpc_1_1simpleHelix.html
  // assume that JBL is same as what this page calls "perigee" parametrisation
  double kappa = -omega;
  double theta = M_PI/2. - atan(tanl);
  double dca = -d0; // do we have to worry about sign?

  return JBLHelix(kappa, phi0, theta, dca, z0);

}


JBLHelix TrackParticleFitObject::getJBLHelix (const ThreeVector & newRefPoint ) {
  updateCache();

  // equations taken from L3 internal note 1666

  // this is the helix wrt the track reference point
  //   mutable double phi0  ;
  //   mutable double omega ; // this is C in note 1666
  //   mutable double tanl  ;
  //   mutable double d0    ; // this is delta in note
  //   mutable double z0    ;

  // first convert the ILD parameters to those at the new reference point
  // taken from L3 internal note 1666  

  double DeltaX = newRefPoint.getX() - trackReferencePoint.getX();
  double DeltaY = newRefPoint.getY() - trackReferencePoint.getY();
  double DeltaZ = newRefPoint.getY() - trackReferencePoint.getY();

  // if new reference point is vert close to existing one, just use that one - it's simpler
  if ( sqrt( DeltaX*DeltaX + DeltaY*DeltaY + DeltaZ*DeltaZ ) < 1e-6 ) 
    return getJBLHelix();

  // now convert

  double R = 1./omega;

  double new_phi0 = atan2 ( sin(phi0) - DeltaX/(R - d0) , cos(phi0) + DeltaY/(R - d0) );
  double new_d0 = R - ( R - d0 )*sqrt( 1 + 
				       ( -2.*DeltaX*sin(phi0) + 2.*DeltaY*cos(phi0) )/( R - d0 ) + 
				       ( pow(DeltaX, 2.) + pow(DeltaY, 2.) )/pow ( R - d0 , 2 ) );

  double DeltaPhi = new_phi0 - phi0;
  // put in region -pi -> pi
  while ( DeltaPhi >  M_PI ) DeltaPhi-=2*M_PI;
  while ( DeltaPhi < -M_PI ) DeltaPhi+=2*M_PI;
  double sincDeltaPhi = fabs(DeltaPhi)>1e-10 ? sin(DeltaPhi)/DeltaPhi : 1.0;
  double s = (DeltaX*cos(phi0) + DeltaY*sin(phi0))/sincDeltaPhi;

  double new_z0 = z0 + s*tanl;

  // convert ILD parameters to those of JBL
  // daniel took this from http://ilcsoft.desy.de/MarlinTPC/current/doc/html/classmarlintpc_1_1simpleHelix.html
  // assume that JBL is same as what this page calls "perigee" parametrisation
  double kappa = -omega;
  double theta = M_PI/2. - atan(tanl);
  // double dca = -d0; // do we have to worry about sign?

  return JBLHelix( kappa , new_phi0, theta, -new_d0, new_z0);

}


/// Set start (i=0) or stop (i=1) vertex to a point as close as possible to given point
void TrackParticleFitObject::setVertex (int ivertex,         ///< Vertex number: 0=start, 1=stop
					const TwoVector& v   ///< Vertex position // DANIEL why is this 2d vertex?
					) {

  //  cout << "hello from TrackParticleFitObject::setVertex " << ivertex << " : " << v.getX() << " " << v.getY() << endl;

  updateCache();

  // find value of s for PCA to vertex v
  // vertex is in absolute position; vertex is wrt reference point
  JBLHelix jblh = getJBLHelix();

  TwoVector relPos( v.getX() - trackReferencePoint.getX() ,  v.getY() - trackReferencePoint.getY() );
  double best_s = jblh.getClosestS (relPos); 

  //  cout << "abs pos " << v.getX() << " " << v.getY() << " rel pos = " << relPos.getX() << " " << relPos.getY() << " bestS = " << best_s << endl;

  if ( ivertex==0 ) { // start vertex
    s_start = best_s;
  } else if (ivertex==1 ) { // end vertex
    s_end = best_s;
  } else {
    cout << "invalid vertex number " << ivertex << " only 0 (start), 1 (end) allowed" << endl;
    assert(0);
  }

  cachevalid=false;

  return;
}

/// Get derivative of vertex w.r.t. parameter ilocal 
ThreeVector TrackParticleFitObject::getVertexDerivative (int ivertex,       ///< vertex number: 0=start, 1=stop
							 int ilocal         ///< Local parameter number
							 ) const { 

  updateCache();

  ThreeVector vtxDer(0,0,0);

  if ( ivertex==0 ) {
    vtxDer.setValues( getTrajectoryStartFirstDerivatives(0, ilocal),
		      getTrajectoryStartFirstDerivatives(1, ilocal),
		      getTrajectoryStartFirstDerivatives(2, ilocal) );
  } else if ( ivertex==1 ) {
    vtxDer.setValues( getTrajectoryEndFirstDerivatives(0, ilocal),
		      getTrajectoryEndFirstDerivatives(1, ilocal),
		      getTrajectoryEndFirstDerivatives(2, ilocal) );
  } else {
    cout << "invalid vertex number " << ivertex << " only 0 (start), 1 (end) allowed" << endl;
    assert( 0 ); 
  }

  return vtxDer;
}

FourVector TrackParticleFitObject::getMomentum (int ivertex) const {
  // get the 4mom at start or end vertex

  //  cout << "hello0 " << ivertex << endl;

  updateCache();

  //  cout << "hello1" << endl;

  ThreeVector threeMom(0,0,0);
  if ( ivertex==0 ) {
    threeMom = momentumAtStart;
  } else if ( ivertex==1 ) {
    threeMom = momentumAtEnd;
  } else {
    cout << "invalid vertex number " << ivertex << " only 0 (start), 1 (end) allowed" << endl;
    assert( 0 );
  }

  //  cout << threeMom << " " << mass << endl;

  return FourVector(threeMom, mass);
}


bool TrackParticleFitObject::fixVertexParam (int ivertex, bool fix) {
  bool bb(false);
  if ( ivertex==0 ) {
    bb=fixParam( iStart, fix );
  } else if ( ivertex==1 ) {
    bb=fixParam( iEnd, fix );
  } else {
    assert(0);
  }
  return bb;
}

bool TrackParticleFitObject::releaseVertexParam (int ivertex) {    ///< Vertex number: 0=start, 1=stop
  assert( ivertex==0 || ivertex==1 );
  return fixVertexParam (ivertex, false);
}

