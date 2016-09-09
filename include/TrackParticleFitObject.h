/*! \file
 *  \brief Declares class TrackParticleFitObject
 *
 *
 */

#ifndef __TRACKPARTICLEFITOBJECT_H
#define __TRACKPARTICLEFITOBJECT_H

#include "ParticleFitObject.h"
#include "TwoVector.h"
#include "JBLHelix.h"

#undef NDEBUG
#include <cassert>

#include "EVENT/TrackState.h"
#include "EVENT/Track.h"

class TrackParticleFitObject : public ParticleFitObject {
 public:

  TrackParticleFitObject( const EVENT::Track*      trk, double m);
  TrackParticleFitObject( const EVENT::TrackState* trk, double m);
  TrackParticleFitObject( const double* _ppars, const double* _cov, double m, const double* refPt_=0);

  virtual ~TrackParticleFitObject();

  TrackParticleFitObject (const TrackParticleFitObject& rhs
			  );

  TrackParticleFitObject& operator= (const TrackParticleFitObject& rhs
				     );

  /// Return a new copy of itself
  virtual TrackParticleFitObject *copy() const;

  /// Assign from anther object, if of same type
  virtual TrackParticleFitObject& assign (const BaseFitObject& source   ///< The source object
                                   );

  /// Get name of parameter ilocal
  virtual const char *getParamName (int ilocal     ///< Local parameter number
                                    ) const;

  /// Read values from global vector, readjust vector; return: significant change
  virtual bool   updateParams (double p[],   ///< The parameter vector
                               int idim      ///< Length of the vector
                               );

  // these depend on actual parametrisation!
  virtual double getDPx(int ilocal) const;
  virtual double getDPy(int ilocal) const;
  virtual double getDPz(int ilocal) const;
  virtual double getDE(int ilocal) const;

  virtual double getFirstDerivative_Meta_Local( int iMeta, int ilocal , int metaSet ) const; // derivative of intermediate variable iMeta wrt local parameter ilocal
  virtual double getSecondDerivative_Meta_Local( int iMeta, int ilocal , int jlocal, int metaSet ) const; // derivative of intermediate variable iMeta wrt local parameter ilocal

  virtual double getChi2 () const;

  virtual int getNPar() const {return NPAR;}

  virtual ThreeVector getTrackPlaneNormal() const;
  virtual ThreeVector getTrackPcaVector() const;

  virtual int getCharge() const;

  /// Set the B field for all tracks
  static double setBfield (double bfield_             ///< New Value of B field (in Tesla)
			   );
  
  /// Get the B field for all tracks (in Tesla)
  inline static double getBfield ();
  
  /// Global B field in Tesla(!)
  static double bfield;

  static const double omega_pt_conv;
  static const double maxpt;



  //--------------------------------------------------
  // extra methods needed by VertrexFitObject
  //--------------------------------------------------

  /// Get start (i=0) or stop (i=1) vertex 
  virtual ThreeVector getVertex (int ivertex) const;

  /// Set start (i=0) or stop (i=1) vertex to a point as close as possible to given point
  virtual void setVertex (int ivertex,         ///< Vertex number: 0=start, 1=stop
			  const TwoVector& v   ///< Vertex position // DANIEL why is this 2d vertex?
			  );

  /// Get momentum at vertex
  virtual FourVector getMomentum (int ivertex) const;

  /// Set parameters such that track passes through a vertex with a given 4-momentum; return=success
  // only needed for unmeasured tracks, I think
  virtual bool  setParameters (int ivertex,                  ///< Vertex number: 0=start, 1=stop
			       const ThreeVector& vertex,    ///< Vertex position
			       const FourVector& momentum,   ///< Four-momentum
			       double charge_                ///< Charge (signed, in units of e)
                               ) {assert(0);}

  /// Get helix that is tangential at a certain arc length s
  //  virtual JBLHelix getTangentialHelix (double s     ///<  Arc length
  //				       ) {assert(0);}

  /// Fix parameter(s) pertaining to a vertex, or release it
  virtual bool fixVertexParam (int ivertex,    ///< Vertex number
			       bool fix=true  ///< fix if true, release if false
			       );

  /// Release parameter(s) pertaining to a vertex
  virtual bool releaseVertexParam (int ivertex    ///< Vertex number: 0=start, 1=stop
				   );


  /// Get derivative of vertex w.r.t. parameter ilocal 
  virtual ThreeVector getVertexDerivative (int ivertex,       ///< vertex number: 0=start, 1=stop
					   int ilocal         ///< Local parameter number
					   ) const;

  // this one probably useful

  /// Get JBLhelix
  virtual JBLHelix getJBLHelix (); // with respect to the reference point of the track
  virtual JBLHelix getJBLHelix (const ThreeVector& newRefPoint); // wrt a general reference point


  //---------------- end extra needed methods for vertexfitobject

  enum {iD0=0, iPhi0, iOmega, iZ0, iTanL, iStart, iEnd, NPAR};
  
 protected:

  static const double parfact[NPAR];

  virtual void initialise( const double* _pars, const double* _cov, double m);

  void updateCache() const;
  void updateMomentumDerivatives() const;
  void updateNormalDerivatives() const;
  void updateTrajectoryDerivatives() const;

  mutable ThreeVector trackReferencePoint;
  mutable ThreeVector trackPlaneNormal;
  mutable ThreeVector trackPcaVector;

  mutable ThreeVector trajectoryPointAtPCA;
  mutable ThreeVector trajectoryPointAtStart;
  mutable ThreeVector trajectoryPointAtEnd;

  mutable ThreeVector momentumAtPCA;
  mutable ThreeVector momentumAtStart;
  mutable ThreeVector momentumAtEnd;

  mutable double momentumFirstDerivatives[4][NPAR];
  mutable double momentumSecondDerivatives[4][NPAR][NPAR];

  mutable double normalFirstDerivatives[3][NPAR];
  mutable double normalSecondDerivatives[3][NPAR][NPAR];

  mutable double trajectoryStartFirstDerivatives[3][NPAR];
  mutable double trajectoryStartSecondDerivatives[3][NPAR][NPAR];

  mutable double trajectoryEndFirstDerivatives[3][NPAR];
  mutable double trajectoryEndSecondDerivatives[3][NPAR][NPAR];

  mutable double phi0  ;
  mutable double omega ;
  mutable double tanl  ;
  mutable double d0    ;
  mutable double z0    ;
  mutable double s_start;
  mutable double s_end  ;

  mutable double chi2;

  void   resetMomentumFirstDerivatives() const;
  void   resetMomentumSecondDerivatives() const;

  void   resetNormalFirstDerivatives() const;
  void   resetNormalSecondDerivatives() const;

  void   resetTrajectoryFirstDerivatives() const;
  void   resetTrajectorySecondDerivatives() const;

  void   setMomentumFirstDerivatives (int iMeta, int jLocal, double x) const;
  void   setMomentumSecondDerivatives(int iMeta, int jLocal, int kLocal, double x) const;

  void   setNormalFirstDerivatives (int iMeta, int jLocal, double x) const;
  void   setNormalSecondDerivatives(int iMeta, int jLocal, int kLocal, double x) const;

  void   setTrajectoryStartFirstDerivatives (int iMeta, int jLocal, double x) const;
  void   setTrajectoryStartSecondDerivatives(int iMeta, int jLocal, int kLocal, double x) const;

  void   setTrajectoryEndFirstDerivatives (int iMeta, int jLocal, double x) const;
  void   setTrajectoryEndSecondDerivatives(int iMeta, int jLocal, int kLocal, double x) const;

  virtual double getMomentumFirstDerivatives (int iMeta, int jLocal) const;
  virtual double getMomentumSecondDerivatives(int iMeta, int jLocal, int kLocal) const;

  virtual double getNormalFirstDerivatives (int iMeta, int jLocal) const;
  virtual double getNormalSecondDerivatives(int iMeta, int jLocal, int kLocal) const;

  virtual double getTrajectoryStartFirstDerivatives (int iMeta, int jLocal) const;
  virtual double getTrajectoryStartSecondDerivatives(int iMeta, int jLocal, int kLocal) const;

  virtual double getTrajectoryEndFirstDerivatives (int iMeta, int jLocal) const;
  virtual double getTrajectoryEndSecondDerivatives(int iMeta, int jLocal, int kLocal) const;

};



#endif // __TRACKPARTICLEFITOBJECT_H

