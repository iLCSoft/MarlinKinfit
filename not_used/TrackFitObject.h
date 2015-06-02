/*! \file 
 *  \brief Declares class TrackFitObject
 *
 * \b Changelog:
 * - 17.11.04 BL First version
 *
 * \b CVS Log messages:
 * - $Log: TrackFitObject.h,v $
 * - Revision 1.5  2009/05/22 12:12:07  blist
 * - New tracers
 * -
 * - Revision 1.4  2008/01/30 09:14:55  blist
 * - Preparations for NewtonFitter
 * -
 * - Revision 1.3  2007/09/17 12:50:15  blist
 * - Some parameters reordered
 * -
 * - Revision 1.2  2007/09/13 13:33:06  blist
 * - Print methods return os
 * -
 *
 */ 

#ifndef __TRACKFITOBJECT_H
#define __TRACKFITOBJECT_H

#include "BaseFitObject.h"

#include <iostream>

class TwoVector;
class ThreeVector;
class FourVector;
class JBLHelix;


// Class TrackFitObject
/// Abstract base class for track objects of kinematic fits
/**
 *
 * A TrackFitObject has usually 5 track parameters, e.g.
 * (kappa, phi0, theta, dca, z0) for a helix,
 * plus one or two parameters for a start and a possible stop
 * arc length, and a mass.
 * Therefore, the base class provides space to store
 * NPARMAX=8 parameters. 
 *
 * A TrackFitObject is able to calculate a point on its trajectory
 * for a given arclength s, a tangent vector at this point,
 * and a momentum. In general, it will need to know the b-field
 * (assumed to be along z) for this. The b-field can be set as 
 * static value for the whole TrackFitObject class.
 *
 * A TrackFitObject can also return a JBLHelix object that represents
 * a tangential helix at a certain arc length s. The helix object is 
 * able to find a point where the helix comes close to another helix,
 * which is needed to determine starting values for fits.
 * 
 * Author:Benno List, Jenny List
 * $Date: 2009/05/22 12:12:07 $
 * $Author: blist $
 *
 * \b Changelog:
 * - 17.11.04 BL First version
 * - 5.1.05 BL: Add getTangentialHelix
 */ 


class TrackFitObject: public BaseFitObject {
  public:
    /// Default constructor
    TrackFitObject(const char *name_ = 0                   ///< Optional name
                  );
    /// Copy Constructor
    TrackFitObject (const TrackFitObject& rhs              ///< Source object
                   );
    /// Assignment
    TrackFitObject& operator= (const TrackFitObject& rhs   ///< Source object
                   );
    /// Return a new copy of itself
    virtual TrackFitObject *copy() const = 0;
    
    /// Assign from anther object, if of same type
    virtual TrackFitObject& assign (const BaseFitObject& source   ///< The source object
                                 ) = 0;
    
    /// Destructor
    virtual ~TrackFitObject();
    
    /// Set value and measured flag of parameter ilocal; return=success
    virtual bool   setParam (int ilocal,         ///< Local parameter number
                             double par_,        ///< New parameter value
                             bool measured_,     ///< New "measured" flag
                             bool fixed_ = false ///< New "fixed" flag
                            );  
    /// Set value of parameter ilocal; return=success
    virtual bool   setParam (int ilocal,         ///< Local parameter number
                             double par_         ///< New parameter value
                            );  

    /// Set parameters such that track passes through a vertex with a given 4-momentum; return=success
    virtual bool  setParameters (int ivertex,                  ///< Vertex number: 0=start, 1=stop
                                const ThreeVector& vertex,    ///< Vertex position
                                const FourVector& momentum,   ///< Four-momentum
                                double charge_                ///< Charge (signed, in units of e)
                               ) = 0; 
    /// Set measured value of parameter ilocal; return=success
    virtual bool   setMParam (int ilocal,         ///< Local parameter number
                              double mpar_        ///< New parameter value
                             );  
    /// Set error of parameter ilocal; return=success
    virtual bool   setError (int ilocal,          ///< Local parameter number
                             double err_          ///< New error value
                             );
    /// Set covariance of parameters ilocal and jlocal; return=success
    virtual bool   setCov (int ilocal,            ///< Local parameter number
                           int jlocal,            ///< Local parameter number
                           double cov_            ///< New error value
                          );
    /// Set number of parameter ilocal in global list
    /// return true signals OK
    virtual bool setGlobalParNum (int ilocal,     ///< Local parameter number
                                  int iglobal     ///< Global parameter number
                                 ); 
    
    /// Fix a parameter (fix=true), or release it (fix=false)
    virtual bool fixParam (int ilocal,            ///< Local parameter number
                           bool fix=true          ///< fix if true, release if false
                          );
    
    /// Get current value of parameter ilocal
    virtual double getParam (int i                ///< Local parameter number
                            ) const;
    /// Get measured value of parameter ilocal
    virtual double getMParam (int i              ///< Local parameter number
                            ) const;
    /// Get error of parameter ilocal
    virtual double getError (int ilocal          ///< Local parameter number
                            ) const;
    /// Get covariance between parameters ilocal and jlocal
    virtual double getCov (int ilocal,           ///< Local parameter number 1
                           int jlocal            ///< Local parameter number 2
                          ) const;
    /// Get measured flag for parameter i
    virtual bool isParamMeasured (int ilocal     ///< Local parameter number
                                  ) const;
    /// Get fixed flag for parameter i
    virtual bool isParamFixed (int ilocal        ///<  Local parameter number
                              ) const;
    /// Get global parameter number of parameter ilocal
    virtual int getGlobalParNum(int ilocal       ///<  Local parameter number
                               ) const;
    /// Get number of parameters of this FitObject
    virtual int getNPar() const = 0;
    
    /// Get object's name
    virtual const char *getName () const;
    
    /// Set object's name
    virtual void setName (const char *name_       ///< Pointer to new name string 
                         );
    
    /// Get point along trajectory into existing 3-vector
    virtual void getTrajectoryPointEx (double s,              ///< Arc length
                                       ThreeVector& p         ///< Target 3-vector
                                      ) const = 0; 
    /// Get point along trajectory
    virtual ThreeVector getTrajectoryPoint (double s             ///<  Arc length
                                           ) const; 
    /// Get start (i=0) or stop (i=1) vertex into existing 3-vector
    virtual void getVertexEx (int ivertex,              ///< Vertex number (0: start, 1: stop)
                              ThreeVector& p            ///< Target 3-vector
                             ) const = 0; 
    /// Get start (i=0) or stop (i=1) vertex 
    virtual ThreeVector getVertex (int ivertex    ///< vertex number: 0=start, 1=stop
                                  ) const; 

    /// Set start (i=0) or stop (i=1) vertex to a point as close as possible to given point
    virtual void setVertex (int ivertex,         ///< Vertex number: 0=start, 1=stop
                            const TwoVector& v   ///< Vertex position
                           ) = 0; 
    
    /// Get derivative of trajectory w.r.t. parameter ilocal into  existing 3-vector
    virtual void getTrajectoryDerivativeEx (double s,            ///< Arc length
                                            int ilocal,          ///< Local parameter number
                                            ThreeVector& p       ///< Target 3-vector
                                           ) const = 0;
    /// Get derivative of trajectory w.r.t. parameter ilocal 
    virtual ThreeVector getTrajectoryDerivative (double s,       ///< Arc length
                                                 int ilocal      ///< Local parameter number
                                                ) const;
    /// Get derivative of vertex w.r.t. parameter ilocal into  existing 3-vector
    virtual void getVertexDerivativeEx (int ivertex,            ///< vertex number: 0=start, 1=stop
                                        int ilocal,             ///< Local parameter number
                                        ThreeVector& p          ///< Target 3-vector
                                       ) const = 0;
    /// Get derivative of vertex w.r.t. parameter ilocal 
    virtual ThreeVector getVertexDerivative (int ivertex,       ///< vertex number: 0=start, 1=stop
                                             int ilocal         ///< Local parameter number
                                            ) const;

    /// Get momentum along trajectory into existing 4-vector
    virtual void getMomentumAtTrajectoryEx (double s,           ///< Arc length
                                            FourVector& p       ///< Target 3-vector
                                           ) const = 0; 
    /// Get momentum along trajectory
    virtual FourVector getMomentumAtTrajectory (double s        ///< Arc length
                                               ) const; 
    /// Get momentum at vertex into existing 4-vector
    virtual void getMomentumEx (int ivertex,                    ///<  Vertex number: 0=start, 1=stop
                                FourVector& p                   ///< Target 3-vector
                               ) const = 0; 
    /// Get momentum at vertex
    virtual FourVector getMomentum (int ivertex                 ///< Vertex number: 0=start, 1=stop
                                   ) const; 

    /// Get charge in units of e
    virtual double getCharge () const = 0; 

    /// Get mass in GeV
    virtual double getMass () const = 0; 

    /// Get momentum along trajectory into existing 4-vector
    virtual void getMomentumDerivativeAtTrajectoryEx (double s,      ///< Arc length
                                                      int ilocal,    ///< Local parameter number 
                                                      FourVector& p  ///< Target 3-vector
                                                     ) const = 0; 
    /// Get momentum along trajectory
    virtual FourVector getMomentumDerivativeAtTrajectory (double s,   ///< Arc length
                                                          int ilocal  ///< Local parameter number 
                                                         ) const; 
    /// Get derivative of momentum w.r.t. parameter ilocal into  existing 4-vector
    virtual void getMomentumDerivativeEx (int ivertex,           ///< Vertex number: 0=start, 1=stop
                                          int ilocal,            ///< Local parameter number 
                                          FourVector& p          ///< Target 3-vector
                                         ) const = 0;
    /// Get derivative of momentum w.r.t. parameter ilocal 
    virtual FourVector getMomentumDerivative (int ivertex,       ///< Vertex number: 0=start, 1=stop
                                              int ilocal         ///< Local parameter number
                                             ) const;
    /// Get s (arclength in r/phi) of vertex ivertex
    virtual double getArcLength (int ivertex                           ///<  Vertex number: 0=start, 1=stop
                                ) const = 0;
    /// Add covariance matrix elements to 
    /// global covariance matrix of size idim x idim
    virtual void addToGlobCov(double *globvcov,                       ///< Covariance matrix
                              int idim                                ///< Dimension of covariance matrix
                              ) const; 
        
    
    /// Get chi squared from measured and fitted parameters
    virtual double getChi2() const;
    /// Get derivative of chi squared w.r.t. parameter ilocal
    virtual double getDChi2DParam(int ilocal                ///<  Local parameter number
                                    ) const = 0;
    /// Get second derivative of chi squared w.r.t. parameters ilocal1 and ilocal2
    virtual double getD2Chi2DParam2(int ilocal,               ///<  Local parameter number 1
                                    int jlocal                ///<  Local parameter number 2
                                   ) const = 0;
    
    /// Add derivatives of chi squared to global covariance matrix
    virtual void addToGlobalChi2DerMatrix (double *M,   ///< Global covariance matrix
                                           int idim     ///< First dimension of global covariance matrix
                                           ) const = 0;
    
    /// Get helix that is tangential at a certain arc length s
    virtual JBLHelix getTangentialHelix (double s     ///<  Arc length
                                        ) = 0;
    /// print object to ostream
    virtual std::ostream& print (std::ostream& os              ///< Output stream
                                ) const;
    
    /// Invalidate the cache
    virtual void invalidateCache();  
    
    
    /// Set the B field for all tracks
    static double setBfield (double bfield_             ///< New Value of B field (in Tesla)
                            );
  
    /// Get the B field for all tracks (in Tesla)
    inline static double getBfield ();
    
    /// Global B field in Tesla(!)
    static double bfield;
    
    /// Fix parameter(s) pertaining to a vertex, or release it
    virtual bool fixVertexParam (int ivertex,    ///< Vertex number: 0=start, 1=stop
                                 bool fix=true  ///< fix if true, release if false
                                 ) = 0;
    /// Release parameter(s) pertaining to a vertex
    virtual bool releaseVertexParam (int ivertex    ///< Vertex number: 0=start, 1=stop
                                    );
    /// Print Covariance matrix
    virtual void printCov (std::ostream& os              ///< Output stream
                          ) const;
  protected:
    /// Init covariance matrix to dummy values
    virtual void initCov();
    /// Check covariance matrix to dummy values
    virtual void checkCov();
    
    /// Calculate the inverse of the covariance matrix
    virtual bool calculateCovInv() const;
    
    /// Calculate chi2
    virtual void calculateChi2() const;
    /// Copy another TrackFitObject
    virtual void copy (const TrackFitObject& source              ///< Source object
                   );
  
  
    /// Maximum number of parameters
    enum {NPARMAX = 8};
    
    /// fit parameters
    double par[NPARMAX];
    /// measured parameters
    double mpar[NPARMAX];
    /// errors
    double err[NPARMAX];
    /// measured flag
    bool measured[NPARMAX];
    /// fixed flag
    bool fixed[NPARMAX];
    /// global paramter number for each parameter
    int globalParNum [NPARMAX];
    /// local covariance matrix
    double cov [NPARMAX][NPARMAX];    
    /// inverse of local covariance matrix
    mutable double covinv [NPARMAX][NPARMAX];   
    /// flag for valid cache
    mutable bool cachevalid; 
    /// flag for valid inverse covariance matrix
    mutable bool covinvvalid; 
    
    // Cache variables
    mutable double chi2;               ///< chi^2
    mutable double resid[NPARMAX];     ///< residuals
    mutable bool   chi2contr[NPARMAX]; ///< contributes to chi2?
    
    char *name;                        ///< object name (name string must be allocated and deleted outside object!)
    
};
    
inline std::ostream& operator<< (std::ostream& os,              ///< Output stream
                                 const TrackFitObject& tfo      ///< Object to be printed
                                ) {
  tfo.print(os);
  return os;
}

double TrackFitObject::getBfield() {
  return bfield;
}

#endif // __TRACKFITOBJECT_H

