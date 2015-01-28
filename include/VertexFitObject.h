/*! \file 
 *  \brief Declares class VertexFitObject
 *
 * \b Changelog:
 * - 6.12.04 BL First version
 *
 * \b CVS Log messages:
 * - $Log: VertexFitObject.h,v $
 * - Revision 1.5  2008/01/30 09:14:55  blist
 * - Preparations for NewtonFitter
 * -
 * - Revision 1.4  2007/09/17 12:50:15  blist
 * - Some parameters reordered
 * -
 * - Revision 1.3  2007/09/14 10:58:42  blist
 * - Better documentation,
 * - added PyConstraint::add1stDerivativesToMatrix,
 * - fixed bug in PzConstraint::add1stDerivativesToMatrix
 * -
 * - Revision 1.2  2007/09/13 13:33:06  blist
 * - Print methods return os
 * -
 *
 */ 

#ifndef __VERTEXFITOBJECT_H
#define __VERTEXFITOBJECT_H

#include "BaseFitObject.h"

#include <iostream>
#include <cassert>
#include <vector>

class ThreeVector;
class FourVector;
class TrackFitObject;
class BaseFitter;
class BaseConstraint;

// Class VertexFitObject
/// Class that represents vertices
/**
 *
 * A VertexFitObject represents a vertex, parametrized by
 * its coordinates (x, y, z).
 *
 * A VertexFitObject is a BaseFitObject, and therefore the vertex coordinates
 * can enter a fit, either as measured values for a measured (primary) vertex,
 * or as unmeasured values for a decay vertex.
 *
 * Additionally, a VertexFitObject keeps a list of outgoing tracks and knows
 * its incoming track (if any), which enables the VertexFitObject to set a list
 * of constraints (VertexConstraint and TrackMomentumConstraint objects)
 * for a fitter object.
 *
 * 
 * Author:Benno List, Jenny List
 * $Date: 2008/01/30 09:14:55 $
 * $Author: blist $
 */ 
class VertexFitObject: public BaseFitObject {
  public:
    /// Constructor
    VertexFitObject(const char *name_,
                    double x,
                    double y,
                    double z
                   );
    /// Copy Constructor
    VertexFitObject (const VertexFitObject& rhs);
    /// assignment
    VertexFitObject& operator= (const VertexFitObject& rhs);
    /// Return a new copy of itself
    virtual VertexFitObject *copy() const ;
    
    /// Assign from anther object, if of same type
    virtual VertexFitObject& assign (const BaseFitObject& source   ///< The source object
                                 ) ;
    /// Destructor
    virtual ~VertexFitObject();
    
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
    /// Set measured value of parameter ilocal; return=success
    virtual bool   setMParam (int ilocal,         ///< Local parameter number
                              double mpar_        ///< Measured parameter value
                             );  
    /// Set error of parameter ilocal
    virtual bool   setError (int ilocal,    ///< Local parameter number
                             double err_    ///< New error value
                             );
    /// Set covariance of parameters ilocal and jlocal
    virtual bool   setCov (int ilocal,    ///< Local parameter number
                           int jlocal,    ///< Local parameter number
                           double cov_    ///< New error value
                          );
    /// Set number of parameter ilocal in global list
    /// return true signals OK
    virtual bool setGlobalParNum (int ilocal, int iglobal); 
    
    /// Fix a parameter (fix=true), or release it (fix=false)
    virtual bool fixParam (int ilocal,    ///< Local parameter number
                           bool fix=true  ///< fix if true, release if false
                          );
    
    /// Get parameter ilocal
    virtual double getParam (int ilocal     ///< Local parameter number
                            ) const;
    /// Get measured value of parameter ilocal
    virtual double getMParam (int ilocal     ///< Local parameter number
                            ) const;
    /// Get error of parameter ilocal
    virtual double getError (int ilocal     ///< Local parameter number
                            ) const;
    /// Get covariance between parameters ilocal and jlocal
    virtual double getCov (int ilocal,    ///< Local parameter number
                           int jlocal     ///< Local parameter number
                          ) const;
    /// Get measured flag for parameter i
    virtual bool isParamMeasured (int ilocal) const;
    /// Get fixed flag for parameter i
    virtual bool isParamFixed (int ilocal) const;
    /// Get global parameter number of parameter ilocal
    virtual int getGlobalParNum(int ilocal) const;
    /// Get number of parameters of this FitObject
    virtual int getNPar() const;
    
    /// Get object's name
    virtual const char *getName () const;
    
    /// Set object's name
    virtual void setName (const char *name_);
    
     /// Get name of parameter ilocal
    virtual const char *getParamName (int ilocal     ///< Local parameter number
                                     ) const;
   
    
    /// Get vertex into  existing 3-vector
    virtual void getVertexEx (ThreeVector& p
                             ) const;
    /// Get vertex 
    virtual ThreeVector getVertex() const; 
    
    /// Get derivative of vertex w.r.t. parameter ilocal into  existing 3-vector
    virtual void getVertexDerivativeEx (int ilocal, 
                                        ThreeVector& p
                                       ) const;
    /// Get derivative of vertex w.r.t. parameter ilocal 
    virtual ThreeVector getVertexDerivative (int ilocal 
                                            ) const;

    /// Add covariance matrix elements to 
    /// global covariance matrix of size idim x idim
    virtual void addToGlobCov(double *cov, int idim) const; 
        
    
    /// Get chi squared from measured and fitted parameters
    virtual double getChi2() const;
    /// Get derivative of chi squared w.r.t. parameter ilocal
    virtual double getDChi2DParam(int ilocal) const 
    {assert (0);}
    /// Get second derivative of chi squared w.r.t. parameters ilocal1 and ilocal2
    virtual double getD2Chi2DParam2(int ilocal1, int ilocal2) const {assert (0);};
    
    /// Add derivatives of chi squared to global covariance matrix
    virtual void addToGlobalChi2DerMatrix (double *M,   ///< Global covariance matrix
                                           int idim     ///< First dimension of global covariance matrix
                                           ) const 
    {assert (0);}
    /// Add derivatives of chi squared to global derivative matrix
    virtual void addToGlobalChi2DerVector (double *y,   ///< Vector of chi2 derivatives
                                           int idim     ///< Vector size 
                                           ) const 
    {assert (0);}
    
    /// Add derivatives to global covariance matrix
    virtual void addToGlobalDerMatrix (int idim, double c, double *M) const {assert (0);};
    
    /// print object to ostream
    virtual std::ostream& print (std::ostream& os   ///< The output stream
                                ) const;
    
    /// get track i
    virtual TrackFitObject *getTrack (int i
                                     ) const {return tracks[i].track;}
        
    /// add track
    virtual void addTrack (TrackFitObject *track,   ///< The track
                           bool inbound,            ///< Is track decaying vertex (inbound=true), or originating?
                           bool measured            ///< Is this a measured track?
                          );
    
    enum  {
      VX = 1, VY = 2, VZ = 4, VXY = 3, VXYZ=7,
      PX = 8, PY = 16, PZ = 32, PXY = 24, PXYZ = 56,
      E = 64, EPXYZ = 120,
      THE_FULL_MONTY = 127
    };
    
    /// Create constraints and add them to a BaseFit object
    virtual void addConstraints (BaseFitter& fitter, 
                                 int mask = THE_FULL_MONTY
                                );
    
    /// Estimate vertex position
    virtual ThreeVector estimatePosition ();
    
    /// Initialize this object and attatched tracks for fit with initial estimates
    virtual void initForFit();
            
  protected:
    /// generate vertex constraints
    virtual void addVertexConstraints (BaseFitter& fitter, int axis);
    /// generate momentum constraints
    virtual void addMomentumConstraint (BaseFitter& fitter, int axis);
  
    /// init covariance matrix to dummy values
    virtual void initCov();
    
    /// Calculate the inverse of the covariance matrix
    virtual bool calculateCovInv() const;
    
    /// Calculate chi2
    virtual void calculateChi2() const;
    /// Copy another VertexFitObject
    virtual void copy (const VertexFitObject& source);
  
  
    /// Number of parameters
    enum {NPAR = 3};
    
    /// fit parameters
    double par[NPAR];
    /// measured parameters
    double mpar[NPAR];
    /// errors
    double err[NPAR];
    /// measured flag
    bool measured[NPAR];
    /// fixed flag
    bool fixed[NPAR];
    /// global paramter number for each parameter
    int globalParNum [NPAR];
    /// local covariance matrix
    double cov [NPAR][NPAR];    
    /// inverse of local covariance matrix
    mutable double covinv [NPAR][NPAR];   
    /// flag for valid inverse covariance matrix
    mutable bool covinvvalid; 
    
    // Cache variables
    mutable double chi2;               ///< chi^2
    mutable double resid[NPAR];     ///< residuals
    mutable bool   chi2contr[NPAR]; ///< contributes to chi2?
    
    char *name;                        ///< object name
    
    struct TrackDescriptor {
      TrackFitObject *track;
      bool inbound;
      bool measured;
      TrackDescriptor (TrackFitObject *track_, bool inbound_, bool measured_)
        : track (track_), inbound (inbound_), measured (measured_) 
        {}
      TrackDescriptor ()
        : track (0), inbound (0), measured (0) 
        {}
    };
    
    typedef std::vector<TrackDescriptor> TContainer;
    typedef TContainer::iterator TIterator;
    TContainer tracks;
    typedef std::vector<BaseConstraint *> CContainer;
    typedef CContainer::iterator CIterator;
    CContainer constraints;
    
};
    
#endif // __VERTEXFITOBJECT_H

