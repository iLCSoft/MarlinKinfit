/*! \file 
 *  \brief Declares class TrackMomentumConstraint
 *
 * \b Changelog:
 * - 23.11.04 BL: First version (taken from PConstraint)
 *
 */ 

#ifndef __TRACKMOMENTUMCONSTRAINT_H
#define __TRACKMOMENTUMCONSTRAINT_H

#include "TrackConstraint.h"
#include "FourVector.h"

class ParticleFitObject;

class TrackMomentumConstraint : public TrackConstraint {
  public:
    TrackMomentumConstraint (double pxfact_, 
                             double pyfact_, 
                             double pzfact_,
                             double efact_=0, 
                             double value_ = 0
                             );
    TrackMomentumConstraint (int axis, 
                             double value_ = 0
                            );                          
    virtual ~TrackMomentumConstraint();
    virtual double getValue() const;
    virtual void getDerivatives (int idim, double der[]) const;
    virtual void add1stDerivativesToMatrix(int idim, double *M) const;
    virtual void add2ndDerivativesToMatrix(int idim, double *M, double lambda) const;
    
    virtual void addToGlobalDerMatrix (double lambda, int idim, double *M) const;
    
    virtual void invalidateCache() const;
  
  protected:
    void updateCache() const;
  
    FourVector factor;
    double value;
    
    mutable bool cachevalid;
    mutable int  nparams;
};

#endif // __TRACKMOMENTUMCONSTRAINT_H
